#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import json

import healpy as hp
import matplotlib.pylab as plt
import numpy as np
import rubin_scheduler
import rubin_scheduler.scheduler.basis_functions as bf
import rubin_scheduler.scheduler.detailers as detailers
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils import iers
from rubin_scheduler.scheduler import sim_runner
from rubin_scheduler.site_models import _read_fields
from rubin_scheduler.scheduler.model_observatory import ModelObservatory
from rubin_scheduler.scheduler.schedulers import CoreScheduler, SimpleFilterSched
from rubin_scheduler.scheduler.surveys import (
    BlobSurvey,
    GreedySurvey,
    LongGapSurvey,
    ScriptedSurvey,
    BaseMarkovSurvey,
    generate_ddf_scheduled_obs,
)
from rubin_scheduler.scheduler.utils import (
    ConstantFootprint,
    EuclidOverlapFootprint,
    make_rolling_footprints,
    scheduled_observation,
    comcam_tessellate,
    gnomonic_project_toxy,
    tsp_convex,
    SimTargetooServer,
    TargetoO,
)
from rubin_scheduler.utils import (
    _hpid2_ra_dec,
    _ra_dec2_hpid,
    _approx_ra_dec2_alt_az,
    _angular_separation,
)

from tabulate import tabulate

# So things don't fail on hyak
iers.conf.auto_download = False
# XXX--note this line probably shouldn't be in production
iers.conf.auto_max_age = None


class ToO_scripted_survey(ScriptedSurvey, BaseMarkovSurvey):
    """If there is a new ToO event, generate a set of scripted observations to try and follow it up.

    Parameters
    ----------
    

    times : list of floats
        The times after the detection that observations should be attempted (hours)

    alt_min : float
        Do not attempt observations below this limit (degrees). Note the telescope alt limit is 20
        degrees, however, slew and filter change time means alt_min here should be set higher
        (otherwise, target will pass altitude check, but then fail to observe by the time the
        telescope gets there).

    """

    def __init__(
        self,
        basis_functions,
        followup_footprint=None,
        nside=32,
        reward_val=1e6,
        times=[0, 1, 2, 4, 24],
        filters_at_times=["gz", "gz", "gz", "gz", "gy"],
        nvis=[1, 1, 1, 1, 6],
        exptime=30.0,
        camera="LSST",
        survey_name="ToO",
        flushtime=2.0,
        mjd_tol=1.0 / 24.0,
        dist_tol=0.5,
        alt_min=25.0,
        alt_max=85.0,
        HA_min=5,
        HA_max=19,
        ignore_obs="dummy",
        dither=True,
        seed=42,
        npositions=7305,
        n_snaps=2,
        n_usnaps=1,
        id_start=1,
        detailers=None,
    ):
        # Figure out what else I need to super here

        self.basis_functions = basis_functions
        self.survey_name = survey_name
        self.followup_footprint = followup_footprint
        self.last_event_id = -1
        self.night = -1
        self.reward_val = reward_val
        self.times = np.array(times) / 24.0  # to days
        self.filters_at_times = filters_at_times
        self.exptime = exptime
        self.nvis = nvis
        self.n_snaps = n_snaps
        self.n_usnaps = n_usnaps
        self.nside = nside
        self.flushtime = flushtime / 24.0
        self.mjd_tol = mjd_tol
        self.dist_tol = np.radians(dist_tol)
        self.alt_min = np.radians(alt_min)
        self.alt_max = np.radians(alt_max)
        self.HA_min = HA_min
        self.HA_max = HA_max
        self.ignore_obs = ignore_obs
        self.extra_features = {}
        self.extra_basis_functions = {}
        self.detailers = []
        self.dither = dither
        self.id_start = id_start
        self.detailers = detailers

        self.camera = camera
        # Load the OpSim field tesselation and map healpix to fields
        if self.camera == "LSST":
            ra, dec = _read_fields()
            self.fields_init = np.empty(
                ra.size, dtype=list(zip(["RA", "dec"], [float, float]))
            )
            self.fields_init["RA"] = ra
            self.fields_init["dec"] = dec
        elif self.camera == "comcam":
            self.fields_init = comcam_tessellate()
        else:
            ValueError('camera %s unknown, should be "LSST" or "comcam"' % camera)
        self.fields = self.fields_init.copy()

        self.hp2fields = np.array([])
        self._hp2fieldsetup(self.fields["RA"], self.fields["dec"])

        # Initialize the list of scripted observations
        self.clear_script()

        # Generate and store rotation positions to use.
        # This way, if different survey objects are seeded the same, they will
        # use the same dither positions each night
        rng = np.random.default_rng(seed)
        self.lon = rng.random(npositions) * np.pi * 2
        # Make sure latitude points spread correctly
        # http://mathworld.wolfram.com/SpherePointPicking.html
        self.lat = np.arccos(2.0 * rng.random(npositions) - 1.0)
        self.lon2 = rng.random(npositions) * np.pi * 2

    def _check_list(self, conditions):
        """Check to see if the current mjd is good"""
        observation = None
        if self.obs_wanted is not None:
            # Scheduled observations that are in the right time window and have not been executed
            in_time_window = np.where(
                (self.mjd_start < conditions.mjd)
                & (self.obs_wanted["flush_by_mjd"] > conditions.mjd)
                & (~self.obs_wanted["observed"])
            )[0]

            if np.size(in_time_window) > 0:
                pass_checks = self._check_alts_ha(
                    self.obs_wanted[in_time_window], conditions
                )
                matches = in_time_window[pass_checks]
            else:
                matches = []

            if np.size(matches) > 0:
                # If we have something in the current filter, do that, otherwise whatever is first
                # in_filt = np.where(self.obs_wanted[matches]['filter'] == conditions.current_filter)[0]
                # if np.size(in_filt) > 0:
                #    indx = matches[in_filt[0]]
                # else:
                #    indx = matches[0]
                # observation = self._slice2obs(self.obs_wanted[indx])
                observation = self._slice2obs(self.obs_wanted[matches[0]])

        return observation

    def flush_script(self, conditions):
        """Remove things from the script that aren't needed anymore"""
        if self.obs_wanted is not None:
            still_relevant = np.where(
                (self.obs_wanted["observed"] == False)
                & (self.obs_wanted["flush_by_mjd"] < conditions.mjd)
            )[0]
            if np.size(still_relevant) > 0:
                observations = self.obs_wanted[still_relevant]
                self.set_script(observations)
            else:
                self.clear_script()

    def _new_event(self, target_o_o, conditions):
        """A new ToO event, generate any observations for followup"""
        # flush out any old observations or ones that have been completed
        self.flush_script(conditions)
        # Check that the event center is in the footprint we want to observe
        hpid_center = _ra_dec2_hpid(
            self.nside, target_o_o.ra_rad_center, target_o_o.dec_rad_center
        )
        if self.followup_footprint[hpid_center] > 0:
            target_area = self.followup_footprint * target_o_o.footprint
            # generate a list of pointings for that area
            hpid_to_observe = np.where(target_area > 0)[0]

            # Check if we should spin the tesselation for the night.
            if self.dither & (conditions.night != self.night):
                self._spin_fields(conditions)
                self.night = conditions.night + 0

            field_ids = np.unique(self.hp2fields[hpid_to_observe])

            # Put the fields in a good order. Skipping dither positions for now.
            better_order = order_observations(
                self.fields["RA"][field_ids], self.fields["dec"][field_ids]
            )
            ras = self.fields["RA"][field_ids[better_order]]
            decs = self.fields["dec"][field_ids[better_order]]

            # Figure out an MJD start time for the object if it is still rising and low.
            alt, az = _approx_ra_dec2_alt_az(
                target_o_o.ra_rad_center,
                target_o_o.dec_rad_center,
                conditions.site.latitude_rad,
                None,
                conditions.mjd,
                lmst=np.max(conditions.lmst),
            )
            HA = np.max(conditions.lmst) - target_o_o.ra_rad_center * 12.0 / np.pi

            if (HA < self.HA_max) & (HA > self.HA_min):
                t_to_rise = (self.HA_max - HA) / 24.0
                mjd0 = conditions.mjd + t_to_rise
            else:
                mjd0 = conditions.mjd + 0.0

            obs_list = []
            for time, filternames, nv in zip(
                self.times, self.filters_at_times, self.nvis
            ):
                for filtername in filternames:
                    # Subsitute y for z if needed
                    if (filtername == "z") & (
                        filtername not in conditions.mounted_filters
                    ):
                        filtername = "y"
                    for i in range(nv):
                        if filtername in conditions.mounted_filters:
                            if filtername == "u":
                                nexp = self.n_usnaps
                            else:
                                nexp = self.n_snaps

                            obs = scheduled_observation(ras.size)
                            obs["RA"] = ras
                            obs["dec"] = decs
                            obs["mjd"] = mjd0 + time
                            obs["flush_by_mjd"] = mjd0 + time + self.flushtime
                            obs["exptime"] = self.exptime
                            obs["nexp"] = nexp
                            obs["filter"] = filtername
                            obs[
                                "rotSkyPos"
                            ] = 0  # XXX--maybe throw a rotation detailer in here
                            obs["mjd_tol"] = self.mjd_tol
                            obs["dist_tol"] = self.dist_tol
                            obs["alt_min"] = self.alt_min
                            obs["alt_max"] = self.alt_max
                            obs["HA_max"] = self.HA_max
                            obs["HA_min"] = self.HA_min

                            obs["note"] = self.survey_name + ", %i_t%i" % (
                                target_o_o.id,
                                time * 24,
                            )
                            obs_list.append(obs)
            observations = np.concatenate(obs_list)
            if self.obs_wanted is not None:
                if np.size(self.obs_wanted) > 0:
                    observations = np.concatenate([self.obs_wanted, observations])
            self.set_script(observations)

    def calc_reward_function(self, conditions):
        """If there is an observation ready to go, execute it, otherwise, -inf"""
        # check if any new event has come in

        if conditions.targets_of_opportunity is not None:
            for target_o_o in conditions.targets_of_opportunity:
                if target_o_o.id > self.last_event_id:
                    self._new_event(target_o_o, conditions)
                    self.last_event_id = target_o_o.id

        observation = self._check_list(conditions)
        if observation is None:
            self.reward = -np.inf
        else:
            self.reward = self.reward_val
        return self.reward




def generate_events(
    nside=32,
    mjd_start=59853.5,
    radius=6.5,
    survey_length=365.25 * 10,
    rate=10.0,
    expires=3.0,
    seed=42,
):
    """Generate a bunch of ToO events

    Parameters
    ----------
    rate : float (10)
        The number of events per year.
    expires : float (3)
        How long to keep broadcasting events as relevant (days)
    """

    np.random.seed(seed=seed)
    ra, dec = _hpid2_ra_dec(nside, np.arange(hp.nside2npix(nside)))
    print(ra, '  ' , dec)
    radius = np.radians(radius)
    print(radius)
    # Use a ceil here so we get at least 1 event even if doing a short run.
    n_events = int(np.ceil(survey_length / 365.25 * rate))
    names = ["mjd_start", "ra", "dec", "expires"]
    types = [float] * 4
    event_table = np.zeros(n_events, dtype=list(zip(names, types)))

    event_table["mjd_start"] = (
        np.sort(np.random.random(n_events)) * survey_length + mjd_start
    )
    event_table["expires"] = event_table["mjd_start"] + expires
    # Make sure latitude points spread correctly
    # http://mathworld.wolfram.com/SpherePointPicking.html
    event_table["ra"] = np.random.rand(n_events) * np.pi * 2
    event_table["dec"] = np.arccos(2.0 * np.random.rand(n_events) - 1.0) - np.pi / 2.0

    events = []
    for i, event_time in enumerate(event_table["mjd_start"]):
        dist = _angular_separation(ra, dec, event_table["ra"][i], event_table["dec"][i])
        good = np.where(dist <= radius)
        footprint = np.zeros(ra.size, dtype=float)
        footprint[good] = 1
        events.append(
            TargetoO(
                i,
                footprint,
                event_time,
                expires,
                ra_rad_center=event_table["ra"][i],
                dec_rad_center=event_table["dec"][i],
            )
        )
        print(ra.size)
        print(footprint)
        
    table = tabulate(event_table, names)
    print('ToO Events: \n',table)

    events = SimTargetooServer(events)
    return events, event_table

def incoming_event(
        #json_file=json_file,
    nside=32,
    mjd_start=59853.5,
    #radius=6.5,
    #survey_length=365.25 * 10,
    #rate=10.0,
    expires=3.0,
    #seed=42,
):
    """Gather Information from the Startegy Code .json file"""

    # Open the json file 
    with open(sys.argv[1], 'r') as f:
        data = json.load(f)

    expTime = []
    RA = []
    dec = []
    filters = []

    # No. of pointings
    n_point = len([ele for ele in data if isinstance(ele, dict)])
    
    for i in range(n_point):
        expTime.append(data[i]['expTime'])
        RA.append(data[i]['RA'])
        dec.append(data[i]['dec'])
        filters.append(data[i]['filter'])

    f.close()

    names = ["mjd_start", "ra", "dec", "expires"]
    types = [float] * 4
    pointing_table = np.zeros(n_point, dtype=list(zip(names, types)))

    pointing_table["mjd_start"] = np.full(n_point, mjd_start)
    pointing_table["expires"] = pointing_table["mjd_start"] + expires
    pointing_table["ra"] = np.array(RA)
    pointing_table["dec"] = np.array(dec)

    table = tabulate(pointing_table, names)
    print("Total Number of Pointings: ", pointing_table.size)
    print('ToO Events: \n',table)
    
    pointings = []
    for i, event_time in enumerate(pointing_table["mjd_start"]):
        footprint = np.array([1.0])
        pointings.append(
            TargetoO(
                i,
                footprint,
                event_time,
                expires,
                ra_rad_center=pointing_table["ra"][i],
                dec_rad_center=pointing_table["dec"][i],
            )
        )
    pointings = SimTargetooServer(pointings)
    
    return pointings, pointing_table, filters


def set_run_info(dbroot=None, file_end="v3.4_", out_dir="."):#, rate=None, ntoo=None):
    """Gather versions of software used to record"""
    extra_info = {}  # Dictionary to store executed command, git hash, file executed and scheduler git hash
    exec_command = ""
    for arg in sys.argv:
        exec_command += " " + arg
    extra_info["exec command"] = exec_command
    try:
        extra_info["git hash"] = subprocess.check_output(["git", "rev-parse", "HEAD"])
    except subprocess.CalledProcessError:
        extra_info["git hash"] = "Not in git repo"

    extra_info["file executed"] = os.path.realpath(__file__)
    try:
        rs_path = rubin_scheduler.__path__[0]
        hash_file = os.path.join(rs_path, "../", ".git/refs/heads/main")
        extra_info["rubin_scheduler git hash"] = subprocess.check_output(
            ["cat", hash_file]
        )
    except subprocess.CalledProcessError:
        pass

    # Use the filename of the script to name the output database
    if dbroot is None:
        fileroot = os.path.basename(sys.argv[0]).replace(".py", "") + "_"   # Use the filename of python script (too_example in this case)
    else:
        fileroot = dbroot 
    #fileroot += "%i_%i" % (rate, ntoo) # Append the too rate and nfollow to the db name
    fileroot = os.path.join(out_dir, fileroot + file_end) # Append the version to the filename
    return fileroot, extra_info




def run_sched(
    scheduler,
    observatory,
    survey_length=365.25,
    nside=32,
    fileroot="baseline_",
    verbose=False,
    extra_info=None,
    illum_limit=40.0,
    mjd_start=60796.0,
    event_table=None,
):
    """Run survey"""
    years = np.round(survey_length / 365.25)
    n_visit_limit = None
    fs = SimpleFilterSched(illum_limit=illum_limit)
    observatory, scheduler, observations = sim_runner(
        observatory,
        scheduler,
        survey_length=survey_length,
        filename=fileroot + "%iyrs.db" % years,
        delete_past=True,
        n_visit_limit=n_visit_limit,
        verbose=verbose,
        extra_info=extra_info,
        filter_scheduler=fs,
        event_table=event_table,
    )

    return observatory, scheduler, observations




def example_scheduler(args):
    survey_length = args.survey_length  # Days
    out_dir = args.out_dir
    verbose = args.verbose
    max_dither = args.maxDither
    illum_limit = args.moon_illum_limit
    nexp = args.nexp
    nslice = args.rolling_nslice
    rolling_scale = args.rolling_strength
    dbroot = args.dbroot
    nights_off = args.nights_off
    neo_night_pattern = args.neo_night_pattern
    neo_filters = args.neo_filters
    neo_repeat = args.neo_repeat
    ddf_season_frac = args.ddf_season_frac
    neo_am = args.neo_am
    neo_elong_req = args.neo_elong_req
    neo_area_req = args.neo_area_req
    nside = args.nside
    too_rate = args.too_rate
    too_filters = args.filters
    too_nfollow = args.nfollow
    #json_file = args.json_file

    mjd_start = args.mjd_start# 60796.0
    per_night = True  # Dither DDF per night

    camera_ddf_rot_limit = 75.0  # degrees
    
    fileroot, extra_info = set_run_info(
        dbroot=dbroot,
        file_end="v3.4_",
        out_dir=out_dir,
        #rate=too_rate,
        #ntoo=too_nfollow,
    )
    
    #  print('fileroot: ',fileroot,'\n Extra Info: ', extra_info)

    sim_ToOs, event_table, filters = incoming_event(nside=nside,mjd_start=mjd_start)

    too_filters = np.unique(np.array(filters))
    print("Filters: ", too_filters)
    #sim_ToOs, event_table = generate_events(
    #    nside=nside, survey_length=survey_length, rate=too_rate, mjd_start=mjd_start
    # )

    #print('Events: ',sim_ToOs,'\n Event Table: ', event_table)


    observatory = ModelObservatory(nside=nside, mjd_start=mjd_start, sim_to_o=sim_ToOs)   # sim_ToOs = events = SimTargetooServer(events)
    conditions = observatory.return_conditions()




    # ToO Survey ----------------------------------------------------------------------------------

    # Modify the footprint
    sky = EuclidOverlapFootprint(nside=nside, smc_radius=4, lmc_radius=6)
    footprints_hp_array, labels = sky.return_maps()

    wfd_indx = np.where(
        (labels == "lowdust") | (labels == "LMC_SMC") | (labels == "virgo")
    )[0]
    wfd_footprint = footprints_hp_array["r"] * 0
    wfd_footprint[wfd_indx] = 1

    footprints_hp = {}
    for key in footprints_hp_array.dtype.names:
        footprints_hp[key] = footprints_hp_array[key]

    footprint_mask = footprints_hp["r"] * 0
    footprint_mask[np.where(footprints_hp["r"] > 0)] = 1


    
    footprints = make_rolling_footprints(
        fp_hp=footprints_hp,
        mjd_start=conditions.mjd_start,
        sun_ra_start=conditions.sun_ra_start,
        nslice=nslice,
        scale=rolling_scale,
        nside=nside,
        wfd_indx=wfd_indx,
        order_roll=1,
        n_cycles=4,
    )

    
    # Let's make a footprint to follow up ToO events
    too_footprint = footprints_hp["r"] * 0 + np.nan
    too_footprint[np.where(footprints_hp["r"] > 0)[0]] = 1.0
 

    # Set up the DDF surveys to dither
    u_detailer = detailers.FilterNexp(filtername="u", nexp=1)
    dither_detailer = detailers.DitherDetailer(
        per_night=per_night, max_dither=max_dither
    )
    details = [
        detailers.CameraRotDetailer(
            min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit
        ),
        dither_detailer,
        u_detailer,
        detailers.Rottep2RotspDesiredDetailer(),
    ]
    euclid_detailers = [
        detailers.CameraRotDetailer(
            min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit
        ),
        detailers.EuclidDitherDetailer(),
        u_detailer,
        detailers.Rottep2RotspDesiredDetailer(),
    ]

    
    # Set up the damn ToO kwargs
    times = [0, 1, 2, 4, 24]
    filters_at_times = [too_filters] * 5
    nvis = [1, 1, 1, 1, 6]
    toos = [
        ToO_scripted_survey(
            [],
            nside=nside,
            followup_footprint=too_footprint,
            times=times[0:too_nfollow],
            filters_at_times=filters_at_times[0:too_nfollow],
            nvis=nvis[0:too_nfollow],
            detailers=details,
        )
    ]

    surveys = [toos]#, ddfs, long_gaps, blobs, twi_blobs, neo, greedy]

    scheduler = CoreScheduler(surveys, nside=nside)

    if args.setup_only:
        return scheduler
    else:
        observatory, scheduler, observations = run_sched(
            scheduler,
            observatory,
            survey_length=survey_length,
            verbose=verbose,
            fileroot=os.path.join(out_dir, fileroot),
            extra_info=extra_info,
            nside=nside,
            illum_limit=illum_limit,
            mjd_start=mjd_start,
            event_table=event_table,
        )
        return observatory, scheduler, observations



def sched_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", dest="verbose", action="store_true")
    parser.set_defaults(verbose=False)
    parser.add_argument("--survey_length", type=float, default= 30)#365.25 * 10)
    parser.add_argument("--out_dir", type=str, default="")
    parser.add_argument(
        "--maxDither", type=float, default=0.7, help="Dither size for DDFs (deg)"
    )
    parser.add_argument(
        "--moon_illum_limit",
        type=float,
        default=40.0,
        help="illumination limit to remove u-band",
    )
    parser.add_argument("--nexp", type=int, default=2)
    parser.add_argument("--rolling_nslice", type=int, default=2)
    parser.add_argument("--rolling_strength", type=float, default=0.9)
    parser.add_argument("--dbroot", type=str)
    parser.add_argument("--ddf_season_frac", type=float, default=0.2)
    parser.add_argument("--nights_off", type=int, default=3, help="For long gaps")
    parser.add_argument("--neo_night_pattern", type=int, default=4)
    parser.add_argument("--neo_filters", type=str, default="riz")
    parser.add_argument("--neo_repeat", type=int, default=4)
    parser.add_argument("--neo_am", type=float, default=2.5)
    parser.add_argument("--neo_elong_req", type=float, default=45.0)
    parser.add_argument("--neo_area_req", type=float, default=0.0)
    parser.add_argument(
        "--setup_only", dest="setup_only", default=False, action="store_true"
    )
    parser.add_argument(
        "--nside",
        type=int,
        default=32,
        help="Nside should be set to default (32) except for tests.",
    )
    parser.add_argument("--too_rate", type=float, default=10, help="N events per year")
    parser.add_argument("--filters", type=str, default="z")
    parser.add_argument("--nfollow", type=int, default=1)
    parser.add_argument("--mjd_start", type=float, default=60796.0)


    parser.add_argument("json", type = argparse.FileType('r'))
    return parser


if __name__ == "__main__":
    parser = sched_argparser()
    args = parser.parse_args()
    example_scheduler(args)
