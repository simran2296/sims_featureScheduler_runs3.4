#!/usr/bin/env python

import argparse
import os
import subprocess
import sys

import healpy as hp
import matplotlib.pylab as plt
import numpy as np
import rubin_scheduler
import rubin_scheduler.scheduler.basis_functions as bf
import rubin_scheduler.scheduler.detailers as detailers
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils import iers
from rubin_scheduler.site_models import _read_fields
from rubin_scheduler.scheduler import sim_runner
from rubin_scheduler.scheduler.model_observatory import ModelObservatory
from rubin_scheduler.scheduler.schedulers import CoreScheduler, SimpleFilterSched
from rubin_scheduler.scheduler.surveys import (
    BlobSurvey,
    GreedySurvey,
    LongGapSurvey,
    ScriptedSurvey,
    generate_ddf_scheduled_obs,
    BaseMarkovSurvey,
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
from rubin_scheduler.site_models import Almanac
from rubin_scheduler.utils import (
    _hpid2_ra_dec,
    _ra_dec2_hpid,
    _approx_ra_dec2_alt_az,
    _angular_separation,
)

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
        times=[1, 2, 4, 24, 48],
        filters_at_times=["gz", "gz", "gz", "gz", "gz", "gz"],
        nvis=[1, 1, 1, 1, 6, 6],
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
        self.last_mjd = -1

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

        # Don't bother with checking if we can run before twilight ends
        self.before_twi_check = False

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
            # Put the fields in a good order.
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

    def generate_observations(self, conditions):
        observations = self.generate_observations_rough(conditions)

        if len(observations) > 0:
            for detailer in self.detailers:
                observations = detailer(observations, conditions)

        return observations


def mean_longitude(longitude):
    """Compute a mean longitude, accounting for wrap around."""
    x = np.cos(longitude)
    y = np.sin(longitude)
    meanx = np.mean(x)
    meany = np.mean(y)
    angle = np.arctan2(meany, meanx)
    radius = np.sqrt(meanx**2 + meany**2)
    mid_longitude = angle % (2.0 * np.pi)
    if radius < 0.1:
        mid_longitude = np.pi
    return mid_longitude


# XXX--should move this to utils generally and refactor it out of base MarkovSurvey
def order_observations(RA, dec):
    """
    Take a list of ra,dec positions and compute a traveling salesman solution through them
    """
    # Let's find a good spot to project the points to a plane
    mid_dec = (np.max(dec) - np.min(dec)) / 2.0
    mid_ra = mean_longitude(RA)
    # Project the coordinates to a plane. Could consider scaling things to represent
    # time between points rather than angular distance.
    pointing_x, pointing_y = gnomonic_project_toxy(RA, dec, mid_ra, mid_dec)
    # Round off positions so that we ensure identical cross-platform performance
    #scale = 1e6
    #pointing_x = np.round(pointing_x * scale).astype(int)
    #pointing_y = np.round(pointing_y * scale).astype(int)
    # Now I have a bunch of x,y pointings. Drop into TSP solver to get an effiencent route
    towns = np.vstack((pointing_x, pointing_y)).T
    # Leaving optimize=False for speed. The optimization step doesn't usually improve much.
    better_order = tsp_convex(towns, optimize=False)
    return better_order


def generate_events(
    nside=32,
    mjd_start=59853.5,
    radius_possible=[2.52, 4.0, 5.64],
    radius_fractions=[0.15, 0.5, 1.],
    survey_length=365.25 * 2,
    rate=10.0,
    expires=3.0,
    seed=42,
    frac_inter=0.8,
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
    # Use a ceil here so we get at least 1 event even if doing a short run.
    n_events = int(np.ceil(survey_length / 365.25 * rate))

    names = ["mjd_start", "ra", "dec", "expires", "radius"]
    types = [float] * 5
    event_table = np.zeros(n_events, dtype=list(zip(names, types)))

    event_table["mjd_start"] = (
        np.sort(np.random.random(n_events)) * survey_length + mjd_start
    )
    
    # Make sure latitude points spread correctly
    # http://mathworld.wolfram.com/SpherePointPicking.html
    event_table["ra"] = np.random.rand(n_events) * np.pi * 2
    event_table["dec"] = np.arccos(2.0 * np.random.rand(n_events) - 1.0) - np.pi / 2.0

    # set a searcha radius for each event
    draw = np.random.uniform(size=n_events)

    indx = np.searchsorted(radius_fractions, draw)

    radius = np.radians(radius_possible)[indx]

    event_table["radius"] = radius

    draw = np.random.uniform(size=n_events)
    no_intur = draw > frac_inter

    # for now, use flat 100 hours
    intur_time = np.random.uniform(size=n_events)*100/24. + event_table["mjd_start"]

    # If we are not inturupting, use expires time
    intur_time[np.where(no_intur)[0]] = expires + event_table["mjd_start"][np.where(no_intur)]

    event_table["expires"] = intur_time

    events = []
    for i, event_time in enumerate(event_table["mjd_start"]):
        dist = _angular_separation(ra, dec, event_table["ra"][i], event_table["dec"][i])
        good = np.where(dist <= radius[i])
        footprint = np.zeros(ra.size, dtype=float)
        footprint[good] = 1
        events.append(
            TargetoO(
                i,
                footprint,
                event_time,
                intur_time[i],
                ra_rad_center=event_table["ra"][i],
                dec_rad_center=event_table["dec"][i],
            )
        )
    events = SimTargetooServer(events)
    return events, event_table


def standard_bf(
    nside,
    filtername="g",
    filtername2="i",
    m5_weight=6.0,
    footprint_weight=1.5,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    template_weight=12.0,
    u_template_weight=50.0,
    g_template_weight=50.0,
    footprints=None,
    n_obs_template=None,
    season=300.0,
    season_start_hour=-4.0,
    season_end_hour=2.0,
    moon_distance=30.0,
    strict=True,
    wind_speed_maximum=20.0,
):
    """Generate the standard basis functions that are shared by blob surveys

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filtername : list of str
        The filternames for the first set
    filtername2 : list of str
        The filter names for the second in the pair (None if unpaired)
    n_obs_template : dict (None)
        The number of observations to take every season in each filter
    season : float (300)
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : float (-4.)
        For weighting how strongly a template image needs to be observed (hours)
    sesason_end_hour : float (2.)
        For weighting how strongly a template image needs to be observed (hours)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    template_weight : float (12.)
        The weight to place on getting image templates every season
    u_template_weight : float (24.)
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a little higher than
        the standard template_weight kwarg.
    g_template_weight : float (24.)
        The weight to place on getting image templates in g-band. Since there
        are so few g-visits, it can be helpful to turn this up a little higher than
        the standard template_weight kwarg.

    Returns
    -------
    basis_functions_weights : `list`
        list of tuple pairs (basis function, weight) that is
        (rubin_scheduler.scheduler.BasisFunction object, float)

    """
    template_weights = {
        "u": u_template_weight,
        "g": g_template_weight,
        "r": template_weight,
        "i": template_weight,
        "z": template_weight,
        "y": template_weight,
    }

    bfs = []

    if filtername2 is not None:
        bfs.append(
            (
                bf.M5DiffBasisFunction(filtername=filtername, nside=nside),
                m5_weight / 2.0,
            )
        )
        bfs.append(
            (
                bf.M5DiffBasisFunction(filtername=filtername2, nside=nside),
                m5_weight / 2.0,
            )
        )

    else:
        bfs.append(
            (bf.M5DiffBasisFunction(filtername=filtername, nside=nside), m5_weight)
        )

    if filtername2 is not None:
        bfs.append(
            (
                bf.FootprintBasisFunction(
                    filtername=filtername,
                    footprint=footprints,
                    out_of_bounds_val=np.nan,
                    nside=nside,
                ),
                footprint_weight / 2.0,
            )
        )
        bfs.append(
            (
                bf.FootprintBasisFunction(
                    filtername=filtername2,
                    footprint=footprints,
                    out_of_bounds_val=np.nan,
                    nside=nside,
                ),
                footprint_weight / 2.0,
            )
        )
    else:
        bfs.append(
            (
                bf.FootprintBasisFunction(
                    filtername=filtername,
                    footprint=footprints,
                    out_of_bounds_val=np.nan,
                    nside=nside,
                ),
                footprint_weight,
            )
        )

    bfs.append(
        (
            bf.SlewtimeBasisFunction(filtername=filtername, nside=nside),
            slewtime_weight,
        )
    )
    if strict:
        bfs.append(
            (bf.StrictFilterBasisFunction(filtername=filtername), stayfilter_weight)
        )
    else:
        bfs.append(
            (bf.FilterChangeBasisFunction(filtername=filtername), stayfilter_weight)
        )

    if n_obs_template is not None:
        if filtername2 is not None:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername),
                        n_obs=n_obs_template[filtername],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weights[filtername] / 2.0,
                )
            )
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername2,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername2),
                        n_obs=n_obs_template[filtername2],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weights[filtername2] / 2.0,
                )
            )
        else:
            bfs.append(
                (
                    bf.NObsPerYearBasisFunction(
                        filtername=filtername,
                        nside=nside,
                        footprint=footprints.get_footprint(filtername),
                        n_obs=n_obs_template[filtername],
                        season=season,
                        season_start_hour=season_start_hour,
                        season_end_hour=season_end_hour,
                    ),
                    template_weights[filtername],
                )
            )

    # The shared masks
    bfs.append(
        (
            bf.MoonAvoidanceBasisFunction(nside=nside, moon_distance=moon_distance),
            0.0,
        )
    )
    bfs.append((bf.AvoidDirectWind(nside=nside, wind_speed_maximum=wind_speed_maximum), 0))
    filternames = [fn for fn in [filtername, filtername2] if fn is not None]
    bfs.append((bf.FilterLoadedBasisFunction(filternames=filternames), 0))
    bfs.append((bf.PlanetMaskBasisFunction(nside=nside), 0.0))

    return bfs


def blob_for_long(
    nside,
    nexp=2,
    exptime=30.0,
    filter1s=["g"],
    filter2s=["i"],
    pair_time=33.0,
    camera_rot_limits=[-80.0, 80.0],
    n_obs_template=None,
    season=300.0,
    season_start_hour=-4.0,
    season_end_hour=2.0,
    shadow_minutes=60.0,
    max_alt=76.0,
    moon_distance=30.0,
    ignore_obs=["DD", "twilight_near_sun"],
    m5_weight=6.0,
    footprint_weight=1.5,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    template_weight=12.0,
    u_template_weight=50.0,
    g_template_weight=50.0,
    footprints=None,
    u_nexp1=True,
    night_pattern=[True, True],
    time_after_twi=30.0,
    HA_min=12,
    HA_max=24 - 3.5,
    blob_names=[],
):
    """
    Generate surveys that take observations in blobs.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filter1s : list of str
        The filternames for the first set
    filter2s : list of str
        The filter names for the second in the pair (None if unpaired)
    pair_time : float (33)
        The ideal time between pairs (minutes)
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    n_obs_template : dict (None)
        The number of observations to take every season in each filter. If None, sets to 3 each.
    season : float (300)
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : float (-4.)
        For weighting how strongly a template image needs to be observed (hours)
    sesason_end_hour : float (2.)
        For weighting how strongly a template image needs to be observed (hours)
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    template_weight : float (12.)
        The weight to place on getting image templates every season
    u_template_weight : float (24.)
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a little higher than
        the standard template_weight kwarg.
    u_nexp1 : bool (True)
        Add a detailer to make sure the number of expossures in a visit is always 1 for u observations.
    """

    BlobSurvey_params = {
        "slew_approx": 7.5,
        "filter_change_approx": 140.0,
        "read_approx": 2.0,
        "min_pair_time": 15.0,
        "search_radius": 30.0,
        "alt_max": 85.0,
        "az_range": None,
        "flush_time": 30.0,
        "smoothing_kernel": None,
        "nside": nside,
        "seed": 42,
        "dither": True,
        "twilight_scale": True,
    }

    surveys = []
    if n_obs_template is None:
        n_obs_template = {"u": 3, "g": 3, "r": 3, "i": 3, "z": 3, "y": 3}

    times_needed = [pair_time, pair_time * 2]
    for filtername, filtername2 in zip(filter1s, filter2s):
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                filtername=filtername,
                filtername2=filtername2,
                m5_weight=m5_weight,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayfilter_weight=stayfilter_weight,
                template_weight=template_weight,
                u_template_weight=u_template_weight,
                g_template_weight=g_template_weight,
                footprints=footprints,
                n_obs_template=n_obs_template,
                season=season,
                season_start_hour=season_start_hour,
                season_end_hour=season_end_hour,
            )
        )

        # Masks, give these 0 weight
        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(
                    nside=nside,
                    shadow_minutes=shadow_minutes,
                    max_alt=max_alt,
                ),
                0.0,
            )
        )
        if filtername2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=time_needed), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))
        bfs.append((bf.AfterEveningTwiBasisFunction(time_after=time_after_twi), 0.0))
        bfs.append((bf.HaMaskBasisFunction(ha_min=HA_min, ha_max=HA_max), 0.0))
        # don't execute every night
        bfs.append((bf.NightModuloBasisFunction(night_pattern), 0.0))
        # only execute one blob per night
        bfs.append((bf.OnceInNightBasisFunction(notes=blob_names), 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if filtername2 is None:
            survey_name = "blob_long, %s" % filtername
        else:
            survey_name = "blob_long, %s%s" % (filtername, filtername2)
        if filtername2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(filtername=filtername2))

        if u_nexp1:
            detailer_list.append(detailers.FilterNexp(filtername="u", nexp=1))
        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                filtername1=filtername,
                filtername2=filtername2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_note=survey_name,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **BlobSurvey_params
            )
        )

    return surveys


def gen_long_gaps_survey(
    footprints,
    nside=32,
    night_pattern=[True, True],
    gap_range=[2, 7],
    HA_min=12,
    HA_max=24 - 3.5,
    time_after_twi=120,
    u_template_weight=50.0,
    g_template_weight=50.0,
):
    """
    Paramterers
    -----------
    HA_min(_max) : float
        The hour angle limits passed to the initial blob scheduler.
    """

    surveys = []
    f1 = ["g", "r", "i"]
    f2 = ["r", "i", "z"]
    # Maybe force scripted to not go in twilight?
    blob_names = []
    for fn1, fn2 in zip(f1, f2):
        for ab in ["a", "b"]:
            blob_names.append("blob_long, %s%s, %s" % (fn1, fn2, ab))
    for filtername1, filtername2 in zip(f1, f2):
        blob = blob_for_long(
            footprints=footprints,
            nside=nside,
            filter1s=[filtername1],
            filter2s=[filtername2],
            night_pattern=night_pattern,
            time_after_twi=time_after_twi,
            HA_min=HA_min,
            HA_max=HA_max,
            u_template_weight=u_template_weight,
            g_template_weight=g_template_weight,
            blob_names=blob_names,
        )
        scripted = ScriptedSurvey(
            [bf.AvoidDirectWind(nside=nside)], nside=nside, ignore_obs=["blob", "DDF", "twi", "pair"]
        )
        surveys.append(
            LongGapSurvey(blob[0], scripted, gap_range=gap_range, avoid_zenith=True)
        )

    return surveys


def gen_greedy_surveys(
    nside=32,
    nexp=2,
    exptime=30.0,
    filters=["r", "i", "z", "y"],
    camera_rot_limits=[-80.0, 80.0],
    shadow_minutes=60.0,
    max_alt=76.0,
    moon_distance=30.0,
    ignore_obs=["DD", "twilight_near_sun"],
    m5_weight=3.0,
    footprint_weight=0.75,
    slewtime_weight=3.0,
    stayfilter_weight=100.0,
    repeat_weight=-1.0,
    footprints=None,
):
    """
    Make a quick set of greedy surveys

    This is a convenience function to generate a list of survey objects that can be used with
    rubin_scheduler.scheduler.schedulers.Core_scheduler.
    To ensure we are robust against changes in the sims_featureScheduler codebase, all kwargs are
    explicitly set.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filters : list of str (['r', 'i', 'z', 'y'])
        Which filters to generate surveys for.
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    """
    # Define the extra parameters that are used in the greedy survey. I
    # think these are fairly set, so no need to promote to utility func kwargs
    greed_survey_params = {
        "block_size": 1,
        "smoothing_kernel": None,
        "seed": 42,
        "camera": "LSST",
        "dither": True,
        "survey_name": "greedy",
    }

    surveys = []
    detailer_list = [
        detailers.CameraRotDetailer(
            min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
        )
    ]
    detailer_list.append(detailers.Rottep2RotspDesiredDetailer())

    for filtername in filters:
        bfs = []
        bfs.extend(
            standard_bf(
                nside,
                filtername=filtername,
                filtername2=None,
                m5_weight=m5_weight,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayfilter_weight=stayfilter_weight,
                template_weight=0,
                u_template_weight=0,
                g_template_weight=0,
                footprints=footprints,
                n_obs_template=None,
                strict=False,
            )
        )

        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=2 * 60.0, filtername=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )
        # Masks, give these 0 weight
        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(nside=nside, shadow_minutes=shadow_minutes, max_alt=max_alt),
                0,
            )
        )
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        surveys.append(
            GreedySurvey(
                basis_functions,
                weights,
                exptime=exptime,
                filtername=filtername,
                nside=nside,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **greed_survey_params
            )
        )

    return surveys


def generate_blobs(
    nside,
    nexp=2,
    exptime=30.0,
    filter1s=["u", "u", "g", "r", "i", "z", "y"],
    filter2s=["g", "r", "r", "i", "z", "y", "y"],
    pair_time=33.0,
    camera_rot_limits=[-80.0, 80.0],
    n_obs_template=None,
    season=300.0,
    season_start_hour=-4.0,
    season_end_hour=2.0,
    shadow_minutes=60.0,
    max_alt=76.0,
    moon_distance=30.0,
    ignore_obs=["DD", "twilight_near_sun"],
    m5_weight=6.0,
    footprint_weight=1.5,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    template_weight=12.0,
    u_template_weight=50.0,
    g_template_weight=50.0,
    footprints=None,
    u_nexp1=True,
    scheduled_respect=45.0,
    good_seeing={"g": 3, "r": 3, "i": 3},
    good_seeing_weight=3.0,
    mjd_start=1,
    repeat_weight=-20,
):
    """
    Generate surveys that take observations in blobs.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filter1s : list of str
        The filternames for the first set
    filter2s : list of str
        The filter names for the second in the pair (None if unpaired)
    pair_time : float (33)
        The ideal time between pairs (minutes)
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    n_obs_template : Dict (None)
        The number of observations to take every season in each filter. If None, sets to 3 each.
    season : float (300)
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : float (-4.)
        For weighting how strongly a template image needs to be observed (hours)
    sesason_end_hour : float (2.)
        For weighting how strongly a template image needs to be observed (hours)
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    template_weight : float (12.)
        The weight to place on getting image templates every season
    u_template_weight : float (24.)
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a little higher than
        the standard template_weight kwarg.
    u_nexp1 : bool (True)
        Add a detailer to make sure the number of expossures in a visit is always 1 for u observations.
    scheduled_respect : float (45)
        How much time to require there be before a pre-scheduled observation (minutes)
    """

    BlobSurvey_params = {
        "slew_approx": 7.5,
        "filter_change_approx": 140.0,
        "read_approx": 2.0,
        "min_pair_time": 15.0,
        "search_radius": 30.0,
        "alt_max": 85.0,
        "az_range": None,
        "flush_time": 30.0,
        "smoothing_kernel": None,
        "nside": nside,
        "seed": 42,
        "dither": True,
        "twilight_scale": False,
    }

    if n_obs_template is None:
        n_obs_template = {"u": 3, "g": 3, "r": 3, "i": 3, "z": 3, "y": 3}

    surveys = []

    times_needed = [pair_time, pair_time * 2]
    for filtername, filtername2 in zip(filter1s, filter2s):
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.Rottep2RotspDesiredDetailer())
        detailer_list.append(detailers.CloseAltDetailer())
        detailer_list.append(detailers.FlushForSchedDetailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                filtername=filtername,
                filtername2=filtername2,
                m5_weight=m5_weight,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayfilter_weight=stayfilter_weight,
                template_weight=template_weight,
                u_template_weight=u_template_weight,
                g_template_weight=g_template_weight,
                footprints=footprints,
                n_obs_template=n_obs_template,
                season=season,
                season_start_hour=season_start_hour,
                season_end_hour=season_end_hour,
            )
        )

        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=3 * 60.0, filtername=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )

        # Insert things for getting good seeing templates
        if filtername2 is not None:
            if filtername in list(good_seeing.keys()):
                bfs.append(
                    (
                        bf.NGoodSeeingBasisFunction(
                            filtername=filtername,
                            nside=nside,
                            mjd_start=mjd_start,
                            footprint=footprints.get_footprint(filtername),
                            n_obs_desired=good_seeing[filtername],
                        ),
                        good_seeing_weight,
                    )
                )
            if filtername2 in list(good_seeing.keys()):
                bfs.append(
                    (
                        bf.NGoodSeeingBasisFunction(
                            filtername=filtername2,
                            nside=nside,
                            mjd_start=mjd_start,
                            footprint=footprints.get_footprint(filtername2),
                            n_obs_desired=good_seeing[filtername2],
                        ),
                        good_seeing_weight,
                    )
                )
        else:
            if filtername in list(good_seeing.keys()):
                bfs.append(
                    (
                        bf.NGoodSeeingBasisFunction(
                            filtername=filtername,
                            nside=nside,
                            mjd_start=mjd_start,
                            footprint=footprints.get_footprint(filtername),
                            n_obs_desired=good_seeing[filtername],
                        ),
                        good_seeing_weight,
                    )
                )
        # Make sure we respect scheduled observations
        bfs.append((bf.TimeToScheduledBasisFunction(time_needed=scheduled_respect), 0))
        # Masks, give these 0 weight
        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(
                    nside=nside,
                    shadow_minutes=shadow_minutes,
                    max_alt=max_alt,
                ),
                0.0,
            )
        )
        if filtername2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append((bf.TimeToTwilightBasisFunction(time_needed=time_needed), 0.0))
        bfs.append((bf.NotTwilightBasisFunction(), 0.0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if filtername2 is None:
            survey_name = "pair_%i, %s" % (pair_time, filtername)
        else:
            survey_name = "pair_%i, %s%s" % (pair_time, filtername, filtername2)
        if filtername2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(filtername=filtername2))

        if u_nexp1:
            detailer_list.append(detailers.FilterNexp(filtername="u", nexp=1))
        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                filtername1=filtername,
                filtername2=filtername2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_note=survey_name,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **BlobSurvey_params
            )
        )

    return surveys


def generate_twi_blobs(
    nside,
    nexp=2,
    exptime=30.0,
    filter1s=["r", "i", "z", "y"],
    filter2s=["i", "z", "y", "y"],
    pair_time=15.0,
    camera_rot_limits=[-80.0, 80.0],
    n_obs_template=None,
    season=300.0,
    season_start_hour=-4.0,
    season_end_hour=2.0,
    shadow_minutes=60.0,
    max_alt=76.0,
    moon_distance=30.0,
    ignore_obs=["DD", "twilight_near_sun"],
    m5_weight=6.0,
    footprint_weight=1.5,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    template_weight=12.0,
    footprints=None,
    repeat_night_weight=None,
    wfd_footprint=None,
    scheduled_respect=15.0,
    repeat_weight=-1.0,
    night_pattern=None,
):
    """
    Generate surveys that take observations in blobs.

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    nexp : int (1)
        The number of exposures to use in a visit.
    exptime : float (30.)
        The exposure time to use per visit (seconds)
    filter1s : list of str
        The filternames for the first set
    filter2s : list of str
        The filter names for the second in the pair (None if unpaired)
    pair_time : float (22)
        The ideal time between pairs (minutes)
    camera_rot_limits : list of float ([-80., 80.])
        The limits to impose when rotationally dithering the camera (degrees).
    n_obs_template : dict (None)
        The number of observations to take every season in each filter. If None, sets to 3 each.
    season : float (300)
        The length of season (i.e., how long before templates expire) (days)
    season_start_hour : float (-4.)
        For weighting how strongly a template image needs to be observed (hours)
    sesason_end_hour : float (2.)
        For weighting how strongly a template image needs to be observed (hours)
    shadow_minutes : float (60.)
        Used to mask regions around zenith (minutes)
    max_alt : float (76.
        The maximium altitude to use when masking zenith (degrees)
    moon_distance : float (30.)
        The mask radius to apply around the moon (degrees)
    ignore_obs : str or list of str ('DD')
        Ignore observations by surveys that include the given substring(s).
    m5_weight : float (3.)
        The weight for the 5-sigma depth difference basis function
    footprint_weight : float (0.3)
        The weight on the survey footprint basis function.
    slewtime_weight : float (3.)
        The weight on the slewtime basis function
    stayfilter_weight : float (3.)
        The weight on basis function that tries to stay avoid filter changes.
    template_weight : float (12.)
        The weight to place on getting image templates every season
    u_template_weight : float (24.)
        The weight to place on getting image templates in u-band. Since there
        are so few u-visits, it can be helpful to turn this up a little higher than
        the standard template_weight kwarg.
    """

    BlobSurvey_params = {
        "slew_approx": 7.5,
        "filter_change_approx": 140.0,
        "read_approx": 2.0,
        "min_pair_time": 10.0,
        "search_radius": 30.0,
        "alt_max": 85.0,
        "az_range": None,
        "flush_time": 30.0,
        "smoothing_kernel": None,
        "nside": nside,
        "seed": 42,
        "dither": True,
        "twilight_scale": False,
        "in_twilight": True,
    }

    surveys = []

    if n_obs_template is None:
        n_obs_template = {"u": 3, "g": 3, "r": 3, "i": 3, "z": 3, "y": 3}

    times_needed = [pair_time, pair_time * 2]
    for filtername, filtername2 in zip(filter1s, filter2s):
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.Rottep2RotspDesiredDetailer())
        detailer_list.append(detailers.CloseAltDetailer())
        detailer_list.append(detailers.FlushForSchedDetailer())
        # List to hold tuples of (basis_function_object, weight)
        bfs = []

        bfs.extend(
            standard_bf(
                nside,
                filtername=filtername,
                filtername2=filtername2,
                m5_weight=m5_weight,
                footprint_weight=footprint_weight,
                slewtime_weight=slewtime_weight,
                stayfilter_weight=stayfilter_weight,
                template_weight=template_weight,
                u_template_weight=0,
                g_template_weight=0,
                footprints=footprints,
                n_obs_template=n_obs_template,
                season=season,
                season_start_hour=season_start_hour,
                season_end_hour=season_end_hour,
            )
        )

        bfs.append(
            (
                bf.VisitRepeatBasisFunction(
                    gap_min=0, gap_max=2 * 60.0, filtername=None, nside=nside, npairs=20
                ),
                repeat_weight,
            )
        )

        if repeat_night_weight is not None:
            bfs.append(
                (
                    bf.AvoidLongGapsBasisFunction(
                        nside=nside,
                        filtername=None,
                        min_gap=0.0,
                        max_gap=10.0 / 24.0,
                        ha_limit=3.5,
                        footprint=wfd_footprint,
                    ),
                    repeat_night_weight,
                )
            )
        # Make sure we respect scheduled observations
        bfs.append((bf.TimeToScheduledBasisFunction(time_needed=scheduled_respect), 0))
        # Masks, give these 0 weight
        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(
                    nside=nside,
                    shadow_minutes=shadow_minutes,
                    max_alt=max_alt,
                ),
                0.0,
            )
        )
        if filtername2 is None:
            time_needed = times_needed[0]
        else:
            time_needed = times_needed[1]
        bfs.append(
            (bf.TimeToTwilightBasisFunction(time_needed=time_needed, alt_limit=12), 0.0)
        )

        # Let's turn off twilight blobs on nights where we are
        # doing NEO hunts
        bfs.append((bf.NightModuloBasisFunction(pattern=night_pattern), 0))

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]
        if filtername2 is None:
            survey_name = "pair_%i, %s" % (pair_time, filtername)
        else:
            survey_name = "pair_%i, %s%s" % (pair_time, filtername, filtername2)
        if filtername2 is not None:
            detailer_list.append(detailers.TakeAsPairsDetailer(filtername=filtername2))
        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                filtername1=filtername,
                filtername2=filtername2,
                exptime=exptime,
                ideal_pair_time=pair_time,
                survey_note=survey_name,
                ignore_obs=ignore_obs,
                nexp=nexp,
                detailers=detailer_list,
                **BlobSurvey_params
            )
        )

    return surveys


def ddf_surveys(detailers=None, season_unobs_frac=0.2, euclid_detailers=None, nside=None):
    obs_array = generate_ddf_scheduled_obs(season_unobs_frac=season_unobs_frac)

    euclid_obs = np.where(
        (obs_array["note"] == "DD:EDFS_b") | (obs_array["note"] == "DD:EDFS_a")
    )[0]
    all_other = np.where(
        (obs_array["note"] != "DD:EDFS_b") & (obs_array["note"] != "DD:EDFS_a")
    )[0]

    survey1 = ScriptedSurvey([bf.AvoidDirectWind(nside=nside)], detailers=detailers)
    survey1.set_script(obs_array[all_other])

    survey2 = ScriptedSurvey([bf.AvoidDirectWind(nside=nside)], detailers=euclid_detailers)
    survey2.set_script(obs_array[euclid_obs])

    return [survey1, survey2]


def ecliptic_target(nside=32, dist_to_eclip=40.0, dec_max=30.0, mask=None):
    """Generate a target_map for the area around the ecliptic

    Parameters
    ----------
    nside : int (32)
        The HEALpix nside to use
    dist_to_eclip : float (40)
        The distance to the ecliptic to constrain to (degrees).
    dec_max : float (30)
        The max declination to alow (degrees).
    mask : np.array (None)
        Any additional mask to apply, should be a HEALpix mask with matching nside.
    """

    ra, dec = _hpid2_ra_dec(nside, np.arange(hp.nside2npix(nside)))
    result = np.zeros(ra.size)
    coord = SkyCoord(ra=ra * u.rad, dec=dec * u.rad)
    eclip_lat = coord.barycentrictrueecliptic.lat.radian
    good = np.where(
        (np.abs(eclip_lat) < np.radians(dist_to_eclip)) & (dec < np.radians(dec_max))
    )
    result[good] += 1

    if mask is not None:
        result *= mask

    return result


def generate_twilight_near_sun(
    nside,
    night_pattern=None,
    nexp=1,
    exptime=15,
    ideal_pair_time=5.0,
    max_airmass=2.0,
    camera_rot_limits=[-80.0, 80.0],
    time_needed=10,
    footprint_mask=None,
    footprint_weight=0.1,
    slewtime_weight=3.0,
    stayfilter_weight=3.0,
    area_required=None,
    filters="riz",
    n_repeat=4,
    sun_alt_limit=-14.8,
    slew_estimate=4.5,
    moon_distance=30.0,
    shadow_minutes=0,
    min_alt=20.,
    max_alt=76.0,
    max_elong=60.0,
    az_range=180.0,
    ignore_obs=["DD", "pair", "long", "blob", "greedy"],
    filter_dist_weight=0.3,
    time_to_12deg=25.0,
):
    """Generate a survey for observing NEO objects in twilight

    Parameters
    ----------
    night_pattern : list of bool (None)
        A list of bools that set when the survey will be active. e.g., [True, False]
        for every-other night, [True, False, False] for every third night.
    nexp : int (1)
        Number of snaps in a visit
    exptime : float (15)
        Exposure time of visits
    ideal_pair_time : float (5)
        Ideal time between repeat visits (minutes).
    max_airmass : float (2)
        Maximum airmass to attempt (unitless).
    camera_rot_limits : list of float ([-80, 80])
        The camera rotation limits to use (degrees).
    time_needed : float (10)
        How much time should be available (e.g., before twilight ends) (minutes).
    footprint_mask : np.array (None)
        Mask to apply to the constructed ecliptic target mask (None).
    footprint_weight : float (0.1)
        Weight for footprint basis function
    slewtime_weight : float (3.)
        Weight for slewtime basis function
    stayfilter_weight : float (3.)
        Weight for staying in the same filter basis function
    area_required : float (None)
        The area that needs to be available before the survey will return observations (sq degrees?)
    filters : str ('riz')
        The filters to use, default 'riz'
    n_repeat : int (4)
        The number of times a blob should be repeated, default 4.
    sun_alt_limit : float (-14.8)
        Do not start unless sun is higher than this limit (degrees)
    slew_estimate : float (4.5)
        An estimate of how long it takes to slew between neighboring fields (seconds).
    time_to_sunrise : float (25.)
        Do not execute if time to sunrise is greater than (minutes).
    """
    survey_name = "twilight_near_sun"
    footprint = ecliptic_target(nside=nside, mask=footprint_mask)
    constant_fp = ConstantFootprint(nside=nside)
    for filtername in filters:
        constant_fp.set_footprint(filtername, footprint)

    surveys = []
    for filtername in filters:
        detailer_list = []
        detailer_list.append(
            detailers.CameraRotDetailer(
                min_rot=np.min(camera_rot_limits), max_rot=np.max(camera_rot_limits)
            )
        )
        detailer_list.append(detailers.CloseAltDetailer())
        # Should put in a detailer so things start at lowest altitude
        detailer_list.append(
            detailers.TwilightTripleDetailer(
                slew_estimate=slew_estimate, n_repeat=n_repeat
            )
        )
        bfs = []

        bfs.append(
            (
                bf.FootprintBasisFunction(
                    filtername=filtername,
                    footprint=constant_fp,
                    out_of_bounds_val=np.nan,
                    nside=nside,
                ),
                footprint_weight,
            )
        )

        bfs.append(
            (
                bf.SlewtimeBasisFunction(filtername=filtername, nside=nside),
                slewtime_weight,
            )
        )
        bfs.append(
            (bf.StrictFilterBasisFunction(filtername=filtername), stayfilter_weight)
        )
        bfs.append(
            (bf.FilterDistBasisFunction(filtername=filtername), filter_dist_weight)
        )
        # Need a toward the sun, reward high airmass, with an airmass cutoff basis function.
        bfs.append(
            (bf.NearSunTwilightBasisFunction(nside=nside, max_airmass=max_airmass), 0)
        )
        bfs.append(
            (
                bf.AltAzShadowMaskBasisFunction(nside=nside, shadow_minutes=shadow_minutes, max_alt=max_alt,
                                                min_alt=min_alt),
                0,
            )
        )
        bfs.append(
            (bf.MoonAvoidanceBasisFunction(nside=nside, moon_distance=moon_distance), 0)
        )
        bfs.append((bf.FilterLoadedBasisFunction(filternames=filtername), 0))
        bfs.append((bf.PlanetMaskBasisFunction(nside=nside), 0))
        bfs.append(
            (
                bf.SolarElongationMaskBasisFunction(
                    min_elong=0.0, max_elong=max_elong, nside=nside
                ),
                0,
            )
        )

        bfs.append((bf.NightModuloBasisFunction(pattern=night_pattern), 0))
        # Do not attempt unless the sun is getting high
        bfs.append(
            (
                (
                    bf.SunHighLimitBasisFunction(
                        sun_alt_limit=sun_alt_limit, time_to_12deg=time_to_12deg
                    )
                ),
                0,
            )
        )

        # unpack the basis functions and weights
        weights = [val[1] for val in bfs]
        basis_functions = [val[0] for val in bfs]

        # Set huge ideal pair time and use the detailer to cut down the list of observations to fit twilight?
        surveys.append(
            BlobSurvey(
                basis_functions,
                weights,
                filtername1=filtername,
                filtername2=None,
                ideal_pair_time=ideal_pair_time,
                nside=nside,
                exptime=exptime,
                survey_note=survey_name,
                ignore_obs=ignore_obs,
                dither=True,
                nexp=nexp,
                detailers=detailer_list,
                az_range=az_range,
                twilight_scale=False,
                area_required=area_required,
            )
        )
    return surveys


def set_run_info(dbroot=None, file_end="v3.4_", out_dir=".", rate=None, ntoo=None, tfilt=""):
    """Gather versions of software used to record
    """
    extra_info = {}
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
        extra_info["rubin_scheduler git hash"] = subprocess.check_output(["cat", hash_file])
    except subprocess.CalledProcessError:
        pass

    # Use the filename of the script to name the output database
    if dbroot is None:
        fileroot = os.path.basename(sys.argv[0]).replace(".py", "") + "_"
    else:
        fileroot = dbroot + "_"
    fileroot += "r%i_d%i_%s_" % (rate, ntoo, tfilt)
    fileroot = os.path.join(out_dir, fileroot + file_end)
    return fileroot, extra_info


def run_sched(
    scheduler,
    survey_length=365.25,
    nside=32,
    filename=None,
    verbose=False,
    extra_info=None,
    illum_limit=40.0,
    mjd_start=60796.0,
    event_table=None,
    sim_to_o=None,
):
    """Run survey
    """
    n_visit_limit = None
    fs = SimpleFilterSched(illum_limit=illum_limit)
    observatory = ModelObservatory(nside=nside, mjd_start=mjd_start, sim_to_o=sim_to_o)
    observatory, scheduler, observations = sim_runner(
        observatory,
        scheduler,
        survey_length=survey_length,
        filename=filename,
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
    too_filters = args.too_filters
    too_days = args.too_days

    # Be sure to also update and regenerate DDF grid save file if changing mjd_start
    mjd_start = 60796.0
    per_night = True  # Dither DDF per night

    camera_ddf_rot_limit = 75.0  # degrees

    fileroot, extra_info = set_run_info(dbroot=dbroot, file_end="v3.4_", out_dir=out_dir,
                                        rate=too_rate, ntoo=too_days, tfilt=too_filters)

    pattern_dict = {
        1: [True],
        2: [True, False],
        3: [True, False, False],
        4: [True, False, False, False],
        # 4 on, 4 off
        5: [True, True, True, True, False, False, False, False],
        # 3 on 4 off
        6: [True, True, True, False, False, False, False],
        7: [True, True, False, False, False, False],
    }
    neo_night_pattern = pattern_dict[neo_night_pattern]
    reverse_neo_night_pattern = [not val for val in neo_night_pattern]

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

    repeat_night_weight = None

    sim_ToOs, event_table = generate_events(
        nside=nside, survey_length=365*2., rate=too_rate*5,
        mjd_start=mjd_start + 365*3
    )

    # Use the Almanac to find the position of the sun at the start of survey
    almanac = Almanac(mjd_start=mjd_start)
    sun_moon_info = almanac.get_sun_moon_positions(mjd_start)
    sun_ra_start = sun_moon_info["sun_RA"].copy()

    footprints = make_rolling_footprints(
        fp_hp=footprints_hp,
        mjd_start=mjd_start,
        sun_ra_start=sun_ra_start,
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

    gaps_night_pattern = [True] + [False] * nights_off

    long_gaps = gen_long_gaps_survey(
        nside=nside,
        footprints=footprints,
        night_pattern=gaps_night_pattern,
    )

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

    ddfs = ddf_surveys(
        detailers=details,
        season_unobs_frac=ddf_season_frac,
        euclid_detailers=euclid_detailers,
        nside=nside,
    )

    greedy = gen_greedy_surveys(nside, nexp=nexp, footprints=footprints)
    neo = generate_twilight_near_sun(
        nside,
        night_pattern=neo_night_pattern,
        filters=neo_filters,
        n_repeat=neo_repeat,
        footprint_mask=footprint_mask,
        max_airmass=neo_am,
        max_elong=neo_elong_req,
        area_required=neo_area_req,
    )
    blobs = generate_blobs(
        nside,
        nexp=nexp,
        footprints=footprints,
        mjd_start=mjd_start,
    )
    twi_blobs = generate_twi_blobs(
        nside,
        nexp=nexp,
        footprints=footprints,
        wfd_footprint=wfd_footprint,
        repeat_night_weight=repeat_night_weight,
        night_pattern=reverse_neo_night_pattern,
    )

    # Set up the damn ToO kwargs
    times = [1, 2, 4, 24, 48]
    filters_at_times = [too_filters] * 3 + ["gz", "gz"]
    nvis = [1, 1, 1, 6, 6]

    if too_days == 2:
        too_nfollow = 4
    elif too_days == 3:
        too_nfollow = 5

    camera_rot_limits = [-80., 80.]
    detailer_list = []
    detailer_list.append(detailers.CameraRotDetailer(min_rot=np.min(camera_rot_limits),
                                                     max_rot=np.max(camera_rot_limits)))
    detailer_list.append(detailers.Rottep2RotspDesiredDetailer())

    toos = [
        ToO_scripted_survey(
            [],
            nside=nside,
            followup_footprint=too_footprint,
            times=times[0:too_nfollow],
            filters_at_times=filters_at_times[0:too_nfollow],
            nvis=nvis[0:too_nfollow],
            detailers=detailer_list,
        )
    ]

    surveys = [toos, ddfs, long_gaps, blobs, twi_blobs, neo, greedy]

    scheduler = CoreScheduler(surveys, nside=nside)

    if args.setup_only:
        return scheduler
    else:
        years = np.round(survey_length / 365.25)
        observatory, scheduler, observations = run_sched(
            scheduler,
            survey_length=survey_length,
            verbose=verbose,
            filename=os.path.join(fileroot + "%iyrs.db" % years),
            extra_info=extra_info,
            nside=nside,
            illum_limit=illum_limit,
            mjd_start=mjd_start,
            event_table=event_table,
            sim_to_o=sim_ToOs,
        )
        return observatory, scheduler, observations


def sched_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", dest="verbose", action="store_true", help="Print more output")
    parser.set_defaults(verbose=False)
    parser.add_argument("--survey_length", type=float, default=365.25 * 10, help="Survey length in days")
    parser.add_argument("--out_dir", type=str, default="", help="Output directory")
    parser.add_argument("--maxDither", type=float, default=0.7, help="Dither size for DDFs (deg)")
    parser.add_argument(
        "--moon_illum_limit",
        type=float,
        default=40.0,
        help="illumination limit to remove u-band",
    )
    parser.add_argument("--nexp", type=int, default=2, help="Number of exposures per visit")
    parser.add_argument("--rolling_nslice", type=int, default=2, help="Number of independent rolling stripes")
    parser.add_argument("--rolling_strength", type=float, default=0.9, help="Rolling strength between 0-1")
    parser.add_argument("--dbroot", type=str, help="Database root")
    parser.add_argument(
        "--ddf_season_frac", type=float, default=0.2, help="How much of season to use for DDFs"
    )
    parser.add_argument("--nights_off", type=int, default=3, help="For long gaps")
    parser.add_argument(
        "--neo_night_pattern", type=int, default=4, help="Which night pattern to use for inner solar system"
    )
    parser.add_argument("--neo_filters", type=str, default="riz", help="Filters for inner solar system")
    parser.add_argument(
        "--neo_repeat", type=int, default=4, help="Number of repeat visits for inner solar system"
    )
    parser.add_argument("--neo_am", type=float, default=2.5, help="Airmass limit for twilight NEO visits")
    parser.add_argument(
        "--neo_elong_req", type=float, default=45.0, help="Solar elongation required for inner solar system"
    )
    parser.add_argument(
        "--neo_area_req",
        type=float,
        default=0.0,
        help="Sky area required before attempting inner solar system",
    )
    parser.add_argument(
        "--setup_only",
        dest="setup_only",
        default=False,
        action="store_true",
        help="Only construct scheduler, do not simulate",
    )
    parser.add_argument(
        "--nside", type=int, default=32, help="Nside should be set to default (32) except for tests."
    )
    
    parser.add_argument("--too_rate", type=float, default=10, help="N events per year")
    parser.add_argument("--too_days", type=int)
    parser.add_argument("--too_filters", type=str, default="gz")

    return parser


if __name__ == "__main__":

    parser = sched_argparser()
    args = parser.parse_args()
    example_scheduler(args)

