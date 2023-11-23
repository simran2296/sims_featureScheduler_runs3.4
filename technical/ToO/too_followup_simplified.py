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
    seed=42,
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
    print('ToO Events: \n',table)


    
    pointings = []
    #for 
    return pointings, pointing_table


def set_run_info(dbroot=None, file_end="v3.4_", out_dir=".", rate=None, ntoo=None):
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
    fileroot += "%i_%i" % (rate, ntoo) # Append the too rate and nfollow to the db name
    fileroot = os.path.join(out_dir, fileroot + file_end) # Append the version to the filename
    return fileroot, extra_info


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

    mjd_start = 60796.0
    per_night = True  # Dither DDF per night

    camera_ddf_rot_limit = 75.0  # degrees
    '''
    fileroot, extra_info = set_run_info(
        dbroot=dbroot,
        file_end="v3.4_",
        out_dir=out_dir,
        rate=too_rate,
        ntoo=too_nfollow,
    )
    '''
    #  print('fileroot: ',fileroot,'\n Extra Info: ', extra_info)

    sim_ToOs, event_table = incoming_event(nside=nside,mjd_start=mjd_start)

    
    #sim_ToOs, event_table = generate_events(
    #    nside=nside, survey_length=survey_length, rate=too_rate, mjd_start=mjd_start
    # )

    #print('Events: ',sim_ToOs,'\n Event Table: ', event_table)




def sched_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", dest="verbose", action="store_true")
    parser.set_defaults(verbose=False)
    parser.add_argument("--survey_length", type=float, default=365.25 * 10)
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
    parser.add_argument("--filters", type=str, default="gz")
    parser.add_argument("--nfollow", type=int, default=1)


    parser.add_argument("json", type = argparse.FileType('r'))
    return parser


if __name__ == "__main__":
    parser = sched_argparser()
    args = parser.parse_args()
    example_scheduler(args)
