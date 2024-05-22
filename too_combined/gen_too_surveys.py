import numpy as np
from rubin_scheduler.site_models import _read_fields
from rubin_scheduler.scheduler.surveys import (
    ScriptedSurvey,
    BaseMarkovSurvey,
)
from rubin_scheduler.scheduler.utils import (
    scheduled_observation,
    comcam_tessellate,
    gnomonic_project_toxy,
    tsp_convex,
)
from rubin_scheduler.utils import (
    _hpid2_ra_dec,
    _ra_dec2_hpid,
    _approx_ra_dec2_alt_az,
    _angular_separation,
)


def gen_too_surveys(too_footprint, detailer_list=[], nside=32):
    surveys = []

    # Let's make a footprint to follow up ToO events
    too_footprint = too_footprint * 0 + np.nan
    too_footprint[np.where(too_footprint > 0)[0]] = 1.0

    # Set up the damn ToO kwargs
    too_filters = "gz"
    times = [1, 2, 4, 24, 48]
    filters_at_times = [too_filters] * 3 + ["gz", "gz"]
    nvis = [1, 1, 1, 6, 6]
    too_nfollow = 5

    s1 = ToO_scripted_survey(
        [],
        nside=nside,
        followup_footprint=too_footprint,
        times=times[0:too_nfollow],
        filters_at_times=filters_at_times[0:too_nfollow],
        nvis=nvis[0:too_nfollow],
        detailers=detailer_list,
    )

    surveys.append(s1)
    return surveys


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
                            obs["rotSkyPos"] = (
                                0  # XXX--maybe throw a rotation detailer in here
                            )
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
    # scale = 1e6
    # pointing_x = np.round(pointing_x * scale).astype(int)
    # pointing_y = np.round(pointing_y * scale).astype(int)
    # Now I have a bunch of x,y pointings. Drop into TSP solver to get an effiencent route
    towns = np.vstack((pointing_x, pointing_y)).T
    # Leaving optimize=False for speed. The optimization step doesn't usually improve much.
    better_order = tsp_convex(towns, optimize=False)
    return better_order


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
