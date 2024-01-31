import numpy as np
from rubin_scheduler.scheduler.utils import restore_scheduler
from baseline import example_scheduler, sched_argparser
from rubin_scheduler.scheduler.model_observatory import ModelObservatory


if __name__ == "__main__":

    parser = sched_argparser()
    args = parser.parse_args()
    args.setup_only = True

    scheduler = example_scheduler(args)

    observation_id = 1359637-4

    #observation_id = 5000

    observatory = ModelObservatory()

    scheduler, observatory = restore_scheduler(observation_id, scheduler,
                                               observatory, 'baseline_v3.4_10yrs.db', filter_sched=None, fast=True)

    conditions = observatory.return_conditions()
    conditions.mounted_filters = ['u', 'g', 'r', 'i', 'z']

    for sl in scheduler.survey_lists:
        for survey in sl:
            print(survey, survey.calc_reward_function(conditions))

    obs = np.concatenate(scheduler.survey_lists[2][4].generate_observations(conditions))
    import pdb ; pdb.set_trace()

    sur = scheduler.survey_lists[2][4]
    