import numpy as np
import sqlite3
import pandas as pd
from initial import sched_argparser, example_scheduler
from rubin_scheduler.scheduler.utils import restore_scheduler
from rubin_scheduler.scheduler.model_observatory import ModelObservatory


if __name__ == "__main__":
    args = sched_argparser()
    args = args.parse_args([])
    args.setup_only = True

    scheduler = example_scheduler(args)

    mo = ModelObservatory(mjd_start=60796.0)
    sched, mo = restore_scheduler(188, scheduler, mo, 'm2_initial_v3.4_0yrs.db')
    conditions = mo.return_conditions()

    for sl in sched.survey_lists:
        for survey in sl:
            print(survey.survey_name, np.nanmax(survey.calc_reward_function(conditions)))

    for bf in sched.survey_lists[2][0].basis_functions:
        print(bf, np.nanmax(bf(conditions)))

    import pdb ; pdb.set_trace()
