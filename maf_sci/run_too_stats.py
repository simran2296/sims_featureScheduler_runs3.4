import argparse
import numpy as np
import pandas as pd
import sqlite3
import os

from rubin_sim.maf.slicers import UserPointsSlicer
import rubin_sim.maf.metrics as metrics
import rubin_sim.maf as maf


def gen_too_slicer(mjd_starts, ra_rad, dec_rad, file_indx=0, distance=50):
    slicer = UserPointsSlicer(np.degrees(ra_rad), np.degrees(dec_rad), lat_lon_deg=True, badval=0)
    slicer.slice_points["peak_time"] = mjd_starts
    # XXX--which file to use?
    slicer.slice_points["file_indx"] = file_indx
    # Distance to use? Probably in Mpc
    slicer.slice_points["distance"] = distance
    return slicer


def gen_table(db_file, out_dir):

    run_name = os.path.basename(db_file).replace('.db', '')
    con = sqlite3.connect(db_file)
    too_events = pd.read_sql('select * from events', con)
    too_obs = pd.read_sql('select note from observations where note like "%ToO%"', con)
    con.close()

    bundle_list = []
    slicer = gen_too_slicer(too_events["mjd_start"].values, too_events["ra"].values, too_events["dec"].values)
    metric = maf.PrestoColorKNePopMetric(mjd0=0, skyregion="extragalactic", metric_name="presto")
    summary = [
            # XXX--maybe this should be a probability sum?
            metrics.SumMetric(metric_name="Total detected"),
            metrics.CountMetric(metric_name="Total lightcurves in footprint"),
            metrics.CountMetric(metric_name="Total lightcurves on sky", mask_val=0),
            metrics.MeanMetric(mask_val=0, metric_name="Fraction detected of total (mean)"),
        ]
    bundle = maf.MetricBundle(
            metric,
            slicer,
            "",
            run_name=run_name,
            summary_metrics=summary,
        )
    bundle_list.append(bundle)
    metric = maf.KNePopMetric(mjd0=0)
    bundle = maf.MetricBundle(
            metric,
            slicer,
            "",
            run_name=run_name,
            summary_metrics=summary,
        )
    bundle_list.append(bundle)

    bd = maf.metricBundles.make_bundles_dict_from_list(bundle_list)
    bg = maf.metricBundles.MetricBundleGroup(bd, db_file, save_early=False)
    bg.run_all()

    counts_df = too_obs.value_counts()
    arr = np.zeros(len(too_events), dtype=list(zip(['t0', 't1', 't2', 't4', 't24', 't48'], [int, int, int, int, int, int])))
    for label, val in zip(counts_df.index, counts_df.values):
        label = label[0].replace('ToO, ', '')
        i, t = label.split('_')
        i = int(i.replace('_', ''))
        arr[t][i] = val

    # Add the number of observations per time to the event table
    for key in arr.dtype.names:
        too_events[key] = arr[key]

    # Add the run name
    too_events["run_name"] = run_name

    keys = ["presto_presto_color_detect", "presto_score_p", 'presto_score_s', 
            'KNePopMetric_blue_color_detect', 'KNePopMetric_multi_color_detect',
            'KNePopMetric_multi_detect', 'KNePopMetric_red_color_detect',
            'KNePopMetric_ztfrest_simple', 'KNePopMetric_ztfrest_simple_blue',
            'KNePopMetric_ztfrest_simple_red']
    for key in keys:
        too_events[key] = bd[key].metric_values

    too_events.to_hdf(os.path.join(out_dir, run_name + '.h5'), 'too')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--db_file",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default=None,
    )

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    gen_table(args.db_file, out_dir=args.out_dir)
