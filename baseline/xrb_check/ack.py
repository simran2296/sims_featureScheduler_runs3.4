import rubin_sim.maf as maf
import matplotlib.pylab as plt
import sqlite3
import pandas as pd
import os



out_dir = 'temp'
results_db = maf.ResultsDb(out_dir=out_dir)

filename = '../baseline_v3.4_10yrs.db'
runName = os.path.basename(filename).replace('.db', '')

con = sqlite3.connect(filename)
mjd0_df = pd.read_sql("select min(observationStartMJD) from observations;", con)
con.close()
mjd0 = mjd0_df.values.min()

bundleList = []
n_events = 10000
xrbslicer = maf.generate_xrb_pop_slicer(n_events=n_events)
metric = maf.XRBPopMetric(output_lc=False, mjd0=mjd0)
xrb_summaryMetrics = [
    maf.SumMetric(metric_name="Total detected"),
    maf.CountMetric(metric_name="Total lightcurves in footprint"),
    maf.CountMetric(metric_name="Total lightcurves on sky", mask_val=0),
    maf.MeanMetric(metric_name="Fraction detected in footprint"),
    maf.MeanMetric(mask_val=0, metric_name="Fraction detected of total"),
    maf.MedianMetric(metric_name="Median"),
    maf.MeanMetric(metric_name="Mean"),
]

bundleList.append(
    maf.MetricBundle(
        metric,
        xrbslicer,
        "",
        run_name=runName,
        summary_metrics=xrb_summaryMetrics,
    )
)
bdict = maf.make_bundles_dict_from_list(bundleList)

group = maf.MetricBundleGroup(bdict, filename, out_dir=out_dir, results_db=results_db,
        )

group.run_all()
group.plot_all(closefigs=False)

print(bdict['XRBPopMetric_early_detect'].summary_values)
