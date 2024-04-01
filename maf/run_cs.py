import argparse
import os
import rubin_sim.maf as maf


def run(filename):

    nside = 64

    run_name = os.path.basename(filename).replace('.db', '')

    out_dir = run_name + '_sci'
    results_db = maf.db.ResultsDb(out_dir=out_dir)

    bundle_list = []

    summary_metrics = [maf.SumMetric()]

    displayDict = {}

    displayDict["group"] = "Variables/Transients"
    displayDict["subgroup"] = "Color and slope"
    displayDict["caption"] = "Number of times a color and slope are measured in a night"

    metric = maf.ColorSlope2NightMetric()
    sql = 'visitExposureTime > 19'
    slicer = maf.slicers.HealpixSlicer(nside=nside)
    plot_dict = {}

    bundle_list.append(maf.MetricBundle(metric, slicer, sql,
                                        run_name=run_name,
                                        plot_dict=plot_dict,
                                        summary_metrics=summary_metrics,
                                        display_dict=displayDict))

    displayDict["group"] = "Variables/Transients"
    displayDict["subgroup"] = "Color and slope"
    displayDict["caption"] = "Number of times a color and slope are measured over 2 nights."

    metric = maf.ColorSlopeMetric()
    bundle_list.append(maf.MetricBundle(metric, slicer, sql,
                                        run_name=run_name,
                                        plot_dict=plot_dict,
                                        summary_metrics=summary_metrics,
                                        display_dict=displayDict))

    bd = maf.metricBundles.make_bundles_dict_from_list(bundle_list)
    bg = maf.metricBundles.MetricBundleGroup(bd, filename, out_dir=out_dir, results_db=results_db)
    bg.run_all()
    bg.plot_all()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=str, default=None)
    args = parser.parse_args()

    run(args.db)
    
