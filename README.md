[![DOI](https://zenodo.org/badge/718290536.svg)](https://zenodo.org/doi/10.5281/zenodo.10126361)

# sims_featureScheduler_runs3.4
Yet another repo for Rubin observing strategy experiements

Completed simulations can be found at: https://s3df.slac.stanford.edu/data/rubin/sim-data/sims_featureScheduler_runs3.4/

Latest analysis of runs often served at:  http://astro-lsst-01.astro.washington.edu:8080/


## Families of simulations

baseline : The latest baseline. No major strategy changes compared to v3.3

ddf_acor : Variations on DDF season length and cadence within season. Intended generally to extend season length without sacrificing mid-season cadence.

ddf_loaded : Executing DDF observations only when a filter is already loaded. Intended to reduce filter changes, and potentially spread DDF observations over more days.

ddf_loaded_half : Like ddf_loaded, but scheduling DDF sequences of half-length to execute twice as often.

early_science : TBD

extra_ss : Running more twilight near-sun observations than the baseline. Mainly intended to show that if we have ToO observations for earth-interior objects the impact of the rest of the survey would be negligable.

initial : An initial simulation before finilizing the new baseline.

mw : Tests what happens if the LMC/SMC region is excluded from rolling.

noroll : A baseline without any rolling cadence.

roll_pause/roll_pause2 : Tests to try and make more uniform surveys at release dates.

too : Updating ToO simulations from earlier verions

too_ignore : Like too sims, but the main survey ignores ToO observations (i.e., ToO observations don't get to count twice.)

too_elab : ToO simulations with a more elaborate followup strategy

too_elab_heavy : Similar to too_elab, but going deeper on first-day ToO observations (4 visits per filter rather than 1)

too_elab_lots : same as too_elab, should be merged added to those in case I forgot.

vary_u_exp_n : variations on u-band exposure time and number of u-band visits. Motivated by decrease in u throughput going to tripple silver coating.

weather : Running the baseline for different weather realizations. Useful for checking expected levels of variation in science metrics.

