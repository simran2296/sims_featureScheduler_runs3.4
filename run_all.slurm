#!/bin/bash

## Job Name

#SBATCH --job-name=big_run3.4

## Allocation Definition

## On mox and ikt, the account and partition options should be the same.
#SBATCH --account=astro
#SBATCH --partition=astro

## Resources

## Nodes

#SBATCH --nodes=3
#SBATCH --ntasks-per-node=28

## Walltime (hours:min:sec) Do not specify a walltime substantially more than your job needs.

#SBATCH --time=24:00:00

## Memory per node. It is important to specify the memory since the default memory is very small.

## For mox, --mem may be more than 100G depending on the memory of your nodes.

## For ikt, --mem may be 58G or more depending on the memory of your nodes.

## See above section on "Specifying memory" for choices for --mem.

#SBATCH --mem=500G

## Specify the working directory for this job

#SBATCH --chdir=/gscratch/scrubbed/yoachim/sims_featureScheduler_runs3.4

##turn on e-mail notification

#SBATCH --mail-type=ALL

#SBATCH --mail-user=yoachim@uw.edu

## export all your environment variables to the batch job session

#SBATCH --export=all

## Set up the evironment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate guro
export OPENBLAS_NUM_THREADS=1

cd /gscratch/scrubbed/yoachim/sims_featureScheduler_runs3.4

## run all the baseline commands in parallel
module load parallel-20170722


# Need to activate for each I guess
cat run_all.sh | xargs -I'{}' echo "cd /gscratch/scrubbed/yoachim/sims_featureScheduler_runs3.4 ; source setup_env ; {}" > final.sh

scontrol show hostnames > list_of_nodes
cpu_count=28

cat final.sh | parallel --sshloginfile list_of_nodes -j $cpu_count

