#!/bin/sh

for sub in $(seq 1 5)
do
    sbatch -p aetkin,owners,normal --array=$sub -N 1 -n 16 --mem=64000 -t 8:00:00 --begin=now+${sub}minute run_fmriprep.sh

#    sbatch -p aetkin,owners,normal --array=$sub -N 1 -n 16 --mem=64000 -t 2:00:00 --begin=now run_fmriprep.sh

   sbatch -p aetkin,owners --array=$sub -N 1 -n 16 --mem=64000 -t 4:00:00 --begin=now+8hour run_fmriprep.sh
done
