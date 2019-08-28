#!/bin/sh

for sub in $(seq 19 19)
do
   # sbatch -p aetkin,owners,normal --array=$sub -N 1 -n 16 --mem=64000 -t 8:00:00 --begin=now+${sub}minute run_fmriprep.sh

   sbatch -p aetkin,owners,normal --array=$sub -N 1 -n 16 --mem=64000 -t 4:00:00 --begin=now run_fmriprep.sh

   # 6 hours was enough to complete freesurfer and T1w processing but not the rest. 
#   sbatch -p aetkin,owners,normal --array=$sub -N 1 -n 16 --mem=64000 -t 4:00:00 --begin=now+8hour run_fmriprep.sh
done
