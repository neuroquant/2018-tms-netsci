#!/bin/sh

SUBJID=$(sed "1p;1d" ../00-bidsify/subjects_tms.txt)
#SUBJID=(NTHC1001)

for subj in ${SUBJID}
do
    # for taskno in $(seq 1 18)
    # do
    # --array=$taskno 
    sbatch -n 1 -t 1:00:00 --mem=32000 -p owners,normal --export=SUBJECT=sub-${subj} run_denoiser_T1w.sh
    # done
done
