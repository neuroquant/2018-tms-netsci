#!/bin/sh

#SUBJID=$(sed -n "54,116p;" ../00-bidsify/subjects.txt)
SUBJID=$(sed "1p;30q" ../01-preprocess/subjects_missing.txt)
#SUBJID=(NTHC1001)

for subj in ${SUBJID}
do
    # for taskno in $(seq 1 18)
    # do
    # --array=$taskno 
    sbatch -n 1 -t 1:00:00 --mem=32000 -p owners,normal --export=SUBJECT=sub-${subj} run_denoiser.sh
    # done
done
