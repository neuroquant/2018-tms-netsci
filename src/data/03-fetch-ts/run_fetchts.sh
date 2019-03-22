#!/bin/sh

export BASEDIR=${PI_SCRATCH}/COMET/CausalConnectome
export DERIVATIVES=${BASEDIR}/derivatives/denoiser
export ROITS_DIR=${BASEDIR}/derivatives/denoiser

export VIRTUALENVWRAPPER_PYTHON=/share/software/user/open/python/3.6.1/bin/python3
alias source_venv='source ~/.local/bin/virtualenvwrapper.sh'
source_venv
workon ni-denoise

SUBJID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../00-bidsify/subjects_tms_part2.txt)
#SUBJID=(NTHC1003 NTHC1009)

for subjno in ${SUBJID}
do
    #subj=sub-${SUBJID[$subjno]}
    #echo $subj
    
    # for taskno in $(seq 1 18)
    # do
    # --array=$taskno 
    #sbatch -n 1 -t 1:00:00 --mem=32000 --export=SUBJECT=sub-${subj} 
    # done
    
    export SUBJECT=$subjno
    python fetch_roi_timeseries.py
done
