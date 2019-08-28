#!/bin/sh

export BASEDIR=${PI_SCRATCH}/COMET/CausalConnectome
export DERIVATIVES=${BASEDIR}/derivatives/fmriprep-fsl/
export WORK_DIR=${BASEDIR}/work/work/fmriprep/fmriprep_wf/

for sub in $(seq 1 53)
do 

    SUBJID=$(sed -n "${sub}p" ../00-bidsify/subjects.txt)
    
done
