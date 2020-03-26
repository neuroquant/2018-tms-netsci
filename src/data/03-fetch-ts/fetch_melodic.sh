#!/bin/sh

export BASEDIR=${PI_SCRATCH}/COMET/CausalConnectome
export DERIVATIVES=${BASEDIR}/derivatives/fmriprep-fsl
export WORK_DIR=${BASEDIR}/work/fmriprep/fmriprep_wf

for sub in $(seq 118 118) #117-159
do 

    SUBJID=$(sed -n "${sub}p" ../00-bidsify/subjects.txt)
    mkdir -p ${DERIVATIVES}/melodic/sub-${SUBJID}
    
    for tdir in ${WORK_DIR}/single_subject_${SUBJID}_wf/func_preproc_ses*
    do 
        taskdir=$(echo ${tdir} | sed -e "s|.*preproc_||g" | sed -e "s|_run.*||g")
        echo "Zipping ${DERIVATIVES}/melodic/sub-${SUBJID}/${taskdir}"
        zip "${DERIVATIVES}/melodic/sub-${SUBJID}/${taskdir}.zip" ${tdir}/ica_aroma_wf/melodic/*
    done
done