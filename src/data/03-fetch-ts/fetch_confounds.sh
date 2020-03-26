#!/bin/sh

export BASEDIR=${PI_SCRATCH}/COMET/CausalConnectome
export DERIVATIVES=${BASEDIR}/derivatives/fmriprep-fsl/fmriprep
export ROITS_DIR=${BASEDIR}/derivatives/fmriprep-fsl/denoiser
export CONFOUND_DIR=${BASEDIR}/derivatives/fmriprep-fsl/confound

mkdir -p ${CONFOUND_DIR}

# rsync --update --list-only ${DERIVATIVES}/sub-*/ses-*/func/*.tsv ${CONFOUND_DIR}/.
rsync --update ${DERIVATIVES}/sub-*/ses-*/func/*.tsv ${CONFOUND_DIR}/.

#cp ${DERIVATIVES}/sub-*/ses-*/func/*.tsv ${CONFOUND_DIR}/.
#cp ${DERIVATIVES}/sub-*/ses-*/func/*.csv ${CONFOUND_DIR}/.
