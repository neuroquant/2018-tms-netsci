#!/bin/sh

STUDY_DIR=CausalConnectome
FPREP_DIR=${PI_SCRATCH}/COMET/${STUDY_DIR}/derivatives/fmriprep-fsl/fmriprep
FSF_DIR=${PI_SCRATCH}/COMET/${STUDY_DIR}/derivatives/fmriprep-fsl/freesurfer
DERIVATIVES=${PI_SCRATCH}/COMET/${STUDY_DIR}/derivatives/fmriprep-fsl/denoiser

#echo "rclone confound derivatives"
#rclone copy --update drive:Data/${STUDY_DIR}/derivatives/fmriprep-v1.2.3/confound ${PI_SCRATCH}/COMET/${STUDY_DIR}/derivatives/fmriprep-fsl/confound


echo "rclone copy --update drive:Data/${STUDY_DIR}/derivatives/fmriprep-v1.2.3/fmriprep ${FPREP_DIR}"
#rclone copy --update drive:Data/${STUDY_DIR}/derivatives/fmriprep-v1.2.3/fmriprep/sub-NTHC1043 ${FPREP_DIR}/sub-NTHC1043
rclone copy --update -v drive:Data/${STUDY_DIR}/derivatives/fmriprep-v1.2.3/fmriprep/sub-TEHC2001 ${FPREP_DIR}/sub-TEHC2001