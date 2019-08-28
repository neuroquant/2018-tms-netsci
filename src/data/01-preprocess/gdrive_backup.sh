#!/bin/sh

STUDY_DIR=CausalConnectome
FPREP_DIR=${PI_SCRATCH}/COMET/${STUDY_DIR}/derivatives/fmriprep-fsl/fmriprep
FSF_DIR=${PI_SCRATCH}/COMET/${STUDY_DIR}/derivatives/fmriprep-fsl/freesurfer
DERIVATIVES=${PI_SCRATCH}/COMET/${STUDY_DIR}/derivatives/fmriprep-fsl/denoiser

# rclone sync --copy-links --update ${BACKUP_DIR}/derivatives drive:Data/${STUDY_DIR}/derivatives
#

# echo "rclone copy -v --update ${DERIVATIVES} drive:Data/${STUDY_DIR}/derivatives/fmriprep-v1.2.3/denoiser"
# rclone copy -v --update ${DERIVATIVES} drive:Data/${STUDY_DIR}/derivatives/fmriprep-v1.2.3/denoiser

echo "rclone copy --update --filter "- *smoothAROMAnonaggr_bold.nii.gz.nii.gz" ${FPREP_DIR} box:Stanford/datasets/CausalConnectome/derivatives/fmriprep-v1.2.3/fmriprep"
#rclone copy -v --update --filter "- *smoothAROMAnonaggr_bold.nii.gz.nii.gz" ${FPREP_DIR} drive:Data/${STUDY_DIR}/derivatives/fmriprep-v1.2.3/fmriprep
rclone copyto -v --update --filter "- *smoothAROMAnonaggr_bold.nii.gz.nii.gz" ${FPREP_DIR} box:Stanford/datasets/CausalConnectome/derivatives/fmriprep-v1.2.3

#echo "rclone copy --update ${FSF_DIR}  drive:Data/${STUDY_DIR}/derivatives/fmriprep-v1.2.3/freesurfer"
#rclone copy -v --update ${FSF_DIR} drive:Data/${STUDY_DIR}/derivatives/fmriprep-v1.2.3/freesurfer

#
# rclone copy --checksum --update ${BACKUP_DIR}/derivatives/cpac drive:Data/${STUDY_DIR}/derivatives/cpac
#RAWBIDS=BIDS
# rclone sync --copy-links --update ${BACKUP_DIR}/${RAWBIDS} box:Stanford/datasets/${STUDY_DIR}/RawBIDS
