#!/bin/sh

STUDY_DIR=CausalConnectome
BACKUP_DIR=${PI_SCRATCH}/COMET/${STUDY_DIR}/derivatives/fmriprep

# rclone sync --copy-links --update ${BACKUP_DIR}/derivatives drive:Data/${STUDY_DIR}/derivatives
#
echo "rclone copy --update ${BACKUP_DIR}/fmriprep drive:Data/${STUDY_DIR}/derivatives/fmriprep"
#rclone copy -v --update ${BACKUP_DIR}/fmriprep drive:Data/${STUDY_DIR}/derivatives/fmriprep
rclone copy -v --update ${PI_SCRATCH}/COMET/${STUDY_DIR}/derivatives/denoiser drive:Data/${STUDY_DIR}/derivatives/denoiser

#
# rclone copy --checksum --update ${BACKUP_DIR}/derivatives/cpac drive:Data/${STUDY_DIR}/derivatives/cpac
#RAWBIDS=BIDS
# rclone sync --copy-links --update ${BACKUP_DIR}/${RAWBIDS} box:Stanford/datasets/${STUDY_DIR}/RawBIDS
