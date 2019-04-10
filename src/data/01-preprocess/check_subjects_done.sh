#!/bin/sh

DATADIR=${OAK}/projects/causcon-bids.v1.0
ALL_SITES=(LaMFG LFp LIFGAnat LIFGBLA LpMFG LVmFp RaMFG RFEF RFp RIFGAnat RIFGBLA RIFJ RIPL RM1 RpMFG RpreSMA RVmFp)

for (( siteno=0; siteno<${#ALL_SITES[@]}; siteno++))
do

    TMS_SITE=${ALL_SITES[$siteno]}

    ls ${DATADIR}/sub-*/ses-d2/func/*.nii.gz | grep ${TMS_SITE} | sed -e "s?.*/sub-??g" | sed -e "s?_ses-d2.*??g" > ${TMS_SITE}_bids.txt

    ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/fmriprep/sub-*/ses-d2/func/*MNI152*preproc_bold.nii.gz | grep ${TMS_SITE} | sed -e "s?.*func/sub-??g" | sed -e "s?_ses-d2.*??g" > ${TMS_SITE}_files.txt

    ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/fmriprep/sub-*/ses-d2/func/*MELODIC*.tsv | grep ${TMS_SITE} | sed -e "s?.*func/sub-??g" | sed -e "s?_ses-d2.*??g" > ${TMS_SITE}_ica.txt

done
