#!/bin/sh

ALL_SITES=(LaMFG LFp LIFGAnat LIFGBLA LpMFG LVmFp RaMFG RFEF RFp RIFGAnat RIFGBLA RIFJ RIPL RM1 RpMFG RpreSMA RVmFp)

for (( siteno=0; siteno<${#ALL_SITES[@]}; siteno++))
do

    TMS_SITE=${ALL_SITES[$siteno]}

    ls /scratch/groups/aetkin/COMET/CausalConnectome/BIDS/sub-*/ses-d2/func/*.nii.gz | grep ${TMS_SITE} | sed -e "s?.*/sub-??g" | sed -e "s?_ses-d2.*??g" > ${TMS_SITE}_bids.txt

    ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep/fmriprep/sub-*/ses-d2/func/*MNI152*preproc_bold.nii.gz | grep ${TMS_SITE} | sed -e "s?.*func/sub-??g" | sed -e "s?_ses-d2.*??g" > ${TMS_SITE}_files.txt

done