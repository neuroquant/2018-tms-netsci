#!/bin/sh

DATADIR=${OAK}/projects/causcon-bids.v1.1
ALL_SITES=(LaMFG LFp LIFGAnat LIFGBLA LpMFG LVmFp RaMFG RFEF RFp RIFGAnat RIFGBLA RIFJ RIPL RM1 RpMFG RpreSMA RVmFp)

ls ${DATADIR}/sub-*/ses-d1/func/*.nii.gz | grep rest | sed -e "s?.*/sub-??g" | sed -e "s?_ses-d1.*??g" > rest_bids.txt

	ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/confound/sub-*MELODIC*.tsv | grep rest | sed -e "s?.*confound/sub-??g" | sed -e "s?_ses-d1.*??g" > rest_ica.txt

    #ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/fmriprep/sub-*/ses-d1/func/*MELODIC*.tsv | grep rest | sed -e "s?.*func/sub-??g" | sed -e "s?_ses-d1.*??g" > rest_ica.txt

    #ls /scratch/groups/aetkin/COMET/CausalConnectome/work/fmriprep/fmriprep_wf/single_subject_*/func_preproc_*/ica_aroma_wf/ica_aroma_confound_extraction/*Confounds.tsv | grep rest | sed -e "s?.*single_subject_??g" | sed -e "s?_wf/func.*??g" > rest_ica.txt
    
grep -v -F -f  rest_ica.txt rest_bids.txt > rest_missing.txt

echo "Rest counted"

ls ${DATADIR}/sub-*/ses-d1/func/*.nii.gz | grep monetary | sed -e "s?.*/sub-??g" | sed -e "s?_ses-d1.*??g" > mid_bids.txt

	ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/confound/sub-*MELODIC*.tsv | grep monetary | sed -e "s?.*confound/sub-??g" | sed -e "s?_ses-d1.*??g" > mid_ica.txt

    #ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/fmriprep/sub-*/ses-d1/func/*MELODIC*.tsv | grep monetary | sed -e "s?.*func/sub-??g" | sed -e "s?_ses-d1.*??g" > mid_ica.txt

    #ls /scratch/groups/aetkin/COMET/CausalConnectome/work/fmriprep/fmriprep_wf/single_subject_*/func_preproc_*/ica_aroma_wf/ica_aroma_confound_extraction/*Confounds.tsv | grep monetary | sed -e "s?.*single_subject_??g" | sed -e "s?_wf/func.*??g" > mid_ica.txt
    
grep -v -F -f  mid_ica.txt mid_bids.txt > mid_missing.txt

echo "MID counted"

ls ${DATADIR}/sub-*/ses-d1/func/*.nii.gz | grep multisource | sed -e "s?.*/sub-??g" | sed -e "s?_ses-d1.*??g" > msit_bids.txt

	ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/confound/sub-*MELODIC*.tsv | grep multisource | sed -e "s?.*confound/sub-??g" | sed -e "s?_ses-d1.*??g" > msit_ica.txt

   # ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/fmriprep/sub-*/ses-d1/func/*MELODIC*.tsv | grep multisource | sed -e "s?.*func/sub-??g" | sed -e "s?_ses-d1.*??g" > msit_ica.txt

   # ls /scratch/groups/aetkin/COMET/CausalConnectome/work/fmriprep/fmriprep_wf/single_subject_*/func_preproc_*/ica_aroma_wf/ica_aroma_confound_extraction/*Confounds.tsv | grep multisource | sed -e "s?.*single_subject_??g" | sed -e "s?_wf/func.*??g" > msit_ica.txt
    
grep -v -F -f  msit_ica.txt msit_bids.txt > msit_missing.txt

echo "MSIT counted"


for (( siteno=0; siteno<${#ALL_SITES[@]}; siteno++))
do

    TMS_SITE=${ALL_SITES[$siteno]}

    ls ${DATADIR}/sub-*/ses-d2/func/*.nii.gz | grep ${TMS_SITE} | sed -e "s?.*/sub-??g" | sed -e "s?_ses-d2.*??g" > ${TMS_SITE}_bids.txt

    ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/fmriprep/sub-*/ses-d2/func/*MNI152*preproc_bold.nii.gz | grep ${TMS_SITE} | sed -e "s?.*func/sub-??g" | sed -e "s?_ses-d2.*??g" > ${TMS_SITE}_files.txt

	ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/confound/sub-*MELODIC*.tsv | grep ${TMS_SITE} | sed -e "s?.*confound/sub-??g" | sed -e "s?_ses-d2.*??g" > ${TMS_SITE}_ica.txt

    #ls /scratch/groups/aetkin/COMET/CausalConnectome/derivatives/fmriprep-fsl/fmriprep/sub-*/ses-d2/func/*MELODIC*.tsv | grep ${TMS_SITE} | sed -e "s?.*func/sub-??g" | sed -e "s?_ses-d2.*??g" > ${TMS_SITE}_ica.txt
    
    #ls /scratch/groups/aetkin/COMET/CausalConnectome/work/fmriprep/fmriprep_wf/single_subject_*/func_preproc_*/ica_aroma_wf/ica_aroma_confound_extraction/*Confounds.tsv | grep ${TMS_SITE} | sed -e "s?.*single_subject_??g" | sed -e "s?_wf/func.*??g" > ${TMS_SITE}_ica.txt

    grep -v -F -f  ${TMS_SITE}_ica.txt ${TMS_SITE}_bids.txt > ${TMS_SITE}_missing.txt
    
    echo "${TMS_SITE} counted."
done

rm all_missing.txt
sort *_missing.txt | uniq -c > all_missing.txt
