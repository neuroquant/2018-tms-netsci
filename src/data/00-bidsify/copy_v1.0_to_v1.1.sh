#!/bin/sh

#FROM_DIR=${OAK}/projects/causcon-bids.v1.0
FROM_DIR=/scratch/groups/aetkin/COMET/CausalConnectome/BIDS
TO_DIR=${OAK}/projects/causcon-bids.v1.1

TASK=multisourceinterference

for sub in ${FROM_DIR}/sub-*
do
    file=`ls ${sub}/ses-d1/func/*${TASK}*.nii.gz`
    #newfile=`echo $file | sed "s?.*bids.v1.0/??g"`
    newfile=`echo $file | sed "s?.*BIDS/??g"`
    
    echo $newfile
    rsync --update --checksum --copy-links ${FROM_DIR}/$newfile ${TO_DIR}/$newfile
done
