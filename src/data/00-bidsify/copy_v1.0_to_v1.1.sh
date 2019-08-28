#!/bin/sh

FROM_DIR=${OAK}/projects/causcon-bids.v1.0
TO_DIR=${OAK}/projects/causcon-bids.v1.1

TASK=multisourceinterference

for sub in ${FROM_DIR}/sub-*
do
    file=`ls ${sub}/ses-d1/func/*${TASK}*.nii.gz`
    newfile=`echo $file | sed "s?.*bids.v1.0/??g"`
    echo $newfile
    cp ${FROM_DIR}/$newfile ${TO_DIR}/$newfile    
done
