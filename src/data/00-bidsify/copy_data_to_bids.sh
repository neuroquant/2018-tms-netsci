#!/bin/sh

############################################
# Specify Input/Output Locations
############################################
BASEDIR=${PI_SCRATCH}/COMET/CausalConnectome/rawdata
STUDYDIR=(structural rest tmsfmri)
OUTPUTDIR=${PI_SCRATCH}/COMET/CausalConnectome
mkdir -p ${OUTPUTDIR}/BIDS
OUTPUTDIR=${PI_SCRATCH}/COMET/CausalConnectome/BIDS
##########################################
# Create list of subjects per TMS task
##########################################
ls ${BASEDIR}/structural/*.nii.gz | sed -e "s?.*CausCon_??g" | sed -e "s?_t1.*??g" > subjects.txt
SUBJID=$(sed "1p;1d" subjects.txt)



#########################################
# BIDSify structural and resting
#########################################

for sub in ${SUBJID[@]}
do
    echo "$sub"

    SUBJDIR=${OUTPUTDIR}/sub-${sub}

    DAY_NO=1
    SES_NO=d${DAY_NO}
    BIDS_SESSION=ses-${SES_NO}

    mkdir -p ${SUBJDIR}
    mkdir -p ${SUBJDIR}/${BIDS_SESSION}/

    # #################################
    # ####### Structural ##########
    # #################################
	IMG_DIR=structural
    IMG_TYPE=t1
    anat_fname=sub-${sub}_ses-${SES_NO}_T1w.nii.gz
    ls ${BASEDIR}/${IMG_DIR}/*${sub}*${IMG_TYPE}.nii.gz
    #echo "Copy to ${SUBJDIR}/${BIDS_SESSION}/anat/$anat_fname"

    # mkdir & link
    mkdir -p ${SUBJDIR}/${BIDS_SESSION}/anat
    unlink ${SUBJDIR}/${BIDS_SESSION}/anat/$anat_fname
    ln -s ${BASEDIR}/${IMG_DIR}/*${sub}*${IMG_TYPE}.nii.gz ${SUBJDIR}/${BIDS_SESSION}/anat/$anat_fname

    #################################
    ###### Resting ##########
    #################################
    TASK=rest
    TASKNAME=rest
    IMG_TYPE=${TASK}
    RUN_NO=1
    rest_fname=sub-${sub}_ses-${SES_NO}_task-${TASKNAME}_run-${RUN_NO}_bold.nii.gz
    ls ${BASEDIR}/${IMG_TYPE}/*${sub}*${TASK}.nii.gz
    #echo "copy to ${SUBJDIR}/${BIDS_SESSION}/func/$rest_fname"
    # mkdir & link
    mkdir -p ${SUBJDIR}/${BIDS_SESSION}/func
    unlink ${SUBJDIR}/${BIDS_SESSION}/func/$rest_fname
    ln -s ${BASEDIR}/${IMG_TYPE}/*${sub}*${TASK}.nii.gz ${SUBJDIR}/${BIDS_SESSION}/func/$rest_fname

done
# Remove broken links
for f in `find -L ${OUTPUTDIR} -maxdepth 4 -type l`; do unlink $f; done


# ###################################
# # BIDSify TMS Tasks
# ###################################

for sub in ${SUBJID[@]}
do

    echo "$sub"
    SUBJDIR=${OUTPUTDIR}/sub-${sub}

	DAY_NO=2
	SES_NO=d${DAY_NO}
	BIDS_SESSION=ses-${SES_NO}

	IMG_TYPE=tmsfmri
	TASK=singlepulse
	mkdir -p ${SUBJDIR}/${BIDS_SESSION}/
	mkdir -p ${SUBJDIR}/${BIDS_SESSION}/func

	FR_TMS_SP_SITES=(L_aMFG L_Fp L_IFG_Anat L_IFG_BLA L_pMFG L_VmFp R_aMFG R_FEF R_Fp R_IFG_Anat R_IFG_BLA R_IFJ R_IPL R_M1 R_pMFG R_preSMA R_VmFp)

	TO_TMS_SP_SITES=(LaMFG LFp LIFGAnat LIFGBLA LpMFG LVmFp RaMFG RFEF RFp RIFGAnat RIFGBLA RIFJ RIPL RM1 RpMFG RpreSMA RVmFp)

	RUN_NO=1
	for (( siteno=0; siteno<${#FR_TMS_SP_SITES[@]}; siteno++))
	do
	    #echo "From ${FR_TMS_SP_SITES[$siteno]} to ${TO_TMS_SP_SITES[$siteno]}"
		ls ${BASEDIR}/${IMG_TYPE}/${FR_TMS_SP_SITES[$siteno]}/*${sub}_*${FR_TMS_SP_SITES[$siteno]}*.nii.gz
		func_fname=sub-${sub}_ses-${SES_NO}_task-${TASK}${TO_TMS_SP_SITES[$siteno]}_run-${RUN_NO}_bold.nii.gz
		echo "Copying to ${SUBJDIR}/${BIDS_SESSION}/func/$func_fname"
		unlink ${SUBJDIR}/${BIDS_SESSION}/func/$func_fname
		ln -s ${BASEDIR}/${IMG_TYPE}/${FR_TMS_SP_SITES[$siteno]}/*${sub}_*${FR_TMS_SP_SITES[$siteno]}*.nii.gz ${SUBJDIR}/${BIDS_SESSION}/func/$func_fname
	done
done
# Remove broken links
for f in `find -L ${OUTPUTDIR} -maxdepth 4 -type l`; do unlink $f; done


###################################
# BIDSify TMS Tasks
# colorID_r1  colorID_r2  emoconflict  MID  MSIT
###################################

for sub in ${SUBJID[@]}
do
	
    echo "$sub"
    SUBJDIR=${OUTPUTDIR}/sub-${sub}
	
	DAY_NO=1
	SES_NO=d${DAY_NO}
	BIDS_SESSION=ses-${SES_NO}
	
	IMG_TYPE=mid
	TASKNAME=monetaryincentivedelay
	TASKS=(colorID_r1 colorID_r2 emoconflict mid msit)
	TASK=mid
	mkdir -p ${SUBJDIR}/${BIDS_SESSION}/
	mkdir -p ${SUBJDIR}/${BIDS_SESSION}/func

	RUN_NO=1
    func_fname=sub-${sub}_ses-${SES_NO}_task-${TASKNAME}_run-${RUN_NO}_bold.nii.gz
    ls ${BASEDIR}/${IMG_TYPE}/*${sub}*${TASK}.nii.gz
	echo "Copying to ${SUBJDIR}/${BIDS_SESSION}/func/$func_fname"

	unlink ${SUBJDIR}/${BIDS_SESSION}/func/$func_fname
	ln -s ${BASEDIR}/${IMG_TYPE}/*${sub}*${TASK}.nii.gz ${SUBJDIR}/${BIDS_SESSION}/func/$func_fname
	
	# for (( taskno=0; taskno<${#TASKS[@]}; taskno++))
	# do
	#     func_fname=sub-${sub}_ses-${SES_NO}_task-${TASKNAME}_run-${RUN_NO}_bold.nii.gz
	#     ls ${BASEDIR}/${IMG_TYPE}/*${sub}*${TASK}.nii.gz
	# 	echo "Copying to ${SUBJDIR}/${BIDS_SESSION}/func/$func_fname"
	#
	# 	ln -s ${BASEDIR}/${IMG_TYPE}/*${sub}*${TASK}.nii.gz ${SUBJDIR}/${BIDS_SESSION}/func/$func_fname
	# done
done
# Remove broken links
for f in `find -L ${OUTPUTDIR} -maxdepth 4 -type l`; do unlink $f; done
echo "find -L ${OUTPUTDIR} -maxdepth 4  -type d -empty -print"
find -L ${OUTPUTDIR} -maxdepth 4  -type d -empty -delete

