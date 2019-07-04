#!/bin/sh
################################
# Input and output directories #
################################
BASEDIR=${PI_SCRATCH}/COMET/CausalConnectome
BIDSDIR=${BASEDIR}/derivatives/fmriprep-fsl/fmriprep
OUTPUTDIR=${BASEDIR}/derivatives/denoiser
##########################
DENOISE="python denoiser/run_denoise.py"
##########################

export WORKON_HOME=${PI_HOME}/software/envs
export VIRTUALENVWRAPPER_PYTHON=/share/software/user/open/python/3.6.1/bin/python3
alias source_venv='source ~/.local/bin/virtualenvwrapper.sh'

source_venv
workon ni-denoise

# TASKS=(rest multisourceinterference singlepulseLaMFG singlepulseLFp singlepulseLIFGAnat singlepulseLIFGBLA singlepulseLpMFG singlepulseLVmFp singlepulseRaMFG singlepulseRFp singlepulseRIFGAnat singlepulseRIFGBLA singlepulseRpMFG singlepulseRVmFp singlepulseRFEF singlepulseRIPL singlepulseRIFJ singlepulseRM1 singlepulseRpreSMA)

TASKS=(singlepulseLIFGAnat singlepulseLIFGBLA singlepulseLpMFG singlepulseLVmFp singlepulseRaMFG singlepulseRFp singlepulseRIFGAnat singlepulseRIFGBLA singlepulseRpMFG singlepulseRVmFp singlepulseRFEF singlepulseRIPL singlepulseRIFJ singlepulseRM1 singlepulseRpreSMA)


CONFOUND_NAMES="csf	white_matter    std_dvars   framewise_displacement   t_comp_cor_00   t_comp_cor_01   t_comp_cor_02   t_comp_cor_03   t_comp_cor_04   t_comp_cor_05  a_comp_cor_00	a_comp_cor_01	a_comp_cor_02	a_comp_cor_03	a_comp_cor_04	a_comp_cor_05	cosine00	cosine01	cosine02	cosine03	cosine04 trans_x	trans_y	trans_z	rot_x	rot_y	rot_z"

for (( taskno=0; taskno<${#TASKS[@]}; taskno++))
do

    TMS_TASK=${TASKS[$taskno]}
    echo ${TMS_TASK}
    
    # thefile=$(find ${BIDSDIR}/${SUBJECT} -maxdepth 4 -type f -name "*T1w*preproc_bold*" -printf "%T@ %p " | sort -k 1n | tail -n 1)
    thefile=$(find ${BIDSDIR}/${SUBJECT} -maxdepth 4 -type f -name "*${TMS_TASK}*T1w*preproc_bold.nii.gz*" | sort -k 1n | tail -n 20)

    # echo "The list of files:
    # $thefile"

    for (( fileno=0; fileno<${#thefile[@]}; fileno++))
    do
        IMG_FILE=${thefile[$fileno]}
        CONFOUND_BASE=$(echo ${IMG_FILE} | sed -e "s?.*func/??g" | sed -e "s?space-T1w.*??g")
       CONFOUND_FILE=$(find ${BIDSDIR}/${SUBJECT} -maxdepth 4 -type f -name "*${CONFOUND_BASE}*confound*" | head -n 1)
      echo "IMG: ${IMG_FILE}, CONFOUND: ${CONFOUND_FILE}"
    #    echo "Running: ${DENOISE} ${IMG_FILE} ${CONFOUND_FILE} ${OUTPUTDIR}"
        mkdir -p ${OUTPUTDIR}/${SUBJECT}
        $DENOISE ${IMG_FILE} ${CONFOUND_FILE} ${OUTPUTDIR}/${SUBJECT} --col_names ${CONFOUND_NAMES} --hp_filter .009
    done
done
