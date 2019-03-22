#!/bin/sh
################################
# Input and output directories #
################################
BASEDIR=${PI_SCRATCH}/COMET/CausalConnectome
BIDSDIR=${BASEDIR}/BIDS
OUTPUTDIR=${BASEDIR}/derivatives/fmriprep
WORKDIR=${BASEDIR}/work/fmriprep
##########################
# Singularity containers #
##########################
# FPREP_IMG=${PI_HOME}/singularity_images/poldracklab_fmriprep_1.2.5-2018-12-04-2ef6b23ede2a.img
FPREP_IMG=${PI_HOME}/singularity_images/bids-fmriprep-1.2.3_latest.sif
##########################
# usage: fmriprep [-h] [--version] [--skip_bids_validation]
#                [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
#                [-t TASK_ID] [--echo-idx ECHO_IDX] [--nthreads NTHREADS]
#                [--omp-nthreads OMP_NTHREADS] [--mem_mb MEM_MB] [--low-mem]
#                [--use-plugin USE_PLUGIN] [--anat-only] [--boilerplate]
#                [--ignore-aroma-denoising-errors] [-v] [--debug]
#                [--ignore {fieldmaps,slicetiming,sbref} [{fieldmaps,slicetiming,sbref} ...]]
#                [--longitudinal] [--t2s-coreg] [--bold2t1w-dof {6,9,12}]
#                [--output-space {T1w,template,fsnative,fsaverage,fsaverage6,fsaverage5} [{T1w,template,fsnative,fsaverage,fsaverage6,fsaverage5} ...]]
#                [--force-bbr] [--force-no-bbr]
#                [--template {MNI152NLin2009cAsym}]
#                [--output-grid-reference OUTPUT_GRID_REFERENCE]
#                [--template-resampling-grid TEMPLATE_RESAMPLING_GRID]
#                [--medial-surface-nan] [--use-aroma]
#                [--aroma-melodic-dimensionality AROMA_MELODIC_DIMENSIONALITY]
#                [--skull-strip-template {OASIS,NKI}]
#                [--skull-strip-fixed-seed] [--fmap-bspline] [--fmap-no-demean]
#                [--use-syn-sdc] [--force-syn] [--fs-license-file PATH]
#                [--no-submm-recon] [--cifti-output | --fs-no-reconall]
#                [-w WORK_DIR] [--resource-monitor] [--reports-only]
#                [--run-uuid RUN_UUID] [--write-graph] [--stop-on-first-crash]
#                [--notrack] [--sloppy]
#                bids_dir output_dir {participant}
#########################
# Recommended resource constraints
# one compute node per subject, 16 cores with 64Gb
#########################
# Preprocessing Options #
#########################
# SLURM_ARRAY_TASK_ID=1
SUBID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../00-bidsify/subjects_tms.txt)
TASKID=(rest)
TASK_TMS=(singlepulseLaMFG)
# singlepulseLaMFG singlepulseLFp singlepulseLIFGAnat singlepulseLIFGBLA singlepulseLpMFG singlepulseLVmFp
# singlepulseRaMFG singlepulseRFp singlepulseRIFGAnat singlepulseRIFGBLA singlepulseRpMFG singlepulseRVmFp
# singlepulseRFEF singlepulseRIPL singlepulseRIFJ singlepulseRM1 singlepulseRpreSMA
N_CPUS=16
MEM_GB=64
MEM_MB=64000
# OUTPUT_SPACE=T1w template
# --template
MNI152_SYM=/home/users/manjarin/ANALYSIS/pipelines/sherlock-preproc/2017-preproc-tms-fmri/configs/resources/MNI152_T1_2mm_brain.nii.gz
##
IGNORE_OPTS=(slicetiming fieldmaps) 
OUTPUT_SPACE=(T1w template)
BOLD2T1DOF=9 # can be 6,9,12. --bold2t1w-dof
OUTPUT_GRID_REFERENCE=${HOME}/ANALYSIS/pipelines/sherlock-preproc/2017-preproc-tms-fmri/configs/resources/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz
# To be added or considered
# --no-skull-strip-ants (not available in v1.0.7)
#--fs-no-reconall --no-freesurfer
# --fs-no-reconall
# --skull-strip-ants 
# --output-grid-reference (Specify craddock or yeo compatible one) –template-resampling-grid
# --use-aroma
# --ignore-aroma-denoising-errors
# --force-no-bbr
#########################
#     Preliminaries     #
#########################
export SINGULARITY_BIND="/scratch,/share"
PYTHONPATH=""
#########################

for sub in ${SUBID}
do
    subname=sub-${sub}
    echo "$subname"
    
    # # OUTPUT SPACE for FMRI
    echo "${SCRATCHDIR}"
    echo "${WORKDIR}"
    
    echo "singularity run ${FPREP_IMG} ${BIDSDIR} ${OUTPUTDIR} participant --boilerplate \
    --skip_bids_validation \
    --participant_label ${sub} \
    --task-id ${TASK_TMS} \
    --bold2t1w-dof ${BOLD2T1DOF} \
    --ignore-aroma-denoising-errors \
    --output-space "${OUTPUT_SPACE}" \
    --write-graph \
    --template-resampling-grid "${OUTPUT_GRID_REFERENCE}" \
    --fs-no-reconall \
    --ignore ${IGNORE_OPTS} fieldmaps \
    -w ${WORKDIR} --n_cpus ${N_CPUS} --omp-nthreads 4 --mem_mb ${MEM_MB}"
    
    singularity run ${FPREP_IMG}  ${BIDSDIR} ${OUTPUTDIR} participant \
	--participant_label ${sub} \
        --bold2t1w-dof ${BOLD2T1DOF} \
        --force-bbr \
        --use-aroma --aroma-melodic-dimensionality -75 \
        --ignore-aroma-denoising-errors \
        --output-space "${OUTPUT_SPACE}" \
        --write-graph \
        --template-resampling-grid "${OUTPUT_GRID_REFERENCE}" \
        --ignore "${IGNORE_OPTS}" \
        --fs-license-file /share/software/user/open/freesurfer/6.0.0/license.txt \
        --longitudinal \
        -w ${WORKDIR} --n_cpus ${N_CPUS} --omp-nthreads 4 --mem_mb ${MEM_MB}

    # for task in ${TASKID[@]}
    # do
    #     echo "Task specific ${task}"
    #     # # OUTPUT SPACE for FMRI
    #     singularity run ${FPREP_IMG} ${BIDSDIR} ${OUTPUTDIR} \
    #         participant --participant_label ${sub} \
    #         -t ${task} --no-freesurfer \
    #         --ignore "${IGNORE_OPTS}" --no-skull-strip-ants \
    #         --ignore ${IGNORE_OPTS} \ 
    #         --ignore-aroma-denoising-errors \
    #         --output-space "${OUTPUT_SPACE}" \
    #    --output-grid-reference ${OUTPUT_GRID_REFERENCE} \   
    #         --use-aroma --write-graph \
    #         -w ${WORKDIR} --n_cpus ${N_CPUS} --mem_mb ${MEM_MB}
    # done
done
