#!/bin/sh

OUTPUTDIR=$PI_SCRATCH/COMET/CausalConnectome/derivatives/fmriprep-fsl/fmriprep

echo "find ${OUTPUTDIR}/sub-NTHC*/figures/* -newermt 'Jan 25 00:00' \!  -newermt 'Jan 28 00:00' -print"
find ${OUTPUTDIR}/sub-TEHC*/figures/* -newermt 'Jan 25 00:00' \!  -newermt 'Jan 28 00:00' -ls
echo "find ${OUTPUTDIR}/sub-TEHC*/figures/* -newermt 'Jan 25 00:00' \!  -newermt 'Jan 28 00:00' -print"
