#!/bin/sh

# SUBID=(6) # 13, 20, 29, 6 16 17 35 37  
# Sept 17 we had to run 16, 35 and 37. The 206x seem to be particularly failing due to mixed July runs?

#for sub in ${SUBID[@]}
for sub in $(seq 1 14)
do
   sbatch -p aetkin,normal,owners --array=$sub -N 1 -n 16 --mem=64000 -t 8:00:00 --begin=now run_fmriprep.sh
   # now+${sub}minute
   #sbatch -p aetkin,owners,normal --array=$sub -N 1 -n 16 --mem=64000 -t 4:00:00 --begin=now run_fmriprep.sh

   # 6 hours was enough to complete freesurfer and T1w processing but not the rest. 
   # sbatch -p aetkin,owners,normal --array=$sub -N 1 -n 16 --mem=64000 -t 4:00:00 --begin=now+9hour run_fmriprep.sh
done
