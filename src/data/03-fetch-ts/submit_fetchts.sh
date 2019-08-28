#!/bin/sh

for subno in $(seq 45 53)
do
    sbatch --mem=32GB -t 1:00:00 --array=$subno -p aetkin,owners,normal run_fetchts.sh
done
