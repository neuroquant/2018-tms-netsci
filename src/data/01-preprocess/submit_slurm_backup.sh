#!/bin/sh

sbatch -N 1 -n 2 -p aetkin,owners,normal -t 8:00:00 --begin=now gdrive_backup.sh
sbatch -N 1 -n 2 -p aetkin,owners,normal -t 6:00:00 --begin=now+10hour gdrive_backup.sh
