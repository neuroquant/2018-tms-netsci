#!/bin/sh

sbatch -N 1 -n 2 -p aetkin,owners,normal -t 4:00:00 --begin=now gdrive_backup.sh
sbatch -N 1 -n 2 -p aetkin,owners,normal -t 6:00:00 --begin=now+4hour gdrive_backup.sh
