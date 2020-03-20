#!/bin/sh

sbatch -N 1 -n 2 -p normal -t 8:00:00 gdrive_fetch.sh
sbatch -N 1 -n 2 -p aetkin,normal -t 6:00:00 --begin=now+10hour gdrive_backup.sh
