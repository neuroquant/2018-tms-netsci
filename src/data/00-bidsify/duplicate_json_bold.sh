#!/bin/sh

TMS_TASKS=(LFp LIFGAnat LIFGBLA LpMFG LVmFp RaMFG RFEF RFp RIFGAnat RIFGBLA RIFJ RIPL RM1 RpMFG RpreSMA RVmFp)

for (( siteno=0; siteno<${#TMS_TASKS[@]}; siteno++))
do
   echo "Duplicating task-singlepulseLaMFG for site ${TMS_TASKS[$siteno]}"
   cp task-singlepulseLaMFG_bold.json task-singlepulse${TMS_TASKS[$siteno]}_bold.json
done
