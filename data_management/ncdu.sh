#!/bin/sh

source /etc/profile.d/lsf.sh 
source /etc/profile.d/modules.sh

#source /etc/profile.d/00-modulepath.sh
#source /etc/profile.d/z00_lmod.sh

bsub -sla miket_sc < /data/talkowski/iwong/src/scripts/ncdu.lsf



