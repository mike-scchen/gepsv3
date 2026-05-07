#!/bin/ksh
#PJM -L "node=1:noncont"
#PJM -L elapse=03:00:00
#PJM -L rscgrp=small
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=1"
#PJM -g sum
#PJM -j
#PJM -g sum
#PJM -N rsm_post
#PJM -o %j.log
#PJM -e %j.err

#set -x

OMP=1
source /users/xa09/sample/setup_mpi+omp.fx1000  $OMP

export RSMDISK=`cd ../..; pwd`
export SDATE=2023010100             # initial date
export RUNDIR=/nwpr/gfs/xb157/CWBSUMc2_cpl/rsm/wrk/$SDATE
export RESTRHR=0                    # restart hour (zero means NOT restart)
export FCSTHR=24                  # forecast hour

RSMPATH=${RSMDISK}/run; cd ${RSMPATH} || exit 8
. ${RSMPATH}/configure ||  exit 8
echo $SHSDIR  
echo $RUNDIR

# cd to the RSM output directory
cd $RUNDIR

# run RSM post processor
for ihr in {${RESTRHR}..${FCSTHR}..${PRTHOUR}};
do
  echo "processing fhour: $ihr"
  /usr/bin/time ${SHSDIR}/run_rpgb.sh  $ihr
  wait
done
echo "Finish RSM rpgb"
