#!/bin/ksh
#PJM -L "node=12:noncont"
#PJM -L elapse=12:00:00
#PJM -L rscgrp=small
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=576"
#PJM -j
#PJM -g sum
#PJM -N rsm_timcom
#PJM -o %j.log
#PJM -e %j.err


set -x

export RMPI=576
OMP=1
source /users/xa09/sample/setup_mpi+omp.fx1000  $OMP

export RSMDISK=`cd ../..; pwd`
export SDATE=2023010100             # initial date
export RESTRHR=0                    # restart hour (zero means NOT restart)
export FCSTHR=24                   # forecast hour
export RUNDIR=$RSMDISK/wrk/${SDATE}   # running directory

RSMPATH=${RSMDISK}/run; cd ${RSMPATH} || exit 8
 sed -i "s/^export SDATE=.*$/export SDATE=$SDATE/" /nwpr/gfs/xb157/CWBSUMc2/components/rsm/run/exp/configure_5km
 sed -i "s/^export FCSTHR=.*$/export FCSTHR=$FCSTHR/" /nwpr/gfs/xb157/CWBSUMc2/components/rsm/run/exp/configure_5km
 RSM_FCST=${RSMPATH}/exe/rsm_cpl.x
. ${RSMPATH}/run_2cpl.sh ${RSMDISK} ${SDATE} ${RESTRHR} ${FCSTHR} ${RUNDIR} ||  exit 8

#
#RSM
#/usr/bin/time -p mpiexec  -stdin rfcstparm.all --of-proc stdout.mpmd -n $RMPI ${RSM_FCST}
# data from GFS (file 11,12,18) are big-endian
/usr/bin/time -p mpiexec  -stdin rfcstparm.all -n $RMPI ${RSM_FCST} -Wl,-T11,-T12,-T18

#/usr/bin/time -p mpiexec -n $RMPI ${RSM_FCST} -Wl,-T
if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!"
fi



