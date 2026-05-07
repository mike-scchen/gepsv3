#!/bin/bash
#PJM -L "node=36:noncont"
#PJM -L rscgrp=large
#PJM -x PJM_CACHE_MODE=4
#PJM -L elapse=15:00:00
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=1728"
#PJM -j
#PJM -N rsm_230423
#PJM -o %j.log
#PJM -e %j.err

OMP=1
source /users/xa09/sample/setup_mpi+omp.fx1000 $OMP

export n_proc_gfs=0
export n_proc_rsm=1152
export n_proc_gocn=0
export n_proc_rocn=576
export MPI=1728

export RSMDISK=/nwpr/gfs/xb124/RSM_2cpl
# initial date
export SDATE=2023042300   
# restart hour (zero means NOT restart)
export RESTRHR=0  
# forecast hour
export FCSTHR=1080 
# running directory
export RUNDIR=$RSMDISK/wrk/$SDATE 
export RSMPATH=$RSMDISK/run
cd ${RSMPATH} 
. ${RSMPATH}/run.sh ${RSMDISK} ${SDATE} ${RESTRHR} ${FCSTHR} ${RUNDIR}

FCT_MODEL=$RSMPATH/exe/rsm_4cpl.x
/usr/bin/time -p mpiexec -stdin ./rfcstparm.all -n ${MPI} ${FCT_MODEL} 000 -Wl,-T11,-T12,-T18

if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!"
fi

echo "Ending at: " $(date)
