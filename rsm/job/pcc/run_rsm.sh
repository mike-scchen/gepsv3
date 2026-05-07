#!/bin/bash
#PBS -l select=9:ncpus=36:mpiprocs=36
#PBS -q research
#PBS -l walltime=06:00:00
#PBS -N RSM
#PBS -j oe
#PBS -ko

set -x

module load intel/18.0.3

machine=pcc

export RMPI=324
OMP=1

export  I_MPI_PIN_CELL=core
# export I_MPI_PIN_DOMAIN=omp
# export KMP_AFFINITY=granularity=thread
export OMP_NUM_THREADS=$OMP

export RSMDISK=${PBS_O_WORKDIR}/../..
export SDATE=2021060200             # initial date
export RESTRHR=0                    # restart hour (zero means NOT restart)
export FCSTHR=24                    # forecast hour
export RUNDIR=$RSMDISK/wrk/$SDATE   # running directory

RSMPATH=${RSMDISK}/run; cd ${RSMPATH} || exit 8
RSM_FCST=${RSMPATH}/exe/rsm.x
. ${RSMPATH}/run.sh ${RSMDISK} ${SDATE} ${RESTRHR} ${FCSTHR} ${RUNDIR} ||  exit 8

#RSM
/usr/bin/time -p mpiexec  -prepend-rank -n  $RMPI ${RSM_FCST}
#/usr/bin/time -p mpiexec -s all -n $RMPI ${RSM_FCST}
if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!"
fi

