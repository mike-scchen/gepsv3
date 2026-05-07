#!/bin/bash
#SBATCH -A MST107192
#SBATCH -p ct560
#SBATCH -t 4-00:00:00
#SBATCH -N 9
#SBATCH --ntasks-per-node=48
#SBATCH --switches=1
#SBATCH -J 20200201_rsm_cpl
#SBATCH -o %j.log
#SBATCH -e %j.err

set -x

module purge
module load compiler/intel/2020u4
module load IntelMPI/2020
module load netcdf-4.8.0-NC4-intel2020-impi
module load pnetcdf-1.8.1-intel2020-impi
module load rcec/ncl/6.6.2
module list

machine=pcc

export RMPI=432
OMP=1

export  I_MPI_PIN_CELL=core
# export I_MPI_PIN_DOMAIN=omp
# export KMP_AFFINITY=granularity=thread
export OMP_NUM_THREADS=$OMP

export RSMDISK=`cd ../..; pwd`
export SDATE=2020030100             # initial date
export RESTRHR=0                    # restart hour (zero means NOT restart)
export FCSTHR=8760                  # forecast hour
export RUNDIR=$RSMDISK/wrk/${SDATE}_dm5d5   # running directory

RSMPATH=${RSMDISK}/run; cd ${RSMPATH} || exit 8
RSM_FCST=${RSMPATH}/exe/rsm_timcom.x
. ${RSMPATH}/run.sh ${RSMDISK} ${SDATE} ${RESTRHR} ${FCSTHR} ${RUNDIR} ||  exit 8

#RSM
/usr/bin/time -p mpiexec  -prepend-rank -n  $RMPI ${RSM_FCST} 000
#/usr/bin/time -p mpiexec -s all -n $RMPI ${RSM_FCST}
if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!"
fi

