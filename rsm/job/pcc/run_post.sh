#!/bin/ksh
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime=00:20:00
#PBS -q research
#PBS -N rpgb
#PBS -j oe
#PBS -ko

set -x

module load intel/18.0.3

export  I_MPI_PIN_CELL=core

export RSMDISK=${PBS_O_WORKDIR}/../..
export SDATE=2021060200             # initial date
export RUNDIR=$RSMDISK/wrk/$SDATE   # running directory

RSMPATH=${RSMDISK}/run; cd ${RSMPATH} || exit 8
. ${RSMPATH}/configure ||  exit 8
echo $SHSDIR  
echo $RUNDIR

# cd to the RSM output directory
cd $RUNDIR

# run RSM post processor
for ihr in {0..${FCSTHR}..${PRTHOUR}};
do
  echo "processing fhour: $ihr"
  /usr/bin/time ${SHSDIR}/run_rpgb.sh  $ihr
  wait
done
echo "Finish RSM rpgb"


# use loop above, or submit array job
#/usr/bin/time ${SHSDIR}/run_rpgb.sh ${PBS_ARRAY_INDEX}

