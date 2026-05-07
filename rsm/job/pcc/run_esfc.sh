#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime=00:20:00
#PBS -q dc20190047
#PBS -P MST107203 
#PBS -N extrsfc
#PBS -j oe
#PBS -ko

set -x

module load intel/2018_u1

 export  I_MPI_PIN_CELL=core

PWDIR=/work1/yj8taiwania/model/tco639l72
#################################################################

#SDATE=19082200
echo $SDATE
echo $RSMVER

 GFSWRK=/project/pm25/pub/team2/msm/wrk_${RSMVER}_${SDATE}
cd $GFSWRK

# send data to datamv05
Pathto=rsmop@datamv05.cwb.gov.tw:/nwpr/rsm/rsmop/data2/RSM5km_archive/cwb_TCoRSM_${SDATE}/
rsync -avhz ${GFSWRK}/r_pgb* $Pathto
chmod 644 ${GFSWRK}/r_ltn*
rsync -avhz ${GFSWRK}/r_ltn* $Pathto

#chmod -R 755 ${GFSWRK}/*

# extract sfc field from grib
/usr/bin/time ${GFSWRK}/run_extrsfc.sh
echo "Finish RSM post processor - extract sfc field"

# use sfc field to plot
mkdir -p $GFSWRK/plot
cp  ${PWDIR}/plot/*  $GFSWRK/plot/
#/usr/bin/time ${GFSWRK}/plot/drawsfc.sh   $GFSWRK
#/usr/bin/time ${GFSWRK}/plot/drawatm.sh   $GFSWRK  $SDATE
#echo "Finish RSM 5km surface field plotting"

/usr/bin/time ${PWDIR}/plot/files_TWN/drawrain.sh $GFSWRK $SDATE
echo "Finish RSM 5km surface field plotting"


#chmod wind turbines output
#chmod 644 /project/pm25/pub/team2/msm/wrk_${RSMVER}_${SDATE}/power.txt
#chmod 644 /project/pm25/pub/team2/msm/wrk_${RSMVER}_${SDATE}/windspeed.txt

#rsync figure to dm05
rsync -avhz /home/rsmop/TCO_RSM_branch/RSM_archive/FIG_GFS_RSM/fig_${SDATE} rsmop@datamv05.cwb.gov.tw:/nwpr/rsm/rsmop/fig_RSM_MSM/TCo383+RSM5km_OP_FIGURE/TCo383_RSM5km/
#sh /home/yj8taiwania/data/sh/power/get_data_power.sh
