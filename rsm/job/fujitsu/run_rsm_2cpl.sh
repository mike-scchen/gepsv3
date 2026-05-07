#!/bin/bash
#PJM -L "node=36:noncont"
#PJM -L rscgrp=large
#PJM -x PJM_CACHE_MODE=4
#PJM -L elapse=15:00:00
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=1728"
#PJM -j
#PJM -N rsm_230101
#PJM -o %j.log
#PJM -e %j.err

OMP=1
source /users/xa09/sample/setup_mpi+omp.fx1000 $OMP

export n_proc_gfs=0
export n_proc_rsm=1152
export n_proc_gocn=0
export n_proc_rocn=576
export MPI=1728

export RSMDISK=/nwpr/gfs/xb157/CWBSUMc2/4cpl/op_work/rsm
# initial date
export SDATE=2023010100   
# restart hour (zero means NOT restart)
export RESTRHR=0  
# forecast hour
export FCSTHR=1080
# running directory
export RUNDIR=$RSMDISK/wrk/${SDATE}_2cpl 
export RSMPATH=$RSMDISK/run

dtgg=`echo ${SDATE} | cut -c3-10`
glb_path=/nwpr/gfs/xb157/CWBSUMc2/4cpl/op_work
TIMCOMr_bcneed=${glb_path}/${dtgg}.ufs/TIMCOM_glb
GFSWRK=${glb_path}/${dtgg}.ufs/wrk
NAMP=${glb_path}/${dtgg}.ufs/rsm
ln -fs ${TIMCOMr_bcneed} ${RUNDIR}/
ln -fs ${GFSWRK}/rsm_data* $RUNDIR/
ln -fs ${GFSWRK}/rsm_idate_${SDATE} $RUNDIR/rsm_idate
cd ${RUNDIR}
rename "_${SDATE}00" "" rsm_data_*.f*
ln -fs  ${GFSWRK}/rsm_xlat_${SDATE}00.txt $RUNDIR/g2r_xlat.txt
ln -fs  ${GFSWRK}/rsm_xlon_${SDATE}00.txt $RUNDIR/g2r_xlon.txt
ln -fs ${NAMP}/tai_grid.nc ${RUNDIR}/.
ln -fs ${NAMP}/tai_lbc.nc ${RUNDIR}/.
ln -fs ${NAMP}/tai000.nc ${RUNDIR}/.

cp -p ${NAMP}/namelist000.inp_tai ${RUNDIR}/.
cp -p ${NAMP}/namelist000.rof ${RUNDIR}/.
cp -p ${NAMP}/namelist000.run ${RUNDIR}/.

rm -rf ${RUNDIR}/tai_runtime.log
cd ${RSMPATH} 
. ${RSMPATH}/run_2cpl.sh ${RSMDISK} ${SDATE} ${RESTRHR} ${FCSTHR} ${RUNDIR}

ln -sf /nwpr/gfs/xb157/CWBSUMc2/4cpl/components/driver/TCoTIMCOM_2cpl $RSMPATH/exe/rsm_2cpl.x
FCT_MODEL=$RSMPATH/exe/rsm_2cpl.x

/usr/bin/time -p mpiexec -stdin ./rfcstparm.all -n ${MPI} ${FCT_MODEL} 000 -Wl,-T11,-T12,-T18

if [ $? != 0 ] ; then
  echo "error occured: fct model fail !!"
fi

echo "Ending at: " $(date)
