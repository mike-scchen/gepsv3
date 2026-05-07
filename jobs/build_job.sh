#/bin/bash
set -x
kshpath=`pwd `
. ${kshpath}/00_setDate.sh

export caldtg="/data/common/gfs/scripts/Caldtg.ksh"
export dtg=$(echo ${op_date}|cut -c 3-10)
export fgdtg=$(${caldtg} ${dtg} -6)
export ocn_restart=.false.
export commonpath=/data/common/gfs/GEPSv3_lib
export MDIR=${kshpath}/${dtg}_${suffix}
export DMSDIR=${kshpath}/dmsdb_${suffix}
export GFSDIR2=${kshpath}/../tco639l72
export GFSDIR=${GFSDIR2}
export RSMDIR=${kshpath}/../rsm
 cd ${kshpath}

export GFSWRK=${MDIR}.ufs/wrk
export GFSOUT=${MDIR}.ufs/TCo$JCAP
export RSMOUT=${MDIR}.ufs/rsm
export GLBOUT=${MDIR}.ufs/TIMCOM_glb
export TAIOUT=${MDIR}.ufs/TIMCOM_tai
export CICEOUT=${MDIR}.ufs/CICE

mkdir -p ${MDIR}.ufs
mkdir -p ${DMSDIR}.ufs/in
mkdir -p ${DMSDIR}.ufs/bck
mkdir -p ${GFSWRK}
mkdir -p ${GFSOUT}
mkdir -p ${RSMOUT}
mkdir -p ${GLBOUT}
mkdir -p ${TAIOUT}
mkdir -p ${CICEOUT}/restart
######################################################################
if [ ${mach} = gpu -a ${struc} = 2cpl ] ;then
  if [ $JCAP = 383  ] ; then
    n_proc_gfs=8
    n_proc_rsm=0
    n_proc_gocn=8
    n_proc_rocn=0
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
    ln -sf ${commonpath}/data/rmp_tco2timcom_xnew.nc ${GFSWRK}/rmp_tco2timcom.nc
    ln -sf ${commonpath}/data/rmp_timcom2tco_xnew.nc ${GFSWRK}/rmp_timcom2tco.nc
  elif [ $JCAP = 199  ] ; then
    n_proc_gfs=8
    n_proc_rsm=0
    n_proc_gocn=8
    n_proc_rocn=0
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
    ln -sf ${commonpath}/data/rmp_tco199_to_timcom_1000.nc ${GFSWRK}/rmp_tco2timcom.nc
    ln -sf ${commonpath}/data/rmp_timcom_to_tco199_1000.nc ${GFSWRK}/rmp_timcom2tco.nc
  fi

fi

if [ ${mach} = fx1000 -a ${struc} = 2cpl ] ;then
  if [ $JCAP = 639  ] ; then
    n_proc_gfs=768
    n_proc_rsm=0
    n_proc_gocn=768
    n_proc_rocn=0
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
  elif [ $JCAP = 383  ] ; then
    n_proc_gfs=256
    n_proc_rsm=0
    n_proc_gocn=240
    n_proc_rocn=0
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
    ln -sf ${commonpath}/data/rmp_tco2timcom_xnew.nc ${GFSWRK}/rmp_tco2timcom.nc
    ln -sf ${commonpath}/data/rmp_timcom2tco_xnew.nc ${GFSWRK}/rmp_timcom2tco.nc
  elif [ $JCAP = 199  ] ; then
    n_proc_gfs=600
    n_proc_rsm=0
    n_proc_gocn=1920
    n_proc_rocn=0
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
    ln -sf ${commonpath}/data/rmp_tco199_to_timcom_1000.nc ${GFSWRK}/rmp_tco2timcom.nc
    ln -sf ${commonpath}/data/rmp_timcom_to_tco199_1000.nc ${GFSWRK}/rmp_timcom2tco.nc
  fi  
fi

if [ ${mach} = fx1000 -a ${struc} = 2cpl_CICE ] ;then
  if [ $JCAP = 639  ] ; then
    n_proc_gfs=768
    n_proc_rsm=0
    n_proc_gocn=768
    n_proc_rocn=0
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
  elif [ $JCAP = 383  ] ; then
    n_proc_gfs=384
    n_proc_rsm=0
    n_proc_gocn=384
    n_proc_rocn=0
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
    ln -sf ${commonpath}/data/rmp_tco2timcom_xnew.nc ${GFSWRK}/rmp_tco2timcom.nc
    ln -sf ${commonpath}/data/rmp_timcom2tco_xnew.nc ${GFSWRK}/rmp_timcom2tco.nc
  elif [ $JCAP = 199  ] ; then
    n_proc_gfs=600
    n_proc_rsm=0
    n_proc_gocn=384 #1152
    n_proc_rocn=0
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
    ln -sf ${commonpath}/data/rmp_tco199_to_timcom_1000.nc ${GFSWRK}/rmp_tco2timcom.nc
    ln -sf ${commonpath}/data/rmp_timcom_to_tco199_1000.nc ${GFSWRK}/rmp_timcom2tco.nc
  fi
fi


if [ ${mach} = fx1000 -a ${struc} = 4cpl ] ;then
  if [ $JCAP = 639  ] ; then
    n_proc_gfs=768
    n_proc_rsm=1152
    n_proc_gocn=768
    n_proc_rocn=576
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
  elif [ $JCAP = 383  ] ; then
    n_proc_gfs=256
    n_proc_rsm=256
    n_proc_gocn=240
    n_proc_rocn=64
    MPI=$((${n_proc_gfs}+${n_proc_rsm}+${n_proc_gocn}+${n_proc_rocn}))
    ln -sf ${commonpath}/data/rmp_tco2timcom_xnew.nc ${GFSWRK}/rmp_tco2timcom.nc
    ln -sf ${commonpath}/data/rmp_timcom2tco_xnew.nc ${GFSWRK}/rmp_timcom2tco.nc
  fi
fi

nodes=$(( (MPI + 47) / 48 ))
if [ $nodes -gt 12 ]; then
   rgrp="large"
else
   rgrp="small"
fi
#######################################################################
if [ "${struc}" = "2cpl_CICE" ]; then
    /usr/bin/bash "${kshpath}/config_gfs_cice_fx1000.sh" "${JCAP}" "${RESTRHR}" "${FCSTGAP}" "${FCSTHR}"
elif [ "${mach}" = "fx1000" ]; then
    if [ "${suffix}" = "000" ]; then
        /usr/bin/bash "${kshpath}/config_gfs_fx1000.sh" "${JCAP}" "${RESTRHR}" "${FCSTGAP}" "${FCSTHR}"
    else
        echo "GFS : in EM case"
        /usr/bin/bash "${kshpath}/config_gfs_fx1000_em.sh" "${JCAP}" "${RESTRHR}" "${FCSTGAP}" "${FCSTHR}"
    fi 
elif [ "${mach}" = "gpu" ]; then
    /usr/bin/bash "${kshpath}/config_gfs.sh" "${JCAP}" "${RESTRHR}" "${FCSTGAP}" "${FCSTHR}"
fi 

if [ ${em_det} = 0  -a ! -f ${commonpath}/data/MASOP_ic/MASOPS_eps${suffix}/*${dtg}*0000/W00100* ] ;then
   INIDMS=${commonpath}/data/MASOP_ic/MASOPS_eps${suffix}_e
else
   INIDMS=${commonpath}/data/MASOP_ic/MASOPS_eps${suffix}
fi

#INIDMS=/data/common/gfs/dms_data/ncep_ana.ufs/TCo${JCAP}l72_${dtg}
#INIDMS=/nwpr/gfs/xb157/data/dmsdb/TCo383L72CFSRn1.ufs/TCo383l72_${dtg}
if [ $JCAP = 199  ] ; then
INIDMS=/data/common/gfs/dms_data/ncep_ana_n2_xnew.ufs/TCo199l72_${dtg}
fi
ln -fs ${INIDMS}/*${dtg}*0000   ${DMSDIR}.ufs/in/
ln -fs ${INIDMS}/*${fgdtg}*0006 ${DMSDIR}.ufs/in/

#BCKDMS=/data/common/gfs/dms_data/bckdms.ufs/BCK_TCo383_GI30S
if [ $JCAP = 199  ] ; then
BCKDMS=/data/common/gfs/dms_data/bckdms.ufs/BCK_TCo${JCAP}_GK30S_1000
ln -fs ${BCKDMS}/* ${DMSDIR}.ufs/bck/
else
BCKDMS=/data/common/gfs/dms_data/bckdms.ufs/BCK_TCo${JCAP}_GI30S_xnew
ln -fs ${BCKDMS}/* ${DMSDIR}.ufs/bck/
fi

if [ "${struc}" = "2cpl_CICE" ]; then
    /usr/bin/bash "${kshpath}/config_glb_cice_fx1000.sh" "${JCAP}" "${RESTRHR}" "${FCSTGAP}" "${FCSTHR}"
elif [ "${mach}" = "fx1000" ]; then
    /usr/bin/bash "${kshpath}/config_glb_fx1000.sh" "${JCAP}" "${RESTRHR}" "${FCSTGAP}" "${FCSTHR}"
elif [ "${mach}" = "gpu" ]; then
    /usr/bin/bash "${kshpath}/config_glb.sh" "${JCAP}" "${RESTRHR}" "${FCSTGAP}" "${FCSTHR}"
fi

ln -fs ${commonpath}/data/timcom_grid_1536x720x55.nc ${GFSWRK}/glb_grid.nc
#-----------------pre 2024/09/04---------
#ln -fs ${commonpath}/data/HYCOM_ic/TIMCOM_g1536_${op_date:0:10}.nc ${GFSWRK}/glb000.nc
#-----------------aft 2024/09/04---------
if [ -f /nwpr/gfs/xb99/GEPSv3/M02new/dtg/glb/${op_date:2:10}/TIMCOM_g1536_${op_date:0:10}.nc ]; then
  ln -fs /nwpr/gfs/xb99/GEPSv3/M02new/dtg/glb/${op_date:2:10}/TIMCOM_g1536_${op_date:0:10}.nc ${GFSWRK}/glb000.nc
else
  ln -fs ${commonpath}/data/HYCOM_ic/TIMCOM_g1536_${op_date:0:10}.nc ${GFSWRK}/glb000.nc
fi

/usr/bin/bash ${kshpath}/config_tai${tai_v}.sh ${FCSTHR}
ln -fs ${commonpath}/data/timcom_grid${tai_v}.nc ${RSMOUT}/tai_grid.nc
#-----------------pre 2024/09/04---------
#ln -fs ${commonpath}/data/HYCOM_ic/tai_rsm5km_${op_date:0:10}${tai_v}.nc ${RSMOUT}/tai000.nc
#-----------------aft 2024/09/04---------
ln -fs ${commonpath}/data/HYCOM_ic_prerun/${dtg}/tai_rsm5km_${op_date:0:10}.nc ${RSMOUT}/tai000.nc
ln -fs ${commonpath}/data/lateral_bc_diff${tai_v}.nc ${RSMOUT}/tai_lbc.nc

glb_path=${kshpath}
TIMCOMr_bcneed=${glb_path}/${dtg}_${suffix}.ufs/TIMCOM_glb
ln -fs ${TIMCOMr_bcneed} ${RSMOUT}/

if [ ${mach} = gpu ] ;then
cat > ${MDIR}.ufs/submit_job.sh <<EOF
#!/bin/bash
#PJM -L vnode=1
#PJM -L vnode-core=32
#PJM -L gpu-share=8
#PJM -L rscunit=rscunit_pg01
#PJM -L rscgrp=gpu-rd-small
#PJM -x PJM_CACHE_MODE=4
#PJM -L elapse=03:00:00
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=${MPI}"
#PJM -j
#PJM -g sum
#PJM -N op_${dtg}
#PJM -o ${MDIR}.ufs/%j.log
#PJM -e ${MDIR}.ufs/%j.err

. /usr/share/Modules/init/bash
module use /package/x86_64/nvidia/hpc_sdk/modulefiles
module load nvhpc-hpcx-cuda12/24.11
module unuse /package/x86_64/nvidia/hpc_sdk/modulefiles

export LD_LIBRARY_PATH=/package/x86_64/nvidia/netcdf-4.9.0/lib:\$LD_LIBRARY_PATH

export n_proc_gfs=${n_proc_gfs}
export n_proc_rsm=${n_proc_rsm}
export n_proc_gocn=${n_proc_gocn}
export n_proc_rocn=${n_proc_rocn}

export MPI=${MPI}
# initial date
export SDATE=${op_date}
# restart hour (zero means NOT restart)
export RESTRHR=${RESTRHR}
# forecast hour
export FCSTHR=${FCSTHR}

# running directory
cd ${GFSWRK}
export GFSWRK=${GFSWRK}
export ANADMS=in@${DMSDIR}
export FCSTDMS=TCo${JCAP}@${MDIR}
export BCKOPS=bck@${DMSDIR}
export NWPETCGLB=${GFSWRK}
export GLB_TYPHINI=/ncs/ncsatyp/TYP/M00/dtg/ty
export FIXDIR=${GFSDIR}/fix
# << Enviroment variables of OpenACC with NVCOMPILER >>
export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC="1"
export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC_POOLSIZE="40G"
export NVCOMPILER_ACC_USE_GRAPH="1"
export NVCOMPILER_ACC_CUDA_NOCOPY="1"


FCT_MODEL=${FCT_MODEL}
/usr/bin/time -p mpiexec -n ${MPI} ${FCT_MODEL} 000
if [ \$? != 0 ] ; then
  echo "error occurred: fct model failed!!"
fi

echo "Ending at: " \$(date)
EOF
fi

if [ ${mach} = fx1000 -a ${struc} = 4cpl ] ;then
cat > ${MDIR}.ufs/submit_job.sh << EOF
#!/bin/bash
#PJM -L "node=${nodes}:noncont"
#PJM -L rscgrp=${rgrp}
#PJM -x PJM_CACHE_MODE=4
#PJM -L elapse=18:00:00
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=${MPI}"
#PJM -j
#PJM -g sum
#PJM -N op_${dtg}_${suffix}
#PJM -o ${MDIR}.ufs/%j.log
#PJM -e ${MDIR}.ufs/%j.err
#PJM -S
#PJM --spath  ${MDIR}.ufs/%n.i%j

OMP=1
source /users/xa09/sample/setup_mpi+omp.fx1000 \$OMP

export n_proc_gfs=${n_proc_gfs}
export n_proc_rsm=${n_proc_rsm}
export n_proc_gocn=${n_proc_gocn}
export n_proc_rocn=${n_proc_rocn}
export MPI=${MPI}

export RSMDISK=${RSMDIR}
# initial date
export SDATE=${op_date}
# restart hour (zero means NOT restart)
export RESTRHR=${RESTRHR}
# forecast hour
export FCSTHR=${FCSTHR}
# running directory
export RUNDIR=${RSMOUT}
export RSMPATH=${RSMDIR}/run
cd \${RSMPATH}
. \${RSMPATH}/run_4cpl.sh \${RSMDISK} \${SDATE} \${RESTRHR} \${FCSTHR} \${RUNDIR}

cd ${GFSWRK}
export GFSWRK=${GFSWRK}
export ANADMS=in@${DMSDIR}
export FCSTDMS=TCo${JCAP}@${MDIR}
export BCKOPS=bck@${DMSDIR}
export NWPETCGLB=${GFSWRK}
export GLB_TYPHINI=/ncs/ncsatyp/TYP/M00/dtg/ty
export FIXDIR=${GFSDIR}/fix

FCT_MODEL=${FCT_MODEL}
/usr/bin/time -p mpiexec -stdin ../rsm/rfcstparm.all --mca coll_base_reduce_commute_safe 1 -n \${MPI} \${FCT_MODEL} 000 -Wl,-T11,-T12,-T18

if [ \$? != 0 ] ; then
  echo "error occured: fct model fail !!"
fi

echo "Ending at: " \$(date)
EOF
fi
if [ ${mach} = fx1000 -a \( "${struc}" = "2cpl" -o "${struc}" = "2cpl_CICE" \) ] ;then
cat > ${MDIR}.ufs/submit_job.sh << EOF
#!/bin/bash
#PJM -L "node=${nodes}:noncont"
#PJM -L rscgrp=${rgrp}
#PJM -x PJM_CACHE_MODE=4
#PJM -L elapse=12:00:00
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=${MPI}"
#PJM -j
#PJM -g sum
#PJM -N op_${dtg}_${suffix}
#PJM -o ${MDIR}.ufs/%j.log
#PJM -e ${MDIR}.ufs/%j.err
#PJM -S
#PJM --spath  ${MDIR}.ufs/%n.i%j

OMP=1
source /users/xa09/sample/setup_mpi+omp.fx1000 \$OMP

export n_proc_gfs=${n_proc_gfs}
export n_proc_rsm=${n_proc_rsm}
export n_proc_gocn=${n_proc_gocn}
export n_proc_rocn=${n_proc_rocn}
export MPI=${MPI}

export RSMDISK=${RSMDIR}
# initial date
export SDATE=${op_date}
# restart hour (zero means NOT restart)
export RESTRHR=${RESTRHR}
# forecast hour
export FCSTHR=${FCSTHR}
# running directory
export RUNDIR=${RSMOUT}
export RSMPATH=${RSMDIR}/run

cd ${GFSWRK}
export GFSWRK=${GFSWRK}
export ANADMS=in@${DMSDIR}
export FCSTDMS=TCo${JCAP}@${MDIR}
export BCKOPS=bck@${DMSDIR}
export NWPETCGLB=${GFSWRK}
export GLB_TYPHINI=/ncs/ncsatyp/TYP/M00/dtg/ty
export FIXDIR=${GFSDIR}/fix

FCT_MODEL=${FCT_MODEL}
/usr/bin/time -p mpiexec --mca coll_base_reduce_commute_safe 1 -n \${MPI} \${FCT_MODEL} 000 -Wl,-T11,-T12,-T18

if [ \$? != 0 ] ; then
  echo "error occured: fct model fail !!"
fi

echo "Ending at: " \$(date)
EOF

fi
