#!/bin/bash

export kshpath=$(pwd)
export caldtg="/data/common/gfs/scripts/Caldtg.ksh"
export op_date=2025072800  #${1} #'2001070500'
export suffix=''
export RESTRHR=0 #${2} #0
export FCSTGAP=6 #${3} #6
export FCSTHR=48 #${4}  #1080
export fyy=${op_date:0:4}
export fmm=${op_date:4:2}
export fdd=${op_date:6:2}

export dtg=${op_date:2:8}
export fgdtg=$(${caldtg} ${dtg} -6)

export ocn_restart=.false.
export MDIR=${kshpath}/${dtg}${suffix}
export GLBOUT=${MDIR}/TIMCOM_glb

mkdir -p ${MDIR}
mkdir -p ${GLBOUT}

ln -fs /data/common/gfs/GEPSv3_lib/data/timcom_grid_1536x720x55.nc ${MDIR}/glb_grid.nc
ln -fs /data/common/gfs/GEPSv3_lib/data/HYCOM_ic/TIMCOM_g1536_${op_date:0:8}00.nc ${MDIR}/glb000.nc

cat > ${MDIR}/namelist000.run << EOF
&general_info
 path_output   = '${GLBOUT}/'
 ocn_grid_num  = 1
 ocn_grid_list = 'glb:NoData'
 atm_grid_num  = 1
 atm_grid_list = 'jra:NoData'
 ice_grid_num  = 1
 ice_grid_list = 'ice:NoData'
 rof_grid_num  = 1
 rof_grid_list = 'rof:NoData'
 daodt_drv = 12
 dtg_init = ${dtg}
 duration = ${FCSTHR}
 duration_unit = 'hour'
/

&restart_info
 restart = ${ocn_restart}
 path_restart = '${MDIR}/'
/
EOF

cat > ${MDIR}/namelist000.inp_glb << EOF
&grid_info
 nx_grid = 1536
 ny_grid = 720
 nz_grid = 55
 daodt = 1920
/

&mpi_info
 npx = ${NPX}
 npy = ${NPY}
 peri_x = .true.
 peri_y = .false.
 symm_np = .true.
/

&io_info
 freq_out = ${FCSTGAP}
 freq_restart = 24
 freq_out_unit = 'hour'
/

&phy_info
 opt_windmix = 1
 opt_arbr_p0 = 0
 dm0 = 5.d5
 opt_da = 0
/

&solver_info
 opt_solver     = 2
 threshold_bicg = 1.d-16
 threshold_p0   = 1.d-6
 max_iter_p0    = 600
 PCSICriterion  = 1.d-24
 pcsi_cuda_graph = .true.
/

EOF


cat > ${MDIR}/namelist000.jra << EOF
&datm_uwnd
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_u10m_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'u10m'
&end

&datm_vwnd
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_v10m_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'v10m'
&end

&datm_temp
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_t10m_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 't10m'
&end

&datm_qhum
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_q10m_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'q10m'
&end

&datm_pslv
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_SLP_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'pslv'
&end

&datm_swup
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_swup_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'swup'
&end

&datm_swdn
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_swdn_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'swdn'
&end

&datm_lwdn
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_lwdn_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'lwdn'
&end

&datm_rain
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_rain_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'rain'
&end

&datm_snow
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_snow_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'snow'
&end
EOF

cat > ${MDIR}/namelist000.ice << EOF
&dice_ifrac
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TIMCOM_dice_frac.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'frac'
&end
EOF

cat > ${MDIR}/namelist000.rof << EOF
&drof_roff
 fpath = '/data/common/gfs/GEPSv3_lib/data/JRA55_forcing/fix/'
 fname = 'GEPSv3_JRA_forcing_roff_${fyy}_${fmm}.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'roff'
&end
EOF
