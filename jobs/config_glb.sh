#!/bin/bash

cat > ${GFSWRK}/namelist000.run << EOF
&general_info
 path_output   = '${GLBOUT}/'
 ocn_grid_num  = 1
 ocn_grid_list = 'glb:NoData'
 atm_grid_num  = 0
 atm_grid_list = 'tco:NoData'
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
 path_restart = '${GFSWRK}/'
/
EOF

cat > ${GFSWRK}/namelist000.inp_glb << EOF
&grid_info
 nx_grid = 1536
 ny_grid = 720
 nz_grid = 55
 daodt = 1920
/

&mpi_info
 npx = 4
 npy = 2
 peri_x = .true.
 peri_y = .false.
 symm_np = .true.
/

&io_info
 freq_out = ${FCSTGAP}
 freq_restart = ${FCSTHR}
 freq_out_unit = 'hour'
/

&phy_info
 opt_windmix = 1
 opt_arbr_p0 = 0
 dm0 = 5.d5
 opt_da = 0
/

&solver_info
 opt_solver = 2
 threshold_bicg = 1.d-16
 threshold_p0   = 1.d-6
 max_iter_p0    = 600
 PCSICriterion  = 1.d-24
 pcsi_cuda_graph = .false.
/
EOF


cat > ${GFSWRK}/namelist000.ice << EOF
&dice_ifrac
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TIMCOM_dice_frac.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'frac'
&end
EOF

cat > ${GFSWRK}/namelist000.rof << EOF
&drof_roff
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TIMCOM_drof_roff.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'roff'
&end
EOF



cat > ${GFSWRK}/namelist000.tco << EOF
&datm_uwnd
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'u10m'
&end

&datm_vwnd
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'v10m'
&end

&datm_temp
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 't02m'
&end

&datm_qhum
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'q02m'
&end

&datm_pslv
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'pslv'
&end

&datm_swup
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'swup'
&end

&datm_swdn
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'swdn'
&end

&datm_lwdn
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'lwdn'
&end

&datm_rain
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'rain'
&end

&datm_snow
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TCo383_2001010100_atmo.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'snow'
&end
EOF

