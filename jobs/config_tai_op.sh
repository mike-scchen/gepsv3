#!/bin/bash

cat > ${RSMOUT}/namelist000.run << EOF
&general_info
 path_output   = '${TAIOUT}/'
 ocn_grid_num  = 1
 ocn_grid_list = 'tai:NoData'
 atm_grid_num  = 0
 atm_grid_list = 'NoData'
 ice_grid_num  = 0
 ice_grid_list = 'NoData'
 rof_grid_num  = 1
 rof_grid_list = 'rof:NoData'
 daodt_drv = 12
 dtg_init = ${dtg}
 duration = ${FCSTHR}
 duration_unit = 'hour'
/

&restart_info
 restart = ${ocn_restart}
 path_restart = '${RSMOUT}/'
/
EOF

cat > ${RSMOUT}/namelist000.inp_tai << EOF
&grid_info
 nx_grid = 288
 ny_grid = 288
 nz_grid = 55
 daodt = 1920
/

&mpi_info
 npx = 4
 npy = 16
 peri_x = .false.
 peri_y = .false.
 symm_np = .false.
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
 threshold_bicg = 1.d-16
 threshold_p0   = 1.d-6
 max_iter_p0    = 100
/
EOF

cat > ${RSMOUT}/namelist000.rof << EOF
&drof_roff
 fpath = '/data/common/gfs/GEPSv3_lib/data/'
 fname = 'TIMCOM_drof_roff_tai_op.nc'
 gridx = 'lon'
 gridy = 'lat'
 time  = 'time'
 vara  = 'roff'
&end
EOF
