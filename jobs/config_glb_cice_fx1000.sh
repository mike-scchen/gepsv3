#!/bin/bash
kshpath=`pwd `
. ${kshpath}/00_setDate.sh

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

if [ $1 = 383  ] ; then
x_grid=1536
y_grid=720
npex=24
npey=16
cice_x=$((x_grid / npex))
cice_y=$((y_grid / npey))
c_np=$((npex * npey))
cat > ${GFSWRK}/namelist000.inp_glb << EOF
&grid_info
 nx_grid = ${x_grid}
 ny_grid = ${y_grid}
 nz_grid = 55
 daodt = 1920
/

&mpi_info
 npx = ${npex}
 npy = ${npey}
 peri_x = .true.
 peri_y = .false.
 symm_np = .true.
/

&io_info
 freq_out = ${FCSTGAP}
 freq_restart = ${FRE_R}
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

&init_info
 opt_ic = 1
 ts3z_file = './hycom_ts3z.nc'
 uv3z_file = './hycom_uv3z.nc'
 ssh_file  = './hycom_ssh.nc'
/
EOF
elif [ $1 = 199 ] ; then
x_grid=1536
y_grid=720
npex=24
npey=16
cice_x=$((x_grid / npex))
cice_y=$((y_grid / npey))
c_np=$((npex * npey))
cat > ${GFSWRK}/namelist000.inp_glb << EOF
&grid_info
 nx_grid = ${x_grid}
 ny_grid = ${y_grid}
 nz_grid = 55
 daodt = 1920
/

&mpi_info
 npx = ${npex}
 npy = ${npey}
 peri_x = .true.
 peri_y = .false.
 symm_np = .true.
/

&io_info
 freq_out = ${FCSTGAP}
 freq_restart = ${FRE_R}
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
fi

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

cat > ${GFSWRK}/ice_in << EOF
&setup_nml
    days_per_year  = 365
    use_leap_years = .true.
    year_init      = ${ice_year}
    month_init     = ${ice_mon}
    day_init       = ${ice_day}
    sec_init       = 0
    istep0         = 0
    dt             = 600.0
    dt_cpl         = 7200.0
    npt_unit       = 'd'
    npt            = ${ft_day}
    ndtd           = 2
    runtype        = 'initial'
    ice_ic         = 'internal'
    restart_ext    = .false.
    use_restart_time = .false.
    restart_format = 'default'
    lcdf64         = .false.
    numin          = 21
    numax          = 89
    restart_dir    = '../CICE/restart/'
    restart_file   = 'iced'
    pointer_file   = '../CICE/ice.restart_file'
    dumpfreq       = 'd','x','x','x','x'
    dumpfreq_n     =  ${ft_day}, 1 , 1 , 1 , 1
    dumpfreq_base  = 'init','init','init','init','init'
    dump_last      = .false.
    bfbflag        = 'off'
    diagfreq       = 1
    diag_type      = 'stdout'
    diag_file      = 'ice_diag.d'
    debug_model    = .false.
    debug_model_step = 0
    debug_model_i  = -1
    debug_model_j  = -1
    debug_model_iblk = -1
    debug_model_task = -1
    debug_forcing  = .false.
    print_global   = .true.
    print_points   = .true.
    timer_stats    = .false.
    memory_stats   = .false.
    conserv_check  = .false.
    latpnt(1)      =  90.
    lonpnt(1)      =  0.
    latpnt(2)      = -65.
    lonpnt(2)      = -45.
    histfreq       = 'd','x','x','x','x'
    histfreq_n     =  1 , 1 , 1 , 1 , 1
    histfreq_base  = 'init','init','init','init','init'
    hist_avg       = .true.,.true.,.true.,.true.,.true.
    history_dir    = '../CICE/'
    history_file   = 'iceh'
    history_precision = 4
    history_format = 'default'
    hist_time_axis = 'end'
    write_ic       = .true.
    incond_dir     = '../CICE/'
    incond_file    = 'iceh_ic'
    version_name   = 'CICE_6.5.0'
/

&grid_nml
    grid_format  = 'nc'
    grid_type    = 'latlon'
    grid_ice     = 'B'
    grid_atm     = 'A'
    grid_ocn     = 'A'
    grid_file    = '/data/common/gfs/GEPSv3_lib/data/CICE_data/grid/0p25/cice_grid_0p25.nc'
    kmt_type     = 'file'
    kmt_file     = '/data/common/gfs/GEPSv3_lib/data/CICE_data/grid/0p25/cice_kmt_0p25.nc'
    bathymetry_file = '/data/common/gfs/GEPSv3_lib/data/CICE_data/grid/0p25/timcom_lev_thick.txt'
    bathymetry_format = 'pop'
    use_bathymetry = .True.
    gridcpl_file = 'unknown_gridcpl_file'
    kcatbound    = 0
    dxrect         = 30.e5
    dyrect         = 30.e5
    lonrefrect    = -156.50d0
    latrefrect    =   71.35d0
    scale_dxdy     = .false.
    dxscale        = 1.d0
    dyscale        = 1.d0
    close_boundaries = .false.
    ncat         = 5
    nfsd         = 1
    nilyr        = 4
    nslyr        = 1
    nblyr        = 1
    orca_halogrid = .false.
/

&tracer_nml
    n_iso        = 0
    n_aero       = 0
    n_zaero      = 0
    n_algae      = 0
    n_doc        = 0
    n_dic        = 0
    n_don        = 0
    n_fed        = 0
    n_fep        = 0
    tr_iage      = .false.
    restart_age  = .false.
    tr_FY        = .false.
    restart_FY   = .false.
    tr_lvl       = .false.
    restart_lvl  = .false.
    tr_pond_topo = .false.
    restart_pond_topo = .false.
    tr_pond_lvl  = .false.
    restart_pond_lvl  = .false.
    tr_snow      = .false.
    restart_snow = .false.
    tr_iso       = .false.
    restart_iso  = .false.
    tr_aero      = .false.
    restart_aero = .false.
    tr_fsd       = .false.
    restart_fsd  = .false.
/

&thermo_nml
    kitd              = 1
    ktherm            = 2
    conduct           = 'MU71'
    ksno              = 0.3d0
    hi_min            = 0.01d0
    a_rapid_mode      =  0.5e-3
    Rac_rapid_mode    =    10.0
    aspect_rapid_mode =     1.0
    dSdt_slow_mode    = -1.5e-7
    phi_c_slow_mode   =    0.05
    phi_i_mushy       =    0.85
    Tliquidus_max     = -0.1d0
    hfrazilmin        = 0.05d0
    floediam          = 300.0d0
/

&dynamics_nml
    kdyn            = 1
    ndte            = 120
    revised_evp     = .false.
    evp_algorithm   = 'standard_2d'
    brlx            = 300.0
    arlx            = 300.0
    advection       = 'remap'
    kstrength       = 1
    krdg_partic     = 1
    krdg_redist     = 1
    mu_rdg          = 4
    Pstar           = 2.75e4
    Cstar           = 20
    Cf              = 17.
    Ktens           = 0.
    e_yieldcurve    = 2.
    e_plasticpot    = 2.
    visc_method     = 'avg_zeta'
    elasticDamp     = 0.36d0
    deltaminEVP     = 1e-11
    deltaminVP      = 2e-9
    capping_method  = 'max'
    seabed_stress   = .false.
    seabed_stress_method = 'LKD'
    k1              = 8.
    k2              = 15.
    alphab          = 20.
    threshold_hw    = 30.
    coriolis        = 'latitude'
    kridge          = 1
    ktransport      = 1
    ssh_stress      = 'coupled'
    maxits_nonlin   = 10
    precond         = 'pgmres'
    dim_fgmres       = 50
    dim_pgmres       = 5
    maxits_fgmres   = 1
    maxits_pgmres   = 1
    monitor_nonlin  = .false.
    monitor_fgmres  = .false.
    monitor_pgmres  = .false.
    ortho_type      = 'mgs'
    reltol_nonlin   = 1e-8
    reltol_fgmres   = 1e-1
    reltol_pgmres   = 1e-6
    algo_nonlin     = 'picard'
    use_mean_vrel   = .true.
/

&shortwave_nml
    shortwave       = 'dEdd'
    snw_ssp_table   = 'test'
    albedo_type     = 'ccsm3'
    albicev         = 0.75
    albicei         = 0.45
    albsnowv        = 0.98
    albsnowi        = 0.73 
    ahmax           = 0.3
    R_ice           = 0.
    R_pnd           = 0.
    R_snw           = 1.5
    dT_mlt          = 1.5
    rsnw_mlt        = 1500.
    kalg            = 0.6
    sw_redist         = .false.
    sw_frac           = 0.9d0
    sw_dtemp          = 0.02d0
/

&ponds_nml
    hp1             = 0.01
    hs0             = 0.
    hs1             = 0.03
    dpscale         = 1.e-3
    frzpnd          = 'hlid'
    rfracmin        = 0.15
    rfracmax        = 1.
    pndaspect       = 0.8
/

&snow_nml
    snwredist       = 'none'
    snwgrain        = .false.
    use_smliq_pnd   = .false.
    rsnw_fall       =  100.0
    rsnw_tmax       = 1500.0
    rhosnew         =  100.0
    rhosmin         =  100.0
    rhosmax         =  450.0
    windmin         =   10.0
    drhosdwind      =   27.3
    snwlvlfac       =    0.3
    snw_aging_table = 'test'
    snw_filename    = 'unknown'
    snw_rhos_fname  = 'unknown'
    snw_Tgrd_fname  = 'unknown'
    snw_T_fname     = 'unknown'
    snw_tau_fname   = 'unknown'
    snw_kappa_fname = 'unknown'
    snw_drdt0_fname = 'unknown'
/

&forcing_nml
    formdrag        = .false.
    atmbndy         = 'similarity'
    rotate_wind     = .false.
    calc_strair     = .true.
    calc_Tsfc       = .true.
    highfreq        = .false.
    natmiter        = 5
    atmiter_conv    = 0.0d0
    ustar_min       = 0.0005
    iceruf          = 0.0005
    calc_dragio     = .false.
    iceruf_ocn      = 0.03
    emissivity      = 0.985
    fbot_xfer_type  = 'constant'
    update_ocn_f    = .false.
    l_mpond_fresh   = .false.
    tfrz_option     = 'mushy'
    saltflux_option = 'constant'
    ice_ref_salinity = 4.0
    oceanmixed_ice  = .false.
    wave_spec_type  = 'none'
    wave_spec_file  = 'unknown_wave_spec_file'
    nfreq           = 25
    restart_coszen  = .true.
    restore_ice     = .false.
    restore_ocn     = .false.
    trestore        =  90
    precip_units    = 'mks'
    default_season  = 'winter'
    atm_data_type   = 'default'
    atm_data_version = ''
    ocn_data_type   = 'default'
    bgc_data_type   = 'default'
    fe_data_type    = 'default'
    ice_data_type   = 'latsst'
    ice_data_conc   = 'parabolic'
    ice_data_dist   = 'uniform'
    fyear_init      = 2005
    ycycle          = 1
    atm_data_format = 'nc'
    atm_data_dir    = '/data/common/gfs/GEPSv3_lib/data/CICE_data/forcing/0p25'
    bgc_data_dir    = '/data/common/gfs/GEPSv3_lib/data/CICE_data/forcing/0p25/WOA/MONTHLY'
    ocn_data_format = 'bin'
    ocn_data_dir    = '/data/common/gfs/GEPSv3_lib/data/CICE_data/forcing/0p25/CESM/MONTHLY'
    oceanmixed_file = 'unknown_oceanmixed_file'
/

&domain_nml
    nprocs = ${c_np}
    nx_global         = 1536
    ny_global         = 720
    block_size_x      = ${cice_x}
    block_size_y      = ${cice_y}
    max_blocks        = 1
    processor_shape   = 'square-ice'
    distribution_type = 'cartesian'
    distribution_wght = 'blockall'
    distribution_wght_file = 'unknown'
    ew_boundary_type  = 'cyclic'
    ns_boundary_type  = 'open'
    maskhalo_dyn      = .true.
    maskhalo_remap    = .true.
    maskhalo_bound    = .true.
    add_mpi_barriers  = .false.
    debug_blocks      = .false.
/

&zbgc_nml
    tr_brine        = .false.
    restart_hbrine  = .false.
    tr_zaero        = .false.
    modal_aero      = .false.
    skl_bgc         = .false.
    z_tracers       = .false.
    dEdd_algae      = .false.
    solve_zbgc      = .false.
    bgc_flux_type   = 'Jin2006'
    restore_bgc     = .false.
    restart_bgc     = .false.
    scale_bgc       = .false.
    solve_zsal      = .false.
    restart_zsal    = .false.
    tr_bgc_Nit      = .false.
    tr_bgc_C        = .false.
    tr_bgc_chl      = .false.
    tr_bgc_Am       = .false.
    tr_bgc_Sil      = .false.
    tr_bgc_DMS      = .false.
    tr_bgc_PON      = .false.
    tr_bgc_hum      = .false.
    tr_bgc_DON      = .false.
    tr_bgc_Fe       = .false. 
    grid_o          = 0.006
    grid_o_t        = 0.006
    l_sk            = 0.024
    grid_oS         = 0.0
    l_skS           = 0.028
    phi_snow        = -0.3
    initbio_frac    = 0.8
    frazil_scav     = 0.8  
    ratio_Si2N_diatoms = 1.8                         
    ratio_Si2N_sp      = 0.0
    ratio_Si2N_phaeo   = 0.0
    ratio_S2N_diatoms  = 0.03  
    ratio_S2N_sp       = 0.03 
    ratio_S2N_phaeo    = 0.03
    ratio_Fe2C_diatoms = 0.0033
    ratio_Fe2C_sp      = 0.0033
    ratio_Fe2C_phaeo   = 0.1
    ratio_Fe2N_diatoms = 0.023 
    ratio_Fe2N_sp      = 0.023
    ratio_Fe2N_phaeo   = 0.7
    ratio_Fe2DON       = 0.023
    ratio_Fe2DOC_s     = 0.1
    ratio_Fe2DOC_l     = 0.033
    fr_resp            = 0.05
    tau_min            = 5200.0
    tau_max            = 173000.0
    algal_vel          = 0.0000000111
    R_dFe2dust         = 0.035
    dustFe_sol         = 0.005
    chlabs_diatoms     = 0.03
    chlabs_sp          = 0.01
    chlabs_phaeo       = 0.05
    alpha2max_low_diatoms = 0.8
    alpha2max_low_sp      = 0.67
    alpha2max_low_phaeo   = 0.67
    beta2max_diatoms   = 0.018
    beta2max_sp        = 0.0025
    beta2max_phaeo     = 0.01
    mu_max_diatoms     = 1.44
    mu_max_sp          = 0.851
    mu_max_phaeo       = 0.851
    grow_Tdep_diatoms  = 0.06
    grow_Tdep_sp       = 0.06
    grow_Tdep_phaeo    = 0.06
    fr_graze_diatoms   = 0.0
    fr_graze_sp        = 0.1
    fr_graze_phaeo     = 0.1
    mort_pre_diatoms   = 0.007
    mort_pre_sp        = 0.007
    mort_pre_phaeo     = 0.007
    mort_Tdep_diatoms  = 0.03
    mort_Tdep_sp       = 0.03
    mort_Tdep_phaeo    = 0.03
    k_exude_diatoms    = 0.0
    k_exude_sp         = 0.0
    k_exude_phaeo      = 0.0
    K_Nit_diatoms      = 1.0
    K_Nit_sp           = 1.0
    K_Nit_phaeo        = 1.0
    K_Am_diatoms       = 0.3
    K_Am_sp            = 0.3
    K_Am_phaeo         = 0.3
    K_Sil_diatoms      = 4.0
    K_Sil_sp           = 0.0
    K_Sil_phaeo        = 0.0
    K_Fe_diatoms       = 1.0
    K_Fe_sp            = 0.2
    K_Fe_phaeo         = 0.1
    f_don_protein      = 0.6
    kn_bac_protein     = 0.03
    f_don_Am_protein   = 0.25
    f_doc_s            = 0.4
    f_doc_l            = 0.4
    f_exude_s          = 1.0
    f_exude_l          = 1.0
    k_bac_s            = 0.03
    k_bac_l            = 0.03
    T_max              = 0.0
    fsal               = 1.0
    op_dep_min         = 0.1
    fr_graze_s         = 0.5
    fr_graze_e         = 0.5
    fr_mort2min        = 0.5
    fr_dFe             = 0.3
    k_nitrif           = 0.0
    t_iron_conv        = 3065.0
    max_loss           = 0.9
    max_dfe_doc1       = 0.2
    fr_resp_s          = 0.75
    y_sk_DMS           = 0.5
    t_sk_conv          = 3.0
    t_sk_ox            = 10.0
    algaltype_diatoms  = 0.0
    algaltype_sp       = 0.5
    algaltype_phaeo    = 0.5
    nitratetype        = -1.0
    ammoniumtype       = 1.0
    silicatetype       = -1.0
    dmspptype          = 0.5
    dmspdtype          = -1.0
    humtype            = 1.0
    doctype_s          = 0.5
    doctype_l          = 0.5
    dontype_protein    = 0.5
    fedtype_1          = 0.5
    feptype_1          = 0.5
    zaerotype_bc1      = 1.0
    zaerotype_bc2      = 1.0
    zaerotype_dust1    = 1.0
    zaerotype_dust2    = 1.0
    zaerotype_dust3    = 1.0
    zaerotype_dust4    = 1.0
    ratio_C2N_diatoms  = 7.0
    ratio_C2N_sp       = 7.0
    ratio_C2N_phaeo    = 7.0
    ratio_chl2N_diatoms= 2.1
    ratio_chl2N_sp     = 1.1
    ratio_chl2N_phaeo  = 0.84
    F_abs_chl_diatoms  = 2.0
    F_abs_chl_sp       = 4.0
    F_abs_chl_phaeo    = 5.0
    ratio_C2N_proteins = 7.0
/

&icefields_nml
    f_tmask        = .true.
    f_umask        = .false.
    f_nmask        = .false.
    f_emask        = .false.
    f_blkmask      = .true.
    f_tarea        = .true.
    f_uarea        = .true.
    f_narea        = .false.
    f_earea        = .false.
    f_dxt          = .false.
    f_dyt          = .false.
    f_dxu          = .false.
    f_dyu          = .false.
    f_dxe          = .false.
    f_dye          = .false.
    f_dxn          = .false.
    f_dyn          = .false.
    f_HTN          = .false.
    f_HTE          = .false.
    f_ANGLE        = .true.
    f_ANGLET       = .true.
    f_NCAT         = .true.
    f_VGRDi        = .false.
    f_VGRDs        = .false.
    f_VGRDb        = .false.
    f_VGRDa        = .true.
    f_bounds       = .false.
    f_aice         = 'd1' 
    f_hi           = 'd1'
    f_hs           = 'd1' 
    f_Tsfc         = 'd1' 
    f_sice         = 'd1' 
    f_uvel         = 'd1' 
    f_vvel         = 'd1' 
    f_uatm         = 'd1' 
    f_vatm         = 'd1' 
    f_fswdn        = 'd1' 
    f_flwdn        = 'd1'
    f_snowfrac     = 'd1'
    f_snow         = 'd1' 
    f_snow_ai      = 'x' 
    f_rain         = 'd1' 
    f_rain_ai      = 'x' 
    f_sst          = 'd1' 
    f_sss          = 'd1' 
    f_uocn         = 'd1' 
    f_vocn         = 'd1' 
    f_frzmlt       = 'd1'
    f_fswfac       = 'd1'
    f_fswint_ai    = 'x'
    f_fswabs       = 'x' 
    f_fswabs_ai    = 'x' 
    f_albsni       = 'x' 
    f_alvdr        = 'x'
    f_alidr        = 'x'
    f_alvdf        = 'x'
    f_alidf        = 'x'
    f_alvdr_ai     = 'x'
    f_alidr_ai     = 'x'
    f_alvdf_ai     = 'x'
    f_alidf_ai     = 'x'
    f_albice       = 'x'
    f_albsno       = 'x'
    f_albpnd       = 'x'
    f_coszen       = 'x'
    f_flat         = 'd' 
    f_flat_ai      = 'x' 
    f_fsens        = 'd' 
    f_fsens_ai     = 'x' 
    f_fswup        = 'd' 
    f_flwup        = 'd' 
    f_flwup_ai     = 'x' 
    f_evap         = 'd' 
    f_evap_ai      = 'x' 
    f_Tair         = 'd' 
    f_Tref         = 'd' 
    f_Qref         = 'd'
    f_congel       = 'd' 
    f_frazil       = 'd' 
    f_snoice       = 'd' 
    f_dsnow        = 'x' 
    f_melts        = 'd'
    f_meltt        = 'd'
    f_meltb        = 'd'
    f_meltl        = 'd'
    f_fresh        = 'd'
    f_fresh_ai     = 'x'
    f_fsalt        = 'd'
    f_fsalt_ai     = 'x'
    f_fbot         = 'd'
    f_fhocn        = 'd' 
    f_fhocn_ai     = 'x' 
    f_fswthru      = 'd' 
    f_fswthru_ai   = 'x' 
    f_fsurf_ai     = 'x'
    f_fcondtop_ai  = 'x'
    f_fmeltt_ai    = 'd' 
    f_strairx      = 'd' 
    f_strairy      = 'd' 
    f_strtltx      = 'd' 
    f_strtlty      = 'd' 
    f_strcorx      = 'x' 
    f_strcory      = 'x' 
    f_strocnx      = 'd' 
    f_strocny      = 'd' 
    f_strintx      = 'x' 
    f_strinty      = 'x'
    f_taubx        = 'x'
    f_tauby        = 'x'
    f_strength     = 'd'
    f_divu         = 'x'
    f_shear        = 'x'
    f_sig1         = 'x' 
    f_sig2         = 'x' 
    f_sigP         = 'x' 
    f_dvidtt       = 'x' 
    f_dvidtd       = 'x' 
    f_daidtt       = 'x'
    f_daidtd       = 'x' 
    f_dagedtt      = 'x'
    f_dagedtd      = 'x' 
    f_mlt_onset    = 'd'
    f_frz_onset    = 'd'
    f_hisnap       = 'x'
    f_aisnap       = 'x'
    f_trsig        = 'd'
    f_icepresent   = 'd'
    f_iage         = 'd'
    f_FY           = 'x'
    f_aicen        = 'x'
    f_vicen        = 'x'
    f_vsnon        = 'x'
    f_snowfracn    = 'x'
    f_keffn_top    = 'x'
    f_Tinz         = 'x'
    f_Sinz         = 'x'
    f_Tsnz         = 'x'
    f_fsurfn_ai    = 'x'
    f_fcondtopn_ai = 'x'
    f_fmelttn_ai   = 'x'
    f_flatn_ai     = 'x'
    f_fsensn_ai    = 'x'
    f_CMIP         = 'x'
/

&icefields_mechred_nml
    f_alvl         = 'm'
    f_vlvl         = 'm'
    f_ardg         = 'm'
    f_vrdg         = 'm'
    f_dardg1dt     = 'x'
    f_dardg2dt     = 'x'
    f_dvirdgdt     = 'x'
    f_opening      = 'x'
    f_ardgn        = 'x'
    f_vrdgn        = 'x'
    f_dardg1ndt    = 'x'
    f_dardg2ndt    = 'x'
    f_dvirdgndt    = 'x'
    f_krdgn        = 'x'
    f_aparticn     = 'x'
    f_aredistn     = 'x'
    f_vredistn     = 'x'
    f_araftn       = 'x'
    f_vraftn       = 'x'
/

&icefields_pond_nml
    f_apondn       = 'x'
    f_apeffn       = 'x'
    f_hpondn       = 'x'
    f_apond        = 'm'
    f_hpond        = 'm'
    f_ipond        = 'm'
    f_apeff        = 'm'
    f_apond_ai     = 'm'
    f_hpond_ai     = 'm'
    f_ipond_ai     = 'm'
    f_apeff_ai     = 'm'
/

&icefields_snow_nml
    f_smassicen    = 'x'
    f_smassliqn    = 'x'
    f_rhos_cmpn    = 'x'
    f_rhos_cntn    = 'x'
    f_rsnwn        = 'x'
    f_smassice     = 'm'
    f_smassliq     = 'm'
    f_rhos_cmp     = 'm'
    f_rhos_cnt     = 'm'
    f_rsnw         = 'm'
    f_meltsliq     = 'm'
    f_fsloss       = 'm'
/

&icefields_bgc_nml
    f_fiso_atm     = 'x'
    f_fiso_ocn     = 'x'
    f_iso          = 'x'
    f_faero_atm    = 'x'
    f_faero_ocn    = 'x'
    f_aero         = 'x'
    f_fbio         = 'm'
    f_fbio_ai      = 'm'
    f_zaero        = 'x'
    f_bgc_S        = 'm'
    f_bgc_N        = 'm'
    f_bgc_C        = 'x'
    f_bgc_DOC      = 'm'
    f_bgc_DIC      = 'x'
    f_bgc_chl      = 'x'
    f_bgc_Nit      = 'm'
    f_bgc_Am       = 'm'
    f_bgc_Sil      = 'm'
    f_bgc_DMSPp    = 'x'
    f_bgc_DMSPd    = 'x'
    f_bgc_DMS      = 'x'
    f_bgc_DON      = 'x'  
    f_bgc_Fe       = 'm'  
    f_bgc_hum      = 'm'   
    f_bgc_PON      = 'm'
    f_bgc_ml       = 'm'
    f_upNO         = 'm'
    f_upNH         = 'm'
    f_bTin         = 'm'
    f_bphi         = 'm' 
    f_iDi          = 'm'
    f_iki          = 'm'
    f_fbri         = 'm'  
    f_hbri         = 'm'
    f_zfswin       = 'm'
    f_bionet       = 'm'
    f_biosnow      = 'm'
    f_grownet      = 'm'
    f_PPnet        = 'm'
    f_algalpeak    = 'm'
    f_zbgc_frac    = 'm'
/

&icefields_drag_nml
    f_drag         = 'x'
    f_Cdn_atm      = 'x'
    f_Cdn_ocn      = 'x'
/

&icefields_fsd_nml
    f_fsdrad       = 'm'
    f_fsdperim     = 'm'
    f_afsd         = 'm'
    f_afsdn        = 'm'
    f_dafsd_newi   = 'x'
    f_dafsd_latg   = 'x'
    f_dafsd_latm   = 'x'
    f_dafsd_wave   = 'x'
    f_dafsd_weld   = 'x'
    f_wave_sig_ht  = 'x'
    f_aice_ww      = 'x'
    f_diam_ww      = 'x'
    f_hice_ww      = 'x'
/
EOF
