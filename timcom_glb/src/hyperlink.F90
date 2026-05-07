module hyperlink
  use timcom_const
  use timcom_grid, only: ocn_panel
  use MPI, only: MPI_OFFSET_KIND
  implicit none

  character(len=8) :: mytag

  integer, pointer :: &
    nx_grid, nx, nxf, &
    ny_grid, ny, nyf, &
    nz_grid, nz, nzf, &
    xu_bgn,  xu_end,  &
    yv_bgn,  yv_end,  daodt, &
    m_comm_cart, myid, np, & 
    m_comm_x, myid_x, npx, &
    m_comm_y, myid_y, npy, &
    daodt_out, daodt_out_unit, &
    daodt_restart,         &
    opt_windmix, opt_arbr_p0, &
    max_iter_p0, opt_da,      &
    opt_solver  ! NTU 20250417 M. Hsieh

  integer, dimension(:), pointer :: &
    nbid, r8type2d, r8type3d, &
    i2type2d, i2type3d, &
    r8type3du, r8type3dv, ndim

  integer(MPI_OFFSET_KIND), dimension(3) :: &
    ns_ijk, nc_ijk, nc_ujk, nc_ivk, nc_ijw

  real(r8), pointer :: &
    dt, odt, &
    total_vol, total_area, &
    avg_glb_t, avg_glb_s, dm0, &
    threshold_bicg, threshold_p0, &
    PCSICriterion  ! NTU 20250417 M. Hsieh

  logical, pointer :: &
    peri_x, peri_y, symm_np, pcsi_cuda_graph

  real(r8), pointer, dimension(:) :: &
    x_grid, y_grid, z_grid, &
    x_face, y_face, z_face, &
    cs, csv, ocs, ocsv, &
    dx, dxu, odx, odxu, &
    dy, dyv, ody, odyv, &
    dz, dzw, odz, odzw, &
    curv_f, tanphi, vbk, hbk

  integer(2), pointer, dimension(:,:) :: &
    kb, iu0, iv0

  integer(2), pointer, dimension(:,:,:) :: &
    in, iu, iv ,iw

  real(r8), pointer, dimension(:,:,:) :: &
      u1, u2, ulf, &
      v1, v2, vlf, &
      t1, t2, tlf, &
      s1, s2, slf, &
      c1, c2, clf, &
      u,   v,   w, &
      p,  ev, add, &
      dmx, dmy, dhx, dhy, &
      rho, sw_trans, volume, &
      stf, smft, vvc, &
      hmean_nudge

  real(r8), pointer, dimension(:,:) :: &
      ar, at, al, ab, ac, &
      cr, ct, cl, cb, cc, &
      cgr, cgrh, cgp, cgph, &
      cgv, cgt,  cgs, cgsh, &
      x,  p0, ssh, u_change, v_change, &
      u_10, v_10, q_10, t_10, pslv, &
      taux, tauy, swup, swdn, &
      lwup, lwdn, lath, senh, &
      qdot, qdot2, melth, snow_f, &
      evap, rain, snow, roff, ioff, &
      ifrc, vice, vsno, melt, salt, area, &
      shf_qsw, tidal_energy, kpp_hblt, &
      qice, aqice, qflux, tgfs, &
      t_nudge, s_nudge, t_da

  real(r8), pointer, dimension(:,:,:,:) :: &
      vdc, kpp_src, hv, t_clim, s_clim

contains

subroutine assign_grid_all(grd)
  implicit none

  type(ocn_panel), target :: grd

  call assign_grid_scalar(grd)
  call assign_grid_vector(grd)

end subroutine assign_grid_all

subroutine assign_grid_scalar(grd)
  implicit none

  type(ocn_panel), target :: grd

  mytag = grd%tag

  nx_grid => grd%nx_grid
  ny_grid => grd%ny_grid
  nz_grid => grd%nz_grid
  nx  => grd%nx
  ny  => grd%ny
  nz  => grd%nz
  nxf => grd%nxf
  nyf => grd%nyf
  nzf => grd%nzf
  xu_bgn => grd%xu_bgn
  xu_end => grd%xu_end
  yv_bgn => grd%yv_bgn
  yv_end => grd%yv_end

  daodt => grd%daodt
  dt    => grd%dt
  odt   => grd%odt
  daodt_out => grd%daodt_out
  daodt_out_unit => grd%daodt_out_unit
  daodt_restart => grd%daodt_restart

  total_vol => grd%total_vol
  total_area => grd%total_area
  avg_glb_t => grd%avg_glb_t
  avg_glb_s => grd%avg_glb_s

  np  => grd%mpicom%np
  npx => grd%mpicom%npx
  npy => grd%mpicom%npy

  peri_x  => grd%mpicom%peri_x
  peri_y  => grd%mpicom%peri_y
  symm_np => grd%mpicom%symm_np

  opt_windmix => grd%opt_windmix
  opt_arbr_p0 => grd%opt_arbr_p0
  max_iter_p0 => grd%max_iter_p0
  opt_da => grd%opt_da
  dm0 => grd%dm0

  threshold_bicg => grd%threshold_bicg
  threshold_p0   => grd%threshold_p0

  opt_solver => grd%opt_solver  ! NTU 20250417 M. Hsieh
  PCSICriterion => grd%PCSICriterion
  pcsi_cuda_graph => grd%pcsi_cuda_graph

end subroutine assign_grid_scalar

subroutine assign_grid_vector(grd)
  implicit none

  type(ocn_panel), target :: grd

  ns_ijk = grd%mpicom%ns_ijk
  nc_ijk = grd%mpicom%nc_ijk
  nc_ujk = grd%mpicom%nc_ujk
  nc_ivk = grd%mpicom%nc_ivk
  nc_ijw = grd%mpicom%nc_ijw

  x_grid => grd%x_grid
  y_grid => grd%y_grid
  z_grid => grd%z_grid
  x_face => grd%x_face
  y_face => grd%y_face
  z_face => grd%z_face

  cs   => grd%cs
  csv  => grd%csv
  ocs  => grd%ocs
  ocsv => grd%ocsv

  dx   => grd%dx
  dxu  => grd%dxu
  odx  => grd%odx
  odxu => grd%odxu

  dy   => grd%dy
  dyv  => grd%dyv
  ody  => grd%ody
  odyv => grd%odyv

  dz   => grd%dz
  dzw  => grd%dzw
  odz  => grd%odz
  odzw => grd%odzw

  curv_f => grd%curv_f
  tanphi => grd%tanphi
  vbk => grd%vbk
  hbk => grd%hbk
  
  kb  => grd%kb
  iu0 => grd%iu0
  iv0 => grd%iv0
  in  => grd%in
  iu  => grd%iu
  iv  => grd%iv
  iw  => grd%iw

  u1  => grd%u1
  u2  => grd%u2
  ulf => grd%ulf
  v1  => grd%v1
  v2  => grd%v2
  vlf => grd%vlf
  t1  => grd%t1
  t2  => grd%t2
  tlf => grd%tlf
  s1  => grd%s1
  s2  => grd%s2
  slf => grd%slf
  c1  => grd%c1
  c2  => grd%c2
  clf => grd%clf
  u   => grd%u
  v   => grd%v
  w   => grd%w
  p   => grd%p
  ev  => grd%ev
  hv  => grd%hv
  add => grd%add
  dmx => grd%dmx
  dmy => grd%dmy
  dhx => grd%dhx
  dhy => grd%dhy
  rho => grd%rho      
  sw_trans => grd%sw_trans

  hmean_nudge => grd%hmean_nudge
  t_clim => grd%t_clim
  s_clim => grd%s_clim
  t_nudge => grd%t_nudge
  s_nudge => grd%s_nudge
  t_da    => grd%t_da

  ar => grd%ar
  at => grd%at
  al => grd%al
  ab => grd%ab
  ac => grd%ac
  cr => grd%cr
  ct => grd%ct
  cl => grd%cl
  cb => grd%cb
  cc => grd%cc
  cgr  => grd%cgr
  cgrh => grd%cgrh
  cgp  => grd%cgp
  cgph => grd%cgph
  cgv  => grd%cgv
  cgt  => grd%cgt
  cgs  => grd%cgs
  cgsh => grd%cgsh
  x  => grd%x
  p0 => grd%p0
  ssh => grd%ssh
  u_change => grd%u_change
  v_change => grd%v_change

  u_10=>grd%u_10
  v_10=>grd%v_10
  q_10=>grd%q_10
  t_10=>grd%t_10
  pslv=>grd%pslv

  taux=>grd%taux
  tauy=>grd%tauy

  swdn=>grd%swdn
  swup=>grd%swup
  lwdn=>grd%lwdn
  lwup=>grd%lwup
  lath=>grd%lath
  senh=>grd%senh
  melth=>grd%melth
  snow_f=>grd%snow_f
  qdot=>grd%qdot
  qdot2=>grd%qdot2

  evap=>grd%evap
  rain=>grd%rain
  snow=>grd%snow
  roff=>grd%roff
  ioff=>grd%ioff
  ifrc=>grd%ifrc
  vice=>grd%vice
  vsno=>grd%vsno
  melt=>grd%melt
  salt=>grd%salt
  qice=>grd%qice
  aqice=>grd%aqice
  qflux=>grd%qflux
  tgfs=>grd%tgfs
  stf=>grd%stf
  smft=>grd%smft
  shf_qsw=>grd%shf_qsw
  tidal_energy=>grd%tidal_energy
  kpp_src=>grd%kpp_src
  kpp_hblt=>grd%kpp_hblt
  vvc=>grd%vvc
  vdc=>grd%vdc

  area=>grd%area
  volume=>grd%volume

  m_comm_cart => grd%mpicom%m_comm_cart
  m_comm_x    => grd%mpicom%m_comm_x
  m_comm_y    => grd%mpicom%m_comm_y

  myid   => grd%mpicom%myid
  myid_x => grd%mpicom%myid_x
  myid_y => grd%mpicom%myid_y
  nbid   => grd%mpicom%nbid
  r8type2d => grd%mpicom%r8type2d
  r8type3d => grd%mpicom%r8type3d
  i2type2d => grd%mpicom%i2type2d
  i2type3d => grd%mpicom%i2type3d
  r8type3du => grd%mpicom%r8type3du
  r8type3dv => grd%mpicom%r8type3dv
  ndim => grd%ndim
end subroutine assign_grid_vector

end module hyperlink
