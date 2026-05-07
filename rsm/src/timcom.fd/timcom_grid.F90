module tai_timcom_grid
  use tai_timcom_const, only: r8
  use tai_timcom_comm,  only: comm_panel
  implicit none

  type :: data_panel
    integer :: ncid, varid
    integer :: nd(3), t(2)
    real(r8), pointer :: time(:), lon(:), lat(:)
    real(r8), pointer :: src(:,:,:), srcm(:,:)
    character(len=256) :: fname
    character(len=16)  :: fvara
    character(len=8)  :: fvadm(3), fdims(3)
  end type data_panel

  type :: ocn_panel
    character(len=8) :: tag
    character(len=128) :: nml_file

    type(comm_panel) :: mpicom

    integer :: &
      nx_grid, nx, nxf, & 
      ny_grid, ny, nyf, &
      nz_grid, nz, nzf, &
      xu_bgn,  xu_end,  &
      yv_bgn,  yv_end,  &
      daodt, ndim(3),   &
      daodt_out,        &
      daodt_out_unit,   &
      daodt_restart,    &
      opt_windmix, opt_arbr_p0, &
      lateral_bc_ntime, &
      max_iter_p0

    real(r8) :: &
      dt, odt, &
      total_vol, total_area, &
      avg_glb_t, avg_glb_s, dm0, &
      threshold_bicg, threshold_p0

    real(r8), pointer, dimension(:) :: &
      x_grid, y_grid, z_grid, &
      x_face, y_face, z_face, &
      cs, csv, ocs, ocsv, &
      dx, dxu, odx, odxu, &
      dy, dyv, ody, odyv, &
      dz, dzw, odz, odzw, &
      curv_f, tanphi, vbk, hbk, &
      lateral_bc_time

    integer(2), pointer, dimension(:,:) :: &
      kb, iu0, iv0

    integer(2), pointer, dimension(:,:,:) :: &
      in, iu, iv ,iw

    integer, pointer, dimension(:,:) :: &
      jwbc, jebc, isbc, inbc

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
      hmean_nudge, t_clim, s_clim, &
      wbc, ebc, sbc, nbc

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
      ifrc, meltw, area, &
      shf_qsw, tidal_energy, kpp_hblt, &
      qice, aqice, qflux, tsea, &
      t_nudge, s_nudge

    real(r8), pointer, dimension(:,:,:,:) :: &
      vdc, kpp_src, hv

  end type ocn_panel

  type :: output_panel
    real(r8), pointer, dimension(:,:,:) :: &
      u2, v2, w, t2, s2
    real(r8), pointer, dimension(:,:) :: &
      p0 
  end type
 
  type(ocn_panel), pointer :: ocn_grid(:)=>null()

contains
subroutine allocate_grid_var(grd)
  implicit none

  type(ocn_panel) :: grd

  call allocate_r8_1d(grd%x_grid, grd%nx)
  call allocate_r8_1d(grd%y_grid, grd%ny)
  call allocate_r8_1d(grd%z_grid, grd%nz)

  call allocate_r8_1d(grd%x_face, grd%nxf)
  call allocate_r8_1d(grd%y_face, grd%nyf)
  call allocate_r8_1d(grd%z_face, grd%nzf)

  call allocate_r8_1d(grd%cs,  grd%ny)
  call allocate_r8_1d(grd%ocs, grd%ny)
  call allocate_r8_1d(grd%csv, grd%nyf)
  call allocate_r8_1d(grd%ocsv,grd%nyf)

  call allocate_r8_1d(grd%dx,  grd%ny)
  call allocate_r8_1d(grd%odx, grd%ny)
  call allocate_r8_1d(grd%dxu, grd%nyf)
  call allocate_r8_1d(grd%odxu,grd%nyf)

  call allocate_r8_1d(grd%dy,  grd%ny)
  call allocate_r8_1d(grd%ody, grd%ny)
  call allocate_r8_1d(grd%dyv, grd%nyf)
  call allocate_r8_1d(grd%odyv,grd%nyf)

  call allocate_r8_1d(grd%dz,  grd%nz)
  call allocate_r8_1d(grd%odz, grd%nz)
  call allocate_r8_1d(grd%dzw, grd%nzf)
  call allocate_r8_1d(grd%odzw,grd%nzf)

  call allocate_r8_1d(grd%curv_f, grd%ny)
  call allocate_r8_1d(grd%tanphi, grd%ny)

  call allocate_i2_2d(grd%kb, grd%nx, grd%ny, 1, 1)
  call allocate_i2_2d(grd%iu0,grd%nxf,grd%ny, 1, 1)
  call allocate_i2_2d(grd%iv0,grd%nx, grd%nyf,1, 1)
  call allocate_i2_3d(grd%in, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_i2_3d(grd%iu, grd%nxf,grd%ny, grd%nz, 1, 1)
  call allocate_i2_3d(grd%iv, grd%nx, grd%nyf,grd%nz, 1, 1)
  call allocate_i2_3d(grd%iw, grd%nx, grd%ny, grd%nzf,1, 1)

  call allocate_r8_2d(grd%area, grd%nx, grd%ny, 0, 0)
  call allocate_r8_3d(grd%volume, grd%nx, grd%ny, grd%nz, 0, 0)

  call allocate_r8_3d(grd%u1, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%u2, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%ulf,grd%nx, grd%ny, grd%nz, 2, 2)

  call allocate_r8_3d(grd%v1, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%v2, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%vlf,grd%nx, grd%ny, grd%nz, 2, 2)

  call allocate_r8_3d(grd%t1, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%t2, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%tlf,grd%nx, grd%ny, grd%nz, 2, 2)

  call allocate_r8_3d(grd%s1, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%s2, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%slf,grd%nx, grd%ny, grd%nz, 2, 2)

  call allocate_r8_3d(grd%c1, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%c2, grd%nx, grd%ny, grd%nz, 2, 2)
  call allocate_r8_3d(grd%clf,grd%nx, grd%ny, grd%nz, 2, 2)
 
  call allocate_r8_3d(grd%u, grd%nxf, grd%ny, grd%nz, 1, 1)
  call allocate_r8_3d(grd%v, grd%nx,  grd%nyf,grd%nz, 1, 1)
  call allocate_r8_3d(grd%w, grd%nx,  grd%ny, grd%nzf,0, 0)
  call allocate_r8_3d(grd%p, grd%nx,  grd%ny, grd%nz, 2, 2)
  call allocate_r8_2d(grd%p0,grd%nx,  grd%ny, 1, 1)
  call allocate_r8_2d(grd%ssh, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%u_change, grd%nxf, grd%ny, 1, 1)
  call allocate_r8_2d(grd%v_change, grd%nx, grd%nyf, 1, 1)

  call allocate_r8_3d(grd%ev,  grd%nx, grd%ny, grd%nz, 0, 0)
  !call allocate_r8_3d(grd%hv,  grd%nx, grd%ny, grd%nz, 0, 0)
  call allocate_r8_1d(grd%hbk, grd%nz)
  call allocate_r8_1d(grd%vbk, grd%nz)
  call allocate_r8_3d(grd%add, grd%nx, grd%ny, grd%nz, 0, 0)

  call allocate_r8_3d(grd%hmean_nudge, grd%nx, grd%ny, grd%nz, 0, 0)
  call allocate_r8_3d(grd%t_clim, grd%nx, grd%ny, 12, 0, 0)
  call allocate_r8_3d(grd%s_clim, grd%nx, grd%ny, 12, 0, 0)
  call allocate_r8_2d(grd%t_nudge, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%s_nudge, grd%nx, grd%ny, 0, 0)

  call allocate_r8_2d(grd%ar,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%at,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%al,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%ab,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%ac,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%cr,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%ct,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%cl,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%cb,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%cc,   grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%x,    grd%nx, grd%ny, 1, 1)
  call allocate_r8_2d(grd%cgr,  grd%nx, grd%ny, 1, 1)
  call allocate_r8_2d(grd%cgrh, grd%nx, grd%ny, 1, 1)
  call allocate_r8_2d(grd%cgp,  grd%nx, grd%ny, 1, 1)
  call allocate_r8_2d(grd%cgph, grd%nx, grd%ny, 1, 1)
  call allocate_r8_2d(grd%cgv,  grd%nx, grd%ny, 1, 1)
  call allocate_r8_2d(grd%cgt,  grd%nx, grd%ny, 1, 1)
  call allocate_r8_2d(grd%cgs,  grd%nx, grd%ny, 1, 1)
  call allocate_r8_2d(grd%cgsh, grd%nx, grd%ny, 1, 1)

  call allocate_r8_2d(grd%u_10, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%v_10, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%q_10, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%t_10, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%pslv, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%taux, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%tauy, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%swup, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%swdn, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%lwup, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%lwdn, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%lath, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%senh, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%qdot, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%qdot2, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%melth, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%snow_f,grd%nx, grd%ny, 0, 0)

  call allocate_r8_2d(grd%evap, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%rain, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%snow, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%roff, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%ioff, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%ifrc, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%meltw,grd%nx, grd%ny, 0, 0)

  call allocate_r8_3d(grd%dhx, grd%nxf, grd%ny, grd%nz, 0, 0)
  call allocate_r8_3d(grd%dhy, grd%nx, grd%nyf, grd%nz, 0, 0)
  call allocate_r8_3d(grd%dmx, grd%nxf, grd%ny, grd%nz, 0, 0)
  call allocate_r8_3d(grd%dmy, grd%nx, grd%nyf, grd%nz, 0, 0)

  call allocate_r8_3d(grd%rho, grd%nx, grd%ny, grd%nz, 0, 0)
  call allocate_r8_3d(grd%sw_trans, grd%nx, grd%ny, grd%nz, 0, 0)

  call allocate_r8_2d(grd%qice,  grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%aqice, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%qflux, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%tsea,  grd%nx, grd%ny, 0, 0)

  call allocate_r8_3d(grd%stf,     grd%nx, grd%ny, 2, 0, 0)
  call allocate_r8_3d(grd%smft,    grd%nx, grd%ny, 2, 0, 0)
  call allocate_r8_2d(grd%shf_qsw, grd%nx, grd%ny, 0, 0)
  call allocate_r8_2d(grd%tidal_energy, grd%nx, grd%ny, 0, 0)


  call allocate_r8_2d(grd%kpp_hblt, grd%nx, grd%ny, 0, 0)  
  call allocate_r8_3d(grd%vvc, grd%nx, grd%ny, grd%nz, 0, 0)

  allocate(grd%hv(grd%nx, grd%ny, grd%nz, 2))
  grd%hv = 0.d0

  allocate(grd%kpp_src(grd%nx, grd%ny, grd%nz, 2))
  grd%kpp_src = 0.d0

  allocate(grd%vdc(grd%nx, grd%ny, 0:grd%nzf, 2))
  grd%vdc = 0.d0

  allocate(grd%wbc(grd%ny, grd%nz, 4))
  grd%wbc = 0.d0
  allocate(grd%ebc(grd%ny, grd%nz, 4))
  grd%ebc = 0.d0
  allocate(grd%sbc(grd%nx, grd%nz, 4))
  grd%sbc = 0.d0
  allocate(grd%nbc(grd%nx, grd%nz, 4))
  grd%nbc = 0.d0

  allocate(grd%jwbc(grd%ny, grd%nz))
  grd%jwbc = 0
  allocate(grd%jebc(grd%ny, grd%nz))
  grd%jebc = 0
  allocate(grd%isbc(grd%nx, grd%nz))
  grd%isbc = 0
  allocate(grd%inbc(grd%nx, grd%nz))
  grd%inbc = 0

end subroutine allocate_grid_var

subroutine allocate_r8_1d(var, ni)
  implicit none

  real(r8), pointer :: var(:)
  integer :: ni

  allocate(var(ni))
  var = 0.d0
end subroutine allocate_r8_1d

subroutine allocate_r8_2d(var, ni, nj, nigh, njgh)
  implicit none

  real(r8), pointer :: var(:,:)
  integer :: ni, nj, nigh, njgh

  allocate(var(1-nigh:ni+nigh, 1-njgh:nj+njgh))
  var = 0.d0
end subroutine allocate_r8_2d

subroutine allocate_r8_3d(var, ni, nj, nk, nigh, njgh)
  implicit none

  real(r8), pointer :: var(:,:,:)
  integer :: ni, nj, nk, nigh, njgh

  allocate(var(1-nigh:ni+nigh, 1-njgh:nj+njgh, nk))
  var = 0.d0
end subroutine allocate_r8_3d

subroutine allocate_i2_2d(var, ni, nj, nigh, njgh)
  implicit none

  integer(2), pointer :: var(:,:)
  integer :: ni, nj, nigh, njgh

  allocate(var(1-nigh:ni+nigh, 1-njgh:nj+njgh))
  var = 0
end subroutine allocate_i2_2d

subroutine allocate_i2_3d(var, ni, nj, nk, nigh, njgh)
  implicit none

  integer(2), pointer :: var(:,:,:)
  integer :: ni, nj, nk, nigh, njgh

  allocate(var(1-nigh:ni+nigh, 1-njgh:nj+njgh, nk))
  var = 0
end subroutine allocate_i2_3d

end module tai_timcom_grid

