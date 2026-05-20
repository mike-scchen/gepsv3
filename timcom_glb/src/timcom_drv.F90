module timcom_drv
  use timcom_const
  use timcom_comm
  use timcom_grid
  use timcom_general
  use hyperlink, only: nx, ny, nz, nxf, nyf, nzf, &
                       xu_bgn, xu_end, yv_bgn, yv_end, &
                       m_comm_cart, myid, r8type3d, &
                       nbid, ndim, symm_np, &
                       kb, in, iu, iv, iw, iu0, iv0, &
                       area, volume, npx, npy, &
                       total_area, total_vol, &
                       avg_glb_t, avg_glb_s, myid_x, myid_y
  implicit none

  real(r8) :: incomp_res

contains

subroutine timcom_check_restart(itf_drv, timer)
  use hyperlink, only: assign_grid_all
  implicit none

  integer, intent(in) :: itf_drv
  type(timer_panel) :: timer
  real(r8) :: dtmp, mjd2

  integer :: itf, n
  character(len=256) :: log_info

  do n = 1, ocn_grid_num
    call assign_grid_all(ocn_grid(n))

    itf =  itf_drv*ocn_grid(n)%daodt/timer%daodt
    dtmp = ocn_grid(n)%dt/day2sec*0.5d0
    mjd2 = dtmp*dble(itf)

    call check_restart(itf, timer%base_mjd+mjd2)

  end do
end subroutine timcom_check_restart

subroutine timcom_exec(itf_drv, timer)
  use hyperlink, only: assign_grid_all
  implicit none

  integer, intent(in) :: itf_drv
  type(timer_panel) :: timer
  real(r8) :: dtmp, mjd1, mjd2, time_stamp, time_spend, dt_cpl

  integer :: itf, it0, mxit, n, ierr
  character(len=256) :: log_info

  if(myid .eq. rootid) then
    write(log_info,'(a,f12.4,a3,f12.4,a,i8)') &
      "[INFO]: DRIVER forecat from:", timer%mjd1, " to", timer%mjd2, ", steps:", itf_drv
    call comm_write_log_info(fid_log, log_info)
  end if

  do n = 1, ocn_grid_num
    call assign_grid_all(ocn_grid(n))
    
    time_stamp = MPI_WTIME()

    it0  = (itf_drv-1)*ocn_grid(n)%daodt/timer%daodt 
    mxit =  itf_drv   *ocn_grid(n)%daodt/timer%daodt
    dtmp = ocn_grid(n)%dt/day2sec*0.5d0
    dt_cpl = day2sec*timer%dt_mjd
    do itf = 1+it0, mxit
      mjd1 = dtmp*dble(itf-1)
      mjd2 = dtmp*dble(itf)  
      
      if(myid .eq. rootid) then
        write(log_info,'(a,f12.4,a3,f12.4,a,i8)') &
          "[INFO]: TIMCOM forecat from:", mjd1, " to", mjd2, ", steps:", itf
        call comm_write_log_info(fid_log, log_info)
      end if
       
      call timcom_fs(itf, it0, timer%syng_mon, mxit, dt_cpl)
      call check_blowup(itf)

#ifdef cpl
  #ifndef owdate
      call check_output_all(itf)
  #else
      call check_output_all(itf, timer%base_mjd+mjd2)
  #endif
#else
      call check_output_all(itf, timer%base_mjd+timer%init_mjd+mjd2)
      call check_restart(itf, timer%base_mjd+timer%init_mjd+mjd2)
      !if(itf .eq. 1800) call comm_finalize
#endif
!      call check_restart(itf, timer%base_mjd+timer%init_mjd+mjd2)

    end do

    time_spend = MPI_WTIME() - time_stamp
    call MPI_ALLREDUCE(time_spend, time_stamp, 1, MPI_REAL8, MPI_MAX, m_comm_cart, ierr)
    if(myid .eq. rootid) then
      write(log_info,'(a,f8.1,a)') &
          "[INFO]: TIMCOM runtime spends ", time_stamp, " sec in a driver loop"
        call comm_write_log_info(fid_log, log_info)
    end if
  end do

end subroutine timcom_exec

subroutine timcom_fs(itf, it0, syng_mon, mxit, dt_cpl)
  use hyperlink, only: u2, v2, t2, s2, p0, u,v,w, u_10, v_10, q_10, t_10, pslv, smft, rho
  implicit none
  integer, intent(in) :: itf, it0, syng_mon, mxit
  real(r8), intent(in) :: dt_cpl
  character(len=256) :: log_info

  call solve_hydrostatic_equation
  call calu_surf_flux
#ifdef cpl
 #ifndef clm_r
  call wind_mixing
 #else
  call wind_mixing(itf)
 #endif
#else
  call wind_mixing(itf)
#endif
  call solve_dynamic_equation
  call impose_surface_flux
  call impose_temp_salt_mixing
  call impose_temp_salt_nudging(syng_mon)
  call impose_bottom_stress
  call impose_trapezoidal_coriolis
  call impose_open_boundary_conditons
  call interp_Agrid_to_Cgrid
  call pressure_solver(itf,it0)
  call solve_continuity_equation
  call check_incompressibility
  call eliminate_arbitrary_pressure
  call interp_Cgrid_to_Agrid
  call modified_filter
  call adjust_seaice_ts(itf,mxit,dt_cpl)

end subroutine timcom_fs

subroutine solve_hydrostatic_equation
  use hyperlink,  only: p, p0, t2, s2, rho, &
                        z_grid, dzw, dz
  use timcom_eos, only: compute_eos_mkcoef_zgrid
  implicit none

  real(r8), pointer :: tp(:,:,:), ts(:,:,:), pp(:,:,:)
  integer(2), pointer :: msk(:,:,:)

  real(r8) :: tmpd, wface(nx,ny,nzf)
  integer :: i, j, k, n, ierr

  tp => t2(1:nx,1:ny,1:nz)
  ts => s2(1:nx,1:ny,1:nz)
  pp =>  p(1:nx,1:ny,1:nz)
  msk => in(1:nx,1:ny,1:nz)

  call compute_eos_mkcoef_zgrid(nx, ny, nz, msk, tp, ts, z_grid, rho) 

  !wface = 0.d0
  !wface(1:nx,1:ny,1) = p0(1:nx,1:ny)

  !do k = 2, nzf
  !  wface(1:nx,1:ny,k) = wface(1:nx,1:ny,k-1) + rho(1:nx,1:ny,k-1)*grav*dz(k-1)
  !end do

  !p(1:nx,1:ny,1:nz) = 0.5d0*(wface(1:nx,1:ny,1:nz) + wface(1:nx,1:ny,2:nzf))

  tmpd = grav*z_grid(1)
  do j = 1, ny
    do i = 1, nx
      p(i,j,1) = p0(i,j) + tmpd*rho(i,j,1)
    end do
  end do

  do k = 2, nz
    tmpd = 0.5d0*grav*dzw(k)
    do j = 1, ny
      do i = 1, nx
        p(i,j,k) = p(i,j,k-1) + tmpd*(rho(i,j,k)+rho(i,j,k-1))
      end do
    end do
  end do

  p = p*in
  call MPI_BARRIER(m_comm_cart, ierr)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, p, symm_np)

end subroutine solve_hydrostatic_equation

subroutine solve_hydrostatic_equation_balance
  use hyperlink,  only: p, p0, t2, s2, rho, &
                        z_grid, dz 
  use timcom_eos, only: compute_eos_wright, &
                        compute_eos_mkcoef, &
                        compute_eos_jacmcd
  implicit none

  real(r8), pointer :: tp(:,:,:), ts(:,:,:), pp(:,:,:)
  integer(2), pointer :: msk(:,:,:)
  real(r8), dimension(nx,ny,nz) :: rho1

  real(r8) :: tmpd, wface(nx,ny), sltd, pisd, res0(1), res1(1), dum0(1), dum1(1)
  integer :: i, j, k, n, opt=2, ierr
  integer, parameter :: &
    opt_wright = 1, &
    opt_mkcoef = 2, &
    opt_JacMcD = 3
  character(len=256) :: log_info

  tp => t2(1:nx,1:ny,1:nz)
  ts => s2(1:nx,1:ny,1:nz)
  pp =>  p(1:nx,1:ny,1:nz)
  msk => in(1:nx,1:ny,1:nz)

  do n = 1, 10

    wface = p0(1:nx,1:ny)
    do k = 1, nz
      pp(1:nx,1:ny,k) = wface + rho(1:nx,1:ny,k)*grav*dz(k)*0.5d0
      wface = wface + rho(1:nx,1:ny,k)*grav*dz(k)
    end do

    select case(opt)
    case(opt_wright)
      call compute_eos_wright(nx, ny, nz, msk, tp, ts, pp, rho1)
    case(opt_mkcoef)
      call compute_eos_mkcoef(nx, ny, nz, msk, tp, ts, pp, rho1)
    case(opt_JacMcd)
      call compute_eos_JacMcd(nx, ny, nz, msk, tp, ts, -1.d-2*z_grid, rho1)
    end select

    res0(1) = sum((rho-dble(msk))**2)
    res1(1) = sum((rho1-rho)**2)

    call MPI_ALLREDUCE(res0, dum0, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
    call MPI_ALLREDUCE(res1, dum1, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

    res1 = dsqrt(dum1)/dsqrt(dum0)

    if(myid .eq. rootid) then
      write(log_info,'(a,i4,f)') "[CHECK RHO]: ", n, res1
      call comm_write_log_info(fid_log, log_info)
    end if
    rho = rho1
    if(res1(1) < 0.001d0) exit
  end do

  p = p*in
  call MPI_BARRIER(m_comm_cart, ierr)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, p , symm_np)

end subroutine solve_hydrostatic_equation_balance 

subroutine calu_surf_flux
  use timcom_phy, only: surface_flux_atmocn
  use hyperlink, only: u2, v2, t2, s2, rho, &
                       u_10, v_10, t_10, q_10, pslv, &
                       swup, swdn, lwup, lwdn, taux, &
                       tauy, senh, lath, evap, rain, &
                       snow, roff, ifrc, melth, melt, salt, &
                       stf, shf_qsw, smft
  implicit none

  integer :: nsize
  integer(2) :: mask(nx,ny)
  real(r8), dimension(nx,ny) :: &
    u_ocn, v_ocn, t_ocn, &
    qdot, wflux, qdot2, snow_f
  integer :: i, j  
  real(r8) :: tmp
  real(r8), parameter :: shr_const_latice  = 3.337d5

  nsize = nx*ny

  mask  = in(1:nx,1:ny,1)
  u_ocn = u2(1:nx,1:ny,1)/100.d0
  v_ocn = v2(1:nx,1:ny,1)/100.d0
  t_ocn = t2(1:nx,1:ny,1) + 273.15d0

  call surface_flux_atmocn(nsize, mask, ifrc, & 
                      u_ocn, v_ocn, t_ocn, &
                      u_10, v_10, t_10, q_10, pslv,  &
                      taux, tauy, senh, lath, lwup, evap)

  snow_f = snow*shr_const_latice
  qdot  = (1.d0-ifrc)*(lwup + lwdn + senh + lath - snow_f) + melth
  wflux = (1.d0-ifrc)*(rain + snow + evap + roff) + melt
  qdot2 = (1.d0-ifrc)*(swup + swdn)

  do j = 1, ny
    do i = 1, nx
      stf(i,j,1)   = dble(in(i,j,1))*qdot(i,j)*hflux_factor
      stf(i,j,2)   = dble(in(i,j,1))*(wflux(i,j)*(-34.7d0)*1.d-4 + salt(i,j)*0.1d0)
      shf_qsw(i,j) = dble(in(i,j,1))*qdot2(i,j)*hflux_factor
      smft(i,j,1)  = dble(in(i,j,1))*taux(i,j)
      smft(i,j,2)  = dble(in(i,j,1))*tauy(i,j)
    end do
  end do 
end subroutine calu_surf_flux

#ifdef cpl
 #ifndef clm_r
  subroutine wind_mixing
 #else
  subroutine wind_mixing(itf)
 #endif
#else
 subroutine wind_mixing(itf)
#endif
  use timcom_kpp, only: kppglo
  use timcom_phy, only: pp82
  use hyperlink, only: w, u2, v2, t2, s2, rho, &
                       sw_trans, stf, shf_qsw, &
                       smft, vbk, hbk, add, vdc, &
                       vvc, ev, hv, kpp_src, kpp_hblt, &
                       opt_windmix
  implicit none

  integer, parameter :: &
    opt_kpp  = 1, &
    opt_pp82 = 2

  integer :: ierr, itf

  select case(opt_windmix)
  case(opt_kpp) 
#ifdef cpl
 #ifndef clm_r
    call kppglo(u2=u2, v2=v2, t2=t2, s2=s2, &
                rho=rho, trans=sw_trans, stf=stf, shf_qsw=shf_qsw, smft=smft, &
                vbk=vbk, hbk=hbk, add=add, vdc=vdc, vvc=vvc, &
                ev=ev, hv=hv, kpp_src=kpp_src, kpp_hblt=kpp_hblt)
 #else
    call kppglo(itf, u2=u2, v2=v2, t2=t2, s2=s2, &                                 
                rho=rho, trans=sw_trans, stf=stf, shf_qsw=shf_qsw, smft=smft, &    
                vbk=vbk, hbk=hbk, add=add, vdc=vdc, vvc=vvc, &                     
                ev=ev, hv=hv, kpp_src=kpp_src, kpp_hblt=kpp_hblt)                  
 #endif
#else
    call kppglo(itf, u2=u2, v2=v2, t2=t2, s2=s2, &
                rho=rho, trans=sw_trans, stf=stf, shf_qsw=shf_qsw, smft=smft, &
                vbk=vbk, hbk=hbk, add=add, vdc=vdc, vvc=vvc, &
                ev=ev, hv=hv, kpp_src=kpp_src, kpp_hblt=kpp_hblt)
#endif  
  case(opt_pp82)
    call pp82(nx, ny, nz, w, u2(1:nx,1:ny,1:nz), v2(1:nx,1:ny,1:nz), rho, hbk, ev, hv)
  
  end select

  call MPI_BARRIER(m_comm_cart, ierr)
end subroutine wind_mixing

subroutine solve_dynamic_equation
  use hyperlink, only: u, v, w, p,  &
                       u1, ulf, u2, &
                       v1, vlf, v2, &
                       t1, tlf, t2, &
                       s1, slf, s2, &
                       dmx, dmy, ev, &
                       dhx, dhy, hv         
  implicit none

  real(r8) :: p_grad(nx,ny,nz)

  call pressure_gradient_x(p, p_grad) 
  call solve_momentun_equation(u, v, w, p_grad, u1, ulf, u2, dmx, dmy, ev) 

  call pressure_gradient_y(p, p_grad)
  call solve_momentun_equation(u, v, w, p_grad, v1, vlf, v2, dmx, dmy, ev)

!  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanTX =', avg_glb_t) 
!  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanSX =', avg_glb_s)
  
  call solve_tracer_equation(u, v, w, t1, tlf, t2, dhx, dhy, hv(:,:,:,1))

  call solve_tracer_equation(u, v, w, s1, slf, s2, dhx, dhy, hv(:,:,:,2))

!  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanT1 =', avg_glb_t)
!  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanS1 =', avg_glb_s)

 !call solve_tracer_equation(u, v, w, c1, clf, c2, dhx, dhy, hv(:,:,:,1))
end subroutine solve_dynamic_equation

subroutine pressure_gradient_x(p, px)
  use hyperlink, only: odx
  implicit none

  real(r8), intent(in) :: &
    p(-1:nx+2,-1:ny+2,1:nz)

  real(r8), intent(out) :: &
    px(nx,ny,nz)

  integer :: i, j, k
  real(r8), dimension(0:nxf+1,0:ny+1) :: xp, msk
  real(r8) :: tmp

  px = 0.d0

  do k = 1, nz
    xp = 0.d0
    msk = dble(iu(:,:,k))
    do j = 1, ny
      do i = xu_bgn-1, xu_end+1
        xp(i,j) = msk(i,j)*(p(i,j,k)-p(i-1,j,k))
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        tmp = o12*msk(i,j)*msk(i+1,j)
        px(i,j,k) = 0.5d0*(xp(i+1,j) + xp(i,j))  &
                   + tmp*(-xp(i-1,j) + xp(i,j)   &
                          +xp(i+1,j) - xp(i+2,j))
      end do
      px(:,j,k) = px(:,j,k)*odx(j)
    end do
  end do

  px = px/rho_sw
end subroutine pressure_gradient_x

subroutine pressure_gradient_y(p, py)
  use hyperlink, only: ocs, ody
  implicit none

  real(r8), intent(in) :: &
    p(-1:nx+2,-1:ny+2,1:nz)

  real(r8), intent(out) :: &
    py(nx,ny,nz)

  integer :: i, j, k
  real(r8), dimension(0:nx+1,0:nyf+1) :: yp, msk
  real(r8) :: tmp

  py = 0.d0

  do k = 1, nz
    yp = 0.d0
    msk = dble(iv(:,:,k))
    do j = yv_bgn-1, yv_end+1
      do i = 1, nx
        yp(i,j) = msk(i,j)*(p(i,j,k)-p(i,j-1,k))
      end do
    end do

    do j = 1, ny
      do i = 1, nx
        tmp = o12*msk(i,j)*msk(i,j+1)
        py(i,j,k) = 0.5d0*(yp(i,j+1) + yp(i,j)) &
                  +  tmp*(-yp(i,j-1) + yp(i,j)  &
                          +yp(i,j+1) - yp(i,j+2))
      end do
      py(:,j,k) = py(:,j,k)*ody(j)
    end do
  end do

  py = py/rho_sw
end subroutine pressure_gradient_y

subroutine solve_momentun_equation(u, v, w, dp, q1, qlf, q2, hmix_x, hmix_y, vmix_z)
  use hyperlink, only: dt 
  implicit none

  real(r8), intent(inout) ::   &
    q2(-1:nx+2,-1:ny+2,1:nz)

  real(r8), intent(in) ::      &
    u(0:nxf+1, 0:ny+1,  1:nz), &
    v(0:nx+1,  0:nyf+1, 1:nz), &
    w(1:nx,    1:ny,    1:nzf),&
    dp ( 1:nx,   1:ny,  1:nz), &
    q1 (-1:nx+2,-1:ny+2,1:nz), &
    qlf(-1:nx+2,-1:ny+2,1:nz), &
    hmix_x(nxf,ny, nz),     &
    hmix_y(nx ,nyf,nz),     &
    vmix_z(nx, ny, nz)

  real(r8), dimension(nx,ny) :: &
    dqx_dx, dqy_dy, dqz_dz, qz
 
  real(r8) :: dtin 
  integer :: i, j, k

  dqx_dx = 0.d0
  dqy_dy = 0.d0
  dqz_dz = 0.d0

  select case(lfsrf)
  case(0)
    qz(1:nx,1:ny) = 0.d0
  case(1)
    qz(1:nx,1:ny) = w(1:nx,1:ny,1)*q2(1:nx,1:ny,1)
  end select

  do k = 1, nz
    call compute_vertical_fluxes(k, w, q1, qlf, q2, vmix_z(:,:,k), qz, dqz_dz)
    call compute_longitudinal_fluxes(k, u, q1(:,:,k), q2(:,:,k), hmix_x(:,:,k), dqx_dx)
    call compute_latitudinal_fluxes (k, v, q1(:,:,k), q2(:,:,k), hmix_y(:,:,k), dqy_dy)

    do j = 1, ny
      do i = 1, nx
        dtin = dt*dble(in(i,j,k))
        q2(i,j,k) = q1(i,j,k) - dtin*(dp(i,j,k) + dqx_dx(i,j) + dqy_dy(i,j) + dqz_dz(i,j)) 
      end do
    end do
  end do

end subroutine solve_momentun_equation

subroutine solve_tracer_equation(u, v, w, q1, qlf, q2, hmix_x, hmix_y, vmix_z)
  use hyperlink, only: dt
  implicit none

  real(r8), intent(inout) :: &
    q2(-1:nx+2,-1:ny+2,1:nz)

  real(r8), intent(in) ::      &
    u(0:nxf+1, 0:ny+1,  1:nz), &
    v(0:nx+1,  0:nyf+1, 1:nz), &
    w(1:nx,    1:ny,    1:nzf),&
    q1 (-1:nx+2,-1:ny+2,1:nz), &
    qlf(-1:nx+2,-1:ny+2,1:nz), &
    hmix_x(nxf,ny, nz),     &
    hmix_y(nx ,nyf,nz),     &
    vmix_z(nx, ny, nz)

  real(r8), dimension(nx,ny) :: &
    dqx_dx, dqy_dy, dqz_dz, qz

  real(r8) :: dtin
  integer :: i, j, k

  dqx_dx = 0.d0
  dqy_dy = 0.d0
  dqz_dz = 0.d0

  select case(lfsrf)
  case(0)
    qz(1:nx,1:ny) = 0.d0
  case(1)
    qz(1:nx,1:ny) = w(1:nx,1:ny,1)*q2(1:nx,1:ny,1)
  end select

  do k = 1, nz
    call compute_vertical_fluxes(k, w, q1, qlf, q2, vmix_z(:,:,k), qz, dqz_dz)
    call compute_longitudinal_fluxes(k, u, q1(:,:,k), q2(:,:,k), hmix_x(:,:,k), dqx_dx)
    call compute_latitudinal_fluxes (k, v, q1(:,:,k), q2(:,:,k), hmix_y(:,:,k), dqy_dy)

    do j = 1, ny
      do i = 1, nx
        dtin = dt*dble(in(i,j,k))
        q2(i,j,k) = q1(i,j,k) - dtin*(dqx_dx(i,j) + dqy_dy(i,j) + dqz_dz(i,j))
      end do
    end do
  end do

end subroutine solve_tracer_equation

subroutine compute_vertical_fluxes(k, w, q1, qlf, q2, vmix_z, qz_t, dqz_dz)
  use hyperlink, only: odz
  implicit none

  integer, intent(in) :: &
    k

  real(r8), intent(in) :: &
    w(1:nx,1:ny,1:nzf),   &
    q1 (-1:nx+2,-1:ny+2,nz), &
    qlf(-1:nx+2,-1:ny+2,nz), &
    q2 (-1:nx+2,-1:ny+2,nz), &
    vmix_z(nx,ny)

  real(r8), intent(inout) :: &
    qz_t(nx,ny)

  real(r8), intent(out) :: &
    dqz_dz(nx,ny)

  real(r8) :: qz_b(nx,ny), scr(nx,ny), tmp

  integer :: i, j, l

  l = k + 1

  if(k .eq. nz) then
    qz_b = 0.d0
  else
    scr(:,:) = 6.d0*(q2(1:nx,1:ny,k) + q2(1:nx,1:ny,l))
    if(k .gt. 1 .and. k .lt. nz-1) then
      do j = 1, ny
        do i = 1, nx
          tmp = dble(iw(i,j,k)*iw(i,j,l+1))
          scr(i,j) = scr(i,j) + tmp*(-qlf(i,j,k-1)+q2(i,j,k)+q2(i,j,l)-q2(i,j,l+1))
        end do
      end do
    end if

    scr = o12*scr

    do j = 1, ny
      do i = 1, nx
        tmp = dble(iw(i,j,l))
        qz_b(i,j) = tmp*(w(i,j,l)*scr(i,j) - vmix_z(i,j)*(q1(i,j,l)-q1(i,j,k)))
      end do
    end do
  end if

  dqz_dz = (qz_b - qz_t)*odz(k)
  qz_t  = qz_b

end subroutine compute_vertical_fluxes

subroutine compute_latitudinal_fluxes(k, v, q1, q2, hmix_y, dqy_dy)
  use hyperlink, only: csv, odyv, ocs, ody
  implicit none

  integer, intent(in) :: &
    k

  real(r8), intent(in) :: &
    v(0:nx+1,0:nyf+1,1:nz), &
    q1(-1:nx+2,-1:ny+2),  &
    q2(-1:nx+2,-1:ny+2),  &
    hmix_y(nx,nyf)

  real(r8), intent(out) :: &
    dqy_dy(nx,ny)

  real(r8) :: qy(nx,nyf), scr(nx,nyf), tmp, odyj
  integer :: i, j
  real(r8) :: msk(0:nx+1,0:nyf+1)

  scr = 0.d0
  qy  = 0.d0

  msk = dble(iv(:,:,k))
  scr(1:nx,yv_bgn:yv_end) = 6.d0*(q2(1:nx,yv_bgn-1:yv_end-1) + q2(1:nx,yv_bgn:yv_end))
  
  do j = yv_bgn, yv_end
    do i = 1, nx
      tmp = msk(i,j-1)*msk(i,j+1)
      scr(i,j) = scr(i,j) + tmp*(-q2(i,j-2) + q2(i,j-1) + q2(i,j) - q2(i,j+1))
    end do
  end do

  scr = o12*scr

  do j = yv_bgn, yv_end
    do i = 1, nx
      tmp = csv(j)*msk(i,j)
      qy(i,j) = tmp*(v(i,j,k)*scr(i,j) - hmix_y(i,j)*odyv(j)*(q1(i,j) - q1(i,j-1)))
    end do
  end do

  do j = 1, ny
    odyj = ocs(j)*ody(j)
    do i = 1, nx
      dqy_dy(i,j) = (qy(i,j+1) - qy(i,j))*odyj
    end do
  end do

end subroutine compute_latitudinal_fluxes

subroutine compute_longitudinal_fluxes(k, u, q1, q2, hmix_x, dqx_dx)
  use hyperlink, only: odx
  implicit none

  integer, intent(in) :: &
    k

  real(r8), intent(in) :: &
    u(0:nxf+1,0:ny+1,1:nz), &
    q1(-1:nx+2,-1:ny+2),  &
    q2(-1:nx+2,-1:ny+2),  &
    hmix_x(nxf,ny)

  real(r8), intent(out) :: &
    dqx_dx(nx,ny)

  real(r8) :: qx(nxf,ny), scr(nxf,ny), tmp
  integer :: i, j
  real(r8) :: msk(0:nxf+1,0:ny+1)

  scr = 0.d0
  qx  = 0.d0

  msk = dble(iu(:,:,k))
  scr(xu_bgn:xu_end,1:ny) = 6.d0*(q2(xu_bgn-1:xu_end-1,1:ny) + q2(xu_bgn:xu_end,1:ny))
  
  do j = 1, ny
    do i = xu_bgn, xu_end
      tmp = msk(i-1,j)*msk(i+1,j)
      scr(i,j) = scr(i,j) + tmp*(-q2(i-2,j) + q2(i-1,j) + q2(i,j) - q2(i+1,j))
    end do
  end do

  scr = o12*scr

  do j = 1, ny
    do i = xu_bgn, xu_end
      tmp = msk(i,j)
      qx(i,j) = tmp*(u(i,j,k)*scr(i,j) - hmix_x(i,j)*odx(j)*(q1(i,j) - q1(i-1,j)))
    end do
  end do

  do j = 1, ny
    do i = 1, nx
      dqx_dx(i,j) = (qx(i+1,j) - qx(i,j))*odx(j)
    end do
  end do
end subroutine compute_longitudinal_fluxes

subroutine impose_surface_flux
  use timcom_phy, only: wind_surf, heat_surf, heat_volume, salt_surf
  use hyperlink, only: dt, odz, u2, v2, t2, s2, &
                       smft, stf, sw_trans, shf_qsw
  implicit none

  real(r8) :: qvol(nx,ny,nz), avg_qvol
  integer :: k
  ! ===============
  ! surface effect
  ! ===============
  ! A.wind stress
  call wind_surf(nx, ny, smft, dt, odz(1), u2(1:nx,1:ny,1), v2(1:nx,1:ny,1))

  ! B.heat fluxes, qdot
  ! (longwave + latent + sensible), qdot (positive downward)
  call heat_surf(nx, ny, stf(:,:,1), dt, odz(1), t2(1:nx,1:ny,1))

  ! C.short wave radiation, qdot2 (positive downward)
  call heat_volume(nx, ny, nz, sw_trans, shf_qsw, dt, odz, t2(1:nx,1:ny,1:nz))

  ! D.evaporation and precipitation (no salt sources)
  call salt_surf(nx, ny, stf(:,:,2), dt, odz(1), s2(1:nx,1:ny,1))

  do k = 1, nz
    qvol(:,:,k) = sw_trans(:,:,k)*shf_qsw
  end do
  avg_qvol = 0.d0

!  call go_check_conservation(qvol, volume, total_vol, 'glob_meanQSW = ', avg_qvol)
!  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanT2 =', avg_glb_t)
!  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanS2 =', avg_glb_s)
end subroutine impose_surface_flux

subroutine impose_temp_salt_mixing
  use timcom_kpp, only: kppsrcglo
  use hyperlink,  only: t2, s2, kpp_src, opt_windmix
  implicit none

  integer, parameter :: &
    opt_kpp  = 1, &
    opt_pp82 = 2

  integer :: ierr

  select case(opt_windmix)
  case(opt_kpp)
    call kppsrcglo(t2(1:nx,1:ny,1:nz), s2(1:nx,1:ny,1:nz), kpp_src(1:nx,1:ny,1:nz,1:2))
  case(opt_pp82)
    continue
  end select 

!  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanT3 =', avg_glb_t)
!  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanS3 =', avg_glb_s)
end subroutine impose_temp_salt_mixing

subroutine impose_temp_salt_nudging(syng_mon)
  use hyperlink, only: t2, s2, t_nudge, s_nudge, t_clim, s_clim, t_da, opt_da
  implicit none

  integer, intent(in) :: syng_mon

  if(opt_da .eq. 1) &
    t2(1:nx,1:ny,1) = t2(1:nx,1:ny,1) + t_nudge(:,:)*(t_da(:,:) - t2(1:nx,1:ny,1))
!  t2(1:nx,1:ny,1) = t2(1:nx,1:ny,1) + t_nudge(:,:)*(t_clim(:,:,syng_mon) - t2(1:nx,1:ny,1))
!  s2(1:ny,1:ny,1) = s2(1:nx,1:ny,1) + s_nudge(:,:)*(s_clim(:,:,syng_mon) - s2(1:nx,1:ny,1))

!  t2(1:nx,1:ny,2:nz) = (1.d0-heman_nudge(1:nx,1:ny,2:nz))*t2(1:nx,1:ny,2:nz) &
!                     + hmean_nudge(1:nx,1:ny,2:nz)*t_clim(1:nx,1:ny,2:nz)
!  s2(1:nx,1:ny,2:nz) = (1.d0-heman_nudge(1:nx,1:ny,2:nz))*s2(1:nx,1:ny,2:nz) &
!                     + hmean_nudge(1:nx,1:ny,2:nz)*s_clim(1:nx,1:ny,2:nz)

end subroutine impose_temp_salt_nudging

subroutine impose_bottom_stress
  use hyperlink, only: dt, odz, u1, u2, v1, v2
  implicit none

  integer :: i, j, k
  real(r8) :: tmp, drg
  real(r8), parameter :: drag = 1.d-3

  tmp = dt*drag

  do j = 1, ny
    do i = 1, nx
      k = kb(i,j)
      drg = in(i,j,k)*tmp*odz(k)*dsqrt(u1(i,j,k)**2+v1(i,j,k)**2)
      drg = min(0.2d0, drg)
      u2(i,j,k) = u2(i,j,k) - drg*u1(i,j,k)
      v2(i,j,k) = v2(i,j,k) - drg*v1(i,j,k)
    end do
  end do
end subroutine impose_bottom_stress

subroutine impose_trapezoidal_coriolis
  use hyperlink, only: dt, curv_f, tanphi, &
                       u1, ulf, u2, v1, v2
  implicit none

  integer i, j, k
  real(r8) :: tmp, temp, qu, qv

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        tmp = 0.5d0*(curv_f(j) + ulf(i,j,k)*tanphi(j))*dt
        temp = 1.d0/(1.d0 + tmp**2)
        qu = u2(i,j,k) + tmp*v1(i,j,k)
        qv = v2(i,j,k) - tmp*u1(i,j,k)
        u2(i,j,k) = temp*(qu + tmp*qv)
        v2(i,j,k) = temp*(qv - tmp*qu)
      end do
    end do
  end do
end subroutine impose_trapezoidal_coriolis

subroutine impose_open_boundary_conditons
  use hyperlink, only: u2, v2, t2, s2, u, v, &
                       ulf, vlf, tlf, slf, dt, &
                       odx, ody, ocs, csv, &
                       myid_x, myid_y, peri_x, peri_y
  implicit none

  integer :: i, j, k, ierr
  real(r8) :: tmp, tmpin, temp, temp1, temp2

  u2(1:nx,1:ny,1:nz) = u2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)
  v2(1:nx,1:ny,1:nz) = v2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)
  t2(1:nx,1:ny,1:nz) = t2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)
  s2(1:nx,1:ny,1:nz) = s2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)

  if(.not. peri_x) then
    tmp = 0.5d0*dt 
    if(myid_x .eq. 0) then
      do k = 1, nz
        do j = 1, ny
          tmpin = in(1,j,k)*tmp*odx(j)
          temp  = dabs(u(1,j,k))
          temp1 = tmpin*(temp + u(1,j,k))
          temp2 = tmpin*(temp - u(1,j,k))
          u2(1,j,k) = u2(1,j,k) + temp1*u2(0,j,k) - temp2*ulf(1,j,k)
          v2(1,j,k) = v2(1,j,k) + temp1*v2(0,j,k) - temp2*vlf(1,j,k)
          t2(1,j,k) = t2(1,j,k) + temp1*t2(0,j,k) - temp2*tlf(1,j,k)
          s2(1,j,k) = s2(1,j,k) + temp1*s2(0,j,k) - temp2*slf(1,j,k)
        end do
      end do
    end if

    if(myid_x .eq. npx-1) then
      do k = 1, nz
        do j = 1, ny
          tmpin = in(nx,j,k)*tmp*odx(j)
          temp  = dabs(u(nxf,j,k))
          temp1 = tmpin*(temp + u(nxf,j,k))
          temp2 = tmpin*(temp - u(nxf,j,k))
          u2(nx,j,k) = u2(nx,j,k) - temp1*ulf(nx,j,k) + temp2*u2(nx+1,j,k)
          v2(nx,j,k) = v2(nx,j,k) - temp1*vlf(nx,j,k) + temp2*v2(nx+1,j,k)
          t2(nx,j,k) = t2(nx,j,k) - temp1*tlf(nx,j,k) + temp2*t2(nx+1,j,k)
          s2(nx,j,k) = s2(nx,j,k) - temp1*slf(nx,j,k) + temp2*s2(nx+1,j,k)
        end do
      end do
    end if
  end if

  if(.not. peri_y) then
    if(myid_y .eq. 0) then
      tmp = 0.5d0*dt*ody(1)*csv(1)*ocs(1)
      do k = 1, nz
        do i = 1, nx
          tmpin = in(i,1,k)*tmp
          temp  = dabs(v(i,1,k))
          temp1 = tmpin*(temp + v(i,1,k))
          temp2 = tmpin*(temp - v(i,1,k))
          u2(i,1,k) = u2(i,1,k) + temp1*u2(i,0,k) - temp2*ulf(i,1,k)
          v2(i,1,k) = v2(i,1,k) + temp1*v2(i,0,k) - temp2*vlf(i,1,k)
          t2(i,1,k) = t2(i,1,k) + temp1*t2(i,0,k) - temp2*tlf(i,1,k)
          s2(i,1,k) = s2(i,1,k) + temp1*s2(i,0,k) - temp2*slf(i,1,k)
        end do
      end do 
    end if
 
    if(myid_y .eq. npy-1 .and. .not. symm_np) then
      tmp = 0.5d0*dt*ody(ny)*csv(nyf)*ocs(ny)
      do k = 1, nz
        do i = 1, nx
          tmpin = in(i,ny,k)*tmp
          temp  = dabs(v(i,nyf,k))
          temp1 = tmpin*(temp + v(i,nyf,k))
          temp2 = tmpin*(temp - v(i,nyf,k))
          u2(i,ny,k) = u2(i,ny,k) - temp1*ulf(i,ny,k) + temp2*u2(i,ny+1,k)
          v2(i,ny,k) = v2(i,ny,k) - temp1*vlf(i,ny,k) + temp2*v2(i,ny+1,k)
          s2(i,ny,k) = s2(i,ny,k) - temp1*slf(i,ny,k) + temp2*s2(i,ny+1,k)
          t2(i,ny,k) = t2(i,ny,k) - temp1*tlf(i,ny,k) + temp2*t2(i,ny+1,k)
         end do
       end do
     end if
  end if

  call MPI_BARRIER(m_comm_cart, ierr)

  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, u2, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, v2, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, t2, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, s2, symm_np)

  if(symm_np) v2(:,ny+1:ny+2,:) = -v2(:,ny+1:ny+2,:)
  
end subroutine impose_open_boundary_conditons

subroutine interp_Agrid_to_Cgrid
  use hyperlink, only: u, u2, v, v2, &
                       m_comm_cart, r8type3du, r8type3dv, &
                       nbid, ndim, symm_np
  implicit none

  integer :: i, j, k, ierr
  real(r8) :: tmp, tmp1
  real(r8) :: msk(0:nxf+1,0:nyf+1)

  msk = 0
  ! LONGITUDINAL DIRECTION
  u(xu_bgn:xu_end,1:ny,1:nz) = 6.0d0*(u2(xu_bgn-1:xu_end-1,1:ny,1:nz) + u2(xu_bgn:xu_end,1:ny,1:nz))
  do k = 1, nz
    msk(0:nxf+1,1:ny) = iu(0:nxf+1,1:ny,k)
    do j = 1, ny
      do i = xu_bgn, xu_end
        tmp1 = msk(i-1,j)*msk(i+1,j)
        tmp  = (-u2(i-2,j,k) + u2(i-1,j,k) + u2(i,j,k) - u2(i+1,j,k))*tmp1
        u(i,j,k) =( u(i,j,k) + tmp)*msk(i,j)
      end do
    end do
  end do
  u(xu_bgn:xu_end,1:ny,1:nz) = o12*u(xu_bgn:xu_end,1:ny,1:nz)

  msk = 0
  v(1:nx,yv_bgn:yv_end,1:nz) = 6.d0*(v2(1:nx,yv_bgn-1:yv_end-1,1:nz) + v2(1:nx,yv_bgn:yv_end,1:nz))
  do k = 1, nz
    msk(1:nx,0:nyf+1) = iv(1:nx,0:nyf+1,k)
    do j = yv_bgn, yv_end
      do i = 1, nx
        tmp1 = msk(i,j-1)*msk(i,j+1)
        tmp  = (-v2(i,j-2,k) + v2(i,j-1,k) + v2(i,j,k) - v2(i,j+1,k))*tmp1
        v(i,j,k) = (v(i,j,k) + tmp)*msk(i,j)
      end do
    end do
  end do
  v(1:nx,yv_bgn:yv_end,1:nz) = o12*v(1:nx,yv_bgn:yv_end,1:nz)

  call MPI_BARRIER(m_comm_cart, ierr)
  call mpi_exch_r8type3du(m_comm_cart, r8type3du, nbid, ndim, u, symm_np)
  call mpi_exch_r8type3dv(m_comm_cart, r8type3dv, nbid, ndim, v, symm_np)
  if(symm_np) v(:,nyf+1,:) = -v(:,nyf+1,:)

end subroutine interp_Agrid_to_Cgrid

subroutine solve_continuity_equation
! ======================================== !
! CALCULATE W BY CONTINUITY EQUATION       !
! ======================================== !
! this completes the advanced time level advection velocity
! having exactly zero barotropic divergence and satisfying all
! kinematic and specified open boundary conditions
! in hydrostatic assumption
  use hyperlink, only: u, v, w, odx, ody, odz, ocs, csv
  implicit none

  integer :: i, j, k
  real(r8) :: tmp ,temp

  w = 0.d0

  select case(lfsrf)
  case(0)
    do k = 1, nz             !(k:layer& top   face for w, k+1:bottom face for w)
      tmp = 1.d0/odz(k)
      do j = 1, ny           !(j:cell & north face for v)
        temp = ocs(j)*ody(j)
        do i = 1, nx         !(i:cell & east  face for u)
          w(i,j,k+1) = (w(i,j,k)  &
                     - ((u(i+1,j,k)-u(i,j,k))*odx(j) &
                     + (csv(j+1)*v(i,j+1,k)-csv(j)*v(i,j,k))*temp)*tmp)
        end do
      end do
    end do

  case(1)
    do j = 1, ny              !(j:cell & north face for v)
      temp = ocs(j)*ody(j)
      do i = 1, nx            !(i:cell & east  face for u)
        do k = kb(i,j), 1, -1
          tmp = 1.d0/odz(k)
          w(i,j,k) = w(i,j,k+1) &
                   + ((u(i+1,j,k)-u(i,j,k))*odx(j) &
                   + (csv(j+1)*v(i,j+1,k)-csv(j)*v(i,j,k))*temp)*tmp
        end do
      end do
    end do
    w(1:nx,1:ny,1:nzf) = w(1:nx,1:ny,1:nzf)*iw(1:nx,1:ny,1:nzf)
  end select
end subroutine solve_continuity_equation

subroutine pressure_solver(itf,it0)

  use timcom_solver, only: p_bicgstab_le, PCSI, &  ! NTU 20250417 M. Hsieh
                           solver_BiCGStab, solver_PCSI

  use hyperlink, only: p0, ssh, odz, u_change, v_change, &
                       ab, al, ac, ar, at, x, &
                       cb, cl, cc, cr, ct,    &
                       cgr, cgrh, cgp, cgv,   &
                       cgs, cgt, cgph, cgsh,  &
                       u, odx, v, odyv, w, odt,  &
                       m_comm_cart, symm_np,  & 
                       r8type3du, r8type3dv,  &
                       nbid, ndim, myid_x, myid_y, &
                       threshold_p0, max_iter_p0, daodt_out_unit, &
                       opt_solver  ! NTU 20250417 M. Hsieh

  implicit none

  integer, intent(in) :: itf, it0

  real(r8) :: s(nx,ny), scr, max_vel, tp, tmp1, tmp2, res_p2, res_p3
  integer :: i, j, k, ierr, itmask, iter
  character(len=256) :: log_info

  u_change = u(:,:,1)
  v_change = v(:,:,1)
  call MPI_BARRIER(m_comm_cart, ierr)
  
  do itmask = 1, max_iter_p0

    call solve_continuity_equation

    select case(lfsrf)
    case(0)
      do j = 1, ny
        do i = 1, nx
          s(i,j) = w(i,j,kb(i,j)+1)
        end do
      end do
    case(1)
      s(:,:) = -w(:,:,1)
    end select

    if ( opt_solver == solver_BiCGStab )         &  ! NTU 20250417 M. Hsieh
    call p_bicgstab_le(ab, al, ac, ar, at, s, x, &
                       cb, cl, cc, cr, ct,       &
                       cgr, cgrh, cgp, cgv,      &
                       cgs, cgt, cgph, cgsh, iter, res_p2, res_p3)

    if ( opt_solver == solver_PCSI )  &  ! NTU 20250417 M. Hsieh
       call PCSI(ab, al, ac, ar, at,  &
                 cb, cl, cc, cr, ct, s, x, iter, res_p2)

    p0 = p0 + odt*x*rho_sw

    max_vel = 0.d0
    do j = 1, ny
      do i = xu_bgn, xu_end
        scr = -(x(i,j)-x(i-1,j))*odx(j)
        u(i,j,1:nz) = u(i,j,1:nz) + scr*iu(i,j,1:nz)
        max_vel = max(max_vel, dble(1-iu(i,j,1))*dabs(scr))
      end do
    end do

    do j = yv_bgn, yv_end
      do i = 1, nx
        scr = -(x(i,j)-x(i,j-1))*odyv(j)
        v(i,j,1:nz) = v(i,j,1:nz) + scr*iv(i,j,1:nz)
        max_vel = max(max_vel, dble(1-iv(i,j,1))*dabs(scr))
      end do
    end do

    call MPI_BARRIER(m_comm_cart, ierr)
    call mpi_exch_r8type3du(m_comm_cart, r8type3du, nbid, ndim, u, symm_np)
    call mpi_exch_r8type3dv(m_comm_cart, r8type3dv, nbid, ndim, v, symm_np)
    if(symm_np) v(:,nyf+1,:) = -v(:,nyf+1,:)

    tp = max_vel
    call MPI_ALLREDUCE(tp, max_vel, 1, MPI_REAL8, MPI_MAX, m_comm_cart, ierr)

    tmp1 = sum(p0(1:nx,1:ny)**2)
    tp   = tmp1
    call MPI_ALLREDUCE(tp, tmp1, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

    tmp2 = sum(x(1:nx,1:ny)**2)
    tp   = tmp2
    call MPI_ALLREDUCE(tp, tmp2, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

    tmp2 = rho_sw*odt*dsqrt(tmp2/tmp1) 
    
!    if(myid_x .eq. 16 .and. myid_y .eq. 17) then
!      write(log_info,'(a,2i6,6f21.16)') &
!        "CHECK RESTART p: ", itf, itmask, w(nx/2,ny,1), u(nx/2,ny,1), u(nx/2+1,ny,1), v(nx/2,ny,1), v(nx/2,ny+1,1), x(nx/2,ny)
!      call comm_write_log_info(fid_log, log_info)
!    end if
    if(tmp2 .lt. threshold_p0 .and. itmask .gt. 2) exit
  end do
 
  if(myid .eq. rootid .and. mod(itf,1) .eq. 0) then
    if ( opt_solver == solver_BiCGStab ) then
      write(log_info,'(a,i8,a,i3,a,1e16.9,a14,f8.4,a4,a, i4, 2e10.3)') &
        '@itf-it0===', itf, ',itmask=', itmask, &
        ', p0 convergence rate=', tmp2, ', vmx on land=', max_vel, 'cm/s', &
        ', BICGSTAB convergence=', iter, res_p2, res_p3
    end if
    if ( opt_solver == solver_PCSI ) then  ! NTU 20250324 M. Hsieh
      write(log_info,*) &
        '@itf-it0===', itf, ',itmask=', itmask, &
        ', p0 convergence rate=', tmp2, ', vmx on land=', max_vel, 'cm/s', &
        ', PCSI iter, res = ', iter, res_p2
    end if
    call comm_write_log_info(fid_log, log_info)
  end if 
  u_change = iu(:,:,1)*(u(:,:,1) - u_change)
  v_change = iv(:,:,1)*(v(:,:,1) - v_change)
  p0(1:nx,1:ny) = p0(1:nx,1:ny)*in(1:nx,1:ny,1)
  ssh(1:nx,1:ny) = p0(1:nx,1:ny)/rho_sw/grav
end subroutine pressure_solver

subroutine check_incompressibility
  use hyperlink, only: u, v, w, odx, csv, ocs, ody, odz
  implicit none

  integer :: i, j, k, ierr
  real(r8) :: mxp, res, tep, tmp(3)
  character(len=256) :: log_info

  mxp = 0.d0
  res = 0.d0

  do j = 1, ny
    do i = 1, nx
      do k = kb(i,j),1,-1
        tmp(1) = (u(i+1,j,k)-u(i,j,k))*odx(j)
        tmp(2) = (csv(j+1)*v(i,j+1,k)-csv(j)*v(i,j,k))*ocs(j)*ody(j)
        tmp(3) = (w(i,j,k+1)-w(i,j,k))*odz(k)
        mxp = mxp + maxval(dabs(tmp))*in(i,j,k)
        res = res +    dabs(sum(tmp))*in(i,j,k)
      end do
    end do
  end do

  tep = res
  call MPI_ALLREDUCE(tep, res, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
  tep = mxp
  call MPI_ALLREDUCE(tep, mxp, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
  res = res/mxp
  incomp_res = res

  if(myid .eq. rootid) then
    write(log_info, '(a,1pe9.2)') ' *** normalized mean incompressibility error = ', res
    call comm_write_log_info(fid_log, log_info)
  end if

end subroutine check_incompressibility

subroutine eliminate_arbitrary_pressure
  use hyperlink, only: p0, ssh, nx_grid, ny_grid, opt_arbr_p0, area, total_area
  implicit none

  integer :: i,j, ierr
  real(r8) :: psm, tmp

  if(opt_arbr_p0 .eq. 1) then
    psm = sum(p0(1:nx,1:ny)*area(1:nx,1:ny))
    tmp = psm
    call MPI_ALLREDUCE(tmp, psm, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
    psm = psm/total_area
    p0(1:nx,1:ny) = (p0(1:nx,1:ny)-psm)*in(1:nx,1:ny,1)
    ssh(1:nx,1:ny) = p0(1:nx,1:ny)/rho_sw/grav
  end if

end subroutine eliminate_arbitrary_pressure

subroutine interp_incompressibility_to_Agrid
  use hyperlink, only: u2, v2, u, v, u_change, v_change, &
                       m_comm_cart, r8type3d, &
                       nbid, ndim, symm_np
  implicit none

  real(r8), dimension(nx,ny) :: incompress_change, high_order_change
  integer :: i, j, k, ierr

  incompress_change = 0.5d0*(u_change(1:nxf-1,1:ny)+u_change(2:nxf,1:ny))
  do k = 1, nz
    high_order_change = o24*iu(0:nxf-1,1:ny,k)*iu(3:nxf+1,1:ny,k) &
                      *(-u_change(0:nxf-1,1:ny) + u_change(1:nxf,  1:ny) &
                        +u_change(2:nxf+1,1:ny) - u_change(3:nxf+2,1:ny))

    u2(1:nx,1:ny,k) = u2(1:nx,1:ny,k) + incompress_change + high_order_change
  end do
  u2(1:nx,1:ny,:) = u2(1:nx,1:ny,:)*in(1:nx,1:ny,:)

  incompress_change = 0.5d0*(v_change(1:nx,1:nyf-1)+v_change(1:nx,2:nyf))
  do k = 1, nz
    high_order_change = o24*iv(1:nx,0:nyf-1,k)*iv(1:nx,3:nyf+1,k) &
                      *(-v_change(1:nx,0:nyf-1) + v_change(1:nx,2:nyf  ) &
                        +v_change(1:nx,2:nyf+1) - v_change(1:nx,3:nyf+2))

    v2(1:nx,1:ny,k) = v2(1:nx,1:ny,k) + incompress_change + high_order_change
  end do
  v2(1:nx,1:ny,:) = v2(1:nx,1:ny,:)*in(1:nx,1:ny,:)

  call MPI_BARRIER(m_comm_cart, ierr)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, u2, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, v2, symm_np)
  if(symm_np) v2(:,ny+1:ny+2,:) = -v2(:,ny+1:ny+2,:)

end subroutine interp_incompressibility_to_Agrid

subroutine interp_Cgrid_to_Agrid
  use hyperlink, only: u2, v2, u, v, &
                       m_comm_cart, r8type3d, &
                       nbid, ndim, symm_np
  implicit none

  real(r8) :: tmp, msk(0:nxf+1,0:nyf+1)
  integer :: i, j, k, ierr

  msk = 0.d0
  u2(1:nx,1:ny,1:nz) = 12.d0*(u(1:nx,1:ny,1:nz) + u(2:nxf,1:ny,1:nz))
  do k = 1, nz
    msk(0:nxf+1,1:ny) = dble(iu(0:nxf+1,1:ny,k))
    do j = 1, ny
      do i = 1, nx
        tmp = msk(i-1,j)*msk(i+2,j)
        u2(i,j,k) = u2(i,j,k) + tmp*(-u(i-1,j,k)+u(i,j,k)+u(i+1,j,k)-u(i+2,j,k))
      end do
    end do
  end do
  u2(1:nx,1:ny,1:nz) = o24*in(1:nx,1:ny,1:nz)*u2(1:nx,1:ny,1:nz)

  msk = 0.d0
  v2(1:nx,1:ny,1:nz) = 12.d0*(v(1:nx,1:ny,1:nz) + v(1:nx, 2:nyf,1:nz))
  do k = 1, nz
    msk(1:nx,0:nyf+1) = dble(iv(1:nx,0:nyf+1,k))
    do j = 1, ny
      do i = 1, nx
        tmp = msk(i,j-1)*msk(i,j+2)
        v2(i,j,k) = v2(i,j,k) + tmp*(-v(i,j-1,k)+v(i,j,k)+v(i,j+1,k)-v(i,j+2,k))
      end do
    end do
  end do
  v2(1:nx,1:ny,1:nz) = o24*in(1:nx,1:ny,1:nz)*v2(1:nx,1:ny,1:nz)

  call MPI_BARRIER(m_comm_cart, ierr)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, u2, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, v2, symm_np)
  if(symm_np) v2(:,ny+1:ny+2,:) = -v2(:,ny+1:ny+2,:)

end subroutine interp_Cgrid_to_Agrid

subroutine adjust_seaice_ts(itf, mxit, dt_cpl)
  use timcom_ice, only: ice_formation
  use hyperlink, only: in, odz, t2, s2, qice, aqice, qflux, ssh 
  implicit none

  integer :: itf, mxit
  real(r8) :: dt_cpl
  logical :: if_calu_qflux = .false.

  if(itf .eq. mxit) then
    if_calu_qflux = .true.
  else
    if_calu_qflux = .false.
  end if

  call ice_formation(nx, ny, nz, kb, odz, ssh, t2, s2, qice, aqice, qflux, if_calu_qflux, dt_cpl)
 
  t2(1:nx,1:ny,1:nz) = t2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)
  s2(1:nx,1:ny,1:nz) = s2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)

  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, t2, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, s2, symm_np)

!  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanT4 =', avg_glb_t)
!  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanS4 =', avg_glb_s)
end subroutine adjust_seaice_ts

subroutine modified_filter
  use hyperlink, only: u1, ulf, u2, &
                       v1, vlf, v2, &
                       t1, tlf, t2, &
                       s1, slf, s2
  implicit none
  
  call modified_filter_q(u1, ulf, u2)
  call modified_filter_q(v1, vlf, v2)
  call modified_filter_q(t1, tlf, t2)
  call modified_filter_q(s1, slf, s2)

  if(symm_np) then
    v1(:,ny+1:ny+2,:)  = -v1(:,ny+1:ny+2,:)
    vlf(:,ny+1:ny+2,:) = -vlf(:,ny+1:ny+2,:)
  end if
end subroutine modified_filter

subroutine modified_filter_q(q1, qlf, q2)
  use hyperlink, only: m_comm_cart, r8type3d, &
                       nbid, ndim, symm_np
  implicit none

  real(r8), intent(inout) ::   &
    q1 (-1:nx+2,-1:ny+2,1:nz), &
    qlf(-1:nx+2,-1:ny+2,1:nz)

  real(r8), intent(in) ::      &
    q2 (-1:nx+2,-1:ny+2,1:nz)

  real(r8), parameter :: &
    fltw = 0.10d0, &
    wraf = 0.53d0

  real(r8) :: tmp
  
  integer :: i, j, k

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        tmp = fltw*(0.5d0*q1(i,j,k) - qlf(i,j,k) + 0.5d0*q2(i,j,k))
        q1(i,j,k) = qlf(i,j,k) + wraf*tmp
        qlf(i,j,k) = q2(i,j,k) + (wraf-1.d0)*tmp
      end do
    end do
  end do

  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, q1,  symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, qlf, symm_np)

end subroutine modified_filter_q

subroutine check_blowup(itf)
  use hyperlink, only: v2, t2, s2
  use timcom_pncio, only: put_blowup_pncio
  use hyperlink, only: mytag, myid
  implicit none

  integer, intent(in) :: itf
  integer :: lmin(2), lmax(2), ierr
  real(r8) :: fmin, fmax, tmp
  character(len=256) :: outfile, log_info

  lmin = minloc(v2(1:nx,1:ny,1), mask=in(1:nx,1:ny,1) .eq. 1)
  fmin = minval(v2(1:nx,1:ny,1), mask=in(1:nx,1:ny,1) .eq. 1)
  if(fmin > 1.d20) fmin = 0.d0

  lmax = maxloc(v2(1:nx,1:ny,1), mask=in(1:nx,1:ny,1) .eq. 1)
  fmax = maxval(v2(1:nx,1:ny,1), mask=in(1:nx,1:ny,1) .eq. 1)
  if(fmax < -1.d20) fmax = 0.d0

  tmp = fmax
  call MPI_ALLREDUCE(tmp, fmax, 1, MPI_REAL8, MPI_MAX, m_comm_cart, ierr)
  tmp = fmin
  call MPI_ALLREDUCE(tmp, fmin, 1, MPI_REAL8, MPI_MIN, m_comm_cart, ierr)

  tmp = fmax - fmin

  if(tmp > 1.d3) then
    if(myid .eq. rootid) then
      write(log_info,*) "[ERROR]: un-expectedly large velocity at itf ", itf, fmax, fmin
      call comm_write_log_info(fid_log, log_info)
    end if

    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',itf,'_blowup.nc'

    call put_blowup_pncio(outfile)
    call comm_finalize
  end if

end subroutine check_blowup

subroutine check_output(itf)
  use timcom_pncio, only: put_flux_pncio, put_timcom_pncio
  use hyperlink, only: mytag, myid, daodt_out, daodt_out_unit
  implicit none

  integer, intent(in) :: itf
  character(len=256) :: outfile, log_info

  outfile = trim(path_output)
  !if(mod(itf,daodt_out) .eq.0) then
  !  write(outfile, '(a,a,i3.3,a,i6.6,a)') &
  !    trim(path_output),trim(mytag),member,'_',itf/daodt_out_unit,'_flux.nc'
  !  call put_flux_pncio(outfile)

    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_c_',itf,'.nc'
    call put_timcom_pncio(outfile)

    if(myid .eq. rootid) then
      write(log_info,*) "[INFO]: output standard file ", trim(outfile)
      call comm_write_log_info(fid_log, log_info)
    end if
  !end if

end subroutine check_output

subroutine check_output_flux(itf)
  use timcom_pncio, only: put_flux_pncio, put_timcom_pncio
  use hyperlink, only: mytag, myid, daodt_out, daodt_out_unit, daodt
  implicit none

  integer, intent(in) :: itf
  character(len=256) :: outfile, log_info
  real :: nitf,cal_ndt,rdt_out_unit,rdaodt
  integer :: ndt_out

  outfile = trim(path_output)
  if(mod(itf,daodt_out) .eq. 0) then
  nitf = real(itf)                 !integer to real
  rdaodt=real(daodt)
  rdt_out_unit=rdaodt/24.d0
  cal_ndt=nitf/rdt_out_unit  !cal real
  ndt_out=cal_ndt            !real to integer
#ifndef owdate
    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',ndt_out,'_flux.nc'
    call put_flux_pncio(outfile)

    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',ndt_out,'.nc'
    call put_timcom_pncio(outfile)
#else
    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',itf/daodt_out_unit,'_flux.nc'
    call put_flux_pncio(outfile)

    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',itf/daodt_out_unit,'.nc'
    call put_timcom_pncio(outfile)
#endif
    if(myid .eq. rootid) then
      write(log_info,*) "[INFO]: output file ", trim(outfile)
      call comm_write_log_info(fid_log, log_info)
    end if
  end if
end subroutine check_output_flux

#ifdef cpl
 #ifndef owdate
  subroutine check_output_all(itf)
  use calendar, only: gregd2_noleap
 #else
  subroutine check_output_all(itf,mjd_check)
  use calendar, only: gregd2_noleap, gregd2
 #endif
#else
  subroutine check_output_all(itf,mjd_check)
  use calendar, only: gregd2_noleap, gregd2 
#endif
  use timcom_pncio, only: put_flux_pncio, &
                          put_windmix_pncio, &
                          put_timcom_pncio
  use hyperlink, only: mytag, myid, daodt, daodt_out, daodt_out_unit
  implicit none

  integer, intent(in) :: itf

#ifdef cpl
  #ifdef owdate
  real(r8), intent(in) :: mjd_check
  #endif
#else
  real(r8), intent(in) :: mjd_check
#endif
  character(len=256) :: outfile, log_info
  integer :: greg(5), dtg
  real :: nitf,cal_ndt,rdt_out_unit,rdaodt
  integer :: ndt_out

#ifdef cpl
  if(mod(itf,daodt_out) .eq. 0 .or. itf .eq. 1) then

  #ifndef owdate
  nitf = real(itf)                 !integer to real
  rdaodt=real(daodt)
  rdt_out_unit=rdaodt/24.d0
  cal_ndt=nitf/rdt_out_unit  !cal real
  ndt_out=int(cal_ndt)            !real to integer

    !write(outfile, '(a,a,i3.3,a,i6.6,a)') &
    !  trim(path_output),trim(mytag),member,'_wndmix_',itf/daodt_out_unit,'.nc'
    !call put_windmix_pncio(outfile)

    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_srflx_',ndt_out,'.nc'
    call put_flux_pncio(outfile)

    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',ndt_out,'.nc'
    call put_timcom_pncio(outfile)

    if(myid .eq. rootid) then
      write(log_info,*) "[INFO]: output file ", trim(outfile)
      call comm_write_log_info(fid_log, log_info)
    end if
  #else
    call gregd2(mjd_check, greg(1), greg(2), greg(3), greg(4), greg(5))
    dtg = greg(1)*1000000 + greg(2)*10000 + greg(3)*100 + greg(4)
    write(outfile, '(a,a,i3.3,a,i10.10,a)') &
      trim(path_output),trim(mytag),member,'_',dtg,'.nc'
    call put_timcom_pncio(outfile)

    if(myid .eq. rootid) then
      write(log_info,*) "[INFO]: output file ", trim(outfile)
      call comm_write_log_info(fid_log, log_info)
    end if
  #endif

  end if

#else
  if(mod(itf,daodt_out) .eq. 0) then

    call gregd2_noleap(mjd_check, greg(1), greg(2), greg(3), greg(4), greg(5))
    dtg = greg(1)*1000000 + greg(2)*10000 + greg(3)*100 + greg(4)

    write(outfile, '(a,a,i3.3,a,i10,a)') &
      trim(path_output),trim(mytag),member,'_',dtg,'.nc'
    call put_timcom_pncio(outfile)

    if(myid .eq. rootid) then
      write(log_info,*) "[INFO]: output file ", trim(outfile)
      call comm_write_log_info(fid_log, log_info)
    end if
  end if
#endif
end subroutine check_output_all

subroutine check_restart(itf, mjd_check)
  use calendar,     only: gregd2_noleap, gregd2
  use timcom_pncio, only: put_restart_pncio
  use hyperlink, only: mytag, myid, daodt_restart
  implicit none

  integer, intent(in) :: itf
  real(r8), intent(in) :: mjd_check
  integer :: greg(5), dtg
  character(len=256) :: outfile, log_info

  if(mod(itf,daodt_restart) .eq. 0) then
#ifdef cpl
    call gregd2(mjd_check, greg(1), greg(2), greg(3), greg(4), greg(5))
#else
    call gregd2_noleap(mjd_check, greg(1), greg(2), greg(3), greg(4), greg(5)) 
#endif
    dtg = greg(1)*1000000 + greg(2)*10000 + greg(3)*100 + greg(4)
    write(outfile, '(a,a,i3.3,a,i10,a)') &
      trim(path_restart),trim(mytag),member,'_restart_',dtg,'.nc'
    call put_restart_pncio(outfile)

    if(myid .eq. rootid) then
      write(log_info,*) "[INFO]: output restart file ", trim(outfile), mjd_check
      call comm_write_log_info(fid_log, log_info)
    end if
  end if
end subroutine check_restart

subroutine go_check_conservation( var_chk, vol_chk, total_vol, msg_chk, t_avg_glb)
  implicit none

  real(r8), intent(inout) :: var_chk(nx,ny,nz)
  real(r8), intent(in)    :: vol_chk(nx,ny,nz), total_vol
  real(r8), intent(inout) :: t_avg_glb
  character(len=*), intent(in) :: msg_chk
  integer :: ierr
  real(r8) :: totalt, dtmp
  character(len=256) :: log_info

  totalt = sum(var_chk*vol_chk)/total_vol
  dtmp   = totalt
  call MPI_ALLREDUCE(dtmp, totalt, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
  if (myid .eq. rootid) then
    write(log_info,'(a40,2e18.3)') trim(msg_chk), totalt, totalt-t_avg_glb
    call comm_write_log_info(fid_log, log_info)
  end if
  t_avg_glb = totalt

end subroutine go_check_conservation
end module timcom_drv
