module tai_timcom_drv
  use tai_timcom_const
  use tai_timcom_comm
  use tai_timcom_grid
  use tai_timcom_general
  use tai_hyperlink, only: nx, ny, nz, nxf, nyf, nzf, &
                       xu_bgn, xu_end, yv_bgn, yv_end, &
                       m_comm_cart, myid, r8type3d, &
                       nbid, ndim, symm_np, &
                       kb, in, iu, iv, iw, iu0, iv0, &
                       area, volume, npx, npy, &
                       total_area, total_vol, &
                       avg_glb_t, avg_glb_s
  implicit none

contains

subroutine timcom_exec(itf_drv, timer)
  use tai_hyperlink, only: assign_grid_all,nx,ny,myid,x_grid,y_grid,u2,v2
  implicit none

  integer, intent(in) :: itf_drv
  type(timer_panel) :: timer
  real(r8) :: dtmp, mjd1, mjd2, time_stamp, time_spend

  integer :: itf, it0, mxit, n, ierr,i,j
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
    do itf = 1+it0, mxit
      mjd1 = dtmp*dble(itf-1)
      mjd2 = dtmp*dble(itf)  
      
      if(myid .eq. rootid) then
        write(log_info,'(a,f12.4,a3,f12.4,a,i8)') &
          "[INFO]: TIMCOM forecat from:", mjd1, " to", mjd2, ", steps:", itf
        call comm_write_log_info(fid_log, log_info)
      end if

      call get_lateral_bc(itf)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'pre timcom_fs check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'pre timcom_fs check v2=', v2(30,7,1)

      call timcom_fs(itf, it0, timer%syng_mon)
      call check_blowup(itf)
     

!       do i = 1, nx
!        do j =1, ny
!         if(x_grid(i) .gt. 128.30 .and. x_grid(i) .lt. 128.31  .and. y_grid(j) .gt. 31.05 .and.  y_grid(j) .lt. 31.06)then
!           write(*,*) 'check u2=', u2(i,j,1) , 'myid=', myid, i, j, itf
!         end if
!        end do
!       end do

!      if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'done check u2=', u2(30,7,1)
!      if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'done check v2=', v2(30,7,1)

      call check_output_all(itf)
!0719
!      call check_restart(itf, timer%base_mjd+mjd2)
      call check_restart(itf, timer%base_mjd+timer%init_mjd+mjd2)
      !if(itf .eq. 2881) call comm_finalize
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

subroutine timcom_fs(itf, it0, syng_mon)
  use tai_hyperlink, only: u2,v2, myid
  implicit none

  integer, intent(in) :: itf, it0, syng_mon
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 0 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 0 check v2=', v2(30,7,1)
  call update_boundary_condition
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 1 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 1 check v2=', v2(30,7,1)
  call solve_hydrostatic_equation
     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 2 check u2=', u2(30,7,1)
     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 2 check v2=', v2(30,7,1)
  call calu_surf_flux
     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 3 check u2=', u2(30,7,1)
     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 3 check v2=', v2(30,7,1)
  call wind_mixing
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 4 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 4 check v2=', v2(30,7,1)
  call solve_dynamic_equation
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 5 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 5 check v2=', v2(30,7,1)
  call impose_surface_flux
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 6 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 6 check v2=', v2(30,7,1)
  call impose_temp_salt_mixing
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 7 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 7 check v2=', v2(30,7,1)
  call impose_temp_salt_nudging(syng_mon)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 8 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 8 check v2=', v2(30,7,1)
  call impose_bottom_stress
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 9 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 9 check v2=', v2(30,7,1)
  call impose_trapezoidal_coriolis
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 10 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 10 check v2=', v2(30,7,1)
  call impose_open_boundary_conditons
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 11 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 11 check v2=', v2(30,7,1)
  call interp_Agrid_to_Cgrid
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 12 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 12 check v2=', v2(30,7,1)
  call pressure_solver(itf,it0)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 13 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 13 check v2=', v2(30,7,1)
  call solve_continuity_equation
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 14 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 14 check v2=', v2(30,7,1)
  call check_incompressibility
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 15 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 15 check v2=', v2(30,7,1)
  call eliminate_arbitrary_pressure
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 16 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 16 check v2=', v2(30,7,1)
  call interp_Cgrid_to_Agrid
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 17 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 17 check v2=', v2(30,7,1)
  call adjust_seaice_ts
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 18 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 18 check v2=', v2(30,7,1)
  call modified_filter
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 19 check u2=', u2(30,7,1)
!     if(myid .eq. 518 .and. itf .eq. 1) write(*,*) 'in timcom_fs 19 check v2=', v2(30,7,1)

end subroutine timcom_fs

subroutine get_lateral_bc(itf)
  use tai_timcom_pncio, only: get_lateral_bc_field, check_lateral_bc
  use tai_hyperlink, only: daodt_out_unit, myid_x, myid_y, &
                       u1, u2, v1, v2, t1, t2, s1, s2, ns_lbc, nc_lbc, &
                       jwbc, jebc, isbc, inbc, wbc, ebc, sbc, nbc
  implicit none

  integer, intent(in) :: itf
  integer :: fhr
  real(r8) :: frac
  character(len=256) :: field_file, log_info

  if(mod(itf, daodt_out_unit*6) .eq. 1) then
    fhr = (itf-1)/daodt_out_unit + 0
    frac = 1.d0/dble(daodt_out_unit*6)
    write(field_file,"(a,i6.6,a)") "./TIMCOM_glb/glb000_", fhr, ".nc"

    if(myid .eq. rootid) then
        write(log_info,'(a,a)') &
          "[INFO]: READ lateral bc file:", trim(field_file)
        call comm_write_log_info(fid_log, log_info)
      end if

    call get_lateral_bc_field(trim(field_file), ns_lbc, nc_lbc, &
                              jwbc, jebc, isbc, inbc, &
                               wbc,  ebc,  sbc,  nbc)

    write(field_file,"(a,i6.6,a)") "tai000_lbc_", fhr, ".nc"
    call check_lateral_bc(trim(field_file), wbc, ebc, sbc, nbc)

    if(myid_x .eq. 0) then
      wbc(:,:,1) = (wbc(:,:,1) - u2(0,1:ny,1:nz))*frac
      wbc(:,:,2) = (wbc(:,:,2) - v2(0,1:ny,1:nz))*frac
      wbc(:,:,3) = (wbc(:,:,3) - t2(0,1:ny,1:nz))*frac
      wbc(:,:,4) = (wbc(:,:,4) - s2(0,1:ny,1:nz))*frac

      u1(0,1:ny,1:nz) = u2(0,1:ny,1:nz)
      v1(0,1:ny,1:nz) = v2(0,1:ny,1:nz)
      t1(0,1:ny,1:nz) = t2(0,1:ny,1:nz)
      s1(0,1:ny,1:nz) = s2(0,1:ny,1:nz)
    end if
    if(myid_x .eq. npx-1) then
      ebc(:,:,1) = (ebc(:,:,1) - u2(nx+1,1:ny,1:nz))*frac
      ebc(:,:,2) = (ebc(:,:,2) - v2(nx+1,1:ny,1:nz))*frac
      ebc(:,:,3) = (ebc(:,:,3) - t2(nx+1,1:ny,1:nz))*frac
      ebc(:,:,4) = (ebc(:,:,4) - s2(nx+1,1:ny,1:nz))*frac
  
      u1(nx+1,1:ny,1:nz) = u2(nx+1,1:ny,1:nz)
      v1(nx+1,1:ny,1:nz) = v2(nx+1,1:ny,1:nz)
      t1(nx+1,1:ny,1:nz) = t2(nx+1,1:ny,1:nz)
      s1(nx+1,1:ny,1:nz) = s2(nx+1,1:ny,1:nz)
    end if
    if(myid_y .eq. 0) then
      sbc(:,:,1) = (sbc(:,:,1) - u2(1:nx,0,1:nz))*frac
      sbc(:,:,2) = (sbc(:,:,2) - v2(1:nx,0,1:nz))*frac
      sbc(:,:,3) = (sbc(:,:,3) - t2(1:nx,0,1:nz))*frac
      sbc(:,:,4) = (sbc(:,:,4) - s2(1:nx,0,1:nz))*frac

      u1(1:nx,0,1:nz) = u2(1:nx,0,1:nz)
      v1(1:nx,0,1:nz) = v2(1:nx,0,1:nz)
      t1(1:nx,0,1:nz) = t2(1:nx,0,1:nz)
      s1(1:nx,0,1:nz) = s2(1:nx,0,1:nz)
    end if
    if(myid_y .eq. npy-1) then
      nbc(:,:,1) = (nbc(:,:,1) - u2(1:nx,ny+1,1:nz))*frac
      nbc(:,:,2) = (nbc(:,:,2) - v2(1:nx,ny+1,1:nz))*frac
      nbc(:,:,3) = (nbc(:,:,3) - t2(1:nx,ny+1,1:nz))*frac
      nbc(:,:,4) = (nbc(:,:,4) - s2(1:nx,ny+1,1:nz))*frac

      u1(1:nx,ny+1,1:nz) = u2(1:nx,ny+1,1:nz)
      v1(1:nx,ny+1,1:nz) = v2(1:nx,ny+1,1:nz)
      t1(1:nx,ny+1,1:nz) = t2(1:nx,ny+1,1:nz)
      s1(1:nx,ny+1,1:nz) = s2(1:nx,ny+1,1:nz)
    end if
  end if
end subroutine get_lateral_bc

subroutine update_boundary_condition
  use tai_hyperlink,  only: t2, s2, u2, v2,         &
                        myid_x, myid_y, u, v,   &
                        wbc, ebc, sbc, nbc, u1, &
                        v1, t1, s1
  implicit none

  integer :: i, j, k

  if(myid_x .eq. 0 )  then
    do k = 1, nz
      do j = 1, ny
        u2(-1:0,j,k) = wbc(j,k,1) + u1(0,j,k)
        v2(-1:0,j,k) = wbc(j,k,2) + v1(0,j,k)
        t2(-1:0,j,k) = wbc(j,k,3) + t1(0,j,k)
        s2(-1:0,j,k) = wbc(j,k,4) + s1(0,j,k)
         u(0:1,j,k) = u2(0,j,k)
      end do
    end do
  end if

  if(myid_x .eq. npx-1) then
    do k = 1, nz
      do j = 1, ny
        u2(nxf:nxf+1,j,k) = ebc(j,k,1) + u1(nxf,j,k)
        v2(nxf:nxf+1,j,k) = ebc(j,k,2) + v1(nxf,j,k)
        t2(nxf:nxf+1,j,k) = ebc(j,k,3) + t1(nxf,j,k)
        s2(nxf:nxf+1,j,k) = ebc(j,k,4) + s1(nxf,j,k)
         u(nxf:nxf+1,j,k) = u2(nxf,j,k)
      end do
    end do
  end if

  if(myid_y .eq. 0 ) then
    do k = 1, nz
      do i = 1, nx
        u2(i,-1:0,k) = sbc(i,k,1) + u1(i,0,k)
        v2(i,-1:0,k) = sbc(i,k,2) + v1(i,0,k)
        t2(i,-1:0,k) = sbc(i,k,3) + t1(i,0,k)
        s2(i,-1:0,k) = sbc(i,k,4) + s1(i,0,k)
         v(i,0:1,k) = v2(i,0,k)
      end do
    end do
  end if

  if(myid_y .eq. npy-1) then
    do k = 1, nz
      do i = 1, nx
        u2(i,nyf:nyf+1,k) = nbc(i,k,1) + u1(i,nyf,k)
        v2(i,nyf:nyf+1,k) = nbc(i,k,2) + v1(i,nyf,k)
        t2(i,nyf:nyf+1,k) = nbc(i,k,3) + t1(i,nyf,k)
        s2(i,nyf:nyf+1,k) = nbc(i,k,4) + s1(i,nyf,k)
         v(i,nyf:nyf+1,k) = v2(i,nyf,k)
      end do
    end do
  end if
end subroutine update_boundary_condition

subroutine solve_hydrostatic_equation
  use tai_hyperlink,  only: p, p0, t2, s2, rho, &
                        z_grid, dzw, dz
  use tai_timcom_eos, only: compute_eos_mkcoef_zgrid
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
  use tai_hyperlink,  only: p, p0, t2, s2, rho, &
                        z_grid, dz 
  use tai_timcom_eos, only: compute_eos_wright, &
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
    if(res1(1) < 1.d-6) exit
  end do

  p = p*in
  call MPI_BARRIER(m_comm_cart, ierr)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, p , symm_np)

end subroutine solve_hydrostatic_equation_balance 

subroutine calu_surf_flux
  use tai_timcom_phy, only: surface_flux_atmocn
  use tai_hyperlink, only: u2, v2, t2, s2, rho, &
                       u_10, v_10, t_10, q_10, pslv, &
                       swup, swdn, lwup, lwdn, taux, &
                       tauy, senh, lath, evap, rain, &
                       snow, roff, ifrc, &
                       stf, shf_qsw, smft, myid
  implicit none

  integer :: nsize
  integer(2) :: mask(nx,ny)
  real(r8), dimension(nx,ny) :: &
    u_ocn, v_ocn, t_ocn, &
    qdot, wflux, qdot2
  integer :: i, j
  real(r8) :: tmp

  nsize = nx*ny

  mask  = in(1:nx,1:ny,1)
  u_ocn = u2(1:nx,1:ny,1)/100.d0
  v_ocn = v2(1:nx,1:ny,1)/100.d0
  t_ocn = t2(1:nx,1:ny,1) + 273.15d0

  call surface_flux_atmocn(nsize, mask, ifrc, & 
                      u_ocn, v_ocn, t_ocn, &
                      u_10, v_10, t_10, q_10, pslv,  &
                      taux, tauy, senh, lath, lwup, evap)


  qdot  = lwup + lwdn + senh + lath
  wflux = rain + snow + evap + roff
  qdot2 = swup + swdn

  do j = 1, ny
    do i = 1, nx
      stf(i,j,1)   = dble(in(i,j,1))*qdot(i,j)*hflux_factor
      stf(i,j,2)   = dble(in(i,j,1))*wflux(i,j)*(-34.7d0)*1.d-4
      shf_qsw(i,j) = dble(in(i,j,1))*qdot2(i,j)*hflux_factor
      smft(i,j,1)  = dble(in(i,j,1))*taux(i,j)
      smft(i,j,2)  = dble(in(i,j,1))*tauy(i,j)
    end do
  end do 

  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check u10=', u_10(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check v10=', v_10(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check t10=', t_10(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check q10=', q_10(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check pslv=', pslv(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check swdn=', swdn(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check swup=', swup(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check lwdn=', lwdn(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check taux=', taux(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check tauy=', tauy(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check lath=', lath(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check senh=', senh(30,7)
  if(myid .eq. 518 ) write(*,*) 'calu_surf_flux check evap=', evap(30,7)

end subroutine calu_surf_flux

subroutine wind_mixing
  use tai_timcom_kpp, only: kppglo
  use tai_timcom_phy, only: pp82
  use tai_hyperlink, only: w, u2, v2, t2, s2, rho, &
                       sw_trans, stf, shf_qsw, &
                       smft, vbk, hbk, add, vdc, &
                       vvc, ev, hv, kpp_src, kpp_hblt, &
                       opt_windmix
  implicit none

  integer, parameter :: &
    opt_kpp  = 1, &
    opt_pp82 = 2

  integer :: ierr

  select case(opt_windmix)
  case(opt_kpp) 
    call kppglo(u2=u2, v2=v2, t2=t2, s2=s2, &
                rho=rho, trans=sw_trans, stf=stf, shf_qsw=shf_qsw, smft=smft, &
                vbk=vbk, hbk=hbk, add=add, vdc=vdc, vvc=vvc, &
                ev=ev, hv=hv, kpp_src=kpp_src, kpp_hblt=kpp_hblt)
  
  case(opt_pp82)
    call pp82(nx, ny, nz, w, u2(1:nx,1:ny,1:nz), v2(1:nx,1:ny,1:nz), rho, hbk, ev, hv)
  
  end select

  call MPI_BARRIER(m_comm_cart, ierr)
end subroutine wind_mixing

subroutine solve_dynamic_equation
  use tai_hyperlink, only: u, v, w, p,  &
                       u1, ulf, u2, &
                       v1, vlf, v2, &
                       t1, tlf, t2, &
                       s1, slf, s2, &
                       dmx, dmy, ev, &
                       dhx, dhy, hv, myid        
  implicit none
  integer :: itf
  real(r8) :: p_grad(nx,ny,nz)
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-1 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-1 check v2=', v2(30,7,1), 'itf=', itf
  call pressure_gradient_x(p, p_grad) 
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2 check v2=', v2(30,7,1), 'itf=', itf  
  call solve_momentun_equation(u, v, w, p_grad, u1, ulf, u2, dmx, dmy, ev) 
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-3 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-3 check v2=', v2(30,7,1), 'itf=', itf

  call pressure_gradient_y(p, p_grad)
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-4 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-4 check v2=', v2(30,7,1), 'itf=', itf
  call solve_momentun_equation(u, v, w, p_grad, v1, vlf, v2, dmx, dmy, ev)
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-5 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-5 check v2=', v2(30,7,1), 'itf=', itf

!  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanTX =', avg_glb_t) 
!  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanSX =', avg_glb_s)
  
  call solve_tracer_equation(u, v, w, t1, tlf, t2, dhx, dhy, hv(:,:,:,1))
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-6 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-6 check v2=', v2(30,7,1), 'itf=', itf

  call solve_tracer_equation(u, v, w, s1, slf, s2, dhx, dhy, hv(:,:,:,2))
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-7 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-7 check v2=', v2(30,7,1), 'itf=', itf

!  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanT1 =', avg_glb_t)
!  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanS1 =', avg_glb_s)

 !call solve_tracer_equation(u, v, w, c1, clf, c2, dhx, dhy, hv(:,:,:,1))
end subroutine solve_dynamic_equation

subroutine pressure_gradient_x(p, px)
  use tai_hyperlink, only: odx
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
  use tai_hyperlink, only: ocs, ody
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
  use tai_hyperlink, only: dt,myid,u2,v2
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
  integer :: i, j, k, itf

  dqx_dx = 0.d0
  dqy_dy = 0.d0
  dqz_dz = 0.d0

  select case(lfsrf)
  case(0)
    qz(1:nx,1:ny) = 0.d0
  case(1)
    qz(1:nx,1:ny) = w(1:nx,1:ny,1)*q2(1:nx,1:ny,1)
  end select

!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-1 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-1 check v2=', v2(30,7,1), 'itf=', itf
  
  do k = 1, nz
    call compute_vertical_fluxes(k, w, q1, qlf, q2, vmix_z(:,:,k), qz, dqz_dz)
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-2 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-2 check v2=', v2(30,7,1), 'itf=', itf
    call compute_longitudinal_fluxes(k, u, q1(:,:,k), q2(:,:,k), hmix_x(:,:,k), dqx_dx)
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-3 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-3 check v2=', v2(30,7,1), 'itf=', itf
    call compute_latitudinal_fluxes (k, v, q1(:,:,k), q2(:,:,k), hmix_y(:,:,k), dqy_dy)
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-4 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-4 check v2=', v2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-4 check q1=', q1(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-4 check dtin=', dtin, 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-4 check dp=', dp(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-4 check dqx_dx=', dqx_dx(30,7), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-4 check dqy_dy=', dqy_dy(30,7), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-4 check dqz_dz=', dqz_dz(30,7), 'itf=', itf

    do j = 1, ny
      do i = 1, nx
        dtin = dt*dble(in(i,j,k))
        q2(i,j,k) = q1(i,j,k) - dtin*(dp(i,j,k) + dqx_dx(i,j) + dqy_dy(i,j) + dqz_dz(i,j)) 
      end do
    end do
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5 check u2=', u2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5 check v2=', v2(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5 check q1=', q1(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5 check dtin=', dtin, 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5 check dp=', dp(30,7,1), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5 check dqx_dx=', dqx_dx(30,7), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5 check dqy_dy=', dqy_dy(30,7), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5 check dqz_dz=', dqz_dz(30,7), 'itf=', itf

  end do

end subroutine solve_momentun_equation

subroutine solve_tracer_equation(u, v, w, q1, qlf, q2, hmix_x, hmix_y, vmix_z)
  use tai_hyperlink, only: dt
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
  use tai_hyperlink, only: odz,myid
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

  integer :: i, j, l, itf

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

!    if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-0 check qz_b=', qz_b(30,7), 'itf=', itf
!    if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-0 check w=', w(30,7,1), 'itf=', itf
!    if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-0 check scr=', scr(30,7), 'itf=', itf
!    if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-0 check vmix_z=', vmix_z(30,7), 'itf=', itf
!    if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-0 check q1=', q1(30,7,1), 'itf=', itf

    do j = 1, ny
      do i = 1, nx
        tmp = dble(iw(i,j,l))
        qz_b(i,j) = tmp*(w(i,j,l)*scr(i,j) - vmix_z(i,j)*(q1(i,j,l)-q1(i,j,k)))
      end do
    end do
  end if

!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-1 check dqz_dz=', dqz_dz(30,7), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-1 check dz_b=', qz_b(30,7), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-1 check dz_t=', qz_t(30,7), 'itf=', itf

  dqz_dz = (qz_b - qz_t)*odz(k)
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-2 check dqz_dz=', dqz_dz(30,7), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-2 check qz_b=', qz_b(30,7), 'itf=', itf
!     if(myid .eq. 518 ) write(*,*) 'in timcom_fs 4-2-5-2 check qz_t=', qz_t(30,7), 'itf=', itf
  qz_t  = qz_b

end subroutine compute_vertical_fluxes

subroutine compute_latitudinal_fluxes(k, v, q1, q2, hmix_y, dqy_dy)
  use tai_hyperlink, only: csv, odyv, ocs, ody
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
  use tai_hyperlink, only: odx
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
  use tai_timcom_phy, only: wind_surf, heat_surf, heat_volume, salt_surf
  use tai_hyperlink, only: dt, odz, u2, v2, t2, s2, &
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
  use tai_timcom_kpp, only: kppsrcglo
  use tai_hyperlink,  only: t2, s2, kpp_src, opt_windmix
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
  use tai_hyperlink, only: t2, s2, t_nudge, s_nudge, t_clim, s_clim
  implicit none

  integer, intent(in) :: syng_mon

  t2(1:nx,1:ny,1) = t2(1:nx,1:ny,1) + t_nudge(:,:)*(t_clim(:,:,syng_mon) - t2(1:nx,1:ny,1))
  s2(1:ny,1:ny,1) = s2(1:nx,1:ny,1) + s_nudge(:,:)*(s_clim(:,:,syng_mon) - s2(1:nx,1:ny,1))

!  t2(1:nx,1:ny,2:nz) = (1.d0-heman_nudge(1:nx,1:ny,2:nz))*t2(1:nx,1:ny,2:nz) &
!                     + hmean_nudge(1:nx,1:ny,2:nz)*t_clim(1:nx,1:ny,2:nz)
!  s2(1:nx,1:ny,2:nz) = (1.d0-heman_nudge(1:nx,1:ny,2:nz))*s2(1:nx,1:ny,2:nz) &
!                     + hmean_nudge(1:nx,1:ny,2:nz)*s_clim(1:nx,1:ny,2:nz)

end subroutine impose_temp_salt_nudging

subroutine impose_bottom_stress
  use tai_hyperlink, only: dt, odz, u1, u2, v1, v2
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
  use tai_hyperlink, only: dt, curv_f, tanphi, &
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
  use tai_hyperlink, only: u2, v2, t2, s2, u, v, &
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
  use tai_hyperlink, only: u, u2, v, v2, &
                       m_comm_cart, r8type3du, r8type3dv, &
                       nbid, ndim, symm_np
  implicit none

  integer :: i, j, k, ierr
  real(r8) :: tmp, tmp1
  real(r8) :: msk(0:nxf+1,0:nyf+1)

  msk = 0
  ! LONGITUDINAL DIRECTION
  u2(1:nx,1:ny,1:nz) = u2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)
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
  v2(1:nx,1:ny,1:nz) = v2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)
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
  call mpi_exch_r8type3du(m_comm_cart, r8type3du, nbid, ndim, u, .false.)
  call mpi_exch_r8type3dv(m_comm_cart, r8type3dv, nbid, ndim, v, .false.)
  !if(symm_np) v(:,nyf+1,:) = -v(:,nyf+1,:)

end subroutine interp_Agrid_to_Cgrid

subroutine solve_continuity_equation
! ======================================== !
! CALCULATE W BY CONTINUITY EQUATION       !
! ======================================== !
! this completes the advanced time level advection velocity
! having exactly zero barotropic divergence and satisfying all
! kinematic and specified open boundary conditions
! in hydrostatic assumption
  use tai_hyperlink, only: u, v, w, odx, ody, odz, ocs, csv
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
  use tai_timcom_solver, only: p_bicgstab_le
  use tai_hyperlink, only: p0, odz, u_change, v_change, &
                       ab, al, ac, ar, at, x, &
                       cb, cl, cc, cr, ct,    &
                       cgr, cgrh, cgp, cgv,   &
                       cgs, cgt, cgph, cgsh,  &
                       u, odx, v, odyv, w, odt,  &
                       m_comm_cart, symm_np,  & 
                       r8type3du, r8type3dv,  &
                       nbid, ndim, myid_x, myid_y, &
                       threshold_p0, max_iter_p0, daodt_out_unit

  implicit none

  integer, intent(in) :: itf, it0

  real(r8) :: s(nx,ny), scr, max_vel, tp, tmp1, tmp2, res_p2, res_p3
  integer :: i, j, k, ierr, itmask, mxmask, iter
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

    call p_bicgstab_le(ab, al, ac, ar, at, s, x, &
                       cb, cl, cc, cr, ct,       &
                       cgr, cgrh, cgp, cgv,      &
                       cgs, cgt, cgph, cgsh, iter, res_p2, res_p3, itf)

    p0 = p0 + odt*x*rho_sw

    max_vel = 0.d0
    do j = 1, ny
      do i = xu_bgn, xu_end
        scr = -(x(i,j)-x(i-1,j))*odx(j)
        u(i,j,:) = u(i,j,:) + scr*iu(i,j,:)
        max_vel = max(max_vel, dble(1-iu(i,j,1))*dabs(scr))
      end do
    end do

    do j = yv_bgn, yv_end
      do i = 1, nx
        scr = -(x(i,j)-x(i,j-1))*odyv(j)
        v(i,j,:) = v(i,j,:) + scr*iv(i,j,:)
        max_vel = max(max_vel, dble(1-iv(i,j,1))*dabs(scr))
      end do
    end do

    call MPI_BARRIER(m_comm_cart, ierr)
    call mpi_exch_r8type3du(m_comm_cart, r8type3du, nbid, ndim, u, .false.)
    call mpi_exch_r8type3dv(m_comm_cart, r8type3dv, nbid, ndim, v, .false.)
    !if(symm_np) v(:,nyf+1,:) = -v(:,nyf+1,:)
    
    tp = max_vel
    call MPI_ALLREDUCE(tp, max_vel, 1, MPI_REAL8, MPI_MAX, m_comm_cart, ierr)

    tmp1 = sum(p0(1:nx,1:ny)**2)
    tp   = tmp1
    call MPI_ALLREDUCE(tp, tmp1, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

    tmp2 = sum(x(1:nx,1:ny)**2)
    tp   = tmp2
    call MPI_ALLREDUCE(tp, tmp2, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

    tmp2 = rho_sw*odt*dsqrt(tmp2/tmp1) 
    
    if(tmp2 .lt. threshold_p0 .and. itmask .gt. 2) exit
  end do
 
  if(myid .eq. rootid .and. mod(itf,1) .eq.0) then
    write(log_info,'(a,i8,a,i8,a,1e16.9,a14,f8.4,a4,a, i4, 2e10.3)') &
      '@itf-it0===', itf, ',itmask=', itmask, &
      ', p0 convergence rate=', tmp2, ', vmx on land=', max_vel, 'cm/s', &
      ', BICGSTAB convergence=', iter, res_p2, res_p3
    call comm_write_log_info(fid_log, log_info)
  end if 
  u_change = iu(:,:,1)*(u(:,:,1) - u_change)
  v_change = iv(:,:,1)*(v(:,:,1) - v_change)
  p0(1:nx,1:ny) = p0(1:nx,1:ny)*in(1:nx,1:ny,1)
end subroutine pressure_solver

subroutine check_incompressibility
  use tai_hyperlink, only: u, v, w, odx, csv, ocs, ody, odz
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

  if(myid .eq. rootid) then
    write(log_info, '(a,1pe9.2)') ' *** normalized mean incompressibility error = ', res
    call comm_write_log_info(fid_log, log_info)
  end if

end subroutine check_incompressibility

subroutine eliminate_arbitrary_pressure
  use tai_hyperlink, only: p0, nx_grid, ny_grid, opt_arbr_p0, area, total_area
  implicit none

  integer :: i,j, ierr
  real(r8) :: psm, tmp

  if(opt_arbr_p0 .eq. 1) then
    psm = sum(p0(1:nx,1:ny)*area(1:nx,1:ny))
    tmp = psm
    call MPI_ALLREDUCE(tmp, psm, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
    psm = psm/total_area
    p0(1:nx,1:ny) = (p0(1:nx,1:ny)-psm)*in(1:nx,1:ny,1)
  end if

end subroutine eliminate_arbitrary_pressure

subroutine interp_incompressibility_to_Agrid
  use tai_hyperlink, only: u2, v2, u, v, u_change, v_change, &
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
  use tai_hyperlink, only: u2, v2, u, v, &
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

subroutine adjust_seaice_ts
  use tai_timcom_ice, only: ice_formation
  use tai_hyperlink, only: in, odz, t2, s2, p0, qice, aqice, qflux 
  implicit none

  call ice_formation(nx, ny, nz, kb, odz, p0, t2, s2, qice, aqice, qflux)
  
  t2(1:nx,1:ny,1:nz) = t2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)
  s2(1:nx,1:ny,1:nz) = s2(1:nx,1:ny,1:nz)*in(1:nx,1:ny,1:nz)

  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, t2, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, s2, symm_np)

!  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanT4 =', avg_glb_t)
!  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanS4 =', avg_glb_s)
end subroutine adjust_seaice_ts

subroutine modified_filter
  use tai_hyperlink, only: u1, ulf, u2, &
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
  use tai_hyperlink, only: m_comm_cart, r8type3d, &
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
  use tai_hyperlink, only: v2, t2, s2
  use tai_timcom_pncio, only: put_blowup_pncio
  use tai_hyperlink, only: mytag, myid
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
      write(log_info,*) "[ERROR]: unix-pectedly large velocity at itf ", itf, fmax, fmin
      call comm_write_log_info(fid_log, log_info)
    end if

    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',itf,'_blowup.nc'

    call put_blowup_pncio(outfile)
    call comm_finalize
  end if

end subroutine check_blowup

subroutine check_output(itf)
  use tai_timcom_pncio, only: put_flux_pncio, put_timcom_pncio
  use tai_hyperlink, only: mytag, myid, daodt_out, daodt_out_unit
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
  use tai_timcom_pncio, only: put_flux_pncio, put_timcom_pncio
  use tai_hyperlink, only: mytag, myid, daodt_out, daodt_out_unit
  implicit none

  integer, intent(in) :: itf
  character(len=256) :: outfile, log_info

  outfile = trim(path_output)
  if(mod(itf,daodt_out) .eq. 0) then
    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',itf/daodt_out_unit,'_flux.nc'
    call put_flux_pncio(outfile)
    
    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',itf/daodt_out_unit,'.nc'
    call put_timcom_pncio(outfile)
    
    if(myid .eq. rootid) then
      write(log_info,*) "[INFO]: output file ", trim(outfile)
      call comm_write_log_info(fid_log, log_info)
    end if
  end if
end subroutine check_output_flux

subroutine check_output_all(itf)
  use tai_timcom_pncio, only: put_flux_pncio, &
                          put_windmix_pncio, &
                          put_timcom_pncio
  use tai_hyperlink, only: mytag, myid, daodt, daodt_out, daodt_out_unit
  implicit none

  integer, intent(in) :: itf
  character(len=256) :: outfile, log_info

  if(mod(itf,daodt_out) .eq. 0 .or. itf .eq. 1) then
    !write(outfile, '(a,a,i3.3,a,i6.6,a)') &
    !  trim(path_output),trim(mytag),member,'_wndmix_',itf/daodt_out_unit,'.nc'
    !call put_windmix_pncio(outfile)

    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_srflx_',itf/daodt_out_unit,'.nc'
    call put_flux_pncio(outfile)

    write(outfile, '(a,a,i3.3,a,i6.6,a)') &
      trim(path_output),trim(mytag),member,'_',itf/daodt_out_unit,'.nc'
    call put_timcom_pncio(outfile)

    if(myid .eq. rootid) then
      write(log_info,*) "[INFO]: output file ", trim(outfile)
      call comm_write_log_info(fid_log, log_info)
    end if
  end if
end subroutine check_output_all

subroutine check_restart(itf, mjd_check)
  use tai_calendar,     only: gregd2
  use tai_timcom_pncio, only: put_restart_pncio
  use tai_hyperlink, only: mytag, myid, daodt_restart
  implicit none

  integer, intent(in) :: itf
  real(r8), intent(in) :: mjd_check
  integer :: greg(5), dtg
  character(len=256) :: outfile, log_info

  if(mod(itf,daodt_restart) .eq. 0) then
    call gregd2(mjd_check, greg(1), greg(2), greg(3), greg(4), greg(5)) 
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
end module tai_timcom_drv
