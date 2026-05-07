!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define NCCLCHECK(ierr) call ocn_nccl_check_helper(ierr, __FILE__, __LINE__)

module timcom_drv_gpu
  use timcom_const
  use timcom_comm
  use timcom_comm_gpu
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

  real(r8) :: incomp_res_gpu
  real(r8), allocatable :: s(:,:)

contains

  subroutine timcom_exec_gpu(itf_drv, timer)
    use hyperlink
    use timcom_kpp, only: zgrid_kpp => zgrid, tidal_coef_kpp => tidal_coef, bckgrnd_vvc_kpp => bckgrnd_vvc, &
                          bckgrnd_vdc_kpp => bckgrnd_vdc, tidal_diff_kpp => tidal_diff, &
                          fcort_kpp => fcort, bolus_sp_kpp => bolus_sp, fstokes_kpp => fstokes, &
                          hmxl_kpp => hmxl, hwide_kpp => hwide, cg_kpp => cg
    use timcom_drv
    implicit none

    integer, intent(in) :: itf_drv
    type(timer_panel) :: timer
    real(r8) :: dtmp, mjd1, mjd2, time_stamp, time_spend, dt_cpl

    integer :: itf, it0, mxit, n, ierr
    character(len=256) :: log_info
    integer :: async_id = 1

    if(myid .eq. rootid) then
      write(log_info,'(a,f12.4,a3,f12.4,a,i8)') &
        "[INFO]: DRIVER forecat from:", timer%mjd1, " to", timer%mjd2, ", steps:", itf_drv
      call comm_write_log_info(fid_log, log_info)
    end if

    do n = 1, ocn_grid_num
      call assign_grid_all(ocn_grid(n))

      time_stamp = MPI_WTIME()

      !$acc update async(async_id) device( &
      !$acc& u, v, p0, iw, odx, ody, odz, ocs, csv, kb, w, x, iu, iv, in, &
      !$acc& ab, al, ac, ar, at, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh, &
      !$acc& nx, ny, nxf, nyf, nz, max_iter_p0, odyv, u_change, v_change, &
      !$acc& xu_bgn, xu_end, yv_bgn, yv_end, odt, nzf, &
      !$acc& incomp_res_gpu, area, total_area, u2, v2, t2, s2, qice, aqice, qflux, &
      !$acc& u1, ulf, v1, vlf, t1, tlf, s1, slf, dt, npx, npy, symm_np, myid_x, myid_y, peri_x, peri_y, &
      !$acc& tanphi, curv_f, t_nudge, t_da, kpp_src, smft, stf, sw_trans, shf_qsw, &
      !$acc& p, dmx, dmy, ev, dhx, dhy, hv, rho, z_grid, dzw, &
      !$acc& taux, tauy, senh, lath, lwup, evap, ifrc, vice, vsno, melt, melth, salt, u_10, v_10, t_10, q_10, pslv, &
      !$acc& lwdn, rain, snow, roff, swup, swdn, vdc, vvc, hbk, vbk, odzw, &
      !$acc& zgrid_kpp, tidal_coef_kpp, bckgrnd_vvc_kpp, bckgrnd_vdc_kpp, tidal_diff_kpp, &
      !$acc& fcort_kpp, bolus_sp_kpp, fstokes_kpp, hmxl_kpp, kpp_hblt, hwide_kpp, cg_kpp, z_face)

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
        
        call timcom_fs_gpu(itf, it0, timer%syng_mon, mxit, dt_cpl)
        call check_blowup_gpu(itf)
  #ifdef cpl
   #ifndef owdate
        call check_output_all_gpu(itf)
   #else
        call check_output_all_gpu(itf, timer%base_mjd+mjd2)
   #endif
  #else
        call check_output_all_gpu(itf, timer%base_mjd+timer%init_mjd+mjd2)
!        call check_restart_gpu(itf, timer%base_mjd+timer%init_mjd+mjd2)
  #endif
        call check_restart_gpu(itf, timer%base_mjd+timer%init_mjd+mjd2)
        !if(itf .eq. 1800) call comm_finalize
      end do

      !$acc update async(async_id) self( &
      !$acc& u, v, p0, iw, odx, ody, odz, ocs, csv, kb, w, x, iu, iv, in, &
      !$acc& ab, al, ac, ar, at, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh, &
      !$acc& nx, ny, nxf, nyf, nz, max_iter_p0, odyv, u_change, v_change, &
      !$acc& xu_bgn, xu_end, yv_bgn, yv_end, odt, nzf, &
      !$acc& incomp_res_gpu, area, total_area, u2, v2, t2, s2, qice, aqice, qflux, &
      !$acc& u1, ulf, v1, vlf, t1, tlf, s1, slf, dt, npx, npy, symm_np, myid_x, myid_y, peri_x, peri_y, &
      !$acc& tanphi, curv_f, t_nudge, t_da, kpp_src, smft, stf, sw_trans, shf_qsw, &
      !$acc& p, dmx, dmy, ev, dhx, dhy, hv, rho, z_grid, dzw, &
      !$acc& taux, tauy, senh, lath, lwup, evap, ifrc, vice, vsno, melt, melth, salt, u_10, v_10, t_10, q_10, pslv, &
      !$acc& lwdn, rain, snow, roff, swup, swdn, vdc, vvc, hbk, vbk, odzw, &
      !$acc& zgrid_kpp, tidal_coef_kpp, bckgrnd_vvc_kpp, bckgrnd_vdc_kpp, tidal_diff_kpp, &
      !$acc& fcort_kpp, bolus_sp_kpp, fstokes_kpp, hmxl_kpp, kpp_hblt, hwide_kpp, cg_kpp, z_face)

      time_spend = MPI_WTIME() - time_stamp
      call MPI_ALLREDUCE(time_spend, time_stamp, 1, MPI_REAL8, MPI_MAX, m_comm_cart, ierr)
      if(myid .eq. rootid) then
        write(log_info,'(a,f8.1,a)') &
            "[INFO]: TIMCOM runtime spends ", time_stamp, " sec in a driver loop"
          call comm_write_log_info(fid_log, log_info)
      end if
    end do

  end subroutine timcom_exec_gpu

  subroutine timcom_fs_gpu(itf, it0, syng_mon, mxit, dt_cpl)
    use hyperlink, only: v2, t2, s2
    implicit none

    integer, intent(in) :: itf, it0, syng_mon, mxit
    real(r8), intent(in) :: dt_cpl
    character(len=256) :: log_info
    integer :: async_id

    async_id = 1

    if (myid .eq. rootid) then
      print *, "timcom_fs", itf, it0
    end if

    call solve_hydrostatic_equation_gpu(async_id)
    call calu_surf_flux_gpu(async_id)
  #ifdef cpl
   #ifndef clm_r
    call wind_mixing_gpu(async_id)
   #else
    call wind_mixing_gpu(itf, async_id)
   #endif
  #else
    call wind_mixing_gpu(itf, async_id)
  #endif
    call solve_dynamic_equation_gpu(async_id)
    call impose_surface_flux_gpu(async_id)
    call impose_temp_salt_mixing_gpu(async_id)
  #ifdef cpl
    call impose_temp_salt_nudging_gpu(syng_mon, async_id)
  #endif
    call impose_bottom_stress_gpu(async_id)
    call impose_trapezoidal_coriolis_gpu(async_id)
    call impose_open_boundary_conditons_gpu(async_id)
    call interp_Agrid_to_Cgrid_gpu(async_id)
    call pressure_solver_gpu(itf, it0, async_id)
    call solve_continuity_equation_gpu(async_id)
    call check_incompressibility_gpu(async_id)
    call eliminate_arbitrary_pressure_gpu(async_id)
    call interp_Cgrid_to_Agrid_gpu(async_id)
    call adjust_seaice_ts_gpu(itf,mxit,dt_cpl,async_id)
    call modified_filter_gpu(async_id)

    !$acc update async(async_id) self(v2, t2, s2)
    !$acc wait(async_id)

  end subroutine timcom_fs_gpu

  subroutine check_blowup_gpu(itf)
    use hyperlink, only: x_grid, y_grid, z_grid, u2, v2, t2, s2, u, v, w, p0
    use timcom_pncio, only: put_blowup_pncio
    use hyperlink, only: mytag, myid
    implicit none

    integer, intent(in) :: itf
    integer :: lmin(2), lmax(2), ierr
    real(r8) :: fmin, fmax, tmp
    character(len=256) :: outfile, log_info
    integer :: async_id

    async_id = 1

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

      !$acc update async(async_id) self(z_grid, u2, v2, t2, s2, u, v, w, p0)
      !$acc wait(async_id)
      call put_blowup_pncio(outfile)
      call comm_finalize
    end if

  end subroutine check_blowup_gpu

  #ifdef cpl
   #ifndef owdate
    subroutine check_output_all_gpu(itf)
   #else
    subroutine check_output_all_gpu(itf,mjd_check)
   #endif
  #else
   subroutine check_output_all_gpu(itf,mjd_check)
  #endif
    use timcom_drv, only: check_output_all
    use hyperlink, only: daodt_out, x_grid, y_grid, z_grid, u_10, &
                         v_10, t_10, q_10, pslv, taux, tauy, swup, swdn, &
                         lwup, lwdn, lath, senh, evap, rain, snow, ifrc, &
                         roff, u2, v2, t2, s2, rho, p0, w
    implicit none

    integer, intent(in) :: itf
  #ifdef cpl
    #ifdef owdate
    real(r8), intent(in) :: mjd_check
    #endif
  #else
    real(r8), intent(in) :: mjd_check
  #endif
   integer :: async_id

    async_id = 1

    #ifdef cpl
      if(mod(itf,daodt_out) .eq. 0 .or. itf .eq. 1) then
        !$acc update async(async_id) self(z_grid, u_10, v_10, t_10, q_10, &
        !$acc& pslv, taux, tauy, swup, swdn, lwup, lwdn, lath, senh, evap, rain, snow, ifrc, roff, &
        !$acc& u2, v2, t2, s2, rho, p0, w)
        !$acc wait(async_id)
      end if
    #else
      if(mod(itf,daodt_out) .eq. 0) then
        !$acc update async(async_id) self(z_grid, u2, v2, t2, s2, rho, p0, w)
        !$acc wait(async_id)
      end if
    #endif

    #ifdef cpl
     #ifndef owdate
      call check_output_all(itf)
     #else
      call check_output_all(itf, mjd_check)
     #endif
    #else
      call check_output_all(itf, mjd_check)
    #endif

  end subroutine check_output_all_gpu

  subroutine check_restart_gpu(itf, mjd_check)
    use calendar,     only: gregd2_noleap, gregd2
    use timcom_drv, only: check_restart
    use hyperlink, only: daodt_restart, x_grid, y_grid, z_grid, &
                         x_face, y_face, z_face, u1, u2, ulf, v1, v2, vlf, &
                         t1, t2, tlf, s1, s2, slf, u, v, w, p0, x, cgr, cgrh, &
                         cgp, cgv, cgs, cgt
    implicit none

    integer, intent(in) :: itf
    real(r8), intent(in) :: mjd_check
    integer :: async_id

    async_id = 1

    if(mod(itf,daodt_restart) .eq. 0) then
      !$acc update async(async_id) self(z_grid, &
      !$acc& z_face, u1, u2, ulf, v1, v2, vlf, t1, t2, tlf, &
      !$acc& s1, s2, slf, u, v, w, p0, x, cgr, cgrh, cgp, cgv, cgs, cgt)
      !$acc wait(async_id)
    end if

    call check_restart(itf, mjd_check)

  end subroutine check_restart_gpu

  subroutine solve_continuity_equation_gpu(async_id)
    use hyperlink, only: u, v, w, odx, ody, odz, ocs, csv
    implicit none

    integer, intent(in) :: async_id
    integer :: i, j, k
    real(r8) :: tmp, temp, sum_w

    select case(lfsrf)
    case(0)
      !$acc parallel loop collapse(2) private(tmp, temp, sum_w) async(async_id)
      do j = 1, ny
        do i = 1, nx
          temp = ocs(j)*ody(j)
          w(i,j,1) = 0.d0
          sum_w = 0.d0
          !$acc loop seq
          do k = 1, nz
            tmp = 1.d0/odz(k)
            sum_w = (sum_w &
                        - ((u(i+1,j,k)-u(i,j,k))*odx(j) &
                        + (csv(j+1)*v(i,j+1,k)-csv(j)*v(i,j,k))*temp)*tmp)
            w(i,j,k+1) = sum_w
          end do
          !$acc loop seq
          do k = nz + 2, nzf
            w(i,j,k) = 0.d0
          end do
        end do
      end do
    case(1)
      !$acc parallel loop collapse(2) private(temp, tmp, sum_w) async(async_id)
      do j = 1, ny
        do i = 1, nx
          temp = ocs(j)*ody(j)
          sum_w = 0.d0
          !$acc loop seq
          do k = kb(i,j), 1, -1
            tmp = 1.d0/odz(k)
            sum_w = sum_w &
                      + ((u(i+1,j,k)-u(i,j,k))*odx(j) &
                      + (csv(j+1)*v(i,j+1,k)-csv(j)*v(i,j,k))*temp)*tmp
            w(i,j,k) = sum_w*iw(i,j,k)
          end do
          !$acc loop seq
          do k = kb(i,j) + 1, nzf
            w(i,j,k) = 0.d0
          end do
        end do
      end do
    end select
  end subroutine solve_continuity_equation_gpu


  subroutine impose_bottom_stress_gpu(async_id)
    use hyperlink, only: dt, odz, u1, u2, v1, v2, kb, in
    implicit none

    integer, intent(in) :: async_id
    integer :: i, j, k
    real(r8) :: tmp, drg
    real(r8), parameter :: drag = 1.d-3

    tmp = dt*drag

    !$acc parallel loop collapse(2) private(i, j , k, drg) async(async_id)
    do j = 1, ny
      do i = 1, nx
        k = kb(i,j)
        drg = in(i,j,k)*tmp*odz(k)*dsqrt(u1(i,j,k)**2+v1(i,j,k)**2)
        drg = min(0.2d0, drg)
        u2(i,j,k) = u2(i,j,k) - drg*u1(i,j,k)
        v2(i,j,k) = v2(i,j,k) - drg*v1(i,j,k)
      end do
    end do
  end subroutine impose_bottom_stress_gpu

  subroutine impose_trapezoidal_coriolis_gpu(async_id)
  use hyperlink, only: dt, curv_f, tanphi, &
                       u1, ulf, u2, v1, v2
  implicit none

  integer, intent(in) :: async_id
  integer i, j, k
  real(r8) :: tmp, temp, qu, qv

  !$acc parallel loop collapse(3) private(i, j, k, tmp, temp, qu, qv) async(async_id)
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
  end subroutine impose_trapezoidal_coriolis_gpu

  subroutine impose_open_boundary_conditons_gpu(async_id)
    use hyperlink, only: u2, v2, t2, s2, u, v, &
                        ulf, vlf, tlf, slf, dt, &
                        odx, ody, ocs, csv, &
                        myid_x, myid_y, peri_x, peri_y
    implicit none

    integer, intent(in) :: async_id
    integer :: i, j, k, ierr
    real(r8) :: tmp, tmpin, temp, temp1, temp2

    !$acc parallel loop collapse(3) &
    !$acc& private(i, j, k, tmp, tmpin, temp, temp1, temp2) async(async_id)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          u2(i,j,k) = u2(i,j,k)*in(i,j,k)
          v2(i,j,k) = v2(i,j,k)*in(i,j,k)
          t2(i,j,k) = t2(i,j,k)*in(i,j,k)
          s2(i,j,k) = s2(i,j,k)*in(i,j,k)
          if (.not. peri_x .and. myid_x .eq. 0 .and. i .eq. 1) then
            tmp = 0.5d0*dt
            tmpin = in(1,j,k)*tmp*odx(j)
            temp  = dabs(u(1,j,k))
            temp1 = tmpin*(temp + u(1,j,k))
            temp2 = tmpin*(temp - u(1,j,k))
            u2(1,j,k) = u2(1,j,k) + temp1*u2(0,j,k) - temp2*ulf(1,j,k)
            v2(1,j,k) = v2(1,j,k) + temp1*v2(0,j,k) - temp2*vlf(1,j,k)
            t2(1,j,k) = t2(1,j,k) + temp1*t2(0,j,k) - temp2*tlf(1,j,k)
            s2(1,j,k) = s2(1,j,k) + temp1*s2(0,j,k) - temp2*slf(1,j,k)
          end if
          if (.not. peri_x .and. myid_x .eq. npx-1 .and. i .eq. nx) then
            tmp = 0.5d0*dt
            tmpin = in(nx,j,k)*tmp*odx(j)
            temp  = dabs(u(nxf,j,k))
            temp1 = tmpin*(temp + u(nxf,j,k))
            temp2 = tmpin*(temp - u(nxf,j,k))
            u2(nx,j,k) = u2(nx,j,k) - temp1*ulf(nx,j,k) + temp2*u2(nx+1,j,k)
            v2(nx,j,k) = v2(nx,j,k) - temp1*vlf(nx,j,k) + temp2*v2(nx+1,j,k)
            t2(nx,j,k) = t2(nx,j,k) - temp1*tlf(nx,j,k) + temp2*t2(nx+1,j,k)
            s2(nx,j,k) = s2(nx,j,k) - temp1*slf(nx,j,k) + temp2*s2(nx+1,j,k)
          end if
          if (.not. peri_y .and. myid_y .eq. 0 .and. j .eq. 1) then
            tmp = 0.5d0*dt*ody(1)*csv(1)*ocs(1)
            tmpin = in(i,1,k)*tmp
            temp  = dabs(v(i,1,k))
            temp1 = tmpin*(temp + v(i,1,k))
            temp2 = tmpin*(temp - v(i,1,k))
            u2(i,1,k) = u2(i,1,k) + temp1*u2(i,0,k) - temp2*ulf(i,1,k)
            v2(i,1,k) = v2(i,1,k) + temp1*v2(i,0,k) - temp2*vlf(i,1,k)
            t2(i,1,k) = t2(i,1,k) + temp1*t2(i,0,k) - temp2*tlf(i,1,k)
            s2(i,1,k) = s2(i,1,k) + temp1*s2(i,0,k) - temp2*slf(i,1,k)
          end if
          if(.not. peri_y .and. myid_y .eq. npy-1 &
             .and. .not. symm_np .and. j .eq. ny) then
            tmp = 0.5d0*dt*ody(ny)*csv(nyf)*ocs(ny)
            tmpin = in(i,ny,k)*tmp
            temp  = dabs(v(i,nyf,k))
            temp1 = tmpin*(temp + v(i,nyf,k))
            temp2 = tmpin*(temp - v(i,nyf,k))
            u2(i,ny,k) = u2(i,ny,k) - temp1*ulf(i,ny,k) + temp2*u2(i,ny+1,k)
            v2(i,ny,k) = v2(i,ny,k) - temp1*vlf(i,ny,k) + temp2*v2(i,ny+1,k)
            s2(i,ny,k) = s2(i,ny,k) - temp1*slf(i,ny,k) + temp2*s2(i,ny+1,k)
            t2(i,ny,k) = t2(i,ny,k) - temp1*tlf(i,ny,k) + temp2*t2(i,ny+1,k)
          end if
        end do
      end do
    end do

    call MPI_BARRIER(m_comm_cart, ierr)

    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, u2, symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, v2, symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, t2, symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, s2, symm_np, async_id)

    if (symm_np) then
      !$acc parallel loop collapse(3) private(i, j, k) async(async_id)
      do k = 1, nz
        do j = ny+1, ny+2
          do i = -1, nx+2
            v2(i,j,k) = -v2(i,j,k)
          end do
        end do
      end do
    end if

  end subroutine impose_open_boundary_conditons_gpu

  subroutine impose_surface_flux_gpu(async_id)
    use timcom_phy, only: wind_surf, heat_surf, heat_volume, salt_surf
    use hyperlink, only: dt, odz, u2, v2, t2, s2, &
                        smft, stf, sw_trans, shf_qsw, &
                        nx, ny, nz
    implicit none

    integer, intent(in) :: async_id
    integer :: i, j, k

    !$acc parallel loop collapse(3) async(async_id)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (k .eq. 1) then
            ! A.wind stress
            u2(i,j,1) = u2(i,j,1) + dt*odz(1)*smft(i,j,1)
            v2(i,j,1) = v2(i,j,1) + dt*odz(1)*smft(i,j,2)
            ! B.heat fluxes, qdot
            ! (longwave + latent + sensible), qdot (positive downward)
            t2(i,j,1) = t2(i,j,1) + dt*odz(1)*stf(i,j,1)
            ! D.evaporation and precipitation (no salt sources)
            s2(i,j,1) = s2(i,j,1) + dt*odz(1)*1000.d0*stf(i,j,2)
          end if
          ! C. short wave radiation, qdot2 (positive downward)
          t2(i,j,k) = t2(i,j,k) + dt*odz(k)*shf_qsw(i,j)*sw_trans(i,j,k)
        end do
      end do
    end do

  end subroutine impose_surface_flux_gpu

  subroutine pressure_solver_gpu(itf,it0, async_id)
    use timcom_solver, only: solver_BiCGStab, solver_PCSI
    use timcom_solver_gpu, only: p_bicgstab_le_gpu, pcsi_gpu, res2, res3
    use timcom_comm_gpu, only: mpi_exch_r8type3du_gpu, mpi_exch_r8type3dv_gpu, nccl_comm
    use timcom_drv
    use hyperlink, only: p0, odz, u_change, v_change, &
                         ab, al, ac, ar, at, x, &
                         cb, cl, cc, cr, ct,    &
                         cgr, cgrh, cgp, cgv,   &
                         cgs, cgt, cgph, cgsh,  &
                         u, odx, v, odyv, w, odt,  &
                         m_comm_cart, symm_np,  &
                         r8type3du, r8type3dv,  &
                         nbid, ndim, myid_x, myid_y, &
                         threshold_p0, max_iter_p0, daodt_out_unit, &
                         opt_solver
    use hyperlink, only: ody, ocs, csv
    use openacc
    use cudafor
    use nccl

    implicit none

    integer, intent(in) :: itf, it0

    real(r8) :: scr, max_vel, tp, tmp1, tmp2, res_p2, res_p3, max_t, sum_t, sum_t2
    integer :: i, j, k, ierr, itmask, iter
    character(len=256) :: log_info
    integer :: async_id
    integer(kind=cuda_stream_kind) :: stream

    stream = acc_get_cuda_stream(async_id)

    !$acc enter data create(max_vel, tmp1, tmp2) async(async_id)

    !$acc parallel loop collapse(2) async(async_id)
    do j = 0, ny + 1
      do i = 0, nxf + 1
        u_change(i, j) = u(i, j, 1)
      end do
    end do

    !$acc parallel loop collapse(2) async(async_id)
    do j = 0, nyf + 1
      do i = 0, nx + 1
        v_change(i, j) = v(i, j, 1)
      end do
    end do

    do itmask = 1, max_iter_p0

      call solve_continuity_equation_gpu(async_id)

      select case(lfsrf)
      case(0)
        !$acc parallel loop collapse(2) async(async_id)
        do j = 1, ny
          do i = 1, nx
            s(i,j) = w(i,j,kb(i,j)+1)
          end do
        end do
      case(1)
        !$acc parallel loop collapse(2) async(async_id)
        do j = 1, ny
          do i = 1, nx
            s(i, j) = -w(i, j, 1)
          end do
        end do
      end select

      select case(opt_solver)
      case(solver_BiCGStab)
        call p_bicgstab_le_gpu(ab, al, ac, ar, at, s, x, &
                               cb, cl, cc, cr, ct,       &
                               cgr, cgrh, cgp, cgv,      &
                               cgs, cgt, cgph, cgsh, iter, async_id)
      case(solver_PCSI)
        call pcsi_gpu(ab, al, ac, ar, at, s, x, iter, res_p2, async_id)
      end select

      !$acc parallel loop collapse(2) async(async_id)
      do j = 1, ny
        do i = 1, nx
          p0(i, j) = p0(i, j) + odt * x(i, j) * rho_sw
        end do
      end do

      !$acc kernels present(max_vel) async(async_id)
      max_vel = 0.d0
      !$acc end kernels

      !$acc parallel loop gang async(async_id) private(max_t)
      do j = 1, ny
        max_t = 0.d0
        !$acc loop vector reduction(max:max_t) private(scr)
        do i = xu_bgn, xu_end
          scr = -(x(i,j)-x(i-1,j))*odx(j)
          !$acc loop seq
          do k = 1, nz
            u(i,j,k) = u(i,j,k) + scr*iu(i,j,k)
          end do
          max_t = max(max_t, dble(1-iu(i,j,1))*dabs(scr))
        end do
        !$acc atomic
        max_vel = max(max_vel, max_t)
      end do

      !$acc parallel loop gang async(async_id) private(max_t)
      do j = yv_bgn, yv_end
        max_t = 0.d0
        !$acc loop vector reduction(max:max_t) private(scr)
        do i = 1, nx
          scr = -(x(i,j)-x(i,j-1))*odyv(j)
          !$acc loop seq
          do k = 1, nz
            v(i,j,k) = v(i,j,k) + scr*iv(i,j,k)
          end do
          max_t = max(max_t, dble(1-iv(i,j,1))*dabs(scr))
        end do
        !$acc atomic
        max_vel = max(max_vel, max_t)
      end do

      call mpi_exch_r8type3du_gpu(m_comm_cart, r8type3du, nbid, ndim, u, symm_np, async_id)
      call mpi_exch_r8type3dv_gpu(m_comm_cart, r8type3dv, nbid, ndim, v, symm_np, async_id)

      if (symm_np) then
        !$acc parallel loop collapse(2) async(async_id)
        do k = 1, nz
          do i = 0, nx + 1
            v(i, nyf + 1, k) = -v(i, nyf + 1, k)
          end do
        end do
      end if

      !$acc host_data use_device(max_vel)
      NCCLCHECK(ncclAllReduce(max_vel, max_vel, 1, ncclFloat64, ncclMax, nccl_comm, stream))
      !$acc end host_data

      !$acc kernels present(tmp1, tmp2) async(async_id)
      tmp1 = 0.d0
      tmp2 = 0.d0
      !$acc end kernels

      !$acc parallel loop gang private(sum_t, sum_t2) async(async_id)
      do j = 1, ny
        sum_t = 0.d0
        sum_t2 = 0.d0
        !$acc loop vector reduction(+:sum_t) reduction(+:sum_t2)
        do i = 1, nx
          sum_t = sum_t + p0(i, j) * p0(i, j)
          sum_t2 = sum_t2 + x(i, j) * x(i, j)
        end do
        !$acc atomic
        tmp1 = tmp1 + sum_t
        !$acc atomic
        tmp2 = tmp2 + sum_t2
      end do

      !$acc host_data use_device(tmp1, tmp2)
      NCCLCHECK(ncclAllReduce(tmp1, tmp1, 1, ncclFloat64, ncclSum, nccl_comm, stream))
      NCCLCHECK(ncclAllReduce(tmp2, tmp2, 1, ncclFloat64, ncclSum, nccl_comm, stream))
      !$acc end host_data

      !$acc kernels present(tmp1, tmp2) async(async_id)
      tmp2 = rho_sw*odt*dsqrt(tmp2/tmp1)
      !$acc end kernels

      !    if(myid_x .eq. 16 .and. myid_y .eq. 17) then
      !      write(log_info,'(a,2i6,6f21.16)') &
      !        "CHECK RESTART p: ", itf, itmask, w(nx/2,ny,1), u(nx/2,ny,1), u(nx/2+1,ny,1), v(nx/2,ny,1), v(nx/2,ny+1,1),
      !        x(nx/2,ny)
      !      call comm_write_log_info(fid_log, log_info)
      !    end if

      !$acc update self(tmp2) async(async_id)
      !$acc wait(async_id)
      if(tmp2 .lt. threshold_p0 .and. itmask .gt. 2) exit
    end do

    if(myid .eq. rootid .and. mod(itf,1) .eq. 0) then
      !$acc update self(max_vel) async(async_id)
      !$acc wait(async_id)
      select case(opt_solver)
      case(solver_BiCGStab)
        write(log_info,'(a,i8,a,i3,a,1e16.9,a14,f8.4,a4,a, i4, 2e10.3)') &
          '@itf-it0===', itf, ',itmask=', itmask, &
          ', p0 convergence rate=', tmp2, ', vmx on land=', max_vel, 'cm/s', &
          ', BICGSTAB convergence=', iter, res_p2, res_p3
      case(solver_PCSI)
        write(log_info,*) &
          '@itf-it0===', itf, ',itmask=', itmask, &
          ', p0 convergence rate=', tmp2, ', vmx on land=', max_vel, 'cm/s', &
          ', PCSI iter, res = ', iter, res_p2
      end select
      call comm_write_log_info(fid_log, log_info)
    end if

    !$acc parallel loop collapse(2) async(async_id)
    do j = 0, ny + 1
      do i = 0, nxf + 1
        u_change(i, j) = iu(i, j, 1)*(u(i, j, 1) - u_change(i, j))
      end do
    end do

    !$acc parallel loop collapse(2) async(async_id)
    do j = 0, nyf + 1
      do i = 0, nx + 1
        v_change(i, j) = iv(i, j, 1)*(v(i, j, 1) - v_change(i, j))
      end do
    end do

    !$acc parallel loop collapse(2) async(async_id)
    do j = 1, ny
      do i = 1, nx
        p0(i, j) = p0(i, j) * in(i, j, 1)
      end do
    end do

    !$acc exit data delete(max_vel, tmp1, tmp2) async(async_id)

  end subroutine pressure_solver_gpu

  subroutine impose_temp_salt_nudging_gpu(syng_mon, async_id)
    use hyperlink, only: t2, s2, t_nudge, s_nudge, t_clim, s_clim, t_da, opt_da, nx, ny
    implicit none

    integer, intent(in) :: syng_mon
    integer, intent(in) :: async_id
    integer :: i, j

    if(opt_da .eq. 1) then
      !$acc parallel loop collapse(2) async(async_id)
      do j = 1, ny
        do i = 1, nx
          t2(i,j,1) = t2(i,j,1) + t_nudge(i,j)*(t_da(i,j) - t2(i,j,1))
        end do
      end do
    end if

  end subroutine impose_temp_salt_nudging_gpu

  subroutine interp_Agrid_to_Cgrid_gpu(async_id)
    use timcom_comm_gpu, only: mpi_exch_r8type3du_gpu, mpi_exch_r8type3dv_gpu
    use hyperlink, only: u, u2, v, v2, &
                         m_comm_cart, r8type3du, r8type3dv, &
                         nbid, ndim, symm_np
    implicit none

    integer, intent(in) :: async_id
    integer :: i, j, k, ierr
    real(r8) :: tmp, tmp1

    ! LONGITUDINAL DIRECTION

    !$acc parallel loop collapse(3) private(i, j, k, tmp, tmp1) async(async_id)
    do k = 1, nz
      do j = 1, yv_end
        do i = 1, xu_end
          if (i >= xu_bgn .and. j <= ny) then
            u(i,j,k) = 6.0d0*(u2(i-1,j,k) + u2(i,j,k))
            tmp1 = 0.0d0
            if (i-1 >= 0 .and. i+1 <= nxf+1) then
              tmp1 = iu(i-1,j,k)*iu(i+1,j,k)
            end if
            tmp  = (-u2(i-2,j,k) + u2(i-1,j,k) + u2(i,j,k) - u2(i+1,j,k))*tmp1
            u(i,j,k) = (u(i,j,k) + tmp) * iu(i,j,k) * o12
          end if

          if (j >= yv_bgn .and. i <= nx) then
            v(i,j,k) = 6.d0*(v2(i,j-1,k) + v2(i,j,k))
            tmp1 = 0.0d0
            if (j-1 >= 0 .and. j+1 <= nyf+1) then
              tmp1 = iv(i,j-1,k)*iv(i,j+1,k)
            end if
            tmp  = (-v2(i,j-2,k) + v2(i,j-1,k) + v2(i,j,k) - v2(i,j+1,k))*tmp1
            v(i,j,k) = (v(i,j,k) + tmp) * iv(i,j,k) * o12
          end if
        end do
      end do
    end do

    call mpi_exch_r8type3du_gpu(m_comm_cart, r8type3du, nbid, ndim, u, symm_np, async_id)
    call mpi_exch_r8type3dv_gpu(m_comm_cart, r8type3dv, nbid, ndim, v, symm_np, async_id)

    if(symm_np) then
      !$acc parallel loop collapse(3) private(i, j, k) async(async_id)
      do k = 1, nz
        do j = nyf+1, nyf+2
          do i = 0, nx+1
            v(i,j,k) = -v(i,j,k)
          end do
        end do
      end do
    end if

  end subroutine interp_Agrid_to_Cgrid_gpu

  subroutine interp_Cgrid_to_Agrid_gpu(async_id)
    use hyperlink, only: u2, v2, u, v, &
                        m_comm_cart, r8type3d, &
                        nbid, ndim, symm_np
    implicit none

    integer, intent(in) :: async_id
    real(r8) :: tmp
    integer :: i, j, k, ierr

    !$acc parallel loop collapse(3) private(i, j, k, tmp) async(async_id)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          tmp = 0.0d0
          if (i-1 >= 0 .and. i+2 <= nxf+1) then
            tmp = dble(iu(i-1,j,k))*dble(iu(i+2,j,k))
          end if
          u2(i,j,k) = 12.d0*(u(i,j,k) + u(i+1,j,k))
          u2(i,j,k) = u2(i,j,k) + tmp*(-u(i-1,j,k)+u(i,j,k)+u(i+1,j,k)-u(i+2,j,k))
          u2(i,j,k) = o24 * in(i,j,k) * u2(i,j,k)

          v2(i,j,k) = 12.d0*(v(i,j,k) + v(i,j+1,k))
          tmp = 0.0d0
          if (j-1 >= 0 .and. j+2 <= nyf+1) then
            tmp = dble(iv(i,j-1,k))*dble(iv(i,j+2,k))
          end if
          v2(i,j,k) = v2(i,j,k) + tmp*(-v(i,j-1,k)+v(i,j,k)+v(i,j+1,k)-v(i,j+2,k))
          v2(i,j,k) = o24 * in(i,j,k) * v2(i,j,k)
        end do
      end do
    end do

    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, u2, symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, v2, symm_np, async_id)

    if (symm_np) then
      !$acc parallel loop collapse(3) private(i, j, k) async(async_id)
      do k = 1, nz
        do j = ny+1, ny+2
          do i = -1, nx+2
            v2(i,j,k) = -v2(i,j,k)
          end do
        end do
      end do
    end if

  end subroutine interp_Cgrid_to_Agrid_gpu

  subroutine adjust_seaice_ts_gpu(itf, mxit, dt_cpl, async_id)
    use timcom_ice_gpu, only: ice_formation_gpu
    use hyperlink, only: in, odz, t2, s2, ssh, qice, aqice, qflux, &
                         m_comm_cart, r8type3d, &
                         nbid, ndim, symm_np, &
                         nx, ny, nz, kb
    implicit none

    integer, intent(in) :: itf, mxit
    real(r8), intent(in) :: dt_cpl
    integer, intent(in) :: async_id
    logical :: if_calu_qflux = .false.
    integer :: i, j, k

    if(itf .eq. mxit) then
      if_calu_qflux = .true.
    else
      if_calu_qflux = .false.
    end if
    call ice_formation_gpu(nx, ny, nz, kb, odz, ssh, t2, s2, qice, aqice, qflux, if_calu_qflux, dt_cpl, async_id)

    !$acc parallel loop collapse(3) async(async_id)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          t2(i,j,k) = t2(i,j,k)*in(i,j,k)
          s2(i,j,k) = s2(i,j,k)*in(i,j,k)
        end do
      end do
    end do

    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, t2, symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, s2, symm_np, async_id)

  end subroutine adjust_seaice_ts_gpu

  subroutine modified_filter_gpu(async_id)
    use hyperlink, only: u1, ulf, u2, &
                       v1, vlf, v2, &
                       t1, tlf, t2, &
                       s1, slf, s2
    use timcom_drv, only: modified_filter_q
    implicit none

    integer, intent(in) :: async_id
    integer :: i, j, k
    real(r8) :: tmp

    real(r8), parameter :: &
      fltw = 0.10d0, &
      wraf = 0.53d0

    !$acc parallel loop collapse(3) async(async_id) private(i, j, k, tmp)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          tmp = fltw*(0.5d0*u1(i,j,k) - ulf(i,j,k) + 0.5d0*u2(i,j,k))
          u1(i,j,k) = ulf(i,j,k) + wraf*tmp
          ulf(i,j,k) = u2(i,j,k) + (wraf-1.d0)*tmp
          tmp = fltw*(0.5d0*v1(i,j,k) - vlf(i,j,k) + 0.5d0*v2(i,j,k))
          v1(i,j,k) = vlf(i,j,k) + wraf*tmp
          vlf(i,j,k) = v2(i,j,k) + (wraf-1.d0)*tmp
          tmp = fltw*(0.5d0*t1(i,j,k) - tlf(i,j,k) + 0.5d0*t2(i,j,k))
          t1(i,j,k) = tlf(i,j,k) + wraf*tmp
          tlf(i,j,k) = t2(i,j,k) + (wraf-1.d0)*tmp
          tmp = fltw*(0.5d0*s1(i,j,k) - slf(i,j,k) + 0.5d0*s2(i,j,k))
          s1(i,j,k) = slf(i,j,k) + wraf*tmp
          slf(i,j,k) = s2(i,j,k) + (wraf-1.d0)*tmp
        end do
      end do
    end do

    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, u1,  symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, ulf, symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, v1,  symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, vlf, symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, t1,  symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, tlf, symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, s1,  symm_np, async_id)
    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, slf, symm_np, async_id)

    if(symm_np) then
      !$acc parallel loop collapse(2) private(i, k) async(async_id)
      do k = 1, nz
        do i = -1, nx+2
          v1(i,ny+1:ny+2,k)  = -v1(i,ny+1:ny+2,k)
          vlf(i,ny+1:ny+2,k) = -vlf(i,ny+1:ny+2,k)
        end do
      end do
    end if

  end subroutine modified_filter_gpu

  subroutine eliminate_arbitrary_pressure_gpu(async_id)
    use hyperlink, only: p0, nx_grid, ny_grid, opt_arbr_p0, area, total_area
    use timcom_comm_gpu
    use openacc
    use nccl
    implicit none

    integer, intent(in) :: async_id
    integer :: i, j, ierr
    integer(kind=cuda_stream_kind) :: stream
    real(r8) :: psm

    if(opt_arbr_p0 .eq. 1) then
      !$acc enter data async(async_id) create(psm)

      !$acc kernels present(psm) async(async_id)
      psm = 0.0_r8
      !$acc end kernels

      !$acc parallel loop collapse(2) reduction(+:psm) present(psm) async(async_id)
      do j = 1, ny
        do i = 1, nx
          psm = psm + p0(i,j)*area(i,j)
        end do
      end do

      stream = acc_get_cuda_stream(async_id)
      !$acc host_data use_device(psm)
      NCCLCHECK(ncclAllReduce(psm, psm, 1, ncclDouble, ncclSum, nccl_comm, stream))
      !$acc end host_data

      !$acc parallel loop collapse(2) private(i, j) present(psm) async(async_id)
      do j = 1, ny
        do i = 1, nx
          p0(i,j) = (p0(i,j)-psm/total_area)*in(i,j,1)
        end do
      end do
      !$acc exit data async(async_id) delete(psm)
    end if

  end subroutine eliminate_arbitrary_pressure_gpu

  subroutine impose_temp_salt_mixing_gpu(async_id)
    use hyperlink,  only: t2, s2, kpp_src, opt_windmix, dt, in, nx, ny, nz
    implicit none

    integer, intent(in) :: async_id
    integer, parameter :: &
      opt_kpp  = 1, &
      opt_pp82 = 2

    integer :: ierr, i, j, k

    select case(opt_windmix)
    case(opt_kpp)
      !$acc parallel loop collapse(3) async(async_id)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            t2(i,j,k) = t2(i,j,k) + dt*kpp_src(i,j,k,1)*in(i,j,k)
            s2(i,j,k) = s2(i,j,k) + dt*kpp_src(i,j,k,2)*in(i,j,k)
          end do
        end do
      end do

    case(opt_pp82)
      continue
    end select

  end subroutine impose_temp_salt_mixing_gpu

  subroutine check_incompressibility_gpu(async_id)
    use nccl
    use openacc
    use hyperlink, only: u, v, w, odx, csv, ocs, ody, odz, in, nx, ny, nz, kb
    implicit none

    integer, intent(in) :: async_id
    integer :: i, j, k, ierr
    real(r8) :: mxp, res, tep, tmp(3)
    character(len=256) :: log_info
    integer(kind=cuda_stream_kind) :: stream

    !$acc enter data async(async_id) create(mxp, res, tmp)

    !$acc kernels present(mxp, res) async(async_id)
    mxp = 0.d0
    res = 0.d0
    !$acc end kernels

    !$acc parallel loop collapse(3) async(async_id) &
    !$acc& private(i, j, k, tmp) &
    !$acc& reduction(+:mxp, res) present(mxp, res)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (k <= kb(i,j)) then
            tmp(1) = (u(i+1,j,k)-u(i,j,k))*odx(j)
            tmp(2) = (csv(j+1)*v(i,j+1,k)-csv(j)*v(i,j,k))*ocs(j)*ody(j)
            tmp(3) = (w(i,j,k+1)-w(i,j,k))*odz(k)
            mxp = mxp + max(dabs(tmp(1)), max(dabs(tmp(2)), dabs(tmp(3))))*in(i,j,k)
            res = res + dabs(tmp(1) + tmp(2) + tmp(3))*in(i,j,k)
          end if
        end do
      end do
    end do

    stream = acc_get_cuda_stream(async_id)
    !$acc host_data use_device(mxp, res)
    NCCLCHECK(ncclAllReduce(mxp, mxp, 1, ncclDouble, ncclSum, nccl_comm, stream))
    NCCLCHECK(ncclAllReduce(res, res, 1, ncclDouble, ncclSum, nccl_comm, stream))
    !$acc end host_data

    !$acc kernels present(res, mxp, incomp_res_gpu) async(async_id)
    res = res/mxp
    incomp_res_gpu = res
    !$acc end kernels
    !$acc update self(res) async(async_id)
    !$acc wait(async_id)

    if(myid .eq. rootid) then
      write(log_info, '(a,1pe9.2)') ' *** normalized mean incompressibility error = ', res
      call comm_write_log_info(fid_log, log_info)
    end if

    !$acc exit data async(async_id) delete(mxp, res, tmp)

  end subroutine check_incompressibility_gpu

  subroutine calu_surf_flux_gpu(async_id)
    use timcom_phy, only: surface_flux_atmocn, opt_bulk, opt_bulk_ncar, opt_bulk_cesm, &
                        shr_const_zvir, shr_const_cpdair, shr_const_cpvir, &
                        shr_const_g, shr_const_karman, shr_const_stebol, &
                        shr_const_latvap
    use hyperlink, only: u2, v2, t2, s2, rho, &
                        u_10, v_10, t_10, q_10, pslv, &
                        swup, swdn, lwup, lwdn, taux, &
                        tauy, senh, lath, evap, rain, &
                        snow, roff, ifrc, melth, melt, salt,&
                        stf, shf_qsw, smft
    implicit none

    integer, intent(in) :: async_id
    integer(2) :: mask(nx,ny)
    real(r8), dimension(nx,ny) :: u_ocn, v_ocn, t_ocn
    real(r8) :: qdot, wflux, qdot2
    integer :: i, j
    real(r8) :: tmp

    ! Variables needed for inlined surface_flux_atmocn_gpu
    real(r8) :: rhovec_a, zatm, work

    ! Variables needed for inlined cesm_shr_flux_atmocn_gpu
    real(r8), parameter :: &
      umin  =  0.5d0, & ! minimum wind speed       (m/s)
      zref  = 10.0d0, & ! reference height           (m)
      ztref =  2.0d0    ! reference height for air T (m)

    real(r8) :: &
      vmag,     & ! surface wind magnitude   (m/s)
      thvbot,   & ! virtual temperature      (K)
      ssq,      & ! sea surface humidity     (kg/kg)
      delt,     & ! potential T difference   (K)
      delq,     & ! humidity difference      (kg/kg)
      stable,   & ! stability factor
      rdn,      & ! sqrt of neutral exchange coeff (momentum)
      rhn,      & ! sqrt of neutral exchange coeff (heat)
      ren,      & ! sqrt of neutral exchange coeff (water)
      rd,       & ! sqrt of exchange coefficient (momentum)
      rh,       & ! sqrt of exchange coefficient (heat)
      re,       & ! sqrt of exchange coefficient (water)
      ustar,    & ! ustar
      qstar,    & ! qstar
      tstar,    & ! tstar
      hol,      & ! H (at zbot) over L
      xsq,      & ! ?
      xqq,      & ! ?
      psimh,    & ! stability function at zbot (momentum)
      psixh,    & ! stability function at zbot (heat and water)
      psix2,    & ! stability function at ztref reference height
      alz,      & ! ln(zbot/zref)
      al2,      & ! ln(zref/ztref)
      u10n,     & ! 10m neutral wind
      tau,      & ! stress at zbot
      cp,       & ! specific heat of moist air
      bn,       & ! exchange coef funct for interpolation
      bh,       & ! exchange coef funct for interpolation
      fac,      & ! vertical interpolation factor
      spval       ! local missing value

    !$acc enter data async(async_id) create(mask, u_ocn, v_ocn, t_ocn)

    !$acc parallel loop collapse(2) private(qdot, wflux, qdot2) async(async_id)
    do j = 1, ny
      do i = 1, nx
        mask(i, j)  = in(i,j,1)
        u_ocn(i, j) = u2(i,j,1)/100.d0
        v_ocn(i, j) = v2(i,j,1)/100.d0
        t_ocn(i, j) = t2(i,j,1) + 273.15d0

        select case(opt_bulk)
        case(opt_bulk_ncar)
          continue
        case(opt_bulk_cesm)
          zatm = 10.d0
          rhovec_a = 0.d0
          if (mask(i, j) .eq. 1) then
            rhovec_a = pslv(i, j)/((1.d0+0.608d0*q_10(i, j))*287.04d0*t_10(i, j))
          end if

          al2 = dlog(zref/ztref)
          taux(i, j) = 0.d0
          tauy(i, j) = 0.d0
          senh(i, j) = 0.d0
          lath(i, j) = 0.d0
          lwup(i, j) = 0.d0
          evap(i, j) = 0.d0

          if(mask(i, j) /= 0) then
            !--- compute some needed quantities ---
            vmag   = max(umin, dsqrt((u_10(i, j)-u_ocn(i, j))**2 + (v_10(i, j)-v_ocn(i, j))**2))
            thvbot = t_10(i, j) * (1.d0 + shr_const_zvir*q_10(i, j)) ! virtual temp (K)
            ssq    = 0.98d0*(640380.d0 / dexp(5107.4d0/t_ocn(i, j)))/rhovec_a ! sea surf hum (kg/kg)
            delt   = t_10(i, j) - t_ocn(i, j)                         ! pot temp diff (K)
            delq   = q_10(i, j) - ssq                               ! spec hum dif (kg/kg)
            alz    = dlog(zatm/zref)
            cp     = shr_const_cpdair*(1.d0 + shr_const_cpvir*ssq)

            !------------------------------------------------------------
            ! first estimate of Z/L and ustar, tstar and qstar
            !------------------------------------------------------------

            !--- neutral coefficients, z/L = 0.0 ---
            stable = 0.5d0 + dsign(0.5d0, delt)
            rdn    = dsqrt(0.0027d0 / vmag + 0.000142d0 + 0.0000764d0 * vmag)
            rhn    = (1.d0-stable)*0.0327d0 + stable*0.018d0
            ren    = 0.0346d0

            !--- ustar, tstar, qstar ---
            ustar = rdn*vmag
            tstar = rhn*delt
            qstar = ren*delq

            !--- compute stability & evaluate all stability functions ---
            hol    = shr_const_karman*shr_const_g*zatm*  &
                      (tstar/thvbot+qstar/(1.d0/shr_const_zvir+q_10(i, j)))/ustar**2
            hol    = dsign( min(dabs(hol),10.d0), hol )
            stable = 0.5d0 + dsign(0.5d0, hol)
            xsq    = max(dsqrt(dabs(1.d0 - 16.d0*hol)), 1.d0)
            xqq    = dsqrt(xsq)
            psimh  = -5.d0*hol*stable + (1.d0-stable)*(dlog((1.d0+xqq*(2.d0+xqq))*(1.d0+xqq*xqq)/8.d0) - 2.d0*datan(xqq) + 1.571d0)
            psixh  = -5.d0*hol*stable + (1.d0-stable)*(2.d0*dlog((1.d0 + xqq*xqq)/2.d0))

            !--- shift wind speed using old coefficient ---
            rd   = rdn / (1.d0 + rdn/shr_const_karman*(alz-psimh))
            u10n = vmag * rd / rdn

            !--- update transfer coeffs at 10m and neutral stability ---
            rdn = dsqrt(0.0027d0 / u10n + 0.000142d0 + 0.0000764d0 * u10n)
            ren = 0.0346d0
            rhn = (1.d0 - stable)*0.0327d0 + stable*0.018d0

            !--- shift all coeffs to measurement height and stability ---
            rd = rdn/(1.d0 + rdn/shr_const_karman*(alz-psimh))
            rh = rhn/(1.d0 + rhn/shr_const_karman*(alz-psixh))
            re = ren/(1.d0 + ren/shr_const_karman*(alz-psixh))

            !--- update ustar, tstar, qstar using updated, shifted coeffs ---
            ustar = rd * vmag
            tstar = rh * delt
            qstar = re * delq

            !------------------------------------------------------------
            ! iterate to converge on Z/L, ustar, tstar and qstar
            !------------------------------------------------------------
            !--- compute stability & evaluate all stability functions ---
            hol    = shr_const_karman*shr_const_g*zatm* &
                      (tstar/thvbot+qstar/(1.d0/shr_const_zvir+q_10(i, j)))/ustar**2
            hol    = dsign( min(dabs(hol),10.d0), hol )
            stable = .5d0 + dsign(.5d0 , hol)
            xsq    = max(dsqrt(dabs(1.d0 - 16.d0*hol)) , 1.d0)
            xqq    = dsqrt(xsq)
            psimh  = -5.d0*hol*stable + (1.d0-stable)*(dlog((1.d0+xqq*(2.d0+xqq))*(1.d0+xqq*xqq)/8.d0) - 2.d0*datan(xqq) + 1.571d0)
            psixh  = -5.d0*hol*stable + (1.d0-stable)*(2.d0*dlog((1.d0 + xqq*xqq)/2.d0))

            !--- shift wind speed using old coeffs ---
            rd   = rdn / (1.d0 + rdn/shr_const_karman*(alz-psimh))
            u10n = vmag*rd/rdn

            !--- update transfer coeffs at 10m and neutral stability ---
            rdn = dsqrt(0.0027d0 / u10n + 0.000142d0 + 0.0000764d0 * u10n)
            ren = 0.0346d0
            rhn = (1.d0 - stable)*0.0327d0 + stable*0.018d0

            !--- shift all coeffs to measurement height and stability ---
            rd = rdn/(1.d0 + rdn/shr_const_karman*(alz-psimh))
            rh = rhn/(1.d0 + rhn/shr_const_karman*(alz-psixh))
            re = ren/(1.d0 + ren/shr_const_karman*(alz-psixh))

            !--- update ustar, tstar, qstar using updated, shifted coeffs ---
            ustar = rd*vmag
            tstar = rh*delt
            qstar = re*delq

            !------------------------------------------------------------
            ! compute the fluxes
            !------------------------------------------------------------
            tau = rhovec_a*ustar*ustar
            !--- momentum flux ---
            taux(i, j) = tau*(u_10(i, j)-u_ocn(i, j))/vmag
            ! tauy
            tauy(i, j) = tau*(v_10(i, j)-v_ocn(i, j))/vmag
            !--- heat flux ---
            ! sensible heat
            senh(i, j) = cp*tau*tstar/ustar
            ! latent heat
            lath(i, j) = shr_const_latvap*tau*qstar/ustar
            ! long wave radiation
            lwup(i, j) = -shr_const_stebol*t_ocn(i, j)**4
            !--- water flux ---
            ! evaporation
            evap(i, j) = lath(i, j)/shr_const_latvap
          end if

          work = mask(i, j)*rhovec_a*1.63d-3*dsqrt(u_10(i, j)**2 + v_10(i, j)**2)
          taux(i, j) = 10.d0*(ifrc(i, j)*work*u_10(i, j) + (1.d0-ifrc(i, j))*taux(i, j))
          tauy(i, j) = 10.d0*(ifrc(i, j)*work*v_10(i, j) + (1.d0-ifrc(i, j))*tauy(i, j))
        end select

        qdot  = (1.d0-ifrc(i,j))*(lwup(i,j) + lwdn(i,j) + senh(i,j) + lath(i,j) - snow(i,j)*3.337d5) + melth(i,j)
        wflux = (1.d0-ifrc(i,j))*(rain(i,j) + snow(i,j) + evap(i,j) + roff(i,j)) + melt(i,j)
        qdot2 = (1.d0-ifrc(i,j))*(swup(i,j) + swdn(i,j))
        stf(i,j,1)   = dble(in(i,j,1))*qdot*hflux_factor
        stf(i,j,2)   = dble(in(i,j,1))*(wflux*(-34.7d0)*1.d-4 + salt(i,j)*0.1d0)
        shf_qsw(i,j) = dble(in(i,j,1))*qdot2*hflux_factor
        smft(i,j,1)  = dble(in(i,j,1))*taux(i,j)
        smft(i,j,2)  = dble(in(i,j,1))*tauy(i,j)
      end do
    end do

    !$acc exit data async(async_id) delete(mask, u_ocn, v_ocn, t_ocn)

  end subroutine calu_surf_flux_gpu

  subroutine solve_hydrostatic_equation_gpu(async_id)
    !$acc routine(eos_mkcoef) seq
    use timcom_eos_gpu, only: eos_mkcoef
    use hyperlink,  only: p, p0, t2, s2, rho, &
                          z_grid, dzw, dz
    implicit none

    integer, intent(in) :: async_id

    real(r8) :: tmpd, wface(nx,ny,nzf), sltd, pisd
    integer :: i, j, k, n, ierr

    !$acc parallel loop collapse(3) async(async_id)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho(i,j,k) = 0.d0
          if(in(i,j,k) .eq. 1) then
            pisd = z_grid(k)*1.d-3
            tmpd = max(-1.8d0, t2(i,j,k))
            tmpd = min(40.d0, tmpd)
            sltd = min(41.d0, s2(i,j,k))
            sltd = max(0.d0, sltd)
            rho(i,j,k) = eos_mkcoef(tmpd, sltd, pisd)
          end if
          if (k == 1) then
            tmpd = grav*z_grid(1)
            p(i,j,k) = p0(i,j) + tmpd*rho(i,j,k)
          end if
        end do
      end do
    end do

    !$acc parallel loop collapse(2) async(async_id)
    do j = 1, ny
      do i = 1, nx
        !$acc loop seq
        do k = 2, nz
          tmpd = 0.5d0*grav*dzw(k)
          p(i,j,k) = p(i,j,k-1) + tmpd*(rho(i,j,k)+rho(i,j,k-1))
        end do
      end do
    end do

    !$acc parallel loop collapse(3) async(async_id)
    do k = 1, nz
      do j = -1, ny+2
        do i = -1, nx+2
          p(i,j,k) = p(i,j,k)*in(i,j,k)
        end do
      end do
    end do

    call ghost_cell_exch_gpu(m_comm_cart, r8type3d, nbid, ndim, p, symm_np, async_id)

  end subroutine solve_hydrostatic_equation_gpu

#ifdef cpl
 #ifndef clm_r
  subroutine wind_mixing_gpu(async_id)
 #else
  subroutine wind_mixing_gpu(itf, async_id)
 #endif
#else
  subroutine wind_mixing_gpu(itf, async_id)
#endif
    use timcom_kpp_gpu, only: kppglo_gpu
    use timcom_phy_gpu, only: pp82_gpu
    use hyperlink, only: w, u2, v2, t2, s2, rho, &
                        sw_trans, stf, shf_qsw, &
                        smft, vbk, hbk, add, vdc, &
                        vvc, ev, hv, kpp_src, kpp_hblt, &
                        opt_windmix, odzw, dt
    implicit none

    integer, intent(in) :: async_id
    integer, parameter :: &
      opt_kpp  = 1, &
      opt_pp82 = 2

    integer :: ierr, itf

    select case(opt_windmix)
    case(opt_kpp)
#ifdef cpl
 #ifndef clm_r
      call kppglo_gpu(u2=u2, v2=v2, t2=t2, s2=s2, &
                  rho=rho, trans=sw_trans, stf=stf, shf_qsw=shf_qsw, smft=smft, &
                  vbk=vbk, hbk=hbk, add=add, vdc=vdc, vvc=vvc, &
                  ev=ev, hv=hv, kpp_src=kpp_src, kpp_hblt=kpp_hblt, async_id=async_id)
 #else
      call kppglo_gpu(itf, u2=u2, v2=v2, t2=t2, s2=s2, &
                  rho=rho, trans=sw_trans, stf=stf, shf_qsw=shf_qsw, smft=smft, &
                  vbk=vbk, hbk=hbk, add=add, vdc=vdc, vvc=vvc, &
                  ev=ev, hv=hv, kpp_src=kpp_src, kpp_hblt=kpp_hblt, async_id=async_id)
 #endif
#else
      call kppglo_gpu(itf, u2=u2, v2=v2, t2=t2, s2=s2, &
                  rho=rho, trans=sw_trans, stf=stf, shf_qsw=shf_qsw, smft=smft, &
                  vbk=vbk, hbk=hbk, add=add, vdc=vdc, vvc=vvc, &
                  ev=ev, hv=hv, kpp_src=kpp_src, kpp_hblt=kpp_hblt, async_id=async_id)
#endif
    case(opt_pp82)
      call pp82_gpu(nx, ny, nz, w, u2, v2, rho, hbk, ev, hv, async_id)

    end select

    call MPI_BARRIER(m_comm_cart, ierr)
  end subroutine wind_mixing_gpu

  subroutine solve_dynamic_equation_gpu(async_id)
    use hyperlink, only: u, v, w, p,  &
                        u1, ulf, u2, &
                        v1, vlf, v2, &
                        t1, tlf, t2, &
                        s1, slf, s2, &
                        dmx, dmy, ev, &
                        dhx, dhy, hv, &
                        nx, ny, nz
    implicit none

    integer, intent(in) :: async_id
    real(r8) :: p_grad(nx,ny,nz)

    !$acc enter data async(async_id) create(p_grad)

    call pressure_gradient_x_gpu(p, p_grad, async_id)

    call solve_momentum_equation_gpu(u, v, w, p_grad, u1, ulf, u2, dmx, dmy, ev, async_id)

    call pressure_gradient_y_gpu(p, p_grad, async_id)

    call solve_momentum_equation_gpu(u, v, w, p_grad, v1, vlf, v2, dmx, dmy, ev, async_id)

    call solve_tracer_equation_gpu(u, v, w, t1, tlf, t2, dhx, dhy, hv(:,:,:,1), async_id)

    call solve_tracer_equation_gpu(u, v, w, s1, slf, s2, dhx, dhy, hv(:,:,:,2), async_id)

    !$acc exit data async(async_id) delete(p_grad)

  end subroutine solve_dynamic_equation_gpu

  subroutine pressure_gradient_x_gpu(p, px, async_id)
    use hyperlink, only: odx
    implicit none

    integer, intent(in) :: async_id
    real(r8), intent(in) :: &
      p(-1:nx+2,-1:ny+2,1:nz)

    real(r8), intent(out) :: &
      px(nx,ny,nz)

    integer :: i, j, k
    real(r8), dimension(0:nxf+1,0:ny+1,1:nz) :: xp, msk
    real(r8) :: tmp

    !$acc enter data async(async_id) create(xp, msk, px)

    !$acc parallel loop collapse(3) async(async_id)
    do k = 1, nz
      do j = 0, ny+1
        do i = 0, nxf+1
          xp(i,j,k) = 0.d0
          msk(i,j,k) = dble(iu(i,j,k))
          if (i >= xu_bgn-1 .and. i <= xu_end+1 .and. j >= 1 .and. j <= ny) then
            xp(i,j,k) = msk(i,j,k)*(p(i,j,k)-p(i-1,j,k))
          end if
        end do
      end do
    end do

    !$acc parallel loop collapse(3) async(async_id)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          tmp = o12*msk(i,j,k)*msk(i+1,j,k)
          px(i,j,k) = 0.5d0*(xp(i+1,j,k) + xp(i,j,k))  &
                    + tmp*(-xp(i-1,j,k) + xp(i,j,k)   &
                            +xp(i+1,j,k) - xp(i+2,j,k))
          px(i,j,k) = px(i,j,k)*odx(j)/rho_sw
        end do
      end do
    end do

    !$acc exit data async(async_id) delete(xp, msk, px)

  end subroutine pressure_gradient_x_gpu

  subroutine pressure_gradient_y_gpu(p, py, async_id)
    use hyperlink, only: ody
    implicit none

    integer, intent(in) :: async_id
    real(r8), intent(in) :: &
      p(-1:nx+2,-1:ny+2,1:nz)

    real(r8), intent(out) :: &
      py(nx,ny,nz)

    integer :: i, j, k
    real(r8), dimension(0:nx+1,0:nyf+1,1:nz) :: yp, msk
    real(r8) :: tmp

    !$acc enter data async(async_id) create(yp, msk)

    !$acc parallel loop collapse(3) async(async_id)
    do k = 1, nz
      do j = 0, nyf+1
        do i = 0, nx+1
          yp(i,j,k) = 0.d0
          msk(i,j,k) = dble(iv(i,j,k))
          if (i >= 1 .and. i <= nx .and. j >= yv_bgn-1 .and. j <= yv_end+1) then
            yp(i,j,k) = msk(i,j,k)*(p(i,j,k)-p(i,j-1,k))
          end if
        end do
      end do
    end do

    !$acc parallel loop collapse(3) async(async_id)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          tmp = o12*msk(i,j,k)*msk(i,j+1,k)
          py(i,j,k) = 0.5d0*(yp(i,j+1,k) + yp(i,j,k))  &
                    + tmp*(-yp(i,j-1,k) + yp(i,j,k)   &
                            +yp(i,j+1,k) - yp(i,j+2,k))
          py(i,j,k) = py(i,j,k)*ody(j)/rho_sw
        end do
      end do
    end do

    !$acc exit data async(async_id) delete(yp, msk)

  end subroutine pressure_gradient_y_gpu

  subroutine solve_momentum_equation_gpu(u, v, w, dp, q1, qlf, q2, hmix_x, hmix_y, vmix_z, async_id)
    use hyperlink, only: dt, odz, csv, odyv, ocs, ody, odx
    implicit none

    integer, intent(in) :: async_id
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

    call solve_momentum_or_tracer_equation_gpu(u, v, w, q1, qlf, q2, hmix_x, hmix_y, vmix_z, async_id, dp)

  end subroutine solve_momentum_equation_gpu

  subroutine solve_tracer_equation_gpu(u, v, w, q1, qlf, q2, hmix_x, hmix_y, vmix_z, async_id)
    use hyperlink, only: dt, odz, csv, odyv, ocs, ody, odx
    implicit none

    integer, intent(in) :: async_id
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

    call solve_momentum_or_tracer_equation_gpu(u, v, w, q1, qlf, q2, hmix_x, hmix_y, vmix_z, async_id)

  end subroutine solve_tracer_equation_gpu

  subroutine solve_momentum_or_tracer_equation_gpu(u, v, w, q1, qlf, q2, hmix_x, hmix_y, vmix_z, async_id, dp)
    use hyperlink, only: dt, odz, csv, odyv, ocs, ody, odx
    implicit none

    integer, intent(in) :: async_id
    real(r8), intent(inout) ::   &
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

    real(r8), intent(in), optional ::  &
      dp ( 1:nx,   1:ny,  1:nz)

    real(r8), dimension(nx,ny) :: &
      dqx_dx, dqy_dy, dqz_dz, qz

    real(r8) ::   &
      q2_out(-1:nx+2,-1:ny+2,1:nz)

    real(r8) :: dtin
    integer :: i, j, k
    logical :: is_momentum

    ! Variables for inlined compute_vertical_fluxes
    real(r8) :: vert_tmp, vert_qz_b_val, vert_scr_val

    ! Variables for inlined compute_latitudinal_fluxes
    real(r8) :: lati_tmp, lati_odyj
    real(r8) :: lati_msk_val, lati_qy_j, lati_qy_jp1
    real(r8) :: lati_scr_j, lati_scr_jp1

    ! Variables for inlined compute_longitudinal_fluxes
    real(r8) :: long_tmp
    real(r8) :: long_msk_val, long_qx_i, long_qx_ip1
    real(r8) :: long_scr_i, long_scr_ip1

    !$acc enter data async(async_id) create(dqx_dx, dqy_dy, dqz_dz, qz, q2_out)

    !$acc parallel loop collapse(2) async(async_id)
    do j = 1, ny
      do i = 1, nx
        select case(lfsrf)
        case(0)
          qz(i,j) = 0.d0
        case(1)
          qz(i,j) = w(i,j,1)*q2(i,j,1)
        end select

        !$acc loop seq
        do k = 1, nz
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute_vertical_fluxes
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if(k .eq. nz) then
            vert_qz_b_val = 0.d0
          else
            vert_scr_val = 6.d0*(q2(i,j,k) + q2(i,j,k+1))
            if(k .gt. 1 .and. k .lt. nz-1) then
              vert_tmp = dble(iw(i,j,k)*iw(i,j,k+2))
              vert_scr_val = vert_scr_val + vert_tmp*(-qlf(i,j,k-1)+q2(i,j,k)+q2(i,j,k+1)-q2(i,j,k+2))
            end if
            vert_scr_val = o12*vert_scr_val
            vert_tmp = dble(iw(i,j,k+1))
            vert_qz_b_val = vert_tmp*(w(i,j,k+1)*vert_scr_val - vmix_z(i,j,k)*(q1(i,j,k+1)-q1(i,j,k)))
          end if
          dqz_dz(i,j) = (vert_qz_b_val - qz(i,j))*odz(k)
          qz(i,j) = vert_qz_b_val

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute_longitudinal_fluxes
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! Calculate scr(i,j)
          if (i >= xu_bgn .and. i <= xu_end) then
            long_scr_i = 6.d0*(q2(i-1,j,k) + q2(i,j,k))
            long_msk_val = dble(iu(i-1,j,k))*dble(iu(i+1,j,k))
            long_scr_i = long_scr_i + long_msk_val*(-q2(i-2,j,k) + q2(i-1,j,k) + q2(i,j,k) - q2(i+1,j,k))
            long_scr_i = o12*long_scr_i
          else
            long_scr_i = 0.d0
          end if

          ! Calculate scr(i+1,j)
          if (i+1 >= xu_bgn .and. i+1 <= xu_end) then
            long_scr_ip1 = 6.d0*(q2(i,j,k) + q2(i+1,j,k))
            long_msk_val = dble(iu(i,j,k))*dble(iu(i+2,j,k))
            long_scr_ip1 = long_scr_ip1 + long_msk_val*(-q2(i-1,j,k) + q2(i,j,k) + q2(i+1,j,k) - q2(i+2,j,k))
            long_scr_ip1 = o12*long_scr_ip1
          else
            long_scr_ip1 = 0.d0
          end if

          ! Calculate qx(i,j)
          long_msk_val = dble(iu(i,j,k))
          long_qx_i = long_msk_val*(u(i,j,k)*long_scr_i - hmix_x(i,j,k)*odx(j)*(q1(i,j,k) - q1(i-1,j,k)))

          ! Calculate qx(i+1,j)
          long_msk_val = dble(iu(i+1,j,k))
          long_qx_ip1 = long_msk_val*(u(i+1,j,k)*long_scr_ip1 - hmix_x(i+1,j,k)*odx(j)*(q1(i+1,j,k) - q1(i,j,k)))

          dqx_dx(i,j) = (long_qx_ip1 - long_qx_i)*odx(j)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! compute_latitudinal_fluxes
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          lati_odyj = ocs(j)*ody(j)

          ! Calculate scr(i,j)
          if (j >= yv_bgn .and. j <= yv_end) then
            lati_scr_j = 6.d0*(q2(i,j-1,k) + q2(i,j,k))
            lati_msk_val = dble(iv(i,j-1,k))*dble(iv(i,j+1,k))
            lati_scr_j = lati_scr_j + lati_msk_val*(-q2(i,j-2,k) + q2(i,j-1,k) + q2(i,j,k) - q2(i,j+1,k))
            lati_scr_j = o12*lati_scr_j
          else
            lati_scr_j = 0.d0
          end if

          ! Calculate scr(i,j+1)
          if (j+1 >= yv_bgn .and. j+1 <= yv_end) then
            lati_scr_jp1 = 6.d0*(q2(i,j,k) + q2(i,j+1,k))
            lati_msk_val = dble(iv(i,j,k))*dble(iv(i,j+2,k))
            lati_scr_jp1 = lati_scr_jp1 + lati_msk_val*(-q2(i,j-1,k) + q2(i,j,k) + q2(i,j+1,k) - q2(i,j+2,k))
            lati_scr_jp1 = o12*lati_scr_jp1
          else
            lati_scr_jp1 = 0.d0
          end if

          ! Calculate qy(i,j)
          lati_msk_val = dble(iv(i,j,k))
          lati_tmp = csv(j)*lati_msk_val
          lati_qy_j = lati_tmp*(v(i,j,k)*lati_scr_j - hmix_y(i,j,k)*odyv(j)*(q1(i,j,k) - q1(i,j-1,k)))

          ! Calculate qy(i,j+1)
          lati_msk_val = dble(iv(i,j+1,k))
          lati_tmp = csv(j+1)*lati_msk_val
          lati_qy_jp1 = lati_tmp*(v(i,j+1,k)*lati_scr_jp1 - hmix_y(i,j+1,k)*odyv(j+1)*(q1(i,j+1,k) - q1(i,j,k)))

          dqy_dy(i,j) = (lati_qy_jp1 - lati_qy_j)*lati_odyj
          dtin = dt*dble(in(i,j,k))

          ! Check if this is a momentum equation (dp is present) or tracer equation
          is_momentum = present(dp)
          if (is_momentum) then
            q2_out(i,j,k) = q1(i,j,k) - dtin*(dp(i,j,k) + dqx_dx(i,j) + dqy_dy(i,j) + dqz_dz(i,j))
          else
            q2_out(i,j,k) = q1(i,j,k) - dtin*(dqx_dx(i,j) + dqy_dy(i,j) + dqz_dz(i,j))
          end if
        end do
      end do
    end do

    !$acc kernels async(async_id)
    q2(1:nx,1:ny,1:nz) = q2_out(1:nx,1:ny,1:nz)
    !$acc end kernels

    !$acc exit data async(async_id) delete(dqx_dx, dqy_dy, dqz_dz, qz, q2_out)

  end subroutine solve_momentum_or_tracer_equation_gpu

end module timcom_drv_gpu

