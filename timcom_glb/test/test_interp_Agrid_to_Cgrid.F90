!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_interp_Agrid_to_Cgrid
  implicit none

  call timcom_initialize_gpu
  call interp_Agrid_to_Cgrid_unit
  call timcom_finalize_gpu

end program test_interp_Agrid_to_Cgrid

subroutine interp_Agrid_to_Cgrid_unit
  use timcom_drv, only: interp_Agrid_to_Cgrid
  use timcom_drv_gpu, only: interp_Agrid_to_Cgrid_gpu
  use hyperlink, only: u, u2, v, v2, &
                       m_comm_cart, r8type3du, r8type3dv, &
                       nbid, ndim, symm_np, &
                       nx, ny, nz, nxf, nyf, nzf, &
                       xu_bgn, xu_end, yv_bgn, yv_end, &
                       iu, iv, o12
  implicit none

  real(8), dimension(:), allocatable :: u_cpu, v_cpu
  real(8), dimension(:), allocatable :: u_gpu, v_gpu
  real(8), dimension(:,:,:), allocatable :: u_initial, v_initial
  real(8), dimension(:,:,:), allocatable :: iu_r8, iv_r8
  real :: time1, time2, ct, gt
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  symm_np = .true.

  ! Initialize variables
  allocate(u_initial(size(u,1), size(u,2), size(u,3)), v_initial(size(v,1), size(v,2), size(v,3)))
  allocate(u_cpu(size(u)), v_cpu(size(v)))
  allocate(u_gpu(size(u)), v_gpu(size(v)))
  allocate(iu_r8(size(iu,1), size(iu,2), size(iu,3)), iv_r8(size(iv,1), size(iv,2), size(iv,3)))

  ! Initialize input parameters
  call random_number(u)
  call random_number(v)
  call random_number(u2)
  call random_number(v2)
  call random_number(iu_r8)
  iu = iu_r8
  call random_number(iv_r8)
  iv = iv_r8

  ! Store initial values
  u_initial = u
  v_initial = v

  ! CPU execution
  ct = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call interp_Agrid_to_Cgrid
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do
  u_cpu = reshape(u, [size(u)])
  v_cpu = reshape(v, [size(v)])

  ! Restore initial values before GPU execution
  u = u_initial
  v = v_initial

  ! GPU execution
  !$acc enter data async(async_id) &
  !$acc& copyin(u, u2, iu, v, v2, iv, &
  !$acc&        xu_bgn, xu_end, yv_bgn, yv_end, &
  !$acc&        nx, ny, nz, nxf, nyf)
  gt = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call interp_Agrid_to_Cgrid_gpu(async_id)
    !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(u, v), &
  !$acc& delete(u2, v2, iu, iv, xu_bgn, xu_end, &
  !$acc&        yv_bgn, yv_end, nx, ny, nz, nxf, nyf)
  !$acc wait(async_id)

  u_gpu = reshape(u, [size(u)])
  v_gpu = reshape(v, [size(v)])

  ! Assertions
  call assert_allclose(u_gpu, size(u_gpu), u_cpu, size(u_cpu), 1e-10, 1e-10, "Array u")
  call assert_allclose(v_gpu, size(v_gpu), v_cpu, size(v_cpu), 1e-10, 1e-10, "Array v")

  ! Print timing results
  write (*, *) 'Timing for CPU & GPU: ', ct, gt
  write (*, *) 'Speedup ratio: ', ct/gt

end subroutine interp_Agrid_to_Cgrid_unit
