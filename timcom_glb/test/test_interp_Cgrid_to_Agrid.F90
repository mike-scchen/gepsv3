!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_interp_Cgrid_to_Agrid
  implicit none

  call timcom_initialize_gpu
  call interp_Cgrid_to_Agrid_unit
  call timcom_finalize_gpu

end program test_interp_Cgrid_to_Agrid

subroutine interp_Cgrid_to_Agrid_unit
  use timcom_drv, only: interp_Cgrid_to_Agrid
  use timcom_drv_gpu, only: interp_Cgrid_to_Agrid_gpu
  use hyperlink, only: u, u2, v, v2, &
                       m_comm_cart, r8type3d, &
                       nbid, ndim, symm_np, &
                       nx, ny, nz, nxf, nyf, nzf, &
                       iu, iv, o24, in, myid
  implicit none

  real(8), dimension(:,:,:), allocatable :: u2_cpu, v2_cpu
  real(8), dimension(:,:,:), allocatable :: u2_gpu, v2_gpu
  real(8), dimension(:,:,:), allocatable :: u2_initial, v2_initial
  real(8), dimension(:,:,:), allocatable :: iu_r8, iv_r8, in_r8
  real :: time1, time2, ct, gt
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i
  integer :: seed_size
  integer, allocatable :: seed(:)

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 12345 + myid
  call random_seed(put=seed)

  ! Initialize variables
  allocate(u2_initial(size(u2,1), size(u2,2), size(u2,3)), v2_initial(size(v2,1), size(v2,2), size(v2,3)))
  allocate(u2_cpu(size(u2,1), size(u2,2), size(u2,3)), v2_cpu(size(v2,1), size(v2,2), size(v2,3)))
  allocate(u2_gpu(size(u2,1), size(u2,2), size(u2,3)), v2_gpu(size(v2,1), size(v2,2), size(v2,3)))
  allocate(iu_r8(size(iu,1), size(iu,2), size(iu,3)))
  allocate(iv_r8(size(iv,1), size(iv,2), size(iv,3)))
  allocate(in_r8(size(in,1), size(in,2), size(in,3)))

  ! Initialize input parameters
  call random_number(u)
  call random_number(v)
  call random_number(u2)
  call random_number(v2)
  call random_number(iu_r8)
  iu = nint(iu_r8)
  call random_number(iv_r8)
  iv = nint(iv_r8)
  call random_number(in_r8)
  in = nint(in_r8)

  ! Store initial values
  u2_initial = u2
  v2_initial = v2

  ! CPU execution
  ct = 0.0
  do i = 1, num_steps
    u2 = u2_initial
    v2 = v2_initial
    call cpu_time(time1)
    call interp_Cgrid_to_Agrid
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do
  u2_cpu = u2
  v2_cpu = v2

  ! Restore initial values before GPU execution
  u2 = u2_initial
  v2 = v2_initial

  ! GPU execution
  gt = 0.0
  !$acc enter data async(async_id) &
  !$acc& copyin(u, u2, iu, v, v2, iv, in, nx, ny, nz, nxf, nyf)
  do i = 1, num_steps
    !$acc update device(u2, v2) async(async_id)
    !$acc wait(async_id)
    call cpu_time(time1)
    call interp_Cgrid_to_Agrid_gpu(async_id)
    !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(u2, v2), &
  !$acc& delete(u, v, iu, iv, in, nx, ny, nz, nxf, nyf, nzf)
  !$acc wait(async_id)

  u2_gpu = u2
  v2_gpu = v2

  ! Assertions
  call assert_allclose_r8type3d(u2_gpu, u2_cpu, nx, ny, nz, 1e-10, 1e-10, "Array u2")
  call assert_allclose_r8type3d(v2_gpu, v2_cpu, nx, ny, nz, 1e-10, 1e-10, "Array v2")

  ! Print timing results
  write (*, *) 'Timing for CPU & GPU: ', ct, gt
  write (*, *) 'Speedup ratio: ', ct/gt

end subroutine interp_Cgrid_to_Agrid_unit
