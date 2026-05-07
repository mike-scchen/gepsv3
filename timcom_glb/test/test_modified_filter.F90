!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_modified_filter
  implicit none

  call timcom_initialize_gpu
  call modified_filter_unit
  call timcom_finalize_gpu

end program test_modified_filter

subroutine modified_filter_unit
  use timcom_drv, only: modified_filter
  use timcom_drv_gpu, only: modified_filter_gpu
  use hyperlink, only: u1, ulf, u2, &
                       v1, vlf, v2, &
                       t1, tlf, t2, &
                       s1, slf, s2, &
                       symm_np, nx, ny, nz, myid
  implicit none

  real(8), dimension(:,:,:), allocatable :: u1_initial, ulf_initial, v1_initial, vlf_initial
  real(8), dimension(:,:,:), allocatable :: t1_initial, tlf_initial, s1_initial, slf_initial
  real(8), dimension(:,:,:), allocatable :: u1_cpu, ulf_cpu, v1_cpu, vlf_cpu, t1_cpu, tlf_cpu, s1_cpu, slf_cpu
  real(8), dimension(:,:,:), allocatable :: u1_gpu, ulf_gpu, v1_gpu, vlf_gpu, t1_gpu, tlf_gpu, s1_gpu, slf_gpu
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
  allocate(u1_initial(size(u1,1), size(u1,2), size(u1,3)), ulf_initial(size(ulf,1), size(ulf,2), size(ulf,3)))
  allocate(v1_initial(size(v1,1), size(v1,2), size(v1,3)), vlf_initial(size(vlf,1), size(vlf,2), size(vlf,3)))
  allocate(t1_initial(size(t1,1), size(t1,2), size(t1,3)), tlf_initial(size(tlf,1), size(tlf,2), size(tlf,3)))
  allocate(s1_initial(size(s1,1), size(s1,2), size(s1,3)), slf_initial(size(slf,1), size(slf,2), size(slf,3)))

  allocate(u1_cpu(size(u1,1), size(u1,2), size(u1,3)), ulf_cpu(size(ulf,1), size(ulf,2), size(ulf,3)))
  allocate(v1_cpu(size(v1,1), size(v1,2), size(v1,3)), vlf_cpu(size(vlf,1), size(vlf,2), size(vlf,3)))
  allocate(t1_cpu(size(t1,1), size(t1,2), size(t1,3)), tlf_cpu(size(tlf,1), size(tlf,2), size(tlf,3)))
  allocate(s1_cpu(size(s1,1), size(s1,2), size(s1,3)), slf_cpu(size(slf,1), size(slf,2), size(slf,3)))
  allocate(u1_gpu(size(u1,1), size(u1,2), size(u1,3)), ulf_gpu(size(ulf,1), size(ulf,2), size(ulf,3)))
  allocate(v1_gpu(size(v1,1), size(v1,2), size(v1,3)), vlf_gpu(size(vlf,1), size(vlf,2), size(vlf,3)))
  allocate(t1_gpu(size(t1,1), size(t1,2), size(t1,3)), tlf_gpu(size(tlf,1), size(tlf,2), size(tlf,3)))
  allocate(s1_gpu(size(s1,1), size(s1,2), size(s1,3)), slf_gpu(size(slf,1), size(slf,2), size(slf,3)))

  ! Initialize input parameters with random values
  call random_number(u1)
  call random_number(ulf)
  call random_number(u2)
  call random_number(v1)
  call random_number(vlf)
  call random_number(v2)
  call random_number(t1)
  call random_number(tlf)
  call random_number(t2)
  call random_number(s1)
  call random_number(slf)
  call random_number(s2)

  ! Store initial values
  u1_initial = u1
  ulf_initial = ulf
  v1_initial = v1
  vlf_initial = vlf
  t1_initial = t1
  tlf_initial = tlf
  s1_initial = s1
  slf_initial = slf

  ! CPU execution
  ct = 0.0
  do i = 1, num_steps
    u1 = u1_initial
    ulf = ulf_initial
    v1 = v1_initial
    vlf = vlf_initial
    t1 = t1_initial
    tlf = tlf_initial
    s1 = s1_initial
    slf = slf_initial
    call cpu_time(time1)
    call modified_filter
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do

  u1_cpu = u1
  ulf_cpu = ulf
  v1_cpu = v1
  vlf_cpu = vlf
  t1_cpu = t1
  tlf_cpu = tlf
  s1_cpu = s1
  slf_cpu = slf

  ! Restore initial values before GPU execution
  u1 = u1_initial
  ulf = ulf_initial
  v1 = v1_initial
  vlf = vlf_initial
  t1 = t1_initial
  tlf = tlf_initial
  s1 = s1_initial
  slf = slf_initial

  ! GPU execution
  !$acc enter data async(async_id) &
  !$acc& copyin(u1, ulf, u2, v1, vlf, v2, t1, tlf, t2, s1, slf, s2, &
  !$acc&        nx, ny, nz)
  gt = 0.0
  do i = 1, num_steps
    !$acc update device(u1, ulf, v1, vlf, t1, tlf, s1, slf) async(async_id)
    !$acc wait(async_id)
    call cpu_time(time1)
    call modified_filter_gpu(async_id)
    !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(u1, ulf, v1, vlf, t1, tlf, s1, slf) &
  !$acc& delete(u2, v2, t2, s2, nx, ny, nz)
  !$acc wait(async_id)

  u1_gpu = u1
  ulf_gpu = ulf
  v1_gpu = v1
  vlf_gpu = vlf
  t1_gpu = t1
  tlf_gpu = tlf
  s1_gpu = s1
  slf_gpu = slf

  ! Assertions
  call assert_allclose_r8type3d(u1_gpu, u1_cpu, nx, ny, nz, 1e-10, 1e-10, "Array u1")
  call assert_allclose_r8type3d(ulf_gpu, ulf_cpu, nx, ny, nz, 1e-10, 1e-10, "Array ulf")
  call assert_allclose_r8type3d(v1_gpu, v1_cpu, nx, ny, nz, 1e-10, 1e-10, "Array v1")
  call assert_allclose_r8type3d(vlf_gpu, vlf_cpu, nx, ny, nz, 1e-10, 1e-10, "Array vlf")
  call assert_allclose_r8type3d(t1_gpu, t1_cpu, nx, ny, nz, 1e-10, 1e-10, "Array t1")
  call assert_allclose_r8type3d(tlf_gpu, tlf_cpu, nx, ny, nz, 1e-10, 1e-10, "Array tlf")
  call assert_allclose_r8type3d(s1_gpu, s1_cpu, nx, ny, nz, 1e-10, 1e-10, "Array s1")
  call assert_allclose_r8type3d(slf_gpu, slf_cpu, nx, ny, nz, 1e-10, 1e-10, "Array slf")

  ! Print timing results
  write (*, *) 'Timing for CPU & GPU: ', ct, gt
  write (*, *) 'Speedup ratio: ', ct/gt

end subroutine modified_filter_unit
