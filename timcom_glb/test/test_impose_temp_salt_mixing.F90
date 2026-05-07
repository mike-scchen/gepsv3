!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_impose_temp_salt_mixing
  implicit none

  call timcom_initialize_gpu
  call impose_temp_salt_mixing_unit
  call timcom_finalize_gpu

end program test_impose_temp_salt_mixing

subroutine impose_temp_salt_mixing_unit
  use timcom_drv, only: impose_temp_salt_mixing
  use timcom_drv_gpu, only: impose_temp_salt_mixing_gpu
  use hyperlink, only: nx, ny, nz, dt, t2, s2, kpp_src, opt_windmix, in
  implicit none

  real(8), dimension(:), allocatable :: t2_cpu, s2_cpu
  real(8), dimension(:), allocatable :: t2_gpu, s2_gpu
  real(8), dimension(:,:,:), allocatable :: t2_initial, s2_initial
  real(8), dimension(:,:,:), allocatable :: in_r8
  real :: time1, time2, ct, gt
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  ! Initialize variables
  allocate(t2_initial(size(t2,1), size(t2,2), size(t2,3)), s2_initial(size(s2,1), size(s2,2), size(s2,3)))
  allocate(t2_cpu(size(t2)), s2_cpu(size(s2)))
  allocate(t2_gpu(size(t2)), s2_gpu(size(s2)))
  allocate(in_r8(size(in,1), size(in,2), size(in,3)))

  call random_number(dt)
  call random_number(t2)
  call random_number(s2)
  call random_number(in_r8)
  in = nint(in_r8)
  call random_number(kpp_src)

  ! Set opt_windmix to use KPP (value 1)
  opt_windmix = 1

  ! Store initial values
  t2_initial = t2
  s2_initial = s2

  ! CPU execution
  ct = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call impose_temp_salt_mixing
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do
  t2_cpu = reshape(t2, [size(t2)])
  s2_cpu = reshape(s2, [size(s2)])

  ! Restore initial values before GPU execution
  t2 = t2_initial
  s2 = s2_initial

  ! GPU execution
  !$acc enter data async(async_id) &
  !$acc& copyin(t2(1:nx,1:ny,1:nz), s2(1:nx,1:ny,1:nz), &
  !$acc&        kpp_src(1:nx,1:ny,1:nz,1:2), in(1:nx,1:ny,1:nz), &
  !$acc&        dt, nx, ny, nz)
  gt = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call impose_temp_salt_mixing_gpu(async_id)
    !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(t2(1:nx,1:ny,1:nz), s2(1:nx,1:ny,1:nz)), &
  !$acc& delete(kpp_src(1:nx,1:ny,1:nz,1:2), in(1:nx,1:ny,1:nz), &
  !$acc&        dt, nx, ny, nz)
  !$acc wait(async_id)

  t2_gpu = reshape(t2, [size(t2)])
  s2_gpu = reshape(s2, [size(s2)])

  ! Assertions
  call assert_allclose(t2_gpu, size(t2_gpu), t2_cpu, size(t2_cpu), 1e-10, 1e-10, "Array t2")
  call assert_allclose(s2_gpu, size(s2_gpu), s2_cpu, size(s2_cpu), 1e-10, 1e-10, "Array s2")

  ! Print timing results
  write (*, *) 'Timing for CPU & GPU: ', ct, gt
  write (*, *) 'Speedup ratio: ', ct/gt

end subroutine impose_temp_salt_mixing_unit

