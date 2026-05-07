!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_impose_temp_salt_nudging
  implicit none

  call timcom_initialize_gpu
  call impose_temp_salt_nudging_unit
  call timcom_finalize_gpu

end program test_impose_temp_salt_nudging

subroutine impose_temp_salt_nudging_unit
  use timcom_drv, only: impose_temp_salt_nudging
  use timcom_drv_gpu, only: impose_temp_salt_nudging_gpu
  use hyperlink, only: nx, ny, nz, t2, t_nudge, t_clim, t_da, opt_da
  implicit none

  real(8), dimension(:), allocatable :: t2_cpu
  real(8), dimension(:), allocatable :: t2_gpu
  real(8), dimension(:,:,:), allocatable :: t2_initial
  real :: time1, time2, ct, gt
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i
  integer :: syng_mon

  ! Initialize variables
  allocate(t2_initial(size(t2,1), size(t2,2), size(t2,3)))
  allocate(t2_cpu(size(t2)))
  allocate(t2_gpu(size(t2)))

  ! Initialize input parameters
  call random_number(t2)
  call random_number(t_nudge)
  call random_number(t_clim)
  call random_number(t_da)

  ! Set opt_da to 1 to test the data assimilation branch
  opt_da = 1

  ! Set a random month for the synoptic month parameter
  syng_mon = 6  ! Example: June

  ! Store initial values
  t2_initial = t2

  ! CPU execution
  ct = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call impose_temp_salt_nudging(syng_mon)
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do
  t2_cpu = reshape(t2, [size(t2)])

  ! Restore initial values before GPU execution
  t2 = t2_initial

  ! GPU execution
  !$acc enter data async(async_id) &
  !$acc& copyin(t2(1:nx,1:ny,1), t_nudge(1:nx,1:ny), t_da(1:nx,1:ny), &
  !$acc&        nx, ny)
  gt = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call impose_temp_salt_nudging_gpu(syng_mon, async_id)
    !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(t2(1:nx,1:ny,1)), &
  !$acc& delete(t_nudge(1:nx,1:ny), t_da(1:nx,1:ny), nx, ny)
  !$acc wait(async_id)

  t2_gpu = reshape(t2, [size(t2)])

  ! Assertions
  call assert_allclose(t2_gpu, size(t2_gpu), t2_cpu, size(t2_cpu), 1e-10, 1e-10, "Array t2")

  ! Print timing results
  write (*, *) 'Timing for CPU & GPU: ', ct, gt
  write (*, *) 'Speedup ratio: ', ct/gt

end subroutine impose_temp_salt_nudging_unit
