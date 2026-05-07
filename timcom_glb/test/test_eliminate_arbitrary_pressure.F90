!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_eliminate_arbitrary_pressure
  implicit none

  call timcom_initialize_gpu
  call eliminate_arbitrary_pressure_unit
  call timcom_finalize_gpu

end program test_eliminate_arbitrary_pressure

subroutine eliminate_arbitrary_pressure_unit
  use hyperlink, only: p0, nx_grid, ny_grid, opt_arbr_p0, area, total_area, in, nx, ny
  use timcom_drv, only: eliminate_arbitrary_pressure
  use timcom_drv_gpu, only: eliminate_arbitrary_pressure_gpu
  use timcom_const
  implicit none

  real(r8), dimension(:,:,:), allocatable :: in_r8
  real(r8), dimension(:), allocatable :: p0_cpu, p0_gpu
  real(r8), dimension(:,:), allocatable :: p0_initial
  real(r8) :: time1, time2, ct, gt
  integer :: i
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16

  ! Initialize variables
  allocate(p0_initial(size(p0,1), size(p0,2)))
  allocate(p0_cpu(size(p0)), p0_gpu(size(p0)))
  allocate(in_r8(size(in,1), size(in,2), size(in,3)))

  ! Initialize test data
  call random_number(p0)
  call random_number(area)
  area = area * 1000.0_r8
  total_area = sum(area)
  call random_number(in_r8)
  in = nint(in_r8)
  opt_arbr_p0 = 1

  ! Store initial values
  p0_initial = p0

  ! CPU execution
  ct = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call eliminate_arbitrary_pressure()
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do

  ! Store CPU results
  p0_cpu = reshape(p0, [size(p0)])

  ! Restore initial values before GPU execution
  p0 = p0_initial

  ! GPU execution
  gt = 0.0
  !$acc enter data async(async_id) &
  !$acc& copyin(p0(1:nx,1:ny), area(1:nx,1:ny), in(1:nx,1:ny,1), total_area, nx, ny)
  do i = 1, num_steps
    call cpu_time(time1)
    call eliminate_arbitrary_pressure_gpu(async_id)
    !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(p0(1:nx,1:ny)), &
  !$acc& delete(in(1:nx,1:ny,1), area(1:nx,1:ny), total_area, nx, ny)
  !$acc wait(async_id)

  p0_gpu = reshape(p0, [size(p0)])

  ! Assertions
  call assert_allclose(p0_gpu, size(p0_gpu), p0_cpu, size(p0_cpu), &
                      1e-10, 1e-10, "Array p0")

  ! Print timing results
  write (*, *) 'Timing for CPU & GPU: ', ct, gt
  write (*, *) 'Speedup ratio: ', ct/gt

  ! Clean up
  deallocate(p0_initial, p0_cpu, p0_gpu)

end subroutine eliminate_arbitrary_pressure_unit
