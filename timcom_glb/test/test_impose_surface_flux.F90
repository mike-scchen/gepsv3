!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_impose_surface_flux
  implicit none

  call timcom_initialize_gpu
  call impose_surface_flux_unit
  call timcom_finalize_gpu

end program test_impose_surface_flux

subroutine impose_surface_flux_unit
  use timcom_drv, only: impose_surface_flux
  use timcom_drv_gpu, only: impose_surface_flux_gpu
  use hyperlink, only: nx, ny, nz, dt, odz, u2, v2, t2, s2, smft, stf, sw_trans, shf_qsw
  implicit none

  real(8), dimension(:), allocatable :: u2_cpu, v2_cpu, t2_cpu, s2_cpu
  real(8), dimension(:), allocatable :: u2_gpu, v2_gpu, t2_gpu, s2_gpu
  real(8), dimension(:,:,:), allocatable :: u2_initial, v2_initial, t2_initial, s2_initial
  real :: time1, time2, ct, gt
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  ! Initialize variables
  allocate(u2_initial(size(u2,1), size(u2,2), size(u2,3)), v2_initial(size(v2,1), size(v2,2), size(v2,3)), &
           t2_initial(size(t2,1), size(t2,2), size(t2,3)), s2_initial(size(s2,1), size(s2,2), size(s2,3)))
  allocate(u2_cpu(size(u2)), v2_cpu(size(v2)), t2_cpu(size(t2)), s2_cpu(size(s2)))
  allocate(u2_gpu(size(u2)), v2_gpu(size(v2)), t2_gpu(size(t2)), s2_gpu(size(s2)))
  call random_number(dt)
  call random_number(odz)
  call random_number(u2)
  call random_number(v2)
  call random_number(t2)
  call random_number(s2)
  call random_number(smft)
  call random_number(stf)
  call random_number(sw_trans)
  call random_number(shf_qsw)

  ! Store initial values
  u2_initial = u2
  v2_initial = v2
  t2_initial = t2
  s2_initial = s2

  ! CPU execution
  ct = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call impose_surface_flux
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do
  u2_cpu = reshape(u2, [size(u2)])
  v2_cpu = reshape(v2, [size(v2)])
  t2_cpu = reshape(t2, [size(t2)])
  s2_cpu = reshape(s2, [size(s2)])

  ! Restore initial values before GPU execution
  u2 = u2_initial
  v2 = v2_initial
  t2 = t2_initial
  s2 = s2_initial

  ! GPU execution
  !$acc enter data async(async_id) &
  !$acc& copyin(u2(1:nx,1:ny,1), v2(1:nx,1:ny,1), &
  !$acc& t2(1:nx,1:ny,1:nz), s2(1:nx,1:ny,1), &
  !$acc& smft(1:nx,1:ny,1:2), stf(1:nx,1:ny,1:2), &
  !$acc& odz, sw_trans, shf_qsw, nx, ny, nz, dt)
  gt = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call impose_surface_flux_gpu(async_id)
    !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(u2(1:nx,1:ny,1), v2(1:nx,1:ny,1), &
  !$acc& t2(1:nx,1:ny,1:nz), s2(1:nx,1:ny,1)), &
  !$acc& delete(nx, ny, nz, dt, odz, smft, stf, sw_trans, shf_qsw)
  !$acc wait(async_id)

  u2_gpu = reshape(u2, [size(u2)])
  v2_gpu = reshape(v2, [size(v2)])
  t2_gpu = reshape(t2, [size(t2)])
  s2_gpu = reshape(s2, [size(s2)])

  ! Assertions
  call assert_allclose(u2_gpu, size(u2_gpu), u2_cpu, size(u2_cpu), 1e-10, 1e-10, "Array u2")
  call assert_allclose(v2_gpu, size(v2_gpu), v2_cpu, size(v2_cpu), 1e-10, 1e-10, "Array v2")
  call assert_allclose(t2_gpu, size(t2_gpu), t2_cpu, size(t2_cpu), 1e-10, 1e-10, "Array t2")
  call assert_allclose(s2_gpu, size(s2_gpu), s2_cpu, size(s2_cpu), 1e-10, 1e-10, "Array s2")

  ! Print timing results
  write (*, *) 'Timing for CPU & GPU: ', ct, gt
  write (*, *) 'Speedup ratio: ', ct/gt

end subroutine impose_surface_flux_unit

