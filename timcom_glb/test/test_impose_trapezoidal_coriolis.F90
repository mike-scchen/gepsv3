program test_impose_trapezoidal_coriolis
  implicit none

  call timcom_initialize_gpu
  call impose_trapezoidal_coriolis_unit
  call timcom_finalize_gpu

end program test_impose_trapezoidal_coriolis

subroutine impose_trapezoidal_coriolis_unit
  use timcom_drv, only: impose_trapezoidal_coriolis
  use timcom_drv_gpu, only: impose_trapezoidal_coriolis_gpu
  use hyperlink, only: u1, u2, v1, v2, dt, curv_f, tanphi, ulf
  use hyperlink, only: nx, ny, nz, nzf
  implicit none

  real(8), dimension(:), allocatable :: u2_cpu, u2_gpu, v2_cpu, v2_gpu
  real(8), dimension(:,:,:), allocatable :: u2_tmp, v2_tmp
  real :: time1, time2, ct, gt
  integer, parameter :: async_id = 1 !UT要給1
!  integer, parameter :: async_id = -1 !測時間時要給-1

  integer, parameter :: num_steps = 16
  integer :: i

  allocate(u2_tmp(nx,ny,nz), v2_tmp(nx,ny,nz))
  call random_number(u1)
  call random_number(v1)
  call random_number(u2_tmp)
  call random_number(v2_tmp)
  call random_number(dt)
  call random_number(curv_f)
  call random_number(tanphi)
  call random_number(ulf)
  allocate(u2_cpu(size(u2)), u2_gpu(size(u2)))
  allocate(v2_cpu(size(v2)), v2_gpu(size(v2)))
  u2=u2_tmp
  v2=v2_tmp
  ct = 0.
  do i = 1, num_steps
    call cpu_time(time1)
    call impose_trapezoidal_coriolis
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do
  u2_cpu = reshape(u2, [size(u2)])
  v2_cpu = reshape(v2, [size(v2)])

  u2=u2_tmp
  v2=v2_tmp
  !$acc data async(async_id) copy(u2(1:nx,1:ny,1:nz), v2(1:nx,1:ny,1:nz)) &
  !$acc& copyin(nx, ny, nz, dt, tanphi(:)) &
  !$acc& copyin(u1(1:nx,1:ny,1:nz),v1(1:nx,1:ny,1:nz), curv_f(:), ulf(1:nx,1:ny,1:nz))
  gt = 0.
  do i = 1, num_steps
    call cpu_time(time1)
    call impose_trapezoidal_coriolis_gpu(async_id)
    !$acc wait(async_id)
    call cpu_time(time2)  
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc end data
  !$acc wait(async_id)
  u2_gpu = reshape(u2, [size(u2)])
  v2_gpu = reshape(v2, [size(v2)])

  call assert_allclose(u2_gpu, size(u2_gpu), u2_cpu, size(u2_cpu), 1e-10, 1e-10, "Array u2")
  call assert_allclose(v2_gpu, size(v2_gpu), v2_cpu, size(v2_cpu), 1e-10, 1e-10, "Array v2")
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt

end subroutine impose_trapezoidal_coriolis_unit

