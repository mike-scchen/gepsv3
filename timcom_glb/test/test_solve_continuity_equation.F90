!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_solve_continuity_equation
  implicit none

  call timcom_initialize_gpu
  call solve_continuity_equation_unit
  call timcom_finalize_gpu

end program test_solve_continuity_equation

subroutine solve_continuity_equation_unit
  use timcom_drv, only: solve_continuity_equation
  use timcom_drv_gpu, only: solve_continuity_equation_gpu
  use hyperlink, only: u, v, w, odx, ody, odz, ocs, csv, iw, kb
  use hyperlink, only: nx, ny, nzf, nxf, nz, nyf
  implicit none

  real(8), dimension(:), allocatable :: w_cpu, w_gpu
  real :: time1, time2, ct, gt
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  call random_number(u)
  call random_number(v)
  call random_number(w)
  call random_number(odx)
  call random_number(ody)
  call random_number(odz)
  call random_number(ocs)
  call random_number(csv)
  allocate(w_cpu(size(w)), w_gpu(size(w)))
  ct = 0.
  do i = 1, num_steps
    call cpu_time(time1)
    call solve_continuity_equation
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do
  w_cpu = reshape(w, [size(w)])

  !$acc data async(async_id) copyout(w(1:nx,1:ny,1:nzf)) copyin(ny, nx, nz, nzf) &
  !$acc& copyin(u(0:nxf+1,0:ny+1,1:nz), v(0:nx+1,0:nyf+1,1:nz), iw(0:nx+1,0:ny+1,1:nzf)) &
  !$acc& copyin(odx(1:ny), ody(1:ny), odz(1:nz), ocs(1:ny), csv(1:nyf), kb(0:nx+1,0:ny+1))
  gt = 0.
  do i = 1, num_steps
    call cpu_time(time1)
    call solve_continuity_equation_gpu(async_id)
  !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc end data
  !$acc wait(async_id)
  w_gpu = reshape(w, [size(w)])

  call assert_allclose(w_gpu, size(w_gpu), w_cpu, size(w_cpu), 1e-10, 1e-10, "Array w")

   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt

end subroutine solve_continuity_equation_unit
