!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_pcsi
  implicit none

  call timcom_initialize_gpu
  call pcsi_unit
  call timcom_finalize_gpu

end program test_pcsi

subroutine pcsi_unit
  use timcom_solver, only: PCSI
  use timcom_solver_gpu, only: pcsi_gpu, a0r, s, r, q, csalpha, csbeta, csy, csomega, one_csy, rr
  use timcom_const, only: r8
  use hyperlink, only: ab, al, ac, ar, at, x, cb, cl, cc, cr, ct, in
  use hyperlink, only: nx, ny, myid
  implicit none

  real(r8), dimension(:, :), allocatable :: x_cpu, x_gpu, b
  real(r8) :: res_cpu, res_gpu
  integer :: iter
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  allocate(x_cpu(0:nx+1,0:ny+1), x_gpu(0:nx+1,0:ny+1), b(nx, ny))
  allocate(a0r(0:nx+1,0:ny+1), s(0:nx+1,0:ny+1), q(0:nx+1,0:ny+1))
  allocate(r(0:nx+1,0:ny+1,2))

  ! Generate random test data
  call random_number(b)
  b = 3.4 * b - 1.5

  ! CPU
  do i = 1, num_steps
    x_cpu = 0.d0
    call PCSI(ab, al, ac, ar, at, cb, cl, cc, cr, ct, b, x_cpu, iter, res_cpu)
  end do
  print '(A, I2, I3, E10.3)', "CPU", myid, iter, res_cpu

  ! GPU
  !$acc enter data async(async_id) &
  !$acc& copyin(ac, ab, al, ar, at, b, in, nx, ny) &
  !$acc& create(x_gpu, a0r, s, r, q, csalpha, csbeta, csy, csomega, one_csy, rr)
  do i = 1, num_steps
    !$acc kernels async(async_id)
    x_gpu = 0.d0
    !$acc end kernels
    call pcsi_gpu(ab, al, ac, ar, at, b, x_gpu, iter, res_gpu, async_id)
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(x_gpu) &
  !$acc& delete(ac, ab, al, ar, at, b, in, nx, ny, a0r, s, r, q, csalpha, csbeta, csy, csomega, one_csy, rr)
  !$acc wait(async_id)

  print '(A, I2, I3, E10.3)', "GPU", myid, iter, res_gpu

  ! Compare CPU and GPU results
  call assert_allclose(x_gpu(1:nx, 1:ny), nx*ny, x_cpu(1:nx, 1:ny), nx*ny, 1e-10, 1e-10, "Array x")
  call assert_allclose(x_gpu(1:nx, 0), nx, x_cpu(1:nx, 0), nx, 1e-10, 1e-10, "Array x north boundary")
  call assert_allclose(x_gpu(1:nx, ny + 1), nx, x_cpu(1:nx, ny + 1), nx, 1e-10, 1e-10, "Array x south boundary")
  call assert_allclose(x_gpu(0, 1:ny), ny, x_cpu(0, 1:ny), ny, 1e-10, 1e-10, "Array x west boundary")
  call assert_allclose(x_gpu(nx + 1, 1:ny), ny, x_cpu(nx + 1, 1:ny), ny, 1e-10, 1e-10, "Array x east boundary")

end subroutine pcsi_unit
