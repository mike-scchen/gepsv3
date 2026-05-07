!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_p_bicgstab_le
  implicit none

  call timcom_initialize_gpu
  call p_bicgstab_le_unit
  call timcom_finalize_gpu

end program test_p_bicgstab_le

subroutine p_bicgstab_le_unit
  use timcom_solver, only: p_bicgstab_le
  use timcom_solver_gpu, only: p_bicgstab_le_gpu, rhn, rho, alpha, beta, w, tmp, tp, res1, res2, res3
  use timcom_const, only: r8
  use hyperlink, only: ab, al, ac, ar, at, x, cb, cl, cc, cr, ct, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh
  use hyperlink, only: nx, ny, myid
  implicit none

  real(r8), dimension(:, :), allocatable :: x_cpu, x_gpu, s
  real(r8), dimension(1) :: res_p2, res_p3
  integer :: iter
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  allocate(x_cpu(0:nx+1,0:ny+1), x_gpu(0:nx+1,0:ny+1), s(nx, ny))

  call random_number(s)
  s = 3.4 * s - 1.5

  do i = 1, num_steps
    x_cpu = 0.d0

    call p_bicgstab_le( &
      ab, al, ac, ar, at, s, x_cpu, &
      cb, cl, cc, cr, ct,       &
      cgr, cgrh, cgp, cgv,      &
      cgs, cgt, cgph, cgsh, iter, res_p2(1), res_p3(1))
  end do
  print '(A, I2, I3, E10.3, E10.3)', "CPU", myid, iter, res_p2(1), res_p3(1)
  
  !$acc enter data copyin(ab, al, ac, ar, at, s, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh, nx, ny) create(x_gpu) async(async_id)
  !$acc enter data create(rhn, rho, alpha, beta, w, tmp, tp, res1, res2, res3) async(async_id)

  do i = 1, num_steps
    !$acc kernels async(async_id)
    x_gpu = 0.d0
    !$acc end kernels

    call p_bicgstab_le_gpu( &
      ab, al, ac, ar, at, s, x_gpu, &
      cb, cl, cc, cr, ct,       &
      cgr, cgrh, cgp, cgv,      &
      cgs, cgt, cgph, cgsh, iter, async_id)
  end do
  !$acc exit data copyout(x_gpu) delete(ab, al, ac, ar, at, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh, s, nx, ny) async(async_id)
  !$acc exit data delete(rhn, rho, alpha, beta, w, tmp, tp, res1, res2, res3) async(async_id)
  !$acc wait(async_id)
  print '(A, I2, I3, E10.3, E10.3)', "GPU", myid, iter, res2, res3
  
  call assert_allclose(x_gpu(1:nx, 1:ny), nx*ny, x_cpu(1:nx, 1:ny), nx*ny, 1e-10, 1e-10, "Array x")
  call assert_allclose(x_gpu(1:nx, 0), nx, x_cpu(1:nx, 0), nx, 1e-10, 1e-10, "Array x north boundary")
  call assert_allclose(x_gpu(1:nx, ny + 1), nx, x_cpu(1:nx, ny + 1), nx, 1e-10, 1e-10, "Array x south boundary")
  call assert_allclose(x_gpu(0, 1:ny), ny, x_cpu(0, 1:ny), ny, 1e-10, 1e-10, "Array x west boundary")
  call assert_allclose(x_gpu(nx + 1, 1:ny), ny, x_cpu(nx + 1, 1:ny), ny, 1e-10, 1e-10, "Array x east boundary")

end subroutine p_bicgstab_le_unit
