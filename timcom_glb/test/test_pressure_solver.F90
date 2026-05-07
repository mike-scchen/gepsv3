!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_pressure_solver
  implicit none

  call timcom_initialize_gpu
  call pressure_solver_unit
  call timcom_finalize_gpu

end program test_pressure_solver

subroutine pressure_solver_unit
  use timcom_const, only: r8
  use timcom_drv, only: pressure_solver
  use timcom_drv_gpu, only: pressure_solver_gpu, s_drv => s
  use timcom_solver_gpu, only: rhn, rho, alpha, beta, w_solver_gpu => w, tmp, tp, res1, res2, res3
  use timcom_solver_gpu, only: pcsi_gpu, a0r, s, r, q, csalpha, csbeta, csy, csomega, one_csy, rr
  use hyperlink, only: u, v, p0, iw, odx, ody, odz, ocs, csv, kb, w, x, iu, iv, in
  use hyperlink, only: ab, al, ac, ar, at, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh
  use hyperlink, only: nx, ny, nxf, nyf, nz, max_iter_p0, odyv, u_change, v_change, xu_bgn, xu_end, yv_bgn, yv_end, odt, nzf
  implicit none

  real(r8), dimension(:, :, :), allocatable :: &
    u_origin, u_cpu, u_gpu, v_origin, v_cpu, v_gpu
  real(r8), dimension(:, :), allocatable :: p0_origin, p0_cpu, p0_gpu, &
    u_change_cpu, u_change_gpu, v_change_cpu, v_change_gpu
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  allocate( &
    u_origin(0:nxf + 1, 0:ny + 1, nz), v_origin(0:nx + 1, 0:nyf + 1, nz), &
    u_cpu(0:nxf + 1, 0:ny + 1, nz), v_cpu(0:nx + 1, 0:nyf + 1, nz), &
    u_gpu(0:nxf + 1, 0:ny + 1, nz), v_gpu(0:nx + 1, 0:nyf + 1, nz), &
    p0_cpu(0:nx + 1, 0:ny + 1), p0_gpu(0:nx + 1, 0:ny + 1), &
    p0_origin(0:nx + 1, 0:ny + 1))
  allocate( &
    u_change_cpu(0:nxf + 1, 0:ny + 1), u_change_gpu(0:nxf + 1, 0:ny + 1), &
    v_change_cpu(0:nx + 1, 0:nyf + 1), v_change_gpu(0:nx + 1, 0:nyf + 1))
  allocate(a0r(0:nx+1,0:ny+1), s(0:nx+1,0:ny+1), q(0:nx+1,0:ny+1), r(0:nx+1,0:ny+1,2))
  allocate(s_drv(nx,ny))

  call random_number(u_origin)
  call random_number(v_origin)
  call random_number(p0_origin)
  max_iter_p0 = 10

  do i = 1, num_steps
    u = u_origin
    v = v_origin
    p0 = p0_origin
    call pressure_solver(1, 1)
  end do
  u_cpu = u
  v_cpu = v
  p0_cpu = p0
  u_change_cpu = u_change
  v_change_cpu = v_change

  !$acc enter data copyin(u, v, iw, odx, ody, odz, ocs, csv, kb, w, x, &
  !$acc& ab, al, ac, ar, at, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh, &
  !$acc& p0, iu, iv, odyv, u_change, v_change, nx, nxf, ny, nyf, nz, xu_bgn, &
  !$acc& xu_end, yv_bgn, yv_end, odt, nzf, in, u_origin, v_origin, p0_origin) async(async_id)
  !$acc enter data create(s_drv, rhn, rho, alpha, beta, w_solver_gpu, tmp, tp, res1, res2, res3, &
  !$acc& a0r, s, r, q, csalpha, csbeta, csy, csomega, one_csy, rr) async(async_id)
  do i = 1, num_steps
    !$acc kernels async(async_id)
    u = u_origin
    v = v_origin
    p0 = p0_origin
    !$acc end kernels
    call pressure_solver_gpu(1, 1, async_id)
  end do
  !$acc exit data copyout(u, v, iw, odx, ody, odz, ocs, csv, kb, w, x, &
  !$acc& ab, al, ac, ar, at, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh, &
  !$acc& p0, iu, iv, odyv, u_change, v_change, nx, nxf, ny, nyf, nz, &
  !$acc& xu_bgn, xu_end, yv_bgn, yv_end, odt, nzf, in, u_origin, v_origin, p0_origin) async(async_id)
  !$acc exit data delete(s_drv, rhn, rho, alpha, beta, w_solver_gpu, tmp, tp, res1, res2, res3, &
  !$acc& a0r, s, r, q, csalpha, csbeta, csy, csomega, one_csy, rr) async(async_id)
  !$acc wait(async_id)
  u_gpu = u
  v_gpu = v
  p0_gpu = p0
  u_change_gpu = u_change
  v_change_gpu = v_change

  call assert_allclose(u_gpu, size(u_gpu), u_cpu, size(u_cpu), 1e-10, 1e-10, "Array u")
  call assert_allclose(v_gpu, size(v_gpu), v_cpu, size(v_cpu), 1e-10, 1e-10, "Array v")
  call assert_allclose(p0_gpu(1:nx, 1:ny), nx*ny, p0_cpu(1:nx, 1:ny), nx*ny, 1e-10, 1e-10, "Array p0")
  call assert_allclose(u_change_gpu, size(u_change_gpu), u_change_cpu, size(u_change_cpu), 1e-10, 1e-10, "Array u_change")
  call assert_allclose(v_change_gpu, size(v_change_gpu), v_change_cpu, size(v_change_cpu), 1e-10, 1e-10, "Array v_change")

end subroutine pressure_solver_unit

