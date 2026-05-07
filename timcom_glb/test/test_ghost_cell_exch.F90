!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_ghost_cell_exch
  implicit none

  call timcom_initialize_gpu
  call mpi_exch_r8type2d_unit
  call timcom_finalize_gpu

end program test_ghost_cell_exch

subroutine mpi_exch_r8type2d_unit
  use timcom_comm
  use timcom_comm_gpu
  use hyperlink, only: m_comm_cart, r8type2d, nbid, ndim, symm_np, nx, ny
  implicit none

  real(r8), dimension(0:nx+1,0:ny+1) :: x, x_cpu, x_gpu
  integer, parameter :: async_id = 1

  call random_number(x)

  x_cpu = x
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, x_cpu, symm_np)

  x_gpu = x
  !$acc enter data copyin(x_gpu) async(async_id)
  call mpi_exch_r8type2d_gpu(m_comm_cart, r8type2d, nbid, ndim, x_gpu, symm_np, async_id)
  !$acc exit data copyout(x_gpu) async(async_id)
  !$acc wait(async_id)

  call assert_allclose(x_gpu(1:nx, 1:ny), nx*ny, x_cpu(1:nx, 1:ny), nx*ny, 1e-10, 1e-10, "Array x")
  call assert_allclose(x_gpu(1:nx, 0), nx, x_cpu(1:nx, 0), nx, 1e-10, 1e-10, "Array x north boundary")
  call assert_allclose(x_gpu(1:nx, ny + 1), nx, x_cpu(1:nx, ny + 1), nx, 1e-10, 1e-10, "Array x south boundary")
  call assert_allclose(x_gpu(0, 1:ny), ny, x_cpu(0, 1:ny), ny, 1e-10, 1e-10, "Array x west boundary")
  call assert_allclose(x_gpu(nx + 1, 1:ny), ny, x_cpu(nx + 1, 1:ny), ny, 1e-10, 1e-10, "Array x east boundary")

end subroutine mpi_exch_r8type2d_unit
