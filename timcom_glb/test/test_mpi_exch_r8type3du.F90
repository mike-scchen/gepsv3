!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_mpi_exch_r8type3du
  implicit none

  call timcom_initialize_gpu
  call mpi_exch_r8type3du_unit
  call timcom_finalize_gpu

end program test_mpi_exch_r8type3du

subroutine mpi_exch_r8type3du_unit
  use timcom_comm
  use timcom_comm_gpu
  use hyperlink, only: m_comm_cart, r8type3du, nbid, ndim, symm_np, nxf, ny, nz
  implicit none

  real(r8), dimension(0:nxf+1,0:ny+1,nz) :: u, u_cpu, u_gpu
  integer, parameter :: num_steps = 16
  integer, parameter :: async_id = 1
  integer :: i

  call random_number(u)

  do i = 1, num_steps
    u_cpu = u
    call mpi_exch_r8type3du(m_comm_cart, r8type3du, nbid, ndim, u_cpu, symm_np)
  end do

  !$acc enter data copyin(u_gpu, u, nxf, ny, nz) async(async_id)
  do i = 1, num_steps
    !$acc kernels async(async_id)
    u_gpu = u
    !$acc end kernels
    call mpi_exch_r8type3du_gpu(m_comm_cart, r8type3du, nbid, ndim, u_gpu, symm_np, async_id)
  end do
  !$acc exit data copyout(u_gpu) delete(u, nxf, ny, nz) async(async_id)
  !$acc wait(async_id)

  call assert_allclose(u_gpu, size(u_gpu), u_cpu, size(u_cpu), 1e-10, 1e-10, "Array u")

end subroutine mpi_exch_r8type3du_unit
