!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_mpi_exch_r8type3dv
  implicit none

  call timcom_initialize_gpu
  call mpi_exch_r8type3dv_unit
  call timcom_finalize_gpu

end program test_mpi_exch_r8type3dv

subroutine mpi_exch_r8type3dv_unit
  use timcom_comm
  use timcom_comm_gpu
  use hyperlink, only: m_comm_cart, r8type3dv, nbid, ndim, symm_np, nx, nyf, nz
  implicit none

  real(r8), dimension(0:nx+1,0:nyf+1,nz) :: v, v_cpu, v_gpu
  integer, parameter :: num_steps = 16
  integer, parameter :: async_id = 1
  integer :: i

  call random_number(v)

  do i = 1, num_steps
    v_cpu = v
    call mpi_exch_r8type3dv(m_comm_cart, r8type3dv, nbid, ndim, v_cpu, symm_np)
  end do

  !$acc enter data copyin(v_gpu, v, nx, nyf, nz) async(async_id)
  do i = 1, num_steps
    !$acc kernels async(async_id)
    v_gpu = v
    !$acc end kernels
    call mpi_exch_r8type3dv_gpu(m_comm_cart, r8type3dv, nbid, ndim, v_gpu, symm_np, async_id)
  end do
  !$acc exit data copyout(v_gpu) delete(v, nx, nyf, nz) async(async_id)
  !$acc wait(async_id)

  call assert_allclose(v_gpu, size(v_gpu), v_cpu, size(v_cpu), 1e-10, 1e-10, "Array v")

end subroutine mpi_exch_r8type3dv_unit
