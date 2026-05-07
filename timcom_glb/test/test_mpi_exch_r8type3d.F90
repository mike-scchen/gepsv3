!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_mpi_exch_r8type3d
  implicit none

  call timcom_initialize_gpu
  call mpi_exch_r8type3d
  call timcom_finalize_gpu

end program test_mpi_exch_r8type3d

subroutine mpi_exch_r8type3d
  use timcom_comm
  use timcom_comm_gpu
  use hyperlink, only: m_comm_cart, r8type3d, nbid, ndim, symm_np, nx, ny, nz, myid
  implicit none

  real(r8), dimension(-1:nx+2,-1:ny+2,nz) :: u2, u2_cpu, u2_gpu
  integer, parameter :: num_steps = 16
  integer, parameter :: async_id = 1
  integer :: i
  integer :: seed_size
  integer, allocatable :: seed(:)

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed = 12345 + myid
  call random_seed(put=seed)

  call random_number(u2)

  do i = 1, num_steps
    u2_cpu = u2
    call mpi_exch_r8type3d(m_comm_cart, r8type3d, nbid, ndim, u2_cpu, symm_np)
  end do

  !$acc enter data copyin(u2_gpu, u2, nx, ny, nz) async(async_id)
  do i = 1, num_steps
    !$acc kernels async(async_id)
    u2_gpu = u2
    !$acc end kernels
    call mpi_exch_r8type3d_gpu(m_comm_cart, r8type3d, nbid, ndim, u2_gpu, symm_np, async_id)
  end do
  !$acc exit data copyout(u2_gpu) delete(u2, nx, ny, nz) async(async_id)
  !$acc wait(async_id)

  call assert_allclose_r8type3d(u2_gpu, u2_cpu, nx, ny, nz, 1e-10, 1e-10, "Array u2")


end subroutine mpi_exch_r8type3d
