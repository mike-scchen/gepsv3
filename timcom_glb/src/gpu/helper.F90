!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define NCCLCHECK(ierr) call ocn_nccl_check_helper(ierr, __FILE__, __LINE__)

subroutine ocn_device_init
  use openacc
  use timcom_comm, only: myid_timcom
  implicit none

  integer :: num_device, device_id

  num_device = acc_get_num_devices(acc_device_nvidia)
  device_id = mod(myid_timcom, num_device)
  call acc_set_device_num(device_id, acc_device_nvidia)

end subroutine ocn_device_init

subroutine ocn_cuda_check_helper(ierr, filename, line)
  use cudafor
  implicit none

  integer, intent(in) :: ierr
  character(len=*), intent(in) :: filename
  integer, intent(in) :: line

  if (ierr .ne. cudaSuccess) then
    print '(5g0)', "Failed, CUDA error ", filename, ":", line
    print '(5g0)', cudaGetErrorString(ierr)
  end if

end subroutine ocn_cuda_check_helper

subroutine ocn_nccl_check_helper(ierr, filename, line)
  use nccl
  implicit none

  type(ncclResult), intent(in) :: ierr
  character(len=*), intent(in) :: filename
  integer, intent(in) :: line

  if (ierr .ne. ncclSuccess) then
    print '(5g0)', "Failed, NCCL error ", filename, ":", line
    print '(5g0)', ncclGetErrorString(ierr)
  end if

end subroutine ocn_nccl_check_helper

subroutine ocn_nccl_init
  use mpi
  use nccl
  use timcom_comm, only: myid_timcom, m_comm_timcom, rootid, np_timcom
  use timcom_comm_gpu, only: nccl_id, nccl_comm
  implicit none

  integer :: ierr

  if (myid_timcom .eq. rootid) then
    NCCLCHECK(ncclGetUniqueId(nccl_id))
  end if
  call MPI_Bcast(nccl_id%internal, size(nccl_id%internal), MPI_BYTE, 0, m_comm_timcom, ierr)
  NCCLCHECK(ncclCommInitRank(nccl_comm, np_timcom, nccl_id, myid_timcom))

end subroutine ocn_nccl_init
