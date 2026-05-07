!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define NCCLCHECK(ierr) call ocn_nccl_check_helper(ierr, __FILE__, __LINE__)

module timcom_comm_gpu
  use timcom_comm
  use nccl
  implicit none

  type(ncclUniqueId) :: nccl_id
  type(ncclComm) :: nccl_comm
  real(r8), allocatable, dimension(:) :: wsend_2d, wrecv_2d, esend_2d, erecv_2d

  interface ghost_cell_exch_gpu
    module procedure mpi_exch_r8type2d_gpu
    module procedure mpi_exch_r8type3d_gpu
  end interface

contains

subroutine timcom_comm_gpu_init
  use openacc
  use hyperlink, only: ny
  implicit none

  allocate(wsend_2d(ny), wrecv_2d(ny), esend_2d(ny), erecv_2d(ny))
  !$acc enter data create(wsend_2d, wrecv_2d, esend_2d, erecv_2d)

end subroutine timcom_comm_gpu_init

subroutine timcom_comm_gpu_finalize
  use openacc
  implicit none

  !$acc exit data delete(wsend_2d, wrecv_2d, esend_2d, erecv_2d)

end subroutine timcom_comm_gpu_finalize

subroutine mpi_exch_r8type2d_gpu(comm_id, datatype, nbid, ndim, scr, symm_np, async_id)
  ! Present on device
  !   Array: scr
  use cudafor
  use hyperlink, only: nx, ny
  use openacc
  use nccl

  implicit none
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  real(r8), intent(inout) :: scr(:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh = 1
  integer :: iw, ie, i, j, async_id
  integer(kind=cuda_stream_kind) :: stream
  logical :: has_s, has_n, has_np

  has_s = (nbid(3) .ge. 0) .and. (nbid(3) .lt. np_timcom)
  has_n = (nbid(4) .ge. 0) .and. (nbid(4) .lt. np_timcom)
  has_np = (nbid(5) .ge. 0) .and. (nbid(5) .lt. np_timcom)

  stream = acc_get_cuda_stream(async_id)

  iw = ngh + 1
  ie = ndim(1) + 1
  !$acc parallel loop async(async_id)
  do j = 1, ny
    wsend_2d(j) = scr(iw, j + 1)
    esend_2d(j) = scr(ie, j + 1)
  end do

  !$acc host_data use_device(scr, wsend_2d, wrecv_2d, esend_2d, erecv_2d)
  NCCLCHECK(ncclGroupStart())
  NCCLCHECK(ncclSend(wsend_2d(1), ny*ngh, ncclFloat64, nbid(1), nccl_comm, stream))
  NCCLCHECK(ncclRecv(wrecv_2d(1), ny*ngh, ncclFloat64, nbid(2), nccl_comm, stream))
  NCCLCHECK(ncclSend(esend_2d(1), ny*ngh, ncclFloat64, nbid(2), nccl_comm, stream))
  NCCLCHECK(ncclRecv(erecv_2d(1), ny*ngh, ncclFloat64, nbid(1), nccl_comm, stream))
  if (has_s) then
    NCCLCHECK(ncclSend(scr(2, ngh + 1), nx*ngh, ncclFloat64, nbid(3), nccl_comm, stream))
    NCCLCHECK(ncclRecv(scr(2, 1), nx*ngh, ncclFloat64, nbid(3), nccl_comm, stream))
  end if
  if (has_n) then
    NCCLCHECK(ncclRecv(scr(2, ndim(2) + ngh + 1), nx*ngh, ncclFloat64, nbid(4), nccl_comm, stream))
    NCCLCHECK(ncclSend(scr(2, ndim(2) + 1), nx*ngh, ncclFloat64, nbid(4), nccl_comm, stream))
  end if
  if (symm_np .and. has_np) then
    NCCLCHECK(ncclSend(scr(2, ndim(2) + 1), nx*ngh, ncclFloat64, nbid(5), nccl_comm, stream))
    NCCLCHECK(ncclRecv(scr(2, ndim(2) + ngh + 1), nx*ngh, ncclFloat64, nbid(5), nccl_comm, stream))
  end if
  NCCLCHECK(ncclGroupEnd())
  !$acc end host_data

  iw = ndim(1) + ngh + 1
  ie = 1
  !$acc parallel loop async(async_id)
  do j = 1, ny
    scr(iw, j + 1) = wrecv_2d(j)
    scr(ie, j + 1) = erecv_2d(j)
  end do
  
end subroutine mpi_exch_r8type2d_gpu

subroutine mpi_exch_r8type3d_gpu(comm_id, datatype, nbid, ndim, scr, symm_np, async_id)
  use cudafor
  use hyperlink, only: nxf, nx, ny, nz
  use nccl
  use openacc

  implicit none
  integer, intent(in) :: async_id
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  real(r8), intent(inout) :: scr(:,:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh=2
  integer :: ierr, istat(MPI_STATUS_SIZE), iw, ie, in, is, inp
  real(r8), dimension(ngh, ny + 4, nz) :: wsend, wrecv, esend, erecv
  real(r8), dimension(nx + 4, ngh, nz) :: nsend, nrecv, ssend, srecv, npsend, nprecv
  integer(kind=cuda_stream_kind) :: stream
  integer :: i, j, k, nxgh, nygh
  logical :: has_s, has_n, has_np

  stream = acc_get_cuda_stream(async_id)

  nxgh = nx+2*ngh
  nygh = ny+2*ngh
  has_s = (nbid(3) .ge. 0) .and. (nbid(3) .lt. np_timcom)
  has_n = (nbid(4) .ge. 0) .and. (nbid(4) .lt. np_timcom)
  has_np = (nbid(5) .ge. 0) .and. (nbid(5) .lt. np_timcom)

  !$acc enter data async(async_id) &
  !$acc& create(wsend, esend, wrecv, erecv) &
  !$acc& create(nsend, ssend, npsend, srecv, nrecv, nprecv)

  iw = ngh + 1
  ie = ndim(1) + 1
  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz
    do j = 1, ny + 4
      do i = 1, ngh
        wsend(i, j, k) = scr(iw + i - 1, j, k)
        esend(i, j, k) = scr(ie + i - 1, j, k)
      end do
    end do
  end do

  !$acc host_data use_device(wsend, wrecv, esend, erecv)
  NCCLCHECK(ncclGroupStart())
  NCCLCHECK(ncclSend(wsend(1, 1, 1), ngh*nygh*nz, ncclFloat64, nbid(1), nccl_comm, stream))
  NCCLCHECK(ncclSend(esend(1, 1, 1), ngh*nygh*nz, ncclFloat64, nbid(2), nccl_comm, stream))
  NCCLCHECK(ncclRecv(wrecv(1, 1, 1), ngh*nygh*nz, ncclFloat64, nbid(2), nccl_comm, stream))
  NCCLCHECK(ncclRecv(erecv(1, 1, 1), ngh*nygh*nz, ncclFloat64, nbid(1), nccl_comm, stream))
  NCCLCHECK(ncclGroupEnd())
  !$acc end host_data

  iw = ndim(1) + ngh + 1
  ie = 1
  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz
    do j = 1, ny + 4
      do i = 1, ngh
        scr(iw + i - 1, j, k) = wrecv(i, j, k)
        scr(ie + i - 1, j, k) = erecv(i, j, k)
      end do
    end do
  end do

  is = ngh + 1
  in = ndim(2) + 1
  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz
    do j = 1, ngh
      do i = 1, nx + 4
        if (has_s) then
          ssend(i, j, k) = scr(i, is + j - 1, k)
        end if
        if (has_n) then
          nsend(i, j, k) = scr(i, in + j - 1, k)
        end if
      end do
    end do
  end do

  !$acc host_data use_device(ssend, srecv, nsend, nrecv)
  NCCLCHECK(ncclGroupStart())
  if (has_s) then
    NCCLCHECK(ncclSend(ssend(1, 1, 1), nz*ngh*nxgh, ncclFloat64, nbid(3), nccl_comm, stream))
    NCCLCHECK(ncclRecv(nrecv(1, 1, 1), nz*ngh*nxgh, ncclFloat64, nbid(3), nccl_comm, stream))
  end if
  if (has_n) then
    NCCLCHECK(ncclSend(nsend(1, 1, 1), nz*ngh*nxgh, ncclFloat64, nbid(4), nccl_comm, stream))
    NCCLCHECK(ncclRecv(srecv(1, 1, 1), nz*ngh*nxgh, ncclFloat64, nbid(4), nccl_comm, stream))
  end if
  NCCLCHECK(ncclGroupEnd())
  !$acc end host_data

  is = ndim(2) + ngh + 1
  in = 1
  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz
    do j = 1, ngh
      do i = 1, nx + 4
        if (has_n) then
          scr(i, is + j - 1, k) = srecv(i, j, k)
        end if
        if (has_s) then
          scr(i, in + j - 1, k) = nrecv(i, j, k)
        end if
      end do
    end do
  end do

  if (symm_np) then
    inp = ndim(2) + 1
    !$acc parallel loop collapse(3) async(async_id)
    do k = 1, nz
      do j = 1, ngh
        do i = 1, nx + 4
          npsend(i, j, k) = scr(i, inp + j - 1, k)
        end do
      end do
    end do

    !$acc host_data use_device(npsend, nprecv)
    NCCLCHECK(ncclGroupStart())
    if (has_np) then
      NCCLCHECK(ncclSend(npsend(1, 1, 1), nz*ngh*nxgh, ncclFloat64, nbid(5), nccl_comm, stream))
      NCCLCHECK(ncclRecv(nprecv(1, 1, 1), nz*ngh*nxgh, ncclFloat64, nbid(5), nccl_comm, stream))
    end if
    NCCLCHECK(ncclGroupEnd())
    !$acc end host_data

    inp = ndim(2) + ngh + 1
    !$acc parallel loop collapse(2) async(async_id)
    do k = 1, nz
      do i = 1, nx + 4
        if (has_np) then
          scr(i, inp, k) = nprecv(i, 2, k)
          scr(i, inp + 1, k) = nprecv(i, 1, k)
        end if
      end do
    end do

  end if

  !$acc exit data async(async_id) &
  !$acc& delete(wsend, esend, wrecv, erecv) &
  !$acc& delete(nsend, ssend, npsend, srecv, nrecv, nprecv)

end subroutine mpi_exch_r8type3d_gpu

subroutine mpi_exch_r8type3du_gpu(comm_id, datatype, nbid, ndim, scr, symm_np, async_id)
  use cudafor
  use hyperlink, only: nxf, ny, nz
  use nccl
  use openacc

  implicit none
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  real(r8), intent(inout) :: scr(:,:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh=1
  integer :: ierr, istat(MPI_STATUS_SIZE), iw, ie, in, is, inp
  real(r8), dimension(ny + 2, nz) :: wsend, wrecv, esend, erecv
  real(r8), dimension(nxf + 2, nz) :: nsend, nrecv, ssend, srecv, npsend, nprecv
  integer(kind=cuda_stream_kind) :: stream
  integer :: async_id, i, j, k, nxgh, nygh
  logical :: has_s, has_n, has_np

  stream = acc_get_cuda_stream(async_id)

  nxgh = 2*ngh+nxf
  nygh = 2*ngh+ny
  has_s = (nbid(3) .ge. 0) .and. (nbid(3) .lt. np_timcom)
  has_n = (nbid(4) .ge. 0) .and. (nbid(4) .lt. np_timcom)
  has_np = (nbid(5) .ge. 0) .and. (nbid(5) .lt. np_timcom)

  !$acc enter data create(wsend, esend, wrecv, erecv) async(async_id)
  !$acc enter data create(nsend, ssend, npsend, srecv, nrecv, nprecv) async(async_id)

  iw = ngh + 2
  ie = ndim(1) + ngh
  !$acc parallel loop collapse(2) async(async_id)
  do k = 1, nz
    do j = 1, ny + 2
      wsend(j, k) = scr(iw, j, k)
      esend(j, k) = scr(ie, j, k)
    end do
  end do

  !$acc host_data use_device(wsend, wrecv, esend, erecv)
  NCCLCHECK(ncclGroupStart())
  NCCLCHECK(ncclSend(wsend(1, 1), nygh*nz, ncclFloat64, nbid(1), nccl_comm, stream))
  NCCLCHECK(ncclRecv(erecv(1, 1), nygh*nz, ncclFloat64, nbid(1), nccl_comm, stream))
  NCCLCHECK(ncclSend(esend(1, 1), nygh*nz, ncclFloat64, nbid(2), nccl_comm, stream))
  NCCLCHECK(ncclRecv(wrecv(1, 1), nygh*nz, ncclFloat64, nbid(2), nccl_comm, stream))
  NCCLCHECK(ncclGroupEnd())
  !$acc end host_data
  
  iw = ndim(1) + ngh + 2
  ie = 1
  !$acc parallel loop collapse(2) async(async_id)
  do k = 1, nz
    do j = 1, ny + 2
      scr(iw, j, k) = wrecv(j, k)
      scr(ie, j, k) = erecv(j, k)
    end do
  end do

  is = ngh + 1
  in = ndim(2) + ngh
  !$acc parallel loop collapse(2) async(async_id)
  do k = 1, nz
    do i = 1, nxf + 2
      if (has_s) then
        ssend(i, k) = scr(i, is, k)
      end if
      if (has_n) then
        nsend(i, k) = scr(i, in, k)
      end if
    end do
  end do

  !$acc host_data use_device(ssend, srecv, nsend, nrecv)
  NCCLCHECK(ncclGroupStart())
  if (has_s) then
    NCCLCHECK(ncclSend(ssend(1, 1), nz*ngh*nxgh, ncclFloat64, nbid(3), nccl_comm, stream))
    NCCLCHECK(ncclRecv(nrecv(1, 1), nz*ngh*nxgh, ncclFloat64, nbid(3), nccl_comm, stream))
  end if
  if (has_n) then
    NCCLCHECK(ncclSend(nsend(1, 1), nz*ngh*nxgh, ncclFloat64, nbid(4), nccl_comm, stream))
    NCCLCHECK(ncclRecv(srecv(1, 1), nz*ngh*nxgh, ncclFloat64, nbid(4), nccl_comm, stream))
  end if
  NCCLCHECK(ncclGroupEnd())
  !$acc end host_data

  is = ndim(2) + ngh + 1
  in = 1
  !$acc parallel loop collapse(2) async(async_id)
  do k = 1, nz
    do i = 1, nxf + 2
      if (has_n) then
        scr(i, is, k) = srecv(i, k)
      end if
      if (has_s) then
        scr(i, in, k) = nrecv(i, k)
      end if
    end do
  end do

  if (symm_np) then
    inp = ndim(2) + ngh
    !$acc parallel loop collapse(2) async(async_id)
    do k = 1, nz
      do i = 1, nxf + 2
        npsend(i, k) = scr(i, inp, k)
      end do
    end do
    !$acc host_data use_device(npsend, nprecv)
    NCCLCHECK(ncclGroupStart())
    if (has_np) then
      NCCLCHECK(ncclSend(npsend(1, 1), nz*ngh*nxgh, ncclFloat64, nbid(5), nccl_comm, stream))
      NCCLCHECK(ncclRecv(nprecv(1, 1), nz*ngh*nxgh, ncclFloat64, nbid(5), nccl_comm, stream))
    end if
    NCCLCHECK(ncclGroupEnd())
    !$acc end host_data
    inp = ndim(2) + 2*ngh
    !$acc parallel loop collapse(2) async(async_id)
    do k = 1, nz
      do i = 1, nxf + 2
        if (has_np) then
          scr(i, inp, k) = nprecv(i, k)
        end if
      end do
    end do
  end if

  !$acc exit data delete(wsend, esend, wrecv, erecv) async(async_id)
  !$acc exit data delete(nsend, ssend, npsend, srecv, nrecv, nprecv) async(async_id)

end subroutine mpi_exch_r8type3du_gpu

subroutine mpi_exch_r8type3dv_gpu(comm_id, datatype, nbid, ndim, scr, symm_np, async_id)
  use cudafor
  use hyperlink, only: nx, nyf, nz
  use nccl
  use openacc

  implicit none
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  real(r8), intent(inout) :: scr(:,:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh=1
  integer :: isend, irecv, ierr, istat(MPI_STATUS_SIZE), iw, ie, is, in, inp
  integer :: async_id, i, j, k, nxgh, nygh
  integer(kind=cuda_stream_kind) :: stream
  real(r8), dimension(nyf + 2, nz) :: wsend, wrecv, esend, erecv
  real(r8), dimension(nx + 2, nz) :: ssend, srecv, nsend, nrecv, npsend, nprecv
  logical :: has_s, has_n, has_np

  stream = acc_get_cuda_stream(async_id)

  nxgh = 2*ngh + nx
  nygh = 2*ngh + nyf
  has_s = (nbid(3) .ge. 0) .and. (nbid(3) .lt. np_timcom)
  has_n = (nbid(4) .ge. 0) .and. (nbid(4) .lt. np_timcom)
  has_np = (nbid(5) .ge. 0) .and. (nbid(5) .lt. np_timcom)

  !$acc enter data create(wsend, wrecv, esend, erecv) async(async_id)
  !$acc enter data create(ssend, srecv, nsend, nrecv, npsend, nprecv) async(async_id)

  iw = ngh + 1
  ie = ndim(1) + ngh
  !$acc parallel loop collapse(2) async(async_id)
  do k = 1, nz
    do j = 1, nyf + 2
      wsend(j, k) = scr(iw, j, k)
      esend(j, k) = scr(ie, j, k)
    end do
  end do

  !$acc host_data use_device(wsend, wrecv, esend, erecv)
  NCCLCHECK(ncclGroupStart())
  NCCLCHECK(ncclSend(wsend(1, 1), nygh*nz, ncclFloat64, nbid(1), nccl_comm, stream))
  NCCLCHECK(ncclRecv(wrecv(1, 1), nygh*nz, ncclFloat64, nbid(2), nccl_comm, stream))
  NCCLCHECK(ncclSend(esend(1, 1), nygh*nz, ncclFloat64, nbid(2), nccl_comm, stream))
  NCCLCHECK(ncclRecv(erecv(1, 1), nygh*nz, ncclFloat64, nbid(1), nccl_comm, stream))
  NCCLCHECK(ncclGroupEnd())
  !$acc end host_data

  iw = ndim(1) + ngh + 1
  ie = 1
  !$acc parallel loop collapse(2) async(async_id)
  do k = 1, nz
    do j = 1, nyf + 2
      scr(iw, j, k) = wrecv(j, k)
      scr(ie, j, k) = erecv(j, k)
    end do
  end do

  is = ngh + 2
  in = ndim(2) + ngh
  !$acc parallel loop collapse(2) async(async_id)
  do k = 1, nz
    do i = 1, nx + 2
      if (has_s) then
        ssend(i, k) = scr(i, is, k)
      end if
      if (has_n) then
        nsend(i, k) = scr(i, in, k)
      end if
    end do
  end do

  !$acc host_data use_device(ssend, srecv, nsend, nrecv)
  NCCLCHECK(ncclGroupStart())
  if (has_s) then
    NCCLCHECK(ncclSend(ssend(1, 1), nxgh*nz, ncclFloat64, nbid(3), nccl_comm, stream))
    NCCLCHECK(ncclRecv(nrecv(1, 1), nxgh*nz, ncclFloat64, nbid(3), nccl_comm, stream))
  end if
  if (has_n) then
    NCCLCHECK(ncclRecv(srecv(1, 1), nxgh*nz, ncclFloat64, nbid(4), nccl_comm, stream))
    NCCLCHECK(ncclSend(nsend(1, 1), nxgh*nz, ncclFloat64, nbid(4), nccl_comm, stream))
  end if
  NCCLCHECK(ncclGroupEnd())
  !$acc end host_data

  is = ndim(2) + ngh + 2
  in = 1
  !$acc parallel loop collapse(2) async(async_id)
  do k = 1, nz
    do i = 1, nx + 2
      if (has_n) then
        scr(i, is, k) = srecv(i, k)
      end if
      if (has_s) then
        scr(i, in, k) = nrecv(i, k)
      end if
    end do
  end do

  if (symm_np) then
    inp = ndim(2) + ngh
    !$acc parallel loop collapse(2) async(async_id)
    do k = 1, nz
      do i = 1, nx + 2
        if (has_np) then
          npsend(i, k) = scr(i, inp, k)
        end if
      end do
    end do

    !$acc host_data use_device(npsend, nprecv)
    NCCLCHECK(ncclGroupStart())
    if (has_np) then
      NCCLCHECK(ncclSend(npsend(1, 1), nxgh*nz, ncclFloat64, nbid(5), nccl_comm, stream))
      NCCLCHECK(ncclRecv(nprecv(1, 1), nxgh*nz, ncclFloat64, nbid(5), nccl_comm, stream))
    end if
    NCCLCHECK(ncclGroupEnd())
    !$acc end host_data
    inp = ndim(2) + ngh + 2
    !$acc parallel loop collapse(2) async(async_id)
    do k = 1, nz
      do i = 1, nx + 2
        if (has_np) then
          scr(i, inp, k) = nprecv(i, k)
        end if
      end do
    end do
  end if

  !$acc exit data delete(wsend, wrecv, esend, erecv) async(async_id)
  !$acc exit data delete(ssend, srecv, nsend, nrecv, npsend, nprecv) async(async_id)

end subroutine mpi_exch_r8type3dv_gpu

end module timcom_comm_gpu
