module timcom_comm
  use timcom_const
  use MPI
  implicit none

  integer, parameter :: &
    rootid = 0

  integer :: &
    m_comm_timcom = MPI_COMM_WORLD,    & ! communication id of timcom model
    myid_timcom   = MPI_PROC_NULL,     & ! process id in the m_comm_timcom
    np_timcom     = -1,                & ! numbers of process in the m_comm_timcom
    fid_log       

  type :: comm_panel
    integer :: &
      m_comm_cart = MPI_COMM_NULL,     & ! communication id of the grid in timcom model 
      m_comm_x    = MPI_COMM_NULL,     & ! communication id in the x-dir of m_comm_cart
      m_comm_y    = MPI_COMM_NULL,     & ! communication id in the y-dir of m_comm_cart
      coords(3)   = MPI_COMM_NULL,     & ! coordinate in m_comm_cart; 1:z, 2:y, 3:x
      myid        = MPI_PROC_NULL,     & ! process id in the m_comm_cart
      myid_x      = MPI_PROC_NULL,     & ! process id in the x-dir of m_comm_cart
      myid_y      = MPI_PROC_NULL,     & ! process id in the y-dri of m_comm_cart
      nbid(5)     = MPI_PROC_NULL,     & ! neighbor process id; 1:w, 2:e, 3:s, 4:n, 5:np
      r8type2d(2) = MPI_DATATYPE_NULL, & ! dble datatype of MPI for 2d variable; 1:x, 2:y
      r8type3d(2) = MPI_DATATYPE_NULL, & ! dble datatype of MPI for 3d variable; 1:x, 2:y
      i2type2d(2) = MPI_DATATYPE_NULL, & ! int2 datatype of MPI for 2d variable; 1:x, 2:y
      i2type3d(2) = MPI_DATATYPE_NULL, & ! int2 datatype of MPI for 3d variable; 1:x, 2:y
      r8type3du(2)= MPI_DATATYPE_NULL, & ! dble datatype of MPI for 3d u variable; 1:x, 2:y
      r8type3dv(2)= MPI_DATATYPE_NULL, & ! dble datatype of MPI for 3d v variable; 1:x, 2:y
      np          = -1,                & ! numbers of process in the m_comm_timcom
      npx         = -1,                & ! numbers of process in the x-dir of m_comm_cart
      npy         = -1                   ! numbers of process in the y-dir of m_comm_cart
    integer(MPI_OFFSET_KIND) :: &
      ns_ijk(3) = -1, &
      nc_ijk(3) = -1, &
      nc_ujk(3) = -1, &
      nc_ivk(3) = -1, &
      nc_ijw(3) = -1
    logical :: &
      peri_x      = .false.,           & ! periodic boundary condition in the x-dir
      peri_y      = .false.,           & ! periodic boundary condition in the  y-dir
      symm_np     = .false.              ! symmetry boundary condition in North Pole
  end type
  
  interface ghost_cell_exch
     module procedure mpi_exch_r8type3d
     module procedure mpi_exch_r8type2d
     module procedure mpi_exch_i2type3d
     module procedure mpi_exch_i2type2d
  end interface

contains

#ifdef cpl
subroutine comm_init(comm_id)
#else
subroutine comm_init
#endif
  implicit none

#ifdef cpl
  integer, intent(in), optional :: comm_id
#endif
  integer ierr

#ifdef cpl
  if(present(comm_id)) then
    m_comm_timcom = comm_id
  else
    call MPI_INIT(ierr)
  end if
#else
  call MPI_INIT(ierr)
#endif

  call MPI_COMM_SIZE(m_comm_timcom,   np_timcom, ierr)
  call MPI_COMM_RANK(m_comm_timcom, myid_timcom, ierr)

#ifdef USE_GPU
  call ocn_device_init
  call ocn_nccl_init
#endif
  
  call MPI_BARRIER(m_comm_timcom, ierr)

  call MPI_FILE_OPEN(m_comm_timcom, 'runtime.log',  &
                     MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                     MPI_INFO_NULL, fid_log, ierr)

end subroutine comm_init

subroutine comm_cart_gen(comm)
  implicit none

  type(comm_panel), intent(inout) :: comm
  integer :: npdim(3), nbid_np(2) = MPI_PROC_NULL, ierr, n
  logical :: rdim(3) = .false., peri(3)
  character(len=256) :: log_info

  comm%np = comm%npx*comm%npy

  npdim(1) = 1
  npdim(2) = comm%npy
  npdim(3) = comm%npx

  peri(1) = .false.
  peri(2) = comm%peri_y
  peri(3) = comm%peri_x

! cartisean remark
! -------------------------------
! |   8   |   9   |   10  |   11  |
! | (0,2) | (1,2) | (2,2) | (3,2) |
! |-------|-------|-------|-------|
! |   4   |   5   |   6   |   7   |
! | (0,1) | (1,1) | (2,1) | (3,1) |
! |-------|-------|-------|-------|
! |   0   |   1   |   2   |   3   |
! | (0,0) | (1,0) | (2,0) | (3,0) |
! -------------------------------

! call MPI_DIMS_CREATE(comm%ng, 3, npdim, ierr) !
  call MPI_CART_CREATE(m_comm_timcom, 3, npdim, peri, .true., comm%m_comm_cart, ierr)
  call MPI_COMM_RANK(comm%m_comm_cart, comm%myid,  ierr)
  call MPI_CART_COORDS(comm%m_comm_cart, comm%myid, 3, comm%coords, ierr)
  call MPI_CART_SHIFT(comm%m_comm_cart, 2, 1, comm%nbid(1), comm%nbid(2), ierr)
  call MPI_CART_SHIFT(comm%m_comm_cart, 1, 1, comm%nbid(3), comm%nbid(4), ierr)

  rdim = .false.
  rdim(3) = .true.
  call MPI_CART_SUB(comm%m_comm_cart, rdim, comm%m_comm_x, ierr)
  call MPI_COMM_RANK(comm%m_comm_x, comm%myid_x, ierr)
  rdim = .false.
  rdim(2) = .true.
  call MPI_CART_SUB(comm%m_comm_cart, rdim, comm%m_comm_y, ierr)
  call MPI_COMM_RANK(comm%m_comm_y, comm%myid_y, ierr)

! North Pole Symmetry Boundary Condition
  if(comm%symm_np) then
    if(mod(comm%npx,2) .eq. 0) then
      call MPI_CART_SHIFT(comm%m_comm_cart, 2, comm%npx/2, nbid_np(1), nbid_np(2), ierr)
      if(comm%myid_y .eq. comm%npy-1) then
        comm%nbid(5) = nbid_np(1) ! nbid_np(1) == nbid_np(2) due to NP symmetry
      else
        comm%symm_np = .false.
        comm%nbid(5) = MPI_PROC_NULL
      end if
    else
      write(log_info,'(a)') '[ERROR]: can not impose symmetry BC in North Pole'
      call comm_write_log_info(fid_log, log_info)
      call comm_finalize
    end if
  end if

!  if(comm%myid .eq. rootid) then
!    write(*,'(a,i6)') "[INFO]: comm_cart_gen done.", MPI_OFFSET_KIND
!    write(*,'(10a8,3a10)') &
!          'myid', 'myid_x', 'myid_y', &
!          'nbid_w', 'nbid_e', 'nbid_s', 'nbid_n', &
!          'nbid_np', 'nbid_ns', 'nbid_nd', & 
!          'comm_cart', 'comm_x', 'comm_y'
!  end if
!  do n = 0, comm%np-1
!    if(comm%myid .eq. n) then
!      write(*,'(10i8,3i10)') &
!            comm%myid, comm%myid_x, comm%myid_y, &
!            comm%nbid(1), comm%nbid(2), comm%nbid(3), comm%nbid(4), &
!            comm%nbid(5), nbid_np(1), nbid_np(2), &
!            comm%m_comm_cart, comm%m_comm_x, comm%m_comm_y
!    end if
!    call MPI_BARRIER(comm%m_comm_cart, ierr)
!  end do
  call MPI_BARRIER(comm%m_comm_cart, ierr)
end subroutine comm_cart_gen

subroutine comm_def_mpi_type(comm, nx, ny, nz)
  implicit none

  type(comm_panel), intent(inout) :: comm
  integer, pointer :: nx, ny, nz
  integer :: ngh, nxgh, nygh, nxf, nyf, nzf

  nxf = nx + 1
  nyf = ny + 1
  nzf = nz + 1

  comm%ns_ijk(1) = comm%myid_x*nx + 1
  comm%ns_ijk(2) = comm%myid_y*ny + 1
  comm%ns_ijk(3) = 1

  comm%nc_ijk(1) = nx
  comm%nc_ijk(2) = ny
  comm%nc_ijk(3) = nz

  comm%nc_ujk = comm%nc_ijk
  comm%nc_ivk = comm%nc_ijk
  comm%nc_ijw = comm%nc_ijk

  comm%nc_ujk(1) = nxf
  comm%nc_ivk(2) = nyf
  comm%nc_ijw(3) = nzf

  !if(myid_timcom .eq. rootid) then
  !  write(*,*) '[INFO]: comm_def_mpi_type '
  !  write(*,*) nx, ny, nz
  !  write(*,*) comm%nc_ijk, comm%ns_ijk
  !end if

  ngh  = 2
  nxgh = nx+2*ngh
  nygh = ny+2*ngh

  call commit_mpi_type(nygh*nz, ngh, nxgh, MPI_REAL8,    comm%r8type3d(1))
  call commit_mpi_type(nygh*nz, ngh, nxgh, MPI_INTEGER2, comm%i2type3d(1))

  call commit_mpi_type( nz, ngh*nxgh, nxgh*nygh, MPI_REAL8,    comm%r8type3d(2))
  call commit_mpi_type( nz, ngh*nxgh, nxgh*nygh, MPI_INTEGER2, comm%i2type3d(2))

  ngh  = 1
  nxgh = nx+2*ngh
  nygh = ny+2*ngh

  call commit_mpi_type(nygh, ngh, nxgh, MPI_REAL8,    comm%r8type2d(1))
  call commit_mpi_type(nygh, ngh, nxgh, MPI_INTEGER2, comm%i2type2d(1))

  call commit_mpi_type(1, ngh*nxgh, 1, MPI_REAL8,    comm%r8type2d(2))
  call commit_mpi_type(1, ngh*nxgh, 1, MPI_INTEGER2, comm%i2type2d(2))

  ngh = 1
  nxgh = 2*ngh+nxf
  nygh = 2*ngh+ny
  call commit_mpi_type(nygh*nz, ngh,      nxgh,      MPI_REAL8, comm%r8type3du(1))
  call commit_mpi_type(     nz, ngh*nxgh, nxgh*nygh, MPI_REAL8, comm%r8type3du(2))

  nxgh = 2*ngh+nx
  nygh = 2*ngh+nyf
  call commit_mpi_type(nygh*nz, ngh,      nxgh,      MPI_REAL8, comm%r8type3dv(1))
  call commit_mpi_type(     nz, ngh*nxgh, nxgh*nygh, MPI_REAL8, comm%r8type3dv(2))
end subroutine comm_def_mpi_type

subroutine commit_mpi_type(nb, nelm, stride, oldtype, newtype)
  implicit none

  integer :: nb, nelm, stride
  integer :: oldtype, newtype
  integer :: ierr

  call MPI_TYPE_VECTOR(nb, nelm, stride, oldtype, newtype, ierr)
  call MPI_TYPE_COMMIT(newtype, ierr)

end subroutine commit_mpi_type

subroutine mpi_exch_r8type3d(comm_id, datatype, nbid, ndim, scr, symm_np)
  implicit none
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  real(r8), intent(inout) :: scr(:,:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh=2
  integer :: isend, irecv, ierr, istat(MPI_STATUS_SIZE)

  isend = ngh+1
  irecv = ndim(1)+isend

  call MPI_SENDRECV(scr( isend, 1, 1), 1, datatype(1), nbid(1), 1, &
                    scr( irecv, 1, 1), 1, datatype(1), nbid(2), 1, &
                    comm_id, istat, ierr)
  isend = ndim(1)+1
  irecv = 1
  call MPI_SENDRECV(scr(isend, 1, 1), 1, datatype(1), nbid(2), 1, &
                    scr(irecv, 1, 1), 1, datatype(1), nbid(1), 1, &
                    comm_id, istat, ierr)
  isend = ngh+1
  irecv = ndim(2)+isend
  call MPI_SENDRECV(scr( 1, isend, 1), 1, datatype(2), nbid(3), 1, &
                    scr( 1, irecv, 1), 1, datatype(2), nbid(4), 1, &
                    comm_id, istat, ierr)
  isend = ndim(2)+1
  irecv = 1
  call MPI_SENDRECV(scr( 1,isend, 1), 1, datatype(2), nbid(4), 1, &
                    scr( 1,irecv, 1), 1, datatype(2), nbid(3), 1, &
                    comm_id, istat, ierr)

  if(symm_np) then
    isend = ndim(2)+1
    irecv = isend + ngh
    call MPI_SENDRECV(scr( 1, isend, 1), 1, datatype(2), nbid(5), 1, &
                      scr( 1, irecv, 1), 1, datatype(2), nbid(5), 1, &
                      comm_id, istat, ierr)
    scr(:,irecv:irecv+1,:) = scr(:,irecv+1:irecv:-1,:)
  endif

end subroutine mpi_exch_r8type3d

subroutine mpi_exch_i2type3d(comm_id, datatype, nbid, ndim, scr, symm_np)
  implicit none
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  integer(2), intent(inout) :: scr(:,:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh=2
  integer :: isend, irecv, ierr, istat(MPI_STATUS_SIZE)

  isend = ngh+1
  irecv = ndim(1)+isend

  call MPI_SENDRECV(scr( isend, 1, 1), 1, datatype(1), nbid(1), 1, &
                    scr( irecv, 1, 1), 1, datatype(1), nbid(2), 1, &
                    comm_id, istat, ierr)
  isend = ndim(1)+1
  irecv = 1
  call MPI_SENDRECV(scr(isend, 1, 1), 1, datatype(1), nbid(2), 1, &
                    scr(irecv, 1, 1), 1, datatype(1), nbid(1), 1, &
                    comm_id, istat, ierr)
  isend = ngh+1
  irecv = ndim(2)+isend
  call MPI_SENDRECV(scr( 1, isend, 1), 1, datatype(2), nbid(3), 1, &
                    scr( 1, irecv, 1), 1, datatype(2), nbid(4), 1, &
                    comm_id, istat, ierr)
  isend = ndim(2)+1
  irecv = 1
  call MPI_SENDRECV(scr( 1,isend, 1), 1, datatype(2), nbid(4), 1, &
                    scr( 1,irecv, 1), 1, datatype(2), nbid(3), 1, &
                    comm_id, istat, ierr)
  if(symm_np) then
    isend = ndim(2)+1
    irecv = isend+ngh
    call MPI_SENDRECV(scr( 1, isend, 1), 1, datatype(2), nbid(5), 1, &
                      scr( 1, irecv, 1), 1, datatype(2), nbid(5), 1, &
                      comm_id, istat, ierr)
    scr(:,irecv:irecv+1,:) = scr(:,irecv+1:irecv:-1,:)
  endif

end subroutine mpi_exch_i2type3d

subroutine mpi_exch_r8type2d(comm_id, datatype, nbid, ndim, scr, symm_np)
  implicit none
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  real(r8), intent(inout) :: scr(:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh=1
  integer :: isend, irecv, ierr, istat(MPI_STATUS_SIZE)

  isend = ngh+1
  irecv = ndim(1)+isend

  call MPI_SENDRECV(scr( isend, 1), 1, datatype(1), nbid(1), 1, &
                    scr( irecv, 1), 1, datatype(1), nbid(2), 1, &
                    comm_id, istat, ierr)
  isend = ndim(1)+1
  irecv = 1
  call MPI_SENDRECV(scr(isend, 1), 1, datatype(1), nbid(2), 1, &
                    scr(irecv, 1), 1, datatype(1), nbid(1), 1, &
                    comm_id, istat, ierr)
  isend = ngh+1
  irecv = ndim(2)+isend
  call MPI_SENDRECV(scr( 1, isend), 1, datatype(2), nbid(3), 1, &
                    scr( 1, irecv), 1, datatype(2), nbid(4), 1, &
                    comm_id, istat, ierr)
  isend = ndim(2)+1
  irecv = 1
  call MPI_SENDRECV(scr( 1,isend), 1, datatype(2), nbid(4), 1, &
                    scr( 1,irecv), 1, datatype(2), nbid(3), 1, &
                    comm_id, istat, ierr)
  if(symm_np) then
    isend = ndim(2)+1
    irecv = isend+ngh
    call MPI_SENDRECV(scr( 1, isend), 1, datatype(2), nbid(5), 1, &
                      scr( 1, irecv), 1, datatype(2), nbid(5), 1, &
                      comm_id, istat, ierr)
  endif

end subroutine mpi_exch_r8type2d

subroutine mpi_exch_i2type2d(comm_id, datatype, nbid, ndim, scr, symm_np)
  implicit none
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  integer(2), intent(inout) :: scr(:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh=1
  integer :: isend, irecv, ierr, istat(MPI_STATUS_SIZE)

  isend = ngh+1
  irecv = ndim(1)+isend

  call MPI_SENDRECV(scr( isend, 1), 1, datatype(1), nbid(1), 1, &
                    scr( irecv, 1), 1, datatype(1), nbid(2), 1, &
                    comm_id, istat, ierr)
  isend = ndim(1)+1
  irecv = 1
  call MPI_SENDRECV(scr(isend, 1), 1, datatype(1), nbid(2), 1, &
                    scr(irecv, 1), 1, datatype(1), nbid(1), 1, &
                    comm_id, istat, ierr)
  isend = ngh+1
  irecv = ndim(2)+isend
  call MPI_SENDRECV(scr( 1, isend), 1, datatype(2), nbid(3), 1, &
                    scr( 1, irecv), 1, datatype(2), nbid(4), 1, &
                    comm_id, istat, ierr)
  isend = ndim(2)+1
  irecv = 1
  call MPI_SENDRECV(scr( 1,isend), 1, datatype(2), nbid(4), 1, &
                    scr( 1,irecv), 1, datatype(2), nbid(3), 1, &
                    comm_id, istat, ierr)
  if(symm_np) then
    isend = ndim(2)+1
    irecv = isend+ngh
    call MPI_SENDRECV(scr( 1, isend), 1, datatype(2), nbid(5), 1, &
                      scr( 1, irecv), 1, datatype(2), nbid(5), 1, &
                      comm_id, istat, ierr)
  endif

end subroutine mpi_exch_i2type2d

subroutine mpi_exch_r8type3du(comm_id, datatype, nbid, ndim, scr, symm_np)
  implicit none
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  real(r8), intent(inout) :: scr(:,:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh=1
  integer :: isend, irecv, ierr, istat(MPI_STATUS_SIZE)

  isend = ngh+2
  irecv = ndim(1)+isend

  call MPI_SENDRECV(scr( isend, 1, 1), 1, datatype(1), nbid(1), 1, &
                    scr( irecv, 1, 1), 1, datatype(1), nbid(2), 1, &
                    comm_id, istat, ierr)
  isend = ndim(1)+ngh
  irecv = 1
  call MPI_SENDRECV(scr(isend, 1, 1), 1, datatype(1), nbid(2), 1, &
                    scr(irecv, 1, 1), 1, datatype(1), nbid(1), 1, &
                    comm_id, istat, ierr)
  isend = ngh+1
  irecv = ndim(2)+isend
  call MPI_SENDRECV(scr( 1, isend, 1), 1, datatype(2), nbid(3), 1, &
                    scr( 1, irecv, 1), 1, datatype(2), nbid(4), 1, &
                    comm_id, istat, ierr)
  isend = ndim(2)+ngh
  irecv = 1
  call MPI_SENDRECV(scr( 1,isend, 1), 1, datatype(2), nbid(4), 1, &
                    scr( 1,irecv, 1), 1, datatype(2), nbid(3), 1, &
                    comm_id, istat, ierr)
  if(symm_np) then
    isend = ndim(2)+ngh
    irecv = isend + ngh
    call MPI_SENDRECV(scr( 1, isend, 1), 1, datatype(2), nbid(5), 1, &
                      scr( 1, irecv, 1), 1, datatype(2), nbid(5), 1, &
                      comm_id, istat, ierr)
  endif

end subroutine mpi_exch_r8type3du

subroutine mpi_exch_r8type3dv(comm_id, datatype, nbid, ndim, scr, symm_np)
  implicit none
  integer, intent(in) :: comm_id, datatype(2), nbid(5), ndim(3)
  real(r8), intent(inout) :: scr(:,:,:)
  logical, intent(in) :: symm_np
  integer, parameter :: ngh=1
  integer :: isend, irecv, ierr, istat(MPI_STATUS_SIZE)

  isend = ngh+1
  irecv = ndim(1)+isend

  call MPI_SENDRECV(scr( isend, 1, 1), 1, datatype(1), nbid(1), 1, &
                    scr( irecv, 1, 1), 1, datatype(1), nbid(2), 1, &
                    comm_id, istat, ierr)
  isend = ndim(1)+ngh
  irecv = 1
  call MPI_SENDRECV(scr(isend, 1, 1), 1, datatype(1), nbid(2), 1, &
                    scr(irecv, 1, 1), 1, datatype(1), nbid(1), 1, &
                    comm_id, istat, ierr)
  isend = ngh+2
  irecv = ndim(2)+isend
  call MPI_SENDRECV(scr( 1, isend, 1), 1, datatype(2), nbid(3), 1, &
                    scr( 1, irecv, 1), 1, datatype(2), nbid(4), 1, &
                    comm_id, istat, ierr)
  isend = ndim(2)+ngh
  irecv = 1
  call MPI_SENDRECV(scr( 1,isend, 1), 1, datatype(2), nbid(4), 1, &
                    scr( 1,irecv, 1), 1, datatype(2), nbid(3), 1, &
                    comm_id, istat, ierr)
  if(symm_np) then
    isend = ndim(2)+ngh
    irecv = isend+2
    call MPI_SENDRECV(scr( 1, isend, 1), 1, datatype(2), nbid(5), 1, &
                      scr( 1, irecv, 1), 1, datatype(2), nbid(5), 1, &
                      comm_id, istat, ierr)
  endif

end subroutine mpi_exch_r8type3dv

subroutine comm_write_log_info(fid, log_info)
  implicit none
  
  integer, intent(in) :: fid
  character(len=*) :: log_info

  integer :: leng, istat(MPI_STATUS_SIZE), ierr

  log_info = trim(log_info)//achar(10)
  leng = len(trim(log_info))

  call MPI_FILE_WRITE_SHARED(fid, log_info, leng, MPI_CHAR, istat, ierr)
end subroutine comm_write_log_info

subroutine comm_finalize
  implicit none

  integer :: ierr

  call MPI_FILE_CLOSE(fid_log, ierr)
  call MPI_FINALIZE(ierr)
  stop
end subroutine comm_finalize

subroutine comm_finalize_cpl
  implicit none

  integer :: ierr

  call MPI_FILE_CLOSE(fid_log, ierr)

end subroutine comm_finalize_cpl

end module timcom_comm
