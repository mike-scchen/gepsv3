module cpl_smat
  use m_SparseMatrix,     only: SparseMatrix
  use m_SparseMatrix,     only: SparseMatrix_init => init
  use m_SparseMatrix,     only: SparseMatrix_importGRowInd => importGlobalRowIndices
  use m_SparseMatrix,     only: SparseMatrix_importGColInd => importGlobalColumnIndices
  use m_SparseMatrix,     only: SparseMatrix_importMatElts => importMatrixElements
  use m_SparseMatrixPlus, only: SparseMatrixPlus
  use m_SparseMatrixPlus, only: SparseMatrixPlus_init  => init
  use m_SparseMatrixPlus, only: SparseMatrixPlus_clean => clean
  use m_SparseMatrixPlus, only: Xonly
  use cpl_rank
  implicit none

  type(SparseMatrix), target ::   gfs2rsm_sMat
  type(SparseMatrix), target ::  gfs2gocn_sMat
  type(SparseMatrix), target ::  gocn2gfs_sMat
  type(SparseMatrix), target ::  rsm2rocn_sMat
  type(SparseMatrix), target ::  rocn2rsm_sMat
  type(SparseMatrix), target :: gocn2rocn_sMat

  type(SparseMatrixPlus), target ::   gfs2rsm_sMatP
  type(SparseMatrixPlus), target ::  gfs2gocn_sMatP
  type(SparseMatrixPlus), target ::  gocn2gfs_sMatP
  type(SparseMatrixPlus), target ::  rsm2rocn_sMatP
  type(SparseMatrixPlus), target ::  rocn2rsm_sMatP
  type(SparseMatrixPlus), target :: gocn2rocn_sMatP

contains
subroutine cpl_smat_init(myrank, rank_root, id_src, id_dst, mpi_comm_mct, SCRIP_file)
  use m_GlobalSegMap, only: GlobalSegMap
  use cpl_gsmap 
  implicit none

  integer, intent(in) :: myrank, rank_root
  integer, intent(in) :: id_src, id_dst
  integer, intent(in) :: mpi_comm_mct
  character(len=*), intent(in) :: SCRIP_file

  type(GlobalSegMap), pointer :: srcGSMap, dstGSMap
  type(SparseMatrix), pointer :: remap_sMat
  type(SparseMatrixPlus), pointer :: remap_sMatPlus

  if(id_src .eq. id_gfs .and. id_dst .eq. id_gocn) then
    srcGSMap=>gfsGSMap
    dstGSMap=>gocn_in_gfsGSMap
    remap_sMat=>gfs2gocn_sMat
    remap_sMatPlus=>gfs2gocn_sMatP

  elseif(id_src .eq. id_gfs .and. id_dst .eq. id_rsm) then
    srcGSMap=>gfsGSMap
    dstGSMap=>rsm_in_gfsGSMap
    remap_sMat=>gfs2rsm_sMat
    remap_sMatPlus=>gfs2rsm_sMatP

  elseif(id_src .eq. id_gocn .and. id_dst .eq. id_gfs) then
    srcGSMap=>gocnGSMap
    dstGSMap=>gfs_in_gocnGSMap
    remap_sMat=>gocn2gfs_sMat
    remap_sMatPlus=>gocn2gfs_sMatP

  elseif(id_src .eq. id_gocn .and. id_dst .eq. id_rocn) then
    srcGSMap=>gocnGSMap
    dstGSMap=>rocn_in_gocnGSMap
    remap_sMat=>gocn2rocn_sMat
    remap_sMatPlus=>gocn2rocn_sMatP

  elseif(id_src .eq. id_rsm .and. id_dst .eq. id_rocn) then
    srcGSMap=>rsmGSMap
    dstGSMap=>rocn_in_rsmGSMap
    remap_sMat=>rsm2rocn_sMat
    remap_sMatPlus=>rsm2rocn_sMatP

  elseif(id_src .eq. id_rocn .and. id_dst .eq. id_rsm) then
    srcGSMap=>rocnGSMap
    dstGSMap=>rsm_in_rocnGSMap
    remap_sMat=>rocn2rsm_sMat
    remap_sMatPlus=>rocn2rsm_sMatP

  end if
!  write(*,*) "yc check smat init:", myrank, id_src, trim(SCRIP_file)
  if(trim(SCRIP_file) .eq. "SameGrid") then
    if(myrank .eq. 0) write(*,*) "[CPL]: The Grid of ", id_src, " and ", id_dst, " is the SAME."
  else
    if(myrank .eq. 0) call read_SCRIP_file(SCRIP_file, remap_sMat, id_src)
    call SparseMatrixPlus_init(remap_sMatPlus, remap_sMat, srcGSMap, dstGSMap, Xonly, rank_root, mpi_comm_mct, id_src) 
  end if
end subroutine cpl_smat_init

subroutine read_SCRIP_file(SCRIP_file, remap_sMat, tag)
  use NETCDF
  implicit none
  
  character(len=*), intent(in) :: SCRIP_file
  type(SparseMatrix), intent(out) :: remap_sMat
  integer, intent(in) :: tag

  integer :: ierr, ncid, varid
  integer :: src_dims(2), dst_dims(2), dimIDs(1)
  integer :: num_elems, nRows, nCols
  real(8), pointer :: sparse_weights(:)
  integer, pointer :: sparse_rows(:), sparse_cols(:)

  ierr = NF90_OPEN(trim(SCRIP_file), NF90_NOWRITE, ncid)
 
  ierr = NF90_INQ_VARID(ncid, "src_grid_dims", varid)
  ierr = NF90_GET_VAR(ncid, varid, src_dims, (/1/), (/2/))

  ierr = NF90_INQ_VARID(ncid, "dst_grid_dims", varid)
  ierr = NF90_GET_VAR(ncid, varid, dst_dims, (/1/), (/2/))

  ierr = NF90_INQ_VARID(ncid, "S", varid)
  ierr = NF90_INQUIRE_VARIABLE(ncid, varid, dimids=dimIDs)
  ierr = NF90_INQUIRE_DIMENSION(ncid, dimIDs(1), len=num_elems)

  allocate( sparse_weights(num_elems) )
  allocate( sparse_rows(num_elems) )
  allocate( sparse_cols(num_elems) )

  ierr = NF90_GET_VAR(ncid, varid, sparse_weights, (/1/), (/num_elems/))
  ierr = NF90_INQ_VARID(ncid, "row", varid)
  ierr = NF90_GET_VAR(ncid, varid, sparse_rows, (/1/), (/num_elems/))
  ierr = NF90_INQ_VARID(ncid, "col", varid)
  ierr = NF90_GET_VAR(ncid, varid, sparse_cols, (/1/), (/num_elems/))  
  ierr = NF90_CLOSE(ncid)

  nRows = dst_dims(1)*dst_dims(2)
  nCols = src_dims(1)*src_dims(2)

  write(*,*) "[CPL]: in read sMat file:", nRows, nCols, tag, num_elems
  call SparseMatrix_init(remap_sMat, nRows, nCols, num_elems)
  write(*,*) "[CPL]: check sMat nRows: ", tag, size(sparse_rows)
  call SparseMatrix_importGRowInd(remap_sMat, sparse_rows,    size(sparse_rows))
  write(*,*) "[CPL]: check sMat nCols: ", tag, size(sparse_cols)
  call SparseMatrix_importGColInd(remap_sMat, sparse_cols,    size(sparse_cols))
  write(*,*) "[CPL]: check sMat nWgts: ", tag, size(sparse_weights)
  call SparseMatrix_importMatElts(remap_sMat, sparse_weights, size(sparse_weights))

  deallocate( sparse_weights )
  deallocate( sparse_rows )
  deallocate( sparse_cols )
end subroutine read_SCRIP_file

end module cpl_smat
