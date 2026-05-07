subroutine rsm_cpl_init(compid)
  use cpl_rank,  only: id_rocn
  use cpl_gsmap, only: cpl_gsmap_init
  use cpl_attr,  only: cpl_attr_init
  use cpl_smat,  only: cpl_smat_init
#include <paramodel.h>
#include <npesi.h>
#include <comrmpi.h>

  integer, intent(in) :: compid
  integer :: j, rank_root, seg_strt(jgrd1p-1), seg_leng(jgrd1p-1), nx, ny

  rank_root = 0
  nx = igrd1p-1
  ny = jgrd1p-1

  do j = 1, ny
    seg_strt(j) = (mycol*nx + (j-1)*igrd + 1) + myrow*ny*igrd
    seg_leng(j) = nx
  end do

  call cpl_gsmap_init(mype, ny, seg_strt, seg_leng, rank_root, mpi_comm_rsm, compid, "drsm GSMap:")
  call cpl_smat_init(mype, rank_root, compid, id_rocn, MPI_COMM_RSM, "SameGrid")
  call cpl_attr_init(compid, id_rocn, MPI_COMM_RSM)
end subroutine rsm_cpl_init
