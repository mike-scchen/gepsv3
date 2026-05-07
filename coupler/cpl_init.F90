subroutine cpl_init(n_comps, mct_comm_world, id_comp, ngrid_comp)
  use m_MCTWorld, only: MCTWorld_init => init
  use cpl_rank, only:myrank_root_all, ngrid_comp_all
  use MPI
  implicit none

  integer, intent(in) :: n_comps, mct_comm_world, id_comp
  integer, intent(in), optional :: ngrid_comp
  integer :: ierr

  call MCTWorld_init(n_comps, MPI_COMM_WORLD, mct_comm_world, id_comp)
  !ngrid_comp_all(id_comp) = ngrid_comp 
  !call MPI_BCAST(ngrid_comp_all(id_comp), 1, MPI_INT, myrank_root_all, MPI_COMM_WORLD, ierr)
  !itag = 1
  !call MPI_SEND( ngrid_comp_all(id_comp), 1, MPI_INTEGER, rank_root_all(id_gfs), 1, MPI_COMM_WORLD, ierr ) 
  !ngrid_comp_all(id_comp) = ngrid_comp
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !write(*,'(a16,5i7)') 'myrank_root_all:',myrank_root_all, ngrid_comp_all
end subroutine cpl_init


