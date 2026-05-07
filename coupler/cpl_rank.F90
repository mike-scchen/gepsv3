module cpl_rank
  implicit none

  integer, parameter :: n_comp_all = 4

  integer :: id_gfs  = 1
  integer :: id_gocn = 3
  integer :: id_rsm  = 2
  integer :: id_rocn = 4

  integer :: nproc_comp_all(n_comp_all) = 0
  integer :: ngrid_comp_all(n_comp_all) = 0
  logical :: l_comp_act(n_comp_all)     = .false.

  integer :: avsize(n_comp_all) = 0
  integer :: myrank_root_all = 0
  integer :: rank_root_all(n_comp_all) = 0
end module cpl_rank
