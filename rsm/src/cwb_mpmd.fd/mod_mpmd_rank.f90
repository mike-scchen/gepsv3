      module mpmd_rank

      implicit none

      integer nsize_all, nsize_gfs, nsize_rsm,  &
              myrank_all,myrank_gfs,myrank_rsm, &
              root_gfs,root_rsm,                &
              MPI_COMM_gfs,MPI_COMM_rsm

      common/mpmd_comm/nsize_all, nsize_gfs, nsize_rsm, &
              myrank_all,myrank_gfs,myrank_rsm,         &
              root_gfs,root_rsm,                        &
              MPI_COMM_gfs,MPI_COMM_rsm   

      end module mpmd_rank
