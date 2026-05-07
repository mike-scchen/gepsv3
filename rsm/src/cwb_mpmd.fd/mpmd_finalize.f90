      subroutine mpmd_finalize
 
      use mpi
 
      call MPI_BARRIER(MPI_COMM_WORLD, IERR )
      call MPI_FINALIZE(IERR)
 
      return
      end
