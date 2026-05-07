      subroutine mpmd_barrier(icomm)

      use mpi

      integer icomm

      call MPI_BARRIER(icomm, IERR )

      if(ierr .ne. 0)then
         print *,"mpmd_barrier fail, error code=", IERR
         stop
      endif
 
      return
      end
