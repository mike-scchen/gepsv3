      subroutine mpmd_irecv(dt,len,isrc,itag,type)

      dimension dt(*)
      integer   len,isrc,itag
      character*1 type

      if(type .eq. 'I')then
        call mpmd_irecv_i(dt,len,isrc,itag)
      elseif(type .eq. 'R')then
        call mpmd_irecv_r(dt,len,isrc,itag)
      else
         print *,"unknown type ---> ",type 
         stop
      endif

      return
      end

!---------------------------------------------------

      subroutine mpmd_irecv_i(idata,len,isrc,itag)

      use mpi

      integer idata(*),len,isrc,itag,istatus(MPI_STATUS_SIZE)

      call MPI_IRECV( idata, len, MPI_INTEGER, isrc, itag, &
                       MPI_COMM_WORLD, irecvreq, IERR )

      CALL MPI_WAIT(irecvreq, istatus, ierr)

      if(ierr .ne. 0)then
         print *,"mpmd_irecv_i fail, error code=", IERR
         stop
      endif
 
      return
      end

!---------------------------------------------------

      subroutine mpmd_irecv_r(rdata,len,rsrc,rtag)

      use mpi

      integer len,rsrc,rtag,istatus(MPI_STATUS_SIZE)
      real*8  rdata(*)

      call MPI_IRECV( rdata, len, MPI_DOUBLE_PRECISION, rsrc, rtag, &
                       MPI_COMM_WORLD, irecvreq, IERR )

      CALL MPI_WAIT(irecvreq, istatus, ierr)

      if(ierr .ne. 0)then
         print *,"mpmd_irecv_r fail, error code=", IERR
         stop
      endif
 
      return
      end
