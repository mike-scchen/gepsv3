      subroutine mpmd_recv(dt,len,isrc,itag,type)

      dimension dt(*)
      integer   len,isrc,itag
      character*1 type

      if(type .eq. 'I')then
        call mpmd_recv_i(dt,len,isrc,itag)
      elseif(type .eq. 'R')then
        call mpmd_recv_r(dt,len,isrc,itag)
      else
         print *,"unknown type ---> ",type 
         stop
      endif

      return
      end

!---------------------------------------------------

      subroutine mpmd_recv_i(idata,len,isrc,itag)

      use mpi

      integer idata(*),len,isrc,itag,istatus(MPI_STATUS_SIZE)

      call MPI_RECV( idata, len, MPI_INTEGER, isrc, itag, &
                       MPI_COMM_WORLD, istatus, IERR )

      if(ierr .ne. 0)then
         print *,"mpmd_recv_i fail, error code=", IERR
         stop
      endif
 
      return
      end

!---------------------------------------------------

      subroutine mpmd_recv_r(rdata,len,rsrc,rtag)

      use mpi

      integer len,rsrc,rtag,istatus(MPI_STATUS_SIZE)
      real*8  rdata(*)

      call MPI_RECV( rdata, len, MPI_DOUBLE_PRECISION, rsrc, rtag, &
                       MPI_COMM_WORLD, istatus, IERR )

      if(ierr .ne. 0)then
         print *,"mpmd_recv_r fail, error code=", IERR
         stop
      endif
 
      return
      end
