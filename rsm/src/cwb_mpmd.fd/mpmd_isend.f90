      subroutine mpmd_isend(dt,len,idst,itag,type)

      dimension dt(*)
      integer len,idst,itag
      character*1 type

      if(type .eq. 'I')then
         call mpmd_isend_i(dt,len,idst,itag)
      elseif(type .eq. 'R')then
         call mpmd_isend_r(dt,len,idst,itag)
      else
         print *,"unknown type ---> ",type 
         stop
      endif

      return
      end

!---------------------------------------------------

      subroutine mpmd_isend_i(idata,len,idst,itag)

      use mpi

      integer idata(*),len,idst,itag

      call MPI_ISEND( idata, len, MPI_INTEGER, idst, itag, &
                       MPI_COMM_WORLD, isendreq, IERR )

      if(ierr .ne. 0)then
         print *,"mpmd_isend_i fail, error code=", IERR
         stop
      endif
 
      return
      end

!-------------------------------------------------

      subroutine mpmd_isend_r(rdata,len,rdst,rtag)

      use mpi

      integer len,rdst,rtag
      real*8  rdata(*)

      call MPI_ISEND( rdata, len, MPI_DOUBLE_PRECISION, rdst, rtag, &
                       MPI_COMM_WORLD, isendreq, IERR )

      if(ierr .ne. 0)then
         print *,"mpmd_isend_r fail, error code=", IERR
         stop
      endif
 
      return
      end
