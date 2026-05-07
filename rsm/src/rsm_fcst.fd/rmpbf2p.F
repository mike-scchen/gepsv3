      subroutine rmpbf2p(a,lngrdfull,                                    &
     &                   b,lngrdpart,ntotal)
!$$$  subprogram documentation block
!
! subprogram: 	rmpbf2p
!            
! prgmmr: ying-ju chen  date:2021-07-19
!
! abstract: transpose (if,jf,kf) to (ip,jp,kf)
!
! program history log:
!    2021-07-19  ying-ju chen	build for basefield IO (for restart)
!
! usage:	call mpgf2p(a,igrd12,jgrd12,b,igrd12p,jgrd12p,ntotal)
!
!    input argument lists:
!	a	- real (igrd12,jgrd12,ntotal) full field 
!	igrd12	- integer total grid in longitude time 2
!	jgrd12	- integer total grid in latitude divide 2
!	ntotal	- integer total set of fields
!
!    output argument list:
!	b	- real (igrd12p,jgrd12p,ntotal) partial field for each pe 
!	igrd12p	- integer partial grid in longitude time 2
!	jgrd12p	- integer partial grid in latitude divide 2
! 
! subprograms called:
!   mpi_scatterv	- scatter message from master pe to all pe
!
! attributes:
!    library: mpi
!    language: fortran 90
!$$$
!
#include <define.h>
#ifdef MP
      implicit none
#include <npesi.h>
#include <comrmpi.h>
      integer lngrdfull,lngrdpart,ntotal                                  & 
     &       ,ii,n,mk,i,k,ierr,lend,strd
      real    a,b
      dimension a(lngrdfull,ntotal),b(lngrdpart,ntotal)
      real(kind=kind_mpi),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::len(:),loc(:)
!
      allocate(tmpsnd(lngrdfull*ntotal))
      allocate(tmprcv(lngrdpart*ntotal))
      allocate(len(0:npes-1))
      allocate(loc(0:npes-1))
!
!
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          loc(n)=mk
          lend=lngrdpart
          strd=lngrdpart*n
          do k=1,ntotal
            do i=1,lend
              ii=i+strd
              mk=mk+1
              tmpsnd(mk)=a(ii,k)
            enddo
          enddo
          len(n)=mk-loc(n)
        enddo
      else
        len(mype)=ntotal*lngrdpart
      endif
!
#if defined(CWB_MPMD) || defined(CWB_CFS)
      call mpi_scatterv(tmpsnd,len(0),loc(0),MPIREAL,                    & 
     &     tmprcv,len(mype),MPIREAL,0,mpi_comm_rsm,ierr)
#else      
      call mpi_scatterv(tmpsnd,len(0),loc(0),MPIREAL,                    & 
     &     tmprcv,len(mype),MPIREAL,0,mpi_comm_world,ierr)
#endif      
!
      mk=0
      lend=lngrdpart
      do k=1,ntotal
        do i=1,lend
          mk=mk+1
          b(i,k)=tmprcv(mk)
        enddo
      enddo
!
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(len)
      deallocate(loc)
!
#endif
      return
      end
