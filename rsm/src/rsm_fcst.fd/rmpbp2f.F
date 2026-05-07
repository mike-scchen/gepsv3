      subroutine rmpbp2f(a,lngrdpart,b,lngrdfull,ntotal)
!$$$  subprogram documentation block
!
! subprogram: 	rmpbp2f
!            
! prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
!
! abstract:  transpose (ip,jp,kf) to (if,jf,kf)
!
! program history log:
!    99-06-27  henry juang 	finish entire test for gsm
!
! usage:	call rmpbp2f(a,igrd1p,jgrd1p,b,igrd1,jgrd1,ntotal)
!
!    input argument lists:
!	b	- real (igrd1p,jgrd1p,ntotal) partial field for each pe 
!	igrd1p	- integer partial grid in longitude time 2
!	jgrd1p	- integer partial grid in latitude divide 2
!	ntotal	- integer total set of fields
!
!    output argument list:
!	a	- real (igrd1,jgrd1,ntotal) full field 
!	igrd1	- integer total grid in longitude time 2
!	jgrd1	- integer total grid in latitude divide 2
! 
! subprograms called:
!   mpi_gatherv	- gather message from all pe to master
!
! attributes:
!    library: mpi
!    language: fortran 90
!$$$
#include <define.h>
#ifdef MP
      implicit none
#include <npesi.h>
#include <comrmpi.h>
      integer lngrdpart,lngrdfull,ntotal                           & 
     &       ,ii,n,i,k,mk,ierr,lend,strd
      real    a,b
      dimension a(lngrdpart,ntotal),b(lngrdfull,ntotal)
!
      real(kind=kind_mpi),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::len(:),loc(:)
!
      allocate(tmpsnd(lngrdpart*ntotal))
      allocate(tmprcv(lngrdfull*ntotal))
      allocate(len(0:npes-1))
      allocate(loc(0:npes-1))
!
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          loc(n)=mk
          len(n)=ntotal*lngrdpart
          mk=loc(n)+len(n)
        enddo
      endif
!      
      lend=lngrdpart
      mk=0
      do k=1,ntotal
        do i=1,lend
          mk=mk+1
          tmpsnd(mk)=a(i,k)
        enddo
      enddo
      len(mype)=mk
!
#if defined(CWB_MPMD) || defined(CWB_CFS)
      call mpi_gatherv(tmpsnd,len(mype),MPIREAL,                         &  
     &     tmprcv,len(0),loc(0),MPIREAL,0,mpi_comm_rsm,ierr)
#else      
      call mpi_gatherv(tmpsnd,len(mype),MPIREAL,                         &  
     &     tmprcv,len(0),loc(0),MPIREAL,0,mpi_comm_world,ierr)
#endif      
!
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          lend=lngrdpart
          strd=lngrdpart*n
          do k=1,ntotal
            do i=1,lend
              ii=i+strd
              mk=mk+1
              b(ii,k)=tmprcv(mk)
            enddo
          enddo
        enddo
!jw initialize for all other processors
      else
        b=0.
      endif
!      print*,'finish routine rmpbp2f'
!
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(len)
      deallocate(loc)
!
#endif
      return
      end
