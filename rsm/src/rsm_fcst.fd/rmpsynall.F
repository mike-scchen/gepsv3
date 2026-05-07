      subroutine rmpsynall
!$$$  subprogram documentation block
!
! subprogram: 	rmpsynall
!            
! prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
!
! abstract:  set mpi barrier for all nodes
!
! program history log:
!    99-06-27  henry juang 	finish entire test for gsm
!
! usage:	call rmpsynall
!
! subprograms called:
!   mpi_barrier  - set barrier for comm
!
! attributes:
!    library: mpi
!    language: fortran 90
!$$$
#include <define.h>
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
#if defined(CWB_MPMD) || defined(CWB_CFS)
      call mpi_barrier(mpi_comm_rsm,ierr)
#else      
      call mpi_barrier(mpi_comm_world,ierr)
#endif      
      if( ierr.ne.0 ) then
        print *,'PE',mype,':***** Error stop in mpsynall ******** '
        print *,'PE',mype,':error code from mpi_barrier = ',ierr
        call rmpabort
      endif
#endif
      return
      end
