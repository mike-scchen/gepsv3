      subroutine rmpinit(ncoli,nrowi,strwtime)
!$$$  subprogram documentation block
!
! subprogram: 	rmpinit
!            
! prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
!
! abstract: initializing the mpi by each pe.
!
! program history log:
!    99-06-27  henry juang 	finish entire test for gsm
!
! usage:	call rmpinit
!
! subprograms called:
!   mpi_init		- to initial mpi call
!   mpi_comm_size	- to get total number of pe for the comm
!   mpi_comm_rank	- to get my own pe number 
!   mpi_comm_split	- to split entire comm into sub comm
!
! attributes:
!    library: mpi
!    language: fortran 90
!$$$
#include <define.h>
#ifdef MP
#ifdef CWB_CFS
      use CMP_COMM, ONLY: process_rank_local
      use RSM_cc,   ONLY: MPI_COMM_cRsm
#endif
#include <npesi.h>
#include <comrmpi.h>

#ifdef CWB_MPMD
! for mpmd mode >
      call mpmd_init(npes,  mype , mpi_comm_rsm, root_gfs, info)
      ierr=info
      jerr=info
! for mpmd mode <
      strwtime=mpi_wtime()
#else
#ifdef CWB_CFS
      call mpi_init(info)
      strwtime=mpi_wtime()
! TYW >
      print*,'RM: back from MPI_INIT'
      call RSM_CMP_START
      print*,'RM: back from RSM_CMP_START',MPI_COMM_cRsm
      MPI_COMM_RSM=MPI_COMM_cRsm
      mype=process_rank_local
      jerr=info
      print*,'RM: ranking number ',mype
! TYW < 
      call mpi_comm_size(MPI_COMM_RSM,npes,ierr)
#else
     !call mpi_init(info)
     strwtime=mpi_wtime()
     info=0
     call mpi_comm_size(mpi_comm_world,npes,ierr)
     call mpi_comm_rank(mpi_comm_world,mype,jerr)
#endif
#endif

      if( ierr.ne.0 .or. jerr.ne.0 .or.                                  &  
     &    info.ne.0 .or. npes.ne.npesi ) then
        print *,'PE',mype,':********* Error stop in rmpinit ********* '
        print *,'PE',mype,':error code from mpi_init = ',info
        print *,'PE',mype,':error code from mpi_comm_size = ',ierr
        print *,'PE',mype,':error code from mpi_comm_rank = ',jerr
        print *,'PE',mype,':npes preset in code = ',npesi
        print *,'PE',mype,':npes request from environment = ',npes
        print *,'PE',mype,':******* End of output for rmpinit ******* '
        call rmpabort
      else
        master=0
        msgtag=0
        ncol=ncoli
        nrow=nrowi
        myrow=mype/ncol
        mycol=mod(mype,ncol)
#if defined(CWB_MPMD) || defined(CWB_CFS)
        call mpi_comm_split(mpi_comm_rsm,  myrow,mycol,comm_row,ierr)
        call mpi_comm_split(mpi_comm_rsm,  mycol,myrow,comm_column,jerr)
#else
        call mpi_comm_split(mpi_comm_world,myrow,mycol,comm_row,ierr)
        call mpi_comm_split(mpi_comm_world,mycol,myrow,comm_column,jerr)
#endif       
        if( ierr.ne.0 .or. jerr.ne.0 ) then
          print *,'PE',mype,':error code for doing comm_row = ',ierr
          print *,'PE',mype,':error code for doing comm_column = ',jerr
          call rmpabort
        else
          print *,'PE',mype,' npes ncol nrow mycol myrow',                 &
     &                        npes,ncol,nrow,mycol,myrow
        endif
      endif
#endif
      return
      end

!     function mpi_wtime( )
!     real*8 mpi_wtime
!     mpi_wtime = timef()
!     return
!     end
