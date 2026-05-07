!***********************************************************************
!
      MODULE RSM_cc

      USE CMP_COMM, ONLY: MPI_COMM_cRsm => COMM_local, Coupler_id, &
     &   component_master_rank_local, process_rank_local, &
!          Note: the latter two are only to compare with each
!          other and thus determine if the process is the local
!          master (root) process. (Comparison of
!          component_master_rank_global with process_rank_global
!          would not work because the former is known only to
!          Coupler process and the local master process itself.)
     &   component_nprocs, kind_REAL,MPI_kind_REAL,MPI_INTEGER,MPI_STATUS_SIZE, &
     &   ibuffer
      USE mpi_def, ONLY: COMM_TILES => MC_COMP

      implicit none

      integer N2D


      logical COMP /.false./

!controls:
!     integer nunit_announce_cc /6/, VerbLev /5/
!     integer nunit_announce_cc /6/, VerbLev /2/
      integer nunit_announce_cc /6/, VerbLev /3/
      save

      END MODULE RSM_cc
!
!***********************************************************************
!
      SUBROUTINE RSM_CMP_START

      USE RSM_cc, ONLY: component_nprocs,VerbLev,ibuffer,Coupler_id

      implicit none

      integer Rsm_id /3/, Rsm_master_rank_local /0/
      character*20 s
!

!        print*,'AM: to call CMP_INIT'
                      !<-id of AM as a component of the coupled system
      call CMP_INIT(Rsm_id,1)
                             !<-"flexibility level"
!        print*,'AM: back from CMP_INIT'
!     if (Coupler_id.ge.0) VerbLev=min(VerbLev,ibuffer(4))
      if (Coupler_id.ge.0) VerbLev=min(VerbLev,2)

        Rsm_master_rank_local=component_nprocs-1
                               !<- this redefinition is to meet the
                               ! requirement of subr. split2d_r used
                               ! in DISASSEMBLE_cc for disassembling
                               ! 2D fields. The requirement seems
                               ! to be that the input argument
                               ! representing a whole grid array be
                               ! defined in process of the largest rank
                               ! which seems to be considered i/o
                               ! process. To use a different value,
                               ! e.g. the conventional 0, split2d_r
                               ! (or DISASSEMBLE_cc) must be rewritten.
                     ! (Strangely, unsplit2d_r does not pose this
                     ! requirement and uses a dummy arg. ioproc to
                     ! identify the process where the whole grid array
                     ! is to be defined. Seemingly.)

      Rsm_master_rank_local=0    ! see above for modifications needed
                                 ! to support this change

      call CMP_INTRO(Rsm_master_rank_local)

      write(s,'(i2)') VerbLev
      call RSM_ANNOUNCE('back from CMP_INTRO, VerbLev='//s,2)

      return
      END
!
!***********************************************************************
!
      SUBROUTINE RSM_CMP_START1

      USE RSM_cc, ONLY: process_rank_local,VerbLev,ibuffer,Coupler_id

      implicit none

      integer Rsm_id /3/
!

                      !<-id of AM as a component of the coupling system
      call CMP_INIT(Rsm_id,1)
                             !<-"flexibility level"

!     if (Coupler_id.ge.0) VerbLev=min(VerbLev,ibuffer(4))
      if (Coupler_id.ge.0) VerbLev=min(VerbLev,2)

!           print*,'AM: back from CMP_INIT, process_rank_local=', &
!           process_rank_local

      return
      END
!
!***********************************************************************
!
      SUBROUTINE RSM_CMP_START2(me)      

      USE RSM_cc, ONLY: VerbLev

      implicit none

      integer me

      character*20 s
!

      if (me .eq. 0) then
        CALL CMP_INTRO_m
      else
        CALL CMP_INTRO_s
      end if

      write(s,'(i2)') VerbLev
      call RSM_ANNOUNCE('back from CMP_INTRO_m, VerbLev='//s,1)

      return
      END
!
!***********************************************************************
      SUBROUTINE RSM_RECVA2Rdate
      use RSM_cc, only: Coupler_id
      use module_gfsdata, only: idate
      implicit none

      integer,parameter :: ilen=1
      integer, dimension(ilen) :: F
      if (Coupler_id.lt.0) return    !   <- standalone mode
      call CMP_INTEGER_RECV(F,ilen)
      call RSM_ANNOUNCE('recv basefield idate=',2)
      idate=F(ilen)

      END subroutine RSM_RECVA2Rdate
!***********************************************************************
!     SUBROUTINE RSM_RECVA2Rdata
!     USE RSM_cc, ONLY: Coupler_id
!     use module_gfsdata,ONLY: clngrd,kmsoil,levs, &
!    & fhour,temp,spfh,clwr,ozon,geop,terr,   &
!    & u,v,tg,smc,snr,stc,cice,slmsk
!
!     implicit none

!     integer :: ilen1,ilen2,ilens,ilen3
!     real :: fhourin(1)
!     ilen1=1
!     ilen2=clngrd
!     ilens=clngrd*kmsoil
!     ilen3=clngrd*levs
!     if (Coupler_id.lt.0) return    !   <- standalone mode
!     call RSM_ANNOUNCE('start recv basefield data=',2)
!     call CMP_RECV(fhourin,ilen1)
!     fhour=fhourin(1)
!     call RSM_ANNOUNCE('recv basefield fhour =',2)
!     call CMP_RECV(temp,ilen3)
!     call CMP_RECV(spfh,ilen3)
!     call CMP_RECV(clwr,ilen3)
!     call CMP_RECV(ozon,ilen3)
!     call CMP_RECV(geop,ilen2)
!     call CMP_RECV(terr,ilen2)
!     call CMP_RECV(u,ilen3)
!     call CMP_RECV(v,ilen3)
!     call CMP_RECV(tg,ilen2)
!     call CMP_RECV(smc,ilens)
!     call CMP_RECV(snr,ilen2)
!     where (snr<1.e-3) snr=0.0
!     call CMP_RECV(stc,ilens)
!     call CMP_RECV(cice,ilen2)
!     call CMP_RECV(slmsk,ilen2)
!     where(slmsk.ge.0.01)
!       slmsk=1.0
!     elsewhere
!       slmsk=0.0
!     end where

!     END subroutine RSM_RECVA2Rdata
!
!***********************************************************************
!
      SUBROUTINE RSM_ANNOUNCE(s,DbgLev)

      USE RSM_cc, ONLY: nunit_announce_cc,VerbLev

      implicit none

      character*(*) s
      integer DbgLev
!
      if (DbgLev.le.VerbLev) &
     &  CALL CMP_ANNOUNCE(nunit_announce_cc,'RM: '//s)
!jwhwu
        call flush(nunit_announce_cc)
!jwhwu

      return
      END
!
!***********************************************************************
!
!
!***********************************************************************
!
       subroutine rsm_maxmin(xdim,ydim,x,s)
  
       USE RSM_cc

       implicit none

       integer xdim,ydim,i,j
       real(kind=kind_REAL) x(xdim,ydim),xmax,xmin
       character(*) s

      xmax=x(1,1)
      xmin=x(1,1)
      do j=1,ydim
      do i=1,xdim
       if ( xmax .lt. x(i,j) ) xmax=x(i,j)
       if ( xmin .gt. x(i,j) ) xmin=x(i,j)
      enddo
      enddo
      print *,s//' in rsm_maxmin,xdim=',xdim,'ydim=',ydim, &
     &   'xmax=',xmax,'xmin=',xmin

      return
      end
! ******************************************************************
       subroutine rsm_maxmin_int(xdim,ydim,x,s)

       USE RSM_cc

       implicit none

       integer xdim,ydim,i,j
       integer x(xdim,ydim),xmax,xmin
       character(*) s

      xmax=x(1,1)
      xmin=x(1,1)
      do j=1,ydim
      do i=1,xdim
       if ( xmax .lt. x(i,j) ) xmax=x(i,j)
       if ( xmin .gt. x(i,j) ) xmin=x(i,j)
      enddo
      enddo
      print *,s//' in rsm_maxmin,xdim=',xdim,'ydim=',ydim, &
     &   'xmax=',xmax,'xmin=',xmin

      return
      end

