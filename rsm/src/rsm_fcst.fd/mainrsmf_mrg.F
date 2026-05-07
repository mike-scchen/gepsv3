#include <defopt.h>
      subroutine mainrsmf_mrg(compid, mpi_comm_mct)
!....
! add option to run from global as well as regional files
!....
      use machine , only : kind_io4
      use cnttime
!     use module_rinpinit,only: deallocate_rinpinitout_array
!yj---use module_rinpbase for 2nd cread----------------      
!     use module_rinpbase, only:  &
!    &  gzo,qo,teo,uuo,vvo,rqo,deallocate_rinpbaseout_array
!yj     &  gz=>gzo,q=>qo,te=>teo,du=>uuo,zv=>vvo,rq=>rqo 
!yj---use module_rinpbase for 2nd cread--- end   
!yj-note: vars. of rinpbase are allocated in 
!yj-note: its mainprogram.
!yj----------------------------------------------------
!yj---use module_gfsdata for 2nd cread----------------      
      use module_gfsdata, only: fhourg=>fhour
!     use module_gfsdata, only:   &
!    &  idateo=>idate,fhouro=>fhour, &
!    &  allocate_gfsdata,deallocate_gfsdata
!yj---- allocate_gfsdata in the recv_data subroutine---
#ifdef do3DRT
      use module_SurfaceAlbedo, only : deallocate_3drtdata
#endif
!yj----------------------------------------------------     
#ifdef wind_fitch
      use module_wind_fitch , only : finalized_wind_fitch_rsm , fi_nt 
#endif
!ccj
#include <paramodel.h>
!
      common /rcomfgi/ n0,n1,n2
#include <rscomio.h>
!
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
#ifdef G2R
#include <commpi.h>
#include <commpildef.h>
#endif
#endif
#include <rscomgpd.h>
!yj---mpmd itag common---
      common /mpmditag/itag
!yj---mpmd itag common---
!
      integer, intent(in) :: compid, mpi_comm_mct
      real(kind=kind_io4) fhours
      integer,dimension(4):: idate
!
      character cfhour*16,cform*40,filename*100
      integer iclim, ierr
      logical ldosfc,lsfcchg
      real rdummy(lngrdp)
      real*8 time_stamp, time_spend
!..........................................
!ty---wall time common---
#ifdef TIMING
      real*8 tm_1,tm_2,tm_use,mpi_wtime
#endif
!ty---wall time common---
#ifdef MP
      mpi_comm_rsm = mpi_comm_mct
      call rmpinit(ncoli,nrowi,strwtime)
      root_gfs = 0
      if(mype.eq.master) write(*,*) "[RSM] check root_gfs: ", root_gfs
      call rsm_cpl_init(compid)
#else
#ifdef IBMSP
      strwtime=timef()
#endif
#endif
!cc
!..........................................
!  
#ifdef MP
#ifdef G2R
       call mpdimset
#endif
       call rmpdimset
#endif
!.........................................
!c
#ifdef NONHYD
      if(mype.eq.master) print 100
100     format (1h0,'ncep nonhydrostatic mesoscale spectral model',      &
     &        '-current version created may 2004')
#else
      if(mype.eq.master) print 100
100   format (1h0,'ncep regional spectral model',                        &  
     &        '-current version created february 2003')
#endif
!jun
#ifdef CNT
       call cnttime_init()
#endif
!jun
      call rgetcon_mrg
#ifdef MP
!tst      if (mype.eq.master) then
#endif
!      print*,'mype=',mype,'at mainrsmf after rgetcon_mrg,rsfcsec=',      &
!     &     rsfcsec
#ifdef MP
!tst      endif
#endif

!
! rsm: run regional model only.
#ifdef MP
         if (mype.eq.master) then
#endif
#ifdef G2R
         print *,'rsm: run regional only with global files'
#endif
#ifdef C2R
         print *,'rsm: run regional with regional files'
#endif
#ifdef MP
         endif
#endif
#ifdef MP
      if (mype.eq.master) then
#endif
!yj ---------------------------------------------        
!yj --- mark the below lins                   ---
!yj --- replace the binary file input by mpmd ---
!yj ---------------------------------------------        
!yj      rewind nrsmi1
!yj      read(nrsmi1)
!yj      read(nrsmi1) fhours
!yj      fhour=fhours
!yj      rewind nrsmi1
!yj      kh=nint(fhour)
!     print*,'at mainrsmf',' fhour= ',fhour,' kh= ',kh
!yj      ndig=max(log10(kh+0.5)+1.,2.)
!yj      write(cform,'("(i",i1,".",i1,")")') ndig,ndig
!yj      write(cfhour,cform) kh
!hmhj open(unit=n1,file='rb_sigf'//cfhour                              & 
!hmhj&    ,status='old',form='unformatted')
!yj      write(filename,"(a,a)") 'rb_sigf',cfhour
!yj      open(unit=n1,file=filename                                       & 
!yj     &    ,status='old',form='unformatted')
!yj ---------------------------------------------        
#ifdef MP
      endif
#endif
#ifdef MP
      if (mype.eq.master) then
#endif
      print*,'at mainrsmf before rsmini,rsfcsec=',rsfcsec
#ifdef MP
      endif
#endif
      call rsmini_mpi(n0,n1)
#ifdef CPL_TIMCOM
      !call tai_timcom_initialize(mpi_comm_rsm)
      if(mype.eq.0) print*, 'check tai_timcom init:', cpl_kdt, cpl_time
#endif
#ifdef MP
      if (mype.eq.master) then
#endif
      print*,'at mainrsmf after rsmini'
!yj      close(n1)
#ifdef MP
      endif
#endif

! time loop
!jw for test
      maxstp=nint(fhseg/fhbas)
!yjprint maxstp-- start --
#ifdef MP
        if (mype.eq.master) then
#endif
          print*,'maxstp=',maxstp
#ifdef MP
        endif
#endif

      do n=1,maxstp

!yj ------call recv_data, recv. gfs data------- start        
!       call allocate_gfsdata
#ifdef MP
      if (mype.eq.master) then
        print*,'--yjcheck-- in mainrsm, tloop(n)=',n
#endif
        fhnxta=phour+fhbas

        print*,'fhnxta=',fhnxta
        print*,'phour=',phour
        print*,'fhbas=',fhbas
#ifdef MP
      endif
#endif

      fhourg=fhnxta
      call rmpbcastr(fhnxta,1)
      call rmpbcastl(lsfcmrg,1)

!notyet#ifdef CWB_CFS
!notyet!     call RSM_RECVA2Rdata
!notyet      stop '----error---'
!notyet#ifdef TIMING
!notyet      tm_1=mpi_wtime()
!notyet#endif
!notyet#endif

      kfhours=nint(fhnxta*3600)
      kdays=mod(kfhours,nint(rsfcsec))      
      if ( lsfcmrg .and. kdays.eq.0 .and. fhnxta .gt.0) then
        if (mype.eq.master) then
          print*,'------rsfcsec =',rsfcsec,'------'
          print*,'update the surface variables at fhour=',fhnxta
        endif
        iclim=0
        ldosfc=.true.
        lsfcchg=.true.
        time_stamp = MPI_WTIME()
#if defined (CWB_MPMD) || defined (CWBSUM)       
        call recv_data(root_gfs,iclim,ldosfc,lsfcchg,           &
     &                 mype,master,nrorg,rdummy,rdummy,         &
     &                 lonlen(mype),latlen(mype))
#else
#ifdef CWB_CFS
        call RSM_RECVA2Rdata(root_gfs,iclim,ldosfc,lsfcchg,     &
     &                 mype,master,nrorg,rdummy,rdummy,         &
     &                 lonlen(mype),latlen(mype))
#else
        call read_data(n2,iclim,ldosfc,lsfcchg,                 &
     &                 mype,master,nrorg,rdummy,rdummy,         &
     &                 lonlen(mype),latlen(mype))
#endif
#endif
        time_spend = MPI_WTIME() - time_stamp
        call MPI_ALLREDUCE(time_spend, time_stamp, 1, MPI_REAL8, MPI_MAX, mpi_comm_rsm, ierr)
        if(mype.eq.master) then
          write(*,'(a,f8.1,a)') '[RSM] recv gfs spends: ', time_stamp, "sec, surf on"
        end if
      else
        if (mype.eq.master) then
          print*,'not update the surface variables at fhour=',fhnxta
        endif
        iclim=0
        ldosfc=.false.
        lsfcchg=.false.
        time_stamp = MPI_WTIME()
#if defined (CWB_MPMD) || defined (CWBSUM)        
        call recv_data(root_gfs,iclim,ldosfc,lsfcchg,           &
     &                 mype,master,nrorg,rdummy,rdummy,         &
     &                 lonlen(mype),latlen(mype))
#else
#ifdef CWB_CFS
        call RSM_RECVA2Rdata(root_gfs,iclim,ldosfc,lsfcchg,     &
     &                 mype,master,nrorg,rdummy,rdummy,         &
     &                 lonlen(mype),latlen(mype))
#else
        call read_data(n2,iclim,ldosfc,lsfcchg,                 &
     &                 mype,master,nrorg,rdummy,rdummy,         &
     &                 lonlen(mype),latlen(mype))
#endif
#endif
        time_spend = MPI_WTIME() - time_stamp
        call MPI_ALLREDUCE(time_spend, time_stamp, 1, MPI_REAL8, MPI_MAX, mpi_comm_rsm, ierr)
        if(mype.eq.master) then
          write(*,'(a,f8.1,a)') '[RSM] recv gfs spends: ', time_stamp, "sec, surf off"
        end if
      endif  

#ifdef MP
      if (mype.eq.master) then
#endif
!yj-check the fhour in recv_data/read_data is the same as fhnxta 
      if (abs(fhourg-fhnxta) .gt. 0.0001) then
        print*,'fhour=',fhour,'fhnxta=',fhnxta
        print*,'Error! fhour should equal to fhnxta'
        stop
      endif  

#ifdef MP
      endif
#endif
!
!
#ifdef G2R
#ifdef MP
        ijmaxf=lnt22
        ijmaxp=lnt22p
#else
        ijmaxf=lnt22
        ijmaxp=lnt22
#endif
        kmax=levs
#endif
#ifdef C2R
        ijmaxf=clngrd
        ijmaxp=clngrd
        kmax=levr
#endif

        fhnxt=phour+fhbas
!yjprint fhnxt -- start --
#ifdef MP
        if (mype.eq.master) then
#endif
          print*,'fhnxt=',fhnxt
          print*,'phour=',phour
          print*,'fhbas=',fhbas
#ifdef MP
        endif
#endif
!yjprint fhnxt -- end --

#ifdef MP
        if (mype.eq.master) then
#endif
        kh=nint(fhnxt)
        ndig=max(log10(kh+0.5)+1.,2.)
        write(cform,'("(i",i1,".",i1,")")') ndig,ndig
        write(cfhour,cform) kh
!hmhj   open(unit=n2,file='rb_sigf'//cfhour                              &
!hmhj&      ,status='old',form='unformatted')
        write(filename,"(a,a)") 'rb_sigf',cfhour
        print*,'open n2, filename=',filename
!yj        open(unit=n2,file=filename                                       & 
!yj     &      ,status='old',form='unformatted')
#ifdef MP
        endif
#endif
!
        call rsmsmf(fhnxt,compid)
!
#ifdef MP
#ifdef TIMING
      if (mype.eq.master) then
        tm_2=mpi_wtime()
        tm_use=tm_2-tm_1
        print 265, phour,fhnxta,tm_use
  265   format( 'RSM intgrating ',f10.3,' to ',f10.3,' Timing=',f10.3,' elapse seconds')
      endif
#endif
#endif
      enddo
!yj======finish time integrate loop=========

#ifdef RKN
#ifdef MP
      call rmpkencollect(itnum,nvrken,nptken,nstken,svdata)
#endif
#endif
      call rsmsav

#ifdef do3DRT
      if (lsphy) then
      call deallocate_3drtdata
      endif
#endif

#ifdef wind_fitch
      if (lsphy) then
      call finalized_wind_fitch_rsm
      endif
#endif

#ifdef RKN
#ifdef MP
      if (mype.eq.master) then
#endif
!j      call rsmkenprof    !add by binbin to convert to eta bufr
#ifdef MP
      endif
#endif
#endif
!
!jun
#ifdef CNT
       call cnttime_final()
#endif
!jun
#ifdef MP
#ifdef NDSL
      call rndslfv_final
#endif
!      call rmpfine(endwtime)
      if (mype.eq.master) then
        print *,' MPI wall time is ',endwtime-strwtime
      endif
#else
#ifdef IBMSP
      endwtime=timef()
      print *,' SGL wall time is ',(endwtime-strwtime)/1000.
#endif
#endif
!c
!c    stop
      end
