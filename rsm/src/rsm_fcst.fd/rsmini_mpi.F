#include <defopt.h>
      subroutine rsmini_mpi(n0,n1)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rsmini
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  prepare all the constant arraies for regional computation
!      before regional forecast and read the input data set.
!
! program history log:
!
! usage:    call rsmini(gfhour,filtai,
!          1                   nin1i,nin2i,
!          2                   nrsmi1i,nrsmi2i,nrflipi,
!          3                   nrsmo1i,nrsmo2i,nrflopi,nrsflxi,nriniti)
!   input argument list:
!     gfhour  - forecast hour of global model
!     nin1i - unit number for global sigma file input 1
!     nin2i - unit number for global sigma file input 2
!     nrsmi1i - unit number for regional sigma file input 1
!     nrsmi2i - unit number for regional sigma file input 2
!     nrflipi - unit number for regional surface file input
!     nrsmo1i - unit number for regional sigma file outout 1
!     nrsmo2i - unit number for regional sigma file outout 2
!     nrflopi - unit number for regional surface file outout
!     nrsflxi - unit number for regional flux file outout
!     nriniti - unit number for regional initialization temperory file
!
!   output argument list:
!
!   common block:
! rscomf  rscomltb  comcon
!   input files: none
!
!   output files: none
!
!   subprograms called:
! setbgd  setg2b  getwt  rftini plnini rltbini  rsminp  rsetsig
! rsetmap  ramtm  epslxy  pint3  rsetsfc  rfixio
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
!.................................................................
!..... this routine setup all the routine and constant for rsm.
!..... and get input data from gsm or others include surface files.
!.....
      use machine , only : kind_io4
      use physcons, pi => con_pi
      use sigmap_cons, only : allocate_sigma_array,sigma
!----------------------------------------------------------------------      
!------- All real variables from gfs are double precision.   ----------
!----------------------------------------------------------------------      
      use module_gfsdata, only:                                        &
     &  idateg=>idate,fhourg=>fhour,                                   &
     &  config_gfsdata,clat,clon,                                      &
     &  allocate_climdata,deallocate_climdata,allocate_gridata,        &
     &  tg3clm,zorclm,albclm,vegclm,vetclm,sotclm,alfclm,              &
     &  vmnclm,vmxclm,slpclm,absclm,fflat,fflon,rslm
#ifdef do3DRT
      use module_SurfaceAlbedo,only: allocate_3drtdata,read_3drtdata
#endif

!ccj
#include <paramodel.h>
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
      real, allocatable :: grid1(:,:,:)
      real, allocatable :: coef(:,:)
#define IGRD1S igrd1p
#define JGRD1S jgrd1p
#define LNWAVS lnwavp
#define LNGRDS lngrdp
#else
#define IGRD1S igrd1
#define JGRD1S jgrd1
#define LNWAVS lnwav
#define LNGRDS lngrd
#endif
#include <rscomfi.h>
#include <rscomspc.h>
#include <rscomgrd.h>
#include <rscomver.h>
#include <rscomio.h>
#include <rscompln.h>
#include <rscomrad.h>
#include <rscomloc.h>
#include <rscomltb.h>
#include <rscomgpd.h>
#include <rscomoz.h>
!---mpmd itag common---
      common /mpmditag/itag
!cc
      dimension gzzb(IGRD1S,JGRD1S),                                     &  
     &          gzxb(IGRD1S,JGRD1S),gzyb(IGRD1S,JGRD1S)                  
!cc Jun's chnage for g2rini
      dimension  gclat(latg)                                             &
     &  , colrad(latg2), wgt(latg2),wgtcs(latg2),rcs2(latg2)
!cc
!
      parameter (km=1)
!
      dimension tmpsp(lnwav,levrp1)
      real(kind=kind_io4), allocatable :: tmps(:,:)
      real(kind=kind_io4) fhours,sis(levr+1),sls(levr)
!
      character cfhour*16,cform*40
!
! restart
      logical ireadbase
! 
      namelist/namloc/                                                  &
     & rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd &
     &,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd &
     &,clat1,clat2,clon1,clon2

! climate data namlist
      character*500          fntg3c,fnzorc,fnalbc,fnvegc,fnvetc,fnsotc, &
     &                 fnmskh,fnvmnc,fnvmxc,fnslpc,fnabsc
      namelist/namclim/iclim,fntg3c,fnzorc,fnalbc,fnvegc,fnvetc,fnsotc, &
     &                 fnmskh,fnvmnc,fnvmxc,fnslpc,fnabsc
      real, allocatable :: dflat(:),dflon(:)

      logical ldosfc,lsfcchg
! base field mapping array
      real,allocatable ::  cmap(:),bmap(:)
! model grid mapping array
!      real,allocatable :: fflat(:),fflon(:)         ! only used by mype=master
      real,allocatable :: cm2(:)                    ! only used by mype=master
      real,allocatable :: frcsln(:),frsnln(:)       ! only used by mype=master
      real,allocatable :: ffm2(:),ffm2x(:),ffm2y(:) ! only used by mype=master
      dimension  rcsln(lngrdp),rsnln(lngrdp)        ! for rotate u v to new projection
!      real,allocatable :: rslm(:)                   ! only used by mype=master
! tmp
      real*4,dimension(:),allocatable:: tmp4

#ifdef KEout
! yj--add output the kinetic energy of each zonal wave number >
      real,dimension(:,:),allocatable::  &
     &      utmp, vtmp, ttmp,ptmp, ke 
      real,dimension(:),allocatable::  &
     &      kezon, kemer
! yj--add output the kinetic energy of each zonal wave number <
#endif


#ifdef MP
      if (mype.eq.master) then
#endif

      print *,' ----- into rsmini ----- '
!
      print *,' ++++++++ header of rsm +++++++++++++++'
      print *,' ++ version with nest to global model +'
      print *,' ++ iwave=', iwav,' jwave=', jwav 
      print *,' ++ levs =', levr,' levh =', levh  
      print *,' ++ igrid=', igrd1,' jgrid=', jgrd1  
      print *,' ++ ncol=', ncol,' nrow=', nrow  
#ifdef MP
      endif
#endif

      if (lrestart) then
        fhour=restrhr
        if (mype.eq.master) print*,'run re-start'
      else
        fhour=0.0
        if (mype.eq.master) print*,'NOT re-start'
      endif

! get sigma-p coeff. of RSM vertical coordinate
      call allocate_sigma_array(levr)
!
! receive idate from GFS
!
      if (mype.eq.master) then
#if defined(CWBSUM) || defined(CWB_MPMD)
        call recv_idate(root_gfs,itag)
#else      
#ifdef CWB_CFS
        call RSM_RECVA2Rdate
#else
        kh=nint(restrhr)
        ndig=max(log10(kh+0.5)+1.,3.)
        write(cform,'("(i",i1,".",i1,")")') ndig,ndig
        write(cfhour,cform) kh
        open(unit=n0,file='rsm_idate',convert='little_endian' &
     &        ,status='unknown',form='unformatted')
        call read_idate(n0)
#endif
#endif
        call newidate(idateg,0.0,idate)
      endif
!
! get GFS data configuration
!   read namlist
      if ( mype.eq.master ) then
        rewind(5)
        read (5,namloc)
        write (6,namloc)
      endif
      call config_gfsdata(mype,master,root_gfs,cigrd1,cjgrd1,itag)
!
! ..... setup base field grid
!
#ifdef MP
      istr=lonstr(mype)
      iend=lonstr(mype)+lonlen(mype)-1
      jstrs=latstr(mype)
      jends=latstr(mype)+latlen(mype)-1
      jstrn=jgrd1-jends+1
      jendn=jgrd1-jstrs+1
      ilen=lonlen(mype)  
#else
      istr=1
      iend=igrd1
      jstrs=1
      jends=jgrd1
      jstrn=jgrd1+1
      jendn=jgrd1
      ilen=igrd1     
#endif
      ib1=(istr-1)/bgf-border+1
      ib2=ibgd1-1+ib1
      jb1=(jstrs-1)/bgf-border+1
      jb2=jbgd1-1+jb1
      jbx=jb2
#ifdef MP
      if (mype.eq.master) print *, 'ib1=',ib1,'jb2=',ib2,                &  
     &   'jb1=',jb1,'jb2=',jb2,'jbx=',jbx,                               &
     &   'bgf=',bgf,'border=',border,'ibgd1=',ibgd1,'jbgd=',jbgd1        &
     &    ,'igrd1p=',igrd1p,'jgrd1p=',jgrd1p
#endif
      call setbgd(blat(1:lngrdb),blon(1:lngrdb),delx,dely,dlamda0,       &
     &            bgf,ib1,ib2,jb1,jb2                                    &
#ifdef MP
     &           ,mype,master                                            &
#endif
     &            )

!
! ..... prepare for transform gsm to rsm
! (1) base field grid setting  (after this, can call cgtobg)
!
      call ctobini(clat,clon,blat,blon                                     &
#ifdef MP
     &           ,mype,master                                            &
#endif
     &            )
      call getwt(border*2-bsmooth,bsmooth,bgf,wsltb,dwsltb)

      ! setup u v rotate for rsm2bgd
      allocate(cmap(lngrdb),bmap(lngrdb))
      call cgtobg(clat,cmap,1)
      ! define cgrid project transformation for u and v
       rad = pi / 180.0
       hfpi= 0.5 * pi
       truth = rtruth
       nproj = rproj
       orient = rorient
       dlamda0 = orient * rad + hfpi
       ncproj = cproj
       cdlamda0 = corient * rad + hfpi
       do n=1,lngrdb
         ! base field mapping
         if( nproj.eq.1 .or. nproj.eq.-1 ) then
           a2 = ( 1.0 + sin( nproj * truth * rad ) )
           bmap(n)=a2/(1.+sin( nproj * blat(n) ))
           rgln = blon(n) - dlamda0
           rgln = nproj * rgln
           csmm = nproj * cos( rgln ) 
           snmm = nproj * sin( rgln )
         else if( nproj.eq.0 ) then
           a2 = cos( truth * rad )
           bmap(n)=a2/cos( blat(n) )
           csmm =  0.0
           snmm = -1.0 
         endif
         ! GFS data mapping
         if( ncproj.eq.0 ) then
           csll =  0.0
           snll = -1.0 
         endif
!!!! some bug in this mapping rotation, skip for test! >>>
!        cmorm = sqrt( cmap(n) )/ bmap(n)
!        bcsln(n) = ( csll * csmm + snll * snmm ) * cmorm
!        bsnln(n) = ( csll * snmm - snll * csmm ) * cmorm
!!!! some bug in this mapping rotation, skip for test! <<<
       enddo
      deallocate(cmap,bmap)
      !
      !(2) model grid setting
      if (mype.eq.master) then
!        allocate(fflat(lngrd),fflon(lngrd))
        call allocate_gridata(lngrd)
        call setlatlon(fflat,fflon,delx,dely,dlamda0)
        call rtorini(clat,clon,fflat,fflon)               !after this, can call cgtorg
      endif

      call rmpgf2p(fflat,igrd1,jgrd1,flat,igrd1p,jgrd1p,1)
      call rmpgf2p(fflon,igrd1,jgrd1,flon,igrd1p,jgrd1p,1)

      !(3) model grid for surface field
      call c2rini_mpi(clat,clon,flat,flon,mype,master)
! ................................................
! define regional map factors
      if ( mype.eq.master ) then

        allocate ( cm2(lngrd))
        allocate ( frcsln(lngrd),frsnln(lngrd))
        allocate ( ffm2(lngrd),ffm2x(lngrd),ffm2y(lngrd) )
        rad = pi / 180.0
        hfpi= 0.5 * pi
        truth = rtruth
        nproj = rproj
        if( nproj.eq.1 .or. nproj.eq.-1 ) then
          orient = rorient
          dlamda0 = orient * rad + hfpi
          a2 = ( 1.0 + sin( nproj * truth * rad ) )
          do n=1,lngrd
            ffm2(n)=(1.+sin( nproj * fflat(n) ))/a2  ! 1/(map factor)
            rgln = fflon(n) - dlamda0
            rgln = nproj * rgln
            frcsln(n) = nproj * cos( rgln ) 
            frsnln(n) = nproj * sin( rgln )
          enddo
        else if( nproj.eq.0 ) then
          a2 = cos( truth * rad )
          do n=1,lngrd
            ffm2(n)=cos( fflat(n) )/a2               ! 1/(map factor)
            frcsln(n) =  0.0
            frsnln(n) = -1.0 
          enddo
        endif
! define cgrid project transformation for u and v
        do n=1,lngrd
!         if (nproj.eq.4) then
            cm2(n)=1.0
            csll = 0.0
            snll = -1.0
!         endif
!!!! some bug in this mapping rotation, skip for test! >>>
!         csmm = frcsln(n)
!         snmm = frsnln(n)
!         cmorm = sqrt( cm2(n) / ffm2(n) )
!         frcsln(n) = ( csll * csmm + snll * snmm ) * cmorm
!         frsnln(n) = ( csll * snmm - snll * csmm ) * cmorm
!!!! some bug in this mapping rotation, skip for test! <<<
        enddo
!
        call todxdy(ffm2,ffm2x,ffm2y,igrd1,jgrd1,rdelx,rdely)
!
        do ij=1,lngrd
        ffm2(ij) = 1.0/ ( ffm2(ij) * ffm2(ij) )
        ffm2x(ij) = -2.0 * ffm2(ij)**1.5 * ffm2x(ij)
        ffm2y(ij) = -2.0 * ffm2(ij)**1.5 * ffm2y(ij)
        enddo
        call maxmin(ffm2,lngrd,1,1,1,'fm2')
        call maxmin(ffm2x,lngrd,1,1,1,'fm2x')
        call maxmin(ffm2y,lngrd,1,1,1,'fm2y')

      endif ! (mype.eq.master)

      call rmpgf2p(frcsln,igrd1,jgrd1,rcsln,igrd1p,jgrd1p,1)
      call rmpgf2p(frsnln,igrd1,jgrd1,rsnln,igrd1p,jgrd1p,1)
      call rmpgf2p(ffm2,igrd1,jgrd1,fm2,igrd1p,jgrd1p,1)
      call rmpgf2p(ffm2x,igrd1,jgrd1,fm2x,igrd1p,jgrd1p,1)
      call rmpgf2p(ffm2y,igrd1,jgrd1,fm2y,igrd1p,jgrd1p,1)

      if ( mype.eq.master ) then
        deallocate ( ffm2,ffm2x,ffm2y )
        deallocate ( frcsln,frsnln )
        deallocate ( cm2 )
      endif
!    
! read(regional land-sea mask [0/1/2]
        if( mype.eq.master ) then
!          allocate(rslm(lngrd),tmp4(lngrd))
          allocate(tmp4(lngrd))
          rewind(nrslm)
          read(nrslm) tmp4
          rslm=tmp4
          deallocate(tmp4)
        endif
        call rmpgf2p(rslm,igrd1,jgrd1,slmsk,igrd1p,jgrd1p,1)
!
! ---------------- read climate data on regional grid--------------------
      if (mype.eq.master) then
        rewind(5)
        read(5,namclim)
      endif
      call rmpbcasti(iclim,1)
      if (iclim.eq.1) then
      if (mype.eq.master) then
        call allocate_climdata(lngrd)
        allocate ( dflat(lngrd),dflon(lngrd) )
        ! degree
        dflat=fflat*180./pi
        dflon=fflon*180./pi
        call rclima(idate(4),idate(2),idate(3),idate(1),fhour,            &
     &             lngrd,lsoil, rslm,dflat,dflon,                         &
     &             fntg3c,fnzorc,fnalbc,fnvegc,fnvetc,fnsotc,fnmskh,      &
     &             fnvmnc,fnvmxc,fnslpc,fnabsc,                           & 
     &             tg3clm,zorclm,albclm,vegclm,vetclm,sotclm,alfclm,      &
     &             vmnclm,vmxclm,slpclm,absclm)
        deallocate ( dflat,dflon )
      endif
        call clima_f2p
      if (mype.eq.master) then
        call deallocate_climdata
      endif

! when read climate data, give default initial value for the following variables
        do i=1,igrd1p
        do j=1,jgrd1p
            cv (i,j)=0.0     ! cloud amount
            cvt(i,j)=0.0     ! cloud top
            cvb(i,j)=100.0   ! cloud base
        enddo
        enddo
        do i=1,igrd1p
        do j=1,jgrd1p
#ifdef doNOAHLSM
          canopy(i,j)=60.
#else
          canopy(i,j)=0.6
#endif
          f10m(i,j)=0.8
          uustar(i,j) = 1.
          ffmm(i,j) = log(30.)
          ffhh(i,j) = log(30.)
        enddo
        enddo
      endif  ! iclim
#ifdef CPL_TIMCOM
      cpl_field = 0.0
#endif
!
! deallocate some arrays
!      if (mype.eq.master) then
!        deallocate(rslm)
!        deallocate(fflat,fflon)
!      endif
!===============================================================
!
!....< set up some constants for model >
!
#ifndef T0SI
#define T0SI 300.0
#endif
      do 30 lev=1,levr
      tov(lev)= T0SI
 30   continue
!
      call rftini
! jun fft
#ifndef MP
      call plnini
#endif
      call rltbini
!
!*************get initial data****************
! start receive data from GFS
!
      fhourg=fhour
      ldosfc=.true.
      lsfcchg=.false.
#if defined(CWBSUM) || defined(CWB_MPMD)
      call recv_data(root_gfs,iclim,ldosfc,lsfcchg,                &
     &               mype,master,nrorg,rcsln,rsnln,                &
     &               lonlen(mype),latlen(mype))
#else
#ifdef CWB_CFS
      call RSM_RECVA2Rdata(root_gfs,iclim,ldosfc,lsfcchg,          &
     &               mype,master,nrorg,rcsln,rsnln,                &
     &               lonlen(mype),latlen(mype))
#else
      call read_data(n1,iclim,ldosfc,lsfcchg,                      &
     &               mype,master,nrorg,rcsln,rsnln,                &
     &               lonlen(mype),latlen(mype))
#endif
#endif
!
      if(fhour.eq.0.0) then
        if (lrestart) stop "Error: don't use restart when fhour=0"
!
        if(ncld.gt.0) then
         if ( nmmiph.eq.4 .or. nmmiph.eq.6 ) then    ! WSM3/WSM6 (num_p3d=3)
          do k=1,levr
          do j=1,JGRD1S
          do i=1,IGRD1S
             phy_f3d(i,j,k,1) = 10.    !re_cloud
             phy_f3d(i,j,k,2) = 50.    !re_ice
             phy_f3d(i,j,k,3) = 250.   !re_snow
          enddo
          enddo
          enddo
         elseif ( nmmiph.eq.11 .or. nmmiph.eq.12 .or. nmmiph.eq.13 ) then     ! GFDL MP v1/v2/v3 (num_p3d=5)
          do k=1,levr
          do j=1,JGRD1S
          do i=1,IGRD1S
             phy_f3d(i,j,k,1) = 10.    !re_cloud
             phy_f3d(i,j,k,2) = 50.    !re_ice
             phy_f3d(i,j,k,3) = 1000.  !re_rain
             phy_f3d(i,j,k,4) = 250.   !re_snow
             phy_f3d(i,j,k,5) = 1000.  !re_graupel
          enddo
          enddo
          enddo
         else   !Zhao-Carr(num_p3d=4),Brad-Ferrier(num_p3d=3)
          do n=1,num_p3d
          do k=1,levr
          do j=1,JGRD1S
          do i=1,IGRD1S
             phy_f3d(i,j,k,n) = 0.0
          enddo
          enddo
          enddo
          enddo
         endif
          do n=1,num_p2d
          do j=1,JGRD1S
          do i=1,IGRD1S
             phy_f2d(i,j,n)=0.0
          enddo
          enddo
          enddo
         endif
        do k=1,levr
        do j=1,JGRD1S
        do i=1,IGRD1S
           cldcov(i,j,k)=0.0
        enddo
        enddo
        enddo
      else ! fhour.eq.0.0
#ifdef MP
      if (mype == master) then
#endif
        print *, '---------- begin rsm restart from file ----------'
        print *, 'file base       =', nb
        print *, 'file unit1      =', nrsmi1
        print *, 'file unit2      =', nrsmi2
        print *, 'iorestartfile   =', lrestart

        open(unit=nb, file='r_sigbase', action='read', &
     &       form='unformatted', status='old')
        open(unit=nrsmi1, file='r_sigf', action='read', &
     &       form='unformatted', status='old')
        open(unit=nrsmi2, file='r_sigftdt', action='read', &
     &       form='unformatted', status='old')
#ifdef MP
      end if
#endif
      if (lrestart) then
! get base field
      call readbase(nb,qb,teb,uub,vvb,rqb,lngrdb,levr,ntotal)
! get init field
      call sreadwav(nrsmi1,fhour,idate,                                 &
     &           flat,flon,fm2,fm2x,fm2y,sl,si,                         &
     &           gz,qm,tem,uum,vvm                                      &
#ifdef NDSL
     &           ,rmgrd                                                 &
#else
     &           ,rm                                                    &
#endif
#ifdef NONHYD
     &                     ,pnm,tnm,onm,tov                             &
#endif
     &                                       )
      call sreadwav(nrsmi2,fhour,idate,                                 &
     &           flat,flon,fm2,fm2x,fm2y,sl,si,                         &
     &           gz,q,te,uu,vv                                          &
#ifdef NDSL
     &           ,rqgrd                                                 &
#else
     &           ,rq                                                    &
#endif
#ifdef NONHYD
     &                     ,pn,tn,on,tov                                &
#endif
     &                                       )
      endif
        if(ncld.gt.0) then
          if ((nmmiph .eq. 2) .or.             & ! Zhao-Carr
         &    (nmmiph .eq. 3) .or.             & ! Brad-Ferrier
#ifdef NEW_WSM3
         &    (nmmiph .eq. 4) .or.             & ! WSM3
#endif
#ifdef NEW_WSM6
         &    (nmmiph .eq. 6) .or.             & ! WSM6
#endif
         &    (nmmiph .eq. 11) .or.            & ! GFDL MP v1
         &    (nmmiph .eq. 12) .or.            & ! GFDL MP v2
         &    (nmmiph .eq. 13) .or.            & ! GFDL MP v3
         &    (nmmiph .eq. 15)     ) then        ! GCE 3ICE
            ! Restart additional microphysics variables
            call sread2(nrsmi2,phy_f3d,phy_f2d,lrestart)
          endif
        endif
      if (mype.eq.master) then
        close(nb)
        close(nrsmi1)
        close(nrsmi2)
      endif
      endif
!yj ------ this endif (fhour.eq.0.0) -------
!
#ifdef MP
      if (mype.eq.master) then
#endif
      print*,'in rsmini_mpi.F, fhour and idate=',fhour,idate
      print 20,fhour,idate
 20   format(1h ,'fhour=',f6.2,2x,' date= ',4(i4,2x))
#ifdef MP
      endif
#endif
!
      if (abs(fhour-restrhr).gt.1.E-5) then
        print*,'STOP:---forecast already finished---'
        print*,'STOP:---forecast hour is',nsseg*deltim/3600.
        print*,'STOP:---restart hour is ',restrhr
        print*,'STOP:---current process hour is ',fhour
        stop
      endif
      kdt=nint(fhour*3600/deltim)
      fhour=kdt*deltim/3600
      phour=fhour
!
!cc
! rsetmap:  pass in latitudes longitudes map factor, and compute the
!            mean map factor or the largest map factor and it derivatives
!            in x and y directions.
!c       call rsetmap(flat(1),flon(1),fm2,fm2x,fm2y,delx,dely)
!c       call rsetmap(flat(n),flon(n),fm2,fm2x,fm2y,delx,dely)
      call rsetmap(flat,flon,fm2,fm2x,fm2y,delx,dely)

#ifndef SIGP
! get si,sl,del
      ps0=101.36
      do k=1,levrp1
        si(k)=(sigma(k,2)+sigma(k,1)*ps0*10.)/(ps0*10.)
      enddo
      do k=1,levr
        sl(k)=(si(k)+si(k+1))/2.
        del(k)=si(k)-si(k+1)
      enddo
!sigma
      call rsetsig(levr, levrm1, levrp1,                                 &  
     &             ci, si, del, sl, cl, rpi)
      sl1=sl(1)
      do 40 k=1,levrm1
      rpirec(k) = 1.e0/rpi(k)
 40   continue
!
      do 50 k=2,levr
      rdel2(k)=0.5e0 / del(k)
 50   continue
      rdel2(1)=0.5e0/del(1)
      rdel2(levrp1)=0.5e0/del(levr)
!
!     do 50 k=1,levr
!     rdel2(k)=0.5e0/del(k)
!50   continue
!      
! ramtm:  prepare some constant arraies for dynamic computation and semi-
!            implicit computation.
      call ramtm(del,rpi,sv,p1,p2,am,hm,tm,si)
      call rbmcm(tov,p1,p2,h1,h2,del,ci,bm,cm,sv,am)
#ifdef NONHYD
      call cdefghm(cmn,dmn,emn,fmn,gmn,hmn,                              &  
     &                   si,sl,tov,del,h1,h2)
#endif
!
#else
!sigma-P
!     tov - reference temperature [unit: K], given above, default=300.0
!     ps0 - reference surface pressure [unit: cb]
      ps0=101.36
      tov(:)=300.0
      call ramhm_sp(ps0,sigma,tov,sv,p1,p2,am,hm)
      call rbmcm_sp(ps0,sigma,tov,sv,p1,p2,am,hm,h1,h2,bm,cm)
#ifdef NONHYD
      stop 'sigma-P for MSM not implemented'
#endif
! use si,sl,del for reference sigma value (for physics parameterization)
      do k=1,levrp1
        si(k)=(sigma(k,2)+sigma(k,1)*ps0*10.)/(ps0*10.)
      enddo
      do k=1,levr
        sl(k)=(si(k)+si(k+1))/2.
        del(k)=si(k)-si(k+1)
      enddo
#endif
!
! epslxy:  prepare derivative and laplacian for regional wave field.
!            for regional spectral model.
!
#ifdef MP
      allocate (coef(lnwav,4))
      call epslxy(lnwav,iwav1,jwav1,igrd,jgrd,                           &  
     &            coef(1,1),coef(1,2),coef(1,3),coef(1,4),delx,dely)
      rnnp1max=coef(lnwav,4)
      epsxmax=coef(lnwav,1)
      epsymax=coef(lnwav,2)
      call rmpsf2p(coef(1,1),lnwav,epsx,lnwavp,1)
      call rmpsf2p(coef(1,2),lnwav,epsy,lnwavp,1)
      call rmpsf2p(coef(1,3),lnwav,snnp1,lnwavp,1)
      call rmpsf2p(coef(1,4),lnwav,rnnp1,lnwavp,1)
      deallocate (coef)
#else
      call epslxy(lnwav,iwav1,jwav1,igrd,jgrd,                           &  
     &            epsx,epsy,snnp1,rnnp1,delx,dely)
      rnnp1max=rnnp1(lnwav)
      epsxmax=epsx(lnwav)
      epsymax=epsy(lnwav)
#endif
!
      do nn=1,LNWAVS
        dpsx(nn) = epsx(nn)*delx
        dpsy(nn) = epsy(nn)*dely
      enddo
!
#ifdef MP
#define LATLENS latlen(mype)
#else
#define LATLENS jgrd1
#endif
      do j=1,LATLENS
         lat=j+jstrs-1
         call pint3(gzb,gzzb(1,j),gzxb(1,j),gzyb(1,j)                    &  
     &             ,wsltb,dwsltb,lat,km,IGRD1S,ib1,ib2                   &
     &             ,jb1,jb2,jbx,istr,ilen                                &
     &             ,bgf,border,jgrd1)
      enddo
      call rsetsfc(gz,gzzb,gzxb,gzyb,epsx,epsy)
!....< initial setup for model physics >
!ccj
!
!c set ozon indices for interpolation
!c
      call setindxoz(LNGRDS,lats18,flat,ddy,jindx1,jindx2)
!c
! read climatologically varying ozone prod and dist rate
!
      call readoz_disprd(ozprdin,ozdisin,poz)
!c
#ifdef do3DRT
      if (lsphy) then
! read SurfaceAlbedoMod variables
!
      call allocate_3drtdata(IGRD1S,JGRD1S)
      call read_3drtdata(igrd1,jgrd1,IGRD1S,JGRD1S)
      endif
#endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccc
      rcl=1.0
      iflws=0
!     read radiation files and some constants
!c      call gcljms(si)
!c      call grlwsw(sl,nco2)
!c      call gaeros(si,sl)
!....< set initial solar hour >
      solhr = fhour + idate(1)
      iday = solhr / 24.e0
      solhr = solhr - iday*24.e0
      if( solhr.lt.0.0 ) solhr=solhr+24.e0
!....< read fixed fields from regional fixfld prog>
#ifdef MP
      if (mype.eq.master) then
#endif
      print *,' -rsmini-read surface from regional.rsfcsec=',rsfcsec
#ifdef MP
      endif
#endif
      call rfixio(fhour,idate,                                           &  
     &   tsea,smc,sheleg,stc,tg3,zorl,                                   &
     &   cv,cvb,cvt,slmsk,f10m,canopy,vfrac,vtype,stype,                 &
     &   shdmin,shdmax,slope,snoalb,                                     &
     &   uustar,ffmm,ffhh,alvsf,alvwf,alnsf,alnwf,facsf,facwf,           &
     &   0,nrflip,rsfcsec,lsfcmrg,lrestart)
!snwdph initial value
        do j = 1, jgrd1p
          do i = 1, igrd1p
            if(slmsk(i,j).eq.1.) then
              snwdph(i,j) = sheleg(i,j) / 100.0
            else if (slmsk(i,j).eq.2) then
              snwdph(i,j) = sheleg(i,j) / 1000.0 * 3.0
            endif
          enddo
        enddo
!
#ifdef NONHYD
      if(fhour.eq.0.0) then
        mid=0.5*LNWAVS
        dzdtsum=0.0
        do k=1,levr+1
          dzdtsum=dzdtsum+abs(on(mid,k))
        enddo
        if( dzdtsum.eq.0.0 ) then
          call rloopw
!
          call rmsgtr( q, pn, tn, on,del,rq)
          do k=1,levr+1
            do i=1,LNWAVS
              onm(i,k)=on(i,k)
            enddo
          enddo
        endif
      endif
#endif

!testjun
     if (mype.eq.master) then
          print *,' in rsmini_mpi,after (nonhydrostatic) '
     endif
!....read mountain variance etc
!
      allocate (tmps(lngrd,nmtvr))
#ifdef MP
      allocate (grid1(igrd1,jgrd1,nmtvr))
#define HPRIMES grid1
#else
#define HPRIMES hprime
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
         rewind(nmtnv)
         read(nmtnv) tmps
         call sgl2ful(tmps,HPRIMES,lngrd*nmtvr)
#ifdef MP
      endif
      call rmpgf2p(grid1,igrd1,jgrd1,                                    &  
     &            hprime,igrd1p,jgrd1p,nmtvr)
#endif
!
      do k=1,nmtvr
        do j = 1, JGRD1S
          do i = 1, IGRD1S
            if(slmsk(i,j).ne.1.) hprime(i,j,k) = 0.
          enddo
        enddo
      enddo
#ifdef MP
      call rmpgp2f(hprime,igrd1p,jgrd1p,grid1,igrd1,jgrd1,nmtvr)
      if( mype.eq.master ) then
#endif
        call maxmin(grid1,lngrd,nmtvr,1,nmtvr,'hprimes full')
!       call maxmin(hprime,IGRD1S*JGRD1S,nmtvr,1,nmtvr,'hprimes')
#ifdef MP
      endif
#endif
!

      deallocate (tmps)
#undef HPRIMES
#ifdef MP
      deallocate (grid1)
#endif
!
      call rzerout
! output the initial condition----sig file---
      if (.not.lrestart) then
! sig
#ifdef MP
      if( mype.eq.master ) then
#endif
        print*,'output initial condition'
      open(unit=nrsmo1, file='r_sigf00', action='write',               &
     &     form='unformatted', status='replace')
#ifdef MP
      endif
#endif
      call swrite(nrsmo1,label,fhour,idate,                              &
     &               flat,flon,fm2,fm2x,fm2y,sl,si,                      &
     &               gz,qm ,tem ,uum ,vvm                                &
#ifdef NDSL
     &                    ,rmgrd                                       &
#else
     &                    ,rm                                          &
#endif
#ifdef NONHYD
     &                    ,pnm,tnm,onm                                 &
#endif
     &                    )
      if(ncld.gt.0) then
      if ((nmmiph .eq. 2) .or.       & ! Zhao-Carr
     &    (nmmiph .eq. 3) .or.       & ! Brad-Ferrier
#ifdef NEW_WSM3
     &    (nmmiph .eq. 4) .or.       & ! WSM3
#endif
#ifdef NEW_WSM6
     &    (nmmiph .eq. 6) .or.       & ! WSM6
#endif
     &    (nmmiph .eq. 11) .or.      & ! GFDL MP v1
     &    (nmmiph .eq. 12) .or.      & ! GFDL MP v2
     &    (nmmiph .eq. 13) .or.      & ! GFDL MP v3
     &    (nmmiph .eq. 15)) then       ! GCE 3ICE
              
! iflag=1 for output f_ice f_rain f_rimsf
        iflag=1
        call swrite2(nrsmo1,phy_f3d,phy_f2d,iflag,.false.)
      endif
      endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      close(nrsmo1)

! sfc
      open(unit=nrflop,file='r_sfcf00',                                  &   
     &        status='unknown',form='unformatted') 
#ifdef MP
      endif
#endif
      call rfixio(fhour,idate,                                           &  
     &   tsea,smc,sheleg,stc,tg3,zorl,                                   &
     &   cv,cvb,cvt,slmsk,f10m,canopy,vfrac,vtype,stype,                 &
     &   shdmin,shdmax,slope,snoalb,                                     &
     &   uustar,ffmm,ffhh,alvsf,alvwf,alnsf,alnwf,facsf,facwf,           &
     &   1,nrflop,rsfcsec,lsfcmrg,.false.)
#ifdef MP
      if( mype.eq.master ) then
#endif
      close(nrflop)
#ifdef MP
      endif
#endif
      endif
! ---end output initial file
! ...............prepare ken points .................
#ifdef RKN
       isave=0
       itnum=0
       if(npoint.ne.0) then
         isave=1
         itnum=1
         call rkenpre(flat,flon,gzzb,nrflip)
       endif
#endif

      return
      end
