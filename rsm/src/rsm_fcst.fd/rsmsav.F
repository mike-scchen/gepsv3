#include <defopt.h>
      subroutine rsmsav
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rsmsav
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  do regional sigma, flux, and surface files output
!
! program history log:
!
! usage:    call rsmsav
!   input argument list:
!       gfhour  - global model forecast hour
!
!   output argument list:
!
!   common block:
!       rscomf  comcon
!
!   input files: none
!
!   output files: none
!
!   subprograms called:
!       swrite  rfixio
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
      use machine , only : kind_io4
!     use module_rinpbase, only: gzo,qo,teo,uuo,vvo,rqo
!ccj

#ifdef wind_fitch
      use module_wind_fitch , only :writ_out_fitch_rsm, wival, wjval, power, fi_nt, spdh
#endif

#include <paramodel.h>
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
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
#include <rscommap.h>
#ifdef RKN
#include <rscomgpd.h>
#endif
#include <defdim.h>

#ifdef wind_fitch
      integer nn,fn
      real spd(fi_nt),pp(fi_nt)
#endif
!cjun
      dimension tmp(lngrd,levr),tmpk1(lngrd)
      character cfhour*16,cform*40
      integer iflag
      character*20 stringn
      character*100 filename
!cc
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,' ---- into rwrtres write restart files  ---- '
!
      print 800,deltim,fhour
 800  format(1h0,'dt=',e10.2,2x,'fcst hour =',f10.2,' h')
!yj------------change name of output------------
!yj      write(stringn,"(A10,I4.4,I2.2,I2.2,I2.2)") &
!yj     &  'rename.sh ',idate(4),idate(2),idate(3),idate(1)
      kh=nint(fhour)
      ndig=max(log10(kh+0.5)+1.,2.)
      write(cform,'("(i",i1,".",i1,")")') ndig,ndig
      write(cfhour,cform) kh

      write(filename,"(a,a)") 'r_sigf',cfhour
      open(unit=nrsmo1, file=filename, action='write', &
     &     form='unformatted', status='replace')
      open(unit=nrsms1, file='r_sigf', action='write', &
     &     form='unformatted', status='replace')
      open(unit=nrsms2, file='r_sigftdt', action='write', &
     &     form='unformatted', status='replace')
      open(unit=nb, file='r_sigbase', action='write', &
     &       form='unformatted', status='replace')

!yj idate:  (4)year, (2)month, (3)date, (1)hour     
!yj      print*,'stringn for rename.sh is:',stringn
!      print*,'*****start call call system shell in rsmsav.F*****'
!yj      call system(stringn)
!      call system('rename.sh '//cfhour)
!      print*,'*****finish call call system shell in rsmsav.F*****'
!yj-note-  !check with rgetcom_mrg.F   
!yj!ln -fs r_sig.${idate}   fort.70
!yj!ln -fs r_sfc.${idate}   fort.72
!yj!ln -fs r_flx.${idate}   fort.73
!yj-----------change name of output------------
#ifdef MP
      endif
#endif
!yj: basefield updated to be the same as fhour at output timestep
!     call rsm2bgd(gzo,qo,teo,uuo,vvo,rqo)
!     do 644 n=1,lngrdb
!       qb(n)= qltb(n)
!644  continue
!     do 655 k=1,levr
!     do 655 n=1,lngrdb
!       teb(n,k)= teltb(n,k)
!       uub(n,k)= uultb(n,k)
!       vvb(n,k)= vvltb(n,k)
!655  continue
!     do 656 nc=1,ntotal
!     do 656 k=1,levr
!     do 656 i=1,lngrdb
!       rqb(i,k,nc)= rqltb(i,k,nc)
!656  continue
!
!..................................................................
      call swrite(nrsmo1,label,fhour,idate,                              &  
     &            flat,flon,fm2,fm2x,fm2y,sl,si,                         &
     &            gz,qm,tem,uum,vvm                                      &
#ifdef NDSL
     &               ,rmgrd                                              &
#else
     &               ,rm                                                 &
#endif
#ifdef NONHYD
     &                     ,pnm,tnm,onm                                  &
#endif
     &                                )

#ifdef MP
      call rmpsynall
#endif

!      call swrite(nrsmo2,label,fhour,idate,                              &
!     &            flat,flon,fm2,fm2x,fm2y,sl,si,                         &
!     &            gz,q ,te ,uu ,vv                                       &
!#ifdef NDSL
!     &               ,rqgrd                                              &
!#else
!     &               ,rq                                                 &
!#endif
!#ifdef NONHYD
!     &                     ,pn ,tn ,on                                   &
!#endif
!     &                                 )
!
!#ifdef MP
!      call rmpsynall
!#endif
      call swritewav(nrsms1,fhour,idate,                                 &  
     &            flat,flon,fm2,fm2x,fm2y,sl,si,                         &
     &            gz,qm,tem,uum,vvm                                      &
#ifdef NDSL
     &               ,rmgrd                                              &
#else
     &               ,rm                                                 &
#endif
#ifdef NONHYD
     &                     ,pnm,tnm,onm,tov                              &
#endif
     &                                )

#ifdef MP
      call rmpsynall
#endif

      call swritewav(nrsms2,fhour,idate,                                 &
     &            flat,flon,fm2,fm2x,fm2y,sl,si,                         &
     &            gz,q ,te ,uu ,vv                                       &
#ifdef NDSL
     &               ,rqgrd                                              &
#else
     &               ,rq                                                 &
#endif
#ifdef NONHYD
     &                     ,pn ,tn ,on,tov                               &
#endif
     &                                 )

#ifdef MP
      call rmpsynall
#endif

!ccj
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
     &    (nmmiph .eq. 15) ) then      ! GCE 3ICE
              
! iflag=1 for output f_ice f_rain f_rimsf
        iflag=1
        call swrite2(nrsmo1,phy_f3d,phy_f2d,iflag,.false.)
#ifdef MP
        call rmpsynall
#endif
        iflag=0
        call swrite2(nrsms2,phy_f3d,phy_f2d,iflag,.true.)
#ifdef MP
        call rmpsynall
#endif
      endif
      endif

      call writebase(nb,qb,teb,uub,vvb,rqb,lngrdb,levr,ntotal)
#ifdef MP
      call rmpsynall
#endif

!
#ifdef MP
      if( mype.eq.master ) then
#endif
      close(nrsmo1)
      close(nrsms1)
      close(nrsms2)
      close(nb)
      write(filename,"(a,a)") 'r_sfcf',cfhour
      open(unit=nrflop,file=filename,                                    &   
     &        status='unknown',form='unformatted') 
#ifdef MP
      endif
#endif
!cc
!
!sfc
!................................................
!yj-note: ioflag=1, the arguments are input and are
!yj-note: passed by rscomgrd.h here!!! 
#ifdef MP
      if( mype.eq.master ) print*,'call rfixio.F in rsmsav.F ' 
#endif
      call rfixio(fhour,idate,                                     &
     &   tsea,smc,sheleg,stc,tg3,zorl,                                   &
     &   cv,cvb,cvt,slmsk,f10m,canopy,vfrac,vtype,stype,                 &
     &   shdmin,shdmax,slope,snoalb,                                     &
     &   uustar,ffmm,ffhh,alvsf,alvwf,alnsf,alnwf,facsf,facwf,           &
     &   1,nrflop,rsfcsec,lsfcmrg,.false.)
#ifdef MP
      call rmpsynall
#endif
      call rfixio(fhour,idate,                                     &
     &   tsea,smc,sheleg,stc,tg3,zorl,                                   &
     &   cv,cvb,cvt,slmsk,f10m,canopy,vfrac,vtype,stype,                 &
     &   shdmin,shdmax,slope,snoalb,                                     &
     &   uustar,ffmm,ffhh,alvsf,alvwf,alnsf,alnwf,facsf,facwf,           &
     &   1,nrflop8,rsfcsec,.false.,.true.)
      ! don't do lsfcmrg twice
#ifdef MP
      call rmpsynall
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      close(nrflop)
      close(nrflop8)

#ifdef MP
      endif
#endif
! flx
      secphy=(fhour-zhour)*3600.
      secswr=max(secphy,fhswr*3600.)
      seclwr=max(secphy,fhlwr*3600.)
      nrflxbs=73

#ifdef MP
      if( mype.eq.master ) then
#endif
      write(filename,"(a,a)") 'r_flxf',cfhour
      open(unit=nrflxbs,file=filename,                                    &   
     &        status='unknown',form='unformatted')
#ifdef MP
      endif
#endif
!
      call rwrtsfc(zhour,fhour,idate,nrflxbs,slmsk,rlat,rlon,            &  
     &            dusfc,dvsfc,dtsfc,dqsfc,tsea,smc,stc,sheleg,           &
     &            secswr,seclwr,dlwsfc,ulwsfc,                           &
     &            geshem,bengsh,gflux,fluxr,                             &
     &            u10m,v10m,gust,t2m,q2m,psurf,zorl,tmpmax,tmpmin,       &
     &            spfhmax,spfhmin,wvuflx,wvvflx,                         &
     &            srunoff,soilm,snwdph,                                  &
     &            runoff,ep,cldwrk,cflash,dugwd,dvgwd,hpbl,pwat,         &
     &            icen,igen,icen2,ienst,iensi)

#ifdef wind_fitch
      if (lsphy) then
       pp = 0.0
      spd = 0.0
      do nn=1,fi_nt
      if ( (wival(nn).ne.-9999) .and. (wjval(nn).ne.-9999) )  then
       pp(nn)=power(nn)/3600.*deltim/1000.
       spd(nn)=spdh(nn)/3600.*deltim
      endif
      enddo
      call writ_out_fitch_rsm(spd,pp)
      power=0.0
      spdh=0.0
      endif
#endif
#ifdef MP
      call rmpsynall
#endif

#ifdef MP
      if( mype.eq.master ) then
#endif
      close(nrflxbs)
#ifdef MP
      endif
#endif
!..............regional ken points ...................
#ifdef RKN
      if( npoint.gt.0 ) then
#ifdef MP
      if( mype.eq.master ) then
#endif
        open(unit=nrkenp,file='r_ken_final',                             &   
     &           status='unknown',form='unformatted')
#ifdef MP
      endif
#endif
        call rkenwrt(nrkenp,label,fhour,idate,si,sl)
#ifdef MP
      if( mype.eq.master ) then
#endif
        close(nrkenp)
#ifdef MP
        call rmpsynall
      endif
#endif
      endif
#endif
!................................................
      return
      end
