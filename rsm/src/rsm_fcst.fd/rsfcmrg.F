!yj: don't merge clima data:  (2019.Sep.3)
!yj:     tg3,zorl,alb*4(alvsf,alvwf,alnsf,alnwf),
!yj:     vfrac,vtype,stype,alf*2(facsf,facwf)
       subroutine rsfcmrg(label,fhour,idate,                             &  
     &    tsea,smc,sheleg,stc,                                           &
     &    cv,cvb,cvt,slmsk,f10m,canopy,                                  &
     &    uustar,ffmm,ffhh                                    )
!yj: ori       
!yj:       subroutine rsfcmrg(label,fhour,idate,                             &  
!yj:     &    tsea,smc,sheleg,stc,tg3,zorl,                                  &
!yj:     &    cv,cvb,cvt,slmsk,f10m,canopy,vfrac,vtype,stype,                &
!yj:     &    uustar,ffmm,ffhh,alvsf,alvwf,alnsf,alnwf,facsf,facwf)
!
! usage:    call rsfcmrg(fhour,tsea,smc,sheleg,stc,tg3,zorl,
!          &    cv,cvb,cvt,albed,slmsk,f10m,canopy,hprime,
!          &    ioflag,nread,nunit,shour)
!   input argument list: or output argument list: depends on ioflag
!     fhour     - forecast hour
!     tsea      - land and sea surface temperature (k)
!     smc       - two layers of soil moisture contents (0.47 - 0.1)
!     sheleg    - snow depth (cm)
!     stc       - two layers of soil temperature (k)
!     tg3       - the lowest soil temperature at 3m (k)
!     zorl      - surface roughness
!     plantr    - planetary resistence
!     cv        - cloud amount
!     cvb       - cloud base (sigma layer number)
!     cvt       - cloud top (sigma layer number)
!     albed     - albedo
!     slmsk     - sea land mask
!     f10m      - 10 m height factor for wind
!     canopy    - surface canopy
!
!$$$
      use machine , only : kind_io4
!yj----use module_rinpbase, only: sfc-field for rfixio----  
      use module_rinpbase, only:             &
     &   tsea_i=>tsea,sheleg_i=>sheleg       &
     &  ,tg3_i=>tg3clm,zorl_i=>zorclm        &
     &  ,slmsk_i=>slmsk  &
     &  ,vfrac_i=>vegclm  &
     &  ,vtype_i=>vetclm,stype_i=>sotclm &
     &  ,albedo_i=>albclm,alfclm
!yj -------------- use module rinpbase --------------- end
#include <paramodel.h>
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
#define LNGRDS lngrdp
#else
#define LNGRDS lngrd
#endif
#ifdef G2R
#define IIMAX lonf
#define JJMAX latg
#else
#define IIMAX cigrd1
#define JJMAX cjgrd1
#endif

      parameter (imax=IIMAX, jmax=JJMAX, ijmax=imax*jmax)

      real(kind=kind_io4) ghours,fhours
      real(kind=kind_io4) ggs(ijmax),ggls(ijmax,lsoil)
      real(kind=kind_io4) gg2s(ijmax,2),gg4s(ijmax,4)
      real(kind=kind_io4) rgs(LNGRDS),rgls(LNGRDS,lsoil)
      real(kind=kind_io4) rg2s(LNGRDS,2),rg4s(LNGRDS,4)
!cc jun write out
      real(kind=kind_io4) grgs(lngrd),grgls(lngrd,lsoil)
      real(kind=kind_io4) grg2s(lngrd,2),grg4s(lngrd,4)
      real(kind=kind_io4) rrgs(lngrd),rrgls(lngrd,lsoil)
      real(kind=kind_io4) rrg2s(lngrd,2),rrg4s(lngrd,4)
!cc
      dimension slg(ijmax),ggl(ijmax,lsoil)                              &  
     &         , gg(ijmax),gg2(ijmax,2)                                  &
     &         ,gg4(ijmax,4)
      dimension ggsea(ijmax),ggland(ijmax),iwork(ijmax)
      dimension fgrgs(lngrd), fgrgls(lngrd,lsoil)                        &
     &         ,fgrg2s(lngrd,2), fgrg4s(lngrd,4)                         
      dimension  gclat(jmax)                                             &
     &  , colrad(latg2), wgt(latg2),wgtcs(latg2),rcs2(latg2)
      dimension idate(4)
      character*8 label1(4), label(4)
      dimension gtsea  (LNGRDS), gsmc   (LNGRDS,lsoil),                  &
     &          gsheleg(LNGRDS), gstc   (LNGRDS,lsoil),                  &
     &          gtg3   (LNGRDS), gvfrac (LNGRDS),                        &
     &          gzorl  (LNGRDS), gslmsk (LNGRDS),                        &
     &          gcv    (LNGRDS), gcvb   (LNGRDS),                        &
     &          gcvt   (LNGRDS),                                         &
     &          gf10m  (LNGRDS), gcanopy(LNGRDS),                        &
     &          gvtype (LNGRDS),                                         &
     &          gstype (LNGRDS), guustar(LNGRDS),                        &
     &          gffmm  (LNGRDS),   gffhh(LNGRDS),                        &
     &          galvsf (LNGRDS),  galvwf(LNGRDS),                        &
     &          galnsf (LNGRDS),  galnwf(LNGRDS),                        &
     &          gfacsf (LNGRDS),  gfacwf(LNGRDS)
      dimension rr     (LNGRDS)
      dimension tsea  (LNGRDS), smc   (LNGRDS,lsoil),                    &
     &          sheleg(LNGRDS), stc   (LNGRDS,lsoil),                    &
     &          tg3   (LNGRDS), vfrac (LNGRDS),                          &
     &          zorl  (LNGRDS), slmsk (LNGRDS),                          &
     &          cv    (LNGRDS), cvb   (LNGRDS),                          &
     &          cvt   (LNGRDS),                                          &
     &          f10m  (LNGRDS), canopy(LNGRDS),                          &
     &          vtype (LNGRDS),                                          &
     &          stype (LNGRDS), uustar(LNGRDS),                          &
     &          ffmm  (LNGRDS),   ffhh(LNGRDS),                          &
     &          alvsf (LNGRDS),  alvwf(LNGRDS),                          &
     &          alnsf (LNGRDS),  alnwf(LNGRDS),                          &
     &          facsf (LNGRDS),  facwf(LNGRDS)
!     &          ,rslmsk (LNGRDS),
      integer,save :: version
      data version/199802/
      dimension slmk(LNGRDS),simk(LNGRDS)                                &
     &         ,flat(LNGRDS),flon(LNGRDS)
      character cfhour*16,cform*40,filename*100
!yj      data  nsfc/12/
!yj      save nsfc

!cc
#ifdef MP
      if( mype.eq.master ) then
#endif
      kh=nint(fhour)
      ndig=max(log10(kh+0.5)+1.,2.)
      write(cform,'("(i",i1,".",i1,")")') ndig,ndig
      write(cfhour,cform) kh
!hmhj open(unit=nsfc,file='rb_sfcf'//cfhour                              &  
!hmhj&    ,status='old',form='unformatted')
      write(filename,"(a,a)") 'rb_sfcf',cfhour
!yj      open(unit=nsfc,file=filename                                       &  
!yj     &    ,status='old',form='unformatted')
!yj      print *,'in s2rinp_mrg, nsfc=',nsfc,'open rb_sfcf=',cfhour
      print *,'in s2rinp_mrg, replace open rb_sfcf=',cfhour,'with mpmd'
!
!yj      rewind nsfc
!
! get slmask first
!yj      read(nsfc) label1
!yj      read(nsfc) ghours,idate
!yj      ghour=ghours
!yj      read(nsfc)               ! tsea
!yj      read(nsfc)               ! smc
!yj      read(nsfc) ggs            ! sheleg
!yj      call sgl2ful(ggs,gg,ijmax)
#ifdef MP
      endif
       gg=sheleg_i
       call rmpbcastr(gg,ijmax)
#endif
#ifdef G2R
      call gg2rg_mpi(gg,gsheleg,1)
#else
      call cg2rg_mpi(gg,gsheleg,1)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
!yj      read(nsfc)               ! stc
!yj      read(nsfc)               ! tg3
!yj      read(nsfc)               ! zorl
!yj      read(nsfc)               ! cv
!yj      read(nsfc)               ! cvb
!yj      read(nsfc)               ! cvt
!yj      read(nsfc)               ! alvsf,alvwf,alnsf,alnwf
!yj      read(nsfc) ggs          ! slmsk
!yj      call sgl2ful(ggs,slg,ijmax)
! get sea and ice index
      do 10 i=1,ijmax
      gg(i) = 0.0
      slg=slmsk_i
      if( slg(i) .eq. 2.0 ) gg(i) = 2.0
 10   continue
#ifdef MP
       endif
       call rmpbcastr(gg,ijmax)
#endif
#ifdef G2R
       call gg2rg_mpi(gg,simk,1)
#else
       call cg2rg_mpi(gg,simk,1)
#endif
      do 20 i=1,LNGRDS
      if( simk(i) .gt. 0.9 ) then
          simk(i) = 2.0
      else
          simk(i) = 0.0
      endif
 20   continue

#ifdef MP
      if( mype.eq.master ) then
#endif
! get sea and land index
      do 30 i=1,ijmax
      gg(i) = 0.0
      if( slg(i) .eq. 1.0 ) gg(i) = 1.0
 30   continue
#ifdef MP
       endif
       call rmpbcastr(gg,ijmax)
#endif
#ifdef G2R
       call gg2rg_mpi(gg,slmk,1)
#else
       call cg2rg_mpi(gg,slmk,1)
#endif
        
      do 40 i=1,LNGRDS
      if( slmk(i) .gt. 0.4 ) then
          slmk(i) = 1.0
      else
          slmk(i) = 0.0
      endif
 40   continue
!       do 45 i=1,LNGRDS
!       gslmsk(i) = max( simk(i), slmk(i) )
!   45  continue
        do 45 i=1,LNGRDS
! add sea ice info
        slmsk(i) = max(slmsk(i),simk(i))
 45     continue
!
! put sea ice point into sea point
#ifdef MP
      if( mype.eq.master ) then
#endif
        do i=1,ijmax 
        if(slg(i).eq.2.0) slg(i)=0.0
        enddo
#ifdef MP
       endif
#endif
!
#ifdef MP
      if( mype.eq.master ) then
#endif
!yj      rewind nsfc

!yj      read(nsfc) label1
! hour idate ----- 2
!yj      read(nsfc) ghours,idate
! tsea ----------- 3  -----merge from Ocean 
!yj      read(nsfc) ggs
!yj      call sgl2ful(ggs,gg,ijmax)
      gg=tsea_i
      call maxmin(gg,ijmax,1,1,1,'tsea')
      call expnd(gg,slg,0,ggsea,iwork,imax,jmax)
      call expnd(gg,slg,1,ggland,iwork,imax,jmax)
#ifdef MP
      endif
      call rmpbcastr(ggsea,ijmax)
      call rmpbcastr(ggland,ijmax)
#endif
#ifdef G2R
      call gg2rg_mpi(ggsea,gtsea,1)
      call gg2rg_mpi(ggland,rr,1)
#else
      call cg2rg_mpi(ggsea,gtsea,1)
      call cg2rg_mpi(ggland,rr,1)
#endif
      do i=1,LNGRDS
         if (slmsk(i).eq.1.0) gtsea(i)=rr(i)
      enddo
      do 50 i=1,LNGRDS
      if( slmsk(i) .eq. 2.0e0 ) then
          gtsea(i)= min( gtsea(i), 271.20e0 )
      else if( slmsk(i).eq.0.0e0 ) then
          gtsea(i)= max( gtsea(i), 271.21e0 )
      else if( slmsk(i).eq.1.0e0 .and.                                &  
     &         gsheleg(i).gt.0.0e0 ) then
          gtsea(i)= min( gtsea(i), 273.16e0 )
      endif
 50   continue

      do 101 i=1,LNGRDS
        if( slmsk(i) .eq. 0.0 ) tsea(i)= gtsea(i)
 101  continue
!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
! soil moisture ----- 4
!yj      read(nsfc) ggls
!!!      call sgl2ful(ggls,ggl,lsoil*ijmax)
!yj#ifdef MP
!yj      endif
!!!      call rmpbcastr(ggl,lsoil*ijmax)
!yj#endif
!!!#ifdef G2R
!!!      call gg2rg_mpi(ggl,gsmc,lsoil)
!!!#else
!!!      call cg2rg_mpi(ggl,gsmc,lsoil)
!!!#endif
!
! snow depth ----- 5
!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
!yj      read(nsfc) ggs
!!!      call sgl2ful(ggs,gg,ijmax)
!!!      call maxmin(gg,ijmax,1,1,1,'snow depth')
!yj#ifdef MP
!yj      endif
!!!      call rmpbcastr(gg,ijmax)
!yj#endif
!!!#ifdef G2R
!!!      call gg2rg_mpi(gg,gsheleg,1)
!!!#else
!!!      call cg2rg_mpi(gg,gsheleg,1)
!!!#endif
!
! stc ------------ 6
!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
!yj       read(nsfc) ggls
!!!      call sgl2ful(ggls,ggl,lsoil*ijmax)
!yj#ifdef MP
!yj      endif
!!!      call rmpbcastr(ggl,lsoil*ijmax)
!yj#endif
!!!#ifdef G2R
!!!      call gg2rg_mpi(ggl,gstc,lsoil)
!!!#else
!!!      call cg2rg_mpi(ggl,gstc,lsoil)
!!!#endif

! tg3 ------------ 7  ---- Merge all from global cycle ----
!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!!yj      read(nsfc)  ggs
!!yj      call sgl2ful(ggs,gg,ijmax)
!      gg=tg3_i
!      call expnd(gg,slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg,slg,1,ggland,iwork,imax,jmax)
!#ifdef MP
!      endif
!      call rmpbcastr(ggsea,ijmax)
!      call rmpbcastr(ggland,ijmax)
!#endif
!#ifdef G2R
!      call gg2rg_mpi(ggsea,gtg3,1)
!      call gg2rg_mpi(ggland,rr,1)
!#else
!      call cg2rg_mpi(ggsea,gtg3,1)
!      call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) gtg3(i)=rr(i)
!      enddo
!      do 60 i=1,LNGRDS
!        tg3(i)= gtg3(i)
! 60   continue
! roughness ----- 8 ---- merge ovre land ----
!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!!yj      read(nsfc) ggs
!!yj      call sgl2ful(ggs,gg,ijmax)
!      gg=zorl_i
!      call expnd(gg,slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg,slg,1,ggland,iwork,imax,jmax)
!      call maxmin(gg,ijmax,1,1,1,'roughness')
!#ifdef MP
!      endif
!      call rmpbcastr(ggsea,ijmax)
!      call rmpbcastr(ggland,ijmax)
!#endif
!#ifdef G2R
!      call gg2rg_mpi(ggsea,gzorl,1)
!      call gg2rg_mpi(ggland,rr,1)
!#else
!      call cg2rg_mpi(ggsea,gzorl,1)
!      call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) gzorl(i)=rr(i) 
!      enddo
!      do 70 i=1,LNGRDS
!        if( slmsk(i) .ge. 1.0 ) zorl(i)= gzorl(i)
! 70   continue
!
!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
! cloud amount -- 9
!yj      read(nsfc) ggs
!!!      call sgl2ful(ggs,gg,ijmax)
!yj#ifdef MP
!yj      endif
!!!      call rmpbcastr(gg,ijmax)
!yj#endif
!!!#ifdef G2R
!!!      call gg2rg_mpi(gg,gcv,1)
!!!#else
!!!      call cg2rg_mpi(gg,gcv,1)
!!!#endif

!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
! cloud base --- 10
!yj      read(nsfc) ggs
!!!      call sgl2ful(ggs,gg,ijmax)
!yj#ifdef MP
!yj      endif
!!!      call rmpbcastr(gg,ijmax)
!yj#endif
!!!#ifdef G2R
!!!      call gg2rg_mpi(gg,gcvb,1)
!!!#else
!!!      call cg2rg_mpi(gg,gcvb,1)
!!!#endif

!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
! cloud top ---- 11
!yj      read(nsfc) ggs
!!!      call sgl2ful(ggs,gg,ijmax)
!!!      call maxmin(gg,ijmax,1,1,1,'cloud top')
!yj#ifdef MP
!yj      endif
!!!      call rmpbcastr(gg,ijmax)
!yj#endif
!!!#ifdef G2R
!!!      call gg2rg_mpi(gg,gcvt,1)
!!!#else
!!!      call cg2rg_mpi(gg,gcvt,1)
!!!#endif

! 4 albedo ------- 12
!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!!yj      read(nsfc) gg4s
!!yj      call sgl2ful(gg4s,gg4,ijmax*4)
!      gg4=albedo_i
!      call expnd(gg4(1,1),slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg4(1,1),slg,1,ggland,iwork,imax,jmax)
!#ifdef MP
!      endif
!      call rmpbcastr(ggsea,ijmax)
!      call rmpbcastr(ggland,ijmax)
!#endif
!#ifdef G2R
!      call gg2rg_mpi(ggsea,galvsf,1)
!      call gg2rg_mpi(ggland,rr,1)
!#else
!      call cg2rg_mpi(ggsea,galvsf,1)
!      call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) galvsf(i)=rr(i)
!      enddo
!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!      call expnd(gg4(1,2),slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg4(1,2),slg,1,ggland,iwork,imax,jmax)
!#ifdef MP
!      endif
!      call rmpbcastr(ggsea,ijmax)
!      call rmpbcastr(ggland,ijmax)
!#endif
!#ifdef G2R
!      call gg2rg_mpi(ggsea,galvwf,1)
!      call gg2rg_mpi(ggland,rr,1)
!#else
!      call cg2rg_mpi(ggsea,galvwf,1)
!      call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) galvwf(i)=rr(i)
!      enddo
!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!      call expnd(gg4(1,3),slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg4(1,3),slg,1,ggland,iwork,imax,jmax)
!#ifdef MP
!      endif
!      call rmpbcastr(ggsea,ijmax)
!      call rmpbcastr(ggland,ijmax)
!#endif
!#ifdef G2R
!      call gg2rg_mpi(ggsea,galnsf,1)
!      call gg2rg_mpi(ggland,rr,1)
!#else
!      call cg2rg_mpi(ggsea,galnsf,1)
!      call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) galnsf(i)=rr(i)
!      enddo
!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!      call expnd(gg4(1,4),slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg4(1,4),slg,1,ggland,iwork,imax,jmax)
!#ifdef MP
!      endif
!      call rmpbcastr(ggsea,ijmax)
!      call rmpbcastr(ggland,ijmax)
!!#endif
!#ifdef G2R
!      call gg2rg_mpi(ggsea,galnwf,1)
!      call gg2rg_mpi(ggland,rr,1)
!#else
!      call cg2rg_mpi(ggsea,galnwf,1)
!      call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) galnwf(i)=rr(i)
!      enddo
!cjun
!
!
!!  for r sfc
!
!      do 2103 i=1,LNGRDS
!          alvsf(i)= galvsf(i)
!2103     continue
!      do 2104 i=1,LNGRDS
!          alvwf(i)= galvwf(i)
!2104     continue
!      do 2105 i=1,LNGRDS
!          alnsf(i)= galnsf(i)
!2105     continue
!      do 2106 i=1,LNGRDS
!          alnwf(i)= galnwf(i)
!2106     continue

!
!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
! sea land mask -- 13
!yj      read(nsfc) ggs
!yj#ifdef MP
!yj      endif
!yj#endif

!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!! vegitation fraction ---- 14
!!yj      read(nsfc) ggs
!!yj      call sgl2ful(ggs,gg,ijmax)
!      gg=vfrac_i
!      call expnd(gg,slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg,slg,1,ggland,iwork,imax,jmax)
!      call maxmin(gg,ijmax,1,1,1,'vegitation fraction')
!#ifdef MP
!      endif
!      call rmpbcastr(ggsea,ijmax)
!      call rmpbcastr(ggland,ijmax)
!#endif
!#ifdef G2R
!      call gg2rg_mpi(ggsea,gvfrac,1)
!      call gg2rg_mpi(ggland,rr,1)
!#else
!      call cg2rg_mpi(ggsea,gvfrac,1)
!      call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) gvfrac(i)=rr(i)
!      enddo
!      do 103 i=1,LNGRDS
!          if( slmsk(i) .ge. 1.0 ) vfrac(i)= gvfrac(i)
!103     continue
!! prepare 15 to end for end of file condition
!      do i=1,LNGRDS
!#ifdef doNOAHLSM
!        gcanopy(i)=60.0
!#else
!        gcanopy(i)=0.6
!#endif
!        gf10m(i)=0.8
!        gvtype(i)=1.0
!        gstype(i)=1.0
!        gfacsf(i)=1.0
!        gfacwf(i)=1.0
!        guustar(i) = 1.
!        gffmm(i) = log(30.)
!        gffhh(i) = log(30.)
!      enddo
!!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
! vegitation canopy ------ 15
!yj      read(nsfc,end=1500) ggs
!!!      call sgl2ful(ggs,gg,ijmax)
!!!      call expnd(gg,slg,0,ggsea,iwork,imax,jmax)
!!!      call expnd(gg,slg,1,ggland,iwork,imax,jmax)
!yj#ifdef MP
!yj      endif
!!!      call rmpbcastr(ggsea,ijmax)
!!!      call rmpbcastr(ggland,ijmax)
!yj#endif
!!!#ifdef G2R
!!!      call gg2rg_mpi(ggsea,gcanopy,1)     
!!!     call gg2rg_mpi(ggland,rr,1)
!!!#else
!!!      call cg2rg_mpi(ggsea,gcanopy,1)
!!!      call cg2rg_mpi(ggland,rr,1)
!!!#endif
!!!      do i=1,LNGRDS
!!!         if (slmsk(i).eq.1.0) gcanopy(i)=rr(i)    
!!!      enddo
!!!        do 104 i=1,LNGRDS
!!!          if( slmsk(i) .ge. 1.0 ) canopy(i)= gcanopy(i)
!!!104     continue
!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
! f10m -------------- 16
!yj      read(nsfc,end=1500)ggs
!!!      call sgl2ful(ggs,gg,ijmax)
!yj#ifdef MP
!yj      endif
!!!      call rmpbcastr(gg,ijmax)
!yj#endif
!!!#ifdef G2R
!!!      call gg2rg_mpi(gg,gf10m,1)
!!!#else
!!!      call cg2rg_mpi(gg,gf10m,1)
!!!#endif

!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!! vegitation type ---- 17 ---merge over land ----
!!yj      read(nsfc,end=1500) ggs
!!yj      call sgl2ful(ggs,gg,ijmax)
!      gg=vtype_i
!      call expnd(gg,slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg,slg,1,ggland,iwork,imax,jmax)
!#ifdef MP
!       endif
!      call rmpbcastr(ggsea,ijmax)
!      call rmpbcastr(ggland,ijmax)
!#endif
!#ifdef G2R
!       call gg2rg_mpi(ggsea,gvtype,1)
!       call gg2rg_mpi(ggland,rr,1)
!#else
!       call cg2rg_mpi(ggsea,gvtype,1)
!       call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) gvtype(i)=rr(i)     
!      enddo
!!cjun
!
!!      call ful2sgl(vtype,rgs,LNGRDS)
!!      call sgl2ful(rgs,vtype,LNGRDS)
!        do 105 i=1,LNGRDS
!          if( slmsk(i) .ge. 1.0 ) vtype(i)= gvtype(i)
!105     continue
!!
!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!! soil type ---------- 18 ----merge over land ----
!!yj      read(nsfc,end=1500) ggs
!!yj      call sgl2ful(ggs,gg,ijmax)
!      gg=stype_i
!      call expnd(gg,slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg,slg,1,ggland,iwork,imax,jmax)
!      call maxmin(gg,ijmax,1,1,1,'soil type 18')
!#ifdef MP
!       endif
!       call rmpbcastr(ggsea,ijmax)
!       call rmpbcastr(ggland,ijmax)
!#endif
!#ifdef G2R
!       call gg2rg_mpi(ggsea,gstype,1)
!       call gg2rg_mpi(ggland,rr,1)
!#else
!       call cg2rg_mpi(ggsea,gstype,1)
!       call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) gstype(i)=rr(i)      
!      enddo
!
!!      call ful2sgl(stype,rgs,LNGRDS)
!!      call sgl2ful(rgs,stype,LNGRDS)
!        do 106 i=1,LNGRDS
!          if( slmsk(i) .ge. 1.0 ) stype(i)= gstype(i)
!106     continue
!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!! zenith angle dependent factor ---- 19 ----merge onver land ----
!!yj      read(nsfc,end=1500) gg2s
!!yj      call sgl2ful(gg2s,gg2,ijmax*2)
!      gg2=alfclm
!      call expnd(gg2(1,1),slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg2(1,1),slg,1,ggland,iwork,imax,jmax)
!      call maxmin(gg,ijmax*2,1,1,1,'zenith angle')
!#ifdef MP
!       endif
!       call rmpbcastr(ggsea,ijmax)
!       call rmpbcastr(ggland,ijmax)
!#endif
!#ifdef G2R
!      call gg2rg_mpi(ggsea,gfacsf,1)
!      call gg2rg_mpi(ggland,rr,1)
!#else
!      call cg2rg_mpi(ggsea,gfacsf,1)
!      call cg2rg_mpi(ggland,rr,1)
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) gfacsf(i)=rr(i)
!      enddo
!#ifdef MP
!      if( mype.eq.master ) then
!#endif
!      call expnd(gg2(1,2),slg,0,ggsea,iwork,imax,jmax)
!      call expnd(gg2(1,2),slg,1,ggland,iwork,imax,jmax)
!#ifdef MP
!       endif
!       call rmpbcastr(ggsea,ijmax)   
!       call rmpbcastr(ggland,ijmax)   
!#endif
!#ifdef G2R
!      call gg2rg_mpi(ggsea,gfacwf,1)   
!      call gg2rg_mpi(ggland,rr,1)          
!#else
!      call cg2rg_mpi(ggsea,gfacwf,1)   
!      call cg2rg_mpi(ggland,rr,1)          
!#endif
!      do i=1,LNGRDS
!         if (slmsk(i).eq.1.0) gfacwf(i)=rr(i)
!      enddo
!
!!  for r sfc
!      do 107 i=1,LNGRDS
!        if( slmsk(i) .ge. 1.0 ) facsf(i)= gfacsf(i)
!        if( slmsk(i) .ge. 1.0 ) facwf(i)= gfacwf(i)
!107   continue
!!yj1500  continue
!!yj#ifdef MP
!!yj      if( mype.eq.master ) then
!!yj#endif
! ustar ------------ 20
!yj      read(nsfc,end=2000) ggs
!!!      call sgl2ful(ggs,gg,ijmax)
!!!      call maxmin(gg,ijmax,1,1,1,'ustar')
!yj#ifdef MP
!yj       endif
!!!       call rmpbcastr(gg,ijmax)
!yj#endif
!!!#ifdef G2R
!!!       call gg2rg_mpi(gg,guustar,1)
!!!#else
!!!       call cg2rg_mpi(gg,guustar,1)
!!!#endif
!
!
!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
!!!      call maxmin(gg,ijmax,1,1,1,'ustar')
! momentum exchange similarity function -- 21
!yj      read(nsfc,end=2000) ggs
!!!      call sgl2ful(ggs,gg,ijmax)
!!!      call maxmin(gg,ijmax,1,1,1,'momentum function')
!yj#ifdef MP
!yj       endif
!!!       call rmpbcastr(gg,ijmax)
!yj#endif
!!!#ifdef G2R
!!!       call gg2rg_mpi(gg,gffmm,1)
!!!#else
!!!       call cg2rg_mpi(gg,gffmm,1)
!!!#endif
!cjun

!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
! heat exchange similarity function -- 22
!yj      read(nsfc,end=2000) ggs
!!!      call sgl2ful(ggs,gg,ijmax)
!!!      call maxmin(gg,ijmax,1,1,1,'heat exchange function')
!yj#ifdef MP
!yj       endif
!!!       call rmpbcastr(gg,ijmax)
!yj#endif
!!!#ifdef G2R
!!!       call gg2rg_mpi(gg,gffhh,1)
!!!#else
!!!       call cg2rg_mpi(gg,gffhh,1)
!!!#endif


!yj 2000  continue
!yj#ifdef MP
!yj      if( mype.eq.master ) then
!yj#endif
!yj       close(nsfc)
!yj#ifdef MP
!yj       endif
!yj#endif

       return
       end
