      subroutine rfixio(fhour,idate,                                     &  
!yj      subroutine rfixio(label,fhour,idate,                            &  
     &    tsea,smc,sheleg,stc,tg3,zorl,                                  &
     &    cv,cvb,cvt,slmsk,f10m,canopy,vfrac,vtype,stype,                &
     &    shdmin,shdmax,slope,snoalb,                                    &
     &    uustar,ffmm,ffhh,alvsf,alvwf,alnsf,alnwf,facsf,facwf,          &
     &    ioflag,nunit,rsfcsec,lsfcmrg,iorestartfile)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rfixio
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  do regional surface input and output
!
! program history log:
!
! usage:    call rfixio(fhour,tsea,smc,sheleg,stc,tg3,zorl,plantr,
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
!     hprime    - mountain variance
!     ioflag    - 0 : for read, so all of above as output
!                 1 : for write, so all of above as input
!     shour     - second for the current sequential integration
!
!
!   input files:
!     nread
!
!   output files:
!     nunit
!
!   subprograms called: 
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
      use machine , only : kind_io4,kind_io8
!yj----use module_rinpinit, only: sfc-field for rfixio----  
!     use module_rinpinit, only:                        &
!    &   tsea_i=>tsea,smc_i=>smc,sheleg_i=>sheleg       &
!    &  ,stc_i=>stc,tg3_i=>tg3clm,zorl_i=>zorclm        &
!    &  ,cv_i=>cv,cvb_i=>cvb,cvt_i=>cvt,slmsk_i=>slmsk  &
!    &  ,f10m_i=>f10m,canopy_i=>canopy,vfrac_i=>vegclm  &
!    &  ,vtype_i=>vetclm,stype_i=>sotclm,uustar_i=>uustar &
!    &  ,shdmin_i=>vmnclm,shdmax_i=>vmxclm              &
!    &  ,slope_i=>slpclm,snoalb_i=>absclm               &
!    &  ,ffmm_i=>ffmm,ffhh_i=>ffhh                      &
!    &  ,albclm,alfclm
!yj     &  ,alvsf_i=>albclm(:,1) &
!yj     &  ,alvwf_i=>albclm(:,2),alnsf_i=>albclm(:,3)      &
!yj     &  ,alnwf_i=>albclm(:,4)   &
!yj     &  ,facsf_i=>alfclm(:,1),facwf_i=>alfclm(:,2)
!yj -------------- use module rinpinit --------------- end
!yj------------------------------------------------------
!yj-note: rfixio in rsmini_mpi get these vars. from -----
!yj------ module_rinpinit, and rfixio in rsmsav do  -----
!yj------ not get vars from module_rinpinit, since  -----
!yj------ they are deallocate in rsmini_mpi.F. In   -----
!yj------ rsmsav.F, vars are input argument.        -----
!yj------------------------------------------------------
#include <paramodel.h>
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
      real, allocatable :: spec(:)
      real, allocatable :: spec2(:)
      real, allocatable :: spec3(:)
      real, allocatable :: grid(:,:)
#define LNGRDS lngrdp
#else
#define LNGRDS lngrd
#endif
!cc
      real(kind=kind_io4) fhours
      real(kind=kind_io4), allocatable :: tmps(:)
      real(kind=kind_io4), allocatable :: tmpsl(:,:)
      real(kind=kind_io4), allocatable :: tmps2(:,:)
      real(kind=kind_io4), allocatable :: tmps4(:,:)
      real(kind=kind_io8) fhoursd
      real(kind=kind_io8), allocatable :: tmpsd(:)
      real(kind=kind_io8), allocatable :: tmpsld(:,:)
      real(kind=kind_io8), allocatable :: tmps2d(:,:)
      real(kind=kind_io8), allocatable :: tmps4d(:,:)
      dimension idate(4)
      character*8 label(4)
      dimension tsea  (LNGRDS), smc   (LNGRDS,lsoil),                    &  
     &          sheleg(LNGRDS), stc   (LNGRDS,lsoil),                    &
     &          tg3   (LNGRDS), vfrac (LNGRDS),                          &
     &          zorl  (LNGRDS), slmsk (LNGRDS),                          &
     &          cv    (LNGRDS), cvb   (LNGRDS),                          &
     &          cvt   (LNGRDS),                                          &
     &          f10m  (LNGRDS), canopy(LNGRDS),                          &
     &          vtype (LNGRDS),                                          &
     &          stype (LNGRDS), uustar(LNGRDS),                          &
     &         shdmin (LNGRDS), shdmax(LNGRDS),                          &
     &          slope (LNGRDS), snoalb(LNGRDS),                          &
     &          ffmm  (LNGRDS),   ffhh(LNGRDS),                          &
     &          alvsf (LNGRDS),  alvwf(LNGRDS),                          &
     &          alnsf (LNGRDS),  alnwf(LNGRDS),                          &
     &          facsf (LNGRDS),  facwf(LNGRDS)
      integer,save :: version
      data version/199802/
      logical lsfcmrg
      logical iorestartfile
#ifdef MP
      allocate (spec(lngrd))
      allocate (spec2(lngrd))
      allocate (spec3(lngrd))
      allocate (grid(lngrd,lsoil))
#endif
      if (.not.iorestartfile) then
      allocate (tmps(lngrd))
      allocate (tmpsl(lngrd,lsoil))
      allocate (tmps2(lngrd,2))
      allocate (tmps4(lngrd,4))
      else
      allocate (tmpsd(lngrd))
      allocate (tmpsld(lngrd,lsoil))
      allocate (tmps2d(lngrd,2))
      allocate (tmps4d(lngrd,4))
      endif
!yj note: nunit is sfc initial data if ioflag = 0      
!
!  ioflag = 0  ...  read fixed field from unit nunit
!  ioflag = 1  ...  write fixed field to unit nunit
!
      if(ioflag.eq.0) then
!
! start read --------------------------------------
!
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        rewind nunit
        read(nunit) label
        read(nunit) fhoursd, idate
        fhour=fhoursd
      endif
99    format(1h ,'fhour, idate=',f6.2,2x,4(1x,i4))
      print *,'fix field read in from unit=',nunit
      print 99,fhour, idate
#ifdef MP
      endif
!yj      call rmpbcastc(label,32)
      call rmpbcastr(fhour,1)
      call rmpbcasti(idate,4)
#endif
!--------------------------------------------1.tsea
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       tsea_i=spec
        call maxmin(spec,lngrd,1,1,1,' tsea')
!     else
!       spec=tsea_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,tsea,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,tsea,lngrd)
#endif
!--------------------------------------------2.smc
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsld
!       call sgl2ful(tmpsl,grid, lsoil*lngrd)
        grid=tmpsld
!       smc_i=grid
        do k=1,lsoil
          call maxmin(grid(1,k),lngrd,1,1,1,' smc')
        enddo
!     else
!       grid=smc_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(grid,igrd1,jgrd1,smc,igrd1p,jgrd1p,lsoil)
      endif
#else
!     call sgl2ful(tmpsl,smc,lsoil*lngrd)
#endif
!--------------------------------------------3.sheleg
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       sheleg_i=spec
        call maxmin(spec,lngrd,1,1,1,' sheleg')
!     else
!       spec=sheleg_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,sheleg,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,sheleg,lngrd)
#endif
!--------------------------------------------4.stc
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsld
!       call sgl2ful(tmpsl,grid, lsoil*lngrd)
        grid=tmpsld
!       stc_i=grid
        do k =1,lsoil
          call maxmin(grid(1,k),lngrd,1,1,1,' stc')
        enddo
!     else
!       grid=stc_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
        call rmpgf2p(grid,igrd1,jgrd1,stc,igrd1p,jgrd1p,lsoil)
      endif
#else
!     call sgl2ful(tmpsl,stc,lsoil*lngrd)
#endif
!--------------------------------------------5.tg3
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       tg3_i=spec
        call maxmin(spec,lngrd,1,1,1,' tg3')
!     else
!       spec=tg3_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
        call rmpgf2p(spec,igrd1,jgrd1,tg3,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,tg3,lngrd)
#endif
!--------------------------------------------6.zorl
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       zorl_i=spec
        call maxmin(spec,lngrd,1,1,1,' zorl')
!     else
!       spec=zorl_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
        call rmpgf2p(spec,igrd1,jgrd1,zorl,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,zorl,lngrd)
#endif
!--------------------------------------------7.cv
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       cv_i=spec
        call maxmin(spec,lngrd,1,1,1,' cv')
!     else
!       spec=cv_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,cv,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,cv,lngrd)
#endif
!--------------------------------------------8.cvb
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       cvb_i=spec
        call maxmin(spec,lngrd,1,1,1,' cvb')
!     else
!       spec=cvb_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,cvb,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,cvb,lngrd)
#endif
!--------------------------------------------9.cvt
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       cvt_i=spec
        call maxmin(spec,lngrd,1,1,1,' cvt')
!     else
!       spec=cvt_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,cvt,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,cvt,lngrd)
#endif
!--------------------------------------------10.albedo
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmps4d
!       call sgl2ful(tmps4(1,1),spec, lngrd)
        spec=tmps4d(:,1)
!       albclm(:,1)=spec
        call maxmin(spec,lngrd,1,1,1,' alvsf')
!     else
!       spec=albclm(:,1)
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,alvsf,igrd1p,jgrd1p,1)
      endif
      if( mype.eq.master ) then
      if (iorestartfile) then
!       call sgl2ful(tmps4(1,2),spec, lngrd)
        spec=tmps4d(:,2)
!       albclm(:,2)=spec
        call maxmin(spec,lngrd,1,1,1,' alvwf')
!     else
!       spec=albclm(:,2)
!       print*,'-skip-'
      endif
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,alvwf,igrd1p,jgrd1p,1)
      endif
      if( mype.eq.master ) then
      if (iorestartfile) then
!        call sgl2ful(tmps4(1,3),spec, lngrd)
         spec=tmps4d(:,3)
!        albclm(:,3)=spec
         call maxmin(spec,lngrd,1,1,1,' alnsf')
!     else
!        spec=albclm(:,3)
!       print*,'-skip-'
      endif
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,alnsf,igrd1p,jgrd1p,1)
      endif
      if( mype.eq.master ) then
      if (iorestartfile) then
!        call sgl2ful(tmps4(1,4),spec, lngrd)
         spec=tmps4d(:,4)
!        albclm(:,4)=spec
         call maxmin(spec,lngrd,1,1,1,' alnwf')
!     else
!        spec=albclm(:,4)
!       print*,'-skip-'
      endif
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,alnwf,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps4(1,1),alvsf,lngrd)
!     call sgl2ful(tmps4(1,2),alvwf,lngrd)
!     call sgl2ful(tmps4(1,3),alnsf,lngrd)
!     call sgl2ful(tmps4(1,4),alnwf,lngrd)
#endif
!--------------------------------------------11.slmsk
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       slmsk_i=spec
        call maxmin(spec,lngrd,1,1,1,' slmsk')
!     else
!       spec=slmsk_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,slmsk,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,slmsk,lngrd)
#endif
!--------------------------------------------12.vfrac
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       vfrac_i=spec
        call maxmin(spec,lngrd,1,1,1,' vfrac')
!     else
!       spec=vfrac_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,vfrac,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,vfrac,lngrd)
#endif
!--------------------------------------------13.canopy
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,err=5000)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       canopy_i=spec
        call maxmin(spec,lngrd,1,1,1,' canopy')
!     else
!       spec=canopy_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,canopy,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,canopy,lngrd)
#endif
!--------------------------------------------14.f10m
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,err=5000)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       f10m_i=spec
        call maxmin(spec,lngrd,1,1,1,' f10m')
!     else
!       spec=f10m_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,f10m,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,f10m,lngrd)
#endif
!--------------------------------------------15.vtype
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,err=5000)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       vtype_i=spec
        call maxmin(spec,lngrd,1,1,1,' vtype')
!     else
!       spec=vtype_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,vtype,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,vtype,lngrd)
#endif
!--------------------------------------------16.stype
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,err=5000)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       stype_i=spec
        call maxmin(spec,lngrd,1,1,1,' stype')
!     else
!       spec=stype_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,stype,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,stype,lngrd)
#endif
!--------------------------------------------17.facswf
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,err=5000)  tmps2d
!       call sgl2ful(tmps2(1,1),spec, lngrd)
        spec=tmps2d(:,1)
!       alfclm(:,1)=spec
        call maxmin(spec,lngrd,1,1,1,' facsf')
!     else
!       spec=alfclm(:,1)
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,facsf,igrd1p,jgrd1p,1)
      endif
      if( mype.eq.master ) then
      if (iorestartfile) then
!       call sgl2ful(tmps2(1,2),spec, lngrd)
        spec=tmps2d(:,2)
!       alfclm(:,2)=spec
        call maxmin(spec,lngrd,1,1,1,' facwf')
!     else
!       spec=alfclm(:,2)
!       print*,'-skip-'
      endif
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,facwf,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps2(1,1),facsf,lngrd)
!     call sgl2ful(tmps2(1,2),facwf,lngrd)
#endif
!--------------------------------
#ifdef MP
      do ij = 1, lngrd
        spec(ij) = 1.
        spec2(ij) = log(30.)
        spec3(ij) = log(30.)
      enddo
#else
      do ij = 1, lngrd
        uustar(ij) = 1.
        ffmm(ij) = log(30.)
        ffhh(ij) = log(30.)
      enddo
#endif
!--------------------------------------------18.uustar
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,end=200)   tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       uustar_i=spec
        call maxmin(spec,lngrd,1,1,1,' uustar')
!     else
!       spec=uustar_i
!       print*,'-skip-'
      endif
#ifdef MP
 200  continue
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,uustar,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,uustar,lngrd)
 200  continue
#endif
!--------------------------------------------19.ffmm
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,end=201)  tmpsd
!       call sgl2ful(tmps,spec2, lngrd)
        spec2=tmpsd
!       ffmm_i=spec2
        call maxmin(spec2,lngrd,1,1,1,' ffmm')
!     else
!       spec2=ffmm_i
!       print*,'-skip-'
      endif
#ifdef MP
 201  continue
      endif
      if (iorestartfile) then
      call rmpgf2p(spec2,igrd1,jgrd1,ffmm,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,ffmm,lngrd)
 201  continue
#endif
!--------------------------------------------20.ffhh
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,end=202)  tmpsd
!       call sgl2ful(tmps,spec3, lngrd)
        spec3=tmpsd
!       ffhh_i=spec3
        call maxmin(spec3,lngrd,1,1,1,' ffhh')
!     else
!       spec3=ffhh_i
!       print*,'-skip-'
      endif
#ifdef MP
 202  continue
      endif
      if (iorestartfile) then
      call rmpgf2p(spec3,igrd1,jgrd1,ffhh,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,ffhh,lngrd)
 202  continue
#endif
!
!--------------------------------------------21.shdmin
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,end=5000)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       shdmin_i=spec
        call maxmin(spec,lngrd,1,1,1,' shdmin')
!     else
!       spec=shdmin_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,shdmin,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,shdmin,lngrd)
#endif
!
!--------------------------------------------22.shdmax
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,end=5000)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       shdmax_i=spec
        call maxmin(spec,lngrd,1,1,1,' shdmax')
!     else
!       spec=shdmax_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,shdmax,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,shdmax,lngrd)
#endif
!
!--------------------------------------------23.slopetyp
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,end=5000)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       slope_i=spec
        call maxmin(spec,lngrd,1,1,1,' slopetyp')
!     else
!       spec=slope_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,slope,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,slope,lngrd)
#endif
!
!--------------------------------------------24.snoalb
#ifdef MP
      if( mype.eq.master ) then
#endif
      if (iorestartfile) then
        read(nunit,end=5000)  tmpsd
!       call sgl2ful(tmps,spec, lngrd)
        spec=tmpsd
!       snoalb_i=spec
        call maxmin(spec,lngrd,1,1,1,' snoalb')
!     else
!       spec=snoalb_i
!       print*,'-skip-'
      endif
#ifdef MP
      endif
      if (iorestartfile) then
      call rmpgf2p(spec,igrd1,jgrd1,snoalb,igrd1p,jgrd1p,1)
      endif
#else
!     call sgl2ful(tmps,snoalb,lngrd)
#endif
!      
#ifdef MP
      if( mype.eq.master ) then
#endif
      print *,'in rfixio, rsfc read finished!'
#ifdef MP
      endif
#endif

      else
!yj-note: below, if ioflag =1
!  start to write -----------------------------------------------
!
!ccc jun change to add merge sfc
!
      kfhours=nint(fhour*3600)
      kdays=mod(kfhours,nint(rsfcsec))      
      if ( lsfcmrg .and. kdays.eq.0 .and. fhour .gt.0) then
#ifdef MP
      if( mype.eq.master ) then
#endif
      print *,'rsfcmrg write from unit=',nunit,                          &  
     &      ' fhour,idate=',fhour,idate,                                 &
     &      'already updated when getting gfs data in mainrsmf_mrg'
#ifdef MP
      endif
#endif
!yj---add rinpbase sfc data---start
!  kfhours= 21600 kdays= 21600 rsfcsec= 86400.00000000000
!  so rsfcmrg will not be called in the mpmd version.
!yj---add rinpbase sfc data---end
!
!yj: don't merge clima data: (2019.Sep.3)
!yj:     tg3,zorl,alb*4(alvsf,alvwf,alnsf,alnwf),
!yj:     vfrac,vtype,stype,alf*2(facsf,facwf)
!          call rsfcmrg(label,fhour,idate,                               &   
!     &   tsea,smc,sheleg,stc,                                            &
!     &   cv,cvb,cvt,slmsk,f10m,canopy,                                   &
!     &   uustar,ffmm,ffhh                                    )
!yj: ori       
!yj:           call rsfcmrg(label,fhour,idate,                               &   
!yj:      &   tsea,smc,sheleg,stc,tg3,zorl,                                   &
!yj:      &   cv,cvb,cvt,slmsk,f10m,canopy,vfrac,vtype,stype,                 &
!yj:      &   uustar,ffmm,ffhh,alvsf,alvwf,alnsf,alnwf,facsf,facwf)

      endif
! endif ( lsfcmrg .and. kdays.eq.0 .and. fhour .gt.0) above      
#ifdef MP
      if( mype.eq.master ) then
!yj      print *,'rsfcmrg finished'
#endif
      rewind nunit
      label(1)='ncep'
      label(2)='rsm'
      label(3)='mpi'
      label(4)='version'
      write(nunit) label
      if (iorestartfile) then
        fhoursd=fhour
        write(nunit) fhoursd, idate, igrd1, jgrd1, version
      else
        fhours=fhour
        write(nunit) fhours, idate, igrd1, jgrd1, version
      endif
#ifdef MP
      endif
#endif
!
!1.tsea
#ifdef MP
      call rmpgp2f(tsea,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(tsea,  tmps,lngrd)
      write(nunit) tmps
#endif
!2.smc
#ifdef MP
      call rmpgp2f(smc,igrd1p,jgrd1p,grid,igrd1,jgrd1,lsoil)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsld=grid
        write(nunit) tmpsld
      else
!       call ful2sgl(grid,  tmpsl,lsoil*lngrd)
        tmpsl=grid
        write(nunit) tmpsl
      endif
      endif
#else
      call ful2sgl(smc,   tmpsl,lsoil*lngrd)
      write(nunit) tmpsl
#endif
!3.sheleg
#ifdef MP
      call rmpgp2f(sheleg,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(sheleg,tmps,lngrd)
      write(nunit) tmps
#endif
!4.stc
#ifdef MP
      call rmpgp2f(stc,igrd1p,jgrd1p,grid,igrd1,jgrd1,lsoil)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsld=grid
        write(nunit) tmpsld
      else
!       call ful2sgl(grid,  tmpsl,lsoil*lngrd)
        tmpsl=grid
        write(nunit) tmpsl
      endif
      endif
#else
      call ful2sgl(stc,   tmpsl,lsoil*lngrd)
      write(nunit) tmpsl
#endif
!5.tg3
#ifdef MP
      call rmpgp2f(tg3,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(tg3,   tmps,lngrd)
      write(nunit) tmps
#endif
!6.zorl
#ifdef MP
      call rmpgp2f(zorl,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(zorl,  tmps,lngrd)
      write(nunit) tmps
#endif
!7.cv
#ifdef MP
      call rmpgp2f(cv,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(cv,    tmps,lngrd)
      write(nunit) tmps
#endif
!8.cvb
#ifdef MP
      call rmpgp2f(cvb,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(cvb,   tmps,lngrd)
      write(nunit) tmps
#endif
!9.cvt
#ifdef MP
      call rmpgp2f(cvt,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(cvt,   tmps,lngrd)
      write(nunit) tmps
#endif
!10.albedo
#ifdef MP
      call rmpgp2f(alvsf,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmps4d(:,1)=spec
      else
!       call ful2sgl(spec,  tmps4(1,1),lngrd)
        tmps4(:,1)=spec
      endif
      endif
      call rmpgp2f(alvwf,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmps4d(:,2)=spec
      else
!       call ful2sgl(spec,  tmps4(1,2),lngrd)
        tmps4(:,2)=spec
      endif
      endif
      call rmpgp2f(alnsf,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmps4d(:,3)=spec
      else
!       call ful2sgl(spec,  tmps4(1,3),lngrd)
        tmps4(:,3)=spec
      endif
      endif
      call rmpgp2f(alnwf,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmps4d(:,4)=spec
        write(nunit) tmps4d
      else
!       call ful2sgl(spec,  tmps4(1,4),lngrd)
        tmps4(:,4)=spec
        write(nunit) tmps4
      endif
      endif
#else
      call ful2sgl(alvsf, tmps4(1,1),lngrd)
      call ful2sgl(alvwf, tmps4(1,2),lngrd)
      call ful2sgl(alnsf, tmps4(1,3),lngrd)
      call ful2sgl(alnwf, tmps4(1,4),lngrd)
      write(nunit) tmps4
#endif
!11.slmsk
#ifdef MP
      call rmpgp2f(slmsk,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(slmsk, tmps,lngrd)
      write(nunit) tmps
#endif
!12.vfrac
#ifdef MP
      call rmpgp2f(vfrac,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(vfrac, tmps,lngrd)
      write(nunit) tmps
#endif
!13.canopy
#ifdef MP
      call rmpgp2f(canopy,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(canopy,tmps,lngrd)
      write(nunit) tmps
#endif
!14.f10m
#ifdef MP
      call rmpgp2f(f10m,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(f10m,  tmps,lngrd)
      write(nunit) tmps
#endif
!15.vtype
#ifdef MP
      call rmpgp2f(vtype,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(vtype, tmps,lngrd)
      write(nunit) tmps
#endif
!16.stype
#ifdef MP
      call rmpgp2f(stype,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(stype, tmps,lngrd)
      write(nunit) tmps
#endif
!17.facswf
#ifdef MP
      call rmpgp2f(facsf,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmps2d(:,1)=spec
      else
!       call ful2sgl(spec,  tmps2(1,1),lngrd)
        tmps2(:,1)=spec
      endif
      endif
      call rmpgp2f(facwf,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmps2d(:,2)=spec
        write(nunit) tmps2d
      else
!       call ful2sgl(spec,  tmps2(1,2),lngrd)
        tmps2(:,2)=spec
        write(nunit) tmps2
      endif
      endif
#else
      call ful2sgl(facsf, tmps2(1,1),lngrd)
      call ful2sgl(facwf, tmps2(1,2),lngrd)
      write(nunit) tmps2
#endif
!18.uustar
#ifdef MP
      call rmpgp2f(uustar,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(uustar,tmps,lngrd)
      write(nunit) tmps
#endif
!19.ffmm
#ifdef MP
      call rmpgp2f(ffmm,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(ffmm,  tmps,lngrd)
      write(nunit) tmps
#endif
!20.ffhh
#ifdef MP
      call rmpgp2f(ffhh,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(ffhh,  tmps,lngrd)
      write(nunit) tmps
#endif
!      
!21.shdmin
#ifdef MP
      call rmpgp2f(shdmin,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(shdmin, tmps,lngrd)
      write(nunit) tmps
#endif
!
!22.shdmax
#ifdef MP
      call rmpgp2f(shdmax,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(shdmax, tmps,lngrd)
      write(nunit) tmps
#endif
!
!23.slope
#ifdef MP
      call rmpgp2f(slope,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(slope, tmps,lngrd)
      write(nunit) tmps
#endif
!
!24.snoalb
#ifdef MP
      call rmpgp2f(snoalb,igrd1p,jgrd1p,spec,igrd1,jgrd1,1)
      if( mype.eq.master ) then
      if (iorestartfile) then
        tmpsd=spec
        write(nunit) tmpsd
      else
!       call ful2sgl(spec,  tmps,lngrd)
        tmps=spec
        write(nunit) tmps
      endif
      endif
#else
      call ful2sgl(snoalb, tmps,lngrd)
      write(nunit) tmps
#endif
!
#ifdef MP
      if( mype.eq.master ) then
#endif
      close(nunit)
      print *,'rfixio write to unit=',nunit,' fhour,idate=',fhour,idate
#ifdef MP
      endif
#endif
!----end if(ioflag)----
      endif
!
#ifdef MP
      deallocate (spec)
      deallocate (spec2)
      deallocate (spec3)
      deallocate (grid)
#endif
      if (.not.iorestartfile) then
      deallocate (tmps)
      deallocate (tmpsl)
      deallocate (tmps2)
      deallocate (tmps4)
      else
      deallocate (tmpsd)
      deallocate (tmpsld)
      deallocate (tmps2d)
      deallocate (tmps4d)
      endif
!
      return
5000  print *,' error read in rfixio '
      call abort
      end
