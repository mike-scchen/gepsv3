!#define NEW_WSM3
!#define NEW_WSM6
#define fix_neg_q
! ---------------------------------------------------------------------!
! 2019.Sep.01 Y.J.CHEN: add noah lsm (lsm=0: OSU lsm; lsm=1: Noah lsm) !
!                       lsm=0: ivegsrc=2, isot=0                       !
!                       lsm=1: ivegsrc=1, isot=1                       !
! ---------------------------------------------------------------------!
#include <define.h>
      subroutine rbphys(dtp,dtf,xkt2,                                    &   
#ifndef NONHYD
     & plamgrs,pphigrs,                                                  &
#endif
     & ugrs,vgrs,pgr,tgrs,qgrs,                                          &
#ifndef NONHYD
     & xgrs,                                                             &
#endif
     & gt0,gq0,gu0,gv0,                                                  &
#ifdef NONHYD
     & pgrs,wgrs,gp0,gw0,tgrh,gh0,                                       &
#endif

#ifdef RAS
     & rannum,flipv, ncrnd,                                              &
#endif
#ifdef RKN
     & psexp,radsl,                                                      &
     & dlwsf1,qss,                                                       &
     & gflx,cd,                                                          &
     & cdq,evap,                                                         &
     & hflx,drain,                                                       &
     & runof,cld1d,                                                      &
     & gamt,gamq,                                                        &
     & dqsfc1,dtsfc1,                                                    &
     & dusfc1,dvsfc1,                                                    &
     & dusfcg,dvsfcg,                                                    &
     & rainc,rainl,                                                      &
     & hsw,vvel,                                                         &
     & stsoil,smsoil,                                                    &
#endif
     & cldtot,                                                           &
     & phy_f3dv,phy_f2dv, me,                                            &
!    & vvel,phy_f3dv,phy_f2dv, me,                                            &
     & lonlens,igrd1s,                                                   &
     & lat)
!x   & lat,ptotj,pwatj)
!fpp$ noconcur r
!
      use machine , only : kind_phys
      use sigmap_cons, only : sig=>sigma
      use funcphys , only : fpkap
      use physcons, rocp => con_rocp, fv => con_fvirt                    & 
     &,             grav => con_g, cp => con_cp, rd => con_rd            &
     &,             rv => con_rv, hvap => con_hvap,  hfus => con_hfus    &
     &,             rerth => con_rerth                                   &
     &,             qcmin => con_qcmin                                   &
#ifdef NONHYD
     &,             ccv => con_cv                                        &
#endif
     &,             pi => con_pi                                         &
     &,             cvap => con_cvap, t0c => con_t0c, rdorv => con_eps   &
     &,             rhoair0 => con_rhoair0, rhoh2o => con_rhoh2o         &
     &,             cliq => con_cliq, cice => con_csol, psat => con_psat
! for WSM3
#ifdef NEW_WSM3
      use rmodule_mp_wsm3,        only: wsm3
#else
      use rmodule_mp_wsm3
#endif
! for WSM5
      use module_mp_wsm5
! for WSM6
#ifdef NEW_WSM6
      use rmodule_mp_wsm6,        only: wsm6
#else
      use rmodule_mp_wsm6
#endif
! for GFDL MP v1
      use module_mp_gfdl_v16, only: gfdl_cloud_microphys_driver,        &
                                    cloud_diagnosis
! for GFDL MP v2
      use rmodule_mp_gfdl_v2,  only: gfdlv2_driver => gfdl_cld_mp_driver
! for GFDL MP v3
      use module_mp_gfdl_v3,  only: gfdlv3_driver => gfdl_cld_mp_driver
! for Goddard (GCE) 3ICE MP
!      use module_mp_gsfcgce,  only: gsfcgce
      use module_mp_gsfcgce_3ice_nuwrf,  only: gsfcgce_3ice_nuwrf
#ifdef wind_fitch
      use module_wind_fitch , only : dragforce, fi_nt, wjval, wi
#endif     
!
!ccj
#include <paramodel.h>
#include <rscomgrd.h>
#include <rscomver.h>
#include <rscomio.h>
#include <rscomloc.h>
!
#include <rscommap.h>
#include <rscomgpd.h>
#include <rscombgt.h>
!
      dimension  xkt2(igrd1s)
!
      integer ncw,nrw,niw,nsw,ngl,ntcw,ntrw,ntiw,ntsw,ntgl,lat,jcapr
!cc
!...................................................................
! warning: gt0,gq0,gu0,gv0 may overlay tgrs,qgrs,ugrs,vgrs respectively.
#ifndef NONHYD
      dimension  plamgrs(igrd1s)
      dimension  pphigrs(igrd1s)
#endif
      dimension     ugrs(igrd1s,levr)
      dimension     vgrs(igrd1s,levr)
      dimension      pgr(igrd1s)
      dimension     tgrs(igrd1s,levr)
      dimension     qgrs(igrd1s,levh)
      real, allocatable,dimension(:,:,:) :: qgrsx, dqdtx
#ifndef NONHYD
      dimension     xgrs(igrd1s,levr)
#endif
      dimension      gt0(igrd1s,levr)
      dimension      gq0(igrd1s,levh)
      dimension      gu0(igrd1s,levr)
      dimension      gv0(igrd1s,levr)
#ifdef NONHYD
      dimension     pgrs(igrd1s,levr)
      dimension     wgrs(igrd1s,levr+1)
      dimension      gp0(igrd1s,levr)
      dimension      gw0(igrd1s,levr+1)
      dimension     tgrh(igrd1s,levr)
      dimension     gh0(igrd1s,levr)
#endif
      dimension      sio(igrd1s,levr+1)
      dimension      slo(igrd1s,levr)
      dimension      delo(igrd1s,levr+1)
      dimension      sik(igrd1s,levr+1)
      dimension      slk(igrd1s,levr)
!ccj
      dimension  phy_f3dv(igrd1s,levr,num_p3d),                          &
     &           phy_f2dv(igrd1s,num_p2d)
!c
      real(kind=kind_phys) prsl(igrd1s,levr),  prsi(igrd1s,levr+1)       &   
     &,    prslk(igrd1s,levr)                                            &
     &,    prsik(igrd1s,levr+1)                                          &
     &,    phii(igrd1s,levr+1), phil(igrd1s,levr)                        &
     &,    delp(igrd1s,levr)
!c
!...................................................................
      real,parameter ::    epsq=1.e-12, hsub=hvap+hfus
      dimension kbot(igrd1s),ktop(igrd1s),kuo(igrd1s),                   &   
     &dudt(igrd1s,levr),dvdt(igrd1s,levr),                               &
     &dtdt(igrd1s,levr),dqdt(igrd1s,levh),                               &
#ifdef NONHYD
     & dwdt(igrd1s,levr+1),                                              &
#endif
     & hsw(igrd1s,levr),gflx(igrd1s),                                    &
     & rain(igrd1s),rainc(igrd1s),rainl(igrd1s),rain1(igrd1s),           &
     & rain2(igrd1s),evapc(igrd1s),wind(igrd1s),kpbl(igrd1s),            &
     & w2(igrd1s,levr),vvel(igrd1s,levr)
!xb110> 
#ifdef NTDK
      real(kind=kind_phys) rphil(igrd1s,levr),rdudt(igrd1s,levr)          &
     &     ,rdvdt(igrd1s,levr),rdtdt(igrd1s,levr),rdqdt(igrd1s,levr)      &
     &     ,rgu0(igrd1s,levr),rgv0(igrd1s,levr),rgt0(igrd1s,levr)         &
     &     ,rgq0(igrd1s,levr),rvvel(igrd1s,levr)                          &
     &     ,rprsl(igrd1s,levr),rprsi(igrd1s,levr+1)                       &
     &     ,rclw(igrd1s,levr,2)
      real(kind=kind_phys) rkbot(igrd1s),rktop(igrd1s)
#endif
!xb110<
! for wsm
      real                 ::  ww(igrd1s,levr)
      real                 ::  zl(igrd1s,levr)
      real                 ::  zi(igrd1s,levr+1)
      real                 ::  delz(igrd1s,levr)
      real                 ::  prsl_pa(igrd1s,levr)
      real                 ::  rainncv(igrd1s)
      real                 ::  snow1(igrd1s,levr)
      real                 ::  graupel1(igrd1s,levr)
      real                 ::  snowncv(igrd1s,levr)
      real                 ::  graupelncv(igrd1s,levr)
! for moninedmf & YSU
      real                 ::  prsi_pa(igrd1s,levr+1)
      real                 ::  delp_pa(igrd1s,levr)
#ifdef LTN
!lightning
      real flash(igrd1s)   !flash density (unit in flashes km^-2 day^-1)
      real snow_flx(igrd1s,levr),ptu(igrd1s,levr),pqu(igrd1s,levr)
      real ud_ptu(igrd1s,levr),ud_pqu(igrd1s,levr),                      &
           ud_ztenh(igrd1s,levr),ud_zqenh(igrd1s,levr),                  &
           ud_cnvw(igrd1s,levr),ud_rho(igrd1s,levr),                     &
           ud_phil(igrd1s,levr),ud_prsl_hpa(igrd1s,levr),                &
           ud_snow_flx(igrd1s,levr)

      integer ud_kbot(igrd1s),ud_ktop(igrd1s)
!      
!.................................................................
#endif
      dimension  psexp(igrd1s),                                          &
     &          snowmt(igrd1s),                                          &
     &             fm1(igrd1s),   fh1(igrd1s),                           &
     &              cd(igrd1s),   cdq(igrd1s),                           &
     &             qss(igrd1s), radsl(igrd1s)
!for noahlsm >     
      dimension  psexp_pa(igrd1s),snwdph_mm(igrd1s)
!for noahlsm <     
      dimension dusfcg(igrd1s),dvsfcg(igrd1s)
      dimension dusfc1(igrd1s),dvsfc1(igrd1s),                           &
     &          dtsfc1(igrd1s),dqsfc1(igrd1s),                           &
     &          dlwsf1(igrd1s),ulwsf1(igrd1s)
      dimension smsoil(igrd1s,lsoil),stsoil(igrd1s,lsoil)
      dimension slsoil(igrd1s,lsoil)
      integer soiltyp(igrd1s),vegtype(igrd1s)
      dimension sigmaf(igrd1s)
      dimension rb(igrd1s),rhscnpy(igrd1s)
      dimension ai(igrd1s,lsoil),bi(igrd1s,lsoil),drain(igrd1s)
      dimension cci(igrd1s,lsoil),rhsmc(igrd1s,lsoil),runof(igrd1s)
      dimension zsoil(igrd1s,lsoil),cld1d(igrd1s)
      dimension evap(igrd1s),hflx(igrd1s),stress(igrd1s)
!c    dimension evap(igrd1s),hflx(igrd1s),rnet(igrd1s)
      dimension t850(igrd1s),uv925(igrd1s),uv850(igrd1s),               &
                ep1d(igrd1s),uvdiff(igrd1s)
      dimension gamt(igrd1s),gamq(igrd1s)
      dimension ustar1(igrd1s)
!c
      integer lonfequiv
      dimension var(igrd1s),oc(igrd1s),oa4(igrd1s,4),clx(igrd1s,4)
      real(kind=kind_phys) theta(igrd1s),gamma(igrd1s)
      real(kind=kind_phys) sigma(igrd1s),elvmax(igrd1s)
      real(kind=kind_phys) qr_col(igrd1s,levr), fc_ice(igrd1s,levr)
      real(kind=kind_phys) work1(igrd1s),work2(igrd1s)
      real(kind=kind_phys) rhc(igrd1s,levr),clstp
      real(kind=kind_phys) sr(igrd1s),xncw(igrd1s)
      real(kind=kind_phys) work3(igrd1s,levr),ptem,ptem1
      real(kind=kind_phys) qmin, rhzbot, rhztop
      real(kind=kind_phys) rhbbot, rhbtop, rhpbl
      real(kind=kind_phys) dxmax, dxmin, dxinv, cb2mb, tempv
!yj-implement noah land model---------------------------------------->>>
!input arguments, files .................
      integer  lsm,ivegsrc,isot
!old:           lsm=0: ivegsrc=2, isot=0 
!new:           lsm=1: ivegsrc=1, isot=1 
#ifdef doNOAHLSM
      dimension qsurf(igrd1s)
      parameter(lsm=1, ivegsrc=1, isot=1)
#else
#ifdef doPROGTM
      parameter(lsm=0, ivegsrc=2, isot=0)
#endif
#endif
      logical redrag
      parameter(redrag=.true.)
!skip nst--->>>
      integer  nstf_name(5)
      data nstf_name /0,0,0,0,0/  !!if turn on 2,1(0),0,0,0
!skip nst---<<<
      real(kind=kind_phys), dimension(igrd1s):: sfcemis 
      real(kind=kind_phys), dimension(igrd1s):: tprcp,srflag
      logical mom4ice
      parameter(mom4ice=.false.)
!........................................
!rad:
      real(kind=kind_phys), dimension(igrd1s)::                        &
     &    gabsbdlw,adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw

!lsm:
      logical, dimension(igrd1s) :: flag_iter, flag_guess
      real(kind=kind_phys), dimension(igrd1s) :: tsurf,evbs,evcw,      &
     &         trans,sbsno,snowc,snohf,zlvl,smcwlt2,smcref2,wet1,      &
     &         fm10,fh2,tseal
     !---input/output:
      real(kind=kind_phys), dimension(igrd1s) :: evbsa,evcwa,transa,   &
     &         sbsnoa,snowca,snohfa

      real(kind=kind_phys), dimension(igrd1s) :: rzice,rcice,rtice
      integer, dimension(igrd1s) :: slopetyp, islmsk
      real(kind=kind_phys), dimension(igrd1s) :: chh, cmm

      real(kind=kind_phys), dimension(igrd1s) :: sncovr
      !---output:
!yj:     real(kind=kind_phys), dimension(igrd1s) ::                       &
!yj:    &         epi,dlwsfci,ulwsfci,uswsfci,dswsfci,gfluxi,t1,q1,u1,v1

!nst:
!skip nst:      real(kind=kind_phys), dimension(igrd1s) ::  &
!skip nst:!      real(kind=kind_phys), dimension(igrd1s), optional,  intent(inout) ::  &
!skip nst:     &      xt, xs, xu, xv, xz, zm, xtts, xzts, d_conv, ifd, dt_cool,   &
!skip nst:     &      qrain, tref, z_c, c_0, c_d, w_0, w_d
!yj-implement noah land model----------------------------------------<<<
#ifdef RAS
      logical flipv
#endif
      parameter (qmin=1.0e-10, rhzbot=0.85, rhztop=0.85)   ! zhao micro
      parameter (rhpbl=0.95,   rhbbot=0.95, rhbtop=0.90)   ! brad micro
!     parameter (              rhbbot=0.95, rhbtop=0.85)   ! brad micro
!     parameter (qmin=1.0e-10, rhcbot=0.85, rhctop=0.85)
!     parameter (qmin1=1.0e-30)
!     parameter (dxmax=log(1.0/7200.0), dxmin=log(1.0/192.0)
!     parameter (dxmax=-8.8818363, dxmin=-5.2574954
!    &,          dxinv=1.0/(dxmax-dxmin))
      parameter (dxmax=-9.5468126, dxmin=-5.2574954                      &  
     &,          dxinv=1.0/(dxmax-dxmin))
      parameter (cb2mb=10.0)
      real(kind=kind_phys) work4, dpshc, p850, p925, utmp, vtmp
      real(kind=kind_phys) frain, qi, qw, qr, wc, f_rain, f_ice, tem
!c
      real(kind=kind_phys) rclx(igrd1s),rcs(igrd1s)
!cc
      integer me, lonlens, igrd1s, ipr
!ccj
#ifdef RAS
      integer ncrnd
      real(kind=kind_phys) rannum(ncrnd), rarea(igrd1s)
#endif
      real(kind=kind_phys) clw(igrd1s,levr,2)
      real(kind=kind_phys) xmu(igrd1s)
!cc
      integer kinver(igrd1s)
      logical ldiag
      logical lprnt, invrsn(igrd1s)
#if defined(SCASAS) || defined(SAMFSAS)     
      dimension      cnvw(igrd1s,levr),cnvc(igrd1s,levr)
      dimension      sigmagfm(igrd1s)
#endif 
#ifdef wind_fitch
      integer nw
#endif     
      real(kind=kind_phys) zvfun(igrd1s),parea(igrd1s)
      logical lseaspray,lheatstrg
#ifdef CHS
! for canopy heat storage (CHS)
      real(kind=kind_phys) hffac(igrd1s)
      real(kind=kind_phys),parameter :: h0facu=0.25_kind_phys,h0facs=1.0_kind_phys
#endif
!GFDLMP>>>
      dimension cldtot(igrd1s,levr)
!      real(kind=kind_phys), parameter :: rainmin=1.0e-13 !(m)
      real(kind=kind_phys), parameter :: qmin_adj=1.0e-15 !(kg/kg)
      real(kind=kind_phys), dimension(:), allocatable ::                &
     &         frland,area,rain0,snow0,ice0,graupel0,water0
      real(kind=kind_phys), dimension(:,:), allocatable ::              &
     &         rho,rew,rei,rer,res,reg
      real(kind=kind_phys), dimension(:,:), allocatable ::              &
     &         rqv,rqc,rqr,rqi,rqs,rqg,rqa,rqn,rt0,rw0,ru0,rv0,rdp,rdz
      real(kind=kind_phys), dimension(:,:), allocatable ::              &
     &         qv_dt,qc_dt,qr_dt,qi_dt,qs_dt,qg_dt,qa_dt,pt_dt,udt,vdt
      logical  hydrostatic,phys_hydrostatic,sedi_w,effr_in
      logical  consv_te,last_step,do_inline_mp
      real(kind=kind_phys), dimension(:), allocatable ::                &
     &         gsize,hs,cond0,dep0,evap0,sub0
      real(kind=kind_phys), dimension(:,:), allocatable ::              &
     &         qnl,qni,q_con,cappa,te,prefluxr,prefluxi,prefluxs,       &
     &         prefluxg,prefluxw
!Goddard : 
!      real(kind=kind_phys) :: rhowater,rhosnow,dx
      real(kind=kind_phys),dimension(:,:),allocatable ::                &
              ht,hail2d,rainnc2d,rainncv2d,snownc2d,snowncv2d,          &
              graupelnc2d,graupelncv2d,hailnc2d,sr2d                    
      real(kind=kind_phys),dimension(:,:,:),allocatable ::              &
              t3d,qq3d,qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,rho3d,pii3d,p3d,   &
              th3d,dz3d,phii3d,rew3d,rer3d,rei3d,res3d,                    &
              reg3d,reh3d,refl_10cm
      logical diagflag 
      integer do_radar_ref
#ifdef EXT_DIAG
      real(kind=kind_phys),dimension(:,:,:),allocatable ::              &
              physc, physe, physd, physs, physm, physf,                 &
              acphysc, acphyse, acphysd, acphyss, acphysm, acphysf,     &
              preci3d, precs3d, precg3d, prech3d, precr3d
#endif
!
!     Scaling factor for GFS version of GWD parameterization scheme
!     "C_d" is scaled by cdmbgwd(1)
!     "C_m" is scaled by cdmbgwd(2)
!
!     1.0 means no scaling is applied
      real(kind=kind_phys), parameter :: &
     &  cdmbgwd(2)  = (/3.0,1.0/)!1.0_kind_phys
!ccj
#ifdef spray
      lseaspray = .true.
#else
      lseaspray = .false.
#endif     
#ifdef CHS
      lheatstrg = .true.
#else
      lheatstrg = .false.
#endif     
      lprnt = .false.
      ipr = 1
      rain1=0.0
      rain2=0.0
      rainc=0.0
!
      if (nmmiph.eq.6 .or. nmmiph.eq.15 .or. &    !WSM6 or GCE_3ice
          nmmiph.eq.11 .or. nmmiph.eq.12 .or. nmmiph.eq.13) then !GFDLMP
        nvdiff = ntotal-3  ! don't diffuse rain,snow,graupel
      elseif (nmmiph.eq.5) then !WSM5
        nvdiff = ntotal-2  ! don't diffuse rain,snow
      elseif (nmmiph.eq.4) then !WSM3
        nvdiff = ntotal-1  ! diffuse only vapor,o3,water+ice
      else
        nvdiff = ntotal         !GFDLMP
      endif

!      
!yj-implement noah land model---------------------------------------->>>
      sfcemis=1.0
      tprcp=0.0 
!yj-implement noah land model----------------------------------------<<<
#ifdef LTN
!lightning
      flash=0.0     !flash density (unit in flashes km^-2 day^-1)
      snow_flx=0.0
      ptu=0.0
      pqu=0.0
#endif
!
!cccccccccccccccccccccccccccccc
!
      do i=1,lonlens
          psexp(i)= exp(pgr(i))   ! in cbar
      enddo
! use si and sl for reference sigma value (for physics parameterization)
      do k=1,levr+1
      do i=1,lonlens
        sio(i,k)=(sig(k,2)+sig(k,1)*psexp(i)*10.)/(psexp(i)*10.)
!       sio(i,k)=si(k)
      enddo
      enddo
      do k=1,levr
      do i=1,lonlens
        slo(i,k)=(sio(i,k)+sio(i,k+1))/2.
!       slo(i,k)=((sio(i,k)**(1.+rocp)-sio(i,k+1)**(1.+rocp))/      &
!    &          ((1.+rocp)*(sio(i,k)-sio(i,k+1))))**(1./rocp)
!       slo(i,k)=sl(k)
      enddo
      enddo
      do k=1,levr
      do i=1,lonlens
        delo(i,k)=sio(i,k)-sio(i,k+1)
      enddo
      enddo
!c
      do k=1,levr+1
      do i=1,lonlens
         sik(i,k)=sio(i,k)**rocp
      enddo
      enddo
      do k=1,levr
      do i=1,lonlens
         slk(i,k)=slo(i,k)**rocp
      enddo
      enddo
!c
      dpshc = 30.0                ! in cbar
!
!     prsl = kPa(cbar) ,prsl_pa = Pa
      do i=1,lonlens
          psexp_pa(i)=psexp(i)*1000.0    ! in pa
          work4          = fpkap(1000.0*psexp(i))
          prsik(i,1)     = work4
          do k=1,levr
            prsl(i,k)    = psexp(i)     * slo(i,k)
            prsl_pa(i,k) = prsl(i,k)*1.e3
            prsi(i,k)    = psexp(i)     * sio(i,k)
            prsi_pa(i,k) = prsi(i,k)*1.e3
            prslk(i,k)   = work4        * slk(i,k)
            prsik(i,k+1) = work4        * sik(i,k+1)
          enddo
          prsi(i,levr+1) = psexp(i) * sio(i,levr+1)
          prsi_pa(i,levr+1) = prsi(i,levr+1)*1.e3
      enddo
      snwdph_mm(:)=snwdph(:,lat)*1000.
      do k=1,levr
        do i=1,lonlens
          phii(i,k) = 0.0
          phil(i,k) = 0.0
        enddo
      enddo
      do i=1,lonlens
        phii(i,levr+1) = 0.0
      enddo
!
!  get dry tep.
!
      do 270 k=1,levr
      do 270 j=1,lonlens
#ifdef T
      i=j+(lat-1)*igrd1s
      wt(i,k,11)=tgrs(j,k)
#endif
      tgrs(j,k)=tgrs(j,k)/(1.+fv*max(qgrs(j,k),qmin))
#ifdef NONHYD
      tgrh(j,k)=tgrh(j,k)/(1.+fv*max(qgrs(j,k),qmin))
#endif
#ifdef T
      wt(i,k,11)=tgrs(j,k)-wt(i,k,11)
#endif
  270 continue
!
#ifndef SIGP
#ifndef NONHYD
      call get_prs(lonlens,igrd1s,levr,rocp,cp,fv,tgrs,qgrs             &
     &,            prsi,prsik,prsl,prslk,phii,phil,delp)
#else
      call get_prs(lonlens,igrd1s,levr,rocp,cp,fv,tgrh,qgrs             &  
     &,            prsi,prsik,prsl,prslk,phii,phil,delp)
#endif
#else
#ifndef NONHYD
      call get_prs_sp(lonlens,igrd1s,levr,rocp,cp,fv,tgrs,qgrs          &
     &,            prsi,prsik,prsl,prslk,phii,phil,delp)
#else
      call get_prs_sp(lonlens,igrd1s,levr,rocp,cp,fv,tgrh,qgrs          &  
     &,            prsi,prsik,prsl,prslk,phii,phil,delp)
#endif
#endif
      do i=1,lonlens
          do k=1,levr
            delp_pa(i,k) = delp(i,k)*1.e3
          enddo
      enddo
!
!ccj
      p850 = 850. * 0.1
      p925 = 925. * 0.1
!
      rbs2 = 1.e0
      jcapr= 126.*106000./ rdelx
!
!  frain is the factor for centered difference scheme correction of rain
!  amount.
!
      frain = dtf/dtp
      ldiag=.false.
!
!  transfer vegetation fraction from global to local location
!  as well as vegetation type and soil type
!
!yj-implement noah land model---------------------------------------->>>
      do j = 1, lonlens
        sigmaf(j) = max(vfrac(j,lat),0.01)
!       sigmaf(j) = max(vfrac(j,lat),.3)
        if (lsm .eq.0) sigmaf(j) = 0.5 + vfrac(j,lat) * 0.5

        islmsk(j)   = nint(slmsk(j,lat))

        if(islmsk(j).eq.2) then
          if (isot.eq.1) then
            soiltyp(j) = 16
          else
            soiltyp(j) = 9
          endif
          if (ivegsrc.eq.1) then
            vegtype(j) = 15
          elseif (ivegsrc.eq.2) then
            vegtype(j) = 13
          endif
          slopetyp(j) = 9  !! clu: qa(clopetyp)
        else
          soiltyp(j) = int(stype(j,lat)+.5)
          vegtype(j) = int(vtype(j,lat)+.5)
          slopetyp(j) = int(slope(j,lat)+0.5 ) !! clu: slope -> slopetyp
        endif

        if(islmsk(j).eq.1) then
          if (isot.eq.1) then
            if(soiltyp(j).eq.0)  soiltyp(j) = 6   !loam
          else
            if(soiltyp(j).eq.0)  soiltyp(j) = 2   !silty clay loam
          endif
          if (ivegsrc.eq.1) then
            if(vegtype(j).eq.0)  vegtype(j) = 10  !grasslands
          elseif (ivegsrc.eq.2) then
            if(vegtype(j).eq.0)  vegtype(j) = 7   !groundcover only (perennial)
          endif
        endif
! ---...xw: transfer ice thickness & concentration from global to local
! variables
!yj: skip sea ice--->>> (simulation domain at subtropical)
!     zice(j) = hice(j)
!     cice(j) = fice(j)
!     tice(j) = tisfc(j)
     rzice(j) = 0.
     rcice(j) = 0.
     if (islmsk(j) == 2) then
       rtice(j) = 271.2
     else
       rtice(j)=tsea(j,lat)
     endif
!yj: skip sea ice---<<<

        fm1(j) = ffmm(j,lat)
        fh1(j) = ffhh(j,lat)
        ustar1(j) = uustar(j,lat)

      enddo
!
!  transfer soil moisture and temperature from global to local variables
!
      do k = 1, lsoil
        do j = 1, lonlens
          smsoil(j,k) = smc(j,lat,k)
          stsoil(j,k) = stc(j,lat,k)
!yj: replace slsoil(unfrozen soil moisture content(fraction)) 
!yj: by smsoil(total soil moisture content (fraction))       
          slsoil(j,k) = smc(j,lat,k)  !! clu: slc -> slsoil
        enddo
      enddo
!
      do 240 k=1,levr
      do 240 j=1,lonlens
      dudt(j,k)=0.e0
      dvdt(j,k)=0.e0
      dtdt(j,k)=0.e0
  240 continue
#ifdef NONHYD
      do k=1,levr+1
        do j=1,lonlens
          dwdt(j,k)=0.e0
        enddo
      enddo
#endif
      do 250 k=1,levh
      do 250 j=1,lonlens
      dqdt(j,k)=0.e0
  250 continue
      rcl=1.0
      do i=1,lonlens
        rclx(i)   = 1.0
        rcs(i)   = sqrt(rclx(i))
        work1(i) = prsik(i,1) / prslk(i,1)
        parea(i) = rdelx * rdely / xm2(i,lat)
#ifdef RAS
!       garea(i) = rerth * rerth * (pi+pi)*pi/(rcs(i)*lonf*latg)
        rarea(i) = rdelx * rdely / xm2(i,lat)
!        if(me.eq.0.and.i.eq.lonlens/2) then
!          print*,' lat= ',lat,' rarea= ',rarea(i),
!     &           ' delx= ',rdelx/sqrt(xm2(i,lat))
!        endif
#endif
      enddo
!yj-implement noah land model----------------------------------------<<<
!.................................................................
!
!  get surface pressure.
!
      do 280 j=1,lonlens
      psurf(j,lat)=psexp(j)
      if(lssav) psmean(j,lat)=psmean(j,lat)+psexp(j)*dtf
  280 continue

#ifdef doDCYC2
!
!  initialize dtdt with heating rate from dcyc2 and get radsl for progtm
!
      call rdcyc2t3(lonlens,igrd1s,levr,solhr,slag,                       &   
     &           sinlar(1,lat),coslar(1,lat),                            &
     &           sdec,cdec,rlon(1,lat),coszer(1,lat),                    &
     &           sfcdlw(1,lat),sfcnsw(1,lat),tgrs(1,1),                  &
     &           tsea(1,lat),tsflw(1,lat),swh(1,1,lat),hlw(1,1,lat),     &
     &           dlwsf1,ulwsf1,radsl,dtdt,xmu)
!
!yj: tmp variables for future implement of new radiation to noah land model>
! sswdnf: net downward sw flux (positive downward)
! sfcnsw: net sw flux (positive upward*)
! dlwsf1: downward lonwave flux (positive downward)
! radsl : ((net sw flux)  - dlwsf1) (positive upward*)
! ulwsf1: upward longwave flux (positive upward)
! note: for Noah LSM, 
!         adjsfcdsw: positive downward
!         adjsfcnsw: positive downward*
!         adjsfcdlw: positive downward
!         adjsfculw: positive upward
      do j = 1, lonlens
      adjsfcdsw(j)=sswdnf(j,lat)*xmu(j)   !sfc: downward short wave(sw) = upward sw - net sw
      adjsfcnsw(j)=-radsl(j)-dlwsf1(j)
      adjsfcdlw(j)=dlwsf1(j)
      adjsfculw(j)=sfcemis(j)*ulwsf1(j)+(1.0-sfcemis(j))*adjsfcdlw(j)
      enddo

!  --- ...  define the downward lw flux absorbed by ground

      do j = 1, lonlens
        gabsbdlw(j) = sfcemis(j) * adjsfcdlw(j)
      enddo
!yj: tmp variables for future implement of new radiation to noah land model<
!
!       print *,'in rbphys, aft dcyc2t3, me=',me
      if(lssav) then
      do 285 j=1,lonlens
        !yj:keep upward positive
      dlwsfc(j,lat)=dlwsfc(j,lat)+adjsfcdlw(j)*dtf
      ulwsfc(j,lat)=ulwsfc(j,lat)+adjsfculw(j)*dtf
!yj      dlwsfc(j,lat)=dlwsfc(j,lat)+dlwsf1(j)*dtf
!yj      ulwsfc(j,lat)=ulwsfc(j,lat)+ulwsf1(j)*dtf
285   continue
      endif
!
      do 286 k=1,levr
      do 286 j=1,lonlens
      hsw(j,k)=dtdt(j,k)-hlw(j,k,lat)
#ifdef T
      i=j+(lat-1)*igrd1s
      wt(i,k, 9)=hsw(j,k)
      wt(i,k,10)=hlw(j,k,lat)
#endif
286   continue
#endif  
! end of doDCYC2
!
!
!cc
!     phy_f2dv(:,num_p2d) = 0.0  ! disable downdraft effect on evap
!       if (me .eq. 0) print *,'--- in rbphys, before progtm'
!       if (me .eq. 0) print *,'-- tsea  ', tsea(:,lat)
!
#ifdef doLSM
#ifdef doNOAHLSM
!yj-implement noah land model---------------------------------------->>>
!
! --- ... lu: initialize flag_guess, flag_iter, tsurf

      do j = 1, lonlens
        tsurf(j)      = tsea(j,lat)
        flag_guess(j) = .false.
        flag_iter(j)  = .true.
        drain(j)      = 0.0
        ep1d(j)       = 0.0
        runof(j)      = 0.0
        hflx(j)       = 0.0
        evap(j)       = 0.0

        evbs(j)       = 0.0
        evcw(j)       = 0.0
        trans(j)      = 0.0
        sbsno(j)      = 0.0
        snowc(j)      = 0.0
        snohf(j)      = 0.0
        zlvl(j)       = phil(j,1) / grav
        smcwlt2(j)    = 0.0
        smcref2(j)    = 0.0
      enddo
!      
!yj: refer CWBGFS's setting
!  estimate t850 for noah land model sfc_drv, uv850-uv925=uvdiff( for gusts)
!
      do j=1,lonlens
        t850(j) = tgrs(j,1)
        uv850(j)=sqrt(ugrs(j,1)*ugrs(j,1)+vgrs(j,1)*vgrs(j,1))
        uv925(j)=sqrt(ugrs(j,1)*ugrs(j,1)+vgrs(j,1)*vgrs(j,1))
      enddo
      do k = 1, levr - 1
        do j=1,lonlens
          if(prsl(j,k) .gt. p850 .and. prsl(j,k+1) .le. p850) then
            t850(j) = tgrs(j,k) - (prsl(j,k)-p850)                        &   
     &              / (prsl(j,k)-prsl(j,k+1)) * (tgrs(j,k)-tgrs(j,k+1))
            utmp = ugrs(j,k) - (prsl(j,k)-p850)                        &   
     &              / (prsl(j,k)-prsl(j,k+1)) * (ugrs(j,k)-vgrs(j,k+1))
            vtmp = vgrs(j,k) - (prsl(j,k)-p850)                        &   
     &              / (prsl(j,k)-prsl(j,k+1)) * (vgrs(j,k)-vgrs(j,k+1))
            uv850(j)=sqrt(utmp*utmp+vtmp*vtmp)
          endif
          if(prsl(j,k) .gt. p925 .and. prsl(j,k+1) .le. p925) then
            utmp = ugrs(j,k) - (prsl(j,k)-p925)                        &   
     &              / (prsl(j,k)-prsl(j,k+1)) * (ugrs(j,k)-vgrs(j,k+1))
            vtmp = vgrs(j,k) - (prsl(j,k)-p925)                        &   
     &              / (prsl(j,k)-prsl(j,k+1)) * (vgrs(j,k)-vgrs(j,k+1))
            uv925(j)=sqrt(utmp*utmp+vtmp*vtmp)
          endif
        enddo
      enddo
!
        do j=1,lonlens
        srflag(j) = 0.
        if(t850(j).le.273.16)then
          srflag(j)=1.
        endif
         uvdiff(j)=uv850(j)-uv925(j)
        enddo


!  --- ...  lu: iter-loop over (sfc_diff,sfc_drv,sfc_ocean,sfc_sice)

      do iter = 1, 2

!  --- ...  surface exchange coefficients
!
!      sheleg unit:[mm], snwdph unit: [m]
      if (lprnt) write(0,*)' tsea=',tsea(ipr,lat),' tsurf=',tsurf(ipr),iter
        call rsfc_diff(lonlens,igrd1s,psexp_pa,ugrs,vgrs,tgrs,qgrs,zlvl, &
     &         snwdph_mm,tsea(1,lat),zorl(1,lat),cd,cdq,rb,             &
     &                prsl_pa,work1,islmsk,                             &
     &                stress,fm1,fh1,zvfun,parea,ssu(1,lat),ssv(1,lat), &
     &                ustar1,wind,phy_f2dv(1,num_p2d),fm10,fh2,         &
     &                sigmaf,vegtype,shdmax(1,lat),ivegsrc,             &
     &                tsurf, flag_iter, redrag)

       if (lprnt) write(0,*)' cdq=',cdq(ipr),' iter=',iter             &
    &,   ' wind=',wind(ipr),'phy_f2dv=',phy_f2dv(ipr,num_p2d),' ugrs='   &
    &,   ugrs(ipr,1),' vgrs=',vgrs(ipr,1)

!  --- ...  lu: update flag_guess

        do j = 1, lonlens
          if (iter == 1 .and. wind(j) < 2.0) then
            flag_guess(j) = .true.
          endif
        enddo

!skip nst:        if ( nstf_name(1) > 0 ) then
!skip nst:
!skip nst:          do j = 1, lonlens
!skip nst:            if ( islmsk(j) == 0 ) then
!skip nst:              tem      = (oro(j)-oro_uf(j)) * rlapse
!skip nst:              tseal(j) = tsea(j)  + tem
!skip nst:              tsurf(j) = tsurf(j) + tem
!skip nst:            endif
!skip nst:          enddo
!skip nst:
!skip nst:!       if (lprnt) write(0,*)' tseaz1=',tsea(ipr),' tref=',tref(ipr)
!skip nst:!    &,   ' dt_cool=',dt_cool(ipr),' dt_warm=',2.0*(xt(ipr)/xz(ipr)
!skip nst:!    &,   ' kdt=',kdt
!skip nst:!    &,   ' tgrs=',tgrs(ipr,1),' prsl=',prsl(ipr,1)
!skip nst:!    &,   ' work3=',work3(ipr),' kdt=',kdt
!skip nst:
!skip nst:!yj: lack xlon, sinlat, some radiation related variables
!skip nst:          call sfc_nst                                                  &
!skip nst:     &       ( lonens,lsoil,psexp_pa,ugrs,vgrs,tgrs,qgrs,tref,cd,cdq,   &
!skip nst:     &         prsl_pa,work1,islmsk,xlon,sinlat,stress,                 &
!skip nst:     &         sfcemis,gabsbdlw,adjsfcnsw,tprcp,dtf,kdt,solhr,xcosz,    &
!skip nst:     &         phy_f2dv(1,num_p2d),flag_iter,flag_guess,nstf_name,      &
!skip nst:     &         lprnt,ipr,                                               &
!skip nst:!  --- Input/output
!skip nst:     &         tseal,tsurf,xt,xs,xu,xv,xz,zm,xtts,xzts,dt_cool,         &
!skip nst:     &         z_c,c_0,c_d,w_0,w_d,d_conv,ifd,qrain,                    &
!skip nst:!  ---  outputs:
!skip nst:     &         qss, gflx, cmm, chh, evap, hflx, ep1d)
!skip nst:
!skip nst:!         if (lprnt) print *,' tseaz2=',tseal(ipr),' tref=',tref(ipr),
!skip nst:!    &     ' dt_cool=',dt_cool(ipr),' dt_warm=',2.0*xt(ipr)/xz(ipr),
!skip nst:!    &     ' kdt=',kdt
!skip nst:
!skip nst:          do i = 1, im
!skip nst:            if ( islmsk(i) == 0 ) then
!skip nst:              tsurf(i) = tsurf(i) - (oro(i)-oro_uf(i)) * rlapse
!skip nst:            endif
!skip nst:          enddo
!skip nst:
!skip nst:!  --- ...  run nsst model  ... ---
!skip nst:
!skip nst:          if ( nstf_name(1) > 1 ) then
!skip nst:            zsea1 = 0.001*real(nstf_name(4))
!skip nst:            zsea2 = 0.001*real(nstf_name(5))
!skip nst:            call get_dtzm_2d(xt,xz,dt_cool,z_c,slmsk,
!skip nst:     &                       zsea1,zsea2,im,1,dtzm)
!skip nst:            do i = 1, im
!skip nst:              if ( islmsk(i) == 0 ) then
!skip nst:              tsea(i) = max(271.2,tref(i) + dtzm(i))
!skip nst:     &                      -(oro(i)-oro_uf(i))*rlapse
!skip nst:              endif
!skip nst:            enddo
!skip nst:          endif
!skip nst:
!skip nst:!         if (lprnt) print *,' tseaz2=',tsea(ipr),' tref=',tref(ipr),   &
!skip nst:!    &    ' dt_cool=',dt_cool(ipr),' dt_warm=',dt_warm(ipr),' kdt=',kdt
!skip nst:
!skip nst:        else

!  --- ...  surface energy balance over ocean

          call rsfc_ocean                                                &
!  ---  inputs:
     &     ( lonlens,igrd1s,psexp_pa,ugrs,vgrs,tgrs,qgrs,tsea(1,lat),   &
     &       cd,cdq,prsl_pa,work1,islmsk,phy_f2dv(1,num_p2d),flag_iter, &
     &       fm10,fm1,lseaspray,ssu(1,lat), ssv(1,lat), &
!  ---  outputs:
     &       qsurf,cmm,chh,gflx,evap,hflx,ep1d                            &
     &     )

!skip nst:        endif       ! if ( nstf_name(1) > 0 ) then

!       if (lprnt) write(0,*)' sfalb=',sfalb(ipr),' ipr=',ipr           &
!    &,   ' weasd=',weasd(ipr),' snwdph=',snwdph(ipr)                   &
!    &,   ' tprcp=',tprcp(ipr),' kdt=',kdt,' iter=',iter
!    &,' tseabefland=',tsea(ipr)

!  --- ...  surface energy balance over land
!
        if (lsm == 1) then                          ! noah lsm call

!     if (lprnt) write(0,*)' tsead=',tsea(ipr),' tsurf=',tsurf(ipr),iter
!    &,' pgr=',pgr(ipr),' sfcemis=',sfcemis(ipr)

          call rsfc_drv                                                  &

!  ---  inputs:
     &     ( lonlens,igrd1s,lsoil,psexp_pa,ugrs,vgrs,tgrs,qgrs,soiltyp, &
     &       vegtype,sigmaf,sfcemis,gabsbdlw,adjsfcdsw,adjsfcnsw,dtf,   &
     &       tg3(1,lat),cd,cdq,prsl_pa,work1,zlvl,islmsk,               &
     &       phy_f2dv(1,num_p2d),slopetyp,shdmin(1,lat),shdmax(1,lat),  &
     &       snoalb(1,lat),sfalb(1,lat),flag_iter,flag_guess,           &
     &       .false.  ,isot,ivegsrc,                                    &
!    &       lheatstrg,isot,ivegsrc,                                    &
!  ---  in/outs:
!yj: use sheleg to replace weasd, use smsoil to replace slsoil
     &       sheleg(1,lat),snwdph_mm,tsea(1,lat),tprcp,srflag,smsoil,   &
     &       stsoil,slsoil,canopy(1,lat),trans,tsurf,zorl(1,lat),       &
!  ---  outputs:
     &       sncovr,qsurf,gflx,drain,evap,hflx,ep1d,runof,              &
     &       cmm,chh,evbs,evcw,sbsno,snowc,soilm(1,lat),snohf,          &
     &       smcwlt2,smcref2,wet1                                       &
     &     )

      if (lprnt) write(0,*)' tseae=',tsea(ipr,lat),' tsurf=',tsurf(ipr),iter &
     &,' phy_f2dv=',phy_f2dv(ipr,num_p2d)

        else                                       ! osu lsm call

          call rsfc_land                                                 &
!  ---  inputs:
     &     ( lonlens,igrd1s,lsoil,psexp_pa,ugrs,vgrs,tgrs,qgrs,smsoil,  &
     &       soiltyp,sigmaf,vegtype,sfcemis,adjsfcdlw,adjsfcnsw,dtf,    &
     &                   tg3(1,lat),cd,cdq,prsl_pa,work1,islmsk,        &
!    &       zorl(1,lat),tg3(1,lat),cd,cdq,prsl(1,1),work3,islmsk,      &
     &       phy_f2dv(1,num_p2d),flag_iter,flag_guess,                  &
!  ---  input/outputs:
!yj: use sheleg to replace weasd
     &       sheleg(1,lat),tsea(1,lat),tprcp,srflag,stsoil,             &
     &       canopy(1,lat),tsurf,                                       &
!  ---  outputs:
     &       qsurf,snowmt,gflx,zsoil,rhscnpy,rhsmc,                     &
     &       ai,bi,cci,drain,evap,hflx,ep1d,cmm,chh,                    &
     &       evbs,evcw,trans,sbsno,snowc,soilm(1,lat),                  &
     &       snohf,smcwlt2,smcref2                                      &
     &     )

        endif

        if (lprnt) write(0,*)' tseabeficemodel =',tsea(ipr,lat),' me=',me   &
     &,   ' kdt=',kdt

!  --- ...  surface energy balance over seaice

!yj        if (lssav_cpl) then
!yj          do i = 1, im
!yj            if (flag_cice(i)) then
!yj               islmsk(i) = islmsk_cice(i)
!yj            endif
!yj          enddo
!yj        endif

        call rsfc_sice                                                   &
!  ---  inputs:
     &     ( lonlens,igrd1s,lsoil,psexp_pa,ugrs,vgrs,tgrs,qgrs,dtf,     &
     &       sfcemis,gabsbdlw,adjsfcnsw,adjsfcdsw,srflag,               &
     &       cd,cdq,prsl_pa,work1,islmsk,phy_f2dv(1,num_p2d),           &
     &       flag_iter,mom4ice,lsm, lprnt,ipr,                          &
!    &       flag_iter,mom4ice,lsm,                                     &
!  ---  input/outputs:
!yj: use sheleg to replace weasd
     &       rzice,rcice,rtice,sheleg(1,lat),tsea(1,lat),tprcp,stsoil,  &
     &       ep1d,                                                      &
!  ---  outputs:
     &       snwdph(1,lat),qsurf,snowmt,gflx,cmm,chh,evap,hflx          &
     &     )

!yj        if (lssav_cpl) then
!yj          do i = 1, im
!yj            if (flag_cice(i)) then
!yj               islmsk(i) = nint(slmsk(i))
!yj            endif
!yj          enddo
!yj
!yj          call sfc_cice                                                 &
!yj!  ---     inputs:
!yj     &     ( im,ugrs,vgrs,tgrs,qgrs,cd,cdq,prsl(1,1),work3,             &
!yj     &       islmsk_cice,phy_f2d(1,num_p2d),flag_iter,                  &
!yj     &       dqsfc_cice,dtsfc_cice,                                     &
!yj!  ---     outputs:
!yj     &       qss,cmm,chh,evap,hflx                                      &
!yj     &     )
!yj        endif

!  --- ...  lu: update flag_iter and flag_guess

        do i = 1, lonlens
          flag_iter(i)  = .false.
          flag_guess(i) = .false.

          if(islmsk(i) == 1 .and. iter == 1) then
            if (wind(i) < 2.0) flag_iter(i) = .true.
          elseif (islmsk(i) == 0 .and. iter == 1                        &
     &                           .and. nstf_name(1) > 0) then
            if (wind(i) < 2.0) flag_iter(i) = .true.
          endif
        enddo

      enddo   ! end iter_loop

!yj:  do i = 1, lonlens
!yj:    epi(i)     = ep1d(i)
!yj:    dlwsfci(i) = adjsfcdlw(i)
!yj:    ulwsfci(i) = adjsfculw(i)
!yj:    uswsfci(i) = adjsfcdsw(i) - adjsfcnsw(i)
!yj:    dswsfci(i) = adjsfcdsw(i)
!yj:    gfluxi(i)  = gflx(i)
!yj:    t1(i)      = tgrs(i,1)
!yj:    q1(i)      = qgrs(i,1)
!yj:    u1(i)      = ugrs(i,1)
!yj:    v1(i)      = vgrs(i,1)
!yj:  enddo

      if (lsm == 0) then                          ! osu lsm call
        do i = 1, lonlens
         sncovr(i) = 0.0
!yj: use sheleg to replace weasd
         if (sheleg(i,lat) > 0.0) sncovr(i) = 1.0
        enddo
      endif

!  --- ...  update near surface fields

      call rsfc_diag(lonlens,igrd1s,psexp_pa,ugrs,vgrs,tgrs,qgrs,qss,       &
     &              tsea(1,lat),qsurf,f10m(1,lat),u10m(1,lat),v10m(1,lat), &
     &              t2m(1,lat),q2m(1,lat),work1,                           &
     &              evap,fm1,fh1,fm10,fh2)

      do i = 1, lonlens
        phy_f2dv(i,num_p2d) = 0.0
      enddo

!yj      if (lssav_cpl) then
!yj        do i = 1, im
!yj          dlwsfci_cpl(i)   = adjsfcdlw(i)
!yj          dswsfci_cpl(i)   = adjsfcdsw(i)
!yj          dlwsfc_cpl(i)    = dlwsfc_cpl(i) + adjsfcdlw(i)*dtf
!yj          dswsfc_cpl(i)    = dswsfc_cpl(i) + adjsfcdsw(i)*dtf
!yj          dnirbmi_cpl(i)   = adjnirbmd(i)
!yj          dnirdfi_cpl(i)   = adjnirdfd(i)
!yj          dvisbmi_cpl(i)   = adjvisbmd(i)
!yj          dvisdfi_cpl(i)   = adjvisdfd(i)
!yj          dnirbm_cpl(i)    = dnirbm_cpl(i) + adjnirbmd(i)*dtf
!yj          dnirdf_cpl(i)    = dnirdf_cpl(i) + adjnirdfd(i)*dtf
!yj          dvisbm_cpl(i)    = dvisbm_cpl(i) + adjvisbmd(i)*dtf
!yj          dvisdf_cpl(i)    = dvisdf_cpl(i) + adjvisdfd(i)*dtf
!yj          nlwsfci_cpl(i)   = adjsfcdlw(i)  - adjsfculw(i)
!yj          nlwsfc_cpl(i)    = nlwsfc_cpl(i) + nlwsfci_cpl(i)*dtf
!yj          t2mi_cpl(i)      = t2m(i)
!yj          q2mi_cpl(i)      = q2m(i)
!yj          u10mi_cpl(i)     = u10m(i)
!yj          v10mi_cpl(i)     = v10m(i)
!yj          tseai_cpl(i)     = tsea(i)
!yj          psurfi_cpl(i)    = pgr(i)
!yj        enddo
!yj
!yj!  ---  estimate mean albedo for ocean point without ice cover and apply
!yj!       them to net SW heat fluxes
!yj
!yj        do i = 1, im
!yj          if (islmsk(i) /= 1) then  ! not a land point
!yj!  ---  compute open water albedo
!yj            xcosz_loc = max( 0.0, min( 1.0, xcosz(i) ))
!yj            ocalnirdf_cpl(i) = 0.06
!yj            ocalnirbm_cpl(i) = max(albdf, 0.026/(xcosz_loc**1.7+0.065)  &
!yj     &                       + 0.15 * (xcosz_loc-0.1) * (xcosz_loc-0.5) &
!yj     &                       * (xcosz_loc-1.0))
!yj            ocalvisdf_cpl(i) = 0.06
!yj            ocalvisbm_cpl(i) = ocalnirbm_cpl(i)
!yj
!yj            nnirbmi_cpl(i) = adjnirbmd(i)-adjnirbmd(i)*ocalnirbm_cpl(i)
!yj            nnirdfi_cpl(i) = adjnirdfd(i)-adjnirdfd(i)*ocalnirdf_cpl(i)
!yj            nvisbmi_cpl(i) = adjvisbmd(i)-adjvisbmd(i)*ocalvisbm_cpl(i)
!yj            nvisdfi_cpl(i) = adjvisdfd(i)-adjvisdfd(i)*ocalvisdf_cpl(i)
!yj          else
!yj            nnirbmi_cpl(i) = adjnirbmd(i) - adjnirbmu(i)
!yj            nnirdfi_cpl(i) = adjnirdfd(i) - adjnirdfu(i)
!yj            nvisbmi_cpl(i) = adjvisbmd(i) - adjvisbmu(i)
!yj            nvisdfi_cpl(i) = adjvisdfd(i) - adjvisdfu(i)
!yj          endif
!yj          nswsfci_cpl(i) = nnirbmi_cpl(i) + nnirdfi_cpl(i)              &
!yj     &                   + nvisbmi_cpl(i) + nvisdfi_cpl(i)
!yj          nswsfc_cpl(i)  = nswsfc_cpl(i)  + nswsfci_cpl(i)*dtf
!yj          nnirbm_cpl(i)  = nnirbm_cpl(i)  + nnirbmi_cpl(i)*dtf
!yj          nnirdf_cpl(i)  = nnirdf_cpl(i)  + nnirdfi_cpl(i)*dtf
!yj          nvisbm_cpl(i)  = nvisbm_cpl(i)  + nvisbmi_cpl(i)*dtf
!yj          nvisdf_cpl(i)  = nvisdf_cpl(i)  + nvisdfi_cpl(i)*dtf
!yj        enddo
!yj      endif

      if (lssav) then
        do i = 1, lonlens  !! 290 do loop
          gflux(i,lat)   = gflux(i,lat)  + gflx(i)  * dtf
          evbsa(i)   = evbsa(i)  + evbs(i)  * dtf
          evcwa(i)   = evcwa(i)  + evcw(i)  * dtf
          transa(i)  = transa(i) + trans(i) * dtf
          sbsnoa(i)  = sbsnoa(i) + sbsno(i) * dtf
          snowca(i)  = snowca(i) + snowc(i) * dtf
          snohfa(i)  = snohfa(i) + snohf(i) * dtf
          ep(i,lat)      = ep(i,lat)     + ep1d(i)  * dtf

          tmpmax(i,lat)  = max(tmpmax(i,lat),t2m(i,lat))
          tmpmin(i,lat)  = min(tmpmin(i,lat),t2m(i,lat))

          spfhmax(i,lat) = max(spfhmax(i,lat),q2m(i,lat))
          spfhmin(i,lat) = min(spfhmin(i,lat),q2m(i,lat))
          uustar(i,lat) = ustar1(i)
          ffmm(i,lat) = fm1(i)
          ffhh(i,lat) = fh1(i)
        enddo
      endif

!!!!!!!!!!!!!!!!!Commented by Moorthi on July 18, 2012 !!!!!!!!!!!!!!!!!!!
!     do i = 1, lonlens
!  --- ...  compute coefficient of evaporation in evapc
!
!       if (evapc(i) > 1.0e0) evapc(i) = 1.0e0
!  --- ...  over snow cover or ice or sea, coef of evap =1.0e0
!       if (weasd(i) > 0.0 .or. slmsk(i) /= 1.0) evapc(i) = 1.0e0
!     enddo
!1874
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
!end noah land model

!do progtm (osu?)
#ifdef doPROGTM
      call progtm(lonlens,igrd1s,lsoil,psexp,                            &   
     & ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),phil(1,1),                &
     & sheleg(1,lat),tsea(1,lat),qss,                                    &
     & smsoil,stsoil,evapc,soiltyp,sigmaf,                               &
     & vegtype,canopy(1,lat),dlwsf1,radsl,                               &
     & snowmt,dtf,zorl(1,lat),                                           &
     & tg3(1,lat),gflx,f10m(1,lat),u10m(1,lat),v10m(1,lat),              &
     & t2m(1,lat),q2m(1,lat),                                            &
     & zsoil,cd,cdq,rb,rhscnpy,rhsmc,ai,bi,cci,                          &
     & rclx,prsl(1,1),work1,slmsk(1,lat),                                &
     & drain,evap,hflx,stress,ep1d,                                      &
     & fm1,fh1,ustar1,wind,phy_f2dv(1,num_p2d))
!
!       if (me .eq. 0) print *,'--- in rbphys, after progtm'
!       if (me .eq. 0) print *,'-- tsea  ', tsea(:,lat)
!       print *,'in rbphys,progtm,me=',me,'lat=',lat
      do i=1,lonlens
        phy_f2dv(i,num_p2d) = 0.0
      enddo
!
      if(lssav) then
      do 290 j=1,lonlens
      gflux(j,lat)=gflux(j,lat)+gflx(j)*dtf
      tmpmax(j,lat) = max(tmpmax(j,lat),t2m(j,lat))
      tmpmin(j,lat) = min(tmpmin(j,lat),t2m(j,lat))
      spfhmax(j,lat) = max(spfhmax(j,lat),q2m(j,lat))
      spfhmin(j,lat) = min(spfhmin(j,lat),q2m(j,lat))
      ep(j,lat) = ep(j,lat) + ep1d(j) * dtf
      uustar(j,lat) = ustar1(j)
      ffmm(j,lat) = fm1(j)
      ffhh(j,lat) = fh1(j)
  290 continue
      endif
!
!     compute coefficient of evaporation in evapc
!
      do 320 j=1,lonlens
      if (evapc(j).gt.1.e0) evapc(j)=1.0e0
  320 continue
!
!     over snow cover or ice or sea, coef of evap =1.0.e0
!
      do 330 j=1,lonlens
      if ((sheleg(j,lat).gt.0.e0) .or.                                 &  
     &     (slmsk(j,lat).ne.1.e0))                                     &
     &      evapc(j)=1.e0
  330 continue
#endif
! end of doPROGTM
#endif
! end of doLSM
!
!ccj
      if(ncld.gt.0) then
        ncw = (1+ntoz)*levr
      else
        ncw = 0
      endif

      if (ncld.eq.6) then !GCE_4ICE
        nrw=1*levr+ncw
        niw=2*levr+ncw
        nsw=3*levr+ncw
        ngl=4*levr+ncw
        nhl=5*levr+ncw 
      elseif (ncld.eq.5 ) then !WSM6 or GFDLMP or GCE_3ICE
        nrw=1*levr+ncw
        niw=2*levr+ncw
        nsw=3*levr+ncw
        ngl=4*levr+ncw
      elseif (ncld.eq.4 ) then !WSM5
        nrw=1*levr+ncw
        niw=2*levr+ncw
        nsw=3*levr+ncw
      elseif (ncld.eq.2 ) then !WSM3
        nrw=1*levr+ncw
      elseif (ncld.eq.1 ) then !Zhao-Carr
        ncw=(1+ntoz)*levr
      else
        if (me.eq.0) print*,'Error: ncld undefined!'
        stop 99
      endif
!      deallocate ( rqv,rqc,rqr,rqi,rqs,rqg )
!
#ifdef fix_neg_q
      if ( ncld .eq. 5 ) then
         ix = igrd1s
         im = lonlens
         call fix_neg_q6(ix, im, levr, levh, ncw, nrw, niw, nsw, ngl,   &
                         delp_pa, tgrs, qgrs, qcmin)
      endif
#endif
!
!  do vertical diffusion
!
      allocate(qgrsx(igrd1s,levr,nvdiff), dqdtx(igrd1s,levr,nvdiff))
      dqdtx(:,:,:)=0.0
      if (nvdiff.eq.7) then     !GFDLMP
        do k=1,levr
        do i=1,lonlens
           qgrsx(i,k,1)=qgrs(i,k)
           qgrsx(i,k,2)=qgrs(i,k+ntoz*levr)
           qgrsx(i,k,3)=qgrs(i,k+ncw)
           qgrsx(i,k,4)=qgrs(i,k+niw)
           qgrsx(i,k,5)=qgrs(i,k+nrw)
           qgrsx(i,k,6)=qgrs(i,k+nsw)
           qgrsx(i,k,7)=qgrs(i,k+ngl)
        enddo
        enddo
      elseif (nvdiff.eq.4) then !WSM6 or WSM5 or GFDLMP
        do k=1,levr
        do i=1,lonlens
           qgrsx(i,k,1)=qgrs(i,k)
           qgrsx(i,k,2)=qgrs(i,k+ntoz*levr)
           qgrsx(i,k,3)=qgrs(i,k+ncw)
           qgrsx(i,k,4)=qgrs(i,k+niw)
        enddo
        enddo
      else !nvdiff.eq.3
        do kc1 = 1,nvdiff
           kc = (kc1-1)*levr
        do k=1,levr
        do i=1,lonlens
           qgrsx(i,k,kc1)=qgrs(i,kc+k)
        enddo
        enddo
        enddo
      endif

#ifdef CHS
!<<CCPPv6.0.0
!  in order to achieve heat storage within canopy layer, in the canopy
!    heat torage parameterization the kinematic sensible heat flux
!    (hflx) as surface boundary forcing to the pbl scheme is
!    reduced in a factor of hffac given as a function of surface roughness &
!    green vegetation fraction (zvfun)
!
        if (lheatstrg) then
          do i=1,lonlens
            if (islmsk(i).eq.1) then
              if(hflx(i) > 0.) then
                hffac(i) = h0facu * zvfun(i)
              else
                hffac(i) = h0facs * zvfun(i)
              endif
              hffac(i) = 1. + hffac(i)
              hflx(i) = hflx(i) / hffac(i)
            endif
          enddo
        endif
!>>CCPPv6.0.0
#endif

#ifdef doYSU
      call rysu2d (ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrsx,                  &
     &     prsl_pa,prsi_pa,prslk,prsik(1,1),                            &
     &     igrd1s,lonlens,levr,nvdiff,ncw,delp_pa,                      &
     &     phii,phil,psexp_pa,                                          &
     &     zorl(1,lat),stress,hpbl(1,lat),kpbl,fm1,fh1,                 &
     &     islmsk,hflx,evap,wind,rb,                                    &
     &     dtf,u10m(1,lat),v10m(1,lat),                                 &
     &     swh(1,1,lat),hlw(1,1,lat),xmu,dudt,dvdt,dtdt,dqdtx,          &
     &     dusfc1,dvsfc1,dtsfc1,dqsfc1,gamt,gamq,gust(1,lat),uvdiff )
#endif
#ifdef doMONI_N
       ntcw=1+ntoz+1 ! qt+oz+ql
       call rmoninedmf(igrd1s,lonlens,levr,nvdiff,                       &
     &                ntcw,dvdt,dudt,dtdt,dqdtx,                        &
     &                ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrsx,              &
     &                swh(1,1,lat),hlw(1,1,lat),xmu,prsik(1,1),         &
     &                rb,zorl(1,lat),u10m(1,lat),v10m(1,lat),fm1,fh1,   &
     &                tsea(1,lat),    hflx,evap,stress,wind,kpbl,       &
!    &                tsea(1,lat),qss,hflx,evap,stress,wind,kpbl,       &
     &                prsi_pa,delp_pa,prsl_pa,prslk,phii,phil,          &
     &                dtf,dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl(1,lat),gamt, &
     &                gamq,gust(1,lat),uvdiff )
#endif
#ifdef doMONINP
      call moninp(lonlens,igrd1s,levr,nvdiff,                            &
     &  dvdt,dudt,dtdt,dqdtx,                                            &
     &  ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrsx,u10m(1,lat),v10m(1,lat),     &
     &  prsik(1,1),rb,fm1,fh1,                                           &
     &  tsea(1,lat),    hflx,evap,stress,wind,kpbl,                      &
!    &  tsea(1,lat),qss,hflx,evap,stress,wind,kpbl,                      &
     &  prsi,delp,prsl,prslk,phii,phil,rclx,dtf,                         &
     &  dusfc1,dvsfc1,dtsfc1,dqsfc1,hpbl(1,lat),gamt,gamq,               &
     &  gust(1,lat),uvdiff)
#endif

      if (nvdiff.eq.7) then     !GFDLMP
        do k=1,levr
        do i=1,lonlens
          dqdt(i,k)           = dqdtx(i,k,1)
          dqdt(i,k+ntoz*levr) = dqdtx(i,k,2)
          dqdt(i,k+ncw)       = dqdtx(i,k,3)
          dqdt(i,k+niw)       = dqdtx(i,k,4)
          dqdt(i,k+nrw)       = dqdtx(i,k,5)
          dqdt(i,k+nsw)       = dqdtx(i,k,6)
          dqdt(i,k+ngl)       = dqdtx(i,k,7)
        enddo
        enddo
      elseif (nvdiff.eq.4) then !WSM6 or WSM5
        do k=1,levr
        do i=1,lonlens
          dqdt(i,k)           = dqdtx(i,k,1)
          dqdt(i,k+ntoz*levr) = dqdtx(i,k,2)
          dqdt(i,k+ncw)       = dqdtx(i,k,3)
          dqdt(i,k+niw)       = dqdtx(i,k,4)
        enddo
        enddo
      else
        do kc1 = 1,nvdiff
           kc = (kc1-1)*levr
        do k=1,levr
        do i=1,lonlens
           dqdt(i,kc+k)=dqdtx(i,k,kc1)
        enddo
        enddo
        enddo
      endif

      deallocate(qgrsx,dqdtx)
!cc
      if(lssav) then
      do 340 j=1,lonlens
      dusfc(j,lat)=dusfc(j,lat)+dusfc1(j)*dtf
      dvsfc(j,lat)=dvsfc(j,lat)+dvsfc1(j)*dtf
      dtsfc(j,lat)=dtsfc(j,lat)+dtsfc1(j)*dtf
      dqsfc(j,lat)=dqsfc(j,lat)+dqsfc1(j)*dtf
340   continue
      endif
#ifdef A
      do k=1,levr
      do j=1,lonlens
        i=j+(lat-1)*igrd1s
#endif
#ifdef T
        wt(i,k,4)=dtdt(j,k)-wt(i,k,9)-wt(i,k,10)
#endif
#ifdef Q
        wq(i,k,3)=dqdt(j,k)
#endif
#ifdef U
        wu(i,k,5)=dudt(j,k)
#endif
#ifdef V
        wv(i,k,5)=dvdt(j,k)
#endif
#ifdef A
      enddo
      enddo
#endif
! end of PBL doMONINP
!
!ccj

#ifdef wind_fitch

      do nw=1,fi_nt
        if ( wjval(nw).eq. lat )  then
        wi(nw)=me
        call dragforce(phii,phil,ugrs(1,1),vgrs(1,1),rdelx,dtf,         &
     &                dudt,dvdt,igrd1s,lonlens,levr,nw)
        endif
      enddo

#endif

#ifdef doGWDPS
!     Dynamically initialize orographic variables from hprime
!     depending on the value of nmtvr
      do i = 1, lonlens
        var(i) = hprime(i,lat,1)
      enddo
!
      oc     = 0.0
      oa4    = 0.0
      clx    = 0.0
      theta  = 0.0
      gamma  = 0.0
      sigma  = 0.0
      elvmax = 0.0

! SCC 10/02/05
      if (nmtvr .eq. 6) then
       do i=1,lonlens
         oc(i) = hprime(i,lat,2)
       enddo
       do k = 1, 4
         do i=1,lonlens
           oa4(i,k) = hprime(i,lat,k+2)
           clx(i,k) = 0.0
         enddo
       enddo
! --- lm mb (*j*) 1997 gwd
      elseif(nmtvr .eq. 10) then
       do i=1,lonlens
         oc(i) = hprime(i,lat,2)
       enddo
       do k = 1, 4
         do i=1,lonlens
           oa4(i,k) = hprime(i,lat,k+2)
           clx(i,k) = hprime(i,lat,k+6)
         enddo
       enddo
! --- lm mb (*j*)
      elseif(nmtvr .eq. 14) then
! --- get the kim fields (until this is changed)
       do i=1,lonlens
         oc(i) = hprime(i,lat,2)
       enddo
       do k = 1, 4
         do i=1,lonlens
           oa4(i,k) = hprime(i,lat,k+2)
           clx(i,k) = hprime(i,lat,k+6)
         enddo
       enddo
       do i=1,lonlens
          theta(i)  = hprime(i,lat,11)
          gamma(i)  = hprime(i,lat,12)
          sigma(i)  = hprime(i,lat,13)
          elvmax(i) = hprime(i,lat,14)
       enddo
      endif
! SCC 10/02/05
!
      lonfequiv=int(192.0*200.e3/rdelx)
!
!     GFS version of GWD parameterization scheme
      call rgwdps(lonlens,igrd1s,igrd1s,levr,dvdt,dudt,dtdt,   &
                 ugrs(1,1),vgrs(1,1),tgrs(1,1),qgrs(1,1),     &
                 kpbl,prsi_pa,delp_pa,prsl_pa,prslk,phii,phil,dtp,kdt, &
                 var,oc,oa4,clx,theta,sigma,gamma,elvmax,     &
                 dusfcg,dvsfcg,grav,cp,rd,rv,lonfequiv,nmtvr,cdmbgwd)

!       print *,'aft gwdps, me=',me,'lat=',lat
!c lonfequiv is an equivalent number of lonf of mrf in rsm
! for mrf, cleff = 0.5e-5 * sqrt(float(imx)/192.0) !  this is inverse of cleff!
! where imx=lonf
!  scale cleff between im=384*2 and 192*2 for t126/t170 and t62
!ccc
!
      if(lssav) then
      do 341 j=1,lonlens
      dugwd(j,lat)=dugwd(j,lat)+dusfcg(j)*dtf
      dvgwd(j,lat)=dvgwd(j,lat)+dvsfcg(j)*dtf
341   continue
      endif
!
#ifdef A
      do k=1,levr
      do j=1,lonlens
        i=j+(lat-1)*igrd1s
#endif
#ifdef U
        wu(i,k,7)=dudt(j,k)-wu(i,k,5)
#endif
#ifdef V
        wv(i,k,7)=dvdt(j,k)-wv(i,k,5)
#endif
#ifdef A
      enddo
      enddo
#endif
#endif
! end of doGWDPS
#ifdef NONHYD
!
! ustr = u / m    vstr = v /m
! dwdt=m**2(d ustr/dt*dzdx + d vstr/dt*dzdy)
! dwdt=m(dudt*dgzdx+dvdt*dgzdy)/g
!
      do j=1,lonlens
      dwdt(j,1)=dudt(j,1)*gzdx(j,lat)+dvdt(j,1)*gzdy(j,lat)
      dwdt(j,1)=dwdt(j,1)*xm(j,lat) / grav
      enddo
#endif
!
!.......................................................
!...      get vertical motion (cb/sec) in vvel
!.......................................................
#ifdef NONHYD
      do k=1,levr
      do j=1,lonlens
        wmean=-0.5*(wgrs(j,k)+wgrs(j,k+1))
#ifdef HYDPHY
!hmhj vvel=-(pg)/(rT) dz/dt  201812

        tempv= tgrh(j,k)*(1.+fv*max(qgrs(j,k),qmin))
        vvel(j,k)=wmean*prsl(j,k)*grav/(rd*tempv)
#else
        vvel(j,k)=wmean*prsl(j,k)*grav/(rd*tgrs(j,k))
#endif
      enddo
      enddo
#else
#ifndef SIGP
      call omegast3(lonlens,igrd1s,levr,                                 &  
     &            pphigrs(1),plamgrs(1),                                 &
     &            ugrs(1,1),vgrs(1,1),                                   &
     &            xgrs(1,1),del,rbs2,vvel,                               &
     &            psexp,slo)
#else
      call omegast3_sp(lonlens,igrd1s,levr,                              &  
     &            pphigrs(1),plamgrs(1),                                 &
     &            ugrs(1,1),vgrs(1,1),                                   &
     &            xgrs(1,1),rbs2,vvel,                                   &
     &            psexp,sio,sig(k,1))
#endif
#endif
!       print *,'aft omegac, me=',me,'lat=',lat
!
!   w2 is used to store wind speed
!
!     do k=1,levr
!      do j=1,lonlens
!         w2(j,k)=rcs*sqrt(ugrs(j,k)**2+vgrs(j,k)**2)
!     enddo
!     enddo
!
!yj: for wsm--------------------------------------------------------->>>
!#ifdef NONHYD
!   do k = 1,levr 
!     do i = 1,lonlens
!        tempv= tgrh(j,k)*(1.+fv*max(qgrs(j,k),qmin))
!        delz(i,k) = (rd*tempv)/grav*(log(si(k))-log(si(k+1)))
!     enddo
!   enddo
!#else
   do k = 1,levr 
     do i = 1,lonlens
       zl(i,k)=phil(i,k)/grav
       zi(i,k)=phii(i,k)/grav
     enddo
   enddo
   do i = 1,lonlens
     zi(i,levr+1)=phii(i,levr+1)/grav
   enddo
   do k = 1,levr 
     do i = 1,lonlens
       delz(i,k) = zi(i,k+1)-zi(i,k)
!        if (delz(i,k).le.1.) then
!          if (me.eq.0) print*,'delz too small,i,k=',delz,i,k
!          call abort
!        endif
     enddo
   enddo
!#endif

   do k = 1,levr 
     do i = 1,lonlens
!#ifdef NONHYD
!       ww(i,k) = 0.5*(wgrs(j,k)+wgrs(j,k+1))
!#else
      !yj: delp=p(z-1)-p(z+1), is positive
       ww(i,k) = -vvel(i,k)*delz(i,k)/delp(i,k)
!yj       ww(i,k) = -xgrs(i,k)*delz(i,k)/delp(i,k)
!#endif
     enddo
   enddo
!yj: for wsm---------------------------------------------------------<<<
!
!testjun
      do 350 k=1,levr
      do 350 j=1,lonlens
      gt0(j,k)=tgrs(j,k)+dtdt(j,k)*dtp
      gu0(j,k)=ugrs(j,k)+dudt(j,k)*dtp
      gv0(j,k)=vgrs(j,k)+dvdt(j,k)*dtp
#ifdef NONHYD
#ifdef HYDPHY
      gh0(j,k)=tgrh(j,k)+dtdt(j,k)*dtp
#endif
#endif
  350 continue

#ifdef NONHYD
      do k=1,levr+1
        do j=1,lonlens
          gw0(j,k)=wgrs(j,k)+dwdt(j,k)*dtp
        enddo
      enddo
#endif
!ccj
!testjun
!      call maxmin(dqdt,lonlens*levh,1,1,1,'dqdt in rbphys')
      do 360 k=1,levh
      do 360 j=1,lonlens
      gq0(j,k)=qgrs(j,k)+dqdt(j,k)*dtp
  360 continue
!      do 360 k=1,levh
!      do 360 j=1,lonlens
!      gq0(j,k)=qgrs(j,k)
!  360 continue
!
!  ozone physics
!
#ifdef doOZPHYS
      if(ntoz .eq. 1)then
        call ozphys(lat,lonlens,igrd1s,dtp,prsl,                         &  
     &              gq0(1,1+levr))
      endif
#endif
!cc
!
!  call to gwater
!.....
!
#ifndef SIGP
#ifndef NONHYD
      call rget_phi(lonlens,igrd1s,levr,cp,fv,gt0(1,1),gq0(1,1),          &
     &             prsik,prslk,phii,phil)
#else
#ifdef HYDPHY
      call rget_phi(lonlens,igrd1s,levr,cp,fv,gh0(1,1),gq0(1,1),          &
     &             prsik,prslk,phii,phil)
#else
      call rget_phi(lonlens,igrd1s,levr,cp,fv,gt0(1,1),gq0(1,1),          &
     &             prsik,prslk,phii,phil)
#endif
#endif
#else
#ifndef NONHYD
      call get_phi_sp(lonlens,igrd1s,levr,rd,fv,gt0(1,1),gq0(1,1),       &
     &             prsi,phii,phil)
#else
#ifdef HYDPHY
      call get_phi_sp(lonlens,igrd1s,levr,rd,fv,gh0(1,1),gq0(1,1),       &
     &             prsi,phii,phil)
#else
      call get_phi_sp(lonlens,igrd1s,levr,rd,fv,gt0(1,1),gq0(1,1),       &
     &             prsi,phii,phil)
#endif
#endif
#endif
!
!ccj
      do k=1,levr
        do i=1,lonlens
          clw(i,k,1) = 0.0
          clw(i,k,2) = -999.9
        enddo
      enddo
#if defined(SCASAS) || defined(SAMFSAS)
         cnvc(:,:)  = 0.0
         cnvw(:,:)  = 0.0
#endif
       kuo(:)   = 0
       ktop(:)  = 1
       kbot(:)  = levr
!
      if(ncld.gt.0) then
        ncw = (1+ntoz)*levr
      else
        ncw = 0
      endif
!cc
#ifdef A
      do k=1,levr
      do j=1,lonlens
        i=j+(lat-1)*igrd1s
#endif
#ifdef T
        wt(i,k,6)=gt0(j,k)
#endif
#ifdef Q
        wq(i,k,5)=gq0(j,k)
#endif
#ifdef A
      enddo
      enddo
#endif

!ccj
!
!     calling convective parameterization
!
      if (ncld .gt. 0) then
!
        if (nmmiph .eq. 3) then    ! call brad ferrier microphysics
!
!***          algorithm to separate different hydrometeor species
!
          do k=1,levr
            kk=k+ncw
            do i=1,lonlens
              wc     = gq0(i,kk)
              qi     = 0.
              qr     = 0.
              qw     = 0.
              f_ice  = max(0.0, min(1.0, phy_f3dv(i,k,1)))
              f_rain = max(0.0, min(1.0, phy_f3dv(i,k,2)))
!
              qi = f_ice*wc
              qw = wc-qi
              if (qw .gt. 0.0) then
                qr = f_rain*qw
                qw = qw-qr
              endif
!
!             if (f_ice .ge. 1.) then
!               qi = wc
!             else if (f_ice .le. 0.) then
!               qw = wc
!             else
!               qi = f_ice*wc
!               qw = wc-qi
!             endif
!
!             if (qw.gt.0. .and. f_rain.gt.0.) then
!               if (f_rain .ge. 1.) then
!                 qr = qw
!                 qw = 0.
!               else
!                 qr = f_rain*qw
!                 qw = qw-qr
!               endif
!             endif
!
              qr_col(i,k) = qr
!             clw(i,k)    = qi + qw
              clw(i,k,1)  = qi
              clw(i,k,2)  = qw
!
!
!***          array to track fraction of "cloud" in the form of ice
!
!             if (qi+qw .gt. epsq) then
!               fc_ice(i,k) = qi / (qi+qw)
!             else
!               fc_ice(i,k) = 0.
!             endif

            enddo
          enddo
        elseif (nmmiph.eq.2) then  !Zhao-Carr MP
            do k=1,levr
              do i=1,lonlens
                clw(i,k,1) = gq0(i,ncw+k)
              enddo
            enddo
        elseif (nmmiph.eq.4) then !WSM3
            do k=1,levr
              do i=1,lonlens
                clw(i,k,1) = gq0(i,ncw+k)
              enddo
            enddo
        elseif (nmmiph.eq.5) then !WSM5
            do k=1,levr
              do i=1,lonlens
                clw(i,k,1) = gq0(i,niw+k)  !ice
                clw(i,k,2) = gq0(i,ncw+k)  !water
              enddo
            enddo
        elseif (nmmiph.eq.6 .or. nmmiph.eq.15 .or.                     & !WSM6
                nmmiph.eq.11 .or. nmmiph.eq.12 .or. nmmiph.eq.13) then !GFDLMP
            do k=1,levr
              do i=1,lonlens
                clw(i,k,1) = gq0(i,niw+k)  !ice
                clw(i,k,2) = gq0(i,ncw+k)  !water
              enddo
            enddo
        endif    ! endif clw defined
      endif      ! endif (ncld .gt. 0)

#ifdef doCONCNV
!!
#ifndef RAS
!4x4
!        if ( me.eq.0 ) then
!         print *,'me=0,lat=',lat,'gq0(1:5,13)=',gq0(1:5,13),             &
!     &                         'gt0(1:5,13)=',gt0(1:5,13),               &
!     &   'delp=',delp(1:5,1),'prsl=',prsl(1:5,1),'psexp=',psexp(1:5),    &
!     &   'phil=',phil(1:5,1),'clw=',clw(1:5,1,1),clw(1:5,1,2),'rcs=',    &
!     &   rcs(1:5),'cld1s=',cld1d(1:5),'rain1=',rain1(1:5),'kbot=',       &
!     &   kbot(1:5),'ktop=',ktop(1:5),'kuo=',kuo(1:5),'slmsk=',           &
!     &   slmsk(1:5,lat),'vvel=',vvel(1:5,1),'xkt2=',xkt2(1:5)            
!       endif
#ifdef SAS
      do i=1,lonlens
!      call sascnv(lonlens,igrd1s,levr,                                   &  
!     &            jcapr,dtp,delp,prsl,psexp,phil,                        &
!     &            clw,gq0(1,1),gt0,gu0,gv0,rcs,cld1d,                    &
!     &            rain1,kbot,ktop,kuo,slmsk(1,lat),                      &
!     &            vvel,xkt2,ncld)
      call rsascnv(1,igrd1s,levr,                                         &  
     &            jcapr,dtp,delp(i,1),prsl(i,1),psexp(i),phil(i,1),      &
     &            clw(i,1,1),gq0(i,1),gt0(i,1),gu0(i,1),gv0(i,1),rcs(i), &
     &            cld1d(i),                                              &
     &            rain1(i),kbot(i),ktop(i),kuo(i),slmsk(i,lat),          &
     &            vvel(i,1),xkt2(i),ncld)
      enddo
#endif      
!end of ifdef SAS
#ifdef NSAS      
      call sascnv_nsas (lonlens,igrd1s,levr ,jcapr,dtp                  &
     &              ,delp    ,prsl  ,psexp,phil ,clw                     &
     &              ,gq0(1,1),gt0   ,gu0  ,gv0  ,rcs                     &
     &              ,cld1d   ,rain1 ,kbot ,ktop ,kuo                     &
     &              ,slmsk(1,lat)  ,vvel  ,xkt2 ,ncld )
#endif
#ifdef SCASAS
!xb110> 20190731
      call mfdeepcnv (lonlens,igrd1s,levr ,dtp ,delp                    &
     &               ,prsl   ,psexp ,phil ,clw ,gq0(1,1)                &
     &               ,gt0    ,gu0   ,gv0  ,vvel,cld1d                   &
     &               ,rain1  ,kbot  ,ktop ,kuo ,islmsk                  &
     &               ,ncld   ,cnvw  ,cnvc ,sigmagfm, rdelx)
!xb110> 20190731
!      call mfdeepcnv_v2 (lonlens,igrd1s,levr ,dtp ,delp                 &
!     &               ,prsl   ,psexp ,phil ,clw ,gq0(1,1)                &
!     &               ,gt0    ,gu0   ,gv0  ,vvel,cld1d                   &
!     &               ,rain1  ,kbot  ,ktop ,kuo ,slmsk(1,lat)            &
!     &               ,ncld   ,cnvw  ,cnvc ,sigmagfm, rdelx)
!xb110< 20190731
#endif
#ifdef SAMFSAS
#ifdef LTN
!lightning, keep temprature before wet physics and invert the vertical 
!           direction (upside down)
      do k=1,levr
        ud_ztenh(:,levr-k+1)=gt0(:,k)
        ud_zqenh(:,levr-k+1)=qgrs(:,k)
      enddo
      do k=1,levr
      do i=1,lonlens
        ud_rho(i,levr-k+1)= prsl_pa(i,k) / (rd *  gt0(i,k) )
      enddo
      enddo
#endif

!
!xb110> 20190911, the scale-awareness SAS scheme in FV3
       call rsamfdeepcnv (lonlens,igrd1s,levr,dtp,delp,prsl,psexp,phil    &
     &    ,clw,gq0(1,1),gt0,gu0,gv0,cld1d,rain1,kbot,ktop,kuo,islmsk     &
     &    ,rdelx,vvel,ncld,cnvw,cnvc                                     &
#ifdef LTN
     &                              ,snow_flx,ptu,pqu                    &
#endif
     &                              )
!xb110< 20190911
#ifdef LTN
!lightning----------------->
!
!need to check:
!! WARNING: in lightning_ec, top layer is 1, bottom layer is levr       !!
!! WARNING: in lightning_ec, flash unit: km-2 day-1 but not sure        !!
!                            let's supposed the unit is km-2 s-1 here   !!
!invert all variables (upside down)
      do k=1,levr
        ud_ptu(:,levr-k+1)=ptu(:,k)
        ud_pqu(:,levr-k+1)=pqu(:,k)
        ud_cnvw(:,levr-k+1)=cnvw(:,k)
        ud_phil(:,levr-k+1)=phil(:,k)
        ud_prsl_hpa(:,levr-k+1)=prsl(:,k)*10.0
        ud_snow_flx(:,levr-k+1)=snow_flx(:,k)
      enddo
      do i=1,lonlens
        ud_kbot(i)=levr-kbot(i)+1
        ud_ktop(i)=levr-ktop(i)+1
      enddo
!
!      call lightning_ec(lonlens,igrd1s,levr,ud_ptu,ud_pqu                &
!                       ,ud_ztenh,ud_zqenh                                &
!                       ,ud_cnvw,ud_rho,ud_kbot,ud_ktop                   &
!                       ,ud_phil,ud_prsl_hpa,ud_snow_flx                  &
!                       ,flash,islmsk,kuo)
! for output
!      if(lssav) then
!        do j = 1, lonlens 
!          cflash(j,lat)=cflash(j,lat)+flash(j)*dtf
!        enddo
!      endif
!lightning-----------------<
#endif
#endif
!xb110> new Tiedtke scheme developed by Bechtold et al.2014
!and further modified to decrease the downdraft entrinment with
!incraseing bouyancy.
#ifdef NTDK
      do k=1,levr !to invert the vertical coordinate for nTDK scheme
         kc=levr-k+1
         do i = 1,lonlens
            rphil(i,kc) = phil(i,k)
            rdudt(i,kc) = dudt(i,k)
            rdvdt(i,kc) = dvdt(i,k)
            rdtdt(i,kc) = dtdt(i,k)
            rdqdt(i,kc) = dqdt(i,k)
            rgu0(i,kc)  = gu0(i,k)
            rgv0(i,kc)  = gv0(i,k)
            rgt0(i,kc)  = gt0(i,k)
            rgq0(i,kc)  = gq0(i,k)
            rvvel(i,kc)  = vvel(i,k)
            rprsl(i,kc) = prsl(i,k)
            if (ncld.le.2) then !Zhao-Carr and WSM3
              do j = 1,2
              rclw(i,kc,j)  = clw(i,k,j)  !cloud water
              enddo
            else                ! WSM5,WSM6,GFDLMP
              do j = 1,2
                rclw(i,kc,j)  = clw(i,k,j) !cloud ice and water
              end do
            endif
         end do !for i loop
      end do !for k loop

      do k = 1,levr+1
         kc=(levr+1)-k+1
         do i = 1,lonlens
            rprsi(i,kc) = prsi(i,k)
         end do
      end do

      call cumastr_driv_n                                               &
     &        (lonlens,igrd1s,levr,dtp,rphil,rdudt,rdvdt,rdtdt,rdqdt    &
     &        ,rgu0,rgv0,rgt0,rgq0,rain1,rvvel,evap,rkbot,rktop,ncld    &
     &        ,rprsl,psexp,islmsk,hflx,kuo,rprsi,rdelx,rclw,hpbl(1,lat))

      do k = 1,levr  ! to invert back to the original vertical axis
         kc = levr-k+1
         do i = 1,lonlens
            gu0(i,k)   = rgu0(i,kc)
            gv0(i,k)   = rgv0(i,kc)
            gt0(i,k)   = rgt0(i,kc)
            gq0(i,k)   = rgq0(i,kc)
            if (ncld.le.2) then ! Zhao-Carr and WSM3
              do j = 1,2
              clw(i,k,j) = rclw(i,kc,j)    !cloud water
              enddo
            else                ! WSM5,WSM6,GFDLMP
              do j = 1,2
                clw(i,k,j) = rclw(i,kc,j) !cloud ice and water
              end do
            endif
         end do !for i loop
       end do !for k loop

       do i = 1,lonlens
         kbot(i) = levr-rkbot(i)+1
         ktop(i) = levr-rktop(i)+1
       end do


#endif
!xb110<

!        if ( me.eq.0 ) then
!         print *,'me=0,lat=',lat,'gq0(1:5,13)=',gq0(1:5,13),             &
!     &    'gt0(1:5,13)=',gt0(1:5,13)   
!        endif
!
#else
      call rascnv(lonlens,igrd1s,levr, dtp, ncrnd, rannum                &
     &,           gt0,    gq0,   gu0,    gv0, clw, 0                     &
     &,           prsi,   prsl,   prsik,  prslk, phil,  phii             &
     &,           kpbl,   cd,     rain1,  kbot,  ktop,  kuo              &
     &,           phy_f2dv(1,num_p2d), flipv, cb2mb                      &
     &,           me, 1, 1, rarea, lprnt)
      cld1d = 0.
#endif
!end of ifndef RAS
#ifdef A
      do k=1,levr
      do j=1,lonlens
        i=j+(lat-1)*igrd1s
#endif
#ifdef T
        wt(i,k,6)=(gt0(j,k)-wt(i,k,6))/dt2
#endif
#ifdef Q
        wq(i,k,5)=(gq0(j,k)-wq(i,k,5))/dt2
#endif
#ifdef T
        wt(i,k,7)=gt0(j,k)
#endif
#ifdef Q
        wq(i,k,6)=gq0(j,k)
#endif
#ifdef A
      enddo
      enddo
#endif
#endif
! end of doCONCNV
!-------yj20191030----end deep conv

#ifdef doSCONCNV
#if defined(RAS) || defined(SAS)     
! This will run for SAS and RAS, but not newSAS(NSAS)
!
      call shalcv(lonlens,igrd1s,levr,dtp,delp,prsi,prsl,                &
     &              prslk,kuo,gq0(1,1),gt0(1,1),dpshc)
#endif     
#ifdef NSAS      
      call shalcv_nsas (lonlens,igrd1s,levr ,jcapr,dtp                   &
     &                 ,delp   ,prsl  ,psexp,phil ,clw                   &
     &                 ,gq0    ,gt0   ,gu0  ,gv0  ,rain2                 &
     &                 ,kbot   ,ktop  ,kuo  ,slmsk(1,lat),vvel           &
     &                 ,ncld   ,hpbl(1,lat)  ,hflx ,evap )
!
#endif      
#ifdef SCASAS
!xb110> 20190731
      call mfshalcnv (lonlens,igrd1s,levr ,dtp                          &
     &               ,delp   ,prsl  ,psexp,phil ,clw                    &
     &               ,gq0    ,gt0   ,gu0  ,gv0  ,rain2                  &
     &               ,kbot   ,ktop  ,kuo  ,islmsk,vvel                  &
     &               ,ncld   ,hpbl(1,lat) ,cnvw ,cnvc ,sigmagfm, rdelx)
!xb110> 20190731
!      call mfshalcnv_v2 (lonlens,igrd1s,levr ,dtp                       &
!     &               ,delp   ,prsl  ,psexp,phil ,clw                    &
!     &               ,gq0    ,gt0   ,gu0  ,gv0  ,rain2                  &
!     &               ,kbot   ,ktop  ,kuo  ,slmsk(1,lat),vvel            &
!     &               ,ncld   ,hpbl(1,lat) ,cnvw ,cnvc ,sigmagfm, rdelx)
!xb110< 20190731
!
#endif
#if defined(SAMFSAS) || defined(NTDK)
!#if defined(SAMFSAS)
!xb110> 20190911, the scale-awareness SAS scheme in FV3
      call rsamfshalcnv(lonlens,igrd1s,levr,dtp,delp,prsl,psexp,phil,clw  &
     &    ,gq0,gt0,gu0,gv0,rain2,kbot,ktop,kuo,islmsk,rdelx,vvel         &
     &    ,ncld,hpbl(1,lat))
!xb110< 20190911
!
#endif
!
#ifdef A
      do k=1,levr
      do j=1,lonlens
        i=j+(lat-1)*igrd1s
#endif
#ifdef T
        wt(i,k,7)=(gt0(j,k)-wt(i,k,7))/dt2
#endif
#ifdef Q
        wq(i,k,6)=(gq0(j,k)-wq(i,k,6))/dt2
#endif
#ifdef T
        wt(i,k,8)=gt0(j,k)
#endif
#ifdef Q
        wq(i,k,7)=gq0(j,k)
#endif
#ifdef A
      enddo
      enddo
#endif
#endif
! end of doSCONCNV
!j
      do k=1,levr
        do i=1,lonlens
          if (clw(i,k,2) .le. -999.0) clw(i,k,2) = 0.0
        enddo
      enddo
!
      if (ncld .gt. 0) then
        if (nmmiph .eq. 3) then    ! call brad ferrier microphysics
!
!***          extract cloud water & ice from fc_ice
!
          do k=1,levr
            kk=k+ncw
            do i=1,lonlens
!             qi = clw(i,k)*fc_ice(i,k)
!             qw = clw(i,k) - qi
!
              qi = clw(i,k,1)
              qw = clw(i,k,2)
!
!***          algorithm to combine different hydrometeor species
!
!             gq0(i,kk) = max(epsq, qi+qw+qr_col(i,k))
              gq0(i,kk) = qi + qw + qr_col(i,k)
              if (qi .le. epsq) then
                phy_f3dv(i,k,1) = 0.
              else
                phy_f3dv(i,k,1) = qi/gq0(i,kk)
              endif
              if (qr_col(i,k) .le. epsq) then
                phy_f3dv(i,k,2) = 0.
              else
                phy_f3dv(i,k,2) = qr_col(i,k) / (qw+qr_col(i,k))
              endif
            enddo
          enddo
        elseif (nmmiph .eq. 2) then  !Zhao-Carr MP
            do k=1,levr
!             kk=k+ncw
              do i=1,lonlens
!               gq0(i,kk) = clw(i,k) + gq0(i,kk)
                gq0(i,ncw+k) = clw(i,k,1) + clw(i,k,2)
              enddo
            enddo
        elseif (nmmiph.eq.4) then   ! WSM3
            do k=1,levr
              do i=1,lonlens
                gq0(i,ncw+k) = clw(i,k,1) + clw(i,k,2)
              enddo
            enddo
        elseif (nmmiph.eq.5) then   ! WSM5
            do k=1,levr
              do i=1,lonlens
                gq0(i,niw+k) = clw(i,k,1) !ice
                gq0(i,ncw+k) = clw(i,k,2) !water
              enddo
            enddo
        elseif (nmmiph.eq.6 .or. nmmiph.eq.15 .or.                      &  ! WSM6
                nmmiph.eq.11 .or. nmmiph.eq.12 .or. nmmiph.eq.13) then  ! GFDLMP
            do k=1,levr
              do i=1,lonlens
                gq0(i,niw+k) = clw(i,k,1) !ice
                gq0(i,ncw+k) = clw(i,k,2) !water
              enddo
            enddo

        endif   !endif nmmiph
      else
        do k=1,levr
          do i=1,lonlens
              clw(i,k,1) = clw(i,k,1) + clw(i,k,2)
          enddo
        enddo
      endif    !endif (ncld .gt. 0)
!
      if(lssav) then
      do  j=1,lonlens
         cldwrk(j,lat) = cldwrk(j,lat) + cld1d(j) * dtf
      enddo
      endif
!
      do 235 j=1,lonlens
      rain1(j) = rain1(j)+rain2(j)
      rainc(j)=frain*rain1(j)
      if(lssav) bengsh(j,lat)=bengsh(j,lat)+rainc(j)
235   continue
!
!ccj
!
      call cnvc90(clstp,lonlens,igrd1s,rainc,kbot,ktop,                  &   
     &           levr,prsi,                                              &
     &           acv(1,lat),acvb(1,lat),acvt(1,lat),                     &
     &           cv(1,lat), cvb(1,lat), cvt(1,lat))
!
!ccj
#ifdef doLRGSCL
      ncw = (1+ntoz)*levr
      if (ncld .eq. 0) then
        call lrgscl(lonlens,igrd1s,levr,dtp,gt0(1,1),gq0(1,1),           &   
     &                 prsl,delp,prslk,rain1,clw(1,1,1))
      elseif (ncld .eq. 1) then
!
!       to call moist convective adjustment above pbl
!
!       call mstcnv(lonlens,igrd1s,levr,dtp,gt0,gq0,prsl,del,prslk,rain1
!    &,                                                 lprnt)
!       rainc = rainc + frain * rain1
!       if(lssav) bengsh(ilon) = bengsh(ilon)
!    &                              + rain1 * frain
!
!       moist convective adjustment over
!
        nlonsequiv=int((78500.0/rdelx)*512.) !nlonsequiv is an equivalence of nlons(i) in mrf
        do i=1,lonlens
          work1(i) = (log(1.0 / (rcs(i)*nlonsequiv)) - dxmin) * dxinv
          work2(i) = 1.0 - work1(i)
        enddo
!
        invrsn = .false.
        kinver = 0
        do k=1,levr/2
          do i=1,lonlens
            if (prsl(i,k) .gt. 60.0 .and. (.not. invrsn(i))) then
               ptem = gt0(i,k+1) - gt0(i,k)
               if (ptem .gt. 0.2 .and. gt0(i,k) .gt. 278.0) then
                 invrsn(i) = .true.
                 kinver(i) = k
               endif
            endif
          enddo
        enddo
!      if ( levr .gt. 100) then
!--------------------------------------------------------------------------
!
      if (nmmiph .eq. 3) then    ! call brad ferrier microphysics
!
        do k=1,levr
          do i=1,lonlens
            kk = kpbl(i)
            if (k .le. kk) then
              ptem  = 0.0
              ptem1 = rhbbot
            else
 !            ptem  = (rhbtop-rhbbot) / (prsik(i,levr+1)-prsik(i,kk))
              ptem  = (rhbtop-rhpbl) / (prsik(i,levr+1)-prsik(i,kk))
              ptem1 = rhpbl
            endif
            rhc(i,k) = ptem1 + ptem * (prslk(i,k)-prslk(i,kk))
            rhc(i,k) = 0.999 * work1(i) + rhc(i,k) * work2(i)
          enddo
        enddo
!
        do k=1,levr/2
          do i=1,lonlens
!           if (k .le. kinver(i) .and. slmsk(i,lat) .eq. 0.0) then
            if (k .le. kinver(i)) then
              rhc(i,k) = min(0.85, rhc(i,k))
!             rhc(i,k) = 0.85
            endif
          enddo
        enddo
!
!       if (me .eq. 0 .and. kdt .eq. 1 .and. ilon .eq. 1)
!    &      print *,' rhc=',rhc
        do i=1,lonlens
!         xncw(i) = 100.0 * work1(i) + 10.0 * work2(i)
!         xncw(i) = 200.0 * work1(i) + 30.0 * work2(i)
!         xncw(i) = 250.0 * work1(i) + 25.0 * work2(i)
!         xncw(i) = 250.0 * work1(i) + 20.0 * work2(i) ! 20021217
          xncw(i) = 100.0 * work1(i) + 20.0 * work2(i)
!!        xncw(i) = 250.0 * work1(i) + 15.0 * work2(i)
!         xncw(i) = 250.0 * work1(i) + 10.0 * work2(i)
!         xncw(i) = 250.0 * work1(i) + 50.0 * work2(i)
        enddo
!       if (me .eq. 0 .and. kdt .eq. 1 .and. ilon .eq. 1)
!    &          print *,' xncw=',xncw
!
        call gsmdrive(lonlens, igrd1s, levr, dtp, prsl, delp             & 
     &,               gt0, gq0(1,1), gq0(1,1+ncw), slmsk(1,lat)          &
     &,               phy_f3dv(1,1,1),  phy_f3dv(1,1,2)                  &
     &,               phy_f3dv(1,1,3), rain1, sr, grav                   &
     &,               hvap, hsub, cp, rhc, xncw, me, lprnt,ipr)
!
      elseif (nmmiph .eq. 2) then ! call zhao/carr/sundqvist microphysics
!
        do k=1,levr
          do i=1,lonlens
            ptem   = (rhztop-rhzbot) / (prslk(i,levr)-prslk(i,1))
            rhc(i,k) = rhzbot + ptem * (prslk(i,k)-prslk(i,1))
            rhc(i,k) = 0.999 * work1(i) + rhc(i,k) * work2(i)
            rhc(i,k) = min(0.999,rhc(i,k))
          enddo
        enddo
!
!   the following commented by moorthi on 11/01/2002
!       do k=1,levr/2
!         do i=1,lonlens
!           if (k .le. kinver(i) .and. slmsk(i,lat) .eq. 0.0) then
!           if (k .le. kinver(i)) then
!             rhc(i,k) = 0.80
!           endif
!         enddo
!       enddo
!       if (me .eq. 0 .and. kdt .eq. 1 .and. ilon .eq. 1)
!    &      print *,' rhc=',rhc
!
#if defined(RAS) || defined(SAS)     
        call rgscond(lonlens, igrd1s, levr,                               &    
     &              dtp, prsl, psexp,                                    &
     &              gq0(1,1), gq0(1,1+ncw), gt0,                         &
     &              num_p3d,num_p2d,phy_f3dv(:,:,1:num_p3d),             &
     &              phy_f2dv(:,1:num_p2d),                               & 
!     &              phy_f3dv(1,1,1), phy_f3dv(1,1,2), phy_f2dv(1,1),     &
!     &              phy_f3dv(1,1,3), phy_f3dv(1,1,4), phy_f2dv(1,2),     &
     &              rhc,lprnt)
        call rprecpd(lonlens, igrd1s, levr,                               &
     &              dtp, delp, prsl, psexp,                              &
     &              gq0(1,1), gq0(1,1+ncw), gt0, rain1, rhc, lprnt)
#endif
#if defined(NSAS)        
        call gscond_nsas(lonlens, igrd1s, levr,                          &    
     &              dtp, prsl, psexp,                                    &
     &              gq0(1,1), gq0(1,1+ncw), gt0,                         &
     &              num_p3d,num_p2d,phy_f3dv(:,:,1:num_p3d),             &
     &              phy_f2dv(:,1:num_p2d),                               & 
!     &              phy_f3dv(1,1,1), phy_f3dv(1,1,2), phy_f2dv(1,1),     &
!     &              phy_f3dv(1,1,3), phy_f3dv(1,1,4), phy_f2dv(1,2),     &
     &              rhc,lprnt)
        call precpd_nsas(lonlens, igrd1s, levr,                           &
     &              dtp, delp, prsl, psexp,                              &
     &              gq0(1,1), gq0(1,1+ncw), gt0, rain1, rhc, lprnt)
#endif        
#if defined(SCASAS) || defined(SAMFSAS) || defined(NTDK)
        call gscond_sca(lonlens, igrd1s, levr,                           &    
     &              dtp, prsl, psexp,                                    &
     &              gq0(1,1), gq0(1,1+ncw), gt0,                         &
     &              num_p3d,num_p2d,phy_f3dv(:,:,1:num_p3d),             &
     &              phy_f2dv(:,1:num_p2d),                               & 
!     &              phy_f3dv(1,1,1), phy_f3dv(1,1,2), phy_f2dv(1,1),     &
!     &              phy_f3dv(1,1,3), phy_f3dv(1,1,4), phy_f2dv(1,2),     &
     &              rhc,lprnt)
        call precpd_sca(lonlens, igrd1s, levr,                           &
     &              dtp, delp, prsl, psexp,                              &
     &              gq0(1,1), gq0(1,1+ncw), gt0, rain1, rhc, lprnt)
#endif
      endif
!ccc
!
     else  !ncld>1
      if (nmmiph .eq. 4) then
#ifdef NEW_WSM3
          im = lonlens
          ix = igrd1s
          ims = 1 ; ime = ix ; kms = 1 ; kme = levr ; its = 1 ; ite = im ; kts = 1 ; kte = levr
          ntcw=1+ncw
          ntrw=1+ncw+1*levr

           call wsm3(gt0(1,1), phii(1:ix,1:levr+1),                              &
                                gq0(1:ix,1:levr),                                &
                                gq0(1:ix,ntcw:ntcw+levr-1),                      &
                                gq0(1:ix,ntrw:ntrw+levr-1),                      &
                                ww(1,1),                                         &
                                ! rain1 will be zero out and calculate in WSM3
                                prsl_pa(1,1), delz, dtp, rain1,                  &
                                sr,                                              &
                                !re_cloud, re_ice, re_snow ! for radiation
                                phy_f3dv(:,:,1),phy_f3dv(:,:,2),phy_f3dv(:,:,3), &
                                ims,ime, kms,kme,                                &
                                its,ite, kts,kte                                 &
#ifdef LTN
                                !xb110> output snow flux for lightning
                               ,snow_flx                                    &
#endif
                                 )
#else
      call phys_mps_wsm3_p(gt0(1,1),gq0(1,1),gq0(1,1+ncw),ncld,                  &
     &             ww(1,1),prsl_pa(1,1),                                         &
     &             delz(1,1),dtp,grav,cp,cvap,rd,rv,t0c,                         &
     &             fv,rdorv,qmin,hsub,hvap,hfus,rhoair0,rhoh2o,                  &
     &             cliq,cice,psat,                                               &
     &             lat, rain1,rainncv,snow1(1,1),snowncv(1,1),                   &
     &             sr,                                                           &
     &             1,igrd1s, 1,1, 1,levr,                                        &
     &             1,igrd1s, 1,1, 1,levr,                                        &
     &             1,lonlens, 1,1, 1,levr) 
!      call cloud3(lonlens,igrd1s, levr,dtp,psexp,                        &
!     &            gt0(1,1),gq0(1,1),gq0(1,ncw+1),ncld,sl,delp,slk,       &
!     &            rain1,lat,vvel,kdt,fhour)
#endif
      elseif (nmmiph .eq. 5) then
      call phys_mps_wsm5_p(gt0(1,1),gq0(1,1),gq0(1,1+ncw),ncld,prsl_pa(1,1),     &
     &             delz(1,1),dtp,grav,cp,cvap,rd,rv,t0c,                         &
     &             fv,rdorv,qmin,hsub,hvap,hfus,rhoair0,rhoh2o,                  &
     &             cliq,cice,psat,                                               &
     &             lat, rain1,rainncv,sr,                                        &
     &             1,igrd1s, 1,1, 1,levr,                                        &
     &             1,igrd1s, 1,1, 1,levr,                                        &
     &             1,lonlens, 1,1, 1,levr,                                       &
     &             snow1(1,1), snowncv(1,1))
!     call cloud5(1, 1, levr,dtp,psexp,
!    1            gt0,gq0(1,1),gq0(1,ncw+1),ncld,sl,delp,slk,
!    2            rain1,lat,vvel,kdt,fhour)

      elseif (nmmiph .eq. 6) then
#ifdef NEW_WSM6
          im = lonlens
          ix = igrd1s
          ims = 1 ; ime = ix ; kms = 1 ; kme = levr ;
          its = 1 ; ite = im ; kts = 1 ; kte = levr ;
          ntcw=1+ncw
          ntrw=1+ncw+1*levr
          ntiw=1+ncw+2*levr
          ntsw=1+ncw+3*levr
          ntgl=1+ncw+4*levr

           call wsm6(gt0(1,1), phii(1:ix,1:levr+1),                              &
                                gq0(1:ix,1:levr),                                &
                                gq0(1:ix,ntcw:ntcw+levr-1),                      &
                                gq0(1:ix,ntrw:ntrw+levr-1),                      &
                                gq0(1:ix,ntiw:ntiw+levr-1),                      &
                                gq0(1:ix,ntsw:ntsw+levr-1),                      &
                                gq0(1:ix,ntgl:ntgl+levr-1),                      &
                                ! rain1 will be zero out and calculate in WSM6
                                prsl_pa(1,1), delz, dtp, rain1,                  &
                                sr,                                              &
                                islmsk,                                          &
                                !re_cloud, re_ice, re_snow ! for radiation
                                phy_f3dv(:,:,1),phy_f3dv(:,:,2),phy_f3dv(:,:,3), &
                                ims,ime, kms,kme,                                &
                                its,ite, kts,kte)

#else
      call phys_mps_wsm6_p(gt0(1,1),gq0(1,1),gq0(1,1+ncw),ncld,prsl_pa(1,1),    &
     &             delz(1,1),dtp,grav,cp,cvap,rd,rv,t0c,                        &
     &             fv,rdorv,qmin,hsub,hvap,hfus,rhoair0,rhoh2o,                 &
     &             cliq,cice,psat,                                              &
     &             lat, rain1,rainncv,sr,                                       &
     &             1,igrd1s, 1,1, 1,levr,                                       &
     &             1,igrd1s, 1,1, 1,levr,                                       &
     &             1,lonlens, 1,1, 1,levr,                                      &
     &             snow1(1,1), snowncv(1,1), graupel1(1,1), graupelncv(1,1))
!     call cloud6(1, 1, levr,dtp,psexp,
!    1            gt0,gq0(1,1),gq0(1,ncw+1),ncld,sl,delp,slk,
!    2            rain1,lat,vvel,kdt,fhour)

#endif
      elseif ( nmmiph.eq.11 .or. nmmiph.eq.12 .or. nmmiph.eq.13) then  ! GFDL MP v1 & v2 & v3
          im = lonlens
          ix = lonlens!shian igrd1s

          hydrostatic = .false.        !flag for hydrostatic solver
          phys_hydrostatic = .true.    !flag for hydrostatic heating from physics
          sedi_w = .false.             !flag for w transfer
          effr_in = .true.             !flag for using effective radii of hydrometeors in radiation process

          consv_te = .false.           !flag for energy conservation
          last_step = .true.           !flag for final clean-up (not sure)
          do_inline_mp = .false.       !flag for inline GFDLMP

          allocate                                                      &
     &     ( area(im),rain0(im),ice0(im),snow0(im),graupel0(im),        &
     &       rqv(im,levr),rqc(im,levr),rqr(im,levr),rqi(im,levr),       &
     &       rqs(im,levr),rqg(im,levr),rqa(im,levr),rqn(im,levr),       &
     &       rt0(im,levr),rw0(im,levr),ru0(im,levr),rv0(im,levr),       &
     &       rdp(im,levr),rdz(im,levr),frland(im),gsize(im),hs(im),     &
     &       water0(im) )
          area  = 0.0 ; frland = 0.0 ; gsize = 0.0 ; hs   = 0.0
          rain0 = 0.0 ; snow0 = 0.0 ; ice0 = 0.0  ; graupel0 = 0.0
          water0 = 0.0

        if ( nmmiph .eq. 11 ) then
          allocate                                                      &
     &     ( qv_dt(im,levr),qc_dt(im,levr),qr_dt(im,levr),              &
     &       qi_dt(im,levr),qs_dt(im,levr),qg_dt(im,levr),              &
     &       qa_dt(im,levr),pt_dt(im,levr),udt(im,levr),vdt(im,levr) )
          if ( effr_in ) allocate ( rho(im,levr),rew(im,levr),          &
     &       rei(im,levr),rer(im,levr),res(im,levr),reg(im,levr) )
          qv_dt = 0.0 ; qc_dt = 0.0 ; qr_dt = 0.0 ; qi_dt = 0.0 
          qs_dt = 0.0 ; qg_dt = 0.0 ; qa_dt = 0.0 ; pt_dt = 0.0
          udt   = 0.0 ; vdt   = 0.0

          if ( effr_in ) then
            rew = 10.
            rei = 50.
            rer = 1000.
            res = 250.
            reg = 1000.
          endif

        elseif ( nmmiph.eq.12 .or. nmmiph.eq.13 ) then
          allocate                                                      &
     &     ( qnl(im,levr),qni(im,levr),q_con(im,levr),                  &
     &       cappa(im,levr),te(im,levr) )
#ifdef EXT_DIAG
          allocate                                                      &
     &     ( cond0(im),dep0(im),evap0(im),sub0(im),prefluxr(im,levr),   &
     &       prefluxi(im,levr),prefluxs(im,levr),prefluxg(im,levr),     &
     &       prefluxw(im,levr) )
          cond0 = 0.0 ; dep0 = 0.0 ; evap0 = 0.0 ; sub0 = 0.0
          prefluxr = 0.0 ; prefluxi = 0.0
          prefluxs = 0.0 ; prefluxg = 0.0 ; prefluxw = 0.0
#endif
          qnl   = 0.0 ; qni  = 0.0 ; q_con = 0.0
          cappa = 0.0 ; te   = 0.0 
        endif

        do i = 1, im
            area(i) = rdelx * rdely / xm2(i,lat)  !area of grid box (m^2)
            gsize(i) = sqrt(area(i))              !length of grid box (m)
            hs(i) = phil(i,1)                     !terrain geopotential (gpm)
            if( islmsk(i).eq.1 ) frland(i) = 1.   !land fraction

            do k = 1, levr
              kc = levr - k + 1

              rqv(i,k) =  gq0(i,    kc)
              rqc(i,k) =  gq0(i,ncw+kc)
              rqr(i,k) =  gq0(i,nrw+kc)
              rqi(i,k) =  gq0(i,niw+kc)
              rqs(i,k) =  gq0(i,nsw+kc)
              rqg(i,k) =  gq0(i,ngl+kc)
              rqa(i,k) =  0.  !layer cloud fraction (should set to zero for do_qa=.false.)
              rqn(i,k) =  0.
              rt0(i,k) =  gt0(i,kc)
              ru0(i,k) =  gu0(i,kc)
              rv0(i,k) =  gv0(i,kc)
              rdp(i,k) =  delp_pa(i,kc)  !difference of interface pressure (Pa)
              rdz(i,k) = -delz(i,kc)     !difference of height (m) , rdz<0

              if ( sedi_w ) then
                rw0(i,k) =  ww(i,kc)
              else
                rw0(i,k) =  0.
              endif
            enddo
        enddo

        if ( nmmiph .eq. 11 ) then
          call gfdl_cloud_microphys_driver                              &
     &          ( rqv, rqc, rqr, rqi, rqs, rqg, rqa, rqn,               &
     &            rt0, rw0, ru0, rv0,                                   &
     &            qv_dt, qc_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,      &
     &            pt_dt, udt, vdt,                                      &
     &            rdz, rdp, area, dtp, frland,                          &
     &            rain0, snow0, ice0, graupel0,                         &
     &            hydrostatic, phys_hydrostatic,                        &
     &            1, im, 1, levr, 1, levr )

          do k = 1, levr
            kc = levr - k + 1
            do i = 1, im
              rqv(i,kc) = rqv(i,kc) + qv_dt(i,kc) * dtp
              rqc(i,kc) = rqc(i,kc) + qc_dt(i,kc) * dtp
              rqr(i,kc) = rqr(i,kc) + qr_dt(i,kc) * dtp
              rqi(i,kc) = rqi(i,kc) + qi_dt(i,kc) * dtp
              rqs(i,kc) = rqs(i,kc) + qs_dt(i,kc) * dtp
              rqg(i,kc) = rqg(i,kc) + qg_dt(i,kc) * dtp
              rqa(i,kc) = rqa(i,kc) + qa_dt(i,kc) * dtp
              rt0(i,kc) = rt0(i,kc) + pt_dt(i,kc) * dtp
              ru0(i,kc) = ru0(i,kc) + udt  (i,kc) * dtp
              rv0(i,kc) = rv0(i,kc) + vdt  (i,kc) * dtp
            enddo
          enddo

          if ( effr_in ) then
            do k = 1, levr
              kc = levr - k + 1
              do i = 1, im
                rho(i,kc) = 0.622*prsl_pa(i,k)                          &
     &                     /( rd*rt0(i,kc)*(rqv(i,kc)+0.622) ) !air density (kg/m^3)
              enddo
            enddo

            call cloud_diagnosis                                        &
     &           ( 1, im, 1, levr, rho, rdp, islmsk,                    &
     &             rqc, rqi, rqr, rqs, rqg, rt0,                        &  !module_mp_gfdl_fv3_v16.F
     &             rew, rei, rer, res, reg )

            do k = 1, levr
              do i = 1, im
                phy_f3dv(i,k,1) = rew(i,k)   !(micron) k=1 at top
                phy_f3dv(i,k,2) = rei(i,k)   !(micron)
                phy_f3dv(i,k,3) = rer(i,k)   !(micron)
                phy_f3dv(i,k,4) = res(i,k)   !(micron)
                phy_f3dv(i,k,5) = reg(i,k)   !(micron)
              enddo
            enddo
          endif  !end of if effr_in

          deallocate                                                    &
           ( qv_dt,qc_dt,qr_dt,qi_dt,qs_dt,qg_dt,qa_dt,pt_dt,udt,vdt )
          if ( effr_in ) deallocate ( rho,rew,rei,rer,res,reg )

        elseif ( nmmiph .eq. 12 ) then
          call gfdlv2_driver                                            &
     &          ( rqv, rqc, rqr, rqi, rqs, rqg, rqa, qnl, qni,          &
     &            rt0, rw0, ru0, rv0, rdz, rdp, gsize, dtp, hs,         &
     &            rain0, snow0, ice0, graupel0, hydrostatic,            &
     &            1, im, 1, levr, q_con, cappa, consv_te, te,           &
#ifdef EXT_DIAG
     &            prefluxr, prefluxi, prefluxs, prefluxg,               &
     &            cond0, dep0, evap0, sub0,                             &
#endif
     &            last_step, do_inline_mp, me )

          deallocate                                                    &
     &     ( qnl,qni,q_con,cappa,te )
#ifdef EXT_DIAG
          deallocate                                                    &
           ( prefluxr,prefluxi,prefluxs,prefluxg,prefluxw,              &
     &       cond0,dep0,evap0,sub0 )
#endif

        elseif ( nmmiph .eq. 13 ) then
          call gfdlv3_driver                                            &
     &          ( rqv, rqc, rqr, rqi, rqs, rqg, rqa, qnl, qni,          &
     &            rt0, rw0, ru0, rv0, rdz, rdp, gsize, dtp, hs,         &
     &            frland,                                               &
     &            water0, rain0, snow0, ice0, graupel0, hydrostatic,    &
     &            1, im, 1, levr, q_con, cappa, consv_te, te,           &
#ifdef EXT_DIAG
     &            prefluxw, prefluxr, prefluxi, prefluxs, prefluxg,     &
     &            cond0, dep0, evap0, sub0,                             &
#endif
     &            last_step, do_inline_mp )

          deallocate                                                    &
     &     ( qnl,qni,q_con,cappa,te )
#ifdef EXT_DIAG
          deallocate                                                    &
           ( prefluxr,prefluxi,prefluxs,prefluxg,prefluxw,              &
     &       cond0,dep0,evap0,sub0 )
#endif
        endif

          do k = 1, levr
            kc = levr - k + 1
            do i = 1, im
              gq0(i,    k) = rqv(i,kc)
              gq0(i,ncw+k) = rqc(i,kc)
              gq0(i,nrw+k) = rqr(i,kc)
              gq0(i,niw+k) = rqi(i,kc)
              gq0(i,nsw+k) = rqs(i,kc)
              gq0(i,ngl+k) = rqg(i,kc)
              gt0(i,k)     = rt0(i,kc)
              gu0(i,k)     = ru0(i,kc)
              gv0(i,k)     = rv0(i,kc)
              cldtot(i,k)  = rqa(i,kc)  !k=1 at sfc
            enddo
          enddo

          do i = 1, im
            rain1(i) = water0(i)+rain0(i)+snow0(i)+ice0(i)+graupel0(i)  !total large scale precipitation (m)
          enddo

          deallocate                                                    &
     &     ( area,frland,gsize,hs,rain0,snow0,ice0,graupel0,water0,     &
     &       rqv,rqc,rqr,rqi,rqs,rqg,rqa,rqn,rt0,rw0,ru0,rv0,rdp,rdz )

        elseif ( nmmiph.eq.15 ) then   ! GCE 3ICE CWARSM V1

          im = lonlens
          ix = igrd1s
          ims = 1 ; ime = ix ; kms = 1 ; kme = levr ;
          its = 1 ; ite = im ; kts = 1 ; kte = levr ;

          allocate                                                      &
           ( t3d(ix,levr,1),qq3d(ix,levr,1),qv3d(ix,levr,1),            &
             qc3d(ix,levr,1), th3d(ix,levr,1),                           &
             qr3d(ix,levr,1),qi3d(ix,levr,1),qs3d(ix,levr,1),           &
             qg3d(ix,levr,1),rho3d(ix,levr,1),                          &
             pii3d(ix,levr,1),p3d(ix,levr,1),phii3d(ix,levr+1,1),            &
             ht(ix,1),dz3d(ix,levr,1),rainnc2d(ix,1),                   &
             sr2d(ix,1),snownc2d(ix,1),graupelnc2d(ix,1),               &
             rainncv2d(ix,1),snowncv2d(ix,1),graupelncv2d(ix,1),        &
             refl_10cm(ix,levr,1) )

          qq3d = 0. ; qv3d = 0. ; qc3d = 0. ; qr3d = 0. ; qi3d = 0. ; qs3d = 0.
          qg3d = 0.
          t3d = 0.  ; pii3d = 0.; th3d = 0. 
          rho3d = 0.
          p3d = 0.  ; phii3d = 0. ; dz3d = 0.
          ht = 0.
          rainnc2d = 0. ; snownc2d = 0. ; graupelnc2d = 0.
          rainncv2d = 0. ; snowncv2d = 0. ; graupelncv2d = 0. 
          sr2d = 0.
!          land2d = 0.
          refl_10cm = 0.

          diagflag = .false.          !if diagflag=true and do_radar_ref=1, call calc_refl10cm
          do_radar_ref = 0


          do i = 1, ix
            ht(i,1) = phii(i,1) / grav      !terrain geopotential height above sea level (m)
          enddo
!          print *,' debug1=', levr,ix,  gt0(1,1), prsl_pa(1,1)

          do k = 1, levr
            do i = 1, ix
              qq3d (i,k,1) = gq0(i,    k)
              qv3d (i,k,1) = qq3d (i,k,1) / (1. - qq3d (i,k,1))
              qc3d (i,k,1) = gq0(i,ncw+k) 
              qr3d (i,k,1) = gq0(i,nrw+k) 
              qi3d (i,k,1) = gq0(i,niw+k) 
              qs3d (i,k,1) = gq0(i,nsw+k) 
              qg3d (i,k,1) = gq0(i,ngl+k) 
              p3d  (i,k,1) = prsl_pa(i,k)                          !layer mean pressure (Pa)
              pii3d(i,k,1) = (p3d  (i,k,1) / 100000.)**(rocp)
              t3d  (i,k,1) = gt0(i,k)                              !temperature (K)
              th3d (i,k,1) = gt0(i,k) / ((prsl_pa(i,k) / 100000.)**(rocp))
              dz3d (i,k,1) = delz(i,k)                             !layer thickness (m)
              rho3d(i,k,1) = 0.622*prsl_pa(i,k)                         &
                            /(rd*gt0(i,k)*(qv3d (i,k,1)+0.622))    !air density (kg/m^3)
            enddo
          enddo

          do k = 1, levr+1
            do i = 1, ix
              phii3d(i,k,1) = phii(i,k)                            !geopotential 
            enddo
          enddo

          call gsfcgce_3ice_nuwrf(th3d(1:ix,1:levr,1:1),                     &
                       qv3d(1:ix,1:levr,1:1), qc3d(1:ix,1:levr,1:1),         &
                       qr3d(1:ix,1:levr,1:1), qi3d(1:ix,1:levr,1:1),         &
                       qs3d(1:ix,1:levr,1:1),                                &
                       rho3d(1:ix,1:levr,1:1), pii3d(1:ix,1:levr,1:1),       &
                       p3d(1:ix,1:levr,1:1), dtp,phii3d(1:ix,1:levr+1,1:1),  &
                       ht(1:ix,1:1), dz3d(1:ix,1:levr,1:1), grav,            &
                       ww(1:ix,1:levr), 1000., 100.,                     &
                       1, islmsk,                                            &
!!                       ids,ide, jds,jde, kds,kde,                          & ! domain dims
!!                      ims,ime, 1,1, kms,kme,                               & ! memory dims
!!                       its,ite, 1,1, kts,kte,                              & ! tile   dims
                       1, ix, 1, 1, 1, levr,                                 & ! memory dims
                       1, im, 1, 1, 1, levr,                                 & ! tile   dims                      
                       rainnc2d(1:ix,1:1), rainncv2d(1:ix,1:1),              &
                       snownc2d(1:ix,1:1), snowncv2d(1:ix,1:1),              &
                       sr2d(1:ix,1:1),                                       &
                       graupelnc2d(1:ix,1:1), graupelncv2d(1:ix,1:1),        &
                       refl_10cm(1:ix,1:levr,1:1), diagflag,                 &
                       do_radar_ref,                                         &
!!                       f_qg, qg, pgold,                                     &
                       .true., qg3d(1:ix,1:levr,1:1),                        &
                       0, 0                                                  &
                                                 )


          do k = 1, levr
            do i = 1, ix
              gq0(i,    k) = qv3d (i,k,1) / (1. + qv3d (i,k,1))
              gq0(i,ncw+k) = qc3d (i,k,1) 
              gq0(i,nrw+k) = qr3d (i,k,1) 
              gq0(i,niw+k) = qi3d (i,k,1) 
              gq0(i,nsw+k) = qs3d (i,k,1) 
              gq0(i,ngl+k) = qg3d (i,k,1) 
              gt0(i,    k) = th3d  (i,k,1) *  (p3d  (i,k,1) / 100000.)**(rd/cp)
            enddo
          enddo

          do i = 1, ix
            rain1(i) = rainnc2d(i,1) / 1000.   !total large scale precipitation (m)
            sr   (i) = sr2d(i,1) / 1000.       !solid precipitation (m)
          enddo

          deallocate                                                    &
           ( t3d,th3d,qq3d,qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,rho3d,pii3d,   &
             p3d,phii3d,ht,dz3d,rainnc2d,sr2d,snownc2d,graupelnc2d,     &
             rainncv2d,snowncv2d,graupelncv2d,refl_10cm )


        endif  !end of if nmmiph
      endif   !end of if ncld

!xb110> lightning
#ifdef LTN
      if (ncld .eq. 2 .or. ncld .eq. 5)then    !for wsm3 or wsm6 or GFDLMP
        if (ncld .eq. 5) snow_flx = 0.
          call lightning_ecn(me,lonlens,igrd1s,levr,gt0,gu0,gv0         &
                            ,ncld,levh,ncw,gq0                          &
                            ,prsl_pa,phii,snow_flx,islmsk,flash)
! for output
        if(lssav) then
          do j = 1, lonlens 
            cflash(j,lat)=cflash(j,lat)+flash(j)*dtf
          enddo
        endif
      end if
#endif
!xb110<
#ifdef fix_neg_q
      if ( ncld .eq. 5 ) then
         ix = igrd1s
         im = lonlens
         call fix_neg_q6(ix, im, levr, levh, ncw, nrw, niw, nsw, ngl,   &
                         delp_pa, gt0, gq0, qcmin)
      endif
#endif

!cc
#ifdef A
      do k=1,levr
      do j=1,lonlens
        i=j+(lat-1)*igrd1s
#endif
#ifdef T
        wt(i,k,8)=(gt0(j,k)-wt(i,k,8))/dt2
#endif
#ifdef Q
        wq(i,k,7)=(gq0(j,k)-wq(i,k,7))/dt2
#endif
#ifdef A
      enddo
      enddo
#endif
!
      do 440 j=1,lonlens
      rainl(j)=frain*rain1(j)
      rain(j)=rainc(j)+rainl(j)
      if(lssav) geshem(j,lat)=geshem(j,lat)+rain(j)
  440 continue
#endif
! end of doLRGSCL
!      if ( levr .gt. 100) then
!--------------------------------------------------------------------------
!-----------------------------------------------------------------------
!      endif
!
!  estimate t850 for rain-snow decision
!
      do i=1,lonlens
        t850(i) = gt0(i,1)
      enddo
      do k = 1, levr - 1
        do i=1,lonlens
          if(prsl(i,k) .gt. p850 .and. prsl(i,k+1) .le. p850) then
            t850(i) = gt0(i,k) - (prsl(i,k)-p850)                        &   
     &              / (prsl(i,k)-prsl(i,k+1)) * (gt0(i,k)-gt0(i,k+1))
          endif
        enddo
      enddo
!
#ifdef doPROGT2
!  factor=weighted mean tep.
      do 450 j=1,lonlens
      if(t850(j).le.273.16) then
        if(slmsk(j,lat).ne.0.) then
          sheleg(j,lat)=sheleg(j,lat)+1.e3*rain(j)
          rain(j)=0.
        endif
      endif
  450 continue
!
!
!  update soil moisture and canopy water after precipitation has been
!  computed
!
      call progt2(lonlens,igrd1s,lsoil,rhscnpy,rhsmc,ai,bi,cci,smsoil,   &  
     &  slmsk(1,lat),canopy(1,lat),rain,runof,snowmt,                    &
     &  zsoil,soiltyp,sigmaf,dtf)
!
!  add srunoff soilm snwdph by zyf
!
      if(lssav) then
      do j = 1, lonlens
         soilm(j,lat) = -1.0*smsoil(j,1)*zsoil(j,1)
         do k=2,lsoil
           soilm(j,lat) = soilm(j,lat) +                                 &
     &                    smsoil(j,k)*(zsoil(j,k-1)-zsoil(j,k))
         enddo
         if (slmsk(j,lat).eq.1) then
           snwdph(j,lat) = sheleg(j,lat) / 100.0
         else if (slmsk(j,lat).eq.2) then
           snwdph(j,lat) = sheleg(j,lat) / 1000.0 * 3.0
         endif
      enddo
      endif
#endif
! end of doPROGT2
!
!  total runoff is composed of drainage into water table and
!  runoff at the surface and is accumulated in unit of meters
!
      if(lssav) then
      do j = 1, lonlens
        runoff(j,lat) = (drain(j) + runof(j)) * dtf / 1000.              &
     &                + runoff(j,lat)
       srunoff(j,lat) = runof(j) * dtf / 1000. + srunoff(j,lat)
      enddo
      endif

!
!  return updated smsoil and stsoil to global arrays
!
      do k = 1, lsoil
        do j = 1, lonlens
          smc(j,lat,k) = smsoil(j,k)
          stc(j,lat,k) = stsoil(j,k)
!yj:      slc(j,lat,k) = slsoil(i,k)
        enddo
      enddo
!
!ccj
!
! calc. integral of moisture in pwat
!
      do k=1,levr
        do i=1,lonlens
          work3(i,k) = 0.0
        enddo
      enddo
!
      if (ncld .gt. 0) then
       ncw = (1+ntoz)*levr
       do ic=ncw,levr*(ntotal-1),levr
         do k=1,levr
           do j=1,lonlens
            work3(j,k) = work3(j,k) + gq0(j,ic+k)
           enddo
         enddo
       enddo
      endif
!ccc
!
#ifdef NONHYD
      cpocv=cp/ccv
#endif
      do 493 j=1,lonlens
      pwat(j,lat)=0.
      do 492 k=1,levr
! restore all drt temperatures back to virtual temperatures
      tgrs(j,k)=tgrs(j,k)*(1.+fv*max(qgrs(j,k),qmin))
#ifdef T
      i=j+(lat-1)*igrd1s
      wt(i,k,11)=-wt(i,k,11)+gt0(j,k)
#endif
      gt0(j,k)=gt0(j,k)*(1.+fv*max(gq0(j,k),qmin))
#ifdef T
      wt(i,k,11)=(gt0(j,k)-wt(i,k,11))/dt2
#endif
#ifdef NONHYD
#ifdef HYDPHY
      tgrh(j,k)=tgrh(j,k)*(1.+fv*max(qgrs(j,k),qmin))
      dtv=gt0(j,k)-tgrs(j,k)
      gh0(j,k)=tgrh(j,k)+dtv
!hmhj isobaric process even for NONHYD
      gp0(j,k)=pgrs(j,k)
#else
      dtv=gt0(j,k)-tgrs(j,k)
      gp0(j,k)=pgrs(j,k)+(cpocv*dtv/tgrs(j,k))
#endif
#endif
      pwat(j,lat)=pwat(j,lat)+delp(j,k)*(gq0(j,k)+work3(j,k))
      wvuflx(j,lat)=wvuflx(j,lat)                                        &
     &              +delp(j,k)*(gq0(j,k)+work3(j,k))*gu0(j,k)            &
     &               *(1.e3/grav)*dtf
      wvvflx(j,lat)=wvvflx(j,lat)                                        &
     &              +delp(j,k)*(gq0(j,k)+work3(j,k))*gv0(j,k)            &
     &               *(1.e3/grav)*dtf
  492 continue
      pwat(j,lat)=pwat(j,lat)*(1.e3/grav)
  493 continue
!
!  compute total mass and water mass over the domain
!x    ptotj=0.
!x    pwatj=0.
!x    do 903 j=1,lonlens
!x    ptotj=ptotj+psexp(j)
!x    pwatj=pwatj+pwat(j,lat)*(grav/1.e3)
!x903 continue
!x    ptotj=ptotj/lonlens
!x    pwatj=pwatj/lonlens
!
!      print *,'end of rbphys,me=',me
      return
      end
