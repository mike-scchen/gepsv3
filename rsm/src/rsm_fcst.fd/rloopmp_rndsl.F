#include <define.h>
       subroutine rloopmp_rndsl
#ifdef MP
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rloopa
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  
!   program  starts with spectral coefficients temp.
!   of x-wind, y-wind, specific humidity, and
!   ln((psfc).  converts them to the gaussian grid at each
!   row and calls rfidi,  for the northern and southern
!   rows at the same time.  after return from rfidi
!   sr.  completes calculation of perturbation tendencies 
!   of temp. wind. and lnps. specific humidity with relaxation. 
!   all input/output  is via commons.
!
! program history log:
!
!   2018Jul05: yj - add ndsl (do physics)
! usage:    call rloopa
!   input argument list:
!
!   output argument list:
!
!   common block:
!       rscomf  rscommap  rscomltb
!
!   input files: none
!
!   output files: none
!
!   subpcrograms called: 
!       dellcc  dellsc  dellcs  sumfcs  sumffc  ffsncs ffscos  pint3
!       rfidi   pint1   ffancs  ffacos  sumgcs  sumggc 
!       
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
!.....
!.....fold north y and south y 
!.....use perturbation method by substraction of full field forcing
!     by large scale forcing.
!.....
!ccj
        use physcons, pi => con_pi, grav => con_g
        use physics_module
        use radiation_module
        use cnttime
!c
#include <paramodel.h>
#include <ncpu.h>
#include <rscomcpu.h>
#include <npesi.h>
#include <comrmpi.h>
!
#include <rscomfi.h>
#include <rscomspc.h>
#include <rscomgrd.h>
#include <rscomver.h>
#include <rscomio.h>
#include <rscompln.h>
#include <rscomrad.h>
#include <rscomloc.h>
!
#include <rscommap.h>
#include <rscomltb.h>
#include <rscomgpd.h>
#include <rscombgt.h>
#include <rscomoz.h>
#include <rdparm.h>
!
!ifdef AFA
!     parameter(ncpus = ncpu )
!endif
#ifdef T
      dimension at(igrd1p,levr),bt(igrd1p,levr),ct(igrd1p,levr)
#endif
#ifdef Q
      dimension aq(igrd1p,levr,ntotal),bq(igrd1p,levr,ntotal)
#endif
#ifdef U
      dimension au(igrd1p,levr),bu(igrd1p,levr),fu(igrd1p,levr)
     &         ,pu(igrd1p,levr)
#endif
#ifdef V
      dimension av(igrd1p,levr),bv(igrd1p,levr),fv(igrd1p,levr)
     &         ,pv(igrd1p,levr)
#endif
#ifdef P
      dimension ap(igrd1p),bp(igrd1p)
#endif
      real, allocatable,dimension(:,:,:) :: dot
!      dimension dot(igrd1p,levr+1,jgrd1p)
!....
      parameter(numcs = 1 + levr*5       , numcs2 = numcs * 2 )
      parameter(numcos = 1 + levr*1       , numcos2 = numcos * 2 )
!                 cos             sin          in x direction
      parameter( kpy=1       , kpx=numcs+kpy,                            &  
     &           kvv=1   +kpy, kuu=numcs+kvv,                            &
     &           kvm=levr+kvv, kum=numcs+kvm,                            &
     &           kty=levr+kvm, ktx=numcs+kty,                            &
     &           kvy=levr+kty, kvx=numcs+kvy,                            &
     &           kux=levr+kvy, kuy=numcs+kux )
      parameter( kpp=numcs2+1,                                           &
     &           ktt=1    +kpp,                                          &
     &           ktm=levr+ktt,                                           &
     &           kpm=levr+ktm)
!
!for mpi
      parameter(numcsp =1+levrp*5      , numcs2p=numcsp*2)
      parameter(numcosp=1+levrp*1      , numcos2p=numcosp*2)
      parameter( kpyp=1         , kpxp=numcsp+kpyp,                      &
     &           kvvp=1    +kpyp, kuup=numcsp+kvvp,                      &
     &           kvmp=levrp+kvvp, kump=numcsp+kvmp,                      &
     &           ktyp=levrp+kvmp, ktxp=numcsp+ktyp,                      &
     &           kvyp=levrp+ktyp, kvxp=numcsp+kvyp,                      &
     &           kuxp=levrp+kvyp, kuyp=numcsp+kuxp )
      parameter( kppp=numcs2p+1,                                         &
     &           kttp=1     +kppp,                                       &
     &           ktmp=levrp+kttp,                                        &
     &           kpmp=levrp+ktmp)
!
      parameter( numbgd=1+levr*3+levh,                                   &
     &           kplb=1,                                                 &
     &           kulb=kplb+1,                                            &
     &           kvlb=kulb+levr,                                         &
     &           ktlb=kvlb+levr,                                         &
     &           kqlb=ktlb+levr )
      parameter( kppb=numbgd+1,                                          &
     &           kuub=kppb+1,                                            &
     &           kvvb=kuub+levr,                                         &
     &           kttb=kvvb+levr,                                         &
     &           kqqb=kttb+levr)
      parameter( kpxb=2*numbgd+1,                                        &
     &           kuxb=kpxb+1,                                            &
     &           kvxb=kuxb+levr,                                         &
     &           ktxb=kvxb+levr,                                         &
     &           kqxb=ktxb+levr )
      parameter( kpyb=3*numbgd+1,                                        &
     &           kuyb=kpyb+1,                                            &
     &           kvyb=kuyb+levr,                                         &
     &           ktyb=kvyb+levr,                                         &
     &           kqyb=ktyb+levr )
!!
#ifdef NONHYD
#ifdef NDSLDYN
      parameter( noncs=  levr*1,noncs2=noncs*2,                          &  
#else
      parameter( noncs=1+levr*3,noncs2=noncs*2,                          &  
#endif
     &           locncs=1,                                               &
     &           knpy=locncs   ,   knpx=noncs+knpy,                      &
#ifndef NDSLDYN
     &           knty=knpy+levr,   kntx=noncs+knty,                      &
     &           knwy=knty+levr,   knwx=noncs+knwy,                      &
#endif
     &           noncos=1+levr*3,noncos2=noncos*2,                       &
     &           locncos=locncs+noncs2,                                  &
     &           knpp=locncos   ,  knpm=noncos+knpp,                    &
     &           kntt=knpp+levr,   kntm=noncos+kntt,                     &
     &           knww=kntt+levr,   knwm=noncos+knww)
      parameter( knppa=1,                                                &
     &           kntta=knppa+levr,                                       &
     &           knwwa=kntta+levr )
      parameter( knppm=1,                                                &
     &           knttm=knppm+levr,                                       &
     &           knwwm=knttm+levr )
!for mpi
#ifdef NDSLDYN
      parameter( noncsp=  levrp*1,noncs2p=noncsp*2,                      &
#else
      parameter( noncsp=1+levrp*3,noncs2p=noncsp*2,                      &
#endif
     &           locncsp=1,                                              &
     &           knpyp=locncsp   ,    knpxp=noncsp+knpyp,                &
#ifndef NDSLDYN
     &           kntyp=knpyp+levrp,   kntxp=noncsp+kntyp,                &
     &           knwyp=kntyp+levrp,   knwxp=noncsp+knwyp,                &
#endif
     &           noncosp=1+levrp*3,noncos2p=noncosp*2,                   &
     &           locncosp=locncsp+noncs2p,                               &
     &           knppp=locncosp   ,   knpmp=noncosp+knppp,               &
     &           knttp=knppp+levrp,   kntmp=noncosp+knttp,               &
     &           knwwp=knttp+levrp,   knwmp=noncosp+knwwp)
      parameter(lotsnnk=noncs2p+noncos2p)
      parameter(lotannk=levrp*3+1)
      real, allocatable,dimension(:,:,:) :: snnxk,snnyk,annxk,annyk
!      dimension snnxk(igrd1p,lotsnnk,jgrd1)
!      dimension snnyk(igrd1,lotsnnk,jgrd1p)
!      dimension annxk(igrd1p,lotannk,jgrd1)
!      dimension annyk(igrd1, lotannk,jgrd1p)

!<<
      parameter(lotsnn=noncs2+noncos2)
      parameter(lotann=levr*3+1)
      real, allocatable,dimension(:,:,:) :: snn,snnm,ann
!      dimension snn(igrd1p,lotsnn,jgrd1p)
!      dimension snnm(igrd1p,lotann,jgrd1p)
!      dimension ann(igrd1p,lotann,jgrd1p)
#endif

!jun
      parameter(kus=1,                                                   &
     &          kvs=levr+kus,                                            &
     &          kts=levr+kvs,                                            &
     &          kqs=levr+kts)
      parameter(kpyr=1,        kvvr=1    +kpyr,                          & !cos
     &          kpxr=levr+kvvr,kuur=1    +kpxr,                          & !sin
     &          kddr=levr+kuur,                                          &
     &          kppr=levr+kddr,                                          &
     &          kttr=1   +kppr,                                          &
     &          kqqr=levr+kttr )
!synm>      
      parameter(kuum=1,                                                  &
     &          kvvm=levr+kuum,                                          &
     &          kttm=levr+kvvm,                                          &
     &          kppm=levr+kttm,                                          & !3D ln(ps)
     &          kqqm=levr+kppm                                           & !tracers
#ifdef NONHYD     
#ifdef NDSLDYN
     &         ,kgzm=levh+kqqm                                           &!geopotential height
#endif
#endif
     &          )
!synm<      
#ifdef NDSLDYN
!NDSL:synmk>
      parameter(kumk=1,                                                  &
     &          kvmk=levrp+kumk,                                         &
     &          ktmk=levrp+kvmk,                                         &
     &          kpmk=levrp+ktmk,                                         & !3D ln(ps)
     &          kqmk=levrp+kpmk                                          & !tracers
#ifdef NONHYD     
     &         ,kgzmk=levhp+kqmk                                         &!geopotential height
     &         ,kpnmk=levrp+kgzmk                                        &
     &         ,ktnmk=levrp+kpnmk                                        &
     &         ,kwmk=levrp+ktnmk                                         &!levrp (start from k=2)
#endif
     &          )
#else
      parameter(kqmk=1)
#endif
!NDSL:synmk<
#ifdef NDSLDYN
!NDSL:synt>
      parameter(kuut=1,                                                 &
     &          kvvt=levr+kuut,                                         &
     &          kttt=levr+kvvt,                                         &
     &          kppt=levr+kttt,                                         &
     &          kqqt=levr+kppt                                          &!kppt to 3D
#ifdef NONHYD     
     &         ,kgzt=levh+kqqt                                          &
     &         ,kpnt=levr+kgzt                                          &
     &         ,ktnt=levr+kpnt                                          &
     &         ,kwwt=levr+ktnt                                          &!levr (start from k=2)
#endif
     &          )
!NDSL:synt<
#else
      parameter(kqqt=1)
#endif

      parameter(                                                         &
     &          kusp=1,                                                  &
     &          kvsp=levrp+kusp,                                         &
     &          ktsp=levrp+kvsp,                                         &
     &          kqsp=levrp+ktsp)
!c
      parameter(lotsyn=numcs2+numcos2)
      parameter(lotsynr=3+4*levr+levh)
#ifdef NONHYD
#ifdef NDSLDYN
      parameter(lotsynm=5*levr+levh)
#else
      parameter(lotsynm=4*levr+levh)
#endif
#else
      parameter(lotsynm=4*levr+levh)
#endif
!NDSL:synmk  >
#ifdef NDSLDYN
#ifdef NONHYD
      parameter(lotsynmk=8*levrp+levhp)
      parameter(nlotsynmk=8+ntotal)
#else
      parameter(lotsynmk=4*levrp+levhp)
      parameter(nlotsynmk=4+ntotal)
#endif
#else
      parameter(lotsynmk=levhp)
      parameter(nlotsynmk=ntotal)
#endif
!NDSL:synmk  <
!NDSL:synt  >
#ifdef NDSLDYN
#ifdef NONHYD
      parameter(lotsynt=8*levr+levh)
#else
      parameter(lotsynt=  4*levr+levh)
#endif
#else
      parameter(lotsynt= levh)
      real, allocatable,dimension(:,:,:) ::utmp,vtmp
!      dimension uvtmp(igrd1p,2*levr,jgrd1p)
#endif
!NDSL:synt  <
      real, allocatable,dimension(:,:,:) :: syn, srq, synr, synm, syndd, bgd
      real, allocatable,dimension(:,:,:) :: syfk, synk                                 

!      dimension syn(igrd1p,lotsyn,jgrd1p)
!      dimension srq(igrd1p,levh,jgrd1p)
!      dimension synr(igrd1p,lotsynr,jgrd1p)
!      dimension synm(igrd1p,lotsynm,jgrd1p)
!      dimension syndd(igrd1p,levr,jgrd1p)
!      dimension bgd(igrd1p,numbgd*4,jgrd1p)
      parameter(lotsynk=numcs2p+numcos2p)
!      dimension syfk(igrd1p,lotsynk,jgrd1)
!      dimension synk(igrd1,lotsynk,jgrd1p)
!NDSL: n-1 fields and n wind      
#ifdef NDSLDYN
#ifdef NONHYD
      real, allocatable,dimension(:,:,:) :: w2nd                                       
!      dimension w2nd(igrd1p,levr,jgrd1p)
#endif
#endif
      real, allocatable,dimension(:,:,:) :: synmk, rindex, synt, pindex                
      real, allocatable,dimension(:,:) :: ppindex                
!      dimension synmk(igrd1,lotsynmk,jgrd1p)
!      dimension rindex(igrd1,levrp,jgrd1p)
!      dimension synt(igrd1p,lotsynt,jgrd1p)
!      dimension pindex(igrd1p,levr,jgrd1p)
      real, allocatable,dimension(:,:,:) :: synuk, synvk        
!      dimension synuk(igrd1,levrp,jgrd1p)
!      dimension synvk(igrd1,levrp,jgrd1p)
      logical llind
      character*4  cfhour
#ifdef RAYLEIFR      
! rayleigh friction      
      integer levtop
      real    windmax
      data levtop/4/
      data windmax/45./
#endif
!
      common/reproa/                                                     &  
     &  tex(lnwavp,levr)                                                 &! cos in y direction
     &, rqx(lnwavp,levr,ntotal)                                          &
     &, vvy(lnwavp,levr)                                                 &
     &, uux(lnwavp,levr)                                                 &
     &, tey(lnwavp,levr)                                                 &! sin in y direction
     &, rqy(lnwavp,levr,ntotal)                                          &
     &, vvx(lnwavp,levr)                                                 &
     &, uuy(lnwavp,levr)                                                 &
     &, texk(llwavp,levrp)                                               &! cos in y direction
     &, rqxk(llwavp,levrp,ntotal)                                        &
     &, vvyk(llwavp,levrp)                                               &
     &, uuxk(llwavp,levrp)                                               &
     &, teyk(llwavp,levrp)                                               &! sin in y direction
     &, rqyk(llwavp,levrp,ntotal)                                        &
     &, vvxk(llwavp,levrp)                                               &
     &, uuyk(llwavp,levrp)
!
#ifdef NONHYD
      common/repronon/                                                   &
     &  pnx(lnwavp,levr)                                                 &
     &, tnx(lnwavp,levr)                                                 &
     &, onx(lnwavp,levrp1)                                               &
     &, pny(lnwavp,levr)                                                 &
     &, tny(lnwavp,levr)                                                 &
     &, ony(lnwavp,levrp1)                                               &!mpi
     &, pnxk(llwavp,levrp)                                               &
     &, tnxk(llwavp,levrp)                                               &
     &, onxk(llwavp,levrp1p)                                             &
     &, pnyk(llwavp,levrp)                                               &
     &, tnyk(llwavp,levrp)                                               &
     &, onyk(llwavp,levrp1p)                                                   
#endif

      parameter(lotanl=levr*3+levh)
      parameter(lotbgd=levr*7+levh+3)
      real, allocatable,dimension(:,:,:) :: anl
!      dimension anl(igrd1p,lotanl, jgrd1p)
!
#ifdef RKN
       dimension                                                         &
     & udx(igrd1p,levr,jgrd1p),udy(igrd1p,levr,jgrd1p),                  &
     & vdx(igrd1p,levr,jgrd1p),vdy(igrd1p,levr,jgrd1p)
!radiation_run
      real (kind=kind_rad) flwupf(igrd1p,jgrd1p),                        &
     &     fswupf(igrd1p,jgrd1p),fswdnf(igrd1p,jgrd1p),                  &
     &     sswupf(igrd1p,jgrd1p),sswdnf(igrd1p,jgrd1p),                  &
     &     slwupf(igrd1p,jgrd1p),slwdnf(igrd1p,jgrd1p),                  &
     &     cldsaf(igrd1p,5,jgrd1p)
      integer mtopaf(igrd1p,3,jgrd1p),mbotaf(igrd1p,3,jgrd1p)
!physics_run
       dimension psexp(igrd1p,jgrd1p),radsl(igrd1p,jgrd1p),              &
     &  dlwsf1(igrd1p,jgrd1p),qss(igrd1p,jgrd1p),                        &
     &  gflx(igrd1p,jgrd1p),cd(igrd1p,jgrd1p),                           &
     &  cdq(igrd1p,jgrd1p),evap(igrd1p,jgrd1p),                          &
     &  hflx(igrd1p,jgrd1p),drain(igrd1p,jgrd1p),                        &
     &  runof(igrd1p,jgrd1p),cld1d(igrd1p,jgrd1p),                       &
     &  gamt(igrd1p,jgrd1p),gamq(igrd1p,jgrd1p),                         &
     &  dqsfc1(igrd1p,jgrd1p),dtsfc1(igrd1p,jgrd1p),                     &
     &  dusfc1(igrd1p,jgrd1p),dvsfc1(igrd1p,jgrd1p),                     &
     &  dusfcg(igrd1p,jgrd1p),dvsfcg(igrd1p,jgrd1p),                     &
     &  rainc(igrd1p,jgrd1p),rainl(igrd1p,jgrd1p),                       &
     &  hsw(igrd1p,levr,jgrd1p),vvel(igrd1p,levr,jgrd1p),                &
     &  stsoil(igrd1p,lsoil,jgrd1p),smsoil(igrd1p,lsoil,jgrd1p),         &
     &  gu0(igrd1p,levr,jgrd1p),gv0(igrd1p,levr,jgrd1p),                 &
     &  gt0(igrd1p,levr,jgrd1p),gq0(igrd1p,levh,jgrd1p)
#endif
!...........................................................
!...
      common/rcomdhc/spdlat(levr,jgrd1p),sm(lnwavp)
!jun
!yj       dimension tmp(lnwav),tmp1(igrd1p,jgrd1p),                      &    
!yj     &   tmpg(lngrd),tmpsk(lnwav,levr+1)
!yj       real(kind=4) atmp(lnwav),atmp1(igrd1p,jgrd1p),                 &
!yj     &   atmpg(lngrd),atmpsk(lnwav,levr)
!yj      dimension tmpu(IGRD1S,jgrd1p),tmpv(IGRD1S,jgrd1p)
!yj      dimension tmpw(IGRD1S,jgrd1p),tmpwork(igrd1,jgrd1)
!jun- radiation
       integer itimsw,itimlw, kk
! zyf test
!yj      dimension zyftmp(lnwavp,levr),zyfgrid(lngrdp,levr)
!yj      dimension zyffull(lngrd,levr)

!yj:add index
!      real(kind=8), allocatable :: indfull(:,:),tmpgrd(:,:)
!      real(kind=8), allocatable :: indfullv(:,:),tmpgrdv(:,:)
!      real(kind=4), allocatable :: inds(:)

!      allocate(indfull(lngrd,levr),inds(lngrd),tmpgrd(lngrdp,levr))
!      allocate(indfullv(lngrd,levr),tmpgrdv(lngrdp,levr))


! yj: allocate local array >-------------------------------------------->
      allocate(dot(igrd1p,levr+1,jgrd1p))
      allocate(syn(igrd1p,lotsyn,jgrd1p), srq(igrd1p,levh,jgrd1p),       &
     &         synr(igrd1p,lotsynr,jgrd1p), synm(igrd1p,lotsynm,jgrd1p), &
     &         syndd(igrd1p,levr,jgrd1p), bgd(igrd1p,numbgd*4,jgrd1p),   &
     &         syfk(igrd1p,lotsynk,jgrd1), synk(igrd1,lotsynk,jgrd1p)  )
      allocate(synmk(igrd1,lotsynmk,jgrd1p), rindex(igrd1,levrp,jgrd1p), &
     &         synt(igrd1p,lotsynt,jgrd1p), pindex(igrd1p,levr,jgrd1p) )
      allocate(ppindex(igrd1p,jgrd1p))
      allocate(synuk(igrd1,levrp,jgrd1p), synvk(igrd1,levrp,jgrd1p),     &
     &         anl(igrd1p,lotanl, jgrd1p) )
#ifdef NONHYD
      allocate(snnxk(igrd1p,lotsnnk,jgrd1), snnyk(igrd1,lotsnnk,jgrd1p), &
     &         annxk(igrd1p,lotannk,jgrd1), annyk(igrd1, lotannk,jgrd1p)) 
      allocate( snn(igrd1p,lotsnn,jgrd1p), snnm(igrd1p,lotann,jgrd1p),   &
     &          ann(igrd1p,lotann,jgrd1p) ) 
#endif
#ifdef NDSLDYN
#ifdef NONHYD
      allocate(w2nd(igrd1p,levr,jgrd1p))
#endif
#endif
#ifndef NDSLDYN 
      allocate( utmp(igrd1p,levr,jgrd1p),vtmp(igrd1p,levr,jgrd1p))
#endif 
! yj: allocate local array <--------------------------------------------<


! yj: set initial values >---------------------------------------------->
#ifdef DOINDEX
      llind=.true.
#else
      llind=.false.
#endif
      syn=0.0
      srq=0.0
      synr=0.0
      synm=0.0
      syndd=0.0
      bgd=0.0
      syfk=0.0
      synk=0.0
      synmk=0.0
      rindex=1.0
      synt=0.0
      pindex=1.0
      ppindex=1.0
      synuk=0.0
      synvk=0.0

#ifdef NDSLDYN
#ifdef NONHYD
      w2nd=0.0
#endif
#endif
! yj: set initial values <----------------------------------------------<
!.................................................................
!
      me=mype
      lonlens=lonlen(mype)
      latlens=latlen(mype)
      lonstrs=lonstr(mype)
      lat1s=1
      lat2s=latlen(mype)
      latts=latstr(mype)-1
!c
!yj  
      dot=0.0
!     if(me.eq.0) print *,' rloopmp start'
!      if(me.eq.0) &
!     &    print *,'rmgrd,rqgrd,rtgrd(2,2,2)=',rmgrd(2,2,2),rqgrd(2,2,2),rtgrd(2,2,2)
!............................................................
! ---------------- q   dqdx    dqdy   d dx  d dy   ---------
      call dellcc(lnwavp,q,dpdlam,dpdphi,epsx,epsy,1)
      call dellsc(lnwavp,uu,uux,uuy,epsx,epsy,levr)
      call dellcs(lnwavp,vv,vvx,vvy,epsx,epsy,levr)
      call dellcc(lnwavp,te,tex,tey,epsx,epsy,levr)
#ifdef NONHYD
#ifdef NDSLDYN
      call dellcc(lnwavp,pn,pnx,pny,epsx,epsy,levr*1)
#else
      call dellcc(lnwavp,pn,pnx,pny,epsx,epsy,levr*3+1)
#endif
#endif
!              tmp=0.
!      call rmpsp2f(on,lnwavp,tmpsk,lnwav,levr)
!      if (mype.eq.master) then
!        call maxmin(tmpsk,lnwav,levr,1,levr,'rloopmp wn sp')
!      endif

!
       if (lsphy) then
       call physics_init(me,lonlens,latlens,lngrdp,                      &  
     &   nmmiph,num_p3d,num_p2d,levr,levh,levrp1,jgrd1p,                 &
     &   igrd1p,jgrd1p,igrd1,jgrd1,fhour,fhswr,phour,lsoil,lsfwd,        &
     &   deltim,lscca,lsswr,                                             &
     &   si(1),idate(1),phy_f3d(1,1,1,1),phy_f2d(1,1,1),                 &
     &   flat,flon)
!
       call radiation_init(me,lonlens,lonstrs,latlens,levr,levrp1,       &
     &   levh,jgrd1p,lat1s,lat2s,latts,l,lp1,lp1v,                       &
     &   igrd1p,jgrd1p,ncpu,ntotal,                                      &
     &   lnwavp,ncld,nvect,num_p3d,nmtvr,                                &
     &   fhour,phour,fhswr,fhlwr,deltim,                                 &
     &   lsswr,lslwr,lssav,lw_type,ntoz,                                 &
     &   si,sl,idate,rlat,rlon,                                          &
     &   sinlar,coslar,itimsw,itimlw)
       endif

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
       if (lsphy) then
        call physics_prep(ozprd(1,1,1),ozdis(1,1,1),                     &
     &    ozprdout(1,1),ozdisout(1,1),solhr,clstp)
!        print *,'in rloopmp,ozprd=',ozprd(1,1,1),ozprd(1,1,lev46)
!
        if(itimsw.ne.0 .or. itimlw.ne.0) then
          call radiation_prep(sdec,cdec,slag,q,uu,te,vv,del,rq,          &
     &      tsea,sheleg,coszer,hprime(1,1,1),stc(1,1,1),                 &
     &      sinlar,coslar)
        endif
       endif
!ccjune
!
      uua=0.
      vva=0.
      tea=0.
      rqa=0.
      xa=0.
      ya=0.
      wa=0.
      rta=0.
      za=0.
#ifdef NONHYD
      pa=0.
      ta=0.
      oa=0.
      pna=0.
      tna=0.
      ona=0.
#endif
!
!           lat loop
! compute lat band limits
!
! first lat loop
!
#ifdef CNT_1
        call cnttime_on("one_loop_spectal")
#endif
!
       call rmpnn2n(q,lnwavp,qa,llwavp,1)
       call rmpnn2n(dpdlam,lnwavp,dpdlama,llwavp,1)
       call rmpnn2n(dpdphi,lnwavp,dpdphia,llwavp,1)
       call rmpnn2n(qm,lnwavp,qma,llwavp,1)
       call rmpnn2nk(uu,lnwavp,levr,uua,llwavp,levrp,2)
       call rmpnn2nk(vv,lnwavp,levr,vva,llwavp,levrp,2)
       call rmpnn2nk(te,lnwavp,levr,tea,llwavp,levrp,1)
       call rmpnn2nk(tem,lnwavp,levr,tema,llwavp,levrp,1)
       call rmpnn2nk(tex,lnwavp,levr,texk,llwavp,levrp,1)
       call rmpnn2nk(vvy,lnwavp,levr,vvyk,llwavp,levrp,3)
       call rmpnn2nk(vvx,lnwavp,levr,vvxk,llwavp,levrp,2)
#ifdef NONHYD
#ifdef NDSLDYN
       call rmpnn2nk(pnx,lnwavp,levr,pnxk,llwavp,levrp,1)
       call rmpnn2nk(pny,lnwavp,levr,pnyk,llwavp,levrp,1)
#else
       call rmpnn2nk(pnx,lnwavp,levr,pnxk,llwavp,levrp,3)
       call rmpnn2n (onx(1,levrp1),lnwavp,onxk(1,levrp1p),llwavp,1)
       call rmpnn2nk(pny,lnwavp,levr,pnyk,llwavp,levrp,3)
       call rmpnn2n (ony(1,levrp1),lnwavp,onyk(1,levrp1p),llwavp,1)
#endif
       call rmpnn2nk(pn,lnwavp,levr,pna,llwavp,levrp,3)
       call rmpnn2n (on(1,levrp1),lnwavp,ona(1,levrp1p),llwavp,1)
       call rmpnn2nk(pnm,lnwavp,levr,pnma,llwavp,levrp,3)
       call rmpnn2n (onm(1,levrp1),lnwavp,onma(1,levrp1p),llwavp,1)
       call rmpnn2nk(p,lnwavp,levr,pa,llwavp,levrp,3)
       call rmpnn2n(o(1,levrp1),lnwavp,oa(1,levrp1p),llwavp,1)
#endif
!
!jun -- define variable
       syfk=0.
! 1.dpdla,dpdphi...uu,vv...uum...vvm     
       call ffsncsy1(dpdlama,dpdphia,syfk,                               &   
     &    kpxp,kpyp,iwav1p,igrd1p,lotsynk,levrp*2+1)
!.2.tex,tey...       ...vvy,vvx...uux,uuy
       call ffsncsy1(texk,teyk,syfk,                                     &
     &    ktxp,ktyp,iwav1p,igrd1p,lotsynk,levrp*1)
       call ffsncsy1(vvyk,vvxk,syfk,                                     &
     &    kvxp,kvyp,iwav1p,igrd1p,lotsynk,levrp*2)
!
#ifdef NONHYD
!.2.NHYD.  pnx,pny...tnx,tny...onx,ony
#ifdef NDSLDYN
      call ffsncsy1(pnxk,pnyk,snnxk,                                     &
     &    knpxp,knpyp,iwav1p,igrd1p,lotsnnk,levrp)
#else
      call ffsncsy1(pnxk,pnyk,snnxk,                                     &
     &    knpxp,knpyp,iwav1p,igrd1p,lotsnnk,levrp*3+1)
#endif
!.2.NHYD.  pna,tna,ona,pnma,tnma,onma
      call ffscosy1(pna,snnxk,knppp,                                      &
     &    iwav1p,igrd1p,lotsnnk,levrp*6+2)
!      
#endif

!.3. q,te,rq,rwm,qm,rm
       call ffscosy1(qa,syfk,kppp,                                       &
     &    iwav1p,igrd1p,lotsynk,levrp+1)
       call ffscosy1(tema,syfk,ktmp,                                       &
     &    iwav1p,igrd1p,lotsynk,levrp+1)

! switch uux uuy and vvy vvx for ffsncs
#ifdef OPENMP
!$omp parallel do private(j,k,i,synhold)
#endif
      do j=1,jgrd1
      do k=1,levrp
      do i=1,lonlens
        synhold            = syfk(i,kuxp-1+k,j)
        syfk(i,kuxp-1+k,j) = syfk(i,kuyp-1+k,j)
        syfk(i,kuyp-1+k,j) = synhold
        synhold            = syfk(i,kvxp-1+k,j)
        syfk(i,kvxp-1+k,j) = syfk(i,kvyp-1+k,j)
        syfk(i,kvyp-1+k,j) = synhold
      enddo
      enddo
      enddo
!     
! second loop
!
      call rmpnl2ny(syfk,igrd1p,jgrd1,synk,igrd1,jgrd1p,                 &   
     &               lotsynk,1,lotsynk)
#ifdef NONHYD
      call rmpnl2ny(snnxk,igrd1p,jgrd1,snnyk,igrd1,jgrd1p,               &
     &               lotsnnk,1,lotsnnk)
#endif
!
      call ffsncsx1(synk,jgrd1p,kpyp,kpxp,lotsynk,numcsp)
      call ffscosx1(synk,jgrd1p,numcs2p+1,lotsynk,numcos2p)
#ifdef NONHYD
      call ffsncsx1(snnyk,jgrd1p,knpyp,knpxp,lotsnnk,noncsp)
      call ffscosx1(snnyk,jgrd1p,locncosp,lotsnnk,noncos2p)
#endif
!
      call rmpx2nx(synk,igrd1,lotsynk,syn,igrd1p,lotsyn,                 & 
     &             jgrd1p,kpyp,kpy,1)
      call rmpnk2nx(synk,igrd1,lotsynk,syn,igrd1p,lotsyn,                &
     &             jgrd1p,levrp,levr,kvvp,kvv,5)
      call rmpx2nx(synk,igrd1,lotsynk,syn,igrd1p,lotsyn,                 &
     &             jgrd1p,kpxp,kpx,1)
      call rmpnk2nx(synk,igrd1,lotsynk,syn,igrd1p,lotsyn,                &
     &             jgrd1p,levrp,levr,kuup,kuu,5)
      call rmpx2nx(synk,igrd1,lotsynk,syn,igrd1p,lotsyn,                 &
     &             jgrd1p,kppp,kpp,1)
      call rmpnk2nx(synk,igrd1,lotsynk,syn,igrd1p,lotsyn,                &
     &             jgrd1p,levrp,levr,kttp,ktt,2)
      call rmpx2nx(synk,igrd1,lotsynk,syn,igrd1p,lotsyn,                 &
     &             jgrd1p,kpmp,kpm,1)
#ifdef NONHYD
#ifdef NDSLDYN
      call rmpnk2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,               &
     &             jgrd1p,levrp,levr,knpyp,knpy,1)
#else
      call rmpnk2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,               &
     &             jgrd1p,levrp,levr,knpyp,knpy,3)
      call rmpx2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,                &
     &             jgrd1p,knwyp+levrp,knwy+levr,1)
#endif
#ifdef NDSLDYN
      call rmpnk2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,               &
     &             jgrd1p,levrp,levr,knpxp,knpx,1)
#else
      call rmpnk2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,               &
     &             jgrd1p,levrp,levr,knpxp,knpx,3)
      call rmpx2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,                &
     &             jgrd1p,knwxp+levrp,knwx+levr,1)
#endif
      call rmpnk2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,               &
     &             jgrd1p,levrp,levr,knppp,knpp,3)
      call rmpx2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,                &
     &             jgrd1p,knwwp+levrp,knww+levr,1)
      call rmpnk2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,               &
     &             jgrd1p,levrp,levr,knpmp,knpm,3)
      call rmpx2nx(snnyk,igrd1,lotsnnk,snn,igrd1p,lotsnn,                &
     &             jgrd1p,knwmp+levrp,knwm+levr,1)
#endif
!
#ifdef CNT_1
        call cnttime_off("one_loop_spectal")
#endif
!
      lat1s=1
      lat2s=latlen(mype)
      latdons=0
      latts=latstr(mype)-1
#ifdef OPENMP
!$omp parallel do private(lat,lan,latt,latj,i,j,k,n,kk)
#endif
      do 2052 lat=lat1s,lat2s
        lan=lat-latdons
        latt=lat+latts
        latj = (lat-1)*igrd1p     
!
! ........obtain base field values
!
      call pint3(qb  ,bgd(1,kppb,lan),                                   &  
     &                bgd(1,kpxb,lan),                                   &  
     &                bgd(1,kpyb,lan),wsltb,dwsltb,latt,numbgd,          &
     &                igrd1p,                                            &
     &                ib1,ib2,jb1,jb2,jbx,istr,ilen,                     &
     &                bgf,border,jgrd1)
!
! .......obtain full field values
!
!       call maxmin(syn(1,ktt,lan),lonlens,levr,1,levr,'aft pint3 tmp')
!       call maxmin(syn(1,ktt,lan),lonlens,levr,1,levr,'aft pint3 tmp')
!       call maxmin(bgd(1,kttb,lan),lonlens,levr,1,levr,'aft pint3 bgd')
!       call maxmin(snn(1,kntt,lan),lonlens,levr,1,levr,'aft pint3 nt')
!       call maxmin(snn(1,knww,lan),lonlens,levr,1,levr,'aft pint3 nw')
!       call maxmin(snn(1,knpp,lan),lonlens,levr,1,levr,'aft pint3 np')

      do 310 i=1,lonlens
      syn(i,kpp,lan)=syn(i,kpp,lan) + bgd(i,kppb,lan)
      syn(i,kpy,lan)=syn(i,kpy,lan) + bgd(i,kpyb,lan)
      syn(i,kpx,lan)=syn(i,kpx,lan) + bgd(i,kpxb,lan)
 310  continue
      do 320 k=1,levr
      do 320 i=1,lonlens
      syn(i,ktt-1+k,lan)=syn(i,ktt-1+k,lan)+bgd(i,kttb-1+k,lan)
      syn(i,kuu-1+k,lan)=syn(i,kuu-1+k,lan)+bgd(i,kuub-1+k,lan)
      syn(i,kvv-1+k,lan)=syn(i,kvv-1+k,lan)+bgd(i,kvvb-1+k,lan)
      syn(i,kux-1+k,lan)=syn(i,kux-1+k,lan)+bgd(i,kuxb-1+k,lan)
      syn(i,kuy-1+k,lan)=syn(i,kuy-1+k,lan)+bgd(i,kuyb-1+k,lan)
      syn(i,kvx-1+k,lan)=syn(i,kvx-1+k,lan)+bgd(i,kvxb-1+k,lan)
      syn(i,kvy-1+k,lan)=syn(i,kvy-1+k,lan)+bgd(i,kvyb-1+k,lan)
      syn(i,ktx-1+k,lan)=syn(i,ktx-1+k,lan)+bgd(i,ktxb-1+k,lan)
      syn(i,kty-1+k,lan)=syn(i,kty-1+k,lan)+bgd(i,ktyb-1+k,lan)
#ifdef RKN
      udx(i,k,lan)=syn(i,kux-1+k,lan)
      udy(i,k,lan)=syn(i,kuy-1+k,lan)
      vdx(i,k,lan)=syn(i,kvx-1+k,lan)
      vdy(i,k,lan)=syn(i,kvy-1+k,lan)
#endif
 320  continue
      do 325 n=1,ntotal
      do 325 k=1,levr
        kk=(n-1)*ntotal+k
      do 325 i=1,lonlens
        srq(i,kk,lan)=rqgrd(i+latj,k,n)+bgd(i,kqqb-1+kk,lan)
 325  continue
!
!jun radiation
      if (lsphy) then

      do j=1,lonlens
          synr(j,kppr,lan) = syn(j,kpp,lan)
          synr(j,kpxr,lan) = xm2(j,lat)*syn(j,kpx,lan)
          synr(j,kpyr,lan) = xm2(j,lat)*syn(j,kpy,lan)
      enddo
      do k=1,levr
      do j=1,lonlens
          synr(j,kttr-1+k,lan)=syn(j,ktt-1+k,lan)
          synr(j,kuur-1+k,lan)=syn(j,kuu-1+k,lan)
          synr(j,kvvr-1+k,lan)=syn(j,kvv-1+k,lan)
          synr(j,kddr-1+k,lan)=xm2(j,lat)*                               &    
     &            (syn(j,kux-1+k,lan)+syn(j,kvy-1+k,lan))
      enddo
      enddo
      do k=1,levh
      do j=1,lonlens
         synr(j,kqqr-1+k,lan)=srq(j,k,lan)
      enddo
      enddo

        if(itimsw.ne.0 .or. itimlw.ne.0) then
         call radiation_run(lat,lan,lotsynr,                             &
     &    kpyr,kvvr,kpxr,kuur,kddr,kppr,kttr,kqqr,                       &  
#ifdef NONHYD
     &    snn(1,knww,lan),                                               &
#endif
     &    synr,phy_f3d,fluxr,slmsk,                                      &
     &    hprime(1,1,1),stc(1,1,1),sheleg,zorl,cv,cvt,cvb,               &
     &    alvsf,alnsf,alvwf,alnwf,facsf,facwf,                           &
     &    swh,hlw,sfcnsw,sfcdlw,coszer,tsflw,cvavg,del                   &
#ifdef RKN
     &    ,flwupf(1,lan),fswupf(1,lan),fswdnf(1,lan),                    &
     &    sswupf(1,lan),sswdnf(1,lan),                                   &
     &    slwupf(1,lan),slwdnf(1,lan),                                   &
     &    cldsaf(1,1,lan),mtopaf(1,1,lan),mbotaf(1,1,lan)                &
#endif
     &    ,nmmiph,cldcov                                                 &
!yj---get sfc downward short wave radiation and snowfree albedo--->>>     
     &    ,sswdnf(1,lan),sfalb(1,lan))
!yj---get sfc downward short wave radiation and snowfree albedo---<<<     
        endif
      endif
! end radiation
 2052 continue
!
!yj---move n-1 full field to here---------
!
#ifdef OPENMP
!$omp parallel do private(lat,lan,latt,latj,i,j,k,n,kk,dt)
#endif
      do 2054 lat=lat1s,lat2s
        lan=lat-latdons
        latt=lat+latts
        latj = (lat-1)*igrd1p

!cc compute n-1 partial total field temp,u,v,q,p 
!yj: step 1: n-1 perturbation + n base field
!
      do i=1,lonlens
      synm(i,kppm,lan) = syn(i,kpm,lan)+bgd(i,kppb,lan)
      enddo
      do k=1,levr
      do i=1,lonlens
      synm(i,kuum-1+k,lan)=syn(i,kum-1+k,lan)+bgd(i,kuub-1+k,lan)
      synm(i,kvvm-1+k,lan)=syn(i,kvm-1+k,lan)+bgd(i,kvvb-1+k,lan)
      synm(i,kttm-1+k,lan)=syn(i,ktm-1+k,lan)+bgd(i,kttb-1+k,lan)
      enddo
      enddo
      do n=1,ntotal
      do  k=1,levr
        kk=(n-1)*levr+k
      do  i=1,lonlens
        synm(i,kqqm-1+kk,lan)=rmgrd(i+latj,k,n)+bgd(i,kqqb-1+kk,lan)
      enddo
      enddo
      enddo
!
! ...... get base field tendencies ....
!
      call pint1(qltb,bgd(1,kplb,lan),wsltb,latt,numbgd,igrd1p,          &
     &           ib1,ib2,jb1,jb2,jbx,istr,ilen,                          &
     &           bgf,border,jgrd1)
!
!
!yj: step 2: (n-1 perturbation + n base field) - (base tend*dt)
!yj:          = n-1 full field
! update total field: add base field increment to nstar
      dt=dt2/2.
!NDSL 3D ln(ps) >      
      do i=1,lonlens
      synm(i,kppm,lan)= -bgd(i,kplb,lan)*dt+synm(i,kppm,lan)
      enddo
      do k=2,levr
      synm(:,kppm-1+k,lan)= synm(:,kppm,lan)
      enddo
!NDSL 3D ln(ps) <      
      do k=1,levr
      do i=1,lonlens
      synm(i,kuum-1+k,lan)=-bgd(i,kulb-1+k,lan)*dt+synm(i,kuum-1+k,lan)
      synm(i,kvvm-1+k,lan)=-bgd(i,kvlb-1+k,lan)*dt+synm(i,kvvm-1+k,lan)
      synm(i,kttm-1+k,lan)=-bgd(i,ktlb-1+k,lan)*dt+synm(i,kttm-1+k,lan)
      enddo
      enddo
              
      do k=1,levh
      do i=1,lonlens
      synm(i,kqqm-1+k,lan)=-bgd(i,kqlb-1+k,lan)*dt+synm(i,kqqm-1+k,lan)
      enddo
      enddo
          
!snnm
#ifdef NONHYD
!get n-1 full field for horizontal NDSL
      do k=1,levr
      do i=1,lonlens
      snnm(i,knppm-1+k,lan)=snn(i,knpm-1+k,lan)+                        &
     &     synm(i,kppm-1+k,lan)+log(sl(k))
      snnm(i,knttm-1+k,lan)=snn(i,kntm-1+k,lan)+synm(i,kttm-1+k,lan)
      enddo
      enddo
!levrp1=levr+1      
      do  k=1,levrp1
      do  i=1,lonlens
        snnm(i,knwwm-1+k,lan)=snn(i,knwm-1+k,lan)
      enddo
      enddo

#ifdef NDSLDYN
!yjadd gz to NDSL >
      do k=1,levr
        do i=1,lonlens
          synm(i,kgzm+k-1,lan)=gzz(i,lan)
        enddo
      enddo
      do k=1,levr
        do kk=1,levr
          do i=1,lonlens
            synm(i,kgzm-1+k,lan)=synm(i,kgzm-1+k,lan)+                 &
     &                           synm(i,kttm-1+kk,lan)*am(k,kk)
          enddo
        enddo
      enddo
!yjadd gz to NDSL <
#endif

#endif
!endif nonhydro.

 2054 continue

!NDSL>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef NDSLDYN
!NDSL_nx2nk >
#ifdef NONHYD
!do ww first, borrow synt to store half layer uu and vv
      synt=0.0
!t=n
#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k)
#endif
      do lat=lat1s,lat2s
        lan=lat-latdons
      do k=1,levr-1
      do i=1,lonlens
        synt(i,kuut+k-1,lan)=  &
     &   0.5*(syn(i,kuu+k-1,lan)+syn(i,kuu+k,lan))*xm2(i,lat)
        synt(i,kvvt+k-1,lan)=  &
     &   0.5*(syn(i,kvv+k-1,lan)+syn(i,kvv+k,lan))*xm2(i,lat)
      enddo
      enddo
      enddo
      call rmpnx2nk(synt,igrd1p,lotsynt,synuk,igrd1,levrp,                &
     &             jgrd1p,levr,levrp,kuut,1,1)
      call rmpnx2nk(synt,igrd1p,lotsynt,synvk,igrd1,levrp,                &
     &             jgrd1p,levr,levrp,kvvt,1,1)
!t=n-1
      w2nd(:,:,:)=snnm(:,knwwm+1:knwwm+levr,:)
      call rmpnx2nk(w2nd,igrd1p,levr,synmk,igrd1,lotsynmk,              &
     &             jgrd1p,levr,levrp,1,kwmk,1)
      k1=kwmk
      k2=kwmk+levrp-1
      call rndslfv_monoadvh2(synmk,levrp,k1,k2,lotsynmk,      &
     &                     synuk,synvk,deltim,rdelx,rdely,difhsec)

!note: synmk(kwmk)=2nd layer of w, only 2~levr+1 of w is calculated 
!note: but only 2~levr of w should be used in nfidi_rloopmp
#endif

!other variables:
!t=n
#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k)
#endif
      do lat=lat1s,lat2s
        lan=lat-latdons
      do k=1,levr
        do i=1,lonlens
          synt(i,kuut+k-1,lan)=syn(i,kuu+k-1,lan)*xm2(i,lat)
          synt(i,kvvt+k-1,lan)=syn(i,kvv+k-1,lan)*xm2(i,lat)
        enddo
      enddo
      enddo
      call rmpnx2nk(synt,igrd1p,lotsynt,synuk,igrd1,levrp,              &
     &             jgrd1p,levr,levrp,kuut,1,1)
      call rmpnx2nk(synt,igrd1p,lotsynt,synvk,igrd1,levrp,              &
     &             jgrd1p,levr,levrp,kvvt,1,1)
      synt=0.0
!t=n-1
      call rmpnx2nk(synm,igrd1p,lotsynm,synmk,igrd1,lotsynmk,                &
     &             jgrd1p,levr,levrp,kuum,kumk,4+ntotal)
#ifdef NONHYD
      call rmpnx2nk(synm,igrd1p,lotsynm,synmk,igrd1,lotsynmk,                &
     &             jgrd1p,levr,levrp,kgzm,kgzmk,1)
      call rmpnx2nk(snnm,igrd1p,lotann,synmk,igrd1,lotsynmk,              &
     &             jgrd1p,levr,levrp,knppm,kpnmk,2)
#endif
!NDSL_nx2nk <
!NDSL-horiadv >
#ifdef NONHYD
      k1=kumk
      k2=ktnmk+levrp-1
      call rndslfv_monoadvh2(synmk,levrp,k1,k2,lotsynmk, &
     &                      synuk,synvk,deltim,rdelx,rdely,difhsec)
#else
      k1=kumk
      k2=kqmk+levhp-1
      call rndslfv_monoadvh2(synmk,levrp,k1,k2,lotsynmk, &
     &                      synuk,synvk,deltim,rdelx,rdely,difhsec)
#endif
!NDSL-horiadv <
!NDSL_nk2nx >
      call rmpnk2nx(synmk,igrd1,lotsynmk,synt,igrd1p,lotsynt,             &
     &             jgrd1p,levrp,levr,kumk,kuut,4+ntotal)
#ifdef NONHYD
      call rmpnk2nx(synmk,igrd1,lotsynmk,synt,igrd1p,lotsynt,             &
     &             jgrd1p,levrp,levr,kgzmk,kgzt,1)
      call rmpnk2nx(synmk,igrd1,lotsynmk,synt,igrd1p,lotsynt,             &
     &             jgrd1p,levrp,levr,kpnmk,kpnt,3)
#endif
!NDSL_nk2nx <
!yj:add index for lateral boundary>
      if (llind) then
      k1=1
      k2=levrp
      rindex=1.0
      call rndslfv_monoadvh2_ind(rindex,levrp,k1,k2,levrp,              &
     &                      synuk,synvk,deltim)
      call rmpnk2nx(rindex,igrd1,levrp,pindex,igrd1p,levr,              &
     &             jgrd1p,levrp,levr,1,1,1)
!      if (mype.eq.master) then
!        call maxmin(pindex,lngrdp,levr,1,levr,'rloopmp_rndsl pindex')
!      endif
      else
        pindex=1.0
      endif
!
!      if (llind) then
!         if (nint(fhour*3600.).le.5.0*deltim) then
!         do j=1,jgrd1p
!           jlat=(j-1)*igrd1p
!           do k=1,levr
!             do i=1,igrd1p
!               ij=i+jlat
!               tmpgrd(ij,k)=pindex(i,k,j)
!             enddo
!           enddo
!         enddo
!         call rmpgp2f(tmpgrd,igrd1p,jgrd1p,indfull,igrd1,jgrd1,levr)
!
!          if( mype.eq.master ) then
!!          call maxmin(indfull,lngrd,levr,1,levr,'rloopmp_rndsl pindex')
!            write(cfhour,'(i4.4)') nint(fhour*3600.)
!            print*,'cfhour for pindex 1st10=',cfhour
!            open(98,file='pindex_1st10_'//cfhour,form='unformatted',   &
!     &           access='direct',recl=lngrd*4)
!            do k=1,levr
!              call ful2sgl(indfull(1,k),inds,lngrd)
!              write(98,rec=k) (inds(i),i=1,lngrd)
!            enddo   
!            close(98)
!          endif
!         endif
!      endif
!yj:add index for lateral boundary<


#else
!if not doing NDSLDYN

!t=n
#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k)
#endif
      do lat=lat1s,lat2s
        lan=lat-latdons
      do k=1,levr
        do i=1,lonlens
          utmp(i,k,lan)=syn(i,kuu+k-1,lan)*xm2(i,lat)
          vtmp(i,k,lan)=syn(i,kvv+k-1,lan)*xm2(i,lat)
        enddo
      enddo
      enddo
      call rmpnx2nk(utmp,igrd1p,levr,synuk,igrd1,levrp,              &
     &             jgrd1p,levr,levrp,1,1,1)
      call rmpnx2nk(vtmp,igrd1p,levr,synvk,igrd1,levrp,              &
     &             jgrd1p,levr,levrp,1,1,1)
!t=n-1
      call rmpnx2nk(synm,igrd1p,lotsynm,synmk,igrd1,lotsynmk,           &
     &             jgrd1p,levr,levrp,kqqm,kqmk,ntotal)
!NDSL-horiadv >
      k1=kqmk
      k2=kqmk+levhp-1
      call rndslfv_monoadvh2(synmk,levrp,k1,k2,lotsynmk, &
     &                      synuk,synvk,deltim,rdelx,rdely,difhsec)
!NDSL-horiadv <
!NDSL_nk2nx >
      call rmpnk2nx(synmk,igrd1,lotsynmk,synt,igrd1p,lotsynt,             &
     &             jgrd1p,levrp,levr,kqmk,kqqt,ntotal)
!NDSL_nk2nx <
!yj:add index for lateral boundary>
      if (llind) then
      k1=1
      k2=levrp
      rindex=1.0
      call rndslfv_monoadvh2_ind(rindex,levrp,k1,k2,levrp,              &
     &                      synuk,synvk,deltim)
      call rmpnk2nx(rindex,igrd1,levrp,pindex,igrd1p,levr,              &
     &             jgrd1p,levrp,levr,1,1,1)
!      if (mype.eq.master) then
!        call maxmin(pindex,lngrdp,levr,1,levr,'rloopmp_rndsl pindex')
!      endif
      else
        pindex=1.0
      endif
!yj:add index for lateral boundary<

#endif
!NDSL<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
!yj2020check>
!     if (mype.eq.0) print*,'t=n,srq(1,1,1)=',srq(1,1,1)
!yj2020check<

!
#ifdef OPENMP
!$omp parallel do private(lat,lan,latt)
#endif
      do 2053 lat=lat1s,lat2s
        lan=lat-latdons
        latt=lat+latts
! ......... compute full field forcing
!
#ifdef NONHYD
      call nfidi_rndsl(syn(1,kuu,lan),syn(1,kvv,lan),                    &
     &           syn(1,ktt,lan),srq(1,  1,lan),                          &
     &           syn(1,kpx,lan),syn(1,kpy,lan),                          &
     &           syn(1,kux,lan),syn(1,kuy,lan),                          &
     &           syn(1,kvx,lan),syn(1,kvy,lan),                          &
     &           syn(1,ktx,lan),syn(1,kty,lan),                          &
     &           syn(1,kpp,lan),spdlat(1,lan),                           &
     &           snn(1,knpp,lan),snn(1,kntt,lan),snn(1,knww,lan),        &
#ifdef NDSLDYN
     &           snn(1,knpx,lan),                                        &
     &           snn(1,knpy,lan),                                        &
#else
     &           snn(1,knpx,lan),snn(1,kntx,lan),snn(1,knwx,lan),        &
     &           snn(1,knpy,lan),snn(1,knty,lan),snn(1,knwy,lan),        &
#endif
     &           snn(1,knwm+levr,lan),deltim,                            &
#ifndef INTDTB     
     &           bgd(1,ktlb,lan),                                        &
#endif
     &           si,sl,p1,p2,am,                                         &
     &           xm2(1,lat),xm2px(1,lat),xm2py(1,lat),                   &
     &           corf(1,lat),corf2(1,lat),gzdx(1,lat),gzdy(1,lat),       &
     &           del,rdel2,lonlens,igrd1p,levr,ntotal                    &
!yj NDSL >
#ifdef NDSLDYN
     &          ,synt(1,kuut,lan),synt(1,kvvt,lan)                       &
     &                           ,synt(1,kppt,lan)                       &
     &          ,synt(1,kqqt,lan),synt(1,kgzt,lan)                       &
     &          ,synt(1,kpnt,lan),synt(1,ktnt,lan)                       &
     &          ,synt(1,kwwt,lan)                                        &
     &          ,dot(1,1,lan)                                            &
#ifdef INTDTB     
     &          ,synt(1,kttt,lan),synm(1,kttm,lan)                       &
#endif
!yj NDSL < 
#else
     &          ,synt(1,kqqt,lan)                                        &
     &          ,dot(1,1,lan)                                            &
#endif
     &    )
#else
      call rfidi_rndsl(syn(1,kuu,lan),syn(1,kvv,lan),                    & 
     &           syn(1,ktt,lan),srq(1,  1,lan),                          &
     &           syn(1,kpx,lan),syn(1,kpy,lan),                          &
     &           syn(1,kux,lan),syn(1,kuy,lan),                          &
     &           syn(1,kvx,lan),syn(1,kvy,lan),                          &
     &           syn(1,ktx,lan),syn(1,kty,lan),                          &
     &           syn(1,kpp,lan),spdlat(1,lan),                           &
     &     p1,p2,am,xm2(1,lat),xm2px(1,lat),xm2py(1,lat),                &
     &     corf(1,lat),gzdx(1,lat),gzdy(1,lat),del,rdel2,                &
     &     lonlens,igrd1p,levr,levrm1,levrp1,ntotal                      &
#ifdef T
     &    ,at,bt,ct                                                      &
#endif
#ifdef Q
     &    ,aq,bq                                                         &
#endif
#ifdef U
     &    ,au,bu,fu,pu                                                   &
#endif
#ifdef V
     &    ,av,bv,fv,pv                                                   &
#endif
#ifdef P
     &    ,ap,bp                                                         &
#endif
!yj NDSL >
#ifdef NDSLDYN
     &          ,synt(1,kuut,lan),synt(1,kvvt,lan)                       &
     &          ,synt(1,kttt,lan),synt(1,kppt,lan)                       &
     &          ,synt(1,kqqt,lan),sl                                     &
     &    ,dot(1,1,lan)                                                  &
#else
     &          ,synt(1,kqqt,lan)                                        &
     &    ,dot(1,1,lan)                                                  &
#endif
!yj NDSL < 
     &    )
#endif
! end radiation
!yj NDSL> 
!yj call maxmin(dot(1,1,lan),igrd1p,levr+1,1,levr+1,'rloopmp_rndsl sdot grid')
!yj NDSL< 
 2053 continue

!yj2020check>
!     if (mype.eq.0) print*,'t=n,tend,srq(1,1,1)=',srq(1,1,1)
!yj2020check<

!
!yj-----------cut to be 2055--------------------------
#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k)
#endif
      do 2055 lat=lat1s,lat2s
        lan=lat-latdons

!update to nstar
      do i=1,lonlens
      syn(i,kpp,lan) = syn(i,kpp,lan)*dt2+synm(i,kppm,lan)
      enddo
      do k=1,levr
      do i=1,lonlens
      syn(i,kuu-1+k,lan)=syn(i,kuu-1+k,lan)*dt2+synm(i,kuum-1+k,lan)
      syn(i,kvv-1+k,lan)=syn(i,kvv-1+k,lan)*dt2+synm(i,kvvm-1+k,lan)
      syn(i,ktt-1+k,lan)=syn(i,ktt-1+k,lan)*dt2+synm(i,kttm-1+k,lan)
      enddo
      enddo
      do  k=1,levh
      do  i=1,lonlens
      srq(i,k,lan)=srq(i,k,lan)*dt2+synm(i,kqqm-1+k,lan)
      enddo
      enddo

#ifdef NONHYD
!for vertical NDSL, tn: full field, pn: nonhyd. perturbation      
      do k=1,levr
      do i=1,lonlens
!yj      snn(i,kntt-1+k,lan)=snn(i,kntt-1+k,lan)*dt2+snnm(i,knttm-1+k,lan)
      snn(i,kntt-1+k,lan)=snn(i,kntt-1+k,lan)*dt2+snn(i,kntm-1+k,lan)   &
     &                   +syn(i,ktt-1+k,lan)
#ifdef NDSLDYN
!yj      snn(i,knpp-1+k,lan)=snn(i,knpp-1+k,lan)*dt2  &
!yj     &           +snnm(i,knppm-1+k,lan) -synm(i,kppm-1+k,lan) -log(sl(k))
      snn(i,knpp-1+k,lan)=snn(i,knpp-1+k,lan)*dt2+snn(i,knpm-1+k,lan)  
!yj-v      snn(i,knpp-1+k,lan)=snn(i,knpp-1+k,lan)*dt2+snnm(i,knppm-1+k,lan)
#else
      snn(i,knpp-1+k,lan)=snn(i,knpp-1+k,lan)*dt2+snn(i,knpm-1+k,lan)   &
     &                   +syn(i,kpp,lan)+log(sl(k))
#endif
      enddo
      enddo
      do k=1,levr+1
      do i=1,lonlens
      snn(i,knww-1+k,lan)=snn(i,knww-1+k,lan)*dt2+snnm(i,knwwm-1+k,lan)
      enddo
      enddo
#endif

#ifdef NDSLDYN
!yj NDSL >
!yj only i,k (j is the lan loop)
!yj: for non-hyd do nonhyd. affected factors, u,v,tracers
!yj: u,v,(t,if hyd),tracers
#ifndef NONHYD
      !u+v+T+tracers 
      call rndslfv_monoadvv(syn(1,kuu,lan),lonlens,    &
     &      si,dot(1,1,lan),deltim)
      call rndslfv_monoadvv(syn(1,kvv,lan),lonlens,    &
     &      si,dot(1,1,lan),deltim)
      call rndslfv_monoadvv(syn(1,ktt,lan),lonlens,    &
     &      si,dot(1,1,lan),deltim)
      do k=1,ntotal
        call rndslfv_monoadvv(srq(1,1+(k-1)*levr,lan),lonlens,  &
     &      si,dot(1,1,lan),deltim)
      enddo
#else
      !u+v+tracers 
      call rndslfv_monoadvv(syn(1,kuu,lan),lonlens,    &
     &      si,dot(1,1,lan),deltim)
      call rndslfv_monoadvv(syn(1,kvv,lan),lonlens,    &
     &      si,dot(1,1,lan),deltim)
      do k=1,ntotal
        call rndslfv_monoadvv(srq(1,1+(k-1)*levr,lan),lonlens,  &
     &      si,dot(1,1,lan),deltim)
      enddo
!yj: pn,tn,ww      
!     do i=1,lonlens
!       snn(i,knww,lan)=(syn(i,kuu,lan)*gzdx(i,lat)   &
!    &   +syn(i,kvv,lan)*gzdy(i,lat))*xm2(i,lat)/grav
!     enddo
      call rndslfv_monoadvv_nonhyd(snn(1,1,lan),lonlens,    &
     &      si,dot(1,1,lan),deltim) 
!yj: return pn full field      
      do k=1,levr
      do i=1,lonlens
      snn(i,knpp-1+k,lan)=snn(i,knpp-1+k,lan)  &
     &        +syn(i,kpp,lan) +log(sl(k))
      enddo
      enddo
!yj NDSL <     
#endif

#else
!if not doing NDSLDYN
      !tracers 
      do k=1,ntotal
      call rndslfv_monoadvv(srq(1,1+(k-1)*levr,lan),lonlens,  &
     &      si,dot(1,1,lan),deltim)
      enddo
#endif

      if (llind) then
!      do k=1,levr
!      do i=1,lonlens
!        if (pindex(i,k,lan).lt.1.0) then
!          pindex(i,k,lan)=0.0
!        else
!          pindex(i,k,lan)=1.0
!        endif
!      enddo
!      enddo
      do i=1,lonlens
        if (sum(pindex(i,1:levr,lan)).lt.(levr*1.0-0.1)) then
          ppindex(i,lan)=0.0
        endif
      enddo
      endif

 2055 continue

!yj2020check>
!     if (mype.eq.0) print*,'t=n+1*,srq(1,1,1)=',srq(1,1,1)
!yj2020check<

!yj> check vertical index
!      if (llind) then
!         do j=1,jgrd1p
!           jlat=(j-1)*igrd1p
!           do k=1,levr
!             do i=1,igrd1p
!               ij=i+jlat
!               tmpgrdv(ij,k)=pindex(i,k,j)
!             enddo
!           enddo
!         enddo
!        call rmpgp2f(tmpgrdv,igrd1p,jgrd1p,indfullv,igrd1,jgrd1,levr)
!        if( mype.eq.master ) then
!          checkdif=0.0
!          do k=1,levr
!            do ij=1,lngrd
!              checkdif=checkdif+abs(indfullv(ij,k))-abs(indfull(ij,k))
!           enddo
!         enddo
!         print*,'at fhour, diff btw ndsl vertical and horizontal=',   &
!     &     cfhour,checkdif
!        endif
!
!      endif
!yj< check vertical index

!
!****************** physics
      if(lsphy) then !yj if physics

#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k)
#endif
      do lat=lat1s,lat2s
        lan=lat-latdons
        do i=1,lonlens
        syn(i,kpx,lan) = xm2(i,lat)*syn(i,kpx,lan)
        syn(i,kpy,lan) = xm2(i,lat)*syn(i,kpy,lan)
        enddo
        do k=1,levr
        do i=1,lonlens
        syn(i,kuu-1+k,lan)=xm (i,lat)*syn(i,kuu-1+k,lan)
        syn(i,kvv-1+k,lan)=xm (i,lat)*syn(i,kvv-1+k,lan)
                                     !yj:virtual wind
        enddo
        enddo
        do k=1,levr
        do i=1,lonlens
        syndd(i,k,lan)=xm2(i,lat)*(syn(i,kux-1+k,lan)+                     &  
     &                               syn(i,kvy-1+k,lan))
        enddo
        enddo
      enddo

#ifdef OPENMP
!$omp parallel do private(lat,lan)
#endif
!       print *,'physics_run is called,me=',me,'lats1=',lat1s,           &
!     &  'lat2s=',lat2s,'latdons=',latdons,'latt=',latt,'latts=',latts

      do 2056 lat=lat1s,lat2s
        lan=lat-latdons

        call physics_run(lan,lat,                                        &
     &    syn(1,kpx,lan),syn(1,kpy,lan),                                 &  
     &    syn(1,kuu,lan),syn(1,kvv,lan),syn(1,kpp,lan),                  &
     &    syn(1,ktt,lan),srq(1,  1,lan),syndd(1,1,lan),                  &
     &    anl(1,kts,lan),anl(1,kqs,lan),anl(1,kus,lan),anl(1,kvs,lan),   &
#ifdef NONHYD
     &    snn(1,knpp,lan),snn(1,knww,lan),snn(1,kntt,lan),               &
     &    ann(1,knppa,lan),ann(1,knwwa,lan),ann(1,kntta,lan),            &
#endif
#ifdef RKN
     &    psexp(1,lan),radsl(1,lan),                                     &
     &    dlwsf1(1,lan),qss(1,lan),                                      &
     &    gflx(1,lan),cd(1,lan),                                         &
     &    cdq(1,lan),evap(1,lan),                                        &
     &    hflx(1,lan),drain(1,lan),                                      &
     &    runof(1,lan),cld1d(1,lan),                                     &
     &    gamt(1,lan),gamq(1,lan),                                       &
     &    dqsfc1(1,lan),dtsfc1(1,lan),                                   &
     &    dusfc1(1,lan),dvsfc1(1,lan),                                   &
     &    dusfcg(1,lan),dvsfcg(1,lan),                                   &
     &    rainc(1,lan),rainl(1,lan),                                     &
     &    hsw(1,1,lan),vvel(1,1,lan),                                    &
     &    stsoil(1,1,lan),smsoil(1,1,lan),                               &
#endif
     &   cldcov,                                                         &
     &   phy_f3d,phy_f2d)

 2056 continue

#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k)
#endif
      do lat=lat1s,lat2s
        lan=lat-latdons
      do k=1,levr
      do i=1,lonlens
      anl(i,kus-1+k,lan)=anl(i,kus-1+k,lan)/xm(i,lat)
      anl(i,kvs-1+k,lan)=anl(i,kvs-1+k,lan)/xm(i,lat)
      !yj:virtual wind
      enddo
      enddo
      enddo
      
      else   !yj else physics

#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k)
#endif
      do lat=lat1s,lat2s
        lan=lat-latdons
        do k=1,levr
        do i=1,lonlens
         anl(i,kus-1+k,lan)=syn(i,kuu-1+k,lan)
         anl(i,kvs-1+k,lan)=syn(i,kvv-1+k,lan)
         anl(i,kts-1+k,lan)=syn(i,ktt-1+k,lan)
#ifdef NONHYD
         ann(i,kntta-1+k,lan)=snn(i,kntt-1+k,lan)
         ann(i,knppa-1+k,lan)=snn(i,knpp-1+k,lan)
#endif
        enddo
        enddo

        do k=1,levh
        do i=1,lonlens
          anl(i,kqs-1+k,lan)=srq(i,k,lan)
        enddo
        enddo

#ifdef NONHYD
        do k=1,levr+1
        do i=1,lonlens
          ann(i,knwwa-1+k,lan)=snn(i,knww-1+k,lan)
        enddo
        enddo
#endif

      enddo 

      endif  !yj-endif physics
!yj....................................................
!
#ifdef RAYLEIFR      
!yj: damp wind with rayleigh friction --->   
#ifdef OPENMP
!$omp parallel do private(k,lat)
#endif
      do k=1,levr
      spdmax(k) = 0.0
       do lat=lat1s,lat2s
        spdmax(k)=max(spdmax(k),spdlat(k,lat))
       enddo
      enddo
      call rmpgetspd(levr,spdmax)

#ifdef OPENMP
!$omp parallel do private(k,lat,lan)
#endif
      do k=levr-levtop,levr
      do lat=lat1s,lat2s
        lan=lat-latdons
        if (spdmax(k) .gt. windmax) then
          if( mype.eq.master ) then
            print*,'call rayleifr, k,spdmax=',k,spdmax(k)
          endif
          call rayleifr(igrd1p,lonlens,levr,levtop                     &
     &              ,anl(1,kus,lan),anl(1,kvs,lan))
        endif
      enddo
      enddo
!yj: damp wind with rayleigh friction ---<   
#endif

!
#ifdef OPENMP
!$omp parallel do private(lat,lan,i,k)
#endif
      do 2058 lat=lat1s,lat2s
        lan=lat-latdons

#ifdef NONHYD
! for pn, tn      
!yj: nonhyd. pert.= total field - hydro.
!yj: nonhyd. pert. tendency = (nonhyd.pert.(n+1)-nonhyd.pert.(n))/dt2
!yj: note: syn(kts), syn(kpp) are not updated by physics,
!yj:       i.e., anl(kts)=syn(ktt), and syn(kpp) n = syn(kpp) n+1
      do k=1,levr
      do i=1,lonlens
      snn(i,kntt-1+k,lan)=(ann(i,kntta-1+k,lan)-                         &
     &                     anl(i,kts+k-1,lan)-                           &
     &                     snn(i,kntm-1+k,lan))/dt2    
      snn(i,knpp-1+k,lan)=(ann(i,knppa-1+k,lan)-                         &
     &                     syn(i,kpp,lan)-log(sl(k))-                    & 
     &                     snn(i,knpm-1+k,lan))/dt2   
      enddo
      enddo
! for w      
      do k=1,levr+1
      do i=1,lonlens
      snn(i,knww-1+k,lan)=(ann(i,knwwa-1+k,lan)-                          &
     &                     snnm(i,knwwm-1+k,lan))/dt2  
      enddo
      enddo
#endif
      
!yj: note: syn(kpp) updated as tendency (Do nonhyd. before update
!yj:       syn(kpp) to tendency.)
      do i=1,lonlens
      syn(i,kpp,lan) = (syn(i,kpp,lan)-synm(i,kppm,lan))/dt2
      enddo


       do k=1,levr
       do i=1,lonlens
      syn(i,kuu-1+k,lan)=(anl(i,kus-1+k,lan)-                            &  
     &                    synm(i,kuum-1+k,lan))/dt2
      syn(i,kvv-1+k,lan)=(anl(i,kvs-1+k,lan)-                            &
     &                    synm(i,kvvm-1+k,lan))/dt2
      syn(i,ktt-1+k,lan)= (anl(i,kts-1+k,lan)-synm(i,kttm-1+k,lan))      &
     &                    /dt2
      enddo
      enddo
!
      do k=1,levh
      do i=1,lonlens
      srq(i,k,lan)=(anl(i,kqs-1+k,lan)-synm(i,kqqm-1+k,lan))      &
     &                  /dt2
      enddo
      enddo
!

#ifdef RKN
      if(lsphy) then
      do k=1,levr
      do i=1,lonlens
        gu0(i,k,lan)=anl(i,kus-1+k,lan)
        gv0(i,k,lan)=anl(i,kvs-1+k,lan)
        gt0(i,k,lan)=anl(i,kts-1+k,lan)
      enddo
      enddo
      do k=1,levh
      do i=1,lonlens
        gq0(i,k,lan)=anl(i,kqs-1+k,lan)
      enddo
      enddo
      endif
#endif
!
!colse phys
 2058 continue

!yj2020check>
!     if (mype.eq.0) print*,'t=n+1*,total forcing,srq(1,1,1)=',srq(1,1,1)
!yj2020check<

#ifdef RKN
      if(lsphy) then

!......................................................
!      if(me.eq.0) print 100,(spdmax(k),k=1,levr)

!  ken point diagnostics
!
!
       if(isave.ne.0.and.npoint.gt.0) then
         lones=lonlen(mype)
         lon1s=lonstr(mype)
         lon2s=lon1s-1+lones
         latt1=latts+lat1s
         latt2=latts+lat2s
!radiation _run
         call rkenputr(igrd1p,jgrd1p,lon1s,lon2s,                        &
     &        lones,                                                     &
     &        latt1,latt2,sdec,                                          &
     &        slmsk(1,lat1s),tser(1,lat1s),                              &
     &        flwupf(1,1),fswupf(1,1),fswdnf(1,1),                       &
     &        sswupf(1,1),sswdnf(1,1),slwupf(1,1),slwdnf(1,1),           &
     &        cldsaf(1,1,1),mtopaf(1,1,1),mbotaf(1,1,1),                 &
     &        itimsw,itimlw)                         
!         call rkenputr(igrd1p,jgrd1p,lon1s,lon2s,                        &
!     &        lones,                                                     &
!     &        itimsw,itimlw)
!
!physics _run
         call rkenput(igrd1p,jgrd1p,lon1s,lon2s,                         &   
     &                    lones,                                         &
     &        latt1,latt2,rcl,                                           &
     &        slmsk(1,lat1s),psexp(1,1),tg3(1,lat1s),sheleg(1,lat1s),    &
     &        radsl(1,1),dlwsf1(1,1),                                    &
     &        tsea(1,lat1s),qss(1,1),gflx(1,1),                          &
     &        zorl(1,lat1s),cd(1,1),cdq(1,1),                            &
     &        evap(1,1),hflx(1,1),stsoil(1,1,1),                         &
     &        canopy(1,lat1s),drain(1,1),smsoil(1,1,1),                  &
     &        runof(1,1),cld1d(1,1),                                     &
     &        u10m(1,lat1s),v10m(1,lat1s),t2m(1,lat1s),q2m(1,lat1s),     &
     &        hpbl(1,lat1s),gamt(1,1),gamq(1,1),                         &
     &        dqsfc1(1,1),dtsfc1(1,1),dusfc1(1,1),dvsfc1(1,1),           &
     &        dusfcg(1,1),dvsfcg(1,1),                                   &
     &        rainc(1,1),rainl(1,1),                                     &
     &        gu0(1,1,1),gv0(1,1,1),                                     &
     &        gt0(1,1,1),gq0(1,1,1),                                     &
     &        hsw(1,1,1),hlw(1,1,lat1s),vvel(1,1,1),sl(1),               & !binbin add sl(1) here
     &        vfrac(1,lat1s),tmpmin(1,lat1s),tmpmax(1,lat1s),            &
     &        bengsh(1,lat1s),geshem(1,lat1s))
!
       endif
!
      endif
#endif

!
      if (rdelx.le.1000) then
       rlxzsec=2400.
      else
       rlxzsec=1800.
      endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,latj,i,j,nc,n,k,kk,rsltb)
#endif
      do 2070 lat=lat1s,lat2s
      lan=lat-latdons
      latj = (lat-1)*igrd1p
!cc
#ifdef A
      do k=1,levr
      do i=1,lonlens
        j=i+latj
#endif
#ifdef T
        wt(j,k,1)=at(i,k)
        wt(j,k,2)=bt(i,k)
        wt(j,k,3)=ct(i,k)
#endif
#ifdef U
        wu(j,k,1)=au(i,k)
        wu(j,k,2)=bu(i,k)
        wu(j,k,3)=fu(i,k)
        wu(j,k,4)=pu(i,k)
#endif
#ifdef V
        wv(j,k,1)=av(i,k)
        wv(j,k,2)=bv(i,k)
        wv(j,k,3)=fv(i,k)
        wv(j,k,4)=pv(i,k)
#endif
#ifdef A
      enddo
      enddo
      do nc=1,ntotal
      do k=1,levr
      do i=1,lonlens
        j=i+latj
#endif
#ifdef Q
        wq(j,k,1)=aq(i,k,nc)
        wq(j,k,2)=bq(i,k,nc)
#endif
#ifdef A
      enddo
      enddo
      enddo
#endif
#ifdef P
      do i=1,lonlens
        j=i+latj
        wp(j,1,1)=ap(i)
        wp(j,1,2)=bp(i)
      enddo
#endif
!ccjun
#ifdef A
      do k=1,levr
      do i=1,lonlens
      j=i+latj
#endif
#ifdef T
      wt(j,k,13)=syn(i,ktt-1+k,lan)-bgd(i,ktlb-1+k,lan)
#endif
#ifdef U
      wu(j,k, 9)=syn(i,kuu-1+k,lan)-bgd(i,kulb-1+k,lan)
#endif
#ifdef V
      wv(j,k, 9)=syn(i,kvv-1+k,lan)-bgd(i,kvlb-1+k,lan)
#endif
#ifdef A
      enddo
      enddo
      do k=1,levh
      do i=1,lonlens
      j=i+latj
#endif
#ifdef Q
      wq(j,k, 8)=srq(i,k,lan)-bgd(i,kqlb-1+k,lan)
#endif
#ifdef A
      enddo
      enddo
#endif
#ifdef P
      do i=1,lonlens
      j=i+latj
      wp(j,1, 4)=syn(i,kpp,lan)-bgd(i,kplb,lan)
      enddo
#endif
!
! ......... adding relaxation and damping to total forcing
! ......... compute perturbation forcing by (full - base) forcing
!

      do i=1,lonlens
        syn(i,kpp,lan)=(syn(i,kpp,lan)-bgd(i,kplb,lan))                  &  
     &                 -syn(i,kpm,lan)*sltb(i+latj)/rlxmsec
      enddo

      do k=1,levr
        if (k.gt.levt) then
          rsltb=1.-cos(pi*(sl(levt)-sl(k))/(sl(levt)-sl(levr)))
          rsltb=rsltb/(2.*rlxzsec)
        else
          rsltb=0.
        endif

      do i=1,lonlens
#ifdef UPDAMP
      syn(i,kuu-1+k,lan)=(syn(i,kuu-1+k,lan)-bgd(i,kulb-1+k,lan))        &
     &                   -syn(i,kum-1+k,lan)*sltb(i+latj)/rlxmsec        &
     &                   -syn(i,kum-1+k,lan)*rsltb
      syn(i,kvv-1+k,lan)=(syn(i,kvv-1+k,lan)-bgd(i,kvlb-1+k,lan))        &
     &                   -syn(i,kvm-1+k,lan)*sltb(i+latj)/rlxmsec        &
     &                   -syn(i,kvm-1+k,lan)*rsltb
#else
      syn(i,kuu-1+k,lan)=(syn(i,kuu-1+k,lan)-bgd(i,kulb-1+k,lan))        &
     &                   -syn(i,kum-1+k,lan)*sltb(i+latj)/rlxmsec
      syn(i,kvv-1+k,lan)=(syn(i,kvv-1+k,lan)-bgd(i,kvlb-1+k,lan))        &
     &                   -syn(i,kvm-1+k,lan)*sltb(i+latj)/rlxmsec
#endif
#ifdef NONHYD     
#ifdef INTDTB
      syn(i,ktt-1+k,lan)=(syn(i,ktt-1+k,lan)-bgd(i,ktlb-1+k,lan))        &
     &                   -syn(i,ktm-1+k,lan)*sltb(i+latj)/rlxhsec
#else
#ifdef UPDAMP
      syn(i,ktt-1+k,lan)=(syn(i,ktt-1+k,lan)-bgd(i,ktlb-1+k,lan))        &
     &                    -syn(i,ktm-1+k,lan)*rsltb 
#else
      syn(i,ktt-1+k,lan)=(syn(i,ktt-1+k,lan)-bgd(i,ktlb-1+k,lan))        
#endif
#endif
#else
      syn(i,ktt-1+k,lan)=(syn(i,ktt-1+k,lan)-bgd(i,ktlb-1+k,lan))        &
     &                   -syn(i,ktm-1+k,lan)*sltb(i+latj)/rlxhsec
#endif
#ifdef NONHYD     
      snn(i,knpp-1+k,lan)= snn(i,knpp-1+k,lan)                           &
     &                   -snn(i,knpm-1+k,lan)*sltb(i+latj)/rlxhsec       &
     &                   -snn(i,knpm-1+k,lan)*rsltb
      snn(i,kntt-1+k,lan)= snn(i,kntt-1+k,lan)                           &
     &                   -snn(i,kntm-1+k,lan)*sltb(i+latj)/rlxhsec       &
     &                   -snn(i,kntm-1+k,lan)*rsltb
#endif

      enddo
      enddo
      do n=1,ntotal
      do k=1,levr
        kk=(n-1)*levr+k
      do i=1,lonlens
      srq(i,kk,lan)=(srq(i,kk,lan)-bgd(i,kqlb-1+kk,lan))     &
     &                   -rmgrd(i+latj,k,n)*sltb(i+latj)/rlxhsec
      enddo
      enddo
      enddo
2070  continue
!yj2020check>
!     if (mype.eq.0) print*,'t=n+1*,relx1,srq(1,1,1)=',srq(1,1,1)
!yj2020check<
!
!      if (llind) then
!        if (mod(nint(fhour*3600.),1*3600).eq.0) then
!         do j=1,jgrd1p
!           jlat=(j-1)*igrd1p
!           do k=1,levr
!             do i=1,igrd1p
!               ij=i+jlat
!               tmpgrd(ij,k)=pindex(i,k,j)
!             enddo
!           enddo
!         enddo
!         call rmpgp2f(tmpgrd,igrd1p,jgrd1p,indfull,igrd1,jgrd1,levr)
!          if( mype.eq.master ) then
!            write(cfhour,'(i4.4)') nint(fhour)
!            print*,'cfhour for pindex last=',cfhour
!            open(99,file='pindex_'//cfhour,form='unformatted',         &
!     &           access='direct',recl=lngrd*4)
!            do k=1,levr
!              call ful2sgl(indfull(1,k),inds,lngrd)
!              write(99,rec=k) (inds(i),i=1,lngrd)
!            enddo   
!            close(99)
!          endif
!        endif
!      endif

#ifdef OPENMP
!$omp parallel do private(lat,lan,latj,n,nkkc,i,j,k,rsltb,rsltn,rsltbm,rsltbh)
#endif

      do 2080 lat=lat1s,lat2s
      lan=lat-latdons
      latj = (lat-1)*igrd1p
      if( iimprlx.eq.1 ) then
        do n=1,ntotal
        nkkc=(n-1)*levr
        do k=1,levr
          do i=1,lonlens
            rsltb=1./(1.+sltb(i+latj)*2.*deltim/rlxhsec)
            rtgrd(i+latj,k,n)= srq(i,nkkc+k,lan)*rsltb*ppindex(i,lan)
            rtgrd(i+latj,k,n)= rtgrd(i+latj,k,n)*dt2+rmgrd(i+latj,k,n)
          enddo
        enddo
        enddo
!yj2020check>
!     if (mype.eq.0) print*,'t=n+1,rtgrd(1,1,1)=',rtgrd(1,1,1)
!yj2020check<
#ifndef NDSLDYN
          ppindex=1.0
#endif
       
!        levt=levr-6
        do i=1,lonlens
          rsltb=1./(1.+sltb(i+latj)*2.*deltim/rlxhsec)
          syn(i,kpp,lan)= syn(i,kpp,lan)*rsltb*ppindex(i,lan)
        enddo
        do k=1,levr
          if (k.gt.levt) then
            rsltn=1.-cos(pi*(sl(levt)-sl(k))/(sl(levt)-sl(levr)))
            rsltn=rsltn/(2.*rlxzsec)
          else
            rsltn=0.
          endif
          rsltn=1./(1.+2.*deltim*rsltn)
          do i=1,lonlens
            rsltbm=1./(1.+sltb(i+latj)*2.*deltim/rlxmsec)
            rsltbh=1./(1.+sltb(i+latj)*2.*deltim/rlxhsec)
#ifdef UPDAMP
            syn(i,kuu-1+k,lan)= syn(i,kuu-1+k,lan)*rsltbm*rsltn*ppindex(i,lan)
            syn(i,kvv-1+k,lan)= syn(i,kvv-1+k,lan)*rsltbm*rsltn*ppindex(i,lan)
#else
            syn(i,kuu-1+k,lan)= syn(i,kuu-1+k,lan)*rsltbm*ppindex(i,lan)
            syn(i,kvv-1+k,lan)= syn(i,kvv-1+k,lan)*rsltbm*ppindex(i,lan)
#endif
#ifdef NONHYD
#ifdef INTDTB
            syn(i,ktt-1+k,lan)= syn(i,ktt-1+k,lan)*rsltbh*ppindex(i,lan)
#endif
            snn(i,knpp-1+k,lan)=                                       &
     &       snn(i,knpp-1+k,lan)*rsltbh*rsltn*ppindex(i,lan)
            snn(i,kntt-1+k,lan)=                                       &
     &       snn(i,kntt-1+k,lan)*rsltbh*rsltn*ppindex(i,lan)
#else
#ifdef UPDAMP
            syn(i,ktt-1+k,lan)= syn(i,ktt-1+k,lan)*rsltbh*rsltn*ppindex(i,lan)
#else
            syn(i,ktt-1+k,lan)= syn(i,ktt-1+k,lan)*rsltbh*ppindex(i,lan)
#endif
#endif
          enddo
        enddo
      endif
!
#ifdef A
      do k=1,levr
      do i=1,lonlens
      j=i+latj
#endif
#ifdef T
      wt(j,k,13)=syn(i,ktt-1+k,lan)-wt(j,k,13)
#endif
#ifdef U
      wu(j,k, 9)=syn(i,kuu-1+k,lan)-wu(j,k, 9)
#endif
#ifdef V
      wv(j,k, 9)=syn(i,kvv-1+k,lan)-wv(j,k, 9)
#endif
#ifdef T
      wt(j,k,14)=-wt(j,k,14)+syn(i,ktt-1+k,lan)
#endif
#ifdef U
      wu(j,k,10)=-wu(j,k,10)+syn(i,kuu-1+k,lan)
#endif
#ifdef V
      wv(j,k,10)=-wv(j,k,10)+syn(i,kvv-1+k,lan)
#endif
#ifdef A
      enddo
      enddo
      do k=1,levh
      do i=1,lonlens
      j=i+latj
#endif
#ifdef Q
      wq(j,k, 8)=srq(i,k,lan)-wq(j,k, 8)
      wq(j,k, 9)=-wq(j,k, 9)+srq(i,k,lan)
#endif
#ifdef A
      enddo
      enddo
#endif
#ifdef P
      do i=1,lonlens
      j=i+latj
      wp(j,1, 4)=syn(i,kpp,lan)-wp(j,1, 4)
      wp(j,1, 5)=-wp(j,1, 5)+syn(i,kpp,lan)
      enddo
#endif
!
2080  continue

!...........................
!cj
#ifdef RKN
!
!  ken point - div and vort
!
      if(ilshrt.lt.1) then
        if(isave.ne.0.and.npoint.gt.0) then
         lones=lonlen(mype)
         lon1s=lonstr(mype)
         lon2s=lon1s-1+lones
         latt1=latts+lat1s
         latt2=latts+lat2s
!
#ifdef OPENMP
!$omp parallel do private(igpt,jlat,jlat1,ilon,k)
#endif
         do igpt=1,npoint
         if(jgpd(igpt).ge.latt1.and.jgpd(igpt).le.latt2) then
         if((igpd(igpt).ge.lon1s.and.igpd(igpt).le.lon2s) ) then
!     &      (igpd(igpt).ge.lon1n.and.igpd(igpt).le.lon2n)) then
           jlat=jgpd(igpt)-latt1+1
           jlat1=jlat
           if(igpd(igpt).ge.lon1s.and.igpd(igpt).le.lon2s) then
             ilon=igpd(igpt)-lon1s+1
           endif

           do k=1,levr
           svdata(k+rslvark+8*levr,igpt,itnum) =                         &  
     &        xm2(ilon,jlat1)*                                           &
     &       (udx(ilon,k,jlat)+vdy(ilon,k,jlat))
           svdata(k+rslvark+9*levr,igpt,itnum) =                         &
     &        xm2(ilon,jlat1)*                                           &
     &       (vdx(ilon,k,jlat)-udy(ilon,k,jlat))
           enddo
         endif
         endif
         enddo
        endif
      endif
#endif
!
#ifdef CNT_1
        call cnttime_on("one_loop_spectal")
#endif
!
      call rmpnx2x(syn,igrd1p,lotsyn,synk,igrd1,lotsynk,                 &
     &             jgrd1p,kpp,kppp,1)
      call rmpnx2nk(syn,igrd1p,lotsyn,synk,igrd1,lotsynk,                &
     &             jgrd1p,levr,levrp,kvv,kvvp,1)
      call rmpnx2nk(syn,igrd1p,lotsyn,synk,igrd1,lotsynk,                &
     &             jgrd1p,levr,levrp,kuu,kuup,1)
      call rmpnx2nk(syn,igrd1p,lotsyn,synk,igrd1,lotsynk,                &
     &             jgrd1p,levr,levrp,ktt,kttp,1)
#ifdef NONHYD
      call rmpnx2nk(snn,igrd1p,lotsnn,snnyk,igrd1,lotsnnk,               &
     &             jgrd1p,levr,levrp,knpp,knppp,3)
      call rmpnx2x(snn,igrd1p,lotsnn,snnyk,igrd1,lotsnnk,                &
     &             jgrd1p,knww+levr,knwwp+levrp,1)
#endif
!
      call ffancsx1(synk,jgrd1p,kvvp,kuup,lotsynk,levrp)
      call ffacosx1(synk,jgrd1p,kppp,lotsynk,levrp+1)
#ifdef NONHYD
      call ffacosx1(snnyk,jgrd1p,knppp,lotsnnk, levrp*3+1 )
#endif
!
      call rmpny2nl(synk,igrd1,jgrd1p,syfk,igrd1p,                       &
     &       jgrd1,lotsynk,1,lotsynk)
#ifdef NONHYD
      call rmpny2nl(snnyk,igrd1,jgrd1p,snnxk,igrd1p,                     &
     &       jgrd1,lotsnnk,1,lotsnnk)
#endif
!
!.1. uu,vv
      call ffancsy1(syfk,xa,wa,kuup,kvvp,igrd1p,iwav1p,lotsynk,levrp)
!.2. z,y,rt
      call ffacosy1(syfk,za,kppp,igrd1p,iwav1p,lotsynk,levrp+1)
#ifdef NONHYD
!.3. 
      call ffacosy1(snnxk,pa,knppp,igrd1p,iwav1p,lotsnnk,levrp*3+1)
#endif
!..................................................................
!
      call rmpnk2nn(ya,llwavp,levrp,y,lnwavp,levr,1)
      call rmpn2nn(za,llwavp,z,lnwavp,1)
      call rmpnk2nn(xa,llwavp,levrp,x,lnwavp,levr,1)
      call rmpnk2nn(wa,llwavp,levrp,w,lnwavp,levr,1)
#ifdef NONHYD
      call rmpnk2nn(pa,llwavp,levrp,p,lnwavp,levr,3)
      call rmpn2nn(oa(1,levrp1p),llwavp,o(1,levrp1),lnwavp,1)
#endif
!
#ifdef CNT_1
        call cnttime_off("one_loop_spectal")
#endif
!
#ifdef T
      call cctogd(y,tmpbgt,levr)
      do k=1,levr
        do i=1,lngrdp
          wt(i,k,14)=tmpbgt(i,k)-wt(i,k,14)
        enddo
      enddo
#endif
!
#ifdef Q
      call cctogd(rt,tmpbgt,levh)
      do k=1,levh
        do i=1,lngrdp
          wq(i,k, 9)=tmpbgt(i,k)-wq(i,k, 9)
        enddo
      enddo
#endif
!
#ifdef U
      call sctogd(x,tmpbgt,levr)
      do k=1,levr
        do i=1,lngrdp
          wu(i,k,10)=tmpbgt(i,k)-wu(i,k,10)
        enddo
      enddo
#endif
!
#ifdef V
      call cstogd(w,tmpbgt,levr)
      do k=1,levr
        do i=1,lngrdp
          wv(i,k,10)=tmpbgt(i,k)-wv(i,k,10)
        enddo
      enddo
#endif
!
#ifdef P
      call cctogd(z,tmpbgt(1,1),1)
      do i=1,lngrdp
        wp(i,1, 5)=tmpbgt(i,1)-wp(i,1, 5)
      enddo
#endif
!
#ifndef RAYLEIFR      
!yj: move to above if do rayleifr --->
#ifdef OPENMP
!$omp parallel do private(k,lat)
#endif
!
      do k=1,levr
      spdmax(k) = 0.0
       do lat=1,latlen(mype)
        spdmax(k)=max(spdmax(k),spdlat(k,lat))
       enddo
      enddo
      call rmpgetspd(levr,spdmax)
!yj: move to above if do rayleifr ---<
#endif
!
! update ww at top for n time step (see nfidi)
!
#ifdef NONHYD
      do j=1,lnwavp
        on(j,levr+1)=onm(j,levr+1)+deltim*o(j,levr+1)
      enddo
#endif

!      deallocate(indfull,inds,tmpgrd)
!      deallocate(indfullv,tmpgrdv)
! yj: deallocate local array >------------------------------------------>
       deallocate( dot)
       deallocate(syn,srq,synr,synm,syndd,bgd,syfk,synk,                 &
     &            synmk,rindex,synt,pindex,ppindex,synuk,synvk,anl)
#ifndef NDSLDYN 
       deallocate( utmp, vtmp)
#endif
#ifdef NONHYD
       deallocate( snnxk,snnyk,annxk,annyk,snn,snnm,ann)
#endif
#ifdef NDSLDYN
#ifdef NONHYD
       deallocate(w2nd)
#endif
#endif
! yj: deallocate local array <------------------------------------------<

!......................................................
      if(me.eq.0) print 100,(spdmax(k),k=1,levr)
!orig        
!100   format(' rsm spdmx(01:10)=',10f5.0,/' rsm spdmx(11:20)=',10f5.0,   &  
!     &     /' rsm spdmx(21:30)=',10f5.0,/' rsm spdmx(31:40)=',10f5.0,    &
!     &     /' rsm spdmx(41:50)=',10f5.0,/' rsm spdmx(51:60)=',10f5.0,    &
!     &     /' rsm spdmx(61:70)=',10f5.0,/' rsm spdmx(71:80)=',10f5.0,    &
!     &     /' rsm spdmx(81:90)=',10f5.0,/' rsm spdmx(91:00)=',10f5.0)
 100   format(' rsm spdmx(01:10)=',10f10.5,/' rsm spdmx(11:20)=',10f10.5,   &  
      &     /' rsm spdmx(21:30)=',10f10.5,/' rsm spdmx(31:40)=',10f10.5,    &
      &     /' rsm spdmx(41:50)=',10f10.5,/' rsm spdmx(51:60)=',10f10.5,    &
      &     /' rsm spdmx(61:70)=',10f10.5,/' rsm spdmx(71:80)=',10f10.5,    &
      &     /' rsm spdmx(81:90)=',10f10.5,/' rsm spdmx(91:00)=',10f10.5)
!
#endif
!

      return
      end

