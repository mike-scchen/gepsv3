!#define NEW_WSM3
!#define NEW_WSM6
#include <define.h>
      subroutine rgetcon_mrg
!....
      use machine , only : kind_io4
      use funcphys
      use module_progtm
      use physcons, rhoair0 => con_rhoair0 , rhoh2o => con_rhoh2o ,        &
                    rhosnow => con_rhosnow , cliq => con_cliq ,            &
                    cvap => con_cvap  
#ifdef NEW_WSM3
      use rmodule_mp_wsm3,        only: wsm3init
#else
      use rmodule_mp_wsm3
#endif
      use module_mp_wsm5
#ifdef NEW_WSM6
      use rmodule_mp_wsm6,        only: wsm6init
#else
      use rmodule_mp_wsm6
#endif
! for GFDL MP v1
      use module_mp_gfdl_v16,    only : gfdl_cloud_microphys_init
! for GFDL MP v2
      use rmodule_mp_gfdl_v2,     only : gfdlv2_init => gfdl_cld_mp_init
! for GFDL MP v3
      use module_mp_gfdl_v3,     only : gfdlv3_init => gfdl_cld_mp_init
! for GCE 3ICE
      use module_mp_gsfcgce_3ice_nuwrf,     only : consat_s 

!ccj
!cccc  jun change
! add isfcmrg in namelist and ilist and rsfcsec in namelist and rlist
!cccc  jun change end
#include <paramodel.h>
!
      common /rcomfgi/ n0,n1,n2
#include <rscomio.h>
#ifdef CRA
#include <rscomcpu.h>
#endif
#include <rscomloc.h>
!
#include <rscomgpd.h>
!cc
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
      dimension ilist(100),rlist(100),clist(100),glist(100)
#endif
      namelist/namrsm/                                                   &  
     &                icen,igen,icen2,ienst,iensi,runid,usrid,           &
     &                deltime,filta,seminpa,ncpus,                       &
     &                fcstsec,prntsec,rswrsec,rlwrsec,                   &
     &                basesec,rdfisec,rkensec,                           &
     &                rlxhsec,rlxmsec,rlxzsec,difhsec,difmsec,rsfcsec,   &
     &                restrhr,readclimhr,                                &
     &                isemimp,iimprlx,idmpjet,imdlphy,ioutnhr,           &
     &                isfcmrg,irclim,cplf_timcom
#ifdef RKN
      namelist/namrkn/                                                   &
     &                npoint,idstat,gpdnam,gpdlat,gpdlon,gpdgzs,         & !binbin add: idstat
     &                fhken,nsken,isshrt,ilshrt
#endif
      namelist/namloc/                                                  &
     & rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd &
     &,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd &
     &,clat1,clat2,clon1,clon2
      real(kind=kind_io4) fhours
      logical, parameter   ::  allowed_to_read = .true.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!shian
      cplf_timcom=7200.
!  default namelist
      filta=0.92
      semimpa=0.8
      deltime=0
      restrhr=0.
      fcstsec=10800.
      prntsec=10800.
      basesec=10800.
      readclimhr=24.
      rlwrsec=3600.
      rswrsec=3600.
      rdfisec=0.
      rkensec=60.
      difhsec=600.
      difmsec=600.
      rlxhsec=1800.
      rlxmsec=1800.
      rlxzsec=1800.
      rsfcsec=86400.
      isemimp=1
      iimprlx=1
      idmpjet=1
      imdlphy=1
      ioutnhr=1
      isfcmrg=0
      irclim=1
      dtphys=0
      icen=7
      igen=rgen
      icen2=0
      ienst=0
      iensi=0
      runid=0
      usrid=0
      ncpus=1	! for opnemp, default is 1
      ndigyr=4
#ifdef RKN
      npoint=0
      nsken=0
      fhken=1
      gpdlat=0
      gpdlon=0
      isshrt=0
      ilshrt=0
#endif
      call gncpus(ncpus)
      ncpus1=ncpus+1
#ifdef MP
      if( mype.eq.master ) then
#endif
      print *,' ncpus = ',ncpus
#ifdef MP
      endif
#endif
!.............................................
! read namelist of forecast parameter
!
#ifdef MP
      if( mype.eq.master ) then
#endif
      read (5,namrsm)
      write(6,namrsm)
#ifdef MP
!--------
      ilist(1)=icen
      ilist(2)=igen
      ilist(3)=icen2
      ilist(4)=ienst
      ilist(5)=iensi
      ilist(6)=ncpus
      ilist(7)=isemimp
      ilist(8)=iimprlx
      ilist(9)=idmpjet
      ilist(10)=imdlphy
      ilist(11)=ioutnhr
      ilist(12)=isfcmrg
      ilist(13)=irclim
!
      rlist(1)=runid
      rlist(2)=usrid
      rlist(3)=deltime
      rlist(4)=filta
      rlist(5)=seminpa
      rlist(6)=fcstsec
      rlist(7)=prntsec
      rlist(8)=rswrsec
      rlist(9)=rlwrsec
      rlist(10)=basesec
      rlist(11)=rdfisec
      rlist(12)=rkensec
      rlist(13)=rlxhsec
      rlist(14)=rlxmsec
      rlist(15)=difhsec
      rlist(16)=difmsec
      rlist(17)=rsfcsec
      rlist(18)=restrhr
      rlist(19)=readclimhr
!shian
      rlist(20)=cplf_timcom
      endif
!
      call rmpbcasti(ilist,13)
      call rmpbcastr(rlist,20)
!
      icen   =ilist(1)
      igen   =ilist(2)
      icen2  =ilist(3)
      ienst  =ilist(4)
      iensi  =ilist(5)
      ncpus  =ilist(6)
      isemimp=ilist(7)
      iimprlx=ilist(8)
      idmpjet=ilist(9)
      imdlphy=ilist(10)
      ioutnhr=ilist(11)
      isfcmrg=ilist(12)
      irclim=ilist(13)
!
      runid  =rlist(1)
      usrid  =rlist(2)
      deltime=rlist(3)
      filta  =rlist(4)
      seminpa=rlist(5)
      fcstsec=rlist(6)
      prntsec=rlist(7)
      rswrsec=rlist(8)
      rlwrsec=rlist(9)
      basesec=rlist(10)
      rdfisec=rlist(11)
      rkensec=rlist(12)
      rlxhsec=rlist(13)
      rlxmsec=rlist(14)
      difhsec=rlist(15)
      difmsec=rlist(16)
      rsfcsec=rlist(17)
      restrhr=rlist(18)
      readclimhr=rlist(19)
!shian
      cplf_timcom=rlist(20)
!------
#endif
!
! read namelist of ken point station location
! gwv  everybody reads the namelists
!
#ifdef RKN
!j      read (5,namrkn)   !standard input can't read larger than 32k file
!jw      open (500,file='station.parm',status='old')
!jw      read (500,namrkn)   !so, change to 500 to read 1300 stations
!jw      do i=1,npoint
!jw        idstat(i)=i
!jw      enddo
! -- done by Jun Du and Binbin on Jan. 7, 2005
      open (500,file='station.parm',status='old')
      do i = 1, nptken
        read (500,501) gpdlat(i),gpdlon(i),gpdnam(i),ztmp
!hmhj   gpdgzs(i) = float(ztmp)
        gpdgzs(i) = ztmp
        gpdlat(i) = 100.0 + gpdlat(i)
        gpdlon(i) = -1.0*gpdlon(i)
#ifdef MP
      if( mype.eq.master ) then
#endif
          print *,gpdlat(i),gpdlon(i),gpdnam(i),gpdgzs(i)
#ifdef MP
      endif
#endif
501     format(6x,f6.2,1x,f7.2,2x,a4,31x,i4)
      end do
#ifdef MP
!jw
      call rmpsynall
#endif
      close(500)

      npoint = nptken
      do i=1,npoint
        idstat(i)=i
      enddo

!
#ifdef MP
      if( mype.eq.master ) then
#endif
        do i=1,npoint
          print *,idstat(i),gpdnam(i),gpdlat(i),gpdlon(i),gpdgzs(i)
        enddo
#ifdef MP
      endif
#endif
#endif

! read namelist of model location
!
#ifdef MP
      if( mype.eq.master ) then
#endif
      read (5,namloc)
      write(6,namloc)
#ifdef MP
!--------
      glist(1)=rproj
      glist(2)=rtruth
      glist(3)=rorient
      glist(4)=rdelx
      glist(5)=rdely
      glist(6)=rcenlat
      glist(7)=rcenlon
      glist(8)=rlftgrd
      glist(9)=rbtmgrd
      clist(1)=cproj
      clist(2)=ctruth
      clist(3)=corient
      clist(4)=cdelx
      clist(5)=cdely
      clist(6)=ccenlat
      clist(7)=ccenlon
      clist(8)=clftgrd
      clist(9)=cbtmgrd
!
      endif
      call rmpbcastr(glist,9)
      call rmpbcastr(clist,9)
!
      rproj  =glist(1)
      rtruth =glist(2)
      rorient=glist(3)
      rdelx  =glist(4)
      rdely  =glist(5)
      rcenlat=glist(6)
      rcenlon=glist(7)
      rlftgrd=glist(8)
      rbtmgrd=glist(9)
      cproj  =clist(1)
      ctruth =clist(2)
      corient=clist(3)
      cdelx  =clist(4)
      cdely  =clist(5)
      ccenlat=clist(6)
      ccenlon=clist(7)
      clftgrd=clist(8)
      cbtmgrd=clist(9)
!------
#endif
!.............................................
!
! set constant and alter namelist if necessary
      filtb=(1.-filta)*0.5
      deltim=deltime
      if(deltim.le.0.0) deltim=rdelx/10000.*60.
      if(difmsec.le.0.0) difmsec=3*deltim
      if(difhsec.le.0.0) difhsec=6*deltim
      nsseg = nint(fcstsec/deltim)
      nsres = nint(restrhr*3600./deltim)
      nsout = nint(prntsec/deltim)
      nsbas = nint(basesec/deltim)
      nsclm = nint(readclimhr*3600./deltim)
      nsswr = nint(rswrsec/deltim)
      nslwr = nint(rlwrsec/deltim)
      nsdfi = nint(rdfisec/deltim)
!     fhseg = nsseg*deltim/3600.
      fhseg = (nsseg-nsres)*deltim/3600.
      fhout = nsout*deltim/3600.
      fhbas = nsbas*deltim/3600.
      fhswr = nsswr*deltim/3600.
      fhlwr = nslwr*deltim/3600.
      fhdfi = nsdfi*deltim/3600.
      lfnhr = ioutnhr.eq.1
      lsfcmrg = isfcmrg.eq.1
      lrestart = nsres.ne.0
      lsphy = imdlphy.eq.1
      lrclm = irclim.eq.1
#ifdef CPL_TIMCOM
!      kdt_cpl = nsout
!      kdt_cpl = kdt_cpl/3
      kdt_cpl = nint(cplf_timcom/deltim)
      odt_cpl = 1.0/(kdt_cpl*deltim)
      if(mype.eq.master) then
        print *,'[RSM] kdt_cpl=',kdt_cpl,', dt_cpl=',kdt_cpl*deltim
      end if 
#endif
#ifdef RKN
      if(npoint.gt.0.and.npoint.le.nptken) then
        if(nsken.le.0) then
           fhken=max(fhken,deltim/3600)
        else
          fhken=nsken*deltim/3600
        endif
        fhken=fhout/max(nint(fhout/fhken),1)
        nsken=nint(fhken*3600/deltim)
        if(nsout.gt.nsken*nstken) npoint=0
      else
        npoint=0
      endif
      nsken=max(nsken,1)
#endif
!
!  define unit numbers
!
!  input
      n0    = 18  !for read rsm_idate
      n1    = 11  !for read rsm_data
      n2    = 12  !for read rsm_data
      if(nsseg.le.0 .or. nsout.le.0 .or. nsbas.le.0 ) then 
#ifdef MP
      if( mype.eq.master ) then
#endif
        print*,'need to read fhour1 and fhour2'
!       rewind n1
!       rewind n2
!       read(n1)
!       read(n2)
!       read(n1) fhours
!       fhour1=fhours
!       read(n2) fhours
!       fhour2=fhours
!       rewind n1
!       rewind n2
#ifdef MP
      endif
      call rmpbcastr(fhour1,1)
      call rmpbcastr(fhour2,1)
#endif
        basesec=(fhour2-fhour1)*3600.
        nsbas = nint(basesec/deltim)
        nsseg = nsbas
        nsout = nsbas
        fhbas = nsbas*deltim/3600.
        fhseg = fhbas
        fhout = fhbas
      endif
      nges  = 13  !yj not used
      nflip = 14  !yj not used
      nco2  = 15  !yj not used (link empty file)
      nradr = 21  !yj not used
      nradf = 22  !yj not used
!     nnmod = 23  !yj used rmtnslm
      nrslm = 23
      nrorg = 24
      nmtnv = 42  !yj used rmtnvar
      nrsmi0= 33  !yj used input rsm_idate
      nrsmi1= 30  !yj used input r_sigi (r_sigfxx)
!      idate from GFS
      nrsmi2= 31  !yj used input r_sigitdt (restart)
!      data(u,v,t....) from GFS
      nrflip= 32  !yj used input r_sfci (r_sfcfxx)
      nb    = 34  !yj used input base field for restart
      nctune= 43  !yj used TUNE1
      nozon = 45  !yj used O3CLIM #48 conflict with gfs
!  output
      n3    = 51
      n4    = 52
      nflop = 53
      nflps = 54
      nsigi = 55
      nsigs = 56
      nznli = 61
      nsfci = 62
      nsfcf = 63
      nznlf = 64
      ndgi  = 65
      ndgf  = 66
      ngpken= 67
      nsfcs = 68
      nznls = 69
      nrsmo1= 70 ! used output r_sig.f
      nrsmo2= 71 ! used output r_sigftdt 
      nrflop= 72 ! used output r_sfc.f
      nrflxi= 74 ! used output r_flx.f00
      nrsmop= 75 ! for intermediate sigma file output
      nrsfcp= 76 ! for intermediate surface file output
      nrkenp= 77 ! for intermediate ken point output
      nrinit= 78
      nrflxp= 79 ! for intermediate flux file output
      nrsms1= 80 ! used output r_sig.f
      nrsms2= 81 ! used output r_sigftdt
      nrflop8= 82 ! used output r_sfc.f
#ifdef A
      nrbgt = 85
!  work files 
      nrbgt1 = 95
      nrbgt2 = 96
#endif
      nr2dda = 98
!c
      call cmpind
!c      call gpvs
!c      call gtdp
!c      call gthe
!c      call gtma
      call gfuncphys
#ifdef G2R
      call gpln2i
      call epsilo(epsi,jcap)
      call ggozrm(epsi)
      call gft_lonf
#endif
      call grddf
      call grdkt
!
      if (allowed_to_read) then
#ifdef NEW_WSM3
        if (nmmiph.eq.4) call wsm3init()
#else
        if (nmmiph.eq.4) call wsm3init(rhoair0, rhoh2o, rhosnow, cliq, cvap,  &
                           allowed_to_read)
#endif
        if (nmmiph.eq.5) call wsm5init(rhoair0, rhoh2o, rhosnow, cliq, cvap,  &
                           allowed_to_read)
#ifdef NEW_WSM6
        if (nmmiph.eq.6) call wsm6init()
#else
        if (nmmiph.eq.6) call wsm6init(rhoair0, rhoh2o, rhosnow, cliq, cvap,  &
                          allowed_to_read)
#endif
        if (nmmiph.eq.11) call gfdl_cloud_microphys_init()
        if (nmmiph.eq.12) call gfdlv2_init()
        if (nmmiph.eq.13) call gfdlv3_init()
      endif
!yj
!yj      print*,'check dimension in rgetcon_mrg.F'
!yj      print*,'rdelx=',rdelx
!yj      print*,'cdelx=',cdelx
!
#ifdef NDSL
      call rndslfv_init
      if( mype.eq.master )  print*,'after call rndslfv_init'
#endif      
      return
      end
