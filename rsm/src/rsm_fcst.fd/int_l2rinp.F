#include <inp_defopt.h>
      subroutine int_l2rinp(nrslm)
!......................................................
!$$$  subprogram documentation block
! abstract:  do regional lon-lat grid to regional grid interpolation
!
! program history log:
!
! yunfei zhang 2009
!
!$$$
!....
      use machine , only : kind_io4
      use physcons, pi => con_pi
!yj add module start -----------------------      
      use module_rinpinit, only :  &      
     &  siin,slin,delin,ciin,clin,rpiin &
     & ,siin3d,slin3d,terr &
!yj     & ,idate &
!comrinp2     
     & ,q,te,uu,vv,rq,gz &
!yj!#ifdef NONHYD
!yj!     &      ,pp,tt,ww &
!yj!#endif
     &      ,flat,flon,fm2,fm2x,fm2y &
!comrinpsfc
     & ,tsea,smc,sheleg,stc,tg3clm,zorclm &
     & ,cv,cvb,cvt &
     & ,albclm,slmsk,vegclm,canopy, f10m  &
     & ,vetclm,sotclm,alfclm              &
     & ,vmnclm,vmxclm,slpclm,absclm       &
     & ,uustar,ffmm,ffhh
!yj add module   end -----------------------      
!yj add module gfsdata start----------------
      use module_gfsdata,only:             &
     & idate,fhour,       &
     & te_gfs=>temp, sp_gfs=>spfh,         &
     & clw_gfs=>clwr, rain_gfs=>rain,      &
     & qice_gfs=>qice, snow_gfs=>snow,      &
     & grpl_gfs=>grpl, o3_gfs=>ozon,        &
     & ge_gfs=>geop, u_gfs=>u, v_gfs=>v,   &
     & terr_gfs=>terr,                     &
     & tg_gfs=>tg, smc_gfs=>smc,           &
     & snr_gfs=>snr, stc_gfs=>stc,         &
     & cice_gfs=>cice,slmk_gfs=>slmsk
!yj add module gfsdata   end----------------
!ccj
#include <inpint_paramodel.h>
#include <inp_rscomloc.h>
!cc
!...
      parameter (nfld=8)
!gfs: t spfh clwmr ozone z u v terrain     
!rsm: z t u v spfh o3 clwmr terrain
!yjgfs      parameter (nsfcfld=13)
      parameter (nsfcfld=8)
!gfs: tg 0.1soilmois 0.4soilmois snr 0.1soiltemp 0.4soiltemp lsmask cice      
      parameter (imax=cigrd1, jmax=cjgrd1)
      character*8 label(4)
      real plevels(levs)
      data plevels/1000,975,950,925,900,875,850,825,800    &
                  ,775,750,725,700,675,650,625,600,575     &
                  ,550,525,500,475,450,425,400,375,350     &
                  ,325,300,275,250,225,200,175,150,125     &
                  ,100,70,50,30,20,10,7,5,3,2,1         /
!
      real(kind=kind_io4) fhours                                        
!     &         ,rgs(lngrd),  rgls(lngrd,lsoil)                          &
!     &         ,rg2s(lngrd,2),  rg4s(lngrd,4)
!yj
      integer :: idatenew(4)
!      integer :: idate(4)
!      real :: sl(levs),si(levs+1)
!      real :: flat(lngrd),flon(lngrd) 
!      real :: fm2(lngrd),fm2x(lngrd),fm2y(lngrd)
!      real :: gz(lngrd), q(lngrd) 
!      real :: te(lngrd,levs), rq(lngrd,levs,ntotal)
!      real :: uu(lngrd,levs), vv(lngrd,levs) 
!yj
!
      real(kind=kind_io4), allocatable :: rgs(:)
      real, allocatable :: slg(:)
      real, allocatable ::  cg(:)  
!      real, allocatable ::  cgg(:,:)
!      real, allocatable ::  cg(:)  , cgl(:,:)
!      real, allocatable :: cg2(:,:), cg4(:,:)
      real, allocatable :: blat(:),blon(:)
!
      real, allocatable :: cgsea(:),cgland(:),rr(:)
      integer, allocatable :: iwork(:)
!      integer, allocatable :: iwork(:),idate(:)
!      real, allocatable :: sl(:),si(:)
!      real, allocatable :: flat(:),flon(:)
      real, allocatable :: cm2(:)
!      real, allocatable :: fm2(:),fm2x(:),fm2y(:)
      real, allocatable :: dflat(:),dflon(:)
!
!      real, allocatable :: gz(:), q(:)
!      real, allocatable :: te(:,:), rq(:,:,:)
!      real, allocatable :: uu(:,:), vv(:,:)
!#ifdef NONHYD
!      real, allocatable :: pp(:,:), tt(:,:)
!      real, allocatable :: ww(:,:)
!#endif
!      real, allocatable :: rain(:),cvrn(:)
      real, allocatable :: rcsln(:),rsnln(:)
!....
!      real, allocatable :: tsea  (:), smc   (:,:)
!      real, allocatable :: sheleg(:), stc   (:,:)
!      real, allocatable :: tg3   (:), vfrac (:)
!      real, allocatable :: zorl  (:)
!      real, allocatable :: zorl  (:), slmsk (:)
!      real, allocatable :: cv    (:), cvb   (:)
!      real, allocatable :: cvt   (:)
!      real, allocatable :: f10m  (:), canopy(:)
!      real, allocatable :: vtype (:)
!      real, allocatable :: stype (:)
!      real, allocatable :: stype (:), uustar(:)
!      real, allocatable :: ffmm  (:),   ffhh(:)
!      real, allocatable :: alvsf (:),  alvwf(:)
!      real, allocatable :: alnsf (:),  alnwf(:)
!      real, allocatable :: facsf (:),  facwf(:)
      real, allocatable :: slmk(:),simk(:)
!      real, allocatable :: simk(:)
!      real, allocatable :: tg3clm(:),zorclm(:)
!      real, allocatable :: albclm(:,:),vegclm(:)
!      real, allocatable :: vetclm(:),sotclm(:),alfclm(:,:)
!yj      real, allocatable :: cfld(:,:,:)
!yj      real, allocatable :: csfcfld(:,:)
!yj      real, allocatable :: spfh(:,:)
!
      character*500          fntg3c,fnzorc,fnalbc,fnvegc,fnvetc,fnsotc, &
     &                 fnmskh,fnvmnc,fnvmxc,fnslpc,fnabsc
      namelist/namclim/iclim,fntg3c,fnzorc,fnalbc,fnvegc,fnvetc,fnsotc, &
     &                 fnmskh,fnvmnc,fnvmxc,fnslpc,fnabsc
!yj      common /hydrid/ siin3d(lngrd,levs+1),slin3d(lngrd,levs)
!yj      common /llgz/ terr(lngrd)
      data fntg3c/'global_tg3clim.2.6x1.5.grb'/
      data fnmskh/'global_slmask.t126.grb'/
      data fnalbc/'global_albedo4.1x1.grb'/
      data fnsotc/'global_soiltype.1x1.grb'/
      data fnvegc/'global_vegfrac.1x1.grb'/
      data fnvetc/'global_vegtype.1x1.grb'/
      data fnzorc/'global_zorclim.1x1.grb'/
      data fnvmnc/'global_shdmin.0.144x0.144.grb'/
      data fnvmxc/'global_shdmax.0.144x0.144.grb'/
      data fnslpc/'global_slope.1x1.grb'/
      data fnabsc/'global_snoalb.1x1.grb'/
!yj      character*500 fncwbpre, fncwbsfc
!yj      data fncwbpre/'input.prel'/
!yj      data fncwbsfc/'input.sfcl'/
!
      print *,' enter int_l2rinp '
!
      allocate ( rgs(lngrd) )
      allocate ( slg(clngrd) )
      allocate (  cg(clngrd) )
!      allocate (  cg(clngrd)  , cgl(clngrd,lsoil) )
!      allocate ( cg2(clngrd,2), cg4(clngrd,4) )
      allocate ( blat(clngrd) ,blon(clngrd) )
!
      allocate ( cgsea(clngrd),cgland(clngrd) )
      allocate ( iwork(clngrd) )
!      allocate ( idate(4),sl(levs),si(levs+1) )
!      allocate ( flat(lngrd),flon(lngrd) )
      allocate ( cm2(lngrd) )
!      allocate ( fm2(lngrd),fm2x(lngrd),fm2y(lngrd) )
      allocate ( dflat(lngrd),dflon(lngrd) )
!
!      allocate ( gz(lngrd), q(lngrd) )
!      allocate ( te(lngrd,levs), rq(lngrd,levs,ntotal) )
!      allocate ( uu(lngrd,levs), vv(lngrd,levs) )
!      allocate ( rain(lngrd),cvrn(lngrd) )
      allocate ( rcsln(lngrd),rsnln(lngrd) )
!....
      allocate ( rr    (lngrd) )
!      allocate ( tsea  (lngrd), smc   (lngrd,lsoil) )
!      allocate ( sheleg(lngrd), stc   (lngrd,lsoil) )
!      allocate ( tg3   (lngrd), vfrac (lngrd) )
!      allocate ( zorl  (lngrd) )
!      allocate ( zorl  (lngrd), slmsk (lngrd) )
!      allocate ( cv    (lngrd), cvb   (lngrd) )
!      allocate ( cvt   (lngrd) )
!      allocate ( f10m  (lngrd), canopy(lngrd) )
!      allocate ( vtype (lngrd) )
!      allocate ( stype (lngrd) )
!      allocate ( stype (lngrd), uustar(lngrd) )
!      allocate ( ffmm  (lngrd),   ffhh(lngrd) )
!      allocate ( alvsf (lngrd),  alvwf(lngrd) )
!      allocate ( alnsf (lngrd),  alnwf(lngrd) )
!      allocate ( facsf (lngrd),  facwf(lngrd) )
       allocate ( slmk(lngrd),simk(lngrd) )
!      allocate ( simk(lngrd) )
!      allocate ( tg3clm(lngrd),zorclm(lngrd) )
!      allocate ( albclm(lngrd,4),vegclm(lngrd) )
!      allocate ( vetclm(lngrd),sotclm(lngrd),alfclm(lngrd,2) )
!yj      allocate ( cfld(clngrd,levs,nfld) )
!yj      allocate ( csfcfld(clngrd,nsfcfld) )
!yj      allocate ( spfh(clngrd,levs) )
!
!.................................................................
! define regional grid
      call int_setlatlon(flat,flon,delx,dely,dlamda0)
      call maxmin(flat,lngrd,1,1,1,' flat ')
      call maxmin(flon,lngrd,1,1,1,' flon ')
!
!!!!!gfs!!!!!      
!yj      if (cwbgfs.eq.1) then
! yj add binary output from gfs
!yj        print*,'replace by cwbgfs --- start reading rsm_atm.out'
!yj         close(npgb)
!yj         open(npgb,file=fncwbpre,status='old', &
!yj             form='unformatted',iostat=ios)
!yj         if( ios.ne.0) then
!yj           print *,' error open rsm_atm.out'
!yj           stop
!yj         endif 
!yj         do n=1,nfld    ! terr is read from rsm-terr.bin
!yj           do k=1,levs
!yj             if (n.eq.nfld) then
!yj               read(npgb,iostat=eof) cfld(:,1,n)
!yj               if (eof.ne.0) print*,'error reading file'
!yj               print*,'read the end of the atm file,eof=',eof
!yj               print*,'replace with terr_grib2 n=nfld=',n
!yj                 open(77,file='terr_grib2.bin',form='unformatted')
!yj                 read(77) cfld(:,1,n)
!yj                 close(77)
!yj               print*,'checking (n,k)=(8,1)?',n,k
!yj               exit
!yj             else  
!yj               print*,'k=',k
!yj               read(npgb,iostat=eof) cfld(:,levs-k+1,n)
!yj               if (eof.ne.0) print*,'error reading file'
!yj             endif  
!yj           enddo
!yj         enddo
!-----------------------------------------------------------
!the order of var in cfld is :
!1.temp 2.spfh 3.cloud water 4.o3 5.geop 6.u 7.v 8.terr
!-----------------------------------------------------------
! yj mpmd recv data and put in cfld!!!! wait for MPMD !!!!!

!yj         open(99,file='date.input.log',status='old',iostat=ios)
!yj         if(ios.ne.0) then
!yj           print*,' error open date.input.log'
!yj           stop
!yj         endif
!yj         read(99,*) idate(1)
!yj         read(99,*) idate(2)
!yj         read(99,*) idate(3)
!yj         read(99,*) idate(4)
!yj         read(99,*) fhour
!yj         close(99)
!yj      endif  !!!end if (cwbgfs.eq.1)
!yj        close(npgb)
!!!!!gfs!!!!!      
      
      call int_setcgrid(blat,blon)

! .....................................................
! define coarse regional to regional grid interpolation
      call int_r2rini(blat,blon,flat,flon)
      call int_rtorini(blat,blon,flat,flon)
! ................................................
! define regional map factors
      rad = pi / 180.0
      hfpi= 0.5 * pi
      truth = rtruth
      nproj = rproj
      if( nproj.eq.1 .or. nproj.eq.-1 ) then
        orient = rorient
        dlamda0 = orient * rad + hfpi
        a2 = ( 1.0 + sin( nproj * truth * rad ) )
        do n=1,lngrd
          fm2(n)=(1.+sin( nproj * flat(n) ))/a2  ! 1/(map factor)
          rgln = flon(n) - dlamda0
          rgln = nproj * rgln
          rcsln(n) = nproj * cos( rgln ) 
          rsnln(n) = nproj * sin( rgln )
        enddo
      else if( nproj.eq.0 ) then
        a2 = cos( truth * rad )
        do n=1,lngrd
          fm2(n)=cos( flat(n) )/a2               ! 1/(map factor)
          rcsln(n) =  0.0
          rsnln(n) = -1.0 
        enddo
      endif
!
      call todxdy(fm2,fm2x,fm2y,igrd1,jgrd1,rdelx,rdely)
!
      do ij=1,lngrd
      fm2(ij) = 1.0/ ( fm2(ij) * fm2(ij) )
      fm2x(ij) = -2.0 * fm2(ij)**1.5 * fm2x(ij)
      fm2y(ij) = -2.0 * fm2(ij)**1.5 * fm2y(ij)
      enddo
      call maxmin(fm2,lngrd,1,1,1,'fm2')
      call maxmin(fm2x,lngrd,1,1,1,'fm2x')
      call maxmin(fm2y,lngrd,1,1,1,'fm2y')
!    
!....................................................
!     iunit=10
! do sig file or not
! ....................................................
        label(1)='GFS GRIB'
        write(6,*) label(1)
        do k=1,levs
          slin(k)=plevels(k)/plevels(1)
        enddo
        siin(1)=1.0
        siin(levs+1)=0.0
        do k=2,levs
          siin(k)=(slin(k-1)+slin(k))/2.0
        enddo
! use siin3d & slin3d to transfer si & sl
        do i=1,lngrd
              siin3d(i,:)=siin
              slin3d(i,:)=slin
        enddo
!
! interpolate real terrain for int_rchgr modify boundary mountain.
!      allocate (  cgg(1440,721) )
!       iaaa=26
!       rewind(iaaa)
!       read(iaaa) cg !!!!! link aaa to unit number of terr file(grib2).
!        read(iaaa) cgg !!!!! link aaa to unit number of terr file(grib2), real8
!        ix1=int((60.-0.)/0.25+1.)
!        ix2=int((180.-0.)/0.25+1.-1.)
!        iy1=int(-20.-(-90.)/0.25+1.)
!        iy2=int(70.-(-90.)/0.25+1.)
!       ii=0
!        do j = iy1,iy2
!          do i = ix1,ix2 
!            ii=ii+1
!            cg(ii)=cgg(i,j)
!          enddo
!        enddo
!      deallocate (  cgg )
!yj        cg=cfld(:,1,8)
        if (abs(rproj).eq.1) then
!       call int_cg2rg(cg,terr,1)
        call int_cg2rg(terr_gfs,terr,1)
        else
!       call int_cgtorg(cg,terr,1)
        call int_cgtorg(terr_gfs,terr,1)
        endif
!
!yj        write(6,210) (idate(i),i=1,4),fhour,                             &
        write(6,210) idate,fhour,                             &
     &            (siin(k),k=1,levs+1),(slin(k),k=1,levs)
210     format(5x,' input date and fcst hour ',i10,f7.1/(2x,g13.6))
!yj210     format(5x,' input date and fcst hour ',4i5,f7.1/(2x,g13.6))
! gz
!        cg=ge_gfs
!yj        cg=cfld(:,1,5)
        if (abs(rproj).eq.1) then
        call int_cg2rg(ge_gfs,gz,1)
        else
        call int_cgtorg(ge_gfs,gz,1)
        endif
! q
        cg=log(plevels(1)/10.0)
        if (abs(rproj).eq.1) then
        call int_cg2rg(cg,q,1)
        else
        call int_cgtorg(cg,q,1)
        endif
! te 
        do k=1,levs
!          cg=te_gfs(:,k)
!yj          cg=cfld(:,k,1)
          if (abs(rproj).eq.1) then
          call int_cg2rg(te_gfs(:,k),te(1,k),1)
          else
          call int_cgtorg(te_gfs(:,k),te(1,k),1)
          endif
        enddo
        do k=1,levs
!          cg=u_gfs(:,k)
!yj          cg=cfld(:,k,6)
          if (abs(rproj).eq.1) then
          call int_cg2rg(u_gfs(:,k),uu(1,k),1)
          else
          call int_cgtorg(u_gfs(:,k),uu(1,k),1)
          endif
!          cg=v_gfs(:,k)
!yj          cg=cfld(:,k,7)
          if (abs(rproj).eq.1) then
          call int_cg2rg(v_gfs(:,k),vv(1,k),1)
          else
          call int_cgtorg(v_gfs(:,k),vv(1,k),1)
          endif
        enddo
! rq
        do k=1,levs
!          cg=sp_gfs(:,k)
!yj            cg=cfld(:,k,2)
          if (abs(rproj).eq.1) then
          call int_cg2rg(sp_gfs(:,k),rq(1,k,1),1)
          else
          call int_cgtorg(sp_gfs(:,k),rq(1,k,1),1)
          endif
        enddo
! convert temperature to virtual temperature
!
        do k=1,levs
        do i=1,lngrd
           te(i,k)=(1.+0.61*rq(i,k,1))*te(i,k)
        enddo
        enddo
!ccj
       if(ntoz.eq.1) then
         do ntt = 1,ntoz
           ktrace = 1+ntt
           do k=1,levs
!o3             
!             cg=o3_gfs(:,k)
!yj             cg=cfld(:,k,4)
             if (abs(rproj).eq.1) then
             call int_cg2rg(o3_gfs(:,k),rq(1,k,ktrace),1)
             else
             call int_cgtorg(o3_gfs(:,k),rq(1,k,ktrace),1)
             endif
           enddo
         enddo
       endif
!
       if(ncldb.ge.1) then
         do nc = 1,ncldb
           kcloud = 1+ntoz+nc
           do k=1,levs
!cloud water             
!yj             cg=cfld(:,k,3)
             if (nc.eq.1) then
               cg(:)=clw_gfs(:,k)
             elseif (nc.eq.2) then
               cg(:)=rain_gfs(:,k)
             elseif (nc.eq.3) then
               cg(:)=qice_gfs(:,k)
             elseif (nc.eq.4) then
               cg(:)=snow_gfs(:,k)
             elseif (nc.eq.5) then
               cg(:)=grpl_gfs(:,k)
             endif
             if (abs(rproj).eq.1) then
             call int_cg2rg(cg,rq(1,k,kcloud),1)
             else
             call int_cgtorg(cg,rq(1,k,kcloud),1)
             endif
           enddo
         enddo
       endif
!
       if(ncldb.eq.0.and.ncld.ge.1) then
         do kc = 1,ncld
           kcloud = 1+ntoz+kc
           do k=1,levs
           do i=1,clngrd
             rq(i,k,kcloud)=0.
           enddo
           enddo
         enddo
       endif
!cc
!yj!#ifdef NONHYD
!yj!! pp
!yj!        do k=1,levs
!yj!         do i=1,lngrd
!yj!          pp(i,k)=q(i)+log(slin(k))
!yj!         enddo
!yj!        enddo
!yj!! tt
!yj!        do k=1,levs
!yj!         do i=1,lngrd
!yj!          tt(i,k)=te(i,k)
!yj!         enddo
!yj!        enddo
!yj!! ww
!yj!        do k=1,levs+1
!yj!         do i=1,lngrd
!yj!          ww(i,k)=0.0
!yj!         enddo
!yj!        enddo
!yj!#endif
!.................................................................
! define cgrid project transformation for u and v
      corient=(clon2-clon1)/2.0
      ctruth=(clat2-clat1)/2.0
      truth=ctruth
      nproj = cproj
      orient = corient
      dlamda0 = orient * rad + hfpi
        do n=1,lngrd
          if( nproj.eq.1 .or. nproj.eq.-1 ) then
            a2 = ( 1.0 + sin( nproj * truth * rad ) )
            cm2(n)=(1.+sin( nproj * flat(n) ))/a2  ! 1/(map factor)
            call inp_ll2xy(flat(n),flon(n),x,y)
            call inp_xy2ll(x,y,clat,clon)
            rgln = clon - dlamda0
            rgln = nproj * rgln
            csll = nproj * cos( rgln )
            snll = nproj * sin( rgln )
          else if( nproj.eq.0 ) then
            a2 = cos( truth * rad )
            cm2(n)=cos( flat(n) )/a2               ! 1/(map factor)
            csll =  0.0
            snll = -1.0 
          else if (nproj.eq.4) then
            cm2(n)=1
            csll = 0.0
            snll = -1.0
          endif
          csmm = rcsln(n)
          snmm = rsnln(n)
          cmorm = sqrt( cm2(n) / fm2(n) )
          rcsln(n) = ( csll * csmm + snll * snmm ) * cmorm
          rsnln(n) = ( csll * snmm - snll * csmm ) * cmorm
        enddo

!
! rotate u v to new projection
        do k=1,levs
          do n=1,lngrd
            hold = uu(n,k)
            uu(n,k) = hold * rcsln(n) - vv(n,k) * rsnln(n)
            vv(n,k) = hold * rsnln(n) + vv(n,k) * rcsln(n)
          enddo
        enddo
!
        call maxmin(gz,lngrd,1,1,1,'gz full')
        call maxmin( q,lngrd,1,1,1,' q full')
        call maxmin(te,lngrd,levs,1,levs,'te full')
        call maxmin(rq(1,1,1),lngrd,levs,1,levs,'rq full')
        if(ntoz.eq.1)                                                    &         
     &   call maxmin(rq(1,1,2),lngrd,levs,1,levs,'ozone full')
        if(ncld.eq.1) then
         call maxmin(rq(1,1,3),lngrd,levs,1,levs,'cloud water full')
        else
         call maxmin(rq(1,1,3),lngrd,levs,1,levs,'cloud water full')
         call maxmin(rq(1,1,4),lngrd,levs,1,levs,'rain full')
         call maxmin(rq(1,1,5),lngrd,levs,1,levs,'cloud ice full')
         call maxmin(rq(1,1,6),lngrd,levs,1,levs,'snow full')
         call maxmin(rq(1,1,7),lngrd,levs,1,levs,'graupel full')
        endif
        call maxmin(uu,lngrd,levs,1,levs,'uu full')
        call maxmin(vv,lngrd,levs,1,levs,'vv full')
!yj!#ifdef NONHYD
!yj!        call maxmin(pp,lngrd,levs,1,levs,'pp full')
!yj!        call maxmin(tt,lngrd,levs,1,levs,'tt full')
!yj!        call maxmin(ww,lngrd,levs+1,1,levs+1,'ww full')
!yj!#endif
        call maxmin(fm2,lngrd,1,1,1,'fm2 full')
        call maxmin(fm2x,lngrd,1,1,1,'fm2x full')
        call maxmin(fm2y,lngrd,1,1,1,'fm2y full')
        call maxmin(flat,lngrd,1,1,1,'flat full')
        call maxmin(flon,lngrd,1,1,1,'flon full')
!
7777    print *,' === write regional data at ',fhour
!
! ...............................................
!
        print *,' end of sig file interpolation. '
!
! ==============================================================
!.............................................................
! do sfc 
!!!!!gfs!!!!!      
!yj      if (cwbgfs.eq.1) then
! yj add binary output from gfs
!yj!!!!!!!!!!!!!!!--MPMD-- surface field input!!!!!!!!!!!!!
!yj        close(npgb)
!yj        print*,'replace by cwbgfs --- start reading rsm_sfc.out'
!yj         open(npgb,file=fncwbsfc,status='old', &
!yj             form='unformatted',iostat=ios2)
!yj         if( ios2.ne.0) then
!yj           print *,' error open rsm_sfc.out'
!yj           stop
!yj         endif 
!yj         do n=1,nsfcfld
!yj             if (n.eq.nsfcfld) then
!yj               print*,'n=nsfcfld=',n
!yj               read(npgb,iostat=eof) csfcfld(:,n)
!yj               if (eof.ne.0) print*,'error reading file'
!yj               print*,'read the end of the sfc file,eof=',eof
!yj               print*,'n=8?',n
!yj!yj               print*,'replace with simk_grib2.bin from grib2 file'
!yj!yj               open(66,file='simk_grib2.bin',form='unformatted')
!yj!yj               read(66) csfcfld(:,n)
!yj!yj               close(66)
!yj               print*,'checking (n)=(8)?',n
!yj             elseif (n.eq.7) then
!yj               read(npgb,iostat=eof) csfcfld(:,n)
!yj               print*,'n=7?',n
!yj               print*,'replace land-sea mask from grib2 file'
!yj               open(78,file='slmk_grib2.bin',form='unformatted')
!yj               read(78) csfcfld(:,n)
!yj               close(78)
!yj!yj             elseif(n.eq.4) then
!yj!yj               print*,'n=4?',n
!yj!yj               print*,'replace sheleg_grib2.bin from grib2 file'
!yj!yj               open(68,file='sheleg_grib2.bin',form='unformatted')
!yj!yj               read(68) csfcfld(:,n)
!yj!yj               close(68)
!yj             else  
!yj               read(npgb,iostat=eof) csfcfld(:,n)
!yj               if (eof.ne.0) print*,'error reading file'
!yj             endif  
!yj         enddo  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!the order of var in csfcfld is :
!1.tg 2.0-0.1m soil mois 3.0.4-1m soil mois 4.snr 5.0-0.1m soil temp 
!6.0.4-1m soil temp 7.land-sea mask 8.cice 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!yj         open(99,file='date.input.log',status='old',iostat=ios)
!yj         if(ios.ne.0) then
!yj           print*,' error open date.input.log'
!yj           stop
!yj         endif
!yj         read(99,*) idate(1)
!yj         read(99,*) idate(2)
!yj         read(99,*) idate(3)
!yj         read(99,*) idate(4)
!yj         read(99,*) fhour
!yj         close(99)
!yj      endif  !!!end if (cwbgfs.eq.1)
!yj        close(npgb)
!!!!!gfs!!!!!      
!

!      cg=snr_gfs
!yj      cg=csfcfld(:,4)            ! sheleg = snr
      call int_cg2rg(snr_gfs,sheleg,1)
!      cg=cice_gfs
!yj      cg=csfcfld(:,8)           ! ice mask = cice
! get sea and ice index
      call int_cg2rg(cice_gfs,simk,1)
      do 20 i=1,lngrd
        if( simk(i) .gt. 0.4 ) then
!yj          print*,'yj check simk'
          simk(i) = 2.0
        else
          simk(i) = 0.0
        endif
 20   continue
! get sea and land index
!      allocate (  cgg(1440,721) )
!      ibbb=27
!      rewind(ibbb)
!      read(ibbb) slg !!! link bbb to unit number of slmsk file
!      read(ibbb) cgg !!! link bbb to unit number of slmsk file
!        ix1=int((60.-0.)/0.25+1.)
!        ix2=int((180.-0.)/0.25+1.-1.)
!        iy1=int(-20.-(-90.)/0.25+1.)
!        iy2=int(70.-(-90.)/0.25+1.)
!       ii=0
!        do j = iy1,iy2
!          do i = ix1,ix2 
!            ii=ii+1
!            slg(ii)=cgg(i,j)
!          enddo
!        enddo
!      deallocate (  cgg )
!yj      cg=csfcfld(:,7)           ! land mask
!      call int_cg2rg(slg,slmk,1)
!      do 40 i=1,lngrd
!        if( slmk(i) .gt. 0.4 ) then
!          slmk(i) = 1.0
!        else
!          slmk(i) = 0.0
!        endif
! 40   continue
!     do 45 i=1,lngrd
!        slmsk(i) = max( simk(i), slmk(i) )
!  45 continue
!
! get regional land-sea mask
        read(nrslm) rgs  
        call sgl2ful(rgs,slmsk,lngrd)
        rewind nrslm
! before add sea ice info
! read climatology data
      call newidate(idate,0.0,idatenew)
!yj      call newidate(idate,fhour,idatenew)
      rewind(5)
      read(5,namclim)
      if (iclim.eq.1) then
        dflat=flat*180./pi
        dflon=flon*180./pi
        call rclima(idatenew(4),idatenew(2),idatenew(3),idatenew(1),      &
     &             fhour,                                                 &
     &             lngrd,lsoil, slmsk,dflat,dflon,                        &
     &             fntg3c,fnzorc,fnalbc,fnvegc,fnvetc,fnsotc,fnmskh,      &
     &             fnvmnc,fnvmxc,fnslpc,fnabsc,                           & 
     &             tg3clm,zorclm,albclm,vegclm,vetclm,sotclm,alfclm,      &
     &             vmnclm,vmxclm,slpclm,absclm)
      else
        write(6,*) 'Error: iclim should be 1 if pgb input'
        call abort
      endif
!
      do 45 i=1,lngrd
! add sea ice info
        slmsk(i) = max(slmsk(i),simk(i))
 45   continue
!
! update climate var sea ice info
      if (iclim.eq.1) then
        do i=1,lngrd
          if (slmsk(i).eq.2.0) then
#ifdef doNOAHLSM
          !for NOAHLSM
            vetclm(i)=15.0
            sotclm(i)=16.0
            zorclm(i)=1.0
#else
          !for OSULSM
            vetclm(i)=13.0
            sotclm(i)=9.0
            zorclm(i)=1.0
#endif
          endif
        enddo
      endif

!
!yj      slg=csfcfld(:,7)
!yj2019
       slg=slmk_gfs
!
!------ read from gsm, interpolate to rsm  -----
! label --- 1
! hour idate --- 2
        fhours=fhour
! tsea 3
!        cg=tg_gfs
!yj        cg=csfcfld(:,1)
!       call int_cg2rg(cg,tsea,1)
        call expnd(tg_gfs,slg,0,cgsea,iwork,imax,jmax)
        call expnd(tg_gfs,slg,1,cgland,iwork,imax,jmax)
        call int_cg2rg(cgsea,tsea,1)
        call int_cg2rg(cgland,rr,1)
        do i=1,lngrd
          if (slmsk(i).eq.1.0) tsea(i)=rr(i)
        enddo
        do 50 i=1,lngrd
          if( slmsk(i) .eq. 2.0e0 ) then
            tsea(i)= min( tsea(i), 271.20e0 )
          else if( slmsk(i).eq.0.0e0 ) then
            tsea(i)= max( tsea(i), 271.21e0 )
          else if( slmsk(i).eq.1.0e0 .and.                               &   
     &         sheleg(i).gt.0.0e0 ) then
            tsea(i)= min( tsea(i), 273.16e0 )
          endif
 50     continue
        call maxmin(  tsea,lngrd,1,1,1,'tsea')
!        call ful2sgl( tsea,rgs,lngrd)
! soil moisture ----- 4
!        cgl(:,1)=smc_gfs(:,1)       
!        cgl(:,2)=smc_gfs(:,2)
!yj        cgl(:,1)=csfcfld(:,2)
!yj        cgl(:,2)=csfcfld(:,3)
!       call int_cg2rg(cgl,smc,lsoil)
        do k=1,lsoil
#ifdef doNOAHLSM
          !for NOAHLSM, lsoil=4, k=1~4
          call expnd(smc_gfs(1,k),slg,0,cgsea,iwork,imax,jmax)
          call expnd(smc_gfs(1,k),slg,1,cgland,iwork,imax,jmax)
#else
          !for OSULSM, lsoil=2, k=1,3
          kk=(k-1)*2+1
          call expnd(smc_gfs(1,kk),slg,0,cgsea,iwork,imax,jmax)
          call expnd(smc_gfs(1,kk),slg,1,cgland,iwork,imax,jmax)
#endif
          call int_cg2rg(cgsea,smc(1,k),1)
          call int_cg2rg(cgland,rr,1)
          do i=1,lngrd
            if (slmsk(i).eq.1.0) then
               smc(i,k)=rr(i)
            else
               smc(i,k)=0.0
            endif
          enddo
        enddo
        call maxmin(   smc,lngrd,lsoil,1,lsoil,'smc')
!        call ful2sgl(  smc,rgls,lngrd*lsoil)
! snow depth ----- 5
!        cg=snr_gfs      
!yj        cg=csfcfld(:,4)
!       call int_cg2rg(cg,sheleg,1)
        call expnd(snr_gfs,slg,0,cgsea,iwork,imax,jmax)
        call expnd(snr_gfs,slg,1,cgland,iwork,imax,jmax)
        call int_cg2rg(cgsea,sheleg,1)
        call int_cg2rg(cgland,rr,1)
        do i=1,lngrd
          if (slmsk(i).eq.1.0) sheleg(i)=rr(i)
        enddo
        call maxmin(  sheleg,lngrd,1,1,1,'sheleg')
!        call ful2sgl( sheleg,rgs,lngrd)
! stc ------------ 6
!        cgl(:,1)=stc_gfs(:,1)       
!        cgl(:,2)=stc_gfs(:,2)       
!yj        cgl(:,1)=csfcfld(:,5)
!yj        cgl(:,2)=csfcfld(:,6)
!       call int_cg2rg(cgl,stc,lsoil)
        do k=1,lsoil
#ifdef doNOAHLSM
          !for NOAHLSM, lsoil=4, k=1~4
          call expnd(stc_gfs(1,k),slg,0,cgsea,iwork,imax,jmax)
          call expnd(stc_gfs(1,k),slg,1,cgland,iwork,imax,jmax)
#else
          !for OSULSM, lsoil=2, k=1,3
          kk=(k-1)*2+1
          call expnd(stc_gfs(1,kk),slg,0,cgsea,iwork,imax,jmax)
          call expnd(stc_gfs(1,kk),slg,1,cgland,iwork,imax,jmax)
#endif
          call int_cg2rg(cgsea,stc(1,k),1)
          call int_cg2rg(cgland,rr,1)
          do i=1,lngrd
            if (slmsk(i).eq.1.0) then
               stc(i,k)=rr(i)
            else
               stc(i,k)=0.0
            endif
          enddo
        enddo
        call maxmin(   stc,lngrd,lsoil,1,lsoil,'stc')
!        call ful2sgl(  stc,rgls,lngrd*lsoil)
! tg3 ------------ 7
!        if (iclim.eq.1) then
!          call ful2sgl( tg3clm,rgs,lngrd)
!        endif
! roughness ----- 8
!        if (iclim.eq.1) then
!          call ful2sgl( zorclm,rgs,lngrd)
!        endif
! cloud amount -- 9
! cloud base --- 10
! cloud top ---- 11
! consistency check among cv cvb and cvt
        do i=1,lngrd
            cv (i)=0.0
            cvt(i)=0.0
            cvb(i)=100.0
        enddo
        call maxmin(  cv, lngrd ,1,1,1,'cv')
!        call ful2sgl( cv,rgs, lngrd )
        call maxmin(  cvb, lngrd ,1,1,1,'cvb')
!        call ful2sgl( cvb,rgs, lngrd )
        call maxmin(  cvt, lngrd ,1,1,1,'cvt')
!        call ful2sgl( cvt,rgs, lngrd )
! four albedoes ------- 12
!        if (iclim.eq.1) then
!          call ful2sgl( albclm(1,1),rg4s(1,1),lngrd)
!          call ful2sgl( albclm(1,2),rg4s(1,2),lngrd)
!          call ful2sgl( albclm(1,3),rg4s(1,3),lngrd)
!          call ful2sgl( albclm(1,4),rg4s(1,4),lngrd)
!        endif
! sea land mask -- 13
        call maxmin(  slmsk,lngrd,1,1,1,'slmsk')
!        call ful2sgl( slmsk,rgs,lngrd)
! vegitation fraction --14
!        if (iclim.eq.1) then
!          call ful2sgl( vegclm,rgs,lngrd)
!        endif
! prepare 15 to end for end of file condition
        do i=1,lngrd
#ifdef doNOAHLSM
          canopy(i)=60.
#else
          canopy(i)=0.6
#endif
          f10m(i)=0.8
!          vtype(i)=1.0
!          stype(i)=1.0
!          facsf(i)=1.0
!          facwf(i)=1.0
          uustar(i) = 1.
          ffmm(i) = log(30.)
          ffhh(i) = log(30.)
        enddo
! vegitation canopy ------ 15
!        call ful2sgl( canopy,rgs,lngrd)
!
! f10m -------------- 16
!        call ful2sgl( f10m,rgs,lngrd)
!
! vegitation type ---- 17
!        if (iclim.eq.1) then
!          call ful2sgl( vetclm,rgs,lngrd)
!        endif
!
! soil type ---------- 18
!        if (iclim.eq.1) then
!          call ful2sgl( sotclm,rgs,lngrd)
!        endif
!
! zenith angle dependent factor ---- 19
!        if (iclim.eq.1) then
!          call ful2sgl( alfclm(1,1),rg2s(1,1),lngrd)
!          call ful2sgl( alfclm(1,2),rg2s(1,2),lngrd)
!        endif
!
! ustar ------------ 20
!        call ful2sgl( uustar,rgs,lngrd)
!
! momentum exchange similarity function -- 21
!        call ful2sgl( ffmm,rgs,lngrd)
!
! heat exchange similarity function -- 22
!        call ful2sgl( ffhh,rgs,lngrd)
!
        print *,' ****** nornal end of int_l2rinp ******* '
!
      deallocate ( rgs )
      deallocate ( slg )
      deallocate (  cg  )
!      deallocate (  cg  , cgl )
!      deallocate ( cg2, cg4 )
      deallocate ( blat ,blon )
!
      deallocate ( cgsea,cgland,rr )
      deallocate ( iwork )
      deallocate ( cm2 )
      deallocate ( dflat,dflon )
!
!      deallocate ( rain,cvrn )
      deallocate ( rcsln,rsnln )
!....
!      deallocate ( tg3   , vfrac  )
!      deallocate ( zorl )
!      deallocate ( vtype  )
!      deallocate ( stype )
!      deallocate ( alvsf ,  alvwf )
!      deallocate ( alnsf ,  alnwf )
!      deallocate ( facsf ,  facwf )
      deallocate ( slmk,simk )
!      deallocate ( simk )
!yj      deallocate ( cfld )
!yj      deallocate ( csfcfld )
!yj      deallocate ( spfh )
    
!
      return
      end
