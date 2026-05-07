      subroutine l2rinp
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
      use module_rinpbase, only :  &      
     &  siin,slin,delin,ciin,clin,rpiin &
     & ,siin3d,slin3d,terr &
!yj     & ,idate &
!comrinp2     
     & ,q,te,uu,vv,rq,gz &
#ifdef NONHYD
     &      ,pp,tt,ww &
#endif
     &      ,flat,flon,fm2,fm2x,fm2y
!yj add module   end -----------------------      
!yj add module_gfsdata start ---------------
      use module_gfsdata,only:             &
     & idate,fhour,                        &
     & te_gfs=>temp, sp_gfs=>spfh,         &
     & clw_gfs=>clwr, o3_gfs=>ozon,        &
     & ge_gfs=>geop, u_gfs=>u, v_gfs=>v
!yj add module_gfsdata   end ---------------

!      
!ccj
#include <inp_paramodel.h>  
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
      integer ios,ios2,cwbgfs,eof
      parameter (cwbgfs=1)
!
      real(kind=kind_io4) fhours,sis(levs+1),sls(levs)            & 
     &         ,cgs(clngrd), cgls(clngrd,lsoil)                         &
     &         ,cg2s(clngrd,2), cg4s(clngrd,4)                          &
     &         ,rgs(lngrd),  rgls(lngrd,lsoil)                          &
     &         ,rg2s(lngrd,2),  rg4s(lngrd,4)
!yj
!      integer :: idate(4)
!      real :: sl(levs),si(levs+1)
!      real :: flat(lngrd),flon(lngrd) 
!      real :: fm2(lngrd),fm2x(lngrd),fm2y(lngrd)
!      real :: gz(lngrd), q(lngrd) 
!      real :: te(lngrd,levs), rq(lngrd,levs,ntotal)
!      real :: uu(lngrd,levs), vv(lngrd,levs) 
!yj
!
      real, allocatable :: slg(:)
      real, allocatable ::  cg(:)  , cgl(:,:)
      real, allocatable :: cg2(:,:), cg4(:,:)
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
      real, allocatable :: rain(:),cvrn(:)
      real, allocatable :: rcsln(:),rsnln(:)
!....
!      real, allocatable :: tsea  (:), smc   (:,:)
!      real, allocatable :: sheleg(:), stc   (:,:)
      real, allocatable :: tg3   (:), vfrac (:)
      real, allocatable :: zorl  (:)
!      real, allocatable :: zorl  (:), slmsk (:)
!      real, allocatable :: cv    (:), cvb   (:)
!      real, allocatable :: cvt   (:)
!      real, allocatable :: f10m  (:), canopy(:)
      real, allocatable :: vtype (:)
      real, allocatable :: stype (:)
!      real, allocatable :: stype (:), uustar(:)
!      real, allocatable :: ffmm  (:),   ffhh(:)
      real, allocatable :: alvsf (:),  alvwf(:)
      real, allocatable :: alnsf (:),  alnwf(:)
      real, allocatable :: facsf (:),  facwf(:)
      real, allocatable :: slmk(:),simk(:)
!      real, allocatable :: tg3clm(:),zorclm(:)
!      real, allocatable :: albclm(:,:),vegclm(:)
!      real, allocatable :: vetclm(:),sotclm(:),alfclm(:,:)
!yj      real, allocatable :: cfld(:,:,:)
!yj      real, allocatable :: csfcfld(:,:)
!yj      real, allocatable :: spfh(:,:)
!
      character*500 fntg3c,fnzorc,fnalbc,fnvegc,fnvetc,fnsotc,fnmskh
      namelist/namclim/iclim,fntg3c,fnzorc,fnalbc,fnvegc,fnvetc,fnsotc, &
     &                 fnmskh
!yj      common /hydrid/ siin3d(lngrd,levs+1),slin3d(lngrd,levs)
!yj      common /llgz/ terr(lngrd)
      data fntg3c/'global_tg3clim.2.6x1.5.grb'/
      data fnmskh/'global_slmask.t126.grb'/
      data fnalbc/'global_albedo4.1x1.grb'/
      data fnsotc/'global_soiltype.1x1.grb'/
      data fnvegc/'global_vegfrac.1x1.grb'/
      data fnvetc/'global_vegtype.1x1.grb'/
      data fnzorc/'global_zorclim.1x1.grb'/
      character*500 fncwbpre, fncwbsfc
      data fncwbpre/'input.prel'/
      data fncwbsfc/'input.sfcl'/
!
      print *,' enter l2rinp '
!
!yj      call allocate_inpver_array
!yj      call allocate_hydrid_array
!yj      call allocate_llgz_array 
!yj      call allocate_comrinp1_array
!yj      call allocate_comrinp2_array

      allocate ( slg(clngrd) )
      allocate (  cg(clngrd)  , cgl(clngrd,lsoil) )
      allocate ( cg2(clngrd,2), cg4(clngrd,4) )
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
      allocate ( rain(lngrd),cvrn(lngrd) )
      allocate ( rcsln(lngrd),rsnln(lngrd) )
!....
      allocate ( rr    (lngrd) )
!      allocate ( tsea  (lngrd), smc   (lngrd,lsoil) )
!      allocate ( sheleg(lngrd), stc   (lngrd,lsoil) )
      allocate ( tg3   (lngrd), vfrac (lngrd) )
      allocate ( zorl  (lngrd) )
!      allocate ( zorl  (lngrd), slmsk (lngrd) )
!      allocate ( cv    (lngrd), cvb   (lngrd) )
!      allocate ( cvt   (lngrd) )
!      allocate ( f10m  (lngrd), canopy(lngrd) )
      allocate ( vtype (lngrd) )
      allocate ( stype (lngrd) )
!      allocate ( stype (lngrd), uustar(lngrd) )
!      allocate ( ffmm  (lngrd),   ffhh(lngrd) )
      allocate ( alvsf (lngrd),  alvwf(lngrd) )
      allocate ( alnsf (lngrd),  alnwf(lngrd) )
      allocate ( facsf (lngrd),  facwf(lngrd) )
      allocate ( slmk(lngrd),simk(lngrd) )
!      allocate ( tg3clm(lngrd),zorclm(lngrd) )
!      allocate ( albclm(lngrd,4),vegclm(lngrd) )
!      allocate ( vetclm(lngrd),sotclm(lngrd),alfclm(lngrd,2) )
!yj      allocate ( cfld(clngrd,levs,nfld) )
!yj      allocate ( csfcfld(clngrd,nsfcfld) )
!yj      allocate ( spfh(clngrd,levs) )
!
!.................................................................
! define regional grid
      call setlatlon(flat,flon,delx,dely,dlamda0)
      call maxmin(flat,lngrd,1,1,1,' flat ')
      call maxmin(flon,lngrd,1,1,1,' flon ')
!
!!!!!gfs!!!!!      
      if (cwbgfs.eq.1) then
! yj add binary output from gfs
        print*,'replace by cwbgfs --- start reading rsm_atm.out'
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
      endif  !!!end if (cwbgfs.eq.1)
!yj        close(npgb)
!!!!!gfs!!!!!      
      
      call setcgrid(blat,blon)

! .....................................................
! define coarse regional to regional grid interpolation
      call r2rini(blat,blon,flat,flon)
      call rtorini(blat,blon,flat,flon)
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
! do sig field 
! ....................................................
        label(1)='GFS GRIB'
        write(6,*) label
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
! interpolate real terrain for rchgr modify boundary mountain.
        iaaa=16  !!! fort.16 for rinpbase terr.file
        rewind(iaaa)
        read(iaaa) cg !!!!! link aaa to unit number of terr file(grib2).
!yj        cg=cfld(:,1,8)
        if (abs(rproj).eq.1) then
        call cg2rg(cg,terr,1)
        else
        call cgtorg(cg,terr,1)
        endif
!
!yj        write(6,210) (idate(i),i=1,4),fhour,                             &
        write(6,210) idate,fhour,                             &
     &            (siin(k),k=1,levs+1),(slin(k),k=1,levs)
210     format(5x,' input date and fcst hour ',i10,f7.1/(2x,g13.6))
!yj210     format(5x,' input date and fcst hour ',4i5,f7.1/(2x,g13.6))
! gz
        cg=ge_gfs
!yj        cg=cfld(:,1,5)
        if (abs(rproj).eq.1) then
        call cg2rg(cg,gz,1)
        else
        call cgtorg(cg,gz,1)
        endif
! q
        cg=log(plevels(1)/10.0)
        if (abs(rproj).eq.1) then
        call cg2rg(cg,q,1)
        else
        call cgtorg(cg,q,1)
        endif
! te 

        do k=1,levs
          cg=te_gfs(:,k)
!yj          cg=cfld(:,k,1)
          if (abs(rproj).eq.1) then
          call cg2rg(cg,te(1,k),1)
          else
          call cgtorg(cg,te(1,k),1)
          endif
        enddo
        do k=1,levs
          cg=u_gfs(:,k)
!yj          cg=cfld(:,k,6)
          if (abs(rproj).eq.1) then
          call cg2rg(cg,uu(1,k),1)
          else
          call cgtorg(cg,uu(1,k),1)
          endif
          cg=v_gfs(:,k)
!yj          cg=cfld(:,k,7)
          if (abs(rproj).eq.1) then
          call cg2rg(cg,vv(1,k),1)
          else
          call cgtorg(cg,vv(1,k),1)
          endif
        enddo
! rq
        do k=1,levs
          cg=sp_gfs(:,k)
!yj            cg=cfld(:,k,2)
          if (abs(rproj).eq.1) then
          call cg2rg(cg,rq(1,k,1),1)
          else
          call cgtorg(cg,rq(1,k,1),1)
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
             cg=o3_gfs(:,k)
!yj             cg=cfld(:,k,4)
             if (abs(rproj).eq.1) then
             call cg2rg(cg,rq(1,k,ktrace),1)
             else
             call cgtorg(cg,rq(1,k,ktrace),1)
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
             cg=clw_gfs(:,k)
!yj             cg=cfld(:,k,3)
             if (abs(rproj).eq.1) then
             call cg2rg(cg,rq(1,k,kcloud),1)
             else
             call cgtorg(cg,rq(1,k,kcloud),1)
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
#ifdef NONHYD
! pp
        do k=1,levs
         do i=1,lngrd
          pp(i,k)=q(i)+log(sl(k))
         enddo
        enddo
! tt
        do k=1,levs
         do i=1,lngrd
          tt(i,k)=te(i,k)
         enddo
        enddo
! ww
        do k=1,levs+1
         do i=1,lngrd
          ww(i,k)=0.0
         enddo
        enddo
#endif
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
        if(ncld.eq.1)                                                    &
     &   call maxmin(rq(1,1,3),lngrd,levs,1,levs,'cloud water full')
        call maxmin(uu,lngrd,levs,1,levs,'uu full')
        call maxmin(vv,lngrd,levs,1,levs,'vv full')
#ifdef NONHYD
        call maxmin(pp,lngrd,levs,1,levs,'pp full')
        call maxmin(tt,lngrd,levs,1,levs,'tt full')
        call maxmin(ww,lngrd,levs+1,1,levs+1,'ww full')
#endif

        call maxmin(fm2,lngrd,1,1,1,'fm2 full')
        call maxmin(fm2x,lngrd,1,1,1,'fm2x full')
        call maxmin(fm2y,lngrd,1,1,1,'fm2y full')
        call maxmin(flat,lngrd,1,1,1,'flat full')
        call maxmin(flon,lngrd,1,1,1,'flon full')
!
7777    print *,' === write regional data at ',fhour
!
!......................................................
        print *,' end of sig file interpolation. '
!
! ==============================================================
!.............................................................
        print *,' ****** nornal end of l2rinp ******* '
!
      deallocate ( slg )
      deallocate (  cg  , cgl )
      deallocate ( cg2, cg4 )
      deallocate ( blat ,blon )
!
      deallocate ( cgsea,cgland,rr )
      deallocate ( iwork )
      deallocate ( cm2 )
      deallocate ( dflat,dflon )
!
      deallocate ( rain,cvrn )
      deallocate ( rcsln,rsnln )
!....
      deallocate ( tg3   , vfrac  )
      deallocate ( zorl )
      deallocate ( vtype  )
      deallocate ( stype )
      deallocate ( alvsf ,  alvwf )
      deallocate ( alnsf ,  alnwf )
      deallocate ( facsf ,  facwf )
      deallocate ( slmk,simk )
!yj      deallocate ( cfld )
!yj      deallocate ( csfcfld )
!yj      deallocate ( spfh )
!
      return
      end
