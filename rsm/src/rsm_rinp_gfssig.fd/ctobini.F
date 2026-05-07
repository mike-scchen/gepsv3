#include <define.h>
      subroutine ctobini(clat,clon,blat,blon                                 &
!yj      subroutine ctobini(clats,clons,blat,blon                                 &
!      subroutine ctobini(inic,blat,blon                                 &
#ifdef MP
     &           ,mype,master                                           &
#endif
     &            )
!
      use physcons, pi => con_pi
!yj add
!     use module_rinpbase, only : clatm => flat, clonm=>flon
!ccj
#include <paramodel.h>
#include <rsparltb.h>
!cc
!     dimension xlonin(cigrd1),xlatin(cjgrd1)
      dimension blat(lngrdb),blon(lngrdb)        &
     &         ,clat(cigrd1,cjgrd1)              &
     &         ,clon(cigrd1,cjgrd1)
!yj      real(kind=kind_io4) clats(cigrd1,cjgrd1),clons(cigrd1,cjgrd1)
!
      common/comctobi/      n11(lngrdb)
      common/comctob/                                                   &
     &          xc0(lngrdb),xc1(lngrdb),xc2(lngrdb),xc3(lngrdb)         &
     &         ,yc0(lngrdb),yc1(lngrdb),yc2(lngrdb),yc3(lngrdb)
!
#ifdef MP
      if (mype.eq.master) then
#endif
      print 1234
 1234 format(' = in routine ctobini for horizontal interpolation = ')
!
      hfpi = pi / 2.0
      twopi = pi * 2.0
!         open(99,file='rsm_xlon_202110010000.txt',status='unknown', &
!             form='formatted',iostat=ierr)
!         if (ierr/=0) stop "rsmout_sigp: open file xlon fail"
!         read(99,*) xlonin
!         close(99)
!         open(98,file='rsm_xlat_202110010000.txt',status='unknown', &
!             form='formatted',iostat=ierr)
!         if (ierr/=0) stop "rsmout_sigp: open file xlat fail"
!         read(98,*) xlatin
!         close(98)
!     do j=1,cjgrd1
!     do i=1,cigrd1
!       clat(i,j)=xlatin(j)/180.*pi
!       clon(i,j)=xlonin(i)/180.*pi
!     enddo
!     enddo
      call maxmin(clat,cigrd1*cjgrd1,1,1,1,'lat from c')
      call maxmin(clon,cigrd1*cjgrd1,1,1,1,'lon from c')
      call maxmin(blat,lngrdb,1,1,1,'lat from b')
      call maxmin(blon,lngrdb,1,1,1,'lon from b')
      x00=clon(1,1); y00=clat(1,1)
      x11=clon(cigrd1,cjgrd1); y11=clat(cigrd1,cjgrd1)
! estimated cdelx and cdely      
      cdelx=(x11-x00)/(cigrd1-1)
      cdely=(y11-y00)/(cjgrd1-1)
#ifdef MP
      endif
!     call rmpbcastr(x00,1)
!     call rmpbcastr(y00,1)
!     call rmpbcastr(x11,1)
!     call rmpbcastr(y11,1)
!     call rmpbcastr(cdelx,1)
!     call rmpbcastr(cdely,1)
!     call rmpbcastr(clat,cigrd1*cjgrd1)
!     call rmpbcastr(clon,cigrd1*cjgrd1)
#endif
!
      n = 0
      do 1000 i=1,lngrdb
        x=blon(i); y=blat(i)
        xlon=(x-x00)/cdelx+1
        xlat=(y-y00)/cdely+1
        lon=xlon
        lat=xlat
        if( lon.lt.1 .or. lon.gt.cigrd1 ) then
           print *,x,' out side of ',x00,' and ',x11
           stop
        else if( lat.lt.1 .or. lat.gt.cjgrd1 ) then
           print *,y,' out side of ',y00,' and ',y11
           stop
        endif
        lonp1 = lon + 1
        latp1 = lat + 1
        if (.not.(((x.ge.clon(lon,lat)).and.(x.le.clon(lonp1,latp1)))   &
     &     .and. ((y.ge.clat(lon,lat)).and.(y.le.clat(lonp1,latp1)))))  &
     &     then
           print*,'not good, x,y=',x,y,clon(lon,lat),clon(lonp1,latp1)  &
                                  ,clat(lon,lat),clat(lonp1,latp1)
           print*,'please find the correct index (lon,lat) so that',    &
     &            'x in between clon(lon,lat) and clon(lonp1,latp1)',   &
     &            'y in between clat(lon,lat) and clat(lonp1,latp1)'
        endif

        dclat = clat(lonp1,latp1)-clat(lon,lat)
        dclon = clon(lonp1,latp1)-clon(lon,lat)
        dblat1 = (y-clat(lon,lat))/dclat
        dblon1 = (x-clon(lon,lat))/dclon
        n = n + 1
!
        n11(n) = lon   + (lat   -1)*cigrd1
        yr = dblat1 - 0.5
        b1 = 0.5*(yr*yr-0.25)
        yc3(n) = b1*(0.5+yr/3.)
        yc2(n) = (0.5+yr)*(1.-b1)
        yc1(n) = (0.5-yr)*(1.-b1)
        yc0(n) = b1*(0.5-yr/3.)
        xr = dblon1 - 0.5
        b1 = 0.5*(xr*xr-0.25)
        xc3(n) = b1*(0.5+xr/3.)
        xc2(n) = (0.5+xr)*(1.-b1)
        xc1(n) = (0.5-xr)*(1.-b1)
        xc0(n) = b1*(0.5-xr/3.)
!
1000  continue
!
      return
      end
