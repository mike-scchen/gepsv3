#include <define.h>
      subroutine c2rini_mpi(clat,clon,rlat,rlon                                &
!yj      subroutine c2rini_mpi(clat,clon,rlat,rlon                                &
!yj      subroutine c2rini_mpi(inic,rlat,rlon                                &
#ifdef MP
     &           ,mype,master                                             &
#endif
     &           )
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  g2rini
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  prepare global grid to regional grid interpolation by storing
!            the interpolation coefficient and grid arrangement.
!
! program history log:
!
! usage:     g2rini(gclat,flat,flon)
!   input argument list:
!     gclat     - gaussian latitude of global in dimension of latg
!     flat      - latitude of regional grids in dimension of lngrd
!     flon      - longitude of regional grids in dimension of lngrd
!
!   output argument list:
!
!   common blocks:
!       comg2r
!
!   input files: none
!
!   output files: none
!
!   subprograms called: none
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
      use machine , only : kind_io4
!yj add
!     use module_rinpbase, only : clatm => flat, clonm=>flon
!ccj
#include <paramodel.h>
#ifdef MP
#define LNGRDS lngrdp
#else
#define LNGRDS lngrd
#endif
!cc
!
      dimension rlat(LNGRDS),rlon(LNGRDS)                                   &
     &         ,clat(cigrd1,cjgrd1)                                       &
     &         ,clon(cigrd1,cjgrd1)
!yj      real(kind=kind_io4) clats(cigrd1,cjgrd1),clons(cigrd1,cjgrd1)
!
      common/comc2r/                                                      &
     &          n00(LNGRDS),n10(LNGRDS)                                     &
     &         ,n11(LNGRDS),n01(LNGRDS)                                     &
     &         ,d00(LNGRDS),d10(LNGRDS)                                     &
     &         ,d11(LNGRDS),d01(LNGRDS)
!
      save icall
      data icall/0/
!
      if( icall .eq. 0 ) then
      icall = 1
#ifdef MP
      if (mype.eq.master) then
#endif
      print 1234
 1234 format(' = in routine c2rini for horizontal interpolation = ')
!yj      rewind inic
!yj      nskip=4+3*levr+(1+ntoz+ncldb)*levr+3
!yj#ifdef N2R
!yj      nskip=nskip+3*levr+1
!yj#endif
!yj      do nc=1,nskip
!yj        read(inic)
!yj      enddo
!yj      read(inic) clats
!yj      read(inic) clons
!yj      call sgl2ful(clats,clat,cigrd1*cjgrd1)
!yj      call sgl2ful(clons,clon,cigrd1*cjgrd1)
!     do j=1,cjgrd1
!     do i=1,cigrd1
!       ii=cigrd1*(j-1)+i
!       clat(i,j)=clatm(ii)
!       clon(i,j)=clonm(ii)
!     enddo
!     enddo
      call maxmin(clat,cigrd1*cjgrd1,1,1,1,'lat from c')
      call maxmin(clon,cigrd1*cjgrd1,1,1,1,'lon from c')
      call maxmin(rlat,LNGRDS,1,1,1,'lat from r')
      call maxmin(rlon,LNGRDS,1,1,1,'lon from r')
      x00=clon(1,1); y00=clat(1,1)
      x11=clon(cigrd1,cjgrd1); y11=clat(cigrd1,cjgrd1)
! estimated cdelx and cdely      
      cdelx=(x11-x00)/(cigrd1-1)
      cdely=(y11-y00)/(cjgrd1-1)
      print *,' c(1,1)=',x00,y00
      print *,' c(im,jm)=',x11,y11
      print *,' c delx dely =',cdelx,cdely
#ifdef MP
      endif
      call rmpbcastr(x00,1)
      call rmpbcastr(y00,1)
      call rmpbcastr(x11,1)
      call rmpbcastr(y11,1)
      call rmpbcastr(cdelx,1)
      call rmpbcastr(cdely,1)
!     call rmpbcastr(clat,cigrd1*cjgrd1)
!     call rmpbcastr(clon,cigrd1*cjgrd1)
#endif
!
      n = 0
      do 1000 i=1,LNGRDS
!       call ll2xy(rlat(i),rlon(i),x,y)
        x=rlon(i); y=rlat(i)
        xlon=(x-x00)/cdelx+1
        xlat=(y-y00)/cdely+1
        lon=xlon
        lat=xlat
        if( lon.lt.1 .or. lon.gt.cigrd1 ) then
           print *,'c2rini_mpi',x,' out side of ',x00,' and ',x11
           stop
        else if( lat.lt.1 .or. lat.gt.cjgrd1 ) then
           print *,'c2rini_mpi',y,' out side of ',y00,' and ',y11
           stop
        endif
!       call ll2xy(clat(lon  ,lat  ),clon(lon  ,lat  ),xc00,yc00)
             xc00=clon(lon,lat); yc00=clat(lon,lat)
!       call ll2xy(clat(lon+1,lat+1),clon(lon+1,lat+1),xc11,yc11)
             xc11=clon(lon+1,lat+1); yc11=clat(lon+1,lat+1)
!       call ll2xy(clat(lon+1,lat  ),clon(lon+1,lat  ),xc10,yc10)
             xc10=clon(lon+1,lat  ); yc10=clat(lon+1,lat  )
!       call ll2xy(clat(lon  ,lat+1),clon(lon  ,lat+1),xc01,yc01)
             xc01=clon(lon  ,lat+1); yc01=clat(lon  ,lat+1)
        r00=(x-xc00)*(x-xc00)+(y-yc00)*(y-yc00)
        r10=(x-xc10)*(x-xc10)+(y-yc10)*(y-yc10)
        r01=(x-xc01)*(x-xc01)+(y-yc01)*(y-yc01)
        r11=(x-xc11)*(x-xc11)+(y-yc11)*(y-yc11)
        if(r00.eq.0.0) then
          rsum=1.0
          r00=1.0
          r01=0.0
          r10=0.0
          r11=0.0
        elseif(r01.eq.0.0) then
          rsum=1.0
          r00=0.0
          r01=1.0
          r10=0.0
          r11=0.0
        elseif(r11.eq.0.0) then
          rsum=1.0
          r00=0.0
          r01=0.0
          r10=0.0
          r11=1.0
        elseif(r10.eq.0.0) then
          rsum=1.0
          r00=0.0
          r01=0.0
          r10=1.0
          r11=0.0
        else
          r00=1./r00
          r01=1./r01
          r10=1./r10
          r11=1./r11
          rsum=r00+r01+r10+r11
        endif
!
      n = n + 1
      d00(n) = r00 / rsum
      d10(n) = r10 / rsum
      d11(n) = r11 / rsum
      d01(n) = r01 / rsum
      n00(n) = lon   + (lat   -1)*cigrd1
      n10(n) = lon+1 + (lat   -1)*cigrd1
      n11(n) = lon+1 + (lat+1 -1)*cigrd1
      n01(n) = lon   + (lat+1 -1)*cigrd1
!
1000  continue
!
!     print *,' d00 to d11 should be all positive.'
!     call maxmin(d00,LNGRDS,1,1,1,'d00')
!     call maxmin(d10,LNGRDS,1,1,1,'d10')
!     call maxmin(d01,LNGRDS,1,1,1,'d01')
!     call maxmin(d11,LNGRDS,1,1,1,'d11')
!
      else
#ifdef MP
      if (mype.eq.master) then
#endif
      print *,' c2rini already called '
#ifdef MP
      endif
#endif
      endif
!
      return
      end
