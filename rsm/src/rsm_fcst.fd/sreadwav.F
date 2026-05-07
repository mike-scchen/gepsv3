#include <define.h>
      subroutine sreadwav(n,fhour,idate,                                 &
     &                 flat,flon,fm2,fm2x,fm2y,sl,si,                    &
     &                 fgz,fq,fte,fuu,fvv                                &
#ifdef NDSL
     &                 ,frqg                                             &
#else
     &                 ,frq                                              &
#endif
#ifdef NONHYD
     &                     ,fpp,ftt,fww,tov                              &
#endif
     &                 )
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  sreadwav
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  read regional sigma surface file and base field to obtain
!            regional perturbation in wave coefficient.
!
! program history log:
!
! usage:    call sreadwav(n,fhour,idate,flat,flon,fm2,fm2x,fm2y,sl,si,
!          1           fgz,fq,fte,fuu,fvv,frq)
!   input argument list:
!     n         - file input unit number
!
!   output argument list:
!     fhour     - forecast hour
!     idate     - integrer date of hour month day year.
!     flat      - latitudes
!     flon      - longitudes
!     fm2       - square of map factor
!     fm2x      - d (fm2) / dx
!     fm2y      - d (fm2) / dy
!     sl        - sigma layers
!     si        - sigma interfaces
!     fgz       - perturbation wave coefficient of terrain height
!     fq        - perturbation wave coefficient of log surface pressure (cb)
!     fte       - perturbation wave coefficient of virtual temperature
!     fuu       - perturbation wave coefficient of x-direction wind 
!     fvv       - perturbation wave coefficient of y-direction wind
!     frq       - perturbation wave coefficient of specific humidity
!
!   input files: 
!       n
!
!   output files: none
!
!   subprograms called:
!       maxmin  getptb  gdtocc gdtocs gdtosc 
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
!.....
!.....developed by hann-ming henry juang (nov 1990, version 0)......
!.....use u v as input, update u v tendency by div vor tendency.....
!.....
!ccj
#include <paramodel.h>
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
#endif
      real, allocatable :: tmp1(:),tmp2(:)
      real, allocatable :: tmpgrd(:)
#include <rscomltb.h>
#ifdef MP
#define LNGRDS lngrdp
#define LNWAVS lnwavp
#else
#define LNGRDS lngrd
#define LNWAVS lnwav
#endif
!cc
      parameter(levmax=100)
      real fhour
      dimension idate(4),si(levrp1),sl(levr),                            &
     & flat(LNGRDS), flon(LNGRDS),                                       &
     &  fm2(LNGRDS), fm2x(LNGRDS), fm2y(LNGRDS),                         &
     & fgz(LNWAVS), fq(LNWAVS),fte(LNWAVS,levr),                         &
     & fuu(LNWAVS,levr),fvv(LNWAVS,levr),                                &
#ifdef NDSL
     & frqg(LNGRDS,levr,ntotal)                                          &
#else
     & frq(LNWAVS,levr,ntotal)                                           &
#endif
#ifdef NONHYD
     &,fpp(LNWAVS,levr),ftt(LNWAVS,levr),fww(LNWAVS,levr+1),tov(levr)    &
#endif
     &,dummy(1)
!cjun
      dimension tmp(lngrd,levr),tmpk1(lngrd)
!
#ifdef MP
      if( mype.eq.master ) then
#endif
      print *,' ---- begin sreadwav -----'
#ifdef MP
      endif
#endif
!
#ifdef MP
      allocate (tmp1(lnwav))
      allocate (tmp2(lnwav))
      allocate (tmpgrd(lngrd))
#endif
!
#ifdef MP
      if( mype.eq.master ) then
#endif
      rewind n
      read(n) fhour,idate,(si(k),k=1,levrp1),(sl(k),k=1,levr)
      print *,' fhour idate si sl ',fhour,idate,si,sl
#ifdef MP
      endif
!yj      call rmpbcastc(label,32)
      call rmpbcastr(fhour,1)
      call rmpbcasti(idate,4)
      call rmpbcastr(si,levrp1)
      call rmpbcastr(sl,levr)
#endif
!-----------------------------------------gz
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
      call rmpsf2p(tmp1,lnwav,fgz,lnwavp,1)
#endif
!---------------------------------------- q
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
      call rmpsf2p(tmp1,lnwav,fq,lnwavp,1)
#endif
!--------------------------------------- te
      do 100 k=1,levr
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
      call rmpsf2p(tmp1,lnwav,fte(1,k),lnwavp,1)
#endif
100   continue
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef NONHYD
      read(n) (tov(k),k=1,levr)
#endif
#ifdef MP
      endif
#ifdef NONHYD
      call rmpbcastr(tov,levr)
#endif
#endif

!---------------------------------------- uu vv
      do 110 k=1,levr
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmp1(i),i=1,lnwav)
      read(n) (tmp2(i),i=1,lnwav)
#ifdef MP
      endif
      call rmpsf2p(tmp1,lnwav,fuu(1,k),lnwavp,1)
      call rmpsf2p(tmp2,lnwav,fvv(1,k),lnwavp,1)
#endif
110   continue
!-------------------------------------------- rq
#ifdef NDSL
      do 121 nc = 1,ntotal
      do 120 k=1,levr
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
      call rmpgf2p(tmpgrd,igrd1,jgrd1,frqg(1,k,nc),igrd1p,jgrd1p,1)
#endif 
120   continue
121   continue
!else not NDSL
#else
      do 121 nc = 1,ntotal
      do 120 k=1,levr
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
      call rmpsf2p(tmp1,lnwav,frq(1,k,nc),lnwavp,1)
#endif 
120   continue
121   continue
!endif NDSL
#endif
!
#ifdef NONHYD
!c--------------------------------------- pn
      do k=1,levr
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
      call rmpsf2p(tmp1,lnwav,fpp(1,k),lnwavp,1)
#endif
      enddo
!--------------------------------------- tn
      do k=1,levr
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmp1(i,k),i=1,lnwav)
#ifdef MP
      endif
      call rmpsf2p(tmp1,lnwav,ftt(1,k),lnwavp,1)
#endif
      enddo
!--------------------------------------- wn
      do k=1,levr+1
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
      call rmpsf2p(tmp1,lnwav,fww(1,k),lnwavp,1)
#endif
      enddo
!
#endif
!-------------------------------------------- fm2
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
      call rmpgf2p(tmpgrd,igrd1,jgrd1,fm2,igrd1p,jgrd1p,1)
#endif
!-------------------------------------------- fm2x
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
      call rmpgf2p(tmpgrd,igrd1,jgrd1,fm2x,igrd1p,jgrd1p,1)
#endif
!-------------------------------------------- fm2y
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
      call rmpgf2p(tmpgrd,igrd1,jgrd1,fm2y,igrd1p,jgrd1p,1)
#endif
!-------------------------------------------- flat
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
      call rmpgf2p(tmpgrd,igrd1,jgrd1,flat,igrd1p,jgrd1p,1)
#endif
!-------------------------------------------- flon
#ifdef MP
      if( mype.eq.master ) then
#endif
      read(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
      call rmpgf2p(tmpgrd,igrd1,jgrd1,flon,igrd1p,jgrd1p,1)
#endif
!
#ifdef MP
      if( mype.eq.master ) then
#endif
      print *,' ------ end sreadwav ----- '
#ifdef MP
      endif
#endif
!
#ifdef MP
      deallocate (tmpgrd)
      deallocate (tmp2)
      deallocate (tmp1)
#endif
!
      return
      end
