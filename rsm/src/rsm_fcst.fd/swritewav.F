#include <define.h>
!=============================================================
      subroutine swritewav(n,fhour,idate,                                &  
     &                  flat,flon,fm2,fm2x,fm2y,sl,si,                   &
     &                  fgz,fq,fte,fuu,fvv                               &
#ifdef NDSL
     &                  ,frqg                                            &
#else
     &                  ,frq                                             &
#endif
#ifdef NONHYD
     &                            ,fpp,ftt,fww,tov                       &
#endif
     &                                        )
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  swritewav
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  write regional sigma surface file by adding base field to
!            regional perturbation from  wave gridficient.
!
! program history log:
!
! usage:    call swritewav(n,fhour,idate,flat,flon,fm2,fm2x,fm2y,sl,si,
!          1            fgz,fq,fte,fuu,fvv,frq)
!   input argument list:
!     n   - file input unit number
!     fhour - forecast hour
!     idate - integrer date of hour month day year.
!     flat  - latitudes
!     flon  - longitudes
!     fm2 - square of map factor
!     fm2x  - d (fm2) / dx
!     fm2y  - d (fm2) / dy
!     sl  - sigma layers
!     si  - sigma interfaces
!     fgz - perturbation wave gridficient of terrain height
!     fq  - perturbation wave gridficient of log surface pressure (cb)
!     fte - perturbation wave gridficient of virtual temperature
!     fuu - perturbation wave gridficient of x-direction wind
!     fvv - perturbation wave gridficient of y-direction wind
!     frq - perturbation wave gridficient of specific humidity
!
!   output argument list:
!
!   input files:
!   n
!
!   output files: none
!
!   subprograms called:
!   maxmin  getful  cctogd  cstogd sctogd
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
#ifdef MP
#define LNGRDS lngrdp
#define LNWAVS lnwavp
#else
#define LNGRDS lngrd
#define LNWAVS lnwav
#endif
      real, allocatable :: tmp1(:),tmp2(:)
      real, allocatable :: tmpgrd(:)
#include <rscomltb.h>
#include <rscomloc.h>
!cc
      parameter(levmax=100)
      real fhour
      dimension idate(4),si(levrp1),sl(levr),                            &
     & flat(LNGRDS), flon(LNGRDS),                                       &
     &   fm2(LNGRDS), fm2x(LNGRDS), fm2y(LNGRDS),                        &
     & fgz(LNWAVS), fq(LNWAVS),fte(LNWAVS,levr),                         &
     & fuu(LNWAVS,levr),fvv(LNWAVS,levr)                                 &
#ifdef NDSL
     & ,frqg(LNGRDS,levr,ntotal)                                         &
#else
     & ,frq(LNWAVS,levr,ntotal)                                          &
#endif
#ifdef NONHYD
     &,fpp(LNWAVS,levr),ftt(LNWAVS,levr),fww(LNWAVS,levr+1),tov(levr)    &
#endif
     &,dummy(1)

      character cfhour*16,cform*40

!
#ifdef MP
      allocate (tmp1(lnwav))
      allocate (tmp2(lnwav))
      allocate (tmpgrd(lngrd))
#endif
!c
#ifdef MP
      if( mype.eq.master ) then
#endif
!c
      write(n) fhour,idate,(si(k),k=1,levrp1),(sl(k),k=1,levr)
!
      print *,' fhour idate si sl ',fhour,idate,si,sl
#ifdef MP
      endif
#endif
#ifdef MP
#endif
!------------------------------------------ gz
#ifdef MP
      call rmpsp2f(fgz,lnwavp,tmp1,lnwav,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
#endif
!-------------------------------------------- q
#ifdef MP
      call rmpsp2f(fq,lnwavp,tmp1,lnwav,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
#endif
!---------------------------------------------- te
      do 300 k=1,levr
#ifdef MP
      call rmpsp2f(fte(1,k),lnwavp,tmp1,lnwav,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
#endif
300   continue
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef NONHYD
      write(n) (tov(k),k=1,levr)
#endif
#ifdef MP
      endif
#endif
!------------------------------------------------ uu vv

      do 310 k=1,levr
#ifdef MP
      call rmpsp2f(fuu(1,k),lnwavp,tmp1,lnwav,1)
      call rmpsp2f(fvv(1,k),lnwavp,tmp2,lnwav,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmp1(i),i=1,lnwav)
      write(n) (tmp2(i),i=1,lnwav)
#ifdef MP
      endif
#endif
310   continue
!------------------------------------------------- rq
!ccj
#ifdef NDSL
      do nc = 1,ntotal
      do 320 k=1,levr
#ifdef MP
         call rmpgp2f(frqg(1,k,nc),igrd1p,jgrd1p,tmpgrd,igrd1,jgrd1,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
#endif
320   continue
      enddo
!else not NDSL
#else
      do nc = 1,ntotal
      do 320 k=1,levr
#ifdef MP
         call rmpsp2f(frq(1,k,nc),lnwavp,tmp1,lnwav,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
#endif
320   continue
      enddo
!endif NDSL
#endif
!cc
#ifdef NONHYD
!c---------------------------------------------- pn
      do 330 k=1,levr
#ifdef MP
      call rmpsp2f(fpp(1,k),lnwavp,tmp1,lnwav,1)  
      if( mype.eq.master ) then
#endif
      write(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
#endif
330   continue
!c---------------------------------------------- tn
      do 340 k=1,levr
#ifdef MP
      call rmpsp2f(ftt(1,k),lnwavp,tmp1,lnwav,1)  
      if( mype.eq.master ) then
#endif
      write(n) (tmp1(i),i=1,lnwav)
#ifdef MP
      endif
#endif
340   continue
!c---------------------------------------------- wn
      do 350 k=1,levr+1
#ifdef MP
      call rmpsp2f(fww(1,k),lnwavp,tmp1,lnwav)  
      if( mype.eq.master ) then
#endif
      write(n) (tmp1(i,k),i=1,lnwav)
#ifdef MP
      endif
#endif
350   continue
!c------------
#endif
!-------------------------------------------------- fm2
#ifdef MP
      call rmpgp2f(fm2,igrd1p,jgrd1p,tmpgrd,igrd1,jgrd1,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!-------------------------------------------------- fm2x
#ifdef MP
      call rmpgp2f(fm2x,igrd1p,jgrd1p,tmpgrd,igrd1,jgrd1,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!-------------------------------------------------- fm2y
#ifdef MP
      call rmpgp2f(fm2y,igrd1p,jgrd1p,tmpgrd,igrd1,jgrd1,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!-------------------------------------------------- flat
#ifdef MP
      call rmpgp2f(flat,igrd1p,jgrd1p,tmpgrd,igrd1,jgrd1,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!-------------------------------------------------- flon
#ifdef MP
      call rmpgp2f(flon,igrd1p,jgrd1p,tmpgrd,igrd1,jgrd1,1)
      if( mype.eq.master ) then
#endif
      write(n) (tmpgrd(i),i=1,lngrd)
#ifdef MP
      print *,' ------ end swritewav ----- '
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
