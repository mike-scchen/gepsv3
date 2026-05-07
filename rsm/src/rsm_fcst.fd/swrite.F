#include <define.h>
!=============================================================
      subroutine swrite(n,label,fhour,idate,                             &  
     &                  flat,flon,fm2,fm2x,fm2y,sl,si,                   &
     &                  fgz,fq,fte,fuu,fvv                               &
#ifdef NDSL
     &                  ,frqg                                            &
#else
     &                  ,frq                                             &
#endif
#ifdef NONHYD
     &                            ,fpp,ftt,fww                           &
#endif
     &                                        )
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  swrite
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  write regional sigma surface file by adding base field to
!            regional perturbation from  wave gridficient.
!
! program history log:
!
! usage:    call swrite(n,fhour,idate,flat,flon,fm2,fm2x,fm2y,sl,si,
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
      use machine , only : kind_io4
      use sigmap_cons , only : aki, bki
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
      real, allocatable ::grid(:,:),grdv(:,:),gridk1(:)
#include <rscomltb.h>
#include <rscomloc.h>
!cc
      parameter(levmax=100)
      character*8 label(4)
      real(kind=kind_io4), allocatable :: tmps(:)
      real(kind=kind_io4) fhours,sis(levr+1),sls(levr)
      real(kind=kind_io4) dummy(2*levmax+1-levrp1-levr)                  &  
     &         ,ext(512-(6+2*levmax))
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
     &,fpp(LNWAVS,levr),ftt(LNWAVS,levr),fww(LNWAVS,levr+1)
      real, allocatable :: tbs(:,:)
      real*4 qs4(lngrd),sl4(levr),tempik                                 &
#endif
     &  ,swritetmp1(2)

      character cfhour*16,cform*40

!
      allocate (grid(LNGRDS,levr+1))
      allocate (grdv(LNGRDS,levr))
      allocate (gridk1(LNGRDS))
      allocate (tmps(lngrd))
      allocate (tmp1(lngrd))
      allocate (tmp2(lngrd))
#ifdef NONHYD
      allocate (tbs(lngrd,levr))
#endif
!c
#ifdef MP
      if( mype.eq.master ) then
#endif
      label(1)='ncep'
      label(2)='rsm'
      label(3)='mpi'
      label(4)='version'
!c
      write(n) label
      print 3000,label,n
3000  format(1h0,'swrite label  ',4a10,' n=',i3)
!
      do i=1,2*levmax+1-levrp1-levr
      dummy(i)=0.0
      enddo
      do i=1,512-(6+2*levmax)
      ext(i)=0.0
      enddo
!
      ext(1)=iwav1
      ext(2)=jwav1
      ext(3)=igrd1
      ext(4)=jgrd1
      ext(5)=levr
      ext(6)=2+levr*4+5
      ext(7)=rproj             !projection flag
      ext(8)=rtruth
      ext(9)=rorient
      ext(10)=rcenlat
      ext(11)=rcenlon
      ext(12)=rlftgrd
      ext(13)=rbtmgrd
      ext(14)=rdelx
      ext(15)=rdely
#ifdef NONHYD
      ext(16)=1.0
#endif
!
      fhours=fhour
      call ful2sgl(si,sis,levrp1)
      call ful2sgl(sl,sls,levr)
      write(n) fhours,idate,(sis(k),k=1,levrp1),(sls(k),k=1,levr)        &
     &         ,dummy,ext
      write(n) real(aki,4),real(bki,4)
!
      print *,' fhour idate si sl ',fhour,idate,si,sl
      print *,' aki ',aki
      print *,' bki ',bki
      print *,' other model experimental constants: '                    & 
     &       ,' wavei wavej gridi gridj xlayers totfld '                 &
     &         ,(ext(i),i=1,6)                                           &
     &       ,' proj truth orient cenlat cenlon '                        &
     &         ,(ext(i),i=7,11)                                          &
     &       ,' grdleft grdbotm delx dely '                              &
     &         ,(ext(i),i=12,15)
#ifdef MP
      endif
#endif
!------------------------------------------ gz
      call cctogd(fgz,gridk1,1)
      call getful(gridk1,gzb,wsltb,1)
#ifdef MP
      call rmpgp2f(gridk1,igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = gridk1(:)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      call ful2sgl(tmp1,tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!-------------------------------------------- q
      call cctogd(fq,gridk1,1)
      call getful(gridk1,qb,wsltb,1)
#ifdef MP
      call rmpgp2f(gridk1,igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = gridk1(:)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef NONHYD
      do i=1,lngrd
        qs4(i) = tmp1(i)
      enddo
#endif
      call ful2sgl(tmp1,tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!---------------------------------------------- te
      call cctogd(fte,grid,levr)
      call getful(grid,teb,wsltb,levr)
      do 300 k=1,levr
#ifdef MP
      call rmpgp2f(grid(1,k),igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = grid(:,k)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef NONHYD
      do i=1,lngrd
        tbs(i,k) = tmp1(i)
      enddo
#endif
      call ful2sgl(tmp1,tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
300   continue
!------------------------------------------------ uu vv

      call sctogd(fuu,grid,levr)
      call cstogd(fvv,grdv,levr)
      call getful(grid,uub,wsltb,levr)
      call getful(grdv,vvb,wsltb,levr)
      do 310 k=1,levr
#ifdef MP
      call rmpgp2f(grid(1,k),igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
      call rmpgp2f(grdv(1,k),igrd1p,jgrd1p,tmp2,igrd1,jgrd1,1)
#else
      tmp1(:) = grid(:,k)
      tmp2(:) = grdv(:,k)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      call ful2sgl(tmp1,tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
      call ful2sgl(tmp2,tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
310   continue
!------------------------------------------------- rq
!ccj
      do nc = 1,ntotal
#ifdef NDSL
      grid(:,1:levr)=frqg(:,:,nc)
#else
      call cctogd(frq(1,1,nc),grid,levr)
#endif
      call getful(grid,rqb(1,1,nc),wsltb,levr)
      do 320 k=1,levr
#ifdef MP
      call rmpgp2f(grid(1,k),igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = grid(:,k)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      call ful2sgl(tmp1,tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
320   continue
      enddo
!cc
#ifdef NONHYD
!c---------------------------------------------- pn
      call cctogd(fpp,grid,levr)
      do 330 k=1,levr
#ifdef MP
      call rmpgp2f(grid(1,k),igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)  
#else
      tmp1(:) = grid(:,k)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      sl4(k) = log( sl(k) )
      do i=1,lngrd
        tempik = qs4(i) + sl4(k)
        tmp1(i) = tmp1(i) + tempik
      enddo
      call ful2sgl(tmp1,tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
330   continue
!c---------------------------------------------- tn
      call cctogd(ftt,grid,levr)
      do 340 k=1,levr
#ifdef MP
      call rmpgp2f(grid(1,k),igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = grid(:,k)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      do i=1,lngrd
        tmp1(i) = tmp1(i) + tbs(i,k)
      enddo
      call ful2sgl(tmp1,tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
340   continue
!c---------------------------------------------- wn
      call cctogd(fww,grid,levr)
      call cctogd(fww(1,levr+1),grid(1,levr+1),1)
      do 350 k=1,levr+1
#ifdef MP
      call rmpgp2f(grid(1,k),igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = grid(:,k)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      call ful2sgl(tmp1,tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
350   continue
!c------------
#endif
!-------------------------------------------------- fm2
#ifdef MP
      call rmpgp2f(fm2,igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = fm2(:)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      call ful2sgl(tmp1, tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!-------------------------------------------------- fm2x
#ifdef MP
      call rmpgp2f(fm2x,igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = fm2x(:)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      call ful2sgl(tmp1, tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!-------------------------------------------------- fm2y
#ifdef MP
      call rmpgp2f(fm2y,igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = fm2y(:)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      call ful2sgl(tmp1, tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!-------------------------------------------------- flat
#ifdef MP
      call rmpgp2f(flat,igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = flat(:)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      call ful2sgl(tmp1, tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      endif
#endif
!-------------------------------------------------- flon
#ifdef MP
      call rmpgp2f(flon,igrd1p,jgrd1p,tmp1,igrd1,jgrd1,1)
#else
      tmp1(:) = flon(:)
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
      call ful2sgl(tmp1, tmps,lngrd)
      write(n) (tmps(i),i=1,lngrd)
#ifdef MP
      print *,' ------ end swrite ----- '
      endif
#endif
!
      deallocate (grid)
      deallocate (grdv)
      deallocate (gridk1)
      deallocate (tmps)
      deallocate (tmp1)
      deallocate (tmp2)
#ifdef NONHYD
      deallocate (tbs)
#endif
!
      return
      end
