      subroutine rupdate
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rupdate
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  update prognostic variables from forcing (or tendency).
!
! program history log:
!
! usage:    call rupdate
!   input argument list:
!
!   output argument list:
!
!   common block:
! rscomf
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
      use machine, only: kind_io4
!ccj
#include <paramodel.h>
#ifdef MP
#define LNWAVS lnwavp
#include <npesi.h>
#include <comrmpi.h>
#else
#define LNWAVS lnwav
#endif
#include <rscomspc.h>
#include <rscomio.h>
!
      dimension bc_mask(lnwav),bc_maskp(LNWAVS)
!      real,allocatable:: tmp(:),tmp1(:),tmpp(:),fullk1(:,:),gridk1(:)
!      real(kind=kind_io4),allocatable::fullk1sgl(:)
!      real :: mean
#ifdef KEout
! yj--add output the kinetic energy of each zonal wave number >
      character cfhour*3
      real,dimension(:,:),allocatable::  &
     &      utmp, vtmp, ke 
      real,dimension(:),allocatable::  &
     &      kezon, kemer
! yj--add output the kinetic energy of each zonal wave number <
#endif
      
!
#ifdef MP
      if( mype.eq.master ) then
#endif
      print *,'enter rupdate,fhour=',fhour
#ifdef MP
      endif
#endif
!cc
!
!------------------------------ ln (ps) -------------
      do 100 i=1,LNWAVS
      z(i) = qm(i) + dt2 * z(i)
100   continue
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef MBC
! mean bias correction
!      print *,'mype=',mype,'lnwavs=',LNWAVS,                       &
!     &    'jwav1=',jwav1
       z(1)=0.
!      remkh=mod(fhour,6.0)
!      if ( remkh .eq. 0) then
!        print *,'mype=',mype,' fhour=',fhour,'rem=',remkh
!        print *,'wave0, k=1,z(0)=',z(1)
!        print *,'wave0i,k=1,z(1)=',z(2)
!        print *,'wave1r,k=1,z(2)=',z(3)
!        print *,'wave1i,k=1,z(3)=',z(1+jwav1)
!        print *,'wave2r,k=1,z(4)=',z(2+jwav1)
!        print *,'wave2i,k=1,z(5)=',z(3+jwav1)
!        print *,'wave3r,k=1,z(6)=',z(1+jwav1*2)
!        print *,'wave3i,k=1,z(7)=',z(2+jwav1*2)
!        print *,'wave4r,k=1,z(8)=',z(3+jwav1*2)
!      endif
#endif
#ifdef MP
      endif
#endif
!
!-----------------------------------------------
#ifdef ORIGIN_THREAD
!$doacross share(tem,uum,vvm,dt2,y,x,w),local(i,k)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(tem,uum,vvm,dt2)
!mic$1 shared(  y,  x,  w)
!mic$1 private(i,k)
#endif
#ifdef OPENMP
!$omp parallel do private(i,k)
#endif
      do 200 k=1,levr
      do 201 i=1,LNWAVS
!-------------------------------------- temp ------
      y(i,k)=  tem(i,k) + dt2 * y(i,k)
!-------------------------------------- uu ---------
      x(i,k) = uum(i,k) + dt2 * x(i,k)
! ------------------------------------- vv ---------
      w(i,k) = vvm(i,k) + dt2 * w(i,k)
! --------------------------------------------------
#ifdef NONHYD
      p(i,k) = pnm(i,k) + dt2 * p(i,k)
      t(i,k) = tnm(i,k) + dt2 * t(i,k)
#endif
201   continue
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef MBC
!      print *,'mype=',mype,'lev=',levr,'lnwavs=',LNWAVS,            &
!     &    'jwav1=',jwav1
        y(1,k)=0.
        x(1,k)=0.
        w(1,k)=0.
!      remkh=mod(fhour,6.0)
!      if ( remkh .eq. 0) then
!        print *,'mype=',mype,' fhour=',fhour,'rem=',remkh
!        print *,'wave0, k=',k,'y(0,k)=',y(1,k)
!        print *,'wave0i,k=',k,'y(1,k)=',y(2,k)
!        print *,'wave1r,k=',k,'y(2,k)=',y(3,k)
!        print *,'wave1i,k=',k,'y(3,k)=',y(1+jwav1,k)
!        print *,'wave2r,k=',k,'y(4,k)=',y(2+jwav1,k)
!        print *,'wave2i,k=',k,'y(5,k)=',y(3+jwav1,k)
!        print *,'wave3r,k=',k,'y(6,k)=',y(1+jwav1*2,k)
!        print *,'wave3i,k=',k,'y(7,k)=',y(2+jwav1*2,k)
!        print *,'wave4r,k=',k,'y(8,k)=',y(3+jwav1*2,k)

!        print *,'mype=',mype,' fhour=',fhour,'rem=',remkh
!        print *,'wave0, k=',k,'x(0,k)=',x(1,k)
!        print *,'wave0i,k=',k,'x(1,k)=',x(2,k)
!        print *,'wave1r,k=',k,'x(2,k)=',x(3,k)
!        print *,'wave1i,k=',k,'x(3,k)=',x(1+jwav1,k)
!        print *,'wave2r,k=',k,'x(4,k)=',x(2+jwav1,k)
!        print *,'wave2i,k=',k,'x(5,k)=',x(3+jwav1,k)
!        print *,'wave3r,k=',k,'x(6,k)=',x(1+jwav1*2,k)
!        print *,'wave3i,k=',k,'x(7,k)=',x(2+jwav1*2,k)
!        print *,'wave4r,k=',k,'x(8,k)=',x(3+jwav1*2,k)

!        print *,'mype=',mype,' fhour=',fhour,'rem=',remkh
!        print *,'wave0, k=',k,'w(0,k)=',w(1,k)
!        print *,'wave0i,k=',k,'w(1,k)=',w(2,k)
!        print *,'wave1r,k=',k,'w(2,k)=',w(3,k)
!        print *,'wave1i,k=',k,'w(3,k)=',w(1+jwav1,k)
!        print *,'wave2r,k=',k,'w(4,k)=',w(2+jwav1,k)
!        print *,'wave2i,k=',k,'w(5,k)=',w(3+jwav1,k)
!        print *,'wave3r,k=',k,'w(6,k)=',w(1+jwav1*2,k)
!        print *,'wave3i,k=',k,'w(7,k)=',w(2+jwav1*2,k)
!        print *,'wave4r,k=',k,'w(8,k)=',w(3+jwav1*2,k)

!      endif
#endif
#ifdef MP
      endif
#endif
200   continue
#ifdef NONHYD
      do k=1,levr+1
      do i=1,LNWAVS
      o(i,k) = onm(i,k) + dt2 * o(i,k)
      enddo
      enddo
#endif
#ifdef CHGTLEV
! hmhj 04272007 force zero perturbations to u,v on model top
      k = levr
      do i=1,lnwavp
!-------------------------------------- uu ---------
      uum(i,k) = 0.0
      uu(i,k) = 0.0
      x(i,k) = 0.0
! ------------------------------------- vv ---------
      vvm(i,k) = 0.0
      vv(i,k) = 0.0
      w(i,k) = 0.0
      enddo
!-----------------------------------------------
#endif
!-----------------------------------------------
#ifdef ORIGIN_THREAD
!$doacross share(rm,dt2,rt),local(i,k)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(rm,dt2)
!mic$1 shared(rt)
!mic$1 private(i,k)
#endif
#ifdef OPENMP
!$omp parallel do private(nc,i,k)
#endif
      do 202 nc=1,ntotal
      do 202 k=1,levr
      do 205 i=1,LNWAVS
! ------------------------------------- mixing ratio
      rt(i,k,nc) = rm(i,k,nc) + dt2 * rt(i,k,nc)
205   continue
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef MBC
!      print *,'mype=',mype,'lev=',levr,'lnwavs=',LNWAVS,                 &
!     &    'jwav1=',jwav1
      if (nc .eq. 1) then
         rt(1,k,nc)=0.
      endif
!      remkh=mod(fhour,6.0)
!      if ( remkh .eq. 0) then
!      print *,'mype=',mype,' fhour=',fhour,'rem=',remkh
!      print *,'wave0, nc=',nc,'k=',k,'rt(0,k)=',rt(1,k,nc)
!      print *,'wave0i,nc=',nc,'k=',k,'rt(1,k)=',rt(2,k,nc)
!      print *,'wave1r,nc=',nc,'k=',k,'rt(2,k)=',rt(3,k,nc)
!      print *,'wave1i,nc=',nc,'k=',k,'rt(3,k)=',rt(1+jwav1,k,nc)
!      print *,'wave2r,nc=',nc,'k=',k,'rt(4,k)=',rt(2+jwav1,k,nc)
!      print *,'wave2i,nc=',nc,'k=',k,'rt(5,k)=',rt(3+jwav1,k,nc)
!      print *,'wave3r,nc=',nc,'k=',k,'rt(6,k)=',rt(1+jwav1*2,k,nc)
!      print *,'wave3i,nc=',nc,'k=',k,'rt(7,k)=',rt(2+jwav1*2,k,nc)
!      print *,'wave4r,nc=',nc,'k=',k,'rt(8,k)=',rt(3+jwav1*2,k,nc)
!      endif
#endif
#ifdef MP
      endif
#endif
202   continue
!


!!yjtest-print the mean perturbation and the waveij=0 at
!!fhour=0.00833333333333333 >
!! yj,close this output part! set fhour to -9999
!      if (abs(fhour-(-9999)) .lt. 0.001) then
!!      if (abs(fhour-0.00833333333333333) .lt. 0.001) then
!      allocate(tmp(lnwav))
!      allocate(tmp1(lnwav))
!      allocate(tmpp(lnwavp))
!      allocate(gridk1(lngrdp))
!      allocate(fullk1(igrd1,jgrd1))
!      allocate(fullk1sgl(lngrd))
!#ifdef MP
!! output the spectral ln(ps) perturbation and wave0,0 of ln(ps) perturbation
!      call rmpsp2f(z,lnwavp,tmp,lnwav,1)
!      if (mype.eq.master) then
!        print*,'size of spec ln(ps)=',iwav1,'x',jwav1,'=',lnwav
!        open(93,file='speclnps',form='unformatted')
!        call ful2sgl(tmp,tmp1,lnwav)
!        write(93) tmp1
!        tmp1=0.0
!        close(93)
!        do i=2,lnwav
!          tmp(i)=0.0
!        enddo  
!      endif        
!!tmp, tmpp, fullk1, gridk1     
!      call rmpsf2p(tmp,lnwav,tmpp,lnwavp,1)
!      call cctogd(tmpp,gridk1,1)
!      call rmpgp2f(gridk1,igrd1p,jgrd1p,fullk1(1,1),igrd1,jgrd1,1)
!      if (mype.eq.master) then
!        print*,'yjtest- wave0,0 of ln(ps)=',fullk1(1,1)
!        call ful2sgl(fullk1(1,1),fullk1sgl,lngrd)
!        open(99,file='meanlnps',form='unformatted')
!        write(99) fullk1sgl
!        close(99)
!      endif  
!! output the full value of ln(ps) perturbation
!      call cctogd(z,gridk1,1)
!      call rmpgp2f(gridk1,igrd1p,jgrd1p,fullk1(1,1),igrd1,jgrd1,1)
!      if (mype.eq.master) then
!!print the domain mean value        
!        mean=0.0
!        do i=1,lngrd
!          mean=mean+fullk1(i,1)
!        enddo  
!        print*,'yjtest- mean of ln(ps)=',mean/float(lngrd)
!!print the zonal mean
!        do j=1,jgrd1
!        mean=0.0
!        do i=1,igrd1
!          mean=mean+fullk1(i,j)
!        enddo  
!        print*,'yjtest- mean of j,ln(ps)=',j,mean/float(igrd1)
!        enddo
!!print the meridinol mean
!        do i=1,igrd1
!        mean=0.0
!        do j=1,jgrd1
!          mean=mean+fullk1(i,j)
!        enddo  
!        print*,'yjtest- mean of i,ln(ps)=',i,mean/float(jgrd1)
!        enddo
!!output file       
!        call ful2sgl(fullk1(1,1),fullk1sgl,lngrd)
!        open(98,file='fulllnps',form='unformatted')
!        write(98) fullk1sgl
!        close(98)
!      endif  
!! output the wave0,: of ln(ps) perturbation
!      call rmpsp2f(z,lnwavp,tmp,lnwav,1)
!      if (mype.eq.master) then
!        do j=jwav1+1,lnwav
!          tmp(j)=0.0
!        enddo  
!      endif  
!      call rmpsf2p(tmp,lnwav,tmpp,lnwavp,1)
!      call cctogd(tmpp,gridk1,1)
!      call rmpgp2f(gridk1,igrd1p,jgrd1p,fullk1(1,1),igrd1,jgrd1,1)
!      if (mype.eq.master) then
!        print*,'yjtest- wave0,10 of ln(ps)=',fullk1(0,10)
!        call ful2sgl(fullk1(1,1),fullk1sgl,lngrd)
!        open(97,file='meanjlnps',form='unformatted')
!        write(97) fullk1sgl
!        close(97)
!      endif  
!! output the wave:,0 of ln(ps) perturbation
!      call rmpsp2f(z,lnwavp,tmp,lnwav,1)
!      call rmpsp2f(z,lnwavp,tmp1,lnwav,1)
!      if (mype.eq.master) then
!        do i=1,iwav1
!          tmp1(1+(i-1)*jwav1)=0.0
!        enddo  
!        tmp=tmp-tmp1
!      endif  
!      call rmpsf2p(tmp,lnwav,tmpp,lnwavp,1)
!      call cctogd(tmpp,gridk1,1)
!      call rmpgp2f(gridk1,igrd1p,jgrd1p,fullk1(1,1),igrd1,jgrd1,1)
!      if (mype.eq.master) then
!        print*,'yjtest- wave10,0 of ln(ps)=',fullk1(10,0)
!        call ful2sgl(fullk1(1,1),fullk1sgl,lngrd)
!        open(96,file='meanilnps',form='unformatted')
!        write(96) fullk1sgl
!        close(96)
!      endif  
!! output the LBC meanij0  of ln(ps) perturbation
!      call rmpsp2f(z,lnwavp,tmp,lnwav,1)
!      call rmpsp2f(z,lnwavp,tmp1,lnwav,1)
!      if (mype.eq.master) then
!        do j=1,jwav1
!          tmp1(j)=0.0
!        enddo  
!        do i=1,iwav1
!          tmp1(1+(i-1)*jwav1)=0.0
!        enddo  
!        tmp=tmp-tmp1
!      endif  
!      call rmpsf2p(tmp,lnwav,tmpp,lnwavp,1)
!      call cctogd(tmpp,gridk1,1)
!      call rmpgp2f(gridk1,igrd1p,jgrd1p,fullk1(1,1),igrd1,jgrd1,1)
!      if (mype.eq.master) then
!        print*,'yjtest- wave10,0 of ln(ps)=',fullk1(10,0)
!        call ful2sgl(fullk1(1,1),fullk1sgl,lngrd)
!        open(95,file='meanijlnps',form='unformatted')
!        write(95) fullk1sgl
!        close(95)
!      print*,'yjtest- before deallocate'
!      endif  
!!output the lnps-meanij0     
!      call rmpsf2p(tmp1,lnwav,tmpp,lnwavp,1)
!      call cctogd(tmpp,gridk1,1)
!      call rmpgp2f(gridk1,igrd1p,jgrd1p,fullk1(1,1),igrd1,jgrd1,1)
!      if (mype.eq.master) then
!        print*,'yjtest- wave10,0 of ln(ps)=',fullk1(10,0)
!        call ful2sgl(fullk1(1,1),fullk1sgl,lngrd)
!        open(94,file='lnps-meanijlnps',form='unformatted')
!        write(94) fullk1sgl
!        close(94)
!      print*,'yjtest- before deallocate'
!      endif  
!#endif      
!      deallocate(tmp,tmp1,tmpp,gridk1,fullk1,fullk1sgl)
!#ifdef MP
!      if (mype.eq.master) then
!        print*,'yjtest-finish output and deallocate'
!      endif  
!#endif      
!      endif  !end if(fhour)
!!yjtest-print the mean perturbation and the waveij=0 at
!!fhour=0.00833333333333333 <




#ifdef LBC
#ifdef MP
      if ( mype.eq.master) then
#endif
      bc_mask=1.0
!yjtest >      
      bc_mask(1)=0.0
      bc_mask(2)=0.0
      bc_mask(1+jwav1)=0.0
!      do j=1,jwav1
!         bc_mask(j)=0.0
!      enddo
!      do i=1,iwav1
!         bc_mask(1+(i-1)*jwav1)=0.0
!      enddo
!yjtest <      
#ifdef MP
      endif
      call rmpsf2p(bc_mask,lnwav,bc_maskp,lnwavp,1)
#else
      bc_maskp = bc_mask
#endif
      do i=1,LNWAVS
         z(i)=z(i)*bc_maskp(i)
      enddo
      do k=1,levr
         do i=1,LNWAVS
            x(i,k)=x(i,k)*bc_maskp(i)
            y(i,k)=y(i,k)*bc_maskp(i)
            w(i,k)=w(i,k)*bc_maskp(i)
            rt(i,k,1)=rt(i,k,1)*bc_maskp(i)
         enddo
      enddo
#endif

! ----------------------------------------------------
! yj--add output the kinetic energy of each zonal wave number >
#ifdef KEout
      if (mod(fhour,6) .lt. 0.0001) then
        write(cfhour,'(I3.3)') nint(fhour)
      allocate(utmp(lnwav,levr))
      allocate(vtmp(lnwav,levr))
      allocate(ke(iwav1,jwav1))
      allocate(kezon(iwav1))
      allocate(kemer(jwav1))
      call rmpsp2f(x,lnwavp,utmp,lnwav,levr)
      call rmpsp2f(w,lnwavp,vtmp,lnwav,levr)
#ifdef MP
      if( mype.eq.master ) then
      print*,'start open the files'
        open(89,file='ke2d_'//cfhour//'.bin',form='unformatted')
        open(88,file='kezon_'//cfhour//'.txt',form='formatted')
        open(87,file='kemer_'//cfhour//'.txt',form='formatted')
      do k=1,levr
        ke=0.0
        kezon=0.0
        kemer=0.0
        do j=1,jwav1
        do i=1,iwav1
          ii=j+(i-1)*jwav1
          ke(i,j)= 0.5 *                 &
     &        ( utmp(ii,k)**2 + vtmp(ii,k)**2)     
        enddo
        enddo
        do i=1,iwav1
        do j=1,jwav1
          kezon(i)=kezon(i) + ke(i,j)
        enddo
        enddo
        do j=1,jwav1
        do i=1,iwav1
          kemer(j)=kemer(j) + ke(i,j)   
        enddo
        enddo
        write(89) ke
        write(88,*) kezon
        write(87,*) kemer
      enddo
        close(89)
        close(88)
        close(87)
      endif
#endif
      deallocate(utmp)
      deallocate(vtmp)
      deallocate(ke)
      deallocate(kezon)
      deallocate(kemer)
      endif
#endif

! yj--add output the kinetic energy of each zonal wave number <
      return
      end
