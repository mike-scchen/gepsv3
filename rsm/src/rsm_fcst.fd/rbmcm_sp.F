! --------------------- rver ---------------------------
! modified by hann-ming henry juang for rsm.
! ------------------------------------------------------
      subroutine rbmcm_sp(ps0,sigma,tov,sv,p1,p2,am,hm,h1,h2,bm,cm)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rbmcm_sp
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!   prgmmr:  ying-ju chen               org: cwb        date: 22-02-23
!
! abstract:  prepare some constant arraies for dynamic computation and semi-
!            implicit computation.
!
! program history log:
!   22-02-23 modified for sigma-P coordinate [yjchen]
!
! usage:    call rbmcm_sp(ps0,sigma,tov,sv,p1,p2,am,hm,h1,h2,bm,cm)
!   input argument list:
!     ps0 - reference surface pressure [unit: cb]
!     sigma - interface values, sigma1 = Bk; sigma2 = Ak; P[mb]=Ak+Bk*Ps[mb]
!           (bottom:1, top:levr+1)
!     tov - reference temperature [unit: K]
!     sv     - m0^2/Ps0* dp  (del time mean of map factor)
!     p1     - inverse of rpi
!     p2     - rpi
!     am     - matrix to compute geopotential height
!     hm     - matrix coeff. of lnPs
!
!   output argument list:
!     h1     - vector for temperature equation
!     h2     - as h1
!     bm     - for semi-implicit computation
!     cm     - for semi-implicit computation
!
!   local variables list:
!     ci     - sigma interface values
!
!   common block :
!     rscommap
!
!   input files: none
!
!   output files: none
!
!   subprograms called:
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
      use physcons, rocp => con_rocp
!ccj
#include <paramodel.h>
#include <rscommap.h>
!cc
      real ps0
      real tov(levr)
      dimension sigma(levrp1,2),sv(levr),am(levr,levr),hm(levr,levr),   &
     &          p1(levr),p2(levr),rpi(levrm1),                          &
     &          p0i(levrp1),dp(levr)
      real rk
      integer k,j
      dimension ci(levrp1),pak(levr),x1(levr),x2(levr)
      dimension bm(levr,levr), cm(levr,levr), h1(levr), h2(levr)
!      print *,' ---- start bmcm ---- '
!
! reference pressure levels p0i [unit: cb]
!     p0i(levrp1)=0.0, sigma(levr+1,:) are zeros by setting
      do k=1,levrp1
        p0i(k)=sigma(k,1)*ps0*10.+sigma(k,2)
        p0i(k)=p0i(k)*0.1
         ci(k)=sigma(k,1)-1.0
      enddo
!
      do k=1,levr
        dp(k)=p0i(k)-p0i(k+1)
        pak(k)=(sigma(k,1)+sigma(k+1,1))/(p0i(k)+p0i(k+1))
      enddo
! bm
      rk = rocp
      do k=1,levrm1
      h1(k+1) = ( p1(k+1) * tov(k) - tov(k+1))/(2.0*dp(k+1))
      enddo
      h1(1)=0.e0

      do k=1,levrm1
      h2(k) = ( tov(k) - p2(k) * tov(k+1))/(2.0*dp(k))
      enddo
      h2(levr)=0.e0

      do k=1,levr
      x1(k) = rk*tov(k)*pak(k)+ci(k)*h1(k)+ci(k+1)*h2(k)
      x2(k) = h1(k)+h2(k)
      enddo

      bm(:,:)=0.e0
      do k=1,levr
      do j=1,levr
      bm(k,j) = x1(k)*dp(j)
      enddo
      enddo

      loop_k: do k=1,levr
      loop_j: do j=1,levr
          if(j.gt.k) exit loop_j
          bm(k,j) = bm(k,j) + x2(k)*dp(j)
        enddo loop_j
      enddo loop_k
      do k=1,levr
      bm(k,k) = bm(k,k) - h1(k) * dp(k)
      enddo

      do j=1,levr
      do i=1,levr
      bm(i,j) = bm(i,j) * xm2m
      enddo
      enddo
! cm
      do i=1,levr
      do j=1,levr
      cm(i,j) = 0.e0
      do k=1,levr
      cm(i,j) = cm(i,j) + am(i,k) * bm(k,j)
      enddo
      cm(i,j) = cm(i,j) + hm(i,1) * sv(j)
      enddo
      enddo
! cm, to co-use rgetcd.F, let cm be -cm
      do i=1,levr
      do j=1,levr
      cm(i,j) = -1.0*cm(i,j)
      enddo
      enddo
!cc   print *,' ---- end bmcm ---- '
      return
      end
