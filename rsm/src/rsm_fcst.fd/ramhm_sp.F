      subroutine ramhm_sp(ps0,sigma,tov,sv,p1,p2,am,hm)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  ramhm_sp
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!   prgmmr:  ying-ju chen               org: cwb        date: 22-02-17
!
! abstract:  prepare some constant arraies for dynamic computation and semi-
!            implicit computation.
!
! program history log:
!   22-02-17 modified for sigma-P coordinate [yjchen]
!
! usage:    call ramhm_sp(ps0,sigma,tov,sv,p1,p2,am,hm)
!   input argument list:
!     ps0 - reference surface pressure [unit: cb]
!     sigma - interface values, sigma1 = Bk; sigma2 = Ak; P[mb]=Ak+Bk*Ps[mb]
!           (bottom:1, top:levr+1)
!     tov - reference temperature [unit: K]
!
!   output argument list:
!     sv     - m0^2/Ps0* dp 
!     p1     - inverse of rpi
!     p2     - rpi
!     am     - matrix to compute geopotential height
!     hm     - matrix coeff. of lnPs
!
!   input files: none
!
!   output files: none
!
!   subprograms called:
! rprnta
! rprntb
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
      use physcons, cp => con_cp, rd => con_rd
!ccj
#include <paramodel.h>
#include <rscommap.h>
!cc
      real ps0
      real tov(levr)
      dimension sigma(levrp1,2),sv(levr),am(levr,levr),hm(levr,levr),   &
     &          p1(levr),p2(levr),rpi(levrm1),                          &
     &          dsi(levr),sl2(levr),p0i(levrp1),dp(levr),rp2(levr)
      real pak,rk
      integer k,le,j

!     print *,' ----- start ramhm_sp ------ '
!     print*,sigma(:,1)
! initialization
      do k=1,levr
        dsi(k)=sigma(k,1)-sigma(k+1,1)
        sl2(k)=sigma(k,1)+sigma(k+1,1)
      enddo
! reference pressure levels p0i [unit: cb]
!     p0i(levrp1)=0.0, sigma(levr+1,:) are zeros by setting
      do k=1,levrp1
        p0i(k)=sigma(k,1)*ps0*10.+sigma(k,2)
        p0i(k)=p0i(k)*0.1
      enddo
!
      do k=1,levr
        dp(k)=p0i(k)-p0i(k+1)
        rp2(k)=1.0/(p0i(k)+p0i(k+1))
      enddo
! sv
      do le=1,levr
      sv(le) = dp(le)/ps0 * xm2m
      enddo
!
! rpi(k) = ( pl(k+1) / pl(k) ) ** kapa
      rk=rd/cp
      do k=1,levrm1
        rpi(k)=( (p0i(k)+p0i(k+1))/(p0i(k+1)+p0i(k+2)) ) ** rk
      enddo
!
      do le=1,levrm1
      p1(le+1) = 1.e0 / rpi(le)
      p2(le) = rpi(le)
      enddo
      p1(1)    = 0.e0
      p2(levr) = 0.e0
! am
      am(:,:)=0.0
      do k=1,levr
      am(k,k) = rd * dp(k)*rp2(k)
      enddo
      do k=2,levr
      do j=1,k-1
      am(k,j) = 2.0 * rd * dp(j)*rp2(j)
      enddo
      enddo
! hm
      do k=1,levr-1
        pak = rd*ps0*tov(k)*rp2(k)*(dsi(k)-dp(k)*sl2(k)*rp2(k))
        hm(k  ,1) = hm(k,1) + pak
        hm(k+1,1) = hm(k,1) + pak
      enddo
      k=levr
        pak = rd*ps0*tov(k)*rp2(k)*(dsi(k)-dp(k)*sl2(k)*rp2(k))
        hm(k  ,1) = hm(k,1) + pak
      do k=1,levr
        hm(k,1) = hm(k,1) + rd*ps0*tov(k)*sl2(k)*rp2(k)
      enddo

      return
      end
