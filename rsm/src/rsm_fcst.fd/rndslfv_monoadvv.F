#include <define.h>
      subroutine rndslfv_monoadvv (grid_gr,lonlens,siin,pdot,deltim   &
#ifdef SIGP
     &                            ,psi                                &
#endif
      )
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log:
! 2011 02 20 : henry juang, initial implemented into nems as NDSL with mass_dp
! 2013 09 30 : henry juang, add option of theta advection, (used later)
!
! 2018 08 09 : ying-ju chen, changed to RSM usage
!
      use machine, only: kind_grid => kind_io8
      use mod_rndsl
#ifdef SIGP
      use sigmap_cons, only : sigma
#endif
      implicit none
#include <paramodel.h>
!
!input:
      integer, intent(in) ::  lonlens
      real(kind=kind_grid), intent(inout) :: grid_gr(igrd1p,levr)
      real(kind=kind_grid), intent(in) :: siin(levr+1)
      real(kind=kind_grid), intent(in) :: pdot(igrd1p,levr+1)
      real,   intent(in):: deltim
!local:
      logical 	lprint
      integer mass,nvars
#ifdef SIGP
      real(kind=kind_grid) :: pi(levr+1)
      real(kind=kind_grid) :: psi(igrd1p)
      real(kind=kind_grid) :: ps(igrd1p)
#endif
      integer i,k
!      integer i,n,k,lon
!
!     lprint = .false.

!     if( lprint ) then
!       print *,' enter ndslfv_advect  with monotonicity '
!     endif
!

      mass  = 0
!
! var. at n+1* 
!
#ifndef SIGP
        nvars=1
        call rvertical_cell_advect(lonlens,igrd1p,levr,nvars,          &
     &            deltim,siin,pdot(1,1),grid_gr(1,1),mass)
#else
!     input psi is surface pressure [unit: cb]
      do i=1,lonlens
        ps(i)=exp(psi(i))
      enddo
      nvars=1
      pi(levr+1)=0.0
      do i = 1, lonlens
        do k = 1, levr
          ! pi [unit: cb]
          pi(k)=0.1*(sigma(k,1)*ps(i)*10.+sigma(k,2))
        enddo
        call rvertical_cell_advect(1,1,levr,nvars,                     &
     &            deltim,pi,pdot(i,1:levr+1),grid_gr(i,1:levr),mass)
      enddo
#endif
!
! var. at n+1

!
!       if( lprint ) then
!       call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' q vertadv')
!       print *,' ------------------------------------------- '
!       print *,' finish updating n+1* in grid_gr at lan=',lan
!       endif
!

! 
! ===============================
!
      return
      end
