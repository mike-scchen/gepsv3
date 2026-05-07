      subroutine rndslfv_monoadvv_nonhyd (grid_gr,lonlens,siin,pdot,deltim)
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
      implicit none
#include <paramodel.h>
!yj-same as rloopmp.f >
      integer noncs, noncs2, locncs, noncos, noncos2, locncos
      integer lotgr
      integer knpy,knpx,knty,kntx,knwy,knwx
      integer knpp,knpm,kntt,kntm,knww,knwm 
#ifdef NDSLDYN
      parameter( noncs=  levr*1,noncs2=noncs*2,                          &  
#else
      parameter( noncs=1+levr*3,noncs2=noncs*2,                          &  
#endif
     &           locncs=1,                                               &
     &           knpy=locncs   ,   knpx=noncs+knpy,                      &
#ifndef NDSLDYN
     &           knty=knpy+levr,   kntx=noncs+knty,                      &
     &           knwy=knty+levr,   knwx=noncs+knwy,                      &
#endif
     &           noncos=1+levr*3,noncos2=noncos*2,                       &
     &           locncos=locncs+noncs2,                                  &
     &           knpp=locncos   ,  knpm=noncos+knpp,                    &
     &           kntt=knpp+levr,   kntm=noncos+kntt,                     &
     &           knww=kntt+levr,   knwm=noncos+knww)
      parameter(lotgr=noncs2+noncos2)
!yj-same as rloopmp.f <
!

      real(kind=kind_grid) grid_gr(igrd1p,lotgr)
      real(kind=kind_grid) pdot(igrd1p,levr+1)
      real(kind=kind_grid) siin(levr+1)
      real,   intent(in):: deltim
      integer, intent(in) :: lonlens
   
      real      qqlon(igrd1p,levr*3) !knpp,kntt,knww

      logical 	lprint

      integer mass,nvars
      integer i,n,k,kk,lon,lon_dim
      integer kpp,ktt,kww
      integer kp,kt,kw
      integer kpg,ktg,kwg
!
!     lprint = .false.

!     if( lprint ) then
!       print *,' enter ndslfv_advect  with monotonicity '
!     endif
!

      mass  = 0
!
      kpp = 1
      ktt = kpp + levr
      kww = ktt + levr
!      
      nvars = 3 
!
!!$omp parallel do schedule(dynamic,1) private(lan)
!!$omp+private(i,k,mass)
!!$omp+private(kpg,ktg,kwg)
!!$omp+private(kp ,kt ,kw )

!
!
! p t w at n+1* 
!$omp parallel do private(k,i,kp,kt,kw,kpg,ktg,kwg)
        do k=1,levr
          kp=kpp+k-1
          kt=ktt+k-1
          kw=kww+k-1
          kpg=knpp+k-1
          ktg=kntt+k-1
          kwg=knww+k-1
          do i=1,lonlens
            qqlon(i,kp) = grid_gr(i,kpg) 
            qqlon(i,kt) = grid_gr(i,ktg) 
            qqlon(i,kw) = 0.5*(grid_gr(i,kwg)+grid_gr(i,kwg+1)) 
          enddo
        enddo
!
        mass=0
        call rvertical_cell_advect (lonlens,igrd1p,levr,nvars,          &
     &            deltim,siin,pdot(1,1),qqlon(1,1),mass)
!
! p t at n+1
!$omp parallel do private(k,i,kp,kt,kpg,ktg)
        do k=1,levr
          kp=kpp+k-1
          kt=ktt+k-1
          kpg=knpp+k-1
          ktg=kntt+k-1
          do i=1,lonlens
            grid_gr(i,kpg) = qqlon(i,kp)
            grid_gr(i,ktg) = qqlon(i,kt)
          enddo
        enddo
! w at n+1, w(1) and w(levr+1) are not updated in this subroutine!
!$omp parallel do private(k,i,kw,kwg)
        do k=1,levr-1
          kw=kww+k
          kwg=knww+k
          do i=1,lonlens
!wwv            grid_gr(i,kwg) = 0.5*(qqlon(i,kw-1)+qqlon(i,kw))
          enddo
        enddo
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
