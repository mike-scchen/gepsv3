      subroutine rayleifr(ix,im,lev,levtop,ut,vt)

      implicit  none
      integer   ix,im,lev

      real      ut(ix,lev),vt(ix,lev),xm(ix)
!      real      wmax(lev)
!
!     integer, parameter :: levtop=6,lev2=3
      integer levtop,lev2

      real     ckdx,day,windmax
      data ckdx/40./, day/1./

!      data windmax/60./

!     logical dofric

      integer k,i   
      real    ckdy,frictime,fac,fricd,cdx,cdy

      lev2=levtop/2
!
!      wmax(1:lev)=0.
!      do 10 k=1,lev
!      do 10 i=1,im
!        wmax(k)= max(wmax(k),xm(i)*sqrt(ut(i,k)**2+vt(i,k)**2))
! 10   continue
!     call mpe_global_max(wmax,lev,mpe_double)
!     dofric=.false.
!     do k=levr-levtop,lev
!       if( wmax(k) .gt. windmax )then
!         dofric=.true.
!         if(myrank.eq.0)print *,' rayleifr k=',k,' wmax=',wmax(k), &
!          ' dofric=',dofric
!       endif
!     enddo
!     if(myrank.eq.0)print *,' rayleifr  dofric=',dofric
!
      ckdy=ckdx
      frictime=1./(day*86400.)
!     if(dofric)then
        do k=lev-levtop,lev
          fac=1.+exp(-1.0*((lev-k)-lev2)/lev2)
          if(fac.gt.3.)fac=3.
          fricd=frictime*fac
          cdx=ckdx*fricd
          cdy=ckdy*fricd
          do i = 1, im 
            ut(i,k)=ut(i,k)-cdx*ut(i,k)
            vt(i,k)=vt(i,k)-cdy*vt(i,k)
          enddo
        enddo
!     endif
!
      return
      end
