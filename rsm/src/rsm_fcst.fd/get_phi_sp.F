      subroutine get_phi_sp(im,ix,levs,rd,nu,t,q,prsi,phii,phil)
!
      use machine , only : kind_phys
      implicit none
!
      integer im, ix, levs
      real(kind=kind_phys) rd, nu
      real(kind=kind_phys) prsi(ix,levs+1)                              &
     &,                    phii(ix,levs+1), phil(ix,levs)               &
     &,                    t(ix,levs),      q(ix,levs)
      real(kind=kind_phys) tem, dphib, dphit
      integer i, l
!
      do i=1,im
        phii(i,1)   = 0.0           ! ignoring topography height here
      enddo
      do l=1,levs
        do i=1,im
          dphib       = -rd * t(i,l) * (1.0 + nu * max(q(i,l),0.0))     &          
     &                * (prsi(i,l+1)-prsi(i,l))/(prsi(i,l+1)+prsi(i,l))
          dphit       = dphib
          phil(i,l)   = phii(i,l) + dphib
          phii(i,l+1) = phil(i,l) + dphit
        enddo
      enddo
!
      return
      end
