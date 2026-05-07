      subroutine rsfc_diag(im,ix,ps,u1,v1,t1,q1,qss,
     &                    tskin,qsurf,f10m,u10m,v10m,t2m,q2m,
     &                    prslki,evap,fm,fh,fm10,fh2)
!
      use rmachine , only : kind_phys
      use rfuncphys, only : fpvs
      use rphyscons, grav => con_g,  cp => con_cp,
     &              eps => con_eps, epsm1 => con_epsm1
      implicit none
!
      integer              im, ix
      real, dimension(ix) :: ps,   u1,   v1,   t1,  q1,  tskin,  qsurf,
     &                       f10m, u10m, v10m, t2m, q2m, prslki, evap,
     &                       fm,   fh,   fm10, fh2,qss
!
!     locals
!
      real (kind=kind_phys), parameter :: qmin=1.0e-8
      integer              k,i
!
      real(kind=kind_phys)        fhi,qssm, wrk
!     real(kind=kind_phys)        fhi, qss, wrk
!     real(kind=kind_phys) sig2k, fhi, qss
!
!     real, parameter :: g=grav
!
!     estimate sigma ** k at 2 m
!
!     sig2k = 1. - 4. * g * 2. / (cp * 280.)
!
!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals
!
!!
      do i = 1, im
        f10m(i) = fm10(i) / fm(i)
!       f10m(i) = min(f10m(i),1.)
        u10m(i) = f10m(i) * u1(i)
        v10m(i) = f10m(i) * v1(i)
        fhi     = fh2(i) / fh(i)
!       t2m(i)  = tskin(i)*(1. - fhi) + t1(i) * prslki(i) * fhi
!       sig2k   = 1. - (grav+grav) / (cp * t2m(i))
!       t2m(i)  = t2m(i) * sig2k
        wrk     = 1.0 - fhi

        t2m(i)  = tskin(i)*wrk + t1(i)*prslki(i)*fhi - (grav+grav)/cp

        if(evap(i) >= 0.) then !  for evaporation>0, use inferred qsurf to deduce q2m
          q2m(i) = qsurf(i)*wrk + max(qmin,q1(i))*fhi
        else                   !  for dew formation, use saturated q at tskin
          qss(i) = fpvs(tskin(i))
          qssm   = eps * qss(i) / (ps(i) + epsm1 * qss(i))
          q2m(i) = qssm*wrk + max(qmin,q1(i))*fhi
        endif
        qss(i) = fpvs(t2m(i))
        qssm   = eps * qss(i) / (ps(i) + epsm1 * qss(i))
        q2m(i) = min(q2m(i),qssm)
      enddo

      return
      end
