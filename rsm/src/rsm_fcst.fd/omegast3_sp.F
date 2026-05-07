      subroutine omegast3_sp(im,ix,nz,                                   &   
     &      dphi,dlam,ug,vg,dg,rcl,vvel,ps,si,bki)
!fpp$ noconcur r
!....   code lifted from post (mcp1840) jun 88--computes vvel (cb/sec)
!....    input ps in cb,output vvel in cb/sec
!....   do loops altered for better vectorization possibilities..k.a.c.
!c
      use machine , only : kind_phys
      implicit none
!!
      integer              i,k,le,nz,ix,im
      real(kind=kind_phys) rcl
!c
      real(kind=kind_phys) dphi(ix),dlam(ix),ps(ix)
      real(kind=kind_phys)                                               &  
     & cg(ix,nz),ug(ix,nz),vg(ix,nz),                                    &
     & dg(ix,nz),si(ix,nz+1)
!...   vvel contains omega in layers on return from subroutine...
      real(kind=kind_phys) vvel(ix,nz)
      real(kind=kind_phys) SUM(ix)
      real(kind=kind_phys) bki(nz+1)
      real(kind=kind_phys) VGRADP(ix),GRADPV(ix)
!!
!...  compute c=v(true)*del(ln(ps)).divide by cos for del cos for v
      do 48 i=1,im
        dphi(i)=dphi(i)*rcl
        dlam(i)=dlam(i)*rcl
   48 continue
      SUM=0.
      do 5 le=nz,1,-1
        do 50 i=1,im
          VGRADP(i)=(ug(i,le)*dlam(i)+vg(i,le)*dphi(i))*(bki(le)-bki(le+1))
          GRADPV(i)=VGRADP(i)+dg(i,le)*(si(i,le)-si(i,le+1))
          SUM(i)=SUM(i)+GRADPV(i)
          vvel(i,le)=( (0.5*dg(i,le)*(si(i,le)-si(i,le+1)))             &
     &            +(bki(le)*VGRADP(i)/(bki(le)-bki(le+1)))                 &
     &            -SUM(i) )*ps(i)
   50 continue
    5 continue
!!
      return
      end
