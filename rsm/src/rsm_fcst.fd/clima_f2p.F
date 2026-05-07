      subroutine clima_f2p
!input
        use module_gfsdata, only :                                      &      
     &             tg3clm,zorclm,albclm,vegclm,vetclm,sotclm,alfclm,     &
     &             vmnclm,vmxclm,slpclm,absclm
!output are stored in rscomgrd.h comm/rcomgd4/
#include <paramodel.h>
#include <rscomgrd.h>
!
      real, allocatable :: spec(:)
      allocate (spec(lngrd))

!----------------------------------------------tg3
      call rmpgf2p(tg3clm,igrd1,jgrd1,tg3,igrd1p,jgrd1p,1)
!----------------------------------------------zorl
      call rmpgf2p(zorclm,igrd1,jgrd1,zorl,igrd1p,jgrd1p,1)
!----------------------------------------------albedo
      call rmpgf2p(albclm(:,1),igrd1,jgrd1,alvsf,igrd1p,jgrd1p,1)
      call rmpgf2p(albclm(:,2),igrd1,jgrd1,alvwf,igrd1p,jgrd1p,1)
      call rmpgf2p(albclm(:,3),igrd1,jgrd1,alnsf,igrd1p,jgrd1p,1)
      call rmpgf2p(albclm(:,4),igrd1,jgrd1,alnwf,igrd1p,jgrd1p,1)
!----------------------------------------------vfrac
      call rmpgf2p(vegclm,igrd1,jgrd1,vfrac,igrd1p,jgrd1p,1)
!----------------------------------------------vfrac
      call rmpgf2p(vetclm,igrd1,jgrd1,vtype,igrd1p,jgrd1p,1)
!----------------------------------------------stype
      call rmpgf2p(sotclm,igrd1,jgrd1,stype,igrd1p,jgrd1p,1)
!----------------------------------------------facswf
      call rmpgf2p(alfclm(:,1),igrd1,jgrd1,facsf,igrd1p,jgrd1p,1)
      call rmpgf2p(alfclm(:,2),igrd1,jgrd1,facwf,igrd1p,jgrd1p,1)
!----------------------------------------------shdmin
      call rmpgf2p(vmnclm,igrd1,jgrd1,shdmin,igrd1p,jgrd1p,1)
!----------------------------------------------shdmax
      call rmpgf2p(vmxclm,igrd1,jgrd1,shdmax,igrd1p,jgrd1p,1)
!----------------------------------------------slopetyp
      call rmpgf2p(slpclm,igrd1,jgrd1,slope,igrd1p,jgrd1p,1)
!----------------------------------------------snoalb
      call rmpgf2p(absclm,igrd1,jgrd1,snoalb,igrd1p,jgrd1p,1)
!---
      return
      end subroutine clima_f2p
