      subroutine readbase(n,qb,teb,uub,vvb,rqb,lngrdb,levs,ntotal)
#include <defopt.h>
#include <npesi.h>
#include <comrmpi.h>
      integer lngrdb,levs
      dimension qb(lngrdb),uub(lngrdb,levs),vvb(lngrdb,levs)           &
     &  ,teb(lngrdb,levs),rqb(lngrdb,levs,ntotal)
      dimension tmp2(lngrdb*npes),tmp3(lngrdb*levs*npes)
      integer len2
      len2=lngrdb*npes

!---------------------------------------- q
      if (mype.eq.master) then
      rewind(n)
      read(n) tmp2
      endif
      call rmpbf2p(tmp2,len2,qb,lngrdb,1)
!--------------------------------------- te
      if (mype.eq.master) then
      read(n) tmp3
      endif
      call rmpbf2p(tmp3,len2,teb,lngrdb,levs)
!--------------------------------------- uu
      if (mype.eq.master) then
      read(n) tmp3
      endif
      call rmpbf2p(tmp3,len2,uub,lngrdb,levs)
!--------------------------------------- vv
      if (mype.eq.master) then
      read(n) tmp3
      endif
      call rmpbf2p(tmp3,len2,vvb,lngrdb,levs)
!--------------------------------------- rq
      do nc=1,ntotal
      if (mype.eq.master) then
      read(n) tmp3
      endif
      call rmpbf2p(tmp3,len2,rqb(1,1,nc),lngrdb,levs)
      enddo

      if (mype.eq.master) print*,'finish readbase'
      return
      end
