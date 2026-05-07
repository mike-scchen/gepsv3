      subroutine writebase(n,qb,teb,uub,vvb,rqb,lngrdb,levs,ntotal)
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
      call rmpbp2f(qb,lngrdb,tmp2,len2,1)
      if (mype.eq.master) then
      write(n) tmp2
      endif
!--------------------------------------- te
      call rmpbp2f(teb,lngrdb,tmp3,len2,levs)
      if (mype.eq.master) then
      write(n) tmp3
      endif
!--------------------------------------- uu
      call rmpbp2f(uub,lngrdb,tmp3,len2,levs)
      if (mype.eq.master) then
      write(n) tmp3
      endif
!--------------------------------------- vv
      call rmpbp2f(vvb,lngrdb,tmp3,len2,levs)
      if (mype.eq.master) then
      write(n) tmp3
      endif
!--------------------------------------- rq
      do nc=1,ntotal
      call rmpbp2f(rqb(1,1,nc),lngrdb,tmp3,len2,levs)
      if (mype.eq.master) then
      write(n) tmp3
      endif
      enddo

      if (mype.eq.master) print*,'finish writebase'
      return
      end
