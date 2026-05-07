      subroutine rinpinit(mype,master)
!!!!!!!!!!!!!!!!!The output varialbes!!!!!!!!!!!!!!!!!!!!!!!!
! Use module_rinpinit, only: (the following vars.)
! Then, call deallocate_rinpinitout_array to deallocate them.
! comrinp1:idate
! comvet(flase):sl,si       
! comrinp2(false): fm2,fm2x,fm2y,flat,flon
! comrinp4: gzo,qo,teo,uuo,vvo,rqo,  &
!           #ifdef NONHYD
!             & ppo,tto,wwo,
!           #endif
! comrinpsfc: 
!     &   tsea(:),smc(:,:),sheleg(:),stc(:,:), &
!     &   tg3clm(:),zorclm(:),cv(:),cvb(:),cvt(:), &
!     &   albclm(:,:),slmskin(:),vegclm(:),canopy(:), &
!     &   f10m(:),vetclm(:),sotclm(:),alfclm(:,:), &
!     &   uustar(:),ffmm(:),ffhh(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! variables are save in the module_rinpinit !!!!!!!!!
!!!!!!!!! use module_rinpinit to get the output of rinpinit!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! interpolate global to regional by spherical transform in dynamic grid
! or interpolate coarse regional to regional 
! by bi-linear interpolation in physic grids.
!
! ---------------- common of rinp output -----------------
! --- comrinp1, comrinp3, comrinp4, comrinpsfc, comvet ---
! ---------------- common of rinp output -----------------
! Note #include <inpint_paramodel.h> and define ijdimi, ijdimo
!
      use funcphys
      use machine , only : kind_io4
!yj add module start -----------------------      
      use module_rinpinit      
!yj add module   end -----------------------      
!ccj
!yj#include <inpint_paramodel.h>
#include <inp_rscomcpu.h>
!
!#include <inp_rscomver.h>
#include <inp_rscomloc.h>
!
      logical         newsst,newsig,newmtn,sig2rg,sfc2rg,pgb2rg
      real plevels(levs)
      namelist/namrin/newsst,newsig,newmtn,sig2rg,sfc2rg,          &
     &                percmtn,pgb2rg
      data pgb2rg/.false./
      data newsig/.false./
      data newmtn/.false./
      data sig2rg/.false./
      data sfc2rg/.false./
      data percmtn/0.0/
      data nrslm,nrorg,nsig/   23,   24,  25/
      namelist /namlocinit/                                              &
     & rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd  &
     &,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd  &
     &,clat1,clat2,clon1,clon2
!yj
!yj      real(kind=kind_io4) fhour
!
!yj      parameter(idimi=igrd1i,jdimi=jgrd1i,ijdimi=idimi*jdimi,            &                 
!yj     &          kdimi=levri,kdimpi=kdimi+1,kdimmi=kdimi-1)
!yj      parameter(idimo=igrd1o,jdimo=jgrd1o,ijdimo=idimo*jdimo,            &
!yj     &          kdimo=levro,kdimpo=kdimo+1,kdimmo=kdimo-1)
!yj!
!yj      common /comrinp1i/idate(4)
!yj!
!yj      common/comrinp3ini/q (ijdimo),te (ijdimo,kdimi),uu (ijdimo,kdimi), &
!yj     &       vv(ijdimo,kdimi),rq(ijdimo,kdimi,ntotal),gz(ijdimo)         &
!yj#ifdef NONHYD
!yj     &      ,pp (ijdimo,kdimi),tt (ijdimo,kdimi)                         &
!yj     &      ,ww (ijdimo,kdimi+1)                                         &
!yj#endif
!yj     &      ,flato(ijdimo),flono(ijdimo)                                 &
!yj     &      ,fm2o(ijdimo),fm2xo(ijdimo),fm2yo(ijdimo)
!yj!
!yj      common/comrinp4ini/qo(ijdimo),teo(ijdimo,kdimo),uuo(ijdimo,kdimo), &
!yj     &       vvo(ijdimo,kdimo),rqo(ijdimo,kdimo,ntotal),gzo(ijdimo)      &
!yj#ifdef NONHYD
!yj     &      ,ppo(ijdimo,kdimo),tto(ijdimo,kdimo)                         &
!yj     &      ,wwo(ijdimo,kdimo+1)                                         &
!yj#endif
!yj!
!yj      common /comrinpsfcini/tsea(ijdimi),smc(ijdimi,lsoil),  &
!yj     & sheleg(ijdimi),stc(ijdimi,lsoil),tg3clm(ijdimi),zorclm(ijdimi), &
!yj     & cv(ijdimi),cvb(ijdimi),cvt(ijdimi), &
!yj     & albclm(ijdimi,4),slmskin(ijdimi),vegclm(ijdimi), &
!yj     & canopy(ijdimi), f10m(ijdimi),  &
!yj     & vetclm(ijdimi),sotclm(ijdimi),alfclm(ijdimi,2), &
!yj     & uustar(ijdimi),ffmm(ijdimi),ffhh(ijdimi)
!yj!
!yj      common/comvetini/si (kdimpo),sl (kdimo),del (kdimo),           &
!yj     &               ci(kdimpo),cl (kdimo),rpi (kdimmo)
!yj      
! enter
!yjadd
       call allocate_rinpinit_array
      if (mype.eq.master) then

      print *,'in mainrinp_init'
!      call inp_gncpus(ncpus)
!      print *,'after mainrinp_init'
!
      print *,'read namelist'
      rewind(5)
      read(5,namrin)
      print *,'after read namrin'
      read(5,namlocinit)
      print *,'after read namlocinit'
      write(6,namlocinit)
!
!
!      call gfuncphys
!
      if ( pgb2rg ) then
        write(*,*) 'levs=',levs
!yj        read(nplevels,*) plevels
        call int_l2rinp(nrslm)
!
      if( newsig .or. newmtn ) then
              print *,'call int_rchgr'
         call int_rchgr(nrslm,nrorg,nsig,      &
     &           newsig,newmtn,pgb2rg,percmtn,      &
     &           newsst)
              print *,'after call int_rchgr'
      endif
      endif
!yjadd-----the following may be marked for mpmd----
      endif   !endif of mype=master
      call deallocate_rinpinit_array
!yj----  broadcast in rmpgf2p in sread in rsmini_mpi ----      
!
      return
      end subroutine
