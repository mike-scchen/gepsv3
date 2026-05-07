      subroutine rinpbase(mype,master)
!!!!!!!!!!!!!!!!!The output varialbes!!!!!!!!!!!!!!!!!!!!!!!!
! Use module_rinpbase, only: (the following vars.)
! Then, call deallocate_rinpbaseout_array to deallocate them.
! comrinp1:idate
! comvet(flase):sl,si       
! comrinp2(false): fm2,fm2x,fm2y,flat,flon
! comrinp4: gzo,qo,teo,uuo,vvo,rqo,  &
!           #ifdef NONHYD
!             & ppo,tto,wwo,
!           #endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! variables are save in the module_rinpbase !!!!!!!!!
!!!!!!!!! use module_rinpbase to get the output of rinpbase!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! interpolate global to regional by spherical transform in dynamic grid
! or interpolate coarse regional to regional 
! by bi-linear interpolation in physic grids.
!
      use funcphys
      use machine , only : kind_io4
!yj add module start -----------------------      
      use module_rinpbase      
!yj add module   end -----------------------      

!ccj
!yj <inp_paramodel.h> is used in the module_rinpbase
!yj#include <inp_paramodel.h>
#include <inp_rscomcpu.h>
#include <inp_rscomloc.h>
!
      logical         newsst,newsig,newmtn,sig2rg,sfc2rg,pgb2rg
      namelist/namrin/newsst,newsig,newmtn,sig2rg,sfc2rg,          &
     &                percmtn,pgb2rg

      data pgb2rg/.false./
      data newsst/.false./
      data newsig/.false./
      data newmtn/.false./
      data sig2rg/.false./
      data sfc2rg/.false./
      data percmtn/0.0/
      data nrorg,nsig/   14,  15/
      namelist /namlocbase/                                              &
     & rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd  &
     &,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd  &
     &,clat1,clat2,clon1,clon2
! enter
!yjadd
      call allocate_rinpbase_array
      if (mype.eq.master) then
!yjadd
!
      print *,'in mainrinp'
      call inp_gncpus(ncpus)
      print *,'after mainrinp'
!
      print *,'read namelist'
      rewind(5)
      read(5,namrin)
      print *,'after read namrin'
      read(5,namlocbase)
      print *,'after read namlocbase'
      write(6,namlocbase)
!
      call gfuncphys
!
      if ( pgb2rg ) then
        write(*,*) 'levs=',levs
!yj surface data are not used for the mpmd
        call l2rinp
!
      if( newsig .or. newmtn ) then
              print *,'call rchgr'
         call rchgr(nrorg,nsig,      &
     &           newsig,newmtn,pgb2rg,percmtn )
              print *,'after call rchgr'
      endif
      endif
!yjadd-----the following may be marked for mpmd----
!endif of mype=master below
      endif
      call deallocate_rinpbase_array
!yj----  broadcast in cread in rsmini_mpi ----      
!yj      call rmpbcasti(idate,4)
!yj      call rmpbcastr(sl,levr)
!yj      call rmpbcastr(si,levr+1)
!yj      call rmpbcastr(fm2,lngrd)
!yj      call rmpbcastr(fm2x,lngrd)
!yj      call rmpbcastr(fm2y,lngrd)
!yj      call rmpbcastr(flat,lngrd)
!yj      call rmpbcastr(flon,lngrd)
!yj      call rmpbcastr(gzo,lngrd)
!yj      call rmpbcastr(qo,lngrd)
!yj      call rmpbcastr(teo,lngrd*levr)
!yj      call rmpbcastr(uuo,lngrd*levr)
!yj      call rmpbcastr(vvo,lngrd*levr)
!yj      call rmpbcastr(rqo,lngrd*levr*ntotal)
!yj#ifdef NONHYD
!yj      call rmpbcastr(ppo,lngrd*levr)
!yj      call rmpbcastr(tto,lngrd*levr)
!yj      call rmpbcastr(wwo,lngrd*levr)
!yj#endif
! comrinp1:idate
! comvet(flase):sl,si       
! comrinp2(false): fm2,fm2x,fm2y,flat,flon
! comrinp4: gzo,qo,teo,uuo,vvo,rqo,  &
!           #ifdef NONHYD
!             & ppo,tto,wwo,
!           #endif
!yjadd
!
      return
      end subroutine
