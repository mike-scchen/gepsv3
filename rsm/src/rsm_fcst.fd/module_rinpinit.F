#include <define.h>
      MODULE module_rinpinit

      use machine, only : kind_io8
      implicit none
      public
#include <inpint_paramodel.h>
!----------------------------------------------------------      
! Module contains the output of rinpinit in form of
! dynamically allocated array
! inpver, comvet, hydrid, llgz, comrinp1, comrinp2, 
! comrinp4, (comrinpsfc)
!----------------------------------------------------------     
!yj----- replace inpver ----
      REAL(kind=kind_io8),DIMENSION(:),ALLOCATABLE ::                 &
     &      siin,slin,delin,ciin,clin,rpiin 
!yj----- replace comvet common -------------
      REAL(kind=kind_io8),DIMENSION(:),ALLOCATABLE ::                 &
     &   si,sl,del,ci,cl,rpi
!yj----- replace hydrid and llgz common ----
      REAL(kind=kind_io8),DIMENSION(:,:),ALLOCATABLE ::                 &
     &      siin3d,slin3d 
      REAL(kind=kind_io8),DIMENSION(:),ALLOCATABLE ::                 &
     &      terr 
!yj----- replace comrinp1 common -------
!yj      INTEGER,DIMENSION(:),ALLOCATABLE ::              &
!yj     &      idate 
!
!      REAL(kind=kind_io4),DIMENSION(:),ALLOCATABLE :: fhour
!yj----- replace comrinp2 common -------
      REAL(kind=kind_io8),ALLOCATABLE ::                    &
     &       q(:),te(:,:),uu(:,:),vv(:,:)    &
     &      ,rq(:,:,:),gz(:)                       & 
#ifdef NONHYD
     &      ,pp(:,:),tt(:,:),ww(:,:)               &
#endif
     &      ,flat(:),flon(:),fm2(:),fm2x(:),fm2y(:)
!yj----- replace comrinp4 common ------
      REAL(kind=kind_io8),ALLOCATABLE ::                 &
     &      qo(:),teo(:,:),uuo(:,:),vvo(:,:),rqo(:,:,:)  & 
#ifdef NONHYD
     &      ,ppo(:,:),tto(:,:),wwo(:,:)                &     
#endif
     &      ,gzo(:)

!yj----- replace comrinpsfc common ---------
      REAL(kind=kind_io8),ALLOCATABLE ::                 &
     &   tsea(:),smc(:,:),sheleg(:),stc(:,:),            &
     &   tg3clm(:),zorclm(:),cv(:),cvb(:),cvt(:),        &
     &   albclm(:,:),slmsk(:),vegclm(:),canopy(:),       &
     &   f10m(:),vetclm(:),sotclm(:),alfclm(:,:),        &
     &   vmnclm(:),vmxclm(:),slpclm(:),absclm(:),        &
     &   uustar(:),ffmm(:),ffhh(:)
!
      save
!yj      
      contains

        subroutine allocate_inpver_array
          integer  ierr
          allocate (siin(levs+1),slin(levs),delin(levs), &
     &     ciin(levs+1),clin(levs),rpiin(levs-1), stat= ierr)
          if (ierr/= 0) then
            write(6,*) 'mod_rinpinit : allocate fail 1 '
            stop
          endif  
          return
        end subroutine  
        subroutine deallocate_inpver_array
          integer  ierr
!          deallocate_inpver_array
          deallocate(siin,slin,delin,ciin,clin,rpiin, stat= ierr)
          if (ierr/=0) &
     &       stop "mod_rinpinit : deallocate inpver fail"
          return
        end subroutine
!
        subroutine allocate_comvet_array
          integer  ierr
          allocate (si(levr+1),sl(levr),del(levr), &
     &     ci(levr+1),cl(levr),rpi(levr-1), stat= ierr)
          if (ierr/= 0) then
            write(6,*) 'mod_rinpinit : allocate fail comvet '
            stop
          endif  
          return
        end subroutine  
        subroutine deallocate_comvet_array(levio)
          logical levio
          integer  ierr
          if (levio) then
          deallocate(si,sl,               stat= ierr)
          else
          deallocate(      del,ci,cl,rpi, stat= ierr)
          endif
          if (ierr/=0) &
     &       stop "mod_rinpinit : deallocate comvet fail"
          return
        end subroutine
!
        subroutine allocate_hydrid_array
          integer  ierr
          allocate (siin3d(lngrd,levs+1),slin3d(lngrd,levs) &
     &     , stat= ierr)
          if (ierr/= 0) &
     &       stop 'mod_rinpinit : allocate fail hydrid'
          return
        end subroutine  
        subroutine deallocate_hydrid_array
          integer  ierr
          deallocate (siin3d,slin3d, stat= ierr)
          if (ierr/= 0) & 
     &       stop 'mod_rinpinit : deallocate fail hydrid'
          return
        end subroutine 
!
        subroutine allocate_llgz_array
          integer  ierr
          allocate (terr(lngrd), stat= ierr)
          if (ierr/= 0) &
     &       stop 'mod_rinpinit : allocate fail llgz'
          return
        end subroutine  
        subroutine deallocate_llgz_array
          integer  ierr
          deallocate (terr, stat= ierr)
          if (ierr/= 0) &
     &       stop 'mod_rinpinit : deallocate fail llgz'
          return
        end subroutine  


!yj        subroutine allocate_comrinp1_array
!yj          integer  ierr
!yj          allocate( idate(4), stat= ierr)
!yj          if (ierr/= 0) &
!yj     &       stop 'mod_rinpinit : allocate fail comrinp1'
!yj          return
!yj        end subroutine  
!yj        subroutine deallocate_comrinp1_array
!yj          integer  ierr
!yj          deallocate( idate, stat= ierr)
!yj          if (ierr/= 0) &
!yj     &       stop 'mod_rinpinit : deallocate fail comrinp1'
!yj          return
!yj        end subroutine  

        subroutine allocate_comrinp2_array
          integer  ierr
          allocate( q(lngrd),te(lngrd,levs),uu(lngrd,levs)  &
     &     ,vv(lngrd,levs),rq(lngrd,levs,ntotal),gz(lngrd)  &
#ifdef NONHYD
     &     ,pp(lngrd,levs),tt(lngrd,levs),ww(lngrd,levs+1)  &
#endif
     &     ,flat(lngrd),flon(lngrd),fm2(lngrd),fm2x(lngrd)  &
     &     ,fm2y(lngrd), stat= ierr)
          if (ierr/= 0) stop 'mod_rinpinit : allocate fail 2'
          return
        end subroutine  
        subroutine deallocate_comrinp2_array(mapio)
          logical mapio
          integer  ierr
          if (mapio) then
          deallocate(   &
     &     flat,flon,fm2,fm2x,fm2y, stat= ierr)
          else
          deallocate(q,te,uu,vv,rq,gz  &
#ifdef NONHYD
     &     ,pp,tt,ww  &
#endif
     &     , stat= ierr)
          endif

          if (ierr/=0)  &
     &       stop "mod_rinpinit : deallocate comrinp2 fail"
          return
        end subroutine


        subroutine allocate_comrinp4_array
          integer  ierr
          allocate( qo(lngrd),teo(lngrd,levr),uuo(lngrd,levr)  &
     &     ,vvo(lngrd,levr),rqo(lngrd,levr,ntotal),gzo(lngrd)  &
#ifdef NONHYD
     &     ,ppo(lngrd,levr),tto(lngrd,levr),wwo(lngrd,levr+1)  &
#endif
     &     , stat= ierr)
          if (ierr/= 0) stop 'mod_rinpinit : allocate fail comrinp4'
          return
        end subroutine  
        subroutine deallocate_comrinp4_array
          integer  ierr
          deallocate(qo,teo,uuo,vvo,rqo,gzo  &
#ifdef NONHYD
     &     ,ppo,tto,wwo  &
#endif
     &     , stat= ierr)
          if (ierr/=0)  &
     &       stop "mod_rinpinit : deallocate comrinp4 fail"
          return
        end subroutine

        subroutine allocate_comrinpsfc_array
          integer  ierr
          allocate( tsea(lngrd),smc(lngrd,lsoil)          &
     &   ,sheleg(lngrd),stc(lngrd,lsoil),tg3clm(lngrd)    &
     &   ,zorclm(lngrd),cv(lngrd),cvb(lngrd),cvt(lngrd)   &
     &   ,albclm(lngrd,4),slmsk(lngrd),vegclm(lngrd)      &
     &   ,canopy(lngrd),f10m(lngrd),vetclm(lngrd)         &
     &   ,sotclm(lngrd),alfclm(lngrd,2)                   &
     &   ,vmnclm(lngrd),vmxclm(lngrd),slpclm(lngrd)       &
     &   ,absclm(lngrd),uustar(lngrd)                     &
     &   ,ffmm(lngrd),ffhh(lngrd), stat=ierr)
          if (ierr/= 0) stop 'mod_rinpinit : allocate fail 3'
          return
        end subroutine  

        subroutine deallocate_comrinpsfc_array
          integer  ierr
!          deallocate_inpver_array
          deallocate(tsea,smc,sheleg,stc,tg3clm,zorclm    &
     &     ,cv,cvb,cvt,albclm,slmsk,vegclm,canopy,f10m    &
     &     ,vetclm,sotclm,alfclm,vmnclm,vmxclm,slpclm     &
     &     ,absclm,uustar                                 &
     &     ,ffmm,ffhh,stat=ierr)
          if (ierr/=0) &
     &       stop "mod_rinpinit : deallocate comrinpsfc fail"
          return
        end subroutine
!        
        subroutine allocate_rinpinit_array
          call allocate_inpver_array
          call allocate_comvet_array
          call allocate_hydrid_array
          call allocate_llgz_array
!yj          call allocate_comrinp1_array
          call allocate_comrinp2_array
          call allocate_comrinp4_array
          call allocate_comrinpsfc_array
          return
        end subroutine  
        subroutine deallocate_rinpinit_array
          call deallocate_inpver_array
          call deallocate_comvet_array(.false.)
          call deallocate_hydrid_array
          call deallocate_llgz_array
          call deallocate_comrinp2_array(.false.)
          return
        end subroutine  
        subroutine deallocate_rinpinitout_array
          call deallocate_comvet_array(.true.)
!yj          call deallocate_comrinp1_array
          call deallocate_comrinp2_array(.true.)
          call deallocate_comrinp4_array
          call deallocate_comrinpsfc_array
          return
        end subroutine  

        end module module_rinpinit
