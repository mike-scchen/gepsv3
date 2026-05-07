      subroutine rndslfv_monoadvh2                                     &
     &  (grid_gr,levrp,k1,k2,lotgr,uulon,vvlon,deltim,rdelx,rdely,difhsec)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
      use mod_rndsl
      implicit none
#include <define.h>
#include <npesi.h>
!!#include <comrmpi.h>

      
!for mpi
!yj-same as rloopmp.f >
!yj      integer kumk,kvmk,ktmk,kpmk,kqmk,kgzmk
!yj#ifdef NONHYD
!yj      integer kpnmk,ktnmk,kwmk
!yj#endif
!NDSL:synmk>
!yj      parameter(kumk=1,                                                  &
!yj     &          kvmk=levrp+kumk,                                         &
!yj     &          ktmk=levrp+kvmk,                                         &
!yj     &          kpmk=levrp+ktmk,                                         &
!yj     &          kqmk=levrp+kpmk,                                         & 
!yj     &          kgzmk=levhp+kqmk                                         &
!yj#ifdef NONHYD     
!yj     &         ,kpnmk=levrp+kgzmk                                        &
!yj     &         ,ktnmk=levrp+kpnmk                                        &
!yj     &         ,kwmk=levrp+ktnmk !kwmt:levrp+1 (midpoint: 1~levrp+1)     &
!yj#endif
!yj     &          )
!NDSL:synmk<
!yj#ifdef NONHYD
!yj      parameter(lotgr=1+8*levrp+levhp)
!yj#else
!yj      parameter(lotgr=  5*levrp+levhp)
!yj#endif
!yj-same as rloopmp.f <
      
      integer, intent(in):: levrp,k1,k2,lotgr
      real, intent(in):: rdelx, rdely, difhsec
      integer:: nkk,nvar

      real            :: grid_gr(lonfull,lotgr,latpart)
      real, intent(in):: deltim
   
      real, intent(in):: uulon(lonfull,levrp,latpart)
      real, intent(in):: vvlon(lonfull,levrp,latpart)
      real,allocatable:: qqlon(:,:,:),rrlon(:,:,:)

      real      vvlat(latfull,levrp,lonpart)
      real,allocatable::qqlat(:,:,:),rrlat(:,:,:)
      real      rma, rm2a, rdt2, rkt, pi, cons0, cons1

      logical   lprint

      integer mono,mass
      integer ilan,i,j,n,k,kk,lon,lan,lat,lons_lat,jlonf,irc
!
       nkk=k2-k1+1
! this will do k1 ~ k2  !(k2 include)       
       nvar=(nkk-1)/levrp+1
      allocate(                         &
     &  qqlon(lonfull,nkk,latpart),   & 
     &  rrlon(lonfull,nkk,latpart),   &
     &  qqlat(latfull,nkk,lonpart),   &
     &  rrlat(latfull,nkk,lonpart)    &
     &  )
!      
      lprint = .false.

!yj      print *,' enter ndslfv_monoadvh '
!yj      if( lprint ) print *,' enter ndslfv_monoadvh '
!
      mono  = 1
      mass  = 0
!
      rdt2 = 0.5 / deltim      
!
! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
!
#ifdef OPENMP
!$omp parallel do private(lan)  &
!$omp schedule(dynamic)
#endif
      do lan=1,mylatlen
!
        qqlon(:,:,lan) = grid_gr(:,k1:k2,lan)

        rrlon(:,:,lan) = qqlon(:,:,lan)
!
! first set positive advection in east-west direction 
!
        call rsopenend_cell_massadv(lonfull,levrp,nvar,deltim   &
     &               ,uulon(1,1,lan),rrlon(1,1,lan),mass,rdelx)
#ifdef dediff
        do k=1,nkk
        call dediffusion(lonfull,rrlon(1,k,lan),rdelx,difhsec)
        enddo
#endif

      enddo  !end do lan=1,mylatlen
#ifdef OPENMP
!$omp end parallel do
#endif

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat

       if( lprint ) print *,' ndslfv_advect transport from we to ns '

       call rpara_we2ns(vvlon,vvlat,levrp)
       call rpara_we2ns(rrlon,rrlat,nkk)
       call rpara_we2ns(qqlon,qqlat,nkk)

       if( lprint ) then
       print *,' ------------ after we2ns ---------------------- '
       do lon=1,mylonlen
        print *,'  lon=',lon
!yj        call mymaxmin(vvlat(1,1  ,lon),latfull,latfull,1,' we2ns v')
!yj        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' we2ns r')
!yj        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' we2ns q')
       enddo
       endif
!
! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

       if( lprint ) then
       print *,' ndslfv_advect adv loop in y '
       print *,' mylonlen=',mylonlen
       endif

#ifdef OPENMP
!$omp parallel do private(lon)  &
!$omp schedule(dynamic) 
#endif
       do lon=1,mylonlen
! 
        if( lprint ) print *,' lon=',lon

! first set advection in north-south direction in great circle through two poles
!
        call rsopenend_cell_massadv(latfull,levrp,nvar,deltim   &
     &               ,vvlat(1,1,lon),rrlat(1,1,lon),mass,rdely)
#ifdef dediff
        do k=1,nkk
        call dediffusion(latfull,rrlat(1,k,lon),rdely,difhsec)
        enddo
#endif

        if( lprint ) then
        print *,' ------------------------------------------- '
!yj        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry u ')
!yj        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry v ')
!yj        call mymaxmin(rrlat(1,ktt,lon),latfull,latfull,1,' advry t ')
!yj        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' advry q ')
!yj        call mymaxmin(rrlat(1,kdp,lon),latfull,latfull,1,' advry dp ')
!yj        call mymaxmin(rrlat(1,kp2,lon),latfull,latfull,1,' advry p2 ')
        endif
!
! second set advection in north-south direction in great circle through two poles
!
        call rsopenend_cell_massadv(latfull,levrp,nvar,deltim   &
     &               ,vvlat(1,1,lon),qqlat(1,1,lon),mass,rdely)
#ifdef dediff
        do k=1,nkk
        call dediffusion(latfull,qqlat(1,k,lon),rdely,difhsec)
        enddo
#endif

        if( lprint ) then
        print *,' ------------------------------------------- '
!yj        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy u ')
!yj        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy v ')
!yj        call mymaxmin(qqlat(1,ktt,lon),latfull,latfull,1,' advqy t ')
!yj        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' advqy q ')
!yj        call mymaxmin(qqlat(1,kdp,lon),latfull,latfull,1,' advqy dp ')
!yj        call mymaxmin(qqlat(1,kp2,lon),latfull,latfull,1,' advqy p2 ')
        print *,' done with y at lon=',lon
        endif

       enddo  !lon=1,mylonlen
#ifdef OPENMP
!$omp end parallel do
#endif
!
! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory 
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon

       if( lprint ) print *,' ndslfv_advect transport from ns to we '

       call rpara_ns2we(rrlat,rrlon,nkk)
       call rpara_ns2we(qqlat,qqlon,nkk)

       if(lprint) then
       print *,' ------------ after ns2we ---------------------- '
!yj       do lan=1,lats_node_a
!yj        print *,'  lan=',lan
!yj        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' ns2we r')
!yj        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' ns2we q')
!yj       enddo
       endif

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!      print *,' ndslfv_advect adv loop in x for last '

#ifdef OPENMP
!$omp parallel do private(lan,i,k)  & 
!$omp schedule(dynamic) 
#endif
      do lan=1,mylatlen

! second set advection in x for the second of the pair
        call rsopenend_cell_massadv(lonfull,levrp,nvar,deltim   &
     &               ,uulon(1,1,lan),qqlon(1,1,lan),mass,rdelx)
#ifdef dediff
        do k=1,nkk
        call dediffusion(lonfull,qqlat(1,k,lan),rdelx,difhsec)
        enddo
#endif

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!yj        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' adv x q ')
!        endif

        do k=1,nkk
          do i=1,lonfull
            qqlon(i,k,lan) = 0.5 * ( qqlon(i,k,lan) + rrlon(i,k,lan) )
          enddo
        enddo

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!yj        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu u ')
!yj        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu v ')
!yj        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' do redu t ')
!yj        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' do redu q ')
!yj        call mymaxmin(qqlon(1,kdp,lan),lons_lat,lonfull,1,' do redu dp')
!yj        call mymaxmin(qqlon(1,kp2,lan),lons_lat,lonfull,1,' do redu p2')
!        print *,' finish horizonatal advection at lan=',lan
!        endif
! tendency at n
        grid_gr(:,k1:k2,lan) = (qqlon(:,:,lan)-grid_gr(:,k1:k2,lan))*rdt2
   
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        print *,' finish horizonatal advection at lan=',lan
!        endif
!
      enddo  ! do lan=1, mylatlen
#ifdef OPENMP
!$omp end parallel do
#endif

! 
! ===============================
      deallocate(qqlon,rrlon,qqlat,rrlat)
!
      return
      end
