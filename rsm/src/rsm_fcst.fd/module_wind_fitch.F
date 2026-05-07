!WRF:MODEL_LAYER:PHYSICS
#include <defopt.h>

MODULE module_wind_fitch
!
!Represents kinetic energy extracted by wind turbines and turbulence
! (TKE) they produce at model levels within the rotor area.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NOTICE
!! The following paper should be cited whenever presenting results using this scheme
!! (using either the original version or any modified versions of the scheme):
!! Fitch, A. C. et al. 2012: Local and Mesoscale Impacts of Wind Farms as Parameterized in a
!! Mesoscale NWP Model. Monthly Weather Review, doi:http://dx.doi.org/10.1175/MWR-D-11-00352.1
!!
!! Anna C. Fitch, National Center for Atmospheric Research (formerly University of Bergen)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! History of changes:
!
! WRFV3.5.1: 
! WRFV3.6:   Modified by Pedro A. Jimenez to include:
!             - Initialize the wind turbines in this module.
!             - Introduce z_at_walls to avoid instabilities due to neglecting
!                the perturbation of the geopotential height.
!             - User friendly interface to introduce the technical characteritics of
!                the wind turbines.
!             - Only uses one set of turbine coefficients using the wind speed at hub height
!             - Two standing coefficients.
!             - Calculates the power produced by the wind turbines.
!
! References:
!
! Fitch, A. C. et al. 2012: Local and Mesoscale Impacts of Wind Farms as Parameterized in a
!    Mesoscale NWP Model. Monthly Weather Review, doi:http://dx.doi.org/10.1175/MWR-D-11-00352.1
! Fitch, A. C. et al. 2013: Mesoscale Influences of Wind Farms Throughout a Diurnal Cycle.
!    Monthly Weather Review, doi:http://dx.doi.org/10.1175/MWR-D-12-00185.1
! Fitch, A. C. et al. 2013: Parameterization of Wind Farms in Climate Models.
!    Journal of Climate, doi:http://dx.doi.org/10.1175/JCLI-D-12-00376.1
! Jimenez, P.A., J. Navarro, A.M. Palomares and J. Dudhia:  Mesoscale modeling of offshore wind turbines
!    wakes at the wind farm resolving scale: a composite-based analysis with the WRF model over Horns Rev.
!    Wind Energy, (In Press.).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  USE module_driver_constants, ONLY : max_domains
!  USE module_model_constants, ONLY :  piconst
!
!  USE module_llxy
!  USE module_dm, ONLY : wrf_dm_min_real
!  USE module_configure, ONLY : grid_config_rec_type

#include <define.h>

  IMPLICIT NONE
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
#endif

  INTEGER, PARAMETER :: MAXVALS  = 100
  INTEGER, PARAMETER :: MAXVALS2 = 100

  real,parameter     :: piconst  = 3.1415926535897931 ! pi
!
  INTEGER           :: fi_nt = 0
  INTEGER, DIMENSION(:), ALLOCATABLE :: NKIND
  INTEGER, DIMENSION(:), ALLOCATABLE :: wival,wjval
  REAL, DIMENSION(:), ALLOCATABLE :: hubheight,diameter,stc,stc2,cutin,cutout,npower,power,spdh
  INTEGER, DIMENSION(:), ALLOCATABLE :: wim,wi
!
  REAL :: turbws(maxvals,maxvals2),turbtc(maxvals,maxvals2),turbpw(maxvals,maxvals2)
!
   logical :: windfarm_initialized = .false.
!
CONTAINS

  SUBROUTINE  dragforce(                               &
       & phii,phil,u,v,dx,dt,du,dv,ix,im,km,kt    )
!       & id,z_at_w,u,v,dx,dz,dt,qke,du,dv        &
!       &,windfarm_opt,power                      &
!       &,ids,ide,jds,jde,kds,kde                 &
!       &,ims,ime,jms,jme,kms,kme                 &
!       &,its,ite,jts,jte,kts,kte                 &
!       &)
!
!
      use physcons, g=>con_g

  implicit none

!
!  INTEGER, INTENT(IN) :: id,windfarm_opt 
!  INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
!  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
!  INTEGER, INTENT(IN) :: ids,ide,jds,jde,kds,kde
  REAL, INTENT(IN) :: dx,dt
  INTEGER, INTENT(IN) :: ix,im,km
  REAL, DIMENSION(ix,km), INTENT(IN) :: phil,u,v
  REAL, DIMENSION(ix,km) :: dz,z_at_w
  REAL, DIMENSION(ix,km+1), INTENT(IN) :: phii
!  REAL, DIMENSION(ix,km), INTENT(INOUT) :: du,dv,qke
  REAL, DIMENSION(ix,km), INTENT(INOUT) :: du,dv
!
! Local
!
  REAL     blade_l_point,blade_u_point,zheightl,zheightu,z1,z2,tarea
  REAL     speed,tkecof,powcof,thrcof,wfdensity,z8w(im,km+1)
  real     conw
!  INTEGER  itf,jtf,ktf
  INTEGER  i,j,k,n
  INTEGER  k_turbine_bot, k_turbine_top

  LOGICAL :: kfound
!
! ... PAJ: more variables ...
!
  REAL :: speedhub,speed1,speed2
  real :: power1,power2,area,ec
  INTEGER :: kbot,ktop,kt

  conw=1./g
  do k = 1,km
  do i = 1,im
      z8w(i,k) = phii(i,k)*conw
   z_at_w(i,k) = phil(i,k)*conw
  enddo
  enddo
  do i = 1,im
      z8w(i,km+1) = phii(i,km+1)*conw
  enddo
  do k = 1,km
  do i = 1,im
   dz(i,k) = z8w(i,k+1)-z8w(i,k)
  enddo
  enddo

!  itf=MIN0(ite,ide-1)
!  jtf=MIN0(jte,jde-1)
!  ktf=MIN0(kte,kde-1)

    wfdensity = 1.0/(dx*dx)   !  per turbine, so numerator is 1

!    DO kt = 1,nt
!      IF ( windfarm_opt .eq. 1 ) THEN
!
! vertical layers cut by turbine blades
!
        k_turbine_bot=0      !bottom level
        k_turbine_top=-1     !top level
        i = wival(kt)
!        j = wjval(kt)
!
         if (i .ne. -9999) then
!         if (i.ne.-9999.and.j.ne.-9999) then
         IF (i .le. im ) then
!        IF (( its .LE. i .AND. i .LE. itf ) .AND. &
!            ( jts .LE. j .AND. j .LE. jtf )  ) THEN
!
          blade_l_point=hubheight(kt)-diameter(kt)/2. ! height of lower blade tip above ground (m)
          blade_u_point=hubheight(kt)+diameter(kt)/2. ! height of upper blade tip above ground (m)
!
          kfound = .false.
          zheightl=0.0
          ! find vertical levels cut by turbine blades
          DO k=1,km
!          DO k=kts,ktf
            IF(.NOT. kfound) THEN
              zheightu = zheightl + dz(i,k) ! increment height

              IF(blade_l_point .GE. zheightl .AND. blade_l_point .LE. zheightu) THEN
                k_turbine_bot=k ! lower blade tip cuts this level
              ENDIF

              IF(blade_u_point .GE. zheightl .AND. blade_u_point .LE. zheightu) THEN
                k_turbine_top=k ! upper blade tip cuts this level
                kfound = .TRUE.
              ENDIF

              zheightl = zheightu
            ENDIF
          ENDDO
          IF ( kfound ) THEN
!
! ... PAJ: Changes introduced to compute only one set of turbine coefficients ...
!          First computes the wind speed at the hub height.
!
          kfound = .false.
          zheightl=0.
          ! find vertical levels (half levels) within the hub height
          DO k=1,km
!          DO k=kts,ktf
            IF(.NOT. kfound) THEN
              z2 = zheightl + 0.5*dz(i,k) 
!
              IF(hubheight(kt) .GE. z2 ) THEN
                kbot=k 
              ELSE
                ktop=k
                kfound = .TRUE.
              ENDIF
!
              if (.NOT. kfound) z1=z2
              zheightl = z2 + 0.5*dz(i,k)
            ENDIF
          ENDDO
!
          speed1=0.
          speed2=0.
          if (ktop.eq.1) then
           speedhub=sqrt(u(i,1)**2.+v(i,1)**2.)*hubheight(kt)/z1
          else
           speed1=sqrt(u(i,kbot)**2.+v(i,kbot)**2.)
           speed2=sqrt(u(i,ktop)**2.+v(i,ktop)**2.)
           speedhub=speed1+((speed2-speed1)/(z2-z1))*(hubheight(kt)-z1)
          endif

          spdh(kt)=spdh(kt)+speedhub

!
! ... calculate TKE, power and thrust coeffs
!
              CALL dragcof(tkecof,powcof,thrcof,               &
                           speedhub,cutin(kt),cutout(kt),   &
                           npower(kt),diameter(kt),stc(kt),stc2(kt),nkind(kt))
!
! ... PAJ: Computation of power generated by the wind turbine ...
!
          area=piconst/4.*diameter(kt)**2.          ! area swept by turbine blades
          power1=0.5*1.23*speedhub**3.*area*powcof
          power(kt)=power1+power(kt)
!          print*,'power:',power(kt),power1
          power2=0.
!
            DO k=k_turbine_bot,k_turbine_top ! loop over turbine blade levels
              z1=z_at_w(i,k)-blade_l_point-z_at_w(i,1)  ! distance between k level and lower blade tip
              z2=z_at_w(i,k+1)-blade_l_point-z_at_w(i,1) ! distance between k+1 level and lower blade tip
              IF(z1 .LT. 0.) z1=0.0 ! k level lower than lower blade tip
              IF(z2 .GT. diameter(kt)) z2=diameter(kt) ! k+1 level higher than turbine upper blade tip
              CALL turbine_area(z1,z2,diameter(kt),wfdensity,tarea)
!
              speed=sqrt(u(i,k)**2.+v(i,k)**2.)
              power2=power2+0.5*powcof*1.23*(speed**3.)*tarea/wfdensity
            ENDDO
!
! ... PAJ: Computes the tendencies of TKE and momentum ...
!
            DO k=k_turbine_bot,k_turbine_top ! loop over turbine blade levels
              z1=z_at_w(i,k)-blade_l_point-z_at_w(i,1)  ! distance between k lev and lower blade tip
              z2=z_at_w(i,k+1)-blade_l_point-z_at_w(i,1) !distance between k+1 lev and lower blade tip
              IF(z1 .LT. 0.) z1=0.0 ! k level lower than lower blade tip
              IF(z2 .GT. diameter(kt)) z2=diameter(kt) ! k+1 level higher than turbine upper blade tip
!
              CALL turbine_area(z1,z2,diameter(kt),wfdensity,tarea)
!
              speed=sqrt(u(i,k)**2.+v(i,k)**2.)
!
! ... PAJ: normalization introduced to conserve energy ...
!
              if (power1.eq.0.or.power2.eq.0) then
              ec=1.
              else
              ec=power1/power2
              endif
!
!              ! output TKE
!              qke(i,k,j) = qke(i,k,j)+speed**3.*tarea*tkecof*dt/dz(i,k,j)*ec
              ! output u tendency
              du(i,k) = du(i,k)-.5*u(i,k)*thrcof*speed*tarea/dz(i,k)*ec
              ! output v tendency
              dv(i,k) = dv(i,k)-.5*v(i,k)*thrcof*speed*tarea/dz(i,k)*ec
            ENDDO
          ENDIF
        ENDIF
        endif
!      ENDIF
!    ENDDO

      return
  END SUBROUTINE dragforce

! This subroutine calculates area of turbine between two vertical levels
! Input variables : 
!            z1 = distance between k level and lower blade tip
!            z2 = distance between k+1 level and lower blade tip
!            wfdensity = wind farm density in m^-2
!     tdiameter = turbine diameter
! Output variable :
!         tarea = area of turbine between two levels * wfdensity
  SUBROUTINE turbine_area(z1,z2,tdiameter,wfdensity,tarea)

  REAL, INTENT(IN) ::tdiameter,wfdensity
  REAL, INTENT(INOUT) ::z1,z2
  REAL, INTENT(OUT):: tarea
  REAL r,zc1,zc2
  REAL a1,a2

  r=tdiameter/2.              !r = turbine radius
  z1=r-z1                   !distance of kth level from turbine center 
  z2=r-z2                   !distance of k+1 th level from turbine center
  zc1=abs(z1)
  zc2=abs(z2)

  IF (zc1.ge.r) THEN
    a1=r*r*asin(1.0)
  ELSE
    a1=zc1*sqrt(r*r-zc1*zc1)+r*r*asin(zc1/r)
  ENDIF

  IF (zc2.ge.r) THEN
    a2=r*r*asin(1.0)
  ELSE
    a2=zc2*sqrt(r*r-zc2*zc2)+r*r*asin(zc2/r)
  ENDIF


  !turbine area between z1 and z2
  IF(z1 .GT. 0. .AND. z2 .GT. 0.) THEN
     tarea=a1-a2
  ELSE IF(z1 .LT. 0. .AND. z2 .LT. 0.) THEN
     tarea=a2-a1
  ELSE
     tarea=a2+a1
  ENDIF
  tarea=tarea*wfdensity      !turbine area * wind farm density 

      return
  END SUBROUTINE turbine_area


  SUBROUTINE dragcof(tkecof,powcof,thrcof,speed,cispeed,cospeed, &
                     tpower,tdiameter,stdthrcoef,stdthrcoef2,nkind)


  REAL, INTENT(IN):: speed, cispeed, cospeed, tpower,tdiameter,stdthrcoef,stdthrcoef2
  REAL, INTENT(OUT):: tkecof,powcof,thrcof
  REAL :: power,area,mspeed,hspeed
!
! ... PAJ ...
!
   INTEGER :: nkind,k,nu,nb
   LOGICAL :: vfound
   REAL :: fac1,fac2

  area=piconst/4.*tdiameter**2.          ! area swept by turbine blades

      vfound=.false.
      DO k=1,maxvals2
            IF(.NOT. vfound) THEN
              IF(turbws(nkind,k).GT.speed) THEN
                nu=k 
                nb=k-1
                vfound=.true.
              ENDIF
            ENDIF
      ENDDO
!
  IF (speed .LE. cispeed) THEN
     thrcof = stdthrcoef
  ELSE
    IF (speed .GE. cospeed) THEN
     thrcof = stdthrcoef2
     ELSE
     thrcof = turbtc(nkind,nb)+(turbtc(nkind,nu)-turbtc(nkind,nb))/(turbws(nkind,nu)-turbws(nkind,nb))*(speed-turbws(nkind,nb))
    ENDIF
  ENDIF
!
! ... power coeficient ...
!
  IF(speed .LE. cispeed .OR. speed .GE. cospeed) THEN
     power=0.
     powcof=0.
  ELSE
      fac1=1000./(0.5*1.23*turbws(nkind,nb)**3.*area)
      fac2=1000./(0.5*1.23*turbws(nkind,nu)**3.*area)
      power = turbpw(nkind,nb)+(turbpw(nkind,nu)-turbpw(nkind,nb))/(turbws(nkind,nu)-turbws(nkind,nb)) &
                               *(speed-turbws(nkind,nb))
      powcof = turbpw(nkind,nb)*fac1+(turbpw(nkind,nu)*fac2-turbpw(nkind,nb)*fac1)/(turbws(nkind,nu)-turbws(nkind,nb)) &
                                     *(speed-turbws(nkind,nb))
  ENDIF
!
  ! tke coefficient calculation 

  tkecof=thrcof-powcof
  IF(tkecof .LT. 0.) tkecof=0.
!
      return
  END SUBROUTINE dragcof
!===xb118
 subroutine init_module_wind_fitch_rsm(igrd1s,jgrd1s,lonlens,latlens,flat,flon,me)

  IMPLICIT NONE
! 
   CHARACTER*256 num,input
   real rlat(igrd1s,jgrd1s),rlon(igrd1s,jgrd1s)
   real flat(igrd1s*jgrd1s),flon(igrd1s*jgrd1s)
   INTEGER i,j,me,nval,k,jlat,ijlat,igrd1s,jgrd1s,lonlens,latlens

  real lat,lon,latb,late,lonb,lone
  integer fn,nn

    if (.not.windfarm_initialized) then
!IF
    power=0.
! open file:'windturbines.txt' to read number of wind turbine.
    open(170,file='windturbines.txt',form='formatted',status='old')
       fi_nt=0
 10    read(170,*,end=100)
       fi_nt=fi_nt+1
       goto 10

 100   continue
       rewind (170)
!ENDIF

! allocate variables
       allocate (nkind(fi_nt),wival(fi_nt),wjval(fi_nt),power(fi_nt),spdh(fi_nt),wim(fi_nt),wi(fi_nt))
       allocate (hubheight(fi_nt),stc(fi_nt),stc2(fi_nt),cutin(fi_nt),cutout(fi_nt),diameter(fi_nt),npower(fi_nt))
       wival=-9999
       wjval=-9999
       power =0.0
       spdh = 0.0
       wi = 0
       wim = 0

!IF
  do k=1,fi_nt
!read lat lon from 'windturbines.txt'
   read (170,*)lat,lon,nkind(k)
   if (me .eq. 0) write(*,*)'WINDFARM Turbine #',k,': Lat, lon = ',lat,lon,'; Type = ',nkind(k)

!            wival(k)=-9999
!            wjval(k)=-9999
!change from lat->i lon->j
      do j = 1,latlens
      jlat = (j-1)*IGRD1S
       do i=1,lonlens
       ijlat = i + jlat
       rlat(i,j) = flat(ijlat)*180./piconst
       rlon(i,j) = flon(ijlat)*180./piconst
       enddo
      enddo

     do j=1,latlens
      do i=1,lonlens

       if (j .eq. 1) then
       latb=rlat(i,j)-(rlat(i,j+1)-rlat(i,j))/2.
       late=rlat(i,j)+(rlat(i,j+1)-rlat(i,j))/2.
       else if (j .eq. latlens) then
       latb=rlat(i,j)-(rlat(i,j)-rlat(i,j-1))/2.
       late=rlat(i,j)+(rlat(i,j)-rlat(i,j-1))/2.
       else
       latb=rlat(i,j)-(rlat(i,j)-rlat(i,j-1))/2.
       late=rlat(i,j)+(rlat(i,j+1)-rlat(i,j))/2.
       endif

       if (i .eq. 1) then
       lonb=rlon(i,j)-(rlon(i+1,j)-rlon(i,j))/2.
       lone=rlon(i,j)+(rlon(i+1,j)-rlon(i,j))/2.
       else if (i .eq. lonlens) then
       lonb=rlon(i,j)-(rlon(i,j)-rlon(i-1,j))/2.
       lone=rlon(i,j)+(rlon(i,j)-rlon(i-1,j))/2.
       else
       lonb=rlon(i,j)-(rlon(i,j)-rlon(i-1,j))/2.
       lone=rlon(i,j)+(rlon(i+1,j)-rlon(i,j))/2.
       endif
      if((( latb .le. lat).and.(lat .le. late)).and.(( lonb .le. lon).and.(lon .le. lone))) then
       wival(k) = i
       wjval(k) = j
      endif

      enddo
     enddo

  enddo
  close(170)
#ifdef MP
      if( mype.eq.master ) then
#endif
       fn=175
      open (fn,file='power.txt',form='formatted',status='unknown')
      print*,'open output file: power.txt'
       fn=176
      open (fn,file='windspeed.txt',form='formatted',status='unknown')
      print*,'open output file: windspeed.txt'
#ifdef MP
       endif
#endif
         turbws=0.
         turbtc=0.
         turbpw=0.
         DO i=1,fi_nt
          write(num,*) nkind(i)
          num=adjustl(num)
          input="wind-turbine-"//trim(num)//".tbl"
          OPEN(file=TRIM(input),unit=49,FORM='FORMATTED',STATUS='OLD')
          READ (49,*,ERR=132)nval
          READ(49,*,ERR=132)hubheight(i),diameter(i),stc(i),npower(i)
            DO k=1,nval
              READ(49,*,ERR=132)turbws(nkind(i),k),turbtc(nkind(i),k),turbpw(nkind(i),k)
            ENDDO
          cutin(i)  = turbws(nkind(i),1)
          cutout(i) = turbws(nkind(i),nval)
          stc2(i) = turbtc(nkind(i),nval)
          close (49)
         ENDDO

 132   continue

      windfarm_initialized=.true.
    endif
!ENDIF

      return
     end subroutine init_module_wind_fitch_rsm

     subroutine writ_out_fitch_rsm(spd,pp)
      integer fn, i, ierr, nn
      logical alive,lprnt
      real spd(fi_nt),pp(fi_nt),spd0(fi_nt),pp0(fi_nt)
      lprnt=.false.
      spd0 = 0.0
       pp0 = 0.0

#ifdef MP
#if defined(CWB_MPMD) || defined(CWB_CFS)
!      call mpi_barrier(mpi_comm_rsm,ierr)
      call mpi_reduce(spd,spd0,fi_nt,MPIREAL,mpi_sum,              &
     &                master,mpi_comm_rsm,ierr)
!      call mpi_barrier(mpi_comm_rsm,ierr)
      call mpi_reduce(pp,pp0,fi_nt,MPIREAL,mpi_sum,              &
     &                master,mpi_comm_rsm,ierr)
!      call mpi_barrier(mpi_comm_rsm,ierr)
      call mpi_reduce(wi,wim,fi_nt,mpi_integer,mpi_sum,              &
     &                master,mpi_comm_rsm,ierr)
#else
!      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_reduce(spd,spd0,fi_nt,MPIREAL,mpi_sum,              &
     &                master,mpi_comm_world,ierr)
!      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_reduce(pp,pp0,fi_nt,MPIREAL,mpi_sum,              &
     &                master,mpi_comm_world,ierr)
!      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_reduce(wi,wim,fi_nt,mpi_integer,mpi_sum,              &
     &                master,mpi_comm_world,ierr)
#endif
#endif

#ifdef MP
      if( mype .eq. master ) then
#endif
       do i=1,fi_nt
         if (pp0(i) .lt. 0.001) pp0(i)=0.0
         if (spd0(i) .lt. 0.001) spd0(i)=0.0
       enddo
       if(lprnt)   print*,'mpi number:',wim
       if(lprnt)   print*,'hourly tot. power:',pp0
       if(lprnt)   print*,'hourly ave. wind speed:',spd0
       fn=175
      inquire(fn,exist=alive)
#ifdef __FUJITSU
       write (fn,*) pp0(:)
#else
       write (fn,"(*(g0,x))") pp0(:)
#endif
       fn=176
      inquire(fn,exist=alive)
#ifdef __FUJITSU
       write (fn,*) spd0(:)
#else
       write (fn,"(*(g0,x))") spd0(:)
#endif
#ifdef MP
      endif
#endif

       return
     end subroutine writ_out_fitch_rsm

     subroutine finalized_wind_fitch_rsm
      integer fn

       deallocate (nkind,wival,wjval,power,spdh,wim,wi)
       deallocate (hubheight,stc,stc2,cutin,cutout,diameter)

#ifdef MP
      if( mype .eq. master ) then
#endif
       fn=175
      close (fn)
       fn=176
      close (fn)
#ifdef MP
      endif
#endif

      return
 end subroutine finalized_wind_fitch_rsm
!===xb118
!
!!!  SUBROUTINE init_module_wind_fitch(id,config_flags,xlong,xlat,windfarm_initialized,&
!!!                                            ims,ime,jms,jme,its,ite,jts,jte,ids,ide,jds,jde)
!!!!
!!!  IMPLICIT NONE
!!!!
!!!   integer ims,ime,jms,jme,ids,ide,jds,jde
!!!   integer its,ite,jts,jte
!!!   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN) :: xlong,xlat
!!!   TYPE (grid_config_rec_type) :: config_flags
!!!   TYPE (PROJ_INFO) :: ts_proj
!!!   logical :: windfarm_initialized
!!!! 
!!!   CHARACTER*256 num,input,message_wind
!!!   real lat,lon,ts_rx,ts_ry
!!!   REAL :: known_lat, known_lon
!!!   INTEGER i,j,nval,k,id
!!!
!!!   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
!!!!
!!!      IF ( wrf_dm_on_monitor() ) THEN
!!!!
!!!! ... PAJ: Opens the file with the location of the wind turbines ...
!!!!
!!!        if ( config_flags%windfarm_ij .eq. 1 ) then
!!!          open(70,file='windturbines-ij.txt',form='formatted',status='old')
!!!        else
!!!          open(70,file='windturbines.txt',form='formatted',status='old')
!!!        end if
!!!!
!!!! ... PAJ: Counts the turbines ...
!!!!
!!!       nt=0
!!! 10    read(70,*,end=100) 
!!!       nt=nt+1
!!!       goto 10
!!!!
!!! 100   continue
!!!       rewind (70)
!!!     END IF
!!!!
!!!     CALL wrf_dm_bcast_integer(nt,1)
!!!!
!!!! ... PAJ: Initializes the configuration of the wind farm(s) ...
!!!!
!!!     if (.not. windfarm_initialized) then
!!!       allocate (nkind(nt),ival(nt,max_domains),jval(nt,max_domains))
!!!       allocate (hubheight(nt),stc(nt),stc2(nt),cutin(nt),cutout(nt),diameter(nt),npower(nt))
!!!       ival=-9999
!!!       jval=-9999
!!!       windfarm_initialized=.true.
!!!     endif
!!!!
!!!     IF ( wrf_dm_on_monitor() ) THEN
!!!     do k=1,nt
!!!       if ( config_flags%windfarm_ij .eq. 1 ) then
!!!         read(70,*) ival(k,id), jval(k,id), nkind(k)
!!!         write(message_wind,*)'WINDFARM Turbine #',k,': I, J = ',ival(k,id), jval(k,id),'; Type = ',nkind(k)
!!!         CALL wrf_message(message_wind)
!!!
!!!       else
!!!
!!!         read(70,*)lat,lon,nkind(k)
!!!         write(message_wind,*)'WINDFARM Turbine #',k,': Lat, lon = ',lat,lon,'; Type = ',nkind(k)
!!!         CALL wrf_message(message_wind)
!!!
!!!         CALL map_init(ts_proj)
!!!
!!!         known_lat = xlat(its,jts)
!!!         known_lon = xlong(its,jts)
!!!
!!!      ! Mercator
!!!      IF (config_flags%map_proj == PROJ_MERC) THEN
!!!         CALL map_set(PROJ_MERC, ts_proj,               &
!!!                      truelat1 = config_flags%truelat1, &
!!!                      lat1     = known_lat,             &
!!!                      lon1     = known_lon,             &
!!!                      knowni   = REAL(its),             &
!!!                      knownj   = REAL(jts),             &
!!!                      dx       = config_flags%dx)
!!!
!!!      ! Lambert conformal
!!!      ELSE IF (config_flags%map_proj == PROJ_LC) THEN
!!!         CALL map_set(PROJ_LC, ts_proj,                  &
!!!                      truelat1 = config_flags%truelat1,  &
!!!                      truelat2 = config_flags%truelat2,  &
!!!                      stdlon   = config_flags%stand_lon, &
!!!                      lat1     = known_lat,              &
!!!                      lon1     = known_lon,              &
!!!                      knowni   = REAL(its),              &
!!!                      knownj   = REAL(jts),              &
!!!                      dx       = config_flags%dx)
!!!!      ! Polar stereographic
!!!      ELSE IF (config_flags%map_proj == PROJ_PS) THEN
!!!         CALL map_set(PROJ_PS, ts_proj,                  &
!!!                      truelat1 = config_flags%truelat1,  &
!!!                      stdlon   = config_flags%stand_lon, &
!!!                      lat1     = known_lat,              &
!!!                      lon1     = known_lon,              &
!!!                      knowni   = REAL(its),              &
!!!                      knownj   = REAL(jts),              &
!!!                      dx       = config_flags%dx)
!!!!#if (EM_CORE == 1)
!!!!      ! Cassini (global ARW)
!!!!      ELSE IF (config_flags%map_proj == PROJ_CASSINI) THEN
!!!!         CALL map_set(PROJ_CASSINI, ts_proj,                            &
!!!!                      latinc   = grid%dy*360.0/(2.0*EARTH_RADIUS_M*PI), &
!!!!                      loninc   = grid%dx*360.0/(2.0*EARTH_RADIUS_M*PI), &
!!!!                      lat1     = known_lat,                             &
!!!!                      lon1     = known_lon,                             &
!!!!                      lat0     = config_flags%pole_lat,                 &
!!!!                      lon0     = config_flags%pole_lon,                 &
!!!!                      knowni   = 1.,                                    &
!!!!                      knownj   = 1.,                                    &
!!!!                      stdlon   = config_flags%stand_lon)
!!!!#endif
!!!!
!!!!      ! Rotated latitude-longitude
!!!!      ELSE IF (config_flags%map_proj == PROJ_ROTLL) THEN
!!!!         CALL map_set(PROJ_ROTLL, ts_proj,                      &
!!!!! I have no idea how this should work for NMM nested domains
!!!!                      ixdim    = grid%e_we-1,                   &
!!!!                      jydim    = grid%e_sn-1,                   &
!!!!                      phi      = real(grid%e_sn-2)*grid%dy/2.0, &
!!!!                      lambda   = real(grid%e_we-2)*grid%dx,     &
!!!!                      lat1     = config_flags%cen_lat,          &
!!!!                      lon1     = config_flags%cen_lon,          &
!!!!                      latinc   = grid%dy,                       &
!!!!                      loninc   = grid%dx,                       &
!!!!                      stagger  = HH)
!!!!
!!!      END IF
!!!!
!!!         CALL latlon_to_ij(ts_proj, lat, lon, ts_rx, ts_ry)
!!!!
!!!          ival(k,id)=nint(ts_rx)
!!!          jval(k,id)=nint(ts_ry)
!!!          if (ival(k,id).lt.ids.and.ival(k,id).gt.ide) then
!!!            ival(k,id)=-9999
!!!            jval(k,id)=-9999
!!!          endif
!!!!
!!!!         write(73,*) k,id,ival(k,id),jval(k,id)
!!!          write(message_wind,*)'WINDFARM Turbine #',k,': Lat, lon = ',lat,lon, &
!!!                               ', (i,j) = (',ival(k,id),',',jval(k,id),'); Type = ',nkind(k)
!!!          CALL wrf_debug(0,message_wind)
!!!!
!!!     end if
!!!!
!!!     enddo
!!!      close(70)
!!!!
!!!! ... PAJ: Read the tables with the turbine's characteristics ...
!!!!
!!!         turbws=0.
!!!         turbtc=0.
!!!         turbpw=0.
!!!         DO i=1,nt
!!!          write(num,*) nkind(i)
!!!          num=adjustl(num)
!!!          input="wind-turbine-"//trim(num)//".tbl"
!!!          OPEN(file=TRIM(input),unit=19,FORM='FORMATTED',STATUS='OLD')
!!!          READ (19,*,ERR=132)nval
!!!          READ(19,*,ERR=132)hubheight(i),diameter(i),stc(i),npower(i)
!!!            DO k=1,nval
!!!              READ(19,*,ERR=132)turbws(nkind(i),k),turbtc(nkind(i),k),turbpw(nkind(i),k)
!!!            ENDDO
!!!          cutin(i)  = turbws(nkind(i),1)
!!!          cutout(i) = turbws(nkind(i),nval)
!!!          stc2(i) = turbtc(nkind(i),nval)
!!!          close (19)
!!!         ENDDO
!!!
!!! 132   continue
!!!!
!!!! ... ...
!!!!
!!!      endif
!!!
!!!        CALL wrf_dm_bcast_integer(ival,nt*max_domains)
!!!        CALL wrf_dm_bcast_integer(jval,nt*max_domains)
!!!        CALL wrf_dm_bcast_real(hubheight,nt)
!!!        CALL wrf_dm_bcast_real(diameter,nt)
!!!        CALL wrf_dm_bcast_real(stc,nt)
!!!        CALL wrf_dm_bcast_real(npower,nt)
!!!        CALL wrf_dm_bcast_real(cutin,nt)
!!!        CALL wrf_dm_bcast_real(cutout,nt)
!!!        CALL wrf_dm_bcast_integer(nkind,nt) 
!!!        CALL wrf_dm_bcast_real(turbws,maxvals*maxvals2) 
!!!        CALL wrf_dm_bcast_real(turbtc,maxvals*maxvals2) 
!!!        CALL wrf_dm_bcast_real(turbpw,maxvals*maxvals2) 
!!!
!!!  END SUBROUTINE init_module_wind_fitch
  
END MODULE module_wind_fitch
