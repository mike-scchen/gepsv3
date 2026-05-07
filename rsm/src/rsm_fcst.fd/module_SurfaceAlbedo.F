module module_SurfaceAlbedo
implicit none
public

integer r8
parameter (r8 = 8)
real(r8), allocatable ::                                                 &
    sky_view(:,:),terrain_config(:,:),sinsl_cosas(:,:),sinsl_sinas(:,:), &
    msk_of_3ddata(:,:)

save
contains
!
  subroutine allocate_3drtdata(igrd1s,jgrd1s)
  integer igrd1s,jgrd1s,ierr
  allocate(sky_view(igrd1s,jgrd1s),terrain_config(igrd1s,jgrd1s),        &
           sinsl_cosas(igrd1s,jgrd1s),sinsl_sinas(igrd1s,jgrd1s),        &
           msk_of_3ddata(igrd1s,jgrd1s),                                 &
           stat=ierr)

  if (ierr/=0) stop "SurfaceAlbedoMod_allocate_3drtdata: fail"
  return
  end subroutine
!
  subroutine deallocate_3drtdata
  integer ierr
  deallocate(sky_view,terrain_config,sinsl_cosas,sinsl_sinas,            &
             msk_of_3ddata,stat=ierr)

  if (ierr/=0) stop "SurfaceAlbedoMod_deallocate_3drtdata: fail"
  return
  end subroutine
!
  subroutine read_3drtdata(igrd1,jgrd1,igrd1s,jgrd1s)
#include <define.h>
  integer igrd1,jgrd1,igrd1s,jgrd1s,ierr,recsize,iunit
  real(r8), allocatable :: tmp_read(:,:)
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
  msk_of_3ddata=1.0
  if( mype.eq.master ) then
#endif

  allocate(tmp_read(igrd1,jgrd1))

  iunit=99
  inquire(IOLENGTH=recsize) tmp_read
  open(iunit,file='rsm_3drt.bin',                                        &
       form='unformatted',access='direct',convert='big_endian',          &
       status='old',recl=recsize,iostat=ierr)
  if (ierr/=0) stop "SurfaceAlbedoMod_read_3drtdata: open file fail"

! 1. topo
  read(iunit,rec=1,iostat=ierr) tmp_read
! print*,'--3d--terr read in=',tmp_read(igrd1/2,jgrd1/2)
  if (ierr/=0) stop "SurfaceAlbedoMod_read_3drtdata: read fail1"
#ifdef MP
  endif
#endif      
!#ifdef MP
!  call rmpgf2p(tmp_read,igrd1,jgrd1,terrain_config,igrd1s,jgrd1s,1)
!#else
!  terrain=tmp_read
!  where(terrain_config .le. -9999.) msk_of_3ddata=0.0
!#endif

! 2. sinsl_cosas
#ifdef MP
  if( mype.eq.master ) then
#endif
  read(iunit,rec=2,iostat=ierr) tmp_read
! print*,'--3d--sinsl_cosas read in=',tmp_read(igrd1/2,jgrd1/2)
  if (ierr/=0) stop "SurfaceAlbedoMod_read_3drtdata: read fail2"
#ifdef MP
  endif
#endif      
#ifdef MP
  call rmpgf2p(tmp_read,igrd1,jgrd1,sinsl_cosas,igrd1s,jgrd1s,1)
#else
  sinsl_cosas=tmp_read
  where(sinsl_cosas .le. -9999.) msk_of_3ddata=0.0
#endif

! 3. sinsl_sinas
#ifdef MP
  if( mype.eq.master ) then
#endif
  read(iunit,rec=3,iostat=ierr) tmp_read
! print*,'--3d--sinsl_sinas read in=',tmp_read(igrd1/2,jgrd1/2)
  if (ierr/=0) stop "SurfaceAlbedoMod_read_3drtdata: read fail3"
#ifdef MP
  endif
#endif      
#ifdef MP
  call rmpgf2p(tmp_read,igrd1,jgrd1,sinsl_sinas,igrd1s,jgrd1s,1)
#else
  sinsl_sinas=tmp_read
  where(sinsl_sinas .le. -9999.) msk_of_3ddata=0.0
#endif
  
! 4. svf
#ifdef MP
  if( mype.eq.master ) then
#endif
  read(iunit,rec=4,iostat=ierr) tmp_read
! print*,'--3d--svf read in=',tmp_read(igrd1/2,jgrd1/2)
  if (ierr/=0) stop "SurfaceAlbedoMod_read_3drtdata: read fail4"
#ifdef MP
  endif
#endif      
#ifdef MP
  call rmpgf2p(tmp_read,igrd1,jgrd1,sky_view,igrd1s,jgrd1s,1)
#else
  sky_view=tmp_read
  where(sky_view .le. -9999.) msk_of_3ddata=0.0
#endif

! 5. tcf
#ifdef MP
  if( mype.eq.master ) then
#endif
 read(iunit,rec=5,iostat=ierr) tmp_read
 if (ierr/=0) stop "SurfaceAlbedoMod_read_3drtdata: read fail5"
#ifdef MP
  endif
#endif      
#ifdef MP
  call rmpgf2p(tmp_read,igrd1,jgrd1,terrain_config,igrd1s,jgrd1s,1)
#else
  terrain_config=tmp_read
  where(terrain_config .le. -9999.) msk_of_3ddata=0.0
#endif
  
#ifdef MP
  if( mype.eq.master ) then
#endif
  close(iunit)
#ifdef MP
  endif
#endif      

  return
  end subroutine
!
!----------------------------------------------------------------------
!  Added by Wei-Liang Lee for 3D topography effect
!----------------------------------------------------------------------
  subroutine Albedo_3D_Adjustment(im,nextsw_cday,lon,lat,albd,albi,      &
                                  sdec,coszen,slmsk,lan)
! !DESCRIPTION:
! Adjust surface albedo by accounting for 3-D topography effect on surface solar radiation
!

! !USES:
!   use clmtype
!   use shr_orb_mod
!   use clm_varctl  , only: iulog

!
! !ARGUMENTS:
    implicit none

    integer, intent(in) :: im                 ! length of the array
    real(r8), intent(in) :: nextsw_cday       ! calendar day at Greenwich (1.00, ..., days/year)
! in-out arguments
!
    real(r8), intent(in) :: lon(im)           ! gridcell longitude (radians)
    real(r8), intent(in) :: lat(im)           ! gridcell latitude (radians)
    integer, parameter :: numrad = 1          ! number of radiation 
    real(r8), intent(inout) :: albd(im,numrad)! surface albedo (direct)
    real(r8), intent(inout) :: albi(im,numrad)! surface albedo (diffuse)
    real(r8), intent(in) :: sdec              ! sin(solar declination angle (radians))  decl=asin(sdec)
    real(r8), intent(in) :: coszen(im)        ! cosine of solar zenith angle
    real(r8), intent(in) :: slmsk(im)         ! sea land mask - sea(0),land(1),ice(2) on fcst model grid
    integer, intent(in) :: lan                ! latitude index
!
! !CALLED FROM:
! subroutine SurfaceAlbedo

! !REVISION HISTORY:
! Author: Wei-Liang Lee, 01/16/2015
! modified for RSM/MSM: Ying-Ju Chen, 07/28/2020
!

!
! local pointers to implicit in arguments
!
!   integer , pointer :: pgridcell(:) ! gridcell of corresponding pft
!   integer , pointer :: pcolumn(:)   ! column of corresponding pft
!    
! get the following data from module >
!   real(r8), pointer :: stdev_elev(:)       ! standard deviation of elevation within a gridcell
!   real(r8), pointer :: sky_view(:)         ! mean of (sky view factor / cos(slope))
!   real(r8), pointer :: terrain_config(:)   ! mean of (terrain configuration factor / cos(slope))
!   real(r8), pointer :: sinsl_cosas(:)      ! sin(slope)*cos(aspect) / cos(slope)
!   real(r8), pointer :: sinsl_sinas(:)      ! sin(slope)*sin(aspect) / cos(slope)
! get the above     data from module <
!
!
! !OTHER LOCAL VARIABLES:
!
    real(r8), parameter :: mpe = 1.e-06_r8 ! prevents overflow for division by zero
    real(r8), parameter :: pi = 3.14159265358979323846_r8 ! pi
    integer  :: fp,fc,ii                   ! indices
    integer  :: ib                         ! band index
    integer  :: izen
    real(r8) :: cosz    ! cosine solar zenith angle for next time step
    real(r8) :: sinz    ! sine of solar zenith angle
    real(r8) :: azi_angle    ! solar azimuth angle
    real(r8) :: next_tod    ! time of day for nextsw_cday in second
    real(r8) :: solar_inc ! (solar incident angle) / cos(slope) / cosz
    real(r8) :: f_dir     ! adjustment factor for direct flux
    real(r8) :: f_rdir, f_rdir_temp    ! adjustment factor for reflected-direct flux
    real(r8) :: f_dif     ! adjustment factor for diffuse flux
    real(r8) :: f_rdif, f_rdif_temp    ! adjustment factor for reflected-diffuse flux
    real(r8) :: fd_prime     ! temp adjustment factor for direct flux
    real(r8) :: fi_prime     ! temp adjustment factor for diffuse flux
    real(r8) :: albd_adjust     ! adjusted albedo for direct flux
    real(r8) :: albi_adjust     ! adjusted albedo for diffuse flux
    real(r8) :: ftemp1, ftemp2, dzen1, dzen2
    real(r8) :: local_timeofday     ! local time of day (second)
    real(r8) :: coeff_dir(3,0:7)  ! regression coefficients for f_dir
    real(r8) :: coeff_rdir(3,0:7) ! regression coefficients for f_rdir
    real(r8) :: coeff_dif(4,0:7)  ! regression coefficients for f_dif
    real(r8) :: coeff_rdif(3,0:7) ! regression coefficients for f_rdif

    data coeff_dir(:,1) /2.045E+1, 6.792E-1, -2.103E+1/
    data coeff_dir(:,2) /1.993E+0, 9.284E-1, -2.911E+0/
    data coeff_dir(:,3) /5.900E-2, 9.863E-1, -1.045E+0/
    data coeff_dir(:,4) /5.270E-3, 9.942E-1, -9.995E-1/
    data coeff_dir(:,5) /2.977E-3, 9.959E-1, -9.990E-1/
    data coeff_dir(:,6) /2.977E-3, 9.959E-1, -9.990E-1/
    data coeff_dir(:,7) /8.347E-3,        0, -8.393E-3/

    data coeff_rdir(:,1) / 2.351E-1, 1.590E-1, -2.332E-1/
    data coeff_rdir(:,2) / 1.368E-1, 1.642E-1, -1.358E-1/
    data coeff_rdir(:,3) / 1.254E-1, 1.653E-1, -1.247E-1/
    data coeff_rdir(:,4) / 1.274E-1, 1.635E-1, -1.267E-1/
    data coeff_rdir(:,5) / 1.314E-1, 1.623E-1, -1.307E-1/
    data coeff_rdir(:,6) / 1.359E-1, 1.620E-1, -1.352E-1/
    data coeff_rdir(:,7) /-4.463E-6, 1.556E-1,  1.287E-3/

    data coeff_dif(:,1) /3.146E-7, 4.385E+0, 6.723E-3, -4.382E+0/
    data coeff_dif(:,2) /6.001E-7, 4.068E+0, 2.456E-2, -4.085E+0/
    data coeff_dif(:,3) /7.436E-7, 3.911E+0, 5.606E-2, -3.960E+0/
    data coeff_dif(:,4) /7.806E-7, 3.763E+0, 1.049E-1, -3.863E+0/
    data coeff_dif(:,5) /7.581E-7, 3.559E+0, 1.734E-1, -3.727E+0/
    data coeff_dif(:,6) /7.015E-7, 3.298E+0, 2.543E-1, -3.547E+0/
    data coeff_dif(:,7) /6.359E-7, 2.984E+0,        0, -2.984E+0/

    data coeff_rdif(:,1) / 1.493E-1, 1.621E-1, -1.483E-1/
    data coeff_rdif(:,2) / 1.462E-1, 1.654E-1, -1.454E-1/
    data coeff_rdif(:,3) / 1.454E-1, 1.673E-1, -1.446E-1/
    data coeff_rdif(:,4) / 1.465E-1, 1.683E-1, -1.457E-1/
    data coeff_rdif(:,5) / 1.443E-1, 1.682E-1, -1.435E-1/
    data coeff_rdif(:,6) / 1.446E-1, 1.686E-1, -1.439E-1/

    data coeff_rdif(:,7) /-3.427E-6, 1.576E-1,  1.199E-3/

    ! Assign local pointers to derived subtypes components (gridcell-level)
!yj, include rscomgrd
!   lat              => grc%lat  !flat
!   lon              => grc%lon  !flon
!    
! read in data >    
!
!   stdev_elev       => grc%stdev_elev  
!   sky_view         => grc%sky_view  !???
!   terrain_config   => grc%terrain_config !???
!   sinsl_cosas      => grc%sinsl_cosas !???
!   sinsl_sinas      => grc%sinsl_sinas !???
!    
! read in data <    
!

    ! Assign local pointers to derived subtypes components (column-level)

!   decl           => cps%decl  !dlt(in astronomy.F,or sdec(sin),cdec(cos) in rscomgrd.h)
!   coszen         => cps%coszen !coszer (in radiation_module.F,rscomgrd.h)

    ! Assign local pointers to derived subtypes components (pft-level)

!   pgridcell => pft%gridcell ! ii?
!   pcolumn   => pft%column   ! jj?
!   albd      => pps%albd !alvbf(vis beam)+alnbf(nir beam) in grrad.F
!   albi      => pps%albi !alvdf(vis diff)+alndf(nir diff) in grrad.F


    coeff_dir(:,0) = coeff_dir(:,1)
    coeff_rdir(:,0) = coeff_rdir(:,1)
    coeff_dif(:,0) = coeff_dif(:,1)
    coeff_rdif(:,0) = coeff_rdif(:,1)

    next_tod = 86400._r8 * (nextsw_cday - int(nextsw_cday))

!   do fp = 1,num_nourbanp
!      g = pgridcell(p)
!      c = pcolumn(p)
    do ii=1,im
       if (( slmsk(ii).eq. 1.0 ).and.(msk_of_3ddata(ii,lan).eq.1.0)) then
       cosz = coszen(ii)

! If the sun is below the horizon, latitude larger than 60 degree, or the grid
! is flat, don't do anything.
!      if (cosz > 0._r8 .and. abs(lat(ii)) < 1.047_r8 .and. stdev_elev(ii) > 0._r8) then
       if (cosz > 0._r8 .and. abs(lat(ii)) < 1.047_r8 ) then
!         !print
!         print*,'-in 3d code- terrain_config=',terrain_config(ii,lan)
!         print*,'-in 3d code- sinsl_cosas=',sinsl_cosas(ii,lan)
!         print*,'-in 3d code- sinsl_sinas=',sinsl_sinas(ii,lan)
!         print*,'-in 3d code- sky_view=',sky_view(ii,lan)
          local_timeofday = next_tod + lon(ii) / pi * 180._r8 * 240._r8

          if (local_timeofday >= 86400._r8) then
             local_timeofday = local_timeofday - 86400._r8
          endif
    
          if (cosz == 1._r8) then
             azi_angle = 0._r8
             sinz = 0._r8
             solar_inc = 1
          else
             sinz = sqrt(1._r8-cosz*cosz)
             azi_angle = (sin(lat(ii))*cosz- sdec ) / (cos(lat(ii))*sinz)
             azi_angle = max(-1._r8,min(1._r8,azi_angle))
             azi_angle = acos(-azi_angle)
             if (local_timeofday >=43200._r8) then
                azi_angle = 2._r8*pi - azi_angle
             endif
             azi_angle = pi / 2._r8 - azi_angle
             solar_inc = 1+(sinz/cosz)*(cos(azi_angle)*sinsl_cosas(ii,lan)+sin(azi_angle)*sinsl_sinas(ii,lan))
          endif

          izen = int((cosz + 0.05_r8) / 0.15_r8)
          dzen1 = (cosz - (izen * 0.15_r8 - 0.05_r8)) / 0.15_r8
          dzen2 = 1._r8 - dzen1

          ftemp1 = coeff_dir(1,izen) * sky_view(ii,lan) + &
                   coeff_dir(2,izen) * solar_inc + coeff_dir(3,izen)
          ftemp2 = coeff_dir(1,izen+1) * sky_view(ii,lan) + &
                   coeff_dir(2,izen+1) * solar_inc + coeff_dir(3,izen+1)
          f_dir = ftemp2 * dzen1 + ftemp1 * dzen2
          f_dir = max(-1._r8,f_dir)

          ftemp1 = coeff_rdir(1,izen) * sky_view(ii,lan) + &
                   coeff_rdir(2,izen) * terrain_config(ii,lan) + coeff_rdir(3,izen)
          ftemp2 = coeff_rdir(1,izen+1) * sky_view(ii,lan) + &
                   coeff_rdir(2,izen+1) * terrain_config(ii,lan) + coeff_rdir(3,izen+1)
          f_rdir_temp = ftemp2 * dzen1 + ftemp1 * dzen2

!         ftemp1 = coeff_dif(1,izen) * stdev_elev(ii) + &
          ftemp1 =                                     &
                   coeff_dif(2,izen) * sky_view(ii,lan) + &
                   coeff_dif(3,izen) * solar_inc + coeff_dif(4,izen)
!         ftemp2 = coeff_dif(1,izen+1) * stdev_elev(ii) + &
          ftemp2 =                                       &
                   coeff_dif(2,izen+1) * sky_view(ii,lan) + &
                   coeff_dif(3,izen+1) * solar_inc + coeff_dif(4,izen+1)
          f_dif = ftemp2 * dzen1 + ftemp1 * dzen2
          f_dif = max(-1._r8,f_dif) 

          ftemp1 = coeff_rdif(1,izen) * sky_view(ii,lan) + &
                   coeff_rdif(2,izen) * terrain_config(ii,lan) + coeff_rdif(3,izen)
          ftemp2 = coeff_rdif(1,izen+1) * sky_view(ii,lan) + &
                   coeff_rdif(2,izen+1) * terrain_config(ii,lan) + coeff_rdif(3,izen+1)
          f_rdif_temp = ftemp2 * dzen1 + ftemp1 * dzen2

          do ib = 1, numrad
             f_rdir = f_rdir_temp * (albd(ii,ib) / 0.1_r8)
             f_rdif = f_rdif_temp * (albi(ii,ib) / 0.1_r8)

             fd_prime = 1._r8 + f_dir + f_rdir
             fi_prime = 1._r8 + f_dif + f_rdif

             albd_adjust = fd_prime * albd(ii,ib) - (fd_prime-1._r8)
             albi_adjust = fi_prime * albi(ii,ib) - (fi_prime-1._r8)

             if (albd_adjust <= 0) then 
!               print*,'lt0=',f_rdir_temp
!               print*,'lt1=',f_dir,f_rdir
!               print*,'lt2=',fd_prime,albd(ii,ib),(fd_prime-1._r8)
                albd_adjust = 0
                fd_prime = 1._r8 / (1._r8 - albd(ii,ib))
             endif

             if (albi_adjust <= 0) then
                albi_adjust = 0
                fi_prime = 1._r8 / (1._r8 - albi(ii,ib))
             endif

             !print-debug>>
!            print*,'new,old albd(ii)=',ii,albd(ii,ib),albd_adjust
             !print-debug<<
             albd(ii,ib) = albd_adjust
             albi(ii,ib) = albi_adjust
          enddo

       endif
       endif !end slmsk
    enddo


  end subroutine Albedo_3D_Adjustment
!--- End Lee ----------------------------------------------------------



end module module_SurfaceAlbedo
