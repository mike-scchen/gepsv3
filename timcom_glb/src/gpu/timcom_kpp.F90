!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module timcom_kpp_gpu
  use timcom_const
  use timcom_comm
  use timcom_comm_gpu
  use hyperlink, only: nx, ny, nz, nzf, kb, in, &
                       dt, iw, dz, odz, dzw, odzw, &
                       z_grid, z_face, x_grid, y_grid,&
                       myid, myid_x, myid_y
  implicit none

  real(r8), parameter :: &
    convect_diff =  20000.d0, &       ! diffusivity to mimic convection
    convect_visc =  20000.d0          ! viscosity   to mimic convection

contains

#ifdef cpl
 #ifndef clm_r
 subroutine kppglo_gpu(u2, v2, t2, s2, rho, trans, &
                      stf, shf_qsw, smft, vbk, hbk, &
                      add, vdc, vvc, ev, hv, kpp_src, kpp_hblt, async_id)
 #else
 subroutine kppglo_gpu(itf, u2, v2, t2, s2, rho, trans, &
                      stf, shf_qsw, smft, vbk, hbk, &
                      add, vdc, vvc, ev, hv, kpp_src, kpp_hblt, async_id)
 #endif

#else
subroutine kppglo_gpu(itf, u2, v2, t2, s2, rho, trans, &
                      stf, shf_qsw, smft, vbk, hbk, &
                      add, vdc, vvc, ev, hv, kpp_src, kpp_hblt, async_id)
#endif
  implicit none
  integer :: itf
  integer, intent(in) :: async_id
  real(r8), dimension(-1:nx+2,-1:ny+2,nz) :: &
    u2, v2, t2, s2

  real(r8), dimension(nx,nx,nz) :: &
    rho, trans

  real(r8), optional :: &
    stf(nx,ny,2), &
    shf_qsw(nx,ny), &
    smft(nx,ny,2)

  real(r8) :: &
    vbk(nz), &
    hbk(nz), &
    add(nx,ny,nz)

  real(r8) :: &
    vdc(nx,ny,0:nzf,2), &
    vvc(nx,ny,nz), &
    ev(nx,ny,nz), &
    hv(nx,ny,nz,2), &
    kpp_src(nx,ny,nz,2), &
    kpp_hblt(nx,ny)

  integer :: i,j,k,l
  real(r8) :: emax, dtmp
  real(r8) :: trcr(nx,ny,nz,2)

  !$acc enter data async(async_id) create(trcr)

  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        trcr(i,j,k,1) = t2(i,j,k)
        trcr(i,j,k,2) = s2(i,j,k)*ppt_to_salt
      end do
    end do
  end do

  #ifdef chk_kpp
    if(myid_x .eq. 7 .and. myid_y .eq. 11) then
      write(*,*) "[KPP CHECK]: 01", itf
      write(*,*) "hv   = ",  hv(nx/2,ny,1,1)
      write(*,*) "vdct = ", vdc(nx/2,ny,1,1)
      write(*,*) "vdcs = ", vdc(nx/2,ny,1,2)
      write(*,*) "ev   = ",  ev(nx/2,ny,1)
      write(*,*) "vvc  = ", vvc(nx/2,ny,1)
      write(*,*) "rho  = ", rho(nx/2,ny,1)
      write(*,*) "hblt = ", kpp_hblt(nx/2,ny)
      write(*,*) "smft = ", smft(nx/2,ny,1:2)
      write(*,*) "t2   = ", trcr(nx/2,ny,1,1)
      write(*,*) "s2   = ", trcr(nx/2,ny,1,2)
      write(*,*) "u2   = ", u2(nx/2,ny,1)
      write(*,*) "v2   = ", v2(nx/2,ny,1)
      write(*,*) "emax = ", 0.5d0/(dt*odzw(2:4)**2)
      write(*,*) "odzw = ", odzw(2:4)*iw(nx/2,ny,2:4)
      write(*,*) achar(10)
    end if
  #endif

  call vmix_coeffs_kpp_gpu(trans, kpp_src, vdc, vvc, trcr, u2, v2, rho,&
                           stf, shf_qsw, kpp_hblt, convect_diff, convect_visc,&
                           smft=smft, async_id=async_id)


  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz-1
    do j = 1, ny
      do i = 1, nx
        l = k+1
        emax = 0.5d0/(dt*odzw(l)**2)
        dtmp = odzw(l)*iw(i,j,l)
        ev(i,j,k) = dtmp*max(min(emax,vvc(i,j,k)),vbk(k)*5.d0)
        hv(i,j,k,1) = dtmp*max(min(emax,vdc(i,j,k,1)),hbk(k)*5.d0)
        hv(i,j,k,2) = dtmp*max(min(emax,vdc(i,j,k,2)),hbk(k)*5.d0)
      end do
    end do
  end do

  #ifdef chk_kpp
    if(myid_x .eq. 7 .and. myid_y .eq. 11) then
      write(*,*) "[KPP CHECK]: 02", itf
      write(*,*) "hv   = ",  hv(nx/2,ny,1,1)
      write(*,*) "vdct = ", vdc(nx/2,ny,1,1)
      write(*,*) "vdcs = ", vdc(nx/2,ny,1,2)
      write(*,*) "ev   = ",  ev(nx/2,ny,1)
      write(*,*) "vvc  = ", vvc(nx/2,ny,1)
      write(*,*) "rho  = ", rho(nx/2,ny,1)
      write(*,*) "hblt = ", kpp_hblt(nx/2,ny)
      write(*,*) "smft = ", smft(nx/2,ny,1:2)
      write(*,*) "t2   = ", trcr(nx/2,ny,1,1)
      write(*,*) "s2   = ", trcr(nx/2,ny,1,2)
      write(*,*) "u2   = ", u2(nx/2,ny,1)
      write(*,*) "v2   = ", v2(nx/2,ny,1)
      write(*,*) "emax = ", 0.5d0/(dt*odzw(2:4)**2)
      write(*,*) "odzw = ", odzw(2:4)*iw(nx/2,ny,2:4)
      write(*,*) achar(10)
    end if
  #endif

  !$acc exit data async(async_id) delete(trcr)

end subroutine kppglo_gpu

subroutine vmix_coeffs_kpp_gpu(trans, kpp_src, vdc, vvc, trcr, uuu, vvv, rhomix,&
                               stf, shf_qsw, kpp_hblt, convect_diff, convect_visc,&
                               smf, smft, async_id)
  ! !DESCRIPTION:
  !  This is the driver routine which calculates the vertical
  !  mixing coefficients for the KPP mixing scheme as outlined in
  !  Large, McWilliams and Doney, Reviews of Geophysics, 32, 363
  !  (November 1994).  The non-local mixing is also computed here, but
  !  is treated as a source term in baroclinic.
  use timcom_kpp, only: hmxl, ldbl_diff, zgrid, BVSQcon, &
                        tidal_diff, tidal_coef, bckgrnd_vvc, &
                        bckgrnd_vdc, zgrid, rrho0, dsfmax, hwide, cg, &
                        fcort, bolus_sp, fstokes
  implicit none
  ! !INPUT PARAMETERS:
  !  integer(2), dimension(nx,ny,nz), intent(in) :: IN

  integer, intent(in) :: async_id

  real(r8), dimension(nx,ny,nz,2), intent(inout) :: &
    trcr                ! tracers at current time

  real(r8), dimension(-1:nx+2,-1:ny+2,nz), intent(in) :: &
    uuu, vvv            ! velocities at current time

  real(r8), dimension(nx,ny,nz), intent(in) :: &
    trans, &            ! transmission (1 to 0) from surface (calculated in kppglo)
    rhomix              ! density at mix time

  real(r8), dimension(nx,ny,2), intent(in) :: &
    stf                 ! surface forcing for all tracers

  real(r8), dimension(nx,ny), intent(in) :: &
    shf_qsw             ! short-wave forcing

  real(r8), dimension(nx,ny,2), intent(in), optional :: &
    smf,               &! surface momentum forcing at U points
    smft                ! surface momentum forcing at T points
                        ! *** either one or the other (not
                        ! *** both) should be passed
  real(r8), intent(in) ::&
    convect_diff,      &! diffusivity to mimic convection
    convect_visc        ! viscosity   to mimic convection

  ! !INPUT/OUTPUT PARAMETERS:

  real(r8), dimension(nx,ny,nz), intent(inout) ::      &
    vvc        ! viscosity for momentum diffusion

  real(r8), dimension(nx,ny,0:nzf,2),intent(inout) :: &
    vdc        ! diffusivity for tracer diffusion

  real(r8), dimension(nx,ny,nz,2), intent(out) ::    &
    kpp_src

  real(r8), dimension(nx,ny), intent(out) :: &
    kpp_hblt

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer  :: &
    k,                 &! vertical level index
    i,j,               &! horizontal loop indices
    n,                 &! tracer index
    mt2                 ! index for separating temp from other trcrs

  integer, dimension(nx,ny) :: &
    kbl                   ! index of first lvl below hbl

  real(r8), dimension(nx,ny) :: &
    ustar,      &! surface friction velocity
    bfsfc,      &! surface buoyancy forcing
    work1,work2,&! temporary storage
    fcon,       &! convection temporary
    stable       ! = 1 for stable forcing; = 0 for unstable forcing

  real(r8), dimension(nx,ny,nz) :: &
    dbloc,      &! buoyancy difference between adjacent levels
    dbsfc,      &! buoyancy difference between level and surface
    ghat         ! non-local mixing coefficient

  real(r8), dimension(nx,ny,0:nzf) :: &
    visc        ! local temp for viscosity
  !-----------------------------------------------------------------------
  !  initialize
  !-----------------------------------------------------------------------
  !  if(lkpp_ini) then
  !    call init_vmix_kpp
  !    lkpp_ini = .false.
  !  end if

  !$acc enter data async(async_id) &
  !$acc& create(ustar, bfsfc, work1, work2, fcon, stable, dbloc, dbsfc, ghat, visc, kbl)

  !-----------------------------------------------------------------------
  !  compute buoyancy differences at each vertical level.
  !  derived variables: dbloc, dbsfc
  !-----------------------------------------------------------------------

  call buoydiff_gpu(dbloc, dbsfc, trcr, async_id)

  !-----------------------------------------------------------------------
  !  compute mixing due to shear instability, internal waves and
  !  convection
  !  derived variables: vdc,visc
  !-----------------------------------------------------------------------

  call ri_iwmix_gpu(dbloc, visc, vdc, uuu, vvv, async_id)

  !-----------------------------------------------------------------------
  !  compute double diffusion if desired
  !-----------------------------------------------------------------------
  if(ldbl_diff)  call ddmix_gpu(vdc,trcr, async_id)

  !-----------------------------------------------------------------------
  !  compute boundary layer depth
  !  derived variables: kpp_hblt, ustar, bfsfc, stable, kbl
  !-----------------------------------------------------------------------

  if(present(smft)) then
    call bldepth_gpu(trans,dbloc, dbsfc, trcr, uuu, vvv, stf, shf_qsw, &
                 kpp_hblt, ustar, bfsfc, stable, kbl, &
                 smfx = smft, is_smft = .true., async_id=async_id)
  else
    call bldepth_gpu(trans,dbloc, dbsfc, trcr, uuu, vvv, stf, shf_qsw, &
                 kpp_hblt, ustar, bfsfc, stable, kbl, &
                 smfx = smf, is_smft = .false., async_id=async_id)
  end if

  !-----------------------------------------------------------------------
  !  compute boundary layer diffusivities
  !  derived variables: visc,vdc,ghat
  !-----------------------------------------------------------------------

  call blmix_gpu(visc, vdc, kpp_hblt, ustar, bfsfc, stable, &
             kbl, ghat, async_id)

  ! #ifdef chk_kpp1
  !   if(myid_x .eq. 7 .and. myid_y .eq. 11) then
  !     write(*,*) "[KPP CHECK vmix]: 01"
  !     write(*,*) "kbl = ", kbl(nx/2,ny)
  !     write(*,*) "hblt = ", kpp_hblt(nx/2,ny)
  !     write(*,*) achar(10)
  !   end if
  ! #endif
  !-----------------------------------------------------------------------
  !
  !  consider interior convection:
  !
  !    compute function of Brunt-Vaisala squared for convection.
  !
  !  use either a smooth
  !
  !    WORK1 = N**2,  FCON is function of N**2
  !    FCON = 0 for N**2 > 0
  !    FCON = [1-(1-WORK1/BVSQcon)**2]**3 for BVSQcon < N**2 < 0
  !    FCON = 1 for N**2 < BVSQcon
  !
  !  or a step function. The smooth function has been used with
  !  BVSQcon = -0.2e-4_dbl_kind.
  !
  !  after convection, average viscous coeffs to U-grid and reset sea
  !  floor values
  !
  !-----------------------------------------------------------------------

  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz-1
    do j = 1, ny
      do i = 1, nx
        work1(i, j) = dbloc(i,j,k)/(zgrid(k) - zgrid(k+1))
        if(BVSQcon /= 0.d0) then
          work2(i, j) = min(1.d0-(max(work1(i, j),BVSQcon))/BVSQcon, 1.d0)
          fcon(i, j)  = (1.d0 - work2(i, j)*work2(i, j))**3
        else
          if (work1(i, j) > 0.d0) then
            fcon(i, j) = 0.d0
          else
            fcon(i, j) = 1.d0
          end if
        end if

        !*** add convection and reset sea floor values to zero
        if(k >= kbl(i,j)) then
          visc(i,j,k)   = visc(i,j,k)  + convect_visc * fcon(i,j)
           vdc(i,j,k,1) = vdc(i,j,k,1) + convect_diff * fcon(i,j)
           vdc(i,j,k,2) = vdc(i,j,k,2) + convect_diff * fcon(i,j)
        end if
        if(k >= kb(i,j)) then
          visc(i,j,k  ) = 0.d0
          vdc (i,j,k,1) = 0.d0
          vdc (i,j,k,2) = 0.d0
        end if

        if (k < kb(i, j)) then
          vvc(i, j, k) = visc(i, j, k)
        else
          vvc(i, j, k) = 0.d0
        end if

        vdc(i,j,nz,:) = 0.d0
        vvc(i,j,nz)   = 0.d0
      end do
    end do
  end do

  ! #ifdef chk_kpp1
  !   if(myid_x .eq. 7 .and. myid_y .eq. 11) then
  !     write(*,*) "[KPP CHECK vmix]: 02"
  !     write(*,*) "kbl = ", kbl(nx/2,ny)
  !     write(*,*) "hblt = ", kpp_hblt(nx/2,ny)
  !     write(*,*) achar(10)
  !   end if
  ! #endif
  !-----------------------------------------------------------------------
  !  add ghatp term from previous computation to right-hand-side
  !  source term on current row
  !-----------------------------------------------------------------------
  !$acc parallel loop collapse(4) async(async_id)
  do n = 1, 2
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          mt2 = min(n,2)
          if (k == 1) then
            kpp_src(i,j,k,n) = stf(i,j,n)*odz(1)              &
                            *(-vdc(i,j,1,mt2)*ghat(i,j,1))
          else
            kpp_src(i,j,k,n) = stf(i,j,n)*odz(k)              &
                            *( vdc(i,j,k-1,mt2)*ghat(i,j,k-1) &
                              -vdc(i,j,k  ,mt2)*ghat(i,j,k  ))
          end if
          if (n == 2) then
            kpp_src(i,j,k,2) = 1000.d0*kpp_src(i,j,k,2) !from msu to psu
          end if
        end do
      end do
      ! #ifdef chk_kpp1
      !   if(myid_x .eq. 7 .and. myid_y .eq. 11) then
      !     write(*,*) "[KPP CHECK vmix]: 02-1", n, k
      !     write(*,*) "kbl = ", kbl(nx/2,ny)
      !     write(*,*) "hblt = ", kpp_hblt(nx/2,ny)
      !     write(*,*) achar(10)
      !   end if
      ! #endif
    end do
  end do

  !-----------------------------------------------------------------------
  !
  !  compute diagnostic mixed layer depth (cm) using a max buoyancy
  !  gradient criterion.  Use USTAR and BFSFC as temps.
  !
  !-----------------------------------------------------------------------
  ! #ifdef chk_kpp1
  !   if(myid_x .eq. 7 .and. myid_y .eq. 11) then
  !     write(*,*) "[KPP CHECK vmix]: 03"
  !     write(*,*) "kbl = ", kbl(nx/2,ny)
  !     write(*,*) "hblt = ", kpp_hblt(nx/2,ny)
  !     write(*,*) achar(10)
  !   end if
  ! #endif

  !$acc parallel loop collapse(2) async(async_id)
  do j = 1, ny
    do i = 1, nx
      ustar(i,j) = 0.d0
      if (kb(i,j) == 1) then
        hmxl(i,j) = z_grid(1)-z_face(1)
      else
        hmxl(i,j) = 0.d0
      end if
      !$acc loop seq
      do k = 2, nz
        if (k <= kb(i,j)) then
          ustar(i,j) = max(dbsfc(i,j,k)/(z_grid(k)-z_face(1)), ustar(i,j))
          hmxl(i,j) = z_grid(k)-z_face(1)
        end if
      end do
    end do
  end do

  !$acc parallel loop collapse(3) async(async_id)
  do j = 1, ny
    do i = 1, nx
      do k = 2, nz
        visc(i,j,1) = 0.d0
        if (ustar(i,j) > 0.d0) then
          visc(i,j,k) = (dbsfc(i,j,k)-dbsfc(i,j,k-1))/ &
                        (z_grid(k) - z_grid(k-1))
        end if
        if (visc(i,j,k) >= ustar(i,j) .and. &
           (visc(i,j,k)-visc(i,j,k-1)) /= 0.d0 .and. &
            ustar(i,j) > 0.d0) then ! avoid divide by zero
          bfsfc(i,j) = (visc(i,j,k)-ustar(i,j))/ &
                       (visc(i,j,k)-visc(i,j,k-1))
          hmxl(i,j) = -0.5d0*(zgrid(k  ) + zgrid(k-1))*(1.0d0-bfsfc(i,j)) &
                      -0.5d0*(zgrid(k-1) + zgrid(k-2))*bfsfc(i,j)
          ustar(i,j) = 0.d0
        end if
      end do
    end do
  end do

  ! #ifdef chk_kpp1
  !   if(myid_x .eq. 7 .and. myid_y .eq. 11) then
  !     write(*,*) "[KPP CHECK vmix]: 04"
  !     write(*,*) "hblt = ", kpp_hblt(nx/2,ny)
  !     write(*,*) achar(10)
  !   end if
  ! #endif

  !$acc exit data async(async_id) &
  !$acc& delete(ustar, bfsfc, work1, work2, fcon, stable, dbloc, dbsfc, ghat, visc, kbl)
end subroutine vmix_coeffs_kpp_gpu

subroutine buoydiff_gpu(dbloc, dbsfc, trcr, async_id)
  !$acc routine(state_gpu) seq
  use timcom_kpp, only: zgrid
  implicit none
  ! !DESCRIPTION:
  !  This routine calculates the buoyancy differences at model levels.
  !
  ! !REVISION HISTORY:
  !  same as module
  !   implicit none
  ! !INPUT PARAMETERS:
  integer, intent(in) :: async_id
  real(r8), dimension(nx,ny,nz,2), intent(inout) :: &
    trcr     ! tracers at current time
  ! !OUTPUT PARAMETERS:
  real(r8), dimension(nx,ny,nz), intent(out) :: &
    dbloc, & ! buoyancy difference between adjacent levels
    dbsfc    ! buoyancy difference between level and surface

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer :: &
    k,                 &! vertical level index
    i,j,               &! horizontal indices
    kprev, klvl, ktmp   ! indices for 2-level TEMPK array
  integer, dimension(nx,ny) :: kmt
  real(r8) , dimension(nx,ny) :: &
    rho1,              &! density of sfc t,s displaced to k
    rhokm,             &! density of t(k-1),s(k-1) displaced to k
    rhok,              &! density at level k
    tempsfc             ! adjusted temperature at surface
  real(r8) , dimension(nx,ny,2) :: &
    tempk               ! temp adjusted for freeze at levels k,k-1
  real(r8), dimension(nx,ny) :: test,test2,outs
  real(r8) :: tmp

  !$acc enter data async(async_id) &
  !$acc& create(rho1, rhokm, rhok, tempsfc, tempk)

  !$acc parallel loop collapse(2) private(kprev, klvl, ktmp) async(async_id)
  do j = 1, ny
    do i = 1, nx
      !-----------------------------------------------------------------------
      !  calculate density and buoyancy differences at surface
      !-----------------------------------------------------------------------
      if (trcr(i,j,1,1) < -2.d0) then
        tempsfc(i,j) = -2.d0
      else
        tempsfc(i,j) = trcr(i,j,1,1)
      endif
      klvl  = 2
      kprev = 1
      tempk(i,j,kprev) = tempsfc(i,j)
      dbsfc(i,j,1) = 0.d0
      !-----------------------------------------------------------------------
      !  calculate DBLOC and DBSFC for all other levels
      !-----------------------------------------------------------------------
      !$acc loop seq
      do k = 2, nz
        if (trcr(i,j,k,1) < -2.0d0) then
          tempk(i,j,klvl) = -2.0d0
        else
          tempk(i,j,klvl) = trcr(i,j,k,1)
        endif

        ! calculate the density from Temperature and Salinity using an equation of state
        ! derived from McDougall, Wright, Jackett and Feistel (hereafter MWJF, 2001 )
        call state_gpu(i, j, k, k, tempsfc(i,j),          trcr(i,j,1,  2), &
          tmp, rho1(i,j), tmp, tmp, .false., .true., .false., .false., zgrid(k), async_id)
        call state_gpu(i, j, k, k, tempk(i,j,kprev), trcr(i,j,k-1,2), &
          tmp, rhokm(i,j), tmp, tmp, .false., .true., .false., .false., zgrid(k), async_id)
        call state_gpu(i, j, k, k, tempk(i,j,klvl),  trcr(i,j,k  ,2), &
          tmp, rhok(i,j), tmp, tmp, .false., .true., .false., .false., zgrid(k), async_id)

        if(rhok(i,j) /= 0.d0) then
          dbsfc(i,j,k)   = grav*(1.0d0 - rho1 (i,j)/rhok(i,j))
          dbloc(i,j,k-1) = grav*(1.0d0 - rhokm(i,j)/rhok(i,j))
        else
          dbsfc(i,j,k)   = 0.d0
          dbloc(i,j,k-1) = 0.d0
        end if
        if ( k-1 >= kb(i,j) ) dbloc(i,j,k-1) = 0.d0
        ktmp  = klvl
        klvl  = kprev
        kprev = ktmp
      end do
      dbloc(i,j,nz) = 0.d0
    end do
  end do

  !$acc exit data async(async_id) &
  !$acc& delete(rho1, rhokm, rhok, tempsfc, tempk)

end subroutine buoydiff_gpu

function pressure(depth)
  !$acc routine seq
! !INPUT PARAMETERS:
  real(r8), intent(in) :: depth    ! depth in meters
! !OUTPUT PARAMETERS:
  real(r8) :: pressure   ! pressure [bars]

  pressure = 0.059808d0*(dexp(-0.025d0*depth) - 1.0d0) &
           + 0.100766d0*depth + 2.28405d-7*depth**2.d0
end function pressure

subroutine state_gpu(i, j, k, kk, tempk, saltk, &
                 rhoout, rhofull, drhodt, drhods, &
                 calc_rhoout, calc_rhofull, calc_drhodt, calc_drhods, zgrid, async_id)
  !$acc routine seq
  !$acc routine(pressure) seq
  implicit none
  ! !DESCRIPTION:
  !  Returns the density of water at level k from equation of state
  !  $\rho = \rho(d,\theta,S)$ where $d$ is depth, $\theta$ is
  !  potential temperature, and $S$ is salinity. the density can be
  !  returned as a perturbation (RHOOUT) or as the full density
  !  (RHOFULL). Note that only the polynomial EOS choice will return
  !  a perturbation density; in other cases the full density is returned
  !  regardless of which argument is requested.
  !
  !  This routine also computes derivatives of density with respect
  !  to temperature and salinity at level k from equation of state
  !  if requested (ie the calc_* flags are set to .true.).
  !
  !  If $k = kk$ are equal the density for level k is returned.
  !  If $k \neq kk$ the density returned is that for a parcel
  !  adiabatically displaced from level k to level kk.

  ! !INPUT PARAMETERS:
  integer, intent(in) :: &
    i, j, &
    k,                   &! depth level index
    kk,                  &! level to which water is adiabatically
                          ! displaced
    async_id              ! GPU stream ID for asynchronous execution
  real(r8), intent(in) :: &
    tempk,               &! temperature at level k
    saltk                 ! salinity    at level k

  logical, intent(in) :: &
    calc_rhoout,        &! flag to calculate perturbation density
    calc_rhofull,       &! flag to calculate full density
    calc_drhodt,        &! flag to calculate temperature derivative
    calc_drhods          ! flag to calculate salinity derivative

  real(r8), intent(in) :: &
    zgrid

  ! !OUTPUT PARAMETERS:
  real(r8), intent(out) :: &
    rhoout,  &! perturbation density of water
    rhofull, &! full density of water
    drhodt,  &! derivative of density with respect to temperature
    drhods    ! derivative of density with respect to salinity
  !EOP
  !BOC
  !-----------------------------------------------------------------------
  !  local variables:
  !-----------------------------------------------------------------------
  integer :: dep
  real(r8) :: &
    TQ, SQ,             &! adjusted T,S
    sqr, denomk,        &! work values
    work1, work2, work3, work4
  real(r8) :: ppz ! temporary pressure scalars
  real(r8) ::      &
    mwjfnums0t0, mwjfnums0t1, mwjfnums0t2, mwjfnums0t3,              &
    mwjfnums1t0, mwjfnums1t1, mwjfnums2t0,                           &
    mwjfdens0t0, mwjfdens0t1, mwjfdens0t2, mwjfdens0t3, mwjfdens0t4, &
    mwjfdens1t0, mwjfdens1t1, mwjfdens1t3,                           &
    mwjfdensqt0, mwjfdensqt2
  !-----------------------------------------------------------------------
  !  valid ranges and pressure as function of depth
  !-----------------------------------------------------------------------
  real(r8), dimension(nz) :: &
    pressz              ! ref pressure (bars) at each level
  !-----------------------------------------------------------------------
  !  MWJF EOS coefficients
  !-----------------------------------------------------------------------
  !*** these constants will be used to construct the numerator
  !*** factor unit change (kg/m^3 -> g/cm^3) into numerator terms
  real(r8), parameter :: p001 = 0.001d0
  real(r8), parameter ::                   &
    mwjfnp0s0t0 =  9.99843699d02*p001, &
    mwjfnp0s0t1 =  7.35212840d+0*p001, &
    mwjfnp0s0t2 = -5.45928211d-2*p001, &
    mwjfnp0s0t3 =  3.98476704d-4*p001, &
    mwjfnp0s1t0 =  2.96938239d+0*p001, &
    mwjfnp0s1t1 = -7.23268813d-3*p001, &
    mwjfnp0s2t0 =  2.12382341d-3*p001, &
    mwjfnp1s0t0 =  1.04004591d-2*p001, &
    mwjfnp1s0t2 =  1.03970529d-7*p001, &
    mwjfnp1s1t0 =  5.18761880d-6*p001, &
    mwjfnp2s0t0 = -3.24041825d-8*p001, &
    mwjfnp2s0t2 = -1.23869360d-11*p001

  !*** these constants will be used to construct the denominator
  real(r8), parameter ::            &
    mwjfdp0s0t0 =   1.0d+0,         &
    mwjfdp0s0t1 =   7.28606739d-3,  &
    mwjfdp0s0t2 =  -4.60835542d-5,  &
    mwjfdp0s0t3 =   3.68390573d-7,  &
    mwjfdp0s0t4 =   1.80809186d-10, &
    mwjfdp0s1t0 =   2.14691708d-3,  &
    mwjfdp0s1t1 =  -9.27062484d-6,  &
    mwjfdp0s1t3 =  -1.78343643d-10, &
    mwjfdp0sqt0 =   4.76534122d-6,  &
    mwjfdp0sqt2 =   1.63410736d-9,  &
    mwjfdp1s0t0 =   5.30848875d-6,  &
    mwjfdp2s0t3 =  -3.03175128d-16, &
    mwjfdp3s0t1 =  -1.27934137d-17

  !*** prevent problems with garbage on land points or ghost cells
  TQ = min(tempk, 999.0d0)
  TQ = max(TQ, -2.0d0)
  SQ = min(saltk, 0.999d0)
  SQ = max(SQ, 0.d0)

  !------------------------------------------------------------
  ppz = 10.d0*pressure(-zgrid*0.01d0) ! cm to m

  SQ = 1000.d0*SQ
  SQR = dsqrt(SQ)

  !***
  !*** first calculate numerator of MWJF density [P_1(S,T,p)]
  !***
  mwjfnums0t0 = mwjfnp0s0t0 + ppz*(mwjfnp1s0t0 + ppz*mwjfnp2s0t0)
  mwjfnums0t1 = mwjfnp0s0t1
  mwjfnums0t2 = mwjfnp0s0t2 + ppz*(mwjfnp1s0t2 + ppz*mwjfnp2s0t2)
  mwjfnums0t3 = mwjfnp0s0t3
  mwjfnums1t0 = mwjfnp0s1t0 + ppz*mwjfnp1s1t0
  mwjfnums1t1 = mwjfnp0s1t1
  mwjfnums2t0 = mwjfnp0s2t0

  WORK1 = mwjfnums0t0 + TQ * (mwjfnums0t1 + TQ * (mwjfnums0t2 + &
          mwjfnums0t3 * TQ)) + SQ * (mwjfnums1t0 +              &
          mwjfnums1t1 * TQ + mwjfnums2t0 * SQ)

  !***
  !*** now calculate denominator of MWJF density [P_2(S,T,p)]
  !***
  mwjfdens0t0 = mwjfdp0s0t0 + ppz*mwjfdp1s0t0
  mwjfdens0t1 = mwjfdp0s0t1 + ppz**3 * mwjfdp3s0t1
  mwjfdens0t2 = mwjfdp0s0t2
  mwjfdens0t3 = mwjfdp0s0t3 + ppz**2 * mwjfdp2s0t3
  mwjfdens0t4 = mwjfdp0s0t4
  mwjfdens1t0 = mwjfdp0s1t0
  mwjfdens1t1 = mwjfdp0s1t1
  mwjfdens1t3 = mwjfdp0s1t3
  mwjfdensqt0 = mwjfdp0sqt0
  mwjfdensqt2 = mwjfdp0sqt2

  WORK2 = mwjfdens0t0 + TQ * (mwjfdens0t1 + TQ * (mwjfdens0t2 +    &
          TQ * (mwjfdens0t3 + mwjfdens0t4 * TQ))) +                   &
          SQ * (mwjfdens1t0 + TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3)+ &
          SQR * (mwjfdensqt0 + TQ*TQ*mwjfdensqt2))

  DENOMK = 1.d0/WORK2

  if(calc_rhoout) then
    RHOOUT = WORK1*DENOMK
  end if

  if(calc_rhofull) then
    RHOFULL = WORK1*DENOMK
  end if

  if(calc_drhodt) then
    WORK3 = mwjfnums0t1 + TQ * (2.d0*mwjfnums0t2 +    &
            3.d0*mwjfnums0t3 * TQ) + mwjfnums1t1 * SQ
    WORK4 = mwjfdens0t1 + SQ * mwjfdens1t1 +               &
            TQ * (2.d0*(mwjfdens0t2 + SQ*SQR*mwjfdensqt2) +  &
            TQ * (3.d0*(mwjfdens0t3 + SQ * mwjfdens1t3) +    &
            TQ *  4.d0*mwjfdens0t4))
    DRHODT = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK
  end if

  if(calc_drhods) then
    WORK3 = mwjfnums1t0 + mwjfnums1t1 * TQ + 2.d0*mwjfnums2t0 * SQ
    WORK4 = mwjfdens1t0 +   &
            TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3) +   &
            1.5d0*SQR*(mwjfdensqt0 + TQ*TQ*mwjfdensqt2)
    DRHODS = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK * 1000.d0
  end if

end subroutine state_gpu

subroutine ri_iwmix_gpu(dbloc, visc, vdc, uuu, vvv, async_id)
  ! !DESCRIPTION:
  !  Computes viscosity and diffusivity coefficients for the interior
  !  ocean due to shear instability (richardson number dependent),
  !  internal wave activity, and to static instability (Ri < 0).
  use timcom_kpp, only: tidal_diff, tidal_coef, bckgrnd_vvc, bckgrnd_vdc, &
                        zgrid, num_v_smooth_Ri, ltidal_mixing, lccsm_control_compatible, &
                        prandtl, tidal_mix_max, lrich, Riinfty, rich_mix
  implicit none
  ! !INPUT PARAMETERS:
  integer, intent(in) :: async_id
  real(r8), dimension(-1:nx+2,-1:ny+2,nz), intent(in) :: &
    uuu, &     ! U velocities at current time
    vvv        ! V velocities at current time
  real(r8), dimension(nx,ny,nz), intent(in) :: &
    dbloc      ! buoyancy difference between adjacent levels
  ! !INPUT/OUTPUT PARAMETERS:
  real(r8), dimension(nx,ny,0:nzf,2), intent(inout) :: &
    vdc        ! diffusivity for tracer diffusion
  ! !OUTPUT PARAMETERS:
  real(r8), dimension(nx,ny,0:nzf), intent(out) :: &
    visc       ! viscosity
  !EOP
  !BOC
  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer :: &
    k,                 &! index for vertical levels
    i,j,               &! horizontal loop indices
    n                   ! vertical smoothing index
  real(r8), dimension(nx,ny) :: &
    vshear,            &! (local velocity shear)^2
    ri_loc,            &! local Richardson number
    fri,               &! function of Ri for shear
    work1
  real(r8), dimension(nx,ny) :: kvmix,kvmix_m
  real(r8) :: temp_val

  !$acc enter data async(async_id) &
  !$acc& create(kvmix, kvmix_m, vshear, ri_loc, fri, work1)

  !$acc parallel loop collapse(2) async(async_id)
  do j = 1, ny
    do i = 1, nx
      !-----------------------------------------------------------------------
      !  compute mixing at each level
      !-----------------------------------------------------------------------
      ! Initialize arrays with loops
      kvmix(i,j) = 0.d0
      kvmix_m(i,j) = 0.d0
      visc(i,j,0) = 0.d0

      !$acc loop seq
      do k = 1, nz
        !-----------------------------------------------------------------------
        !     compute velocity shear squared and average to T points:
        !     VSHEAR = (UUU(k)-UUU(k+1))**2+(VVV(k)-VVV(k+1))**2
        !     Use FRI here as a temporary.
        !-----------------------------------------------------------------------
        if(k < nz) then
          fri(i,j) = (uuu(i,j,k)-uuu(i,j,k+1))**2 + &
                     (vvv(i,j,k)-vvv(i,j,k+1))**2
          vshear(i,j) = fri(i,j)
        else
          vshear(i,j) = 0.d0
        end if
        !-----------------------------------------------------------------------
        !     compute local richardson number
        !     use visc array as temporary Ri storage to be smoothed
        !-----------------------------------------------------------------------
        ri_loc(i,j) = dbloc(i,j,k)*(zgrid(k)-zgrid(k+1))/(vshear(i,j) + eps)
        if (k <= int(kb(i,j))) then
          visc(i,j,k) = ri_loc(i,j)
        else
          visc(i,j,k) = visc(i,j,k-1)
        end if
      end do

      !-----------------------------------------------------------------------
      !  vertically smooth Ri num_v_smooth_Ri times with 1-2-1 weighting
      !  result again stored temporarily in VISC and use RI_LOC and FRI
      !  as temps
      !----------------------------------------------------------------------
      !$acc loop seq
      do n = 1, num_v_smooth_Ri
        fri(i,j) = 0.25d0*visc(i,j,1)
        visc(i,j,nz+1) = visc(i,j,nz)

        !$acc loop seq
        do k = 1, nz
          ri_loc(i,j) = visc(i,j,k)
          if(int(kb(i,j)) >= 2 ) then !modified from KMT, land=0
              visc(i,j,k) = fri(i,j) + 0.5d0*ri_loc(i,j) &
                                      + 0.25d0*visc(i,j,k+1)
          end if
          fri(i,j) = 0.25d0*ri_loc(i,j)
        end do
      end do
      !-----------------------------------------------------------------------
      !  now that we have a smoothed Ri field, finish computing coeffs
      !  at each level
      !-----------------------------------------------------------------------
      !$acc loop seq
      do k = 1, nz
        if(ltidal_mixing) then
          tidal_diff(i,j,k) = 0.d0
          !-----------------------------------------------------------------------
          !
          !     if Ri-number mixing requested,
          !     evaluate function of Ri for shear instability:
          !       for 0 < Ri < Riinfty, function = (1 - (Ri/Riinfty)**2)**3
          !       for     Ri > Riinfty, function = 0
          !       for     Ri < 0      , function = 1
          !     compute contribution due to shear instability
          !     VISC holds smoothed Ri at k, but replaced by real VISC
          !
          !     otherwise only use iw
          !     convection is added later
          !
          !-----------------------------------------------------------------------
          !-----------------------------------------------------------------------
          !
          !  consider the internal wave mixing first. rich_mix is used as the
          !  upper limit for internal wave mixing coefficient. bckgrnd_vvc
          !  was already multiplied by Prandtl.
          !
          !  NOTE: no partial_bottom_cell implementation at this time
          !
          !-----------------------------------------------------------------------
          work1(i,j) = dbloc(i,j,k)/(zgrid(k) - zgrid(k+1))
          if (work1(i,j) > 0.d0) then
            tidal_diff(i,j,k) = tidal_coef(i,j,k)/work1(i,j)
          end if

          ! Notes:
          ! (1) this step breaks backwards compatibility
          ! (2) check for k>2 was added to if statement to
          ! avoid out of bounds access
          if(.not. lccsm_control_compatible) then   ! this step breaks backwards compatibility
            if(k > 2) then
              if(k == int(kb(i,j))-1 .or. k == int(kb(i,j))-2) then
                tidal_diff(i,j,k) = max(tidal_diff(i,j,k), tidal_diff(i,j,k-1))
              end if
            end if
          end if

          temp_val = bckgrnd_vvc(i,j,k)/Prandtl + tidal_diff(i,j,k)
          work1(i,j) = Prandtl * min(temp_val, tidal_mix_max)

          if(k < nz) then
            kvmix_m(i,j) = work1(i,j)
          end if

          if(k < nz) then
            vdc(i,j,k,2) = min(bckgrnd_vdc(i,j,k) + tidal_diff(i,j,k), tidal_mix_max)
            kvmix(i,j) = vdc(i,j,k,2)
          end if

          if(lrich) then
            temp_val = max(visc(i,j,k), 0.d0)/Riinfty
            fri(i,j) = min(temp_val, 1.d0)
            visc(i,j,k) = work1(i,j) + rich_mix*(1.d0 - fri(i,j)*fri(i,j))**3

            if(k < nz) then
              vdc(i,j,k,2) = vdc(i,j,k,2) + rich_mix*(1.d0 - fri(i,j)*fri(i,j))**3
              vdc(i,j,k,1) = vdc(i,j,k,2)
            end if
          else
            visc(i,j,k) = work1(i,j)

            if(k < nz) then
              vdc(i,j,k,1) = vdc(i,j,k,2)
            end if
          end if
        else ! .not. ltidal_mixing
          if(k < nz) then
            kvmix(i,j) = bckgrnd_vdc(i,j,k)
            kvmix_m(i,j) = bckgrnd_vvc(i,j,k)
          end if
          if(lrich) then
            temp_val = max(visc(i,j,k), 0.d0)/Riinfty
            fri(i,j) = min(temp_val, 1.0d0)
            visc(i,j,k) = bckgrnd_vvc(i,j,k) + rich_mix*(1.0d0 - fri(i,j)*fri(i,j))**3
            if(k < nz) then
              vdc(i,j,k,2) = bckgrnd_vdc(i,j,k) + rich_mix*(1.0d0 - fri(i,j)*fri(i,j))**3
              vdc(i,j,k,1) = vdc(i,j,k,2)
            end if
          else
            visc(i,j,k) = bckgrnd_vvc(i,j,k)
            if(k < nz) then
              vdc(i,j,k,2) = bckgrnd_vdc(i,j,k)
              vdc(i,j,k,1) = vdc(i,j,k,2)
            end if
          end if
        end if ! ltidal_mixing

        !-----------------------------------------------------------------------
        !     set seafloor values to zero
        !-----------------------------------------------------------------------
        if(k >= int(kb(i,j))) then
          visc(i,j,k  ) = 0.d0
          vdc (i,j,k,1) = 0.d0
          vdc (i,j,k,2) = 0.d0
        end if
      end do
      !-----------------------------------------------------------------------
      !  fill extra coefficients for blmix
      !-----------------------------------------------------------------------
      visc(i,j,0) = 0.d0
      vdc(i,j,0,1) = 0.d0
      vdc(i,j,0,2) = 0.d0
      visc(i,j,nzf) = 0.d0
      vdc(i,j,nzf,1) = 0.d0
      vdc(i,j,nzf,2) = 0.d0
    end do
  end do

  !$acc exit data async(async_id) &
  !$acc& delete(kvmix, kvmix_m, vshear, ri_loc, fri, work1)

end subroutine ri_iwmix_gpu

subroutine ddmix_gpu(vdc, trcr, async_id)
  ! !DESCRIPTION:
  !  $R_\rho$ dependent interior flux parameterization.
  !  Add double-diffusion diffusivities to Ri-mix values at blending
  !  interface and below.
  !$acc routine(state_gpu) seq
  use timcom_kpp, only: rrho0, dsfmax, zgrid
  implicit none
  !!INPUT PARAMETERS:
  integer, intent(in) :: async_id
  real(r8), dimension(nx,ny,nz,2), intent(in) :: &
    trcr                ! tracers at current time
  ! !INPUT/OUTPUT PARAMETERS:
  real(r8), dimension(nx,ny,0:nzf,2),intent(inout) :: &
    vdc        ! diffusivity for tracer diffusion
  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer ::  k, i, j, kup, knxt
  real(r8), dimension(nx,ny) :: &
    alphadt,           &! alpha*DT  across interfaces
    betads,            &! beta *DS  across interfaces
    rrho,              &! dd density ratio
    diffdd,            &! dd diffusivity scale
    prtl_tmp             ! prandtl number
  real(r8), dimension(nx,ny,2) :: &
    talpha,            &! temperature expansion coeff
    sbeta               ! salinity    expansion coeff
  real(r8) :: tmp

  !-----------------------------------------------------------------------
  !  compute alpha*DT and beta*DS at interfaces.  use RRHO and
  !  PRANDTL for temporary storage for call to state
  !-----------------------------------------------------------------------

  !$acc enter data async(async_id) &
  !$acc& create(alphadt, betads, rrho, diffdd, prtl_tmp, talpha, sbeta, tmp)

  !$acc parallel loop collapse(2) private(kup, knxt) async(async_id)
  do j = 1, ny
    do i = 1, nx
      kup = 1
      knxt = 2
      prtl_tmp(i,j) = merge(-2.0d0, trcr(i,j,1,1), trcr(i,j,1,1) < -2.0d0)
      call state_gpu(i, j, 1, 1, prtl_tmp(i,j), trcr(i,j,1,2), &
        tmp, rrho(i,j), talpha(i,j,kup), sbeta(i,j,kup), .false., .true., .true., .true., zgrid(1), async_id)

      !$acc loop seq
      do k = 1, nz
        if(k < nz) then
          prtl_tmp(i,j) = merge(-2.0d0, trcr(i,j,k+1,1), trcr(i,j,k+1,1) < -2.0d0)
          call state_gpu(i, j, k+1, k+1, prtl_tmp(i,j), trcr(i,j,k+1,2), &
            tmp, rrho(i,j), talpha(i,j,knxt), sbeta(i,j,knxt), .false., .true., .true., .true., zgrid(k+1), async_id)
          alphadt(i,j) = -0.5d0*(talpha(i,j,kup) + talpha(i,j,knxt)) &
                               *(trcr(i,j,k,1) - trcr(i,j,k+1,1))

          betads(i,j) = 0.5d0*(sbeta(i,j,kup) + sbeta(i,j,knxt)) &
                              *(trcr(i,j,k,2) - trcr(i,j,k+1,2))
          kup  = knxt
          knxt = 3 - kup
        else
          alphadt(i,j) = 0.d0
          betads(i,j) = 0.d0
        end if
        !-----------------------------------------------------------------------
        !     salt fingering case
        !-----------------------------------------------------------------------
        if (alphadt(i,j) > betads(i,j) .and. betads(i,j) > 0.d0) then
          rrho(i,j) = min(alphadt(i,j)/betads(i,j), rrho0)
          diffdd(i,j) = dsfmax*(1.0d0-(rrho(i,j)-1.0d0)/(rrho0-1.0d0))**3
          vdc(i,j,k,1) = vdc(i,j,k,1) + 0.7d0*diffdd(i,j)
          vdc(i,j,k,2) = vdc(i,j,k,2) + diffdd(i,j)
        end if
        !-----------------------------------------------------------------------
        !     diffusive convection
        !-----------------------------------------------------------------------
        if (alphadt(i,j) < 0.d0 .and. betads(i,j) < 0.d0 .and. alphadt(i,j) > betads(i,j)) then
          rrho(i,j) = alphadt(i,j) / betads(i,j)
          diffdd(i,j) = 1.5d-2*0.909d0* &
              dexp(4.6d0*dexp(-0.54d0*(1.0d0/rrho(i,j)-1.0d0)))
          prtl_tmp(i,j) = 0.15d0*rrho(i,j)
        else
          rrho(i,j) = 0.d0
          diffdd(i,j) = 0.d0
          prtl_tmp(i,j) = 0.d0
        end if

        if (rrho(i,j) > 0.5d0) then
          prtl_tmp(i,j) = (1.85d0 - 0.85d0/rrho(i,j))*rrho(i,j)
        end if

        vdc(i,j,k,1) = vdc(i,j,k,1) + diffdd(i,j)
        vdc(i,j,k,2) = vdc(i,j,k,2) + prtl_tmp(i,j)*diffdd(i,j)
      end do
    end do
  end do

  !$acc exit data async(async_id) &
  !$acc& delete(alphadt, betads, rrho, diffdd, prtl_tmp, talpha, sbeta, tmp)

end subroutine ddmix_gpu

subroutine bldepth_gpu(trans, dbloc, dbsfc, trcr, uuu, vvv, stf, shf_qsw,&
                    hblt, ustar, bfsfc, stable, kbl, smfx, is_smft, async_id)
  !$acc routine(wscale_gpu) seq
  !$acc routine(state_gpu) seq
  ! !DESCRIPTION:
  !  This routine computes the ocean boundary layer depth defined as
  !  the shallowest depth where the bulk Richardson number is equal to
  !  a critical value, Ricr.
  !
  !  NOTE: bulk richardson numbers are evaluated by computing
  !        differences between values at zgrid(kl) $< 0$ and surface
  !        reference values. currently, the reference values are equal
  !        to the values in the surface layer.  when using higher
  !        vertical grid resolution, these reference values should be
  !        computed as the vertical averages from the surface down to
  !        epssfc*zgrid(kl).
  !
  !  This routine also computes where surface forcing is stable
  !  or unstable (STABLE)
  use timcom_kpp, only: bolus_sp, zgrid, lcheckekmo, lshort_wave, cmonob, &
                        vonkar, cekman, fcort, epssfc, vtc, ricr, concv, &
                        linertial, fstokes, llangmuir, hmxl
  implicit none
  ! !INPUT PARAMETERS:
  integer, intent(in) :: async_id
  real(r8), dimension(nx,ny,nz,2), intent(in) :: &
    trcr                ! tracers at current time

  real(r8), dimension(-1:nx+2,-1:ny+2,nz), intent(in) :: &
    uuu,vvv         ! velocities at current time

  real(r8), dimension(nx,ny,nz), intent(in) :: &
    dbloc,         &! buoyancy difference between adjacent levels
    dbsfc,         &! buoyancy difference between level and surface
    trans           ! transmission (1 to 0) from surface (calculated in kppglo)

  real(r8), dimension(nx,ny,2), intent(in) :: &
    stf             ! surface forcing for all tracers

  real(r8), dimension(nx,ny), intent(in) :: &
    shf_qsw         ! short-wave forcing

  real(r8), dimension(nx,ny,2), intent(in) :: &
    smfx            ! surface momentum forcing at U points or T points

  logical, intent(in) :: is_smft

  ! !OUTPUT PARAMETERS:
  integer, dimension(nx,ny), intent(out) :: &
    kbl                    ! index of first lvl below hbl

  real(r8), dimension(nx,ny), intent(out) :: &
    hblt,               &! boundary layer depth
    bfsfc,              &! Bo+radiation absorbed to d
    stable,             &! =1 stable forcing; =0 unstab
    ustar                ! surface friction velocity
  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer :: &
    i,j,k,kkk,               &! loop indices
    kupper, kup, kdn, ktmp, kl  ! vertical level indices

  real(r8), dimension(nx,ny) :: &
    vshear,            &! (velocity shear re sfc)^2
    sigma,             &! d/hbl
    wm, ws,            &! turb vel scale functions
    bo,                &! surface buoyancy forcing
    bosol,             &! radiative buoyancy forcing
    talpha,            &! temperature expansion coeff
    sbeta,             &! salinity    expansion coeff
    rho1,              &! density at the surface
    work,              &! temp array
    zkl,               &! depth at current z level
    b_frqncy,          &! buoyancy frequency
    rsh_hblt,          &! resolved shear contribution to HBLT (fraction)
    hlangm,            &! Langmuir depth
    hekman,            &! Eckman depth limit
    hlimit,test,       &! limit to mixed-layer depth (= min(HEKMAN,HMONOB))
    tmparr              ! temporary array

  real(r8), dimension(nx,ny,3) :: &
    ri_bulk,           &! Bulk Ri number at 3 lvls
    hmonob              ! Monin-Obukhov depth limit

  real(r8) ::          &
    absorb_frac,       &! shortwave absorption frac
    sqrt_arg,          &! dummy sqrt argument
    z_upper, z_up       ! upper depths for RI_BULK interpolation

  real(r8) :: &
    a_co, b_co, c_co    ! coefficients of the quadratic equation
                        ! $(a_{co}z^2+b_{co}|z|+c_{co}=Ri_b) used to
                        ! find the boundary layer depth. when
                        ! finding the roots, c_co = c_co - Ricr
  real(r8) :: &
    slope_up            ! slope of the above quadratic equation
                        ! at zup. this is used as a boundary
                        ! condition to determine the coefficients.
  real(r8) :: tmp, work1_0, work1_1, work1_2, work1_3
  !-----------------------------------------------------------------------
  !  compute friction velocity USTAR.
  !-----------------------------------------------------------------------

  !$acc enter data async(async_id) &
  !$acc& create(vshear, sigma, wm, ws, bo, bosol, talpha, sbeta, rho1, work, zkl, &
  !$acc&        b_frqncy, rsh_hblt, hlangm, hekman, hlimit, tmparr, ri_bulk, hmonob)

  !$acc parallel loop collapse(2) async(async_id) &
  !$acc& private(tmp, work1_0, work1_1, work1_2, work1_3, kupper, kup, kdn, ktmp, kl, &
  !$acc&         z_upper, z_up, slope_up, a_co, b_co, c_co, sqrt_arg)
  do j = 1, ny
    do i = 1, nx
      if(is_smft) then
        ustar(i,j) = dsqrt(dsqrt(smfx(i,j,1)**2 + smfx(i,j,2)**2))
      else
        work(i,j) = dsqrt(dsqrt(smfx(i,j,1)**2 + smfx(i,j,2)**2))
        ustar(i,j) = work(i,j)
      end if
      !-----------------------------------------------------------------------
      !  compute density and expansion coefficients at surface
      !-----------------------------------------------------------------------
      if (trcr(i,j,1,1) < -2.0d0) then
        work(i,j) = -2.0d0
      else
        work(i,j) = trcr(i,j,1,1)
      end if
      call state_gpu(i, j, 1, 1, work(i,j), trcr(i,j,1,2), &
        tmp, rho1(i,j), talpha(i,j), sbeta(i,j), .false., .true., .true., .true., zgrid(1), async_id)
      ! #ifdef chk_kpp1
      !   if(myid_x .eq. 7 .and. myid_y .eq. 11 .and. i .eq. nx/2 .and. j .eq. ny) then
      !     write(*,*) "[KPP CHECK bldepth]: 01"
      !     write(*,*) "ustar = ", ustar(nx/2,ny)
      !     write(*,*) "rho1 = ", rho1(nx/2,ny)
      !     write(*,*) "t2   = ", trcr(nx/2,ny,1,1)
      !     write(*,*) "s2   = ", trcr(nx/2,ny,1,2)
      !     write(*,*) achar(10)
      !   end if
      ! #endif
      !-----------------------------------------------------------------------
      !  compute turbulent and radiative sfc buoyancy forcing
      !-----------------------------------------------------------------------
      if(rho1(i,j) /= 0.d0) then
        bo(i,j) = grav*(-talpha(i,j)*stf(i,j,1) - &
                          sbeta(i,j)*stf(i,j,2))/rho1(i,j)

        bosol(i,j) = -grav*talpha(i,j)*shf_qsw(i,j)/rho1(i,j)
      else
        bo(i,j) = 0.d0
        bosol(i,j) = 0.d0
      end if
      !-----------------------------------------------------------------------
      !  Find bulk Richardson number at every grid level until > Ricr
      !  max values when Ricr never satisfied are KBL = KMT and
      !  HBLT = -zgrid(KMT)
      !
      !  NOTE: the reference depth is -epssfc/2.*zgrid(i,k), but the
      !        reference u,v,t,s values are simply the surface layer
      !        values and not the averaged values from 0 to 2*ref.depth,
      !        which is necessary for very fine grids(top layer < 2m
      !        thickness)
      !
      !
      !  Initialize hbl and kbl to bottomed out values
      !  Initialize HEKMAN and HLIMIT (= HMONOB until reset) to model bottom
      !  Initialize Monin Obukhov depth to value at z_up
      !  Set HMONOB=-zgrid(km) if unstable
      !-----------------------------------------------------------------------
      kupper = 1
      kup = 2
      kdn = 3
      z_upper = 0.d0
      z_up = zgrid(1)
      ri_bulk(i,j,kupper) = 0.d0
      ri_bulk(i,j,kup) = 0.d0
      if (int(kb(i,j)) > 1) then
        kbl(i,j) = int(kb(i,j))
      else
        kbl(i,j) = 1
      end if
      hlangm(i,j) = 0.d0

      !$acc loop seq
      do kl = 1, nz
        ! zkl = -zgrid(kl)
        if(kl == kbl(i,j)) hblt(i,j) = -zgrid(kl)
      end do

      ! #ifdef chk_kpp1
      !   if(myid_x .eq. 7 .and. myid_y .eq. 11) then
      !     write(*,*) "[KPP CHECK bldepth]: 02"
      !     write(*,*) "kbl  = ", kbl(nx/2,ny)
      !     write(*,*) "hblt = ", hblt(nx/2,ny)
      !     write(*,*) achar(10)
      !   end if
      ! #endif

      if(lcheckekmo) then
        hekman(i,j) = -zgrid(nz) + eps
        hlimit(i,j) = -zgrid(nz) + eps

        if(lshort_wave) then
          bfsfc(i,j) = bo(i,j) + bosol(i,j)*(1.d0-trans(i,j,1))
        else
          bfsfc(i,j) = bo(i,j)
        end if

        if (bfsfc(i,j) >= 0.d0) then
          stable(i,j) = 1.d0
        else
          stable(i,j) = 0.d0
        end if

        bfsfc(i,j) = bfsfc(i,j) + stable(i,j)*eps

        if (stable(i,j) > 0.5d0) then
          work(i,j) = cmonob*ustar(i,j)*ustar(i,j)*ustar(i,j)/vonkar/bfsfc(i,j)
        else
          work(i,j) = zgrid(nz)
        end if

        if (work(i,j) <= -z_up) then
          hmonob(i,j,kup) = -z_up+eps
        else
          hmonob(i,j,kup) = work(i,j)
        end if
      end if
      rsh_hblt(i,j) = 0.d0

      !-----------------------------------------------------------------------
      !  compute velocity shear squared on U-grid and use the maximum
      !  of the four surrounding U-grid values for the T-grid.
      !-----------------------------------------------------------------------
      !$acc loop seq
      do kl = 2, nz
        zkl(i,j) = -zgrid(kl)
        work1_0 = (uuu(i,j,1)-uuu(i,j,kl))**2 + &
                    (vvv(i,j,1)-vvv(i,j,kl))**2
        work1_1 = (uuu(i-1,j,1)-uuu(i-1,j,kl))**2 + &
                      (vvv(i-1,j,1)-vvv(i-1,j,kl))**2
        work1_2 = (uuu(i,j-1,1)-uuu(i,j-1,kl))**2 + &
                      (vvv(i,j-1,1)-vvv(i,j-1,kl))**2
        work1_3 = (uuu(i-1,j-1,1)-uuu(i-1,j-1,kl))**2 + &
                      (vvv(i-1,j-1,1)-vvv(i-1,j-1,kl))**2
        vshear(i,j) = max(work1_0, work1_1, work1_2, work1_3)
        !-----------------------------------------------------------------------
        !     compute bfsfc= Bo + radiative contribution down to hbf * hbl
        !     add epsilon to BFSFC to ensure never = 0
        !-----------------------------------------------------------------------
        if(lshort_wave) then
          bfsfc(i,j) = bo(i,j) + bosol(i,j)*(1.d0-trans(i,j,kl-1))
        else
          bfsfc(i,j) = bo(i,j)
        end if
        if (bfsfc(i,j) >= 0.d0) then
          stable(i,j) = 1.d0
        else
          stable(i,j) = 0.d0
        end if
        bfsfc(i,j) = bfsfc(i,j) + stable(i,j)*eps
        !-----------------------------------------------------------------------
        !     compute the Ekman and Monin Obukhov depths using above stability
        !-----------------------------------------------------------------------
        if(lcheckekmo) then
          if(stable(i,j) > 0.5d0 .and. hekman(i,j) >= -zgrid(nz) ) then
             hekman(i,j) = max(zkl(i,j), &
                         cekman*ustar(i,j)/(dabs(fcort(j))+eps))
          end if
          hmonob(i,j,kdn) = stable(i,j)*cmonob*ustar(i,j)*ustar(i,j)*ustar(i,j)/vonkar/bfsfc(i,j) + &
                            (stable(i,j)-1.0d0)*zgrid(nz)

          if(hmonob(i,j,kdn) <= zkl(i,j) .and. &
             hmonob(i,j,kup) >  -z_up) then
            work(i,j) = (hmonob(i,j,kdn) - hmonob(i,j,kup))/ &
                        (z_up + zkl(i,j))
            hlimit(i,j) = (hmonob(i,j,kdn) - work(i,j)*zkl(i,j))/ &
                          (1.0d0 - work(i,j))
          end if
        end if
        !-----------------------------------------------------------------------
        !     compute velocity scales at sigma, for hbl = -zgrid(kl)
        !-----------------------------------------------------------------------
        sigma(i,j) = epssfc
        call wscale_gpu(sigma(i,j), zkl(i,j), ustar(i,j), bfsfc(i,j), 2, wm(i,j), ws(i,j))
        !-----------------------------------------------------------------------
        !     compute the turbulent shear contribution to RI_BULK and store
        !     in WM.
        !-----------------------------------------------------------------------
        b_frqncy(i,j) = dsqrt( &
               0.5d0*(dbloc(i,j,kl) + dabs(dbloc(i,j,kl)) + eps2)/  &
               (zgrid(kl)-zgrid(kl+1)) )

        wm(i,j) = zkl(i,j)*ws(i,j)*b_frqncy(i,j)* &
                 ((vtc/ricr)*max(2.1d0 - 200.d0*b_frqncy(i,j),concv) )
        !-----------------------------------------------------------------------
        !     compute bulk Richardson number at new level
        !-----------------------------------------------------------------------
        if (kb(i,j) >= kl) then
          work(i,j) = (zgrid(1)-zgrid(kl))*dbsfc(i,j,kl)
        else
          work(i,j) = 0.d0
        end if

        if(linertial) then
          ri_bulk(i,j,kdn) = work(i,j)/(vshear(i,j)+wm(i,j)+ustar(i,j)*bolus_sp(i,j)+eps)
        else
          ri_bulk(i,j,kdn) = work(i,j)/(vshear(i,j)+wm(i,j)+eps)
        end if
        !-----------------------------------------------------------------------
        !       find hbl where Rib = Ricr. if possible, use a quadratic
        !       interpolation. if not, linearly interpolate. the quadratic
        !       equation coefficients are determined using the slope and
        !       Ri_bulk at z_up and Ri_bulk at zgrid(kl). the slope at
        !       z_up is computed linearly between z_upper and z_up.
        !       compute Langmuir depth always
        !-----------------------------------------------------------------------
        if(kbl(i,j) == kb(i,j) .and. ri_bulk(i,j,kdn) > ricr) then
          slope_up = (ri_bulk(i,j,kupper) - ri_bulk(i,j,kup))/(z_up - z_upper)

          a_co = (ri_bulk(i,j,kdn) - ri_bulk(i,j,kup) -         &
                  slope_up*(zkl(i,j) + z_up) )/(z_up + zkl(i,j))**2

          b_co = slope_up + 2.d0 * a_co * z_up

          c_co = ri_bulk(i,j,kup) + z_up*(a_co*z_up + slope_up) - ricr

          sqrt_arg = b_co**2 - 4.d0*a_co*c_co

          if((dabs(b_co) > eps .and. dabs(a_co)/dabs(b_co) <= eps ) &
              .or. sqrt_arg <= 0.d0 ) then

            hblt(i,j) = -z_up + (z_up + zkl(i,j)) *               &
                        (ricr             - ri_bulk(i,j,kup))/    &
                        (ri_bulk(i,j,kdn) - ri_bulk(i,j,kup))

          else
            hblt(i,j) = (-b_co + dsqrt(sqrt_arg)) / (2.d0*a_co)
          end if
          kbl(i,j) = kl
          rsh_hblt(i,j) =  (vshear(i,j)*ricr/ &
                           (dbsfc(i,j,kl)+eps))/hblt(i,j)

          hlangm(i,j) = ustar(i,j) * dsqrt(fstokes(i,j)*zkl(i,j)/(dbsfc(i,j,kl)+eps))/0.9d0
        end if
        ! #ifdef chk_kpp1
        !   if(myid_x .eq. 7 .and. myid_y .eq. 11) then
        !     write(*,*) "[KPP CHECK bldepth]: 03"
        !     write(*,*) "level k =", kl
        !     write(*,*) "kbl  = ", kbl(nx/2,ny)
        !     write(*,*) "hblt = ", hblt(nx/2,ny)
        !     write(*,*) achar(10)
        !   end if
        ! #endif
        !-----------------------------------------------------------------------
        !     swap klevel indices and move to next level
        !-----------------------------------------------------------------------
        ktmp   = kupper
        kupper = kup
        kup    = kdn
        kdn    = ktmp
        z_upper = z_up
        z_up    = zgrid(kl)
      end do

      !$acc loop seq
      do kl = nz,2,-1
        !-----------------------------------------------------------------------
        !     apply Langmuir parameterization if requested
        !-----------------------------------------------------------------------
        if(llangmuir) then
          if (hlangm(i,j) > hblt(i,j)          .and. &
              hlangm(i,j) >  -zgrid(kl-1) .and. &
              hlangm(i,j) <= zkl(i,j)) then
            hblt(i,j) = hlangm(i,j)
            kbl(i,j)  = kl
          end if
        end if
      end do
      !-----------------------------------------------------------------------
      !  first combine Ekman and Monin-Obukhov depth limits. then apply
      !  these restrictions to HBLT. note that HLIMIT is set to -zgrid(km)
      !  in unstable forcing.
      !-----------------------------------------------------------------------
      if(lcheckekmo) then
        if (hekman(i,j) < hlimit(i,j)) then
          hlimit(i,j) = hekman(i,j)
        end if
        do kl = 2, nz
          if (hlimit(i,j) < hblt(i,j)         .and. &
              hlimit(i,j) >  -zgrid(kl-1) .and. &
              hlimit(i,j) <= zkl(i,j)) then
            hblt(i,j) = hlimit(i,j)
            kbl(i,j) = kl
          end if
        end do
      end if
    end do
  end do

  !-----------------------------------------------------------------------
  !  apply a Gaussian filter
  !-----------------------------------------------------------------------
  ! #ifdef chk_kpp1
  !   if(myid_x .eq. 7 .and. myid_y .eq. 11) then
  !     write(*,*) "[KPP CHECK bldepth]: 04"
  !     write(*,*) "kbl  = ", kbl(nx/2,ny)
  !     write(*,*) "hblt = ", hblt(nx/2,ny)
  !     write(*,*) achar(10)
  !   end if
  ! #endif

  call smooth_hblt_gpu(.true., .false., hblt, kbl, tmparr, async_id)

  ! #ifdef chk_kpp1
  !   if(myid_x .eq. 7 .and. myid_y .eq. 11) then
  !     write(*,*) "[KPP CHECK bldepth]: 05"
  !     write(*,*) "kbl  = ", kbl(nx/2,ny)
  !     write(*,*) "hblt = ", hblt(nx/2,ny)
  !     write(*,*) achar(10)
  !   end if
  ! #endif
  !-----------------------------------------------------------------------
  !  correct stability and buoyancy forcing for SW up to boundary layer
  !-----------------------------------------------------------------------
  !$acc parallel loop collapse(2) async(async_id)
  do j = 1, ny
    do i = 1, nx
      if(lshort_wave) then
        bfsfc(i,j) = bo(i,j) + bosol(i,j)*(1.0d0-trans(i,j,1))
      end if
      if (bfsfc(i,j) >= 0.d0) then
        stable(i,j) = 1.0d0
      else
        stable(i,j) = 0.d0
      end if
      bfsfc(i,j) = bfsfc(i,j) + stable(i,j) * eps ! ensures bfsfc never=0
    end do
  end do

  !$acc exit data async(async_id) &
  !$acc& delete(vshear, sigma, wm, ws, bo, bosol, talpha, sbeta, rho1, work, zkl, &
  !$acc&        b_frqncy, rsh_hblt, hlangm, hekman, hlimit, tmparr, ri_bulk, hmonob)

end subroutine bldepth_gpu

subroutine wscale_gpu(sigma, hbl, ustar, bfsfc, m_or_s, wm, ws)
  !$acc routine seq
  ! !DESCRIPTION:
  !  Computes turbulent velocity scales.
  !
  !  For $\zeta \geq 0,
  !    w_m = w_s = \kappa U^\star/(1+5\zeta)$
  !
  !  For $\zeta_m \leq \zeta < 0,
  !    w_m = \kappa U^\star (1-16\zeta)^{1\over 4}$
  !
  !  For $\zeta_s \leq \zeta < 0,
  !    w_s = \kappa U^\star (1-16\zeta)^{1\over 2}$
  !
  !  For $\zeta < \zeta_m,
  !    w_m = \kappa U^\star (a_m - c_m\zeta)^{1\over 3}$
  !
  !  For $\zeta < \zeta_s,
  !    w_s = \kappa U^\star (a_s - c_s\zeta)^{1\over 3}$
  !
  !  where $\kappa$ is the von Karman constant.
  use timcom_kpp, only: vonkar, zeta_m, zeta_s, a_m, c_m, a_s, c_s
  implicit none
  ! !INPUT PARAMETERS:
  integer, intent(in) :: &
    m_or_s              ! flag =1 for wm only, 2 for ws, 3 for both

  real(r8), intent(in) :: &
    sigma,             &! normalized depth (d/hbl)
    hbl,               &! boundary layer depth
    bfsfc,             &! surface buoyancy forcing
    ustar               ! surface friction velocity

  ! !OUTPUT PARAMETERS:
  real(r8), intent(out) :: &
    wm,                &! turb velocity scales: momentum
    ws                  ! turb velocity scales: tracer

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  real(r8) :: &
    zeta,           &! d/L or sigma*hbl/L(monin-obk)
    zetah            ! sigma*hbl*vonkar*BFSFC or ZETA = ZETAH/USTAR**3

  !-----------------------------------------------------------------------
  !  compute zetah and zeta - surface layer is special case
  !-----------------------------------------------------------------------
  zetah = sigma*hbl*vonkar*bfsfc
  zeta = zetah/(ustar**3 + eps)
  !-----------------------------------------------------------------------
  !  compute velocity scales for momentum
  !-----------------------------------------------------------------------
  if(m_or_s == 1 .or. m_or_s == 3) then
    if(zeta >= 0.d0) then ! stable region
      wm = vonkar*ustar/(1.0d0 + 5.0d0*zeta)
    else if (zeta >= zeta_m) then
      wm = vonkar*ustar*(1.0d0 - 16.0d0*zeta)**0.25d0
    else
      wm = vonkar*(a_m*(ustar**3)-c_m*zetah)**0.33d0
    end if
  end if
  !-----------------------------------------------------------------------
  !  compute velocity scales for tracers
  !-----------------------------------------------------------------------
  if(m_or_s == 2 .or. m_or_s == 3) then
    if(zeta >= 0.d0) then
      ws = vonkar*ustar/(1.0d0 + 5.0d0*zeta)
    else if (zeta >= zeta_s) then
      ws = vonkar*ustar*dsqrt(1.0d0 - 16.0d0*zeta)
    else
      ws = vonkar*(a_s*(ustar**3)-c_s*zetah)**0.33d0
    end if
  end if

end subroutine wscale_gpu

subroutine smooth_hblt_gpu(overwrite_hblt, use_hmxl, &
                       hblt, kbl, smooth_out, async_id)
  ! !DESCRIPTION:
  !  This subroutine uses a 1-1-4-1-1 Laplacian filter one time
  !  on HBLT or HMXL to reduce any horizontal two-grid-point noise.
  !  If HBLT is overwritten, KBL is adjusted after smoothing.
  use timcom_kpp, only: hmxl, zgrid
  use hyperlink, only: m_comm_cart, r8type2d, nbid, ndim, symm_np
  implicit none
  ! !INPUT PARAMETERS:
  integer, intent(in) :: async_id
  logical, intent(in) :: &
    overwrite_hblt,   &    ! if .true.,  HBLT is overwritten
                             ! if .false., the result is returned in
                             !  a dummy array
    use_hmxl               ! if .true., smooth HMXL
                             ! if .false., smooth HBLT
  ! !INPUT/OUTPUT PARAMETERS:

  real(r8), dimension(nx,ny), intent(inout) :: &
    hblt                   ! boundary layer depth

  integer, dimension(nx,ny), intent(inout) :: &
    kbl                    ! index of first lvl below hbl

  ! !OUTPUT PARAMETERS:

  real(r8), dimension(nx,ny), intent(out) ::  &
    smooth_out              ! optional output array containing the
                              !  smoothened field if overwrite_hblt is false
  !-----------------------------------------------------------------------
  !     local variables
  !-----------------------------------------------------------------------
  integer :: &
    i, j, k

  real(r8) ::  &
    work1(0:nx+1,0:ny+1), work2(nx,ny)

  real(r8) ::  &
    ccc, ccw, cce, ccn, ccs, &  ! averaging weights, * change to ccc and ccs
    ztmp                   ! temp for level depth

  !$acc enter data create(work1, work2) async(async_id)

  !-----------------------------------------------------------------------
  !     perform one smoothing pass since we cannot do the necessary
  !     boundary updates for multiple passes.
  !-----------------------------------------------------------------------

  !$acc kernels async(async_id)
  if(use_hmxl) then
    work2 = hmxl
  else
    work2 = hblt
  endif
  work1 = 0.d0
  work1(1:nx,1:ny) = work2
  !$acc end kernels

  call ghost_cell_exch_gpu(m_comm_cart, r8type2d, nbid, ndim, work1, symm_np, async_id)

  !$acc parallel loop collapse(2) private(ccw, cce, ccn, ccs, ccc) async(async_id)
  do j = 1, ny
    do i = 1, nx
      if(in(i,j,1)*kb(i,j) /= 0 ) then
        ccw = .125d0
        cce = .125d0
        ccn = .125d0
        ccs = .125d0
        ccc = .5d0
        if(in(i-1,j,1)*kb(i-1,j) == 0) then
          ccc = ccc + ccw
          ccw = 0.d0
        end if
        if(in(i+1,j,1)*kb(i+1,j) == 0) then
          ccc = ccc + cce
          cce = 0.d0
        end if
        if(in(i,j-1,1)*kb(i,j-1) == 0 ) then
          ccc = ccc + ccs
          ccs = 0.d0
        end if
        if(in(i,j+1,1)*kb(i,j+1) == 0 ) then
          ccc = ccc + ccn
          ccn = 0.d0
        end if
        work2(i,j) =  ccw * work1(i-1,j)   &
                    + cce * work1(i+1,j)   &
                    + ccs * work1(i,j-1)   &
                    + ccn * work1(i,j+1)   &
                    + ccc * work1(i,j)
      end if
    end do
  end do

  !$acc parallel loop collapse(2) private(ztmp) async(async_id)
  do j = 1, ny
    do i = 1, nx
      !$acc loop vector
      do k = 1, nz
        ztmp = -zgrid(k)
        if(k == kb(i,j) .and. work2(i,j) > ztmp ) then
          work2(i,j) = ztmp
        end if
      end do
      !$acc loop seq
      do k = 1, nz
        if(overwrite_hblt .and. .not.use_hmxl) then
          hblt(i,j) = work2(i,j)
          ztmp = -zgrid(k)
          if(in(i,j,1)*kb(i,j) /= 0  .and. &
            (hblt(i,j)> -zgrid(k-1)) .and. &
            (hblt(i,j)<= ztmp      )) kbl(i,j) = k
        else
          smooth_out(i,j) = work2(i,j)
        end if
      end do
    end do
  end do

  !$acc exit data delete(work1, work2) async(async_id)

end subroutine smooth_hblt_gpu

subroutine blmix_gpu(visc, vdc, hblt, ustar, bfsfc, stable, &
                  kbl, ghat, async_id)
  !$acc routine(wscale_gpu) seq
  ! !DESCRIPTION:
  !  This routine computes mixing coefficients within boundary layer
  !  which depend on surface forcing and the magnitude and gradient
  !  of interior mixing below the boundary layer (matching).  These
  !  quantities have been computed in other routines.
  !
  !  Caution: if mixing bottoms out at hbl = -zgrid(km) then
  !  fictitious layer at km+1 is needed with small but finite width
  !  hwide(km+1).
  use timcom_kpp, only: epssfc, zgrid, hwide, cg
  implicit none
  ! !INPUT/OUTPUT PARAMETERS:
  integer, intent(in) :: async_id
  real(r8), dimension(nx,ny,0:nzf), intent(inout) :: &
    visc               ! interior mixing coeff on input
                         ! combined interior/bndy layer coeff output

  real(r8), dimension(nx,ny,0:nzf,2), intent(inout) :: &
    vdc        ! diffusivity for tracer diffusion

  ! !INPUT PARAMETERS:

  integer, dimension(nx,ny), intent(in) ::  &
    kbl                    ! index of first lvl below hbl

  real(r8), dimension(nx,ny), intent(in) ::  &
    hblt,                 & ! boundary layer depth
    bfsfc,                & ! surface buoyancy forcing
    stable,               & ! =1 stable forcing; =0 unstab
    ustar                   ! surface friction velocity
  ! !OUTPUT PARAMETERS:

  real(r8), dimension(nx,ny,nz),intent(out) :: &
    ghat                ! non-local mixing coefficient
  !-----------------------------------------------------------------------
  !
  !  local variables
  !
  !-----------------------------------------------------------------------
  integer :: &
    k,kp1,             &! dummy k level index
    i,j,               &! horizontal indices
    n

  integer, dimension(nx,ny) :: &
    kn                  ! klvl closest to HBLT

  real(r8), dimension(nx,ny,nz,3) :: &
    blmc                ! bndy layer mixing coefs

  real(r8), dimension(nx,ny,3) :: &
    gat1,              &! shape function at sigma=1
    dat1,              &! derivative of shape function
    dkm1                ! bndy layer difs at kbl-1 lvl

  real(r8), dimension(nx,ny) :: &
    wm,ws,             &! turbulent velocity scales
    casea,             &! =1 in case A, =0 in case B
    sigma,             &! normalized depth (d/hbl)
    visch,             &! viscosity at hbl
    difth,             &! temp diffusivity at hbl
    difsh,             &! tracer diffusivity at hbl
    delhat, r, dvdzup, dvdzdn, &
    viscp, diftp, difsp, f1, &
    work1,work2

  !$acc enter data async(async_id) &
  !$acc& create(kn, blmc, gat1, dat1, dkm1, wm, ws, casea, sigma, visch, difth, difsh, &
  !$acc&        delhat, r, dvdzup, dvdzdn, viscp, diftp, difsp, f1, work1, work2)

  !$acc parallel loop collapse(2) private(k) async(async_id)
  do j = 1, ny
    do i = 1, nx
      !-----------------------------------------------------------------------
      !  compute velocity scales at hbl
      !-----------------------------------------------------------------------
      sigma(i,j) = epssfc
      call wscale_gpu(sigma(i,j), hblt(i,j), ustar(i,j), bfsfc(i,j), 3, wm(i,j), ws(i,j))
      !-----------------------------------------------------------------------
      !  determine caseA = 0 if closer to KBL than KBL-1
      !  KN is then the closest klevel to HBLT
      !----------------------------------------------------------------------
      k = kbl(i,j)
      casea(i,j)  = .5d0 + dsign(.5d0, -zgrid(k)-.5d0*hwide(k)-hblt(i,j))
      kn(i,j) = nint(casea(i,j))*(kbl(i,j)-1) + (1-nint(casea(i,j)))*kbl(i,j)
      !-------------------------------------------------------------------
      !  find the interior viscosities and derivatives at hbl by
      !  interpolating derivative values at vertical interfaces.  compute
      !  matching conditions for shape function.
      !-----------------------------------------------------------------------
      f1(i,j) = stable(i,j)*5.d0*bfsfc(i,j)/(ustar(i,j)**4+eps)
    end do
  end do

  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if(k == kn(i,j)) then
          delhat(i,j) = .5d0*hwide(k) - zgrid(k) - hblt(i,j)
          r     (i,j) = 1.d0 - delhat(i,j) / hwide(k)

          dvdzup(i,j) = (visc(i,j,k-1) - visc(i,j,k  ))/hwide(k)
          dvdzdn(i,j) = (visc(i,j,k  ) - visc(i,j,k+1))/hwide(k+1)
          viscp (i,j) = .5d0*( (1.d0-r(i,j))* &
                               (dvdzup(i,j) + dabs(dvdzup(i,j))) + &
                                   r(i,j) * &
                               (dvdzdn(i,j) + dabs(dvdzdn(i,j))) )

          dvdzup(i,j) = (vdc(i,j,k-1,2) - vdc(i,j,k  ,2))/hwide(k)
          dvdzdn(i,j) = (vdc(i,j,k  ,2) - vdc(i,j,k+1,2))/hwide(k+1)
          difsp (i,j) = .5d0*( (1.d0-r(i,j))* &
                               (dvdzup(i,j) + dabs(dvdzup(i,j))) + &
                                   r(i,j) * &
                               (dvdzdn(i,j) + dabs(dvdzdn(i,j))) )
          dvdzup(i,j) = (vdc(i,j,k-1,1) - vdc(i,j,k  ,1))/hwide(k)
          dvdzdn(i,j) = (vdc(i,j,k  ,1) - vdc(i,j,k+1,1))/hwide(k+1)
          diftp (i,j) = .5d0*( (1.d0-r(i,j))* &
                               (dvdzup(i,j) + dabs(dvdzup(i,j))) + &
                                   r(i,j) * &
                               (dvdzdn(i,j) + dabs(dvdzdn(i,j))) )
          visch(i,j) = visc(i,j,k)  + viscp(i,j)*delhat(i,j)
          difsh(i,j) = vdc(i,j,k,2) + difsp(i,j)*delhat(i,j)
          difth(i,j) = vdc(i,j,k,1) + diftp(i,j)*delhat(i,j)
          gat1(i,j,1) = visch(i,j) / hblt(i,j) /(wm(i,j)+eps)
          dat1(i,j,1) = -viscp(i,j)/(wm(i,j)+eps) + &
                           f1(i,j)*visch(i,j)

          gat1(i,j,2) = difsh(i,j) / hblt(i,j) /(ws(i,j)+eps)
          dat1(i,j,2) = -difsp(i,j)/(ws(i,j)+eps) + &
                           f1(i,j)*difsh(i,j)

          gat1(i,j,3) = difth(i,j) / hblt(i,j) /(ws(i,j)+eps)
          dat1(i,j,3) = -diftp(i,j)/(ws(i,j)+eps) + &
                           f1(i,j)*difth(i,j)
        end if
      end do
    end do
  end do

  !$acc parallel loop collapse(3) async(async_id)
  do j = 1, ny
    do i = 1, nx
      do n = 1, 3
        dat1(i,j,n) = min(dat1(i,j,n),0.d0)
      end do
    end do
  end do

  !-----------------------------------------------------------------------
  !  compute the dimensionless shape functions and diffusivities
  !  at the grid interfaces.  also compute function for non-local
  !  transport term (GHAT).
  !-----------------------------------------------------------------------

  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        sigma(i,j) = (-zgrid(k) + 0.5d0*hwide(k))/hblt(i,j)
        f1(i,j) = min(sigma(i,j),epssfc)

        call wscale_gpu(f1(i,j), hblt(i,j), ustar(i,j), bfsfc(i,j), 3, wm(i,j), ws(i,j))

        blmc(i,j,k,1) = hblt(i,j)*wm(i,j)*sigma(i,j)*       &
                         (1.d0 + sigma(i,j)*((sigma(i,j)-2.d0) + &
                         (3.d0-2.d0*sigma(i,j))*gat1(i,j,1) +    &
                         (sigma(i,j)-1.d0)*dat1(i,j,1)))

        blmc(i,j,k,2) = hblt(i,j)*ws(i,j)*sigma(i,j)*       &
                         (1.d0 + sigma(i,j)*((sigma(i,j)-2.d0) + &
                         (3.d0-2.d0*sigma(i,j))*gat1(i,j,2) +    &
                         (sigma(i,j)-1.d0)*dat1(i,j,2)))

        blmc(i,j,k,3) = hblt(i,j)*ws(i,j)*sigma(i,j)*       &
                         (1.d0 + sigma(i,j)*((sigma(i,j)-2.d0) + &
                         (3.d0-2.d0*sigma(i,j))*gat1(i,j,3) +    &
                         (sigma(i,j)-1.d0)*dat1(i,j,3)))
        ghat(i,j,k) = (1.d0-stable(i,j))* cg/(ws(i,j)*hblt(i,j) +eps)
      end do
    end do
  end do

  !-----------------------------------------------------------------------
  !  find diffusivities at kbl-1 grid level
  !-----------------------------------------------------------------------

  !$acc parallel loop collapse(2) async(async_id)
  do j = 1, ny
    do i = 1, nx
      k = kbl(i,j) - 1
      sigma(i,j) = -zgrid(k)/hblt(i,j)
      f1(i,j) = min(sigma(i,j),epssfc)
      call wscale_gpu(f1(i,j), hblt(i,j), ustar(i,j), bfsfc(i,j), 3, wm(i,j), ws(i,j))

      dkm1(i,j,1) = hblt(i,j)*wm(i,j)*sigma(i,j)*     &
                    (1.d0+sigma(i,j)*((sigma(i,j)-2.d0) + &
                    (3.d0-2.d0*sigma(i,j))*gat1(i,j,1)  + &
                    (sigma(i,j)-1.d0)*dat1(i,j,1)))

      dkm1(i,j,2) = hblt(i,j)*ws(i,j)*sigma(i,j)*     &
                    (1.d0+sigma(i,j)*((sigma(i,j)-2.d0) + &
                    (3.d0-2.d0*sigma(i,j))*gat1(i,j,2)  + &
                    (sigma(i,j)-1.d0)*dat1(i,j,2)))

      dkm1(i,j,3) = hblt(i,j)*ws(i,j)*sigma(i,j)*     &
                    (1.d0+sigma(i,j)*((sigma(i,j)-2.d0) + &
                    (3.d0-2.d0*sigma(i,j))*gat1(i,j,3)  + &
                    (sigma(i,j)-1.d0)*dat1(i,j,3)))

    end do
  end do

  !$acc parallel loop collapse(3) async(async_id)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        !-----------------------------------------------------------------------
        !  compute the enhanced mixing
        !-----------------------------------------------------------------------
        if(k < nz .and. k == (kbl(i,j) - 1)) then
          delhat(i,j) = (hblt(i,j) + zgrid(k))/(zgrid(k)-zgrid(k+1))

          blmc(i,j,k,1) = (1.d0-delhat(i,j))*visc(i,j,k) +           &
                        delhat(i,j) *(                                 &
                        (1.d0-delhat(i,j))**2*dkm1(i,j,1) +            &
                        delhat(i,j)**2*(casea(i,j)*visc(i,j,k)+        &
                        (1.d0-casea(i,j))*blmc(i,j,k,1)))

          blmc(i,j,k,2) = (1.d0-delhat(i,j))*vdc(i,j,k,2) +          &
                        delhat(i,j)*(                                  &
                        (1.d0-delhat(i,j))**2*dkm1(i,j,2) +            &
                        delhat(i,j)**2*(casea(i,j)*vdc(i,j,k,2)+       &
                        (1.d0-casea(i,j))*blmc(i,j,k,2)))

          blmc(i,j,k,3) = (1.d0-delhat(i,j))*vdc(i,j,k,1) +          &
                        delhat(i,j) *(                                 &
                        (1.d0-delhat(i,j))**2*dkm1(i,j,3) +            &
                        delhat(i,j)**2*(casea(i,j)*vdc(i,j,k,1)+       &
                        (1.d0-casea(i,j))*blmc(i,j,k,3)))

          ghat(i,j,k) = (1.d0-casea(i,j)) * ghat(i,j,k)

        end if
        !-----------------------------------------------------------------------
        !  combine interior and boundary layer coefficients and nonlocal term
        !-----------------------------------------------------------------------
        if(k < kbl(i,j)) then
          visc(i,j,k)  = blmc(i,j,k,1)
          vdc(i,j,k,2) = blmc(i,j,k,2)
          vdc(i,j,k,1) = blmc(i,j,k,3)
        else
          ghat(i,j,k) = 0.d0
        end if
      end do
    end do
  end do

  !$acc exit data async(async_id) &
  !$acc& delete(kn, blmc, gat1, dat1, dkm1, wm, ws, casea, sigma, visch, difth, difsh, &
  !$acc&        delhat, r, dvdzup, dvdzdn, viscp, diftp, difsp, f1, work1, work2)

end subroutine blmix_gpu

end module timcom_kpp_gpu
