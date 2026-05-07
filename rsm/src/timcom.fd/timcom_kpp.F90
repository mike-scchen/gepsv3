module tai_timcom_kpp
  use tai_timcom_const
  use tai_timcom_comm
  use tai_hyperlink, only: nx, ny, nz, nzf, kb, in, &
                       dt, iw, dz, odz, dzw, odzw, &
                       z_grid, z_face, x_grid, y_grid,&
                       myid, myid_x, myid_y

  implicit none

  integer, parameter :: num_v_smooth_Ri = 1
  logical, parameter :: &
    ldbl_diff     = .true., & ! dbl diffusion
    lshort_wave   = .true., & ! computing short-wave forcing
    ltidal_mixing = .true., & ! tidal mixing
    lccsm_control_compatible = .false., &
    lhoriz_varying_bckgrnd = .true.,    & ! horizontal varying backgrd vdc,vvc
    lcheckekmo = .false., & !check Ekman, Monin-Obhukov depth limit
    llangmuir = .false., &  !! Langmuir parameterization
    linertial = .false., & ! inertial mixing parameterization
    lrich = .true., &   ! compute Ri-dependent mixing
    lkpp_ini = .true.

  real(r8) :: &
    cg,       & !coefficient for counter-gradient term
    Vtc         !resolution and buoyancy independent part of the
                !turbulent velocity shear coefficient (for bulk Ri no)
  real(r8), parameter :: &
    convect_diff =  20000.d0,     &! diffusivity to mimic convection
    convect_visc =  20000.d0, &       ! viscosity   to mimic convection
    hmix_tracer_type_gm = 0.d0, &
    Prandtl = 10.d0, &
    rich_mix = 50.d0,      &! coefficient for rich number term
    Riinfty = .8d0,        &! Rich. no. limit for shear instability
    BVSQcon = .0d0, &          ! Brunt-Vaisala square cutoff(s**-2)
    epssfc = .1d0, &
    vonkar = .4d0,         &! von Karman constant
    cstar = 10.d0, &       ! coeff for nonlocal transport
    Ricr = 0.3d0, &
    cmonob = 1.d0,         &! coefficient for Monin-Obukhov depth
    cekman = .7d0,         &! coefficient for Ekman depth
    concv = 1.7d0, &
    Rrho0  = 2.55d0,     &! limit for double-diff density ratio
    dsfmax = 1.d0, &        ! max diffusivity for salt fingering
    !** parameters for velocity scale function (from Large et al.)
    zeta_m = -0.2d0,      &
    zeta_s = -1.d0,      &
    c_m    =  8.38d0,     &
    c_s    =  98.96d0,    &
    a_m    =  1.26d0,     &
    a_s    = -28.86d0,    &
    tidal_mix_max = 100.d0
  real(r8), pointer, dimension(:) :: &
    zgrid, & !(0:nzf)
    hwide, & !(0:nzf)
    fcort    !(ny)
  real(r8), pointer, dimension(:,:,:) :: &
    bckgrnd_vvc, & !(nx,ny,nz)
    bckgrnd_vdc, & !(nx,ny,nz)
    tidal_diff,  & !(nx,ny,nz)
    tidal_coef     !(nx,ny,nz)
  real(r8), pointer, dimension(:,:) :: &
    fstokes, &   !ratio of stokes velocity to ustar
    hmxl, &
    bolus_sp

  type :: kpp_panel
    real(r8), pointer, dimension(:) :: &
      zgrid, & !(0:nzf)
      hwide, & !(0:nzf)
      fcort    !(ny)
    real(r8), pointer, dimension(:,:,:) :: &
      bckgrnd_vvc, & !(nx,ny,nz)
      bckgrnd_vdc, & !(nx,ny,nz)
      tidal_diff,  & !(nx,ny,nz)
      tidal_coef     !(nx,ny,nz)
    real(r8), pointer, dimension(:,:) :: &
      fstokes, &   !ratio of stokes velocity to ustar
      hmxl, &
      bolus_sp
  end type kpp_panel

  type(kpp_panel) :: kpp_var
contains

subroutine init_vmix_kpp(tidal_energy_flux)
  implicit none
!-----------------------------------------------------------------------
! local coefficients
!-----------------------------------------------------------------------
! type(kpp_panel) :: kpp_var
  real(r8) :: tidal_energy_flux(nx,ny)
  integer :: i,j,k
  real(r8) :: bckgrnd_vdc_psis  ! PSI diffusivity in northern hemisphere
  real(r8) :: bckgrnd_vdc_psin  ! PSI diffusivity in southern hemisphere
  real(r8), dimension(nz) :: zw ! vert dist from sfc to bottom of layer
  real(r8), parameter :: bckgrnd_vdc1 = 0.16d0 ! background diffusivity (Ledwell)
  real(r8), parameter :: bckgrnd_vdc2 = 0.d0   ! variation in diffusivity
  real(r8), parameter :: bckgrnd_vdc_eq   = 0.01d0    ! equatorial diffusivity (Gregg)
  real(r8), parameter :: bckgrnd_vdc_psim = 0.13d0    ! Max. PSI induced diffusivity (MacKinnon)
  real(r8), parameter :: bckgrnd_vdc_ban  = 1.0d0     ! Banda Sea diffusivity (Gordon)
  real(r8), parameter :: bckgrnd_vdc_dpth = 1000.0d02 !2500 ! depth at which diff equals vdc1
  real(r8), parameter :: bckgrnd_vdc_linv = 4.5d-05   ! inverse length for transition region
! tidal mixing
  real(r8), dimension(nx,ny) :: &
    vertical_func,              &
    work1,HT
  real(r8) :: coef, rho_fw
  real(r8), parameter ::             &
    local_mixing_fraction = 0.33d0,  &
    mixing_efficiency = 0.2d0,       &
    vertical_decay_scale = 500.0d02, &
    tidal_mix_max = 100.d0

  allocate(kpp_var%zgrid(0:nzf))
  zgrid => kpp_var%zgrid
  zgrid = 0.d0
  
  allocate(kpp_var%hwide(0:nzf))
  hwide => kpp_var%hwide
  hwide = 0.d0

  allocate(kpp_var%fcort(ny))
  fcort => kpp_var%fcort
  fcort = 0.d0

  allocate(kpp_var%bckgrnd_vvc(nx,ny,nz))
  bckgrnd_vvc => kpp_var%bckgrnd_vvc
  bckgrnd_vvc = 0.d0

  allocate(kpp_var%bckgrnd_vdc(nx,ny,nz))
  bckgrnd_vdc => kpp_var%bckgrnd_vdc
  bckgrnd_vdc = 0.d0

  allocate(kpp_var%tidal_diff(nx,ny,nz))
  tidal_diff => kpp_var%tidal_diff
  tidal_diff = 0.d0

  allocate(kpp_var%tidal_coef(nx,ny,nz))
  tidal_coef => kpp_var%tidal_coef
  tidal_coef = 0.d0

  allocate(kpp_var%fstokes(nx,ny))
  fstokes => kpp_var%fstokes
  fstokes = 0.d0

  allocate(kpp_var%hmxl(nx,ny))
  hmxl => kpp_var%hmxl
  hmxl = 0.d0

  allocate(kpp_var%bolus_sp(nx,ny))
  bolus_sp => kpp_var%bolus_sp
  bolus_sp = 0.d0

  if(lhoriz_varying_bckgrnd) then
    k = 1
    do j = 1, ny
      do i =1, nx
        bckgrnd_vdc_psis   = bckgrnd_vdc_psim*exp(-(0.4d0*(y_grid(j)+28.9d0))**2.)
        bckgrnd_vdc_psin   = bckgrnd_vdc_psim*exp(-(0.4d0*(y_grid(j)-28.9d0))**2.)
        bckgrnd_vdc(i,j,k) = bckgrnd_vdc_eq+bckgrnd_vdc_psin+bckgrnd_vdc_psis
     
        if ( y_grid(j) .lt. -10.d0 ) then
          bckgrnd_vdc(i,j,k) = bckgrnd_vdc(i,j,k) + bckgrnd_vdc1
        else if ( y_grid(j) .le. 10.d0 ) then
          bckgrnd_vdc(i,j,k) = bckgrnd_vdc(i,j,k) + bckgrnd_vdc1*(y_grid(j)/10.d0)**2.
        else
          bckgrnd_vdc(i,j,k) = bckgrnd_vdc(i,j,k) + bckgrnd_vdc1
        end if
      !----------------
      ! North Banda Sea
      !----------------
        if((y_grid(j) .lt. -1.0d0)  .and. (y_grid(j) .gt. -4.0d0)  .and.  &
           (x_grid(i) .gt. 103.0d0) .and. (x_grid(i) .lt. 134.0d0)) then
          bckgrnd_vdc(i,j,k) = bckgrnd_vdc_ban
        end if
      !-----------------
      ! Middle Banda Sea
      !-----------------
        if((y_grid(j) .le. -4.0d0)  .and. (y_grid(j) .gt. -7.0d0)  .and.  &
           (x_grid(i) .gt. 106.0d0) .and. (x_grid(i) .lt. 140.d0)) then
          bckgrnd_vdc(i,j,k) = bckgrnd_vdc_ban
        end if
      !----------------
      ! South Banda Sea
      !----------------
        if((y_grid(j) .le. -7.0d0)  .and. (y_grid(j) .gt. -8.3d0)  .and.  &
           (x_grid(i) .gt. 111.0d0) .and. (x_grid(i) .lt. 142.0d0)) then
          bckgrnd_vdc(i,j,k) = bckgrnd_vdc_ban
        endif

        bckgrnd_vvc(i,j,k) = Prandtl*bckgrnd_vdc(i,j,k)
      end do ! nx
    end do ! ny

    do k = 2, nz
      bckgrnd_vdc(:,:,k) = bckgrnd_vdc(:,:,1)
      bckgrnd_vvc(:,:,k) = bckgrnd_vvc(:,:,1)
    end do
  else ! .not. lhoriz_varying_bckgrnd
    do k = 1, nz
      zw(k) = z_face(k+1)
      bckgrnd_vdc(:,:,k) = bckgrnd_vdc1 + bckgrnd_vdc2* &
                      datan(bckgrnd_vdc_linv*           &
                            (zw(k)-bckgrnd_vdc_dpth))
      bckgrnd_vvc(:,:,k) = Prandtl*bckgrnd_vdc(:,:,k)
    end do
  endif ! lhoriz_varying_bckgrnd

! initialize grid info
  zgrid(0) = eps
  zgrid(1:nz) = -z_grid(1:nz)
  zgrid(nzf) = -z_face(nzf)

  hwide(0) = eps
  hwide(1:nz) = dz(1:nz)
  hwide(nzf) = eps

! initialize cg and Vtc
  Vtc = dsqrt(0.2d0/c_s/epssfc)/vonkar**2
  cg  = cstar*vonkar*(c_s*vonkar*epssfc)**0.33d0

! stokes ratio
  do j = 1, ny
    fstokes(:,ny) = 11.0d0 - max(5.0d0*dcos(3.0d0*(y_grid(j)/180.d0*pi)),0.d0)
  end do

! coriolis coeff
  do j = 1, ny
    fcort(j) = dsin(y_grid(j)*d2r)*pi/2.16d4
  end do

  if(ltidal_mixing) then
    tidal_diff = 0.d0
    tidal_coef = 0.d0
    rho_fw     = 1.0d0

    HT = 0.d0
    do j = 1, ny
      do i = 1, nx
        k = int(kb(i,j))
        HT(i,j) = z_face(k+1)
      end do
    end do
!! tidal_energy now move to ./src/driver.F90 for model initialization.
!! (see. set windmix)
!-----------------------------------------------------------------------
!  convert TIDAL_ENERGY_FLUX from W/m^2 to gr/s^3.
!-----------------------------------------------------------------------
!     tidal_energy_flux(:,:) = 1000.d0 * tidal_energy_flux(:,:)
!-----------------------------------------------------------------------
!  compute the time independent part of the tidal mixing coefficients
!-----------------------------------------------------------------------

    work1 = 0.d0
    do k = 1, nz
      where ( k < int(kb(1:nx,1:ny)) )
        work1(:,:) = work1(:,:) +  &
        dexp(-(HT(:,:) - z_face(k+1))/vertical_decay_scale)*dzw(k+1) !dzw(k)
      endwhere
    end do ! k

    tidal_coef(:,:,:) = 0.d0

    coef = local_mixing_fraction * mixing_efficiency / rho_fw

    do k = 1, nz
      where ( k <= int(kb(1:nx,1:ny)) )
        vertical_func(:,:) =                        &
             exp(-(HT(:,:) - z_face(k+1))           &
             /vertical_decay_scale) / work1(:,:)
        tidal_coef(:,:,k) = coef *  &
                            tidal_energy_flux(:,:) * vertical_func(:,:)
      elsewhere
        tidal_coef(:,:,k) = 0.d0
      endwhere
    end do ! k
    !call mpi_exch_3d_r8(tidal_coef)
 end if ! ltidal_mixing
end subroutine init_vmix_kpp

subroutine kppsrcglo(t,s,kpp_src)
  implicit none

  real(r8), dimension(nx,ny,nz), intent(inout) :: t, s
  real(r8), dimension(nx,ny,nz,2), intent(in) :: kpp_src
! KPP_SRC unit: T (degree C/s) S (g/kg/S)
! SO ONLY NEED DT CONVERSION
! NEED TO CHECK THE UNIT OF STF/DZ*VDC again

  t = t + dt*kpp_src(:,:,:,1)*in(1:nx,1:ny,1:nz)
  s = s + dt*kpp_src(:,:,:,2)*in(1:nx,1:ny,1:nz) !add the KPP NON-LOCAL SOURCE

end subroutine kppsrcglo

subroutine kppglo(u2, v2, t2, s2, rho, trans, &
                  stf, shf_qsw, smft, vbk, hbk, &
                  add, vdc, vvc, ev, hv, kpp_src, kpp_hblt)
  implicit none

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
!!DESCRIPTION:
! This is the main driver routine which calculates EV,HV and KPP_SRC
! The KPP mixing scheme is outlined in Large et al, Reviews of Geophysics, 32, 363 (November 1994).
! KPP_SRC part is treated in ==surface effect==

  trcr(:,:,:,1) = t2(1:nx,1:ny,1:nz)
  trcr(:,:,:,2) = s2(1:nx,1:ny,1:nz)*ppt_to_salt ! change to pop2 salt unit (g/g)
#ifdef chk_kpp
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK]: 01"
    write(*,*) "hv   = ",  hv(48,84,1:3)
    write(*,*) "vdct = ", vdc(48,84,1:3,1)
    write(*,*) "vdcs = ", vdc(48,84,1:3,2)
    write(*,*) "ev   = ",  ev(48,84,1:3)
    write(*,*) "vvc  = ", vvc(48,84,1:3)
    write(*,*) "rho  = ", rho(48,84,1:3)
    write(*,*) "hblt = ", kpp_hblt(48,84)
    write(*,*) "smft = ", smft(48,84,:)
    write(*,*) "t2   = ", trcr(48,84,1:3,1)
    write(*,*) "s2   = ", trcr(48,84,1:3,2)
    write(*,*) "u2   = ", u2(48,84,1:3)
    write(*,*) "v2   = ", v2(48,84,1:3)
    write(*,*) "emax = ", 0.5d0/(dt*odzw(2:4)**2)
    write(*,*) "odzw = ", odzw(2:4)*iw(48,84,2:4)
  end if
#endif
  call vmix_coeffs_kpp(trans, kpp_src, vdc, vvc, trcr, u2, v2, rho,&
                       stf, shf_qsw, kpp_hblt, convect_diff, convect_visc,&
                       smft=smft)
  do k = 1, nz-1
    l = k+1
    emax = 0.5d0/(dt*odzw(l)**2)
!TS            emax=MIN(1200.d0,emax)
            !emax=500.d0
    do j = 1, ny
      do i = 1, nx
!       if(k .ge. kbl(i+1,j+1))then
!         emax = 0.5d0/(dt*odzw(l)**2)
!         emax = MIN(6000.d0,emax)
!       else
!         emax=600.d0
!       end if
        dtmp = odzw(l)*iw(i,j,l)
        !! limit the value not to over 900 or .5/(DT*ODZW(L)**2)
        ev(i,j,k) = dtmp*MAX(MIN(emax,vvc(i,j,k)),vbk(k)*5.d0)
        !ev(i,j,k) = dtmp*hbk(k)*5.d0
        !if(hmix_tracer_type_gm < 0.d0) then
        hv(i,j,k,1) = dtmp*MAX(MIN(emax,vdc(i,j,k,1)),hbk(k)*5.d0)
        hv(i,j,k,2) = dtmp*MAX(MIN(emax,vdc(i,j,k,2)),hbk(k)*5.d0)
        !else
        !  hv(i,j,k) = dtmp*(vdc(i,j,k,1) + add(i,j,k))
       ! end if
      end do
    end do
  end do
#ifdef chk_kpp
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK]: 02"
    write(*,*) "hv   = ",  hv(48,84,1:3)
    write(*,*) "vdct = ", vdc(48,84,1:3,1)
    write(*,*) "vdcs = ", vdc(48,84,1:3,2)
    write(*,*) "ev   = ",  ev(48,84,1:3)
    write(*,*) "vvc  = ", vvc(48,84,1:3)
    write(*,*) "rho  = ", rho(48,84,1:3)
    write(*,*) "hblt = ", kpp_hblt(48,84)
    write(*,*) "smft = ", smft(48,84,:)
    write(*,*) "t2   = ", trcr(48,84,1:3,1)
    write(*,*) "s2   = ", trcr(48,84,1:3,2)
    write(*,*) "u2   = ", u2(48,84,1:3)
    write(*,*) "v2   = ", v2(48,84,1:3)
    write(*,*) "emax = ", 0.5d0/(dt*odzw(1:3)**2)
    write(*,*) "odzw = ", odzw(2:4)*iw(48,84,2:4)
  end if
#endif
end subroutine kppglo

subroutine vmix_coeffs_kpp(trans, kpp_src, vdc, vvc, trcr, uuu, vvv, rhomix,&
                           stf, shf_qsw, kpp_hblt, convect_diff, convect_visc,&
                           smf, smft)
! !DESCRIPTION:
!  This is the driver routine which calculates the vertical
!  mixing coefficients for the KPP mixing scheme as outlined in
!  Large, McWilliams and Doney, Reviews of Geophysics, 32, 363
!  (November 1994).  The non-local mixing is also computed here, but
!  is treated as a source term in baroclinic.
  implicit none
! !INPUT PARAMETERS:
!  integer(2), dimension(nx,ny,nz), intent(in) :: IN

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
!-----------------------------------------------------------------------
!  compute buoyancy differences at each vertical level.
!  derived variables: dbloc, dbsfc
!-----------------------------------------------------------------------
  call buoydiff(dbloc, dbsfc, trcr)
!-----------------------------------------------------------------------
!  compute mixing due to shear instability, internal waves and
!  convection
!  derived variables: vdc,visc
!-----------------------------------------------------------------------
  call ri_iwmix(dbloc, visc, vdc, uuu, vvv)
!-----------------------------------------------------------------------
!  compute double diffusion if desired
!-----------------------------------------------------------------------
  if(ldbl_diff)  call ddmix(vdc,trcr)
!-----------------------------------------------------------------------
!  compute boundary layer depth
!  derived variables: kpp_hblt, ustar, bfsfc, stable, kbl
!-----------------------------------------------------------------------
  if(present(smft)) then
    call bldepth(trans,dbloc, dbsfc, trcr, uuu, vvv, stf, shf_qsw, &
                 kpp_hblt, ustar, bfsfc, stable, kbl, &
                 smft = smft)
  else
    call bldepth(trans,dbloc, dbsfc, trcr, uuu, vvv, stf, shf_qsw, &
                 kpp_hblt, ustar, bfsfc, stable, kbl, &
                 smf = smf)
  end if
!-----------------------------------------------------------------------
!  compute boundary layer diffusivities
!  derived variables: visc,vdc,ghat
!-----------------------------------------------------------------------
  call blmix(visc, vdc, kpp_hblt, ustar, bfsfc, stable, &
             kbl, ghat)

#ifdef chk_kpp1
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK vmix]: 01"
    write(*,*) "kbl = ", kbl(48,84)
    write(*,*) "hblt = ", kpp_hblt(48,84)
    write(*,*) achar(10)
  end if
#endif
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
  do k = 1, nz-1
    work1 = dbloc(:,:,k)/(zgrid(k) - zgrid(k+1))
    if(BVSQcon /= 0.d0) then
      work2 = min(1.d0-(max(work1,BVSQcon))/BVSQcon, 1.d0)
      fcon  = (1.d0 - work2*work2)**3
    else
      where (work1 > 0.d0)
        fcon = 0.d0
      elsewhere
        fcon = 1.d0
      end where
    end if

    !*** add convection and reset sea floor values to zero
    do j = 1, ny
      do i = 1, nx
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
      end do
    end do

    vvc(:,:,k) = merge(visc(:,:,k), 0.d0, ( k < kb(1:nx,1:ny) ) )
  end do

  VDC(:,:,nz,:) = 0.d0
  VVC(:,:,nz)   = 0.d0
#ifdef chk_kpp1
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK vmix]: 02"
    write(*,*) "kbl = ", kbl(48,84)
    write(*,*) "hblt = ", kpp_hblt(48,84)
    write(*,*) achar(10)
  end if
#endif
!-----------------------------------------------------------------------
!  add ghatp term from previous computation to right-hand-side
!  source term on current row
!-----------------------------------------------------------------------  
  do n = 1, 2
    mt2 = min(n,2)
    KPP_SRC(:,:,1,n) = STF(:,:,n)*odz(1)                &
                     *(-VDC(:,:,1,mt2)*GHAT(:,:,1))
    do k = 2, nz
      KPP_SRC(:,:,k,n) = STF(:,:,n)*odz(k)              &
                      *( VDC(:,:,k-1,mt2)*GHAT(:,:,k-1) &
                        -VDC(:,:,k  ,mt2)*GHAT(:,:,k  ))
#ifdef chk_kpp1
       if(myid_x .eq. 8 .and. myid_y .eq. 3) then
         write(*,*) "[KPP CHECK vmix]: 02-1", n, k
         write(*,*) "kbl = ", kbl(48,84)
         write(*,*) "hblt = ", kpp_hblt(48,84)
         write(*,*) achar(10)
       end if
#endif
    end do
  end do
  KPP_SRC(:,:,:,2) = 1000.d0*KPP_SRC(:,:,:,2) !from msu to psu

!-----------------------------------------------------------------------
!
!  compute diagnostic mixed layer depth (cm) using a max buoyancy
!  gradient criterion.  Use USTAR and BFSFC as temps.
!
!-----------------------------------------------------------------------
#ifdef chk_kpp1
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK vmix]: 03"
    write(*,*) "kbl = ", kbl(48,84)
    write(*,*) "hblt = ", kpp_hblt(48,84)
    write(*,*) achar(10)
  end if
#endif
  ustar = 0.d0
  where (kb(1:nx,1:ny) == 1)
    hmxl(:,:) = z_grid(1)-z_face(1)
  elsewhere
    hmxl(:,:) = 0.d0
  endwhere

  do k = 2, nz
    where(k <= kb(1:nx,1:ny))
      ustar = max(dbsfc(:,:,k)/(z_grid(k)-z_face(1)),ustar)
      hmxl(:,:) = z_grid(k)-z_face(1)
    endwhere
  end do

  visc(:,:,1) = 0.d0
  do k = 2, nz
    where(ustar > 0.d0 )
      visc(:,:,k) = (dbsfc(:,:,k)-dbsfc(:,:,k-1))/ &
                         (z_grid(k) - z_grid(k-1))
    end where
    where(visc(:,:,k) >= ustar .and.              &
         (visc(:,:,k)-visc(:,:,k-1)) /= 0.d0 .and.  &
          ustar > 0.d0 )   ! avoid divide by zero
      bfsfc = (visc(:,:,k)-ustar)/ &
              (visc(:,:,k)-visc(:,:,k-1))
      hmxl(:,:) = -0.5d0*(zgrid(k  ) + zgrid(k-1))*(1.0d0-bfsfc) &
                  -0.5d0*(zgrid(k-1) + zgrid(k-2))*bfsfc
      ustar(:,:) = 0.d0
    endwhere
  end do
#ifdef chk_kpp1
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK vmix]: 04"
    write(*,*) "hblt = ", kpp_hblt(48,84)
    write(*,*) achar(10)
  end if
#endif
end subroutine vmix_coeffs_kpp

!--------------------------------------------------------------
subroutine buoydiff(dbloc, dbsfc, trcr)
  implicit none
! !DESCRIPTION:
!  This routine calculates the buoyancy differences at model levels.
!
! !REVISION HISTORY:
!  same as module
!   implicit none
! !INPUT PARAMETERS:
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

!-----------------------------------------------------------------------
!  calculate density and buoyancy differences at surface
!-----------------------------------------------------------------------
  tempsfc = merge(-2.d0,trcr(:,:,1,1),trcr(:,:,1,1) < -2.d0)
  klvl  = 2
  kprev = 1

  tempk(:,:,kprev) = tempsfc
  dbsfc(:,:,1) = 0.d0

!-----------------------------------------------------------------------
!  calculate DBLOC and DBSFC for all other levels
!-----------------------------------------------------------------------\
  do k = 2, nz
    tempk(:,:,klvl) = merge(-2.0d0,trcr(:,:,k,1),trcr(:,:,k,1) < -2.0d0)

! calculate the density from Temperature and Salinity using an equation of state
! derived from McDougall, Wright, Jackett and Feistel (hereafter MWJF, 2001 )
    call state(k, k, tempsfc,          trcr(:,:,1,  2), rhofull=rho1)
    call state(k, k, tempk(:,:,kprev), trcr(:,:,k-1,2), rhofull=rhokm)
    call state(k, k, tempk(:,:,klvl),  trcr(:,:,k  ,2), rhofull=rhok)
    do j = 1, ny
      do i = 1, nx
        if(rhok(i,j) /= 0.d0) then
          dbsfc(i,j,k)   = grav*(1.0d0 - rho1 (i,j)/rhok(i,j))
          dbloc(i,j,k-1) = grav*(1.0d0 - rhokm(i,j)/rhok(i,j))
        else
          dbsfc(i,j,k)   = 0.d0
          dbloc(i,j,k-1) = 0.d0
        end if
        if ( k-1 >= kb(i,j) ) dbloc(i,j,k-1) = 0.d0
      end do
    end do
    ktmp  = klvl
    klvl  = kprev
    kprev = ktmp
  end do

  dbloc(:,:,nz) = 0.d0

end subroutine buoydiff

subroutine state(k, kk, tempk, saltk, &
                 rhoout, rhofull, drhodt, drhods)
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
!  if requested (ie the optional arguments are present).
!
!  If $k = kk$ are equal the density for level k is returned.
!  If $k \neq kk$ the density returned is that for a parcel
!  adiabatically displaced from level k to level kk.

! !INPUT PARAMETERS:
  integer, intent(in) :: &
    k,                   &! depth level index
    kk                    ! level to which water is adiabatically
                          ! displaced
  real(r8), dimension(nx,ny), intent(in) :: &
    tempk,               &! temperature at level k
    saltk                 ! salinity    at level k
! !OUTPUT PARAMETERS:

  real(r8), dimension(nx,ny), optional, intent(out) :: &
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
  real(r8), dimension(nx,ny) :: &
    TQ, SQ,             &! adjusted T,S
    sqr, denomk,        &! work arrays
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
  TQ = max(TQ,   -2.0d0)
  SQ = min(saltk, 0.999d0)
  SQ = max(SQ, 0.d0)

!------------------------------------------------------------
  do dep = 1, nz
    pressz(dep) = pressure(-zgrid(dep)*0.01d0) ! cm to m
  end do
  ppz = 10.d0*pressz(kk)
  SQ  = 1000.d0*SQ
!#ifdef CCSMCOUPLED
! call shr_vmath_sqrt(SQ, SQR, nx_block*ny_block)
!#else
  SQR = sqrt(SQ)
!#endif

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

  if(present(RHOOUT)) then
    RHOOUT = WORK1*DENOMK
  end if

  if(present(RHOFULL)) then
    RHOFULL = WORK1*DENOMK
  end if

  if(present(DRHODT)) then
    WORK3 = &! dP_1/dT
            mwjfnums0t1 + TQ * (2.d0*mwjfnums0t2 +    &
            3.d0*mwjfnums0t3 * TQ) + mwjfnums1t1 * SQ
    WORK4 = &! dP_2/dT
            mwjfdens0t1 + SQ * mwjfdens1t1 +               &
            TQ * (2.d0*(mwjfdens0t2 + SQ*SQR*mwjfdensqt2) +  &
            TQ * (3.d0*(mwjfdens0t3 + SQ * mwjfdens1t3) +    &
            TQ *  4.d0*mwjfdens0t4))

    DRHODT = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK
  end if

  if(present(DRHODS)) then
    WORK3 = &! dP_1/dS
            mwjfnums1t0 + mwjfnums1t1 * TQ + 2.d0*mwjfnums2t0 * SQ

    WORK4 = mwjfdens1t0 +   &! dP_2/dS
            TQ * (mwjfdens1t1 + TQ*TQ*mwjfdens1t3) +   &
            1.5d0*SQR*(mwjfdensqt0 + TQ*TQ*mwjfdensqt2)

    DRHODS = (WORK3 - WORK1*DENOMK*WORK4)*DENOMK * 1000.d0 !unit converison
  end if

end subroutine state

function pressure(depth)
! !INPUT PARAMETERS:
  real(r8), intent(in) :: depth    ! depth in meters
! !OUTPUT PARAMETERS:
  real(r8) :: pressure   ! pressure [bars]

  pressure = 0.059808d0*(dexp(-0.025d0*depth) - 1.0d0) &
           + 0.100766d0*depth + 2.28405d-7*depth**2
end function pressure

subroutine ri_iwmix(dbloc, visc, vdc, uuu, vvv)
! !DESCRIPTION:
!  Computes viscosity and diffusivity coefficients for the interior
!  ocean due to shear instability (richardson number dependent),
!  internal wave activity, and to static instability (Ri < 0).
  implicit none
! !INPUT PARAMETERS:
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
!-----------------------------------------------------------------------
!  compute mixing at each level
!-----------------------------------------------------------------------
  kvmix       = 0.d0
  kvmix_m     = 0.d0
  visc(:,:,0) = 0.d0

  do k = 1, nz
!-----------------------------------------------------------------------
!     compute velocity shear squared and average to T points:
!     VSHEAR = (UUU(k)-UUU(k+1))**2+(VVV(k)-VVV(k+1))**2
!     Use FRI here as a temporary.
!-----------------------------------------------------------------------
    if(k < nz) then
      fri = (uuu(1:nx,1:ny,k)-uuu(1:nx,1:ny,k+1))**2 + &
            (vvv(1:nx,1:ny,k)-vvv(1:nx,1:ny,k+1))**2
      vshear = fri
    else
      vshear = 0.d0
    end if
!-----------------------------------------------------------------------
!     compute local richardson number
!     use visc array as temporary Ri storage to be smoothed
!-----------------------------------------------------------------------
    ri_loc = dbloc(:,:,k)*(zgrid(k)-zgrid(k+1))/(vshear + eps)
    visc(:,:,k)   = merge( ri_loc, visc(:,:,k-1), k <= int(kb(1:nx,1:ny)) )
  end do

!-----------------------------------------------------------------------
!  vertically smooth Ri num_v_smooth_Ri times with 1-2-1 weighting
!  result again stored temporarily in VISC and use RI_LOC and FRI
!  as temps
!----------------------------------------------------------------------
  do n = 1, num_v_smooth_Ri
    fri            = 0.25d0*visc(:,:,1)
    visc(:,:,nz+1) = visc(:,:,nz)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ri_loc(i,j) = visc(i,j,k)
          if(int(kb(i,j)) >= 2 ) then !modified from KMT, land=0
             visc(i,j,k) = fri(i,j) + 0.5d0*ri_loc(i,j) &
                                      + 0.25d0*visc(i,j,k+1)
          end if
          fri(i,j) = 0.25d0*ri_loc(i,j)
        end do
      end do
    end do
  end do
!-----------------------------------------------------------------------
!  now that we have a smoothed Ri field, finish computing coeffs
!  at each level
!-----------------------------------------------------------------------
  if(ltidal_mixing) tidal_diff(:,:,:) = 0.d0
  do k = 1, nz
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
    if(ltidal_mixing ) then

!-----------------------------------------------------------------------
!
!  consider the internal wave mixing first. rich_mix is used as the
!  upper limit for internal wave mixing coefficient. bckgrnd_vvc
!  was already multiplied by Prandtl.
!
!  NOTE: no partial_bottom_cell implementation at this time
!
!-----------------------------------------------------------------------

      work1 = dbloc(:,:,k)/(zgrid(k) - zgrid(k+1))
      where (work1 > 0.d0)
        tidal_diff(:,:,k) = tidal_coef(:,:,k)/work1
      endwhere
      ! Notes:
      ! (1) this step breaks backwards compatibility
      ! (2) check for k>2 was added to if statement to
      ! avoid out of bounds access
      if(.not. lccsm_control_compatible) then   ! this step breaks backwards compatibility
        if(k > 2) then
          where((k ==kb(1:nx,1:ny)-1 .or. k == kb(1:nx,1:ny)-2))
            tidal_diff(:,:,k) = max( tidal_diff(:,:,k),  &
                                     tidal_diff(:,:,k-1))
          endwhere
        end if
      end if

      work1 = Prandtl*min(bckgrnd_vvc(:,:,k)/Prandtl  &
                         + tidal_diff(:,:,k), tidal_mix_max)
      if(k < nz) then
        kvmix_m(:,:) = work1(:,:)
      end if
      if(k < nz) then
        vdc(:,:,k,2) = min(bckgrnd_vdc(:,:,k) + tidal_diff(:,:,k),  &
                             tidal_mix_max)
        kvmix(:,:) = vdc(:,:,k,2)
      end if

      if(lrich) then
        fri = min((max(visc(:,:,k),0.d0))/Riinfty, 1.d0)
        visc(:,:,k) = work1 + rich_mix*(1.d0 - fri*fri)**3
        if(k < nz) then
          vdc(:,:,k,2) = vdc(:,:,k,2) + rich_mix*(1.d0 - fri*fri)**3
          vdc(:,:,k,1) = vdc(:,:,k,2)
        end if
      else
        visc(:,:,k) = work1
        if(k < nz) then
          vdc(:,:,k,1) = vdc(:,:,k,2)
        end if
      end if
    else ! .not. ltidal_mixing
      if(k < nz) then
        kvmix(:,:) = bckgrnd_vdc(:,:,k)
        kvmix_m(:,:) = bckgrnd_vvc(:,:,k)
      end if

      if(lrich) then
        fri = min((max(visc(:,:,k),0.d0))/Riinfty, 1.0d0)
        visc(:,:,k) = bckgrnd_vvc(:,:,k) + &
                      rich_mix*(1.0d0 - fri*fri)**3
        if(k < nz) then
          vdc(:,:,k,2) = bckgrnd_vdc(:,:,k) + &
                              rich_mix*(1.0d0 - fri*fri)**3
          vdc(:,:,k,1) = vdc(:,:,k,2)
        end if
      else
        visc(:,:,k) = bckgrnd_vvc(:,:,k)

        if(k < nz) then
          vdc(:,:,k,2) = bckgrnd_vdc(:,:,k)
          vdc(:,:,k,1) = vdc(:,:,k,2)
        end if
      end if
    end if ! ltidal_mixing
!-----------------------------------------------------------------------
!     set seafloor values to zero
!-----------------------------------------------------------------------
    do j = 1, ny
      do i = 1, nx
        if(k >= int(kb(i,j))) then
          visc(i,j,k  ) = 0.d0
          vdc (i,j,k,1) = 0.d0
          vdc (i,j,k,2) = 0.d0
        end if
      end do
    end do
  end do
!-----------------------------------------------------------------------
!  fill extra coefficients for blmix
!-----------------------------------------------------------------------
  visc(:,:,0  ) = 0.d0
  vdc (:,:,0,:) = 0.d0
  visc(:,:,nzf  ) = 0.d0
  vdc (:,:,nzf,:) = 0.d0

end subroutine ri_iwmix

subroutine ddmix(VDC, TRCR)
! !DESCRIPTION:
!  $R_\rho$ dependent interior flux parameterization.
!  Add double-diffusion diffusivities to Ri-mix values at blending
!  interface and below.
  implicit none
!!INPUT PARAMETERS:
  real(r8), dimension(nx,ny,nz,2), intent(in) :: &
    TRCR                ! tracers at current time
! !INPUT/OUTPUT PARAMETERS:
  real(r8), dimension(nx,ny,0:nzf,2),intent(inout) :: &
    VDC        ! diffusivity for tracer diffusion
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
  integer ::  k, kup, knxt
  real(r8), dimension(nx,ny) :: &
    ALPHADT,           &! alpha*DT  across interfaces
    BETADS,            &! beta *DS  across interfaces
    RRHO,              &! dd density ratio
    DIFFDD,            &! dd diffusivity scale
    prtl_tmp             ! prandtl number
  real(r8), dimension(nx,ny,2) :: &
    TALPHA,            &! temperature expansion coeff
    SBETA               ! salinity    expansion coeff
!-----------------------------------------------------------------------
!  compute alpha*DT and beta*DS at interfaces.  use RRHO and
!  PRANDTL for temporary storage for call to state
!-----------------------------------------------------------------------
  kup  = 1
  knxt = 2

  prtl_tmp = merge(-2.0d0,TRCR(:,:,1,1),TRCR(:,:,1,1) < -2.0d0)
  call state(1, 1, prtl_tmp, TRCR(:,:,1,2),RHOFULL=RRHO, &
                   DRHODT=TALPHA(:,:,kup), DRHODS=SBETA(:,:,kup))

  do k = 1, nz
    if(k < nz) then
      prtl_tmp = merge(-2.0d0,TRCR(:,:,k+1,1),TRCR(:,:,k+1,1) < -2.0d0)
      call state(k+1, k+1, prtl_tmp, TRCR(:,:,k+1,2),  &
                           RHOFULL=RRHO, DRHODT=TALPHA(:,:,knxt), &
                           DRHODS= SBETA(:,:,knxt))

      ALPHADT = -0.5d0*(TALPHA(:,:,kup) + TALPHA(:,:,knxt)) &
                      *(TRCR(:,:,k,1) - TRCR(:,:,k+1,1))

      BETADS  = 0.5d0*( SBETA(:,:,kup) +  SBETA(:,:,knxt)) &
                     *(TRCR(:,:,k,2) - TRCR(:,:,k+1,2))
      kup  = knxt
      knxt = 3 - kup
    else
      ALPHADT = 0.d0
      BETADS  = 0.d0
    end if

!-----------------------------------------------------------------------
!     salt fingering case
!-----------------------------------------------------------------------
    where( ALPHADT > BETADS .and. BETADS > 0.d0 )
      RRHO = MIN(ALPHADT/BETADS, Rrho0)
      DIFFDD = dsfmax*(1.0d0-(RRHO-1.0d0)/(Rrho0-1.0d0))**3
      VDC(:,:,k,1) = VDC(:,:,k,1) + 0.7d0*DIFFDD
      VDC(:,:,k,2) = VDC(:,:,k,2) + DIFFDD
    endwhere

!-----------------------------------------------------------------------
!     diffusive convection
!-----------------------------------------------------------------------
    where( ALPHADT < 0.d0 .and. BETADS < 0.d0 .and. ALPHADT > BETADS )
      RRHO   = ALPHADT / BETADS
      DIFFDD = 1.5d-2*0.909d0* &
          dexp(4.6d0*dexp(-0.54d0*(1.0d0/RRHO-1.0d0)))
         prtl_tmp = 0.15d0*RRHO
    elsewhere
      RRHO     = 0.d0
      DIFFDD   = 0.d0
      prtl_tmp = 0.d0
    endwhere

    where(RRHO > 0.5d0) prtl_tmp = (1.85d0 - 0.85d0/RRHO)*RRHO

    VDC(:,:,k,1) = VDC(:,:,k,1) + DIFFDD
    VDC(:,:,k,2) = VDC(:,:,k,2) + prtl_tmp*DIFFDD
  end do

end subroutine ddmix

subroutine bldepth(trans, dbloc, dbsfc, trcr, uuu, vvv, stf, shf_qsw,&
                    hblt, ustar, bfsfc, stable, kbl, smf, smft)
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
  implicit none
! !INPUT PARAMETERS:
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

  real(r8), dimension(nx,ny,2), intent(in), optional :: &
    smf,           &! surface momentum forcing at U points
    smft            ! surface momentum forcing at T points
                    ! *** either one or the other (not
                    ! *** both) should be passed
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
    hlimit,test              ! limit to mixed-layer depth
                          ! (= min(HEKMAN,HMONOB))

  real(r8), dimension(nx,ny,3) :: &
    ri_bulk,           &! Bulk Ri number at 3 lvls
    hmonob              ! Monin-Obukhov depth limit

  real(r8) ::          &
    absorb_frac,       &! shortwave absorption frac
    sqrt_arg,          &! dummy sqrt argument
    z_upper, z_up,     &! upper depths for RI_BULK interpolation
    work1(-1:nx+2,-1:ny+2)

  real(r8) :: &
    a_co, b_co, c_co    ! coefficients of the quadratic equation
                        ! $(a_{co}z^2+b_{co}|z|+c_{co}=Ri_b) used to
                        ! find the boundary layer depth. when
                        ! finding the roots, c_co = c_co - Ricr
  real(r8) :: &
    slope_up            ! slope of the above quadratic equation
                        ! at zup. this is used as a boundary
                        ! condition to determine the coefficients.
!-----------------------------------------------------------------------
!  compute friction velocity USTAR.
!-----------------------------------------------------------------------
  if(present(smft)) then
    ustar = sqrt(sqrt(smft(:,:,1)**2 + smft(:,:,2)**2))
  else
    work = sqrt(sqrt(smf(:,:,1)**2 + smf(:,:,2)**2))
    ustar = work
  end if
!-----------------------------------------------------------------------
!  compute density and expansion coefficients at surface
!-----------------------------------------------------------------------
  work = merge(-2.0d0,trcr(:,:,1,1),trcr(:,:,1,1) < -2.0d0)
  call state(1,1,work,trcr(:,:,1,2), &
             rhofull=rho1, drhodt=talpha, drhods=sbeta)
#ifdef chk_kpp1
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK bldepth]: 01"
    write(*,*) "ustar = ", ustar(48,84)
    write(*,*) "rho1 = ", rho1(48,84)
    write(*,*) "t2   = ", trcr(48,84,1:3,1)
    write(*,*) "s2   = ", trcr(48,84,1:3,2)
    write(*,*) achar(10)
  end if
#endif
!-----------------------------------------------------------------------
!  compute turbulent and radiative sfc buoyancy forcing
!-----------------------------------------------------------------------
  do j = 1, ny
    do i = 1, nx
      if(rho1(i,j) /= 0.d0) then
        bo(i,j) = grav*(-talpha(i,j)*stf(i,j,1) - &
                          sbeta(i,j)*stf(i,j,2))/rho1(i,j)

        bosol(i,j) = -grav*talpha(i,j)*shf_qsw(i,j)/rho1(i,j)
      else
        bo(i,j) = 0.d0
        bosol(i,j) = 0.d0
      end if
    end do
  end do
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

  ri_bulk(:,:,kupper) = 0.d0
  ri_bulk(:,:,kup) = 0.d0
  kbl = merge( int(kb(1:nx,1:ny)), 1, (int(kb(1:nx,1:ny)) > 1))
  hlangm = 0.d0

  do kl = 1, nz
    zkl = -zgrid(kl)
    do j = 1, ny
      do i = 1, nx
        if(kl == kbl(i,j)) hblt(i,j) = zkl(i,j)
      end do
    end do
  end do
#ifdef chk_kpp1
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK bldepth]: 02"
    write(*,*) "kbl  = ", kbl(48,84)
    write(*,*) "hblt = ", hblt(48,84)
    write(*,*) achar(10)
  end if
#endif
  if(lcheckekmo) then
    hekman = -zgrid(nz) + eps
    hlimit = -zgrid(nz) + eps
    if(lshort_wave) then
      bfsfc = bo + bosol*(1.d0-trans(:,:,1))
    else
      bfsfc = bo
    end if

    stable = merge(1.d0, 0.d0, bfsfc >= 0.d0)

    bfsfc  = bfsfc + stable*eps

    work   =   stable * cmonob*ustar*ustar*ustar/vonkar/bfsfc &
            + (stable -1.0d0)*zgrid(nz)
    hmonob(:,:,kup) = merge( -z_up+eps, work, work <= -z_up )
  end if
  rsh_hblt = 0.d0

!-----------------------------------------------------------------------
!  compute velocity shear squared on U-grid and use the maximum
!  of the four surrounding U-grid values for the T-grid.
!-----------------------------------------------------------------------
  do kl = 2, nz ! loop, 1846
    work1 = (uuu(:,:,1)-uuu(:,:,kl))**2 + &
            (vvv(:,:,1)-vvv(:,:,kl))**2
    ZKL = -zgrid(kl)

    do j = 1, ny
      do i = 1, nx
        vshear(i,j) = max(work1(i,j  ), work1(i-1,j  ),   &
                          work1(i,j-1), work1(i-1,j-1))
      end do
    end do
!-----------------------------------------------------------------------
!     compute bfsfc= Bo + radiative contribution down to hbf * hbl
!     add epsilon to BFSFC to ensure never = 0
!-----------------------------------------------------------------------
    if(lshort_wave) then
      bfsfc = bo + bosol*(1.d0-trans(:,:,kl-1))
    else
      bfsfc = bo
    end if
    
    stable = merge(1.d0, 0.d0, bfsfc >= 0.d0)
    bfsfc  = bfsfc + stable*eps
!-----------------------------------------------------------------------
!     compute the Ekman and Monin Obukhov depths using above stability
!-----------------------------------------------------------------------
    if(lcheckekmo) then
      do j = 1, ny
        do i = 1, nx
          if(stable(i,j) > 0.5d0 .and. hekman(i,j) >= -zgrid(nz) ) then
             hekman(i,j) = max(zkl(i,j), &
                         cekman*ustar(i,j)/(dabs(fcort(j))+eps))
          end if
        end do
      end do

      hmonob(:,:,kdn) = stable*cmonob*ustar*ustar*ustar/vonkar/bfsfc + &
                       (stable-1.0d0)*zgrid(nz)
      do j = 1, ny
        do i = 1, nx
          if(hmonob(i,j,kdn) <= zkl(i,j) .and. &
             hmonob(i,j,kup) >  -z_up) then
            work(i,j) = (hmonob(i,j,kdn) - hmonob(i,j,kup))/ &
                        (z_up + zkl(i,j))
            hlimit(i,j) = (hmonob(i,j,kdn) - work(i,j)*zkl(i,j))/ &
                          (1.0d0 - work(i,j))
          end if
        end do
      end do
    end if
!-----------------------------------------------------------------------
!     compute velocity scales at sigma, for hbl = -zgrid(kl)
!-----------------------------------------------------------------------
    sigma = epssfc
    call wscale(sigma, zkl, ustar, bfsfc, 2, wm, ws)
!-----------------------------------------------------------------------
!     compute the turbulent shear contribution to RI_BULK and store
!     in WM.
!-----------------------------------------------------------------------
    b_frqncy = dsqrt( &
               0.5d0*(dbloc(:,:,kl) + abs(dbloc(:,:,kl)) + eps2)/  &
               (zgrid(kl)-zgrid(kl+1)) )

    wm = zkl*ws*b_frqncy* &
        ((vtc/Ricr)*max(2.1d0 - 200.d0*b_frqncy,concv) )
!-----------------------------------------------------------------------
!     compute bulk Richardson number at new level
!-----------------------------------------------------------------------
    work = merge((zgrid(1)-zgrid(kl))*dbsfc(:,:,kl), &
                0.d0, kb(1:nx,1:ny) >= kl)
    if(linertial) then
      ri_bulk(:,:,kdn) = work/(vshear+wm+ustar*bolus_sp(:,:)+eps)
    else
      ri_bulk(:,:,kdn) = work/(vshear+wm+eps)
    end if
!-----------------------------------------------------------------------
!       find hbl where Rib = Ricr. if possible, use a quadratic
!       interpolation. if not, linearly interpolate. the quadratic
!       equation coefficients are determined using the slope and
!       Ri_bulk at z_up and Ri_bulk at zgrid(kl). the slope at
!       z_up is computed linearly between z_upper and z_up.
!       compute Langmuir depth always
!-----------------------------------------------------------------------
    do j = 1, ny
      do i = 1, nx
        if(kbl(i,j) == kb(i,j) .and. ri_bulk(i,j,kdn) > Ricr) then
          slope_up = (RI_BULK(i,j,kupper) - RI_BULK(i,j,kup))/(z_up - z_upper)

          a_co = (RI_BULK(i,j,kdn) - RI_BULK(i,j,kup) -         &
                  slope_up*(ZKL(i,j) + z_up) )/(z_up + ZKL(i,j))**2

          b_co = slope_up + 2.d0 * a_co * z_up

          c_co = RI_BULK(i,j,kup) + z_up*(a_co*z_up + slope_up) - Ricr

          sqrt_arg = b_co**2 - 4.d0*a_co*c_co

          if((dabs(b_co) > eps .and. dabs(a_co)/dabs(b_co) <= eps ) &
              .or. sqrt_arg <= 0.d0 ) then

            HBLT(i,j) = -z_up + (z_up + ZKL(i,j)) *               &
                        (Ricr             - RI_BULK(i,j,kup))/    &
                        (RI_BULK(i,j,kdn) - RI_BULK(i,j,kup))

          else
            HBLT(i,j) = (-b_co + dsqrt(sqrt_arg)) / (2.d0*a_co)
          end if
          kbl(i,j) = kl
          RSH_HBLT(i,j) =  (VSHEAR(i,j)*Ricr/ &
                           (DBSFC(i,j,kl)+eps))/HBLT(i,j)

          HLANGM(i,j) = USTAR(i,j) * DSQRT(FSTOKES(i,j)*ZKL(i,j)/(DBSFC(i,j,kl)+eps))/0.9d0
        end if
      end do
    end do
#ifdef chk_kpp1
    if(myid_x .eq. 8 .and. myid_y .eq. 3) then
      write(*,*) "[KPP CHECK bldepth]: 03"
      write(*,*) "level k =", kl
      write(*,*) "kbl  = ", kbl(48,84)
      write(*,*) "hblt = ", hblt(48,84)
    end if
#endif
!-----------------------------------------------------------------------
!     swap klevel indices and move to next level
!-----------------------------------------------------------------------
    ktmp   = kupper
    kupper = kup
    kup    = kdn
    kdn    = ktmp
    z_upper = z_up
    z_up    = zgrid(kl)

  end do ! end loop, 1846
!-----------------------------------------------------------------------
!     apply Langmuir parameterization if requested
!-----------------------------------------------------------------------
  if(llangmuir) then
    do kl = nz,2,-1
      where(hlangm > hblt          .and. &
            hlangm >  -zgrid(kl-1) .and. &
            hlangm <= zkl                )
        hblt = hlangm
        kbl  = kl
      end where
    end do
  end if
!-----------------------------------------------------------------------
!  first combine Ekman and Monin-Obukhov depth limits. then apply
!  these restrictions to HBLT. note that HLIMIT is set to -zgrid(km)
!  in unstable forcing.
!-----------------------------------------------------------------------
  if(lcheckekmo) then
    where(hekman < hlimit) hlimit = hekman
    do kl = 2, nz
      where(hlimit < hblt         .and. &
            hlimit >  -zgrid(kl-1).and. &
            hlimit <= zkl               )
        hblt = hlimit
        kbl = kl
      end where
    end do
  end if
!-----------------------------------------------------------------------
!  apply a Gaussian filter
!-----------------------------------------------------------------------
#ifdef chk_kpp1
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK bldepth]: 04"
    write(*,*) "kbl  = ", kbl(48,84)
    write(*,*) "hblt = ", hblt(48,84)
    write(*,*) achar(10)
  end if
#endif
  call smooth_hblt(.true., .false., hblt=hblt, kbl=kbl)
#ifdef chk_kpp1
  if(myid_x .eq. 8 .and. myid_y .eq. 3) then
    write(*,*) "[KPP CHECK bldepth]: 05"
    write(*,*) "kbl  = ", kbl(48,84)
    write(*,*) "hblt = ", hblt(48,84)
    write(*,*) achar(10)
  end if
#endif
!-----------------------------------------------------------------------
!  correct stability and buoyancy forcing for SW up to boundary layer
!-----------------------------------------------------------------------
  if(lshort_wave)  bfsfc = bo + bosol*(1.0d0-trans(:,:,1))
  
  stable = merge(1.0d0, 0.d0, bfsfc >= 0.d0)
  bfsfc  = bfsfc + stable * eps ! ensures bfsfc never=0

end subroutine bldepth

subroutine wscale(sigma, hbl, ustar, bfsfc, m_or_s, wm, ws)
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
  implicit none
! !INPUT PARAMETERS:
  integer, intent(in) :: &
    m_or_s              ! flag =1 for wm only, 2 for ws, 3 for both

  real(r8), dimension(nx,ny), intent(in) :: &
    sigma,             &! normalized depth (d/hbl)
    hbl,               &! boundary layer depth
    bfsfc,             &! surface buoyancy forcing
    ustar               ! surface friction velocity

! !OUTPUT PARAMETERS:
  real(r8), dimension(nx,ny), intent(out) :: &
    wm,                &! turb velocity scales: momentum
    ws                  ! turb velocity scales: tracer

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------
  integer :: i,j  ! dummy loop indices

  real(r8), dimension(nx,ny) :: &
    zeta,           &! d/L or sigma*hbl/L(monin-obk)
    zetah            ! sigma*hbl*vonkar*BFSFC or ZETA = ZETAH/USTAR**3

!-----------------------------------------------------------------------
!  compute zetah and zeta - surface layer is special case
!-----------------------------------------------------------------------
  zetah = sigma*hbl*vonkar*bfsfc
  zeta  = zetah/(ustar**3 + eps)

!-----------------------------------------------------------------------
!  compute velocity scales for momentum
!-----------------------------------------------------------------------
  if(m_or_s == 1 .or. m_or_s == 3) then
    do j = 1, ny
      do i = 1, nx
        if(zeta(i,j) >= 0.d0) then ! stable region
          wm(i,j) = vonkar*ustar(i,j)/(1.0d0 + 5.0d0*zeta(i,j))
        else if (zeta(i,j) >= zeta_m) then
          wm(i,j) = vonkar*ustar(i,j)*(1.0d0 - 16.0d0*zeta(i,j))**0.25d0
        else
          wm(i,j) = vonkar*(a_m*(ustar(i,j)**3)-c_m*zetah(i,j))**0.33d0
        end if
      end do
    end do
  end if
!-----------------------------------------------------------------------
!  compute velocity scales for tracers
!-----------------------------------------------------------------------
  if(m_or_s == 2 .or. m_or_s == 3) then
    do j = 1, ny
      do i = 1, nx
        if(zeta(i,j) >= 0.d0) then
          ws(i,j) = vonkar*ustar(i,j)/(1.0d0 + 5.0d0*zeta(i,j))
        else if (zeta(i,j) >= zeta_s) then
          ws(i,j) = vonkar*ustar(i,j)*dsqrt(1.0d0 - 16.0d0*zeta(i,j))
        else
          ws(i,j) = vonkar*(a_s*(ustar(i,j)**3)-c_s*zetah(i,j))**0.33d0
        end if
      end do
    end do
  end if

end subroutine wscale

subroutine smooth_hblt(overwrite_hblt, use_hmxl, &
                       hblt, kbl, smooth_out)
! !DESCRIPTION:
!  This subroutine uses a 1-1-4-1-1 Laplacian filter one time
!  on HBLT or HMXL to reduce any horizontal two-grid-point noise.
!  If HBLT is overwritten, KBL is adjusted after smoothing.
  use tai_hyperlink, only: m_comm_cart, r8type2d, nbid, ndim, symm_np
  implicit none
! !INPUT PARAMETERS:
  logical, intent(in) :: &
    overwrite_hblt,   &    ! if .true.,  HBLT is overwritten
                             ! if .false., the result is returned in
                             !  a dummy array
    use_hmxl               ! if .true., smooth HMXL
                             ! if .false., smooth HBLT
! !INPUT/OUTPUT PARAMETERS:

  real(r8), dimension(nx,ny), optional, intent(inout) :: &
    hblt                   ! boundary layer depth

  integer, dimension(nx,ny), optional, intent(inout) :: &
    kbl                    ! index of first lvl below hbl

! !OUTPUT PARAMETERS:

  real(r8), dimension(nx,ny), optional, intent(out) ::  &
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


!-----------------------------------------------------------------------
!     perform one smoothing pass since we cannot do the necessary
!     boundary updates for multiple passes.
!-----------------------------------------------------------------------
  if(use_hmxl) then
    work2 = hmxl
  else
    work2 = hblt
  endif
  work1 = 0.d0
  work1(1:nx,1:ny) = work2
 
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, work1, symm_np) 
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

  do k = 1, nz
    do j = 1, ny
      do i= 1, nx
        ztmp = -zgrid(k)
        if(k == kb(i,j) .and. work2(i,j) > ztmp ) then
          work2(i,j) = ztmp
        end if
      end do
    end do
  end do

  if(overwrite_hblt .and. .not.use_hmxl ) then
    hblt = work2
    do k = 1, nz
      do j = 1, ny
        do i= 1, nx
          ztmp = -zgrid(k)
          if(in(i,j,1)*kb(i,j) /= 0  .and. &
            (hblt(i,j)> -zgrid(k-1)) .and. &
            (hblt(i,j)<= ztmp      )) kbl(i,j) = k

        end do
      end do
    end do
  else
    smooth_out = work2
  end if

end subroutine smooth_hblt

subroutine blmix(visc, vdc, hblt, ustar, bfsfc, stable, &
                  kbl, ghat)
! !DESCRIPTION:
!  This routine computes mixing coefficients within boundary layer
!  which depend on surface forcing and the magnitude and gradient
!  of interior mixing below the boundary layer (matching).  These
!  quantities have been computed in other routines.
!
!  Caution: if mixing bottoms out at hbl = -zgrid(km) then
!  fictitious layer at km+1 is needed with small but finite width
!  hwide(km+1).
  implicit none
! !INPUT/OUTPUT PARAMETERS:

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
    i,j                 ! horizontal indices


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

!-----------------------------------------------------------------------
!  compute velocity scales at hbl
!-----------------------------------------------------------------------
  sigma = epssfc

  call wscale(sigma, hblt, ustar, bfsfc, 3, wm, ws)

!-----------------------------------------------------------------------
!  determine caseA = 0 if closer to KBL than KBL-1
!  KN is then the closest klevel to HBLT
!----------------------------------------------------------------------
  do j = 1, ny
    do i = 1, nx
      k = kbl(i,j)
      casea(i,j)  = .5d0 + SIGN(.5d0, -zgrid(k)-.5d0*hwide(k)-hblt(i,j))
    end do
  end do

  kn = NINT(casea)*(kbl-1) + (1-NINT(casea))*kbl

!-------------------------------------------------------------------
!  find the interior viscosities and derivatives at hbl by
!  interpolating derivative values at vertical interfaces.  compute
!  matching conditions for shape function.
!-----------------------------------------------------------------------
  f1 = stable*5.d0*bfsfc/(ustar**4+eps)

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if(k == kn(i,j)) then

          DELHAT(i,j) = .5d0*hwide(k) - zgrid(k) - HBLT(i,j)
          R     (i,j) = 1.d0 - DELHAT(i,j) / hwide(k)

          DVDZUP(i,j) = (VISC(i,j,k-1) - VISC(i,j,k  ))/hwide(k)
          DVDZDN(i,j) = (VISC(i,j,k  ) - VISC(i,j,k+1))/hwide(k+1)
          VISCP (i,j) = .5d0*( (1.d0-R(i,j))* &
                               (DVDZUP(i,j) + dabs(DVDZUP(i,j))) + &
                                   R(i,j) * &
                               (DVDZDN(i,j) + dabs(DVDZDN(i,j))) )

          DVDZUP(i,j) = (VDC(i,j,k-1,2) - VDC(i,j,k  ,2))/hwide(k)
          DVDZDN(i,j) = (VDC(i,j,k  ,2) - VDC(i,j,k+1,2))/hwide(k+1)
          DIFSP (i,j) = .5d0*( (1.d0-R(i,j))* &
                               (DVDZUP(i,j) + dabs(DVDZUP(i,j))) + &
                                   R(i,j) * &
                               (DVDZDN(i,j) + dabs(DVDZDN(i,j))) )

          DVDZUP(i,j) = (VDC(i,j,k-1,1) - VDC(i,j,k  ,1))/hwide(k)
          DVDZDN(i,j) = (VDC(i,j,k  ,1) - VDC(i,j,k+1,1))/hwide(k+1)
          DIFTP (i,j) = .5d0*( (1.d0-R(i,j))* &
                               (DVDZUP(i,j) + dabs(DVDZUP(i,j))) + &
                                   R(i,j) * &
                               (DVDZDN(i,j) + dabs(DVDZDN(i,j))) )
          VISCH(i,j) = VISC(i,j,k)  + VISCP(i,j)*DELHAT(i,j)
          DIFSH(i,j) = VDC(i,j,k,2) + DIFSP(i,j)*DELHAT(i,j)
          DIFTH(i,j) = VDC(i,j,k,1) + DIFTP(i,j)*DELHAT(i,j)

          GAT1(i,j,1) = VISCH(i,j) / HBLT(i,j) /(WM(i,j)+eps)
          DAT1(i,j,1) = -VISCP(i,j)/(WM(i,j)+eps) + &
                           F1(i,j)*VISCH(i,j)

          GAT1(i,j,2) = DIFSH(i,j) / HBLT(i,j) /(WS(i,j)+eps)
          DAT1(i,j,2) = -DIFSP(i,j)/(WS(i,j)+eps) + &
                           F1(i,j)*DIFSH(i,j)

          GAT1(i,j,3) = DIFTH(i,j) / HBLT(i,j) /(WS(i,j)+eps)
          DAT1(i,j,3) = -DIFTP(i,j)/(WS(i,j)+eps) + &
                           F1(i,j)*DIFTH(i,j)

        endif
      end do
    end do
  end do

  dat1 = min(dat1,0.d0)

!-----------------------------------------------------------------------
!  compute the dimensionless shape functions and diffusivities
!  at the grid interfaces.  also compute function for non-local
!  transport term (GHAT).
!-----------------------------------------------------------------------
  do k = 1, nz
    sigma = (-zgrid(k) + 0.5d0*hwide(k))/hblt
    f1 = min(sigma,epssfc)

    call wscale(f1, hblt, ustar, bfsfc, 3, wm, ws)

    do j = 1, ny
      do i = 1, nx
        BLMC(i,j,k,1) = HBLT(i,j)*WM(i,j)*SIGMA(i,j)*       &
                         (1.d0 + SIGMA(i,j)*((SIGMA(i,j)-2.d0) + &
                         (3.d0-2.d0*SIGMA(i,j))*GAT1(i,j,1) +    &
                         (SIGMA(i,j)-1.d0)*DAT1(i,j,1)))

        BLMC(i,j,k,2) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*       &
                         (1.d0 + SIGMA(i,j)*((SIGMA(i,j)-2.d0) + &
                         (3.d0-2.d0*SIGMA(i,j))*GAT1(i,j,2) +    &
                         (SIGMA(i,j)-1.d0)*DAT1(i,j,2)))

        BLMC(i,j,k,3) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*       &
                         (1.d0 + SIGMA(i,j)*((SIGMA(i,j)-2.d0) + &
                         (3.d0-2.d0*SIGMA(i,j))*GAT1(i,j,3) +    &
                         (SIGMA(i,j)-1.d0)*DAT1(i,j,3)))
        GHAT(i,j,k) = (1.d0-STABLE(i,j))* cg/(WS(i,j)*HBLT(i,j) +eps)
      end do
    end do
  end do

!-----------------------------------------------------------------------
!  find diffusivities at kbl-1 grid level
!-----------------------------------------------------------------------
  do j = 1, ny
    do i = 1, nx
      k = kbl(i,j) - 1
      sigma(i,j) = -zgrid(k)/hblt(i,j)
    end do
  end do

  f1 = min(sigma,epssfc)
  call wscale(f1, hblt, ustar, bfsfc, 3, wm, ws)

  do j = 1, ny
    do i = 1, nx

      DKM1(i,j,1) = HBLT(i,j)*WM(i,j)*SIGMA(i,j)*     &
                    (1.d0+SIGMA(i,j)*((SIGMA(i,j)-2.d0) + &
                    (3.d0-2.d0*SIGMA(i,j))*GAT1(i,j,1)  + &
                    (SIGMA(i,j)-1.d0)*DAT1(i,j,1)))

      DKM1(i,j,2) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*     &
                    (1.d0+SIGMA(i,j)*((SIGMA(i,j)-2.d0) + &
                    (3.d0-2.d0*SIGMA(i,j))*GAT1(i,j,2)  + &
                    (SIGMA(i,j)-1.d0)*DAT1(i,j,2)))

      DKM1(i,j,3) = HBLT(i,j)*WS(i,j)*SIGMA(i,j)*     &
                    (1.d0+SIGMA(i,j)*((SIGMA(i,j)-2.d0) + &
                    (3.d0-2.d0*SIGMA(i,j))*GAT1(i,j,3)  + &
                    (SIGMA(i,j)-1.d0)*DAT1(i,j,3)))

    end do
  end do

!-----------------------------------------------------------------------
!  compute the enhanced mixing
!-----------------------------------------------------------------------
  do k = 1, nz-1
    where (k == (KBL - 1)) DELHAT = (HBLT + zgrid(k))/(zgrid(k)-zgrid(k+1))
    do j = 1, ny
      do i = 1, nx
        if(k == (kbl(i,j) - 1)) then

          BLMC(i,j,k,1) = (1.d0-DELHAT(i,j))*VISC(i,j,k) +           &
                        DELHAT(i,j) *(                                 &
                        (1.d0-DELHAT(i,j))**2*DKM1(i,j,1) +            &
                        DELHAT(i,j)**2*(CASEA(i,j)*VISC(i,j,k)+        &
                        (1.d0-CASEA(i,j))*BLMC(i,j,k,1)))

          BLMC(i,j,k,2) = (1.d0-DELHAT(i,j))*VDC(i,j,k,2) +          &
                        DELHAT(i,j)*(                                  &
                        (1.d0-DELHAT(i,j))**2*DKM1(i,j,2) +            &
                        DELHAT(i,j)**2*(CASEA(i,j)*VDC(i,j,k,2)+       &
                        (1.d0-CASEA(i,j))*BLMC(i,j,k,2)))

          BLMC(i,j,k,3) = (1.d0-DELHAT(i,j))*VDC(i,j,k,1) +          &
                        DELHAT(i,j) *(                                 &
                        (1.d0-DELHAT(i,j))**2*DKM1(i,j,3) +            &
                        DELHAT(i,j)**2*(CASEA(i,j)*VDC(i,j,k,1)+       &
                        (1.d0-CASEA(i,j))*BLMC(i,j,k,3)))

          GHAT(i,j,k) = (1.d0-CASEA(i,j)) * GHAT(i,j,k)

        end if
      end do
    end do
  end do

!-----------------------------------------------------------------------
!  combine interior and boundary layer coefficients and nonlocal term
!-----------------------------------------------------------------------
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
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

end subroutine blmix

end module tai_timcom_kpp
