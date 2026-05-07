module tai_timcom_phy
  use tai_timcom_const
  use tai_timcom_comm
  implicit none

  integer, parameter  ::    &
    opt_bulk_ncar = 1,      &! ncar_ocean_fluxes_mod
    opt_bulk_cesm = 2        ! cesmtimcom_shr_flux_mod

  integer, parameter  ::    &
    opt_bulk = opt_bulk_cesm   ! option for deriving turbulent flux

  real(kind=r8), parameter ::       &
    avgalb           = 0.066d0,     &
    shr_const_g      = 9.80616d0,   &
    shr_const_karman = 0.4d0,       &
    shr_const_boltz  = 1.38065d-23, &! Boltzmann's constant ~ J/K/molecule
    shr_const_avogad = 6.02214d26,  &! Avogadro's number ~ molecules/kmole
    shr_const_mwdair = 28.966d0,    &! molecular weight dry air ~ kg/kmole
    shr_const_mwwv   = 18.016d0,    &! molecular weight water vapor
    shr_const_cpdair = 1.00464d3,   &! specific heat of dry air   ~ J/kg/K
    shr_const_cpwv   = 1.810d3,     &! specific heat of water vap ~ J/kg/K
    shr_const_latvap = 2.501d6,     &! latent heat of evaporation ~ J/kg
    shr_const_latice = 3.337d5,     &! latent heat of fusion      ~ J/kg
    shr_const_stebol = 5.67d-8,     &! Stefan-Boltzmann constant ~ W/m^2/K^4
    shr_const_rgas   = shr_const_avogad*shr_const_boltz,     &! Universal gas constant ~ J/K/kmole
    shr_const_rdair  = shr_const_rgas/shr_const_mwdair,      &! Dry air gas constant     ~ J/K/kg
    shr_const_rwv    = shr_const_rgas/shr_const_mwwv,        &! Water vapor gas constant ~ J/K/kg
    shr_const_zvir   = (shr_const_rwv/shr_const_rdair)-1.d0, &! RWV/RDAIR - 1.0
    shr_const_cpvir  = (shr_const_cpwv/shr_const_cpdair)-1.d0 ! CPWV/CPDAIR - 1.0

contains

subroutine wind_surf(nx, ny, smft, dt, odz, ssu, ssv)
  implicit none

  integer,  intent(in) :: nx, ny
  real(r8), intent(in) :: smft(nx,ny,2)
  real(r8), intent(in) :: dt, odz
  real(r8), intent(inout) :: ssu(nx,ny), ssv(nx,ny)
  real(r8) :: tmp
! the velocity arrays u,v have a perimeter "ghost zone".
! thus, we set wind forcing only in the interior zones.
! taux,tauy are surface wind stress components.
! taux,tauy units are force per unit area (i.e., energy per unit volume).
! tmp=odz(1)/rho
! all units are cgs, so we use rho=1.
  tmp = dt*odz
  ssu = ssu + tmp*smft(:,:,1)
  ssv = ssv + tmp*smft(:,:,2)

end subroutine wind_surf

subroutine heat_surf(nx, ny, stf_heat, dt, odz, sst)
  implicit none

  integer,  intent(in) :: nx, ny
  real(r8), intent(in) :: stf_heat(nx,ny), dt, odz
  real(r8), intent(inout) :: sst(nx,ny)

! the temperature array t has a permiter "ghost zone".
! thus, we set atmospheri! heat exchange only in the interior zones.
! typical latitudinal range: -50 to 50 watts per square meter.
! 1 watt per square meter = 1 erg per second per equare cm.
! example: surface heating qdot = 50 watts per square meter
! qdot=50.
! tmp is conversion factor from watts per square meter
! to deg ! per time step in top model layer.
! tmp=dt*odz(1)/rho/cp
! all units are cgs, so we use rho=1          !dog -> 1.026
! cp for water is 2.5e4 ergs/gram/deg c.
! Cp for water is 3.996
!POP_CpSW                = 3.996e7_POPd0  ! specific heat salt water
!POP_CpSW                = SHR_CONST_CPSW*10000.0_POPd0 ! erg/g/K

  sst = sst + dt*odz*stf_heat

end subroutine heat_surf

subroutine heat_volume(nx, ny, nz, qsw_trans, qsw, dt, odz, t2)
  implicit none

  integer, intent(in) :: nx, ny, nz
  real(r8), intent(in) :: qsw_trans(nx,ny,nz), qsw(nx,ny)
  real(r8), intent(in) :: dt, odz(nz)
  real(r8), intent(inout) :: t2(nx,ny,nz)
  integer :: k
 
! tmp is conversion factor from watts per square meter
! to deg ! per time step in top model layer.
! tmp=dt*odz(1)/rho/cp
! all units are cgs, so we use rho=1          !dog -> 1.026
! cp for water is 2.5e4 ergs/gram/deg c.
! Cp for water is 3.996
  do k = 1, nz
    t2(:,:,k) = t2(:,:,k) + dt*odz(k)*qsw*qsw_trans(:,:,k)
  end do
end subroutine heat_volume

subroutine salt_surf(nx, ny, stf_mass, dt, odz, sss)
  implicit none

  integer,  intent(in) :: nx, ny
  real(r8), intent(in) :: stf_mass(nx,ny), dt, odz
  real(r8), intent(inout) :: sss(nx,ny)

! Salinity= Salinity + [ tmp * sdot ]     (sdot==evapo)
! tmp converts precipitation and evaporation rate (kg/m2/s) to psu
! hence tmp is DT/DZ /waterdensity * salinity
! tmp =        dt*odz/1000 * psu
! from kg/m**2/s -> msu cm/s -> psu cm/s : -34.7*1.e-4*1.e3
! EVAPO<0 (negative up); PREC>0 (positive down); - (negative EVAPO+PREC) -->
! SALINITY INCREASE
! unit of SALT is kg(psu)/m^2/s

  sss = sss + dt*odz*1000.d0*stf_mass

end subroutine salt_surf

subroutine surface_flux_atmocn(nsize, mask, ifrac, uvec_o, vvec_o, tvec_o, &
                               uvec_a, vvec_a, tvec_a, qvec_a, pvec_a,  &
                               taux, tauy, senh, lath, lwup, evap)
  implicit none
  integer, intent(in) :: &
    nsize

  integer(2), intent(in) :: &
    mask(nsize)

  real(r8), dimension(nsize), intent(in) :: &
    uvec_o, vvec_o, tvec_o,  &
    uvec_a, vvec_a, tvec_a, qvec_a, pvec_a, ifrac

  real(r8), dimension(nsize), intent(out) :: &
    taux, tauy, senh, lath, lwup, evap

  real(r8) :: &
    rhovec_a(nsize), zatm(nsize), work(nsize)

  zatm = 10.d0
  select case(opt_bulk)
  case(opt_bulk_ncar)
  case(opt_bulk_cesm)

    rhovec_a = 0.d0
    where(mask .eq. 1)
      rhovec_a = pvec_a/((1.d0+0.608d0*qvec_a)*287.04d0*tvec_a)
    end where

    call cesm_shr_flux_atmocn(nsize, mask, zatm, uvec_a, vvec_a, tvec_a, &
                              qvec_a, rhovec_a, uvec_o, vvec_o, tvec_o, &
                              taux, tauy, senh, lath, lwup, evap)
 
    work = mask*rhovec_a*1.63d-3*dsqrt(uvec_a**2 + vvec_a**2)
    taux = ifrac*work*uvec_a + (1.d0-ifrac)*taux
    tauy = ifrac*work*vvec_a + (1.d0-ifrac)*tauy

    taux = 10.d0*taux
    tauy = 10.d0*tauy
    senh = (1.d0-ifrac)*senh
    lath = (1.d0-ifrac)*lath
    lwup = (1.d0-ifrac)*lwup
    evap = (1.d0-ifrac)*evap
  end select 
end subroutine surface_flux_atmocn

subroutine cesm_shr_flux_atmocn(nMax, mask, zbot, ubot, vbot, thbot, &
                                qbot, rbot, us, vs, ts, & 
                                taux, tauy, senh, lath, lwup, evap)
  implicit none

  integer, intent(in) :: &
    nMax

  integer(2), intent(in) :: &
    mask(nMax)

  real(r8), dimension(nMax), intent(in) :: &
    zbot, ubot, vbot, thbot, &
    qbot, rbot, us, vs, ts

  real(r8), dimension(nMax), intent(out) :: &
    taux, tauy, senh, lath, lwup, evap

  real(r8), parameter :: &
    umin  =  0.5d0, & ! minimum wind speed       (m/s)
    zref  = 10.0d0, & ! reference height           (m)
    ztref =  2.0d0    ! reference height for air T (m)

  real(r8) :: &
    vmag,     & ! surface wind magnitude   (m/s)
    thvbot,   & ! virtual temperature      (K)
    ssq,      & ! sea surface humidity     (kg/kg)
    delt,     & ! potential T difference   (K)
    delq,     & ! humidity difference      (kg/kg)
    stable,   & ! stability factor
    rdn,      & ! sqrt of neutral exchange coeff (momentum)
    rhn,      & ! sqrt of neutral exchange coeff (heat)
    ren,      & ! sqrt of neutral exchange coeff (water)
    rd,       & ! sqrt of exchange coefficient (momentum)
    rh,       & ! sqrt of exchange coefficient (heat)
    re,       & ! sqrt of exchange coefficient (water)
    ustar,    & ! ustar
    qstar,    & ! qstar
    tstar,    & ! tstar
    hol,      & ! H (at zbot) over L
    xsq,      & ! ?
    xqq,      & ! ?
    psimh,    & ! stability function at zbot (momentum)
    psixh,    & ! stability function at zbot (heat and water)
    psix2,    & ! stability function at ztref reference height
    alz,      & ! ln(zbot/zref)
    al2,      & ! ln(zref/ztref)
    u10n,     & ! 10m neutral wind
    tau,      & ! stress at zbot
    cp,       & ! specific heat of moist air
    bn,       & ! exchange coef funct for interpolation
    bh,       & ! exchange coef funct for interpolation
    fac,      & ! vertical interpolation factor
    spval       ! local missing value

  integer :: n
! --- local functions --------------------------------
  real(r8) :: qsat   ! function: the saturation humididty of air (kg/m^3)
  real(r8) :: cdn    ! function: neutral drag coeff at 10m
  real(r8) :: psimhu ! function: unstable part of psimh
  real(r8) :: psixhu ! function: unstable part of psimx
  real(r8) :: Umps   ! dummy arg ~ wind velocity (m/s)
  real(r8) :: Tk     ! dummy arg ~ temperature (K)
  real(r8) :: xd     ! dummy arg ~ ?

  qsat(Tk)   = 640380.d0 / dexp(5107.4d0/Tk)
  cdn(Umps)  = 0.0027d0 / Umps + 0.000142d0 + 0.0000764d0 * Umps
  psimhu(xd) = dlog((1.d0+xd*(2.d0+xd))*(1.d0+xd*xd)/8.d0) - 2.d0*datan(xd) + 1.571d0
  psixhu(xd) = 2.d0*dlog((1.d0 + xd*xd)/2.d0)
! -------------------------------------------------------------------------------
! PURPOSE:
!   computes atm/ocn surface fluxes
!
! NOTES:
!   o all fluxes are positive downward
!   o net heat flux = net sw + lw up + lw down + sen + lat
!   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
!
! ASSUMPTIONS:
!   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
!   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
!                                 ctn = .0180 sqrt(cdn), stable
!   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
! ------------------------------------------------------------------------------
  al2 = dlog(zref/ztref)

  taux = 0.d0
  tauy = 0.d0
  senh = 0.d0
  lath = 0.d0
  lwup = 0.d0
  evap = 0.d0

  do n=1, nMax
    if(mask(n) /= 0) then
    !--- compute some needed quantities ---
      vmag   = max(umin, dsqrt((ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2))
      thvbot = thbot(n) * (1.d0 + shr_const_zvir*qbot(n)) ! virtual temp (K)
      ssq    = 0.98d0*qsat(ts(n))/rbot(n)                 ! sea surf hum (kg/kg)
      delt   = thbot(n) - ts(n)                           ! pot temp diff (K)
      delq   = qbot(n) - ssq                              ! spec hum dif (kg/kg)
      alz    = dlog(zbot(n)/zref)
      cp     = shr_const_cpdair*(1.d0 + shr_const_cpvir*ssq)

    !------------------------------------------------------------
    ! first estimate of Z/L and ustar, tstar and qstar
    !------------------------------------------------------------

    !--- neutral coefficients, z/L = 0.0 ---
      stable = 0.5d0 + dsign(0.5d0, delt)
      rdn    = dsqrt(cdn(vmag))
      rhn    = (1.d0-stable)*0.0327d0 + stable*0.018d0
      ren    = 0.0346d0

    !--- ustar, tstar, qstar ---
      ustar = rdn*vmag
      tstar = rhn*delt
      qstar = ren*delq

    !--- compute stability & evaluate all stability functions ---
      hol    = shr_const_karman*shr_const_g*zbot(n)*  &
               (tstar/thbot(n)+qstar/(1.d0/shr_const_zvir+qbot(n)))/ustar**2
      hol    = dsign( min(dabs(hol),10.d0), hol )
      stable = 0.5d0 + dsign(0.5d0, hol)
      xsq    = max(sqrt(dabs(1.d0 - 16.d0*hol)), 1.d0)
      xqq    = sqrt(xsq)
      psimh  = -5.d0*hol*stable + (1.d0-stable)*psimhu(xqq)
      psixh  = -5.d0*hol*stable + (1.d0-stable)*psixhu(xqq)

    !--- shift wind speed using old coefficient ---
      rd   = rdn / (1.d0 + rdn/shr_const_karman*(alz-psimh))
      u10n = vmag * rd / rdn

    !--- update transfer coeffs at 10m and neutral stability ---
      rdn = sqrt(cdn(u10n))
      ren = 0.0346d0
      rhn = (1.d0-stable)*0.0327d0 + stable*0.018d0

    !--- shift all coeffs to measurement height and stability ---
      rd = rdn/(1.d0 + rdn/shr_const_karman*(alz-psimh))
      rh = rhn/(1.d0 + rhn/shr_const_karman*(alz-psixh))
      re = ren/(1.d0 + ren/shr_const_karman*(alz-psixh))

    !--- update ustar, tstar, qstar using updated, shifted coeffs --
      ustar = rd * vmag
      tstar = rh * delt
      qstar = re * delq

    !------------------------------------------------------------
    ! iterate to converge on Z/L, ustar, tstar and qstar
    !------------------------------------------------------------
    !--- compute stability & evaluate all stability functions ---
      hol    = shr_const_karman*shr_const_g*zbot(n)* &
               (tstar/thbot(n)+qstar/(1.d0/shr_const_zvir+qbot(n)))/ustar**2
      hol    = dsign( min(dabs(hol),10.d0), hol )
      stable = .5d0 + dsign(.5d0 , hol)
      xsq    = max(sqrt(dabs(1.d0 - 16.d0*hol)) , 1.d0)
      xqq    = sqrt(xsq)
      psimh  = -5.d0*hol*stable + (1.d0-stable)*psimhu(xqq)
      psixh  = -5.d0*hol*stable + (1.d0-stable)*psixhu(xqq)

    !--- shift wind speed using old coeffs ---
      rd   = rdn / (1.d0 + rdn/shr_const_karman*(alz-psimh))
      u10n = vmag*rd/rdn

    !--- update transfer coeffs at 10m and neutral stability ---
      rdn = sqrt(cdn(u10n))
      ren = 0.0346d0
      rhn = (1.d0 - stable)*0.0327d0 + stable*0.018d0

    !--- shift all coeffs to measurement height and stability ---
      rd = rdn/(1.d0 + rdn/shr_const_karman*(alz-psimh))
      rh = rhn/(1.d0 + rhn/shr_const_karman*(alz-psixh))
      re = ren/(1.d0 + ren/shr_const_karman*(alz-psixh))

    !--- update ustar, tstar, qstar using updated, shifted coeffs ---
      ustar = rd*vmag
      tstar = rh*delt
      qstar = re*delq

    !------------------------------------------------------------
    ! compute the fluxes
    !------------------------------------------------------------
      tau = rbot(n)*ustar*ustar
    !--- momentum flux ---
    ! taux
      taux(n) = tau*(ubot(n)-us(n))/vmag
    ! tauy
      tauy(n) = tau*(vbot(n)-vs(n))/vmag
    !--- heat flux ---
    ! sensible heat
      senh(n) = cp*tau*tstar/ustar
    ! latent heat
      lath(n) = shr_const_latvap*tau*qstar/ustar
    ! long wave radiation
      lwup(n) = -shr_const_stebol*ts(n)**4
    !--- water flux ---
    ! evaporation
      evap(n) = lath(n)/shr_const_latvap
    end if
  end do

end subroutine cesm_shr_flux_atmocn

subroutine pp82(nx, ny, nz, w, u2, v2, rho, hbk, ev, hv)
! ----------------------------------------------------------------------
! Calculate gradient Ri based vertical mixing coefficents
! so called "eddy viscosity and diffusivity"
! as per Pacanowski and Philander(1982)

! In the loop below, EV,HV units are cm-cm/s
! but EV,HV normalization by DZ is done
  use tai_hyperlink, only: iw, odzw, dt
  implicit none

  integer :: nx, ny, nz
  real(r8), intent(in) :: &
    w(nx,ny,nz+1), &
    u2(nx,ny,nz),  &
    v2(nx,ny,nz),  &
    rho(nx,ny,nz), &
    hbk(nz)

  real(r8), intent(out) :: &
    ev(nx,ny,nz), &
    hv(nx,ny,nz,2)

  integer :: i, j, k, l
  real(r8) :: tmpw, hbk0, tmp, temp, Ri, evisc_local, emax
  real(r8), parameter :: rzmx = 50.d0, ri_max = -0.9d0, evisc = 5.d0

  ev = 0.d0
  hv = 0.d0
  do k = 1, nz-1
    l = k + 1
    tmpw = 1.d0/(rzmx*odzw(l))
    hbk0 = hbk(k)
    do j = 1, ny
      do i = 1, nx
        temp = tmpw*dabs(w(i,j,l)) !units of cm-cm/s
        Ri = grav*(rho(i,j,l)-rho(i,j,k))*odzw(l)/(odzw(l)**2 &
                 *(0.001d0+(u2(i,j,l)-u2(i,j,k))**2+(v2(i,j,l)-v2(i,j,k))**2))
        Ri = max(Ri_max, Ri)
        tmp = 1.d0/(1.d0 + Ri)
        temp = tmp*temp

        ! we add ODZW factor & apply explicit stability limit
        evisc_local = min(evisc*tmp**2, 100.d0)
        hv(i,j,k,1) = (evisc_local*tmp + temp + hbk0)
        ev(i,j,k)   = (evisc_local     + temp )!+ add(i,j,k)) ! ADD includes VBK0
        !if(myid.eq.5.and.i.eq.i1.and.j.eq.46) then
        !  write(*,*) "yc check ev:"
        !  write(*,"(i6,3f16.5)") k, rho(i1,46,k:l), Ri
        !end if
      end do
    end do
  end do

  do k = 1, nz-1
    l = k + 1
    emax = 0.2d0/(dt*odzw(l)**2)
    do j = 1, ny
      do i = 1, nx
        tmp = dble(iw(i,j,l))*odzw(l)
        hv(i,j,k,1) = tmp*min(emax, hv(i,j,k,1))
        ev(i,j,k)   = tmp*min(emax, ev(i,j,k))
      end do
    end do
  end do

  hv(:,:,:,2) = hv(:,:,:,1)*0.1d0
end subroutine pp82
end module tai_timcom_phy

