module timcom_const
  implicit none

  integer, parameter ::             &
    r4 = selected_real_kind(6,37),  &
    r8 = selected_real_kind(15,307),&
    lfsrf = 1

  real(r8), parameter ::   &
    o12 = 1.d0/12.d0,      &
    o24 = 1.d0/24.d0,      & 
    grav = 980.d0,         &
    pi = 4.d0*datan(1.d0), & 
    r0 = 6.4d8,          & ! earth radius[cm], 6.371d8
    d2r = pi/180.d0,       & ! degree to radian
    r2d = 180.d0/pi,       & ! radian to degree
    tlz_depth = 5.5d5,     & ! maxium ocean depth[cm]
    mn2sec  = 60.d0,       & ! minute to second
    hr2mn   = 60.d0,       & ! hour to minute
    hr2sec  = 3600.d0,     & ! hour to second 
    day2hr  = 24.d0,       & ! day to hour
    day2mn  = 1440.d0,     & ! day to minute
    day2sec = 86400.d0,    & ! day to second
    eps     = 1.0d-10, &
    eps2    = 1.0d-20, &
    ppt_to_salt = 1.d-3,   &              ! from psu to msu. namely, ppt to g/g.
    salt_to_ppt = 1000.d0, &              ! from msu to psu
    sea_ice_salinity = 4.0d0, &           ! salinity of sea ice formed (psu)
    ocn_ref_salinity = 34.7d0,&           ! ocean reference salinity (psu)
    rho_sw = 4.1d0/3.996d0,   &           ! density of seawater in (g/cm^3)
    cp_sw = 3.996d7, &                    ! specific heat salt water
    hflux_factor = 1000.d0/(rho_sw*cp_sw)
end module timcom_const
