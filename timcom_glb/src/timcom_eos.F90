module timcom_eos
  use timcom_const
  use timcom_comm
  implicit none

contains

subroutine compute_eos_Wright(nx, ny, nz, msk, tp, ts, pp, rho)
  implicit none

  integer, intent(in) :: &
    nx, ny, nz

  integer(2), intent(in) :: &
    msk(nx,ny,nz)

  real(r8), intent(in) :: &
    tp(nx,ny,nz), &
    ts(nx,ny,nz), &
    pp(nx,ny,nz)

  real(r8), intent(out) :: &
    rho(nx,ny,nz)

  integer :: i, j, k
  real(r8) :: tmpd, sltd, pisd

  rho = 0.d0
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if(msk(i,j,k) .eq. 1) then
          tmpd = max(-1.8d0, tp(i,j,k))
          tmpd = min(40.d0, tmpd)
          sltd = min(41.d0, ts(i,j,k))
          sltd = max(0.d0, sltd)
          pisd = pp(i,j,k)*1.d-6
          rho(i,j,k) = eos_Wright(tmpd, sltd, pisd)
        end if
      end do
    end do
  end do
end subroutine compute_eos_Wright

subroutine compute_eos_mkcoef(nx, ny, nz, msk, tp, ts, pp, rho)
  implicit none

  integer, intent(in) :: &
    nx, ny, nz

  integer(2), intent(in) :: &
    msk(nx,ny,nz)

  real(r8), intent(in) :: &
    tp(nx,ny,nz), &
    ts(nx,ny,nz), &
    pp(nx,ny,nz)

  real(r8), intent(out) :: &
    rho(nx,ny,nz)

  integer :: i, j, k
  real(r8) :: tmpd, sltd, pisd

  rho = 0.d0
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if(msk(i,j,k) .eq. 1) then
          tmpd = max(-1.8d0, tp(i,j,k))
          tmpd = min(40.d0, tmpd)
          sltd = min(41.d0, ts(i,j,k))
          sltd = max(0.d0, sltd)
          pisd = pp(i,j,k)*1.d-6
          rho(i,j,k) = eos_mkcoef(tmpd, sltd, pisd)
        end if
      end do
    end do
  end do
end subroutine compute_eos_mkcoef

subroutine compute_eos_mkcoef_zgrid(nx, ny, nz, msk, tp, ts, z_grid, rho)
  implicit none

  integer, intent(in) :: &
    nx, ny, nz

  integer(2), intent(in) :: &
    msk(nx,ny,nz)

  real(r8), intent(in) :: &
    tp(nx,ny,nz), &
    ts(nx,ny,nz), &
    z_grid(nz)

  real(r8), intent(out) :: &
    rho(nx,ny,nz)

  integer :: i, j, k
  real(r8) :: tmpd, sltd, pisd

  rho = 0.d0
  do k = 1, nz
    pisd = z_grid(k)*1.d-3
    do j = 1, ny
      do i = 1, nx
        if(msk(i,j,k) .eq. 1) then
          tmpd = max(-1.8d0, tp(i,j,k))
          tmpd = min(40.d0, tmpd)
          sltd = min(41.d0, ts(i,j,k))
          sltd = max(0.d0, sltd)
          rho(i,j,k) = eos_mkcoef(tmpd, sltd, pisd)
        end if
      end do
    end do
  end do
end subroutine compute_eos_mkcoef_zgrid

subroutine compute_eos_JacMcd(nx, ny, nz, msk, tp, ts, depth, rho)
  implicit none

  integer, intent(in) :: &
    nx, ny, nz

  integer(2), intent(in) :: &
    msk(nx,ny,nz)

  real(r8), intent(in) :: &
    tp(nx,ny,nz), &
    ts(nx,ny,nz), &
    depth(nz)

  real(r8), intent(out) :: &
    rho(nx,ny,nz)

  integer :: i, j, k
  real(r8) :: tmpd, sltd

  rho = 0.d0

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        if(msk(i,j,k) .eq. 1) then
          tmpd = max(-1.8d0, tp(i,j,k))
          tmpd = min(40.d0, tmpd)
          sltd = min(41.d0, ts(i,j,k))
          sltd = max(0.d0, sltd)
          rho(i,j,k) = eos_JacMcd(tmpd, sltd, depth(k))
        end if
      end do
    end do
  end do
end subroutine compute_eos_JacMcd

function eos_Wright(tp, s, p)
  implicit none
! Dan Wright's full equation of state
! TP potential temperature [degrees Celcius]
! S  salinity [ppt, parts per thousand]
! P  pressure .001/cm depth change [bars, 10^5 Pascals]
! Wright_density [grams/cm^3] multiply by 1000 to get kg/m^3
  real(r8) :: tp, s, p
  real(r8) :: p01, p02, eos_Wright

  p02 = 1747.4508988d0 + TP*(11.51588d0 - 0.046331033d0*TP) &
                       - S*(3.85429655d0 + 0.01353985d0*TP)
  p01 = P + 5884.81703666d0 + &
        TP*(39.803732d0 + TP*(-0.3191477d0 + TP*0.0004291133d0)) &
        + 2.6126277d0*S

  eos_Wright = p01/(p02 + 0.7028423*p01)

end function eos_Wright

function eos_JacMcD(Tt, Ts, z_r)
  implicit none
! Computes density via Equation Of State (EOS) for seawater.
! If so prescribed, non-linear EOS of Jackett and McDougall (1995)
! is used.
!
! Tt potential temperature [deg Celsius].
! Ts salinity [PSU].
! Tz pressure/depth, [depth in meters and negative].
!
! K0, K1 and K2 are the pressure polynomial coefficients for secant
! bulk modulus, so that
!
!               bulk = K0 - K1 * z + K2 * z**2 ;
!
! while rho1 is sea-water density [kg/m^3] at standard pressure
! of 1 Atm, so that the density anomaly at in-sity pressure is
!
!               rho = rho1 / (1 + z / bulk) - 1000
!
! << This equation of state formulation has been derived by Jackett
!    and McDougall (1992), unpublished manuscript, CSIRO, Australia.
!    It computes in-situ density anomaly as a function of potential
!    temperature (Celsius) relative to the surface, salinity (PSU),
!    and depth (meters).  It assumes  no  pressure  variation along
!    geopotential  surfaces,  that  is,  depth  and  pressure  are
!    interchangeable. >>
!                                          John Wilkin, 29 July 92
!
! Check Values: T=3 C S=35.5 PSU Z=-5000 m rho=1050.3639165364
!
!
!  Further Information:
!  http://www.brest.ird.fr/Roms_tools/

  real(r8) :: Tt, Ts, z_r
  real(r8) :: sqrtTs, K0, K1, K2, K3, rho1
  real(r8) :: eos_JacMcd


  real(r8), parameter :: &
    A00 = 19092.56d0, A01 = 209.8925d0, A02 =-3.041638d0 ,A03 =-1.852732d-3, A04 =-1.361629d-5, &
    A10 = 104.4077d0, A11 =-6.500517d0, A12 = 0.155319d0 ,A13 = 2.326469d-4, &
    AS0 =-5.587545d0, AS1 = 0.7390729d0,AS2 =-1.909078d-2, &
    B00 = 4.721788d-1,B01 = 1.028859d-2,B02 =-2.512549d-4,B03 =-5.939910d-7, &
    B10 =-1.571896d-2,B11 =-2.598241d-4,B12 = 7.267926d-6,BS1 = 2.042967d-3, &
    E00 = 1.045941d-5,E01 =-5.782165d-10,E02= 1.296821d-7, &
    E10 =-2.595994d-7,E11 =-1.248266d-9,E12 =-3.508914d-9, &
    QR  =999.842594d0,Q01 = 6.793952d-2,Q02 =-9.095290d-3, &
    Q03 = 1.001685d-4,Q04 =-1.120083d-6,Q05 = 6.536332d-9, &
    Q10 = 0.824493d0, Q11 =-4.08990d-3, Q12 = 7.64380d-5, Q13 =-8.24670d-7, Q14 = 5.38750d-9, &
    QS0 =-5.72466d-3, QS1 = 1.02270d-4, QS2 =-1.65460d-6, Q20 = 4.8314d-4

  sqrtTs = sqrt(Ts);

  K0 = A00 + Tt*(A01 + Tt*(A02 + Tt*(A03 + Tt*A04))) &
           + Ts*(A10 + Tt*(A11 + Tt*(A12 + Tt*A13))  &
                 + sqrtTs*(AS0 + Tt*(AS1 + Tt*AS2)))

  K1 = B00 + Tt*(B01 + Tt*(B02 + Tt*B03)) &
           + Ts*(B10 + Tt*(B11 + Tt*B12) + sqrtTs*BS1)

  K2 = E00 + Tt*(E01 + Tt*E02) + Ts*(E10 + Tt*(E11 + Tt*E12))

  rho1 = QR + Tt*(Q01 + Tt*(Q02 + Tt*(Q03 + Tt*(Q04 + Tt*Q05)))) &
            + Ts*(Q10 + Tt*(Q11 + Tt*(Q12 + Tt*(Q13 + Tt*Q14))) &
                  + sqrtTs*(QS0 + Tt*(QS1 + Tt*QS2))+ Ts*Q20)

  eos_JacMcd = rho1/(1.d0 + 0.1d0*z_r/(K0 - z_r*(K1 - z_r*K2)))

  eos_JacMcd = eos_JacMcd*1.d-3 ! unit from kg/m^3 to g/cm^3
end function eos_JacMcd

function eos_mkcoef(TP,SS,PP)
  implicit none
! Calculates the density from local coefficients.
! TP    potential temperature [degrees, Celcius]
! SS     salinity              [psu, effectively ppt]
! PP    pressure              [bars, 10^5 Pascals]
! crho  local density when TP=S=0
! cTP   gradient of local density wrt TP [kg/m^3/C]
! cS    gradient of local density wrt S  [kg/m^3/ppt]
! cTPTP .5*second derivative of local density wrt TP [kg/m^3/C/C]
! cSTP  cross-derivative of local density wrt S and TP [kg/m^3/C/ppt]
! Operation count = 8
  real(r8) :: TP, SS, PP
  real(r8) :: crho,cTP,cS,cTPTP,cSTP
  real(r8) :: eos_mkcoef

  call make_local_coefficients(TP, SS, PP, crho, cTP, cS, cTPTP, cSTP)
  eos_mkcoef = crho + (cTP + cTPTP*tp + cSTP*ss)*tp + cS*ss

  eos_mkcoef = eos_mkcoef*1.d-3 + 1.d0
end function eos_mkcoef

subroutine make_local_coefficients(TP, SS, PP, crho, cTP, cS, cTPTP, cSTP)
  implicit none
! Calculates local coefficients
! TP    potential temperature [degrees, Celcius]
! SS    salinity              [psu, effectively ppt]
! PP    pressure              [bars, 10^5 Pascals]
! crho  local density when TP=S=0
! cTP   gradient of local density wrt TP
! cS    gradient of local density wrt S
! cTPTP .5*second derivative of local density wrt TP
! cSTP  cross-derivative of local density wrt S and TP [kg/m^3/C/ppt]
! USE UNESCO_density, ONLY: density_T_S_TT_TS, temperature, &
!                           potential_temperature_T, potential_temperature_TT, &
!                           potential_temperature_TS
  real(r8), intent(in) :: TP,SS,PP
  real(r8), intent(out) :: cTP,cS,cTPTP,cSTP,crho
  real(r8) dum,TT
! REAL*8 temperature,potential_temperature_T
! REAL*8 potential_temperature_TS,potential_temperature_TT
! find coeffs for rho = crho + cTP*(TP-TPr) + cS*(S-Sr)
!                            + cTPTP*(TP-TPr)**2 + cSTP*(S-Sr)*(TP-TPr)
  TT = temperature(TP,SS,PP)
  call density_T_S_TT_TS(TT,SS,PP,crho,cTP,cS,cTPTP,cSTP)
  dum = 1.d0/potential_temperature_T(TT,SS,PP)
  cTP = cTP*dum
  cSTP = dum*(cSTP - cTP*potential_temperature_TS(PP))
  cTPTP= 0.5d0*dum*dum*(cTPTP-cTP*potential_temperature_TT(TT,PP))
! modify coeffs for rho = crho + cTP*TP + cS*S +cTPTP*TP**2 + cSTP*S*TP
  crho = crho-(cTP-cTPTP*TP)*TP-cS*SS+cSTP*SS*TP
  cTP = cTP-2.d0*cTPTP*TP-cSTP*SS
  cS = cS - cSTP*TP

end subroutine make_local_coefficients

function temperature(TP, SS, PP)
  implicit none
! Calculates insitu temperature from potential temperature
! TP potential temperature [degrees Celcius]
! SS  salinity [ppt]
! PP  pressure [bars, 10^5 pascals]
!
! The UNESCO equation for potential_temperature is inverted iteratively.
!
! See also: potential_temperature, density

  real(r8) :: TP, SS, PP, TT!,potential_temperature
  real(r8) :: temperature
  integer(2) :: iterate
  integer :: ierr
  character(len=256) :: log_info

  TT = TP
  temperature = TT     ! first `guess'
! under relax to obtain temperature from potential temperatuer
  do iterate =1, 100
    TT = TT + 0.999d0*(TP-potential_temperature(TT,SS,PP))
!  T=T+(TP-potential_temperature(T,S,P))
    if(dabs(temperature-TT)<1.d-5) then
      exit
    end if
    temperature = TT
  end do

  if(iterate >= 100) then
    if(myid_timcom .eq. rootid) then
      write(log_info,*) &
        '[ERROR]: temperature did not converge in 100 iterations, '// &
        'Program stopped by FUNCTION temperature'
      call comm_write_log_info(fid_log, log_info)
    end if
    call comm_finalize
  end if

end function temperature

function potential_temperature(TT,SS,PP)
  implicit none
! Calculates the potential temperature from
! T  temperature  [degrees Celcius]
! S  salinity in [parts per thousand (ppt)]
! P  pressure in [bars (10^5 Pascals)]
! EXAMPLE:
!        t=potential_temperature(T=10,S=25,P=1000)
!        gives t=8.46785160000000
! References:
! Gill, A. E. (1981) Atmosphere-Ocean Dynamics. Academic Press
! Bryden, H. L. (1973) New polynomials for thermal expansion,
!        adiabatic temperature gradient and geopotential
!        temperature gradient of sea water. Deep-Sea Res., 20, 401-408.
  real(r8) :: TT,SS,PP
  real(r8) :: potential_temperature

  potential_temperature =  &
      TT-PP*(3.6504D-4 + TT*(8.3198D-5 - TT*(5.4065D-7-TT*4.0274D-9))) &
      -PP*(SS-35.d0)*(1.7439D-5 - TT*2.9778D-7)                        &
      -PP*PP*(8.9309D-7 - TT*(3.1628D-8 - TT*2.1987D-10))              &
      +4.1057D-9*(SS-35.d0)*PP*PP                                      &
      -PP*PP*PP*(-1.6056D-10 + TT*5.0484D-12)
end function potential_temperature

function potential_temperature_T(TT,SS,PP)
! Calculates the derivative of potential temperature wrt insitu temperature T
! TT  temperature in degrees Celcius
! SS  salinity in parts per thousand
! PP  pressure in bars (10^5 Pascals)
!
! Note: the derivative of T wrt potential temperature is
!       1/potential_temperature_T
!
! References:
! Gill, A. E. (1981) Atmosphere-Ocean Dynamics. Academic Press
! Bryden, H. L. (1973) New polynomials for thermal expansion,
!        adiabatic temperature gradient and geopotential
!        temperature gradient of sea water. Deep-Sea Res., 20, 401-408.
  real(r8) :: TT,SS,PP
  real(r8) :: potential_temperature_T

  potential_temperature_T=  &
      1.d0+PP*(-8.3198d-5+TT*(1.0813d-6-TT*1.20822d-8) &
      +(SS-35.d0)*2.9778d-7+PP*((3.1628d-8-TT*4.3974d-10)-PP*5.0484d-12))

end function potential_temperature_T

function potential_temperature_TT(TT,PP)
  implicit none
! Calculates the 2nd derivative of potential temperature wrt temperature T
! TT  temperature in degrees Celcius
! PP  pressure in bars (10^5 Pascals)
!
! References:
! Gill, A. E. (1981) Atmosphere-Ocean Dynamics. Academic Press
! Bryden, H. L. (1973) New polynomials for thermal expansion,
!        adiabatic temperature gradient and geopotential
!        temperature gradient of sea water. Deep-Sea Res., 20, 401-408.
  real(r8) :: TT,PP
  real(r8) :: potential_temperature_TT

  potential_temperature_TT = &
      PP*( 1.0813d-6 - TT*2.41644d-8 - PP*4.3974d-10 )
end function potential_temperature_TT

function potential_temperature_TS(PP)
  implicit none
! Calculates the second derivative of potential temperature wrt
! insitu temperature T and salinity S
! P  pressure in bars (10^5 Pascals)
!
! References:
! Gill, A. E. (1981) Atmosphere-Ocean Dynamics. Academic Press
! Bryden, H. L. (1973) New polynomials for thermal expansion,
!        adiabatic temperature gradient and geopotential
!        temperature gradient of sea water. Deep-Sea Res., 20, 401-408.

  real(r8) :: PP
  real(r8) :: potential_temperature_TS

  potential_temperature_TS = PP*2.9778d-7
end function potential_temperature_TS

subroutine density_T_S_TT_TS(TT,SS,PP,rho,drho_dT,drho_dS, &
       drho_dT_dT,drho_dT_dS)
  implicit none
! Uses equation (A3.1) of Gill 1982 Atmosphere Ocean Dynamics to
! calculate the density of fresh water and then uses (A3.2) to
! calculate the density of saline water at one standard atmosphere.
! Pressure effects then accounted for using (A3.3).
!
! The derivatives of density with respect to temperature T and
! salinity S are also calculated. The second derivative wrt T
! and the TS cross-derivative are calculated.
!
! INPUT VARIABLES:
! TT   temperature, [degrees Celcius]
! SS   salinity, [psu, which effectively equals ppt]
! PP   pressure, [bars, 10^5 Pascals]
!
! OUTPUT VARIABLES:
! rho         density  at T,S,P            [kg/m^3]
! drho_dT     derivative of density wrt T  [kg/m^3/C]
! drho_dS     derivative of density wrt S  [kg/m^3/ppt]
!
! OPTIONAL OUTPUT VARIABLES:
! drho_dT_dT  2nd derivative of density wrt T  [kg/m^3/C/C]
! drho_dT_dS  cross derivative of density wrt T and S  [kg/m^3/C/ppt]
!
! test values are rho(T=5, S=0, P=0)      =  999.96675 - 1000
!                 rho(T=5, S=35,P=0)      = 1027.67547 - 1000
!                 rho(T=25, S=35, P=1000) = 1062.53817 - 1000
!
! This algorithm is modified so it works well at single precision.
! This is done by extracting an offset density of 1000 kg/m^3
! USE global_parameters, ONLY: ks
  real(r8) ::  TT, SS, PP, rho, drho_dT, drho_dS, drho_dT_dT, drho_dT_dS
  real(r8) ::  K, dum, dKdT, dKdS, r, dKdTdT, rdT, rdS, dKdSdT

! calculate the density of fresh water at temperature TT
  rho = -0.1574060d0 + (6.793952d-2 - (9.095290d-3-(1.001685d-4 &
      - (1.120083d-6 - 6.536332d-9*TT)*TT)*TT)*TT)*TT

! calculate the derivative wrt T of density of fresh water at temperature T
  drho_dT =  6.793952d-2 - (1.819058d-2 - (3.005055d-4 &
      -(4.480332d-6 -3.268166d-8*TT)*TT)*TT)*TT

! calculate the 2nd derivative wrt T of density of fresh water at temperature T
!CW      IF(PRESENT(drho_dT_dT))
  drho_dT_dT = - (1.819058d-2 - (6.01011d-4 &
      - (1.3440996d-5 - 1.3072664d-7*TT)*TT)*TT)

  dum = DSQRT(SS)

! calculate the density of saline water at 1 atmosphere
  rho = rho + SS*(0.824493d0 +4.8314d-4*SS - (4.0899d-3 - (7.6438d-5 &
      - (8.2467d-7 - 5.3875d-9*TT)*TT)*TT)*TT)                         &
      + SS*dum*(-5.72466d-3 + (1.0227d-4 - 1.6546d-6*TT)*TT)

! calculate the derivative wrt T of density of saline water at 1 atmosphere
  drho_dT = drho_dT &
      +SS*(-4.0899d-3 +(1.52876d-4-(2.47401d-6 - 2.155d-8*TT)*TT)*TT) &
      +SS*dum*(1.0227d-4-3.3092d-6*TT)

! calculate the 2nd derivative wrt T of density of saline water at 1 atmosphere
!CW      IF(PRESENT(drho_dT_dT))
  drho_dT_dT = drho_dT_dT &
      +SS*(1.52876d-4 - (4.94802d-6-6.465d-8*TT)*TT)-SS*dum*3.3092d-6

! calculate the derivative wrt S of density of saline water at 1 atmosphere
  drho_dS = 0.824493d0-(4.0899d-3 - (7.6438d-5                    &
      -(8.2467d0-7-5.3875d-9*TT)*TT)*TT)*TT                         &
      +dum*(-8.58699d-3+(1.53405d-4-2.4819d-6*TT)*TT)+9.6628d-4*SS

! calculate the density cross-derivative wrt S and T of saline water at 1 atmos
!CW       IF(PRESENT(drho_dT_dS))
  drho_dT_dS = - (4.0899d-3 - (1.52876d-4 &
      -(2.47401d-6-2.155d-8*TT)*TT)*TT)+dum*(1.53405d-4-4.9638d-6*TT)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%  PRESSURE EFFECTS via secant bulk modulus  %%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! calculate the secant bulk modulus of pure water
  K = 19652.21d0 + TT*(148.4206d0 - TT*(2.327105d0 - TT*(1.360477d-2 &
      - TT*5.155288d-5)))

! calculate the derivative wrt T of the secant bulk modulus of pure water
  dKdT = 148.4206d0 - TT*(4.65421d0 - TT*(4.081431d-2 -TT*2.0621152d-4))

! calculate the 2nd derivative wrt T of the secant bulk modulus of pure water
  dKdTdT =  - (4.65421d0 - TT*(8.162862d-2 - TT*6.1863456d-4))

! calculate the bulk modulus at one atmosphere
  K = K+SS*(54.6746d0 - TT*(0.603459d0 -TT*(1.09987d-2 - TT*6.1670d-5))) &
      +SS*dum*(7.944d-2 + TT*(1.6483d-2 - TT*5.3009d-4))

! calculate derivative wrt T of the bulk modulus at one atmosphere
  dKdT = dKdT +SS*(- 0.603459d0 + TT*(2.19974d-2 - TT*1.8501d-4)) &
      +SS*dum*(1.6483d-2 - 1.06018d-3*TT)

! calculate 2nd derivative wrt T of the bulk modulus at one atmosphere
!CW      IF(PRESENT(drho_dT_dT))
  dKdTdT = dKdTdT +SS*(2.19974d-2 - TT*3.7002d-4 - dum*1.06018d-3)

! calculate the derivative wrt S of the bulk modulus at one atmosphere
  dKdS = 54.6746d0 - TT*(0.603459d0 - TT*(1.09987d-2 - TT*6.1670d-5)) &
      +dum*(1.1916d-1 + TT*(2.47245d-2 - TT*7.95135d-4))

! calculate the cross derivative wrt S and T
! of the bulk modulus at one atmosphere
!CW      IF(PRESENT(drho_dT_dS))
  dKdSdT = - (.603459d0 - TT*(2.19974E-2 - TT*1.8501E-4)) &
      + dum*(2.47245d-2 - TT*1.59027E-3)
! correct the bulk modulus for pressure p
  K =K+PP*((3.239908d0+TT*(1.43713d-3+TT*(1.16092d-4-TT*5.77905d-7))) &
      + SS*(2.2838d-3 - TT*(1.0981d-5+TT*1.6078d-6))+1.91075d-4*SS*dum  &
      + PP*(   (8.50935d-5 - TT*(6.12293d-6 - TT*5.2787d-8))            &
      + SS*(-9.9348d-7 + TT*(2.0816d-8 + TT*9.1697d-10))))

! correct the derivative wrt T of the bulk modulus for pressure p
  dKdT=dKdT + PP*((1.43713d-3 + TT*(2.32184d-4 - TT*1.733715d-6))   &
      + SS*(- 1.0981d-5 - 3.2156d-6*TT)                                 &
      + PP*(  (- 6.12293d-6 + 1.05574d-7*TT)                            &
      + SS*(2.0816d-8 + 1.83394d-9*TT)))

! correct the 2nd derivative wrt T of the bulk modulus for pressure p
!CW      IF (PRESENT(drho_dT_dT))
  dKdTdT = dKdTdT + PP*(  (2.32184d-4 - TT*3.46743d-6)  &
      - SS*3.2156d-6 + PP*(1.05574d-7+ SS*1.83394d-9))

! correct the derivative wrt S of the bulk modulus for pressure p
  dKdS = dKdS + PP*((2.2838d-3 - TT*(1.0981d-5 + TT*1.6078d-6)) &
      +2.866125d-4*dum+PP*(-9.9348d-7+TT*(2.0816d-8 + TT*9.1697d-10)))

! correct the cross derivative wrt S and T of the bulk modulus for pressure p
!CW      IF(PRESENT(drho_dT_dS))
  dKdSdT = dKdSdT - PP*(1.0981d-5 + TT*3.2156d-6 &
      - PP*(2.0816d-8 + TT*1.83394d-9))

! correct the density for pressure effects
  r = rho
  rdT = drho_dT
  rdS = drho_dS

  dum = 1.d0/(K-PP)
  rho = (K*rho+1000.d0*PP)*dum
  drho_dT = (K*drho_dT + (r - rho)*dKdT)*dum
  drho_dS = (K*drho_dS + (r - rho)*dKdS)*dum

!CW      IF(PRESENT(drho_dT_dT))
  drho_dT_dT = dum*(r*dKdTdT+2.d0*rdT*dKdT+K*drho_dT_dT &
      + rho*(2.d0*dum*dKdT*dKdT-dKdTdT)- 2.d0*dum*(r*dKdT+K*rdT)*dKdT)

!CW      IF(PRESENT(drho_dT_dS))
  drho_dT_dS = dum*( dKdS*rdT+K*drho_dT_dS+rdS*dKdT &
      +r*dKdSdT-rho*dKdSdT-drho_dS*dKdT -drho_dT*dKdS)

end subroutine density_T_S_TT_TS

subroutine compute_eos_stability(nx, ny, nz, msk, z_grid, tp, ts, rho)
  implicit none

  integer :: nx, ny, nz

  integer(2), intent(in) :: &
    msk(nx,ny,nz)

  real(r8), intent(in) :: &
    z_grid(nz)

  real(r8), intent(inout) :: &
    tp(nx,ny,nz), &
    ts(nx,ny,nz), &
    rho(nx,ny,nz)

  integer :: i, j, k, n
  real(r8) :: tmpd, sltd, pisd

  call compute_eos_JacMcd(nx, ny, nz, msk, tp, ts, -1.d-2*z_grid, rho)

  do k = 2, nz
    do j = 1, ny
      do i = 1, nx
        if(rho(i,j,k) < rho(i,j,k-1) .and. msk(i,j,k) .eq. 1) then
          tp(i,j,k) = tp(i,j,k-1)
          ts(i,j,k) = ts(i,j,k-1)
          tmpd = max(-1.8d0, tp(i,j,k))
          tmpd = min(40.d0, tmpd)
          sltd = min(41.d0, ts(i,j,k))
          sltd = max(0.d0, sltd)
          rho(i,j,k) = eos_JacMcd(tmpd, sltd, -1.d-2*z_grid(k))
        end if
      end do
    end do
  end do
end subroutine compute_eos_stability

end module timcom_eos

