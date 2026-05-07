MODULE mo_constants

  use machine, only: kind_phys

  IMPLICIT NONE

  ! ----------------------------------------------------------------
  !
  ! module *mo_constants* basic universal constants  and derived constants.
  !
  ! ----------------------------------------------------------------

  REAL(kind=kind_phys), PARAMETER :: dayl = 86400.  ! length of the day (in seconds).

  REAL(kind=kind_phys) :: api     ! 2.*arcsin(1.).
  REAL(kind=kind_phys) :: a       ! radius of the earth.
  REAL(kind=kind_phys) :: omega   ! solid rotation velocity of the earth.
  REAL(kind=kind_phys) :: g       ! gravity acceleration.
  REAL(kind=kind_phys) :: cpd     ! specific heat at constant pressure (dry air).
  REAL(kind=kind_phys) :: cpv     !            idem               (water vapour).
  REAL(kind=kind_phys) :: rd      ! gas constant for dry air.
  REAL(kind=kind_phys) :: rv      !    idem      for water vapour.
  REAL(kind=kind_phys) :: rcpd    ! rcpd=1./cpd.
  REAL(kind=kind_phys) :: vtmpc1  ! vtmpc1=rv/rd-1.
  REAL(kind=kind_phys) :: vtmpc2  ! vtmpc2=cpv/cpd-1.
  REAL(kind=kind_phys) :: rhoh2o  ! density of liquid water.
  REAL(kind=kind_phys) :: alv     ! latent heat for vaporisation.
  REAL(kind=kind_phys) :: als     ! latent heat for sublimation.
  REAL(kind=kind_phys) :: alf     ! latent heat for fusion.
  REAL(kind=kind_phys) :: clw     ! specific heat for liquid water.
  REAL(kind=kind_phys) :: tmelt   ! temperature of fusion of ice.
  REAL(kind=kind_phys) :: solc    ! solar constant.
  REAL(kind=kind_phys) :: stbo    ! stephan boltzmann constant.
  REAL(kind=kind_phys) :: yearl   ! length of the year (in days).

  ! constants used for computation of saturation mixing ratio
  !   over liquid water(*c_les*) or ice(*c_ies*).

  REAL(kind=kind_phys) :: c1es    ! 610.78
  REAL(kind=kind_phys) :: c2es    ! 1es*rd/rv
  REAL(kind=kind_phys) :: c3les   ! 17.269
  REAL(kind=kind_phys) :: c3ies   ! 21.875
  REAL(kind=kind_phys) :: c4les   ! 35.86
  REAL(kind=kind_phys) :: c4ies   !  7.66
  REAL(kind=kind_phys) :: c5les   ! c3les*(tmelt-c4les)
  REAL(kind=kind_phys) :: c5ies   ! c3ies*(tmelt-c4ies)
  REAL(kind=kind_phys) :: c5alvcp ! c5les*alv/cpd
  REAL(kind=kind_phys) :: c5alscp ! c5ies*als/cpd
  REAL(kind=kind_phys) :: alvdcp  ! alv/cpd
  REAL(kind=kind_phys) :: alsdcp  ! als/cpd
  REAL(kind=kind_phys) :: zrg
  REAL(kind=kind_phys) :: r5alvcp
  REAL(kind=kind_phys) :: r5alscp
  REAL(kind=kind_phys) :: ralvdcp 
  REAL(kind=kind_phys) :: ralsdcp
  REAL(kind=kind_phys) :: ralfdcp 
  REAL(kind=kind_phys) :: rtwat
  REAL(kind=kind_phys) :: rtber
  REAL(kind=kind_phys) :: rtice

CONTAINS

  SUBROUTINE inicon

    ! Description:
    ! Preset constants in mo_constants.
    !
    ! Method:
    !
    ! *inicon* is called from *setdyn*.
    !
    ! Authors:
    !
    ! M. Jarraud, ECMWF, December 1982, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! H.-S. Bauer, MPI, Jul 1998, changed
    ! A. Rhodin, MPI, Jan 1999, subroutine inicon put into module mo_constants
    !
    ! for more details see file AUTHORS
    !

    IMPLICIT NONE

    !  Intrinsic functions 
    INTRINSIC ASIN

    !  Executable statements 

    !-- 1. Preset constants

    api    = 2.*ASIN(1.)
    a      = 6371000.
    omega  = .7292E-4
    g      = 9.80665
    cpd    = 1005.46
    cpv    = 1869.46
    rd     = 287.05
    rv     = 461.51

    rcpd   = 1./cpd
    vtmpc1 = rv/rd - 1.
    vtmpc2 = cpv/cpd - 1.

    rhoh2o = 1000.
    alv    = 2.5008E6
    als    = 2.8345E6
    alf    = als - alv      
!    alf    = 3.50E5          !new Tiedtke (MPAS)

    clw    = 4186.84
    tmelt  = 273.16

    solc   = 1365.
    stbo   = 5.67E-8

    c1es    = 610.78
    c2es    = c1es*rd/rv
    c3les   = 17.269
    c3ies   = 21.875
    c4les   = 35.86
    c4ies   =  7.66
    c5les   = c3les*(tmelt-c4les)
    c5ies   = c3ies*(tmelt-c4ies)
    c5alvcp = c5les*alv/cpd
    c5alscp = c5ies*als/cpd
    alvdcp  = alv/cpd
    alsdcp  = als/cpd
!> xb110
    zrg = 1.0/g
    r5alvcp = c5les*alv*rcpd
    r5alscp = c5ies*als*rcpd
    ralvdcp=alv*rcpd
    ralsdcp=als*rcpd 
    ralfdcp=alf*rcpd
    rtwat=tmelt
    rtber=tmelt-5.
    rtice=tmelt-23.
!< xb110

  END SUBROUTINE inicon

END MODULE mo_constants
