      module rphyscons
      use machine,only:kind_phys
!  physical constants as set in nmc handbook from smithsonian tables.
!  physical constants are given to 5 places.
!  1990/04/30: g and rd are made consistent with nws usage.
!  2001/10/22: g made consistent with si usage.
!  2009/12/28: add some microphysics constants by zyf
!  math constants
      real(kind=kind_phys),parameter:: con_pi      =3.1415926535897931 ! pi
      real(kind=kind_phys),parameter:: con_sqrt2   =1.414214e+0 ! square root of 2
      real(kind=kind_phys),parameter:: con_sqrt3   =1.732051e+0 ! square root of 3
!  primary constants
      real(kind=kind_phys),parameter:: con_rerth   =6.3712e+6 ! radius of earth     (m)
      real(kind=kind_phys),parameter:: con_g       =9.80665e+0! gravity             (m/s2)
      real(kind=kind_phys),parameter:: con_omega   =7.2921e-5 ! ang vel of earth    (1/s)
!> std atms pressure (pa)
       real(kind=kind_phys),parameter:: con_p0     =1.01325e5 
!      real(kind=kind_phys),parameter:: con_solr   =1.36822e+3     ! solar constant    (W/m2)-aer(2001)
!> solar constant (\f$W/m^{2}\f$)-liu(2002)
       real(kind=kind_phys),parameter:: con_solr_old =1.3660e+3  
!> solar constant (\f$W/m^{2}\f$)-nasa-sorce Tim(2008)
       real(kind=kind_phys),parameter:: con_solr   =1.3608e+3    
!      real(kind=kind_phys),parameter:: con_solr   =1.36742732e+3  ! solar constant    (W/m2)-gfdl(1989) - OPR as of Jan 2006

!> molar gas constant (\f$J/mol/K\f$)
      real(kind=kind_phys),parameter:: con_rgas   =8.314472      
      real(kind=kind_phys),parameter:: con_rd      =2.8705e+2 ! gas constant air    (j/kg/k)
      real(kind=kind_phys),parameter:: con_rv      =4.6150e+2 ! gas constant h2o    (j/kg/k)
      real(kind=kind_phys),parameter:: con_cp      =1.0046e+3 ! spec heat air @p    (j/kg/k)
      real(kind=kind_phys),parameter:: con_cv      =7.1760e+2 ! spec heat air @v    (j/kg/k)
      real(kind=kind_phys),parameter:: con_cvap    =1.8460e+3 ! spec heat h2o gas   (j/kg/k)
      real(kind=kind_phys),parameter:: con_cliq    =4.1855e+3 ! spec heat h2o liq   (j/kg/k)
      real(kind=kind_phys),parameter:: con_csol    =2.1060e+3 ! spec heat h2o ice   (j/kg/k)
      real(kind=kind_phys),parameter:: con_hvap    =2.5000e+6 ! lat heat h2o cond   (j/kg)
      real(kind=kind_phys),parameter:: con_hfus    =3.3358e+5 ! lat heat h2o fusion (j/kg)
      real(kind=kind_phys),parameter:: con_psat    =6.1078e+2 ! pres at h2o 3pt     (pa)  
!     real(kind=kind_phys),parameter:: con_sbc     =5.6730e-8 ! stefan-boltzmann    (w/m2/k4)
!yj      real(kind=kind_phys),parameter:: con_solr    =1.3533e+3 ! solar constant      (w/m2)
      real(kind=kind_phys),parameter:: con_t0c     =2.7315e+2 ! temp at 0c          (k)
      real(kind=kind_phys),parameter:: con_ttp     =2.7316e+2 ! temp at h2o 3pt     (k)
!> temp freezing sea (K)
      real(kind=kind_phys),parameter:: con_tice   =2.7120e+2 
      real(kind=kind_phys),parameter:: con_jcal    =4.1855e+0 ! joules per calorie  ()
!> sea water reference density (\f$kg/m^{3}\f$)
      real(kind=kind_phys),parameter:: con_rhw0   =1022.0     
!> min q for computing precip type
      real(kind=kind_phys),parameter:: con_epsq   =1.0E-12    
!> min q for condensate mass mixing ratio
      real(kind=kind_phys),parameter:: con_qcmin  =1.0E-20
! add by zyf
      real(kind=kind_phys),parameter:: con_rhoair0 = 1.28      ! air density        (kg/m3)
      real(kind=kind_phys),parameter:: con_rhoh2o  = 1000.0    ! water density      (kg/m3)
      real(kind=kind_phys),parameter:: con_rhosnow = 100.0     ! snow density       (kg/m3)

!  secondary constants
      real(kind=kind_phys),parameter:: con_rocp    =con_rd/con_cp
      real(kind=kind_phys),parameter:: con_cpor    =con_cp/con_rd
      real(kind=kind_phys),parameter:: con_rog     =con_rd/con_g
      real(kind=kind_phys),parameter:: con_fvirt   =con_rv/con_rd-1.
      real(kind=kind_phys),parameter:: con_eps     =con_rd/con_rv
      real(kind=kind_phys),parameter:: con_epsm1   =con_rd/con_rv-1.
      real(kind=kind_phys),parameter:: con_dldt    =con_cvap-con_cliq
      real(kind=kind_phys),parameter:: con_xpona   =-con_dldt/con_rv
      real(kind=kind_phys),parameter:: con_xponb   =                    &
     &                       -con_dldt/con_rv+con_hvap/(con_rv*con_ttp)

!> \name Other Physics/Chemistry constants (source: 2002 CODATA)

!> speed of light (\f$m/s\f$)
      real(kind=kind_phys),parameter:: con_c      =2.99792458e+8  
!> planck constant (\f$J/s\f$)
      real(kind=kind_phys),parameter:: con_plnk   =6.6260693e-34 
!> boltzmann constant (\f$J/K\f$)
      real(kind=kind_phys),parameter:: con_boltz  =1.3806505e-23 
!> stefan-boltzmann (\f$W/m^{2}/K^{4}\f$)
       real(kind=kind_phys),parameter:: con_sbc    =5.670400e-8  
!> avogadro constant (\f$mol^{-1}\f$)
      real(kind=kind_phys),parameter:: con_avgd   =6.0221415e23 
!> vol of ideal gas at 273.15K, 101.325kPa (\f$m^{3}/mol\f$)
      real(kind=kind_phys),parameter:: con_gasv   =22413.996e-6  
!     real(kind=kind_phys),parameter:: con_amd    =28.970         ! molecular wght of dry air (g/mol)
!> molecular wght of dry air (\f$g/mol\f$)
      real(kind=kind_phys),parameter:: con_amd    =28.9644     
!> molecular wght of water vapor (\f$g/mol\f$)
      real(kind=kind_phys),parameter:: con_amw    =18.0154   
!> molecular wght of o3 (\f$g/mol\f$)
      real(kind=kind_phys),parameter:: con_amo3   =47.9982 
!     real(kind=kind_phys),parameter:: con_amo3   =48.0           ! molecular wght of o3  (g/mol)
!> molecular wght of co2 (\f$g/mol\f$)
      real(kind=kind_phys),parameter:: con_amco2  =44.011     
!> molecular wght of o2 (\f$g/mol\f$)
      real(kind=kind_phys),parameter:: con_amo2   =31.9999  
!> molecular wght of ch4 (\f$g/mol\f$)
      real(kind=kind_phys),parameter:: con_amch4  =16.043  
!> molecular wght of n2o (\f$g/mol\f$)
      real(kind=kind_phys),parameter:: con_amn2o  =44.013 
!> temperature the H.G.Nuc. ice starts
      real(kind=kind_phys), parameter:: con_thgni  =-38.15   

!> \name Miscellaneous physics related constants (Moorthi - Jul 2014)

! integer, parameter :: max_lon=16000, max_lat=8000, min_lon=192, min_lat=94
! integer, parameter :: max_lon=5000,  max_lat=2500, min_lon=192, min_lat=94 ! current opr
! integer, parameter :: max_lon=5000,  max_lat=2000, min_lon=192, min_lat=94 ! current opr
! integer, parameter :: max_lon=8000,  max_lat=4000, min_lon=192, min_lat=94 ! current opr
! real(kind=kind_phys), parameter:: rlapse  = 0.65e-2, rhc_max = 0.9999      ! current opr
! real(kind=kind_phys), parameter:: rlapse  = 0.65e-2, rhc_max = 0.9999999   ! new
! real(kind=kind_phys), parameter:: rlapse  = 0.65e-2, rhc_max = 0.9900

  real(kind=kind_phys), parameter:: rlapse  = 0.65e-2
! real(kind=kind_phys), parameter:: cb2mb   = 10.0, pa2mb   = 0.01
! for wsm6
  real(kind=kind_phys),parameter:: rhowater   = 1000.         ! density of water (kg/m^3)
! real(kind=kind_phys),parameter:: rhosnow    = 100.          ! density of snow (kg/m^3)
  real(kind=kind_phys),parameter:: rhoair     = 1.28          ! density of air near surface (kg/m^3)

! real(kind=kind_phys) :: dxmax, dxmin, dxinv, rhc_max

!........................................!
      end module rphyscons                !
!========================================!
!! @}
