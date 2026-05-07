      subroutine lightning_ecn                                           &
                 (me,nxj,klon,klev,t,u,v,             &
                 ncld,levh,ncw,qall,p,pgeo,snow_flx,lndj,ft)

!reference:
!Lopez P., 2016: A lightning parameterization for the ECMWF Integrated Forecasting System. Mon. Wea. Rev., 144, 3057-3075.

!input-------------------------------------------
!     nxj:
!     klon:
!     klev:
!     t: temperature (K)
!     u: zonal wind (m/s)
!     v: meridional wind (m/s)
!     ncld:
!     levh:
!     ncw:
!     qall: all water species (kg/kg)
!     p: pressure field (Pa)
!     pgeo: geopotential (m^2s^-2)
!     snow_flx: snow flux from mp (kgm^-2s^-1), only from wsm3
!     lndj: land sea index (land=1;sea=2) 
!output------------------------------------------
!     ft: flash density (km^2day^-1)

      use physcons, only:vtmpc1 => con_fvirt, g => con_g,               &
                         rd => con_rd, eps1 => con_fvirt

      implicit none
!------------------------------------------------------------------------
!input:
      integer, intent (in) :: klev,klon,nxj,me,ncld,levh,ncw
      real,    intent (in) :: t(klon,klev), u(klon,klev), v(klon,klev)   &
                             ,p(klon,klev), qall(klon,levh)   &
                             ,pgeo(klon,klev), snow_flx(klon,klev)  
      integer, intent (in) :: lndj(klon)
!output:
      real, intent (out)   :: ft(klon) !lightning flash density (flashes km^-2 day^-1)
!------------------------------------------------------------------------
      integer jl,jk,ntcw,ntrw,ntiw,ntsw,ntgl
      real    qgraup(klon,klev),qsnow(klon,klev),pap(klon,klev),den(klon,klev)         &
             ,q(klon,klev)                               & !water vaper mixing ratio
             ,qv(klon,klev)                              & !specific humidity
             ,qc(klon,klev)                              & !cloud water mixing ratio
             ,qi(klon,klev)                              & !ice mixing ratio
             ,qr(klon,klev)                              & !rain mixing ratio
             ,qs(klon,klev)                              & !snow mixing ratio
             ,qg(klon,klev)                              & !graupel mixing ratio
             ,qcond(klon,klev)                             !condensate amount
      real    charg(klon)                  !charging rate
      real    rh(klev),p1d(klev),t1d(klev),qv1d(klev),hgt1d(klev),      &
              cape1d(klev),zlcl1d(klev)
      real    tmp,ncape,ncin,zlcl,plcl,lidx
      real,parameter::  vgraup= 3.0, &     ! typical fall speed for graupel (ms^-1)
                        vsnow = 0.5, &     ! typical fall speed for snow (ms^-1)
                        alpha = 32.4
      real    beta,zdz
!-------------------------------------------------------------------
      den    = 0.
      qgraup = 0.
      qsnow  = 0.
      charg  = 0.
      ft     = 0.
      cape1d = 0.
      zlcl1d = 0.

      ntcw=1+ncw
      ntrw=1+ncw+1*klev
      ntiw=1+ncw+2*klev
      ntsw=1+ncw+3*klev
      ntgl=1+ncw+4*klev

         qv(:,:) = qall(:,1:klev)
         qc(:,:) = qall(:,ntcw:ntcw+klev-1)    !qc=qc+qi in wsm3
         qr(:,:) = qall(:,ntrw:ntrw+klev-1) !qr=rain+snow in wsm3
         do jl = 1,nxj
          do jk = 1,klev
            if (t(jl,jk) .lt. 273.15) qr(jl,jk) = 0.
          end do
         end do
         qcond(:,:) = qc(:,:) + qr(:,:)    ! water species without snow and graupel

         if (ncld .eq. 5)then
           qr(:,:) = qall(:,ntrw:ntrw+klev-1)
           qi(:,:) = qall(:,ntiw:ntiw+klev-1)
           qs(:,:) = qall(:,ntsw:ntsw+klev-1)
           qg(:,:) = qall(:,ntgl:ntgl+klev-1)
           qcond(:,:) = qc(:,:) + qr(:,:) + qi(:,:)
         end if

      do jk = 1,klev
       do jl = 1,nxj

         pap(jl,jk) = p(jl,jk)
         q(jl,jk) = qv(jl,jk)/(1.-qv(jl,jk))                !convert to mixing ratio
         tmp = MAX(1.E-10, q(jl,jk))                        !qv should be positive
         den(jl,jk) = pap(jl,jk)/(rd*t(jl,jk)*(1.+eps1*tmp))  !air density 

         if (ncld .eq. 2)then           !for wsm3
           if (lndj(jl) .eq. 1) beta = 0.7    !over land
           if (lndj(jl) .ne. 1) beta = 0.45    !over ocean or ice
           qgraup(jl,jk) = max(0.0,beta*snow_flx(jl,jk)/(den(jl,jk)*vgraup)) !eq.1
           qsnow(jl,jk) = max(0.0,(1-beta)*snow_flx(jl,jk)/(den(jl,jk)*vsnow)) !eq.2
         else if (ncld .eq. 5)then     !for wsm6
           qgraup(jl,jk) = qg(jl,jk)
           qsnow(jl,jk) = qs(jl,jk)
         end if

        end do
       end do

       do jl = 1,nxj
        do jk = 1,klev
          if (t(jl,jk) .ge. 248.15 .and. t(jl,jk) .le. 273.15) then
            zdz = (pgeo(jl,jk+1)-pgeo(jl,jk))/g             
            charg(jl) = charg(jl) + qgraup(jl,jk)*                      &
                        (qcond(jl,jk) + qsnow(jl,jk))*den(jl,jk)*zdz        !eq.3
          end if
        end do
       end do
!
!from WRF V4.0 module_diag_functions.F

       do jl = 1,nxj
        do jk = 1,klev
          p1d(jk) = p(jl,jk)
          t1d(jk) = t(jl,jk)
          qv1d(jk) = qv(jl,jk)
          hgt1d(jk) = pgeo(jl,jk) 
        end do
        
        call calc_rh(klev,p1d,t1d,qv1d,rh) 
        call Buoyancy(klev,t1d,rh,p1d,hgt1d,1,ncape,ncin,zlcl,plcl,lidx,1) 
!        
        cape1d(jl) = ncape
        zlcl1d(jl) = max(0.0,zlcl)

        ft(jl) = alpha*charg(jl)*sqrt(cape1d(jl))*(min(zlcl1d(jl)/g*1.0e-3,1.8))**2      !eq.5

       end do

!      if (me .eq. 0)then
!        print*,"MaxMin RH = ",maxval(rh)," ",minval(rh)
!        print*,"MaxMin den = ",maxval(den)," ",minval(den)
!        print*,"MaxMin cape = ",maxval(cape1d)," ",minval(cape1d)
!        print*,"MaxMin zlcl = ",maxval(zlcl1d/g)," ",minval(zlcl1d/g)
!        print*,"MaxMin qgraup = ",maxval(qgraup)," ",minval(qgraup)
!        print*,"MaxMin qsnow = ",maxval(qsnow)," ",minval(qsnow)
!        print*,"MaxMin charg = ",maxval(charg)," ",minval(charg)
!        print*,"MaxMin flash = ",maxval(ft)," ",minval(ft)
!      end if

      return
      end subroutine lightning_ecn

  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  !~                                                                          ~!
  !~ Name:                                                                    ~!
  !~    Buoyancy                                                              ~!
  !~                                                                          ~!
  !~ Description:                                                             ~!
  !~    This function computes Convective Available Potential Energy (CAPE)   ~!
  !~    with inhibition as a result of water loading given the data required  ~!
  !~    to run up a sounding.                                                 ~!
  !~                                                                          ~!
  !~    Additionally, since we are running up a sounding anyways, this        ~!
  !~    function returns the height of the Level of Free Convection (LFC) and ~!
  !~    the pressure at the LFC.  That-a-ways, we don't have to run up a      ~!
  !~    sounding later, saving a relatively computationally expensive         ~!
  !~    routine.                                                              ~!
  !~                                                                          ~!
  !~ Usage:                                                                   ~!
  !~    ostat = Buoyancy ( tK, rh, p, hgt, sfc, CAPE, ZLFC, PLFC, parcel )    ~!
  !~                                                                          ~!
  !~ Where:                                                                   ~!
  !~                                                                          ~!
  !~    IN                                                                    ~!
  !~    --                                                                    ~!
  !~    tK   = Temperature ( K )                                              ~!
  !~    rh   = Relative Humidity ( % )                                        ~!
  !~    p    = Pressure ( Pa )                                                ~!
  !~    hgt  = Geopotential heights ( m )                                     ~!
  !~    sfc  = integer rank within submitted arrays that represents the       ~!
  !~           surface                                                        ~!
  !~                                                                          ~!
  !~    OUT                                                                   ~!
  !~    ---                                                                   ~!
  !~    ostat         INTEGER return status. Nonzero is bad.                  ~!
  !~    CAPE ( J/kg ) Convective Available Potential Energy                   ~!
  !~    ZLFC ( gpm )  Height at the LFC                                       ~!
  !~    PLFC ( Pa )   Pressure at the LFC                                     ~!
  !~                                                                          ~!
  !~    tK, rh, p, and hgt are all REAL arrays, arranged from lower levels    ~!
  !~    to higher levels.                                                     ~!
  !~                                                                          ~!
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
    subroutine Buoyancy ( nz, tk, rh, p, hgt, sfc, cape, cin, zlfc, plfc, lidx,  &
                        parcel)
  
      IMPLICIT NONE
  
      INTEGER, INTENT ( IN )  :: nz          !~ Number of vertical levels
      INTEGER, INTENT ( IN )  :: sfc         !~ Surface level in the profile
      REAL,    INTENT ( IN )  :: tk   ( nz ) !~ Temperature profile ( K )
      REAL,    INTENT ( IN )  :: rh   ( nz ) !~ Relative Humidity profile ( % )
      REAL,    INTENT ( IN )  :: p    ( nz ) !~ Pressure profile ( Pa )
      REAL,    INTENT ( IN )  :: hgt  ( nz ) !~ Height profile ( gpm )
      REAL,    INTENT ( OUT ) :: cape        !~ CAPE ( J kg^-1 )
      REAL,    INTENT ( OUT ) :: cin         !~ CIN  ( J kg^-1 )
      REAL,    INTENT ( OUT ) :: zlfc        !~ LFC Height ( gpm )
      REAL,    INTENT ( OUT ) :: plfc        !~ LFC Pressure ( Pa )
      REAL,    INTENT ( OUT ) :: lidx        !~ Lifted index
      INTEGER                 :: ostat       !~ Function return status
                                             !~ Nonzero is bad.

      INTEGER, INTENT ( IN  ) :: parcel      !~ Most Unstable = 1 (default)
                                             !~ Mean layer    = 2
                                             !~ Surface based = 3
  
      !~ Derived profile variables
      !  -------------------------
      REAL                    :: ws   ( nz ) !~ Saturation mixing ratio
      REAL                    :: w    ( nz ) !~ Mixing ratio
      REAL                    :: dTvK ( nz ) !~ Parcel / ambient Tv difference
      REAL                    :: buoy ( nz ) !~ Buoyancy
      REAL                    :: tlclK       !~ LCL temperature ( K )
      REAL                    :: plcl        !~ LCL pressure ( Pa )
      REAL                    :: nbuoy       !~ Negative buoyancy
      REAL                    :: pbuoy       !~ Positive buoyancy
  
      !~ Source parcel information
      !  -------------------------
      REAL                    :: srctK       !~ Source parcel temperature ( K )
      REAL                    :: srcrh       !~ Source parcel rh ( % )
      REAL                    :: srcws       !~ Source parcel sat. mixing ratio
      REAL                    :: srcw        !~ Source parcel mixing ratio
      REAL                    :: srcp        !~ Source parcel pressure ( Pa )
      REAL                    :: srctheta    !~ Source parcel theta ( K )
      REAL                    :: srcthetaeK  !~ Source parcel theta-e ( K )
      INTEGER                 :: srclev      !~ Level of the source parcel
      REAL                    :: spdiff      !~ Pressure difference
   
      !~ Parcel variables
      !  ----------------
      REAL                    :: ptK        !~ Parcel temperature ( K )
      REAL                    :: ptvK       !~ Parcel virtual temperature ( K )
      REAL                    :: tvK        !~ Ambient virtual temperature ( K )
      REAL                    :: pw         !~ Parcel mixing ratio
  
      !~ Other utility variables
      !  -----------------------
      INTEGER                 :: i, j, k    !~ Dummy iterator
      INTEGER                 :: lfclev     !~ Level of LFC
      INTEGER                 :: prcl       !~ Internal parcel type indicator
      INTEGER                 :: mlev       !~ Level for ML calculation
      INTEGER                 :: lyrcnt     !~ Number of layers in mean layer
      LOGICAL                 :: flag       !~ Dummy flag
      LOGICAL                 :: wflag      !~ Saturation flag
      REAL                    :: freeze     !~ Water loading multiplier
      REAL                    :: pdiff      !~ Pressure difference between levs 
      REAL                    :: pm, pu, pd !~ Middle, upper, lower pressures
      REAL                    :: lidxu      !~ Lifted index at upper level
      REAL                    :: lidxd      !~ Lifted index at lower level
  
      !~ Thermo / dynamical constants
      !  ----------------------------
      REAL                    :: Rd         !~ Dry gas constant
         PARAMETER ( Rd = 287.058 )         !~ J deg^-1 kg^-1
      REAL                    :: Cp         !~ Specific heat constant pressure
         PARAMETER ( Cp = 1004.67 )         !~ J deg^-1 kg^-1
      REAL                    :: g          !~ Acceleration due to gravity
         PARAMETER ( g  = 9.80665 )         !~ m s^-2
      REAL                    :: RUNDEF
         PARAMETER ( RUNDEF = -9.999E30 )

      REAL :: TLCL, The2T, VirtualTemperature, SaturationMixingRatio, Theta, &
              Thetae
  
      !~ Initialize variables
      !  --------------------
      ostat  = 0
      CAPE   = REAL ( 0 )
      CIN    = REAL ( 0 )
      ZLFC   = RUNDEF
      PLFC   = RUNDEF
  
      !~ Look for submitted parcel definition
      !~ 1 = Most unstable
      !~ 2 = Mean layer
      !~ 3 = Surface based
      !  -------------------------------------
      IF ( parcel > 3 .or. parcel < 1 ) THEN
         prcl = 1
      ELSE
         prcl =  parcel
      END IF
  
      !~ Initalize our parcel to be (sort of) surface based.  Because of
      !~ issues we've been observing in the WRF model, specifically with
      !~ excessive surface moisture values at the surface, using a true
      !~ surface based parcel is resulting a more unstable environment
      !~ than is actually occuring.  To address this, our surface parcel
      !~ is now going to be defined as the parcel between 25-50 hPa
      !~ above the surface. UPDATE - now that this routine is in WRF,
      !~ going to trust surface info. GAC 20140415
      !  ----------------------------------------------------------------
  
      !~ Compute mixing ratio values for the layer
      !  -----------------------------------------
      DO k = sfc, nz
        ws  ( k )   = SaturationMixingRatio ( tK(k), p(k) )
        w   ( k )   = ( rh(k)/100.0 ) * ws ( k )
      END DO
  
      srclev      = sfc
      srctK       = tK    ( sfc )
      srcrh       = rh    ( sfc )
      srcp        = p     ( sfc )
      srcws       = ws    ( sfc )
      srcw        = w     ( sfc )
      srctheta    = Theta ( tK(sfc), p(sfc)/100.0 )
   
      !~ Compute the profile mixing ratio.  If the parcel is the MU parcel,
      !~ define our parcel to be the most unstable parcel within the lowest
      !~ 180 mb.
      !  -------------------------------------------------------------------
      mlev = sfc + 1
      DO k = sfc + 1, nz
   
         !~ Identify the last layer within 100 hPa of the surface
         !  -----------------------------------------------------
         pdiff = ( p (sfc) - p (k) ) / REAL ( 100 )
         IF ( pdiff <= REAL (100) ) mlev = k

         !~ If we've made it past the lowest 180 hPa, exit the loop
         !  -------------------------------------------------------
         IF ( pdiff >= REAL (180) ) EXIT

         IF ( prcl == 1 ) THEN
            !IF ( (p(k) > 70000.0) .and. (w(k) > srcw) ) THEN
            IF ( (w(k) > srcw) ) THEN
               srctheta = Theta ( tK(k), p(k)/100.0 )
               srcw = w ( k )
               srclev  = k
               srctK   = tK ( k )
               srcrh   = rh ( k )
               srcp    = p  ( k )
            END IF
         END IF
   
      END DO
   
      !~ If we want the mean layer parcel, compute the mean values in the
      !~ lowest 100 hPa.
      !  ----------------------------------------------------------------
      lyrcnt =  mlev - sfc + 1
      IF ( prcl == 2 ) THEN
   
         srclev   = sfc
         srctK    = SUM ( tK (sfc:mlev) ) / REAL ( lyrcnt )
         srcw     = SUM ( w  (sfc:mlev) ) / REAL ( lyrcnt )
         srcrh    = SUM ( rh (sfc:mlev) ) / REAL ( lyrcnt )
         srcp     = SUM ( p  (sfc:mlev) ) / REAL ( lyrcnt )
         srctheta = Theta ( srctK, srcp/100. )
   
      END IF
   
      srcthetaeK = Thetae ( srctK, srcp/100.0, srcrh, srcw )
   
      !~ Calculate temperature and pressure of the LCL
      !  ---------------------------------------------
      tlclK = TLCL ( tK(srclev), rh(srclev) )
      plcl  = p(srclev) * ( (tlclK/tK(srclev))**(Cp/Rd) )
   
      !~ Now lift the parcel
      !  -------------------
   
      buoy  = REAL ( 0 )
      pw    = srcw
      wflag = .false.
      DO k  = srclev, nz
         IF ( p (k) <= plcl ) THEN
   
            !~ The first level after we pass the LCL, we're still going to
            !~ lift the parcel dry adiabatically, as we haven't added the
            !~ the required code to switch between the dry adiabatic and moist
            !~ adiabatic cooling.  Since the dry version results in a greater
            !~ temperature loss, doing that for the first step so we don't over
            !~ guesstimate the instability.
            !  ----------------------------------------------------------------
   
            IF ( wflag ) THEN
               flag  = .false.
   
               !~ Above the LCL, our parcel is now undergoing moist adiabatic
               !~ cooling.  Because of the latent heating being undergone as
               !~ the parcel rises above the LFC, must iterative solve for the
               !~ parcel temperature using equivalant potential temperature,
               !~ which is conserved during both dry adiabatic and
               !~ pseudoadiabatic displacements.
               !  --------------------------------------------------------------
               ptK   = The2T ( srcthetaeK, p(k), flag )
   
               !~ Calculate the parcel mixing ratio, which is now changing
               !~ as we condense moisture out of the parcel, and is equivalent
               !~ to the saturation mixing ratio, since we are, in theory, at
               !~ saturation.
               !  ------------------------------------------------------------
               pw = SaturationMixingRatio ( ptK, p(k) )
   
               !~ Now we can calculate the virtual temperature of the parcel
               !~ and the surrounding environment to assess the buoyancy.
               !  ----------------------------------------------------------
               ptvK  = VirtualTemperature ( ptK, pw )
               tvK   = VirtualTemperature ( tK (k), w (k) )
   
               !~ Modification to account for water loading
               !  -----------------------------------------
               freeze = 0.033 * ( 263.15 - pTvK )
               IF ( freeze > 1.0 ) freeze = 1.0
               IF ( freeze < 0.0 ) freeze = 0.0
   
               !~ Approximate how much of the water vapor has condensed out
               !~ of the parcel at this level
               !  ---------------------------------------------------------
               freeze = freeze * 333700.0 * ( srcw - pw ) / 1005.7
   
               pTvK = pTvK - pTvK * ( srcw - pw ) + freeze
               dTvK ( k ) = ptvK - tvK
               buoy ( k ) = g * ( dTvK ( k ) / tvK )
   
            ELSE
   
               !~ Since the theta remains constant whilst undergoing dry
               !~ adiabatic processes, can back out the parcel temperature
               !~ from potential temperature below the LCL
               !  --------------------------------------------------------
               ptK   = srctheta / ( 100000.0/p(k) )**(Rd/Cp)
   
               !~ Grab the parcel virtual temperture, can use the source
               !~ mixing ratio since we are undergoing dry adiabatic cooling
               !  ----------------------------------------------------------
               ptvK  = VirtualTemperature ( ptK, srcw )
   
               !~ Virtual temperature of the environment
               !  --------------------------------------
               tvK   = VirtualTemperature ( tK (k), w (k) )
   
               !~ Buoyancy at this level
               !  ----------------------
               dTvK ( k ) = ptvK - tvK
               buoy ( k ) = g * ( dtvK ( k ) / tvK )
   
               wflag = .true.
   
            END IF
   
         ELSE
   
            !~ Since the theta remains constant whilst undergoing dry
            !~ adiabatic processes, can back out the parcel temperature
            !~ from potential temperature below the LCL
            !  --------------------------------------------------------
            ptK   = srctheta / ( 100000.0/p(k) )**(Rd/Cp)
   
            !~ Grab the parcel virtual temperture, can use the source
            !~ mixing ratio since we are undergoing dry adiabatic cooling
            !  ----------------------------------------------------------
            ptvK  = VirtualTemperature ( ptK, srcw )
   
            !~ Virtual temperature of the environment
            !  --------------------------------------
            tvK   = VirtualTemperature ( tK (k), w (k) )
   
            !~ Buoyancy at this level
            !  ---------------------
            dTvK ( k ) = ptvK - tvK
            buoy ( k ) = g * ( dtvK ( k ) / tvK )
   
         END IF

         !~ Chirp
         !  -----
  !          WRITE ( *,'(I15,6F15.3)' )k,p(k)/100.,ptK,pw*1000.,ptvK,tvK,buoy(k)
   
      END DO
   
      !~ Add up the buoyancies, find the LFC
      !  -----------------------------------
      flag   = .false.
      lfclev = -1
      nbuoy  = REAL ( 0 )
      pbuoy = REAL ( 0 )
      DO k = sfc + 1, nz
         IF ( tK (k) < 253.15 ) EXIT
         CAPE = CAPE + MAX ( buoy (k), 0.0 ) * ( hgt (k) - hgt (k-1) )
         CIN  = CIN  + MIN ( buoy (k), 0.0 ) * ( hgt (k) - hgt (k-1) )
   
         !~ If we've already passed the LFC
         !  -------------------------------
         IF ( flag .and. buoy (k) > REAL (0) ) THEN
            pbuoy = pbuoy + buoy (k)
         END IF
   
         !~ We are buoyant now - passed the LFC
         !  -----------------------------------
         IF ( .not. flag .and. buoy (k) > REAL (0) .and. p (k) < plcl ) THEN
            flag = .true.
            pbuoy = pbuoy + buoy (k)
            lfclev = k
         END IF
   
         !~ If we think we've passed the LFC, but encounter a negative layer
         !~ start adding it up.
         !  ----------------------------------------------------------------
         IF ( flag .and. buoy (k) < REAL (0) ) THEN
            nbuoy = nbuoy + buoy (k)

            !~ If the accumulated negative buoyancy is greater than the
            !~ positive buoyancy, then we are capped off.  Got to go higher
            !~ to find the LFC. Reset positive and negative buoyancy summations
            !  ----------------------------------------------------------------
            IF ( ABS (nbuoy) > pbuoy ) THEN
               flag   = .false.
               nbuoy  = REAL ( 0 )
               pbuoy  = REAL ( 0 )
               lfclev = -1
            END IF
         END IF

      END DO

      !~ Calculate lifted index by interpolating difference between
      !~ parcel and ambient Tv to 500mb.
      !  ----------------------------------------------------------
      DO k = sfc + 1, nz

         pm = 50000.
         pu = p ( k )
         pd = p ( k - 1 )

         !~ If we're already above 500mb just set lifted index to 0.
         !~ --------------------------------------------------------
         IF ( pd .le. pm ) THEN
            lidx = 0.
            EXIT
   
         ELSEIF ( pu .le. pm .and. pd .gt. pm) THEN

            !~ Found trapping pressure: up, middle, down.
            !~ We are doing first order interpolation.  
            !  ------------------------------------------
            lidxu = -dTvK ( k ) * ( pu / 100000. ) ** (Rd/Cp)
            lidxd = -dTvK ( k-1 ) * ( pd / 100000. ) ** (Rd/Cp)
            lidx = ( lidxu * (pm-pd) + lidxd * (pu-pm) ) / (pu-pd)
            EXIT

         ENDIF

      END DO
   
      !~ Assuming the the LFC is at a pressure level for now
      !  ---------------------------------------------------
      IF ( lfclev > 0 ) THEN
         PLFC = p   ( lfclev )
         ZLFC = hgt ( lfclev )
      END IF
   
      IF ( PLFC /= PLFC .OR. PLFC < REAL (0) ) THEN
         PLFC = REAL ( -1 )
         ZLFC = REAL ( -1 )
      END IF
   
      IF ( CAPE /= CAPE ) cape = REAL ( 0 )
   
      IF ( CIN  /= CIN  ) cin  = REAL ( 0 )

      !~ Chirp
      !  -----
  !       WRITE ( *,* ) ' CAPE: ', cape, ' CIN:  ', cin
  !       WRITE ( *,* ) ' LFC:  ', ZLFC, ' PLFC: ', PLFC
  !       WRITE ( *,* ) ''
  !       WRITE ( *,* ) ' Exiting buoyancy.'
  !       WRITE ( *,* ) ' ==================================== '
  !       WRITE ( *,* ) ''
   
  end subroutine Buoyancy 
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  !~ 
  !~ Name:
  !~    calc_rh
  !~
  !~ Description:
  !~    This function calculates relative humidity given pressure(Pa), 
  !~    temperature(K), and specific humidity(kg/kg).
  !~ 
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  subroutine calc_rh (nz, p, t, qv, rh )
    
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: nz !Number of vertical levels 
    REAL, INTENT(IN)  :: p(nz), t(nz), qv(nz)
    REAL, INTENT(OUT) :: rh(nz)

    ! Local
    ! -----
    REAL, PARAMETER :: pq0=379.90516
    REAL, PARAMETER :: a2=17.2693882
    REAL, PARAMETER :: a3=273.16
    REAL, PARAMETER :: a4=35.86
    REAL, PARAMETER :: rhmin=1.
    REAL :: q(nz), qs(nz)
    INTEGER :: i,j,k
  
    ! Following algorithms adapted from WRFPOST
    ! May want to substitute with another later
    ! -----------------------------------------
     do k = 1,nz
      q(k)=qv(k)
      qs(k)=pq0/p(k)*exp(a2*(t(k)-a3)/(t(k)-a4))
      rh(k)=100.*q(k)/qs(k)
      IF (rh(k) .gt. 100.) THEN
        rh(k)=100.
      ELSE IF (rh(k) .lt. rhmin) THEN
        rh(k)=rhmin
      ENDIF
     end do

  end subroutine calc_rh
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  !~
  !~ Name:
  !~    SaturationMixingRatio
  !~
  !~ Description:
  !~    This function calculates saturation mixing ratio given the
  !~    temperature ( K ) and the ambient pressure ( Pa ).  Uses 
  !~    approximation of saturation vapor pressure.
  !~
  !~ References:
  !~    Bolton (1980), Monthly Weather Review, pg. 1047, Eq. 10
  !~
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  FUNCTION SaturationMixingRatio ( tK, p ) result ( ws )

    IMPLICIT NONE

    REAL, INTENT ( IN ) :: tK
    REAL, INTENT ( IN ) :: p
    REAL                :: ws

    REAL :: es

    es = 6.122 * exp ( (17.67*(tK-273.15))/ (tK-29.66) )
    ws = ( 0.622*es ) / ( (p/100.0)-es )

  END FUNCTION SaturationMixingRatio
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  !~
  !~ Name:
  !~    Theta
  !~
  !~ Description:
  !~    This function calculates potential temperature as defined by
  !~    Poisson's equation, given temperature and pressure ( hPa ).
  !~
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  FUNCTION Theta ( t, p )
  IMPLICIT NONE

     !~ Variable declaration
     !  --------------------
     REAL, INTENT ( IN ) :: t
     REAL, INTENT ( IN ) :: p
     REAL                :: theta

     REAL :: Rd ! Dry gas constant
     REAL :: Cp ! Specific heat of dry air at constant pressure
     REAL :: p0 ! Standard pressure ( 1000 hPa )
  
     Rd =  287.04
     Cp = 1004.67
     p0 = 1000.00

     !~ Poisson's equation
     !  ------------------
     theta = t * ( (p0/p)**(Rd/Cp) )
  
  END FUNCTION Theta



  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  !~
  !~ Name:
  !~    Thetae
  !~
  !~ Description:
  !~    This function returns equivalent potential temperature using the 
  !~    method described in Bolton 1980, Monthly Weather Review, equation 43.
  !~
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  FUNCTION Thetae ( tK, p, rh, mixr )
  IMPLICIT NONE

     !~ Variable Declarations
     !  ---------------------
     REAL :: tK        ! Temperature ( K )
     REAL :: p         ! Pressure ( hPa )
     REAL :: rh        ! Relative humidity
     REAL :: mixr      ! Mixing Ratio ( kg kg^-1)
     REAL :: te        ! Equivalent temperature ( K )
     REAL :: thetae    ! Equivalent potential temperature
  
     REAL, PARAMETER :: R  = 287.04         ! Universal gas constant (J/deg kg)
     REAL, PARAMETER :: P0 = 1000.0         ! Standard pressure at surface (hPa)
     REAL, PARAMETER :: lv = 2.54*(10**6)   ! Latent heat of vaporization
                                            ! (J kg^-1)
     REAL, PARAMETER :: cp = 1004.67        ! Specific heat of dry air constant
                                            ! at pressure (J/deg kg)
     REAL :: tlc,TLCL                            ! LCL temperature
  
     !~ Calculate the temperature of the LCL
     !  ------------------------------------
     tlc = TLCL ( tK, rh )
  
     !~ Calculate theta-e
     !  -----------------
     thetae = (tK * (p0/p)**( (R/Cp)*(1.- ( (.28E-3)*mixr*1000.) ) ) )* &
                 exp( (((3.376/tlc)-.00254))*&
                    (mixr*1000.*(1.+(.81E-3)*mixr*1000.)) )
  
  END FUNCTION Thetae
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  !~                                                                     
  !~ Name:                                                                
  !~    tlcl                                                               
  !~                                                                        
  !~ Description:                                                            
  !~    This function calculates the temperature of a parcel of air would have
  !~    if lifed dry adiabatically to it's lifting condensation level (lcl).  
  !~                                                                          
  !~ References:                                                              
  !~    Bolton (1980), Monthly Weather Review, pg. 1048, Eq. 22
  !~                                                                          
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  FUNCTION TLCL ( tk, rh )
    
    IMPLICIT NONE
 
    REAL, INTENT ( IN ) :: tK   !~ Temperature ( K )
    REAL, INTENT ( IN ) :: rh   !~ Relative Humidity ( % )
    REAL                :: tlcl
    
    REAL :: denom, term1, term2

    term1 = 1.0 / ( tK - 55.0 )
    IF ( rh > REAL (0) ) THEN
      term2 = ( LOG (rh/100.0)  / 2840.0 )
    ELSE
      term2 = ( LOG (0.001/1.0) / 2840.0 )
    END IF
    denom = term1 - term2
    tlcl = ( 1.0 / denom ) + REAL ( 55 ) 

  END FUNCTION TLCL
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  !~
  !~ Name:
  !~    The2T.f90
  !~
  !~ Description:
  !~    This function returns the temperature at any pressure level along a
  !~    saturation adiabat by iteratively solving for it from the parcel
  !~    thetae.
  !~
  !~ Dependencies:
  !~    function thetae.f90
  !~
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  FUNCTION The2T ( thetaeK, pres, flag ) result ( tparcel )
  IMPLICIT NONE
  
     !~ Variable Declaration
     !  --------------------
     REAL,    INTENT     ( IN ) :: thetaeK
     REAL,    INTENT     ( IN ) :: pres
     LOGICAL, INTENT ( INOUT )  :: flag
     REAL                       :: tparcel
  
     REAL :: thetaK
     REAL :: tovtheta
     REAL :: tcheck
     REAL :: svpr, svpr2
     REAL :: smixr, smixr2
     REAL :: thetae_check, thetae_check2
     REAL :: tguess_2, correction
  
     LOGICAL :: found
     INTEGER :: iter
  
     REAL :: R     ! Dry gas constant
     REAL :: Cp    ! Specific heat for dry air
     REAL :: kappa ! Rd / Cp
     REAL :: Lv    ! Latent heat of vaporization at 0 deg. C
     REAL :: Thetae
  
     R     = 287.04
     Cp    = 1004.67
     Kappa = R/Cp
     Lv    = 2.500E+6

     !~ Make initial guess for temperature of the parcel
     !  ------------------------------------------------
     tovtheta = (pres/100000.0)**(r/cp)
     tparcel  = thetaeK/exp(lv*.012/(cp*295.))*tovtheta

     iter = 1
     found = .false.
     flag = .false.

     DO
        IF ( iter > 105 ) EXIT

        tguess_2 = tparcel + REAL ( 1 )

        svpr   = 6.122 * exp ( (17.67*(tparcel-273.15)) / (tparcel-29.66) )
        smixr  = ( 0.622*svpr ) / ( (pres/100.0)-svpr )
        svpr2  = 6.122 * exp ( (17.67*(tguess_2-273.15)) / (tguess_2-29.66) )
        smixr2 = ( 0.622*svpr2 ) / ( (pres/100.0)-svpr2 )

        !  ------------------------------------------------------------------ ~!
        !~ When this function was orinially written, the final parcel         ~!
        !~ temperature check was based off of the parcel temperature and      ~!
        !~ not the theta-e it produced.  As there are multiple temperature-   ~!
        !~ mixing ratio combinations that can produce a single theta-e value, ~!
        !~ we change the check to be based off of the resultant theta-e       ~!
        !~ value.  This seems to be the most accurate way of backing out      ~!
        !~ temperature from theta-e.                                          ~!
        !~                                                                    ~!
        !~ Rentschler, April 2010                                             ~!
        !  ------------------------------------------------------------------  !

        !~ Old way...
        !thetaK = thetaeK / EXP (lv * smixr  /(cp*tparcel) )
        !tcheck = thetaK * tovtheta

        !~ New way
        thetae_check  = Thetae ( tparcel,  pres/100., 100., smixr  )
        thetae_check2 = Thetae ( tguess_2, pres/100., 100., smixr2 )

        !~ Whew doggies - that there is some accuracy...
        !IF ( ABS (tparcel-tcheck) < .05) THEN
        IF ( ABS (thetaeK-thetae_check) < .001) THEN
           found = .true.
           flag  = .true.
           EXIT
        END IF

        !~ Old
        !tparcel = tparcel + (tcheck - tparcel)*.3

        !~ New
        correction = ( thetaeK-thetae_check ) / ( thetae_check2-thetae_check )
        tparcel = tparcel + correction

        iter = iter + 1
     END DO

     !IF ( .not. found ) THEN
     !   print*, "Warning! Thetae to temperature calculation did not converge!"
     !   print*, "Thetae ", thetaeK, "Pressure ", pres
     !END IF

  END FUNCTION The2T
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  !~
  !~ Name:
  !~    VirtualTemperature
  !~
  !~ Description:
  !~    This function returns virtual temperature given temperature ( K )
  !~    and mixing ratio.
  !~
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  FUNCTION VirtualTemperature ( tK, w ) result ( Tv )
  IMPLICIT NONE

     !~ Variable declaration
     real, intent ( in ) :: tK !~ Temperature
     real, intent ( in ) :: w  !~ Mixing ratio ( kg kg^-1 )
     real                :: Tv !~ Virtual temperature

     Tv = tK * ( 1.0 + (w/0.622) ) / ( 1.0 + w )

  END FUNCTION VirtualTemperature
