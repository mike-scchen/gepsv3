 SUBROUTINE GRIBIT(F,LBM,IDRT,IGRID,IM,JM,MXBIT,COLAT1,             &
&                  ILPDS,IPTV,ICEN,IGEN,IBMS,                       &
&                  IPU,ITL,IL1,IL2,                                 &
&                  IYR,IMO,IDY,IHR,IFTU,IP1,IP2,ITR,                &
&                  INA,INM,ICEN2,IDS,IENS,                          &
&                  XLAT1,XLON1,XLAT2,XLON2,DELX,DELY,ORI,TRU,PROJ,  &
&                  GRIB,LGRIB,IERR)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:    GRIBIT      CREATE GRIB MESSAGE
!   PRGMMR: IREDELL          ORG: W/NMC23    DATE: 92-10-31
!
! ABSTRACT: CREATE A GRIB MESSAGE FROM A FULL FIELD.
!   AT PRESENT, ONLY GLOBAL LATLON GRIDS AND GAUSSIAN GRIDS
!   AND REGIONAL POLAR PROJECTIONS ARE ALLOWED.
!
! PROGRAM HISTORY LOG:
!   92-10-31  IREDELL
!   94-05-04  JUANG (FOR GSM AND RSM USE)
!
! USAGE:    CALL GRIBIT(F,LBM,IDRT,IGRID,IM,JM,MXBIT,COLAT1,
!    &                  ILPDS,IPTV,ICEN,IGEN,IBMS,IPU,ITL,IL1,IL2,
!    &                  IYR,IMO,IDY,IHR,IFTU,IP1,IP2,ITR,
!    &                  INA,INM,ICEN2,IDS,IENS,
!    &                  XLAT1,XLON1,DELX,DELY,ORI,TRU,PROJ,
!    &                  GRIB,LGRIB,IERR)
!   INPUT ARGUMENT LIST:
!     F        - REAL (IM*JM) FIELD DATA TO PACK INTO GRIB MESSAGE
!     LBM      - LOGICAL (IM*JM) BITMAP TO USE IF IBMS=1
!     IDRT     - INTEGER DATA REPRESENTATION TYPE
!                (0 FOR LATLON OR 4 FOR GAUSSIAN OR 5 FOR POLAR)
!     IGRID    - grid number
!     IM       - INTEGER LONGITUDINAL DIMENSION
!     JM       - INTEGER LATITUDINAL DIMENSION
!     MXBIT    - INTEGER MAXIMUM NUMBER OF BITS TO USE (0 FOR NO LIMIT)
!     COLAT1   - REAL FIRST COLATITUDE OF GRID IF IDRT=4 (RADIANS)
!     ILPDS    - INTEGER LENGTH OF THE PDS (USUALLY 28)
!     IPTV     - INTEGER PARAMETER TABLE VERSION (USUALLY 1)
!     ICEN     - INTEGER FORECAST CENTER (USUALLY 7)
!     IGEN     - INTEGER MODEL GENERATING CODE
!     IBMS     - INTEGER BITMAP FLAG (0 FOR NO BITMAP)
!     IPU      - INTEGER PARAMETER AND UNIT INDICATOR
!     ITL      - INTEGER TYPE OF LEVEL INDICATOR
!     IL1      - INTEGER FIRST LEVEL VALUE (0 FOR SINGLE LEVEL)
!     IL2      - INTEGER SECOND LEVEL VALUE
!     IYR      - INTEGER YEAR
!     IMO      - INTEGER MONTH
!     IDY      - INTEGER DAY
!     IHR      - INTEGER HOUR
!     IFTU     - INTEGER FORECAST TIME UNIT (1 FOR HOUR)
!     IP1      - INTEGER FIRST TIME PERIOD
!     IP2      - INTEGER SECOND TIME PERIOD (0 FOR SINGLE PERIOD)
!     ITR      - INTEGER TIME RANGE INDICATOR (10 FOR SINGLE PERIOD)
!     INA      - INTEGER NUMBER INCLUDED IN AVERAGE
!     INM      - INTEGER NUMBER MISSING FROM AVERAGE
!     ICEN2    - INTEGER FORECAST SUBCENTER
!                (USUALLY 0 BUT 1 FOR REANAL OR 2 FOR ENSEMBLE)
!     IDS      - INTEGER DECIMAL SCALING
!     IENS     - INTEGER (5) ENSEMBLE EXTENDED PDS VALUES
!                (APPLICATION,TYPE,IDENTIFICATION,PRODUCT,SMOOTHING)
!                (USED ONLY IF ICEN2=2 AND ILPDS>=45)
!     XLAT1    - REAL FIRST POINT OF REGIONAL LATITUDE (RADIANS)
!     XLON1    - REAL FIRST POINT OF REGIONAL LONGITUDE (RADIANS)
!     XLAT2    - REAL LAST  POINT OF REGIONAL LATITUDE (RADIANS)
!     XLON2    - REAL LAST  POINT OF REGIONAL LONGITUDE (RADIANS)
!     DELX     - REAL DX ON 60N FOR REGIONAL (M)
!     DELY     - REAL DY ON 60N FOR REGIONAL (M)
!     PROJ     - REAL POLAR PROJECTION FLAG 1 FOR NORTH -1 FOR SOUTH
!                     MERCATER PROJECTION 0
!     ORI      - REAL ORIENTATION OF REGIONAL POLAR PROJECTION OR
!     TRU      - TRUTH FOR REGIONAL MERCATER AND LAMBERT PROJECTION
!
!   OUTPUT ARGUMENT LIST:
!     GRIB     - CHARACTER (LGRIB) GRIB MESSAGE
!     LGRIB    - INTEGER LENGTH OF GRIB MESSAGE
!                (NO MORE THAN 100+ILPDS+IM*JM*(MXBIT+1)/8)
!     IERR     - INTEGER ERROR CODE (0 FOR SUCCESS)
!
! SUBPROGRAMS CALLED:
!   GTBITS     - COMPUTE NUMBER OF BITS AND ROUND DATA APPROPRIATELY
!   W3FI72     - ENGRIB DATA INTO A GRIB1 MESSAGE
!
! ATTRIBUTES:
!   LANGUAGE: CRAY FORTRAN
!
!$$$
      use grib_mod
      INTEGER IENS(5)
      REAL F(IM*JM)
      LOGICAL LBM(IM*JM)
      CHARACTER GRIB(*)
      INTEGER IBM(IM*JM*IBMS+1-IBMS),IPDS(100),IBDS(100)
      integer::IGDS(5)
      REAL FR(IM*JM)
      CHARACTER PDS(ILPDS)
!grib2
!setction 0 & 1
integer::listsec0(2)=(/0,2/)
integer::listsec1(13)
!

integer*4,parameter :: igdstmplen=22
integer*4 :: igdstmpl(igdstmplen)
integer*4 :: ideflist=0,idefnum=1,ipdsnum

integer*4,parameter :: ipdstmplen=29,ipdstmplenacc=32
integer*4::ipdstmpl(ipdstmplen),ipdstmplacc(ipdstmplenacc)
integer*4,parameter :: numcoord=1
real*4 :: coordlist(numcoord)
integer*4::idrsnum3=3
integer*4,parameter :: idrstmplen3=17
integer*4::idrstmpl3(idrstmplen3)
!GRIB2 SECTION 6
character,allocatable,save :: cgrib(:)*1
integer*4:: lcgrib=1*1e8 ,lengrib
integer*4 :: ibmap=255
logical*1,allocatable,save :: bmap(:)
integer::ngrdpts,NF
integer*8::idtg,idtg2
character::cdtg*12,cdtg2*12
integer::p0,p1,p2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  DETERMINE GRID PARAMETERS
 PI=ACOS(-1.)
 NF=IM*JM
 if( .not. allocated( bmap ) )allocate(bmap(NF))
 IGDS( 1)=0
 IGDS( 2)=NF   !number of grid point
 IGDS( 3)=0
 IGDS( 4)=0
!------ section 3
 igdstmpl(01)=6          !Shpape of the Earth( See Code talbe 3.2) 
 igdstmpl(02)=0          !
 igdstmpl(03)=6370000    !
 igdstmpl(04)=0          !
 igdstmpl(05)=0          !
 igdstmpl(06)=0          !
 igdstmpl(07)=0          !
 igdstmpl(08)=IM         ! EAST-WEST POINTS
 igdstmpl(09)=JM         ! NORTH-SOUTH POINTS
 IF(IDRT.EQ.0) THEN    !lat-lon mode
   IGDS( 5)=0
!   IF(IM.EQ.144.AND.JM.EQ.73) THEN
!     IGRID=2
!   ELSEIF(IM.EQ.360.AND.JM.EQ.181) THEN
!     IGRID=3
!!  ELSE
!!    IGRID=255
!   ENDIF
   IRESFL=128
   ISCAN=0
   IGDS09=-NINT(90.E3)
   IGDS10=0
   IGDS11=NINT(180.E3/(JM-1))
   IGDS12=NINT(360.E3/IM)
   !change to regional latlon
   LAT1=NINT(180.E3/ACOS(-1.)*XLAT1)
   LON1=NINT(180.E3/ACOS(-1.)*XLON1)
   LAT2=NINT(180.E3/ACOS(-1.)*XLAT2)
   LON2=NINT(180.E3/ACOS(-1.)*XLON2)
   IGDS11=(IGDS09-LAT1)/(JM-1)
   IGDS12=(IGDS10-LON1)/(IM-1)
   IGDS13=64
! 
   igdstmpl(10) = 0   ! Basic Angle of init projection (not important to us)
   igdstmpl(11) = 0   ! Subdivision of basic angle
   igdstmpl(12) = LAT1*1e+3   ! LATITUDE OF first grid point
   igdstmpl(13) = LON1*1e+3   ! LONGITUDE OF first grid point
   igdstmpl(14) = IRESFL  ! Resolution and component flags
   igdstmpl(15) = LAT2*1e+3 !latitude of last grid point
   igdstmpl(16) = LON2*1e+3 !longitude of last grid point 
   igdstmpl(17) = IGDS12*1e3 ! i-direction increment in micro degs
   igdstmpl(18) = IGDS11*1e3 ! j-direction increment in micro degs
   igdstmpl(19) = 64   ! Scanning mode
 ELSEIF(IDRT.EQ.4) THEN
   IF(IM.EQ.192.AND.JM.EQ.94) THEN
     IGRID=98
   ELSEIF(IM.EQ.384.AND.JM.EQ.190) THEN
     IGRID=126
!  ELSE
!    IGRID=255
   ENDIF
   IRESFL=128
   ISCAN=0
   LAT1=NINT(90.E3-180.E3/PI*COLAT1)
   LON1=0
   LATI=JM/2
   LONI=NINT(360.E3/IM)
   IGDS09=-LAT1
   IGDS10=-LONI
   IGDS11=LATI
   IGDS12=LONI
   IGDS13=ISCAN
   IGDS14=0
   IGDS15=0
   IGDS16=0
   IGDS17=0
   IGDS18=0
 ELSEIF(IDRT.EQ.5) THEN    ! POLAR PROJECTION
!  print*,'POLAR PROJECTION'
!  print*,'PROJ=',PROJ
!  IGRID=255
   IGDS( 5)=20
   LAT1=NINT(180.E3/ACOS(-1.) * XLAT1)
   LON1=NINT(180.E3/ACOS(-1.) * XLON1)
   IRESFL=8
   IGDS09=NINT(ORI*1.E3)
   IGDS10=DELX
   IGDS11=DELY
   IGDS12=0
   IF( NINT(PROJ).EQ.1  ) IGDS12=0         ! NORTH POLAR PROJ
   IF( NINT(PROJ).EQ.-1 ) IGDS12=128       ! SOUTH POLAT PROJ
   ISCAN=64
   IGDS13=ISCAN
   IGDS14=0
   IGDS15=0
   IGDS16=0
   IGDS17=0
   IGDS18=0
!??????????
   center_lat=23.5
   central_lon=120.5
!
   igdstmpl(10) = LAT1*1e+3
   igdstmpl(11) = LON1*1e+3
   igdstmpl(12) = 128 ! Resolution and component flag
   igdstmpl(13) = LAT1*1e+3 ! Latitude where Dx and Dy are specified
   igdstmpl(14) = central_lon
   igdstmpl(15) = IGDS10*1e+03  ! x-dimension grid-spacing  in units of m^-3
   igdstmpl(16) = IGDS11*1e+03
   if (center_lat .lt. 0) then
      igdstmpl(17) = 1
   else
      igdstmpl(17) = 0
   endif
   igdstmpl(18) = 64   ! Scanning mode
 ELSEIF(IDRT.EQ.1) THEN    ! MERCATER PROJECTION
   IGDS( 5)=10
   LAT1=NINT(180.E3/ACOS(-1.) * XLAT1)
   LON1=NINT(180.E3/ACOS(-1.) * XLON1)
   IRESFL=128
   LAT2=NINT(180.E3/ACOS(-1.) * XLAT2)
   LON2=NINT(180.E3/ACOS(-1.) * XLON2)
   IGDS11=DELX
   IGDS12=DELY
   IGDS13=NINT(TRU*1.E3)
   ISCAN=64
   igdstmpl(10)=LAT1*1e+3   ! LATITUDE OF first grid point
   igdstmpl(11)=LON1*1e+3   ! LONGITUDE OF first grid point
   igdstmpl(12)=IRESFL     ! RESOLUTION FLAG
   igdstmpl(13)=LAT1*1e3   !Mercator projection intersects the Earth
   igdstmpl(14)=LAT2*1e3 ! LATITUDE OF END OR ORIENTATION
   igdstmpl(15)=LON2*1e3 ! LONGITUDE OF END OR DX IN METER ON 60N
   igdstmpl(16)=ISCAN      !Scanning mode
   igdstmpl(17)=0          !Orientation of the grid between i-direction and equator
   igdstmpl(18)=IGDS11*1e3 !longitudinal direction grid length
   igdstmpl(19)=IGDS12*1e3 !latitudinal direction grid length
 ELSEIF(IDRT.EQ.3) THEN    ! LAMBERT PROJECTION
   IGDS( 5)=30
   LAT1=NINT(180.E3/ACOS(-1.)*XLAT1)
   LON1=NINT(180.E3/ACOS(-1.)*XLON1)
   LAT2=NINT(180.E3/ACOS(-1.)*XLAT2)
   LON2=NINT(180.E3/ACOS(-1.)*XLON2)
!  print*,'LAT1,LON1=',LAT1,LON1
   IRESFL=8
   IGDS09=NINT(ORI*1.E3)
!  print*,'IGDS09=',IGDS09
   IGDS10=DELX
   IGDS11=DELY
!  print*,'IGDS10,IGDS11=',IGDS10,IGDS11
   IGDS12=0
   IF( NINT(PROJ).EQ.2  ) IGDS12=0         ! NORTH LAMBERT PROJ
   IF( NINT(PROJ).EQ.-2 ) IGDS12=128       ! SOUTH LAMBERT PROJ
!  print*,'IGDS12=',IGDS12
   ISCAN=64
   IGDS13=ISCAN
   IGDS14=0
   IGDS15=NINT(TRU*1.E3)
   IGDS16=IGDS15
!  print*,'IGDS15,IGDS16=',IGDS15,IGDS16
!  IGDS17=-90000    !lat of South Pole (MILLIDEGREES)
   IGDS17=0         !???????
   IGDS18=0         !lon of South Pole (MILLIDEGREES)
!  ???????????????????
   central_lon=120.5 !?????
   center_lat=23.5 !?????
!
   igdstmpl(10) = LAT1*1e+3   ! LATITUDE OF first grid point
   igdstmpl(11) = LON1*1e+3   ! LONGITUDE OF first grid point
   igdstmpl(12) = 128 ! Resolution and component flag
   igdstmpl(13) = latin1 ! Latitude where Dx and Dy are specified
   igdstmpl(14) = central_lon 
   igdstmpl(15) = IGDS10*1e+3  ! x-dimension grid-spacing  in units of m^-3
   igdstmpl(16) = IGDS11*1e+3
   if (center_lat .lt. 0) then
      igdstmpl(17) = 1
   else
      igdstmpl(17) = 0
   endif
   igdstmpl(18) = 64   ! Scanning mode
   igdstmpl(19) = lat1*1e+03
   igdstmpl(20) = lat2*1e+03
   igdstmpl(21) = -90*1e+06
   igdstmpl(22) = central_lon
 ELSE
   IERR=40
   RETURN
 ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  RESET TIME RANGE PARAMETER IN CASE OF OVERFLOW
      IF(ITR.GE.2.AND.ITR.LE.5.AND.IP2.GE.256) THEN
        JP1=IP2
        JP2=0
        JTR=10
      ELSE
        JP1=IP1
        JP2=IP2
        JTR=ITR
      ENDIF
!  FIX YEAR AND CENTURY
      IYR4=IYR
      IF(IYR.LE.100) IYR4=2050-MOD(2050-IYR,100)
      IYC=MOD(IYR4-1,100)+1
      ICY=(IYR4-1)/100+1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!grib2 section 1
 listsec1( 1)=ICEN   !center ID
 listsec1( 2)=ICEN2  !sub center ID
 listsec1( 3)=2
 listsec1( 4)=1
 listsec1( 5)=1 
 listsec1( 6)=IYR
 listsec1( 7)=IMO
 listsec1( 8)=IDY
 listsec1( 9)=IHR
 listsec1(10)=0
 listsec1(11)=0
 listsec1(12)=0
 if(IP1==0)then
   listsec1(13)=0  !:anl:
 else
   listsec1(13)=1  !:HH hours fcst:
 endif
!
!grib2 section 4
!print*,'IP1=',ip1,'IP2=',ip2
!0:normal 8:acc_time
coordlist=1
call grb1to2table(IPU,p0,p1,p2)
listsec0(1)=p0

if(IPTV==133)then
  !print*,'IPTV==133'
  select case(IPU)
    case(204) ;p1=1  ;  p2=219
    case(205) ;p1=1  ;  p2=220
  end select 
endif
ipdstmpl( 1)= p1
ipdstmpl( 2)= p2
if(IP1==0)then
  ipdstmpl( 3)=0
else
  ipdstmpl( 3)=2
endif
ipdstmpl( 4)=0
ipdstmpl( 5)=ICEN2
ipdstmpl( 6)=0
ipdstmpl( 7)=0
ipdstmpl( 8)=1
if(ip2==0)then
  ipdsnum=0 
  ipdstmpl( 9)=IP1
else
  ipdsnum=8
  !ipdstmpl( 9)=IP2
  ipdstmpl( 9)=IP1
endif
ipdstmpl(10)=ITL   !layer type
!some fixed for diffrent between grib1 and grib2
if(ITL==103)ipdstmpl(10)=102 !above mean sea level
if(ITL==105)ipdstmpl(10)=103 !above ground 
if(ITL==108)ipdstmpl(10)=104 !sigma layer
if(ITL==112)ipdstmpl(10)=106 !Depth Below Land Surface (m)
if(ITL==116)ipdstmpl(10)=108 !Level at Specified Pressure Difference from Ground to Level
if(ITL==200)ipdstmpl(10)=10  !Entire Atmosphere
!print*,'ITL=',ITL
!layer set
!print*,'IL1=',IL1,'IL2=',IL2
if(ipdstmpl(10)==100.or.ipdstmpl(10)==108)then
  ipdstmpl(11)=-2 !change hPa to Pa for pressure level
elseif(ipdstmpl(10)==106)then
  ipdstmpl(11)=2  !cm to m
elseif(ipdstmpl(10)==104)then
  ipdstmpl(11)=2  !sigma layer value 72 to 0.72
else
  ipdstmpl(11)=0     !scale factor
endif
if(IL1==0)then
 ipdstmpl(12)=IL2   !value
 ipdstmpl(13)=255
 ipdstmpl(14)=0
 ipdstmpl(15)=0
else
  !print *,'IL1,IL2',ITL,IL1,IL2
  ipdstmpl(12)=IL1   !value
  ipdstmpl(13)=ipdstmpl(10)
  ipdstmpl(14)=ipdstmpl(11)
  ipdstmpl(15)=IL2
endif
if(ipdsnum==8)then  !accum precipitation
  write(cdtg,'(I4.4,I2.2,I2.2,I2.2,A2)')IYR,IMO,IDY,IHR,'00'
  call dtgfix12(cdtg,cdtg2 ,IP2 )
  read(cdtg2(1:10),'(i4,i2,i2,i2)')ipdstmpl(16:19)!yyyy mm dd hh
  ipdstmpl(20)=0 !t20 ! minute
  ipdstmpl(21)=0 !t21 ! Second | Time of end of overall time interval
  ipdstmpl(22)=1 !t22
  ipdstmpl(23)=0 !t23
  ipdstmpl(24)=0 !t24  0 ave 1 accu
  ipdstmpl(25)=2 !t25 code 4.11
  ipdstmpl(26)=1 !t26 ! unit of time for time range. 1:hour, 2:day
  ipdstmpl(27)=IP2-IP1  ! time interval value
  ipdstmpl(28)=0  !t28 fcst,dt unit 0:minut 1:hour 2:day
  ipdstmpl(29)=0  !t29
endif
!grib2 section 5
!print*,'IDS=',IDS
 idrstmpl3=(/0,0,IDS,0,0,0,0,0,0,0,0,0,0,0,0,0,2/)


!!not use
!!  FILL PDS PARAMETERS
!      IPDS(01)=ILPDS    ! LENGTH OF PDS
!      IPDS(02)=IPTV     ! PARAMETER TABLE VERSION ID
!      IPDS(03)=ICEN     ! CENTER ID
!      IPDS(04)=IGEN     ! GENERATING MODEL ID
!      IPDS(05)=IGRID    ! GRID ID
!      IPDS(06)=1        ! GDS FLAG, USUALLY 1
!      IPDS(07)=IBMS     ! BMS FLAG
!      IPDS(08)=IPU      ! PARAMETER UNIT ID
!      IPDS(09)=ITL      ! TYPE OF LEVEL ID
!      IPDS(10)=IL1      ! LEVEL 1 OR 0
!      IPDS(11)=IL2      ! LEVEL 2
!      IPDS(12)=IYC      ! YEAR
!      IPDS(13)=IMO      ! MONTH
!      IPDS(14)=IDY      ! DAY
!      IPDS(15)=IHR      ! HOUR
!      IPDS(16)=0        ! MINUTE
!      IPDS(17)=IFTU     ! FORECAST TIME UNIT ID
!      IPDS(18)=JP1      ! TIME PERIOD 1
!      IPDS(19)=JP2      ! TIME PERIOD 2 OR 0
!      IPDS(20)=JTR      ! TIME RANGE INDICATOR
!      IPDS(21)=INA      ! NUMBER IN AVERAGE
!      IPDS(22)=INM      ! NUMBER MISSING
!      IPDS(23)=ICY      ! CENTURY
!      IPDS(24)=ICEN2    ! FORECAST SUBCENTER
!      IPDS(25)=IDS      ! DECIMAL SCALING
!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! IBDS(1:9)=0       ! BDS FLAGS
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  FILL BITMAP AND COUNT VALID DATA.  RESET BITMAP FLAG IF ALL VALID.
!      NBM=NF
!      IF(IBMS.NE.0) THEN
!        NBM=0
!        DO I=1,NF
!          IF(LBM(I)) THEN
!            IBM(I)=1
!            NBM=NBM+1
!          ELSE
!            IBM(I)=0
!          ENDIF
!        ENDDO
!        IF(NBM.EQ.NF) IPDS(7)=0
!      ENDIF
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  ROUND DATA AND DETERMINE NUMBER OF BITS
!!      IF(NBM.EQ.0) THEN
!!        DO I=1,NF
!!          FR(I)=0.
!!        ENDDO
!!        NBIT=0
!!      ELSE
!!        CALL GTBITS(IPDS(7),IDS,NF,IBM,F,FR,FMIN,FMAX,NBIT)
!!!       WRITE(0,'("GTBITS:",4I4,4X,2I4,4X,2G16.6)')
!!!    &   IPU,ITL,IL1,IL2,IDS,NBIT,FMIN,FMAX
!!        IF(MXBIT.GT.0) NBIT=MIN(NBIT,MXBIT)
!!      ENDIF
!
!! FR(i) is packed.
!!     write(6,*) (FR(i),i=1,10)
!!     write(6,*) 'NBIT= ',NBIT
!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  CREATE PRODUCT DEFINITION SECTION
!      CALL W3FI68(IPDS,PDS)
!      IF(ICEN2.EQ.2.AND.ILPDS.GE.45) THEN
!        ILAST=45
!        CALL PDSENS(IENS,KPROB,XPROB,KCLUST,KMEMBR,ILAST,PDS)
!      ENDIF
!      CALL W3FI72(0,FR,0,NBIT,1,IPDS,PDS,                                &
!     &            1,255,IGDS,0,0,IBM,NF,IBDS,                            &
!     &            NFO,GRIB,LGRIB,IERR)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call gribcreate(GRIB,lcgrib,listsec0,listsec1,ierr)
call addgrid(GRIB,lcgrib,igds,igdstmpl,igdstmplen,ideflist,idefnum,ierr)
call addfield(GRIB,lcgrib,ipdsnum,ipdstmpl,ipdstmplen, &
     coordlist,numcoord,idrsnum3,idrstmpl3,idrstmplen3, &
     F,NF,ibmap,bmap,ierr)
!call gribend(cgrib,lcgrib,lengrib,ierr)
call gribend(GRIB,lcgrib,LGRIB,ierr)

!pass to output
!turn off output layer type 116
!if(ITL==116)LGRIB=0  
!turn off output variable
!if(IPU==140)LGRIB=0 !CRAIN
!if(IPU==141)LGRIB=0 !CFRZR
!if(IPU==142)LGRIB=0 !CICEP
!if(IPU==143)LGRIB=0 !CSNOW
!if(IPU==222)LGRIB=0 !5WAVH
!
!if(IPTV==131)LGRIB=0 !WVUFLX WVVFLX
RETURN
END


