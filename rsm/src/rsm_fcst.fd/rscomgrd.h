#ifdef MP
#define IGRD1S igrd1p
#define JGRD1S jgrd1p
#define LNGRDS lngrdp
#else
#define IGRD1S igrd1
#define JGRD1S jgrd1
#define LNGRDS lngrd
#endif
!.................................................................
!....begin regional grid common....................
!
      common /rcomgd1/                                           &   
     & flat(LNGRDS), flon(LNGRDS),                               &
     &  fm2(LNGRDS), fm2x(LNGRDS), fm2y(LNGRDS)
!....
      common /rcomgd2/                                           &
     &  slmsk(IGRD1S,JGRD1S),                                    &
     & hprime(IGRD1S,JGRD1S,14),                                 &
     &    swh(IGRD1S,levr,JGRD1S),hlw(IGRD1S,levr,JGRD1S),       &
     & sfcnsw(IGRD1S,JGRD1S),   sfcdlw(IGRD1S,JGRD1S),           &
     & sinlar(IGRD1S,JGRD1S),   coslar(IGRD1S,JGRD1S),           &
     &                            coszer(IGRD1S,JGRD1S),         &
     &    acv(IGRD1S,JGRD1S),     acvt(IGRD1S,JGRD1S),           &
     &     cv(IGRD1S,JGRD1S),      cvt(IGRD1S,JGRD1S),           &
     &    cvb(IGRD1S,JGRD1S),     acvb(IGRD1S,JGRD1S),           &
     &  tsflw(IGRD1S,JGRD1S),     f10m(IGRD1S,JGRD1S),           &
     &  sdec,cdec,slag,solhr,clstp,                              &
!yj---get sfc downward short wave radiation and snowfree albedo--->>>     
!radiation_run
     & sswdnf(IGRD1S,JGRD1S),    sfalb(IGRD1S,JGRD1S)
!yj---get sfc downward short wave radiation and snowfree albedo---<<<     
!....
      common /rcomgd3/                                           &
     &  dusfc(IGRD1S,JGRD1S),    dvsfc(IGRD1S,JGRD1S),           &
     &  dtsfc(IGRD1S,JGRD1S),    dqsfc(IGRD1S,JGRD1S),           &
     & dlwsfc(IGRD1S,JGRD1S),   ulwsfc(IGRD1S,JGRD1S),           &
     & geshem(IGRD1S,JGRD1S),     tsea(IGRD1S,JGRD1S),           &
     &    ssu(IGRD1S,JGRD1S),      ssv(IGRD1S,JGRD1S),           &
     &  dugwd(IGRD1S,JGRD1S),    dvgwd(IGRD1S,JGRD1S),           &
     &   u10m(IGRD1S,JGRD1S),     v10m(IGRD1S,JGRD1S),           &
     &    t2m(IGRD1S,JGRD1S),      q2m(IGRD1S,JGRD1S),           &
     &  psurf(IGRD1S,JGRD1S),   psmean(IGRD1S,JGRD1S),           &
     & cflash(IGRD1S,JGRD1S)
!....
      common /rcomgd4/                                           &
     &    tg3(IGRD1S,JGRD1S),     zorl(IGRD1S,JGRD1S),           &
     & sheleg(IGRD1S,JGRD1S),   bengsh(IGRD1S,JGRD1S),           &
     &  gflux(IGRD1S,JGRD1S),    slrad(IGRD1S),                  &
     &    smc(IGRD1S,JGRD1S,lsoil),                              &
     &    stc(IGRD1S,JGRD1S,lsoil),                              &
     & canopy(IGRD1S,JGRD1S),   runoff(IGRD1S,JGRD1S),           &
     &srunoff(IGRD1S,JGRD1S),    soilm(IGRD1S,JGRD1S),           &
     & snwdph(IGRD1S,JGRD1S),                                    &
     & tmpmax(IGRD1S,JGRD1S),   tmpmin(IGRD1S,JGRD1S),           &
     & spfhmax(IGRD1S,JGRD1S), spfhmin(IGRD1S,JGRD1S),           &
     &  wvuflx(IGRD1S,JGRD1S),  wvvflx(IGRD1S,JGRD1S),           &
     &     ep(IGRD1S,JGRD1S),   cldwrk(IGRD1S,JGRD1S),           &
     &   hpbl(IGRD1S,JGRD1S),     pwat(IGRD1S,JGRD1S),           &
     &   gust(IGRD1S,JGRD1S),                                    &
     &  alvsf(IGRD1S,JGRD1S),    alvwf(IGRD1S,JGRD1S),           &
     &  alnsf(IGRD1S,JGRD1S),    alnwf(IGRD1S,JGRD1S),           &
     &  facsf(IGRD1S,JGRD1S),    facwf(IGRD1S,JGRD1S),           &    
     &  vfrac(IGRD1S,JGRD1S),    vtype(IGRD1S,JGRD1S),           &  
     &  stype(IGRD1S,JGRD1S),   uustar(IGRD1S,JGRD1S),           &    
     &   ffmm(IGRD1S,JGRD1S),     ffhh(IGRD1S,JGRD1S),           &
!yj-implement noah land model---------------------------------------->>>
     &  slope(IGRD1S,JGRD1S),   shdmin(IGRD1S,JGRD1S),           &
     &  shdmax(IGRD1S,JGRD1S),  snoalb(IGRD1S,JGRD1S),           &
!skip nst:     &     oro(IGRD1S,JGRD1S),  oro_uf(IGRD1S,JGRD1S)
!yj-implement noah land model----------------------------------------<<<
!
!yj-implement 3DRT package   ---------------------------------------->>>
     &  sky_view(IGRD1S,JGRD1S),terrain_config(IGRD1S,JGRD1S),   &
     &  sinsl_cosas(IGRD1S,JGRD1S),sinsl_sinas(IGRD1S,JGRD1S)
!yj-implement 3DRT package   ----------------------------------------<<<
!.... end of rcomgdn
      common /rcomgd5/                                           &   
     & rtgrd(LNGRDS,levr,ntotal), rqgrd(LNGRDS,levr,ntotal),     &
     & rmgrd(LNGRDS,levr,ntotal) 
#ifdef CPL_TIMCOM
      common /rcomgd6/                                           &
     & cpl_field(IGRD1S,JGRD1S,9)
#endif
