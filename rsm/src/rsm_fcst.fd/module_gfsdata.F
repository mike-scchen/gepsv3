#include <defopt.h>
      module module_gfsdata
      use physcons, pi => con_pi
      
      implicit none

!#include <paramodel.h>

      public

! GFS data send/write tag/record index
      integer   kfh,knt,ktt,ksgeo,ksfcp,kuu,kvv
      integer   kqt,kcw,krw,kiw,ksw,kgl,koz
      integer   ktg,ksmc,ksnr,kstc,ksimk,kslmk
      integer   ktrace,kend
      integer   klon,klat

! GFS data vertical coordinate and horizontal lat/lon info
      integer   glev
      parameter(glev=72)
! for new NoahLSM (NOAH MP), 4 soil layers are needed, but for OSU LSM only 2
! soil layers are need. For convenient, GFS always gives RSM 4 soil
! layer(kmsoil), but lsoil will be 4 for NOAHLSM and will be 2 for OSU.
      integer kmsoil
      parameter(kmsoil=4)
      real(kind=8), dimension(:)  , allocatable, save  :: gaki,gbki
      real, allocatable:: clat(:,:),clon(:,:)
      
      integer:: idate
      real:: fhour

!yj2019: add terr and slmsk from gfs
      real, allocatable::  &
     &   sgeo(:,:), psfc(:,:), temp(:,:,:), u(:,:), v(:,:), &
     &   rqgb(:,:,:,:)
!    &   spfh(:,:), clwr(:,:),                        &
!    &   rain(:,:), qice(:,:), snow(:,:), grpl(:,:),  &
!    &   ozon(:,:)
      real, allocatable:: simk(:,:)
!    &   slmsk(:), simk(:), tg(:),smc(:,:),snr(:),stc(:,:)

      real, allocatable::                                               &
     &   tg3clm(:),zorclm(:),albclm(:,:),vegclm(:),vetclm(:),           &
     &   sotclm(:),alfclm(:,:),vmnclm(:),vmxclm(:),slpclm(:),absclm(:)

!ty2022 : add lat lon and slmsk from gfs
      real, allocatable::                                               &
     &   fflat(:),fflon(:),rslm(:)

      save

      contains

        subroutine config_gfsdata(mype,master,root_gfs,cigrd1,cjgrd1,itag)
        integer mype,master,cigrd1,cjgrd1
        integer ierr,i,j
#ifdef CWB_CFS
        real, allocatable:: xlonin(:,:),xlatin(:,:)
#else
        real, allocatable:: xlonin(:),xlatin(:)
#endif
!for defined CWB_MPMD or CWB_CFS
        integer itag,root_gfs
        !
        ! GFS data configuration
        ! allocate GFS horizontal coordinate setting
        allocate (clat(cigrd1,cjgrd1),clon(cigrd1,cjgrd1),stat=ierr)
      if (ierr/=0) stop "module_gfsdata: allocate horizontal vars error"
        !
        ! GFS data send/write tag/record order
          ktrace=-999
        ! |___lat/lon
          klon=2
          klat=3
        ! |___sig
          kfh=1
          knt=kfh+1
          ksgeo=knt+1
          ksfcp=ksgeo+1
          ktt=ksfcp+1
          kuu=ktt+glev
          kvv=kuu+glev
          kqt=kvv+glev
          kcw=kqt+glev
          koz=kcw+glev
        ! |___sfc
          kslmk=koz+glev
          ksnr=kslmk+1
          ksimk=ksnr+1
          ktg=ksimk+1
          ksmc=ktg+1
          kstc=ksmc+kmsoil
        ! |___sig-tracers
          krw=kstc+kmsoil
          kiw=krw+glev
          ksw=kiw+glev
          kgl=ksw+glev
          kend=kgl+glev

        if ( mype.eq.master ) then
          print*,'klon,klat='
          print*,klon,klat
          print*,'ktrace,kfh,knt,ktt,ksgeo,ksfcp,kuu,kvv,kqt,kcw,koz='
          print*,ktrace,kfh,knt,ktt,ksgeo,ksfcp,kuu,kvv,kqt,kcw,koz
          print*,'ktg,ksmc,ksnr,kstc,ksimk,kslmk,krw,kiw,ksw,kgl,kend='
          print*,ktg,ksmc,ksnr,kstc,ksimk,kslmk,krw,kiw,ksw,kgl,kend
        endif

        ! get info of GFS data lat/lon
        if ( mype.eq.master ) then

#if defined (CWB_MPMD) || defined (CWBSUM)
          allocate (xlonin(cigrd1),xlatin(cjgrd1),stat=ierr)
          if (ierr/=0) stop "module_gfsdata: allocate xlon/ylat error"
          itag=klon
          call mpmd_recv(xlonin,cigrd1,root_gfs,itag,'R')
          print*,'klon itag=',itag
          itag=klat
          call mpmd_recv(xlatin,cjgrd1,root_gfs,itag,'R')
          print*,'klat itag=',itag
#else
#ifdef CWB_CFS
          print*,'TYW config_gfsdata, before recv xlon,cigrd1,cjgrd1=' &
     &          ,cigrd1,cjgrd1
          allocate (xlonin(cigrd1,1),xlatin(1,cjgrd1),stat=ierr)
          if (ierr/=0) stop "module_gfsdata: allocate xlon/ylat error"
          call CMP_RECV(xlonin,cigrd1)
          print*,'recv klon, cigrd1=',cigrd1
          print*,'and klon=',xlonin
          call CMP_RECV(xlatin,cjgrd1)
          print*,'recv klat, cjgrd1=',cjgrd1
          print*,'and klat=',xlatin
#else
          allocate (xlonin(cigrd1),xlatin(cjgrd1),stat=ierr)
          if (ierr/=0) stop "module_gfsdata: allocate xlon/ylat error"
          open(99,file='g2r_xlon.txt',status='old', &
              form='formatted',iostat=ierr)
          if (ierr/=0) stop "module_gfsdata: open file xlon fail"
          read(99,*) xlonin
          close(99)
          open(98,file='g2r_xlat.txt',status='old', &
              form='formatted',iostat=ierr)
          if (ierr/=0) stop "module_gfsdata: open file xlat fail"
          read(98,*) xlatin
          close(98)
#endif
#endif
          do j=1,cjgrd1
          do i=1,cigrd1
#ifdef CWB_CFS
            clat(i,j)=xlatin(1,j)/180.*pi
            clon(i,j)=xlonin(i,1)/180.*pi
#else
            clat(i,j)=xlatin(j)/180.*pi
            clon(i,j)=xlonin(i)/180.*pi
#endif
          enddo
          enddo
          deallocate (xlonin,xlatin,stat=ierr)
          if (ierr/=0) stop "module_gfsdata: deallocate xlon/ylat error"

        endif

        call rmpbcastr(clat,cigrd1*cjgrd1)
        call rmpbcastr(clon,cigrd1*cjgrd1)
        !
        ! allocate GFS vertical coordinate setting
        allocate (gaki(glev+1),gbki(glev+1),stat=ierr)
        if (ierr/=0) stop "module_gfsdata: allocate vertical vars error"
        if (glev==60) then
! L60: hybrid coordinate
         gaki=(/ 0.00000,   .00000,   .02193,   .26557,   .97701,  2.26873 &
            ,    4.22363,  6.91977, 10.42730, 14.80485, 20.09558, 26.32304 &
            ,   33.48706, 41.55990, 50.48301, 60.16483, 70.47979, 81.26901 &
            ,   92.34269,103.48445,114.45724,125.01073,134.88958,143.84211 &
            ,  151.62865,158.02912,162.84922,165.92496,167.12546,166.35379 &
            ,  163.55255,158.81248,152.35664,144.43608,135.32210,125.30212 &
            ,  114.67568,103.75080, 92.84074, 82.26131, 72.32859, 63.32831 &
            ,   55.29569, 48.14158, 41.78099, 36.13468, 31.12936, 26.69766 &
            ,   22.77809, 19.31478, 16.25721, 13.55986, 11.18186,  9.08662 &
            ,    7.24145,  5.61725,  4.18811,  2.93106,  1.82570,   .85400 &
            ,    0.00000/)

!
         gbki=(/1.0000000,.99058760,.98124505,.96996497,.95697164,.94205290 &
            ,  .92502870,.90572942,.88400268,.85972113,.83279100,.80316112  &
            ,  .77083188,.73586353,.69838324,.65858995,.61675643,.57322786  &
            ,  .52841650,.48279246,.43687069,.39119502,.34632012,.30279276  &
            ,  .26113359,.22182071,.18527624,.15185637,.12184549,.09545413  &
            ,  .07281429,.05387320,.03841285,.02615862,.01679191,.00995912  &
            ,  .00528030,.00235711,.00077978,.00013304,.00000104,.00000000  &
            ,  .00000000,.00000000,.00000000,.00000000,.00000000,.00000000  &
            ,  .00000000,.00000000,.00000000,.00000000,.00000000,.00000000  &
            ,  .00000000,.00000000,.00000000,.00000000,.00000000,.00000000  &
            ,  0.00000000/)

        elseif (glev==72) then
        gaki=(/   0.00000,  0.00000,  0.00691,  0.06934,  0.25125,  0.60919,  &
                  1.19905,  2.07840,  3.30576,  4.93965,  7.03738,  9.65364,  &
                 12.83884, 16.63717, 21.08455, 26.20638, 32.01524, 38.50863,  &
                 45.66689, 53.45145, 61.80350, 70.64337, 79.87055, 89.36473,  &
                 98.98758,108.58553,117.99331,127.03812,135.54422,143.33774,  &
                150.25131,156.12838,160.82692,164.22232,166.20951,166.70406,  &
                165.64246,162.99021,158.82621,153.30553,146.59965,138.89340,  &
                130.38263,121.27200,111.77279,102.10102, 92.47566, 83.11719,  &
                 74.24635, 66.07746, 58.68863, 52.02149, 46.01449, 40.60953,  &
                 35.75213, 31.39157, 27.48082, 23.97656, 20.83899, 18.03171,  &
                 15.52153, 13.27825, 11.27451,  9.48552,  7.88892,  6.46451,  &
                  5.19414,  4.06147,  3.05182,  2.15204,  1.35033,  0.63614,  &
                  0.00000 /)
         gbki=(/ 1.0000000,.99472349,.98885067,.98227564,.97487364,.96652416, &
                 .95710466,.94648889,.93454851,.92115521,.90618329,.88951283, &
                 .87103331,.85064780,.82827754,.80386683,.77738809,.74884679, &
                 .71828604,.68579044,.65148892,.61555624,.57821275,.53972232, &
                 .50038827,.46054736,.42056203,.38081121,.34168025,.30355047, &
                 .26678885,.23173872,.19871153,.16798052,.13977611,.11428347, &
                 .09164181,.07193703,.05511859,.04103882,.02951509,.02033548, &
                 .01326388,.00804493,.00440868,.00207486,.00075667,.00016402, &
                 .00000624,.00000000,.00000000,.00000000,.00000000,.00000000, &
                 .00000000,.00000000,.00000000,.00000000,.00000000,.00000000, &
                 .00000000,.00000000,.00000000,.00000000,.00000000,.00000000, &
                 .00000000,.00000000,.00000000,.00000000,.00000000,.00000000, &
                 .00000000 /)
        else
          print*,'Error: glev not supported, glev=',glev
          stop 88
        endif
! 
        return
        end subroutine config_gfsdata
!
        subroutine allocate_gfsdata(igrd1p,jgrd1p,ntotal)
          integer igrd1p,jgrd1p,ntotal
          integer lngrdp,ierr
          lngrdp=igrd1p*jgrd1p
          allocate(sgeo(igrd1p,jgrd1p),psfc(igrd1p,jgrd1p),temp(igrd1p,jgrd1p,glev),        &
     &       u(lngrdp,glev),v(lngrdp,glev),rqgb(igrd1p,jgrd1p,glev,ntotal),   &
     &       simk(igrd1p,jgrd1p),stat=ierr) 
!>> defined locally in read_data.F <<
!    &       spfh(lngrdp,glev),     &
!    &       clwr(lngrdp,glev),ozon(lngrdp,glev),                      &
!    &       slmsk(lngrdp),   &
!    &       tsea(lngrdp),smc(lngrdp,kmsoil),          &
!    &       sheleg(lngrdp),stc(lngrdp,kmsoil),stat=ierr) 
!>> defined locally in read_data.F <<
          if (ierr/=0) stop "module_gfsdata: allocate fail"

!>> defined together in rqgb(:,:,1:ntotal) <<
!         if (ntotal.gt.3) then
!           allocate(rain(lngrdp,glev),qice(lngrdp,glev),         &
!    &       snow(lngrdp,glev),grpl(lngrdp,glev),stat=ierr) 
!            rain=0.0
!            qice=0.0
!            snow=0.0
!            grpl=0.0
!         if (ierr/=0) stop "module_gfsdata: allocate fail"
!         endif
        return
        end subroutine

        subroutine deallocate_gfsdata
          integer ierr
          deallocate(sgeo,psfc,temp,u,v,rqgb,simk,stat=ierr)
!         deallocate(sgeo,psfc,temp,u,v,spfh,clwr,ozon,simk,stat=ierr)
!                    tsea,smc,sheleg,stc,slmsk,stat=ierr)
          if (ierr/=0) stop "module_gfsdata: deallocate fail"
!         if (ntotal.gt.3) then
!           deallocate(rain,qice,snow,grpl,stat=ierr)
!           if (ierr/=0) stop "module_gfsdata: deallocate fail"
!         endif
        return
        end subroutine

        subroutine allocate_climdata(lngrd)
          integer lngrd
          integer ierr
          allocate(tg3clm(lngrd),zorclm(lngrd),albclm(lngrd,4),         &
     &       vegclm(lngrd),vetclm(lngrd),sotclm(lngrd),alfclm(lngrd,2), &
     &       vmnclm(lngrd),vmxclm(lngrd),slpclm(lngrd),absclm(lngrd),   &
     &       stat=ierr) 
          if (ierr/=0) stop "module_gfsdata: allocate clima fail"
        return
        end subroutine

        subroutine deallocate_climdata
          integer ierr
          deallocate(tg3clm,zorclm,albclm,                              &
     &       vegclm,vetclm,sotclm,alfclm,                               &
     &       vmnclm,vmxclm,slpclm,absclm,                               &
     &       stat=ierr) 
          if (ierr/=0) stop "module_gfsdata: deallocate clima fail"
        return
        end subroutine

        subroutine allocate_gridata(lngrd)
          integer lngrd
          integer ierr
          allocate(fflat(lngrd),fflon(lngrd),rslm(lngrd),stat=ierr)
          if (ierr/=0) stop "module_gridata: allocate grid fail"
        return
        end subroutine

        subroutine deallocate_gridata
          integer ierr
          deallocate(fflat,fflon,rslm,stat=ierr)
          if (ierr/=0) stop "module_gridata: deallocate grid fail"
        return
        end subroutine

      end module module_gfsdata  
