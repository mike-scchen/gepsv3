      subroutine rclima(iy,im,id,ih,fh,len,lsoil,slmask,outlat,outlon,     & 
     &                 fntg3c,fnzorc,fnalbc,fnvegc,fnvetc,fnsotc,fnmskh,   & 
     &                 fnvmnc,fnvmxc,fnslpc,fnabsc,                        & 
     &                 tg3clm,zorclm,albclm,vegclm,vetclm,sotclm,alfclm,   &
     &                 vmnclm,vmxclm,slpclm,absclm)
!
      use machine , only : kind_io8,kind_io4
      use sfccyc_module
#include <defopt.h>
      implicit none
      real (kind=kind_io8) rjday,wei1x,wei2x,rjdayh,wei2m,wei1m,wei1s, &
     &                     wei2s,fh,stcmon1s,blto,blno,deltsfc,rjdayh2
      real (kind=kind_io8) wei1y,wei2y
      integer jdoy,jday,jh,jdow,mmm,mmp,mm,iret,monend,i,k,jm,jd,iy4, &
     &        jy,mon1,is2,isx,kpd9,is1,l,nn,mon2,mon,is, &
     &        im,id, &
     &        lugb,iy,len,lsoil,ih,iprnt,me,m1,m2,k1,k2, &
     &        jmsk,imsk,j,ialb,isot &
     &,       landice_cat
      integer nvet,nsot,nveg,nrec
!yj: use original albedo
      parameter(ialb=0)
#ifdef doNOAHLSM
!yj: use original soiltype isot=0, use new soiltype isot=1
      parameter(isot=1)
      parameter(nvet=97,nsot=98,nveg=99)
#else
#ifdef doPROGTM
      parameter(isot=0)
#endif
#endif
!
      real      slplmx
      parameter(slplmx=9.0)
!
      real (kind=kind_io8) zsca 
!
      character*500 fntsfc,fnwetc,fnsnoc,fnzorc,fnalbc,fnaisc, &
     &             fntg3c,fnscvc,fnsmcc,fnstcc,fnacnc,fnvegc, &
     &             fnvetc,fnsotc,fnalbc2 &
     &,            fnvmnc,fnvmxc,fnslpc,fnabsc,fnmskh
      real (kind=kind_io8) tsfclm(len),tsfcl2(len), &
     &     wetclm(len),snoclm(len), &
     &     zorclm(len),albclm(len,4),aisclm(len), &
     &     tg3clm(len),acnclm(len), &
     &     cvclm (len),cvbclm(len),cvtclm(len), &
     &     cnpclm(len), &
     &     smcclm(len,lsoil),stcclm(len,lsoil), &
     &     sliclm(len),scvclm(len),vegclm(len), &
     &     vetclm(len),sotclm(len),alfclm(len,2) &
     &,    vmnclm(len),vmxclm(len),slpclm(len),absclm(len)
!     real (kind=kind_io8) slmskh(5000*2500)
      real (kind=kind_io8) outlat(len), outlon(len)
!
      real (kind=kind_io8) slmask(len), tsfcl0(len)
      real (kind=kind_io8), allocatable :: slmask_noice(:)
!      
      real (kind=kind_io8) zorlmx,zorlmn,zoromx,zoromn,zorsmx,             &
     &                     zorsmn,zorimx,zorimn,zorjmx,zorjmn,             &
     &                     alblmx,alblmn,albomx,albomn,albsmx,             &
     &                     albsmn,albimx,albimn,albjmx,albjmn
      real (kind=kind_io8) abslmx,abssmx,abssmn,abslmn

      parameter(zorlmx=300.,zorlmn=2.,zoromx=1.0,zoromn=1.e-05,            &
     &          zorsmx=300.,zorsmn=2.,zorimx=1.0,zorimn=1.0,               &
     &          zorjmx=1.0,zorjmn=1.0)
      parameter(alblmx=0.80,alblmn=0.06,albomx=0.06,albomn=0.06,           &
     &          albsmx=0.80,albsmn=0.06,albimx=0.80,albimn=0.06,           &
     &          albjmx=0.80,albjmn=0.06)
!
      logical lanom, gaus, first
!
! set z0 based on sib vegetation type
      real (kind=kind_io8) z0_sib(13)
      data z0_sib /2.653, 0.826, 0.563, 1.089, 0.854, 0.856, &
     &             0.035, 0.238, 0.065, 0.076, 0.011, 0.125, &
     &             0.011 /
! set z0 based on igbp vegetation type
      real (kind=kind_io8) z0_igbp_min(20), z0_igbp_max(20)
      real (kind=kind_io8) z0_season(12)
      data z0_igbp_min /1.089, 2.653, 0.854, 0.826, 0.800, 0.050, &
     &                  0.030, 0.856, 0.856, 0.150, 0.040, 0.130, &
!    &                  1.000, 0.250, 0.011, 0.011, 0.001, 0.076, &
     &                  2.100, 0.250, 0.011, 0.011, 0.001, 0.076, &
     &                  0.050, 0.030/
      data z0_igbp_max /1.089, 2.653, 0.854, 0.826, 0.800, 0.050, &
     &                  0.030, 0.856, 0.856, 0.150, 0.040, 0.130, &
!    &                  1.000, 0.250, 0.011, 0.011, 0.001, 0.076, &
     &                  2.100, 0.250, 0.011, 0.011, 0.001, 0.076, &
     &                  0.050, 0.030/
!
! dayhf : julian day of the middle of each month
!
      real (kind=kind_io8) dayhf(13)
      data dayhf/ 15.5, 45.0, 74.5,105.0,135.5,166.0, &
     &           196.5,227.5,258.0,288.5,319.0,349.5,380.5/
!
      real (kind=kind_io8) fha(5)
      real(4) fha4(5)
      integer w3kindreal,w3kindint
      integer ida(8),jda(8),ivtyp, kpd7
!
      real (kind=kind_io8), allocatable :: tsf(:,:),sno(:,:), &
     &                     zor(:,:),wet(:,:), &
     &                     ais(:,:), acn(:,:),   scv(:,:), smc(:,:,:), &
     &                     tg3(:),   alb(:,:,:), alf(:,:), &
     &                     vet(:),   sot(:),     tsf2(:), &
     &                     veg(:,:), stc(:,:,:) &
     &,                    vmn(:), vmx(:),  slp(:), absm(:) &
     &,                    slmskh(:)
      integer mon1s, mon2s, sea1s, sea2s, sea1, sea2, hyr1, hyr2
      data first/.true./
      data mon1s/0/, mon2s/0/, sea1s/0/, sea2s/0/
!
      save first, tsf, sno, zor, wet,  ais, acn, scv, smc, tg3, &
     &     alb,   alf, vet, sot, tsf2, veg, stc, &
     &     vmn,   vmx, slp, absm, &
     &     mon1s, mon2s, sea1s, sea2s, dayhf, k1, k2, m1, m2, &
     &     landice_cat, slmskh
!
      logical lprnt
!
!yj-------------- >>
! deltsfc: frequency for surface data cycling (hour)
      deltsfc=24.0  !not involved in the model
!yj---------------<<      
      if (ialb.eq.1) then
        kpdabs = kpdabs_1
        kpdalb = kpdalb_1
!        alblmx = .99
!        albsmx = .99
!        alblmn = .01
!        albsmn = .01
        abslmx = 1.0
        abssmx = 1.0
        abssmn = .01
        abslmn = .01
      else
        kpdabs = kpdabs_0
        kpdalb = kpdalb_0
!        alblmx = .80
!        albsmx = .80
!        alblmn = .06
!        albsmn = .06
        abslmx = .80
        abssmx = .80
        abslmn = .01
        abssmn = .01
      endif
      if (isot.eq.1) then
        kpdsot=kpdsot_1
      else
        kpdsot=kpdsot_0
      endif
!        
      fntsfc=''
      fnwetc=''
      fnsnoc=''
      fnaisc=''
      fnscvc=''
      fnsmcc=''
      fnstcc=''
      fnacnc=''
!     fnvmnc=''
!     fnvmxc=''
!     fnslpc=''
!     fnabsc=''
!
! move from sfccycle
      lugb=101
      me=0
      if (first) then
      allocate (slmskh(5000*2500))
      call hmskrd(lugb,imsk,jmsk,fnmskh,                                   & 
     &              kpdmsk,slmskh,gaus,blno,blto,me)
      endif
!
      do i=1,len
        tsfclm(i) = 0.0
        tsfcl2(i) = 0.0
        snoclm(i) = 0.0
        wetclm(i) = 0.0
        zorclm(i) = 0.0
        aisclm(i) = 0.0
        tg3clm(i) = 0.0
        acnclm(i) = 0.0
        cvclm(i)  = 0.0
        cvbclm(i) = 0.0
        cvtclm(i) = 0.0
        cnpclm(i) = 0.0
        sliclm(i) = 0.0
        scvclm(i) = 0.0
        vmnclm(i) = 0.0
        vmxclm(i) = 0.0
        slpclm(i) = 0.0
        absclm(i) = 0.0
      enddo
      do k=1,lsoil
        do i=1,len
          smcclm(i,k) = 0.0
          stcclm(i,k) = 0.0
        enddo
      enddo
      do k=1,4
        do i=1,len
          albclm(i,k) = 0.0
        enddo
      enddo
      do k=1,2
        do i=1,len
          alfclm(i,k) = 0.0
        enddo
      enddo
!
      iret   = 0
      monend = 9999
!
      if (first) then
!
!    allocate variables to be saved
!
       allocate (tsf(len,2), sno(len,2),      zor(len,2), &
     &           wet(len,2), ais(len,2),      acn(len,2), &
     &           scv(len,2), smc(len,lsoil,2), &
     &           tg3(len),   alb(len,4,2),    alf(len,2), &
     &           vet(len),   sot(len), tsf2(len),  &
!clu [+1l] add vmn, vmx, slp, abs
     &           vmn(len),   vmx(len), slp(len), absm(len), &
     &           veg(len,2), stc(len,lsoil,2))
     
!
!     get tsf climatology for the begining of the forecast
!
        if (fh .gt. 0.0) then
!cbosu
          if (me == 0) print*,'bosu fh gt 0'

          iy4=iy
          if(iy.lt.101) iy4=1900+iy4
          fha=0
          ida=0
          jda=0
!         fha(2)=nint(fh)
          ida(1)=iy
          ida(2)=im
          ida(3)=id
          ida(5)=ih
          call w3kind(w3kindreal,w3kindint)
          if(w3kindreal == 4) then
            fha4=fha
            call w3movdat(fha4,ida,jda)
          else
            call w3movdat(fha,ida,jda)
          endif
          jy=jda(1)
          jm=jda(2)
          jd=jda(3)
          jh=jda(5)
          if (me .eq. 0) write(6,*) ' forecast jy,jm,jd,jh', &
     &                   jy,jm,jd,jh
          jdow = 0
          jdoy = 0
          jday = 0
          call w3doxdat(jda,jdow,jdoy,jday)
          rjday=jdoy+jda(5)/24.
          if(rjday.lt.dayhf(1)) rjday=rjday+365.
!
          if (me .eq. 0) write(6,*) 'forecast jy,jm,jd,jh=',jy,jm,jd,jh
!
!         for monthly mean climatology
!
          monend = 12
          do mm=1,monend
            mmm=mm
            mmp=mm+1
            if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
               mon1=mmm
               mon2=mmp
               go to 10
            endif
          enddo
          print *,'wrong rjday',rjday
          call abort
   10     continue
          wei1m = (dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
          wei2m = (rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
          if(mon2.eq.13) mon2=1
          if (me .eq. 0) print *,'rjday,mon1,mon2,wei1m,wei2m=', &
     &                   rjday,mon1,mon2,wei1m,wei2m
!
!       read monthly mean climatology of tsf
!
          kpd7 = -1
          do nn=1,2
            mon = mon1
            if (nn .eq. 2) mon = mon2
!!!            call fixrdc(lugb,fntsfc,kpdtsf,kpd7,mon,slmask, &
!!!     &                 tsf(1,nn),len,iret &
!!!     &,                imsk, jmsk, slmskh, gaus,blno, blto &
!!!     &,                outlat, outlon, me)
          enddo
!
!  tsf at the begining of forecast i.e. fh=0
!
!         do i=1,len
!           tsfcl0(i) = wei1m * tsf(i,1) + wei2m * tsf(i,2)
!         enddo
        endif
      endif
!
!  compute current jy,jm,jd,jh of forecast and the day of the year
!
      iy4=iy
      if(iy.lt.101) iy4=1900+iy4
      fha    = 0
      ida    = 0
      jda    = 0
      fha(2) = nint(fh)
      ida(1) = iy
      ida(2) = im
      ida(3) = id
      ida(5) = ih
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal==4) then
        fha4=fha
        call w3movdat(fha4,ida,jda)
      else
        call w3movdat(fha,ida,jda)
      endif
      jy     = jda(1)
      jm     = jda(2)
      jd     = jda(3)
      jh     = jda(5)
!     if (me .eq. 0) write(6,*) ' forecast jy,jm,jd,jh,rjday=', &
!    &               jy,jm,jd,jh,rjday
      jdow   = 0
      jdoy   = 0
      jday   = 0
      call w3doxdat(jda,jdow,jdoy,jday)
      rjday  = jdoy+jda(5)/24.
      if(rjday.lt.dayhf(1)) rjday=rjday+365.

      if (me .eq. 0) write(6,*) ' forecast jy,jm,jd,jh,rjday=', &
     &               jy,jm,jd,jh,rjday
!
      if (me .eq. 0) write(6,*) 'forecast jy,jm,jd,jh=',jy,jm,jd,jh
!
!     for monthly mean climatology
!
      monend = 12
      do mm=1,monend
         mmm=mm
         mmp=mm+1
         if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
            mon1=mmm
            mon2=mmp
            go to 20
         endif
      enddo
      print *,'wrong rjday',rjday
      call abort
   20 continue
      wei1m=(dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
      wei2m=(rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
      if(mon2.eq.13) mon2=1
      if (me .eq. 0) print *,'rjday,mon1,mon2,wei1m,wei2m=', &
     &               rjday,mon1,mon2,wei1m,wei2m
!
!     for seasonal mean climatology
!
      monend = 4
      is     = im/3 + 1
      if (is.eq.5) is = 1
      do mm=1,monend
        mmm = mm*3 - 2
        mmp = (mm+1)*3 - 2
        if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
          sea1 = mmm
          sea2 = mmp
          go to 30
        endif
      enddo
      print *,'wrong rjday',rjday
      call abort
   30 continue
      wei1s = (dayhf(sea2)-rjday)/(dayhf(sea2)-dayhf(sea1))
      wei2s = (rjday-dayhf(sea1))/(dayhf(sea2)-dayhf(sea1))
      if(sea2.eq.13) sea2=1
      if (me .eq. 0) print *,'rjday,sea1,sea2,wei1s,wei2s=', &
     &               rjday,sea1,sea2,wei1s,wei2s
!
!     for summer and winter values (maximum and minimum).
!
      monend = 2
      is     = im/6 + 1
      if (is.eq.3) is = 1
      do mm=1,monend
        mmm = mm*6 - 5
        mmp = (mm+1)*6 - 5
        if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
          hyr1 = mmm
          hyr2 = mmp
          go to 31
        endif
      enddo
      print *,'wrong rjday',rjday
      call abort
   31 continue
      wei1y = (dayhf(hyr2)-rjday)/(dayhf(hyr2)-dayhf(hyr1))
      wei2y = (rjday-dayhf(hyr1))/(dayhf(hyr2)-dayhf(hyr1))
      if(hyr2.eq.13) hyr2=1
      if (me .eq. 0) print *,'rjday,hyr1,hyr2,wei1y,wei2y=', &
     &               rjday,hyr1,hyr2,wei1y,wei2y
!
!  start reading in climatology and interpolate to the date
!
      first_time : if (first) then
!cbosu
        if (me == 0) print*,'bosu first time thru'
!
!     annual mean climatology
!
!  fraction of vegetation field for albedo --  there are two
!  fraction fields in this version: strong zenith angle dependent
!  and weak zenith angle dependent
!
        kpd9 = -1
!cjfe
        alf=0.
!cjfe

        kpd7=-1
        if (ialb.eq.1) then
!cbosu    still need facsf and facwf.  read them from the production
!cbosu    file
!cbosu
          call fixrdc(lugb,fnalbc2,kpdalf(1),kpd7,kpd9,slmask &
     &,               alf,len,iret &
     &,               imsk, jmsk, slmskh, gaus,blno, blto &
     &,               outlat, outlon, me)
        else
          call fixrdc(lugb,fnalbc,kpdalf(1),kpd7,kpd9,slmask &
     &,               alf,len,iret &
     &,               imsk, jmsk, slmskh, gaus,blno, blto &
     &,               outlat, outlon, me)
        endif
        do i = 1, len
          if(slmask(i).eq.1.) then
            alf(i,2) = 100. - alf(i,1)
          endif
        enddo
!
!  deep soil temperature
!
        if(fntg3c(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fntg3c,kpdtg3,kpd7,kpd9,slmask, &
     &                tg3,len,iret &
     &,               imsk, jmsk, slmskh, gaus,blno, blto &
     &,               outlat, outlon, me)
        endif
!
!  vegetation type
!
!  when using the new gldas soil moisture climatology, a veg type
!  dataset must be selected.
!
        if(fnvetc(1:8).ne.'        ') then
          if(fnvetc(1:15).eq.'out_landuse.bin')then
          if (me .eq. 0) write(6,*) 'climatological vegetation', &
     &                         '(reproduced data) type read in.' 
          inquire(iolength=nrec) vet
          open(nvet,file=fnvetc(1:15),access='direct',form='unformatted', &
              status='unknown',convert='big_endian',recl=nrec)

          read(nvet,rec=1) vet
          close(nvet)
          else
          kpd7=-1
          call fixrdc(lugb,fnvetc,kpdvet,kpd7,kpd9,slmask, &
     &                vet,len,iret &
     &,               imsk, jmsk, slmskh, gaus,blno, blto &
     &,               outlat, outlon, me)

          if (me .eq. 0) write(6,*) 'climatological vegetation', &
     &                              ' type read in.'
          endif
          landice_cat=13
          if (maxval(vet)> 13.0) landice_cat=15
        elseif(index(fnsmcc,'soilmgldas').ne.0) then ! new soil moisture climo
          if (me .eq. 0) write(6,*) 'fatal error: must choose'
          if (me .eq. 0) write(6,*) 'climatological veg type when'
          if (me .eq. 0) write(6,*) 'using new gldas soil moisture.'
          call abort
        endif
!
!  soil type
!
        if(fnsotc(1:8).ne.'        ') then
          if(fnsotc(1:16).eq.'out_soiltype.bin')then
          if (me .eq. 0) write(6,*) 'climatological soil ', &
     &                         '(reproduced data) type read in.' 
          inquire(iolength=nrec) sot
          open(nsot,file=fnsotc(1:16),access='direct',form='unformatted', &
              status='unknown',convert='big_endian',recl=nrec)
          read(nsot,rec=1) sot
          close(nsot)
          else
          kpd7=-1
          call fixrdc(lugb,fnsotc,kpdsot,kpd7,kpd9,slmask, &
     &                sot,len,iret &
     &,               imsk, jmsk, slmskh, gaus,blno, blto &
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological soil type read in.'
          endif
        endif

!
!  min vegetation cover
!
        if(fnvmnc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnvmnc,kpdvmn,kpd7,kpd9,slmask, &
     &                vmn,len,iret &
     &,               imsk, jmsk, slmskh, gaus,blno, blto &
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological shdmin read in.'
        endif
!
!  max vegetation cover
!
        if(fnvmxc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnvmxc,kpdvmx,kpd7,kpd9,slmask, &
     &                vmx,len,iret &
     &,               imsk, jmsk, slmskh, gaus,blno, blto &
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological shdmax read in.'
        endif
!
!  slope type
!
        if(fnslpc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnslpc,kpdslp,kpd7,kpd9,slmask, &
     &                slp,len,iret &
     &,               imsk, jmsk, slmskh, gaus,blno, blto &
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological slope read in.'
        endif
!
!  max snow albeod
!
        if(fnabsc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnabsc,kpdabs,kpd7,kpd9,slmask, &
     &                absm,len,iret &
     &,               imsk, jmsk, slmskh, gaus,blno, blto &
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological snoalb read in.'
        endif
!clu ----------------------------------------------------------------------
!
        is1 = sea1/3 + 1
        is2 = sea2/3 + 1
        if (is1 .eq. 5) is1 = 1
        if (is2 .eq. 5) is2 = 1
        do nn=1,2
!
!     seasonal mean climatology
          if(nn.eq.1) then
             isx=is1
          else
             isx=is2
          endif
          if(isx.eq.1) kpd9 = 12
          if(isx.eq.2) kpd9 = 3
          if(isx.eq.3) kpd9 = 6
          if(isx.eq.4) kpd9 = 9
!
!         seasonal mean climatology
!
!  albedo
!  there are four albedo fields in this version:
!  two for strong zeneith angle dependent (visible and near ir)
!  and two for weak zeneith angle dependent (vis ans nir)
!
          if (ialb.eq.0) then
            kpd7=-1
            do k = 1, 4
              call fixrdc(lugb,fnalbc,kpdalb(k),kpd7,kpd9,slmask, &
     &                    alb(1,k,nn),len,iret &
     &,                   imsk, jmsk, slmskh, gaus,blno, blto &
     &,                   outlat, outlon, me)
            enddo
          endif
!
!         monthly mean climatology
!
          mon = mon1
          if (nn .eq. 2) mon = mon2
!cbosu
!cbosu  new snowfree albedo database is monthly.
          if (ialb.eq.1) then
            kpd7=-1
            do k = 1, 4
              call fixrdc(lugb,fnalbc,kpdalb(k),kpd7,mon,slmask, &
     &                    alb(1,k,nn),len,iret &
     &,                   imsk, jmsk, slmskh, gaus,blno, blto &
     &,                   outlat, outlon, me)
            enddo
          endif

!     if(lprnt) print *,' mon1=',mon1,' mon2=',mon2
!
!  tsf at the current time t
!
!!!          kpd7=-1
!!!          call fixrdc(lugb,fntsfc,kpdtsf,kpd7,mon,slmask, &
!!!     &               tsf(1,nn),len,iret &
!!!     &,              imsk, jmsk, slmskh, gaus,blno, blto &
!!!     &,              outlat, outlon, me)
!     if(lprnt) print *,' tsf=',tsf(iprnt,nn),' nn=',nn
!
!  tsf...at time t-deltsfc
!
!     fh2 = fh - deltsfc
!     if (fh2 .gt. 0.0) then
!       call fixrd(lugb,fntsfc,kpdtsf,lclim,slmask, &
!    &             iy,im,id,ih,fh2,tsfcl2,len,iret &
!    &,            imsk, jmsk, slmskh, gaus,blno, blto &
!    &,            outlat, outlon, me)
!     else
!       do i=1,len
!         tsfcl2(i) = tsfclm(i)
!       enddo
!     endif
!
!  soil wetness
!
          if(fnwetc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnwetc,kpdwet,kpd7,mon,slmask, &
     &                  wet(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
          elseif(fnsmcc(1:8).ne.'        ') then
            if (index(fnsmcc,'global_soilmcpc.1x1.grb').ne.0) then ! the old
!climo data
              kpd7=-1
              call fixrdc(lugb,fnsmcc,kpdsmc,kpd7,mon,slmask, &
     &                    smc(1,lsoil,nn),len,iret &
     &,                   imsk, jmsk, slmskh, gaus,blno, blto &
     &,                   outlat, outlon, me)
              do l=1,lsoil-1
                do i = 1, len
                 smc(i,l,nn) = smc(i,lsoil,nn)
                enddo
              enddo
            else  ! the new gldas data.  it does not have data defined at landice
                  ! points.  so for efficiency, don't have fixrdc try to
                  ! find a value at landice points as defined by the vet type
                  ! (vet).
              allocate(slmask_noice(len))
              slmask_noice=1.0
              do i = 1, len
                if (nint(vet(i)) < 1 .or. &
     &              nint(vet(i)) == landice_cat) then
                  slmask_noice(i) = 0.0
                endif
              enddo
              do k = 1, lsoil
                if (k==1) kpd7=10     ! 0_10 cm    (pds octs 11 and 12)
                if (k==2) kpd7=2600   ! 10_40 cm
                if (k==3) kpd7=10340  ! 40_100 cm
                if (k==4) kpd7=25800  ! 100_200 cm
                call fixrdc(lugb,fnsmcc,kpdsmc,kpd7,mon,slmask_noice, &
     &                      smc(1,k,nn),len,iret &
     &,                     imsk, jmsk, slmskh, gaus,blno, blto &
     &,                     outlat, outlon, me)
              enddo
              deallocate(slmask_noice)
            endif
          else
            write(6,*) 'climatological soil wetness file not given'
!!!            call abort
          endif
!
!  soil temperature
!
          if(fnstcc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnstcc,kpdstc,kpd7,mon,slmask, &
     &                  stc(1,lsoil,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
            do l=1,lsoil-1
              do i = 1, len
               stc(i,l,nn) = stc(i,lsoil,nn)
              enddo
            enddo
          endif
!
!  sea ice
!
          kpd7=-1
          if(fnacnc(1:8).ne.'        ') then
            call fixrdc(lugb,fnacnc,kpdacn,kpd7,mon,slmask, &
     &                  acn(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
          elseif(fnaisc(1:8).ne.'        ') then
            call fixrdc(lugb,fnaisc,kpdais,kpd7,mon,slmask, &
     &                  ais(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
          else
            write(6,*) 'climatological ice cover file not given'
!!!            call abort
          endif
!
!  snow depth
!
          kpd7=-1
!!!          call fixrdc(lugb,fnsnoc,kpdsno,kpd7,mon,slmask, &
!!!     &                sno(1,nn),len,iret &
!!!     &,               imsk, jmsk, slmskh, gaus,blno, blto &
!!!     &,               outlat, outlon, me)
!
!  snow cover
!
          if(fnscvc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnscvc,kpdscv,kpd7,mon,slmask, &
     &                  scv(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
            write(6,*) 'climatological snow cover read in.'
          endif
!
!  surface roughness
!
      if(fnzorc(1:3) == 'sib') then
        if (me == 0) then
          write(6,*) 'roughness length to be set from sib veg type'
        endif
      elseif(fnzorc(1:4) == 'igbp') then
        if (me == 0) then
          write(6,*) 'roughness length to be set from igbp veg type'
        endif
      else
        kpd7=-1
        call fixrdc(lugb,fnzorc,kpdzor,kpd7,mon,slmask, &
     &              zor(1,nn),len,iret &
     &,             imsk, jmsk, slmskh, gaus,blno, blto &
     &,             outlat, outlon, me)
      endif
!
          do i = 1, len
!                           set clouds climatology to zero
            cvclm (i) = 0.
            cvbclm(i) = 0.
            cvtclm(i) = 0.
!
            cnpclm(i) = 0.  !set canopy water content climatology to zero
          enddo
!
!  vegetation cover
!
          if(fnvegc(1:8).ne.'        ') then
            if(fnvegc(1:15).eq.'out_vegfrac.bin')then
              if (me .eq. 0) write(6,*) 'climatological vegetation', &
     &               '(reproduced data) cover read in for mon=',mon
              inquire(iolength=nrec) veg(:,1)
              open(nveg,file=fnvegc(1:15),access='direct',form='unformatted', &
                 status='unknown',convert='big_endian',recl=nrec)
              read(nveg,rec=mon) veg(:,nn)
              close(nveg)
            else
            kpd7=-1
            call fixrdc(lugb,fnvegc,kpdveg,kpd7,mon,slmask, &
     &                  veg(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
            if (me .eq. 0) write(6,*) 'climatological vegetation', &
     &                                ' cover read in for mon=',mon
            endif
          endif

        enddo
!
      mon1s = mon1 ; mon2s = mon2 ; sea1s = sea1 ; sea2s = sea2
!
      if (me .eq. 0) print *,' mon1s=',mon1s,' mon2s=',mon2s &
     &,' sea1s=',sea1s,' sea2s=',sea2s
!
        k1  = 1 ; k2 = 2
        m1  = 1 ; m2 = 2
!
        first = .false.
      endif first_time
!
!     to get tsf climatology at the previous call to sfccycle
!
!     if (fh-deltsfc >= 0.0) then
        rjdayh = rjday - deltsfc/24.0
!     else
!       rjdayh = rjday
!     endif
!     if(lprnt) print *,' rjdayh=',rjdayh,' mon1=',mon1,' mon2=' &
!    &,mon2,' mon1s=',mon1s,' mon2s=',mon2s,' k1=',k1,' k2=',k2
      if (rjdayh .ge. dayhf(mon1)) then
        if (mon2 .eq. 1) mon2 = 13
        wei1x = (dayhf(mon2)-rjdayh)/(dayhf(mon2)-dayhf(mon1))
        wei2x = 1.0 - wei1x
        if (mon2 .eq. 13) mon2 = 1
      else
        rjdayh2 = rjdayh
        if (rjdayh .lt. dayhf(1)) rjdayh2 = rjdayh2 + 365.0
        if (mon1s .eq. mon1) then
          mon1s = mon1 - 1
          if (mon1s .eq. 0) mon1s = 12
          k2  = k1
          k1  = mod(k2,2) + 1
          mon = mon1s
!!!          kpd7=-1
!!!          call fixrdc(lugb,fntsfc,kpdtsf,kpd7,mon,slmask, &
!!!     &               tsf(1,k1),len,iret &
!!!     &,              imsk, jmsk, slmskh, gaus,blno, blto &
!!!     &,              outlat, outlon, me)
        endif
        mon2s = mon1s + 1
!       if (mon2s .eq. 1) mon2s = 13
        wei1x = (dayhf(mon2s)-rjdayh2)/(dayhf(mon2s)-dayhf(mon1s))
        wei2x = 1.0 - wei1x
        if (mon2s .eq. 13) mon2s = 1
!       do i=1,len
!         tsf2(i) = wei1x * tsf(i,k1) + wei2x * tsf(i,k2)
!       enddo
      endif
!
!cbosu new albedo is monthly
      if (sea1 .ne. sea1s) then
         sea1s = sea1
         sea2s = sea2
         m1    = mod(m1,2) + 1
         m2    = mod(m1,2) + 1
!
!     seasonal mean climatology
!
         isx   = sea2/3 + 1
         if (isx .eq. 5) isx = 1
         if(isx.eq.1) kpd9 = 12
         if(isx.eq.2) kpd9 = 3
         if(isx.eq.3) kpd9 = 6
         if(isx.eq.4) kpd9 = 9
!
!  albedo
!  there are four albedo fields in this version:
!  two for strong zeneith angle dependent (visible and near ir)
!  and two for weak zeneith angle dependent (vis ans nir)
!
!cbosu
        if (ialb.eq.0) then
           kpd7=-1
           do k = 1, 4
             call fixrdc(lugb,fnalbc,kpdalb(k),kpd7,kpd9,slmask &
     &,                 alb(1,k,m2),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
           enddo
        endif

      endif

      if (mon1 .ne. mon1s) then

         mon1s = mon1
         mon2s = mon2
         k1    = mod(k1,2) + 1
         k2    = mod(k1,2) + 1
!
!     monthly mean climatology
!
          mon = mon2
          nn  = k2
!cbosu
          if (ialb.eq.1) then
            if (me == 0) print*,'bosu 2nd time in clima for month ', &
     &                   mon, k1,k2
            kpd7=-1
            do k = 1, 4
              call fixrdc(lugb,fnalbc,kpdalb(k),kpd7,mon,slmask, &
     &                    alb(1,k,nn),len,iret &
     &,                   imsk, jmsk, slmskh, gaus,blno, blto &
     &,                   outlat, outlon, me)
            enddo
          endif
!
!  tsf at the current time t
!
!!!          kpd7=-1
!!!          call fixrdc(lugb,fntsfc,kpdtsf,kpd7,mon,slmask, &
!!!     &               tsf(1,nn),len,iret &
!!!     &,              imsk, jmsk, slmskh, gaus,blno, blto &
!!!     &,              outlat, outlon, me)
!
!  soil wetness
!
          if(fnwetc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnwetc,kpdwet,kpd7,mon,slmask, &
     &                  wet(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
          elseif(fnsmcc(1:8).ne.'        ') then
            if (index(fnsmcc,'global_soilmcpc.1x1.grb') /= 0) then ! the old
!climo data
              kpd7=-1
              call fixrdc(lugb,fnsmcc,kpdsmc,kpd7,mon,slmask, &
     &                    smc(1,lsoil,nn),len,iret &
     &,                   imsk, jmsk, slmskh, gaus,blno, blto &
     &,                   outlat, outlon, me)
              do l=1,lsoil-1
                do i = 1, len
                 smc(i,l,nn) = smc(i,lsoil,nn)
                enddo
              enddo
            else  ! the new gldas data.  it does not have data defined at landice
                  ! points.  so for efficiency, don't have fixrdc try to
                  ! find a value at landice points as defined by the vet type
                  ! (vet).
              allocate(slmask_noice(len))
              slmask_noice=1.0
              do i = 1, len
                if (nint(vet(i)) < 1 .or. &
     &              nint(vet(i)) == landice_cat) then
                  slmask_noice(i) = 0.0
                endif
              enddo
              do k = 1, lsoil
                if (k==1) kpd7=10     ! 0_10 cm   (pds octs 11 and 12)
                if (k==2) kpd7=2600   ! 10_40 cm
                if (k==3) kpd7=10340  ! 40_100 cm
                if (k==4) kpd7=25800  ! 100_200 cm
                call fixrdc(lugb,fnsmcc,kpdsmc,kpd7,mon,slmask_noice, &
     &                      smc(1,k,nn),len,iret &
     &,                     imsk, jmsk, slmskh, gaus,blno, blto &
     &,                     outlat, outlon, me)
              enddo
              deallocate(slmask_noice)
            endif
          else
            write(6,*) 'climatological soil wetness file not given'
!!!            call abort
          endif
!
!  sea ice
!
          kpd7=-1
          if(fnacnc(1:8).ne.'        ') then
            call fixrdc(lugb,fnacnc,kpdacn,kpd7,mon,slmask, &
     &                  acn(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
          elseif(fnaisc(1:8).ne.'        ') then
            call fixrdc(lugb,fnaisc,kpdais,kpd7,mon,slmask, &
     &                  ais(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
          else
            write(6,*) 'climatological ice cover file not given'
!!!            call abort
          endif
!
!  snow depth
!
!!!          kpd7=-1
!!!          call fixrdc(lugb,fnsnoc,kpdsno,kpd7,mon,slmask, &
!!!     &                sno(1,nn),len,iret &
!!!     &,               imsk, jmsk, slmskh, gaus,blno, blto &
!!!     &,               outlat, outlon, me)
!
!  snow cover
!
          if(fnscvc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnscvc,kpdscv,kpd7,mon,slmask, &
     &                  scv(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
            write(6,*) 'climatological snow cover read in.'
          endif
!
!  surface roughness
!
      if(fnzorc(1:3) == 'sib') then
        if (me == 0) then
          write(6,*) 'roughness length to be set from sib veg type'
        endif
      elseif(fnzorc(1:4) == 'igbp') then
        if (me == 0) then
          write(6,*) 'roughness length to be set from igbp veg type'
        endif
      else
        kpd7=-1
        call fixrdc(lugb,fnzorc,kpdzor,kpd7,mon,slmask, &
     &              zor(1,nn),len,iret &
     &,             imsk, jmsk, slmskh, gaus,blno, blto &
     &,             outlat, outlon, me)
      endif
!
!  vegetation cover
!
          if(fnvegc(1:8).ne.'        ') then
            if(fnvegc(1:15).eq.'out_vegfrac.bin')then
              if (me .eq. 0) write(6,*) 'climatological vegetation', &
     &               '(reproduced data) cover read in for mon=',mon
              inquire(iolength=nrec) veg(:,1)
              open(nveg,file=fnvegc(1:15),access='direct',form='unformatted', &
                 status='unknown',convert='big_endian',recl=nrec)
              read(nveg,rec=mon) veg(:,nn)
              close(nveg)
            else
            kpd7=-1
            call fixrdc(lugb,fnvegc,kpdveg,kpd7,mon,slmask, &
     &                  veg(1,nn),len,iret &
     &,                 imsk, jmsk, slmskh, gaus,blno, blto &
     &,                 outlat, outlon, me)
!           if (me .eq. 0) write(6,*) 'climatological vegetation', &
!    &                                ' cover read in for mon=',mon
            endif
          endif
!
      endif
!
!     now perform the time interpolation
!
! when chosen, set the z0 based on the vegetation type.
! for this option to work, namelist variable fnvetc must be
! set to point at the proper vegetation type file.
      if(fnzorc(1:3) == 'sib') then
        if(fnvetc(1:4) == '   ') then
          if (me==0) write(6,*) "must choose sib veg type climo file"
          call abort
        endif
        zorclm = 0.0
        do i=1,len
          ivtyp=nint(vet(i))
          if (ivtyp >= 1 .and. ivtyp <= 13) then
            zorclm(i) = z0_sib(ivtyp)
          endif
        enddo
      elseif(fnzorc(1:4) == 'igbp') then
        if(fnvetc(1:4) == '   ') then
          if (me==0) write(6,*) "must choose igbp veg type climo file"
          call abort
        endif
        zorclm = 0.0
        do i=1,len
          ivtyp=nint(vet(i))
          if (ivtyp >= 1 .and. ivtyp <= 20) then
            z0_season(1) = z0_igbp_min(ivtyp)
            z0_season(7) = z0_igbp_max(ivtyp)
            if(outlat(i) < 0.0)then
              zorclm(i) = wei1y * z0_season(hyr2) + &
     &                    wei2y *z0_season(hyr1)
             else
              zorclm(i) = wei1y * z0_season(hyr1) + &
     &                    wei2y *z0_season(hyr2)
           endif
          endif
        enddo
      else
        do i=1,len
          zorclm(i) = wei1m * zor(i,k1) + wei2m * zor(i,k2)
        enddo
      endif
!
      do i=1,len
!       tsfclm(i) = wei1m * tsf(i,k1) + wei2m * tsf(i,k2)
!       snoclm(i) = wei1m * sno(i,k1) + wei2m * sno(i,k2)
        cvclm(i)  = 0.0
        cvbclm(i) = 0.0
        cvtclm(i) = 0.0
        cnpclm(i) = 0.0
!       tsfcl2(i) = tsf2(i)
      enddo
!     if(lprnt) print *,' tsfclm=',tsfclm(iprnt),' wei1m=',wei1m &
!    &,' wei2m=',wei2m,' tsfk12=',tsf(iprnt,k1),tsf(iprnt,k2)
!
!     if (fh .eq. 0.0) then
!       do i=1,len
!         tsfcl0(i) = tsfclm(i)
!       enddo
!     endif
!     if (rjdayh .ge. dayhf(mon1)) then
!       do i=1,len
!         tsf2(i)   = wei1x * tsf(i,k1) + wei2x * tsf(i,k2)
!         tsfcl2(i) = tsf2(i)
!       enddo
!     endif
!     if(lprnt) print *,' tsf2=',tsf2(iprnt),' wei1x=',wei1x &
!    &,' wei2x=',wei2x,' tsfk12=',tsf(iprnt,k1),tsf(iprnt,k2) &
!    &,' mon1s=',mon1s,' mon2s=',mon2s &
!    &,' slmask=',slmask(iprnt)
!
      if(fnacnc(1:8).ne.'        ') then
        do i=1,len
          acnclm(i) = wei1m * acn(i,k1) + wei2m * acn(i,k2)
        enddo
      elseif(fnaisc(1:8).ne.'        ') then
        do i=1,len
          aisclm(i) = wei1m * ais(i,k1) + wei2m * ais(i,k2)
        enddo
      endif
!
      if(fnwetc(1:8).ne.'        ') then
        do i=1,len
          wetclm(i) = wei1m * wet(i,k1) + wei2m * wet(i,k2)
        enddo
      elseif(fnsmcc(1:8).ne.'        ') then
        do k=1,lsoil
          do i=1,len
            smcclm(i,k) = wei1m * smc(i,k,k1) + wei2m * smc(i,k,k2)
          enddo
        enddo
      endif
!
      if(fnscvc(1:8).ne.'        ') then
        do i=1,len
          scvclm(i) = wei1m * scv(i,k1) + wei2m * scv(i,k2)
        enddo
      endif
!
      if(fntg3c(1:8).ne.'        ') then
        do i=1,len
          tg3clm(i) =         tg3(i)
        enddo
      elseif(fnstcc(1:8).ne.'        ') then
        do k=1,lsoil
          do i=1,len
            stcclm(i,k) = wei1m * stc(i,k,k1) + wei2m * stc(i,k,k2)
          enddo
        enddo
      endif
!
      if(fnvegc(1:8).ne.'        ') then
        do i=1,len
          vegclm(i) = wei1m * veg(i,k1) + wei2m * veg(i,k2)
        enddo
      endif
!
      if(fnvetc(1:8).ne.'        ') then
        do i=1,len
          vetclm(i) =         vet(i)
        enddo
      endif
!
      if(fnsotc(1:8).ne.'        ') then
        do i=1,len
          sotclm(i) =         sot(i)
        enddo
      endif


!clu ----------------------------------------------------------------------
!
      if(fnvmnc(1:8).ne.'        ') then
        do i=1,len
          vmnclm(i) = vmn(i)
        enddo
      endif
!
      if(fnvmxc(1:8).ne.'        ') then
        do i=1,len
          vmxclm(i) = vmx(i)
        enddo
      endif
!
      if(fnslpc(1:8).ne.'        ') then
        do i=1,len
          slpclm(i) = slp(i)
        enddo
      endif
!
      if(fnabsc(1:8).ne.'        ') then
        do i=1,len
          absclm(i) = absm(i)
        enddo
      endif
!clu ----------------------------------------------------------------------
!
!cbosu  diagnostic print
      if (me == 0) print*,'monthly albedo weights are ', &
     &             wei1m,' for k', k1, wei2m, ' for k', k2

      if (ialb.eq.1) then
        do k=1,4
          do i=1,len
            albclm(i,k) = wei1m * alb(i,k,k1) + wei2m * alb(i,k,k2)
          enddo
        enddo
      else
        do k=1,4
          do i=1,len
            albclm(i,k) = wei1s * alb(i,k,m1) + wei2s * alb(i,k,m2)
          enddo
        enddo
      endif
!
      do k=1,2
        do i=1,len
          alfclm(i,k) = alf(i,k)
        enddo
      enddo
!
!  end of climatology reads
!
!!! move from sfccycle--->>>
!
!  scale surface roughness and albedo to model required units
!
      zsca=100.
      call scale(zorclm,len,zsca)
      zsca=0.01
      call scale(albclm,len,zsca)
      call scale(albclm(1,2),len,zsca)
      call scale(albclm(1,3),len,zsca)
      call scale(albclm(1,4),len,zsca)
      call scale(alfclm,len,zsca)
      call scale(alfclm(1,2),len,zsca)
!clu [+4l] scale vmn, vmx, abs from percent to fraction
      zsca=0.01
      call scale(vmnclm,len,zsca)
      call scale(vmxclm,len,zsca)
      call scale(absclm,len,zsca)

!
!  set albedo over ocean to albomx
!
      call albocn(albclm,slmask,albomx,len)
!
!  make sure vegetation type and soil type are non zero over land
!
      call landtyp(vetclm,sotclm,slpclm,slmask,len)
!!!
!!      albclm=albclm*0.01
!!      alfclm=alfclm*0.01
!!      zorclm=zorclm*100
!!!
!!      vmnclm=vmnclm*0.01
!!      vmxclm=vmxclm*0.01
!!      absclm=absclm*0.01
      do i=1,len
         if (slmask(i).eq.1) then
            if (zorclm(i).gt.zorlmx) zorclm(i)=zorlmx
            if (zorclm(i).lt.zorlmn) zorclm(i)=zorlmn
            do j=1,4
            if (albclm(i,j).gt.alblmx) albclm(i,j)=alblmx
            if (albclm(i,j).lt.alblmn) albclm(i,j)=alblmn
            enddo
         else
            if (zorclm(i).gt.zoromx) zorclm(i)=zoromx
            if (zorclm(i).lt.zoromn) zorclm(i)=zoromn
            do j=1,4
            if (albclm(i,j).gt.albomx) albclm(i,j)=albomx
            if (albclm(i,j).lt.albomn) albclm(i,j)=albomn
            enddo
         endif
         if (slmask(i).eq.2) then
           slpclm(i)=slplmx
         endif
      enddo

!!! move from sfccycle---<<<
!!!
      return
      end subroutine rclima
