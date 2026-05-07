#include <inp_defopt.h>
      subroutine int_rchgr(nslmsk,norog,nsig,          &   
     &                 lsig,lmtn,lp2r,percmtn,         &
     &                 lsst)  
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  int_rchgr
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  regional change resolution for sigma file and surface file
!
! program history log:
!
! usage:    call int_rchgr(nin,nsfcin,nslmsk,norog,mtnvar,
!          1           nou,nsfcou,lsig,lmtn,lhor)
!   input argument list:
!     nin       - unit number for sigma file input
!     nsfcin    - unit number for surface file input
!     nslmsk    - unit number for sea-land mask file input
!     norog     - unit number for terrain file input
!     mtnvar    - unit number for mountain variance file input
!     nou       - unit number for sigma file output
!     nsfcou    - unit number for surface file output
!     lsig      - logical index of using new sigma layes
!     lmtn      - logical index of using new mountain
!     lhor      - logical index of horizontal interpolation
!
!   output argument list:
!
!   input files: see input argument
!
!   output files: see input argument
!
!   subprograms called: 
!      int_rread  rsetsig  maxmin  int_newsig  int_extrap  rwrite  
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
      use machine , only : kind_io4
!yj add module start -----------------------      
      use module_rinpinit, only :  &
!comrinp1      
!yj     &  idate &
!comrinp2     
     & q,te,uu,vv,rq,gz &
!yj!#ifdef NONHYD 
!yj!     & ,pp,tt,ww &
!yj!#endif
     & ,flat,flon,fm2,fm2x,fm2y    &
!comrinp4
     & ,qo,teo,uuo,vvo,rqo,gzo                    & 
#ifdef NONHYD
     &  ,ppo,tto,wwo                               &
#endif
!inpver     
     &  ,siin,slin,delin,ciin,clin,rpiin &
!comvet     
     &  ,si,sl,del,ci,cl,rpi    &
!llgz     
     &  ,terr  &
!comrinpsfc     
     &  ,tsea  &
     &  ,slmskin=>slmsk

!yj add module   end -----------------------      

!ccj
#include <inpint_paramodel.h>
!cc
!
! change resolution in regional domain
!
      parameter(idimi=igrd1i,jdimi=jgrd1i,ijdimi=idimi*jdimi,            &                 
     &          kdimi=levri,kdimpi=kdimi+1,kdimmi=kdimi-1)
      parameter(idimo=igrd1o,jdimo=jgrd1o,ijdimo=idimo*jdimo,            &
     &          kdimo=levro,kdimpo=kdimo+1,kdimmo=kdimo-1)
!
      dimension ggridn(ijdimo),slmsk(ijdimo)
      dimension ggrido1(ijdimo)
      dimension pgz(ijdimo),fgz(lnwav)
!
!      real(kind=kind_io4) sgridil(ijdimi,lsoil),sgridol(ijdimo,lsoil)
      real(kind=kind_io4) sgridn(ijdimo)
!      real(kind=kind_io4) sgridi1(ijdimi),sgrido1(ijdimo)
      real(kind=kind_io4) sgrido1(ijdimo)
!
!yj      integer version
!
      logical lsig,lmtn,lp2r,lsst
!yj      logical lsig,lmtn,lhor,lc2r,ls2r,lp2r,lsst
!yj                         X  , yes, yes,    , no
!
      print *,' +++++++++++++ start regional chgr ++++++++++ '
!yj      if(lhor) then
!yj        print *,' +++ change horizontal resolution ++ '
!yj        print *,' +++ have not yet implemented!!!! ++ '
!yj        print *,' +++ abort!!!!!!                  ++ '
!yj        call abort
!yj      endif
!
!  read in original initial data
!
!yj      if( lc2r ) then
! ----------------------------------------------------------------------
!yj      print 1101,fhour,idate
!yj 1101 format(1h0,1x,'original initial data read in',/,                   &
!yj     &           3x,'fhour=',f6.1,' ihour=',i2,' month=',i2,' day=',i2,  &
!yj     &              ' year=',i4)
      print 1102,idimi,jdimi,ijdimi,kdimi
 1102 format(1x,'output file resolution =',6i10)
      print 1102,idimo,jdimo,ijdimo,kdimo
 1103 format(1x,'output file resolution =',6i10)
!
!      call sigset(ciin, siin, delin, slin, clin, rpiin, kdimi)
!
      print *,' --- check input --- '
      call maxmin(uu,ijdimi,kdimi,1,kdimi,' uu ')
      call maxmin(vv,ijdimi,kdimi,1,kdimi,' vv ')
      call maxmin(te,ijdimi,kdimi,1,kdimi,' te ')
      call maxmin(rq(1,1,1),ijdimi,kdimi,1,kdimi,' rq ')
      if(ntoz.eq.1)                                                      &                       
     & call maxmin(rq(1,1,2),ijdimi,kdimi,1,kdimi,' oz ')       
      if(ncld.eq.1)                                                      &
     & call maxmin(rq(1,1,3),ijdimi,kdimi,1,kdimi,' cw ')
      call maxmin( q,ijdimi,    1,1,    1,'  q ')
      call maxmin(fm2,ijdimi,    1,1,    1,'  fm2 ')
!
!  define new sigma layer/level ?
!
      if( lsig ) then
         call int_newsig(nsig,ci,si,del,sl,cl,rpi)
      else
         do 1110 k=1,kdimo
           sl(k) = slin(k)
           si(k) = siin(k)
1110     continue
         si(kdimo+1) = siin(kdimi+1)
         call sigset(ci, si, del, sl, cl, rpi, kdimo)
      endif
!
      if(lmtn) then
        rewind(norog)
        read(norog) sgrido1
        call sgl2ful(sgrido1,gzo,ijdimo)
        print *,'gz of output resolution read in'
        call maxmin(gzo,ijdimo,1,1,1,'  gzo ')
        call maxmin(terr,ijdimo,1,1,1,'  terr ')
        if(percmtn.gt.0.0) then
! modify lateral boundary mountain
          xc = float(igrd)/2.0 + 1.0  
          yc = float(jgrd)/2.0 + 1.0 
          xs = percmtn*xc
          ys = percmtn*yc
          xc1= xc - xs - 1.0            
          yc1= yc - ys - 1.0 
          pow=0.5
!         pow=1.0
          pi=2.*asin(1.0)
          do lat = 1,jgrd1
            yls= min( 1.0, max( abs(lat-yc)-yc1,0.0 )/ys )
            yls= (1.+sin(pi*(yls-0.5)))/2.
            jlats=(lat-1)*igrd1
            do i=1,igrd1      
              xl = min( 1.0, max( abs(i   -xc)-xc1,0.0 )/xs )
              xl = (1.+sin(pi*(xl-0.5)))/2.
              coes = max(xl,yls)
              if( coes.lt.1.e-6 ) then
                coes = 0.0
              else
                coes = coes ** pow
              endif
              if (lp2r) then
              gzo(i+jlats)=terr(i+jlats)*coes+gzo(i+jlats)*(1.-coes)
              else
              gzo(i+jlats)=gz(i+jlats)*coes+gzo(i+jlats)*(1.-coes)
              endif
            enddo
          enddo
! modify gzo with regional spectral smooth
! init spectral constant
          call inpint_plnini
          call inpint_rftini
! get ptb
          do ij=1,lngrd
            if (lp2r) then
            pgz(ij)=gzo(ij)-terr(ij)
            else 
            pgz(ij)=gzo(ij)-gz(ij)
            endif
          enddo
          print *,'after pgz'
          call maxmin(pgz,ijdimo,1,1,1,' before inpint_gdtocc, pgz ')
! grid to coef back to grid
          call inpint_gdtocc(pgz,fgz,1)
          call inpint_cctogd(fgz,pgz,1)
          print *,'after cc gd pgz'
          call maxmin(pgz,ijdimo,1,1,1,' after inpint_cctogd, pgz ')
! get full
          do ij=1,lngrd
            if (lp2r) then
            gzo(ij)=pgz(ij)+terr(ij)
            else
            gzo(ij)=pgz(ij)+gz(ij)
            endif
          enddo
          print *,'gz of lateral boundar modified'
          call maxmin(gzo,ijdimo,1,1,1,'  gzo ')
        endif
!
      else
        do 3000 n=1,ijdimo
        gzo(n) = gz(n)
3000    continue
      endif
!ccccccccc
!
      print *,' ----- check grid value before sg2sg ---- '
      call maxmin(uu,ijdimo,kdimi,1,kdimi,' uu ')
      call maxmin(vv,ijdimo,kdimi,1,kdimi,' vv ')
      call maxmin(te,ijdimo,kdimi,1,kdimi,' te ')
      call maxmin(rq(1,1,1),ijdimo,kdimi,1,kdimi,' rq ')
      if(ntoz.eq.1)                                                      &                       
     & call maxmin(rq(1,1,2),ijdimo,kdimi,1,kdimi,' oz ')
      if(ncld.eq.1)                                                      &
     & call maxmin(rq(1,1,3),ijdimo,kdimi,1,kdimi,' cw ')
      call maxmin( q,ijdimo,    1,1,    1,'  q ')
!yj!#ifdef NONHYD
!yj!      call maxmin(pp,ijdimo,kdimi,1,kdimi,' pp ')
!yj!      call maxmin(tt,ijdimo,kdimi,1,kdimi,' tt ')
!yj!      call maxmin(ww,ijdimo,kdimi+1,1,kdimi+1,' ww ')
!yj!#endif
!
      if( lsig .or. lmtn ) then
          call int_extrap(lsig,lmtn)
      else
        do 3010 n=1,ijdimo
        qo(n) = q(n)
3010    continue
        do 3020 k=1,kdimo
        do 3020 n=1,ijdimo
        uuo(n,k) = uu(n,k)
        vvo(n,k) = vv(n,k)
        teo(n,k) = te(n,k)
#ifdef NONHYD
!yj!        ppo(n,k) = pp(n,k)
!yj!        tto(n,k) = tt(n,k)
        ppo(n,k) = qo(n)+log(sl(k))
        tto(n,k) = teo(n,k)
#endif
3020    continue
#ifdef NONHYD
        do k=1,kdimo+1
          do n=1,ijdimo
            wwo(n,k) = 0.0
          enddo
        enddo
#endif
        do 3030 nc=1,ntotal
        do 3030 k=1,kdimo
        do 3030 i=1,ijdimo
        rqo(i,k,nc) = rq(i,k,nc)
3030    continue
      endif
        where (rqo<0.0) rqo=0.0
!
      print *,' --- check output --- '
      call maxmin(uuo,ijdimo,kdimo,1,kdimo,' uuo ')
      call maxmin(vvo,ijdimo,kdimo,1,kdimo,' vvo ')
      call maxmin(teo,ijdimo,kdimo,1,kdimo,' teo ')
      call maxmin(rqo(1,1,1),ijdimo,kdimo,1,kdimo,' rqo ')
      if(ntoz.eq.1)                                                      &                  
     & call maxmin(rqo(1,1,2),ijdimo,kdimo,1,kdimo,' ozo ')
!      if(ncld.eq.1)                                                      &
!     & call maxmin(rqo(1,1,3),ijdimo,kdimo,1,kdimo,' cwo ')
      if(ncld.eq.1) then
       call maxmin(rqo(1,1,3),ijdimo,kdimo,1,kdimo,' cwo ')
      else
       call maxmin(rqo(1,1,3),ijdimo,kdimo,1,kdimo,' cwo ')
       call maxmin(rqo(1,1,4),ijdimo,kdimo,1,kdimo,' rwo ')
       call maxmin(rqo(1,1,5),ijdimo,kdimo,1,kdimo,' iwo ')
       call maxmin(rqo(1,1,6),ijdimo,kdimo,1,kdimo,' swo ')
       call maxmin(rqo(1,1,7),ijdimo,kdimo,1,kdimo,' gwo ')
      endif
      call maxmin( qo,ijdimo,    1,1,    1,'  qo ')
      call maxmin( flat,ijdimo,    1,1,    1,'  flat ')
      call maxmin( flon,ijdimo,    1,1,    1,'  flon ')
#ifdef NONHYD
      call maxmin(ppo,ijdimo,kdimo,1,kdimo,' ppo ')
      call maxmin(tto,ijdimo,kdimo,1,kdimo,' tto ')
      call maxmin(wwo,ijdimo,kdimo+1,1,kdimo+1,' wwo ')
#endif
!
! --------------------------------------------------------------------
!yj      endif  !(lc2r)
!
!  addition for the change resolution of surface file
!
!yj      if( ls2r ) then
! ********************************************************************
      if (lsst) then
! read regional land-sea mask
         rewind nslmsk
         read(nslmsk) sgridn
         call sgl2ful(sgridn,slmsk,ijdimo)
         rewind nslmsk
!
! use new sst on sea points
!
         call ana_sst(ggrido1,ijdimo,flat,flon)
         do n=1,ijdimo
           if( slmsk(n).eq.0.0 ) then
             tsea(n)= max( ggrido1(n), 271.21e0 )
           endif
         enddo
      endif   
         rewind nslmsk
         read (nslmsk,end=5115) sgridn
         rewind nslmsk
         call sgl2ful(sgridn,ggridn,ijdimo)
         print *,' read input grid of slmask'
!         call maxmin(ggridn,ijdimo,1,1,1,' slmask')
         do 5113 n=1,ijdimo
           if( slmskin(n).le.1.0 ) slmskin(n)=0.0
           slmskin(n) = max(ggridn(n),slmskin(n))
5113     continue
5115     continue
! ******************************************************************
!endif (ls2r) below
!yj      endif
!
!
      print *,' ++++++ normal exit of regional chgr ++++ '
      return
!
      end
