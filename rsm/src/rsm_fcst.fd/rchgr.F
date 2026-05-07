#include <inp_defopt.h>
      subroutine rchgr(nslmsk,norog,nsig,          &   
     &            lsig,lmtn,lp2r,percmtn,          &
     &            lsst,dosfc)
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rchgr
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  regional change resolution for sigma file and surface file
!
! program history log:
!
! usage:    call rchgr(nin,nsfcin,nslmsk,norog,mtnvar,
!          1           nou,nsfcou,lsig,lmtn)
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
!
!   output argument list:
!
!   input files: see input argument
!
!   output files: see input argument
!
!   subprograms called: 
!      rread  rsetsig  maxmin  newsig  extrap  rwrite  
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
      use module_rinpbase, only :  &
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
!yj!#ifdef NONHYD
!yj!     &  ,ppo,tto,wwo                               &
!yj!#endif
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
#include <inp_paramodel.h>
!cc
      dimension ggridn(lngrd),slmsk(lngrd)
      dimension ggrido1(lngrd)
      dimension pgz(lngrd),fgz(lnwav)
!
!      real(kind=kind_io4) sgridil(lngrd,lsoil),sgridol(lngrd,lsoil)
      real(kind=kind_io4) sgridn(lngrd)
!      real(kind=kind_io4) sgridi1(lngrd),sgrido1(lngrd)
      real(kind=kind_io4) sgrido1(lngrd)
!
      logical lsig,lmtn,lp2r,lsst,dosfc
!
      print *,' +++++++++++++ start regional chgr ++++++++++ '
!  read in original initial data
!
!yj      print 1101,idate
!yj 1101 format(1h0,1x,'original initial data read in',/,                   &
!yj     &           3x,' ihour=',i2,' month=',i2,' day=',i2,  &
!yj     &              ' year=',i4)
!yj      print 1102,igrd1,jgrd1,lngrd,levs
 1102 format(1x,'output file resolution =',6i10)
      print 1102,igrd1,jgrd1,lngrd,levr
 1103 format(1x,'output file resolution =',6i10)
!
      print *,' --- check input --- '
      call maxmin(uu,lngrd,levs,1,levs,' uu ')
      call maxmin(vv,lngrd,levs,1,levs,' vv ')
      call maxmin(te,lngrd,levs,1,levs,' te ')
      call maxmin(rq(1,1,1),lngrd,levs,1,levs,' rq ')
      if(ntoz.eq.1)                                                      &                       
     & call maxmin(rq(1,1,2),lngrd,levs,1,levs,' oz ')       
      if(ncld.eq.1)                                                      &
     & call maxmin(rq(1,1,3),lngrd,levs,1,levs,' cw ')
      call maxmin( q,lngrd,    1,1,    1,'  q ')
      call maxmin(fm2,lngrd,    1,1,    1,'  fm2 ')
!
!  define new sigma layer/level ?
!
      if( lsig ) then
         call newsig(nsig,ci,si,del,sl,cl,rpi)
      else
         do 1110 k=1,levr
           sl(k) = slin(k)
           si(k) = siin(k)
1110     continue
         si(levr+1) = siin(levs+1)
         call sigset(ci, si, del, sl, cl, rpi, levr)
      endif
!
      if(lmtn) then
        rewind(norog)
        read(norog) sgrido1
        call sgl2ful(sgrido1,gzo,lngrd)
        print *,'gz of output resolution read in'
        call maxmin(gzo,lngrd,1,1,1,'  gzo ')
        call maxmin(terr,lngrd,1,1,1,'  terr ')
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
          call inp_plnini
          call inp_rftini
! get ptb
          do ij=1,lngrd
            if (lp2r) then
            pgz(ij)=gzo(ij)-terr(ij)
            else 
            pgz(ij)=gzo(ij)-gz(ij)
            endif
          enddo
          print *,'after pgz'
          call maxmin(pgz,lngrd,1,1,1,' before inp_gdtocc, pgz ')
! grid to coef back to grid
          call inp_gdtocc(pgz,fgz,1)
          call inp_cctogd(fgz,pgz,1)
          print *,'after cc gd pgz'
          call maxmin(pgz,lngrd,1,1,1,' after inp_cctogd, pgz ')
! get full
          do ij=1,lngrd
            if (lp2r) then
            gzo(ij)=pgz(ij)+terr(ij)
            else
            gzo(ij)=pgz(ij)+gz(ij)
            endif
          enddo
          print *,'gz of lateral boundar modified'
          call maxmin(gzo,lngrd,1,1,1,'  gzo ')
        endif
!
      else
        do 3000 n=1,lngrd
        gzo(n) = gz(n)
3000    continue
      endif
!ccccccccc
!
      print *,' ----- check grid value before sg2sg ---- '
      call maxmin(uu,lngrd,levs,1,levs,' uu ')
      call maxmin(vv,lngrd,levs,1,levs,' vv ')
      call maxmin(te,lngrd,levs,1,levs,' te ')
      call maxmin(rq(1,1,1),lngrd,levs,1,levs,' rq ')
      if(ntoz.eq.1)                                                      &                       
     & call maxmin(rq(1,1,2),lngrd,levs,1,levs,' oz ')
      if(ncld.eq.1)                                                      &
     & call maxmin(rq(1,1,3),lngrd,levs,1,levs,' cw ')
      call maxmin( q,lngrd,    1,1,    1,'  q ')
!yj!#ifdef NONHYD
!yj!      call maxmin(pp,lngrd,levs,1,levs,' pp ')
!yj!      call maxmin(tt,lngrd,levs,1,levs,' tt ')
!yj!      call maxmin(ww,lngrd,levs+1,1,levs+1,' ww ')
!yj!#endif
!
      if( lsig .or. lmtn ) then
          call extrap(lsig,lmtn)
      else
        do 3010 n=1,lngrd
        qo(n) = q(n)
3010    continue
        do 3020 k=1,levr
        do 3020 n=1,lngrd
        uuo(n,k) = uu(n,k)
        vvo(n,k) = vv(n,k)
        teo(n,k) = te(n,k)
!yj!#ifdef NONHYD
!yj!        ppo(n,k) = pp(n,k)
!yj!        tto(n,k) = tt(n,k)
!yj!#endif
3020    continue
!yj!#ifdef NONHYD
!yj!        do k=1,levr
!yj!          do n=1,lngrd
!yj!            wwo(n,k) = ww(n,k)
!yj!          enddo
!yj!        enddo
!yj!#endif
        do 3030 nc=1,ntotal
        do 3030 k=1,levr
        do 3030 i=1,lngrd
        rqo(i,k,nc) = rq(i,k,nc)
3030    continue
      endif
!
      print *,' --- check output --- '
      call maxmin(uuo,lngrd,levr,1,levr,' uuo ')
      call maxmin(vvo,lngrd,levr,1,levr,' vvo ')
      call maxmin(teo,lngrd,levr,1,levr,' teo ')
      call maxmin(rqo(1,1,1),lngrd,levr,1,levr,' rqo ')
      if(ntoz.eq.1)                                                      &                  
     & call maxmin(rqo(1,1,2),lngrd,levr,1,levr,' ozo ')
!      if(ncld.eq.1)                                                      &
!     & call maxmin(rqo(1,1,3),lngrd,levr,1,levr,' cwo ')
      if(ncld.eq.1) then
       call maxmin(rqo(1,1,3),lngrd,levr,1,levr,' cwo ')
      else
       call maxmin(rqo(1,1,3),lngrd,levr,1,levr,' cwo ')
       call maxmin(rqo(1,1,4),lngrd,levr,1,levr,' rwo ')
       call maxmin(rqo(1,1,5),lngrd,levr,1,levr,' iwo ')
       call maxmin(rqo(1,1,6),lngrd,levr,1,levr,' swo ')
       call maxmin(rqo(1,1,7),lngrd,levr,1,levr,' gwo ')
      endif
      call maxmin( qo,lngrd,    1,1,    1,'  qo ')
      call maxmin( flat,lngrd,    1,1,    1,'  flat ')
      call maxmin( flon,lngrd,    1,1,    1,'  flon ')
!yj!#ifdef NONHYD
!yj!      call maxmin(ppo,lngrd,kdimo,1,levr,' ppo ')
!yj!      call maxmin(tto,lngrd,kdimo,1,levr,' tto ')
!yj!      call maxmin(wwo,lngrd,kdimo+1,1,levr+1,' wwo ')
!yj!#endif
!
!  addition for the change resolution of surface file
!
!yj      if( ls2r ) then
! ********************************************************************
      if (dosfc) then
      if (lsst) then
! read regional land-sea mask
         rewind nslmsk
         read(nslmsk) sgridn
         call sgl2ful(sgridn,slmsk,lngrd)
         rewind nslmsk
!
! use new sst on sea points
!
         call ana_sst(ggrido1,lngrd,flat,flon)
         do n=1,lngrd
           if( slmsk(n).eq.0.0 ) then
             tsea(n)= max( ggrido1(n), 271.21e0 )
           endif
         enddo
      endif   
         rewind nslmsk
         read (nslmsk,end=5115) sgridn
         rewind nslmsk
         call sgl2ful(sgridn,ggridn,lngrd)
         print *,' read input grid of slmask'
!         call maxmin(ggridn,lngrd,1,1,1,' slmask')
         do 5113 n=1,lngrd
           if( slmskin(n).le.1.0 ) slmskin(n)=0.0
           slmskin(n) = max(ggridn(n),slmskin(n))
5113     continue
5115     continue
      endif
! ******************************************************************
!endif (ls2r) below
!yj      endif
!
!
      print *,' ++++++ normal exit of regional chgr ++++ '
      return
!
      end
