      subroutine rdeldif
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rdeldif
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  horizontal diffusion of regional perturbation waves
!            for regional spectral model.
!
! program history log:
!
! usage:    call rdeldif
!   input argument list:
!
!   output argument list:
!
!   common block :
!       rscomf
!
!   input files: none
!
!   output files: none
!
!   subprograms called: none
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
!.....

      use physcons,only: con_pi
!ccj
#include <paramodel.h>
#ifdef MP
#define LNWAVS lnwavp
#else
#define LNWAVS lnwav
#endif
#include <rscomspc.h>
#include <rscomio.h>
#include <rscomver.h>

#ifdef DIFDZ
!yj
      common /repron/                                                    &
     &  tmp(LNWAVS,levr+1),tmp1(LNWAVS,levr+1)                           &
     &, dim(LNWAVS,levr)
#endif
!cc
      dimension xkm(LNWAVS),xkh(LNWAVS)
      dimension rkm(LNWAVS),rkh(LNWAVS)
#ifdef DIVDIF
      dimension xkd(LNWAVS), rkd(LNWAVS)
#endif

#ifdef DIVDIF
#ifndef TSDIVDIF
#define TSDIVDIF 3.0
#endif
#endif

#ifdef UPDIFF
#ifndef NHDIFORD
#define NHDIFORD 2
#endif
#ifndef LEVTOP
#define LEVTOP 2
#endif
!yj:top 3 levels (<18hPa) use 2nd diffusion, others 4th
      dimension rkmt(LNWAVS),rkht(LNWAVS),rkmt2(LNWAVS),rkht2(LNWAVS)
!      &    , xkmtt(LNWAVS), xkhtt(LNWAVS) &
!     &         ,rkmtt(LNWAVS), rkhtt(LNWAVS) 
#endif
      dimension rtrd(levr)
!x    dimension tt(LNWAVS),ct(levr),cq(levr)
#ifdef UPDIFF
!      save ddt,xkm,xkh,rkm,rkh,rtrd
#else
      save ddt,xkm,xkh,rkm,rkh,rtrd
#ifdef DIVDIF
      save xkd, rkd
#endif
#endif
!.................................................................
!
#ifdef UPDIFF
      ddt = deltim
#else
      data ddt/0/
      if(ddt.ne.deltim) then
         ddt = deltim
         coefmax = (epsxmax**2+epsymax**2)**2
 !---------------------------------------------------------------->>>>
         do 20 i=1,LNWAVS
         coef=(epsx(i)**2+epsy(i)**2)**2
         coef = coef/coefmax
#ifdef DIFFSW              
              ii=min(mod(i,iwav1),1)
              jj=((i-1)/iwav1)+1
              if ((ii.le.0.5*iwav1).and.(jj.le.0.5*jwav1)) then
                xkh(i)=0.0
                xkm(i)=0.0
#ifdef DIVDIF
                xkd(i)=0.0
#endif
              else
                xkh(i) = coef / difhsec
                xkm(i) = coef / difmsec
#ifdef DIVDIF
                xkd(i) = coef /(difmsec/TSDIVDIF)
#endif
              endif
#else
                xkh(i) = coef / difhsec
                xkm(i) = coef / difmsec
#ifdef DIVDIF
                xkd(i) = coef /(difmsec/TSDIVDIF)
#endif
#endif
         rkm(i) = 1.0 / ( 1.0 + 2.0 * ddt * xkm(i) )
         rkh(i) = 1.0 / ( 1.0 + 2.0 * ddt * xkh(i) )
#ifdef DIVDIF
         rkd(i) = 1.0 / ( 1.0 + 2.0 * ddt * xkd(i) )
#endif
 20      continue
 !----------------------------------------------------------------<<<<
!
#ifdef CHGTLEV
!  sigma level at which to begin rayleigh damping
        slrd0=0.012
!  reciprocal of time scale per scale height
!  above beginning sigma level for rayleigh damping
!       rtrd1=1./(5.*86400.)
        rtrd1=1./(1.*3600.)
!
        do k=1,levr
          if(sl(k).lt.slrd0) then
            rtrd(k)=rtrd1*log(slrd0/sl(k))
          else
            rtrd(k)=0
          endif
        enddo
#endif
!
!x       do 40 k=1,levr
!x       ku = min0 ( k+1, levr )
!x       kd = max0 ( k-1, 1     )
!x       sldsl = sl(k) / (sl(ku) - sl(kd)) / 4.0
!x       ct(k) = ( te(1,ku) - te(1,kd) ) * sldsl
!x       cq(k) = ( rq(1,ku) - rq(1,kd) ) * sldsl
!40      continue
!
      endif
#endif
!
! time splitting implicit diffusion
!
#ifdef ORIGIN_THREAD
!$doacross share(x,w,y,rkm,rkh,rkmt,rkht),local(i,k)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(x,w,y)
!mic$1 shared(rkm,rkh,rkmt,rkht)
!mic$1 private(i,k)
#endif
#ifdef OPENMP
!$omp parallel do private(i,k,rkmik)
#endif
        do 1100 k=1,levr
        do 1100 i=1,LNWAVS
#ifdef UPDIFF
!yj:top 3 levels (<18hPa) use 2nd diffusion, others 4th
            coefmaxt = (epsxmax**2+epsymax**2)**1
            coefmax = (epsxmax**2+epsymax**2)**NHDIFORD
            coefmaxh = (epsxmax**2+epsymax**2)**NHDIFORD
            if (k .ge. levr-LEVTOP) then
              difhsec2=1./15.
              difmsec2=difhsec2
              coef=(epsx(i)**2+epsy(i)**2)**1
              coef = coef/coefmaxt
              xkh(i) = coef * difhsec2
              xkm(i) = coef * difmsec2
#ifdef DIVDIF
              xkd(i) = xkm (i)
#endif
!            elseif (k.le.20) then
            elseif (k .le. 10) then
              coef=(epsx(i)**2+epsy(i)**2)**NHDIFORD
              coef = coef/coefmax
              xkh(i) = coef / difhsec
              xkm(i) = coef / difmsec
#ifdef DIVDIF
              xkd(i) = xkm (i)
#endif
            else
              coef=(epsx(i)**2+epsy(i)**2)**NHDIFORD
              coef = coef/coefmax
              coefh=(epsx(i)**2+epsy(i)**2)**NHDIFORD
              coefh = coefh/coefmaxh
#ifdef DIFFSW              
              ii=min(mod(i,iwav1),1)
              jj=((i-1)/iwav1)+1
              if ((ii.le.0.5*iwav1).and.(jj.le.0.5*jwav1)) then
                xkh(i)=0.0
                xkm(i)=0.0
#ifdef DIVDIF
                xkd(i)=0.0
#endif
              else
                xkh(i) = coefh / difhsec
                xkm(i) = coef / difmsec
#ifdef DIVDIF
                if (spdmax(k).ge.80.0) then
                xkd(i) = coef /(difmsec/TSDIVDIF)
                else
                xkd(i) = coef / difmsec
                endif
#endif
              endif
#else
                xkh(i) = coefh / difhsec
                xkm(i) = coef / difmsec
#ifdef DIVDIF
                if (spdmax(k).ge.80.0) then
                xkd(i) = coef /(difmsec/TSDIVDIF)
                else
                xkd(i) = coef / difmsec
                endif
#endif
#endif
!            else
!              mu=30                 !exp1=42; exp2=38
!              sig=3                 !exp1=9.; exp2=3
!              gfac=1.3              !exp1=1.5;exp2=1.3
!              sig22=2.*sig**2.
!              gcoef=1./sqrt(sig22*con_pi)
!              gexp=exp(-1.*((float(k-mu))**2./sig22))
!!              difhsec2=min(max(gcoef*gexp*gfac,1./difhsec),1./15.)
!              difhsec2=max(gcoef*gexp*gfac,1./difhsec)
!              difmsec2=difhsec2    !exp1 & exp2 =difhsec2*2./3.
!              coef=(epsx(i)**2+epsy(i)**2)**2
!              coef = coef/coefmax
!              xkh(i) = coef * difhsec2
!              xkm(i) = coef * difmsec2
            endif
            rkm(i) = 1.0 / ( 1.0 + 2.0 * ddt * xkm(i) )
            rkh(i) = 1.0 / ( 1.0 + 2.0 * ddt * xkh(i) )
#ifdef DIVDIF
            rkd(i) = 1.0 / ( 1.0 + 2.0 * ddt * xkd(i) )
#endif
#ifdef CHGTLEV
        rkmik = 1.0 / ( 1.0 + 2.0 * ddt *xkm(i) + ddt*rtrd(k) )
        x (i,k) = x (i,k) * rkmik
        w (i,k) = w (i,k) * rkmik
#else
        x (i,k) = x (i,k) * rkm(i)
!#ifdef DIFDZ        
!        w (i,k) = w (i,k) * rkh(i)
!#else        
#ifdef DIVDIF
        w (i,k) = w (i,k) * rkd(i)
#else
        w (i,k) = w (i,k) * rkm(i)
#endif
!#endif        
#endif
        y (i,k) = y (i,k) * rkh(i)
#ifdef NONHYD
        p (i,k) = p (i,k) * rkh(i)
        t (i,k) = t (i,k) * rkh(i)
        o (i,k) = o (i,k) * rkm(i)
#endif

#else
 !else  UPDIFF

#ifdef CHGTLEV
        rkmik = 1.0 / ( 1.0 + 2.0 * ddt *xkm(i) + ddt*rtrd(k) )
        x (i,k) = x (i,k) * rkmik
        w (i,k) = w (i,k) * rkmik
#else
        x (i,k) = x (i,k) * rkm(i)
!#ifdef DIFDZ        
!        w (i,k) = w (i,k) * rkh(i)
!#else        
#ifdef DIVDIF
        w (i,k) = w (i,k) * rkd(i)
#else
        w (i,k) = w (i,k) * rkm(i)
#endif
!#endif        
#endif
        y (i,k) = y (i,k) * rkh(i)
#ifdef NONHYD
        p (i,k) = p (i,k) * rkh(i)
        t (i,k) = t (i,k) * rkh(i)
        o (i,k) = o (i,k) * rkm(i)
#endif

#endif
 !endif UPDIFF
1100    continue



#ifdef NONHYD
        k=levr+1
        do i=1,LNWAVS
        o (i,k) = o (i,k) * rkm(i)
        enddo
#endif
!
#ifndef NDSL
#ifdef ORIGIN_THREAD
!$doacross share(rt,rkh),local(i,k)
#endif
#ifdef CRAY_THREAD
!mic$ do all
!mic$1 shared(rt)
!mic$1 shared(rkh)
!mic$1 private(i,k)
#endif
#ifdef OPENMP
!$omp parallel do private(nc,i,k)
#endif
        do 1200 nc=1,ntotal
        do 1200 k=1,levr
        do 1200 i=1,LNWAVS
#ifdef UPDIFF
!yj:top 3 levels (<18hPa) use 2nd diffusion, others 4th
            coefmaxt = (epsxmax**2+epsymax**2)**1
            coefmax = (epsxmax**2+epsymax**2)**NHDIFORD
            if (k .ge. levr-LEVTOP) then
              difhsec2=1./15.
              difmsec2=difhsec2
              coef=(epsx(i)**2+epsy(i)**2)**1
              coef = coef/coefmaxt
              xkh(i) = coef * difhsec2
!            elseif (k.le.20) then
            else
              coef=(epsx(i)**2+epsy(i)**2)**NHDIFORD
              coef = coef/coefmax
              xkh(i) = coef / difhsec
            endif
            rkh(i) = 1.0 / ( 1.0 + 2.0 * ddt * xkh(i) )
#endif
        rt(i,k,nc) = rt(i,k,nc) * rkh(i)
1200    continue
#endif
!
#ifdef DIFDZ
!yj!-------------------change dz forcing to uv---------
      call dz2uv(LNWAVS,x,w,tmp,tmp1,epsx,epsy,rnnp1,levr)
!!mic$ do all
!!mic$1 shared(tmp,tmp1,x,w)
!!mic$1 private(i)
      do 870 k=1,levr
      do 870 i=1,LNWAVS
      x(i,k) = tmp (i,k)
      w(i,k) = tmp1(i,k)
 870  continue
!yj!-------------------change dz forcing to uv---------
#endif
      return
      end
