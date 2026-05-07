      subroutine gfstopo(nx,my,jtrun,mlmax,nn,xt,topog)
c     subroutine gfstopo(nx,my,jtrun,mlmax,nn,xt,topog,xtm)
c
c  based on tom's code from c. s. liou       
c  purpose : to get necessary topo datas ,including silhouette terrain
c            , topo spectral coefs. and trunated-filtered terrain,
c            from  5" topo dataset.
c
c!!!!!!!!!!!!!!!!!!!!
c while change dimension nx,my ,something else must be noticed:
c you must change 1.   im,jm           in subroutine trans, transr
c                 2.   common /fft/
c!!!!!!!!!!!!!!!!!!!!
c
c      implicit none
c      include './bck_dim.h'
c     parameter (nx=720, my=nx/2)
      integer   nx,my,nxg,myg,err,jm2,i,j,nn,istat,idms,nz
      integer   jtrun,mlmax
      real      pi,hmax,hmin,gs,ahh
      character*1   grd
      parameter (nxg=360*120, myg=nxg/2)
c     character ifile*34,flagn*4,flag*4,flago*4
c     character*72 type_w,argument
      dimension  xt(nx,my),xtm(nx,my),sinl(my),weight(my),
     &           topog(nx,my),zdvg(nx,my),xtg(nx,my),
     &           poly(mlmax,my/2),dpoly(mlmax,my/2),
     &           mlsort(jtrun,jtrun),msort(mlmax),lsort(mlmax),
     &           filt(mlmax),topos(mlmax,2),xts(mlmax,2),
     &           zdvs(mlmax,2)

c------------------------------------
c       call dmsmsg("all",istat)
c
c       type_w='WORDER'//char(0)
c       argument='34'//char(0)
c       call dmscfg(type_w,argument,istat_w)
c
c     ifile='BCK'
c
c     call dmsopn(ifile,"w",istat)
c     if(istat.eq.0)print*,'open dmsfile ok, file=',ifile
      
            
      print*,'jtrun ,  mlmax = ',jtrun,mlmax

c
c  build pointer arrays for locating zonal and total wavenumber
c  values in the one-dimensional spherical harmonic arrays.
c
      call sortml(jtrun,mlmax,msort,lsort,mlsort)      
c
      print*,'1 step ok'
c
c  gaussian quadrature weights and latitudes
c
      call gausl3(my,-1.0,1.0,weight,sinl)
c
c
      print*,'2 step ok'
c  associated legendre polynomials and their derivatives
c
      jm2=my/2
      call lgndr(jm2,jtrun,mlmax,mlsort,sinl,poly,dpoly)
c
      print*,'3 step ok'
c
c  difine how much area a grid point characters
c (delete)
      pi=4.0*atan(1.0)
c
c  find max and min height of topo
c
      hmax=0.
      hmin=0.
      do 250 j=1,my
      do 250 i=1,nx
      hmax=max(hmax,xt(i,j))
      hmin=min(hmin,xt(i,j))
 250  continue

      print*,'max height:',hmax
      print*,'min height:',hmin
c
c  put datd into dmsfile
c
c      do i=1,nx*my
c      tabls(i,1)=float(ils(i,1))
c      end do
c
clzl      idms=1
      idms=1
c
c      if(idms .eq.1)then
c      call crtdms34(ifile,nx*my,xt ,'s00060','gbkf',flag)
c      call crtdms34(ifile,nx*my,xtm,'s00062','gbkf',flag)
c      call crtdms34(ifile,nx*my,tabls,'s00070','gbck',flag)
c     call crtdms(ifile,nx*my,xt,'s06','gg',99)
c     call crtdms(ifile,nx*my,xtm,'s08','gg',99)
c     call crtdms(ifile,nx*my,tabls,'s07','gg',00)
c     call crtdms(ifile,nx*my,facls,'s09','gg',00)
c      end if
c
c      print*,'4 step ok'
c just do land-sea table
ccc   go to 611
c
c  transfer meter to geopotential (m**2/s**2)
c
      gs=9.80616
      do j=1,my
        do i=1,nx
          xtg(i,j)=xt(i,j)*gs
cx          xtm(i,1)=xtm(i,1)*gs
        end do
      end do
c
c     grid to spectral
c
      call tranrs(jtrun,mlmax,nx,my,poly,weight,xtg,xts)
c      call flush(6)
cx      call tranrs(jtrun,mlmax,nx,my,poly,weight,xtm,zdvs)
      print*,'4 step ok'
c
c  spectral to grid for topo standard deviation
c
cx      call transr(jtrun,mlmax,nx,my,poly,zdvs,zdvg)
cx      do i=1,nx*my
cx        zdvg(i,1)=zdvg(i,1)/gs
cx      end do
cx      if(idms .eq. 1)call crtdms(ifile,nx*my,zdvg,'s08','gg',0)
c
c      do 610 nn=1,2
c
c  lanczos filter for terrain smoothing
c
        ahh=nn
        nz=0
       print*,'before filter ahh =',ahh
       call filter(ahh,lsort,mlmax,jtrun,filt,xts,topos)
c
c      call flush(6)
c   spectral to grid
c
      call transr(jtrun,mlmax,nx,my,poly,topos,topog)
      print*,'5-',nn,' step ok'
c      call flush(6)
c
c  transfer geopotential to meter
c
      do j=1,my
        do i=1,nx
          topog(i,j)=topog(i,j)/gs
        end do
      end do
c
c     if(idms .eq. 1) then
c       flagn='gbkf'
c       flago='gbkf'
c       write(flagn(4:4),'(i1.1)')nn
c       call getenv("flag",flag)
c       write(flago(4:4),'(i1.1)')nz
c       print*,'flagn=',flagn,'flago=',flago,' flag=',flag
c
c       call crtdms34(ifile,nx*my,topog,'s00060',flagn,flag)
c       call crtdms34(ifile,nx*my,xtm,'s00062',flago,flag)
c      call crtdms34(ifile,mlmax*2,topos,'s00060',flagn,flags)
c     end if

c     if(idms .eq. 1) then
c     call crtdms(ifile,nx*my,topog,'s06','gg',nn)
c     call crtdms(ifile,mlmax*2,topos,'s06','gs',nn)
c     end if
c     print*,'6-3 step ok'
c     call flush(6)
c
cx...............................
c
clzl      iav=2*((1+4320/nx)/2)
cOKlzl      iav=2*((1+nxg/nx)/2)
c
c      do 2 jq=1,my
c      alat(jq)=asin(sinl(jq))
c     rjj=1.5+2160.*(alat(jq)/pi+0.5)
clzl      rjj=1.0+2160.*(alat(jq)/pi+0.5)
c      rjj=1.0+myg*(alat(jq)/pi+0.5)
c      jj=rjj-0.5*iav
cfong iavi=iav/cos(alat(jq))
c      iavi=iav
c      iavi=2*(iavi/2)
cfong iavi=min(50,iavi)
c
c      do 42 j=1,iav
c      jp=jj+j-1
clzl      do 142 i=1,4320
c      do 142 i=1,nxg
c      ibits(i,j)=ils1(i,jp)
c      x(i,j)=itopo(i,jp)
c      if(ibits(i,j) .eq. 0) x(i,j)=0.0
c  142 continue
c   42 continue

c
c  get characteristic topo value, topo standard deviation, adn land-
c  sea bit within the area
c  the characteristic topo value gets from silhouette method
c
c      call ciloet2(jj,iav,iavi,nx,x,topog(1,jq),xtm(1,jq))

c    2 continue
c
c      do 450 j=1,my
c      do 450 i=1,nx
c      if(.not.bils(i,j))xtm(i,j)=0.0
c  450 continue
c
c      if(idms .eq. 1) then
c        call crtdms34(ifile,nx*my,xtm,'s00062',flagn,flag)
c     call crtdms(ifile,nx*my,xtm,'s08','gg',nn-1)
c      end if
c
c      print*,'6-4 step ok'
cc      call flush(6)
cx...............................
c
c
c  find max and min height of topo
c
      hmax=0.
      hmin=0.
      do 500 j=1,my
      do 500 i=1,nx
      hmax=max(hmax,topog(i,j))
      hmin=min(hmin,topog(i,j))
 500  continue

      print*,'max height:',hmax
      print*,'min height:',hmin
c 610  continue
c 611  continue
c
c     call dmscls(ifile,istat)
c      call dmsexit(0)
c
c
      return
      end
c
      subroutine filter (azero,lsort,mlmax,jtrun,filt,sgeo,sgeox)
      dimension filt(mlmax),sgeo(mlmax,2),lsort(mlmax),sgeox(mlmax,2)
c
      filt=0.0
      pi= 4.0*atan(1.0)
      jj= jtrun-1
c
      print*,' in filter azero = ',azero
      print*,' in filter jtrun = ',jtrun
      print*,' in filter jj = ',jj
c
      rj= 1.0/(jj)
c
      do 1 l=2,jj
c
      fac= pi*(l-1.0)*rj
      filt(l)= sin(fac)/fac
clzl+add--
      if( azero .eq. 0.0 )then 
        filt(l)= filt(l)
      else
      filt(l)= filt(l)**azero
      endif
c
    1 continue
      filt(1)= 1.0
      do 2 ml=1,mlmax
      l= lsort(ml)
      sgeox(ml,1)= sgeo(ml,1)*filt(l)
      sgeox(ml,2)= sgeo(ml,2)*filt(l)
    2 continue
c
      return
      end

      subroutine gausl3 (n,xa,xb,wt,ab)
c
c weights and abscissas for nth order gaussian quadrature on (xa,xb).
c input arguments
c
c n  -the order desired
c xa -the left endpoint of the interval of integration
c xb -the right endpoint of the interval of integration
c output arguments
c ab -the n calculated abscissas
c wt -the n calculated weights
c
      implicit double precision (a-h,o-z)
c
      real  ab(n) ,wt(n),xa,xb
c
c machine dependent constants---
c  tol - convergence criterion for double precision iteration
c  pi  - given to 15 significant digits
c  c1  -  1/8                     these are coefficients in mcmahon"s
c  c2  -  -31/(2*3*8**2)          expansions of the kth zero of the
c  c3  -  3779/(2*3*5*8**3)       bessel function j0(x) (cf. abramowitz,
c  c4  -  -6277237/(3*5*7*8**5)   handbook of mathematical functions).
c  u   -  (1-(2/pi)**2)/4
c
      data tol/1.d-14/,pi/3.14159265358979/,u/.148678816357662/
      data c1,c2,c3,c4/.125,-.080729166666667,.246028645833333,
     1                -1.82443876720609 /
c
c maximum number of iterations before giving up on convergence
c
      data maxit /5/
c
c arithmetic statement function for converting integer to double
c
      dbli(i) = dble(float(i))
c
      ddif = .5d0*(dble(xb)-dble(xa))
      dsum = .5d0*(dble(xb)+dble(xa))
      if (n .gt. 1) go to 101
      ab(1) = 0.
      wt(1) = 2.*ddif
      go to 107
  101 continue
      nnp1 = n*(n+1)
      cond = 1./sqrt((.5+float(n))**2+u)
      lim = n/2
c
      do 105 k=1,lim
	 b = (float(k)-.25)*pi
	 bisq = 1./(b*b)
c
c rootbf approximates the kth zero of the bessel function j0(x)
c
	 rootbf = b*(1.+bisq*(c1+bisq*(c2+bisq*(c3+bisq*c4))))
c
c      initial guess for kth root of legendre poly p-sub-n(x)
c
	 dzero = cos(rootbf*cond)
	 do 103 i=1,maxit
c
	    dpm2 = 1.d0
	    dpm1 = dzero
c
c       recursion relation for legendre polynomials
c
	    do 102 nn=2,n
		dp = (dbli(2*nn-1)*dzero*dpm1-dbli(nn-1)*dpm2)/dbli(nn)
		dpm2 = dpm1
		dpm1 = dp
  102       continue
	    dtmp = 1.d0/(1.d0-dzero*dzero)
	    dppr = dbli(n)*(dpm2-dzero*dp)*dtmp
	    dp2pri = (2.d0*dzero*dppr-dbli(nnp1)*dp)*dtmp
	    drat = dp/dppr
c
c       cubically-convergent iterative improvement of root
c
	    dzeri = dzero-drat*(1.d0+drat*dp2pri/(2.d0*dppr))
	    ddum= dabs(dzeri-dzero)
	 if (ddum .le. tol) go to 104
	    dzero = dzeri
  103    continue
	 print 504
  504    format(1x,' in gausl3, convergence failed')
  104    continue
	 ddifx = ddif*dzero
	 ab(k) = dsum-ddifx
	 wt(k) = 2.d0*(1.d0-dzero*dzero)/(dbli(n)*dpm2)**2*ddif
	 i = n-k+1
	 ab(i) = dsum+ddifx
	 wt(i) = wt(k)
  105 continue
c
      if (mod(n,2) .eq. 0) go to 107
      ab(lim+1) = dsum
      nm1 = n-1
      dprod = n
      do 106 k=1,nm1,2
	 dprod = dbli(nm1-k)*dprod/dbli(n-k)
  106 continue
      wt(lim+1) = 2.d0/dprod**2*ddif
  107 return
      end
      subroutine sortml(jtrun,mlmax,msort,lsort,mlsort)
c
c  sortml builds pointer arrays for functional dependency between
c  1-d spectral index and zonal and total wavenumber indices.  this
c  subroutine reflects the coefficient storage strategy used in the
c  model
c
c ***input***
c
c  jtrun:  zonal and total wavenumber limit
c  mlmax:  total number of 1-d spectral index for triangular trunc
c
c  ***output***
c
c  msort:  zonal wavenumber as function of 1-d spectral index
c  lsort:  total wavenumber as function of 1-d spectral index
c  mlsort: total wavenumber index as function of zonal and total
c          wavenumber
c
c *******************************************************************
c
      dimension msort(mlmax),lsort(mlmax),mlsort(jtrun,jtrun)
c
      mlx= (jtrun/2)*((jtrun+1)/2)
      ml= 0
      do 1 k=1,jtrun-1,2
      do 1 m=1,jtrun-k
      ml= ml+1
      mlp= ml+mlx
      mlsort(m,m+k)= mlp
      mlsort(m,m+k-1)= ml
      msort(ml)= m
      lsort(ml)= m+k-1
      msort(mlp)= m
      lsort(mlp)= m+k
    1 continue
c
      ml= mlp
      do 2 m=2,jtrun,2
      ml= ml+1
      mlsort(m,jtrun)= ml
      msort(ml)= m
      lsort(ml)= jtrun
    2 continue
      return
      end
c
      subroutine lgndr(jm2,jtrun,mlmax,mlsort,sinl,poly,dpoly)
c
c  generate legendre polynomials and their derivatives on the
c  gaussian latitudes
c
c ***input***
c
c  jm2:  number of gaussian latitudes from south pole and equator
c  jtrun:  zonal wavenumber truncation limit
c  mlmax: total number of triangular truncation spherical harmonics
c  mlsort: pointer array of 1-d indexs at functions of zonal and
c          total wavenumbers
c  sinl: sin of gaussian latitudes
c
c  ***output***
c
c  poly: associated legendre coefficients
c  dpoly: d(poly)/d(sinl)
c
c ******************************************************************
c
c ref= belousov, s. l., 1962= tables of normalized associated
c        legendre polynomials. pergamon press, new york
c
      dimension poly(mlmax,jm2),dpoly(mlmax,jm2),sinl(jm2)
     *, mlsort(jtrun,jtrun)
c
c CWB orig      dimension pnm(jtrun+1,jtrun),dpnm(jtrun+1,jtrun)
      dimension pnm(jtrun+1,jtrun+1),dpnm(jtrun+1,jtrun+1)
c
c sinl is sin(latitude) = cos(colatitude)
c pnm(np,mp) is legendre polynomial p(n,m) with np=n+1, mp=m+1
c pnm(mp,np+1) is x derivative of p(n,m) with np=n+1, mp=m+1
c
      jtrunp= jtrun+1
      do 1001 j=1,jm2
      xx= sinl(j)
      sn= sqrt(1.0-xx*xx)
       sn2i = 1.0/(1.0 - xx*xx)
      rt2= sqrt(2.0)
       c1 = rt2
c
       pnm(1,1) = 1.0/rt2
      theta=-atan(xx/sqrt(1.0-xx*xx))+2.0*atan(1.0)
c
      do 20 n=1,jtrun
       np = n + 1
      fn=n
       fn2 = fn + fn
       fn2s = fn2*fn2
c eq 22
      c1= c1*sqrt(1.0-1.0/fn2s)
      c3= c1/sqrt(fn*(fn+1.0))
       ang = fn*theta
       s1 = 0.0
       s2 = 0.0
       c4 = 1.0
       c5 = fn
       a = -1.0
       b = 0.0
c
      do 27 kp=1,np,2
       k = kp - 1
      s2= s2+c5*sin(ang)*c4
      if (k.eq.n) c4 = 0.5*c4
      s1= s1+c4*cos(ang)
       a = a + 2.0
       b = b + 1.0
      fk=k
       ang = theta*(fn - fk - 2.0)
       c4 = (a*(fn - b + 1.0)/(b*(fn2 - a)))*c4
       c5 = c5 - 2.0
   27 continue
c eq 19
       pnm(np,1) = s1*c1
c eq 21
       pnm(np,2) = s2*c3
   20 continue
c
      do 4 mp=3,jtrunp
       m = mp - 1
      fm= m
       fm1 = fm - 1.0
       fm2 = fm - 2.0
       fm3 = fm - 3.0
      c6= sqrt(1.0+1.0/(fm+fm))
c eq 23
       pnm(mp,mp) = c6*sn*pnm(m,m)
      if (mp - jtrunp) 3,4,4
    3 continue
       nps = mp + 1
c
      do 41 np=nps,jtrunp
       n = np - 1
      fn= n
       fn2 = fn + fn
       c7 = (fn2 + 1.0)/(fn2 - 1.0)
       c8 = (fm1 + fn)/((fm + fn)*(fm2 + fn))
      c= sqrt((fn2+1.0)*c8*(fm3+fn)/(fn2-3.0))
      d= -sqrt(c7*c8*(fn-fm1))
      e= sqrt(c7*(fn-fm)/(fn+fm))
c eq 17
       pnm(np,mp) = c*pnm(np-2,mp-2)
     1            + xx*(d*pnm(np-1,mp-2) + e*pnm(np - 1,mp))
   41 continue
    4 continue
c
      do 50 mp=1,jtrun
      fm= mp-1.0
       fms = fm*fm
      do 50 np=mp,jtrun
      fnp= np
       fnp2 = fnp + fnp
       cf = (fnp*fnp - fms)*(fnp2 - 1.0)/(fnp2 + 1.0)
      cf= sqrt(cf)
c der
      dpnm(np,mp)   = -sn2i*(cf*pnm(np+1,mp) - fnp*xx*pnm(np,mp))
   50 continue
c
      do 71 m=1,jtrun
      do 71 l=m,jtrun
      ml= mlsort(m,l)
      poly(ml,j)= pnm(l,m)
      dpoly(ml,j)=dpnm(l,m)
   71 continue
      dpoly(1,j)= 0.0
 1001 continue
      return
      end
c
      subroutine tranrs(jtrun,mlmax,nx,my,poly,w,r,s)
c
c  subroutine to transform a scalar grid point field to spectral
c  coefficient
c
c *** input ***
c
c  jtrun: zonal wavenumber truncation limit
c  mlmax: total number of spherical harmonic coeff  (horizontal field)
c  nx: e-w dimension no.
c  my: n-s dimension no.
c  ll: number of vertical levels to transform
c  poly: legendre polynomials
c  w: gaussian quadrature weights
c  r: 3-dim input grid pt. field to be transformed
c
c *** output ***
c
c  s: spectral coefficient fields
c
c  **********************************
c
      dimension poly(mlmax,my/2),s(mlmax,2),r(nx,my),w(my)
c     include 'fftcom.h'             
ct180 common/fft/ trigs(1024),ifax(19)
c     common/fft/ trigs(2048),ifax(19)
c      common/fft/ trigs(12288),ifax(19)
      dimension trigs(4096),ifax(19)
c-lzl-orig      dimension cc(nx+3,my),work(nx*my,2)
      dimension cc(nx+2,my),work(nx*my,2)
      mlx= (jtrun/2)*((jtrun+1)/2)
c
c     do 30 k=1,ll
c
c  put grid point fields into two dimensional horizontal array
c
      do 23 j=1,my
      do 23 i=1,nx
      cc(i,j)= r(i,j)
   23 continue
c
c  fft for each guassian latitude of 2-d field
c
      call fftfax(nx,ifax,trigs)
      call rfftmlt(cc,work,trigs,ifax,1,nx+2,nx,my,-1)
c-lzl-orig      call rfftmlt(cc,work,trigs,ifax,1,nx+3,nx,my,-1)
c
c  to start quadrature integral we compute contribution without
c  adding to existing sum
c
      m1= 0
      do 62 l=jtrun-1,1,-2
cdir@    ivdep
	 do 63 m = 1, l
	 mm= 2*m-1
	 mp= mm+1
	 ml= m+m1
	 mk= ml+mlx
	    s(ml,1) = w(1)*poly(ml,1)*(cc(mm,1)+cc(mm,my))
	    s(ml,2) = w(1)*poly(ml,1)*(cc(mp,1)+cc(mp,my))
	    s(mk,1) = w(1)*poly(mk,1)*(cc(mm,1)-cc(mm,my))
	    s(mk,2) = w(1)*poly(mk,1)*(cc(mp,1)-cc(mp,my))
   63    continue
      m1= m1+l
   62 continue
c
      ml= mlx*2
cdir@ ivdep
      do 64 m=2,jtrun,2
      ml=ml+1
      mm= 2*m-1
      mp= mm+1
      s(ml,1)= w(1)*poly(ml,1)*(cc(mm,1)+cc(mm,my))
      s(ml,2)= w(1)*poly(ml,1)*(cc(mp,1)+cc(mp,my))
   64 continue
c
c  for rest of quadrature integral we add contribution to
c  existing sum
c
      do 70 j=2,my/2
      jj= my-j+1
      m1= 0
      do 72 l=jtrun-1,1,-2
cdir@    ivdep
	 do 73 m = 1, l
	 mm= 2*m-1
	 mp= mm+1
	 ml= m+m1
	 mk= ml+mlx
	    s(ml,1) = s(ml,1)+w(j)*poly(ml,j)*(cc(mm,j)+cc(mm,jj))
	    s(ml,2) = s(ml,2)+w(j)*poly(ml,j)*(cc(mp,j)+cc(mp,jj))
	    s(mk,1) = s(mk,1)+w(j)*poly(mk,j)*(cc(mm,j)-cc(mm,jj))
	    s(mk,2) = s(mk,2)+w(j)*poly(mk,j)*(cc(mp,j)-cc(mp,jj))
   73    continue
      m1= m1+l
   72 continue
c
      ml= mlx*2
cdir@ ivdep
      do 65 m=2,jtrun,2
      ml=ml+1
      mm= 2*m-1
      mp= mm+1
      s(ml,1)= s(ml,1)+w(j)*poly(ml,j)*(cc(mm,j)+cc(mm,jj))
      s(ml,2)= s(ml,2)+w(j)*poly(ml,j)*(cc(mp,j)+cc(mp,jj))
   65 continue
   70 continue
c  30 continue
c
      return
      end
c
      subroutine transr(jtrun,mlmax,nx,my,poly,s,r)
c
c  subroutine to transform a spectral coefficient field to
c  grid point form
c
c *** input ***
c
c  jtrun: zonal wavenumber resolution limit
c  mlmax: number of spectral coefficients (horizontal field)
c  nx: e-w dimension no.
c  my: n-s dimension no.
c  ll: number of levels to transform
c  poly: legendre polynomials
c  s: spectral coefficient array to transform
c
c *** output ***
c
c  r: 3-d output grid point fields
c
c  **************************************
c
      dimension poly(mlmax,my/2),s(mlmax,2),r(nx,my)
c     include 'fftcom.h'              
ct180 common/fft/ trigs(1024),ifax(19)
c     common/fft/ trigs(2048),ifax(19)
c      common/fft/ trigs(12288),ifax(19)
      dimension trigs(4096),ifax(19)
c-lzl-orig      dimension cc(nx+3,my),work(nx*my,2)
      dimension cc(nx+2,my),work(nx*my,2)
c
      mlx= (jtrun/2)*((jtrun+1)/2)
c     do 20 k=1,ll
clzl      do 55 m=1,(nx+3)*my/2
      do 55 m=1,(nx+2)*my/2
      cc(m,1)= 0.0
      cc(m,my/2+1)= 0.0
   55 continue
c
      do 5 j=1,my/2
      jj= my+1-j
      ml= 2*mlx
cdir@ ivdep
      do 3 m=2,jtrun,2
      ml= ml+1
      mm= 2*m-1
      mp= mm+1
      cc(mm,j)= poly(ml,j)*s(ml,1)
      cc(mp,j)= poly(ml,j)*s(ml,2)
      cc(mm,jj)= cc(mm,j)
      cc(mp,jj)= cc(mp,j)
    3 continue
c
      m1= 0
      do 5 l=jtrun-1,1,-2
cdir@ ivdep
      do 6 m=1,l
      mm= 2*m-1
      mp= mm+1
      ml= m+m1
      mk= ml+mlx
      cc(mm,j)= cc(mm,j)+poly(ml,j)*s(ml,1)+poly(mk,j)*s(mk,1)
      cc(mm,jj)=cc(mm,jj)+poly(ml,j)*s(ml,1)-poly(mk,j)*s(mk,1)
      cc(mp,j)= cc(mp,j)+poly(ml,j)*s(ml,2)+poly(mk,j)*s(mk,2)
      cc(mp,jj)=cc(mp,jj)+poly(ml,j)*s(ml,2)-poly(mk,j)*s(mk,2)
    6 continue
      m1= m1+l
    5 continue
c
      call fftfax(nx,ifax,trigs)
      call rfftmlt(cc,work,trigs,ifax,1,nx+2,nx,my,1)
c-lzl-orig      call rfftmlt(cc,work,trigs,ifax,1,nx+3,nx,my,1)
c
      do 22 j=1,my
      do 22 i=1,nx
      r(i,j)= cc(i,j)
   22 continue
c
   20 continue
c
      return
      end

