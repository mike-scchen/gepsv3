!----------------------------------------------------------------------------
      subroutine rndslfv_init
!yj      subroutine ndslfv_init(lonf,latg,ntrac,
!yj     &                       coslat,colrad,wgt,lats_nodes_a)
!
! initialize some constant variable/arrays for ndsl advection -hmhj
!
! program log
! 2011 02 20 : henry jaung, initiated for ndsl advection
! 2013 06 20 : henry jaung, enhance cyclic condition for extra strong wind in wam
! 2013 09 30 : henry jaung, add option to advect between regular and
!                           reduced grids.
! 2018 07 02 : YING-JU CHEN, modified for RSM
!
!yj      use gfs_dyn_physcons
!yj      use gfs_dyn_layout1
         use mod_rndsl
!
      implicit none
#include <paramodel.h>
#include <npesi.h>      
#include <comrmpi.h>

!yj      real   ,intent(in):: coslat(latg)
!yj      real   ,intent(in):: colrad(latg/2)
!yj      real   ,intent(in):: wgt   (latg/2)
!yj      integer,intent(in):: lats_nodes_a(nodes), lonf, latg, ntrac
!yj      integer,intent(in):: ntrac = ntotal in paramodel.h
!
!yj      integer	jm2,jm,jmh,i,j
!yj      real 	pi,hfpi,twopi,dlat,dlon
!yj      real 	y0,y1,y2,y3,y4,dyregg

      logical   lprint
!
      integer	n,nm,nr,lonp,latp
!
      lprint = .false.
!yj      pi = con_pi
!yj      hfpi  = pi * 0.5
!yj      twopi = pi * 2.0
!yj      jm  = latg
!yj      jm2 = jm*2
!yj      jmh = jm/2
! ------------------------------------------------------------
! ------------------- gaussian latitude 
! cosglat   cos(lat)
! gglat     gaussian latitude count from 0 to 2 pi
!           from north pole to south pole back to north pole
! gglati    gaussian latitude at edges for gglat
! gslati    cos(lat) at gglati
! ggfact    interpolation factor from gglat to gglati for wind
!
!yj      allocate ( cosglat(jm) )
!yj      allocate ( gglat(jm2), gglati(jm2+1) )
!yj      allocate ( gslati(jm2+1), ggfact(jm2,4) )
!
!yj      do j=1,jm
!yj        cosglat(j)=coslat(j)
!yj      enddo
!
!yj      do j=1,jmh
!yj        gglat(      j) =         colrad(j)
!yj        gglat(jm +1-j) =    pi - colrad(j)
!yj        gglat(jm +  j) =    pi + colrad(j)
!yj        gglat(jm2+1-j) = twopi - colrad(j)
!yj      enddo
!
!yj      gslati(    1)=0.0
!yj      gslati(jm +1)=2.0
!yj      gslati(jm2+1)=4.0
!yj      do j=1,jmh
!yj        gslati(    j+1) = gslati(    j  ) + wgt(j)
!yj        gslati(jm -j+1) = gslati(jm -j+2) - wgt(j)
!yj        gslati(jm +j+1) = gslati(jm +j  ) + wgt(j)
!yj        gslati(jm2-j+1) = gslati(jm2-j+2) - wgt(j)
!yj      enddo

!yj      if( me == 0 .and. lprint ) then
!yj      print *,' cell_edge=1    gslati=',gslati(1)
!yj      do j=1,jm2
!yj      print *,'    cell_mean=',j,'  dgslat= ',gslati(j+1)-gslati(j)
!yj      print *,' cell_edge=',j+1,' gslati= ',gslati(j+1)
!yj      enddo
!yj      endif
!
! real latitude values first at edge, use temporary ggfact
!yj      gglati(1) = hfpi
!     print *,'        gglati =',gglati(1)
!yj      do j=1,jmh
!yj        gglati(j+1) = asin( sin(gglati(j)) - wgt(j) )
!yj        ggfact(j,1) = gglati(j) - gglati(j+1)
!       print *,' j dgglat ',j,ggfact(j,1)
!       print *,'        gglati =',gglati(j+1)
!yj      enddo
!
! then relative latitude from 0 to 2pi
!yj      gglati(    1)=0.0
!yj      gglati(jm +1)=pi
!yj      gglati(jm2+1)=twopi
!yj      do j=1,jmh
!yj        gglati(    j+1) = gglati(    j  ) + ggfact(j,1)
!yj        gglati(jm -j+1) = gglati(jm -j+2) - ggfact(j,1)
!yj        gglati(jm +j+1) = gglati(jm +j  ) + ggfact(j,1)
!yj        gglati(jm2-j+1) = gglati(jm2-j+2) - ggfact(j,1)
!yj      enddo
!
! -------------------------
!yj      if( me == 0 .and. lprint ) then
!yj      print *,' cell_edge=1    gglati=',gglati(1)
!yj      do j=1,jm2
!yj      print *,'  dg1 =',gglat(j)-gglati(j)
!yj      print *,'                   cell_mean=',j,'  gglat= ',gglat(j)
!yj      print *,'  dg2 =',gglati(j+1)-gglat(j)
!yj      print *,' cell_edge=',j+1,' gglati= ',gglati(j+1)
!yj      enddo
!yj      endif
!
! ------------- interpolation factor ---------------
!yj      do j=1,jm2
!yj        y0 = gglati(j)
!yj        if( j == 1 ) then
!yj          y1 = -gglat(2)
!yj          y2 = -gglat(1)
!yj          y3 =  gglat(1)
!yj          y4 =  gglat(2)
!yj        else if( j.eq.2 ) then
!yj          y1 = -gglat(1)
!yj          y2 =  gglat(1)
!yj          y3 =  gglat(2)
!yj          y4 =  gglat(3)
!yj        else if( j == jm2 ) then
!yj          y1 = gglat(jm2-2)
!yj          y2 = gglat(jm2-1)
!yj          y3 = gglat(jm2)
!yj          y4 = twopi+gglat(1)
!yj        else
!yj          y1 = gglat(j-2)
!yj          y2 = gglat(j-1)
!yj          y3 = gglat(j)
!yj          y4 = gglat(j+1)
!yj        endif
!yj        ggfact(j,4) = (y0-y1)*(y0-y2)*(y0-y3)/((y4-y1)*(y4-y2)*(y4-y3))
!yj        ggfact(j,3) = (y0-y1)*(y0-y2)*(y0-y4)/((y3-y1)*(y3-y2)*(y3-y4))
!yj        ggfact(j,2) = (y0-y1)*(y0-y3)*(y0-y4)/((y2-y1)*(y2-y3)*(y2-y4))
!yj        ggfact(j,1) = (y0-y2)*(y0-y3)*(y0-y4)/((y1-y2)*(y1-y3)*(y1-y4))
!yj      enddo
!
! -------------------------------------------------------------
!  determin the longitude with full grid
!
!yj      allocate( gglon(lonf), ggloni(lonf+1) )
!yj      dlon = pi * 2. / lonf
!yj      do i=1,lonf
!yj        gglon(i)=(i-1)*dlon
!yj      enddo
!yj      do i=2,lonf
!yj        ggloni(i)=0.5*(gglon(i-1)+gglon(i))
!yj      enddo
!yj      ggloni(     1)=ggloni(   2)-dlon
!yj      ggloni(lonf+1)=ggloni(lonf)+dlon
!
!yj      if( me == 0 .and. lprint ) then
!yj      print *,'  dlon=',dlon
!yj      print *,' edge number=1    ggloni=',ggloni(1)
!yj      do i=1,lonf
!yj      print *,'               cell number=',i,'  gglon= ',gglon(i)
!yj      print *,' edge number=',i+1,' ggloni= ',ggloni(i+1)
!yj      enddo
!yj      endif
!
! --------------------- for parallel --------------------
!
      allocate( lonstrndsl(nrow), lonlenndsl(nrow) )
      allocate( latstrndsl(nrow), latlenndsl(nrow) )
!yj      allocate( lonstrndsl(nodes), lonlenndsl(nodes) )
!yj      allocate( latstrndsl(nodes), latlenndsl(nodes) )
!     print*,'enter rndsl_init,mype=',mype
!
      lonfull = igrd1  !!yj longfull=igrd1
!yj      lonfull = lonf  !!yj longfull=igrd1
!yj      lonhalf = lonf / 2	! lonf has to be even !!yj no need
      lonpart = (igrd1-1)/nrow+1  !!yj lonpart=(igrd1-1)/nrow+1
!yj      lonpart = (lonhalf-1)/nodes+1  !!yj lonpart=(igrd1-1)/nrow+1
!
      lonp = igrd1 / nrow  !!yj lonp=igrd1/nrow
!yj      lonp = lonhalf / nodes  !!yj lonp=igrd1/nrow
!yj      do n=1,nodes
      do n=1,nrow
        lonlenndsl(n)=lonp
      enddo
      nr=mod(igrd1,nrow)  !!yj nr=mod(igrd1/nrow)
!yj      nr=mod(lonhalf,nodes)  !!yj nr=mod(igrd1/nrow)
      if( nr.ne.0 ) then
        do n=1,nr
          lonlenndsl(n)=lonlenndsl(n)+1
        enddo
      endif                        
      nm=1
      do n=1,nrow
!yj      do n=1,nodes
        lonstrndsl(n) = nm
        nm = nm + lonlenndsl(n)
!yj        print *,' node lonstrndsl lonlenndsl ',n,lonstrndsl(n),lonlenndsl(n)
      enddo
!! 
      latfull = jgrd1     !!yj latfull=jgrd1
!yj      latfull = latg * 2     !!yj latfull=jgrd1
!yj      lathalf = latg         !!yj no need
      latpart = (jgrd1-1)/nrow+1            !!yj (jgrd1-1)/nrow+1
      latp = jgrd1 / nrow  
      do n=1,nrow
        latlenndsl(n)=latp
      enddo
      nr=mod(jgrd1,nrow)
      if( nr.ne.0 ) then
        do n=1,nr
          latlenndsl(n)=latlenndsl(n)+1
        enddo
      endif                        
!yj      latpart = 0            !!yj (jgrd1-1)/nrow+1
!yj      do n=1,nodes  !!yj no need
!yj        latpart=max(latpart,lats_nodes_a(n))     
!yj        latlenndsl(n) = lats_nodes_a(n)
!yj      enddo

      nm=1
!yj      do n=1,nodes  !!yj n=1,nrow
      do n=1,nrow  !!yj n=1,nrow
        latstrndsl(n) = nm 
        nm = nm + latlenndsl(n)
!yj        print *,' node latstrndsl latlenndsl ',n,latstrndsl(n),latlenndsl(n)
      enddo 
!
      lonlenmax=0
      latlenmax=0
!yj      do n=1,nodes
      do n=1,nrow
        lonlenmax = max(lonlenmax,lonlenndsl(n))
        latlenmax = max(latlenmax,latlenndsl(n))
      enddo
!
      mylonlen = lonlenndsl(myrow+1)
      mylatlen = latlenndsl(myrow+1)
!      if (mype.eq.master) print*,'mylonlen=',mylonlen
!      if (mype.eq.master) print*,'mylatlen=',mylatlen
!      print*,'mype,myrow,mylonlen,mylatlen=',    &
!     &        mype,myrow,mylonlen,mylatlen
!yj      mylonlen = lonlenndsl(me+1)
!yj      mylatlen = latlenndsl(me+1)
!
!      ndslhvar = 5 + ntotal 	! u,v,t,p,dp,tracers
!      ndslvvar = 3 + ntotal	! u,v,t,tracers
!     print *,' ndslfv_init: total variables for ndsl hadv =',ndslhvar
!     print *,' ndslfv_init: total variables for ndsl vadv =',ndslvvar
!!
      return
      end subroutine rndslfv_init

!yj2020------------------------------------------------------------------------>
      subroutine rsopenend_cell_massadv(im,lev,nvars,delt,uc,qq,mass,rdelx)
#include <define.h>
!---------------------------------------------------------------------      
!This is semi-open end cell mass advy, the boundary is exterpolated
!to outer 2 points and then do the fixend cell adv.
!   RSM input im = igrd1 (full grid points)    
!   RSM input lev= levrp (partial grid levels)   
!   nvars = num of variables to be advected
!   delt = timestep (s)
!   uc = wind (u) at cell center (m/s)
!   qq = array of all fields to be advected
!   mass = if 1, do mass conservation, doesn't work in RSM/MSM.
!   rdelx = grid spacing of dx or dy
!
! compute local positive advection with mass conserving
! qq will be advect by uc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
! modified for RSM by YING-Ju CHEN 2018
!
!yj: add extra 2 points outside 1~im+1
!yj:  -1, 0, 1, ..., im+1, p1, p2
!
      use physcons, pi => con_pi
      implicit none
!
      integer   im,lev,nvars,mass
      real	delt
      real	uc(im,lev) 
      real	qq(im,lev,nvars)
      real  rdelx

! index for cell center extrapolation      
      integer m1,m0,p1,p2
      logical loutflow(2)

! array after extrapolation
      real, allocatable :: ue(:),qe(:),de(:)
      real, allocatable :: past(:,:),da(:,:),next(:,:)
      real, allocatable :: var(:), dxfact(:)     &
     &      ,xpast(:),xnext(:),dist(:),interfi(:)
!
      real	step(10), dist_step
      integer  	n,k,i,nst,nstep

!yj ------------v---------------------------v------      
!yj |-o-|-o-|-o-|-o-|-o-|-o-|-o-|- ... -|-o-|-o-|-o-|-o-|
!yj -2  -1  0   1   2   3   4   5  ...  im  +1  +2  +3  +4 <-interfi
!yj   -2  -1  0   1   2   3   4           im  +1  +2  +3   <-c
!yj ------------^---------------------------^------      
!yj   uc  uc  uc  uc  uc  uc  uc   ...    uc
!yj   qq  qq  qq  qq  qq  qq  qq   ...    qq
!yj -----------------------------------------------      

!
! preparations ---------------------------
!yj      print*,'enter rsopenend_cell_massadv'
      m1=-1
      m0=0
      p1=im+1
      p2=im+2
      allocate(ue(m1:p2),qe(m1:p2),de(m1:p2))
      allocate(var(m1:p2),past(m1:p2,nvars),     &
     &    da(m1:p2,nvars),next(m1:p2,nvars),     &
     &  	dxfact(m1:p2),                         &
     &    xpast(m1:p2+1),xnext(m1:p2+1),         &
     &    dist (m1:p2+1),interfi(m1:p2+1))
      ue=0.0
      qe=0.0
      var=0.0
      past=0.0
      da=0.0
      next=0.0
      dxfact=0.0
      xpast=0.0
      xnext=0.0
      dist=0.0
      interfi=0.0
!
      de   = rdelx  !if ds not equal, need to be modified
!
      do k=1,lev
        ue(1:im) = uc(1:im,k)
!
        do i=1,im
          var(i)     = ue(i) * delt !yj: departure points at c
        enddo
!yj:check interior region or exterior region >
! loutflow(1): exterpolation at index=1
! loutflow(2): exterpolation at index=im
        loutflow = .false.
        if (ue(1)  .le. 0.0)   loutflow(1)=.true.
        if (ue(im) .ge. 0.0)   loutflow(2)=.true.
!        if (var(1)  .lt. -0.5*de(1))   loutflow(1)=.true.
!        if (var(im) .gt.  0.5*de(im))  loutflow(2)=.true.
!yj:check interior region or exterior region <
        if (loutflow(1)) then
          var(m0)=2*var(1)-var(2)               !outflow: negative
          var(m1)=2*var(m0)-var(1)               !outflow: negative
        else
          var(m0)=min(var(1),var(2))  !inflow: positive
          var(m1)=var(m0)
        endif
        if (loutflow(2)) then
          var(p1)=2*var(im)-var(im-1)               !outflow: positive
          var(p2)=2*var(p1)-var(im)               !outflow: positive
        else
          var(p1)=max(var(im),var(im-1)) !inflow: negative
          var(p2)=var(p1)
        endif

        do i=m0,p2
          !yj: (un)equal grid, moving distance of interface
          dist(i)=(var(i)*de(i-1)+var(i-1)*de(i)) / (de(i-1)+de(i)) 
        enddo
        dist(m1)=0.0
        dist(p2+1)=0.0
!yj       call rdef_cfl_step(im+5,dist(m1),de(m1),step,nstep)
        call rdef_cfl_step(im+3,dist(m0),de(m0),step,nstep)


        do nst = 1, nstep

        interfi(m1) = 0.0
        xpast(m1)   = interfi(m1)
        xnext(m1)   = interfi(m1)
        do i=m0,p2+1
          interfi(i) = interfi(i-1)+de(i-1)
          dist_step  = dist(i)*step(nst)
          xpast(i)   = interfi(i) - dist_step
          xnext(i)   = interfi(i) + dist_step
        enddo
#ifdef NDSLPPM        
!yj for ppm        
        interfi(:)=-interfi(:)
        xpast(:)=-xpast(:)
        xnext(:)=-xnext(:)
#endif
!        
        if( mass.eq.1 ) then
         do i=m1,p2
          dxfact(i) = (xpast(i+1)-xpast(i)) / (xnext(i+1)-xnext(i))
         enddo
        endif
        do n=1,nvars
          qe(1:im) = qq(1:im,k,n)
          qe(m0)=qe(1)
          qe(m1)=qe(1)
          qe(p1)=qe(im)
          qe(p2)=qe(im)
          past(m1:p2,n) = qe(m1:p2)
        enddo

#ifdef NDSLPPM        
        call rvertical_cell_ppm_intp(interfi(m1),past(m1,1),    &
     &   xpast(m1),da(m1,1),im+4,nvars)
#else
        call rfixend_cell_plm_intp(interfi(m1),past(m1,1),     &
     &   xpast(m1),da(m1,1),im+4,nvars)
#endif

        if( mass.eq.1 ) then
          do n=1,nvars
            da(m1:p2,n) = da(m1:p2,n) * dxfact(m1:p2)
          enddo
        endif

#ifdef NDSLPPM        
        call rvertical_cell_ppm_intp(xnext(m1),da(m1,1),        &
     &   interfi(m1),next(m1,1),im+4,nvars)
#else
        call rfixend_cell_plm_intp(xnext(m1),da(m1,1),         &
     &   interfi(m1),next(m1,1),im+4,nvars)	
#endif

        do n=1,nvars
          qq(1:im,k,n) = next(1:im,n)
        enddo

        enddo   ! for nst


      enddo ! for k

      deallocate(ue,qe,de)
      deallocate(var,past,da,next,    &
     &  	dxfact,xpast,xnext,         &
     &    dist,interfi)
      return
      end subroutine rsopenend_cell_massadv
!-------------------------------------------------------------------
!yj2020------------------------------------------------------------------------<
! ------------------------------------------------------------------------------
      subroutine rsopenend_cell_massadvx(im,lev,nvars,delt,uc,qq,mass)
!yj: This is semi-open end cell mass advy, the boundary is exterpolated
!yj: to outer 2 points and then do the fixend cell adv.
!yj:RSM input im = igrd1      
!yj:RSM input lev = partial lev = levrp      
!yj:nvars = num of variables to be advected
!yj:delt = timestep
!yj:uc = wind (u) at cell center
!yj:qq = array of all fields to be advected
!yj:mass = if 1, do mass conservation, doesn't work in RSM/MSM.
!
! compute local positive advection with mass conserving
! qq will be advect by uc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
! modified for RSM by YING-Ju CHEN 2018
!
!yj: add extra 3 points outside 1~im+1
!yj: -2, -1, 0, 1, ..., im+1, p1, p2, p3
!
      use physcons, pi => con_pi
      implicit none
#include <rscomloc.h>
!
      integer   im,lev,nvars,mass
      real	delt
      real	uc(im,lev) 
      real	qq(im,lev,nvars)
! index for cell center extrapolation      
      integer m4,m3,m2,m1,m0
      integer p1,p2,p3,p4,p5
      logical loutflow(2)
! array after extrapolation
      real, allocatable :: ue(:),qe(:),de(:)
      real, allocatable :: past(:,:),da(:,:),next(:,:)
      real, allocatable :: var(:), dxfact(:)     &
     &      ,xpast(:),xnext(:),dist(:),interfi(:)
!
      real	step(10), dist_step
      integer  	n,k,i,nst,nstep

!yj ------------v---------------------------v------      
!yj |-o-|-o-|-o-|-o-|-o-|-o-|-o-|- ... -|-o-|-o-|-o-|-o-|
!yj -2  -1  0   1   2   3   4   5  ...  im  +1  +2  +3  +4 <-interfi
!yj   -2  -1  0   1   2   3   4           im  +1  +2  +3   <-c
!yj ------------^---------------------------^------      
!yj   uc  uc  uc  uc  uc  uc  uc   ...    uc
!yj   qq  qq  qq  qq  qq  qq  qq   ...    qq
!yj -----------------------------------------------      

!
! preparations ---------------------------
!yj      print*,'enter rsopenend_cell_massadvx'
      m2=-2
      m1=-1
      m0=0
      p1=im+1
      p2=im+2
      p3=im+3
      allocate(ue(m1:p2),qe(m1:p2),de(m1:p2))
      allocate(var(m1:p2),past(m1:p2,nvars),     &
     &    da(m1:p2,nvars),next(m1:p2,nvars),     &
     &  	dxfact(m1:p2),                         &
     &    xpast(m1:p2+1),xnext(m1:p2+1),         &
     &    dist (m1:p2+1),interfi(m1:p2+1))
      ue=0.0
      qe=0.0
      var=0.0
      past=0.0
      da=0.0
      next=0.0
      dxfact=0.0
      xpast=0.0
      xnext=0.0
      dist=0.0
      interfi=0.0
!
      de   = rdelx  !if ds not equal, need to be modified
!
      do k=1,lev
        ue(1:im) = uc(1:im,k)
!
        do i=1,im
          var(i)     = ue(i) * delt !yj: departure points at c
        enddo
!yj:check interior region or exterior region >
! loutflow(1): exterpolation at index=1
! loutflow(2): exterpolation at index=im
        loutflow = .false.
!        if (ue(1)  .lt. 0.0)   loutflow(1)=.true.
!        if (ue(im) .gt. 0.0)   loutflow(2)=.true.
        if (var(1)  .lt. -0.5*de(1))   loutflow(1)=.true.
        if (var(im) .gt.  0.5*de(im))  loutflow(2)=.true.
!yj:check interior region or exterior region <
        if (loutflow(1)) then
          var(m0)=2*var(1)-var(2)               !outflow: negative
!          var(m0)=min(var(1),2*de(1)*var(1)-var(2))  !outflow: negative
          var(m1)=2*var(m0)-var(1)               !outflow: negative
        else
          var(m0)=min(var(1),var(2))  !inflow: positive
          var(m1)=var(m0)
        endif
        if (loutflow(2)) then
          var(p1)=2*var(im)-var(im-1)               !outflow: positive
!          var(p1)=max(var(im),2*de(im)*var(im)-var(im-1)) !outflow: positive
          var(p2)=2*var(p1)-var(im)               !outflow: positive
        else
          var(p1)=max(var(im),var(im-1)) !inflow: negative
          var(p2)=var(p1)
        endif

        do i=m0,p2
          !yj: (un)equal grid, moving distance of interface
          dist(i)=(var(i)*de(i-1)+var(i-1)*de(i)) / (de(i-1)+de(i)) 
        enddo
        dist(m1)=0.0
        dist(p2+1)=0.0
        call rdef_cfl_step(im+5,dist(m1),de(m1),step,nstep)

        do nst = 1, nstep

        interfi(m1) = 0.0
        xpast(m1)   = interfi(m1)
        xnext(m1)   = interfi(m1)
        do i=m0,p2+1
          interfi(i) = interfi(i-1)+de(i-1)
          dist_step  = dist(i)*step(nst)
          xpast(i)   = interfi(i) - dist_step
          xnext(i)   = interfi(i) + dist_step
        enddo
#ifdef NDSLPPM        
!yj for ppm        
        interfi(:)=-interfi(:)
        xpast(:)=-xpast(:)
        xnext(:)=-xnext(:)
#endif
!        
        if( mass.eq.1 ) then
         do i=m1,p2
          dxfact(i) = (xpast(i+1)-xpast(i)) / (xnext(i+1)-xnext(i))
         enddo
        endif
        do n=1,nvars
          qe(1:im) = qq(1:im,k,n)
          qe(m0)=qe(1)
          qe(m1)=qe(1)
          qe(p1)=qe(im)
          qe(p2)=qe(im)
          past(m1:p2,n) = qe(m1:p2)
        enddo

#ifdef NDSLPPM        
        call rvertical_cell_ppm_intp(interfi(m1),past(m1,1),    &
     &   xpast(m1),da(m1,1),im+4,nvars)
#else
        call rfixend_cell_plm_intp(interfi(m1),past(m1,1),     &
     &   xpast(m1),da(m1,1),im+4,nvars)
#endif

        if( mass.eq.1 ) then
          do n=1,nvars
            da(m1:p2,n) = da(m1:p2,n) * dxfact(m1:p2)
          enddo
        endif

#ifdef NDSLPPM        
        call rvertical_cell_ppm_intp(xnext(m1),da(m1,1),        &
     &   interfi(m1),next(m1,1),im+4,nvars)
#else
        call rfixend_cell_plm_intp(xnext(m1),da(m1,1),         &
     &   interfi(m1),next(m1,1),im+4,nvars)	
#endif

        do n=1,nvars
          qq(1:im,k,n) = next(1:im,n)
        enddo

        enddo   ! for nst

      enddo ! for k

      deallocate(ue,qe,de)
      deallocate(var,past,da,next,    &
     &  	dxfact,xpast,xnext,         &
     &    dist,interfi)
      return
      end subroutine rsopenend_cell_massadvx
!-------------------------------------------------------------------
      subroutine rsopenend_cell_massadvy(jm,lev,nvars,delt,vc,qq,mass)
!yj: This is semi-open end cell mass advy, the boundary is exterpolated
!yj: to outer 3 points and then do the fixend cell adv.
!yj      subroutine fixend_cell_massadvy(jm,jmh,lev,nvars,delt,vc,qq,mass)
!yj:RSM input jm = jgrd1      
!yj:RSM input lev = partial lev = levrp      
!
! compute local positive advection with mass conserving
! qq will be advect by vc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
!
      use physcons, pi => con_pi
!
      implicit none
#include <rscomloc.h>
!
      integer   jm,lev,nvars,mass
      real	delt
      real	vc(jm,lev)  !!yj  vc(jgrd1,levpartial)
      real	qq(jm,lev,nvars)
! index for cell center extrapolation      
      integer m4,m3,m2,m1,m0
      integer p1,p2,p3,p4,p5
      logical loutflow(2)
! array after extrapolation
      real, allocatable :: ve(:),qe(:),de(:)
      real, allocatable :: past(:,:),da(:,:),next(:,:)
      real, allocatable :: var(:), dyfact(:)     &
     &      ,ypast(:),ynext(:),dist(:),interfi(:)
!
      real	step(10), dist_step
      integer  	n,k,j,nst,nstep
!
! preparations ---------------------------
!yj      print*,'enter rsopenend_cell_massadvy'
      m2=-2
      m1=-1
      m0=0
      p1=jm+1
      p2=jm+2
      p3=jm+3
      allocate(ve(m1:p2),qe(m1:p2),de(m1:p2))
      allocate(var(m1:p2),past(m1:p2,nvars),     &
     &    da(m1:p2,nvars),next(m1:p2,nvars),     &
     &  	dyfact(m1:p2),                         &
     &    ypast(m1:p2+1),ynext(m1:p2+1),         &
     &    dist (m1:p2+1),interfi(m1:p2+1))
      ve=0.0
      qe=0.0
      var=0.0
      past=0.0
      da=0.0
      next=0.0
      dyfact=0.0
      ypast=0.0
      ynext=0.0
      dist=0.0
      interfi=0.0
!
      de   = rdely  !if ds not equal, need to be modified
!   
!yj: c:cell center
!yj: i:cell interface
      do k=1,lev
        ve(1:jm) = vc(1:jm,k)
!        
        do j=1,jm
          var(j)     = ve(j) * delt  !yj at c
        enddo
!yj:check interior region or exterior region >
! loutflow(1): exterpolation at index=1
! loutflow(2): exterpolation at index=jm
        loutflow = .false.
!        if (var(1)  .lt. 0.0)  loutflow(1)=.true.
!        if (var(jm) .gt. 0.0) loutflow(2)=.true.
        if (var(1)  .lt. -0.5*de(1))   loutflow(1)=.true.
        if (var(jm) .gt.  0.5*de(jm))  loutflow(2)=.true.
!yj:check interior region or exterior region <
        if (loutflow(1)) then
          var(m0)=2*var(1)-var(2)               !outflow: negative
          var(m1)=2*var(m0)-var(1)               !outflow: negative
!          var(m1)=var(m0)
        else
          var(m0)=min(var(1),var(2))  !inflow: positive
          var(m1)=var(m0)           
        endif
        if (loutflow(2)) then
          var(p1)=2*var(jm)-var(jm-1)               !outflow: positive
          var(p2)=2*var(p1)-var(jm)               !outflow: positive
!          var(p2)=var(p1)
        else
          var(p1)=max(var(jm),var(jm-1)) !inflow: negative
          var(p2)=var(p1)
        endif

        do j=m0,p2
          !yj: (un)equal grid, moving distance of interface
          dist(j)=(var(j)*de(j-1)+var(j-1)*de(j)) / (de(j-1)+de(j))
        enddo
        dist(m1)=0.0
        dist(p2+1)=0.0
!
        call rdef_cfl_step(jm+5,dist(m1),de(m1),step,nstep)

        do nst = 1, nstep !yj: loop smaller distance

          interfi(m1) = 0.0
          ypast(m1)   = interfi(m1)
          ynext(m1)   = interfi(m1)
        do j=m0,p2+1
          interfi(j) = interfi(j-1)+de(j-1)
          dist_step  = dist(j)*step(nst) !yj: get smaller distance
          ypast(j)   = interfi(j) - dist_step !yj:gslati(latinterface)->j-1
          ynext(j)   = interfi(j) + dist_step
        enddo
#ifdef NDSLPPM        
!yj for ppm        
        interfi(:)=-interfi(:)
        ypast(:)=-ypast(:)
        ynext(:)=-ynext(:)
#endif
!        
        if( mass.eq.1 ) then
         do j=m1,p2
          dyfact(j) = (ypast(j+1)-ypast(j)) / (ynext(j+1)-ynext(j))
         enddo
        endif
        do n=1,nvars
          qe(1:jm) = qq(1:jm,k,n)
          qe(m0)=qe(1)
          qe(m1)=qe(1)
          qe(p1)=qe(jm)
          qe(p2)=qe(jm)
          past(m1:p2,n) = qe(m1:p2) !yj at c
        enddo
#ifdef NDSLPPM        
        call rvertical_cell_ppm_intp(interfi(m1),past(m1,1),    &
     &   ypast(m1),da(m1,1),jm+4,nvars)
#else
        call rfixend_cell_plm_intp(interfi(m1),past(m1,1),          &
     &                       ypast(m1),da(m1,1),jm+4,nvars)
#endif
        if( mass.eq.1 ) then
          do n=1,nvars
            da(m1:p2,n) = da(m1:p2,n) * dyfact(m1:p2)
          enddo
        endif
#ifdef NDSLPPM        
        call rvertical_cell_ppm_intp(ynext(m1),da(m1,1),        &
     &   interfi(m1),next(m1,1),jm+4,nvars)
#else
        call rfixend_cell_plm_intp(ynext(m1),da(m1,1),              &
     &   interfi(m1),next(m1,1),jm+4,nvars)	
#endif
        do n=1,nvars
          qq(1:jm,k,n) = next(1:jm,n)
        enddo

        enddo   ! for nst

      enddo ! for k

      deallocate(ve,qe,de)
      deallocate(var,past,da,next,    &
     &  	dyfact,ypast,ynext,         &
     &    dist,interfi)
      return
      end subroutine rsopenend_cell_massadvy
!
!-------------------------------------------------------------------------------
      subroutine rfixend_cell_plm_intp(pp,qq,pn,qn,lons,nv)
!
! linear quick version of the same two end points cell interpolation
! so pp(1) = pn(1) and pp(lonp+1) = pn(lonn+1)
!
! pp    location at interfac point as input
! qq    quantity at averaged-cell as input
! pn    location at interface of new grid structure as input
! qn    quantity at averaged-cell as output
! lons  numer of cells for dimension
!
! author : henry.juang@noaa.gov
!
!yj      use gfs_dyn_physcons
!yj      use gfs_dyn_layout1
!
      implicit none
!
      real      pp(lons+1)
      real      qq(lons  ,nv)
      real      pn(lons+1)
      real      qn(lons  ,nv)
      integer   lons,nv
!
      real      hfdp(lons),rdsi(lons)
      real      locf(lons+lons),df(lons+lons-1),dt(lons+lons-1)
      real      qc(lons),dq(lons),qmi(lons),qpi(lons)
      real      valf(lons+lons-1)
      real      rdd(lons)
      real      dd,check,ss
      integer   js(lons+1),i4j(lons+lons)
      integer   i,j,n,ip,in
!
! arrange input array cover output location with cyclic boundary
! condition
!
      locf(1) = pp(1)
      i4j(1)=1
      js(1)=1
      ip=2
      in=2
      do j=2,lons+lons-1 
        if( pp(ip).le.pn(in) ) then
          locf(j) = pp(ip)
          ip = ip + 1
        else
          locf(j) = pn(in)
          js(in) = j
          in = in + 1
        endif
        i4j(j) = ip - 1
      enddo
      locf(lons+lons) = pp(lons+1)
      js(lons+1) = lons+lons
!
! interpolation coefficient
!
      do i=1,lons
        hfdp(i) = 0.5*(pp(i+1)-pp(i))
      enddo
      do i=2,lons
        rdsi(i) = 1.0 / (hfdp(i)+hfdp(i-1))
      enddo
      do i=1,lons
        dd = 0.0
        do j=js(i),js(i+1)-1
          df(j) = locf(j+1) - locf(j)
          dd    = dd    + df(j)
        enddo
        rdd(i) = 1.0 / dd
      enddo
      do j=1,lons+lons-1
        i = i4j(j)
        dt(j) = (locf(j) + 0.5*df(j)) - (pp(i) + hfdp(i))
      enddo
!
! start interpolation by integral of ppm
!
      dq(1) = 0.0
      dq(lons) = 0.0

      do n=1,nv

        do i=1,lons
          qc(i) = qq(i,n)
        enddo
        do i=2,lons
          qmi(i) = (qc(i)-qc(i-1))*rdsi(i)
        enddo
        do i=2,lons-1
          qpi(i) = qmi(i+1)
        enddo
        do i=2,lons-1
          check = qmi(i)*qpi(i)
          if( check.lt.0.0 ) then
            dq(i) = 0.0
          else
            dq(i)  = 0.5*(qmi(i)+qpi(i))
          endif
        enddo

        do j=1,lons+lons-1
          i = i4j(j)
          valf(j) = qc(i) + dq(i)*dt(j)
        enddo

        do i=1,lons
          ss = 0.0
          do j=js(i),js(i+1)-1
            ss = ss + valf(j)*df(j)
          enddo
          qn(i,n) = ss * rdd(i)
        enddo

      enddo
!
      return
      end subroutine rfixend_cell_plm_intp
!
! ------------------------------------------------------------------------
      subroutine rvertical_cell_advect(lons,londim,levs,nvars,           &
     &                                deltim,ssi,wwi,qql,mass)
!
      implicit none

      integer 	londim,levs,nvars,lons,mass
      real 	deltim
!yj      real	ssi(londim,levs+1)
      real	ssi(levs+1)
      real	wwi(londim,levs+1)
      real	qql(londim,levs,nvars)

      real ssii(levs+1)
      real ssid(levs+1),ssia(levs+1)
      real dd(levs+1),ds(levs),step(10),dd_step
      real dsfact(levs), sstmp, dpdt, check
      real rqmm(levs,nvars),rqnn(levs,nvars),rqda(levs,nvars)
      integer km,i,k,n,nst,nstep

        do k=1,levs
          ds(k)=ssi(k)-ssi(k+1)
        enddo
      do i=1,lons                           
!yj        
!        do k=1,levs
!          ds(k)=ssi(i,k)-ssi(i,k+1)
!        enddo
        dd(1) = 0.0
        do k=2,levs
          dd(k) = wwi(i,k) * deltim
        enddo
        dd(levs+1) = 0.0
!cflx   call def_cfl_max (levs+1,dd,ds,nstep)
        call rdef_cfl_step(levs+1,dd,ds,step,nstep)
!
       do nst = 1, nstep
!
!       do k=1,levs+1
!         ssii(k)=-ssi(i,k)
!         ssid(k)=-ssi(i,k)+dd(k)
!         ssia(k)=-ssi(i,k)-dd(k)
!       enddo
!hmhj give direction for value larger with k larger
        do k=1,levs+1
          dd_step= dd(k)*step(nst)
!yj for plm          
          ssii(k)=-ssi(k)
          ssid(k)=-ssi(k)-dd_step
          ssia(k)=-ssi(k)+dd_step
!yj          ssid(k)=-ssi(i,k)+dd_step
!yj          ssia(k)=-ssi(i,k)-dd_step
#ifdef NDSLPPMV        
!yj for ppm          
          ssii(k)=ssi(k)
          ssid(k)=ssi(k)+dd_step
          ssia(k)=ssi(k)-dd_step
#endif
!          
        enddo
        if( mass.eq.1 ) then
          do k=1,levs
            dsfact(k)=(ssid(k)-ssid(k+1))/(ssia(k)-ssia(k+1))
          enddo
        endif
!
        do n=1,nvars
          do k=1,levs
            rqmm(k,n) = qql(i,k,n)
          enddo
        enddo
#ifdef NDSLPPMV        
        call rvertical_cell_ppm_intp(ssii,rqmm,ssid,rqda,levs,nvars)
#else
        call rfixend_cell_plm_intp(ssii,rqmm,ssid,rqda,levs,nvars)
#endif
        if( mass.eq.1 ) then
          do n=1,nvars
            do k=1,levs
              rqda(k,n) = rqda(k,n) * dsfact(k)
            enddo
          enddo
        endif
#ifdef NDSLPPMV        
        call rvertical_cell_ppm_intp(ssia,rqda,ssii,rqnn,levs,nvars)
#else
        call rfixend_cell_plm_intp(ssia,rqda,ssii,rqnn,levs,nvars)
#endif
        do n=1,nvars
          do k=1,levs
            qql(i,k,n)=rqnn(k,n)
          enddo
        enddo
!
       enddo
!
      enddo

!yj      print*,'finish vertical_cell_advect'
      return
      end subroutine rvertical_cell_advect
! 
! ------------------------------------------------------------------------
      subroutine rvertical_cell_ppm_intp(pp,qq,pn,qn,levs,nvars)
!
! mass conservation in vertical interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in vertical with mass conservation is under the condition that
! pp(1)=pn(1), pp(levs+1)=pn(levs+1)
!
! pp    pressure at interfac level as input
! qq    quantity at layer as input
! pn    pressure at interface of new grid structure as input
! qn    quantity at layer as output
! levs  numer of verical layers
!
! author : henry.juang@noaa.gov
!
      implicit none
!
      real      pp(levs+1)
      real      qq(levs,nvars)
      real      pn(levs+1)
      real      qn(levs,nvars)
      integer   levs,nvars
!
      real      massm,massc,massp,massbot,masstop
      real      qmi(levs,nvars),qpi(levs,nvars)
      real      dql(nvars),dqh(nvars)
      real      hh(levs)
      real      dqi,dqimax,dqimin,dqmono(levs,nvars)
      real      tl,tl2,tl3,tlp,tlm,tlc
      real      th,th2,th3,thp,thm,thc
      real      dpp,dqq,c1,c2
      integer   i,k, kl, kh, kk, kkl, kkh,n
      integer, parameter :: mono=1
!
      if( pp(1).ne.pn(1) .or. pp(levs+1).ne.pn(levs+1) ) then
        print *,' Error in vertical_cell_ppm_intp for domain values '
        print *,' i pp1 pn1 ppt pnt ',i,     &
     &          pp(1),pn(1),pp(levs+1),pn(levs+1)
        call abort
      endif
!
! prepare thickness for grid
!
      do k=1,levs
        hh(k) = pp(k+1)-pp(k)
      enddo
!
! prepare location with monotonic concerns
!
      do n=1,nvars
        massbot=(3.*hh(1)+hh(2))*qq(1,n)-2.*hh(1)*qq(2,n)
        massm = massbot/(hh(1)+hh(2))
        massc = qq(1  ,n)
        massp = qq(1+1,n)
        dqi = 0.25*(massp-massm)
        dqimax = max(massm,massc,massp) - massc
        dqimin = massc - min(massm,massc,massp)
        dqmono(1,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
        do k=2,levs-1
          massp = qq(k+1,n)
          massc = qq(k  ,n)
          massm = qq(k-1,n)
          dqi = 0.25*(massp-massm)
          dqimax = max(massm,massc,massp) - massc
          dqimin = massc - min(massm,massc,massp)
          dqmono(k,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
        enddo
        masstop=(3.*hh(levs)+hh(levs-1))*qq(levs,n)  &
     &             -2.*hh(levs)*qq(levs-1,n)
        massp = masstop/(hh(levs)+hh(levs-1))
        massc = qq(levs  ,n)
        massm = qq(levs-1,n)
        dqi = 0.25*(massp-massm)
        dqimax = max(massm,massc,massp) - massc
        dqimin = massc - min(massm,massc,massp)
        dqmono(levs,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
!
! compute value at interface with momotone
!
        do k=2,levs
          qmi(k,n)=(qq(k-1,n)*hh(k)+qq(k,n)*hh(k-1))/(hh(k)+hh(k-1)) &
     &       +(dqmono(k-1,n)-dqmono(k,n))/3.0
        enddo
        do k=1,levs-1
          qpi(k,n)=qmi(k+1,n)
        enddo
        qmi(1,n)=qq(1,n)
        qpi(1,n)=qq(1,n)
        qmi(levs,n)=qq(levs,n)
        qpi(levs,n)=qq(levs,n)
      enddo
!
! do monotonicity
!
      if( mono.eq.1 ) then
        do n=1,nvars
        do k=1,levs
          c1=qpi(k,n)-qq(k,n)
          c2=qq(k,n)-qmi(k,n)
          if( c1*c2.le.0.0 ) then
            qmi(k,n)=qq(k,n)
            qpi(k,n)=qq(k,n)
          endif
        enddo
        do k=1,levs
          c1=(qpi(k,n)-qmi(k,n))*(qq(k,n)-0.5*(qpi(k,n)+qmi(k,n)))
          c2=(qpi(k,n)-qmi(k,n))*(qpi(k,n)-qmi(k,n))/6.
          if( c1.gt.c2 ) then
            qmi(k,n)=3.*qq(k,n)-2.*qpi(k,n)
          else if( c1.lt.-c2 ) then
            qpi(k,n)=3.*qq(k,n)-2.*qmi(k,n)
          endif
        enddo
        enddo
      endif
!
! start interpolation by integral of ppm spline
!
      kkl = 1
      tl=0
      do n=1,nvars
        dql(n)=0.0
      enddo

      do k=1,levs

        kl = k
        kh = k + 1
! find kkh
        do kk=kkl+1,levs+1
          if( pn(kh).ge.pp(kk) ) then
            kkh = kk-1
            go to 100
          endif
        enddo
        print *,' Error in vertical_cell_ppm_intp for no lev found '
        print *,' i kh kl ',i,kh,kl
        print *,' pn ',(pn(kk),kk=1,levs+1)
        print *,' pp ',(pp(kk),kk=1,levs+1)
        call abort
 100    continue
        th=(pn(kh)-pp(kkh))/hh(kkh)
        th2=th*th
        th3=th2*th
        thp = th3-th2
        thm = th3-2.*th2+th
        thc = -2.*th3+3.*th2
        do n=1,nvars
          dqh(n)=thp*qpi(kkh,n)+thm*qmi(kkh,n)+thc*qq(kkh,n)
        enddo
! mass interpolate
        if( kkh.eq.kkl ) then
          do n=1,nvars
            qn(k,n) = (dqh(n)-dql(n))/(th-tl)
          enddo
        else if( kkh.gt.kkl ) then
          dpp  = (1.-tl)*hh(kkl) + th*hh(kkh)
          do kk=kkl+1,kkh-1
            dpp = dpp + hh(kk)
          enddo
          do n=1,nvars
            dql(n) = qq(kkl,n)-dql(n)
            dqq  = dql(n)*hh(kkl) + dqh(n)*hh(kkh)
            do kk=kkl+1,kkh-1
              dqq = dqq + qq(kk,n)*hh(kk)
            enddo
            qn(k,n) = dqq / dpp
          enddo
        else
          print *,' Error in vertical_cell_ppm_intp for lev messed up '
          print *,' i kh kl ',i,kh,kl
          print *,' pn ',(pn(kk),kk=1,levs+1)
          print *,' pp ',(pp(kk),kk=1,levs+1)
          call abort
        endif
! next one
        kkl = kkh
        tl  = th
        do n=1,nvars
          dql(n) = dqh(n)
        enddo

      enddo	! end of k loop
!
      return
      end subroutine rvertical_cell_ppm_intp
!!
! ------------------------------------------------------------------------
!!
      subroutine def_cfl_max (im,dist,del,nstep)
!
! compute the deformation cfl condition
! select the maxima value of the deformation CFL and provide step to
! avoid it.
!
!
      implicit none
!
      integer im,nstep
      real dist(im),del(im-1)

!yj: im is all interface (igrd1+1 or jgrd1+1)      
!yj: dist is the moving distance of interface 
!yj: del is delx or dely
!
! local
      integer   n
      real      rstep,check,check_max,check_point
!
      check_point=1.00
      nstep = 1
      check_max = 0.0
      check_loop: do n=1,im-1
        check = abs ( (dist(n+1)-dist(n))/del(n) )
        if( check.lt.check_point ) then
          cycle check_loop
        else
          check_max = max( check_max, check )
        endif
      enddo check_loop
      if(check_max.ge.check_point) then
        nstep = int(check_max/check_point) + 1
        print *,' max def_cfl ',check_max,' needs ',nstep,'steps in ',im
        rstep = 1. / nstep
        do n=1,im
          dist(n) = dist(n) * rstep
        enddo
      endif
!
      return
      end subroutine def_cfl_max

!!
! ------------------------------------------------------------------------
!!
      subroutine rdef_cfl_step (im,dist,del,step,nstep)
!
! compute the deformation cfl condition
! select the maxima value of the deformation CFL and provide step to
! avoid it.
!
!yj: im is all interface (igrd1+1 or jgrd1+1)      
!yj: dist is the moving distance of interface 
!yj: del is delx or dely
!
      implicit none
!
      integer im,nstep
      real dist(im),del(im-1),step(10)
! local
      integer   n
      real      rstep,check,check_max,check_point
      real      safe_step,last_step
!
      check_point=1.00
      safe_step=0.99
      nstep = 1
      step(1) = 1.0

      check_max = 0.0
      check_loop: do n=1,im-1
        check = abs ( (dist(n+1)-dist(n))/del(n) )
        if( check.lt.check_point ) then
          cycle check_loop
        else
          check_max = max( check_max, check )
          print*,'n,dist(n+1),dist(n),del(n),check=', &
     &     n,dist(n+1),dist(n),del(n),check
        endif
      enddo check_loop
      if(check_max.ge.check_point) then
        nstep = int(check_max/safe_step) + 1
        print *,' max def_cfl ',check_max,' needs ',nstep,'steps in ',im
        rstep =  safe_step / check_max 
        do n=1,nstep-1
          step(n) = rstep
        enddo
        last_step = 1. - ( nstep - 1 ) * rstep
        step(nstep) = last_step
      endif
!
      return
      end subroutine rdef_cfl_step

!!
! ------------------------------------------------------------------------
      subroutine rmymaxmin(a,im,ix,kx,ch)
        implicit none
        real :: fmax
        real :: fmin
        integer :: i
        integer :: im
        integer :: ix
        integer :: k
        integer :: kx

      real a(ix,kx)
      character*(*) ch
      do k=1,kx
        fmin=a(1,k)
        fmax=a(1,k)
        do i=1,im
          fmin=min(fmin,a(i,k))
          fmax=max(fmax,a(i,k))
        enddo
        print *,' max=',fmax,' min=',fmin,' at k=',k,' for ',ch
      enddo
      return
      end subroutine rmymaxmin
! ------------------------------------------------------------------------

      subroutine extp4_back( qi, ds, qo)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! YJ CHEN 2018 Jul 2
!    backward exterpolation 
! Arguments:
!    input
!      qi: real array to be exterpolated (ip)
!      ip: integer, involved points  =4
!      ds: integer array, distance between every q point
!    output
!      qo: real array, after exterpolated (2*ip-1)
!
! Examples:
!   ip=4
!
!                     V this is the boundary inner cell of qi
!      di-2  di-1  di0   di1   di2   di3
!   ds-2  ds-1  ds0   ds1   ds2   ds3   ds4
!   X     X     X     *     O     O     O      <-qo and qi
!   -2    -1    0     1     2     3     4
!      w     w     w     w     w     w         <-q1 (tendency of q, q')
!      -2    -1    0     1     2     3
!         w     w     w     w     w            <-q2 (q'')
!         -2    -1    0     1     2
!            w     w     w     w               <-q3 (q''')
!            -2    -1    0     1
!               w     w     w                  <-q4 (q'''')
!               -2    -1    0
!                        w                                  
!                        -1  
!                     <------  exterpolation direction 
!
!   *: boundary point of q
!   x: points needed to be exterpolated
!   o: inner points of q
!   w: tendency of previous line 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none
!ip=4
      integer ip, m0, m1, m2
      parameter(ip=4,m0=0,m1=-1,m2=-2)
      integer qloc
      real qi(1:ip)  !(1:4)  !(x,x,x,1,2,3,4)!(1,2,3,4,0,-1,-2)
      real qo(m2:ip) !(-2:4)
      real ds(m2:ip) !(-2,4)
      real di(m2:ip-1) !(-2:3)

! local array for exterpolation
      real, allocatable :: qten(:,:)
! local loop index      
      integer i,j

      allocate( qten(m2:ip,0:ip) ) 

! start
      qten(1:ip,0)=qi(1:ip)
      di(m2:ip-1)=0.5*(ds(m1:ip)+ds(m2:ip-1))

! get tendency q1 
      qten(1:3,1) = (qten(2:4,0) - qten(1:3,0) )/ di(1:3)
! get tendancy q2
      qten(1:2,2) = (qten(2:3,1) - qten(1:2,1))/ ds(2:3)

! exterpolation

!for q(0)
      qloc=0
      qten(qloc+1,3)=(qten(qloc+2,2)-qten(qloc+1,2))/di(2)
      qten(qloc,2)=qten(qloc+1,2)-(qten(qloc+1,3)*di(1))
      qten(qloc,1)=qten(qloc+1,1)-qten(qloc,2)*ds(1)
      qten(qloc,0)=qten(qloc+1,0)-qten(qloc,1)*di(0)

      qten(qloc,3)=(qten(qloc+1,2)-qten(qloc,2))/di(1)
!for q(-1)
      qloc=m1
      qten(qloc+1,4)=(qten(qloc+2,3)-qten(qloc+1,3))/ds(2)
      qten(qloc,3)=qten(qloc+1,3)-qten(qloc+1,4)*ds(1)
      qten(qloc,2)=qten(qloc+1,2)-qten(qloc,3)*di(0)
      qten(qloc,1)=qten(qloc+1,1)-qten(qloc,2)*ds(0)
      qten(qloc,0)=qten(qloc+1,0)-qten(qloc,1)*di(-1)

      qten(qloc,4)=(qten(qloc+1,3)-qten(qloc,3))/ds(1)
!for q(-2)
      qloc=m2
      qten(qloc,4)=qten(qloc+1,4)-(qten(qloc+2,4)-qten(qloc+1,4))/di(1)*di(0)
      qten(qloc,3)=qten(qloc+1,3)-qten(qloc,4)*ds(0)
      qten(qloc,2)=qten(qloc+1,2)-qten(qloc,3)*di(-1)
      qten(qloc,1)=qten(qloc+1,1)-qten(qloc,2)*ds(-1)
      qten(qloc,0)=qten(qloc+1,0)-qten(qloc,1)*di(-2)


      qo(m2:ip) = qten(m2:ip,0)

      deallocate( qten ) 

      return
      end subroutine extp4_back
! ------------------------------------------------------------------------

      subroutine extp4_forw( qi, ds, qo)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! YJ CHEN 2018 Jul 2
!    forward exterpolation 
! Arguments:
!    input
!      qi: real array to be exterpolated (ip)
!      ip: integer, involved points  =4
!      ds: integer array, distance between every q point
!    output
!      qo: real array, after exterpolated (2*ip-1)
!
! Examples:
!   ip=4
!
!                     V this is the boundary point of qi
!      di1   di2   di3   di4   dip1  dip2  
!   ds1   ds2   ds3   ds4   dsp1  dsp2  dsp3
!   X     X     X     *     O     O     O      <-qo and qi
!   1     2     3     4     p1    p2    p3
!      w     w     w     w     w     w         <-q1 (tendency of q, q')
!      2     3     4     p1    p2    p3
!         w     w     w     w     w            <-q2 (q'')
!         3     4     p1    p2    p3
!            w     w     w     w               <-q3 (q''')
!            4     p1    p2    p3
!               w     w     w                  <-q4 (q'''')
!               p1    p2    p3
!                  w                                        
!                  p2        
! exterpolation direction ------->
!
!   *: boundary point of q
!   x: points needed to be exterpolated
!   o: inner points of q
!   w: tendency of previous line 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none
!ip=4
      integer ip, p1, p2, p3
      parameter(ip=4,p1=5,p2=6,p3=7)
      integer qloc
      real qi(1:ip)  
      real qo(1:p3) 
      real ds(1:p3)
      real di(1:p2) 

! local array for exterpolation
      real, allocatable :: qten(:,:)
! local loop index      
      integer i,j

      allocate( qten(1:p3,0:ip) ) 

! start
      qten(1:ip,0)=qi(1:ip)
      di(1:p2)=0.5*(ds(2:p3)+ds(1:p2))

! get tendency q1 
      qten(2:4,1) = (qten(2:4,0) - qten(1:3,0) )/ di(1:3)
! get tendancy q2
      qten(3:4,2) = (qten(3:4,1) - qten(2:3,1))/ ds(2:3)

! exterpolation

!for q(p1)
      qloc=p1
      qten(qloc-1,3)=(qten(qloc-1,2)-qten(qloc-2,2))/di(2)

      qten(qloc,2)=qten(qloc-1,2)+(qten(qloc-1,3)*di(3))
      qten(qloc,1)=qten(qloc-1,1)+qten(qloc,2)*ds(4)
      qten(qloc,0)=qten(qloc-1,0)+qten(qloc,1)*di(4)

      qten(qloc,3)=(qten(qloc,2)-qten(qloc-1,2))/di(3)
!for q(p2)
      qloc=p2
      qten(qloc-1,4)=(qten(qloc-1,3)-qten(qloc-2,3))/ds(3)
      qten(qloc,3)=qten(qloc-1,3)+qten(qloc-1,4)*ds(4)
      qten(qloc,2)=qten(qloc-1,2)+qten(qloc,3)*di(4)
      qten(qloc,1)=qten(qloc-1,1)+qten(qloc,2)*ds(p1)
      qten(qloc,0)=qten(qloc-1,0)+qten(qloc,1)*di(p1)

      qten(qloc,4)=(qten(qloc,3)-qten(qloc-1,3))/ds(4)
!for q(p3)
      qloc=p3
      qten(qloc,4)=qten(qloc-1,4)+(qten(qloc-1,4)-qten(qloc-2,4))/di(3)*di(4)
      qten(qloc,3)=qten(qloc-1,3)+qten(qloc,4)*ds(p1)
      qten(qloc,2)=qten(qloc-1,2)+qten(qloc,3)*di(p1)
      qten(qloc,1)=qten(qloc-1,1)+qten(qloc,2)*ds(p2)
      qten(qloc,0)=qten(qloc-1,0)+qten(qloc,1)*di(p2)


      qo(1:p3) = qten(1:p3,0)

      deallocate( qten ) 

      return
      end subroutine extp4_forw

! ------------------------------------------------------------------------
!yj      subroutine glob_bdy_back()
!yj      return
!yj      end subroutine glob_bdy_back
! ------------------------------------------------------------------------
!yj      subroutine glob_bdy_forw()
!yj      return
!yj      end subroutine glob_bdy_forw
! ------------------------------------------------------------------------
      subroutine rndslfv_final
      use mod_rndsl
        deallocate ( lonstrndsl, lonlenndsl )
        deallocate( latstrndsl, latlenndsl )
      return
      end subroutine rndslfv_final




!-----for index (yj 2019.02.13)------------------------------------------      
      subroutine rndslfv_monoadvh2_ind(grid_gr,levrp,k1,k2,lotgr,uulon,vvlon,deltim)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
      use machine, only: kind_grid => kind_io8
      use mod_rndsl
      implicit none
#include <define.h>
#include <npesi.h>
#include <comrmpi.h>

      
!for mpi
!yj-same as rloopmp.f >
!yj      integer kumk,kvmk,ktmk,kpmk,kqmk,kgzmk
!yj#ifdef NONHYD
!yj      integer kpnmk,ktnmk,kwmk
!yj#endif
!NDSL:synmk>
!yj      parameter(kumk=1,                                                  &
!yj     &          kvmk=levrp+kumk,                                         &
!yj     &          ktmk=levrp+kvmk,                                         &
!yj     &          kpmk=levrp+ktmk,                                         &
!yj     &          kqmk=levrp+kpmk,                                         & 
!yj     &          kgzmk=levhp+kqmk                                         &
!yj#ifdef NONHYD     
!yj     &         ,kpnmk=levrp+kgzmk                                        &
!yj     &         ,ktnmk=levrp+kpnmk                                        &
!yj     &         ,kwmk=levrp+ktnmk !kwmt:levrp+1 (midpoint: 1~levrp+1)     &
!yj#endif
!yj     &          )
!NDSL:synmk<
!yj#ifdef NONHYD
!yj      parameter(lotgr=1+8*levrp+levhp)
!yj#else
!yj      parameter(lotgr=  5*levrp+levhp)
!yj#endif
!yj-same as rloopmp.f <
      
      integer, intent(in):: levrp,k1,k2,lotgr
      integer:: nkk,nvar

      real(kind=kind_grid)::grid_gr(lonfull,lotgr,latpart)
      real,   intent(in):: deltim
   
      real, intent(in):: uulon(lonfull,levrp,latpart)
      real, intent(in):: vvlon(lonfull,levrp,latpart)
      real,allocatable::qqlon(:,:,:),rrlon(:,:,:)

      real      vvlat(latfull,levrp,lonpart)
      real,allocatable::qqlat(:,:,:),rrlat(:,:,:)
      real      rma, rm2a, rdt2, rkt, pi, cons0, cons1

      logical   lprint

      integer mono,mass
      integer ilan,i,j,n,k,kk,lon,lan,lat,lons_lat,jlonf,irc
!
       nkk=k2-k1+1
! this will do k1 ~ k2  !(k2 include)       
       nvar=(nkk-1)/levrp+1
      allocate(                         &
     &  qqlon(lonfull,nkk,latpart),   & 
     &  rrlon(lonfull,nkk,latpart),   &
     &  qqlat(latfull,nkk,lonpart),   &
     &  rrlat(latfull,nkk,lonpart)    &
     &  )
!      
      lprint = .false.

!yj      print *,' enter ndslfv_monoadvh '
!yj      if( lprint ) print *,' enter ndslfv_monoadvh '
!
      mono  = 1
      mass  = 0
!
      rdt2 = 0.5 / deltim      
!
! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
!
#ifdef OPENMP
!$omp parallel do schedule(dynamic,1) private(lan)
#endif
      do lan=1,mylatlen
!
        qqlon(:,:,lan) = grid_gr(:,k1:k2,lan)

        rrlon(:,:,lan) = qqlon(:,:,lan)
!
! first set positive advection in east-west direction 
!
        call rsopenend_cell_massadvx_ind(lonfull,levrp,nvar,deltim   &
     &                    ,uulon(1,1,lan),rrlon(1,1,lan),mass)

      enddo  !end do lan=1,mylatlen


! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat

       if( lprint ) print *,' ndslfv_advect transport from we to ns '

       call rpara_we2ns(vvlon,vvlat,levrp)
       call rpara_we2ns(rrlon,rrlat,nkk)
       call rpara_we2ns(qqlon,qqlat,nkk)

       if( lprint ) then
       print *,' ------------ after we2ns ---------------------- '
       do lon=1,mylonlen
        print *,'  lon=',lon
!yj        call mymaxmin(vvlat(1,1  ,lon),latfull,latfull,1,' we2ns v')
!yj        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' we2ns r')
!yj        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' we2ns q')
       enddo
       endif
!
! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

       if( lprint ) then
       print *,' ndslfv_advect adv loop in y '
       print *,' mylonlen=',mylonlen
       endif

#ifdef OPENMP
!$omp parallel do schedule(dynamic,1) private(lon)
#endif
       do lon=1,mylonlen
! 
        if( lprint ) print *,' lon=',lon

! first set advection in north-south direction in great circle through two poles
!
        call rsopenend_cell_massadvy_ind(latfull,levrp,nvar,deltim,    &
     &                   vvlat(1,1,lon),rrlat(1,1,lon),mass)

        if( lprint ) then
        print *,' ------------------------------------------- '
!yj        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry u ')
!yj        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry v ')
!yj        call mymaxmin(rrlat(1,ktt,lon),latfull,latfull,1,' advry t ')
!yj        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' advry q ')
!yj        call mymaxmin(rrlat(1,kdp,lon),latfull,latfull,1,' advry dp ')
!yj        call mymaxmin(rrlat(1,kp2,lon),latfull,latfull,1,' advry p2 ')
        endif
!
! second set advection in north-south direction in great circle through two poles
!
        call rsopenend_cell_massadvy_ind(latfull,levrp,nvar,deltim,    &
     &                   vvlat(1,1,lon),qqlat(1,1,lon),mass)

        if( lprint ) then
        print *,' ------------------------------------------- '
!yj        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy u ')
!yj        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy v ')
!yj        call mymaxmin(qqlat(1,ktt,lon),latfull,latfull,1,' advqy t ')
!yj        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' advqy q ')
!yj        call mymaxmin(qqlat(1,kdp,lon),latfull,latfull,1,' advqy dp ')
!yj        call mymaxmin(qqlat(1,kp2,lon),latfull,latfull,1,' advqy p2 ')
        print *,' done with y at lon=',lon
        endif

       enddo  !lon=1,mylonlen
!
! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory 
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon

       if( lprint ) print *,' ndslfv_advect transport from ns to we '

       call rpara_ns2we(rrlat,rrlon,nkk)
       call rpara_ns2we(qqlat,qqlon,nkk)

       if(lprint) then
       print *,' ------------ after ns2we ---------------------- '
!yj       do lan=1,lats_node_a
!yj        print *,'  lan=',lan
!yj        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' ns2we r')
!yj        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' ns2we q')
!yj       enddo
       endif

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!      print *,' ndslfv_advect adv loop in x for last '

#ifdef OPENMP
!$omp parallel do schedule(dynamic,1) private(lan), &
!$omp&             private(i,k)
#endif

      do lan=1,mylatlen

! second set advection in x for the second of the pair
        call rsopenend_cell_massadvx_ind(lonfull,levrp,nvar,deltim,       &
     &                   uulon(1,1,lan),qqlon(1,1,lan),mass)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!yj        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' adv x q ')
!        endif


        do k=1,nkk
          do i=1,lonfull
            qqlon(i,k,lan) = 0.5 * ( qqlon(i,k,lan) + rrlon(i,k,lan) )
          enddo
        enddo

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!yj        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu u ')
!yj        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu v ')
!yj        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' do redu t ')
!yj        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' do redu q ')
!yj        call mymaxmin(qqlon(1,kdp,lan),lons_lat,lonfull,1,' do redu dp')
!yj        call mymaxmin(qqlon(1,kp2,lan),lons_lat,lonfull,1,' do redu p2')
!        print *,' finish horizonatal advection at lan=',lan
!        endif
!yj180829>
!yj      enddo  ! do lan=1, mylatlen
!yj        print*,'after o3=',qqlon(:,5*levrp+1,:)
!yj      do lan=1,mylatlen
!yj180829<
! tendency at n
!yj        grid_gr(:,k1:k2,lan) = (qqlon(:,:,lan)-grid_gr(:,k1:k2,lan))*rdt2
        grid_gr(:,k1:k2,lan) = qqlon(:,:,lan)
   
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        print *,' finish horizonatal advection at lan=',lan
!        endif
!
      enddo  ! do lan=1, mylatlen
!yj       print*,'in ndsl_monoadvh2 tendency    &
!yj     &  u(8,1:levrp,8)=',grid_gr(8,1:levrp,8)  
!yj      print*,'in ndsl_monoadvh2 tendency  &
!yj     & Qs(8,1:levrp,8)=',grid_gr(8,levrp*3+1:levrp*4,8)
!yj      print*,'in ndsl_monoadvh2 tendency  &
!yj     & spfh(8,1:levrp,8)=',grid_gr(8,levrp*4+1:levrp*5,8)

! 
! ===============================
      deallocate(qqlon,rrlon,qqlat,rrlat)
!
      return
      end



!!yj----------------------------------------------------------------      
      subroutine rsopenend_cell_massadvx_ind(im,lev,nvars,delt,uc,qq,mass)
!yj: This is semi-open end cell mass advy, the boundary is exterpolated
!yj: to outer 3 points and then do the fixend cell adv.
!yj:RSM input im = igrd1      
!yj:RSM input lev = partial lev = levrp      
!yj:nvars = num of variables to be advected
!yj:delt = timestep
!yj:uc = wind (u) at cell center
!yj:qq = array of all fields to be advected
!yj:mass = if 1, do mass conservation?
!
! compute local positive advection with mass conserving
! qq will be advect by uc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
! modified for RSM by YING-Ju CHEN 2018
!
!yj: add extra 3 points outside 1~im+1
!yj: -2, -1, 0, 1, ..., im+1, p1, p2, p3
!
      use physcons, pi => con_pi
      implicit none
#include <rscomloc.h>
!
      integer   im,lev,nvars,mass
      real	delt
      real	uc(im,lev) 
      real	qq(im,lev,nvars)
! index for cell center extrapolation      
      integer m4,m3,m2,m1,m0
      integer p1,p2,p3,p4,p5
      logical loutflow(2)
! array after extrapolation
      real, allocatable :: ue(:),qe(:),de(:)
      real, allocatable :: past(:,:),da(:,:),next(:,:)
      real, allocatable :: var(:), ds(:), dxfact(:)     &
     &      ,xpast(:),xnext(:),dist(:),interfi(:)
!
      real	step(10), dist_step

!yj ------------v---------------------------v------      
!yj |-o-|-o-|-o-|-o-|-o-|-o-|-o-|- ... -|-o-|-o-|-o-|-o-|
!yj -2  -1  0   1   2   3   4   5  ...  im  +1  +2  +3  +4 <-interfi
!yj   -2  -1  0   1   2   3   4           im  +1  +2  +3   <-c
!yj ------------^---------------------------^------      
!yj   uc  uc  uc  uc  uc  uc  uc   ...    uc
!yj   qq  qq  qq  qq  qq  qq  qq   ...    qq
!yj -----------------------------------------------      

      integer  	n,k,i,nst,nstep
!
! preparations ---------------------------
!yj      print*,'enter rsopenend_cell_massadvx'
      m2=-2
      m1=-1
      m0=0
      p1=im+1
      p2=im+2
      p3=im+3
      allocate(ue(m1:p2),qe(m1:p2),de(m1:p2))
      allocate(var(m1:p2),past(m1:p2,nvars),     &
     &    da(m1:p2,nvars),next(m1:p2,nvars),     &
     &  	dxfact(m1:p2),                         &
     &    xpast(m1:p2+1),xnext(m1:p2+1),         &
     &    dist (m1:p2+1),interfi(m1:p2+1))
      ue=0.0
      qe=0.0
      var=0.0
      past=0.0
      da=0.0
      next=0.0
      dxfact=0.0
      xpast=0.0
      xnext=0.0
      dist=0.0
      interfi=0.0
!
      de   = rdelx  !if ds not equal, need to be modified
!
      do k=1,lev
! extrapolation for left and right domain 
        ue(1:im) = uc(1:im,k)
!
        do i=1,im
          var(i)     = ue(i) * delt !yj: departure points at c
        enddo
!yj:check interior region or exterior region >
! loutflow(1): exterpolation at index=1
! loutflow(2): exterpolation at index=im
        loutflow = .false.
!        if (var(1)  .lt. 0.0)   loutflow(1)=.true.
!        if (var(im) .gt. 0.0)   loutflow(2)=.true.
        if (var(1)  .lt. -0.5*de(1))   loutflow(1)=.true.
        if (var(im) .gt.  0.5*de(im))  loutflow(2)=.true.
!yj:check interior region or exterior region <
        if (loutflow(1)) then
          var(m0)=2*var(1)-var(2)               !outflow: negative
          var(m1)=2*var(m0)-var(1)               !outflow: negative
        else
          var(m0)=min(var(1),var(2))  !inflow: positive
          var(m1)=var(m0)
        endif
        if (loutflow(2)) then
          var(p1)=2*var(im)-var(im-1)               !outflow: positive
          var(p2)=2*var(p1)-var(im)               !outflow: positive
        else
          var(p1)=max(var(im),var(im-1))
          var(p2)=var(p1)
        endif

        do i=m0,p2
          !yj:moving dist. of interface
          dist(i)=(var(i)*de(i-1)+var(i-1)*de(i)) / (de(i-1)+de(i)) 
        enddo
        dist(m1)=0.0
        dist(p2+1)=0.0
!        
        call rdef_cfl_step(im+5,dist(m1),de(m1),step,nstep)

        do nst = 1, nstep

        interfi(m1) = 0.0
        xpast(m1)   = interfi(m1)
        xnext(m1)   = interfi(m1)
        do i=m0,p2+1
          interfi(i) = interfi(i-1)+de(i-1)
          dist_step  = dist(i)*step(nst)
          xpast(i)   = interfi(i) - dist_step
          xnext(i)   = interfi(i) + dist_step
        enddo
#ifdef NDSLPPM        
!yj for ppm        
        interfi(:)=-interfi(:)
        xpast(:)=-xpast(:)
        xnext(:)=-xnext(:)
#endif
        if( mass.eq.1 ) then
         do i=m1,p2
          dxfact(i) = (xpast(i+1)-xpast(i)) / (xnext(i+1)-xnext(i))
         enddo
        endif
        do n=1,nvars
          qe(1:im) = qq(1:im,k,n)
          qe(m0)=0.0
          qe(m1)=0.0
          qe(p1)=0.0
          qe(p2)=0.0
          past(m1:p2,n) = qe(m1:p2)
        enddo

#ifdef NDSLPPM        
!        call vertical_cell_ppm_intp_ind(interfi(m1),past(m1,1),    &
        call rvertical_cell_ppm_intp(interfi(m1),past(m1,1),    &
     &   xpast(m1),da(m1,1),im+4,nvars)
#else
        call rfixend_cell_plm_intp_ind(interfi(m1),past(m1,1),     &
     &   xpast(m1),da(m1,1),im+4,nvars)
#endif

        if( mass.eq.1 ) then
          do n=1,nvars
            da(m1:p2,n) = da(m1:p2,n) * dxfact(m1:p2)
          enddo
        endif

#ifdef NDSLPPM        
!        call vertical_cell_ppm_intp_ind(xnext(m1),da(m1,1),        &
        call rvertical_cell_ppm_intp(xnext(m1),da(m1,1),        &
     &   interfi(m1),next(m1,1),im+4,nvars)
#else
        call rfixend_cell_plm_intp_ind(xnext(m1),da(m1,1),         &
     &   interfi(m1),next(m1,1),im+4,nvars)	
#endif

        do n=1,nvars
          qq(1:im,k,n) = next(1:im,n)
        enddo

        enddo   ! for nst

      enddo ! for k

      deallocate(ue,qe,de)
      deallocate(var,past,da,next,    &
     &  	dxfact,xpast,xnext,         &
     &    dist,interfi)
!yj      print*,'finish rsopenend_cell_massadvx'
      return
      end subroutine rsopenend_cell_massadvx_ind

!-------------------------------------------------------------------
      subroutine rsopenend_cell_massadvy_ind(jm,lev,nvars,delt,vc,qq,mass)
!yj: This is semi-open end cell mass advy, the boundary is exterpolated
!yj: to outer 3 points and then do the fixend cell adv.
!yj      subroutine fixend_cell_massadvy(jm,jmh,lev,nvars,delt,vc,qq,mass)
!yj:RSM input jm = jgrd1      
!yj:RSM input lev = partial lev = levrp      
!
! compute local positive advection with mass conserving
! qq will be advect by vc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
!
!yj      use gfs_dyn_physcons
!yj      use gfs_dyn_layout1
      use physcons, pi => con_pi
!
      implicit none
!yj#include <paramodel.h>
!yj#include <npesi.h>
#include <rscomloc.h>
!
      integer   jm,lev,nvars,mass
!yj      integer   jm,jmh,lev,nvars,mass
      real	delt
      real	vc(jm,lev)  !!yj  vc(jgrd1,levpartial)
!!yj pass in vc(is:ie)-> larger domain for boundary (maybe skip this)
      real	qq(jm,lev,nvars)
! index for cell center extrapolation      
      integer m4,m3,m2,m1,m0
      integer p1,p2,p3,p4,p5
      logical loutflow(2)
! array after extrapolation
      real, allocatable :: ve(:),qe(:),de(:)
      real, allocatable :: past(:,:),da(:,:),next(:,:)
      real, allocatable :: var(:), dyfact(:)     &
     &      ,ypast(:),ynext(:),dist(:),interfi(:)
!
      real	step(10), dist_step
      integer  	n,k,j,nst,nstep
!
! preparations ---------------------------
!yj      print*,'enter rsopenend_cell_massadvy'
      m2=-2
      m1=-1
      m0=0
      p1=jm+1
      p2=jm+2
      p3=jm+3
      allocate(ve(m1:p2),qe(m1:p2),de(m1:p2))
      allocate(var(m1:p2),past(m1:p2,nvars),     &
     &    da(m1:p2,nvars),next(m1:p2,nvars),     &
     &  	dyfact(m1:p2),                         &
     &    ypast(m1:p2+1),ynext(m1:p2+1),         &
     &    dist (m1:p2+1),interfi(m1:p2+1))
      ve=0.0
      qe=0.0
      var=0.0
      past=0.0
      da=0.0
      next=0.0
      dyfact=0.0
      ypast=0.0
      ynext=0.0
      dist=0.0
      interfi=0.0
!
      de   = rdely  !if ds not equal, need to be modified
!   
!yj: c:cell center
!yj: i:cell interface
      do k=1,lev

! extrapolation for left and right domain 
        ve(1:jm) = vc(1:jm,k)
!yj        ve(jm:1:-1) = vc(1:jm,k)
!yj        call extp4_back(ve(1),de(m2),ve(m2))
!yj        call extp4_forw(ve(jm-3),de(jm-3),ve(jm-3))
!
!yj        do j=m2,p3
!yj          var(j)     = -ve(j) * delt  !yj at c
!yj        enddo
        do j=1,jm
          var(j)     = ve(j) * delt  !yj at c
        enddo
!yj:check interior region or exterior region >
! loutflow(1): exterpolation at index=1
! loutflow(2): exterpolation at index=jm
        loutflow = .false.
!        if (var(1)  .lt. 0.0)  loutflow(1)=.true.
!        if (var(jm) .gt. 0.0)  loutflow(2)=.true.
        if (var(1)  .lt. -0.5*de(1))   loutflow(1)=.true.
        if (var(jm) .gt.  0.5*de(jm))  loutflow(2)=.true.
!yj:check interior region or exterior region <
        if (loutflow(1)) then
          var(m0)=2*var(1)-var(2)               !outflow: negative
          var(m1)=2*var(m0)-var(1)               !outflow: negative
        else
          var(m0)=min(var(1),var(2))  !inflow: positive
          var(m1)=var(m0)           
        endif
        if (loutflow(2)) then
          var(p1)=2*var(jm)-var(jm-1)               !outflow: positive
          var(p2)=2*var(p1)-var(jm)               !outflow: positive
        else
          var(p1)=max(var(jm),var(jm-1)) !inflow: negative
          var(p2)=var(p1)
        endif

        do j=m0,p2
          !yj: (un)equal grid, moving distance of interface
          dist(j)=(var(j)*de(j-1)+var(j-1)*de(j)) / (de(j-1)+de(j))
        enddo
        dist(m1)=0.0
        dist(p2+1)=0.0
!
        call rdef_cfl_step(jm+5,dist(m1),de(m1),step,nstep)

        do nst = 1, nstep !yj: loop smaller distance

          interfi(m1) = 0.0
          ypast(m1)   = interfi(m1)
          ynext(m1)   = interfi(m1)
        do j=m0,p2+1
          interfi(j) = interfi(j-1)+de(j-1)
          dist_step  = dist(j)*step(nst) !yj: get smaller distance
          ypast(j)   = interfi(j) - dist_step !yj:gslati(latinterface)->j-1
          ynext(j)   = interfi(j) + dist_step
        enddo
#ifdef NDSLPPM        
!yj for ppm        
        interfi(:)=-interfi(:)
        ypast(:)=-ypast(:)
        ynext(:)=-ynext(:)
#endif
        if( mass.eq.1 ) then
         do j=m1,p2
          dyfact(j) = (ypast(j+1)-ypast(j)) / (ynext(j+1)-ynext(j))
         enddo
        endif
        do n=1,nvars
          qe(1:jm) = qq(1:jm,k,n)
          qe(m0)=0.0
          qe(m1)=0.0
          qe(p1)=0.0
          qe(p2)=0.0
          past(m1:p2,n) = qe(m1:p2) !yj at c
        enddo
#ifdef NDSLPPM        
!        call vertical_cell_ppm_intp_ind(interfi(m1),past(m1,1),    &
        call rvertical_cell_ppm_intp(interfi(m1),past(m1,1),    &
     &   ypast(m1),da(m1,1),jm+4,nvars)
#else
        call rfixend_cell_plm_intp_ind(interfi(m1),past(m1,1),          &
     &                       ypast(m1),da(m1,1),jm+4,nvars)
#endif
        if( mass.eq.1 ) then
          do n=1,nvars
            da(m1:p2,n) = da(m1:p2,n) * dyfact(m1:p2)
          enddo
        endif
#ifdef NDSLPPM        
!        call vertical_cell_ppm_intp_ind(ynext(m1),da(m1,1),        &
        call rvertical_cell_ppm_intp(ynext(m1),da(m1,1),        &
     &   interfi(m1),next(m1,1),jm+4,nvars)
#else
        call rfixend_cell_plm_intp_ind(ynext(m1),da(m1,1),              &
     &   interfi(m1),next(m1,1),jm+4,nvars)	
#endif
        do n=1,nvars
          qq(1:jm,k,n) = next(1:jm,n)
!yj          qq(1:jm,k,n) = next(jm:1:-1,n)
        enddo

        enddo   ! for nst

      enddo ! for k

      deallocate(ve,qe,de)
      deallocate(var,past,da,next,    &
     &  	dyfact,ypast,ynext,         &
     &    dist,interfi)

!yj      print*,'finish rsopenend_cell_massadvy'
      return
      end subroutine rsopenend_cell_massadvy_ind

!
!-------------------------------------------------------------------------------
      subroutine rfixend_cell_plm_intp_ind(pp,qq,pn,qn,lons,nv)
!
! linear quick version of the same two end points cell interpolation
! so pp(1) = pn(1) and pp(lonp+1) = pn(lonn+1)
!
! pp    location at interfac point as input
! qq    quantity at averaged-cell as input
! pn    location at interface of new grid structure as input
! qn    quantity at averaged-cell as output
! lons  numer of cells for dimension
!
! author : henry.juang@noaa.gov
!
!yj      use gfs_dyn_physcons
!yj      use gfs_dyn_layout1
!
      implicit none
!
      real      pp(lons+1)
      real      qq(lons  ,nv)
      real      pn(lons+1)
      real      qn(lons  ,nv)
      integer   lons,nv
!
      real      hfdp(lons),rdsi(lons)
      real      locf(lons+lons),df(lons+lons-1),dt(lons+lons-1)
      real      qc(lons),dq(lons),qmi(lons),qpi(lons)
      real      valf(lons+lons-1)
      real      rdd(lons)
      real      dd,check,ss
      integer   js(lons+1),i4j(lons+lons)
      integer   i,j,n,ip,in
!
! arrange input array cover output location with cyclic boundary
! condition
!
      locf(1) = pp(1)
      i4j(1)=1
      js(1)=1
      ip=2
      in=2
      do j=2,lons+lons-1 
        if( pp(ip).le.pn(in) ) then
          locf(j) = pp(ip)
          ip = ip + 1
        else
          locf(j) = pn(in)
          js(in) = j
          in = in + 1
        endif
        i4j(j) = ip - 1
      enddo
      locf(lons+lons) = pp(lons+1)
      js(lons+1) = lons+lons
!
! interpolation coefficient
!
      do i=1,lons
        hfdp(i) = 0.5*(pp(i+1)-pp(i))
      enddo
      do i=2,lons
        rdsi(i) = 1.0 / (hfdp(i)+hfdp(i-1))
      enddo
      do i=1,lons
        dd = 0.0
        do j=js(i),js(i+1)-1
          df(j) = locf(j+1) - locf(j)
          dd    = dd    + df(j)
        enddo
        rdd(i) = 1.0 / dd
      enddo
      do j=1,lons+lons-1
        i = i4j(j)
        dt(j) = (locf(j) + 0.5*df(j)) - (pp(i) + hfdp(i))
      enddo
!
! start interpolation by integral of ppm
!
      dq(1) = 0.0
      dq(lons) = 0.0

      do n=1,nv

        do i=1,lons
          qc(i) = qq(i,n)
        enddo
        do i=2,lons
          qmi(i) = (qc(i)-qc(i-1))*rdsi(i)
        enddo
        do i=2,lons-1
          qpi(i) = qmi(i+1)
        enddo
        do i=2,lons-1
          check = qmi(i)*qpi(i)
          if( check.lt.0.0 ) then
            dq(i) = 0.0
          else
            dq(i)  = 0.5*(qmi(i)+qpi(i))
          endif
        enddo

        do j=1,lons+lons-1
          i = i4j(j)
!          valf(j) = qc(i) + dq(i)*dt(j)
          if (qc(i).eq.0.0) then
            valf(j)=0.0
          else
            valf(j)=1.0
          endif
        enddo

        do i=1,lons
          ss = 1.0
          do j=js(i),js(i+1)-1
!            ss = ss + valf(j)*df(j)
            if (valf(j).eq.0.0) then
              ss=0.0
            endif
          enddo
!          qn(i,n) = ss * rdd(i)
          qn(i,n) = ss 
        enddo

      enddo
!
      return
      end subroutine rfixend_cell_plm_intp_ind
! ------------------------------------------------------------------------
      subroutine vertical_cell_ppm_intp_ind(pp,qq,pn,qn,levs,nvars)
!
! mass conservation in vertical interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in vertical with mass conservation is under the condition that
! pp(1)=pn(1), pp(levs+1)=pn(levs+1)
!
! pp    pressure at interfac level as input
! qq    quantity at layer as input
! pn    pressure at interface of new grid structure as input
! qn    quantity at layer as output
! levs  numer of verical layers
!
! author : henry.juang@noaa.gov
!
      implicit none
!
      real      pp(levs+1)
      real      qq(levs,nvars)
      real      pn(levs+1)
      real      qn(levs,nvars)
      integer   levs,nvars
!
      real      massm,massc,massp,massbot,masstop
      real      qmi(levs,nvars),qpi(levs,nvars)
      real      dql(nvars),dqh(nvars)
      real      hh(levs)
      real      dqi,dqimax,dqimin,dqmono(levs,nvars)
      real      tl,tl2,tl3,tlp,tlm,tlc
      real      th,th2,th3,thp,thm,thc
      real      dpp,dqq,c1,c2
      integer   i,k, kl, kh, kk, kkl, kkh,n
      integer, parameter :: mono=1
!
      if( pp(1).ne.pn(1) .or. pp(levs+1).ne.pn(levs+1) ) then
        print *,' Error in vertical_cell_ppm_intp_ind for domain values '
        print *,' i pp1 pn1 ppt pnt ',i,     &
     &          pp(1),pn(1),pp(levs+1),pn(levs+1)
        call abort
      endif
!
! prepare thickness for grid
!
      do k=1,levs
        hh(k) = pp(k+1)-pp(k)
      enddo
!
! prepare location with monotonic concerns
!
      do n=1,nvars
        massbot=(3.*hh(1)+hh(2))*qq(1,n)-2.*hh(1)*qq(2,n)
        massm = massbot/(hh(1)+hh(2))
        massc = qq(1  ,n)
        massp = qq(1+1,n)
        dqi = 0.25*(massp-massm)
        dqimax = max(massm,massc,massp) - massc
        dqimin = massc - min(massm,massc,massp)
        dqmono(1,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
        do k=2,levs-1
          massp = qq(k+1,n)
          massc = qq(k  ,n)
          massm = qq(k-1,n)
          dqi = 0.25*(massp-massm)
          dqimax = max(massm,massc,massp) - massc
          dqimin = massc - min(massm,massc,massp)
          dqmono(k,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
        enddo
        masstop=(3.*hh(levs)+hh(levs-1))*qq(levs,n)  &
     &             -2.*hh(levs)*qq(levs-1,n)
        massp = masstop/(hh(levs)+hh(levs-1))
        massc = qq(levs  ,n)
        massm = qq(levs-1,n)
        dqi = 0.25*(massp-massm)
        dqimax = max(massm,massc,massp) - massc
        dqimin = massc - min(massm,massc,massp)
        dqmono(levs,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
!
! compute value at interface with momotone
!
        do k=2,levs
          qmi(k,n)=(qq(k-1,n)*hh(k)+qq(k,n)*hh(k-1))/(hh(k)+hh(k-1)) &
     &       +(dqmono(k-1,n)-dqmono(k,n))/3.0
        enddo
        do k=1,levs-1
          qpi(k,n)=qmi(k+1,n)
        enddo
        qmi(1,n)=qq(1,n)
        qpi(1,n)=qq(1,n)
        qmi(levs,n)=qq(levs,n)
        qpi(levs,n)=qq(levs,n)
      enddo
!
! do monotonicity
!
      if( mono.eq.1 ) then
        do n=1,nvars
        do k=1,levs
          c1=qpi(k,n)-qq(k,n)
          c2=qq(k,n)-qmi(k,n)
          if( c1*c2.le.0.0 ) then
            qmi(k,n)=qq(k,n)
            qpi(k,n)=qq(k,n)
          endif
        enddo
        do k=1,levs
          c1=(qpi(k,n)-qmi(k,n))*(qq(k,n)-0.5*(qpi(k,n)+qmi(k,n)))
          c2=(qpi(k,n)-qmi(k,n))*(qpi(k,n)-qmi(k,n))/6.
          if( c1.gt.c2 ) then
            qmi(k,n)=3.*qq(k,n)-2.*qpi(k,n)
          else if( c1.lt.-c2 ) then
            qpi(k,n)=3.*qq(k,n)-2.*qmi(k,n)
          endif
        enddo
        enddo
      endif
!
! start interpolation by integral of ppm spline
!
      kkl = 1
      tl=0
      do n=1,nvars
        dql(n)=0.0
      enddo

      do k=1,levs

        kl = k
        kh = k + 1
! find kkh
        do kk=kkl+1,levs+1
          if( pn(kh).ge.pp(kk) ) then
            kkh = kk-1
            go to 100
          endif
        enddo
        print *,' Error in vertical_cell_ppm_intp_ind for no lev found '
        print *,' i kh kl ',i,kh,kl
        print *,' pn ',(pn(kk),kk=1,levs+1)
        print *,' pp ',(pp(kk),kk=1,levs+1)
        call abort
 100    continue
        th=(pn(kh)-pp(kkh))/hh(kkh)
        th2=th*th
        th3=th2*th
        thp = th3-th2
        thm = th3-2.*th2+th
        thc = -2.*th3+3.*th2
        do n=1,nvars
          dqh(n)=thp*qpi(kkh,n)+thm*qmi(kkh,n)+thc*qq(kkh,n)
        enddo
! mass interpolate
        if( kkh.eq.kkl ) then
          do n=1,nvars
            qn(k,n) = (dqh(n)-dql(n))/(th-tl)
          enddo
        else if( kkh.gt.kkl ) then
          dpp  = (1.-tl)*hh(kkl) + th*hh(kkh)
          do kk=kkl+1,kkh-1
            dpp = dpp + hh(kk)
          enddo
          do n=1,nvars
            dql(n) = qq(kkl,n)-dql(n)
            dqq  = dql(n)*hh(kkl) + dqh(n)*hh(kkh)
            do kk=kkl+1,kkh-1
              dqq = dqq + qq(kk,n)*hh(kk)
            enddo
            qn(k,n) = dqq / dpp
          enddo
        else
          print *,' Error in vertical_cell_ppm_intp_ind for lev messed up '
          print *,' i kh kl ',i,kh,kl
          print *,' pn ',(pn(kk),kk=1,levs+1)
          print *,' pp ',(pp(kk),kk=1,levs+1)
          call abort
        endif
! next one
        kkl = kkh
        tl  = th
        do n=1,nvars
          dql(n) = dqh(n)
        enddo

      enddo	! end of k loop
!
      return
      end subroutine vertical_cell_ppm_intp_ind
!!
! ------------------------------------------------------------------------
