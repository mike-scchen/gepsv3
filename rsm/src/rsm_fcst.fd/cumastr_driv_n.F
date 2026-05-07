      subroutine cumastr_driv_n  &
                   (nxj,nx,lev,dt,phi,pvom,pvol,ptte,pqte,ut,vt,tt,qt,  &
                    rcup,sd,evap,kcbot,kctop,ncld,plt,pt,lndj,hfx,kcnv,&
                    paphp1,dx,qlt,hpbl)

!c#######################################################################
!c                     subroutine description
!c
!c 1. parameter specification
!c
!c     nx : dimension of horizontal-direction
!c    lev : dimension of z-direction
!c     dt : time step for ut, vt, tt, qt updates                (s)
!c      g : gravity                                             (m/s2)
!c      r : air gas constant                                    (j/kg)
!c     cp : air specific heat constant                          (j/kg)
!c    phi : odd level geopotential (nx,lev)                     (m2/s2)
!c      u : x-direction velocity   (nx,lev) at current time lvl (m/s)
!c      v : y-direction velocity   (nx,lev) at current time lvl (m/s)
!c      t : potent temp on odd lvl (nx,lev) at current time lvl (k)
!c      q : moisture    on odd lvl (nx,lev) at current time lvl (kg/kg)
!c     ut : x-direction velocity   (nx,lev) at past/future  lvl (m/s)
!c     vt : y-direction velocity   (nx,lev) at past/future  lvl (m/s)
!c     tt : real temp.  on odd lvl (nx,lev) at past/future  lvl (k)
!c     qt : moisture    on odd lvl (nx,lev) at past/future  lvl (kg/kg)
!c     qlt: liquit water on odd lvl at past/futrue lvl (kg/kg)
!c     sd : vertical velocity (nx,lev) mb/s
!c   evap: evaporation      (nx)               (w/m2)
!c#####################################################################

      use mo_constants, only:vtmpc1, inicon
      USE mo_cumulus_flux, only : cuparam
      use physcons, g=>con_g, cp =>con_cp , r => con_rd

      implicit none
!c input & output variable
      integer nx,nxj,lev,ncld,j,jj
      integer kcbot(nx),kctop(nx)
      real phi(nx,lev)                                     &
!         , u(nx,lev),v(nx,lev),t(nx,lev),q(nx,lev*ncld)             &
          , ut(nx,lev),vt(nx,lev),tt(nx,lev),qt(nx,lev)         &
          , evap(nx),sd(nx,lev)                                     &
          , pk(nx,lev),pk2(nx,lev)                                   &
          , plt(nx,lev),pt(nx),qlt(nx,lev,2)
!c  local work arrays
      real pkxmb(nx,lev)
      real pk2x(nx,lev),pkx(nx,lev)
!ccumastr
      real papp1(nx,lev),paphp1(nx,lev+1),pgeo(nx,lev)              &
           ,zqsat(nx,lev),pverv(nx,lev)                              &
           ,ztp1(nx,lev),zqp1(nx,lev),zup1(nx,lev),zvp1(nx,lev)      &
           ,ptte(nx,lev),pqte(nx,lev),pvom(nx,lev),pvol(nx,lev)      &
           ,zxp1(nx,lev),pxtec(nx,lev),ptop,dth,dt       &
           ,ztu(nx,lev),zqu(nx,lev),zmfu(nx,lev),zmfd(nx,lev)        &
           ,pcte(nx,lev)
      real pqhfl(nx),prsfc(nx),pssfc(nx),rcup(nx),zrain(nx)        &
            ,phhfl(nx),hfx(nx)
!
      real rhoh2o,dx,d2r,xlatj,xlat,tt1,rho
!c
      integer klevp1,klevm1,k,i,kc
      integer lndj(nx)
!xb110>
      real zew,zqs,zcor,foeewm,mdlon,tem1
      real hpbl(nx)
      logical locum(nx)
      integer kcnv(nx)
      real,parameter :: tf=233.16,tcr=263.16, tcrf=1.0/(tcr-tf)
!xb110<

      call cuparam
      call inicon 

      rhoh2o=1000.
!
      do k=1,lev
      do i=1,nxj
        papp1(i,k)=plt(i,k)*1000.  ! from cb to pa
        paphp1(i,k)= paphp1(i,k)*1000. ! from cb to pa
        pgeo(i,k) = phi(i,k) 
        pverv(i,k)=sd(i,k)*1000.  ! from cb to pa
      enddo
      enddo
      do i=1,nxj
        paphp1(i,lev+1)= paphp1(i,lev+1)*1000. !from cb to pa
      enddo

      do  i=1,nxj
      zrain(i) = 0.0
      locum(i) =.false.

      rho = paphp1(i,lev+1)/&
     &         (r*(tt(i,lev)*(1.+vtmpc1*qt(i,lev))))

      phhfl(i) = hfx(i)*rho*cp !xb110, sensible heat flux, convert unit from m/s to w/m2
      pqhfl(i) = evap(i)*rho   !xb110, upward surface moisture flux (kgm^-2s^-1)
      enddo

      dth = dt

      do k=1,lev
      kc=k
      do i=1,nxj
        ztp1(i,kc)=tt(i,k)
        zqp1(i,kc)=qt(i,k)
        zup1(i,kc)=ut(i,k)
        zvp1(i,kc)=vt(i,k)
        pcte(i,kc)=0.0
        tt1=tt(i,k)

        zew  = foeewm(tt1)
        zqs  = zew/papp1(i,kc)
        zqs  = min(0.5,zqs)
        zcor = 1./(1.-vtmpc1*zqs)
        zqsat(i,kc)=zqs*zcor

      enddo
      enddo
!============================================================
      klevp1=lev+1
      klevm1=lev-1

      call cumastr_n  &
             (nxj,nx,lev,klevp1,klevm1,ztp1,zqp1,zup1,zvp1,pverv   &
             ,zqsat,pqhfl,papp1,paphp1,pgeo,ptte,pqte,pxtec,pvom,pvol   &
             ,prsfc,pssfc,kcbot,kctop,dth,ztu,zqu,zmfu,zmfd,zrain,pcte  &
             ,phhfl,lndj,locum,kcnv,dx,hpbl)

      do k=1,lev
      kc=k
      do i=1,nxj
        tt(i,k)=tt(i,k)+ptte(i,kc)*dt
        qt(i,k)=qt(i,k)+pqte(i,kc)*dt
        ut(i,k)=ut(i,k)+pvom(i,kc)*dt
        vt(i,k)=vt(i,k)+pvol(i,kc)*dt
!xb110> output cloud ice for wsm6 MP
        if (qlt(i,k,2)>-999.0)then
          tem1 =max(0.0, min(1.0, (tcr-tt(i,k))*tcrf))
          qlt(i,k,1) = max(qlt(i,k,1) + pxtec(i,kc)*dt,0.) * tem1 !cloud ice
          qlt(i,k,2) = max(qlt(i,k,2) + pxtec(i,kc)*dt,0.) * (1.0-tem1) !cloud water
        else
          qlt(i,k,1)=max(qlt(i,k,1)+pxtec(i,kc)*dt,0.) !cloud water
        end if
!xb110<
      enddo
      enddo

      do i=1,nxj
        rcup(i)=(prsfc(i)+pssfc(i))*dt/rhoh2o
      enddo

      return
      end
