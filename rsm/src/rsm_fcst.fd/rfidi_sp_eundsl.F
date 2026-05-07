#include <define.h>
      subroutine rfidi_sp_eundsl(ug,vg,tg,rqg,dlam,dphi,udx,udy,vdx,vdy, &    
     &                 tdx,tdy,          dlax,spdlat,                    &
     &                 xm2,xm2px,xm2py,corf,gzdx,gzdy,                   &
     &                 lonlens,igrd1s,levr,levrm1,levrp1,ntotal,         &
     &                 sigma,                                            &
     &                 uadv,vadv,tadv,qsadv,qadv,dot                     &
     &           )
!fpp$ noconcur r
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rfidi_sp_eundsl
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  do regional dynamic computation (full tendency)
!
! program history log:
!   add sigma-pressure coordinate with def p = sigma1 * psfc + sigma2
!       --- ying-ju chen (org: cwb; date: 2022.Jan.25)
!
! usage:    call rfidi_sp(ug,vg,tg,rqg,dlam,dphi,udx,udy,vdx,vdy,
!          1           tdx,tdy,          dlax,spdlat,
!          2           xm2,xm2px,xm2py,corf,gzdx,gzdy,
!          3           lonlens,igrd1s,levr,levrm1,levrp1,ntotal,
!          4           sigma,
!          5           uadv,vadv,tadv,qsadv,qadv,dot)
!
!   input argument list:
!     ug        - full grid field of u
!     vg        - full grid field of v
!     tg        - full grid field of virtual temperaturev
!     rqg       - full grid field of specific humidity
!     dlam      - x-derivative of log (surface pressure) 
!     dphi      - y-derivative of log (surface pressure) 
!     udx       - x-derivative of u
!     udy       - y-derivative of u
!     vdx       - x-derivative of v
!     vdy       - y-derivative of v
!     tdx       - x-derivative of virtual temperature
!     tdy       - y-derivative of virtual temperature
!     xm2       - map facter square
!     xm2px     - x-derivative of map factor square
!     xm2py     - y-derivative of map factor square
!     corf      - coriolis force
!     gzdx      - x-derivative of terrain height
!     gzdy      - y-derivative of terrain height
!     sigma     - interface values, sigma1 = Bk; sigma2 = Ak; P[mb]=Ak+Bk*Ps[mb]
!                 (bottom:1, top:levr+1)
!     uadv      - horizonal advection of u
!     vadv      - horizonal advection of v
!     tadv      - horizonal advection of t
!     qsadv     - horizonal advection of log(ps)
!     qadv      - horizonal advection of tracers
!     dot       - sigma-p vertical flux
!
!   output argument list:
!     (vertical advection hasn't done yet in this program)
!     ug        - full tendency grid field of u
!     vg        - full tendency grid field of v
!     tg        - full tendency grid field of virtual temperaturev
!     rqg       - full tendency grid field of specific humidity
!     dlax      - full tendency grid field of log surface pressure
!
!   local variables:
!     p1        - inverse of rpi
!     p2        - rpi
!     rdelp2    - inverse of pressure layer thickness
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
      use physcons, rocp => con_rocp, rd => con_rd
!.....
!.....developed by hann-ming henry juang (dec 1990, version 0)......
!..... this version of total forcing computation including
!..... linear and nonlinear in equations.
!.....
!
      dimension   sigma(levr+1,2),sia(levr+1),sib(levr+1),               &
     &   dsi(levr),sl2(levr),pi(igrd1s,levr+1),ps(igrd1s),               &
     &   dp(igrd1s,levr),rp2(igrd1s,levr),rdelp2(igrd1s,levr),           &
     &   pak(igrd1s)
      real pl
      dimension   spdlat(levr),                                          &
     &   tg(igrd1s,levr),tdx(igrd1s,levr),tdy(igrd1s,levr),              &
     &   ug(igrd1s,levr),udx(igrd1s,levr),udy(igrd1s,levr),              &
     &   vg(igrd1s,levr),vdx(igrd1s,levr),vdy(igrd1s,levr),              &
     &   rqg(igrd1s,levr,ntotal),                                        &
!    &   rqdx(igrd1s,levr,ntotal),rqdy(igrd1s,levr,ntotal),              &
     &   dphi(igrd1s),dlam(igrd1s)
      dimension  dlax(igrd1s),                                           &
     &    b(igrd1s,levr),  a(igrd1s,levr),                               &
     &  tau(igrd1s,levr),rtg(igrd1s,levr,ntotal)
      dimension vadv(igrd1s,levr),uadv(igrd1s,levr),                     &
     &  qsadv(igrd1s,levr),tadv(igrd1s,levr),                            &
     &  qadv(igrd1s,levr,ntotal)
      dimension  p1(levr),p2(levr),                                      &
     &                xm2(igrd1s),xm2px(igrd1s),                         &
     &  xm2py(igrd1s),corf(igrd1s),                                      &
     &  gzdx(igrd1s),gzdy(igrd1s)
!
      dimension dg(igrd1s,levr), ek(igrd1s,levr),spd(igrd1s,levr),       &
     &     cg (igrd1s,levr), db(igrd1s,levr),cb(igrd1s,levr),            &
     &     dot(igrd1s,levrp1),dup(igrd1s,levr),dvp(igrd1s,levr),         &
     &     dum(igrd1s,levr ),dvm(igrd1s,levr)
!    &     dum(igrd1s,levr ),dvm(igrd1s,levr),vvel(igrd1s,levr)
!
! --- initialize arrays ---
!
      sia=0.e0
      sib=0.e0
      dsi=0.e0
      sl2=0.e0
      ps=0.e0
      pi=0.e0
      dp=0.e0
      rp2=0.e0
      rdelp2=0.e0
      dot=0.e0
      dvm=0.e0
      dum=0.e0
      dvp=0.e0
      dup=0.e0
      a=0.e0
      b=0.e0
      tau=0.e0
      rtg=0.e0
!
! --- prepare coordinate related values ---
!
!     sigma(:,1)=bki, sigma(:,2)=aki, and p=aki+ps*bki
!     sigma(:,2)=aki in the unit mb
!     ps,pi in the unit cb
!     vertical index: 1 -> bottom;
!                     levr(levr+1) -> top layer(level)
!     use sia<-aki; sib<-bki
!
      do k=1,levr+1
          sia(k) =sigma(k,2)
          sib(k) =sigma(k,1)
      enddo
!
      do k=1,levr
        dsi(k)=sib(k)-sib(k+1)
        sl2(k)=sib(k)+sib(k+1)
      enddo
!
!     input dlax is surface pressure [unit: cb]
      do j=1,lonlens
        ps(j)=exp(dlax(j))
      enddo
      pi(:,levr+1)=0.0
      do k=1,levr
        do j=1,lonlens
          pi(j,k)=sib(k)*ps(j)*10.+sia(k)
          pi(j,k)=pi(j,k)*0.1
        enddo
      enddo
!
      do k=1,levr
        do j=1,lonlens
           dp(j,k)=pi(j,k)-pi(j,k+1)
          rp2(j,k)=1./(pi(j,k)+pi(j,k+1))
        enddo
      enddo
!
      do k=1,levr
        do j=1,lonlens
          rdelp2(j,k)=1./(2.*dp(j,k))
        enddo
      enddo
!
      do 20 k=1,levr
      do 20 j=1,lonlens
      dg(j,k)=udx(j,k)+vdy(j,k)
      ek(j,k)=(ug(j,k)*ug(j,k)+vg(j,k)*vg(j,k))*0.5e0
      spd(j,k)=sqrt(2.*ek(j,k))*xm2(j)
  20  continue
      do k=1,levr
      spdlat(k) = 0.0e0
      do j=1,lonlens
        spdlat(k) = max( spdlat(k),spd(j,k))
      enddo
      enddo
!
! --- compute pressure gradient force ---
!
! pressure gradient force (PGF) is
!     -(R T div(Qk) + div(phi_k)
!        = -R T Bki*Ps/Pk *grad(Qs) .....................(A)
!        + grad(phi_s) ..................................(B)
!        + grad(Sum_1,k-1 {R (Pi-Pi+1)/(Pi+Pi+1) Ti}) ...(C)
!        + grad(Sum_1,k {R (Pi-Pi+1)/(Pi+Pi+1) Ti}) .....(D)
!
! PGF: (B) is ready as -gzdx and -gzdy
!
! PGF: expand and calculate (C)+(D); store in dup and dvp
!
      dvm=0.e0
      dum=0.e0
      dvp=0.e0
      dup=0.e0
      do k=1,levr-1
        do j=1,lonlens
          pak(j)=dsi(k)-dp(j,k)*rp2(j,k)*sl2(k)

          dum(j,k)=dp(j,k)*tdx(j,k)+tg(j,k)*ps(j)*dlam(j)*pak(j)
          dum(j,k)=rd*rp2(j,k)*dum(j,k)
          dup(j,k  )=dup(j,k)+dum(j,k)
          dup(j,k+1)=dup(j,k)+dum(j,k)

          dvm(j,k)=dp(j,k)*tdy(j,k)+tg(j,k)*ps(j)*dphi(j)*pak(j)
          dvm(j,k)=rd*rp2(j,k)*dvm(j,k)
          dvp(j,k  )=dvp(j,k)+dvm(j,k)
          dvp(j,k+1)=dvp(j,k)+dvm(j,k)
        enddo
      enddo
      k=levr
      do j=1,lonlens
        pak(j)=dsi(k)-dp(j,k)*rp2(j,k)*sl2(k)
        dum(j,k)=dp(j,k)*tdx(j,k)+tg(j,k)*ps(j)*dlam(j)*pak(j)
        dum(j,k)=rd*rp2(j,k)*dum(j,k)
        dup(j,k)=dup(j,k)+dum(j,k)
        dvm(j,k)=dp(j,k)*tdy(j,k)+tg(j,k)*ps(j)*dphi(j)*pak(j)
        dvm(j,k)=rd*rp2(j,k)*dvm(j,k)
        dvp(j,k)=dvp(j,k)+dvm(j,k)
      enddo
!
! PGF: calculate (A), then store in dum and dvm
!
      do k=1,levr
        do j=1,lonlens
          pak(j  )=rd*tg(j,k)*ps(j)*sl2(k)*rp2(j,k)
          dum(j,k)=pak(j)*dlam(j)
          dvm(j,k)=pak(j)*dphi(j)
        enddo
      enddo
!
! a=dv/dt  b=du/dt
      do 100 k=1,levr
      do 100 j=1,lonlens
#ifdef NDSLDYN
       a(j,k)= vadv(j,k) -ek(j,k)*xm2py(j)
       b(j,k)= uadv(j,k) -ek(j,k)*xm2px(j)
#else
       a(j,k)= ug(j,k) * vdx(j,k) + vg(j,k) * vdy(j,k)
       b(j,k)= ug(j,k) * udx(j,k) + vg(j,k) * udy(j,k)
       a(j,k)=-a(j,k)*xm2(j)-ek(j,k)*xm2py(j)
       b(j,k)=-b(j,k)*xm2(j)-ek(j,k)*xm2px(j)
#endif
       a(j,k)=a(j,k)-gzdy(j)-dvp(j,k)-dvm(j,k)
       b(j,k)=b(j,k)-gzdx(j)-dup(j,k)-dum(j,k)
       a(j,k)=a(j,k)-corf(j)*ug(j,k)
       b(j,k)=b(j,k)+corf(j)*vg(j,k)
  100 continue

      do 120 k=1,levr
      do 120 j=1,lonlens
#ifdef NDSLDYN
      tau(j,k)=tadv(j,k)
#else
      tau(j,k)= ug(j,k) * tdx(j,k) + vg(j,k) * tdy(j,k)
      tau(j,k)=-tau(j,k)*xm2(j)
#endif
  120 continue
!
!     do 140 nc=1,ntotal
!     do 140 k=1,levr
!     do 140 j=1,lonlens
!     rtg(j,k,nc)= ug(j,k) * rqdx(j,k,nc) + vg(j,k) * rqdy(j,k,nc)
!     rtg(j,k,nc)=-rtg(j,k,nc)*xm2(j)
! 140 continue
!
! Next, want vertical flux zeta_dot(partial P/ partial zeta) at levels
! compute from continuity equation, and we have
! 1/ps zeta_dot(partial P/ partial zeta)
!   = -m^2 Sum_k,K{(Bi-Bi+1)V.grad(Qs)+(Pi-Pi+1)/Ps Div}
!   + Bk m^2 Sum_1,K{(Bi-Bi+1)V.grad(Qs)+(Pi-Pi+1)/Ps Div}
!
! cg <- V.grad(Qs)
!
      do 180 k=1,levr
      do 180 j=1,lonlens
#ifdef NDSLDYN
      cg(j,k)= -qsadv(j,k)
#else
      cg(j,k)=ug(j,k)*dlam(j)+ vg(j,k)*dphi(j)
#endif
  180 continue
!
      do 190 j=1,lonlens
      db(j,1)=dp(j,1)*dg(j,1)/ps(j)
      cb(j,1)=dsi(1)*cg(j,1)
  190 continue
      do 210 k=1,levrm1
      do 210 j=1,lonlens
      db(j,k+1)=db(j,k)+dp(j,k+1)*dg(j,k+1)/ps(j)
      cb(j,k+1)=cb(j,k)+dsi(k+1)*cg(j,k+1)
  210 continue
!
!   store forcing of lnp in dlax
!
      do 220 j=1,lonlens
#ifdef NDSLDYN
      dlax(j)= -cb(j,levr)-db(j,levr)*xm2(j)
#else
      dlax(j)= -(cb(j,levr)+db(j,levr))*xm2(j)
#endif
  220 continue
!
!   sigma dot computed only at interior interfaces.
!
      dvm=0.e0
      dum=0.e0
      dvp=0.e0
      dup=0.e0
      do 230 j=1,lonlens
      dot(j,1)=0.e0
      dvm(j,1)=0.e0
      dum(j,1)=0.e0
      dot(j,levrp1)=0.e0
      dvp(j,levr )=0.e0
      dup(j,levr )=0.e0
  230 continue
!
! dot is negative vertical flux
!     - zeta_dot (partial P)/(partial zeta) at levels
!     do 240 k=2,levr
!     do 240 j=1,lonlens
!     dot(j,k)=-ps(j)*xm2(j)*(                                           &
!    &                 (sib(k)-1) * (db(j,levr)+cb(j,levr))              &
!    &                 + db(j,k-1) + cb(j,k-1)                           &
!    &                      )
! 240 continue
      do 240 k=1,levrm1
      do 240 j=1,lonlens
      dot(j,k+1)=dot(j,k)+                                               &
     &        ( dsi(k)*(db(j,levr)+cb(j,levr)-cg(j,k))                   &
     &         -dp(j,k)*dg(j,k)/ps(j) )*xm2(j)

  240 continue
! change dot without psfc
      do k=2,levr
      do j=1,lonlens
        dot(j,k) = dot(j,k)*ps(j)
      enddo
      enddo

!      call maxmin(dot,lonlens,levrm1,1,levrm1,'in rfidi,dot')
!
! pak is (1/p) (dp/dt) at layers
! vvel = (dp/dt) for physics
!
      rk= rocp
!     do 250 k=2,levr
!     do 250 j=1,lonlens
!       vvel(j,k)=0.5*xm2(j)*ps(j)*(                                     &
!    &                  sl2(k)*cg(j,k)                                   &
!    &                 -(db(j,levr)-db(j,k-1)+cb(j,levr)-cb(j,k-1))      &
!    &                 -(db(j,levr)-db(j,k  )+cb(j,levr)-cb(j,k  ))      &
!    &                                )
!       pak(j)=vvel(j,k)*2.0*rp2(j,k)
!       tau(j,k)=tau(j,k)+rk*tg(j,k)*pak(j)
! 250 continue
! for k=1, vvel and tau
!     do j=1,lonlens
!       vvel(j,1)=0.5*xm2(j)*ps(j)*(                                     &
!    &                  sl2(1)*cg(j,1)                                   &
!    &                 -(db(j,levr)          +cb(j,levr)          )      &
!    &                 -(db(j,levr)-db(j,1  )+cb(j,levr)-cb(j,1  ))      &
!    &                                )
!       pak(j)=vvel(j,1)*2.0*rp2(j,1)
!       tau(j,1)=tau(j,1)+rk*tg(j,1)*pak(j)
!     enddo
!
#ifdef NDSLDYN
      do k=1,levr
        do j=1,lonlens
          pak(j  )=rk*tg(j,k)*ps(j)*sl2(k)*rp2(j,k)
          pl=0.5*(pi(j,k)+pi(j,k+1))
          tadv(j,k)=tau(j,k)                                             &
     &     +pak(j)*xm2(j)*(cg(j,k)-cb(j,levr)-db(j,levr))                &
     &     -rk*tg(j,k)*0.5*(dot(j,k+1)+dot(j,k))/pl
        enddo
      enddo
#else
      do k=1,levr
        do j=1,lonlens
          pak(j  )=rk*tg(j,k)*ps(j)*sl2(k)*rp2(j,k)
          tau(j,k)=tau(j,k)                                              &
     &     +pak(j)*xm2(j)*(cg(j,k)-cb(j,levr)-db(j,levr))
        enddo
      enddo
#endif
!
#ifndef NDSLDYN
      dvm=0.e0
      dum=0.e0
      dvp=0.e0
      dup=0.e0
      do 260 k=1,levrm1
      do 260 j=1,lonlens
      dvp(j,k  )=vg(j,k+1)-vg(j,k)
      dup(j,k  )=ug(j,k+1)-ug(j,k)
      dvm(j,k+1)=vg(j,k+1)-vg(j,k)
      dum(j,k+1)=ug(j,k+1)-ug(j,k)
  260 continue
#endif
!ccccccccccccccccccccccccccccccccc

      do 270 k=1,levr
      do 270 j=1,lonlens
#ifdef NDSLDYN
      vadv(j,k)=a(j,k)      
#else
      vg(j,k)=a(j,k)                                                     &
     &           -( dot(j,k+1)*dvp(j,k)                                  &
     &             +dot(j,k  )*dvm(j,k) )*rdelp2(j,k)
      vadv(j,k)=vg(j,k)      
#endif
#ifdef NDSLDYN
      uadv(j,k)=b(j,k)     
#else
      ug(j,k)=b(j,k)                                                     &
     &           -( dot(j,k+1)*dup(j,k)                                  &
     &             +dot(j,k  )*dum(j,k) )*rdelp2(j,k)
      uadv(j,k)=ug(j,k)     
#endif
  270 continue

#ifndef NDSLDYN
      dum=0.e0
      dup=0.e0
      p1(1)    = 0.e0
      p2(levr) = 0.e0
      do 280 k=1,levrm1
      do 280 j=1,lonlens
!     dup(j,k  )=tg(j,k+1)-tg(j,k)
!     dum(j,k+1)=tg(j,k+1)-tg(j,k)
! if vertical flux of T is combined with vertical component of conti.equ.,
! then the previous calculated pak should not involved vertical part of 
! conti.equ.
! rpi(k) = ( pl(k+1) / pl(k) ) ** kapa
      p2(k)=( (pi(j,k)+pi(j,k+1))/(pi(j,k+1)+pi(j,k+2)) ) ** rk
      p1(k+1) = 1.e0 / p2(k)
      dup(j,k  )=p2(k  )*tg(j,k+1)-tg(j,k)
      dum(j,k+1)=tg(j,k+1)-p1(k+1)*tg(j,k)
  280 continue
!
      do 290 k=1,levr
      do 290 j=1,lonlens
      tg (j,k)=tau(j,k)                                                  &  
     &               -( dot(j,k+1)*dup(j,k)                              &
     &                 +dot(j,k  )*dum(j,k) )*rdelp2(j,k)
      tadv(j,k)=tg(j,k)
  290 continue
#endif
!
!      dum=0.e0
!      dup=0.e0
!      do nc=1,ntotal
!      do 300 k=1,levrm1
!      do 300 j=1,lonlens
!      dup(j,k  )=rqg(j,k+1,nc)-rqg(j,k,nc)
!      dum(j,k+1)=rqg(j,k+1,nc)-rqg(j,k,nc)
!  300 continue
!!     do 310 j=1,lonlens
!!     dup(j,levr)=0.e0
!! 310 continue
!      do 320 k=1,levr
!      do 320 j=1,lonlens
!      rqg(j,k,nc)=rtg(j,k,nc)                                            &  
!     &               -( dot(j,k+1)*dup(j,k)                              &
!     &                 +dot(j,k  )*dum(j,k) )*rdelp2(j,k)
!  320 continue
!      enddo
!
      return
      end
