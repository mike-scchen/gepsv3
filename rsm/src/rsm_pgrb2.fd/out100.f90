      subroutine out100(nx,my,lev,t2,q2,u10,v10,ps,ts,                  &
     &                  pt,ut,vt,tt,qt,fx100m)
!    Q2 T2 U10 V10 TSFC (igrd1*jgrd1)    => q2 t2 u10 v10 ts
!    FXS(igrd1,ksps,jgrd1): surface pressure  => ps
!    FXS(igrd1,kspn,jgrd1): pressure(on sig) [kPa]  => pt
!    FXS(igrd1,kstn,jgrd1): temperature(on sig) => tt
!    FXS(igrd1,ksu,jgrd1): U wind(on sig)  => ut
!    FXS(igrd1,ksv,jgrd1): V wind(on sig)  => vt
!    FXS(igrd1,ksq,jgrd1): specific humidity(on sig) => qt
!    fx100m(igrd1,5,jgrd1): output P T Q U V at 100m above terrain surface
!
      implicit  none

      integer   nx,my,lev,nxmy
      real, parameter ::rgas=2.8705E+2,grav=9.8

      real     u10(nx*my),v10(nx*my),t2(nx*my),q2(nx*my),ps(nx*my),     &
     &         ts(nx*my),pt(nx,lev,my),ut(nx,lev,my),vt(nx,lev,my),     &
     &         tt(nx,lev,my),qt(nx,lev,my),fx100m(nx,5,my)
      integer  i,j,n,ll,mm
!
      integer,parameter :: m= 1
      real psp(nx*my),ptp(nx,lev,my)
      real   avett,pp,p2,p10,hm(m),temp
      data hm/100.0/
!
      nxmy=nx*my
!----------------------------------------------------------------------
      do mm=1,m    ! 1: 100m.
!-----------------------------------------------------------------------
! calculate PQUVT at new layer.
      n=0
!      print*,'t2:',t2
!      print*,'ps:',ps
      do j=1,my
       do i=1,nx
       n=n+1
      psp(n)=ps(n)*1000.
        do ll=1,lev
        ptp(i,ll,j)=pt(i,ll,j)*1000.
        enddo
       enddo
      enddo
      n=0
      do j=1,my
       do i=1,nx
       n=n+1
! ave. temp. between surface and 1st sig-P layer
          temp=((psp(n)-ptp(i,1,j))/2.)/(ptp(i,1,j)-ptp(i,2,j))+1.0
          avett=(ts(n)+tt(i,1,j))/2.0 *                                 &
     &          (1.0+(0.608*(qt(i,2,j)+temp*(qt(i,1,j)-qt(i,2,j)))))
! calculate 40m&100m P, 2m P, 10m P.
        pp=psp(n)*(exp(-(hm(mm)*grav)/(rgas*avett)))
        p2=psp(n)*(exp(-(2.0*grav)/(rgas*avett)))
        p10=psp(n)*(exp(-(10.0*grav)/(rgas*avett)))
        fx100m(i,1,j)=pp

! find pp pressure position
         if(pp.gt.ptp(i,1,j))then
          temp=(ptp(i,1,j)-pp)/(ptp(i,1,j)-p2)
          fx100m(i,2,j)=t2(n)*temp+tt(i,1,j)*(1.0-temp)
          fx100m(i,3,j)=max(q2(n)*temp+qt(i,1,j)*(1.0-temp),1.e-8)
          temp=(ptp(i,1,j)-pp)/(ptp(i,1,j)-p10)
          fx100m(i,4,j)=u10(n)*temp+ut(i,1,j)*(1.0-temp)
          fx100m(i,5,j)=v10(n)*temp+vt(i,1,j)*(1.0-temp)
          goto 99
         endif

         do ll=1,10
          if((pp.le.ptp(i,ll,j)).and.(pp.gt.ptp(i,ll+1,j)))then
           temp=(ptp(i,ll,j)-pp)/(ptp(i,ll,j)-ptp(i,ll+1,j))
           fx100m(i,2,j)=tt(i,ll+1,j)*temp+tt(i,ll,j)*(1.0-temp)
           fx100m(i,3,j)=                                               &
     &              max(qt(i,ll+1,j)*temp+qt(i,ll,j)*(1.0-temp),1.e-8)
           fx100m(i,4,j)=ut(i,ll+1,j)*temp+ut(i,ll,j)*(1.0-temp)
           fx100m(i,5,j)=vt(i,ll+1,j)*temp+vt(i,ll,j)*(1.0-temp)
           goto 99
          endif
         enddo

 99   continue

       enddo  ! end (i)
      enddo  ! end (j)

!-----------------------------------------------------------------------
      enddo  ! end (mm)
!=======================================================================
      return
      end
