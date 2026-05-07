#include <defopt.h>
      subroutine RSM_RECVA2Rdata(root_gfs,iclim,ldosfc,lsfcchg,  &
     &                     mype,master,nrorg,rcsln,rsnln,  &
     &                     lonlens,latlens)
!
! for initialize: ldosfc=t, lsfcchg=f
! for base field: ldosfc=f, lsfcchg=-(doesn't matter)
! for every 24hr: ldosfc=t, lsfcchg=t
! Note: lsfcchg only change land-sea-ice mask and tsea (sea-surface temperature)
!
      USE RSM_cc, ONLY: Coupler_id
      use module_gfsdata, only:                                        &
     &  kfh,knt,ktt,ksgeo,ksfcp,kuu,kvv,                               &
     &  kqt,kcw,krw,kiw,ksw,kgl,koz,                                   &
     &  ktg,ksmc,ksnr,kstc,ksimk,kslmk,                                &
     &  ktrace,kend,                                                   &
     &  glev,kmsoil,gaki,gbki,fhour,                                   &
     &  sgeo,psfc,temp,u,v,rqgb,simk,                                  &
     &  allocate_gfsdata,deallocate_gfsdata
      use physcons, only : qcmin => con_qcmin
     ! RSM sig value
      use sigmap_cons, only : aki, bki
!
#include <paramodel.h>
! base field
#include <rscomltb.h>
#include <rscomgrd.h>
#include <rscomspc.h>
!
      integer root_gfs
      integer k,len2,lens,len3
      integer lonlens,latlens
      real :: fhourin(1)
!
      logical ldosfc,lsfcchg
      integer iclim
! model grid field      
      real,allocatable ::  gr(:)
      real hold
      dimension  rcsln(lngrdp),rsnln(lngrdp)   ! for rotate u v to new projection
! model terrain file      
      integer nrorg
      real*4,dimension(:),allocatable:: tmp4
      real,allocatable ::  rorg(:),rorgp(:,:),qo(:,:)
      real,allocatable ::  pgz(:),fgz(:)
      real,allocatable ::  tern(:,:)
! array for vertical interpolation
      real,allocatable ::  s1(:,:,:),s1l(:,:,:)
      real,allocatable ::  si(:,:,:),sl(:,:,:)
! base field vertical interpolation
      real,allocatable ::  gteltb(:,:),guultb(:,:),gvvltb(:,:),grqltb(:,:,:)
! local array for read, write file      
      real,allocatable ::  g2recv(:,:)
! local array for GFS data (surface field)
      dimension gslmk(cigrd1,cjgrd1),gsimk(cigrd1,cjgrd1)
      dimension  gsnr(cigrd1,cjgrd1)
      dimension cgsea(cigrd1,cjgrd1),cgland(cigrd1,cjgrd1)
      dimension rr     (igrd1p,jgrd1p)
      dimension rkeep  (igrd1p,jgrd1p)
! tracer index from global model
      integer   ngindex(7)
      integer   ngtotal, ngoz, ngcw, ngrw, ngiw, ngsw, nggl
!     data ngtotal,ngoz,ngcw,ngrw,ngiw,ngsw,nggl/2,3,2,3,4,5,6/
! local variables
      integer nto3,ntcw,ntrw,ntiw,ntsw,ntgl
      real rdummy

! tracer index of rsm
      nto3 = 1+ntoz
      ntcw = 1+ntoz+1
      if (ncld.eq.5 ) then !WSM6 or GFDLMP
        ntrw=1+ntcw
        ntiw=2+ntcw
        ntsw=3+ntcw
        ntgl=4+ntcw
      elseif (ncld.eq.4 ) then !WSM5
        ntrw=1+ntcw
        ntiw=2+ntcw
        ntsw=3+ntcw
      elseif (ncld.eq.2 ) then !WSM3
        ntrw=1+ntcw
      elseif (ncld.eq.1 ) then !Zhao-Carr
        ntcw=ntcw
      else
        stop "Error: ncld undefined!"
      endif

      len2=clngrd
      lens=len2*kmsoil
      len3=clngrd*glev

! allocate variables
      call allocate_gfsdata(igrd1p,jgrd1p,ntotal)

      allocate(g2recv(cigrd1,cjgrd1),gr(lngrd),stat=ierr)
      if (ierr/=0) stop "recv_data: allocate g2recv fail"

      if (Coupler_id.lt.0) return    !   <- standalone mode
      call RSM_ANNOUNCE('start recv basefield data=',2)

!for atmospheric fields        
! fhour
      if (mype.eq.master) then
        call CMP_RECV(fhourin,1)
        fhour=fhourin(1)
        call RSM_ANNOUNCE('recv basefield fhour =',2)
        print*,'recv fhour=',fhour
      endif
      call rmpbcastr(fhour,1)
! tracer index from gfs
      if (mype.eq.master) then
        call CMP_INTEGER_RECV(ngindex,7)
        print*,'recv ngindex=',ngindex
      endif
      call rmpbcasti(ngindex,7)
! list tracer index
      ngtotal=ngindex(1)
      ngoz=ngindex(2)
      ngcw=ngindex(3)
      ngrw=ngindex(4)
      ngiw=ngindex(5)
      ngsw=ngindex(6)
      nggl=ngindex(7)

! sgeo (surface terrain)
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        print*,'recv_sgeo 1'
        call maxmin(g2recv,len2,1,1,1,'recv_sgeo')
        ! model grid-full2d
        call cgtorg(g2recv,gr, 1)
      endif
        ! model grid-partial2d
        call rmpgf2p(gr,igrd1,jgrd1,sgeo,igrd1p,jgrd1p,1)
        ! base field grid
        call rmpbcastr(g2recv,cigrd1*cjgrd1)
        call cgtobg(g2recv,gzb,1)

! psfc [hPa] (surface pressure)
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        print*,'recv_psfc 1'
        call maxmin(g2recv,len2,1,1,1,'recv_psfc')
        ! model grid-full2d
        call cgtorg(g2recv,gr, 1)
      endif
        ! model grid-partial2d
        call rmpgf2p(gr,igrd1,jgrd1,psfc,igrd1p,jgrd1p,1)
        ! base field grid
        call rmpbcastr(g2recv,cigrd1*cjgrd1)
        call cgtobg(g2recv,qltb,1)
        ! change unit to log([cb]) for base field
        qltb(:)=log(qltb(:)/10.)

        if (glev.ne.levr) then
          if (mype.eq.master) then
            print*,'GFS lev=',glev
            print*,'RSM lev=',levr
            print*,'ERROR: GFS and RSM must have same vertical layers'
            print*,'This version does NOT support changing vertical layers'
          endif
          stop 88
        endif
        if (glev.ne.levr) then
          allocate(gteltb(lngrdb,glev),guultb(lngrdb,glev),     &
     &     gvvltb(lngrdb,glev),grqltb(lngrdb,glev,ntotal))
        endif

! temperature
      do k=1,glev
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        call maxmin(g2recv,len2,1,1,1,'recv_temp')
      endif
        if (mype.eq.master) then
          ! model grid-full2d
          call cgtorg(g2recv,gr, 1)
        endif
        ! model grid-partial2d
        call rmpgf2p(gr,igrd1,jgrd1,temp(:,:,k),igrd1p,jgrd1p,1)
        ! base field grid
        call rmpbcastr(g2recv,cigrd1*cjgrd1)
        if (glev.ne.levr) then
          call cgtobg(g2recv,gteltb(:,k),1)
        else
          call cgtobg(g2recv,teltb(:,k),1)
        endif
      enddo
! u
      do k=1,glev
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        call maxmin(g2recv,len2,1,1,1,'recv_u')
      endif
        if (mype.eq.master) then
          ! model grid-full2d
          call cgtorg(g2recv,gr, 1)
        endif
        ! model grid-partial2d
        call rmpgf2p(gr,igrd1,jgrd1,u(:,k),igrd1p,jgrd1p,1)
        ! base field grid
        call rmpbcastr(g2recv,cigrd1*cjgrd1)
        if (glev.ne.levr) then
          call cgtobg(g2recv,guultb(:,k),1)
        else
          call cgtobg(g2recv,uultb(:,k),1)
        endif
      enddo
!v
      do k=1,glev
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        call maxmin(g2recv,len2,1,1,1,'recv_v')
      endif
        if (mype.eq.master) then
          ! model grid-full2d
          call cgtorg(g2recv,gr, 1)
        endif
        ! model grid-partial2d
        call rmpgf2p(gr,igrd1,jgrd1,v(:,k),igrd1p,jgrd1p,1)
        ! base field grid
        call rmpbcastr(g2recv,cigrd1*cjgrd1)
        if (glev.ne.levr) then
          call cgtobg(g2recv,gvvltb(:,k),1)
        else
          call cgtobg(g2recv,vvltb(:,k),1)
        endif
      enddo

! uv mapping rotation
      if (.false.) then !!!! some bug in this mapping rotation, skip for test!
!     base field
      do k=1,glev
        do n=1,lngrdb
          hold = uultb(n,k)
          uultb(n,k) = hold * bcsln(n) - vvltb(n,k) * bsnln(n)
          vvltb(n,k) = hold * bsnln(n) + vvltb(n,k) * bcsln(n)
        enddo
      enddo
!     init field
      if (fhour.eq.0.0) then
        do k=1,levr
          do n=1,lngrdp
            hold = u(n,k)
            u(n,k) = hold * rcsln(n) - v(n,k) * rsnln(n)
            v(n,k) = hold * rsnln(n) + v(n,k) * rcsln(n)
          enddo
        enddo
      endif
      endif

! spfh: rqgb(1)
      do k=1,glev
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        call maxmin(g2recv,len2,1,1,1,'recv_spfh')
      endif
        if (mype.eq.master) then
          ! model grid-full2d
          call cgtorg(g2recv,gr, 1)
        endif
        ! model grid-partial2d
        call rmpgf2p(gr,igrd1,jgrd1,rqgb(:,:,k,1),igrd1p,jgrd1p,1)
        ! base field grid
        call rmpbcastr(g2recv,cigrd1*cjgrd1)
        if (glev.ne.levr) then
          call cgtobg(g2recv,grqltb(:,k,1),1)
        else
          call cgtobg(g2recv,rqltb(:,k,1),1)
        endif
      enddo
!
! convert temperature to virtual temperature
! |__ base field
      if (glev.ne.levr) then
        do k=1,glev
        do i=1,lngrdb
           gteltb(i,k)=(1.+0.61*grqltb(i,k,1))*gteltb(i,k)
        enddo
        enddo
      else
        do k=1,levr
        do i=1,lngrdb
           teltb(i,k)=(1.+0.61*rqltb(i,k,1))*teltb(i,k)
        enddo
        enddo
      endif
! |__ init field
        do k=1,glev
        do j=1,latlens
        do i=1,lonlens
           temp(i,j,k)=(1.+0.61*rqgb(i,j,k,1))*temp(i,j,k)
        enddo
        enddo
        enddo
!

! other tracers: rqgb(...)
      do ntrac=2,ngtotal
        if(ntrac.eq.ngcw)then
          ktrace=kcw                                ! cloud liquid water content
        else if(ntrac.eq.ngiw)then
          ktrace=kiw                                ! cloud ice content
        else if(ntrac.eq.ngrw)then
          ktrace=krw                                ! rain
        else if(ntrac.eq.ngsw)then
          ktrace=ksw                                ! snow 
        else if(ntrac.eq.nggl)then
          ktrace=kgl                                ! graupel
        else
          stop "ERROR: recv_data: wrong tracer index found !"
        endif

        do k=1,glev
        if (mype.eq.master) then
          call CMP_RECV(g2recv,len2)
          call maxmin(g2recv,len2,1,1,1,'recv_tracer')
        endif
          if (mype.eq.master) then
            ! model grid-full2d
            call cgtorg(g2recv,gr, 1)
          endif
          if(ntrac.eq.ngcw)then
            ! model grid-partial2d
            call rmpgf2p(gr,igrd1,jgrd1,rqgb(:,:,k,ntcw),igrd1p,jgrd1p,1)
            ! base field grid
            call rmpbcastr(g2recv,cigrd1*cjgrd1)
            if (glev.ne.levr) then
              call cgtobg(g2recv,grqltb(:,k,ntcw),1)
            else
              call cgtobg(g2recv,rqltb(:,k,ntcw),1)
            endif
          else if(ntrac.eq.ngiw)then
            ! model grid-partial2d
            call rmpgf2p(gr,igrd1,jgrd1,rqgb(:,:,k,ntiw),igrd1p,jgrd1p,1)
            ! base field grid
            call rmpbcastr(g2recv,cigrd1*cjgrd1)
            if (glev.ne.levr) then
              call cgtobg(g2recv,grqltb(:,k,ntiw),1)
            else
              call cgtobg(g2recv,rqltb(:,k,ntiw),1)
            endif
          else if(ntrac.eq.ngrw)then
            ! model grid-partial2d
            call rmpgf2p(gr,igrd1,jgrd1,rqgb(:,:,k,ntrw),igrd1p,jgrd1p,1)
            ! base field grid
            call rmpbcastr(g2recv,cigrd1*cjgrd1)
            if (glev.ne.levr) then
              call cgtobg(g2recv,grqltb(:,k,ntrw),1)
            else
              call cgtobg(g2recv,rqltb(:,k,ntrw),1)
            endif
          else if(ntrac.eq.ngsw)then
            ! model grid-partial2d
            call rmpgf2p(gr,igrd1,jgrd1,rqgb(:,:,k,ntsw),igrd1p,jgrd1p,1)
            ! base field grid
            call rmpbcastr(g2recv,cigrd1*cjgrd1)
            if (glev.ne.levr) then
              call cgtobg(g2recv,grqltb(:,k,ntsw),1)
            else
              call cgtobg(g2recv,rqltb(:,k,ntsw),1)
            endif
          else if(ntrac.eq.nggl)then
            ! model grid-partial2d
            call rmpgf2p(gr,igrd1,jgrd1,rqgb(:,:,k,ntgl),igrd1p,jgrd1p,1)
            ! base field grid
            call rmpbcastr(g2recv,cigrd1*cjgrd1)
            if (glev.ne.levr) then
              call cgtobg(g2recv,grqltb(:,k,ntgl),1)
            else
              call cgtobg(g2recv,rqltb(:,k,ntgl),1)
            endif
          endif
        enddo
      enddo
! ozone
      do k=1,glev
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        call maxmin(g2recv,len2,1,1,1,'recv_ozon')
      endif
        if (mype.eq.master) then
          ! model grid-full2d
          call cgtorg(g2recv,gr, 1)
        endif
        ! model grid-partial2d
        call rmpgf2p(gr,igrd1,jgrd1,rqgb(:,:,k,nto3),igrd1p,jgrd1p,1)
        ! base field grid
        call rmpbcastr(g2recv,cigrd1*cjgrd1)
        if (glev.ne.levr) then
          call cgtobg(g2recv,grqltb(:,k,nto3),1)
        else
          call cgtobg(g2recv,rqltb(:,k,nto3),1)
        endif
      enddo
        

      if (mype.eq.master) print*,'ok-recv part 1'
!-----------------------------------------------------------------------
! 2D surface variables (interplation is different)
! (only horizontal interpolation)
      if (ldosfc) then
! need to get the variables (1)~(4) first
! (1) land-sea mask of GFS data
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        print*,'recv_slmsk 1'
        call maxmin(g2recv,len2,1,1,1,'recv_slmsk')
        where(g2recv.ge.0.01)
          gslmk=1.0
        elsewhere
          gslmk=0.0
        end where
      endif
! (2) sheleg from GFS (for RSM land-sea-ice mask)
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        where (gsnr<1.e-3) gsnr=0.0
        print*,'recv_snr 1'
        call maxmin(gsnr,len2,1,1,1,'recv_snr')
      endif
! (3) sea and ice index from GFS (for RSM land-sea-ice mask)
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        print*,'recv_simk 1'
        call maxmin(gsimk,len2,1,1,1,'recv_simk')
      endif
!
      call rmpbcastr(gsnr,clngrd)
      call rmpbcastr(gsimk,clngrd)
      if (.not.lsfcchg) then
        call cg2rg_mpi(gsnr,sheleg,1)
      endif
      call cg2rg_mpi(gsimk,simk,1)
      do j=1,latlens
      do i=1,lonlens
        if( simk(i,j) .gt. 0.4 ) then
          simk(i,j) = 2.0
        else
          simk(i,j) = 0.0
        endif
      enddo
      enddo
! (4) add ice mask info to regional land-sea mask
      ! add sea ice info [0/1/2]
      do j=1,latlens
      do i=1,lonlens
        slmsk(i,j) = max(slmsk(i,j),simk(i,j))
      enddo
      enddo

! update climate var sea ice info
      if (iclim.eq.1) then
        do j=1,latlens
        do i=1,lonlens
          if (slmsk(i,j).eq.2.0) then
#ifdef doNOAHLSM
          !for NOAHLSM
            vtype(i,j)=15.0
            stype(i,j)=16.0
            zorl(i,j)=1.0
#else
          !for OSULSM
            vtype(i,j)=13.0
            stype(i,j)=9.0
            zorl(i,j)=1.0
#endif
          endif
        enddo
        enddo
      endif

! tsea (tg, ground temperature)
      if (lsfcchg) then
        rkeep=tsea
      endif
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        print*,'recv_tsea 1'
        call maxmin(g2recv,len2,1,1,1,'recv_tsea')

        call expnd(g2recv,gslmk,0,cgsea,cigrd1,cjgrd1)
        call expnd(g2recv,gslmk,1,cgland,cigrd1,cjgrd1)
      endif
      call rmpbcastr(cgsea,clngrd)
      call rmpbcastr(cgland,clngrd)
      call cg2rg_mpi(cgsea,tsea, 1)
      call cg2rg_mpi(cgland,rr,1)
      do j=1,latlens
      do i=1,lonlens
        if (slmsk(i,j).eq.1.0) tsea(i,j)=rr(i,j)
      enddo
      enddo
      do j=1,latlens
      do i=1,lonlens
        if( slmsk(i,j) .eq. 2.0e0 ) then
          tsea(i,j)= min( tsea(i,j), 271.20e0 )
        else if( slmsk(i,j).eq.0.0e0 ) then
          tsea(i,j)= max( tsea(i,j), 271.21e0 )
        else if( slmsk(i,j).eq.1.0e0 .and.                               &   
   &         sheleg(i,j).gt.0.0e0 ) then
          tsea(i,j)= min( tsea(i,j), 273.16e0 )
        endif
      enddo
      enddo
      call maxmin(  tsea,lngrdp,1,1,1,'tsea')
      if (lsfcchg) then
        ! replace only sea surface temperature
        do j=1,latlens
        do i=1,lonlens
          if( slmsk(i,j) .ne. 0.0 ) then
            tsea(i,j)= rkeep(i,j)
          endif
        enddo
        enddo
      endif

! smc
      do k=1,kmsoil
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        call maxmin(g2recv,len2,1,1,1,'recv_smc')
      endif
      if (.not.lsfcchg) then
        if (mype.eq.master) then
          call expnd(g2recv,gslmk,0,cgsea,cigrd1,cjgrd1)
          call expnd(g2recv,gslmk,1,cgland,cigrd1,cjgrd1)
        endif
        call rmpbcastr(cgsea,clngrd)
        call rmpbcastr(cgland,clngrd)
        call cg2rg_mpi(cgsea,smc(1,1,k), 1)
        call cg2rg_mpi(cgland,rr,1)
        do j=1,latlens
        do i=1,lonlens
          if (slmsk(i,j).eq.1.0) then
             smc(i,j,k)=rr(i,j)
          else
             smc(i,j,k)=0.0
          endif
        enddo
        enddo
        call maxmin(   smc(1,1,k),lngrdp,1,1,1,'smc')
      endif
      enddo

! stc
      do k=1,kmsoil
      if (mype.eq.master) then
        call CMP_RECV(g2recv,len2)
        call maxmin(g2recv,len2,1,1,1,'recv_stc')
      endif
      if (.not.lsfcchg) then
        if (mype.eq.master) then
          call expnd(g2recv,gslmk,0,cgsea,cigrd1,cjgrd1)
          call expnd(g2recv,gslmk,1,cgland,cigrd1,cjgrd1)
        endif
        call rmpbcastr(cgsea,clngrd)
        call rmpbcastr(cgland,clngrd)
        call cg2rg_mpi(cgsea,stc(1,1,k), 1)
        call cg2rg_mpi(cgland,rr,1)
        do j=1,latlens
        do i=1,lonlens
          if (slmsk(i,j).eq.1.0) then
             stc(i,j,k)=rr(i,j)
          else
             stc(i,j,k)=0.0
          endif
        enddo
        enddo
        call maxmin(   stc(1,1,k),lngrdp,1,1,1,'stc')
      endif
      enddo
! ... 
!
      if (mype.eq.master) print*,'ok-recv part 2'
      endif  ! ldosfc
!-----------------------------------------------------------------------
!
! vertical interpolation for base field if GFS and RSM have different layers
!
      if (glev.ne.levr) then
        if (mype.eq.master) then
          print*,'GFS lev=',glev
          print*,'RSM lev=',levr
          print*,'do vertical interpolation for base field'
        endif
! allocate var for vertical interpolation
      allocate(s1(lngrdb,1,glev+1),s1l(lngrdb,1,glev)) ! base sig value
      do k=1,glev+1
        do j=1,1
        do i=1,lngrdb
          rdummy=exp(qltb(i))*10.
          ! pressure at levels (gfs aki,bki)
          s1(i,j,k)=gaki(k)+gbki(k)*rdummy
          ! sigma value
          s1(i,j,k)=s1(i,j,k)/rdummy
        enddo
        enddo
      enddo
!     print*,'mype,psfc,s1=',mype,psfc(1,1),s1(1,1,:)
      do k=1,glev
        do j=1,1
        do i=1,lngrdb
          s1l(i,j,k)=(s1(i,j,k)+s1(i,j,k+1))/2.
        enddo
        enddo
      enddo
! get new sigma levels and layers based on new surface pressure
      allocate(si(lngrdb,1,levr+1),sl(lngrdb,1,levr)) !RSM sig value
      do k=1,levr+1
        do j=1,1
        do i=1,lngrdb
          rdummy=exp(qltb(i))*10.
          ! rsm aki,bki
          si(i,j,k)=(aki(k)+bki(k)*rdummy)/rdummy
        enddo
        enddo
      enddo
      do k=1,levr
        do j=1,1
        do i=1,lngrdb
          sl(i,j,k)=(si(i,j,k+1)+si(i,j,k))/2.
        enddo
        enddo
      enddo

      allocate(tern(lngrdb,levr))
      in=0
      tensn=10.0
      iuv=1
        call sg2sg(qltb,s1l(:,1,:),gteltb(1,1:glev),qltb,sl(:,1,:),tern,lngrdb,1,glev,levr,    &
     &             in,tensn,iuv,lngrdb,1)
        teltb(:,:)=tern(:,:)
        call sg2sg(qltb,s1l(:,1,:),guultb(1,1:glev),qltb,sl(:,1,:),tern,lngrdb,1,glev,levr,    &
     &             in,tensn,iuv,lngrdb,1)
        uultb(:,:)=tern(:,:)
        call sg2sg(qltb,s1l(:,1,:),gvvltb(1,1:glev),qltb,sl(:,1,:),tern,lngrdb,1,glev,levr,    &
     &             in,tensn,iuv,lngrdb,1)
        vvltb(:,:)=tern(:,:)
      do nc=1,ntotal
        call sg2sg(qltb,s1l(:,1,:),grqltb(1,1:glev,nc),qltb,sl(:,1,:),tern,lngrdb,1,glev,levr,    &
     &             in,tensn,iuv,lngrdb,1)
        rqltb(:,:,nc)=tern(:,:)
      enddo

      deallocate(tern)
      deallocate(si,sl)  ! base field sig value
      deallocate(s1,s1l) ! base field sig value
      !
      deallocate(gteltb,guultb,gvvltb,grqltb)
      !
      endif ! (glev.ne.levr)
!-----------------------------------------------------------------------
!
      if (fhour.eq.0.0) then

! fhour=0, get base field
      if (mype.eq.master) print*,'get base field at +0 hr'
      do n=1,lngrdb
        qb(n)= qltb(n)
      enddo
      do k=1,levr
      do n=1,lngrdb
        teb(n,k)= teltb(n,k)
        uub(n,k)= uultb(n,k)
        vvb(n,k)= vvltb(n,k)
      enddo
      enddo
      do nc=1,ntotal
      do k=1,levr
      do i=1,lngrdb
        rqb(i,k,nc)= rqltb(i,k,nc)
      enddo
      enddo
      enddo
! after horizontal interpolation, do vertical interpolation on model grid
!
! get regional terrain
! read(regional terrain)
        if( mype.eq.master ) then
          allocate(rorg(lngrd))
          allocate(tmp4(lngrd))
          rewind(nrorg)
          read(nrorg) tmp4
          rorg=tmp4
          deallocate(tmp4)
        endif
!>>> if terrain blending between GFS and RSM is wanted, add here <<<
        call rmpgp2f(sgeo,igrd1p,jgrd1p,gr,igrd1,jgrd1,1)

        if( mype.eq.master ) then
          allocate(pgz(lngrd))
          allocate(fgz(lnwav))
        percmtn=0.1
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
              rorg(i+jlats)=gr(i+jlats)*coes+rorg(i+jlats)*(1.-coes)
            enddo
          enddo
! modify gzo with regional spectral smooth
! init spectral constant
          call inpint_plnini
          call inpint_rftini
! get ptb
          do ij=1,lngrd
            pgz(ij)=rorg(ij)-gr(ij)
          enddo
          print *,'after pgz'
          call maxmin(pgz,lngrd,1,1,1,' before inpint_gdtocc, pgz ')
! grid to coef back to grid
          call inpint_gdtocc(pgz,fgz,1)
          call inpint_cctogd(fgz,pgz,1)
          print *,'after cc gd pgz'
          call maxmin(pgz,lngrd,1,1,1,' after inpint_cctogd, pgz ')
! get full
          do ij=1,lngrd
            rorg(ij)=pgz(ij)+gr(ij)
          enddo
          print *,'gz of lateral boundar modified'
          call maxmin(rorg,lngrd,1,1,1,'  rorg ')
        endif  !endif percmtn.gt.0.0
          deallocate(pgz)
          deallocate(fgz)
!
      endif  !endif mype.eq.master
!<<<
        allocate(rorgp(igrd1p,jgrd1p))
        call rmpgf2p(rorg,igrd1,jgrd1,rorgp,igrd1p,jgrd1p,1)
        if( mype.eq.master ) then
          deallocate(rorg)
        endif
!
! allocate var for vertical interpolation
      allocate(s1(igrd1p,jgrd1p,glev+1),s1l(igrd1p,jgrd1p,glev)) !GFS sig value
      do k=1,glev+1
        do j=1,latlens
        do i=1,lonlens
          ! pressure at levels (gfs aki,bki)
          s1(i,j,k)=gaki(k)+gbki(k)*psfc(i,j)
          ! sigma value
          s1(i,j,k)=s1(i,j,k)/psfc(i,j)
        enddo
        enddo
      enddo
!     print*,'mype,psfc,s1=',mype,psfc(1,1),s1(1,1,:)
      do k=1,glev
        do j=1,latlens
        do i=1,lonlens
          s1l(i,j,k)=(s1(i,j,k)+s1(i,j,k+1))/2.
        enddo
        enddo
      enddo
      ! change unit to log([cb])
      psfc(:,:)=log(psfc(:,:)/10.)

! get surface pressure over regional terrain
      if (mype.eq.master) call maxmin(rorgp,lngrdp,1,1,1,'  rorgp ')
      allocate(qo(igrd1p,jgrd1p))

      do j=1,latlens
        call newps(igrd1p,lonlens,glev,temp(1:igrd1p,j,1:glev),psfc(:,j),sgeo(:,j),  &
     &             rorgp(:,j),qo(:,j),s1(1:igrd1p,j,1:glev+1))
      enddo
      if (mype.eq.master) call maxmin(qo,lngrdp,1,1,1,'  qo ')

! get new sigma levels and layers based on new surface pressure
      allocate(si(igrd1p,jgrd1p,levr+1),sl(igrd1p,jgrd1p,levr)) !RSM sig value
      do k=1,levr+1
        do j=1,latlens
        do i=1,lonlens
          rdummy=exp(qo(i,j))*10.
          ! rsm aki,bki
          si(i,j,k)=(aki(k)+bki(k)*rdummy)/rdummy
        enddo
        enddo
      enddo
      do k=1,levr
        do j=1,latlens
        do i=1,lonlens
          sl(i,j,k)=(si(i,j,k+1)+si(i,j,k))/2.
        enddo
        enddo
      enddo
        
! vertical interpolation (output tern)     
      allocate(tern(lngrdp,levr))
      in=0
      tensn=10.0
      iuv=1
! te
! te vertical interpolation
      call sg2sg(psfc,s1l,temp,qo,sl,tern,igrd1p,jgrd1p,glev,levr,      &
     &           in,tensn,iuv,lonlens,latlens)
! te spectral coef.
      call getptb(tern,teb,wsltb,levr)
      call gdtocc(tern,te,levr)
! uu
! uu vertical interpolation
      call sg2sg(psfc,s1l,u,qo,sl,tern,igrd1p,jgrd1p,glev,levr,      &
     &           in,tensn,iuv,lonlens,latlens)
! uu spectral coef.
      call getptb(tern,uub,wsltb,levr)
      call gdtosc(tern,uu,levr)

! vv
! vv vertical interpolation
      call sg2sg(psfc,s1l,v,qo,sl,tern,igrd1p,jgrd1p,glev,levr,      &
     &           in,tensn,iuv,lonlens,latlens)
! vv spectral coef.
      call getptb(tern,vvb,wsltb,levr)
      call gdtocs(tern,vv,levr)
! rq
! rq vertical interpolation
! spfh, ozon, clwr, rain, qice, snow, grpl
!   1     2     3     4     5     6     7 
      do nc = 1,ntotal
      call sg2sg(psfc,s1l,rqgb(:,:,:,nc),qo,sl,tern,igrd1p,jgrd1p,glev,levr,      &
     &           in,tensn,iuv,lonlens,latlens)
      call getptb(tern,rqb(1,1,nc),wsltb,levr)
#ifdef NDSL
      rqgrd(:,:,nc)=tern(:,1:levr)
      ! keep the calculation of clwr,rain,qice,snow,grpl always at "full field"
      ! spfh and ozon still at perturbation
      if ( nc.gt.2 ) then
        call getful(tern,rqb(1,1,nc),wsltb,levr)
        rqgrd(:,:,nc) = tern(:,1:levr)
        rqb(:,:,nc) = qcmin
        ! make sure the full field of clwr,rain,qice,snow,grpl be positive
        do k=1,levr
          do i=1,lngrdp
            rqgrd(i,k,nc) = max(rqgrd(i,k,nc),qcmin)
          enddo
        enddo
      endif
#else
      call gdtocc(tern,rq(1,1,nc),levr)
#endif
      enddo

! get gz and q spectral coef.
! gz
      call getptb(rorgp,gzb,wsltb,1)
      call gdtocc(rorgp,gz,1)
! q
      call getptb(qo,qb,wsltb,1)
      call gdtocc(qo,q,1)
!
        do k=1,levr
          do i=1,lnwavp
            uum(i,k)=uu(i,k)
            vvm(i,k)=vv(i,k)
            tem(i,k)=te(i,k)
#ifdef NONHYD
            pnm(i,k)=pn(i,k)
            tnm(i,k)=tn(i,k)
            onm(i,k)=on(i,k)
#endif
          enddo
        enddo
#ifdef NONHYD
        k=levr+1
        do i=1,lnwavp
          onm(i,k)=on(i,k)
        enddo
#endif
        do nc=1,ntotal
          do k=1,levr
#ifdef NDSL
          do i=1,lngrdp
            rmgrd(i,k,nc)=rqgrd(i,k,nc)
          enddo
#else
          do i=1,lnwavp
            rm(i,k,nc)=rq(i,k,nc)
          enddo
#endif
          enddo
        enddo
        do i=1,lnwavp
          qm(i)=q(i)
        enddo
      endif   ! (fhour.eq.0.0)

      if (mype.eq.master) print*,'ok-vertical interpolation'
!-----------------------------------------------------------------------
! deallocate variables
      if (fhour.eq.0.0) then
        deallocate(rorgp)
        deallocate(s1,s1l)
        deallocate(qo)
        deallocate(si,sl)
        deallocate(tern)
      endif
      deallocate(g2recv,gr)
      call deallocate_gfsdata
!
      return
      end
