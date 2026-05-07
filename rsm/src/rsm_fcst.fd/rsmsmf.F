      subroutine rsmsmf(fhnxt,compid)
!     subroutine rsmsmf(fhnxt,ggz,gq,gte,gdu,gzv,grq)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rsmsmf
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  do regional model forecast.
!
! program history log:
!
! usage:    call rsmsmf(gfhour,gshour,ggz,gq,gte,gdi,gze,grq)
!   input argument list:
!     gfhour    - global forecast hour
!     gshour    - global current time in second
!     ggz       - global terrain coefficient
!     gq        - global log surface pressure coefficient
!     gte       - global virtual temperature coefficient
!     gdi       - global divergence coefficient
!     gze       - global vorticity coefficient
!     grq       - global specific humidity coefficient
!
!   output argument list:
!
!   common block:
!     rscomf  rscomltb  comcon
!
!   input files: none
!
!   output files: none
!
!   subprograms called: 
!       rsmltb  rdfini  rmsgtr  rloopr  rbmcm  rloopa  rloopb  rdeldif
!       rsicdif  rlatbnd  rupdate  rfilt1  rfilt2
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
!.................................................................
! regional spectral model developed by hann-ming henry juang
!       version with nested to global spectral model
!.................................................................
!.....
!ccj
#include <paramodel.h>
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
#endif
#include <rscomfi.h>
#include <rscomspc.h>
#include <rscomgrd.h>
#include <rscomver.h>
#include <rscomio.h>
#include <rscompln.h>
#include <rscomrad.h>
#include <rscomloc.h>
#include <rscomltb.h>
#include <rscomgpd.h>
!cc
#ifdef MP
#define LNT22S lnt22p
#else
#define LNT22S lnt22
#endif
!
#ifdef G2R
      dimension                                                          &  
     &    ggz(LNT22S), gq(LNT22S)                                        &
     &  , gte(LNT22S,levs), grq(LNT22S,levs,ntotal)                      &
     &  , gdu(LNT22S,levs), gzv(LNT22S,levs)
#endif
#ifdef C2R
      dimension                                                          &
     &    ggz(clngrd), gq(clngrd)                                        &
     &  , gte(clngrd,levr), grq(clngrd,levr,ntotal)                      &
     &  , gdu(clngrd,levr), gzv(clngrd,levr)
#endif
!
!test
      dimension tmpsp(lnwav,levrp1)
      logical lsouf
      integer :: compid, ncpl, ierr
      real :: time_stamp, time_spend
!
#ifdef MP
        if (mype.eq.master) then
#endif
      print *,' ------ into rsmsmf ----- '
#ifdef MP
        endif
#endif
!
! ....... run regional per certain seconds
      secnow=phour*3600.
      secnxt=fhnxt*3600.
      if( secnxt.le.secnow ) return
#ifdef MP
        if (mype.eq.master) then
#endif
      print *,' run rsm from ',secnow,' to ',secnxt
#ifdef MP
        endif
#endif
!
! ....... setup the target lateral boundary for rsm
!     call rsmltb(secnow,secnxt,ggz,gq,gte,gdu,gzv,grq)
      call rsmltb(secnow,secnxt)
!
! ................................................
! ....... do step 1 if fhour=0 
#ifdef CPL_TIMCOM
  if(kdt .eq. 0) then
    if(mype.eq.master) print *,'[RSM] recv init data from TIMCOM Tai'
    call rsm_cpl_recv4rocn(compid, slmsk, tsea, ssu, ssv)
  end if
#endif
      if( fhour.eq.0.0 ) then
#ifndef MP
!yj fhdfi=0 by default        
        if( fhdfi.gt.0. ) then
          print*,'call rsmdfi'
          call rsmdfi
        endif
#endif
        kdt=kdt+1
        fhour=kdt*deltim/3600
        call rfwdstep
        phour=fhour
      endif
!
!....
!...............................................................
!...      smooth start time loop .............................
!
      lsimp=isemimp.eq.1
      lsdmp=idmpjet.eq.1
!     lsphy=imdlphy.eq.1
#ifdef A
      lsbgt=.true.
#endif
      lsfwd=.false.
      lsftr=.true.     ! the same as lsgen=.false.
      lssav=.true.
!
      dtltb=deltim
      dt2=2.*deltim
      idt=kdt+1
      mdt=nint(secnxt/deltim)
#ifdef A
      call rsetbgt(dt2,idt,mdt,nrbgt,nrbgt1,nrbgt2)
#endif
#ifdef CPL_TIMCOM
      if(kdt.eq.1) then
        if(mype.eq.master) print *,'[RSM] runtime step: ', kdt, '/', nsseg
        cpl_field(:,:,1) = u10m
        cpl_field(:,:,2) = v10m
        cpl_field(:,:,3) = t2m
        cpl_field(:,:,4) = q2m
        cpl_field(:,:,5) = psurf*1.d3
        cpl_field(:,:,6) =-fluxr(:,:,3)/3600.d0
        cpl_field(:,:,7) = fluxr(:,:,4)/3600.d0
        cpl_field(:,:,8) = dlwsfc/3600.d0
        cpl_field(:,:,9) = geshem*1.d3/deltim
        call rsm_cpl_send2rocn(compid, cpl_field, tsea) 
        cpl_field(:,:,5) = psurf
        cpl_field(:,:,6:9) = 0.d0
      end if
      time_stamp = MPI_WTIME() 
#endif
!
!     print *,' run step from ',idt,' to ',mdt
!
      do 10000 jdt=idt,mdt
        kdt=kdt+1
        fhour=kdt*deltim/3600
        lsout=mod(kdt,nsout).eq.0
        lscca=mod(kdt,nsswr).eq.0
        lsswr=mod(kdt,nsswr).eq.1
        lslwr=mod(kdt,nslwr).eq.1
        lsclm=mod(kdt,nsclm).eq.1
#ifdef RKN
        lsken=npoint.gt.0.and.mod(kdt,nsken).eq.0
        isave = 0
        if(lsken) isave = 1
#endif
!....
        call rsmstep
#ifdef CPL_TIMCOM
        if(mype.eq.master) print *,'[RSM] runtime step: ', kdt, '/', nsseg

        cpl_field(:,:,1) = cpl_field(:,:,1) + u10m
        cpl_field(:,:,2) = cpl_field(:,:,2) + v10m
        cpl_field(:,:,3) = cpl_field(:,:,3) + t2m
        cpl_field(:,:,4) = cpl_field(:,:,4) + q2m
        cpl_field(:,:,5) = cpl_field(:,:,5) + psurf

        if(mod(kdt,kdt_cpl).eq.0) then
          time_spend = MPI_WTIME() - time_stamp
          call MPI_ALLREDUCE(time_spend, time_stamp, 1, MPI_REAL8, MPI_MAX, mpi_comm_rsm, ierr)
          if(mype.eq.master) then
            write(*,'(a,f8.1,a)') '[RSM] runtime spends: ', time_stamp, "sec in a cpl loop"
          end if
   
          ncpl = kdt/kdt_cpl
          if(mype.eq.master) print *,'[RSM] couple step: ', kdt, '/', nsseg

          cpl_field(:,:,1:5) = cpl_field(:,:,1:5)/dble(kdt_cpl)
          cpl_field(:,:,5)   = cpl_field(:,:,5)*1.d3  ! Pa to hPa
          cpl_field(:,:,6)   =-odt_cpl*(fluxr(:,:,3) - cpl_field(:,:,6))
          cpl_field(:,:,7)   = odt_cpl*(fluxr(:,:,4) - cpl_field(:,:,7))
          cpl_field(:,:,8)   = odt_cpl*(dlwsfc - cpl_field(:,:,8))
          cpl_field(:,:,9)   = odt_cpl*(geshem - cpl_field(:,:,9))*1.d3

          call rsm_cpl_recv4rocn(compid, slmsk, tsea, ssu, ssv)
          call rsm_cpl_send2rocn(compid, cpl_field, tsea)

          if(lsout) then
            cpl_field = 0.d0
          else
            cpl_field(:,:,1:5) = 0.d0
            cpl_field(:,:,6)   = fluxr(:,:,3)
            cpl_field(:,:,7)   = fluxr(:,:,4)
            cpl_field(:,:,8)   = dlwsfc
            cpl_field(:,:,9)   = geshem
          end if
          
          time_stamp = MPI_WTIME()
        end if
#endif
!>>-------------------------debug
!       if (abs(fhour-1.2083).le.1.e-4) then
!         call rwrtout
!         stop 88
!       endif
!<<-------------------------debug
!.....
! intermedium output
        if(lsout) then
        lsouf=mod(kdt,nsseg).eq.0
       if(mype.eq.master) print *,'in rsmsmf lsout=',lsout,'kdt=',kdt       &   
     &  ,'nsseg=',nsseg,'nsout=',nsout,'lsouf=',lsouf
!
        if(.not.lsouf) then
#ifdef RKN
#ifdef MP
          call rmpkencollect(itnum,nvrken,nptken,nstken,svdata)
#endif
#endif
          call rwrtout
          call rzerout
#ifdef RKN
!
! reset itnum after writing point sounding data.
!
          itnum=1
          do k=1,nstken
            do j=1,nptken
              do i=1,nvrken
                svdata(i,j,k) = 0.
              enddo
            enddo
          enddo
#endif
!
        endif
        endif
        phour=fhour
10000 continue
!................................................
!          time loop end
!................................................
#ifdef A
      call rwrtbgt(fhour,idate)
#endif
!................................................
!
      return
      end
