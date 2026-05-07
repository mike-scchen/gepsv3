      program mainlsmdata
!     implicit none
!c
      use physcons, rerth => con_rerth, pi => con_pi
      use machine , only : kind_io8,kind_io4
!
#include <paramodel.h>
#include <rscomloc.h>
!
! declaration:
! for terrain, inpt is 1206 :: istr=-2, iend=1203, datain is integer(kind=2)
!              the first and last 3 points is redundent.
! for BNU+NCU_top tile, inpt=1200 :: istr=1, iend=1200, data is integer
!
! input data
!      integer, parameter:: istr=1, iend=1200  ! BNU+NCU_top tile
!      integer, parameter:: inpt=iend-istr+1
!      integer, parameter:: itilex=1200,itiley=1200 ! pts of data used
!      real(kind=kind_io8), parameter:: dlon=0.00833333, dlat=0.00833333
!      real(kind=kind_io8), parameter:: pwest=-179.99583, psouth=-89.995833
      integer :: istr, iend,ar,xx,yy
      integer :: tile_x,tile_y,tile_z
      integer :: itilex,itiley,itilez  ! pts of data used
      integer :: recsize,recsize21,recsize22,recsize31,recsize88
      real(kind=kind_io8) :: dlon, dlat
      real(kind=kind_io8) :: pwest, psouth
      real(kind=kind_io8) :: c_min, c_max
      character :: projt*10
      character :: units*10
      real(kind=kind_io8) :: scale_fractor,missing_value
      integer :: wordsize
      real(kind=kind_io8) :: mask_n,NAN_val
      integer :: indx1,indx2,indy1,indy2
      integer :: ix1,ix2,iy1,iy2 
      integer :: ntile,ntilx,ntily
      integer :: pntilx,pntily,pntile
      character :: tiln*23
      character :: pathn*100,pathn2*100
      character :: filen*500,filen2*500
! array for reading input data
!      integer(kind=1) :: datain(istr:iend,istr:iend)
!      real(kind=kind_io8) :: slmskin(istr:iend,istr:iend)
      integer(kind=1),allocatable :: datain(:,:)
      integer(kind=1),allocatable :: slmskin(:,:)
! array for storing source data used
      real(kind=kind_io8),allocatable :: datar(:,:),datar2(:,:),slmskr(:,:)
      real(kind=kind_io8),allocatable :: rinlon(:),rinlat(:)
! dummy
      integer :: i,j,k,ii,jj,jz
      integer :: ntx,nty,ntz
      real(kind=kind_io8) :: aa,bb,cc,dd,ee,a0,b0,c0,d0,e0
      character :: outf*500,outf1*500,inlu*500
!      
!      integer,parameter :: inttyp=1    !inttyp=1  take the closest point value
      integer :: inttyp  !inttyp=1  take the closest point value
!      logical,parameter :: lmask=.false.
      logical :: lmask
      logical :: checkNAN
      logical :: ck0
! array for storing interpolation data and location
      integer, parameter:: nrslm=23
      real(kind=kind_io8), allocatable :: dataout(:),dataout2(:),flon(:),flat(:),slmskout(:),smsk(:,:)
      real(kind=kind_io8), allocatable :: flatr(:),slmskout2(:)
    
      real(kind=kind_io8) :: wlon, elon, nrlat, srlat
      real(kind=kind_io8) :: wlonin, nrlatin
      real(kind=kind_io8), parameter :: rad=pi / 180.

! array for grads output
      real(kind=kind_io4), allocatable :: datar4(:,:),dataout4(:)
      real(kind=kind_io4), allocatable :: slmskout4(:),slmskout24(:)

! MPI, no MPI currently    
      integer,parameter :: me=0
!-------------------
      character*120 condir
      character*7 dotype
!
      namelist /namloc/                                                      &
     & rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd      &
     &,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd      &
! This following line is redundent, but since rmtn, rlsmdata, and 3drtdata share
! the same namlist, the following 4 variables are kept in this rlsmdata code.
     &,clat1,clat2,clon1,clon2
!
      namelist /namcondir/ condir,pathn,pathn2
!      
      namelist /naminfo/ dotype,dlon,dlat,pwest,psouth,c_min,c_max,    &
                         projt,wordsize,tile_x,tile_y,tile_z,units,    &
                         scale_fractor,missing_value
!
! define regional grid from stdin and save to comloc
!
      read(5,namloc)
      read(5,namcondir)
      read(5,naminfo)

      print namloc 
      print namcondir
      print naminfo
!
! check input file information
!
      if(wordsize.ne.1) stop 'Err:wordsize of input data not 1!'
      if(projt.ne.'regular_ll') stop 'Err:projt of input data not ll!'
      do i=1,120
         if(condir(i:i).eq.' ') then
            mn=i
            exit
         endif
      enddo
! 
! info of landuse, soiltyp, vegfrac
!
      if (dotype.eq.'landuse') then
        istr=1
        iend=tile_x
        itilex=tile_x
        itiley=tile_y
        itilez=tile_z
        outf='out_landuse.bin'
        outf1='in_landuse.bin'
        if (units.eq.'category') then
          inttyp=1
          lmask=.true.
          mask_n=17.
        else
          stop 'error reading info of landuse data'
        endif
!       pathn='modis_landuse_20class_15s'
        print*,'pathn=',pathn
        NAN_val=missing_value
        ! 17 is water
        checkNAN=.false.
      elseif (dotype.eq.'soiltyp') then
        istr=1
        iend=tile_x
        itilex=tile_x
        itiley=tile_y
        itilez=tile_z
        outf='out_soiltype.bin'
        outf1='in_soiltype.bin'
        if (units.eq.'category') then
          inttyp=1
          lmask=.true.
          mask_n=14.
        else
          stop 'error reading info of soiltyp data'
        endif
!       pathn='bnu_soiltype_top+NCU'
        print*,'pathn=',pathn
        NAN_val=missing_value
        checkNAN=.false.
      elseif (dotype.eq.'vegfrac') then
        istr=1
        iend=tile_x
        itilex=tile_x
        itiley=tile_y
        itilez=tile_z
        outf='out_vegfrac.bin'
        outf1='in_vegfrac.bin'
        if (units.eq.'category') then
          stop 'error reading info of vegfrac data'
        else
          inttyp=0
          lmask=.true.
          mask_n=14.
        endif
!       pathn='greenfrac_fpar_modis'
        print*,'pathn=',pathn
        NAN_val=missing_value
        checkNAN=.true.
        ! maybe missing_value is ocean?
      else
        print*,'Error: type of surface data missing!'
        stop 8
      endif
      
          NAN_val=NAN_val*scale_fractor
          c_max=c_max*scale_fractor
          c_min=c_min*scale_fractor
      print*,'istr=',istr
      print*,'iend=',iend
!
!-------------start-------------------------------
! 
      do ntz=1,tile_z
! 
      allocate(dataout(lngrd),dataout2(lngrd),flon(lngrd),flat(lngrd),slmskout(lngrd))
      allocate(flatr(lngrd),slmskout2(lngrd))
      allocate(dataout4(lngrd),slmskout24(lngrd))
      slmskout=0.0
!      
! get output location
      print *,'+++++ define model grid lat/lon +++++'
      call setlatlon(flat,flon,delx,dely,dlamda0)
      flat(:)=flat(:)/rad
      flon(:)=flon(:)/rad
!      
! open file for lat information
      open(90,file='flatr.txt',form='formatted',status='unknown')
      open(91,file='flon.txt',form='formatted',status='unknown')
      write(91,*) flon(1:igrd1)
      close(91)
! reverse flat      
      do j = 1,jgrd1
      lats = j
      ijlats = (lats-1)*igrd1
      ijlatso = (lats)*igrd1
      jilats = (jgrd1-lats)*igrd1
      do i=1,igrd1
        flatr(ijlats+i)=flat(jilats+i)
      enddo
      write(90,*) flatr(ijlatso)
      enddo
      close(90)
      print *,'flatr=',(flatr(i),i=1,20)

!      
      wlon=flon(1)
      elon=flon(lngrd ) 
      nrlat=flatr(1) 
      srlat=flatr(lngrd)
      print*,'wlon=',wlon
      print*,'elon=',elon
      print*,'nrlat=',nrlat
      print*,'srlat=',srlat
! decide which tiles of source data are included, +-1 in case near tile boundary
      indx1=nint((wlon-pwest)/dlon)-1
      indx2=nint((elon-pwest)/dlon)+1
      indy1=nint((srlat-psouth)/dlat)-1
      indy2=nint((nrlat-psouth)/dlat)+1
      print*,'indx1,indx2,indy1,indy2=',   &
              indx1,indx2,indy1,indy2
      ntilx=int((indx2-indx1)/itilex)+2
      ntily=int((indy2-indy1)/itiley)+2
      ntile=ntilx*ntily
      pntilx=ntilx*itilex
      pntily=ntily*itiley
      pntile=ntile*itilex*itiley

      allocate(datar(pntilx,pntily),datar2(pntilx,pntily))
      allocate(datar4(pntilx,pntily))
      allocate(slmskr(pntilx,pntily))
      allocate(rinlon(pntilx))
      allocate(rinlat(pntily))

      slmskr=0.0
!      
!input file name = 36001-37200.13201-14400
!                   ix1   ix2   iy1   iy2
!

      print*,'recl=',ntilx,itilex,ntily,itiley,4
      print*,'recl=',ntilx*itilex*ntily*itiley*4

!get array size for file reclength      
      inquire(iolength=recsize21) dataout(1:lngrd)
      inquire(iolength=recsize22) dataout4(1:lngrd)
      inquire(iolength=recsize31) datar4(1:pntilx,1:pntily)
      inquire(iolength=recsize88) slmskout24(1:lngrd)
      print*,'recsize21=',recsize21
      print*,'recsize22=',recsize22
      print*,'recsize31=',recsize31
      print*,'recsize88=',recsize88


! open file for output data
      open(21,file=trim(outf),access='direct',form='unformatted',     &
             status='unknown',convert='big_endian',recl=recsize21)
      open(22,file=trim(outf)//'_4',access='direct',form='unformatted',     &
              status='unknown',convert='big_endian',recl=recsize22)
      open(31,file=trim(outf1),access='direct',form='unformatted',    &
             status='unknown',convert='big_endian',recl=recsize31)
      open(88,file='slmskout2.bin',access='direct',form='unformatted', &
              status='unknown',convert='big_endian',recl=recsize88)
!        
! --- loop for opening source tiles ---
      do nty=1,ntily
      do ntx=1,ntilx

        allocate(datain(istr:iend,istr:iend))
        inquire(iolength=recsize) datain(istr:iend,1)
        print*,'recsize=',recsize

        ix1=int(indx1/itilex)*itilex+1+itilex*(ntx-1)
        ix2=ix1-1+itilex
        iy1=int(indy2/itiley)*itiley+1-itiley*(nty-1)
        iy2=iy1-1+itiley
        write(tiln,'(I5.5,A1,I5.5,A1,I5.5,A1,I5.5)')                   &
                     ix1,'-',ix2,'.',iy1,'-',iy2
        print*,ntx,nty,tiln
        filen=condir(1:mn-1)//'/'//trim(pathn)//'/'//tiln
        print*,trim(filen)
        open(11,file=trim(filen),access='direct',form='unformatted',   &
             status='unknown',convert='big_endian',recl=recsize)
!       inquire(unit11,recl=recsize)
         if (dotype.eq.'vegfrac') then
         if (pathn.eq."veg_frac_lsa") then
             lmask=.false.
         else
        filen2=condir(1:mn-1)//'/'//trim(pathn2)//'/'//tiln
        print*,trim(filen2)
        open(12,file=trim(filen2),access='direct',form='unformatted',  &
             status='unknown',convert='big_endian',recl=recsize)
        allocate(slmskin(istr:iend,istr:iend),smsk(pntilx,pntily))
         endif
         endif

! read input data
        do j=1,tile_y
        jz=(ntz-1)*tile_y+j
          read(11,rec=jz)(datain(istr-1+i,istr-1+j),i=1,tile_x)
         if ((dotype.eq.'vegfrac') .and. (pathn.ne."veg_frac_lsa")) then
          read(12,rec=j)(slmskin(istr-1+i,istr-1+j),i=1,tile_x)
         endif
        enddo

        close(11)
        if((dotype.eq.'vegfrac').and.(pathn.ne."veg_frac_lsa"))&
             close(12)

! turn integer to real
        do i=1,itilex
          do j=1,itiley
          ii=(ntx-1)*itilex+i
          jj=(nty-1)*itiley+j
          datar(ii,jj)=DBLE(datain(i,iend+1-j))*scale_fractor
        if ((dotype.eq.'vegfrac') .and. (pathn.ne."veg_frac_lsa")) then
          smsk(ii,jj)=DBLE(slmskin(i,iend+1-j))
        endif
! sea land mask          
          if (lmask) then
          if (dotype.eq.'vegfrac') then
            if (smsk(ii,jj).eq.mask_n) then
              slmskr(ii,jj)=0.0
              datar(ii,jj)=0.0
            else
              slmskr(ii,jj)=1.0
             if(checkNAN) then
             if ((datar(ii,jj).gt.c_max).or.(datar(ii,jj).lt.c_min)) then
              datar(ii,jj)=NAN_val
             endif
             endif
            endif
          else
            if (datar(ii,jj).eq.mask_n ) then
              slmskr(ii,jj)=0.0
            else
              slmskr(ii,jj)=1.0
             if(checkNAN) then
             if ((datar(ii,jj).gt.c_max).or.(datar(ii,jj).lt.c_min)) then
              datar(ii,jj)=NAN_val
             endif
             endif
            endif
          endif
          endif
!          if(datar(ii,jj).le.-1.0*scale_fractor) then
!            datar(ii,jj)=NAN_val
!!            slmskr(ii,jj)=0.0
!!          else
!!            slmskr(ii,jj)=1.0
!          endif
          enddo
        enddo

!          
! check NAN values
        datar2=datar
        if(checkNAN) then
        do i=1,itilex
          do j=1,itiley
          ii=(ntx-1)*itilex+i
          jj=(nty-1)*itiley+j
            if(slmskr(ii,jj).eq.1.0 .and. datar(ii,jj).eq.NAN_val) then
!              print*,'missing value found in source data'
!              print*,'interpolation from nearest non-NAN 4 points'
              aa=0.
              bb=0.
              cc=0.
              dd=0.
              if(ii.gt.1)then
                 if(datar(ii-1,jj).ne.NAN_val) aa=1.
              endif
              if(ii.lt.itilex) then
                 if(datar(ii+1,jj).ne.NAN_val) bb=1.
              endif
              if(jj.gt.1)then
                 if(datar(ii,jj-1).ne.NAN_val) cc=1.
              endif
              if(jj.lt.itiley)then
                 if(datar(ii,jj+1).ne.NAN_val) dd=1.
              endif
              ee=aa+bb+cc+dd
!              a0=0.
!              b0=0.
!              c0=0.
!              d0=0.
!              if(datar(ii-1,jj).eq.0.0) a0=1.
!              if(datar(ii+1,jj).eq.0.0) b0=1.
!              if(datar(ii,jj-1).eq.0.0) c0=1.
!              if(datar(ii,jj+1).eq.0.0) d0=1.
!              e0=a0+b0+c0+d0
              if( ee .ne. 0. ) then
!inside                      
                if((ii.gt.1).and.(jj.gt.1).and.(ii.lt.itilex).and.(jj.lt.itiley)) then
                result=datar(ii-1,jj)*datar(ii+1,jj)*datar(ii,jj-1)*datar(ii,jj+1)
                if (result.eq.0.0) then
                ee=aa+bb+cc+dd
                datar2(ii,jj)=(datar(ii-1,jj)*aa+datar(ii+1,jj)*bb+     &
                              datar(ii,jj-1)*cc+datar(ii,jj+1)*dd)/ee
                endif
!corner                
                elseif ((ii+jj.eq.2)) then
                result=datar(ii+1,jj)*datar(ii,jj+1)
                if (result.eq.0.0) then
                ee=bb+dd
                datar2(ii,jj)=(datar(ii+1,jj)*bb+     &
                              datar(ii,jj+1)*dd)/ee
                endif
                elseif ((ii+jj.eq.1+itiley)) then
                result=datar(ii+1,jj)*datar(ii,jj-1)
                if (result.eq.0.0) then
                ee=bb+cc
                datar2(ii,jj)=(datar(ii+1,jj)*bb+     &
                              datar(ii,jj-1)*cc)/ee
                endif
                elseif ((ii+jj.eq.itilex+1)) then
                result=datar(ii-1,jj)*datar(ii,jj+1)
                if (result.eq.0.0) then
                ee=aa+dd
                datar2(ii,jj)=(datar(ii-1,jj)*aa+     &
                              datar(ii,jj+1)*dd)/ee
                endif
                elseif ((ii+jj.eq.itilex+itiley)) then
                result=datar(ii-1,jj)*datar(ii,jj-1)
                if (result.eq.0.0) then
                ee=aa+cc
                datar2(ii,jj)=(datar(ii-1,jj)*aa+     &
                              datar(ii,jj-1)*cc)/ee
                endif
!i=1                
                elseif ((ii.eq.1).and.(ii+jj.ne.2).and.(ii+jj.ne.1+itiley)) then
                result=datar(ii+1,jj)*datar(ii,jj-1)*datar(ii,jj+1)
                if (result.eq.0.0) then
                ee=bb+cc+dd
                datar2(ii,jj)=(datar(ii+1,jj)*bb+     &
                              datar(ii,jj-1)*cc+datar(ii,jj+1)*dd)/ee
                endif
!i=itilex         
                elseif ((ii.eq.itilex).and.(ii+jj.ne.itilex+1).and.(ii+jj.ne.itilex+itiley)) then
                result=datar(ii-1,jj)*datar(ii,jj-1)*datar(ii,jj+1)
                if (result.eq.0.0) then
                ee=aa+cc+dd
                datar2(ii,jj)=(datar(ii-1,jj)*aa+     &
                              datar(ii,jj-1)*cc+datar(ii,jj+1)*dd)/ee
                endif
!j=1                
                elseif ((jj.eq.1).and.(ii+jj.ne.2).and.(ii+jj.ne.1+itilex)) then
                result=datar(ii-1,jj)*datar(ii+1,jj)*datar(ii,jj+1)
                if (result.eq.0.0) then
                ee=aa+bb+dd
                datar2(ii,jj)=(datar(ii-1,jj)*aa+datar(ii+1,jj)*bb+     &
                              datar(ii,jj+1)*dd)/ee
                endif
!j=itiley           
                elseif ((jj.eq.itiley).and.(ii+jj.ne.1+itiley).and.(ii+jj.ne.itilex+itiley)) then
                result=datar(ii-1,jj)*datar(ii+1,jj)*datar(ii,jj-1)
                if (result.eq.0.0) then
                ee=aa+bb+cc
                datar2(ii,jj)=(datar(ii-1,jj)*aa+datar(ii+1,jj)*bb+     &
                              datar(ii,jj-1)*cc)/ee
                endif
                endif
              endif
              !ee\=0
            endif
            !slmsk,datar
          enddo
        enddo
        print*,'finish first part of NAN refill'

       datar=datar2
!       datar2=datar
!     ck0=.true.
      ck0=.false.
       if (.true.) then
!      do a=1,10
         jj=0
        do i=1,itilex/2
          ii=(ntx-1)*itilex+i
          jj=jj+1

          yy=jj
          do xx=ii,(itilex-ii+1)
            if((xx.gt.0).and.(yy.gt.0))  &
            datar2(xx,yy)=pt8_intp(datar,pntilx,pntily,xx,yy,NAN_val,slmskr,ck0)
          enddo

          xx=(itilex-ii+1)
          do yy=(jj+1),(itiley-jj+1)
            if((xx.gt.0).and.(yy.gt.0))  &
            datar2(xx,yy)=pt8_intp(datar,pntilx,pntily,xx,yy,NAN_val,slmskr,ck0)
          enddo

          yy=(itiley-jj+1)
          do xx=(itilex-ii),ii,-1
            if((xx.gt.0).and.(yy.gt.0))  &
            datar2(xx,yy)=pt8_intp(datar,pntilx,pntily,xx,yy,NAN_val,slmskr,ck0)
          enddo

          if (mod(itiley,2).eq.0 .and. jj.eq.(itiley/2))goto 100
          xx=ii
          do yy=(itiley-jj),(jj+1),-1
            if((xx.gt.0).and.(yy.gt.0))  &
            datar2(xx,yy)=pt8_intp(datar,pntilx,pntily,xx,yy,NAN_val,slmskr,ck0)
          enddo

        enddo
       if (mod(itilex,2).eq.1) then
           xx=itilex/2+1
           yy=itiley/2+1
            if((xx.gt.0).and.(yy.gt.0))  &
           datar2(xx,yy)=pt8_intp(datar,pntilx,pntily,xx,yy,NAN_val,slmskr,ck0)
       endif
 100   continue
!       enddo
       endif
       datar=datar2
       print*,'finish second part of NAN refill'

        endif
! endif check NAN values            
       print*,'finish check NAN values'
        if (dotype.eq.'vegfrac' .and. (pathn.ne."veg_frac_lsa"))      &
             deallocate(slmskin,smsk)
        deallocate(datain)

! get input location
        wlonin=pwest+dlon*DBLE(ix1)
        nrlatin=psouth+dlat*DBLE(iy2)
        print*,'wlonin=',wlonin,'nrlatin=',nrlatin
        do i=1,itilex
          ii=(ntx-1)*itilex+i
          rinlon(ii)=wlonin+dlon*DBLE(i-1)
        enddo
        do j=1,itiley
          jj=(nty-1)*itiley+j
          rinlat(jj)=nrlatin-dlat*DBLE(j-1)
        enddo

      enddo ! nty=1,ntily
      enddo ! ntx=1,ntilx
      print*,'finish all tiles'
!      
! --- end loop for opening source tiles ---
!
!      if (lmask) then
! get regional land-sea mask
        allocate(slmskout4(lngrd))
        read(nrslm) slmskout4
        slmskout=slmskout4
        do j=1,jgrd1
        do i=1,igrd1
          jj=(j-1)*igrd1+i
          jr=(jgrd1-j)*igrd1+i
          slmskout4(jj)=slmskout(jr)
        enddo
        enddo
        slmskout=0.0
        call sgl2ful(slmskout4,slmskout,lngrd)
        rewind nrslm
        deallocate(slmskout4)
!      endif

! do interpolation
       do j=1,pntily
       do i=1,pntilx
         if (datar(i,j).lt.0.0 .or.datar(i,j).eq.NAN_val ) slmskr(i,j)=0.0
       enddo
       enddo
       ar=0
       call la2ga(datar,pntilx,pntily,rinlon,rinlat,wlon,nrlat,        &
                  inttyp,dataout,lngrd,lmask,slmskr,slmskout,          &
                  flatr,flon,me,ar)

! mask used by slmskout
        do i=1,lngrd
        if (slmskout(i).eq.0.0) then
        if (dotype.eq.'vegfrac') then
            dataout(i)=0.0
        else
            dataout(i)=mask_n
        endif
        endif
        enddo

! print
      print*,'ntilx=',ntilx,'ntily=',ntily
      print*,'itilex=',itilex,'itiley=',itiley
      print*,'rinlon=',rinlon(1),'~',rinlon(ntilx*itilex)
      print*,'rinlat=',rinlat(1),'~',rinlat(ntily*itiley)
      print*,'flon=',flon(1),'~',flon(lngrd)
      print*,'flatr=',flatr(1),'~',flatr(lngrd)
! print
      print*,'input array=',ntilx*itilex,'x',ntily*itiley
      print*,'output array=',igrd1,'x',jgrd1
!      
! output data
      datar4(:,:)=datar(:,:)
      dataout4(:)=dataout(:)

      slmskout24(:)=slmskout(:)
      write(31,rec=ntz) datar4(1:pntilx,1:pntily)
!      write(21,rec=ntz) dataout4
!      write(21,rec=ntz) dataout
! reverse output data
      write(21,rec=ntz)(dataout((j-1)*igrd1+1:j*igrd1),j=jgrd1,1,-1)
      write(22,rec=ntz)(dataout4((j-1)*igrd1+1:j*igrd1),j=jgrd1,1,-1)
      write(88,rec=ntz)(slmskout24((j-1)*igrd1+1:j*igrd1),j=jgrd1,1,-1)

! close file
      close(21)
      close(22)
      close(31)
      close(88)
     
! deallocation of input data
      deallocate(rinlat)
      deallocate(rinlon)
      deallocate(slmskr)
      deallocate(datar4)
      deallocate(datar,datar2)
! deallocation of output data
      deallocate(dataout4,slmskout24)
      deallocate(flatr,slmskout2)
      deallocate(dataout,dataout2,flon,flat,slmskout)

      print *,' !---finish ntimes=',ntz,'---!'
      enddo !---end of ntz

      print *,' !!!! Normal exit !!!!'
      stop 0
      end program mainlsmdata

      function pt4_intp(ind,ix,iy,ii,jj,nan,slm,ck0)
      integer ii,jj,ix,iy
      real(kind=8) :: aa,bb,cc,dd,ee, nan
      real(kind=8) :: ind(ix,iy),slm(ix,iy)
      logical :: ck0
!              print*,'interpolation from nearest non-NAN 4 points'
         if (slm(ii,jj).eq.1.0 .and. ind(ii,jj).eq.nan)then
              aa=0.
              bb=0.
              cc=0.
              dd=0.
              if(ind(ii-1,jj).ne.nan .and. slm(ii-1,jj).eq.1.0) aa=1.
              if(ind(ii+1,jj).ne.nan .and. slm(ii+1,jj).eq.1.0) bb=1.
              if(ind(ii,jj-1).ne.nan .and. slm(ii,jj-1).eq.1.0) cc=1.
              if(ind(ii,jj+1).ne.nan .and. slm(ii,jj+1).eq.1.0) dd=1.
              ee=aa+bb+cc+dd
              if( ee .ne. 0. ) then

! check all 4 points are .gt. 0.0, if anynoe .eq. 0, then output .eq. 0
              if (ck0) then
              if (ind(ii-1,jj)*ind(ii+1,jj)*ind(ii,jj-1)*ind(ii,jj+1).eq.0.0) &
                pt4_intp=0.0
              else
                pt4_intp=(ind(ii-1,jj)*aa+ind(ii+1,jj)*bb+         &
                              ind(ii,jj-1)*cc+ind(ii,jj+1)*dd)/ee
              endif
              else
                pt4_intp=(ind(ii-1,jj)*aa+ind(ii+1,jj)*bb+         &
                              ind(ii,jj-1)*cc+ind(ii,jj+1)*dd)/ee
              endif
         else
            pt4_intp=ind(ii,jj)
         endif
      return
      end function

      function pt8_intp(ind,ix,iy,ii,jj,nan,slm,ck0)
      integer ii,jj,ix,iy
      real(kind=8) :: aa,bb,cc,dd,ee,ff,gg,hh,ss, nan
      real(kind=8) :: ind(ix,iy),slm(ix,iy)
      logical :: ck0
!              print*,'interpolation from nearest non-NAN 4 points'
         if (slm(ii,jj).eq.1.0 .and. &
             (ind(ii,jj).eq.nan .or. (ind(ii,jj).gt.0.0 .and. ind(ii,jj).le.0.1)))then
              aa=0.
              bb=0.
              cc=0.
              dd=0.
              ee=0.
              ff=0.
              gg=0.
              hh=0.
              ss=0.
              if (ii.ne.1) then
              if(ind(ii-1,jj).ne.nan .and. slm(ii-1,jj).eq.1.0) aa=1.
              endif
              if (ii.ne.ix) then
              if(ind(ii+1,jj).ne.nan .and. slm(ii+1,jj).eq.1.0) bb=1.
              endif
              if (jj.ne.1) then
              if(ind(ii,jj-1).ne.nan .and. slm(ii,jj-1).eq.1.0) cc=1.
              endif
              if (jj.ne.iy) then
              if(ind(ii,jj+1).ne.nan .and. slm(ii,jj+1).eq.1.0) dd=1.
              endif
              if ((ii.ne.1).and.(jj.ne.1)) then
              if(ind(ii-1,jj-1).ne.nan .and. slm(ii-1,jj-1).eq.1.0) ee=1.
              endif
              if ((ii.ne.ix).and.(jj.ne.iy)) then
              if(ind(ii+1,jj+1).ne.nan .and. slm(ii+1,jj+1).eq.1.0) ff=1.
              endif
              if ((ii.ne.ix).and.(jj.ne.1)) then
              if(ind(ii+1,jj-1).ne.nan .and. slm(ii+1,jj-1).eq.1.0) gg=1.
              endif
              if ((ii.ne.1).and.(jj.ne.iy)) then
              if(ind(ii-1,jj+1).ne.nan .and. slm(ii-1,jj+1).eq.1.0) hh=1.
              endif
              if (ii.eq.1) then
                 aa=0.
                 ee=0.
                 hh=0.
              elseif (ii.eq.ix) then
                 bb=0.
                 ff=0.
                 gg=0.
              endif
              if (jj.eq.1) then
                 cc=0.
                 ee=0.
                 gg=0.
              elseif (jj.eq.iy) then
                 dd=0.
                 ff=0.
                 hh=0.
              endif

              ss=aa+bb+cc+dd+ee+ff+gg+hh
              if( ss .ne. 0. ) then

! check all 8 points are .gt. 0.0, if anynoe .eq. 0, then output .eq. 0
!inside                      
               if((ii.ne.1).and.(jj.ne.1).and.(ii.ne.ix).and.(jj.ne.iy)) then
                 result= &
                   ind(ii-1,jj)*ind(ii+1,jj)*ind(ii,jj-1)*ind(ii,jj+1)* &
                   ind(ii-1,jj-1)*ind(ii+1,jj+1)*ind(ii+1,jj-1)*        &
                   ind(ii-1,jj+1)
               if (ck0.and.(result.eq.0.0))  then
                 pt8_intp=0.0
               else
                 ss=aa+bb+cc+dd+ee+ff+gg+hh
                 pt8_intp=(ind(ii-1,jj)*aa+ind(ii+1,jj)*bb+             &
                           ind(ii,jj-1)*cc+ind(ii,jj+1)*dd+             &
                           ind(ii-1,jj-1)*ee+ind(ii+1,jj+1)*ff+         &
                           ind(ii+1,jj-1)*gg+ind(ii-1,jj+1)*hh)/ss
               endif
!corner
               elseif((ii.eq.1).and.(jj.eq.1)) then
                 result= &
                   ind(ii+1,jj)*ind(ii,jj+1)*ind(ii+1,jj+1)
               if (ck0.and.(result.eq.0.0))  then
                 pt8_intp=0.0
               else
                 ss=bb+dd+ff
                 pt8_intp=(ind(ii+1,jj)*bb+             &
                           ind(ii,jj+1)*dd+             &
                           ind(ii+1,jj+1)*ff)/ss
               endif
               elseif((ii.eq.1).and.(jj.eq.iy)) then
                 result= &
                   ind(ii+1,jj)*ind(ii,jj-1)*ind(ii+1,jj-1)
               if (ck0.and.(result.eq.0.0))  then
                 pt8_intp=0.0
               else
                 ss=bb+cc+gg
                 pt8_intp=(ind(ii+1,jj)*bb+             &
                           ind(ii,jj-1)*cc+             &
                           ind(ii+1,jj-1)*gg)/ss
               endif
               elseif((ii.eq.ix).and.(jj.eq.1)) then
                 result= &
                   ind(ii-1,jj)*ind(ii,jj+1)*ind(ii-1,jj+1)
               if (ck0.and.(result.eq.0.0))  then
                 pt8_intp=0.0
               else
                 ss=aa+dd+hh
                 pt8_intp=(ind(ii-1,jj)*aa+             &
                           ind(ii,jj+1)*dd+             &
                           ind(ii-1,jj+1)*hh)/ss
               endif
               elseif((ii.eq.ix).and.(jj.eq.iy)) then
                 result= &
                   ind(ii-1,jj)*ind(ii,jj-1)*ind(ii-1,jj-1)
               if (ck0.and.(result.eq.0.0))  then
                 pt8_intp=0.0
               else
                 ss=aa+bb+cc+dd+ee+ff+gg+hh
                 pt8_intp=(ind(ii-1,jj)*aa+             &
                           ind(ii,jj-1)*cc+             &
                           ind(ii-1,jj-1)*ee)/ss
               endif
!ii=1
               elseif((ii.eq.1).and.(jj.ne.1).and.(jj.ne.iy)) then
                 result= &
                   ind(ii+1,jj)*ind(ii,jj-1)*ind(ii,jj+1)* &
                   ind(ii+1,jj+1)*ind(ii+1,jj-1)
               if (ck0.and.(result.eq.0.0))  then
                 pt8_intp=0.0
               else
                 ss=bb+cc+dd+ff+gg
                 pt8_intp=(ind(ii+1,jj)*bb+             &
                           ind(ii,jj-1)*cc+             &
                           ind(ii,jj+1)*dd+             &
                           ind(ii+1,jj+1)*ff+           &
                           ind(ii+1,jj-1)*gg)/ss
               endif
!ii=ix
               elseif((ii.eq.ix).and.(jj.ne.1).and.(jj.ne.iy)) then
                 result= &
                   ind(ii-1,jj)*ind(ii,jj-1)*ind(ii,jj+1)* &
                   ind(ii-1,jj-1)*ind(ii-1,jj+1)
               if (ck0.and.(result.eq.0.0))  then
                 pt8_intp=0.0
               else
                 ss=aa+cc+dd+ee+hh
                 pt8_intp=(ind(ii-1,jj)*aa+             &
                           ind(ii,jj-1)*cc+             &
                           ind(ii,jj+1)*dd+             &
                           ind(ii-1,jj-1)*ee+           &
                           ind(ii-1,jj+1)*hh)/ss
               endif
!jj=1
               elseif((jj.eq.1).and.(ii.ne.1).and.(ii.ne.ix)) then
                 result= &
                   ind(ii-1,jj)*ind(ii+1,jj)*ind(ii,jj+1)* &
                   ind(ii+1,jj+1)*ind(ii-1,jj+1)
               if (ck0.and.(result.eq.0.0))  then
                 pt8_intp=0.0
               else
                 ss=aa+bb+dd+ff+hh
                 pt8_intp=(ind(ii-1,jj)*aa+             &
                           ind(ii+1,jj)*bb+             &
                           ind(ii,jj+1)*dd+             &
                           ind(ii+1,jj+1)*ff+           &
                           ind(ii-1,jj+1)*hh)/ss
               endif
!jj=iy
               elseif((jj.eq.iy).and.(ii.ne.1).and.(ii.ne.ix)) then
                 result= &
                   ind(ii-1,jj)*ind(ii+1,jj)*ind(ii,jj-1)* &
                   ind(ii-1,jj-1)*ind(ii+1,jj-1)
               if (ck0.and.(result.eq.0.0))  then
                 pt8_intp=0.0
               else
                 ss=aa+bb+cc+dd+ee+ff+gg+hh
                 pt8_intp=(ind(ii-1,jj)*aa+             &
                           ind(ii+1,jj)*bb+             &
                           ind(ii,jj-1)*cc+             &
                           ind(ii-1,jj-1)*ee+         &
                           ind(ii+1,jj-1)*gg)/ss
               endif
               endif

              else
                pt8_intp=nan

              endif
         else
            pt8_intp=ind(ii,jj)
         endif
      return
      end function



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>>----modis_landuse_20class_15s info----->>      
!type=categorical
!category_min=1  (c_min)
!category_max=20 (c_max)
!projection=regular_ll (projt)
!dx=0.00416667 (dlon)
!dy=0.00416667 (dlat)
!known_x=1.0
!known_y=1.0
!known_lat=-89.9979167  (psouth)
!known_lon=-179.9979167 (pwest)
!wordsize=1   (wordsize)
!tile_x=2400  
!tile_y=2400
!tile_z=1
!units="category"
!description="MODIS modified-IGBP landuse - 500 meter"
!mminlu="MODIFIED_IGBP_MODIS_NOAH"
!iswater=17
!isice=15
!isurban=13
!<<----modis_landuse_20class_15s info-----<<      
