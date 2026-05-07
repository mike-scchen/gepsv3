      subroutine rsm_interp(ilen,jlen,nvar,NAN_val)
!c
      use physcons, pi => con_pi
      use machine , only : kind_io8,kind_io4
      implicit none
!
#include <paramodel.h>
#include <rscomloc.h>
!
! ---------------------------------------------------------------------
!
!   "rsm_interp.F":
!     read "tile_3drt.bin" and RSM/MSM's domain setting and interpolate
!     the data to RSM/MSM's grid.
!
! 
! Declaration:
      real(kind=kind_io8), parameter :: rad=pi / 180.
! input data
      integer :: ilen, jlen, nvar
      real    :: NAN_val 
!>      
      real    :: valmax 
!<      
!      
! local var
! array for storing source data used
      real(kind=kind_io4),allocatable ::                                &
     &     datar(:,:,:),slmskr(:,:),rinlon(:),rinlat(:)
      real(kind=kind_io8),allocatable ::                                &
     &     datar8(:,:,:),slmskr8(:,:),rinlon8(:),rinlat8(:)
! array for storing interpolation data and location
      integer, parameter:: nrslm=23
      real(kind=kind_io8), allocatable ::                               &
     &     dataout(:),flon(:),flat(:),slmskout(:),flatr(:)
      real(kind=kind_io4), allocatable ::                               &
     &     dataout4(:),slmskout4(:)
! logical     
      logical,parameter :: lmask=.true.
! interpolation flag
      integer :: inttyp  !inttyp=1  take the closest point value
! real      
      real(kind=kind_io8) :: wlon, elon, nrlat, srlat
! read, mapping argument for setlatlon      
      real(kind=kind_io8) :: dlamda0
! integer, loop index
      integer i,j,jj,jr,nn
      integer lats,ijlats,ijlatso,jilats
! integer, size of the file rec length
      integer mn,ireclen, jreclen, irecsize, recsizeo, recsizeo4
! integer, file units
      integer iunit1,iunit2,iunit3,iunit4,iunit5,iunit6
!-------------start-------------------------------
! all 5 variables are continuous, no need to use closest point interpolation
! 5 vars: 1.topo, 2.sinsl_cosas, 3.sinsl_sinas, 4.svf, 5.tcf
      inttyp=0

      iunit1=21
      iunit2=22
      iunit3=23
      iunit4=24
      iunit5=25
      iunit6=26
! 
! allocation 
! array for storing source data used
      allocate(datar(ilen,jlen,nvar),slmskr(ilen,jlen))
      allocate(rinlon(ilen),rinlat(jlen))
      allocate(datar8(ilen,jlen,nvar),slmskr8(ilen,jlen))
      allocate(rinlon8(ilen),rinlat8(jlen))

! array for storing interpolation data and location
      allocate(dataout(lngrd),flon(lngrd),flat(lngrd),slmskout(lngrd))
      allocate(dataout4(lngrd))
      allocate(flatr(lngrd))
!      
! get output location (RSM's grid location)
      print *,'+++++ define model grid lat/lon +++++'
      print*,' do mercater'
      dlamda0 = 0.0
      call setlatlon(flat,flon,rdelx,rdely,dlamda0)
      flat(:)=flat(:)/rad
      flon(:)=flon(:)/rad
!      
! open file to output lat information
      do i=1,120
         if(condir(i:i).eq.' ') then
            mn=i
            exit
         endif
      enddo
      open(99,file=condir(1:mn-1)//'/'//'flon.txt',form='formatted',    &
     &     status='unknown')
      write(99,*) flon(1:igrd1)
      close(99)
! open file to output lat information
      open(98,file=condir(1:mn-1)//'/'//'flatr.txt',form='formatted',   &
     &     status='unknown')
! reverse flat      
      do j = 1,jgrd1
      lats = j
      ijlats = (lats-1)*igrd1
      ijlatso = (lats)*igrd1
      jilats = (jgrd1-lats)*igrd1
      do i=1,igrd1
        flatr(ijlats+i)=flat(jilats+i)
      enddo
      write(98,*) flatr(ijlatso)
      enddo
      close(98)
      print *,'flatr=',(flatr(i),i=1,20)

!      
      wlon=flon(1)
      elon=flon(lngrd ) 
      nrlat=flatr(1) 
      srlat=flatr(lngrd)
      print*,'wlon=',wlon,'elon=',elon,'nrlat=',nrlat,'srlat=',srlat

!      
!get array size for file reclength      
      inquire(iolength=recsizeo) dataout(1:lngrd)
      inquire(iolength=recsizeo4) dataout4(1:lngrd)
      print*,'recsizeo=',recsizeo,'recsizeo4=',recsizeo4,'lngrd=',lngrd

! open file for output data
      open(iunit1,file=condir(1:mn-1)//'/'//'rsm_3drt.bin',             &
     &     form='unformatted',access='direct',convert='big_endian',     &
     &     status='unknown',recl=recsizeo)
      open(iunit5,file=condir(1:mn-1)//'/'//'rsm_3drt_4.bin',           &
     &     form='unformatted',access='direct',convert='big_endian',     &
     &     status='unknown',recl=recsizeo4)


! open file of input data and read data
      inquire(iolength=irecsize)  datar(1:ilen,1:jlen,1)
      open(iunit2,file=condir(1:mn-1)//'/'//'tile_3drt.bin',            &
     &     form='unformatted',access='direct',convert='big_endian',     &
     &     status='old',recl=irecsize)
! output source data land sea mask
      open(iunit6,file=condir(1:mn-1)//'/'//'tile_3drt_lsmsk_4.bin',    &
     &     form='unformatted',access='direct',convert='big_endian',     &
     &     status='unknown',recl=irecsize)

      do nn=1,nvar
        read(iunit2,rec=nn) datar(:,:,nn)
      enddo
      close(iunit2)
!>
      valmax=-99999.0
      do j=1,jlen
      do i=1,ilen
          valmax=max(datar(i,j,1),valmax)
      enddo
      enddo
      print*,'valmax of topo=',valmax
!<

! get input location
      inquire(iolength=ireclen)  rinlon(:)
      inquire(iolength=jreclen)  rinlat(:)
      open(iunit3,file=condir(1:mn-1)//'/'//'tile_lon.bin',              &
     &     form='unformatted',access='direct',convert='big_endian',     &
     &     status='old',recl=ireclen)
      open(iunit4,file=condir(1:mn-1)//'/'//'tile_lat.bin',              &
     &     form='unformatted',access='direct',convert='big_endian',     &
     &     status='old',recl=jreclen)
      read(iunit3,rec=1) rinlon
      read(iunit4,rec=1) rinlat
      close(iunit3)
      close(iunit4)

!      
! get regional land-sea mask and reverse the lat direction
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

!
      !let it be real 8 and reverse the lat direction
      rinlon8=rinlon
      do j=1,jlen
        rinlat8(jlen+1-j)=rinlat(j)
      enddo
      do nn=1,nvar
      do j=1,jlen
      do i=1,ilen
        datar8(i,jlen+1-j,nn)=datar(i,j,nn)
      enddo
      enddo
      enddo
!>
      valmax=-99999.0
      do j=1,jlen
      do i=1,ilen
          valmax=max(datar8(i,j,1),valmax)
      enddo
      enddo
      print*,'valmax of topo=',valmax
!<

      print*,'nvar=',nvar
      do nn=1,nvar

      do j=1,jlen
      do i=1,ilen
        if (datar8(i,j,nn).eq.NAN_val) then
            slmskr8(i,j)=0.0
        else
            slmskr8(i,j)=1.0
        endif
      enddo
      enddo

! do interpolation
      call la2ga(datar8(1,1,nn),ilen,jlen,rinlon8,rinlat8,wlon,         &
     &           nrlat,inttyp,dataout,lngrd,lmask,slmskr8,slmskout,     &
     &           flatr,flon)

! mask used by slmskout
!       do i=1,lngrd
!       if (slmskout(i).eq.0.0) then
!           dataout(i)=0.0 !check different variables sea value
!       endif
!       enddo

!        
! reverse and output the data
      write(iunit1,rec=nn)(dataout((j-1)*igrd1+1:j*igrd1),j=jgrd1,1,-1)
      dataout4(:)=dataout(:)
      write(iunit5,rec=nn)(dataout4((j-1)*igrd1+1:j*igrd1),j=jgrd1,1,-1)

      print*,'end of nn=',nn
! enddo of nn=1,nvar 
      enddo
!>
      valmax=-99999.0
      do i=1,lngrd
          valmax=max(dataout(i),valmax)
      enddo
      print*,'valmax of topo=',valmax
!<
!
! write land sea mask of the source data
    
      slmskr=slmskr8
      write(iunit6,rec=1)(slmskr(:,j),j=jlen,1,-1)

! close file
      close(iunit1)
      close(iunit5)
      close(iunit6)
     
! deallocation 
      deallocate(flatr)
      deallocate(dataout4)
      deallocate(dataout,flon,flat,slmskout)
      deallocate(rinlon8,rinlat8)
      deallocate(datar8,slmskr8)
      deallocate(rinlon,rinlat)
      deallocate(datar,slmskr)


      print *,' !!!! Normal exit !!!!'
      return
      end subroutine rsm_interp
