      program main_3drt
!     implicit none
!c
!
#include <paramodel.h>
#include <rscomloc.h>
!
! Discription of the Param_3DRT data:
! ---------------------------------------------------------------------
! In Param_3DRT data, 5 variables are used.
! They are in order of 
!     1.topo, 2.sinsl_cosas, 3.sinsl_sinas, 4.svf, 5.tcf
! Since the global netcdf4 data is too large, we cut it into a small 
! tile for MSM 1km Taiwan region [ about 18N~29N, 115E~126E ] via
! the code "cut_tile.F". The exact lat/lon is record in tile_lat.bin 
! and tile_lon.bin, and the data of this small tile is tile_3drt.bin.
! They are written in the following format:
!     open(iunit1,file='tile_3drt.bin',form='unformatted',
!    +     access='direct',convert='big_endian',recl=isize)
!     open(iunit2,file='tile_lon.bin',form='unformatted',
!    +     access='direct',convert='big_endian',recl=ilen))
!     open(iunit3,file='tile_lat.bin',form='unformatted',
!    +     access='direct',convert='big_endian',recl=jlen)
! ---------------------------------------------------------------------
!
! Subroutine:
!   "cut_tile.F":
!     cut a small tile and save it to "tile_3drt.bin" for MSM1km.
!     If the preexisting "tile_3drt.bin" is large enough, no need to run
!     this part again.
!
!   "rsm_interp.F":
!     read "tile_3drt.bin" and RSM/MSM's domain setting and interpolate
!     the data to RSM/MSM's grid. The interpolated data on RSM/MSM's grid
!     is "rsm_3drt.bin".
! 
! Declaration:
!-------------------
!
!     real    TS_LAT, TE_LAT, TS_LON, TE_LON
!     parameter (TS_LAT = 18.0, TE_LAT = 29.0 )
!     parameter (TS_LON = 115.0,TE_LON = 126.0)
      integer ilen, jlen

! The _FillValue (or say missing value is -9999.f in the netcdf file, check
! "ncdump -h" of the file)       
      real NAN_val
      parameter (NAN_val=-9999.0)
! There are 5 varibles we want to retrieve from the netcdf file, check "ncdump
! -h" of the file)
      integer nvar
      parameter (nvar=5)

      logical lcut
      parameter (lcut=.true.)
!
      namelist /namloc/                                                      &
     & rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd      &
     &,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd      &
! This following line is redundent, but since rmtn, rlsmdata, and 3drtdata share
! the same namlist, the following 4 variables are kept in this 3drtdata code.
     &,clat1,clat2,clon1,clon2
!      
      namelist /naminfo/ TS_LAT,TE_LAT,TS_LON,TE_LON,condir,FILE_NAME
!
! define regional grid from stdin and save to common (rscomloc.h)
!
      read(5,namloc)
      read(5,naminfo)

      print namloc 
      print naminfo
!
!-------------------
!

      print*,'--run cut_tile--'
      call cut_tile(ilen,jlen,lcut)

      print*,'--run rsm_interp--'
      call rsm_interp(ilen, jlen, nvar, NAN_val)

      stop 0 
      end program main_3drt
