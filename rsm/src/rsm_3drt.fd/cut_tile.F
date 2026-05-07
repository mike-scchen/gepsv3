!     This is part of the netCDF package.
!     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!     See COPYRIGHT file for conditions of use.

!     This is an example which reads some surface pressure and
!     temperatures. The data file read by this program is produced
!     comapnion program sfc_pres_temp_wr.f. It is intended to illustrate
!     the use of the netCDF fortran 77 API.

!     This program is part of the netCDF tutorial:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

!     Full documentation of the netCDF Fortran 77 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77

!     $Id: sfc_pres_temp_rd.f,v 1.8 2007/01/24 19:45:09 russ Exp $

!-----------------------------------------------------------------------------
      Subroutine cut_tile(ilen, jlen, lcut)
!
!     INPUT (from rscomloc.h):
!       TS_LAT, TE_LAT, TS_LON, TE_LON is the domain of interest. 
!       For MSM 1km, [ 18N~29N, 115E~126E ] is enough.
!     OUTPUT:
!       ilen: the length of index of zonal/x direction
!       jlen: the length of index of meridional/y direction
!
      use machine , only : kind_io8,kind_io4
!
      implicit none
#include <rscomloc.h>
      include 'netcdf.inc'

!     Logical variable for cutting data or not
      logical lcut

!     This is the name of the data file we will read.
!     character*(*) FILE_NAME
!     parameter (FILE_NAME='nc3.nc')
      character filen*500
      integer ncid

!     We are reading 2D data, a 17280 x 43200 lat-lon grid.
      integer NDIMS
      parameter (NDIMS=2)
      integer NLATS, NLONS
      parameter (NLATS = 17280, NLONS = 43200)
      character*(*) LAT_NAME, LON_NAME
      parameter (LAT_NAME='lat', LON_NAME='lon')
      integer lat_dimid, lon_dimid

!     For the lat lon coordinate netCDF variables.
      real(kind=kind_io4) lats(NLATS), lons(NLONS)
      integer lat_varid, lon_varid

!     We will read topography, sinsl_cosas, sinsl_sinas_in, svf and tcf
      character*(*) TOPO_NAME, SINSL_COSAS_NAME, SINSL_SINAS_NAME,      &
     &              SVF_NAME, TCF_NAME
      parameter (TOPO_NAME='topo')
      parameter (SINSL_COSAS_NAME='sinsl_cosas')
      parameter (SINSL_SINAS_NAME='sinsl_sinas')
      parameter (SVF_NAME='svf')
      parameter (TCF_NAME='tcf')
      integer topo_varid, sinsl_cosas_varid, sinsl_sinas_varid,         &
     &        svf_varid, tcf_varid
      integer dimids(NDIMS)

!     SKIP: To check the units attributes.
!     character*(*) UNITS
!     parameter (UNITS = 'units')
!     character*(*) TOPO_UNITS, SVF_UNITS, LAT_UNITS, LON_UNITS
!     parameter (TOPO_UNITS = ' ', SVF_UNITS = ' ')
!     parameter (LAT_UNITS = 'degrees_north')
!     parameter (LON_UNITS = 'degrees_east')
!     integer MAX_ATT_LEN
!     parameter (MAX_ATT_LEN = 80)
!     character*(MAX_ATT_LEN) topo_units_in, svf_units_in
!     character*(MAX_ATT_LEN) lat_units_in, lon_units_in
!     integer att_len

!     Read the data into these arrays.
      real(kind=kind_io4) :: topo_in(NLONS, NLATS)
      real(kind=kind_io4) :: sinsl_cosas_in(NLONS, NLATS)
      real(kind=kind_io4) :: sinsl_sinas_in(NLONS, NLATS)
      real(kind=kind_io4) :: svf_in(NLONS, NLATS)
      real(kind=kind_io4) :: tcf_in(NLONS, NLATS)
      real(kind=kind_io4),allocatable ::                                &
    &      datain(:,:),tmplons(:),tmplats(:)

!     Domain setting of the tile of interest. 
!     For MSM 1km, [ 18N~29N, 115E~126E ] is enough.
!     real TS_LAT, TE_LAT, TS_LON, TE_LON
!     parameter (TS_LAT = 18.0, TE_LAT = 29.0 )
!     parameter (TS_LON = 115.0,TE_LON = 126.0)
!     Also need indices of this small domain.
      integer ist, iend, jst, jend, ilen, jlen, isize
      integer ireclen, jreclen, irecsize
      integer iunit1, iunit2, iunit3

!     We will learn about the data file and store results in these
!     program variables.
      integer ndims_in, nvars_in, ngatts_in, unlimdimid_in

!     Loop indices
      integer lat, lon, i, mn

!     Error handling
      integer retval


!     Open the netcdf file. 
      do i=1,120
         if(condir(i:i).eq.' ') then
            mn=i
            exit
         endif
      enddo
      filen=condir(1:mn-1)//'/'//trim(FILE_NAME)
      retval = nf_open(trim(filen), nf_nowrite, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     There are a number of inquiry functions in netCDF which can be
!     used to learn about an unknown netCDF file. NF_INQ tells how many
!     netCDF variables, dimensions, and global attributes are in the
!     file; also the dimension id of the unlimited dimension, if there
!     is one.
      retval = nf_inq(ncid, ndims_in, nvars_in, ngatts_in,              &
     &     unlimdimid_in)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     In this case we know that there are 2 netCDF dimensions, 4 netCDF
!     variables, no global attributes, and no unlimited dimension.
!     if (ndims_in .ne. 2 .or. nvars_in .ne. 4 .or. ngatts_in .ne. 0    &
!    &     .or. unlimdimid_in .ne. -1) stop 2

!     Get the varids of the latitude and longitude coordinate variables.
      retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_inq_varid(ncid, LON_NAME, lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Read the latitude and longitude data.
      retval = nf_get_var_real(ncid, lat_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_var_real(ncid, lon_varid, lons)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Check to make sure we got what we expected.
!     do lat = 1, NLATS
!        if (lats(lat) .ne. START_LAT + (lat - 1) * 5.0) stop 2
!     end do
!     do lon = 1, NLONS
!        if (lons(lon) .ne. START_LON + (lon - 1) * 5.0) stop 2
!     end do
!     print*, lats(1)


!     Calculate the indices of the small tile
      ist =nint( (TS_LON-lons(1)) / (lons(2)-lons(1)) )
      iend=nint( (TE_LON-lons(1)) / (lons(2)-lons(1)) )
      jst =nint( (TS_LAT-lats(1)) / (lats(2)-lats(1)) )
      jend=nint( (TE_LAT-lats(1)) / (lats(2)-lats(1)) )
      ilen = iend-ist+1
      jlen = jend-jst+1
      isize= ilen*jlen
      print*,'isize=',isize
      print*,'ilen=',ilen
      print*,'jlen=',jlen
      inquire(iolength=irecsize) topo_in(ist:iend,jst:jend)
      inquire(iolength=ireclen) lons(ist:iend)
      inquire(iolength=jreclen) lats(jst:jend)
      print*,'irecsize=',irecsize
      print*,'ireclen=',ireclen
      print*,'jreclen=',jreclen

!     RETURN, if the preexisting "tile_3drt.bin" is large enough
      if (.not.lcut) return

!     Open the binary file for the small tile.
      iunit1=21
      iunit2=22
      iunit3=23
      open(iunit1,file=condir(1:mn-1)//'/'//'tile_3drt.bin',form='unformatted',              &
     &     access='direct',convert='big_endian',recl=irecsize)
      print*,'open 2nd file'
      open(iunit2,file=condir(1:mn-1)//'/'//'tile_lon.bin',form='unformatted',               &
     &     access='direct',convert='big_endian',recl=ireclen)
      print*,'open 3rd file'
      open(iunit3,file=condir(1:mn-1)//'/'//'tile_lat.bin',form='unformatted',               &
     &     access='direct',convert='big_endian',recl=jreclen)

      write(iunit2,rec=1) lons(ist:iend)
      write(iunit3,rec=1) lats(jst:jend)

!     print*,'lons=',lons(ist:iend)
!     print*,'lats=',lats(jst:jend)

!     Get the varids of the pressure and temperature netCDF variables.
      retval = nf_inq_varid(ncid, TOPO_NAME, topo_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, SINSL_COSAS_NAME, sinsl_cosas_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, SINSL_SINAS_NAME, sinsl_sinas_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, SVF_NAME, svf_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, TCF_NAME, tcf_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Read the surface pressure and temperature data from the file.
!     Since we know the contents of the file we know that the data
!     arrays in this program are the correct size to hold all the data.
      retval = nf_get_var_real(ncid, topo_varid, topo_in)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_get_var_real(ncid, sinsl_cosas_varid, sinsl_cosas_in)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_get_var_real(ncid, sinsl_sinas_varid, sinsl_sinas_in)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_get_var_real(ncid, svf_varid, svf_in)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_get_var_real(ncid, tcf_varid, tcf_in)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Check the data. It should be the same as the data we wrote.
!     do lon = 1, NLONS
!        do lat = 1, NLATS
!            if (pres_in(lon, lat) .ne. SAMPLE_PRESSURE +               &
!    &           (lon - 1) * NLATS + (lat - 1)) stop 2
!            if (temp_in(lon, lat) .ne. SAMPLE_TEMP +                   &
!    &           .25 * ((lon - 1) * NLATS + (lat - 1))) stop 2
!        end do
!     end do
!     print*,'topo(1,1)=',topo_in(1,1)
!     print*,'svf(1:10,1)=',svf_in(1:10,1)

!     Write data of the small tile to the file
!     1. topo
      write(iunit1,rec=1) topo_in(ist:iend,jst:jend)
!     2. sinsl_cosas
      write(iunit1,rec=2) sinsl_cosas_in(ist:iend,jst:jend)
!     3. sinsl_sinas
      write(iunit1,rec=3) sinsl_sinas_in(ist:iend,jst:jend)
!     4. svf
      write(iunit1,rec=4) svf_in(ist:iend,jst:jend)
!     5. tcf
      write(iunit1,rec=5) tcf_in(ist:iend,jst:jend)

!     Each of the netCDF variables has a "units" attribute. Let's read
!     them and check them.

!     retval = nf_get_att_text(ncid, lat_varid, UNITS, lat_units_in)
!     if (retval .ne. nf_noerr) call handle_err(retval)
!     retval = nf_inq_attlen(ncid, lat_varid, UNITS, att_len)
!     if (retval .ne. nf_noerr) call handle_err(retval)
!     if (lat_units_in(1:att_len) .ne. LAT_UNITS) stop 2
 
!     retval = nf_get_att_text(ncid, lon_varid, UNITS, lon_units_in)
!     if (retval .ne. nf_noerr) call handle_err(retval)
!     retval = nf_inq_attlen(ncid, lon_varid, UNITS, att_len)
!     if (retval .ne. nf_noerr) call handle_err(retval)
!     if (lon_units_in(1:att_len) .ne. LON_UNITS) stop 2

!     retval = nf_get_att_text(ncid, pres_varid, UNITS, pres_units_in)
!     if (retval .ne. nf_noerr) call handle_err(retval)
!     retval = nf_inq_attlen(ncid, pres_varid, UNITS, att_len)
!     if (retval .ne. nf_noerr) call handle_err(retval)
!     if (pres_units_in(1:att_len) .ne. PRES_UNITS) stop 2

!     retval = nf_get_att_text(ncid, temp_varid, UNITS, temp_units_in)
!     if (retval .ne. nf_noerr) call handle_err(retval)
!     retval = nf_inq_attlen(ncid, temp_varid, UNITS, att_len)
!     if (retval .ne. nf_noerr) call handle_err(retval)
!     if (temp_units_in(1:att_len) .ne. TEMP_UNITS) stop 2

!     Close the nc file. This frees up any internal netCDF resources
!     associated with the file.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Close the file of the small tile
      close(iunit1)
      close(iunit2)
      close(iunit3)

!     If we got this far, everything worked as expected. Yipee!
      print *,'*** SUCCESS reading Param_3DRT file Param_3drt.nc!'
      return
      end subroutine cut_tile

      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end subroutine handle_err
