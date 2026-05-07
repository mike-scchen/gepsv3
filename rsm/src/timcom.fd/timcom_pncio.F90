module tai_timcom_pncio
  use tai_timcom_const
  use tai_timcom_comm
  use tai_timcom_grid, only: ocn_panel
  use PNETCDF
  implicit none

contains
subroutine get_lateral_bc_all(gmpi, pncfile, wbc, ebc, sbc, nbc)
  use tai_hyperlink, only: myid
  implicit none

  type(comm_panel), intent(in) :: gmpi
  character(len=*) :: pncfile
  real(r8), dimension(:,:,:,:), pointer :: wbc, ebc, sbc, nbc
  integer :: nlon, nlat, nlev, ntime, ncid, ierr
  integer(kind=MPI_OFFSET_KIND) :: ns_lbc(3),nc_lbc(3)

  nlat = size(wbc,1)
  nlon = size(sbc,1)
  nlev = size(wbc,2)
  ntime = size(sbc,3)

  ierr = NFMPI_OPEN(gmpi%m_comm_cart, trim(pncfile), NF_NOWRITE, MPI_INFO_NULL, ncid)

  ns_lbc(1) = gmpi%myid_y*nlat + 1
  ns_lbc(2) = 1
  ns_lbc(3) = 1
  nc_lbc(1) = nlat
  nc_lbc(2) = nlev
  nc_lbc(3) = ntime

  call get_nfmpi_vara_double(ncid, 'U1_BC_E', ns_lbc, nc_lbc, ebc, 1)
  call get_nfmpi_vara_double(ncid, 'V1_BC_E', ns_lbc, nc_lbc, ebc, 2)
  call get_nfmpi_vara_double(ncid, 'T1_BC_E', ns_lbc, nc_lbc, ebc, 3)
!  if(myid .eq. 1) write(*,*) 'get_nfmpi_vara_double ebc=', ebc(5,1,1,3)
  call get_nfmpi_vara_double(ncid, 'S1_BC_E', ns_lbc, nc_lbc, ebc, 4)
  call get_nfmpi_vara_double(ncid, 'U1_BC_W', ns_lbc, nc_lbc, wbc, 1)
  call get_nfmpi_vara_double(ncid, 'V1_BC_W', ns_lbc, nc_lbc, wbc, 2)
  call get_nfmpi_vara_double(ncid, 'T1_BC_W', ns_lbc, nc_lbc, wbc, 3)
  call get_nfmpi_vara_double(ncid, 'S1_BC_W', ns_lbc, nc_lbc, wbc, 4)

  ns_lbc(1) = gmpi%myid_x*nlon + 1
  nc_lbc(1) = nlon

  call get_nfmpi_vara_double(ncid, 'U1_BC_N', ns_lbc, nc_lbc, nbc, 1)
  call get_nfmpi_vara_double(ncid, 'V1_BC_N', ns_lbc, nc_lbc, nbc, 2)
  call get_nfmpi_vara_double(ncid, 'T1_BC_N', ns_lbc, nc_lbc, nbc, 3)
  call get_nfmpi_vara_double(ncid, 'S1_BC_N', ns_lbc, nc_lbc, nbc, 4)
  call get_nfmpi_vara_double(ncid, 'U1_BC_S', ns_lbc, nc_lbc, sbc, 1)
  call get_nfmpi_vara_double(ncid, 'V1_BC_S', ns_lbc, nc_lbc, sbc, 2)
  call get_nfmpi_vara_double(ncid, 'T1_BC_S', ns_lbc, nc_lbc, sbc, 3)
  call get_nfmpi_vara_double(ncid, 'S1_BC_S', ns_lbc, nc_lbc, sbc, 4)

  ierr = NFMPI_CLOSE(ncid)

end subroutine get_lateral_bc_all 

subroutine get_nfmpi_vara_double(ncid, varname, ns_lbc, nc_lbc, lbc_all, ind)
  implicit none

  integer :: ncid, ind
  character(len=*) :: varname
  integer(kind=MPI_OFFSET_KIND) :: ns_lbc(3),nc_lbc(3)
  real(r8), dimension(:,:,:,:), pointer :: lbc_all
  real(r8) :: lbc(nc_lbc(1),nc_lbc(2),nc_lbc(3))
  integer :: ierr, varid

  ierr = NFMPI_INQ_VARID(ncid, trim(varname), varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_lbc, nc_lbc, lbc)
  lbc(:,1,:) = lbc(:,2,:)
  lbc(:,nc_lbc(2),:) = lbc(:,nc_lbc(2)-1,:)
  lbc_all(:,:,:,ind) = lbc

end subroutine get_nfmpi_vara_double

subroutine get_indata_pncio(pncfile)
  use tai_hyperlink
  implicit none

  character(len=*) :: pncfile
  integer :: ncid, ierr, varid

  ierr = NFMPI_OPEN(m_comm_cart, trim(pncfile), NF_NOWRITE, MPI_INFO_NULL, ncid)
  ierr = NFMPI_INQ_VARID(ncid, 'x_grid', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1), nc_ijk(1), x_grid)
  ierr = NFMPI_INQ_VARID(ncid, 'y_grid', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ijk(2), y_grid)
  ierr = NFMPI_INQ_VARID(ncid, 'z_grid', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(3), nc_ijk(3), z_grid)
  ierr = NFMPI_INQ_VARID(ncid, 'x_face', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1), nc_ujk(1), x_face)
  ierr = NFMPI_INQ_VARID(ncid, 'y_face', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ivk(2), y_face)
  ierr = NFMPI_INQ_VARID(ncid, 'z_face', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(3), nc_ijw(3), z_face)
  ierr = NFMPI_INQ_VARID(ncid, 'cs', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ijk(2), cs)
  ierr = NFMPI_INQ_VARID(ncid, 'csv', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ivk(2), csv)
  ierr = NFMPI_INQ_VARID(ncid, 'ocs', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ijk(2), ocs)
  ierr = NFMPI_INQ_VARID(ncid, 'ocsv', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ivk(2), ocsv)
  ierr = NFMPI_INQ_VARID(ncid, 'dx', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ijk(2), dx)
  ierr = NFMPI_INQ_VARID(ncid, 'dxu', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ivk(2), dxu)
  ierr = NFMPI_INQ_VARID(ncid, 'odx', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ijk(2), odx)
  ierr = NFMPI_INQ_VARID(ncid, 'odxu', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ivk(2), odxu)
  ierr = NFMPI_INQ_VARID(ncid, 'dy', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ijk(2), dy)
  ierr = NFMPI_INQ_VARID(ncid, 'dyv', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ivk(2), dyv)
  ierr = NFMPI_INQ_VARID(ncid, 'ody', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ijk(2), ody)
  ierr = NFMPI_INQ_VARID(ncid, 'odyv', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(2), nc_ivk(2), odyv)
  ierr = NFMPI_INQ_VARID(ncid, 'odz', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(3), nc_ijk(3), odz)
  ierr = NFMPI_INQ_VARID(ncid, 'odzw', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(3), nc_ivk(3), odzw)
  !ierr = NFMPI_INQ_VARID(ncid, 'area', varid)
  !ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), area)
  !ierr = NFMPI_INQ_VARID(ncid, 'volume', varid)
  !ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:3), nc_ijk(1:3), volume)
  ierr = NFMPI_INQ_VARID(ncid, 'kb', varid)
  ierr = NFMPI_GET_VARA_INT2_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), kb(1:nx,1:ny))
  ierr = NFMPI_INQ_VARID(ncid, 'in', varid)
  ierr = NFMPI_GET_VARA_INT2_ALL(ncid, varid, ns_ijk(1:3), nc_ijk(1:3), in(1:nx,1:ny,1:nz))
  !ierr = NFMPI_INQ_VARID(ncid, 'iu', varid)
  !ierr = NFMPI_GET_VARA_INT2_ALL(ncid, varid, ns_ijk(1:3), nc_ujk(1:3), iu(1:nxf,1:ny,1:nz))
  !ierr = NFMPI_INQ_VARID(ncid, 'iv', varid)
  !ierr = NFMPI_GET_VARA_INT2_ALL(ncid, varid, ns_ijk(1:3), nc_ivk(1:3), iv(1:nx,1:nyf,1:nz))
  !ierr = NFMPI_INQ_VARID(ncid, 'iw', varid)
  !ierr = NFMPI_GET_VARA_INT2_ALL(ncid, varid, ns_ijk(1:3), nc_ijw(1:3), iw(1:nx,1:ny,1:nzf))

  ierr = NFMPI_CLOSE(ncid)
end subroutine get_indata_pncio

subroutine get_ssh_pncio(ssh_file)
  use tai_hyperlink, only: m_comm_cart, nx, ny, ns_ijk, nc_ijk, ssh
  implicit none

  character(len=*) :: ssh_file
  integer :: ncid, ierr, varid
  real(r4) :: var2d(nx,ny)

  var2d = 0.0
  ierr = NFMPI_OPEN(m_comm_cart, trim(ssh_file), NF_NOWRITE, MPI_INFO_NULL, ncid)
  ierr = NFMPI_INQ_VARID(ncid, 'ssh', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), var2d)
  ssh(1:nx,1:ny) = dble(var2d)

  ierr = NFMPI_CLOSE(ncid)
end subroutine get_ssh_pncio

subroutine get_ic_pncio(ic_file)
  use tai_hyperlink, only: m_comm_cart, nx, ny, nz, ns_ijk, nc_ijk, &
                       u1, v1, t1, s1, ssh
  implicit none
  
  character(len=*) :: ic_file
  integer :: ncid, ierr, varid
  real(r4) :: var3d(nx,ny,nz)

  ierr = NFMPI_OPEN(m_comm_cart, trim(ic_file), NF_NOWRITE, MPI_INFO_NULL, ncid)
  ierr = NFMPI_INQ_VARID(ncid, 'U1', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_ijk, nc_ijk, var3d)
  u1(1:nx,1:ny,1:nz) = dble(var3d)
  ierr = NFMPI_INQ_VARID(ncid, 'V1', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_ijk, nc_ijk, var3d)
  v1(1:nx,1:ny,1:nz) = dble(var3d)
  ierr = NFMPI_INQ_VARID(ncid, 'T1', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_ijk, nc_ijk, var3d)
  t1(1:nx,1:ny,1:nz) = dble(var3d)
  ierr = NFMPI_INQ_VARID(ncid, 'S1', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_ijk, nc_ijk, var3d)
  s1(1:nx,1:ny,1:nz) = dble(var3d)
  ierr = NFMPI_INQ_VARID(ncid, 'SSH', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), var3d(:,:,1))
  ssh(1:nx,1:ny) = dble(var3d(:,:,1))
  ierr = NFMPI_CLOSE(ncid)
end subroutine get_ic_pncio

subroutine put_ic_pncio(pncfile)
  use tai_hyperlink, only: m_comm_cart, nx, ny, nz, ns_ijk, nc_ijk, &
                       nx_grid, ny_grid, nz_grid, &
                        x_grid,  y_grid,  z_grid, &
                       u1, v1, t1, s1, rho, p0, in, u, v, w
  implicit none

  character(len=*) :: pncfile
  integer(MPI_OFFSET_KIND) :: ngdim(3), length
  integer :: ncid, iid, jid, kid, nid(3), ierr
  integer :: id_lon, id_lat, id_lev, &
             id_u1, id_v1, id_t1, id_s1, &
             id_im, id_kb, id_rho, id_p0, id_u, id_v, id_w
  real(r4) :: lon(nx), lat(ny), dep(nz), tmp3d(nx,ny,nz), tmp2d(nx,ny), zero(1)

  zero = 0.d0
  ngdim(1) = nx_grid
  ngdim(2) = ny_grid
  ngdim(3) = nz_grid

  ierr = NFMPI_CREATE(m_comm_cart, trim(pncfile), NF_64BIT_OFFSET, MPI_INFO_NULL, ncid)

  ierr = NFMPI_DEF_DIM(ncid, "lon", ngdim(1), nid(1))
  ierr = NFMPI_DEF_DIM(ncid, "lat", ngdim(2), nid(2))
  ierr = NFMPI_DEF_DIM(ncid, "dep", ngdim(3), nid(3))

  ierr = NFMPI_DEF_VAR(ncid, "lon", NF_REAL, 1, nid(1), id_lon)
  length = len("degrees_east")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lon, 'units', length, "degrees_east")
  ierr = NFMPI_DEF_VAR(ncid, "lat", NF_REAL, 1, nid(2), id_lat)
  length = len("degrees_north")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lat, 'units', length, "degrees_north")
  ierr = NFMPI_DEF_VAR(ncid, "dep", NF_REAL, 1, nid(3), id_lev)
  length = len("meters")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lev, 'units', length, "meters")

  ierr = NFMPI_DEF_VAR(ncid, "U1", NF_REAL, 3, nid, id_u1)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u1, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u1, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u1, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V1", NF_REAL, 3, nid, id_v1)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v1, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v1, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v1, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "T1", NF_REAL, 3, nid, id_t1)
  length = len("C")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t1, 'units', length, "C")
  length = len("Water Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t1, 'long_name', length, "Water Temperature")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_t1, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "S1", NF_REAL, 3, nid, id_s1)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s1, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s1, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_s1, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "RHO", NF_REAL, 3, nid, id_rho)
  length = len("g/cm^3")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rho, 'units', length, "g/cm^3")
  length = len("density")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rho, 'long_name', length, "density")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_rho, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "P0", NF_REAL, 2, nid, id_p0)
  length = len("Pa")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rho, 'units', length, "Pa")
  length = len("Surface Pressure")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rho, 'long_name', length, "Surface Pressure")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_rho, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "U", NF_REAL, 3, nid, id_u)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V", NF_REAL, 3, nid, id_v)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "W", NF_REAL, 3, nid, id_w)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_w, 'units', length, "cm s-1")
  length = len("Vertical Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_w, 'long_name', length, "Vertical Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_w, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "IN", NF_SHORT, 3, nid, id_im)
  length = len("IN")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_im, 'units', length, "IN")

  ierr = NFMPI_ENDDEF(ncid)

  lon = real(x_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lon, ns_ijk(1), nc_ijk(1), lon)
  lat = real(y_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lat, ns_ijk(2), nc_ijk(2), lat)
  dep = real(z_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lev, ns_ijk(3), nc_ijk(3), dep)
  tmp3d = real(u1(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_u1,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(v1(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_v1,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(t1(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_t1,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(s1(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_s1,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(rho(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_rho, ns_ijk,    nc_ijk, tmp3d)
  tmp2d = real(p0(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_p0,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp3d = real(u(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_u,   ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(v(1:nx,2:ny+1,1:nz)) 
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_v,   ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(w(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_w,   ns_ijk,    nc_ijk, tmp3d)
  ierr = NFMPI_PUT_VARA_INT2_ALL(ncid, id_im,  ns_ijk,    nc_ijk, in(1:nx,1:ny,1:nz))
  
  ierr = NFMPI_CLOSE(ncid)
end subroutine put_ic_pncio

subroutine put_flux_pncio(pncfile)
  use tai_hyperlink, only: m_comm_cart, nx, ny, nz, ns_ijk, nc_ijk, &
                       nx_grid, ny_grid, nz_grid, &
                        x_grid,  y_grid,  z_grid, &
                       u_10, v_10, t_10, q_10, pslv, &
                       taux, tauy, swup, swdn, &
                       lwup, lwdn, lath, senh, &
                       evap, rain, snow, ifrc, roff
  implicit none

  character(len=*) :: pncfile
  integer(MPI_OFFSET_KIND) :: ngdim(2), length
  integer :: ncid, iid, jid, kid, nid(2), ierr
  integer :: id_lon, id_lat, id_lev, id_u10, id_v10, &
             id_t10, id_q10, id_slp, id_swu, id_swd, &
             id_lwd, id_ran, id_snw, id_ice, id_rof, &
             id_tux, id_tuy, id_lwu, id_lth, id_snh, id_evp
  real(r4) :: lon(nx), lat(ny), tmp2d(nx,ny), zero(1)

  zero = 0.d0
  ngdim(1) = nx_grid
  ngdim(2) = ny_grid

  ierr = NFMPI_CREATE(m_comm_cart, trim(pncfile), NF_64BIT_OFFSET, MPI_INFO_NULL, ncid)

  ierr = NFMPI_DEF_DIM(ncid, "lon", ngdim(1), nid(1))
  ierr = NFMPI_DEF_DIM(ncid, "lat", ngdim(2), nid(2))

  ierr = NFMPI_DEF_VAR(ncid, "lon", NF_REAL, 1, nid(1), id_lon)
  length = len("degrees_east")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lon, 'units', length, "degrees_east")

  ierr = NFMPI_DEF_VAR(ncid, "lat", NF_REAL, 1, nid(2), id_lat)
  length = len("degrees_north")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lat, 'units', length, "degrees_north")

  ierr = NFMPI_DEF_VAR(ncid, "U_10", NF_REAL, 2, nid, id_u10)
  length = len("m/s")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u10, 'units', length, "m/s")
  length = len("10 m Eastward wind velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u10, 'long_name', length, "10 m Eastward wind velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u10, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V_10", NF_REAL, 2, nid, id_v10)
  length = len("m/s")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v10, 'units', length, "m/s")
  length = len("10 m Northward wind velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v10, 'long_name', length, "10 m Northward wind velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v10, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "T_10", NF_REAL, 2, nid, id_t10)
  length = len("degree")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t10, 'units', length, "degree")
  length = len("10 m height air temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t10, 'long_name', length, "10 m height air temperature")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_t10, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "Q_10", NF_REAL, 2, nid, id_q10)
  length = len("%/%")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_q10, 'units', length, "%/%")
  length = len("10 m height specific humidity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_q10, 'long_name', length, "10 m height specific humidity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_q10, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "PSLV", NF_REAL, 2, nid, id_slp)
  length = len("Pa")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_slp, 'units', length, "Pa")
  length = len("sea level pressure")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_slp, 'long_name', length, "sea level pressure")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_slp, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "SWUP", NF_REAL, 2, nid, id_swu)
  length = len("W/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_swu, 'units', length, "W/m2")
  length = len("Upward short-wave radiation")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_swu, 'long_name', length, "Upward short-wave radiation")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_swu, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "SWDN", NF_REAL, 2, nid, id_swd)
  length = len("W/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_swd, 'units', length, "W/m2")
  length = len("Downward short-wave radiation")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_swd, 'long_name', length, "Downward short-wave radiation")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_swd, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "LWDN", NF_REAL, 2, nid, id_lwd)
  length = len("W/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lwd, 'units', length, "W/m2")
  length = len("Downward long-wave radiation")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lwd, 'long_name', length, "Downward long-wave radiation")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_lwd, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "RAIN", NF_REAL, 2, nid, id_ran)
  length = len("kg/sec/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ran, 'units', length, "kg/sec/m^2")
  length = len("rain")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ran, 'long_name', length, "rain")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_ran, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "SNOW", NF_REAL, 2, nid, id_snw)
  length = len("kg/sec/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_snw, 'units', length, "kg/sec/m^2")
  length = len("snow")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_snw, 'long_name', length, "snow")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_snw, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "ifrac", NF_REAL, 2, nid, id_ice)
  length = len("[0-1]")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ice, 'units', length, "[0-1]")
  length = len("ice fraction")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ice, 'long_name', length, "ice fraction")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_ice, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "roff", NF_REAL, 2, nid, id_rof)
  length = len("kg/sce/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rof, 'units', length, "kg/sec/m^2")
  length = len("runoff")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rof, 'long_name', length, "runoff")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_rof, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "taux", NF_REAL, 2, nid, id_tux)
  length = len("dyne/cm^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_tux, 'units', length, "dyne/cm^2")
  length = len("longitude wind stress")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_tux, 'long_name', length, "longitude wind stress")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_tux, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "tauy", NF_REAL, 2, nid, id_tuy)
  length = len("dyne/cm^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_tuy, 'units', length, "dyne/cm^2")
  length = len("latitude wind stress")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_tuy, 'long_name', length, "latitude wind stress")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_tuy, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "LWUP", NF_REAL, 2, nid, id_lwu)
  length = len("Watt/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lwu, 'units', length, "Watt/m^2")
  length = len("upward long wave radiation")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lwu, 'long_name', length, "upward long wave radiation")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_lwu, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "LATH", NF_REAL, 2, nid, id_lth)
  length = len("Watt/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lth, 'units', length, "Watt/m^2")
  length = len("Latent heat flux")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lth, 'long_name', length, "Latent heat flux")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_lth, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "SENH", NF_REAL, 2, nid, id_snh)
  length = len("Watt/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_snh, 'units', length, "Watt/m^2")
  length = len("Sensible heat flux")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_snh, 'long_name', length, "Sensible heat flux")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_snh, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "EVAP", NF_REAL, 2, nid, id_evp)
  length = len("Kg/sec/m^2")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_evp, 'units', length, "Kg/sec/m^2")
  length = len("Evaporation")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_evp, 'long_name', length, "Evaporation")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_evp, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_ENDDEF(ncid)

  lon = real(x_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lon, ns_ijk(1), nc_ijk(1), lon)
  lat = real(y_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lat, ns_ijk(2), nc_ijk(2), lat)
  tmp2d = real(u_10)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_u10,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(v_10)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_v10,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(t_10)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_t10,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(q_10)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_q10,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(pslv)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_slp,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(swup)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_swu,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(swdn)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_swd,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(lwdn)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lwd,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(rain)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_ran,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(snow)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_snw,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(ifrc)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_ice,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(roff)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_rof,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(taux)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_tux,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(tauy)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_tuy,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(lwup)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lwu,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(lath)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lth,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(senh)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_snh,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp2d = real(evap)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_evp,  ns_ijk(1:2), nc_ijk(1:2), tmp2d)


  ierr = NFMPI_CLOSE(ncid)
end subroutine put_flux_pncio

subroutine put_timcom_pncio(pncfile)
  use tai_hyperlink, only: m_comm_cart, nx, ny, nz, ns_ijk, nc_ijk, &
                       nx_grid, ny_grid, nz_grid, &
                        x_grid,  y_grid,  z_grid, &
                       u2, v2, t2, s2, rho, p, p0, u, v, w
  implicit none

  character(len=*) :: pncfile
  integer(MPI_OFFSET_KIND) :: ngdim(3), length
  integer :: ncid, iid, jid, kid, nid(3), ierr
  integer :: id_lon, id_lat, id_lev, &
             id_u2, id_v2, id_t2, id_s2, id_rho, id_p, id_p0, id_u, id_v, id_w
  real(r4) :: lon(nx), lat(ny), dep(nz), tmp3d(nx,ny,nz), tmp2d(nx,ny), zero(1)

  zero = 0.d0
  ngdim(1) = nx_grid
  ngdim(2) = ny_grid
  ngdim(3) = nz_grid

  ierr = NFMPI_CREATE(m_comm_cart, trim(pncfile), NF_64BIT_OFFSET, MPI_INFO_NULL, ncid)

  ierr = NFMPI_DEF_DIM(ncid, "lon", ngdim(1), nid(1))
  ierr = NFMPI_DEF_DIM(ncid, "lat", ngdim(2), nid(2))
  ierr = NFMPI_DEF_DIM(ncid, "dep", ngdim(3), nid(3))

  ierr = NFMPI_DEF_VAR(ncid, "lon", NF_REAL, 1, nid(1), id_lon)
  length = len("degrees_east")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lon, 'units', length, "degrees_east")
  ierr = NFMPI_DEF_VAR(ncid, "lat", NF_REAL, 1, nid(2), id_lat)
  length = len("degrees_north")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lat, 'units', length, "degrees_north")
  ierr = NFMPI_DEF_VAR(ncid, "dep", NF_REAL, 1, nid(3), id_lev)
  length = len("cm")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lev, 'units', length, "cm")

  ierr = NFMPI_DEF_VAR(ncid, "U2", NF_REAL, 3, nid, id_u2)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u2, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u2, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V2", NF_REAL, 3, nid, id_v2)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v2, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v2, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "T2", NF_REAL, 3, nid, id_t2)
  length = len("C")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t2, 'units', length, "C")
  length = len("Water Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t2, 'long_name', length, "Water Temperature")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_t2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "S2", NF_REAL, 3, nid, id_s2)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s2, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s2, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_s2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "RHO", NF_REAL, 3, nid, id_rho)
  length = len("g/cm^3")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rho, 'units', length, "g/cm^3")
  length = len("density")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rho, 'long_name', length, "density")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_rho, '_FillValue', NF_REAL, length, zero)

  !ierr = NFMPI_DEF_VAR(ncid, "P", NF_REAL, 3, nid, id_p)
  !length = len("Pa")
  !ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rho, 'units', length, "Pa")
  !length = len("Pressure")
  !ierr = NFMPI_PUT_ATT_TEXT(ncid, id_rho, 'long_name', length, "Pressure")
  !length = 1
  !ierr = NFMPI_PUT_ATT_REAL(ncid, id_rho, '_FillValue', NF_REAL, length, zero)

!  ierr = NFMPI_DEF_VAR(ncid, "U", NF_REAL, 3, nid, id_u)
!  length = len("cm s-1")
!  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u, 'units', length, "cm s-1")
!  length = len("Eastward Water Velocity")
!  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u, 'long_name', length, "Eastward Water Velocity")
!  length = 1
!  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u, '_FillValue', NF_REAL, length, zero)

!  ierr = NFMPI_DEF_VAR(ncid, "V", NF_REAL, 3, nid, id_v)
!  length = len("cm s-1")
!  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v,  'units', length, "cm s-1")
!  length = len("Northward Water Velocity")
!  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v,  'long_name', length, "Northward Water Velocity")
!  length = 1
!  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v,  '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "w", NF_REAL, 3, nid, id_w)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_w, 'units', length, "cm s-1")
  length = len("Vertical Velocity at face")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_w, 'long_name', length, "Vertical Velocity at face")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_w, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "p0", NF_REAL, 2, nid(1:2), id_p0)
  length = len("Pa")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_p0, 'units', length, "Pa")
  length = len("Surface Pressure")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_p0, 'long_name', length, "Surface Pressure")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_p0, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_ENDDEF(ncid)

  lon = real(x_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lon, ns_ijk(1), nc_ijk(1), lon)
  lat = real(y_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lat, ns_ijk(2), nc_ijk(2), lat)
  dep = real(z_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lev, ns_ijk(3), nc_ijk(3), dep)
  tmp3d = real(u2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_u2,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(v2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_v2,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(t2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_t2,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(s2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_s2,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(rho(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_rho, ns_ijk,    nc_ijk, tmp3d)
  !tmp3d = real(p(1:nx,1:ny,1:nz))
  !ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_p,   ns_ijk,    nc_ijk, tmp3d)
  !tmp3d = real(u(2:nx+1,1:ny,1:nz))
  !ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_u,   ns_ijk,    nc_ijk, tmp3d)
  !tmp3d = real(v(1:nx,2:ny+1,1:nz))
  !ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_v,   ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(w(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_w,   ns_ijk,    nc_ijk, tmp3d)
  tmp2d = real(p0(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_p0,  ns_ijk(1:2),nc_ijk(1:2), tmp2d)

  ierr = NFMPI_CLOSE(ncid)
end subroutine put_timcom_pncio

subroutine put_windmix_pncio(pncfile)
  use tai_hyperlink, only: m_comm_cart, nx, ny, nz, nxf, nyf, ns_ijk, nc_ijk, &
                       nx_grid, ny_grid, nz_grid, &
                        x_grid,  y_grid,  z_grid, &
                       sw_trans, ev, hv, vdc, vvc,&
                       kpp_src, kpp_hblt, dmx, dmy
  implicit none

  character(len=*) :: pncfile
  integer(MPI_OFFSET_KIND) :: ngdim(3), length
  integer :: ncid, iid, jid, kid, nid(3), ierr
  integer :: id_lon, id_lat, id_lev, &
             id_tran, id_ev, id_hv, id_vdct, id_vdcs, id_vvc, id_hblt, &
             id_srcs, id_srct, id_dmx, id_dmy
  real(r4) :: lon(nx), lat(ny), dep(nz), tmp3d(nx,ny,nz), tmp2d(nx,ny), zero(1)

  zero = 0.d0
  ngdim(1) = nx_grid
  ngdim(2) = ny_grid
  ngdim(3) = nz_grid

  ierr = NFMPI_CREATE(m_comm_cart, trim(pncfile), NF_64BIT_OFFSET, MPI_INFO_NULL, ncid)

  ierr = NFMPI_DEF_DIM(ncid, "lon", ngdim(1), nid(1))
  ierr = NFMPI_DEF_DIM(ncid, "lat", ngdim(2), nid(2))
  ierr = NFMPI_DEF_DIM(ncid, "dep", ngdim(3), nid(3))

  ierr = NFMPI_DEF_VAR(ncid, "lon", NF_REAL, 1, nid(1), id_lon)
  length = len("degrees_east")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lon, 'units', length, "degrees_east")
  ierr = NFMPI_DEF_VAR(ncid, "lat", NF_REAL, 1, nid(2), id_lat)
  length = len("degrees_north")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lat, 'units', length, "degrees_north")
  ierr = NFMPI_DEF_VAR(ncid, "dep", NF_REAL, 1, nid(3), id_lev)
  length = len("meters")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lev, 'units', length, "meters")

  ierr = NFMPI_DEF_VAR(ncid, "transmit", NF_REAL, 3, nid, id_tran)
  length = len("[0-1]")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_tran, 'units', length, "[0-1]")
  length = len("short-wave radiation transmmit rate")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_tran, 'long_name', length, "short-wave radiation transmmit rate")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_tran, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "ev", NF_REAL, 3, nid, id_ev)
  length = len("cm/sec")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ev, 'units', length, "cm/sec")
  length = len("vertical eddy viscosity rate")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ev, 'long_name', length, "vertical eddy viscosity rate")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_ev, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "hv", NF_REAL, 3, nid, id_hv)
  length = len("cm/sec")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_hv, 'units', length, "cm/sec")
  length = len("vertical eddy diffusivity rate")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_hv, 'long_name', length, "vertical eddy diffusivity rate")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_hv, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "vdc_t", NF_REAL, 3, nid, id_vdct)
  length = len("cm^2/sec")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_vdct, 'units', length, "PSU")
  length = len("vertical eddy diffusivity for temp")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_vdct, 'long_name', length, "vertical eddy diffusivity for temp")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_vdct, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "vdc_s", NF_REAL, 3, nid, id_vdcs)
  length = len("cm^2/sec")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_vdcs, 'units', length, "g/cm^3")
  length = len("vertical eddy diffusivity for salt")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_vdcs, 'long_name', length, "vertical eddy diffusivity for salt")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_vdcs, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "vvc", NF_REAL, 3, nid, id_vvc)
  length = len("cm^2/sec")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_vvc, 'units', length, "cm^2/sec")
  length = len("vertical eddy viscosity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_vvc, 'long_name', length, "vertical eddy viscosity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_vvc, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "kpp_srcs", NF_REAL, 3, nid, id_srcs)
  length = len("1d-5*psu/sec")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_srcs, 'units', length, "1d-5*psu/sec")
  length = len("kpp_src salt")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_srcs, 'long_name', length, "kpp_src salt")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_srcs, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "kpp_srct", NF_REAL, 3, nid, id_srct)
  length = len("1d-5*degC/sec")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_srct, 'units', length, "1d-5*degC/sec")
  length = len("kpp_src temp")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_srct, 'long_name', length, "kpp_src temp")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_srct, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "kpp_hblt", NF_REAL, 2, nid(1:2), id_hblt)
  length = len("cm")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_hblt, 'units', length, "cm")
  length = len("kpp boundary layer depth")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_hblt, 'long_name', length, "kpp boundary layer depth")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_hblt, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "dmx", NF_REAL, 3, nid, id_dmx)
  length = len("cm/sec")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_dmx, 'units', length, "cm/sec")
  length = len("horizental eddy viscosity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_dmx, 'long_name', length, "horizental eddy viscosity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_dmx, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "dmy", NF_REAL, 3, nid, id_dmy)
  length = len("cm/sec")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_dmy, 'units', length, "cm/sec")
  length = len("horizental eddy viscosity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_dmy, 'long_name', length, "horizental eddy viscosity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_dmy, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_ENDDEF(ncid)

  lon = real(x_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lon, ns_ijk(1), nc_ijk(1), lon)
  lat = real(y_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lat, ns_ijk(2), nc_ijk(2), lat)
  dep = real(z_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lev, ns_ijk(3), nc_ijk(3), dep)
  tmp3d = real(sw_trans)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_tran, ns_ijk, nc_ijk, tmp3d)
  tmp3d = real(ev)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_ev,   ns_ijk, nc_ijk, tmp3d)
  tmp3d = real(hv(:,:,:,1))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_hv,   ns_ijk, nc_ijk, tmp3d)
  tmp3d = real(vdc(:,:,1:nz,1))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_vdct, ns_ijk, nc_ijk, tmp3d)
  tmp3d = real(vdc(:,:,1:nz,2))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_vdcs, ns_ijk, nc_ijk, tmp3d)
  tmp3d = real(vvc)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_vvc,  ns_ijk, nc_ijk, tmp3d)
  tmp3d = real(kpp_src(:,:,:,2)*1.d5)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_srcs, ns_ijk, nc_ijk, tmp3d)
  tmp3d = real(kpp_src(:,:,:,1)*1.d5)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_srct, ns_ijk, nc_ijk, tmp3d)
  tmp2d = real(kpp_hblt)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_hblt, ns_ijk(1:2), nc_ijk(1:2), tmp2d)
  tmp3d = real(dmx(2:nxf,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_dmx,  ns_ijk, nc_ijk, tmp3d)
  tmp3d = real(dmy(1:nx,2:nyf,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_dmy,  ns_ijk, nc_ijk, tmp3d)

  ierr = NFMPI_CLOSE(ncid)
end subroutine put_windmix_pncio

subroutine get_restart_pncio(pncfile)
  use tai_hyperlink, only: m_comm_cart, ns_ijk, nc_ijk, &
                       nc_ujk, nc_ivk, nc_ijw,      &
                       nx, ny, nz, nxf, nyf, nzf,   &
                       nx_grid, ny_grid, nz_grid,   &
                        x_grid,  y_grid,  z_grid,   &
                        x_face,  y_face,  z_face,   &
                       u1, u2, ulf, v1, v2, vlf,    &
                       t1, t2, tlf, s1, s2, slf,    &
                       u, v, w, p0, x, cgr, cgrh,   &
                       cgp, cgv, cgs, cgt
  implicit none

  character(len=*) :: pncfile
  integer :: ierr, ncid, varid

  ierr = NFMPI_OPEN(m_comm_cart, trim(pncfile), NF_NOWRITE, MPI_INFO_NULL, ncid)
  ierr = NFMPI_INQ_VARID(ncid, 'u1', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, u1(1:nx,1:ny,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 'u2', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, u2(1:nx,1:ny,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 'ulf', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, ulf(1:nx,1:ny,1:nz))

  ierr = NFMPI_INQ_VARID(ncid, 'v1', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, v1(1:nx,1:ny,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 'v2', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, v2(1:nx,1:ny,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 'vlf', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, vlf(1:nx,1:ny,1:nz))

  ierr = NFMPI_INQ_VARID(ncid, 't1', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, t1(1:nx,1:ny,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 't2', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, t2(1:nx,1:ny,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 'tlf', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, tlf(1:nx,1:ny,1:nz))

  ierr = NFMPI_INQ_VARID(ncid, 's1', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, s1(1:nx,1:ny,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 's2', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, s2(1:nx,1:ny,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 'slf', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijk, slf(1:nx,1:ny,1:nz))

  ierr = NFMPI_INQ_VARID(ncid, 'u', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ujk, u(1:nxf,1:ny,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 'v', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ivk, v(1:nx,1:nyf,1:nz))
  ierr = NFMPI_INQ_VARID(ncid, 'w', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk, nc_ijw, w(1:nx,1:ny,1:nzf))

  ierr = NFMPI_INQ_VARID(ncid, 'p0', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), p0(1:nx,1:ny))
  ierr = NFMPI_INQ_VARID(ncid, 'x', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), x(1:nx,1:ny))
  ierr = NFMPI_INQ_VARID(ncid, 'cgr', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), cgr(1:nx,1:ny))
  ierr = NFMPI_INQ_VARID(ncid, 'cgrh', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), cgrh(1:nx,1:ny))
  ierr = NFMPI_INQ_VARID(ncid, 'cgp', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), cgp(1:nx,1:ny))
  ierr = NFMPI_INQ_VARID(ncid, 'cgv', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), cgv(1:nx,1:ny))
  ierr = NFMPI_INQ_VARID(ncid, 'cgs', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), cgs(1:nx,1:ny))
  ierr = NFMPI_INQ_VARID(ncid, 'cgt', varid)
  ierr = NFMPI_GET_VARA_DOUBLE_ALL(ncid, varid, ns_ijk(1:2), nc_ijk(1:2), cgt(1:nx,1:ny))

  ierr = NFMPI_CLOSE(ncid)

end subroutine get_restart_pncio

subroutine put_restart_pncio(pncfile)
  use tai_hyperlink, only: m_comm_cart, ns_ijk, nc_ijk, &
                       nc_ujk, nc_ivk, nc_ijw,      &
                       nx, ny, nz, nxf, nyf, nzf,   &
                       nx_grid, ny_grid, nz_grid,   &
                        x_grid,  y_grid,  z_grid,   &
                        x_face,  y_face,  z_face,   &
                       u1, u2, ulf, v1, v2, vlf,    &
                       t1, t2, tlf, s1, s2, slf,    &
                       u, v, w, p0, x, cgr, cgrh,   &
                       cgp, cgv, cgs, cgt
  implicit none

  character(len=*) :: pncfile
  integer(MPI_OFFSET_KIND) :: ngdim(6), length
  integer :: ncid, iid, jid, kid, nid(6), ierr
  integer :: id_lon, id_lat, id_lev, id_lonf, id_latf, id_levf, &
             id_u2, id_v2, id_t2, id_s2, id_u, id_v, id_w, id_p0, &
             id_u1, id_v1, id_t1, id_s1, id_ulf, id_vlf, id_tlf, id_slf, &
             id_x, id_cgr, id_cgrh, id_cgp, id_cgv, id_cgs, id_cgt
  real(r4) :: lon(nx), lat(ny), dep(nz), tmp3d(nx,ny,nz), tmp2d(nx,ny), zero(1)

  zero = 0.d0
  ngdim(1) = nx_grid
  ngdim(2) = ny_grid
  ngdim(3) = nz_grid
  ngdim(4) = nx_grid+1
  ngdim(5) = ny_grid+1
  ngdim(6) = nz_grid+1

  ierr = NFMPI_CREATE(m_comm_cart, trim(pncfile), NF_64BIT_OFFSET, MPI_INFO_NULL, ncid)

  ierr = NFMPI_DEF_DIM(ncid, "lon_c", ngdim(1), nid(1))
  ierr = NFMPI_DEF_DIM(ncid, "lat_c", ngdim(2), nid(2))
  ierr = NFMPI_DEF_DIM(ncid, "dep_c", ngdim(3), nid(3))
  ierr = NFMPI_DEF_DIM(ncid, "lon_f", ngdim(4), nid(4))
  ierr = NFMPI_DEF_DIM(ncid, "lat_f", ngdim(5), nid(5))
  ierr = NFMPI_DEF_DIM(ncid, "dep_f", ngdim(6), nid(6))

  ierr = NFMPI_DEF_VAR(ncid, "lon_c", NF_DOUBLE, 1, nid(1), id_lon)
  length = len("degrees_east")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lon, 'units', length, "degrees_east")
  ierr = NFMPI_DEF_VAR(ncid, "lat_c", NF_DOUBLE, 1, nid(2), id_lat)
  length = len("degrees_north")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lat, 'units', length, "degrees_north")
  ierr = NFMPI_DEF_VAR(ncid, "dep_c", NF_DOUBLE, 1, nid(3), id_lev)
  length = len("cm")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lev, 'units', length, "cm")
  ierr = NFMPI_DEF_VAR(ncid, "lon_f", NF_DOUBLE, 1, nid(4), id_lon)
  length = len("degrees_east")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lon, 'units', length, "degrees_east")
  ierr = NFMPI_DEF_VAR(ncid, "lat_f", NF_DOUBLE, 1, nid(5), id_lat)
  length = len("degrees_north")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lat, 'units', length, "degrees_north")
  ierr = NFMPI_DEF_VAR(ncid, "dep_f", NF_DOUBLE, 1, nid(6), id_lev)
  length = len("cm")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lev, 'units', length, "cm")

  ierr = NFMPI_DEF_VAR(ncid, "u2", NF_DOUBLE, 3, nid(1:3), id_u2)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u2, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u2, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "u1", NF_DOUBLE, 3, nid(1:3), id_u1)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u1, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u1, 'long_name', length, "Eastward Water Velocity")
  length = 1            
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u1, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "ulf", NF_DOUBLE, 3, nid(1:3), id_ulf)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ulf, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ulf, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_ulf, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "v2", NF_DOUBLE, 3, nid(1:3), id_v2)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v2, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v2, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "v1", NF_DOUBLE, 3, nid(1:3), id_v1)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v1, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v1, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v1, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "vlf", NF_DOUBLE, 3, nid(1:3), id_vlf)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_vlf, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_vlf, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_vlf, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "t2", NF_DOUBLE, 3, nid(1:3), id_t2)
  length = len("C")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t2, 'units', length, "C")
  length = len("Water Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t2, 'long_name', length, "Water Temperature")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_t2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "t1", NF_DOUBLE, 3, nid(1:3), id_t1)
  length = len("C")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t1, 'units', length, "C")
  length = len("Water Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t1, 'long_name', length, "Water Temperature")
  length = 1           
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_t1, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "tlf", NF_DOUBLE, 3, nid(1:3), id_tlf)
  length = len("C")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_tlf, 'units', length, "C")
  length = len("Water Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_tlf, 'long_name', length, "Water Temperature")
  length = 1           
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_tlf, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "s2", NF_DOUBLE, 3, nid(1:3), id_s2)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s2, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s2, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_s2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "s1", NF_DOUBLE, 3, nid(1:3), id_s1)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s1, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s1, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_s1, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "slf", NF_DOUBLE, 3, nid(1:3), id_slf)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_slf, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_slf, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_slf, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "u", NF_DOUBLE, 3, (/nid(4),nid(2),nid(3)/), id_u)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "v", NF_DOUBLE, 3, (/nid(1),nid(5),nid(3)/), id_v)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v,  'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v,  'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v,  '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "w", NF_DOUBLE, 3, (/nid(1),nid(2),nid(6)/), id_w)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_w, 'units', length, "cm s-1")
  length = len("Vertical Velocity at face")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_w, 'long_name', length, "Vertical Velocity at face")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_w, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "p0", NF_DOUBLE, 2, nid(1:2), id_p0)
  length = len("Pa")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_p0, 'units', length, "Pa")
  length = len("Surface Pressure")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_p0, 'long_name', length, "Surface Pressure")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_p0, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "x", NF_DOUBLE, 2, nid(1:2), id_x)
  length = len("Pressure solver residual x")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_x, 'long_name', length, "Pressure solver residual x")

  ierr = NFMPI_DEF_VAR(ncid, "cgr", NF_DOUBLE, 2, nid(1:2), id_cgr)
  length = len("Pressure solver cgr")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_cgr, 'lng_name', length, "Pressure solver cgr")

  ierr = NFMPI_DEF_VAR(ncid, "cgrh", NF_DOUBLE, 2, nid(1:2), id_cgrh)
  length = len("Pressure solver cgr")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_cgrh, 'lng_name', length, "Pressure solver cgrh")

  ierr = NFMPI_DEF_VAR(ncid, "cgp", NF_DOUBLE, 2, nid(1:2), id_cgp)
  length = len("Pressure solver cgp")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_cgp, 'lng_name', length, "Pressure solver cgp")

  ierr = NFMPI_DEF_VAR(ncid, "cgv", NF_DOUBLE, 2, nid(1:2), id_cgv)
  length = len("Pressure solver cgv")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_cgv, 'lng_name', length, "Pressure solver cgv")

  ierr = NFMPI_DEF_VAR(ncid, "cgs", NF_DOUBLE, 2, nid(1:2), id_cgs)
  length = len("Pressure solver cgs")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_cgs, 'lng_name', length, "Pressure solver cgs")

  ierr = NFMPI_DEF_VAR(ncid, "cgt", NF_DOUBLE, 2, nid(1:2), id_cgt)
  length = len("Pressure solver cgt")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_cgt, 'lng_name', length, "Pressure solver cgt")

  ierr = NFMPI_ENDDEF(ncid)

  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_lon, ns_ijk(1), nc_ijk(1), x_grid)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_lat, ns_ijk(2), nc_ijk(2), y_grid)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_lev, ns_ijk(3), nc_ijk(3), z_grid)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_lon, ns_ijk(1), (/nc_ijk(1)+1/), x_face)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_lat, ns_ijk(2), (/nc_ijk(2)+1/), y_face)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_lev, ns_ijk(3), (/nc_ijk(3)+1/), z_face)

  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_u2,  ns_ijk, nc_ijk, u2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_u1,  ns_ijk, nc_ijk, u1(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_ulf, ns_ijk, nc_ijk, ulf(1:nx,1:ny,1:nz))

  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_v2,  ns_ijk, nc_ijk, v2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_v1,  ns_ijk, nc_ijk, v1(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_vlf, ns_ijk, nc_ijk, vlf(1:nx,1:ny,1:nz))

  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_t2,  ns_ijk, nc_ijk, t2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_t1,  ns_ijk, nc_ijk, t1(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_tlf, ns_ijk, nc_ijk, tlf(1:nx,1:ny,1:nz))

  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_s2,  ns_ijk, nc_ijk, s2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_s1,  ns_ijk, nc_ijk, s1(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_slf, ns_ijk, nc_ijk, slf(1:nx,1:ny,1:nz))

  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_u,   ns_ijk, nc_ujk, u(1:nxf,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_v,   ns_ijk, nc_ivk, v(1:nx,1:nyf,1:nz))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_w,   ns_ijk, nc_ijw, w(1:nx,1:ny,1:nzf))

  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_p0,  ns_ijk(1:2), nc_ijk(1:2),  p0(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_x,   ns_ijk(1:2), nc_ijk(1:2),   x(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_cgr, ns_ijk(1:2), nc_ijk(1:2), cgr(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_cgrh,ns_ijk(1:2), nc_ijk(1:2),cgrh(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_cgp, ns_ijk(1:2), nc_ijk(1:2), cgp(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_cgv, ns_ijk(1:2), nc_ijk(1:2), cgv(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_cgs, ns_ijk(1:2), nc_ijk(1:2), cgs(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_cgt, ns_ijk(1:2), nc_ijk(1:2), cgt(1:nx,1:ny))

  ierr = NFMPI_CLOSE(ncid)
end subroutine put_restart_pncio

subroutine put_blowup_pncio(pncfile)
  use tai_hyperlink, only: m_comm_cart, nx, ny, nz, ns_ijk, nc_ijk, &
                       nx_grid, ny_grid, nz_grid, &
                        x_grid,  y_grid,  z_grid, &
                       u2, v2, t2, s2, u, v, w, p0
  implicit none

  character(len=*) :: pncfile
  integer(MPI_OFFSET_KIND) :: ngdim(3), length
  integer :: ncid, iid, jid, kid, nid(3), ierr
  integer :: id_lon, id_lat, id_lev, &
             id_u2, id_v2, id_t2, id_s2, id_u, id_v, id_w, id_p0
  real(r4) :: lon(nx), lat(ny), dep(nz), tmp3d(nx,ny,nz), tmp2d(nx,ny), zero(1)

  zero = 0.d0
  ngdim(1) = nx_grid
  ngdim(2) = ny_grid
  ngdim(3) = nz_grid

  ierr = NFMPI_CREATE(m_comm_cart, trim(pncfile), NF_64BIT_OFFSET, MPI_INFO_NULL, ncid)

  ierr = NFMPI_DEF_DIM(ncid, "lon", ngdim(1), nid(1))
  ierr = NFMPI_DEF_DIM(ncid, "lat", ngdim(2), nid(2))
  ierr = NFMPI_DEF_DIM(ncid, "dep", ngdim(3), nid(3))

  ierr = NFMPI_DEF_VAR(ncid, "lon", NF_DOUBLE, 1, nid(1), id_lon)
  length = len("degrees_east")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lon, 'units', length, "degrees_east")
  ierr = NFMPI_DEF_VAR(ncid, "lat", NF_REAL, 1, nid(2), id_lat)
  length = len("degrees_north")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lat, 'units', length, "degrees_north")
  ierr = NFMPI_DEF_VAR(ncid, "dep", NF_REAL, 1, nid(3), id_lev)
  length = len("cm")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lev, 'units', length, "cm")

  ierr = NFMPI_DEF_VAR(ncid, "U2", NF_REAL, 3, nid, id_u2)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u2, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u2, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V2", NF_REAL, 3, nid, id_v2)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v2, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v2, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "T2", NF_REAL, 3, nid, id_t2)
  length = len("C")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t2, 'units', length, "C")
  length = len("Water Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_t2, 'long_name', length, "Water Temperature")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_t2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "S2", NF_REAL, 3, nid, id_s2)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s2, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_s2, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_s2, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "U", NF_REAL, 3, nid, id_u)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_u, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_u, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V", NF_REAL, 3, nid, id_v)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v,  'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_v,  'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_v,  '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "w", NF_REAL, 3, nid, id_w)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_w, 'units', length, "cm s-1")
  length = len("Vertical Velocity at face")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_w, 'long_name', length, "Vertical Velocity at face")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_w, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "p0", NF_REAL, 2, nid(1:2), id_p0)
  length = len("Pa")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_p0, 'units', length, "Pa")
  length = len("Surface Pressure")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_p0, 'long_name', length, "Surface Pressure")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_p0, '_FillValue', NF_REAL, length, zero)

  ierr = NFMPI_ENDDEF(ncid)

  lon = real(x_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lon, ns_ijk(1), nc_ijk(1), lon)
  lat = real(y_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lat, ns_ijk(2), nc_ijk(2), lat)
  dep = real(z_grid)
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_lev, ns_ijk(3), nc_ijk(3), dep)
  tmp3d = real(u2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_u2,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(v2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_v2,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(t2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_t2,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(s2(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_s2,  ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(u(2:nx+1,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_u,   ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(v(1:nx,2:ny+1,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_v,   ns_ijk,    nc_ijk, tmp3d)
  tmp3d = real(w(1:nx,1:ny,1:nz))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_w ,  ns_ijk,    nc_ijk, tmp3d)
  tmp2d = real(p0(1:nx,1:ny))
  ierr = NFMPI_PUT_VARA_REAL_ALL(ncid, id_p0,  ns_ijk(1:2),nc_ijk(1:2), tmp2d)

  ierr = NFMPI_CLOSE(ncid)

end subroutine put_blowup_pncio

subroutine get_ts_climate(nx, ny, nz, t_clim, s_clim)
  use tai_hyperlink, only: m_comm_cart, ns_ijk, nc_ijk
  implicit none

  integer, intent(in) :: nx, ny, nz
  real(r8), dimension(nx,ny,12), intent(out) :: t_clim, s_clim

  character(len=256) :: pncfile
  integer :: ncid, ierr, varid
  real(r4) :: var4d(nx,ny,nz,12)
  integer(MPI_OFFSET_KIND) :: ns_ijkm(4), nc_ijkm(4)

  pncfile = "/data/common/gfs/GEPSv3_lib/data/TIMCOM_docn_ts_tai.nc"

  ns_ijkm(1:3) = ns_ijk
  nc_ijkm(1:3) = nc_ijk

  ns_ijkm(4) = 1
  nc_ijkm(4) = 12

  ierr = NFMPI_OPEN(m_comm_cart, trim(pncfile), NF_NOWRITE, MPI_INFO_NULL, ncid)
  ierr = NFMPI_INQ_VARID(ncid, 'tclim', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_ijkm, nc_ijkm, var4d)
  t_clim(1:nx,1:ny,1:12) = dble(var4d(1:nx,1:ny,1,1:12))
  ierr = NFMPI_INQ_VARID(ncid, 'sclim', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_ijkm, nc_ijkm, var4d)
  s_clim(1:nx,1:ny,1:12) = dble(var4d(1:nx,1:ny,1,1:12))
  ierr = NFMPI_CLOSE(ncid)
end subroutine get_ts_climate

subroutine get_lateral_bc_info(jwbc, jebc, isbc, inbc, ns_lbc, nc_lbc)
  use tai_hyperlink, only: m_comm_cart, nx, ny, nz, myid_x, myid_y, &
                       ns_ijk, nc_ijk, x_grid, y_grid, z_grid, npx
  implicit none

  integer, intent(out) :: jwbc(ny,nz), jebc(ny,nz), isbc(nx,nz), inbc(nx,nz)
  integer(MPI_OFFSET_KIND), intent(out) :: ns_lbc(3), nc_lbc(3)
  integer :: i, j, k, e, ierr, ncid, dimid, varid, xbnd(2), ybnd(2)
  integer(MPI_OFFSET_KIND) :: nds(3), ns_k(2), nc_k(2)
  character(len=6) :: fdims(3)
  character(len=256) :: pncfile
  real(r8), allocatable :: xgrd2(:), ygrd2(:), zgrd2(:)

  pncfile = "/data/common/gfs/GEPSv3_lib/data/timcom_grid_1536x720x55.nc"
  fdims(1) = 'x_face'
  fdims(2) = 'y_face'
  fdims(3) = 'z_grid'

  ierr = NFMPI_OPEN(m_comm_cart, trim(pncfile), NF_NOWRITE, MPI_INFO_NULL, ncid)
  do i = 1, 3
    ierr = NFMPI_INQ_DIMID(ncid, trim(fdims(i)), dimid)
    ierr = NFMPI_INQ_DIMLEN(ncid, dimid, nds(i))
  end do

  allocate(xgrd2(nds(1)))
  allocate(ygrd2(nds(2)))
  allocate(zgrd2(nds(3)))

  ierr = NFMPI_INQ_VARID(ncid, 'x_face', varid)
  ierr = NFMPI_GET_VAR_ALL(ncid, varid, xgrd2, nds(1), MPI_REAL8)
  ierr = NFMPI_INQ_VARID(ncid, 'y_face', varid)
  ierr = NFMPI_GET_VAR_ALL(ncid, varid, ygrd2, nds(2), MPI_REAL8)
  ierr = NFMPI_INQ_VARID(ncid, 'z_grid', varid)
  ierr = NFMPI_GET_VAR_ALL(ncid, varid, zgrd2, nds(3), MPI_REAL8)

  ierr = NFMPI_CLOSE(ncid)

  xbnd = 0
  ybnd = 0

  e = nx
  do i = 2, nds(1)
    if(xgrd2(i-1).le.x_grid(1) .and. xgrd2(i).gt.x_grid(1)) then
      xbnd(1) = i-1
    end if
    if(xgrd2(i-1).le.x_grid(e) .and. xgrd2(i).gt.x_grid(e)) then
      xbnd(2) = i-1
      exit
    end if
  end do

  e = ny
  do j = 2, nds(2)
    if(ygrd2(j-1).le.y_grid(1) .and. ygrd2(j).gt.y_grid(1)) then
      ybnd(1) = j-1
    end if
    if(ygrd2(j-1).le.y_grid(e) .and. ygrd2(j).gt.y_grid(e)) then
      ybnd(2) = j-1
      exit
    end if
  end do

  ns_lbc = 0
  nc_lbc = 0
  ns_lbc = (/xbnd(1), ybnd(1), 1/)
  nc_lbc = (/xbnd(2)-xbnd(1)+1, ybnd(2)-ybnd(1)+1, nz/)

  pncfile = "tai_lbc.nc"
  ierr = NFMPI_OPEN(m_comm_cart, trim(pncfile), NF_NOWRITE, MPI_INFO_NULL, ncid)

  ns_k = (/ns_ijk(1), ns_ijk(3)/)
  nc_k = (/nc_ijk(1), nc_ijk(3)/)

  ierr = NFMPI_INQ_VARID(ncid, 'sbc', varid)
  ierr = NFMPI_GET_VARA_INT_ALL(ncid, varid, ns_k, nc_k, isbc)
  ierr = NFMPI_INQ_VARID(ncid, 'nbc', varid)
  ierr = NFMPI_GET_VARA_INT_ALL(ncid, varid, ns_k, nc_k, inbc)

  ns_k = (/ns_ijk(2), ns_ijk(3)/)
  nc_k = (/nc_ijk(2), nc_ijk(3)/)

  ierr = NFMPI_INQ_VARID(ncid, 'wbc', varid)
  ierr = NFMPI_GET_VARA_INT_ALL(ncid, varid, ns_k, nc_k, jwbc)
  ierr = NFMPI_INQ_VARID(ncid, 'ebc', varid)
  ierr = NFMPI_GET_VARA_INT_ALL(ncid, varid, ns_k, nc_k, jebc)

  ierr = NFMPI_CLOSE(ncid)

  !if(myid_x .eq. npx-1 .and. myid_y .eq. 9) then
  !  write(*,"(a,10i4)") "[CHECK EBC0: ", ns_lbc, nc_lbc, xbnd, ybnd
  !  write(*,*) "[CHECK GRD]:", xgrd2(1:5), x_grid(1:5), ygrd2(1:5), y_grid(1:5)
  !  write(*,*) "[CHECK JEBC0: ", jebc
  !end if
  !if(isbc(1,1) .lt. ns_lbc(1))  write(*,'(9i5)') '[ERROR] check xbnd: ', myid_x, myid_y, xbnd, isbc(1,1), isbc(nx,1), ns_lbc
  !if(jebc(1,1) .lt. ns_lbc(2))  write(*,'(9i5)') '[ERROR] check ybnd: ', myid_x, myid_y, ybnd, jebc(1,1), jebc(ny,1), ns_lbc
  !if(nc_lbc(1) .lt. 1) write(*,"(10i5)") "[ERROR] check xlbc: ", myid_x, myid_y, xbnd, ns_lbc, nc_lbc
  !if(nc_lbc(2) .lt. 1) write(*,"(10i5)") "[ERROR] check ylbc: ", myid_x, myid_y, ybnd, ns_lbc, nc_lbc

  isbc = isbc - ns_lbc(1) + 1
  inbc = inbc - ns_lbc(1) + 1
  jwbc = jwbc - ns_lbc(2) + 1
  jebc = jebc - ns_lbc(2) + 1

  deallocate(xgrd2)
  deallocate(ygrd2)
  deallocate(zgrd2)
end subroutine get_lateral_bc_info

subroutine get_lateral_bc_field(pncfile, ns_lbc, nc_lbc, &
                                jwbc, jebc, isbc, inbc,  &
                                wbc, ebc, sbc, nbc)
  use tai_hyperlink, only: m_comm_cart, myid_x, myid_y, nx, ny, nz, in, npx, npy
  implicit none

  character(len=*), intent(in) :: pncfile
  integer(MPI_OFFSET_KIND) :: ns_lbc(3), nc_lbc(3)
  integer, intent(in) :: jwbc(ny,nz), jebc(ny,nz), isbc(nx,nz), inbc(nx,nz)
  real(r8), intent(out) :: wbc(ny,nz,4), ebc(ny,nz,4), sbc(nx,nz,4), nbc(nx,nz,4)
  integer :: i, j, k, ierr, ncid, varid
  real(r4), allocatable :: var3d(:,:,:,:)

  allocate(var3d(nc_lbc(1),nc_lbc(2),nc_lbc(3),4))

  var3d = 0.0
  ierr = NFMPI_OPEN(m_comm_cart, trim(pncfile), NF_NOWRITE, MPI_INFO_NULL, ncid)
  ierr = NFMPI_INQ_VARID(ncid, 'U2', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_lbc, nc_lbc, var3d(:,:,:,1))
  ierr = NFMPI_INQ_VARID(ncid, 'V2', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_lbc, nc_lbc, var3d(:,:,:,2))
  ierr = NFMPI_INQ_VARID(ncid, 'T2', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_lbc, nc_lbc, var3d(:,:,:,3))
  ierr = NFMPI_INQ_VARID(ncid, 'S2', varid)
  ierr = NFMPI_GET_VARA_REAL_ALL(ncid, varid, ns_lbc, nc_lbc, var3d(:,:,:,4))
  ierr = NFMPI_CLOSE(ncid)

!  if(myid_x .eq. npx-1 .and. myid_y .eq. 9) then
!    write(*,"(a,6i4,6D12.3)") "[CHECK EBC]: ", ns_lbc, nc_lbc, var3d(:,:,:,3) 
!    write(*,*) "[CHECK JEBC]: ", jebc
!  end if

  sbc = 0.d0
  nbc = 0.d0
  wbc = 0.d0
  ebc = 0.d0

  do k = 1, nz
    do i = 1, nx
      if(isbc(i,k).gt.0) sbc(i,k,:) = dble(var3d(isbc(i,k),1,k,:))
      if(inbc(i,k).gt.0) nbc(i,k,:) = dble(var3d(inbc(i,k),nc_lbc(2),k,:))
    end do
    do j = 1, ny
      if(jwbc(j,k).gt.0) wbc(j,k,:) = dble(var3d(1,jwbc(j,k),k,:))
      if(jebc(j,k).gt.0) ebc(j,k,:) = dble(var3d(nc_lbc(1),jebc(j,k),k,:))
    end do
  end do

  do k = 2, nz
    do i = 1, nx
      if(in(i, 1,k).eq.1 .and. sbc(i,k,4).lt.1.d0) sbc(i,k,:) = sbc(i,k-1,:)
      if(in(i,ny,k).eq.1 .and. nbc(i,k,4).lt.1.d0) nbc(i,k,:) = nbc(i,k-1,:)
    end do
    do j = 1, ny
      if(in( 1,j,k).eq.1 .and. wbc(j,k,4).lt.1.d0) wbc(j,k,:) = wbc(j,k-1,:)
      if(in(nx,j,k).eq.1 .and. ebc(j,k,4).lt.1.d0) ebc(j,k,:) = ebc(j,k-1,:)
    end do
  end do

  deallocate(var3d)

end subroutine get_lateral_bc_field

subroutine check_lateral_bc(field_file, wbc, ebc, sbc, nbc)
  use tai_hyperlink, only: m_comm_cart, myid_x, myid_y, npx, npy, nx_grid, ny_grid, &
                       x_grid, y_grid, z_grid,  nx, ny, nz, ns_ijk, nc_ijk
  implicit none

  character(len=*) :: field_file
  real(r8) :: wbc(ny,nz,4), ebc(ny,nz,4), sbc(nx,nz,4), nbc(nx,nz,4)

  integer(MPI_OFFSET_KIND) :: length, ns_bc(2), nc_bc(2), ngdim(3)
  integer :: ierr, nid(3), nid_we(2), nid_sn(2), ncid, id_lon, id_lat, id_lev, &
             id_wbc_u, id_wbc_v, id_wbc_t, id_wbc_s, id_ebc_u, id_ebc_v, id_ebc_t, id_ebc_s, &
             id_sbc_u, id_sbc_v, id_sbc_t, id_sbc_s, id_nbc_u, id_nbc_v, id_nbc_t, id_nbc_s
  real(r4) :: zero(1)

  zero = 0.d0

  ngdim(1) = nx_grid
  ngdim(2) = ny_grid
  ngdim(3) = nz

  ierr = NFMPI_CREATE(m_comm_cart, trim(field_file), NF_64BIT_OFFSET, MPI_INFO_NULL, ncid)

  ierr = NFMPI_DEF_DIM(ncid, "lon", ngdim(1), nid(1))
  ierr = NFMPI_DEF_DIM(ncid, "lat", ngdim(2), nid(2))
  ierr = NFMPI_DEF_DIM(ncid, "dep", ngdim(3), nid(3))

  ierr = NFMPI_DEF_VAR(ncid, "lon", NF_DOUBLE, 1, nid(1), id_lon)
  length = len("degrees_east")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lon, 'units', length, "degrees_east")
  ierr = NFMPI_DEF_VAR(ncid, "lat", NF_DOUBLE, 1, nid(2), id_lat)
  length = len("degrees_north")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lat, 'units', length, "degrees_north")
  ierr = NFMPI_DEF_VAR(ncid, "dep", NF_DOUBLE, 1, nid(3), id_lev)
  length = len("centimeter")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_lev, 'units', length, "centimeter")

  nid_we(1) = nid(2)
  nid_we(2) = nid(3)
  ierr = NFMPI_DEF_VAR(ncid, "U1_BC_W", NF_DOUBLE, 2, nid_we, id_wbc_u)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_wbc_u, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_wbc_u, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_wbc_u, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V1_BC_W", NF_DOUBLE, 2, nid_we, id_wbc_v)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_wbc_v, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_wbc_v, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_wbc_v, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "T1_BC_W", NF_DOUBLE, 2, nid_we, id_wbc_t)
  length = len("K")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_wbc_t, 'units', length, "K")
  length = len("Potential Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_wbc_t, 'long_name', length, "Potential Temperature")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_wbc_t, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "S1_BC_W", NF_DOUBLE, 2, nid_we, id_wbc_s)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_wbc_s, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_wbc_s, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_wbc_s, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "U1_BC_E", NF_DOUBLE, 2, nid_we, id_ebc_u)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ebc_u, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ebc_u, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_ebc_u, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V1_BC_E", NF_DOUBLE, 2, nid_we, id_ebc_v)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ebc_v, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ebc_v, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_ebc_v, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "T1_BC_E", NF_DOUBLE, 2, nid_we, id_ebc_t)
  length = len("K")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ebc_t, 'units', length, "K")
  length = len("Potential Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ebc_t, 'long_name', length, "Potential Temperature")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_ebc_t, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "S1_BC_E", NF_DOUBLE, 2, nid_we, id_ebc_s)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ebc_s, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_ebc_s, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_ebc_s, '_FillValue', NF_DOUBLE, length, zero)

  nid_sn(1) = nid(1)
  nid_sn(2) = nid(3)

  ierr = NFMPI_DEF_VAR(ncid, "U1_BC_S", NF_DOUBLE, 2, nid_sn, id_sbc_u)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_sbc_u, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_sbc_u, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_sbc_u, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V1_BC_S", NF_DOUBLE, 2, nid_sn, id_sbc_v)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_sbc_v, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_sbc_v, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_sbc_v, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "T1_BC_S", NF_DOUBLE, 2, nid_sn, id_sbc_t)
  length = len("K")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_sbc_t, 'units', length, "K")
  length = len("Potential Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_sbc_t, 'long_name', length, "Potential Temperature")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_sbc_t, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "S1_BC_S", NF_DOUBLE, 2, nid_sn, id_sbc_s)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_sbc_s, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_sbc_s, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_sbc_s, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "U1_BC_N", NF_DOUBLE, 2, nid_sn, id_nbc_u)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_nbc_u, 'units', length, "cm s-1")
  length = len("Eastward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_nbc_u, 'long_name', length, "Eastward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_nbc_u, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "V1_BC_N", NF_DOUBLE, 2, nid_sn, id_nbc_v)
  length = len("cm s-1")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_nbc_v, 'units', length, "cm s-1")
  length = len("Northward Water Velocity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_nbc_v, 'long_name', length, "Northward Water Velocity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_nbc_v, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "T1_BC_N", NF_DOUBLE, 2, nid_sn, id_nbc_t)
  length = len("K")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_nbc_t, 'units', length, "K")
  length = len("Potential Temperature")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_nbc_t, 'long_name', length, "Potential Temperature")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_nbc_t, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_DEF_VAR(ncid, "S1_BC_N", NF_DOUBLE, 2, nid_sn, id_nbc_s)
  length = len("PSU")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_nbc_s, 'units', length, "PSU")
  length = len("Salinity")
  ierr = NFMPI_PUT_ATT_TEXT(ncid, id_nbc_s, 'long_name', length, "Salinity")
  length = 1
  ierr = NFMPI_PUT_ATT_REAL(ncid, id_nbc_s, '_FillValue', NF_DOUBLE, length, zero)

  ierr = NFMPI_ENDDEF(ncid)

  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_lon, ns_ijk(1), nc_ijk(1), x_grid)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_lat, ns_ijk(2), nc_ijk(2), y_grid)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_lev, ns_ijk(3), nc_ijk(3), z_grid)

  ns_bc = (/ns_ijk(2), ns_ijk(3)/)
  nc_bc = 0
  if(myid_x .eq. 0) nc_bc = (/nc_ijk(2), nc_ijk(3)/)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_wbc_u, ns_bc, nc_bc, wbc(:,:,1))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_wbc_v, ns_bc, nc_bc, wbc(:,:,2))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_wbc_t, ns_bc, nc_bc, wbc(:,:,3))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_wbc_s, ns_bc, nc_bc, wbc(:,:,4))

  nc_bc = 0
  if(myid_x .eq. npx-1) nc_bc = (/nc_ijk(2), nc_ijk(3)/)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_ebc_u, ns_bc, nc_bc, ebc(:,:,1))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_ebc_v, ns_bc, nc_bc, ebc(:,:,2))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_ebc_t, ns_bc, nc_bc, ebc(:,:,3))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_ebc_s, ns_bc, nc_bc, ebc(:,:,4))

  ns_bc = (/ns_ijk(1), ns_ijk(3)/)
  nc_bc = 0
  if(myid_y .eq. 0) nc_bc = (/nc_ijk(1), nc_ijk(3)/)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_sbc_u, ns_bc, nc_bc, sbc(:,:,1))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_sbc_v, ns_bc, nc_bc, sbc(:,:,2))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_sbc_t, ns_bc, nc_bc, sbc(:,:,3))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_sbc_s, ns_bc, nc_bc, sbc(:,:,4))

  nc_bc = 0
  if(myid_y .eq. npy-1) nc_bc = (/nc_ijk(1), nc_ijk(3)/)
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_nbc_u, ns_bc, nc_bc, nbc(:,:,1))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_nbc_v, ns_bc, nc_bc, nbc(:,:,2))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_nbc_t, ns_bc, nc_bc, nbc(:,:,3))
  ierr = NFMPI_PUT_VARA_DOUBLE_ALL(ncid, id_nbc_s, ns_bc, nc_bc, nbc(:,:,4))

  ierr = NFMPI_CLOSE(ncid)
end subroutine check_lateral_bc

end module tai_timcom_pncio
