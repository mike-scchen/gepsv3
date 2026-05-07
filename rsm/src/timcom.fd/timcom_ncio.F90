module tai_timcom_ncio
  use tai_timcom_const
  use tai_timcom_grid, only: data_panel
  use NETCDF
  implicit none

contains
subroutine get_ncdata_dims(ncdata)
  implicit none

  type(data_panel) :: ncdata
  integer :: i, ierr, dimid, vadmid, ndtmp
  character(len=8) :: dname

  ierr = NF90_OPEN(trim(ncdata%fname), NF90_NOWRITE, ncdata%ncid)

  do i = 1, 3
    ierr = NF90_INQ_DIMID(ncdata%ncid, trim(ncdata%fdims(i)), dimid)
    ierr = NF90_INQUIRE_DIMENSION(ncdata%ncid, dimid, dname, ncdata%nd(i))
  end do

  allocate(ncdata%time(ncdata%nd(1)))
  allocate(ncdata%lon(ncdata%nd(2)))
  allocate(ncdata%lat(ncdata%nd(3)))
  allocate(ncdata%src(ncdata%nd(2),ncdata%nd(3),2))
  allocate(ncdata%srcm(ncdata%nd(2),ncdata%nd(3)))
 
  ierr = NF90_INQ_VARID(ncdata%ncid, trim(ncdata%fvadm(1)), vadmid)
  ierr = NF90_GET_VAR(ncdata%ncid, vadmid, ncdata%time, (/1/), (/ncdata%nd(1)/))
!  ncdata%time = ncdata%time-1.d0  ! should be removed after restart testing

  ierr = NF90_INQ_VARID(ncdata%ncid, trim(ncdata%fvadm(2)), vadmid)
  ierr = NF90_INQUIRE_VARIABLE(ncdata%ncid, vadmid, ndims=ndtmp)
  if(ndtmp .eq. 1) then
    ierr = NF90_GET_VAR(ncdata%ncid, vadmid, ncdata%lon, (/1/), (/ncdata%nd(2)/))
  else
    ierr = NF90_GET_VAR(ncdata%ncid, vadmid, ncdata%lon, (/1,1/), (/ncdata%nd(2),1/))
  end if

  ierr = NF90_INQ_VARID(ncdata%ncid, trim(ncdata%fvadm(3)), vadmid)
  ierr = NF90_INQUIRE_VARIABLE(ncdata%ncid, vadmid, ndims=ndtmp)
  if(ndtmp .eq. 1) then
    ierr = NF90_GET_VAR(ncdata%ncid, vadmid, ncdata%lat, (/1/), (/ncdata%nd(3)/))
  else
    ierr = NF90_GET_VAR(ncdata%ncid, vadmid, ncdata%lat, (/1,1/), (/ncdata%nd(3),1/))
  end if

  ierr = NF90_INQ_VARID(ncdata%ncid, trim(ncdata%fvara), ncdata%varid)
  ierr = NF90_CLOSE(ncdata%ncid)
end subroutine get_ncdata_dims

subroutine get_ncdata_var(ncdata, t_index, src)
  implicit none

  type(data_panel) :: ncdata
  integer, intent(in) :: t_index
  real(r8), dimension(:,:) :: src
  integer :: ncid, nlon, nlat, ierr

  nlon = ncdata%nd(2)
  nlat = ncdata%nd(3)

  ierr = NF90_OPEN(trim(ncdata%fname), NF90_NOWRITE, ncid)
  ierr = NF90_GET_VAR(ncid, ncdata%varid, src, (/1,1,t_index/), (/nlon,nlat,1/))
  ierr = NF90_CLOSE(ncid)

end subroutine get_ncdata_var

subroutine get_lateral_bc_time(ncfile, lateral_bc_ntime, lateral_bc_time)
  implicit none

  character(len=*) :: ncfile
  integer :: lateral_bc_ntime
  real(r8), dimension(:), pointer :: lateral_bc_time
  integer :: ierr, ncid, dimid, ntime, vadmid

  ierr = NF90_OPEN(trim(ncfile), NF90_NOWRITE, ncid)
  ierr = NF90_INQ_DIMID(ncid, 'time', dimid)
  ierr = NF90_INQUIRE_DIMENSION(ncid, dimid, len=ntime)

  allocate(lateral_bc_time(ntime))

  lateral_bc_ntime = ntime

  ierr = NF90_INQ_VARID(ncid, 'time', vadmid)
  ierr = NF90_GET_VAR(ncid, vadmid, lateral_bc_time, (/1/), (/ntime/))

  ierr = NF90_CLOSE(ncid)

end subroutine get_lateral_bc_time
end module tai_timcom_ncio
