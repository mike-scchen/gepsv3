module timcom_datm
  use timcom_const
  use timcom_comm
  use timcom_grid, only: data_panel
  use timcom_general, only: atm_grid_num,  timer_panel
  implicit none

  integer, parameter ::& 
    idx_a2x_uwnd =  1, &
    idx_a2x_vwnd =  2, &
    idx_a2x_temp =  3, &
    idx_a2x_qhum =  4, &
    idx_a2x_pslv =  5, &
    idx_a2x_swup =  6, &
    idx_a2x_swdn =  7, &
    idx_a2x_lwup =  8, &
    idx_a2x_lwdn =  9, &
    idx_a2x_lath = 10, &
    idx_a2x_senh = 11, &
    idx_a2x_rain = 12, &
    idx_a2x_snow = 13, &
    n_a2x_flux   = 13

  character(len=4), dimension(n_a2x_flux) :: &
    aflux_var = (/"uwnd", "vwnd", "temp", "qhum", "pslv", &
                  "swup", "swdn", "lwup", "lwdn",         &
                  "lath", "senh", "rain", "snow" /)

  type :: atm_panel
    character(len=8) :: tag
    character(len=128) :: nml_file
    integer :: nvar = 0
    type(data_panel) :: flx(n_a2x_flux)
    integer :: nlon, nlat
    real(r8), pointer :: lon(:), lat(:)
  end type

  type(atm_panel), pointer :: atm_grid(:)=>null()
contains

subroutine datm_init(grid_list)
  use timcom_general, only: str_split, member, atm_grid_num
  implicit none

  character(len=*) :: grid_list

  integer :: i, ierr, n
  logical :: set_atm_grid = .true.
  character(len=256) :: log_info

  if(atm_grid_num .gt. 0) then
    allocate(atm_grid(atm_grid_num))
    do i = 1, atm_grid_num
      call str_split(trim(grid_list), ":", atm_grid(i)%tag, ierr)
      write(atm_grid(i)%nml_file, '(a8,i3.3,a1,a)') 'namelist',member,'.',trim(atm_grid(i)%tag)
      if(myid_timcom .eq. rootid) then
        write(log_info,'(5a)') &
          '[INFO]: Register atm grid: ', trim(atm_grid(i)%tag), achar(10), &
          '        Read atm grid namelist: ', trim(atm_grid(i)%nml_file)
        call comm_write_log_info(fid_log, log_info)
      end if
      do n = 1, n_a2x_flux
        call get_datm_var_nml(trim(atm_grid(i)%nml_file), atm_grid(i)%flx(n), n)
        if(set_atm_grid .and. associated(atm_grid(i)%flx(n)%src)) then
          atm_grid(i)%nlon = atm_grid(i)%flx(n)%nd(2)
          atm_grid(i)%nlat = atm_grid(i)%flx(n)%nd(3)
          atm_grid(i)%lon => atm_grid(i)%flx(n)%lon
          atm_grid(i)%lat => atm_grid(i)%flx(n)%lat
          set_atm_grid = .false.
        end if
      end do
    end do
  else
    if(myid_timcom .eq. rootid) then
      write(log_info,'(a,a)') '[INFO]: There is no assigned atm grid. ', trim(grid_list)
      call comm_write_log_info(fid_log, log_info)
    end if
  end if

  !if(myid_timcom .eq. rootid) then
  !  write(*,'(a,2i6)') "[INFO]: atm grid dimension: ", atm_grid(1)%nlon, atm_grid(1)%nlat
  !  write(*,'(a)') "[INFO]: atm grid lon: "
  !  write(*,'(10f12.6)') atm_grid(1)%lon
  !  write(*,'(a)') "[INFO]: atm grid lat: "
  !  write(*,'(10f12.6)') atm_grid(1)%lat
  !  write(*,'(a)') "[INFO]: atm grid time: "
  !  write(*,'(10f12.6)') atm_grid(1)%flx(1)%time
  !end if
end subroutine datm_init

subroutine datm_exec(timer)
  use timcom_general, only: integral_var2d
  implicit none

  type(timer_panel) :: timer
  integer :: i, ierr, n

  do i = 1, atm_grid_num
    do n = 1, n_a2x_flux
      if(associated(atm_grid(i)%flx(n)%src)) then
#ifdef cpl
 #ifndef clm_r
  call integral_var2d(atm_grid(i)%flx(n), timer%mjd1, timer%mjd2,1)
 #else
  call integral_var2d(atm_grid(i)%flx(n), timer%mjd1, timer%mjd2, timer%shift_mjd, 1)
!  call integral_var2d(atm_grid(i)%flx(n), timer%syng_mjd, timer%syng_mjd+timer%dt_mjd, timer%shift_mjd, 1)
 #endif
#else
  call integral_var2d(atm_grid(i)%flx(n), timer%mjd1, timer%mjd2, timer%shift_mjd, 1)
#endif

      end if
    end do
  end do
end subroutine datm_exec

subroutine get_datm_var_nml(nml_file, datm_data, nid)
  use timcom_ncio, only: get_ncdata_dims
  implicit none

  character(len=*), intent(in) :: nml_file
  type(data_panel) :: datm_data
  integer, intent(in) :: nid

  character(len=128) :: fpath
  character(len=128) :: fname
  character(len=8)   :: gridx, gridy, time, vara
  integer :: fid = 11, ierr = -1
  character(len=256) :: log_info

  namelist /datm_uwnd/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_vwnd/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_temp/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_qhum/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_pslv/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_qhum/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_swup/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_swdn/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_lwup/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_lwdn/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_lath/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_senh/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_rain/ fpath, fname, gridx, gridy, time, vara
  namelist /datm_snow/ fpath, fname, gridx, gridy, time, vara

  open(unit=fid, file=trim(nml_file), status='OLD')
  select case(nid)
  case(idx_a2x_uwnd)
    read(unit=fid, nml=datm_uwnd, iostat=ierr)
  case(idx_a2x_vwnd)
    read(unit=fid, nml=datm_vwnd, iostat=ierr)
  case(idx_a2x_temp)
    read(unit=fid, nml=datm_temp, iostat=ierr)
  case(idx_a2x_qhum)
    read(unit=fid, nml=datm_qhum, iostat=ierr)
  case(idx_a2x_pslv)
    read(unit=fid, nml=datm_pslv, iostat=ierr)
  case(idx_a2x_swup)
    read(unit=fid, nml=datm_swup, iostat=ierr)
  case(idx_a2x_swdn)
    read(unit=fid, nml=datm_swdn, iostat=ierr)
  case(idx_a2x_lwup)
    read(unit=fid, nml=datm_lwup, iostat=ierr)
  case(idx_a2x_lwdn)
    read(unit=fid, nml=datm_lwdn, iostat=ierr)
  case(idx_a2x_lath)
    read(unit=fid, nml=datm_lath, iostat=ierr)
  case(idx_a2x_senh)
    read(unit=fid, nml=datm_senh, iostat=ierr)
  case(idx_a2x_rain)
    read(unit=fid, nml=datm_rain, iostat=ierr)
  case(idx_a2x_snow)
    read(unit=fid, nml=datm_snow, iostat=ierr)
  case default
    continue
  end select

  close(fid)

  if(ierr .eq. 0) then
    datm_data%fname = trim(fpath)//trim(fname)
    datm_data%fvara = trim(vara)
    datm_data%fdims(1) = trim(time)
    datm_data%fdims(2) = trim(gridx)
    datm_data%fdims(3) = trim(gridy)
    datm_data%fvadm(1) = trim(time)
    datm_data%fvadm(2) = trim(gridx)
    datm_data%fvadm(3) = trim(gridy)

    if(myid_timcom .eq. rootid ) then
      write(log_info,'(10a)') &
        'find ', trim(datm_data%fvara), &
        '(',  trim(datm_data%fdims(2)), &
        ', ', trim(datm_data%fdims(3)), &
        ', ', trim(datm_data%fdims(1)), &
        ') in the file ', trim(datm_data%fname)
      call comm_write_log_info(fid_log, log_info)
    end if
    call get_ncdata_dims(datm_data)
    datm_data%t(:) = -1
  end if
end subroutine get_datm_var_nml

subroutine datm_to_timcom
  use timcom_grid, only: ocn_grid
  implicit none

  integer :: i, j, n
  integer :: x1, x2, y1, y2
  real(r8) :: avg_alb

  x1 = ocn_grid(i)%mpicom%ns_ijk(1)
  x2 = x1+ocn_grid(i)%nx - 1

  y1 = ocn_grid(i)%mpicom%ns_ijk(2)
  y2 = y1+ocn_grid(i)%ny - 1

  i = 1
#ifdef cpl
  ocn_grid(i)%u_10 = atm_grid(i)%flx(idx_a2x_uwnd)%srcm(x1:x2,y1:y2)
  ocn_grid(i)%v_10 = atm_grid(i)%flx(idx_a2x_vwnd)%srcm(x1:x2,y1:y2)
  ocn_grid(i)%t_10 = atm_grid(i)%flx(idx_a2x_temp)%srcm(x1:x2,y1:y2)
  ocn_grid(i)%q_10 = atm_grid(i)%flx(idx_a2x_qhum)%srcm(x1:x2,y1:y2)
  ocn_grid(i)%pslv = atm_grid(i)%flx(idx_a2x_pslv)%srcm(x1:x2,y1:y2)
  ocn_grid(i)%swup = -atm_grid(i)%flx(idx_a2x_swup)%srcm(x1:x2,y1:y2)
  ocn_grid(i)%swdn = atm_grid(i)%flx(idx_a2x_swdn)%srcm(x1:x2,y1:y2)
  ocn_grid(i)%lwdn = atm_grid(i)%flx(idx_a2x_lwdn)%srcm(x1:x2,y1:y2)
  ocn_grid(i)%rain = atm_grid(i)%flx(idx_a2x_rain)%srcm(x1:x2,y1:y2)
  ocn_grid(i)%snow = 0.d0 !atm_grid(i)%flx(idx_a2x_snow)%srcm(x1:x2,y1:y2)
#else
  ocn_grid(i)%u_10 = dble(atm_grid(i)%flx(idx_a2x_uwnd)%srcm(x1:x2,y1:y2))
  ocn_grid(i)%v_10 = dble(atm_grid(i)%flx(idx_a2x_vwnd)%srcm(x1:x2,y1:y2))
  ocn_grid(i)%t_10 = dble(atm_grid(i)%flx(idx_a2x_temp)%srcm(x1:x2,y1:y2))
  ocn_grid(i)%q_10 = dble(atm_grid(i)%flx(idx_a2x_qhum)%srcm(x1:x2,y1:y2))
  ocn_grid(i)%pslv = dble(atm_grid(i)%flx(idx_a2x_pslv)%srcm(x1:x2,y1:y2))
!  ocn_grid(i)%swup = dble(-atm_grid(i)%flx(idx_a2x_swup)%srcm(x1:x2,y1:y2))
  ocn_grid(i)%swdn = dble(atm_grid(i)%flx(idx_a2x_swdn)%srcm(x1:x2,y1:y2))
  ocn_grid(i)%lwdn = dble(atm_grid(i)%flx(idx_a2x_lwdn)%srcm(x1:x2,y1:y2))
  ocn_grid(i)%rain = dble(atm_grid(i)%flx(idx_a2x_rain)%srcm(x1:x2,y1:y2))
  ocn_grid(i)%snow = 0.d0 !atm_grid(i)%flx(idx_a2x_snow)%srcm(x1:x2,y1:y2)

  n = 1
  do j = 1, ocn_grid(n)%ny
      avg_alb = 0.069 - 0.11*dcos(2.d0*ocn_grid(n)%y_grid(j)*d2r)
      do i = 1, ocn_grid(n)%nx
        ocn_grid(n)%swup(i,j) = -ocn_grid(n)%swdn(i,j)*avg_alb
      end do
  end do

  do j = 1, ocn_grid(n)%ny
      do i = 1, ocn_grid(n)%nx
        if(ocn_grid(n)%t_10(i,j) .lt. 273.15d0) then
           ocn_grid(n)%snow(i,j) = ocn_grid(n)%rain(i,j)
           ocn_grid(n)%rain(i,j) = 0.d0
        end if
      end do
  end do
#endif
 
end subroutine datm_to_timcom

end module timcom_datm

