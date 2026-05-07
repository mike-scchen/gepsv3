module tai_timcom_dice
  use tai_timcom_const
  use tai_timcom_comm
  use tai_timcom_grid, only: data_panel
  use tai_timcom_general, only: ice_grid_num, timer_panel
  implicit none

  integer, parameter ::& 
    idx_i2x_ifrc =  1, &
    n_i2x_flux   =  1

  character(len=4), dimension(n_i2x_flux) :: &
    iflux_var = (/"ifrc"/)

  type :: ice_panel
    character(len=8) :: tag
    character(len=128) :: nml_file
    integer :: nvar = 0
    type(data_panel) :: flx(n_i2x_flux)
    integer :: nlon, nlat
    real(r8), pointer :: lon(:), lat(:)
  end type

  type(ice_panel), pointer :: ice_grid(:)=>null()
contains

subroutine dice_init(grid_list)
  use tai_timcom_general, only: str_split, member
  implicit none

  character(len=*) :: grid_list

  integer :: i, ierr, n
  logical :: set_ice_grid = .true.
  character(len=256) :: log_info

  if(ice_grid_num .gt. 0) then
    allocate(ice_grid(ice_grid_num))
    do i = 1, ice_grid_num
      call str_split(trim(grid_list), ":", ice_grid(i)%tag, ierr)
      write(ice_grid(i)%nml_file, '(a8,i3.3,a1,a3)') 'namelist',member,'.',trim(ice_grid(i)%tag)
      if(myid_timcom .eq. rootid) then
        write(log_info,'(5a)') &
          '[INFO]: Register ice grid: ', trim(ice_grid(i)%tag), achar(10), & 
          '        Read ice grid namelist: ', trim(ice_grid(i)%nml_file)
        call comm_write_log_info(fid_log, log_info)
      end if

      do n = 1, n_i2x_flux
        call get_dice_var_nml(trim(ice_grid(i)%nml_file), ice_grid(i)%flx(n), n)
        if(set_ice_grid .and. associated(ice_grid(i)%flx(n)%src)) then
          ice_grid(i)%nlon = ice_grid(i)%flx(n)%nd(2)
          ice_grid(i)%nlat = ice_grid(i)%flx(n)%nd(3)
          ice_grid(i)%lon => ice_grid(i)%flx(n)%lon
          ice_grid(i)%lat => ice_grid(i)%flx(n)%lat
          set_ice_grid = .false.
        end if
      end do
    end do
  else
    if(myid_timcom .eq. rootid) then
      write(log_info,'(a,a)') '[INFO]: There is no assigned ice grid. ', trim(grid_list)
      call comm_write_log_info(fid_log, log_info)
    end if
  end if

  !if(myid_timcom .eq. rootid) then
  !  write(*,'(a,2i6)') "[INFO]: ice grid dimension: ", ice_grid(1)%nlon, ice_grid(1)%nlat
  !  write(*,'(a)') "[INFO]: ice grid lat: "
  !  write(*,'(10f12.6)') ice_grid(1)%lon
  !  write(*,'(a)') "[INFO]: ice grid lon: "
  !  write(*,'(10f12.6)') ice_grid(1)%lat
  !  write(*,'(a)') "[INFO]: ice grid time: "
  !  write(*,'(10f12.6)') ice_grid(1)%flx(1)%time
  ! end if
end subroutine dice_init

subroutine dice_exec(timer)
  use tai_timcom_general, only: integral_var2d
  implicit none

  type(timer_panel) :: timer
  integer :: i, ierr, n

  do i = 1, ice_grid_num
    do n = 1, n_i2x_flux
      if(associated(ice_grid(i)%flx(n)%src)) then
        call integral_var2d(ice_grid(i)%flx(n), timer%mjd1, timer%mjd2,2)
      end if
    end do
  end do
end subroutine dice_exec

subroutine get_dice_var_nml(nml_file, dice_data, nid)
  use tai_timcom_ncio, only: get_ncdata_dims
  implicit none

  character(len=*), intent(in) :: nml_file
  type(data_panel) :: dice_data
  integer, intent(in) :: nid

  character(len=128) :: fpath
  character(len=32)  :: fname
  character(len=8)   :: gridx, gridy, time, vara
  integer :: fid = 12, ierr = -1
  character(len=256) :: log_info

  namelist /dice_ifrac/ fpath, fname, gridx, gridy, time, vara

  open(unit=fid, file=trim(nml_file), status='OLD')
  select case(nid)
  case(idx_i2x_ifrc)
    read(unit=fid, nml=dice_ifrac, iostat=ierr)
  case default
    continue
  end select

  close(fid)

  if(ierr .eq. 0) then
    dice_data%fname = trim(fpath)//trim(fname)
    dice_data%fvara = trim(vara)
    dice_data%fdims(1) = trim(time)
    dice_data%fdims(2) = trim(gridx)
    dice_data%fdims(3) = trim(gridy)
    dice_data%fvadm(1) = trim(time)
    dice_data%fvadm(2) = trim(gridx)
    dice_data%fvadm(3) = trim(gridy)

    if(myid_timcom .eq. rootid ) then
      write(log_info,'(10a)') &
        'find ', trim(dice_data%fvara), &
        '(',  trim(dice_data%fdims(2)), &
        ', ', trim(dice_data%fdims(3)), &
        ', ', trim(dice_data%fdims(1)), &
        ') in the file ', trim(dice_data%fname)
      call comm_write_log_info(fid_log, log_info)
    end if
    call get_ncdata_dims(dice_data)
    dice_data%t(:) = -1
  end if
end subroutine get_dice_var_nml

subroutine dice_to_timcom
  use tai_timcom_grid, only: ocn_grid
  implicit none

  integer, parameter :: i = 1
  integer :: x1, x2, y1, y2

  x1 = ocn_grid(i)%mpicom%ns_ijk(1)
  x2 = x1+ocn_grid(i)%nx - 1

  y1 = ocn_grid(i)%mpicom%ns_ijk(2)
  y2 = y1+ocn_grid(i)%ny - 1

  ocn_grid(i)%ifrc = ice_grid(i)%flx(idx_i2x_ifrc)%srcm(x1:x2,y1:y2)

end subroutine dice_to_timcom

end module tai_timcom_dice

