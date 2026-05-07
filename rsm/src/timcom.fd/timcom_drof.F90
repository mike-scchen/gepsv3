module tai_timcom_drof
  use tai_timcom_const
  use tai_timcom_comm
  use tai_timcom_grid, only: data_panel
  use tai_timcom_general, only: rof_grid_num, timer_panel
  implicit none

  integer, parameter ::& 
    idx_r2x_roff =  1, &
    n_r2x_flux   =  1

  character(len=4), dimension(n_r2x_flux) :: &
    rflux_var = (/"roff"/)

  type :: rof_panel
    character(len=8) :: tag
    character(len=128) :: nml_file
    integer :: nvar = 0
    type(data_panel) :: flx(n_r2x_flux)
    integer :: nlon, nlat
    real(r8), pointer :: lon(:), lat(:)
  end type

  type(rof_panel), pointer :: rof_grid(:)=>null()
contains

subroutine drof_init(grid_list)
  use tai_timcom_general, only: str_split, member
  implicit none

  character(len=*) :: grid_list

  integer :: i, ierr, n
  logical :: set_rof_grid = .true.
  character(len=256) :: log_info

  if(rof_grid_num .gt. 0) then
    allocate(rof_grid(rof_grid_num))
    do i = 1, rof_grid_num
      call str_split(trim(grid_list), ":", rof_grid(i)%tag, ierr)
      write(rof_grid(i)%nml_file, '(a8,i3.3,a1,a3)') 'namelist',member,'.',trim(rof_grid(i)%tag)
      if(myid_timcom .eq. rootid) then
        write(log_info,'(5a)') &
          '[INFO]: Register rof grid: ', trim(rof_grid(i)%tag), achar(10), &
          '        Read rof grid namelist: ', trim(rof_grid(i)%nml_file)
        call comm_write_log_info(fid_log, log_info)
      end if
      do n = 1, n_r2x_flux
        call get_drof_var_nml(trim(rof_grid(i)%nml_file), rof_grid(i)%flx(n), n)
        if(set_rof_grid .and. associated(rof_grid(i)%flx(n)%src)) then
          rof_grid(i)%nlon = rof_grid(i)%flx(n)%nd(2)
          rof_grid(i)%nlat = rof_grid(i)%flx(n)%nd(3)
          rof_grid(i)%lon => rof_grid(i)%flx(n)%lon
          rof_grid(i)%lat => rof_grid(i)%flx(n)%lat
          set_rof_grid = .false.
        end if
      end do
    end do
  else
    if(myid_timcom .eq. rootid) then
      write(log_info,'(a,a)') '[INFO]: There is no assigned rof grid. ', trim(grid_list)
      call comm_write_log_info(fid_log, log_info)
    end if
  end if

  !if(myid_timcom .eq. rootid) then
  !  write(*,'(a,2i6)') "[INFO]: rof grid dimension: ", rof_grid(1)%nlon, rof_grid(1)%nlat
  !  write(*,'(a)') "[INFO]: rof grid lat: "
  !  write(*,'(10f12.6)') rof_grid(1)%lon
  !  write(*,'(a)') "[INFO]: rof grid lon: "
  !  write(*,'(10f12.6)') rof_grid(1)%lat
  !  write(*,'(a)') "[INFO]: rof grid time: "
  !  write(*,'(10f12.6)') rof_grid(1)%flx(1)%time
  ! end if
end subroutine drof_init

subroutine drof_exec(timer)
  use tai_timcom_general, only: integral_var2d
  implicit none

  type(timer_panel) :: timer
  integer :: i, ierr, n

  do i = 1, rof_grid_num
    do n = 1, n_r2x_flux
      if(associated(rof_grid(i)%flx(n)%src)) then
        call integral_var2d(rof_grid(i)%flx(n), timer%mjd1, timer%mjd2, 4)
      end if
    end do
  end do
end subroutine drof_exec

subroutine get_drof_var_nml(nml_file, drof_data, nid)
  use tai_timcom_ncio, only: get_ncdata_dims
  implicit none

  character(len=*), intent(in) :: nml_file
  type(data_panel) :: drof_data
  integer, intent(in) :: nid

  character(len=128) :: fpath
  character(len=32)  :: fname
  character(len=8)   :: gridx, gridy, time, vara
  integer :: fid = 13, ierr = -1
  character(len=256) :: log_info

  namelist /drof_roff/ fpath, fname, gridx, gridy, time, vara

  open(unit=fid, file=trim(nml_file), status='OLD')
  select case(nid)
  case(idx_r2x_roff)
    read(unit=fid, nml=drof_roff, iostat=ierr)
  case default
    continue
  end select

  close(fid)

  if(ierr .eq. 0) then
    drof_data%fname = trim(fpath)//trim(fname)
    drof_data%fvara = trim(vara)
    drof_data%fdims(1) = trim(time)
    drof_data%fdims(2) = trim(gridx)
    drof_data%fdims(3) = trim(gridy)
    drof_data%fvadm(1) = trim(time)
    drof_data%fvadm(2) = trim(gridx)
    drof_data%fvadm(3) = trim(gridy)

    if(myid_timcom .eq. rootid ) then
      write(log_info,'(10a)') &
        'find ', trim(drof_data%fvara), &
        '(',  trim(drof_data%fdims(2)), &
        ', ', trim(drof_data%fdims(3)), &
        ', ', trim(drof_data%fdims(1)), &
        ') in the file ', trim(drof_data%fname)
      call comm_write_log_info(fid_log, log_info)
    end if
    call get_ncdata_dims(drof_data)
    drof_data%t(:) = -1
  end if
end subroutine get_drof_var_nml

subroutine drof_to_timcom
  use tai_timcom_grid, only: ocn_grid
  implicit none

  integer, parameter :: i = 1
  integer :: x1, x2, y1, y2

  x1 = ocn_grid(i)%mpicom%ns_ijk(1)
  x2 = x1+ocn_grid(i)%nx - 1

  y1 = ocn_grid(i)%mpicom%ns_ijk(2)
  y2 = y1+ocn_grid(i)%ny - 1

  ocn_grid(i)%roff = 0.1d0*rof_grid(i)%flx(idx_r2x_roff)%srcm(x1:x2,y1:y2)

end subroutine drof_to_timcom

end module tai_timcom_drof

