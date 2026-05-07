module tai_timcom_general
  use tai_timcom_const
  use tai_timcom_comm
  use tai_timcom_grid, only: data_panel
  implicit none

  character(len=256) :: path_output, path_restart
  logical :: restart
  integer :: member = 0

  integer ::          & 
    ocn_grid_num = 0, &
    atm_grid_num = 0, &
    rof_grid_num = 0, &
    ice_grid_num = 0, &
    wav_grid_num = 0

  type :: timer_panel
    integer :: daodt, it0, mxit
    integer :: greg1(4), greg2(4), syng_mon
    real(r8) :: mjd1, mjd2, dt_mjd, syng_mjd
    real(r8) :: base_mjd, init_mjd, dest_mjd
    logical :: stopsig = .false., &
               end_day = .false., &
               end_mn  = .false., &
               end_yr  = .false., &
               chk_out = .false.
    character(len=32) :: mjd_unit
  end type timer_panel

  type(timer_panel) :: timer_drv
contains

subroutine general_setup(grid_list)
  implicit none

  character(len=*)  :: grid_list(4)
  integer :: ierr
!  write(*,*) 'pre-get_member_from_args' 
  call get_member_from_args
!  write(*,*) 'pre-get_namelist_run', grid_list
  call get_namelist_run(grid_list)

end subroutine general_setup

subroutine get_member_from_args
  implicit none

  character(len=8) :: arg
  character(len=256) :: log_info
  integer :: istat, iargc
  
!  write(*,*) 'pre-getarg', 'arg = ', arg
  call getarg(1, arg)

!  write(*,*) 'iargc', iargc(), 'myid_timcom', myid_timcom, 'rootid', rootid

  if(iargc() .eq. 0) then
    if(myid_timcom .eq. rootid) then
      write(log_info,'(a)') "[ERROR]: There is no member number in the command line argument."
      call comm_write_log_info(fid_log, log_info)
    end if
!    write(*,*) 'pre-comm_finalize1'
    call comm_finalize
  else
    read(arg, *, iostat=istat) member

 !  write(*,*) 'istat', istat, 'myid_timcom', myid_timcom, 'rootid', rootid, 'member', member

    if(istat /= 0 .or. member > 999) then
      if(myid_timcom .eq. rootid) then
        write(log_info,'(a)') "[ERROR]: the datatype of member should be a integer from 0 to 999"
        call comm_write_log_info(fid_log, log_info)
      end if
!      write(*,*) 'pre-comm_finalize2'
      call comm_finalize
    else
      if(myid_timcom .eq. rootid) then
        write(log_info,'(a,i4.3)') "[INFO]: TIMCOM MEMBER", member
        call comm_write_log_info(fid_log, log_info)
      end if
    end if
  end if
end subroutine get_member_from_args

subroutine get_namelist_run(grid_list)
  use tai_calendar
  implicit none

  character(len=*) :: grid_list(4)
  character(len=32) :: filename
  character(len=256) :: &
    ocn_grid_list, atm_grid_list, &
    rof_grid_list, ice_grid_list, &
    log_info
  integer :: daodt_drv,    duration
  character(len=8) :: duration_unit
  integer :: fid = 10, ierr
  integer :: dtg_init

  namelist /general_info/ &
    path_output, &
    ocn_grid_num, ocn_grid_list, &
    atm_grid_num, atm_grid_list, &
    rof_grid_num, rof_grid_list, &
    ice_grid_num, ice_grid_list, &
    daodt_drv,    dtg_init,     &
    duration,     duration_unit

  namelist /restart_info/ &
    restart, path_restart

  write(filename, '(A8,i3.3,A4)') 'namelist',member,'.run'
  open(unit=fid, file=trim(filename), status='OLD', iostat=ierr)
  if(ierr /= 0) then
    if(myid_timcom .eq. rootid) then
       write(log_info,'(a,a)') '[ERROR]: There is no ', trim(filename)
       call comm_write_log_info(fid_log, log_info)
    end if
    close(fid)
    call comm_finalize
  else
    read(unit=fid, nml=general_info)
    read(unit=fid, nml=restart_info)
    close(fid)
    if(myid_timcom .eq. rootid) then
      write(log_info,'(a,a)') '[INFO]: Reading ', trim(filename)
      call comm_write_log_info(fid_log, log_info)
    end if
  end if

  grid_list(1) = trim(ocn_grid_list)
  grid_list(2) = trim(atm_grid_list)
  grid_list(3) = trim(ice_grid_list)
  grid_list(4) = trim(rof_grid_list)

  call dtg2greg(dtg_init, timer_drv%greg1)
  write(timer_drv%mjd_unit, &
    '("days since ",i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)', iostat=ierr) &
     timer_drv%greg1, 0, 0

  timer_drv%greg2(:) = timer_drv%greg1(:)
  select case(trim(duration_unit))
  case('days')
    timer_drv%greg2(3) = timer_drv%greg1(3) + duration
  case('hour')
    timer_drv%greg2(4) = timer_drv%greg1(4) + duration
  case default
    if(myid_timcom .eq. rootid) then
      write(log_info,'(a)') '[ERROR]: duration unit: hour or days'
      call comm_write_log_info(fid_log, log_info)
    end if
    call comm_finalize
  end select

!0719
!  timer_drv%base_mjd = mjd2(timer_drv%greg1(1), timer_drv%greg1(2), timer_drv%greg1(3), timer_drv%greg1(4), 0)
  timer_drv%base_mjd = mjd2(timer_drv%greg1(1), 1, 1, 0, 0)
  timer_drv%init_mjd = mjd2(timer_drv%greg1(1), timer_drv%greg1(2), timer_drv%greg1(3), timer_drv%greg1(4), 0) - timer_drv%base_mjd
  timer_drv%dest_mjd = mjd2(timer_drv%greg2(1), timer_drv%greg2(2), timer_drv%greg2(3), timer_drv%greg2(4), 0) - timer_drv%base_mjd
  timer_drv%mjd1 = timer_drv%init_mjd
  timer_drv%mjd2 = timer_drv%init_mjd

  timer_drv%daodt = daodt_drv
  timer_drv%dt_mjd = 1.d0/dble(daodt_drv)
  timer_drv%it0 = 0
  timer_drv%mxit = idint(timer_drv%dest_mjd-timer_drv%init_mjd)*timer_drv%daodt
 
  if(myid_timcom .eq. rootid) then
    write(log_info,'(a,i4)') 'daodt_drv: ', timer_drv%daodt
    call comm_write_log_info(fid_log, log_info)
    write(log_info,'(a,f12.4,2a)') 'init mjd: ' , timer_drv%init_mjd, ' ', trim(timer_drv%mjd_unit)
    call comm_write_log_info(fid_log, log_info)
    write(log_info,'(a,f12.4,2a)') 'dest mjd: ' , timer_drv%dest_mjd, ' ', trim(timer_drv%mjd_unit)
    call comm_write_log_info(fid_log, log_info)
    write(log_info,'(a,i6, a)') 'duration: ' , duration, ' days'
    call comm_write_log_info(fid_log, log_info)
    write(log_info,'(a,4i4)') 'init greg: ', timer_drv%greg1
    call comm_write_log_info(fid_log, log_info)
    write(log_info,'(a,4i4)') 'dest greg: ', timer_drv%greg2
    call comm_write_log_info(fid_log, log_info)
    write(log_info,'(a,i8)')  'init itf: ',  timer_drv%it0
    call comm_write_log_info(fid_log, log_info)
    write(log_info,'(a,i8)')  'dest itf: ',  timer_drv%mxit
    call comm_write_log_info(fid_log, log_info)
  end if

end subroutine get_namelist_run

subroutine update_timer(itf, timer)
  implicit none

  integer, intent(in) :: itf
  type(timer_panel) :: timer
  character(len=256) :: log_info

!0719
!  timer%mjd1 = timer%dt_mjd*dble(itf-1)
!  timer%mjd2 = timer%dt_mjd*dble(itf)
  timer%mjd1 = timer%init_mjd + timer%dt_mjd*dble(itf-1)
  timer%mjd2 = timer%init_mjd + timer%dt_mjd*dble(itf)
  timer%syng_mjd = dmod(timer%mjd2,365.d0)
  timer%syng_mon = idint(12.d0*timer%syng_mjd/365.d0) + 1

  if(itf .gt. timer%mxit) then
    if(myid_timcom .eq. rootid) then
      write(log_info,'(a,2f12.3)') "[ERROR]: out of duration: ", timer%mjd2+timer%init_mjd, timer%dest_mjd
      call comm_write_log_info(fid_log, log_info)
    end if
  else
    if(myid_timcom .eq. rootid) then
      write(log_info,'(a,f12.3,a3,f12.3,a,i8)') "[INFO]: forecast days:", timer%mjd1, " to", timer%mjd2, ", steps:", itf
      call comm_write_log_info(fid_log, log_info)
    end if
  endif
end subroutine update_timer

subroutine str_split(majstr, delimiter, substr, ierr)
  implicit none

  character(len=256) :: majstr
  character(len=1)   :: delimiter
  character(len=8)   :: substr
  integer, intent(out) :: ierr

  integer :: loc

  ierr = 0
  loc = index(majstr, delimiter)
  if(loc .eq. 0) then
    ierr = 1
  else
    substr = majstr(1:loc-1)
    majstr = majstr(loc+1:256)
  end if
end subroutine str_split

subroutine dtg2greg(dtg, greg)
  implicit none

  integer, intent(in)  :: dtg
  integer, intent(out) :: greg(4)
  integer :: dtg_tmp

  dtg_tmp = dtg

  greg(1) = dtg_tmp/1000000
  dtg_tmp = dtg_tmp - greg(1)*1000000
  greg(2) = dtg_tmp/10000
  dtg_tmp = dtg_tmp - greg(2)*10000
  greg(3) = dtg_tmp/100
  dtg_tmp = dtg_tmp - greg(3)*100
  greg(4) = dtg_tmp

  if(greg(1) .lt. 100) then
    if(greg(1) .ge. 48) then
      greg(1) = greg(1) + 1900
    else
      greg(1) = greg(1) + 2000
    end if
  end if
end subroutine dtg2greg

subroutine integral_var2d(var2d, mjd1_drv, mjd2_drv, opt)
  use tai_timcom_ncio, only: get_ncdata_var
  implicit none

  type(data_panel) :: var2d
  real(r8) :: mjd1_drv, mjd2_drv, mjd1, mjd2
  integer, optional :: opt
  real(r8) :: dt_var2d, tmp1, tmp2
  integer :: t, mjdt1(2), mjdt2(2)
  integer :: len_data, ierr
  character(len=256) :: log_info

  mjd1 = dmod(mjd1_drv, 365.d0)
  mjd2 = mjd1 + (mjd2_drv - mjd1_drv)

  mjdt1 = -1
  mjdt2 = -1
  dt_var2d = 0
  var2d%srcm = 0.d0

! if(myid_timcom .eq. rootid) then
  do t = 1, var2d%nd(1)-1
    if(mjd1 .ge. var2d%time(t) .and. &
       mjd1 .lt. var2d%time(t+1)) then
       mjdt1(1) = t
       mjdt1(2) = t + 1
       exit
    end if
  end do
  do t = 1, var2d%nd(1)-1
    if(mjd2 .gt. var2d%time(t) .and. &
       mjd2 .le. var2d%time(t+1)) then
       mjdt2(1) = t
       mjdt2(2) = t + 1
       exit
    end if
  end do

  if(mjd1 .lt. var2d%time(1)) then
    mjdt1(1:2) = 1
  end if
  if(mjd1 .ge. var2d%time(var2d%nd(1))) then
    mjdt1(1:2) = var2d%nd(1)
  end if
  if(mjd2 .le. var2d%time(1)) then
    mjdt2(1:2) = 1
  end if
  if(mjd2 .gt. var2d%time(var2d%nd(1))) then
    mjdt2(1:2) = var2d%nd(1)
  end if

  if(mjdt1(1) .eq. -1 .or. mjdt2(1) .eq. -1) then
    if(myid_timcom .eq. rootid) then
      write(log_info,*) "[ERROR]: data ", trim(var2d%fvara), &
              " out of time!!", mjd1, mjd2, var2d%time(1)
      call comm_write_log_info(fid_log, log_info)
    end if
    call comm_finalize
  end if
  
  if(mjdt1(1) .ne. var2d%t(1)) then
    if(mjdt1(1) .eq. var2d%t(2)) then
      var2d%t(1) = var2d%t(2)
      var2d%src(:,:,1) = var2d%src(:,:,2)
    else
      var2d%t(1) = mjdt1(1)
      call get_ncdata_var(var2d, mjdt1(1), var2d%src(:,:,1))
    end if
    var2d%t(2) = mjdt1(2)
    call get_ncdata_var(var2d, mjdt1(2), var2d%src(:,:,2))
  end if

  tmp1 = (var2d%time(mjdt1(2)) - mjd1)
  if(mjdt1(1) .lt. mjdt1(2)) then
    dt_var2d = var2d%time(mjdt1(2)) - var2d%time(mjdt1(1))
    tmp1 = tmp1/dt_var2d
    tmp2 = 2.d0 - tmp1
    var2d%srcm = (tmp1*var2d%src(:,:,1)+tmp2*var2d%src(:,:,2))*tmp1*0.5d0*dt_var2d
  else
    var2d%srcm = var2d%src(:,:,1)*tmp1
  end if
  
  if(mjdt2(1) .eq. mjdt2(2) .and. mjdt2(2) .eq. var2d%nd(1)) then
    tmp1 = (mjd2 - var2d%time(mjdt2(2)))
    var2d%srcm = var2d%srcm + var2d%src(:,:,2)*tmp1
    var2d%srcm = var2d%srcm/(mjd2-mjd1)
  else
    if(mjdt1(2) .le. mjdt2(1) .and. mjdt2(1) .ne. mjdt2(2)) then
      do t = mjdt1(2), mjdt2(1)
        var2d%t(1) = var2d%t(2)
        var2d%src(:,:,1) = var2d%src(:,:,2)
        var2d%t(2) = t + 1
        call get_ncdata_var(var2d, t+1, var2d%src(:,:,2))

        dt_var2d = var2d%time(t+1) - var2d%time(t)
        var2d%srcm = var2d%srcm + sum(var2d%src,3)*dt_var2d*0.5d0
      end do
    end if
  
    tmp1 = (var2d%time(mjdt2(2)) - mjd2)
    if(mjdt2(1) .lt. mjdt2(2)) then
      dt_var2d = var2d%time(mjdt2(2)) - var2d%time(mjdt2(1))
      tmp1 = tmp1/dt_var2d
      tmp2 = 2.d0 - tmp1
      var2d%srcm = var2d%srcm - &
                     (tmp1*var2d%src(:,:,1)+tmp2*var2d%src(:,:,2))*tmp1*0.5d0*dt_var2d
      var2d%srcm = var2d%srcm/(mjd2-mjd1)
    else
      var2d%srcm = var2d%srcm - var2d%src(:,:,2)*tmp1
      var2d%srcm = var2d%srcm/(mjd2-mjd1)
    end if
  end if
  
!  end if
!  len_data = var2d%nd(2)*var2d%nd(3)
!  call MPI_BCAST(var2d%srcm,  len_data, MPI_REAL8, 0, m_comm_timcom, ierr)
!  call MPI_BCAST(var2d%src, 2*len_data, MPI_REAL8, 0, m_comm_timcom, ierr)
!  call MPI_BCAST(mjdt1, 2, MPI_INT, 0, m_comm_timcom, ierr)
!  call MPI_BCAST(mjdt2, 2, MPI_INT, 0, m_comm_timcom, ierr)

  if(myid_timcom .eq. rootid .and. present(opt)) then
     write(log_info,'(a,a,a,f12.4,a3,f12.4,a4)') &
       "[INFO]: integral ", trim(var2d%fvara), " from:", mjd1, " to", mjd2, " days"
     call comm_write_log_info(fid_log, log_info)
  end if
end subroutine integral_var2d

end module tai_timcom_general
