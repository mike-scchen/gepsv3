subroutine tai_timcom_main(compid, mpi_comm_mct)
  implicit none

  integer, intent(in) :: compid, mpi_comm_mct

  call tai_timcom_initialize(compid, mpi_comm_mct)
  call tai_timcom_run(compid)
  call tai_timcom_finalize

end subroutine tai_timcom_main

subroutine tai_timcom_initialize(compid, mpi_comm_mct)
  use tai_timcom_comm,    only: comm_init
  use tai_timcom_general, only: general_setup
  use tai_timcom_init, only: grid_init
  use tai_timcom_datm, only: datm_init
  use tai_timcom_dice, only: dice_init
  use tai_timcom_drof, only: drof_init
  implicit none

  integer, intent(in) :: compid, mpi_comm_mct
  character(len=256) :: grid_list(4)

  call comm_init(mpi_comm_mct)
  call general_setup(grid_list)
  call grid_init(grid_list(1), compid)
  call datm_init(grid_list(2))
  call dice_init(grid_list(3))
  call drof_init(grid_list(4))
  
  return
end subroutine tai_timcom_initialize

subroutine tai_timcom_run(compid)
  use tai_timcom_general, only: timer_drv, update_timer
  use tai_timcom_datm,    only: datm_exec, datm_to_timcom
  use tai_timcom_dice,    only: dice_exec, dice_to_timcom
  use tai_timcom_drof,    only: drof_exec, drof_to_timcom
  use tai_timcom_drv,     only: timcom_exec
  use tai_timcom_cplmct,  only: coupler_data_exch
  implicit none

  integer, intent(in) :: compid
  integer :: itf, it0, mxit

  it0  = timer_drv%it0
  mxit = timer_drv%mxit

  do itf = 1+it0, mxit
    call update_timer(itf, timer_drv)

    !call datm_exec(timer_drv)
    !call datm_to_timcom
    if(itf .eq. 1) call coupler_data_exch(compid, timer_drv)

    !call dice_exec(timer_drv)
    !call dice_to_timcom

    call drof_exec(timer_drv)
    call drof_to_timcom
  
    call timcom_exec(itf, timer_drv)

    call coupler_data_exch(compid, timer_drv)
  end do

  return
end subroutine tai_timcom_run

subroutine timcom_run_rsm_cpl(itf, igrd, jgrd, rsm_field, tsea, ssu, ssv)
  use tai_timcom_const,   only: r8
  use tai_timcom_general, only: timer_drv, update_timer
  use tai_timcom_drv,     only: timcom_exec
  use timcom_rsm_cpl, only: import_rsm_grid, export_rsm_grid
  implicit none

  integer :: itf, igrd, jgrd
  real(r8) :: rsm_field(igrd,jgrd,9), tsea(igrd,jgrd), ssu(igrd,jgrd), ssv(igrd,jgrd)

  call import_rsm_grid(itf, igrd, jgrd, rsm_field, tsea)

  call update_timer(itf, timer_drv)
  call timcom_exec(itf, timer_drv)

  call export_rsm_grid(igrd, jgrd, tsea, ssu, ssv)

end subroutine timcom_run_rsm_cpl

subroutine tai_timcom_finalize
  use tai_timcom_comm, only: comm_finalize_cpl
  implicit none

  call comm_finalize_cpl

end subroutine tai_timcom_finalize
