#ifdef cpl

subroutine timcom_main(compid, mpi_comm_mct)
  implicit none

  integer, intent(in) :: compid, mpi_comm_mct

  call timcom_initialize(compid, mpi_comm_mct)
  call timcom_run(compid)
  call timcom_finalize

end subroutine timcom_main

#else

subroutine timcom_main
  implicit none

  call timcom_initialize
  call timcom_run
  call timcom_finalize

end subroutine timcom_main

#endif

#ifdef cpl
subroutine timcom_initialize(compid, mpi_comm_mct)
#else
subroutine timcom_initialize
#endif

  use timcom_comm,    only: comm_init
  use timcom_general, only: general_setup
  use timcom_init, only: grid_init
  use timcom_datm, only: datm_init
  use timcom_dice, only: dice_init
#ifdef cpl_cice
  use CICE_InitMod, only: CICE_Initialize
#endif
  use timcom_drof, only: drof_init
  implicit none
#ifdef cpl
  integer, intent(in) :: compid, mpi_comm_mct
#endif
  character(len=256) :: grid_list(4)

#ifdef cpl
  call comm_init(mpi_comm_mct)
#else
  call comm_init
#endif
  call general_setup(grid_list)

#ifdef cpl
  call grid_init(grid_list(1), compid)
#else
  call grid_init(grid_list(1))
#endif

  call datm_init(grid_list(2))
  call dice_init(grid_list(3))
  call drof_init(grid_list(4))
#ifdef cpl_cice
  call CICE_Initialize(mpi_comm_mct)
#endif  
  return
end subroutine timcom_initialize

#ifdef cpl
subroutine timcom_run(compid)
#else
subroutine timcom_run
#endif

  use timcom_general, only: timer_drv, update_timer, restart
  use timcom_datm,    only: datm_exec, datm_to_timcom
  use timcom_dice,    only: dice_exec, dice_to_timcom
  use timcom_drof,    only: drof_exec, drof_to_timcom
  use timcom_drv,     only: timcom_exec, timcom_check_restart
  use timcom_cplmct,  only: coupler_data_exch
#ifdef cpl_rocn
  use timcom_cplmct,  only: coupler_glb2rocn_exch
#endif
#ifdef cpl_cice
  use CICE_RunMod
  use timcom_cplmct,  only: coupler_cice_import, coupler_cice_export
#endif
  implicit none
#ifdef cpl
  integer, intent(in) :: compid
#endif
  integer :: itf, it0, mxit

  it0  = timer_drv%it0
  mxit = timer_drv%mxit

  do itf = 1+it0, mxit
    call update_timer(itf, timer_drv)

    call dice_exec(timer_drv)
    call dice_to_timcom
#ifdef cpl
    if(itf .eq. 1 .and. .not. restart) call coupler_data_exch(compid, timer_drv)
#else
    call datm_exec(timer_drv)
    call datm_to_timcom
#endif

#ifdef cpl_cice
    call coupler_cice_import(timer_drv)
    call CICE_Run
    call coupler_cice_export(timer_drv)
#endif
    call drof_exec(timer_drv)
    call drof_to_timcom
   
    call timcom_exec(itf, timer_drv)
#ifdef cpl
    call coupler_data_exch(compid, timer_drv)
#endif
    call timcom_check_restart(itf, timer_drv)

#ifdef cpl_rocn
    if(mod(itf,3).eq.0) call coupler_glb2rocn_exch(compid, timer_drv)
#endif

  end do

  return
end subroutine timcom_run

subroutine timcom_finalize
#ifdef cpl
  use timcom_comm, only: comm_finalize_cpl
#else
  use timcom_comm, only: comm_finalize
#endif

  implicit none

#ifdef cpl
  call comm_finalize_cpl
#else
  call comm_finalize
#endif

end subroutine timcom_finalize
