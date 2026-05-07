#ifdef cpl

subroutine timcom_main_gpu(compid, mpi_comm_mct)
  implicit none

  integer, intent(in) :: compid, mpi_comm_mct

  call timcom_initialize_gpu(compid, mpi_comm_mct)
  call timcom_run_gpu(compid)
  call timcom_finalize_gpu

end subroutine timcom_main_gpu

#else

subroutine timcom_main_gpu
  implicit none

  call timcom_initialize_gpu
  call timcom_run_gpu
  call timcom_finalize_gpu

end subroutine timcom_main_gpu

#endif

#ifdef cpl
subroutine timcom_initialize_gpu(compid, mpi_comm_mct)
#else
subroutine timcom_initialize_gpu
#endif

  use timcom_comm_gpu, only: timcom_comm_gpu_init
  implicit none

#ifdef cpl
  integer, intent(in) :: compid, mpi_comm_mct
  call timcom_initialize(compid, mpi_comm_mct)
#else
  call timcom_initialize
#endif

  call timcom_comm_gpu_init

  return
end subroutine timcom_initialize_gpu

subroutine timcom_finalize_gpu

  use timcom_comm_gpu, only: timcom_comm_gpu_finalize
  implicit none

  call timcom_finalize
  call timcom_comm_gpu_finalize

end subroutine timcom_finalize_gpu

#ifdef cpl
subroutine timcom_run_gpu(compid)
#else
subroutine timcom_run_gpu
#endif

  use timcom_general, only: timer_drv, update_timer, restart
  use timcom_datm,    only: datm_exec, datm_to_timcom
  use timcom_dice,    only: dice_exec, dice_to_timcom
  use timcom_drof,    only: drof_exec, drof_to_timcom
  use timcom_drv_gpu, only: timcom_exec_gpu, incomp_res_gpu, s_drv => s
  use timcom_drv,     only: timcom_check_restart
  use timcom_cplmct,  only: coupler_data_exch
#ifdef cpl_rocn
  use timcom_cplmct,  only: coupler_glb2rocn_exch
#endif
#ifdef cpl_cice
  use CICE_RunMod
  use timcom_cplmct,  only: coupler_cice_import, coupler_cice_export
#endif
  use hyperlink
  use timcom_solver_gpu, only: rhn_solver => rhn, rho_solver => rho, alpha_solver => alpha, beta_solver => beta, &
                               w_solver => w, tmp_solver => tmp, tp_solver => tp, &
                               res1_solver => res1, res2_solver => res2, res3_solver => res3, &
                               a0r_solver => a0r, s_solver => s, r_solver => r, q_solver => q, &
                               csalpha_solver => csalpha, csbeta_solver => csbeta, csy_solver => csy, &
                               csomega_solver => csomega, one_csy_solver => one_csy, rr_solver => rr
  use timcom_kpp, only: zgrid_kpp => zgrid, tidal_coef_kpp => tidal_coef, bckgrnd_vvc_kpp => bckgrnd_vvc, &
                        bckgrnd_vdc_kpp => bckgrnd_vdc, tidal_diff_kpp => tidal_diff, &
                        fcort_kpp => fcort, bolus_sp_kpp => bolus_sp, fstokes_kpp => fstokes, &
                        hmxl_kpp => hmxl, hwide_kpp => hwide, cg_kpp => cg
  use timcom_comm_gpu, only: timcom_comm_gpu_init, timcom_comm_gpu_finalize
  implicit none
#ifdef cpl
  integer, intent(in) :: compid
#endif
  integer :: itf, it0, mxit, async_id

  async_id = 1

  allocate(a0r_solver(0:nx+1,0:ny+1), s_solver(0:nx+1,0:ny+1), q_solver(0:nx+1,0:ny+1), r_solver(0:nx+1,0:ny+1,2))
  allocate(s_drv(nx,ny))

  !$acc enter data async(async_id) create(s_drv, rhn_solver, rho_solver, &
  !$acc& alpha_solver, beta_solver, w_solver, tmp_solver, tp_solver, &
  !$acc& res1_solver, res2_solver, res3_solver, &
  !$acc& a0r_solver, s_solver, r_solver, q_solver, csalpha_solver, csbeta_solver, csy_solver, &
  !$acc& csomega_solver, one_csy_solver, rr_solver, &
  !$acc& u, v, p0, iw, odx, ody, odz, ocs, csv, kb, w, x, iu, iv, in, &
  !$acc& ab, al, ac, ar, at, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh, &
  !$acc& nx, ny, nxf, nyf, nz, max_iter_p0, odyv, u_change, v_change, &
  !$acc& xu_bgn, xu_end, yv_bgn, yv_end, odt, nzf, &
  !$acc& incomp_res_gpu, area, total_area, u2, v2, t2, s2, qice, aqice, qflux, &
  !$acc& u1, ulf, v1, vlf, t1, tlf, s1, slf, dt, npx, npy, symm_np, myid_x, myid_y, peri_x, peri_y, &
  !$acc& tanphi, curv_f, t_nudge, t_da, kpp_src, smft, stf, sw_trans, shf_qsw, &
  !$acc& p, dmx, dmy, ev, dhx, dhy, hv, rho, z_grid, dzw, &
  !$acc& taux, tauy, senh, lath, lwup, evap, ifrc, vice, vsno, melt, melth, salt, u_10, v_10, t_10, q_10, pslv, &
  !$acc& lwdn, rain, snow, roff, swup, swdn, vdc, vvc, hbk, vbk, odzw, &
  !$acc& zgrid_kpp, tidal_coef_kpp, bckgrnd_vvc_kpp, bckgrnd_vdc_kpp, tidal_diff_kpp, &
  !$acc& fcort_kpp, bolus_sp_kpp, fstokes_kpp, hmxl_kpp, kpp_hblt, hwide_kpp, cg_kpp, z_face)

  !$acc wait(async_id)

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
   
    call timcom_exec_gpu(itf, timer_drv)
#ifdef cpl
    call coupler_data_exch(compid, timer_drv)
#endif
    call timcom_check_restart(itf, timer_drv)

#ifdef cpl_rocn
    if(mod(itf,3).eq.0) call coupler_glb2rocn_exch(compid, timer_drv)
#endif

  end do

  !$acc exit data async(async_id) delete(s_drv, rhn_solver, rho_solver, &
  !$acc& alpha_solver, beta_solver, w_solver, tmp_solver, tp_solver, &
  !$acc& res1_solver, res2_solver, res3_solver, &
  !$acc& a0r_solver, s_solver, r_solver, q_solver, csalpha_solver, csbeta_solver, csy_solver, &
  !$acc& csomega_solver, one_csy_solver, rr_solver, &
  !$acc& u, v, p0, iw, odx, ody, odz, ocs, csv, kb, w, x, iu, iv, in, &
  !$acc& ab, al, ac, ar, at, cgr, cgrh, cgp, cgv, cgs, cgt, cgph, cgsh, &
  !$acc& nx, ny, nxf, nyf, nz, max_iter_p0, odyv, u_change, v_change, &
  !$acc& xu_bgn, xu_end, yv_bgn, yv_end, odt, nzf, &
  !$acc& incomp_res_gpu, area, total_area, u2, v2, t2, s2, qice, aqice, qflux, &
  !$acc& u1, ulf, v1, vlf, t1, tlf, s1, slf, dt, npx, npy, symm_np, myid_x, myid_y, peri_x, peri_y, &
  !$acc& tanphi, curv_f, t_nudge, t_da, kpp_src, smft, stf, sw_trans, shf_qsw, &
  !$acc& p, dmx, dmy, ev, dhx, dhy, hv, rho, z_grid, dzw, &
  !$acc& taux, tauy, senh, lath, lwup, evap, ifrc, vice, vsno, melt, melth, salt, u_10, v_10, t_10, q_10, pslv, &
  !$acc& lwdn, rain, snow, roff, swup, swdn, vdc, vvc, hbk, vbk, odzw, &
  !$acc& zgrid_kpp, tidal_coef_kpp, bckgrnd_vvc_kpp, bckgrnd_vdc_kpp, tidal_diff_kpp, &
  !$acc& fcort_kpp, bolus_sp_kpp, fstokes_kpp, hmxl_kpp, kpp_hblt, hwide_kpp, cg_kpp, z_face)

  !$acc wait(async_id)

  return
end subroutine timcom_run_gpu
