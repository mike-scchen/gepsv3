!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_wind_mixing
  implicit none

  call timcom_initialize_gpu

  write (*, *) '=== Testing KPP option (opt_windmix = 1) ==='
  call wind_mixing_unit(1)

  write (*, *) '=== Testing PP82 option (opt_windmix = 2) ==='
  call wind_mixing_unit(2)

  call timcom_finalize_gpu

end program test_wind_mixing

subroutine wind_mixing_unit(test_opt_windmix)
  use timcom_drv, only: wind_mixing
  use timcom_drv_gpu, only: wind_mixing_gpu
  use timcom_kpp, only: zgrid, tidal_coef, bckgrnd_vvc, bckgrnd_vdc, tidal_diff, &
                        fcort, bolus_sp, fstokes, hmxl, hwide, cg
  use hyperlink, only: nx, ny, nz, w, u2, v2, t2, s2, rho, &
                       sw_trans, stf, shf_qsw, smft, vbk, hbk, add, vdc, &
                       vvc, ev, hv, kpp_src, kpp_hblt, opt_windmix, z_grid, z_face, odz, odzw, iw, kb, nzf, in, dt
  implicit none

  integer, intent(in) :: test_opt_windmix

  real(8), dimension(:), allocatable :: w_cpu, u2_cpu, v2_cpu, rho_cpu
  real(8), dimension(:), allocatable :: vbk_cpu, hbk_cpu, add_cpu, vdc_cpu
  real(8), dimension(:), allocatable :: vvc_cpu, ev_cpu, hv_cpu, kpp_src_cpu
  real(8), dimension(:), allocatable :: kpp_hblt_cpu, t2_cpu, s2_cpu
  real(8), dimension(:), allocatable :: sw_trans_cpu, stf_cpu, shf_qsw_cpu, smft_cpu

  real(8), dimension(:), allocatable :: w_gpu, u2_gpu, v2_gpu, rho_gpu
  real(8), dimension(:), allocatable :: vbk_gpu, hbk_gpu, add_gpu, vdc_gpu
  real(8), dimension(:), allocatable :: vvc_gpu, ev_gpu, hv_gpu, kpp_src_gpu
  real(8), dimension(:), allocatable :: kpp_hblt_gpu, t2_gpu, s2_gpu
  real(8), dimension(:), allocatable :: sw_trans_gpu, stf_gpu, shf_qsw_gpu, smft_gpu

  real(8), dimension(:,:,:), allocatable :: w_initial, u2_initial, v2_initial, rho_initial
  real(8), dimension(:), allocatable :: vbk_initial, hbk_initial
  real(8), dimension(:,:,:), allocatable :: add_initial, ev_initial
  real(8), dimension(:,:,:,:), allocatable :: hv_initial, vdc_initial, kpp_src_initial
  real(8), dimension(:,:,:), allocatable :: vvc_initial, t2_initial, s2_initial
  real(8), dimension(:,:,:), allocatable :: sw_trans_initial
  real(8), dimension(:,:), allocatable :: kpp_hblt_initial
  real(8), dimension(:,:,:), allocatable :: stf_initial, smft_initial
  real(8), dimension(:,:), allocatable :: shf_qsw_initial

  integer, parameter :: &
    opt_kpp  = 1, &
    opt_pp82 = 2
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i, itf

  ! Set the opt_windmix to the test value
  opt_windmix = test_opt_windmix

  ! Initialize variables
  allocate(w_initial(size(w,1), size(w,2), size(w,3)))
  allocate(u2_initial(size(u2,1), size(u2,2), size(u2,3)))
  allocate(v2_initial(size(v2,1), size(v2,2), size(v2,3)))
  allocate(t2_initial(size(t2,1), size(t2,2), size(t2,3)))
  allocate(s2_initial(size(s2,1), size(s2,2), size(s2,3)))
  allocate(rho_initial(size(rho,1), size(rho,2), size(rho,3)))
  allocate(sw_trans_initial(size(sw_trans,1), size(sw_trans,2), size(sw_trans,3)))
  allocate(stf_initial(size(stf,1), size(stf,2), size(stf,3)))
  allocate(smft_initial(size(smft,1), size(smft,2), size(smft,3)))
  allocate(vbk_initial(size(vbk)))
  allocate(hbk_initial(size(hbk)))
  allocate(add_initial(size(add,1), size(add,2), size(add,3)))
  allocate(vdc_initial(size(vdc,1), size(vdc,2), size(vdc,3), size(vdc,4)))
  allocate(vvc_initial(size(vvc,1), size(vvc,2), size(vvc,3)))
  allocate(ev_initial(size(ev,1), size(ev,2), size(ev,3)))
  allocate(hv_initial(size(hv,1), size(hv,2), size(hv,3), size(hv,4)))
  allocate(kpp_src_initial(size(kpp_src,1), size(kpp_src,2), size(kpp_src,3), size(kpp_src,4)))
  allocate(kpp_hblt_initial(size(kpp_hblt,1), size(kpp_hblt,2)))

  allocate(w_cpu(size(w)), u2_cpu(size(u2)), v2_cpu(size(v2)), rho_cpu(size(rho)))
  allocate(vbk_cpu(size(vbk)), hbk_cpu(size(hbk)), add_cpu(size(add)), vdc_cpu(size(vdc)))
  allocate(vvc_cpu(size(vvc)), ev_cpu(size(ev)), hv_cpu(size(hv)), kpp_src_cpu(size(kpp_src)))
  allocate(kpp_hblt_cpu(size(kpp_hblt)), t2_cpu(size(t2)), s2_cpu(size(s2)))
  allocate(sw_trans_cpu(size(sw_trans)), stf_cpu(size(stf)))
  allocate(shf_qsw_cpu(size(shf_qsw)), smft_cpu(size(smft)))

  allocate(w_gpu(size(w)), u2_gpu(size(u2)), v2_gpu(size(v2)), rho_gpu(size(rho)))
  allocate(vbk_gpu(size(vbk)), hbk_gpu(size(hbk)), add_gpu(size(add)), vdc_gpu(size(vdc)))
  allocate(vvc_gpu(size(vvc)), ev_gpu(size(ev)), hv_gpu(size(hv)), kpp_src_gpu(size(kpp_src)))
  allocate(kpp_hblt_gpu(size(kpp_hblt)), t2_gpu(size(t2)), s2_gpu(size(s2)))
  allocate(sw_trans_gpu(size(sw_trans)), stf_gpu(size(stf)))
  allocate(shf_qsw_gpu(size(shf_qsw)), smft_gpu(size(smft)))

  ! Initialize with random values
  call random_number(w)
  call random_number(u2)
  call random_number(v2)
  call random_number(t2)
  call random_number(s2)
  call random_number(rho)
  call random_number(sw_trans)
  call random_number(stf)
  call random_number(shf_qsw)
  call random_number(smft)
  call random_number(vbk)
  call random_number(hbk)
  call random_number(add)
  call random_number(vdc)
  call random_number(vvc)
  call random_number(ev)
  call random_number(hv)
  call random_number(kpp_src)
  call random_number(kpp_hblt)

  ! Store initial values
  w_initial = w
  u2_initial = u2
  v2_initial = v2
  t2_initial = t2
  s2_initial = s2
  rho_initial = rho
  sw_trans_initial = sw_trans
  stf_initial = stf
  smft_initial = smft
  vbk_initial = vbk
  hbk_initial = hbk
  add_initial = add
  vdc_initial = vdc
  vvc_initial = vvc
  ev_initial = ev
  hv_initial = hv
  kpp_src_initial = kpp_src
  kpp_hblt_initial = kpp_hblt
  shf_qsw_initial = shf_qsw

  ! Set itf for non-coupled case
#ifndef cpl
  itf = 1
#endif

  ! CPU execution
  do i = 1, num_steps
#ifdef cpl
    call wind_mixing
#else
    call wind_mixing(itf)
#endif
  end do

  ! Store CPU results
  w_cpu = reshape(w, [size(w)])
  u2_cpu = reshape(u2, [size(u2)])
  v2_cpu = reshape(v2, [size(v2)])
  t2_cpu = reshape(t2, [size(t2)])
  s2_cpu = reshape(s2, [size(s2)])
  rho_cpu = reshape(rho, [size(rho)])
  sw_trans_cpu = reshape(sw_trans, [size(sw_trans)])
  stf_cpu = reshape(stf, [size(stf)])
  smft_cpu = reshape(smft, [size(smft)])
  vbk_cpu = reshape(vbk, [size(vbk)])
  hbk_cpu = reshape(hbk, [size(hbk)])
  add_cpu = reshape(add, [size(add)])
  vdc_cpu = reshape(vdc, [size(vdc)])
  vvc_cpu = reshape(vvc, [size(vvc)])
  ev_cpu = reshape(ev, [size(ev)])
  hv_cpu = reshape(hv, [size(hv)])
  kpp_src_cpu = reshape(kpp_src, [size(kpp_src)])
  kpp_hblt_cpu = reshape(kpp_hblt, [size(kpp_hblt)])
  shf_qsw_cpu = reshape(shf_qsw, [size(shf_qsw)])

  ! Restore initial values before GPU execution
  w = w_initial
  u2 = u2_initial
  v2 = v2_initial
  t2 = t2_initial
  s2 = s2_initial
  rho = rho_initial
  sw_trans = sw_trans_initial
  stf = stf_initial
  smft = smft_initial
  vbk = vbk_initial
  hbk = hbk_initial
  add = add_initial
  vdc = vdc_initial
  vvc = vvc_initial
  ev = ev_initial
  hv = hv_initial
  kpp_src = kpp_src_initial
  kpp_hblt = kpp_hblt_initial
  shf_qsw = shf_qsw_initial

  ! GPU execution

  select case(opt_windmix)
  case(opt_kpp)
    !$acc enter data async(async_id) &
    !$acc& copyin(t2, s2, hv, ev, vdc, vvc, hbk, vbk, odzw, iw, ny, nx, nz, u2, v2, &
    !$acc&        zgrid, kb, tidal_coef, bckgrnd_vvc, bckgrnd_vdc, nzf, tidal_diff, &
    !$acc&        fcort, bolus_sp, fstokes, sw_trans, in, hmxl, kpp_hblt, hwide, cg, &
    !$acc&        stf, odz, z_grid, z_face, kpp_src, shf_qsw, smft, dt)
  case(opt_pp82)
    !$acc enter data async(async_id) &
    !$acc& copyin(nx, ny, nz, w, u2, v2, rho, hbk, odzw, dt, iw) &
    !$acc& create(ev, hv)
  end select

  do i = 1, num_steps
#ifdef cpl
    call wind_mixing_gpu(async_id)
#else
    call wind_mixing_gpu(itf, async_id)
#endif
  end do

  select case(opt_windmix)
  case(opt_kpp)
    !$acc exit data async(async_id) &
    !$acc& copyout(t2, s2, hv, ev, vdc, vvc, hbk, vbk, odzw, iw, ny, nx, nz, u2, v2, &
    !$acc&         zgrid, kb, tidal_coef, bckgrnd_vvc, bckgrnd_vdc, nzf, tidal_diff, &
    !$acc&         fcort, bolus_sp, fstokes, sw_trans, in, hmxl, kpp_hblt, hwide, cg, &
    !$acc&         stf, odz, z_grid, z_face, kpp_src, shf_qsw, smft, dt)
  case(opt_pp82)
    !$acc exit data async(async_id) &
    !$acc& copyout(ev, hv) &
    !$acc& delete(nx, ny, nz, w, u2, v2, rho, hbk, odzw, dt, iw)
  end select
  !$acc wait(async_id)

  ! Store GPU results
  w_gpu = reshape(w, [size(w)])
  u2_gpu = reshape(u2, [size(u2)])
  v2_gpu = reshape(v2, [size(v2)])
  t2_gpu = reshape(t2, [size(t2)])
  s2_gpu = reshape(s2, [size(s2)])
  rho_gpu = reshape(rho, [size(rho)])
  sw_trans_gpu = reshape(sw_trans, [size(sw_trans)])
  stf_gpu = reshape(stf, [size(stf)])
  smft_gpu = reshape(smft, [size(smft)])
  vbk_gpu = reshape(vbk, [size(vbk)])
  hbk_gpu = reshape(hbk, [size(hbk)])
  add_gpu = reshape(add, [size(add)])
  vdc_gpu = reshape(vdc, [size(vdc)])
  vvc_gpu = reshape(vvc, [size(vvc)])
  ev_gpu = reshape(ev, [size(ev)])
  hv_gpu = reshape(hv, [size(hv)])
  kpp_src_gpu = reshape(kpp_src, [size(kpp_src)])
  kpp_hblt_gpu = reshape(kpp_hblt, [size(kpp_hblt)])
  shf_qsw_gpu = reshape(shf_qsw, [size(shf_qsw)])

  ! Assertions based on the option being tested
  if (test_opt_windmix == 1) then
    ! Assertions for KPP option
    call assert_allclose(w_gpu, size(w_gpu), w_cpu, size(w_cpu), 1e-10, 1e-10, "Array w (KPP)")
    call assert_allclose(u2_gpu, size(u2_gpu), u2_cpu, size(u2_cpu), 1e-10, 1e-10, "Array u2 (KPP)")
    call assert_allclose(v2_gpu, size(v2_gpu), v2_cpu, size(v2_cpu), 1e-10, 1e-10, "Array v2 (KPP)")
    call assert_allclose(t2_gpu, size(t2_gpu), t2_cpu, size(t2_cpu), 1e-10, 1e-10, "Array t2 (KPP)")
    call assert_allclose(s2_gpu, size(s2_gpu), s2_cpu, size(s2_cpu), 1e-10, 1e-10, "Array s2 (KPP)")
    call assert_allclose(rho_gpu, size(rho_gpu), rho_cpu, size(rho_cpu), 1e-10, 1e-10, "Array rho (KPP)")
    call assert_allclose(sw_trans_gpu, size(sw_trans_gpu), sw_trans_cpu, size(sw_trans_cpu), 1e-10, 1e-10, "Array sw_trans (KPP)")
    call assert_allclose(stf_gpu, size(stf_gpu), stf_cpu, size(stf_cpu), 1e-10, 1e-10, "Array stf (KPP)")
    call assert_allclose(shf_qsw_gpu, size(shf_qsw_gpu), shf_qsw_cpu, size(shf_qsw_cpu), 1e-10, 1e-10, "Array shf_qsw (KPP)")
    call assert_allclose(smft_gpu, size(smft_gpu), smft_cpu, size(smft_cpu), 1e-10, 1e-10, "Array smft (KPP)")
    call assert_allclose(vbk_gpu, size(vbk_gpu), vbk_cpu, size(vbk_cpu), 1e-10, 1e-10, "Array vbk (KPP)")
    call assert_allclose(hbk_gpu, size(hbk_gpu), hbk_cpu, size(hbk_cpu), 1e-10, 1e-10, "Array hbk (KPP)")
    call assert_allclose(add_gpu, size(add_gpu), add_cpu, size(add_cpu), 1e-10, 1e-10, "Array add (KPP)")
    call assert_allclose(vdc_gpu, size(vdc_gpu), vdc_cpu, size(vdc_cpu), 1e-10, 1e-10, "Array vdc (KPP)")
    call assert_allclose(vvc_gpu, size(vvc_gpu), vvc_cpu, size(vvc_cpu), 1e-10, 1e-10, "Array vvc (KPP)")
    call assert_allclose(ev_gpu, size(ev_gpu), ev_cpu, size(ev_cpu), 1e-10, 1e-10, "Array ev (KPP)")
    call assert_allclose(hv_gpu, size(hv_gpu), hv_cpu, size(hv_cpu), 1e-10, 1e-10, "Array hv (KPP)")
    call assert_allclose(kpp_src_gpu, size(kpp_src_gpu), kpp_src_cpu, size(kpp_src_cpu), 1e-10, 1e-10, "Array kpp_src (KPP)")
    call assert_allclose(kpp_hblt_gpu, size(kpp_hblt_gpu), kpp_hblt_cpu, size(kpp_hblt_cpu), 1e-10, 1e-10, "Array kpp_hblt (KPP)")

  else if (test_opt_windmix == 2) then
    ! Assertions for PP82 option
    call assert_allclose(w_gpu, size(w_gpu), w_cpu, size(w_cpu), 1e-10, 1e-10, "Array w (PP82)")
    call assert_allclose(u2_gpu, size(u2_gpu), u2_cpu, size(u2_cpu), 1e-10, 1e-10, "Array u2 (PP82)")
    call assert_allclose(v2_gpu, size(v2_gpu), v2_cpu, size(v2_cpu), 1e-10, 1e-10, "Array v2 (PP82)")
    call assert_allclose(rho_gpu, size(rho_gpu), rho_cpu, size(rho_cpu), 1e-10, 1e-10, "Array rho (PP82)")
    call assert_allclose(hbk_gpu, size(hbk_gpu), hbk_cpu, size(hbk_cpu), 1e-10, 1e-10, "Array hbk (PP82)")
    call assert_allclose(ev_gpu, size(ev_gpu), ev_cpu, size(ev_cpu), 1e-10, 1e-10, "Array ev (PP82)")
    call assert_allclose(hv_gpu, size(hv_gpu), hv_cpu, size(hv_cpu), 1e-10, 1e-10, "Array hv (PP82)")
  end if

  write (*, *) 'Test completed'

end subroutine wind_mixing_unit
