!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_calu_surf_flux
  implicit none

  call timcom_initialize_gpu
  call calu_surf_flux_unit
  call timcom_finalize_gpu

end program test_calu_surf_flux

subroutine calu_surf_flux_unit
  use timcom_drv, only: calu_surf_flux
  use timcom_drv_gpu, only: calu_surf_flux_gpu
  use timcom_phy, only: opt_bulk, opt_bulk_cesm, opt_bulk_ncar
  use hyperlink, only: in, nx, ny, u2, v2, t2, s2, rho, &
                       u_10, v_10, t_10, q_10, pslv, &
                       swup, swdn, lwup, lwdn, taux, &
                       tauy, senh, lath, evap, rain, &
                       snow, roff, ifrc, &
                       stf, shf_qsw, smft, &
                       m_comm_cart, r8type3d, &
                       nbid, ndim, symm_np, hflux_factor
  implicit none

  real(8), dimension(:), allocatable :: stf_cpu, shf_qsw_cpu, smft_cpu
  real(8), dimension(:), allocatable :: taux_cpu, tauy_cpu, senh_cpu, lath_cpu, lwup_cpu, evap_cpu
  real(8), dimension(:), allocatable :: stf_gpu, shf_qsw_gpu, smft_gpu
  real(8), dimension(:), allocatable :: taux_gpu, tauy_gpu, senh_gpu, lath_gpu, lwup_gpu, evap_gpu
  real(8), dimension(:,:,:), allocatable :: stf_initial, smft_initial
  real(8), dimension(:,:), allocatable :: shf_qsw_initial
  real(8), dimension(:,:), allocatable :: taux_initial, tauy_initial, senh_initial, lath_initial
  real(8), dimension(:,:), allocatable :: lwup_initial, evap_initial
  real(8), dimension(:,:,:), allocatable :: in_r8
  real(8), dimension(:,:), allocatable :: u2_initial, v2_initial, t2_initial
  real(8), dimension(:,:), allocatable :: u_10_initial, v_10_initial, t_10_initial, q_10_initial, pslv_initial
  real(8), dimension(:,:), allocatable :: swup_initial, swdn_initial, lwdn_initial
  real(8), dimension(:,:), allocatable :: rain_initial, snow_initial, roff_initial, ifrc_initial
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  ! opt_bulk = opt_bulk_cesm

  ! Initialize variables
  allocate(stf_initial(size(stf,1), size(stf,2), size(stf,3)), smft_initial(size(smft,1), size(smft,2), size(smft,3)))
  allocate(shf_qsw_initial(size(shf_qsw,1), size(shf_qsw,2)))
  allocate(taux_initial(size(taux,1), size(taux,2)), tauy_initial(size(tauy,1), size(tauy,2)))
  allocate(senh_initial(size(senh,1), size(senh,2)), lath_initial(size(lath,1), size(lath,2)))
  allocate(lwup_initial(size(lwup,1), size(lwup,2)), evap_initial(size(evap,1), size(evap,2)))

  allocate(stf_cpu(size(stf)), shf_qsw_cpu(size(shf_qsw)), smft_cpu(size(smft)))
  allocate(taux_cpu(size(taux)), tauy_cpu(size(tauy)), senh_cpu(size(senh)))
  allocate(lath_cpu(size(lath)), lwup_cpu(size(lwup)), evap_cpu(size(evap)))
  allocate(stf_gpu(size(stf)), shf_qsw_gpu(size(shf_qsw)), smft_gpu(size(smft)))
  allocate(taux_gpu(size(taux)), tauy_gpu(size(tauy)), senh_gpu(size(senh)))
  allocate(lath_gpu(size(lath)), lwup_gpu(size(lwup)), evap_gpu(size(evap)))

  allocate(in_r8(size(in,1), size(in,2), size(in,3)))
  allocate(u2_initial(size(u2,1), size(u2,2)), v2_initial(size(v2,1), size(v2,2)), t2_initial(size(t2,1), size(t2,2)))
  allocate(u_10_initial(size(u_10,1), size(u_10,2)), v_10_initial(size(v_10,1), size(v_10,2)))
  allocate(t_10_initial(size(t_10,1), size(t_10,2)), q_10_initial(size(q_10,1), size(q_10,2)))
  allocate(pslv_initial(size(pslv,1), size(pslv,2)))
  allocate(swup_initial(size(swup,1), size(swup,2)), swdn_initial(size(swdn,1), size(swdn,2)))
  allocate(lwdn_initial(size(lwdn,1), size(lwdn,2)))
  allocate(rain_initial(size(rain,1), size(rain,2)))
  allocate(snow_initial(size(snow,1), size(snow,2)), roff_initial(size(roff,1), size(roff,2)))
  allocate(ifrc_initial(size(ifrc,1), size(ifrc,2)))

  ! Initialize input parameters with random values
  call random_number(in_r8)
  in = nint(in_r8)
  call random_number(u2_initial)
  u2(1:nx,1:ny,1) = u2_initial * 100.0d0
  call random_number(v2_initial)
  v2(1:nx,1:ny,1) = v2_initial * 100.0d0
  call random_number(t2_initial)
  t2(1:nx,1:ny,1) = t2_initial * 30.0d0
  call random_number(u_10_initial)
  u_10 = u_10_initial * 10.0d0
  call random_number(v_10_initial)
  v_10 = v_10_initial * 10.0d0
  call random_number(t_10_initial)
  t_10 = t_10_initial * 30.0d0
  call random_number(q_10_initial)
  q_10 = q_10_initial * 0.02d0
  call random_number(pslv_initial)
  pslv = pslv_initial * 1000.0d0 + 100000.0d0
  call random_number(swup_initial)
  swup = swup_initial * 200.0d0
  call random_number(swdn_initial)
  swdn = swdn_initial * 1000.0d0
  call random_number(lwdn_initial)
  lwdn = lwdn_initial * 400.0d0
  call random_number(rain_initial)
  rain = rain_initial * 0.001d0
  call random_number(snow_initial)
  snow = snow_initial * 0.0001d0
  call random_number(roff_initial)
  roff = roff_initial * 0.0001d0
  call random_number(ifrc_initial)
  ifrc = ifrc_initial * 0.5d0

  ! Store initial values of output arrays
  stf_initial = stf
  smft_initial = smft
  shf_qsw_initial = shf_qsw
  taux_initial = taux
  tauy_initial = tauy
  senh_initial = senh
  lath_initial = lath
  lwup_initial = lwup
  evap_initial = evap

  ! CPU execution
  do i = 1, num_steps
    call calu_surf_flux
  end do
  stf_cpu = reshape(stf, [size(stf)])
  smft_cpu = reshape(smft, [size(smft)])
  shf_qsw_cpu = reshape(shf_qsw, [size(shf_qsw)])
  taux_cpu = reshape(taux, [size(taux)])
  tauy_cpu = reshape(tauy, [size(tauy)])
  senh_cpu = reshape(senh, [size(senh)])
  lath_cpu = reshape(lath, [size(lath)])
  lwup_cpu = reshape(lwup, [size(lwup)])
  evap_cpu = reshape(evap, [size(evap)])

  ! Restore initial values before GPU execution
  stf = stf_initial
  smft = smft_initial
  shf_qsw = shf_qsw_initial
  taux = taux_initial
  tauy = tauy_initial
  senh = senh_initial
  lath = lath_initial
  lwup = lwup_initial
  evap = evap_initial

  ! GPU execution

  !$acc enter data async(async_id) &
  !$acc& create(taux, tauy, senh, lath, lwup, evap) &
  !$acc& copyin(ifrc, u_10, v_10, t_10, q_10, pslv) &
  !$acc& copyin(u2, v2, t2, in, nx, ny) &
  !$acc& create(stf, shf_qsw, smft) &
  !$acc& copyin(lwdn, rain, snow, roff, swup, swdn)
  do i = 1, num_steps
    call calu_surf_flux_gpu(async_id)
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(stf, shf_qsw, smft, taux, tauy, senh, lath, lwup, evap) &
  !$acc& delete(u2, v2, t2) &
  !$acc& delete(u_10, v_10, t_10, q_10, pslv) &
  !$acc& delete(ifrc, lwdn, rain, snow, roff, swup, swdn) &
  !$acc& delete(in, nx, ny)
  !$acc wait(async_id)

  stf_gpu = reshape(stf, [size(stf)])
  smft_gpu = reshape(smft, [size(smft)])
  shf_qsw_gpu = reshape(shf_qsw, [size(shf_qsw)])
  taux_gpu = reshape(taux, [size(taux)])
  tauy_gpu = reshape(tauy, [size(tauy)])
  senh_gpu = reshape(senh, [size(senh)])
  lath_gpu = reshape(lath, [size(lath)])
  lwup_gpu = reshape(lwup, [size(lwup)])
  evap_gpu = reshape(evap, [size(evap)])

  ! Assertions
  call assert_allclose(stf_gpu, size(stf_gpu), stf_cpu, size(stf_cpu), 1e-10, 1e-10, "Array stf")
  call assert_allclose(smft_gpu, size(smft_gpu), smft_cpu, size(smft_cpu), 1e-10, 1e-10, "Array smft")
  call assert_allclose(shf_qsw_gpu, size(shf_qsw_gpu), shf_qsw_cpu, size(shf_qsw_cpu), 1e-10, 1e-10, "Array shf_qsw")
  call assert_allclose(taux_gpu, size(taux_gpu), taux_cpu, size(taux_cpu), 1e-10, 1e-10, "Array taux")
  call assert_allclose(tauy_gpu, size(tauy_gpu), tauy_cpu, size(tauy_cpu), 1e-10, 1e-10, "Array tauy")
  call assert_allclose(senh_gpu, size(senh_gpu), senh_cpu, size(senh_cpu), 1e-10, 1e-10, "Array senh")
  call assert_allclose(lath_gpu, size(lath_gpu), lath_cpu, size(lath_cpu), 1e-10, 1e-10, "Array lath")
  call assert_allclose(lwup_gpu, size(lwup_gpu), lwup_cpu, size(lwup_cpu), 1e-10, 1e-10, "Array lwup")
  call assert_allclose(evap_gpu, size(evap_gpu), evap_cpu, size(evap_cpu), 1e-10, 1e-10, "Array evap")

  write (*, *) 'Test complete'

  ! Deallocate arrays
  deallocate(stf_initial, smft_initial, shf_qsw_initial)
  deallocate(taux_initial, tauy_initial, senh_initial, lath_initial, lwup_initial, evap_initial)
  deallocate(stf_cpu, shf_qsw_cpu, smft_cpu)
  deallocate(taux_cpu, tauy_cpu, senh_cpu, lath_cpu, lwup_cpu, evap_cpu)
  deallocate(stf_gpu, shf_qsw_gpu, smft_gpu)
  deallocate(taux_gpu, tauy_gpu, senh_gpu, lath_gpu, lwup_gpu, evap_gpu)
  deallocate(in_r8)
  deallocate(u2_initial, v2_initial, t2_initial)
  deallocate(u_10_initial, v_10_initial, t_10_initial, q_10_initial, pslv_initial)
  deallocate(swup_initial, swdn_initial, lwdn_initial)
  deallocate(rain_initial, snow_initial, roff_initial, ifrc_initial)

end subroutine calu_surf_flux_unit
