!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_solve_dynamic_equation
  implicit none

  call timcom_initialize_gpu
  call solve_dynamic_equation_unit
  call timcom_finalize_gpu

end program test_solve_dynamic_equation

subroutine solve_dynamic_equation_unit
  use timcom_drv, only: solve_dynamic_equation
  use timcom_drv_gpu, only: solve_dynamic_equation_gpu
  use hyperlink, only: u, v, w, p, &
                       u1, ulf, u2, &
                       v1, vlf, v2, &
                       t1, tlf, t2, &
                       s1, slf, s2, &
                       dmx, dmy, ev, &
                       dhx, dhy, hv, &
                       nx, ny, nz, &
                       in, iv, iw, iu, &
                       dt, odz, csv, odyv, ocs, ody, odx, &
                       xu_bgn, xu_end, yv_bgn, yv_end, nxf, nyf
  implicit none

  ! Arrays to store CPU and GPU results
  real(8), dimension(:), allocatable :: u2_cpu, u2_gpu
  real(8), dimension(:), allocatable :: v2_cpu, v2_gpu
  real(8), dimension(:), allocatable :: t2_cpu, t2_gpu
  real(8), dimension(:), allocatable :: s2_cpu, s2_gpu

  ! Arrays to store initial values
  real(8), dimension(:,:,:), allocatable :: u_initial, v_initial, w_initial, p_initial
  real(8), dimension(:,:,:), allocatable :: u1_initial, ulf_initial, u2_initial
  real(8), dimension(:,:,:), allocatable :: v1_initial, vlf_initial, v2_initial
  real(8), dimension(:,:,:), allocatable :: t1_initial, tlf_initial, t2_initial
  real(8), dimension(:,:,:), allocatable :: s1_initial, slf_initial, s2_initial
  real(8), dimension(:,:,:), allocatable :: dmx_initial, dmy_initial
  real(8), dimension(:,:,:), allocatable :: dhx_initial, dhy_initial
  real(8), dimension(:,:,:), allocatable :: ev_initial
  real(8), dimension(:,:,:,:), allocatable :: hv_initial

  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  ! Allocate memory for arrays
  allocate(u_initial(size(u,1), size(u,2), size(u,3)))
  allocate(v_initial(size(v,1), size(v,2), size(v,3)))
  allocate(w_initial(size(w,1), size(w,2), size(w,3)))
  allocate(p_initial(size(p,1), size(p,2), size(p,3)))

  allocate(u1_initial(size(u1,1), size(u1,2), size(u1,3)))
  allocate(ulf_initial(size(ulf,1), size(ulf,2), size(ulf,3)))
  allocate(u2_initial(size(u2,1), size(u2,2), size(u2,3)))

  allocate(v1_initial(size(v1,1), size(v1,2), size(v1,3)))
  allocate(vlf_initial(size(vlf,1), size(vlf,2), size(vlf,3)))
  allocate(v2_initial(size(v2,1), size(v2,2), size(v2,3)))

  allocate(t1_initial(size(t1,1), size(t1,2), size(t1,3)))
  allocate(tlf_initial(size(tlf,1), size(tlf,2), size(tlf,3)))
  allocate(t2_initial(size(t2,1), size(t2,2), size(t2,3)))

  allocate(s1_initial(size(s1,1), size(s1,2), size(s1,3)))
  allocate(slf_initial(size(slf,1), size(slf,2), size(slf,3)))
  allocate(s2_initial(size(s2,1), size(s2,2), size(s2,3)))

  allocate(dmx_initial(size(dmx,1), size(dmx,2), size(dmx,3)))
  allocate(dmy_initial(size(dmy,1), size(dmy,2), size(dmy,3)))

  allocate(dhx_initial(size(dhx,1), size(dhx,2), size(dhx,3)))
  allocate(dhy_initial(size(dhy,1), size(dhy,2), size(dhy,3)))

  allocate(ev_initial(size(ev,1), size(ev,2), size(ev,3)))

  allocate(hv_initial(size(hv,1), size(hv,2), size(hv,3), size(hv,4)))

  allocate(u2_cpu(size(u2)), u2_gpu(size(u2)))
  allocate(v2_cpu(size(v2)), v2_gpu(size(v2)))
  allocate(t2_cpu(size(t2)), t2_gpu(size(t2)))
  allocate(s2_cpu(size(s2)), s2_gpu(size(s2)))

  ! Initialize input parameters with random values
  call random_number(u)
  call random_number(v)
  call random_number(w)
  call random_number(p)

  call random_number(u1)
  call random_number(ulf)
  call random_number(u2)

  call random_number(v1)
  call random_number(vlf)
  call random_number(v2)

  call random_number(t1)
  call random_number(tlf)
  call random_number(t2)

  call random_number(s1)
  call random_number(slf)
  call random_number(s2)

  call random_number(dmx)
  call random_number(dmy)

  call random_number(dhx)
  call random_number(dhy)

  call random_number(ev)
  call random_number(hv)

  ! Store initial values
  u_initial = u
  v_initial = v
  w_initial = w
  p_initial = p

  u1_initial = u1
  ulf_initial = ulf
  u2_initial = u2

  v1_initial = v1
  vlf_initial = vlf
  v2_initial = v2

  t1_initial = t1
  tlf_initial = tlf
  t2_initial = t2

  s1_initial = s1
  slf_initial = slf
  s2_initial = s2

  dmx_initial = dmx
  dmy_initial = dmy

  dhx_initial = dhx
  dhy_initial = dhy

  ev_initial = ev
  hv_initial = hv

  ! CPU execution
  do i = 1, num_steps
    call solve_dynamic_equation
  end do

  ! Store CPU results
  u2_cpu = reshape(u2, [size(u2)])
  v2_cpu = reshape(v2, [size(v2)])
  t2_cpu = reshape(t2, [size(t2)])
  s2_cpu = reshape(s2, [size(s2)])

  ! Restore initial values before GPU execution
  u = u_initial
  v = v_initial
  w = w_initial
  p = p_initial

  u1 = u1_initial
  ulf = ulf_initial
  u2 = u2_initial

  v1 = v1_initial
  vlf = vlf_initial
  v2 = v2_initial

  t1 = t1_initial
  tlf = tlf_initial
  t2 = t2_initial

  s1 = s1_initial
  slf = slf_initial
  s2 = s2_initial

  dmx = dmx_initial
  dmy = dmy_initial

  dhx = dhx_initial
  dhy = dhy_initial

  ev = ev_initial
  hv = hv_initial

  ! GPU execution
  !$acc enter data async(async_id) &
  !$acc& copyin(u, v, w, p, u1, v1, t1, s1, u2, v2, t2, s2, &
  !$acc&        ulf, vlf, tlf, slf, dmx, dmy, ev, dhx, dhy, hv, &
  !$acc&        dt, odz, csv, odyv, ocs, ody, odx, in, iv, iw, iu, &
  !$acc&        nx, ny, nz, xu_bgn, xu_end, yv_bgn, yv_end, nxf, nyf)
  do i = 1, num_steps
    call solve_dynamic_equation_gpu(async_id)
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(u2, v2, t2, s2) &
  !$acc& delete(u, v, w, p, u1, v1, t1, s1, ulf, vlf, tlf, slf, &
  !$acc&        dmx, dmy, ev, dhx, dhy, hv, dt, &
  !$acc&        odz, csv, odyv, ocs, ody, odx, in, iv, iw, iu, &
  !$acc&        nx, ny, nz, xu_bgn, xu_end, yv_bgn, yv_end, nxf, nyf)
  !$acc wait(async_id)

  ! Store GPU results
  u2_gpu = reshape(u2, [size(u2)])
  v2_gpu = reshape(v2, [size(v2)])
  t2_gpu = reshape(t2, [size(t2)])
  s2_gpu = reshape(s2, [size(s2)])

  ! Assertions
  call assert_allclose(u2_gpu, size(u2_gpu), u2_cpu, size(u2_cpu), 1e-10, 1e-10, "Array u2")
  call assert_allclose(v2_gpu, size(v2_gpu), v2_cpu, size(v2_cpu), 1e-10, 1e-10, "Array v2")
  call assert_allclose(t2_gpu, size(t2_gpu), t2_cpu, size(t2_cpu), 1e-10, 1e-10, "Array t2")
  call assert_allclose(s2_gpu, size(s2_gpu), s2_cpu, size(s2_cpu), 1e-10, 1e-10, "Array s2")

  write (*, *) 'Test completed'

end subroutine solve_dynamic_equation_unit
