!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_adjust_seaice_ts
  implicit none

  call timcom_initialize_gpu
  call adjust_seaice_ts_unit
  call timcom_finalize_gpu

end program test_adjust_seaice_ts

subroutine adjust_seaice_ts_unit
  use timcom_drv, only: adjust_seaice_ts
  use timcom_drv_gpu, only: adjust_seaice_ts_gpu
  use hyperlink, only: in, odz, t2, s2, ssh, qice, aqice, qflux, &
                       m_comm_cart, r8type3d, &
                       nbid, ndim, symm_np, &
                       nx, ny, nz, kb
  implicit none

  real(8), dimension(:), allocatable :: t2_cpu, s2_cpu
  real(8), dimension(:), allocatable :: t2_gpu, s2_gpu
  real(8), dimension(:,:,:), allocatable :: t2_initial, s2_initial
  real(8), dimension(:,:,:), allocatable :: in_r8
  real(8), dimension(:,:), allocatable :: kb_r8
  real :: time1, time2, ct, gt
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  ! Initialize variables
  allocate(t2_initial(size(t2,1), size(t2,2), size(t2,3)), s2_initial(size(s2,1), size(s2,2), size(s2,3)))
  allocate(t2_cpu(size(t2)), s2_cpu(size(s2)))
  allocate(t2_gpu(size(t2)), s2_gpu(size(s2)))
  allocate(in_r8(size(in,1), size(in,2), size(in,3)), kb_r8(size(kb,1), size(kb,2)))

  ! Initialize input parameters
  call random_number(in_r8)
  in = nint(in_r8)
  call random_number(odz)
  call random_number(t2)
  call random_number(s2)
  call random_number(ssh)
  call random_number(qice)
  call random_number(aqice)
  call random_number(qflux)
  call random_number(kb_r8)
  kb = nint(kb_r8)

  ! Store initial values
  t2_initial = t2
  s2_initial = s2

  ! CPU execution
  ct = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call adjust_seaice_ts(1,num_steps,ct)
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do
  t2_cpu = reshape(t2, [size(t2)])
  s2_cpu = reshape(s2, [size(s2)])

  ! Restore initial values before GPU execution
  t2 = t2_initial
  s2 = s2_initial

  ! GPU execution
  !$acc enter data async(async_id) &
  !$acc& copyin(t2, s2, in, qice, aqice, qflux, kb, odz, ssh, &
  !$acc&        nx, ny, nz)
  gt = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call adjust_seaice_ts_gpu(1,num_steps,gt,async_id)
    !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc exit data async(async_id) &
  !$acc copyout(t2, s2) &
  !$acc delete(kb, odz, ssh, in, qice, aqice, qflux, nx, ny, nz)
  !$acc wait(async_id)

  t2_gpu = reshape(t2, [size(t2)])
  s2_gpu = reshape(s2, [size(s2)])

  ! Assertions
  call assert_allclose(t2_gpu, size(t2_gpu), t2_cpu, size(t2_cpu), 1e-10, 1e-10, "Array t2")
  call assert_allclose(s2_gpu, size(s2_gpu), s2_cpu, size(s2_cpu), 1e-10, 1e-10, "Array s2")

  ! Print timing results
  write (*, *) 'Timing for CPU & GPU: ', ct, gt
  write (*, *) 'Speedup ratio: ', ct/gt

end subroutine adjust_seaice_ts_unit
