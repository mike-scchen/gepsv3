!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_check_incompressibility
  implicit none

  call timcom_initialize_gpu
  call check_incompressibility_unit
  call timcom_finalize_gpu

end program test_check_incompressibility

subroutine check_incompressibility_unit
  use timcom_drv, only: check_incompressibility, incomp_res
  use timcom_drv_gpu, only: check_incompressibility_gpu, incomp_res_gpu
  use hyperlink, only: u, v, w, odx, ody, odz, ocs, csv, in, kb
  use hyperlink, only: nx, ny, nz, nxf, nyf, nzf, m_comm_cart, myid
  implicit none

  ! Variables for storing initial values
  real(8), dimension(:,:,:), allocatable :: u_initial, v_initial, w_initial
  real(8), dimension(:), allocatable :: odx_initial, ody_initial, odz_initial, ocs_initial, csv_initial
  integer(2), dimension(:,:), allocatable :: kb_initial
  integer(2), dimension(:,:,:), allocatable :: in_initial

  real(8) :: incomp_res_cpu

  ! Timing variables
  real(8) :: time1, time2, ct, gt

  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  ! Allocate memory for arrays
  allocate(u_initial(0:nxf+1, 0:ny+1, nz))
  allocate(v_initial(0:nx+1, 0:nyf+1, nz))
  allocate(w_initial(0:nx+1, 0:ny+1, 0:nzf))
  allocate(odx_initial(ny))
  allocate(ody_initial(ny))
  allocate(odz_initial(nz))
  allocate(ocs_initial(ny))
  allocate(csv_initial(nyf))
  allocate(kb_initial(0:nx+1, 0:ny+1))
  allocate(in_initial(0:nx+1, 0:ny+1, 0:nz))

  ! Initialize arrays with random values
  call random_number(u)
  call random_number(v)
  call random_number(w)
  call random_number(odx)
  call random_number(ody)
  call random_number(odz)
  call random_number(ocs)
  call random_number(csv)

  ! Store initial values
  u_initial = u
  v_initial = v
  w_initial = w
  odx_initial = odx
  ody_initial = ody
  odz_initial = odz
  ocs_initial = ocs
  csv_initial = csv
  kb_initial = kb
  in_initial = in

  ! CPU execution
  ct = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call check_incompressibility
    call cpu_time(time2)
    if (i .ne. 1) ct = ct + time2 - time1
  end do

  ! Store CPU results
  incomp_res_cpu = incomp_res

  ! Restore initial values for GPU execution
  u = u_initial
  v = v_initial
  w = w_initial
  odx = odx_initial
  ody = ody_initial
  odz = odz_initial
  ocs = ocs_initial
  csv = csv_initial
  kb = kb_initial
  in = in_initial

  ! GPU execution
  !$acc enter data async(async_id) &
  !$acc& create(incomp_res_gpu) &
  !$acc& copyin(u, v, w, odx, csv, ocs, ody, odz, in, nx, ny, nz, kb)
  gt = 0.0
  do i = 1, num_steps
    call cpu_time(time1)
    call check_incompressibility_gpu(async_id)
    !$acc wait(async_id)
    call cpu_time(time2)
    if (i .ne. 1) gt = gt + time2 - time1
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(incomp_res_gpu) &
  !$acc& delete(u, v, w, odx, csv, ocs, ody, odz, in, nx, ny, nz, kb)
  !$acc wait(async_id)

  ! Compare CPU and GPU results for incomp_res
  call assert_allclose([incomp_res_gpu], 1, [incomp_res_cpu], 1, &
                      1e-10, 1e-10, "incomp_res")

  ! Print performance information
  if (myid == 0) then
    write (*, *) 'Timing for CPU & GPU: ', ct, gt
    write (*, *) 'Speedup ratio: ', ct/gt
    write (*, *) 'CPU incomp_res: ', incomp_res_cpu
    write (*, *) 'GPU incomp_res: ', incomp_res_gpu
  end if

  ! Deallocate memory
  deallocate(u_initial, v_initial, w_initial)
  deallocate(odx_initial, ody_initial, odz_initial, ocs_initial, csv_initial)
  deallocate(kb_initial, in_initial)

end subroutine check_incompressibility_unit
