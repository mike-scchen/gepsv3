!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_solve_hydrostatic_equation
  implicit none

  call timcom_initialize_gpu
  call solve_hydrostatic_equation_unit
  call timcom_finalize_gpu

end program test_solve_hydrostatic_equation

subroutine solve_hydrostatic_equation_unit
  use timcom_drv, only: solve_hydrostatic_equation
  use timcom_drv_gpu, only: solve_hydrostatic_equation_gpu
  use hyperlink, only: p, p0, t2, s2, rho, &
                       z_grid, dzw, dz, in, grav, &
                       m_comm_cart, r8type3d, &
                       nbid, ndim, symm_np, &
                       nx, ny, nz, nzf
  implicit none

  real(8), dimension(:), allocatable :: p_cpu, p_gpu
  real(8), dimension(:), allocatable :: rho_cpu, rho_gpu
  real(8), dimension(:,:,:), allocatable :: p_initial
  real(8), dimension(:,:,:), allocatable :: rho_initial
  real(8), dimension(:,:), allocatable :: p0_initial
  real(8), dimension(:,:,:), allocatable :: t2_initial, s2_initial
  real(8), dimension(:,:,:), allocatable :: in_r8
  integer, parameter :: async_id = 1
  integer, parameter :: num_steps = 16
  integer :: i

  ! Initialize variables
  allocate(p_initial(size(p,1), size(p,2), size(p,3)))
  allocate(rho_initial(size(rho,1), size(rho,2), size(rho,3)))
  allocate(p0_initial(size(p0,1), size(p0,2)))
  allocate(t2_initial(size(t2,1), size(t2,2), size(t2,3)))
  allocate(s2_initial(size(s2,1), size(s2,2), size(s2,3)))
  allocate(p_cpu(size(p)), p_gpu(size(p)))
  allocate(rho_cpu(size(rho)), rho_gpu(size(rho)))
  allocate(in_r8(size(in,1), size(in,2), size(in,3)))

  ! Initialize input parameters with random values
  call random_number(in_r8)
  in = nint(in_r8)
  call random_number(t2)
  call random_number(s2)
  call random_number(p0)
  call random_number(rho)
  call random_number(p)

  ! Store initial values
  p_initial = p
  rho_initial = rho
  p0_initial = p0
  t2_initial = t2
  s2_initial = s2

  ! CPU execution
  do i = 1, num_steps
    call solve_hydrostatic_equation
  end do
  p_cpu = reshape(p, [size(p)])
  rho_cpu = reshape(rho, [size(rho)])

  ! Restore initial values before GPU execution
  p = p_initial
  rho = rho_initial
  p0 = p0_initial
  t2 = t2_initial
  s2 = s2_initial

  ! GPU execution
  !$acc enter data async(async_id) &
  !$acc& copyin(rho, p) &
  !$acc& copyin(in, nx, ny, nz, t2, s2, p0, z_grid, dzw)
  do i = 1, num_steps
    call solve_hydrostatic_equation_gpu(async_id)
  end do
  !$acc exit data async(async_id) &
  !$acc& copyout(rho, p) &
  !$acc& delete(in, nx, ny, nz, t2, s2, p0, z_grid, dzw)
  !$acc wait(async_id)

  p_gpu = reshape(p, [size(p)])
  rho_gpu = reshape(rho, [size(rho)])

  ! Assertions
  call assert_allclose(p_gpu, size(p_gpu), p_cpu, size(p_cpu), 1e-10, 1e-10, "Array p")
  call assert_allclose(rho_gpu, size(rho_gpu), rho_cpu, size(rho_cpu), 1e-10, 1e-10, "Array rho")

  write (*, *) 'Test complete'

end subroutine solve_hydrostatic_equation_unit
