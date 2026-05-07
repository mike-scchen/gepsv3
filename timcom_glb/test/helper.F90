!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assert_allclose(actual, n_actual, desired, n_desired, rtol, atol, err_msg)
  implicit none

  integer, intent(in) :: n_actual, n_desired
  real(8), dimension(n_actual), intent(in) :: actual
  real(8), dimension(n_desired), intent(in) :: desired
  real(8), intent(in) :: rtol, atol
  character(len=*), intent(in), optional :: err_msg
  real(8), parameter :: eps = nearest(real(0.0, 8), 1.0)
  real(8) :: rel_diff, abs_diff
  logical :: equal
  integer :: i, err_count

  if (n_actual .ne. n_desired) then
    print *, "Arrays have different lengths!"
    if (present(err_msg)) print *, err_msg
    call exit(1)
  end if

  equal = .true.
  err_count = 0
  rel_diff = 0.0
  abs_diff = 0.0
  do i = 1, n_actual
    rel_diff = max(rel_diff, abs((actual(i) - desired(i))/(desired(i) + eps)))
    abs_diff = max(abs_diff, abs(actual(i) - desired(i)))
    if (abs(actual(i) - desired(i)) > atol + rtol*abs(desired(i))) then
      equal = .false.
      if (err_count .lt. 5) then
        print '(5g0)', "Desired = ", desired(i), " but get actual = ", actual(i)
        print '(5g0)', "Absolute difference = ", abs(actual(i) - desired(i)), "; Tolerance = ", atol + rtol*abs(desired(i))
      end if
      err_count = err_count + 1
    end if
  end do

  if (.not. equal) then
    write (*, '(2(1A, 1pe12.4))') &
      "Arrays are not close within tolerance rtol =", rtol, &
      ", atol =", atol
    write (*, '(1A, 1pe15.7)') "Max relative difference = ", rel_diff
    write (*, '(1A, 1pe15.7)') "Max absolute difference = ", abs_diff
    if (present(err_msg)) then
      print '(5g0)', "In array ", err_msg, ","
      print '(5g0)', err_count, "/", n_actual, " elements are not close within tolerance rtol = ", rtol, " atol =", atol
    end if
    call exit(1)
  end if

end subroutine assert_allclose

subroutine assert_allclose_r8type3d(gpu_array, cpu_array, nx, ny, nz, rtol, atol, err_msg)
  implicit none

  integer, intent(in) :: nx, ny, nz
  real(8), dimension(-1:nx+2, -1:ny+2, nz), intent(in) :: gpu_array, cpu_array
  real(8), intent(in) :: rtol, atol
  character(len=*), intent(in), optional :: err_msg

  ! Check main region
  call assert_allclose(gpu_array(1:nx, 1:ny, 1:nz), nx*ny*nz, &
                      cpu_array(1:nx, 1:ny, 1:nz), nx*ny*nz, &
                      rtol, atol, trim(err_msg)//" main region")

  ! Check north boundary (excluding corners)
  call assert_allclose(gpu_array(1:nx, -1:0, 1:nz), 2*nx*nz, &
                      cpu_array(1:nx, -1:0, 1:nz), 2*nx*nz, &
                      rtol, atol, trim(err_msg)//" north boundary")

  ! Check south boundary (excluding corners)
  call assert_allclose(gpu_array(1:nx, ny+1:ny+2, 1:nz), 2*nx*nz, &
                      cpu_array(1:nx, ny+1:ny+2, 1:nz), 2*nx*nz, &
                      rtol, atol, trim(err_msg)//" south boundary")

  ! Check west boundary (excluding corners)
  call assert_allclose(gpu_array(-1:0, 1:ny, 1:nz), 2*ny*nz, &
                      cpu_array(-1:0, 1:ny, 1:nz), 2*ny*nz, &
                      rtol, atol, trim(err_msg)//" west boundary")

  ! Check east boundary (excluding corners)
  call assert_allclose(gpu_array(nx+1:nx+2, 1:ny, 1:nz), 2*ny*nz, &
                      cpu_array(nx+1:nx+2, 1:ny, 1:nz), 2*ny*nz, &
                      rtol, atol, trim(err_msg)//" east boundary")

end subroutine assert_allclose_r8type3d

