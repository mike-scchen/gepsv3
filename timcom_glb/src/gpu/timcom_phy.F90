!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module timcom_phy_gpu
  use timcom_const
  use timcom_comm
  implicit none

contains

subroutine pp82_gpu(nx, ny, nz, w, u2, v2, rho, hbk, ev, hv, async_id)
! ----------------------------------------------------------------------
! Calculate gradient Ri based vertical mixing coefficents
! so called "eddy viscosity and diffusivity"
! as per Pacanowski and Philander(1982)

! In the loop below, EV,HV units are cm-cm/s
! but EV,HV normalization by DZ is done
  use hyperlink, only: iw, odzw, dt
  implicit none

  integer, intent(in) :: async_id
  integer :: nx, ny, nz
  real(r8), intent(in) :: &
    w(nx,ny,nz+1), &
    u2(-1:nx+2,-1:ny+2,nz),  &
    v2(-1:nx+2,-1:ny+2,nz),  &
    rho(nx,ny,nz), &
    hbk(nz)

  real(r8), intent(out) :: &
    ev(nx,ny,nz), &
    hv(nx,ny,nz,2)

  integer :: i, j, k, l
  real(r8) :: tmpw, hbk0, tmp, temp, Ri, evisc_local, emax
  real(r8), parameter :: rzmx = 50.d0, ri_max = -0.9d0, evisc = 5.d0

  !$acc parallel loop collapse(3) async(async_id) &
  !$acc& private(tmpw, hbk0, tmp, temp, Ri, evisc_local, emax, l)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ev(i,j,k) = 0.d0
        hv(i,j,k,1) = 0.d0
        hv(i,j,k,2) = 0.d0

        if (k < nz) then
          l = k + 1
          tmpw = 1.d0/(rzmx*odzw(l))
          hbk0 = hbk(k)
          temp = tmpw*dabs(w(i,j,l)) !units of cm-cm/s
          Ri = grav*(rho(i,j,l)-rho(i,j,k))*odzw(l)/(odzw(l)**2 &
                   *(0.001d0+(u2(i,j,l)-u2(i,j,k))**2+(v2(i,j,l)-v2(i,j,k))**2))
          Ri = max(Ri_max, Ri)
          tmp = 1.d0/(1.d0 + Ri)
          temp = tmp*temp

          evisc_local = min(evisc*tmp**2, 100.d0)
          hv(i,j,k,1) = (evisc_local*tmp + temp + hbk0)
          ev(i,j,k)   = (evisc_local     + temp)

          emax = 0.2d0/(dt*odzw(l)**2)
          tmp = dble(iw(i,j,l))*odzw(l)
          hv(i,j,k,1) = tmp*min(emax, hv(i,j,k,1))
          ev(i,j,k)   = tmp*min(emax, ev(i,j,k))
        endif

        hv(i,j,k,2) = hv(i,j,k,1)*0.1d0
      end do
    end do
  end do

end subroutine pp82_gpu

end module timcom_phy_gpu
