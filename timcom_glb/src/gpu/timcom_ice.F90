!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module timcom_ice_gpu
  use timcom_const
  use timcom_comm
  implicit none

contains

subroutine ice_formation_gpu(nx, ny, nz, kb, odz, ssh, t2, s2, qice, aqice, qflux, if_calu_qflux, dt_cpl, async_id)
  implicit none

  integer, intent(in) ::    &
    nx, ny, nz, async_id

  integer(2), intent(in) :: &
    kb(0:nx+1,0:ny+1)

  real(r8), intent(in) ::   &
    odz(nz),                &
    ssh(nx,ny), dt_cpl

  real(r8), intent(inout) :: &
    t2(-1:nx+2,-1:ny+2,nz),  &
    s2(-1:nx+2,-1:ny+2,nz),  &
    qice(nx,ny),  &!tot column cooling from ice form (in C*cm)
    aqice(nx,ny), &!sum of accumulated ice heat flux since tlast
    qflux(nx,ny)   !internal ocn heat flux due to ice formation

  logical :: if_calu_qflux

  integer, parameter :: &
    kmxice = 1    !lowest level from which to integrate ice formation

  real(r8), parameter :: &
    time_weight = 1.d0,  &
    salice = sea_ice_salinity*ppt_to_salt, & !sea ice salinity in msu
    salref = ocn_ref_salinity*ppt_to_salt, & !ocean ref salinity in msu
    ref_val = salref - salice, & !tracer reference value
    cp_over_lhfusion = &
      rho_sw*cp_sw/(3.34d9*1.0d0) !cp_sw/latent_heat_fusion = rho_sw*cp_sw/(latent_heat_fusion*rho_fw)

  real(r8) :: &
    fw_freeze(nx,ny),  & !water flux at T points due to frazil ice formation
    tfrz(nx,ny),       & !freezing temperature of water in deg C
    potice(nx,ny),     & !potential amount of ice formation
    work1(nx,ny),      & !working variable
    work2(nx,ny)

  real(r8) :: trcr1, trcr2

  integer :: i, j, k

  !$acc enter data create(potice, tfrz, ssh, work1) async(async_id)
  !$acc parallel loop collapse(2) private(i,j,k, trcr1, trcr2) async(async_id)
  do j = 1, ny
    do i = 1, nx
      !$acc loop seq
      do k = 1, 1
        qice(i,j) = 0.d0
        aqice(i,j) = 0.d0
        potice(i,j) = 0.d0

        trcr1 = t2(i,j,k)
        trcr2 = s2(i,j,k)*ppt_to_salt

        if (k >= 2 .and. k <= kmxice) then
          tfrz(i,j) = -1.8d0 ! call tfreez(nx,ny,tfrz,trcr(i,j,k,2))
          if ( k <= int(kb(i,j)) ) then
            potice(i,j) = (tfrz(i,j) - trcr1)/odz(k)
          end if
          potice(i,j) = max(potice(i,j), qice(i,j))
          trcr1 = trcr1 + potice(i,j)*odz(k)
          if(ref_val /= 0.0)  then
            trcr2 = trcr2 + ref_val*potice(i,j)*cp_over_lhfusion*odz(k)
          end if
          qice(i,j) = qice(i,j) - potice(i,j)
        endif

        if (k .eq. 1) then
          tfrz(i,j) = -1.8d0 ! call tfreez(nx,ny,tfrz,trcr(i,j,k,2))
          work1(i,j) = 1.d0/odz(k) + ssh(i,j)
          if ( k <= int(kb(i,j)) ) then
            potice(i,j) = (tfrz(i,j) - trcr1)*work1(i,j)
          end if
          potice(i,j) = max(potice(i,j), qice(i,j))
          trcr1 = trcr1 + potice(i,j)/work1(i,j)
          if (ref_val /= 0.0)  then
            trcr2 = trcr2 + ref_val*potice(i,j)*cp_over_lhfusion/work1(i,j)
          end if
          qice(i,j) = qice(i,j) - potice(i,j)
          aqice(i,j) = aqice(i,j) + time_weight*qice(i,j)
          tfrz(i,j) = -1.8d0 ! call tfreez(nx,ny,tfrz,trcr(i,j,k,2))
          if ( k <= int(kb(i,j)) ) then
            potice(i,j) = (tfrz(i,j) - trcr1)*work1(i,j)
          end if
          potice(i,j) = max(potice(i,j), aqice(i,j))
          trcr1 = trcr1 + potice(i,j)/work1(i,j)
          if (ref_val /= 0.0)  then
            trcr2 = trcr2 + ref_val*potice(i,j)*cp_over_lhfusion/work1(i,j)
          end if
          aqice(i,j) = aqice(i,j) - time_weight*potice(i,j)
        end if

        t2(i,j,k) = trcr1
        s2(i,j,k) = trcr2/ppt_to_salt
      end do
    end do
  end do
  !$acc exit data delete(potice, tfrz, ssh, work1) async(async_id)

end subroutine ice_formation_gpu

end module timcom_ice_gpu

