module timcom_ice
  use timcom_const
  use timcom_comm
  implicit none

contains

subroutine ice_formation(nx, ny, nz, kb, odz, ssh, t2, s2, qice, aqice, qflux, if_calu_qflux, dt_cpl)
  implicit none

  integer, intent(in) ::    &
    nx, ny, nz

  integer(2), intent(in) :: &
    kb(0:nx+1,0:ny+1)

  real(r8), intent(in) ::   &
    odz(nz),                &
    ssh(nx,ny)

  real(r8), intent(inout) :: &
    t2(-1:nx+2,-1:ny+2,nz),  &
    s2(-1:nx+2,-1:ny+2,nz),  &
    qice(nx,ny),  &!tot column cooling from ice form (in C*cm)
    aqice(nx,ny), &!sum of accumulated ice heat flux since tlast
    qflux(nx,ny), &!internal ocn heat flux due to ice formation
    dt_cpl
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
    work2(nx,ny),      &     
    trcr(nx,ny,nz,2)     

  integer :: k
  qice = 0.d0
  aqice = 0.d0
  potice = 0.d0

!===================================================================================
! ice formation
!===================================================================================
!  compute frazil ice formation for sub-surface layers. if ice
!  forms in lower layers but layers above are warm - the heat is
!  used to melt the ice. the ice formation occurs at salinity, Si.
!  this volume is replaced with an equal volume at the salinity of
!  the layer above. the total ice heat flux is accumulated.
!
!  WARNING: unless a monotone advection scheme is in place,
!  advective errors could lead to temps that are far below freezing
!  in some locations and this scheme will form lots of ice.
!  ice formation should be limited to the top layer (kmxice=1)
!  if the advection scheme is not monotone.
!
!  modified from pop2 ice.F90, partial_bottom_cells == .False.
!                              lfw_as_salt_flx == .True.
!                              sfc_layer_type == sfc_layer_varthick
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
! assign working variable and change the salinity unit to (g/g)
!------------------------------------------------------------------------
  trcr(:,:,:,1) = t2(1:nx,1:ny,1:nz)
  trcr(:,:,:,2) = s2(1:nx,1:ny,1:nz)*ppt_to_salt

  do k = kmxice, 2, -1
  !***
  !*** potice is the potential amount of ice formation
  !*** (potice>0) or melting (potice<0) in layer k
  !***

    call tfreez(nx,ny,tfrz,trcr(:,:,k,2))
    where ( k <= int(kb(1:nx,1:ny)) ) &
    potice = (tfrz - trcr(:,:,k,1))/odz(k)

  !***
  !*** if potice < 0, use the heat to melt any ice
  !*** from lower layers
  !*** if potice > 0, keep on freezing (QICE < 0)
  !***

    potice = max(potice,qice(:,:))

  !***
  !*** adjust tracer values based on freeze/melt (tfrz or ori-val)
  !***

    trcr(:,:,k,1) = trcr(:,:,k,1) + potice*odz(k)

   !ref_val = salref - salice
    if(ref_val /= 0.0)  then
      trcr(:,:,k,2) = trcr(:,:,k,2) + ref_val*potice*cp_over_lhfusion*odz(k)
    end if

  !*** accumulate freezing potential
    qice(:,:) = qice(:,:) - potice

  end do

!-----------------------------------------------------------------------
!  now repeat the above algorithm for the surface layer. when fresh
!  water flux formulation is used, the surface layer does not get
!  any salt from other layers. instead, its volume changes.
!-----------------------------------------------------------------------

   k = 1

   call tfreez(nx,ny,tfrz,trcr(:,:,k,2))
   work1 = 1.d0/odz(k) + ssh

   where ( k <= int(kb(1:nx,1:ny)) )
     potice = (tfrz - trcr(:,:,k,1))*work1
   end where

   potice = max(potice, qice)

   trcr(:,:,k,1) = trcr(:,:,k,1) + potice/work1

   if (ref_val /= 0.0)  then
      trcr(:,:,k,2) = trcr(:,:,k,2) + ref_val*potice*cp_over_lhfusion/work1
   end if

   qice = qice - potice

!-----------------------------------------------------------------------
!  let any residual heat in the upper layer melt previously formed ice
!-----------------------------------------------------------------------
   aqice = aqice + time_weight*qice

!-----------------------------------------------------------------------
!  recalculate freezing potential based on adjusted T.
!  only interested in melt potential now (POTICE < 0) - use this
!  melt to offset any accumulated freezing (AQICE < 0) and
!  adjust T and S to reflect this melting. when freshwater flux
!  formulation, compute the associated freshwater flux instead of
!  adjusting S.
!-----------------------------------------------------------------------

   call tfreez(nx,ny,tfrz,trcr(:,:,k,2))

   where (k <= int(kb(1:nx,1:ny)) )
     potice = (tfrz - trcr(:,:,k,1)) * work1
   endwhere

   potice = max(potice, aqice)

   trcr(:,:,k,1) = trcr(:,:,k,1) + potice/work1

  !ref_val = salref - salice
   if (ref_val /= 0.0)  then
      trcr(:,:,k,2) = trcr(:,:,k,2) + ref_val*potice*cp_over_lhfusion/work1
   end if

   aqice = aqice - time_weight*potice

!---------------------------------------------------------------------------
! update T2, S2 and change the salinity unit to (g/kg)
!---------------------------------------------------------------------------
   t2(1:nx,1:ny,1:nz) = trcr(:,:,:,1)
   s2(1:nx,1:ny,1:nz) = trcr(:,:,:,2)/ppt_to_salt

   if(if_calu_qflux) &
     call ice_flx_to_coupler(nx, ny, nz, kb, odz, ssh, dt_cpl, t2, s2, qice, aqice, qflux)

end subroutine ice_formation

subroutine ice_flx_to_coupler(nx, ny, nz, kb, odz, ssh, dt_cpl, t2, s2, qice, aqice, qflux)
  implicit none

  integer,    intent(in) :: nx, ny, nz
  integer(2), intent(in) :: kb(0:nx+1,0:ny+1)
  real(r8),   intent(in) :: odz(nz), ssh(nx,ny), dt_cpl
  real(r8),   intent(inout) :: &
    t2(-1:nx+2,-1:ny+2,nz),  &
    s2(-1:nx+2,-1:ny+2,nz),  &
    qice(nx,ny),  &!tot column cooling from ice form (in C*cm)
    aqice(nx,ny), &!sum of accumulated ice heat flux since tlast
    qflux(nx,ny)   !internal ocn heat flux due to ice formation

  real(r8) :: trcr(nx,ny,nz,2), work1(nx,ny), work2(nx,ny), tfrz(nx,ny)
  integer :: k

  trcr(:,:,:,1) = t2(1:nx,1:ny,1:nz)
  trcr(:,:,:,2) = s2(1:nx,1:ny,1:nz)*ppt_to_salt

  k = 1

  call tfreez(nx,ny,tfrz,trcr(:,:,k,2))
  
  work1 = 1.d0/odz(k) + ssh
  work2 = 0.d0

  where ( k <= int(kb(1:nx,1:ny)) )
    work2 = (tfrz - trcr(:,:,k,1))*work1
  end where

  aqice = aqice*0.5d0

  where( aqice .lt. 0.d0)
    qice = -aqice
  elsewhere
    qice = work2
  end where

  qflux = qice/dt_cpl/hflux_factor
end subroutine

subroutine tfreez(nx, ny, TFRZ, SALT)
  implicit none

! !DESCRIPTION:
!  This function computes the freezing point of salt water.

  integer, intent(in) :: nx, ny

  real(r8), dimension(nx,ny), intent(in) :: &
    SALT    !salinity in model units (g/g)

  real(r8), dimension(nx,ny), intent(out) :: &
    TFRZ    !freezing temperature of water in deg C

  !TFRZ = -0.0544d0*SALT*salt_to_ppt
  TFRZ = -1.8d0

end subroutine tfreez

end module timcom_ice

