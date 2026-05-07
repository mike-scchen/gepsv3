module timcom_cpl
  use timcom_const
  use timcom_comm, only: myid_timcom
  use timcom_general
  implicit none

  type :: cpl_data_2d
    integer :: nd(3)
    real(kind=r8), pointer :: lon(:), lat(:)
    real(kind=r8), pointer :: rAttr(:,:,:)
  end type

  type :: smat_panel
    integer, pointer :: ij(:,:,:)
    real(r8), pointer :: wgt(:,:,:)
  end type 

  type(cpl_data_2d), pointer :: a2x_a(:), x2o_o(:), o2x_x(:), a2o_x(:)
  type(smat_panel),  pointer :: a2o_smat(:) !ocn_grid_num
  
contains

subroutine init_timcom_cpl
  use grids, only: ocn_grid, n_x2o_flux, n_o2x_flux
  use datm,  only: atm_grid, n_a2x_flux
  implicit none

  integer :: i 

  allocate(x2o_o(ocn_grid_num))
  allocate(o2x_x(ocn_grid_num))
  allocate(a2o_x(ocn_grid_num))

  do i = 1, ocn_grid_num
    call init_cpl_data(x2o_o(i),  &
                     n_x2o_flux,  &
                 ocn_grid(i)%nx,  &
                 ocn_grid(i)%ny,  &
                 ocn_grid(i)%x_grid, &
                 ocn_grid(i)%y_grid)

    call init_cpl_data(o2x_x(i),  &
                     n_o2x_flux,  &
                 ocn_grid(i)%nx,  &
                 ocn_grid(i)%ny,  &
                 ocn_grid(i)%x_grid, &
                 ocn_grid(i)%y_grid)

    call init_cpl_data(a2o_x(i),  &
                     n_a2x_flux,  &
                 ocn_grid(i)%nx,  &
                 ocn_grid(i)%ny,  &
                 ocn_grid(i)%x_grid, &
                 ocn_grid(i)%y_grid)
  end do

  allocate(a2x_a(atm_grid_num))
  do i = 1, atm_grid_num
    call init_cpl_data(a2x_a(i), &
                     n_a2x_flux, &
               atm_grid(i)%nlon, & 
               atm_grid(i)%nlat, &
               atm_grid(i)%lon,  &
               atm_grid(i)%lat)
  end do

  allocate(a2o_smat(ocn_grid_num))
  do i = 1, ocn_grid_num
    call init_cpl_smat(a2x_a(1), a2o_x(i), a2o_smat(i))
  end do

  call ocn_export_cpl
 
end subroutine init_timcom_cpl

subroutine init_cpl_data(cpl_data, nflux, nlon, nlat, lon, lat)
  implicit none

  type(cpl_data_2d) :: cpl_data
  integer, intent(in) :: nflux, nlon, nlat
  real(r8) :: lon(nlon), lat(nlat)

  allocate(cpl_data%lon(nlon))
  allocate(cpl_data%lat(nlat))
  allocate(cpl_data%rAttr(nflux,nlon,nlat)) 

  cpl_data%lon = lon
  cpl_data%lat = lat
  cpl_data%rAttr = 0.d0
  cpl_data%nd(1) = nflux
  cpl_data%nd(2) = nlon
  cpl_data%nd(3) = nlat

end subroutine init_cpl_data

subroutine init_cpl_smat(src_data, dst_data, smat)
  implicit none

  type(cpl_data_2d) :: src_data, dst_data
  type(smat_panel) :: smat

  integer :: i, j, ip, jp, ii, jj, n
  real(r8) :: odydat, wgt_lat, odxdat, wgt_lon
  
  allocate(smat%ij(2,dst_data%nd(2),dst_data%nd(3)))
  allocate(smat%wgt(2,dst_data%nd(2),dst_data%nd(3)))

  do jj = 1, dst_data%nd(3)
    j = 0
    do n = 1, src_data%nd(3)
      if(src_data%lat(n)-dst_data%lat(jj) .ge. 0.d0) then
        j = n-1
        exit
      end if
    end do
    jp = j + 1
    odydat = 1.d0/(src_data%lat(j+1)-src_data%lat(j))
    wgt_lat = (dst_data%lat(jj)-src_data%lat(j))*odydat
    
    do ii = 1, dst_data%nd(2)
      i = 0
      do n = 1, src_data%nd(2)
        if(src_data%lon(n)-dst_data%lon(ii) .ge. 0.d0) then
          i = n-1
          exit
        end if
      end do
      if(i .eq. 0 .or. i .eq. src_data%nd(2)) then
        i  = src_data%nd(2)
        ip = 1
      else
        ip = i + 1
      end if
      odxdat = 1.d0/dmod(src_data%lon(ip)+360.d0-src_data%lon(i),360.d0)
      wgt_lon = dmod(dst_data%lon(ii)-src_data%lon(i)+360.d0,360.d0)*odxdat
 
      smat%ij(1,ii,jj)  = i
      smat%ij(2,ii,jj)  = j
      smat%wgt(1,ii,jj) = wgt_lon
      smat%wgt(2,ii,jj) = wgt_lat
    end do
  end do
end subroutine init_cpl_smat

subroutine smat_remapping(src_data, dst_data, smat)
  implicit none

  type(cpl_data_2d) :: src_data, dst_data
  type(smat_panel) :: smat

  integer :: i, j, ip, jp, ii, jj
  real(r8), dimension(src_data%nd(1)) :: &
    dat_sw, dat_se, dat_nw, dat_ne, tmp1, tmp2
  real(r8) :: wgt_lon, wgt_lat

  do jj = 1, dst_data%nd(3)
    do ii = 1, dst_data%nd(2)
      i = smat%ij(1,ii,jj)
      j = smat%ij(2,ii,jj)
      ip = i + 1
      jp = j + 1
      wgt_lon = smat%wgt(1,ii,jj)
      wgt_lat = smat%wgt(2,ii,jj)
     
      dat_sw = src_data%rAttr(:, i , j)
      dat_se = src_data%rAttr(:, ip, j)
      dat_nw = src_data%rAttr(:, i,  jp)
      dat_ne = src_data%rAttr(:, ip, jp)

      tmp1 = (dat_se-dat_sw)*wgt_lon + dat_sw
      tmp2 = (dat_ne-dat_nw)*wgt_lon + dat_nw

      dst_data%rAttr(:,ii,jj) = (tmp2-tmp1)*wgt_lat + tmp1 
    end do
  end do
end subroutine smat_remapping

subroutine datm_export_cpl
  use datm,  only: atm_grid, n_a2x_flux
  implicit none

  integer :: i, n

  do i = 1, atm_grid_num
    do n = 1, n_a2x_flux
      if(associated(atm_grid(i)%flx(n)%srcm)) then
        a2x_a(i)%rAttr(n,:,:) = atm_grid(i)%flx(n)%srcm
!        if(n.eq.1) write(*,*) "check a2x expt:", myid, n, a2x_a(1)%rAttr(1,60,55)
      else
        a2x_a(i)%rAttr(n,:,:) = 0.d0
      end if
    end do 
  end do
!  write(*,*) "check a2x expt:", myid, a2x_a(1)%rAttr(1,60,55)
end subroutine datm_export_cpl

subroutine a2o_remap_cpl
  implicit none

  integer i

  do i = 1, ocn_grid_num
!    write(*,*) "check a2x remp:", myid, a2x_a(i)%rAttr(1,60,55)
!    write(*,*) "check a2o smat:", myid, a2o_smat(i)%ij(1,100,100), &
!                                        a2o_smat(i)%ij(2,100,100), &
!                                        a2o_smat(i)%wgt(1,100,100),&
!                                        a2o_smat(i)%wgt(2,100,100)
    call smat_remapping(a2x_a(1), a2o_x(i), a2o_smat(i))
!    write(*,*) "check a2o remp:", myid, a2o_x(i)%rAttr(1,100,100), &
!                                        a2o_x(i)%rAttr(2,100,100)
  end do
end subroutine a2o_remap_cpl

subroutine ocn_export_cpl
  use grids
  implicit none

  integer :: i
  
  do i = 1, ocn_grid_num
    o2x_x(i)%rAttr(idx_o2x_ssu,:,:) = ocn_grid(i)%u2(:,:,1)/100.d0
    o2x_x(i)%rAttr(idx_o2x_ssv,:,:) = ocn_grid(i)%v2(:,:,1)/100.d0
    o2x_x(i)%rAttr(idx_o2x_sst,:,:) = ocn_grid(i)%t2(:,:,1)+273.15d0
    o2x_x(i)%rAttr(idx_o2x_sss,:,:) = ocn_grid(i)%s2(:,:,1)
    o2x_x(i)%rAttr(idx_o2x_msk,:,:) = dble(ocn_grid(i)%in(:,:,1))
  end do
end subroutine ocn_export_cpl

subroutine calu_atmocn_flux
  use grids
  use datm
  use surface_flux, only: cesm_atmocn_fluxes
  implicit none

  integer :: i

  do i = 1, ocn_grid_num
    call cesm_atmocn_fluxes(ocn_grid(i)%i2*ocn_grid(i)%j2, &
                             a2o_x(1)%rAttr(idx_a2x_u10m,:,:), &
                             a2o_x(1)%rAttr(idx_a2x_v10m,:,:), &
                             a2o_x(1)%rAttr(idx_a2x_t02m,:,:), &
                             a2o_x(1)%rAttr(idx_a2x_q02m,:,:), &
                             a2o_x(1)%rAttr(idx_a2x_pslv,:,:), &
                             o2x_x(i)%rAttr(idx_o2x_ssu, :,:), &
                             o2x_x(i)%rAttr(idx_o2x_ssv, :,:), &         
                             o2x_x(i)%rAttr(idx_o2x_sst, :,:), &    
                             o2x_x(i)%rAttr(idx_o2x_msk, :,:), &           
                             x2o_o(i)%rAttr(idx_x2o_taux,:,:), &
                             x2o_o(i)%rAttr(idx_x2o_tauy,:,:), &
                             x2o_o(i)%rAttr(idx_x2o_evap,:,:), &
                             x2o_o(i)%rAttr(idx_x2o_lath,:,:), &
                             x2o_o(i)%rAttr(idx_x2o_senh,:,:), &
                             x2o_o(i)%rAttr(idx_x2o_lwup,:,:))

   x2o_o(i)%rAttr(idx_x2o_lwdn,:,:) = a2o_x(1)%rAttr(idx_a2x_lwdn,:,:)
   x2o_o(i)%rAttr(idx_x2o_swup,:,:) = a2o_x(1)%rAttr(idx_a2x_swup,:,:)
   x2o_o(i)%rAttr(idx_x2o_swdn,:,:) = a2o_x(1)%rAttr(idx_a2x_swdn,:,:)
   x2o_o(i)%rAttr(idx_x2o_rain,:,:) = a2o_x(1)%rAttr(idx_a2x_rain,:,:)
   x2o_o(i)%rAttr(idx_x2o_snow,:,:) = a2o_x(1)%rAttr(idx_a2x_snow,:,:)

  end do

end subroutine calu_atmocn_flux

subroutine ocn_import_cpl(grd, n)
  use grids
  implicit none

  type(ocn_panel) :: grd
  integer, intent(in) :: n

  grd%taux = x2o_o(n)%rAttr(idx_x2o_taux,:,:)*10.d0
  grd%tauy = x2o_o(n)%rAttr(idx_x2o_tauy,:,:)*10.d0
  grd%swdn = x2o_o(n)%rAttr(idx_x2o_swdn,:,:)
  grd%swup = x2o_o(n)%rAttr(idx_x2o_swup,:,:)
  grd%lwdn = x2o_o(n)%rAttr(idx_x2o_lwdn,:,:)
  grd%lwup = x2o_o(n)%rAttr(idx_x2o_lwup,:,:)
  grd%lath = x2o_o(n)%rAttr(idx_x2o_lath,:,:)
  grd%senh = x2o_o(n)%rAttr(idx_x2o_senh,:,:)
  grd%rain = x2o_o(n)%rAttr(idx_x2o_rain,:,:)
  grd%snow = x2o_o(n)%rAttr(idx_x2o_snow,:,:)


end subroutine

end module timcom_cpl
