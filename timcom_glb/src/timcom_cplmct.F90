module timcom_cplmct
  use timcom_const
  use timcom_comm
  implicit none
contains

subroutine ocn_cpl_init(compid, grd)
  use hyperlink, only: nx_grid, nx, ny, &
                       myid_x, myid_y, myid, m_comm_cart
  use timcom_grid, only: ocn_panel
  use cpl_rank,  only: id_gfs, id_rocn
  use cpl_gsmap, only: cpl_gsmap_init
  use cpl_attr,  only: cpl_attr_init
  use cpl_smat,  only: cpl_smat_init
  implicit none

  integer, intent(in) :: compid
  type(ocn_panel) :: grd
  integer :: j, rank_root, seg_strt(ny), seg_leng(ny), ierr

  rank_root = 0
  do j = 1, ny
    seg_strt(j) = (myid_x*nx + (j-1)*nx_grid + 1) + myid_y*ny*nx_grid
    seg_leng(j) = nx
  end do
  !write(*,*) "check docn gsmap size:", myid, seg_leng(1:5), seg_strt(1:5)
  call cpl_gsmap_init(myid, ny, seg_strt, seg_leng, rank_root, m_comm_cart, compid, "docn GSMap:")
  call cpl_smat_init(myid, rank_root, compid, id_gfs, m_comm_cart, "./rmp_timcom2tco.nc")
  call cpl_attr_init(compid, id_gfs, m_comm_cart)

end subroutine ocn_cpl_init

subroutine coupler_data_exch(compid, timer)
  use timcom_general, only: ocn_grid_num, timer_panel
  use timcom_grid,    only: ocn_grid
  use hyperlink,      only: myid
  implicit none

  integer, intent(in) :: compid
  type(timer_panel) :: timer

  real(r8), pointer, dimension(:,:) :: &
    sst, ssu, ssv, &
    u10m, v10m, t02m, q02m, pslv, &
    swup, swdn, lwdn, rain, snow, tgfs, &
    ifrc, vice, vsno

  integer(2), pointer, dimension(:,:) :: &
    msk

  integer :: n, nx, ny

  if(myid .eq. 0) write(*,*) "[TIMCOM] Glb to coupler: ", timer%mjd1, timer%mjd2
 
  do n = 1, ocn_grid_num
    nx = ocn_grid(n)%nx
    ny = ocn_grid(n)%ny
    sst  => ocn_grid(n)%t2(1:nx,1:ny,1)
    ssu  => ocn_grid(n)%u2(1:nx,1:ny,1)
    ssv  => ocn_grid(n)%v2(1:nx,1:ny,1)
    msk  => ocn_grid(n)%in(1:nx,1:ny,1)
    u10m => ocn_grid(n)%u_10
    v10m => ocn_grid(n)%v_10
    t02m => ocn_grid(n)%t_10
    q02m => ocn_grid(n)%q_10
    pslv => ocn_grid(n)%pslv
    swup => ocn_grid(n)%swup
    swdn => ocn_grid(n)%swdn
    lwdn => ocn_grid(n)%lwdn
    rain => ocn_grid(n)%rain
    snow => ocn_grid(n)%snow
    tgfs => ocn_grid(n)%tgfs
    ifrc => ocn_grid(n)%ifrc
    vice => ocn_grid(n)%vice
    vsno => ocn_grid(n)%vsno

    call ocn_cpl_send2gfs(compid, sst, ssu, ssv, ifrc, vice, vsno, msk)
    call ocn_cpl_recv4gfs(compid, u10m, v10m, t02m, q02m, pslv, &
                            swup, swdn, lwdn, rain, snow, tgfs) 
  end do 

end subroutine coupler_data_exch

subroutine ocn_cpl_send2gfs(compid, SST, SSU, SSV, ifrc, vice, vsno, msk)
  use hyperlink, only: nx, ny, myid
  use cpl_rank, only: id_gfs
  use cpl_attr, only: AttrVect_importRAttr, export_AV=>gocn_export_Attr 
  use cpl_sendrecv, only: cpl_send
  implicit none

  integer, intent(in) :: compid
  real(r8), dimension(nx*ny), target, optional, intent(in)  :: SST, SSU, SSV, ifrc, vice, vsno
  integer(2), dimension(nx*ny), intent(in) :: msk
  real(r8), dimension(:), pointer :: send_ptr
  integer :: temp

  if(.not. associated(export_AV%rAttr)) write(*,*) "[TIMCOM ERR]: gocn_export_Atrr Error"
  
  send_ptr=>SST
  call AttrVect_importRAttr(export_AV, "sst", send_ptr)
  
  send_ptr=>SSU
  call AttrVect_importRAttr(export_AV, "ssu", send_ptr)

  send_ptr=>SSV
  call AttrVect_importRAttr(export_AV, "ssv", send_ptr)

  send_ptr=>ifrc
  call AttrVect_importRAttr(export_AV, "ifrc", send_ptr)

  send_ptr=>vice
  call AttrVect_importRAttr(export_AV, "vice", send_ptr)

  send_ptr=>vsno
  call AttrVect_importRAttr(export_AV, "vsno", send_ptr)

  export_AV%rAttr(1,:) = export_AV%rAttr(1,:) + 273.15d0
  export_AV%rAttr(2,:) = export_AV%rAttr(2,:)*0.01d0
  export_AV%rAttr(3,:) = export_AV%rAttr(3,:)*0.01d0
  export_AV%rAttr(4,:) = export_AV%rAttr(4,:)
  export_AV%rAttr(5,:) = export_AV%rAttr(5,:)*0.001d0
  export_AV%rAttr(6,:) = export_AV%rAttr(6,:)*0.001d0

  call cpl_send(compid, id_gfs, export_AV, myid)
end subroutine ocn_cpl_send2gfs

subroutine ocn_cpl_recv4gfs(compid,u10m,v10m,t02m,q02m,pslv,swup,swdn,lwdn,rain,snow,tgfs)
  use hyperlink, only: nx, ny, myid
  use cpl_rank,  only: id_gfs
  use cpl_attr,  only: gocn_recv_gfs_Attr
  use cpl_sendrecv, only: cpl_recv
  implicit none

  integer, intent(in) :: compid
  real(r8), dimension(nx*ny), intent(out) :: u10m, v10m, t02m, q02m, &
                                                 pslv, swup, swdn, lwdn, &
                                                 rain, snow, tgfs
  integer :: n
  if(.not. associated(gocn_recv_gfs_Attr%rAttr)) write(*,*) "[TIMCOM ERR]: gocn_recv_gfs_Atrr Error"
  call cpl_recv(compid, id_gfs, gocn_recv_gfs_Attr, myid)
  u10m = gocn_recv_gfs_Attr%rAttr(1,:)
  v10m = gocn_recv_gfs_Attr%rAttr(2,:)
  t02m = gocn_recv_gfs_Attr%rAttr(3,:)
  q02m = gocn_recv_gfs_Attr%rAttr(4,:)
  pslv = gocn_recv_gfs_Attr%rAttr(5,:)
  swup = gocn_recv_gfs_Attr%rAttr(6,:)
  swdn = gocn_recv_gfs_Attr%rAttr(7,:)
  lwdn = gocn_recv_gfs_Attr%rAttr(8,:)
  rain = gocn_recv_gfs_Attr%rAttr(9,:)
  snow = gocn_recv_gfs_Attr%rAttr(10,:)
  tgfs = gocn_recv_gfs_Attr%rAttr(11,:)

  do n = 1, nx*ny
    if(t02m(n) .lt. 273.15d0) then
        snow(n) = rain(n)
        rain(n) = 0.d0
     end if
  end do
end subroutine ocn_cpl_recv4gfs

subroutine coupler_glb2rocn_exch(compid, timer)
  use timcom_general, only: ocn_grid_num, timer_panel
  use timcom_grid,    only: ocn_grid
  use hyperlink,      only: myid
  implicit none

  integer, intent(in) :: compid
  type(timer_panel) :: timer

  real(r8), pointer, dimension(:,:) :: &
    src, ocn_u, ocn_v, ocn_t, ocn_s

  integer(2), pointer, dimension(:,:) :: &
    var2d, msk

  integer :: n, k, nx, ny, nz

  if(myid .eq. 0) write(*,*) "[TIMCOM] Glb to regional coupler: ", timer%mjd1, timer%mjd2

  do n = 1, ocn_grid_num
    nx = ocn_grid(n)%nx
    ny = ocn_grid(n)%ny
    nz = ocn_grid(n)%nz
    allocate(ocn_u(nx,ny))
    allocate(ocn_v(nx,ny))
    allocate(ocn_t(nx,ny))
    allocate(ocn_s(nx,ny))

    do k = 1, nz
      msk => ocn_grid(n)%in(:,:,k)

      src => ocn_grid(n)%u2(:,:,k)
      call extrapolation(nx, ny, src, msk, ocn_u)

      src => ocn_grid(n)%v2(:,:,k)
      call extrapolation(nx, ny, src, msk, ocn_v)

      src => ocn_grid(n)%t2(:,:,k)
      call extrapolation(nx, ny, src, msk, ocn_t)

      src => ocn_grid(n)%s2(:,:,k)
      call extrapolation(nx, ny, src, msk, ocn_s)

      call ocn_cpl_send2rocn(compid, ocn_u, ocn_v, ocn_t, ocn_s)
    end do
    deallocate(ocn_u)
    deallocate(ocn_v)
    deallocate(ocn_t)
    deallocate(ocn_s)
  end do
end subroutine coupler_glb2rocn_exch

subroutine ocn_cpl_send2rocn(compid, u2, v2, t2, s2)
  use hyperlink, only: nx, ny, myid
  use cpl_rank, only: id_rocn
  use cpl_attr, only: AttrVect_importRAttr, nested_AV=>gocn_nested_Attr
  use cpl_sendrecv, only: cpl_send
  implicit none

  integer, intent(in) :: compid
  real(r8), dimension(nx*ny), target, optional, intent(in)  :: u2, v2, t2, s2
  real(r8), dimension(:), pointer :: send_ptr

  if(.not. associated(nested_AV%rAttr)) write(*,*) "[TIMCOM ERR]: gocn_nested_Atrr Error"
  send_ptr=>u2
  call AttrVect_importRAttr(nested_AV, "ocn_u", send_ptr)

  send_ptr=>v2
  call AttrVect_importRAttr(nested_AV, "ocn_v", send_ptr)

  send_ptr=>t2
  call AttrVect_importRAttr(nested_AV, "ocn_t", send_ptr)

  send_ptr=>s2
  call AttrVect_importRAttr(nested_AV, "ocn_s", send_ptr)

  call cpl_send(compid, id_rocn, nested_AV, myid)
end subroutine ocn_cpl_send2rocn

subroutine extrapolation(nx, ny, src_in, msk_in, dst)
  use hyperlink, only: m_comm_cart, r8type2d, nbid, ndim, symm_np
  use timcom_comm
  implicit none

  integer, intent(in) :: nx, ny
  real(r8), intent(in) :: src_in(-1:nx+2,-1:ny+2)
  integer(2), intent(in) :: msk_in(-1:nx+2,-1:ny+2)
  real(r8), intent(out) :: dst(nx,ny)

  real(r8) :: var(0:nx+1,0:ny+1), otmp_sum
  integer(2) :: msk(0:nx+1,0:ny+1), tmp_sum
  integer :: n, i, j

  msk = msk_in(0:nx+1,0:ny+1)
  var = src_in(0:nx+1,0:ny+1)

  do n = 1, nx/2
    do j = 1, ny
      do i = 1, nx
        tmp_sum = msk(i+1,j) + msk(i-1,j)
        otmp_sum = 1.d0/dble(tmp_sum)
        if(msk(i,j) .eq. 0 .and. tmp_sum .gt. 0) then
           var(i,j) = otmp_sum*(var(i+1,j)*msk(i+1,j) + var(i-1,j)*msk(i-1,j))
           msk(i,j) = 1
        end if
      end do
    end do
    call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, var, symm_np)
  end do

  dst = var(1:nx,1:ny)

end subroutine extrapolation

#ifdef cpl_cice
subroutine coupler_cice_import(timer)
  use timcom_general, only: ocn_grid_num, timer_panel
  use timcom_grid,    only: ocn_grid
  use hyperlink

  use ice_blocks, only: nx_block, ny_block
  use ice_flux, only: zlvl, zlvs, uatm, vatm, wind, potT, Tair, Qa, &
                      rhoa, swvdr, swvdf, swidr, swidf, fsw, flw, &
                      frain, fsnow, uocn, vocn, ss_tltx, ss_tlty, &
                      sst, sss, frzmlt
  use ice_state, only: aice
  implicit none

  type(timer_panel) :: timer

  integer :: i, j, ii, jj, k, iblk

  real(r8), parameter :: &
    frcvdr = 0.28d0, & ! frac of incoming sw in vis direct band
    frcvdf = 0.24d0, & ! frac of incoming sw in vis diffuse band
    frcidr = 0.31d0, & ! frac of incoming sw in near IR direct band
    frcidf = 0.17d0    ! frac of incoming sw in near IR diffuse band

  real(r8), tmpdata(0:nx+1,0:ny+1), pslv_10(0:nx+1,0:ny+1)

  if(myid .eq. 0) write(*,*) "[TIMCOM] Glb to cice coupler: ", timer%mjd1, timer%mjd2, nx_block, ny_block 

  iblk = 1

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%swdn(1:nx,1:ny) + ocn_grid(1)%swup(1:nx,1:ny)
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  fsw(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)
  swvdr = fsw*frcvdr
  swidr = fsw*frcidr
  swvdf = fsw*frcvdf
  swidf = fsw*frcidf

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%lwdn(1:nx,1:ny)
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  flw(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%u_10(1:nx,1:ny)
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  uatm(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%v_10(1:nx,1:ny)
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  vatm(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%t_10(1:nx,1:ny)
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  do j = 1, ny_block
    do i = 1, nx_block
      if (aice(i,j,iblk) > 0.1d0) Tair(i,j,iblk) = min(Tair(i,j,iblk), 273.25d0)
    end do
  end do
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  Tair(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%q_10(1:nx,1:ny)
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  Qa(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  where(in(1:nx,1:ny,1) .eq. 1)
    tmpdata(1:nx,1:ny) = ocn_grid(1)%pslv(1:nx,1:ny)/((1.d0     \
                        + 0.608d0*ocn_grid(1)%q_10(1:nx,1:ny))  \
                        * 287.04d0*ocn_grid(1)%t_10(1:nx,1:ny))
  end where
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  rhoa(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  pslv_10 = 0.d0
  where(in(1:nx,1:ny,1) .eq. 1)
    pslv_10(1:nx,1:ny) = ocn_grid(1)%pslv(1:nx,1:ny)*dexp(-10.d0  \
                        / (29.3d0*ocn_grid(1)%t_10(1:nx,1:ny)     \
                        * (1.d0+0.608d0*ocn_grid(1)%q_10(1:nx,1:ny))))
    tmpdata(1:nx,1:ny) = ocn_grid(1)%t_10(1:nx,1:ny)*(100000.d0/pslv_10(1:nx,1:ny))**(2.d0/7.d0)
  end where
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  potT(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%snow(1:nx,1:ny)
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  fsnow(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%rain(1:nx,1:ny)
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  frain(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%qflux(1:nx,1:ny)
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  frzmlt(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%u2(1:nx,1:ny,1)*0.01d0
  tmpdata(1:nx,1:ny) = tmpdata(1:nx,1:ny)*ocn_grid(1)%in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  uocn(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%v2(1:nx,1:ny,1)*0.01d0
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  if(symm_np) then
    tmpdata(:,ny+1) = -tmpdata(:,ny+1)
  end if
  vocn(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%t2(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  sst(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  tmpdata(1:nx,1:ny) = ocn_grid(1)%s2(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  sss(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, p0, symm_np)
  tmpdata = 0.d0
  do i = 1, nx
     tmpdata(i,1:ny) = 0.5d0*ocn_grid(1)%odx(1:ny)*(ocn_grid(1)%p0(i+1,1:ny) - ocn_grid(1)%p0(i-1,1:ny))/981.d0
  end do
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  ss_tltx(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)

  tmpdata = 0.d0
  do j = 1, ny
     tmpdata(1:nx,j) = 0.5d0*ocn_grid(1)%ody(j)*(ocn_grid(1)%p0(1:nx,j+1) - ocn_grid(1)%p0(1:nx,j-1))/981.d0
  end do
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, tmpdata, symm_np)
  ss_tlty(:,:,iblk) = tmpdata(0:nx+1,0:ny+1)
   
end subroutine coupler_cice_import

subroutine coupler_cice_export(timer)
  use timcom_general, only: ocn_grid_num, timer_panel
  use timcom_grid,    only: ocn_grid
  use hyperlink,      only: myid, nx, ny

  use ice_blocks, only: nx_block, ny_block
  use ice_state,  only: aice, vice, vsno
  use ice_flux,   only: fresh, fsalt, fhocn

  implicit none

  type(timer_panel) :: timer

  integer :: i, j, k, iblk

  if(myid .eq. 0) write(*,*) "[TIMCOM] CICE to glb coupler: ", timer%mjd1, timer%mjd2

  k = 1
  iblk = 1
  do j = 1, ny
    do i = 1, nx
     ocn_grid(1)%ifrc(i,j)  = min(aice(i+1,j+1,iblk), 1.d0)
     ocn_grid(1)%melt(i,j)  = fresh(i+1,j+1,iblk)*ocn_grid(1)%ifrc(i,j)
     ocn_grid(1)%salt(i,j)  = fsalt(i+1,j+1,iblk)*ocn_grid(1)%ifrc(i,j)
     ocn_grid(1)%melth(i,j) = fhocn(i+1,j+1,iblk)*ocn_grid(1)%ifrc(i,j)
     ocn_grid(1)%vice(i,j)  = vice(i+1,j+1,iblk)
     ocn_grid(1)%vsno(i,j)  = vsno(i+1,j+1,iblk)
    end do
  end do

  ocn_grid(1)%aqice = 0.d0
  ocn_grid(1)%qflux = 0.d0
end subroutine coupler_cice_export
#endif
end module timcom_cplmct
