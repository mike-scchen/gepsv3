module tai_timcom_cplmct
  use tai_timcom_const
  implicit none
contains

subroutine ocn_cpl_init(compid, grd)
  use tai_hyperlink, only: nx_grid, nx, ny, &
                       myid_x, myid_y, myid, m_comm_cart
  use tai_timcom_grid, only: ocn_panel
  use cpl_rank,  only: id_rsm
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
  call cpl_gsmap_init(myid, ny, seg_strt, seg_leng, rank_root, m_comm_cart, compid, "Tai GSMap:")
  call cpl_smat_init(myid, rank_root, compid, id_rsm, m_comm_cart, "SameGrid")
  call cpl_attr_init(compid, id_rsm, m_comm_cart)

end subroutine ocn_cpl_init

subroutine coupler_data_exch(compid, timer)
  use tai_timcom_general, only: ocn_grid_num, timer_panel
  use tai_timcom_grid,    only: ocn_grid
  use tai_hyperlink,      only: myid
  implicit none

  integer, intent(in) :: compid
  type(timer_panel) :: timer

  real(r8), pointer, dimension(:,:) :: &
    sst, ssu, ssv, &
    u10m, v10m, t02m, q02m, pslv, &
    swup, swdn, lwdn, rain, snow, tsea

  integer(2), pointer, dimension(:,:) :: &
    msk

  integer :: n, nx, ny

  if(myid .eq. 0) write(*,*) "[TIMCOM]: Tai to coupler: ", timer%mjd1, timer%mjd2
 
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
    tsea => ocn_grid(n)%tsea
  
    call ocn_cpl_send2rsm(compid, sst, ssu, ssv, msk)
    call ocn_cpl_recv4rsm(compid, u10m, v10m, t02m, q02m, pslv, &
                            swup, swdn, lwdn, rain, snow, tsea) 
  end do 

end subroutine coupler_data_exch

subroutine ocn_cpl_send2rsm(compid, SST, SSU, SSV, msk)
  use tai_hyperlink, only: nx, ny, myid
  use cpl_rank, only: id_rsm
  use cpl_attr, only: AttrVect_importRAttr, rocn_export_attr, rocnExport_rList 
  use cpl_sendrecv, only: cpl_send
  implicit none

  integer, intent(in) :: compid
  real(r8), dimension(nx*ny), target, optional, intent(in)  :: SST, SSU, SSV
  integer(2), dimension(nx*ny), intent(in) :: msk
  real(r8), dimension(:), pointer :: send_ptr
  integer :: temp

  if(present(SST)) then
    send_ptr=>SST
    call AttrVect_importRAttr(rocn_export_attr, "sst", send_ptr)
  end if
  if(present(SSU)) then
    send_ptr=>SSU
    call AttrVect_importRAttr(rocn_export_attr, "ssu", send_ptr)
  end if
  if(present(SSV)) then
    send_ptr=>SSV
    call AttrVect_importRAttr(rocn_export_attr, "ssv", send_ptr)
  end if
  rocn_export_attr%rAttr(1,:) = rocn_export_attr%rAttr(1,:) + 273.15d0
  rocn_export_attr%rAttr(2,:) = rocn_export_attr%rAttr(2,:)/100.d0
  rocn_export_attr%rAttr(3,:) = rocn_export_attr%rAttr(3,:)/100.d0
  temp = size(SST)
  !write(*,*) 'yc check cpl_send@ocn:', myid, temp
  !if(myid .eq. 192) write(*,*) "yc check ocn send in TIMCOM:", SST(100), SSU(100), SSV(100), msk(100)
  call cpl_send(compid, id_rsm, rocn_export_attr, myid)

end subroutine ocn_cpl_send2rsm

subroutine ocn_cpl_recv4rsm(compid,u10m,v10m,t02m,q02m,pslv,swup,swdn,lwdn,rain,snow,tsea)
  use tai_hyperlink, only: nx, ny, myid
  use cpl_rank,  only: id_rsm
  use cpl_attr,  only: rocn_recv_rsm_Attr
  use cpl_sendrecv, only: cpl_recv
  implicit none

  integer, intent(in) :: compid
  real(r8), dimension(nx*ny), intent(out) :: u10m, v10m, t02m, q02m, &
                                                 pslv, swup, swdn, lwdn, &
                                                 rain, snow, tsea

  call cpl_recv(compid, id_rsm, rocn_recv_rsm_Attr, myid)
  u10m = rocn_recv_rsm_Attr%rAttr(1,:)
  v10m = rocn_recv_rsm_Attr%rAttr(2,:)
  t02m = rocn_recv_rsm_Attr%rAttr(3,:)
  q02m = rocn_recv_rsm_Attr%rAttr(4,:)
  pslv = rocn_recv_rsm_Attr%rAttr(5,:)
  swup = rocn_recv_rsm_Attr%rAttr(6,:)
  swdn = rocn_recv_rsm_Attr%rAttr(7,:)
  lwdn = rocn_recv_rsm_Attr%rAttr(8,:)
  rain = rocn_recv_rsm_Attr%rAttr(9,:)
  tsea = rocn_recv_rsm_Attr%rAttr(10,:)
  snow = 0.d0

end subroutine ocn_cpl_recv4rsm
end module tai_timcom_cplmct
