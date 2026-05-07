module cpl_attr
  use m_AttrVect, only: AttrVect
  use m_AttrVect, only: AttrVect_init => init
  use m_AttrVect, only: AttrVect_importRAttr => importRAttr
  use m_AttrVect, only: mct_aVect_indexIA => indexIA
  use m_AttrVect, only: mct_aVect_indexRA => indexRA
  use m_Router,   only: Router
  use m_Router,   only: Router_init => init
  use cpl_rank
  implicit none
                  
  type(AttrVect), target :: gfs_export_Attr 
  type(AttrVect), target :: gfs_send_gocn_Attr
  type(AttrVect), target :: gfs_recv_gocn_Attr
  type(AttrVect), target :: gfs_nested_Attr
  type(AttrVect), target :: gfs_send_rsm_Attr

  type(AttrVect), target :: rsm_export_Attr 
  type(AttrVect), target :: rsm_send_rocn_Attr
  type(AttrVect), target :: rsm_recv_rocn_Attr
  type(AttrVect), target :: rsm_nested_Attr
  type(AttrVect), target :: rsm_recv_gfs_Attr

  type(AttrVect), target :: gocn_export_Attr
  type(AttrVect), target :: gocn_send_gfs_Attr
  type(AttrVect), target :: gocn_recv_gfs_Attr
  type(AttrVect), target :: gocn_nested_Attr
  type(AttrVect), target :: gocn_send_rocn_Attr

  type(AttrVect), target :: rocn_export_Attr
  type(AttrVect), target :: rocn_send_rsm_Attr
  type(AttrVect), target :: rocn_recv_rsm_Attr
  type(AttrVect), target :: rocn_nested_Attr
  type(AttrVect), target :: rocn_recv_gocn_Attr

  type(Router), target :: gfs_send_gocn_Rout
  type(Router), target :: gfs_recv_gocn_Rout
  type(Router), target :: gfs_send_rsm_Rout

  type(Router), target :: rsm_send_rocn_Rout
  type(Router), target :: rsm_recv_rocn_Rout
  type(Router), target :: rsm_recv_gfs_Rout

  type(Router), target :: gocn_send_gfs_Rout
  type(Router), target :: gocn_recv_gfs_Rout
  type(Router), target :: gocn_send_rocn_Rout

  type(Router), target :: rocn_send_rsm_Rout
  type(Router), target :: rocn_recv_rsm_Rout
  type(Router), target :: rocn_recv_gocn_Rout

  character(len=256), target ::  gfsExport_rList = "u10m:v10m:t02m:q02m:pslv:swup:swdn:lwdn:rain:snow:tg"
  character(len=256), target ::  rsmExport_rList = "u10m:v10m:t02m:q02m:pslv:swup:swdn:lwdn:rain:tsea"
  character(len=256), target :: gocnExport_rList = "sst:ssu:ssv"
  character(len=256), target :: rocnExport_rList = "sst:ssu:ssv"

  character(len=256), target ::  gfsNested_rList = "field1:field2"
  character(len=256), target ::  rsmNested_rList = "field1:field2"
  character(len=256), target :: gocnNested_rList = "ocn_u:ocn_v:ocn_t:ocn_s"
  character(len=256), target :: rocnNested_rList = "ocn_u:ocn_v:ocn_t:ocn_s"


  character(len=256), target ::   gfs2gocn_rList = "u10m:v10m:t02m:q02m:pslv:swup:swdn:lwdn:rain:snow:tg"
  character(len=256), target ::   gocn2gfs_rList = "sst:ssu:ssv"
  character(len=256), target ::    gfs2rsm_rList = "field1:field2"
  character(len=256), target ::  gocn2rocn_rList = "ocn_u:ocn_v:ocn_t:ocn_s"
  character(len=256), target ::   rsm2rocn_rList = "u10m:v10m:t02m:q02m:pslv:swup:swdn:lwdn:rain:tsea"
  character(len=256), target ::   rocn2rsm_rList = "sst:ssu:ssv"
contains

subroutine cpl_attr_init(id_src, id_dst, mpi_comm_mct)
!  use m_GlobalSegMap,only: GlobalSegMap
  use cpl_gsmap 
  implicit none

  integer, intent(in) :: id_src 
  integer, intent(in) :: id_dst
  integer, intent(in) :: mpi_comm_mct
  integer :: send_size, recv_size, export_size
  character(len=256), pointer :: send_list, recv_list, export_list
  type(GlobalSegMap), pointer :: dstGSMap, srcGSMap
  type(AttrVect), pointer :: send_attr, recv_attr, export_attr
  type(Router), pointer :: send_rout, recv_rout

  if(id_src .eq. id_gfs .and. id_dst .eq. id_gocn) then
    srcGSMap=>gfsGSMap
    dstGSMap=>gocn_in_gfsGSMap

    export_attr=>gfs_export_Attr
    export_list=>gfsExport_rList
    export_size=avsize(id_src)

    send_attr=>gfs_send_gocn_Attr
    send_list=>gfs2gocn_rList
    send_size=avsize(id_dst)
    send_rout=>gfs_send_gocn_Rout

    recv_attr=>gfs_recv_gocn_Attr
    recv_list=>gocn2gfs_rList
    recv_size=avsize(id_src)
    recv_rout=>gfs_recv_gocn_Rout

    call Router_init(id_dst, dstGSMap, mpi_comm_mct, send_Rout)  
    !write(*,'(a14,8i8)') 'AttrVect send:', id_src, id_dst, send_size, &
    !                  send_Rout%comp1id, send_Rout%comp2id, &
    !                  send_Rout%nprocs, send_Rout%maxsize, send_Rout%lavsize
    call Router_init(id_dst, srcGSMap, mpi_comm_mct, recv_Rout)
    !write(*,'(a14,8i8)') 'AttrVect recv:', id_src, id_dst, recv_size, &
    !                  recv_Rout%comp1id, recv_Rout%comp2id, &
    !                  recv_Rout%nprocs, recv_Rout%maxsize, recv_Rout%lavsize
 
  elseif(id_src .eq. id_gocn .and. id_dst .eq. id_gfs) then
    srcGSMap=>gocnGSMap
    dstGSMap=>gfs_in_gocnGSMap

    export_attr=>gocn_export_Attr
    export_list=>gocnExport_rList
    export_size=avsize(id_src)

    send_attr=>gocn_send_gfs_Attr
    send_list=>gocn2gfs_rList
    send_size=avsize(id_dst)
    send_rout=>gocn_send_gfs_Rout

    recv_attr=>gocn_recv_gfs_Attr
    recv_list=>gfs2gocn_rList
    recv_size=avsize(id_src)
    recv_rout=>gocn_recv_gfs_Rout

    call Router_init(id_dst, srcGSMap, mpi_comm_mct, recv_Rout)

    call Router_init(id_dst, dstGSMap, mpi_comm_mct, send_Rout)

  elseif(id_src .eq. id_gfs .and. id_dst .eq. id_rsm) then
    srcGSMap=>gfsGSMap
    dstGSMap=>rsm_in_gfsGSMap

    export_attr=>gfs_nested_Attr
    export_list=>gfsNested_rList
    export_size=avsize(id_src)

    send_attr=>gfs_send_rsm_Attr
    send_list=>gfs2rsm_rList
    send_size=avsize(id_dst)
    send_rout=>gfs_send_rsm_Rout

    call Router_init(id_dst, dstGSMap, mpi_comm_mct, send_Rout)
  
  elseif(id_src .eq. id_rsm .and. id_dst .eq. id_gfs) then
    srcGSMap=>rsmGSMap
    dstGSMap=>gfs_in_rsmGSMap

    export_attr=>rsm_nested_Attr
    export_list=>rsmNested_rList
    export_size=avsize(id_src)

    recv_attr=>rsm_recv_gfs_Attr
    recv_list=>gfs2rsm_rList
    recv_size=avsize(id_src)
    recv_rout=>rsm_recv_gfs_Rout

    call Router_init(id_dst, srcGSMap, mpi_comm_mct, recv_Rout)
  
  elseif(id_src .eq. id_gocn .and. id_dst .eq. id_rocn) then
    srcGSMap=>gocnGSMap
    dstGSMap=>rocn_in_gocnGSMap

    export_attr=>gocn_nested_Attr
    export_list=>gocnNested_rList
    export_size=avsize(id_src)

    send_attr=>gocn_send_rocn_Attr
    send_list=>gocn2rocn_rList
    send_size=avsize(id_dst)
    send_rout=>gocn_send_rocn_Rout

    call Router_init(id_dst, dstGSMap, mpi_comm_mct, send_Rout)

  elseif(id_src .eq. id_rocn .and. id_dst .eq. id_gocn) then
    srcGSMap=>rocnGSMap
    dstGSMap=>gocn_in_rocnGSMap

    export_attr=>rocn_nested_Attr
    export_list=>rocnNested_rList
    export_size=avsize(id_src)

    recv_attr=>rocn_recv_gocn_Attr
    recv_list=>gocn2rocn_rList
    recv_size=avsize(id_src)
    recv_rout=>rocn_recv_gocn_Rout

    call Router_init(id_dst, srcGSMap, mpi_comm_mct, recv_Rout)

  elseif(id_src .eq. id_rsm .and. id_dst .eq. id_rocn) then
    srcGSMap=>rsmGSMap
    dstGSMap=>rocn_in_rsmGSMap

    export_attr=>rsm_export_Attr
    export_list=>rsmExport_rList
    export_size=avsize(id_src)

    send_attr=>rsm_send_rocn_Attr
    send_list=>rsm2rocn_rList
    send_size=avsize(id_dst)
    send_rout=>rsm_send_rocn_Rout

    recv_attr=>rsm_recv_rocn_Attr
    recv_list=>rocn2rsm_rList
    recv_size=avsize(id_src)
    recv_rout=>rsm_recv_rocn_Rout

    call Router_init(id_dst, dstGSMap, mpi_comm_mct, send_Rout)
    call Router_init(id_dst, srcGSMap, mpi_comm_mct, recv_Rout)

  elseif(id_src .eq. id_rocn .and. id_dst .eq. id_rsm) then
    srcGSMap=>rocnGSMap
    dstGSMap=>rsm_in_rocnGSMap

    export_attr=>rocn_export_Attr
    export_list=>rocnExport_rList
    export_size=avsize(id_src)

    send_attr=>rocn_send_rsm_Attr
    send_list=>rocn2rsm_rList
    send_size=avsize(id_dst)
    send_rout=>rocn_send_rsm_Rout

    recv_attr=>rocn_recv_rsm_Attr
    recv_list=>rsm2rocn_rList
    recv_size=avsize(id_src)
    recv_rout=>rocn_recv_rsm_Rout
  
    call Router_init(id_dst, srcGSMap, mpi_comm_mct, recv_Rout)
    call Router_init(id_dst, dstGSMap, mpi_comm_mct, send_Rout)

  end if

  if(.not. associated(export_attr)) write(*,*) "MCT ATTR ERROR01:", id_src, id_dst
  call AttrVect_init(export_attr, rList=export_list, lsize=export_size)
  if(associated(send_Attr)) then
    if(.not. associated(export_attr)) write(*,*) "MCT ATTR ERROR02:", id_src, id_dst
    call AttrVect_init(send_attr, rList=send_list, lsize=send_size)
  end if
  if(associated(recv_Attr)) then
    if(.not. associated(export_attr)) write(*,*) "MCT ATTR ERROR03:", id_src, id_dst
    call AttrVect_init(recv_attr, rList=recv_list, lsize=recv_size)
  end if
end subroutine cpl_attr_init

end module cpl_attr
