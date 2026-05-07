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
  character(len=256), target :: gocnExport_rList = "sst:ssu:ssv:ifrc:vice:vsno"
  character(len=256), target :: rocnExport_rList = "sst:ssu:ssv"

  character(len=256), target ::  gfsNested_rList = "field1:field2"
  character(len=256), target ::  rsmNested_rList = "field1:field2"
  character(len=256), target :: gocnNested_rList = "ocn_u:ocn_v:ocn_t:ocn_s"
  character(len=256), target :: rocnNested_rList = "ocn_u:ocn_v:ocn_t:ocn_s"


  character(len=256), target ::   gfs2gocn_rList = "u10m:v10m:t02m:q02m:pslv:swup:swdn:lwdn:rain:snow:tg"
  character(len=256), target ::   gocn2gfs_rList = "sst:ssu:ssv:ifrc:vice:vsno"
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

    call AttrVect_init(gfs_export_Attr,    rList=gfsExport_rList, lsize=avsize(id_src))
    call AttrVect_init(gfs_send_gocn_Attr, rList=gfs2gocn_rList,  lsize=avsize(id_dst))
    call AttrVect_init(gfs_recv_gocn_Attr, rList=gocn2gfs_rList,  lsize=avsize(id_src))
    call Router_init(id_dst, gocn_in_gfsGSMap, mpi_comm_mct, gfs_send_gocn_Rout)  
    call Router_init(id_dst,         gfsGSMap, mpi_comm_mct, gfs_recv_gocn_Rout)
 
  elseif(id_src .eq. id_gocn .and. id_dst .eq. id_gfs) then

    call AttrVect_init(gocn_export_Attr,   rList=gocnExport_rList, lsize=avsize(id_src))
    call AttrVect_init(gocn_send_gfs_Attr, rList=gocn2gfs_rList,   lsize=avsize(id_dst))
    call AttrVect_init(gocn_recv_gfs_Attr, rList=gfs2gocn_rList,   lsize=avsize(id_src))
    call Router_init(id_dst,        gocnGSMap, mpi_comm_mct, gocn_recv_gfs_Rout)
    call Router_init(id_dst, gfs_in_gocnGSMap, mpi_comm_mct, gocn_send_gfs_Rout)

  elseif(id_src .eq. id_gfs .and. id_dst .eq. id_rsm) then

    call AttrVect_init(gfs_nested_Attr,    rList=gfsNested_rList, lsize=avsize(id_src))
    call AttrVect_init(gfs_send_rsm_Attr,  rList=gfs2rsm_rList,   lsize=avsize(id_dst))
    call Router_init(id_dst, rsm_in_gfsGSMap, mpi_comm_mct, gfs_send_rsm_Rout)
  
  elseif(id_src .eq. id_rsm .and. id_dst .eq. id_gfs) then

    call AttrVect_init(rsm_nested_Attr,    rList=rsmNested_rList, lsize=avsize(id_src))
    call AttrVect_init(rsm_recv_gfs_Attr,  rList=gfs2rsm_rList,   lsize=avsize(id_src))
    call Router_init(id_dst, rsmGSMap, mpi_comm_mct, rsm_recv_gfs_Rout)

  elseif(id_src .eq. id_gocn .and. id_dst .eq. id_rocn) then

    call AttrVect_init(gocn_nested_Attr,    rList=gocnNested_rList, lsize=avsize(id_src))
    call AttrVect_init(gocn_send_rocn_Attr, rList=gocn2rocn_rList,  lsize=avsize(id_dst))
    call Router_init(id_dst, rocn_in_gocnGSMap, mpi_comm_mct, gocn_send_rocn_Rout)

  elseif(id_src .eq. id_rocn .and. id_dst .eq. id_gocn) then

    call AttrVect_init(rocn_nested_Attr,    rList=rocnNested_rList, lsize=avsize(id_src))
    call AttrVect_init(rocn_recv_gocn_Attr, rList=gocn2rocn_rList,  lsize=avsize(id_src))
    call Router_init(id_dst, rocnGSMap, mpi_comm_mct, rocn_recv_gocn_Rout)

  elseif(id_src .eq. id_rsm .and. id_dst .eq. id_rocn) then

    call AttrVect_init(rsm_export_Attr,    rList=rsmExport_rList, lsize=avsize(id_src))
    call AttrVect_init(rsm_send_rocn_Attr, rList=rsm2rocn_rList,  lsize=avsize(id_dst))
    call AttrVect_init(rsm_recv_rocn_Attr, rList=rocn2rsm_rList,  lsize=avsize(id_src))
    call Router_init(id_dst, rocn_in_rsmGSMap, mpi_comm_mct, rsm_send_rocn_Rout)
    call Router_init(id_dst,         rsmGSMap, mpi_comm_mct, rsm_recv_rocn_Rout)

  elseif(id_src .eq. id_rocn .and. id_dst .eq. id_rsm) then

    call AttrVect_init(rocn_export_Attr,   rList=rocnExport_rList, lsize=avsize(id_src))
    call AttrVect_init(rocn_send_rsm_Attr, rList=rocn2rsm_rList,   lsize=avsize(id_dst))
    call AttrVect_init(rocn_recv_rsm_Attr, rList=rsm2rocn_rList,   lsize=avsize(id_src))
    call Router_init(id_dst,        rocnGSMap, mpi_comm_mct, rocn_recv_rsm_Rout)
    call Router_init(id_dst, rsm_in_rocnGSMap, mpi_comm_mct, rocn_send_rsm_Rout)
  
  end if

end subroutine cpl_attr_init

end module cpl_attr
