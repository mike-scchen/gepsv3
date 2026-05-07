module cpl_sendrecv
  use cpl_rank
  use cpl_gsmap
  use cpl_attr
  use cpl_smat
  implicit none
contains
  subroutine cpl_send(id_src, id_dst, remap_Attr, myid)
    use m_MatAttrVectMul, only: MCT_MatVecMul => sMatAvMult
    use m_Transfer,only : MCT_Send => send
    implicit none

    integer, intent(in) :: id_src, id_dst, myid
    type(AttrVect), intent(in), target :: remap_Attr 
    type(AttrVect), pointer :: send_Attr 
    type(SparseMatrixPlus), pointer ::  remp_MatPlus
    type(Router), pointer :: send_Rout
    character(len=256), pointer :: rList

    if(id_src .eq. id_gfs .and. id_dst .eq. id_gocn) then
      remp_MatPlus=>gfs2gocn_sMatP
      send_Attr=>gfs_send_gocn_Attr
      send_Rout=>gfs_send_gocn_Rout

      call MCT_MatVecMul(remap_Attr, remp_MatPlus, send_Attr)
    else if(id_src .eq. id_gocn .and. id_dst .eq. id_gfs) then
      remp_MatPlus=>gocn2gfs_sMatP
      send_Attr=>gocn_send_gfs_Attr
      send_Rout=>gocn_send_gfs_Rout

      call MCT_MatVecMul(remap_Attr, remp_MatPlus, send_Attr)
    else if(id_src .eq. id_rsm .and. id_dst .eq. id_rocn) then
      rsm_send_rocn_Attr%rAttr=remap_Attr%rAttr
      send_Attr=>rsm_send_rocn_Attr
      send_Rout=>rsm_send_rocn_Rout

    else if(id_src .eq. id_rocn .and. id_dst .eq. id_rsm) then
      rocn_send_rsm_Attr%rAttr=remap_Attr%rAttr
      send_Attr=>rocn_send_rsm_Attr
      send_Rout=>rocn_send_rsm_Rout

    else if(id_src .eq. id_gocn .and. id_dst .eq. id_rocn) then
      remp_MatPlus=>gocn2rocn_sMatP
      send_Attr=>gocn_send_rocn_Attr
      send_Rout=>gocn_send_rocn_Rout

      call MCT_MatVecMul(remap_Attr, remp_MatPlus, send_Attr)
    end if

    call MCT_Send(send_Attr, send_Rout)
  end subroutine cpl_send

  subroutine cpl_recv(id_src, id_dst, recv_Attr, myid)
    use m_Transfer,only : MCT_Recv => recv 
    implicit none

    integer, intent(in) :: id_src, id_dst, myid
    type(AttrVect), intent(out) :: recv_Attr
    type(Router), pointer :: recv_Rout

    if(id_src .eq. id_gfs .and. id_dst .eq. id_gocn) then
      recv_Rout=>gfs_recv_gocn_Rout
    else if(id_src .eq. id_gocn .and. id_dst .eq. id_gfs) then
      recv_Rout=>gocn_recv_gfs_Rout
    else if(id_src .eq. id_rsm .and. id_dst .eq. id_rocn) then
      recv_Rout=>rsm_recv_rocn_Rout
    else if(id_src .eq. id_rocn .and. id_dst .eq. id_rsm) then
      recv_Rout=>rocn_recv_rsm_Rout
    else if(id_src .eq. id_rocn .and. id_dst .eq. id_gocn) then
      recv_Rout=>rocn_recv_gocn_Rout
    end if
  
    call MCT_Recv(recv_Attr, recv_Rout)
  end subroutine cpl_recv
end module cpl_sendrecv
