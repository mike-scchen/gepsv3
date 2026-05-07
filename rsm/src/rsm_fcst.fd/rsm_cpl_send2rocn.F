subroutine rsm_cpl_send2rocn(compid, cpl_field, tsea)
  use cpl_rank,     only: id_rocn
  use cpl_attr,     only: AttrVect_importRAttr, rsm_export_attr, rsmExport_rList
  use cpl_sendrecv, only: cpl_send
#include <paramodel.h>
#include <npesi.h>
#include <comrmpi.h>

  integer, intent(in) :: compid
  real*8 :: cpl_field(igrd1p,jgrd1p,9), tsea(igrd1p,jgrd1p)
  integer :: cnt, i, ibgn, iend, j, jbgn, jend 

  cnt = 0
  ibgn = 1
  iend = igrd1p-1
  jbgn = 1
  jend = jgrd1p-1

  if(mycol .eq. 0) then
    ibgn = 2
    iend = igrd1p
  end if
  if(myrow .eq. 0) then
    jbgn = 2
    jend = jgrd1p
  end if


  do j = jbgn, jend
    do i = ibgn, iend
      cnt = cnt + 1
      rsm_export_attr%rAttr( 1,cnt) = cpl_field(i,j,1)
      rsm_export_attr%rAttr( 2,cnt) = cpl_field(i,j,2)
      rsm_export_attr%rAttr( 3,cnt) = cpl_field(i,j,3)
      rsm_export_attr%rAttr( 4,cnt) = cpl_field(i,j,4)
      rsm_export_attr%rAttr( 5,cnt) = cpl_field(i,j,5)
      rsm_export_attr%rAttr( 6,cnt) = cpl_field(i,j,6)
      rsm_export_attr%rAttr( 7,cnt) = cpl_field(i,j,7)
      rsm_export_attr%rAttr( 8,cnt) = cpl_field(i,j,8)
      rsm_export_attr%rAttr( 9,cnt) = cpl_field(i,j,9)
      rsm_export_attr%rAttr(10,cnt) = tsea(i,j)
    end do
  enddo

  call cpl_send(compid, id_rocn, rsm_export_attr, mype)
end subroutine rsm_cpl_send2rocn

