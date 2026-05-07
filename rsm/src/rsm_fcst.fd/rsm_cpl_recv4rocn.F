subroutine rsm_cpl_recv4rocn(compid, slmsk, sst, ssu, ssv)
  use cpl_rank,     only: id_rocn
  use cpl_attr,     only: recv_AV=>rsm_recv_rocn_Attr
  use cpl_sendrecv, only: cpl_recv
#include <paramodel.h>
#include <npesi.h>
#include <comrmpi.h>

  integer :: compid
  real*8, dimension(igrd1p,jgrd1p) :: slmsk, sst, ssu, ssv
  integer :: cnt, i, ibgn, iend, j, jbgn, jend
  
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

  call cpl_recv(compid, id_rocn, recv_AV, mype) 

  ssu = 0.d0
  ssv = 0.d0
  cnt = 0
  do j = jbgn, jend
    do i = ibgn, iend
      cnt = cnt + 1
      if(slmsk(i,j) .le. 0.1d0) then
        sst(i,j) = recv_AV%rAttr(1,cnt)
        ssu(i,j) = recv_AV%rAttr(2,cnt)
        ssv(i,j) = recv_AV%rAttr(3,cnt)
      end if
    end do
  end do

end subroutine rsm_cpl_recv4rocn
