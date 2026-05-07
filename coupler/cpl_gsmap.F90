module cpl_gsmap
  use m_GlobalSegMap, only: GlobalSegMap
  use m_GlobalSegMap, only: GlobalSegMap_init => init
  use m_GlobalSegMap, only: GlobalSegMap_lsize => lsize
  use m_GlobalSegMap, only: GlobalSegMap_clean => clean
  use cpl_rank
  implicit none 
 
  type(GlobalSegMap), target ::         gfsGSMap
  type(GlobalSegMap), target :: gocn_in_gfsGSMap
  type(GlobalSegMap), target ::  rsm_in_gfsGSMap

  type(GlobalSegMap), target ::         rsmGSMap
  type(GlobalSegMap), target :: rocn_in_rsmGSMap
  type(GlobalSegMap), target ::  gfs_in_rsmGSMap

  type(GlobalSegMap), target ::         gocnGSMap
  type(GlobalSegMap), target ::  gfs_in_gocnGSMap
  type(GlobalSegMap), target :: rocn_in_gocnGSMap

  type(GlobalSegMap), target ::         rocnGSMap
  type(GlobalSegMap), target ::  rsm_in_rocnGSMap
  type(GlobalSegMap), target :: gocn_in_rocnGSMap

contains
subroutine cpl_gsmap_init(myrank, seg_size, seg_strt, seg_leng, rank_root, mct_comm_world, id_comp, msg)
  implicit none

  integer, intent(in) :: myrank, seg_size
  integer, intent(in) :: seg_strt(seg_size), seg_leng(seg_size)
  integer, intent(in) :: rank_root
  integer, intent(in) :: mct_comm_world
  integer, intent(in) :: id_comp
  character(len=*) :: msg
  integer :: seg_strt1(1), seg_leng1(1), res

  if(id_comp.eq.id_gfs) then
    call GlobalSegMap_init(gfsGSMap, seg_strt, seg_leng, rank_root, mct_comm_world, id_comp)
    avsize(id_comp) = GlobalSegMap_lsize(gfsGSMap, mct_comm_world)

    if(l_comp_act(id_gocn)) then
      res = mod(ngrid_comp_all(id_gocn), nproc_comp_all(id_gfs))
      seg_leng1(1) = (ngrid_comp_all(id_gocn)-res)/nproc_comp_all(id_gfs)
      if(myrank .lt. res) then
        seg_leng1(1) = seg_leng1(1) + 1
        seg_strt1(1) = (myrank*seg_leng1(1)) + 1
      else
        seg_strt1(1) = res*(seg_leng1(1)+1) + (myrank - res)*seg_leng1(1) + 1
      end if
      call GlobalSegMap_init(gocn_in_gfsGSMap, seg_strt1, seg_leng1, rank_root, mct_comm_world, id_comp)
      avsize(id_gocn) = GlobalSegMap_lsize(gocn_in_gfsGSMap, mct_comm_world)
    end if

!    if(l_comp_act(id_rsm)) then
!      seg_leng1(1) = ngrid_comp_all(id_rsm)/nproc_comp_all(id_gfs)
!      seg_strt1(1) = (myrank*seg_leng1(1)) + 1
!      call GlobalSegMap_init(rsm_in_gfsGSMap, seg_strt1, seg_leng1, rank_root, mct_comm_world, id_comp)
!      avsize(id_rsm) = GlobalSegMap_lsize(rsm_in_gfsGSMap, mct_comm_world)
!      write(*,*) trim(msg), ", myrank:", myrank, ", build  rsm in gfs GSMap", avsize(id_rsm)
!    end if

  elseif(id_comp.eq.id_rsm) then
    call GlobalSegMap_init(rsmGSMap, seg_strt, seg_leng, rank_root, mct_comm_world, id_comp)
    avsize(id_comp) = GlobalSegMap_lsize(rsmGSMap, mct_comm_world)

!    if(l_comp_act(id_gfs)) then
!      seg_leng1(1) = ngrid_comp_all(id_gfs)/nproc_comp_all(id_rsm)
!      seg_strt1(1) = (myrank*seg_leng1(1)) + 1
!      call GlobalSegMap_init(gfs_in_rsmGSMap, seg_strt1, seg_leng1, rank_root, mct_comm_world, id_comp)
!      avsize(id_gfs) = GlobalSegMap_lsize(gfs_in_rsmGSMap, mct_comm_world)
!      write(*,*) trim(msg), ", myrank:", myrank, ", build  gfs in rsm GSMap", avsize(id_gfs)
!    end if

    if(l_comp_act(id_rocn)) then
      res = mod(ngrid_comp_all(id_rocn), nproc_comp_all(id_rsm))
      seg_leng1(1) = (ngrid_comp_all(id_rocn)-res)/nproc_comp_all(id_rsm)
      if(myrank .lt. res) then
        seg_leng1(1) = seg_leng1(1) + 1
        seg_strt1(1) = (myrank*seg_leng1(1)) + 1
      else
        seg_strt1(1) = res*(seg_leng1(1)+1) + (myrank - res)*seg_leng1(1) + 1
      end if      
      call GlobalSegMap_init(rocn_in_rsmGSMap, seg_strt, seg_leng, rank_root, mct_comm_world, id_comp)
      avsize(id_rocn) = GlobalSegMap_lsize(rocn_in_rsmGSMap, mct_comm_world)
    end if

  elseif(id_comp.eq.id_gocn) then
    call GlobalSegMap_init(gocnGSMap, seg_strt, seg_leng, rank_root, mct_comm_world, id_comp)
    avsize(id_comp) = GlobalSegMap_lsize(gocnGSMap, mct_comm_world)

    if(l_comp_act(id_gfs)) then
      res = mod(ngrid_comp_all(id_gfs), nproc_comp_all(id_gocn))
      seg_leng1(1) = (ngrid_comp_all(id_gfs)-res)/nproc_comp_all(id_gocn)
      if(myrank .lt. res) then
        seg_leng1(1) = seg_leng1(1) + 1
        seg_strt1(1) = (myrank*seg_leng1(1)) + 1
      else
        seg_strt1(1) = res*(seg_leng1(1)+1) + (myrank - res)*seg_leng1(1) + 1
      end if
      call GlobalSegMap_init(gfs_in_gocnGSMap, seg_strt1, seg_leng1, rank_root, mct_comm_world, id_comp)
      avsize(id_gfs) = GlobalSegMap_lsize(gfs_in_gocnGSMap, mct_comm_world)
    end if

    if(l_comp_act(id_rocn)) then
      res = mod(ngrid_comp_all(id_rocn), nproc_comp_all(id_gocn))
      seg_leng1(1) = (ngrid_comp_all(id_rocn)-res)/nproc_comp_all(id_gocn)
      if(myrank .lt. res) then
        seg_leng1(1) = seg_leng1(1) + 1
        seg_strt1(1) = (myrank*seg_leng1(1)) + 1
      else
        seg_strt1(1) = res*(seg_leng1(1)+1) + (myrank - res)*seg_leng1(1) + 1
      end if

      call GlobalSegMap_init(rocn_in_gocnGSMap, seg_strt1, seg_leng1, rank_root, mct_comm_world, id_comp)
      avsize(id_rocn) = GlobalSegMap_lsize(rocn_in_gocnGSMap, mct_comm_world)
!      write(*,*) trim(msg), ", myrank:", myrank, ", build rocn in gocn GSMap", avsize(id_rocn), seg_leng1, seg_strt1, ngrid_comp_all(id_rocn), nproc_comp_all(id_gocn)
    end if

  elseif(id_comp.eq.id_rocn) then
    call GlobalSegMap_init(rocnGSMap, seg_strt, seg_leng, rank_root, mct_comm_world, id_comp)
    avsize(id_comp) = GlobalSegMap_lsize(rocnGSMap, mct_comm_world)

!    if(l_comp_act(id_gocn)) then
!      seg_leng1(1) = ngrid_comp_all(id_gocn)/nproc_comp_all(id_rocn)
!      seg_strt1(1) = (myrank*seg_leng1(1)) + 1
!      call GlobalSegMap_init(gocn_in_rocnGSMap, seg_strt1, seg_leng1, rank_root, mct_comm_world, id_comp)
!      avsize(id_gocn) = GlobalSegMap_lsize(gocn_in_rocnGSMap, mct_comm_world)
!      write(*,*) trim(msg), ", myrank:", myrank, ", build gocn in rocn GSMap", avsize(id_gocn)
!    end if

    if(l_comp_act(id_rsm)) then
      res = mod(ngrid_comp_all(id_rsm), nproc_comp_all(id_rocn))
      seg_leng1(1) = (ngrid_comp_all(id_rsm)-res)/nproc_comp_all(id_rocn)
      if(myrank .lt. res) then
        seg_leng1(1) = seg_leng1(1) + 1
        seg_strt1(1) = (myrank*seg_leng1(1)) + 1
      else
        seg_strt1(1) = res*(seg_leng1(1)+1) + (myrank - res)*seg_leng1(1) + 1
      end if
      call GlobalSegMap_init(rsm_in_rocnGSMap, seg_strt, seg_leng, rank_root, mct_comm_world, id_comp)
      avsize(id_rsm) = GlobalSegMap_lsize(rsm_in_rocnGSMap, mct_comm_world)
    end if

  end if
end subroutine cpl_gsmap_init

end module cpl_gsmap
