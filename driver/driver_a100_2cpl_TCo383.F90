PROGRAM driver
  USE cpl_rank, only: n_comp_all, id_gfs, id_rsm, id_gocn, id_rocn, &
                      l_comp_act, nproc_comp_all, myrank_root_all, &
                      ngrid_comp_all
  USE MPI
  implicit none

  integer :: nprocs, myrank, id_comm
  integer :: nprocs1, myrank1, comm1, comm2
  integer :: n_comps, id_comp(n_comp_all), nproc_comp(0:n_comp_all)
  integer :: n, ierr, n0_cnt, n1_cnt,  split_comm
  character(len=512) :: env_tmp
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

  id_comp(id_gfs)  = id_gfs
  id_comp(id_rsm)  = id_rsm
  id_comp(id_gocn) = id_gocn
  id_comp(id_rocn) = id_rocn
  nproc_comp = -1
 
  if(myrank .eq. 0) then 
    call get_env_info('n_proc_gfs',  nproc_comp(id_gfs))
    call get_env_info('n_proc_rsm',  nproc_comp(id_rsm))
    call get_env_info('n_proc_gocn', nproc_comp(id_gocn))
    call get_env_info('n_proc_rocn', nproc_comp(id_rocn))
  end if
  call MPI_BCAST(nproc_comp, 5, MPI_INT, 0, MPI_COMM_WORLD, ierr)

  if(sum(nproc_comp(1:n_comp_all)) .ne. nprocs) then
    if(myrank .eq. 0) &
      write(*,*) 'ERROR(driver): process inbalance', nprocs, sum(nproc_comp(1:n_comp_all))
    call MPI_FINALIZE(ierr)
    stop
  else
    if(myrank .eq. 0) write(*,'(a19,i4)') "total process used:", nprocs
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
  end if
  n_comps = count(nproc_comp.gt.0)
  if(myrank.eq.0) write(*,*) "before:", id_comp, nproc_comp(1:n_comp_all)
  id_comm = -1
  n0_cnt = 1
  n1_cnt = 1
  do n = 1, n_comp_all
    if(nproc_comp(n) .eq. 0) then
      id_comp(n) = n_comps + n0_cnt
      n0_cnt = n0_cnt + 1
    else
      id_comp(n) = n1_cnt
      n1_cnt = n1_cnt + 1
    end if
    nproc_comp_all(id_comp(n)) = nproc_comp(n)
    nproc_comp(n) = nproc_comp(n) + nproc_comp(n-1)
    if(myrank.gt.nproc_comp(n-1) .and. myrank.le.nproc_comp(n)) then
      id_comm = id_comp(n)
      myrank_root_all = nproc_comp(n-1)+1
    end if
  end do

  if(myrank.eq.0) write(*,*) "after:", id_comp, nproc_comp_all
  id_gfs  = id_comp(id_gfs)
  id_rsm  = id_comp(id_rsm)
  id_gocn = id_comp(id_gocn)
  id_rocn = id_comp(id_rocn)

  l_comp_act = .false. 
  do n = 1, n_comps
    l_comp_act(n) = .true.   
  end do
  ngrid_comp_all(id_gfs)  = 1552*768
  ngrid_comp_all(id_gocn) = 1536*720
  ngrid_comp_all(id_rsm)  = 768*432
  ngrid_comp_all(id_rocn) = 768*432

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
 
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, id_comm, id_comp(id_comm), split_comm, ierr)
  call MPI_COMM_SIZE(split_comm, nprocs1, ierr)
  call MPI_COMM_RANK(split_comm, myrank1, ierr)
  if(myrank1 .eq. nproc_comp_all(id_comm)-1) then
    write(*,*) "check id_comm: ", n_comps, id_comm, myrank, nproc_comp_all(id_comm), ngrid_comp_all(id_comm), split_comm
  end if
#ifdef check_info
  write(*,'(a7,3i4,a9,2i4,a7,2i4,a10,i4,a13,2i4)') &
          "myrank:", myrank, myrank1, myrank_root_all, ", nprocs:", nprocs, nprocs1, ", comm:", MPI_COMM_WORLD, split_comm, &
       ", id_comp:", id_comm, ", rank group:", nproc_comp(n-1)+1, nproc_comp(n)
#endif
  if(id_comm .eq. id_gfs) then
    call cpl_init(n_comps, split_comm, id_gfs)
    call MPI_BARRIER(split_comm, ierr)
    if(myrank1 .eq. 0) write(*,*) 'run TCo at: ', MPI_WTIME()
    write(*,*) 'check split_comm=',split_comm
    call gfcst(id_gfs, split_comm)
    if(myrank1 .eq. 0) write(*,*) 'gfs model done'

  elseif(id_comm .eq. id_rsm) then
    call cpl_init(n_comps, split_comm, id_rsm)
    call getenv('RUNDIR', env_tmp)
    call chdir(trim(env_tmp))
    call MPI_BARRIER(split_comm, ierr)
!    if(myrank1 .eq. 0) write(*,*) 'run RSM at: ', MPI_WTIME()
!    call mainrsmf_mrg(id_rsm, split_comm) 
!    if(myrank1 .eq. 0) write(*,*) 'rsm model done'

  elseif(id_comm .eq. id_gocn) then
    call cpl_init(n_comps, split_comm, id_gocn)
    call MPI_BARRIER(split_comm, ierr)
    if(myrank1 .eq. 0) write(*,*) 'run GLB at: ', MPI_WTIME()
    call timcom_main_gpu(id_gocn, split_comm)
    if(myrank1 .eq. 0) write(*,*) 'global oceann model done'

  elseif(id_comm .eq. id_rocn) then
    call cpl_init(n_comps, split_comm, id_rocn)
    call getenv('RUNDIR', env_tmp)
    call chdir(trim(env_tmp))
    call MPI_BARRIER(split_comm, ierr)
!    if(myrank1 .eq. 0) write(*,*) 'run TAI at: ', MPI_WTIME()
!    call tai_timcom_main(id_rocn, split_comm)
!    if(myrank1 .eq. 0) write(*,*) 'regional ocean model done'

  else
    write(*,'(a21,i4,a9,i4)') "color error, myrank =", myrank, ", color =", id_comm
  end if

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_FINALIZE(ierr) 
  stop
END PROGRAM driver

subroutine get_env_info(env_var, envinfo)
  implicit none

  character(len=*),intent(in) :: env_var
  integer, intent(out) :: envinfo

  character(len=4) :: env_tmp
  
  envinfo = -1
  call getenv(trim(env_var), env_tmp)
  read(env_tmp,'(i4)') envinfo

end subroutine get_env_info
