      subroutine mpmd_init(nsize,myrank,mycomm,who,istat)

      use mpmd_rank
      use mpi

      implicit none

      integer nsize,myrank,mycomm,who
      integer i,ierr,istat,iworld,igfs,irsm,NGFS,NRSM,NALL
      character*4 cgfs,crsm

      integer,allocatable ::  ranks_gfs(:) ,ranks_rsm(:)

! the whole group, (gfs + rsm)
      call MPI_INIT( istat )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myrank_all, ierr )
      istat=istat+ierr
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nsize_all, ierr )
      istat=istat+ierr

      call MPI_COMM_GROUP( MPI_COMM_WORLD, iworld, IERR )

      call getenv('GMPI',cgfs)
      read(cgfs,'(i4)')NGFS
      print *,'NGFS=',NGFS
      call getenv('RMPI',crsm)
      read(crsm,'(i4)')NRSM
      print *,'NRSM=',NRSM

      NALL=NGFS+NRSM
      if(NALL .ne. nsize_all)print*,'NGFS,NRSM,NALL=',NGFS,NRSM,NALL
      if(NALL .ne. nsize_all)stop'NGFS+NRSM .ne. nsize'

      allocate (ranks_gfs(NGFS),ranks_rsm(NRSM),stat=ierr)
      if (ierr/= 0) then
         write(6,*) 'mpmd_init : allocate fail '
      end if

      do i=1,NGFS
         ranks_gfs(i)=i-1
      enddo
      do i=1,NRSM
         ranks_rsm(i)=(i-1)+NGFS
      enddo

      call MPI_GROUP_excl( iworld, NRSM ,ranks_rsm, igfs, IERR )
      call MPI_GROUP_excl( iworld, NGFS ,ranks_gfs, irsm, IERR )

! create the sub_group(gfs)
      call MPI_COMM_create( MPI_COMM_WORLD,igfs,MPI_COMM_gfs,IERR )
! create the sub_group(rsm)
      call MPI_COMM_create( MPI_COMM_WORLD,irsm,MPI_COMM_rsm,IERR )

      if(myrank_all .le. NGFS-1)then

! gfs group goes here
         call MPI_COMM_RANK( MPI_COMM_gfs, myrank_gfs, ierr )
         call MPI_COMM_SIZE( MPI_COMM_gfs, nsize_gfs, ierr )
         myrank=myrank_gfs
         nsize=nsize_gfs
         mycomm=MPI_COMM_gfs
         who=NGFS
        print *, 'in mpmd init : myrank_gfs= ', myrank_gfs, ' nsize_gfs= ', nsize_gfs
         if(nsize.ne.NGFS)stop'mpmd_init error : NGFS inconsistent with mpiexec setting !'

      else
! rsm group goes here
         call MPI_COMM_RANK( MPI_COMM_rsm, myrank_rsm, ierr )
         call MPI_COMM_SIZE( MPI_COMM_rsm, nsize_rsm, ierr )
         myrank=myrank_rsm
         nsize=nsize_rsm
         mycomm=MPI_COMM_rsm
         who=0
        print *, 'in mpmd init : myrank_rsm= ', myrank_rsm, ' nsize_rsm= ', nsize_rsm
         if(nsize.ne.NRSM)stop'mpmd_init error : NRSM inconsistent with mpiexec setting !'

      endif

      root_gfs=0
      root_rsm=NGFS

      call MPI_GROUP_free( iworld, IERR )
      call MPI_GROUP_free( igfs ,  IERR )
      call MPI_GROUP_free( irsm ,  IERR )

! from now on, gfs and rsm groups can
!              communicate with whole group through MPI_COMM_WORLD, or
!              communicate with gfs   group through MPI_COMM_gfs,   or
!              communicate with rsm   group through MPI_COMM_rsm      

      deallocate (ranks_gfs,ranks_rsm)

      return
      end
