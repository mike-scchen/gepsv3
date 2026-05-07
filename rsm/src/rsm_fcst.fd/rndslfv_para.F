      subroutine rpara_we2ns(a,b,lev)
!yj: RSM input lev = partial lev = levrp
!yj      subroutine rpara_we2ns(a,b,lev,global_lats_a,latg)
!
! mpi transport from full dimension of west-east to full dimension of 
! north-south with latitude shuffl. created by hann-ming henry juang
! 
! program log
! 2011 02 20 : henry juang, created for ndsl advection
!
!
!yj      use gfs_dyn_mpi_def
!yj      use gfs_dyn_layout1
      use mod_rndsl
      implicit none
!yj#include <paramodel.h>
#include <define.h>
#include <npesi.h>
#include <comrmpi.h>
!
!yj      integer lev,latg
      integer lev
      real a(lonfull,lev,latpart)
      real b(latfull,lev,lonpart)
!yj      integer global_lats_a(latg)
!
      real (kind=kind_mpi) works(1,lev,lonlenmax*latlenmax,nrow) 
      real (kind=kind_mpi) workr(1,lev,lonlenmax*latlenmax,nrow) 
!yj      real (kind=kind_mpi) works(1,lev,lonlenmax*latlenmax,nrow) 
!yj      real (kind=kind_mpi) workr(1,lev,lonlenmax*latlenmax,nrow) 
      integer lensend(nrow),lenrecv(nrow)
      integer locsend(nrow),locrecv(nrow)
      integer i,j,k,n,mn,jj,ierr
!yj      integer i,j,k,n,mn,jj,lat1,lat2,ierr
!
!yj      print *,'mype, enter rndslfv_we2ns ',mype

#ifdef OPENMP
!$omp parallel do private(n,mn,i,j,k) &
!$omp schedule(dynamic)
#endif
      do n=1,nrow 
        mn=0
        do j=1,mylatlen
          do i=lonstrndsl(n),lonstrndsl(n)+lonlenndsl(n)-1
            mn = mn + 1
            do k=1,lev
              works(1,k,mn,n)=a(i        ,k,j)
            enddo
          enddo
        enddo
        lensend(n) = mn * 1 * lev 
        locsend(n) = (n-1)*lonlenmax*latlenmax*1*lev
        lenrecv(n)=latlenndsl(n)*mylonlen*1*lev
        locrecv(n)=locsend(n)
!yj        if (mype.eq.2) then
!yj        print*,'mype,n,lensend,lenrecv=',mype,n,lensend(n),lenrecv(n)
!yj        print*,'mype,n,latlenndsl(n)=',mype,n,latlenndsl(n)
!yj        print*,'mype,n,mylonlen=',mype,n,mylonlen
!yj        print*,'mype,n,locsend(n),locrecv(n)=',mype,n,locsend(n),locrecv(n)
!yj        print*,'mype,n,mn,lev=',mype,n,mn,lev
!yj        endif
      enddo
#ifdef OPENMP
!$omp end parallel do
#endif
!
!yj      print*,'before mpi_barrier'
      call mpi_barrier (comm_column,ierr)  
!yj      print*,'after mpi_barrier'
!yj      call mpi_barrier (mc_comp,ierr)  !!yj mc_comp->com_row?
      
!     call mpi_alltoallv(works,lensend,locsend,mpi_r_mpi_r,
!    &                   workr,lenrecv,locrecv,mpi_r_mpi_r,
!    &                   mc_comp,ierr)
      call mpi_alltoallv(works,lensend,locsend,MPIREAL,   &
     &                   workr,lenrecv,locrecv,MPIREAL,   &
     &                   comm_column,ierr)
!yj        print*,'mype, myrow,finish mpi_alltoallv,ierr=',  &
!yj     &       mype,myrow,ierr
!yj     &                   mc_comp,ierr)
     
!
#ifdef OPENMP
!$omp parallel do private(n,mn,i,j,jj,k) &
!$omp schedule(dynamic)
#endif
!yj      do n=1,nodes
      do n=1,nrow
        mn=0
        do j=1,latlenndsl(n)
          jj = latstrndsl(n) + j - 1
!yj          lat1=global_lats_a(jj) !!yj no need
!yj          lat2=latfull+1-lat1    !!yj no need
          do i=1,mylonlen
            mn = mn + 1
            do k=1,lev
              b(jj,k,i) = workr(1,k,mn,n)  !!yj lat1->jj
!yj              b(lat1,k,i) = workr(1,k,mn,n)  !!yj lat1->jj
            enddo
          enddo
        enddo
      enddo
#ifdef OPENMP
!$omp end parallel do
#endif
!
!yj      print *,'mype, end of ndslfv_we2ns ',mype

      return
      end

! ======================================================================
      subroutine rpara_ns2we(a,b,lev) !yj keepa,b,lev
!yj      subroutine para_ns2we(a,b,lev,global_lats_a,latg) !yj keepa,b,lev
!
! mpi transport from full dimension of west-east to full dimension of 
! north-south with latitude shuffl.
!
!yj      use gfs_dyn_mpi_def
!yj      use gfs_dyn_layout1
      use mod_rndsl
      implicit none
!yj#include <paramodel.h>
#include <define.h>
#include <npesi.h>
#include <comrmpi.h>
!
!yj      integer lev,latg
      integer lev
      real a(latfull,lev,lonpart)
      real b(lonfull,lev,latpart)
!yj      integer global_lats_a(latg)
!
!     real (kind=kind_mpi_r) works(2*lev*mylonlen*lathalf)
!     real (kind=kind_mpi_r) workr(2*lev*mylatlen*lonhalf)
      real (kind=kind_mpi) works(1,lev,lonlenmax*latlenmax,nrow)
      real (kind=kind_mpi) workr(1,lev,lonlenmax*latlenmax,nrow)
!yj      real (kind=kind_mpi) works(2,lev,lonlenmax*latlenmax,nodes)
!yj      real (kind=kind_mpi) workr(2,lev,lonlenmax*latlenmax,nodes)
      integer lensend(nrow),lenrecv(nrow)
      integer locsend(nrow),locrecv(nrow)
!yj      integer lensend(nodes),lenrecv(nodes)
!yj      integer locsend(nodes),locrecv(nodes)
      integer i,j,k,n,mn,jj,ierr
!yj      integer i,j,k,n,mn,jj,lat1,lat2,ierr
!
!yj      print *,'mype, enter rndslfv_ns2we ',mype

#ifdef OPENMP
!$omp parallel do private(n,mn,i,j,jj,k) &
!$omp schedule(dynamic)
#endif
      do n=1,nrow
!yj      do n=1,nodes
        mn=0
        do j=1,latlenndsl(n)
          jj=latstrndsl(n) + j -1
!yj          lat1=global_lats_a(jj)
!yj          lat2=latfull+1-lat1
          do i=1,mylonlen
            mn = mn + 1
            do k=1,lev
              works(1,k,mn,n) = a(jj,k,i)
!yj              works(1,k,mn,n) = a(lat1,k,i)
!yj              works(2,k,mn,n) = a(lat2,k,i)
            enddo
          enddo
        enddo
        lensend(n)=mn*1*lev
        locsend(n)=(n-1)*lonlenmax*latlenmax*1*lev
        lenrecv(n)=mylatlen*lonlenndsl(n)*1*lev
!yj        lensend(n)=mn*2*lev
!yj        locsend(n)=(n-1)*lonlenmax*latlenmax*2*lev
!yj        lenrecv(n)=mylatlen*lonlenndsl(n)*2*lev
        locrecv(n)=locsend(n)
      enddo
#ifdef OPENMP
!$omp end parallel do
#endif
!
!yj      print*,'before mpi_barrier'
      call mpi_barrier (comm_column,ierr)
!yj      print*,'after mpi_barrier'
!yj      call mpi_barrier (mc_comp,ierr)
!     call mpi_alltoallv(works,lensend,locsend,mpi_r_mpi_r,
!    &                   workr,lenrecv,locrecv,mpi_r_mpi_r,
!    &                   mc_comp,ierr)
      call mpi_alltoallv(works,lensend,locsend,MPIREAL,   &
     &                   workr,lenrecv,locrecv,MPIREAL,   &
     &                   comm_column,ierr)
!yj     &                   mc_comp,ierr)
!
#ifdef OPENMP
!$omp parallel do private(n,mn,i,j,k) &
!$omp schedule(dynamic)
#endif
!yj      do n=1,nodes 
      do n=1,nrow 
        mn=0
        do j=1,mylatlen
          do i=lonstrndsl(n),lonstrndsl(n)+lonlenndsl(n)-1
            mn = mn + 1
            do k=1,lev
              b(i        ,k,j) = workr(1,k,mn,n)
!yj              b(i+lonhalf,k,j) = workr(2,k,mn,n)
            enddo
          enddo
        enddo
      enddo
#ifdef OPENMP
!$omp end parallel do
#endif
!
!yj      print *,'mype, end of rndslfv_ns2we ',mype
!
      return
      end
