      subroutine rfftmlt(a,work,trigs,ifax,inc,jump,n,m,isign)

!CWB2018 for gfs by CH Lee, using fftw-3.3.7

      implicit none

      include "fftw3.f"

      real*8 a(jump,m),work(*),trigs(*),scale
      integer ifax(*),inc,jump,n,m,isign,i,j,nn

      integer*8 plan_c2r,plan_r2c
      integer   fftwm,fftwn,fftwtype

!      common/cwbfftw/plan_c2r,plan_r2c, &
!                fftwm,fftwn,fftwtype

!     logical,save :: first_c2r,first_r2c
!     data    first_c2r,first_r2c/.true. , .true./

!$omp threadprivate(/cwbfftw/)

      nn=jump/2
      if(mod(jump,2).eq.0) then  !jump is even 
         if(isign.eq.1) then
!           if(isign.ne.fftwtype.or.m.ne.fftwm.or.n.ne.fftwn) then
!              if(.not.first_c2r) call dfftw_destroy_plan(plan_c2r)
!              first_c2r=.false.

! still got memory leak problem
!              if(.not.first_c2r) then
!                 call dfftw_destroy_plan(plan_c2r)
!                 call dfftw_cleanup
!                 first_c2r=.false.
!              endif

               call dfftw_plan_many_dft_c2r(plan_c2r,1,n,m,a,n,inc,nn
     &              ,a,n,inc,jump,FFTW_ESTIMATE)

!              fftwn=n
!              fftwm=m
!              fftwtype=isign
!           endif

            call dfftw_execute(plan_c2r)

         else
            scale=1.0/dfloat(n)

!           if(isign.ne.fftwtype.or.m.ne.fftwm.or.n.ne.fftwn) then
!              if(.not.first_r2c) call dfftw_destroy_plan(plan_r2c)
!              first_r2c=.false.

! still got memory leak problem
!              if(.not.first_r2c) then
!                 call dfftw_destroy_plan(plan_r2c)
!                 call dfftw_cleanup
!                 first_r2c=.false.
!              endif

               call dfftw_plan_many_dft_r2c(plan_r2c,1,n,m,a,n,inc,jump
     &              ,a,n,inc,nn,FFTW_ESTIMATE)

!              fftwn=n
!              fftwm=m
!              fftwtype=isign
!           endif

            call dfftw_execute(plan_r2c)

            do j=1,m
            do i = 1, n+2
               a(i,j)=a(i,j)*scale
            enddo
            enddo
         endif
      else    !jump is odd
         print *,'CWB obsolete, jump=',jump
      endif

      return
      end

!-------------------------------------------------
      subroutine fftfax (n,ifax,trigs)

      real*8  trigs(*)
      integer ifax(*),n

      character*80 wisdom_file

      call cwb_fftw_init_threads( )

      call cwb_fftw_make_planner_thread_safe()

      call cwb_fftw_plan_with_nthreads( )

!     wisdom_file="/users/xa09/fftw/wisdom.gfsT512L60"//char(0)
!     call cwb_fftw_import_wisdom_from_filename(wisdom_file)

      return
      end

!-------------------------------------------------
      subroutine export_wisdom

      character*80 wisdom_file

      wisdom_file="/users/xa09/fftw/wisdom.gfsT512L60"//char(0)
      call cwb_fftw_export_wisdom_to_filename(wisdom_file)

      return
      end
