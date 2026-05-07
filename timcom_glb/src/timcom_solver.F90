module timcom_solver
  use timcom_const
  use timcom_comm
  implicit none

  integer, parameter :: solver_BiCGStab = 1,  & ! NTU 20250320 M. Hsieh
                            solver_PCSI = 2,  &
                         MaxLanczosStep = 20, & ! Max lanczos steps to get eigenvalues
!                        PCSICheckStart = 60, & ! start checking from given steps
!                         PCSICheckFreq = 10, & ! check convergence every freq steps
                         PCSICheckStart = 1,  &
                          PCSICheckFreq = 2,  &
                          maxIterations = 1000  ! max number of PCSI iterations

  real(r8), parameter :: &
!           PCSICriterion = 1.0d-12, & ! error criterion for PCSI
            LanczosCriterion = 0.1     ! Lanczos convergence criterion

  logical, parameter :: usePreconditioner = .true.

  real(r8) :: gamma_tmp,   & ! correction term to central coefficient
              PcsiMaxEig,  & ! Largest eigenvalue for PCSI method
              PcsiMinEig     ! smallest eigenvalue for PCSI method

contains
subroutine sipinc(al,ab,ac,ar,at,cl,cb,cc,cr,ct,alpha)
  use hyperlink, only: nx, ny
  implicit none

  real(r8), dimension(nx,ny), intent(in)  :: al, ab, ac, ar, at
  real(r8), intent(in) :: alpha
  real(r8), dimension(nx,ny), intent(out) :: cl, cb, cr, cc, ct
  real(r8), dimension(nx,ny) :: act
  integer :: i, j

  act = ac

if (usePreconditioner) then
  cl(1,1) = al(1,1)
  cb(1,1) = ab(1,1)
  cc(1,1) = 1.d0/ac(1,1)
  ct(1,1) = at(1,1)*cc(1,1)
  cr(1,1) = ar(1,1)*cc(1,1)

  do i=2,nx
    cl(i,1) = al(i,1)
    cb(i,1) = ab(i,1)/(1+alpha*cr(i-1,1))
    cc(i,1) = 1.d0/(act(i,1)+alpha*(cb(i,1)*cr(i-1,1))-cb(i,1)*ct(i-1,1))
    ct(i,1) = at(i,1)*cc(i,1)
    cr(i,1) = (ar(i,1)-alpha*cb(i,1)*ar(i-1,1))*cc(i,1)
  end do

  do j=2,ny
    cl(1,j) = al(1,j)/(1+alpha*ct(1,j-1))
    cb(1,j) = ab(1,j)
    cc(1,j) = 1.d0/(act(1,j)+alpha*(cl(1,j)*ct(1,j-1))-cl(1,j)*cr(1,j-1))
    ct(1,j) = (at(1,j)-alpha*cl(1,j)*ct(1,j-1))*cc(1,j)
    cr(1,j) = ar(1,j)*cc(1,j)
    do i=2,nx
    cl(i,j) = al(i,j)/(1+alpha*ct(i,j-1))
      cb(i,j) = ab(i,j)/(1+alpha*cr(i-1,j))
      cc(i,j) = 1.d0/(act(i,j)+alpha*(cl(i,j)*ct(i,j-1)+cb(i,j)*cr(i-1,j))-cl(i,j)*cr(i,j-1)-cb(i,j)*ct(i-1,j))
      ct(i,j) = (at(i,j)-alpha*cl(i,j)*ct(i,j-1))*cc(i,j)
      cr(i,j) = (ar(i,j)-alpha*cb(i,j)*cr(i-1,j))*cc(i,j)
    end do
  end do
else
  cl = 0.0
  cb = 0.0
  ct = 0.0
  cr = 0.0
  do j = 1, ny
  do i = 1, nx
     if ( abs(ac(i,j)) > 1.0d-18 ) then
        cc(i,j) = 1.0 / ac(i,j)
     else
        cc(i,j) = 0.0
     end if
  end do
  end do

end if
end subroutine sipinc

subroutine p_bicgstab_le(ab,al,ac,ar,at,b,x,cb,cl,cc,cr,ct,r,rh,p,v,s,t,ph,sh,iter,res2,res3)
  use hyperlink, only: nx, ny, m_comm_cart, r8type2d, nbid, ndim, myid, symm_np, threshold_bicg
  implicit none

  real(r8),dimension(nx,ny) :: &
    ab, al, ac, ar, at, &
    cb, cl, cc, cr, ct, b

  real(r8),dimension(0:nx+1,0:ny+1) :: &
    x, r, rh, p, v, s, t, ph, sh

  real(r8) :: res1, res2, res3, tp, tmp, rho, rhn, alpha, beta, w
  integer :: ierr, istat(MPI_STATUS_SIZE)
  integer :: iter, i, j, i0j0
  character(len=256) :: log_info

  i0j0 = (nx+2)*(ny+2)
  iter = 0
  rho = 1.d0
  alpha = 1.d0
  w = 1.d0
  
  call p_sqprod(ab,al,ac,ar,at,x,r)

  do j = 1, ny
    do i = 1, nx
      r(i,j) = b(i,j)-r(i,j)
    end do
  end do

  call DCOPY(i0j0,r,1,rh,1)
  call p_inprod(r,r,tp)

  res1 = dsqrt(tp)

  do iter = 1, 100
    call p_inprod(r,rh,rhn)
    beta = (rhn/rho)*(alpha/w)

    call DAXPY(i0j0, -w, v, 1, p, 1)
    call DSCAL(i0j0, beta, p, 1)
    call DAXPY(i0j0, 1.d0, r, 1, p, 1)

    call p_sipsqelim(cb,cl,cc,cr,ct,p,ph)
    call p_sqprod(ab,al,ac,ar,at,ph,v)
    call p_inprod(rh,v,tmp)
    alpha = rhn/tmp

    call DCOPY(i0j0, r, 1, s, 1)
    call DAXPY(i0j0,-alpha, v, 1, s, 1)

    call p_sipsqelim(cb,cl,cc,cr,ct,s,sh)
    call p_sqprod(ab,al,ac,ar,at,sh,t)
    call p_inprod(t,t,tmp)
    call p_inprod(t,s,w)
    w = w/tmp

    tmp = 0.d0

    call DAXPY(i0j0, alpha, ph, 1, x, 1)
    call DAXPY(i0j0,     w, sh, 1, x, 1)
    call DCOPY(i0j0,  s, 1, r, 1)
    call DAXPY(i0j0,    -w, t, 1, r, 1)

    call p_inprod(r,r,tp)
     
    rho = rhn

    res2 = dsqrt(tp)
    res3 = res2/res1
    if(res2 .lt. threshold_bicg .or. res3 .lt. 1.d-6 ) then
      call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, x, symm_np)
    !  if(myid .eq. rootid) then
    !    write(log_info,*) "bicg pressure solver convergence:", iter, res2, res3
    !    call comm_write_log_info(fid_log, log_info)
    !  end if
      return
    end if
    res1 = res2
  end do

  if(myid .eq. rootid) then
    write(log_info,*) "bicg pressure solver not convergence:", iter, res2, res3
    call comm_write_log_info(fid_log, log_info)
  end if
  call comm_finalize
end subroutine p_bicgstab_le

subroutine p_sqprod(ab,al,ac,ar,at,x,b)
  use hyperlink, only: m_comm_cart, nx, ny, r8type2d, nbid, ndim, symm_np
                     
  implicit none
  real(r8), dimension(nx,ny) :: ab, al, ac, ar, at
  real(r8), dimension(0:nx+1,0:ny+1) :: x, b
  integer :: ierr, istat(MPI_STATUS_SIZE)

  call MPI_BARRIER(m_comm_cart, ierr)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, x, symm_np)

  b(1:nx,1:ny) = ab*x(1:nx,  0:ny-1) &
               + al*x(0:nx-1,1:ny)   &
               + ac*x(1:nx,  1:ny)   &
               + ar*x(2:nx+1,1:ny)   &
               + at*x(1:nx,  2:ny+1)

#ifndef cpl
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, b, symm_np)
#endif

end subroutine p_sqprod

subroutine p_inprod(a1,a2,pd)
  use hyperlink, only: nx, ny, m_comm_cart
  implicit none

  real(r8), dimension(0:nx+1,0:ny+1) :: a1, a2
  real(r8), intent(out) :: pd
  real(r8) :: tmp
  integer :: ierr, i, j

  tmp = 0.d0

  do j = 1, ny
    do i = 1, nx
      tmp = tmp + a1(i,j)*a2(i,j)
    end do
  end do

  call MPI_ALLREDUCE(tmp, pd, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
end subroutine p_inprod

subroutine p_sipsqelim(cb,cl,cc,cr,ct,b,bh)
  use hyperlink, only: nx, ny
  implicit none

  real(r8), dimension(nx,ny), intent(in) :: &
    cb, cl, cc, cr, ct
  real(r8), intent(in)  :: b(0:nx+1,0:ny+1)
  real(r8), intent(out) :: bh(0:nx+1,0:ny+1) 

  integer :: i0j0, i, j

  i0j0 = (nx+2)*(ny+2)

  call DCOPY(i0j0, b, 1, bh, 1)

  do j = 1, ny-1
    do i = 1, nx-1
      bh(i,j)   = bh(i,j)*cc(i,j)
      bh(i+1,j) = bh(i+1,j) - cl(i,j)*bh(i,j)
      bh(i,j+1) = bh(i,j+1) - cb(i,j)*bh(i,j)
    end do
    bh(nx, j)   = bh(nx,j  )*cc(nx,j)
    bh(nx, j+1) = bh(nx,j+1) - cb(nx,j)*bh(nx,j)
  end do
  
  do i = 1, nx-1
    bh(i,  ny) = bh(i,ny)*cc(i,ny)
    bh(i+1,ny) = bh(i+1,ny) - cl(i,ny)*bh(i,ny)
  end do 
  bh(nx,ny) = bh(nx,ny)*cc(nx,ny)  

  do j = ny, 2, -1
    do i = nx, 2, -1
      bh(i-1,j) = bh(i-1,j) - cr(i,j)*bh(i,j)
      bh(i,j-1) = bh(i,j-1) - ct(i,j)*bh(i,j)
    end do
    bh(1,j-1) = bh(1,j-1) - ct(1,j)*bh(1,j)
  end do

  do i = nx, 2, -1
    bh(i-1,1) = bh(i-1,1) - cr(i,1)*bh(i,1)
  end do

end subroutine p_sipsqelim

! NTU 20250320 M. Hsieh

!subroutine PcsiLanczos(nx, ny, al, ab, ac, ar, at, &
subroutine PcsiLanczos(nx, ny, al, ab, ac, ar, at, &
                               cl, cb, cc, cr, ct, &
                               u, v, errorCode)

   use hyperlink, only: m_comm_cart, r8type2d, nbid, ndim, myid, symm_np, in
   implicit none

! Ported from POP in CESM

! !DESCRIPTION:
!  This routine implements the Lanczos based method to solve
!  the maximum and minimum eigenvalues of the coefficient matrix $A$.
!  Or $AM$ if a preconditioner $M^{-1}$ is provided.
!
!  References:
!     Hu, Y., Huang, X., Wang, X., Fu, H., Xu, S., Ruan, H., ... & Yang, G.
!     (2013). A
!        scalable barotropic mode solver for the parallel ocean program. In
!        Euro-Par 2013
!        Parallel Processing (pp. 739-750). Springer Berlin Heidelberg.
!     Paige, C. C. (1980). Accuracy and effectiveness of the Lanczos algorithm
!     for the
!        symmetric eigenproblem. Linear algebra and its applications, 34,
!        235-258.
!
! !REVISION HISTORY:
!  this routine implemented by Yong Hu, et al., Tsinghua University

 ! !INPUT PARAMETERS:

   integer, intent(in) :: nx, ny           ! block size

   real(r8), dimension(nx,ny), intent(in) :: al, ab, ac, ar, at, &
                                             cl, cb, cc, cr, ct

   real(r8), intent(inout) :: u, v         ! max and min eigenvalues

   integer, intent(out) :: errorCode       ! returned error code

 ! LOCAL VARIABLES

   real(r8), parameter :: &
      eps_mineig = 1.0d-8  ! convergence tolerance for the smallest eigenvalue

   real(r8), dimension(nx,ny) :: WORK  ! temporary vectors

   real(r8), dimension(0:nx+1,0:ny+1) ::  &
      R,                  & ! r
      S,                  & ! preconditioned r : Mr
      Q,                  & ! normalized r : r/||r||_A
      Q1,                 & ! previous step q
      P,                  & ! preconditioned q : Mq
      A0R,                & ! inverse of diagonal of A
      WORK1                 ! temporary vectors

   real(r8), dimension(MaxLanczosStep) :: &
      vcsa                  ! diagonal elements of tridiagonal matrix $T$
   real(r8), dimension(0:MaxLanczosStep) :: &
      vcsb                  ! off-diagonal elements of tridiagonal matrix $T$

   real(r8), dimension(:), allocatable :: &
      mcsa, mcsb                 ! temporary vectors saving

   integer :: i, j, m,      & ! local iteration counter
              info, istat,  & ! local flags
              ierr

   real(r8) :: csa,  &! inner product alpha = (p,r)
               csb,  &! inner product beta = (r,s)
               csc,  &! vector magnitude  ||r0||_A
               mineig, tmp, u0

   character(len=256) :: log_info

   if(myid_timcom .eq. rootid) then
      write(log_info,*) "(PcsiLanczos) nx,ny, initial u,v = ", nx,ny,u,v
      call comm_write_log_info(fid_log, log_info)
   end if

   errorCode = 0

   ! compute diagonal inverse for diagonal preconditioning

   if ( .not. usePreconditioner ) then
      do j = 1, ny
      do i = 1, nx
         if ( abs(ac(i,j)) > 1.0d-18 ) then
            A0R(i,j) = 1.0d0 / ac(i,j)
         else
            A0R(i,j) = 0.0d0
         end if
      end do
      end do
   end if

   ! set random vector R, initialize temporary vectors P, Q

    R(:,:) = 1.0
    Q(:,:) = 0.0
   Q1(:,:) = 0.0

   ! preconditioning s = M^{-1}r

   if (usePreconditioner) then
      call p_sipsqelim(cb,cl,cc,cr,ct,R,S)
   else
      do j = 1, ny
      do i = 1, nx
         S(i,j) = R(i,j)*A0R(i,j)
      end do
      end do
   endif
   do j = 1, ny
   do i = 1, nx
      WORK(i,j) = S(i,j)*R(i,j)*dble(in(i,j,1))
   end do
   end do
   ! (R,S)  $A$ and $M$  is negative defined for POP
   !                        positive definite(?) for TIMCOM so use a*0 instead

!  tmp = -sum(WORK)
   tmp = sum(WORK)

   call MPI_ALLREDUCE(tmp, csc, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

   ! q = R/||R||_{AM}  normalize R based on measurement of $AM^{-1}$

   if (csc .gt. 0.0) then
      tmp = 1.0/sqrt(csc)
      do j = 1, ny
      do i = 1, nx
         Q(i,j) = tmp*R(i,j)
      end do
      end do
   else      ! either random vector R == 0 or A is sigular
      if ( myid_timcom == rootid ) then
         write(log_info,*) &
            'PcsiLanczos: error1 estimating in lanczos: b == 0 !!! csc = ', csc
         call comm_write_log_info(fid_log, log_info)
      end if
      return
   endif
   call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, Q, symm_np)

      csb = 0.0
        u = 0.0
        v = 0.0
   mineig = 1.0
!-----------------------------------------------------------------------
!  lanczos iterations
!-----------------------------------------------------------------------
        u0 = -1.0e5
   vcsb(0) = 0.0
   lanczos_iter: do m = 1, MaxLanczosStep

   !  preconditioning q

      if (usePreconditioner) then
         call p_sipsqelim(cb,cl,cc,cr,ct,Q,P)
      else
         do j = 1, ny
         do i = 1, nx
            P(i,j) = Q(i,j)*A0R(i,j)
         end do
         end do
      endif
      call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, P, symm_np)

   ! r = Mp - beta_{j-1} q_{j-1}

      call p_sqprod(ab,al,ac,ar,at,P,WORK1)

      do j = 1, ny
      do i = 1, nx
         R(i,j) =  WORK1(i,j) - csb*Q1(i,j)
         WORK(i,j) = P(i,j)*R(i,j)*dble(in(i,j,1))
      end do
      end do

!     tmp = -sum(WORK)  ! (r,r)
      tmp = sum(WORK)
      call MPI_ALLREDUCE(tmp, csa, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

   ! r = r - alpha*q

      do j = 1, ny
      do i = 1, nx
         R(i,j) = R(i,j) - csa*Q(i,j)
      end do
      end do
      call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, R, symm_np)

   ! preconditioning R

      if (usePreconditioner) then
         call p_sipsqelim(cb,cl,cc,cr,ct,R,S)
      else
         do j = 1, ny
         do i = 1, nx
             S(i,j) = R(i,j)*A0R(i,j)
         end do
         end do
      endif
      do j = 1, ny
      do i = 1, nx
         WORK(i,j) = S(i,j)*R(i,j)*dble(in(i,j,1))
      end do
      end do

!     tmp = -sum(WORK)  ! (r,s)  $A$ and $M$  is negative defined for POP
      tmp = sum(WORK)
      call MPI_ALLREDUCE(tmp, csc, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

      if (csc .lt. 0.0) then
         if ( myid_timcom == rootid ) then
            write(log_info,*) "PcsiLanczos: negative csc = ", csc
            call comm_write_log_info(fid_log, log_info)
         end if
         return
      endif

      csb = sqrt(csc)

      vcsa(m) = csa
      vcsb(m) = csb

   ! Gershgorin circle theorem to estimate the largest eigenvalue

      tmp = vcsa(m) + vcsb(m) + vcsb(m-1)
      u = 0.5 * ( tmp + u0 + abs(tmp-u0) )
      u0 = u

   ! normalize r  : q = r/||r||_A

      if ( abs(csb) > 1.0d-18 ) then
         tmp = 1.0/csb
         do j = 1, ny
         do i = 1, nx
            Q1(i,j) = Q(i,j)
             Q(i,j) = tmp*R(i,j)
         end do
         end do
      else
         if ( myid_timcom == rootid ) then
            write(log_info,'(a)') &
               'PcsiLanczos: error2 estimating in lanczos: b == 0 !!!'
            call comm_write_log_info(fid_log, log_info)
         end if
         return
      endif

!-----------------------------------------------------------------------
!  compute eigenvalues of the tridiagonal matrix T every 10 steps
!-----------------------------------------------------------------------

      if ( ( mod(m,10) == 0 ) .or. ( m == MaxLanczosStep ) ) then

         if ( myid_timcom == rootid ) then
            allocate(mcsa(m), mcsb(m), stat = istat)
            if ( istat > 0 ) then
               write(log_info,'(a)') &
                  'PcsiLanczos: error allocating Lanczos Tridiagonal'
               call comm_write_log_info(fid_log, log_info)
               return
            endif

            do i = 1, m-1
                mcsa(i)   = vcsa(i)
                mcsb(i+1) = vcsb(i)
            end do
            mcsa(m) = vcsa(m)
            mcsb(1) = 0.0

            info = 0
            ! compute smallest eigenvalue
            call ratQR(m, eps_mineig, mcsa,mcsb,v,info)
            if ( info /= 0 ) then
               write(log_info,'(a)') &
                  'PcsiLanczos: error estimating smallest eigenvalue !!!'
               call comm_write_log_info(fid_log, log_info)
               return
            endif

            write(log_info,'(a18,i3,a8,2e15.7)') &
               "Lanczos steps ", m, " eigs:", v, u

            deallocate(mcsa, mcsb, stat = istat)
            if ( istat > 0 ) then
               write(log_info,'(a)') &
                  'PcsiLanczos: error deallocating Lanczos Tridiagonal'
               call comm_write_log_info(fid_log, log_info)
               return
            endif
         endif
         call MPI_BCAST(v, 1, MPI_REAL8, 0, m_comm_cart, ierr)

         ! check convergence of the smallest eigenvalue
         if ( abs(1.0-v/mineig ) <  LanczosCriterion ) then
            exit lanczos_iter
         endif
         mineig = v  ! save previous value

      endif

   enddo lanczos_iter

end subroutine PcsiLanczos

subroutine ratQR(n, eps1, d, e, mineig, ierr)

! !DESCRIPTION:
!  This subroutine finds the algebraically smallest
!  eigenvalue of a positive definite symmetric tridiagonal matrix by the
!  rational QR method with Newton corrections.
!  Follow EISPACK lib subroutine RATQR
!  Ref: https://www.netlib.org/eispack/ratqr.f
!
! ! INPUT VARIABLES
   integer,  intent(in) :: n       ! matrix order
   real(r8), intent(in) :: d(n), & ! diagonal elements
                           e(n), & ! off-diagonal elements, e(1) is arbitrary
                           eps1    ! convergence tolerance

!  ! OUTPUT VARIABLES
   integer,  intent(inout) :: ierr    ! status flag
   real(r8), intent(inout) :: mineig  ! smallest eigevalue

!  ! LOCAL VARIABLES
   real(r8), dimension(n) :: e2, bd, w
   real(r8) :: f, ep, delta, ferr
   integer i, ii ! local counters
   real(r8) p, q, qp, r, s, tot

   ierr = 0
   w(1:n) = d(1:n)
   ferr = 0.0
   s = 0.0
   tot = w(1)
   q = 0.0

   do i = 1, n
      p = q
      bd(i) = e(i)*e(i)
      q = 0.0
      if ( i /= n ) then
        q = abs( e(i+1) )
      endif
      tot = min( w(i) - p - q, tot )
   end do
   bd(1) = 0.0

   if ( tot < 0.0 ) then
     tot = 0.0
   else
     w(1:n) = w(1:n) - tot
   endif

!  QR transformation.

   do while ( .true. )

      tot = tot + s
      delta = w(n) - s
      i = n
      if ( delta <= eps1 ) exit

      f = bd(n) / delta
      qp = delta + f
      p = 1.0

      do ii = 1, n - 1
         i = n - ii
         q = w(i) - s - f
         r = q / qp
         p = p * r + 1.0
         ep = f * r
         w(i+1) = qp + ep
         delta = q - ep

         ! check convergence
         if ( delta <= eps1 ) exit

         f = bd(i) / q
         qp = delta + f
         bd(i+1) = qp * ep
      end do

      ! check convergence
      if ( delta <= eps1 ) exit

      w(1) = qp
      s = qp / p

      !  Set error: irregular end of iteration.
      if ( tot + s <= tot ) then
         ierr = 1
         return
      endif

   end do

   w(1) = tot
   ferr = ferr + abs(delta)
   bd(1) = ferr
   mineig = w(1)

end subroutine ratQR

subroutine PCSI(ab,al,ac,ar,at,cb,cl,cc,cr,ct,B,X,numIterations,rmsResidual)

  use hyperlink, only: nx, ny, in, &
                       m_comm_cart, r8type2d, nbid, ndim, myid, symm_np, &
                       total_area, PCSICriterion  ! NTU 20250324 M. Hsieh
  implicit none

! !DESCRIPTION:
!  This routine implements the Preconditioned Classical Stiefel Iteration
!  (PCSI)  solver for solving the linear system $Ax=b$.
!  It uses the two extreme eigenvalues of $A$ instead of the norm of
!  residual. Thus, PCSI eliminates global reductions in each iteration.
!  The eigenvalues are estimated by routine PcsiLanczos.
!  PCSI supports all kinds of preconditioners supported by PCG/ChronGear,
!  including diagonal and EVP preconditioning.
!
!  References:
!     Stiefel, E. L. (1958). Kernel polynomial in linear algebra and their
!        numerical applications, in: Further contributions to the
!        determination of eigenvalues. NBS Applied Math. Ser., 49, 1-22.
!     Hu, Y., Huang, X., Wang, X., Fu, H., Xu, S., Ruan, H., Xue, W. and Yang,
!     G. (2013).
!        A scalable barotropic mode solver for the parallel ocean program.
!        In Euro-Par 2013  Parallel Processing (pp. 739-750) Springer Berlin
!        Heidelberg.
!
! !REVISION HISTORY:
!  this routine implemented by Yong Hu et al., Tsinghua University

! !INPUT PARAMETERS:

   real(r8), dimension(nx,ny), intent(in) :: &
      ab, al, ac, ar, at, &
      cb, cl, cc, cr, ct, &
      B    ! right hand side of linear system

   real(r8), dimension(0:nx+1,0:ny+1), intent(inout) :: &
      X    ! on input,  an initial guess for the solution
           ! on output, solution of the linear system

! OUTPUT PARAMETERS:

   integer, intent(out) :: &
      numIterations  ! accumulated no of iterations (diagnostic)

   real(r8), intent(out) :: rmsResidual  ! residual (diagnostic)

! local variables

   integer  :: i,j,m  ! local iteration counter

   real(r8) :: csalpha,csbeta,csy,csomega,one_csy,cc1, &
               tmp,rr,       &  ! scalar inner product results
               residualNorm     ! residual normalization

   real(r8), dimension(nx,ny) :: WORK0

   real(r8), dimension(0:nx+1,0:ny+1) :: &
      R,   &  ! residual (b-Ax)
      S,   &  ! conjugate direction vector
      Q,   &  ! various cg intermediate results
      A0R     ! not used w/o preconditioning

   integer :: ierr, istat(MPI_STATUS_SIZE)
   character(len=256) :: log_info

   if ( .not. usePreconditioner ) then
      !--- diagonal preconditioner if preconditioner not specified
      do j = 1, ny
      do i = 1, nx
         if ( abs(ac(i,j)) > 1.0d-18 ) then
            A0R(i,j) = 1.0d0 / ac(i,j)
         else
            A0R(i,j) = 0.0d0
         end if
      end do
      end do
   end if
!-----------------------------------------------------------------------
!  step 1 : compute iteration parameters by eigenvalues
!  $\alpha =\frac{2}{\mu -\nu}$, $ \beta = \frac{\mu +\nu}{\mu -\nu}$,
!  $\gamma = \frac{\beta}{\alpha}$, $\omega_0 =\frac{ 2}{\gamma}$
!-----------------------------------------------------------------------

   csalpha = 2.0d0/(PcsiMaxEig-PcsiMinEig)
   csbeta = (PcsiMaxEig+PcsiMinEig)/(PcsiMaxEig-PcsiMinEig)
   csy = csbeta/csalpha
   csomega = 2.0d0/csy
   one_csy = 1.0d0/csy
   residualNorm = 1.0d0/total_area

!  if ( myid_timcom == rootid ) then
!     write(log_info,*) "(PCSI) csy,resNorm = ", csy, residualNorm
!     call comm_write_log_info(fid_log, log_info)
!  end if
!-----------------------------------------------------------------------
! step 2 : compute initial residual and initialize X
! $\textbf{r}_0 = \textbf{b}-\textbf{B}\textbf{x}_0$;
! $\textbf{x}_1 =\textbf{x}_0 -\gamma^{-1}\textbf{M}^{-1}\textbf{r}_0$;
! $\textbf{r}_1 =\textbf{b} -\textbf{B}\textbf{x}_1$;
!-----------------------------------------------------------------------

   call p_sqprod(ab,al,ac,ar,at,X,S)
   R(1:nx,1:ny) = B(1:nx,1:ny) - S(1:nx,1:ny)
!  do j = 1, ny
!  do i = 1, nx
!     R(i,j) = B(i,j) - S(i,j)
!  end do
!  end do

   if (usePreconditioner) then
      call p_sipsqelim(cb,cl,cc,cr,ct,R,S)
      R(1:nx,1:ny) = S(1:nx,1:ny)
   else
      R(1:nx,1:ny) = R(1:nx,1:ny) * A0R(1:nx,1:ny)
   end if

   Q(1:nx,1:ny) = one_csy * R(1:nx,1:ny)
   call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, Q, symm_np)

   X(1:nx,1:ny) = X(1:nx,1:ny) + Q(1:nx,1:ny)

   call p_sqprod(ab,al,ac,ar,at,X,S)
   R(1:nx,1:ny) = B(1:nx,1:ny) - S(1:nx,1:ny)
   call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, R, symm_np)

   numIterations = maxIterations

!  if ( myid_timcom == rootid ) then
!     write(log_info,*) "(PCSI) b4 iter, B, R = ", B(2,2), R(2,2)
!     call comm_write_log_info(fid_log, log_info)
!  end if

   iterationLoop: do m = 1, maxIterations

!-----------------------------------------------------------------------
! step 3 : update iteration parameter
! $\omega_k = 1/(\gamma - \frac{1}{4\alpha^2}\omega_{k-1})$
!-----------------------------------------------------------------------
      csomega = 1.0d0/(csy-csomega/(4.0d0*csalpha*csalpha))
!  if ( myid_timcom == rootid ) then
!     write(log_info,*) "(PCSI) ", m, "th iter, csomega = ", csomega
!     call comm_write_log_info(fid_log, log_info)
!  end if
!-----------------------------------------------------------------------
! step 4 : preconditioning
! $\textbf{r}'_{k-1} =\textbf{M}^{-1}\textbf{r}_{k-1}$
!-----------------------------------------------------------------------
      if (usePreconditioner) then
         call p_sipsqelim(cb,cl,cc,cr,ct,R,S)
         R(1:nx,1:ny) = S(1:nx,1:ny)
      else
         R(1:nx,1:ny) = R(1:nx,1:ny) * A0R(1:nx,1:ny)
      end if
      call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, R, symm_np)
!-----------------------------------------------------------------------
! step 5 : compute X increment and update X and R
! $\Delta \textbf{x}_{k} =\omega_k\textbf{r}'_{k-1}+(\gamma \omega_k-1)\Delta
! \textbf{x}_{k-1}$;
! $\textbf{x}_{k} =\textbf{x}_{k-1}+\Delta \textbf{x}_{k-1}$;
! $\textbf{r}_{k} =b- \textbf{B}\textbf{x}_{k}$;
!-----------------------------------------------------------------------
      cc1 = csy*csomega-1.0d0
      Q(1:nx,1:ny) = csomega*R(1:nx,1:ny) + cc1*Q(1:nx,1:ny)

      X(1:nx,1:ny) = X(1:nx,1:ny) + Q(1:nx,1:ny)

      call p_sqprod(ab,al,ac,ar,at,X,S)
      R(1:nx,1:ny) = B(1:nx,1:ny) - S(1:nx,1:ny)

      if ( ( m >= PCSICheckStart ) .and. &
           ( mod(m,PCSICheckFreq) == 0 ) ) then
         WORK0(1:nx,1:ny) = R(1:nx,1:ny)*R(1:nx,1:ny)*dble(in(1:nx,1:ny,1))      ! (r,r)
         tmp = sum(WORK0)
         call MPI_ALLREDUCE(tmp, rr, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

!        if ( myid_timcom == rootid ) then
!           write(log_info,*) "(PCSI) ", m, "th iter, tmp, sum(tmp) = ", tmp, rr
!           call comm_write_log_info(fid_log, log_info)
!        end if

         if ( rr < PCSICriterion ) then
            numIterations = m
            exit iterationLoop
         end if
      end if
   enddo iterationLoop

   rmsResidual = sqrt(rr*residualNorm)

end subroutine PCSI

end module timcom_solver
