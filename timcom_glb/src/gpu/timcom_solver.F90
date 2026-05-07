!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define CUDACHECK(ierr) call ocn_cuda_check_helper(ierr, __FILE__, __LINE__)
#define NCCLCHECK(ierr) call ocn_nccl_check_helper(ierr, __FILE__, __LINE__)

module timcom_solver_gpu
  use timcom_const
  use timcom_comm
  use timcom_comm_gpu
  use cudafor
  use timcom_solver, only: usePreconditioner, PcsiMaxEig, PcsiMinEig, &
                           maxIterations, PCSICheckFreq, PCSICheckStart, &
                           p_sqprod
  use hyperlink, only: pcsi_cuda_graph
  implicit none

  real(r8) :: rhn, rho, alpha, beta, w, tmp, tp, res1, res2, res3, &
              csalpha, csbeta, csy, csomega, one_csy, rr
  real(r8), allocatable, dimension(:,:,:) :: &
    r      ! residual (b-Ax)
  real(r8), allocatable, dimension(:,:) :: &
    s,   &  ! conjugate direction vector
    q,   &  ! various cg intermediate results
    a0r     ! not used w/o preconditioning
  type(c_ptr) :: h_ptr_rr
  type(c_devptr) :: d_ptr_rr
  real(r8), pointer :: zc_rr
  real(r8), device, pointer :: d_zc_rr
  type(cudaEvent) :: event
  logical :: bicgstab_cg_created = .false., pcsi_init_cg_created = .false., &
             pcsi_stage1_cg_created(2) = .false., pcsi_stage2_cg_created(2) = .false., &
             pcsi_cuda_init = .false.
  type(cudaGraph) :: bicgstab_graph, pcsi_init_graph, pcsi_stage1_graph(2), pcsi_stage2_graph(2)
  type(cudaGraphExec) :: bicgstab_graph_exec, pcsi_init_graph_exec, pcsi_stage1_graph_exec(2), pcsi_stage2_graph_exec(2)

contains

subroutine p_bicgstab_le_gpu(ab,al,ac,ar,at,b,x,cb,cl,cc,cr,ct,r,rh,p,v,s,t,ph,sh,iter, async_id)
  ! Present on device
  !  Array: ab, al, ac, ar, at, b, x, r, rh, p, v, s, t, ph, sh
  use hyperlink, only: nx, ny, m_comm_cart, r8type2d, nbid, ndim, myid, symm_np, threshold_bicg
  use nccl
  use openacc
  use cudafor
  implicit none

  real(r8),dimension(nx,ny) :: &
    ab, al, ac, ar, at, &
    cb, cl, cc, cr, ct, b

  real(r8),dimension(0:nx+1,0:ny+1) :: &
    x, r, rh, p, v, s, t, ph, sh

  integer :: iter, i, j, n, l
  character(len=256) :: log_info
  integer :: async_id
  logical :: converged
  real(r8) :: r_t, sum_t, sum2_t
  integer(kind=cuda_stream_kind) :: stream
  integer, parameter :: tiling_size = 128

  stream = acc_get_cuda_stream(async_id)

  converged = .false.

  !$acc kernels present(rho, alpha, w, tp) async(async_id)
  rho = 1.d0
  alpha = 1.d0
  w = 1.d0
  tp = 0.d0
  rhn = 0.d0
  !$acc end kernels

  ! sqprod(A, x, r); r = b - r; rh = r; inprod(r, r, tp)
  call ghost_cell_exch_gpu(m_comm_cart, r8type2d, nbid, ndim, x, symm_np, async_id)

  !$acc parallel loop gang vector_length(tiling_size) async(async_id) private(sum_t)
  do n = 0, nx*ny - 1, tiling_size
    sum_t = 0.d0
    !$acc loop vector reduction(+:sum_t) private(r_t)
    do l = 0, tiling_size - 1
      i = mod(n + l, nx) + 1
      j = (n + l) / nx + 1
      r_t = ab(i, j)*x(i, j - 1) &
        + al(i, j)*x(i - 1, j) &
        + ac(i, j)*x(i, j) &
        + ar(i, j)*x(i + 1, j) &
        + at(i, j)*x(i, j + 1)
      r_t = b(i, j) - r_t
      r(i, j) = r_t
      rh(i, j) = r_t
      sum_t = sum_t + r_t * r_t
    end do
    !$acc atomic
    tp = tp + sum_t
  end do

  !$acc host_data use_device(tp)
  NCCLCHECK(ncclAllReduce(tp, tp, 1, ncclFloat64, ncclSum, nccl_comm, stream))
  !$acc end host_data

  !$acc kernels present(res2, tp) async(async_id)
  res2 = dsqrt(tp)
  !$acc end kernels

  do iter = 1, 100
    
    ! inprod(r, rh, rhn)
    !$acc parallel loop gang vector_length(tiling_size) async(async_id) private(sum_t)
    do n = 0, nx*ny - 1, tiling_size
      sum_t = 0.d0
      !$acc loop vector reduction(+:sum_t) private(r_t)
      do l = 0, tiling_size - 1
        i = mod(n + l, nx) + 1
        j = (n + l) / nx + 1
        sum_t = sum_t + r(i,j) * rh(i,j)
      end do
      !$acc atomic
      rhn = rhn + sum_t
    end do

    !$acc host_data use_device(rhn)
    NCCLCHECK(ncclAllReduce(rhn, rhn, 1, ncclFloat64, ncclSum, nccl_comm, stream))
    !$acc end host_data

    !$acc kernels present(beta, rhn, rho, alpha, w, tmp, res1, res2) async(async_id)
    beta = (rhn/rho)*(alpha/w)
    tmp = 0.d0
    res1 = res2
    !$acc end kernels

    ! p = r + beta * (-w * v + p); precondition(p, ph); sqprod(a, ph, v); inprod(rh, v, tmp)
    !$acc parallel loop collapse(2) async(async_id) present(beta, w) private(r_t)
    do j = 1, ny
      do i = 1, nx
        r_t = r(i, j) + beta*(-w*v(i, j) + p(i, j))
        p(i, j) = r_t
        ph(i, j) = r_t / ac(i, j)
      end do
    end do

    call ghost_cell_exch_gpu(m_comm_cart, r8type2d, nbid, ndim, ph, symm_np, async_id)

    !$acc parallel loop gang vector_length(tiling_size) async(async_id) private(sum_t)
    do n = 0, nx*ny - 1, tiling_size
      sum_t = 0.d0
      !$acc loop vector reduction(+:sum_t) private(r_t)
      do l = 0, tiling_size - 1
        i = mod(n + l, nx) + 1
        j = (n + l) / nx + 1
        r_t = ab(i, j)*ph(i, j - 1) &
          + al(i, j)*ph(i - 1, j) &
          + ac(i, j)*ph(i, j) &
          + ar(i, j)*ph(i + 1, j) &
          + at(i, j)*ph(i, j + 1)
        v(i, j) = r_t
        sum_t = sum_t + rh(i, j) * r_t
      end do
      !$acc atomic
      tmp = tmp + sum_t
    end do

    !$acc host_data use_device(tmp)
    NCCLCHECK(ncclAllReduce(tmp, tmp, 1, ncclFloat64, ncclSum, nccl_comm, stream))
    !$acc end host_data

    !$acc kernels async(async_id) present(alpha, rhn, tmp, w)
    alpha = rhn/tmp
    tmp = 0.d0
    w = 0.d0
    !$acc end kernels

    ! s = r; precondition(sh, s); sqprod(a, sh, t); inprod(t, t, tmp); inprod(t, s, w)
    !$acc parallel loop collapse(2) async(async_id) present(alpha) private(r_t)
    do j = 1, ny
      do i = 1, nx
        r_t = -alpha*v(i, j) + r(i, j)
        s(i, j) = r_t
        sh(i, j) = r_t / ac(i, j)
      end do
    end do

    call ghost_cell_exch_gpu(m_comm_cart, r8type2d, nbid, ndim, sh, symm_np, async_id)

    !$acc parallel loop gang vector_length(tiling_size) async(async_id) private(sum_t, sum2_t)
    do n = 0, nx*ny - 1, tiling_size
      sum_t = 0.d0
      sum2_t = 0.d0
      !$acc loop vector reduction(+:sum_t) reduction(+:sum2_t) private(r_t)
      do l = 0, tiling_size - 1
        i = mod(n + l, nx) + 1
        j = (n + l) / nx + 1
        r_t = ab(i, j)*sh(i, j - 1) &
          + al(i, j)*sh(i - 1, j) &
          + ac(i, j)*sh(i, j) &
          + ar(i, j)*sh(i + 1, j) &
          + at(i, j)*sh(i, j + 1)
        t(i, j) = r_t
        sum_t = sum_t + r_t * r_t
        sum2_t = sum2_t + r_t * s(i, j)
      end do
      !$acc atomic
      tmp = tmp + sum_t
      !$acc atomic
      w = w + sum2_t
    end do

    !$acc host_data use_device(tmp, w)
    NCCLCHECK(ncclGroupStart())
    NCCLCHECK(ncclAllReduce(tmp, tmp, 1, ncclFloat64, ncclSum, nccl_comm, stream))
    NCCLCHECK(ncclAllReduce(w, w, 1, ncclFloat64, ncclSum, nccl_comm, stream))
    NCCLCHECK(ncclGroupEnd())
    !$acc end host_data

    !$acc kernels async(async_id) present(w, tmp, tp)
    w = w/tmp
    tp = 0.d0
    !$acc end kernels
    
    ! x = w * sh + alpha * ph + x; r = -w * t + s; inprod(r, r, tp)
    !$acc parallel loop gang vector_length(tiling_size) async(async_id) private(sum_t) present(w, alpha)
    do n = 0, nx*ny - 1, tiling_size
      sum_t = 0.d0
      !$acc loop vector reduction(+:sum_t) private(r_t)
      do l = 0, tiling_size - 1
        i = mod(n + l, nx) + 1
        j = (n + l) / nx + 1
        x(i, j) = w*sh(i, j) + alpha*ph(i, j) + x(i, j)
        r_t = -w*t(i, j) + s(i, j)
        r(i, j) = r_t
        sum_t = sum_t + r_t * r_t
      end do
      !$acc atomic
      tp = tp + sum_t
    end do
    
    !$acc host_data use_device(tp)
    NCCLCHECK(ncclAllReduce(tp, tp, 1, ncclFloat64, ncclSum, nccl_comm, stream))
    !$acc end host_data
    
    !$acc kernels async(async_id) present(rho, rhn, res2, tp, res3, res1)
    rho = rhn
    res2 = dsqrt(tp)
    res3 = res2 / res1
    rhn = 0.d0
    !$acc end kernels

    !$acc update self(res2, res3) async(async_id)
    !$acc wait(async_id)
    
    if((res2 .lt. threshold_bicg) .or. (res3 .lt. 1.d-6)) then
      call ghost_cell_exch_gpu(m_comm_cart, r8type2d, nbid, ndim, x, symm_np, async_id)
      ! if(myid .eq. rootid) then
      !   write(log_info,*) "bicg pressure solver convergence:", iter, res2, res3
      !   call comm_write_log_info(fid_log, log_info)
      ! end if
      converged = .true.
      exit
    end if
  end do

  if (.not. converged) then
    if(myid .eq. rootid) then
      write(log_info,*) "bicg pressure solver not convergence:", iter, res2, res3
      call comm_write_log_info(fid_log, log_info)
    end if
    call comm_finalize
  end if

end subroutine p_bicgstab_le_gpu

subroutine pcsi_gpu(ab,al,ac,ar,at,b,x,numiterations,rmsresidual, async_id)

  use openacc
  use hyperlink, only: nx, ny, in, &
                       m_comm_cart, r8type2d, nbid, ndim, myid, symm_np, &
                       total_area, pcsicriterion  ! NTU 20250324 M. Hsieh
  use cudafor
  use iso_c_binding
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

  integer, intent(in) :: async_id
  real(r8), dimension(nx,ny), intent(in) :: &
    ab, al, ac, ar, at, &
    b    ! right hand side of linear system

  real(r8), dimension(0:nx+1,0:ny+1), intent(inout) :: &
    x    ! on input,  an initial guess for the solution
         ! on output, solution of the linear system

  ! OUTPUT PARAMETERS:

  integer, intent(out) :: &
    numiterations  ! accumulated no of iterations (diagnostic)

  real(r8), intent(out) :: rmsresidual  ! residual (diagnostic)

  ! local variables
  integer(kind=cuda_stream_kind) :: stream, check_stream
  integer  :: i,j,m,n,l  ! local iteration counter
  integer  :: curr, next
  logical :: converged

  real(r8) :: cc1, tmp, sum_t, s_t, r_t, q_t, &
              residualnorm     ! residual normalization

  integer :: ierr, istat(mpi_status_size)
  character(len=256) :: log_info

  stream = acc_get_cuda_stream(async_id)
  check_stream = acc_get_cuda_stream(async_id + 1)

  if (.not. pcsi_cuda_init) then
    CUDACHECK(cudaEventCreate(event))
    CUDACHECK(cudaHostAlloc(h_ptr_rr, sizeof(1.0d0), cudaHostAllocMapped))
    call c_f_pointer(h_ptr_rr, zc_rr)
    CUDACHECK(cudaHostGetDevicePointer(d_ptr_rr, h_ptr_rr, 0))
    call c_f_pointer(d_ptr_rr, d_zc_rr)
    pcsi_cuda_init = .true.
  end if

  converged = .false.
  numiterations = maxiterations

  call pcsi_init(ab, al, ac, ar, at, b, x, async_id)
  call pcsi_stage1(1, ab, al, ac, ar, at, b, x, async_id)

  CUDACHECK(cudaEventRecord(event, stream))

  iterationloop: do m = PCSICheckFreq, maxiterations, PCSICheckFreq

    CUDACHECK(cudaStreamWaitEvent(check_stream, event, 0))
    call pcsi_stage2(mod(m/PCSICheckFreq, 2) + 1, async_id + 1)
    call pcsi_stage1(mod(m/PCSICheckFreq, 2) + 1, ab, al, ac, ar, at, b, x, async_id)
    CUDACHECK(cudaEventRecord(event, stream))

    !$acc wait(async_id + 1)
    call MPI_ALLREDUCE(zc_rr, rr, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
    if (rr < pcsicriterion) then
      numiterations = m
      converged = .true.
      exit iterationloop
    end if

  end do iterationloop

  if (.not. converged) then
    if(myid .eq. rootid) then
      write (*,*) "pcsi solver not convergence:", m, rr
      write(log_info,*) "pcsi solver not convergence:", m, rr
      call comm_write_log_info(fid_log, log_info)
    end if
    call comm_finalize
  end if

  residualnorm = 1.0d0/total_area
  rmsresidual = sqrt(rr*residualnorm)

end subroutine pcsi_gpu

subroutine pcsi_init(ab, al, ac, ar, at, b, x, async_id)
  use openacc
  use hyperlink, only: nx, ny, in, &
                       m_comm_cart, r8type2d, nbid, ndim, myid, symm_np, &
                       total_area, pcsicriterion
  use cudafor
  use iso_c_binding
  implicit none

  integer, intent(in) :: async_id
  real(r8), dimension(nx,ny), intent(in) :: ab, al, ac, ar, at, b
  real(r8), dimension(0:nx+1,0:ny+1), intent(inout) :: x

  integer(kind=cuda_stream_kind) :: stream
  integer :: i, j

  stream = acc_get_cuda_stream(async_id)

  if (.not. pcsi_init_cg_created) then
    if (pcsi_cuda_graph) then
      CUDACHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal))
    end if

    !$acc serial async(async_id) present(csalpha, csbeta, csy, csomega, one_csy)
    csalpha = 2.0d0/(pcsimaxeig-pcsimineig)
    csbeta = (pcsimaxeig+pcsimineig)/(pcsimaxeig-pcsimineig)
    csy = csbeta/csalpha
    csomega = 2.0d0/csy
    one_csy = 1.0d0/csy
    rr = 0.0d0
    !$acc end serial

    !$acc parallel loop collapse(2) async(async_id)
    do j = 1, ny
      do i = 1, nx
        !--- diagonal preconditioner if preconditioner not specified
        if ( .not. usepreconditioner ) then
          if ( abs(ac(i,j)) > 1.0d-18 ) then
            a0r(i,j) = 1.0d0 / ac(i,j)
          else
            a0r(i,j) = 0.0d0
          end if
        end if
      end do
    end do

    call ghost_cell_exch_gpu(m_comm_cart, r8type2d, nbid, ndim, x, symm_np, async_id)

    !$acc parallel loop collapse(2) async(async_id) present(one_csy)
    do j = 1, ny
      do i = 1, nx
        s(i,j) = ab(i,j)*x(i,j-1) &
                + al(i,j)*x(i-1,j) &
                + ac(i,j)*x(i,j)   &
                + ar(i,j)*x(i+1,j) &
                + at(i,j)*x(i,j+1)
        r(i,j,1) = b(i,j) - s(i,j)
        if (usepreconditioner) then
          s(i, j) = r(i, j, 1) / ac(i, j)
          r(i,j,1) = s(i,j)
        else
          r(i,j,1) = r(i,j,1) * a0r(i,j)
        end if
        q(i,j) = one_csy * r(i,j,1)
      end do
    end do

    !$acc parallel loop collapse(2) async(async_id)
    do j = 1, ny
      do i = 1, nx
        x(i,j) = x(i,j) + q(i,j)
      end do
    end do

    call ghost_cell_exch_gpu(m_comm_cart, r8type2d, nbid, ndim, x, symm_np, async_id)

    !$acc parallel loop collapse(2) async(async_id)
    do j = 1, ny
      do i = 1, nx
        s(i,j) = ab(i,j)*x(i,j-1) &
                + al(i,j)*x(i-1,j) &
                + ac(i,j)*x(i,j)   &
                + ar(i,j)*x(i+1,j) &
                + at(i,j)*x(i,j+1)
        r(i,j,1) = b(i,j) - s(i,j)
      end do
    end do

    if (pcsi_cuda_graph) then
      CUDACHECK(cudaStreamEndCapture(stream, pcsi_init_graph))
      CUDACHECK(cudaGraphInstantiate(pcsi_init_graph_exec, pcsi_init_graph, 0))
      pcsi_init_cg_created = .true.
    end if
  end if

  if (pcsi_cuda_graph) then
    CUDACHECK(cudaGraphLaunch(pcsi_init_graph_exec, stream))
  end if

end subroutine pcsi_init

subroutine pcsi_stage1(curr, ab, al, ac, ar, at, b, x, async_id)
  use openacc
  use hyperlink, only: nx, ny, in, &
                       m_comm_cart, r8type2d, nbid, ndim, myid, symm_np, &
                       total_area, pcsicriterion
  use cudafor
  use iso_c_binding
  implicit none

  integer, intent(in) :: curr, async_id
  real(r8), dimension(nx,ny), intent(in) :: ab, al, ac, ar, at, b
  real(r8), dimension(0:nx+1,0:ny+1), intent(inout) :: x

  integer(kind=cuda_stream_kind) :: stream
  integer :: i, j, m, next
  real(r8) :: s_t, r_t, q_t, cc1

  stream = acc_get_cuda_stream(async_id)

  if (.not. pcsi_stage1_cg_created(curr)) then
    next = curr

    if (pcsi_cuda_graph) then
      CUDACHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal))
    end if

    do m = 1, PCSICheckFreq
      !$acc serial async(async_id) present(csomega, csy, csalpha)
      csomega = 1.0d0/(csy-csomega/(4.0d0*csalpha*csalpha))
      !$acc end serial

      !$acc parallel loop gang vector_length(128) async(async_id) &
      !$acc& present(csalpha, csomega, csy, r)
      do j = 1, ny
        !$acc loop vector private(s_t, r_t, q_t, cc1)
        do i = 1, nx
          if (usepreconditioner) then
            s_t = r(i, j, curr) / ac(i, j)
            r_t = s_t
          else
            r_t = r(i, j, curr) * a0r(i,j)
          end if
          cc1 = csy*csomega-1.0d0
          q_t = csomega*r_t + cc1*q(i,j)
          q(i,j) = q_t
          x(i,j) = x(i,j) + q_t
        end do
      end do

      call ghost_cell_exch_gpu(m_comm_cart, r8type2d, nbid, ndim, x, symm_np, async_id)

      if (m == PCSICheckFreq) then
        next = 3 - curr
      end if

      !$acc parallel loop collapse(2) async(async_id) private(s_t) present(r)
      do j = 1, ny
        do i = 1, nx
          s_t = ab(i,j)*x(i,j-1) &
                + al(i,j)*x(i-1,j) &
                + ac(i,j)*x(i,j) &
                + ar(i,j)*x(i+1,j) &
                + at(i,j)*x(i,j+1)
          s(i,j) = s_t
          r(i,j,next) = b(i,j) - s_t
        end do
      end do
    end do

    if (pcsi_cuda_graph) then
      CUDACHECK(cudaStreamEndCapture(stream, pcsi_stage1_graph(curr)))
      CUDACHECK(cudaGraphInstantiate(pcsi_stage1_graph_exec(curr), pcsi_stage1_graph(curr), 0))
      pcsi_stage1_cg_created(curr) = .true.
    end if
  end if

  if (pcsi_cuda_graph) then
    CUDACHECK(cudaGraphLaunch(pcsi_stage1_graph_exec(curr), stream))
  end if

end subroutine pcsi_stage1

subroutine pcsi_stage2(curr, async_id)
  use openacc
  use hyperlink, only: nx, ny, in, &
                       m_comm_cart, r8type2d, nbid, ndim, myid, symm_np, &
                       total_area, pcsicriterion
  use cudafor
  use iso_c_binding
  implicit none

  integer, intent(in) :: curr, async_id

  integer(kind=cuda_stream_kind) :: stream
  integer :: i, j, n, l
  real(r8) :: sum_t, r_t

  integer, parameter :: tiling_size = 128

  stream = acc_get_cuda_stream(async_id)

  if (.not. pcsi_stage2_cg_created(curr)) then
    if (pcsi_cuda_graph) then
      CUDACHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal))
    end if

    !$acc parallel loop gang vector_length(tiling_size) async(async_id) &
    !$acc& private(sum_t)
    do n = 0, nx*ny - 1, tiling_size
      sum_t = 0.d0
      !$acc loop vector reduction(+:sum_t) private(r_t)
      do l = 0, tiling_size - 1
        i = mod(n + l, nx) + 1
        j = (n + l) / nx + 1
        r_t = r(i,j,curr)
        sum_t = sum_t + r_t*r_t*dble(in(i,j,1))
      end do
      !$acc atomic
      rr = rr + sum_t
    end do

    !$acc serial async(async_id) present(rr) deviceptr(d_zc_rr)
    d_zc_rr = rr
    rr = 0.0d0
    !$acc end serial

    if (pcsi_cuda_graph) then
      CUDACHECK(cudaStreamEndCapture(stream, pcsi_stage2_graph(curr)))
      CUDACHECK(cudaGraphInstantiate(pcsi_stage2_graph_exec(curr), pcsi_stage2_graph(curr), 0))
      pcsi_stage2_cg_created(curr) = .true.
    end if
  end if

  if (pcsi_cuda_graph) then
    CUDACHECK(cudaGraphLaunch(pcsi_stage2_graph_exec(curr), stream))
  end if

end subroutine pcsi_stage2

end module timcom_solver_gpu
