module tai_timcom_solver
  use tai_timcom_const
  use tai_timcom_comm
  implicit none
contains
subroutine sipinc(al,ab,ac,ar,at,cl,cb,cc,cr,ct,alpha)
  use tai_hyperlink, only: nx, ny
  implicit none

  real(r8), dimension(nx,ny), intent(in)  :: al, ab, ac, ar, at
  real(r8), intent(in) :: alpha
  real(r8), dimension(nx,ny), intent(out) :: cl, cb, cr, cc, ct
  real(r8), dimension(nx,ny) :: act
  integer :: i, j

  act = ac

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
end subroutine sipinc

subroutine p_bicgstab_le(ab,al,ac,ar,at,b,x,cb,cl,cc,cr,ct,r,rh,p,v,s,t,ph,sh,iter,res2,res3,itf)
  use tai_hyperlink, only: nx, ny, m_comm_cart, r8type2d, nbid, ndim, myid, symm_np, threshold_bicg
  implicit none

  real(r8),dimension(nx,ny) :: &
    ab, al, ac, ar, at, &
    cb, cl, cc, cr, ct, b

  real(r8),dimension(0:nx+1,0:ny+1) :: &
    x, r, rh, p, v, s, t, ph, sh

  real(r8) :: res1, res2, res3, tp, tmp, rho, rhn, alpha, beta, w
  integer :: ierr, istat(MPI_STATUS_SIZE), itf
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
      call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, x, .false.)
    !  if(myid .eq. rootid .and. itf .eq. 2881) then
    !    write(log_info,*) "bicg pressure solver convergence:", iter, res2, res3
    !    call comm_write_log_info(fid_log, log_info)
    !   end if
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
  use tai_hyperlink, only: m_comm_cart, nx, ny, r8type2d, nbid, ndim, symm_np
                     
  implicit none
  real(r8), dimension(nx,ny) :: ab, al, ac, ar, at
  real(r8), dimension(0:nx+1,0:ny+1) :: x, b
  integer :: ierr, istat(MPI_STATUS_SIZE)

  call MPI_BARRIER(m_comm_cart, ierr)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, x, .false.)

  b(1:nx,1:ny) = ab*x(1:nx,  0:ny-1) &
               + al*x(0:nx-1,1:ny)   &
               + ac*x(1:nx,  1:ny)   &
               + ar*x(2:nx+1,1:ny)   &
               + at*x(1:nx,  2:ny+1) 

end subroutine p_sqprod

subroutine p_inprod(a1,a2,pd)
  use tai_hyperlink, only: nx, ny, m_comm_cart
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
  use tai_hyperlink, only: nx, ny
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

end module tai_timcom_solver

