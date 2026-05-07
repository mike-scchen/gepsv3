module timcom_rsm_cpl
  use tai_timcom_const
  use tai_timcom_comm
  use tai_timcom_grid
  use tai_timcom_general
  implicit none
 
contains

subroutine import_rsm_grid(itf, igrd, jgrd, rsm_field, tsea)
  use tai_hyperlink, only: assign_grid_all, nx, ny, myid_x, myid_y
  implicit none

  integer :: n, itf, igrd, jgrd
  real(r8) :: rsm_field(igrd,jgrd,9), tsea(igrd,jgrd)
  integer :: ibgn, iend, jbgn, jend

  n = 1
  call assign_grid_all(ocn_grid(n)) 
  ibgn = 1
  iend = igrd-1
  jbgn = 1
  jend = jgrd-1

  if(myid_x .eq. 0) then
    ibgn = 2
    iend = igrd
  end if
  if(myid_y .eq. 0) then
    jbgn = 2
    jend = jgrd
  end if

  ocn_grid(n)%u_10 = rsm_field(ibgn:iend, jbgn:jend, 1)
  ocn_grid(n)%v_10 = rsm_field(ibgn:iend, jbgn:jend, 2)
  ocn_grid(n)%t_10 = rsm_field(ibgn:iend, jbgn:jend, 3)
  ocn_grid(n)%q_10 = rsm_field(ibgn:iend, jbgn:jend, 4)
  ocn_grid(n)%pslv = rsm_field(ibgn:iend, jbgn:jend, 5)
  ocn_grid(n)%swup = rsm_field(ibgn:iend, jbgn:jend, 6)
  ocn_grid(n)%swdn = rsm_field(ibgn:iend, jbgn:jend, 7)
  ocn_grid(n)%lwdn = rsm_field(ibgn:iend, jbgn:jend, 8)
  ocn_grid(n)%rain = rsm_field(ibgn:iend, jbgn:jend, 9)
  ocn_grid(n)%snow = 0.d0 

end subroutine import_rsm_grid

subroutine export_rsm_grid(igrd, jgrd, tsea, ssu, ssv)
  use tai_hyperlink, only: assign_grid_all, myid_x, myid_y, nx, ny
  implicit none

  integer :: n, igrd, jgrd
  real(r8), dimension(igrd, jgrd) :: tsea, ssu, ssv
  integer :: ibgn, iend, jbgn, jend, i, j

  ssu = 0.d0
  ssv = 0.d0

  n = 1
  call assign_grid_all(ocn_grid(n))
  ibgn = 1
  iend = igrd-1
  jbgn = 1
  jend = jgrd-1

  if(myid_x .eq. 0) then
    ibgn = 2
    iend = igrd
  end if
  if(myid_y .eq. 0) then
    jbgn = 2
    jend = jgrd
  end if

  do j = 1, ny
    do i = 1, nx
      if(ocn_grid(n)%in(i,j,1) .eq. 1) then
        tsea(ibgn+i-1,jbgn+j-1) = ocn_grid(n)%t2(i,j,1) + 273.15d0
         ssu(ibgn+i-1,jbgn+j-1) = ocn_grid(n)%u2(i,j,1)*1.d-2
         ssv(ibgn+i-1,jbgn+j-1) = ocn_grid(n)%v2(i,j,1)*1.d-2
      end if
    end do
  end do

end subroutine export_rsm_grid

end module timcom_rsm_cpl
