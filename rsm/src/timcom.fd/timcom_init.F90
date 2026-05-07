module tai_timcom_init
  use tai_timcom_const
  use tai_timcom_comm
  use tai_timcom_grid, only: ocn_grid, ocn_panel
  use tai_hyperlink
  implicit none

contains

subroutine grid_init(grid_list, compid)
  use tai_timcom_general, only: member, ocn_grid_num, str_split
  use tai_timcom_grid,    only: allocate_grid_var
  use tai_timcom_kpp,     only: init_vmix_kpp
  use tai_timcom_cplmct,  only: ocn_cpl_init
  implicit none

  character(len=*) :: grid_list
  integer, intent(in), optional :: compid
  character(len=256) :: log_info

  integer :: i, ierr, n

  if(ocn_grid_num .gt. 0) then
    allocate(ocn_grid(ocn_grid_num))
    do i = 1, ocn_grid_num
      call str_split(trim(grid_list), ":", ocn_grid(i)%tag, ierr)
      write(ocn_grid(i)%nml_file, '(a8,i3.3,a5,a)') 'namelist',member,'.inp_',trim(ocn_grid(i)%tag)

      if(myid_timcom .eq. rootid) then
        write(log_info,'(5a)') &
          '[INFO]: Register ocn grid: ', trim(ocn_grid(i)%tag), achar(10), &
          '        Read TIMCOM grid namelist: ', trim(ocn_grid(i)%nml_file)
        call comm_write_log_info(fid_log, log_info)
      end if

      call assign_grid_scalar(ocn_grid(i))
      call get_grid_inp_nml(trim(ocn_grid(i)%nml_file))
      call comm_cart_gen(ocn_grid(i)%mpicom)
      call comm_def_mpi_type(ocn_grid(i)%mpicom, nx, ny, nz)
      call allocate_grid_var(ocn_grid(i))
      call assign_grid_vector(ocn_grid(i))
      call get_grid_indata
      call set_grid_solver
      call init_lateral_bc(ocn_grid(i))
      call set_initial_condition
      call set_ts_nudging
      call init_vmix_kpp(ocn_grid(i)%tidal_energy)
      if(present(compid)) call ocn_cpl_init(compid, ocn_grid(i))
    end do
  else
    if(myid_timcom .eq. rootid) then
      write(log_info,'(a,a)') '[INFO]: There is no assigned ocn grid.', trim(grid_list)
      call comm_write_log_info(fid_log, log_info)
    end if
  end if
end subroutine grid_init

subroutine get_grid_inp_nml(nml_file)
  implicit none

  character(len=*) :: nml_file
  integer :: fid, ierr, freq_out, freq_restart
  character(len=4) :: freq_out_unit
  character(len=256) :: log_info

  namelist /grid_info/   nx_grid, ny_grid, nz_grid, daodt
  namelist /mpi_info/    npx, npy, peri_x, peri_y, symm_np
  namelist /io_info/     freq_out, freq_restart, freq_out_unit
  namelist /phy_info/    opt_windmix, opt_arbr_p0, dm0
  namelist /solver_info/ threshold_bicg, threshold_p0, max_iter_p0

  fid = myid_timcom + 20
  open(unit=fid, file=trim(nml_file), status='OLD')
  read(unit=fid, nml=grid_info,   iostat=ierr)
  read(unit=fid, nml=mpi_info,    iostat=ierr)
  read(unit=fid, nml=io_info,     iostat=ierr)
  read(unit=fid, nml=phy_info,    iostat=ierr)
  read(unit=fid, nml=solver_info, iostat=ierr)
  close(fid)

  nx = nx_grid/npx
  ny = ny_grid/npy
  nz = nz_grid

  nxf = nx+1
  nyf = ny+1 
  nzf = nz+1

  dt  = 2.d0*day2sec/dble(daodt)
  odt = 1.d0/dt

  select case(trim(freq_out_unit))
  case('days')
    daodt_out = daodt*freq_out
    daodt_restart = daodt*freq_restart
    daodt_out_unit = daodt
  case('hour')
    daodt_out = int(dble(daodt*freq_out)/24.d0)
    daodt_restart = int(dble(daodt*freq_restart)/24.d0)
    daodt_out_unit = daodt/24
  end select
 
  if(myid_timcom .eq. rootid) then
    write(log_info, *) '[INFO]: get_grid_inp_nml: ', trim(mytag)
    call comm_write_log_info(fid_log, log_info)
    write(log_info, *) '        grid dims: ', nx_grid, ny_grid, nz_grid
    call comm_write_log_info(fid_log, log_info)
    write(log_info, *) '        subgrid dims: ', nx, ny, nz
    call comm_write_log_info(fid_log, log_info)
    write(log_info, *) '        proces. dims: ', npx, npy
    call comm_write_log_info(fid_log, log_info)
    write(log_info, *) '        daodt: ', daodt
    call comm_write_log_info(fid_log, log_info)
    write(log_info, *) '        daodt restart: ', daodt_restart
    call comm_write_log_info(fid_log, log_info)
    write(log_info, *) '        dt and odt: ', dt, odt
    call comm_write_log_info(fid_log, log_info)
  end if
end subroutine get_grid_inp_nml

subroutine get_grid_indata
  use tai_timcom_pncio, only: get_indata_pncio
  implicit none

  character(len=128) :: pncfile
  real(r8) :: tmp, total_trans, R, xi1, xi2, myZ 
  integer :: i, j, k, ierr, ncid, varid
  character(len=256) :: log_info

  write(pncfile, '(a,a8)') trim(mytag),'_grid.nc'
  if(myid .eq. rootid) then
    write(log_info, '(2a)') "[INFO]: read indata file: ", trim(pncfile)
    call comm_write_log_info(fid_log, log_info)
  end if

  ndim(1) = nx
  ndim(2) = ny
  ndim(3) = nz

  in = 0
  kb = 0
  call get_indata_pncio(trim(pncfile))

  call ghost_cell_exch(m_comm_cart, i2type3d, nbid, ndim, in, .false.)
  call ghost_cell_exch(m_comm_cart, i2type2d, nbid, ndim, kb, .false.)

  dz = 1.d0/odz
  dzw(2:nz) = z_grid(2:nz) - z_grid(1:nz-1)
  dzw(1) = z_grid(1) - z_face(1)
  dzw(nzf) = z_face(nzf) - z_grid(nz)
  odzw = 1.d0/dzw
  tanphi = dtan(y_grid*d2r)/r0
  curv_f = dsin(y_grid*d2r)*pi/2.16d4

  iu = 0
  iu(0:nxf+1,0:ny+1 ,1:nz) = in(-1:nx+1, 0:ny+1,1:nz  )*in(0:nx+2,0:ny+1,1:nz)
 
  iv = 0
  iv(0:nx+1 ,0:nyf+1,1:nz) = in( 0:nx+1,-1:ny+1,1:nz  )*in(0:nx+1,0:ny+2,1:nz)
  
  iw = 0
  iw(1:nx,1:ny,2:nz) = in(1:nx,1:ny,1:nz-1)*in(1:nx,1:ny,2:nz)

  if(lfsrf .eq. 1) iw(1:nx,1:ny,1) = in(1:nx,1:ny,1) 

  xu_bgn = 1
  xu_end = nxf
  if(.not. peri_x) then
    if(myid_x .eq. 0)     xu_bgn = 2
    if(myid_x .eq. npx-1) xu_end = nxf-1
  end if

  yv_bgn = 1
  yv_end = nyf
  if(.not. peri_y) then
    if(myid_y .eq. 0)     yv_bgn = 2
    if(myid_y .eq. npy-1) yv_end = nyf-1
  end if

  do j = 1, ny
    area(:,j) = dx(j)*dy(j)
    do k = 1, nz
      volume(:,j,k) = area(:,j)*dz(k)
    end do
  end do

  volume = volume*dble(in(1:nx,1:ny,1:nz))
  total_vol = sum(volume)
  tmp = total_vol
  call MPI_ALLREDUCE(tmp, total_vol, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

  area = area*dble(in(1:nx,1:ny,1))
  total_area = sum(area)
  tmp = total_area
  call MPI_ALLREDUCE(tmp, total_area, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

  avg_glb_t = 0.d0
  avg_glb_s = 0.d0

  R   = 0.58d0
  xi1 = 0.35d0
  xi2 =  23.d0
  do j = 1, ny
    do i = 1, nx
      total_trans=0.d0
      do k = 1, nz
        myZ = -z_grid(k)/100.d0
        sw_trans(i,j,k)=(R*dexp(myZ/xi1)+(1-R)*dexp(myZ/xi2))* &
                               in(i,j,k)*odz(k)
        total_trans = total_trans + sw_trans(i,j,k)
      end do
      do k = 1, nz
        if( total_trans .gt. 1.d-6) then
          sw_trans(i,j,k) = sw_trans(i,j,k)/total_trans
        else
          sw_trans(i,j,k) = 0.d0
        end if
      end do
    end do
  end do

  sw_trans = sw_trans*dble(in(1:nx,1:ny,1:nz))
! ====================================================
! set background vertical eddy viscosity & diffusivity
! ====================================================
! augment molecular viscosity & diffusivity by parameterized synopti! wind
! events, and by breaking (u,v,t,s, near surface only) and non-breaking
! (u,v only, at all rmetss) internal waves.
! synopti! wind forced mixing having 20m e-folding scale.
! bigger momentum mixing due to pressure-xfers
! associated with internal waves that have no counterpart in t,s xfers.
! after yr 7
! 10m scale height for t,s mixing due to breaking internal waves
! deeper tmpuv enhancement is motivated by findings of hurlburt & hogan
! that deep bathymetry and deep transient eddies influences gs path.
! tmpts is not enhanced (layer model has none across layer interface!)
! 1000m scale height for u,v mixing (due to pressure-xfers even without
! internal wave breaking). this enhances deep bathymetry effe!cts on
! upper level flow (e.g., gs separation, kuroshio meanders)
  do k = 1, nz-1
    vbk(k) = 1.d-2 + 5.d-2*dexp(-1.d-5*(z_face(k+1) - z_face(1)))
    hbk(k) = 2.d-3 + 2.d-2*dexp(-1.d-3*(z_face(k+1) - z_face(1)))
  end do 

  call set_lateral_mixing
  dhx = dmx
  dhy = dmy
end subroutine get_grid_indata

subroutine set_lateral_mixing
  implicit none

  integer :: i, j, k, it, jt
  real(r8) :: temp, scr(-1:nx+2,-1:ny+2,nz), scr2(-1:nx+2,-1:ny+2,nz)

  scr = 0.d0

  it = myid_x*nx
  jt = myid_y*ny

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        temp = dm0*0.2d0*(dexp(-1.d-3*(2-j-jt)**2)+dexp(-1.d-3*(ny_grid-jt-j)**2))      ! north & south poles
        scr(i,j,k) = dm0*0.2d0*dexp(-1.d-3*(y_grid(j)**2)) + temp
      end do
    end do
  end do

  if(y_grid(1) .gt. 80.d0) scr = scr*100.d0
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, scr, symm_np)

  scr = merge(scr*1.d3, scr, in.eq.0)
  scr2 = scr

  do k = 1, 3
    scr(1:nx,1:ny,1:nz) = 0.5d0*scr2(1:nx,1:ny,1:nz) &
                        + 0.125d0*(scr2(0:nx-1,1:ny,1:nz) + scr2(2:nx+1,1:ny,1:nz) &
                                  +scr2(1:nx,0:ny-1,1:nz) + scr2(1:nx,2:ny+1,1:nz))
    call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, scr, symm_np)
    scr2 = scr
  end do

  dmx = scr(1:nxf,1:ny,1:nz)
  dmy = scr(1:nx,1:nyf,1:nz)
!  do j = 1, ny
!    do i = 1, nxf
!      dmx(i,j,:) = 0.375d0*(scr(i-1,j,:)+scr(i,j,:)) + 0.125d0*(scr(i-2,j,:)+scr(i+1,j,:))
!    end do
!  end do

!  do j = 1, nyf
!    do i = 1, nx
!      dmy(i,j,:) = 0.375d0*(scr(i,j,:)+scr(i,j-1,:)) + 0.125d0*(scr(i,j-2,:)+scr(i,j+1,:))
!    end do
!  end do
end subroutine set_lateral_mixing

subroutine set_grid_solver
  use tai_timcom_solver
  implicit none
  
  integer :: i, j, k, ii, jj
  real(r8) :: gamma_tmp, dzx, dzyb, dzyt 
 
  k = 1
  do j = 1, ny
    jj = j+1
    dzx  = dz(k)*odx(j)**2
    dzyb = csv(j )*ocs(j)*odyv(j )*ody(j)*dz(k)
    dzyt = csv(jj)*ocs(j)*odyv(jj)*ody(j)*dz(k)
    do i = 1, nx
      al(i,j) = dzx 
      ar(i,j) = dzx 
      ab(i,j) = dzyb
      at(i,j) = dzyt
    end do
  end do

  do k = 2, nz
    do j = 1, ny
      jj = j+1
      dzx  = dz(k)*odx(j)**2
      dzyb = csv(j )*ocs(j)*odyv(j )*ody(j)*dz(k)
      dzyt = csv(jj)*ocs(j)*odyv(jj)*ody(j)*dz(k)
      do i = 1, nx
        ii = i+1
        al(i,j) = al(i,j) + dzx *iu( i, j,k)
        ar(i,j) = ar(i,j) + dzx *iu(ii, j,k)
        ab(i,j) = ab(i,j) + dzyb*iv( i, j,k)
        at(i,j) = at(i,j) + dzyt*iv( i,jj,k)
      end do
    end do
  end do

  if(.not. peri_x) then
    if(myid_x .eq. 0) al(1,1:ny) = 0.d0
    if(myid_x .eq. npx-1) ar(nx,1:ny) = 0.d0
  end if

  if(.not. peri_y) then
    if(myid_y .eq. 0) ab(1:nx,1) = 0.d0
    if(myid_y .eq. npy-1) at(1:nx,ny) = 0.d0
  end if

  gamma_tmp = 0.d0
  if (lfsrf .eq. 1) gamma_tmp = 1.d0/(0.5d0*dt**2)
  do j = 1, ny
    do i = 1, nx
      al(i,j) = -al(i,j)
      ar(i,j) = -ar(i,j)
      ab(i,j) = -ab(i,j)
      at(i,j) = -at(i,j)
      ac(i,j) = (-al(i,j)-ar(i,j)-ab(i,j)-at(i,j)) + gamma_tmp
    end do
  end do

  call sipinc(al, ab, ac, ar, at, cl, cb, cc, cr, ct, 0.5d0)
end subroutine set_grid_solver

subroutine init_lateral_bc(grd)
  use tai_timcom_pncio, only: get_lateral_bc_info, get_lateral_bc_field, check_lateral_bc
  use tai_timcom_general, only: member
  implicit none

  type(ocn_panel) :: grd
  character(len=256) :: field_file

  call get_lateral_bc_info(grd%jwbc, grd%jebc, grd%isbc, grd%inbc, grd%mpicom%ns_lbc, grd%mpicom%nc_lbc)

  field_file = "./TIMCOM_glb/glb000_000000.nc"

  call get_lateral_bc_field(trim(field_file), grd%mpicom%ns_lbc, grd%mpicom%nc_lbc, &
                            grd%jwbc, grd%jebc, grd%isbc, grd%inbc, &
                            grd%wbc, grd%ebc, grd%sbc, grd%nbc)
!  field_file = "tai000_lbc_000000.nc"
!  call check_lateral_bc(trim(field_file), grd%wbc, grd%ebc, grd%sbc, grd%nbc)

end subroutine init_lateral_bc

subroutine set_initial_condition
  use tai_timcom_general, only: restart
  implicit none

  if(restart) then
    call set_initial_condition_restart
  else
    call set_initial_condition_newrun
  end if
end subroutine set_initial_condition

subroutine set_initial_condition_restart
  use tai_timcom_general, only: member, path_restart, timer_drv
  use tai_timcom_pncio, only: get_restart_pncio, put_ic_pncio
  implicit none

  character(len=256) :: restart_file
  integer :: i, j, k
  real(r8) :: tmpd

  write(restart_file,'(a,a,i3.3,a,i4.4,i2.2,i2.2,i2.2,a)') &
    trim(path_restart), trim(mytag), member,'_restart_',timer_drv%greg1(1:4),'.nc'
 
  call get_restart_pncio(restart_file)

  u1 = u1*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, u1,  symm_np)
  u2 = u2*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, u2,  symm_np)
  ulf=ulf*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, ulf, symm_np)

  v1 = v1*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, v1,  symm_np)
  v2 = v2*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, v2,  symm_np)
  vlf=vlf*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, vlf, symm_np)

  t1 = t1*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, t1,  symm_np)
  t2 = t2*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, t2,  symm_np)
  tlf=tlf*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, tlf, symm_np)

  s1 = s1*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, s1,  symm_np)
  s2 = s2*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, s2,  symm_np)
  slf=slf*in
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, slf, symm_np)

  u(xu_bgn:xu_end,1:ny,1:nz) = u(xu_bgn:xu_end,1:ny,1:nz)*iu(xu_bgn:xu_end,1:ny,1:nz)
  call mpi_exch_r8type3du(m_comm_cart, r8type3du, nbid, ndim, u, .false.)
  v(1:nx,yv_bgn:yv_end,1:nz) = v(1:nx,yv_bgn:yv_end,1:nz)*iv(1:nx,yv_bgn:yv_end,1:nz)
  call mpi_exch_r8type3dv(m_comm_cart, r8type3dv, nbid, ndim, v, .false.)
  w(1:nx,1:ny,1:nzf) = w(1:nx,1:ny,1:nzf)*iw(1:nx,1:ny,1:nzf)

  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, p0,  .false.)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, x,   .false.)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, cgr, .false.)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, cgrh,.false.)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, cgp, .false.)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, cgv, .false.)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, cgs, .false.)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, cgt, .false.)

  if(.not. peri_x) then
    if(myid_x .eq. 0) then
      u2(0,1:ny,1:nz) = 2.d0*u(1,1:ny,1:nz) - u2(1,1:ny,1:nz)
      u1(0,1:ny,1:nz) = u2(0,1:ny,1:nz)
      ulf(0,1:ny,1:nz) = u2(0,1:ny,1:nz)

      v2(0,1:ny,1:nz) = v2(1,1:ny,1:nz)
      v1(0,1:ny,1:nz) = v2(0,1:ny,1:nz)
      vlf(0,1:ny,1:nz) = v2(0,1:ny,1:nz)

      t2(0,1:ny,1:nz) = t2(1,1:ny,1:nz)
      t1(0,1:ny,1:nz) = t2(0,1:ny,1:nz)
      tlf(0,1:ny,1:nz) = t2(0,1:ny,1:nz)

      s2(0,1:ny,1:nz) = s2(1,1:ny,1:nz)
      s1(0,1:ny,1:nz) = s2(0,1:ny,1:nz)
      slf(0,1:ny,1:nz) = s2(0,1:ny,1:nz)
    end if

    if(myid_x .eq. npx-1) then
      u2(nx+1,1:ny,1:nz) = 2.d0*u(nxf,1:ny,1:nz) - u2(nx,1:ny,1:nz)
      u1(nx+1,1:ny,1:nz) = u2(nx+1,1:ny,1:nz)
      ulf(nx+1,1:ny,1:nz) = u2(nx+1,1:ny,1:nz)

      v2(nx+1,1:ny,1:nz) = v2(nx,1:ny,1:nz)
      v1(nx+1,1:ny,1:nz) = v2(nx+1,1:ny,1:nz)
      vlf(nx+1,1:ny,1:nz) = v2(nx+1,1:ny,1:nz)

      t2(nx+1,1:ny,1:nz) = t2(nx,1:ny,1:nz)
      t1(nx+1,1:ny,1:nz) = t2(nx+1,1:ny,1:nz)
      tlf(nx+1,1:ny,1:nz) = t2(nx+1,1:ny,1:nz)
  
      s2(nx+1,1:ny,1:nz) = s2(nx,1:ny,1:nz)
      s1(nx+1,1:ny,1:nz) = s2(nx+1,1:ny,1:nz)
      slf(nx+1,1:ny,1:nz) = s2(nx+1,1:ny,1:nz)
    end if
  end if

  if(.not. peri_y) then
    if(myid_y .eq. 0) then
      u2(1:nx,0,1:nz) = u2(1:nx,1,1:nz)
      u1(1:nx,0,1:nz) = u2(1:nx,0,1:nz)
      ulf(1:nx,0,1:nz) = u2(1:nx,0,1:nz)

      v2(1:nx,0,1:nz) = 2.d0*v(1:nx,1,1:nz) - v2(1:nx,1,1:nz)
      v1(1:nx,0,1:nz) = v2(1:nx,0,1:nz)
      vlf(1:nx,0,1:nz) = v2(1:nx,0,1:nz)

      t2(1:nx,0,1:nz) = t2(1:nx,1,1:nz)
      t1(1:nx,0,1:nz) = t2(1:nx,0,1:nz)
      tlf(1:nx,0,1:nz) = t2(1:nx,0,1:nz)

      s2(1:nx,0,1:nz) = s2(1:nx,1,1:nz)
      s1(1:nx,0,1:nz) = s2(1:nx,0,1:nz)
      slf(1:nx,0,1:nz) = s2(1:nx,0,1:nz)
    end if

    if(myid_y .eq. npy-1 .and. .not. symm_np) then
      u2(1:nx,ny+1,1:nz) = u2(1:nx,ny,1:nz)
      u1(1:nx,ny+1,1:nz) = u2(1:nx,ny,1:nz)
      ulf(1:nx,ny+1,1:nz) = u2(1:nx,ny,1:nz)

      v2(1:nx,ny+1,1:nz) = 2.d0*v(1:nx,nyf,1:nz) - v2(1:nx,ny,1:nz)
      v1(1:nx,ny+1,1:nz) = v2(1:nx,ny+1,1:nz)
      vlf(1:nx,ny+1,1:nz) = v2(1:nx,ny+1,1:nz)

      t2(1:nx,ny+1,1:nz) = t2(1:nx,ny,1:nz)
      t1(1:nx,ny+1,1:nz) = t2(1:nx,ny+1,1:nz)
      tlf(1:nx,ny+1,1:nz) = t2(1:nx,ny+1,1:nz)

      s2(1:nx,ny+1,1:nz) = s2(1:nx,ny,1:nz)
      s1(1:nx,ny+1,1:nz) = s2(1:nx,ny+1,1:nz)
      slf(1:nx,ny+1,1:nz) = s2(1:nx,ny+1,1:nz)
    end if
  end if

  call put_ic_pncio("check_restart_grid.nc")

  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanT0 =', avg_glb_t)
  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanS0 =', avg_glb_s)
end subroutine set_initial_condition_restart

subroutine set_initial_condition_newrun
  use tai_timcom_general, only: member
  use tai_timcom_pncio, only: get_ssh_pncio, get_ic_pncio, put_ic_pncio
  implicit none

  character(len=128) :: ic_file, ic_ssh_file
  character(len=256) :: log_info
  integer :: i, j, k, ierr
  real(r8) :: tmpd, tmp, psm

  write(ic_file, '(a,i3.3,a3)') trim(mytag), member, '.nc'
  write(ic_ssh_file, '(a,i3.3,a7)') trim(mytag), member, '_ssh.nc'
  if(myid .eq. rootid) then
    write(log_info, '(2a)') "[INFO]: read IC file: ", trim(ic_file)
    call comm_write_log_info(fid_log, log_info)
  end if

  call get_ic_pncio(trim(ic_file))
  u1 = u1*in
  v1 = v1*in
  t1 = t1*in
  s1 = s1*in

  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, u1, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, v1, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, t1, symm_np)
  call ghost_cell_exch(m_comm_cart, r8type3d, nbid, ndim, s1, symm_np)

  p0(1:nx,1:ny) = ssh*rho_sw*grav*in(1:nx,1:ny,1)
  !psm = sum(p0(1:nx,1:ny)*area(1:nx,1:ny))
  !tmp = psm
  !call MPI_ALLREDUCE(tmp, psm, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
  !psm = psm/total_area
  !p0(1:nx,1:ny) = (p0(1:nx,1:ny)-psm)*in(1:nx,1:ny,1)
  call ghost_cell_exch(m_comm_cart, r8type2d, nbid, ndim, p0, symm_np)

  if(myid_x .eq. 0) then
    u1(0,:,:) = u1(1,:,:)
    v1(0,:,:) = v1(1,:,:)
    t1(0,:,:) = t1(1,:,:)
    s1(0,:,:) = s1(1,:,:)
    p0(0,:)   = p0(1,:)
  end if

  if(myid_x .eq. npx-1) then
    u1(nx+1,:,:) = u1(nx,:,:)
    v1(nx+1,:,:) = v1(nx,:,:)
    t1(nx+1,:,:) = t1(nx,:,:)
    s1(nx+1,:,:) = s1(nx,:,:)
    p0(nx+1,:)   = p0(nx,:)
  end if

  if(myid_y .eq. 0) then
    u1(:,0,:) = u1(:,1,:)
    v1(:,0,:) = v1(:,1,:)
    t1(:,0,:) = t1(:,1,:)
    s1(:,0,:) = s1(:,1,:)  
    p0(:,0)   = p0(:,1)
  end if

  if(myid_y .eq. npy-1) then
    u1(:,ny+1,:) = u1(:,ny,:)
    v1(:,ny+1,:) = v1(:,ny,:)
    t1(:,ny+1,:) = t1(:,ny,:)
    s1(:,ny+1,:) = s1(:,ny,:)
    p0(:,ny+1)   = p0(:,ny)
  end if

  call pre_def_uvw
  call mpi_exch_r8type3du(m_comm_cart, r8type3du, nbid, ndim, u, symm_np)
  call mpi_exch_r8type3dv(m_comm_cart, r8type3dv, nbid, ndim, v, symm_np)

  u2  = u1
  ulf = u2
  v2  = v1
  vlf = v2
  t2  = t1
  tlf = t2
  s2  = s1
  slf = s2

  call put_ic_pncio("check_ic_grid.nc")

!  call go_check_conservation(t2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanT0 =', avg_glb_t)
!  call go_check_conservation(s2(1:nx,1:ny,1:nz), volume, total_vol, 'glob_meanS0 =', avg_glb_s)
end subroutine set_initial_condition_newrun

subroutine set_ts_nudging
  use tai_timcom_pncio, only: get_ts_climate
  implicit none

  real(r8) :: &
    rate,    n_tmp, &
    tmpp_s,  tmpp_sc, tmpp_shs, &
    tmpp_t,  tmpp_tc, tmpp_ths, &
    tmpp_o,  tmpp_oc, tmpp_hi, tmpp_ei, &
    acc_deg, arc_deg, equ_deg,  &
    acc_rag, arc_rag, equ_rag,  &
    nudge_scale, hmean_scale

  integer :: i, j

  rate = 5.d0

! salinity
  n_tmp = -1.d0
  tmpp_s = max(1.d0/(n_tmp*daodt),0.d0)
! temperature
  n_tmp = -1.d0
  tmpp_t = max(1.d0/(n_tmp*daodt),0.d0)
! salinity coast
  n_tmp = 3.d0
  tmpp_sc = max(1.d0/(n_tmp*daodt),0.d0)
! temperature coast
  n_tmp = -1.d0
  tmpp_tc = max(1.d0/(n_tmp*daodt),0.d0)
! salinity acc and arc
  n_tmp = -1.d0
  tmpp_shs = max(1.d0/(n_tmp*daodt),0.d0)
! temperature acc and arc
  n_tmp = -1.d0
  tmpp_ths = max(1.d0/(n_tmp*daodt),0.d0)

! hmean acc and arc
  n_tmp = -1.d0
  tmpp_hi = max(1.d0/(n_tmp*daodt),0.d0)
! hmean equator
  n_tmp = -1.d0
  tmpp_ei = max(1.d0/(n_tmp*daodt),0.d0)
! hmean other
  n_tmp = -1.d0
  tmpp_o = max(1.d0/(n_tmp*daodt),0.d0)
! hmean coast
  n_tmp = -1.d0
  tmpp_oc = max(1.d0/(n_tmp*daodt),0.d0)

  acc_deg = -65.d0
  acc_rag = 2.d0

  arc_deg = 70.d0
  arc_rag = 2.d0

  equ_deg = 0.d0
  equ_rag = 5.d0

  do j = 1, ny
    nudge_scale = 0.d0
    hmean_scale = 0.d0
    if(y_grid(j) .ge. arc_deg .or. y_grid(j) .le. acc_deg) then
      nudge_scale = 1.d0
      hmean_scale = 1.d0
    elseif(y_grid(j) .ge. arc_deg - arc_rag) then
      nudge_scale = dexp(-rate*dabs((y_grid(j)-arc_deg)/arc_rag))
      hmean_scale = nudge_scale
    elseif(y_grid(j) .le. acc_deg + acc_rag) then
      nudge_scale = dexp(-rate*dabs((y_grid(j)-acc_deg)/acc_rag))
      hmean_scale = nudge_scale
    elseif(y_grid(j) .ge. equ_deg - equ_rag .and. y_grid(j) .le. equ_deg + equ_rag) then
      hmean_scale = dexp(-rate*2.5d0*dabs((y_grid(j)-equ_deg)/equ_rag))
      tmpp_hi = tmpp_ei
    end if

    t_nudge(:,j) = (tmpp_ths-tmpp_t)*nudge_scale + tmpp_t
    s_nudge(:,j) = (tmpp_shs-tmpp_s)*nudge_scale + tmpp_s
    hmean_nudge(:,j,:) = (tmpp_hi-tmpp_o)*hmean_scale + tmpp_o

    do i = 1, nx
      if(kb(i,j).le.1 .or. &
         kb(i-1,j) .le.1 .or. kb(i+1,j).le.1 .or. &
         kb(i,j-1) .le.1 .or. kb(i,j+1).le.1) then
        t_nudge(i,j) = tmpp_tc
        s_nudge(i,j) = tmpp_sc
        hmean_nudge(i,j,:) = tmpp_oc
      end if
    end do
  end do

  call get_ts_climate(nx, ny, nz, t_clim, s_clim)
end subroutine set_ts_nudging

subroutine pre_def_pressure
  use mpi
  use tai_timcom_eos, only: compute_eos_wright, &
                        compute_eos_mkcoef, &
                        compute_eos_jacmcd, &
                        compute_eos_stability
  implicit none

  real(r8), pointer :: tp(:,:,:), ts(:,:,:), pp(:,:,:)
  integer(2), pointer :: msk(:,:,:)
  real(r8), dimension(nx,ny,nz) :: rho1

  real(r8) :: tmpd, res0(1), res1(1), dum0(1), dum1(1)
  integer :: i, j, k, n, opt=2, ierr
  integer, parameter :: &
    opt_wright = 1, &
    opt_mkcoef = 2, &
    opt_JacMcD = 3
  character(len=256) :: log_info

  tp => t1(1:nx,1:ny,1:nz)
  ts => s1(1:nx,1:ny,1:nz)
  pp =>  p(1:nx,1:ny,1:nz)
  msk => in(1:nx,1:ny,1:nz)

  call compute_eos_stability(nx, ny, nz, msk, z_grid, tp, ts, rho)

  do n = 1, 10
    tmpd = grav*z_grid(1)
    do j = 1, ny
      do i = 1, nx
        pp(i,j,1) = tmpd*rho(i,j,1)
      end do
    end do

    do k = 2, nz
      tmpd = 0.5d0*grav*dzw(k)
      do j = 1, ny
        do i = 1, nx
          pp(i,j,k) = pp(i,j,k-1) + tmpd*(rho(i,j,k)+rho(i,j,k-1))
        end do
      end do
    end do

    select case(opt)
    case(opt_wright)
      call compute_eos_wright(nx, ny, nz, msk, tp, ts, pp, rho1)
    case(opt_mkcoef)
      call compute_eos_mkcoef(nx, ny, nz, msk, tp, ts, pp, rho1)
    case(opt_JacMcd)
      call compute_eos_JacMcd(nx, ny, nz, msk, tp, ts, -1.d-2*z_grid, rho1)
    end select

    res0(1) = sum((rho-dble(msk))**2)
    res1(1) = sum((rho1-rho)**2)

    call MPI_ALLREDUCE(res0, dum0, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
    call MPI_ALLREDUCE(res1, dum1, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)

    res1 = dsqrt(dum1)/dsqrt(dum0)

    if(myid .eq. rootid) then
      write(log_info,'(a,i4,f)') "[CHECK RHO]: ", n, res1
      call comm_write_log_info(fid_log, log_info)
    end if
    rho = rho1
    if(res1(1) < 1.d-6) exit
  end do

end subroutine pre_def_pressure

subroutine pre_def_velocity
  implicit none

  real(r8) :: &
    dpx(nx,ny,nz), &
    dpy(nx,ny,nz)

  integer :: i, j, k
  real(r8) :: xp(0:nxf+1,1:ny), yp(1:nx,0:nyf+1), temp, tmp, dtin

  dpx = 0.d0
  do k = 1, nz
    xp = 0.d0
    do j = 1, ny
      do i = xu_bgn-1, xu_end+1
        xp(i,j) = iu(i,j,k)*(p(i,j,k)-p(i-1,j,k))
      end do
    end do
    do j = 1, ny
      do i = 1, nx
        tmp = o12*dble(iu(i,j,k)*iu(i+1,j,k))
        dpx(i,j,k) = 0.5d0*(xp(i+1,j) + xp(i,j))  &
                   + tmp*(-xp(i-1,j) + xp(i,j)   &
                          +xp(i+1,j) - xp(i+2,j))
      end do
      dpx(:,j,k) = dpx(:,j,k)*odx(j)
    end do
  end do

  dpy = 0.d0
  do k = 1, nz
    yp = 0.d0
    do j = yv_bgn-1, yv_end+1
      do i = 1, nx
        yp(i,j) = iv(i,j,k)*(p(i,j,k)-p(i,j-1,k))
      end do
    end do
    do j = 1, ny
      do i = 1, nx
        tmp = o12*dble(iv(i,j,k)*iv(i,j+1,k))
        dpy(i,j,k) = 0.5d0*(yp(i,j+1) + yp(i,j)) &
                  +  tmp*(-yp(i,j-1) + yp(i,j)  &
                          +yp(i,j+1) - yp(i,j+2))
      end do
      dpy(:,j,k) = dpy(:,j,k)*ody(j)
    end do
  end do

  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        dtin = dt*in(i,j,k)*0.d0
        u1(i,j,k) = -dtin*dpx(i,j,k)
        v1(i,j,k) = -dtin*dpy(i,j,k)
      end do
    end do
  end do
end subroutine pre_def_velocity

subroutine pre_def_uvw
  implicit none

  integer :: i, j, k
  real(r8) :: temp, tmp

  u(1:nxf,1:ny,:) = 0.5d0*(u1(0:nx,1:ny,:) + u1(1:nx+1,1:ny,:))
  u(xu_bgn:xu_end,1:ny,1:nz) = u(xu_bgn:xu_end,1:ny,1:nz)*iu(xu_bgn:xu_end,1:ny,1:nz)

  v(1:nx,1:nyf,:) = 0.5d0*(v1(1:nx,0:ny,:) + v1(1:nx,1:ny+1,:))
  v(1:nx,yv_bgn:yv_end,1:nz) = v(1:nx,yv_bgn:yv_end,1:nz)*iv(1:nx,yv_bgn:yv_end,1:nz)

  w = 0.d0
  do j = 1, ny              !(j:cell & north face for v)
    temp = ocs(j)*ody(j)
    do i = 1, nx            !(i:cell & east  face for u)
      do k = kb(i,j), 1, -1
        tmp = 1.d0/odz(k)
        w(i,j,k) = (w(i,j,k+1) &
                 + ((u(i+1,j,k)-u(i,j,k))*odx(j) &
                 + (csv(j+1)*v(i,j+1,k)-csv(j)*v(i,j,k))*temp)*tmp)
      end do
    end do
  end do
end subroutine pre_def_uvw

subroutine go_check_conservation( var_chk, vol_chk, total_vol, msg_chk, t_avg_glb)
  implicit none

  real(r8), intent(in) :: var_chk(nx,ny,nz)
  real(r8), intent(in) :: vol_chk(nx,ny,nz), total_vol
  real(r8), intent(inout) :: t_avg_glb
  character(len=*), intent(in) :: msg_chk
  character(len=256) :: log_info
  integer :: ierr
  real(r8) :: totalt, dtmp


  totalt = sum(var_chk*vol_chk)/total_vol
  dtmp   = totalt
  call MPI_ALLREDUCE(dtmp, totalt, 1, MPI_REAL8, MPI_SUM, m_comm_cart, ierr)
  if (myid .eq. rootid) then
    write(log_info,'(a40,2f18.12)') trim(msg_chk), totalt, totalt-t_avg_glb
    call comm_write_log_info(fid_log, log_info)
  end if

  t_avg_glb = totalt
end subroutine go_check_conservation
end module tai_timcom_init
