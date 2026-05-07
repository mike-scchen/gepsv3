module timcom_output
  implicit none

  type :: output_dim_panel
    character(len=8) :: dname(4)
    integer(MPI_OFFSET_KIND) :: ndim(4)
    character(len=64) :: long_name(4), units(4)
  end type output_dim_panel

  type :: output_head_panel
    character(len=64) :: vname
  end type output_head_panel
contains

subroutine def_output_var(varname)
  implicit none

  select case(trim(varname))
  case("")
  case("")
  case default
    continue
  end select

end subroutine def_output_var

end module timcom_output

