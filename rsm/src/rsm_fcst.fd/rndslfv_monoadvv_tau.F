      subroutine rndslfv_monoadvv_tau(grid_gr,lonlens,siin,igrd1p,levr,pdot,deltim)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log:
! 2011 02 20 : henry juang, initial implemented into nems as NDSL with mass_dp
! 2013 09 30 : henry juang, add option of theta advection, (used later)
!
! 2018 08 09 : ying-ju chen, changed to RSM usage
!
      use machine, only: kind_grid => kind_io8
      implicit none
      integer, intent(in):: lonlens,igrd1p,levr
      real,   intent(in):: deltim
!
      real(kind=kind_grid) grid_gr(igrd1p,levr)
      real(kind=kind_grid) pdot(igrd1p,levr+1)
      real(kind=kind_grid) siin(levr+1)
   
      real      qqlon(igrd1p,levr) !t
!      real      rdt2

      integer mass,nvars
      integer i,k
!
      mass  = 0
      nvars = 1
!      rdt2  = 0.5/deltim
!
! t at n+1* 
        do k=1,levr
          do i=1,lonlens
            qqlon(i,k) = grid_gr(i,k) 
          enddo
        enddo
!
        mass=0
        call rvertical_cell_advect (lonlens,igrd1p,levr,nvars,          &
     &            deltim,siin,pdot(1,1),qqlon(1,1),mass)
!
        do k=1,levr
          do i=1,lonlens
! t n+1
            grid_gr(i,k) = qqlon(i,k)
! t tendency at n
!yj            grid_gr(i,k) = (qqlon(i,k)-grid_gr(i,k))*rdt2
          enddo
        enddo
!
! ===============================
!
      return
      end
