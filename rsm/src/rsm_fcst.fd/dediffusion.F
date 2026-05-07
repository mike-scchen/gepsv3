      subroutine dediffusion(nx,ui,dx,tau)
!      
! this is de-diffusion for adding short waves enery back     
!
! vars:
! nx  dimension of input data (ui) and output data (uo)
! ui  input data, data used to do diffusion
!     output data, data after de-diffusion (= ui + coef * ui_xxxx)

      integer,intent(in):: nx
      real, intent(inout)  :: ui(nx)
      real,intent(in):: dx,tau
!local var
      real :: taui, coef
      real :: tmp(nx)
      integer :: i

      taui=1.0/tau
      ! fourth order diffusion
!     coef=0.5*taui*(dx/2.)**4.0
!test coef
      coef=16.0*0.5*taui*(dx/2.)**4.0
!     coef=1.5e14

! calculate 4th order differential
      do i=1,nx
        if (i.ge.3 .and. i.le.nx-2) then
          !central finite difference
          call cfd4(ui(i-2),ui(i-1),ui(i),ui(i+1),ui(i+2),          &
                   tmp(i),dx)
        elseif (i.lt.3) then
          !forward finite difference
          call ffd4(ui(i),ui(i+1),ui(i+2),ui(i+3),ui(i+4),ui(i+5),  &
                   tmp(i),dx)
        else !(i.gt.nx-2)
          !backward finite difference
          call bfd4(ui(i-5),ui(i-4),ui(i-3),ui(i-2),ui(i-1),ui(i),  &
                   tmp(i),dx)
        endif
      enddo


! add diffusion to ui
      do i=1,nx
        ui(i) = ui(i) + coef*tmp(i)
      enddo

      return
      end subroutine dediffusion

!-----central finite difference-----
      subroutine cfd4(m2,m1,p0,p1,p2,df,dx)
      real, intent(in)  :: m2,m1,p0,p1,p2,dx
      real, intent(out) :: df

      df=(1./(dx**4.))*(m2-4.*m1+6.*p0-4.*p1+p2)

      return
      end subroutine cfd4
! 
!-----forward finite difference-----
      subroutine ffd4(p0,p1,p2,p3,p4,p5,df,dx)
      real, intent(in)  :: p0,p1,p2,p3,p4,p5,dx
      real, intent(out) :: df

      df=(1./(dx**4.))*(3.*p0-14.*p1+26.*p2-24.*p3+11.*p4-2.*p5)

      return
      end subroutine ffd4
!       
!-----forward finite difference-----
      subroutine bfd4(m5,m4,m3,m2,m1,p0,df,dx)
      real, intent(in)  :: p0,m1,m2,m3,m4,m5,dx
      real, intent(out) :: df

      df=(1./(dx**4.))*(3.*p0-14.*m1+26.*m2-24.*m3+11.*m4-2.*m5)

      return
      end subroutine bfd4
