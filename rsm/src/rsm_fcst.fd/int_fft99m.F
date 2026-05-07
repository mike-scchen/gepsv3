!IBM-TWN   subroutine fft991(a,work,trigs,ifax,inc,jump,n,m,isign)
!ch   subroutine rfftmlt(a,work,trigs,ifax,inc,jump,n,m,isign)
      subroutine int_fft99m(a,work,trigs,ifax,inc,jump,n,m,isign)

!     *******************************************************
!     *****   Sept 1995 - Zaphiris Christidis, IBM Res  *****
!     *******************************************************
!     *******************************************************
!     ****   Modified by Zhan Zhang, Met. FSU, Aug. 1996  ***
!     *******************************************************
!
!     *******************************************************
!     ****   Modified by IBM-TWN for CWB Aug. 2011        ***
!     *******************************************************
!

      Real*8 a,trigs,work,scale
      DIMENSION A(*),WORK(N+2,M),TRIGS(*),IFAX(1)
!
!  ESSL fft wrapper routine.
!  Common containing fourier transform paramters for ESSL routines
!
!        parameter (naux1=2156,naux2=772)
!CWB2014
!       parameter (naux1=22000,naux2=20000)
        parameter (naux1=44002,naux2=20000)
        real*8 aux1(naux1),aux2(naux2),aux3(1)
        common/ibmfft_init/aux1,aux2,aux3,ibmn,ibmm,ibmtype
!CWB    save /ibmfft_init/
!  2012-11-02 >>>
!$omp threadprivate(/ibmfft_init/)
!  2012-11-02 <<<
!
!--   Stride (jump) is a even number .......
!
      if(mod(jump,2).eq.0) then
         jump2=jump/2
         if(isign.eq.1) then
            scale=1.0
            ibmsign=-1
            if(isign.ne.ibmtype.or.m.ne.ibmm.or.n.ne.ibmn) then
!
!     compute twiddle factors and save them. Never call it again.
!
               call dcrft(1,a,jump2,a,jump,n,m,ibmsign,scale,   &
                          aux1,naux1,aux2,naux2,aux3,0)
               ibmn=n
               ibmm=m
               ibmtype=isign
            end if
            call dcrft(0,a,jump2,a,jump,n,m,ibmsign,scale,      &
                       aux1,naux1,aux2,naux2,aux3,0)
!
!--    Initialize data for consistency with fft991
!
            ib=n*inc+1
            do i=1,m
               a(ib)=0.0
               a(ib+inc)=0.0
               ib=ib+jump
            end do
         else
            scale=1.0/dfloat(n)
            ibmsign=1
            if(isign.ne.ibmtype.or.m.ne.ibmm.or.n.ne.ibmn) then
!
!     compute twiddle factors and save them. Never call it again.
!
               call drcft(1,a,jump,a,jump2,n,m,ibmsign,scale,    &
                          aux1,naux1,aux2,naux2,aux3,0)
               ibmn=n
               ibmm=m
               ibmtype=isign
            end if

            call drcft(0,a,jump,a,jump2,n,m,ibmsign,scale,       &
                       aux1,naux1,aux2,naux2,aux3,0)
         end if
!
!--   Stride (jump) is an odd number .......
!
      else
         jump2=(n/2)+1

         if(isign.eq.1) then
            scale=1.0
            ibmsign=-1
!
!---  Swap the data from a to work, while jump is an odd number
!
            i2=1
            do i=1,m
               do j=1,n+2
                  work(j,i)=a(i2)
                  i2=i2+1
               end do
               i2=jump*i+1
            end do
!
            if(isign.ne.ibmtype.or.m.ne.ibmm.or.n.ne.ibmn) then
!
!     compute twiddle factors and save them. Never call it again.
!
               call dcrft(1,work,jump2,a,jump,n,m,ibmsign,scale,    &
                          aux1,naux1,aux2,naux2,aux3,0)
               ibmn=n
               ibmm=m
               ibmtype=isign
            end if

            call dcrft(0,work,jump2,a,jump,n,m,ibmsign,scale,       &
                       aux1,naux1,aux2,naux2,aux3,0)

         else
            scale=1.0/dfloat(n)
            ibmsign=+1
            if(isign.ne.ibmtype.or.m.ne.ibmm.or.n.ne.ibmn) then
!
!     compute twiddle factors and save them. Never call it again.
!
               call drcft(1,a,jump,work,jump2,n,m,ibmsign,scale,    &
                          aux1,naux1,aux2,naux2,aux3,0)
               ibmn=n
               ibmm=m
               ibmtype=isign
            end if

            call drcft(0,a,jump,work,jump2,n,m,ibmsign,scale,       &
                       aux1,naux1,aux2,naux2,aux3,0)
!
!---  Swap the data from work to a
!
            i2=1
            do i=1,m
               do j=1,n+2
                  a(i2)=work(j,i)
                  i2=i2+1
               end do
               i2=jump*i+1
            end do
         end if
      end if

      return
      end
! =====================================================================
