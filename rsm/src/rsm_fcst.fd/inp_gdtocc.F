      subroutine inp_gdtocc(grid,coef,km)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  inp_gdtocc
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  do regional grid to cos-cos wave transformation.
!
! program history log:
!
! usage:    call  inp_gdtocc(grid,coef,km)
!   input argument list:
!     grid  - grid field with dimension of (lngrd,km)
!     km  - the second dimension of grid and ff
!
!   output argument list:
!     coef     - wave coefficent with dimension of (lnwav,km)
!
!   common block:
! repro  - a temperory common block.
!
!   input files: none
!
!   output files: none
!
!   subprograms called:
! ffacos  sumggc
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
#include <inp_paramodel.h>
#include <rscompln.h>
#define LNGRDS lngrd
#define LNWAVS lnwav
      dimension grid(LNGRDS,km),coef(LNWAVS,km)
      dimension tmp(lngrd,km),tmpsp(lnwav,km)
      real, allocatable :: syn(:,:)
!
      allocate (syn(igrd1,km))
      syn=0.0
      do k=1,km
       do i=1,lnwav
        coef(i,k) = 0.0
       enddo
      enddo
!      print *,' in inp_gdtocc,j=',jgrd1,'k=',km
!      call maxmin(grid,lngrd,km,1,km,' in inp_gdtocc, grid')
      call inp_ffacosx1(grid, jgrd1,1,1,1)
      call inp_ffacosy1(grid,coef,1,igrd1,iwav1,1,1)
!      do j=1,jgrd1
!        jlat=(j-1)*igrd1
!        do k=1,km
!         do i=1,igrd1
!          ij=i+jlat
!          syn(i,k) = grid(ij,k)
!         enddo
!        enddo
!        call inp_ffacosx1(syn, km,1,1,1)
!        call sumggc1(syn,coef,gcosc(1,j),
!     1             lnwav,iwav1,jwav1,igrd1,km)
!      enddo
!      call maxmin(coef,lnwav,km,1,km,' in inp_gdtocc, coef')
      deallocate (syn)
      return
      end
