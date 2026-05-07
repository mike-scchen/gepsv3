      subroutine inpint_cctogd(coef,grid,km)                                             
!$$$  subprogram documentation block                                            
!                .      .    .                                       .          
! subprogram:  inpint_cctogd                                                           
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 94-02-06          
!                                                                               
! abstract:  transform regional cos-cos wave coefficient                        
!            to regional grids.                                                 
!	     for regional spectral model.                                             
!                                                                               
! program history log:                                                          
!                                                                               
! usage:    call inpint_cctogd(coef,grid,km)                                             
!   input argument list:                                                        
!     coef	- wave coefficent with dimension of lnwav
!     km	- vertical dimension.                                                  
!                                                                               
!   output argument list:                                                       
!     grid     - grid field with dimension of lngrd
!                                                                               
!   input files: none                                                           
!                                                                               
!   output files: none                                                          
!                                                                               
!   subprograms called:                                                         
!	sumffc	- summation coefficient of cos wave in y direction                     
!	ffscos	- fast fourier coefficient transform of cos wave in                    
!                 x direction.                                                  
!                                                                               
!   remark: none                                                                
!                                                                               
! attributes:                                                                   
!   language: fortran 77.                                                       
!   machine:  cray c90.                                                         
!                                                                               
!$$$                                                                            
#include <inpint_paramodel.h>
#include <rscompln.h>
#define LNGRDS lngrd
#define LNWAVS lnwav
      dimension coef(LNWAVS,km),grid(LNGRDS,km)
      real, allocatable :: syn(:,:)
!
!      allocate (syn(igrd1,km))
!      call maxmin(coef,lnwav,km,1,km,' in inpint_gdtocc, coef')
      call inpint_ffscosy1(coef,grid,1,iwav1,igrd1,1,1)
      call inpint_ffscosx1(grid,jgrd1,1,1,1)
!      call maxmin(grid,lngrd,km,1,km,' in inpint_gdtocc, grid')
!      do lat=1,jgrd1                                                          
!        call sumffc1(coef,syn,ccosg(1,lat),
!     1             lnwav,iwav1,jwav1,igrd1,km)                                     
!        call inpint_ffscosx1(syn,km,1,1,1)                                                  
!      call maxmin(syn,igrd1,1,1,1,' in inpint_gdtocc, syn')
!        jlat=(lat-1)*igrd1                                                    
!        do k=1,km                                                               
!         do i=1,igrd1                                                         
!          j=i+jlat                                                              
!          grid(j,k) = syn(i,k)                                                  
!         enddo                                                                  
!        enddo                                                                   
!      enddo
!      deallocate (syn)         
!
      return                                                                    
      end                                                                       
