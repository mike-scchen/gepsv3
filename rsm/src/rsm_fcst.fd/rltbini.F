      subroutine rltbini
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  rltbini
!   prgmmr:  hann-ming henry juang      org: w/nmc20    date: 92-02-06
!
! abstract:  prepare constant arries for lateral boundary relaxation
!
! program history log:
!
! usage:    call rltbini
!   input argument list:
!
!   output argument list:
!
!   common block:
! rscomver  rscomltb
!
!   input files: none
!
!   output files: none
!
!   subprograms called:
! intrinsic function : asin  abs  float  amax1  amin1
!
!   remark: none
!
! attributes:
!   language: fortran 77.
!   machine:  cray ymp.
!
!$$$
#include <paramodel.h>
#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
#define RLTBS rltbse
#define SLTBS sltbse
      real, allocatable :: rltbse(:)
      real, allocatable :: sltbse(:)
!yj>      
!yj<      
#else
#define RLTBS rltb
#define SLTBS sltb
!yj>      
!yj<      
#endif
!....
!.....developed by hann-ming henry juang (nov 1990, version 0)......
!.....use u v as input, update u v tendency ....
!.....
!....
#include <rscomio.h>
#include <rscomltb.h>
!cc
!....
      save icall
!                                                         
      data icall/0/                                      
!
#ifdef MP
      allocate (rltbse(lngrd))
      allocate (sltbse(lngrd))
#endif
!
      if( icall.eq.0 ) then                            
        icall = 1
        pi = asin(1.0) * 2.0                         
!                                                   
#ifdef STPLTB
! steep slope of lateral boundary relaxation coefficients
        pow=51.0
#else
        pow=15.0
#endif
#ifdef MP
      if (mype.eq.master) then
#endif  
        print *,' in latbnd   factor=',pow
#ifdef MP
      endif
#endif
!                                             
        xc = float(igrd)/2.0 + 1.0  
        yc = float(jgrd)/2.0 + 1.0 
        xc1= xc - 1            
        yc1= yc - 1           
!
       do 10 lat = 1,jgrd1
        lats=lat
        jlats=(lats-1)*igrd1
!                             
        do 4 i=1,igrd1      
        xl = abs(i-xc)      
!
        yl = abs(lats-yc)   
!yjtest>
#ifdef NEWLTB
        grcoeout= ( min( max(xl/xc1,yl/yc1) , 1.) )
        grcoeout= exp(-0.05*(grcoeout/0.7)**25.0)
!yj        coeout = 1.0-exp(-0.05*(coeout/0.5)**8.0)
        if( grcoeout.lt.0.01 ) grcoeout=0.0
        RLTBS(i+jlats)=grcoeout
        SLTBS(i+jlats)=1.-grcoeout
!!NEWLTB2
!        grcoeout= ( min( max(xl/xc1,yl/yc1) , 1.) )
!        grcoeout= exp(-0.05*((grcoeout-0.15)/0.7)**25.0)
!!yj        coeout = 1.0-exp(-0.05*(coeout/0.5)**8.0)
!        if( grcoeout.lt.0.01 ) grcoeout=0.0
!        RLTBS(i+jlats)=grcoeout
!        SLTBS(i+jlats)=1.-grcoeout
!yj<          
!NEWLTB3
!       grcoeout= ( min( max(xl/xc1,yl/yc1) , 1.) )
!       grcoeout= exp(-0.95*((grcoeout-0.17)/0.8)**50.0)
!       if( grcoeout.lt.0.01 ) grcoeout=0.0
!       RLTBS(i+jlats)=grcoeout
!       SLTBS(i+jlats)=1.-grcoeout
!yj<
!--tanh
!        grcoeout= ( min( max(xl/xc1,yl/yc1) , 1.) )
!        grcoeout= 1.-( tanh(1.5*(1-(grcoeout)**1.5)*pi) )**2.
!        if( grcoeout.lt.0.01 ) grcoeout=0.0
!        RLTBS(i+jlats)=1.-grcoeout
!        SLTBS(i+jlats)=grcoeout
!yjtest<
#else
        coeout = ( min( max(xl/xc1,yl/yc1) , 1.) )**pow
        if( coeout.lt.0.02 ) coeout=0.0
        RLTBS(i+jlats)=1.-coeout
        SLTBS(i+jlats)=coeout
#endif
!
 4      continue
!
 10    continue
!
#ifdef MP
      call rmpgf2p(rltbse,igrd1,jgrd1,rltb,igrd1p,jgrd1p,1)
      call rmpgf2p(sltbse,igrd1,jgrd1,sltb,igrd1p,jgrd1p,1)
#endif
#ifdef MP
      if (mype.eq.master) then
#endif
            print *,' ------ initiate latbnd finished ---- '
#ifdef MP
      endif
#endif
      endif
! ------------------------------
#ifdef MP
      deallocate (rltbse)
      deallocate (sltbse)
#endif
!
      return
      end
