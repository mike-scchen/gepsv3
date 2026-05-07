      module mod_rndsl
      implicit none
      save

      integer  lonfull,lonpart,lonlenmax,mylonlen
      integer  latfull,latpart,latlenmax,mylatlen
!      integer  ndslhvar,ndslvvar
      integer, allocatable :: lonstrndsl(:),lonlenndsl(:)
      integer, allocatable :: latstrndsl(:),latlenndsl(:)

      end module mod_rndsl
