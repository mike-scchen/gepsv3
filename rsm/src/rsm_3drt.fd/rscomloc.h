      real                                                               &
     & rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd	 &
     &,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd
      common /comrloc/                                                   &  
     & rproj,rtruth,rorient,rdelx,rdely,rcenlat,rcenlon,rlftgrd,rbtmgrd	 &
     &,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd
!      
      real TS_LAT, TE_LAT, TS_LON, TE_LON
      character :: condir*200, FILE_NAME*200
      common /cominfo/ TS_LAT, TE_LAT, TS_LON, TE_LON, condir, FILE_NAME 
!      
