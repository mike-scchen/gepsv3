      subroutine recv_idate(root_gfs,itag)
!
      use module_gfsdata, only: idate
      integer root_gfs,itag
!

        itag=1
        print*,'recv gfs idate itag=',itag, root_gfs
        call mpmd_recv(idate,1,root_gfs,itag,'I')
        

        print*,'recv gfs idate =',idate
!
      return
      end
