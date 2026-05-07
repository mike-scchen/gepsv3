      subroutine read_idate(ifile)
!
      use module_gfsdata, only: idate
!
      rewind(ifile)
      read(ifile) idate

      print*,'read gfs idate=',idate
      close(ifile)
!
      return
      end
