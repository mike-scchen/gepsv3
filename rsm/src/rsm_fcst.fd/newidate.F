      subroutine newidate(idate,ph,idateo)
! INPUT
!  idate="(I4,I2,I2,I2)" iy,im,id,ih    
!  iy : initial year
!  im : initial month
!  id : initial day
!  ih : initial hour
!  ph : add hour
!
! OUTPUT
!idate has the form idate(1):hh,idate(2):mm,idate(3):dd,idate(4):yy
!  idateo(4)=yy : new year
!  idateo(2)=mm : new month
!  idateo(3)=dd : new day
!  idateo(1)=hh : new hour

      integer idate
      character*10 cdate
      integer,dimension(4):: idateo
      integer iy,im,id,ih,yy,mm,dd,hh
      real ph
      integer md(12)
      data md/31,28,31,30,31,30,31,31,30,31,30,31/

      write(cdate,"(I10)") idate
      read(cdate,"(I4,I2,I2,I2)") iy,im,id,ih

!      logical more

      if (mod(iy,4) .eq. 0) then
       md(2)=29
      else
       md(2)=28
      endif

      hh=ih+ph

      if (hh .ge. 24) then
      dd=id+hh/24
      hh=mod(hh,24)
!       more=.true.
 800 continue
!      do while (more)
        if (dd .gt. md(im)) then
        dd=dd-md(im)
        mm=im+1
          if (mm .gt. 12) then
          yy=iy+1
          mm=mm-12
          else
          yy=iy
          endif
          if (mod(yy,4) .eq. 0) md(2)=29
        else
        mm=im
        yy=iy
        endif
        if (dd .gt. md(mm)) goto 800
!         if (dd .le. md(im)) more=.false.
!      end do

      elseif (hh .lt. 0) then
       dd=id+(hh/24)-1
       hh=24+mod(hh,24)
!       more=.true.
 900 continue
!      do while (more)
        if (dd .le. 0) then
        mm=im-1
          if (mm .le. 0) then
          yy=iy-1
          mm=mm+12
          else
          yy=iy
          endif
          if (mod(iy,4) .eq. 0) md(2)=29
        dd=md(mm)+dd
        else
        mm=im
        yy=iy
        endif
          if (dd .gt. md(mm)) goto 900
!         if (dd .le. md(im)) more=.false.
!       end do

      else
        dd=id
        mm=im
        yy=iy

      endif
!idateo has the form idateo(1):hh,idateo(2):mm,idateo(3):dd,idateo(4):yy
      idateo(1)=hh
      idateo(2)=mm
      idateo(3)=dd
      idateo(4)=yy

      return
      end
