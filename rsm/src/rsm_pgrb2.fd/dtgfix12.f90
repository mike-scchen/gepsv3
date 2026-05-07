subroutine dtgfix12(cdtg1,cdtg2,chg)
!integer*8 idtg1,idtg2
integer   chg,yyyy,mm,dd,hh,ii,dt,diff,rem
integer   month(12)
character(len=12):: cdtg ,cdtg1 ,cdtg2
data      month/31,28,31,30,31,30,31,31,30,31,30,31/
!===start===
!write(cdtg,'(i12)')idtg1
cdtg=cdtg1
read(cdtg,'(i4,i2,i2,i2,i2)')yyyy,mm,dd,hh,ii

if( ((mod(yyyy,4  ) .eq. 0) .and. (mod(yyyy,100) .ne. 0))   .or.  &
    ((mod(yyyy,100) .eq. 0) .and. (mod(yyyy,400) .eq. 0)) ) then
    month(2)=29
endif

dt = hh + chg  
if(dt .lt. 0) then
  diff=chg*(-1)/24
  rem =mod(chg*(-1),24)
  hh  =hh-rem
  if (hh .lt. 0)then
     diff=diff+1
     hh=hh+24
  endif
  dd=dd-diff
  do while(dd .lt. 1)
    mm=mm-1
    if(mm .gt. 0) then
        dd=dd+month(mm)
    else
        yyyy=yyyy-1
        mm=mm+12
        dd=dd+month(mm)
        if(((mod(yyyy,4  ) .eq. 0) .and. (mod(yyyy,100) .ne. 0))  .or.  &
          ( (mod(yyyy,100) .eq. 0) .and. (mod(yyyy,400) .eq. 0))) then
          month(2)=29
        else
          month(2)=28
        endif
    endif
  enddo

elseif (dt .gt. 23)then
  diff=chg/24
  rem=mod(chg,24)
  hh=hh+rem
  if(hh .gt. 23)then
    hh=hh-24
    diff=diff+1
  endif
  dd=dd+diff
  do while (dd.gt.month(mm))
    dd=dd-month(mm)
    if(mm .eq. 12)then
      yyyy=yyyy+1
      mm=1
      if(((mod(yyyy,4  ) .eq. 0) .and. (mod(yyyy,100) .ne. 0))  .or.  &
         ((mod(yyyy,100) .eq. 0) .and. (mod(yyyy,400) .eq. 0))) then
        month(2)=29
      else
        month(2)=28
      endif
    else
      mm=mm+1
    endif
  enddo
else 
  hh=dt 
endif
write(cdtg,'(i4,i2.2,i2.2,i2.2,i2.2)')yyyy,mm,dd,hh,ii
! read(cdtg,'(i12)')idtg2
 cdtg2=cdtg

return

100    print *,'dtgfix12 error, idtg1=',idtg1
return
end subroutine
!
