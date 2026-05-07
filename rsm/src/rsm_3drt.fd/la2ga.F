      subroutine la2ga(regin,imxin,jmxin,rinlon,rinlat,rlon,rlat,inttyp,         &
     &                 gauout,len,lmask,rslmsk,slmask                            &
     &,                outlat, outlon)
      use machine , only : kind_io8,kind_io4
      implicit none
!      integer kind_io8,kind_io4
!      parameter (kind_io8=8,kind_io4=4)
      real (kind=kind_io8) wei4,wei3,wei2,sum2,sum1,sum3,wei1,sum4,              &
     &                     wsum,tem,wsumiv,sums,sumn,wi2j2,x,y,wi1j1,            &
     &                     wi1j2,wi2j1,rlat,rlon,aphi,                           &
     &                     rnume,alamd,denom
      integer jy,ifill,ix,len,inttyp,me,i,j,jmxin,imxin,jq,jx,j1,j2,             &
     &        ii,i1,i2,kmami,it
      integer nx,kxs,kxt,imxnx
!
!  interpolation from lat/lon or gaussian grid to other lat/lon grid
!
      real (kind=kind_io8) outlon(len),outlat(len),gauout(len),                  &
     &                     slmask(len)
      real (kind=kind_io8) regin (imxin,jmxin),rslmsk(imxin,jmxin)
!
      real (kind=kind_io8)    rinlat(jmxin),  rinlon(imxin)
      integer iindx1(len),    iindx2(len)
      integer jindx1(len),    jindx2(len)
      real (kind=kind_io8)    ddx(len),       ddy(len),   wrk(len)
!
      logical lmask
!
      logical first
      integer   num_threads
      data first /.true./
      save num_threads, first
!
      integer len_thread_m, len_thread, i1_t, i2_t,ttt
!
      if (first) then
!#ifdef IBMSP
!         num_threads    = num_parthds()
!#else
         num_threads    = 1
!#endif
         first = .false.
      endif

      ifill = 0
!
!
!     if(me .eq. 0) then
!     print *,'rlon=',rlon,' me=',me
!     print *,'rlat=',rlat,' me=',me,' imxin=',imxin,' jmxin=',jmxin
!     endif
!
!     do j=1,jmxin
!       if(rlat.gt.0.) then
!         rinlat(j) = rlat - float(j-1)*dlain
!       else
!         rinlat(j) = rlat + float(j-1)*dlain
!       endif
!     enddo
!
!     if (me .eq. 0) then
!       print *,'rinlat='
!       print *,(rinlat(j),j=1,jmxin)
!       print *,'rinlon='
!       print *,(rinlon(i),i=1,imxin)
!
!       print *,'outlat='
!       print *,(outlat(j),j=1,len)
!       print *,(outlon(j),j=1,len)
!     endif
!
!     do i=1,imxin
!       rinlon(i) = rlon + float(i-1)*dloin
!     enddo
!
!     print *,'rinlon='
!     print *,(rinlon(i),i=1,imxin)
!
      len_thread_m  = (len+num_threads-1) / num_threads
!
!!$omp parallel do private(i1_t,i2_t,len_thread,it,i,ii,i1,i2)
!!$omp+private(j,j1,j2,jq,ix,jy,nx,kxs,kxt,imxnx,kmami)
!!$omp+private(alamd,denom,rnume,aphi,x,y,wsum,wsumiv,sum1,sum2)
!!$omp+private(sum3,sum4,wi1j1,wi2j1,wi1j2,wi2j2,wei1,wei2,wei3,wei4)
!!$omp+private(sumn,sums)
!!$omp+shared(imxin,jmxin,ifill)
!!$omp+shared(outlon,outlat,wrk,iindx1,rinlon,jindx1,rinlat,ddx,ddy)
!!$omp+shared(rlon,rlat,regin,gauout)
!
      do it=1,num_threads   ! start of threaded loop ...................
        i1_t       = (it-1)*len_thread_m+1
        i2_t       = min(i1_t+len_thread_m-1,len)
        len_thread = i2_t-i1_t+1
!
!       find i-index for interpolation
!
        do i=i1_t, i2_t
          alamd = outlon(i)
          if (alamd .lt. rlon)   alamd = alamd + 360.0
          if (alamd .gt. 360.0+rlon) alamd = alamd - 360.0
          wrk(i)    = alamd
          iindx1(i) = imxin
        enddo
        do i=i1_t,i2_t
          do ii=1,imxin
            if(wrk(i) .ge. rinlon(ii)) iindx1(i) = ii
          enddo
        enddo
        do i=i1_t,i2_t
          i1 = iindx1(i)
          if (i1 .lt. 1) i1 = imxin
          i2 = i1 + 1
          if (i2 .gt. imxin) i2 = 1
          iindx1(i) = i1
          iindx2(i) = i2
          denom     = rinlon(i2) - rinlon(i1)
          if(denom.lt.0.) denom = denom + 360.
          rnume = wrk(i) - rinlon(i1)
          if(rnume.lt.0.) rnume = rnume + 360.
          ddx(i) = rnume / denom
        enddo
!
!  find j-index for interplation
!
        if(rlat.gt.0.) then
          do j=i1_t,i2_t
            jindx1(j)=0
          enddo
          do jx=1,jmxin
            do j=i1_t,i2_t
              if(outlat(j).le.rinlat(jx)) jindx1(j) = jx
            enddo
          enddo
          do j=i1_t,i2_t
            jq = jindx1(j)
            aphi=outlat(j)
            if(jq.ge.1 .and. jq .lt. jmxin) then
              j2=jq+1
              j1=jq
             ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
            elseif (jq .eq. 0) then
              j2=1
              j1=1
              if(abs(90.-rinlat(j1)).gt.0.001) then
                ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            else
              j2=jmxin
              j1=jmxin
              if(abs(-90.-rinlat(j1)).gt.0.001) then
                ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            endif
            jindx1(j)=j1
            jindx2(j)=j2
          enddo
        else
          do j=i1_t,i2_t
            jindx1(j) = jmxin+1
          enddo
          do jx=jmxin,1,-1
            do j=i1_t,i2_t
              if(outlat(j).le.rinlat(jx)) jindx1(j) = jx
            enddo
          enddo
          do j=i1_t,i2_t
            jq = jindx1(j)
            aphi=outlat(j)
            if(jq.gt.1 .and. jq .le. jmxin) then
              j2=jq
              j1=jq-1
              ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
            elseif (jq .eq. 1) then
              j2=1
              j1=1
              if(abs(-90.-rinlat(j1)).gt.0.001) then
                 ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
              else
                 ddy(j)=0.0
              endif
            else
              j2=jmxin
              j1=jmxin
              if(abs(90.-rinlat(j1)).gt.0.001) then
                 ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
              else
                 ddy(j)=0.0
              endif
            endif
            jindx1(j)=j1
            jindx2(j)=j2
          enddo
        endif
!
!     if (me .eq. 0 .and. inttyp .eq. 1) then
!       print *,'la2ga'
!       print *,'iindx1'
!       print *,(iindx1(n),n=1,len)
!       print *,'iindx2'
!       print *,(iindx2(n),n=1,len)
!       print *,'jindx1'
!       print *,(jindx1(n),n=1,len)
!       print *,'jindx2'
!       print *,(jindx2(n),n=1,len)
!       print *,'ddy'
!       print *,(ddy(n),n=1,len)
!       print *,'ddx'
!       print *,(ddx(n),n=1,len)
!     endif
!
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        if (lmask) then
          wei1 = 0.
          wei2 = 0.
          wei3 = 0.
          wei4 = 0.
          do i=1,imxin
            sum1 = sum1 + regin(i,1) * rslmsk(i,1)
            sum2 = sum2 + regin(i,jmxin) * rslmsk(i,jmxin)
            wei1 = wei1 + rslmsk(i,1)
            wei2 = wei2 + rslmsk(i,jmxin)
!
            sum3 = sum3 + regin(i,1) * (1.0-rslmsk(i,1))
            sum4 = sum4 + regin(i,jmxin) * (1.0-rslmsk(i,jmxin))
            wei3 = wei3 + (1.0-rslmsk(i,1))
            wei4 = wei4 + (1.0-rslmsk(i,jmxin))
          enddo
!
          if(wei1.gt.0.) then
            sum1 = sum1 / wei1
          else
            sum1 = 0.
          endif
          if(wei2.gt.0.) then
            sum2 = sum2 / wei2
          else
            sum2 = 0.
          endif
          if(wei3.gt.0.) then
            sum3 = sum3 / wei3
          else
            sum3 = 0.
          endif
          if(wei4.gt.0.) then
            sum4 = sum4 / wei4
          else
            sum4 = 0.
          endif
        else
          do i=1,imxin
            sum1 = sum1 + regin(i,1)
            sum2 = sum2 + regin(i,jmxin)
          enddo
          sum1 = sum1 / imxin
          sum2 = sum2 / imxin
          sum3 = sum1
          sum4 = sum2
        endif
!
!     print *,' sum1=',sum1,' sum2=',sum2
!    *,' sum3=',sum3,' sum4=',sum4
!     print *,' rslmsk=',(rslmsk(i,1),i=1,imxin)
!     print *,' slmask=',(slmask(i),i=1,imxout)
!    *,' j1=',jindx1(1),' j2=',jindx2(1)
!
!
!  inttyp=1  take the closest point value
!
        if(inttyp.eq.1) then
        print*,'do closest point interpolation'

          do i=i1_t,i2_t
            jy = jindx1(i)
            if(ddy(i) .ge. 0.5) jy = jindx2(i)
            ix = iindx1(i)
            if(ddx(i) .ge. 0.5) ix = iindx2(i)
!
!cggg start
!
            if (.not. lmask) then

              gauout(i) = regin(ix,jy)

            else

              if(slmask(i).eq.rslmsk(ix,jy)) then

                gauout(i) = regin(ix,jy)

              else

                i1 = ix
                j1 = jy

! spiral around until matching mask is found.
                ttt=jmxin*imxin/2
                do nx=1,ttt
!                do nx=1,jmxin*imxin/2
                  kxs=sqrt(4*nx-2.5)
                  kxt=nx-int(kxs**2/4+1)
                  select case(mod(kxs,4))
                  case(1)
                    ix=i1-kxs/4+kxt
                    jx=j1-kxs/4
                  case(2)
                    ix=i1+1+kxs/4
                    jx=j1-kxs/4+kxt
                  case(3)
                    ix=i1+1+kxs/4-kxt
                    jx=j1+1+kxs/4
                  case default
                    ix=i1-kxs/4
                    jx=j1+kxs/4-kxt
                  end select
                  if(jx.lt.1) then
                    ix=ix+imxin/2
                    jx=2-jx
                  elseif(jx.gt.jmxin) then
                    ix=ix+imxin/2
                    jx=2*jmxin-jx
                  endif
                  ix=modulo(ix-1,imxin)+1
                  if(slmask(i).eq.rslmsk(ix,jx)) then
                    gauout(i) = regin(ix,jx)
                    go to 81
                  endif
                enddo

!cggg here, set the gauout value to be 0, and let's sarah's land
!cggg routine assign a default.

!             print*,'no matching mask found ',i,i1,j1,ix,jx
!             print*,'set to default value.'
              gauout(i) = 0.0


   81  continue

              end if

            end if

!cggg end

          enddo
          kmami=1
!          if (me .eq. 0) call maxmin(gauout(i1_t),len_thread,kmami,         &
!     &                        1,1,'gauout')
          cycle
        endif  ! nearest neighbor interpolation

!
!  quasi-bilinear interpolation
!
        ifill = 0
        imxnx = 0
        do i=i1_t,i2_t
          y  = ddy(i)
          j1 = jindx1(i)
          j2 = jindx2(i)
          x  = ddx(i)
          i1 = iindx1(i)
          i2 = iindx2(i)
!
          wi1j1 = (1.-x) * (1.-y)
          wi2j1 =     x  *( 1.-y)
          wi1j2 = (1.-x) *      y
          wi2j2 =     x  *      y
!
          tem = 4.*slmask(i) - rslmsk(i1,j1) - rslmsk(i2,j1)                 &
     &                       - rslmsk(i1,j2) - rslmsk(i2,j2)
          if(lmask .and. abs(tem) .gt. 0.01) then
            if(slmask(i).eq.1.) then
                wi1j1 = wi1j1 * rslmsk(i1,j1)
                wi2j1 = wi2j1 * rslmsk(i2,j1)
                wi1j2 = wi1j2 * rslmsk(i1,j2)
                wi2j2 = wi2j2 * rslmsk(i2,j2)
            else
                wi1j1 = wi1j1 * (1.0-rslmsk(i1,j1))
                wi2j1 = wi2j1 * (1.0-rslmsk(i2,j1))
                wi1j2 = wi1j2 * (1.0-rslmsk(i1,j2))
                wi2j2 = wi2j2 * (1.0-rslmsk(i2,j2))
            endif
          endif
!
          wsum   = wi1j1 + wi2j1 + wi1j2 + wi2j2
          wrk(i) = wsum
          if(wsum.ne.0.) then
            wsumiv = 1./wsum
!
            if(j1.ne.j2) then
              gauout(i) = (wi1j1*regin(i1,j1) + wi2j1*regin(i2,j1) +          &
     &                     wi1j2*regin(i1,j2) + wi2j2*regin(i2,j2))           &
     &                  *wsumiv
            else
!
              if (rlat .gt. 0.0) then
                if (slmask(i) .eq. 1.0) then
                  sumn = sum1
                  sums = sum2
                else
                  sumn = sum3
                  sums = sum4
                endif
                if( j1 .eq. 1) then
                  gauout(i) = (wi1j1*sumn        +wi2j1*sumn        +          &
     &                         wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))          &
     &                      * wsumiv
                elseif (j1 .eq. jmxin) then
                  gauout(i) = (wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+          &
     &                         wi1j2*sums        +wi2j2*sums        )          &
     &                      * wsumiv
                endif
!     print *,' slmask=',slmask(i),' sums=',sums,' sumn=',sumn
!    &,' regin=',regin(i1,j2),regin(i2,j2),' j1=',j1,' j2=',j2
!    &,' wij=',wi1j1, wi2j1, wi1j2, wi2j2,wsumiv
              else
                if (slmask(i) .eq. 1.0) then
                  sums = sum1
                  sumn = sum2
                else
                  sums = sum3
                  sumn = sum4
                endif
                if( j1 .eq. 1) then
                  gauout(i) = (wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+          &
     &                         wi1j2*sums        +wi2j2*sums        )          &
     &                      * wsumiv
                elseif (j1 .eq. jmxin) then
                  gauout(i) = (wi1j1*sumn        +wi2j1*sumn        +          &
     &                         wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))          &
     &                      * wsumiv
                endif
              endif
            endif            ! if j1 .ne. j2
          endif
        enddo
        do i=i1_t,i2_t
          j1 = jindx1(i)
          j2 = jindx2(i)
          i1 = iindx1(i)
          i2 = iindx2(i)
          if(wrk(i) .eq. 0.0) then
            if(.not.lmask) then
              write(6,*) ' la2ga called with lmask=.true. but bad',            &
     &                   ' rslmsk or slmask given'
              call abort
            endif
            ifill = ifill + 1
            if(ifill.le.2) then
              if (me .eq. 0) then
                write(6,*) 'i1,i2,j1,j2=',i1,i2,j1,j2
                write(6,*) 'rslmsk=',rslmsk(i1,j1),rslmsk(i1,j2),               &
     &                               rslmsk(i2,j1),rslmsk(i2,j2)
!               write(6,*) 'i,j=',i,j,' slmask(i)=',slmask(i)
                write(6,*) 'i=',i,' slmask(i)=',slmask(i)                       &
     &,         ' outlon=',outlon(i),' outlat=',outlat(i)
              endif
            endif
! spiral around until matching mask is found.
            do nx=1,jmxin*imxin/2
              kxs=sqrt(4*nx-2.5)
              kxt=nx-int(kxs**2/4+1)
              select case(mod(kxs,4))
              case(1)
                ix=i1-kxs/4+kxt
                jx=j1-kxs/4
              case(2)
                ix=i1+1+kxs/4
                jx=j1-kxs/4+kxt
              case(3)
                ix=i1+1+kxs/4-kxt
                jx=j1+1+kxs/4
              case default
                ix=i1-kxs/4
                jx=j1+kxs/4-kxt
              end select
              if(jx.lt.1) then
                ix=ix+imxin/2
                jx=2-jx
              elseif(jx.gt.jmxin) then
                ix=ix+imxin/2
                jx=2*jmxin-jx
              endif
              ix=modulo(ix-1,imxin)+1
              if(slmask(i).eq.rslmsk(ix,jx)) then
                gauout(i) = regin(ix,jx)
                imxnx=max(imxnx,nx)
                go to 71
              endif
            enddo
!
            write(6,*) ' error!!! no filling value found in la2ga'
            call abort
!
   71       continue
          endif
!
        enddo
      enddo            ! end of threaded loop ...................
!!$omp end parallel do
!
      if(ifill.gt.1) then
        if (me .eq. 0) then
        write(6,*) ' unable to interpolate.  filled with nearest',               &
     &             ' point value at ',ifill,' points  imxnx=',imxnx
        endif
      endif
!
      kmami=1
!      if (me .eq. 0) call maxmin(gauout,len,kmami,1,1,'gauout')
!
      return
      end subroutine la2ga
