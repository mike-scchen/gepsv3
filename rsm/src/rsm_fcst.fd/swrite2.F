      subroutine swrite2(n, phy_f3d, phy_f2d, iflag, iorestartfile)
!
! Output 3-D and 2-D microphysics variables
!
      use machine, only: kind_io4, kind_io8

      implicit none

#include <paramodel.h>

#ifdef MP
#include <npesi.h>
#include <comrmpi.h>
#define IGRD1S igrd1p
#define JGRD1S jgrd1p
#else
#define IGRD1S igrd1
#define JGRD1S jgrd1
#endif

      integer :: n
      real :: phy_f3d(IGRD1S, JGRD1S, levr, num_p3d),                  &
     &        phy_f2d(IGRD1S, JGRD1S, num_p2d)
      integer :: iflag
      logical :: iorestartfile

      integer :: i, j, k, m
      real, allocatable :: grid2d(:, :)

#ifdef MP
      if (mype == master) then
#endif
        print '(a)',     '---------- begin swrite2 ----------'
        print '(a, i0)', 'swrite2: file unit       = ', n
        print '(a, i0)', 'swrite2: iflag           = ', iflag
        print '(a, l1)', 'swrite2: iorestartfile   = ', iorestartfile
#ifdef MP
      end if
#endif

      allocate(grid2d(igrd1, jgrd1))

      do m = 1, num_p3d
      do k = 1, levr
#ifdef MP
        call rmpgp2f(phy_f3d(1, 1, k, m), igrd1p, jgrd1p,              &
     &               grid2d(1, 1), igrd1, jgrd1, 1)

        if (mype == master) then
#else
        grid2d(:, :) = phy_f3d(:, :, k, m)
#endif

        do j = 1, jgrd1
        do i = 1, igrd1
! Rationales:
! If nmmiph=2 (Zhao-Carr), output phy_f3d = -999.0
! If nmmiph=3 (Brad-Ferrier), output "f_ice" "f_rain" "f_rimef" only
! If nmmiph=4 (NEW_WSM3), output "rew" "rei" "res"
! If nmmiph=6 (NEW_WSM6), output "rew" "rei" "res"
! If nmmiph=11(GFDL MP v1), output "rew" "rei" "rer" "res" "reg"
! If nmmiph=12(GFDL MP v2), output "rew" "rei" "rer" "res" "reg"
! If nmmiph=13(GFDL MP v3), output "rew" "rei" "rer" "res" "reg"
          if (nmmiph == 2 .and. iflag == 1 .and. m /= 3) then
            grid2d(i, j) = max(0.0, min(1.0, grid2d(i, j)))
          else if (nmmiph == 3 .and. iflag == 1) then
            grid2d(i, j) = -999.0
          else
            grid2d(i, j) = max(0.0, grid2d(i, j))
          end if
        end do
        end do

        if (iflag == 0 .or. m <= 3) then
          print '(4(a, i0), 2(a, es13.6e2))',                          &
     &      'swrite2: 3-D variable ', m, ' of ', num_p3d,              &
     &      ', layer ', k, ' of ', levr,                               &
     &      ', max = ', maxval(grid2d),                                &
     &      ', min = ', minval(grid2d)

          if (iorestartfile) then
            write(n) real(grid2d, kind_io8)
          else
            write(n) real(grid2d, kind_io4)
          end if
        end if
#ifdef MP
        end if ! mype == master
#endif
      end do ! levr
      end do ! num_p3d

      do m = 1, num_p2d
#ifdef MP
        call rmpgp2f(phy_f2d(1, 1, m), igrd1p, jgrd1p,                 &
     &               grid2d(1, 1), igrd1, jgrd1, 1)

        if (mype == master) then
#else
        grid2d(:, :) = phy_f2d(:, :, m)
#endif
        if (iflag == 0) then
          print '(2(a, i0), 2(a, es13.6e2))',                          &
     &      'swrite2: 2-D variable ', m, ' of ', num_p2d,              &
     &      ', max = ', maxval(grid2d),                                &
     &      ', min = ', minval(grid2d)

          if (iorestartfile) then
            write(n) real(grid2d, kind_io8)
          else
            write(n) real(grid2d, kind_io4)
          end if
        end if
#ifdef MP
        end if ! mype == master
#endif
      end do ! num_p2d

      deallocate(grid2d)

#ifdef MP
      if (mype == master) then
#endif
        print '(a)', '---------- end swrite2 ----------'
#ifdef MP
      end if
#endif
      end subroutine swrite2
