      subroutine sread2(n, phy_f3d, phy_f2d, iorestartfile)
!
! Input 3-D and 2-D microphysics variables
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
      logical :: iorestartfile

      integer :: k, m
      real, allocatable :: grid2d(:, :)
      real(kind_io4), allocatable :: buffersp(:)
      real(kind_io8), allocatable :: bufferdp(:)

#ifdef MP
      if (mype == master) then
#endif
        print '(a)',     '---------- begin sread2 ----------'
        print '(a, i0)', 'sread2: file unit       = ', n
        print '(a, l1)', 'sread2: iorestartfile   = ', iorestartfile
#ifdef MP
      end if
#endif

      allocate(grid2d(igrd1, jgrd1))

#ifdef MP
      if (mype == master) then
#endif
        if (iorestartfile) then
          allocate(bufferdp(size(grid2d)))
        else
          allocate(buffersp(size(grid2d)))
        end if
#ifdef MP
      end if
#endif

      do m = 1, num_p3d
      do k = 1, levr
#ifdef MP
        if (mype == master) then
#endif
          if (iorestartfile) then
            read(n) bufferdp

            grid2d(:, :) = reshape(real(bufferdp, kind(grid2d)),       &
     &                             shape(grid2d))
          else
            read(n) buffersp

            grid2d(:, :) = reshape(real(buffersp, kind(grid2d)),       &
     &                             shape(grid2d))
          end if

          print '(4(a, i0), 2(a, es13.6e2))',                          &
     &      'sread2: 3-D variable ', m, ' of ', num_p3d,               &
     &      ', layer ', k, ' of ', levr,                               &
     &      ', max = ', maxval(grid2d),                                &
     &      ', min = ', minval(grid2d)
#ifdef MP
        end if
#endif

#ifdef MP
        call rmpgf2p(grid2d(1, 1), igrd1, jgrd1,                       &
     &               phy_f3d(1, 1, k, m), igrd1p, jgrd1p, 1)
#else
        phy_f3d(:, :, k, m) = grid2d(:, :)
#endif
      end do ! levr
      end do ! num_p3d

      do m = 1, num_p2d
#ifdef MP
        if (mype == master) then
#endif
          if (iorestartfile) then
            read(n) bufferdp

            grid2d(:, :) = reshape(real(bufferdp, kind(grid2d)),       &
     &                             shape(grid2d))
          else
            read(n) buffersp

            grid2d(:, :) = reshape(real(buffersp, kind(grid2d)),       &
     &                             shape(grid2d))
          end if

          print '(2(a, i0), 2(a, es13.6e2))',                          &
     &      'sread2: 2-D variable ', m, ' of ', num_p2d,               &
     &      ', max = ', maxval(grid2d),                                &
     &      ', min = ', minval(grid2d)
#ifdef MP
        end if
#endif

#ifdef MP
        call rmpgf2p(grid2d(1, 1), igrd1, jgrd1,                       &
     &               phy_f2d(1, 1, m), igrd1p, jgrd1p, 1)
#else
        phy_f2d(:, :, m) = grid2d(:, :)
#endif
      end do ! num_p2d

#ifdef MP
      if (mype == master) then
#endif
        if (iorestartfile) then
          deallocate(bufferdp)
        else
          deallocate(buffersp)
        end if
#ifdef MP
      end if
#endif

      deallocate(grid2d)

#ifdef MP
      if (mype == master) then
#endif
        print '(a)', '---------- end sread2 ----------'
#ifdef MP
      end if
#endif
      end subroutine sread2
