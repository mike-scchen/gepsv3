      module rmachine

      implicit none
      save
!  machine dependant constants
      integer kind_io4,kind_io8,kind_phys,kind_rad,kind_REAL,kind_INTEGER
      parameter (kind_rad = selected_real_kind(13,60)) ! the '60' maps to 64-bit real
      parameter (kind_phys = selected_real_kind(13,60)) ! the '60' maps to 64-bit real
      parameter (kind_io4 = 4)
      parameter (kind_io8 = 8)
      parameter (kind_REAL = 8)
      parameter (kind_INTEGER = 4)

      end module rmachine
