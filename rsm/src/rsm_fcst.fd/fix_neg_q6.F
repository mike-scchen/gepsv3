!#define compensate_q_heat
! =======================================================================
! fix negative water species
! this is designed for 6 - class micro - physics schemes
! code from GFDL MP v2
! =======================================================================

subroutine fix_neg_q6 (ix, im, km, kh, ncw, nrw, niw, nsw, ngl, &
                       dp, tt, qt, qmin_adj)
    
    use machine , only : kind_phys
    use physcons, only : con_rd, con_cp, con_cv, con_cvap, con_cliq, con_csol, &
                         con_hvap, con_rv

    integer, intent (in) :: ix, im, km, kh
    integer, intent (in) :: ncw, nrw, niw, nsw, ngl
    real (kind = kind_phys), intent (in) :: qmin_adj
    real (kind = kind_phys), intent (in), dimension (ix,km) :: dp
    real (kind = kind_phys), intent (inout), dimension (ix,km) :: tt
    real (kind = kind_phys), intent (inout), dimension (ix,kh) :: qt
    
    real (kind = kind_phys), dimension (km) :: qv, ql, qr, qi, qs, qg
    real (kind = kind_phys), dimension (km) :: lcpk, icpk

    real (kind = kind_phys) :: c_vap, c_air, c_liq, c_ice
    real (kind = kind_phys) :: c1_vap, c1_liq, c1_ice, d1_vap, d1_ice
    real (kind = kind_phys) :: lv00, li00
    real (kind = kind_phys) :: tice
    
    real (kind = kind_phys) :: dq, cvm

    integer :: i, k

    ! -----------------------------------------------------------------------
    ! define parameters
    ! -----------------------------------------------------------------------
    c_vap  = 3.0 * con_rv  ! =1384.5, heat capacity of water vapor at constant volume
    c_air  = con_cv        ! =717.6, heat capacity of dry air at constant volume
    c_liq  = 4.218e3       ! heat capacity of water at 0 deg c
!    c_liq = 4219.9         ! at 0.01 deg c
!    c_liq = 4195.5         ! at 10 deg c
!    c_liq = 4184.4         ! at 20 deg c
!    c_liq = 4180.1         ! at 30 deg c
!    c_liq = 4179.6         ! at 40 deg c
!    c_liq  = con_cliq      ! =4185.5
    c_ice  = con_csol      ! heat capacity of ice at 0 deg c
    t_ice = 273.16         ! freezing temperature

    c1_vap = c_vap / c_air
    c1_liq = c_liq / c_air
    c1_ice = c_ice / c_air

    d1_vap = ( c_vap - c_liq ) / c_air
    d1_ice = ( c_liq - c_ice ) / c_air

    lv00 = ( con_hvap - ( c_vap - c_liq ) * t_ice ) / c_air
    li00 = ( con_hvap - ( c_liq - c_ice ) * t_ice ) / c_air
    
    do i = 1, im

    ! -----------------------------------------------------------------------
    ! split six hydrometers
    ! -----------------------------------------------------------------------
    do k = 1, km
        qv (k) = qt (i,    k)
        ql (k) = qt (i,ncw+k)
        qr (k) = qt (i,nrw+k)
        qi (k) = qt (i,niw+k)
        qs (k) = qt (i,nsw+k)
        qg (k) = qt (i,ngl+k)
    enddo

    ! -----------------------------------------------------------------------
    ! define heat capacity and latent heat coefficient
    ! -----------------------------------------------------------------------
    
    do k = 1, km
        cvm = 1. + qv (k) * c1_vap + (qr (k) + ql (k)) * c1_liq + (qi (k) + qs (k) + qg (k)) * c1_ice
        lcpk (k) = (lv00 + d1_vap * tt (i,k)) / cvm
        icpk (k) = (li00 + d1_ice * tt (i,k)) / cvm
    enddo

    do k = 1, km
        
        ! -----------------------------------------------------------------------
        ! ice phase:
        ! -----------------------------------------------------------------------
        
        ! if cloud ice < qmin_adj, borrow from snow
        if (qi (k) < qmin_adj) then
            dq = qmin_adj - qi (k)
#ifdef compensate_q_heat
            qs (k) = qs (k) - dq
#endif
            qi (k) = qmin_adj
        endif
        ! if snow < qmin_adj, borrow from graupel
        if (qs (k) < qmin_adj) then
            dq = qmin_adj - qs (k)
#ifdef compensate_q_heat
            qg (k) = qg (k) - dq
#endif
            qs (k) = qmin_adj
        endif
        ! if graupel < qmin_adj, borrow from rain
        if (qg (k) < qmin_adj) then
            dq = qmin_adj - qg (k)
#ifdef compensate_q_heat
            tt (i,k) = tt (i,k) + dq * icpk (k) ! heating
            qr (k) = qr (k) - dq
#endif
            qg (k) = qmin_adj
        endif
        
        ! -----------------------------------------------------------------------
        ! liquid phase:
        ! -----------------------------------------------------------------------
        
        ! if rain < qmin_adj, borrow from cloud water
        if (qr (k) < qmin_adj) then
            dq = qmin_adj - qr (k)
#ifdef compensate_q_heat
            ql (k) = ql (k) - dq
#endif
            qr (k) = qmin_adj
        endif
        ! if cloud water < qmin_adj, borrow from water vapor
        if (ql (k) < qmin_adj) then
            dq = qmin_adj - ql (k)
#ifdef compensate_q_heat
            tt (i,k) = tt (i,k) + dq * lcpk (k) ! heating
            qv (k) = qv (k) - dq
#endif
            ql (k) = qmin_adj
        endif
        
    enddo
    
    ! -----------------------------------------------------------------------
    ! fix water vapor; borrow from below
    ! -----------------------------------------------------------------------
    
!    do k = 2, km
    do k = km, 2, -1
        if (qv (k) < qmin_adj) then
            dq = ( qmin_adj - qv (k) ) * dp (i, k) / dp (i, k - 1)
#ifdef compensate_q_heat
!qv has much larger negative values then other species
!better not be compensated from surroundings
!            qv (k - 1) = qv (k - 1) - dq
#endif
            qv (k) = qmin_adj
        endif
    enddo
    
    ! -----------------------------------------------------------------------
    ! bottom layer; borrow from above
    ! -----------------------------------------------------------------------
    
    if (qv (1) < qmin_adj .and. qv (2) > qmin_adj) then
!        dq = min ( - ( qv (1) - qmin_adj ) * dp (i, 1), qv (2) * dp (i, 2))
        dq = min ( qmin_adj - qv (1) , qv (2) * dp (i, 2) / dp (i, 1) )
#ifdef compensate_q_heat
!qv has much larger negative values then other species
!better not be compensated from surroundings
!        qv (2) = qv (2) - dq
#endif
        qv (1) = qmin_adj
!xb141 >>>
    elseif ( qv (1) < qmin_adj .and. qv (2) < qmin_adj ) then
        qv (2) = qmin_adj
        qv (1) = qmin_adj
!xb141 <<<
    endif

    ! -----------------------------------------------------------------------
    ! update six hydrometers
    ! -----------------------------------------------------------------------
    do k = 1, km
        qt (i,    k) = qv (k)
        qt (i,ncw+k) = ql (k)
        qt (i,nrw+k) = qr (k)
        qt (i,niw+k) = qi (k)
        qt (i,nsw+k) = qs (k)
        qt (i,ngl+k) = qg (k)
    enddo

    enddo
    
end subroutine fix_neg_q6
