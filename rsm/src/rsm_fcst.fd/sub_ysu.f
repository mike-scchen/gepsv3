!-------------------------------------------------------------------------------
      subroutine tridi1n_ysu(cl,cm,cu,r1,r2,au,f1,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------------------
!
      integer, intent(in )      ::     its,ite, kts,kte, nt
!
      real, dimension( its:ite, kts+1:kte+1 ),intent(in   )  ::       cl
!
      real, dimension( its:ite, kts:kte ),intent(in   )  ::       cm, r1
      real, dimension( its:ite, kts:kte,nt ),intent(in   )  ::        r2
!
      real, dimension( its:ite, kts:kte ),intent(inout)  ::   au, cu, f1
      real, dimension( its:ite, kts:kte,nt ),intent(inout)  ::        f2
!
      real     fk
      integer  i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
      l = ite
      n = kte
!
      do i = its,l
        fk = 1./cm(i,1)
        au(i,1) = fk*cu(i,1)
        f1(i,1) = fk*r1(i,1)
      enddo
!
      do it = 1,nt
        do i = its,l
          fk = 1./cm(i,1)
          f2(i,1,it) = fk*r2(i,1,it)
        enddo
      enddo
!
      do k = kts+1,n-1
        do i = its,l
          fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          f1(i,k) = fk*(r1(i,k)-cl(i,k)*f1(i,k-1))
        enddo
      enddo
!
      do it = 1,nt
        do k = kts+1,n-1
          do i = its,l
            fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
            f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
          enddo
        enddo
      enddo
!
      do i = its,l
        fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        f1(i,n) = fk*(r1(i,n)-cl(i,n)*f1(i,n-1))
      enddo
!
      do it = 1,nt
        do i = its,l
          fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
          f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
        enddo
      enddo
!
      do k = n-1,kts,-1
        do i = its,l
          f1(i,k) = f1(i,k)-au(i,k)*f1(i,k+1)
        enddo
      enddo
!
      do it = 1,nt
        do k = n-1,kts,-1
          do i = its,l
            f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
          enddo
        enddo
      enddo
!
      end subroutine tridi1n_ysu
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
      subroutine rtridin_ysu(cl,cm,cu,r2,au,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------------------
!
        integer, intent(in )      ::     its,ite, kts,kte, nt
!
      real, dimension( its:ite, kts+1:kte+1 ),intent(in   )  ::       cl
!
      real, dimension( its:ite, kts:kte ),intent(in   )  ::           cm
      real, dimension( its:ite, kts:kte,nt ),intent(in   )  ::        r2
!
      real, dimension( its:ite, kts:kte ),intent(inout)  ::       au, cu
      real, dimension( its:ite, kts:kte,nt ),intent(inout)  ::        f2
!
      real     fk
      integer  i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
      l = ite
      n = kte
!
      do it = 1,nt
        do i = its,l
          fk = 1./cm(i,1)
          au(i,1) = fk*cu(i,1)
          f2(i,1,it) = fk*r2(i,1,it)
        enddo
      enddo
!
      do it = 1,nt
        do k = kts+1,n-1
          do i = its,l
            fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
            au(i,k) = fk*cu(i,k)
            f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
          enddo
        enddo
      enddo
!
      do it = 1,nt
        do i = its,l
          fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
          f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
        enddo
      enddo
!
      do it = 1,nt
        do k = n-1,kts,-1
          do i = its,l
            f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
          enddo
        enddo
      enddo
!
        end subroutine rtridin_ysu
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
      subroutine rysuinit(rublten,rvblten,rthblten,rqvblten,
     &                   rqcblten,rqiblten,p_qi,p_first_scalar,
     &                   restart, allowed_to_read,
     &                   ids, ide, jds, jde, kds, kde,
     &                   ims, ime, jms, jme, ktes, ktee,
     &                   its, ite, jts, jte, kts, kte                 )
!-------------------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------------------
!
      logical , intent(in)          :: restart, allowed_to_read
      integer , intent(in)          ::  ids, ide, jds, jde, kds, kde,
     &                                  ims, ime, jms, jme, ktes, ktee,
     &                                  its, ite, jts, jte, kts, kte
      integer , intent(in)          ::  p_qi,p_first_scalar
      real,dimension(ims:ime,ktes:ktee,jms:jme),intent(out) ::
     &              rublten,rvblten,rthblten,rqvblten,rqcblten,rqiblten
      integer  i, j, k, itf, jtf, ktf
!
        jtf = min0(jte,jde-1)
        ktf = min0(kte,kde-1)
        itf = min0(ite,ide-1)
!
        if(.not.restart)then
          do j = jts,jtf
            do k = kts,ktf
              do i = its,itf
                 rublten(i,k,j) = 0.
                 rvblten(i,k,j) = 0.
                 rthblten(i,k,j) = 0.
                 rqvblten(i,k,j) = 0.
                 rqcblten(i,k,j) = 0.
              enddo
            enddo
          enddo
        endif
!
        if (p_qi .ge. p_first_scalar .and. .not.restart) then
          do j = jts,jtf
            do k = kts,ktf
              do i = its,itf
                rqiblten(i,k,j) = 0.
              enddo
            enddo
          enddo
        endif
!
        end subroutine rysuinit
!-------------------------------------------------------------------------------
! ==================================================================
      SUBROUTINE RGET_PBLH(KTS,KTE,zi,thetav1D,qke1D,zw1D,dz1D,landsea)
! Copied from MYNN PBL

      !---------------------------------------------------------------
      !             NOTES ON THE PBLH FORMULATION
      !
      !The 1.5-theta-increase method defines PBL heights as the level at
      !which the potential temperature first exceeds the minimum potential
      !temperature within the boundary layer by 1.5 K. When applied to
      !observed temperatures, this method has been shown to produce PBL-
      !height estimates that are unbiased relative to profiler-based
      !estimates (Nielsen-Gammon et al. 2008). However, their study did not
      !include LLJs. Banta and Pichugina (2008) show that a TKE-based
      !threshold is a good estimate of the PBL height in LLJs. Therefore,
      !a hybrid definition is implemented that uses both methods, weighting
      !the TKE-method more during stable conditions (PBLH < 400 m).
      !A variable tke threshold (TKEeps) is used since no hard-wired
      !value could be found to work best in all conditions.
      !---------------------------------------------------------------

      INTEGER,INTENT(IN) :: KTS,KTE
      REAL, INTENT(OUT) :: zi
      REAL, INTENT(IN) :: landsea
      REAL, DIMENSION(KTS:KTE), INTENT(IN) :: thetav1D, qke1D, dz1D
      REAL, DIMENSION(KTS:KTE+1), INTENT(IN) :: zw1D
      !LOCAL VARS
      REAL ::  PBLH_TKE,qtke,qtkem1,wt,maxqke,TKEeps,minthv
      REAL :: delt_thv   !delta theta-v; dependent on land/sea point
      REAL, PARAMETER :: sbl_lim  = 200. !Theta-v PBL lower limit of trust (m).
      REAL, PARAMETER :: sbl_damp = 400. !Damping range for averaging with TKE-based PBLH (m).
      INTEGER :: I,J,K,kthv,ktke

      !FIND MAX TKE AND MIN THETAV IN THE LOWEST 500 M
      k = kts+1
      kthv = 1
      ktke = 1
      maxqke = 0.
      minthv = 9.E9

      DO WHILE (zw1D(k) .LE. 500.)
        qtke  =MAX(Qke1D(k),0.)   ! maximum QKE
         IF (maxqke < qtke) then
            maxqke = qtke
            ktke = k
         ENDIF
         IF (minthv > thetav1D(k)) then
             minthv = thetav1D(k)
             kthv = k
         ENDIF
         k = k+1
      ENDDO
      !TKEeps = maxtke/20. = maxqke/40.
      TKEeps = maxqke/40.
      TKEeps = MAX(TKEeps,0.025)
      TKEeps = MIN(TKEeps,0.25)

      !FIND THETAV-BASED PBLH (BEST FOR DAYTIME).
      zi=0.
      k = kthv+1
      IF((landsea-1.5).GE.0)THEN
      ! WATER
          delt_thv = 0.75
      ELSE
      ! LAND
          delt_thv = 1.5
      ENDIF

      zi=0.
      k = kthv+1
      DO WHILE (zi .EQ. 0.)
         IF (thetav1D(k) .GE. (minthv + delt_thv))THEN
            zi = zw1D(k) - dz1D(k-1)*
     &          MIN((thetav1D(k)-(minthv + delt_thv))/
     &          MAX(thetav1D(k)-thetav1D(k-1),1E-6),1.0)
         ENDIF
        k = k+1
         IF (k .EQ. kte-1) zi = zw1D(kts+1) !EXIT SAFEGUARD
      ENDDO

      !print*,"IN GET_PBLH:",thsfc,zi
      !FOR STABLE BOUNDARY LAYERS, USE TKE METHOD TO COMPLEMENT THE
      !THETAV-BASED DEFINITION (WHEN THE THETA-V BASED PBLH IS BELOW ~0.5 KM).
      !THE TANH WEIGHTING FUNCTION WILL MAKE THE TKE-BASED DEFINITION NEGLIGIBLE
      !WHEN THE THETA-V-BASED DEFINITION IS ABOVE ~1 KM.
      !FIND TKE-BASED PBLH (BEST FOR NOCTURNAL/STABLE CONDITIONS).

      PBLH_TKE=0.
      k = ktke+1
      DO WHILE (PBLH_TKE .EQ. 0.)
        !QKE CAN BE NEGATIVE (IF CKmod == 0)... MAKE TKE NON-NEGATIVE.
         qtke  =MAX(Qke1D(k)/2.,0.)      ! maximum TKE
         qtkem1=MAX(Qke1D(k-1)/2.,0.)
         IF (qtke .LE. TKEeps) THEN
               PBLH_TKE = zw1D(k) - dz1D(k-1)*
     &           MIN((TKEeps-qtke)/MAX(qtkem1-qtke, 1E-6), 1.0)
             !IN CASE OF NEAR ZERO TKE, SET PBLH = LOWEST LEVEL.
             PBLH_TKE = MAX(PBLH_TKE,zw1D(kts+1))
             !print *,"PBLH_TKE:",i,j,PBLH_TKE, Qke1D(k)/2., zw1D(kts+1)
         ENDIF
         k = k+1
         IF (k .EQ. kte-1) PBLH_TKE = zw1D(kts+1) !EXIT SAFEGUARD
      ENDDO

    !BLEND THE TWO PBLH TYPES HERE:

      wt=.5*TANH((zi - sbl_lim)/sbl_damp) + .5
      zi=PBLH_TKE*(1.-wt) + zi*wt

      END SUBROUTINE RGET_PBLH
! ==================================================================

!-------------------------------------------------------------------------------
