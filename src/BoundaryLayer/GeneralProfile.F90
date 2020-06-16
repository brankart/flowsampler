!---------------------------------------------------------------------
! Copyright: CNRS - Université de Grenoble Alpes
!
! Contributors : Jean-Michel Brankart
!
! Jean-Michel.Brankart@univ-grenoble-alpes.fr
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
!---------------------------------------------------------------------
!
!                        MODULE GENERALPROFILE
!
!---------------------------------------------------------------------
! Boundary layer general profile (Poiseuille and Couette flows)
! by Jean-Michel Brankart, January 2020
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ComputeProfile : compute velocity profile with given parameters
! ----------------------------------------------------------------------
MODULE GeneralProfile
  IMPLICIT NONE
  PRIVATE

  PUBLIC ComputeProfile, Integrate

  ! Public parameters 
  INTEGER, PUBLIC, SAVE :: nprint = 4   ! amount of screen notifications
  INTEGER, PUBLIC, SAVE :: nres = 1000  ! grid resolution factor
  REAL(KIND=8), PUBLIC, SAVE :: accuracy = 1e-6  ! required accuracy
  INTEGER, PUBLIC, SAVE :: search = 2  ! search factor to bracket roots
  CHARACTER(len=16), PUBLIC, SAVE :: flowtype = 'Poiseuille' ! flow type
  CHARACTER(len=16), PUBLIC, SAVE :: integtype = 'Explicit'  ! integration type
  LOGICAL, PUBLIC, SAVE :: usebeta = .FALSE. ! use a beta term in the equations

  ! Public arrays: size (n), coordinate (x), gradient (g) and velocity (u)
  INTEGER, PUBLIC, SAVE :: n
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC, SAVE :: x, g, u
  REAL(KIND=8), PUBLIC, SAVE :: Qparam ! Q profile parameter
  REAL(KIND=8), PUBLIC, SAVE :: g0save ! Wall velocity gradient

  ! Profile parameters
  REAL(KIND=8), SAVE :: gam2, lnu, psi, bet2

  ! Other parameters
  INTEGER, SAVE :: ninteg=0 ! total number of integrations

  CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE ComputeProfile(gamma,l)
!----------------------------------------------------------------------
! ** Purpose : compute velocity profile with given parameters
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( in ) :: gamma, l

    REAL(KIND=8) :: g0low, g0up, g0, umaxlow, umaxup, umax, gamcrit

    gam2 = gamma*gamma ; lnu = l ; psi = l ; bet2 = 0.

    ! Check parameters
    IF (accuracy.LE.0.) STOP 'Bad parameter: accuracy'
    IF (search.LE.1) STOP 'Bad parameter: search'

    IF (nprint.GT.0) THEN
      print *, '---------------------------'
      print *, 'General profile computation'
      print *, '---------------------------'
      print *, 'Pysical parameters:'
      print *, '  Flow type:',flowtype
      print *, '  Gamma:',sqrt(gam2)
      print *, '  Lnu:',lnu
      print *, 'Numerical parameters:'
      print *, '  Integration type:',integtype
      print *, '  Grid resolution factor:',nres
      print *, '  Accuracy:',accuracy
      print *, '  Search:',search
    ENDIF

    ! Define grid coordinates (x)
    CALL DefineGrid()

    IF (nprint.GT.0) print *, '  Grid size:',n

    ! Starting lower bound for wall slope g0
    SELECT CASE(flowtype)
    CASE('Poiseuille')
      IF (lnu<3.) THEN
        gamcrit = 1.
      ELSEIF (lnu<6.) THEN
        gamcrit = 4.
      ELSEIF (lnu<9.) THEN
        gamcrit = 5.
      ELSE
        gamcrit = 7.
      ENDIF
      IF (gam2<gamcrit*gamcrit) THEN
        g0low = 4._8  ;  g0up = g0low
      ELSE
        g0low = exp(sqrt(gam2))/sqrt(gam2) ;  g0up = g0low
      ENDIF
    CASE('Couette')
      IF (gam2<1.) THEN
        g0low = 2._8  ;  g0up = g0low
      ELSE
        g0low = exp(sqrt(gam2))/sqrt(gam2) ;  g0up = g0low
      ENDIF
    CASE DEFAULT
      STOP 'Bad flow type'
    END SELECT

    ! Initialize number of integrations
    ninteg = 0

    ! First profile computation with starting wall slope
    CALL ComputeU(g0low,umaxlow)
    umaxup = umaxlow
    g0 = g0low ; umax = umaxlow

    ! Iterate to find upper bound to g0
    DO WHILE(umaxup<1._8)
      g0low = g0up ; umaxlow = umaxup
      g0up = g0up*search
      CALL ComputeU(g0up,umaxup)
      g0 = g0up ; umax = umaxup
    ENDDO
    
    ! Iterate to find lower bound to g0
    DO WHILE(umaxlow>1._8)
      g0up = g0low ; umaxup = umaxlow
      g0low = g0low/search
      CALL ComputeU(g0low,umaxlow)
      g0 = g0low ; umax = umaxlow
    ENDDO

    ! Iterate to find g0 with the required accuracy
    DO WHILE(abs(umax-1._8)>accuracy)
      g0 = g0low + (g0up-g0low) * (1._8-umaxlow) / (umaxup-umaxlow)
      CALL ComputeU(g0,umax)
      IF (umax>1._8) THEN
        g0up = g0 ; umaxup = umax
      ELSE
        g0low = g0 ; umaxlow = umax
      ENDIF
    ENDDO

    g0save = g0

    IF (nprint.GT.0) THEN
      print *, 'Final diagnostics:'
      print *, '  Wall slope:',g0save
      print *, '  Q parameter:',Qparam
      print *, '  Number of iterations:',ninteg
    ENDIF

    END SUBROUTINE ComputeProfile
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE ComputeU(g0,umax)
!----------------------------------------------------------------------
! ** Purpose : compute velocity profile with given wall slope
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( in ) :: g0
    REAL(KIND=8), INTENT( out ) :: umax

    INTEGER :: k

    IF (nprint.GT.1) print *, 'Wall slope:',g0

    ! Compute velocity gradient
    CALL ComputeG(g0)

    ! Integrate to compute u from g
    u(1) = 0.
    DO k=1,n-1
      u(k+1) = u(k) + (x(k+1)-x(k)) * (g(k+1)+g(k)) / 2
    ENDDO

    ! Find maximum velocity
    umax = MAXVAL(u)

    IF (nprint.GT.1) print *, 'Maximum velocity',umax

    END SUBROUTINE ComputeU
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE ComputeG(g0)
!----------------------------------------------------------------------
! ** Purpose : compute velocity gradient with given wall slope
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( in ) :: g0

    REAL(KIND=8) :: Qlow, Qup, Q, gouterrlow, gouterrup, gouterr, gamcrit

    ! Starting lower and upper bound for Q
    SELECT CASE(flowtype)
    CASE('Poiseuille')
      IF (lnu<3.) THEN
        gamcrit = 1.
      ELSEIF (lnu<6.) THEN
        gamcrit = 4.
      ELSEIF (lnu<9.) THEN
        gamcrit = 5.
      ELSE
        gamcrit = 7.
      ENDIF
      IF (gam2<gamcrit*gamcrit) THEN
        Qup = 64._8 ; Qlow = Qup
      ELSE
        Qup = 1. ; Qlow = Qup
      ENDIF
    CASE('Couette')
      IF (gam2==0.) THEN
        Qup = 0. ; Qlow = Qup
      ELSE
        IF (gam2<1.) THEN
          Qup = 64._8 ; Qlow = Qup
        ELSE
          Qup = 1. ; Qlow = Qup
        ENDIF
      ENDIF
    CASE DEFAULT
      STOP 'Bad flow type'
    END SELECT
    IF (nprint.GT.2) print *, '  Starting value of Q:',Qlow

    gouterrlow = Integrate(g0,Qlow)
    gouterrup = gouterrlow
    Q = Qlow ; gouterr = gouterrlow

    ! Iterate to find lower bound to Q
    DO WHILE (gouterrlow<0.)
      Qup = Qlow ; gouterrup = gouterrlow
      Qlow = Qlow/search
      gouterrlow = Integrate(g0,Qlow)
      Q = Qlow ; gouterr = gouterrlow
    ENDDO

    ! Iterate to find upper bound to Q
    DO WHILE (gouterrup>0.)
      Qlow = Qup ; gouterrlow = gouterrup
      Qup = Qup*search
      gouterrup = Integrate(g0,Qup)
      Q = Qup ; gouterr = gouterrup
    ENDDO

    IF (nprint.GT.2) print *, '  Q range:',Qlow,Qup

    ! Iterate to find Q with the required accuracy
    DO WHILE (abs(gouterr)>accuracy)
      Q = Qlow - gouterrlow * (Qup-Qlow) / (gouterrup-gouterrlow)
      gouterr = Integrate(g0,Q)
      IF (gouterr>0.) THEN
        Qlow = Q ; gouterrlow = gouterr
      ELSE
        Qup = Q ; gouterrup = gouterr
      ENDIF
    ENDDO

    IF (nprint.GT.2) print *, '  Q final:',Q

    Qparam = Q

    END SUBROUTINE ComputeG
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    FUNCTION Integrate(g0,Q)
!----------------------------------------------------------------------
! ** Purpose : Integrate differential equation with given parameters
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( in ) :: g0,Q
    REAL(KIND=8) :: Integrate

    INTEGER :: k, klin
    REAL(KIND=8) :: QQ, err, Clinref, Clin
    LOGICAL :: itersense, validp, validc

    IF (nprint.GT.3) print *, '    Q constant:',Q

    SELECT CASE(flowtype)
    CASE('Poiseuille')
      QQ = Q
    CASE('Couette')
      QQ = -Q
    CASE DEFAULT
      STOP 'Bad flow type'
    END SELECT

    ! Select lnu option
    IF (psi.gt.0.5) THEN
      IF (usebeta) THEN
        lnu = 0.
        bet2 = psi * sqrt(gam2) / g0
        bet2 = bet2 * bet2
      ELSE
        IF (gam2==0.) THEN
          lnu = 0.
        ELSE
          lnu = psi / (g0*sqrt(gam2))
        ENDIF
      ENDIF
    ENDIF

    ! Find thickness of viscous sublayer
    IF (lnu>=x(n)) THEN
      klin=n
    ELSE
      klin=1
      DO WHILE (x(klin)<lnu)
        klin=klin+1
      ENDDO
    ENDIF

    IF (klin>1) THEN

      ! Find slope of viscous sublayer
      itersense = .TRUE.
      IF (QQ>0.) THEN
        Clinref = 0.5 * (SQRT(QQ) + g0 / x(klin))
      ELSE
        Clinref = 0.5 * g0 / x(klin)
      ENDIF

      validc = .TRUE.
      validp = Clinref*Clinref>=QQ
      IF (validp) validc = g0>SQRT(SQRT((Clinref*Clinref-QQ)/gam2))

      IF (validp.AND.validc) THEN
        Clin = ( g0 - SQRT(SQRT((Clinref*Clinref-QQ)/gam2)) ) / x(klin)
      ELSE
        itersense = .FALSE.
        IF (QQ+gam2*(g0-Clinref*x(klin)**4)<0.) STOP 'Impossible case'
        Clin = SQRT(QQ+gam2*(g0-Clinref*x(klin)**4))
      ENDIF
      err = ABS((Clin-Clinref)/Clin)

      DO WHILE (err>accuracy)
        Clinref = Clin
        IF (itersense) THEN
          validc = .TRUE.
          validp = Clinref*Clinref>=QQ
          IF (validp) validc = g0>SQRT(SQRT((Clinref*Clinref-QQ)/gam2))

          IF (validp.AND.validc) THEN
            Clin = ( g0 - SQRT(SQRT((Clinref*Clinref-QQ)/gam2)) ) / x(klin)
          ELSE
            itersense = .FALSE.
            Clin = SQRT(QQ+gam2*(g0-Clinref*x(klin)**4))
          ENDIF
        ELSE
          !IF (Clinref>0.) THEN
            Clin = SQRT(QQ+gam2*(g0-Clinref*x(klin)**4))
          !ELSE
          !  itersense = .TRUE.
          !  Clin = ( g0 - SQRT(SQRT((Clinref*Clinref-QQ)/gam2)) ) / x(klin)
          !ENDIF
        ENDIF
        err = ABS((Clin-Clinref)/Clin)
        !print *, 'Clin',Clinref,Clin
      ENDDO

      ! Integrate in viscous sublayer
      DO k=1,klin
        g(k) = g0 - Clin * x(k)
      ENDDO

    ELSE

      g(1) = g0

    ENDIF

    ! Integrate nonlinear profile to compute g
    DO k=klin,n-1
      g(k+1) = g(k) + Dg( g(k) , x(k+1)-x(k) , QQ )
    ENDDO

    ! Diagnose discrepancy to target final slope
    SELECT CASE(flowtype)
    CASE('Poiseuille')
      Integrate = g(n)
    CASE('Couette')
      IF (gam2==0.) THEN
        Integrate = g(n) - g0
      ELSE
        Integrate = g(n) - ( Q / gam2 ) ** 0.25_8
      ENDIF
    CASE DEFAULT
      STOP 'Bad flow type'
    END SELECT

    IF (nprint.GT.3) print *, '    Final slope error:',Integrate

    ninteg = ninteg +1

    END FUNCTION Integrate
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    FUNCTION Dg(g,dx,Q)
!----------------------------------------------------------------------
! ** Purpose : Integrate differential equation with given parameters
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( in ) :: g,dx,Q
    REAL(KIND=8) :: Dg

    REAL(KIND=8) :: Dgbefore

    SELECT CASE(integtype)
    CASE('Explicit')
      Dg = - dx * f(Q,g)
    CASE('Implicit')
      Dg = 0. ; Dgbefore = 1.
      DO WHILE (ABS((Dg-Dgbefore)/Dgbefore)>accuracy)
        Dgbefore = Dg
        Dg = - dx * f(Q,g+Dgbefore/2)
      ENDDO
    CASE DEFAULT
      STOP 'Bad integration type'
    END SELECT

    CONTAINS
      FUNCTION f(Q,g)
        REAL(KIND=8), INTENT( in ) :: Q,g
        REAL(KIND=8) :: f

        REAL(KIND=8) :: gsqr

        gsqr = g * g

        !f = Q + gam2 * gsqr * gsqr - bet2 * gsqr * gsqr * gsqr
        f = Q + gam2 * gsqr * gsqr

        IF (f>=0.) f = SQRT(f)
        IF (f<0.) f = SQRT(-f)

      END FUNCTION f
    END FUNCTION Dg
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE DefineGrid()

    REAL(KIND=8) :: l, xx, dx
    INTEGER :: k

    ! Compute minimal length scale
    IF ((lnu.gt.0.).AND.(lnu.lt.0.5)) THEN
        l = sqrt(gam2)*exp(-sqrt(gam2))/nres
    ELSE
      IF (gam2.gt.2.) THEN
        l = sqrt(gam2)*exp(-sqrt(gam2))/nres
      ELSE
        l = 0.5_8 /nres
      ENDIF
    ENDIF

    ! Compute size of the grid
    xx = 0.5_8 ; n = 1
    DO WHILE (xx>l)
      dx = xx/nres ; xx = xx-dx ; n = n+1
    ENDDO
    DO WHILE (xx>dx)
      xx = xx-dx ; n = n+1
    ENDDO
    xx = 0. ; n = n+1

    ! Allocate arrays
    if (allocated(x)) deallocate(x)
    if (allocated(g)) deallocate(g)
    if (allocated(u)) deallocate(u)
    allocate(x(n),g(n),u(n))

    ! Fill grid array x
    k = n ; x(n) = 0.5_8 ; xx = 0.5_8
    DO WHILE (xx>l)
      dx = xx/nres ; xx = xx-dx ; k = k-1
      x(k) = xx
    ENDDO
    DO WHILE (xx>dx)
      xx = xx-dx ; k = k-1
      x(k) = xx
    ENDDO
    k = k-1 ; x(k) = 0.

    END SUBROUTINE DefineGrid
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE GeneralProfile
