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
!                        MODULE INTERNALPROFILE
!
!---------------------------------------------------------------------
! Boundary layer general profile
! by Jean-Michel Brankart, January 2020
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! ComputeInternalProfile : compute velocity profile with given parameters
! ----------------------------------------------------------------------
MODULE InternalProfile
  use GeneralProfile
  IMPLICIT NONE
  PRIVATE

  PUBLIC ComputeInternalProfile, ComputeExternalProfile, ComputeLogProfile

  ! Public arrays: external profile (uext), internal profile (uint),
  !                log profile (ulog)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC, SAVE :: uext
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC, SAVE :: uint
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC, SAVE :: ulog

  CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE ComputeInternalProfile(gamma,g0)
!----------------------------------------------------------------------
! ** Purpose : compute velocity profile with given parameters
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( in ) :: gamma,g0

    REAL(KIND=8) :: umax

    ! Check parameters
    IF (accuracy.LE.0.) STOP 'Bad parameter: accuracy'
    IF (search.LE.1) STOP 'Bad parameter: search'

    IF (nprint.GT.0) THEN
      print *, '----------------------------'
      print *, 'Internal profile computation'
      print *, '----------------------------'
      print *, 'Pysical parameters:'
      print *, '  Flow type:',flowtype
      print *, '  Gamma:',gamma
      print *, '  Wall slope:',g0
      print *, 'Numerical parameters:'
      print *, '  Integration type:',integtype
      print *, '  Grid resolution factor:',nres
      print *, '  Accuracy:',accuracy
      print *, '  Search:',search
    ENDIF

    if (allocated(uint)) deallocate(uint)
    allocate(uint(n))

    IF (gamma>1.) THEN
      CALL ComputeU(g0,umax)
    ELSE
      uint=0.
    ENDIF

    END SUBROUTINE ComputeInternalProfile
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE ComputeExternalProfile(gamma)
!----------------------------------------------------------------------
! ** Purpose : compute velocity profile with given parameters
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( in ) :: gamma

    REAL(KIND=8), PARAMETER :: pi=3.1415926535897932384626
    REAL(KIND=8) :: th2, z, arglog
    INTEGER :: k

    IF (nprint.GT.0) THEN
      print *, '----------------------------'
      print *, 'External profile computation'
      print *, '----------------------------'
      print *, 'Pysical parameters:'
      print *, '  Gamma:',gamma
    ENDIF

    if (allocated(uext)) deallocate(uext)
    allocate(uext(n))

    SELECT CASE(flowtype)
    CASE('Poiseuille')
      IF (gamma>0.) THEN
        th2 = gamma * SQRT(Qparam)
        th2 = SQRT(th2)
        DO k=1,n
          z = 0.5_8-x(k)
          arglog = COS(th2*z)
          IF (arglog>0.) THEN
            uext(k) = 1._8 + LOG(arglog) / gamma
          ELSE
            uext(k) = 0.
          ENDIF
        ENDDO
      ELSE
        uext = 1.
      ENDIF
    CASE('Couette')
      IF (gamma>0.) THEN
        th2 = gamma * SQRT(Qparam)
        !th2 = SQRT(th2/2._8)
        th2 = SQRT(th2)/(2._8*sqrt(2._8))
        DO k=1,n
          z = 0.5_8-x(k)
          arglog = TAN( pi/4._8 - th2*z )
          IF (arglog>0.) THEN
            uext(k) = 1._8 + sqrt(2._8) * LOG(arglog) / gamma
          ELSE
            uext(k) = 0.
          ENDIF
        ENDDO
      ELSE
        uext = 1.
      ENDIF
    CASE DEFAULT
      STOP 'Bad flow type'
    END SELECT

    END SUBROUTINE ComputeExternalProfile
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE ComputeLogProfile(gamma,g0)
!----------------------------------------------------------------------
! ** Purpose : compute velocity profile with given parameters
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( in ) :: gamma,g0

    INTEGER :: k

    IF (nprint.GT.0) THEN
      print *, '-------------------------------'
      print *, 'Logarithmic profile computation'
      print *, '-------------------------------'
      print *, 'Pysical parameters:'
      print *, '  Gamma:',gamma
      print *, '  Wall slope:',g0
    ENDIF

    if (allocated(ulog)) deallocate(ulog)
    allocate(ulog(n))

    IF (gamma>1.) THEN
      DO k=1,n
        !ulog(k) = LOG(gamma*g0*x(k)+1.0_8) / gamma
        ulog(k) = LOG(gamma*g0*x(k)) / gamma
      ENDDO
    ELSE
      ulog = 0.
    ENDIF

    END SUBROUTINE ComputeLogProfile
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
    uint(1) = 0.
    DO k=1,n-1
      uint(k+1) = uint(k) + (x(k+1)-x(k)) * (g(k+1)+g(k)) / 2
    ENDDO

    ! Find maximum velocity
    umax = MAXVAL(uint)

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

    REAL(KIND=8) :: Q, gouterr

    Q = 0.
    gouterr = Integrate(g0,Q)

    END SUBROUTINE ComputeG
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE InternalProfile
