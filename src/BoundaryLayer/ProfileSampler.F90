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
!                        MODULE PROFILESAMPLER
!
!---------------------------------------------------------------------
! Boundary layer general profile
! by Jean-Michel Brankart, February 2020
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! SampleProfile : sample velocity profile with given parameters
! ----------------------------------------------------------------------
MODULE ProfileSampler
  use GeneralProfile
  use ensdam_storng
  IMPLICIT NONE
  PRIVATE

  PUBLIC SampleProfile

  ! Public arrays: sample profile (usmp), sample gradient
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC, SAVE :: usmp
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC, SAVE :: gsmp
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC, SAVE :: gtst

  ! Sampling parameters
  LOGICAL, PUBLIC, SAVE :: initmax = .FALSE.
  INTEGER, PUBLIC, SAVE :: maxiter = 1000
  INTEGER, PUBLIC, SAVE :: nharmonics = 50
  REAL(KIND=8), PUBLIC, SAVE :: costfac = 0.25
  REAL(KIND=8), PUBLIC, SAVE :: pertstd = 0.1

  ! Parameter of maximum probability profile previously computed
  REAL(KIND=8), SAVE :: saved_gamma = -huge(saved_gamma)
  REAL(KIND=8), SAVE :: saved_lnu = -huge(saved_lnu)

  ! Profile parameters
  REAL(KIND=8), SAVE :: gam2, lnu, psi, cost0, costmin

  ! Number of evaluation of the cost function
  INTEGER, SAVE :: ncost

  CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE SampleProfile(gamma,l)
!----------------------------------------------------------------------
! ** Purpose : compute velocity profile with given parameters
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( in ) :: gamma,l

    REAL(KIND=8) :: umax

    gam2 = gamma*gamma ; lnu = l ; psi = l


    IF ((gamma.NE.saved_gamma).OR.(lnu.NE.saved_lnu)) THEN
      ! Compute maximum probability profile (if not done)
      CALL ComputeProfile(gamma,lnu)
      saved_gamma = gamma ; saved_lnu = lnu
    ENDIF

    if (allocated(usmp)) deallocate(usmp,gsmp,gtst)
    allocate(usmp(2*n-1),gsmp(2*n-1),gtst(2*n-1))

    ! Compute associated minimum cost function
    gtst(1:n) = g(:)
    gtst(2*n-1:n:-1) = g(:)
    costmin = CostFunction()

    print *, 'costmin',costmin

    ! Compute S0 parameter
    cost0 = costfac * costmin
    
    ! Check parameters
    IF (accuracy.LE.0.) STOP 'Bad parameter: accuracy'

    IF (nprint.GT.0) THEN
      print *, '....Sampling new profile'
    ENDIF

    ! Sample gradient
    CALL SampleG()
    ! Compute velocity from gradient
    CALL ComputeU(umax)

    ! Rescale to normalize profilr
    !gsmp(:) = gsmp(:) / umax
    !usmp(:) = usmp(:) / umax

    END SUBROUTINE SampleProfile
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE ComputeU(umax)
!----------------------------------------------------------------------
! ** Purpose : compute velocity profile with given wall slope
!----------------------------------------------------------------------
    REAL(KIND=8), INTENT( out ) :: umax

    INTEGER :: k
    REAL(KIND=8) :: uext

    IF (nprint.GT.1) print *, '......Wall slope:',gsmp(1)

    SELECT CASE(flowtype)
    CASE('Poiseuille')
      uext = 0.
    CASE('Couette')
      uext = 2._8
    CASE DEFAULT
      STOP 'Bad flow type'
    END SELECT

    ! Integrate to compute u from g (from x=0 boundary)
    usmp(1) = 0.
    DO k=1,n-1
      usmp(k+1) = usmp(k) + (x(k+1)-x(k)) * (gsmp(k+1)+gsmp(k)) / 2
    ENDDO

    ! Integrate to compute u from g (from x=1 boundary)
    usmp(2*n-1) = uext
    DO k=1,n-1
      usmp(2*n-k-1) = usmp(2*n-k) - (x(k+1)-x(k)) * (gsmp(2*n-k-1)+gsmp(2*n-k)) / 2
    ENDDO

    ! Find maximum velocity
    umax = MAXVAL(usmp)

    IF (nprint.GT.1) print *, '......Maximum velocity',umax

    END SUBROUTINE ComputeU
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE SampleG()
!----------------------------------------------------------------------
! ** Purpose : sample gradient with Metropolis/Hastings algorithù
!----------------------------------------------------------------------
    REAL(KIND=8) :: cost, costprev, prob, uran, redcost, alpha, beta
    INTEGER :: iter, naccepted, k
    LOGICAL :: accepted

    ! Initialize gradient
    SELECT CASE(flowtype)
    CASE('Poiseuille')

      IF (initmax) THEN
        ! initialize with maximuma probability solution
        gsmp(1:n) = g(:)
        gsmp(2*n-1:n:-1) = - gsmp(1:n)
      ELSE
        ! initialize with constant curvature so that g(1/2)=0
        gsmp(1:n) = g(1) * ( 1.0_8 - 2.0_8 * x(:) )
        gsmp(2*n-1:n:-1) = - gsmp(1:n)
      ENDIF

    CASE('Couette')

      IF (initmax) THEN
        ! initialize with maximuma probability solution
        gsmp(1:n) = g(:)
        gsmp(2*n-1:n:-1) = gsmp(1:n)
      ELSE
        ! initialize with quadratic gradient
        beta = ( 6._8 - g(1) ) * 0.5_8
        alpha = 4._8 * ( 2._8 - beta )
        gsmp(1:n) = 3._8 * alpha * ( x(:) - 0.5_8 ) **2 + beta
        gsmp(2*n-1:n:-1) = gsmp(1:n)
      ENDIF

    CASE DEFAULT
      STOP 'Bad flow type'
    END SELECT


    ! Initialize cost function
    gtst(:) = gsmp(:)
    costprev = CostFunction()

    print *, 'initial cost',costprev
    print *, 'reduced initial cost',(costprev-costmin)/cost0

    ! Metropolis/Hastings iteration
    naccepted = 0
    DO iter=1,maxiter
      ! Sample perturbation from proposal density
      CALL SamplePerturbation()
      ! Evaluate cost function
      cost = CostFunction()
      ! Compute acceptance probability
      prob = EXP( - (cost - costprev) / cost0 )
      ! Decide whether to accept new draw
      CALL kiss_uniform(uran)
      accepted = uran.LT.prob
      ! update gradient if accepted
      IF (accepted) THEN
        gsmp(:) = gtst(:)
        costprev = cost
        naccepted = naccepted + 1
        ! Monitor iteration
        IF (MOD(naccepted,10)==0) THEN
          redcost = (cost-costmin)/cost0
          !print *, '..........iterations:',iter
          !print *, '..........accepted:',naccepted
          !print *, '..........cost:',cost
          !print *, '..........reduced cost:',redcost
          !IF (redcost<1.0) EXIT
        ENDIF
      ENDIF
    ENDDO

    print *, '..........accepted:',naccepted
    print *, '..........cost:',costprev
    print *, '..........reduced cost:',redcost

    END SUBROUTINE SampleG
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    FUNCTION CostFunction()
!----------------------------------------------------------------------
! ** Purpose : Compute cost function
!----------------------------------------------------------------------
    REAL(KIND=8) :: CostFunction

    INTEGER :: k, klin
    REAL(KIND=8) :: S1, S2, gint, dgdx

    ! Select lnu option
    IF (psi.gt.0.5) THEN
      IF (gam2==0.) THEN
        lnu = 0.
      ELSE
        lnu = psi / (gtst(1)*sqrt(gam2))
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

    ! Evaluate viscous part of the cost function (S1)
    S1 = 0.
    DO k=1,n-1
      dgdx = ( gtst(k+1) - gtst(k) ) / ( x(k+1) - x(k) )
      S1 = S1 + ( x(k+1) - x(k) ) * dgdx * dgdx
      dgdx = ( gtst(2*n-k-1) - gtst(2*n-k) ) / ( x(k+1) - x(k) )
      S1 = S1 + ( x(k+1) - x(k) ) * dgdx * dgdx
    ENDDO

    ! Evaluate nonlinear part of the cost function (S2)
    S2 = 0.
    DO k=klin,n-1
      gint = 0.5_8 * ( gtst(k+1) + gtst(k) )
      S2 = S2 + ( x(k+1) - x(k) ) * gint ** 4
      gint = 0.5_8 * ( gtst(2*n-k-1) + gtst(2*n-k) )
      S2 = S2 + ( x(k+1) - x(k) ) * gint ** 4
    ENDDO

    !print *, 'S1S2',S1,S2
    ! Compute the sum of the two parts
    CostFunction = S1 + gam2 * S2
    !CostFunction = gam2 * S2
    !CostFunction = S1

    ncost = ncost +1

    END FUNCTION CostFunction
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE SamplePerturbation()
!----------------------------------------------------------------------
! ** Purpose : sample perturbation from proposal density
!----------------------------------------------------------------------
    REAL(KIND=8), PARAMETER :: pi=3.1415926535897932384626

    REAL(KIND=8) :: kappa, gran, spct, norm, int1, int2, facp, facm, gg
    INTEGER :: i, k, kx

    ! Initialize perturbation to zero
    gtst(:) = 0.

    ! Accumulate Gaussian perturbations along harmonic functions
    ! (with zero condition at the boundaries)
    norm = 0.
    DO i=1,nharmonics
      kappa = REAL(i,8) * pi
      spct = 1./(kappa*kappa)
      norm = norm + spct * spct
      CALL kiss_gaussian(gran)
      gran = gran * spct
      DO k=1,n
        gtst(k) = gtst(k) + gran * SIN(kappa*x(k))
      ENDDO
      DO k = n+1,2*n-1
        kx = 2*n-k
        gtst(k) = gtst(k) + gran * SIN(kappa*(1.-x(kx)))
      ENDDO
    ENDDO

    ! Rescale to the rquired standard deviation
    gtst(:) = gtst(:) * pertstd / SQRT(norm)

    !print *, 'perturb',minval(gtst),maxval(gtst)

    ! Impose the constraint of a zero integral perturbation
    int1 = 0. ; int2 = 0.
    DO k=1,n-1
      gg = gtst(k+1)+gtst(k)
      IF (gg>0.) THEN
        int1 = int1 + (x(k+1)-x(k)) * gg / 2
      ELSE
        int2 = int2 + (x(k+1)-x(k)) * gg / 2
      ENDIF
    ENDDO
    DO k=1,n-1
      gg = gtst(2*n-k-1)+gtst(2*n-k)
      IF (gg>0.) THEN
        int1 = int1 + (x(k+1)-x(k)) * gg / 2
      ELSE
        int2 = int2 + (x(k+1)-x(k)) * gg / 2
      ENDIF
    ENDDO

    facp = - int2 / SQRT(int1*int1+int2*int2)
    facm = int1 / SQRT(int1*int1+int2*int2)

    WHERE(gtst>0.) gtst = gtst * facp
    WHERE(gtst<0.) gtst = gtst * facm

    ! Add expected value of proposal density
    gtst(:) = gsmp(:) + gtst(:)

    END SUBROUTINE SamplePerturbation
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE ProfileSampler
