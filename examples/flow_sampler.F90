program flow_sampler
  use flowsampler_grid
  use flowsampler_kinematics
  use flowsampler_dynamics
  use ensdam_mcmc_update
  use ensdam_storng
  use ensdam_storfg
  use ensdam_stoutil
  use ensdam_sphylm
  use ensdam_ensaugm
  implicit none

  ! Include program parameters
  include 'flow_sampler_parameters.h90'

  ! Size of vectors
  integer :: n ! size of state vector (for local subdomain if parallel computing)
  integer :: smax ! maximum degree of the spherical harmonics

  ! Main program arrays: grid, ensemble
  real(kind=8), dimension(:), allocatable :: x, y ! coordinates of state variable
  real(kind=8), dimension(:), allocatable :: a ! area associated to each grid point
  real(kind=8), dimension(:,:,:), allocatable :: proposal ! proposal distribution
  real(kind=8), dimension(:,:), allocatable :: sample ! output sample
  integer, dimension(:,:), allocatable :: indx ! index of grid location
  real(kind=8) :: total_area
  integer, dimension(:), allocatable :: tmp_sample ! temporary storage

  ! Full grid arrays
  real(kind=8), dimension(:,:), allocatable :: psi, u, v, unorm, zeta
  real(kind=8), dimension(:,:), allocatable :: u0, v0, zeta0, zetaforc
  real(kind=8), dimension(:,:), allocatable :: advresidual, roughness, slope_penalty
  real(kind=8) :: uvrms, zetarms

  ! Arrays for scale separation
  real(kind=8), dimension(:,:), allocatable :: spct ! spectrum in the basis of the spherical harmonics

  ! Temporary storage
  real(kind=8) :: uran, amax, tmp ! temporary storage
  real(kind=8) :: mean, std, misfit, misfit2 ! temporary storage
  real(kind=8), dimension(:), allocatable :: tmp_state ! temporary storage
  real(kind=8), dimension(:), allocatable :: tmp_state2 ! temporary storage
  integer(KIND=8) :: zseed1, zseed2, zseed3, zseed4 ! seeds for random number generator
  character(len=256) :: outfile, infile
  integer :: namunit, ios

  ! Vector indices
  integer :: jm ! index of ensemble member
  integer :: jn ! index of state variable
  integer :: js ! index of scale in multiple scale ensemble
  integer :: kt ! time step
  character(len=4) :: tagkt

  ! Parallel computation variables
  integer :: nproc=1  ! Number of processors
  integer :: iproc=0  ! Current processor index
  integer :: jiproc   ! Processor index

  ! Scores
  real(kind=8) :: score, score_prev
  integer :: iteration_index
  integer :: ktsave = -1
  logical :: conv_test

  ! Diagnostic
  integer :: njo=0  ! number of evaluations of the cost function

#if defined MPI
  include "mpif.h"
  integer, save :: mpi_code
#endif

  ! Initialize parallel computation
#if defined MPI
  call mpi_init(mpi_code)
  call mpi_comm_size(mpi_comm_world,nproc,mpi_code)
  call mpi_comm_rank(mpi_comm_world,iproc,mpi_code)
#endif

  ! Read parameters from namelist
  namunit=10
  open(namunit,FILE='namelist',FORM='FORMATTED')
  read(namunit,experiment,iostat=ios)
  if (ios /= 0) stop 'error reading namelist: experiment'
  read(namunit,grid,iostat=ios)
  if (ios /= 0) stop 'error reading namelist: grid'
  read(namunit,proposal_par,iostat=ios)
  if (ios /= 0) stop 'error reading namelist: proposal_par'
  read(namunit,options,iostat=ios)
  if (ios /= 0) stop 'error reading namelist: options'
  read(namunit,physics,iostat=ios)
  if (ios /= 0) stop 'error reading namelist: physics'
  read(namunit,constraints,iostat=ios)
  if (ios /= 0) stop 'error reading namelist: constraints'
  close(namunit)

  ! Define grid
  nlon = nlongitude
  nlat = nlatitude
  lonmin = longitude_min
  lonmax = longitude_max
  latmin = latitude_min
  latmax = latitude_max
  call defgrid()

  ! Seed random number generator
  call kiss_reset()
  do jiproc = 0, iproc + seedloop
    zseed1 = kiss() ; zseed2 = kiss() ; zseed3 = kiss() ; zseed4 = kiss()
  enddo
  call kiss_seed( zseed1, zseed2, zseed3, zseed4 )

  ! Allocate physical arrays
  allocate(psi(nlon,nlat))
  allocate(u(nlon,nlat))
  allocate(v(nlon,nlat))
  allocate(unorm(nlon,nlat))
  allocate(zeta(nlon,nlat))
  allocate(advresidual(nlon,nlat))
  allocate(roughness(nlon,nlat))
  allocate(slope_penalty(nlon,nlat))

  if (ktend>0) then
    allocate(u0(nlon,nlat))
    allocate(v0(nlon,nlat))
    allocate(zeta0(nlon,nlat))
    if (tau>0.) then
      allocate(zetaforc(nlon,nlat))
    endif
  endif

  ! Problem definition
  ! ==================
  ! Definition of the grid locations (x,y) and
  ! domain decomposition for parallel computing
  if (iproc.eq.0) print *, 'Domain decomposition'
  allocate(indx(nlon,nlat)) ; indx = 0
  call domain_decomposition()

  allocate(tmp_state(n))
  tmp_state = 0.
  allocate(tmp_state2(n))
  tmp_state2 = 0.

  ! Diagnose total area of the domain
  !outfile=trim(expt)//'_area.nc'
  !call create_diag_file(outfile,x,y)
  !call write_diag_variable(outfile,'area',a)
  total_area = sum(a)

#if defined MPI
      call MPI_ALLREDUCE (MPI_IN_PLACE, total_area, 1, MPI_DOUBLE_PRECISION,  &
                    &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

  if (iproc.eq.0) print *, 'Total area:',total_area
  
  ! rescale s0 with total_area
  s0 = s0 * total_area


  allocate(proposal(n,m,s))

  if (read_inputs) then

    if (iproc.eq.0) print *, 'Read proposal distribution'

    infile=trim(expt)//'_proposal.nc'
    call read_ensemble(infile,proposal(:,:,1))

  else

    ! Generate proposal distribution (proposition directions)
    ! by sampling Gaussian random field with the specified spectrum 
    if (iproc.eq.0) print *, 'Generate proposal distribution'

    storfg_ylm_resolution=dlat/real(plegres,8) ! resolution for the Legendre polynomials
    do jm = 1, m
      call gen_field_2s(proposal(:,jm,1),x,y,pow_spectrum,0,lmax)
    enddo

    ! Write proposal distribution in output files (in NetCDF)
    outfile=trim(expt)//'_proposal.nc'
    call write_ensemble(outfile,x,y,proposal(:,:,1))

  endif

  ! Check proposal distribution
  if (ktend==0) then
    outfile=trim(expt)//'_proposal_diag.nc'
    !print *, 'minmaxprop',minval(proposal(:,1,1)),maxval(proposal(:,1,1))
    call state_diag(outfile,proposal(:,1,1))
  endif

  ! Construct multiple scale proposal
  ! =================================
  if (iproc.eq.0) print *, 'Construct multiple scale proposal'
  smax = maxval(scales(:))  ! maximum degree of the spherical harmonics
  allocate(spct(0:smax,-smax:smax))

  external_vector_decomposition=.TRUE. ! parallelization is ruled from outside sphylm
  call init_ylm( smax, 0, latitude_min, latitude_max, dlat/real(plegres,8) )

  do jm=1,m  ! loop on ensemble members

    tmp_state(:) = proposal(:,jm,1) * a(:)
    call proj_ylm( spct, tmp_state, x, y )
#if defined MPI
    call mpi_allreduce (MPI_IN_PLACE, spct, (smax+1)*(2*smax+1),  &
     &                  MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm_world,mpi_code)
#endif

    do js=2,s ! loop on additional scales to include in the proposal
      proposal(:,jm,js) = 0.
      call back_ylm_loc( spct, proposal(:,jm,js), x, y, 0, scales(js) )
    enddo

  enddo

  ! Renormalize std of all additional scales to 1
  do js=2,s
    do jn=1,n
      std = 0. ; mean = 0. ! temporary storage to store ensemble mean
      do jm=1,m
        misfit = proposal(jn,jm,js) - mean
        mean = mean + misfit / jm
        misfit2 = proposal(jn,jm,js) - mean
        std = std + misfit * misfit2
      enddo
      std = SQRT( std / (m-1) )
      proposal(jn,:,js) = proposal(jn,:,js) / std
    enddo
  enddo

  ! Write files (in NetCDF)
  outfile=trim(expt)//'_proposal_js=2.nc'
  call write_ensemble(outfile,x,y,proposal(:,:,2))

  ! Check Schur products
  allocate(sample(n,m))
  allocate( tmp_sample(SUM(multiplicity(:))) )
  do jm=1,m
    CALL newproduct( tmp_state, proposal, multiplicity, tmp_sample )
    sample(:,jm) = tmp_state(:)
  enddo
  outfile=trim(expt)//'_multiple_products.nc'
  call write_ensemble(outfile,x,y,sample(:,:))
  deallocate(sample)
  deallocate( tmp_sample )

  ! Draw sample using MCMC algorithm
  ! ================================
  if (iproc.eq.0) print *, 'Draw sample using MCMC algorithm'

  if (perform_sampling) then

    ! Set number of iterations between convergence checks
    mcmc_convergence_check = convergence_check
    mcmc_control_print=1
    mcmc_proposal=.not.prior
    mcmc_proposal_std=proposal_std

    ! Allocate sample array
    if (mup.ne.1) STOP 'mup>1 not implemented'
    allocate(sample(n,mup))

    ! Time iteration
    do kt=ktini,ktend

      if (iproc.eq.0) print *, 'Running time step kt=',kt

      ! Initialize iteration
      if (kt==0) then
        if (proposal_init) then
          do jm=1,mup
            sample(:,jm)=proposal(:,1+MOD(jm-1,m),1)
          enddo
        else
          sample(:,:)=0.
        endif
      elseif (kt==ktini) then
        if (tau>0.) then
          ! read forcing file
          infile=trim(expt)//'_forcing.nc'
          call read_ensemble(infile,sample(:,:))
          if (iproc.eq.0) print *, 'sample:',minval(sample),maxval(sample)
          ! compute forcing
          call recompose_state(sample(:,1),psi)
          call velocity(psi,u0,v0)
          call vorticity(u0,v0,zetaforc)
          if (iproc.eq.0) print *, 'forcing:',tau, minval(zetaforc),maxval(zetaforc)
          ! forcing diagnostic
          outfile=trim(expt)//'_forcing_diag.nc'
          call state_diag(outfile,sample(:,1))
        endif
        if ((tau>0.).and.(kt==1))  then
          ! start from a zero initial condition
          sample(:,:)=0. ; u0(:,:)=0. ; v0(:,:)=0. ; zeta0(:,:)=0.
        else
          ! read restart file
          write(tagkt,'(i4.4)') kt-1
          infile=trim(expt)//'_'//tagkt//'_sample.nc'
          call read_ensemble(infile,sample(:,:))
          ! initialize u0, v0, zeta0
          call recompose_state(sample(:,1),psi)
          call velocity(psi,u0,v0)
          call vorticity(u0,v0,zeta0)
          ! initial diagnostic
          write(tagkt,'(i4.4)') kt
          outfile=trim(expt)//'_'//tagkt//'_ini_diag.nc'
          call state_diag(outfile,sample(:,1))
        endif
      else
        u0 = u ; v0 = v ; zeta0 = zeta
      endif

      ! Call to MCMC iteration to perform sampling
      if (iproc.eq.0) print *, 'Starting MCMC iteration kt=',kt
      mcmc_index=1
      mcmc_zero_start = .FALSE.
      call mcmc_iteration( mcmc_maxiter, sample, proposal, multiplicity, cost_jo, convergence_test )

      ! Diagnostics of the sample
      if (iproc.eq.0) print *, 'Number of evaluations of the cost function:', njo

      ! Normalize initial condition
      if (kt==0) then
        ! compute u, v, zeta
        call recompose_state(sample(:,1),psi)
        call velocity(psi,u,v)
        call vorticity(u,v,zeta)
        ! compute root mean square zeta
        zetarms = sqrt(sum(zeta*zeta)/(nlon*nlat))
        if (iproc.eq.0) print *, 'Initial renormalization by:',zetarms
        ! renormamize so that root mean square zeta is equal to 1
        sample(:,:) = sample(:,:) / zetarms
      endif

      ! Write sample (in NetCDF)
      write(tagkt,'(i4.4)') kt
      outfile=trim(expt)//'_'//tagkt//'_sample.nc'
      call write_ensemble(outfile,x,y,sample(:,:))

      ! Write physical diagnostic of final sample (in NetCDF)
      write(tagkt,'(i4.4)') kt
      outfile=trim(expt)//'_'//tagkt//'_diag.nc'
      call state_diag(outfile,sample(:,1))

    enddo

  endif

contains

  ! ======================================================================

  ! Callback routine to check the convergence of iterations
  ! Intermediate diagnostics can also be computed here
  function convergence_test(upens,upxens)
  use ensdam_score_crps
  implicit none
  real(kind=8), dimension(:,:), intent(in) :: upens
  real(kind=8), dimension(:,:), intent(in), optional :: upxens
  logical :: convergence_test

  character(len=6) :: tagiter
  real(kind=8) :: cost
  logical :: diag

  convergence_test = .FALSE.

  ! Initialize previous score at first call for this time step
  if (kt.ne.ktsave) then
    ktsave = kt
    iteration_index = convergence_check
    score_prev=tiny(score_prev)
    if (kt>ktini) close(20)
    write(tagkt,'(i4.4)') kt
    outfile=trim(expt)//'_'//tagkt//'_constraint.txt'
    open(unit=20,file=outfile,form='FORMATTED')
  else
    iteration_index = iteration_index + convergence_check
    call flush(20)
  endif

  ! Compute convergence score
  if (check_convergence) then
    if (iproc.eq.0) print *, 'SCORE:',iteration_index,score
    convergence_test = abs( score - score_prev ) / score_prev < convergence_epsilon
    score_prev=score
    if (convergence_test) then
      if (iproc.eq.0) print *, 'Convergence reached !'
    endif
  endif

  ! Write physical diagnostic of first mmember of current sample (in NetCDF)
  if (iteration_diag) then
    write(tagiter,'(i6.6)') iteration_index
    outfile=trim(expt)//'_'//tagkt//'_'//tagiter//'_diag.nc'
    call state_diag(outfile,upens(:,1))
  endif

  ! Diagnostic of constraints (first member)
  if (constraint_diag) then
    diag = .TRUE.
    cost = constraint(upens(:,1),diag)
  endif

  end function convergence_test

  ! ======================================================================

  ! Callback routine to compute cost function
  ! Jo = -log p(yo|hx), computed using state vector as argument
  function cost_jo(state)
  implicit none
  real(kind=8), dimension(:), intent(in) :: state
  real(kind=8) :: cost_jo
  
  real(kind=8), dimension(:), allocatable :: tmpstate
  logical :: diag

  cost_jo = 0.

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE, cost_jo, 1, MPI_DOUBLE_PRECISION,  &
                        &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

  diag = .FALSE.
  if (apply_constraint) cost_jo = cost_jo + constraint(state,diag)

  njo = njo + 1

  end function cost_jo

  ! ======================================================================

  ! Dynamical constraint on state vector
  function constraint(state,diag)
  implicit none
  real(kind=8), dimension(:), intent(in) :: state
  logical, intent(in) :: diag
  real(kind=8) :: constraint

  real(kind=8) :: constraint_velocity, constraint_advection, constraint_diffusion 
  real(kind=8) :: constraint_vorticity, constraint_turbulent_dissipation
  real(kind=8) :: constraint_boundary, area_boundary
  real(kind=8) :: maxu, maxv, zeta0ms, zeta1ms

  integer :: i
  logical :: zero_velocity, compute_velocity, compute_vorticity

  constraint_velocity = 0.
  constraint_vorticity = 0.
  constraint_advection = 0.
  constraint_boundary = 0.
  constraint_diffusion = 0.
  constraint_turbulent_dissipation = 0.


  compute_vorticity = constrain_advection .or. constrain_diffusion .or. constrain_turbulent_dissipation
  compute_vorticity = compute_vorticity .or. constrain_vorticity
  compute_vorticity = compute_vorticity.and.(kt>0)

  compute_velocity = compute_vorticity .or. constrain_velocity

  ! recompose full state
  call recompose_state(state,psi)

  ! compute velocity
  if (compute_velocity.or.diag) then
    call velocity(psi,u,v)
    maxu = maxval(abs(u)) ; maxv = maxval(abs(v))
    zero_velocity = (maxu.eq.0.).and.(maxv.eq.0.)
  endif

  ! constrain velocity
  if (constrain_velocity.or.diag) then
    call decompose_state(u,tmp_state)
    constraint_velocity = sum(tmp_state*tmp_state*a)
    call decompose_state(v,tmp_state)
    constraint_velocity = constraint_velocity + sum(tmp_state*tmp_state*a)

#if defined MPI
    call MPI_ALLREDUCE (MPI_IN_PLACE, constraint_velocity, 1, MPI_DOUBLE_PRECISION,  &
                        &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

    uvrms =sqrt(constraint_velocity)
    if (iproc.eq.0) print *, 'uvrms',uvrms

    if (zero_velocity) then
     constraint_velocity = huge(constraint_velocity)
    else
     constraint_velocity = ( log ( constraint_velocity ) / velvar ) **2
    endif

  endif

  ! compute vorticity
  if (compute_vorticity.or.diag) call vorticity(u,v,zeta)

  ! constrain vorticity
  if ((constrain_vorticity.or.diag).and.(kt>0)) then

    call decompose_state(zeta,tmp_state)
    zeta1ms = sum(tmp_state*tmp_state*a)

#if defined MPI
    call MPI_ALLREDUCE (MPI_IN_PLACE, zeta1ms, 1, MPI_DOUBLE_PRECISION,  &
                        &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

    constraint_vorticity = - mu * mu * zeta1ms / total_area

  endif

  ! compute advection constraint penalty
  if ((constrain_advection.or.diag).and.(kt>0)) then

    if (zero_velocity) then
      constraint_advection = huge(constraint_advection)
      constraint_boundary = huge(constraint_boundary)
    else
      ! compute material derivative of potential vorticity
      call advection(advresidual,zeta0,u0,v0,zeta,u,v,rossby,dt)

      ! apply limitation in unconstrained boundary region
      call decompose_state(advresidual,tmp_state)
      call decompose_state(zeta,tmp_state2)
      constraint_boundary = sum(tmp_state2*tmp_state2*a,mask=tmp_state.eq.0._8)
      area_boundary = sum(a,mask=tmp_state.eq.0._8)

      ! modification of boundary limitation in presence of diffusion
      if (constrain_diffusion) then
        call decompose_state(zeta0,tmp_state)
        zeta0ms = sum(tmp_state*tmp_state*a)
      endif

      ! apply forcing
      if (tau>0.) then
        advresidual = advresidual - zetaforc / tau
      endif

      ! compute advection constraint
      call decompose_state(advresidual,tmp_state)
      constraint_advection = sum(tmp_state*tmp_state*a)

#if defined MPI
      call MPI_ALLREDUCE (MPI_IN_PLACE, constraint_advection, 1, MPI_DOUBLE_PRECISION,  &
                    &     MPI_SUM,mpi_comm_world,mpi_code)
      call MPI_ALLREDUCE (MPI_IN_PLACE, constraint_boundary, 1, MPI_DOUBLE_PRECISION,  &
                    &     MPI_SUM,mpi_comm_world,mpi_code)
      call MPI_ALLREDUCE (MPI_IN_PLACE, area_boundary, 1, MPI_DOUBLE_PRECISION,  &
                    &     MPI_SUM,mpi_comm_world,mpi_code)
      if (constrain_diffusion) then
        call MPI_ALLREDUCE (MPI_IN_PLACE, zeta0ms, 1, MPI_DOUBLE_PRECISION,  &
                           &     MPI_SUM,mpi_comm_world,mpi_code)
      endif
#endif

      constraint_advection = constraint_advection / s0
      constraint_boundary = constraint_boundary / area_boundary

      if ((constrain_diffusion).and.(zeta0ms>0.)) then
        zeta0ms = zeta0ms / total_area
        constraint_boundary = constraint_boundary / min(zeta0ms,1._8)
        constraint_boundary = constraint_boundary * 250
      endif

    endif

  endif

  ! compute diffusion constraint penalty
  if (constrain_diffusion.or.diag) then

    if (zero_velocity) then
      constraint_diffusion = huge(constraint_diffusion)
    else
      call diffusion(roughness,zeta)
      call decompose_state(roughness,tmp_state)
      constraint_diffusion = sum(tmp_state*a)*(nu*nu)

#if defined MPI
      call MPI_ALLREDUCE (MPI_IN_PLACE, constraint_diffusion, 1, MPI_DOUBLE_PRECISION,  &
                    &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

      constraint_diffusion = constraint_diffusion / s0
    endif

  endif

  ! compute turbulent dissipation constraint penalty
  if (constrain_turbulent_dissipation.or.diag) then

    if (zero_velocity) then
      constraint_turbulent_dissipation = huge(constraint_turbulent_dissipation)
    else
      call turbulent_dissipation(slope_penalty,zeta)
      call decompose_state(slope_penalty,tmp_state)
      constraint_turbulent_dissipation = sum(tmp_state*a)*(gamma*gamma)*(nu*nu)

#if defined MPI
      call MPI_ALLREDUCE (MPI_IN_PLACE, constraint_turbulent_dissipation, 1, MPI_DOUBLE_PRECISION,  &
                    &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

      constraint_turbulent_dissipation = constraint_turbulent_dissipation / s0
    endif

  endif

  constraint = 0.
  if (constrain_velocity) constraint = constraint + constraint_velocity
  if (constrain_vorticity.and.(kt>0)) constraint = constraint + constraint_vorticity
  if (constrain_advection.and.(kt>0)) constraint = constraint + constraint_advection
  if (constrain_advection.and.(kt>0)) constraint = constraint + constraint_boundary
  if (constrain_diffusion.and.(kt>0)) constraint = constraint + constraint_diffusion
  if (constrain_turbulent_dissipation.and.(kt>0)) constraint = constraint + constraint_turbulent_dissipation

  if ( (iproc.eq.0) .and. diag ) then
    write(20,'(i6.6,6e18.10)') iteration_index, constraint, &
         &                     constraint_vorticity, &
         &                     constraint_boundary, &
         &                     constraint_advection, &
         &                     constraint_diffusion, &
         &                     constraint_turbulent_dissipation
  endif

  constraint = constraint * (nlon*nlat) * 0.5

  end function constraint

  ! ======================================================================

  ! Power spectrum in the basis of the spherical harmonics
  ! (as a function of the degree l)
  function pow_spectrum(l,m)
  implicit none
  integer, intent(in) :: l,m
  real(kind=8) :: pow_spectrum

  integer :: ll, mm
  real(kind=8) :: norm, rm, rl

  ! Power spectrum
  rl = (real(l,8)-lm)/lc
  pow_spectrum = 1. / ( 1. + rl**lexp )

  ! Normalize the spectrum
  norm = 0.
  do ll=0,lmax
    rl = (real(ll,8)-lm)/lc
    norm = norm + 1. / ( 1. + rl**lexp )
  enddo
  pow_spectrum = pow_spectrum / norm

  ! Scale to account for the multiplicity of each degree
  pow_spectrum = pow_spectrum / ( 1. + 2. * real(l,8) )

  end function pow_spectrum

  ! ======================================================================

  ! Domain decomposition for parallel computing
  subroutine domain_decomposition()
  use ensdam_spharea
  implicit none

  real(kind=8), dimension(:,:), allocatable :: lon, lat, area
  integer :: i, j, k

  ! allocate grid on the full domain
  allocate(lon(0:nlon,0:nlat),lat(0:nlon,0:nlat),area(nlon,nlat))

  ! location of mesh vertices (F points)
  do i=0,nlon
    lon(i,0) = lonmin + i * dlon - dlon / 2
  enddo
  do j=0,nlat
    lat(0,j) = latmin + j * dlat - dlat / 2
  enddo
  lat(0,0) = max(lat(0,0),-90.)
  lat(0,nlat) = min(lat(0,nlat),90.)

  do i=0,nlon
    lon(i,1:nlat) = lon(i,0)
  enddo
  do j=0,nlat
    lat(1:nlon,j) = lat(0,j)
  enddo

  ! compute area associated to each grid mesh
  call mesh_area (lon,lat,area)

  ! location (T points) of state variables (0 indices are not used anymore)
  lon(:,:) = lon(:,:) - dlon / 2
  lat(:,:) = lat(:,:) - dlat / 2

  ! determine size of subdomain for current processor
  n = nlon * nlat / nproc
  if (iproc.lt.mod(nlon*nlat,nproc)) n = n + 1

  ! allocate local grid arrays
  allocate(x(n),y(n),a(n))

  ! provide a subdomain to each processor
  k=0
  do j=1,nlat
  do i=1,nlon
    if (mod(k-iproc,nproc).eq.0) then
      x(1+k/nproc) = lon(i,j)
      y(1+k/nproc) = lat(i,j)
      a(1+k/nproc) = area(i,j)
      indx(i,j) = k+1
    endif
    k = k + 1
  enddo
  enddo

  ! deallocate global grid
  deallocate(lon,lat)

  ! build global index array
#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE,indx,nlon*nlat,MPI_INTEGER, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

  end subroutine domain_decomposition

  ! ======================================================================

  subroutine recompose_state(state,fullstate)

  implicit none
  real(kind=8), dimension(:), intent(in) :: state
  real(kind=8), dimension(:,:), intent(out) :: fullstate

  integer :: i, j, k

  fullstate = 0.

  ! merge processors contributions
  k=0
  do j=1,nlat
  do i=1,nlon
    if (mod(k-iproc,nproc).eq.0) then
      fullstate(i,j) = state(1+k/nproc)
    endif
    k = k + 1
  enddo
  enddo

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE,fullstate,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

  end subroutine recompose_state

  ! ======================================================================

  subroutine decompose_state(fullstate,state)

  implicit none
  real(kind=8), dimension(:,:), intent(in) :: fullstate
  real(kind=8), dimension(:), intent(out) :: state

  integer :: i, j, k

  state = 0.

  ! merge processors contributions
  k=0
  do j=1,nlat
  do i=1,nlon
    if (mod(k-iproc,nproc).eq.0) then
      state(1+k/nproc) = fullstate(i,j)
    endif
    k = k + 1
  enddo
  enddo

  end subroutine decompose_state

  ! ======================================================================
  ! Create diag file
  subroutine create_diag_file(filename,x,y)
  use netcdf
  implicit none

  character(len=256), intent(in) :: filename
  real(kind=8), dimension(:), intent(in) :: x, y

  real(kind=8), dimension(:,:), allocatable :: lon, lat

  integer :: i, j, k
  integer :: is, ncid, idx, idy, idt, idlon, idlat

  ! allocate full grid and full state
  allocate(lon(nlon,nlat),lat(nlon,nlat))
  lon = 0. ; lat = 0.

  ! merge processors contributions
  k=0
  do j=1,nlat
  do i=1,nlon
    if (mod(k-iproc,nproc).eq.0) then
      lon(i,j) = x(1+k/nproc)
      lat(i,j) = y(1+k/nproc)
    endif
    k = k + 1
  enddo
  enddo

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE,lon, nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
  call MPI_ALLREDUCE (MPI_IN_PLACE,lat, nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

  ! write variables in output file
  if (iproc.eq.0) then
    is = NF90_CREATE(filename,NF90_CLOBBER,ncid)
    is = NF90_DEF_DIM(ncid,'lon',nlon,idx)
    is = NF90_DEF_DIM(ncid,'lat',nlat,idy)
    is = NF90_DEF_DIM(ncid,'time',nf90_unlimited,idt)
    is = NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,(/idx,idy/),idlon)
    is = NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,(/idx,idy/),idlat)
    is = NF90_ENDDEF(ncid)
    is = NF90_PUT_VAR(ncid,idlon,lon)
    is = NF90_PUT_VAR(ncid,idlat,lat)
    is = NF90_CLOSE(ncid)
  endif

  deallocate(lon,lat)

  end subroutine create_diag_file

 ! ======================================================================

  ! Write diagnostic variable in output file
  subroutine write_diag_variable(filename,varname,state)
  use netcdf
  implicit none

  character(len=*), intent(in) :: filename, varname
  real(kind=8), dimension(:), intent(in) :: state

  real(kind=8), dimension(:,:), allocatable :: fullstate
  integer :: i, j, k

  integer :: is, ncid, idx, idy, idt, idv

  ! allocate full grid and full state
  allocate(fullstate(nlon,nlat))
  fullstate = 0.

  ! merge processors contributions
  k=0
  do j=1,nlat
  do i=1,nlon
    if (mod(k-iproc,nproc).eq.0) then
      fullstate(i,j) = state(1+k/nproc)
    endif
    k = k + 1
  enddo
  enddo

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE,fullstate, nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

  ! write full variable in output file
  if (iproc.eq.0) then
    is = NF90_OPEN(filename,NF90_WRITE,ncid)
    is = NF90_INQ_DIMID(ncid,'lon',idx)
    is = NF90_INQ_DIMID(ncid,'lat',idy)
    is = NF90_INQ_DIMID(ncid,'time',idt)
    is = NF90_REDEF(ncid)
    is = NF90_DEF_VAR(ncid,varname,NF90_FLOAT,(/idx,idy,idt/),idv)
    is = NF90_ENDDEF(ncid)
    is = NF90_PUT_VAR(ncid,idv,fullstate)
    is = NF90_CLOSE(ncid)
  endif

  deallocate(fullstate)

  end subroutine write_diag_variable

  ! ======================================================================

  ! Write ensemble in file
  subroutine write_ensemble(filename,x,y,ensemble)
  use netcdf
  implicit none

  character(len=256), intent(in) :: filename
  real(kind=8), dimension(:), intent(in) :: x, y
  real(kind=8), dimension(:,:), intent(in) :: ensemble

  real(kind=8), dimension(:,:), allocatable :: lon, lat, fullstate
  integer :: i, j, k, imem, nmem

  integer :: is, ncid, idx, idy, idm, idv, idlon, idlat
  integer, dimension(3) :: vstart

  nmem = SIZE(ensemble,2)

  if (iproc.eq.0) then
    is = NF90_CREATE(filename,NF90_CLOBBER,ncid)
    is = NF90_DEF_DIM(ncid,'lon',nlon,idx)
    is = NF90_DEF_DIM(ncid,'lat',nlat,idy)
    is = NF90_DEF_DIM(ncid,'member',nmem,idm)
    is = NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,(/idx,idy/),idlon)
    is = NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,(/idx,idy/),idlat)
    is = NF90_DEF_VAR(ncid,'ensemble',NF90_FLOAT,(/idx,idy,idm/),idv)
    is = NF90_ENDDEF(ncid)
  endif

  ! allocate full grid
  allocate(lon(nlon,nlat),lat(nlon,nlat))
  lon = 0. ; lat = 0.

  ! merge location arrays
  k=0
  do j=1,nlat
  do i=1,nlon
    if (mod(k-iproc,nproc).eq.0) then
      lon(i,j) = x(1+k/nproc)
      lat(i,j) = y(1+k/nproc)
    endif
    k = k + 1
  enddo
  enddo

#if defined MPI
    call MPI_ALLREDUCE (MPI_IN_PLACE,lon, nlon*nlat,MPI_DOUBLE_PRECISION, &
       &                MPI_SUM,mpi_comm_world,mpi_code)
    call MPI_ALLREDUCE (MPI_IN_PLACE,lat, nlon*nlat,MPI_DOUBLE_PRECISION, &
       &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

  if (iproc.eq.0) then
    is = NF90_PUT_VAR(ncid,idlon,lon)
    is = NF90_PUT_VAR(ncid,idlat,lat)
  endif

  ! allocate full state
  deallocate(lon,lat)
  allocate(fullstate(nlon,nlat))

  DO imem=1,nmem

    fullstate = 0.

    ! merge processors contributions
    k=0
    do j=1,nlat
    do i=1,nlon
      if (mod(k-iproc,nproc).eq.0) then
        fullstate(i,j) = ensemble(1+k/nproc,imem)
      endif
      k = k + 1
    enddo
    enddo

#if defined MPI
    call MPI_ALLREDUCE (MPI_IN_PLACE,fullstate, nlon*nlat,MPI_DOUBLE_PRECISION, &
       &                MPI_SUM,mpi_comm_world,mpi_code)
#endif

    ! write full variable in output file
    if (iproc.eq.0) then
      vstart=1 ; vstart(3)=imem
      is = NF90_PUT_VAR(ncid,idv,fullstate,start=vstart)
    endif

  enddo

  if (iproc.eq.0) then
    is = NF90_CLOSE(ncid)
  endif

  deallocate(fullstate)

  end subroutine write_ensemble

  ! ======================================================================

  ! Read ensemble from file
  subroutine read_ensemble(filename,ensemble)
  use netcdf
  implicit none

  character(len=256), intent(in) :: filename
  real(kind=8), dimension(:,:), intent(out) :: ensemble

  real(kind=8), dimension(:,:), allocatable :: fullstate
  integer :: i, j, k, imem, nmem

  integer :: is, ncid, idv
  integer, dimension(3) :: vstart

  nmem = SIZE(ensemble,2)

  if (iproc.eq.0) print *, 'Reading file:',trim(filename)
  is = NF90_OPEN(filename,NF90_NOWRITE,ncid)
  if (is.ne.0) stop 'error opening NetCDF file'
  is = NF90_INQ_VARID(ncid,'ensemble',idv)

  allocate(fullstate(nlon,nlat))

  DO imem=1,nmem

    ! read full variable in output file
    vstart=1 ; vstart(3)=imem
    is = NF90_GET_VAR(ncid,idv,fullstate,start=vstart)

    ! extract processor contributions
    k=0
    do j=1,nlat
    do i=1,nlon
      if (mod(k-iproc,nproc).eq.0) then
        ensemble(1+k/nproc,imem) = fullstate(i,j)
      endif
      k = k + 1
    enddo
    enddo

  enddo

  is = NF90_CLOSE(ncid)

  deallocate(fullstate)

  end subroutine read_ensemble

  ! ======================================================================

  subroutine state_diag(filename,state)

  implicit none
  character(len=256), intent(in) :: filename
  real(kind=8), dimension(:), intent(in) :: state

  character(len=256) :: outdiag

  outdiag=filename

  call create_diag_file(outdiag,x,y)
  call write_diag_variable(outdiag,'psi',state)

  call recompose_state(state(:),psi)
  call velocity(psi,u,v)

  call decompose_state(u,tmp_state)
  uvrms = sum(tmp_state*tmp_state*a)
  call decompose_state(v,tmp_state)
  uvrms = uvrms + sum(tmp_state*tmp_state*a)

#if defined MPI
  call MPI_ALLREDUCE (MPI_IN_PLACE, uvrms, 1, MPI_DOUBLE_PRECISION,  &
                        &     MPI_SUM,mpi_comm_world,mpi_code)
#endif

  uvrms =sqrt(uvrms)
  if (iproc.eq.0) print *, 'uvrms',uvrms

  call decompose_state(u,tmp_state)
  call write_diag_variable(outdiag,'u',tmp_state)
  call decompose_state(v,tmp_state)
  call write_diag_variable(outdiag,'v',tmp_state)

  call velocity_norm(u,v,unorm)

  call decompose_state(unorm,tmp_state)
  call write_diag_variable(outdiag,'unorm',tmp_state)

  call vorticity(u,v,zeta)

  call decompose_state(zeta,tmp_state)
  call write_diag_variable(outdiag,'zeta',tmp_state)

  call diffusion(roughness,zeta)

  call decompose_state(roughness,tmp_state)
  call write_diag_variable(outdiag,'roughness',tmp_state)

  call turbulent_dissipation(slope_penalty,zeta)

  call decompose_state(slope_penalty,tmp_state)
  call write_diag_variable(outdiag,'slope_penalty',tmp_state)

  if (kt>0) then
    call advection(advresidual,zeta0,u0,v0,zeta,u,v,rossby,dt)

    call decompose_state(advresidual,tmp_state)
    call write_diag_variable(outdiag,'advresidual',tmp_state)
  endif

  end subroutine state_diag

  ! ======================================================================

end program flow_sampler
