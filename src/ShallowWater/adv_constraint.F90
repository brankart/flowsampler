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
!                        MODULE ADV_CONSTRAINT
!
!---------------------------------------------------------------------
! Dynamic operators on the sphere
! by Jean-Michel Brankart, June 2019
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! constraint_init : constraint initialization
! constraint_cost : evaluate constraint cost function
! ----------------------------------------------------------------------
MODULE flowsampler_adv_constraint
  use flowsampler_grid
  use flowsampler_units
  IMPLICIT NONE
  PRIVATE

  PUBLIC constraint_init, constraint_cost

  ! Private arrays with data fields
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: u0, v0, zeta0
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: u1, v1, zeta1
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: xblock0, xblock1

  ! Private arrays with metric/scaling factors
  REAL(KIND=8), DIMENSION(:), SAVE, ALLOCATABLE :: invdx, invdy
  REAL(KIND=8), DIMENSION(:), SAVE, ALLOCATABLE :: invdxfac, invdyfac, zetafac
  REAL(KIND=8), DIMENSION(:), SAVE, ALLOCATABLE :: dtfacx, dtfacy
  REAL(KIND=8), SAVE :: zeta_planetary

#if defined MPI
  ! Definitions for MPI
  include "mpif.h"
  INTEGER, PUBLIC, SAVE  :: mpi_comm_flow_sampler=mpi_comm_world  ! definition of module global communicator
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpi_status_msg
  INTEGER, save :: mpi_code
#endif

  ! Definitions for MPI
  INTEGER, save :: iproc=0, nproc=1
  INTEGER, save :: nproct, iproct, iproc_prev, iproc_next
  INTEGER, save :: mpi_comm_timestep
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ncommt  ! communicator for every timestep


  CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE constraint_init(nt)
!----------------------------------------------------------------------
! ** Purpose : constraint initialization
!----------------------------------------------------------------------
    INTEGER, INTENT( in ) :: nt

    INTEGER :: allocok, iiproct
    INTEGER :: igrp_global
    INTEGER, DIMENSION(:), ALLOCATABLE ::   irankt ! list of processors for current timestep
    INTEGER, DIMENSION(:), ALLOCATABLE ::   igrpt  ! group ID for every timestep
    REAL(KIND=8) :: sinlat, coslat, dlonratio, dlatratio, gravityoverf

#if defined MPI
    CALL mpi_comm_size(mpi_comm_flow_sampler,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flow_sampler,iproc,mpi_code)
    IF (MOD(nproc,nt).NE.0) STOP 'Error in flowsampler/constraint_init: bad nproc'
    nproct = nproc/nt
    iproct = MOD(iproc,nproct)

    ALLOCATE(ncommt(nt),stat=allocok)     ! communicator for timesteps
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    ALLOCATE(irankt(nproct),stat=allocok) ! temporary array
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    ALLOCATE(igrpt(nt),stat=allocok)      ! temporary array
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'

    ! Create global group
    CALL mpi_comm_group(mpi_comm_flow_sampler,igrp_global,mpi_code)

    ! Create one communicator for each timestep
    DO jt = 1, nt
      irankt(:) = (/ (iiproct, iiproct=0,nproct-1) /)
      irankt(:) = irankt(:) + ( jt - 1 ) * nproct
      ! Create the group for current timestep from the global group
      CALL mpi_group_incl( igrp_global, nproct, irankt, igrpt(jt), mpi_code )
      ! Create communicator for jt
      CALL mpi_comm_create( mpi_comm_flow_sampler, igrpt(jt), ncommt(jt), mpi_code )
    ENDDO

    ! Deallocate arrays
    DEALLOCATE( irankt , igrpt )

    ! Identify communicator for current time step
    mpi_comm_timestep = ncommt(1+iproc/nproct)

    ! Identify corresponding processor in previous and next timestep
    ! (index in global communicator)
    iproc_prev = iproc-nproct
    iproc_next = iproc+nproct
#else
    STOP 'Error in flowsampler/constraint_init: MPI is required'
#endif

    ! block size
    nk = nlon * nlat / nproct
    ! block bounds
    k0 = iproct * nk
    k1 = k0 + nk - 1

    ! allocation of blocks of data to be computed by current processor
    allocate(xblock0(nk),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(xblock1(nk),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'

    ! allocation of 2D arrays for one time step
    allocate(u0(nlon,nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(v0(nlon,nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(zeta0(nlon,nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(u1(nlon,nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(v1(nlon,nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(zeta1(nlon,nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'

    ! allocation of metric arrays
    allocate(invdx(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(invdy(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(invdxfac(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(invdyfac(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(zetafac(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(dtfacx(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(dtfacy(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'

    ! compute metric/rescaling arrays
    DO j=1,nlat
      coslat = coslatitude(j)
      sinlat = sinlatitude(j)
      CALL distance_ratio(j,dlonratio,dlatratio)
      CALL gravityoverf_ratio(j,gravityoverf)
      invdx(j) = invdlon / ( dlonratio * coslat )
      invdy(j) = invdlat / dlatratio
      invdxfac(j) = gravityoverf * invdlon / ( dlonratio * coslat )
      invdyfac(j) = gravityoverf * invdlat / dlatratio
      zetafac(j) = 0.5_8 * sinlat / ( dlonratio * coslat )
      dtfacx(j) = dt * 0.5_8 / (dlonratio*deg2rad)
      dtfacy(j) = dt * 0.5_8 / (dlatratio*deg2rad)
    ENDDO

    ! Compute constant factors
    zeta_planetary = 2._8 * earthvorticity

    END SUBROUTINE constraint_init
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    FUNCTION constraint_cost(adt)
!----------------------------------------------------------------------
! ** Purpose : evaluate constraint cost function
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: adt
    REAL(KIND=8) :: constraint_cost

    INTEGER :: i,j,k
    REAL(KIND=8) :: zeta1, zeta2, zeta3, dx, dy
    LOGICAL :: inmask

    ! Compute zonal velocity for this block of current time step
    DO k=k0,k1
      i = 1 + MOD(k,nlon)
      j = 1 + k/nlon
      kk = k - k0 + 1

      inmask = j.LT.nlat
      IF (inmask) THEN
        xblock1(kk) = - ( adt(i,j+1) - adt(i,j) ) * invdyfac(j)
      ELSE
        xblock1(kk) = 0.
      ENDIF
    ENDDO

    ! Get zonal velocity for this block of previous time step
    CALL get_prev_timestep(xblock0,xblock1)
    ! Rebuild previous and current timesteps from individual blocks
    CALL rebuild_timestep(xblock0,u0)
    CALL rebuild_timestep(xblock1,u1)

    ! Compute meridional velocity for this block of current time step
    DO k=k0,k1
      i = 1 + MOD(k,nlon)
      j = 1 + k/nlon
      kk = k - k0 + 1

      inmask = i.LT.nlon
      IF (inmask) THEN
        xblock1(kk) = ( psi(i+1,j) - psi(i,j) ) * invdxfac(j)
      ELSE
        xblock1(kk) = 0.
      ENDIF
    ENDDO

    ! Get meridonal velocity for this block of previous time step
    CALL get_prev_timestep(xblock0,xblock1)
    ! Rebuild previous and current timesteps from individual blocks
    CALL rebuild_timestep(xblock0,v0)
    CALL rebuild_timestep(xblock1,v1)

    ! Compute relative vorticity for this block of current time step
    DO k=k0,k1
      i = 1 + MOD(k,nlon)
      j = 1 + k/nlon
      kk = k - k0 + 1

      inmask = (i.GT.1).AND.(i.LT.nlon).AND.(j.GT.1).AND.(j.LT.nlat)
      IF (inmask) THEN
        zeta1 = ( u(i,j) + u(i,j-1) ) * zetafac(j)
        zeta2 = - ( u(i,j) - u(i,j-1) ) * invdy(j)
        zeta3 = ( v(i,j) - v(i-1,j) ) * zetafac(j)
        xblock1(kk)= zeta1 + zeta2 + zeta3
      ELSE
        xblock1(kk) = 0.
      ENDIF
    ENDDO

    ! Get relative vorticity for this block of previous time step
    CALL get_prev_timestep(xblock0,xblock1)
    ! Rebuild previous and current timesteps from individual blocks
    CALL rebuild_timestep(xblock0,zeta0)
    CALL rebuild_timestep(xblock1,zeta1)

    ! Compute advection residual for this block of current time step
    DO k=k0,k1
      i = 1 + MOD(k,nlon)
      j = 1 + k/nlon

      inmask = (i.GT.1).AND.(i.LT.nlon).AND.(j.GT.1).AND.(j.LT.nlat)
      IF (inmask) THEN
        ! advect potential vorticity from past situation
        dx = - u0(i,j) * dtfacx(j)
        dy = - v0(i,j) * dtfacy(j)
        call get_location_deg(i,j,dx,dy,lon,lat)
        call grid_interp(zeta0,lon,lat,zeta_past)
        IF (zeta_past.NE.spval) zeta_past = zeta_past + zeta_planetary * SIN(lat*deg2rad)
        ! advect potential vorticity from past situation
        dx = u1(i,j) * dtfacx(j)
        dy = v1(i,j) * dtfacy(j)
        call get_location_deg(i,j,dx,dy,lon,lat)
        call grid_interp(zeta1,lon,lat,zeta_future)
        IF (zeta_future.NE.spval) zeta_future = zeta_future + zeta_planetary * SIN(lat*deg2rad)

        xi = zeta_future - zeta_past
        IF (normalize_residual) THEN
          xi = 2._8 * xi / ( ABS(zeta_future) + ABS(zeta_past) )
        ENDIF

        constraint_cost = constraint_cost + xi * xi

      ENDIF
    ENDDO

    END FUNCTION constraint_cost
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE get_prev_timestep(xb0,xb1)
    REAL(KIND=8), DIMENSION(:), INTENT( in )  :: xb0
    REAL(KIND=8), DIMENSION(:), INTENT( out ) :: xb1

#if defined MPI
    ! Send block to next timestep
    IF (iproc_next.LE.nproc) THEN
      CALL MPI_SEND(xb1,nk,MPI_DOUBLE_PRECISION, &
      &             iproc_next,mpi_comm_flowsampler,mpi_code)
    ENDIF
    ! Receive block from previous timestep
    IF (iproc_prev.GE.0) THEN
      CALL MPI_RECV(xb0,nk,MPI_DOUBLE_PRECISION, &
      &             iproc_prev,mpi_comm_flowsampler,mpi_status_msg,mpi_code)
    ENDIF
#endif

    END SUBROUTINE get_prev_timestep
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE rebuild_timestep(xb,field)
    REAL(KIND=8), DIMENSION(:), INTENT( in )    :: xb
    REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: field

#if defined MPI
    ! Gather field for distributed blocks
    CALL MPI_ALLGATHER(xb,   nk,MPI_DOUBLE_PRECISION, &
      &                field,nk,MPI_DOUBLE_PRECISION, &
      &                mpi_comm_timestep,mpi_code)
#endif

    END SUBROUTINE rebuild_timestep
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE flowsampler_adv_constrint
