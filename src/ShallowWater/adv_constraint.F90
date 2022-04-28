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
  INTEGER, SAVE :: nk, k0, k1 ! size and bounds of block of data
  REAL(KIND=8), DIMENSION(:,:), PUBLIC, ALLOCATABLE :: u0, v0, zeta0
  REAL(KIND=8), DIMENSION(:,:), PUBLIC, ALLOCATABLE :: u1, v1, zeta1
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: xblock0, xblock1

  ! Private arrays with metric/scaling factors
  REAL(KIND=8), DIMENSION(:), SAVE, ALLOCATABLE :: invdxfac, invdyfac
  REAL(KIND=8), DIMENSION(:), SAVE, ALLOCATABLE :: zetafacu, zetafacv
  REAL(KIND=8), DIMENSION(:), SAVE, ALLOCATABLE :: dtfacx, dtfacy
  REAL(KIND=8), SAVE :: zeta_planetary

#if defined MPI
  ! Definitions for MPI
  include "mpif.h"
  INTEGER, PUBLIC, SAVE  :: mpi_comm_flowsampler=mpi_comm_world  ! definition of module global communicator
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpi_status_msg
  INTEGER, save :: mpi_code
#endif

  ! Definitions for MPI
  INTEGER, save :: iproc=0, nproc=1
  INTEGER, save :: nproct, iproct, iproc_prev, iproc_next
  INTEGER, PUBLIC, save :: mpi_comm_timestep
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ncommt  ! communicator for every timestep

  CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE constraint_init(nt,dt)
!----------------------------------------------------------------------
! ** Purpose : constraint initialization
!----------------------------------------------------------------------
    INTEGER, INTENT( in ) :: nt
    REAL(KIND=8), INTENT( in ) :: dt

    INTEGER :: allocok, iiproct, jt, j
    INTEGER :: igrp_global
    INTEGER, DIMENSION(:), ALLOCATABLE ::   irankt ! list of processors for current timestep
    INTEGER, DIMENSION(:), ALLOCATABLE ::   igrpt  ! group ID for every timestep
    REAL(KIND=8) :: coslat, coslatu, dlonratio, dlatratio, gravityoverf
    REAL(KIND=8) :: invdx, invdy

#if defined MPI
    CALL mpi_comm_size(mpi_comm_flowsampler,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flowsampler,iproc,mpi_code)
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
    CALL mpi_comm_group( mpi_comm_flowsampler, igrp_global, mpi_code )

    ! Create one communicator for each timestep
    DO jt = 1, nt
      irankt(:) = (/ (iiproct, iiproct=0,nproct-1) /)
      irankt(:) = irankt(:) + ( jt - 1 ) * nproct
      ! Create the group for current timestep from the global group
      CALL mpi_group_incl( igrp_global, nproct, irankt, igrpt(jt), mpi_code )
      ! Create communicator for jt
      CALL mpi_comm_create( mpi_comm_flowsampler, igrpt(jt), ncommt(jt), mpi_code )
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
    allocate(invdxfac(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(invdyfac(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(zetafacu(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(zetafacv(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(dtfacx(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'
    allocate(dtfacy(nlat),stat=allocok)
    IF (allocok.NE.0) STOP 'Allocation error in flowsampler: constraint_init'

    ! compute metric/rescaling arrays
    DO j=1,nlat
      coslat = coslatitude(j)
      coslatu = coslatitude_u(j)
      CALL distance_ratio(j,dlonratio,dlatratio)
      CALL gravityoverf_ratio(j,gravityoverf)
      ! Inverse grid size from gradian^(-1) to meters^(-1)
      invdx = invdlon / ( dlonratio * coslat )
      invdy = invdlat / dlatratio
      ! Factors to obtain velocity from dynamic toprography differences
      invdxfac(j) = gravityoverf * invdx
      invdyfac(j) = gravityoverf * invdy
      ! Factors to obtain relative vorticity from velocity differences
          !invdx(j) = invdlon / ( dlonratio * coslat )
          !invdy(j) = invdlat / dlatratio
          !zetafac(j) = 0.5_8 * sinlat / ( dlonratio * coslat )
      zetafacu(j) = coslatu * invdy / ( invdx * coslat )
      zetafacv(j) = invdx
      ! Factor to obtain displacement from velocity in one time step
      IF (spherical_delta) THEN
        ! WARNING: In this case: dx, dy in radian along a great circle
        dtfacx(j) = dt * 0.25_8 / dlonratio
        dtfacy(j) = dt * 0.25_8 / dlatratio
      ELSE
        ! WARNING: In this case: dx, dy in degrees along coordinates isolines
        dtfacx(j) = dt * 0.25_8 * rad2deg / ( dlonratio * coslat )
        dtfacy(j) = dt * 0.25_8 * rad2deg / dlatratio
      ENDIF
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

    INTEGER :: i,j,k, kk
    REAL(KIND=8) :: zetau1, zetau2, zetav, dx, dy, lon, lat
    REAL(KIND=8) :: zeta_past, zeta_future, xi
    LOGICAL :: inmask, ingrid

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
        xblock1(kk) = ( adt(i+1,j) - adt(i,j) ) * invdxfac(j)
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
            !zetau1 = ( u1(i,j) + u1(i,j-1) ) * zetafac(j)
            !zetau2 = - ( u1(i,j) - u1(i,j-1) ) * invdy(j)
            !zetav = ( v1(i,j) - v1(i-1,j) ) * invdx(j)
            !xblock1(kk)= zeta1 + zeta2 + zeta3
        zetau1 = - u1(i,j) * zetafacu(j)
        zetau2 = u1(i,j-1) * zetafacu(j-1)
        zetav = v1(i,j) - v1(i-1,j)
        xblock1(kk)= ( zetau1 + zetau2 + zetav ) * zetafacv(j)
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
      inmask = inmask .AND. (iproc_prev.GE.0) ! exclude first timestep
      IF (inmask) THEN
        ! advect potential vorticity from past situation
        dx = - ( u0(i,j) + u0(i,j-1) ) * dtfacx(j)
        dy = - ( v0(i,j) + v0(i-1,j) ) * dtfacy(j)
        call get_location(i,j,dx,dy,lon,lat)
        call grid_interp(zeta0,lon,lat,zeta_past)
        ingrid = zeta_past.NE.spval
        IF (ingrid) zeta_past = zeta_past + zeta_planetary * SIN(lat*deg2rad)
        ! advect potential vorticity from past situation
        dx = ( u1(i,j) + u1(i,j-1) ) * dtfacx(j)
        dy = ( v1(i,j) + v1(i-1,j) ) * dtfacy(j)
        call get_location(i,j,dx,dy,lon,lat)
        call grid_interp(zeta1,lon,lat,zeta_future)
        ingrid = ingrid .AND. (zeta_future.NE.spval)
        IF (ingrid) zeta_future = zeta_future + zeta_planetary * SIN(lat*deg2rad)

        IF (ingrid) THEN
          xi = zeta_future - zeta_past
          IF (normalize_residual) THEN
            xi = 2._8 * xi / ( ABS(zeta_future) + ABS(zeta_past) )
          ENDIF

          constraint_cost = constraint_cost + xi * xi

        ENDIF
      ENDIF

      ! For debug
      !kk = k - k0 + 1
      !xblock1(kk)= zeta0(i,j)
      !xblock1(kk)= zeta_past
      !IF (inmask.AND.ingrid) THEN
      !  xblock1(kk)= xi
      !ELSE
      !  xblock1(kk)= 0.
      !ENDIF

    ENDDO

    ! For debug
    !CALL rebuild_timestep(xblock1,zeta1)

    END FUNCTION constraint_cost
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE get_prev_timestep(xb0,xb1)
    REAL(KIND=8), DIMENSION(:), INTENT( in )  :: xb0
    REAL(KIND=8), DIMENSION(:), INTENT( out ) :: xb1

    !INTEGER :: ip, np
    INTEGER, parameter :: tag=110

#if defined MPI
    !print *, 'in send:,',iproc
    !CALL mpi_comm_size(mpi_comm_flowsampler,np,mpi_code)
    !CALL mpi_comm_rank(mpi_comm_flowsampler,ip,mpi_code)
    !print *, 'before send:,',np,ip,iproc_prev,iproc_next

    ! Send block to next timestep
    IF (iproc_next.LT.nproc) THEN
      CALL MPI_SEND(xb1,nk,MPI_DOUBLE_PRECISION, &
      &             iproc_next,tag,mpi_comm_flowsampler,mpi_code)
    ENDIF
    ! Receive block from previous timestep
    IF (iproc_prev.GE.0) THEN
      CALL MPI_RECV(xb0,nk,MPI_DOUBLE_PRECISION, &
      &             iproc_prev,tag,mpi_comm_flowsampler,mpi_status_msg,mpi_code)
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
END MODULE flowsampler_adv_constraint
