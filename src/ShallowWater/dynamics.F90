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
!                        MODULE DYNAMICS
!
!---------------------------------------------------------------------
! Dynamic operators on the sphere
! by Jean-Michel Brankart, June 2019
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! advection : advection of potential vorticity
! diffusion : laplacian diffusion operator
! turbulent_dissipation : nonlinear turbulent dissipation operator
! ----------------------------------------------------------------------
MODULE flowsampler_dynamics
  use flowsampler_grid
  use flowsampler_units
  IMPLICIT NONE
  PRIVATE

  PUBLIC advection, diffusion, turbulent_dissipation, advection_rate

#if defined MPI
  ! Public definitions for MPI
  include "mpif.h"
  INTEGER, PUBLIC, SAVE  :: mpi_comm_flow_sampler_dynamics=mpi_comm_world   ! definition of module global communicator
  INTEGER, save :: mpi_code
  INTEGER, save :: iproc=0, nproc=1
#endif

  CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE advection(advresidual,zeta0,u0,v0,zeta1,u1,v1,rossby,dt)
!----------------------------------------------------------------------
! ** Purpose :   Compute residual with respect to advection of potential vorticity
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: zeta0, u0, v0
    REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: zeta1, u1, v1
    REAL(KIND=8), DIMENSION(:,:), INTENT( out) :: advresidual
    REAL(KIND=8), INTENT( in ) :: rossby, dt

    INTEGER :: i, j, k
    REAL(KIND=8) :: zeta_past, zeta_future, zeta_planetary
    REAL(KIND=8) :: dx, dy, lon, lat, invdx, invdy, dlonratio, dlatratio

    ! Check size of input vectors
    IF (SIZE(zeta0,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(zeta0,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(u0,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(u0,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(v0,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(v0,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(zeta1,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(zeta1,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(u1,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(u1,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(v1,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(v1,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(advresidual,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(advresidual,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'

    ! Initialization
    advresidual = 0.

    invdx = 1._8 ; invdy = 1._8

    IF (physical_units) THEN
      zeta_planetary = 2._8 * earthvorticity
    ELSE
      zeta_planetary = 1._8 / rossby
    ENDIF

#if defined MPI

    CALL mpi_comm_size(mpi_comm_flow_sampler_dynamics,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flow_sampler_dynamics,iproc,mpi_code)

    DO k=iproc,nlon*nlat-1,nproc
      i = 1 + MOD(k,nlon)
      j = 1 + k/nlon
      IF (physical_units) THEN
        CALL distance_ratio(j,dlonratio,dlatratio)
        invdx = 1._8 / dlonratio
        invdy = 1._8 / dlatratio
      ENDIF
      ! advect potential vorticity from past situation
      dx = - u0(i,j) * dt * 0.5_8 * invdx
      dy = - v0(i,j) * dt * 0.5_8 * invdy
      call get_location(i,j,dx,dy,lon,lat)
      call grid_interp(zeta0,lon,lat,zeta_past)
      IF (zeta_past.NE.spval) zeta_past = zeta_past + zeta_planetary * SIN(lat*deg2rad)
      ! advect potential vorticity from past situation
      dx = u1(i,j) * dt * 0.5_8 * invdx
      dy = v1(i,j) * dt * 0.5_8 * invdy
      call get_location(i,j,dx,dy,lon,lat)
      call grid_interp(zeta1,lon,lat,zeta_future)
      IF (zeta_future.NE.spval) zeta_future = zeta_future + zeta_planetary * SIN(lat*deg2rad)
      ! compute misfit between future and past potential vorticity
      IF ( (zeta_past.NE.spval) .AND. (zeta_future.NE.spval) ) THEN
        IF (normalize_residual) THEN
          advresidual(i,j) = 2._8 * ( zeta_future - zeta_past ) / &
                  & ( ABS(zeta_future) + ABS(zeta_past) )
        ELSE
          advresidual(i,j) = ( zeta_future - zeta_past ) / dt
        ENDIF
      ENDIF
    ENDDO

    call MPI_ALLREDUCE (MPI_IN_PLACE,advresidual,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                  MPI_SUM,mpi_comm_flow_sampler_dynamics,mpi_code)

#else

    DO j=1,nlat
      IF (physical_units) THEN
        CALL distance_ratio(j,dlonratio,dlatratio)
        invdx = 1._8 / dlonratio
        invdy = 1._8 / dlatratio
      ENDIF
      DO i=1,nlon
        ! advect potential vorticity from past situation
        dx = - u0(i,j) * dt * 0.5_8 * invdx
        dy = - v0(i,j) * dt * 0.5_8 * invdy
        call get_location(i,j,dx,dy,lon,lat)
        call grid_interp(zeta0,lon,lat,zeta_past)
        IF (zeta_past.NE.spval) zeta_past = zeta_past + zeta_planetary * SIN(lat*deg2rad)
        ! advect potential vorticity from past situation
        dx = u1(i,j) * dt * 0.5_8 * invdx
        dy = v1(i,j) * dt * 0.5_8 * invdy
        call get_location(i,j,dx,dy,lon,lat)
        call grid_interp(zeta1,lon,lat,zeta_future)
        IF (zeta_future.NE.spval) zeta_future = zeta_future + zeta_planetary * SIN(lat*deg2rad)
        ! compute misfit between future and past potential vorticity
        IF ( (zeta_past.NE.spval) .AND. (zeta_future.NE.spval) ) THEN
          IF (normalize_residual) THEN
            advresidual(i,j) = 2._8 * ( zeta_future - zeta_past ) / &
                  & ( ABS(zeta_future) + ABS(zeta_past) )
          ELSE
            advresidual(i,j) = ( zeta_future - zeta_past ) / dt
          ENDIF
        ENDIF
      ENDDO
    ENDDO

#endif

    END SUBROUTINE advection
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE advection_rate(advrate,zeta,u,v,rossby)
!----------------------------------------------------------------------
! ** Purpose :   Compute advection rate of potential vorticity
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: zeta, u, v
    REAL(KIND=8), DIMENSION(:,:), INTENT( out) :: advrate
    REAL(KIND=8), INTENT( in ) :: rossby

    INTEGER :: i, j, k, ip1, im1
    REAL(KIND=8) :: coslat, sinlatp1, sinlatm1
    REAL(KIND=8) :: zetai, zetaj, umean, vmean, zeta_planetary

    advrate = 0.

    ! Check size of input vectors
    IF (SIZE(zeta,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(zeta,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(u,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(u,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(v,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(v,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(advrate,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(advrate,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'

    zeta_planetary = 1._8 / rossby

#if defined MPI

    CALL mpi_comm_size(mpi_comm_flow_sampler_dynamics,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flow_sampler_dynamics,iproc,mpi_code)


    DO k=iproc,nlon*(nlat-2)-1,nproc
      i = 1 + MOD(k,nlon)
      j = 2 + k/nlon
      ip1 = 1 + MOD(i,nlon)
      im1 = 1 + MOD(nlon+i-2,nlon)
      coslat = coslatitude(j)
      zetai = ( zeta(ip1,j) - zeta(im1,j) ) / coslat * invdlon
      zetaj = ( zeta(i,j+1) - zeta(i,j-1) ) * invdlat
      sinlatm1 = sinlatitude(j-1)
      sinlatp1 = sinlatitude(j+1)
      zetaj = zetaj + ( sinlatp1 - sinlatm1 ) * invdlat * zeta_planetary
      umean = ( u(i,j+1) + u(i,j) ) * 0.5_8
      vmean = ( v(im1,j) + v(i,j) ) * 0.5_8
      advrate(i,j) = zetai * umean + zetaj * vmean
      !advrate(i,j) = advrate(i,j) / sqrt( (zetai**2+zetaj**2) * (umean**2+vmean**2) )
    ENDDO

    call MPI_ALLREDUCE (MPI_IN_PLACE,advrate,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                  MPI_SUM,mpi_comm_flow_sampler_dynamics,mpi_code)

#else

    DO j=2,nlat-1
      coslat = coslatitude(j)
      DO i=1,nlon
        ip1 = 1 + MOD(i,nlon)
        im1 = 1 + MOD(nlon+i-2,nlon)
        zetai = ( zeta(ip1,j) - zeta(im1,j) ) / coslat * invdlon
        zetaj = ( zeta(i,j+1) - zeta(i,j-1) ) * invdlat
        sinlatm1 = sinlatitude(j-1)
        sinlatp1 = sinlatitude(j+1)
        zetaj = zetaj + ( sinlatp1 - sinlatm1 ) * invdlat * zeta_planetary
        umean = ( u(i,j+1) + u(i,j) ) * 0.5_8
        vmean = ( v(im1,j) + v(i,j) ) * 0.5_8
        advrate(i,j) = zetai * umean + zetaj * vmean
      ENDDO
    ENDDO

#endif

    IF (.NOT.periodic) advrate(1,:) = 0.
    IF (.NOT.periodic) advrate(nlon,:) = 0.

    END SUBROUTINE advection_rate
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE diffusion(roughness,zeta)
!----------------------------------------------------------------------
! ** Purpose :   Compute diffusion of relative vorticity
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: zeta
    REAL(KIND=8), DIMENSION(:,:), INTENT( out) :: roughness

    INTEGER :: i, j, k, ip1, im1
    REAL(KIND=8) :: coslat, sinlat
    REAL(KIND=8) :: zetaj, zetaii, zetajj, laplacian

    roughness = 0.

    ! Check size of input vectors
    IF (SIZE(zeta,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(zeta,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(roughness,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(roughness,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'

#if defined MPI

    CALL mpi_comm_size(mpi_comm_flow_sampler_dynamics,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flow_sampler_dynamics,iproc,mpi_code)


    DO k=iproc,nlon*(nlat-2)-1,nproc
      i = 1 + MOD(k,nlon)
      j = 2 + k/nlon
      ip1 = 1 + MOD(i,nlon)
      im1 = 1 + MOD(nlon+i-2,nlon)

      sinlat = sinlatitude(j)
      coslat = coslatitude(j)

      zetajj = ( zeta(i,j+1) + zeta(i,j-1) - 2 * zeta(i,j) ) * invdlat * invdlat

      zetaj = ( zeta(i,j+1) - zeta(i,j-1) ) * invdlat * 0.5_8
      zetaii = ( zeta(ip1,j) + zeta(im1,j) - 2 * zeta(i,j) ) * invdlon * invdlon
      zetaii = zetaii / (coslat * coslat) - zetaj * sinlat / coslat

      laplacian = zetaii + zetajj
      roughness(i,j) = laplacian * laplacian
    ENDDO

    call MPI_ALLREDUCE (MPI_IN_PLACE,roughness,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                  MPI_SUM,mpi_comm_flow_sampler_dynamics,mpi_code)

#else

    DO j=2,nlat-1
      sinlat = sinlatitude(j)
      coslat = coslatitude(j)
      DO i=1,nlon
        ip1 = 1 + MOD(i,nlon)
        im1 = 1 + MOD(nlon+i-2,nlon)

        zetajj = ( zeta(i,j+1) + zeta(i,j-1) - 2 * zeta(i,j) ) * invdlat * invdlat

        zetaj = ( zeta(i,j+1) - zeta(i,j-1) ) * invdlat * 0.5_8
        zetaii = ( zeta(ip1,j) + zeta(im1,j) - 2 * zeta(i,j) ) * invdlon * invdlon
        zetaii = zetaii / (coslat * coslat) - zetaj * sinlat / coslat

        laplacian = zetaii + zetajj
        roughness(i,j) = laplacian * laplacian
      ENDDO
    ENDDO

#endif

    roughness(:,2) = 0.
    roughness(:,nlat-1) = 0.
    IF (.NOT.periodic) roughness(1:2,:) = 0.
    IF (.NOT.periodic) roughness(nlon-1:nlon,:) = 0.

    IF (mask_poles) THEN
      roughness(:,1:mask_spole+1) = 0.
      roughness(:,nlat-mask_npole:nlat) = 0.
    ENDIF

    END SUBROUTINE diffusion
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE turbulent_dissipation(slope_penalty,zeta)
!----------------------------------------------------------------------
! ** Purpose :   Compute turbulent dissipation of relative vorticity
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: zeta
    REAL(KIND=8), DIMENSION(:,:), INTENT( out) :: slope_penalty

    INTEGER :: i, j, k, ip1, im1
    REAL(KIND=8) :: coslat
    REAL(KIND=8) :: zetaip, zetaim, zetajp, zetajm, slope_norm

    slope_penalty = 0.

    ! Check size of input vectors
    IF (SIZE(zeta,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(zeta,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(slope_penalty,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(slope_penalty,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'

#if defined MPI

    CALL mpi_comm_size(mpi_comm_flow_sampler_dynamics,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flow_sampler_dynamics,iproc,mpi_code)

    DO k=iproc,nlon*(nlat-2)-1,nproc
      i = 1 + MOD(k,nlon)
      j = 2 + k/nlon
      ip1 = 1 + MOD(i,nlon)
      im1 = 1 + MOD(nlon+i-2,nlon)

      coslat = coslatitude(j)

      zetaip = ( zeta(ip1,j) - zeta(i,j) ) * invdlon / coslat
      zetaim = ( zeta(i,j) - zeta(im1,j) ) * invdlon / coslat
      zetajp = ( zeta(i,j+1) - zeta(i,j) ) * invdlat
      zetajm = ( zeta(i,j) - zeta(i,j-1) ) * invdlat

      slope_norm = zetaip*zetaip + zetaim*zetaim + zetajp*zetajp + zetajm*zetajm
      slope_penalty(i,j) = 0.25_8 * slope_norm * slope_norm
    ENDDO

    call MPI_ALLREDUCE (MPI_IN_PLACE,slope_penalty,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                  MPI_SUM,mpi_comm_flow_sampler_dynamics,mpi_code)

#else

    DO j=2,nlat-1
      coslat = coslatitude(j)
      DO i=1,nlon
        ip1 = 1 + MOD(i,nlon)
        im1 = 1 + MOD(nlon+i-2,nlon)

        zetaip = ( zeta(ip1,j) - zeta(i,j) ) * invdlon / coslat
        zetaim = ( zeta(i,j) - zeta(im1,j) ) * invdlon / coslat
        zetajp = ( zeta(i,j+1) - zeta(i,j) ) * invdlat
        zetajm = ( zeta(i,j) - zeta(i,j-1) ) * invdlat

        slope_norm = zetaip*zetaip + zetaim*zetaim + zetajp*zetajp + zetajm*zetajm
        slope_penalty(i,j) = 0.25_8 * slope_norm * slope_norm
      ENDDO
    ENDDO

#endif

    slope_penalty(:,2) = 0.
    slope_penalty(:,nlat-1) = 0.
    IF (.NOT.periodic) slope_penalty(1:2,:) = 0.
    IF (.NOT.periodic) slope_penalty(nlon-1:nlon,:) = 0.

    IF (mask_poles) THEN
      slope_penalty(:,1:mask_spole+1) = 0.
      slope_penalty(:,nlat-mask_npole:nlat) = 0.
    ENDIF

    END SUBROUTINE turbulent_dissipation
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE diffusion_tensor(roughness,zeta)
!----------------------------------------------------------------------
! ** Purpose :   Compute diffusion of relative vorticity
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: zeta
    REAL(KIND=8), DIMENSION(:,:), INTENT( out) :: roughness

    INTEGER :: i, j, k, ip1, im1
    REAL(KIND=8) :: coslat, sinlat
    REAL(KIND=8) :: zetaim1, zetaip1, zetai, zetaj, zetaii, zetajj, zetaij

    roughness = 0.

    ! Check size of input vectors
    IF (SIZE(zeta,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(zeta,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(roughness,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_dynamics'
    IF (SIZE(roughness,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_dynamics'

#if defined MPI

    CALL mpi_comm_size(mpi_comm_flow_sampler_dynamics,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flow_sampler_dynamics,iproc,mpi_code)


    DO k=iproc,nlon*(nlat-2)-1,nproc
      i = 1 + MOD(k,nlon)
      j = 2 + k/nlon
      ip1 = 1 + MOD(i,nlon)
      im1 = 1 + MOD(nlon+i-2,nlon)

      sinlat = sinlatitude(j)
      coslat = coslatitude(j)

      zetajj = ( zeta(i,j+1) + zeta(i,j-1) - 2 * zeta(i,j) ) * invdlat * invdlat

      zetaip1 = ( zeta(ip1,j+1) - zeta(ip1,j-1) ) * invdlat * 0.5_8 
      zetaim1 = ( zeta(im1,j+1) - zeta(im1,j-1) ) * invdlat * 0.5_8 
      zetaij = ( zetaip1 - zetaim1 ) * invdlon * 0.5_8
      zetai = ( zeta(ip1,j) - zeta(im1,j) ) * invdlon * 0.5_8
      zetaij = - zetaij / coslat - zetai * sinlat / (coslat*coslat)

      zetaj = ( zeta(i,j+1) - zeta(i,j-1) ) * invdlat * 0.5_8
      zetaii = ( zeta(ip1,j) + zeta(im1,j) - 2 * zeta(i,j) ) * invdlon * invdlon
      zetaii = zetaii / (coslat * coslat) - zetaj * sinlat / coslat

      roughness(i,j) = zetajj * zetajj + 2 * zetaij * zetaij + zetaii * zetaii
    ENDDO

    call MPI_ALLREDUCE (MPI_IN_PLACE,roughness,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                  MPI_SUM,mpi_comm_flow_sampler_dynamics,mpi_code)

#else

    DO j=2,nlat-1
      sinlat = sinlatitude(j)
      coslat = coslatitude(j)
      DO i=1,nlon
        ip1 = 1 + MOD(i,nlon)
        im1 = 1 + MOD(nlon+i-2,nlon)

        zetajj = ( zeta(i,j+1) + zeta(i,j-1) - 2 * zeta(i,j) ) * invdlat * invdlat

        zetaip1 = ( zeta(ip1,j+1) - zeta(ip1,j-1) ) * invdlat * 0.5_8
        zetaim1 = ( zeta(im1,j+1) - zeta(im1,j-1) ) * invdlat * 0.5_8
        zetaij = ( zetaip1 - zetaim1 ) * invdlon * 0.5_8
        zetai = ( zeta(ip1,j) - zeta(im1,j) ) * invdlon * 0.5_8
        zetaij = - zetaij / coslat - zetai * sinlat / (coslat*coslat)

        zetaj = ( zeta(i,j+1) - zeta(i,j-1) ) * invdlat * 0.5_8
        zetaii = ( zeta(ip1,j) + zeta(im1,j) - 2 * zeta(i,j) ) * invdlon * invdlon
        zetaii = zetaii / (coslat * coslat) - zetaj * sinlat / coslat

        roughness(i,j) = zetajj * zetajj + 2 * zetaij * zetaij + zetaii * zetaii
      ENDDO
    ENDDO

#endif

    IF (.NOT.periodic) roughness(1,:) = 0.
    IF (.NOT.periodic) roughness(nlon,:) = 0.

    IF (mask_poles) THEN
      roughness(:,1:mask_spole+1) = 0.
      roughness(:,nlat-mask_npole:nlat) = 0.
    ENDIF

    END SUBROUTINE diffusion_tensor
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE flowsampler_dynamics
