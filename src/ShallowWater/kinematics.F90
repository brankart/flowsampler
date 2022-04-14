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
!                        MODULE KINEMATICS
!
!---------------------------------------------------------------------
! Kinematic operators on the sphere
! by Jean-Michel Brankart, June 2019
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! velocity : compute velocity from stream function
! vorticity : compute relative vorticity from velocity
! ----------------------------------------------------------------------
MODULE flowsampler_kinematics
  use flowsampler_grid
  use flowsampler_units
  IMPLICIT NONE
  PRIVATE

  PUBLIC velocity, vorticity, velocity_norm

#if defined MPI
  ! Public definitions for MPI
  include "mpif.h"
  INTEGER, PUBLIC, SAVE  :: mpi_comm_flow_sampler_kinematics=mpi_comm_world   ! definition of module global communicator
  INTEGER, save :: mpi_code
  INTEGER, save :: iproc=0, nproc=1
#endif

  CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE velocity(psi,u,v)
!----------------------------------------------------------------------
! ** Purpose :   Compute velocity (u,v) from stream function (psi)
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: psi
    REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: u, v

    INTEGER :: i, j, k, ip1, jp1
    REAL(KIND=8) :: coslat, invdx, invdy, dlonratio, dlatratio

    ! Check size of input vectors
    IF (SIZE(psi,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(psi,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(u,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(u,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(v,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(v,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_kinematics'

    u=0. ; v=0.

    invdx = invdlon ; invdy = invdlat

#if defined MPI

    CALL mpi_comm_size(mpi_comm_flow_sampler_kinematics,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flow_sampler_kinematics,iproc,mpi_code)

    DO k=iproc,nlon*nlat-1,nproc
      i = 1 + MOD(k,nlon)
      j = 1 + k/nlon
      ip1 = 1 + MOD(i,nlon)
      jp1 = 1 + MOD(j,nlat)
      coslat = coslatitude(j)
      IF (physical_units) THEN
        CALL distance_ratio(j,dlonratio,dlatratio)
        invdx = invdlon / dlonratio
        invdy = invdlat / dlatratio
      ENDIF
      u(i,j) = - ( psi(i,jp1) - psi(i,j) ) * invdy
      v(i,j) = ( psi(ip1,j) - psi(i,j) ) * invdx / coslat
    ENDDO

    call MPI_ALLREDUCE (MPI_IN_PLACE,u,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                  MPI_SUM,mpi_comm_flow_sampler_kinematics,mpi_code)
    call MPI_ALLREDUCE (MPI_IN_PLACE,v,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                  MPI_SUM,mpi_comm_flow_sampler_kinematics,mpi_code)

#else

    DO j=1,nlatœ
      jp1 = 1 + MOD(j,nlat)
      coslat = coslatitude(j)
      IF (physical_units) THEN
        CALL distance_ratio(j,dlonratio,dlatratio)
        invdx = invdlon / dlonratio
        invdy = invdlat / dlatratio
      ENDIF
      DO i=1,nlon
        ip1 = 1 + MOD(i,nlon)
        u(i,j) = - ( psi(i,jp1) - psi(i,j) ) * invdy
        v(i,j) = ( psi(ip1,j) - psi(i,j) ) * invdx / coslat
      ENDDO
    ENDDO

#endif

    u(:,nlat) = 0.
    IF (.NOT.periodic) v(nlon,:) = 0.
    IF (npole) v(:,nlat) = 0.
    IF (spole) v(:,1) = 0.

    IF (physical_units) CALL uv2geo(u,v)

    END SUBROUTINE velocity
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE velocity_norm(u,v,unorm)
!----------------------------------------------------------------------
! ** Purpose :   Compute vorticity(zeta) from velocity (u,v)
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: u, v
    REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: unorm

    INTEGER :: i, j, k, im1
    REAL(KIND=8) :: umean, vmean

    unorm = 0.

    ! Check size of input vectors
    IF (SIZE(unorm,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(unorm,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(u,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(u,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(v,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(v,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_kinematics'

#if defined MPI

    CALL mpi_comm_size(mpi_comm_flow_sampler_kinematics,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flow_sampler_kinematics,iproc,mpi_code)


    DO k=iproc,nlon*(nlat-2)-1,nproc
      i = 1 + MOD(k,nlon)
      j = 2 + k/nlon
      im1 = 1 + MOD(nlon+i-2,nlon)
      umean = ( u(i,j-1) + u(i,j) ) * 0.5_8
      vmean = ( v(im1,j) + v(i,j) ) * 0.5_8
      unorm(i,j) = SQRT(umean * umean + vmean * vmean)
    ENDDO

    call MPI_ALLREDUCE (MPI_IN_PLACE,unorm,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                  MPI_SUM,mpi_comm_flow_sampler_kinematics,mpi_code)

#else

    DO j=2,nlat-1
      DO i=1,nlon
        im1 = 1 + MOD(nlon+i-2,nlon)
        umean = ( u(i,j-1) + u(i,j) ) * 0.5_8
        vmean = ( v(im1,j) + v(i,j) ) * 0.5_8
        unorm(i,j) = SQRT(umean * umean + vmean * vmean)
      ENDDO
    ENDDO

#endif

    IF (.NOT.periodic) unorm(1,:) = 0.
    IF (.NOT.periodic) unorm(nlon,:) = 0.

    END SUBROUTINE velocity_norm
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE vorticity(u,v,zeta)
!----------------------------------------------------------------------
! ** Purpose :   Compute vorticity(zeta) from velocity (u,v)
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: u, v
    REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: zeta

    INTEGER :: i, j, k, im1
    REAL(KIND=8) :: coslat, sinlat, zetan, zetas
    REAL(KIND=8) :: invdx, invdy, dlonratio, dlatratio

    ! Check size of input vectors
    IF (SIZE(zeta,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(zeta,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(u,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(u,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(v,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_kinematics'
    IF (SIZE(v,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_kinematics'

    zeta=0.

    invdx = invdlon ; invdy = invdlat ; dlonratio =1._8

#if defined MPI

    CALL mpi_comm_size(mpi_comm_flow_sampler_kinematics,nproc,mpi_code)
    CALL mpi_comm_rank(mpi_comm_flow_sampler_kinematics,iproc,mpi_code)

    DO k=iproc,nlon*(nlat-2)-1,nproc
      i = 1 + MOD(k,nlon)
      j = 2 + k/nlon
      im1 = 1 + MOD(nlon+i-2,nlon)
      coslat = coslatitude(j)
      sinlat = sinlatitude(j)
      IF (physical_units) THEN
        CALL distance_ratio(j,dlonratio,dlatratio)
        invdx = invdlon / dlonratio
        invdy = invdlat / dlatratio
      ENDIF
      zeta(i,j) = 0.5_8 * ( u(i,j) + u(i,j-1) ) * sinlat / ( coslat * dlonratio )
      zeta(i,j) = zeta(i,j) - ( u(i,j) - u(i,j-1) ) * invdy
      zeta(i,j) = zeta(i,j) + ( v(i,j) - v(im1,j) ) * invdx / coslat
    ENDDO

    call MPI_ALLREDUCE (MPI_IN_PLACE,zeta,nlon*nlat,MPI_DOUBLE_PRECISION, &
     &                  MPI_SUM,mpi_comm_flow_sampler_kinematics,mpi_code)

#else

    DO j=2,nlat-1
      coslat = coslatitude(j)
      sinlat = sinlatitude(j)
      IF (physical_units) THEN
        CALL distance_ratio(j,dlonratio,dlatratio)
        invdx = invdlon / dlonratio
        invdy = invdlat / dlatratio
      ENDIF
      DO i=1,nlon
        im1 = 1 + MOD(i-2,nlon)
        zeta(i,j) = 0.5_8 * ( u(i,j) + u(i,j-1) ) * sinlat / ( coslat * dlonratio )
        zeta(i,j) = zeta(i,j) - ( u(i,j) - u(i,j-1) ) * invdy
        zeta(i,j) = zeta(i,j) + ( v(i,j) - v(im1,j) ) * invdx / coslat
      ENDDO
    ENDDO

#endif

    IF (.NOT.periodic) zeta(1,:) = 0.
    IF (.NOT.periodic) zeta(nlon,:) = 0.

    ! Compute vorticity at the poles using Stokes' theorem
    IF (periodic) THEN
      IF (latmin==-90.) THEN
        zetas = 0.
        DO i=1,nlon
          zetas = zetas - u(i,1)
        ENDDO
        zeta(:,1) = 4._8 * zetas * invdy / nlon
      ENDIF

      IF (latmax==90.) THEN
        zetan = 0.
        DO i=1,nlon
          zetan = zetan + u(i,nlat-1)
        ENDDO
        zeta(:,nlat) = 4._8 * zetan * invdy / nlon
      ENDIF
    ENDIF

    END SUBROUTINE vorticity
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE flowsampler_kinematics
