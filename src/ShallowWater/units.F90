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
!                        MODULE UNITS
!
!---------------------------------------------------------------------
! Kinematic operators on the sphere
! by Jean-Michel Brankart, April 2022
! ----------------------------------------------------------------------
! List of routines/functions :
! ----------------------------
! uv2geo : transform to geostrophic velocities
! ----------------------------------------------------------------------
MODULE flowsampler_units
  use flowsampler_grid
  IMPLICIT NONE
  PRIVATE

  PUBLIC distance_ratio, uv2geo, gravityoverf_ratio

#if defined MPI
  ! Public definitions for MPI
  include "mpif.h"
  INTEGER, PUBLIC, SAVE  :: mpi_comm_flow_sampler_units=mpi_comm_world   ! definition of module global communicator
  INTEGER, save :: mpi_code
  INTEGER, save :: iproc=0, nproc=1
#endif

  LOGICAL, PUBLIC, SAVE :: physical_units=.FALSE. ! use of physical units
  LOGICAL, PUBLIC, SAVE :: normalize_residual = .FALSE.

  ! Use QG/PV model ?
  LOGICAL, PUBLIC, SAVE :: qg_model=.FALSE.
  LOGICAL, PUBLIC, SAVE :: pv_model=.FALSE.
  REAL(KIND=8), PUBLIC, SAVE :: rossby_radius=3.d4 ! Rossby radius in meter
  REAL(KIND=8), PUBLIC, SAVE :: dissip_rate=0. ! Dissipation rate 
  REAL(KIND=8), PUBLIC, SAVE :: adt_ref=-1. ! Reference topography level

  ! Earth geometry and angular velocity
  REAL(KIND=8), PUBLIC, PARAMETER :: earthradius=6.371d6                     ! meters
  REAL(KIND=8), PUBLIC, PARAMETER :: earthsideralday=86164.1                 ! seconds
  REAL(KIND=8), PUBLIC, PARAMETER :: earthvorticity=2._8*pi/earthsideralday  ! rad/second
  ! Ellipsoid correction
  LOGICAL, PUBLIC :: ellipsoid_correction=.FALSE.
  REAL(KIND=8), PUBLIC, PARAMETER :: earthradiusmin=6356752.3                ! meters
  REAL(KIND=8), PUBLIC, PARAMETER :: earthradiusmax=6378137.0                ! meters
  REAL(KIND=8), PUBLIC, PARAMETER :: eccentricity=SQRT(1._8-(earthradiusmin/earthradiusmax)**2)
  REAL(KIND=8), PUBLIC, PARAMETER :: earthradiusratio=earthradiusmax/earthradius
  ! Earth surface gravity: International gravity formula 1980
  REAL(KIND=8), PUBLIC, PARAMETER :: earthgravity=9.81                       ! m/s2
  REAL(KIND=8), PUBLIC, PARAMETER :: earthgravref=9.780327
  REAL(KIND=8), PUBLIC, PARAMETER :: earthgravcoef1=5.2792d-3
  REAL(KIND=8), PUBLIC, PARAMETER :: earthgravcoef2=-2.32d-5

  CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE uv2geo(u,v)
!----------------------------------------------------------------------
! ** Purpose :   transform into geostrophic velocities
!                assuming that the input components are computed
!                using dynamic topography as stream function
!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT( inout ) :: u, v

    INTEGER :: i, j, k
    REAL(KIND=8) :: gravityoverf

    ! Check size of input vectors
    IF (SIZE(u,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_units'
    IF (SIZE(u,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_units'
    IF (SIZE(v,1).NE.nlon) STOP 'Inconsistent size in flow_sampler_units'
    IF (SIZE(v,2).NE.nlat) STOP 'Inconsistent size in flow_sampler_units'

#if defined MPI

    DO k=iproc,nlon*nlat-1,nproc
      i = 1 + MOD(k,nlon)
      j = 1 + k/nlon
      CALL gravityoverf_ratio(j,gravityoverf)
      u(i,j) = u(i,j) * gravityoverf
      v(i,j) = v(i,j) * gravityoverf
    ENDDO

#else

    DO j=1,nlat
      CALL gravityoverf_ratio(j,gravityoverf)
      DO i=1,nlon
        u(i,j) = u(i,j) * gravityoverf
        v(i,j) = v(i,j) * gravityoverf
      ENDDO
    ENDDO

#endif

    END SUBROUTINE uv2geo
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE distance_ratio(j,dlonratio,dlatratio)

    INTEGER, INTENT( in ) :: j
    REAL(KIND=8), INTENT( out ) :: dlonratio,dlatratio

    REAL(KIND=8) :: sinlat, coslat, e2, e2sin2

    IF (ellipsoid_correction) THEN
      sinlat = sinlatitude(j)
      coslat = coslatitude(j)
      e2 = eccentricity * eccentricity
      e2sin2 = e2 * sinlat * sinlat
      dlonratio = earthradiusmax / SQRT(1._8-e2sin2)
      dlatratio = earthradiusmax * (1._8-e2) / (1._8-e2sin2)**1.5_8
    ELSE
      dlonratio = earthradius
      dlatratio = earthradius
    ENDIF

    END SUBROUTINE distance_ratio
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE gravityoverf_ratio(j,gravityoverf)

    INTEGER, INTENT( in ) :: j
    REAL(KIND=8), INTENT( out ) :: gravityoverf

    REAL(KIND=8) :: sinlat, sin2lat, sin4lat

    IF (ellipsoid_correction) THEN
      sinlat = sinlatitude(j)
      sin2lat = sinlat * sinlat
      sin4lat = sin2lat * sin2lat
      gravityoverf = 1._8 + earthgravcoef1 * sin2lat + earthgravcoef2 * sin4lat
      gravityoverf = gravityoverf * earthgravref
      gravityoverf = gravityoverf / ( 2._8 * earthvorticity * sinlat )
    ELSE
      sinlat = sinlatitude(j)
      gravityoverf = earthgravity / ( 2._8 * earthvorticity * sinlat )
    ENDIF

    END SUBROUTINE gravityoverf_ratio
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE flowsampler_units
