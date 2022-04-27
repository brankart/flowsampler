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
!                        MODULE GRID
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
MODULE flowsampler_grid
  IMPLICIT NONE
  PRIVATE

  PUBLIC defgrid, coslatitude, coslatitude_u, sinlatitude
  PUBLIC get_location, grid_interp

  INTEGER, PUBLIC, SAVE :: nlon, nlat
  REAL(KIND=8), PUBLIC, SAVE :: lonmin, latmin
  REAL(KIND=8), PUBLIC, SAVE :: lonmax, latmax
  LOGICAL, PUBLIC, SAVE :: periodic=.FALSE.

  LOGICAL, PUBLIC, SAVE :: spherical_delta=.TRUE.

  LOGICAL, PUBLIC, SAVE :: npole = .FALSE.
  LOGICAL, PUBLIC, SAVE :: spole = .FALSE.
  LOGICAL, PUBLIC, SAVE :: mask_poles = .FALSE.
  INTEGER, PUBLIC, SAVE :: mask_spole = 6
  INTEGER, PUBLIC, SAVE :: mask_npole = 6

  REAL(KIND=8), PUBLIC, SAVE :: dlon, dlat
  REAL(KIND=8), PUBLIC, SAVE :: invdlon, invdlat
  REAL(KIND=8), SAVE :: latminrad, dlatrad
  REAL(KIND=8), SAVE :: lonminrad, dlonrad

  REAL(KIND=8), PUBLIC, PARAMETER :: pi=3.1415926535897932384626
  REAL(KIND=8), PUBLIC, PARAMETER :: deg2rad=pi/180._8
  REAL(KIND=8), PUBLIC, PARAMETER :: rad2deg=1._8/deg2rad
  REAL(KIND=8), PUBLIC, PARAMETER :: spval = -9999.

  CONTAINS
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE defgrid()

    ! Grid steps
    IF (periodic) THEN
      dlon = (lonmax-lonmin) / REAL(nlon,8)
    ELSE
      dlon = (lonmax-lonmin) / REAL(nlon-1,8)
    ENDIF
    dlat = (latmax-latmin) / REAL(nlat-1,8)

    ! Inverse grid steps in radians
    invdlon = 1._8 /(dlon*deg2rad)
    invdlat = 1._8 /(dlat*deg2rad)

    ! Longitude and latitude origin and step in radian
    lonminrad = lonmin * deg2rad
    dlonrad = dlon * deg2rad
    latminrad = latmin * deg2rad
    dlatrad = dlat * deg2rad

    END SUBROUTINE defgrid
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    FUNCTION coslatitude(j)

    INTEGER, INTENT( in ) :: j
    REAL(KIND=8) :: coslatitude

    coslatitude = COS( latminrad + (j-1) * dlatrad )

    END FUNCTION coslatitude
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    FUNCTION coslatitude_u(j)

    INTEGER, INTENT( in ) :: j
    REAL(KIND=8) :: coslatitude_u

    coslatitude_u = COS( latminrad + (j-1) * dlatrad + 0.5_8 * dlatrad )

    END FUNCTION coslatitude_u
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    FUNCTION sinlatitude(j)

    INTEGER, INTENT( in ) :: j
    REAL(KIND=8) :: sinlatitude

    sinlatitude = SIN( latminrad + (j-1) * dlatrad )

    END FUNCTION sinlatitude
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE get_location(i,j,dx,dy,lon,lat)

    INTEGER, INTENT( in ) :: i, j
    REAL(KIND=8) , INTENT( in ) :: dx,dy
    REAL(KIND=8) , INTENT( out ) :: lon, lat

    REAL(KIND=8) :: lonref, latref, delta, azimuth, cosdlon

    IF (spherical_delta) THEN
      ! WARNING: In this case: dx, dy in radian along a great circle

      ! compute reference coordinates in radian
      lonref = lonminrad + (i-1) * dlonrad
      latref = latminrad + (j-1) * dlatrad

      ! compute distance and azimuth in radian
      delta = SQRT( dx*dx + dy*dy )
      if (dy.eq.0.) then
        if (dx.gt.0) then
          azimuth = pi / 2._8    ! eastward
        else
          azimuth = - pi / 2._8  ! westward
        endif
      else
        azimuth = ATAN(dx/dy)  ! dx=0 -> northward (azimuth=0)
      endif

      ! cope with negative dy with same azimuth but negative delta
      ! -> delta > 0 means northward, delta < 0 means southward
      if (dy.lt.0.) then
        delta = - delta
      endif

      ! compute new coordinates
      lat = asin( sin(latref)*cos(delta) + cos(latref)*sin(delta)*cos(azimuth) )
      if (cos(latref).eq.0.) then
        ! reference is at one pole
        lon = lonref
      elseif (cos(lat).eq.0.) then
        ! new location is at one pole
        lon = lonref
      else
        ! both are elsewhere
        cosdlon = (cos(delta) - sin(latref)*sin(lat)) / ( cos(latref)*cos(lat) )
        cosdlon = MAX(MIN(cosdlon,1._8),-1._8)
        lon = lonref + sign(1._8,dx) * acos( cosdlon )
      endif

      ! compute new coordinates in degrees
      lon = lon * rad2deg
      lat = lat * rad2deg

    ELSE
      ! WARNING: In this case: dx, dy in degrees along coordinates isolines

      lon = lonmin + (i-1) * dlon + dx
      lat = latmin + (j-1) * dlat + dy

    ENDIF

    END SUBROUTINE get_location
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    SUBROUTINE grid_interp(field,lon,lat,locvalue)

    REAL(KIND=8), DIMENSION(:,:), INTENT( in ) :: field
    REAL(KIND=8) , INTENT( inout ) :: lon, lat
    REAL(KIND=8) , INTENT( out ) :: locvalue

    LOGICAL :: outside
    INTEGER :: i, j, ip1
    REAL(KIND=8) :: rlon, rlat, val1, val2

    outside = lat .LT. latmin+dlat
    outside = outside .OR. (lat .GT. latmax-dlat)
    IF (periodic) THEN
      DO WHILE (lon.LT.lonmin)
        lon = lon + 360.
      ENDDO
      DO WHILE (lon.GT.lonmax)
        lon = lon - 360.
      ENDDO
    ELSE
      outside = outside .OR. (lon .LT. lonmin+dlon)
      outside = outside .OR. (lon .GT. lonmax-dlon)
    ENDIF

    IF (outside) THEN
      locvalue = spval
    ELSE
      rlon = 1._8 + (lon-lonmin)/dlon
      rlat = 1._8 + (lat-latmin)/dlat
      IF (periodic) THEN
        i = MIN(INT(rlon),nlon-1)
      ELSE
        i = MIN(INT(rlon),nlon)
      ENDIF
      j = MAX(1,MIN(INT(rlat),nlat-1))
      rlon = rlon - REAL(i,8)
      rlat = rlat - REAL(j,8)
      ip1 = 1 + MOD(i,nlon)
      val1 = (1._8 - rlon) * field(i,j) + rlon * field(ip1,j)
      val2 = (1._8 - rlon) * field(i,j+1) + rlon * field(ip1,j+1)
      locvalue = (1._8 - rlat) * val1 + rlat * val2
    ENDIF

    END SUBROUTINE grid_interp
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
END MODULE flowsampler_grid
