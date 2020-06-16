program bl_profile
  use GeneralProfile
  use InternalProfile
  implicit none

  real(kind=8), parameter :: ngamma = 1 ! number of gamma parameters
  real(kind=8), parameter :: gammamin = 5. ! minimum gamma
  real(kind=8), parameter :: gammainc = 1. ! gamma increment
  real(kind=8), parameter :: gammaincs = 0.1 ! small gamma increment (if gamma<1)
  real(kind=8), parameter :: nlnu = 1 ! number of gamma parameters
  real(kind=8), parameter :: lnuin = 6. ! viscous sublayer thickness
  real(kind=8), parameter :: lnuinc = 3. ! gamma increment
  logical, parameter :: couette = .TRUE.
  logical, parameter :: poiseuille = .FALSE.

  real(kind=8), parameter :: chi = 0.41 ! von Karman constant
  real(kind=8), parameter :: c = 5. ! integration constant

  real(kind=8) :: gamma, lnu, xfac, ufac, re, cd, cdexp, cdlin, reexp, relin
  integer :: k, i, j,idx
  character(len=256) :: filename, filecdp, filecdc

  nprint = 4
  integtype = 'Explicit'

  ! Lopp over gamma and lnu values
  do j=1,nlnu
    if (j==1) then
      lnu = lnuin
    else
      lnu = lnu + lnuinc
    endif

  write(filecdp,'(a3,1E9.2,a4)') 'cdp',-lnu,'.txt'
  idx=scan(trim(filecdp),' ')
  do while(idx.ne.0)
    filecdp(idx:idx)='-'
    idx=scan(trim(filecdp),' ')
  enddo
  open(unit=21,file=filecdp)

  write(filecdc,'(a3,1E9.2,a4)') 'cdc',-lnu,'.txt'
  idx=scan(trim(filecdc),' ')
  do while(idx.ne.0)
    filecdc(idx:idx)='-'
    idx=scan(trim(filecdc),' ')
  enddo
  open(unit=22,file=filecdc)

  do i=1,ngamma
    ! Compute gamma
    if (i==1) then
      gamma = gammamin
    else
      if (gamma<1.) then
      !if ((gamma>=5.).and.(gamma<7.)) then
      !if ((gamma>5.).and.(gamma<=7.)) then
        gamma = gamma + gammaincs
      else
        gamma = gamma + gammainc
      endif
    endif

    if (poiseuille) then

    ! Compute Poiseuille profile
    flowtype = 'Poiseuille'
    call ComputeProfile(gamma,lnu)
    !call ComputeInternalProfile(gamma,g(1))
    call ComputeExternalProfile(gamma)
    call ComputeLogProfile(gamma,g(1))

    ! save profile in normal coordinates
    write(filename,'(a1,2E9.2,a4)') 'p',-gamma,-lnu,'.txt'
    idx=scan(trim(filename),' ')
    do while(idx.ne.0)
      filename(idx:idx)='-'
      idx=scan(trim(filename),' ')
    enddo

    open(unit=10,file=filename)
    do k = 1,n
      write(10,'(6E18.10)') x(k),g(k),u(k),g(1)*x(k),ulog(k),uext(k)
    enddo
    do k = n-1,1,-1
      write(10,'(6E18.10)') 1.-x(k),-g(k),u(k),g(1)*x(k),ulog(k),uext(k)
    enddo
    close(10)

    ! save Cd and Reynolds
    if (gamma>0.) then
      re = g(1)*gamma*gamma/(chi*chi)
      cd = 2*chi*chi/(gamma*gamma)
      cdlin = 8./re
      reexp = re * 0.5
      relin = 4.*gamma*gamma/(chi*chi)
      cdexp = cdfct(reexp)
      write(21,'(7E18.10)') gamma,g(1),cd,re,cdexp,cdlin,relin
    endif

    if (gamma>1.) then
      ! save profile in reduced coordinates
      write(filename,'(a2,2E9.2,a4)') 'pr',-gamma,-lnu,'.txt'
      idx=scan(trim(filename),' ')
      do while(idx.ne.0)
        filename(idx:idx)='-'
        idx=scan(trim(filename),' ')
      enddo

      ufac = gamma / chi
      xfac = g(1) * gamma / chi

      x(:) = x(:) * xfac
      u(:) = u(:) * ufac
      ulog(:) = ulog(:) * ufac

      open(unit=10,file=filename)
      do k = 1,n
        write(10,'(4E18.10)') x(k),u(k),ulog(k),log(x(k))/chi+c
      enddo
      close(10)
    endif

    endif

    if (couette) then

    ! Compute Couette profile
    flowtype = 'Couette'
    call ComputeProfile(gamma,lnu)
    !call ComputeInternalProfile(gamma,g(1))
    call ComputeExternalProfile(gamma)
    call ComputeLogProfile(gamma,g(1))

    write(filename,'(a1,2E9.2,a4)') 'c',-gamma,-lnu,'.txt'
    idx=scan(trim(filename),' ')
    do while(idx.ne.0)
      filename(idx:idx)='-'
      idx=scan(trim(filename),' ')
    enddo

    ! save profile in normal coordinates
    open(unit=10,file=filename)
    do k = 1,n
      write(10,'(6E18.10)') x(k),g(k),u(k)-1.,g(1)*x(k)-1.,ulog(k)-1.,uext(k)-1.
    enddo
    do k = n-1,1,-1
      write(10,'(6E18.10)') 1.-x(k),g(k),1.-u(k),1.-g(1)*x(k),1.-ulog(k),1.-uext(k)
    enddo
    close(10)

    ! save Cd and Reynolds
    if (gamma>0.) then
      re = g(1)*gamma*gamma/(chi*chi)
      cd = chi*chi/(gamma*gamma)
      cdexp = 0.048*(re**-0.25)
      cdlin = 4./re
      write(22,'(6E18.10)') gamma,g(1),cd,re,cdexp,cdlin
    endif

    if (gamma>1.) then
      ! save profile in reduced coordinates
      write(filename,'(a2,2E9.2,a4)') 'cr',-gamma,-lnu,'.txt'
      idx=scan(trim(filename),' ')
      do while(idx.ne.0)
        filename(idx:idx)='-'
        idx=scan(trim(filename),' ')
      enddo

      ufac = gamma / chi
      xfac = g(1) * gamma / chi

      x(:) = x(:) * xfac
      u(:) = u(:) * ufac
      ulog(:) = ulog(:) * ufac

      open(unit=10,file=filename)
      do k = 1,n
        write(10,'(4E18.10)') x(k),u(k),ulog(k),log(x(k))/chi+c
      enddo
      close(10)
    endif

    endif

  enddo

  close(21)
  close(22)

  enddo

  contains

    function cdfct(reynolds)
    real(kind=8), intent(in) :: reynolds
    real(kind=8) :: cdfct

    real(kind=8) :: cc, ccp

    ccp = 0.
    cc = 0.048_8 * (reynolds ** -0.25_8)
    
    do while((abs(cc-ccp)/cc)>accuracy)
      ccp = cc
      cc = reynolds * sqrt( ccp * 0.5_8 )
      cc = c + log(cc) / chi
      cc = 1._8 / ( cc * cc * 0.5_8 )
    enddo

    cdfct = cc

    end function cdfct

end program bl_profile
