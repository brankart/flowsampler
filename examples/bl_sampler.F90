program bl_sampler
  use GeneralProfile
  use ProfileSampler
  use ensdam_storng
  implicit none

  real(kind=8), parameter :: ngamma = 20 ! number of gamma parameters
  real(kind=8), parameter :: gammamin = 5. ! minimum gamma
  real(kind=8), parameter :: lnu = 6. ! viscous sublayer thickness
  logical, parameter :: couette = .TRUE.
  logical, parameter :: poiseuille = .FALSE.

  real(kind=8), parameter :: chi = 0.41 ! von Karman constant
  real(kind=8), parameter :: c = 5. ! integration constant

  real(kind=8) :: gamma, xfac, ufac, re, cd, cdexp, cdlin, reexp
  integer :: k, i, j, idx, kx
  character(len=256) :: filename, filecdp, filecdc

  nprint = 2
  nres = 100
  maxiter = 1000
  nharmonics = 80
  !initmax=.TRUE.
  integtype = 'Explicit'
  costfac = 0.1
  pertstd = 0.5

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

  ! Lopp over gamma values
  do i=1,ngamma
    ! Compute gamma
    gamma = gammamin

    if (poiseuille) then

    ! Compute Poiseuille profile
    flowtype = 'Poiseuille'
    call SampleProfile(gamma,lnu)

    ! save profile in normal coordinates
    write(filename,'(a1,i4.4,2E9.2,a4)') 'p',i,-gamma,-lnu,'.txt'
    idx=scan(trim(filename),' ')
    do while(idx.ne.0)
      filename(idx:idx)='-'
      idx=scan(trim(filename),' ')
    enddo

    open(unit=10,file=filename)
    do k = 1,n
      write(10,'(3E18.10)') x(k),gsmp(k),usmp(k)
    enddo
    do k = n+1,2*n-1
      kx = 2*n-k
      write(10,'(3E18.10)') 1.-x(kx),gsmp(k),usmp(k)
    enddo
    close(10)

    ! save Cd and Reynolds
    if (gamma>0.) then
      re = gsmp(1)*gamma*gamma/(chi*chi)
      cd = chi*chi/(gamma*gamma)
      cdlin = 4./re
      reexp = re * 0.5
      cdexp = 0.048*(reexp**-0.25)
      cdexp = 0.5 * cdexp
      write(21,'(6E18.10)') gamma,gsmp(1),cd,re,cdexp,cdlin
    endif

    if (gamma>1.) then
      ! save profile in reduced coordinates
      write(filename,'(a2,i4.4,2E9.2,a4)') 'pr',i,-gamma,-lnu,'.txt'
      idx=scan(trim(filename),' ')
      do while(idx.ne.0)
        filename(idx:idx)='-'
        idx=scan(trim(filename),' ')
      enddo

      ufac = gamma / chi
      xfac = gsmp(1) * gamma / chi

      open(unit=10,file=filename)
      do k = 1,n
        write(10,'(4E18.10)') x(k)*xfac,usmp(k)*ufac
      enddo
      close(10)
    endif

    endif

    if (couette) then

    ! Compute Couette profile
    flowtype = 'Couette'
    call SampleProfile(gamma,lnu)

    write(filename,'(a1,i4.4,2E9.2,a4)') 'c',i,-gamma,-lnu,'.txt'
    idx=scan(trim(filename),' ')
    do while(idx.ne.0)
      filename(idx:idx)='-'
      idx=scan(trim(filename),' ')
    enddo

    ! save profile in normal coordinates
    open(unit=10,file=filename)
    do k = 1,n
      write(10,'(3E18.10)') x(k),gsmp(k),usmp(k)-1.
    enddo
    do k = n+1,2*n-1
      kx = 2*n-k
      write(10,'(3E18.10)') 1.-x(kx),gsmp(k),usmp(k)-1.
    enddo
    close(10)

    ! save Cd and Reynolds
    if (gamma>0.) then
      re = gsmp(1)*gamma*gamma/(chi*chi)
      cd = chi*chi/(gamma*gamma)
      cdexp = 0.048*(re**-0.25)
      cdlin = 4./re
      write(22,'(6E18.10)') gamma,gsmp(1),cd,re,cdexp,cdlin
    endif

    if (gamma>1.) then
      ! save profile in reduced coordinates
      write(filename,'(a2,i4.4,2E9.2,a4)') 'cr',i,-gamma,-lnu,'.txt'
      idx=scan(trim(filename),' ')
      do while(idx.ne.0)
        filename(idx:idx)='-'
        idx=scan(trim(filename),' ')
      enddo

      ufac = gamma / chi
      xfac = gsmp(1) * gamma / chi

      open(unit=10,file=filename)
      do k = 1,n
        write(10,'(4E18.10)') x(k)*xfac,usmp(k)*ufac
      enddo
      close(10)
    endif

    endif

  enddo

  close(21)
  close(22)

end program bl_sampler
