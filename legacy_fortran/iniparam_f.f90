!#######################################################################
!
! Read MGE physical parameters from files.
!
! This Fortran 90 code is compatible with the F language subset
! http://www.fortran.com/imagine1
!
! HISTORY:
!
! V1.0 Written and tested triaxial extension of original axisymmetric
!   version by M. Cappellari. G. van de Ven, Leiden, JUNE/2004
! V2.0 Refactored and simplified by Remco van den Bosch JULY/2004
! V2.1 Added generic dark halo parameters. RvdB March/2010
!
!#######################################################################

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The code related to the bar/disk decomposition is in the subroutine iniparam_bar.
! The observational parameters (from MGE) and intrinsic parameters (p,q,u)
! are sorted and calculated for disk and bar separately, 
! then it is superposed to build the whole galaxy.
!
! The first line in parameters.in is changed to be included 4 parameters
! e.g (12   1   5   7)
! The first one is the number of all Gaussians,
! the second one is to choose whether you need one component.
! deprojection (=0  it is the old version of code) or you need 
! two-component deprojection (=1  it is the new deprojection for barred galaxies).
! In the case of two-component deprojection, the third and fourth parameters 
! are the number of Gaussians of the barred bulge and the disk, respectively.
! Put these two zero if you choose old deeprojection. 
! The last line in parameters.in is a new added parameter which is pattern speed (in unit of km/s/kpc) 
! Note that if you use a positive Omega kinematic maps should be set to be counter-clockwise
!
! adapted from Behzad Tahmasebzadeh's code, July 2023
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


module initial_parameters

    use numeric_kinds

    implicit none

    integer (kind=i4b), public :: ngauss_mge, decmposed, ngaus_bulge, ngaus_disk
    real (kind=dp), dimension(:), allocatable, public :: surf_km,   psi_obs,   sigobs_km,   qobs
    real (kind=dp), dimension(:), allocatable, public :: surf_km_d, psi_obs_d, sigobs_km_d, qobs_d
    real (kind=dp), dimension(:), allocatable, public :: surf_km_b, psi_obs_b, sigobs_km_b, qobs_b


    real(kind=dp), public :: theta_view, phi_view, xmbh, softl_km, psi_view
    real(kind=dp), public :: conversion_factor, totalmass

    real (kind=dp), public :: Omega = 0.0_dp ! set to 0.0_dp by default so that non-bar model is assumed

    ! DM parameters
    integer(kind=i4b), public :: n_dmparam, dm_profile_type
    real(kind=dp), dimension(:), allocatable, public :: dmparam
    real(kind=dp), public :: gamma_var

    ! Global parameters for the orbital starting points
    ! nEner = # energies, nI2 = # I2, nI3 = # I3
    ! rLogMin and rLogMax : log10 of min and max radius in KM
    ! orbit_dithering is the amount of dithering
    integer(kind=i4b), public :: nEner, nI2, nI3, orbit_dithering
    real(kind=dp), public :: rLogMin, rLogMax
    !1 Solar mass = 1.98892 x 10^30 kg [Wikipedia]
    !1 AU = 1.4959787068d8 km [IAU 1976]
    !1 pc = 1.4959787068d8*(648d3/!dpi) km = 3.0856776e+13 km

    ! the NIST recommended constant of gravity:
    !G = 6.67428 (+/- 0.00067) x 10^-11 m^3 kg^-1 s^-2
    ![http://physics.nist.gov/cgi-bin/cuu/Value?bg|search_for=newtonian]

    !           G = 1.33381d11 km^3/(s^2 Msun)
    !  critical density rho_crit = 3H^2/8piG
    real(kind=dp), parameter, public :: &
        grav_const_km = 6.67428e-11_dp*1.98892e30_dp/1e9_dp, &
        parsec_km = 1.4959787068d8*(648d3/pi_d), &
        rho_crit = (3.0_dp*(7.3d-5/parsec_km)**2)/(8.0_dp*pi_d*grav_const_km)

    private ! default private
    public :: iniparam, iniparam_bar

contains

    subroutine iniparam()

        character(len=256) :: infil
        real(kind=dp), dimension(:), allocatable :: surf_pc, sigobs_arcsec
        real(kind=dp) :: distance, upsilon, softl_arcsec, dm_fraction, concentration
        real(kind=dp) :: darkmass, dum, dumdum
        integer(kind=i4b) :: j

        print *, "Gravitational Constant in km^3/(s^2 Msun)", grav_const_km
        print *, "parsec in km", parsec_km
        print *, ""
        print *, "Give the name of the file with input parameters"
        read (unit=*, fmt="(a256)") infil
        open (unit=13, FILE=infil, status="old", action="read", position="rewind")
        ! read number of Gaussians
        read (unit=13, fmt=*) ngauss_mge

        allocate (surf_pc(ngauss_mge), sigobs_arcsec(ngauss_mge), &
                  qobs(ngauss_mge), surf_km(ngauss_mge), sigobs_km(ngauss_mge) &
                  , psi_obs(ngauss_mge))

        ! The MGE observed parameters of the gaussians.
        ! For every gaussian reads:
        ! 1) Central surface brightness: I_j      (in units of L_sun/pc^2)
        ! 2) Dispersion:                 sigma'_j (in units of arcsec)
        ! 3) Observed flattening:        q'_j     (dimenstionless)
        ! 4) Offset Psi viewing angle    psi_j    (unit of degrees)
        ! NB: The gaussians should be arranged in order of increasing sigma
        do j = 1, ngauss_mge
            read (unit=13, fmt=*) surf_pc(j), sigobs_arcsec(j), qobs(j), psi_obs(j)
        end do

        ! read the distance (in units of Mpc)
        read (unit=13, fmt=*) distance
        ! read the viewing angles (in units of degrees)
        read (unit=13, fmt=*) theta_view, phi_view, psi_view
        ! read the mass-to-light ratio (in units of M_sun/L_sun)
        ! NB: in the same photometric band as the surface brightness
        read (unit=13, fmt=*) upsilon
        ! read the black hole mass (in units of M_sun)
        read (unit=13, fmt=*) xmbh
        ! read the softening length of the black hole (in units of arcsec)
        read (unit=13, fmt=*) softl_arcsec
        ! read the # of energies, rlogmin(arcsec) and rlogmax(arcsec)
        read (unit=13, fmt=*) nEner, rlogmin, rlogmax
        ! read the number of orbital range in I2
        read (unit=13, fmt=*) nI2
        ! read the number of orbital range in I3
        read (unit=13, fmt=*) nI3
        ! The number of orbital dithering
        read (unit=13, fmt=*) orbit_dithering

        !added by AW   ! The parameters of the NFW halo: rho_c in solarmass/pc^3 and r_c in arcsec
        !added by JJA, type=5 for gNFW with three inputs: c, DM virial mass in solarmass, and gamma
        read (unit=13, fmt=*) dm_profile_type, n_dmparam
        allocate (dmparam(n_dmparam))
        read (unit=13, fmt=*) dmparam(1:n_dmparam)

        close (unit=13)

        ! apply dithering to the number of orbits
        Nener = Nener*orbit_dithering
        nI2 = nI2*orbit_dithering
        nI3 = nI3*orbit_dithering

        ! conversion factor from arcsec to km
        ! distance is in mpc
        !conversion_factor = distance*1.0e6_dp*1.49598e8_dp
        conversion_factor = distance*1.0e6_dp*tan(pi_d/(648d3))*parsec_km

        print *, '  * Conversion from arcsec to km:', conversion_factor
        print *, 'The softening lenght of the black hole is ', &
            softl_arcsec*conversion_factor/2.95_dp/xmbh, ' Schwarschild radii'

        ! surface brightness in L_sun/km^2
        surf_km(:) = surf_pc(:)/parsec_km**2
        ! from surface brightness to surface density in M_sun/km^2
        ! NB: assuming constant mass-to-light ratio
        surf_km(:) = surf_km(:)*upsilon
        ! dispersion in km
        sigobs_km(:) = sigobs_arcsec(:)*conversion_factor
        ! softening length in km
        softl_km = softl_arcsec*conversion_factor

        ! Apply offset to psi_obs
        psi_obs(:) = psi_obs(:) + psi_view

        ! viewing angles in radians
        theta_view = theta_view*(pi_d/180.0_dp)
        phi_view = phi_view*(pi_d/180.0_dp)
        psi_obs = psi_obs*(pi_d/180.0_dp)
        psi_view = psi_view*(pi_d/180.0_dp)
        ! convert rlog* from log(arcsec) to log(km)
        rLogMin = rLogMin + log10(conversion_factor)
        rLogMax = rLogMax + log10(conversion_factor)

        ! Calculate the total mass of the model
        totalmass = twopi_d*sum(surf_km(:)*qobs(:)*sigobs_km(:)**2)

    end subroutine iniparam
    !     print*,'  * read in initial parameters finished    *'

    subroutine iniparam_bar()

        character(len=256) :: infil
        real (kind=dp), dimension(:), allocatable :: surf_pc, sigobs_arcsec        ! (BT)
        real (kind=dp), dimension(:), allocatable :: surf_pc_d, sigobs_arcsec_d    ! (BT)
        real (kind=dp), dimension(:), allocatable :: surf_pc_b, sigobs_arcsec_b    ! (BT)
        real(kind=dp) :: distance, upsilon, softl_arcsec, dm_fraction, concentration
        real(kind=dp) :: darkmass, dum, dumdum
        integer(kind=i4b) :: j

        print *, "Gravitational Constant in km^3/(s^2 Msun)", grav_const_km
        print *, "parsec in km", parsec_km
        print *, ""
        print *, "Give the name of the file with input parameters"
        read (unit=*, fmt="(a256)") infil
        open (unit=13, FILE=infil, status="old", action="read", position="rewind")
        ! read number of Gaussians and bulge/disk params
        read (unit=13, fmt=*) ngauss_mge, decmposed, ngaus_bulge, ngaus_disk                  ! (BT)
        print*, "Gaussians(n_totall, n_bulge, n_disk) =",ngauss_mge, ngaus_bulge, ngaus_disk  ! (BT)

        allocate ( surf_pc(ngauss_mge)   , sigobs_arcsec(ngauss_mge) ,&
        qobs(ngauss_mge)   ,       surf_km(ngauss_mge) ,&
        psi_obs(ngauss_mge)   ,     sigobs_km(ngauss_mge) ,&
        surf_pc_d(ngaus_disk)   , sigobs_arcsec_d(ngaus_disk) ,&
        qobs_d(ngaus_disk)   ,       surf_km_d(ngaus_disk) ,&
        psi_obs_d(ngaus_disk)   ,     sigobs_km_d(ngaus_disk) ,&
        surf_pc_b(ngaus_bulge)   ,  sigobs_arcsec_b(ngaus_bulge) ,&
        qobs_b(ngaus_bulge)   ,        surf_km_b(ngaus_bulge) ,&
        psi_obs_b(ngaus_bulge)   ,      sigobs_km_b(ngaus_bulge) )     ! (BT)

        ! The MGE observed parameters of the gaussians.
        ! For every gaussian reads:
        ! 1) Central surface brightness: I_j      (in units of L_sun/pc^2)
        ! 2) Dispersion:                 sigma'_j (in units of arcsec)
        ! 3) Observed flattening:        q'_j     (dimenstionless)
        ! 4) Offset Psi viewing angle    psi_j    (unit of degrees)
        ! NB: The gaussians should be arranged in order of increasing sigma
        if (decmposed == 0) then      ! (BT)  Deprojection same as before
           do j = 1, ngauss_mge
              read (unit=13, fmt=*) surf_pc(j), sigobs_arcsec(j), qobs(j), psi_obs(j)
           end do
        else if (decmposed == 1) then ! (BT) Deprojection of barred galaxy
           do j=1 ,ngaus_bulge
              read (unit=13, fmt=*) surf_pc_b(j),sigobs_arcsec_b(j),qobs_b(j),psi_obs_b(j)
              print*, "bulge",surf_pc_b(j),sigobs_arcsec_b(j),qobs_b(j),psi_obs_b(j)
           end do
           
           do j=1, ngaus_disk          ! (BT)
              read (unit=13, fmt=*) surf_pc_d(j),sigobs_arcsec_d(j),qobs_d(j),psi_obs_d(j)
              print*, "disk", surf_pc_d(j),sigobs_arcsec_d(j),qobs_d(j),psi_obs_d(j)
           end do
           
        else
           print*, "choose 0 or 1 for decomposed parameter"
        end if

        ! read the distance (in units of Mpc)
        read (unit=13, fmt=*) distance
        ! read the viewing angles (in units of degrees)
        read (unit=13, fmt=*) theta_view, phi_view, psi_view
        ! read the mass-to-light ratio (in units of M_sun/L_sun)
        ! NB: in the same photometric band as the surface brightness
        read (unit=13, fmt=*) upsilon
        ! read the black hole mass (in units of M_sun)
        read (unit=13, fmt=*) xmbh
        ! read the softening length of the black hole (in units of arcsec)
        read (unit=13, fmt=*) softl_arcsec
        ! read the # of energies, rlogmin(arcsec) and rlogmax(arcsec)
        read (unit=13, fmt=*) nEner, rlogmin, rlogmax
        ! read the number of orbital range in I2
        read (unit=13, fmt=*) nI2
        ! read the number of orbital range in I3
        read (unit=13, fmt=*) nI3
        ! The number of orbital dithering
        read (unit=13, fmt=*) orbit_dithering

        !added by AW   ! The parameters of the NFW halo: rho_c in solarmass/pc^3 and r_c in arcsec
        !added by JJA, type=5 for gNFW with three inputs: c, DM virial mass in solarmass, and gamma
        read (unit=13, fmt=*) dm_profile_type, n_dmparam
        allocate (dmparam(n_dmparam))
        read (unit=13, fmt=*) dmparam(1:n_dmparam)
        read (unit=13, fmt=*) Omega   ! (BT) reading pattern speed
        close (unit=13)

        ! apply dithering to the number of orbits
        Nener = Nener*orbit_dithering
        nI2 = nI2*orbit_dithering
        nI3 = nI3*orbit_dithering

        ! conversion factor from arcsec to km
        ! distance is in mpc
        !conversion_factor = distance*1.0e6_dp*1.49598e8_dp
        conversion_factor = distance*1.0e6_dp*tan(pi_d/(648d3))*parsec_km

        print *, '  * Conversion from arcsec to km:', conversion_factor
        print *, 'The softening lenght of the black hole is ', &
            softl_arcsec*conversion_factor/2.95_dp/xmbh, ' Schwarschild radii'

        ! Omega km/s/kpc to km/s/km
        Omega =  Omega  * (1.0_dp/30856776e9_dp)

        ! surface brightness in L_sun/km^2
        if (decmposed == 0) then          ! (BT)
           surf_km(:) = surf_pc(:)/parsec_km**2      
        else if (decmposed == 1) then
           surf_km_d(:) = surf_pc_d(:)/parsec_km**2                           
           surf_km_b(:) = surf_pc_b(:)/parsec_km**2
        else
           print*, "choose 0 or 1 for decomposed parameter"
        end if

        ! from surface brightness to surface density in M_sun/km^2
        ! NB: assuming constant mass-to-light ratio
        if (decmposed == 0) then           ! (BT)
           surf_km(:) = surf_km(:) * upsilon      
        else if (decmposed == 1) then
           surf_km_d(:) = surf_km_d(:) * upsilon                                 
           surf_km_b(:) = surf_km_b(:) * upsilon
        else
           print*, "choose 0 or 1 for decomposed parameter"
        end if


        ! dispersion in km
        if (decmposed == 0) then            ! (BT)
           sigobs_km(:) = sigobs_arcsec(:)*conversion_factor     
        else if (decmposed == 1) then
           sigobs_km_d(:) = sigobs_arcsec_d(:)*conversion_factor
           sigobs_km_b(:) = sigobs_arcsec_b(:)*conversion_factor
        else
           print*, "choose 0 or 1 for decomposed parameter"
        end if

        ! softening length in km
        softl_km = softl_arcsec*conversion_factor

        ! Apply offset to psi_obs
        if (decmposed == 0) then              ! (BT)
           psi_obs(:) = psi_obs(:) +  psi_view     
        else if (decmposed == 1) then
           psi_obs_d(:) = psi_obs_d(:)  + 90     ! psi_view_disk = pi/2
           psi_obs_b(:) = psi_obs_b(:)  + psi_view
        else
           print*, "choose 0 or 1 for decomposed parameter"
        end if

        ! viewing angles in radians
        theta_view = theta_view*(pi_d/180.0_dp)
        phi_view = phi_view*(pi_d/180.0_dp)

        if (decmposed == 0) then         ! (BT)
           psi_obs  = psi_obs    * (pi_d/180.0_dp)     
        else if (decmposed == 1) then
           psi_obs_d  = psi_obs_d    * (pi_d/180.0_dp)    ! (BT) psi_view_disk = pi/2. the bar is aligned in the disk
           psi_obs_b  = psi_obs_b    * (pi_d/180.0_dp)
        else
           print*, "choose 0 or 1 for decomposed parameter"
        end if

        psi_view = psi_view*(pi_d/180.0_dp)
        ! convert rlog* from log(arcsec) to log(km)
        rLogMin = rLogMin + log10(conversion_factor)
        rLogMax = rLogMax + log10(conversion_factor)


        if (decmposed == 1) then     ! (BT) composition of the bar and the disk to make whole the galaxy
           surf_km    = [surf_km_d , surf_km_b]                              
           qobs       = [qobs_d , qobs_b]                                    
           sigobs_km  = [sigobs_km_d , sigobs_km_b]                          
           psi_obs    = [psi_obs_d , psi_obs_b]                              
        end if


        ! Calculate the total mass of the model
        totalmass = twopi_d*sum(surf_km(:)*qobs(:)*sigobs_km(:)**2)

    end subroutine iniparam_bar


end module initial_parameters

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
