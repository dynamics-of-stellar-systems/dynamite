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

module initial_parameters

    use numeric_kinds

    implicit none

    integer (kind=i4b), public :: ngauss_mge
    real (kind=dp), dimension(:), allocatable, public :: surf_km,psi_obs
    real (kind=dp), dimension(:), allocatable, public :: sigobs_km, qobs
    real (kind=dp), public :: theta_view, phi_view, xmbh, softl_km,psi_view
    real (kind=dp), public :: conversion_factor, totalmass
  
    ! DM parameters
    integer (kind=i4b), public :: n_dmparam,dm_profile_type
    real (kind=dp), dimension(:), allocatable, public :: dmparam
    real (kind=dp), public :: gamma_var

    ! Global parameters for the orbital starting points
    ! nEner = # energies, nI2 = # I2, nI3 = # I3 
    ! rLogMin and rLogMax : log10 of min and max radius in KM
    ! orbit_dithering is the amount of dithering
    integer(kind=i4b), public :: nEner,nI2,nI3,orbit_dithering
    real   (kind=dp ), public :: rLogMin,rLogMax
    !1 Solar mass = 1.98892 x 10^30 kg [Wikipedia]
    !1 AU = 1.4959787068d8 km [IAU 1976]
    !1 pc = 1.4959787068d8*(648d3/!dpi) km = 3.0856776e+13 km

    ! the NIST recommended constant of gravity:
    !G = 6.67428 (+/- 0.00067) x 10^-11 m^3 kg^-1 s^-2
    ![http://physics.nist.gov/cgi-bin/cuu/Value?bg|search_for=newtonian]

    !           G = 1.33381d11 km^3/(s^2 Msun)
    !  critical density rho_crit = 3H^2/8piG
    real (kind=dp), parameter, public :: &
         grav_const_km=6.67428e-11_dp*1.98892e30_dp/1e9_dp, &
         parsec_km = 1.4959787068d8*(648d3/pi_d), &
	     rho_crit = (3.0_dp*(7.3d-5/parsec_km)**2)/(8.0_dp*pi_d*grav_const_km)


    private ! default private
    public :: iniparam
 
contains

    subroutine iniparam()

        character (len=256) :: infil
        real (kind=dp), dimension(:), allocatable :: surf_pc, sigobs_arcsec
        real (kind=dp) :: distance, upsilon, softl_arcsec,dm_fraction, concentration
        real (kind=dp) :: darkmass,dum,dumdum
        integer (kind=i4b) :: j

        print*,"Gravitational Constant in km^3/(s^2 Msun)",grav_const_km
        print*,"parsec in km",parsec_km
        print*,""
        print*, "Give the name of the file with input parameters"
        read (unit=*, fmt="(a256)") infil
        open (unit=13,FILE=infil,status="old",action="read",position="rewind")
        ! read number of Gaussians
        read (unit=13, fmt=*) ngauss_mge

        allocate (surf_pc(ngauss_mge), sigobs_arcsec(ngauss_mge), &
             qobs(ngauss_mge), surf_km(ngauss_mge), sigobs_km(ngauss_mge)&
             ,psi_obs(ngauss_mge))

        ! The MGE observed parameters of the gaussians.
        ! For every gaussian reads:
        ! 1) Central surface brightness: I_j      (in units of L_sun/pc^2)
        ! 2) Dispersion:                 sigma'_j (in units of arcsec)
        ! 3) Observed flattening:        q'_j     (dimenstionless)
        ! 4) Offset Psi viewing angle    psi_j    (unit of degrees)
        ! NB: The gaussians should be arranged in order of increasing sigma
        do j=1,ngauss_mge
          read (unit=13, fmt=*) surf_pc(j),sigobs_arcsec(j),qobs(j),psi_obs(j)
        end do

        ! read the distance (in units of Mpc)
        read (unit=13, fmt=*) distance
        ! read the viewing angles (in units of degrees)
        read (unit=13, fmt=*) theta_view,phi_view,psi_view
        ! read the mass-to-light ratio (in units of M_sun/L_sun)
        ! NB: in the same photometric band as the surface brightness
        read (unit=13, fmt=*) upsilon
        ! read the black hole mass (in units of M_sun)
        read (unit=13, fmt=*) xmbh
        ! read the softening length of the black hole (in units of arcsec)
        read (unit=13, fmt=*) softl_arcsec
        ! read the # of energies, rlogmin(arcsec) and rlogmax(arcsec)
        read (unit=13, fmt=*) nEner,rlogmin,rlogmax
        ! read the number of orbital range in I2
        read (unit=13, fmt=*) nI2
        ! read the number of orbital range in I3
        read (unit=13, fmt=*) nI3
        ! The number of orbital dithering
        read (unit=13, fmt=*) orbit_dithering

		!added by AW   ! The parameters of the NFW halo: rho_c in solarmass/pc^3 and r_c in arcsec
		!added by JJA, type=5 for gNFW with three inputs: c, DM virial mass in solarmass, and gamma 
        read (unit=13, fmt=*) dm_profile_type,n_dmparam
        allocate(dmparam(n_dmparam))
        read (unit=13, fmt=*) dmparam(1:n_dmparam)


        close (unit=13)

        ! apply dithering to the number of orbits
        Nener = Nener * orbit_dithering
        nI2   = nI2   * orbit_dithering
        nI3   = nI3   * orbit_dithering

        ! conversion factor from arcsec to km
        ! distance is in mpc
        !conversion_factor = distance*1.0e6_dp*1.49598e8_dp
        conversion_factor = distance*1.0e6_dp*tan(pi_d/(648d3))*parsec_km

        print*,'  * Conversion from arcsec to km:',conversion_factor    
	    print*,'The softening lenght of the black hole is ', & 
	           softl_arcsec*conversion_factor/2.95_dp/xmbh,' Schwarschild radii'

     
        ! surface brightness in L_sun/km^2
        surf_km(:) = surf_pc(:)/parsec_km**2
        ! from surface brightness to surface density in M_sun/km^2
        ! NB: assuming constant mass-to-light ratio
        surf_km(:) = surf_km(:) * upsilon
        ! dispersion in km
        sigobs_km(:) = sigobs_arcsec(:)*conversion_factor
        ! softening length in km
        softl_km = softl_arcsec*conversion_factor

        ! Apply offset to psi_obs
        psi_obs(:) = psi_obs(:) + psi_view

        ! viewing angles in radians 
        theta_view = theta_view * (pi_d/180.0_dp)
        phi_view   = phi_view   * (pi_d/180.0_dp)
        psi_obs    = psi_obs    * (pi_d/180.0_dp)
        psi_view   = psi_view   * (pi_d/180.0_dp)
        ! convert rlog* from log(arcsec) to log(km)
        rLogMin = rLogMin + log10(conversion_factor) 
        rLogMax = rLogMax + log10(conversion_factor)  
        
        ! Calculate the total mass of the model
        totalmass = twopi_d*sum(surf_km(:)*qobs(:)*sigobs_km(:)**2)
         
      end subroutine iniparam
 !     print*,'  * read in initial parameters finished    *'
end module initial_parameters

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
