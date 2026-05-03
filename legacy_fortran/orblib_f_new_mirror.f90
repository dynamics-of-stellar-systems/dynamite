!######################################################################
!
! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! Sterrewacht Leiden, The Netherlands
!
! HISTORY:
!
! V1.4: MC Verification of intrinsic velocity moments calculation.
!     Fixed significant bug in mergrid_store. Slightly revised
!     implementation
!     to make it easier to check. Independent test of the results
!      still needed
! V2.0: RvdB. Fork for the Triaxial orbit library code.
! V2.0.1 : Changed quadrant grid bins to have more bins.
! V2.0.2 : fix 10**rlogmax in bin setup of qgrid_setup
!          add internal moments
! V2.0.3 : change the zeroth moment grid
! V2.0.4 : add zero psf routine
!          fixed projection
! V3.0.0 : Make intrinisic grid size recipe to be able to numerically
!          Fit the intrinsic mass.
!     RvdB, Leiden, oktober/2005
!
!######################################################################
! $Id: orblib_f.f90,v 1.3 2011/10/25 08:48:45 bosch Exp $

! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! Sterrewacht Leiden, The Netherlands

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Note on rotating bar code:
! The modified or added lines are commented by ! (BT)
! In this file, the modified parts are:
! 1- Integration of orbits in a rotating frame in case of (Omega != 0).
! 2- Applying 4-fold symmetry instead of 8-fold symmetry mirroring and symmetrizing if (Omega != 0).
! 3- if (Omega != 0), Sorting information of each orbit e.g. Circularity and ... during integration.
! created file called  _orb_info.out for each library
!
! adapted from code originally by Behzad Tahmasebzadeh, July 2023
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Added proper motions 2d histograms
! Removed aperture_polygon
! Thomas I. Maindl, November 2025
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module random_gauss_generator
    use numeric_kinds
    implicit none
    private

    ! Seeds the NR random generator
    public :: random_gauss_seed

    ! F version of the gaussian random generator
    ! generate gaussians for a n*m*o array with a width
    public :: random_gauss

    !Generate one 2d gaussian deviate with sigma "width"
    private :: gaussdev

contains

! adapted from NR2.
! Internal computation in SP for speed
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine gaussdev(x)
        real(kind=dp), intent(out), dimension(:) :: x
        !----------------------------------------------------------------------
        real(kind=sp), dimension(2) :: v
        real(kind=sp) :: rsq
        real(kind=dp) :: ran1

        do
            v(1) = ran1(1)
            v(2) = ran1(1)
            ! call random_number(v)
            v = 2.0_sp*v - 1.0_sp
            rsq = sum(v**2)
            if (rsq > 0.0_sp .and. rsq < 1.0_sp) exit
        end do
        x = v*sqrt(-2.0_sp*log(rsq)/rsq)

    end subroutine gaussdev

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine random_gauss_seed()
        !----------------------------------------------------------------------
        logical, save :: initialized = .false.

        print *, "  * Seeding native Random generator"
        if (.not. initialized) then
            ! START reproducible orbit library
            ! uncomment the following line for stochastic orbit library creation
            !call random_seed()
            ! END reproducible orbit library
            print *, "  * Internal Compiler random functions needs to be checked."
            initialized = .true.
        end if

    end subroutine random_gauss_seed

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine random_gauss(t)
        real(kind=dp), dimension(:, :), intent(out) :: t
        !----------------------------------------------------------------------
        integer(kind=i4b)                              :: k

        do k = 1, size(t, 1)
            call gaussdev(t(k, :))
        end do

    end subroutine random_gauss

end module random_gauss_generator

!######################################################################
!######################################################################
!######################################################################

! $Id: orblib_f.f90,v 1.3 2011/10/25 08:48:45 bosch Exp $
! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! Sterrewacht Leiden, The Netherlands

module integrator
    ! Integrator module. Integrates a point in the MGE potential
    use numeric_kinds
    implicit none
    private

    public  :: integrator_integrate

    public  :: integrator_setup, integrator_set_current, integrator_setup_bar

    public  :: integrator_stop, integrator_find_orbtype

    public  :: integrator_setup_write, integrator_write

    private :: real_integrator

    ! current orbit
    integer(kind=i4b), public :: integrator_current

    ! contains the number of points generated by the integrator
    integer(kind=i4b), public :: integrator_points

    ! Starting point to begin with
    integer(kind=i4b), private :: integrator_start

    ! Number of different orbits
    integer(kind=i4b), private :: integrator_number

    ! Number of orbits to integrate for each orbit
    real(kind=dp), private :: integrator_n_orbits

    ! accuracy of integrator
    real(kind=dp), private :: integrator_accuracy

    ! Number of orbit ditherings
    integer(kind=i4b), public :: integrator_dithering

    ! is the current orbital set not regularizble?
    integer(kind=i4b), public :: totalnotregularizable

    ! temporary pos and velocity array for the dense output of the integrator
    real(kind=dp), private, allocatable, dimension(:, :) :: vel_t, pos_t

    ! functions and variables copied from original orblib
    private :: integrator_whichorbit, derivs, SOLOUT
    public :: ini_integ

    ! information about the initial conditions for the orbit integration,
    ! and the grid in integral space.
    real(kind=dp), private, allocatable, dimension(:) :: xini, yini, zini
    real(kind=dp), private, allocatable, dimension(:) :: vxini, vzini
    real(kind=dp), public, allocatable, dimension(:) :: vyini, rcirc
    real(kind=dp), private, allocatable, dimension(:) :: vcirc, tcirc
    integer(kind=i4b), private, allocatable, dimension(:) :: regurizable

    integer(kind=i4b), private, allocatable, dimension(:) :: gEner, gI2, gI3

    integer(kind=i4b), public ::  nEner, nI2, nI3

    integer(kind=i4b), dimension(:), allocatable, private :: integrator_orbittypes
    real(kind=dp), private, allocatable, dimension(:, :)  ::  integrator_moments

    ! store the previously integrated orbit
    real(kind=dp), private, allocatable, dimension(:, :) :: pos_old, vel_old

contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine integrator_set_current(orbit)
        integer(kind=i4b), intent(in):: orbit
        !----------------------------------------------------------------------
        integrator_current = orbit
        Print *, "  * Starting at orbit", orbit + 1
        if (orbit < 0) stop " Can't start at that orbit"
        if (orbit > (nEner*nI2*nI3/integrator_dithering**3) - 1) &
            stop " Not so many orbits"

    end subroutine integrator_set_current

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine integrator_setup()
        use initial_parameters, only: iniparam, orbit_dithering
        !use triaxpotent, only : tp_setup
        use interpolpot, only: ip_setup
        !----------------------------------------------------------------------
        integer(kind=i4b) :: ndith3
        print *, "  ** Setting up integrator module"
        print *, "  * Calling MGE setup"
        call iniparam()
        print *, "  * Calling ip_setup"
        call ip_setup()
        call ini_integ()
        print *, "  * How many orbits should be integrated?"
        read *, integrator_n_orbits
        print *, "    ", integrator_n_orbits
        if (integrator_n_orbits < 1) stop " Too few orbits"
        print *, "  * How many points should be generated per starting point?"
        read *, integrator_points
        print *, "    ", integrator_points
        integrator_dithering = orbit_dithering
        ndith3 = integrator_dithering**3
        allocate (integrator_orbittypes(ndith3))
        allocate (integrator_moments(5, ndith3))
        if (integrator_points < 1) stop " Too few points"
        print *, "  * At which starting point should be started?"
        read *, integrator_start
        print *, "    ", integrator_start
        call integrator_set_current(integrator_start - 1)
        print *, "  * How many starting points should be integrated?"
        read *, integrator_number
        print *, "    ", integrator_number
        if (integrator_number == -1) integrator_number = &
            (nEner*nI2*nI3/integrator_dithering**3)

        if (integrator_number < 1) stop " To few starting points"
        if (integrator_number > &
            (nEner*nI2*nI3/integrator_dithering**3)) &
            & stop " Too many orbits in total"
        print *, "  * How great should te accuracy be of the integrator?"
        read *, integrator_accuracy
        print *, "    ", integrator_accuracy
        if (integrator_accuracy < 0 .or. 0.5 < integrator_accuracy) &
          & stop " wrong accuracy"

        print *, "  ** integrator module setup finished"

    end subroutine integrator_setup

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine integrator_setup_bar()
        use initial_parameters, only: iniparam_bar, orbit_dithering
        !use triaxpotent, only : tp_setup
        use interpolpot, only: ip_setup_bar
        !----------------------------------------------------------------------
        integer(kind=i4b) :: ndith3
        print *, "  ** Setting up integrator module"
        print *, "  * Calling MGE setup"
        call iniparam_bar()
        print *, "  * Calling ip_setup_bar"
        call ip_setup_bar()
        call ini_integ()
        print *, "  * How many orbits should be integrated?"
        read *, integrator_n_orbits
        print *, "    ", integrator_n_orbits
        if (integrator_n_orbits < 1) stop " Too few orbits"
        print *, "  * How many points should be generated per starting point?"
        read *, integrator_points
        print *, "    ", integrator_points
        integrator_dithering = orbit_dithering
        ndith3 = integrator_dithering**3
        allocate (integrator_orbittypes(ndith3))
        allocate (integrator_moments(5, ndith3))
        if (integrator_points < 1) stop " Too few points"
        print *, "  * At which starting point should be started?"
        read *, integrator_start
        print *, "    ", integrator_start
        call integrator_set_current(integrator_start - 1)
        print *, "  * How many starting points should be integrated?"
        read *, integrator_number
        print *, "    ", integrator_number
        if (integrator_number == -1) integrator_number = &
            (nEner*nI2*nI3/integrator_dithering**3)

        if (integrator_number < 1) stop " To few starting points"
        if (integrator_number > &
            (nEner*nI2*nI3/integrator_dithering**3)) &
            & stop " Too many orbits in total"
        print *, "  * How great should te accuracy be of the integrator?"
        read *, integrator_accuracy
        print *, "    ", integrator_accuracy
        if (integrator_accuracy < 0 .or. 0.5 < integrator_accuracy) &
          & stop " wrong accuracy"

        print *, "  ** integrator module setup finished"

    end subroutine integrator_setup_bar

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine integrator_stop()
        !----------------------------------------------------------------------
        if (allocated(xini)) then
            deallocate (xini, yini, zini, vxini, vyini, vzini, gener, gi2, gi3, &
                        vcirc, tcirc, rcirc, regurizable, vel_old, pos_old)
        end if

        if (allocated(pos_old)) then
            deallocate (pos_old, vel_old)
            deallocate (pos_t, vel_t)
        end if

        close (unit=30)
    end subroutine integrator_stop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine integrator_integrate(pos, vel, otype, done, first, alldone)
      use initial_parameters    , only : conversion_factor , Omega     ! (BT)
        logical, intent(out):: done, alldone
        logical, intent(in):: first
        real(kind=dp), intent(out), dimension(integrator_points, 3) :: pos
        real(kind=dp), intent(out), dimension(integrator_points, 3) :: vel
        integer(kind=i4b), intent(out) :: otype
        integer(kind=i4b), save :: dith = 0
        integer(kind=i4b)    :: temporbit
        real    (kind=dp ), dimension(6) :: YY ! (BT)
        real    (kind=dp )   :: Enjo,r_mean,lz2,vca,Svr,Svt,Svz ! (BT)
        real(kind=dp), dimension(5) :: moments
        real    (kind=dp ),dimension(3) :: moments2 ! (BT)
        !----------------------------------------------------------------------
        alldone = .false.
        if (first) then
            integrator_current = integrator_current + 1
            dith = 0
            totalnotregularizable = 0
            print *, "  * Starting integrating :", integrator_current
        end if
        dith = dith + 1
        if (dith <= integrator_dithering**3 .and. &
            integrator_current <= integrator_number) then

            call integrator_whichorbit(integrator_current, dith, temporbit)

            ! check if there is an unregurizable compenent in this set
            if (regurizable(temporbit) == 1) totalnotregularizable = 1
            !print*,"  * Starting integrating :",integrator_current,dith,temporbit
            call real_integrator(temporbit, pos, vel)
            call integrator_find_orbtype(otype, moments, moments2, pos, vel) ! (BT)
            integrator_orbittypes(dith) = otype
            integrator_moments(:, dith) = moments
            done = .false.

            ! write orbit information  in case of figure rotation (BT)
            if (Omega /= 0.0_dp ) then
               r_mean = moments(4)*conversion_factor**(-1.0_dp)
               lz2    = moments(3)*conversion_factor**(-1.0_dp)
               !cir = moments(3)/( SQRT(moments(5))*moments(4) )
               vca  = moments(5)
               Svr  = moments2(1)
               Svt  = moments2(2)
               Svz  = moments2(3)

               YY(1) =  pos(1,1)
               YY(2) =  pos(1,2)
               YY(3) =  pos(1,3)
               YY(4) =  vel(1,1)
               YY(5) =  vel(1,2)
               YY(6) =  vel(1,3)

               ! Compute and store start energy at begin point
               call computer_energy(YY,Enjo)

               write (unit=32, fmt="(25es13.5)")  Enjo, r_mean, lz2, vca, Svr, Svt, Svz
            endif
            ! end write orbit information

        else
            ! integrating done. Set 'done' to true and return.
            pos(:, :) = 0.0_dp
            vel(:, :) = 0.0_dp
            done = .true.
        end if

        if (integrator_current > integrator_number) alldone = .true.

    end subroutine integrator_integrate

!  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  subroutine integrator_integrate_nodith(pos,vel,otype,done,first,alldone)
!    logical           ,intent(out):: done,alldone
!    logical           ,intent(in):: first
!    real    (kind=dp ),intent(out),dimension(integrator_points,3) :: pos
!    real    (kind=dp ),intent(out),dimension(integrator_points,3) :: vel
!    integer (kind=i4b),intent(out) :: otype
!    integer (kind=i4b),save :: dith=0
!    integer (kind=i4b)    :: temporbit
!  !----------------------------------------------------------------------
!    if (first) then
!       integrator_current=integrator_current+1
!       dith=(integrator_dithering**3)/2+1
!       totalnotregularizable=0
!       print*,"  * Starting integrating :",integrator_current
!    endif
!    if (first .and. &
!         integrator_current <= integrator_number) then
!
!       call integrator_whichorbit(integrator_current,dith,temporbit)
!
!       ! check if there is an unregurizable compenent in this set
!       if (regurizable(temporbit) == 1) totalnotregularizable=1
!       print*,"  * Starting integrating :",integrator_current,dith,temporbit
!       call real_integrator (temporbit,pos,vel)
!       call integrator_find_orbtype(otype,pos,vel)
!       integrator_orbittypes(dith)=otype
!       done=.false.
!    else
!       ! integrating done. Set 'done' to true and return.
!       pos(:,:)=0.0_dp
!       vel(:,:)=0.0_dp
!       done=.true.
!    end if
!
!    if (integrator_current > integrator_number) alldone=.true.
!
!  end subroutine integrator_integrate_nodith

! This routines finds the orbit dither number
    subroutine integrator_whichorbit(orbit, dith, dithorbit)
        integer(kind=i4b), intent(in)      :: orbit, dith
        integer(kind=i4b), intent(out)     :: dithorbit
        integer(kind=i4b) :: E1, I2, I3, nd
        integer(kind=i4b) :: d1, d2, d3, DO1, DO2, DO3

        nd = integrator_dithering
        ! Reconstruct E, I2 and I3 number from orbit number:
        I3 = modulo((orbit - 1), nI3/nd) + 1
        I2 = modulo(((orbit - 1)/(nI3/nd)), nI2/nd) + 1
        E1 = ((orbit - 1)/(nI3*nI2/(nd*nd))) + 1
        ! Find which orbit dither we are currently at:
        d3 = modulo((dith - 1), nd) + 1
        d2 = modulo(((dith - 1)/nd), nd) + 1
        d1 = ((dith - 1)/(nd*nd)) + 1
        ! combine previous results to find the undithered orbit number
        DO3 = (I3 - 1)*nd + d3
        DO2 = (I2 - 1)*nd + d2
        DO1 = (E1 - 1)*nd + d1
        dithorbit = DO3 + ((DO2 - 1)*nI3) + ((do1 - 1)*nI3*nI2)
    end subroutine integrator_whichorbit

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine real_integrator(orbit, pos, vel)
        integer(kind=i4b), intent(in)                                :: orbit
        real(kind=dp), intent(out), dimension(integrator_points, 3) :: pos
        real(kind=dp), intent(out), dimension(integrator_points, 3) :: vel
        !----------------------------------------------------------------------
        integer(kind=i4b) :: IOUT, IDID, ITOL
        integer(kind=i4b), parameter :: N = 6, NRDENS = 6, LWORK = 11*N + 8*NRDENS + 21
        integer(kind=i4b), parameter :: LIWORK = NRDENS + 21
        real(kind=dp) :: X, Xend, RTOL, ATOL
        real(kind=dp), dimension(6) :: Y
        real(kind=dp), dimension(lwork) :: WORK
        integer(kind=i4b), dimension(liwork) :: IWORK
        real(kind=dp), dimension(2) :: RPAR
        integer(kind=i4b), dimension(1) :: IPAR
        real(kind=dp), save :: stepsize = 0.0_dp
        real(kind=dp) :: Ebeg, Eend
        integer(kind=i4b), save :: stored_orbit = 0

        if (.NOT. allocated(pos_old)) then
            ! setting up the array for storing the points.
            allocate (pos_old(integrator_points, 3), vel_old(integrator_points, 3))
            allocate (pos_t(integrator_points, 3), vel_t(integrator_points, 3))

        else
            if (integrator_points /= size(pos_old, 1)) then
                deallocate (pos_old, vel_old, pos_t, vel_t)
                allocate (pos_old(integrator_points, 3), vel_old(integrator_points, 3))
                allocate (pos_t(integrator_points, 3), vel_t(integrator_points, 3))
                ! Reset the flag to indicate that we do not have an orbit stored.
                stored_orbit = 0
            end if
        end if

        ! --- REQUIRED TOLERANCE
        ! The integrator keeps the local error on Y(I)
        ! below RTOL(I)*ABS(Y(I)) + ATOL(I).
        ITOL = 0 ! Tolerances are scalars
        RTOL = integrator_accuracy
        ATOL = integrator_accuracy!/10000.0_dp
        ! The Absolutie Tolerance used to be set-up as a vector
        ! ATOL= 1e-6* (/ rcirc(orbit), rcirc(orbit), rcirc(orbit), &
        !                   vcirc(orbit), vcirc(orbit), vcirc(orbit) /)

        Y(1) = xini(orbit)  !  x(0)
        Y(2) = yini(orbit)  !  y(0)
        Y(3) = zini(orbit)  !  z(0)
        Y(4) = vxini(orbit)  ! vx(0)
        Y(5) = vyini(orbit)  ! vy(0)
        Y(6) = vzini(orbit)  ! vz(0)

        ! Compute and store start energy at begin point
        call computer_energy(Y, Ebeg)

        DO  ! While integrating has not succeeded
            ! setting up the array for storing the points.
            ! --- DIMENSION OF THE SYSTEM
            IDID = 0
            ! --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
            IOUT = 2
            ! --- INITIAL VALUES
            X = 0.0

            Y(1) = xini(orbit)  !  x(0)
            Y(2) = yini(orbit)  !  y(0)
            Y(3) = zini(orbit)  !  z(0)
            Y(4) = vxini(orbit)  ! vx(0)
            Y(5) = vyini(orbit)  ! vy(0)
            Y(6) = vzini(orbit)  ! vz(0)

            ! --- ENDPOINT OF INTEGRATION
            XEND = integrator_n_orbits*tcirc(orbit)
            ! Stepsize. Make sure enough steps are in the integration by adding
            ! room for a couple of  extra steps. RvdB 19/12/04
            RPAR(2) = XEND/(integrator_points + 4)
            RPAR(1) = 0.0 ! Unused

            ! --- VALUES FOR PARAMETERS
            ! Default values are used when IWORK or WORK are zero
            IWORK(:) = 0
            WORK(:) = 0.0_dp
            ! Maximum of allowed steps (just really high) (100000)
            IWORK(1) = max(floor(min(dble(huge(1_i4b) - 1), 1.0_dp/RTOL*integrator_n_orbits/200), kind=i4b), 100000_i4b)
            !number of dense components needed. (all in our case)
            IWORK(5) = NRDENS
            !stiffness detection (negative --> do not try to detect)
            IWORK(4) = -1
            !give guess of stepsize
            WORK(7) = stepsize

            !CALL OF THE SUBROUTINE DOPRI8 ( The dop853 integrator.)
            CALL DOP853(N, derivs, X, Y, XEND, RTOL, ATOL, ITOL, SOLOUT, IOUT, &
                        WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)

!   IPAR(1)    COUNT   NUMBER OF STORED INTEGRATION STEPS IN pos_t AND vel_t
!   IWORK(17)  NFCN    NUMBER OF FUNCTION EVALUATIONS
!   IWORK(18)  NSTEP   NUMBER OF COMPUTED STEPS
!   IWORK(19)  NACCPT  NUMBER OF ACCEPTED STEPS
!   IWORK(20)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST),
!                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED)
!   print*,"  * Integrator output"
!   print*,iwork(17:20)

            ! compute final energy
            call computer_energy(Y, Eend)

            select case (IDID)
            case (-1)
                stop "integrator:  INPUT IS NOT CONSISTENT,"
            case (-2)
                print *, "integrator:  LARGER NMAX IS NEEDED,"
            case (-3)
                print *, "integrator:  STEP SIZE BECOMES TOO SMALL."
            case (-4)
                stop "integrator:  PROBLEM IS PROBABLY STIFF (INTERRUPTED)."
            case default
                ! Integrating went ok!

                if (abs((Ebeg - Eend)/Ebeg) .lt. 0.01) then
                    ! If Integration conserved energy within 1 percent then do
                    !   - Store pos and vel as failsafe for next orbit integration
                    !   - Swap the points to the correct array.
                    pos(:, :) = pos_t(:, :)
                    pos_old(:, :) = pos_t(:, :)
                    vel_old(:, :) = vel_t(:, :)
                    vel(:, :) = vel_t(:, :)
                    stored_orbit = 1
                    stepsize = work(7)  ! Store integrater stepsize for next itegration
                    EXIT  ! integration was succesfull
                end if
                print *, "Energy conserved to ", (Ebeg - Eend)/Ebeg*100.0_dp, ", Increasing integrator accuracy"
            end select

            if (RTOL .lt. 1e-12_dp) then
                print *, '  * orbit ', orbit, ' failed. Energy conserved up to ', (Ebeg - Eend)/Ebeg*100.0_dp
                ! Orbit integration unsuccesfull even at higher accuracy
                if (stored_orbit .eq. 0) stop 'Abort, No backup orbit stored'
                pos_t(:, :) = pos_old(:, :)
                vel_t(:, :) = vel_old(:, :)
                EXIT ! Integration not succesfull
            end if

            ! Increase accuracy and try again.
            RTOL = RTOL*0.1_dp
            ATOL = ATOL*0.1_dp

            print *, '  * Retrying orbit', orbit
        end do

    end subroutine real_integrator

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine integrator_setup_write(handle)
        integer(kind=i4b), intent(in) :: handle
        integer(kind=i4b)            :: t1, t2, t3, t4
        !----------------------------------------------------------------------
        print *, "  * Writing integrator output header"
        t1 = (nEner*nI2*nI3/integrator_dithering**3)
        t2 = nEner/integrator_dithering
        t3 = nI2/integrator_dithering
        t4 = nI3/integrator_dithering
        write (unit=handle) t1, t2, t3, t4, integrator_dithering

        write (unit=30, fmt=*) Nener*ni2*ni3, integrator_dithering**3

    end subroutine integrator_setup_write

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine integrator_write(hdl)
        integer(kind=i4b), intent(in) :: hdl
        !----------------------------------------------------------------------
        integer(kind=i4b)            :: nd, orbit, I3, I2, E1
        orbit = integrator_current

        nd = integrator_dithering
        ! Reconstruct E, I2 and I3 number from orbit number:
        I3 = modulo((orbit - 1), nI3/nd) + 1
        I2 = modulo(((orbit - 1)/(nI3/nd)), nI2/nd) + 1
        E1 = ((orbit - 1)/(nI3*nI2/(nd*nd))) + 1

        ! write information about the orbit
        write (unit=hdl) orbit, E1, I2, I3, totalnotregularizable
        write (unit=hdl) integrator_orbittypes(:)

        write (unit=30, fmt="(25es13.5)") integrator_moments(:, :)
    end subroutine integrator_write

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine ini_integ()
        !----------------------------------------------------------------------
        character(len=30) :: infil
        integer(kind=i4b) :: i

        print *, "  * Give name of file with begin conditions"
        read (unit=*, fmt="(a30)") infil
        print *, "   ", infil
        open (unit=31, file=infil, status="OLD", action="read", position="rewind")

        read (unit=31, fmt=*) nEner, nI2, nI3

        i = nEner*nI2*nI3

        print *, " * Orbits in the input file"
        print *, nEner, nI2, nI3, i

        IF (.NOT. ALLOCATED(xini)) &
            allocate (xini(i), yini(i), zini(i), vxini(i), vyini(i), vzini(i), rcirc(i), &
                      tcirc(i), vcirc(i), gEner(i), gI2(i), gI3(i), regurizable(i))

        read (unit=31, fmt="(3I5,9ES30.10,I4)") (gener(i), gi2(i), gi3(i), xini(i), &
                                                 yini(i), zini(i), vxini(i), vyini(i), &
                                                 vzini(i), rcirc(i), tcirc(i), &
                                                 vcirc(i), regurizable(i), &
                                                 i=1, nEner*nI2*nI3)

        close (unit=31)

    end subroutine ini_integ

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine derivs(N, xin, yin, dydx, RPAR, IPAR)
      use initial_parameters, only: Omega        ! (BT)
        use interpolpot, only: ip_accel
        !use triaxpotent, only : tp_accel
        integer, intent(in)              :: N
        real(kind=dp), intent(in)              :: xin
        real(kind=dp), intent(in), dimension(6) :: yin
        real(kind=dp), intent(out), dimension(6) :: dydx
        real(kind=dp), intent(inout), dimension(2) :: RPAR
        integer(kind=i4b), intent(inout), dimension(1) :: IPAR
        !----------------------------------------------------------------------

        ! subroutine which returns the right-hand side derivatives.
        !   x      = t
        !   yin(1) = x
        !   yin(2) = y
        !   yin(3) = z
        !   yin(4) = dx/dt
        !   yin(5) = dy/dt
        !   yin(6) = dz/dt

        if (Omega == 0.0_dp ) then
           ! First calculate the true accelerations at the given position
           dydx(1:3) = yin(4:6)
           call ip_accel(yin(1), yin(2), yin(3), dydx(4), dydx(5), dydx(6))
        else
           ! Rotating frame with Omega pattern speed (BT)
           dydx(1) =  yin(4) + Omega*yin(2)  ! extra terms for integration in co-rotating frame
           dydx(2) =  yin(5) - Omega*yin(1)
           dydx(3) =  yin(6)

           call ip_accel(yin(1),yin(2),yin(3),dydx(4),dydx(5),dydx(6))

           dydx(4) =  dydx(4) + Omega*yin(5)
           dydx(5) =  dydx(5) - Omega*yin(4)
           dydx(6) =  dydx(6)
        end if

    end subroutine derivs

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE SOLOUT(NR, XOLD, X, Y, N, CON, ICOMP, ND, RPAR, IPAR, IRTRN)
        ! --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
        ! --- BY USING "CONTD8", THE CONTINUOUS COLLOCATION SOLUTION
        integer, intent(in)    :: NR, N, ND
        real(kind=dp), intent(in)    :: XOLD, X
        integer, intent(inout) :: IRTRN
        real(kind=dp), intent(in), dimension(N)    :: Y
        real(kind=dp), intent(in), dimension(8*ND) :: CON
        integer(kind=i4b), intent(in), dimension(ND)   :: ICOMP
        real(kind=dp), intent(inout), dimension(2)    :: RPAR
        integer(kind=i4b), intent(inout), dimension(1)    :: IPAR
        !----------------------------------------------------------------------
        real(kind=dp), save          :: XOUT
        integer(kind=i4b), save          :: count = 0
        real(kind=dp)               :: contd8, step, rnd
        real(kind=dp)               :: ran1
        step = RPAR(2)
        IF (NR == 1) THEN
            ! Start storing the orbit after 1+? steps to avoid aliasing
            rnd = ran1(1)
            ! call random_number(rnd) ! 0 < rnd < 1
            XOUT = X + step*(1.0_dp + rnd)
            count = 0
            IPAR(1) = 0
        ELSE
            do
                ! Make sure count > integrator_points to make sure we do not
                ! go out of bounds on pos_t(:,:). RvdB, DK 16/06/03
                IF (X < XOUT .or. count >= integrator_points) exit
                count = count + 1
                IPAR(1) = count
                pos_t(count, 1) = CONTD8(1, XOUT, CON, ICOMP, ND)
                pos_t(count, 2) = CONTD8(2, XOUT, CON, ICOMP, ND)
                pos_t(count, 3) = CONTD8(3, XOUT, CON, ICOMP, ND)
                vel_t(count, 1) = CONTD8(4, XOUT, CON, ICOMP, ND)
                vel_t(count, 2) = CONTD8(5, XOUT, CON, ICOMP, ND)
                vel_t(count, 3) = CONTD8(6, XOUT, CON, ICOMP, ND)
                XOUT = XOUT + step
            end do
        END IF
    END SUBROUTINE SOLOUT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine integrator_find_orbtype(type, moments, moments2, pos, vel) ! (BT)
        integer(kind=i4b), intent(out) :: type
        real(kind=dp), intent(out), dimension(:) :: moments, moments2 ! (BT)
        real(kind=dp), intent(in), dimension(:, :) :: pos
        real(kind=dp), intent(in), dimension(size(pos, 1), 3) :: vel
        real (kind=dp),dimension(size(pos,1)) :: vr, vt, vz ! (BT)
        real (kind=dp) :: sd_vr, sd_vt, sd_vz, mean_vr, mean_vt, mean_vz, n ! (BT)
        !----------------------------------------------------------------------
        real(kind=dp) :: lxc, lyc, lzc
        real(kind=dp), dimension(size(pos, 1)) :: t
        real(kind=dp), parameter :: nul = 0.0_dp
        !  Lx = y*Vz-z*Vy
        t = pos(:, 2)*vel(:, 3) - pos(:, 3)*vel(:, 2)
        lxc = maxval(t)*minval(t)
        moments(1) = sum(t)/size(pos, 1)

        !  Ly = z*Vx-x*Vz
        t = pos(:, 3)*vel(:, 1) - pos(:, 1)*vel(:, 3)
        lyc = maxval(t)*minval(t)
        moments(2) = sum(t)/size(pos, 1)

        !  Lz = x*Vy-y*Vx
        t = pos(:, 1)*vel(:, 2) - pos(:, 2)*vel(:, 1)
        lzc = maxval(t)*minval(t)
        moments(3) = sum(t)/size(pos, 1)

        ! assume orbit is chaotic, unless proven otherwise.
        type = 5

        ! X tube
        if (lxc > nul .and. lyc < nul .and. lzc < nul) type = 1
        ! Y tube
        if (lxc < nul .and. lyc > nul .and. lzc < nul) type = 2
        ! Z tube
        if (lxc < nul .and. lyc < nul .and. lzc > nul) type = 3
        ! Box
        if (lxc < nul .and. lyc < nul .and. lzc < nul) type = 4

        ! vector lenghts should be equal, otherwise DIM  is wrong
        !print*,size(sum(pos(:,:)**2,dim=2)) ,  size(pos,1)

        ! mean radius
        moments(4) = sum(sqrt(sum(pos(:, :)**2, dim=2)))/size(pos, 1)

        ! second moment = vxx+vyy+vzz + 2vxy + 2vyz +2vzx
        moments(5) = sum(sum(vel(:, :)**2, dim=2) &
                         + 2*(vel(:, 1)*vel(:, 2) &
                              + vel(:, 2)*vel(:, 3) &
                              + vel(:, 3)*vel(:, 1)))/size(pos, 1)

        !print*,moments(3),moments(5),moments(3)/moments(4)/sqrt(moments(5))

        ! convert velocity to cylinderical  (anisotropy) (BT)
        vr = ( pos(:,1) * vel(:,1) +  pos(:,2) * vel(:,2) ) / sqrt(pos(:,1)**2+pos(:,2)**2)
        vt = ( pos(:,1) * vel(:,2) +  pos(:,2) * vel(:,1) ) / sqrt(pos(:,1)**2+pos(:,2)**2)
        vz = vel(:,3)

        N = size(pos,1)
        mean_vr = sum(vr) / N
        sd_vr = sqrt(sum((vr-mean_vr)**2) / N)

        mean_vt = sum(vt) / N
        sd_vt = sqrt(sum((vt-mean_vt)**2) / N)

        mean_vz = sum(vz) / N
        sd_vz = sqrt(sum((vz-mean_vz)**2) / N)

        moments2(1) =  sd_vr
        moments2(2) =  sd_vt
        moments2(3) =  sd_vz

    end subroutine integrator_find_orbtype

    subroutine computer_energy(Y, E)
      use initial_parameters, only: Omega ! (BT)
        ! Compute the energy of a particle.
        use interpolpot, only: ip_potent
        real(kind=dp), intent(in), dimension(6) :: Y
        real(kind=dp), intent(out) :: E
        real(kind=dp) :: ep
        call ip_potent(y(1), y(2), y(3), Ep)
        if (Omega /= 0.0_dp ) then
           ep = ep + Omega*(y(1)*y(5)-y(2)*y(4)) ! check the conservation of effective potential (BT)
        end if
        E = ep - 0.5_dp*(y(4)**2 + y(5)**2 + y(6)**2)
    end subroutine computer_energy

end module integrator

!######################################################################
!######################################################################
!######################################################################

! $Id: orblib_f.f90,v 1.3 2011/10/25 08:48:45 bosch Exp $

! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! july 2002 Sterrewacht Leiden.

module projection
    ! module doing the circular projection of the orbit.
    use numeric_kinds
    implicit none
    private

    ! Number of projections around axis
    integer(kind=i4b), private :: proj_number

    ! Readin the number of rotation to be done.
    public :: projection_setup, projection_change_direction

    ! stop projection module
    public :: projection_stop

    ! Project pos(3,:),vel(3,:) to proj(2,:),losvel(:),velx(:),vely(:)
    ! Using the n'th projection
    private :: project_n

    ! Project pos(3,:),vel(3,:) to proj(2,:),losvel(:),velx(:),vely(:)
    !Done is set to .true. if all projections are finished
    public :: project

    ! amount of symmetry multipication. ( for triaxial galaxies there are
    ! 8 symmetries per orbit, but we do one (1) at a time.)
    integer(kind=i4b), public, parameter :: projection_symmetry = 1

    real(kind=dp), public :: theta_proj, phi_proj, psi_proj

contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine projection_setup()
        use initial_parameters, only: theta_view, phi_view, psi_view
        !----------------------------------------------------------------------
        print *, "  ** Projection setup"
        print *, "  * 8 Projections for Triaxial model"
        proj_number = 8
        print *, "  * Inclination of the model is (theta,phi): ", theta_view, phi_view, psi_view
        theta_proj = theta_view
        phi_proj = phi_view
        psi_proj = psi_view
        print *, theta_proj, phi_proj, psi_proj
        print *, "  ** Projection setup finished"

    end subroutine projection_setup

    subroutine projection_change_direction()
        !----------------------------------------------------------------------
        real(kind=dp)              :: t1, t2, t3
        print *, "  ** Projection change direction"
        print *, "  * deprojection angles are (theta,phi,psi): ", &
            theta_proj/(pi_d/180.0_dp), phi_proj/(pi_d/180.0_dp), &
            psi_proj/(pi_d/180.0_dp)
        print *, "  * Give new angles (theta, phi, psi):"
        print *, "  * Anwser -501 0 0 to keep current values"
        read *, t1, t2, t3
        print *, t1, t2, t3
        if (t1 > -500) then
            theta_proj = t1*(pi_d/180.0_dp)
            phi_proj = t2*(pi_d/180.0_dp)
            psi_proj = t3*(pi_d/180.0_dp)
        end if
        print *, "  * New deprojection angles are (theta,phi,psi): ", &
            theta_proj/(pi_d/180.0_dp), phi_proj/(pi_d/180.0_dp), &
            psi_proj/(pi_d/180.0_dp)

        print *, theta_proj, phi_proj, psi_proj

    end subroutine projection_change_direction

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine projection_stop()
        !----------------------------------------------------------------------
        ! empty function

    end subroutine projection_stop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine project_n(type, pos, vel, proj, losvel, vel2d, n)
      use initial_parameters, only :  Omega ! (BT)
        ! use initial_parameters, only : theta_view, phi_view
        ! pos( :, (r,z) )
        real(kind=dp), intent(in), dimension(:, :)             :: pos
        ! vel (:, (r,z,theta))
        real(kind=dp), intent(in), dimension(size(pos, 1), 3)  :: vel
        ! proj(:,(x',y'))
        real(kind=dp), intent(out), dimension(size(pos, 1), 2) :: proj
        ! losvd (:)
        real(kind=dp), intent(out), dimension(size(pos, 1))    :: losvel
        ! velx(:), vely(:)
        real(kind=dp), intent(out), dimension(size(pos, 1), 2) :: vel2d
        integer(kind=i4b), intent(in)                          :: type, n
        !----------------------------------------------------------------------
        real(kind=dp)                          :: t1, t2, t3, theta, phi

        real (kind=dp),dimension(3,8,5) ::vsgn          ! (BT)
        real (kind=dp),dimension(3,8) ::psgn            ! (BT)

        ! Signs of the (vx,vy,vz) for each Projection and type of Orbit
        real(kind=dp), dimension(3, 8, 5), &
            parameter :: vsgn1 = reshape((/ &
                                        ! X tubes
                                        1, 1, 1, -1, 1, 1, 1, 1, -1, -1, 1, -1, &
                                        -1, -1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, &
                                        ! Y tubes
                                        1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, &
                                        -1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, &
                                        ! Z tubes
                                        1, 1, 1, 1, -1, -1, -1, -1, 1, -1, 1, -1, &
                                        1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, 1, &
                                        ! Boxed
                                        1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, &
                                        1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, &
                                        ! Stochastic
                                        1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, &
                                        1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1/), (/3, 8, 5/))

        !Signs of the x,y,z for each projection  :psgn( [x,y,z], project )
        real(kind=dp), dimension(3, 8), &
            parameter :: psgn1 = reshape((/ &
                                        1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, &
                                        1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1/), (/3, 8/))


        ! Signs of the (vx,vy,vz) for each Projection and type of Orbit
        ! (BT) 4-fold symmetry same as before
        real (kind=dp),dimension(3,8,5),parameter :: vsgn2= reshape((/  &
             ! X tubes
             1 , 1 , 1    ,1 , 1 , 1  ,  -1 , -1 ,1  ,   -1 , -1 , 1 , &
             1 ,1 , -1    ,1 ,1 , -1  , -1 ,-1 ,-1  ,  -1 ,-1 , -1 , &
             ! Y tubes
             1 , 1 , 1    , 1 , 1 ,1  ,  1 , 1 ,-1  ,  1 , 1 , -1 , &
             -1 , -1 , 1   ,-1 , -1 ,1  , -1 ,-1 ,-1  , -1 ,-1 , -1 , &
             ! Z tubes
             1 , 1 , 1    , 1 ,1 , 1  , -1 ,-1 , 1  , -1 , -1 , 1 , &
             1 , 1 ,-1    , 1 ,1 ,-1  , -1 ,-1 ,-1  , -1 , -1 ,-1 , &
             ! Boxed
             1 , 1 , 1    ,1 , 1 , 1  , -1 ,-1 , 1  ,  -1 ,-1 , 1 , &
             1 , 1 ,-1    ,1 , 1 ,-1  , -1 ,-1 ,-1  ,  -1 ,-1 ,-1 , &
             ! Stochastic
             1 , 1 , 1    ,1 , 1 , 1  , -1 ,-1 , 1  ,  -1 ,-1 , 1 , &
             1 , 1 ,-1    ,1 , 1 ,-1  , -1 ,-1 ,-1  ,  -1 ,-1 ,-1 /),(/3,8,5/))
        !Signs of the x,y,z for each projection  :psgn( [x,y,z], project )
        real (kind=dp),dimension(3,8),parameter :: psgn2= reshape((/  &
             1 , 1 , 1   , 1 , 1 , 1   , -1 , -1 , 1 ,  -1 , -1 , 1 , &
             1 , 1 ,-1   , 1 , 1 ,-1  , -1 , -1 ,-1 ,  -1 , -1 ,-1 /),(/3,8/))

        ! Use 8-fold for non-rotating, but 4-fold for rotating (BT)
        vsgn=vsgn1
        psgn=psgn1
        if (Omega /= 0.0_dp ) then
           vsgn=vsgn2
           psgn=psgn2
        endif

        theta = theta_proj
        phi = phi_proj

        ! check orbit type
        if (type > 5 .or. type < 1) stop "project_n: Wrong orbit type"

        ! Use sign matrix for the symmetries.
        ! Using the inverse (transpose) of the projection (eq. 4) of Thesis Ellen.

        ! x'
        t1 = -sin(phi)*psgn(1, n)
        t2 = cos(phi)*psgn(2, n)
        proj(:, 1) = t1*pos(:, 1) + t2*pos(:, 2)

        ! y'
        t1 = -cos(theta)*cos(phi)*psgn(1, n)
        t2 = -cos(theta)*sin(phi)*psgn(2, n)
        t3 = sin(theta)*psgn(3, n)
        proj(:, 2) = t1*pos(:, 1) + t2*pos(:, 2) + t3*pos(:, 3)

        ! v_LOS
        t1 = sin(theta)*cos(phi)*vsgn(1, n, type)
        t2 = sin(theta)*sin(phi)*vsgn(2, n, type)
        t3 = cos(theta)*vsgn(3, n, type)
        losvel(:) = t1*vel(:, 1) + t2*vel(:, 2) + t3*vel(:, 3)

        ! pm_x
        t1 = -sin(phi)*vsgn(1, n, type)
        t2 = cos(phi)*vsgn(2, n, type)
        vel2d(:, 1) = t1*vel(:, 1) + t2*vel(:, 2)

        ! pm_y
        t1 = -cos(theta)*cos(phi)*vsgn(1, n, type)
        t2 = -cos(theta)*sin(phi)*vsgn(2, n, type)
        t3 = sin(theta)*vsgn(3, n, type)
        vel2d(:, 2) = t1*vel(:, 1) + t2*vel(:, 2) + t3*vel(:, 3)

        ! ! v_radial and v_tangential START
        ! ! r'
        ! projr = sqrt(proj(:, 1)**2 + proj(:, 2)**2)

        ! ! pm_radial
        ! velr = ( proj(:, 1)*velx(:) + proj(:, 2)*vely(:) ) / projr

        ! ! pm_tangential
        ! velt = ( proj(:, 1)*vely(:) + proj(:, 2)*velx(:) ) / projr
        ! ! v_radial and v_tangential END

        !xaa = (-sin(phi)*x+cos(phi)*y)*sin(psi)-(-cos(theta)*cos(phi)*x-cos(theta)*sin(phi)*y+sin(theta)*z)*cos(psi);
        !yaa = (-sin(phi)*x+cos(phi)*y)*cos(psi)+(-cos(theta)*cos(phi)*x-cos(theta)*sin(phi)*y+sin(theta)*z)*sin(psi);

        !  t1 = sin(phi)
        !  t3 = cos(phi)
        !  t5 = -t1*x+t3*y
        !  t6 = sin(psi)
        !  t8 = cos(theta)
        !  t13 = sin(theta)
        !  t15 = -t8*t3*x-t8*t1*y+t13*z
        !  t16 = cos(psi)
        !  v(1) = t5*t6-t15*t16
        !  v(2) = t5*t16+t15*t6

        !  v(1) = x*sin(psi)-y*cos(psi)
        !  v(2) = x*cos(psi)+y*sin(psi)

    end subroutine project_n

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine project(type, pos, vel, proj, losvel, vel2d, done, first)
        integer(kind=i4b), intent(in)                        :: type
        real(kind=dp), intent(in), dimension(:, :) :: pos
        real(kind=dp), intent(in), dimension(size(pos, 1), 3) :: vel
        real(kind=dp), intent(out), dimension(size(pos, 1)*projection_symmetry, 2) &
             & :: proj
        real(kind=dp), intent(out), dimension(size(pos, 1)*projection_symmetry) &
             & :: losvel
        real(kind=dp), intent(out), dimension(size(pos, 1)*projection_symmetry, 2) &
             & :: vel2d
        logical, intent(out)                          :: done
        logical, intent(in)                          :: first
        !----------------------------------------------------------------------
        integer(kind=i4b), save :: count = 0

        ! reset counter if this is the first projection for this orbit
        if (first) count = 0
        count = count + 1
        done = .false.

        if (count <= proj_number) then
            call project_n(type, pos, vel, proj, losvel, vel2d, count)
        else
            count = proj_number + 1
            done = .true.
            proj(:, :) = 0.0_dp
            losvel(:) = 0.0_dp
            vel2d(:, :) = 0.0_dp
        end if

    end subroutine project

end module projection

!#########################################################
!#########################################################
!#########################################################

! $Id: orblib_f.f90,v 1.3 2011/10/25 08:48:45 bosch Exp $

! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! Sterrewacht Leiden, The Netherlands

module psf
    use numeric_kinds
    implicit none
    private

    ! * Module for the PSF generation.
    ! The way PSF are done in this program is quite simple. We just
    ! take the original point and modify it with a configurable
    ! random gaussian offset (psf_size).

    ! how many psf are there?
    integer(kind=i4b), public                           :: psf_n
    ! what is the histogram dimension of the psf(psf_n)?
    integer(kind=i4b), public, allocatable, dimension(:) :: psf_hist_dim
    ! kind of psf(psf_n)
    integer(kind=i4b), private, allocatable, dimension(:) :: psf_kind
    ! size of psf for (n,psf)
    real(kind=dp), private, allocatable, dimension(:, :) :: psf_sigma
    ! intensity of the psf(n,psf)
    real(kind=dp), private, allocatable, dimension(:, :) :: psf_iten
    ! (i,j,pf) contains a sigma's in random order for psf pf
    real(kind=dp), private, allocatable, dimension(:, :) :: psf_randomsigma
    ! setupup of psf variables
    public :: psf_setup
    ! generate gaussian psf points of input array.
    public :: psf_gaussian

    public :: psf_stop

    ! Find the sigma of the psf number #.
    public :: psf_cal_sigma

    ! Generates an array with proportionals sigma's of a MGE-PSF
    private:: psf_sigma_map

contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine psf_stop()
        !----------------------------------------------------------------------
        if (allocated(psf_kind)) then
            deallocate (psf_kind, psf_sigma, psf_iten, psf_hist_dim)
            deallocate (psf_randomsigma)
        end if

    end subroutine psf_stop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine psf_setup()
        use initial_parameters, only: conversion_factor
        use random_gauss_generator, only: random_gauss_seed
        !----------------------------------------------------------------------
        integer(kind=i4b) :: i, j

        print *, "  ** Setting up PSF module"
        print *, "  * How many different psf's?"
        read *, psf_n
        print *, "   ", psf_n
        allocate (psf_kind(psf_n), psf_hist_dim(psf_n))

        do i = 1, psf_n
            print *, "  * How many gaussians does the ", i, "psf consist of?"
            read *, psf_kind(i)
            print *, psf_kind(i)
            if (psf_kind(i) < 1) stop "gaussian value too low"
        end do

        allocate (psf_sigma(maxval(psf_kind(:)), psf_n))
        allocate (psf_iten(maxval(psf_kind(:)), psf_n))

        do i = 1, psf_n
            print *, "   Intensity, sigma of the gauss for PSF ", i
            do j = 1, psf_kind(i)
                read *, psf_iten(j, i), psf_sigma(j, i)
                print *, psf_iten(j, i), psf_sigma(j, i)
            end do
        end do

        ! convert sizes arcsec to km
        psf_sigma(:, :) = psf_sigma(:, :)*conversion_factor
        call random_gauss_seed()
        call psf_sigma_map()

        print *, "  ** PSF module setup Finished"

    end subroutine psf_setup

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine psf_gaussian(pf, vec, gaus)
        use random_gauss_generator
        integer(kind=i4b), intent(in)               :: pf
        ! input vectors (n,2)
        real(kind=dp), dimension(:, :), intent(in) :: vec
        ! output vectors (n,2)
        real(kind=dp), dimension(:, :), intent(out)::gaus
        !----------------------------------------------------------------------
        real(kind=sp), dimension(size(vec, 1)) :: t
        integer(kind=i4b), dimension(size(vec, 1)) :: ind
        integer(kind=i4b) :: j
        real(kind=dp) :: ran1

        if (psf_kind(pf) == 1) then
            ! One gaussian in this psf
            if (psf_sigma(1, pf) > 1.0_dp) then
                call random_gauss(gaus(:, :))
                gaus(:, :) = vec(:, :) + gaus(:, :)*psf_sigma(1, pf)
            else
                ! psf size is tiny, so no convolution is done.
                gaus(:, :) = vec(:, :)
            end if
        else
            ! MGE PSF. Use the randomsigma to convolve these points
            ! Each sigma has a chance of being used proprotional
            ! to the weight of the corresponding Gaussian component.
            ! M. Cappellari, 14 January 2003
            call random_gauss(gaus(:, :))
            do j = 1, size(vec, 1) ! no forall, want this to be serialized...
                t(j) = ran1(1)
            end do
            ! call random_number(t(:))
            ind = t*(size(vec, 1) - 1) + 1 ! n=size(vec,1) random integers in [1,n]
            forall (j=1:2)
                gaus(:, j) = vec(:, j) + gaus(:, j)*psf_randomsigma(ind(:), pf)
            end forall

        end if

    end subroutine psf_gaussian

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine psf_sigma_map()
        use integrator, only: integrator_points
        use projection, only: projection_symmetry
        !----------------------------------------------------------------------

        integer(kind=i4b)                          :: pf
        ! input vectors (n,2)
        integer(kind=i4b)                          :: i, j, sizex
        real(kind=dp), dimension(:), allocatable :: weightfl
        integer(kind=i4b), dimension(:), allocatable :: weightint
        ! random sigma's

        print *, "  * Making vectors filled of sigmas for psf convolution."
        sizex = integrator_points*projection_symmetry
        allocate (psf_randomsigma(sizex, psf_n))
        psf_randomsigma(:, :) = 0.0_dp

        do pf = 1, psf_n
            allocate (weightfl(psf_kind(pf)), weightint(psf_kind(pf) + 1))

            ! The weight of each PSF gaussian
            do i = 1, psf_kind(pf)
                weightfl(i) = abs(psf_iten(i, pf))
            end do
            do i = 1, psf_kind(pf)
                ! normalized cumulative sum
                weightint(i + 1) = nint(sum(weightfl(1:i))*((sizex - 1)/sum(weightfl(:)))) + 1
            end do
            ! range [1,sizex]
            weightint(1) = 1_i4b
            weightint(psf_kind(pf) + 1) = sizex

            ! Now we generate an array with the sigmas. Each sigma occurs a
            ! relative weighted amount of times in the array.
            do i = 1, psf_kind(pf)
                do j = weightint(i), weightint(i + 1)
                    psf_randomsigma(j, pf) = psf_sigma(i, pf)
                end do
            end do
            print *, 'Weight divided for psf', pf, ':'
            print *, weightint(:)
            deallocate (weightfl, weightint)
        end do
    end subroutine psf_sigma_map

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine psf_cal_sigma(pf, sigma)
        integer(kind=i4b), intent(in) :: pf
        real(kind=dp), intent(out):: sigma
        !----------------------------------------------------------------------
        sigma = maxval(psf_sigma(:, pf))

    end subroutine psf_cal_sigma

end module psf

!######################################################################
!######################################################################
!######################################################################

! $Id: orblib_f.f90,v 1.3 2011/10/25 08:48:45 bosch Exp $

! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! July 2002 Sterrewacht Leiden.

module aperture
    use numeric_kinds
    implicit none
    private

    ! Total number of apertures
    integer(kind=i4b), public                          :: aperture_n
    ! Number of apertures with 0d histograms (mass only)
    integer(kind=i4b), public                          :: ap_hist0d_n
    ! Number of apertures with 2d histograms
    integer(kind=i4b), public                          :: ap_hist2d_n
    ! histogram dimension for each aperture (0=0D, 1=1D, 2=2D)
    integer(kind=i4b), public, allocatable, dimension(:) :: ap_hist_dim

    ! number of bins in aperture
    integer(kind=i4b), public, allocatable, dimension(:) :: aperture_size
    ! Starting point of the aperture in flat array.
    integer(kind=i4b), public, allocatable, dimension(:) :: aperture_start
    ! To which psf does this aperture belong?
    integer(kind=i4b), public, allocatable, dimension(:) :: aperture_psf
    public :: aper_stop

contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine aper_stop()
        !--------------------------------------------------------------
        print *, "  * Stopping aperture module"
        if (allocated(aperture_size)) then
            deallocate (aperture_size)
            deallocate (aperture_start)
            deallocate (ap_hist_dim)
        end if
        print *, "  * Aperture module stopped"

    end subroutine aper_stop

end module aperture

!######################################################################
!######################################################################
!######################################################################

module aperture_boxed
    !contains the aperture functions for the square pixel boxed apertures
    use numeric_kinds
    implicit none
    private

    real(kind=dp), private, dimension(:, :), allocatable :: ap_box_size, ap_box_begin
    real(kind=dp), private, dimension(:), allocatable :: ap_box_idx, ap_box_idy
    real(kind=dp), private, dimension(:), allocatable :: ap_box_rot
    integer(kind=i4b), private, dimension(:), allocatable :: ap_box_bx

    ! read in boxed aperture file
    public :: aperture_boxed_readfile

    !figure out in which aperture the points fit.
    public :: aperture_boxed_find

    public :: aper_boxed_stop

    public :: aperture_boxed_field

contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine aper_boxed_stop()
        !----------------------------------------------------------------------
        if (allocated(ap_box_size)) then
            deallocate (ap_box_size, ap_box_begin, ap_box_idx, ap_box_idy, ap_box_rot)
            deallocate (ap_box_bx)
        end if

    end subroutine aper_boxed_stop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine aperture_boxed_readfile(handle, aper_n)
        use initial_parameters, only: conversion_factor
        use aperture, only: aperture_size, aperture_start
        use file_tools, only: next_content_line
        integer(kind=i4b), intent(in) :: handle, aper_n
        !----------------------------------------------------------------------
        integer(kind=i4b), save       :: amount = 0
        integer(kind=i4b)             :: biny
        character(len=80)             :: string
        !temporary array's
        real(kind=dp), Dimension(:, :), allocatable     :: tr
        integer(kind=i4b), Dimension(:, :), allocatable :: ti

        print *, "  * Reading boxed aperture file."

        amount = amount + 1

        if (allocated(ap_box_size)) then
            allocate (tr(amount, 7), ti(amount, 1))
            tr(:, :) = 0
            ti(:, :) = 0
            tr(1:amount - 1, 1:2) = ap_box_size(:, 1:2)
            tr(1:amount - 1, 3:4) = ap_box_begin(:, 1:2)
            tr(1:amount - 1, 5) = ap_box_idx(:)
            tr(1:amount - 1, 6) = ap_box_idy(:)
            tr(1:amount - 1, 7) = ap_box_rot(:)
            ti(1:amount - 1, 1) = ap_box_bx(:)

            deallocate (ap_box_size, ap_box_begin, ap_box_bx)
            deallocate (ap_box_idx, ap_box_idy, ap_box_rot)

            allocate (ap_box_size(amount, 2), ap_box_begin(amount, 2))
            allocate (ap_box_bx(amount), ap_box_idx(amount))
            allocate (ap_box_idy(amount), ap_box_rot(amount))

            ap_box_size(:, 1:2) = tr(:, 1:2)
            ap_box_begin(:, 1:2) = tr(:, 3:4)
            ap_box_idx(:) = tr(:, 5)
            ap_box_idy(:) = tr(:, 6)
            ap_box_rot(:) = tr(:, 7)
            ap_box_bx(:) = ti(:, 1)

            deallocate (tr, ti)
        else
            allocate (ap_box_size(1, 2), ap_box_begin(1, 2), ap_box_bx(1))
            allocate (ap_box_idx(1), ap_box_idy(1), ap_box_rot(1))
        end if

        print *, "  *  Reading box info"
        print *, "  *  Order: begin(x,y)"
        string = next_content_line(handle)
        read (string, fmt=*) ap_box_begin(amount, 1:2)
        print *, "      size(x,y) "
        string = next_content_line(handle)
        read (string, fmt=*) ap_box_size(amount, 1:2)
        print *, "      rotation"
        string = next_content_line(handle)
        read (string, fmt=*) ap_box_rot(amount)
        ap_box_rot(amount) = ap_box_rot(amount)*(pi_d/180.0_dp)
        print *, "      bin(x,y)"
        string = next_content_line(handle)
        read (string, fmt=*) ap_box_bx(amount), biny

        ! convert arcsec into km
        ap_box_begin(amount, :) = ap_box_begin(amount, :)*conversion_factor
        ap_box_size(amount, :) = ap_box_size(amount, :)*conversion_factor

        ap_box_idx(amount) = (ap_box_bx(amount)/ap_box_size(amount, 1))
        ap_box_idy(amount) = (biny/ap_box_size(amount, 2))

        aperture_start(aper_n) = amount
        aperture_size(aper_n) = ap_box_bx(amount)*biny

        print *, "   Total bins ", aperture_size(aper_n)
        print *, "   begin      ", ap_box_begin(amount, :)
        print *, "   size       ", ap_box_size(amount, :)
        print *, "   rotation   ", ap_box_rot(amount)
        print *, "   binx       ", ap_box_bx(amount)
        print *, "   idx,y      ", ap_box_idx(amount), ap_box_idy(amount)
        print *, " "
        print *, "  * Finished reading aperture"

    end subroutine aperture_boxed_readfile

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine aperture_boxed_find(ap, vec, res)
        ! vec is a n*m*2 matrix with vectors.
        ! res is an n*m matrix with has the resulting pixel of each vector
        ! ap is the aperture number.
        use aperture, only: aperture_start
        !use initial_parameters, only : psi_view
        use projection, only: psi_proj
        integer(kind=i4b), intent(in)                          :: ap
        real(kind=dp), dimension(:, :), intent(in)              :: vec
        integer(kind=i4b), dimension(size(vec, 1)), intent(out) :: res
        !----------------------------------------------------------------------
        integer(kind=i4b) :: n, j, bx
        real(kind=dp) :: r1, r2, b1, b2, idx, idy, sx, sy, x, y, t, q
        !real (kind=dp), dimension(size(vec,1)) :: t, q, x, y

        ! The number of this aperture in the memory
        n = aperture_start(ap)

        r1 = cos(-ap_box_rot(n) + pio2_d - psi_proj)
        r2 = sin(-ap_box_rot(n) + pio2_d - psi_proj)

        b1 = ap_box_size(n, 1)
        b2 = ap_box_size(n, 2)
        idx = ap_box_idx(n)
        idy = ap_box_idy(n)
        bx = ap_box_bx(n)
        sx = ap_box_begin(n, 1)
        sy = ap_box_begin(n, 2)

        ! Perform shift after rotation MC, 19/APR/2004
        ! Meanning of ap_box_begin has changed!

        do j = 1, size(vec, 1)
            t = vec(j, 1)
            q = vec(j, 2)
            x = t*r1 - q*r2 - sx
            res(j) = 0!_i4b
            if (x > 0.0_dp .and. x < b1) then
                y = t*r2 + q*r1 - sy
                if (y > 0.0_dp .and. y < b2) then
                    res(j) = int(x*idx) + int(y*idy)*bx + 1!_i4b
                end if
            end if
        end do

    end subroutine aperture_boxed_find

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine aperture_boxed_field(ap, x, y)
        use aperture, only: aperture_start
        integer(kind=i4b), intent(in) :: ap
        real(kind=dp), intent(out), dimension(:) :: x, y
        !----------------------------------------------------------------------
        integer(kind=i4b) :: i, j, n
        real(kind=dp) :: d, e, f, g, r1, r2

        ! The number of this aperture in the memory
        n = aperture_start(ap)
        r1 = cos(ap_box_rot(n))
        r2 = sin(ap_box_rot(n))
        x(:) = 0.0_dp
        y(:) = 0.0_dp
        do i = 1, 2
            f = (i - 1)*ap_box_size(n, 1) + ap_box_begin(n, 1)
            do j = 1, 2
                g = (j - 1)*ap_box_size(n, 2) + ap_box_begin(n, 2)
                d = r1*f - r2*g
                e = r2*f + r1*g
                x(1) = min(d, x(1))
                y(1) = min(e, y(1))
                x(2) = max(d, x(2))
                y(2) = max(e, y(2))
            end do
        end do

    end subroutine aperture_boxed_field

end module aperture_boxed

!######################################################################
!######################################################################
!######################################################################

module aperture_routines
    ! Basic aperture routines
    ! This module is the overhead to call the aperture_* functions
    use numeric_kinds
    use aperture
    implicit none
    private

    public :: aperture_setup

    public :: aperture_stop

    ! Finds the boundaries (on the sky) of an aperture.
    public :: aperture_field

contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine aperture_setup()
        use aperture_boxed, only: aperture_boxed_readfile
        use psf, only: psf_n
        !----------------------------------------------------------------------
        integer(kind=i4b)  :: i, handle = 11
        character(len=80) :: file, string
        print *, "  **Aperture setup module"
        print *, "  * How many different apertures?  :"
        read *, aperture_n

        allocate (aperture_size(aperture_n))
        allocate (aperture_start(aperture_n))
        allocate (aperture_psf(aperture_n))
        allocate (ap_hist_dim(aperture_n))
        print *, "  * using ", aperture_n, " aperture(s)"

        ap_hist0d_n = 0
        ap_hist2d_n = 0
        do i = 1, aperture_n
            print *, "  * What's the filename of the ", i, " aperture file ? :"
            read *, file
            print *, "  * Reading ", file

            open (unit=handle, file=file, action="read", status="old"&
                 &, position="rewind")
            string = "#counter_rotation_boxed_aperturefile_version_2"
            print *, "  * Assuming type ", string
            call aperture_boxed_readfile(handle, i)
            close (unit=handle)
            print *, "  * To which psf does this aperture belong?"
            read *, aperture_psf(i)
            if (aperture_psf(i) < 1 .or. aperture_psf(i) > psf_n) then
                stop " That PSF does not exist!"
            end if

            print *, "  * Histogram dimensions for this aperture (0, 1, or 2)?"
            read *, ap_hist_dim(i)
            print *, "  * The histograms are ", ap_hist_dim(i), " dimensional."
            if (ap_hist_dim(i) < 0 .or. ap_hist_dim(i) > 2) then
                stop "  Histogram dimension must be 0, 1, or 2!"
            end if
            if (ap_hist_dim(i) == 0) ap_hist0d_n = ap_hist0d_n + 1
            if (ap_hist_dim(i) == 2) ap_hist2d_n = ap_hist2d_n + 1
        end do
        print *, "  ** aperture setup finished"

    end subroutine aperture_setup

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine aperture_stop()
        use aperture, only: aper_stop
        use aperture_boxed, only: aper_boxed_stop
        !----------------------------------------------------------------------
        call aper_stop()
        call aper_boxed_stop()

    end subroutine aperture_stop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine aperture_field(pf, minx, maxx, miny, maxy)
        use aperture, only: aperture_n, aperture_psf
        use aperture_boxed, only: aperture_boxed_field
        use psf, only: psf_n
        integer(kind=i4b), intent(in) :: pf
        real(kind=dp), intent(out) :: maxx, minx, miny, maxy
        !----------------------------------------------------------------------
        logical, save :: initialized = .false.
        real(kind=dp), dimension(:, :), allocatable, save :: f_x, f_y
        integer(kind=i4b) :: i, pfn
        real(kind=dp), dimension(2) :: x, y

        if (.not. initialized) then
            allocate (f_x(psf_n, 2), f_y(psf_n, 2))
            f_x(:, :) = 0.0_dp
            f_y(:, :) = 0.0_dp
            do i = 1, aperture_n
                pfn = aperture_psf(i)
                call aperture_boxed_field(i, x, y)
                f_x(pfn, 1) = min(f_x(pfn, 1), x(1))
                f_x(pfn, 2) = max(f_x(pfn, 2), x(2))
                f_y(pfn, 1) = min(f_y(pfn, 1), y(1))
                f_y(pfn, 2) = max(f_y(pfn, 2), y(2))
            end do
            initialized = .true.
        end if

        minx = f_x(pf, 1)
        miny = f_y(pf, 1)
        maxx = f_x(pf, 2)
        maxy = f_y(pf, 2)

    end subroutine aperture_field

end module aperture_routines

!######################################################################
!######################################################################
!######################################################################

! $Id: orblib_f.f90,v 1.3 2011/10/25 08:48:45 bosch Exp $

! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! Sterrenwacht Leiden, The Netherlands

! This Module provides the binning of seperate histgrams.
! Type 0 : Do nothing.
! Type 1 : Sum all the bins with the same order number.

! This module was written to accomodate the Sauron pixel binning.

module binning
    ! Extension to the histogram module.
    ! This module takes care of possible binning of the histogram.
    use numeric_kinds
    implicit none
    private

    ! set's up the binning array's
    public :: binning_setup

    ! deallocate memory
    public :: binning_stop

    ! bin the aperture
    public :: binning_bin
    public :: binning_bin_h2d

    ! function for binning of type 1
    private :: binning_add_it_up
    private :: binning_add_it_up_h2d

    ! Type of binning. (0=no binning) (1=simple binning)
    integer(kind=i4b), private, allocatable, dimension(:)   ::  bin_type

    ! The way the boxes should be binned (order,ap)
    integer(kind=i4b), private, allocatable, dimension(:, :) ::  bin_order

    ! The total amount of boxes in the binned boxes
    ! ( actually maxval(binning_order(ap)) )
    integer(kind=i4b), public, allocatable, dimension(:) :: bin_max

    ! size(bin_order(:,ap),1)
    integer(kind=i4b), private, allocatable, dimension(:)   ::  bin_size

contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine binning_stop()
        !----------------------------------------------------------------------
        if (allocated(bin_type)) then
            deallocate (bin_type)
            deallocate (bin_order)
            deallocate (bin_max)
            deallocate (bin_size)
        end if

    end subroutine binning_stop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine binning_setup()
        use aperture, only: aperture_n
        use file_tools, only: next_content_line
        !----------------------------------------------------------------------
        integer(kind=i4b) :: i
        character(len=80) :: string
        print *, "  * Starting Binning setup"
        allocate (bin_type(aperture_n))
        allocate (bin_max(aperture_n))
        allocate (bin_size(aperture_n))
        bin_max(:) = 0
        bin_size(:) = 0

        do i = 1, aperture_n
            print *, "  * What kind of binning for aperture ", i
            print *, "    0=none 1=added up"
            do
                read *, bin_type(i)
                if (bin_type(i) == 1 .or. bin_type(i) == 0) exit
                print *, "   - Input incorrect try again!"
            end do
            print *, "  * Type:", bin_type(i)
        end do

        print *, "  * Reading binning files"

        do i = 1, aperture_n
            if (bin_type(i) == 1) then
                print *, "  * Aperture: ", i
                print *, "  * Give the filename of the binning file."
                read *, string
                print *, "  * Opening: ", string
                open (unit=30 + i, file=string, action="read", status="old"&
                  &, position="rewind")
                string = next_content_line(30 + i)  ! skip comment lines
                read (string, *) bin_size(i)
                print *, "  * bins in this aperture:", bin_size(i)
            end if
        end do

        allocate (bin_order(maxval(bin_size(:)), aperture_n))

        bin_order(:, :) = 0

        do i = 1, aperture_n
            if (bin_type(i) == 1) then
                print *, "  * Reading data of aperture:", i, bin_size(i)
                read (unit=30 + i, fmt=*) bin_order(1:bin_size(i), i)
                close (unit=30 + i)
            end if
        end do

        do i = 1, aperture_n
            if (bin_type(i) == 1) then
                bin_max(i) = maxval(bin_order(:, i))
            end if
        end do
        print *, "  ** Binning module setup finished."

    end subroutine binning_setup

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine binning_bin(ap, h, newsize)
        integer(kind=i4b), intent(in)                     :: ap
        real(kind=dp), intent(in out), dimension(:, :) :: h
        integer(kind=i4b), intent(out)                    :: newsize
        !----------------------------------------------------------------------
        if (bin_type(ap) == 1) then
            call binning_add_it_up(ap, h, newsize)
        else
            newsize = size(h, 1)
        end if

    end subroutine binning_bin

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine binning_bin_h2d(ap, h, newsize)  ! version for 2d histograms
        integer(kind=i4b), intent(in)                     :: ap
        real(kind=dp), intent(in out), dimension(:, :, :) :: h
        integer(kind=i4b), intent(out)                    :: newsize
        !----------------------------------------------------------------------
        if (bin_type(ap) == 1) then
            call binning_add_it_up_h2d(ap, h, newsize)
        else
            newsize = size(h, 1)
        end if

    end subroutine binning_bin_h2d

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine binning_add_it_up(ap, h, newsize)
        integer(kind=i4b), intent(in)                        :: ap
        integer(kind=i4b), intent(out)                       :: newsize
        real(kind=dp), intent(in out), dimension(:, :)    :: h
        !----------------------------------------------------------------------
        real(kind=dp), dimension(0:bin_max(ap), size(h, 2)) :: t
        integer(kind=i4b)                                    :: i

        newsize = bin_max(ap)
        t(:, :) = 0.0_dp
        ! check boundaries
        if (newsize > size(h, 1)) stop "Error: binning_add_it_up: new bin&
             &s are bigger then the original "
        if (size(h, 1) /= bin_size(ap)) stop " Wrong number of bins in a bin"

        do i = 1, size(h, 1)
            t(bin_order(i, ap), :) = t(bin_order(i, ap), :) + h(i, :)
        end do

        ! If you assume nothing, there is no way to do this without
        ! copying it back.
        h(1:newsize, :) = t(1:newsize, :)

    end subroutine binning_add_it_up

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine binning_add_it_up_h2d(ap, h, newsize)  ! version for 2d histograms
        integer(kind=i4b), intent(in)                        :: ap
        integer(kind=i4b), intent(out)                       :: newsize
        real(kind=dp), intent(in out), dimension(:, :, :)    :: h
        !----------------------------------------------------------------------
        real(kind=dp), dimension(0:bin_max(ap), size(h, 2), size(h, 3)) :: t
        integer(kind=i4b)                                    :: i

        newsize = bin_max(ap)
        t(:, :, :) = 0.0_dp
        ! check boundaries
        if (newsize > size(h, 1)) stop "Error: binning_add_it_up: new bin&
             &s are bigger then the original "
        if (size(h, 1) /= bin_size(ap)) stop " Wrong number of bins in a bin"

        do i = 1, size(h, 1)
            t(bin_order(i, ap), :, :) = t(bin_order(i, ap), :, :) + h(i, :, :)
        end do

        ! If you assume nothing, there is no way to do this without
        ! copying it back.
        h(1:newsize, :, :) = t(1:newsize, :, :)

    end subroutine binning_add_it_up_h2d

end module binning

!######################################################################
!######################################################################
!######################################################################

module histograms
! Routines for histogram manipulation
    use numeric_kinds
    implicit none
    private

    ! histogram data (aperture,vel)
    real(kind=dp), Dimension(:, :), private, allocatable  :: histogram
    ! 2d histogram data (aperture,vel1,vel2)
    real(kind=dp), Dimension(:, :, :), private, allocatable  :: histogram2d
    ! hist_basic(n,i) n=aperture number, i=width,center,#bins, dim
    real(kind=dp), Dimension(:, :, :), public, allocatable  :: hist_basic
    ! Are the velocity bins all the same?
    logical, public                                      :: hist_thesame
    !h_beg(n,dim),h_end(n,dim) n=aperture_n: begin/end of histogram
    !h_bin(n,dim),h_width(n,dim) : amount of / width of histogram pixels
    real(kind=dp), Dimension(:, :), private, allocatable  :: h_beg, h_end, h_width
    integer(kind=i4b), Dimension(:, :), private, allocatable:: h_bin
    ! h_start(n)  :  where start the first histogram of aperture n
    integer(kind=i4b), Dimension(:), private, allocatable:: h_start
    ! number of polygons/bins in each histogram for each aperture
    integer(kind=i4b), Dimension(:), private, allocatable:: h_blocks
    ! number of histograms (set = aperture_n)
    integer(kind=i4b), private                           :: h_n
    ! number of points stored in histogram ( Used in normalising. )
    real(kind=dp), Dimension(:), private, allocatable     :: h_n_stored
    ! total number of histograms/constraints  after binning
    integer(kind=i4b), private                           :: h_nconstr
    ! maximum histogram dimension
    integer(kind=i4b), public                            :: h_maxdim

    ! routines for writing histogram part of output files
    public :: histogram_write, histogram_setup_write
    public :: histogram_write_compat_sparse, histogram_setup_write_mass

    ! Store velocities in the histogram(n).
    public :: histogram_store

    ! Calculate the velocity bin from the velocity distribution
    public :: histogram_velbin

    ! function to reset the histogram for the next orbit
    public :: histogram_reset

    public :: histogram_stop

    public :: histogram_setup

contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram_reset()
        use aperture, only: aperture_n, ap_hist2d_n
        !----------------------------------------------------------------------
        if (aperture_n - ap_hist2d_n > 0) histogram(:, :) = 0.0_dp
        if (ap_hist2d_n > 0) histogram2d(:, :, :) = 0.0_dp
        h_n_stored(:) = 0.0_dp

    end subroutine histogram_reset

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram_setup_write(handle)
        integer(kind=i4b), intent(in) :: handle
        !----------------------------------------------------------------------
        integer(kind=i4b) :: t1
        !write information about the kinematical constraints and velocity histogram
        ! original names: nconstr,nvcube,dvcube
        t1 = hist_basic(1, 3, 1)/2.0_sp ! corrected by Remco 20/JAN/2003
        write (unit=handle) h_nconstr, t1, hist_basic(1, 1, 1)/hist_basic(1, 3, 1)

    end subroutine histogram_setup_write

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram_setup_write_mass(handle)
        integer(kind=i4b), intent(in) :: handle
        !----------------------------------------------------------------------
        write (unit=handle, fmt="(i5)") h_nconstr

    end subroutine histogram_setup_write_mass

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram_write(handle, handle_pops, handle_pm)
        use aperture, only: ap_hist_dim
        use binning, only: binning_bin, binning_bin_h2d
        integer(kind=i4b), intent(in) :: handle, handle_pops, handle_pm
        !----------------------------------------------------------------------
        integer(kind=i4b)            :: i, bg, ed
        print *, "  * Normalising and Writing histogram data."

        where (h_n_stored(:) > 0.0_dp)
            h_n_stored(:) = 1.0_dp/h_n_stored(:)
        elsewhere
            h_n_stored(:) = 0.0_dp
        end where

        do i = 1, h_n
            bg = h_start(i)
            ed = h_blocks(i) + h_start(i) - 1
            if (ap_hist_dim(i) < 2) then  ! 0d or 1d histograms
                call binning_bin(i, histogram(bg:ed, 1:h_bin(i, 1)), ed)
                ed = h_start(i) - 1 + ed
                !conversion normalizing
                histogram(bg:ed, 1:h_bin(i, 1)) = &
                    &h_n_stored(i)*histogram(bg:ed, 1:h_bin(i, 1))
                if (ap_hist_dim(i) == 0) then
                    if (h_bin(i, 1) /= 1) stop " 0d histogram must have 1 bin only."
                    write (unit=handle_pops) histogram(bg:ed, 1:h_bin(i, 1))
                else  ! 1d histograms
                    call histogram_write_compat_sparse(handle, i, histogram(bg:ed, 1:h_bin(i, 1)))
                end if
            else  ! 2d histograms
                call binning_bin_h2d(i, histogram2d(bg:ed, 1:h_bin(i, 1), 1:h_bin(i, 2)), ed)
                ed = h_start(i) - 1 + ed
                !conversion normalizing
                histogram2d(bg:ed, 1:h_bin(i, 1), 1:h_bin(i, 2)) = &
                    &h_n_stored(i)*histogram2d(bg:ed, 1:h_bin(i, 1), 1:h_bin(i, 2))
                call histogram2d_write_compat_sparse(handle_pm, i, &
                    &histogram2d(bg:ed, 1:h_bin(i, 1), 1:h_bin(i, 2)))
            end if
        end do
    end subroutine histogram_write

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram_write_compat_sparse(handle, i_hist, t)
        integer(kind=i4b), intent(in) :: handle, i_hist
        real(kind=dp), dimension(:, :), intent(in) :: t
        !----------------------------------------------------------------------
        integer(kind=i4b) :: ap, b, e, i, k, bout, eout
        do ap = 1, size(t, 1)
            b = 2*hist_basic(i_hist, 3, 1)
            e = -2*hist_basic(i_hist, 3, 1)
            do i = 1, size(t, 2)
                if (t(ap, i) > 0.0_dp) then
                    b = min(b, i)
                    e = max(e, i)
                end if
            end do

            ! write the relevant information for all velocity histograms to file
            k = hist_basic(i_hist, 3, 1)/2.0_sp + 1.0_sp
            bout = b - k
            eout = e - k
            write (unit=handle) bout, eout
            if (b <= e) write (unit=handle) t(ap, b:e)
        end do

    end subroutine histogram_write_compat_sparse

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram2d_write_compat_sparse(handle, aperture, t)
        integer(kind=i4b), intent(in) :: handle, aperture
        real(kind=dp), dimension(:, :, :), intent(in) :: t
        !----------------------------------------------------------------------
        integer(kind=i4b) :: ap, i, k
        integer(kind=i4b), dimension(2) :: min_idx, max_idx, center_bin
        integer(kind=i4b), dimension(2) :: begin_out, end_out
        do ap = 1, size(t, 1)
            min_idx = h_bin(aperture, :) + 1
            max_idx = -1
            do i = 1, size(t, 2)
                do k = 1, size(t, 3)
                    if (t(ap, i, k) > 0.0_dp) then
                        min_idx = min(min_idx, [i, k])
                        max_idx = max(max_idx, [i, k])
                    end if
                end do
            end do
            ! write the relevant information for all velocity histograms to file
            center_bin = h_bin(aperture, :) / 2 + 1
            begin_out = min_idx - center_bin
            end_out = max_idx - center_bin
            write (unit=handle) begin_out, end_out
            if (all(min_idx <= max_idx)) then
                write (unit=handle) t(ap, min_idx(1):max_idx(1), min_idx(2):max_idx(2))
            end if
        end do

    end subroutine histogram2d_write_compat_sparse

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram_velbin(pf, vel1d, vel2d, bin1d, bin2d)
        use aperture, only: aperture_psf, ap_hist_dim
        integer(kind=i4b), intent(in) :: pf
        real(kind=dp), dimension(:), intent(in) :: vel1d
        real(kind=dp), dimension(:, :), intent(in) :: vel2d
        integer(kind=i4b), dimension(size(vel1d)), intent(out) :: bin1d
        integer(kind=i4b), dimension(size(vel2d, 1), size(vel2d, 2)), intent(out) :: bin2d
        !----------------------------------------------------------------------
        integer(kind=i4b) :: i, ap, bins, dim
        integer(kind=i4b), dimension(size(vel2d, 2)) :: bins2
        real(kind=dp) :: v, beg, width, hend
        real(kind=dp), dimension(size(vel2d, 2)) :: v2, beg2, width2, hend2
        ! find an aperture which is in this pf
        do i = 1, h_n
            if (aperture_psf(i) == pf) ap = i
        end do

        if (ap_hist_dim(ap) < 2) then  ! 0d or 1d histograms
            beg = h_beg(ap, 1)
            hend = h_end(ap, 1)
            width = h_width(ap, 1)
            bins = h_bin(ap, 1)
            do i = 1, size(vel1d)
                v = vel1d(i)
                if (v > beg) then
                    if (v < hend) then
                        ! photon lies within the velocity range
                        bin1d(i) = int(((v - beg)/width)) + 1
                    else
                        ! photon lies above the range
                        ! Assign photon to the last velocity bin.
                        bin1d(i) = bins
                    end if
                else
                    ! photon lies below the velocity range
                    ! assign to first bin
                    bin1d(i) = 1
                end if
            end do
        else  ! 2d histograms  FIXME: can be coded more efficiently, together with 0d, 1d
            beg2 = h_beg(ap, :)
            hend2 = h_end(ap, :)
            width2 = h_width(ap, :)
            bins2 = h_bin(ap, :)
            do i = 1, size(vel2d, 1)
                v2 = vel2d(i, :)
                do dim = 1, size(vel2d, 2)
                    if (v2(dim) > beg2(dim)) then
                        if (v2(dim) < hend2(dim)) then
                            ! photon lies within the velocity range
                            bin2d(i, dim) = int(((v2(dim) - beg2(dim))/width2(dim))) + 1
                        else
                            ! photon lies above the range
                            ! Assign photon to the last velocity bin.
                            bin2d(i, dim) = bins2(dim)
                        end if
                    else
                        ! photon lies below the velocity range
                        ! assign to first bin
                        bin2d(i, dim) = 1
                    end if
                end do
            end do
        end if

    end subroutine histogram_velbin

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram_store(ap, n, velb, velb2d, tot)
        use aperture, only: ap_hist_dim
        integer(kind=i4b), intent(in)  :: ap
        integer(kind=i4b)                                   :: i, k, v
        integer(kind=i4b), dimension(:), intent(in)  :: n
        integer(kind=i4b), dimension(size(n, 1)), intent(in)  :: velb
        integer(kind=i4b), dimension(size(n, 1), 2), intent(in)  :: velb2d
        integer(kind=i4b), intent(in)  :: tot
        integer(kind=i4b), dimension(size(velb2d, 2))       :: v2d
        !----------------------------------------------------------------------
        !update number of points stored (including points not stored)
        ! For normalising.
        h_n_stored(ap) = h_n_stored(ap) + tot

        if (ap_hist_dim(ap) < 2) then  ! 0d or 1d histograms
            do i = 1, size(n, 1)
                k = n(i)
                if (k /= 0) then
                    k = k + h_start(ap) - 1
                    v = velb(i)
                    histogram(k, v) = histogram(k, v) + 1.0_dp
                end if
            end do
        else  ! 2d histograms
            do i = 1, size(n, 1)
                k = n(i)
                if (k /= 0) then
                    k = k + h_start(ap) - 1
                    v2d = velb2d(i, :)
                    histogram2d(k, v2d(1), v2d(2)) = histogram2d(k, v2d(1), v2d(2)) + 1.0_dp
                end if
            end do
        end if

    end subroutine histogram_store

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram_stop()
        use binning, only: binning_stop
        !----------------------------------------------------------------------
        if (allocated(h_bin)) then
            deallocate (hist_basic, h_beg, h_end, h_bin, h_width, h_start)
            deallocate (h_n_stored, h_blocks)
            if (allocated(histogram)) deallocate (histogram)
            if (allocated(histogram2d)) deallocate (histogram2d)
        end if
        call binning_stop()

    end subroutine histogram_stop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine histogram_setup()
        use aperture, only: aperture_n, aperture_size, aperture_psf, ap_hist_dim, &
                            ap_hist2d_n
        use binning, only: binning_setup, bin_max
        use psf, only: psf_n, psf_hist_dim
        !----------------------------------------------------------------------
        integer(kind=i4b)  :: i, j, ap, h_bin_max, ap_size
        integer(kind=i4b), dimension(2)      :: h_bin_max2d
        integer(kind=i4b), dimension(psf_n)  :: psf_aperture
        real(kind=dp)  :: width, center, bins

        print *, "  * Starting Histogram module"
        h_n = aperture_n
        ! for each psf, determine associated aperture
        ! intended for use with DYNAMITE: assuming 1-1 relation
        do j = 1, aperture_n
            psf_aperture(aperture_psf(j)) = j
            psf_hist_dim(aperture_psf(j)) = ap_hist_dim(j)
        end do

        ! allocate memory for histograms
        ! h_start, h_n_stored, h_blocks are independent of dimension
        allocate (h_start(h_n), h_n_stored(h_n), h_blocks(h_n))
        ! hist_basic, h_beg, h_end, h_bin, h_width are dimension-specific
        h_maxdim = maxval(psf_hist_dim)
        allocate(hist_basic(h_n, 3, h_maxdim))
        allocate(h_beg(h_n, h_maxdim), h_end(h_n, h_maxdim))
        allocate(h_bin(h_n, h_maxdim), h_width(h_n, h_maxdim))

        ! read histogram metadata
        do i = 1, psf_n
            print *, "  * PSF ", i, " has ", psf_hist_dim(i), "d histograms"
            if (psf_hist_dim(i) == 0) then
                print *, "  * Setting width, center, #bins to 1., 0., 1"
                hist_basic(psf_aperture(i), 1, 1) = 1.0_dp
                hist_basic(psf_aperture(i), 2, 1) = 0.0_dp
                hist_basic(psf_aperture(i), 3, 1) = 1.0_dp
            else if (psf_hist_dim(i) <= 2) then
                do j = 1, psf_hist_dim(i)
                    print *, "  * Give for psf ", i, " the histogram width, center and"
                    print *, "    amount of bins for dimension ", j
                    read *, width, center, bins
                    print *, width, center, bins
                    if (width <= 0) stop " Width to small"
                    if (bins < 1) stop " Too few bins"
                    hist_basic(psf_aperture(i), 1, j) = width
                    hist_basic(psf_aperture(i), 2, j) = center
                    hist_basic(psf_aperture(i), 3, j) = bins
                end do
            else
                stop "  * Only 0d, 1d, and 2d histograms are supported"
            end if
        end do

        h_beg(:, :) = hist_basic(:, 2, :) - (0.5_dp*hist_basic(:, 1, :))
        h_end(:, :) = hist_basic(:, 2, :) + (0.5_dp*hist_basic(:, 1, :))
        h_bin(:, :) = hist_basic(:, 3, :)
        h_width(:, :) = hist_basic(:, 1, :)/hist_basic(:, 3, :)

        ! allocate memory for 0d and 1d histograms
        if (aperture_n - ap_hist2d_n > 0) then
            h_bin_max = 0
            ap_size = 0
            do ap = 1, h_n
                if (ap_hist_dim(ap) < 2) then
                    h_bin_max = max(h_bin_max, h_bin(ap, 1))
                    ap_size = ap_size + aperture_size(ap)
                end if
            end do
            ! FIXME: sum(aperture_size(:)) is too large
            ! allocate (histogram(sum(aperture_size(:)), h_bin_max))
            allocate (histogram(ap_size, h_bin_max))
            print *, "  * Histogram size : ", size(histogram), "=", size(histogram, 1), "*",&
                & size(histogram, 2)
        end if
        ! allocate memory for 2d histograms
        if (ap_hist2d_n > 0) then
            h_bin_max2d = 0
            ap_size = 0
            do ap = 1, h_n
                if (ap_hist_dim(ap) == 2) then
                    h_bin_max2d(:) = max(h_bin_max2d(:), h_bin(ap, :))
                    ap_size = ap_size + aperture_size(ap)
                end if
            end do
            ! FIXME: sum(aperture_size(:)) is too large
            ! allocate (histogram2d(sum(aperture_size(:)), h_bin_max2d(1), h_bin_max2d(2)))
            allocate (histogram2d(ap_size, h_bin_max2d(1), h_bin_max2d(2)))
            print *, "  * 2d Histogram size : ", size(histogram2d), "=", size(histogram2d, 1), "*",&
                & size(histogram2d, 2), "*", size(histogram2d, 3)
        end if

        h_blocks(:) = aperture_size(:)
        ! i = 1
        ! do ap = 1, h_n
        !     h_start(ap) = i
        !     i = i + h_blocks(ap)
        ! end do
        ! IMPORTANT
        !   0d and 1d histograms:
        !     h_start(ap) points into histogram(individual_aperture, velbin)
        !   2d histograms:
        !     h_start(ap) points into histogram2d(individual_aperture, velbin1, velbin2)
        i = 1  ! 1d histogram counter
        j = 1  ! 2d histogram counter
        do ap = 1, h_n
            if (ap_hist_dim(ap) < 2) then  ! 0d or 1d histograms
                h_start(ap) = i
                i = i + h_blocks(ap)
            else  ! 2d histograms
                h_start(ap) = j
                j = j + h_blocks(ap)
            end if
        end do

        call histogram_reset()
        call binning_setup()

        ! Figure out how many histograms there are.
        ! adapted to DYNAMITE: only 1d histograms are counted
        ! (required by LegacyWeightSolver)
        h_nconstr = 0
        do ap = 1, aperture_n
            if (ap_hist_dim(ap) == 1) then
                if (bin_max(ap) == 0) then
                    h_nconstr = h_nconstr + aperture_size(ap)
                else
                    h_nconstr = h_nconstr + bin_max(ap)
                end if
            end if
        end do

        ! Figure out if all the 0d and 1d velocity histograms are the same.
        hist_thesame = .true.
        do j = 1, aperture_n  ! find first 0d or 1d histogram
            if (ap_hist_dim(j) <= 1) then
                width = hist_basic(j, 1, 1)
                center = hist_basic(j, 2, 1)
                bins = hist_basic(j, 3, 1)
                ap = j
                exit
            end if
        end do
        do j = ap + 1, aperture_n  ! compare with the rest
            if (ap_hist_dim(j) <= 1) then
                if (width /= hist_basic(j, 1, 1)) hist_thesame = .false.
                if (center /= hist_basic(j, 2, 1)) hist_thesame = .false.
                if (bins /= hist_basic(j, 3, 1)) hist_thesame = .false.
                exit
            end if
        end do
        if (hist_thesame) then
            print *, "  * All LOSVD velocity-bins and 0d hist bins are the same"
        else
            print *, "  * LOSVD velocity-bins are not the same. The standard NNLS will not"
            print *, "  * understand the ouput correctly (exception: pops data "
            print *, "  * velocity-bins may differ)."
        end if
        print *, "  ** Histogram module setup finished"

    end subroutine histogram_setup

end module histograms

!######################################################################
!######################################################################
!######################################################################

! $Id: orblib_f.f90,v 1.3 2011/10/25 08:48:45 bosch Exp $

! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! Sterrewacht Leiden, The Netherlands

module quadrantgrid
    use numeric_kinds
    implicit none
    private

    real(kind=dp), private, allocatable, dimension(:, :, :, :) :: quadrant_light
    real(kind=dp), private, allocatable, dimension(:) :: quad_lr, quad_lr2
    real(kind=dp), private, allocatable, dimension(:) :: quad_lth, quad_ltan2th
    real(kind=dp), private, allocatable, dimension(:) :: quad_lph, quad_ltanph

    public  :: qgrid_stop
    public  :: qgrid_write
    public  :: qgrid_setup_write
    public  :: qgrid_reset
    ! Store points in the grid.
    public  :: qgrid_store
    public  :: qgrid_setup
contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine qgrid_stop()
        !----------------------------------------------------------------------
        if (allocated(quadrant_light)) then
            deallocate (quadrant_light)
            deallocate (quad_lr, quad_lth, quad_lph)
            deallocate (quad_lr2, quad_ltan2th, quad_ltanph)
        end if

    end subroutine qgrid_stop

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine qgrid_reset()
        !----------------------------------------------------------------------
        quadrant_light(:, :, :, :) = 0.0_dp
    end subroutine qgrid_reset

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine qgrid_setup()
        use initial_parameters, only: rLogMin, rLogMax, sigobs_km &
                                      , quad_nr, quad_nth, quad_nph
        !----------------------------------------------------------------------
        integer(kind=i4b) :: i
        print *, "  ** Octant grid module setup"

        print *, "  ** Grid dimension:"
        print *, quad_nr, quad_nth, quad_nph

        allocate (quadrant_light(16, quad_nph, quad_nth, quad_nr))
        allocate (quad_lr(quad_nr + 1), quad_lr2(quad_nr + 1))
        allocate (quad_lth(quad_nth + 1), quad_ltan2th(quad_nth + 1))
        allocate (quad_lph(quad_nph + 1), quad_ltanph(quad_nph + 1))

        ! Define a grid in such a way that the boundaries define all possible bins
        ! This also means that there are N+1 boundaries for N bins.

        do i = 2, quad_nr
            quad_lr(i) = 10.0_dp**(rlogmin + (rLogMax - rlogmin + alog10(0.5))*(i - 1.0) &
                                   /(quad_nr - 0.0))
        end do
        quad_lr(1) = 0.0_dp
        quad_lr(quad_nr + 1) = max(10.0_dp**rLogMax*100.0_dp, maxval(sigobs_km)*10.0_dp)

        ! make a lr_squared array for quick computation
        quad_lr2(:) = quad_lr(:)**2_dp

        ! Define the angular bins
        do i = 2, quad_nth
            quad_lth(i) = pio2_d*(i - 1.0_dp)/(quad_nth)
        end do
        quad_lth(1) = 0.0_dp
        quad_lth(quad_nth + 1) = pio2_d

        ! define the angular bins
        do i = 2, quad_nph
            quad_lph(i) = pio2_d*(i - 1.0_dp)/(quad_nph)
        end do
        quad_lph(1) = 0.0_dp
        quad_lph(quad_nph + 1) = pio2_d

        ! make a lr_squared and tan arrays for quick computation
        quad_lr2(:) = quad_lr(:)**2.0_dp
        quad_ltanph(:) = tan(quad_lph(:))
        quad_ltan2th(:) = tan(quad_lth(:))**2.0_dp

        call qgrid_reset()

    end subroutine qgrid_setup

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine qgrid_store(proj, vel, type)
      use initial_parameters, only :  Omega ! (BT)
        ! proj (n, (x,y,z) )
        real(kind=dp), dimension(:, :), intent(in) :: proj, vel
        integer(kind=i4b), intent(in)                       :: type
        !----------------------------------------------------------------------
        real(kind=dp)      :: r2, theta, phi, x, y, z, vx, vy, vz
        integer(kind=i4b) :: i, j, n1, n2, n3, store_type
        integer(kind=i4b), save ::ir = 1, ith = 1, iph = 1

        real (kind=dp),dimension(3,8,5) ::vsgn
        real (kind=dp),dimension(3,8) ::psgn


        ! Signs of the (vx,vy,vz) for each Projection and type of Orbit
        ! (BT) 8-fold symmetry same as before
        real(kind=dp), dimension(3, 8, 5), &
            parameter :: vsgn1 = reshape((/ &
                                        ! X tubes
                                        1, 1, 1, -1, 1, 1, 1, 1, -1, -1, 1, -1, &
                                        -1, -1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, &
                                        ! Y tubes
                                        1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, &
                                        -1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, &
                                        ! Z tubes
                                        1, 1, 1, 1, -1, -1, -1, -1, 1, -1, 1, -1, &
                                        1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, 1, &
                                        ! Boxed
                                        1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, &
                                        1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, &
                                        ! Stochastic
                                        1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, &
                                        1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1/), (/3, 8, 5/))

        !Signs of the x,y,z for each projection  :psgn( [x,y,z], project )
        real(kind=dp), dimension(3, 8), &
            parameter :: psgn1 = reshape((/ &
                                        1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, &
                                        1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1/), (/3, 8/))

        ! Signs of the (vx,vy,vz) for each Projection and type of Orbit
        ! (BT) 4-fold symmetry in case of rotating frame
        real (kind=dp),dimension(3,8,5),parameter :: vsgn2= reshape((/  &
             ! X tubes
             1 , 1 , 1    ,1 , 1 , 1  ,  -1 , -1 ,1  ,   -1 , -1 , 1 , &
             1 ,1 , -1    ,1 ,1 , -1  , -1 ,-1 ,-1  ,  -1 ,-1 , -1 , &
             ! Y tubes
             1 , 1 , 1    , 1 , 1 ,1  ,  1 , 1 ,-1  ,  1 , 1 , -1 , &
             -1 , -1 , 1   ,-1 , -1 ,1  , -1 ,-1 ,-1  , -1 ,-1 , -1 , &
             ! Z tubes
             1 , 1 , 1    , 1 ,1 , 1  , -1 ,-1 , 1  , -1 , -1 , 1 , &
             1 , 1 ,-1    , 1 ,1 ,-1  , -1 ,-1 ,-1  , -1 , -1 ,-1 , &
             ! Boxed
             1 , 1 , 1    ,1 , 1 , 1  , -1 ,-1 , 1  ,  -1 ,-1 , 1 , &
             1 , 1 ,-1    ,1 , 1 ,-1  , -1 ,-1 ,-1  ,  -1 ,-1 ,-1 , &
             ! Stochastic
             1 , 1 , 1    ,1 , 1 , 1  , -1 ,-1 , 1  ,  -1 ,-1 , 1 , &
             1 , 1 ,-1    ,1 , 1 ,-1  , -1 ,-1 ,-1  ,  -1 ,-1 ,-1 /),(/3,8,5/))

        !Signs of the x,y,z for each projection  :psgn( [x,y,z], project )
        real (kind=dp),dimension(3,8),parameter :: psgn2= reshape((/  &
             1 , 1 , 1   , 1 , 1 , 1   , -1 , -1 , 1 ,  -1 , -1 , 1 , &
             1 , 1 ,-1   , 1 , 1 ,-1  , -1 , -1 ,-1 ,  -1 , -1 ,-1 /),(/3,8/))

        ! Use 8-fold for non-rotating, but 4-fold for rotating (BT)
        vsgn=vsgn1
        psgn=psgn1
        if (Omega /= 0.0_dp ) then
           vsgn=vsgn2
           psgn=psgn2
        endif

        ! Hunt assumes open boundaries, but our boundaries are closed
        ! So we dont give the outer boundaries to hunt.

        n1 = size(quad_lr) - 1
        n2 = size(quad_lth) - 1
        n3 = size(quad_lph) - 1

        select case (type)
        case (1)
            store_type = 0
        case (3)
            store_type = 1
        case default
            store_type = 2
        end select

        do i = 1, size(proj, 1) ! loop over photons

            do j = 1, 8            ! loop over projections symmetries
                ! Find the one projection that is in the positive octant.
                ! FIXME: Doing this with a loop is stupid.

                x = proj(i, 1)*psgn(1, j)
                y = proj(i, 2)*psgn(2, j)
                z = proj(i, 3)*psgn(3, j)

                ! only store when the photon is in the positive octant
                ! x==0 and z==0 will throw a divide by zero error.
                if (x > 0.0_dp .and. y >= 0.0_dp .and. z > 0.0_dp) then
                    ! this if is only passed once for every photon.

                    vx = vel(i, 1)*vsgn(1, j, type)
                    vy = vel(i, 2)*vsgn(2, j, type)
                    vz = vel(i, 3)*vsgn(3, j, type)

                    r2 = (x*x + y*y + z*z)
                    theta = (x*x + y*y)/(z*z) ! sqrt atan
                    phi = y/x               ! atan

                    call hunt(quad_lr2(2:n1), n1 - 1, r2, ir)
                    call hunt(quad_ltan2th(2:n2), n2 - 1, theta, ith)
                    call hunt(quad_ltanph(2:n3), n3 - 1, phi, iph)

                    ! store properties of the photon in the grid
                    quadrant_light(1:13, iph + 1, ith + 1, ir + 1) = &
                        quadrant_light(1:13, iph + 1, ith + 1, ir + 1) + &
                        (/1.0_dp, x, y, z, vx, vy, vz, vx*vx, vy*vy, vz*vz, vx*vy, vy*vz, vz*vx/)

                    ! store orbit type
                    quadrant_light(14 + store_type, iph + 1, ith + 1, ir + 1) = &
                        quadrant_light(14 + store_type, iph + 1, ith + 1, ir + 1) + 1
                end if
            end do
        end do

    end subroutine qgrid_store

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine qgrid_setup_write(hdl)
        integer(kind=i4b), intent(in) :: hdl
        !----------------------------------------------------------------------
        ! Write the information about the meridional plane grid. .

        ! remember that N bins have N+1 boundaries
        write (unit=hdl) size(quadrant_light, 1), size(quad_lph) - 1, size(quad_lth) - 1, size(quad_lr) - 1
        write (unit=hdl) quad_lr(:)
        write (unit=hdl) quad_lth(:)
        write (unit=hdl) quad_lph(:)

    end subroutine qgrid_setup_write

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine qgrid_write(hdl)
        integer(kind=i4b), intent(in):: hdl
        real(kind=dp) :: norm
        !----------------------------------------------------------------------

        print *, "  * Writing intrisic moment octant"

        where (quadrant_light(1, :, :, :) /= 0.0_dp)
            quadrant_light(2, :, :, :) = quadrant_light(2, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(3, :, :, :) = quadrant_light(3, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(4, :, :, :) = quadrant_light(4, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(5, :, :, :) = quadrant_light(5, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(6, :, :, :) = quadrant_light(6, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(7, :, :, :) = quadrant_light(7, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(8, :, :, :) = quadrant_light(8, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(9, :, :, :) = quadrant_light(9, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(10, :, :, :) = quadrant_light(10, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(11, :, :, :) = quadrant_light(11, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(12, :, :, :) = quadrant_light(12, :, :, :)/quadrant_light(1, :, :, :)
            quadrant_light(13, :, :, :) = quadrant_light(13, :, :, :)/quadrant_light(1, :, :, :)
        end where

        ! Properly normalize the light
        ! by dividing by the total number of photons
        ! contributing to each grid element.

        norm = sum(quadrant_light(1, :, :, :))

        where (quadrant_light(1, :, :, :) /= 0.0_dp)
            quadrant_light(1, :, :, :) = quadrant_light(1, :, :, :)/norm
            quadrant_light(14, :, :, :) = quadrant_light(14, :, :, :)/norm ! orbtype
            quadrant_light(15, :, :, :) = quadrant_light(15, :, :, :)/norm ! orbtype
            quadrant_light(16, :, :, :) = quadrant_light(16, :, :, :)/norm ! orbtype
        end where

        ! write the light quadrant information
        write (unit=hdl) quadrant_light(:, :, :, :)

    end subroutine qgrid_write

end module quadrantgrid

!######################################################################
!######################################################################
!######################################################################

! $Id: orblib_f.f90,v 1.3 2011/10/25 08:48:45 bosch Exp $

! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! July 2002 Sterrewacht Leiden.

module output
! Module doing all the output of the program
    use numeric_kinds
    implicit none
    private

    integer(kind=i4b), private :: out_handle = 0_i4b
    character(len=80), public  :: out_file_qgrid, out_file_pm, out_file_pops &
                                  , out_file_losvd, out_file_orbclass
    character(len=84), private :: out_tmp_file

    public :: output_setup

    public :: output_close

    public :: output_write

contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine output_setup()
        use integrator, only: integrator_setup_write, integrator_set_current,&
             &                   integrator_current
        use aperture, only: aperture_n, ap_hist0d_n, ap_hist2d_n
        use histograms, only: histogram_setup_write
        use histograms, only: histogram_setup_write_mass
        use quadrantgrid, only: qgrid_setup_write
        !----------------------------------------------------------------------
        character(len=8)  :: d
        character(len=10)  :: t
        character(len=5)  :: g
        integer(kind=i4b)  :: error, tmp

        print *, "  ** Setting up output module"
        print *, "  * Give the name of the qgrid outputfile:"
        ! read (unit=*, fmt="(a80)"), out_file
        read *, out_file_qgrid
        print *, out_file_qgrid

        if (ap_hist0d_n > 0) then
            print *, "  * Give the name of the pops '0d histogram' outputfile:"
            read *, out_file_pops
            print *, out_file_pops
        end if

        if (aperture_n - ap_hist0d_n - ap_hist2d_n > 0) then
            print *, "  * Give the name of the 1d losvd histogram outputfile:"
            read *, out_file_losvd
            print *, out_file_losvd
        end if

        if (ap_hist2d_n > 0) then
            print *, "  * Give the name of the proper motions 2d histogram outputfile:"
            read *, out_file_pm
            print *, out_file_pm
        end if

        print *, "  * Give the name of the orbit classification outputfile:"
        read *, out_file_orbclass
        print *, out_file_orbclass

        ! out_file = adjustl(out_file)
        ! print *, out_file

        out_tmp_file = out_file_qgrid
        out_tmp_file(len_trim(out_file_qgrid) + 1:len_trim(out_file_qgrid) + 4) = ".tmp"
        print *, out_tmp_file

        call date_and_time(date=d, time=t, zone=g)
        print *, "  * Date : ", d, " ", t, " ", g

        out_handle = 50
        error = 0
        ! Check status and setup files
        open (unit=out_handle + 1, iostat=error, file=out_tmp_file, action="write", &
             & status="new", position="rewind")
        if (error == 0) then

            if (error /= 0) stop "  Error opening file."
            ! Write orbit library header in *binary* (typically orblib.dat)
            open (unit=out_handle, iostat=error, file=out_file_qgrid, action="write", &
                  status="new", form="unformatted")
            call integrator_setup_write(out_handle)
            call qgrid_setup_write(out_handle)
            close (unit=out_handle, iostat=error)
            if (error /= 0) stop "  Error closing qgrid file."
            if (ap_hist0d_n > 0) then
                open (unit=out_handle, iostat=error, file=out_file_pops, action="write", &
                    status="new", form="unformatted")
                close (unit=out_handle, iostat=error)  ! no setup, just create file
                if (error /= 0) stop "  Error closing pops file."
            end if
            if (aperture_n - ap_hist0d_n - ap_hist2d_n > 0) then
                open (unit=out_handle, iostat=error, file=out_file_losvd, action="write", &
                    status="new", form="unformatted")
                call histogram_setup_write(out_handle)  ! for LegacyWeightSolver only
                close (unit=out_handle, iostat=error)
                if (error /= 0) stop "  Error closing losvd file."
            end if
            if (ap_hist2d_n > 0) then
                open (unit=out_handle, iostat=error, file=out_file_pm, action="write", &
                  status="new", form="unformatted")
                close (unit=out_handle, iostat=error)  ! no setup, just create file
                if (error /= 0) stop "  Error closing proper motions file."
            end if

            ! Write status file
            write (unit=out_handle + 1, fmt=*, iostat=error) integrator_current
            if (error /= 0) stop "  Error writing to status file."
            close (unit=out_handle + 1, iostat=error)
            if (error /= 0) stop "  Error closing status file."
        else
            print *, "  * Trying to resume previous calculations"
            ! Try to read the status file
            open (unit=out_handle + 1, iostat=error, file=out_tmp_file, action="read", &
                 & status="old", position="rewind")
            if (error /= 0) stop "  Error: Inconsistent status file."
            read (unit=out_handle + 1, fmt=*) tmp
            if (tmp == -1) stop " Error: Orbit library already finished or orbit &
                 & library in inconsistent state"
            call integrator_set_current(tmp)
            close (unit=out_handle + 1, iostat=error)
            if (error /= 0) stop "  Error closing status file."

            ! Checking if orbit library file exists
            open (unit=out_handle, iostat=error, file=out_file_qgrid, action="write", &
                 & status="old", position="append", form="unformatted")
            if (error /= 0) stop "  Error opening library file. Does it exist?"
            close (unit=out_handle, iostat=error)
            if (error /= 0) stop "  Error closing library file."
            print *, "  * Resuming with orbit :", tmp + 1
        end if

        open (unit=30, file=out_file_orbclass, status="replace", action="write")

        print *, "  ** Output file setup finished."

    end subroutine output_setup

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine output_close()
        use integrator, only: integrator_current
        use aperture, only: aperture_n, ap_hist0d_n, ap_hist2d_n
        !----------------------------------------------------------------------
        integer :: error
        print *, "  * Closing files and stopping output module"
        if (out_handle /= 0) then
            open (unit=out_handle, iostat=error, file=out_file_qgrid, action="write", &
                 & status="old", position="append", form="unformatted")
            if (error /= 0) stop "  Error opening qgrid file."
            write (unit=out_handle, iostat=error) " "
            if (error /= 0) stop "  Error writing to qgrid file. Disk full?"
            close (unit=out_handle, iostat=error)
            if (error /= 0) stop "  Error closing qgrid file."
            if (ap_hist0d_n > 0) then
                open (unit=out_handle, iostat=error, file=out_file_pops, action="write", &
                    & status="old", position="append", form="unformatted")
                if (error /= 0) stop "  Error opening pops file."
                write (unit=out_handle, iostat=error) " "
                if (error /= 0) stop "  Error writing to pops file. Disk full?"
                close (unit=out_handle, iostat=error)
                if (error /= 0) stop "  Error closing pops file."
            end if
            if (aperture_n - ap_hist0d_n - ap_hist2d_n > 0) then
                open (unit=out_handle, iostat=error, file=out_file_losvd, action="write", &
                    & status="old", position="append", form="unformatted")
                if (error /= 0) stop "  Error opening losvd file."
                write (unit=out_handle, iostat=error) " "
                if (error /= 0) stop "  Error writing to losvd file. Disk full?"
                close (unit=out_handle, iostat=error)
                if (error /= 0) stop "  Error closing losvd file."
            end if
            if (ap_hist2d_n > 0) then
                open (unit=out_handle, iostat=error, file=out_file_pm, action="write", &
                    & status="old", position="append", form="unformatted")
                if (error /= 0) stop "  Error opening proper motions file."
                write (unit=out_handle, iostat=error) " "
                if (error /= 0) stop "  Error writing to proper motions file. Disk full?"
                close (unit=out_handle, iostat=error)
                if (error /= 0) stop "  Error closing proper motions file."
            end if
        end if

        ! Update the temp file to finished status
        open (unit=out_handle + 1, iostat=error, file=out_tmp_file, action="write", &
             & status="old", position="rewind")
        if (error /= 0) stop "  Error opening status file."
        write (unit=out_handle + 1, fmt=*, iostat=error) - 1_i4b, "orbit library &
             & finished ", integrator_current
        if (error /= 0) stop "  Error writing to status file."
        close (unit=out_handle + 1, iostat=error)
        if (error /= 0) stop "  Error closing status file."

        close(unit=32,iostat=error)                           ! for orbit info (BT)
        if (error/=0) stop "  Error closing status file."

        print *, " * Finished closing files"

    end subroutine output_close

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine output_write()
        use histograms, only: histogram_write
        use aperture, only: aperture_n, ap_hist0d_n, ap_hist2d_n
        use quadrantgrid, only: qgrid_write
        use integrator, only: integrator_write, integrator_current
        !----------------------------------------------------------------------
        integer :: error, out_handle_pops, out_handle_pm
        ! Update the temp file to writing status
        open (unit=out_handle + 1, iostat=error, file=out_tmp_file, action="write", &
             & status="old", position="rewind")
        if (error /= 0) stop "  Error opening status file."
        write (unit=out_handle + 1, fmt=*, iostat=error) - 1_i4b, "Writing orbit: ", &
             & integrator_current - 1
        if (error /= 0) stop "  Error writing to status file."
        close (unit=out_handle + 1, iostat=error)
        if (error /= 0) stop "  Error closing status file."

        ! Write the orbit to the *binary* output files
        ! (typically orblib_qgrid.dat, orblib_pops.dat, orblib_losvd_hist.dat,
        ! orblib_pm.dat).

        ! open and write to the qgrid file
        open (unit=out_handle, iostat=error, file=out_file_qgrid, action="write", &
             & status="old", position="append", form="unformatted")
        if (error /= 0) stop "  Error opening qgrid file."
        call integrator_write(out_handle)
        call qgrid_write(out_handle)
        close (unit=out_handle, iostat=error)
        if (error /= 0) stop "  Error closing qgrid file."

        ! open the pops file if 0d histograms exist
        if (ap_hist0d_n > 0) then
            out_handle_pops = out_handle + 10
            open (unit=out_handle_pops, iostat=error, file=out_file_pops, action="write", &
                & status="old", position="append", form="unformatted")
            if (error /= 0) stop "  Error opening pops file."
        else
            out_handle_pops = 0
        end if
        ! open the losvd file if 1d histograms exist
        if (aperture_n - ap_hist0d_n - ap_hist2d_n > 0) then
            open (unit=out_handle, iostat=error, file=out_file_losvd, action="write", &
                & status="old", position="append", form="unformatted")
            if (error /= 0) stop "  Error opening losvd file."
        end if
        ! open the proper motions file if 2d histograms exist
        if (ap_hist2d_n > 0) then
            out_handle_pm = out_handle + 20
            open (unit=out_handle_pm, iostat=error, file=out_file_pm, action="write", &
                & status="old", position="append", form="unformatted")
            if (error /= 0) stop "  Error opening proper motions file."
        else
            out_handle_pm = 0
        end if
        ! write 0d, 1d, and 2d (if existing) histograms
        call histogram_write(out_handle, out_handle_pops, out_handle_pm)
        ! close the pops file if 0d histograms exist
        if (ap_hist0d_n > 0) then
            close (unit=out_handle_pops, iostat=error)
            if (error /= 0) stop "  Error closing pops file."
        end if
        ! close the losvd file if 1d histograms exist
        if (aperture_n - ap_hist0d_n - ap_hist2d_n > 0) then
            close (unit=out_handle, iostat=error)
            if (error /= 0) stop "  Error closing losvd file."
        end if
        ! close the proper motions file if 2d histograms exist
        if (ap_hist2d_n > 0) then
            close (unit=out_handle_pm, iostat=error)
            if (error /= 0) stop "  Error closing proper motions file."
        end if

        ! Update the temp file to intermediate status
        open (unit=out_handle + 1, iostat=error, file=out_tmp_file, action="write", &
             & status="old", position="rewind")
        if (error /= 0) stop "  Error opening status file."
        write (unit=out_handle + 1, fmt=*, iostat=error) integrator_current
        if (error /= 0) stop "  Error writing to status file."
        close (unit=out_handle + 1, iostat=error)
        if (error /= 0) stop "  Error closing status file."

    end subroutine output_write

end module output

!######################################################################
!######################################################################
!######################################################################

! $Id: orblib_f.f90,v 1.3 2011/10/25 08:48:45 bosch Exp $

! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! July 2002 Sterrenwacht Leiden.

module high_level
    use numeric_kinds
    implicit none
    private

    ! setup/run/stop the program.
    public :: setup, setup_bar, run, stob

contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine setup()
        use integrator, only: integrator_setup
        use projection, only: projection_setup
        use quadrantgrid, only: qgrid_setup
        use aperture_routines, only: aperture_setup
        use histograms, only: histogram_setup
        use psf, only: psf_setup
        use output, only: output_setup
        !----------------------------------------------------------------------
        print *, "  ** Start Setup"
        call integrator_setup()
        call projection_setup()
        call qgrid_setup()
        call psf_setup()
        call aperture_setup()
        call histogram_setup()
        call output_setup()
        print *, "  ** Setup Finished"

    end subroutine setup

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine setup_bar()
        use integrator, only: integrator_setup_bar
        use projection, only: projection_setup
        use quadrantgrid, only: qgrid_setup
        use aperture_routines, only: aperture_setup
        use histograms, only: histogram_setup
        use psf, only: psf_setup
        use output, only: output_setup
        !----------------------------------------------------------------------
        print *, "  ** Start Setup"
        call integrator_setup_bar()
        call projection_setup()
        call qgrid_setup()
        call psf_setup()
        call aperture_setup()
        call histogram_setup()
        call output_setup()
        print *, "  ** Setup Finished"

    end subroutine setup_bar

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine run()
      use initial_parameters, only: Omega
        use histograms, only: histogram_reset, hist_thesame, &
                              histogram_velbin, histogram_store, h_maxdim
        use projection, only: project, projection_symmetry
        use integrator, only: integrator_integrate, integrator_points
        use output, only: output_write
        use quadrantgrid, only: qgrid_reset, qgrid_store
        use psf, only: psf_n, psf_gaussian
        use aperture, only: aperture_n, aperture_psf
        use aperture_boxed, only: aperture_boxed_find

        !----------------------------------------------------------------------
        logical :: done, first, alldone
        real(kind=dp), dimension(integrator_points, 3) :: pos
        real(kind=dp), dimension(integrator_points, 3) :: vel
        real(kind=dp), dimension(integrator_points*projection_symmetry, 2):: proj, vec_gauss
        real(kind=dp), dimension(integrator_points*projection_symmetry):: losvel
        real(kind=dp), dimension(integrator_points*projection_symmetry, 2):: vel2d
        integer(kind=i4b), dimension(integrator_points*projection_symmetry):: velb, poly
        integer(kind=i4b), dimension(integrator_points*projection_symmetry, 2):: velb2d
        integer(kind=i4b)                                          :: ap, i

        integer(kind=i4b) :: type
        real(kind=dp) :: t1, t2
        alldone = .false.
        print *, "  ** Starting Orbit Calculations"
        if (Omega /= 0.0_dp) print*,"Pattern speed ==================== ", Omega ! (BT)
        do  ! for each orbit

            call cpu_time(t1)

            call histogram_reset()
            call qgrid_reset()
            first = .true.
            do ! for all dithers
                call integrator_integrate(pos, vel, type, done, first, alldone)
                first = .false.
                if (done .or. alldone) exit

                call qgrid_store(pos(:, :), vel(:, :), type)
                first = .true.
                do ! for all projections
                    call project(type, pos, vel, proj, losvel, vel2d, done, first)
                    if (done) exit
                    first = .false.

                    if (hist_thesame .and. h_maxdim <= 1) then
                        call histogram_velbin(1, losvel, vel2d, velb, velb2d)
                    end if
                    do i = 1, psf_n
                        if (.not. hist_thesame .or. h_maxdim >= 2) then
                            call histogram_velbin(i, losvel, vel2d, velb, velb2d)
                        end if
                        call psf_gaussian(i, proj, vec_gauss)
                        do ap = 1, aperture_n
                            if (i == aperture_psf(ap)) then
                                call aperture_boxed_find(ap, vec_gauss, poly)
                                call histogram_store(ap, poly, velb, velb2d, size(proj, 1))
                            end if
                        end do
                    end do
                end do
            end do
            if (alldone) exit
            call output_write()
            call cpu_time(t2)
            print *, "  * Time spent one orbit:", t2 - t1, " seconds"
        end do
        print *, "  ** Finished Orbit Calculations"

    end subroutine run

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine stob()
        use integrator, only: integrator_stop
        use projection, only: projection_stop
        use histograms, only: histogram_stop
        use psf, only: psf_stop
        use aperture_routines, only: aperture_stop
        use output, only: output_close
        !----------------------------------------------------------------------
        call output_close()
        call integrator_stop()
        call projection_stop()
        call aperture_stop()
        call psf_stop()
        call histogram_stop()

    end subroutine stob

end module high_level
