! THIS VERSION: GALAHAD 2.2 - 27/05/2008 AT 14:00 GMT.

!-*-*-*-*-*-*-*-  G A L A H A D _ G L T R  double  M O D U L E  *-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released pre GALAHAD Version 1.0. April 21st 1997
!   update released with GALAHAD Version 2.0. February 28th 2006

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_GLTR_double

!      -----------------------------------------------
!      |                                             |
!      | Solve the quadratic program                 |
!      |                                             |
!      |    minimize     1/2 <x, H x> + <c, x> + f0  |
!      |    subject to   < x, M x> <= radius^2       |
!      |       or                                    |
!      |    subject to   < x, M x>  = radius^2       |
!      |                                             |
!      | using a generalized Lanczos method          |
!      |                                             |
!      -----------------------------------------------

      USE GALAHAD_SYMBOLS
      USE GALAHAD_SPACE_double
      USE GALAHAD_RAND_double
      USE GALAHAD_ROOTS_double, ONLY: ROOTS_quadratic
      USE GALAHAD_SPECFILE_double
      USE GALAHAD_NORMS_double, ONLY: TWO_NORM

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: GLTR_initialize, GLTR_read_specfile, GLTR_solve,               &
                GLTR_terminate, GLTR_leftmost_eigenvalue,                      &
                GLTR_leftmost_eigenvector

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: point1 = 0.1_wp
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: point9 = 0.9_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: two = 2.0_wp
      REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp
      REAL ( KIND = wp ), PARAMETER :: roots_tol = ten ** ( - 12 )

!--------------------------
!  Derived type definitions
!--------------------------

      TYPE, PUBLIC :: GLTR_control_type
        INTEGER :: error, out, print_level, itmax, Lanczos_itmax, extra_vectors
        REAL ( KIND = wp ) :: stop_relative, stop_absolute, fraction_opt, f_0
        REAL ( KIND = wp ) :: rminvr_zero
        LOGICAL :: unitm, steihaug_toint, boundary, equality_problem
        LOGICAL :: space_critical, deallocate_error_fatal
        CHARACTER ( LEN = 30 ) :: prefix
      END TYPE

      TYPE, PUBLIC :: GLTR_info_type
        INTEGER :: status, alloc_status, iter, iter_pass2
        REAL ( KIND = wp ) :: multiplier, mnormx, piv, curv, rayleigh
        LOGICAL :: negative_curvature
       CHARACTER ( LEN = 80 ) :: bad_alloc
      END TYPE

      TYPE, PUBLIC :: GLTR_data_type
        PRIVATE
        TYPE ( RAND_seed ) :: seed
        INTEGER :: iter, itm1, itmax, dim_sub, switch, titmax, tinfo, titer
        INTEGER :: Lanczos_itmax, extra_vectors
        REAL ( KIND = wp ) :: alpha, beta, rminvr, rminvr_old
        REAL ( KIND = wp ) :: stop, normp, x_last
        REAL ( KIND = wp ) :: diag, offdiag, rtol,  pgnorm, radius2
        REAL ( KIND = wp ) :: xmx, xmp, pmp, old_leftmost, tau
        LOGICAL :: printi, printd, interior, use_old, try_warm
        LOGICAL :: switch_to_Lanczos, prev_steihaug_toint, save_vectors
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: P
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: D
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: OFFD
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: ALPHAS
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RMINVRS
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: MIN_f
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: LAMBDA
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: D_fact
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: OFFD_fact
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: C_sub
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X_sub
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Z
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: W
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: R_extra
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: P_extra
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: V_extra
      END TYPE

!---------------------------------
!   I n t e r f a c e  B l o c k s
!---------------------------------

      INTERFACE PTTRF

        SUBROUTINE SPTTRF( n, D, E, info )
        INTEGER, INTENT( IN ) :: n
        INTEGER, INTENT( OUT ) :: info
        REAL, INTENT( INOUT ) :: D( n ), E( n - 1 )
        END SUBROUTINE SPTTRF

        SUBROUTINE DPTTRF( n, D, E, info )
        INTEGER, INTENT( IN ) :: n
        INTEGER, INTENT( OUT ) :: info
        DOUBLE PRECISION, INTENT( INOUT ) :: D( n ), E( n - 1 )
        END SUBROUTINE DPTTRF

      END INTERFACE 

    CONTAINS

!-*-*-*-*-*-  G L T R _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE GLTR_initialize( data, control )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  .  Set initial values for the GLTR control parameters  .
!   
!  Argument:
!  =========
!
!   data     private internal data
!   control  a structure containing control information. Components are -
!            stop_relative      the iteration stops successfully when the 
!            stop_absolute      gradient in the M(inverse) norm is less than 
!                        max( stop_relative * initial M(inverse)
!                             gradient norm, stop_absolute )
!            fraction_opt a direction which gives at least fraction_opt times 
!                        the optimal objective value will be found
!            unitm  .TRUE. if the matrix M is the identity
!            steihaug_toint  .TRUE. if we should stop when the Trust-region 
!                            is first encountered
!            space_critical  .TRUE. if allocated arrays are required to be 
!                            exactly the right size 
!            deallocate_error_fatal. .TRUE. if a deallocation error should
!                                    terminate execution
!            boundary .TRUE. if the solution is thought to lie on the boundary
!            equality_problem .TRUE. if the solution is REQUIRED to lie on the
!                             boundary (i.e., the constraint is an equality )
!            extra_vectors  the number of extra work vectors of length n used
!            out         output unit
!            print_level print level. > 0 for output
!            itmax       maximum number of iterations allowed
!            Lanczos_itmax   maximum number of Lanczos iterations allowed
!            prefix  output lines will be prefixed by 
!                        prefix(2:LEN(TRIM(%prefix))-1)
!               where prefix contains the required string enclosed in 
!               quotes, e.g. "string" or 'string'
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t
!-----------------------------------------------

      TYPE ( GLTR_DATA_TYPE ), INTENT( OUT ) :: data
      TYPE ( GLTR_CONTROL_TYPE ), INTENT( OUT ) :: control        

!  Initalize random number seed

      CALL RAND_initialize( data%seed )

!  Set initial control parameter values

      control%stop_relative = SQRT( EPSILON( one ) )
      control%stop_absolute = zero
      control%fraction_opt = one
      control%f_0 = zero
      control%rminvr_zero = ten * EPSILON( one )
      control%unitm = .TRUE.
      control%steihaug_toint = .FALSE.
      control%boundary = .FALSE.
      control%equality_problem = .FALSE.
      control%space_critical = .FALSE.
      control%deallocate_error_fatal  = .FALSE.
      control%error = 6
      control%out = - 1
      control%print_level = 0
      control%extra_vectors = 0
      control%itmax = - 1
      control%Lanczos_itmax = - 1
      control%prefix = '""     '

      data%prev_steihaug_toint = .TRUE.

      RETURN

!  End of subroutine GLTR_initialize

      END SUBROUTINE GLTR_initialize

!-*-*-*-   G L T R _ R E A D _ S P E C F I L E  S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE GLTR_read_specfile( control, device, alt_specname )

!  Reads the content of a specification file, and performs the assignment of 
!  values associated with given keywords to the corresponding control parameters

!  The defauly values as given by GLTR_initialize could (roughly) 
!  have been set as:

!  BEGIN GLTR SPECIFICATIONS (DEFAULT)
!   error-printout-device                           6
!   printout-device                                 6
!   print-level                                     0
!   maximum-number-of-iterations                    -1           
!   maximum-number-of-Lanczos-iterations            -1
!   number-extra-n-vectors-used                     0
!   relative-accuracy-required                      1.0E-8 
!   absolute-accuracy-required                      0.0
!   fraction-optimality-required                    1.0
!   constant-term-in-objective                      0.0
!   zero-gradient-tolerance                         2.0E-15
!   two-norm-trust-region                           T
!   stop-as-soon-as-boundary-encountered            F
!   solution-is-likely-on-boundary                  F
!   equality-problem                                F
!   space-critical                                  F
!   deallocate-error-fatal                          F
!   output-line-prefix                              ""
!  END GLTR SPECIFICATIONS

!  Dummy arguments

      TYPE ( GLTR_control_type ), INTENT( INOUT ) :: control        
      INTEGER, INTENT( IN ) :: device
      CHARACTER( LEN = 16 ), OPTIONAL :: alt_specname

!  Programming: Nick Gould and Ph. Toint, January 2002.

!  Local variables

      INTEGER, PARAMETER :: lspec = 34
      CHARACTER( LEN = 16 ), PARAMETER :: specname = 'GLTR            '
      TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec

!  Define the keywords

     spec%keyword = ''

!  Integer key-words

      spec(  1 )%keyword = 'error-printout-device'
      spec(  2 )%keyword = 'printout-device'
      spec(  3 )%keyword = 'print-level' 
      spec(  4 )%keyword = 'maximum-number-of-iterations'
      spec(  5 )%keyword = 'maximum-number-of-Lanczos-iterations'
      spec( 18 )%keyword = 'number-extra-n-vectors-used'

!  Real key-words

      spec(  6 )%keyword = 'relative-accuracy-required'
      spec(  7 )%keyword = 'absolute-accuracy-required'
      spec(  8 )%keyword = 'fraction-optimality-required'
      spec(  9 )%keyword = 'constant-term-in-objective'
      spec( 17 )%keyword = 'zero-gradient-tolerance'

!  Logical key-words

      spec( 10 )%keyword = 'two-norm-trust-region'
      spec( 11 )%keyword = 'stop-as-soon-as-boundary-encountered'
      spec( 12 )%keyword = 'solution-is-likely-on-boundary'
      spec( 13 )%keyword = 'equality-problem'
      spec( 14 )%keyword = 'space-critical'
      spec( 15 )%keyword = 'deallocate-error-fatal'

!  Character key-words

!     spec( 16 )%keyword = 'output-line-prefix'

!  Read the specfile

      IF ( PRESENT( alt_specname ) ) THEN
        CALL SPECFILE_read( device, alt_specname, spec, lspec, control%error )
      ELSE
        CALL SPECFILE_read( device, specname, spec, lspec, control%error )
      END IF

!  Interpret the result

!  Set integer values

      CALL SPECFILE_assign_value( spec( 1 ), control%error,                    &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 2 ), control%out,                      &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 3 ), control%print_level,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 4 ), control%itmax,                    &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 5 ), control%Lanczos_itmax,            &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 18 ), control%extra_vectors,           &
                                  control%error )

!  Set real values

      CALL SPECFILE_assign_value( spec( 6 ), control%stop_relative,            &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 7 ), control%stop_absolute,            &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 8 ), control%fraction_opt,             &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 9 ), control%f_0,                      &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 17 ), control%rminvr_zero,             &
                                  control%error )
     
!  Set logical values

      CALL SPECFILE_assign_value( spec( 10 ), control%unitm,                   &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 11 ), control%steihaug_toint,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 12 ), control%boundary,                &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 13 ), control%equality_problem,        &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 14 ), control%space_critical,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 15 ),                                  &
                                  control%deallocate_error_fatal,              &
                                  control%error )
!  Set charcter values

!     CALL SPECFILE_assign_value( spec( 16 ), control%prefix,                  &
!                                 control%error )

      RETURN

      END SUBROUTINE GLTR_read_specfile

!-*-*-*-*-*-*-*-*-*-*  G L T R _ S O L V E   S U B R O U T I N E  -*-*-*-*-*-*-*

      SUBROUTINE GLTR_solve( n, radius, f, X, R, VECTOR, data, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!  =========
!
!   n        number of unknowns
!   radius   trust-region radius
!   f        the value of the quadratic function at the current point.
!            Need not be set on entry. On exit it will contain the value at 
!            the best point found
!   X        the vector of unknowns. Need not be set on entry.
!            On exit, the best value found so far
!   R        the residual vector H x + c. On entry this must contain c
!   VECTOR   see inform%status = 2-3 and 6-7
!   data     private internal data
!   control  a structure containing control information. See GLTR_initialize
!   inform   a structure containing information. Components are
!             status     the input/output status. This must be set to 1 on
!                        initial entry or 4 on a re-entry when only radius has
!                        been reduced since the last entry. Other values are
!               2,6 on exit => the inverse of M must be applied to 
!                 VECTOR with the result returned in VECTOR and the subroutine 
!                 re-entered. This will only happen if unitm is .TRUE.
!               3,7 on exit => the product H * VECTOR must be formed, with
!                 the result returned in VECTOR and the subroutine re-entered
!               4 The iteration is to be restarted with a smaller radius but
!                 with all other data unchanged. Set R to c for this entry.
!               5 The iteration will be restarted. Reset R to c and re-enter.
!                 This exit will only occur if control%steihaug_toint is 
!                 .FALSE. and the solution lies on the trust-region boundary
!               0 the solution has been found
!              -1 an array allocation has failed
!              -2 an array deallocation has failed
!              -3 n and/or radius is not positive
!              -15 the matrix M appears to be indefinite
!              -18 the iteration limit has been exceeded
!              -36 the trust-region has been encountered in Steihaug-Toint mode
!             iter        the number of iterations performed
!             iter_pass2  the number of additional iterations performed to find
!                         the solution when it lies on the trust-region boundary
!             multiplier  the Lagrange multiplier corresponding to the 
!                         trust-region constraint
!             mnormx      the M-norm of x
!             negative_curvature  .TRUE. if negative curvature is found
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ) :: radius
      REAL ( KIND = wp ), INTENT( INOUT ) :: f
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: X, R, VECTOR
      TYPE ( GLTR_data_type ), INTENT( INOUT ) :: data
      TYPE ( GLTR_control_type ), INTENT( IN ) :: control        
      TYPE ( GLTR_info_type ), INTENT( INOUT ) :: inform

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: dim_sub, it, itp1, nroots
      REAL ( KIND = wp ) :: alpha, f_tol, other_root, xmx_trial
      CHARACTER ( LEN = 80 ) :: array_name

!  prefix for all output 

      CHARACTER ( LEN = LEN( TRIM( control%prefix ) ) - 2 ) :: prefix
      prefix = control%prefix( 2 : LEN( TRIM( control%prefix ) ) - 1 )

!  Branch to different sections of the code depending on input status

      SELECT CASE ( inform%status )
      CASE ( 1 ) 
        GO TO 100
      CASE ( 2 )  
        GO TO 200
      CASE ( 3 )  
        GO TO 300
      CASE ( 4 )  
        GO TO 400
      CASE ( 5 )  
        GO TO 500
      CASE ( 6 )  
        GO TO 600
      CASE ( 7 )  
        GO TO 700
      END SELECT

!  On initial entry, set constants

  100 CONTINUE 

!  Check for obvious errors

      IF ( n <= 0 ) GO TO 940
      IF ( radius <= zero ) GO TO 950

      data%iter = 0
      data%itmax = control%itmax ; IF ( data%itmax < 0 ) data%itmax = n
      data%Lanczos_itmax = control%Lanczos_itmax
      data%switch_to_Lanczos = .FALSE.
      data%prev_steihaug_toint = control%steihaug_toint
      IF ( data%Lanczos_itmax < 0 ) data%Lanczos_itmax = n
      data%printi = control%out > 0 .AND. control%print_level >= 1
      data%printd = control%out > 0 .AND. control%print_level >= 2
      inform%alloc_status = 0 ; inform%bad_alloc = ''
      inform%iter_pass2 = 0 ; inform%negative_curvature = .FALSE.
      inform%curv = HUGE( one ) ; inform%piv = HUGE( one )
      inform%rayleigh = HUGE( one )
      inform%mnormx = radius ; inform%multiplier = zero
      IF ( control%equality_problem ) THEN
        data%interior = .FALSE. ; data%switch = 0
      ELSE
        data%interior = .TRUE.
      END IF
      data%use_old = .FALSE. ; data%try_warm = .FALSE.
      data%radius2 = radius * radius ; data%rtol = roots_tol
      X = zero ; f = control%f_0 ; data%old_leftmost = zero

!  =====================
!  Array (re)allocations
!  =====================

!  Allocate P

      array_name = 'gltr: P'
      CALL SPACE_resize_array( n, data%P,                                      &
          inform%status, inform%alloc_status, array_name = array_name,         &
          deallocate_error_fatal = control%deallocate_error_fatal,             &
          exact_size = control%space_critical,                                 &
          bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 960

      IF ( .NOT. control%steihaug_toint ) THEN
        data%titmax = 100

!  Allocate space for the Lanczos tridiagonal

        array_name = 'gltr: D'
        CALL SPACE_resize_array( 0, data%itmax + 1, data%D,                    &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        array_name = 'gltr: OFFD'
        CALL SPACE_resize_array( data%itmax + 1, data%OFFD,                    &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  Allocate space for the factors of the Lanczos tridiagonal

        array_name = 'gltr: D_fact'
        CALL SPACE_resize_array( 0, data%itmax + 1, data%D_fact,               &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        array_name = 'gltr: OFFD_fact'
        CALL SPACE_resize_array( data%itmax + 1, data%OFFD_fact,               &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  Allocate space for the RHS and solution for the Lanczos subproblem

        array_name = 'gltr: C_sub'
        CALL SPACE_resize_array( 0, data%itmax + 1, data%C_sub,                &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960
        data%C_sub( 0 ) = one ; data%C_sub( 1 : data%itmax + 1 ) = zero

        array_name = 'gltr: X_sub'
        CALL SPACE_resize_array( 0, data%itmax + 1, data%X_sub,                &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  Allocate space for the history of Lagrange multipliers

        array_name = 'gltr: LAMBDA'
        CALL SPACE_resize_array( 0, data%itmax + 1, data%LAMBDA,               &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960
        data%LAMBDA( 0 ) = zero

!  Allocate space for workspace associated with the Lanczos subproblem

        array_name = 'gltr: Z'
        CALL SPACE_resize_array( 0, data%itmax + 1, data%Z,                    &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        array_name = 'gltr: W'
        CALL SPACE_resize_array( 0, data%itmax + 1, data%W,                    &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  Allocate space to store the alphas and rminvr for more efficient processing
!  in the second pass

        array_name = 'gltr: ALPHAS'
        CALL SPACE_resize_array( data%itmax + 1, data%ALPHAS,                  &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        array_name = 'gltr: RMINVRS'
        CALL SPACE_resize_array( data%itmax + 1, data%RMINVRS,                 &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  Allocate workspace for the sequence of smallest function values

        array_name = 'gltr: MIN_f'
        CALL SPACE_resize_array( 0, data%itmax + 1, data%MIN_f,                &
            inform%status, inform%alloc_status, array_name = array_name,       &
            deallocate_error_fatal = control%deallocate_error_fatal,           &
            exact_size = control%space_critical,                               &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  If required allocate extra space to store R, P and V in case of a second pass

        data%extra_vectors = control%extra_vectors - 2
        IF ( data%extra_vectors > 0 ) THEN
          array_name = 'gltr: R_extra'
          CALL SPACE_resize_array( n, data%R_extra,                            &
              inform%status, inform%alloc_status, array_name = array_name,     &
              deallocate_error_fatal = control%deallocate_error_fatal,         &
              exact_size = control%space_critical,                             &
              bad_alloc = inform%bad_alloc, out = control%error )
          IF ( inform%status == 0 ) THEN
            array_name = 'gltr: P_extra'
            CALL SPACE_resize_array( n, data%P_extra,                          &
                inform%status, inform%alloc_status, array_name = array_name,   &
                deallocate_error_fatal = control%deallocate_error_fatal,       &
                exact_size = control%space_critical,                           &
                bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status == 0 ) THEN
              array_name = 'gltr: V_extra'
              CALL SPACE_resize_array( n, data%extra_vectors, data%V_extra,    &
                  inform%status, inform%alloc_status, array_name = array_name, &
                  deallocate_error_fatal = control%deallocate_error_fatal,     &
                  exact_size = control%space_critical,                         &
                  bad_alloc = inform%bad_alloc, out = control%error )
              IF ( inform%status /= 0 ) THEN
                array_name = 'gltr: R_extra'
                CALL SPACE_dealloc_array( data%R_extra,                        &
                  inform%status, inform%alloc_status, array_name = array_name, &
                  bad_alloc = inform%bad_alloc, out = control%error )
                IF ( inform%status /= 0 ) GO TO 960
                array_name = 'gltr: P_extra'
                CALL SPACE_dealloc_array( data%P_extra,                        &
                  inform%status, inform%alloc_status, array_name = array_name, &
                  bad_alloc = inform%bad_alloc, out = control%error )
                IF ( inform%status /= 0 ) GO TO 960
                data%extra_vectors = 0
              END IF
            ELSE
              array_name = 'gltr: R_extra'
              CALL SPACE_dealloc_array( data%R_extra,                          &
                inform%status, inform%alloc_status, array_name = array_name,   &
                bad_alloc = inform%bad_alloc, out = control%error )
              IF ( inform%status /= 0 ) GO TO 960
              data%extra_vectors = 0
            END IF
          ELSE
            data%extra_vectors = 0
          END IF
        END IF
      ELSE
        data%extra_vectors = 0
      END IF
      data%save_vectors = data%extra_vectors > 0

!  ===========================
!  Start of the main iteration
!  ===========================

  190 CONTINUE

!  ----------------------------------
!  Obtain the preconditioned residual
!  ----------------------------------

      VECTOR = R
      IF ( .NOT. control%unitm ) THEN
        inform%status = 2 ; inform%iter = data%iter
        RETURN
      END IF

!  Obtain the scaled norm of the residual

  200 CONTINUE

      data%rminvr = DOT_PRODUCT( R, VECTOR )
      IF ( ABS( data%rminvr ) < control%rminvr_zero ) data%rminvr = zero
      IF ( data%rminvr < zero ) GO TO 930
      data%pgnorm = SIGN( SQRT( ABS( data%rminvr ) ), data%rminvr ) 

!  If the user has asked to save vectors, save VECTOR, R and P

      IF ( data%save_vectors ) THEN
        IF ( inform%iter < data%extra_vectors )                                &
          data%V_extra( : n , inform%iter + 1 ) = VECTOR
        IF ( inform%iter == data%extra_vectors ) THEN
          data%R_extra( : n ) = R
          data%P_extra( : n ) = data%P( : n )
        END IF
      END IF

      IF ( data%iter > 0 ) THEN
        data%beta = data%rminvr / data%rminvr_old
        data%diag = data%beta / data%alpha
!       data%offdiag = - SQRT( data%beta ) / data%alpha
        data%offdiag = SQRT( data%beta ) / ABS( data%alpha )
      ELSE

!  Compute the stopping tolerance

        data%diag = zero ; data%xmx = zero
        data%stop = MAX( control%stop_relative * data%pgnorm,                  &
                         control%stop_absolute )
        IF ( data%printi )                                                     &
          WRITE( control%out, "( /, A, ' stopping tolerence = ',               &
         &         ES12.4, ', radius = ', ES12.4 )" ) prefix, data%stop, radius
        IF ( .NOT. control%steihaug_toint ) data%C_sub( 0 ) = data%pgnorm
      END IF

!  Print details of the latest iteration

      IF ( data%printi ) THEN
        IF ( MOD( data%iter, 25 ) == 0 .OR. data%printd ) THEN
          IF ( data%interior ) THEN
            WRITE( control%out, "( /, A, '   Iter        f        ',          &
           &       ' pgnorm    step    norm p   norm x     curv ' )" ) prefix
          ELSE
            WRITE( control%out, 2000 ) prefix
          END IF
        END IF 

        IF ( data%interior ) THEN
          IF ( data%iter /= 0 ) THEN
            WRITE( control%out, "( A, I7, ES16.8, 4ES9.2, ES10.2 )" )          &
                   prefix, data%iter, f, data%pgnorm, data%alpha,              &
                   data%normp, SQRT( data%xmx ), inform%rayleigh
          ELSE
            WRITE( control%out, "( A, I7, ES16.8, ES9.2, 4X, '-',              &
           &      2( 8X, '-' ), 9X, '-' )" ) prefix, data%iter, f, data%pgnorm
          END IF
        ELSE
          IF ( data%iter /= 0 ) THEN
            IF ( data%printd ) WRITE( control%out, 2000 ) prefix
            WRITE( control%out, "( A, I7, ES16.8, ES9.2, ES15.8, 2I5 )" )      &
                   prefix, data%iter, data%MIN_f( data%itm1 ) + control%f_0,   &
                   ABS( data%x_last * data%offdiag ),                          &
                   data%LAMBDA( data%itm1 ), data%titer, data%tinfo
          END IF
        END IF
      END IF

!  Test for an interior approximate solution

      IF ( data%interior .AND. data%pgnorm <= data%stop ) THEN
        IF ( data%printi ) WRITE( control%out,                                 &
          "( A, ' pgnorm ', ES12.4, ' < ', ES12.4 )" )                         &
             prefix, data%pgnorm, data%stop
        inform%iter = data%iter ; inform%mnormx = SQRT( data%xmx )
        data%dim_sub = data%iter
        IF ( .NOT. control%steihaug_toint .AND. data%dim_sub > 0 )             &
          data%LAMBDA( data%dim_sub - 1 ) = zero
        GO TO 900
      END IF

      IF ( data%iter > 0 ) THEN

!  Test to see that iteration limit has not been exceeded

         IF ( data%iter >= data%itmax .AND. data%interior ) THEN
           inform%mnormx = SQRT( data%xmx )
           data%dim_sub = data%iter
           GO TO 910
         END IF

!  Obtain the search direction P

         data%xmp = data%beta * ( data%xmp + data%alpha * data%pmp )
         data%pmp = data%rminvr + data%pmp * data%beta * data%beta
         data%P( : n ) = - VECTOR + data%beta * data%P( : n )

!  If required, continue accumulating the Lanczos tridiagonal

         IF ( .NOT. control%steihaug_toint ) THEN
           data%D( data%iter ) = data%diag
           data%OFFD( data%iter ) = data%offdiag
         END IF

!  Check for convergence on the trust-region boundary

         IF ( data%iter >= data%itmax .OR. ( .NOT. data%interior .AND.         &
              ABS( data%offdiag * data%x_last ) <= data%stop ) ) THEN
   
!  Convergence on the trust-region boundary has occured. Determine at which
!  iteration a fraction, fraction_opt, of the optimal solution was found

           IF ( control%fraction_opt < one ) THEN
             f_tol = data%MIN_f( data%itm1 ) * control%fraction_opt
             DO dim_sub = 1, data%iter
                IF ( data%MIN_f( dim_sub - 1 ) <= f_tol ) EXIT
             END DO
           ELSE
             dim_sub = data%iter
           END IF
           data%dim_sub = dim_sub
   
!  Special case: the required fraction of f was achieved by an interior point

           IF ( dim_sub <= data%switch ) THEN
             inform%mnormx = SQRT( data%xmx )
             GO TO 900
           END IF

!          IF ( data%printi ) WRITE( control%out, 2020 )                       &
!            WRITE( control%out, 2020 ) data%MIN_f( dim_sub - 1 ), dim_sub, iter

!  Restore the solution to the Lanczos TR subproblem for this iteration

           data%use_old = .FALSE.
           IF ( dim_sub > 1 + data%switch ) THEN
             data%LAMBDA( dim_sub - 1 ) = data%LAMBDA( dim_sub - 2 )
           ELSE
             data%LAMBDA( dim_sub - 1 ) = zero
             data%try_warm = .FALSE.
           END IF

           CALL GLTR_ttrs( dim_sub, data%D( : dim_sub - 1 ),                   &
                      data%OFFD( : dim_sub - 1 ),                              &
                      data%D_fact( : dim_sub - 1 ),                            &
                      data%OFFD_fact( : dim_sub - 1 ),                         &
                      data%C_sub( : dim_sub - 1 ), radius, data%rtol,          &
                      data%interior, control%equality_problem,                 & 
                      data%titmax, data%try_warm, data%use_old,                &
                      data%old_leftmost, data%LAMBDA( dim_sub - 1 ),           &
                      data%MIN_f( dim_sub - 1 ),                               &
                      data%X_sub( : dim_sub - 1 ), data%tinfo, data%titer,     &
                      data%Z( : dim_sub - 1 ), data%W( : dim_sub - 1 ),        &
                      data%seed, data%printd, control%out, prefix )
  
!          IF ( data%printi ) WRITE( control%out, 2020 )                       &
!          WRITE( control%out, 2020 ) data%MIN_f(dim_sub-1), dim_sub, data%iter

!  Record the optimal objective function value and prepare to recover the
!  approximate solution

          f = data%MIN_f( dim_sub - 1 ) + control%f_0
          inform%multiplier = data%LAMBDA( dim_sub - 1 )
          data%tau = one
          inform%iter = data%iter ; data%iter = 0
          IF ( data%save_vectors ) GO TO 390
          inform%status = 5
          RETURN
        END IF
      ELSE

!  Special case for the first iteration

        data%P( : n ) = - VECTOR
        data%pmp = data%rminvr ; data%xmp = zero
        IF ( .NOT. control%steihaug_toint ) data%D( 0 ) = data%diag
      END IF

      IF ( .NOT. control%steihaug_toint )                                      &
        data%RMINVRS( data%iter + 1 ) = data%rminvr
      data%rminvr_old = data%rminvr

!  Compute the 2-norm of the search direction

      data%normp = TWO_NORM( data%P( : n ) )

!  Test for convergence

      IF ( data%interior .AND. data%normp <= data%rtol ) THEN
        IF ( data%printi ) WRITE( control%out,                                 &
             "( A, ' pnorm ', ES12.4, ' < ', ES12.4 )" )                       &
               prefix, data%normp, data%rtol
        inform%iter = data%iter ; inform%mnormx = SQRT( data%xmx )
        data%dim_sub = data%iter
        GO TO 900
      END IF

      data%itm1 = data%iter ; data%iter = data%iter + 1

!  ------------------------------ 
!  Obtain the product of H with p
!  ------------------------------ 

      VECTOR = data%P( : n ) ; inform%status = 3  ; inform%iter = data%iter
      RETURN

!  Obtain the curvature

  300 CONTINUE 

      inform%curv = DOT_PRODUCT( VECTOR, data%P( : n ) )
      inform%rayleigh = inform%curv / data%pmp

!  Obtain the stepsize and the new diagonal of the Lanczos tridiagonal

      IF ( inform%curv /= zero ) THEN
        data%alpha = data%rminvr / inform%curv
        data%diag = data%diag + one / data%alpha
      ELSE
        data%alpha = HUGE( one ) ** 0.25
        data%switch_to_Lanczos = .TRUE.
      END IF

!  Check that the Lanczos tridiagonal is still positive definite

      IF ( .NOT. inform%negative_curvature ) THEN
        IF ( data%iter > 1 ) THEN
          inform%piv = data%diag - ( data%offdiag / inform%piv ) * data%offdiag
        ELSE
          inform%piv = data%diag
        END IF
        inform%negative_curvature = inform%piv <= zero
      END IF

!  The matrix is indefinite

      IF ( data%interior .AND. inform%negative_curvature ) THEN

!  Find the appropriate point on the boundary

        CALL ROOTS_quadratic(                                                  &
                  data%xmx - data%radius2, two * data%xmp, data%pmp,           &
                  roots_tol, nroots, other_root, alpha )
        data%xmx = data%xmx + alpha * ( two * data%xmp + alpha * data%pmp )
        X = X + alpha * data%P
        f = f + alpha * ( half * alpha * inform%curv - data%rminvr )

!  If the Steihaug-Toint strategy is to be used, find the appropriate
!  point on the boundary and stop

        IF ( control%steihaug_toint .OR. data%switch_to_Lanczos ) THEN
          data%alpha = alpha
          data%dim_sub = data%iter
          IF ( .NOT. control%steihaug_toint )                                 &
            data%LAMBDA( data%dim_sub - 1 ) = zero
          GO TO 920

!  If a more accurate solution is required, switch modes

        ELSE
          data%interior = .FALSE.
          data%itmax = min( data%itmax, data%iter + data%Lanczos_itmax )
          data%switch = data%itm1
          IF ( data%printi ) THEN
             WRITE( control%out, 2010 ) prefix
             IF ( .NOT. data%printd ) WRITE( control%out, 2000 ) prefix
          END IF
        END IF
      END IF

!  If the current estimate of the solution is interior, see if the new point
!  is also interior

      IF ( data%interior ) THEN
         
        xmx_trial = data%xmx + data%alpha * ( data%xmp + data%xmp +            &
                                              data%alpha * data%pmp )

!  The new point is interior

        IF (  xmx_trial <= data%radius2 ) THEN
          data%xmx = xmx_trial
          X = X + data%alpha * data%P( : n )
          f = f - half * data%alpha * data%alpha * inform%curv
          IF ( .NOT. control%steihaug_toint )                                 &
            data%MIN_f( data%itm1 ) = f - control%f_0

!  The new point is outside the trust region

        ELSE

!  Find the appropriate point on the boundary

           CALL ROOTS_quadratic(                                              &
                     data%xmx - data%radius2, two * data%xmp, data%pmp,       &
                     roots_tol, nroots, other_root, alpha )
           data%xmx = data%xmx + alpha * ( two * data%xmp + alpha * data%pmp )
           X = X + alpha * data%P( : n )
           f = f + alpha * ( half * alpha * inform%curv - data%rminvr )

!  If the Steihaug-Toint strategy is to be used, find the appropriate
!  point on the boundary and stop

           IF ( control%steihaug_toint ) THEN
             data%alpha = alpha
             data%dim_sub = data%iter
             GO TO 920

!  If a more accurate solution is required, switch modes

           ELSE
             data%interior = .FALSE.
             data%itmax = min( data%itmax, data%iter + data%Lanczos_itmax )
             data%switch = data%itm1
             IF ( data%printi ) THEN
               WRITE( control%out, 2010 ) prefix
               IF ( .NOT. data%printd ) WRITE( control%out, 2000 ) prefix
             END IF
           END IF
        END IF
      END IF

!  Complete the new diagonal of the Lanczos tridiagonal matrix

      IF ( .NOT. control%steihaug_toint ) THEN
         data%D( data%itm1 ) = data%diag ; data%ALPHAS( data%iter ) = data%alpha
      END IF

      IF ( .NOT. ( control%steihaug_toint .OR. data%interior ) ) THEN

!  Solve the subproblem

        CALL GLTR_ttrs( data%iter, data%D( : data%itm1 ),                      &
             data%OFFD( : data%itm1 ), data%D_fact( : data%itm1 ),             &
             data%OFFD_fact( : data%itm1 ), data%C_sub( : data%itm1 ),         &
             radius, data%rtol, data%interior, control%equality_problem,       &
             data%titmax, data%try_warm,                                       &
             data%use_old, data%old_leftmost, data%LAMBDA( data%itm1 ),        &
             data%MIN_f( data%itm1 ), data%X_sub( : data%itm1 ),               &
             data%tinfo, data%titer, data%Z( : data%itm1 ),                    &
             data%W( : data%itm1 ), data%seed, data%printd, control%out,       &
             prefix )
        data%try_warm = .TRUE.

        data%LAMBDA( data%iter ) = data%LAMBDA( data%itm1 )
        data%use_old = data%old_leftmost < zero
        data%x_last = data%X_sub( data%itm1 )
 
      END IF
         
!  Update the residual

      R = R + data%alpha * VECTOR

!  =========================
!  End of the main iteration
!  =========================

      GO TO 190

!  ===================================
!  Use saved vectors to start 2nd pass
!  ===================================

  390 CONTINUE

      inform%iter_pass2 = MIN( data%dim_sub, data%extra_vectors )
      DO it = 0, inform%iter_pass2 - 1
        itp1 = it + 1
        data%X_sub( it ) = data%tau *                                          &
          ( data%X_sub( it ) / SQRT( data%RMINVRS( itp1 ) ) )
        data%tau = - SIGN( one, data%ALPHAS( itp1 ) ) * data%tau
      END DO

!  Update the solution estimate using the saved vectors

      X = MATMUL( data%V_extra( : n, : inform%iter_pass2 ),                    &
                  data%X_sub( 0 : inform%iter_pass2 - 1 ) )

      IF ( inform%iter_pass2 == data%dim_sub ) GO TO 900
      R = data%R_extra( : n )
      data%P( : n ) = data%P_extra( : n )
      data%iter = inform%iter_pass2
      data%rminvr_old = data%RMINVRS( data%iter )

      GO TO 500

!  ======================================================
!  Re-entry for solution with smaller trust-region radius
!  ======================================================

  400 CONTINUE 
      IF ( control%steihaug_toint ) GO TO 100
      IF ( data%prev_steihaug_toint ) GO TO 100

      X = zero ; f = control%f_0
      inform%iter = 0 ; inform%iter_pass2 = 0
      data%interior = .NOT. control%boundary
      data%use_old = .FALSE. ; data%try_warm = .TRUE.

!  Find the solution to the Lanczos TR subproblem with this radius

      CALL GLTR_ttrs( data%dim_sub, data%D( : data%dim_sub - 1 ),              &
                 data%OFFD( : data%dim_sub - 1 ),                              &
                 data%D_fact( : data%dim_sub - 1 ),                            &
                 data%OFFD_fact( : data%dim_sub - 1 ),                         &
                 data%C_sub( : data%dim_sub - 1 ), radius, data%rtol,          &
                 data%interior, control%equality_problem,                      &
                 data%titmax, data%try_warm, data%use_old,                     &
                 data%old_leftmost, data%LAMBDA( data%dim_sub - 1 ),           &
                 data%MIN_f( data%dim_sub - 1 ),                               &
                 data%X_sub( : data%dim_sub - 1 ), data%tinfo, data%titer,     &
                 data%Z( : data%dim_sub - 1 ), data%W( : data%dim_sub - 1 ),   &
                 data%seed, data%printd, control%out, prefix )      
      
!     IF ( data%printi ) WRITE( control%out, 2020 )                            &
!        WRITE( control%out, 2020 ) data%MIN_f( data%dim_sub - 1 ),            &
!                                   data%dim_sub, data%iter

!  Record the optimal objective function value and prepare to recover the
!  approximate minimizer

      inform%mnormx = radius 
      inform%multiplier = data%LAMBDA( data%dim_sub - 1 )
      f = data%MIN_f( data%dim_sub - 1 ) + control%f_0 
      data%iter = 0 ; data%tau = one

!  --------------------------------------
!  Special part of the code to obtain the 
!  approximate minimizer
!  --------------------------------------

  500 CONTINUE

!  ----------------------------------
!  Obtain the preconditioned residual
!  ----------------------------------

      VECTOR = R
      IF ( .NOT. control%unitm ) THEN
        inform%status = 6 ; inform%iter_pass2 = data%iter
        RETURN
      END IF

!  Obtain the scaled norm of the residual

  600 CONTINUE 

      itp1 = data%iter + 1

!  Update the solution estimate

      data%rminvr = data%RMINVRS( itp1 )
      IF ( data%iter /= 0 ) THEN
        X = X + data%tau                                                       &
            * ( data%X_sub( data%iter ) / SQRT( data%rminvr ) ) * VECTOR
      ELSE
        X = data%tau                                                           &
            * ( data%X_sub( data%iter ) / SQRT( data%rminvr ) ) * VECTOR
      END IF

!  If the approximate minimizer is complete, exit

      IF ( itp1 == data%dim_sub ) THEN
        inform%iter_pass2 = data%iter
        IF ( data%iter >= data%itmax ) GO TO 910
        GO TO 900
      END IF

      IF ( data%iter > 0 ) THEN
        data%beta = data%rminvr / data%rminvr_old
        data%P( : n ) = - VECTOR + data%beta * data%P( : n )
      ELSE

!  Special case for the first iteration

        data%P( : n ) = - VECTOR
      END IF
      data%rminvr_old = data%rminvr ; data%iter = itp1

!  ------------------------------ 
!  Obtain the product of H with p
!  ------------------------------ 

      VECTOR = data%P( : n ) 
      inform%status = 7  ; inform%iter_pass2 = data%iter
      RETURN

!  Obtain the curvature

  700 CONTINUE 

!  Retreive the stepsize

      data%alpha = data%ALPHAS( data%iter )
         
!  Update the residual

      R = R + data%alpha * VECTOR

      data%tau = - SIGN( one, data%alpha ) * data%tau

      GO TO 500

!  ===============
!  Exit conditions
!  ===============

!  Successful returns

  900 CONTINUE
      inform%status = 0
      RETURN

!  Too many iterations

  910 CONTINUE
      IF ( data%printi )                                                       &
        WRITE( control%out, "( /, A, ' Iteration limit exceeded ' ) " ) prefix
      inform%status = GALAHAD_error_max_iterations ; inform%iter = data%iter
      RETURN

!  Boundary encountered in Steihaug-Toint method

  920 CONTINUE

!  Find the gradient at the appropriate point on the boundary

      R = R + data%alpha * VECTOR

      IF ( data%printi ) WRITE( control%out, "( A, I7, ES16.8, 4X, '-', 4X,    &
     &     3ES9.2, ES10.2, //, A,                                              &
     &     ' Now leaving trust region (Steihaug-Toint)' )" )                   &
          prefix, data%iter, f, data%alpha, data%normp, SQRT( data%xmx ),      &
          inform%rayleigh, prefix
      inform%status = GALAHAD_warning_on_boundary ; inform%iter = data%iter
      RETURN

!  Unsuccessful returns

  930 CONTINUE
      IF ( control%error > 0 .AND. control%print_level > 0 )                   &
        WRITE( control%error,                                                  &
         "( A, ' The matrix M appears to be indefinite. Inner product = ',     &
     &      ES12.4  )" ) prefix, data%rminvr
      inform%status = GALAHAD_error_preconditioner ; inform%iter = data%iter
      RETURN

  940 CONTINUE
      IF ( control%error > 0 .AND. control%print_level > 0 )                   &
        WRITE( control%error,                                                  &
         "( A, ' n = ', I6, ' is not positive ' )" ) prefix, n
      inform%status = GALAHAD_error_restrictions ; inform%iter = data%iter
      RETURN

  950 CONTINUE
      IF ( control%error > 0 .AND. control%print_level > 0 )                   &
        WRITE( control%error,                                                  &
         "( A, ' The radius ', ES12.4 , ' is not positive ' )" ) prefix, radius
      inform%status = GALAHAD_error_restrictions ; inform%iter = data%iter
      RETURN

!  Allocation or deallocation error

  960 CONTINUE
      inform%iter = data%iter
      RETURN

!  Non-executable statements

 2000 FORMAT( /, A, '  Iter        f         pgnorm      lambda    ',          &
                    ' tr it info' )
 2010 FORMAT( /, A, ' Boundary encountered. Switching to Lanczos mode ' )
!2020 FORMAT( /, ' MIN_f, it_exit, it_total ', ES22.14, 2I6 )

!  End of subroutine GLTR_solve

      END SUBROUTINE GLTR_solve

!-*-*-*-*-*-  G L T R _ T E R M I N A T E   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE GLTR_terminate( data, control, inform )

!  ..............................................
!  .                                            .
!  .  Deallocate arrays at end of GLTR_solve    .
!  .                                            .
!  ..............................................

!  Arguments:
!  =========
!
!   data    private internal data
!   control see Subroutine GLTR_initialize
!   inform    see Subroutine GLTR_solve

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      TYPE ( GLTR_data_type ), INTENT( INOUT ) :: data
      TYPE ( GLTR_control_type ), INTENT( IN ) :: control        
      TYPE ( GLTR_info_type ), INTENT( INOUT ) :: inform

!-----------------------------------------------
!   L o c a l   V a r i a b l e
!-----------------------------------------------

      CHARACTER ( LEN = 80 ) :: array_name

!  Deallocate all internal arrays

      array_name = 'gltr: P'
      CALL SPACE_dealloc_array( data%P,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: D'
      CALL SPACE_dealloc_array( data%D,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: OFFD'
      CALL SPACE_dealloc_array( data%OFFD,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: ALPHAS'
      CALL SPACE_dealloc_array( data%ALPHAS,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: RMINVRS'
      CALL SPACE_dealloc_array( data%RMINVRS,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: MIN_f'
      CALL SPACE_dealloc_array( data%MIN_f,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: LAMBDA'
      CALL SPACE_dealloc_array( data%LAMBDA,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: D_fact'
      CALL SPACE_dealloc_array( data%D_fact,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: OFFD_fact'
      CALL SPACE_dealloc_array( data%OFFD_fact,                                &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: C_sub'
      CALL SPACE_dealloc_array( data%C_sub,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: X_sub'
      CALL SPACE_dealloc_array( data%X_sub,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: Z'
      CALL SPACE_dealloc_array( data%Z,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'gltr: W'
      CALL SPACE_dealloc_array( data%W,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'glrm: R_extra'
      CALL SPACE_dealloc_array( data%R_extra,                                   &
         inform%status, inform%alloc_status, array_name = array_name,           &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'glrm: P_extra'
      CALL SPACE_dealloc_array( data%P_extra,                                   &
         inform%status, inform%alloc_status, array_name = array_name,           &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'glrm: V_extra'
      CALL SPACE_dealloc_array( data%V_extra,                                   &
         inform%status, inform%alloc_status, array_name = array_name,           &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      RETURN

!  End of subroutine GLTR_terminate

      END SUBROUTINE GLTR_terminate

!-*-*-*-*-*-*-*-*-*-*-  G L T R _ t t r s  S U B R O U T I N E -*-*-*-*-*-*-*-*

      SUBROUTINE GLTR_ttrs( n, D, OFFD, D_fact, OFFD_fact, C, radius, rtol,    &
                           interior, equality, itmax, try_warm, use_old,       &
                           old_leftmost, lambda, f, X, info, iter, Z, W,       &
                           seed, debug, out, prefix ) 

!  Subroutine GLTR_ttrs
!  ====================

!  Given an n by n symmetric tridiagonal matrix T, an n-vector c, and a
!  positive number radius, this subroutine determines a vector
!  x which approximately minimizes the quadratic function

!     f(x) = (1/2)*x'*T*x + c'*x

!  subject to the Euclidean norm constraint

!     norm(x) <= radius

!  (or if equality = .TRUE.   norm(x) = radius).

!  This subroutine computes an approximation x and a Lagrange
!  multiplier lambda such that either lambda is zero and

!      norm(x) <= (1+rtol)*radius,

!  or lambda is positive and

!      abs(norm(x) - radius) <= rtol*radius.

!  If xsol is the solution to the problem, the approximation x
!  satisfies

!      f(x) <= ((1 - rtol)**2)*f(xsol)

!  Dummy arguments
!  ===============

!    n is an integer variable.
!      On entry n is the order of T.
!      On exit n is unchanged.

!    D is a real array of dimension (n).
!      On entry D must contain the diagonal of T
!      Unchanged on exit.

!    OFFD is a real array of dimension (n-1).
!      On entry D must contain the subdiagonal of T
!      Unchanged on exit.

!    D_fact is a real array of dimension (n).
!      On entry D_fact need not be specified.
!      On exit D_fact contains the D part of the LDL(transpose) 
!      factorization of T + lambda I

!    OFFD_fact is a real array of dimension (n-1).
!      On entry OFFD_fact need not be specified.
!      On exit OFFD_fact contains the subdiagonal part of the L factor 
!      from the LDL(transpose) factorization of T + lambda I

!    C is an real array of dimension n.
!      On entry C specifies the linear term in the quadratic.
!      On exit C is unchanged.

!    radius is a real variable.
!      On entry radius is a bound on the Euclidean norm of x.
!      On exit radius is unchanged.

!    rtol is a real variable.
!      On entry rtol is the relative accuracy desired in the
!         solution. Convergence occurs if

!      abs(norm(x) - radius) <= rtol*radius.

!      On exit rtol is unchanged.

!    interior is a logical variable.
!      On entry, interior should be set to .TRUE. if an interior solution
!      is possible, and .FALSE. otherwise.

!      On exit it will be set to .TRUE. if an interior solution has been
!      found and to .FALSE. otherwise.

!    equality is a logical variable.
!      On entry, equality should be set to .TRUE. if a solution on the
!      trust-region boundary is required, and .FALSE. otherwise.

!    itmax is an integer variable.
!      On entry itmax specifies the maximum number of iterations.
!      On exit itmax is unchanged.

!    try_warm is a logical variable.
!      On entry try_warm is .TRUE. if the input value lambda is to be 
!       tried before any other estimate.
!      On exit use_old is unchanged.

!    use_old is a logical variable.
!      On entry use_old is .TRUE. if the leftmost eigenvalue of the leading
!       n-1 by n-1 block is given
!      On exit use_old is unchanged.

!    old_leftmost is a real variable.
!      On entry old_leftmost gives the leftmost eigenvalue of the leading
!       n-1 by n-1 block. Only required if use_old is .TRUE.
!      On exit gives the leftmost eigenvalue of T if T is indefinite.

!    lambda is a real variable.
!      On entry lambda is an initial estimate of the Lagrange
!         multiplier for the constraint norm(x) <= radius.
!      On exit lambda contains the final estimate of the multiplier.

!    f is a real variable.
!      On entry f need not be specified.
!      On exit f is set to f(x) at the output x.

!    x is a real array of dimension n.
!      On entry x need not be specified.
!      On exit x is set to the final estimate of the solution.

!    info is an integer variable.
!      On entry info need not be specified.
!      On exit info is set as follows:

!         info = 1  The function value f(x) has the relative
!                   accuracy specified by rtol.

!         info = 2  The Newton search direction is too small to make
!                   further progress

!         info = 4  Failure to converge after itmax iterations.
!                   On exit x is the best available approximation.

!    iter is an integer variable.
!      On entry iter need not be specified.
!      On exit iter gives the total number of iterations required.

!    Z is a real work array of dimension n.

!    W is a real work array of dimension n.

!    debug is a logical variable.
!      On entry debug should be .TRUE. if debug printing is required
!      On exit debug is unchanged.

!    out is an integer variable.
!      On entry the unit for output if required
!      On exit out is unchanged.

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n, itmax, out
      INTEGER, INTENT( OUT ) :: info, iter 
      LOGICAL, INTENT( IN ) :: debug, use_old, try_warm
      LOGICAL, INTENT( INOUT ) :: interior
      LOGICAL, INTENT( IN ) :: equality
      REAL ( KIND = wp ), INTENT( IN ) :: radius, rtol
      REAL ( KIND = wp ), INTENT( INOUT ) :: lambda, old_leftmost
      REAL ( KIND = wp ), INTENT( OUT ) :: f
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n - 1 ) :: OFFD
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: c, D
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n - 1 ) :: OFFD_fact
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: D_fact, X, Z, W
      TYPE ( RAND_seed ), INTENT( INOUT ) :: seed
      CHARACTER ( LEN = * ), INTENT( IN ) :: prefix

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: indef, i, it
      REAL ( KIND = wp ) :: alpha, lambda_pert, ztx, rxnorm2, distx, tol,      &
                            xnorm, leftmost, delta_lambda, macheps, pert_l

!  Initialization

      iter = 1
      macheps = EPSILON( one ) ; pert_l = macheps ** 0.5 ; tol = macheps ** 0.66

!  First, try a warm start
!  =======================

      IF ( try_warm ) THEN

!  Compute T + lambda*I

        OFFD_fact = OFFD
        D_fact = D + lambda

!  Find the Cholesky factors of T

        CALL PTTRF( n, D_fact, OFFD_fact, indef )
      
!  If T is positive definite, solve  T x = - c

        IF ( indef == 0 ) THEN 
           W = - C
           DO i = 1, n - 1
             W( i + 1 ) = W( i + 1 ) - OFFD_fact( i ) * W( i )
           END DO
           X = W / D_fact
           rxnorm2 = DOT_PRODUCT( W, X )
           DO i = n - 1, 1, - 1
             X( i ) = X( i ) - OFFD_fact( i ) * X( i + 1 )
           END DO

!  If the solution lies outside the trust-region, it provides a good initial
!  estimate of the solution to the TR problem

          xnorm = TWO_NORM( X ) 
          IF ( debug ) THEN
            WRITE( out, "( A, 8X, 'lambda', 13X, 'xnorm', 15X, 'radius' )" )  &
              prefix
            WRITE( out, "( A, 3ES20.12 )" ) prefix, lambda, xnorm, radius
          END IF
          IF ( ABS( xnorm - radius ) <= rtol * radius ) THEN
            f = - half * ( rxnorm2 + lambda * xnorm ** 2 ) 
            info = 1
            RETURN
          END IF

          IF ( xnorm > radius ) GO TO 20

        END IF
      END IF

!  If the warm start fails, check for an unconstrained solution
!  ============================================================

!  10 CONTINUE

      IF ( interior ) THEN

!  Attempt the Cholesky factorization of T

        OFFD_fact = OFFD ; D_fact = D
        CALL PTTRF( n, D_fact, OFFD_fact, indef )

!  If T is positive definite, solve  T x = - c

        IF ( indef == 0 ) THEN 
          W = - C
          DO i = 1, n - 1
            W( i + 1 ) = W( i + 1 ) - OFFD_fact( i ) * W( i )
          END DO
          X = W / D_fact
          DO i = n - 1, 1, - 1
            X( i ) = X( i ) - OFFD_fact( i ) * X( i + 1 )
          END DO

!  If the solution lies within the trust-region, it provides an interior 
!  solution to the TR problem

          xnorm = TWO_NORM( X ) 
          IF ( debug ) THEN
            WRITE( out, "( A, 8X, 'lambda', 13X, 'xnorm', 15X, 'radius' )" )  &
              prefix
            WRITE( out, "( A, 3ES20.12 )" ) prefix, zero, xnorm, radius
          END IF

!  Interior solution case
!  ======================

          IF ( xnorm <= radius ) THEN
            lambda = zero               
            f = half * DOT_PRODUCT( C, X )
            info = 1
            RETURN 

!  The Lagrange multiplier will be positive and lambda is in L

          ELSE
             lambda = zero
          END IF
        ELSE
          interior = .FALSE.
        END IF

      ELSE
        indef = 1
      END IF

!  The solution is not interior. Compute the leftmost eigenvalue
!  =============================================================

      IF ( indef > 0 ) THEN
        leftmost = GLTR_leftmost_eigenvalue( n, D, OFFD, tol, use_old,        &
                     old_leftmost, it, debug, out, prefix )

        IF ( .NOT. equality ) leftmost = MIN( leftmost, zero )
        
        IF ( debug ) WRITE( out, "( A, ' iteration ', I6,                     &
       &     ' leftmost eigenvalue = ', ES22.14 )") prefix, it, leftmost

        old_leftmost = leftmost ;  tol = point1 * tol

!  Special check to see if we are in the hard case

!       IF ( leftmost <= zero ) THEN
          IF ( leftmost <= zero ) THEN
            lambda_pert = - leftmost * ( one + pert_l ) + pert_l
          ELSE
            lambda_pert = - leftmost * ( one - pert_l ) + pert_l
          END IF

!  Compute T + lambda*I

          DO 
            OFFD_fact = OFFD
            D_fact = D + lambda_pert

!  Find the Cholesky factors of T

            CALL PTTRF( n, D_fact, OFFD_fact, indef )

!  Make sure that T really is numerically positive definite

            IF ( indef == 0 ) EXIT

!  H is still numerically indefinite and must be perturbed a bit more

            pert_l = pert_l + pert_l
            IF ( leftmost <= zero ) THEN
              lambda_pert = lambda_pert * ( one + pert_l ) + pert_l
            ELSE
              lambda_pert = lambda_pert * ( one - pert_l ) + pert_l
            END IF
          END DO

!  Solve T x = - c

          W = C
          DO i = 1, n - 1
             W( i + 1 ) = W( i + 1 ) - OFFD_fact( i ) * W( i )
          END DO
          X = - W / D_fact
          rxnorm2 = - DOT_PRODUCT( W, X )
          DO i = n - 1, 1, - 1
             X( i ) = X( i ) - OFFD_fact( i ) * X( i + 1 )
          END DO

          xnorm = TWO_NORM( X ) 
          IF ( debug ) THEN
            WRITE( out, "( A, 8X, 'lambda', 13X, 'xnorm', 15X, 'radius' )" )   &
              prefix
            WRITE( out, "( A, 3ES20.12 )" ) prefix, lambda_pert, xnorm, radius
          END IF

!  Hard case
!  =========

!  If the norm of X is smaller than radius, we are in the hard case

          IF ( xnorm < radius .AND. .NOT. equality ) THEN
            lambda = - leftmost

!  Compute a leftmost eigenvector

            CALL GLTR_leftmost_eigenvector( n, leftmost, D, OFFD, D_fact,      &
                                            OFFD_fact, Z, it, seed )
            IF ( debug ) WRITE( out, "( A, ' iteration ', I6,                  &
           &  ' leftmost eigenvector found ' )" ) prefix, it

!  Compute the step alpha so that X + alpha Z lies on the trust-region
!  boundary and gives the smaller value of q

            ztx = DOT_PRODUCT( Z, X ) / radius 
            distx = ( radius - xnorm ) * ( ( radius + xnorm ) / radius ) 
            alpha = sign( distx / ( abs( ztx ) +                            &
                          sqrt( ztx ** 2 + distx / radius ) ), ztx )

!  Record the optimal values

            X = X + alpha * Z
            f = - half * ( rxnorm2 + lambda * radius ** 2 ) 
            info = 1

            RETURN

!  The Lagrange multiplier will be positive and lambda is in L

          ELSE
             lambda = lambda_pert
          END IF
!       ELSE
!         interior = .TRUE.
!         GO TO 10
!       END IF
      ELSE
        old_leftmost = zero
        IF ( interior .AND. debug )                                            &
          WRITE( out, "( A, 8X,'lambda',13X,'xnorm',15X,'radius' )" ) prefix
      END IF

!  It is now simply a matter of applying Newton's method starting from lambda

!  Main Newton iteration
!  =====================

   20 CONTINUE

      DO iter = 2, itmax 

!  Compute the Newton correction

        W = X / xnorm
        DO i = 1, n - 1
          W( i + 1 ) = W( i + 1 ) - OFFD_fact( i ) * W( i )
        END DO

        delta_lambda = ( ( xnorm - radius ) / radius ) /                       &
                          DOT_PRODUCT( W, W / D_fact )
!       IF ( delta_lambda < zero ) THEN
!          WRITE( 6, * ) 'dlambda', delta_lambda, &
!                   DOT_PRODUCT( W, W / D_fact )
!          write(6,"( ( 2ES12.4 ) )" ) ( W( i ), D_fact( i ), i = 1, n )
!       END IF

!  Check that the Newton correction is significant

!       WRITE( out, "( ' lambda, delta ', 2ES22.14 )" ) lambda, delta_lambda

        IF ( ABS( delta_lambda ) < macheps * ABS( lambda ) ) THEN
          f = - half * ( rxnorm2 + lambda * xnorm ** 2 ) 
          info = 2
          RETURN
        END IF

!  Compute the new estimate of lambda

        lambda = lambda + delta_lambda

!  Find the Cholesky factorization of T + lambda*I

        OFFD_fact = OFFD ; D_fact = D + lambda
        CALL PTTRF( n, D_fact, OFFD_fact, indef )

!       IF ( indef /= 0 ) THEN
!         WRITE( out, * ) ' ... but T+lambdaI should be positive definite ...'
!         STOP
!       END IF

!  Solve the equation (T + lambda*I) x = - c

        W = C
        DO i = 1, n - 1
          W( i + 1 ) = W( i + 1 ) - OFFD_fact( i ) * W( i )
        END DO
        X = - W / D_fact
        rxnorm2 = - DOT_PRODUCT( W, X )
        DO i = n - 1, 1, - 1
          X( i ) = X( i ) - OFFD_fact( i ) * X( i + 1 )
        END DO
        xnorm = TWO_NORM( X ) 
        IF ( debug ) WRITE( out, "( A, 3ES20.12 )" )                          &
          prefix, lambda, xnorm, radius

!  Test for convergence

!       WRITE( out, * ) ' diff, tol ', xnorm - radius, rtol * radius
        IF ( xnorm - radius <= rtol * radius ) THEN
!         f = DOT_PRODUCT( C, X ) + half * DOT_PRODUCT( X, D * X )  &
!            + DOT_PRODUCT( X( : n - 1 ) * OFFD, X( 2 : ) )
!         WRITE( out, "( ' real f = ', ES22.14 )" ) f
          f = - half * ( rxnorm2 + lambda * xnorm ** 2 ) 
          info = 1
          RETURN
        END IF

      END DO

!  Test for termination

      info = 4 
      f = - half * ( rxnorm2 + lambda * xnorm ** 2 ) 
      RETURN

!  End of subroutine GLTR_ttrs

      END SUBROUTINE GLTR_ttrs

!-*-*-*-*-*-*-*-*-  GLTR_leftmost_eigenvalue F U N C T I O N  -*-*-*-*-*-*-*-*-*

      FUNCTION GLTR_leftmost_eigenvalue( n, D, OFFD, tol, use_old,             &
         old_leftmost, iter, debug, out, prefix  )

      REAL ( KIND = wp ) :: GLTR_leftmost_eigenvalue

!  Compute the leftmost eigenvalue, e_1, of a symmetric tridiagonal matrix

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n, out
      INTEGER, INTENT( OUT ) :: iter
      REAL ( KIND = wp ), INTENT( IN ) :: tol, old_leftmost
      LOGICAL, INTENT( IN ) :: debug, use_old
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: D
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n - 1 ) :: OFFD
      CHARACTER ( LEN = * ), INTENT( IN ) :: prefix

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, neg, nroots
      REAL ( KIND = wp ) :: l, u, e, piv, dpiv, e_trial, infinity
      REAL ( KIND = wp ) :: root1, root2, tol_interval, sum, prod
      REAL ( KIND = wp ), PARAMETER :: pert = ten ** ( - 6 )

!  Special case: n = 1

      iter = 0
      IF ( n == 1 ) THEN
        GLTR_leftmost_eigenvalue = D( 1 )
        RETURN
      END IF

!  Initialize lower and upper bounds and eigenvalue estimates using Gersgorin
!  bounds

      l = MIN( D( 1 ) - ABS( OFFD( 1 ) ) , D( n ) - ABS( OFFD( n - 1 ) ),      &
               MINVAL( D( 2 : n - 1 ) - ABS( OFFD( : n - 2 ) )                 &
                                      - ABS( OFFD( 2 : ) ) ) )
      u = MAX( D( 1 ) + ABS( OFFD( 1 ) ) , D( n ) + ABS( OFFD( n - 1 ) ),      &
               MAXVAL( D( 2 : n - 1 ) + ABS( OFFD( : n - 2 ) )                 &
                                      + ABS( OFFD( 2 : ) ) ) )

      infinity = one + u - l
      IF ( use_old ) THEN
        u = MIN( u, old_leftmost )
!       e = half * ( l + u )
        e = u - pert
      ELSE
        e = l
      END IF
!     tol_interval = tol * ( one + half * ABS( l + u ) )
      tol_interval = tol * ( one + half * ABS( l ) + ABS( u ) )

      IF ( debug ) THEN
        WRITE( out, "( A, '  it neg        l               e               u', &
       &               '              piv' )" ) prefix
!       WRITE( out, "( A, 2I4, 3ES16.8, '        -' )") prefix, iter, 0, l, e, u
      END IF
   
!  Main iteration loop

      ITER_LOOP : DO
        iter = iter + 1

!  Compute the inertia of T - e I by implicitly factoring the matrix

        neg = 0
        piv = D( 1 ) - e
        dpiv = - one

!  If a zero pivot is encountered, reset the upper bound

        IF ( piv == zero ) THEN
          u = e
          e = half * ( l + u )
          CYCLE ITER_LOOP

!  If a negative pivot is encountered, exit 

        ELSE IF ( piv < zero ) THEN
          neg = 1
        END IF

        DO i = 2, n

!  Update the pivot

          dpiv = - one + dpiv * ( OFFD( i - 1 ) / piv ) ** 2
          piv = ( D( i ) - ( OFFD( i - 1 ) ** 2 ) / piv ) - e

!  If a zero last pivot is encountered, return with e as the required eigenvalue

          IF ( piv == zero ) THEN
            IF ( neg == 0 .AND. i == n ) THEN
              GLTR_leftmost_eigenvalue = e
              RETURN

!  If a zero pivot is encountered, reset the upper bound

            ELSE
              u = e
              e = half * ( l + u )
              CYCLE ITER_LOOP
            END IF

!  If more than one negative pivot is encountered, exit 

          ELSE IF ( piv < zero ) THEN
            neg = neg + 1
            IF ( neg > 1 ) THEN
              piv = infinity
              dpiv = one
              EXIT
            END IF
          END IF

        END DO

        IF ( debug ) WRITE( out, "( A, 2I4, 4ES16.8 )" )                       & 
          prefix, iter, neg, l, e, u, piv

!  Increase the lower bound

        IF ( neg == 0 ) THEN
          l = e

!  Reduce the upper bound

        ELSE
          u = e
        END IF

!  Test for convergence

        IF ( ABS( piv ) < tol .OR. u - l < tol_interval ) THEN
          GLTR_leftmost_eigenvalue = e
          IF ( debug ) WRITE( out,                                             &
            "(/, A, ' leftmost eigenvalue = ', ES22.14 )")                     &
                              prefix, GLTR_leftmost_eigenvalue
          RETURN
        END IF

!  Compute the Newton step

        IF ( use_old ) THEN
          sum = two * e + piv + ( e - old_leftmost ) * dpiv
          prod = - ( e - old_leftmost ) * piv + e * sum - e * e
          CALL ROOTS_quadratic( prod, - sum, one, roots_tol, nroots,           &
                                root1, root2 )
           e_trial = root1
        ELSE 
           e_trial = e  - piv / dpiv
        END IF

!  If the estimate lies in the interval (l,e_2) and the Newton step 
!  continues lies in [l,u], use the Newton step as the next estimate
!  Otherwise bisect the bounds to get the new eigenvalue estimate

        e = GLTR_new_eig_est( neg, l, u, e_trial )

      END DO ITER_LOOP

      CONTAINS

        FUNCTION GLTR_new_eig_est( neg, l, u, e_trial )

!  Compute the next iterate. This function is PURELY here to stop an 
!  optimization "bug" with certain compilers!

        REAL ( KIND = wp ) :: GLTR_new_eig_est
        INTEGER, INTENT( IN ) :: neg
        REAL ( KIND = wp ), INTENT( IN ) :: l, u, e_trial

!  If the estimate lies in the interval (l,e_2) and the Newton step 
!  continues to lie in [l,u], use the Newton step as the next estimate

        IF ( neg <= 1 .AND. ( e_trial > l .AND. e_trial < u ) ) THEN
           GLTR_new_eig_est = e_trial
        ELSE

!  Otherwise bisect the bounds to get the new eigenvalue estimate

           GLTR_new_eig_est = half * ( l + u )
        END IF
     
        RETURN
        END FUNCTION GLTR_new_eig_est

!  End of subroutine GLTR_leftmost_eigenvalue

      END FUNCTION GLTR_leftmost_eigenvalue

!-*-*-*-*-*-*-*-  GLTR_leftmost_eigenvector S U B R O U T I N E  -*-*-*-*-*-*-*

      SUBROUTINE GLTR_leftmost_eigenvector( n, est, D, OFFD, D_fact,           &
                                            OFFD_fact, U, iter, seed )

!  Compute an eigenvector u corresponding to the leftmost eigenvalue of a 
!  symmetric tridiagonal matrix using inverse iteration

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( OUT ) :: iter
      REAL ( KIND = wp ), INTENT( IN ) :: est
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: D
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n - 1 ) :: OFFD
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: D_fact, U
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n - 1 ) :: OFFD_fact
      TYPE ( RAND_seed ), INTENT( INOUT ) :: seed

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, indef
      REAL ( KIND = wp ) :: e, wnorm, perturb
      REAL ( KIND = wp ), PARAMETER :: pert = ten ** ( - 6 )
      REAL ( KIND = wp ), PARAMETER :: conv = ten ** ( - 8 )

!  Factorize T - est I

      perturb = pert * ( one - est )
      DO 
        e = est - perturb
        D_fact = D - e
        OFFD_fact = OFFD

!  Attempt the Cholesky factorization of T

        CALL PTTRF( n, D_fact, OFFD_fact, indef )

        IF ( indef == 0 ) EXIT
        perturb = ten * perturb
      END DO

!  Initialize a random initial estimate of U

!     CALL RANDOM_NUMBER( U )
      DO i = 1, n
        CALL RAND_random_real( seed, .TRUE., U( i ) )
      END DO

!  Inverse iteration

      DO iter = 1, 5

!  Solve ( T - est I ) w = u, overwriting u with the solution

        DO i = 1, n - 1
          U( i + 1 ) = U( i + 1 ) - OFFD_fact( i ) * U( i )
        END DO

        U = U / D_fact

        DO i = n - 1, 1, - 1
          U( i ) = U( i ) - OFFD_fact( i ) * U( i + 1 )
        END DO

!  Normalize w

        wnorm = one / TWO_NORM( U )
        U = U * wnorm

        IF ( ABS( wnorm - perturb ) <= conv ) EXIT

      END DO
      RETURN

!  end of subroutine GLTR_leftmost_eigenvector

      END SUBROUTINE GLTR_leftmost_eigenvector

!-*-*-*-*-*-  End of G A L A H A D _ G L T R  double  M O D U L E  *-*-*-*-*-*-

   END MODULE GALAHAD_GLTR_double
