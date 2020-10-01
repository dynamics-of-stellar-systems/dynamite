! THIS VERSION: GALAHAD 2.2 - 27/05/2008 AT 16:45 GMT.

!-*-*-*-*-*-*-*-*-*-  G A L A H A D _ L S R T   M O D U L E  *-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.1, November 24th, 2007

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_LSRT_double

!      -------------------------------------------------------------
!      |                                                           |
!      | Solve the regularised least-squares problem               |
!      |                                                           |
!      |    minimize  1/2 || A x - b ||^2 + 1/p sigma || x ||^p    |
!      |                                                           |
!      | for given sigma and p >= 2 using a                        |
!      ! Lanczos bi-diagonalisation method                         |
!      |                                                           |
!      -------------------------------------------------------------

      USE GALAHAD_SYMBOLS
      USE GALAHAD_SPACE_double
      USE GALAHAD_ROOTS_double, ONLY : ROOTS_quadratic
      USE GALAHAD_SPECFILE_double
      USE GALAHAD_NORMS_double, ONLY: TWO_NORM
      USE GALAHAD_LSTR_double, ONLY :                                           &
        LSRT_transform_bidiagonal => LSTR_transform_bidiagonal,                 &
        LSRT_backsolve_bidiagonal => LSTR_backsolve_bidiagonal

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: LSRT_initialize, LSRT_read_specfile, LSRT_solve, LSRT_terminate

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
      REAL ( KIND = wp ), PARAMETER :: three = 3.0_wp
      REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp

      REAL ( KIND = wp ), PARAMETER :: roots_tol = ten ** ( - 12 )
      REAL ( KIND = wp ), PARAMETER :: error_tol = ten ** ( - 12 )
      REAL ( KIND = wp ), PARAMETER :: lambda_start = ten ** ( - 12 )

!--------------------------
!  Derived type definitions
!--------------------------

      TYPE, PUBLIC :: LSRT_control_type
        INTEGER :: error, out, print_level, itmax
        INTEGER :: extra_vectors, stopping_rule, freq
        REAL ( KIND = wp ) :: stop_relative, stop_absolute, fraction_opt
        LOGICAL :: space_critical, deallocate_error_fatal
        CHARACTER ( LEN = 30 ) :: prefix
      END TYPE

      TYPE, PUBLIC :: LSRT_inform_type
        INTEGER :: status, alloc_status, iter, iter_pass2
        INTEGER :: biters, biter_min, biter_max
        REAL ( KIND = wp ) :: obj, multiplier, x_norm, r_norm, Atr_norm
        REAL ( KIND = wp ) :: biter_mean
        CHARACTER ( LEN = 80 ) :: bad_alloc
      END TYPE

      TYPE, PUBLIC :: LSRT_data_type
        PRIVATE
        INTEGER :: iter, itmax, end_pass2, switch, bitmax, bstatus, biter
        INTEGER :: branch, freq, extra_vectors
        REAL ( KIND = wp ) :: obj_0, alpha_kp1, beta_kp1, beta_1, g_norm2
        REAL ( KIND = wp ) :: rho_k, rho_bar, phi_k, phi_bar, theta_kp1
        REAL ( KIND = wp ) :: s, c, s_w, c_w, zeta_bar, gamma, z_norm
        REAL ( KIND = wp ) :: eta_bar, lambda_km1, lambda_bar, ww, sigma2
        REAL ( KIND = wp ) :: error_tol, stop, decrease_st, omega
        LOGICAL :: printi, printd, header, use_old, try_warm, save_vectors
        LOGICAL :: one_pass
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: B_diag
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: B_offdiag
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DECREASE
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: LAMBDA
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: R_diag
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: R_offdiag
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X_sub
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: W
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: F
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: G
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: U_extra
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: V_extra
      END TYPE

!---------------------------------
!   I n t e r f a c e  B l o c k s
!---------------------------------

     INTERFACE ROTG

       SUBROUTINE SROTG( a, b, c, s )
       REAL :: a, b, c, s
       END SUBROUTINE SROTG

       SUBROUTINE DROTG( a, b, c, s )
       DOUBLE PRECISION :: a, b, c, s
       END SUBROUTINE DROTG

     END INTERFACE 

    CONTAINS

!-*-*-*-*-*-  L S R T _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE LSRT_initialize( data, control )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  .  Set initial values for the LSRT control parameters  .
!   
!  Argument:
!  =========
!
!   data     private internal data
!   control  a structure containing control information. Components are -
!            stopping_rule   the stopping rule used (see below) 
!               (0=1.0, 1=norm step, 2=norm step/sigma)
!            stop_relative      the iteration stops successfully when the 
!            stop_absolute      ||A^Tr|| is less than 
!                        max( stop_relative * ||A^Tr initial ||
!                             stop_absolute )
!            fraction_opt a direction which gives at least fraction_opt times 
!                        the optimal objective value will be found
!            space_critical  .TRUE. if allocated arrays are required to be 
!                            exactly the right size 
!            deallocate_error_fatal. .TRUE. if a deallocation error should
!                                    terminate execution
!            extra_vectors    the number of extra work vectors of length n used
!            out         output unit
!            print_level print level. > 0 for output > 1 for debug
!            itmax       maximum number of iterations allowed
!            freq        frequency for solving the bi-diagonal problem
!            prefix  output lines will be prefixed by 
!                        prefix(2:LEN(TRIM(%prefix))-1)
!               where prefix contains the required string enclosed in 
!               quotes, e.g. "string" or 'string'
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t
!-----------------------------------------------

      TYPE ( LSRT_DATA_TYPE ), INTENT( OUT ) :: data
      TYPE ( LSRT_CONTROL_TYPE ), INTENT( OUT ) :: control        

!  Set initial control parameter values

      control%stopping_rule = 1
      control%stop_relative = SQRT( EPSILON( one ) )
      control%stop_absolute = zero
      control%fraction_opt = one
      control%space_critical = .FALSE.
      control%deallocate_error_fatal  = .FALSE.
      control%error = 6
      control%out = 6
      control%print_level = 0
      control%extra_vectors = 0
      control%freq = 1
      control%itmax = - 1
      control%prefix = '""     '

      data%branch = 1

      RETURN

!  End of subroutine LSRT_initialize

      END SUBROUTINE LSRT_initialize

!-*-*-*-*-   L S R T _ R E A D _ S P E C F I L E  S U B R O U T I N E   -*-*-*-

      SUBROUTINE LSRT_read_specfile( control, device, alt_specname )

!  Reads the content of a specification file, and performs the assignment of 
!  values associated with given keywords to the corresponding control parameters

!  The default values as given by LSRT_initialize could (roughly) 
!  have been set as:

!  BEGIN LSRT SPECIFICATIONS (DEFAULT)
!   error-printout-device                           6
!   printout-device                                 6
!   print-level                                     0
!   maximum-number-of-iterations                    -1           
!   bi-diagonal-solve-frequency                     1
!   stopping-rule                                   1
!   number-extra-n-vectors-used                     0
!   relative-accuracy-required                      1.0E-8 
!   absolute-accuracy-required                      0.0
!   fraction-optimality-required                    1.0
!   stop-as-soon-as-boundary-encountered            T
!   space-critical                                  F
!   deallocate-error-fatal                          F
!   output-line-prefix                              ""
!  END LSRT SPECIFICATIONS

!  Dummy arguments

      TYPE ( LSRT_control_type ), INTENT( INOUT ) :: control        
      INTEGER, INTENT( IN ) :: device
      CHARACTER( LEN = 16 ), OPTIONAL :: alt_specname

!  Programming: Nick Gould and Ph. Toint, January 2002.

!  Local variables

      INTEGER, PARAMETER :: lspec = 34
      CHARACTER( LEN = 16 ), PARAMETER :: specname = 'LSRT            '
      TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec

!  Define the keywords

      spec%keyword = ''

!  Integer key-words

      spec(  1 )%keyword = 'error-printout-device'
      spec(  2 )%keyword = 'printout-device'
      spec(  3 )%keyword = 'print-level' 
      spec(  4 )%keyword = 'maximum-number-of-iterations'
      spec(  5 )%keyword = 'stopping-rule'
      spec(  9 )%keyword = 'number-extra-n-vectors-used'
      spec( 10 )%keyword = 'bi-diagonal-solve-frequency'

!  Real key-words

      spec(  6 )%keyword = 'relative-accuracy-required'
      spec(  7 )%keyword = 'absolute-accuracy-required'
      spec(  8 )%keyword = 'fraction-optimality-required'

!  Logical key-words

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

      CALL SPECFILE_assign_value( spec( 1 ), control%error,                     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 2 ), control%out,                       &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 3 ), control%print_level,               &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 4 ), control%itmax,                     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 5 ), control%stopping_rule,             &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 9 ), control%extra_vectors,             &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 10 ), control%freq,                     &
                                  control%error )

!  Set real values

      CALL SPECFILE_assign_value( spec( 6 ), control%stop_relative,             &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 7 ), control%stop_absolute,             &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 8 ), control%fraction_opt,              &
                                  control%error )
     
!  Set logical values

      CALL SPECFILE_assign_value( spec( 14 ), control%space_critical,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 15 ),                                   &
                                  control%deallocate_error_fatal,               &
                                  control%error )
!  Set charcter values

!     CALL SPECFILE_assign_value( spec( 16 ), control%prefix,                   &
!                                 control%error )

      RETURN

      END SUBROUTINE LSRT_read_specfile

!-*-*-*-*-*-*-*-*-*-*  L S R T _ S O L V E   S U B R O U T I N E  -*-*-*-*-*-*-*

      SUBROUTINE LSRT_solve( m, n, p, sigma, X, U, V, data, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!  =========
!
!   m        number of equations (= number of rows of A) > 0
!   n        number of unknowns (= number of columns of A) > 0
!   p        order of regularisation >= 2
!   sigma    regularisation weight
!   X        the vector of unknowns. Need not be set on entry.
!            On exit, the best value found so far
!   U        see inform%status = 2-4. On initial entry this must contain b
!   V        see inform%status = 2-4
!   data     private internal data
!   control  a structure containing control information. See LSRT_initialize
!   inform   a structure containing information. Components are
!             status     the input/output status. This must be set to 1 on
!                        initial entry or 5 on a re-entry when only sigma has
!                        been reduced since the last entry. Other values are
!               2 on exit, the product A with V must be added to U with
!                   the result kept in U, and the subroutine re-entered.
!               3 on exit, the product A^T with U must be added to V with
!                   the result kept in V, and the subroutine re-entered.
!               4 The iteration will be restarted. Reset R to c and re-enter.
!               0 the solution has been found
!              -1 an array allocation has failed
!              -2 an array deallocation has failed
!              -3 m <= 0 and/or n <= 0 and/or sigma < 0 and/or p < 2
!              -18 the iteration limit has been exceeded
!              -25 status is not > 0 on entry
!             iter        the number of bi-diagonalisation iterations performed
!             iter_pass2  the number of additional iterations performed to
!                         recover the solution 
!             multiplier  the multiplier sigma ||x||^(p-2)
!             x_norm      the norm of x
!             r_norm      the norm of Ax-b
!             ATr_norm    the norm of A^T (Ax-b) + multiplier * x
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  A description of LSQR and the notation used here is based on the papers
!
!    C. C. Paige and M. A. Saunders,
!    "LSQR: an algorithm for sparse linear equations and least squares problems",
!    ACM Transactions on Mathematical Software, Vol 8, No 1 (1982), pp 43-71,
!
!  and
!
!    C. C. Paige and M. A. Saunders,
!    "Algorithm 583. LSQR: sparse linear equations and least squares problems",
!    ACM Transactions on Mathematical Software, Vol 8, No 2 (1982), pp 195-209.
!
!  The bi-diagonalisation is due to Golub and Kahan (SIMUM 2, 1965, pp 205-224) 
!  (normalisation ||u_k|| = 1 = ||v_k||, 2-norm throughout).
!
!    Bi-diagonalisation initialization -
!
!       beta_1 u_1 = b   and
!      alpha_1 v_1 = A^T u_1
!
!    Bi-diagonalisation iteration (k = 1,2, ... ) -
!
!       beta_k+1 u_k+1 = A v_k - alpha_k u_k  and
!      alpha_k+1 v_k+1 = A^T u_k+1 - beta_k+1 v_k
!
!  This leads to
!
!    U_k+1^T A V_k = B_k  and  U_k+1^T b = beta_1 e_1 = ||b|| e_1 ,
!
!  where U_k = ( u_1 u_2 ... u_k ), V_k = ( v_1 v_2 ... v_k ) and
!
!         ( alpha_1                              )
!         (  beta_2  alpha_2                     )
!   B_k = (             .        .               )
!         (                    beta_k   alpha_k  )
!         (                             beta_k+1 )
!
!  To solve min 1/2|| A x - b ||^2 + 1/p sigma ||x||^p, find an 
!  approximation in the expanding subspace x = V_k y 
!  => find y_k = arg min || B_k y - beta_1 e_1 || 1/p sigma ||y||^p
!  and then recover x_k = V_k y_k in a second pass
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!---------------------------------
!   D u m m y   A r g u m e n t s
!---------------------------------

      INTEGER, INTENT( IN ) :: m, n
      REAL ( KIND = wp ), INTENT( IN ) :: p, sigma
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: X, V
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: U
      TYPE ( LSRT_data_type ), INTENT( INOUT ) :: data
      TYPE ( LSRT_control_type ), INTENT( IN ) :: control        
      TYPE ( LSRT_inform_type ), INTENT( INOUT ) :: inform

!---------------------------------
!   L o c a l   V a r i a b l e s
!---------------------------------

      INTEGER :: i
      REAL ( KIND = wp ) :: alpha_1, dec_tol, d, g, rho_tilde, phi_tilde
      REAL ( KIND = wp ) :: num, rat, x_kp1_norm, zeta, twoobj, rbiters
      CHARACTER ( LEN = 80 ) :: array_name

!  prefix for all output 

      CHARACTER ( LEN = LEN( TRIM( control%prefix ) ) - 2 ) :: prefix
      prefix = control%prefix( 2 : LEN( TRIM( control%prefix ) ) - 1 )

!  Branch to different sections of the code depending on input status

      IF ( inform%status < 1 ) GO TO 930
      IF ( inform%status == 1 ) data%branch = 1

      SELECT CASE ( data%branch )
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

      IF ( n <= 0 .OR. m < 0 ) GO TO 940
      IF ( sigma < zero ) GO TO 950
      IF ( p < two ) GO TO 980
      data%one_pass = sigma == zero .OR. p == two

!  Set iteration limits

      IF ( control%itmax >= 0 ) THEN
        data%itmax = control%itmax
      ELSE
        data%itmax = MAX( m, n ) + 1
      END IF

!  =====================
!  Array (re)allocations
!  =====================

      IF ( .NOT. data%one_pass ) THEN

!  Allocate space for the solution, y, in the subspace x = V_k y

        array_name = 'gltr: Y'
        CALL SPACE_resize_array( data%itmax, data%Y,                            &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  Allocate space for the Lanczos bidiagonal, B_k

        array_name = 'gltr: B_diag'
        CALL SPACE_resize_array( data%itmax + 1, data%B_diag,                   &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        array_name = 'gltr: B_offdiag'
        CALL SPACE_resize_array( data%itmax, data%B_offdiag,                    &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  Allocate space for the Cholesky factor R_k of B_k^T B_k

        array_name = 'gltr: R_diag'
        CALL SPACE_resize_array( data%itmax, data%R_diag,                       &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        array_name = 'gltr: R_offdiag'
        CALL SPACE_resize_array( data%itmax - 1, data%R_offdiag,                &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  Allocate worspace for the bi-diagonal solves

        array_name = 'gltr: F'
        CALL SPACE_resize_array( data%itmax, data%F,                            &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        array_name = 'gltr: G'
        CALL SPACE_resize_array( 0, data%itmax, data%G,                         &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        array_name = 'gltr: H'
        CALL SPACE_resize_array( data%itmax, data%H,                            &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  Allocate workspace for the sequence of smallest function values & multiplis

        array_name = 'gltr: DECREASE'
        CALL SPACE_resize_array( data%itmax, data%DECREASE,                     &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        array_name = 'gltr: LAMBDA'
        CALL SPACE_resize_array( data%itmax, data%LAMBDA,                       &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

!  If required allocate extra space to store U and V in case of a second pass

        data%extra_vectors = control%extra_vectors
        IF ( data%extra_vectors > 0 ) THEN
          array_name = 'gltr: U_extra'
          CALL SPACE_resize_array( m, data%U_extra,                             &
              inform%status, inform%alloc_status, array_name = array_name,      &
              deallocate_error_fatal = control%deallocate_error_fatal,          &
              exact_size = control%space_critical,                              &
              bad_alloc = inform%bad_alloc, out = control%error )
          IF ( inform%status == 0 ) THEN
            array_name = 'gltr: V_extra'
            CALL SPACE_resize_array( n, data%extra_vectors, data%V_extra,       &
                inform%status, inform%alloc_status, array_name = array_name,    &
                deallocate_error_fatal = control%deallocate_error_fatal,        &
                exact_size = control%space_critical,                            &
                bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= 0 ) THEN
              array_name = 'lsrt: U_extra'
              CALL SPACE_dealloc_array( data%U_extra,                           &
                inform%status, inform%alloc_status, array_name = array_name,    &
                bad_alloc = inform%bad_alloc, out = control%error )
              IF ( inform%status /= 0 ) GO TO 960
              data%extra_vectors = 0
            END IF
          ELSE
            data%extra_vectors = 0
          END IF
        END IF

!  Special case for the case p = 2

      ELSE

!  Allocate space to store the search direction W

        array_name = 'gltr: W'
        CALL SPACE_resize_array( n, data%W,                                     &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 960

        data%extra_vectors = 0
        data%omega = SQRT( sigma )
        data%g_norm2 = zero
        X = zero
      END IF

!  initialization

      inform%iter = 0 ; inform%iter_pass2 = 0
      inform%x_norm = zero ; data%z_norm = zero
      inform%multiplier = zero
      inform%biters = 0
      inform%biter_min = 0 ; inform%biter_max = 0 ; inform%biter_mean = zero

      data%printi = control%out > 0 .AND. control%print_level >= 1
      data%printd = control%out > 0 .AND. control%print_level >= 2
      data%error_tol = error_tol
      data%try_warm = .FALSE.
      data%freq = MAX( 1, control%freq ) ; data%bitmax = 10
      data%save_vectors = data%extra_vectors > 0

!  Recur the bi-diagonalisation vectors in u (= u_k) and v (= v_k)
!
!  Bi-diagonalisation initialization: beta_1 u_1 = b

      V = zero

!  Normalize u_1

      data%beta_1 = TWO_NORM( U )
      IF ( data%beta_1 > zero ) THEN
        U = U / data%beta_1

!  ** Return ** to obtain the product v <- v + A^T * u

        data%branch = 2 ; inform%status = 3
        RETURN
      END IF

!  Continue bi-diagonalisation initialization: alpha_1 v_1 = A^T u_1

  200 CONTINUE 

!  Normalize v_1 and initialize w

      alpha_1 = TWO_NORM( V )
      IF ( .NOT. data%one_pass ) data%B_diag( 1 ) = alpha_1
      IF ( alpha_1 > zero ) THEN
        V = V / alpha_1
        IF ( data%one_pass ) data%W( : n ) = V
      END IF
      data%ww = one

!  Initialize ||r_1|| and ||A^T r_1|| , where the residual r_1 = b - A x_1

      inform%r_norm  = data%beta_1
      inform%ATr_norm = alpha_1 * data%beta_1
      data%obj_0 = half * inform%r_norm ** 2 ; inform%obj = data%obj_0

!  Check for convergence

      IF ( ABS( inform%ATr_norm ) <= control%stop_absolute ) GO TO 900

!   Compute the stopping tolerance

      data%stop = MAX( control%stop_relative * inform%ATr_norm,                 &
                       control%stop_absolute )

!  Initialize R_1 and f_1

      data%rho_bar = alpha_1
      data%phi_bar = data%beta_1

!  Print initial information

      IF ( data%printi ) THEN
        IF ( data%one_pass ) THEN
          WRITE( control%out, 2010 ) prefix
          WRITE( control%out, "( A, I7, 2ES16.8, 2ES9.2 )" ) prefix,            &
             inform%iter, inform%obj, inform%r_norm, inform%ATr_norm,           &
             inform%x_norm
        ELSE 
          WRITE( control%out, 2000 ) prefix
          WRITE( control%out, "( A, I7, 2ES16.8, 2ES9.2 )" )                    &
             prefix, inform%iter, inform%obj, inform%r_norm, inform%ATr_norm,   &
             inform%x_norm
        END IF
      END IF

!  ===========================
!  Start of the main iteration
!  ===========================

      data%alpha_kp1 = alpha_1

  290   CONTINUE

!  Check that the iteration bound has not been exceeded

        IF ( inform%iter + 1 > data%itmax ) GO TO 490
        inform%iter = inform%iter + 1
        data%header = data%printd .OR. MOD( inform%iter, 25 ) == 0

!  If the user has asked to save vectors, save v and possibly u

        IF ( data%save_vectors ) THEN
          IF ( inform%iter <= data%extra_vectors )                              &
            data%V_extra( : n , inform%iter ) = V
          IF ( inform%iter == data%extra_vectors + 1 )                          &
            data%U_extra( : m ) = U
        END IF

!  Bi-diagonalisation iteration: 
!    beta_k+1 u_k+1 = A v_k - alpha_k u_k

!  ** Return ** to obtain the product u <- u + A * v

        U = - data%alpha_kp1 * U
        data%branch = 3 ; inform%status = 2
        RETURN

!  Normalize u

  300   CONTINUE 
        data%beta_kp1 = TWO_NORM( U )
        IF ( .NOT. data%one_pass ) data%B_offdiag( inform%iter ) = data%beta_kp1

        IF ( data%beta_kp1 > zero ) THEN
          U = U / data%beta_kp1

!  Continue bi-diagonalisation iteration: 
!    alpha_k+1 v_k+1 = A^T u_k+1 - beta_k+1 v_k

!  ** Return ** to obtain the product v <- v + A^T * u

          V = - data%beta_kp1 * V
          data%branch = 4 ; inform%status = 3
          RETURN
        END IF

!  Normalize v

  400   CONTINUE
        data%alpha_kp1 = TWO_NORM( V )
        IF ( data%alpha_kp1 > zero ) THEN
          V = V / data%alpha_kp1
        END IF
        
!  general case

        IF ( .NOT. data%one_pass ) THEN
          data%B_diag( inform%iter + 1 ) = data%alpha_kp1

!  Find the solution on the boundary in the subspace x = V_k y

          IF ( data%try_warm )                                                  &
            data%LAMBDA( inform%iter ) = data%LAMBDA( inform%iter - 1 )
          CALL LSRT_solve_bidiagonal( inform%iter,                              &
            data%B_diag( : inform%iter ), data%B_offdiag( : inform%iter ),      &
            data%beta_1, p, sigma, data%LAMBDA( inform%iter ),                  &
            data%Y( : inform%iter ), data%H( : inform%iter ), inform%x_norm,    &
            data%R_diag( : inform%iter ), data%R_offdiag( : inform%iter - 1 ),  &
            data%F( : inform%iter ), data%G( : inform%iter ), data%error_tol,   &
            data%try_warm, control%print_level, control%out, prefix,            &
            data%bitmax, data%biter, data%bstatus )

!  Record statistics about the number of inner iterations performed 

          IF ( INFORM%biters > 0 ) THEN
            rbiters = INFORM%biters
            inform%biters = inform%biters + 1
            inform%biter_min = MIN( inform%biter_min, data%biter )
            inform%biter_max = MAX( inform%biter_max, data%biter )
            inform%biter_mean = ( rbiters / ( rbiters + one ) ) *               &
              ( inform%biter_mean + data%biter / rbiters )
          ELSE
            inform%biters = 1
            inform%biter_min = data%biter
            inform%biter_max = data%biter
            inform%biter_mean = data%biter
          END IF

!  Compute the norms ||r_k|| and ||A^T r_k||

          data%try_warm = .TRUE.
          inform%r_norm = SQRT( ABS( DOT_PRODUCT( data%G( : inform%iter ),      &
                                                  data%G( : inform%iter ) ) -   &
                          data%LAMBDA( inform%iter ) * ( inform%x_norm ** 2 ) ) )
          inform%ATr_norm = ABS( data%alpha_kp1 * data%B_offdiag( inform%iter ) &
                                 * data%Y( inform%iter ) )
          inform%obj = half * inform%r_norm ** 2 +                              &
                         ( sigma / p ) * inform%x_norm ** p
          data%DECREASE( inform%iter ) = data%obj_0 - inform%obj

!  Print progress

          IF ( data%printi ) THEN
            IF ( data%header ) WRITE( control%out, 2000 ) prefix
            WRITE( control%out, "( A, I7, 2ES16.8, 3ES9.2, I7 )" ) prefix,      &
               inform%iter, inform%obj, inform%r_norm, inform%ATr_norm, sigma,  &
               data%LAMBDA( inform%iter ), data%biter
          END IF

!  Check for termination

          IF ( inform%ATr_norm < data%stop ) GO TO 490

!  p = 2 case

        ELSE

!  Columns k, (k+1) and (n+1) contain entries

!       (  rho_bar     0      phi_bar )
!       ( beta_k+1 alpha_k+1     0    )
!       (  omega       0         0    )

!  in rows k, (k+1) and (n+1). Proceed as in Paige and Saunders, 
!  ACM TOMS, Vol 8, No 2 (1982), pp 195-209, Section 2

!  Construct a plane rotation P to eliminate omega from ( rho_bar omega )
!  to create ( rho_tilde 0 )

          IF ( sigma > zero ) THEN
            d = data%omega
            CALL ROTG( data%rho_bar, d, data%c, data%s )
            rho_tilde = data%rho_bar

!  Apply the rotation P to ( phi_bar 0 ) to create ( phi_tilde  g )

            phi_tilde = data%c * data%phi_bar
            g = data%s * data%phi_bar
            data%g_norm2 = data%g_norm2 + g ** 2
          ELSE
            rho_tilde = data%rho_bar
            phi_tilde = data%phi_bar
          END IF

!  Construct a plane rotation Q to eliminate b from ( r_diag_tilde b )
!  to create ( r_diag 0 )

          CALL ROTG( rho_tilde, data%beta_kp1, data%c, data%s )
          data%rho_k = rho_tilde

!  Apply the rotation Q to ( phi_tilde 0 ) to create ( f_k phi_bar )

          data%phi_k = data%c * phi_tilde
          data%phi_bar = data%s * phi_tilde

!  Construct a plane rotation W to remove the super-diagonal entry theta_k from
!  R_k (and thus incrementally turn R_k into a lower-bidiagonal matrix L_k)

          IF ( inform%iter > 1 ) THEN
            data%lambda_km1 = data%lambda_bar
            CALL ROTG( data%lambda_km1, data%theta_kp1, data%c_w, data%s_w )
            zeta = data%zeta_bar * ( data%lambda_bar ) / data%lambda_km1
            data%z_norm = ROOT_SUM_SQUARES( data%z_norm, zeta )

!  Apply the rotation W to the vector ( 0 rho_k )

            data%gamma = data%s_w * data%rho_k
            data%lambda_bar = data%c_w * data%rho_k

!  Compute the norm of ||x_k|| - see Paige and Saunders ACM TOMS 8(1) 1982, 
!  Sec 5.2 for a brief summary

            num = data%phi_k - data%gamma * zeta
            data%zeta_bar = num / data%lambda_bar
            x_kp1_norm = ROOT_SUM_SQUARES( data%z_norm, data%zeta_bar ) 
          ELSE
            data%lambda_bar = data%rho_k
            data%z_norm = zero
            data%zeta_bar = data%phi_k / data%lambda_bar
            x_kp1_norm = ABS( data%zeta_bar )
          END IF
          inform%x_norm = x_kp1_norm

!  Apply the rotation Q to the forthcoming column ( 0 alpha_k+1 ) 
!  in B to create ( r_offdiag r_diag_bar )

          data%theta_kp1 = data%s * data%alpha_kp1
          data%rho_bar = - data%c * data%alpha_kp1

!  Update the approximate solution x and the correction vector w

          X = X + ( data%phi_k / data%rho_k ) * data%W( : n )
          rat =  data%theta_kp1 / data%rho_k
          data%W( : n ) = V - rat * data%W( : n )

!  Update the square of the norm of w

          data%ww = one + data%ww * rat ** 2

!  Compute the norms ||r_k|| and ||A^T r_k|| - see Paige and Saunders 
!  ACM TOMS 8(1) 1982, Sec 5.1

          twoobj = data%phi_bar ** 2 + data%g_norm2
          inform%obj = half * twoobj
          inform%r_norm = SQRT( twoobj - sigma * inform%x_norm ** 2 )
          inform%ATr_norm = ABS( data%alpha_kp1 * data%c * data%phi_bar )

!  Print progress

          IF ( data%printi ) THEN
            IF ( data%header ) WRITE( control%out, 2010 ) prefix
            WRITE( control%out, "( A, I7, 2ES16.8, 2ES9.2 )" ) prefix,          &
              inform%iter, inform%obj, inform%r_norm, inform%ATr_norm,          &
              inform%x_norm
          END IF

!  Check for termination

          IF ( inform%ATr_norm < data%stop ) GO TO 900
        END IF

!  =========================
!  End of the main iteration
!  =========================

      GO TO 290

!  Termination has occured. Determine at which iteration a fraction, 
!  fraction_opt, of the optimal solution was found

  490 CONTINUE

      IF ( inform%iter == 0 ) THEN
        X = zero
        GO TO 920
      END IF  

      IF ( control%fraction_opt < one ) THEN
        dec_tol = data%DECREASE( inform%iter ) * control%fraction_opt

!  Examine subsequent iterates

        DO i = 1, inform%iter
          data%end_pass2 = i
          IF ( data%DECREASE( data%end_pass2 ) >= dec_tol ) EXIT
        END DO
      ELSE
        data%end_pass2 = inform%iter
      END IF

!  Recover the solution on the boundary in the subspace x = V_pass2 y

      CALL LSRT_solve_bidiagonal( data%end_pass2,                               &
        data%B_diag( : data%end_pass2 ),                                        &
        data%B_offdiag( : data%end_pass2 ),                                     &
        data%beta_1, p, sigma, data%LAMBDA( data%end_pass2 ),                   &
        data%Y( : data%end_pass2 ), data%H( : data%end_pass2 ),                 &
        inform%x_norm, data%R_diag( : data%end_pass2 ),                         &
        data%R_offdiag( : data%end_pass2 - 1 ),                                 &
        data%F( : data%end_pass2 ), data%G( : data%end_pass2 ),                 &
        data%error_tol, .TRUE., control%print_level, control%out, prefix,       &
        data%bitmax, data%biter, data%bstatus )

!  Compute the norms ||r_pass2|| and ||A^T r_pass2|| and the Lagrange multiplier

      inform%multiplier = data%Y( data%end_pass2 )
      inform%r_norm = SQRT( ABS( DOT_PRODUCT( data%G( : data%end_pass2 ),       &
                                              data%G( : data%end_pass2 ) ) -    &
                      data%LAMBDA( data%end_pass2 ) * ( inform%x_norm ** 2 ) ) )
      inform%ATr_norm = ABS( data%B_diag( data%end_pass2 + 1 ) *                &
              data%B_offdiag( data%end_pass2 ) * data%Y( data%end_pass2 ) )
      inform%obj = data%obj_0 - data%DECREASE( data%end_pass2 )

!  Restart to compute x_pass2

      IF ( data%extra_vectors == 0 ) THEN
        data%branch = 5 ; inform%status = 4
        RETURN
      END IF

!  =================================================
!  Second pass to recover the solution x_k = V_k y_k 
!  from the vectors V_k = ( v_1 : ... : v_k )
!  =================================================

  500 CONTINUE

      IF ( data%extra_vectors == 0 ) THEN
        X = zero
        inform%iter_pass2 = 0

!  Reentry with no saved vectors: Normalize u_1 and initialize v_1

        U = U / data%beta_1
        V = zero
      ELSE

!  Reentry with  data%extra_vectors saved vectors: loop over saved v_i to 
!  update x_k = sum_i v_i y_i

        inform%iter_pass2 = MIN( data%end_pass2, data%extra_vectors )
        X = MATMUL( data%V_extra( : n, : inform%iter_pass2 ),                   &
                    data%Y( : inform%iter_pass2 ) )

!  Check to see if the solution has been recovered

        IF ( inform%iter_pass2 >= data%end_pass2 ) GO TO 900

!  Initialize u_k and - beta_k v_k-1

        U = data%U_extra( : m )
        V = - data%B_offdiag( inform%iter_pass2 ) *                             &
                data%V_extra( : n, inform%iter_pass2 )
      END IF

!  ===============================
!  Iterate over the columns of V_k
!  ===============================

  590   CONTINUE

!  Bi-diagonalisation iteration: 
!    alpha_k v_k = A^T u_k - beta_k v_k-1

!  ** Return ** to obtain the product v <- v + A^T * u

        data%branch = 6 ; inform%status = 3
        RETURN

!  Normalize v_k

  600   CONTINUE 
        inform%iter_pass2 = inform%iter_pass2 + 1
        V = V / data%B_diag( inform%iter_pass2 )

!  Update x_k

        X = X + data%Y( inform%iter_pass2 ) * V

!  Check to see if the solution has been recovered

        IF ( inform%iter_pass2 >= data%end_pass2 ) GO TO 900

!  Bi-diagonalisation iteration: 
!    beta_k+1 u_k+1 = A v_k - alpha_k u_k

!  ** Return ** to obtain the product u <- u + A * v

        U = - data%B_diag( inform%iter_pass2 ) * U
        data%branch = 7 ; inform%status = 2
        RETURN

!  Normalize u_k+1 and scale v_k

  700   CONTINUE 
        U = U / data%B_offdiag( inform%iter_pass2 )
        V = - data%B_offdiag( inform%iter_pass2 ) * V

!  ======================
!  End of the second pass
!  ======================

      GO TO 590

!  Successful returns

  900 CONTINUE
      inform%status = 0
      RETURN

!  Too many iterations

  920 CONTINUE
      IF ( data%printi )                                                        &
        WRITE( control%out, "( /, A, ' Iteration limit exceeded ' ) " ) prefix
      inform%status = GALAHAD_error_max_iterations
      RETURN

!  Inappropriate entry status

  930 CONTINUE
      inform%status = GALAHAD_error_input_status
      RETURN

!  Inappropriate dimension

  940 CONTINUE
      IF ( control%error > 0 .AND. control%print_level > 0 )                    &
        WRITE( control%error,                                                   &
         "( A, ' n = ', I0, ' or m = ', I0, ' is too small ' )" ) prefix, n, m
      inform%status = GALAHAD_error_restrictions
      RETURN

!  Negative regularisation weight

  950 CONTINUE
      IF ( control%error > 0 .AND. control%print_level > 0 )                    &
        WRITE( control%error,                                                   &
         "( A, ' sigma ', ES12.4 , ' is negative ' )" ) prefix, sigma
      inform%status = GALAHAD_error_restrictions
      RETURN

!  Allocation or deallocation errors

  960 CONTINUE
      RETURN

!  Regularisation order to small

  980 CONTINUE
      IF ( control%error > 0 .AND. control%print_level > 0 )                    &
        WRITE( control%error,                                                   &
         "( A, ' p ', ES12.4 , ' is smaller than 2.0 ' )" ) prefix, p
      inform%status = GALAHAD_error_restrictions
      RETURN

!  Non-executable statement

 2000 FORMAT( /, A, '   Iter        obj          ||Ax-b||  ',                 &
                    '  ||A^Tr||   ||x||   lambda  Newton' )
 2010 FORMAT( /, A, '   Iter        obj          ||Ax-b||    ||A^Tr||   ||x|| ' )

      CONTAINS

! -*-*- R O O T _ S U M _ S Q U A R E S    I N T E R N A L   F U N C T I O N -*-

        FUNCTION ROOT_SUM_SQUARES( a, b )

!  Evaluates the square root of the sum of squares of a and b

!  Dummy arguments

        REAL ( KIND = wp ) :: ROOT_SUM_SQUARES
        REAL ( KIND = wp ), INTENT( IN ) :: a, b

!  Local variable

        REAL ( KIND = wp ) :: s

!  Take precautions to try to prevent overflow and underflow

        s = ABS( a ) + ABS( b )
        IF ( s == zero ) THEN
           ROOT_SUM_SQUARES = zero
        ELSE
           ROOT_SUM_SQUARES = s * SQRT( ( a / s ) ** 2 + ( b / s ) ** 2 )
        END IF

!  End of function ROOT_SUM_SQUARES

        RETURN
        END FUNCTION ROOT_SUM_SQUARES

!  End of subroutine LSRT_solve

      END SUBROUTINE LSRT_solve

!-*-*-*-*-*-  L S R T _ T E R M I N A T E   S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE LSRT_terminate( data, control, inform )

!  ..............................................
!  .                                            .
!  .  Deallocate arrays at end of LSRT_solve    .
!  .                                            .
!  ..............................................

!  Arguments:
!  =========
!
!   data    private internal data
!   control see Subroutine LSRT_initialize
!   inform    see Subroutine LSRT_solve

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      TYPE ( LSRT_data_type ), INTENT( INOUT ) :: data
      TYPE ( LSRT_control_type ), INTENT( IN ) :: control        
      TYPE ( LSRT_inform_type ), INTENT( INOUT ) :: inform

!-----------------------------------------------
!   L o c a l   V a r i a b l e
!-----------------------------------------------

      CHARACTER ( LEN = 80 ) :: array_name

!  Deallocate all internal arrays

      array_name = 'lsrt: Y'
      CALL SPACE_dealloc_array( data%Y,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: V'
      CALL SPACE_dealloc_array( data%H,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: F'
      CALL SPACE_dealloc_array( data%F,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: LAMBDA'
      CALL SPACE_dealloc_array( data%LAMBDA,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: B_diag'
      CALL SPACE_dealloc_array( data%B_diag,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: B_offdiag'
      CALL SPACE_dealloc_array( data%B_offdiag,                                &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: DECREASE'
      CALL SPACE_dealloc_array( data%DECREASE,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: R_diag'
      CALL SPACE_dealloc_array( data%R_offdiag,                                &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: R_offdiag'
      CALL SPACE_dealloc_array( data%R_offdiag,                                &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: U_extra'
      CALL SPACE_dealloc_array( data%U_extra,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'lsrt: V_extra'
      CALL SPACE_dealloc_array( data%V_extra,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      RETURN

!  End of subroutine LSRT_terminate

      END SUBROUTINE LSRT_terminate

!-*-*-*-  L S R T _ S O L V E _ B I D I A G O N A L   S U B R O U T I N E  *-*-*-

      SUBROUTINE LSRT_solve_bidiagonal( n, B_diag, B_offdiag, beta, p, sigma,   &
                                        lambda, Y, H, y_norm, R_diag, R_offdiag,&
                                        F, G, error_tol, try_warm, print_level, &
                                        out, prefix, itmax, iter,  status )

! ---------------------------------------------------------------------------
!
!  Solve min 1/2 || B y - beta e_1 ||^2 + 1/p sigma || y ||^p
!  where
!           (  diag_1                              )
!           ( offdiag_1  diag_2                    )
!     B   = (              .       .               )  is (n+1) by n
!           (                   offdiag_n  diag_n  )
!           (                            offdiag_n )
!
!  The optimality conditions are that y(lambda) solves the weighted problem
!
!     min || B y( lambda ) - beta e_1 ||^2 + lambda || y(lambda) ||^2
!
!  or equivalently
!
!     min || (        B       ) y( lambda ) - beta e_1 ||^2            (*)
!         || ( srqt(lambda) I )                        ||
!
!  and where 
!
!              || y(lambda) ||^(p-2) = lambda / sigma                  (**) 
!
!  So pick lambda to define y(lambda) via (*) and adjust it to enforce 
!  the scalar equation (**) or variants using Newton's method. This 
!  requires that we can solve (*), which we do by reducing
!
!      (        B         beta e_1 ) -> ( R  f )
!      ( srqt(lambda) I      0     )    ( 0  g )
!
!  to n by n upper bi-diagonal form by pre-multiplying by plane rotations.
!  In this case R y(lambda) = f
!
!  Alternative secular functions (**) considered
!    theta(lambda) = || y(lambda) ||^(p-2) - lambda/sigma
!    phi(lambda) = || y(lambda) ||^(2-p) - sigma/lambda
!    zeta(lambda) = lambda * || y(lambda) ||^(2-p) - sigma
!
! ---------------------------------------------------------------------------

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n, itmax, out, print_level
      INTEGER, INTENT( OUT ) :: status, iter 
      LOGICAL, INTENT( IN ) :: try_warm
      REAL ( KIND = wp ), INTENT( IN ) :: p, sigma, beta, error_tol
      REAL ( KIND = wp ), INTENT( INOUT ) :: lambda, y_norm
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: B_diag, B_offdiag
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n - 1 ) :: R_offdiag
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: R_diag, Y, H, F
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( 0 : n ) :: G
      CHARACTER ( LEN = * ), INTENT( IN ) :: prefix

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: it1, nroots
      LOGICAL :: printm, printd
      REAL ( KIND = wp ) :: h_norm, v2oy2, root1, root2, delta, error
      REAL ( KIND = wp ) :: dl_theta, dl_phi, dl_zeta, dl_phi_c
      REAL ( KIND = wp ) :: omega, omega_prime, phi, phi_prime
      REAL ( KIND = wp ) :: zeta, zeta_prime, pi, pi_prime, theta, theta_prime

      printm = out > 0 .AND. print_level >= 2
      printd = out > 0 .AND. print_level >= 3

      IF ( printm ) THEN
         WRITE( out, "( /, A, '   Bi-diagonal subproblem - ' )" ) prefix
         WRITE( out, 2000 ) prefix
      END IF

!  Choose starting lambda ... a warm start

      IF ( try_warm ) THEN
        it1 = 2

!  Transform the bi-diagonal subproblem to upper-traingular form R(lambda) 
!  for the current lambda
        
        CALL LSRT_transform_bidiagonal( n, B_diag, B_offdiag, beta,             &
                                        SQRT( lambda ), R_diag, R_offdiag, F, G )

!  Compute the solution y(lambda) to the bi-diagonal subproblem

        CALL LSRT_backsolve_bidiagonal( n, R_diag, R_offdiag, F, Y, .FALSE. )
        y_norm = TWO_NORM( Y )
        delta = lambda / sigma
        error = y_norm ** ( p - two ) - delta
        IF ( printm ) WRITE( out, "( A, I7, ES22.14, 2ES12.4 )" )               &
          prefix, 1, error, lambda, y_norm

!  Test for convergence 

        IF ( ABS( error ) < error_tol ) THEN
          iter = 1 ; status = 0
          RETURN
        END IF

!  Is this a suitable starting point?

        IF ( error > zero ) THEN

!  Compute the solution R^T(lambda) v = y

          CALL LSRT_backsolve_bidiagonal( n, R_diag, R_offdiag, Y, H, .TRUE. )
          h_norm = TWO_NORM( H )

!  Compute theta = || y(lambda) || and its derivative

          omega = y_norm
          omega_prime = - ( h_norm / y_norm ) * h_norm

!  A suitable correction for lambda is the positive root of
!    sigma/lambda = pi + ( lambda - lambda_current ) pi' ;
!  here the right-hand side is the linearization of pi = 1/||x(lambda)||^(p-2)

          pi = y_norm ** ( two - p )
          pi_prime = ( two - p ) * ( omega ** ( one - p ) ) * omega_prime

          v2oy2 = pi_prime / pi
          CALL ROOTS_quadratic( - sigma / pi, one - lambda * v2oy2, v2oy2,  &
                                roots_tol, nroots, root1, root2 )
          IF ( nroots == 2 ) THEN
            dl_phi_c = root2 - lambda
          ELSE  
            dl_phi_c = root1 - lambda
          END IF

!  Alternatives are Newton corrections to the functions given in the ACO paper

          IF ( printd ) THEN

!  theta correction

            theta = y_norm ** ( p - two ) - lambda / sigma
            theta_prime = ( p - two ) * ( omega ** ( p - three ) ) *            &
                  omega_prime - one / sigma
            dl_theta = - theta / theta_prime 

!  phi correction

            IF ( lambda /= zero ) THEN
              phi = pi - sigma / lambda
              phi_prime = pi_prime + sigma / ( lambda ** 2 )
              dl_phi = - phi / phi_prime 
            ELSE
              dl_phi = zero
            END IF

!  zeta correction

            zeta = lambda * pi - sigma
            zeta_prime = pi + lambda * pi_prime 
            dl_zeta = - zeta / zeta_prime 

!write(6,*) theta, phi, zeta
          
            WRITE( out, 2010 )                                                  &
              prefix, dl_phi_c, prefix, dl_theta, dl_phi, dl_zeta
          END IF

!  Compute the Newton-like update

          lambda = lambda + dl_phi_c
!         lambda = lambda + dl_theta

!  The warm start failed - revert to a cautious starting guess

        ELSE
          lambda = lambda_start
        END IF

!  Choose starting lambda ... a cold start

      ELSE
        it1 = 1
        lambda = lambda_start
      END IF

!  Newton iteration

      DO iter = it1, itmax

!  Transform the bi-diagonal subproblem to upper-traingular form R(lambda) 
!  for the current lambda
        
        CALL LSRT_transform_bidiagonal( n, B_diag, B_offdiag, beta,             &
                                        SQRT( lambda ), R_diag, R_offdiag, F, G )

!  Compute the solution y(lambda) to the bi-diagonal subproblem

        CALL LSRT_backsolve_bidiagonal( n, R_diag, R_offdiag, F, Y, .FALSE. )
        y_norm = TWO_NORM( Y )
        delta = lambda / sigma
        error = y_norm ** ( p - two ) - delta
        IF ( printd ) WRITE( out, 2000 ) prefix
        IF ( printm ) WRITE( out, "( A, I7, ES22.14, 2ES12.4 )" )               &
          prefix, iter, error, lambda, y_norm

!  Test for convergence 

        IF ( ABS( error ) < error_tol ) THEN
          status = 0 ; RETURN
        END IF

!  Compute the solution R^T(lambda) v = y

        CALL LSRT_backsolve_bidiagonal( n, R_diag, R_offdiag, Y, H, .TRUE. )
        h_norm = TWO_NORM( H )

!  Compute theta = || y(lambda) || and its derivative

        omega = y_norm
        omega_prime = - ( h_norm / y_norm ) * h_norm

!  A suitable correction for lambda is the positive root of
!    sigma/lambda = pi + ( lambda - lambda_current ) pi' ;
!  here the right-hand side is the linearization of pi = 1/||x(lambda)||^(p-2)

        pi = y_norm ** ( two - p )
        pi_prime = ( two - p ) * ( omega ** ( one - p ) ) * omega_prime

        v2oy2 = pi_prime / pi
        CALL ROOTS_quadratic( - sigma / pi, one - lambda * v2oy2, v2oy2,       &
                              roots_tol, nroots, root1, root2 )
        IF ( nroots == 2 ) THEN
          dl_phi_c = root2 - lambda
        ELSE  
          dl_phi_c = root1 - lambda
        END IF

!  Alternatives are Newton corrections to the functions given in the ACO paper

        IF ( printd ) THEN

!  theta correction

          theta = y_norm ** ( p - two ) - lambda / sigma
          theta_prime = ( p - two ) * ( omega ** ( p - three ) ) *              &
                omega_prime - one / sigma
          dl_theta = - theta / theta_prime 

!  phi correction

          IF ( lambda /= zero ) THEN
            phi = pi - sigma / lambda
            phi_prime = pi_prime + sigma / ( lambda ** 2 )
            dl_phi = - phi / phi_prime 
          ELSE
            dl_phi = zero
          END IF

!  zeta correction

          zeta = lambda * pi - sigma
          zeta_prime = pi + lambda * pi_prime 
          dl_zeta = - zeta / zeta_prime 
        
!write(6,*) theta, phi, zeta
          WRITE( out, 2010 )                                                    &
            prefix, dl_phi_c, prefix, dl_theta, dl_phi, dl_zeta

        END IF

!  Compute the Newton-like update

        lambda = lambda + dl_phi_c
      END DO

      status = 1
      RETURN

!  Non-executable statement

 2000 FORMAT( /, A, '   Iter    ||y|| - sigma        lambda       ||y||' )
 2010 FORMAT( A, ' correction ', ES16.8, ' alternatives',/,          &
              A, ' theta      ', ES16.8, ' phi ', ES16.8, ' zeta ', ES16.8 )

!  End of subroutine LSRT_solve_bidiagonal

      END SUBROUTINE LSRT_solve_bidiagonal

!-*-*-*-*-*-  End of G A L A H A D _ L S R T  double  M O D U L E  *-*-*-*-*-*-

   END MODULE GALAHAD_LSRT_double

