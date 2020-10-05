! THIS VERSION: GALAHAD 2.2 - 22/04/2008 AT 13:00 GMT.

!-*-*-*-*-*-*-*-*-*-  G A L A H A D _ L S Q P    M O D U L E  -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   started life as GALAHAD_CLS ~ 2000
!   originally released pre GALAHAD Version 1.0. April 10th 2001
!   update released with GALAHAD Version 2.0. November 1st 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_LSQP_double

!     ------------------------------------------------
!     |                                              |
!     | Minimize the linear/separable objective      |
!     |                                              |
!     |  1/2 || W * ( x - x^0 ) ||_2^2 + g^T x + f   |
!     |                                              |
!     | subject to the linear constraints and bounds |
!     |                                              |
!     |          c_l <= A x <= c_u                   |
!     |          x_l <=  x <= x_u                    |
!     |                                              |
!     | for some (possibly zero) diagonal matrix W,  |
!     | using an infeasible-point primal-dual method |
!     |                                              |
!     ------------------------------------------------

!NOT95USE GALAHAD_CPU_time
      USE GALAHAD_SYMBOLS
      USE GALAHAD_SILS_double
      USE GALAHAD_SMT_double
      USE GALAHAD_QPT_double
      USE GALAHAD_SPECFILE_double
      USE GALAHAD_QPP_double, LSQP_dims_type => QPP_dims_type
      USE GALAHAD_QPD_double, LSQP_data_type => QPD_data_type,                 &
                              LSQP_AX => QPD_AX

      USE GALAHAD_SORT_double, ONLY: SORT_heapsort_build,                      &
         SORT_heapsort_smallest, SORT_inverse_permute
      USE GALAHAD_FDC_double
      USE GALAHAD_ROOTS_double

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: LSQP_initialize, LSQP_read_specfile, LSQP_solve,               &
                LSQP_terminate, QPT_problem_type, SMT_type, SMT_put, SMT_get,  &
                LSQP_form_Schur_complement, LSQP_A_by_cols, LSQP_Ax,           &
                LSQP_data_type, LSQP_dims_type, LSQP_indicators

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

      TYPE, PUBLIC :: LSQP_time_type
        REAL :: total, preprocess, find_dependent, analyse, factorize, solve
      END TYPE

      TYPE, PUBLIC :: LSQP_control_type
        INTEGER :: error, out, print_level, start_print, stop_print, maxit
        INTEGER :: factor, max_col, indmin, valmin, itref_max, infeas_max
        INTEGER :: muzero_fixed, indicator_type, restore_problem
        REAL ( KIND = wp ) :: infinity, stop_p, stop_d, stop_c, prfeas, dufeas
        REAL ( KIND = wp ) :: muzero, reduce_infeas, potential_unbounded
        REAL ( KIND = wp ) :: pivot_tol, pivot_tol_for_dependencies, zero_pivot
        REAL ( KIND = wp ) :: identical_bounds_tol, indicator_tol_tapia
        REAL ( KIND = wp ) :: indicator_tol_p, indicator_tol_pd, cpu_time_limit
        LOGICAL :: remove_dependencies, treat_zero_bounds_as_general
        LOGICAL :: just_feasible, getdua, feasol, balance_initial_complentarity
        LOGICAL :: use_corrector, array_syntax_worse_than_do_loop
        CHARACTER ( LEN = 30 ) :: prefix
      END TYPE

      TYPE, PUBLIC :: LSQP_inform_type
        INTEGER :: status, alloc_status, iter, factorization_status
        INTEGER :: factorization_integer, factorization_real, nfacts, nbacts
        REAL ( KIND = wp ) :: obj, potential, non_negligible_pivot
        LOGICAL :: feasible
        CHARACTER ( LEN = 80 ) :: bad_alloc
        TYPE ( LSQP_time_type ) :: time
      END TYPE

!----------------------
!   P a r a m e t e r s
!----------------------

      INTEGER, PARAMETER :: max_sc = 200
      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: point1 = 0.1_wp
      REAL ( KIND = wp ), PARAMETER :: point01 = 0.01_wp
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: two = 2.0_wp
      REAL ( KIND = wp ), PARAMETER :: three = 3.0_wp
      REAL ( KIND = wp ), PARAMETER :: four = 4.0_wp
      REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp
      REAL ( KIND = wp ), PARAMETER :: thousand = 1000.0_wp
      REAL ( KIND = wp ), PARAMETER :: tenm4 = ten ** ( - 4 )
      REAL ( KIND = wp ), PARAMETER :: tenm5 = ten ** ( - 5 )
      REAL ( KIND = wp ), PARAMETER :: tenm7 = ten ** ( - 7 )
      REAL ( KIND = wp ), PARAMETER :: tenm10 = ten ** ( - 10 )
      REAL ( KIND = wp ), PARAMETER :: ten5 = ten ** 5
      REAL ( KIND = wp ), PARAMETER :: infinity = HUGE( one )
      REAL ( KIND = wp ), PARAMETER :: epsmch = EPSILON( one )

   CONTAINS

!-*-*-*-*-*-   L S Q P _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE LSQP_initialize( data, control )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Default control data for LSQP. This routine should be called before
!  LSQP_primal_dual
! 
!  --------------------------------------------------------------------
!
!  Arguments:
!
!  data     private internal data
!  control  a structure containing control information. Components are -
!
!  INTEGER control parameters:
!
!   error. Error and warning diagnostics occur on stream error 
!   
!   out. General output occurs on stream out
!   
!   print_level. The level of output required is specified by print_level
!   
!   maxit. At most maxit inner iterations are allowed 
!   
!   start_print. Any printing will start on this iteration
!
!   stop_print. Any printing will stop on this iteration
!
!   factor. The factorization to be used.
!    Possible values are
!
!      0  automatic 
!      1  Schur-complement factorization
!      2  augmented-system factorization
!
!   max_col. The maximum number of nonzeros in a column of A which is permitted
!    with the Schur-complement factorization
!
!   indmin. An initial guess as to the integer workspace required by SILS
!
!   valmin. An initial guess as to the real workspace required by SILS
! 
!   itref_max. The maximum number of iterative refinements allowed
!
!   infeas_max. The number of iterations for which the overall infeasibility
!     of the problem is not reduced by at least a factor control%reduce_infeas
!     before the problem is flagged as infeasible (see reduce_infeas)
!
!   muzero_fixed. The initial value of the barrier parameter will not be
!    changed for the first muzero_fixed iterations
!
!   restore_problem. Indicates whether and how much of the input problem
!    should be restored on output. Possible values are
!
!      0 nothing restored
!      1 scalar and vector parameters
!      2 all parameters
!
!   indicator_type. Specifies the type of indicator function used.
!    Pssible values are
!
!     1 primal indicator: constraint active <=> distance to nearest bound 
!         <= %indicator_p_tol
!     2 primal-dual indicator: constraint active <=> distance to nearest bound 
!        <= %indicator_tol_pd * size of corresponding multiplier
!     3 primal-dual indicator: constraint active <=> distance to nearest bound 
!        <= %indicator_tol_tapia * distance to same bound at previous iteration
!
!  REAL control parameters:
!
!   infinity. Any bound larger than infinity in modulus will be regarded as 
!    infinite 
!   
!   stop_p. The required accuracy for the primal infeasibility
!   
!   stop_d. The required accuracy for the dual infeasibility
!   
!   stop_c. The required accuracy for the complementarity
!   
!   prfeas. The initial primal variables will not be closer than prfeas 
!    from their bounds 
!   
!   dufeas. The initial dual variables will not be closer than dufeas from 
!    their bounds 
!   
!   muzero. The initial value of the barrier parameter. If muzero is
!    not positive, it will be reset to an appropriate vale
!           
!   reduce_infeas. If the overall infeasibility of the problem is not reduced 
!    by at least a factor reduce_infeas over control%infeas_max iterations,
!    the problem is flagged as infeasible (see infeas_max)
!
!   potential_unbounded. If W=0 and the potential function value is smaller 
!    than potential_unbounded * number of one-sided bounds, the analytic center 
!    will be flagged as unbounded
!
!   pivot_tol. The threshold pivot used by the matrix factorization.
!    See the documentation for SILS for details
!
!   pivot_tol_for_dependencies. The threshold pivot used by the matrix 
!    factorization when attempting to detect linearly dependent constraints.
!    See the documentation for SILS for details
!
!   zero_pivot. Any pivots smaller than zero_pivot in absolute value will 
!    be regarded to be zero when attempting to detect linearly dependent 
!    constraints
!
!   identical_bounds_tol. Any pair of constraint bounds (c_l,c_u) or (x_l,x_u)
!    that are closer than identical_bounds_tol will be reset to the average
!    of their values
!
!   indicator_tol_p. If %indicator_type = 1, a constraint/bound will be 
!    deemed to be active <=> distance to nearest bound <= %indicator_p_tol
!
!   indicator_tol_pd. If %indicator_type = 2, a constraint/bound will be 
!    deemed to be active <=> distance to nearest bound 
!        <= %indicator_tol_pd * size of corresponding multiplier
!
!   indicator_tol_tapia. If %indicator_type = 3, a constraint/bound will be 
!    deemed to be active <=> distance to nearest bound 
!        <= %indicator_tol_tapia * distance to same bound at previous iteration
! 
!   cpu_time_limit. The maximum CPU time allowed
!
!  LOGICAL control parameters:
!
!   remove_dependencies. If true, the equality constraints will be preprocessed
!    to remove any linear dependencies
!
!   treat_zero_bounds_as_general. If true, any problem bound with the value
!    zero will be treated as if it were a general value
!
!   just_feasible. If just_feasible is .TRUE., the algorithm will stop
!    as soon as a feasible point is found. Otherwise, the optimal solution
!    to the problem will be found
!
!   getdua. If getdua, is true, advanced initial values are obtained for the 
!    dual variables
!   
!   feasol. If feasol is true, the final solution obtained will be perturbed 
!    so that variables close to their bounds are moved onto these bounds
!
!   balance_initial_complentarity is .true. if the initial complemetarity
!    is required to be balanced
!
!   ues_corrector is .true. if a corrector step is to be combined with the
!    (predictor) search direction

!   array_syntax_worse_than_do_loop. If array_syntax_worse_than_do_loop is
!    true, f77-style do loops will be used rather than 
!    f90-style array syntax for vector operations
!
!  CHARACTER control parameters:
!
!  prefix (len=30). All output lines will be prefixed by 
!    %prefix(2:LEN(TRIM(%prefix))-1)
!   where %prefix contains the required string enclosed in 
!   quotes, e.g. "string" or 'string'
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      TYPE ( LSQP_data_type ), INTENT( OUT ) :: data
      TYPE ( LSQP_control_type ), INTENT( OUT ) :: control        

!  Initalize SILS components

      CALL SILS_initialize( FACTORS = data%FACTORS, control = data%CNTL )
      data%CNTL%ordering = 3
!57V2 data%CNTL%ordering = 2
!57V3 data%CNTL%ordering = 5
!57V2 data%CNTL%scaling = 0
!57V2 data%CNTL%static_tolerance = zero
!57V2 data%CNTL%static_level = zero
!     data%CNTL%thresh = 50

!  Set control parameters

!  Integer parameters

      control%error = 6
      control%out = 6
      control%print_level = 0
      control%maxit = 1000
      control%start_print = - 1
      control%stop_print = - 1
      control%factor = 0
      control%max_col = 35
      control%indmin = 1000
      control%valmin = 1000
      control%itref_max = 1
      control%infeas_max = 200
      control%muzero_fixed = 1
      control%restore_problem = 2
      control%indicator_type = 3

!  Real parameters

      control%infinity = ten ** 19
      control%stop_p = epsmch ** 0.33
      control%stop_c = epsmch ** 0.33
      control%stop_d = epsmch ** 0.33
!     control%prfeas = ten ** 3
!     control%dufeas = ten ** 3
      control%prfeas = one
      control%dufeas = one
      control%muzero = - one
!     control%pivot_tol = data%CNTL%u
      control%pivot_tol = epsmch ** 0.75
      control%pivot_tol_for_dependencies = half
      control%zero_pivot = epsmch ** 0.75
      control%identical_bounds_tol = epsmch
      control%reduce_infeas = one - point01
      control%potential_unbounded = - 10.0_wp
      control%indicator_tol_p = control%stop_p
      control%indicator_tol_pd = 1.0_wp
      control%indicator_tol_tapia = 0.9_wp
      control%cpu_time_limit = - one

!  Logical parameters

      control%remove_dependencies = .TRUE.
      control%treat_zero_bounds_as_general = .FALSE.
      control%just_feasible = .FALSE.
      control%getdua = .FALSE.
!     control%feasol = .TRUE.
      control%feasol = .FALSE.
      control%balance_initial_complentarity = .FALSE.
      control%use_corrector = .FALSE.
      control%array_syntax_worse_than_do_loop = .FALSE.

!  Character parameters

      control%prefix = '""                            '

      data%trans = 0 ; data%tried_to_remove_deps = .FALSE.
      data%save_structure = .TRUE.

      RETURN  

!  End of LSQP_initialize

      END SUBROUTINE LSQP_initialize

!-*-*-*-*-   L S Q P _ R E A D _ S P E C F I L E  S U B R O U T I N E   -*-*-*-

      SUBROUTINE LSQP_read_specfile( control, device, alt_specname )

!  Reads the content of a specification file, and performs the assignment of 
!  values associated with given keywords to the corresponding control parameters

!  The defauly values as given by LSQP_initialize could (roughly) 
!  have been set as:

! BEGIN LSQP SPECIFICATIONS (DEFAULT)
!  error-printout-device                             6
!  printout-device                                   6
!  print-level                                       0
!  maximum-number-of-iterations                      1000
!  start-print                                       -1
!  stop-print                                        -1
!  factorization-used                                0
!  maximum-column-nonzeros-in-schur-complement       35
!  initial-integer-workspace                         1000
!  initial-real-workspace                            1000
!  maximum-refinements                               1
!  maximum-poor-iterations-before-infeasible         200
!  barrier-fixed-until-iteration                     1
!  indicator-type-used                               3
!  restore-problem-on-output                         0
!  infinity-value                                    1.0D+19
!  primal-accuracy-required                          1.0D-5
!  dual-accuracy-required                            1.0D-5
!  complementary-slackness-accuracy-required         1.0D-5
!  mininum-initial-primal-feasibility                1000.0
!  mininum-initial-dual-feasibility                  1000.0
!  initial-barrier-parameter                         -1.0
!  poor-iteration-tolerance                          0.98
!  minimum-potential-before-unbounded                -10.0
!  pivot-tolerance-used                              1.0D-12
!  pivot-tolerance-used-for-dependencies             0.5
!  zero-pivot-tolerance                              1.0D-12
!  identical-bounds-tolerance                        1.0D-15
!  primal-indicator-tolerence                        1.0D-5
!  primal-dual-indicator-tolerence                   1.0
!  tapia-indicator-tolerence                         0.9
!  maximum-cpu-time-limit                            -1.0
!  remove-linear-dependencies                        T
!  treat-zero-bounds-as-general                      F
!  just-find-feasible-point                          F
!  fix-barrier-parameter-throughout                  F
!  use-corrector-step                                F
!  balance-initial-complentarity                     F
!  get-advanced-dual-variables                       F
!  move-final-solution-onto-bound                    F
!  array-syntax-worse-than-do-loop                   F 
! END LSQP SPECIFICATIONS (DEFAULT)

!  Dummy arguments

      TYPE ( LSQP_control_type ), INTENT( INOUT ) :: control        
      INTEGER, INTENT( IN ) :: device
      CHARACTER( LEN = 16 ), OPTIONAL :: alt_specname

!  Programming: Nick Gould and Ph. Toint, January 2002.

!  Local variables

      INTEGER, PARAMETER :: lspec = 41
      CHARACTER( LEN = 16 ), PARAMETER :: specname = 'LSQP            '
      TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec

!  Define the keywords

!  Integer key-words

      spec(  1 )%keyword = 'error-printout-device'
      spec(  2 )%keyword = 'printout-device'
      spec(  3 )%keyword = 'print-level' 
      spec(  4 )%keyword = 'maximum-number-of-iterations'
      spec(  5 )%keyword = 'start-print'
      spec(  6 )%keyword = 'stop-print'
      spec(  7 )%keyword = 'factorization-used'
      spec(  8 )%keyword = 'maximum-column-nonzeros-in-schur-complement'
      spec(  9 )%keyword = 'initial-integer-workspace'
      spec( 10 )%keyword = 'initial-real-workspace'
      spec( 11 )%keyword = 'maximum-refinements'
      spec( 12 )%keyword = 'maximum-poor-iterations-before-infeasible'
      spec( 41 )%keyword = 'barrier-fixed-until-iteration'
      spec( 13 )%keyword = 'restore-problem-on-output'
      spec( 36 )%keyword = 'indicator-type-used'

!  Real key-words

      spec( 14 )%keyword = 'infinity-value'
      spec( 15 )%keyword = 'primal-accuracy-required'
      spec( 16 )%keyword = 'dual-accuracy-required'
      spec( 17 )%keyword = 'complementary-slackness-accuracy-required'
      spec( 18 )%keyword = 'mininum-initial-primal-feasibility'
      spec( 19 )%keyword = 'mininum-initial-dual-feasibility'
      spec( 20 )%keyword = 'initial-barrier-parameter'
      spec( 21 )%keyword = 'poor-iteration-tolerance'
      spec( 22 )%keyword = 'minimum-potential-before-unbounded'
      spec( 23 )%keyword = 'pivot-tolerance-used'
      spec( 24 )%keyword = 'pivot-tolerance-used-for-dependencies'
      spec( 25 )%keyword = 'zero-pivot-tolerance'
      spec( 26 )%keyword = 'identical-bounds-tolerance'
      spec( 37 )%keyword = 'primal-indicator-tolerence'
      spec( 38 )%keyword = 'primal-dual-indicator-tolerence'
      spec( 39 )%keyword = 'tapia-indicator-tolerence'
      spec( 40 )%keyword = 'maximum-cpu-time-limit'

!  Logical key-words

      spec( 27 )%keyword = 'remove-linear-dependencies'
      spec( 28 )%keyword = 'treat-zero-bounds-as-general'
      spec( 29 )%keyword = 'just-find-feasible-point'
      spec( 30 )%keyword = 'get-advanced-dual-variables'
      spec( 31 )%keyword = 'move-final-solution-onto-bound'
      spec( 32 )%keyword = ''
      spec( 33 )%keyword = 'balance-initial-complentarity'
      spec( 34 )%keyword = 'use-corrector-step'
      spec( 35 )%keyword = 'array-syntax-worse-than-do-loop'

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
      CALL SPECFILE_assign_value( spec( 4 ), control%maxit,                    &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 5 ), control%start_print,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 6 ), control%stop_print,               &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 7 ), control%factor,                   &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 8 ), control%max_col,                  &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 9 ), control%indmin,                   &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 10 ), control%valmin,                  &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 11 ), control%itref_max,               &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 12 ), control%infeas_max,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 41 ), control%muzero_fixed,            &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 13 ), control%restore_problem,         &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 36 ), control%indicator_type,          &
                                  control%error )

!  Set real values


      CALL SPECFILE_assign_value( spec( 14 ), control%infinity,                &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 15 ), control%stop_p,                  &
                                  control%error )     
      CALL SPECFILE_assign_value( spec( 16 ), control%stop_d,                  &
                                  control%error )     
      CALL SPECFILE_assign_value( spec( 17 ), control%stop_c,                  &
                                  control%error )     
      CALL SPECFILE_assign_value( spec( 18 ), control%prfeas,                  &
                                  control%error )     
      CALL SPECFILE_assign_value( spec( 19 ), control%dufeas,                  &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 20 ), control%muzero,                  &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 21 ), control%reduce_infeas,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 22 ), control%potential_unbounded,     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 23 ), control%pivot_tol,               &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 24 ),                                  &
                                  control%pivot_tol_for_dependencies,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 25 ), control%zero_pivot,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 26 ), control%identical_bounds_tol,    &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 37 ), control%indicator_tol_p,         &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 38 ), control%indicator_tol_pd,        &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 39 ), control%indicator_tol_tapia,     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 40 ), control%cpu_time_limit,          &
                                  control%error )

!  Set logical values

      CALL SPECFILE_assign_value( spec( 27 ), control%remove_dependencies,     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 28 ),                                  &
                                  control%treat_zero_bounds_as_general,        &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 29 ), control%just_feasible,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 30 ), control%getdua,                  &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 31 ), control%feasol,                  &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 33 ),                                  &
                                  control%balance_initial_complentarity,       &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 34 ), control%use_corrector,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 35 ),                                  &
                                  control%array_syntax_worse_than_do_loop,     &
                                  control%error )

      RETURN

      END SUBROUTINE LSQP_read_specfile

!-*-*-*-*-*-*-*-*-*-   L S Q P _ S O L V E  S U B R O U T I N E   -*-*-*-*-*-*-*

      SUBROUTINE LSQP_solve( prob, data, control, inform, C_stat, B_stat )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Minimize the linear/separable objective
!
!        1/2 || W * ( x - x^0 ) ||_2^2 + g^T x + f  
!
!  where
!
!               (c_l)_i <= (Ax)_i <= (c_u)_i , i = 1, .... , m,
!
!  and        (x_l)_i <=   x_i  <= (x_u)_i , i = 1, .... , n,
!
!  where x is a vector of n components ( x_1, .... , x_n ),
!  A is an m by n matrix, and any of the bounds (c_l)_i, (c_u)_i
!  (x_l)_i, (x_u)_i may be infinite, using a primal-dual method.
!  The subroutine is particularly appropriate when A is sparse.
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!
!  prob is a structure of type QPT_problem_type, whose components hold 
!   information about the problem on input, and its solution on output.
!   The following components must be set:
!
!   %new_problem_structure is a LOGICAL variable, which must be set to 
!    .TRUE. by the user if this is the first problem with this "structure"
!    to be solved since the last call to LSQP_initialize, and .FALSE. if
!    a previous call to a problem with the same "structure" (but different
!    numerical data) was made.
!
!   %n is an INTEGER variable, which must be set by the user to the
!    number of optimization parameters, n.  RESTRICTION: %n >= 1
!                 
!   %m is an INTEGER variable, which must be set by the user to the
!    number of general linear constraints, m. RESTRICTION: %m >= 0
!                 
!   %gradient_kind is an INTEGER variable which defines the type of linear
!    term of the objective function to be used. Possible values are
!
!     0  the linear term g will be zero, and the analytic centre of the 
!        feasible region will be found if in addition %Hessian_kind is 0.
!        %G (see below) need not be set
!
!     1  each component of the linear terms g will be one. 
!        %G (see below) need not be set
!
!     any other value - the gradients will be those given by %G (see below)
!
!   %Hessian_kind is an INTEGER variable which defines the type of objective
!    function to be used. Possible values are
!
!     0  all the weights will be zero, and the analytic centre of the 
!        feasible region will be found. %WEIGHT (see below) need not be set
!
!     1  all the weights will be one. %WEIGHT (see below) need not be set
!
!     any other value - the weights will be those given by %WEIGHT (see below)
!
!   %WEIGHT is a REAL array, which need only be set if %Hessian_kind is not 0 
!    or 1. If this is so, it must be of length at least %n, and contain the
!    weights W for the objective function. 
!  
!   %X0 is a REAL array, which need only be set if %Hessian_kind is not 0.
!    If this is so, it must be of length at least %n, and contain the
!    weights X^0 for the objective function. 
!  
!   %A is a structure of type SMT_type used to hold the matrix A. 
!    Three storage formats are permitted:
!
!    i) sparse, co-ordinate
!
!       In this case, the following must be set:
!
!       A%type( 1 : 10 ) = TRANSFER( 'COORDINATE', A%type )
!       A%val( : )   the values of the components of A
!       A%row( : )   the row indices of the components of A
!       A%col( : )   the column indices of the components of A
!       A%ne         the number of nonzeros used to store A
!
!    ii) sparse, by rows
!
!       In this case, the following must be set:
!
!       A%type( 1 : 14 ) = TRANSFER( 'SPARSE_BY_ROWS', A%type )
!       A%val( : )   the values of the components of A, stored row by row
!       A%col( : )   the column indices of the components of A
!       A%ptr( : )   pointers to the start of each row, and past the end of
!                    the last row
!
!    iii) dense, by rows
!
!       In this case, the following must be set:
!
!       A%type( 1 : 5 ) = TRANSFER( 'DENSE', A%type )
!       A%val( : )   the values of the components of A, stored row by row,
!                    with each the entries in each row in order of 
!                    increasing column indicies.
!
!    On exit, the components will most likely have been reordered.
!    The output  matrix will be stored by rows, according to scheme (ii) above.
!    However, if scheme (i) is used for input, the output A%row will contain
!    the row numbers corresponding to the values in A%val, and thus in this
!    case the output matrix will be available in both formats (i) and (ii).   
! 
!   %G is a REAL array, which need only be set if %gradient_kind is not 0 
!    or 1. If this is so, it must be of length at least %n, and contain the
!    linear terms g for the objective function. 
!  
!   %f is a REAL variable, which must be set by the user to the value of
!    the constant term f in the objective function. On exit, it may have
!    been changed to reflect variables which have been fixed.
!
!   %C is a REAL array of length %m, which is used to store the values of 
!    A x. It need not be set on entry. On exit, it will have been filled 
!    with appropriate values.
!
!   %X is a REAL array of length %n, which must be set by the user
!    to estimaes of the solution, x. On successful exit, it will contain 
!    the required solution, x.
!
!   %C_l, %C_u are REAL arrays of length %n, which must be set by the user
!    to the values of the arrays c_l and c_u of lower and upper bounds on A x.
!    Any bound c_l_i or c_u_i larger than or equal to control%infinity in 
!    absolute value will be regarded as being infinite (see the entry 
!    control%infinity). Thus, an infinite lower bound may be specified by 
!    setting the appropriate component of %C_l to a value smaller than 
!    -control%infinity, while an infinite upper bound can be specified by 
!    setting the appropriate element of %C_u to a value larger than 
!    control%infinity. On exit, %C_l and %C_u will most likely have been 
!    reordered.
!   
!   %Y is a REAL array of length %m, which must be set by the user to
!    appropriate estimates of the values of the Lagrange multipliers 
!    corresponding to the general constraints c_l <= A x <= c_u. 
!    On successful exit, it will contain the required vector of Lagrange 
!    multipliers.
!
!   %X_l, %X_u are REAL arrays of length %n, which must be set by the user
!    to the values of the arrays x_l and x_u of lower and upper bounds on x.
!    Any bound x_l_i or x_u_i larger than or equal to control%infinity in 
!    absolute value will be regarded as being infinite (see the entry 
!    control%infinity). Thus, an infinite lower bound may be specified by 
!    setting the appropriate component of %X_l to a value smaller than 
!    -control%infinity, while an infinite upper bound can be specified by 
!    setting the appropriate element of %X_u to a value larger than 
!    control%infinity. On exit, %X_l and %X_u will most likely have been 
!    reordered.
!   
!   %Z is a REAL array of length %n, which must be set by the user to
!    appropriate estimates of the values of the dual variables 
!    (Lagrange multipliers corresponding to the simple bound constraints 
!    x_l <= x <= x_u). On successful exit, it will contain
!   the required vector of dual variables. 
!
!  data is a structure of type LSQP_data_type which holds private internal data
!
!  control is a structure of type LSQP_control_type that controls the 
!   execution of the subroutine and must be set by the user. Default values for
!   the elements may be set by a call to LSQP_initialize. See LSQP_initialize 
!   for details
!
!  inform is a structure of type LSQP_inform_type that provides 
!    information on exit from LSQP_solve. The component status 
!    has possible values:
!  
!     0 Normal termination with a locally optimal solution.
!
!    -1 An allocation error occured; the status is given in the component
!       alloc_status.
!
!    -2 A deallocation error occured; the status is given in the component
!       alloc_status.
!
!   - 3 one of the restrictions 
!        prob%n     >=  1
!        prob%m     >=  0
!        prob%A%type in { 'DENSE', 'SPARSE_BY_ROWS', 'COORDINATE' }
!       has been violated.
!
!    -4 The constraints are inconsistent.
!
!    -5 The constraints appear to have no feasible point.
!
!    -7 The objective function appears to be unbounded from below on the
!       feasible set.
!
!    -8 The analytic center appears to be unbounded.
!
!    -9 The analysis phase of the factorization failed; the return status 
!       from the factorization package is given in the component factor_status.
!      
!   -10 The factorization failed; the return status from the factorization
!       package is given in the component factor_status.
!      
!   -11 The solve of a required linear system failed; the return status from 
!       the factorization package is given in the component factor_status.
!      
!   -16 The problem is so ill-conditoned that further progress is impossible.  
!
!   -17 The step is too small to make further impact.
!
!   -18 Too many iterations have been performed. This may happen if
!       control%maxit is too small, but may also be symptomatic of 
!       a badly scaled problem.
!
!   -19 Too much CPU time has passed. This may happen if control%cpu_time_limit 
!       is too small, but may also be symptomatic of a badly scaled problem.
!
!  On exit from LSQP_solve, other components of inform give the 
!  following:
!
!     alloc_status = The status of the last attempted allocation/deallocation 
!     iter   = The total number of iterations required.
!     nbacts = The number of backtracks required in the linesearch.
!     factorization_integer = The total integer workspace required by the 
!              factorization.
!     factorization_real = The total real workspace required by the 
!              factorization.
!     nfacts = The total number of factorizations performed.
!     nbacts = The total number of "wasted" function evaluations during the 
!              linesearch.
!     factorization_status = the return status from the matrix factorization
!              package.   
!     obj = the value of the objective function ||W*(x-x^0)||_2.
!     potential = the value of the logarithmic potential function 
!                 sum -log(distance to constraint boundary)
!     non_negligible_pivot = the smallest pivot which was not judged to be
!       zero when detecting linearly dependent constraints
!     bad_alloc = the name of the array for which an allocation/deallocation
!       error ocurred
!     time%total = the total time spent in the package.
!     time%preprocess = the time spent preprocessing the problem.
!     time%find_dependent = the time spent detecting linear dependencies
!     time%analyse = the time spent analysing the required matrices prior to
!                  factorization.
!     time%factorize = the time spent factorizing the required matrices.
!     time%solve = the time spent computing the search direction.
!
!  C_stat is an optional INTEGER array of length m, which if present will be 
!   set on exit to indicate the likely ultimate status of the constraints. 
!   Possible values are 
!   C_stat( i ) < 0, the i-th constraint is likely in the active set, 
!                    on its lower bound, 
!               > 0, the i-th constraint is likely in the active set
!                    on its upper bound, and
!               = 0, the i-th constraint is likely not in the active set
!
!  B_stat is an optional  INTEGER array of length m, which if present will be 
!   set on exit to indicate the likely ultimate status of the simple bound 
!   constraints. Possible values are 
!   B_stat( i ) < 0, the i-th bound constraint is likely in the active set, 
!                    on its lower bound, 
!               > 0, the i-th bound constraint is likely in the active set
!                    on its upper bound, and
!               = 0, the i-th bound constraint is likely not in the active set
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      TYPE ( LSQP_data_type ), INTENT( INOUT ) :: data
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      TYPE ( LSQP_inform_type ), INTENT( OUT ) :: inform
      INTEGER, INTENT( OUT ), OPTIONAL, DIMENSION( prob%m ) :: C_stat
      INTEGER, INTENT( OUT ), OPTIONAL, DIMENSION( prob%n ) :: B_stat

!  Local variables

      INTEGER :: i, j, a_ne, n_depen, lbreak, lbnds, nzc
      INTEGER :: dy_l_lower, dy_l_upper, dy_u_lower, dy_u_upper
      INTEGER :: dz_l_lower, dz_l_upper, dz_u_lower, dz_u_upper
      REAL :: time, time_start, dum
      REAL ( KIND = wp ) :: fixed_sum, av_bnd
      LOGICAL :: reallocate, printi, remap_freed, reset_bnd, stat_required
      TYPE ( FDC_control_type ) :: FDC_control        
      TYPE ( FDC_inform_type ) :: FDC_inform

      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' entering LSQP_solve ' )" )

!  Initialize time

      CALL CPU_TIME( time_start )

!  Set initial timing breakdowns

      inform%time%total = 0.0 ; inform%time%analyse = 0.0
      inform%time%factorize = 0.0 ; inform%time%solve = 0.0
      inform%time%preprocess = 0.0 ; inform%time%find_dependent = 0.0

!  Initialize counts

      inform%status = GALAHAD_ok
      inform%alloc_status = 0 ; inform%bad_alloc = ''
      inform%factorization_status = 0
      inform%iter = - 1 ; inform%nfacts = - 1 ; inform%nbacts = 0
      inform%factorization_integer = - 1 ; inform%factorization_real = - 1 
      inform%obj = - one ; inform%potential = infinity
      inform%non_negligible_pivot = zero
      inform%feasible = .FALSE.
      stat_required = PRESENT( C_stat ) .AND. PRESENT( B_stat )

!  Basic single line of output per iteration

      printi = control%out > 0 .AND. control%print_level >= 1 

!  Ensure that input parameters are within allowed ranges

      IF ( prob%n < 1 .OR. prob%m < 0 .OR.                                     &
           .NOT. QPT_keyword_A( prob%A%type ) ) THEN
        inform%status = GALAHAD_error_restrictions
        IF ( control%error > 0 .AND. control%print_level > 0 )                 &
          WRITE( control%error, 2010 ) inform%status 
        CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
        GO TO 800 
      END IF 

!  If required, write out problem 

      IF ( control%out > 0 .AND. control%print_level >= 20 ) THEN
        WRITE( control%out, "( ' n, m = ', I0, 1X, I0 )" ) prob%n, prob%m
        WRITE( control%out, "( ' f = ', ES12.4 )" ) prob%f
        IF ( prob%gradient_kind == 0 ) THEN
          WRITE( control%out, "( ' G = zeros' )" )
        ELSE IF ( prob%gradient_kind == 1 ) THEN
          WRITE( control%out, "( ' G = ones' )" )
        ELSE
          WRITE( control%out, "( ' G = ', /, ( 5ES12.4 ) )" )                  &
            prob%G( : prob%n )
        END IF
        IF ( prob%Hessian_kind == 0 ) THEN
          WRITE( control%out, "( ' W = zeros' )" )
        ELSE IF ( prob%Hessian_kind == 1 ) THEN
          WRITE( control%out, "( ' W = ones ' )" )
          WRITE( control%out, "( ' X0 = ', /, ( 5ES12.4 ) )" )                 &
            prob%X0( : prob%n )
        ELSE
          WRITE( control%out, "( ' W = ', /, ( 5ES12.4 ) )" )                  &
            prob%WEIGHT( : prob%n )
          WRITE( control%out, "( ' X0 = ', /, ( 5ES12.4 ) )" )                 &
            prob%X0( : prob%n )
        END IF
        WRITE( control%out, "( ' X_l = ', /, ( 5ES12.4 ) )" )                  &
          prob%X_l( : prob%n )
        WRITE( control%out, "( ' X_u = ', /, ( 5ES12.4 ) )" )                  &
          prob%X_u( : prob%n )
        IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
          WRITE( control%out, "( ' A (dense) = ', /, ( 5ES12.4 ) )" )          &
            prob%A%val( : prob%n * prob%m )
        ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
          WRITE( control%out, "( ' A (row-wise) = ' )" )
          DO i = 1, prob%m
            WRITE( control%out, "( ( 2( 2I8, ES12.4 ) ) )" )                   &
              ( i, prob%A%col( j ), prob%A%val( j ),                           &
                j = prob%A%ptr( i ), prob%A%ptr( i + 1 ) - 1 )
          END DO
        ELSE
          WRITE( control%out, "( ' A (co-ordinate) = ' )" )
          WRITE( control%out, "( ( 2( 2I8, ES12.4 ) ) )" )                     &
          ( prob%A%row( i ), prob%A%col( i ), prob%A%val( i ), i = 1, prob%A%ne )
        END IF
        WRITE( control%out, "( ' C_l = ', /, ( 5ES12.4 ) )" )                  &
          prob%C_l( : prob%m )
        WRITE( control%out, "( ' C_u = ', /, ( 5ES12.4 ) )" )                  &
          prob%C_u( : prob%m )
      END IF

!  Check that problem bounds are consistent; reassign any pair of bounds
!  that are "essentially" the same

      reset_bnd = .FALSE.
      DO i = 1, prob%n
        IF ( prob%X_l( i ) - prob%X_u( i ) > control%identical_bounds_tol ) THEN
          inform%status = GALAHAD_error_bad_bounds
          IF ( control%error > 0 .AND. control%print_level > 0 )               &
            WRITE( control%error, 2010 ) inform%status 
          CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
          GO TO 800 
        ELSE IF ( prob%X_u( i ) == prob%X_l( i )  ) THEN
        ELSE IF ( prob%X_u( i ) - prob%X_l( i )                                &
                  <= control%identical_bounds_tol ) THEN
          av_bnd = half * ( prob%X_l( i ) + prob%X_u( i ) )
          prob%X_l( i ) = av_bnd ; prob%X_u( i ) = av_bnd
          reset_bnd = .TRUE.
        END IF
      END DO   
      IF ( reset_bnd .AND. printi ) WRITE( control%out,                        &
        "( ' ', /, '   **  Warning: one or more variable bounds reset ' )" )

      reset_bnd = .FALSE.
      DO i = 1, prob%m
        IF ( prob%C_l( i ) - prob%C_u( i ) > control%identical_bounds_tol ) THEN
          inform%status = GALAHAD_error_bad_bounds
          IF ( control%error > 0 .AND. control%print_level > 0 )               &
            WRITE( control%error, 2010 ) inform%status 
          CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
          GO TO 800 
        ELSE IF ( prob%C_u( i ) == prob%C_l( i ) ) THEN
        ELSE IF ( prob%C_u( i ) - prob%C_l( i )                                &
                  <= control%identical_bounds_tol ) THEN
          av_bnd = half * ( prob%C_l( i ) + prob%C_u( i ) )
          prob%C_l( i ) = av_bnd ; prob%C_u( i ) = av_bnd
          reset_bnd = .TRUE.
        END IF
      END DO   
      IF ( reset_bnd .AND. printi ) WRITE( control%out,                        &
        "( ' ', /, '   **  Warning: one or more constraint bounds reset ' )" )

!  Record the objective function value for any fixed variables

      fixed_sum = zero        
      IF ( prob%Hessian_kind == 1 ) THEN
        DO i = 1, prob%n
          IF ( prob%X_l( i ) == prob%X_u( i ) ) fixed_sum = fixed_sum +        &
               ( prob%X_l( i ) - prob%X0( i ) ) ** 2
        END DO
      ELSE IF ( prob%Hessian_kind == 2 ) THEN
        DO i = 1, prob%n
          IF ( prob%X_l( i ) == prob%X_u( i ) ) fixed_sum = fixed_sum +        &
               ( prob%WEIGHT( i ) * ( prob%X_l( i ) - prob%X0( i ) ) ) ** 2
        END DO
      END IF

!  Allocate sufficient workspace to hold a null Hessian (if needed)

      reallocate = .TRUE.
      IF ( ALLOCATED( data%RES_x ) ) THEN
        IF ( SIZE( data%RES_x ) < prob%n ) THEN ; DEALLOCATE( data%RES_x ) 
         ELSE ; reallocate = .FALSE. 
         END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%RES_x( prob%n ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%RES_x' ; GO TO 900
        END IF
      END IF
      
      reallocate = .TRUE.
      IF ( ALLOCATED( data%IW ) ) THEN
        IF ( SIZE( data%IW ) < prob%n + 1 ) THEN ; DEALLOCATE( data%IW )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%IW(  prob%n + 1 ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%IW' ; GO TO 900
        END IF
      END IF
      
!  ===========================
!  Preprocess the problem data
!  ===========================

      IF ( data%save_structure ) THEN
        data%new_problem_structure = prob%new_problem_structure
        data%save_structure = .FALSE.
      END IF

      IF ( prob%new_problem_structure ) THEN
        CALL QPP_initialize( data%QPP_map, data%QPP_control )
        data%QPP_control%infinity = control%infinity
        data%QPP_control%treat_zero_bounds_as_general =                       &
          control%treat_zero_bounds_as_general

!  Store the problem dimensions

        IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
          a_ne = prob%m * prob%n
        ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
          a_ne = prob%A%ptr( prob%m + 1 ) - 1
        ELSE
          a_ne = prob%A%ne 
        END IF

        IF ( printi ) WRITE( control%out,                                      &
               "( /, ' problem dimensions before preprocessing: ', /,          &
     &         ' n = ', I8, ' m = ', I8, ' a_ne = ', I8 )" )                   &
               prob%n, prob%m, a_ne

!  Perform the preprocessing. 

        CALL CPU_TIME( time )
        CALL QPP_reorder( data%QPP_map, data%QPP_control,                      &
                          data%QPP_inform, data%dims, prob,                    &
                          .FALSE., .FALSE., .FALSE. )
        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum
  
!  Test for satisfactory termination

        IF ( data%QPP_inform%status /= GALAHAD_ok ) THEN
          inform%status = data%QPP_inform%status
          IF ( control%out > 0 .AND. control%print_level >= 5 )                &
            WRITE( control%out, "( ' status ', I3, ' after QPP_reorder ')" )   &
             data%QPP_inform%status
          IF ( control%error > 0 .AND. control%print_level > 0 )               &
            WRITE( control%error, 2010 ) inform%status 
          CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )
          CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
          GO TO 800 
        END IF 

!  Record array lengths

        IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
          a_ne = prob%m * prob%n
        ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
          a_ne = prob%A%ptr( prob%m + 1 ) - 1
        ELSE
          a_ne = prob%A%ne 
        END IF

        IF ( printi ) WRITE( control%out,                                      &
               "( /, ' problem dimensions after preprocessing: ', /,           &
     &         ' n = ', I8, ' m = ', I8, ' a_ne = ', I8 )" )                   &
               prob%n, prob%m, a_ne

        prob%new_problem_structure = .FALSE.
        data%trans = 1

!  Recover the problem dimensions after preprocessing

      ELSE
        IF ( data%trans == 0 ) THEN
          CALL CPU_TIME( time )
          CALL QPP_apply( data%QPP_map, data%QPP_inform,                       &
                          prob, get_f = .TRUE., get_g = .TRUE.,                &
                          get_x_bounds = .TRUE., get_c_bounds = .TRUE.,        &
                          get_x = .TRUE., get_y = .TRUE., get_z = .TRUE.,      &
                          get_c = .TRUE., get_A = .TRUE. )
          CALL CPU_TIME( dum ) ; dum = dum - time
          inform%time%preprocess = inform%time%preprocess + dum

!  Test for satisfactory termination

          IF ( data%QPP_inform%status /= GALAHAD_ok ) THEN
            inform%status = data%QPP_inform%status
            IF ( control%out > 0 .AND. control%print_level >= 5 )              &
              WRITE( control%out, "( ' status ', I3, ' after QPP_apply ')" )   &
               data%QPP_inform%status
            IF ( control%error > 0 .AND. control%print_level > 0 )             &
              WRITE( control%error, 2010 ) inform%status 
            CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )
            CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
            GO TO 800 
          END IF 
        END IF 
        data%trans = data%trans + 1
      END IF

!  =================================================================
!  Check to see if the equality constraints are linearly independent
!  =================================================================

      IF ( prob%m > 0 .AND.                                                    &
           ( .NOT. data%tried_to_remove_deps .AND.                             &
              control%remove_dependencies ) ) THEN

        CALL CPU_TIME( time ) 

        IF ( control%out > 0 .AND. control%print_level >= 1 )                  &
          WRITE( control%out,                                                  &
            "( /, 1X, I0, ' equalities from ', I0, ' constraints ' )" )        &
            data%dims%c_equality, prob%m

!  Set control parameters

        CALL FDC_initialize( FDC_control )
        FDC_control%error = control%error
        FDC_control%out = control%out
        FDC_control%print_level = control%print_level
        FDC_control%indmin = control%indmin
        FDC_control%valmin = control%valmin
        FDC_control%zero_pivot = control%zero_pivot
        FDC_control%pivot_tol = control%pivot_tol_for_dependencies
        FDC_control%max_infeas = control%stop_p
        FDC_control%CNTL = data%CNTL
        FDC_control%prefix = '" - FDC:"                     '

!  Find any dependent rows

        nzc = prob%A%ptr( data%dims%c_equality + 1 ) - 1 
        CALL FDC_find_dependent( prob%n, data%dims%c_equality,                 &
                                 prob%A%val( : nzc ),                          &
                                 prob%A%col( : nzc ),                          &
                                 prob%A%ptr( : data%dims%c_equality + 1 ),     &
                                 prob%C_l, data%K, data%FACTORS,               &
                                 n_depen, data%Index_C_freed,                  &
                                 FDC_control, FDC_inform )
        inform%status = FDC_inform%status
        inform%non_negligible_pivot = FDC_inform%non_negligible_pivot

!  Record output parameters

        inform%alloc_status = FDC_inform%alloc_status
        inform%factorization_status = FDC_inform%factorization_status
        inform%factorization_integer = FDC_inform%factorization_integer
        inform%factorization_real = FDC_inform%factorization_real
        inform%bad_alloc = FDC_inform%bad_alloc
        inform%time%find_dependent = FDC_inform%time%total
        inform%time%analyse = inform%time%analyse + FDC_inform%time%analyse
        inform%time%factorize = inform%time%factorize + FDC_inform%time%factorize
        inform%nfacts = 1

        CALL CPU_TIME( dum )
        IF ( control%cpu_time_limit >= zero .AND.                              &
             dum - time_start > control%cpu_time_limit ) THEN
          inform%status = GALAHAD_error_cpu_limit
          IF ( control%error > 0 .AND. control%print_level > 0 )               &
            WRITE( control%error, 2010 ) inform%status 
          inform%time%total = dum - time_start ; GO TO 800 
        END IF 
        inform%time%preprocess = inform%time%preprocess + dum - time
        IF ( printi .AND. inform%non_negligible_pivot <                        &
             thousand * control%zero_pivot ) WRITE( control%out, "(            &
       &  /, 1X, 26 ( '*' ), ' WARNING ', 26 ( '*' ), /,                       & 
       &  ' ***  smallest allowed pivot =', ES11.4,' may be too small ***',    &
       &  /, ' ***  perhaps increase control%zero_pivot from', ES11.4,'  ***', &
       &  /, 1X, 26 ( '*' ), ' WARNING ', 26 ( '*' ) )" )                      &
           inform%non_negligible_pivot, control%zero_pivot

!  Check for error exits

        IF ( inform%status /= 0 ) THEN

!  Print details of the error exit

          IF ( control%error > 0 .AND. control%print_level >= 1 ) THEN
            WRITE( control%out, "( ' ' )" )
            IF ( inform%status /= GALAHAD_ok )                                 &
              WRITE( control%error, 2040 ) inform%status, 'LSQP_dependent'
          END IF
          GO TO 700
        END IF

        IF ( control%out > 0 .AND. control%print_level >= 2 .AND. n_depen > 0 )&
          WRITE( control%out, "(/, ' The following ',I7,' constraints appear', &
       &         ' to be dependent', /, ( 8I8 ) )" ) n_depen, data%Index_C_freed

        remap_freed = n_depen > 0 .AND. prob%n > 0

!  Special case: no free variables

        IF ( prob%n == 0 ) THEN
          prob%Y( : prob%m ) = zero
          prob%Z( : prob%n ) = zero
          prob%C( : prob%m ) = zero
          CALL LSQP_AX( prob%m, prob%C( : prob%m ), prob%m,                    &
                        prob%A%ptr( prob%m + 1 ) - 1, prob%A%val,              &
                        prob%A%col, prob%A%ptr, prob%n, prob%X, '+ ')
          GO TO 700
        END IF
        data%tried_to_remove_deps = .TRUE.
      ELSE
        remap_freed = .FALSE.
      END IF

      IF ( remap_freed ) THEN

!  Some of the current constraints will be removed by freeing them

        IF ( control%error > 0 .AND. control%print_level >= 1 )                &
          WRITE( control%out, "( /, ' -> ', I0, ' constraints are',            &
         & ' dependent and will be temporarily removed' )" ) n_depen

!  Allocate arrays to indicate which constraints have been freed

        reallocate = .TRUE.
        IF ( ALLOCATED( data%C_freed ) ) THEN
          IF ( SIZE( data%C_freed ) < n_depen ) THEN
            DEALLOCATE( data%C_freed ) ; ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%C_freed( n_depen ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%C_freed' ; GO TO 900
          END IF
        END IF
        
!  Free the constraint bounds as required

        DO i = 1, n_depen
          j = data%Index_c_freed( i )
          data%C_freed( i ) = prob%C_l( j )
          prob%C_l( j ) = - control%infinity
          prob%C_u( j ) = control%infinity
          prob%Y( j ) = zero
        END DO

        CALL QPP_initialize( data%QPP_map_freed, data%QPP_control )
        data%QPP_control%infinity = control%infinity
        data%QPP_control%treat_zero_bounds_as_general =                        &
          control%treat_zero_bounds_as_general

!  Store the problem dimensions

        data%dims_save_freed = data%dims
        a_ne = prob%A%ne 

        IF ( printi ) WRITE( control%out,                                      &
               "( /, ' problem dimensions before removal of dependecies: ', /, &
     &         ' n = ', I8, ' m = ', I8, ' a_ne = ', I8 )" )                   &
               prob%n, prob%m, a_ne

!  Perform the preprocessing

        CALL CPU_TIME( time ) 
        CALL QPP_reorder( data%QPP_map_freed, data%QPP_control,                &
                          data%QPP_inform, data%dims, prob,                    &
                          .FALSE., .FALSE., .FALSE. )
        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum
  
        data%dims%nc = data%dims%c_u_end - data%dims%c_l_start + 1 
        data%dims%x_s = 1 ; data%dims%x_e = prob%n
        data%dims%c_s = data%dims%x_e + 1 
        data%dims%c_e = data%dims%x_e + data%dims%nc  
        data%dims%c_b = data%dims%c_e - prob%m
        data%dims%y_s = data%dims%c_e + 1 
        data%dims%y_e = data%dims%c_e + prob%m
        data%dims%y_i = data%dims%c_s + prob%m 
        data%dims%v_e = data%dims%y_e

!  Test for satisfactory termination

        IF ( data%QPP_inform%status /= GALAHAD_ok ) THEN
          inform%status = data%QPP_inform%status
          IF ( control%out > 0 .AND. control%print_level >= 5 )                &
            WRITE( control%out, "( ' status ', I3, ' after QPP_reorder ')" )   &
             data%QPP_inform%status
          IF ( control%error > 0 .AND. control%print_level > 0 )               &
            WRITE( control%error, 2010 ) inform%status 
          CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
          CALL QPP_terminate( data%QPP_map_freed, data%QPP_control,            &
                              data%QPP_inform )
          CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )
          GO TO 800 
        END IF 

!  Record revised array lengths

        IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
          a_ne = prob%m * prob%n
        ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
          a_ne = prob%A%ptr( prob%m + 1 ) - 1
        ELSE
          a_ne = prob%A%ne 
        END IF

        IF ( printi ) WRITE( control%out,                                      &
               "( /, ' problem dimensions after removal of dependencies: ', /, &
     &         ' n = ', I8, ' m = ', I8, ' a_ne = ', I8 )" )                   &
               prob%n, prob%m, a_ne

      END IF

!  Compute the dimension of the KKT system 

      data%dims%nc = data%dims%c_u_end - data%dims%c_l_start + 1 

!  Arrays containing data relating to the composite vector ( x  c  y )
!  are partitioned as follows:

!   <---------- n --------->  <---- nc ------>  <-------- m --------->
!                             <-------- m --------->
!                        <-------- m --------->
!   -------------------------------------------------------------------
!   |                   |    |                 |    |                 |
!   |         x              |       c         |          y           |
!   |                   |    |                 |    |                 |
!   -------------------------------------------------------------------
!    ^                 ^    ^ ^               ^ ^    ^               ^
!    |                 |    | |               | |    |               |
!   x_s                |    |c_s              |y_s  y_i             y_e = v_e
!                      |    |                 |                     
!                     c_b  x_e               c_e

      data%dims%x_s = 1 ; data%dims%x_e = prob%n
      data%dims%c_s = data%dims%x_e + 1 
      data%dims%c_e = data%dims%x_e + data%dims%nc  
      data%dims%c_b = data%dims%c_e - prob%m
      data%dims%y_s = data%dims%c_e + 1 
      data%dims%y_e = data%dims%c_e + prob%m
      data%dims%y_i = data%dims%c_s + prob%m 
      data%dims%v_e = data%dims%y_e

!  Allocate real workspace

      reallocate = .TRUE.
      IF ( ALLOCATED( data%SOL_y ) ) THEN
        IF ( SIZE( data%SOL_y ) < prob%m ) THEN ; DEALLOCATE( data%SOL_y )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%SOL_y( prob%m ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%SOL_y' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%RES_y ) ) THEN
        IF ( SIZE( data%RES_y ) < prob%m ) THEN ; DEALLOCATE( data%RES_y )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%RES_y( prob%m ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%RES_y' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%BEST_y ) ) THEN
        IF ( SIZE( data%BEST_y ) < prob%m ) THEN 
          DEALLOCATE( data%BEST_y ) ; ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%BEST_y( prob%m ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%BEST_y' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%SOL ) ) THEN
        IF ( SIZE( data%SOL ) < data%dims%v_e ) THEN ; DEALLOCATE( data%SOL )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%SOL( data%dims%v_e ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%SOL' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%RES ) ) THEN
        IF ( SIZE( data%RES ) < data%dims%v_e ) THEN ; DEALLOCATE( data%RES )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%RES( data%dims%v_e ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%RES' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%BEST ) ) THEN
        IF ( SIZE( data%BEST ) < data%dims%v_e ) THEN ; DEALLOCATE( data%BEST )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%BEST( data%dims%v_e ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%BEST' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%HX ) ) THEN
        IF ( SIZE( data%HX ) /= data%dims%v_e ) THEN ; DEALLOCATE( data%HX )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%HX( data%dims%v_e ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%HX' ; GO TO 900 ; END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%GRAD_L ) ) THEN
        IF ( SIZE( data%GRAD_L ) /= data%dims%c_e ) THEN 
          DEALLOCATE( data%GRAD_L ) ; ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%GRAD_L( data%dims%c_e ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN ; 
          inform%bad_alloc = 'lsqp: data%GRAD_L' ; GO TO 900 ; END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%DIST_X_l ) ) THEN
        IF ( LBOUND( data%DIST_X_l, 1 ) /= data%dims%x_l_start .OR.            &
             UBOUND( data%DIST_X_l, 1 ) /= data%dims%x_l_end ) THEN 
          DEALLOCATE( data%DIST_X_l )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%DIST_X_l( data%dims%x_l_start : data%dims%x_l_end ),    &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%DIST_X_l' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%DIST_X_u ) ) THEN
        IF ( LBOUND( data%DIST_X_u, 1 ) /= data%dims%x_u_start .OR.            &
             UBOUND( data%DIST_X_u, 1 ) /= data%dims%x_u_end ) THEN 
          DEALLOCATE( data%DIST_X_u )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%DIST_X_u( data%dims%x_u_start : data%dims%x_u_end ),    &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%DIST_X_u' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%Z_l ) ) THEN
        IF ( LBOUND( data%Z_l, 1 ) /= data%dims%x_free + 1 .OR.                &
             UBOUND( data%Z_l, 1 ) /= data%dims%x_l_end ) THEN 
          DEALLOCATE( data%Z_l )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%Z_l( data%dims%x_free + 1 : data%dims%x_l_end ),        &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%Z_l' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%Z_u ) ) THEN
        IF ( LBOUND( data%Z_u, 1 ) /= data%dims%x_u_start .OR.                 &
             UBOUND( data%Z_u, 1 ) /= prob%n ) THEN 
          DEALLOCATE( data%Z_u ) ; ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%Z_u( data%dims%x_u_start : prob%n ),               &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%Z_u' ; GO TO 900
        END IF
      END IF
      
      reallocate = .TRUE.
      IF ( ALLOCATED( data%BARRIER_X ) ) THEN
        IF ( LBOUND( data%BARRIER_X, 1 ) /= data%dims%x_free + 1 .OR.          &
             UBOUND( data%BARRIER_X, 1 ) /= prob%n ) THEN 
          DEALLOCATE( data%BARRIER_X )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%BARRIER_X( data%dims%x_free + 1 : prob%n ),        &
                  STAT = inform%alloc_status)
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%BARRIER_X' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%Y_l ) ) THEN
        IF ( LBOUND( data%Y_l, 1 ) /= data%dims%c_l_start .OR.                 &
             UBOUND( data%Y_l, 1 ) /= data%dims%c_l_end ) THEN
          DEALLOCATE( data%Y_l )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%Y_l( data%dims%c_l_start : data%dims%c_l_end ),         &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%Y_l' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%DY_l ) ) THEN
        IF ( LBOUND( data%DY_l, 1 ) /= data%dims%c_l_start .OR.                &
             UBOUND( data%DY_l, 1 ) /= data%dims%c_l_end ) THEN  
          DEALLOCATE( data%DY_l )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%DY_l( data%dims%c_l_start : data%dims%c_l_end ),        &
                  STAT = inform%alloc_status)
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%DY_l' ; GO TO 900 ; END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%DIST_C_l ) ) THEN
        IF ( LBOUND( data%DIST_C_l, 1 ) /= data%dims%c_l_start .OR.            &
             UBOUND( data%DIST_C_l, 1 ) /= data%dims%c_l_end ) THEN  
          DEALLOCATE( data%DIST_C_l )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%DIST_C_l( data%dims%c_l_start : data%dims%c_l_end ),    &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%DIST_C_l' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%Y_u ) ) THEN
        IF ( LBOUND( data%Y_u, 1 ) /= data%dims%c_u_start .OR.                 &
             UBOUND( data%Y_u, 1 ) /= data%dims%c_u_end ) THEN
          DEALLOCATE( data%Y_u )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%Y_u( data%dims%c_u_start : data%dims%c_u_end ),         &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%Y_u' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%DY_u ) ) THEN
        IF ( LBOUND( data%DY_u, 1 ) /= data%dims%c_u_start .OR.                &
             UBOUND( data%DY_u, 1 ) /= data%dims%c_u_end ) THEN  
          DEALLOCATE( data%DY_u )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%DY_u( data%dims%c_u_start : data%dims%c_u_end ),        &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%DY_u' ; GO TO 900 ; END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%DIST_C_u ) ) THEN
        IF ( LBOUND( data%DIST_C_u, 1 ) /= data%dims%c_u_start .OR.            &
             UBOUND( data%DIST_C_u, 1 ) /= data%dims%c_u_end ) THEN  
          DEALLOCATE( data%DIST_C_u )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%DIST_C_u( data%dims%c_u_start : data%dims%c_u_end ),    &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%DIST_C_u' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%C ) ) THEN
        IF ( LBOUND( data%C, 1 ) /= data%dims%c_l_start .OR.                   &
             UBOUND( data%C, 1 ) /= data%dims%c_u_end ) THEN 
          DEALLOCATE( data%C ) ; ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%C( data%dims%c_l_start : data%dims%c_u_end ),           &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%bad_alloc = 'lsqp: data%C' ; GO TO 900 ; END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%BARRIER_C ) ) THEN
        IF ( LBOUND( data%BARRIER_C, 1 ) /= data%dims%c_l_start .OR.           &
             UBOUND( data%BARRIER_C, 1 ) /= data%dims%c_u_end ) THEN 
          DEALLOCATE( data%BARRIER_C )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%BARRIER_C( data%dims%c_l_start : data%dims%c_u_end ),   &
                  STAT = inform%alloc_status)
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%BARRIER_C' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%SCALE_C ) ) THEN
        IF ( LBOUND( data%SCALE_C, 1 ) /= data%dims%c_l_start .OR.             &
             UBOUND( data%SCALE_C, 1 ) /= data%dims%c_u_end ) THEN 
          DEALLOCATE( data%SCALE_C )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%SCALE_C( data%dims%c_l_start : data%dims%c_u_end ),     &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%SCALE_C' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%DELTA ) ) THEN
        IF ( SIZE( data%DELTA ) /= data%dims%v_e ) THEN 
          DEALLOCATE( data%DELTA ) ; ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%DELTA( data%dims%v_e ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%DELTA' ; GO TO 900 ; END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%RHS ) ) THEN
        IF ( SIZE( data%RHS ) /= data%dims%v_e ) THEN ; DEALLOCATE( data%RHS )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%RHS( data%dims%v_e ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%RHS' ; GO TO 900 ; END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%DZ_l ) ) THEN
        IF ( LBOUND( data%DZ_l, 1 ) /= data%dims%x_free + 1 .OR.               &
             UBOUND( data%DZ_l, 1 ) /= data%dims%x_l_end ) THEN 
          DEALLOCATE( data%DZ_l )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%DZ_l( data%dims%x_free + 1 : data%dims%x_l_end ),       &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%DZ_l' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%DZ_u ) ) THEN
        IF ( LBOUND( data%DZ_u, 1 ) /= data%dims%x_u_start .OR.                &
             UBOUND( data%DZ_u, 1 ) /= prob%n ) THEN 
          DEALLOCATE( data%DZ_u ) ; ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%DZ_u( data%dims%x_u_start : prob%n ),                   &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%DZ_u' ; GO TO 900
        END IF
      END IF

!  Allocate optional extra arrays

      lbnds = 0
      lbreak = 0

      reallocate = .TRUE.
      IF ( ALLOCATED( data%COEF0 ) ) THEN
        IF ( SIZE( data%COEF0 ) < lbnds ) THEN ; DEALLOCATE( data%COEF0 )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%COEF0( lbnds ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%COEF0' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%COEF1 ) ) THEN
        IF ( SIZE( data%COEF1 ) < lbnds ) THEN ; DEALLOCATE( data%COEF1 )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%COEF1( lbnds ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%COEF1' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%COEF2 ) ) THEN
        IF ( SIZE( data%COEF2 ) < lbnds ) THEN ; DEALLOCATE( data%COEF2 )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%COEF2( lbnds ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%COEF2' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%COEF3 ) ) THEN
        IF ( SIZE( data%COEF3 ) < lbnds ) THEN ; DEALLOCATE( data%COEF3 )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%COEF3( lbnds ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%COEF3' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%COEF4 ) ) THEN
        IF ( SIZE( data%COEF4 ) < lbnds ) THEN ; DEALLOCATE( data%COEF4 )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%COEF4( lbnds ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%COEF4' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%BREAKP ) ) THEN
        IF ( SIZE( data%BREAKP ) < lbreak ) THEN ; DEALLOCATE( data%BREAKP )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%BREAKP( lbreak ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%BREAKP' ; GO TO 900
        END IF
      END IF

      reallocate = .TRUE.
      IF ( ALLOCATED( data%IBREAK ) ) THEN
        IF ( SIZE( data%IBREAK ) < lbreak ) THEN ; DEALLOCATE( data%IBREAK )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( data%IBREAK( lbreak ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          inform%bad_alloc = 'lsqp: data%IBREAK' ; GO TO 900
        END IF
      END IF

      IF ( control%use_corrector ) THEN
        reallocate = .TRUE.
        IF ( ALLOCATED( data%DELTA_cor ) ) THEN
          IF ( SIZE( data%DELTA_cor ) /= data%dims%v_e ) THEN 
            DEALLOCATE( data%DELTA_cor ) ; ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%DELTA_cor( data%dims%v_e ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DELTA_cor' ; GO TO 900 ; END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( data%DY_cor_l ) ) THEN
          IF ( LBOUND( data%DY_cor_l, 1 ) /= data%dims%c_l_start .OR.           &
               UBOUND( data%DY_cor_l, 1 ) /= data%dims%c_l_end ) THEN  
            DEALLOCATE( data%DY_cor_l )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%DY_cor_l( data%dims%c_l_start : data%dims%c_l_end ),   &
                    STAT = inform%alloc_status)
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DY_cor_l' ; GO TO 900 ; END IF
        END IF
        dy_l_lower = data%dims%c_l_start
        dy_l_upper = data%dims%c_l_end

        reallocate = .TRUE.
        IF ( ALLOCATED( data%DY_cor_u ) ) THEN
          IF ( LBOUND( data%DY_cor_u, 1 ) /= data%dims%c_u_start .OR.           &
               UBOUND( data%DY_cor_u, 1 ) /= data%dims%c_u_end ) THEN  
            DEALLOCATE( data%DY_cor_u )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%DY_cor_u( data%dims%c_u_start : data%dims%c_u_end ),   &
                    STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DY_cor_u' ; GO TO 900 ; END IF
        END IF
        dy_u_lower = data%dims%c_u_start
        dy_u_upper = data%dims%c_u_end

        reallocate = .TRUE.
        IF ( ALLOCATED( data%DZ_cor_l ) ) THEN
          IF ( LBOUND( data%DZ_cor_l, 1 ) /= data%dims%x_free + 1 .OR.          &
               UBOUND( data%DZ_cor_l, 1 ) /= data%dims%x_l_end ) THEN 
            DEALLOCATE( data%DZ_cor_l )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%DZ_cor_l( data%dims%x_free + 1 : data%dims%x_l_end ),  &
                    STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DZ_cor_l' ; GO TO 900
          END IF
        END IF
        dz_l_lower = data%dims%x_free + 1
        dz_l_upper = data%dims%x_l_end

        reallocate = .TRUE.
        IF ( ALLOCATED( data%DZ_cor_u ) ) THEN
          IF ( LBOUND( data%DZ_cor_u, 1 ) /= data%dims%x_u_start .OR.           &
               UBOUND( data%DZ_cor_u, 1 ) /= prob%n ) THEN 
            DEALLOCATE( data%DZ_cor_u ) ; ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%DZ_cor_u( data%dims%x_u_start : prob%n ),             &
                    STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DZ_cor_u' ; GO TO 900
          END IF
        END IF
        dz_u_lower = data%dims%x_u_start
        dz_u_upper = prob%n

      ELSE
        IF ( .NOT. ALLOCATED( data%DELTA_cor ) ) THEN
          ALLOCATE( data%DELTA_cor( 0 ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DELTA_cor' ; GO TO 900 ; END IF
        END IF

        IF ( .NOT. ALLOCATED( data%DY_cor_l ) ) THEN
          ALLOCATE( data%DY_cor_l( 0 ),  STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DY_cor_l' ; GO TO 900
          END IF
        END IF
        dy_l_lower = 1
        dy_l_upper = 0

        IF ( .NOT. ALLOCATED( data%DY_cor_u ) ) THEN
          ALLOCATE( data%DY_cor_u( 0 ),  STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DY_cor_u' ; GO TO 900
          END IF
        END IF
        dy_u_lower = 1
        dy_u_upper = 0

        IF ( .NOT. ALLOCATED( data%DZ_cor_l ) ) THEN
          ALLOCATE( data%DZ_cor_l( 0 ),  STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DZ_cor_l' ; GO TO 900
          END IF
        END IF
        dz_l_lower = 1
        dz_l_upper = 0

        IF ( .NOT. ALLOCATED( data%DZ_cor_u ) ) THEN
          ALLOCATE( data%DZ_cor_u( 0 ),  STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DZ_cor_u' ; GO TO 900
          END IF
        END IF
        dz_u_lower = 1
        dz_u_upper = 0
      END IF

      IF ( stat_required ) THEN
        reallocate = .TRUE.
        IF ( ALLOCATED( data%H_s ) ) THEN
          IF ( SIZE( data%H_s ) < prob%n ) THEN ; DEALLOCATE( data%H_s )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%H_s( prob%n ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%H_s' ; GO TO 900
          END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( data%A_s ) ) THEN
          IF ( SIZE( data%A_s ) < prob%m ) THEN ; DEALLOCATE( data%A_s )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%A_s( prob%m ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%A_s' ; GO TO 900
          END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( data%Y_last ) ) THEN
          IF ( SIZE( data%Y_last ) < prob%m ) THEN ; DEALLOCATE( data%Y_last )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%Y_last( prob%m ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%Y_last' ; GO TO 900
          END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( data%Z_last ) ) THEN
          IF ( SIZE( data%Z_last ) < prob%n ) THEN ; DEALLOCATE( data%Z_last )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( data%Z_last( prob%n ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%Z_last' ; GO TO 900
          END IF
        END IF
      END IF

!  =================
!  Solve the problem
!  =================

!  constraint/variable exit ststus required

      IF ( stat_required ) THEN
        IF ( prob%Hessian_kind == 0 ) THEN
          IF ( prob%gradient_kind == 0 .OR. prob%gradient_kind == 1 ) THEN
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  C_last = data%A_s, X_last = data%H_s,        &
                                  Y_last = data%Y_last, Z_last = data%Z_last,  &
                                  C_stat = C_stat, B_Stat = B_Stat )
          ELSE
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  G = prob%G,                                  &
                                  C_last = data%A_s, X_last = data%H_s,        &
                                  Y_last = data%Y_last, Z_last = data%Z_last,  &
                                  C_stat = C_stat, B_Stat = B_Stat )
          END IF
        ELSE IF ( prob%Hessian_kind == 1 ) THEN
          IF ( prob%gradient_kind == 0 .OR. prob%gradient_kind == 1 ) THEN
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  X0 = prob%X0,                                &
                                  C_last = data%A_s, X_last = data%H_s,        &
                                  Y_last = data%Y_last, Z_last = data%Z_last,  &
                                  C_stat = C_stat, B_Stat = B_Stat )
          ELSE
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  X0 = prob%X0, G = prob%G,                    &
                                  C_last = data%A_s, X_last = data%H_s,        &
                                  Y_last = data%Y_last, Z_last = data%Z_last,  &
                                  C_stat = C_stat, B_Stat = B_Stat )
          END IF
        ELSE
          IF ( prob%gradient_kind == 0 .OR. prob%gradient_kind == 1 ) THEN
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  WEIGHT = prob%WEIGHT, X0 = prob%X0,          &
                                  C_last = data%A_s, X_last = data%H_s,        &
                                  Y_last = data%Y_last, Z_last = data%Z_last,  &
                                  C_stat = C_stat, B_Stat = B_Stat )
          ELSE
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  WEIGHT = prob%WEIGHT, X0 = prob%X0,          &
                                  G = prob%G,                                  &
                                  C_last = data%A_s, X_last = data%H_s,        &
                                  Y_last = data%Y_last, Z_last = data%Z_last,  &
                                  C_stat = C_stat, B_Stat = B_Stat )
          END IF
        END IF  

!  constraint/variable exit ststus not required

      ELSE
        IF ( prob%Hessian_kind == 0 ) THEN
          IF ( prob%gradient_kind == 0 .OR. prob%gradient_kind == 1 ) THEN
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind )
          ELSE
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  G = prob%G )
          END IF
        ELSE IF ( prob%Hessian_kind == 1 ) THEN
          IF ( prob%gradient_kind == 0 .OR. prob%gradient_kind == 1 ) THEN
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  X0 = prob%X0 )
          ELSE
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  X0 = prob%X0, G = prob%G )
          END IF
        ELSE
          IF ( prob%gradient_kind == 0 .OR. prob%gradient_kind == 1 ) THEN
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  WEIGHT = prob%WEIGHT, X0 = prob%X0 )
          ELSE
            CALL LSQP_solve_main( data%dims, prob%n, prob%m,                   &
                                  prob%A%val, prob%A%col, prob%A%ptr,          &
                                  prob%C_l, prob%C_u, prob%X_l, prob%X_u,      &
                                  prob%C, prob%X, prob%Y, prob%Z, data%RES_x,  &
                                  data%SOL_y, data%RES_y, data%BEST_y,         &
                                  data%SOL, data%RES, data%BEST, data%HX,      &
                                  data%GRAD_L, data%DIST_X_l, data%DIST_X_u,   &
                                  data%Z_l, data%Z_u, data%BARRIER_X,          &
                                  data%Y_l, data%DY_l,                         &
                                  data%DIST_C_l, data%Y_u, data%DY_u,          &
                                  data%DIST_C_u, data%C, data%BARRIER_C,       &
                                  data%SCALE_C, data%DELTA, data%RHS,          &
                                  data%DZ_l, data%DZ_u,                        &
                                  data%Abycol_val, data%DIAG_X, data%DIAG_C,   &
                                  data%IW, data%K_colptr, data%Abycol_ptr,     &
                                  data%Abycol_row, data%K, data%FACTORS,       &
                                  prob%f, data%COEF0, data%COEF1, data%COEF2,  &
                                  data%COEF3, data%COEF4, data%DELTA_cor,      &
                                  data%DY_cor_l, dy_l_lower, dy_l_upper,       &
                                  data%DY_cor_u, dy_u_lower, dy_u_upper,       &
                                  data%DZ_cor_l, dz_l_lower, dz_l_upper,       &
                                  data%DZ_cor_u, dz_u_lower, dz_u_upper,       &
                                  data%BREAKP, data%IBREAK, data%CNTL,         &
                                  control, inform,                             &
                                  prob%Hessian_kind, prob%gradient_kind,       &
                                  WEIGHT = prob%WEIGHT, X0 = prob%X0,          &
                                  G = prob%G )
          END IF
        END IF  
      END IF  

      inform%obj = inform%obj + half * fixed_sum

!     write(6,*) ' c_stat ', C_stat( : prob%m )

!  If some of the constraints were freed during the computation, refix them now

      IF ( remap_freed ) THEN

        CALL CPU_TIME( time )
        IF ( stat_required ) THEN
          C_stat( prob%m + 1 : data%QPP_map_freed%m ) = 0
          CALL SORT_inverse_permute( data%QPP_map_freed%m,                     &
                                     data%QPP_map_freed%c_map,                 &
                                     IX = C_stat( : data%QPP_map_freed%m ) )
          B_stat( prob%n + 1 : data%QPP_map_freed%n ) = - 1
          CALL SORT_inverse_permute( data%QPP_map_freed%n,                     &
                                     data%QPP_map_freed%x_map,                 &
                                     IX = B_stat( : data%QPP_map_freed%n ) )
        END IF
        CALL QPP_restore( data%QPP_map_freed, data%QPP_inform, prob,           &
                          get_all = .TRUE.)
        CALL QPP_terminate( data%QPP_map_freed, data%QPP_control,              &
                            data%QPP_inform )
        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum
        data%dims = data%dims_save_freed

!  Fix the temporarily freed constraint bounds

        DO i = 1, n_depen
          j = data%Index_c_freed( i )
          prob%C_l( j ) = data%C_freed( i )
          prob%C_u( j ) = data%C_freed( i )
        END DO
      END IF
      data%tried_to_remove_deps = .FALSE.

!  Retore the problem to its original form

  700 CONTINUE 
      data%trans = data%trans - 1
      IF ( data%trans == 0 ) THEN
        data%IW( : prob%n + 1 ) = 0
        data%RES_x( : prob%n ) = zero
        CALL CPU_TIME( time )
        IF ( stat_required ) THEN
          C_stat( prob%m + 1 : data%QPP_map%m ) = 0
          CALL SORT_inverse_permute( data%QPP_map%m, data%QPP_map%c_map,       &
                                     IX = C_stat( : data%QPP_map%m ) )
          B_stat( prob%n + 1 : data%QPP_map%n ) = - 1
          CALL SORT_inverse_permute( data%QPP_map%n, data%QPP_map%x_map,       &
                                     IX = B_stat( : data%QPP_map%n ) )
        END IF

!  Full restore

        IF ( control%restore_problem >= 2 ) THEN
          CALL QPP_restore( data%QPP_map, data%QPP_inform,                     &
                            prob, get_f = .TRUE., get_g = .TRUE.,              &
                            get_x_bounds = .TRUE., get_c_bounds = .TRUE.,      &
                            get_x = .TRUE., get_y = .TRUE., get_z = .TRUE.,    &
                            get_c = .TRUE., get_A = .TRUE., get_H = .TRUE. )

!  Restore vectors and scalars

        ELSE IF ( control%restore_problem == 1 ) THEN
          CALL QPP_restore( data%QPP_map, data%QPP_inform, prob,               &
                            get_f = .TRUE., get_g = .TRUE.,                    &
                            get_x = .TRUE., get_x_bounds = .TRUE.,             &
                            get_y = .TRUE., get_z = .TRUE.,                    &
                            get_c = .TRUE., get_c_bounds = .TRUE. )

!  Recover solution

        ELSE
          CALL QPP_restore( data%QPP_map, data%QPP_inform, prob,               &
                            get_x = .TRUE., get_y = .TRUE.,                    &
                            get_z = .TRUE., get_c = .TRUE. )
        END IF
        CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )

        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum
        prob%new_problem_structure = data%new_problem_structure
        data%save_structure = .TRUE.
      END IF

!  Compute total time

      CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
      IF ( printi ) WRITE( control%out, 2000 ) inform%time%total,              &
        inform%time%preprocess, inform%time%analyse, inform%time%factorize,    &
        inform%time%solve

  800 CONTINUE 
      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' leaving LSQP_solve ' )" )

      RETURN  

!  Allocation error

  900 CONTINUE 
      inform%status = GALAHAD_error_allocate
      CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
      IF ( printi ) WRITE( control%out, 2900 )                                 &
        inform%bad_alloc, inform%alloc_status
      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' leaving LSQP_solve ' )" )

      RETURN  

!  Non-executable statements

 2000 FORMAT( /, 14X, ' =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=',&
              /, 14X, ' =                  LSQP total time                  =',&
              /, 14X, ' =', 16X, 0P, F12.2, 23x, '='                           &
              /, 14X, ' =    preprocess    analyse    factorize     solve   =',&
              /, 14X, ' =', 4F12.2, 3x, '=',                                   &
              /, 14X, ' =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=' )
 2010 FORMAT( ' ', /, '   **  Error return ',I3,' from LSQP ' ) 
 2040 FORMAT( '   **  Error return ', I6, ' from ', A15 ) 
 2900 FORMAT( ' ** Message from -LSQP_solve-', /,                              &
              ' Allocation error, for ', A, /, ' status = ', I6 ) 

!  End of LSQP_solve

      END SUBROUTINE LSQP_solve

!-*-*-*-*-*-   L S Q P _ S O L V E _ M A I N   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE LSQP_solve_main( dims, n, m, A_val, A_col, A_ptr,             &
                                  C_l, C_u, X_l, X_u, C_RES, X, Y, Z,          &
                                  RES_x, SOL_y, RES_y, BEST_y,                 &
                                  SOL, RES, BEST, HX, GRAD_L,                  &
                                  DIST_X_l, DIST_X_u, Z_l, Z_u, BARRIER_X,     &
                                  Y_l, DY_l, DIST_C_l, Y_u, DY_u, DIST_C_u,    &
                                  C, BARRIER_C, SCALE_C, DELTA, RHS, DZ_l,     &
                                  DZ_u, Abycol_val, DIAG_X, DIAG_C, IW,        &
                                  K_colptr, Abycol_ptr, Abycol_row, K,         &
                                  FACTORS, f, COEF0, COEF1, COEF2, COEF3,      &
                                  COEF4, DELTA_cor,                            &
                                  DY_cor_l, dy_l_lower, dy_l_upper,            &
                                  DY_cor_u, dy_u_lower, dy_u_upper,            &
                                  DZ_cor_l, dz_l_lower, dz_l_upper,            &
                                  DZ_cor_u, dz_u_lower, dz_u_upper,            &
                                  BREAKP, IBREAK, CNTL, control, inform,       &
                                  Hessian_kind, gradient_kind, WEIGHT, X0, G,  &
                                  C_last, X_last, Y_last, Z_last,              &
                                  C_stat, B_Stat )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Minimizes the linear/separable quadratic objective function
!
!    1/2 || W * ( x - x^0 ) ||_2^2 + g^T x + f 
!
!  subject to the constraints
!
!               (c_l)_i <= (Ax)_i <= (c_u)_i , i = 1, .... , m,
!
!    and        (x_l)_i <=   x_i  <= (x_u)_i , i = 1, .... , n,
!
!  where x is a vector of n components ( x_1, .... , x_n ),
!  A is an m by n matrix, and any of the bounds (c_l)_i, (c_u)_i
!  (x_l)_i, (x_u)_i may be infinite, using a primal-dual method.
!  The subroutine is particularly appropriate when A is sparse.
!
!  In order that many of the internal computations may be performed
!  efficiently, it is required that   
!  
!  * the variables are ordered so that their bounds appear in the order
!
!    free                      x
!    non-negativity      0  <= x
!    lower              x_l <= x
!    range              x_l <= x <= x_u   (x_l < x_u)
!    upper                     x <= x_u
!    non-positivity            x <=  0
!
!    Fixed variables are not permitted (ie, x_l < x_u for range variables). 
!
!  * the constraints are ordered so that their bounds appear in the order
!
!    equality           c_l  = A x
!    lower              c_l <= A x
!    range              c_l <= A x <= c_u
!    upper                     A x <= c_u
!
!    Free constraints are not permitted (ie, at least one of c_l and c_u
!    must be finite). Bounds with the value zero are not treated separately.
!
!  These transformations may be effected, in place, using the module
!  GALAHAD_QPP. The same module may subsequently used to recover the solution.
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!
!  dims is a structure of type LSQP_data_type, whose components hold SCALAR
!   information about the problem on input. The components will be unaltered
!   on exit. The following components must be set:
!
!   %x_free is an INTEGER variable, which must be set by the user to the
!    number of free variables. RESTRICTION: %x_free >= 0
!                 
!   %x_l_start is an INTEGER variable, which must be set by the user to the
!    index of the first variable with a nonzero lower (or lower range) bound.
!    RESTRICTION: %x_l_start >= %x_free + 1
!                 
!   %x_l_end is an INTEGER variable, which must be set by the user to the
!    index of the last variable with a nonzero lower (or lower range) bound.
!    RESTRICTION: %x_l_end >= %x_l_start
!                 
!   %x_u_start is an INTEGER variable, which must be set by the user to the
!    index of the first variable with a nonzero upper (or upper range) bound. 
!    RESTRICTION: %x_u_start >= %x_l_start
!                 
!   %x_u_end is an INTEGER variable, which must be set by the user to the
!    index of the last variable with a nonzero upper (or upper range) bound. 
!    RESTRICTION: %x_u_end >= %x_u_start
!                 
!   %c_equality is an INTEGER variable, which must be set by the user to the
!    number of equality constraints, m. RESTRICTION: %c_equality >= 0
!                 
!   %c_l_start is an INTEGER variable, which must be set by the user to the
!    index of the first inequality constraint with a lower (or lower range) 
!    bound. RESTRICTION: %c_l_start = %c_equality + 1 
!    (strictly, this information is redundant!)
!                 
!   %c_l_end is an INTEGER variable, which must be set by the user to the
!    index of the last inequality constraint with a lower (or lower range) 
!    bound. RESTRICTION: %c_l_end >= %c_l_start
!                 
!   %c_u_start is an INTEGER variable, which must be set by the user to the
!    index of the first inequality constraint with an upper (or upper range) 
!    bound. RESTRICTION: %c_u_start >= %c_l_start
!    (strictly, this information is redundant!)
!                 
!   %c_u_end is an INTEGER variable, which must be set by the user to the
!    index of the last inequality constraint with an upper (or upper range) 
!    bound. RESTRICTION: %c_u_end = %m
!    (strictly, this information is redundant!)
!
!   %nc is an INTEGER variable, which must be set by the user to the
!    value dims%c_u_end - dims%c_l_start + 1
!
!   %x_s is an INTEGER variable, which must be set by the user to the
!    value 1
!
!   %x_e is an INTEGER variable, which must be set by the user to the
!    value n
!
!   %c_s is an INTEGER variable, which must be set by the user to the
!    value dims%x_e + 1 
!
!   %c_e is an INTEGER variable, which must be set by the user to the
!    value dims%x_e + dims%nc
!
!   %c_b is an INTEGER variable, which must be set by the user to the
!    value dims%c_e - m
!
!   %y_s is an INTEGER variable, which must be set by the user to the
!    value dims%c_e + 1
!
!   %y_i is an INTEGER variable, which must be set by the user to the
!    value dims%c_s + m
!
!   %y_e is an INTEGER variable, which must be set by the user to the
!    value dims%c_e + m
!
!   %v_e is an INTEGER variable, which must be set by the user to the
!    value dims%y_e
!
!  A_* is used to hold the matrix A by rows. In particular:
!      A_col( : )   the column indices of the components of A
!      A_ptr( : )   pointers to the start of each row, and past the end of
!                   the last row. 
!      A_val( : )   the values of the components of A
!
!  C_l, C_u are REAL arrays of length m, which must be set by the user to 
!   the values of the arrays x_l and x_u of lower and upper bounds on x, ordered
!   as described above (strictly only C_l( dims%c_l_start : dims%c_l_end )
!   and C_u( dims%c_u_start : dims%c_u_end ) need be set, as the other
!   components are ignored!).
!  
!  X_l, X_u are REAL arrays of length n, which must be set by the user to 
!   the values of the arrays x_l and x_u of lower and upper bounds on x, ordered
!   as described above (strictly only X_l( dims%x_l_start : dims%x_l_end )
!   and X_u( dims%x_u_start : dims%x_u_end ) need be set, as the other
!   components are ignored!).
!  
!  C_RES is a REAL array of length m, which need not be set on entry. On exit,
!   the i-th component of C_RES will contain (A * x)_i, for i = 1, .... , m. 
!
!  X is a REAL array of length n, which must be set by
!   the user on entry to LSQP_solve to give an initial estimate of the 
!   optimization parameters, x. The i-th component of X should contain 
!   the initial estimate of x_i, for i = 1, .... , n.  The estimate need 
!   not satisfy the simple bound constraints and may be perturbed by 
!   LSQP_solve prior to the start of the minimization.  Any estimate which is 
!   closer to one of its bounds than control%prfeas may be reset to try to
!   ensure that it is at least control%prfeas from its bounds. On exit from 
!   LSQP_solve, X will contain the best estimate of the optimization 
!   parameters found
!  
!  Y is a REAL array of length m, which must be set by the user
!   on entry to LSQP_solve to give an initial estimates of the
!   optimal Lagrange multipiers, y. The i-th component of Y 
!   should contain the initial estimate of y_i, for i = 1, .... , m.  
!   Any estimate which is smaller than control%dufeas may be 
!   reset to control%dufeas. The dual variable for any variable with both
!   On exit from LSQP_solve, Y will contain the best estimate of
!   the Lagrange multipliers found
!  
!  Z, is a REAL array of length n, which must be set by
!   on entry to LSQP_solve to hold the values of the the dual variables 
!   associated with the simple bound constraints. 
!   Any estimate which is smaller than control%dufeas may be 
!   reset to control%dufeas. The dual variable for any variable with both
!   infinite lower and upper bounds need not be set. On exit from
!   LSQP_solve, Z will contain the best estimates obtained
!  
!  control and inform are exactly as for LSQP_solve
!
!  Hessian_kind is an INTEGER variable which defines the type of objective
!    function to be used. Possible values are
!
!     0  all the weights will be zero, and the analytic centre of the 
!        feasible region will be found. WEIGHT (see below) need not be set
!
!     1  all the weights will be one. WEIGHT (see below) need not be set
!
!     any other value - the weights will be those given by WEIGHT (see below)
!
!   WEIGHT is an optional REAL array, which need only be included if 
!    Hessian_kind is not 0 or 1. If this is so, it must be of length at least 
!    n, and contain the weights W for the objective function. 
!  
!   X0 is an optional REAL array, which need only be included if 
!    Hessian_kind is not 0. If this is so, it must be of length at least 
!    n, and contain the shifts X^0 for the objective function. 
!  
!   gradient_kind is an INTEGER variable which defines the type of linear
!    term of the objective function to be used. Possible values are
!
!     0  the linear term will be zero, and the analytic centre of the 
!        feasible region will be found if in addition Hessian_kind is 0.
!        G (see below) need not be set
!
!     1  each component of the linear terms g will be one. 
!        G (see below) need not be set
!
!     any other value - the gradients will be those given by G (see below)
!
!   G is an optional REAL array, which need only be included if 
!    gradient_kind is not 0 or 1. If this is so, it must be of length at least 
!    n, and contain the gradient term g for the objective function. 
!  
!  The remaining arguments are used as internal workspace, and need not be 
!  set on entry
!  
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, Hessian_kind, gradient_kind
      INTEGER, INTENT( IN ) :: dy_l_lower, dy_l_upper, dy_u_lower, dy_u_upper
      INTEGER, INTENT( IN ) :: dz_l_lower, dz_l_upper, dz_u_lower, dz_u_upper
      REAL ( KIND = wp ), INTENT( IN ) :: f
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      INTEGER, INTENT( IN ), DIMENSION( A_ptr( m + 1 ) - 1 ) :: A_col
      INTEGER, INTENT( OUT ), OPTIONAL, DIMENSION( m ) :: C_stat
      INTEGER, INTENT( OUT ), OPTIONAL, DIMENSION( n ) :: B_stat
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: C_l, C_u
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X_l, X_u
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: Z
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ), OPTIONAL :: WEIGHT, X0
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ), OPTIONAL :: G
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( A_ptr( m + 1 ) - 1 ) :: A_val
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: RES_x
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
                          DIMENSION( m ) :: C_RES, SOL_y, BEST_y, RES_y
      REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL, DIMENSION( n ) :: X_last
      REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL, DIMENSION( m ) :: C_last
      REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL, DIMENSION( m ) :: Y_last
      REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL, DIMENSION( n ) :: Z_last
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
             DIMENSION( dims%v_e ) :: SOL, RES, BEST
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%v_e ) :: DELTA, RHS
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( dims%v_e ) :: HX
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( dims%c_e ) :: GRAD_L
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
                          DIMENSION( dims%x_l_start : dims%x_l_end ) :: DIST_X_l
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
                          DIMENSION( dims%x_u_start : dims%x_u_end ) :: DIST_X_u
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%x_free + 1 : dims%x_l_end ) ::  Z_l, DZ_l
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%x_u_start : n ) :: Z_u, DZ_u
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
                          DIMENSION( dims%x_free + 1 : n ) :: BARRIER_X
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%c_l_start : dims%c_l_end ) :: Y_l, DY_l, DIST_C_l
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%c_u_start : dims%c_u_end ) :: Y_u, DY_u, DIST_C_u
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%c_l_start : dims%c_u_end ) :: C, BARRIER_C, SCALE_C

      INTEGER, ALLOCATABLE, DIMENSION( : ) :: IW, K_colptr, Abycol_ptr,        &
                                              Abycol_row
      INTEGER, DIMENSION( : ) :: IBREAK
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIAG_X, DIAG_C,       &
                                                         Abycol_val
      REAL ( KIND = wp ), DIMENSION( : ) :: COEF0, COEF1, COEF2, COEF3, COEF4
      REAL ( KIND = wp ), DIMENSION( : ) :: BREAKP
      REAL ( KIND = wp ), DIMENSION( dy_l_lower : dy_l_upper ) :: DY_cor_l
      REAL ( KIND = wp ), DIMENSION( dy_u_lower : dy_u_upper ) :: DY_cor_u
      REAL ( KIND = wp ), DIMENSION( dz_l_lower : dz_l_upper ) :: DZ_cor_l
      REAL ( KIND = wp ), DIMENSION( dz_u_lower : dz_u_upper ) :: DZ_cor_u
      REAL ( KIND = wp ), DIMENSION( : ) :: DELTA_cor

      TYPE ( SMT_type ), INTENT( INOUT ) :: K
      TYPE ( SILS_factors ), INTENT( INOUT ) :: FACTORS
      TYPE ( SILS_control ), INTENT( INOUT ) :: CNTL
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      TYPE ( LSQP_inform_type ), INTENT( INOUT ) :: inform

!  Parameters

      REAL ( KIND = wp ), PARAMETER :: eta = tenm4
      REAL ( KIND = wp ), PARAMETER :: sigma_max = point01
!     REAL ( KIND = wp ), PARAMETER :: sigma_max = point1
      REAL ( KIND = wp ), PARAMETER :: gamma_b0 = tenm5
      REAL ( KIND = wp ), PARAMETER :: gamma_f0 = tenm5
!     REAL ( KIND = wp ), PARAMETER :: gamma_b0 = 0.9
!     REAL ( KIND = wp ), PARAMETER :: gamma_f0 = 0.9
      REAL ( KIND = wp ), PARAMETER :: degen_tol = tenm5

!  Local variables

      INTEGER :: A_ne, i, l, lk, liw, ierr, start_print, stop_print, print_level
      INTEGER :: ldiag_x, nbnds, nbnds_x, nbnds_c, muzero_fixed
      INTEGER :: nbact, itref_max, ldiag_c_l, ldiag_c_u
      INTEGER :: out, error, factor, nnzks, it_best, infeas_max
      REAL :: time, dum, time_start
      REAL ( KIND = wp ) :: pjgnrm, errorg, mu, pmax, amax, gamma_f, nu
      REAL ( KIND = wp ) :: alpha, alpha_b, alpha_f, sigma, slope, gamma_b
      REAL ( KIND = wp ) :: gi, merit, merit_model, merit_trial
      REAL ( KIND = wp ) :: res_prim, res_dual, slknes, slkmin, potential_trial
      REAL ( KIND = wp ) :: slknes_x, slknes_c, slkmax_x, slkmax_c, slknes_req
      REAL ( KIND = wp ) :: slkmin_x, slkmin_c, merit_best, reduce_infeas
      REAL ( KIND = wp ) :: prfeas, dufeas, p_min, p_max, d_min, d_max
      REAL ( KIND = wp ) :: step, errorc, one_minus_alpha, pivot_tol, balance
      REAL ( KIND = wp ) :: pmax_cor
      LOGICAL :: set_printt, set_printi, set_printw, set_printd, set_printe
      LOGICAL :: set_printp, printt, printi, printp, printe, printd, printw
      LOGICAL :: get_factors, refact, use_corrector, maxpiv, stat_required
      LOGICAL :: get_stat, use_scale_c = .FALSE.
      CHARACTER ( LEN = 1 ) :: re, mo, co
      TYPE ( SILS_finfo ) :: FINFO
      INTEGER :: sif = 50
!     LOGICAL :: generate_sif = .TRUE.
      LOGICAL :: generate_sif = .FALSE.

      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' entering LSQP_solve_main ' )" )

! -------------------------------------------------------------------
!  If desired, generate a SIF file for problem passed 

      IF ( generate_sif .AND. PRESENT( G ) ) THEN
        WRITE( sif, "( 'NAME          LSQPB_OUT', //, 'VARIABLES', / )" )
        DO i = 1, n
          WRITE( sif, "( '    X', I8 )" ) i
        END DO

        WRITE( sif, "( /, 'GROUPS', / )" )
        DO i = 1, n
          IF ( G( i ) /= zero )                                                &
            WRITE( sif, "( ' N  OBJ      ', ' X', I8, ' ', ES12.5 )" ) i, G( i )
        END DO
        DO i = 1, dims%c_l_start - 1
          DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
            WRITE( sif, "( ' E  C', I8, ' X', I8, ' ', ES12.5 )" )             &
              i, A_col( l ), A_val( l )
          END DO
        END DO
        DO i = dims%c_l_start, dims%c_l_end
          DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
            WRITE( sif, "( ' G  C', I8, ' X', I8, ' ', ES12.5 )" )             &
              i, A_col( l ), A_val( l )
          END DO
        END DO
        DO i = dims%c_l_end + 1, dims%c_u_end
          DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
            WRITE( sif, "( ' L  C', I8, ' X', I8, ' ', ES12.5 )" )             &
              i, A_col( l ), A_val( l )
          END DO
        END DO

        WRITE( sif, "( /, 'CONSTANTS', / )" )
        DO i = 1, dims%c_l_end
          IF ( C_l( i ) /= zero )                                              &
          WRITE( sif, "( '    RHS      ', ' C', I8, ' ', ES12.5 )" ) i, C_l( i )
        END DO
        DO i = dims%c_l_end + 1, dims%c_u_end
          IF ( C_u( i ) /= zero )                                              &
          WRITE( sif, "( '    RHS      ', ' C', I8, ' ', ES12.5 )" ) i, C_u( i )
        END DO

        IF ( dims%c_u_start <= dims%c_l_end ) THEN
          WRITE( sif, "( /, 'RANGES', / )" )
          DO i = dims%c_u_start, dims%c_l_end
            WRITE( sif, "( '    RANGE    ', ' C', I8, ' ', ES12.5 )" )        &
              i, C_u( i ) - C_l( i )
          END DO
        END IF

        IF ( dims%x_free /= 0 .OR. dims%x_l_start <= n ) THEN
          WRITE( sif, "( /, 'BOUNDS', /, ' FR BND       ''DEFAULT''' )" )
          DO i = dims%x_free + 1, dims%x_l_start - 1
            WRITE( sif, "( ' LO BND       X', I8, ' ', ES12.5 )" ) i, zero
          END DO
          DO i = dims%x_l_start, dims%x_l_end
            WRITE( sif, "( ' LO BND       X', I8, ' ', ES12.5 )" ) i, X_l( i )
          END DO
          DO i = dims%x_u_start, dims%x_u_end
            WRITE( sif, "( ' UP BND       X', I8, ' ', ES12.5 )" ) i, X_u( i )
          END DO
          DO i = dims%x_u_end + 1, n
            WRITE( sif, "( ' UP BND       X', I8, ' ', ES12.5 )" ) i, zero
          END DO
        END IF

        WRITE( sif, "( /, 'START POINT', / )" )
        DO i = 1, n
          IF ( X( i ) /= zero )                                                &
            WRITE( sif, "( ' V  START    ', ' X', I8, ' ', ES12.5 )" ) i, X( i )
        END DO

        WRITE( sif, "( /, 'ENDATA' )" )
      END IF

!  SIF file generated
! -------------------------------------------------------------------

!  Initialize time

      CALL CPU_TIME( time_start )

!  ===========================
!  Control the output printing
!  ===========================

      print_level = 0
      IF ( control%start_print < 0 ) THEN
        start_print = - 1
      ELSE
        start_print = control%start_print
      END IF

      IF ( control%stop_print < 0 ) THEN
        stop_print = control%maxit
      ELSE
        stop_print = control%stop_print
      END IF

      error = control%error ; out = control%out 

      set_printe = error > 0 .AND. control%print_level >= 1

!  Basic single line of output per iteration

      set_printi = out > 0 .AND. control%print_level >= 1 

!  As per printi, but with additional timings for various operations

      set_printt = out > 0 .AND. control%print_level >= 2 

!  As per printt but also with an indication of where in the code we are

      set_printp = out > 0 .AND. control%print_level >= 3

!  As per printp but also with details of innner iterations

      set_printw = out > 0 .AND. control%print_level >= 4

!  Full debugging printing with significant arrays printed

      set_printd = out > 0 .AND. control%print_level >= 5

!  Start setting control parameters

      IF ( inform%iter >= start_print .AND. inform%iter < stop_print ) THEN
        printe = set_printe ; printi = set_printi ; printt = set_printt
        printp = set_printp ; printw = set_printw ; printd = set_printd
        print_level = control%print_level
      ELSE
        printe = .FALSE. ; printi = .FALSE. ; printt = .FALSE.
        printp = .FALSE. ; printw = .FALSE. ; printd = .FALSE.
        print_level = 0
      END IF

      factor = control%factor
      IF ( factor < 0 .OR. factor > 2 ) THEN
        IF ( printi ) WRITE( out,                                             &
          "( ' factor = ', I6, ' out of range [0,2]. Reset to 0' )" ) factor
        factor = 0
      END IF

!  If there are no variables, exit

      IF ( n == 0 ) THEN 
        i = COUNT( ABS( C_l( : dims%c_equality ) ) > control%stop_p ) +        &
            COUNT( C_l( dims%c_l_start : dims%c_l_end ) > control%stop_p ) +   &
            COUNT( C_u( dims%c_u_start : dims%c_u_end ) < - control%stop_p )
        IF ( i == 0 ) THEN
          inform%status = GALAHAD_ok
        ELSE
          inform%status = GALAHAD_error_primal_infeasible
        END IF
        IF ( control%array_syntax_worse_than_do_loop ) THEN
          DO i = 1, m ; C_RES( i ) = zero ; Y( i ) = zero ; END DO
        ELSE
          C_RES = zero ; Y = zero
        END IF
        inform%obj = zero
        GO TO 810
      END IF 

!  Record array size

      A_ne = A_ptr( m + 1 ) - 1

!  Set control parameters

      muzero_fixed = control%muzero_fixed
      itref_max = control%itref_max
      prfeas = MAX( control%prfeas, epsmch )
      dufeas = MAX( control%dufeas, epsmch )
      reduce_infeas = MAX( epsmch,                                             &
                           MIN( control%reduce_infeas ** 2, one - epsmch ) )
      infeas_max = MAX( 0, control%infeas_max )
      stat_required = PRESENT( C_stat ) .AND. PRESENT( B_stat )
      IF ( stat_required ) THEN
        B_stat  = 0
        C_stat( : dims%c_equality ) = - 1
        C_stat( dims%c_equality + 1 : ) = 0
      END IF
      get_stat = .FALSE.

      CNTL%u = control%pivot_tol
      CNTL%lp = - 1 ; CNTL%mp = - 1 ; CNTL%wp = - 1
      IF ( print_level > 3 ) THEN
        IF ( print_level < 10 ) THEN
          CNTL%ldiag = 1
        ELSE
          CNTL%ldiag = 2
        END IF
      END IF

!  If required, write out the problem

      IF ( printd ) WRITE( out, 2150 ) ' a ', ( ( i, A_col( l ), A_val( l ),   &
                           l = A_ptr( i ), A_ptr( i + 1 ) - 1 ), i = 1, m )

      IF ( control%balance_initial_complentarity ) THEN
        IF ( control%muzero <= zero ) THEN
          balance = one
        ELSE
          balance = control%muzero
        END IF
      END IF

!  Record the initial point, move the starting point away from any bounds, 
!  and move that for dual variables away from zero

      nbnds_x = 0

!  The variable is free

      IF ( printd ) THEN
        WRITE( out,                                                            &
        "( /, 5X, 'i', 6x, 'x', 10X, 'x_l', 9X, 'x_u', 9X, 'z_l', 9X, 'z_u')" )
        DO i = 1, dims%x_free
          WRITE( out, "( I6, ES12.4, 4( '      -     ') )" ) i, X( i )
        END DO
      END IF

!  The variable is a non-negativity

      DO i = dims%x_free + 1, dims%x_l_start - 1
        nbnds_x = nbnds_x + 1
        X( i ) = MAX( X( i ), prfeas )
        IF ( control%balance_initial_complentarity ) THEN
          Z_l( i ) = balance / X( i )
        ELSE
          Z_l( i ) = MAX( ABS( Z( i ) ), dufeas )
        END IF
        IF ( printd ) WRITE( out, "( I6, 2ES12.4, '      -     ', ES12.4,      &
       &  '      -     ' )" ) i, X( i ), zero, Z_l( i )
      END DO

!  The variable has just a lower bound

      DO i = dims%x_l_start, dims%x_u_start - 1
        nbnds_x = nbnds_x + 1
        X( i ) = MAX( X( i ), X_l( i ) + prfeas )
        DIST_X_l( i ) = X( i ) - X_l( i )
        IF ( control%balance_initial_complentarity ) THEN
          Z_l( i ) = balance / DIST_X_l( i )
        ELSE
          Z_l( i ) = MAX( ABS( Z( i ) ), dufeas )
        END IF
        IF ( printd ) WRITE( out, "( I6, 2ES12.4, '      -     ', ES12.4,      &
       &  '      -     ' )" ) i, X( i ), X_l( i ), Z_l( i )
      END DO

!  The variable has both lower and upper bounds

      DO i = dims%x_u_start, dims%x_l_end

!  Check that range constraints are not simply fixed variables,
!  and that the upper bounds are larger than the corresponing lower bounds

        IF ( X_u( i ) - X_l( i ) <= epsmch ) THEN 
          inform%status = GALAHAD_error_bad_bounds
          GO TO 700
        END IF
        nbnds_x = nbnds_x + 2
        IF ( X_l( i ) + prfeas >= X_u( i ) - prfeas ) THEN 
          X( i ) = half * ( X_l( i ) + X_u( i ) ) 
        ELSE 
          X( i ) = MIN( MAX( X( i ), X_l( i ) + prfeas ), X_u( i ) - prfeas ) 
        END IF 
        DIST_X_l( i ) = X( i ) - X_l( i ) ; DIST_X_u( i ) = X_u( i ) - X( i )
        IF ( control%balance_initial_complentarity ) THEN
          Z_l( i ) = balance / DIST_X_l( i )
          Z_u( i ) = - balance / DIST_X_u( i )
        ELSE
          Z_l( i ) = MAX(   ABS( Z( i ) ),   dufeas )  
          Z_u( i ) = MIN( - ABS( Z( i ) ), - dufeas )
        END IF
        IF ( printd ) WRITE( out, "( I6, 5ES12.4 )" )                          &
             i, X( i ), X_l( i ), X_u( i ), Z_l( i ), Z_u( i )
      END DO

!  The variable has just an upper bound

      DO i = dims%x_l_end + 1, dims%x_u_end
        nbnds_x = nbnds_x + 1
        X( i ) = MIN( X( i ), X_u( i ) - prfeas )
        DIST_X_u( i ) = X_u( i ) - X( i )
        IF ( control%balance_initial_complentarity ) THEN
          Z_u( i ) = - balance / DIST_X_u( i )
        ELSE
          Z_u( i ) = MIN( - ABS( Z( i ) ), - dufeas ) 
        END IF
        IF ( printd ) WRITE( out, "( I6, ES12.4, '      -     ', ES12.4,       &
       &  '      -     ', ES12.4 )" ) i, X( i ), X_u( i ), Z_u( i )
      END DO

!  The variable is a non-positivity

      DO i = dims%x_u_end + 1, n
        nbnds_x = nbnds_x + 1
        X( i ) = MIN( X( i ), - prfeas )
        IF ( control%balance_initial_complentarity ) THEN
          Z_u( i ) = balance / X( i )
        ELSE
          Z_u( i ) = MIN( - ABS( Z( i ) ), - dufeas ) 
        END IF
        IF ( printd ) WRITE( out, "( I6, ES12.4, '      -     ', ES12.4,       &
       &  '      -     ',  ES12.4 )" ) i, X( i ), zero, Z_u( i )
      END DO

!  Compute the value of the constraint, and their residuals

      nbnds_c = 0
      IF ( m > 0 ) THEN

        IF ( control%array_syntax_worse_than_do_loop ) THEN
          DO i = 1, dims%c_equality ; C_RES( i ) = - C_l( i ) ; END DO
          DO i = dims%c_l_start, dims%c_u_end ; C_RES( i ) = zero ; END DO
        ELSE
          C_RES( : dims%c_equality ) = - C_l( : dims%c_equality )
          C_RES( dims%c_l_start : dims%c_u_end ) = zero
        END IF
        CALL LSQP_AX( m, C_RES, m, A_ne, A_val, A_col, A_ptr,        &
                      n, X, '+ ' )
        IF ( printd ) THEN
          WRITE( out,                                                          &
          "( /, 5X,'i', 6x, 'c', 10X, 'c_l', 9X, 'c_u', 9X, 'y_l', 9X, 'y_u' )")
          DO i = 1, dims%c_l_start - 1
            WRITE( out, "( I6, 3ES12.4 )" ) i, C_RES( i ), C_l( i ), C_u( i )
          END DO
        END IF

!  The constraint has just a lower bound

        DO i = dims%c_l_start, dims%c_u_start - 1
          nbnds_c = nbnds_c + 1

!  Compute an appropriate scale factor

          IF ( use_scale_c ) THEN
            SCALE_C( i ) = MAX( one, ABS( C_RES( i ) ) )
          ELSE
            SCALE_C( i ) = one
          END IF

!  Scale the bounds

          C_l( i ) = C_l( i ) / SCALE_C( i )

!  Compute an appropriate initial value for the slack variable

          C( i ) = MAX( C_RES( i ) / SCALE_C( i ), C_l( i ) + prfeas ) 
          DIST_C_l( i ) = C( i ) - C_l( i )
          C_RES( i ) = C_RES( i ) - SCALE_C( i ) * C( i )
          IF ( control%balance_initial_complentarity ) THEN
            Y_l( i ) = balance / DIST_C_l( i )
          ELSE
            Y_l( i ) = MAX( ABS( SCALE_C( i ) * Y( i ) ),  dufeas )
          END IF
          IF ( printd ) WRITE( out,  "( I6, 2ES12.4, '      -     ', ES12.4,   &
         &  '      -     ' )" ) i, C_RES( i ), C_l( i ), Y_l( i )
        END DO

!  The constraint has both lower and upper bounds

        DO i = dims%c_u_start, dims%c_l_end

!  Check that range constraints are not simply fixed variables,
!  and that the upper bounds are larger than the corresponing lower bounds

          IF ( C_u( i ) - C_l( i ) <= epsmch ) THEN 
            inform%status = GALAHAD_error_bad_bounds
            GO TO 700
          END IF
          nbnds_c = nbnds_c + 2

!  Compute an appropriate scale factor

          IF ( use_scale_c ) THEN
            SCALE_C( i ) = MAX( one, ABS( C_RES( i ) ) )
          ELSE
            SCALE_C( i ) = one
          END IF

!  Scale the bounds

          C_l( i ) = C_l( i ) / SCALE_C( i ) 
          C_u( i ) = C_u( i ) / SCALE_C( i )

!  Compute an appropriate initial value for the slack variable

          IF ( C_l( i ) + prfeas >= C_u( i ) - prfeas ) THEN 
            C( i ) = half * ( C_l( i ) + C_u( i ) ) 
          ELSE 
            C( i ) = MIN( MAX( C_RES( i ) / SCALE_C( i ), C_l( i ) + prfeas ), &
                               C_u( i ) - prfeas ) 
          END IF 
          DIST_C_l( i ) = C( i ) - C_l( i ) 
          DIST_C_u( i ) = C_u( i ) - C( i )
          C_RES( i ) = C_RES( i ) - SCALE_C( i ) * C( i )
          IF ( control%balance_initial_complentarity ) THEN
            Y_l( i ) = balance / DIST_C_l( i )
            Y_u( i ) = - balance / DIST_C_u( i )
          ELSE
            Y_l( i ) = MAX(   ABS( SCALE_C( i ) * Y( i ) ),   dufeas )
            Y_u( i ) = MIN( - ABS( SCALE_C( i ) * Y( i ) ), - dufeas )
          END IF
          IF ( printd ) WRITE( out, "( I6, 5ES12.4 )" )                        &
            i, C_RES( i ), C_l( i ), C_u( i ), Y_l( i ), Y_u( i )
        END DO

!  The constraint has just an upper bound

        DO i = dims%c_l_end + 1, dims%c_u_end
          nbnds_c = nbnds_c + 1

!  Compute an appropriate scale factor

          IF ( use_scale_c ) THEN
            SCALE_C( i ) = MAX( one, ABS( C_RES( i ) ) )
          ELSE
            SCALE_C( i ) = one
          END IF

!  Scale the bounds

          C_u( i ) = C_u( i ) / SCALE_C( i )

!  Compute an appropriate initial value for the slack variable

          C( i ) = MIN( C_RES( i ) / SCALE_C( i ), C_u( i ) - prfeas ) 
          DIST_C_u( i ) = C_u( i ) - C( i )
          C_RES( i ) = C_RES( i ) - SCALE_C( i ) * C( i )
          IF ( control%balance_initial_complentarity ) THEN
            Y_u( i ) = - balance / DIST_C_u( i )
          ELSE
            Y_u( i ) = MIN( - ABS( SCALE_C( i ) * Y( i ) ), - dufeas ) 
          END IF
          IF ( printd ) WRITE( out, "( I6, ES12.4, '      -     ', ES12.4,     &
         &  '      -     ', ES12.4 )" ) i, C_RES( i ), C_u( i ), Y_u( i )
        END DO
        res_prim = MAXVAL( ABS( C_RES ) )
      ELSE
        res_prim = zero
      END IF

!  Set useful pointers to stop Compaq optimization inefficiencies

!      ptr_GRAD_L => GRAD_L( dims%x_s : dims%x_e )
!      ptr_HX => HX( : n )
!      ptr_HX_v => HX( : dims%v_e )
!      ptr_DELTA_x => DELTA( dims%x_s : dims%x_e )
!     ptr_DELTA_c => DELTA( dims%c_s : dims%c_e )
!      ptr_DELTA_y => DELTA( dims%y_s : dims%y_e )
!      ptr_RHS_x => RHS( dims%x_s : dims%x_e )
!      ptr_RHS_c => RHS( dims%c_s : dims%c_e )
!      ptr_RHS_y => RHS( dims%y_s : dims%y_e )
      IF ( control%use_corrector ) THEN
!        ptr_DELTA_cor_x => DELTA_cor( dims%x_s : dims%x_e )
!        ptr_DELTA_cor_c => DELTA_cor( dims%c_s : dims%c_e )
!        ptr_DELTA_cor_y => DELTA_cor( dims%y_s : dims%y_e )
      END IF

!  Find the max-norm of the residual

      nbnds = nbnds_x + nbnds_c
      IF ( printi )                                                            &
        WRITE( out, "( /, ' >>>>> factorization package SILS used <<<<<', / )" )
      IF ( printi .AND. m > 0 .AND. dims%c_l_start <= dims%c_u_end )           &
        WRITE( out, "( ' largest/smallest scale factor ', 2ES12.4 )" )         &
          MAXVAL( SCALE_C ), MINVAL( SCALE_C )

!  Compute the complementary slackness

!       DO i = dims%x_free + 1, dims%x_l_start - 1
!         write(6,"(I6, ' x lower', 2ES12.4)" ) i, X( i ), Z_l( i )
!       END DO
!       DO i = dims%x_l_start, dims%x_l_end
!         write(6,"(I6, ' x lower', 2ES12.4)" ) i, DIST_X_l( i ), Z_l( i )
!       END DO
!       DO i = dims%x_u_start, dims%x_u_end
!         write(6,"(I6, ' x upper', 2ES12.4)" ) i, - DIST_X_u( i ), Z_u( i )
!       END DO
!       DO i = dims%x_u_end + 1, n
!         write(6,"(I6, ' x upper', 2ES12.4)" ) i, X( i ), Z_u( i )
!       END DO

!       DO i = dims%c_l_start, dims%c_l_end
!         write(6,"(I6, ' c lower', 2ES12.4)" ) i, DIST_C_l( i ), Y_l( i )
!       END DO
!       DO i = dims%c_u_start, dims%c_u_end
!         write(6,"(I6, ' c upper', 2ES12.4)" ) i, - DIST_C_u( i ), Y_u( i )
!       END DO

      slknes_x = DOT_PRODUCT( X( dims%x_free + 1 : dims%x_l_start - 1 ),       &
                              Z_l( dims%x_free + 1 : dims%x_l_start - 1 ) ) +  &
                 DOT_PRODUCT( DIST_X_l( dims%x_l_start : dims%x_l_end ),       &
                              Z_l( dims%x_l_start : dims%x_l_end ) ) -         &
                 DOT_PRODUCT( DIST_X_u( dims%x_u_start : dims%x_u_end ),       &
                              Z_u( dims%x_u_start : dims%x_u_end ) ) +         &
                 DOT_PRODUCT( X( dims%x_u_end + 1 : n ),                       &
                              Z_u( dims%x_u_end + 1 : n ) )
      slknes_c = DOT_PRODUCT( DIST_C_l( dims%c_l_start : dims%c_l_end ),       &
                              Y_l( dims%c_l_start : dims%c_l_end ) ) -         &
                 DOT_PRODUCT( DIST_C_u( dims%c_u_start : dims%c_u_end ),       &
                              Y_u( dims%c_u_start : dims%c_u_end ) )
      slknes = slknes_x + slknes_c

      slkmin_x = MIN( MINVAL( X( dims%x_free + 1 : dims%x_l_start - 1 ) *      &
                              Z_l( dims%x_free + 1 : dims%x_l_start - 1 ) ),   &
                      MINVAL( DIST_X_l( dims%x_l_start : dims%x_l_end ) *      &
                              Z_l( dims%x_l_start : dims%x_l_end ) ),          &
                      MINVAL( - DIST_X_u( dims%x_u_start : dims%x_u_end ) *    &
                              Z_u( dims%x_u_start : dims%x_u_end ) ),          &
                      MINVAL( X( dims%x_u_end + 1 : n ) *                      &
                              Z_u( dims%x_u_end + 1 : n ) ) )
      slkmin_c = MIN( MINVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) *      &
                              Y_l( dims%c_l_start : dims%c_l_end ) ),          &
                      MINVAL( - DIST_C_u( dims%c_u_start : dims%c_u_end ) *    &
                              Y_u( dims%c_u_start : dims%c_u_end ) ) )
      slkmin = MIN( slkmin_x, slkmin_c )

      slkmax_x = MAX( MAXVAL( X( dims%x_free + 1 : dims%x_l_start - 1 ) *      &
                              Z_l( dims%x_free + 1 : dims%x_l_start - 1 ) ),   &
                      MAXVAL( DIST_X_l( dims%x_l_start : dims%x_l_end ) *      &
                              Z_l( dims%x_l_start : dims%x_l_end ) ),          &
                      MAXVAL( - DIST_X_u( dims%x_u_start : dims%x_u_end ) *    &
                              Z_u( dims%x_u_start : dims%x_u_end ) ),          &
                      MAXVAL( X( dims%x_u_end + 1 : n ) *                      &
                              Z_u( dims%x_u_end + 1 : n ) ) )
      slkmax_c = MAX( MAXVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) *      &
                              Y_l( dims%c_l_start : dims%c_l_end ) ),          &
                      MAXVAL( - DIST_C_u( dims%c_u_start : dims%c_u_end ) *    &
                              Y_u( dims%c_u_start : dims%c_u_end ) ) )

!     WRITE(6,2100) ' >0 ', X( dims%x_free + 1 : dims%x_l_start - 1 )
!     WRITE(6,2100) ' >l ', DIST_X_l( dims%x_l_start : dims%x_l_end )
!     WRITE(6,2100) ' <u ', DIST_X_u( dims%x_u_start : dims%x_u_end )
!     WRITE(6,2100) ' <0 ', - X( dims%x_u_end + 1 : n )
!     stop

    

      p_min = MIN( MINVAL( X( dims%x_free + 1 : dims%x_l_start - 1 ) ),        &
                   MINVAL( DIST_X_l( dims%x_l_start : dims%x_l_end ) ),        &
                   MINVAL( DIST_X_u( dims%x_u_start : dims%x_u_end ) ),        &
                   MINVAL( - X( dims%x_u_end + 1 : n ) ),                      &
                   MINVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) ),        &
                   MINVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) ) )

      p_max = MAX( MAXVAL( X( dims%x_free + 1 : dims%x_l_start - 1 ) ),        &
                   MAXVAL( DIST_X_l( dims%x_l_start : dims%x_l_end ) ),        &
                   MAXVAL( DIST_X_u( dims%x_u_start : dims%x_u_end ) ),        &
                   MAXVAL( - X( dims%x_u_end + 1 : n ) ),                      &
                   MAXVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) ),        &
                   MAXVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) ) )

      d_min = MIN( MINVAL(   Z_l( dims%x_free + 1 : dims%x_l_end ) ),          &
                   MINVAL( - Z_u( dims%x_u_start : n ) ),                      &
                   MINVAL(   Y_l( dims%c_l_start : dims%c_l_end ) ),           &
                   MINVAL( - Y_u( dims%c_u_start : dims%c_u_end ) ) )

      d_max = MAX( MAXVAL(   Z_l( dims%x_free + 1 : dims%x_l_end ) ),          &
                   MAXVAL( - Z_u( dims%x_u_start : n ) ),                      &
                   MAXVAL(   Y_l( dims%c_l_start : dims%c_l_end ) ),           &
                   MAXVAL( - Y_u( dims%c_u_start : dims%c_u_end ) ) )

!  Record the slackness and the deviation from the central path

      IF ( nbnds_x > 0 ) THEN
        slknes_x = slknes_x / nbnds_x
      ELSE
        slknes_x = zero
      END IF

      IF ( nbnds_c > 0 ) THEN
        slknes_c = slknes_c / nbnds_c
      ELSE
        slknes_c = zero
      END IF

      IF ( nbnds > 0 ) THEN
        gamma_f = gamma_f0 * slknes ; slknes = slknes / nbnds
        gamma_b = gamma_b0 * slkmin / slknes
      ELSE
        gamma_f = zero ; slknes = zero ; gamma_b = zero
      END IF

      IF ( printt .AND. nbnds > 0 ) WRITE( out, 2130 )                         &
        slknes, slknes_x, slknes_c, slkmin_x, slkmax_x, slkmin_c, slkmax_c,    &
        p_min, p_max, d_min, d_max

!  Compute the initial objective value

      IF ( Hessian_kind == 0 ) THEN
        inform%obj = f
      ELSE IF ( Hessian_kind == 1 ) THEN
        inform%obj = f + half * SUM( ( X - X0 ) ** 2 )
      ELSE
        inform%obj = f + half * SUM( ( WEIGHT * ( X - X0 ) ) ** 2 )
      END IF
      
      IF ( gradient_kind == 1 ) THEN
        inform%obj = inform%obj + SUM( X )
      ELSE IF ( gradient_kind /= 0 ) THEN
        inform%obj = inform%obj + DOT_PRODUCT( G, X )
      END IF

!  Find the largest components of A

      IF ( A_ne > 0 ) THEN
        amax = MAXVAL( ABS( A_val( : A_ne ) ) )
      ELSE
        amax = zero
      END IF

      IF ( printi ) WRITE( out, "( '  maximum element of A = ', ES12.4 )" ) amax

!  Test to see if we are feasible

      inform%feasible = res_prim <= control%stop_p
      pjgnrm = infinity

      IF ( inform%feasible ) THEN
        IF ( printi ) WRITE( out, 2070 )
        IF ( control%just_feasible ) THEN
          inform%status = GALAHAD_ok
          GO TO 500
        END IF
        IF ( Hessian_kind == 0 .AND. gradient_kind == 0 )                      &
          inform%potential = LSQP_potential_value( dims, n,                    &
                                   X, DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u )
      END IF

!  Set the initial barrier parameter

      sigma = sigma_max ; nu = one
      IF ( control%muzero < zero ) THEN
        mu = sigma * slknes
      ELSE
        mu = control%muzero
      END IF
      slknes_req = slknes

!  Compute the gradient of the Lagrangian function.

      CALL LSQP_Lagrangian_gradient( dims, n, m, X, Y, Y_l, Y_u, Z_l, Z_u,     &
                                     A_ne, A_val, A_col, A_ptr,                &
                                     DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u,   &
                                     GRAD_L( dims%x_s : dims%x_e ),            &
                                     control%getdua, dufeas,                   &
                                     Hessian_kind, gradient_kind, control,     &
                                     WEIGHT = WEIGHT, X0 = X0, G = G )

!  Evaluate the merit function

      merit = LSQP_merit_value( dims, n, m, X, Y, Y_l, Y_u, Z_l, Z_u,          &
                                 DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u,       &
                                 GRAD_L( dims%x_s : dims%x_e ), C_RES, res_dual )

!  Prepare for the major iteration

      inform%iter = 0 ; inform%nfacts = 0
      IF ( printt ) WRITE( out, "( ' ',/,' merit function value = ',           &
     &     ES12.4 )" ) merit 

      IF ( n == 0 ) THEN 
        inform%status = GALAHAD_ok ; GO TO 600 
      END IF 
      merit_best = merit ; it_best = 0

!  Test for convergence

      IF ( res_prim <= control%stop_p .AND. res_dual <= control%stop_d .AND.   &
           slknes_req <= control%stop_c ) THEN
        inform%status = GALAHAD_ok ; GO TO 600 
      END IF 

!  ===================================================
!  Analyse the sparsity pattern of the required matrix
!  ===================================================

      CALL CPU_TIME( time ) 
      CALL LSQP_analyse( dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C,       &
                         factor, nnzks, lk, liw, ldiag_x, ldiag_c_l,           &
                         ldiag_c_u, IW, Abycol_val, Abycol_row, Abycol_ptr,    &
                         K_colptr, DIAG_X, DIAG_C, K, FACTORS, CNTL,           &
                         print_level, control, inform )
      CALL CPU_TIME( dum ) ; dum = dum - time
      IF ( printt ) WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
      inform%time%analyse = inform%time%analyse + dum

      refact = .TRUE. ; re = ' ' ; mo = ' ' ;  nbact = 0
      pivot_tol = control%pivot_tol
      maxpiv = pivot_tol >= half

      IF ( printi ) WRITE( out,                                                &
           "(  /,'  Primal    convergence tolerence = ', ES12.4,               &
          &    /,'  Dual      convergence tolerence = ', ES12.4,               &
          &    /,'  Slackness convergence tolerence = ', ES12.4 )" )           &
           control%stop_p, control%stop_d, control%stop_c

!  ---------------------------------------------------------------------
!  ---------------------- Start of Major Iteration ---------------------
!  ---------------------------------------------------------------------

      DO

!  =====================================================================
!  -*-*-*-*-*-*-*-*-*-*-*-   Test for Optimality   -*-*-*-*-*-*-*-*-*-*-
!  =====================================================================

!  Print a summary of the iteration

        CALL CPU_TIME( time ) ; inform%time%total = time - time_start
  
        IF ( printi ) THEN 

          IF ( inform%iter > 0 ) THEN 
            IF ( printt .OR.                                                   &
              ( printi .AND. inform%iter == start_print ) ) WRITE( out, 2000 ) 
            WRITE( out, 2030 ) inform%iter, re, res_prim, res_dual,            &
              slknes_req, inform%obj, alpha, co, mo, mu, nbact, inform%time%total
          ELSE 
            WRITE( out, 2000 ) 
            WRITE( out, 2020 ) inform%iter, re, res_prim, res_dual,            &
              slknes_req, inform%obj, mu, inform%time%total
          END IF 

          IF ( printd ) THEN 
            WRITE( out, 2100 ) ' X ', X
            WRITE( out, 2100 ) ' Z_l ', Z_l
            WRITE( out, 2100 ) ' Z_u ', Z_u
          END IF 
        END IF 

!  Test for optimality

        IF ( res_prim <= control%stop_p .AND. res_dual <= control%stop_d .AND. &
             slknes_req <= control%stop_c ) THEN
          inform%status = GALAHAD_ok ; GO TO 600 
        END IF 

!  Test to see if more than maxit iterations have been performed

        inform%iter = inform%iter + 1 
        IF ( inform%iter > control%maxit ) THEN 
          inform%status = GALAHAD_error_max_iterations ; GO TO 600 
        END IF 

!  Check that the CPU time limit has not been reached

        IF ( control%cpu_time_limit >= zero .AND.                              &
             inform%time%total > control%cpu_time_limit ) THEN
          inform%status = GALAHAD_error_cpu_limit ; GO TO 600
        END IF 

        IF ( inform%iter == start_print ) THEN
          printe = set_printe ; printi = set_printi ; printt = set_printt
          printw = set_printw ; printd = set_printd
          print_level = control%print_level
        END IF

        IF ( inform%iter == stop_print + 1 ) THEN
          printe = .FALSE. ; printi = .FALSE. ; printt = .FALSE.
          printw = .FALSE. ; printd = .FALSE.
          print_level = 0
        END IF

!       WRITE( 6, "( ' start, stop print, iter ', 3I8 )" )                     &
!         start_print, stop_print, inform%iter

!  Test to see whether the method has stalled

        IF ( merit <= reduce_infeas * merit_best ) THEN
          merit_best = merit
          it_best = 0
        ELSE
          it_best = it_best + 1
          IF ( it_best > infeas_max ) THEN
            IF ( inform%feasible ) THEN
              inform%status = GALAHAD_error_no_center ; GO TO 600
            ELSE 
              IF ( printi ) WRITE( out, "( /, ' =============== the problem',  &
             &  ' appears to be infeasible ====================== ', / )" )
              inform%status = GALAHAD_error_primal_infeasible ; GO TO 600 
            END IF
          END IF
        END IF  

!  Test to see if the potential function appears to be unbounded from below

!       IF ( nbnds == 0 .AND. inform%feasible .AND. Hessian_kind == 0 .AND.    &
        IF ( inform%feasible .AND. Hessian_kind == 0 .AND.                     &
             gradient_kind == 0 ) THEN
          IF ( inform%potential < control%potential_unbounded *                &
               ( ( dims%x_l_end - dims%x_free - 1 ) +                        &
               ( n -  dims%x_u_start ) +                                       &
               ( dims%c_l_end - dims%c_l_start ) +                             &
               ( dims%c_u_end - dims%c_u_start ) ) ) THEN 
            inform%status = GALAHAD_error_no_center ; GO TO 600 
          END IF 

!  Compute the Hessian matrix of the barrier terms

!  Special case for the analytic center

!  problem variables:

          DO i = dims%x_free + 1, dims%x_l_start - 1
            BARRIER_X( i ) =  mu / X( i ) ** 2
          END DO
          DO i = dims%x_l_start, dims%x_u_start - 1
            BARRIER_X( i ) = mu / DIST_X_l( i ) ** 2
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            BARRIER_X( i ) = mu / DIST_X_l( i ) ** 2                           &
                             + mu / DIST_X_u( i ) ** 2
          END DO
          DO i = dims%x_l_end + 1, dims%x_u_end
            BARRIER_X( i ) = mu / DIST_X_u( i ) ** 2
          END DO
          DO i = dims%x_u_end + 1, n
            BARRIER_X( i ) = mu / X( i ) ** 2
          END DO

!  slack variables:

          BARRIER_C( dims%c_l_start : dims%c_u_end ) = zero
          DO i = dims%c_l_start, dims%c_u_start - 1
            BARRIER_C( i ) = mu / DIST_C_l( i ) ** 2
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            BARRIER_C( i ) = mu / DIST_C_l( i ) ** 2                           &
                             + mu / DIST_C_u( i ) ** 2
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            BARRIER_C( i ) = mu / DIST_C_u( i ) ** 2
          END DO

!  General case

        ELSE

!  problem variables:

          DO i = dims%x_free + 1, dims%x_l_start - 1
            IF ( ABS( X( i ) ) <= degen_tol .AND. printw )                     &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, X( i ), Z_l( i )
            BARRIER_X( i ) = Z_l( i ) / X( i )
          END DO
          DO i = dims%x_l_start, dims%x_u_start - 1
            IF ( ABS( DIST_X_l( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, DIST_X_l( i ), Z_l( i )
            BARRIER_X( i ) = Z_l( i ) / DIST_X_l( i ) 
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            IF ( ABS( DIST_X_l( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, DIST_X_l( i ), Z_l( i )
            IF ( ABS( DIST_X_u( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, DIST_X_u( i ), Z_u( i )
            BARRIER_X( i ) = Z_l( i ) / DIST_X_l( i ) - Z_u( i ) / DIST_X_u( i )
          END DO
          DO i = dims%x_l_end + 1, dims%x_u_end
            IF ( ABS( DIST_X_u( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, DIST_X_u( i ), Z_u( i )
            BARRIER_X( i ) = - Z_u( i ) / DIST_X_u( i )
          END DO
          DO i = dims%x_u_end + 1, n
            IF ( ABS( X( i ) ) <= degen_tol .AND. printw )                     &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, X( i ), Z_u( i )
            BARRIER_X( i ) = Z_u( i ) / X( i )
          END DO

!  slack variables:

          BARRIER_C( dims%c_l_start : dims%c_u_end ) = zero
          DO i = dims%c_l_start, dims%c_u_start - 1
            IF ( ABS( DIST_C_l( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST C, Y ', 2ES12.4 )" )            &
                i, DIST_C_l( i ), Y_l( i )
            BARRIER_C( i ) = Y_l( i ) / DIST_C_l( i ) 
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            IF ( ABS( DIST_C_l( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST C, Y ', 2ES12.4 )" )            &
                i, DIST_C_l( i ), Y_l( i )
            IF ( ABS( DIST_C_u( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST C, Y ', 2ES12.4 )" )            &
                i, DIST_C_u( i ), Y_u( i )
            BARRIER_C( i ) = Y_l( i ) / DIST_C_l( i ) - Y_u( i ) / DIST_C_u( i )
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            IF ( ABS( DIST_C_u( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST C, Y ', 2ES12.4 )" )            &
                i, DIST_C_u( i ), Y_u( i )
            BARRIER_C( i ) = - Y_u( i ) / DIST_C_u( i )
          END DO
        END IF

!  =====================================================================
!  -*-*-*-*-*-*-*-*-*-*-*-*-      Factorization      -*-*-*-*-*-*-*-*-*-
!  =====================================================================

        mo = ' '
        IF ( refact ) THEN

!  Only refactorize if B has changed

          re = 'r' 
          get_factors = .TRUE.
      
          CNTL%liw = MAX( 2 * inform%factorization_integer, control%indmin )
          CNTL%la = MAX( 2 * inform%factorization_real, control%valmin )

          CALL CPU_TIME( time ) 

!  For the Schur complement matrix
!  ===============================

          IF ( factor == 0 .OR. factor == 1 ) THEN
      
            IF ( control%array_syntax_worse_than_do_loop ) THEN
              IF ( Hessian_kind == 0 ) THEN
                DO i = 1, dims%x_free ; DIAG_X( i ) = zero ; END DO
                DO i = dims%x_free + 1, n 
                  DIAG_X( i ) = BARRIER_X( i ) ; END DO
              ELSE IF ( Hessian_kind == 1 ) THEN
                DO i = 1, dims%x_free ; DIAG_X( i ) = one ; END DO
                DO i = dims%x_free + 1, n 
                  DIAG_X( i ) = one + BARRIER_X( i ) ; END DO
              ELSE
                DO i = 1, dims%x_free
                  DIAG_X( i ) = WEIGHT( i ) ** 2 ; END DO
                DO i = dims%x_free + 1, n 
                  DIAG_X( i ) = WEIGHT( i ) ** 2 + BARRIER_X( i ) ; END DO
              END IF
              DO i = dims%c_l_start, dims%c_u_end
                DIAG_C( i ) = BARRIER_C( i ) ; END DO
            ELSE
              IF ( Hessian_kind == 0 ) THEN
                DIAG_X( : dims%x_free ) = zero
                DIAG_X( dims%x_free + 1 : ) = BARRIER_X
              ELSE IF ( Hessian_kind == 1 ) THEN
                DIAG_X( : dims%x_free ) = one
                DIAG_X( dims%x_free + 1 : ) = one + BARRIER_X
              ELSE
                DIAG_X( : dims%x_free ) = WEIGHT( : dims%x_free ) ** 2
                DIAG_X( dims%x_free + 1 : ) =                                  &
                  WEIGHT( dims%x_free + 1 : ) ** 2 + BARRIER_X
              END IF
              DIAG_C = BARRIER_C
            END IF
      
            IF ( printd ) THEN 
              WRITE( out, "( ' DIAG_X = ', /, ( 6ES10.2 ) )" ) DIAG_X
              WRITE( out, "( ' DIAG_C = ', /, ( 6ES10.2 ) )" ) DIAG_C
            END IF

!  Form the Schur complement matrix

            IF ( printw ) WRITE( out,                                          &
               "( ' ...... factorization of Schur complement  .......... ' )" )

            IF ( m > 0 ) THEN

!              ptr_K_val => K%val
!              ptr_K_row => K%row
!              ptr_IW_n => IW( : n )
!              ptr_IW_m => IW( n + 1 : )

              CALL LSQP_form_Schur_complement(                                 &
                        dims, n, m, A_ne, Abycol_val, Abycol_row, Abycol_ptr,  &
                        DIAG_X, SCALE_C, DIAG_C, lk, K%val, K%row,     &
                        K_colptr, K%ne, ierr, IW( : n ), IW( n + 1 : ), .TRUE.,      &
                        control%error, print_level )
!             inform%nfacts = inform%nfacts + 1 
            ELSE
              re = 'r' 
              FINFO%flag = 0 
              get_factors = .FALSE.
            END IF
            IF ( printw ) WRITE( out,                                          &
              "( ' ............... end of factorization ............... ' )" )

!  For the KKT matrix
!  ==================

          ELSE

!  Include the values of the barrier terms
 
            IF ( control%array_syntax_worse_than_do_loop ) THEN
              IF ( Hessian_kind == 0 ) THEN
                DO i = 1, dims%x_free ; K%val( nnzks + i ) = zero ; END DO
                DO i = dims%x_free + 1, n 
                  K%val( nnzks + i ) = BARRIER_X( i ) ; END DO
              ELSE IF ( Hessian_kind == 1 ) THEN
                DO i = 1, dims%x_free ; K%val( nnzks + i ) = one ; END DO
                DO i = dims%x_free + 1, n 
                  K%val( nnzks + i ) = one + BARRIER_X( i ) ; END DO
              ELSE
                DO i = 1, dims%x_free
                  K%val( nnzks + i ) = WEIGHT( i ) ** 2 ; END DO
                DO i = dims%x_free + 1, n 
                  K%val( nnzks + i ) = WEIGHT( i ) ** 2 + BARRIER_X( i ) 
                END DO
              END IF
              DO i = 0, dims%nc - 1
                 K%val( nnzks + dims%c_s + i ) = BARRIER_C( dims%c_l_start + i )
              END DO
            ELSE
              IF ( Hessian_kind == 0 ) THEN
                K%val( nnzks + 1 : nnzks + dims%x_free ) = zero
                K%val( nnzks + dims%x_free + 1 : nnzks + n ) = BARRIER_X
              ELSE IF ( Hessian_kind == 1 ) THEN
                K%val( nnzks + 1 : nnzks + dims%x_free ) = one
                K%val( nnzks + dims%x_free + 1 : nnzks + n ) =                 &
                    one + BARRIER_X
              ELSE
                K%val( nnzks + 1 : nnzks + dims%x_free ) =                     &
                  WEIGHT( : dims%x_free ) ** 2
                K%val( nnzks + dims%x_free + 1 : nnzks + n ) =                 &
                  WEIGHT( dims%x_free + 1 : ) ** 2 + BARRIER_X
              END IF
              K%val( nnzks + dims%c_s : nnzks + dims%c_e ) = BARRIER_C
            END IF
          END IF

! ::::::::::::::::::::::::::::::
!  Factorize the required matrix
! ::::::::::::::::::::::::::::::

          IF ( get_factors ) THEN

!           WRITE( 22, "( 2I6 )" ) K%n, K%ne
!           WRITE( 22, "( ( 10I6 ) )" ) K%row( : K%ne )
!           WRITE( 22, "( ( 10I6 ) )" ) K%col( : K%ne ) 
!           WRITE( 22, "( ( 3ES24.16 ) )" ) K%val( : K%ne ) 
!           WRITE( 22, "( ( 3ES24.16 ) )" ) SOL( : K%n )
!           STOP    

            IF ( printw ) WRITE( out,                                          &
               "( ' ......... factorization of KKT matrix ............... ' )" )
            CALL SILS_factorize( K, FACTORS, CNTL, FINFO )
            inform%nfacts = inform%nfacts + 1 
            IF ( printw ) WRITE( out,                                          &
               "( ' ............... end of factorization ............... ' )" )

!  Test that the factorization succeeded

            inform%factorization_status = FINFO%flag
            IF ( inform%factorization_status < 0 ) THEN
              IF ( printe ) WRITE( error, 2040 ) FINFO%flag, 'SILS_factorize'

!  It didn't. We might have run out of options

              IF ( factor == 2 .AND. maxpiv ) THEN
                inform%status = GALAHAD_error_factorization ; GO TO 700
              
!  ... or we may change the method

              ELSE IF ( factor < 2 .AND. maxpiv ) THEN
                factor = 2
                pivot_tol = control%pivot_tol
                maxpiv = pivot_tol >= half
                IF ( printi )                                                  &
                  WRITE( out, "( ' Switching to augmented system method' )" )

!  Re-analyse the sparsity pattern of the required matrix

                CALL CPU_TIME( time ) 
                CALL LSQP_analyse(                                             &
                  dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C, factor,      &
                  nnzks,lk, liw, ldiag_x, ldiag_c_l, ldiag_c_u, IW,            &
                  Abycol_val, Abycol_row, Abycol_ptr, K_colptr, DIAG_X,        &
                  DIAG_C, K, FACTORS, CNTL, print_level, control, inform )
                CALL CPU_TIME( dum ) ; dum = dum - time
                IF ( printt )                                                  &
                  WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
                inform%time%analyse = inform%time%analyse + dum

!  ... or we can increase the pivot tolerance

              ELSE
                maxpiv = .TRUE.
                pivot_tol = half
                IF ( printi ) WRITE( out, "( ' Pivot tolerance increased ' )" )
              END IF
              alpha = zero ; nbact = 0
              inform%factorization_integer = - 1
              inform%factorization_real = - 1
              CYCLE

!  Record warning conditions

            ELSE 
              IF ( FINFO%flag > 0 ) THEN
                IF ( printt ) THEN
                  WRITE( out, 2050 ) FINFO%flag, 'SILS_factorize'
                  IF ( FINFO%flag == 4 ) WRITE( out, "( ' ** Matrix has ', I7, &
               &     ' zero eigenvalues ' )" ) K%n - finfo%rank
                END IF 
              END IF 

!  Record the storage required

              inform%factorization_integer = FINFO%nirbdu 
              inform%factorization_real = FINFO%nrlbdu 

            END IF 
  
            IF ( printt ) WRITE( out,                                          &
              "( ' real/integer space used for factors ', 2I10 )" )            &
                FINFO%nrlbdu, FINFO%nirbdu
  
          ELSE
            inform%factorization_integer = 0 
            inform%factorization_real = 0
          END IF
          CALL CPU_TIME( dum ) ; dum = dum - time
  
          IF ( printt ) THEN
            WRITE( out, "( ' ** factorize time = ', F10.2 ) " ) dum
            WRITE( out, 2060 ) inform%factorization_integer,                   &
                               inform%factorization_real
          END IF 
          inform%time%factorize = inform%time%factorize + dum
  
        ELSE 
          re = ' ' 
        END IF 

!  =======
!  STEP 1:
!  =======

!  =======================================================================
!  -*-*-*-   Obtain the Primal-Dual (Predictor) Search Direction -*-*-*-*-
!  =======================================================================

!  :::::::::::::::::::::::::::::::::
!  Set up the right-hand-side vector
!  :::::::::::::::::::::::::::::::::

        IF ( printd ) WRITE( out, 2100 ) ' GRAD_L ', GRAD_L( dims%x_s : dims%x_e )

!  Problem variables:

        RHS( : dims%x_free ) = - GRAD_L( : dims%x_free )
        DO i = dims%x_free + 1, dims%x_l_start - 1
          RHS( i ) = - GRAD_L( i ) + mu / X( i )
        END DO
        DO i = dims%x_l_start, dims%x_u_start - 1
          RHS( i ) = - GRAD_L( i ) + mu / DIST_X_l( i )
        END DO
        DO i = dims%x_u_start, dims%x_l_end
          RHS( i ) = - GRAD_L( i ) + mu / DIST_X_l( i ) - mu / DIST_X_u( i )
        END DO
        DO i = dims%x_l_end + 1, dims%x_u_end
          RHS( i ) = - GRAD_L( i ) - mu / DIST_X_u( i )
        END DO
        DO i = dims%x_u_end + 1, n
          RHS( i ) = - GRAD_L( i ) + mu / X( i )
        END DO

!  Slack variables:

        DO i = dims%c_l_start, dims%c_u_start - 1
          RHS( dims%c_b + i ) = - Y( i ) + mu / DIST_C_l( i )
        END DO
        DO i = dims%c_u_start, dims%c_l_end
          RHS( dims%c_b + i ) = - Y( i ) + mu / DIST_C_l( i )                  &
                                       - mu / DIST_C_u( i )
        END DO
        DO i = dims%c_l_end + 1, dims%c_u_end
          RHS( dims%c_b + i ) = - Y( i ) - mu / DIST_C_u( i )
        END DO

!  Include the constraint infeasibilities

        IF ( control%array_syntax_worse_than_do_loop ) THEN
          DO i = 1, m
            RHS( dims%c_e + i ) = - C_RES( i ) ; END DO
          DO i = 1, dims%v_e ; DELTA( i ) = RHS( i ) ; END DO
        ELSE
          RHS( dims%y_s : dims%y_e ) = - C_RES
          DELTA = RHS
        END IF
  
        IF ( printd ) THEN 
          WRITE( out, 2100 ) ' RHS_x ', RHS( dims%x_s : dims%x_e )
          IF ( m > 0 )                                                         &
            WRITE( out, 2100 ) ' RHS_y ', RHS( dims%y_s : dims%y_e )
        END IF 

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Obtain the primal-dual direction for the primal variables
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!  Solve  ( H  A^T ) ( Dx^pd ) = - ( grad b )
!         ( A   0  ) ( Dy^pd )     (   r    )

        IF ( printw ) WRITE( out,                                              &
             "( ' ............... compute step  ............... ' )" )

!  Use a direct method

        CALL CPU_TIME( time )
        CALL LSQP_iterative_refinement(                                        &
                   dims, n, m, K, FACTORS, CNTL, A_ne, A_val, A_col, A_ptr,    &
                   BARRIER_X, BARRIER_C, DELTA, factor, DIAG_X, ldiag_x,       &
                   SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,                      &
                   SOL, RES, BEST, SOL_y, BEST_y, RES_y, RES_x, itref_max,     &
                   print_level, control, inform, Hessian_kind, WEIGHT = WEIGHT )
  
        IF ( inform%status /= GALAHAD_ok ) GO TO 700
        CALL CPU_TIME( dum ) ; dum = dum - time
  
        IF ( printt ) WRITE( out, "( ' ** solve time = ', F10.2 ) " ) dum
        inform%time%solve = inform%time%solve + dum

!  Compute the residual of the linear system

        CALL LSQP_residual( dims, n, m, dims%v_e, A_ne, A_val, A_col, A_ptr,   &
                            DELTA( dims%x_s : dims%x_e ),                      &
                            DELTA( dims%c_s : dims%c_e ),                      &
                            DELTA( dims%y_s : dims%y_e ),                      &
                            RHS( dims%x_s : dims%x_e ),                        &
                            RHS( dims%c_s : dims%c_e ),                        &
                            RHS( dims%y_s : dims%y_e ), HX( : dims%v_e ),      &
                            BARRIER_X, BARRIER_C, SCALE_C, errorg, errorc,     &
                            print_level, control, Hessian_kind, WEIGHT = WEIGHT )

!  Check to see if the problem is unbounded from below

        IF ( inform%feasible .AND. FINFO%flag == 4 .AND.                       &
             MAXVAL( ABS( HX( : n ) - RHS( dims%x_s : dims%x_e ) ) ) >         &
               epsmch ** 0.5 ) THEN
          inform%status = GALAHAD_error_unbounded ; GO TO 600 
        END IF

!  If the residual of the linear system is larger than the current 
!  optimality residual, no further progress is likely. Exit

        IF ( SQRT( SUM( ( HX( : dims%v_e ) - RHS ) ** 2 ) ) > merit ) THEN

!  It didn't. We might have run out of options

          IF ( factor == 2 .AND. maxpiv ) THEN
            inform%status = GALAHAD_error_ill_conditioned ; GO TO 600 
              
!  ... or we may change the method

          ELSE IF ( factor < 2 .AND. maxpiv ) THEN
            factor = 2
            pivot_tol = control%pivot_tol
            maxpiv = pivot_tol >= half
            IF ( printi )                                                      &
              WRITE( out, "( ' Switching to augmented system method' )" )

!  Re-analyse the sparsity pattern of the required matrix

            CALL CPU_TIME( time ) 
            CALL LSQP_analyse(                                                 &
              dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C, factor, nnzks,   &
              lk, liw, ldiag_x, ldiag_c_l, ldiag_c_u, IW, Abycol_val,          &
              Abycol_row, Abycol_ptr, K_colptr, DIAG_X, DIAG_C, K,             &
              FACTORS, CNTL, print_level, control, inform )
            CALL CPU_TIME( dum ) ; dum = dum - time
            IF ( printt )                                                      &
              WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
            inform%time%analyse = inform%time%analyse + dum

!  ... or we can increase the pivot tolerance

          ELSE
            maxpiv = .TRUE.
            pivot_tol = half
            IF ( printi ) WRITE( out, "( ' Pivot tolerance increased ' )" )
          END IF
          alpha = zero ; nbact = 0
          CYCLE
        END IF
  
        IF ( printw ) WRITE( out,                                              &
             "( ' ............... step computed ............... ' )" )

        IF ( printd ) THEN 
          WRITE( out, 2120 ) mu 
          WRITE( out, 2100 ) ' DX ', DELTA( dims%x_s : dims%x_e )
          IF ( m > 0 )                                                        &
            WRITE( out, 2100 ) ' DY ', DELTA( dims%y_s : dims%y_e )
        END IF 

!  =======
!  STEP 2:
!  =======

! ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Obtain the search directions for the dual variables
! ::::::::::::::::::::::::::::::::::::::::::::::::::::

!  Problem variables:

!       l = 0
        DO i = dims%x_free + 1, dims%x_l_start - 1
          DZ_l( i ) =   ( mu - Z_l( i ) * ( X( i ) + DELTA( i ) ) ) / X( i )
!         IF ( ABS( one + DELTA( i ) / X( i ) ) < 0.001 ) l = l + 1
        END DO
        DO i = dims%x_l_start, dims%x_l_end
          DZ_l( i ) =   ( mu - Z_l( i ) *                                      &
                           ( DIST_X_l( i ) + DELTA( i ) ) ) / DIST_X_l( i )
!         IF ( ABS( one + DELTA( i ) / DIST_X_l( i ) ) < 0.001 ) l = l + 1
        END DO
  
        DO i = dims%x_u_start, dims%x_u_end
          DZ_u( i ) = - ( mu + Z_u( i ) *                                      &
                          ( DIST_X_u( i ) - DELTA( i ) ) ) / DIST_X_u( i )
!         IF ( ABS( one - DELTA( i ) / DIST_X_u( i ) ) < 0.001 ) l = l + 1
        END DO
  
        DO i = dims%x_u_end + 1, n
          DZ_u( i ) =   ( mu - Z_u( i ) * ( X( i ) + DELTA( i ) ) ) / X( i )
!         IF ( ABS( one + DELTA( i ) / X( i ) )  < 0.001 ) l = l + 1
        END DO
!       write(6,*) l, ' degenerate variable(s)'

!  Slack variables:

        DO i = dims%c_l_start, dims%c_l_end
          DY_l( i ) =   ( mu - Y_l( i ) *                                      &
                    ( DIST_C_l( i ) + DELTA( dims%c_b + i ) ) ) / DIST_C_l( i )
        END DO
  
        DO i = dims%c_u_start, dims%c_u_end
          DY_u( i ) = - ( mu + Y_u( i ) *                                      &
                    ( DIST_C_u( i ) - DELTA( dims%c_b + i ) ) ) / DIST_C_u( i )
        END DO
  
        IF ( printd ) THEN 
          WRITE( out, 2100 ) ' DZ_l ', DZ_l( dims%x_free + 1 : dims%x_l_end )
          WRITE( out, 2100 ) ' DZ_u ', DZ_u( dims%x_u_start : n )
        END IF 

!  Calculate the norm of the search direction

        pmax = MAX( MAXVAL( ABS( DELTA( dims%x_s : dims%x_e ) ) ),             &
                    MAXVAL( ABS( DELTA( dims%c_s : dims%c_e ) ) ),             &
                    MAXVAL( ABS( DELTA( dims%y_s : dims%y_e ) ) ),             &
                    MAXVAL( ABS( DZ_l( dims%x_free + 1 : dims%x_l_end ) ) ),   &
                    MAXVAL( ABS( DZ_u( dims%x_u_start  : n ) ) ),              &
                    MAXVAL( ABS( DY_l( dims%c_l_start  : dims%c_l_end ) ) ),   &
                    MAXVAL( ABS( DY_u( dims%c_u_start  : dims%c_u_end ) ) ) )
  
        IF ( printp ) WRITE( out, 2140 ) pmax 

!  ========
!  STEP 1b:
!  ========

        IF ( control%use_corrector ) THEN

!  =======================================================================
!  -*-*-*-   Obtain the Primal-Dual (Corrector) Search Direction -*-*-*-*-
!  =======================================================================

!  :::::::::::::::::::::::::::::::::
!  Set up the right-hand-side vector
!  :::::::::::::::::::::::::::::::::

!  Problem variables:

          RHS( : dims%x_free ) = zero
          DO i = dims%x_free + 1, dims%x_l_start - 1
            RHS( i ) = - DELTA( i ) * DZ_l( i ) / X( i )
          END DO
          DO i = dims%x_l_start, dims%x_u_start - 1
            RHS( i ) = - DELTA( i ) * DZ_l( i ) / DIST_X_l( i )
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            RHS( i ) = - DELTA( i ) * DZ_l( i ) / DIST_X_l( i )                &
                       + DELTA( i ) * DZ_u( i ) / DIST_X_u( i )
          END DO
          DO i = dims%x_l_end + 1, dims%x_u_end
            RHS( i ) =   DELTA( i ) * DZ_u( i ) / DIST_X_u( i )
          END DO
          DO i = dims%x_u_end + 1, n
            RHS( i ) = - DELTA( i ) * DZ_u( i ) / X( i )
          END DO

!  Slack variables:

          DO i = dims%c_l_start, dims%c_u_start - 1
            RHS( dims%c_b + i ) =                                              &
              - DELTA( dims%c_b + i ) * DY_l( i ) / DIST_C_l( i )
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            RHS( dims%c_b + i ) =                                              &
              - DELTA( dims%c_b + i ) * DY_l( i ) / DIST_C_l( i )              &
              + DELTA( dims%c_b + i ) * DY_u( i ) / DIST_C_u( i )
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            RHS( dims%c_b + i ) =                                              &
                DELTA( dims%c_b + i ) * DY_u( i ) / DIST_C_u( i )
          END DO

!  Include the constraint infeasibilities

          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, m
              RHS( dims%c_e + i ) = zero ; END DO
            DO i = 1, dims%v_e ; DELTA_cor( i ) = RHS( i ) ; END DO
          ELSE
            RHS( dims%y_s : dims%y_e ) = zero
            DELTA_cor = RHS
          END IF
    
          IF ( printd ) THEN 
            WRITE( out, 2100 ) ' RHS_cor_x ', RHS( dims%x_s : dims%x_e )
            IF ( m > 0 )                                                      &
              WRITE( out, 2100 ) ' RHS_cor_y ', RHS( dims%y_s : dims%y_e )
          END IF 

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Obtain the corrector direction for the primal variables
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

          IF ( printw ) WRITE( out,                                            &
               "( ' ............... compute step  ............... ' )" )

!  Use a direct method

          CALL CPU_TIME( time )
          CALL LSQP_iterative_refinement(                                      &
                     dims, n, m, K, FACTORS, CNTL, A_ne, A_val, A_col, A_ptr,  &
                     BARRIER_X, BARRIER_C, DELTA_cor, factor, DIAG_X, ldiag_x, &
                     SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,                    &
                     SOL, RES, BEST, SOL_y, BEST_y, RES_y, RES_x, itref_max,   &
                     print_level, control, inform, Hessian_kind,               &
                     WEIGHT = WEIGHT )
  
          IF ( inform%status /= GALAHAD_ok ) GO TO 700
          CALL CPU_TIME( dum ) ; dum = dum - time
  
          IF ( printt ) WRITE( out, "( ' ** solve time = ', F10.2 ) " ) dum
          inform%time%solve = inform%time%solve + dum

!  Compute the residual of the linear system

          CALL LSQP_residual( dims, n, m, dims%v_e, A_ne, A_val, A_col, A_ptr, &
                              DELTA_cor( dims%x_s : dims%x_e ),                &
                              DELTA_cor( dims%c_s : dims%c_e ),                &
                              DELTA_cor( dims%y_s : dims%y_e ),                &
                              RHS( dims%x_s : dims%x_e ),                      &
                              RHS( dims%c_s : dims%c_e ),                      &
                              RHS( dims%y_s : dims%y_e ), HX( : dims%v_e ),    &
                              BARRIER_X, BARRIER_C,                            &
                              SCALE_C, errorg, errorc, print_level, control,   &
                              Hessian_kind, WEIGHT = WEIGHT )

!  Check to see if the problem is unbounded from below

          IF ( inform%feasible .AND. FINFO%flag == 4 .AND.                     &
               MAXVAL( ABS( HX( : n ) - RHS( dims%x_s : dims%x_e ) ) )         &
                 > epsmch ** 0.5 ) THEN
            inform%status = GALAHAD_error_unbounded ; GO TO 600 
          END IF

!  If the residual of the linear system is larger than the current 
!  optimality residual, no further progress is likely. Exit

          IF ( SQRT( SUM( ( HX( : dims%v_e ) - RHS ) ** 2 ) ) > merit ) THEN

!  It didn't. We might have run out of options

            IF ( factor == 2 .AND. maxpiv ) THEN
              inform%status = GALAHAD_error_ill_conditioned ; GO TO 600 
              
!  ... or we may change the method

            ELSE IF ( factor < 2 .AND. maxpiv ) THEN
              factor = 2
              pivot_tol = control%pivot_tol
              maxpiv = pivot_tol >= half
              IF ( printi )                                                    &
                WRITE( out, "( ' Switching to augmented system method' )" )

!  Re-analyse the sparsity pattern of the required matrix

              CALL CPU_TIME( time ) 
              CALL LSQP_analyse(                                               &
                dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C, factor, nnzks, &
                lk, liw, ldiag_x, ldiag_c_l, ldiag_c_u, IW, Abycol_val,        &
                Abycol_row, Abycol_ptr, K_colptr, DIAG_X, DIAG_C, K,           &
                FACTORS, CNTL, print_level, control, inform )
              CALL CPU_TIME( dum ) ; dum = dum - time
              IF ( printt )                                                    &
                WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
              inform%time%analyse = inform%time%analyse + dum

!  ... or we can increase the pivot tolerance

            ELSE
              maxpiv = .TRUE.
              pivot_tol = half
              IF ( printi ) WRITE( out, "( ' Pivot tolerance increased ' )" )
            END IF
            alpha = zero ; nbact = 0
            CYCLE
          END IF
    
          IF ( printw ) WRITE( out,                                            &
               "( ' ............... step computed ............... ' )" )

          IF ( printd ) THEN 
            WRITE( out, 2120 ) mu 
            WRITE( out, 2100 ) ' DX_cor ', DELTA_cor( dims%x_s : dims%x_e )
            IF ( m > 0 )                                                       &
              WRITE( out, 2100 ) ' DY_cor ', DELTA_cor( dims%y_s : dims%y_e )
          END IF 

!  ========
!  STEP 2b:
!  ========

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Obtain the corrector search directions for the dual variables
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!  Problem variables:

          DO i = dims%x_free + 1, dims%x_l_start - 1
            DZ_cor_l( i ) = - ( DZ_l( i ) * DELTA( i ) +                       &
                                 Z_l( i ) * DELTA_cor( i ) ) / X( i )
          END DO
          DO i = dims%x_l_start, dims%x_l_end
            DZ_cor_l( i ) = - ( DZ_l( i ) * DELTA( i ) +                       &
                                 Z_l( i ) * DELTA_cor( i ) ) / DIST_X_l( i )
          END DO
    
          DO i = dims%x_u_start, dims%x_u_end
            DZ_cor_u( i ) =   ( DZ_u( i ) * DELTA( i ) +                       &
                                 Z_u( i ) * DELTA_cor( i ) ) / DIST_X_u( i )
          END DO
    
          DO i = dims%x_u_end + 1, n
            DZ_cor_u( i ) = - ( DZ_u( i ) * DELTA( i ) +                       &
                                 Z_u( i ) * DELTA_cor( i ) ) / X( i )
          END DO

!  Slack variables:

          DO i = dims%c_l_start, dims%c_l_end
            DY_cor_l( i ) = - ( DY_l( i ) * DELTA( dims%c_b + i ) +            &
              Y_l( i ) * DELTA_cor( dims%c_b + i ) ) / DIST_C_l( i )
          END DO
    
          DO i = dims%c_u_start, dims%c_u_end
            DY_cor_u( i ) =   ( DY_u( i ) * DELTA( dims%c_b + i ) +            &
              Y_u( i ) * DELTA_cor( dims%c_b + i ) ) / DIST_C_u( i )
          END DO
    
          IF ( printd ) THEN 
            WRITE( out, 2100 ) ' DZ_cor_l ',                                   &
              DZ_cor_l( dims%x_free + 1 : dims%x_l_end )
            WRITE( out, 2100 ) ' DZ_cor_u ', DZ_cor_u( dims%x_u_start : n )
          END IF 

!  Calculate the norm of the search direction

          pmax_cor = MAX( MAXVAL( ABS( DELTA_cor( dims%x_s : dims%x_e ) ) ),   &
            MAXVAL( ABS( DELTA_cor( dims%c_s : dims%c_e ) ) ),                 &
            MAXVAL( ABS( DELTA_cor( dims%y_s : dims%y_e ) ) ),                 &
            MAXVAL( ABS( DZ_cor_l( dims%x_free + 1 : dims%x_l_end ) ) ),       &
            MAXVAL( ABS( DZ_cor_u( dims%x_u_start  : n ) ) ),                  &
            MAXVAL( ABS( DY_cor_l( dims%c_l_start  : dims%c_l_end ) ) ),       &
            MAXVAL( ABS( DY_cor_u( dims%c_u_start  : dims%c_u_end ) ) ) )
    
          IF ( printp ) WRITE( out, 2160 ) pmax_cor
        END IF

!  Check to see whether to use a corrector step based on the relative
!  sizes of the predictor and corrector

        IF ( control%use_corrector ) THEN
!         IF ( pmax_cor < ten * pmax ) THEN
            use_corrector = .TRUE.
!         ELSE
!           use_corrector = .FALSE.
!         END IF
        ELSE
          use_corrector = .FALSE.
        END IF

        IF ( use_corrector ) THEN
          co = 'c'
          IF ( printp ) WRITE( out, "( ' ** corrector used' )" )
        ELSE
          co = ' '
        END IF

!  =======
!  STEP 3:
!  =======

!  =====================================================================
!  -*-*-*-*-*-*-*-*-*-*-*-   Line search   -*-*-*-*-*-*-*-*-*-*-*-*-*-*-
!  =====================================================================

!  Perform a line-search to find a point X + alpha DX that
!  sufficiently reduces the merit function

        IF ( printw ) WRITE( out,                                              &
             "( ' .............. get steplength  .............. ' )" )

!  Form the vector H dx + A(trans) dy

        IF ( control%array_syntax_worse_than_do_loop ) THEN
          IF ( Hessian_kind == 0 ) THEN
            DO i = dims%x_s, dims%x_e ; HX( i ) = zero ; END DO
          ELSE IF ( Hessian_kind == 1 ) THEN
            DO i = dims%x_s, dims%x_e ; HX( i ) =  DELTA( i ) ; END DO
          ELSE
            DO i = dims%x_s, dims%x_e
              HX( i ) =  DELTA( i ) * WEIGHT( i ) ** 2 ; END DO
          END IF
        ELSE
          IF ( Hessian_kind == 0 ) THEN
            HX( dims%x_s : dims%x_e ) = zero
          ELSE IF ( Hessian_kind == 1 ) THEN
            HX( dims%x_s : dims%x_e ) = DELTA( dims%x_s : dims%x_e )
          ELSE
            HX( dims%x_s : dims%x_e ) =                                        &
              DELTA( dims%x_s : dims%x_e ) * WEIGHT( dims%x_s : dims%x_e ) ** 2
          END IF
        END IF

        CALL LSQP_AX( n, HX( : n ), m, A_ne, A_val,                            &
                      A_col, A_ptr, m, DELTA( dims%y_s : dims%y_e ), '+T' )

!  Perform a line-search to find a point X + alpha DX which
!  sufficiently reduces the measure of potential

        IF ( inform%feasible .AND. Hessian_kind == 0 .AND.                     &
             gradient_kind == 0 ) THEN

!  Find the largest possible feasible stepsize

          alpha = infinity
          DO i = dims%x_free + 1, dims%x_l_start - 1
            IF ( DELTA( i ) < zero ) alpha = MIN( alpha, - X( i ) / DELTA( i ) )
          END DO

          DO i = dims%x_l_start, dims%x_l_end
            IF ( DELTA( i ) < zero )                                           &
              alpha = MIN( alpha, - DIST_X_l( i ) / DELTA( i ) )
          END DO 
  
          DO i = dims%x_u_start, dims%x_u_end
            IF ( DELTA( i ) > zero )                                           &
              alpha = MIN( alpha, DIST_X_u( i ) / DELTA( i ) )
          END DO 

          DO i = dims%x_u_end + 1, n
            IF ( DELTA( i ) > zero ) alpha = MIN( alpha, - X( i ) / DELTA( i ) )
          END DO

          DO i = dims%c_l_start, dims%c_l_end
            IF ( DELTA(  dims%c_b + i ) < zero )                               &
              alpha = MIN( alpha, - DIST_C_l( i ) / DELTA( dims%c_b + i ) )
          END DO 
    
          DO i = dims%c_u_start, dims%c_u_end
            IF ( DELTA(  dims%c_b + i ) > zero )                               &
              alpha = MIN( alpha, DIST_C_u( i ) / DELTA( dims%c_b + i ) )
          END DO 

!  A step of no larger than one will be attempted

          alpha = MIN( one, 0.9999_wp * alpha )

          IF ( printt ) WRITE( out, "( ' ',/,'       ***  Linesearch      ',   &
         &                             '  step trial centr model centr ' )" )

          nbact = 0 ; step = alpha
          DO

!  Calculate the distances to the bounds at the trial point

            IF ( control%array_syntax_worse_than_do_loop ) THEN
              DO i = 1, n ; X( i ) = X( i ) + step * DELTA( i ) ; END DO
            ELSE
              X = X + step * DELTA( dims%x_s : dims%x_e )
            END IF

            DO i = dims%x_l_start, dims%x_l_end
              DIST_X_l( i ) = DIST_X_l( i ) + step * DELTA( i ) 
            END DO 
    
            DO i = dims%x_u_start, dims%x_u_end
              DIST_X_u( i ) = DIST_X_u( i ) - step * DELTA( i ) 
            END DO 
    
!  Do the same for the slacks

            DO i = dims%c_l_start, dims%c_l_end
              DIST_C_l( i ) = DIST_C_l( i ) + step * DELTA( dims%c_b + i ) 
            END DO 
    
            DO i = dims%c_u_start, dims%c_u_end
              DIST_C_u( i ) = DIST_C_u( i ) - step * DELTA( dims%c_b + i ) 
            END DO 

!  Ensure that the measure of potential has decreased

            one_minus_alpha = one - alpha
            potential_trial = LSQP_potential_value( dims, n,                   &
                                   X, DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u )

!  Check to see if the Amijo criterion is satisfied. If not, halve the 
!  steplength

            IF ( printt ) WRITE( out, "( 22X, 3ES12.4 )" )                     &
              alpha, potential_trial, inform%potential

            IF ( potential_trial <= inform%potential ) EXIT
            alpha = alpha * half ;  step = - alpha ; nbact = nbact + 1 
            IF ( alpha < epsmch ) THEN
              inform%status = GALAHAD_error_tiny_step
              GO TO 500 
            END IF
          END DO
    
!  Update the Lagrange multipliers

          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, m 
              Y( i ) = Y( i ) - alpha * DELTA( dims%c_e + i ) ; END DO
          ELSE
            Y = Y - alpha * DELTA( dims%y_s : dims%y_e )
          END IF

!  Calculate the new dual variables
          
          DO i = dims%x_free + 1, dims%x_l_start - 1
            Z_l( i ) = mu / X( i )
          END DO

          DO i = dims%x_l_start, dims%x_l_end
            Z_l( i ) = mu / DIST_X_l( i )
          END DO 
  
          DO i = dims%x_u_start, dims%x_u_end
            Z_u( i ) = - mu / DIST_X_u( i )
          END DO 

          DO i = dims%x_u_end + 1, n
            Z_u( i ) = mu / X( i )
          END DO

!  Do the same for the Lagrange multipliers

          DO i = dims%c_l_start, dims%c_l_end
            Y_l( i ) = mu / DIST_C_l( i )
          END DO 
    
          DO i = dims%c_u_start, dims%c_u_end
            Y_u( i ) = - mu / DIST_C_u( i )
          END DO 

          inform%potential = potential_trial
          inform%nbacts = inform%nbacts + nbact 

!  Perform a line-search to find a point X + alpha DX that
!  sufficiently reduces the merit function

        ELSE

!  Find the maximum step which keeps the complementarity and feasibilty balanced

          CALL LSQP_compute_maxstep( dims, n, m, Z_l, Z_u, DZ_l, DZ_u,         &
                                     X, DIST_X_l, DIST_X_u,                    &
                                     DELTA( dims%x_s : dims%x_e ),             &
                                     Y_l, Y_u, DY_l, DY_u, DIST_C_l, DIST_C_u, &
                                     DELTA( dims%c_s : dims%c_e ), gamma_b,    &
                                     gamma_f, nu, nbnds,                       &
                                     alpha_b, alpha_f, print_level,            &
                                     control, inform )

          IF ( inform%status /= GALAHAD_ok ) GO TO 500
          IF ( printt ) WRITE( out, "( ' alpha_b, alpha_f = ', 2ES12.4 )" )    &
                                         alpha_b, alpha_f

!  A step of no larger than one will be attempted

          alpha = MIN( one, alpha_b, alpha_f )

! ::::::::::::::::::::::::::::::::::::::::::::
!  Record the slope along the search direction
! ::::::::::::::::::::::::::::::::::::::::::::

          slope = - ( merit - mu * nbnds )
          IF ( printt ) WRITE( out, "( '  Value and slope = ',1P,2D12.4 )" )   &
            merit, slope

! ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Use a backtracking line-search, starting from alpha
! ::::::::::::::::::::::::::::::::::::::::::::::::::::

          IF ( printt ) WRITE( out, "( ' ',/,'       ***  Linesearch    ',     &
         &                             '    step trial value model value ' )" )

          nbact = 0 ; step = alpha
          DO

!  The sqaure of the norm of the new residual should be smaller than a 
!  linear model

            merit_model = merit + alpha * eta * slope 
 
!  Calculate the distances to the bounds and the dual variables at the
!  trial point

            IF ( control%array_syntax_worse_than_do_loop ) THEN
              DO i = 1, n ; X( i ) = X( i ) + step * DELTA( i ) ; END DO
              DO i = 1, m 
                Y( i ) = Y( i ) - step * DELTA( dims%c_e + i ) ; END DO
            ELSE
              X = X + step * DELTA( dims%x_s : dims%x_e )
              Y = Y - step * DELTA( dims%y_s : dims%y_e )
            END IF

            DO i = dims%x_free + 1, dims%x_l_start - 1
              Z_l( i ) = Z_l( i ) + step * DZ_l( i )
            END DO 
  
            DO i = dims%x_l_start, dims%x_l_end
              DIST_X_l( i ) = DIST_X_l( i ) + step * DELTA( i ) 
              Z_l( i ) = Z_l( i ) + step * DZ_l( i )
            END DO 
    
            DO i = dims%x_u_start, dims%x_u_end
              DIST_X_u( i ) = DIST_X_u( i ) - step * DELTA( i ) 
              Z_u( i ) = Z_u( i ) + step * DZ_u( i )
            END DO 
    
            DO i = dims%x_u_end + 1, n
              Z_u( i ) = Z_u( i ) + step * DZ_u( i )
            END DO 

!  Do the same for the slacks and their duals

            DO i = dims%c_l_start, dims%c_l_end
              DIST_C_l( i ) = DIST_C_l( i ) + step * DELTA( dims%c_b + i ) 
              Y_l( i ) = Y_l( i ) + step * DY_l( i )
            END DO 
    
            DO i = dims%c_u_start, dims%c_u_end
              DIST_C_u( i ) = DIST_C_u( i ) - step * DELTA( dims%c_b + i ) 
              Y_u( i ) = Y_u( i ) + step * DY_u( i )
            END DO 

!  Evaluate the merit function at the new point

            one_minus_alpha = one - alpha
            merit_trial = LSQP_merit_value(                                    &
                            dims, n, m, X, Y, Y_l, Y_u, Z_l, Z_u,              &
                            DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u,            &
                            GRAD_L( dims%x_s : dims%x_e ) + alpha * HX( : n ), &
                            one_minus_alpha * C_RES, res_dual )
            IF ( printt ) WRITE( out, "( 22X, 3ES12.4 )" )                     &
              alpha, merit_trial, merit_model 

!  Check to see if the Amijo criterion is satisfied. If not, halve the 
!  steplength

            IF ( merit_trial <= merit_model ) EXIT
            alpha = alpha * half ;  step = - alpha ; nbact = nbact + 1 
            IF ( alpha < epsmch ) THEN
              IF ( inform%iter - 1 > muzero_fixed ) THEN
                inform%status = GALAHAD_error_tiny_step
                GO TO 500 
              ELSE
                muzero_fixed = inform%iter - 2
                EXIT
              END IF
            END IF
          END DO
          merit = merit_trial 
    
          inform%nbacts = inform%nbacts + nbact 
        END IF

!  Update the slack variables 

        IF ( use_corrector ) THEN
          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 0, dims%nc - 1
              C( dims%c_l_start + i ) = C( dims%c_l_start + i ) +              &
                alpha * ( DELTA( dims%c_s + i ) +                              &
                  alpha * DELTA_cor( dims%c_s + i ) )
            END DO
          ELSE
            C = C + alpha * ( DELTA( dims%c_s : dims%c_e ) +                   &
                  alpha * DELTA_cor( dims%c_s : dims%c_e ) )
          END IF
        ELSE
          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 0, dims%nc - 1
              C( dims%c_l_start + i ) = C( dims%c_l_start + i ) +              &
                alpha * DELTA( dims%c_s + i ) ; END DO
          ELSE
            C = C + alpha * DELTA( dims%c_s : dims%c_e )
          END IF
        END IF

!  Update the values of the merit function, the gradient of the Lagrangian, 
!  and the constraint residuals

        GRAD_L( dims%x_s : dims%x_e ) = GRAD_L( dims%x_s : dims%x_e ) +        &
          alpha * HX( dims%x_s : dims%x_e )

        IF ( use_corrector ) THEN
          DO i = dims%x_free + 1, dims%x_l_end
            GRAD_L( i ) = GRAD_L( i ) + alpha * alpha * DZ_cor_l( i )
          END DO 
          DO i = dims%x_u_start, n
            GRAD_L( i ) = GRAD_L( i ) + alpha * alpha * DZ_cor_u( i )
          END DO 
        END IF

        C_RES = one_minus_alpha * C_RES

!       WRITE( 6, "(' updated  grad_l ', ES12.4 )" )                           &
!         MAXVAL( ABS( GRAD_L( dims%x_s : dims%x_e ) ) )

!       CALL LSQP_Lagrangian_gradient( dims, n, m, X, Y, Y_l, Y_u, Z_l, Z_u,   &
!                                    A_ne, A_val, A_col, A_ptr,                &
!                                    DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u,   &
!                                    GRAD_L( dims%x_s : dims%x_e ),            &
!                                    control%getdua, dufeas,                   &
!                                    Hessian_kind, gradient_kind, control,     &
!                                    WEIGHT = WEIGHT, X0 = X0, G = G )

!        WRITE( 6, "(' computed grad_l ', ES12.4 )" )                          &
!         MAXVAL( ABS( GRAD_L( dims%x_s : dims%x_e ) ) )

!  Update the norm of the constraint residual

        res_prim = one_minus_alpha * res_prim
        nu = one_minus_alpha * nu

!  Evaluate the merit function if not already done

        IF ( inform%feasible .AND. Hessian_kind == 0 .AND.                     &
                  gradient_kind == 0 ) THEN
          merit = LSQP_merit_value( dims, n, m, X, Y, Y_l, Y_u, Z_l, Z_u,      &
                            DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u,            &
                            GRAD_L( dims%x_s : dims%x_e ), C_RES, res_dual )
        END IF

!  Compute the objective function value

        IF ( Hessian_kind == 0 ) THEN
          inform%obj = f
        ELSE IF ( Hessian_kind == 1 ) THEN
          inform%obj = f + half * SUM( ( X - X0 ) ** 2 )
        ELSE
          inform%obj = f + half * SUM( ( WEIGHT * ( X - X0 ) ) ** 2 )
        END IF
  
        IF ( gradient_kind == 1 ) THEN
          inform%obj = inform%obj + SUM( X )
        ELSE IF ( gradient_kind /= 0 ) THEN
          inform%obj = inform%obj + DOT_PRODUCT( G, X )
        END IF

        IF ( m > 0 ) THEN 
          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, dims%c_equality ; C_RES( i ) = - C_l( i ) ; END DO
            DO i = dims%c_l_start, dims%c_u_end 
              C_RES( i ) = - SCALE_C( i ) * C( i ) ; END DO
          ELSE
            C_RES( : dims%c_equality ) = - C_l( : dims%c_equality ) 
            C_RES( dims%c_l_start : dims%c_u_end ) = - SCALE_C * C
          END IF
          CALL LSQP_AX( m, C_RES, m, A_ne, A_val, A_col, A_ptr,      &
                        n, X, '+ ' )
          IF ( printt ) WRITE( out, "( '  Constraint residual ', ES12.4 )" )   &
            MAXVAL( ABS( C_RES ) )
!         WRITE( 6, "( ' rec, cal cres = ', 2ES12.4 )" )                       &
!           res_prim, MAXVAL( ABS( C_RES ) )
          IF ( res_prim < MAXVAL( ABS( C_RES ) ) ) THEN
            res_prim = MAXVAL( ABS( C_RES ) )
            IF ( res_prim < one .AND. itref_max < control%itref_max + 2 ) THEN
              IF ( printt ) WRITE( out, "( ' Increasing # refinements .. ' )" )
              itref_max = itref_max + 1
            END IF
          END IF
        END IF

!  Compute the complementary slackness, and the min/max components
!  of the primal/dual infeasibilities

!       DO i = dims%x_free + 1, dims%x_l_start - 1
!         write(6,"(I6, ' x lower', 2ES12.4)" ) i, X( i ), Z_l( i )
!       END DO
!       DO i = dims%x_l_start, dims%x_l_end
!         write(6,"(I6, ' x lower', 2ES12.4)" ) i, DIST_X_l( i ), Z_l( i )
!       END DO
!       DO i = dims%x_u_start, dims%x_u_end
!         write(6,"(I6, ' x upper', 2ES12.4)" ) i, - DIST_X_u( i ), Z_u( i )
!       END DO
!       DO i = dims%x_u_end + 1, n
!         write(6,"(I6, ' x upper', 2ES12.4)" ) i, X( i ), Z_u( i )
!       END DO

!       DO i = dims%c_l_start, dims%c_l_end
!         write(6,"(I6, ' c lower', 2ES12.4)" ) i, DIST_C_l( i ), Y_l( i )
!       END DO
!       DO i = dims%c_u_start, dims%c_u_end
!         write(6,"(I6, ' c upper', 2ES12.4)" ) i, - DIST_C_u( i ), Y_u( i )
!       END DO

        slknes_x = DOT_PRODUCT( X( dims%x_free + 1 : dims%x_l_start - 1 ),     &
                                Z_l( dims%x_free + 1 : dims%x_l_start - 1 ) ) +&
                   DOT_PRODUCT( DIST_X_l( dims%x_l_start : dims%x_l_end ),     &
                                Z_l( dims%x_l_start : dims%x_l_end ) ) -       &
                   DOT_PRODUCT( DIST_X_u( dims%x_u_start : dims%x_u_end ),     &
                                Z_u( dims%x_u_start : dims%x_u_end ) ) +       &
                   DOT_PRODUCT( X( dims%x_u_end + 1 : n ),                     &
                              Z_u( dims%x_u_end + 1 : n ) )
        slknes_c = DOT_PRODUCT( DIST_C_l( dims%c_l_start : dims%c_l_end ),     &
                                Y_l( dims%c_l_start : dims%c_l_end ) ) -       &
                   DOT_PRODUCT( DIST_C_u( dims%c_u_start : dims%c_u_end ),     &
                                Y_u( dims%c_u_start : dims%c_u_end ) )
        slknes = slknes_x + slknes_c
  
        slkmin_x = MIN( MINVAL( X( dims%x_free + 1 : dims%x_l_start - 1 ) *    &
                                Z_l( dims%x_free + 1 : dims%x_l_start - 1 ) ), &
                        MINVAL( DIST_X_l( dims%x_l_start : dims%x_l_end ) *    &
                                Z_l( dims%x_l_start : dims%x_l_end ) ),        &
                        MINVAL( - DIST_X_u( dims%x_u_start : dims%x_u_end ) *  &
                                Z_u( dims%x_u_start : dims%x_u_end ) ),        &
                        MINVAL( X( dims%x_u_end + 1 : n ) *                    &
                                Z_u( dims%x_u_end + 1 : n ) ) )
        slkmin_c = MIN( MINVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) *    &
                                Y_l( dims%c_l_start : dims%c_l_end ) ),        &
                        MINVAL( - DIST_C_u( dims%c_u_start : dims%c_u_end ) *  &
                                Y_u( dims%c_u_start : dims%c_u_end ) ) )
        slkmin = MIN( slkmin_x, slkmin_c )
  
        slkmax_x = MAX( MAXVAL( X( dims%x_free + 1 : dims%x_l_start - 1 ) *    &
                                Z_l( dims%x_free + 1 : dims%x_l_start - 1 ) ), &
                        MAXVAL( DIST_X_l( dims%x_l_start : dims%x_l_end ) *    &
                                Z_l( dims%x_l_start : dims%x_l_end ) ),        &
                        MAXVAL( - DIST_X_u( dims%x_u_start : dims%x_u_end ) *  &
                                Z_u( dims%x_u_start : dims%x_u_end ) ),        &
                        MAXVAL( X( dims%x_u_end + 1 : n ) *                    &
                                Z_u( dims%x_u_end + 1 : n ) ) )
        slkmax_c = MAX( MAXVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) *    &
                                Y_l( dims%c_l_start : dims%c_l_end ) ),        &
                        MAXVAL( - DIST_C_u( dims%c_u_start : dims%c_u_end ) *  &
                                Y_u( dims%c_u_start : dims%c_u_end ) ) )
  
        p_min = MIN( MINVAL( X( dims%x_free + 1 : dims%x_l_start - 1 ) ),      &
                     MINVAL( DIST_X_l( dims%x_l_start : dims%x_l_end ) ),      &
                     MINVAL( DIST_X_u( dims%x_u_start : dims%x_u_end ) ),      &
                     MINVAL( - X( dims%x_u_end + 1 : n ) ),                    &
                     MINVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) ),      &
                     MINVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) ) )
  
        p_max = MAX( MAXVAL( X( dims%x_free + 1 : dims%x_l_start - 1 ) ),      &
                     MAXVAL( DIST_X_l( dims%x_l_start : dims%x_l_end ) ),      &
                     MAXVAL( DIST_X_u( dims%x_u_start : dims%x_u_end ) ),      &
                     MAXVAL( - X( dims%x_u_end + 1 : n ) ),                    &
                     MAXVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) ),      &
                     MAXVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) ) )
  
        d_min = MIN( MINVAL(   Z_l( dims%x_free + 1 : dims%x_l_end ) ),        &
                     MINVAL( - Z_u( dims%x_u_start : n ) ),                    &
                     MINVAL(   Y_l( dims%c_l_start : dims%c_l_end ) ),         &
                     MINVAL( - Y_u( dims%c_u_start : dims%c_u_end ) ) )
  
        d_max = MAX( MAXVAL(   Z_l( dims%x_free + 1 : dims%x_l_end ) ),        &
                     MAXVAL( - Z_u( dims%x_u_start : n ) ),                    &
                     MAXVAL(   Y_l( dims%c_l_start : dims%c_l_end ) ),         &
                     MAXVAL( - Y_u( dims%c_u_start : dims%c_u_end ) ) )
  
        IF ( nbnds_x > 0 ) THEN
          slknes_x = slknes_x / nbnds_x
        ELSE
          slknes_x = zero
        END IF
  
        IF ( nbnds_c > 0 ) THEN
          slknes_c = slknes_c / nbnds_c
        ELSE
          slknes_c = zero
        END IF
        IF ( nbnds > 0 ) THEN
          slknes = slknes / nbnds
          slknes_req = slknes
        ELSE
          slknes = zero
        END IF

        IF ( printt .AND. nbnds > 0 ) WRITE( out, 2130 )                       &
          slknes, slknes_x, slknes_c, slkmin_x, slkmax_x, slkmin_c, slkmax_c,  &
          p_min, p_max, d_min, d_max

        IF ( printd ) THEN
          WRITE( out, "( ' primal-dual -vs- primal dual variables ' )" )
          WRITE( out, "( ' lower ', /, ( 2( I6, 2ES12.4 ) ) )" )               &
            ( i, Z_l( i ), mu / X( i ),                                        &
              i = dims%x_free + 1, dims%x_l_start - 1 ),                       &
            ( i, Z_l( i ), mu / DIST_X_l( i ),                                 &
              i =  dims%x_l_start, dims%x_l_end ) 
          WRITE( out, "( ' upper ', /, ( 2( I6, 2ES12.4 ) ) )" )               &
            ( i, Z_u( i ),  - mu / DIST_X_u( i ),                              &
              i = dims%x_u_start, dims%x_u_end ),                              &
            ( i, Z_u( i ), mu / X( i ),                                        &
              i = dims%x_u_end + 1, n )
        END IF

!  Test to see if we are feasible

        IF ( res_prim <= control%stop_p ) THEN
          IF ( control%just_feasible ) THEN
            inform%status = GALAHAD_ok
            inform%feasible = .TRUE.
            IF ( printi ) THEN
              CALL CPU_TIME( time ) ; inform%time%total = time - time_start
              WRITE( out, 2070 )
              WRITE( out, 2030 ) inform%iter, re, res_prim, res_dual,          &
                slknes_req, zero, alpha, co, mo, mu, nbact, time 
              IF ( printt ) WRITE( out, 2000 ) 
            END IF
            GO TO 500
          END IF

          IF ( .NOT. inform%feasible ) THEN
            IF ( printi ) WRITE( out, 2070 )
            inform%feasible = .TRUE.
            IF ( Hessian_kind == 0 .AND. gradient_kind == 0 ) THEN
              IF ( slkmin_x >= epsmch .AND. slkmin_c >= epsmch ) THEN
                inform%potential = LSQP_potential_value( dims, n, X, DIST_X_l, &
                                                  DIST_X_u, DIST_C_l, DIST_C_u )
              ELSE
                inform%potential = infinity
              END IF
            END IF
          END IF
  
          res_prim = zero ; nu = zero
          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, m ; C_RES( i ) = zero ; END DO
          ELSE
            C_RES = zero
          END IF
        END IF

!  =======
!  STEP 5:
!  =======

        IF ( get_stat ) THEN

!  Estimate the variable and constraint exit status

          CALL LSQP_indicators( dims, n, m, C_l, C_u, C_last, C,               &
                                DIST_C_l, DIST_C_u, X_l, X_u, X_last, X,       &
                                DIST_X_l, DIST_X_u, Y_l, Y_u, Z_l, Z_u,        &
                                Y_last, Z_last,                                &
                                control, C_stat = C_stat, B_stat = B_stat )

!  Count the number of active constraints/bounds

          IF ( printt )                                                        &
            WRITE( out, "( ' indicators: n_active/n, m_active/m ', 4I7 )" )    &
               COUNT( B_stat /= 0 ), n, COUNT( C_stat /= 0 ), m
        END IF

!  Compute the new penalty parameter

        sigma = sigma_max
!       mu = sigma * slknes
        IF ( inform%iter > muzero_fixed ) mu = sigma * slknes
        IF ( mu < one .AND. stat_required ) THEN
          get_stat = .TRUE. 
          C_last( dims%c_l_start : dims%c_u_end )                            &
            = C( dims%c_l_start : dims%c_u_end )
          X_last = X

          DO i = dims%c_l_start, dims%c_u_start - 1
            Y_last( i ) = Y_l( i )
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            IF ( DIST_C_l( i ) <= DIST_C_u( i ) ) THEN
              Y_last( i ) = Y_l( i )
            ELSE
              Y_last( i ) = Y_u( i )
            END IF
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            Y_last( i ) = Y_u( i )
          END DO

          Z_last( : dims%x_free ) = zero
          DO i = dims%x_free + 1, dims%x_u_start - 1
            Z_last( i ) = Z_l( i )
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            IF ( DIST_X_l( i ) <= DIST_X_u( i ) ) THEN
              Z_last( i ) = Z_l( i )
            ELSE
              Z_last( i ) = Z_u( i )
            END IF
          END DO
          DO i = dims%x_l_end + 1, n
            Z_last( i ) = Z_u( i )
          END DO
        END IF       

!  Compute the projected gradient of the Lagrangian function

        pjgnrm = zero 
        DO i = 1, n 
          gi = GRAD_L( i ) 
          IF ( gi < zero ) THEN 
            gi = - MIN( ABS( X_u( i ) - X( i ) ), - gi ) 
          ELSE 
            gi = MIN( ABS( X_l( i ) - X( i ) ), gi ) 
          END IF 
          pjgnrm = MAX( pjgnrm, ABS( gi ) ) 
        END DO 
  
        IF ( printd ) THEN 
          WRITE( out, 2100 ) ' DIST_X_l ',                                     &
            X( dims%x_free + 1 : dims%x_l_start - 1 ), DIST_X_l
          WRITE( out, 2100 ) ' DIST_X_u ',                                     &
            DIST_X_u, - X( dims%x_u_end + 1 : n )
          WRITE( out, "( ' ' )" ) 
        END IF 
  
        IF ( printd ) WRITE( out, 2110 ) pjgnrm, res_prim
      END DO 

!  ---------------------------------------------------------------------
!  ---------------------- End of Major Iteration -----------------------
!  ---------------------------------------------------------------------

  500 CONTINUE 

!  Print details of the solution obtained

  600 CONTINUE 

!  Compute the final objective function value

      IF ( Hessian_kind == 0 ) THEN
        inform%obj = f
      ELSE IF ( Hessian_kind == 1 ) THEN
        inform%obj = f + half * SUM( ( X - X0 ) ** 2 )
      ELSE
        inform%obj = f + half * SUM( ( WEIGHT * ( X - X0 ) ) ** 2 )
      END IF

      IF ( gradient_kind == 1 ) THEN
        inform%obj = inform%obj + SUM( X )
      ELSE IF ( gradient_kind /= 0 ) THEN
        inform%obj = inform%obj + DOT_PRODUCT( G, X )
      END IF

      IF ( printi ) THEN
        WRITE( out, 2010 ) inform%obj, inform%iter, inform%nbacts
        WRITE( out, 2110 ) pjgnrm, res_prim
        IF ( control%getdua ) WRITE( out, 2090 ) 

        IF ( factor == 0 .OR. factor == 1 ) THEN
          WRITE( out, "( /, '  Schur-complement method used ' )" )
        ELSE
          WRITE( out, "( /, '  Augmented system method used ' )" )
        END IF
      END IF

!  If required, make the solution exactly complementary

      IF ( control%feasol ) THEN 
        DO i = dims%x_free + 1, dims%x_l_start - 1
          IF ( ABS( Z_l( i ) ) < ABS( X( i ) ) ) THEN 
            Z_l( i ) = zero 
          ELSE 
            X( i ) = X_l( i ) 
          END IF 
        END DO

        DO i = dims%x_l_start, dims%x_l_end
          IF ( ABS( Z_l( i ) ) < ABS( DIST_X_l( i ) ) ) THEN 
            Z_l( i ) = zero 
          ELSE 
            X( i ) = X_l( i ) 
          END IF 
        END DO

        DO i = dims%x_u_start, dims%x_u_end
          IF ( ABS( Z_u( i ) ) < ABS( DIST_X_u( i ) ) ) THEN 
            Z_u( i ) = zero 
          ELSE 
            X( i ) = X_u( i ) 
          END IF 
        END DO

        DO i = dims%x_u_end + 1, n
          IF ( ABS( Z_u( i ) ) < ABS( X( i ) ) ) THEN 
            Z_u( i ) = zero 
          ELSE 
            X( i ) = X_u( i ) 
          END IF 
        END DO
      END IF 

!  Exit

 700  CONTINUE

!  Set the dual variables

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, dims%x_free ; Z( i ) = zero ; END DO
      ELSE
        Z( : dims%x_free ) = zero
      END IF

      DO i = dims%x_free + 1, dims%x_u_start - 1
        Z( i ) = Z_l( i )
      END DO

      DO i = dims%x_u_start, dims%x_l_end
        IF ( ABS( Z_l( i ) ) <= ABS( Z_u( i ) ) ) THEN
          Z( i ) = Z_u( i )
        ELSE
          Z( i ) = Z_l( i )
        END IF
      END DO

      DO i = dims%x_l_end + 1, n
        Z( i ) = Z_u( i )
      END DO

!  Unscale the constraint bounds

      DO i = dims%c_l_start, dims%c_l_end
        C_l( i ) = C_l( i ) * SCALE_C( i )
      END DO

      DO i = dims%c_u_start, dims%c_u_end
        C_u( i ) = C_u( i ) * SCALE_C( i )
      END DO

!  Compute the values of the constraints

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, m ; C_RES( i ) = zero ; END DO
      ELSE
        C_RES( : m ) = zero
      END IF
      CALL LSQP_AX( m, C_RES( : m ), m, A_ne, A_val, A_col,                    &
                    A_ptr, n, X, '+ ')
      IF ( printi .AND. m > 0 ) THEN 
        WRITE( out, "( '  Constraint residual ', ES12.4 )" )                   &
             MAX( zero, MAXVAL( ABS( C_l( : dims%c_equality ) -                &
                                     C_RES(: dims%c_equality ) ) ),            &
                        MAXVAL( C_l(  dims%c_l_start : dims%c_l_end ) -        &
                                C_RES(  dims%c_l_start : dims%c_l_end ) ),     &
                        MAXVAL( C_RES( dims%c_u_start : dims%c_u_end ) -       &
                                C_u( dims%c_u_start : dims%c_u_end ) ) )     
      END IF

!  If necessary, print warning messages

  810 CONTINUE
      IF ( printi ) then

        SELECT CASE( inform%status )
          CASE( - 1  ) ; WRITE( out, 2210 ) 
          CASE( GALAHAD_error_no_center ) ; WRITE( out, 2240 ) 
          CASE( GALAHAD_error_bad_bounds ) ; WRITE( out, 2250 ) 
          CASE( GALAHAD_error_primal_infeasible ) ; WRITE( out, 2260 ) 
          CASE( GALAHAD_error_factorization ) ; WRITE( out, 2270 ) 
          CASE( GALAHAD_error_ill_conditioned ) ; WRITE( out, 2280 ) 
          CASE( GALAHAD_error_tiny_step ) ; WRITE( out, 2290 ) 
          CASE( GALAHAD_error_max_iterations ) ; WRITE( out, 2300 ) 
          CASE( GALAHAD_error_unbounded ) ; WRITE( out, 2310 ) 
        END SELECT

      END IF
      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' leaving LSQP_solve_main ' )" )

      RETURN  

!  Non-executable statements

 2000 FORMAT( /,' Iter   p-feas  d-feas com-slk   obj    ',                    &
                '  step      mu    bac    time' ) 
 2010 FORMAT( //, '  Final objective function value ', ES22.14,                &
              /,  '  Total number of iterations = ', I6,                       &
              /,  '  Total number of backtracks = ', I6 )
 2020 FORMAT( I5, A1, 3ES8.1, ES9.1, '     -    ', ES7.1,                      &
            '   -', 0P, F9.2 ) 
 2030 FORMAT( I5, A1, 3ES8.1, ES9.1, ES8.1, 2A1, ES7.1, I4, 0P, F9.2 ) 
 2040 FORMAT( '   **  Error return ', I6, ' from ', A15 ) 
 2050 FORMAT( '   **  Warning ', I6, ' from ', A15 ) 
 2060 FORMAT( I8, ' integer and ', I8, ' real words needed for factorization' )
 2070 FORMAT( /, ' ====================== feasible point found',               &
                 ' ======================= ', / )
 2090 FORMAT( /, ' Advanced starting point used for dual variables ' ) 
 2100 FORMAT( A10, 7ES10.2, /, ( 10X, 7ES10.2 ) ) 
 2110 FORMAT( /,'  Norm of projected gradient is ', ES12.4, /,                 &
                '  Norm of infeasibility is      ', ES12.4 ) 
 2120 FORMAT( ' Penalty parameter = ', ES12.4 ) 
 2130 FORMAT( 21X, ' == >  mu estimated   = ', ES9.1, /,                       &
              21X, '       mu_x estimated = ', ES9.1, /,                       &
              21X, '       mu_c estimated = ', ES9.1, /,                       &
              21X, '  min/max slackness_x = ', 2ES12.4, /,                     &
              21X, '  min/max slackness_c = ', 2ES12.4, /,                     &
              21X, '  min/max primal feasibility = ', 2ES12.4, /,              &
              21X, '  min/max dual   feasibility = ', 2ES12.4 ) 
 2140 FORMAT( /, '  Norm of (predictor) search direction = ', ES12.4 ) 
 2150 FORMAT( A6, /, ( 4( 2I5, ES10.2 ) ) )
 2160 FORMAT( /, '  Norm of (corrector) search direction = ', ES12.4 ) 
 2210 FORMAT( /, '  Warning - input paramters incorrect ' ) 
 2240 FORMAT( /, '  Warning - the analytic centre appears to be unbounded ' ) 
 2250 FORMAT( /, '  Warning - the constraints are inconsistent ' ) 
 2260 FORMAT( /, '  Warning - the constraints appear to be inconsistent ' ) 
 2270 FORMAT( /, '  Warning - factorization failure ' ) 
 2280 FORMAT( /, '  Warning - no further progress possible ' ) 
 2290 FORMAT( /, '  Warning - step too small to make further progress ' ) 
 2300 FORMAT( /, '  Warning - iteration bound exceeded ' ) 
 2310 FORMAT( /, '  Warning - problem unbounded from below ' ) 

!  End of LSQP_solve_main

      END SUBROUTINE LSQP_solve_main

!-*-*-*-*-*-*-   L S Q P _ T E R M I N A T E   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE LSQP_terminate( data, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!      ..............................................
!      .                                            .
!      .  Deallocate internal arrays at the end     .
!      .  of the computation                        .
!      .                                            .
!      ..............................................

!  Arguments:
!
!   data    see Subroutine LSQP_initialize
!   control see Subroutine LSQP_initialize
!   inform  see Subroutine LSQP_solve

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_data_type ), INTENT( INOUT ) :: data
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      TYPE ( LSQP_inform_type ), INTENT( INOUT ) :: inform
 
      inform%status = GALAHAD_ok

!  Deallocate all allocated integer arrays

      IF ( ALLOCATED( data%INDEX_C_freed ) ) THEN
        DEALLOCATE( data%INDEX_C_freed, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%INDEX_C_freed'
          IF ( control%error > 0 ) WRITE( control%error, 2900 )                &
            inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%RES_x ) ) THEN
        DEALLOCATE( data%RES_x, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%RES_x'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%SOL_y ) ) THEN
        DEALLOCATE( data%SOL_y, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%SOL_y'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%RES_y ) ) THEN
        DEALLOCATE( data%RES_y, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%RES_y'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%BEST_y ) ) THEN
        DEALLOCATE( data%BEST_y, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%BEST_y'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%SOL ) ) THEN
        DEALLOCATE( data%SOL, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%SOL'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%RES ) ) THEN
        DEALLOCATE( data%RES, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%RES'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%BEST ) ) THEN
        DEALLOCATE( data%BEST, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%BEST'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%HX ) ) THEN
        DEALLOCATE( data%HX, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%HX'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%GRAD_L ) ) THEN
        DEALLOCATE( data%GRAD_L, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%GRAD_L'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DIST_X_l ) ) THEN
        DEALLOCATE( data%DIST_X_l, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DIST_X_l'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DIST_X_u ) ) THEN
        DEALLOCATE( data%DIST_X_u, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DIST_X_u'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%Z_l ) ) THEN
        DEALLOCATE( data%Z_l, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%Z_l'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%Z_u ) ) THEN
        DEALLOCATE( data%Z_u, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%Z_u'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%BARRIER_X ) ) THEN
        DEALLOCATE( data%BARRIER_X, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%BARRIER_X'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%Y_l ) ) THEN
        DEALLOCATE( data%Y_l, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%Y_l'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DY_l ) ) THEN
        DEALLOCATE( data%DY_l, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DY_l'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DIST_C_l ) ) THEN
        DEALLOCATE( data%DIST_C_l, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DIST_C_l'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%Y_u ) ) THEN
        DEALLOCATE( data%Y_u, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%Y_u'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DY_u ) ) THEN
        DEALLOCATE( data%DY_u, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DY_u'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DIST_C_u ) ) THEN
        DEALLOCATE( data%DIST_C_u, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DIST_C_u'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%C ) ) THEN
        DEALLOCATE( data%C, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%C'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%BARRIER_C ) ) THEN
        DEALLOCATE( data%BARRIER_C, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%BARRIER_C'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%SCALE_C ) ) THEN
        DEALLOCATE( data%SCALE_C, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%SCALE_C'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DELTA ) ) THEN
        DEALLOCATE( data%DELTA, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DELTA'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%RHS ) ) THEN
        DEALLOCATE( data%RHS, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%RHS'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DZ_l ) ) THEN
        DEALLOCATE( data%DZ_l, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DZ_l'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DZ_u ) ) THEN
        DEALLOCATE( data%DZ_u, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DZ_u'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

!     IF ( ALLOCATED( data%C_last ) ) THEN
!       DEALLOCATE( data%C_last, STAT = inform%alloc_status )
!       IF ( inform%alloc_status /= 0 ) THEN
!         inform%status = GALAHAD_error_deallocate
!         inform%bad_alloc = 'lsqp: data%C_last'
!         IF ( control%error > 0 )                                             &
!           WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
!       END IF
!     END IF

!     IF ( ALLOCATED( data%X_last ) ) THEN
!       DEALLOCATE( data%X_last, STAT = inform%alloc_status )
!       IF ( inform%alloc_status /= 0 ) THEN
!         inform%status = GALAHAD_error_deallocate
!         inform%bad_alloc = 'lsqp: data%X_last'
!         IF ( control%error > 0 )                                             &
!           WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
!       END IF
!     END IF

      IF ( ALLOCATED( data%Y_last ) ) THEN
        DEALLOCATE( data%Y_last, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%Y_last'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%Z_last ) ) THEN
        DEALLOCATE( data%Z_last, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%Z_last'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%Abycol_val ) ) THEN
        DEALLOCATE( data%Abycol_val, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%Abycol_val'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DIAG_X ) ) THEN
        DEALLOCATE( data%DIAG_X, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DIAG_X'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DIAG_C ) ) THEN
        DEALLOCATE( data%DIAG_C, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DIAG_C'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%IW ) ) THEN
        DEALLOCATE( data%IW, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%IW'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%K_colptr ) ) THEN
        DEALLOCATE( data%K_colptr, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%K_colptr'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%Abycol_ptr ) ) THEN
        DEALLOCATE( data%Abycol_ptr, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%Abycol_ptr'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF
      IF ( ALLOCATED( data%Abycol_row ) ) THEN
        DEALLOCATE( data%Abycol_row, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%Abycol_row'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%K%row ) ) THEN
        DEALLOCATE( data%K%row, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%K%row'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%K%col ) ) THEN
        DEALLOCATE( data%K%col, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%K%col'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%K%val ) ) THEN
        DEALLOCATE( data%K%val, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%K%val'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%COEF0 ) ) THEN
        DEALLOCATE( data%COEF0, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%COEF0'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%COEF1 ) ) THEN
        DEALLOCATE( data%COEF1, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%COEF1'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%COEF2 ) ) THEN
        DEALLOCATE( data%COEF2, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%COEF2'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%COEF2 ) ) THEN
        DEALLOCATE( data%COEF2, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%COEF2'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%COEF3 ) ) THEN
        DEALLOCATE( data%COEF3, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%COEF3'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%COEF4 ) ) THEN
        DEALLOCATE( data%COEF4, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%COEF4'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DELTA_cor ) ) THEN
        DEALLOCATE( data%DELTA_cor, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DELTA_cor'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DY_cor_l ) ) THEN
        DEALLOCATE( data%DY_cor_l, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DY_cor_l'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DY_cor_u ) ) THEN
        DEALLOCATE( data%DY_cor_u, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DY_cor_u'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DZ_cor_l ) ) THEN
        DEALLOCATE( data%DZ_cor_l, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DZ_cor_l'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%DZ_cor_u ) ) THEN
        DEALLOCATE( data%DZ_cor_u, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%DZ_cor_u'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%BREAKP ) ) THEN
        DEALLOCATE( data%BREAKP, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%BREAKP'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( data%IBREAK ) ) THEN
        DEALLOCATE( data%IBREAK, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_deallocate
          inform%bad_alloc = 'lsqp: data%IBREAK'
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
        END IF
      END IF

!  Deallocate all arrays allocated within SILS

      CALL SILS_finalize( data%FACTORS, data%CNTL, inform%alloc_status )
      IF ( inform%alloc_status /= 0 ) THEN
        inform%status = GALAHAD_error_deallocate
        inform%bad_alloc = 'lsqp: data%FACTORS'
        IF ( control%error > 0 )                                               &
          WRITE( control%error, 2900 ) inform%bad_alloc, inform%alloc_status
      END IF

!  Deallocate all arrays allocated for the preprocessing stage

      CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )
      IF ( data%QPP_inform%status /= 0 ) THEN
        inform%status = GALAHAD_error_deallocate
        inform%alloc_status = data%QPP_inform%alloc_status
        inform%bad_alloc = 'lsqp: data%QPP'
      END IF

      RETURN

!  Non-executable statement

 2900 FORMAT( ' ** Message from -LSQP_terminate-', /,                          &
              ' Deallocation error for ', A, /, ' status = ', I6 ) 

!  End of subroutine LSQP_terminate

      END SUBROUTINE LSQP_terminate

!-*-  L S Q P _ L A G R A N G I A N _ G R A D I E N T   S U B R O U T I N E  -*-

      SUBROUTINE LSQP_Lagrangian_gradient( dims, n, m, X, Y, Y_l, Y_u,         &
                                           Z_l, Z_u, A_ne, A_val, A_col, A_ptr,&
                                           DIST_X_l, DIST_X_u, DIST_C_l,       &
                                           DIST_C_u, GRAD_L, getdua, dufeas,   &
                                           Hessian_kind, gradient_kind,        &
                                           control, WEIGHT, X0, G )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute the gradient of the Lagrangian function 
!
!  GRAD_L = W*W*( x - x0 ) - A(transpose) y
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, Hessian_kind, gradient_kind
      REAL ( KIND = wp ), INTENT( IN ) :: dufeas
      LOGICAL, INTENT( IN ) :: getdua 
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: GRAD_L
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_l_start : dims%x_l_end ) :: DIST_X_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_u_start : dims%x_u_end ) :: DIST_X_u
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
             DIMENSION( dims%x_free + 1 : dims%x_l_end ) :: Z_l
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
             DIMENSION( dims%x_u_start : n ) :: Z_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_l_start : dims%c_l_end ) :: DIST_C_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_u_start : dims%c_u_end ) :: DIST_C_u
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
             DIMENSION( dims%c_l_start : dims%c_l_end ) :: Y_l
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
             DIMENSION( dims%c_u_start : dims%c_u_end ) :: Y_u
      INTEGER, INTENT( IN ) :: A_ne
      INTEGER, INTENT( IN ), DIMENSION( A_ne ) :: A_col
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( A_ne ) :: A_val
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ), OPTIONAL :: WEIGHT, X0
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ), OPTIONAL :: G

!  Local variables

      INTEGER :: i
      REAL ( KIND = wp ) :: gi 

!  Add the product A( transpose ) y to the gradient of the quadratic

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        IF ( Hessian_kind == 0 ) THEN
          DO i = 1,  n ; GRAD_L( i ) = zero ; END DO
        ELSE IF ( Hessian_kind == 1 ) THEN
          DO i = 1,  n ; GRAD_L( i ) = X( i ) - X0( i ) ; END DO
        ELSE
          DO i = 1,  n
            GRAD_L( i ) = WEIGHT( i ) ** 2 * ( X( i ) - X0( i ) ) ; END DO
        END IF
        IF ( gradient_kind == 1 ) THEN
          DO i = 1,  n ; GRAD_L( i )  = GRAD_L( i )  + one ; END DO
        ELSE IF ( gradient_kind /= 0 ) THEN
          DO i = 1,  n ; GRAD_L( i )  = GRAD_L( i )  + G( i )  ; END DO
        END IF
      ELSE
        IF ( Hessian_kind == 0 ) THEN
          GRAD_L = zero
        ELSE IF ( Hessian_kind == 1 ) THEN
          GRAD_L = X - X0
        ELSE
          GRAD_L = ( WEIGHT ** 2 ) * ( X - X0 )
        END IF
        IF ( gradient_kind == 1 ) THEN
          GRAD_L = GRAD_L + one
        ELSE IF ( gradient_kind /= 0 ) THEN
          GRAD_L = GRAD_L + G
        END IF
      END IF      

      CALL LSQP_AX( n, GRAD_L, m, A_ne, A_val, A_col, A_ptr, m, Y, '-T' )

!  If required, obtain suitable "good" starting values for the dual
!  variables ( see paper )

      IF ( getdua ) THEN 

!  Problem variables:

!  The variable is a non-negativity

        DO i = dims%x_free + 1, dims%x_l_start - 1
          Z_l( i ) = MAX( dufeas, GRAD_L( i ) / ( one + X( i ) ** 2 ) ) 
        END DO

!  The variable has just a lower bound

        DO i = dims%x_l_start, dims%x_u_start - 1
          Z_l( i ) = MAX( dufeas, GRAD_L( i ) / ( one + DIST_X_l( i ) ** 2 ) ) 
        END DO

!  The variable has both lower and upper bounds

        DO i = dims%x_u_start, dims%x_l_end
          gi = GRAD_L( i ) 
          IF ( ABS( gi ) <= dufeas ) THEN 
            Z_l( i ) = dufeas ; Z_u( i ) = - dufeas 
          ELSE IF ( gi > dufeas ) THEN 
            Z_l( i ) = ( gi + dufeas ) / ( one + DIST_X_l( i ) ** 2 )
            Z_u( i ) = - dufeas 
          ELSE 
            Z_l( i ) = dufeas 
            Z_u( i ) = ( gi - dufeas ) / ( one + DIST_X_u( i ) ** 2 )
          END IF 
        END DO

!  The variable has just an upper bound

        DO i = dims%x_l_end + 1, dims%x_u_end
          Z_u( i ) = MIN( - dufeas, GRAD_L( i ) / ( one + DIST_X_u( i ) ** 2 ) )
        END DO

!  The variable is a non-positivity

        DO i = dims%x_u_end + 1, n
          Z_u( i ) = MIN( - dufeas, GRAD_L( i ) / ( one + X( i ) ** 2 ) )
        END DO

!  Slack variables:

!  The variable has just a lower bound

        DO i = dims%c_l_start, dims%c_u_start - 1
          Y_l( i ) = MAX( dufeas, - Y( i ) / ( one + DIST_C_l( i ) ** 2 ) ) 
        END DO

!  The variable has both lower and upper bounds

        DO i = dims%c_u_start, dims%c_l_end
          gi = - Y( i ) 
          IF ( ABS( gi ) <= dufeas ) THEN 
            Y_l( i ) = dufeas ; Y_u( i ) = - dufeas 
          ELSE IF ( gi > dufeas ) THEN 
            Y_l( i ) = ( gi + dufeas ) / ( one + DIST_C_l( i ) ** 2 )
            Y_u( i ) = - dufeas 
          ELSE 
            Y_l( i ) = dufeas 
            Y_u( i ) = ( gi - dufeas ) / ( one + DIST_C_u( i ) ** 2 )
          END IF 
        END DO

!  The variable has just an upper bound

        DO i = dims%c_l_end + 1, dims%c_u_end
          Y_u( i ) = MIN( - dufeas, - Y( i ) / ( one + DIST_C_u( i ) ** 2 ) )
        END DO
      END IF 

      RETURN  

!  End of LSQP_Lagrangian_gradient

      END SUBROUTINE LSQP_Lagrangian_gradient

!-*-*-*-  L S Q P _ P O T E N T I A L _ V A L U E   S U B R O U T I N E  -*-*-*-

      FUNCTION LSQP_potential_value( dims, n, X, DIST_X_l, DIST_X_u,          &
                                     DIST_C_l, DIST_C_u )
      REAL ( KIND = wp ) LSQP_potential_value

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute the value of the potential function
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_l_start : dims%x_l_end ) :: DIST_X_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_u_start : dims%x_u_end ) :: DIST_X_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_l_start : dims%c_l_end ) :: DIST_C_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_u_start : dims%c_u_end ) :: DIST_C_u

! Compute the potential terms

      LSQP_potential_value =                                                   &
        - SUM( LOG( X( dims%x_free + 1 : dims%x_l_start - 1 ) ) )              &
        - SUM( LOG( DIST_X_l ) ) - SUM( LOG( DIST_X_u ) )                      &
        - SUM( LOG( - X( dims%x_u_end + 1 : n ) ) )                            &
        - SUM( LOG( DIST_C_l ) ) - SUM( LOG( DIST_C_u ) )

      RETURN  

!  End of LSQP_potential_value

      END FUNCTION LSQP_potential_value
 
!-*-*-*-*-*-   L S Q P _ M E R I T _ V A L U E   F U N C T I O N   -*-*-*-*-*-*-

      FUNCTION LSQP_merit_value( dims, n, m, X, Y, Y_l, Y_u, Z_l, Z_u,         &
                                 DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u,       &
                                 GRAD_L, C_RES, res_dual )
      REAL ( KIND = wp ) LSQP_merit_value

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute the value of the merit function
!
!     | < z_l . ( x - x_l ) > +  < z_u . ( x_u - x ) > +
!       < y_l . ( c - c_l ) > +  < y_u . ( c_u - c ) > | + 
!       || ( GRAD_L - z_l - z_u ) ||
!       || (   y - y_l - y_u    ) ||
!       || (  A x - SCALE_c * c ) || 
!
!  where GRAD_L = W*W*( x - x0 ) - A(transpose) y is the gradient 
!  of the Lagrangian
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m
      REAL ( KIND = wp ), INTENT( OUT ) :: res_dual
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X, GRAD_L
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_l_start : dims%x_l_end ) :: DIST_x_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_free + 1 : dims%x_l_end ) :: Z_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_u_start : dims%x_u_end ) :: DIST_X_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_u_start : n ) :: Z_u
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y, C_RES 
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_l_start : dims%c_l_end ) :: Y_l, DIST_C_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_u_start : dims%c_u_end ) :: Y_u, DIST_C_u

!  Local variables

      INTEGER :: i
      REAL ( KIND = wp ) :: res_cs

!  Compute in the l_2-norm

      res_dual = SUM( GRAD_L( : dims%x_free ) ** 2 ) ; res_cs = zero

!  Problem variables:

      DO i = dims%x_free + 1, dims%x_l_start - 1
        res_dual = res_dual + ( GRAD_L( i ) - Z_l( i ) ) ** 2 
        res_cs = res_cs + Z_l( i ) * X( i )
      END DO
      DO i = dims%x_l_start, dims%x_u_start - 1
        res_dual = res_dual + ( GRAD_L( i ) - Z_l( i ) ) ** 2 
        res_cs = res_cs + Z_l( i ) * DIST_X_l( i )
      END DO
      DO i = dims%x_u_start, dims%x_l_end
        res_dual = res_dual + ( GRAD_L( i ) - Z_l( i ) - Z_u( i ) ) ** 2 
        res_cs = res_cs + Z_l( i ) * DIST_X_l( i ) - Z_u( i ) * DIST_X_u( i ) 
      END DO
      DO i = dims%x_l_end + 1, dims%x_u_end
        res_dual = res_dual + ( GRAD_L( i ) - Z_u( i ) ) ** 2 
        res_cs = res_cs - Z_u( i ) * DIST_X_u( i )
      END DO
      DO i = dims%x_u_end + 1, n
        res_dual = res_dual + ( GRAD_L( i ) - Z_u( i ) ) ** 2 
        res_cs = res_cs + Z_u( i ) * X( i )
      END DO

!  Slack variables:

      DO i = dims%c_l_start, dims%c_u_start - 1
        res_dual = res_dual + ( Y( i ) - Y_l( i ) ) ** 2 
        res_cs = res_cs + Y_l( i ) * DIST_C_l( i )
      END DO
      DO i = dims%c_u_start, dims%c_l_end
        res_dual = res_dual + ( Y( i ) - Y_l( i ) - Y_u( i ) ) ** 2 
        res_cs = res_cs + Y_l( i ) * DIST_C_l( i ) - Y_u( i ) * DIST_C_u( i )
      END DO
      DO i = dims%c_l_end + 1, dims%c_u_end
        res_dual = res_dual + ( Y( i ) - Y_u( i ) ) ** 2 
        res_cs = res_cs - Y_u( i ) * DIST_C_u( i )
      END DO

      LSQP_merit_value = ABS( res_cs ) + SQRT( res_dual + SUM( C_RES ** 2 ) )
      res_dual = SQRT( res_dual )

      RETURN  

!  End of function LSQP_merit_value

      END FUNCTION LSQP_merit_value

!-*-  L S Q P _ I T E R A T I V E _ R E F I N E M E N T  S U B R O U T I N E -*-

      SUBROUTINE LSQP_iterative_refinement(                                    &
                 dims, n, m, K, FACTORS, CNTL, A_ne, A_val, A_col, A_ptr,      &
                 BARRIER_X, BARRIER_C, RHS, factor, DIAG_X, ldiag_x,           &
                 SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,                        &
                 SOL, RES, BEST, SOL_y, BEST_y, RES_y, RES_x,                  &
                 itref_max, print_level, control, inform, Hessian_kind, WEIGHT )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Solve the block linear system

!     ( DIAG_X           A(trans)  )
!     (         DIAG_C  - SCALE_C  ) ( sol ) = ( rhs )
!     (   A   - SCALE_C            )

!  using iterative refinement, returning sol in rhs

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, A_ne, ldiag_x, ldiag_c_l, ldiag_c_u
      INTEGER, INTENT( IN ) :: itref_max, factor, Hessian_kind, print_level
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_free + 1 : n ) :: BARRIER_X
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_l_start : dims%c_u_end ) :: BARRIER_C, SCALE_C
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( ldiag_x ) :: DIAG_X
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( ldiag_c_l : ldiag_c_u ) :: DIAG_C
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( dims%v_e ) :: RHS 
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
                          DIMENSION( dims%v_e ) :: SOL, RES
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( dims%v_e ) :: BEST
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
                          DIMENSION( m ) :: SOL_y, BEST_y, RES_y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: RES_x
      INTEGER, INTENT( IN ), DIMENSION( A_ne ) :: A_col
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( A_ne ) :: A_val
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ), OPTIONAL :: WEIGHT
      TYPE ( SMT_type ), INTENT( INOUT ) :: K
      TYPE ( SILS_factors ), INTENT( INOUT ) :: FACTORS
      TYPE ( SILS_control ), INTENT( INOUT ) :: CNTL
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      TYPE ( LSQP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: i, it, itref_max_block
      REAL ( KIND = wp ) :: res_norm, old_res
      TYPE ( SILS_sinfo ) :: SINFO

      IF ( control%out > 0 .AND. print_level >= 3 ) THEN
        IF ( dims%nc > 0 ) THEN
          WRITE( control%out, 2000 )                                           &
             MAXVAL( ABS( RHS( dims%x_s : dims%x_e ) ) ),                      &
             MAXVAL( ABS( RHS( dims%c_s : dims%c_e ) ) ),                      &
             MAXVAL( ABS( RHS( dims%y_s : dims%y_e ) ) )
        ELSE IF ( m > 0 ) THEN
          WRITE( control%out, 2010 )                                           &
             MAXVAL( ABS( RHS( dims%x_s : dims%x_e ) ) ),                      &
             MAXVAL( ABS( RHS( dims%y_s : dims%y_e ) ) )
        ELSE
          WRITE( control%out, 2020 ) MAXVAL( ABS( RHS( dims%x_s : dims%x_e ) ) )
        END IF
      END IF

      old_res = MAXVAL( ABS( RHS ) )
      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, dims%v_e ; SOL( i ) = RHS( i ) ; END DO
      ELSE
        SOL = RHS
      END IF

      IF ( factor == 0 .OR. factor == 1 ) THEN
        itref_max_block = itref_max - 1
        CALL LSQP_block_solve( dims, n, m, SOL( dims%x_s : dims%x_e ),         &
                               SOL( dims%c_s : dims%c_e ),                     &
                               SOL( dims%y_s : dims%y_e ),                     &
                               A_ne, A_val, A_col, A_ptr, DIAG_X, ldiag_x,     &
                               SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,          &
                               SOL_y, BEST_y, RES_y, RES_x, itref_max_block,   &
                               K, FACTORS, CNTL, print_level, control, inform )
      ELSE
        CALL SILS_solve( K, FACTORS, SOL, CNTL, SINFO )
        inform%factorization_status = SINFO%flag
        IF ( SINFO%flag /= 0 ) THEN
          IF ( control%error > 0 .AND. print_level >= 1 )                      &
              WRITE( control%error, 2100 ) SINFO%flag
          inform%status = GALAHAD_error_solve ; RETURN
        END IF
      END IF

      IF ( control%out > 0 .AND. print_level >= 3 ) THEN
        IF ( dims%nc > 0 ) THEN
          WRITE( control%out, 2030 )                                           &
            MAXVAL( ABS( SOL( dims%x_s : dims%x_e ) ) ),                       &
            MAXVAL( ABS( SOL( dims%c_s : dims%c_e ) ) ),                       &
            MAXVAL( ABS( SOL( dims%y_s : dims%y_e ) ) )
        ELSE IF ( m > 0 ) THEN
          WRITE( control%out, 2040 )                                           &
            MAXVAL( ABS( SOL( dims%x_s : dims%x_e ) ) ),                       &
            MAXVAL( ABS( SOL( dims%y_s : dims%y_e ) ) )
        ELSE
          WRITE( control%out, 2050 ) MAXVAL( ABS( SOL( dims%x_s : dims%x_e ) ) )
        END IF
      END IF

      DO it = 1, itref_max

!  Obtain the residual:
!  Remember the barrier terms, and any diagonal perturbations

        IF ( control%array_syntax_worse_than_do_loop ) THEN
          IF ( Hessian_kind == 0 ) THEN
            DO i = dims%x_s, dims%x_free ; RES( i ) = RHS( i ) ; END DO
            DO i = dims%x_free + 1, dims%x_e 
              RES( i ) = RHS( i ) - BARRIER_X( i ) * SOL( i ) ; END DO
          ELSE IF ( Hessian_kind == 1 ) THEN
            DO i = dims%x_s, dims%x_free
              RES( i ) = RHS( i ) - SOL( i ) ; END DO
            DO i = dims%x_free + 1, dims%x_e 
              RES( i ) = RHS( i ) - ( one + BARRIER_X( i ) ) * SOL( i ) ; END DO
          ELSE
            DO i = dims%x_s, dims%x_free
              RES( i ) = RHS( i ) - ( WEIGHT( i ) ** 2 ) * SOL( i ) ; END DO
            DO i = dims%x_free + 1, dims%x_e 
              RES( i ) = RHS( i ) - ( WEIGHT( i ) ** 2 + BARRIER_X( i ) ) *    &
               SOL( i ) ; END DO
          END IF
          DO i = dims%y_s, dims%y_e ; RES( i ) = RHS( i ) ; END DO
          DO i = 0, dims%nc - 1
            RES( dims%c_s + i ) = RHS( dims%c_s + i ) -                        &
              BARRIER_C( dims%c_l_start + i ) * SOL( dims%c_s + i ) ; END DO
        ELSE
          IF ( Hessian_kind == 0 ) THEN
            RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free )
            RES( dims%x_free + 1 : dims%x_e ) =                                &
              RHS( dims%x_free + 1 : dims%x_e ) -                              &
                 BARRIER_X * SOL( dims%x_free + 1 : dims%x_e )
          ELSE IF ( Hessian_kind == 1 ) THEN
            RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free ) -    &
                                            SOL( : dims%x_free )
            RES( dims%x_free + 1 : dims%x_e ) =                                &
              RHS( dims%x_free + 1 : dims%x_e ) -                              &
                ( one + BARRIER_X ) * SOL( dims%x_free + 1 : dims%x_e )
          ELSE
            RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free ) -    &
                              ( WEIGHT( dims%x_s : dims%x_free ) ** 2 ) *      &
                                 SOL( : dims%x_free )
            RES( dims%x_free + 1 : dims%x_e ) =                                &
              RHS( dims%x_free + 1 : dims%x_e ) -                              &
               ( WEIGHT( dims%x_free + 1 : dims%x_e ) ** 2 +                   &
                 BARRIER_X ) * SOL( dims%x_free + 1 : dims%x_e )
          END IF
          RES( dims%y_s : dims%y_e ) = RHS( dims%y_s : dims%y_e )
          RES( dims%c_s : dims%c_e ) = RHS( dims%c_s : dims%c_e ) -            &
            BARRIER_C * SOL( dims%c_s : dims%c_e )
        END IF

!  Include the contribution from A

        CALL LSQP_AX( n, RES( dims%x_s : dims%x_e ), m, A_ne, A_val, A_col,    &
                      A_ptr, m, SOL( dims%y_s : dims%y_e ), '-T' )
        CALL LSQP_AX( m, RES( dims%y_s : dims%y_e ), m, A_ne, A_val, A_col,    &
                      A_ptr, n, SOL( dims%x_s : dims%x_e ), '- ' )


!  Include the contribution from the slack variables

        IF ( control%array_syntax_worse_than_do_loop ) THEN
          DO i = 0, dims%nc - 1
            RES( dims%c_s + i ) = RES( dims%c_s + i ) +                        &
                SCALE_C( dims%c_l_start + i ) * SOL( dims%y_i + i )
            RES( dims%y_i + i ) = RES( dims%y_i + i ) +                        &
                SCALE_C( dims%c_l_start + i ) * SOL( dims%c_s + i )  
          END DO
        ELSE
          RES( dims%c_s : dims%c_e ) =                                         &
            RES( dims%c_s : dims%c_e ) + SCALE_C * SOL( dims%y_i : dims%y_e )
          RES( dims%y_i : dims%y_e ) =                                         &
            RES( dims%y_i : dims%y_e ) + SCALE_C * SOL( dims%c_s : dims%c_e )
        END IF

        IF ( control%out > 0 .AND. print_level >= 3 ) THEN
          IF ( dims%nc > 0 ) THEN
            WRITE( control%out, 2000 )                                         &
              MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) ),                     &
              MAXVAL( ABS( RES( dims%c_s : dims%c_e ) ) ),                     &
              MAXVAL( ABS( RES( dims%y_s : dims%y_e ) ) )
          ELSE IF ( m > 0 ) THEN
            WRITE( control%out, 2010 )                                         &
              MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) ),                     &
              MAXVAL( ABS( RES( dims%y_s : dims%y_e ) ) )
          ELSE
            WRITE( control%out, 2020 )                                         &
              MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) )
          END IF
        END IF

        res_norm = MAXVAL( ABS( RES ) )
        IF ( it > 1 ) THEN
          IF ( res_norm < old_res ) THEN
            old_res = res_norm
            IF ( control%array_syntax_worse_than_do_loop ) THEN
              DO i = 1, dims%v_e ; BEST( i ) = SOL( i ) ; END DO
            ELSE
              BEST = SOL
            END IF
          ELSE
            IF ( control%array_syntax_worse_than_do_loop ) THEN
              DO i = 1, dims%v_e ; SOL( i ) = BEST( i ) ; END DO
            ELSE
              SOL = BEST
            END IF
            EXIT
          END IF
        ELSE
          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, dims%v_e ; BEST( i ) = SOL( i ) ; END DO
          ELSE
            BEST = SOL
          END IF
        END IF
   
        IF ( factor == 0 .OR. factor == 1 ) THEN
          CALL LSQP_block_solve( dims, n, m, RES( dims%x_s : dims%x_e ),       &
                                 RES( dims%c_s : dims%c_e ),                   &
                                 RES( dims%y_s : dims%y_e ),                   &
                                 A_ne, A_val, A_col, A_ptr,                    &
                                 DIAG_X, ldiag_x, SCALE_C, DIAG_C, ldiag_c_l,  &
                                 ldiag_c_u, SOL_y,                             &
                                 BEST_y, RES_y, RES_x,                         &
                                 itref_max_block, K, FACTORS, CNTL,            &
                                 print_level, control, inform )

        ELSE
          CALL SILS_solve( K, FACTORS, RES, CNTL, SINFO )
          inform%factorization_status = SINFO%flag
          IF ( SINFO%flag /= 0 ) THEN
            IF ( control%error > 0 .AND. print_level >= 1 )                    &
               WRITE( control%error, 2100 ) SINFO%flag
            inform%status = GALAHAD_error_solve ; RETURN
          END IF
        END IF

        IF ( control%out > 0 .AND. print_level >= 3 ) THEN
          IF ( dims%nc > 0 ) THEN
            WRITE( control%out, 2060 )                                         &
              MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) ),                     &
              MAXVAL( ABS( RES( dims%c_s : dims%c_e ) ) ),                     &
              MAXVAL( ABS( RES( dims%y_s : dims%y_e ) ) )
          ELSE IF ( m > 0 ) THEN
            WRITE( control%out, 2070 )                                         &
              MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) ),                     &
              MAXVAL( ABS( RES( dims%y_s : dims%y_e ) ) )
          ELSE
            WRITE( control%out, 2080 )                                         &
              MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) )
          END IF
        END IF

        IF ( control%array_syntax_worse_than_do_loop ) THEN
          DO i = 1, dims%v_e ; SOL( i ) = SOL( i ) + RES( i ) ; END DO
        ELSE
          SOL = SOL + RES
        END IF

!  Print final residuals if required

        IF ( control%out > 0 .AND. print_level >= 3                            &
             .AND. it == itref_max ) THEN
   
!  Remember the barrier terms, and any diagonal perturbations

          IF ( control%array_syntax_worse_than_do_loop ) THEN
            IF ( Hessian_kind == 0 ) THEN
              DO i = dims%x_s, dims%x_free ; RES( i ) = RHS( i ) ; END DO
              DO i = dims%x_free + 1, dims%x_e 
                RES( i ) = RHS( i ) - BARRIER_X( i ) * SOL( i ) ; END DO
            ELSE IF ( Hessian_kind == 1 ) THEN
              DO i = dims%x_s, dims%x_free
                RES( i ) = RHS( i ) - SOL( i ) ; END DO
              DO i = dims%x_free + 1, dims%x_e 
                RES( i ) = RHS( i ) - ( one + BARRIER_X( i ) ) * SOL( i ) 
              END DO
            ELSE
              DO i = dims%x_s, dims%x_free
                RES( i ) = RHS( i ) - ( WEIGHT( i ) ** 2 ) * SOL( i ) ; END DO
              DO i = dims%x_free + 1, dims%x_e 
                RES( i ) = RHS( i ) - ( WEIGHT( i ) ** 2 + BARRIER_X( i ) ) *  &
                 SOL( i ) ; END DO
            END IF
            DO i = dims%y_s, dims%y_e ; RES( i ) = RHS( i ) ; END DO
            DO i = 0, dims%nc - 1
              RES( dims%c_s + i ) = RHS( dims%c_s + i ) -                      &
                BARRIER_C( dims%c_l_start + i ) * SOL( dims%c_s + i ) ; END DO
          ELSE
            IF ( Hessian_kind == 0 ) THEN
              RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free )
              RES( dims%x_free + 1 : dims%x_e ) =                              &
                RHS( dims%x_free + 1 : dims%x_e ) -                            &
                   BARRIER_X * SOL( dims%x_free + 1 : dims%x_e )
            ELSE IF ( Hessian_kind == 1 ) THEN
              RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free ) -  &
                                              SOL( : dims%x_free )
              RES( dims%x_free + 1 : dims%x_e ) =                              &
                RHS( dims%x_free + 1 : dims%x_e ) -                            &
                  ( one + BARRIER_X ) * SOL( dims%x_free + 1 : dims%x_e )
            ELSE
              RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free ) -  &
                                ( WEIGHT( dims%x_s : dims%x_free ) ** 2 ) *    &
                                   SOL( : dims%x_free )
              RES( dims%x_free + 1 : dims%x_e ) =                              &
                RHS( dims%x_free + 1 : dims%x_e ) -                            &
                 ( WEIGHT( dims%x_free + 1 : dims%x_e ) ** 2 +                 &
                   BARRIER_X ) * SOL( dims%x_free + 1 : dims%x_e )
            END IF
            RES( dims%y_s : dims%y_e ) = RHS( dims%y_s : dims%y_e )
            RES( dims%c_s : dims%c_e ) = RHS( dims%c_s : dims%c_e ) -          &
              BARRIER_C * SOL( dims%c_s : dims%c_e )
          END IF

!  Include the contribution from A

          CALL LSQP_AX( n, RES( dims%x_s : dims%x_e ), m, A_ne, A_val, A_col,  &
                        A_ptr, m, SOL( dims%y_s : dims%y_e ), '-T' )
          CALL LSQP_AX( m, RES( dims%y_s : dims%y_e ), m, A_ne, A_val, A_col,  &
                        A_ptr, n, SOL( dims%x_s : dims%x_e ), '- ' )

!  Include the contribution from the slack variables

          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 0, dims%nc - 1
              RES( dims%c_s + i ) = RES( dims%c_s + i ) +                      &
                  SCALE_C( dims%c_l_start + i ) * SOL( dims%y_i + i )
              RES( dims%y_i + i ) = RES( dims%y_i + i ) +                      &
                  SCALE_C( dims%c_l_start + i ) * SOL( dims%c_s + i )  
            END DO
          ELSE
            RES( dims%c_s : dims%c_e ) =                                       &
              RES( dims%c_s : dims%c_e ) + SCALE_C * SOL( dims%y_i : dims%y_e )
            RES( dims%y_i : dims%y_e ) =                                       &
              RES( dims%y_i : dims%y_e ) + SCALE_C * SOL( dims%c_s : dims%c_e )
          END IF

          IF ( dims%nc > 0 ) THEN
            WRITE( control%out, 2000 )                                         &
              MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) ),                     &
              MAXVAL( ABS( RES( dims%c_s : dims%c_e ) ) ),                     &
              MAXVAL( ABS( RES( dims%y_s : dims%y_e ) ) )
          ELSE IF ( m > 0 ) THEN
            WRITE( control%out, 2010 )                                         &
              MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) ),                     &
              MAXVAL( ABS( RES( dims%y_s : dims%y_e ) ) )
          ELSE
            WRITE( control%out, 2020 )                                         &
              MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) )
          END IF

        END IF

      END DO

!  Print final residuals if required

      IF ( control%out > 0 .AND. print_level >= 3                              &
           .AND. itref_max == 0 ) THEN
   
!  Remember the barrier terms, and any diagonal perturbations

        IF ( control%array_syntax_worse_than_do_loop ) THEN
          IF ( Hessian_kind == 0 ) THEN
            DO i = dims%x_s, dims%x_free ; RES( i ) = RHS( i ) ; END DO
            DO i = dims%x_free + 1, dims%x_e 
              RES( i ) = RHS( i ) - BARRIER_X( i ) * SOL( i ) ; END DO
          ELSE IF ( Hessian_kind == 1 ) THEN
            DO i = dims%x_s, dims%x_free
              RES( i ) = RHS( i ) - SOL( i ) ; END DO
            DO i = dims%x_free + 1, dims%x_e 
              RES( i ) = RHS( i ) - ( one + BARRIER_X( i ) ) * SOL( i ) ; END DO
          ELSE
            DO i = dims%x_s, dims%x_free
              RES( i ) = RHS( i ) - ( WEIGHT( i ) ** 2 ) * SOL( i ) ; END DO
            DO i = dims%x_free + 1, dims%x_e 
              RES( i ) = RHS( i ) - ( WEIGHT( i ) ** 2 + BARRIER_X( i ) ) *    &
               SOL( i ) ; END DO
          END IF
          DO i = dims%y_s, dims%y_e ; RES( i ) = RHS( i ) ; END DO
          DO i = 0, dims%nc - 1
            RES( dims%c_s + i ) = RHS( dims%c_s + i ) -                        &
              BARRIER_C( dims%c_l_start + i ) * SOL( dims%c_s + i ) ; END DO
        ELSE
          IF ( Hessian_kind == 0 ) THEN
            RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free )
            RES( dims%x_free + 1 : dims%x_e ) =                                &
              RHS( dims%x_free + 1 : dims%x_e ) -                              &
                 BARRIER_X * SOL( dims%x_free + 1 : dims%x_e )
          ELSE IF ( Hessian_kind == 1 ) THEN
            RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free ) -    &
                                            SOL( : dims%x_free )
            RES( dims%x_free + 1 : dims%x_e ) =                                &
              RHS( dims%x_free + 1 : dims%x_e ) -                              &
                ( one + BARRIER_X ) * SOL( dims%x_free + 1 : dims%x_e )
          ELSE
            RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free ) -    &
                              ( WEIGHT( dims%x_s : dims%x_free ) ** 2 ) *      &
                                 SOL( : dims%x_free )
            RES( dims%x_free + 1 : dims%x_e ) =                                &
              RHS( dims%x_free + 1 : dims%x_e ) -                              &
               ( WEIGHT( dims%x_free + 1 : dims%x_e ) ** 2 +                   &
                 BARRIER_X ) * SOL( dims%x_free + 1 : dims%x_e )
          END IF
          RES( dims%y_s : dims%y_e ) = RHS( dims%y_s : dims%y_e )
          RES( dims%c_s : dims%c_e ) = RHS( dims%c_s : dims%c_e ) -            &
            BARRIER_C * SOL( dims%c_s : dims%c_e )
        END IF

!  Include the contribution from A

        CALL LSQP_AX( n, RES( dims%x_s : dims%x_e ), m, A_ne, A_val, A_col,    &
                      A_ptr, m, SOL( dims%y_s : dims%y_e ), '-T' )
        CALL LSQP_AX( m, RES( dims%y_s : dims%y_e ), m, A_ne, A_val, A_col,    &
                              A_ptr, n, SOL( dims%x_s : dims%x_e ), '- ' )

!  Include the contribution from the slack variables

        IF ( control%array_syntax_worse_than_do_loop ) THEN
          DO i = 0, dims%nc - 1
            RES( dims%c_s + i ) = RES( dims%c_s + i ) +                        &
                SCALE_C( dims%c_l_start + i ) * SOL( dims%y_i + i )
            RES( dims%y_i + i ) = RES( dims%y_i + i ) +                        &
                SCALE_C( dims%c_l_start + i ) * SOL( dims%c_s + i )  
          END DO
        ELSE
          RES( dims%c_s : dims%c_e ) =                                         &
            RES( dims%c_s : dims%c_e ) + SCALE_C * SOL( dims%y_i : dims%y_e )
          RES( dims%y_i : dims%y_e ) =                                         &
            RES( dims%y_i : dims%y_e ) + SCALE_C * SOL( dims%c_s : dims%c_e )
        END IF

        IF ( dims%nc > 0 ) THEN
          WRITE( control%out, 2000 )                                           &
            MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) ),                       &
            MAXVAL( ABS( RES( dims%c_s : dims%c_e ) ) ),                       &
            MAXVAL( ABS( RES( dims%y_s : dims%y_e ) ) )
        ELSE IF ( m > 0 ) THEN
          WRITE( control%out, 2010 )                                           &
            MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) ),                       &
            MAXVAL( ABS( RES( dims%y_s : dims%y_e ) ) )
        ELSE
          WRITE( control%out, 2020 ) MAXVAL( ABS( RES( dims%x_s : dims%x_e ) ) )
        END IF

      END IF

      IF ( control%out > 0 .AND. print_level >= 3 ) THEN
        IF ( dims%nc > 0 ) THEN
          WRITE( control%out, 2030 )                                           &
            MAXVAL( ABS( SOL( dims%x_s : dims%x_e ) ) ),                       &
            MAXVAL( ABS( SOL( dims%c_s : dims%c_e ) ) ),                       &
            MAXVAL( ABS( SOL( dims%y_s : dims%y_e ) ) )
        ELSE IF ( m > 0 ) THEN
          WRITE( control%out, 2040 )                                           &
            MAXVAL( ABS( SOL( dims%x_s : dims%x_e ) ) ),                       &
            MAXVAL( ABS( SOL( dims%y_s : dims%y_e ) ) )
        ELSE
          WRITE( control%out, 2050 ) MAXVAL( ABS( SOL( dims%x_s : dims%x_e ) ) )
        END IF
      END IF

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, dims%v_e
          RHS( i ) = SOL( i )
        END DO
      ELSE
        RHS = SOL
      END IF 
      inform%status = GALAHAD_ok

      RETURN

!  Non-executable statements

 2000 FORMAT( '  res(dual) ', ES12.4, '  res(slack ) ', ES12.4,                &
              '  res(primal) ', ES12.4 )
 2010 FORMAT( '  res(dual) ', ES12.4, '  res(primal) ', ES12.4 )
 2020 FORMAT( '  res(dual) ', ES12.4 )
 2030 FORMAT( '  sol(x   ) ', ES12.4, '  sol(slack ) ', ES12.4,                &
              '  sol(y     ) ', ES12.4 )
 2040 FORMAT( '  sol(x   ) ', ES12.4, '  sol(y     ) ', ES12.4 )
 2050 FORMAT( '  sol(x   ) ', ES12.4 )
 2060 FORMAT( ' dsol(x   ) ', ES12.4, ' dsol(slack ) ', ES12.4,                &
              ' dsol(y     ) ', ES12.4 )
 2070 FORMAT( ' dsol(x   ) ', ES12.4, ' dsol(y     ) ', ES12.4 )
 2080 FORMAT( ' dsol(x   ) ', ES12.4 )
 2100 FORMAT( '   **  Error return ', I3, ' from SILS_solve ' ) 

!  End of LSQP_iterative_refinement

      END SUBROUTINE LSQP_iterative_refinement

!-*-*-*-*-*-   L S Q P _ B L O C K _ S O L V E   S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE LSQP_block_solve( dims, n, m, RHS_x, RHS_c, RHS_y, A_ne,      &
                                   A_val, A_col, A_ptr, DIAG_X, ldiag_x,       &
                                   SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,      &
                                   SOL_y, BEST_y, RES_y, RES_x,                &
                                   itref_max, K, FACTORS, CNTL,                &
                                   print_level, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Solve the block system 

!     ( DIAG_X           A(trans)  ) (sol_x) = (rhs_x)
!     (         DIAG_C  - SCALE_C  ) (sol_c)   (rhs_c)
!     (   A   - SCALE_C            ) (sol_y)   (rhs_y)

!  returning ( sol_x , sol_c, sol_y ) in ( rhs_x , rhs_c, rhs_y )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!   Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, A_ne, ldiag_x, ldiag_c_l, ldiag_c_u
      INTEGER, INTENT( IN ) :: itref_max, print_level
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( dims%c_l_start : dims%c_u_end ) :: SCALE_C
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( ldiag_x ) :: DIAG_X
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( ldiag_c_l : ldiag_c_u ) :: DIAG_C
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: RHS_x
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: RHS_y
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
                          DIMENSION( dims%c_l_start : m ) :: RHS_c
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
                          DIMENSION( m ) :: SOL_y, BEST_y, RES_y
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: RES_x
      INTEGER, INTENT( IN ), DIMENSION( A_ne ) :: A_col
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( A_ne ) :: A_val
      TYPE ( SMT_type ), INTENT( INOUT ) :: K
      TYPE ( SILS_factors ), INTENT( INOUT ) :: FACTORS
      TYPE ( SILS_control ), INTENT( INOUT ) :: CNTL
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      TYPE ( LSQP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: i

!  Obtain DIAG_X(inv) * rhs_x ...

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, n ; RHS_x( i ) = RHS_x( i ) / DIAG_X( i ) ; END DO
      ELSE
        RHS_x = RHS_x / DIAG_X
      END IF
      IF (  m == 0 ) RETURN

!  ... and DIAG_C(inv) * rhs_c

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = dims%c_l_start, m
          RHS_c( i ) = RHS_c( i ) / DIAG_C( i ) ; END DO
      ELSE
        RHS_c = RHS_c / DIAG_C
      END IF

!  Now find A * DIAG_X(inv) * rhs_x - rhs_y - SCALE_C * DIAG_C(inv) * rhs_c

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, dims%c_equality ; RHS_y( i ) = - RHS_y( i ) ; END DO
        DO i = dims%c_l_start,  dims%c_u_end
          RHS_y( i ) = - RHS_y( i ) - SCALE_C( i ) * RHS_c( i ) ; END DO
      ELSE
        RHS_y( : dims%c_equality ) = - RHS_y( : dims%c_equality )
        RHS_y( dims%c_l_start : ) = - RHS_y( dims%c_l_start : ) - SCALE_C *RHS_c
      END IF
      CALL LSQP_AX( m, RHS_y, m, A_ne, A_val, A_col, A_ptr,          &
                    n, RHS_x, '+ ' )

!  Next solve K sol_y = A * DIAG_X(inv) * rhs_x - rhs_y - DIAG_X(inv) * rhs_c, 
!  where K = A * DIAG_X(inv) * A(trans) + SCALE_C * DIAG_C(inv) * SCALE_C

      CALL LSQP_block_iterative_refinement(                                    &
                dims, n, m, K, FACTORS, CNTL, RHS_y, A_ne, A_val, A_col,       &
                A_ptr, DIAG_X, ldiag_x, SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u, &
                SOL_y, BEST_y, RES_y, RES_x, itref_max, print_level, control,  &
                inform )

!  Finally, recover rhs_x = DIAG_X(inv) ( rhs_x - A(trans) sol_y ) ...

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, n ; RHS_x( i ) = RHS_x( i ) * DIAG_X( i ) ; END DO
      ELSE
        RHS_x = RHS_x * DIAG_X
      END IF

      CALL LSQP_AX( n, RHS_x, m, A_ne, A_val, A_col, A_ptr,          &
                    m, RHS_y, '-T' )

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, n ; RHS_x( i ) = RHS_x( i ) / DIAG_X( i ) ; END DO
      ELSE
        RHS_x = RHS_x / DIAG_X
      END IF

!  ... and rhs_c = DIAG_C(inv) ( rhs_c + SCALE_C * sol_y )

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = dims%c_l_start,  dims%c_u_end
          RHS_c( i ) = RHS_c( i ) + ( SCALE_C( i ) / DIAG_C( i ) ) * RHS_y( i )
        END DO
      ELSE
        RHS_c = RHS_c + ( SCALE_C / DIAG_C ) * RHS_y( dims%c_l_start : )
      END IF

      RETURN

!  End of subroutine LSQP_block_solve

      END SUBROUTINE LSQP_block_solve

!-*-*-*-*-*-*-   L S Q P _ B L O C K _ I-R   S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE LSQP_block_iterative_refinement( dims, n, m,                  &
                            K, FACTORS, CNTL, RHS_y,                           &
                            A_ne, A_val, A_col, A_ptr, DIAG_X, ldiag_x,        &
                            SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,             &
                            SOL_y, BEST_y, RES_y, RES_x,                       &
                            itref_max, print_level, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Solve the system 

!   (  A * DIAG_X(inv) * A(trans) + SCALE_C * DIAG_C(inv) * SCALE_C ) *
!     sol_y = rhs_y

!  using iterative refinement, and returning sol_y in rhs_y

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!   Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, A_ne, ldiag_x, ldiag_c_l, ldiag_c_u
      INTEGER, INTENT( IN ) :: itref_max, print_level
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( ldiag_x ) :: DIAG_X
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( dims%c_l_start : dims%c_u_end ) :: SCALE_C
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( ldiag_c_l : ldiag_c_u ) :: DIAG_C
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: RHS_y
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
                          DIMENSION( m ) :: SOL_y, BEST_y, RES_y
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: RES_x
      INTEGER, INTENT( IN ), DIMENSION( A_ne ) :: A_col
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( A_ne ) :: A_val
      TYPE ( SMT_type ), INTENT( IN ) :: K
      TYPE ( SILS_factors ), INTENT( INOUT ) :: FACTORS
      TYPE ( SILS_control ), INTENT( IN ) :: CNTL
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      TYPE ( LSQP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: i, it
      REAL ( KIND = wp ) :: res_norm, old_res
      TYPE ( SILS_sinfo ) :: SINFO

      old_res = MAXVAL( ABS( RHS_y ) )
      IF ( control%out > 0 .AND. print_level >= 3 )                            &
        WRITE( control%out, 2000 ) old_res

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, m ; SOL_y( i ) = RHS_y( i ) ; END DO
      ELSE
        SOL_y = RHS_y
      END IF

      CALL SILS_solve( K, FACTORS, SOL_y, CNTL, SINFO )
      inform%factorization_status = SINFO%flag
      IF ( SINFO%flag /= 0 ) THEN
        IF ( control%error > 0 .AND. print_level >= 1 )                        &
          WRITE( control%error, 2040 ) SINFO%flag
        inform%status = GALAHAD_error_solve
        RETURN
      END IF

      IF ( control%out > 0 .AND. print_level >= 3 )                            &
        WRITE( control%out, 2010 ) MAXVAL( ABS( SOL_y ) )

      IF ( itref_max > 0 ) THEN
        DO it = 1, itref_max

!  Form res_x = DIAG_X(inv) * A(trans) sol_y

          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, n ; RES_x( i ) = zero ; END DO
          ELSE
            RES_x = zero
          END IF
          CALL LSQP_AX( n, RES_x, m, A_ne,                         &
                        A_val, A_col,  A_ptr, m, SOL_y, '+T' )
          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, n ; RES_x( i ) = RES_x( i ) / DIAG_X( i ) ; END DO
          ELSE
            RES_x = RES_x / DIAG_X
          END IF

!  Form res_y = rhs_y - A * res_x - SCALE_C * DIAG_C(inv) * SCALE_C * sol_y

          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, dims%c_equality ; RES_y( i ) = RHS_y( i ) ; END DO
            DO i = dims%c_l_start,  dims%c_u_end
              RES_y( i ) = RHS_y( i ) - ( SCALE_C( i ) / DIAG_C( i ) ) *       &
                             SCALE_C( i ) * SOL_y( i )
            END DO
          ELSE
            RES_y( : dims%c_equality ) = RHS_y( : dims%c_equality )
            RES_y( dims%c_l_start : ) = RHS_y( dims%c_l_start : )              &
                   - ( SCALE_C / DIAG_C ) * SCALE_C * SOL_y( dims%c_l_start : )
          END IF
          CALL LSQP_AX( m, RES_y, m, A_ne, A_val, A_col, A_ptr,      &
                       n, RES_x, '- ' )
  
          IF ( control%out > 0 .AND. print_level >= 3 )                        &
            WRITE( control%out, 2000 ) MAXVAL( ABS( RES_y ) )
  
          res_norm = MAXVAL( ABS( RES_y ) )
          IF ( it > 1 ) THEN
            IF ( res_norm < old_res ) THEN
              old_res = res_norm
              IF ( control%array_syntax_worse_than_do_loop ) THEN
                DO i = 1, m ; BEST_y( i ) = SOL_y( i ) ; END DO
              ELSE
                BEST_y = SOL_y
              END IF
            ELSE
              IF ( control%array_syntax_worse_than_do_loop ) THEN
                DO i = 1, m ; SOL_y( i ) = BEST_y( i ) ; END DO
              ELSE
                SOL_y = BEST_y
              END IF
              EXIT
            END IF
          ELSE
            IF ( control%array_syntax_worse_than_do_loop ) THEN
              DO i = 1, m ; BEST_y( i ) = SOL_y( i ) ; END DO
            ELSE
              BEST_y = SOL_y
            END IF
          END IF
     
          CALL SILS_solve( K, FACTORS, RES_y, CNTL, SINFO )
          inform%factorization_status = SINFO%flag
          IF ( SINFO%flag /= 0 ) THEN
            IF ( control%error > 0 .AND. print_level >= 1 )                    &
              WRITE( control%error, 2040 ) SINFO%flag
            inform%status = GALAHAD_error_solve
            RETURN
          END IF
          IF ( control%out > 0 .AND. print_level >= 3 )                        &
            WRITE( control%out, 2020 ) MAXVAL( ABS( RES_y ) )
          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, m ; SOL_y( i ) = SOL_y( i ) + RES_y( i ) ; END DO
          ELSE
            SOL_y = SOL_y + RES_y
          END IF

!  Print final residuals if required

          IF ( control%out > 0 .AND. print_level >= 3                          &
              .AND. it == itref_max ) THEN
   
!  Form res_x = DIAG_X(inv) * A(trans) sol_y

            IF ( control%array_syntax_worse_than_do_loop ) THEN
              DO i = 1, n ; RES_x( i ) = zero ; END DO
            ELSE
              RES_x = zero
            END IF
            CALL LSQP_AX( n, RES_x, m, A_ne, A_val, A_col, A_ptr,    &
                          m, SOL_y, '+T' )
            IF ( control%array_syntax_worse_than_do_loop ) THEN
              DO i = 1, n ; RES_x( i ) = RES_x( i ) / DIAG_X( i ) ; END DO
            ELSE
              RES_x = RES_x / DIAG_X
            END IF

!  Form res_y = rhs_y - A * res_x - SCALE_C * DIAG_C(inv) * SCALE_C * sol_y

            IF ( control%array_syntax_worse_than_do_loop ) THEN
              DO i = 1, dims%c_equality ; RES_y( i ) = RHS_y( i ) ; END DO
              DO i = dims%c_l_start,  dims%c_u_end
                RES_y( i ) = RHS_y( i ) - ( SCALE_C( i ) / DIAG_C( i ) ) *     &
                               SCALE_C( i ) * SOL_y( i )
              END DO
            ELSE
              RES_y( : dims%c_equality ) = RHS_y( : dims%c_equality )
              RES_y( dims%c_l_start : ) = RHS_y( dims%c_l_start : )            &
                    - ( SCALE_C / DIAG_C ) * SCALE_C * SOL_y( dims%c_l_start : )
            END IF
            CALL LSQP_AX( m, RES_y, m, A_ne, A_val, A_col, A_ptr,    &
                          n, RES_x, '- ' )

            WRITE( control%out, 2000 ) MAXVAL( ABS( RES_y ) )
          END IF
        END DO
      ELSE

!  Print residuals if required

        IF ( control%out > 0 .AND. print_level >= 3 ) THEN
   
!  Form res_x = DIAG_X(inv) * A(trans) sol_y

          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, n ; RES_x( i ) = zero ; END DO
          ELSE
            RES_x = zero
          END IF
          CALL LSQP_AX( n, RES_x, m, A_ne, A_val, A_col, A_ptr,      &
                       m, SOL_y, '+T' )
          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, n ; RES_x( i ) = RES_x( i ) / DIAG_X( i ) ; END DO
          ELSE
            RES_x = RES_x / DIAG_X
          END IF

!  Form res_y = rhs_y - A * res_x - DIAG_C(inv) sol_y

          IF ( control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, dims%c_equality ; RES_y( i ) = RHS_y( i ) ; END DO
            DO i = dims%c_l_start,  dims%c_u_end
              RES_y( i ) = RHS_y( i ) - ( SCALE_C( i ) / DIAG_C( i ) ) *       &
                             SCALE_C( i ) * SOL_y( i )
            END DO
          ELSE
            RES_y( : dims%c_equality ) = RHS_y( : dims%c_equality )
            RES_y( dims%c_l_start : ) = RHS_y( dims%c_l_start : )              &
                  - ( SCALE_C / DIAG_C ) * SCALE_C * SOL_y( dims%c_l_start : )
          END IF
          CALL LSQP_AX( m, RES_y, m, A_ne, A_val, A_col, A_ptr,                &
                        n, RES_x, '- ' )
      
          WRITE( control%out, 2000 ) MAXVAL( ABS( RES_y ) )
        END IF
      END IF

      IF ( control%out > 0 .AND. print_level >= 3 )                            &
        WRITE( control%out, 2010 ) MAXVAL( ABS( SOL_y ) )

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 1, m ; RHS_y( i ) = SOL_y( i ) ; END DO
      ELSE
        RHS_y = SOL_y
      END IF
      inform%status = GALAHAD_ok
      RETURN

!  Non-executable statements

 2000 FORMAT( '  res(y) ', ES12.4 )
 2010 FORMAT( '  sol(y) ', ES12.4 )
 2020 FORMAT( ' dsol(y) ', ES12.4 )
 2040 FORMAT( '   **  Error return ', I3, ' from SILS_solve ' ) 

!  End of LSQP_block_iterative_refinement

      END SUBROUTINE LSQP_block_iterative_refinement

!-*-*-*-*-*-*-   L S Q P _ R E S I D U A L   S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE LSQP_residual( dims, n, m, l_res, A_ne, A_val, A_col, A_ptr,  &
                                DX, DC, DY, RHS_x, RHS_c, RHS_y, RES,          &
                                BARRIER_X, BARRIER_C, SCALE_C, errorg, errorc, &
                                print_level, control, Hessian_kind, WEIGHT )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Compute the residual of the linear system

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, A_ne, l_res, Hessian_kind, print_level
      REAL( KIND = wp ), INTENT( OUT ) :: errorg, errorc
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( l_res ) :: RES
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: RHS_x, DX
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_free + 1 : n ) :: BARRIER_X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: RHS_y, DY
      REAL ( KIND = wp ), INTENT( IN ),                                        &
           DIMENSION( dims%c_l_start : m ) :: RHS_c, DC, BARRIER_C, SCALE_C
      INTEGER, INTENT( IN ), DIMENSION( A_ne ) :: A_col
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( A_ne ) :: A_val
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ), OPTIONAL :: WEIGHT
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        

!  Local variables

      INTEGER :: i
      REAL ( KIND = wp ) :: res_tol

      res_tol = epsmch ** 0.5

!  Initalize RES as the zero vector

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = dims%y_s, dims%y_e ; RES( i ) = zero ; END DO
      ELSE
        RES( dims%y_s : dims%y_e ) = zero
      END IF 

!  Remember the barrier terms, and any diagonal perturbations

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        IF ( Hessian_kind == 0 ) THEN
          DO i = 1, dims%x_free ; RES( i ) = zero ; END DO
          DO i = dims%x_free + 1, dims%x_e 
            RES( i ) = BARRIER_X( i ) * DX( i ) ; END DO
        ELSE IF ( Hessian_kind == 1 ) THEN
          DO i = 1, dims%x_free
            RES( i ) = DX( i ) ; END DO
          DO i = dims%x_free + 1, dims%x_e 
            RES( i ) = ( one + BARRIER_X( i ) ) * DX( i ) ; END DO
        ELSE
          DO i = 1, dims%x_free
            RES( i ) = ( WEIGHT( i ) ** 2 ) * DX( i ) ; END DO
          DO i = dims%x_free + 1, dims%x_e 
            RES( i ) = ( WEIGHT( i ) ** 2 + BARRIER_X( i ) ) * DX( i ) ; END DO
        END IF
        DO i = 0, dims%nc - 1
          RES( dims%c_s + i ) = BARRIER_C( dims%c_l_start + i ) *              &
            DC( dims%c_l_start + i ) ; END DO
      ELSE
        IF ( Hessian_kind == 0 ) THEN
          RES( : dims%x_free ) = zero
          RES( dims%x_free + 1 : dims%x_e ) =                                  &
               BARRIER_X * DX( dims%x_free + 1 : )
        ELSE IF ( Hessian_kind == 1 ) THEN
          RES( : dims%x_free ) = DX( : dims%x_free )
          RES( dims%x_free + 1 : dims%x_e ) =                                  &
              ( one + BARRIER_X ) * DX( dims%x_free + 1 : )
        ELSE
          RES( : dims%x_free ) =                                               &
             ( WEIGHT( : dims%x_free ) ** 2 ) * DX( : dims%x_free )
          RES( dims%x_free + 1 : dims%x_e ) =                                  &
             ( WEIGHT( dims%x_free + 1 : dims%x_e ) ** 2 +                     &
               BARRIER_X ) * DX( dims%x_free + 1 : )
        END IF
        RES( dims%c_s : dims%c_e ) = BARRIER_C * DC
      END IF

!  Include the contribution from A and A^T

      CALL LSQP_AX( n, RES( dims%x_s : dims%x_e ), m, A_ne, A_val, A_col,      &
                    A_ptr, m, DY, '+T' )
      CALL LSQP_AX( m, RES( dims%y_s : dims%y_e ), m, A_ne, A_val, A_col,      &
                    A_ptr, n, DX, '+ ' )

!  Include the contribution from the slack variables

      IF ( control%array_syntax_worse_than_do_loop ) THEN
        DO i = 0, dims%nc - 1
          RES( dims%c_s + i ) = RES( dims%c_s + i ) -                          &
              SCALE_C( dims%c_l_start + i ) * DY( dims%c_l_start + i )
          RES( dims%y_i + i ) = RES( dims%y_i + i ) -                          &
              SCALE_C( dims%c_l_start + i ) * DC( dims%c_l_start + i )  
        END DO
      ELSE
        RES( dims%c_s : dims%c_e ) =                                           &
          RES( dims%c_s : dims%c_e ) - SCALE_C * DY( dims%c_l_start : m )
        RES( dims%y_i : dims%y_e ) =                                           &
          RES( dims%y_i : dims%y_e ) - SCALE_C * DC
      END IF

!  Find the largest residual and component of the search direction

      IF ( control%out > 0 .AND. print_level >= 2 ) THEN
        errorg = MAX( MAXVAL( ABS( RES( dims%x_s : dims%x_e ) - RHS_x ) ),     &
                      MAXVAL( ABS( RES( dims%c_s : dims%c_e ) - RHS_c ) ) )
        IF ( print_level >= 4 ) THEN
          DO i = 1, n 
            IF ( ABS( RES( i ) - RHS_x( i ) ) > res_tol )                      &
              WRITE( control%out, 2010 ) 'X', i, RES( i ), RHS_x( i ) 
          END DO
          DO i = dims%c_l_start, dims%c_u_end
            IF ( ABS( RES( dims%c_b + i ) - RHS_c( i ) ) > res_tol )           &
              WRITE( control%out, 2010 ) 'C', i, RES( dims%c_b + i ), RHS_c( i )
          END DO
        END IF
        IF ( m > 0 ) THEN
          errorc = MAXVAL( ABS( RES( dims%y_s : dims%y_e ) - RHS_y ) )
        ELSE
          errorc = zero
        END IF
        IF ( print_level >= 4 ) THEN
          DO i = 1, m 
            IF ( ABS( RES( dims%y_s + i - 1 ) - RHS_y( i ) ) > res_tol )       &
              WRITE( control%out, 2010 )                                       &
                'C', I, RES( dims%y_s + i - 1 ), RHS_y( i )
          END DO 
        END IF
        WRITE( control%out, 2000 ) errorg, errorc, MAX( MAXVAL( ABS( DX ) ),   &
                                                        MAXVAL( ABS( DC ) ),   &
                                                        MAXVAL( ABS( DY ) ) )
      END IF
      RETURN

!  Non-executable statements

 2000 FORMAT( ' ',                                                            &
          /, '         ***  Max component of gradient  residuals = ', ES12.4, &
          /, '         ***  Max component of contraint residuals = ', ES12.4, &
          /, '         ***  Max component of search direction    = ', ES12.4 )
 2010 FORMAT( ' ', A1, '-residual', I6, ' lhs = ', ES12.4,' rhs = ', ES12.4 ) 

!  End of subroutine LSQP_residual

      END SUBROUTINE LSQP_residual

!-*-*-*-  L S Q P _ C O M P U T E _ M A X S T E P   S U B R O U T I N E  -*-*-*-

      SUBROUTINE LSQP_compute_maxstep( dims, n, m, Z_l, Z_u, DZ_l, DZ_u,       &
                                       X, DIST_X_l, DIST_X_u, DX,              &
                                       Y_l, Y_u, DY_l, DY_u,                   &
                                       DIST_C_l, DIST_C_u, DC,                 &
                                       gamma_b, gamma_f, nu, nbnds,            &
                                       alpha_max_b, alpha_max_f,               &
                                       print_level, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the maximum allowable stepsizes alpha_max_b, which balances the 
!  complementarity ie, such that
!
!      min (x-l)_i(z_l)_i - (gamma_b / nbds)( <x-l,z_l> + <x-u,z_u> ) >= 0
!       i                                           
!  and
!      min (x-u)_i(z_u)_i - (gamma_b / nbds)( <x-l,z_l> + <x-u,z_u> ) >= 0 ,
!       i                                           
!
!  and alpha_max_f, which favours feasibility over complementarity, 
!  ie, such that
!
!      <x-l,z_l> + <x-u,z_u> >= nu * gamma_f
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, nbnds, print_level
      REAL ( KIND = wp ), INTENT( IN ) :: gamma_b, gamma_f, nu
      REAL ( KIND = wp ), INTENT( OUT ) :: alpha_max_b, alpha_max_f 
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X, DX 
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_free + 1 : dims%x_l_end ) :: Z_l, DZ_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_l_start : dims%x_l_end ) :: DIST_X_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_u_start : n ) :: Z_u, DZ_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_u_start : dims%x_u_end ) :: DIST_X_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_l_start : dims%c_l_end ) :: Y_l, DY_l, DIST_C_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_u_start : dims%c_u_end ) :: Y_u, DY_u, DIST_C_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( dims%c_l_start : m ) :: DC 
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      TYPE ( LSQP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: i, nroots 

!  Local variables

      REAL ( KIND = wp ) :: compc, compl, compq, coef0, coef1, coef2
      REAL ( KIND = wp ) :: coef0_f, coef1_f, coef2_f
      REAL ( KIND = wp ) :: root1, root2, tol, alpha, alp, nu_gamma_f

      alpha_max_b = infinity ; alpha_max_f = infinity 
      inform%status = GALAHAD_ok
      IF ( nbnds == 0 ) RETURN
      tol = epsmch ** 0.75 

!  ================================================
!             part to compute alpha_max_b
!  ================================================

!  Compute the coefficients for the quadratic expression
!  for the overall complementarity

      coef0_f = zero ; coef1_f = zero ; coef2_f = zero 
      DO i = dims%x_free + 1, dims%x_l_start - 1
        coef0_f = coef0_f + X( i ) * Z_l( i ) 
        coef1_f = coef1_f + X( i ) * DZ_l( i ) + DX( i ) * Z_l( i ) 
        coef2_f = coef2_f + DX( i ) * DZ_l( i ) 
      END DO 
      DO i = dims%x_l_start, dims%x_l_end
        coef0_f = coef0_f + DIST_X_l( i ) * Z_l( i ) 
        coef1_f = coef1_f + DIST_X_l( i ) * DZ_l( i ) + DX( i ) * Z_l( i ) 
        coef2_f = coef2_f + DX( i ) * DZ_l( i ) 
      END DO 
      DO i = dims%x_u_start, dims%x_u_end
        coef0_f = coef0_f - DIST_X_u( i ) * Z_u( i ) 
        coef1_f = coef1_f - DIST_X_u( i ) * DZ_u( i ) + DX( i ) * Z_u( i ) 
        coef2_f = coef2_f + DX( i ) * DZ_u( i ) 
      END DO 
      DO i = dims%x_u_end + 1, n
        coef0_f = coef0_f + X( i ) * Z_u( i ) 
        coef1_f = coef1_f + X( i ) * DZ_u( i ) + DX( i ) * Z_u( i ) 
        coef2_f = coef2_f + DX( i ) * DZ_u( i ) 
      END DO 
      DO i = dims%c_l_start, dims%c_l_end
        coef0_f = coef0_f + DIST_C_l( i ) * Y_l( i ) 
        coef1_f = coef1_f + DIST_C_l( i ) * DY_l( i ) + DC( i ) * Y_l( i ) 
        coef2_f = coef2_f + DC( i ) * DY_l( i ) 
      END DO 
      DO i = dims%c_u_start, dims%c_u_end
        coef0_f = coef0_f - DIST_C_u( i ) * Y_u( i ) 
        coef1_f = coef1_f - DIST_C_u( i ) * DY_u( i ) + DC( i ) * Y_u( i ) 
        coef2_f = coef2_f + DC( i ) * DY_u( i ) 
      END DO 

!  Scale these coefficients

      compc = - gamma_b * coef0_f / nbnds ; compl = - gamma_b * coef1_f / nbnds
      compq = - gamma_b * coef2_f / nbnds
!     write(6,"( ' gamma_b ', ES12.4, I6 )" ) gamma_b, nbnds
!     write( 6, "( 3ES10.2 )" )  compq, compl, compc

!  Compute the coefficients for the quadratic expression
!  for the individual complementarity

      DO i = dims%x_free + 1, dims%x_l_start - 1
        coef0 = compc + X( i ) * Z_l( i ) 
        coef1 = compl + X( i ) * DZ_l( i ) + DX( i ) * Z_l( i ) 
        coef2 = compq + DX( i ) * DZ_l( i ) 
        coef0 = MAX( coef0, zero )
        CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
        IF ( nroots == 2 ) THEN 
          IF ( coef2 > zero ) THEN 
            IF ( root2 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = root2 
          END IF 
        ELSE IF ( nroots == 1 ) THEN 
          IF ( root1 > zero ) THEN 
            alpha = root1 
          ELSE 
            alpha = infinity 
          END IF 
        ELSE 
          alpha = infinity 
        END IF 
        IF ( alpha < alpha_max_b ) alpha_max_b = alpha
      END DO

      DO i = dims%x_l_start, dims%x_l_end
        coef0 = compc + DIST_X_l( i ) * Z_l( i ) 
        coef1 = compl + DIST_X_l( i ) * DZ_l( i ) + DX( i ) * Z_l( i ) 
        coef2 = compq + DX( i ) * DZ_l( i ) 
        coef0 = MAX( coef0, zero )
        CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
        IF ( nroots == 2 ) THEN 
          IF ( coef2 > zero ) THEN 
            IF ( root2 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = root2 
          END IF 
        ELSE IF ( nroots == 1 ) THEN 
          IF ( root1 > zero ) THEN 
            alpha = root1 
          ELSE 
            alpha = infinity 
          END IF 
        ELSE 
          alpha = infinity 
        END IF 
        IF ( alpha < alpha_max_b ) alpha_max_b = alpha 
      END DO

      DO i = dims%x_u_start, dims%x_u_end
        coef0 = compc - DIST_X_u( i ) * Z_u( i ) 
        coef1 = compl - DIST_X_u( i ) * DZ_u( i ) + DX( i ) * Z_u( i ) 
        coef2 = compq + DX( i ) * DZ_u( i ) 
        coef0 = MAX( coef0, zero )
        CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
        IF ( nroots == 2 ) THEN 
          IF ( coef2 > zero ) THEN 
            IF ( root2 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = root2 
          END IF 
        ELSE IF ( nroots == 1 ) THEN 
          IF ( root1 > zero ) THEN 
            alpha = root1 
          ELSE 
            alpha = infinity 
          END IF 
        ELSE 
          alpha = infinity 
        END IF 
        IF ( alpha < alpha_max_b ) alpha_max_b = alpha 
      END DO 

      DO i = dims%x_u_end + 1, n
        coef0 = compc + X( i ) * Z_u( i ) 
        coef1 = compl + X( i ) * DZ_u( i ) + DX( i ) * Z_u( i ) 
        coef2 = compq + DX( i ) * DZ_u( i ) 
        coef0 = MAX( coef0, zero )
        CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
        IF ( nroots == 2 ) THEN 
          IF ( coef2 > zero ) THEN 
            IF ( root2 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = root2 
          END IF 
        ELSE IF ( nroots == 1 ) THEN 
          IF ( root1 > zero ) THEN 
            alpha = root1 
          ELSE 
            alpha = infinity 
          END IF 
        ELSE 
          alpha = infinity 
        END IF 
        IF ( alpha < alpha_max_b ) alpha_max_b = alpha 
      END DO 

      DO i = dims%c_l_start, dims%c_l_end
        coef0 = compc + DIST_C_l( i ) * Y_l( i ) 
        coef1 = compl + DIST_C_l( i ) * DY_l( i ) + DC( i ) * Y_l( i ) 
        coef2 = compq + DC( i ) * DY_l( i ) 
        coef0 = MAX( coef0, zero )
        CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
        IF ( nroots == 2 ) THEN 
          IF ( coef2 > zero ) THEN 
            IF ( root2 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = root2 
          END IF 
        ELSE IF ( nroots == 1 ) THEN 
          IF ( root1 > zero ) THEN 
            alpha = root1 
          ELSE 
            alpha = infinity 
          END IF 
        ELSE 
          alpha = infinity 
        END IF 
        IF ( alpha < alpha_max_b ) alpha_max_b = alpha 
      END DO

      DO i = dims%c_u_start, dims%c_u_end
        coef0 = compc - DIST_C_u( i ) * Y_u( i ) 
        coef1 = compl - DIST_C_u( i ) * DY_u( i ) + DC( i ) * Y_u( i ) 
        coef2 = compq + DC( i ) * DY_u( i ) 
        coef0 = MAX( coef0, zero )
        CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
!       write( 6, "( 3ES10.2, 2ES22.14 )" )  coef2, coef1, coef0, root1, root2
        IF ( nroots == 2 ) THEN 
          IF ( coef2 > zero ) THEN 
            IF ( root2 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = root2 
          END IF 
        ELSE IF ( nroots == 1 ) THEN 
          IF ( root1 > zero ) THEN 
            alpha = root1 
          ELSE 
            alpha = infinity 
          END IF 
        ELSE 
          alpha = infinity 
        END IF 
        IF ( alpha < alpha_max_b ) alpha_max_b = alpha 
      END DO 

      IF ( - compc <= epsmch ** 0.75 ) alpha_max_b = 0.99_wp * alpha_max_b

!  An error has occured. Investigate

      IF ( alpha_max_b <= zero ) THEN 
        IF ( control%out > 0 .AND. print_level >= 2 )                          &
          WRITE( control%out, 2020 ) alpha_max_b
        DO i = dims%x_free + 1, dims%x_l_start - 1
          coef0 = compc + X( i ) * Z_l( i ) 
          coef1 = compl + X( i ) * DZ_l( i ) + DX( i ) * Z_l( i ) 
          coef2 = compq + DX( i ) * DZ_l( i ) 
          CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
          IF ( nroots == 2 ) THEN 
            IF ( coef2 > zero ) THEN 
               IF ( root2 > zero ) THEN 
                  alpha = root1 
               ELSE 
                  alpha = infinity 
               END IF 
            ELSE 
               alpha = root2 
            END IF 
          ELSE IF ( nroots == 1 ) THEN 
            IF ( root1 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = infinity 
          END IF 
          IF ( alpha == alpha_max_b ) THEN
            IF ( control%out > 0 .AND. print_level >= 2 ) THEN
               IF ( nroots == 2 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'X', i, 'L', coef0, coef1, coef2, root1, root2
               ELSE IF ( nroots == 1 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                  'X', i, 'L', coef0, coef1, coef2, root1 
               ELSE 
                 WRITE( control%out, 2000 ) 'X', i, 'L', coef0, coef1, coef2 
               END IF 
               WRITE( control%out, 2010 ) 'X', i, 'L', alpha 
            END IF 
          END IF 
        END DO
        DO i = dims%x_l_start, dims%x_l_end
          coef0 = compc + DIST_X_l( i ) * Z_l( i ) 
          coef1 = compl + DIST_X_l( i ) * DZ_l( i ) + DX( i ) * Z_l( i ) 
          coef2 = compq + DX( i ) * DZ_l( i ) 
          CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
          IF ( nroots == 2 ) THEN 
            IF ( coef2 > zero ) THEN 
               IF ( root2 > zero ) THEN 
                  alpha = root1 
               ELSE 
                  alpha = infinity 
               END IF 
            ELSE 
               alpha = root2 
            END IF 
          ELSE IF ( nroots == 1 ) THEN 
            IF ( root1 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = infinity 
          END IF 
          IF ( alpha == alpha_max_b ) THEN
            IF ( control%out > 0 .AND. print_level >= 2 ) THEN
               IF ( nroots == 2 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'X', i, 'L', coef0, coef1, coef2, root1, root2
               ELSE IF ( nroots == 1 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'X', i, 'L', coef0, coef1, coef2, root1 
               ELSE 
                 WRITE( control%out, 2000 ) 'X', i, 'L', coef0, coef1, coef2 
               END IF 
               WRITE( control%out, 2010 ) 'X', i, 'L', alpha 
            END IF 
          END IF 
        END DO
        DO i = dims%x_u_start, dims%x_u_end
          coef0 = compc - DIST_X_u( i ) * Z_u( i ) 
          coef1 = compl - DIST_X_u( i ) * DZ_u( i ) + DX( i ) * Z_u( i ) 
          coef2 = compq + DX( i ) * DZ_u( i ) 
          CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
          IF ( nroots == 2 ) THEN 
            IF ( coef2 > zero ) THEN 
               IF ( root2 > zero ) THEN 
                  alpha = root1 
               ELSE 
                  alpha = infinity 
               END IF 
            ELSE 
               alpha = root2 
            END IF 
          ELSE IF ( nroots == 1 ) THEN 
            IF ( root1 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = infinity 
          END IF 
          IF ( alpha == alpha_max_b ) THEN
            IF ( control%out > 0 .AND. print_level >= 2 ) THEN
               IF ( nroots == 2 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'X', i, 'U', coef0, coef1, coef2, root1, root2
               ELSE IF ( nroots == 1 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'X', i, 'U', coef0, coef1, coef2, root1 
               ELSE 
                 WRITE( control%out, 2000 ) 'X', i, 'U', coef0, coef1, coef2 
               END IF 
               WRITE( control%out, 2010 ) 'X', i, 'U', alpha 
            END IF 
          END IF 
        END DO 
        DO i = dims%x_u_end + 1, n
          coef0 = compc + X( i ) * Z_u( i ) 
          coef1 = compl + X( i ) * DZ_u( i ) + DX( i ) * Z_u( i ) 
          coef2 = compq + DX( i ) * DZ_u( i ) 
          CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
          IF ( nroots == 2 ) THEN 
            IF ( coef2 > zero ) THEN 
               IF ( root2 > zero ) THEN 
                  alpha = root1 
               ELSE 
                  alpha = infinity 
               END IF 
            ELSE 
               alpha = root2 
            END IF 
          ELSE IF ( nroots == 1 ) THEN 
            IF ( root1 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = infinity 
          END IF 
          IF ( alpha == alpha_max_b ) THEN
            IF ( control%out > 0 .AND. print_level >= 2 ) THEN
               IF ( nroots == 2 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'X', i, 'U', coef0, coef1, coef2, root1, root2
               ELSE IF ( nroots == 1 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'X', i, 'U', coef0, coef1, coef2, root1 
               ELSE 
                 WRITE( control%out, 2000 ) 'X', i, 'U', coef0, coef1, coef2 
               END IF 
               WRITE( control%out, 2010 ) 'X', i, 'U', alpha 
            END IF 
          END IF 
        END DO 
        DO i = dims%c_l_start, dims%c_l_end
          coef0 = compc + DIST_C_l( i ) * Y_l( i ) 
          coef1 = compl + DIST_C_l( i ) * DY_l( i ) + DC( i ) * Y_l( i ) 
          coef2 = compq + DC( i ) * DY_l( i ) 
          CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
          IF ( nroots == 2 ) THEN 
            IF ( coef2 > zero ) THEN 
               IF ( root2 > zero ) THEN 
                  alpha = root1 
               ELSE 
                  alpha = infinity 
               END IF 
            ELSE 
               alpha = root2 
            END IF 
          ELSE IF ( nroots == 1 ) THEN 
            IF ( root1 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = infinity 
          END IF 
          IF ( alpha == alpha_max_b ) THEN
            IF ( control%out > 0 .AND. print_level >= 2 ) THEN
               IF ( nroots == 2 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'C', i, 'L', coef0, coef1, coef2, root1, root2
               ELSE IF ( nroots == 1 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'C', i, 'L', coef0, coef1, coef2, root1 
               ELSE 
                 WRITE( control%out, 2000 ) 'C', i, 'L', coef0, coef1, coef2 
               END IF 
               WRITE( control%out, 2010 ) 'C', i, 'L', alpha 
            END IF 
          END IF 
        END DO
        DO i = dims%c_u_start, dims%c_u_end
          coef0 = compc - DIST_C_u( i ) * Y_u( i ) 
          coef1 = compl - DIST_C_u( i ) * DY_u( i ) + DC( i ) * Y_u( i ) 
          coef2 = compq + DC( i ) * DY_u( i ) 
          CALL ROOTS_quadratic( coef0, coef1, coef2, tol, nroots, root1, root2 ) 
          IF ( nroots == 2 ) THEN 
            IF ( coef2 > zero ) THEN 
               IF ( root2 > zero ) THEN 
                  alpha = root1 
               ELSE 
                  alpha = infinity 
               END IF 
            ELSE 
               alpha = root2 
            END IF 
          ELSE IF ( nroots == 1 ) THEN 
            IF ( root1 > zero ) THEN 
               alpha = root1 
            ELSE 
               alpha = infinity 
            END IF 
          ELSE 
            alpha = infinity 
          END IF 
          IF ( alpha == alpha_max_b ) THEN
            IF ( control%out > 0 .AND. print_level >= 2 ) THEN
               IF ( nroots == 2 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'C', i, 'U', coef0, coef1, coef2, root1, root2
               ELSE IF ( nroots == 1 ) THEN 
                 WRITE( control%out, 2000 )                                    &
                   'C', i, 'U', coef0, coef1, coef2, root1 
               ELSE 
                 WRITE( control%out, 2000 ) 'C', i, 'U', coef0, coef1, coef2 
               END IF 
               WRITE( control%out, 2010 ) 'C', i, 'U', alpha 
            END IF 
          END IF 
        END DO 

        DO i = dims%x_free + 1, dims%x_l_start - 1
          coef0 = X( i ) * Z_l( i ) 
          coef1 = DX( i ) * Z_l( i ) + X( i ) * DZ_l( i )
          coef2 = DX( i ) * DZ_l( i ) 
          alp = alpha_max_b ; alpha = coef0 + alp * ( coef1 + alp * coef2 ) 
          IF ( control%out > 0 .AND. print_level >= 2 )                        &
            WRITE( control%out, 2030 ) 'X', i, 'L', alp, alpha 
        END DO
        DO i = dims%x_l_start, dims%x_l_end
          coef0 = DIST_X_l( i ) * Z_l( i ) 
          coef1 = DX( i ) * Z_l( i ) + DIST_X_l( i ) * DZ_l( i )
          coef2 = DX( i ) * DZ_l( i ) 
          alp = alpha_max_b ; alpha = coef0 + alp * ( coef1 + alp * coef2 ) 
          IF ( control%out > 0 .AND. print_level >= 2 )                        &
            WRITE( control%out, 2030 ) 'X', i, 'L', alp, alpha 
        END DO
        DO i = dims%x_u_start, dims%x_u_end
          coef0 = - DIST_X_u( i ) * Z_u( i ) 
          coef1 = DX( i ) * Z_u( i ) - DIST_X_u( i ) * DZ_u( i )
          coef2 = DX( i ) * DZ_u( i ) 
          alp = alpha_max_b ; alpha = coef0 + alp * ( coef1 + alp * coef2 ) 
          IF ( control%out > 0 .AND. print_level >= 2 )                        &
            WRITE( control%out, 2030 ) 'X', i, 'U', alp, alpha 
        END DO 
        DO i = dims%x_u_end + 1, n
          coef0 = X( i ) * Z_u( i ) 
          coef1 = DX( i ) * Z_u( i ) + X( i ) * DZ_u( i )
          coef2 = DX( i ) * DZ_u( i ) 
          alp = alpha_max_b ; alpha = coef0 + alp * ( coef1 + alp * coef2 ) 
          IF ( control%out > 0 .AND. print_level >= 2 )                        &
            WRITE( control%out, 2030 ) 'X', i, 'U', alp, alpha 
        END DO 
        DO i = dims%c_l_start, dims%c_l_end
          coef0 = DIST_C_l( i ) * Y_l( i ) 
          coef1 = DC( i ) * Y_l( i ) + DIST_C_l( i ) * DY_l( i )
          coef2 = DC( i ) * DY_l( i ) 
          alp = alpha_max_b ; alpha = coef0 + alp * ( coef1 + alp * coef2 ) 
          IF ( control%out > 0 .AND. print_level >= 2 )                        &
            WRITE( control%out, 2030 ) 'C', i, 'L', alp, alpha 
        END DO
        DO i = dims%c_u_start, dims%c_u_end
          coef0 = - DIST_C_u( i ) * Y_u( i ) 
          coef1 = DC( i ) * Y_u( i ) - DIST_C_u( i ) * DY_u( i )
          coef2 = DC( i ) * DY_u( i ) 
          alp = alpha_max_b ; alpha = coef0 + alp * ( coef1 + alp * coef2 ) 
          IF ( control%out > 0 .AND. print_level >= 2 )                        &
            WRITE( control%out, 2030 ) 'C', i, 'U', alp, alpha 
        END DO 
        alp = alpha_max_b ; alpha = compc + alp * ( compl + alp * compq ) 
        IF ( control%out > 0 .AND. print_level >= 2 ) THEN
          WRITE( control%out, 2040 ) alpha 
          WRITE( control%out, 2020 ) alpha_max_b 
        END IF
        WRITE( control%out, "( ' -ve step, no further progress possible ' )" )
        inform%status = GALAHAD_error_tiny_step
        RETURN
      END IF 

!  ================================================
!             part to compute alpha_max_f
!  ================================================

      nu_gamma_f = nu * gamma_f

!  Compute the coefficients for the quadratic expression
!  for the overall complementarity, remembering to first 
!  subtract the term for the feasibility

      coef0_f = coef0_f - nu_gamma_f
      coef1_f = coef1_f + nu_gamma_f

!  Compute the coefficients for the quadratic expression
!  for the individual complementarity
!
      CALL ROOTS_quadratic( coef0_f, coef1_f, coef2_f, tol,                    &
                            nroots, root1, root2 )
      IF ( nroots == 2 ) THEN 
        IF ( coef2_f > zero ) THEN 
          IF ( root2 > zero ) THEN 
            alpha = root1 
          ELSE 
            alpha = infinity 
          END IF 
        ELSE 
          alpha = root2 
        END IF 
      ELSE IF ( nroots == 1 ) THEN 
        IF ( root1 > zero ) THEN 
          alpha = root1 
        ELSE 
          alpha = infinity 
        END IF 
      ELSE 
        alpha = infinity 
      END IF 
      IF ( alpha < alpha_max_f ) alpha_max_f = alpha 
      IF ( - compc <= epsmch ** 0.75 ) alpha_max_f = 0.99_wp * alpha_max_f

      RETURN
  
!  Non-executable statements

 2000 FORMAT( A1, I6, A1,' coefs', 3ES12.4,' roots', 2ES12.4 ) 
 2010 FORMAT( A1, I6, A1,' alpha', ES12.4 ) 
 2020 FORMAT( ' alpha_min ', ES12.4 ) 
 2030 FORMAT( A1, I6, A1,' value at ', ES12.4,' = ', ES12.4 ) 
 2040 FORMAT( ' .vs. ', ES12.4 ) 

!  End of subroutine LSQP_compute_maxstep

      END SUBROUTINE LSQP_compute_maxstep

!-*-*-*-*-*-*-*-   L S Q P _ F O R M _ S _ C   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE LSQP_form_Schur_complement(                                   &
                          dims, n, m, Abycol_ne, Abycol_val, Abycol_row,       &
                          Abycol_ptr, DIAG_X, SCALE_C, DIAG_C, ls, S_val,      &
                          S_row, S_colptr, nes, ierr, col_count, row_count,    &
                          reals, control_error, control_print_level )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  This subroutine computes the Schur-complement matrix 
!    S = A DIAG_X(inv) A(trans) + SCALE_C * DIAG_C(inv) * SCALE_C
!  given the matrix A and vector of weights DIAG_X and DIAG_C

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Arguments:

!  m is the number of rows in matrix A.
!  n is the number of columns in matrix A, and the number of rows and
!    columns in matrix S.
!  Abycol_ne is the number of entries in matrix A. It is the length of 
!     arrays Abycol_val and Abycol_row.
!  Abycol_val holds the values of the entries in matrix A.
!  Abycol_row holds the row indices of the corresponding entries in A
!  Abycol_ptr points to the start of each column of A.
!  DIAG_X, DIAG_C are vectors of weights
!  ls is the length of arrays S and S_row.
!  S_val is set to the values of the entries in the lower triangle of S, 
!   ordered by columns
!  S_row is set to the row numbers of the entries in the lower triangle of S
!  S_colptr is set to point to the start of each column of S.
!  ierr is set by the routine to the code specifying the type of input
!       error:
!            = 0    no input errors;
!            = 1    n less than or equal to zero;
!            = 2    m less than or equal to zero;
!            = 3    ls less than Abycol_ne+m;
!            = 4    ls less than the number of entries in S (the
!                   number required is in nes);
!  nes is the actual number of entries in S
!  COL_count is the number of entries in each column of A
!  ROW_count is the number of entries in each row in A and, later,
!            the number of entries in each row of the lower triangle of S

!  Based on Iain Duff's routine MC35

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, ls, Abycol_ne
      INTEGER, INTENT( IN ) :: control_error, control_print_level
      INTEGER, INTENT( OUT ) :: nes, ierr
      LOGICAL, INTENT( IN ) :: reals
      INTEGER, INTENT( INOUT ), DIMENSION( Abycol_ne ) :: Abycol_row
      INTEGER, INTENT( INOUT ), DIMENSION( n + 1 ) :: Abycol_ptr
      INTEGER, INTENT( OUT ), DIMENSION( ls )  :: S_row
      INTEGER, INTENT( OUT ), DIMENSION( m + 1 )  :: S_colptr
      INTEGER, INTENT( OUT ), DIMENSION( n )  :: COL_count
      INTEGER, INTENT( OUT ), DIMENSION( m )  :: ROW_count
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n )  :: DIAG_X
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( dims%c_l_start : dims%c_u_end ) :: SCALE_C
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( dims%c_l_start : dims%c_u_end ) :: DIAG_C
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( Abycol_ne )  :: Abycol_val
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( ls )  :: S_val

!  Local variables

      INTEGER :: i, i1, i2, ii, ij, ipos, j, j1, j2, jj
      INTEGER :: icount, idist, k, k1, k2, kk
      REAL ( KIND = wp ) :: amult
      LOGICAL :: error

!  Check for input errors

      IF ( n <= 0 ) THEN
        ierr = 1
        IF ( control_error > 0 .AND. control_print_level >= 2 )                &
          WRITE( control_error, 2010 ) n
        GO TO 900
      END IF

      IF ( m <= 0 ) THEN
        ierr = 2
        IF ( control_error > 0 .AND. control_print_level >= 2 )                &
          WRITE( control_error, 2020 ) m
        GO TO 900
      END IF

      IF ( ls < Abycol_ne + m ) THEN
        ierr = 3
        IF ( control_error > 0 .AND. control_print_level >= 2 )                &
          WRITE( control_error, 2030 ) ls, Abycol_ne + m
        GO TO 900
      END IF

      ierr = 0

!  Set up column counts

!     COL_count( : n ) =                                                  &
!       Abycol_ptr( 2 : n + 1 ) - Abycol_ptr( : n )
      DO i = 1, n
        COL_count( i ) = Abycol_ptr( i + 1 ) - Abycol_ptr( i )
      END DO

!  Set up row counts

      ROW_count( : m ) = 0
      DO i = 1, Abycol_ne
        j = Abycol_row( i )
        ROW_count( j ) = ROW_count( j ) + 1
      END DO

!  Set up column pointer

      S_colptr( 1 ) = ls - Abycol_ne + 1
      DO i = 1, m
        S_colptr( i + 1 ) = S_colptr( i ) + ROW_count( i )
      END DO

!  Generate structure of A by columns in S/S_row

      j1 = 1
      DO i = 1, n
        j2 = j1 + COL_count( i ) - 1
        DO jj = j1, j2
          j = Abycol_row( jj )
          ipos = S_colptr( j )
          S_val( ipos ) = Abycol_val( jj )
          S_row( ipos ) = i
          S_colptr( j ) = ipos + 1
        END DO
        j1 = j2 + 1
      END DO

!  Regenerate structure of A by columns but in row order within each column.
!  Order the row indices in Abycol_row using S_val/S_row to create 
!  Abycol_val/Abycol_row

      i1 = ls - Abycol_ne + 1
      DO j = 1, m
        ij = m - j + 1

!  Reset S_colptr after the last loop altered it

        IF ( ij > 1 ) S_colptr( ij ) = S_colptr( ij - 1 )
        i2 = i1 + ROW_count( j ) - 1
        DO ii = i1, i2
          i = S_row( ii )
          ipos = Abycol_ptr( i )
          Abycol_row( ipos ) = j
          Abycol_val( ipos ) = S_val( ii )
          Abycol_ptr( i ) = ipos + 1
        END DO
        i1 = i2 + 1
      END DO
      S_colptr( 1 ) = 1 + ls - Abycol_ne

!  Reset Abycol_ptr since it has been altered

      DO j = n, 2, - 1
        Abycol_ptr( j ) = Abycol_ptr( j - 1 )
      END DO
      Abycol_ptr( 1 ) = 1

!  Now assemble S

      error = .FALSE.
      nes = 0
      icount = 0
      idist = 0
      k1 = S_colptr( 1 )

!  Form row i of S from a linear combination of rows of A

      DO i = 1, m

!  Leave room for possible diagonal if row i of A is empty

        IF ( ROW_count( i ) == 0 ) THEN
          IF ( i >= dims%c_l_start ) THEN
            IF ( error ) THEN
              icount = icount + 1
            ELSE
              ROW_count( i ) = 1
              nes = nes + 1
              IF ( nes <= ls ) THEN
                S_row( nes ) = i
                IF ( reals ) S_val( nes ) = zero
              ELSE
                error = .TRUE.
                icount = 1
              END IF
            END IF
          END IF
          CYCLE
        END IF

!  row i of A is non-empty

        k2 = k1 + ROW_count( i ) - 1
        ROW_count( i ) = 0

!  Scan row i of A(trans) to find which rows of A contribute to row i of S

        DO kk = k1, k2
          k = S_row( kk )
          IF ( reals ) amult = S_val( kk ) / DIAG_X( k )
          j1 = Abycol_ptr( k )
          j2 = j1 + COL_count( k ) - 1

!  The column indices for the entries in each row are in order

          Abycol_ptr( k ) = j1 + 1
          COL_count( k ) = j2 - j1

!  Scan row k of A and add multiple of it to row i of S

          DO jj = j1, j2
            j = Abycol_row( jj )
            IF ( error ) THEN

!  Only count the size of arrays needed because of the error

              IF ( S_colptr( j ) < 0 ) CYCLE
              ROW_count( i ) = ROW_count( i ) + 1
              ipos = ROW_count( i )
              S_colptr( j ) = - ROW_count( i )
              IF ( ipos + idist > kk ) THEN
               icount = icount + 1
               idist = idist - 1
              END IF
  
              S_row( ipos ) = j
              CYCLE

            END IF

            IF ( S_colptr( j ) >= 0 ) THEN

!  First contribution to row i of S

              ROW_count( i ) = ROW_count( i ) + 1

!  S_colptr is used to obtain offset for column positions in present row

              S_colptr( j ) = - ROW_count( i )
              ipos = nes + ROW_count( i )
              IF ( ipos > kk ) THEN

!  Insufficient room in S to continue calculation, so only calculate
!  the size required

                error = .TRUE.

!  Overwrite the beginning of S_row since it is not now needed.

                ij = ROW_count( i )
                S_row( 1 : ij - 1 ) = S_row( nes + 1 : nes + ij - 1 )
                nes = 0

!  idist holds the difference between the current position in S_row
!  and the position that it would be if S_row was long enough

                idist = kk - ij

!  icount is the amount by which S_row is too short

                icount = 1
                ipos = ij
                S_row( ipos ) = j
                CYCLE

              END IF

              S_row( ipos ) = j
              IF ( reals ) S_val( ipos ) = Abycol_val( jj ) * amult
              CYCLE

            END IF
            ipos = nes - S_colptr( j )
            IF ( reals ) S_val( ipos ) = S_val( ipos ) + Abycol_val( jj) * amult
          END DO
        END DO

!  Reset S_colptr

        IF ( error ) nes = 0
        j1 = nes + 1
        nes = nes + ROW_count( i )
        S_colptr( S_row( j1 : nes ) ) = 0

        IF ( error ) THEN
         idist = idist + ROW_count( i )
         nes = 0
        END IF

        k1 = k2 + 1
      END DO

!  Set pointer arrays for the start of each column of A and S

      IF ( .NOT. error ) THEN
        S_colptr( 1 ) = 1
        DO i = 1, m
          S_colptr( i + 1 ) = S_colptr( i ) + ROW_count( i )
        END DO
      END IF

      DO j = n, 2, - 1
        Abycol_ptr( j ) = Abycol_ptr( j - 1 )
      END DO
      Abycol_ptr( 1 ) = 1

      IF ( icount > 0 ) THEN
        nes = icount + ls
        ierr = 4
        IF ( control_error > 0 .AND. control_print_level >= 2 )                &
          WRITE( control_error, 2040 ) ls, nes
      END IF

!  Find the diagonal entries in columns dims%c_l_start to dims%c_u_end

      IF ( reals ) THEN
        DO j = dims%c_l_start, dims%c_u_end
          i = S_colptr( j )
          IF ( j == S_row( i ) ) THEN
            S_val( i ) = S_val( i ) +                                          &
                           ( SCALE_C( j ) / DIAG_C( j ) ) * SCALE_C( j )
          ELSE
            WRITE( 6, "( ' j = ', i6, ' i = ', /, ( 10i6 ) )" ) &
              j, S_row( S_colptr( j ) : S_colptr( j + 1 ) - 1 )
            WRITE( 6, "( ' This should be impossible: '/,                      &
           &             ' please report to n.gould@rl.ac.uk' )" )
            STOP
          END IF
        END DO
      END IF

  900 CONTINUE
      RETURN

!  Non-executable statements

 2010 FORMAT( ' ** Error return from LSQP_form_Schur_complement', /,           &
              '    Value of n (number of rows of A) is set to ', I10, /,       &
              '    but n must be at least 1 ' )
 2020 FORMAT( ' ** Error return from LSQP_form_Schur_complement', /,           &
              '    Value of m (number of columns of A) is set to ', I10, /,    &
              '    but m must be at least 1 ' )
 2030 FORMAT( ' ** Error return from LSQP_form_Schur_complement', /,           &
              '    Value of ls = ', I10, ' too small. ', /,                    &
              '    ls must be at least', I10 )
 2040 FORMAT( ' ** Error return from LSQP_form_Schur_complement', /,           &
              '    Arrays S_row and S are too small ', /, '    ls is ', I10, /,&
              '    but must be at least ', I10 )

!  End of LSQP_form_Schur_complement

      END SUBROUTINE LSQP_form_Schur_complement

!-*-*-*-*-*-*-   L S Q P _ A _ B Y _ C O L S   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE LSQP_A_by_cols( n, m, A_ne, A_val, A_col, A_ptr, B_val,       &
                                 B_row, B_colptr )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Takes a matrix A stored by co-ordinates, and returns the same matrix
!  as B stored by columns

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( IN ) :: n, m, A_ne
      INTEGER, INTENT( IN ), DIMENSION( A_ne ) :: A_col
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      INTEGER, INTENT( OUT ), DIMENSION( A_ne ) :: B_row
      INTEGER, INTENT( OUT ), DIMENSION( n + 1 ) :: B_colptr
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( A_ne ) :: A_val
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( A_ne ) :: B_val

!  Local variables

      INTEGER :: i, j, k, l

!  count the number of nonzeros in each column

      B_colptr( : n ) = 0
      DO l = 1, A_ne
        B_colptr( A_col( l ) ) = B_colptr( A_col( l ) ) + 1
      END DO

!  set the starting addresses for each column

      j = 1
      DO i = 1, n
        l = j
        j = j + B_colptr( i )
        B_colptr( i ) = l
      END DO

!  move the entries from A to B

      DO i = 1, m
        DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
          j = A_col( l )
          k = B_colptr( j )
          B_val( k ) = A_val( l )
          B_row( k ) = i
          B_colptr( j ) = B_colptr( j ) + 1
        END DO
      END DO

!  reset the starting addresses

      DO i = n, 1, - 1
        B_colptr( i + 1 ) = B_colptr( i )
      END DO
      B_colptr( 1 ) = 1

      RETURN

!  End of LSQP_A_by_cols

      END SUBROUTINE LSQP_A_by_cols

!-*-*-*-*-*-*-*-   L S Q P _ A N A L Y S E  S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE LSQP_analyse( dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C, &
                               factor, nnzks, lk, liw, ldiag_x, ldiag_c_l,     &
                               ldiag_c_u, IW, Abycol_val, Abycol_row,          &
                               Abycol_ptr, K_colptr, DIAG_X, DIAG_C,           &
                               K, FACTORS, CNTL, print_level, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Analyse the sparsity pattern of the block matrix

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, A_ne, print_level
      INTEGER, INTENT( INOUT ) :: factor
      INTEGER, INTENT( OUT ) :: nnzks, lk, liw, ldiag_x, ldiag_c_l, ldiag_c_u
      INTEGER, INTENT( IN ), DIMENSION( A_ne ) :: A_col
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( A_ne ) :: A_val
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( dims%c_l_start : dims%c_u_end ) :: SCALE_C
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: IW, K_colptr, Abycol_ptr, Abycol_row
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIAG_X, DIAG_C, Abycol_val
      TYPE ( SMT_type ), INTENT( INOUT ) :: K
      TYPE ( SILS_factors ), INTENT( INOUT ) :: FACTORS
      TYPE ( SILS_control ), INTENT( INOUT ) :: CNTL
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      TYPE ( LSQP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: i, ierr, l, max_col, max_len, y_ii
      LOGICAL :: printi, printt, printe, reallocate, get_factors
      TYPE ( SILS_ainfo ) :: AINFO
!      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ptr_K_row, ptr_IW_n, ptr_IW_m
!      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: K%val

      inform%status = GALAHAD_ok

!  Set useful values

      printe = control%error > 0 .AND. print_level >= 1
      printt = control%out > 0 .AND. print_level >= 2 
      printi = control%out > 0 .AND. print_level >= 1 
      max_col = control%max_col
      IF ( max_col < 0 ) max_col = n

!  Print a header indicating the method selected

   10 CONTINUE

      IF ( printi ) THEN
        IF ( factor == 0 .OR. factor == 1 ) THEN
          WRITE( control%out, "( /, '  Schur-complement method used ' )" )
        ELSE
          WRITE( control%out, "( /, '  Augmented system method used ' )" )
        END IF
      END IF

      get_factors = .TRUE.

!  For the Schur complement matrix
!  ===============================

      IF ( factor == 0 .OR. factor == 1 ) THEN

!  Check to see if there are any free variables

        IF ( dims%x_free /= 0 ) THEN

!  The i-th variable is free. Abandon the Schur-complement factorization
!  in favour of the augmented system

          IF ( printi ) WRITE( control%out, 2020 )
          factor = 2
          GO TO 10
        END IF

!  Allocate the arrays for the analysis phase

        liw = MAX( n + 1, n + m )
        reallocate = .TRUE.
        IF ( ALLOCATED( IW ) ) THEN
          IF ( SIZE( IW ) < liw ) THEN ; DEALLOCATE( IW )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( IW( liw ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%IW' ; GO TO 900 ; END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( Abycol_val ) ) THEN
          IF ( SIZE( Abycol_val ) /= A_ne ) THEN 
            DEALLOCATE( Abycol_val )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( Abycol_val( A_ne ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%Abycol_val' ; GO TO 900
          END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( Abycol_row ) ) THEN
          IF ( SIZE( Abycol_row ) /= A_ne ) THEN 
            DEALLOCATE( Abycol_row )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( Abycol_row( A_ne ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%Abycol_row' ; GO TO 900
          END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( Abycol_ptr ) ) THEN
          IF ( SIZE( Abycol_ptr ) /= n + 1 ) THEN 
            DEALLOCATE( Abycol_ptr )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( Abycol_ptr( n + 1 ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%Abycol_ptr' ; GO TO 900
          END IF
        END IF

!  Reorder A so that its entries are ordered by columns

        CALL LSQP_A_by_cols( n, m, A_ne, A_val, A_col, A_ptr,                  &
                             Abycol_val, Abycol_row, Abycol_ptr )

!  Compute the length of the largest column as well as the average length

        max_len = MAXVAL( Abycol_ptr( 2 : ) - Abycol_ptr( : n ) )

        IF ( printi ) WRITE( control%out, 2000 )                               &
          max_len, float( ( Abycol_ptr( n + 1 ) - 1 ) ) / float( n ),          &
          max_col, COUNT( Abycol_ptr( 2 : ) -                                  &
                         Abycol_ptr( : n ) > max_col )

!  Check that the largest column is not too long

        IF ( factor == 0 .AND. max_len > max_col .AND. m > max_sc ) THEN
          IF ( printi ) WRITE( control%out, 2030 ) max_col
          factor = 2
          DEALLOCATE( Abycol_val, Abycol_row, Abycol_ptr )
          GO TO 10
        END IF

!  Continue allocating the arrays for the analysis phase

        ldiag_x = n
        reallocate = .TRUE.
        IF ( ALLOCATED( DIAG_X ) ) THEN
          IF ( SIZE( DIAG_X ) /= ldiag_x ) THEN ; DEALLOCATE( DIAG_X )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( DIAG_X( ldiag_x ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DIAG_X' ; GO TO 900
          END IF
        END IF

        ldiag_c_l = dims%c_l_start
        ldiag_c_u = dims%c_u_end
        reallocate = .TRUE.
        IF ( ALLOCATED( DIAG_C ) ) THEN
          IF ( LBOUND( DIAG_C, 1 ) /= ldiag_c_l .OR.                           &
               UBOUND( DIAG_C, 1 ) /= ldiag_c_u ) THEN ; DEALLOCATE( DIAG_C )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( DIAG_C( ldiag_c_l : ldiag_c_u ),                           &
                    STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DIAG_C' ; GO TO 900
          END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( K_colptr ) ) THEN
          IF ( SIZE( K_colptr ) /= m + 1 ) THEN ; DEALLOCATE( K_colptr )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( K_colptr( m + 1 ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%K_colptr' ; GO TO 900
          END IF
        END IF

        lk = max( A_ne + m, 2 * A_ne, control%valmin )
        reallocate = .TRUE.
        IF ( ALLOCATED( K%row ) ) THEN
          IF ( SIZE( K%row ) < lk ) THEN ; DEALLOCATE( K%row )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( K%row( lk ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%K%row' ; GO TO 900 ; END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( K%col ) ) THEN
          IF ( SIZE( K%col ) < lk ) THEN ; DEALLOCATE( K%col )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( K%col( lk ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%K%col' ; GO TO 900 ; END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( K%val ) ) THEN
          IF ( SIZE( K%val ) < lk ) THEN ; DEALLOCATE( K%val )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( K%val( lk ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%K%val' ; GO TO 900 ; END IF
        END IF

!  Permute A so that its entries appear by columns, with the entries in
!  each column sorted by increasing row number. Also form the data
!  structure required to hold K = A DIAG_X(inv) A(transpose) + DIAG_C(inv)

   20   CONTINUE

        IF ( m > 0 ) THEN

!          ptr_K_val => K%val
!          ptr_K_row => K%row
!          ptr_IW_n => IW( : n )
!          ptr_IW_m => IW( n + 1 : )

          CALL LSQP_form_Schur_complement(                                     &
                  dims, n, m, A_ne, Abycol_val, Abycol_row, Abycol_ptr,        &
                  DIAG_X, SCALE_C, DIAG_C, lk, K%val, K%row, K_colptr,         &
                  K%ne, ierr, IW( : n ), IW( n + 1 : ), .FALSE.,               &
                  control%error, print_level )
!         inform%nfacts = inform%nfacts + 1 

!  Check for error returns

          IF ( ierr == 3 ) THEN

!  Insufficient space. Allocate more and retry

            DEALLOCATE( K%row, K%col, K%val )
            lk = 2 * lk
            ALLOCATE( K%row( lk ), K%col( lk ), K%val( lk ),    &
                      STAT = inform%alloc_status )
            IF ( inform%alloc_status /= 0 ) THEN 
              inform%bad_alloc = 'lsqp: data%K' ; GO TO 900 ; END IF
  
            GO TO 20
          END IF
  
          IF ( ierr == 4 ) THEN

!  Insufficient space. Allocate more and retry

            DEALLOCATE( K%row, K%col, K%val )
            lk = K%ne
            ALLOCATE( K%row( lk ), K%col( lk ), K%val( lk ),    &
                      STAT = inform%alloc_status )
            IF ( inform%alloc_status /= 0 ) THEN 
              inform%bad_alloc = 'lsqp: data%K' ; GO TO 900 ; END IF
  
            GO TO 20
          END IF

!  Record the column numbers of each entry in A A^T

          DO i = 1, m
            K%col( K_colptr( i ) : K_colptr( i + 1 ) - 1 ) = i
          END DO
        ELSE
          get_factors = .FALSE.
          K%ne = 0
        END IF

!  Now find a good ordering of the rows and columns of A A^T for
!  the sparse factorization. Set the dimensions of K

        K%n = m
        CNTL%pivoting = 2

!  For the KKT matrix
!  ==================

      ELSE

!  Allocate the arrays for the analysis phase

        ldiag_x = 0
        reallocate = .TRUE.
        IF ( ALLOCATED( DIAG_X ) ) THEN
          IF ( SIZE( DIAG_X ) /= ldiag_x ) THEN ; DEALLOCATE( DIAG_X )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( DIAG_X( ldiag_x ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DIAG_X' ; GO TO 900
          END IF
        END IF

        ldiag_c_l = 0
        ldiag_c_u = 0
        reallocate = .TRUE.
        IF ( ALLOCATED( DIAG_C ) ) THEN
          IF ( LBOUND( DIAG_C, 1 ) /= ldiag_c_l .OR.                           &
               UBOUND( DIAG_C, 1 ) /= ldiag_c_u ) THEN ; DEALLOCATE( DIAG_C )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( DIAG_C( ldiag_c_l : ldiag_c_u ),                           &
                    STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%DIAG_C' ; GO TO 900
          END IF
        END IF

        lk = A_ne + n + 2 * dims%nc
        reallocate = .TRUE.
        IF ( ALLOCATED( K%row ) ) THEN
          IF ( SIZE( K%row ) < lk ) THEN ; DEALLOCATE( K%row )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( K%row( lk ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%K%row' ; GO TO 900 ; END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( K%col ) ) THEN
          IF ( SIZE( K%col ) < lk ) THEN ; DEALLOCATE( K%col )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( K%col( lk ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%K%col' ; GO TO 900 ; END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( K%val ) ) THEN
          IF ( SIZE( K%val ) < lk ) THEN ; DEALLOCATE( K%val )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( K%val( lk ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            inform%bad_alloc = 'lsqp: data%K%val'; GO TO 900
          END IF
        END IF

!  Set the coordinates of A in K

        DO i = 1, m
          y_ii = dims%y_s + i - 1
          DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
            K%row( l ) = y_ii
          END DO
        END DO
        K%col( : A_ne ) = A_col( : A_ne ) 
        K%val( : A_ne ) = A_val( : A_ne )

!  Set the coodinates corresponding to the slack variables

        DO i = 1, dims%nc
          K%row( A_ne + i ) = dims%y_i + i - 1
          K%col( A_ne + i ) = n + i
          K%val( A_ne + i ) = - SCALE_C( dims%c_equality + i )
        END DO

!  Set the coordinates of diagonal in K

        nnzks = A_ne + dims%nc

        DO i = dims%x_s, dims%c_e
          K%row( nnzks + i ) = i ; K%col( nnzks + i ) = i
        END DO

!  Set the dimensions of K

        K%n = dims%v_e ; K%ne = nnzks + dims%c_e
        CNTL%pivoting = 1

      END IF
!     WRITE( control%out, "( ' nnzk ', I10 )" ) K%ne

!  Analyse the sparsity pattern of the matrix

      IF ( get_factors ) THEN
        CALL SILS_analyse( K, FACTORS, CNTL, AINFO )

!  Record the storage requested

        inform%factorization_integer = AINFO%nirnec 
        inform%factorization_real = AINFO%nrlnec

!  Check for error returns

        inform%factorization_status = AINFO%flag
        IF ( AINFO%flag < 0 ) THEN
          IF ( printe ) WRITE( control%error, 2040 ) AINFO%flag, 'SILS_analyse'
          inform%status = GALAHAD_error_analysis ; RETURN
        ELSE IF ( AINFO%flag > 0 ) THEN 
          IF ( printt ) WRITE( control%out, 2050 ) AINFO%flag, 'SILS_analyse'
        END IF
        IF ( printt ) WRITE( control%out,                                      &
          "( ' real/integer space required for factors ', 2I10 )" )            &
            AINFO%nrladu, AINFO%niradu

      ELSE
        inform%factorization_integer = 0 ; inform%factorization_real = 0
      END IF

      RETURN

!  Unsuccessful returns

  900 CONTINUE
      inform%status = GALAHAD_error_allocate
      WRITE( control%out, 2900 ) inform%bad_alloc, inform%alloc_status
      RETURN

!  Non-executable statements

 2000 FORMAT( '  maximum, average column lengths ', I7, 0P, F10.1, /,          &
              '  number of columns longer than   ', I7, ' is', I7 )
 2020 FORMAT( /, ' ** There are free variables - abandon the Schur-complement',&
              /, '    factorization in favour of one of the augmented matrix' )
 2030 FORMAT( ' ** The maximum column length in A is larger than', /,          &
              '    max_col =', I7, ' - abandon the Schur-complement', /,       &
              '    factorization in favour of one of the augmented matrix', / )
 2040 FORMAT( '   **  Error return ', I6, ' from ', A15 ) 
 2050 FORMAT( '   **  Warning ', I6, ' from ', A15 ) 
 2900 FORMAT( ' ** Message from -LSQP_analyse-', /,                            &
              ' Allocation error, for ', A, /, ' status = ', I6 )

      END SUBROUTINE LSQP_analyse

!-*-*-*-*-*-   L S Q P _ I N D I C A T O R S   S U B R O U T I N E   -*-*-*-*-*-

     SUBROUTINE LSQP_indicators( dims, n, m, C_l, C_u, C_last, C,               &
                                 DIST_C_l, DIST_C_u, X_l, X_u, X_last, X,       &
                                 DIST_X_l, DIST_X_u, Y_l, Y_u, Z_l, Z_u,        &
                                 Y_last, Z_last, control, C_stat, B_stat )

!  ---------------------------------------------------------------------------

!  Compute indicatirs for active simnple bounds and general constraints

!  C_stat is an INTEGER array of length m, which if present will be 
!   set on exit to indicate the likely ultimate status of the constraints. 
!   Possible values are 
!   C_stat( i ) < 0, the i-th constraint is likely in the active set, 
!                    on its lower bound, 
!               > 0, the i-th constraint is likely in the active set
!                    on its upper bound, and
!               = 0, the i-th constraint is likely not in the active set

!  B_stat is an INTEGER array of length m, which if present will be 
!   set on exit to indicate the likely ultimate status of the simple bound 
!   constraints. Possible values are 
!   B_stat( i ) < 0, the i-th bound constraint is likely in the active set, 
!                    on its lower bound, 
!               > 0, the i-th bound constraint is likely in the active set
!                    on its upper bound, and
!               = 0, the i-th bound constraint is likely not in the active set

!  ---------------------------------------------------------------------------

!  Dummy arguments

      INTEGER, INTENT( IN ) :: n, m
      TYPE ( LSQP_dims_type ), INTENT( IN ) :: dims
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X_l, X_u, X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X_last, Z_last
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: C_l, C_u
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: C_last, Y_last
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%x_l_start : dims%x_l_end ) :: DIST_X_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%x_u_start : dims%x_u_end ) :: DIST_X_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%c_l_start : dims%c_l_end ) :: DIST_C_l, Y_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%c_u_start : dims%c_u_end ) :: DIST_C_u, Y_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%c_l_start : dims%c_u_end ) :: C
      REAL ( KIND = wp ), INTENT(IN ),                                        &
             DIMENSION( dims%x_free + 1 : dims%x_l_end ) ::  Z_l
      REAL ( KIND = wp ), INTENT( IN ),                                       &
             DIMENSION( dims%x_u_start : n ) :: Z_u
      TYPE ( LSQP_control_type ), INTENT( IN ) :: control        
      INTEGER, INTENT( INOUT ), OPTIONAL, DIMENSION( m ) :: C_stat
      INTEGER, INTENT( INOUT ), OPTIONAL, DIMENSION( n ) :: B_stat

!  Local variables

      INTEGER :: i

!     IF ( printd ) WRITE(  control%out,                                     &
!       "( /, ' Constraints : ', /, '                   ',                   &
!    &   '        <------ Bounds ------> ', /                                &
!    &   '      # name       state      Lower       Upper     Multiplier' )" )
!     DO i = dims%c_equality + 1, m 
!       IF ( printd ) WRITE(  control%out,"( 2I7, 4ES12.4 )" ) i,            &
!         C_stat( i ), C( i ), C_l( i ), C_u( i ), Y( i )
!     END DO 

!     IF ( printd ) WRITE(  control%out,                                     &
!        "( /, ' Solution : ', /,'                    ',                     &
!       &    '        <------ Bounds ------> ', /                            &
!       &    '      # name       state      Lower       Upper       Dual' )" )
!     DO i = dims%x_free + 1, n 
!       IF ( printd ) WRITE(  control%out,"( 2I7, 4ES12.4 )" ) i,            &
!         B_stat( i ), X( i ), X_l( i ), X_u( i ), Z( i )
!     END DO 

!  equality constraints

      C_stat( : dims%c_equality ) = - 1

!  free variables

      B_stat( : dims%x_free ) = 0

!  Compute the required indicator

!  ----------------------------------
!  Type 1 ("primal") indicator used: 
!  ----------------------------------

!    a variable/constraint will be "inactive" if
!        distance to nearest bound > indicator_p_tol
!    for some constant indicator_p_tol close-ish to zero

      IF ( control%indicator_type == 1 ) THEN
        DO i = dims%c_equality + 1, m 
          IF ( ABS( C( i ) - C_l( i ) ) < control%indicator_tol_p ) THEN
            C_stat( i ) = - 1
          ELSE IF ( ABS( C( i ) - C_u( i ) ) < control%indicator_tol_p ) THEN
            C_stat( i ) = 1
          ELSE
            C_stat( i ) = 0
          END IF
        END DO
        DO i = dims%x_free + 1, n 
          IF ( ABS( X( i ) - X_l( i ) ) < control%indicator_tol_p ) THEN
            B_stat( i ) = - 1
          ELSE IF ( ABS( X( i ) - X_u( i ) ) < control%indicator_tol_p ) THEN
            B_stat( i ) = 1
          ELSE
            B_stat( i ) = 0
          END IF
        END DO 

!  --------------------------------------
!  Type 2 ("primal-dual") indicator used: 
!  --------------------------------------

!    a variable/constraint will be "inactive" if
!        distance to nearest bound 
!          > indicator_tol_pd * size of corresponding multiplier
!    for some constant indicator_tol_pd close-ish to one.

      ELSE IF ( control%indicator_type == 2 ) THEN

!  constraints with lower bounds

        DO i = dims%c_l_start, dims%c_u_start - 1
          IF ( DIST_C_l( i ) > control%indicator_tol_pd * Y_l( i ) ) THEN
            C_stat( i ) = 0
          ELSE
            C_stat( i ) = - 1
          END IF
        END DO

!  constraints with both lower and upper bounds

        DO i = dims%c_u_start, dims%c_l_end
          IF ( DIST_C_l( i ) <= DIST_C_u( i ) ) THEN
            IF ( DIST_C_l( i ) > control%indicator_tol_pd * Y_l( i ) ) THEN
              C_stat( i ) = 0
            ELSE
              C_stat( i ) = - 1
            END IF
          ELSE
            IF ( DIST_C_u( i ) > - control%indicator_tol_pd * Y_u( i ) ) THEN
              C_stat( i ) = 0
            ELSE
              C_stat( i ) = 1
            END IF
          END IF
        END DO

!  constraints with upper bounds

        DO i = dims%c_l_end + 1, m
          IF ( DIST_C_u( i ) > - control%indicator_tol_pd * Y_u( i ) ) THEN
            C_stat( i ) = 0
          ELSE
            C_stat( i ) = 1
          END IF
        END DO

!  simple non-negativity

        DO i = dims%x_free + 1, dims%x_l_start - 1
          IF ( X( i ) > control%indicator_tol_pd * Z_l( i ) ) THEN
            B_stat( i ) = 0
          ELSE
            B_stat( i ) = - 1
          END IF
        END DO

!  simple bound from below

        DO i = dims%x_l_start, dims%x_u_start - 1
          IF ( DIST_X_l( i ) > control%indicator_tol_pd * Z_l( i ) ) THEN
            B_stat( i ) = 0
          ELSE
            B_stat( i ) = - 1
          END IF
        END DO

!  simple bound from below and above

        DO i = dims%x_u_start, dims%x_l_end
          IF ( DIST_X_l( i ) <= DIST_X_u( i ) ) THEN
            IF ( DIST_X_l( i ) > control%indicator_tol_pd * Z_l( i ) ) THEN
              B_stat( i ) = 0
            ELSE
              B_stat( i ) = - 1
            END IF
          ELSE
            IF ( DIST_x_u( i ) > - control%indicator_tol_pd * Z_u( i ) ) THEN
              B_stat( i ) = 0
            ELSE
              B_stat( i ) = 1
            END IF
          END IF
        END DO

!  simple bound from above

        DO i = dims%x_l_end + 1, dims%x_u_end
          IF ( DIST_x_u( i ) > - control%indicator_tol_pd * Z_u( i ) ) THEN
            B_stat( i ) = 0
          ELSE
            B_stat( i ) = 1
          END IF
        END DO

!  simple non-positivity

        DO i = dims%x_u_end + 1, n
          IF ( - X( i ) > - control%indicator_tol_pd * Z_u( i ) ) THEN
            B_stat( i ) = 0
          ELSE
            B_stat( i ) = 1
          END IF
        END DO

!  --------------------------------
!  Type 3 ("Tapia") indicator used: 
!  --------------------------------

!    a variable/constraint will be "inactive" if
!        distance to nearest bound now 
!          > indicator_tol_tapia * distance to same bound at previous iteration
!    for some constant indicator_tol_tapia close-ish to one.

      ELSE IF ( control%indicator_type == 3 ) THEN

!  constraints with lower bounds

        DO i = dims%c_l_start, dims%c_u_start - 1
!WRITE( 6, "( 'i,dc,dc-,ratio', I6,3ES12.4)" )  i, C( i ) - C_l( i ), &
!C_last( i ) - C_l( i ) , ( C( i ) - C_l( i ) ) / (  C_last( i ) - C_l( i ) )
          IF ( ABS( C( i ) - C_l( i ) ) > control%indicator_tol_tapia *        &
               ABS( C_last( i ) - C_l( i ) ) ) THEN
            C_stat( i ) = 0
          ELSE
            IF ( ABS( Y_l( i ) / Y_last( i ) )                                 &
                 > control%indicator_tol_tapia ) THEN
              C_stat( i ) = - 1
            ELSE
!             write(6,*) i, ABS( Y_l( i ) / Y_last( i ) )
              C_stat( i ) = - 2
            END IF
          END IF
        END DO

!  constraints with both lower and upper bounds

        DO i = dims%c_u_start, dims%c_l_end
          IF ( DIST_C_l( i ) <= DIST_C_u( i ) ) THEN
            IF ( ABS( C( i ) - C_l( i ) ) > control%indicator_tol_tapia *      &
                 ABS( C_last( i ) - C_l( i ) ) ) THEN
              C_stat( i ) = 0
            ELSE
              IF ( ABS( Y_l( i ) / Y_last( i ) )                               &
                   > control%indicator_tol_tapia ) THEN
                C_stat( i ) = - 1
              ELSE
!               write(6,*) i, ABS( Y_l( i ) / Y_last( i ) )
                C_stat( i ) = - 2
              END IF
            END IF
          ELSE
            IF ( ABS( C( i ) - C_u( i ) ) > control%indicator_tol_tapia *      &
                 ABS( C_last( i ) - C_u( i ) ) )  THEN
              C_stat( i ) = 0
            ELSE
              IF ( ABS( Y_u( i ) / Y_last( i ) )                               &
                   > control%indicator_tol_tapia ) THEN
                C_stat( i ) = 1
              ELSE
!               write(6,*) i, ABS( Y_u( i ) / Y_last( i ) )
                C_stat( i ) = 2
              END IF
            END IF
          END IF
        END DO

!  constraints with upper bounds

        DO i = dims%c_l_end + 1, m
!WRITE( 6, "( 'i,dc,dc-,ratio', I6,3ES12.4)" )  i, C_u( i ) - C( i ), &
!C_u( i ) - C_last( i ), ( C_u( i ) - C( i ) ) / ( C_u( i ) - C_last( i ) )
          IF ( ABS( C( i ) - C_u( i ) ) > control%indicator_tol_tapia *        &
               ABS( C_last( i ) - C_u( i ) ) )  THEN
            C_stat( i ) = 0
          ELSE
            IF ( ABS( Y_u( i ) / Y_last( i ) )                                 &
                 > control%indicator_tol_tapia ) THEN
              C_stat( i ) = 1
            ELSE
!             write(6,*) i, ABS( Y_u( i ) / Y_last( i ) )
              C_stat( i ) = 2
            END IF
          END IF
        END DO

!  simple non-negativity

        DO i = dims%x_free + 1, dims%x_l_start - 1
          IF ( ABS( X( i ) ) > control%indicator_tol_tapia *                   &
               ABS( X_last( i ) ) ) THEN
            B_stat( i ) = 0
          ELSE
            IF ( ABS( Z_l( i ) / Z_last( i ) )                                 &
                 > control%indicator_tol_tapia ) THEN
              B_stat( i ) = - 1
            ELSE
!             write(6,*) i, ABS( Z_l( i ) / Z_last( i ) )
              B_stat( i ) = - 2
            END IF
          END IF
        END DO

!  simple bound from below

        DO i = dims%x_l_start, dims%x_u_start - 1
          IF ( ABS( X( i ) - X_l( i ) ) > control%indicator_tol_tapia *        &
               ABS( X_last( i ) - X_l( i ) ) ) THEN
            B_stat( i ) = 0
          ELSE
            IF ( ABS( Z_l( i ) / Z_last( i ) )                                 &
                 > control%indicator_tol_tapia ) THEN
              B_stat( i ) = - 1
            ELSE
!             write(6,*) i, ABS( Z_l( i ) / Z_last( i ) )
              B_stat( i ) = - 2
            END IF
          END IF
        END DO

!  simple bound from below and above

        DO i = dims%x_u_start, dims%x_l_end
          IF ( DIST_X_l( i ) <= DIST_X_u( i ) ) THEN
            IF ( ABS( X( i ) - X_l( i ) ) > control%indicator_tol_tapia *      &
                 ABS( X_last( i ) - X_l( i ) ) ) THEN
              B_stat( i ) = 0
            ELSE
              IF ( ABS( Z_l( i ) / Z_last( i ) )                               &
                   > control%indicator_tol_tapia ) THEN
                B_stat( i ) = - 1
              ELSE
!               write(6,*) i, ABS( Z_l( i ) / Z_last( i ) )
                B_stat( i ) = - 2
              END IF
            END IF
          ELSE
            IF ( ABS( X( i ) - X_u( i ) ) > control%indicator_tol_tapia *      &
                 ABS( X_last( i ) - X_u( i ) ) ) THEN
              B_stat( i ) = 0
            ELSE
              IF ( ABS( Z_u( i ) / Z_last( i ) )                               &
                   > control%indicator_tol_tapia ) THEN
                B_stat( i ) = 1
              ELSE
!               write(6,*) i, ABS( Z_u( i ) / Z_last( i ) )
                B_stat( i ) = 2
              END IF
            END IF
          END IF
        END DO

!  simple bound from above

        DO i = dims%x_l_end + 1, dims%x_u_end
            IF ( ABS( X( i ) - X_u( i ) ) > control%indicator_tol_tapia *      &
                 ABS( X_last( i ) - X_u( i ) ) ) THEN
            B_stat( i ) = 0
          ELSE
            IF ( ABS( Z_u( i ) / Z_last( i ) )                                 &
                 > control%indicator_tol_tapia ) THEN
              B_stat( i ) = 1
            ELSE
!             write(6,*) i, ABS( Z_u( i ) / Z_last( i ) )
              B_stat( i ) = 2
            END IF
          END IF
        END DO

!  simple non-positivity

        DO i = dims%x_u_end + 1, n
          IF ( ABS( X( i ) ) > control%indicator_tol_tapia *                   &
               ABS( X_last( i ) ) ) THEN
            B_stat( i ) = 0
          ELSE
            IF ( ABS( Z_u( i ) / Z_last( i ) )                                 &
                 > control%indicator_tol_tapia ) THEN
              B_stat( i ) = 1
            ELSE
!             write(6,*) i, ABS( Z_u( i ) / Z_last( i ) )
              B_stat( i ) = 2
            END IF
          END IF
        END DO
      ELSE
      END IF

!  End of LSQP_indicators

      END SUBROUTINE LSQP_indicators

!  End of module LSQP

   END MODULE GALAHAD_LSQP_double



