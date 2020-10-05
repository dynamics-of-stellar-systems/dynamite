! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-*-*-  G A L A H A D _ W C P    M O D U L E  -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   started life as part of GALAHAD_CLS ~2000, which mutated into part of
!   GALAHAD LSQP, released pre GALAHAD Version 1.0. April 10th 2001.
!   update extracted and released with GALAHAD Version 2.0.November 1st 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_WCP_double

!     --------------------------------------------------
!     |                                                |
!     | Find a well-centered point within the polytope |
!     |                                                |
!     |          c_l <= A x <= c_u                     |
!     |          x_l <=  x <= x_u                      |
!     |                                                |
!     | using an infeasible-point primal-dual method   |
!     |                                                |
!     --------------------------------------------------

!NOT95USE GALAHAD_CPU_time
      USE GALAHAD_SPACE_double
      USE GALAHAD_SILS_double
      USE GALAHAD_QPT_double
      USE GALAHAD_SPECFILE_double
      USE GALAHAD_QPP_double, WCP_dims_type => QPP_dims_type
      USE GALAHAD_QPD_double, WCP_data_type => QPD_data_type,                  &
                              WCP_AX => QPD_AX

      USE GALAHAD_SORT_double, ONLY: SORT_heapsort_build,                      &
         SORT_heapsort_smallest, SORT_inverse_permute
      USE GALAHAD_ROOTS_double
      USE GALAHAD_FDC_double
      USE GALAHAD_STRING_double

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: WCP_initialize, WCP_read_specfile, WCP_solve,                  &
                WCP_terminate, QPT_problem_type, SMT_type, SMT_put, SMT_get,   &
                WCP_form_Schur_complement, WCP_A_by_cols, WCP_Ax,              &
                WCP_data_type, WCP_dims_type

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

      TYPE, PUBLIC :: WCP_time_type
        REAL :: total, preprocess, find_dependent, analyse, factorize, solve
      END TYPE

      TYPE, PUBLIC :: WCP_control_type
        INTEGER :: error, out, print_level, start_print, stop_print, maxit
        INTEGER :: initial_point
        INTEGER :: factor, max_col, indmin, valmin, itref_max, infeas_max
        INTEGER :: perturbation_strategy, restore_problem
        REAL ( KIND = wp ) :: infinity, stop_p, stop_d, stop_c, prfeas, dufeas
        REAL ( KIND = wp ) :: mu_target, required_infeas_reduction, implicit_tol
        REAL ( KIND = wp ) :: pivot_tol, pivot_tol_for_dependencies, zero_pivot
        REAL ( KIND = wp ) :: perturb_start, alpha_scale, identical_bounds_tol
        REAL ( KIND = wp ) :: reduce_perturb_factor, reduce_perturb_multiplier
        REAL ( KIND = wp ) :: mu_accept_fraction, mu_increase_factor
        REAL ( KIND = wp ) :: insufficiently_feasible, perturbation_small
        REAL ( KIND = wp ) :: cpu_time_limit
        LOGICAL :: remove_dependencies, treat_zero_bounds_as_general
        LOGICAL :: just_feasible, balance_initial_complementarity
        LOGICAL :: use_corrector
        LOGICAL :: space_critical, deallocate_error_fatal
        LOGICAL :: record_x_status, record_c_status
        CHARACTER ( LEN = 30 ) :: prefix
      END TYPE

      TYPE, PUBLIC :: WCP_inform_type
        INTEGER :: status, alloc_status, iter, factorization_status
        INTEGER :: factorization_integer, factorization_real, nfacts
        INTEGER :: c_implicit, x_implicit, y_implicit, z_implicit
        REAL ( KIND = wp ) :: obj, non_negligible_pivot
        LOGICAL :: feasible
        CHARACTER ( LEN = 80 ) :: bad_alloc
        TYPE ( WCP_time_type ) :: time
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: X_status
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: C_status
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
      REAL ( KIND = wp ), PARAMETER :: tenm1 = ten ** ( - 1 )
      REAL ( KIND = wp ), PARAMETER :: tenm2 = ten ** ( - 2 )
      REAL ( KIND = wp ), PARAMETER :: tenm3 = ten ** ( - 3 )
      REAL ( KIND = wp ), PARAMETER :: tenm4 = ten ** ( - 4 )
      REAL ( KIND = wp ), PARAMETER :: tenm5 = ten ** ( - 5 )
      REAL ( KIND = wp ), PARAMETER :: tenm7 = ten ** ( - 7 )
      REAL ( KIND = wp ), PARAMETER :: tenm8 = ten ** ( - 8 )
      REAL ( KIND = wp ), PARAMETER :: tenm10 = ten ** ( - 10 )
      REAL ( KIND = wp ), PARAMETER :: ten2 = ten ** 2
      REAL ( KIND = wp ), PARAMETER :: ten3 = ten ** 3
      REAL ( KIND = wp ), PARAMETER :: ten4 = ten ** 4
      REAL ( KIND = wp ), PARAMETER :: ten5 = ten ** 5
      REAL ( KIND = wp ), PARAMETER :: infinity = HUGE( one )
      REAL ( KIND = wp ), PARAMETER :: epsmch = EPSILON( one )
      REAL ( KIND = wp ), PARAMETER :: mu_tol = ten ** 2

   CONTAINS

!-*-*-*-*-*-   W C P _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE WCP_initialize( data, control )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Default control data for WCP. This routine should be called before
!  WCP_primal_dual
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
!   initial_point. How to choose the initial point.
!     Possible values are
!
!      0  the values input in X, shifted to be at least prfeas from
!         their nearest bound, will be used
!      1  the nearest point to the "bound average" 0.5(X_l+X_u) that satisfies
!          the linear constraints will be used
!
!   factor. The factorization to be used.
!     Possible values are
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
!     of the problem is not reduced by at least a factor 
!     control%required_infeas_reduction before the problem is flagged 
!     as infeasible (see required_infeas_reduction)
!
!   perturbation_strategy. The strategy used to reduce relaxed constraint bounds.
!     Possible values are
!
!      0 do not perturb the constraints
!      1 reduce all perturbations by the same amount with linear reduction
!      2 reduce each perturbation as much as possible with linear reduction
!      3 reduce all perturbations by the same amount with superlinear reduction
!      4 reduce each perturbation as much as possible with superlinear reduction

!   restore_problem. Indicates whether and how much of the input problem
!    should be restored on output. Possible values are
!
!      0 nothing restored
!      1 scalar and vector parameters
!      2 all parameters
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
!   mu_target. The target value of the barrier parameter. If mu_target is
!    not positive, it will be reset to an appropriate value
!           
!   mu_accept_fraction. The complemtary slackness x_i.z_i will be judged to
!     lie within an acceptable margin around its target value mu as soon
!     as mu_accept_fraction * mu <= x_i.z_i <= ( 1 / mu_accept_fraction ) * mu;
!     the perturbations will be reduced as soon as all of the complemtary 
!     slacknesses x_i.z_i lie within acceptable bounds. mu_accept_fraction
!     will be reset to ensure that it lies in the interval (0,1]
!
!   mu_increase_factor. The target value of the barrier parameter will be 
!    increased by mu_increase_factor for infeasible constraints every time 
!    the perturbations are adjusted
!
!   required_infeas_reduction. If the overall infeasibility of the problem is 
!    not reduced by at least a factor required_infeas_reduction over 
!    control%infeas_max iterations, the problem is flagged as infeasible 
!    (see infeas_max)
!
!   alpha_scale. The test for rank defficiency will be to factorize
!     ( alpha_scale I  A^T ; A  0 )
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
!   perturb_start. The constraint bounds will initially be relaxed by 
!     perturb_start. This perturbation will subsequently be reduced to zero. 
!     If perturb_start < 0, the amount by which the bounds are relaxed will 
!     be computed automatically
!
!   reduce_perturb_factor, reduce_perturb_multiplier and insufficiently_feasible.
!    The constraint perturbation will be reduced as follows:
!
!    * if the variable lies outside a bound, the corresponding perturbation
!      will be reduced to 
!        reduce_perturb_factor * current pertubation 
!           + ( 1 - reduce_perturb_factor ) * violation
!    * otherwise, if the variable lies within insufficiently_feasible of its 
!      bound the pertubation will be reduced to
!        reduce_perturb_multiplier * current pertubation
!    * otherwise if will be set to zero
!
!   implicit_tol. Any primal or dual variable that is less feasible than
!    implicit_tol will be regarded as defining an implicit constraint
!
!   perturbation_small. If the maximum constraint pertubation is smaller than
!     perturbation_small and the violation is smaller than implicit_tol, the
!     method will deduce that there is a feasible point but no interior.
!
!   identical_bounds_tol. Any pair of constraint bounds (c_l,c_u) or (x_l,x_u)
!    that are closer than identical_bounds_tol will be reset to the average
!    of their values
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
!   just_feasible. If just_feasible is .TRUE., the algorithm will stop as
!    soon as a feasible interior point is found. Otherwise, a well-centered
!    interior point will be sought
!
!   balance_initial_complementarity is .true. if the initial complemetarity
!    is required to be balanced
!
!   ues_corrector is .true. if a corrector step is to be combined with the
!    (predictor) search direction
!
!   space_critical. If true, every effort will be made to use as little
!     space as possible. This may result in longer computation times
!
!   deallocate_error_fatal. If true, any array/pointer deallocation error
!     will terminate execution. Otherwise, computation will continue
!
!   record_x_status. If true, the array inform%X_status will be allocated
!     and the status of the bound constraints will be reported on exit.
!     In this case, possible values of inform%X_status(i) are as follows:
!        0  the variable lies between its bounds
!       -1  the variable lies on its lower bound for all feasible points
!        1  the variable lies on its upper bound for all feasible points
!       -2  the variable never lies on its lower bound at any feasible point
!        2  the variable never lies on its upper bound at any feasible point
!        3  the bounds are equal, and the variable takes this value for 
!           all feasible points
!       -3  the variable never lies on either bound at any feasible point
!
!   record_c_status. If true, the array inform%C_status will be allocated
!     and the status of the general constraints will be reported on exit.
!     In this case, possible values of inform%C_status(i) are as follows:
!        0  the constraint lies between its bounds
!       -1  the constraint lies on its lower bound for all feasible points
!           and may be fixed at this value and removed from the problem 
!        1  the constraint lies on its upper bound for all feasible points
!           and may be fixed at this value and removed from the problem 
!       -2  the constraint never lies on its lower bound at any feasible point
!           and the bound may be removed from the problem 
!        2  the constraint never lies on its upper bound at any feasible point
!           and the bound may be removed from the problem 
!        3  the bounds are equal, and the constraint takes this value for 
!           all feasible points
!       -3  the constraint never lies on either bound at any feasible point
!           and the constraint may be removed from the problem 
!        4  the constraint is implied by the others and may be removed
!           from the problem
!
!  CHARACTER control parameters:
!
!  prefix (len=30). All output lines will be prefixed by 
!    %prefix(2:LEN(TRIM(%prefix))-1)
!   where %prefix contains the required string enclosed in 
!   quotes, e.g. "string" or 'string'
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      TYPE ( WCP_data_type ), INTENT( OUT ) :: data
      TYPE ( WCP_control_type ), INTENT( OUT ) :: control        

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
      control%initial_point = 0
      control%factor = 0
      control%max_col = 35
      control%indmin = 10000
      control%valmin = 10000
      control%itref_max = 1
      control%infeas_max = 200
      control%restore_problem = 2
      control%perturbation_strategy = 2

!  Real parameters

      control%infinity = ten ** 19
      control%stop_p = epsmch ** 0.33
      control%stop_c = epsmch ** 0.33
      control%stop_d = epsmch ** 0.33
      control%perturb_start = - one
      control%implicit_tol = epsmch ** 0.33
      control%prfeas = one
      control%dufeas = one
      control%mu_target = - one
      control%mu_accept_fraction = one
      control%mu_increase_factor = two
!     control%pivot_tol = data%CNTL%u
      control%pivot_tol = epsmch ** 0.75
      control%pivot_tol_for_dependencies = half
      control%alpha_scale = point01
      control%zero_pivot = epsmch ** 0.75
      control%identical_bounds_tol = epsmch
      control%reduce_perturb_factor = 0.25_wp
      control%reduce_perturb_multiplier = point01
      control%insufficiently_feasible = epsmch ** 0.25
      control%perturbation_small = - one
      control%required_infeas_reduction = one - point01
      control%cpu_time_limit = - one

!  Logical parameters

      control%remove_dependencies = .TRUE.
      control%treat_zero_bounds_as_general = .FALSE.
      control%just_feasible = .FALSE.
      control%balance_initial_complementarity = .FALSE.
      control%use_corrector = .FALSE.
      control%space_critical = .FALSE.
      control%deallocate_error_fatal  = .FALSE.
      control%record_x_status = .TRUE.
      control%record_c_status = .TRUE.

!  Character parameters

      control%prefix = '""                            '

      data%trans = 0 ; data%tried_to_remove_deps = .FALSE.
      data%save_structure = .TRUE.

      RETURN  

!  End of WCP_initialize

      END SUBROUTINE WCP_initialize

!-*-*-*-*-   W C P _ R E A D _ S P E C F I L E  S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE WCP_read_specfile( control, device, alt_specname )

!  Reads the content of a specification file, and performs the assignment of 
!  values associated with given keywords to the corresponding control parameters

!  The defauly values as given by WCP_initialize could (roughly) 
!  have been set as:

! BEGIN WCP SPECIFICATIONS (DEFAULT)
!  error-printout-device                             6
!  printout-device                                   6
!  print-level                                       1
!  maximum-number-of-iterations                      1000
!  start-print                                       -1
!  stop-print                                        -1
!  initial-point-used                                1
!  factorization-used                                0
!  maximum-column-nonzeros-in-schur-complement       35
!  initial-integer-workspace                         10000
!  initial-real-workspace                            10000
!  maximum-refinements                               1
!  maximum-poor-iterations-before-infeasible         200
!  perturbation-strategy                             2
!  restore-problem-on-output                         2
!  infinity-value                                    1.0D+10
!  primal-accuracy-required                          1.0D-5
!  dual-accuracy-required                            1.0D-5
!  complementary-slackness-accuracy-required         1.0D-5
!  initial-bound-perturbation                        -1.0
!  perturbation-small                                -1.0
!  reduce-perturbation-factor                        0.25
!  reduce-perturbation-multiplier                    0.01
!  insufficiently-feasible-tolerance                 0.01
!  implicit-variable-tolerance                       1.0D-5
!  mininum-initial-primal-feasibility                1.0
!  mininum-initial-dual-feasibility                  1.0
!  target-barrier-parameter                          -1.0
!  target-barrier-accept-fraction                    1.0
!  increase-barrier-parameter-by                     2.0
!  required-infeasibility-reduction                  0.99
!  pivot-tolerance-used                              1.0D-12
!  pivot-tolerance-used-for-dependencies             0.5
!  zero-pivot-tolerance                              1.0D-12
!  alpha-scaling-tolerance                           1.0D-2
!  identical-bounds-tolerance                        1.0D-15
!  maximum-cpu-time-limit                            -1.0
!  remove-linear-dependencies                        T
!  treat-zero-bounds-as-general                      F
!  just-find-feasible-point                          F
!  use-corrector-step                                F
!  balance-initial-complementarity                   F
!  record-x-status                                   T
!  record-c-status                                   T
! END WCP SPECIFICATIONS (DEFAULT)

!  Dummy arguments

      TYPE ( WCP_control_type ), INTENT( INOUT ) :: control        
      INTEGER, INTENT( IN ) :: device
      CHARACTER( LEN = 16 ), OPTIONAL :: alt_specname

!  Programming: Nick Gould and Ph. Toint, January 2002.

!  Local variables

      INTEGER, PARAMETER :: lspec = 47
      CHARACTER( LEN = 16 ), PARAMETER :: specname = 'WCP            '
      TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec

!  Define the keywords

!  Integer key-words

      spec(  1 )%keyword = 'error-printout-device'
      spec(  2 )%keyword = 'printout-device'
      spec(  3 )%keyword = 'print-level' 
      spec(  4 )%keyword = 'maximum-number-of-iterations'
      spec(  5 )%keyword = 'start-print'
      spec(  6 )%keyword = 'stop-print'
      spec( 32 )%keyword = 'initial-point-used'
      spec(  7 )%keyword = 'factorization-used'
      spec(  8 )%keyword = 'maximum-column-nonzeros-in-schur-complement'
      spec(  9 )%keyword = 'initial-integer-workspace'
      spec( 10 )%keyword = 'initial-real-workspace'
      spec( 11 )%keyword = 'maximum-refinements'
      spec( 12 )%keyword = 'maximum-poor-iterations-before-infeasible'
      spec( 38 )%keyword = 'perturbation-strategy'
      spec( 13 )%keyword = 'restore-problem-on-output'

!  Real key-words

      spec( 14 )%keyword = 'infinity-value'
      spec( 15 )%keyword = 'primal-accuracy-required'
      spec( 16 )%keyword = 'dual-accuracy-required'
      spec( 17 )%keyword = 'complementary-slackness-accuracy-required'
      spec( 18 )%keyword = 'mininum-initial-primal-feasibility'
      spec( 19 )%keyword = 'mininum-initial-dual-feasibility'
      spec( 20 )%keyword = 'target-barrier-parameter'
      spec( 21 )%keyword = 'poor-iteration-tolerance'
      spec( 22 )%keyword = 'initial-bound-perturbation'
      spec( 23 )%keyword = 'pivot-tolerance-used'
      spec( 24 )%keyword = 'pivot-tolerance-used-for-dependencies'
      spec( 25 )%keyword = 'zero-pivot-tolerance'
      spec( 30 )%keyword = 'alpha-scaling-tolerance'
      spec( 31 )%keyword = ''
      spec( 39 )%keyword = 'implicit-variable-tolerance'
      spec( 26 )%keyword = 'identical-bounds-tolerance'
      spec( 40 )%keyword = 'insufficiently-feasible-tolerance'
      spec( 41 )%keyword = 'reduce-perturbation-factor'
      spec( 42 )%keyword = 'reduce-perturbation-multiplier'
      spec( 43 )%keyword = 'perturbation-small'
      spec( 46 )%keyword = 'target-barrier-accept-fraction'
      spec( 47 )%keyword = 'increase-barrier-parameter-by'
      spec( 35 )%keyword = 'maximum-cpu-time-limit'

!  Logical key-words

      spec( 27 )%keyword = 'remove-linear-dependencies'
      spec( 28 )%keyword = 'treat-zero-bounds-as-general'
      spec( 29 )%keyword = 'just-find-feasible-point'
      spec( 33 )%keyword = 'balance-initial-complementarity'
      spec( 34 )%keyword = 'use-corrector-step'
      spec( 36 )%keyword = 'space-critical'
      spec( 37 )%keyword = 'deallocate-error-fatal'
      spec( 44 )%keyword = 'record-x-status'
      spec( 45 )%keyword = 'record-c-status'

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
      CALL SPECFILE_assign_value( spec( 32 ), control%initial_point,           &
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
      CALL SPECFILE_assign_value( spec( 38 ), control%perturbation_strategy,   &
                                   control%error )
      CALL SPECFILE_assign_value( spec( 13 ), control%restore_problem,         &
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
      CALL SPECFILE_assign_value( spec( 20 ), control%mu_target,               &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 21 ),                                  &
                                  control%required_infeas_reduction,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 22 ), control%perturb_start,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 23 ), control%pivot_tol,               &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 24 ),                                  &
                                  control%pivot_tol_for_dependencies,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 25 ), control%zero_pivot,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 30 ), control%alpha_scale,             &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 39 ), control%implicit_tol,            &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 40 ), control%insufficiently_feasible, &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 41 ), control%reduce_perturb_factor,   &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 42 ),                                  &
                                  control%reduce_perturb_multiplier,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 26 ), control%identical_bounds_tol,    &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 43 ), control%perturbation_small,      &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 46 ), control%mu_accept_fraction,      &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 47 ), control%mu_increase_factor,      &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 35 ), control%cpu_time_limit,          &
                                  control%error )

!  Set logical values

      CALL SPECFILE_assign_value( spec( 27 ), control%remove_dependencies,     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 28 ),                                  &
                                  control%treat_zero_bounds_as_general,        &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 29 ), control%just_feasible,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 33 ),                                  &
                                  control%balance_initial_complementarity,     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 34 ), control%use_corrector,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 36 ), control%space_critical,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 37 ),                                  &
                                  control%deallocate_error_fatal,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 44 ), control%record_x_status,         &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 45 ), control%record_c_status,         &
                                  control%error )

      RETURN

      END SUBROUTINE WCP_read_specfile

!-*-*-*-*-*-*-*-*-*-   W C P _ S O L V E  S U B R O U T I N E   -*-*-*-*-*-*-*

      SUBROUTINE WCP_solve( prob, data, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Finds a well-centred feasible point for the system
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
!    to be solved since the last call to WCP_initialize, and .FALSE. if
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
!        feasible region will be found.
!        %G (see below) need not be set
!
!     1  each component of the linear terms g will be one. 
!        %G (see below) need not be set
!
!     any other value - the gradients will be those given by %G (see below)
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
!   %Y_l, %Y_u are REAL arrays of length %m, which must be set by the user to
!    appropriate estimates of the values of the Lagrange multipliers 
!    corresponding to the general constraints c_l <= A x and A x  <= c_u
!    respectively. Y_l should be positive and Y_u should be negative.
!    On successful exit, they will contain the required vectors of Lagrange 
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
!   %Z_l, %Z_u are REAL arrays of length %m, which must be set by the user to
!    appropriate estimates of the values of the dual variables
!    corresponding to the simple bounds x_l <= x and x  <= x_u
!    respectively. Z_l should be positive and Z_u should be negative.
!    On successful exit, they will contain the required vectors of dual
!    variables.
!
!  data is a structure of type WCP_data_type which holds private internal data
!
!  control is a structure of type WCP_control_type that controls the 
!   execution of the subroutine and must be set by the user. Default values for
!   the elements may be set by a call to WCP_initialize. See WCP_initialize 
!   for details
!
!  inform is a structure of type WCP_inform_type that provides 
!    information on exit from WCP_solve. The component status 
!    has possible values:
!  
!     0 Normal termination with a locally optimal solution.
!
!   - 1 one of the restrictions 
!          prob%n     >=  1
!          prob%m     >=  0
!          prob%A%type in { 'DENSE', 'SPARSE_BY_ROWS', 'COORDINATE' }
!       has been violated.
!
!    -2 An allocation error occured; the status is given in the component
!       alloc_status.
!
!    -3 A deallocation error occured; the status is given in the component
!       alloc_status.
!
!    -4 Too many iterations have been performed. This may happen if
!       control%maxit is too small, but may also be symptomatic of 
!       a badly scaled problem.
!
!    -5 The constraints are inconsistent.
!
!    -6 The constraints appear to have no feasible point.
!
!    -7 The factorization failed; the return status from the factorization
!       package is given in the component factor_status.
!      
!    -8 The problem is so ill-conditoned that further progress is impossible.  
!
!  On exit from WCP_solve, other components of inform give the 
!  following:
!
!     alloc_status = The status of the last attempted allocation/deallocation 
!     iter   = The total number of iterations required.
!     factorization_integer = The total integer workspace required by the 
!              factorization.
!     factorization_real = The total real workspace required by the 
!              factorization.
!     nfacts = The total number of factorizations performed.
!     factorization_status = the return status from the matrix factorization
!              package.   
!     obj = the value of the objective function ||W*(x-x^0)||_2.
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
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      TYPE ( WCP_data_type ), INTENT( INOUT ) :: data
      TYPE ( WCP_control_type ), INTENT( IN ) :: control
      TYPE ( WCP_inform_type ), INTENT( OUT ) :: inform

!  Local variables

      INTEGER :: i, j, a_ne, n_depen, nbnds, lbreak, lbnds, nzc
      INTEGER :: dy_l_lower, dy_l_upper, dy_u_lower, dy_u_upper
      INTEGER :: dz_l_lower, dz_l_upper, dz_u_lower, dz_u_upper
      REAL :: time, time_start, dum
      REAL ( KIND = wp ) :: fixed_sum, av_bnd
      LOGICAL :: printi, remap_freed, reset_bnd, implicit
      CHARACTER ( LEN = 80 ) :: array_name
      TYPE ( FDC_control_type ) :: FDC_control        
      TYPE ( FDC_inform_type ) :: FDC_inform

      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' entering WCP_solve ' )" )

!  Initialize time

      CALL CPU_TIME( time_start )

!  Set initial timing breakdowns

      inform%time%total = 0.0 ; inform%time%analyse = 0.0
      inform%time%factorize = 0.0 ; inform%time%solve = 0.0
      inform%time%preprocess = 0.0 ; inform%time%find_dependent = 0.0

!  Initialize counts

      inform%status = 0 ; inform%alloc_status = 0 ; inform%bad_alloc = ''
      inform%iter = - 1 ; inform%nfacts = - 1
      inform%factorization_integer = - 1 ; inform%factorization_real = - 1 
      inform%obj = - one ; inform%non_negligible_pivot = zero
      inform%feasible = .FALSE. ; inform%factorization_status = 0

!  Basic single line of output per iteration

      printi = control%out > 0 .AND. control%print_level >= 1 

!  Ensure that input parameters are within allowed ranges

      IF ( prob%n < 1 .OR. prob%m < 0 .OR.                                     &
           .NOT. QPT_keyword_A( prob%A%type ) ) THEN
        inform%status = - 1
        IF ( control%error > 0 .AND. control%print_level > 0 )                 &
          WRITE( control%error, 2010 ) inform%status 
        CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
        GO TO 800 
      END IF 
      prob%Hessian_kind = 0

!  If required, write out problem 

      IF ( control%out > 0 .AND. control%print_level >= 20 ) THEN
        WRITE( control%out, "( ' n, m = ', 2I8 )" ) prob%n, prob%m
        WRITE( control%out, "( ' f = ', ES12.4 )" ) prob%f
        IF ( prob%gradient_kind == 0 ) THEN
          WRITE( control%out, "( ' G = zeros' )" )
        ELSE IF ( prob%gradient_kind == 1 ) THEN
          WRITE( control%out, "( ' G = ones' )" )
        ELSE
          WRITE( control%out, "( ' G = ', /, ( 5ES12.4 ) )" )                  &
            prob%G( : prob%n )
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
          inform%status = - 5
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
          inform%status = - 5
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

!  Allocate additional workspace

      array_name = 'wcp: data%RES_x'
      CALL SPACE_resize_array( prob%n, data%RES_x, inform%status,              &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%IW'
      CALL SPACE_resize_array( prob%n + 1, data%IW, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900
      
      IF ( control%record_x_status ) THEN
        array_name = 'wcp: inform%X_status'
        CALL SPACE_resize_array( prob%n, inform%X_status, inform%status,       &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900
      END IF

      IF ( control%record_c_status ) THEN
        array_name = 'wcp: inform%C_status'
        CALL SPACE_resize_array( prob%m, inform%C_status, inform%status,       &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900
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
        data%QPP_control%treat_zero_bounds_as_general =                        &
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

        IF ( data%QPP_inform%status /= 0 ) THEN
          inform%status = data%QPP_inform%status
          IF ( control%out > 0 .AND. control%print_level >= 5 )                &
            WRITE( control%out, "( ' status ', I3, ' after QPP_reorder ')" )   &
             data%QPP_inform%status
          IF ( control%error > 0 .AND. control%print_level > 0 )               &
            WRITE( control%error, 2010 ) inform%status 
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

          IF ( data%QPP_inform%status /= 0 ) THEN
            inform%status = data%QPP_inform%status
            IF ( control%out > 0 .AND. control%print_level >= 5 )              &
              WRITE( control%out, "( ' status ', I3, ' after QPP_apply ')" )   &
               data%QPP_inform%status
            IF ( control%error > 0 .AND. control%print_level > 0 )             &
              WRITE( control%error, 2010 ) inform%status 
            CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
            GO TO 800 
          END IF 
        END IF 
        data%trans = data%trans + 1
      END IF

!  Special case: no free variables

      IF ( prob%n == 0 ) THEN
        prob%Y_l( : prob%m ) = zero ; prob%Y_u( : prob%m ) = zero
        prob%Z_l( : prob%n ) = zero ; prob%Z_u( : prob%n ) = zero
        prob%C( : prob%m ) = zero
        CALL WCP_AX( prob%m, prob%C( : prob%m ), prob%m,                       &
                     prob%A%ptr( prob%m + 1 ) - 1, prob%A%val,                 &
                     prob%A%col, prob%A%ptr, prob%n, prob%X, '+ ')
        GO TO 700
      END IF

!  =================================================================
!  Check to see if the equality constraints are linearly independent
!  =================================================================

      IF ( .NOT. data%tried_to_remove_deps .AND. control%remove_dependencies ) &
        THEN

        CALL CPU_TIME( time ) 

        IF ( control%out > 0 .AND. control%print_level >= 1 )                  &
          WRITE( control%out,                                                  &
            "( /, 1X, I0, ' equalities from ', I0, ' constraints ' )" )        &
            data%dims%c_equality, prob%m

!  Set control parameters

        CALL FDC_initialize( FDC_control )
        FDC_control%error = 6
        FDC_control%out = 6
        FDC_control%CNTL = data%CNTL
        FDC_control%zero_pivot = control%zero_pivot
        FDC_control%pivot_tol = control%pivot_tol_for_dependencies
        FDC_control%max_infeas = control%stop_p

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
          inform%status = - 9
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
            IF ( inform%status /= 0 )                                          &
              WRITE( control%error, 2020 ) inform%status, 'WCP_dependent'
          END IF
          GO TO 700
        END IF

        IF ( control%out > 0 .AND. control%print_level >= 2 .AND. n_depen > 0 )&
          WRITE( control%out, "(/, ' The following ',I0,' constraints appear', &
       &         ' to be dependent', /, ( 8I8 ) )" ) n_depen, data%Index_C_freed

        remap_freed = n_depen > 0 .AND. prob%n > 0

!  Special case: no free variables

        IF ( prob%n == 0 ) THEN
          prob%Y_l( : prob%m ) = zero ; prob%Y_u( : prob%m ) = zero
          prob%Z_l( : prob%n ) = zero ; prob%Z_u( : prob%n ) = zero
          prob%C( : prob%m ) = zero
          CALL WCP_AX( prob%m, prob%C( : prob%m ), prob%m,                     &
                       prob%A%ptr( prob%m + 1 ) - 1, prob%A%val,               &
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

      array_name = 'wcp: data%C_freed'
      CALL SPACE_resize_array( n_depen, data%C_freed, inform%status,           &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900
        
!  Free the constraint bounds as required

        DO i = 1, n_depen
          j = data%Index_c_freed( i )
          data%C_freed( i ) = prob%C_l( j )
          prob%C_l( j ) = - control%infinity
          prob%C_u( j ) = control%infinity
          prob%Y_l( j ) = zero ; prob%Y_u( j ) = zero
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

        IF ( data%QPP_inform%status /= 0 ) THEN
          inform%status = data%QPP_inform%status
          IF ( control%out > 0 .AND. control%print_level >= 5 )                &
            WRITE( control%out, "( ' status ', I3, ' after QPP_reorder ')" )   &
             data%QPP_inform%status
          IF ( control%error > 0 .AND. control%print_level > 0 )               &
            WRITE( control%error, 2010 ) inform%status 
          CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
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

      array_name = 'wcp: data%Y'
      CALL SPACE_resize_array( prob%m, data%Y, inform%status,                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%SOL_y'
      CALL SPACE_resize_array( prob%m, data%SOL_y, inform%status,              &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%RES_y'
      CALL SPACE_resize_array( prob%m, data%RES_y, inform%status,              &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%BEST_y'
      CALL SPACE_resize_array( prob%m, data%BEST_y, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%SOL'
      CALL SPACE_resize_array( data%dims%v_e, data%SOL, inform%status,         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%RES'
      CALL SPACE_resize_array( data%dims%v_e, data%RES, inform%status,         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%BEST'
      CALL SPACE_resize_array( data%dims%v_e, data%BEST, inform%status,        &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%HX'
      CALL SPACE_resize_array( data%dims%v_e, data%HX, inform%status,          &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%GRAD_L'
      CALL SPACE_resize_array( data%dims%c_e, data%GRAD_L, inform%status,      &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DIST_X_l'
      CALL SPACE_resize_array( data%dims%x_free + 1, data%dims%x_l_end,        &
             data%DIST_X_l, inform%status,                                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DIST_Z_l'
      CALL SPACE_resize_array( data%dims%x_free + 1, data%dims%x_l_end,        &
             data%DIST_Z_l, inform%status,                                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%PERTURB_X_l'
      CALL SPACE_resize_array( data%dims%x_free + 1, data%dims%x_l_end,        &
             data%PERTURB_X_l, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%PERTURB_Z_l'
      CALL SPACE_resize_array( data%dims%x_free + 1, data%dims%x_l_end,        &
             data%PERTURB_Z_l, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DIST_X_u'
      CALL SPACE_resize_array( data%dims%x_u_start, prob%n,                    &
             data%DIST_X_u, inform%status,                                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DIST_Z_u'
      CALL SPACE_resize_array( data%dims%x_u_start, prob%n,                    &
             data%DIST_Z_u, inform%status,                                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%PERTURB_X_u'
      CALL SPACE_resize_array( data%dims%x_u_start, prob%n,                    &
             data%PERTURB_X_u, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900
      
      array_name = 'wcp: data%PERTURB_Z_u'
      CALL SPACE_resize_array( data%dims%x_u_start, prob%n,                    &
             data%PERTURB_Z_u, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900
      
      array_name = 'wcp: data%BARRIER_X'
      CALL SPACE_resize_array( data%dims%x_free + 1, prob%n,                   &
             data%BARRIER_X, inform%status,                                    &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%MU_X_l'
      CALL SPACE_resize_array( data%dims%x_free + 1, data%dims%x_l_end,        &
             data%MU_X_l, inform%status,                                       &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%MU_X_u'
      CALL SPACE_resize_array( data%dims%x_u_start, prob%n,                    &
             data%MU_X_u, inform%status,                                       &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DY_l'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_l_end,         &
             data%DY_l, inform%status,                                         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DIST_C_l'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_l_end,         &
             data%DIST_C_l, inform%status,                                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DIST_Y_l'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_l_end,         &
             data%DIST_Y_l, inform%status,                                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%PERTURB_C_l'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_l_end,         &
             data%PERTURB_C_l, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%PERTURB_Y_l'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_l_end,         &
             data%PERTURB_Y_l, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DY_u'
      CALL SPACE_resize_array( data%dims%c_u_start, data%dims%c_u_end,         &
             data%DY_u, inform%status,                                         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DIST_C_u'
      CALL SPACE_resize_array( data%dims%c_u_start, data%dims%c_u_end,         &
             data%DIST_C_u, inform%status,                                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DIST_Y_u'
      CALL SPACE_resize_array( data%dims%c_u_start, data%dims%c_u_end,         &
             data%DIST_Y_u, inform%status,                                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%PERTURB_C_u'
      CALL SPACE_resize_array( data%dims%c_u_start, data%dims%c_u_end,         &
             data%PERTURB_C_u, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%PERTURB_Y_u'
      CALL SPACE_resize_array( data%dims%c_u_start, data%dims%c_u_end,         &
             data%PERTURB_Y_u, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%C'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_u_end,         &
             data%C, inform%status,                                            &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%BARRIER_C'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_u_end,         &
             data%BARRIER_C, inform%status,                                    &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%MU_C_l'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_l_end,         &
             data%MU_C_l, inform%status,                                       &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%MU_C_u'
      CALL SPACE_resize_array( data%dims%c_u_start, data%dims%c_u_end,         &
             data%MU_C_u, inform%status,                                       &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%SCALE_C'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_u_end,         &
             data%SCALE_C, inform%status,                                      &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DELTA'
      CALL SPACE_resize_array( data%dims%v_e, data%DELTA, inform%status,       &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%RHS'
      CALL SPACE_resize_array( data%dims%v_e, data%RHS, inform%status,         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DZ_l'
      CALL SPACE_resize_array( data%dims%x_free + 1, data%dims%x_l_end,        &
             data%DZ_l, inform%status,                                         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%DZ_u'
      CALL SPACE_resize_array( data%dims%x_u_start, prob%n,                    &
             data%DZ_u, inform%status,                                         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

!  Allocate optional extra arrays

      nbnds = prob%n + data%dims%x_l_end - data%dims%x_free                    &
              - data%dims%x_u_start + data%dims%c_l_end - data%dims%c_l_start  & 
              + data%dims%c_u_end - data%dims%c_u_start + 3

      IF ( nbnds > 0 ) THEN
        lbnds = nbnds
        IF ( control%use_corrector ) THEN
          lbreak = 4 * nbnds
        ELSE
          lbreak = 2 * nbnds
        END IF
      ELSE
        lbnds = 0
        lbreak = 0
      END IF

      array_name = 'wcp: data%MU'
      CALL SPACE_resize_array( lbnds, data%MU, inform%status,                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%COEF0'
      CALL SPACE_resize_array( lbnds, data%COEF0, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%COEF1'
      CALL SPACE_resize_array( lbnds, data%COEF1, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%COEF2'
      CALL SPACE_resize_array( lbnds, data%COEF2, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%COEF3'
      CALL SPACE_resize_array( lbnds, data%COEF3, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%COEF4'
      CALL SPACE_resize_array( lbnds, data%COEF4, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%BREAKP'
      CALL SPACE_resize_array( lbreak, data%BREAKP, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      array_name = 'wcp: data%IBREAK'
      CALL SPACE_resize_array( lbreak, data%IBREAK, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 900

      IF ( control%use_corrector ) THEN

        array_name = 'wcp: data%DELTA_cor'
        CALL SPACE_resize_array( data%dims%v_e, data%DELTA_cor,                &
               inform%status, inform%alloc_status, array_name = array_name,    &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        array_name = 'wcp: data%DY_cor_l'
        CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_l_end,       &
               data%DY_cor_l, inform%status,                                   &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        dy_l_lower = data%dims%c_l_start
        dy_l_upper = data%dims%c_l_end

        array_name = 'wcp: data%DY_cor_u'
        CALL SPACE_resize_array( data%dims%c_u_start, data%dims%c_u_end,       &
               data%DY_cor_u, inform%status,                                   &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        dy_u_lower = data%dims%c_u_start
        dy_u_upper = data%dims%c_u_end

        array_name = 'wcp: data%DZ_cor_l'
        CALL SPACE_resize_array( data%dims%x_free + 1, data%dims%x_l_end,      &
               data%DZ_cor_l, inform%status,                                   &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        dz_l_lower = data%dims%x_free + 1
        dz_l_upper = data%dims%x_l_end

        array_name = 'wcp: data%DZ_cor_u'
        CALL SPACE_resize_array( data%dims%x_u_start, prob%n,                  &
               data%DZ_cor_u, inform%status,                                   &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        dz_u_lower = data%dims%x_u_start
        dz_u_upper = prob%n

      ELSE

        array_name = 'wcp: data%DELTA_cor'
        CALL SPACE_resize_array( 0, data%DELTA_cor, inform%status,             &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        array_name = 'wcp: data%DY_cor_l'
        CALL SPACE_resize_array( 0, data%DY_cor_l, inform%status,              &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        dy_l_lower = 1
        dy_l_upper = 0

        array_name = 'wcp: data%DY_cor_u'
        CALL SPACE_resize_array( 0, data%DY_cor_u, inform%status,              &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        dy_u_lower = 1
        dy_u_upper = 0

        array_name = 'wcp: data%DZ_cor_l'
        CALL SPACE_resize_array( 0, data%DZ_cor_l, inform%status,              &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        dz_l_lower = 1
        dz_l_upper = 0

        array_name = 'wcp: data%DZ_cor_u'
        CALL SPACE_resize_array( 0, data%DZ_cor_u, inform%status,              &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 900

        dz_u_lower = 1
        dz_u_upper = 0

      END IF

!  =================
!  Solve the problem
!  =================

      IF ( printi ) WRITE( control%out,                                        &
           "( /, ' <------ variable bounds ------>',                           &
        &        ' <----- constraint bounds ----->',                           &
        &     /, '    free   below    both   above',                           &
        &        '   equal   below    both   above',                           &
        &     /,  8I8 )" )                                                     &
            data%dims%x_free, data%dims%x_u_start - data%dims%x_free - 1,      &
            data%dims%x_l_end - data%dims%x_u_start + 1,                       &
            prob%n - data%dims%x_l_end, data%dims%c_equality,                  &
            data%dims%c_u_start - data%dims%c_equality - 1,                    &
            data%dims%c_l_end - data%dims%c_u_start + 1,                       &
            prob%m - data%dims%c_l_end

      IF ( prob%gradient_kind == 0 .OR. prob%gradient_kind == 1 ) THEN
        CALL WCP_solve_main( data%dims, prob%n, prob%m,                        &
                             prob%A%val, prob%A%col, prob%A%ptr,               &
                             prob%C_l, prob%C_u, prob%X_l, prob%X_u,           &
                             prob%C, prob%X, prob%Y_l, prob%Y_u,               &
                             prob%Z_l, prob%Z_u, data%RES_x,                   &
                             data%Y, data%SOL_y, data%RES_y, data%BEST_y,      &
                             data%SOL, data%RES, data%BEST, data%HX,           &
                             data%GRAD_L, data%DIST_X_l, data%DIST_Z_l,        &
                             data%DIST_X_u, data%DIST_Z_u,                     &
                             data%BARRIER_X, data%MU_X_l, data%MU_X_u,         &
                             data%DY_l, data%DIST_C_l, data%DIST_Y_l,          &
                             data%DY_u, data%DIST_C_u, data%DIST_Y_u,          &
                             data%C, data%BARRIER_C,                           &
                             data%MU_C_l, data%MU_C_u,                         &
                             data%SCALE_C, data%DELTA, data%RHS,               &
                             data%DZ_l, data%DZ_u,                             &
                             data%PERTURB_X_l, data%PERTURB_X_u,               &
                             data%PERTURB_Y_l, data%PERTURB_Y_u,               &
                             data%PERTURB_Z_l, data%PERTURB_Z_u,               &
                             data%PERTURB_C_l, data%PERTURB_C_u,               &
                             data%Abycol_val, data%DIAG_X, data%DIAG_C,        &
                             data%IW, data%K_colptr, data%Abycol_ptr,          &
                             data%Abycol_row, data%K, data%FACTORS, prob%f,    &
                             data%MU, data%COEF0, data%COEF1, data%COEF2,      &
                             data%COEF3, data%COEF4, data%DELTA_cor,           &
                             data%DY_cor_l, dy_l_lower, dy_l_upper,            &
                             data%DY_cor_u, dy_u_lower, dy_u_upper,            &
                             data%DZ_cor_l, dz_l_lower, dz_l_upper,            &
                             data%DZ_cor_u, dz_u_lower, dz_u_upper,            &
                             data%BREAKP, data%IBREAK, data%CNTL,              &
                             control, inform, prob%gradient_kind )
      ELSE
        CALL WCP_solve_main( data%dims, prob%n, prob%m,                        &
                             prob%A%val, prob%A%col, prob%A%ptr,               &
                             prob%C_l, prob%C_u, prob%X_l, prob%X_u,           &
                             prob%C, prob%X, prob%Y_l, prob%Y_u,               &
                             prob%Z_l, prob%Z_u, data%RES_x,                   &
                             data%Y, data%SOL_y, data%RES_y, data%BEST_y,      &
                             data%SOL, data%RES, data%BEST, data%HX,           &
                             data%GRAD_L, data%DIST_X_l, data%DIST_Z_l,        &
                             data%DIST_X_u, data%DIST_Z_u,                     &
                             data%BARRIER_X, data%MU_X_l, data%MU_X_u,         &
                             data%DY_l, data%DIST_C_l, data%DIST_Y_l,          &
                             data%DY_u, data%DIST_C_u, data%DIST_Y_u,          &
                             data%C, data%BARRIER_C,                           &
                             data%MU_C_l, data%MU_C_u,                         &
                             data%SCALE_C, data%DELTA, data%RHS,               &
                             data%DZ_l, data%DZ_u,                             &
                             data%PERTURB_X_l, data%PERTURB_X_u,               &
                             data%PERTURB_Y_l, data%PERTURB_Y_u,               &
                             data%PERTURB_Z_l, data%PERTURB_Z_u,               &
                             data%PERTURB_C_l, data%PERTURB_C_u,               &
                             data%Abycol_val, data%DIAG_X, data%DIAG_C,        &
                             data%IW, data%K_colptr, data%Abycol_ptr,          &
                             data%Abycol_row, data%K, data%FACTORS, prob%f,    &
                             data%MU, data%COEF0, data%COEF1, data%COEF2,      &
                             data%COEF3, data%COEF4, data%DELTA_cor,           &
                             data%DY_cor_l, dy_l_lower, dy_l_upper,            &
                             data%DY_cor_u, dy_u_lower, dy_u_upper,            &
                             data%DZ_cor_l, dz_l_lower, dz_l_upper,            &
                             data%DZ_cor_u, dz_u_lower, dz_u_upper,            &
                             data%BREAKP, data%IBREAK, data%CNTL,              &
                             control, inform, prob%gradient_kind, G = prob%G )
      END IF

      inform%obj = inform%obj + half * fixed_sum

!  If some of the constraints were freed during the computation, refix them now

      IF ( remap_freed ) THEN
        CALL CPU_TIME( time )
        IF ( control%record_x_status )                                         &
          CALL SORT_inverse_permute( data%QPP_map_freed%n,                     &
            data%QPP_map_freed%x_map,                                          &
            IX = inform%X_status( : data%QPP_map_freed%n ) )
        IF ( control%record_c_status ) THEN 
          inform%C_status( prob%m + 1 : data%QPP_map_freed%m ) = 4
          CALL SORT_inverse_permute( data%QPP_map_freed%m,                     &
            data%QPP_map_freed%c_map,                                          &
            IX = inform%C_status( : data%QPP_map_freed%m ) )
        END IF
        CALL QPP_restore( data%QPP_map_freed, data%QPP_inform, prob,           &
                          get_all = .TRUE. )
        CALL QPP_terminate( data%QPP_map_freed, data%QPP_control,              &
                            data%QPP_inform )
        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum
        data%dims = data%dims_save_freed
        inform%c_implicit = inform%c_implicit + n_depen

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
        IF ( control%record_x_status )                                         &
          CALL SORT_inverse_permute( data%QPP_map%n,                           &
            data%QPP_map%x_map, IX = inform%X_status( : data%QPP_map%n ) )
        IF ( control%record_c_status )                                         & 
          CALL SORT_inverse_permute( data%QPP_map%m,                           &
            data%QPP_map%c_map, IX = inform%C_status( : data%QPP_map%m ) )

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

      IF ( control%print_level > 1 .AND. control%record_x_status ) THEN
        implicit = .FALSE.
        DO i = 1, prob%n 
          IF ( inform%X_status( i ) == - 1 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2030 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'LOWER', prob%X( i ), prob%X_l( i ), &
              prob%X_u( i ), prob%Z_l( i ) + prob%Z_u( i )
          ELSE IF ( inform%X_status( i ) == 1 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2030 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'UPPER', prob%X( i ), prob%X_l( i ), &
              prob%X_u( i ), prob%Z_l( i ) + prob%Z_u( i )
          ELSE IF ( inform%X_status( i ) == - 2 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2030 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'LFREE', prob%X( i ), prob%X_l( i ), &
              prob%X_u( i ), prob%Z_l( i ) + prob%Z_u( i )
          ELSE IF ( inform%X_status( i ) == 2 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2030 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'UFREE', prob%X( i ), prob%X_l( i ), &
              prob%X_u( i ), prob%Z_l( i ) + prob%Z_u( i )
          ELSE IF ( inform%X_status( i ) == - 3 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2030 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'FREE ', prob%X( i ), prob%X_l( i ), &
              prob%X_u( i ), prob%Z_l( i ) + prob%Z_u( i )
          END IF
        END DO
      END IF

      IF ( control%print_level > 1 .AND. control%record_c_status ) THEN
        implicit = .FALSE.
        DO i = 1, prob%m
          IF ( inform%C_status( i ) == - 1 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2040 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'LOWER', prob%C( i ), prob%C_l( i ), &
              prob%C_u( i ), prob%Y_l( i ) + prob%Y_u( i )
          ELSE IF ( inform%C_status( i ) == 1 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2040 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'UPPER', prob%C( i ), prob%C_l( i ), &
              prob%C_u( i ), prob%Y_l( i ) + prob%Y_u( i )
          ELSE IF ( inform%C_status( i ) == - 2 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2040 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'LFREE', prob%C( i ), prob%C_l( i ), &
              prob%C_u( i ), prob%Y_l( i ) + prob%Y_u( i )
          ELSE IF ( inform%C_status( i ) == 2 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2040 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'UFREE', prob%C( i ), prob%C_l( i ), &
              prob%C_u( i ), prob%Y_l( i ) + prob%Y_u( i )
         ELSE IF ( inform%C_status( i ) == - 3 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2040 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'FREE ', prob%C( i ), prob%C_l( i ), &
              prob%C_u( i ), prob%Y_l( i ) + prob%Y_u( i )
          ELSE IF ( inform%C_status( i ) == 4 ) THEN
            IF ( .NOT. implicit ) THEN
              WRITE( control%out, 2040 ) ; implicit = .TRUE. ; END IF
            WRITE( control%out, 2050 ) i, 'DEPEN', prob%C( i ), prob%C_l( i ), &
              prob%C_u( i ), prob%Y_l( i ) + prob%Y_u( i )
          END IF
        END DO
      END IF


!  Compute total time

      CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
      IF ( printi ) WRITE( control%out, 2000 ) inform%time%total,              &
        inform%time%preprocess, inform%time%analyse, inform%time%factorize,    &
        inform%time%solve

  800 CONTINUE 
      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' leaving WCP_solve ' )" )

      RETURN  

!  Allocation error

  900 CONTINUE
!     inform%status = - 2
      CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
!     IF ( printi ) WRITE( control%out, 2900 ) bad_alloc, inform%alloc_status
      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' leaving WCP_solve ' )" )

      RETURN  

!  Non-executable statements

 2000 FORMAT( /, 14X, ' =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=',&
              /, 14X, ' =                  WCP total time                   =',&
              /, 14X, ' =', 16X, 0P, F12.2, 23x, '='                           &
              /, 14X, ' =    preprocess    analyse    factorize     solve   =',&
              /, 14X, ' =', 4F12.2, 3x, '=',                                   &
              /, 14X, ' =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=' )
 2010 FORMAT( ' ',  /, '   **  Error return ', I0, ' from WCP ' ) 
 2020 FORMAT( '   **  Error return ', I0, ' from ', A ) 
 2030 FORMAT( /, ' Implicit bounds: ', /, '                    ',              &
                 '        <------ Bounds ------> ', /                          &
                 '      #  state    value   ',                                 &
                 '    Lower       Upper       Dual ' ) 
 2040 FORMAT( /, ' Implicit constraints: ', /, '                   ',          &
                 '        <------ Bounds ------> ', /                          &
                 '      #  state    value   ',                                 &
                 '    Lower       Upper     Multiplier ' ) 
 2050 FORMAT( I7, 1X, A6, 4ES12.4 ) 

!  End of WCP_solve

      END SUBROUTINE WCP_solve

!-*-*-*-*-*-   W C P _ S O L V E _ M A I N   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE WCP_solve_main( dims, n, m, A_val, A_col, A_ptr,              &
                                 C_l, C_u, X_l, X_u, C_RES, X, Y_l, Y_u,       &
                                 Z_l, Z_u, RES_x, Y, SOL_y, RES_y, BEST_y,     &
                                 SOL, RES, BEST, HX, GRAD_L, DIST_X_l,         &
                                 DIST_Z_l, DIST_X_u, DIST_Z_u,                 &
                                 BARRIER_X, MU_X_l, MU_X_u, DY_l,              &
                                 DIST_C_l, DIST_Y_l, DY_u, DIST_C_u,           &
                                 DIST_Y_u, C, BARRIER_C, MU_C_l, MU_C_u,       &
                                 SCALE_C, DELTA, RHS, DZ_l, DZ_u,              &
                                 PERTURB_X_l, PERTURB_X_u,                     &
                                 PERTURB_Y_l, PERTURB_Y_u,                     &
                                 PERTURB_Z_l, PERTURB_Z_u,                     &
                                 PERTURB_C_l, PERTURB_C_u,                     &
                                 Abycol_val, DIAG_X, DIAG_C, IW,               &
                                 K_colptr, Abycol_ptr, Abycol_row, K,          &
                                 FACTORS, f, MU, COEF0, COEF1, COEF2, COEF3,   &
                                 COEF4, DELTA_cor,                             &
                                 DY_cor_l, dy_l_lower, dy_l_upper,             &
                                 DY_cor_u, dy_u_lower, dy_u_upper,             &
                                 DZ_cor_l, dz_l_lower, dz_l_upper,             &
                                 DZ_cor_u, dz_u_lower, dz_u_upper,             &
                                 BREAKP, IBREAK, CNTL, control, inform,        &
                                 gradient_kind, G )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Finds a well-centered point within the polytope
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
!  dims is a structure of type WCP_data_type, whose components hold SCALAR
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
!  n, m, ..., Z_l, Z_u exactly as for prob% in WCP_solve
!
!  control and inform are exactly as for WCP_solve
!
!  The remaining arguments are used as internal workspace, and need not be 
!  set on entry
!  
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, gradient_kind
      INTEGER, INTENT( IN ) :: dy_l_lower, dy_l_upper, dy_u_lower, dy_u_upper
      INTEGER, INTENT( IN ) :: dz_l_lower, dz_l_upper, dz_u_lower, dz_u_upper
      REAL ( KIND = wp ), INTENT( IN ) :: f
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      INTEGER, INTENT( IN ), DIMENSION( A_ptr( m + 1 ) - 1 ) :: A_col
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: C_l, C_u
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X_l, X_u
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: Y_l, Y_u
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: Z_l, Z_u
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ), OPTIONAL :: G
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( A_ptr( m + 1 ) - 1 ) :: A_val
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: RES_x
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
                          DIMENSION( m ) :: C_RES, Y, SOL_y, BEST_y, RES_y
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
             DIMENSION( dims%v_e ) :: SOL, RES, BEST
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%v_e ) :: DELTA, RHS
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( dims%v_e ) :: HX
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( dims%c_e ) :: GRAD_L
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%x_free + 1 : dims%x_l_end ) :: DZ_l,         &
                           DIST_X_l, DIST_Z_l, PERTURB_X_l, PERTURB_Z_l, MU_X_l
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( dims%x_u_start : n ) ::    &
                     DZ_u, DIST_X_u, DIST_Z_u, PERTURB_X_u, PERTURB_Z_u, MU_X_u
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
                          DIMENSION( dims%x_free + 1 : n ) :: BARRIER_X
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%c_l_start : dims%c_l_end ) :: DY_l,               &
                           DIST_C_l, DIST_Y_l, PERTURB_C_l, PERTURB_Y_l, MU_C_l
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%c_u_start : dims%c_u_end ) :: DY_u,               &
                           DIST_C_u, DIST_Y_u, PERTURB_C_u, PERTURB_Y_u, MU_C_u
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
             DIMENSION( dims%c_l_start : dims%c_u_end ) :: C, BARRIER_C, SCALE_C

!  pointers (which will be allocatable arrays in f2000!!)

      INTEGER, ALLOCATABLE, DIMENSION( : ) :: IW, K_colptr
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: Abycol_ptr, Abycol_row
      INTEGER, DIMENSION( : ) :: IBREAK
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIAG_X, DIAG_C
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Abycol_val
      REAL ( KIND = wp ), DIMENSION( : ) :: COEF0, COEF1, COEF2, COEF3, COEF4
      REAL ( KIND = wp ), DIMENSION( : ) :: MU, BREAKP
      REAL ( KIND = wp ), DIMENSION( dy_l_lower : dy_l_upper ) :: DY_cor_l
      REAL ( KIND = wp ), DIMENSION( dy_u_lower : dy_u_upper ) :: DY_cor_u
      REAL ( KIND = wp ), DIMENSION( dz_l_lower : dz_l_upper ) :: DZ_cor_l
      REAL ( KIND = wp ), DIMENSION( dz_u_lower : dz_u_upper ) :: DZ_cor_u
      REAL ( KIND = wp ), DIMENSION( : ) :: DELTA_cor

      TYPE ( SMT_type ), INTENT( INOUT ) :: K
      TYPE ( SILS_factors ), INTENT( INOUT ) :: FACTORS
      TYPE ( SILS_control ), INTENT( INOUT ) :: CNTL
      TYPE ( WCP_control_type ), INTENT( IN ) :: control        
      TYPE ( WCP_inform_type ), INTENT( INOUT ) :: inform

!  Parameters

      REAL ( KIND = wp ), PARAMETER :: eta = tenm4
      REAL ( KIND = wp ), PARAMETER :: sigma_max = point01
      REAL ( KIND = wp ), PARAMETER :: degen_tol = tenm5

!  Local variables

      INTEGER :: A_ne, i, l, lk, liw, ierr, start_print, stop_print, print_level
      INTEGER :: ldiag_x, nbnds, nbnds_x, nbnds_c
      INTEGER :: itref_max, ldiag_c_l, ldiag_c_u, cs_bad
      INTEGER :: out, error, factor, nnzks, it_best, infeas_max
      REAL :: time, dum, time_start
      REAL ( KIND = wp ) :: pjgnrm, errorg, mu_target, pmax, amax, alpha
      REAL ( KIND = wp ) :: gi, merit, res_prim_dual, max_zr, max_yr, slkmax
      REAL ( KIND = wp ) :: res_prim, res_dual, slknes, slkmin, dist_bound
      REAL ( KIND = wp ) :: slknes_x, slknes_c, slkmax_x, slkmax_c, slknes_req
      REAL ( KIND = wp ) :: slkmin_x, slkmin_c, merit_best
      REAL ( KIND = wp ) :: required_infeas_reduction
      REAL ( KIND = wp ) :: prfeas, dufeas, p_min, p_max, d_min, d_max, balance
      REAL ( KIND = wp ) :: errorc, one_minus_alpha, pivot_tol, alpha_max, gmax
      REAL ( KIND = wp ) :: pmax_cor, omega_l, omega_u, perturb, max_r, mu_scale
      REAL ( KIND = wp ) :: reduce_perturb_factor, one_minus_red_pert_fac
      REAL ( KIND = wp ) :: bound_average, bound_length, perturb_max, alpha_est
      REAL ( KIND = wp ) :: perturb_x, perturb_c, perturb_y, perturb_z, too_close
      REAL ( KIND = wp ) :: mu_l, mu_u, max_xr, max_cr, perturbation_small
!     REAL ( KIND = wp ) :: old_merit, min_mu, nu
      LOGICAL :: set_printt, set_printi, set_printw, set_printd, set_printe
      LOGICAL :: set_printp, printt, printi, printp, printe, printd, printw
      LOGICAL :: get_factors, refact, use_corrector, maxpiv, reset_mu
      LOGICAL :: now_interior, now_feasible, start_major
      LOGICAL :: use_scale_c = .FALSE.
!     LOGICAL :: use_scale_c = .TRUE.
      CHARACTER ( LEN = 1 ) :: re, co, al
      CHARACTER ( LEN = 2 ) :: coal
      TYPE ( SILS_finfo ) :: FINFO
      INTEGER :: sif = 50
!     LOGICAL :: generate_sif = .TRUE.
      LOGICAL :: generate_sif = .FALSE.

!     reset_mu = .TRUE.
      reset_mu = .FALSE.

      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' entering WCP_solve_main ' )" )

! -------------------------------------------------------------------
!  If desired, generate a SIF file for problem passed 

      IF ( generate_sif .AND. PRESENT( G ) ) THEN
        WRITE( sif, "( 'NAME          WCPB_OUT', //, 'VARIABLES', / )" )
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
            WRITE( sif, "( '    RANGE    ', ' C', I8, ' ', ES12.5 )" )         &
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

!  ==================
!  Input error checks
!  ==================

!  If there are no variables, exit

      IF ( n == 0 ) THEN 
        i = COUNT( ABS( C_l( : dims%c_equality ) ) > control%stop_p ) +        &
            COUNT( C_l( dims%c_l_start : dims%c_l_end ) > control%stop_p ) +   &
            COUNT( C_u( dims%c_u_start : dims%c_u_end ) < - control%stop_p )
        IF ( i == 0 ) THEN
          inform%status = 0
        ELSE
          inform%status = - 6
        END IF
        C_RES = zero ; Y_l = zero ; Y_u = zero
        inform%obj = zero
        GO TO 810
      END IF 

!  Check that range constraints are not simply fixed variables,
!  and that the upper bounds are larger than the corresponing lower bounds

      DO i = dims%x_u_start, dims%x_l_end
        IF ( X_u( i ) - X_l( i ) <= epsmch ) THEN 
          inform%status = - 5 ; GO TO 700 ; END IF
      END DO

      DO i = dims%c_u_start, dims%c_l_end
        IF ( C_u( i ) - C_l( i ) <= epsmch ) THEN 
          inform%status = - 5 ; GO TO 700 ; END IF
      END DO

!  Set control parameters

      itref_max = control%itref_max
      prfeas = MAX( control%prfeas, epsmch )
      dufeas = MAX( control%dufeas, epsmch )
      required_infeas_reduction = MAX( epsmch,                                 &
        MIN( control%required_infeas_reduction ** 2, one - epsmch ) )
      infeas_max = MAX( 0, control%infeas_max )

      CNTL%u = control%pivot_tol
      CNTL%lp = - 1 ; CNTL%mp = - 1 ; CNTL%wp = - 1
      IF ( print_level > 3 ) THEN
        IF ( print_level < 10 ) THEN
          CNTL%ldiag = 1
        ELSE
          CNTL%ldiag = 2
        END IF
      END IF

!  Record array size

      A_ne = A_ptr( m + 1 ) - 1

!  If required, write out the data matrix for the problem

      IF ( printd ) WRITE( out, 2150 ) ' a ', ( ( i, A_col( l ), A_val( l ),   &
                           l = A_ptr( i ), A_ptr( i + 1 ) - 1 ), i = 1, m )

      IF ( control%balance_initial_complementarity ) THEN
        IF ( control%mu_target <= zero ) THEN
          balance = ten
        ELSE
          balance = control%mu_target
        END IF
      END IF

!  Set the default simple bound perturbations, perturb

      IF ( control%perturbation_strategy <= 0 ) THEN
        perturb = zero
      ELSE
        IF ( control%perturb_start < zero ) THEN
          perturb = MIN( control%stop_p, control%stop_d, control%stop_c )
        ELSE
          perturb = control%perturb_start
        END IF
      END IF
      reduce_perturb_factor = control%reduce_perturb_factor
!     too_close = point01
      too_close = zero

!  =============================
!  Find a suitable initial point
!  =============================

      IF ( control%initial_point >= 1 ) THEN

!  Solve the problem of minimizing 
!     1/2||x-x_m||^2_(X_D^-2) + 1/2||c-c_m||^2_(C_D^-2) + g^T x
!  subject to A x - c = 0, where x_m and c_m are "bound average" 
!  values for x and c and X_d/C_d are "bound lengths"; 

!  for two-sided bounds x_l <= x <= x_u,
!    x_m = 1/2 (x_l + x_u ) and x_d = x_u - x_l
!  for one-sided bounds x_l <= x, 
!    x_m = x_l + pr_feas and x_d = pr_feas
!  for free variables
!    x_m = 0 and x_d = infinite

!  Equivalently, solve the linear system

!   ( X_D^-2   0     A^T ) ( x )   ( X_D^-2 x_m )
!   (   0    C_D^-2  -I  ) ( c ) = ( C_D^-2 c_m )
!   (   A     -I      0  ) ( y )   (     0      )

! ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Analyse the sparsity pattern of the required matrix
! ::::::::::::::::::::::::::::::::::::::::::::::::::::

        SCALE_C = one
        CALL CPU_TIME( time ) 
        CALL WCP_analyse( dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C,      &
                          factor, nnzks, lk, liw, ldiag_x, ldiag_c_l,          &
                          ldiag_c_u, IW, Abycol_val, Abycol_row, Abycol_ptr,   &
                          K_colptr, DIAG_X, DIAG_C, K, FACTORS, CNTL,          &
                          print_level, control, inform )
        CALL CPU_TIME( dum ) ; dum = dum - time
        IF ( printt ) WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
        inform%time%analyse = inform%time%analyse + dum

!  Set up the bound averages and bound lengths; also set up the
!  weights (BARRIER) and right-hand sides (RHS) for the problem
       
!       dist_bound = prfeas
        dist_bound = ten

  200   CONTINUE
        RHS( : dims%x_free ) = zero
        DO i = dims%x_free + 1, dims%x_l_start - 1
          bound_average = dist_bound
          bound_length = dist_bound
          BARRIER_X( i ) =  one / bound_length ** 2
          RHS( i ) = bound_average / bound_length ** 2
        END DO
        DO i = dims%x_l_start, dims%x_u_start - 1
          bound_average = X_l( i ) + dist_bound
          bound_length = dist_bound
          BARRIER_X( i ) = one / bound_length ** 2
          RHS( i ) = bound_average / bound_length ** 2
        END DO
        DO i = dims%x_u_start, dims%x_l_end
          bound_average = half * ( X_l( i ) + X_u( i ) )
          bound_length = X_u( i ) - X_l( i )
          BARRIER_X( i ) =  one / bound_length ** 2
          RHS( i ) = bound_average / bound_length ** 2
        END DO
        DO i = dims%x_l_end + 1, dims%x_u_end
          bound_average = X_u( i ) - dist_bound
          bound_length = dist_bound
          BARRIER_X( i ) =  one / bound_length ** 2
          RHS( i ) = bound_average / bound_length ** 2
        END DO
        DO i = dims%x_u_end + 1, n
          bound_average = - dist_bound
          bound_length = dist_bound
          BARRIER_X( i ) =  one / bound_length ** 2
          RHS( i ) = bound_average / bound_length ** 2
        END DO

        IF ( gradient_kind == 1 ) THEN
          RHS( : n ) = RHS( : n ) - one
        ELSE IF ( gradient_kind /= 0 ) THEN
          RHS( : n ) = RHS( : n ) - G
        END IF

        DO i = dims%c_l_start, dims%c_u_start - 1
          bound_average = C_l( i ) + dist_bound
          bound_length = dist_bound
          BARRIER_C( i ) =  one / bound_length ** 2
          RHS( dims%c_b + i ) = bound_average / bound_length ** 2
        END DO
        DO i = dims%c_u_start, dims%c_l_end
          bound_average = half * ( C_l( i ) + C_u( i ) )
          bound_length = C_u( i ) - C_l( i )
          BARRIER_C( i ) = one / bound_length ** 2
          RHS( dims%c_b + i ) = bound_average / bound_length ** 2
        END DO
        DO i = dims%c_l_end + 1, dims%c_u_end
          bound_average = C_u( i ) - dist_bound
          bound_length = dist_bound
          BARRIER_C( i ) = one / bound_length ** 2
          RHS( dims%c_b + i ) = bound_average / bound_length ** 2
        END DO
        RHS( dims%y_s : dims%y_i - 1 ) = C_l( : dims%c_equality )
        RHS( dims%y_i : dims%y_e ) = zero

!  Factorize the required matrix
        
        get_factors = .TRUE. ; re = 'r' 
        CNTL%liw = MAX( 2 * inform%factorization_integer, control%indmin )
        CNTL%la = MAX( 2 * inform%factorization_real, control%valmin )

        CALL CPU_TIME( time ) 

!  Form the Schur complement matrix

        IF ( factor == 0 .OR. factor == 1 ) THEN
          DIAG_X( : dims%x_free ) = zero
          DIAG_X( dims%x_free + 1 : ) = BARRIER_X
          DIAG_C = BARRIER_C

!         WRITE( 6, "( ' barrier_x:', /, ( 5ES12.4 ) )" ) BARRIER_X
!         WRITE( 6, "( ' barrier_c:', /, ( 5ES12.4 ) )" ) BARRIER_C

          IF ( printw ) WRITE( out,                                            &
             "( ' ...... formation of Schur complement  .......... ' )" )

          IF ( m > 0 ) THEN
            CALL WCP_form_Schur_complement(                                    &
                     dims, n, m, A_ne, Abycol_val, Abycol_row, Abycol_ptr,     &
                     DIAG_X, SCALE_C, DIAG_C, lk, K%val, K%row,        &
                     K_colptr, K%ne, ierr, IW( : n ), IW( n + 1 : ), .TRUE.,         &
                     control%error, print_level )
!           inform%nfacts = inform%nfacts + 1 
          ELSE
            FINFO%flag = 0 ; get_factors = .FALSE.
          END IF
          IF ( printw ) WRITE( out,                                            &
            "( ' ............ Schur complement formed ............ ' )" )
        ELSE
          K%val( nnzks + 1 : nnzks + dims%x_free ) = zero
          K%val( nnzks + dims%x_free + 1 : nnzks + n ) = BARRIER_X
          K%val( nnzks + dims%c_s : nnzks + dims%c_e ) = BARRIER_C
        END IF

! ::::::::::::::::::::::::::::::
!  Factorize the required matrix
! ::::::::::::::::::::::::::::::

        IF ( get_factors ) THEN

!         WRITE( 6, "( ' K: n, nnz ', 2I4 )" ) K%n, K%ne
!         WRITE( 6, "( 3 ( 2I7, ES12.4 ) )" )                                  &
!           ( K%row( i ), K%col( i ), K%val( i ), i = 1, K%ne )
!         WRITE( 22, "( 2I6 )" ) K%n, K%ne
!         WRITE( 22, "( ( 10I6 ) )" ) K%row( : K%ne )
!         WRITE( 22, "( ( 10I6 ) )" ) K%col( : K%ne ) 
!         WRITE( 22, "( ( 3ES24.16 ) )" ) K%val( : K%ne ) 

          IF ( printw ) WRITE( out,                                            &
             "( ' ......... factorization of KKT matrix ...............' )" )
          CALL SILS_factorize( K, FACTORS, CNTL, FINFO )
          IF ( printw ) WRITE( out,                                            &
             "( ' ............... end of factorization ...............' )" )

!  Record the storage required

          inform%nfacts = inform%nfacts + 1 
          inform%factorization_integer = FINFO%nirbdu 
          inform%factorization_real = FINFO%nrlbdu 

!  Test that the factorization succeeded

          inform%factorization_status = FINFO%flag
          IF ( FINFO%flag < 0 ) THEN
            IF ( printe ) WRITE( error, 2040 ) FINFO%flag, 'SILS_factorize'

!  It didn't. We might have run out of options

            IF ( factor == 2 .AND. maxpiv ) THEN
              inform%status = - 7 ; GO TO 700
              
!  ... or we may change the method

            ELSE IF ( factor < 2 .AND. maxpiv ) THEN
              factor = 2
              pivot_tol = control%pivot_tol
              maxpiv = pivot_tol >= half
              IF ( printi )                                                    &
                WRITE( out, "( '       Switching to augmented system method' )" )

!  Re-analyse the sparsity pattern of the required matrix

              CALL CPU_TIME( time ) 
              CALL WCP_analyse(                                                &
                dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C, factor,        &
                nnzks,lk, liw, ldiag_x, ldiag_c_l, ldiag_c_u, IW,              &
                Abycol_val, Abycol_row, Abycol_ptr, K_colptr, DIAG_X,          &
                DIAG_C, K, FACTORS, CNTL, print_level, control, inform )
              CALL CPU_TIME( dum ) ; dum = dum - time
              IF ( printt )                                                    &
                WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
              inform%time%analyse = inform%time%analyse + dum

!  ... or we can increase the pivot tolerance

            ELSE
              maxpiv = .TRUE.
              pivot_tol = half
              IF ( printi )                                                    &
                WRITE( out, "( '       Pivot tolerance increased ' )" )
            END IF
            alpha = zero
            GO TO 200

!  Record warning conditions

          ELSE IF ( FINFO%flag > 0 ) THEN
            IF ( printt ) THEN
              WRITE( out, 2050 ) FINFO%flag, 'SILS_factorize'
              IF ( FINFO%flag == 4 ) WRITE( out, "( ' ** Matrix has ', I7, &
             &   ' zero eigenvalues ' )" ) K%n - finfo%rank
            END IF 
          END IF 
  
          IF ( printt ) WRITE( out,                                            &
            "( ' real/integer space used for factors ', 2I10 )" )              &
              FINFO%nrlbdu, FINFO%nirbdu
  
        ELSE
          inform%factorization_integer = 0 
          inform%factorization_real = 0
        END IF
        CALL CPU_TIME( dum ) ; dum = dum - time
  
        IF ( printt ) THEN
          WRITE( out, "( ' ** factorize time = ', F10.2 ) " ) dum
          WRITE( out, 2060 ) inform%factorization_integer,                     &
                             inform%factorization_real
        END IF 
        inform%time%factorize = inform%time%factorize + dum
 
!  :::::::::::::::::::::::::::::::::::::::::::::::::
!  Solve the linear system to find the initial point
!  :::::::::::::::::::::::::::::::::::::::::::::::::

        IF ( printw ) WRITE( out,                                              &
             "( ' ............... compute initial point ............... ' )" )

!  Use a direct method

        DELTA = RHS

        CALL CPU_TIME( time )
        CALL WCP_iterative_refinement(                                         &
                 dims, n, m, K, FACTORS, CNTL, A_ne, A_val, A_col, A_ptr,      &
                 zero, BARRIER_X, BARRIER_C, DELTA, factor,                    &
                 DIAG_X, ldiag_x, SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,       &
                 SOL, RES, BEST, SOL_y, BEST_y, RES_y, RES_x, itref_max,       &
                 print_level, control, inform )
  
        IF ( inform%status /= 0 ) GO TO 700
        CALL CPU_TIME( dum ) ; dum = dum - time
        X = DELTA( dims%x_s : dims%x_e )
        C = DELTA( dims%c_s : dims%c_e ) 
        Y = DELTA( dims%y_s : dims%y_e )
    
        IF ( printt ) WRITE( out, "( ' ** solve time = ', F10.2 ) " ) dum
        inform%time%solve = inform%time%solve + dum

!  Compute the residual of the linear system

        CALL WCP_residual( dims, n, m, dims%v_e, A_ne, A_val, A_col, A_ptr,    &
                           X, C, Y, RHS( dims%x_s : dims%x_e ),                &
                           RHS( dims%c_s : dims%c_e ),                         &
                           RHS( dims%y_s : dims%y_e ), HX( : dims%v_e ),       &
                           zero, BARRIER_X, BARRIER_C, SCALE_C,                &
                           errorg, errorc, print_level, control )

!  If the residual of the linear system is larger than the current 
!  optimality residual, no further progress is likely.

        IF ( SQRT( SUM( ( HX( : dims%v_e ) - RHS ) ** 2 ) ) > tenm8 ) THEN

!  It wasn't. We might have run out of options ...

          IF ( factor == 2 .AND. maxpiv ) THEN
            inform%status = - 8 ; GO TO 600 
              
!  ... or we may change the method

          ELSE IF ( factor < 2 .AND. maxpiv ) THEN
            factor = 2
            pivot_tol = control%pivot_tol
            maxpiv = pivot_tol >= half
            IF ( printi )                                                      &
              WRITE( out, "( '       Switching to augmented system method' )" )

!  Re-analyse the sparsity pattern of the required matrix

            CALL CPU_TIME( time ) 
            CALL WCP_analyse(                                                  &
                     dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C, factor,   &
                     nnzks, lk, liw, ldiag_x, ldiag_c_l, ldiag_c_u, IW,        &
                     Abycol_val, Abycol_row, Abycol_ptr, K_colptr, DIAG_X,     &
                     DIAG_C, K, FACTORS, CNTL, print_level, control, inform )
            CALL CPU_TIME( dum ) ; dum = dum - time
            IF ( printt )                                                      &
              WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
            inform%time%analyse = inform%time%analyse + dum

!  ... or we can increase the pivot tolerance

          ELSE
            maxpiv = .TRUE.
            pivot_tol = half
            IF ( printi )                                                     &
              WRITE( out, "( '       Pivot tolerance increased ' )" )
          END IF
          alpha = zero
          GO TO 200
        END IF
  
        IF ( printw ) WRITE( out,                                              &
             "( ' ........... initial point computed ........... ' )" )

        IF ( printd ) THEN 
          WRITE( out, 2120 ) ' DX ', DELTA( dims%x_s : dims%x_e )
          IF ( m > 0 ) WRITE( out, 2120 ) ' DY ', DELTA( dims%y_s : dims%y_e )
        END IF 

        Y = - Y

!  Compute the gradient of the Lagrangian function.

        CALL WCP_Lagrangian_gradient( n, m, Y, A_ne, A_val, A_col, A_ptr,      &
                                      GRAD_L( dims%x_s : dims%x_e ),           &
                                      gradient_kind, G = G )

!  Compute individual shifts

        IF ( control%initial_point == 1 ) THEN

!  shifts for the problem variables

          DO i = dims%x_free + 1, dims%x_l_start - 1
            IF ( X( i ) >= prfeas ) THEN
              PERTURB_X_l( i ) = zero
            ELSE
              PERTURB_X_l( i ) = prfeas - X( i )
            END IF
          END DO
          DO i = dims%x_l_start, dims%x_l_end
            IF ( X( i ) >= X_l( i ) + prfeas ) THEN
              PERTURB_X_l( i ) = zero
            ELSE
              PERTURB_X_l( i ) = X_l( i ) + prfeas - X( i )
            END IF
          END DO
          DO i = dims%x_u_start, dims%x_u_end
            IF ( X( i ) <= X_u( i ) - prfeas ) THEN
              PERTURB_X_u( i ) = zero
            ELSE
              PERTURB_X_u( i ) = X( i ) + prfeas - X_u( i )
            END IF
          END DO
          DO i = dims%x_u_end + 1, n
            IF ( X( i ) <= - prfeas ) THEN
              PERTURB_X_u( i ) = zero
            ELSE
              PERTURB_X_u( i ) = X( i ) + prfeas
            END IF
          END DO

!  shifts for the slack variables

          DO i = dims%c_l_start, dims%c_l_end
            IF ( C( i ) >= C_l( i ) + prfeas ) THEN
              PERTURB_C_l( i ) = zero
            ELSE
              PERTURB_C_l( i ) = C_l( i ) + prfeas - C( i )
            END IF
          END DO
          DO i = dims%c_u_start,  dims%c_u_end
            IF ( C( i ) <= C_u( i ) - prfeas ) THEN
              PERTURB_C_u( i ) = zero
            ELSE
              PERTURB_C_u( i ) = C( i ) + prfeas - C_u( i )
            END IF
          END DO

!  values and shifts for the dual problem variables

          DO i = dims%x_free + 1, dims%x_u_start - 1
            Z_l( i ) = GRAD_L( i )
            IF ( Z_l( i ) >= dufeas ) THEN
              PERTURB_Z_l( i ) = zero
            ELSE
              PERTURB_Z_l( i ) = - GRAD_L( i ) + dufeas
            END IF
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            IF ( GRAD_L( i ) >= zero ) THEN
              Z_l( i ) = GRAD_L( i ) + dufeas
              Z_u( i ) = - dufeas
            ELSE
              Z_l( i ) = dufeas 
              Z_u( i ) = GRAD_L( i ) - dufeas
            END IF
            IF ( Z_l( i ) >= dufeas ) THEN
              PERTURB_Z_l( i ) = zero
            ELSE
              PERTURB_Z_l( i ) = dufeas
            END IF
            IF ( Z_u( i ) <= - dufeas ) THEN
              PERTURB_Z_u( i ) = zero
            ELSE
              PERTURB_Z_u( i ) = dufeas
            END IF
          END DO
          DO i = dims%x_l_end + 1, n
            Z_u( i ) = GRAD_L( i )
            IF ( Z_u( i ) <= - dufeas ) THEN
              PERTURB_Z_u( i ) = zero
            ELSE
              PERTURB_Z_u( i ) = GRAD_L( i ) + dufeas
            END IF
          END DO

!  values and shifts for the dual slack variables

          DO i = dims%c_l_start, dims%c_u_start - 1
            Y_l( i ) = Y( i )
            IF ( Y_l( i ) >= dufeas ) THEN
              PERTURB_Y_l( i ) = zero
            ELSE
              PERTURB_Y_l( i ) = - Y( i ) + dufeas
            END IF
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            IF ( Y( i ) >= zero ) THEN
              Y_l( i ) = Y( i ) + dufeas
              Y_u( i ) = - dufeas
            ELSE
              Y_l( i ) = dufeas 
              Y_u( i ) = Y( i ) - dufeas
            END IF
            IF ( Y_l( i ) >= dufeas ) THEN
              PERTURB_Y_l( i ) = zero
            ELSE
              PERTURB_Y_l( i ) = dufeas
            END IF
            IF ( Y_u( i ) <= - dufeas ) THEN
              PERTURB_Y_u( i ) = zero
            ELSE
              PERTURB_Y_u( i ) = dufeas
            END IF
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            Y_u( i ) = Y( i )
            IF ( Y_u( i ) <= - dufeas ) THEN
              PERTURB_Y_u( i ) = zero
            ELSE
              PERTURB_Y_u( i ) = Y( i ) + dufeas
            END IF
          END DO

!  Compute overall shifts

        ELSE

!  shifts for the problem variables

          perturb_x = zero
          DO i = dims%x_free + 1, dims%x_l_start - 1
            IF ( X( i ) < prfeas )                                            &
              perturb_x = MAX( perturb_x, prfeas - X( i ) )
          END DO
          DO i = dims%x_l_start, dims%x_l_end
            IF ( X( i ) < X_l( i ) + prfeas )                                 &
              perturb_x = MAX( perturb_x, X_l( i ) + prfeas - X( i ) )
          END DO
          DO i = dims%x_u_start, dims%x_u_end
            IF ( X( i ) > X_u( i ) - prfeas )                                 &
              perturb_x = MAX( perturb_x, X( i ) + prfeas - X_u( i ) )
          END DO
          DO i = dims%x_u_end + 1, n
            IF ( X( i ) > - prfeas )                                          &
              perturb_x = MAX( perturb_x, X( i ) + prfeas )
          END DO

!  shifts for the slack variables

          perturb_c = zero
          DO i = dims%c_l_start, dims%c_l_end
            IF ( C( i ) < C_l( i ) + prfeas )                                 &
              perturb_c = MAX( perturb_c, C_l( i ) + prfeas - C( i ) )
          END DO
          DO i = dims%c_u_start,  dims%c_u_end
            IF ( C( i ) > C_u( i ) - prfeas )                                 &
              perturb_c = MAX( perturb_c, C( i ) + prfeas - C_u( i ) )
          END DO

!  values and shifts for the dual problem variables

          perturb_z = zero
          DO i = dims%x_free + 1, dims%x_u_start - 1
            Z_l( i ) = GRAD_L( i )
            IF ( GRAD_L( i ) < dufeas )                                       &
              perturb_z = MAX( perturb_z, dufeas - GRAD_L( i ) )
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            IF ( GRAD_L( i ) >= zero ) THEN
              Z_l( i ) = GRAD_L( i ) + dufeas
              Z_u( i ) = - dufeas
            ELSE
              Z_l( i ) = dufeas 
              Z_u( i ) = GRAD_L( i ) - dufeas
            END IF
          END DO
          DO i = dims%x_l_end + 1, n
            Z_u( i ) = GRAD_L( i )
            IF ( GRAD_L( i ) > - dufeas )                                     &
              perturb_z = MAX( perturb_z, GRAD_L( i ) + dufeas )
          END DO

!  values and shifts for the dual slack variables

          perturb_y = zero
          DO i = dims%c_l_start, dims%c_u_start - 1
            Y_l( i ) = Y( i )
            IF ( Y( i ) < dufeas )                                            &
              perturb_y = MAX( perturb_y, dufeas - Y( i ) )
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            IF ( Y( i ) >= zero ) THEN
              Y_l( i ) = Y( i ) + dufeas
              Y_u( i ) = - dufeas
            ELSE
              Y_l( i ) = dufeas 
              Y_u( i ) = Y( i ) - dufeas
            END IF
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            Y_u( i ) = Y( i )
            IF ( Y( i ) > - dufeas )                                           &
              perturb_y = MAX( perturb_y, Y( i ) + dufeas )
          END DO

!  Assign the shifts

!  Shift problem and slack variable bounds by the same amount

          IF ( control%initial_point == 2 ) THEN
            perturb = MAX( perturb_x, perturb_c )
            PERTURB_X_l = perturb ; PERTURB_X_u = perturb
            PERTURB_C_l = perturb ; PERTURB_C_u = perturb
            perturb = MAX( perturb_z, perturb_y )
            PERTURB_Z_l = perturb ; PERTURB_Z_u = perturb
            PERTURB_Y_l = perturb ; PERTURB_Y_u = perturb

!  Shift problem and slack variable bounds by their own overall shifts

          ELSE
            PERTURB_X_l = perturb_x ; PERTURB_X_u = perturb_x
            PERTURB_C_l = perturb_c ; PERTURB_C_u = perturb_c    
            PERTURB_Z_l = perturb_z ; PERTURB_Z_u = perturb_z    
            PERTURB_Y_l = perturb_y ; PERTURB_Y_u = perturb_y    
          END IF
        END IF

      ELSE

!  Set the simple bound perturbations perturb

        PERTURB_X_l = perturb ; PERTURB_X_u = perturb
        PERTURB_C_l = perturb ; PERTURB_C_u = perturb    
        PERTURB_Z_l = perturb ; PERTURB_Z_u = perturb
        PERTURB_Y_l = perturb ; PERTURB_Y_u = perturb    

!  Move the input starting point away from any bounds, 

!  The variable is a non-negativity

        DO i = dims%x_free + 1, dims%x_l_start - 1
          X( i ) = MAX( X( i ), prfeas )
          IF ( control%balance_initial_complementarity ) THEN
            Z_l( i ) = balance / ( X( i ) + PERTURB_X_l( i ) )
          ELSE
            Z_l( i ) = MAX( ABS( Z_l( i ) ), dufeas )
          END IF
          Z_u( i ) = zero
        END DO

!  The variable has just a lower bound

        DO i = dims%x_l_start, dims%x_u_start - 1
          X( i ) = MAX( X( i ), X_l( i ) + prfeas )
          IF ( control%balance_initial_complementarity ) THEN
            Z_l( i ) = balance / ( X( i ) - X_l( i ) + PERTURB_X_l( i ) )
          ELSE
            Z_l( i ) = MAX( ABS( Z_l( i ) ), dufeas )
          END IF
          Z_u( i ) = zero
        END DO

!  The variable has both lower and upper bounds

        DO i = dims%x_u_start, dims%x_l_end
          IF ( X_l( i ) + prfeas >= X_u( i ) - prfeas ) THEN 
            X( i ) = half * ( X_l( i ) + X_u( i ) ) 
          ELSE 
            X( i ) = MIN( MAX( X( i ), X_l( i ) + prfeas ), X_u( i ) - prfeas ) 
          END IF 
          IF ( control%balance_initial_complementarity ) THEN
            Z_l( i ) = balance / ( X( i ) - X_l( i ) + PERTURB_X_l( i ) )
            Z_u( i ) = - balance / ( X_u( i ) + PERTURB_X_u( i ) - X( i ) )
          ELSE
            Z_l( i ) = MAX(   ABS( Z_l( i ) ),   dufeas )  
            Z_u( i ) = MIN( - ABS( Z_u( i ) ), - dufeas )
          END IF
        END DO

!  The variable has just an upper bound

        DO i = dims%x_l_end + 1, dims%x_u_end
          X( i ) = MIN( X( i ), X_u( i ) - prfeas )
          IF ( control%balance_initial_complementarity ) THEN
            Z_u( i ) = - balance / ( X_u( i ) + PERTURB_X_u( i ) - X( i ) )
          ELSE
            Z_u( i ) = MIN( - ABS( Z_u( i ) ), - dufeas ) 
          END IF
          Z_l( i ) = zero
        END DO

!  The variable is a non-positivity

        DO i = dims%x_u_end + 1, n
          X( i ) = MIN( X( i ), - prfeas )
          IF ( control%balance_initial_complementarity ) THEN
            Z_u( i ) = balance / ( PERTURB_X_u( i ) - X( i ) )
          ELSE
            Z_u( i ) = MIN( - ABS( Z_u( i ) ), - dufeas ) 
          END IF
          Z_l( i ) = zero
        END DO

!  Compute the value of the constraint, and their residuals

        nbnds_c = 0
        IF ( m > 0 ) THEN
          C_RES( : dims%c_equality ) = - C_l( : dims%c_equality )
          C_RES( dims%c_l_start : dims%c_u_end ) = zero
          CALL WCP_AX( m, C_RES, m, A_ne, A_val, A_col, A_ptr, n, X, '+ ' )
          IF ( printd ) THEN
            WRITE( out,                                                        &
            "( /, 5X,'i', 6x, 'c', 10X, 'c_l', 9X, 'c_u', 9X, 'y_l', 9X, 'y_u')")
            DO i = 1, dims%c_l_start - 1
              WRITE( out, "( I6, 3ES12.4 )" ) i, C_RES( i ), C_l( i ), C_u( i )
            END DO
          END IF

!  The constraint has just a lower bound

          DO i = dims%c_l_start, dims%c_u_start - 1

!  Compute an appropriate scale factor

            IF ( use_scale_c ) THEN
              SCALE_C( i ) = MAX( one, ABS( C_RES( i ) ) )
            ELSE
              SCALE_C( i ) = one
            END IF

!  Scale the bounds

            C_l( i ) = C_l( i ) / SCALE_C( i )
            C( i ) = MAX( C_RES( i ) / SCALE_C( i ), C_l( i ) + prfeas ) 
            IF ( control%balance_initial_complementarity ) THEN
              Y_l( i ) = balance / ( C( i ) - C_l( i ) + PERTURB_C_l( i ) )
            ELSE
              Y_l( i ) = MAX( ABS( SCALE_C( i ) * Y_l( i ) ),  dufeas )
            END IF
          END DO

!  The constraint has both lower and upper bounds

          DO i = dims%c_u_start, dims%c_l_end

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
              C( i ) = MIN( MAX( C_RES( i ) / SCALE_C( i ),                    &
                                 C_l( i ) + prfeas), C_u( i ) - prfeas ) 
            END IF 
            IF ( control%balance_initial_complementarity ) THEN
              Y_l( i ) = balance / ( C( i ) - C_l( i ) + PERTURB_C_l( i ) )
              Y_u( i ) = - balance / ( C_u( i ) + PERTURB_C_u( i ) - C( i ) )
            ELSE
              Y_l( i ) = MAX(   ABS( SCALE_C( i ) * Y_l( i ) ),   dufeas )
              Y_u( i ) = MIN( - ABS( SCALE_C( i ) * Y_u( i ) ), - dufeas )
            END IF
          END DO

!  The constraint has just an upper bound

          DO i = dims%c_l_end + 1, dims%c_u_end

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
            IF ( control%balance_initial_complementarity ) THEN
              Y_u( i ) = - balance / ( C_u( i ) + PERTURB_C_u( i ) - C( i ) )
            ELSE
              Y_u( i ) = MIN( - ABS( SCALE_C( i ) * Y_u( i ) ), - dufeas ) 
            END IF
          END DO
        END IF
      END IF

!  ==========================
!  End of initial point phase
!  ==========================

!  Now shift the bounds and dual variables appropriately

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
        DIST_X_l( i ) = X( i ) + PERTURB_X_l( i )
        IF ( DIST_X_l( i ) < too_close ) THEN
          DIST_X_l( i ) = DIST_X_l( i ) + too_close
          PERTURB_X_l( i ) = PERTURB_X_l( i ) + too_close
        END IF
        DIST_Z_l( i ) = Z_l( i ) + PERTURB_Z_l( i )
        IF ( DIST_Z_l( i ) < too_close ) THEN
          DIST_Z_l( i ) = DIST_Z_l( i ) + too_close
          PERTURB_Z_l( i ) = PERTURB_Z_l( i ) + too_close
        END IF
        IF ( printd ) WRITE( out, "( I6, 2ES12.4, '      -     ', ES12.4,      &
       &  '      -     ' )" ) i, X( i ), zero, Z_l( i )
      END DO

!  The variable has just a lower bound

      DO i = dims%x_l_start, dims%x_u_start - 1
        nbnds_x = nbnds_x + 1
        DIST_X_l( i ) = X( i ) - X_l( i ) + PERTURB_X_l( i )
        IF ( DIST_X_l( i ) < too_close ) THEN
          DIST_X_l( i ) = DIST_X_l( i ) + too_close
          PERTURB_X_l( i ) = PERTURB_X_l( i ) + too_close
        END IF
        DIST_Z_l( i ) = Z_l( i ) + PERTURB_Z_l( i )
        IF ( DIST_Z_l( i ) < too_close ) THEN
          DIST_Z_l( i ) = DIST_Z_l( i ) + too_close
          PERTURB_Z_l( i ) = PERTURB_Z_l( i ) + too_close
        END IF
        IF ( printd ) WRITE( out, "( I6, 2ES12.4, '      -     ', ES12.4,      &
       &  '      -     ' )" ) i, X( i ), X_l( i ), Z_l( i )
      END DO

!  The variable has both lower and upper bounds

      DO i = dims%x_u_start, dims%x_l_end
        nbnds_x = nbnds_x + 2
        DIST_X_l( i ) = X( i ) - X_l( i ) + PERTURB_X_l( i ) 
        IF ( DIST_X_l( i ) < too_close ) THEN
          DIST_X_l( i ) = DIST_X_l( i ) + too_close
          PERTURB_X_l( i ) = PERTURB_X_l( i ) + too_close
        END IF
        DIST_X_u( i ) = X_u( i ) + PERTURB_X_u( i ) - X( i ) 
        IF ( DIST_X_u( i ) < too_close ) THEN
          DIST_X_u( i ) = DIST_X_u( i ) + too_close
          PERTURB_X_u( i ) = PERTURB_X_u( i ) + too_close
        END IF
        DIST_Z_l( i ) = Z_l( i ) + PERTURB_Z_l( i )
        IF ( DIST_Z_l( i ) < too_close ) THEN
          DIST_Z_l( i ) = DIST_Z_l( i ) + too_close
          PERTURB_Z_l( i ) = PERTURB_Z_l( i ) + too_close
        END IF
        DIST_Z_u( i ) = - Z_u( i ) + PERTURB_Z_u( i )
        IF ( DIST_Z_u( i ) < too_close ) THEN
          DIST_Z_u( i ) = DIST_Z_u( i ) + too_close
          PERTURB_Z_u( i ) = PERTURB_Z_u( i ) + too_close
        END IF
        IF ( printd ) WRITE( out, "( I6, 5ES12.4 )" )                          &
             i, X( i ), X_l( i ), X_u( i ), Z_l( i ), Z_u( i )
      END DO

!  The variable has just an upper bound

      DO i = dims%x_l_end + 1, dims%x_u_end
        nbnds_x = nbnds_x + 1
        DIST_X_u( i ) = X_u( i ) + PERTURB_X_u( i ) - X( i ) 
        IF ( DIST_X_u( i ) < too_close ) THEN
          DIST_X_u( i ) = DIST_X_u( i ) + too_close
          PERTURB_X_u( i ) = PERTURB_X_u( i ) + too_close
        END IF
        DIST_Z_u( i ) = - Z_u( i ) + PERTURB_Z_u( i )
        IF ( DIST_Z_u( i ) < too_close ) THEN
          DIST_Z_u( i ) = DIST_Z_u( i ) + too_close
          PERTURB_Z_u( i ) = PERTURB_Z_u( i ) + too_close
        END IF
        IF ( printd ) WRITE( out, "( I6, ES12.4, '      -     ', ES12.4,       &
       &  '      -     ', ES12.4 )" ) i, X( i ), X_u( i ), Z_u( i )
      END DO

!  The variable is a non-positivity

      DO i = dims%x_u_end + 1, n
        nbnds_x = nbnds_x + 1
        DIST_X_u( i ) = PERTURB_X_u( i ) - X( i ) 
        IF ( DIST_X_u( i ) < too_close ) THEN
          DIST_X_u( i ) = DIST_X_u( i ) + too_close
          PERTURB_X_u( i ) = PERTURB_X_u( i ) + too_close
        END IF
        DIST_Z_u( i ) = - Z_u( i ) + PERTURB_Z_u( i )
        IF ( DIST_Z_u( i ) < too_close ) THEN
          DIST_Z_u( i ) = DIST_Z_u( i ) + too_close
          PERTURB_Z_u( i ) = PERTURB_Z_u( i ) + too_close
        END IF
        IF ( printd ) WRITE( out, "( I6, ES12.4, '      -     ', ES12.4,       &
       &  '      -     ',  ES12.4 )" ) i, X( i ), zero, Z_u( i )
      END DO

!  Compute the value of the constraint, and their residuals

      nbnds_c = 0
      IF ( m > 0 ) THEN
        C_RES( : dims%c_equality ) = - C_l( : dims%c_equality )
        C_RES( dims%c_l_start : dims%c_u_end ) = zero
        CALL WCP_AX( m, C_RES, m, A_ne, A_val, A_col, A_ptr, n, X, '+ ' )
        IF ( printd ) THEN
          WRITE( out,                                                          &
          "( /, 5X,'i', 6x, 'c', 10X, 'c_l', 9X, 'c_u', 9X, 'y_l', 9X, 'y_u')")
          DO i = 1, dims%c_l_start - 1
            WRITE( out, "( I6, 3ES12.4 )" ) i, C_RES( i ), C_l( i ), C_u( i )
          END DO
        END IF

!  The constraint has just a lower bound

        DO i = dims%c_l_start, dims%c_u_start - 1
          nbnds_c = nbnds_c + 1
          C_RES( i ) = C_RES( i ) - SCALE_C( i ) * C( i )
          DIST_C_l( i ) = C( i ) - C_l( i ) + PERTURB_C_l( i ) 
          IF ( DIST_C_l( i ) < too_close ) THEN
            DIST_C_l( i ) = DIST_C_l( i ) + too_close
            PERTURB_C_l( i ) = PERTURB_C_l( i ) + too_close
          END IF
          DIST_Y_l( i ) = Y_l( i ) + PERTURB_Y_l( i ) 
          IF ( DIST_Y_l( i ) < too_close ) THEN
            DIST_Y_l( i ) = DIST_Y_l( i ) + too_close
            PERTURB_Y_l( i ) = PERTURB_Y_l( i ) + too_close
          END IF
          IF ( printd ) WRITE( out,  "( I6, 2ES12.4, '      -     ', ES12.4,   &
         &  '      -     ' )" ) i, C_RES( i ), C_l( i ), Y_l( i )
        END DO

!  The constraint has both lower and upper bounds

        DO i = dims%c_u_start, dims%c_l_end
          nbnds_c = nbnds_c + 2
          C_RES( i ) = C_RES( i ) - SCALE_C( i ) * C( i )
          DIST_C_l( i ) = C( i ) - C_l( i ) + PERTURB_C_l( i ) 
          IF ( DIST_C_l( i ) < too_close ) THEN
            DIST_C_l( i ) = DIST_C_l( i ) + too_close
            PERTURB_C_l( i ) = PERTURB_C_l( i ) + too_close
          END IF
          DIST_C_u( i ) = C_u( i ) + PERTURB_C_u( i ) - C( i ) 
          IF ( DIST_C_u( i ) < too_close ) THEN
            DIST_C_u( i ) = DIST_C_u( i ) + too_close
            PERTURB_C_u( i ) = PERTURB_C_u( i ) + too_close
          END IF
          DIST_Y_l( i ) = Y_l( i ) + PERTURB_Y_l( i ) 
          IF ( DIST_Y_l( i ) < too_close ) THEN
            DIST_Y_l( i ) = DIST_Y_l( i ) + too_close
            PERTURB_Y_l( i ) = PERTURB_Y_l( i ) + too_close
          END IF
          DIST_Y_u( i ) = - Y_u( i ) + PERTURB_Y_u( i ) 
          IF ( DIST_Y_u( i ) < too_close ) THEN
            DIST_Y_u( i ) = DIST_Y_u( i ) + too_close
            PERTURB_Y_u( i ) = PERTURB_Y_u( i ) + too_close
          END IF
          IF ( DIST_Y_u( i ) < too_close ) THEN
            DIST_Y_u( i ) = DIST_Y_u( i ) + too_close
            PERTURB_Y_u( i ) = PERTURB_Y_u( i ) + too_close
          END IF
          IF ( printd ) WRITE( out, "( I6, 5ES12.4 )" )                        &
            i, C_RES( i ), C_l( i ), C_u( i ), Y_l( i ), Y_u( i )
        END DO

!  The constraint has just an upper bound

        DO i = dims%c_l_end + 1, dims%c_u_end
          nbnds_c = nbnds_c + 1

!  Compute an appropriate initial value for the slack variable

          C_RES( i ) = C_RES( i ) - SCALE_C( i ) * C( i )
          DIST_C_u( i ) = C_u( i ) + PERTURB_C_u( i ) - C( i ) 
          IF ( DIST_C_u( i ) < too_close ) THEN
            DIST_C_u( i ) = DIST_C_u( i ) + too_close
            PERTURB_C_u( i ) = PERTURB_C_u( i ) + too_close
          END IF
          DIST_Y_u( i ) = - Y_u( i ) + PERTURB_Y_u( i ) 
          IF ( DIST_Y_u( i ) < too_close ) THEN
            DIST_Y_u( i ) = DIST_Y_u( i ) + too_close
            PERTURB_Y_u( i ) = PERTURB_Y_u( i ) + too_close
          END IF
          IF ( printd ) WRITE( out, "( I6, ES12.4, '      -     ', ES12.4,     &
         &  '      -     ', ES12.4 )" ) i, C_RES( i ), C_u( i ), Y_u( i )
        END DO
        res_prim = MAXVAL( ABS( C_RES ) )
      ELSE
        res_prim = zero
      END IF

!  Find the max-norm of the residual

      nbnds = nbnds_x + nbnds_c
      IF ( printi )                                                            &
        WRITE( out, "( /, ' >>>>> factorization package SILS used <<<<<', / )" )
      IF ( printi .AND. m > 0 .AND. dims%c_l_start <= dims%c_u_end )           &
        WRITE( out, "( ' largest/smallest scale factor ', 2ES12.4 )" )         &
          MAXVAL( SCALE_C ), MINVAL( SCALE_C )

!  Compute the complementary slackness

!     DO i = dims%x_free + 1, dims%x_l_start - 1
!       write(6,"(I6, ' x lower', 2ES12.4)" ) i, X( i ), Z_l( i )
!     END DO
!     DO i = dims%x_l_start, dims%x_l_end
!       write(6,"(I6, ' x lower', 2ES12.4)" ) i, DIST_X_l( i ), DIST_Z_l( i )
!     END DO
!     DO i = dims%x_u_start, dims%x_u_end
!       write(6,"(I6, ' x upper', 2ES12.4)" ) i, DIST_X_u( i ), DIST_Z_u( i )
!     END DO
!     DO i = dims%x_u_end + 1, n
!       write(6,"(I6, ' x upper', 2ES12.4)" ) i, X( i ), Z_u( i )
!     END DO

!     DO i = dims%c_l_start, dims%c_l_end
!       write(6,"(I6, ' c lower', 2ES12.4)" ) i, DIST_C_l( i ), DIST_Y_l( i )
!     END DO
!     DO i = dims%c_u_start, dims%c_u_end
!       write(6,"(I6, ' c upper', 2ES12.4)" ) i, DIST_C_u( i ), DIST_Y_u( i )
!     END DO

      slknes_x = DOT_PRODUCT( DIST_X_l( dims%x_free + 1 : dims%x_l_end ),      &
                              DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) ) +   &
                 DOT_PRODUCT( DIST_X_u( dims%x_u_start : n ),                  &
                              DIST_Z_u( dims%x_u_start : n ) )
      slknes_c = DOT_PRODUCT( DIST_C_l( dims%c_l_start : dims%c_l_end ),       &
                              DIST_Y_l( dims%c_l_start : dims%c_l_end ) ) +    &
                 DOT_PRODUCT( DIST_C_u( dims%c_u_start : dims%c_u_end ),       &
                              DIST_Y_u( dims%c_u_start : dims%c_u_end ) )
      slknes = slknes_x + slknes_c

      slkmin_x = MIN( MINVAL( DIST_X_l( dims%x_free + 1 : dims%x_l_end ) *     &
                              DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) ),    &
                      MINVAL( DIST_X_u( dims%x_u_start : n ) *                 &
                              DIST_Z_u( dims%x_u_start : n ) ) )
      slkmin_c = MIN( MINVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) *      &
                              DIST_Y_l( dims%c_l_start : dims%c_l_end ) ),     &
                      MINVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) *      &
                              DIST_Y_u( dims%c_u_start : dims%c_u_end ) ) )
      slkmin = MIN( slkmin_x, slkmin_c )

      slkmax_x = MAX( MAXVAL( DIST_X_l( dims%x_free + 1 : dims%x_l_end ) *     &
                              DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) ),    &
                      MAXVAL( DIST_X_u( dims%x_u_start : n ) *                 &
                              DIST_Z_u( dims%x_u_start : n ) ) )
      slkmax_c = MAX( MAXVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) *      &
                              DIST_Y_l( dims%c_l_start : dims%c_l_end ) ),     &
                      MAXVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) *      &
                              DIST_Y_u( dims%c_u_start : dims%c_u_end ) ) )
      slkmax = MAX( slkmax_x, slkmax_c )

!     WRITE(6,2120) ' >0 ', X( dims%x_free + 1 : dims%x_l_start - 1 )
!     WRITE(6,2120) ' >l ', DIST_X_l( dims%x_l_start : dims%x_l_end )
!     WRITE(6,2120) ' <u ', DIST_X_u( dims%x_u_start : dims%x_u_end )
!     WRITE(6,2120) ' <0 ', - X( dims%x_u_end + 1 : n )

      p_min = MIN( MINVAL( DIST_X_l( dims%x_free + 1 : dims%x_l_end ) ),       &
                   MINVAL( DIST_X_u( dims%x_u_start : n ) ),                   &
                   MINVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) ),        &
                   MINVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) ) )

      p_max = MAX( MAXVAL( DIST_X_l( dims%x_free + 1 : dims%x_l_end ) ),       &
                   MAXVAL( DIST_X_u( dims%x_u_start : dims%x_u_end ) ),        &
                   MAXVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) ),        &
                   MAXVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) ) )

      d_min = MIN( MINVAL( DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) ),       &
                   MINVAL( DIST_Z_u( dims%x_u_start : n ) ),                   &
                   MINVAL( DIST_Y_l( dims%c_l_start : dims%c_l_end ) ),        &
                   MINVAL( DIST_Y_u( dims%c_u_start : dims%c_u_end ) ) )

      d_max = MAX( MAXVAL( DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) ),       &
                   MAXVAL( DIST_Z_u( dims%x_u_start : n ) ),                   &
                   MAXVAL( DIST_Y_l( dims%c_l_start : dims%c_l_end ) ),        &
                   MAXVAL( DIST_Y_u( dims%c_u_start : dims%c_u_end ) ) )

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
        slknes = slknes / nbnds
      ELSE
        slknes = zero
      END IF

      IF ( printt .AND. nbnds > 0 ) WRITE( out, 2130 )                         &
        slknes, slknes_x, slknes_c, slkmin_x, slkmax_x, slkmin_c, slkmax_c,    &
        p_min, p_max, d_min, d_max

!  Compute the initial objective value

      inform%obj = f
      
      IF ( gradient_kind == 1 ) THEN
        inform%obj = inform%obj + SUM( X )
        gmax = one
      ELSE IF ( gradient_kind /= 0 ) THEN
        inform%obj = inform%obj + DOT_PRODUCT( G, X )
        gmax = MAXVAL( ABS( G ) )
      ELSE
        gmax = zero
      END IF

!  Find the largest components of A

      IF ( A_ne > 0 ) THEN
        amax = MAXVAL( ABS( A_val( : A_ne ) ) )
      ELSE
        amax = zero
      END IF

      IF ( printi ) WRITE( out, "( '  maximum element of A = ', ES12.4 )" ) amax
      IF ( printi ) WRITE( out, "( '  maximum element of g = ', ES12.4 )" ) gmax

!  Set the target barrier parameters

      IF ( control%initial_point >= 1 ) THEN
!       mu_scale = 1.0_wp
!       mu_scale = 1.1_wp
        mu_scale = 10.0_wp
        DO i = dims%x_free + 1, dims%x_l_end
          MU_X_l( i ) = mu_scale * DIST_X_l( i ) * DIST_Z_l( i )
         END DO
        DO i = dims%x_u_start, n
          MU_X_u( i ) = mu_scale * DIST_X_u( i ) * DIST_Z_u( i )
        END DO
        DO i = dims%c_l_start, dims%c_l_end 
          MU_C_l( i ) = mu_scale * DIST_C_l( i ) * DIST_Y_l( i )
        END DO
        DO i = dims%c_u_start, dims%c_u_end
          MU_C_u( i ) = mu_scale * DIST_C_u( i ) * DIST_Y_u( i )
        END DO
        mu_target = mu_scale * slkmax
      ELSE
        IF ( control%mu_target < zero ) THEN
          mu_target = sigma_max * slknes
        ELSE
          mu_target = control%mu_target
        END IF
        MU_X_l = mu_target ; MU_X_u = mu_target
        DO i = dims%x_u_start, dims%x_l_end
          mu_scale = MIN( mu_tol, MAX( one, X_u( i ) - X_l( i ) ) )
          MU_X_l( i ) = MU_X_l( i ) * mu_scale
          MU_X_u( i ) = MU_X_u( i ) * mu_scale
        END DO
        MU_C_l = mu_target ; MU_C_u = mu_target
        DO i = dims%c_u_start, dims%c_l_end
          mu_scale = MIN( mu_tol, MAX( one, C_u( i ) - C_l( i ) ) )
          MU_C_l( i ) = MU_C_l( i ) * mu_scale
          MU_C_u( i ) = MU_C_u( i ) * mu_scale
        END DO
      END IF

!  Compute the error in the complementry slackness

      slknes_req = SUM( ABS( DIST_X_l( dims%x_free + 1 : dims%x_l_end ) *      &
                             DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) -      &
                             MU_X_l( dims%x_free + 1 : dims%x_l_end ) ) )      &
                 + SUM( ABS( DIST_X_u( dims%x_u_start : n ) *                  &
                             DIST_Z_u( dims%x_u_start : n ) -                  &
                             MU_X_u( dims%x_u_start : n ) ) )                  &
                 + SUM( ABS( DIST_C_l( dims%c_l_start : dims%c_l_end ) *       &
                             DIST_Y_l( dims%c_l_start : dims%c_l_end ) -       &
                             MU_C_l( dims%c_l_start : dims%c_l_end ) ) )       &
                 + SUM( ABS( DIST_C_u( dims%c_u_start : dims%c_u_end ) *       &
                             DIST_Y_u( dims%c_u_start : dims%c_u_end ) -       &
                             MU_C_u( dims%c_u_start : dims%c_u_end ) ) )

      IF ( nbnds > 0 ) THEN
        omega_l = MIN( tenm4, tenm4 * slkmin / mu_target )
!       omega_l = MIN( tenm10, tenm10 * slkmin / mu_target )
!       omega_l = MIN( point1, point1 * slkmin / mu_target )
        omega_u = ten ** 10
!       omega_u = one / omega_l
!       write(6,"( ' omega_l, omega_u ', 2ES12.4 )" ) omega_l, omega_u
      END IF

!  Compute the gradient of the Lagrangian function.

      Y = Y_l + Y_u
      CALL WCP_Lagrangian_gradient( n, m, Y, A_ne, A_val, A_col, A_ptr,        &
                                    GRAD_L( dims%x_s : dims%x_e ),             &
                                    gradient_kind, G = G )

!  Evaluate the merit function

      merit = WCP_merit_value( dims, n, m, Y, Y_l, DIST_Y_l, Y_u, DIST_Y_u,    &
                               Z_l, DIST_Z_l, Z_u, DIST_Z_u,                   &
                               DIST_X_l, DIST_X_u, DIST_C_l,                   &
                               DIST_C_u, GRAD_L( dims%x_s : dims%x_e ),        &
                               C_RES, res_dual,                                &
                               MU_X_l, MU_X_u, MU_C_l, MU_C_u )
      res_prim_dual = SUM( ABS( C_RES ) ) + res_dual
!     old_merit = merit

!  Test to see if we are feasible

      inform%feasible =                                                        &
        res_prim <= control%stop_p .AND.res_dual <= control%stop_d
      pjgnrm = infinity

      IF ( inform%feasible ) THEN
        IF ( printi ) WRITE( out, 2070 )
        inform%x_implicit =                                                    &
               COUNT( X( dims%x_free + 1 : dims%x_l_start - 1 ) < zero )       &
             + COUNT( X_l( dims%x_l_start : dims%x_l_end ) >                   &
                      X( dims%x_l_start : dims%x_l_end ) )                     &
             + COUNT( X( dims%x_u_start : dims%x_u_end ) >                     &
                      X_u( dims%x_u_start : dims%x_u_end ) )                   &
             + COUNT( X( dims%x_u_end + 1: n ) > zero )
        inform%z_implicit =                                                    &
               COUNT( Z_l( dims%x_free + 1 : dims%x_l_end ) < zero )           &
             + COUNT( Z_u( dims%x_u_start : n ) > zero )
        inform%c_implicit =                                                    &
               COUNT( C_l( dims%c_l_start : dims%c_l_end ) >                   &
                      C_RES( dims%c_l_start : dims%c_l_end ) )                 &
             + COUNT( C_RES( dims%c_u_start : dims%c_u_end ) >                 &
                      C_u( dims%c_u_start : dims%c_u_end ) )
        inform%y_implicit =                                                    &
               COUNT( Y_l( dims%c_l_start : dims%c_l_end ) < zero )            &
             + COUNT( Y_u( dims%c_u_start : dims%c_u_end ) > zero )

        IF ( MAX( inform%x_implicit, inform%c_implicit,                        &
             inform%y_implicit, inform%z_implicit ) == 0 .AND.                 &
             control%just_feasible ) THEN
          inform%status = 0
          GO TO 500
        END IF
      END IF

!  Prepare for the major iteration

      inform%iter = 0 ; inform%nfacts = 0
      IF ( printt ) WRITE( out, "( ' ',/,' merit function value = ',           &
     &     ES12.4 )" ) merit 

      IF ( n == 0 ) THEN 
        inform%status = 0 ; GO TO 600 
      END IF 
      merit_best = merit ; it_best = 0

!  Test for convergence

!     IF ( res_prim <= control%stop_p .AND. res_dual <= control%stop_d .AND.   &
!          slknes_req <= control%stop_c ) THEN
!       inform%status = 0 ; GO TO 600 
!     END IF 

!  ===================================================
!  Analyse the sparsity pattern of the required matrix
!  ===================================================

      CALL CPU_TIME( time ) 
      CALL WCP_analyse( dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C,        &
                        factor, nnzks, lk, liw, ldiag_x, ldiag_c_l,            &
                        ldiag_c_u, IW, Abycol_val, Abycol_row, Abycol_ptr,     &
                        K_colptr, DIAG_X, DIAG_C, K, FACTORS, CNTL,            &
                        print_level, control, inform )
      CALL CPU_TIME( dum ) ; dum = dum - time
      IF ( printt ) WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
      inform%time%analyse = inform%time%analyse + dum

      refact = .TRUE. ; re = ' '  ;  co = ' ' ; al = ' '

      pivot_tol = control%pivot_tol
      maxpiv = pivot_tol >= half
      now_feasible = .FALSE.
      IF ( control%mu_accept_fraction <= zero .OR.                             &
           control%mu_accept_fraction > one ) THEN
        mu_l = one
      ELSE
        mu_l = control%mu_accept_fraction
      END IF
      mu_u = one / mu_l

      IF ( printi ) WRITE( out,                                                &
           "(  /,'  Primal    convergence tolerence = ', ES12.4,               &
          &    /,'  Dual      convergence tolerence = ', ES12.4,               &
          &    /,'  Slackness convergence tolerence = ', ES12.4 )" )           &
           control%stop_p, control%stop_d, control%stop_c

      IF (  control%perturbation_small > zero ) THEN
        perturbation_small = control%perturbation_small
      ELSE
        perturbation_small = MIN( control%stop_p, control%stop_d )
      END IF

!  ---------------------------------------------------------------------
!  ------------ Start of Perturbation Reduction Loop ------------------
!  ---------------------------------------------------------------------

      DO

         perturb_max =                                                         &
           MAX( MAXVAL( ABS( PERTURB_X_l( dims%x_free + 1 : dims%x_l_end ) ) ),&
                MAXVAL( ABS( PERTURB_X_u( dims%x_u_start : n ) ) ),            &
                MAXVAL( ABS( PERTURB_Z_l( dims%x_free + 1 : dims%x_l_end ) ) ),&
                MAXVAL( ABS( PERTURB_Z_u( dims%x_u_start : n ) ) ),            &
                MAXVAL( ABS( PERTURB_C_l( dims%c_l_start : dims%c_l_end ) ) ), &
                MAXVAL( ABS( PERTURB_C_u( dims%c_u_start : dims%c_u_end ) ) ), &
                MAXVAL( ABS( PERTURB_Y_l( dims%c_l_start : dims%c_l_end ) ) ), &
                MAXVAL( ABS( PERTURB_Y_u( dims%c_u_start : dims%c_u_end ) ) ) )
         now_interior = perturb_max <= perturbation_small

         IF ( printi )                                                         &
           WRITE( out, "( /, '  -*-*-*- perturbations (#>0,max) = ', I0, ', ', &
        &       ES10.4, ' -*-*-*-' )" )                                        &
           COUNT( PERTURB_X_l( dims%x_free + 1 : dims%x_l_end ) > zero ) +     &
           COUNT( PERTURB_X_u( dims%x_u_start : n ) > zero ) +                 &
           COUNT( PERTURB_Z_l( dims%x_free + 1 : dims%x_l_end ) > zero ) +     &
           COUNT( PERTURB_Z_u( dims%x_u_start : n ) > zero ) +                 &
           COUNT( PERTURB_C_l( dims%c_l_start : dims%c_l_end ) > zero ) +      &
           COUNT( PERTURB_C_u( dims%c_u_start : dims%c_u_end ) > zero ) +      &
           COUNT( PERTURB_Y_l( dims%c_l_start : dims%c_l_end ) > zero ) +      &
           COUNT( PERTURB_Y_u( dims%c_u_start : dims%c_u_end ) > zero ),       &
!          MIN( MINVAL( ABS( PERTURB_X_l( dims%x_free + 1 : dims%x_l_end ) ) ),&
!               MINVAL( ABS( PERTURB_X_u( dims%x_u_start : n ) ) ),            &
!               MINVAL( ABS( PERTURB_C_l( dims%c_l_start : dims%c_l_end ) ) ), &
!               MINVAL( ABS( PERTURB_C_u( dims%c_u_start : dims%c_u_end ) ) )),&
           perturb_max

           start_print = inform%iter
           start_major = .TRUE.

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

            IF ( start_major ) THEN 
              WRITE( out, 2000 ) 
              WRITE( out, 2020 ) inform%iter, re, res_prim, res_dual,          &
                                 slknes_req, merit, mu_target,                 &
                                 inform%time%total
              start_major = .FALSE.
            ELSE 
              IF ( printt .OR.                                                 &
                ( printi .AND. inform%iter == start_print ) ) WRITE( out, 2000 ) 
              coal = '  ' ; coal = TRIM( co ) // TRIM( al )
              WRITE( out, 2030 ) inform%iter, re, res_prim, res_dual,          &
                                 slknes_req, merit, alpha, coal,               &
                                 mu_target, inform%time%total
            END IF 

            IF ( printd ) THEN 
              WRITE( out, 2120 ) ' X ', X
              WRITE( out, 2120 ) ' perturb_x_l ', perturb_x_l
!             WRITE( out, 2120 ) ' DIST_x_l ', DIST_x_l
              WRITE( out, 2120 ) ' perturb_x_u ', perturb_x_u
!             WRITE( out, 2120 ) ' DIST_x_u ', DIST_x_u
              WRITE( out, 2120 ) ' Z_l ', Z_l
              WRITE( out, 2120 ) ' perturb_z_l ', perturb_z_l
!             WRITE( out, 2120 ) ' DIST_z_l ', DIST_z_l
              WRITE( out, 2120 ) ' Z_u ', Z_u
              WRITE( out, 2120 ) ' perturb_z_u ', perturb_z_u
!             WRITE( out, 2120 ) ' DIST_z_u ', DIST_z_u
            END IF 
          END IF 

!  Test for optimality

          IF ( res_prim <= control%stop_p .AND. res_dual <= control%stop_d ) THEN
            IF ( slknes_req <= control%stop_c ) THEN
              inform%status = 0 ; GO TO 490
            END IF 

!  Test for quasi-optimality          

            cs_bad = 0
            DO i = dims%x_free + 1, dims%x_l_end
              gi = DIST_X_l( i ) * DIST_Z_l( i )
              IF ( gi < mu_l * MU_X_l( i ) .OR. gi > mu_u * MU_X_l( i ) )      &
                cs_bad = cs_bad + 1
            END DO
            DO i = dims%x_u_start, n
              gi = DIST_X_u( i ) * DIST_Z_u( i )
              IF ( gi < mu_l * MU_X_u( i ) .OR. gi > mu_u * MU_X_u( i ) )      &
                cs_bad = cs_bad + 1
            END DO
            DO i = dims%c_l_start, dims%c_l_end 
              gi = DIST_C_l( i ) * DIST_Y_l( i )
              IF ( gi < mu_l * MU_C_l( i ) .OR. gi > mu_u * MU_C_l( i ) )      &
                cs_bad = cs_bad + 1
            END DO
            DO i = dims%c_u_start, dims%c_u_end
              gi = DIST_C_u( i ) * DIST_Y_u( i )
              IF ( gi < mu_l * MU_C_u( i ) .OR. gi > mu_u * MU_C_u( i ) )      &
                cs_bad = cs_bad + 1
            END DO
            IF ( cs_bad == 0 ) THEN
              WRITE( out, "( /, ' Quasi-optimal point ' )" )
              inform%status = 0 ; GO TO 490
            END IF 
          END IF 

!  Test to see if more than maxit iterations have been performed

          inform%iter = inform%iter + 1 
          IF ( inform%iter > control%maxit ) THEN 
            inform%status = - 4 ; GO TO 600 
          END IF 

!  Check that the CPU time limit has not been reached

          IF ( control%cpu_time_limit >= zero .AND.                            &
               time > control%cpu_time_limit ) THEN
            inform%status = - 9 ; GO TO 600
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

!         WRITE( 6, "( ' start, stop print, iter ', 3I8 )" )                   &
!           start_print, stop_print, inform%iter

!  Test to see whether the method has stalled

          IF ( merit <= required_infeas_reduction * merit_best ) THEN
            merit_best = merit
            it_best = 0
          ELSE
            it_best = it_best + 1
            IF ( it_best > infeas_max ) THEN
              IF ( printi ) WRITE( out, "( /, ' ============ the problem',  &
             &  ' appears to be infeasible  =============' )" )
              inform%status = - 6 ; GO TO 600 
            END IF
          END IF  

!  Compute the (Hessian) barrier terms -

!  problem variables:

          DO i = dims%x_free + 1, dims%x_u_start - 1
            IF ( ABS( DIST_X_l( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, DIST_X_l( i ), DIST_Z_l( i )
            BARRIER_X( i ) = DIST_Z_l( i ) / DIST_X_l( i ) 
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            IF ( ABS( DIST_X_l( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, DIST_X_l( i ), DIST_Z_l( i )
            IF ( ABS( DIST_X_u( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, DIST_X_u( i ), DIST_Z_u( i )
            BARRIER_X( i ) = DIST_Z_l( i ) / DIST_X_l( i ) +                   &
                             DIST_Z_u( i ) / DIST_X_u( i )
          END DO
          DO i = dims%x_l_end + 1, n
            IF ( ABS( DIST_X_u( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST X, Z ', 2ES12.4 )" )            &
                i, DIST_X_u( i ), DIST_Z_u( i )
            BARRIER_X( i ) = DIST_Z_u( i ) / DIST_X_u( i )
          END DO

!  slack variables:

          BARRIER_C( dims%c_l_start : dims%c_u_end ) = zero
          DO i = dims%c_l_start, dims%c_u_start - 1
            IF ( ABS( DIST_C_l( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST C, Y ', 2ES12.4 )" )            &
                i, DIST_C_l( i ), DIST_Y_l( i )
            BARRIER_C( i ) = DIST_Y_l( i ) / DIST_C_l( i ) 
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            IF ( ABS( DIST_C_l( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST C, Y ', 2ES12.4 )" )            &
                i, DIST_C_l( i ), DIST_Y_l( i )
            IF ( ABS( DIST_C_u( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST C, Y ', 2ES12.4 )" )            &
                i, DIST_C_u( i ), DIST_Y_u( i )
            BARRIER_C( i ) = DIST_Y_l( i ) / DIST_C_l( i ) +                   &
                             DIST_Y_u( i ) / DIST_C_u( i )
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            IF ( ABS( DIST_C_u( i ) ) <= degen_tol .AND. printw )              &
              WRITE( 6, "( ' i = ', i6, ' DIST C, Y ', 2ES12.4 )" )            &
                i, DIST_C_u( i ), DIST_Y_u( i )
            BARRIER_C( i ) = DIST_Y_u( i ) / DIST_C_u( i )
          END DO

!  =====================================================================
!  -*-*-*-*-*-*-*-*-*-*-*-*-      Factorization      -*-*-*-*-*-*-*-*-*-
!  =====================================================================

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
        
              DIAG_X( : dims%x_free ) = zero
              DIAG_X( dims%x_free + 1 : ) = BARRIER_X
              DIAG_C = BARRIER_C
        
              IF ( printd ) THEN 
                WRITE( out, "( ' DIAG_X = ', /, ( 6ES10.2 ) )" ) DIAG_X
                WRITE( out, "( ' DIAG_C = ', /, ( 6ES10.2 ) )" ) DIAG_C
              END IF

              IF ( printt ) WRITE( out, "( ' min, max DIAG_X =', 2ES11.4, /,   &
               &                           ' min, max DIAG_C =', 2ES11.4 )" )  &
                  MINVAL( ABS( DIAG_X ) ), MAXVAL( ABS( DIAG_X ) ),            &
                  MINVAL( ABS( DIAG_C ) ), MAXVAL( ABS( DIAG_C ) )

!  Form the Schur complement matrix

              IF ( printw ) WRITE( out,                                        &
                 "( ' ...... factorization of Schur complement  .......... ' )" )

              IF ( m > 0 ) THEN

                CALL WCP_form_Schur_complement(                                &
                         dims, n, m, A_ne, Abycol_val, Abycol_row, Abycol_ptr, &
                         DIAG_X, SCALE_C, DIAG_C, lk, K%val, K%row,            &
                         K_colptr, K%ne, ierr, IW( : n ), IW( n + 1 : ),       &
                         .TRUE., control%error, print_level )
!               inform%nfacts = inform%nfacts + 1 
              ELSE
                re = 'r' 
                FINFO%flag = 0 
                get_factors = .FALSE.
              END IF
              IF ( printw ) WRITE( out,                                        &
                "( ' ............... end of factorization ............... ' )" )

!  For the KKT matrix
!  ==================

            ELSE

!  Include the values of the barrier terms
 
              K%val( nnzks + 1 : nnzks + dims%x_free ) = zero
              K%val( nnzks + dims%x_free + 1 : nnzks + n ) = BARRIER_X
              K%val( nnzks + dims%c_s : nnzks + dims%c_e ) = BARRIER_C
            END IF

! ::::::::::::::::::::::::::::::
!  Factorize the required matrix
! ::::::::::::::::::::::::::::::

            IF ( get_factors ) THEN
!             WRITE( 6, "( ' K: n, nnz ', 2I4 )" ) K%n, K%ne
!             WRITE( 6, "( 3 ( 2I7, ES12.4 ) )" )                              &
!               ( K%row( i ), K%col( i ), K%val( i ), i = 1, K%ne )

!             i = 6
!             WRITE( i, "( 2I6 )" ) K%n, K%ne
!             WRITE( i, "( ( 10I6 ) )" ) K%row( : K%ne )
!             WRITE( i, "( ( 10I6 ) )" ) K%col( : K%ne ) 
!             WRITE( i, "( ( 3ES24.16 ) )" ) K%val( : K%ne ) 
!             WRITE( i, "( ( 3ES24.16 ) )" ) SOL( : K%n )
!             STOP    

              IF ( printw ) WRITE( out,                                        &
                 "( ' ......... factorization of KKT matrix ...............' )" )
              CALL SILS_factorize( K, FACTORS, CNTL, FINFO )
              IF ( printw ) WRITE( out,                                        &
                 "( ' ............... end of factorization ...............' )" )

!  Record the storage required

              inform%nfacts = inform%nfacts + 1 
              inform%factorization_integer = FINFO%nirbdu 
              inform%factorization_real = FINFO%nrlbdu 

!  Test that the factorization succeeded

              inform%factorization_status = FINFO%flag
              IF ( FINFO%flag < 0 ) THEN
                IF ( printe ) WRITE( error, 2040 ) FINFO%flag, 'SILS_factorize'

!  It didn't. We might have run out of options

                IF ( factor == 2 .AND. maxpiv ) THEN
                  inform%status = - 7 ; GO TO 700
              
!  ... or we may change the method

                ELSE IF ( factor < 2 .AND. maxpiv ) THEN
                  factor = 2
                  pivot_tol = control%pivot_tol
                  maxpiv = pivot_tol >= half
                  IF ( printi ) WRITE( out,                                    &
                    "( '       Switching to augmented system method' )" )

!  Re-analyse the sparsity pattern of the required matrix

                  CALL CPU_TIME( time ) 
                  CALL WCP_analyse(                                            &
                    dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C, factor,    &
                    nnzks,lk, liw, ldiag_x, ldiag_c_l, ldiag_c_u, IW,          &
                    Abycol_val, Abycol_row, Abycol_ptr, K_colptr, DIAG_X,      &
                    DIAG_C, K, FACTORS, CNTL, print_level, control, inform )
                  CALL CPU_TIME( dum ) ; dum = dum - time
                  IF ( printt )                                                &
                    WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
                  inform%time%analyse = inform%time%analyse + dum

!  ... or we can increase the pivot tolerance

                ELSE
                  maxpiv = .TRUE.
                  pivot_tol = half
                  IF ( printi )                                                &
                    WRITE( out, "( '       Pivot tolerance increased ' )" )
                END IF
                alpha = zero
                IF ( printi ) WRITE( out, 2000 ) 
                CYCLE

!  Record warning conditions

              ELSE IF ( FINFO%flag > 0 ) THEN
                IF ( printt ) THEN
                  WRITE( out, 2050 ) FINFO%flag, 'SILS_factorize'
                  IF ( FINFO%flag == 4 ) WRITE( out, "( ' ** Matrix has ', I7, &
                 &   ' zero eigenvalues ' )" ) K%n - finfo%rank
                END IF 
              END IF 
  
              IF ( printt ) WRITE( out,                                        &
                "( ' real/integer space used for factors ', 2I10 )" )          &
                  FINFO%nrlbdu, FINFO%nirbdu
  
            ELSE
              inform%factorization_integer = 0 
              inform%factorization_real = 0
            END IF
            CALL CPU_TIME( dum ) ; dum = dum - time
  
            IF ( printt ) THEN
              WRITE( out, "( ' ** factorize time = ', F10.2 ) " ) dum
              WRITE( out, 2060 ) inform%factorization_integer,                 &
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

          IF ( printd ) WRITE( out, 2120 ) ' GRAD_L ',                         &
            GRAD_L( dims%x_s : dims%x_e )


!       A^T u + z - g
!       u - y
! 0 =   Ax - c
!      (x + perturb_x).(z + perturb_z) - mu_x
!      (c + perturb_c).(y + perturb_y) - mu_c

! DX = diag(x + perturb_x) (etc)

! Newton

!  (      I  A^T         ) ( dx )   ( g - A^T u - z  )
!  (          I   -I     ) ( dz )   ( y - u          )
!  ( A                -I ) ( du ) = ( c - Ax         )
!  ( DZ  DX              ) ( dy )   ( mu_x e - DX DZ e )
!  (              DC  DY ) ( dc )   ( mu_c e - DC DY e )

!  remove dz and dy

! DX dz = mu_x e - DZ DX e -    DZ dx 
! dz = DX(inv) ( mu_x e- DZ DX e - DZ dx )
! dy = DC(inv) ( mu_c e - DY DC e - DY dc )

! => 

! A^T du + DX(inv) ( mu_x e - DZ DX e - DZ dx ) = g - A^T u - z 
! =>
! - A^T du + DX(inv) DZ dx = - ( g - A^T u - mu_x DX(inv) e + perturb_z )

! du - DC(inv) ( mu_c e - DY DC e - DY dc ) = y - u
! =>
!  du + DC(inv) DY dc = - ( u + mu_c DC(inv) e + perturb_y )


!  ( DX(inv) DZ            A^T ) ( dx )     (g-A^T u - mu_x DX(inv)e + perturb_z)
!  (            DC(inv) DY  -I ) ( dc ) = - (    u + mu_c DC(inv)e + perturb_y  )
!  (      A        -I          ) (-du )     (          A x - c                  )





!  Problem variables:

          RHS( : dims%x_free ) = - GRAD_L( : dims%x_free )
          DO i = dims%x_free + 1, dims%x_u_start - 1
            RHS( i ) = - GRAD_L( i ) + MU_X_l( i ) / DIST_X_l( i )             &
                       - PERTURB_Z_l( i )
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            RHS( i ) = - GRAD_L( i ) + MU_X_l( i ) / DIST_X_l( i )             &
                       - MU_X_u( i ) / DIST_X_u( i )                           &
                       - PERTURB_Z_l( i ) + PERTURB_Z_u( i )
          END DO
          DO i = dims%x_l_end + 1, n
            RHS( i ) = - GRAD_L( i ) - MU_X_u( i ) / DIST_X_u( i )             &
                       + PERTURB_Z_u( i )
          END DO

!  Slack variables:

          DO i = dims%c_l_start, dims%c_u_start - 1
            RHS( dims%c_b + i ) = - Y( i ) + MU_C_l( i ) / DIST_C_l( i )       &
                                  - PERTURB_Y_l( i )
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            RHS( dims%c_b + i ) = - Y( i ) + MU_C_l( i ) / DIST_C_l( i )       &
                                           - MU_C_u( i ) / DIST_C_u( i )       &
                                           - PERTURB_Y_l( i )  + PERTURB_Y_u( i )
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            RHS( dims%c_b + i ) = - Y( i ) - MU_C_u( i ) / DIST_C_u( i )       &
                                           + PERTURB_Y_u( i )
          END DO

!  Include the constraint infeasibilities

          RHS( dims%y_s : dims%y_e ) = - C_RES
          DELTA = RHS
!IF ( printt ) WRITE( out, "( '  c_res ', ES12.4 )" )  MAXVAL( ABS( C_RES ) )
    
          IF ( printd ) THEN 
            WRITE( out, 2120 ) ' RHS_x ', RHS( dims%x_s : dims%c_e )
            IF ( m > 0 ) WRITE( out, 2120 ) ' RHS_y ', RHS( dims%y_s : dims%y_e )
          END IF 

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Obtain the primal-dual direction for the primal variables
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!  Solve  ( H  A^T ) ( Dx^pd ) = - ( grad b )
!         ( A   0  ) ( Dy^pd )     (   r    )

          IF ( printw ) WRITE( out,                                            &
               "( ' ............... compute step  ............... ' )" )

!  Use a direct method

          CALL CPU_TIME( time )
          CALL WCP_iterative_refinement(                                       &
                   dims, n, m, K, FACTORS, CNTL, A_ne, A_val, A_col, A_ptr,    &
                   zero, BARRIER_X, BARRIER_C, DELTA, factor,                  &
                   DIAG_X, ldiag_x, SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,     &
                   SOL, RES, BEST, SOL_y, BEST_y, RES_y, RES_x, itref_max,     &
                   print_level, control, inform )
    
          IF ( inform%status /= 0 ) GO TO 700
          CALL CPU_TIME( dum ) ; dum = dum - time
    
          IF ( printd ) THEN 
            WRITE( out, 2120 ) ' SOL_x ', SOL( dims%x_s : dims%c_e )
            IF ( m > 0 ) WRITE( out, 2120 ) ' SOL_y ', SOL( dims%y_s : dims%y_e )
          END IF 

          IF ( printt ) WRITE( out, "( ' ** solve time = ', F10.2 ) " ) dum
          inform%time%solve = inform%time%solve + dum

!C_RES = zero
!CALL WCP_AX( m, C_RES, m, A_ne, A_val, A_col, A_ptr, n,  &
!             DELTA( dims%x_s : dims%x_e ), '+ ' )
!IF ( printt ) WRITE( out, "( '  A dx ', ES12.4 )" )  MAXVAL( ABS( C_RES ) )

!  Compute the residual of the linear system

          CALL WCP_residual( dims, n, m, dims%v_e, A_ne, A_val, A_col, A_ptr,  &
                             DELTA( dims%x_s : dims%x_e ),                     &
                             DELTA( dims%c_s : dims%c_e ),                     &
                             DELTA( dims%y_s : dims%y_e ),                     &
                             RHS( dims%x_s : dims%x_e ),                       &
                             RHS( dims%c_s : dims%c_e ),                       &
                             RHS( dims%y_s : dims%y_e ),                       &
                             HX( : dims%v_e ),                                 &
                             zero, BARRIER_X, BARRIER_C, SCALE_C,              &
                             errorg, errorc, print_level, control )

!  If the residual of the linear system is larger than the current 
!  optimality residual, no further progress is likely.

          IF ( SQRT( SUM( ( HX( : dims%v_e ) - RHS ) ** 2 ) ) > merit ) THEN

!  It wasn't. We might have run out of options ...

            IF ( factor == 2 .AND. maxpiv ) THEN
              inform%status = - 8 ; GO TO 600 
              
!  ... or we may change the method

            ELSE IF ( factor < 2 .AND. maxpiv ) THEN
              factor = 2
              pivot_tol = control%pivot_tol
              maxpiv = pivot_tol >= half
              IF ( printi )                                                    &
                WRITE( out, "( '       Switching to augmented system method' )" )

!  Re-analyse the sparsity pattern of the required matrix

              CALL CPU_TIME( time ) 
              CALL WCP_analyse(                                                &
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
              IF ( printi )                                                    &
                 WRITE( out, "( '       Pivot tolerance increased ' )" )
            END IF
            alpha = zero
            CYCLE
          END IF
    
          IF ( printw ) WRITE( out,                                            &
               "( ' ............... step computed ............... ' )" )

          IF ( printd ) THEN 
            WRITE( out, 2120 ) ' DX ', DELTA( dims%x_s : dims%x_e )
            IF ( m > 0 ) WRITE( out, 2120 ) ' DY ', DELTA( dims%y_s : dims%y_e )
          END IF 

!  =======
!  STEP 2:
!  =======

! ::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Obtain the search directions for the dual variables
! ::::::::::::::::::::::::::::::::::::::::::::::::::::

!  Problem variables:

!         l = 0
          DO i = dims%x_free + 1, dims%x_l_end
            DZ_l( i ) =   ( MU_X_l( i ) - DIST_Z_l( i ) *                      &
                             ( DIST_X_l( i ) + DELTA( i ) ) ) / DIST_X_l( i )
!           IF ( ABS( one + DELTA( i ) / DIST_X_l( i ) ) < 0.001 ) l = l + 1
          END DO
  
          DO i = dims%x_u_start, n
            DZ_u( i ) = - ( MU_X_u( i ) - DIST_Z_u( i ) *                      &
                             ( DIST_X_u( i ) - DELTA( i ) ) ) / DIST_X_u( i )
!           IF ( ABS( one - DELTA( i ) / DIST_X_u( i ) ) < 0.001 ) l = l + 1
          END DO
  
!         write(6,*) l, ' degenerate variable(s)'

!  Slack variables:

          DO i = dims%c_l_start, dims%c_l_end
            DY_l( i ) =   ( MU_C_l( i ) - DIST_Y_l( i ) *                      &
                      ( DIST_C_l( i ) + DELTA( dims%c_b + i ) ) ) / DIST_C_l( i )
          END DO
    
          DO i = dims%c_u_start, dims%c_u_end
            DY_u( i ) = - ( MU_C_u( i ) - DIST_Y_u( i ) *                      &
                      ( DIST_C_u( i ) - DELTA( dims%c_b + i ) ) ) / DIST_C_u( i )
          END DO
    
          IF ( printd ) THEN 
            WRITE( out, 2120 ) ' DZ_l ', DZ_l( dims%x_free + 1 : dims%x_l_end )
            WRITE( out, 2120 ) ' DZ_u ', DZ_u( dims%x_u_start : n )
          END IF 

!  Calculate the norm of the search direction

          pmax = MAX( MAXVAL( ABS( DELTA( dims%x_s : dims%x_e ) ) ),           &
                      MAXVAL( ABS( DELTA( dims%c_s : dims%c_e ) ) ),           &
                      MAXVAL( ABS( DELTA( dims%y_s : dims%y_e ) ) ),           &
                      MAXVAL( ABS( DZ_l( dims%x_free + 1 : dims%x_l_end ) ) ), &
                      MAXVAL( ABS( DZ_u( dims%x_u_start  : n ) ) ),            &
                      MAXVAL( ABS( DY_l( dims%c_l_start  : dims%c_l_end ) ) ), &
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
            DO i = dims%x_free + 1, dims%x_u_start - 1
              RHS( i ) = - DELTA( i ) * DZ_l( i ) / DIST_X_l( i )
            END DO
            DO i = dims%x_u_start, dims%x_l_end
              RHS( i ) = - DELTA( i ) * DZ_l( i ) / DIST_X_l( i )              &
                         + DELTA( i ) * DZ_u( i ) / DIST_X_u( i )
            END DO
            DO i = dims%x_l_end + 1, n
              RHS( i ) =   DELTA( i ) * DZ_u( i ) / DIST_X_u( i )
            END DO

!  Slack variables:

            DO i = dims%c_l_start, dims%c_u_start - 1
              RHS( dims%c_b + i ) =                                            &
                - DELTA( dims%c_b + i ) * DY_l( i ) / DIST_C_l( i )
            END DO
            DO i = dims%c_u_start, dims%c_l_end
              RHS( dims%c_b + i ) =                                            &
                - DELTA( dims%c_b + i ) * DY_l( i ) / DIST_C_l( i )            &
                + DELTA( dims%c_b + i ) * DY_u( i ) / DIST_C_u( i )
            END DO
            DO i = dims%c_l_end + 1, dims%c_u_end
              RHS( dims%c_b + i ) =                                            &
                  DELTA( dims%c_b + i ) * DY_u( i ) / DIST_C_u( i )
            END DO

!  Include the constraint infeasibilities

            RHS( dims%y_s : dims%y_e ) = zero
            DELTA_cor = RHS
      
            IF ( printd ) THEN 
              WRITE( out, 2120 ) ' RHS_cor_x ', RHS( dims%x_s : dims%x_e )
              IF ( m > 0 )                                                     &
                WRITE( out, 2120 ) ' RHS_cor_y ', RHS( dims%y_s : dims%y_e )
            END IF 

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Obtain the corrector direction for the primal variables
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            IF ( printw ) WRITE( out,                                          &
                 "( ' ............... compute step  ............... ' )" )

!  Use a direct method

            CALL CPU_TIME( time )
            CALL WCP_iterative_refinement(                                     &
                     dims, n, m, K, FACTORS, CNTL, A_ne, A_val, A_col, A_ptr,  &
                     zero, BARRIER_X, BARRIER_C, DELTA_cor, factor,            &
                     DIAG_X, ldiag_x, SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,   &
                     SOL, RES, BEST, SOL_y, BEST_y, RES_y, RES_x, itref_max,   &
                     print_level, control, inform )
    
            IF ( inform%status /= 0 ) GO TO 700
            CALL CPU_TIME( dum ) ; dum = dum - time
    
            IF ( printt ) WRITE( out, "( ' ** solve time = ', F10.2 ) " ) dum
            inform%time%solve = inform%time%solve + dum

!  Compute the residual of the linear system

            CALL WCP_residual( dims, n, m, dims%v_e, A_ne, A_val, A_col, A_ptr,&
                               DELTA_cor( dims%x_s : dims%x_e ),               &
                               DELTA_cor( dims%c_s : dims%c_e ),               &
                               DELTA_cor( dims%y_s : dims%y_e ),               &
                               RHS( dims%x_s : dims%x_e ),                     &
                               RHS( dims%c_s : dims%c_e ),                     &
                               RHS( dims%y_s : dims%y_e ),                     &
                               HX( : dims%v_e ), zero, BARRIER_X, BARRIER_C,   &
                               SCALE_C, errorg, errorc, print_level, control )

!  If the residual of the linear system is larger than the current 
!  optimality residual, no further progress is likely. Exit

            IF ( SQRT( SUM( ( HX( : dims%v_e ) - RHS ) ** 2 ) ) > merit ) THEN

!  It didn't. We might have run out of options

              IF ( factor == 2 .AND. maxpiv ) THEN
                inform%status = - 8 ; GO TO 600 
              
!  ... or we may change the method

              ELSE IF ( factor < 2 .AND. maxpiv ) THEN
                factor = 2
                pivot_tol = control%pivot_tol
                maxpiv = pivot_tol >= half
                IF ( printi ) WRITE( out,                                      &
                  "( ' Switching to augmented system method' )" )

!  Re-analyse the sparsity pattern of the required matrix

                CALL CPU_TIME( time ) 
                CALL WCP_analyse(                                              &
                  dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C, factor,nnzks,&
                  lk, liw, ldiag_x, ldiag_c_l, ldiag_c_u, IW, Abycol_val,      &
                  Abycol_row, Abycol_ptr, K_colptr, DIAG_X, DIAG_C, K,         &
                  FACTORS, CNTL, print_level, control, inform )
                CALL CPU_TIME( dum ) ; dum = dum - time
                IF ( printt )                                                  &
                  WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) dum
                inform%time%analyse = inform%time%analyse + dum

!  ... or we can increase the pivot tolerance

              ELSE
                maxpiv = .TRUE.
                pivot_tol = half
                IF ( printi )                                                  &
                  WRITE( out, "( '       Pivot tolerance increased ' )" )
              END IF
              alpha = zero
              CYCLE
            END IF
      
            IF ( printw ) WRITE( out,                                          &
                 "( ' ............... step computed ............... ' )" )

            IF ( printd ) THEN 
              WRITE( out, 2120 ) ' DX_cor ', DELTA_cor( dims%x_s : dims%x_e )
              IF ( m > 0 )                                                     &
                WRITE( out, 2120 ) ' DY_cor ', DELTA_cor( dims%y_s : dims%y_e )
            END IF 

!  ========
!  STEP 2b:
!  ========

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Obtain the corrector search directions for the dual variables
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!  Problem variables:

            DO i = dims%x_free + 1, dims%x_l_end
              DZ_cor_l( i ) = - ( DZ_l( i ) * DELTA( i ) +                     &
                              DIST_Z_l( i ) * DELTA_cor( i ) ) / DIST_X_l( i )
            END DO
      
            DO i = dims%x_u_start, n
              DZ_cor_u( i ) =   ( DZ_u( i ) * DELTA( i ) -                     &
                              DIST_Z_u( i ) * DELTA_cor( i ) ) / DIST_X_u( i )
            END DO
    
!  Slack variables:

            DO i = dims%c_l_start, dims%c_l_end
              DY_cor_l( i ) = - ( DY_l( i ) * DELTA( dims%c_b + i ) +          &
                DIST_Y_l( i ) * DELTA_cor( dims%c_b + i ) ) / DIST_C_l( i )
            END DO
      
            DO i = dims%c_u_start, dims%c_u_end
              DY_cor_u( i ) =   ( DY_u( i ) * DELTA( dims%c_b + i ) -          &
                DIST_Y_u( i ) * DELTA_cor( dims%c_b + i ) ) / DIST_C_u( i )
            END DO
      
            IF ( printd ) THEN 
              WRITE( out, 2120 ) ' DZ_cor_l ',                                 &
                DZ_cor_l( dims%x_free + 1 : dims%x_l_end )
              WRITE( out, 2120 ) ' DZ_cor_u ', DZ_cor_u( dims%x_u_start : n )
            END IF 

!  Calculate the norm of the search direction

            pmax_cor = MAX( MAXVAL( ABS( DELTA_cor( dims%x_s : dims%x_e ) ) ), &
              MAXVAL( ABS( DELTA_cor( dims%c_s : dims%c_e ) ) ),               &
              MAXVAL( ABS( DELTA_cor( dims%y_s : dims%y_e ) ) ),               &
              MAXVAL( ABS( DZ_cor_l( dims%x_free + 1 : dims%x_l_end ) ) ),     &
              MAXVAL( ABS( DZ_cor_u( dims%x_u_start  : n ) ) ),                &
              MAXVAL( ABS( DY_cor_l( dims%c_l_start  : dims%c_l_end ) ) ),     &
              MAXVAL( ABS( DY_cor_u( dims%c_u_start  : dims%c_u_end ) ) ) )
    
            IF ( printp ) WRITE( out, 2160 ) pmax_cor
          END IF

!  Check to see whether to use a corrector step based on the relative
!  sizes of the predictor and corrector

          IF ( control%use_corrector ) THEN
!           IF ( pmax_cor < ten * pmax ) THEN
              use_corrector = .TRUE.
!           ELSE
!             use_corrector = .FALSE.
!           END IF
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

          IF ( printw ) WRITE( out,                                            &
               "( ' .............. get steplength  .............. ' )" )

!  Form the vector H dx + A(trans) dy

          HX( dims%x_s : dims%x_e ) = zero

          CALL WCP_AX( n, HX( : n ), m, A_ne, A_val, A_col, A_ptr, m,          &
                       DELTA( dims%y_s : dims%y_e ), '+T' )

          IF ( use_corrector ) THEN

!  If a fixed point on the central path is required, perform a safeguarded 
!  global linesearch

            CALL WCP_min_piecewise_quartic( dims, n, m, nbnds, DIST_Z_l,       &
                                            DIST_Z_u, DZ_l, DZ_u, DZ_cor_l,    &
                                            DZ_cor_u, DIST_X_l, DIST_X_u,      &
                                            DELTA( dims%x_s : dims%x_e ),      &
                                            DELTA_cor( dims%x_s : dims%x_e ),  &
                                            DIST_Y_l, DIST_Y_u, DY_l, DY_u,    &
                                            DY_cor_l,                          &
                                            DY_cor_u, DIST_C_l, DIST_C_u,      &
                                            DELTA( dims%c_s : dims%c_e ),      &
                                            DELTA_cor( dims%c_s : dims%c_e ),  &
                                            res_prim_dual,                     &
                                            MU_X_l, MU_X_u, MU_C_l, MU_C_u,    &
                                            MU, omega_l, omega_u,              &
                                            alpha, alpha_max,                  &
                                            COEF4, COEF3, COEF1, COEF0,        &
                                            BREAKP, IBREAK, print_level,       &
                                            control )
            one_minus_alpha = one - alpha

!  Calculate the distances to the bounds and the dual variables at the
!  new point

            X = X + alpha * ( DELTA( dims%x_s : dims%x_e ) +                   &
                              alpha * DELTA_cor( dims%x_s : dims%x_e ) ) 
            Y = Y - alpha * ( DELTA( dims%y_s : dims%y_e ) +                   &
                              alpha * DELTA_cor( dims%y_s : dims%y_e ) )

            DO i = dims%x_free + 1, dims%x_l_end
              DIST_X_l( i ) = DIST_X_l( i ) +                                  &
                               alpha * ( DELTA( i ) + alpha * DELTA_cor( i ) )
              Z_l( i ) = Z_l( i ) + alpha * ( DZ_l( i ) + alpha * DZ_cor_l( i ) )
              DIST_Z_l( i ) = DIST_Z_l( i ) +                                  &
                               alpha * ( DZ_l( i ) + alpha * DZ_cor_l( i ) )
            END DO 
      
            DO i = dims%x_u_start, n
              DIST_X_u( i ) = DIST_X_u( i ) -                                  &
                                alpha * ( DELTA( i ) + alpha * DELTA_cor( i ) )
              Z_u( i ) = Z_u( i ) + alpha * ( DZ_u( i ) + alpha * DZ_cor_u( i ) )
              DIST_Z_u( i ) = DIST_Z_u( i ) -                                  &
                                alpha * ( DZ_u( i ) + alpha * DZ_cor_u( i ) )
            END DO 
    
!  Do the same for the slacks and their duals

            DO i = dims%c_l_start, dims%c_l_end
              DIST_C_l( i ) = DIST_C_l( i ) + alpha * ( DELTA( dims%c_b + i )  &
                               + alpha * DELTA_cor( dims%c_b + i ) )
              Y_l( i ) = Y_l( i ) + alpha * ( DY_l( i ) + alpha * DY_cor_l( i ) )
              DIST_Y_l( i ) = DIST_Y_l( i ) +                                  &
                                 alpha * ( DY_l( i ) + alpha * DY_cor_l( i ) )
            END DO 
       
            DO i = dims%c_u_start, dims%c_u_end
              DIST_C_u( i ) = DIST_C_u( i ) - alpha * ( DELTA( dims%c_b + i )  &
                              + alpha * DELTA_cor( dims%c_b + i ) )
              Y_u( i ) = Y_u( i ) + alpha * ( DY_u( i ) + alpha * DY_cor_u( i ) )
              DIST_Y_u( i ) = DIST_Y_u( i ) -                                  &
                                alpha * ( DY_u( i ) + alpha * DY_cor_u( i ) )
            END DO 
          ELSE

!  Perform a safeguarded global linesearch to find the step length

            CALL WCP_min_piecewise_quadratic( dims, n, m, nbnds, DIST_Z_l,     &
                                              DIST_Z_u, DZ_l, DZ_u, DIST_X_l,  &
                                              DIST_X_u,                        &
                                              DELTA( dims%x_s : dims%x_e ),    &
                                              DIST_Y_l, DIST_Y_u, DY_l, DY_u,  &
                                              DIST_C_l, DIST_C_u,              &
                                              DELTA( dims%c_s : dims%c_e ),    &
                                              res_prim_dual,                   &
                                              MU_X_l, MU_X_u, MU_C_l, MU_C_u,  &
                                              MU, omega_l, omega_u,            &
                                              alpha, alpha_est, alpha_max,     &
                                              COEF2, COEF1, COEF0, BREAKP,     &
                                              IBREAK, print_level, control )
            one_minus_alpha = one - alpha

!  Calculate the distances to the bounds and the dual variables at the
!  new point

            X = X + alpha * DELTA( dims%x_s : dims%x_e )
            Y = Y - alpha * DELTA( dims%y_s : dims%y_e )

            DO i = dims%x_free + 1, dims%x_l_end
              DIST_X_l( i ) = DIST_X_l( i ) + alpha * DELTA( i ) 
              Z_l( i ) = Z_l( i ) + alpha * DZ_l( i )
              DIST_Z_l( i ) = DIST_Z_l( i ) + alpha * DZ_l( i )
            END DO 
      
            DO i = dims%x_u_start, n
              DIST_X_u( i ) = DIST_X_u( i ) - alpha * DELTA( i ) 
              Z_u( i ) = Z_u( i ) + alpha * DZ_u( i )
              DIST_Z_u( i ) = DIST_Z_u( i ) - alpha * DZ_u( i )
            END DO 
    
!  Do the same for the slacks and their duals

            DO i = dims%c_l_start, dims%c_l_end
              DIST_C_l( i ) = DIST_C_l( i ) + alpha * DELTA( dims%c_b + i ) 
              Y_l( i ) = Y_l( i ) + alpha * DY_l( i )
              DIST_Y_l( i ) = DIST_Y_l( i ) + alpha * DY_l( i )
            END DO 
       
            DO i = dims%c_u_start, dims%c_u_end
              DIST_C_u( i ) = DIST_C_u( i ) - alpha * DELTA( dims%c_b + i ) 
              Y_u( i ) = Y_u( i ) + alpha * DY_u( i )
              DIST_Y_u( i ) = DIST_Y_u( i ) - alpha * DY_u( i )
            END DO 
          END IF

          IF ( alpha == alpha_max ) THEN
            al = 'b'
          ELSE
            al = ' '
          END IF

!  Update the slack variables 

          IF ( use_corrector ) THEN
            C = C + alpha * ( DELTA( dims%c_s : dims%c_e ) +                   &
                  alpha * DELTA_cor( dims%c_s : dims%c_e ) )
          ELSE
            C = C + alpha * DELTA( dims%c_s : dims%c_e )
          END IF

!  Update the values of the merit function, the gradient of the Lagrangian, 
!  and the constraint residuals

          GRAD_L( dims%x_s : dims%x_e ) = GRAD_L( dims%x_s : dims%x_e ) +      &
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

!         WRITE( 6, "(' updated  grad_l ', ES12.4 )" )                         &
!           MAXVAL( ABS( GRAD_L( dims%x_s : dims%x_e ) ) )
!         CALL WCP_Lagrangian_gradient( n, m, Y, A_ne, A_val, A_col,           &
!                                       A_ptr, GRAD_L( dims%x_s : dims%x_e ),                     &
!                                       gradient_kind, G = G )
!          WRITE( 6, "(' computed grad_l ', ES12.4 )" )                        &
!           MAXVAL( ABS( GRAD_L( dims%x_s : dims%x_e ) ) )

!  Update the norm of the constraint residual

          res_prim = one_minus_alpha * res_prim
!         nu = one_minus_alpha * nu

!  Evaluate the merit function if not already done

          merit = WCP_merit_value( dims, n, m, Y, Y_l, DIST_Y_l, Y_u,          &
                                   DIST_Y_u, Z_l, DIST_Z_l, Z_u, DIST_Z_u,     &
                                   DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u,     &
                                   GRAD_L( dims%x_s : dims%x_e ), C_RES, res_dual,                &
                                   MU_X_l, MU_X_u, MU_C_l, MU_C_u )

!         write(6,"(' res_pd ', ES12.4 )") res_prim_dual * one_minus_alpha
          res_prim_dual = SUM( ABS( C_RES ) ) + res_dual

!         write(6,"( ' res_prim , ABS(C) ', 2ES12.4 )" )                       &
!           res_prim, SUM( ABS( C_RES ) )
!         write(6,"(' res_pd,, res_cs ', ES12.4 )") res_prim_dual
          IF ( printw ) write( out,"( ' step, merit = ', 2ES12.4 )" )          &
            alpha, merit

!         min_mu = MIN( MINVAL( MU_X_l( dims%x_free + 1 : dims%x_l_end ) ),    &
!                       MINVAL( MU_X_u( dims%x_u_start : n ) ),                &
!                       MINVAL( MU_C_l( dims%c_l_start : dims%c_l_end ) ),     &
!                       MINVAL( MU_C_u( dims%c_u_start : dims%c_u_end ) ) )
!         WRITE( 6, "( ' dm, bound ', 2ES12.4 )" ) old_merit - merit,          &
!           half * omega_l * min_mu * MIN( half,                               &
!             ( one - omega_l ) * min_mu / old_merit )
!         old_merit = merit

!  Compute the objective function value

          inform%obj = f
  
          IF ( gradient_kind == 1 ) THEN
            inform%obj = inform%obj + SUM( X )
          ELSE IF ( gradient_kind /= 0 ) THEN
            inform%obj = inform%obj + DOT_PRODUCT( G, X )
          END IF

!         write(6,"( 5ES12.4 )" ) C_RES( 1 : 5 )
          IF ( m > 0 ) THEN 

            C_RES( : dims%c_equality ) = - C_l( : dims%c_equality ) 
            C_RES( dims%c_l_start : dims%c_u_end ) = - SCALE_C * C
            CALL WCP_AX( m, C_RES, m, A_ne, A_val, A_col, A_ptr, n, X, '+ ' )
            IF ( printt ) WRITE( out, "( '  Constraint residual ', ES12.4 )" )  &
              MAXVAL( ABS( C_RES ) )
!!           WRITE( 6, "( ' rec, cal cres = ', 2ES12.4 )" )                     &
!!             res_prim, MAXVAL( ABS( C_RES ) )
!            IF ( res_prim < MAXVAL( ABS( C_RES ) ) ) THEN
              res_prim = MAXVAL( ABS( C_RES ) )
!             IF ( res_prim < one .AND. itref_max < control%itref_max + 2 ) THEN
!               IF ( printt ) WRITE( out, "( ' Increasing # refinements .. ' )" )
!               itref_max = itref_max + 1
!             END IF
!           END IF
          END IF
!         write(6,"( 5ES12.4 )" ) C_RES( 1 : 5 )

!         DO i = dims%x_free + 1, dims%x_l_end
!           write(6,"(I6, ' x lower', 2ES12.4)" ) i, DIST_X_l( i ), DIST_Z_l( i )
!         END DO
!         DO i = dims%x_u_start, n
!           write(6,"(I6, ' x upper', 2ES12.4)" ) i, DIST_X_u( i ), DIST_Z_u( i )
!         END DO
!         DO i = dims%c_l_start, dims%c_l_end
!           write(6,"(I6, ' c lower', 2ES12.4)" ) i, DIST_C_l( i ), DIST_Y_l( i )
!         END DO
!         DO i = dims%c_u_start, dims%c_u_end
!           write(6,"(I6, ' c upper', 2ES12.4)" ) i, DIST_C_u( i ), DIST_Y_u( i )
!         END DO

!  Compute the complementary slackness, and the min/max components
!  of the primal/dual infeasibilities

!         IF ( res_prim <= control%stop_p .AND. res_dual <= control%stop_d     &
          IF ( res_prim <= tenm4 .AND. res_dual <= tenm3    &
               .AND. reset_mu .AND. cs_bad == 0 ) THEN
            write(6,*) ' resetting mu '
            reset_mu = .FALSE.
            DO i = dims%x_free + 1, dims%x_l_end
              MU_X_l( i ) = DIST_X_l( i ) * DIST_Z_l( i )
            END DO
            DO i = dims%x_u_start, n
              MU_X_u( i ) = DIST_X_u( i ) * DIST_Z_u( i )
            END DO
            DO i = dims%c_l_start, dims%c_l_end 
              MU_C_l( i ) = DIST_C_l( i ) * DIST_Y_l( i )
            END DO
            DO i = dims%c_u_start, dims%c_u_end
              MU_C_u( i ) = DIST_C_u( i ) * DIST_Y_u( i )
            END DO
          END IF

          IF ( printt .AND. nbnds > 0 ) THEN
            slknes_x = DOT_PRODUCT( DIST_X_l( dims%x_free + 1 : dims%x_l_end ),&
                                    DIST_Z_l( dims%x_free + 1 : dims%x_l_end ))&
                     + DOT_PRODUCT( DIST_X_u( dims%x_u_start : n ),            &
                                    DIST_Z_u( dims%x_u_start : n ) )
            slknes_c = DOT_PRODUCT( DIST_C_l( dims%c_l_start : dims%c_l_end ), &
                                    DIST_Y_l( dims%c_l_start : dims%c_l_end ) )&
                     + DOT_PRODUCT( DIST_C_u( dims%c_u_start : dims%c_u_end ), &
                                    DIST_Y_u( dims%c_u_start : dims%c_u_end ) )
            slknes = slknes_x + slknes_c
      
            slkmin_x = MIN( MINVAL( DIST_X_l( dims%x_free + 1 : dims%x_l_end )*&
                                    DIST_Z_l( dims%x_free + 1 : dims%x_l_end )),&
                            MINVAL( DIST_X_u( dims%x_u_start : n ) *           &
                                    DIST_Z_u( dims%x_u_start : n ) ) )
            slkmin_c = MIN( MINVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) *&
                                    DIST_Y_l( dims%c_l_start : dims%c_l_end )),&
                            MINVAL( DIST_C_u( dims%c_u_start : dims%c_u_end)*  &
                                    DIST_Y_u( dims%c_u_start : dims%c_u_end ) ) )
            slkmin = MIN( slkmin_x, slkmin_c )
      
            slkmax_x = MAX( MAXVAL( DIST_X_l( dims%x_free + 1 : dims%x_l_end )*&
                                    DIST_Z_l( dims%x_free + 1 : dims%x_l_end )),&
                            MAXVAL( DIST_X_u( dims%x_u_start : n ) *           &
                                    DIST_Z_u( dims%x_u_start : n ) ) )
            slkmax_c = MAX( MAXVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) *&
                                    DIST_Y_l( dims%c_l_start : dims%c_l_end )),&
                            MAXVAL( DIST_C_u( dims%c_u_start : dims%c_u_end)*  &
                                    DIST_Y_u( dims%c_u_start : dims%c_u_end ) ) )
      
            p_min = MIN( MINVAL( DIST_X_l( dims%x_free + 1 : dims%x_l_end ) ), &
                         MINVAL( DIST_X_u( dims%x_u_start : n ) ),             &
                         MINVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) ),  &
                         MINVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) ) )
      
            p_max = MAX( MAXVAL( DIST_X_l( dims%x_free + 1 : dims%x_l_end ) ), &
                         MAXVAL( DIST_X_u( dims%x_u_start : n ) ),             &
                         MAXVAL( DIST_C_l( dims%c_l_start : dims%c_l_end ) ),  &
                         MAXVAL( DIST_C_u( dims%c_u_start : dims%c_u_end ) ) )
      
            d_min = MIN( MINVAL( DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) ), &
                         MINVAL( DIST_Z_u( dims%x_u_start : n ) ),             &
                         MINVAL( DIST_Y_l( dims%c_l_start : dims%c_l_end ) ),  &
                         MINVAL( DIST_Y_u( dims%c_u_start : dims%c_u_end ) ) )
      
            d_max = MAX( MAXVAL( DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) ), &
                         MAXVAL( DIST_Z_u( dims%x_u_start : n ) ),             &
                         MAXVAL( DIST_Y_l( dims%c_l_start : dims%c_l_end ) ),  &
                         MAXVAL( DIST_Y_u( dims%c_u_start : dims%c_u_end ) ) )
      
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
            ELSE
              slknes = zero
            END IF

            WRITE( out, 2130 ) slknes, slknes_x, slknes_c, slkmin_x, slkmax_x, &
              slkmin_c, slkmax_c, p_min, p_max, d_min, d_max
          END IF

          slknes_req = SUM( ABS( DIST_X_l( dims%x_free + 1 : dims%x_l_end ) *  &
                                 DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) -  &
                                 MU_X_l( dims%x_free + 1 : dims%x_l_end ) ) )  &
                     + SUM( ABS( DIST_X_u( dims%x_u_start : n ) *              &
                                 DIST_Z_u( dims%x_u_start : n ) -              &
                                 MU_X_u( dims%x_u_start : n ) ) )              &
                     + SUM( ABS( DIST_C_l( dims%c_l_start : dims%c_l_end ) *   &
                                 DIST_Y_l( dims%c_l_start : dims%c_l_end ) -   &
                                 MU_C_l( dims%c_l_start : dims%c_l_end ) ) )   &
                     + SUM( ABS( DIST_C_u( dims%c_u_start : dims%c_u_end ) *   &
                                 DIST_Y_u( dims%c_u_start : dims%c_u_end ) -   &
                                 MU_C_u( dims%c_u_start : dims%c_u_end ) ) )

          IF ( printd ) THEN
            WRITE( out, "( ' primal-dual -vs- primal dual variables ' )" )
            WRITE( out, "( ' lower ', /, ( 2( I6, 2ES12.4 ) ) )" )             &
              ( i, DIST_Z_l( i ), MU_X_l( i ) / DIST_X_l( i ),                 &
                i =  dims%x_free + 1, dims%x_l_end ) 
            WRITE( out, "( ' upper ', /, ( 2( I6, 2ES12.4 ) ) )" )             &
              ( i, DIST_Z_u( i ), - MU_X_u( i ) / DIST_X_u( i ),               &
                i = dims%x_u_start, n )
          END IF

!  Test to see if we are feasible

          IF ( res_prim <= control%stop_p ) THEN
            IF ( control%just_feasible ) THEN
              inform%status = 0
!             inform%feasible = res_dual <= control%stop_d
              IF ( printi ) THEN
                CALL CPU_TIME( time ) ; inform%time%total = time - time_start
                WRITE( out, 2070 )
                coal = '  ' ; coal = TRIM( co ) // TRIM( al )
                WRITE( out, 2030 ) inform%iter, re, res_prim, res_dual,        &
                   slknes_req, zero, alpha, coal, mu_target, time 
                IF ( printt ) WRITE( out, 2000 ) 
              END IF
              GO TO 490
            END IF

            IF ( .NOT. inform%feasible .AND. res_dual <= control%stop_d ) THEN
              IF ( printi ) WRITE( out, 2070 )
              inform%feasible = .TRUE.
            END IF
    
!  Check to see if we are feasible

          END IF
          IF ( control%just_feasible .AND. inform%feasible .AND.               &
            MAX( inform%x_implicit, inform%c_implicit,                         &
                 inform%y_implicit, inform%z_implicit ) == 0  ) THEN
            inform%status = 0
            GO TO 500
          END IF

!  =======
!  STEP 5:
!  =======

!  Compute the projected gradient of the Lagrangian function

          pjgnrm = zero 

          DO i = 1, dims%x_free
            pjgnrm = MAX( pjgnrm, ABS(  GRAD_L( i ) ) )
          END DO

          DO i = dims%x_free + 1, dims%x_u_start - 1
            gi = GRAD_L( i ) 
            IF ( gi > zero )                                                   &
              gi = MIN( ABS( X_l( i ) - PERTURB_X_l( i ) - X( i ) ), gi ) 
            pjgnrm = MAX( pjgnrm, ABS( gi ) ) 
          END DO

          DO i = dims%x_u_start, dims%x_l_end
            gi = GRAD_L( i ) 
            IF ( gi < zero ) THEN 
              gi = - MIN( ABS( X_u( i ) + PERTURB_X_u( i ) - X( i ) ), - gi ) 
            ELSE 
              gi = MIN( ABS( X_l( i ) - PERTURB_X_l( i ) - X( i ) ), gi ) 
            END IF 
            pjgnrm = MAX( pjgnrm, ABS( gi ) ) 
          END DO

          DO i = dims%x_l_end + 1, n
            gi = GRAD_L( i ) 
            IF ( gi < zero )                                                   &
              gi = - MIN( ABS( X_u( i ) + PERTURB_X_u( i ) - X( i ) ), - gi ) 
            pjgnrm = MAX( pjgnrm, ABS( gi ) ) 
          END DO

          IF ( printd ) THEN 
            WRITE( out, 2120 ) ' DIST_X_l ', DIST_X_l
            WRITE( out, 2120 ) ' DIST_X_u ', DIST_X_u
            WRITE( out, "( ' ' )" ) 
          END IF 
  
          IF ( printd ) WRITE( out, 2110 ) pjgnrm, res_prim

        END DO 

!  ---------------------------------------------------------------------
!  ---------------------- End of Major Iteration -----------------------
!  ---------------------------------------------------------------------

  490   CONTINUE 

!  Compute the values of the constraints

        C_RES( : m ) = zero
        CALL WCP_AX( m, C_RES( : m ), m, A_ne, A_val, A_col, A_ptr, n, X, '+ ' )

!  Compute the number of, and largest, constraint violations

        inform%x_implicit =                                                    &
               COUNT( X( dims%x_free + 1 : dims%x_l_start - 1 ) < zero )       &
             + COUNT( X_l( dims%x_l_start : dims%x_l_end ) >                   &
                      X( dims%x_l_start : dims%x_l_end ) )                     &
             + COUNT( X( dims%x_u_start : dims%x_u_end ) >                     &
                      X_u( dims%x_u_start : dims%x_u_end ) )                   &
             + COUNT( X( dims%x_u_end + 1: n ) > zero )
        inform%z_implicit =                                                    &
               COUNT( Z_l( dims%x_free + 1 : dims%x_l_end ) < zero )           &
             + COUNT( Z_u( dims%x_u_start : n ) > zero )
        inform%c_implicit =                                                    &
               COUNT( C_l( dims%c_l_start : dims%c_l_end ) >                   &
                      C_RES( dims%c_l_start : dims%c_l_end ) )                 &
             + COUNT( C_RES( dims%c_u_start : dims%c_u_end ) >                 &
                      C_u( dims%c_u_start : dims%c_u_end ) )
        inform%y_implicit =                                                    &
               COUNT( Y_l( dims%c_l_start : dims%c_l_end ) < zero )            &
             + COUNT( Y_u( dims%c_u_start : dims%c_u_end ) > zero )

        max_xr = MAX( zero,                                                    &
                      MAXVAL( - X( dims%x_free + 1 : dims%x_l_start - 1 ) ),   &
                      MAXVAL( X_l( dims%x_l_start : dims%x_l_end ) -           &
                              X( dims%x_l_start : dims%x_l_end ) ),            &
                      MAXVAL( X( dims%x_u_start : dims%x_u_end ) -             &
                              X_u( dims%x_u_start : dims%x_u_end ) ),          &
                      MAXVAL( X( dims%x_u_end + 1 : n ) ) )
        max_zr = MAX( zero,                                                    &
                      MAXVAL( - Z_l( dims%x_free + 1 : dims%x_l_end ) ),       &
                      MAXVAL( Z_u( dims%x_u_start : n ) ) )
        max_cr = MAX( zero, MAXVAL( C_l( dims%c_l_start : dims%c_l_end ) -     &
                                    C_RES( dims%c_l_start : dims%c_l_end ) ),  &
                            MAXVAL( C_RES( dims%c_u_start : dims%c_u_end ) -   &
                                    C_u( dims%c_u_start : dims%c_u_end ) ) )   
        max_yr = MAX( zero, MAXVAL( - Y_l( dims%c_l_start : dims%c_l_end ) ),  &
                            MAXVAL( Y_u( dims%c_u_start : dims%c_u_end ) ) )

!       DO i = dims%x_free + 1, dims%x_l_end
!         IF ( - Z_l( i ) == max_zr ) WRITE( out, "( I0, ' lower: dual = ',    &
!        &   ES16.8, ' perturb = ', ES16.8 )" ) i, Z_l( i ), PERTURB_Z_l( i )
!         IF ( - Z_l( i ) == max_zr ) WRITE( out, "( I0, ' lower: prim = ',    &
!        &   ES16.8, ' perturb = ', ES16.8 )" ) i, X( i ) - X_l( i ),          &
!          PERTURB_X_l( i )
!       END DO
!       DO i = dims%x_u_start, n
!         IF ( Z_u( i ) == max_zr ) WRITE( out, "( I0, ' upper: dual = ',      &
!        &   ES16.8, ' perturb = ', ES16.8 )" ) i, Z_u( i ), PERTURB_Z_u( i )
!       END DO

        max_r = MAX( max_xr, max_zr, max_cr, max_yr )
!       IF ( MAX( max_xr, max_cr ) <= control%stop_p .AND.                     &
!            MAX( max_zr, max_yr ) <= control%stop_d ) THEN

        IF ( max_r == zero ) THEN
          IF ( now_interior ) THEN
            IF ( printi ) WRITE( out,                                          &
              "( /, ' =============== well-centered interior point found',     &
           &        ' ==============' )" )
            EXIT
          ELSE
            IF ( .NOT. now_feasible ) THEN
              now_feasible = .TRUE.
              IF ( printi ) WRITE( out,                                        &
                "( /, ' ====================== interior point found',          &
             &        ' =====================' )" )
            END IF
          END IF
        ELSE IF ( max_r <= control%implicit_tol .AND.                          &
                  perturb_max <= perturbation_small ) THEN
           IF ( printi ) WRITE( out,                                           &
              "( /, ' ============= feasible but not interior point found',    &
           &        ' =============' )" )
            EXIT
        END IF

        IF ( printi ) THEN 
          IF ( max_xr > zero .OR. max_zr > zero .OR.                           &
               max_yr > zero .OR. max_yr > zero ) WRITE( out, "( '' )" )
          IF ( max_xr > zero )                                                 &
            WRITE( out, "( '  # infeasible variables = ', I0,                  &
           &   ', infeasibility =', ES11.4 )" ) inform%x_implicit, max_xr
          IF ( max_zr > zero )                                                 &
            WRITE( out, "( '  # infeasible duals = ', I0,                      &
           &   ', infeasibility =', ES11.4 )" ) inform%z_implicit, max_zr
          IF ( max_cr > 0 )                                                    &
            WRITE( out, "( '  # infeasible constraints = ', I0,                &
           &   ', infeasibility =', ES11.4 )" ) inform%c_implicit, max_cr
          IF ( max_yr > 0 )                                                    &
            WRITE( out, "( '  # infeasible multipliers = ', I0,                &
           &   ', infeasibility =', ES11.4 )" ) inform%y_implicit, max_yr
        END IF

        one_minus_red_pert_fac = one - reduce_perturb_factor

        IF ( control%perturbation_strategy <= 0 ) THEN
          EXIT

!  Adjust perturb by reducing each relaxtion by the same amount

        ELSE IF ( control%perturbation_strategy == 1 .OR.                      &
                  control%perturbation_strategy == 3 ) THEN
          perturb = reduce_perturb_factor * perturb +                          &
            one_minus_red_pert_fac * max_r
!         write( 6, "( ' max r = ', ES12.4 )" ) max_r
!         write( 6, "( ' new perturb = ', ES12.4 )" ) perturb
!         write( 6, "( ' X = ', / ( 6ES12.4 ) )" ) X( : n )
!         write( 6, "( ' Z = ', / ( 6ES12.4 ) )" ) Z_l( : n )
          PERTURB_X_l = perturb ; PERTURB_X_u = perturb
          PERTURB_C_l = perturb ; PERTURB_C_u = perturb    
          PERTURB_Z_l = perturb ; PERTURB_Z_u = perturb
          PERTURB_Y_l = perturb ; PERTURB_Y_u = perturb    

!  The variable is a non-negativity

          DO i = dims%x_free + 1, dims%x_l_start - 1
            DIST_X_l( i ) = X( i ) + PERTURB_X_l( i )
            DIST_Z_l( i ) = Z_l( i ) + PERTURB_Z_l( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_l( i ), DIST_X_l( i )
          END DO

!  The variable has a lower bound

          DO i = dims%x_l_start, dims%x_l_end
            DIST_X_l( i ) = X( i ) - X_l( i ) + PERTURB_X_l( i )
            DIST_Z_l( i ) = Z_l( i ) + PERTURB_Z_l( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_l( i ), DIST_X_l( i )
          END DO

!  The variable has an upper bound

          DO i = dims%x_u_start, dims%x_u_end
            DIST_X_u( i ) = X_u( i ) + PERTURB_X_u( i ) - X( i ) 
            DIST_Z_u( i ) = - Z_u( i ) + PERTURB_Z_u( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_u( i ), DIST_X_u( i )
          END DO

!  The variable is a non-positivity

          DO i = dims%x_u_end + 1, n
            DIST_X_u( i ) = PERTURB_X_u( i ) - X( i ) 
            DIST_Z_u( i ) = - Z_u( i ) + PERTURB_Z_u( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_u( i ), DIST_X_u( i )
          END DO

!  The constraint has a lower bound

          DO i = dims%c_l_start, dims%c_l_end
            DIST_C_l( i ) = C( i ) - C_l( i ) + PERTURB_C_l( i ) 
            DIST_Y_l( i ) = Y_l( i ) + PERTURB_Y_l( i ) 
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_C_l( i ), DIST_C_l( i )
          END DO

!  The constraint has an upper bound

          DO i = dims%c_u_start, dims%c_u_end
            DIST_C_u( i ) = C_u( i ) + PERTURB_C_u( i ) - C( i ) 
            DIST_Y_u( i ) = - Y_u( i ) + PERTURB_Y_u( i ) 
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_C_u( i ), DIST_C_u( i )
          END DO

!  Adjust perturb by reducing each relaxtion as much as possible

        ELSE IF ( control%perturbation_strategy == 5 ) THEN

!  The variable is a non-negativity

          DO i = dims%x_free + 1, dims%x_l_start - 1
            PERTURB_X_l( i ) =                                                 &
              MAX( zero, reduce_perturb_factor * PERTURB_X_l( i ) -            &
                   one_minus_red_pert_fac * X( i ) )
            DIST_X_l( i ) = X( i ) + PERTURB_X_l( i )
            PERTURB_Z_l( i ) =                                                 &
              MAX( zero, reduce_perturb_factor * PERTURB_Z_l( i ) -            &
                   one_minus_red_pert_fac * Z_l( i ) )
            DIST_Z_l( i ) = Z_l( i ) + PERTURB_Z_l( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_l( i ), DIST_X_l( i )
          END DO

!  The variable has a lower bound

          DO i = dims%x_l_start, dims%x_l_end
            PERTURB_X_l( i ) =                                                 &
              MAX( zero, reduce_perturb_factor * PERTURB_X_l( i ) -            &
                   one_minus_red_pert_fac * ( X( i ) - X_l( i ) ) )
            DIST_X_l( i ) = X( i ) - X_l( i ) + PERTURB_X_l( i )
            PERTURB_Z_l( i ) =                                                 &
              MAX( zero, reduce_perturb_factor * PERTURB_Z_l( i ) -            &
                   one_minus_red_pert_fac * Z_l( i ) )
            DIST_Z_l( i ) = Z_l( i ) + PERTURB_Z_l( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_l( i ), DIST_X_l( i )
          END DO

!  The variable has an upper bound

          DO i = dims%x_u_start, dims%x_u_end
            PERTURB_X_u( i ) =                                                &
              MAX( zero, reduce_perturb_factor * PERTURB_X_u( i ) -           &
                   one_minus_red_pert_fac * ( X_u( i ) - X( i ) ) )
            DIST_X_u( i ) = X_u( i ) + PERTURB_X_u( i ) - X( i ) 
            PERTURB_Z_u( i ) =                                                &
              MAX( zero, reduce_perturb_factor * PERTURB_Z_u( i ) +           &
                   one_minus_red_pert_fac * Z_u( i ) )
            DIST_Z_u( i ) = - Z_u( i ) + PERTURB_Z_u( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_u( i ), DIST_X_u( i )
          END DO

!  The variable is a non-positivity

          DO i = dims%x_u_end + 1, n
            PERTURB_X_u( i ) =                                                &
              MAX( zero, reduce_perturb_factor * PERTURB_X_u( i ) -           &
                   one_minus_red_pert_fac * ( - X( i ) ) )
            DIST_X_u( i ) = PERTURB_X_u( i ) - X( i ) 
            PERTURB_Z_u( i ) =                                                &
              MAX( zero, reduce_perturb_factor * PERTURB_Z_u( i ) +           &
                                  one_minus_red_pert_fac * Z_u( i ) )
            DIST_Z_u( i ) = - Z_u( i ) + PERTURB_Z_u( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_u( i ), DIST_X_u( i )
          END DO

!  The constraint has a lower bound

          DO i = dims%c_l_start, dims%c_l_end
            PERTURB_C_l( i ) =                                                &
              MAX( zero, reduce_perturb_factor * PERTURB_C_l( i ) -           &
                   one_minus_red_pert_fac * ( C( i ) - C_l( i ) ) )
            DIST_C_l( i ) = C( i ) - C_l( i ) + PERTURB_C_l( i ) 
            PERTURB_Y_l( i ) =                                                &
              MAX( zero, reduce_perturb_factor * PERTURB_Y_l( i ) -           &
                   one_minus_red_pert_fac * Y_l( i ) )
            DIST_Y_l( i ) = Y_l( i ) + PERTURB_Y_l( i ) 
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_C_l( i ), DIST_C_l( i )
          END DO

!  The constraint has an upper bound

          DO i = dims%c_u_start, dims%c_u_end
            PERTURB_C_u( i ) =                                                &
               MAX( zero, reduce_perturb_factor * PERTURB_C_u( i ) -          &
                    one_minus_red_pert_fac * ( C_u( i ) - C( i ) ) )
            DIST_C_u( i ) = C_u( i ) + PERTURB_C_u( i ) - C( i ) 
            PERTURB_Y_u( i ) =                                                &
              MAX( zero, reduce_perturb_factor * PERTURB_Y_u( i ) +           &
                   one_minus_red_pert_fac * Y_u( i ) )
            DIST_Y_u( i ) = - Y_u( i ) + PERTURB_Y_u( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_C_u( i ), DIST_C_u( i )
          END DO

        ELSE

!  The variable is a non-negativity

          DO i = dims%x_free + 1, dims%x_l_start - 1
            IF ( X( i ) <= zero ) THEN
              PERTURB_X_l( i ) = reduce_perturb_factor * PERTURB_X_l( i ) -    &
                one_minus_red_pert_fac * X( i )
            ELSE
              IF ( X( i ) <= control%insufficiently_feasible ) THEN
                PERTURB_X_l( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_X_l( i )
              ELSE
                PERTURB_X_l( i ) = zero
              END IF 
            END IF 
            DIST_X_l( i ) = X( i ) + PERTURB_X_l( i )
            IF ( Z_l( i ) <= zero ) THEN
              PERTURB_Z_l( i ) = reduce_perturb_factor * PERTURB_Z_l( i )      &
                - one_minus_red_pert_fac * Z_l( i )
            ELSE
              IF ( Z_l( i ) <= control%insufficiently_feasible ) THEN
                PERTURB_Z_l( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_Z_l( i )
              ELSE
                PERTURB_Z_l( i ) = zero
              END IF 
            END IF 
            DIST_Z_l( i ) = Z_l( i ) + PERTURB_Z_l( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_l( i ), DIST_X_l( i )
          END DO

!  The variable has a lower bound

          DO i = dims%x_l_start, dims%x_l_end
            IF ( X( i ) <= X_l( i ) ) THEN
              PERTURB_X_l( i ) = reduce_perturb_factor * PERTURB_X_l( i ) -    &
                one_minus_red_pert_fac * ( X( i ) - X_l( i ) )
            ELSE
              IF ( X( i ) - X_l( i ) <= control%insufficiently_feasible ) THEN
                PERTURB_X_l( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_X_l( i )
              ELSE
                PERTURB_X_l( i ) = zero
              END IF 
            END IF 
            DIST_X_l( i ) = X( i ) - X_l( i ) + PERTURB_X_l( i )
            IF ( Z_l( i ) <= zero ) THEN
              PERTURB_Z_l( i ) = reduce_perturb_factor * PERTURB_Z_l( i )      &
                - one_minus_red_pert_fac * Z_l( i )
            ELSE
              IF ( Z_l( i ) <= control%insufficiently_feasible ) THEN
                PERTURB_Z_l( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_Z_l( i )
              ELSE
                PERTURB_Z_l( i ) = zero
              END IF 
            END IF 
            DIST_Z_l( i ) = Z_l( i ) + PERTURB_Z_l( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_l( i ), DIST_X_l( i )
          END DO

!  The variable has an upper bound

          DO i = dims%x_u_start, dims%x_u_end
            IF ( X( i ) >= X_u( i ) ) THEN
              PERTURB_X_u( i ) = reduce_perturb_factor * PERTURB_X_u( i ) -    &
                one_minus_red_pert_fac * ( X_u( i ) - X( i ) )
            ELSE
              IF ( X_u( i ) - X( i ) <= control%insufficiently_feasible ) THEN
                PERTURB_X_u( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_X_u( i )
              ELSE
                PERTURB_X_u( i ) = zero
              END IF 
            END IF
            DIST_X_u( i ) = X_u( i ) + PERTURB_X_u( i ) - X( i ) 
            IF ( Z_u( i ) >= zero ) THEN
              PERTURB_Z_u( i ) = reduce_perturb_factor * PERTURB_Z_u( i )      &
                + one_minus_red_pert_fac * Z_u( i )
            ELSE
              IF ( Z_u( i ) >= - control%insufficiently_feasible ) THEN
                PERTURB_Z_u( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_Z_u( i )
              ELSE
                PERTURB_Z_u( i ) = zero
              END IF 
            END IF
            DIST_Z_u( i ) = - Z_u( i ) + PERTURB_Z_u( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_u( i ), DIST_X_u( i )
          END DO

!  The variable is a non-positivity

          DO i = dims%x_u_end + 1, n
            IF ( X( i ) >= zero ) THEN
              PERTURB_X_u( i ) = reduce_perturb_factor * PERTURB_X_u( i ) +    &
                one_minus_red_pert_fac * X( i )
            ELSE
              IF ( X( i ) >= - control%insufficiently_feasible ) THEN
                PERTURB_X_u( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_X_u( i )
              ELSE
                PERTURB_X_u( i ) = zero
              END IF 
            END IF
            DIST_X_u( i ) = PERTURB_X_u( i ) - X( i ) 
            IF ( Z_u( i ) >= zero ) THEN
              PERTURB_Z_u( i ) = reduce_perturb_factor * PERTURB_Z_u( i )      &
                + one_minus_red_pert_fac * Z_u( i )
            ELSE
              IF ( Z_u( i ) >= - control%insufficiently_feasible ) THEN
                PERTURB_Z_u( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_Z_u( i )
              ELSE
                PERTURB_Z_u( i ) = zero
              END IF 
            END IF
            DIST_Z_u( i ) = - Z_u( i ) + PERTURB_Z_u( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_X_u( i ), DIST_X_u( i )
          END DO

!  The constraint has a lower bound

          DO i = dims%c_l_start, dims%c_l_end
            IF ( C( i ) <= C_l( i ) ) THEN
              PERTURB_C_l( i ) = reduce_perturb_factor * PERTURB_C_l( i ) -    &
                one_minus_red_pert_fac * ( C( i ) - C_l( i ) )
            ELSE
              IF ( C( i ) - C_l( i ) <= control%insufficiently_feasible ) THEN
                PERTURB_C_l( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_C_l( i )
              ELSE
                PERTURB_C_l( i ) = zero
              END IF 
            END IF 
            DIST_C_l( i ) = C( i ) - C_l( i ) + PERTURB_C_l( i ) 
            IF ( Y_l( i ) <= zero ) THEN
              PERTURB_Y_l( i ) = reduce_perturb_factor * PERTURB_Y_l( i )      &
                - one_minus_red_pert_fac * Y_l( i )
            ELSE
              IF ( Y_l( i ) <= control%insufficiently_feasible ) THEN
                PERTURB_Y_l( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_Y_l( i )
              ELSE
                PERTURB_Y_l( i ) = zero
              END IF 
            END IF 
            DIST_Y_l( i ) = Y_l( i ) + PERTURB_Y_l( i ) 
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_C_l( i ), DIST_C_l( i )
          END DO

!  The constraint has an upper bound

          DO i = dims%c_u_start, dims%c_u_end
            IF ( C( i ) >= C_u( i ) ) THEN
              PERTURB_C_u( i ) = reduce_perturb_factor * PERTURB_C_u( i ) -    &
                one_minus_red_pert_fac * ( C_u( i ) - C( i ) )
            ELSE
              IF ( C_u( i ) - C( i ) <= control%insufficiently_feasible ) THEN
                PERTURB_C_u( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_C_u( i )
              ELSE
                PERTURB_C_u( i ) = zero
              END IF 
            END IF
            DIST_C_u( i ) = C_u( i ) + PERTURB_C_u( i ) - C( i ) 
            IF ( Y_u( i ) >= zero ) THEN
              PERTURB_Y_u( i ) = reduce_perturb_factor * PERTURB_Y_u( i )      &
                + one_minus_red_pert_fac * Y_u( i )
            ELSE
              IF ( Y_u( i ) >= - control%insufficiently_feasible ) THEN
                PERTURB_Y_u( i ) =                                             &
                  control%reduce_perturb_multiplier * PERTURB_Y_u( i )
              ELSE
                PERTURB_Y_u( i ) = zero
              END IF 
            END IF
            DIST_Y_u( i ) = - Y_u( i ) + PERTURB_Y_u( i )
!           write(48,"( i8, 2es12.4 )" ) i, PERTURB_C_u( i ), DIST_C_u( i )
          END DO
        END IF

!  Adjust mu to try to reduce infeasibilities

!        IF ( control%mu_increase_factor > one ) THEN
!          mu_target = control%mu_increase_factor * mu_target
!          DO i = dims%x_free + 1, dims%x_l_start - 1
!           IF ( X( i ) <= zero .OR. Z_l( i ) <= zero )                         &
!             MU_X_l( i ) = control%mu_increase_factor * MU_X_l( i ) 
!          END DO
!          DO i = dims%x_l_start, dims%x_l_end
!           IF ( X( i ) <= X_l( i ) .OR. Z_l( i ) <= zero )                     &
!             MU_X_l( i ) = control%mu_increase_factor * MU_X_l( i )
!          END DO
!          DO i = dims%x_u_start, dims%x_u_end
!           IF ( X( i ) >= X_u( i ) .OR. Z_u( i ) >= zero )                     &
!             MU_X_u( i ) = control%mu_increase_factor * MU_X_u( i )
!          END DO
!          DO i = dims%x_u_end + 1, n
!           IF ( X( i ) >= zero .OR. Z_u( i ) >= zero )                         &
!             MU_X_u( i ) = control%mu_increase_factor * MU_X_u( i )
!          END DO
!          DO i = dims%c_l_start, dims%c_l_end
!           IF ( C( i ) <= C_l( i ) .OR. Y_l( i ) <= zero )                     &
!             MU_C_l( i ) = control%mu_increase_factor * MU_C_l( i )
!          END DO
!          DO i = dims%c_u_start, dims%c_u_end
!           IF ( C( i ) >= C_u( i ) .OR. Y_u( i ) >= zero )                     &
!             MU_C_u( i ) = control%mu_increase_factor * MU_C_u( i )
!          END DO
!        END IF

!       IF ( .TRUE. ) THEN
!         DO i = dims%x_free + 1, dims%x_l_end
!           MU_X_l( i ) = mu_target * MAX( tenm2, MIN( ten2,                  &
!                               one / MIN( DIST_X_l( i ), DIST_Z_l( i )  ) ) )
!         END DO
!         DO i = dims%x_u_start, n
!           MU_X_u( i ) = mu_target * MAX( tenm2, MIN( ten2,                  &
!                               one / MIN( DIST_X_u( i ), DIST_Z_u( i )  ) ) )
!         END DO
!         DO i = dims%c_l_start, dims%c_l_end
!           MU_C_l( i ) = mu_target * MAX( tenm2, MIN( ten2,                 &
!                               one / MIN( DIST_C_l( i ), DIST_Y_l( i )  ) ) )
!         END DO
!         DO i = dims%c_u_start, dims%c_u_end
!           MU_C_u( i ) = mu_target * MAX( tenm2, MIN( ten2,                 &
!                               one / MIN( DIST_C_u( i ), DIST_Y_u( i )  ) ) )
!         END DO
!       END IF

        IF ( control%mu_increase_factor > one ) THEN
          DO i = dims%x_free + 1, dims%x_l_end
            IF ( PERTURB_X_l( i ) > zero .OR. PERTURB_Z_l( i ) > zero )       &
              MU_X_l( i ) = control%mu_increase_factor * MU_X_l( i ) 
          END DO
          DO i = dims%x_u_start, n
            IF ( PERTURB_X_u( i ) > zero .OR. PERTURB_Z_u( i ) > zero )       &
              MU_X_u( i ) = control%mu_increase_factor * MU_X_u( i ) 
          END DO
          DO i = dims%c_l_start, dims%c_l_end
            IF ( PERTURB_C_l( i ) > zero .OR. PERTURB_C_l( i ) > zero )       &
              MU_C_l( i ) = control%mu_increase_factor * MU_C_l( i ) 
          END DO
          DO i = dims%c_u_start, dims%c_u_end
            IF ( PERTURB_C_u( i ) > zero .OR. PERTURB_C_u( i ) > zero )       &
              MU_C_u( i ) = control%mu_increase_factor * MU_C_u( i ) 
          END DO
        END IF

!  Calculate the complementarity

        slknes_req = SUM( ABS( DIST_X_l( dims%x_free + 1 : dims%x_l_end ) *    &
                               DIST_Z_l( dims%x_free + 1 : dims%x_l_end ) -    &
                               MU_X_l( dims%x_free + 1 : dims%x_l_end ) ) )    &
                   + SUM( ABS( DIST_X_u( dims%x_u_start : n ) *                &
                               DIST_Z_u( dims%x_u_start : n ) -                &
                               MU_X_u( dims%x_u_start : n ) ) )                &
                   + SUM( ABS( DIST_C_l( dims%c_l_start : dims%c_l_end ) *     &
                               DIST_Y_l( dims%c_l_start : dims%c_l_end ) -     &
                               MU_C_l( dims%c_l_start : dims%c_l_end ) ) )     &
                   + SUM( ABS( DIST_C_u( dims%c_u_start : dims%c_u_end ) *     &
                               DIST_Y_u( dims%c_u_start : dims%c_u_end ) -     &
                               MU_C_u( dims%c_u_start : dims%c_u_end ) ) )

!  Compute the constraint residuals

        C_RES( : dims%c_equality ) =                                           &
          C_RES( : dims%c_equality ) - C_l( : dims%c_equality ) 
        C_RES( dims%c_l_start : dims%c_u_end ) =                               &
          C_RES( dims%c_l_start : dims%c_u_end ) - SCALE_C * C

!  Evaluate the merit function if not already done

        merit = WCP_merit_value( dims, n, m, Y, Y_l, DIST_Y_l, Y_u,            &
                                 DIST_Y_u, Z_l, DIST_Z_l, Z_u, DIST_Z_u,       &
                                 DIST_X_l, DIST_X_u, DIST_C_l, DIST_C_u,       &
                                 GRAD_L( dims%x_s : dims%x_e ), C_RES, res_dual,                  &
                                 MU_X_l, MU_X_u, MU_C_l, MU_C_u )
!       old_merit = merit
        IF ( control%perturbation_strategy > 2 )                               &
          reduce_perturb_factor = 0.25_wp * reduce_perturb_factor

!  ---------------------------------------------------------------------
!  -------------- End of Perturbation Reduction Loop -------------------
!  ---------------------------------------------------------------------

      END DO 

  500 CONTINUE 

!  Print details of the solution obtained

  600 CONTINUE 

!  Compute the final objective function value

      inform%obj = f
      IF ( gradient_kind == 1 ) THEN
        inform%obj = inform%obj + SUM( X )
      ELSE IF ( gradient_kind /= 0 ) THEN
        inform%obj = inform%obj + DOT_PRODUCT( G, X )
      END IF

      IF ( printi ) THEN
        WRITE( out, 2010 ) inform%obj, inform%iter
        WRITE( out, 2110 ) pjgnrm, res_prim

        IF ( factor == 0 .OR. factor == 1 ) THEN
          WRITE( out, "( /, '  Schur-complement method used ' )" )
        ELSE
          WRITE( out, "( /, '  Augmented system method used ' )" )
        END IF
      END IF

!  Exit

 700  CONTINUE

!  Unscale the constraint bounds

      DO i = dims%c_l_start, dims%c_l_end
        C_l( i ) = C_l( i ) * SCALE_C( i )
      END DO

      DO i = dims%c_u_start, dims%c_u_end
        C_u( i ) = C_u( i ) * SCALE_C( i )
      END DO

!  Compute the values of the constraints

      C_RES( : m ) = zero
      CALL WCP_AX( m, C_RES( : m ), m, A_ne, A_val, A_col, A_ptr, n, X, '+ ' )

!  Compute the number of constraint violations

!     inform%x_implicit =                                                      &
!            COUNT( X( dims%x_free + 1 : dims%x_l_start - 1 )                  &
!                   < control%implicit_tol )                                   &
!          + COUNT( X( dims%x_l_start : dims%x_l_end ) -                       &
!                   X_l( dims%x_l_start : dims%x_l_end )                       &
!                   < control%implicit_tol )                                   &
!          + COUNT( X( dims%x_u_start : dims%x_u_end ) -                       &
!                   X_u( dims%x_u_start : dims%x_u_end )                       &
!                   > - control%implicit_tol )                                 &
!          + COUNT( X( dims%x_u_end + 1: n ) > - control%implicit_tol )
!     inform%z_implicit =                                                      &
!            COUNT( Z_l( dims%x_free + 1 : dims%x_l_end )                      &
!                   < control%implicit_tol )                                   &
!          + COUNT( Z_u( dims%x_u_start : n ) > - control%implicit_tol )
!     inform%c_implicit =                                                      &
!            COUNT( C_RES( dims%c_l_start : dims%c_l_end ) -                   &
!                   C_l( dims%c_l_start : dims%c_l_end )                       &
!                   < control%implicit_tol )                                   &
!          + COUNT( C_RES( dims%c_u_start : dims%c_u_end ) -                   &
!                   C_u( dims%c_u_start : dims%c_u_end )                       &
!                   > - control%implicit_tol )
!     inform%y_implicit =                                                      &
!       COUNT( Y_l( dims%c_l_start : dims%c_l_end ) < control%implicit_tol )   &
!       + COUNT( Y_u( dims%c_u_start : dims%c_u_end ) > - control%implicit_tol )

      inform%x_implicit = 0 ; inform%z_implicit = 0
      inform%c_implicit = 0 ; inform%y_implicit = 0

      DO i = 1, n 
        IF ( ABS( X_l( i ) - X_u( i ) ) < control%implicit_tol ) THEN
          IF ( control%record_x_status ) inform%X_status( i ) = 3
        ELSE IF ( X( i ) - X_l( i ) < control%implicit_tol ) THEN
          inform%x_implicit = inform%x_implicit + 1
          IF ( control%record_x_status ) inform%X_status( i ) = - 1
        ELSE IF ( X_u( i ) - X( i ) < control%implicit_tol ) THEN
          inform%x_implicit = inform%x_implicit + 1
          IF ( control%record_x_status ) inform%X_status( i ) = 1
        ELSE IF ( i < dims%x_u_start .AND.                                     &
                  ABS( Z_l( i ) ) < control%implicit_tol ) THEN
          inform%z_implicit = inform%z_implicit + 1
          IF ( control%record_x_status ) inform%X_status( i ) = - 2
        ELSE IF ( i > dims%x_l_end .AND.                                       &
                  ABS( Z_u( i ) ) < control%implicit_tol ) THEN
          inform%z_implicit = inform%z_implicit + 1
          IF ( control%record_x_status ) inform%X_status( i ) = 2
        ELSE IF ( ( i >= dims%x_u_start .AND. i <= dims%x_l_end ) .AND.        &
                  ( ABS( Z_l( i ) ) < control%implicit_tol .OR.                &
                    ABS( Z_u( i ) ) < control%implicit_tol ) ) THEN
          inform%z_implicit = inform%z_implicit + 1
          IF ( ABS( Z_l( i ) ) < control%implicit_tol .AND.                    &
               ABS( Z_u( i ) ) < control%implicit_tol )                        &
            inform%z_implicit = inform%z_implicit + 1
          IF ( control%record_x_status ) THEN
            IF ( ABS( Z_l( i ) ) < control%implicit_tol .AND.                  &
                 ABS( Z_u( i ) ) < control%implicit_tol ) THEN
              inform%X_status( i ) = - 3
            ELSE IF ( ABS( Z_l( i ) ) < control%implicit_tol ) THEN
              inform%X_status( i ) = - 2
            ELSE
              inform%X_status( i ) = 2
            END IF
          END IF
        ELSE
          IF ( control%record_x_status ) inform%X_status( i ) = 0
        END IF
      END DO 

      IF ( m > 0 ) THEN 
        DO i = 1, m 
          IF ( ABS( C_l( i ) - C_u( i ) ) < control%implicit_tol ) THEN
            IF ( control%record_c_status ) inform%C_status( i ) = 3
          ELSE IF ( C( I ) - C_l( i ) < control%implicit_tol ) THEN
            inform%c_implicit = inform%c_implicit + 1
            IF ( control%record_c_status ) inform%C_status( i ) = - 1
          ELSE IF ( C_u( I ) - C( i ) < control%implicit_tol ) THEN
            inform%c_implicit = inform%c_implicit + 1
            IF ( control%record_c_status ) inform%C_status( i ) = 1
          ELSE IF ( i < dims%c_u_start .AND.                                   &
                    ABS( Y_l( i ) ) < control%implicit_tol ) THEN
            inform%z_implicit = inform%y_implicit + 1
            IF ( control%record_c_status ) inform%C_status( i ) = - 2
          ELSE IF ( i > dims%c_l_end .AND.                                     &
                    ABS( Y_u( i ) ) < control%implicit_tol ) THEN
            inform%y_implicit = inform%y_implicit + 1
            IF ( control%record_c_status ) inform%C_status( i ) = 2
          ELSE IF ( ( i >= dims%c_u_start .AND. i <= dims%c_l_end ) .AND.      &
                    ( ABS( Y_l( i ) ) < control%implicit_tol .OR.              &
                      ABS( Y_u( i ) ) < control%implicit_tol ) ) THEN
            inform%y_implicit = inform%y_implicit + 1
            IF ( ABS( Y_l( i ) ) < control%implicit_tol .AND.                  &
                 ABS( Y_u( i ) ) < control%implicit_tol )                      &
              inform%y_implicit = inform%y_implicit + 1
            IF ( control%record_c_status ) THEN
              IF ( ABS( Y_l( i ) ) < control%implicit_tol .AND.                &
                   ABS( Y_u( i ) ) < control%implicit_tol ) THEN
                inform%C_status( i ) = - 3
              ELSE IF ( ABS( Y_l( i ) ) < control%implicit_tol ) THEN
                inform%C_status( i ) = - 2
              ELSE
                inform%C_status( i ) = 2
              END IF
            END IF
          ELSE
            IF ( control%record_c_status ) inform%C_status( i ) = 0
          END IF
        END DO 
      END IF 

      inform%feasible =                                                        &
        res_prim <= control%stop_p .AND. res_dual <= control%stop_d .AND.      &
        MAX( inform%x_implicit, inform%c_implicit,                             &
             inform%y_implicit, inform%z_implicit ) == 0

!  If required, print the numbers and maximum sizes of violations

      IF ( printi ) THEN 
        max_xr = MAX( zero,                                                    &
                      MAXVAL( - X( dims%x_free + 1 : dims%x_l_start - 1 ) ),   &
                      MAXVAL( X_l( dims%x_l_start : dims%x_l_end ) -           &
                              X( dims%x_l_start : dims%x_l_end ) ),            &
                      MAXVAL( X( dims%x_u_start : dims%x_u_end ) -             &
                              X_u( dims%x_u_start : dims%x_u_end ) ),          &
                      MAXVAL( X( dims%x_u_end + 1 : n ) ) )
        max_zr = MAX( zero,                                                    &
                      MAXVAL( - Z_l( dims%x_free + 1 : dims%x_l_end ) ),       &
                      MAXVAL( Z_u( dims%x_u_start : n ) ) )
        max_cr = MAX( zero, MAXVAL( C_l( dims%c_l_start : dims%c_l_end ) -     &
                                    C_RES( dims%c_l_start : dims%c_l_end ) ),  &
                            MAXVAL( C_RES( dims%c_u_start : dims%c_u_end ) -   &
                                    C_u( dims%c_u_start : dims%c_u_end ) ) )   
        max_yr = MAX( zero, MAXVAL( - Y_l( dims%c_l_start : dims%c_l_end ) ),  &
                            MAXVAL( Y_u( dims%c_u_start : dims%c_u_end ) ) )

        IF ( max_xr > zero .OR. max_zr > zero .OR.                             &
             max_yr > zero .OR. max_yr > zero ) WRITE( out, "( '' )" )
        IF ( max_xr > zero )                                                   &
          WRITE( out, "( '  # infeasible variables = ', I0,                    &
         &   ', infeasibility =', ES11.4 )" ) inform%x_implicit, max_xr
        IF ( max_zr > zero )                                                   &
          WRITE( out, "( '  # infeasible duals = ', I0,                        &
         &   ', infeasibility =', ES11.4 )" ) inform%z_implicit, max_zr
        IF ( max_cr > 0 )                                                      &
          WRITE( out, "( '  # infeasible constraints = ', I0,                  &
         &   ', infeasibility =', ES11.4 )" ) inform%c_implicit, max_cr
        IF ( max_yr > 0 )                                                      &
          WRITE( out, "( '  # infeasible multipliers = ', I0,                  &
         &   ', infeasibility =', ES11.4 )" ) inform%y_implicit, max_yr
      END IF

!     WRITE( 79, "( ' variables ', /, '       i  violation ' )" )
!     DO i = dims%x_free + 1, dims%x_l_start - 1
!       WRITE( 79, 2991 ) i, MAX( zero, - X( i ) )
!     END DO
!     DO i = dims%x_l_start, dims%x_l_end
!       WRITE( 79, 2991 ) i, MAX( zero, X_l( i ) - X( i ) )
!     END DO
!     DO i = dims%x_u_start, dims%x_u_end
!       WRITE( 79, 2991 ) i, MAX( zero, X( i ) - X_u( i ) )
!     END DO
!     DO i = dims%x_u_end + 1, n
!       WRITE( 79, 2991 ) i, MAX( zero, X( i ) )
!     END DO
!     WRITE( 79, "( /, ' constraints ', /, '       i  violation ' )" )
!     DO i = dims%c_l_start, dims%c_l_end
!       WRITE( 79, 2991 ) i, MAX( zero, C_l( i ) - C_RES( i ) )
!     END DO
!     DO i = dims%c_u_start, dims%c_u_end
!       WRITE( 79, 2991 ) i, MAX( zero, C_RES( i ) - C_u( i ) )
!     END DO
!2991 FORMAT( I8, ES12.4 )

!  If necessary, print warning messages

  810 CONTINUE
      IF ( printi ) then
        SELECT CASE( inform%status )
          CASE( - 1 ) ; WRITE( out,                                            &
            "( /, '  Warning - input paramters incorrect ' )" )
          CASE( - 4 ) ; WRITE( out,                                            &
            "( /, '  Warning - iteration bound exceeded ' ) " )
          CASE( - 5 ) ; WRITE( out,                                            &
            "( /, '  Warning - the constraints are inconsistent ' )" )
          CASE( - 6 ) ; WRITE( out,                                            &
            "( /, '  Warning - the constraints appear to be inconsistent ' )" )
          CASE( - 7 ) ; WRITE( out,                                            &
            "( /, '  Warning - factorization failure ' )" )
          CASE( - 8 ) ; WRITE( out,                                            &
            "( /, '  Warning - residuals too large for further progress ' )" )
        END SELECT
      END IF
      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' leaving WCP_solve_main ' )" )

      RETURN  

!  Non-executable statements

 2000 FORMAT( /,' Iter   p-feas  d-feas com-slk   merit  ',                    &
                '  step        mu       time' ) 
 2010 FORMAT( //, '  Final objective function value ', ES22.14,                &
              /,  '  Total number of iterations = ', I0 )
 2020 FORMAT( I5, A1, 3ES8.1, ES9.1, '     -      ', ES7.1, 0P, F9.2 ) 
 2030 FORMAT( I5, A1, 3ES8.1, ES9.1, ES9.2, A2, 1X, ES7.1, 0P, F9.2 )
 2040 FORMAT( '   **  Error return ', I0, ' from ', A ) 
 2050 FORMAT( '   **  Warning ', I0, ' from ', A ) 
 2060 FORMAT( I8, ' integer and ', I8, ' real words needed for factorization' )
 2070 FORMAT( /, ' ================ point satisfying equations found',         &
                 ' =============== ', / )
 2120 FORMAT( A10, 7ES10.2, /, ( 10X, 7ES10.2 ) ) 
 2110 FORMAT( /,'  Norm of projected gradient is ', ES12.4, /,                 &
                '  Norm of infeasibility is      ', ES12.4 ) 
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

!  End of WCP_solve_main

      END SUBROUTINE WCP_solve_main

!-*-*-*-*-*-*-   W C P _ T E R M I N A T E   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE WCP_terminate( data, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!      ..............................................
!      .                                            .
!      .  Deallocate internal arrays at the end     .
!      .  of the computation                        .
!      .                                            .
!      ..............................................

!  Arguments:
!
!   data    see Subroutine WCP_initialize
!   control see Subroutine WCP_initialize
!   inform  see Subroutine WCP_solve

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( WCP_data_type ), INTENT( INOUT ) :: data
      TYPE ( WCP_control_type ), INTENT( IN ) :: control        
      TYPE ( WCP_inform_type ), INTENT( INOUT ) :: inform
 
!  Local variables

      CHARACTER ( LEN = 80 ) :: array_name

!  Deallocate all arrays allocated within SILS

      CALL SILS_finalize( data%FACTORS, data%CNTL, inform%alloc_status )
      IF ( inform%alloc_status /= 0 ) THEN
        inform%status = - 3
        inform%bad_alloc = ''
        IF ( control%deallocate_error_fatal ) RETURN
      END IF

!  Deallocate all arrays allocated for the preprocessing stage

      CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )
      IF ( data%QPP_inform%status /= 0 ) THEN
        inform%status = - 3
        inform%alloc_status = data%QPP_inform%alloc_status
        inform%bad_alloc = ''
      END IF

!  Deallocate all remaining allocated arrays

      array_name = 'wcp: data%RES_x'
      CALL SPACE_dealloc_array( data%RES_x,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%IW'
      CALL SPACE_dealloc_array( data%IW,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%C_freed'
      CALL SPACE_dealloc_array( data%C_freed,                                &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%Y'
      CALL SPACE_dealloc_array( data%Y,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%SOL_y'
      CALL SPACE_dealloc_array( data%SOL_y,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%RES_y'
      CALL SPACE_dealloc_array( data%RES_y,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%BEST_y'
      CALL SPACE_dealloc_array( data%BEST_y,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%SOL'
      CALL SPACE_dealloc_array( data%SOL,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%RES'
      CALL SPACE_dealloc_array( data%RES,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%BEST'
      CALL SPACE_dealloc_array( data%BEST,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%HX'
      CALL SPACE_dealloc_array( data%HX,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%GRAD_L'
      CALL SPACE_dealloc_array( data%GRAD_L,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIST_X_l'
      CALL SPACE_dealloc_array( data%DIST_X_l,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIST_Z_l'
      CALL SPACE_dealloc_array( data%DIST_Z_l,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%Z_l'
      CALL SPACE_dealloc_array( data%Z_l,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%PERTURB_X_l'
      CALL SPACE_dealloc_array( data%PERTURB_X_l,                            &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIST_X_u'
      CALL SPACE_dealloc_array( data%DIST_X_u,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIST_Z_u'
      CALL SPACE_dealloc_array( data%DIST_Z_u,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%PERTURB_X_u'
      CALL SPACE_dealloc_array( data%PERTURB_X_u,                            &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%PERTURB_Z_u'
      CALL SPACE_dealloc_array( data%PERTURB_Z_u,                            &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%BARRIER_X'
      CALL SPACE_dealloc_array( data%BARRIER_X,                              &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DY_l'
      CALL SPACE_dealloc_array( data%DY_l,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIST_C_l'
      CALL SPACE_dealloc_array( data%DIST_C_l,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIST_Y_l'
      CALL SPACE_dealloc_array( data%DIST_Y_l,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%PERTURB_C_l'
      CALL SPACE_dealloc_array( data%PERTURB_C_l,                            &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%PERTURB_Y_l'
      CALL SPACE_dealloc_array( data%PERTURB_Y_l,                            &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DY_u'
      CALL SPACE_dealloc_array( data%DY_u,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIST_C_u'
      CALL SPACE_dealloc_array( data%DIST_C_u,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIST_Y_u'
      CALL SPACE_dealloc_array( data%DIST_Y_u,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%PERTURB_C_u'
      CALL SPACE_dealloc_array( data%PERTURB_C_u,                            &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%PERTURB_Y_u'
      CALL SPACE_dealloc_array( data%PERTURB_Y_u,                            &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%C'
      CALL SPACE_dealloc_array( data%C,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%BARRIER_C'
      CALL SPACE_dealloc_array( data%BARRIER_C,                              &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%SCALE_C'
      CALL SPACE_dealloc_array( data%SCALE_C,                                &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DELTA'
      CALL SPACE_dealloc_array( data%DELTA,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%RHS'
      CALL SPACE_dealloc_array( data%RHS,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DZ_l'
      CALL SPACE_dealloc_array( data%DZ_l,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DZ_u'
      CALL SPACE_dealloc_array( data%DZ_u,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%MU_X_l'
      CALL SPACE_dealloc_array( data%MU_X_l,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%MU_X_u'
      CALL SPACE_dealloc_array( data%MU_X_u,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%MU_C_l'
      CALL SPACE_dealloc_array( data%MU_C_l,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%MU_C_u'
      CALL SPACE_dealloc_array( data%MU_C_u,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%MU'
      CALL SPACE_dealloc_array( data%MU,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%COEF0'
      CALL SPACE_dealloc_array( data%COEF0,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%COEF1'
      CALL SPACE_dealloc_array( data%COEF1,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%COEF2'
      CALL SPACE_dealloc_array( data%COEF2,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%COEF3'
      CALL SPACE_dealloc_array( data%COEF3,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%COEF4'
      CALL SPACE_dealloc_array( data%COEF4,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%BREAKP'
      CALL SPACE_dealloc_array( data%BREAKP,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%IBREAK'
      CALL SPACE_dealloc_array( data%IBREAK,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DELTA_cor'
      CALL SPACE_dealloc_array( data%DELTA_cor,                              &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DY_cor_l'
      CALL SPACE_dealloc_array( data%DY_cor_l,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DY_cor_u'
      CALL SPACE_dealloc_array( data%DY_cor_u,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DZ_cor_l'
      CALL SPACE_dealloc_array( data%DZ_cor_l,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DZ_cor_u'
      CALL SPACE_dealloc_array( data%DZ_cor_u,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%IW'
      CALL SPACE_dealloc_array( data%IW,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%Abycol_val'
      CALL SPACE_dealloc_array( data%Abycol_val,                             &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%Abycol_row'
      CALL SPACE_dealloc_array( data%Abycol_row,                             &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%Abycol_ptr'
      CALL SPACE_dealloc_array( data%Abycol_ptr,                             &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIAG_X'
      CALL SPACE_dealloc_array( data%DIAG_X,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIAG_C'
      CALL SPACE_dealloc_array( data%DIAG_C,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%K_colptr'
      CALL SPACE_dealloc_array( data%K_colptr,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%K%row'
      CALL SPACE_dealloc_array( data%K%row,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%K%col'
      CALL SPACE_dealloc_array( data%K%col,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%K%val'
      CALL SPACE_dealloc_array( data%K%val,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIAG_X'
      CALL SPACE_dealloc_array( data%DIAG_X,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%DIAG_C'
      CALL SPACE_dealloc_array( data%DIAG_C,                                 &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: data%Index_C_freed'
      CALL SPACE_dealloc_array( data%Index_C_freed,                          &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: inform%X_status'
      CALL SPACE_dealloc_array( inform%X_status,                             &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'wcp: inform%C_status'
      CALL SPACE_dealloc_array( inform%C_status,                             &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      RETURN

!  End of subroutine WCP_terminate

      END SUBROUTINE WCP_terminate 

!-*-  W C P _ L A G R A N G I A N _ G R A D I E N T   S U B R O U T I N E  -*-

      SUBROUTINE WCP_Lagrangian_gradient( n, m, Y, A_ne, A_val, A_col, A_ptr,  &
                                          GRAD_l, gradient_kind, G )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute the gradient of the Lagrangian function 
!
!  GRAD_L = g - A(transpose) y
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( IN ) :: n, m, A_ne, gradient_kind
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: GRAD_L
      INTEGER, INTENT( IN ), DIMENSION( A_ne ) :: A_col
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( A_ne ) :: A_val
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ), OPTIONAL :: G

!  Add the product A( transpose ) y to the gradient of the quadratic

      GRAD_L = zero
      IF ( gradient_kind == 1 ) THEN
        GRAD_L = GRAD_L + one
      ELSE IF ( gradient_kind /= 0 ) THEN
        GRAD_L = GRAD_L + G
      END IF

      CALL WCP_AX( n, GRAD_L, m, A_ne, A_val, A_col, A_ptr, m, Y, '-T' )

      RETURN  

!  End of WCP_Lagrangian_gradient

      END SUBROUTINE WCP_Lagrangian_gradient

!-*-*-*-*-*-   W C P _ M E R I T _ V A L U E   F U N C T I O N   -*-*-*-*-*-*-

      FUNCTION WCP_merit_value( dims, n, m, Y, Y_l, DIST_Y_l, Y_u, DIST_Y_u,   &
                                Z_l, DIST_Z_l, Z_u, DIST_Z_u,                  &
                                DIST_X_l, DIST_X_u, DIST_C_l,                  &
                                DIST_C_u, GRAD_L, C_RES, res_dual,             &
                                MU_X_l, MU_X_u, MU_C_l, MU_C_u )
      REAL ( KIND = wp ) WCP_merit_value

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Compute the value of the merit function
!
!     | < dz_l . dx_l > - mu_x_l | + | < dz_u . dx_u > - mu_x_u | +
!     | < dy_l . dc_l > - mu_c_l | + | < dy_u . dc_u > - mu_c_u | + res_dual
!
!  where 
!
!               || GRAD_L - z_l - z_u ||
!   res_dual =  ||   y - y_l - y_u    ||
!               ||      C_RES         || 
!
!   GRAD_L = g - A(transpose) y = gradient of the Lagrangian
!   C_RES =  A x - SCALE_c * c
!   dx_l = x - x_l + omega_x_l
!   dx_u = x_u - x + omega_x_u
!   dc_l = c - c_l + omega_c_l
!   dc_u = c_u - c + omega_c_u
!   dz_l = z_l + omega_z_l
!   dz_u = - z_u + omega_z_u
!   dy_l = y_l + omega_y_l
!   dy_u = - y_u + omega_y_u
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m
      REAL ( KIND = wp ), INTENT( OUT ) :: res_dual
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: GRAD_L, Z_l, Z_u
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y_l, Y_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_free + 1 : dims%x_l_end ) :: DIST_X_l,          &
                                                            DIST_Z_l, MU_X_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%x_u_start : n ) :: DIST_X_u, DIST_Z_u, MU_X_u
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y, C_RES 
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_l_start : dims%c_l_end ) :: DIST_C_l,           &
                                                           DIST_Y_l, MU_C_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
             DIMENSION( dims%c_u_start : dims%c_u_end ) :: DIST_C_u,           &
                                                           DIST_Y_u, MU_C_u

!  Local variables

      INTEGER :: i
      REAL ( KIND = wp ) :: res_cs

!  Compute in the l_1-norm

!     DO i = 1, dims%x_free
!       write(6,*) ' f ', ABS( GRAD_L( : dims%x_free ) )
!     END DO
      res_dual = SUM( ABS( GRAD_L( : dims%x_free ) ) ) ; res_cs = zero
      
!  Problem variables:

      DO i = dims%x_free + 1, dims%x_u_start - 1
!       write(6,*) ' l ', ABS( GRAD_L( i ) - Z_l( i ) )
        res_dual = res_dual + ABS( GRAD_L( i ) - Z_l( i ) )
        res_cs = res_cs + ABS( DIST_Z_l( i ) * DIST_X_l( i ) - MU_X_l( i ) )
      END DO
      DO i = dims%x_u_start, dims%x_l_end
!       write(6,*) 'lu', ABS( GRAD_L( i ) - Z_l( i ) - Z_u( i ) )
        res_dual = res_dual + ABS( GRAD_L( i ) - Z_l( i ) - Z_u( i ) )
        res_cs = res_cs + ABS( DIST_Z_l( i ) * DIST_X_l( i ) - MU_X_l( i ) )   &
                        + ABS( DIST_Z_u( i ) * DIST_X_u( i ) - MU_X_u( i ) )
      END DO
      DO i = dims%x_l_end + 1, n
!       write(6,*) ' u ', ABS( GRAD_L( i ) - Z_u( i ) )
        res_dual = res_dual + ABS( GRAD_L( i ) - Z_u( i ) )
        res_cs = res_cs + ABS( DIST_Z_u( i ) * DIST_X_u( i ) - MU_X_u( i ) )
      END DO

!  Slack variables:

      DO i = dims%c_l_start, dims%c_u_start - 1
!       write(6,*) 'dl', ABS( Y( i ) - Y_l( i ) )
        res_dual = res_dual + ABS( Y( i ) - Y_l( i ) )
        res_cs = res_cs + ABS( DIST_Y_l( i ) * DIST_C_l( i ) - MU_C_l( i ) )
      END DO
      DO i = dims%c_u_start, dims%c_l_end
!       write(6,*) 'dlu', ABS( Y( i ) - Y_l( i ) - Y_u( i ) )
        res_dual = res_dual + ABS( Y( i ) - Y_l( i ) - Y_u( i ) )
        res_cs = res_cs + ABS( DIST_Y_l( i ) * DIST_C_l( i ) - MU_C_l( i ) )   &
                        + ABS( DIST_Y_u( i ) * DIST_C_u( i ) - MU_C_u( i ) )
      END DO
      DO i = dims%c_l_end + 1, dims%c_u_end
!       write(6,*) 'du', ABS( Y( i ) - Y_u( i ) )
        res_dual = res_dual + ABS( Y( i ) - Y_u( i ) )
        res_cs = res_cs + ABS( DIST_Y_u( i ) * DIST_C_u( i ) - MU_C_u( i ) )
      END DO

      WCP_merit_value = SUM( ABS( C_RES ) ) + res_dual + res_cs
!     write(6,"(' res_cs ', ES12.4 )") res_cs

      RETURN  

!  End of function WCP_merit_value

      END FUNCTION WCP_merit_value
 
!-*-  W C P _ I T E R A T I V E _ R E F I N E M E N T  S U B R O U T I N E -*-

      SUBROUTINE WCP_iterative_refinement(                                     &
                 dims, n, m, K, FACTORS, CNTL, A_ne, A_val, A_col, A_ptr,      &
                 barrier_free, BARRIER_X, BARRIER_C, RHS, factor,              &
                 DIAG_X, ldiag_x, SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,       &
                 SOL, RES, BEST, SOL_y, BEST_y, RES_y, RES_x,                  &
                 itref_max, print_level, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Solve the block linear system

!     ( DIAG_X           A(trans)  )
!     (         DIAG_C  - SCALE_C  ) ( sol ) = ( rhs )
!     (   A   - SCALE_C            )

!  using iterative refinement, returning sol in rhs

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, A_ne, ldiag_x, ldiag_c_l, ldiag_c_u
      INTEGER, INTENT( IN ) :: itref_max, factor, print_level
      REAL ( KIND = wp ), INTENT( IN ) :: barrier_free
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
      TYPE ( SMT_type ), INTENT( INOUT ) :: K
      TYPE ( SILS_factors ), INTENT( INOUT ) :: FACTORS
      TYPE ( SILS_control ), INTENT( INOUT ) :: CNTL
      TYPE ( WCP_control_type ), INTENT( IN ) :: control        
      TYPE ( WCP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: it, itref_max_block
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
      SOL = RHS

      IF ( factor == 0 .OR. factor == 1 ) THEN
        itref_max_block = itref_max - 1
        CALL WCP_block_solve( dims, n, m, SOL( dims%x_s : dims%x_e ),          &
                              SOL( dims%c_s : dims%c_e ),                      &
                              SOL( dims%y_s : dims%y_e ),                      &
                              A_ne, A_val, A_col, A_ptr, DIAG_X, ldiag_x,      &
                              SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,           &
                              SOL_y, BEST_y, RES_y, RES_x, itref_max_block,    &
                              K, FACTORS, CNTL, print_level, control, inform )
      ELSE
!       write(6,*) ' sol ', sol
        CALL SILS_solve( K, FACTORS, SOL, CNTL, SINFO )
!       write(6,*) ' sol ', sol
        inform%factorization_status = SINFO%flag
        IF ( SINFO%flag /= 0 ) THEN
          IF ( control%error > 0 .AND. print_level >= 1 )                      &
              WRITE( control%error, 2090 ) SINFO%flag
          inform%status = - 6 ; RETURN
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

        IF ( barrier_free == zero ) THEN
          RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free )
        ELSE
          RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free ) -      &
            barrier_free * SOL( : dims%x_free )
        END IF
        RES( dims%x_free + 1 : dims%x_e ) =                                    &
          RHS( dims%x_free + 1 : dims%x_e ) -                                  &
             BARRIER_X * SOL( dims%x_free + 1 : dims%x_e )
        RES( dims%y_s : dims%y_e ) = RHS( dims%y_s : dims%y_e )
        RES( dims%c_s : dims%c_e ) = RHS( dims%c_s : dims%c_e ) -              &
          BARRIER_C * SOL( dims%c_s : dims%c_e )

!  Include the contribution from A

        CALL WCP_AX( n, RES( dims%x_s : dims%x_e ), m, A_ne, A_val, A_col,     &
                     A_ptr, m, SOL( dims%y_s : dims%y_e ), '-T' )
        CALL WCP_AX( m, RES( dims%y_s : dims%y_e ), m, A_ne, A_val, A_col,     &
                     A_ptr, n, SOL( dims%x_s : dims%x_e ), '- ' )

!  Include the contribution from the slack variables

        RES( dims%c_s : dims%c_e ) =                                           &
          RES( dims%c_s : dims%c_e ) + SCALE_C * SOL( dims%y_i : dims%y_e )
        RES( dims%y_i : dims%y_e ) =                                           &
          RES( dims%y_i : dims%y_e ) + SCALE_C * SOL( dims%c_s : dims%c_e )

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
            BEST = SOL
          ELSE
            SOL = BEST
            EXIT
          END IF
        ELSE
          BEST = SOL
        END IF
   
        IF ( factor == 0 .OR. factor == 1 ) THEN
          CALL WCP_block_solve( dims, n, m, RES( dims%x_s : dims%x_e ),        &
                                RES( dims%c_s : dims%c_e ),                    &
                                RES( dims%y_s : dims%y_e ),                    &
                                A_ne, A_val, A_col, A_ptr,                     &
                                DIAG_X, ldiag_x, SCALE_C, DIAG_C, ldiag_c_l,   &
                                ldiag_c_u, SOL_y,                              &
                                BEST_y, RES_y, RES_x,                          &
                                itref_max_block, K, FACTORS, CNTL,             &
                                print_level, control, inform )

        ELSE
          CALL SILS_solve( K, FACTORS, RES, CNTL, SINFO )
          inform%factorization_status = SINFO%flag
          IF ( SINFO%flag /= 0 ) THEN
            IF ( control%error > 0 .AND. print_level >= 1 )                    &
               WRITE( control%error, 2090 ) SINFO%flag
            inform%status = - 6 ; RETURN
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

        SOL = SOL + RES

!  Print final residuals if required

        IF ( control%out > 0 .AND. print_level >= 3                            &
             .AND. it == itref_max ) THEN
   
!  Remember the barrier terms, and any diagonal perturbations

          IF ( barrier_free == zero ) THEN
            RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free )
          ELSE
            RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free ) -    &
              barrier_free * SOL( : dims%x_free )
          END IF
          RES( dims%x_free + 1 : dims%x_e ) =                                  &
            RHS( dims%x_free + 1 : dims%x_e ) -                                &
               BARRIER_X * SOL( dims%x_free + 1 : dims%x_e )
          RES( dims%y_s : dims%y_e ) = RHS( dims%y_s : dims%y_e )
          RES( dims%c_s : dims%c_e ) = RHS( dims%c_s : dims%c_e ) -            &
            BARRIER_C * SOL( dims%c_s : dims%c_e )

!  Include the contribution from A

          CALL WCP_AX( n, RES( dims%x_s : dims%x_e ), m, A_ne, A_val, A_col,   &
                       A_ptr, m, SOL( dims%y_s : dims%y_e ), '-T' )
          CALL WCP_AX( m, RES( dims%y_s : dims%y_e ), m, A_ne, A_val, A_col,   &
                       A_ptr, n, SOL( dims%x_s : dims%x_e ), '- ' )

!  Include the contribution from the slack variables

          RES( dims%c_s : dims%c_e ) =                                         &
            RES( dims%c_s : dims%c_e ) + SCALE_C * SOL( dims%y_i : dims%y_e )
          RES( dims%y_i : dims%y_e ) =                                         &
            RES( dims%y_i : dims%y_e ) + SCALE_C * SOL( dims%c_s : dims%c_e )

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

        IF ( barrier_free == zero ) THEN
          RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free )
        ELSE
          RES( dims%x_s : dims%x_free ) = RHS( dims%x_s : dims%x_free ) -      &
            barrier_free * SOL( : dims%x_free )
        END IF
        RES( dims%x_free + 1 : dims%x_e ) =                                    &
          RHS( dims%x_free + 1 : dims%x_e ) -                                  &
               BARRIER_X * SOL( dims%x_free + 1 : dims%x_e )
        RES( dims%y_s : dims%y_e ) = RHS( dims%y_s : dims%y_e )
        RES( dims%c_s : dims%c_e ) = RHS( dims%c_s : dims%c_e ) -              &
          BARRIER_C * SOL( dims%c_s : dims%c_e )

!  Include the contribution from A

        CALL WCP_AX( n, RES( dims%x_s : dims%x_e ), m, A_ne, A_val, A_col,     &
                     A_ptr, m, SOL( dims%y_s : dims%y_e ), '-T' )
        CALL WCP_AX( m, RES( dims%y_s : dims%y_e ), m, A_ne, A_val, A_col,     &
                     A_ptr, n, SOL( dims%x_s : dims%x_e ), '- ' )

!  Include the contribution from the slack variables

        RES( dims%c_s : dims%c_e ) =                                           &
          RES( dims%c_s : dims%c_e ) + SCALE_C * SOL( dims%y_i : dims%y_e )
        RES( dims%y_i : dims%y_e ) =                                           &
          RES( dims%y_i : dims%y_e ) + SCALE_C * SOL( dims%c_s : dims%c_e )

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

      RHS = SOL
      inform%status = 0

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
 2090 FORMAT( '   **  Error return ', I0, ' from SILS_solve ' ) 

!  End of WCP_iterative_refinement

      END SUBROUTINE WCP_iterative_refinement

!-*-*-*-*-*-   W C P _ B L O C K _ S O L V E   S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE WCP_block_solve( dims, n, m, RHS_x, RHS_c, RHS_y, A_ne,       &
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

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
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
      TYPE ( WCP_control_type ), INTENT( IN ) :: control        
      TYPE ( WCP_inform_type ), INTENT( INOUT ) :: inform

!  Obtain DIAG_X(inv) * rhs_x ...

      RHS_x = RHS_x / DIAG_X
      IF (  m == 0 ) RETURN

!  ... and DIAG_C(inv) * rhs_c

      RHS_c = RHS_c / DIAG_C

!  Now find A * DIAG_X(inv) * rhs_x - rhs_y - SCALE_C * DIAG_C(inv) * rhs_c

      RHS_y( : dims%c_equality ) = - RHS_y( : dims%c_equality )
      RHS_y( dims%c_l_start : ) = - RHS_y( dims%c_l_start : ) - SCALE_C *RHS_c
      CALL WCP_AX( m, RHS_y, m, A_ne, A_val, A_col, A_ptr, n, RHS_x, '+ ' )

!  Next solve K sol_y = A * DIAG_X(inv) * rhs_x - rhs_y - DIAG_X(inv) * rhs_c, 
!  where K = A * DIAG_X(inv) * A(trans) + SCALE_C * DIAG_C(inv) * SCALE_C

      CALL WCP_block_iterative_refinement(                                     &
               dims, n, m, K, FACTORS, CNTL, RHS_y, A_ne, A_val, A_col,        &
               A_ptr, DIAG_X, ldiag_x, SCALE_C, DIAG_C, ldiag_c_l, ldiag_c_u,  &
               SOL_y, BEST_y, RES_y, RES_x, itref_max, print_level, control,   &
               inform )

!  Finally, recover rhs_x = DIAG_X(inv) ( rhs_x - A(trans) sol_y ) ...

      RHS_x = RHS_x * DIAG_X

      CALL WCP_AX( n, RHS_x, m, A_ne, A_val, A_col, A_ptr, m, RHS_y, '-T' )

      RHS_x = RHS_x / DIAG_X

!  ... and rhs_c = DIAG_C(inv) ( rhs_c + SCALE_C * sol_y )

      RHS_c = RHS_c + ( SCALE_C / DIAG_C ) * RHS_y( dims%c_l_start : )

      RETURN

!  End of subroutine WCP_block_solve

      END SUBROUTINE WCP_block_solve

!-*-*-*-*-*-*-   W C P _ B L O C K _ I-R   S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE WCP_block_iterative_refinement(                               &
                            dims, n, m, K, FACTORS, CNTL, RHS_y,               &
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

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
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
      TYPE ( WCP_control_type ), INTENT( IN ) :: control        
      TYPE ( WCP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: it
      REAL ( KIND = wp ) :: res_norm, old_res
      TYPE ( SILS_sinfo ) :: SINFO

      old_res = MAXVAL( ABS( RHS_y ) )
      IF ( control%out > 0 .AND. print_level >= 3 )                            &
        WRITE( control%out, 2000 ) old_res

      SOL_y = RHS_y

      CALL SILS_solve( K, FACTORS, SOL_y, CNTL, SINFO )
      inform%factorization_status = SINFO%flag
      IF ( SINFO%flag /= 0 ) THEN
        IF ( control%error > 0 .AND. print_level >= 1 )                        &
          WRITE( control%error, 2040 ) SINFO%flag
        inform%status = - 6
        RETURN
      END IF

      IF ( control%out > 0 .AND. print_level >= 3 )                            &
        WRITE( control%out, 2010 ) MAXVAL( ABS( SOL_y ) )

      IF ( itref_max > 0 ) THEN
        DO it = 1, itref_max

!  Form res_x = DIAG_X(inv) * A(trans) sol_y

          RES_x = zero
          CALL WCP_AX( n, RES_x, m, A_ne, A_val, A_col,  A_ptr, m, SOL_y, '+T' )
          RES_x = RES_x / DIAG_X

!  Form res_y = rhs_y - A * res_x - SCALE_C * DIAG_C(inv) * SCALE_C * sol_y

          RES_y( : dims%c_equality ) = RHS_y( : dims%c_equality )
          RES_y( dims%c_l_start : ) = RHS_y( dims%c_l_start : )                &
                 - ( SCALE_C / DIAG_C ) * SCALE_C * SOL_y( dims%c_l_start : )
          CALL WCP_AX( m, RES_y, m, A_ne, A_val, A_col, A_ptr, n, RES_x, '- ' )
  
          IF ( control%out > 0 .AND. print_level >= 3 )                        &
            WRITE( control%out, 2000 ) MAXVAL( ABS( RES_y ) )
  
          res_norm = MAXVAL( ABS( RES_y ) )
          IF ( it > 1 ) THEN
            IF ( res_norm < old_res ) THEN
              old_res = res_norm
              BEST_y = SOL_y
            ELSE
              SOL_y = BEST_y
              EXIT
            END IF
          ELSE
            BEST_y = SOL_y
          END IF
     
          CALL SILS_solve( K, FACTORS, RES_y, CNTL, SINFO )
          inform%factorization_status = SINFO%flag
          IF ( SINFO%flag /= 0 ) THEN
            IF ( control%error > 0 .AND. print_level >= 1 )                    &
              WRITE( control%error, 2040 ) SINFO%flag
            inform%status = - 6
            RETURN
          END IF
          IF ( control%out > 0 .AND. print_level >= 3 )                        &
            WRITE( control%out, 2020 ) MAXVAL( ABS( RES_y ) )
          SOL_y = SOL_y + RES_y

!  Print final residuals if required

          IF ( control%out > 0 .AND. print_level >= 3                          &
              .AND. it == itref_max ) THEN
   
!  Form res_x = DIAG_X(inv) * A(trans) sol_y

            RES_x = zero
            CALL WCP_AX( n, RES_x, m, A_ne, A_val, A_col, A_ptr,               &
                         m, SOL_y, '+T' )
            RES_x = RES_x / DIAG_X

!  Form res_y = rhs_y - A * res_x - SCALE_C * DIAG_C(inv) * SCALE_C * sol_y

            RES_y( : dims%c_equality ) = RHS_y( : dims%c_equality )
            RES_y( dims%c_l_start : ) = RHS_y( dims%c_l_start : )            &
                  - ( SCALE_C / DIAG_C ) * SCALE_C * SOL_y( dims%c_l_start : )
            CALL WCP_AX( m, RES_y, m, A_ne, A_val, A_col, A_ptr, n, RES_x, '- ' )

            WRITE( control%out, 2000 ) MAXVAL( ABS( RES_y ) )
          END IF
        END DO
      ELSE

!  Print residuals if required

        IF ( control%out > 0 .AND. print_level >= 3 ) THEN
   
!  Form res_x = DIAG_X(inv) * A(trans) sol_y

          RES_x = zero
          CALL WCP_AX( n, RES_x, m, A_ne, A_val, A_col, A_ptr, m, SOL_y, '+T' )
          RES_x = RES_x / DIAG_X

!  Form res_y = rhs_y - A * res_x - DIAG_C(inv) sol_y

          RES_y( : dims%c_equality ) = RHS_y( : dims%c_equality )
          RES_y( dims%c_l_start : ) = RHS_y( dims%c_l_start : )              &
                - ( SCALE_C / DIAG_C ) * SCALE_C * SOL_y( dims%c_l_start : )
          CALL WCP_AX( m, RES_y, m, A_ne, A_val, A_col, A_ptr, n, RES_x, '- ' )
      
          WRITE( control%out, 2000 ) MAXVAL( ABS( RES_y ) )
        END IF
      END IF

      IF ( control%out > 0 .AND. print_level >= 3 )                           &
        WRITE( control%out, 2010 ) MAXVAL( ABS( SOL_y ) )

      RHS_y = SOL_y
      inform%status = 0
      RETURN

!  Non-executable statements

 2000 FORMAT( '  res(y) ', ES12.4 )
 2010 FORMAT( '  sol(y) ', ES12.4 )
 2020 FORMAT( ' dsol(y) ', ES12.4 )
 2040 FORMAT( '   **  Error return ', I0, ' from SILS_solve ' ) 

!  End of WCP_block_iterative_refinement

      END SUBROUTINE WCP_block_iterative_refinement

!-*-*-*-*-*-*-   W C P _ R E S I D U A L   S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE WCP_residual( dims, n, m, l_res, A_ne, A_val, A_col, A_ptr,   &
                               DX, DC, DY, RHS_x, RHS_c, RHS_y, RES,           &
                               barrier_free, BARRIER_X, BARRIER_C, SCALE_C,    &
                               errorg, errorc, print_level, control )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Compute the residual of the linear system

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, A_ne, l_res, print_level
      REAL ( KIND = wp ), INTENT( IN ) :: barrier_free
      REAL(  KIND = wp ), INTENT( OUT ) :: errorg, errorc
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
      TYPE ( WCP_control_type ), INTENT( IN ) :: control        

!  Local variables

      INTEGER :: i
      REAL ( KIND = wp ) :: res_tol

      res_tol = epsmch ** 0.5

!  Initalize RES as the zero vector

      RES( dims%y_s : dims%y_e ) = zero

!  Remember the barrier terms, and any diagonal perturbations

      IF ( barrier_free == zero ) THEN
        RES( : dims%x_free ) = zero
      ELSE
        RES( : dims%x_free ) = barrier_free * DX( : dims%x_free )
      END IF
      RES( dims%x_free + 1 : dims%x_e ) = BARRIER_X * DX( dims%x_free + 1 : )
      RES( dims%c_s : dims%c_e ) = BARRIER_C * DC

!  Include the contribution from A and A^T

      CALL WCP_AX( n, RES( dims%x_s : dims%x_e ), m,                           &
                   A_ne, A_val, A_col, A_ptr, m, DY, '+T' )
      CALL WCP_AX( m, RES( dims%y_s : dims%y_e ), m,                           &
                   A_ne, A_val, A_col, A_ptr, n, DX, '+ ' )

!  Include the contribution from the slack variables

      RES( dims%c_s : dims%c_e ) =                                             &
        RES( dims%c_s : dims%c_e ) - SCALE_C * DY( dims%c_l_start : m )
      RES( dims%y_i : dims%y_e ) =                                             &
        RES( dims%y_i : dims%y_e ) - SCALE_C * DC

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

!  End of subroutine WCP_residual

      END SUBROUTINE WCP_residual

!- W C P _ m i n _ p i e c e w i s e _ q u a d r a t i c  S U B R O U T I N E -

      SUBROUTINE WCP_min_piecewise_quadratic( dims, n, m, nbnds, DIST_Z_l,     &
                                              DIST_Z_u, DZ_l, DZ_u, DIST_X_l,  &
                                              DIST_X_u, DX, DIST_Y_l,          &
                                              DIST_Y_u, DY_l,                  &
                                              DY_u, DIST_C_l, DIST_C_u, DC,    &
                                              res_prim_dual, MU_X_l, MU_X_u,   &
                                              MU_C_l, MU_C_u, MU,              &
                                              omega_l, omega_u,                &
                                              alpha, alpha_est, alpha_max,     &
                                              COEF2, COEF1, COEF0, BREAKP,     &
                                              IBREAK, print_level, control )
!                                             inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find a global minimizer of the piecewise quadratic function 
!
!  q(alpha) = 
!    sum_i=1^n,j=1,2 | q_ij alpha^2 + l_ij alpha + c_ij - mu_i | +
!    sum_i=1^m,j=3,4 | q_ij alpha^2 + l_ij alpha + c_ij - mu_i | + 
!    (1 - alpha) res_prim_dual
!
! for all alpha in [0,1] for which additionally
!
!  q_ij alpha^2 + l_ij alpha + c_ij >= omega_l mu_i  (i=1,n,j=1,2 & i=1,m,j=3,4)
!
! and
!
!  q_ij alpha^2 + l_ij alpha + c_ij <= omega_u mu_i  (i=1,n,j=1,2 & i=1,m,j=3,4)
!
!  for omega_l in (0,1), omega_u > 1
!
!  Here:
!    j      q_ij                     lij                          cij
!    1  dx_i dz^L_i   (x_i - x^L_i) dz^L_i + dx_i z^L_i    (x_i - x^L_i) z^L_i 
!    2  dx_i dz^U_i   (x_i - x^U_i) dz^U_i + dx_i z^U_i    (x_i - x^U_i) z^U_i 
!    3  dc_i dy^L_i   (c_i - c^L_i) dy^L_i + dc_i y^L_i    (c_i - c^L_i) y^L_i 
!    4  dc_i dy^U_i   (c_i - c^U_i) dy^U_i + dc_i y^U_i    (c_i - c^U_i) y^U_i 

!  and res_prim_dual = || residual primal/dual optimality ||
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, nbnds, print_level
      REAL ( KIND = wp ), INTENT( IN ) :: res_prim_dual, omega_l, omega_u
      REAL ( KIND = wp ), INTENT( OUT ) :: alpha, alpha_est, alpha_max
      INTEGER, INTENT( OUT ), DIMENSION( 2 * nbnds ) :: IBREAK
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: DX 
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%x_free + 1 : dims%x_l_end ) :: DIST_Z_l, DZ_l,         &
                                                       DIST_X_l, MU_X_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%x_u_start : n ) :: DIST_Z_u, DZ_u, DIST_X_u, MU_X_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%c_l_start : dims%c_l_end ) :: DIST_Y_l, DY_l,          &
                                                      DIST_C_l, MU_C_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%c_u_start : dims%c_u_end ) :: DIST_Y_u, DY_u,          &
                                                      DIST_C_u, MU_C_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%c_l_start : m ) :: DC 
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
        DIMENSION( nbnds ) :: COEF2, COEF1, COEF0, MU
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( 2 * nbnds ) :: BREAKP
      TYPE ( WCP_control_type ), INTENT( IN ) :: control        
!     TYPE ( WCP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: i, l, nroots, inheap, nbreak, cluster
      REAL ( KIND = wp ) :: c0, c1, c2, c0mu, root( 2 )
      REAL ( KIND = wp ) :: tol, close_tol, res_cs, alpha_min, q_alpha_est
      REAL ( KIND = wp ) :: alpha_l, alpha_u, alpha_m, q_l, q_u, q_m, q_min
      LOGICAL :: details, debug

      details = print_level >= 4
      debug = print_level >= 5
      tol = epsmch ** 0.75 
      close_tol = tol

!  Compute the coefficients c2_ij, c1_ij and c0_ij

      l = 0
      DO i = dims%x_free + 1, dims%x_l_end
        l = l + 1
        COEF2( l ) = DX( i ) * DZ_l( i ) 
        COEF1( l ) = DIST_X_l( i ) * DZ_l( i ) + DX( i ) * DIST_Z_l( i ) 
        COEF0( l ) = DIST_X_l( i ) * DIST_Z_l( i ) 
        MU( l ) = MU_X_l( i )
      END DO 
      DO i = dims%x_u_start, n
        l = l + 1
        COEF2( l ) = DX( i ) * DZ_u( i ) 
        COEF1( l ) = - DIST_X_u( i ) * DZ_u( i ) - DX( i ) * DIST_Z_u( i ) 
        COEF0( l ) = DIST_X_u( i ) * DIST_Z_u( i ) 
        MU( l ) = MU_X_u( i )
      END DO 
      DO i = dims%c_l_start, dims%c_l_end
        l = l + 1
        COEF2( l ) = DC( i ) * DY_l( i ) 
        COEF1( l ) = DIST_C_l( i ) * DY_l( i ) + DC( i ) * DIST_Y_l( i ) 
        COEF0( l ) = DIST_C_l( i ) * DIST_Y_l( i ) 
        MU( l ) = MU_C_l( i )
      END DO 
      DO i = dims%c_u_start, dims%c_u_end
        l = l + 1
        COEF2( l ) = DC( i ) * DY_u( i ) 
        COEF1( l ) = - DIST_C_u( i ) * DY_u( i ) - DC( i ) * DIST_Y_u( i ) 
        COEF0( l ) = DIST_C_u( i ) * DIST_Y_u( i ) 
        MU( l ) = MU_C_u( i )
      END DO 

!  Given omega_l in (0,1), find alpha_max, the largest value in [0,1] for which
!
!  c2_ij alpha^2 + c1_ij alpha + c0_ij >= omega_l mu_i 
!     (i=1,n,j=1,2 & i=1,m,j=3,4)
!
!  for all alpha <= alpha_max 

!     WRITE( control%out, "( ' low ' )" )   
      alpha_max = one
      DO l = 1, nbnds
        c2 = COEF2( l ) ; c1 = COEF1( l ) ; c0 = COEF0( l ) - omega_l * MU( l )
        CALL ROOTS_quadratic( c0, c1, c2, tol, nroots, root( 1 ), root( 2 ) ) 
        DO i = 1, nroots
          IF ( root( i ) > tenm7 ) alpha_max = MIN( alpha_max, root( i ) )
        END DO
!       WRITE( control%out, "( ' l, c0, c1, c2 a ', I6, 4ES12.4 )" )           &
!         l, c0, c1, c2, alpha_max
      END DO
!     WRITE( control%out, "( ' alpha_max ', ES12.4 )" ) alpha_max

!  Given omega_u > 1), find alpha_max, the largest value in [0,1] for which
!
!  c2_ij alpha^2 + c1_ij alpha + c0_ij >= omega_u mu_i  
!     (i=1,n,j=1,2 & i=1,m,j=3,4)
!
!  for all alpha <= alpha_max 

!     WRITE( control%out, "( ' high ' )" )   
      DO l = 1, nbnds
        c2 = COEF2( l ) ; c1 = COEF1( l ) ; c0 = COEF0( l ) - omega_u * MU( l )
        CALL ROOTS_quadratic( c0, c1, c2, tol, nroots, root( 1 ), root( 2 ) ) 
        DO i = 1, nroots
          IF ( root( i ) > tenm7 ) alpha_max = MIN( alpha_max, root( i ) )
        END DO
!       WRITE( control%out, "( ' l, c0, c1, c2 a ', I6, 4ES12.4 )" )           &
!         l, c0, c1, c2, alpha_max
      END DO

!  Compute alpha_est, the global minimizer of the overestimator
!
!  o(alpha) = 
!    alpha^2 [ sum_i=1^n,j=1,2 | q_ij | + sum_i=1^m,j=3,4 | q_ij | ] +
!    (1-alpha) [ sum_i=1^n,j=1,2 | c_ij - mu_i | + 
!                sum_i=1^m,j=3,4 | c_ij - mu_i | ] +
!    (1 - alpha) res_prim_dual
!
! for all alpha in [0,alpha_max]

      c0 = res_prim_dual ; c2 = zero
      DO l = 1, nbnds
        c0 = c0 + ABS( COEF0( l ) - MU( l ) )
        c2 = c2 + ABS( COEF2( l ) )
      END DO

      IF ( c2 /= zero ) THEN
        alpha_est = MIN( alpha_max, half * c0 / c2 )
      ELSE
        alpha_est = alpha_max
      END IF

      q_alpha_est = ( one - alpha_est ) * res_prim_dual
      DO l = 1, nbnds
        q_alpha_est = q_alpha_est + ABS( COEF0( l ) - MU( l ) +               &
          alpha_est * ( COEF1( l ) + alpha_est * COEF2( l ) ) )
      END DO

!  Next find all of the roots (the breakpoints) of  
!     c2_ij alpha^2 + c1_ij alpha + c0_ij - mu_i =0
!  (i=1,n,j=1,2 & i=1,m,j=3,4) that lie in [0,alpha_max].

      nbreak = 0
      DO l = 1, nbnds
        c2 = COEF2( l ) ; c1 = COEF1( l ) ; c0 = COEF0( l ) - MU( l )
        CALL ROOTS_quadratic( c0, c1, c2, tol, nroots, root( 1 ), root( 2 ) ) 
        IF ( nroots >= 1 ) THEN 
          IF ( root( 1 ) > zero .AND. root( 1 ) < alpha_max ) THEN
            nbreak = nbreak + 1
            BREAKP( nbreak ) = root( 1 )
            IBREAK( nbreak ) = l
          END IF
        END IF 
        IF ( nroots == 2 ) THEN 
          IF ( root( 2 ) > zero .AND. root( 2 ) < alpha_max ) THEN
            nbreak = nbreak + 1
            BREAKP( nbreak ) = root( 2 )
            IBREAK( nbreak ) = l
          END IF
        END IF 
      END DO

!  Now use heapsort to arrange the breakpoints in increasing order.
      
      CALL SORT_heapsort_build( nbreak, BREAKP, inheap, INDA = IBREAK )
      IF ( details ) WRITE( control%out,                                       &
       "(  I8, ' breakpoints, interval is [0,', ES12.4, ']' )" )               &
          nbreak, alpha_max

!  Build the coefficients of the quadratic q(alpha) at the start of 
!  the initial interval

      c0 = res_prim_dual ; c1 = - res_prim_dual ; c2 = zero
      DO l = 1, nbnds
        c0mu = COEF0( l ) - MU( l )
        IF ( c0mu < zero .OR.                                                  &
             ( c0mu == zero .AND.( COEF1( l ) < zero .OR.                      &
               ( COEF1( l ) == zero .AND. COEF2( l ) < zero ) ) ) ) THEN
          COEF0( l ) = - c0mu
          COEF1( l ) = - COEF1( l ) ; COEF2( l ) = - COEF2( l )
        ELSE
          COEF0( l ) = c0mu
        END IF
        c0 = c0 + COEF0( l ) ; c1 = c1 + COEF1( l ) ; c2 = c2 + COEF2( l )
      END DO
!     WRITE(6,"( ' c0, c1, c2 ', 3ES12.4 )" ) c0, c1, c2

!  At each stage, consider the piecewise quadratic function between the
!  current breakpoint and the next.

      alpha_l = zero
      q_l = c0
      q_min = q_l
      alpha_min = alpha_l      

      cluster = 0
      IF ( details ) WRITE( control%out, 2000 )
      DO
        IF ( nbreak > 0 ) THEN
          alpha_u = BREAKP( 1 )
          l = IBREAK( 1 )

!  Check that the interval upper bound should not really be considered
!  as a cluster of nearly identical points. If so, simply add the point
!  to the cluster, adjust the coefficient of the quadratic function and 
!  find the next break point.

          IF ( alpha_u - alpha_l < close_tol ) THEN
            COEF0( l ) = - COEF0( l ) ; COEF1( l ) = - COEF1( l )
            COEF2( l ) = - COEF2( l )
            c0 = c0 + two * COEF0( l ) ; c1 = c1 + two * COEF1( l )
            c2 = c2 + two * COEF2( l )
            CALL SORT_heapsort_smallest( nbreak, BREAKP, inheap, INDA = IBREAK )
            nbreak = nbreak - 1
            CYCLE
          END IF
        ELSE
          alpha_u = alpha_max
        END IF
        cluster = cluster + 1

!  Compute the minimizer of q(alpha) in [alpha_l,alpha_u], and 
!  see if this is the new candidate incumbent global minimizer

        q_u = c0 + alpha_u * ( c1 + alpha_u * c2 )

        IF ( debug ) THEN
          WRITE(  control%out, "( ' c0, c1, c2 ', 3ES12.4 )" ) c0, c1, c2
          res_cs = ( 1 - alpha_u ) * res_prim_dual
          DO i = dims%x_free + 1, dims%x_u_start - 1
            res_cs = res_cs + ABS( ( DIST_Z_l( i ) + alpha_u * DZ_l( i ) ) *   &
                         ( DIST_X_l( i ) + alpha_u * DX( i ) ) - MU_X_l( i ) )
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            res_cs = res_cs + ABS( ( DIST_Z_l( i ) + alpha_u * DZ_l( i ) ) *   &
                         ( DIST_X_l( i ) + alpha_u * DX( i ) ) - MU_X_l( i ) ) &
                            + ABS( ( DIST_Z_u( i ) - alpha_u * DZ_u( i ) ) *   &
                         ( DIST_X_u( i ) - alpha_u * DX( i ) ) - MU_X_u( i ) )
          END DO
          DO i = dims%x_l_end + 1, n
            res_cs = res_cs + ABS( ( DIST_Z_u( i ) - alpha_u * DZ_u( i ) ) *   &
                         ( DIST_X_u( i ) - alpha_u * DX( i ) ) - MU_X_u( i ) )
          END DO

!  Slack variables:

          DO i = dims%c_l_start, dims%c_u_start - 1
            res_cs = res_cs + ABS( ( DIST_Y_l( i ) + alpha_u * DY_l( i ) ) *   &
                         ( DIST_C_l( i ) + alpha_u * DC( i ) ) - MU_C_l( i ) )
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            res_cs = res_cs + ABS( ( DIST_Y_l( i ) + alpha_u * DY_l( i ) ) *   &
                         ( DIST_C_l( i ) + alpha_u * DC( i ) ) - MU_C_l( i ) ) &
                            + ABS( ( DIST_Y_u( i ) - alpha_u * DY_u( i ) ) *   &
                         ( DIST_C_u( i ) - alpha_u * DC( i ) ) - MU_C_u( i ) )
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            res_cs = res_cs + ABS( ( DIST_Y_u( i ) - alpha_u * DY_u( i ) ) *   &
                         ( DIST_C_u( i ) - alpha_u * DC( i ) ) - MU_C_u( i ) )
          END DO
          WRITE( control%out, "( ' q_upper = ', 56X, ES12.4 )" ) res_cs
        END IF

        IF ( c2 > zero ) THEN
          alpha_m = - half * c1 / c2
          IF ( alpha_m > alpha_l .AND. alpha_m < alpha_u ) THEN
            q_m = c0 + alpha_m * ( c1 + alpha_m * c2 )
            IF ( q_m < q_min ) THEN
              q_min = q_m
              alpha_min = alpha_m
            END IF
!           WRITE( control%out, 2010 )                                         &
            IF ( details ) WRITE( control%out, 2010 )                          &
              cluster, alpha_l, alpha_m, alpha_u, q_l, q_m, q_u 
          ELSE
            IF ( q_u < q_min ) THEN
              q_min = q_u
              alpha_min = alpha_u
            END IF
!           WRITE( control%out, 2020 )                                         &
            IF ( details ) WRITE( control%out, 2020 )                          &
              cluster, alpha_l, alpha_u, q_l, q_u 
          END IF
        ELSE
          IF ( q_u < q_min ) THEN
            q_min = q_u
            alpha_min = alpha_u
          END IF
!         WRITE( control%out, 2020 )                                           &
          IF ( details ) WRITE( control%out, 2020 )                            &
            cluster, alpha_l, alpha_u, q_l, q_u 
        END IF

!  Move to the end of the interval

        IF ( nbreak == 0 ) EXIT
        alpha_l = alpha_u
        q_l = q_u

!  Adjust the coefficients of the quadratic to reflect the values in the 
!  next interval

        COEF0( l ) = - COEF0( l ) ; COEF1( l ) = - COEF1( l )
        COEF2( l ) = - COEF2( l )
        c0 = c0 + two * COEF0( l ) ; c1 = c1 + two * COEF1( l )
        c2 = c2 + two * COEF2( l )
        CALL SORT_heapsort_smallest( nbreak, BREAKP, inheap, INDA = IBREAK )
        nbreak = nbreak - 1
      END DO

      alpha = alpha_min
      IF ( details ) THEN
        WRITE( control%out, "( /, '   alpha_max       alpha       alpha_est ',  &
       &  '      q(alpha)     q(alpha_est) ', /, 3ES14.7, 2ES15.7 )" )          &
          alpha_max, alpha, alpha_est, q_min, q_alpha_est
      END IF

      RETURN

! Non-executable statements

 2000 FORMAT( ' cluster    alpha_l     alpha_m     alpha_u     ',              &
              '    q_l         q_m         q_u')
 2010 FORMAT( I7, 6ES12.4 )
 2020 FORMAT( I7, ES12.4, '      -     ', 2ES12.4, '      -     ', ES12.4  )

!  End of subroutine WCP_min_piecewise_quadratic

      END SUBROUTINE WCP_min_piecewise_quadratic

!-*- W C P _ m i n _ p i e c e w i s e _ q u a r t i c  S U B R O U T I N E -*-

      SUBROUTINE WCP_min_piecewise_quartic( dims, n, m, nbnds, DIST_Z_l,       &
                                            DIST_Z_u, DZ_l, DZ_u, DZ_c_l,      &
                                            DZ_c_u, DIST_X_l, DIST_X_u, DX,    &
                                            DX_c, DIST_Y_l, DIST_Y_u, DY_l,    &
                                            DY_u, DY_c_l, DY_c_u, DIST_C_l,    &
                                            DIST_C_u, DC, DC_c, res_prim_dual, &
                                            MU_X_l, MU_X_u, MU_C_l, MU_C_u,    &
                                            MU, omega_l, omega_u,              &
                                            alpha, alpha_max,                  &
                                            COEF4, COEF3, COEF1, COEF0,        &
                                            BREAKP, IBREAK,                    &
                                            print_level, control )
!                                           inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find a global minimizer of the piecewise quadratic function 
!
!  q(alpha) = 
!  sum_i=1^n,j=1,2 | c4_ij alpha^4 + c3_ij alpha^3 
!                    + c1_ij alpha + c0_ij - mu_i | +
!  sum_i=1^m,j=3,4 | c4_ij alpha^4 + c3_ij alpha^3 
!                   + c1_ij alpha + c0_ij - mu_i | +
!    (1 - alpha) res_prim_dual
!
! for all alpha in [0,1] for which additionally
!
!  c4_ij alpha^4 + c3_ij alpha^3 + c1_ij alpha + c0_ij >= omega_l mu_i
!
!  and
!
!  c4_ij alpha^4 + c3_ij alpha^3 + c1_ij alpha + c0_ij <= omega_u mu_i  
!
!    (i=1,n,j=1,2 & i=1,m,j=3,4)
!
!  for omega_l in (0,1), omega_u > 1
!
!  Here:
!    j      c4__ij                       c3_ij            
!    1  dx_c_i dz_c^L_i     dx_c_i dz^L_i + dx_i dz_c^L_i 
!    2  dx_c_i dz_c^U_i     dx_c_i dz^U_i + dx_i dz_c^U_i 
!    3  dc_c_i dy_c^L_i     dc_c_i dy^L_i + dc_i dy_c^L_i 
!    4  dc_c_i dy_c^U_i     dc_c_i dy^U_i + dc_i dy_c^U_i 

!    j                   c1_ij                       c0_ij
!    1  (x_i - x^L_i) dz^L_i + dx_i z^L_i    (x_i - x^L_i) z^L_i 
!    2  (x_i - x^U_i) dz^U_i + dx_i z^U_i    (x_i - x^U_i) z^U_i 
!    3  (c_i - c^L_i) dy^L_i + dc_i y^L_i    (c_i - c^L_i) y^L_i 
!    4  (c_i - c^U_i) dy^U_i + dc_i y^U_i    (c_i - c^U_i) y^U_i 

!  and res_prim_dual = || residual primal/dual optimality ||
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, m, nbnds, print_level
      REAL ( KIND = wp ), INTENT( IN ) :: res_prim_dual, omega_l, omega_u
      REAL ( KIND = wp ), INTENT( OUT ) :: alpha, alpha_max
      INTEGER, INTENT( OUT ), DIMENSION( 4 * nbnds ) :: IBREAK
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: DX, DX_c 
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%x_free + 1 : dims%x_l_end ) :: DIST_Z_l, DZ_l, DZ_c_l, &
                                                       DIST_X_l, MU_X_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%x_u_start : n ) :: DIST_Z_u, DZ_u, DZ_c_u,             &
                                           DIST_X_u, MU_X_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%c_l_start : dims%c_l_end ) :: DIST_Y_l, DY_l, DY_c_l,  &
                                                      DIST_C_l, MU_C_l
      REAL ( KIND = wp ), INTENT( IN ),                                        &
        DIMENSION( dims%c_u_start : dims%c_u_end ) :: DIST_Y_u, DY_u, DY_c_u,  &
                                                      DIST_C_u, MU_C_u
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( dims%c_l_start : m ) :: DC, DC_c
      REAL ( KIND = wp ), INTENT( OUT ),                                       &
                          DIMENSION( nbnds ) :: COEF4, COEF3, COEF1, COEF0, MU
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( 4 * nbnds ) :: BREAKP
      TYPE ( WCP_control_type ), INTENT( IN ) :: control        
!     TYPE ( WCP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: i, l, nroots, inheap, nbreak, cluster
      REAL ( KIND = wp ) :: c0, c1, c2, c3, c4, c0mu, root( 4 ), tiny4
      REAL ( KIND = wp ) :: tol, close_tol, res_cs, alpha_min
      REAL ( KIND = wp ) :: alpha_l, alpha_u, alpha_m, q_l, q_u, q_m, q_min
      LOGICAL :: debug, details

      details = print_level >= 4
      debug = print_level >= 5
      tol = epsmch ** 0.75 
      close_tol = tol
!     tiny4 = ten ** ( - 20 )
      tiny4 = ten ** ( - 16 )

!  Compute the coefficients c4_ij, c3_ij, c1_ij and c0_ij

      l = 0
      DO i = dims%x_free + 1, dims%x_l_end
        l = l + 1
        COEF4( l ) = DX_c( i ) * DZ_c_l( i ) 
        COEF3( l ) = DX_c( i ) * DZ_l( i ) + DX( i ) * DZ_c_l( i ) 
        IF ( ABS( COEF4( l ) ) < tiny4 * ABS( COEF3( l ) ) ) COEF4( l ) = zero
!       WRITE(6,*) DX_c( i ) * DIST_Z_l( i ) +  DX( i ) * DZ_l( i )            &
!                  + DIST_X_L( i ) * DZ_c_l( i ) 
        COEF1( l ) = DIST_X_l( i ) * DZ_l( i ) + DX( i ) * DIST_Z_l( i ) 
        COEF0( l ) = DIST_X_l( i ) * DIST_Z_l( i ) 
        MU( l ) = MU_X_l( i )
      END DO 
      DO i = dims%x_u_start, n
        l = l + 1
        COEF4( l ) = DX_c( i ) * DZ_c_u( i ) 
        COEF3( l ) = DX_c( i ) * DZ_u( i ) + DX( i ) * DZ_c_u( i ) 
        IF ( ABS( COEF4( l ) ) < tiny4 * ABS( COEF3( l ) ) ) COEF4( l ) = zero
!      WRITE(6,*) - DX_c( i ) * DIST_Z_u( i ) +  DX( i ) * DZ_u( i )           &
!                 - DIST_X_u( i ) * DZ_c_u( i ) 
        COEF1( l ) = - DIST_X_u( i ) * DZ_u( i ) - DX( i ) * DIST_Z_u( i ) 
        COEF0( l ) = DIST_X_u( i ) * DIST_Z_u( i ) 
        MU( l ) = MU_X_u( i )
      END DO 
      DO i = dims%c_l_start, dims%c_l_end
        l = l + 1
        COEF4( l ) = DC_c( i ) * DY_c_l( i ) 
        COEF3( l ) = DC_c( i ) * DY_l( i ) + DC( i ) * DY_c_l( i ) 
        IF ( ABS( COEF4( l ) ) < tiny4 * ABS( COEF3( l ) ) ) COEF4( l ) = zero
!       WRITE(6,*) DC_c( i ) * Y_l( i ) +  DC( i ) * DY_l( i )                 &
!                  + DIST_C_L( i ) * DY_c_l( i ) 
        COEF1( l ) = DIST_C_l( i ) * DY_l( i ) + DC( i ) * DIST_Y_l( i ) 
        COEF0( l ) = DIST_C_l( i ) * DIST_Y_l( i ) 
        MU( l ) = MU_C_l( i )
      END DO 
      DO i = dims%c_u_start, dims%c_u_end
        l = l + 1
        COEF4( l ) = DC_c( i ) * DY_c_u( i ) 
        COEF3( l ) = DC_c( i ) * DY_u( i ) + DC( i ) * DY_c_u( i ) 
        IF ( ABS( COEF4( l ) ) < tiny4 * ABS( COEF3( l ) ) ) COEF4( l ) = zero
!        WRITE(6,*) DC_c( i ) * Y_u( i ) +  DC( i ) * DY_u( i )                &
!                   - DIST_C_u( i ) * DY_c_u( i ) 
        COEF1( l ) = - DIST_C_u( i ) * DY_u( i ) - DC( i ) * DIST_Y_u( i ) 
        COEF0( l ) = DIST_C_u( i ) * DIST_Y_u( i ) 
        MU( l ) = MU_C_u( i )
      END DO 

!  Given omega_l in (0,1), find alpha_max, the largest value in [0,1] for which
!
!  c4_ij alpha^4 + c3_ij alpha^3 + c1_ij alpha + c0_ij 
!    >= omega_l mu_i    (i=1,n,j=1,2 & i=1,m,j=3,4)
!
!  for all alpha <= alpha_max 

      c2 = zero
      alpha_max = one
      DO l = 1, nbnds
        c4 = COEF4( l ) ; c3 = COEF3( l )
        c1 = COEF1( l ) ; c0 = COEF0( l ) - omega_l * MU( l )
        CALL ROOTS_quartic( c0, c1, c2, c3, c4, tol, nroots,                   &
                            root( 1 ), root( 2 ), root( 3 ), root( 4 ) ) 
        DO i = 1, nroots
          IF ( root( i ) > tenm10 ) alpha_max = MIN( alpha_max, root( i ) )
        END DO
      END DO

!  Given omega_u > 1, find alpha_max, the largest value in [0,1] for which
!
!  c4_ij alpha^4 + c3_ij alpha^3 + c1_ij alpha + c0_ij 
!    <= omega_u mu_i    (i=1,n,j=1,2 & i=1,m,j=3,4)
!
!  for all alpha <= alpha_max 

      DO l = 1, nbnds
        c4 = COEF4( l ) ; c3 = COEF3( l )
        c1 = COEF1( l ) ; c0 = COEF0( l ) - omega_u * MU( l )
        CALL ROOTS_quartic( c0, c1, c2, c3, c4, tol, nroots,                   &
                            root( 1 ), root( 2 ), root( 3 ), root( 4 ) ) 
        DO i = 1, nroots
          IF ( root( i ) > tenm10 ) alpha_max = MIN( alpha_max, root( i ) )
        END DO
      END DO

!  Next find all of the roots (the breakpoints) of  
!     c4_ij alpha^4 + c3_ij alpha^3 + c1_ij alpha + c0_ij - mu_i = 0
!  (i=1,n,j=1,2 & i=1,m,j=3,4) that lie in [0,alpha_max].

      nbreak = 0
      DO l = 1, nbnds
        c4 = COEF4( l ) ; c3 = COEF3( l )
        c1 = COEF1( l ) ; c0 = COEF0( l ) - MU( l )
        CALL ROOTS_quartic( c0, c1, c2, c3, c4, tol, nroots,                   &
                            root( 1 ), root( 2 ), root( 3 ), root( 4 ) ) 
        DO i = 1, nroots
          IF ( root( i ) > zero .AND. root( i ) < alpha_max ) THEN
            nbreak = nbreak + 1
            BREAKP( nbreak ) = root( i )
            IBREAK( nbreak ) = l
          END IF
        END DO
      END DO

!  Now use heapsort to arrange the breakpoints in increasing order.
      
      CALL SORT_heapsort_build( nbreak, BREAKP, inheap, INDA = IBREAK )
      IF ( details ) WRITE( control%out,                                       &
       "(  I8, ' breakpoints, interval is [0,', ES12.4, ']' )" )               &
          nbreak, alpha_max

!  Build the coefficients of the quadratic q(alpha) at the start of 
!  the initial interval

      c0 = res_prim_dual ; c1 = - res_prim_dual ; c3 = zero ; c4 = zero
      DO l = 1, nbnds
        c0mu = COEF0( l ) - MU( l )
        IF ( c0mu < zero .OR.                                                  &
             ( c0mu == zero .AND.( COEF1( l ) < zero .OR.                      &
               ( COEF1( l ) == zero .AND. ( COEF3( l ) < zero .OR.             &
                 ( COEF3( l ) == zero .AND. COEF4( l ) < zero ) ) ) ) ) ) THEN
          COEF0( l ) = - c0mu ; COEF1( l ) = - COEF1( l )
          COEF3( l ) = - COEF3( l ) ; COEF4( l ) = - COEF4( l )
        ELSE
          COEF0( l ) = c0mu
        END IF
        c0 = c0 + COEF0( l ) ; c1 = c1 + COEF1( l )
        c3 = c3 + COEF3( l ) ; c4 = c4 + COEF4( l )
      END DO
!     WRITE(6,"( ' c0, c1, c3, c4 ', 4ES12.4 )" ) c0, c1, c3, c4

!  At each stage, consider the piecewise quadratic function between the
!  current breakpoint and the next.

      alpha_l = zero
      q_l = c0
      q_min = q_l
      alpha_min = alpha_l      

      cluster = 0
      IF ( details ) WRITE( control%out, 2000 )
      DO
        IF ( nbreak > 0 ) THEN
          alpha_u = BREAKP( 1 )
          l = IBREAK( 1 )

!  Check that the interval upper bound should not really be considered
!  as a cluster of nearly identical points. If so, simply add the point
!  to the cluster, adjust the coefficient of the quadratic function and 
!  find the next break point.

          IF ( alpha_u - alpha_l < close_tol ) THEN
            COEF0( l ) = - COEF0( l ) ; COEF1( l ) = - COEF1( l )
            COEF3( l ) = - COEF3( l ) ; COEF4( l ) = - COEF4( l ) ; 
            c0 = c0 + two * COEF0( l ) ; c1 = c1 + two * COEF1( l )
            c3 = c3 + two * COEF3( l ) ; c4 = c4 + two * COEF4( l )
            CALL SORT_heapsort_smallest( nbreak, BREAKP, inheap, INDA = IBREAK )
            nbreak = nbreak - 1
            CYCLE
          END IF
        ELSE
          alpha_u = alpha_max
        END IF
        cluster = cluster + 1

!  Compute the minimizer of q(alpha) in [alpha_l,alpha_u], and 
!  see if this is the new candidate incumbent global minimizer

!       q_u = c0 + alpha_l * ( c1 + alpha_l ** 2 * ( c3 + alpha_l * c4 ) )
!       write(6,"( ' q at start ', ES12.4 )" ) q_u
        q_u = c0 + alpha_u * ( c1 + alpha_u ** 2 * ( c3 + alpha_u * c4 ) )
!       write(6,"( ' res_cs ', ES12.4 )" ) q_u - res_prim_dual * ( 1 - alpha_u )

        IF ( q_u < q_min ) THEN
          q_min = q_u
          alpha_min = alpha_u
        END IF

        IF ( debug ) THEN
          WRITE(  control%out, "( ' c0, c1, c3, c4 ', 4ES12.4 )" ) c0, c1, c3, c4
          res_cs = ( 1 - alpha_u ) * res_prim_dual
          DO i = dims%x_free + 1, dims%x_u_start - 1
            res_cs = res_cs + ABS(                                             &
          ( DIST_Z_l( i ) + alpha_u * ( DZ_l( i ) + alpha_u * DZ_c_l( i ) ) ) *&
          ( DIST_X_l( i ) + alpha_u * ( DX( i ) + alpha_u * DX_c( i ) ) )      &
                - MU_X_l( i ) )
          END DO
          DO i = dims%x_u_start, dims%x_l_end
            res_cs = res_cs + ABS(                                             &
          ( DIST_Z_l( i ) + alpha_u * ( DZ_l( i ) + alpha_u * DZ_c_l( i ) ) ) *&
          ( DIST_X_l( i ) + alpha_u * ( DX( i ) + alpha_u * DX_c( i ) ) )      &
                - MU_X_l( i ) ) + ABS(                                         &
          ( DIST_Z_u( i ) - alpha_u * ( DZ_u( i ) + alpha_u * DZ_c_u( i ) ) ) *&
          ( DIST_X_u( i ) - alpha_u * ( DX( i ) + alpha_u * DX_c( i ) ) )      &
               - MU_X_u( i ) )
          END DO
          DO i = dims%x_l_end + 1, n
            res_cs = res_cs + ABS(                                             &
          ( DIST_Z_u( i ) - alpha_u * ( DZ_u( i ) + alpha_u * DZ_c_u( i ) ) ) *&
          ( DIST_X_u( i ) - alpha_u * ( DX( i ) + alpha_u * DX_c( i ) ) )      &
               - MU_X_u( i ) )
          END DO

!  Slack variables:

          DO i = dims%c_l_start, dims%c_u_start - 1
            res_cs = res_cs + ABS(                                             &
          ( DIST_Y_l( i ) + alpha_u * ( DY_l( i ) + alpha_u * DY_c_l( i ) ) ) *&
          ( DIST_C_l( i ) + alpha_u * ( DC( i ) + alpha_u * DC_c( i ) ) )      &
                - MU_C_l( i) )
          END DO
          DO i = dims%c_u_start, dims%c_l_end
            res_cs = res_cs + ABS(                                             &
          ( DIST_Y_l( i ) + alpha_u * ( DY_l( i ) + alpha_u * DY_c_l( i ) ) ) *&
          ( DIST_C_l( i ) + alpha_u * ( DC( i ) + alpha_u * DC_c( i ) ) )      &
                - MU_C_l( i ) ) + ABS(                                         &
          ( DIST_Y_u( i ) - alpha_u * ( DY_u( i ) + alpha_u * DY_c_u( i ) ) ) *&
          ( DIST_C_u( i ) - alpha_u * ( DC( i ) + alpha_u * DC_c( i ) ) )      &
               - MU_C_u( i ) )
          END DO
          DO i = dims%c_l_end + 1, dims%c_u_end
            res_cs = res_cs + ABS(                                             &
          ( DIST_Y_u( i ) - alpha_u * ( DY_u( i ) + alpha_u * DY_c_u( i ) ) ) *&
          ( DIST_C_u( i ) - alpha_u * ( DC( i ) + alpha_u * DC_c( i ) ) )      &
               - MU_C_u( i ) )
          END DO
          WRITE( control%out, "( ' q_upper = ', 56X, ES12.4 )" ) res_cs
        END IF

!  Find the stationary points of q(alpha) (if any)

        CALL ROOTS_cubic( c1, two * c2, three * c3, four * c4, tol, nroots,    &
                          root( 1 ), root( 2 ), root( 3 ) ) 

!  If the stationary points lie within the interval, check to see whether
!  any of them gives a lower function value

        DO i = 1, nroots
          alpha_m = root( i )
          IF ( alpha_m > alpha_l .AND. alpha_m < alpha_u ) THEN
            q_m = c0 + alpha_m * ( c1 + alpha_m ** 2 * ( c3 + alpha_m * c4 ) )
            IF ( q_m < q_min ) THEN
              q_min = q_m
              alpha_min = alpha_m
            END IF
          END IF
        END DO

        IF ( alpha_min > alpha_l .AND. alpha_min < alpha_u ) THEN
          IF ( details ) WRITE( control%out, 2010 )                            &
            cluster, alpha_l, alpha_min, alpha_u, q_l, q_min, q_u 
        ELSE
          IF ( details ) WRITE( control%out, 2020   )                          &
            cluster, alpha_l, alpha_u, q_l, q_u 
        END IF

!  Move to the end of the interval

        IF ( nbreak == 0 ) EXIT
        alpha_l = alpha_u
        q_l = q_u

!  Adjust the coefficients of the quadratic to reflect the values in the 
!  next interval

        COEF0( l ) = - COEF0( l ) ; COEF1( l ) = - COEF1( l )
        COEF3( l ) = - COEF3( l ) ; COEF4( l ) = - COEF4( l )
        c0 = c0 + two * COEF0( l ) ; c1 = c1 + two * COEF1( l )
        c3 = c3 + two * COEF3( l ) ; c4 = c4 + two * COEF4( l )
        CALL SORT_heapsort_smallest( nbreak, BREAKP, inheap, INDA = IBREAK )
        nbreak = nbreak - 1

      END DO
      alpha = alpha_min
      IF ( details )                                                           &
        WRITE( control%out, "( ' final step, value ', 2ES12.4 )" ) alpha, q_min

      RETURN

! Non-executable statements

 2000 FORMAT( ' cluster    alpha_l     alpha_m     alpha_u     ',              &
              '    q_l         q_m         q_u')
 2010 FORMAT( I7, 6ES12.4 )
 2020 FORMAT( I7, ES12.4, '      -     ', 2ES12.4, '      -     ', ES12.4  )

!  End of subroutine WCP_min_piecewise_quartic

      END SUBROUTINE WCP_min_piecewise_quartic

!-*-*-*-*-*-*-*-   W C P _ F O R M _ S _ C   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE WCP_form_Schur_complement(                                    &
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

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
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

!     COL_count( : n ) = Abycol_ptr( 2 : n + 1 ) - Abycol_ptr( : n )
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
!         ELSE
!           WRITE( 6, "( ' j = ', i6, ' i = ', /, ( 10i6 ) )" ) &
!             j, S_row( S_colptr( j ) : S_colptr( j + 1 ) - 1 )
!           WRITE( 6, "( ' This should be impossible: '/,                      &
!          &             ' please report to n.gould@rl.ac.uk' )" )
!           STOP
          END IF
        END DO
      END IF

  900 CONTINUE
      RETURN

!  Non-executable statements

 2010 FORMAT( ' ** Error return from WCP_form_Schur_complement', /,            &
              '    Value of n (number of rows of A) is set to ', I0, /,       &
              '    but n must be at least 1 ' )
 2020 FORMAT( ' ** Error return from WCP_form_Schur_complement', /,            &
              '    Value of m (number of columns of A) is set to ', I0, /,    &
              '    but m must be at least 1 ' )
 2030 FORMAT( ' ** Error return from WCP_form_Schur_complement', /,            &
              '    Value of ls = ', I0, ' too small. ', /,                    &
              '    ls must be at least', I10 )
 2040 FORMAT( ' ** Error return from WCP_form_Schur_complement', /,            &
              '    Arrays S_row and S are too small ', /, '    ls is ', I0, /,&
              '    but must be at least ', I0 )

!  End of WCP_form_Schur_complement

      END SUBROUTINE WCP_form_Schur_complement

!-*-*-*-*-*-*-   W C P _ A _ B Y _ C O L S   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE WCP_A_by_cols( n, m, A_ne, A_val, A_col, A_ptr, B_val,       &
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

!  End of WCP_A_by_cols

      END SUBROUTINE WCP_A_by_cols

!-*-*-*-*-*-*-*-   W C P _ A N A L Y S E  S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE WCP_analyse( dims, n, m, A_ne, A_val, A_col, A_ptr, SCALE_C,  &
                              factor, nnzks, lk, liw, ldiag_x, ldiag_c_l,      &
                              ldiag_c_u, IW, Abycol_val, Abycol_row,           &
                              Abycol_ptr, K_colptr, DIAG_X, DIAG_C,            &
                              K, FACTORS, CNTL, print_level, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Analyse the sparsity pattern of the block matrix

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( WCP_dims_type ), INTENT( IN ) :: dims
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
      TYPE ( WCP_control_type ), INTENT( IN ) :: control        
      TYPE ( WCP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: i, ierr, l, max_col, max_len, y_ii
      LOGICAL :: printi, printt, printe, get_factors
      CHARACTER ( LEN = 80 ) :: array_name
      TYPE ( SILS_ainfo ) :: AINFO

      inform%status = 0

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

        array_name = 'wcp: data%IW'
        CALL SPACE_resize_array( liw, IW, inform%status,                       &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        array_name = 'wcp: data%Abycol_val'
        CALL SPACE_resize_array( A_ne, Abycol_val, inform%status,              &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        array_name = 'wcp: data%Abycol_row'
        CALL SPACE_resize_array( A_ne, Abycol_row, inform%status,              &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        array_name = 'wcp: data%Abycol_ptr'
        CALL SPACE_resize_array( n + 1, Abycol_ptr, inform%status,             &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

!  Reorder A so that its entries are ordered by columns

        CALL WCP_A_by_cols( n, m, A_ne, A_val, A_col, A_ptr,                   &
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

        array_name = 'wcp: data%DIAG_X'
        CALL SPACE_resize_array( ldiag_x, DIAG_X, inform%status,               &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        ldiag_c_l = dims%c_l_start
        ldiag_c_u = dims%c_u_end

        array_name = 'wcp: data%DIAG_C'
        CALL SPACE_resize_array( ldiag_c_l, ldiag_c_u,                         &
               DIAG_C, inform%status,                                          &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        array_name = 'wcp: data%K_colptr'
        CALL SPACE_resize_array( m + 1, K_colptr, inform%status,               &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        lk = max( A_ne + m, 2 * A_ne, control%valmin )

        array_name = 'wcp: data%K%row'
        CALL SPACE_resize_array( lk, K%row, inform%status,                     &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        array_name = 'wcp: data%K%col'
        CALL SPACE_resize_array( lk, K%col, inform%status,                     &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        array_name = 'wcp: data%K%val'
        CALL SPACE_resize_array( lk, K%val, inform%status,                     &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

!  Permute A so that its entries appear by columns, with the entries in
!  each column sorted by increasing row number. Also form the data
!  structure required to hold K = A DIAG_X(inv) A(transpose) + DIAG_C(inv)

   20   CONTINUE

        IF ( m > 0 ) THEN
          CALL WCP_form_Schur_complement(                                      &
                  dims, n, m, A_ne, Abycol_val, Abycol_row, Abycol_ptr,        &
                  DIAG_X, SCALE_C, DIAG_C, lk, K%val, K%row, K_colptr,         &
                  K%ne, ierr, IW( : n ), IW( n + 1 : ), .FALSE.,               &
                  control%error, print_level )
!         inform%nfacts = inform%nfacts + 1 

!  Check for error returns

          IF ( ierr == 3 .OR. ierr == 4 ) THEN

!  Insufficient space. Allocate more and retry

            IF ( ierr == 3 ) THEN
              lk = 2 * lk
            ELSE
              lk = K%ne
            END IF

            array_name = 'wcp: data%K%row'
            CALL SPACE_resize_array( lk, K%row, inform%status,                 &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= 0 ) RETURN

            array_name = 'wcp: data%K%col'
            CALL SPACE_resize_array( lk, K%col, inform%status,                 &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= 0 ) RETURN

            array_name = 'wcp: data%K%val'
            CALL SPACE_resize_array( lk, K%val, inform%status,                 &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= 0 ) RETURN
  
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

        array_name = 'wcp: data%DIAG_X'
        CALL SPACE_resize_array( ldiag_x, DIAG_X, inform%status,               &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        ldiag_c_l = 0
        ldiag_c_u = 0

        array_name = 'wcp: data%DIAG_C'
        CALL SPACE_resize_array( ldiag_c_l, ldiag_c_u,                         &
               DIAG_C, inform%status,                                          &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        lk = A_ne + n + 2 * dims%nc

        array_name = 'wcp: data%K%row'
        CALL SPACE_resize_array( lk, K%row, inform%status,                     &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        array_name = 'wcp: data%K%col'
        CALL SPACE_resize_array( lk, K%col, inform%status,                     &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

        array_name = 'wcp: data%K%val'
        CALL SPACE_resize_array( lk, K%val, inform%status,                     &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) RETURN

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
          inform%status = - 6 ; RETURN
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

!  Non-executable statements

 2000 FORMAT( '  maximum, average column lengths ', I7, 0P, F10.1, /,          &
              '  number of columns longer than   ', I7, ' is', I7 )
 2020 FORMAT( /, ' ** There are free variables - abandon the Schur-complement',&
              /, '    factorization in favour of one of the augmented matrix' )
 2030 FORMAT( ' ** The maximum column length in A is larger than', /,          &
              '    max_col =', I7, ' - abandon the Schur-complement', /,       &
              '    factorization in favour of one of the augmented matrix', / )
 2040 FORMAT( '   **  Error return ', I0, ' from ', A ) 
 2050 FORMAT( '   **  Warning ', I0, ' from ', A ) 
!2900 FORMAT( ' ** Message from -WCP_analyse-', /,                             &
!             ' Allocation error, for ', A20, /, ' status = ', I0 )

      END SUBROUTINE WCP_analyse

!  End of module WCP

   END MODULE GALAHAD_WCP_double



