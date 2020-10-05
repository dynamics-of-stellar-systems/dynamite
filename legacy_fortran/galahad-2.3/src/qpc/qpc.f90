! THIS VERSION: GALAHAD 2.2 - 21/04/2008 AT 16:00 GMT.

!-*-*-*-*-*-*-*-*-*- G A L A H A D _ Q P C   M O D U L E -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released pre GALAHAD Version 1.0. October 17th 1997
!   update released with GALAHAD Version 2.0. August 11th 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_QPC_double

!               ---------------------------------------------
!               |                                           |
!               | Solve the quadratic program               |
!               |                                           |
!               |    minimize     1/2 x(T) H x + g(T) x + f |
!               |    subject to     c_l <= A x <= c_u       |
!               |                   x_l <=  x  <= x_u       |
!               |                                           |
!               | using an interior-point trust-region      |
!               | approach to find an approximate solution, |
!               ! then crossing over to a working-set       |
!               | method to obtain an accurate solution     |
!               |                                           |
!               ---------------------------------------------

!NOT95USE GALAHAD_CPU_time
      USE GALAHAD_SYMBOLS
      USE GALAHAD_SPACE_double
      USE GALAHAD_SORT_double
      USE GALAHAD_SILS_double, ONLY: SILS_initialize
      USE GALAHAD_QPT_double
      USE GALAHAD_QPP_double
      USE GALAHAD_QPD_double, QPC_data_type => QPD_data_type,                  &
                              QPC_HX => QPD_HX, QPC_AX => QPD_AX
      USE GALAHAD_LSQP_double
      USE GALAHAD_QPA_double
      USE GALAHAD_QPB_double
      USE GALAHAD_SPECFILE_double
      USE GALAHAD_FDC_double

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: QPC_initialize, QPC_read_specfile, QPC_solve, QPC_terminate,   &
                QPC_data_type, QPT_problem_type, SMT_type, SMT_put, SMT_get

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

      TYPE, PUBLIC :: QPC_time_type
        REAL :: total, preprocess, find_dependent, analyse, factorize, solve
      END TYPE

      TYPE, PUBLIC :: QPC_control_type
        INTEGER :: error, out, print_level, indmin, valmin, restore_problem
        REAL ( KIND = wp ) :: infinity, identical_bounds_tol, rho_g, rho_b
        REAL ( KIND = wp ) :: pivot_tol_for_dependencies, zero_pivot
        REAL ( KIND = wp ) :: cpu_time_limit, on_bound_tol
        LOGICAL :: treat_zero_bounds_as_general, array_syntax_worse_than_do_loop
        LOGICAL :: space_critical, deallocate_error_fatal, no_qpa, no_qpb
        LOGICAL :: qpb_or_qpa
        CHARACTER ( LEN = 30 ) :: prefix
        TYPE ( QPA_control_type ) :: QPA_control
        TYPE ( QPB_control_type ) :: QPB_control
      END TYPE

      TYPE, PUBLIC :: QPC_inform_type
        INTEGER :: status, alloc_status, factorization_status
        INTEGER :: factorization_integer, factorization_real, nfacts, nmods
        LOGICAL :: p_found
        REAL ( KIND = wp ) :: obj, non_negligible_pivot
        CHARACTER ( LEN = 80 ) :: bad_alloc
        TYPE ( QPC_time_type ) :: time
        TYPE ( QPA_inform_type ) :: QPA_inform
        TYPE ( QPB_inform_type ) :: QPB_inform
      END TYPE

!----------------------
!   P a r a m e t e r s
!----------------------

      INTEGER, PARAMETER :: max_sc = 200
      INTEGER, PARAMETER :: max_real_store_ratio = 100
      INTEGER, PARAMETER :: max_integer_store_ratio = 100
      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: point01 = 0.01_wp
      REAL ( KIND = wp ), PARAMETER :: point1 = 0.1_wp
      REAL ( KIND = wp ), PARAMETER :: point9 = 0.9_wp
      REAL ( KIND = wp ), PARAMETER :: point99 = 0.99_wp
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: two = 2.0_wp
      REAL ( KIND = wp ), PARAMETER :: four = 4.0_wp
      REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp
      REAL ( KIND = wp ), PARAMETER :: hundred = 100.0_wp
      REAL ( KIND = wp ), PARAMETER :: thousand = 1000.0_wp
      REAL ( KIND = wp ), PARAMETER :: tenm2 = ten ** ( - 2 )
      REAL ( KIND = wp ), PARAMETER :: tenm4 = ten ** ( - 4 )
      REAL ( KIND = wp ), PARAMETER :: k_diag = one
      REAL ( KIND = wp ), PARAMETER :: infinity = HUGE( one )
      REAL ( KIND = wp ), PARAMETER :: epsmch = EPSILON( one )
      REAL ( KIND = wp ), PARAMETER :: teneps = ten * epsmch
      REAL ( KIND = wp ), PARAMETER :: res_large = one
      REAL ( KIND = wp ), PARAMETER :: remote = ten ** 10
      REAL ( KIND = wp ), PARAMETER :: bar_min = zero
      REAL ( KIND = wp ), PARAMETER :: z_min = ten ** ( - 12 )
      LOGICAL :: primal_dual_indicators = .FALSE.

!-------------------------------
!   I n t e r f a c e  B l o c k
!-------------------------------

      INTERFACE TWO_NORM

        FUNCTION SNRM2( n, X, incx )
        REAL :: SNRM2
        INTEGER, INTENT( IN ) :: n, incx
        REAL, INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
!       REAL, INTENT( IN ), DIMENSION( : ) :: X
        END FUNCTION SNRM2

        FUNCTION DNRM2( n, X, incx )
        DOUBLE PRECISION :: DNRM2
        INTEGER, INTENT( IN ) :: n, incx
        DOUBLE PRECISION, INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
!       DOUBLE PRECISION, INTENT( IN ), DIMENSION( : ) :: X
        END FUNCTION DNRM2
        
      END INTERFACE 

   CONTAINS

!-*-*-*-*-*-   Q P C _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE QPC_initialize( data, control )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Default control data for QPC. This routine should be called before
!  QPC_solve
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
!   indmin. An initial guess as to the integer workspace required by SILS
!
!   valmin. An initial guess as to the real workspace required by SILS
! 
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
!   identical_bounds_tol. Any pair of constraint bounds (c_l,c_u) or (x_l,x_u)
!    that are closer than identical_bounds_tol will be reset to the average
!    of their values
!
!   rho_g. The initial value of the penalty parameter used by QPA to penalize
!    general constraints. A non-positive value will be reset to 2 * infinity 
!    norm of the Lagrange multipliers found by QPB or, if QPB has not been 
!    used, 2 * m
!
!   rho_b. The initial value of the penalty parameter used by QPA to penalize
!    simple bound constraints. A non-positive value will be reset to 2 * 
!    infinity norm of the dual variables found by QPB or, if QPB has not been 
!    used, 2 * n
!
!   pivot_tol_for_dependencies. The threshold pivot used by the matrix 
!    factorization when attempting to detect linearly dependent constraints.
!    See the documentation for SILS for details
!
!   zero_pivot. Any pivots smaller than zero_pivot in absolute value will 
!    be regarded to be zero when attempting to detect linearly dependent 
!    constraints
!
!   cpu_time_limit. The maximum CPU time allowed
!
!  LOGICAL control parameters:
!
!   treat_zero_bounds_as_general. If true, any problem bound with the value
!    zero will be treated as if it were a general value
!
!   array_syntax_worse_than_do_loop. If array_syntax_worse_than_do_loop is
!    true, f77-style do loops will be used rather than 
!    f90-style array syntax for vector operations
!
!   space_critical. If true, every effort will be made to use as little
!     space as possible. This may result in longer computation times
!
!   deallocate_error_fatal. If true, any array/pointer deallocation error
!     will terminate execution. Otherwise, computation will continue
!
!   no_qpa. If true, the crossover phase to get an exact solution using
!    the working set solver QPA will be skipped
!
!   no_qpb. If true, the interior-point phase using QPB will be skipped, and 
!    the solution will be found using the working set solver QPA
!
!   qpb_or_qpa. If true the user wishes to use the interior-point phase
!    of the computation, but only to follow this with the working-set phase
!    if the former is unsuccessful. Otherwise both phases are to be used
!    (subject to the requests made in no_qpa and no_qpb).
!
!  CHARACTER control parameters:
!
!  prefix (len=30). All output lines will be prefixed by 
!    %prefix(2:LEN(TRIM(%prefix))-1)
!   where %prefix contains the required string enclosed in 
!   quotes, e.g. "string" or 'string'
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      TYPE ( QPC_data_type ), INTENT( OUT ) :: data
      TYPE ( QPC_control_type ), INTENT( OUT ) :: control        

!  Set control parameters

      CALL QPA_initialize( data, control%QPA_control )
      control%QPA_control%prefix = '" - QPA:"                     '
      CALL QPB_initialize( data, control%QPB_control )
      control%QPB_control%prefix = '" - QPB:"                     '

!  Initalize SILS components (need to do this as data for QPA/B_initiliaze 
!  is INTENT( out ) and thus call to QPB_initialize wipes out QPA-specific data

      CALL SILS_initialize( FACTORS = data%FACTORS, control = data%CNTLA )
      data%CNTLA%ordering = 3
!57V2 data%CNTLA%ordering = 2
!57V3 data%CNTLA%ordering = 5
!57V2 data%CNTLA%scaling = 0
!57V2 data%CNTLA%static_tolerance = zero
!57V2 data%CNTLA%static_level = zero
!     data%CNTLA%fratio = 1.0

!  Integer parameters

      control%error  = 6
      control%out  = 6
      control%print_level = 0
      control%indmin = 1000
      control%valmin = 1000
      control%restore_problem = 2

!  Real parameter

      control%infinity = ten ** 19
      control%identical_bounds_tol = epsmch
      control%rho_g = - one
      control%rho_b = - one
      control%pivot_tol_for_dependencies = half
      control%zero_pivot = epsmch ** 0.75
      control%cpu_time_limit = - one
      control%on_bound_tol = epsmch

!  Logical parameters

      control%treat_zero_bounds_as_general = .FALSE.
      control%array_syntax_worse_than_do_loop = .FALSE.
      control%space_critical = .FALSE.
      control%deallocate_error_fatal  = .FALSE.
      control%no_qpa = .FALSE.
      control%no_qpb = .FALSE.
      control%qpb_or_qpa = .FALSE.

!  Character parameters

      control%prefix = '""                            '

      RETURN

!  End of QPC_initialize

      END SUBROUTINE QPC_initialize

!-*-*-*-*-   Q P C _ R E A D _ S P E C F I L E  S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE QPC_read_specfile( control, device, alt_specname )

!  Reads the content of a specification file, and performs the assignment of 
!  values associated with given keywords to the corresponding control parameters

!  The defauly values as given by QPC_initialize could (roughly) 
!  have been set as:

! BEGIN QPC SPECIFICATIONS (DEFAULT)
!  error-printout-device                             6
!  printout-device                                   6
!  print-level                                       0
!  initial-integer-workspace                         1000
!  initial-real-workspace                            1000
!  restore-problem-on-output                         2
!  infinity-value                                    1.0D+19
!  identical-bounds-tolerance                        1.0D-15
!  initial-rho-g                                     -1.0
!  initial-rho-b                                     -1.0
!  pivot-tolerance-used-for-dependencies             0.5
!  zero-pivot-tolerance                              1.0D-12
!  on-bound-tolerance                                1.0D-15
!  maximum-cpu-time-limit                            -1.0
!  treat-zero-bounds-as-general                      F
!  array-syntax-worse-than-do-loop                   F
!  space-critical                                    F
!  deallocate-error-fatal                            F
!  no-qpa-phase                                      F
!  no-qpb-phase                                      F
!  qpb-or-qpa                                        F
! END QPC SPECIFICATIONS

!  Dummy arguments

      TYPE ( QPC_control_type ), INTENT( INOUT ) :: control        
      INTEGER, INTENT( IN ) :: device
      CHARACTER( LEN = 16 ), OPTIONAL :: alt_specname

!  Programming: Nick Gould and Ph. Toint, January 2002.

!  Local variables

      INTEGER, PARAMETER :: lspec = 21
      CHARACTER( LEN = 16 ), PARAMETER :: specname = 'QPC             '
      TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec

!  Read the specfiles for QPA and QPB

      CALL QPA_read_specfile( control%QPA_control, device )
      CALL QPB_read_specfile( control%QPB_control, device )

!  Define the keywords

     spec%keyword = ''

!  Integer key-words

      spec(  1 )%keyword = 'error-printout-device'
      spec(  2 )%keyword = 'printout-device'
      spec(  3 )%keyword = 'print-level' 
      spec(  4 )%keyword = 'initial-integer-workspace'
      spec(  5 )%keyword = 'initial-real-workspace'
      spec(  6 )%keyword = 'restore-problem-on-output'

!  Real key-words

      spec(  7 )%keyword = 'infinity-value'
      spec(  8 )%keyword = 'identical-bounds-tolerance'
      spec(  9 )%keyword = 'initial-rho-g'
      spec( 10 )%keyword = 'initial-rho-b'
      spec( 11 )%keyword = 'pivot-tolerance-used-for-dependencies'
      spec( 12 )%keyword = 'zero-pivot-tolerance'
      spec( 13 )%keyword = 'maximum-cpu-time-limit'
      spec( 20 )%keyword = 'on-bound-tolerance'

!  Logical key-words

      spec( 14 )%keyword = 'treat-zero-bounds-as-general'
      spec( 15 )%keyword = 'array-syntax-worse-than-do-loop'
      spec( 16 )%keyword = 'space-critical'
      spec( 17 )%keyword = 'deallocate-error-fatal'
      spec( 18 )%keyword = 'no-qpa-phase'
      spec( 19 )%keyword = 'no-qpb-phase'
      spec( 21 )%keyword = 'qpb-or-qpa'

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
      CALL SPECFILE_assign_value( spec( 4 ), control%indmin,                   &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 5 ), control%valmin,                   &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 6 ), control%restore_problem,          &
                                  control%error )

!  Set real value

      CALL SPECFILE_assign_value( spec( 7 ), control%infinity,                 &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 8 ), control%identical_bounds_tol,     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 9 ), control%rho_g,                    &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 10 ), control%rho_b,                   &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 11 ),                                  &
                                  control%pivot_tol_for_dependencies,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 12 ), control%zero_pivot,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 13 ), control%cpu_time_limit,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 20 ), control%on_bound_tol,            &
                                  control%error )

!  Set logical values

      CALL SPECFILE_assign_value( spec( 14 ),                                  &
                                  control%treat_zero_bounds_as_general,        &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 15 ),                                  &
                                  control%array_syntax_worse_than_do_loop,     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 16 ), control%space_critical,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 17 ),                                  &
                                  control%deallocate_error_fatal,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 18 ), control%no_qpa,                  &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 19 ), control%no_qpb,                  &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 21 ), control%qpb_or_qpa,              &
                                  control%error )

      RETURN

      END SUBROUTINE QPC_read_specfile

!-*-*-*-*-*-*-*-*-   Q P C _ S O L V E  S U B R O U T I N E   -*-*-*-*-*-*-*-*-

      SUBROUTINE QPC_solve( prob, C_stat, B_stat, data, control, inform,      &
                            G_p, X_p, Y_p, Z_p )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Solve the quadratic program
!
!     minimize     q(x) = 1/2 x(T) H x + g(T) x + f
!
!     subject to    (c_l)_i <= (Ax)_i <= (c_u)_i , i = 1, .... , m,
!
!        and        (x_l)_i <=   x_i  <= (x_u)_i , i = 1, .... , n,
!
!  where x is a vector of n components ( x_1, .... , x_n ), const is a
!  constant, g is an n-vector, H is a symmetric matrix, 
!  A is an m by n matrix, and any of the bounds (c_l)_i, (c_u)_i
!  (x_l)_i, (x_u)_i may be infinite, using a primal-dual method.
!  The subroutine is particularly appropriate when A and H are sparse
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
!    to be solved since the last call to QPC_initialize, and .FALSE. if
!    a previous call to a problem with the same "structure" (but different
!    numerical data) was made.
!
!   %n is an INTEGER variable, which must be set by the user to the
!    number of optimization parameters, n.  RESTRICTION: %n >= 1
!                 
!   %m is an INTEGER variable, which must be set by the user to the
!    number of general linear constraints, m. RESTRICTION: %m >= 0
!                 
!   %H is a structure of type SMT_type used to hold the LOWER TRIANGULAR part 
!    of H. Four storage formats are permitted:
!
!    i) sparse, co-ordinate
!
!       In this case, the following must be set:
!
!       H%type( 1 : 10 ) = TRANSFER( 'COORDINATE', H%type )
!       H%val( : )   the values of the components of H
!       H%row( : )   the row indices of the components of H
!       H%col( : )   the column indices of the components of H
!       H%ne         the number of nonzeros used to store 
!                    the LOWER TRIANGULAR part of H
!
!    ii) sparse, by rows
!
!       In this case, the following must be set:
!
!       H%type( 1 : 14 ) = TRANSFER( 'SPARSE_BY_ROWS', H%type )
!       H%val( : )   the values of the components of H, stored row by row
!       H%col( : )   the column indices of the components of H
!       H%ptr( : )   pointers to the start of each row, and past the end of
!                    the last row
!
!    iii) dense, by rows
!
!       In this case, the following must be set:
!
!       H%type( 1 : 5 ) = TRANSFER( 'DENSE', H%type )
!       H%val( : )   the values of the components of H, stored row by row,
!                    with each the entries in each row in order of 
!                    increasing column indicies.
!
!    iv) diagonal
!
!       In this case, the following must be set:
!
!       H%type( 1 : 8 ) = TRANSFER( 'DIAGONAL', H%type )
!       H%val( : )   the values of the diagonals of H, stored in order
!
!    On exit, the components will most likely have been reordered.
!    The output  matrix will be stored by rows, according to scheme (ii) above.
!    However, if scheme (i) is used for input, the output H%row will contain
!    the row numbers corresponding to the values in H%val, and thus in this
!    case the output matrix will be available in both formats (i) and (ii).
!    
!   %G is a REAL array of length %n, which must be set by
!    the user to the value of the gradient, g, of the linear term of the
!    quadratic objective function. The i-th component of G, i = 1, ....,
!    n, should contain the value of g_i.  
!    On exit, G will most likely have been reordered.
!   
!   %f is a REAL variable, which must be set by the user to the value of
!    the constant term f in the objective function. On exit, it may have
!    been changed to reflect variables which have been fixed.
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
!   %C is a REAL array of length %m, which is used to store the values of 
!    A x. It need not be set on entry. On exit, it will have been filled 
!    with appropriate values.
!
!   %X is a REAL array of length %n, which must be set by the user
!    to an estimate of the solution x. On successful exit, it will contain
!    the required solution.
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
!  C_stat is a INTEGER array of length m, which may be set by the user
!   on entry to QPA_solve to indicate which of the constraints are to
!   be included in the initial working set. If this facility is required,
!   the component control%cold_start must be set to 0 on entry; C_stat
!   need not be set if control%cold_start is nonzero. On exit,
!   C_stat will indicate which constraints are in the final working set.
!   Possible entry/exit values are 
!   C_stat( i ) < 0, the i-th constraint is in the working set, 
!                    on its lower bound, 
!               > 0, the i-th constraint is in the working set
!                    on its upper bound, and
!               = 0, the i-th constraint is not in the working set
!
!  B_stat is a INTEGER array of length n, which may be set by the user
!   on entry to QPA_solve to indicate which of the simple bound constraints 
!   are to be included in the initial working set. If this facility is required,
!   the component control%cold_start must be set to 0 on entry; B_stat
!   need not be set if control%cold_start is nonzero. On exit,
!   B_stat will indicate which constraints are in the final working set.
!   Possible entry/exit values are 
!   B_stat( i ) < 0, the i-th bound constraint is in the working set, 
!                    on its lower bound, 
!               > 0, the i-th bound constraint is in the working set
!                    on its upper bound, and
!               = 0, the i-th bound constraint is not in the working set
!
!  data is a structure of type QPC_data_type which holds private internal data
!
!  control is a structure of type QPC_control_type that controls the 
!   execution of the subroutine and must be set by the user. Default values for
!   the elements may be set by a call to QPC_initialize. See QPC_initialize 
!   for details
!
!  inform is a structure of type QPC_inform_type that provides 
!    information on exit from QPC_solve. The component status 
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
!        prob%H%type in { 'DENSE', 'SPARSE_BY_ROWS', 'COORDINATE', , 'DIAGONAL' }
!       has been violated.
!
!    -4 The bound constraints are inconsistent.
!
!    -5 The constraints appear to have no feasible point.
!
!    -7 The objective function appears to be unbounded from below on the
!       feasible set.
!
!    -9 The factorization failed; the return status from the factorization
!       package is given in the component factorization_status.
!      
!    -13 The problem is so ill-conditoned that further progress is impossible.  
!
!    -16 The step is too small to make further impact.
!
!    -17 Too many iterations have been performed. This may happen if
!       control%maxit is too small, but may also be symptomatic of 
!       a badly scaled problem.
!
!    -18 Too much CPU time has passed. This may happen if control%cpu_time_limit 
!        is too small, but may also be symptomatic of a badly scaled problem.
!
!    -23 an entry from the strict upper triangle of H has been input.
!
!  On exit from QPC_solve, other components of inform give the 
!  following:
!
!     alloc_status = The status of the last attempted allocation/deallocation 
!     factorization_integer = The total integer workspace required for the 
!       factorization.
!     factorization_real = The total real workspace required for the 
!       factorization.
!     nfacts = The total number of factorizations performed.
!     nmods = The total number of factorizations which were modified to 
!       ensure that the matrix was an appropriate preconditioner. 
!     factorization_status = the return status from the matrix factorization
!       package.   
!     obj = the value of the objective function at the best estimate of the 
!       solution determined by QPC_solve.
!     non_negligible_pivot = the smallest pivot which was not judged to be
!       zero when detecting linearly dependent constraints
!     bad_alloc = the name of the array for which an allocation/deallocation
!       error ocurred
!     time%total = the total time spent in the package.
!     time%preprocess = the time spent preprocessing the problem.
!     time%find_dependent = the time spent detecting linear dependencies
!     time%analyse = the time spent analysing the required matrices prior to
!       factorization.
!     time%factorize = the time spent factorizing the required matrices.
!     time%solve = the time spent computing the search direction.
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      INTEGER, INTENT( INOUT ), DIMENSION( prob%m ) :: C_stat
      INTEGER, INTENT( INOUT ), DIMENSION( prob%n ) :: B_stat
      TYPE ( QPC_data_type ), INTENT( INOUT ) :: data
      TYPE ( QPC_control_type ), INTENT( INOUT ) :: control
      TYPE ( QPC_inform_type ), INTENT( OUT ) :: inform

      REAL ( KIND = wp ), INTENT( IN ), OPTIONAL, DIMENSION( prob%n ) :: G_p
      REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL, DIMENSION( prob%n ) :: X_p
      REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL, DIMENSION( prob%n ) :: Z_p
      REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL, DIMENSION( prob%m ) :: Y_p

!  Local variables

      INTEGER :: a_ne, h_ne, i, j, l, tiny_x, tiny_c, n_depen, n_more_depen, nzc
      INTEGER :: ii, lbreak, k_n_max, lbd, m_link
      INTEGER :: hd_start, hd_end, hnd_start, hnd_end, type, n_pcg
      REAL :: dum, time, time_start, time_qpa_start
      REAL ( KIND = wp ) :: tol, f, av_bnd, a_x, a_norms, best_obj
      LOGICAL :: printi, printd, printt, first_pass, center, reset_bnd, lsqp
      LOGICAL :: remap_fixed, remap_freed, remap_more_freed
      LOGICAL :: diagonal_qp, convex_diagonal_qp, gotsol
      CHARACTER ( LEN = 80 ) :: array_name
      TYPE ( FDC_control_type ) :: FDC_control        
      TYPE ( FDC_inform_type ) :: FDC_inform
      TYPE ( LSQP_control_type ) :: LSQP_control
      TYPE ( LSQP_inform_type ) :: LSQP_inform

!  prefix for all output 

      CHARACTER ( LEN = LEN( TRIM( control%prefix ) ) - 2 ) :: prefix
      prefix = control%prefix( 2 : LEN( TRIM( control%prefix ) ) - 1 )

      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( A, ' -- entering QPC_solve ' )" ) prefix

!  Initialize time

      CALL CPU_TIME( time_start )

!  Set initial timing breakdowns

      inform%time%total = 0.0 ; inform%time%preprocess = 0.0 
      inform%time%analyse = 0.0 ; inform%time%factorize = 0.0 
      inform%time%solve = 0.0

      inform%QPA_inform%time%total = 0.0
      inform%QPA_inform%time%preprocess = 0.0
      inform%QPA_inform%time%analyse = 0.0
      inform%QPA_inform%time%factorize = 0.0 
      inform%QPA_inform%time%solve = 0.0

      inform%QPB_inform%time%total = 0.0 
      inform%QPB_inform%time%preprocess = 0.0
      inform%QPB_inform%time%analyse = 0.0
      inform%QPB_inform%time%factorize = 0.0 
      inform%QPB_inform%time%solve = 0.0
      inform%QPB_inform%time%phase1_total = 0.0
      inform%QPB_inform%time%phase1_analyse = 0.0
      inform%QPB_inform%time%phase1_factorize = 0.0 
      inform%QPB_inform%time%phase1_solve = 0.0

!  Initialize counts

      inform%status = - 31 ; inform%alloc_status = 0 ; inform%bad_alloc = ''
      inform%nfacts = 0 ; inform%nmods = 0
      inform%factorization_integer = 0 ; inform%factorization_real = 0
      inform%obj = infinity
      inform%p_found = .FALSE.
      gotsol = .FALSE.

      inform%QPA_inform%status = - 31 ; inform%QPA_inform%major_iter = 0
      inform%QPA_inform%iter = 0 ; inform%QPA_inform%nfacts = 0
      inform%QPA_inform%cg_iter = 0; inform%QPA_inform%bad_alloc = ''
      inform%QPA_inform%nmods = 0 ; inform%QPA_inform%alloc_status = 0
      inform%QPA_inform%factorization_integer = 0 
      inform%QPA_inform%factorization_real = 0
      inform%QPA_inform%obj = infinity ; inform%QPA_inform%merit = infinity
      inform%QPA_inform%num_g_infeas = - 1
      inform%QPA_inform%num_b_infeas = - 1
      inform%QPA_inform%infeas_g = infinity
      inform%QPA_inform%infeas_b = infinity

      inform%QPB_inform%status = - 31 ; inform%QPB_inform%bad_alloc = ''
      inform%QPB_inform%iter = 0 ; inform%QPB_inform%nfacts = 0
      inform%QPB_inform%cg_iter = 0 ; inform%QPB_inform%nbacts = 0 
      inform%QPB_inform%nmods = 0 ; inform%QPB_inform%alloc_status = 0
      inform%QPB_inform%factorization_integer = 0 
      inform%QPB_inform%factorization_real = 0
      inform%QPB_inform%obj = infinity
      inform%QPB_inform%feasible = .FALSE.

!  Basic single line of output per iteration

      printi = control%out > 0 .AND. control%print_level >= 1 
      printt = control%out > 0 .AND. control%print_level >= 2
      printd = control%out > 0 .AND. control%print_level >= 11

!  Make sure that QPA, QPB and LSQP control parameters follow those from QPC

      IF ( control%qpb_or_qpa .OR. .NOT. control%no_qpa ) THEN
        control%QPA_control%indmin = control%indmin
        control%QPA_control%valmin = control%valmin
        control%QPA_control%restore_problem = control%restore_problem
        control%QPA_control%infinity = control%infinity 
        control%QPA_control%treat_zero_bounds_as_general =                     &
          control%treat_zero_bounds_as_general 
        control%QPA_control%array_syntax_worse_than_do_loop =                  &
          control%array_syntax_worse_than_do_loop 
        control%QPA_control%solve_qp = .TRUE.
        IF ( control%QPA_control%print_level <= 0 ) THEN
          data%CNTLA%ldiag = 0 ; data%CNTLA%sp = - 1 
          data%CNTLA%lp = - 1 ; data%CNTLA%wp = - 1 ; data%CNTLA%mp = - 1
        END IF

        IF ( control%cpu_time_limit >= zero ) THEN
          IF ( control%QPA_control%cpu_time_limit >= zero ) THEN
            control%QPA_control%cpu_time_limit = MIN( control%cpu_time_limit,  &
               control%QPA_control%cpu_time_limit )
          ELSE
            control%QPA_control%cpu_time_limit = control%cpu_time_limit
          END IF
        END IF
      END IF

      IF ( control%qpb_or_qpa .OR. .NOT. control%no_qpb ) THEN
        control%QPB_control%indmin = control%indmin
        control%QPB_control%valmin = control%valmin
        control%QPB_control%restore_problem = control%restore_problem
        control%QPB_control%infinity = control%infinity 
        control%QPB_control%identical_bounds_tol =                             &
          control%identical_bounds_tol
        control%QPB_control%treat_zero_bounds_as_general =                     &
          control%treat_zero_bounds_as_general 
        control%QPB_control%array_syntax_worse_than_do_loop =                  &
          control%array_syntax_worse_than_do_loop 

        control%QPB_control%LSQP_control%indmin = control%indmin
        control%QPB_control%LSQP_control%valmin = control%valmin
        control%QPB_control%LSQP_control%restore_problem =                     &
          control%restore_problem
        control%QPB_control%LSQP_control%infinity = control%infinity 
        control%QPB_control%LSQP_control%identical_bounds_tol =                &
          control%identical_bounds_tol
        control%QPB_control%LSQP_control%treat_zero_bounds_as_general =        &
          control%treat_zero_bounds_as_general 
        control%QPB_control%LSQP_control%array_syntax_worse_than_do_loop =     &
          control%array_syntax_worse_than_do_loop 
        IF ( control%QPB_control%print_level <= 0 ) THEN
          data%CNTL%ldiag = 0 ; data%CNTL%sp = - 1 
          data%CNTL%lp = - 1 ; data%CNTL%wp = - 1 ; data%CNTL%mp = - 1
        END IF

        IF ( control%cpu_time_limit >= zero ) THEN
          IF ( control%QPB_control%cpu_time_limit >= zero ) THEN
            control%QPB_control%cpu_time_limit = MIN( control%cpu_time_limit,  &
               control%QPB_control%cpu_time_limit )
          ELSE
            control%QPB_control%cpu_time_limit = control%cpu_time_limit
          END IF
          IF ( control%QPB_control%LSQP_control%cpu_time_limit >= zero ) THEN
            control%QPB_control%LSQP_control%cpu_time_limit =                  &
              MIN( control%cpu_time_limit,                                     &
                   control%QPB_control%LSQP_control%cpu_time_limit )
          ELSE
            control%QPB_control%LSQP_control%cpu_time_limit =                  &
              control%cpu_time_limit
          END IF
        END IF
      END IF

!  Ensure that input parameters are within allowed ranges

      IF ( prob%n <= 0 .OR. prob%m < 0 .OR.                                    &
           .NOT. QPT_keyword_H( prob%H%type ) .OR.                             &
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
        WRITE( control%out, "( ' G = ', /, ( 5ES12.4 ) )" ) prob%G( : prob%n )
        IF ( SMT_get( prob%H%type ) == 'DIAGONAL' ) THEN
          WRITE( control%out, "( ' H (diagonal) = ', /, ( 5ES12.4 ) )" )       &
            prob%H%val( : prob%n )
        ELSE IF ( SMT_get( prob%H%type ) == 'DENSE' ) THEN
          WRITE( control%out, "( ' H (dense) = ', /, ( 5ES12.4 ) )" )          &
            prob%H%val( : prob%n * ( prob%n + 1 ) / 2 )
        ELSE IF ( SMT_get( prob%H%type ) == 'SPARSE_BY_ROWS' ) THEN
          WRITE( control%out, "( ' H (row-wise) = ' )" )
          DO i = 1, prob%m
            WRITE( control%out, "( ( 2( 2I8, ES12.4 ) ) )" )                   &
              ( i, prob%H%col( j ), prob%H%val( j ),                           &
                j = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1 )
          END DO
        ELSE
          WRITE( control%out, "( ' H (co-ordinate) = ' )" )
          WRITE( control%out, "( ( 2( 2I8, ES12.4 ) ) )" )                     &
          ( prob%H%row( i ), prob%H%col( i ), prob%H%val( i ), i = 1, prob%H%ne )
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
        ELSE IF ( prob%X_u( i ) == prob%X_l( i ) ) THEN
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

!  Set Hessian and gradient types to generic - this may change in future

      prob%Hessian_kind = - 1 ; prob%gradient_kind = - 1

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
          control%QPB_control%treat_zero_bounds_as_general

!  Store the problem dimensions

        IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
          a_ne = prob%m * prob%n
        ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
          a_ne = prob%A%ptr( prob%m + 1 ) - 1
        ELSE
          a_ne = prob%A%ne 
        END IF

        IF ( SMT_get( prob%H%type ) == 'DIAGONAL' ) THEN
          h_ne = prob%n
        ELSE IF ( SMT_get( prob%H%type ) == 'DENSE' ) THEN
          h_ne = ( prob%n * ( prob%n + 1 ) ) / 2
        ELSE IF ( SMT_get( prob%H%type ) == 'SPARSE_BY_ROWS' ) THEN
          h_ne = prob%H%ptr( prob%n + 1 ) - 1
        ELSE
          h_ne = prob%H%ne 
        END IF

        IF ( printi ) WRITE( control%out,                                      &
               "( /, A, ' problem dimensions before preprocessing: ', /,  A,   &
     &         ' n = ', I8, ' m = ', I8, ' a_ne = ', I8, ' h_ne = ', I8 )" )   &
               prefix, prefix, prob%n, prob%m, a_ne, h_ne

!  Perform the preprocessing

        CALL CPU_TIME( time ) 
        CALL QPP_reorder( data%QPP_map, data%QPP_control,                      &
                          data%QPP_inform, data%dims, prob,                    &
                          .FALSE., .FALSE., .FALSE. )
        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum

!  Test for satisfactory termination

        IF ( data%QPP_inform%status /= GALAHAD_ok ) THEN
          inform%status = data%QPP_inform%status
          IF ( control%out > 0 .AND. control%print_level >= 1 )                &
            WRITE( control%out, "( A, ' status ', I3, ' after QPP_reorder')" ) &
             prefix, data%QPP_inform%status
          IF ( control%error > 0 .AND. control%print_level > 0 )               &
            WRITE( control%error, 2010 ) inform%status 
          IF ( control%out > 0 .AND. control%print_level > 0 .AND.             &
               inform%status == GALAHAD_error_upper_entry )                    &
            WRITE( control%error, 2240 ) 
          CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )
          CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
          IF ( data%QPP_inform%status == GALAHAD_error_primal_infeasible ) THEN
            inform%status = GALAHAD_error_primal_infeasible
            GO TO 800 
          ELSE
            GO TO 800 
          END IF 
        END IF 

!  Record array lengths

        IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
          a_ne = prob%m * prob%n
        ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
          a_ne = prob%A%ptr( prob%m + 1 ) - 1
        ELSE
          a_ne = prob%A%ne 
        END IF

        IF ( SMT_get( prob%H%type ) == 'DIAGONAL' ) THEN
          h_ne = prob%n
        ELSE IF ( SMT_get( prob%H%type ) == 'DENSE' ) THEN
          h_ne = ( prob%n * ( prob%n + 1 ) ) / 2
        ELSE IF ( SMT_get( prob%H%type ) == 'SPARSE_BY_ROWS' ) THEN
          h_ne = prob%H%ptr( prob%n + 1 ) - 1
        ELSE
          h_ne = prob%H%ne 
        END IF

        IF ( printi ) WRITE( control%out,                                      &
               "(  A, ' problem dimensions after preprocessing: ', /,  A,      &
     &         ' n = ', I8, ' m = ', I8, ' a_ne = ', I8, ' h_ne = ', I8 )" )   &
               prefix, prefix, prob%n, prob%m, a_ne, h_ne

        prob%new_problem_structure = .FALSE.
        data%trans = 1

!  Recover the problem dimensions after preprocessing

      ELSE
        IF ( data%trans == 0 ) THEN
          CALL CPU_TIME( time ) 
          CALL QPP_apply( data%QPP_map, data%QPP_inform, prob,                 &
                          get_all = .TRUE. )
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
            IF ( data%QPP_inform%status == GALAHAD_error_primal_infeasible ) THEN
              inform%status = GALAHAD_error_primal_infeasible
              GO TO 800 
            ELSE
              GO TO 800 
            END IF 
          END IF 
        END IF 
        data%trans = data%trans + 1

!  Record array lengths

        IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
          a_ne = prob%m * prob%n
        ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
          a_ne = prob%A%ptr( prob%m + 1 ) - 1
        ELSE
          a_ne = prob%A%ne 
        END IF

        IF ( SMT_get( prob%H%type ) == 'DIAGONAL' ) THEN
          h_ne = prob%n
        ELSE IF ( SMT_get( prob%H%type ) == 'DENSE' ) THEN
          h_ne = ( prob%n * ( prob%n + 1 ) ) / 2
        ELSE IF ( SMT_get( prob%H%type ) == 'SPARSE_BY_ROWS' ) THEN
          h_ne = prob%H%ptr( prob%n + 1 ) - 1
        ELSE
          h_ne = prob%H%ne 
        END IF
      END IF

!  Permute initial working sets if provided

      IF ( control%QPA_control%cold_start == 0 ) THEN
        CALL SORT_inplace_permute( data%QPP_map%m, data%QPP_map%c_map,         &
                                   IX = C_stat( : data%QPP_map%m ) )
        CALL SORT_inplace_permute( data%QPP_map%n, data%QPP_map%x_map,         &
                                   IX = B_stat( : data%QPP_map%n ) )
      END IF

      remap_freed = .FALSE. ; remap_more_freed = .FALSE. ; remap_fixed = .FALSE.
      IF ( .NOT. control%qpb_or_qpa .AND. control%no_qpb ) GO TO 150

!  =================================================================
!  Check to see if the equality constraints are linearly independent
!  =================================================================

      IF ( prob%m > 0 .AND.                                                    &
           ( .NOT. data%tried_to_remove_deps .AND.                             &
              control%QPB_control%remove_dependencies ) ) THEN

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
        FDC_control%max_infeas = control%QPB_control%stop_p
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

        IF ( FDC_inform%status /= GALAHAD_ok ) THEN

!  Allocate arrays to hold the matrix vector product

          array_name = 'qpc: data%HX'
          CALL SPACE_resize_array( prob%n, data%HX, inform%status,             &
                 inform%alloc_status, array_name = array_name,                 &
                 deallocate_error_fatal = control%deallocate_error_fatal,      &
                 exact_size = control%space_critical,                          &
                 bad_alloc = inform%bad_alloc, out = control%error )
          IF ( inform%status /= GALAHAD_ok ) GO TO 900

!  On error exit, compute the current objective function value

          data%HX( : prob%n ) = zero
          CALL QPC_HX( data%dims, prob%n, data%HX( : prob%n ),                 &
                       prob%H%ptr( prob%n + 1 ) - 1, prob%H%val,               &
                       prob%H%col, prob%H%ptr, prob%X( : prob%n ), '+' )
          inform%obj = half * DOT_PRODUCT( prob%X( : prob%n ),                 &
                                           data%HX( : prob%n ) )               &
                       + DOT_PRODUCT( prob%X( : prob%n ),                      &
                                      prob%G( : prob%n ) ) + prob%f

!  Print details of the error exit

          IF ( FDC_inform%status == GALAHAD_error_analysis ) THEN
            inform%status = GALAHAD_error_analysis
            inform%factorization_status = GALAHAD_error_upper_entry 
          ELSE IF ( FDC_inform%status == GALAHAD_error_factorization ) THEN
            inform%status = GALAHAD_error_factorization
            inform%factorization_status = GALAHAD_error_bad_bounds
          ELSE
            inform%status = FDC_inform%status
          END IF
          IF ( control%error > 0 .AND. control%print_level >= 1 ) THEN
            WRITE( control%out, "( ' ' )" )
            WRITE( control%error, 2040 ) inform%status, 'FDC_find_dependent'
          END IF
          GO TO 750
        END IF

        IF ( control%out > 0 .AND. control%print_level >= 2 .AND. n_depen > 0 )&
          WRITE( control%out, "(/, ' The following ',I0,' constraints appear', &
       &         ' to be dependent', /, ( 8I8 ) )" ) n_depen, data%Index_C_freed

        remap_freed = n_depen > 0 .AND. prob%n > 0

!  Special case: no free variables

        IF ( prob%n == 0 ) THEN
          IF ( control%QPB_control%array_syntax_worse_than_do_loop ) THEN
            DO i = 1, prob%m ; prob%Y( i ) = zero ; END DO
            DO i = 1, prob%n ; prob%Z( i ) = prob%G( i ) ; END DO
          ELSE
            prob%Y( : prob%m ) = zero
            prob%Z( : prob%n ) = prob%G( : prob%n )
          END IF
          CALL QPC_HX( data%dims, prob%n, prob%Z( : prob%n ),                  &
                       prob%H%ptr( prob%n + 1 ) - 1, prob%H%val,               &
                       prob%H%col, prob%H%ptr, prob%X( : prob%n ), '+' )
          prob%C( : prob%m ) = zero
          CALL QPC_AX( prob%m, prob%C( : prob%m ), prob%m,                     &
                       prob%A%ptr( prob%m + 1 ) - 1, prob%A%val,               &
                       prob%A%col, prob%A%ptr, prob%n, prob%X, '+ ')
          remap_more_freed = .FALSE. ; remap_fixed = .FALSE.
          inform%obj = prob%f
          inform%status = GALAHAD_ok
          GO TO 700
        END IF
        data%tried_to_remove_deps = .TRUE.
      ELSE
        remap_freed = .FALSE.
      END IF

      IF ( remap_freed ) THEN

!  Some of the current constraints will be removed by freeing them

        IF ( control%error > 0 .AND. control%print_level >= 1 )                &
          WRITE( control%out, "( /, ' -> ', i6, ' constraints are',            &
         & ' dependent and will be temporarily removed' )" ) n_depen

!  Allocate arrays to indicate which constraints have been freed

        array_name = 'qpc: data%C_freed'
        CALL SPACE_resize_array( n_depen, data%C_freed, inform%status,         &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= GALAHAD_ok ) GO TO 900

!  Free the constraint bounds as required

        DO i = 1, n_depen
          j = data%Index_C_freed( i )
          data%C_freed( i ) = prob%C_l( j )
          prob%C_l( j ) = - control%infinity
          prob%C_u( j ) = control%infinity
          prob%Y( j ) = zero
        END DO

        CALL QPP_initialize( data%QPP_map_freed, data%QPP_control )
        data%QPP_control%infinity = control%infinity
        data%QPP_control%treat_zero_bounds_as_general =                       &
          control%QPB_control%treat_zero_bounds_as_general

!  Store the problem dimensions

        data%dims_save_freed = data%dims
        a_ne = prob%A%ne 
        h_ne = prob%H%ne 

        IF ( printi ) WRITE( control%out,                                      &
               "( /, A, ' problem dimensions before removal of dependencies:', &
     &            /, A,                                                        &
     &         ' n = ', I8, ' m = ', I8, ' a_ne = ', I8, ' h_ne = ', I8 )" )   &
               prefix, prefix, prob%n, prob%m, a_ne, h_ne

!  Perform the preprocessing

        CALL CPU_TIME( time ) 
        CALL QPP_reorder( data%QPP_map_freed, data%QPP_control,                &
                          data%QPP_inform, data%dims, prob,                    &
                          .FALSE., .FALSE., .FALSE. )
        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum
        inform%QPB_inform%time%preprocess = inform%time%preprocess
        inform%QPA_inform%time%preprocess = inform%time%preprocess
  
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
          IF ( control%out > 0 .AND. control%print_level >= 1 )                &
            WRITE( control%out, "( A, ' status ', I3, ' after QPP_reorder')" ) &
             prefix, data%QPP_inform%status
          IF ( control%error > 0 .AND. control%print_level > 0 )               &
            WRITE( control%error, 2010 ) inform%status 
          IF ( control%out > 0 .AND. control%print_level > 0 .AND.             &
               inform%status == GALAHAD_error_upper_entry )                    &
            WRITE( control%error, 2240 ) 
          CALL QPP_terminate( data%QPP_map_freed, data%QPP_control,            &
                              data%QPP_inform )
          CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )
          CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
          IF ( data%QPP_inform%status == - GALAHAD_error_primal_infeasible ) THEN
            inform%status = GALAHAD_error_primal_infeasible
            GO TO 800 
          ELSE
            GO TO 800 
          END IF 
        END IF 

!  Record revised array lengths

        IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
          a_ne = prob%m * prob%n
        ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
          a_ne = prob%A%ptr( prob%m + 1 ) - 1
        ELSE
          a_ne = prob%A%ne 
        END IF

        IF ( SMT_get( prob%H%type ) == 'DIAGONAL' ) THEN
          h_ne = prob%n
        ELSE IF ( SMT_get( prob%H%type ) == 'DENSE' ) THEN
          h_ne = ( prob%n * ( prob%n + 1 ) ) / 2
        ELSE IF ( SMT_get( prob%H%type ) == 'SPARSE_BY_ROWS' ) THEN
          h_ne = prob%H%ptr( prob%n + 1 ) - 1
        ELSE
          h_ne = prob%H%ne 
        END IF

        IF ( printi ) WRITE( control%out,                                      &
          "( A, ' problem dimensions after removal of dependencies:', /, A,    &
     &          ' n = ', I8, ' m = ', I8, ' a_ne = ', I8, ' h_ne = ', I8 )" )  &
             prefix, prefix, prob%n, prob%m, a_ne, h_ne

!  Further permute initial working sets if provided

        IF ( control%QPA_control%cold_start == 0 ) THEN
          CALL SORT_inplace_permute( data%QPP_map_freed%m,                     &
                                     data%QPP_map_freed%c_map,                 &
                                     IX = C_stat( : data%QPP_map_freed%m ) )
          CALL SORT_inplace_permute( data%QPP_map_freed%n,                     &
                                     data%QPP_map_freed%x_map,                 &
                                     IX = B_stat( : data%QPP_map_freed%n ) )
        END IF
      END IF

!  Special case: Bound-constrained QP

      IF ( a_ne == 0 .AND. h_ne /= 0 ) THEN

!  Check to see if the Hessian is diagonal

        diagonal_qp = .TRUE.
        SELECT CASE ( SMT_get( prob%H%type ) )
        CASE ( 'DIAGONAL' ) 
        CASE ( 'DENSE' ) 
          l = 0
          DO i = 1, prob%n
            DO j = 1, i
              l = l + 1
              IF ( i /= j .AND. prob%H%val( l ) /= zero ) THEN
                diagonal_qp = .FALSE. ; GO TO 3
              END IF
            END DO
          END DO
        CASE ( 'SPARSE_BY_ROWS' )
          DO i = 1, prob%n
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              j = prob%H%col( l )
              IF ( i /= j .AND. prob%H%val( l ) /= zero ) THEN
                diagonal_qp = .FALSE. ; GO TO 3
              END IF
            END DO
          END DO
        CASE ( 'COORDINATE' )
          DO l = 1, prob%H%ne
            i = prob%H%row( l ) ; j = prob%H%col( l )
            IF ( i /= j .AND. prob%H%val( l ) /= zero ) THEN
              diagonal_qp = .FALSE. ; GO TO 3
            END IF
          END DO
        END SELECT
  3     CONTINUE

        remap_more_freed = .FALSE. ; remap_fixed = .FALSE.
        IF ( diagonal_qp ) THEN
          IF ( printi ) WRITE( control%out,                                    &
            "( /, A, ' solving separable bound-constrained QP ' )" ) prefix
          CALL QPB_optimal_for_SBQP( prob, control%QPB_control,                &
                                     inform%QPB_inform,                        &
                                     B_stat = B_stat( : prob%n ) )
          IF ( printi ) WRITE( control%out,                                    &
              "( A, ' on exit from QPB_optimal_for_SBQP: status = ', I0,       &
           &   ', time = ', F0.2, /, A, ' objective value =', ES12.4,          &
           &   /, A, ' active bounds: ', I0, ' from ', I0 )" )                 &
              prefix, inform%QPB_inform%status, inform%QPB_inform%time%total,  &
              prefix, inform%QPB_inform%obj, prefix,                           &
              COUNT( B_stat( : prob%n ) /= 0 ), prob%n
          inform%alloc_status = inform%QPB_inform%alloc_status
          inform%nfacts = inform%nfacts + inform%QPB_inform%nfacts
          inform%nmods = inform%nmods + inform%QPB_inform%nmods
          inform%obj = inform%QPB_inform%obj
          inform%time%total = inform%time%total + inform%QPB_inform%time%total
          inform%status = inform%QPB_inform%status
          GO TO 700
        ELSE
          CALL QPB_feasible_for_BQP( prob, data, control%QPB_control,          &
                                     inform%QPB_inform )
        END IF

!  General case: QP or LP

      ELSE
        first_pass = .TRUE.
        center = control%QPB_control%center
        f = prob%f
    
!  Check to see if the Hessian is diagonal and positive semi-definite

        convex_diagonal_qp = .TRUE.
        SELECT CASE ( SMT_get( prob%H%type ) )
        CASE ( 'DIAGONAL' ) 
          IF ( COUNT( prob%H%val( : prob%n ) < zero ) > 0 )                    &
            convex_diagonal_qp = .FALSE.
        CASE ( 'DENSE' ) 
          l = 0
          DO i = 1, prob%n
            DO j = 1, i
              l = l + 1
              IF ( ( i /= j .AND. prob%H%val( l ) /= zero ) .OR.               &
                   ( i == j .AND. prob%H%val( l ) <  zero ) ) THEN
                convex_diagonal_qp = .FALSE. ; GO TO 5
              END IF
            END DO
          END DO
        CASE ( 'SPARSE_BY_ROWS' )
          DO i = 1, prob%n
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              j = prob%H%col( l )
              IF ( ( i /= j .AND. prob%H%val( l ) /= zero ) .OR.               &
                   ( i == j .AND. prob%H%val( l ) <  zero ) ) THEN
                convex_diagonal_qp = .FALSE. ; GO TO 5
              END IF
            END DO
          END DO
        CASE ( 'COORDINATE' )
          DO l = 1, prob%H%ne
            i = prob%H%row( l ) ; j = prob%H%col( l )
            IF ( ( i /= j .AND. prob%H%val( l ) /= zero ) .OR.                 &
                 ( i == j .AND. prob%H%val( l ) <  zero ) ) THEN
              convex_diagonal_qp = .FALSE. ; GO TO 5
            END IF
          END DO
        END SELECT
  5     CONTINUE

        IF ( h_ne == 0 .OR. convex_diagonal_qp ) THEN
          prob%gradient_kind = 2
          IF ( TWO_NORM( prob%n, prob%G, 1 ) <= epsmch ) THEN
            prob%gradient_kind = 0
          ELSE IF ( TWO_NORM( prob%n, prob%G - one, 1 ) <= epsmch ) THEN
            prob%gradient_kind = 1
          END IF
        ELSE
          prob%gradient_kind = 0
        END IF

!  ==============================
!  Find an initial feasible point
!  ==============================

 10     CONTINUE

!  Set appropraiate control parameters for the phase 1

        LSQP_control = control%QPB_control%LSQP_control
        LSQP_control%print_level = control%QPB_control%print_level
        LSQP_control%out = control%QPB_control%out
        LSQP_control%error = control%QPB_control%error
        LSQP_control%factor = control%QPB_control%factor
        LSQP_control%max_col = control%QPB_control%max_col
        LSQP_control%itref_max = control%QPB_control%itref_max
        LSQP_control%infeas_max = control%QPB_control%infeas_max
        LSQP_control%restore_problem = control%QPB_control%restore_problem 
        LSQP_control%infinity = control%infinity
        LSQP_control%stop_p = control%QPB_control%stop_p
        LSQP_control%stop_c = control%QPB_control%stop_c
        LSQP_control%stop_d = control%QPB_control%stop_d
        LSQP_control%muzero = control%QPB_control%muzero
        LSQP_control%reduce_infeas = control%QPB_control%reduce_infeas
        LSQP_control%pivot_tol = control%QPB_control%pivot_tol
        LSQP_control%pivot_tol_for_dependencies =                              &
          control%QPB_control%pivot_tol_for_dependencies
        LSQP_control%zero_pivot = control%QPB_control%zero_pivot
        LSQP_control%identical_bounds_tol = control%identical_bounds_tol
        LSQP_control%remove_dependencies =                                     &
          control%QPB_control%remove_dependencies
        LSQP_control%treat_zero_bounds_as_general =                            &
          control%treat_zero_bounds_as_general
        LSQP_control%feasol = .FALSE.
        LSQP_control%array_syntax_worse_than_do_loop =                         &
          control%array_syntax_worse_than_do_loop
  
!  Either find the solution to the LP (if there is no Hessian term) ...

        IF ( h_ne == 0 ) THEN
          lsqp = .TRUE.
          LSQP_control%just_feasible = .FALSE.
          LSQP_control%maxit = control%QPB_control%maxit
          LSQP_control%prfeas =                                                &
            MAX( LSQP_control%prfeas, control%QPB_control%prfeas )
          LSQP_control%dufeas =                                                &
            MAX( LSQP_control%dufeas, control%QPB_control%dufeas )
          prob%Hessian_kind = 0

!  .. or find the solution to the Diagonal QP (if the Hessian is diagonal) ...

        ELSE IF ( convex_diagonal_qp ) THEN
          lsqp = .TRUE.
          LSQP_control%just_feasible = .FALSE.
          LSQP_control%maxit = control%QPB_control%maxit
          LSQP_control%prfeas =                                                &
            MAX( LSQP_control%prfeas, control%QPB_control%prfeas )
          LSQP_control%dufeas =                                                &
            MAX( LSQP_control%dufeas, control%QPB_control%dufeas )
          prob%Hessian_kind = 2

          array_name = 'qpc: prob%X0'
          CALL SPACE_resize_array( prob%n, prob%X0, inform%status,             &
                 inform%alloc_status, array_name = array_name,                 &
                 deallocate_error_fatal = control%deallocate_error_fatal,      &
                 exact_size = control%space_critical,                          &
                 bad_alloc = inform%bad_alloc, out = control%error )
          IF ( inform%status /= GALAHAD_ok ) GO TO 900

          array_name = 'qpc: prob%WEIGHT'
          CALL SPACE_resize_array( prob%n, prob%WEIGHT, inform%status,         &
                 inform%alloc_status, array_name = array_name,                 &
                 deallocate_error_fatal = control%deallocate_error_fatal,      &
                 exact_size = control%space_critical,                          &
                 bad_alloc = inform%bad_alloc, out = control%error )
          IF ( inform%status /= GALAHAD_ok ) GO TO 900

          SELECT CASE ( SMT_get( prob%H%type ) )
          CASE ( 'DIAGONAL' ) 
            prob%WEIGHT( : prob%n ) = prob%H%val( : prob%n )
          CASE ( 'DENSE' ) 
            l = 0
            DO i = 1, prob%n
              DO j = 1, i
                l = l + 1
                IF ( i == j ) prob%WEIGHT( i ) = prob%H%val( l )
              END DO
            END DO
          CASE ( 'SPARSE_BY_ROWS' )
            prob%WEIGHT( : prob%n ) = zero
            DO i = 1, prob%n
              DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
                j = prob%H%col( l )
                IF ( i == j )                                                  &
                  prob%WEIGHT( i ) = prob%WEIGHT( i ) + prob%H%val( l )
              END DO
            END DO
          CASE ( 'COORDINATE' )
            prob%WEIGHT( : prob%n ) = zero
            DO l = 1, prob%H%ne
              i = prob%H%row( l ) ; j = prob%H%col( l )
              IF ( i == j )                                                    &
                prob%WEIGHT( i ) = prob%WEIGHT( i ) + prob%H%val( l )
            END DO
          END SELECT
          prob%WEIGHT( : prob%n ) = SQRT( prob%WEIGHT( : prob%n ) )
          prob%X0( : prob%n ) = zero

!  .. or find a centered feasible point ...

        ELSE IF ( center ) THEN
          lsqp = .FALSE.
!         LSQP_control%potential_unbounded = -100.0_wp
          LSQP_control%just_feasible = .FALSE.
          LSQP_control%prfeas = control%QPB_control%prfeas
          LSQP_control%dufeas = control%QPB_control%dufeas
          prob%Hessian_kind = 0
          prob%f = zero

!  .. or minimize the distance to the nearest feasible point

        ELSE
          lsqp = .FALSE.
          LSQP_control%just_feasible = .TRUE.
          LSQP_control%prfeas = control%QPB_control%prfeas
          LSQP_control%dufeas = control%QPB_control%dufeas
          prob%Hessian_kind = 1

          array_name = 'qpc: prob%X0'
          CALL SPACE_resize_array( prob%n, prob%X0, inform%status,             &
                 inform%alloc_status, array_name = array_name,                 &
                 deallocate_error_fatal = control%deallocate_error_fatal,      &
                 exact_size = control%space_critical,                          &
                 bad_alloc = inform%bad_alloc, out = control%error )
          IF ( inform%status /= GALAHAD_ok ) GO TO 900

          array_name = 'qpc: data%X0'
          CALL SPACE_resize_array( prob%n, data%X0, inform%status,             &
                 inform%alloc_status, array_name = array_name,                 &
                 deallocate_error_fatal = control%deallocate_error_fatal,      &
                 exact_size = control%space_critical,                          &
                 bad_alloc = inform%bad_alloc, out = control%error )
          IF ( inform%status /= GALAHAD_ok ) GO TO 900

          data%X0( : prob%n ) = prob%X( : prob%n )
          prob%X0( : prob%n ) = data%X0( : prob%n )
          prob%f = zero
        END IF

!  Solve the LSQP or phase-1 problem

        IF ( lsqp ) THEN
          IF ( printi ) WRITE( control%out, "( /, A, ' entering LSQP ' )" )    &
            prefix
          LSQP_control%indicator_type = control%QPB_control%indicator_type
          LSQP_control%indicator_tol_p = control%QPB_control%indicator_tol_p
          LSQP_control%indicator_tol_pd = control%QPB_control%indicator_tol_pd
          LSQP_control%indicator_tol_tapia =                                   &
            control%QPB_control%indicator_tol_tapia
          CALL LSQP_solve( prob, data, LSQP_control, LSQP_inform,              &
                           C_stat = C_stat( : prob%m ),                        &
                           B_stat = B_stat( : prob%n ) )
!         write(6,*) ' lsqp - status ', LSQP_inform%status
!         write(6,*) ' c_stat ', C_stat( : prob%m )
        ELSE
          IF ( printi )                                                        &
            WRITE( control%out, "( /, A, ' entering LSQP -- phase-one' )" )    &
              prefix
          CALL LSQP_solve( prob, data, LSQP_control, LSQP_inform,              &
                           C_stat = C_stat( : prob%m ),                        &
                           B_stat = B_stat( : prob%n ) )
        END IF

!  Record the exit status
     
        inform%status = LSQP_inform%status 
        inform%QPB_inform%status = LSQP_inform%status 
        inform%QPB_inform%feasible = LSQP_inform%feasible

        IF ( printi ) THEN
          IF ( LSQP_control%out > 0 .AND. LSQP_control%print_level > 0 )       &
            WRITE( control%out, "( /, A, ' returned from LSQP' )" ) prefix
          WRITE( control%out, "( A, ' on exit from LSQP: status = ', I0,       &
         &   ', iterations = ', I0, ', time = ', F0.2 )" ) prefix,             &
            LSQP_inform%status, LSQP_inform%iter, LSQP_inform%time%total
          IF ( lsqp ) WRITE( control%out, "( A, ' objective value =', ES12.4,  &
         &   /, A, ' # active counstraints: ', I0, ' from ', I0,  & 
         &   ', bounds: ', I0, ' from ', I0 )" ) prefix, LSQP_inform%obj,      &
            prefix, COUNT( C_stat( : prob%m ) /= 0 ),  prob%m,                 &
            COUNT( B_stat( : prob%n ) /= 0 ),  prob%n
        END IF

!  If the analytic center appears to be unbounded, have another attempt 
!  at getting feasible

        IF ( inform%QPB_inform%status == GALAHAD_error_upper_entry .OR.        &
             inform%QPB_inform%status == GALAHAD_error_factorization .OR.      &
             inform%QPB_inform%status == GALAHAD_error_ill_conditioned .OR.    &
             inform%QPB_inform%status == GALAHAD_error_tiny_step .OR.          &
             inform%QPB_inform%status == GALAHAD_error_max_iterations .OR.     &
           ( inform%QPB_inform%status == GALAHAD_error_unbounded               &
              .AND. .NOT. lsqp ) )  THEN
          IF ( inform%QPB_inform%feasible ) THEN
            inform%status = GALAHAD_ok
          ELSE
            IF ( first_pass .AND. .NOT. lsqp ) THEN
              center = .NOT. center
              first_pass = .FALSE.
              IF ( printi ) WRITE( control%out,                                &
              "( /, ' .... have a second attempt at getting feasible .... ' )" )
              GO TO 10
            END IF
          END IF
        END IF
        prob%f = f
!       IF ( .NOT. lsqp .AND. .NOT. center ) NULLIFY( prob%X0 )
        prob%Hessian_kind = - 1 ; prob%gradient_kind = - 1

!  Record times for LSQP/phase-1

        CALL CPU_TIME( dum )
        inform%QPB_inform%time%phase1_total = LSQP_inform%time%total
        inform%QPB_inform%time%phase1_analyse = LSQP_inform%time%analyse
        inform%QPB_inform%time%phase1_factorize = LSQP_inform%time%factorize
        inform%QPB_inform%time%phase1_solve = LSQP_inform%time%solve
        inform%QPB_inform%time%total = inform%QPB_inform%time%phase1_total
        inform%QPB_inform%time%analyse = inform%QPB_inform%time%phase1_analyse
        inform%QPB_inform%time%factorize =                                     &
          inform%QPB_inform%time%phase1_factorize
        inform%QPB_inform%time%solve = inform%QPB_inform%time%phase1_solve
   
        IF ( printi )                                                          &
          WRITE( control%out, "( A, ' ', I0, ' integer and ', I0,              &
        &  ' real words required for factors in LSQP' )" ) prefix,             &
           LSQP_inform%factorization_integer, LSQP_inform%factorization_real

!  Record output information from LSQP

        inform%QPB_inform%alloc_status = LSQP_inform%alloc_status 
        inform%QPB_inform%iter = LSQP_inform%iter
        inform%QPB_inform%factorization_status =                               &
          LSQP_inform%factorization_status 
        inform%QPB_inform%factorization_integer =                              &
          LSQP_inform%factorization_integer 
        inform%QPB_inform%factorization_real =                                 &
          LSQP_inform%factorization_real 
        inform%QPB_inform%nfacts = LSQP_inform%nfacts
        inform%QPB_inform%nbacts = LSQP_inform%nbacts
        inform%QPB_inform%nmods = 0
        inform%QPB_inform%obj = LSQP_inform%obj

        inform%alloc_status = inform%QPB_inform%alloc_status
        inform%nfacts = inform%nfacts + inform%QPB_inform%nfacts
        inform%nmods = inform%nmods + inform%QPB_inform%nmods
        inform%obj = inform%QPB_inform%obj
        inform%time%total = inform%time%total + inform%QPB_inform%time%total
        inform%time%analyse =                                                  &
          inform%time%analyse + inform%QPB_inform%time%analyse
        inform%time%factorize =                                                &
          inform%time%factorize + inform%QPB_inform%time%factorize
        inform%time%solve =inform%time%solve + inform%QPB_inform%time%solve
        inform%time%preprocess =                                               &
          inform%time%preprocess + inform%QPB_inform%time%preprocess

!  Check for error exits

        IF ( inform%status /= GALAHAD_ok .AND.                                 &
             inform%status /= GALAHAD_error_primal_infeasible .AND.            & 
             inform%status /= GALAHAD_error_ill_conditioned .AND.              &
             inform%status /= GALAHAD_error_tiny_step .AND.                    &
             inform%status /= GALAHAD_error_max_iterations ) THEN

!  On error exit, compute the current objective function value

          IF ( inform%status /= GALAHAD_error_allocate .AND.                   &
               inform%status /= GALAHAD_error_deallocate ) THEN
            data%HX( : prob%n ) = zero
            CALL QPC_HX( data%dims, prob%n, data%HX( : prob%n ),               &
                         prob%H%ptr( prob%n + 1 ) - 1, prob%H%val,             &
                         prob%H%col, prob%H%ptr, prob%X( : prob%n ), '+' )
            inform%obj = half * DOT_PRODUCT( prob%X( : prob%n ),               &
                                             data%HX( : prob%n ) )             &
                              + DOT_PRODUCT( prob%X( : prob%n ),               &
                                        prob%G( : prob%n ) ) + prob%f
          END IF
          remap_more_freed = .FALSE. ; remap_fixed = .FALSE.
          GO TO 700
        END IF

!  Move to crossover or exit (as appropriate) if the problem is an LSQP

        IF ( lsqp ) THEN
          remap_more_freed = .FALSE. ; remap_fixed = .FALSE.
          IF ( ( control%qpb_or_qpa  .AND. inform%status == GALAHAD_ok ) .OR.  &
                 control%no_qpa ) THEN
!         IF ( control%qpb_or_qpa .OR. control%no_qpa ) THEN
!         IF ( .NOT. control%qpb_or_qpa .AND. control%no_qpa ) THEN
!         IF ( control%no_qpa .AND.                                            &
!             inform%status /= GALAHAD_error_ill_conditioned .AND.             &
!             inform%status /= GALAHAD_error_tiny_step ) THEN
            GO TO 700
          ELSE

!  Save the solution in case QPA fails

            array_name = 'qpc: data%X_trial'
            CALL SPACE_resize_array( prob%n, data%X_trial, inform%status,      &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= GALAHAD_ok ) GO TO 900

            array_name = 'qpc: data%Y_last'
            CALL SPACE_resize_array( prob%m, data%Y_last, inform%status,       &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= GALAHAD_ok ) GO TO 900

            array_name = 'qpc: data%Z_last'
            CALL SPACE_resize_array( prob%n, data%Z_last, inform%status,       &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )

            IF ( inform%status /= GALAHAD_ok ) GO TO 900
            array_name = 'qpc: data%C'
            CALL SPACE_resize_array( prob%m, data%C, inform%status,            &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= GALAHAD_ok ) GO TO 900

            gotsol = .TRUE.
            data%X_trial( : prob%n ) = prob%X( : prob%n )
            data%Y_last( : prob%m ) = prob%Y( : prob%m )
            data%Z_last( : prob%n ) = prob%Z( : prob%n )
            data%C( : prob%m ) = prob%C( : prob%m )
            best_obj = inform%obj

            GO TO 150
          END IF
        END IF



!  ============================
!  Initial feasible point found
!  ============================

!  Check to see if any variables/constraints are flagged as being fixed

        tol = MIN( control%QPB_control%stop_p / ten, SQRT( epsmch ) )
        tiny_x = COUNT( prob%X( data%dims%x_free + 1 : data%dims%x_l_start - 1)&
                        < tol ) +                                              &
                 COUNT( prob%X( data%dims%x_l_start: data%dims%x_l_end ) -     &
                        prob%X_l( data%dims%x_l_start : data%dims%x_l_end )    &
                        < tol ) +                                              &
                 COUNT( prob%X( data%dims%x_u_start: data%dims%x_u_end ) -     &
                        prob%X_u( data%dims%x_u_start : data%dims%x_u_end )    &
                        > - tol ) +                                            &
                 COUNT( prob%X( data%dims%x_u_end + 1 : prob%n )          &
                        > - tol )
  
        tiny_c = COUNT( prob%C( data%dims%c_l_start: data%dims%c_l_end ) -     &
                        prob%C_l( data%dims%c_l_start : data%dims%c_l_end )    &
                        < tol ) +                                              &
                 COUNT( prob%C( data%dims%c_u_start: data%dims%c_u_end ) -     &
                        prob%C_u( data%dims%c_u_start : data%dims%c_u_end )    &
                        > - tol )
  
        remap_fixed = tiny_x > 0 .OR. tiny_c > 0
        IF ( remap_fixed ) THEN

!  Some of the current variables/constraints will be fixed

          IF ( control%error > 0 .AND. control%print_level >= 1 )              &
            WRITE( control%out, "( /, ' -> ', I0, ' further variables and ',   &
           &       I0, ' further constraints will be fixed' )" ) tiny_x, tiny_c

!  Allocate arrays to record the bounds which will be altered

          IF ( tiny_x > 0 ) THEN
            array_name = 'qpc: data%X_fixed'
            CALL SPACE_resize_array( tiny_x, data%X_fixed, inform%status,      &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= GALAHAD_ok ) GO TO 900

            array_name = 'qpc: data%Index_X_fixed'
            CALL SPACE_resize_array( tiny_x,                                   &
                   data%Index_X_fixed, inform%status,                          &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= GALAHAD_ok ) GO TO 900
          END IF
  
          IF ( tiny_c > 0 ) THEN
            array_name = 'qpc: data%C_fixed'
            CALL SPACE_resize_array( tiny_c, data%C_fixed, inform%status,      &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= GALAHAD_ok ) GO TO 900

            array_name = 'qpc: data%Index_C_fixed'
            CALL SPACE_resize_array( tiny_c,                                   &
                   data%Index_C_fixed, inform%status,                          &
                   inform%alloc_status, array_name = array_name,               &
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= GALAHAD_ok ) GO TO 900
          END IF

!  Fix the problem bounds as required

          IF ( tiny_x > 0 ) THEN
            tiny_x = 0
    
            DO i = data%dims%x_free + 1, data%dims%x_l_start - 1
!             write(6,"( I6, A1, ES12.4 )" ) i, 'l', prob%X( i )
              IF ( prob%X( i ) < tol ) THEN
                tiny_x = tiny_x + 1
                data%X_fixed( tiny_x ) = prob%X_u( i )
                data%Index_X_fixed( tiny_x ) = - i
                prob%X_u( i ) = zero
              END IF
            END DO
    
            DO i = data%dims%x_l_start, data%dims%x_u_start - 1
!             write(6,"( I6, A1, ES12.4 )" ) i, 'l', prob%X( i ) - prob%X_l( i )
              IF ( prob%X( i ) - prob%X_l( i ) < tol ) THEN
                tiny_x = tiny_x + 1
                data%X_fixed( tiny_x ) = prob%X_u( i )
                data%Index_X_fixed( tiny_x ) = - i
                prob%X_u( i ) =  prob%X_l( i )
              END IF
            END DO
    
            DO i = data%dims%x_u_start, data%dims%x_l_end
!             write(6,"( I6, A1, ES12.4 )" ) i, 'l', prob%X( i ) - prob%X_l( i )
!             write(6,"( I6, A1, ES12.4 )" ) i, 'u', prob%X_u( i ) - prob%X( i )
              IF ( prob%X( i ) - prob%X_l( i ) < tol ) THEN
                tiny_x = tiny_x + 1
                data%X_fixed( tiny_x ) = prob%X_u( i )
                data%Index_X_fixed( tiny_x ) = - i
                prob%X_u( i ) =  prob%X_l( i )
              ELSE IF ( prob%X( i ) - prob%X_u( i ) > - tol ) THEN
                tiny_x = tiny_x + 1
                data%X_fixed( tiny_x ) = prob%X_l( i )
                data%Index_X_fixed( tiny_x ) = i
                prob%X_l( i ) =  prob%X_u( i )
              END IF
            END DO
  
            DO i = data%dims%x_l_end + 1, data%dims%x_u_end
!             write(6,"( I6, A1, ES12.4 )" ) i, 'u', prob%X_u( i ) - prob%X( i )
              IF ( prob%X( i ) - prob%X_u( i ) > - tol ) THEN
                tiny_x = tiny_x + 1
                data%X_fixed( tiny_x ) = prob%X_l( i )
                data%Index_X_fixed( tiny_x ) = i
                prob%X_l( i ) =  prob%X_u( i )
              END IF
            END DO
    
            DO i = data%dims%x_u_end + 1, prob%n
!             write(6,"( I6, A1, ES12.4 )" ) i, 'u', - prob%X( i )
              IF ( prob%X( i ) > - tol ) THEN
                tiny_x = tiny_x + 1
                data%X_fixed( tiny_x ) = prob%X_l( i )
                data%Index_X_fixed( tiny_x ) = i
                prob%X_l( i ) = zero
              END IF
            END DO
          END IF
       
!  Do the same for the constraint bounds

          IF ( tiny_c > 0 ) THEN
            tiny_c = 0
    
            DO i = data%dims%c_l_start, data%dims%c_u_start - 1
              IF ( prob%C( i ) - prob%C_l( i ) < tol ) THEN
                tiny_c = tiny_c + 1
                data%C_fixed( tiny_c ) = prob%C_u( i )
                data%Index_C_fixed( tiny_c ) = - i
                prob%C_u( i ) =  prob%C_l( i )
              END IF
            END DO
    
            DO i = data%dims%c_u_start, data%dims%c_l_end
              IF ( prob%C( i ) - prob%C_l( i ) < tol ) THEN
                tiny_c = tiny_c + 1
                data%C_fixed( tiny_c ) = prob%C_u( i )
                data%Index_C_fixed( tiny_c ) = - i
                prob%C_u( i ) =  prob%C_l( i )
              ELSE IF ( prob%C( i ) - prob%C_u( i ) > - tol ) THEN
                tiny_c = tiny_c + 1
                data%C_fixed( tiny_c ) = prob%C_l( i )
                data%Index_C_fixed( tiny_c ) = i
                prob%C_l( i ) =  prob%C_u( i )
              END IF
            END DO
    
            DO i = data%dims%c_l_end + 1, data%dims%c_u_end
              IF ( prob%C( i ) - prob%C_u( i ) > - tol ) THEN
                tiny_c = tiny_c + 1
                data%C_fixed( tiny_c ) = prob%C_l( i )
                data%Index_C_fixed( tiny_c ) = i
                prob%C_l( i ) =  prob%C_u( i )
              END IF
            END DO
          END IF

! write(6,"( ' tiny_x = ', I0, /, ( 10I5 ))" ) tiny_x, data%Index_X_fixed
! write(6,"( ' tiny_c = ', I0, /, ( 10I5 ))" ) tiny_c, data%Index_C_fixed
          CALL QPP_initialize( data%QPP_map_fixed, data%QPP_control )
          data%QPP_control%infinity = control%infinity
          data%QPP_control%treat_zero_bounds_as_general =                     &
            control%QPB_control%treat_zero_bounds_as_general

!  Store the problem dimensions

          data%dims_save_fixed = data%dims
          a_ne = prob%A%ne 
          h_ne = prob%H%ne 
  
          IF ( printi ) WRITE( control%out,                                    &
                 "( /, A, ' problem dimensions before preprocessing: ', /,  A, &
     &           ' n = ', I8, ' m = ', I8, ' a_ne = ', I8, ' h_ne = ', I8 )" ) &
                 prefix, prefix, prob%n, prob%m, a_ne, h_ne

!  Perform the preprocessing

          CALL CPU_TIME( time ) 
          CALL QPP_reorder( data%QPP_map_fixed, data%QPP_control,              &
                             data%QPP_inform, data%dims, prob,                 &
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
            IF ( control%out > 0 .AND. control%print_level >= 1 )              &
              WRITE( control%out, "( A, ' status ', I3, ' after QPP_reorder')")&
               prefix, data%QPP_inform%status
            IF ( control%error > 0 .AND. control%print_level > 0 )             &
              WRITE( control%error, 2010 ) inform%status 
            IF ( control%out > 0 .AND. control%print_level > 0 .AND.           &
                 inform%status == GALAHAD_error_upper_entry )                  &
              WRITE( control%error, 2240 ) 
            CALL QPP_terminate( data%QPP_map_fixed, data%QPP_control,          &
                                data%QPP_inform )
            CALL QPP_terminate( data%QPP_map_freed, data%QPP_control,          &
                                data%QPP_inform )
            CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )
            CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
            IF ( data%QPP_inform%status == GALAHAD_error_primal_infeasible ) THEN
              inform%status = GALAHAD_error_primal_infeasible
              GO TO 800 
            ELSE
              GO TO 800 
            END IF 
          END IF 

!  Record revised array lengths

          IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
            a_ne = prob%m * prob%n
          ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
            a_ne = prob%A%ptr( prob%m + 1 ) - 1
          ELSE
            a_ne = prob%A%ne 
          END IF

          IF ( SMT_get( prob%H%type ) == 'DIAGONAL' ) THEN
            h_ne = prob%n
          ELSE IF ( SMT_get( prob%H%type ) == 'DENSE' ) THEN
            h_ne = ( prob%n * ( prob%n + 1 ) ) / 2
          ELSE IF ( SMT_get( prob%H%type ) == 'SPARSE_BY_ROWS' ) THEN
            h_ne = prob%H%ptr( prob%n + 1 ) - 1
          ELSE
            h_ne = prob%H%ne 
          END IF

          IF ( printi ) WRITE( control%out,                                    &
                 "( A, ' problem dimensions after preprocessing: ', /,  A,     &
     &           ' n = ', I8, ' m = ', I8, ' a_ne = ', I8, ' h_ne = ', I8 )" ) &
                 prefix, prefix, prob%n, prob%m, a_ne, h_ne

!  Further permute initial working sets if provided

          IF ( control%QPA_control%cold_start == 0 ) THEN
            CALL SORT_inplace_permute( data%QPP_map_fixed%m,                   &
                                       data%QPP_map_fixed%c_map,               &
                                       IX = C_stat( : data%QPP_map_fixed%m ) )
            CALL SORT_inplace_permute( data%QPP_map_fixed%n,                   &
                                       data%QPP_map_fixed%x_map,               &
                                       IX = B_stat( : data%QPP_map_fixed%n ) )
          END IF

!  If all the variables have now been fixed, the solution has been found

          IF ( prob%n == 0 ) THEN

!  Check that the solution is feasible

            DO i = 1, prob%m
              IF ( prob%C_l( i ) > control%QPB_control%stop_p .OR.             &
                   prob%C_u( i ) < - control%QPB_control%stop_p ) THEN
                inform%status = GALAHAD_error_primal_infeasible
                GO TO 800 
              END IF
            END DO
            prob%C( : prob%m ) = zero
            prob%Y( : prob%m ) = zero
            inform%status = 0
            GO TO 720
          END IF

!  ====================================================================
!  Check to see if the equality constraints remain linearly independent
!  ====================================================================

          IF ( prob%m > 0 .AND. control%QPB_control%remove_dependencies ) THEN
    
            CALL CPU_TIME( time ) 
            IF ( control%out > 0 .AND. control%print_level >= 1 )              &
              WRITE( control%out,                                              &
                "( /, A, 1X, I0, ' equalities from ', I0, ' constraints' )" )  &
                prefix, data%dims%c_equality, prob%m

!  Set control parameters

            CALL FDC_initialize( FDC_control )
            FDC_control%error = control%error
            FDC_control%out = control%out
            FDC_control%print_level = control%print_level
            FDC_control%indmin = control%indmin
            FDC_control%valmin = control%valmin
            FDC_control%zero_pivot = control%zero_pivot
            FDC_control%pivot_tol = control%pivot_tol_for_dependencies
            FDC_control%max_infeas = control%QPB_control%stop_p
            FDC_control%CNTL = data%CNTL
            FDC_control%prefix = '" - FDC:"                     '

!  Find any dependent rows

            nzc = prob%A%ptr( data%dims%c_equality + 1 ) - 1 
            CALL FDC_find_dependent( prob%n, data%dims%c_equality,             &
                                     prob%A%val( : nzc ),                      &
                                     prob%A%col( : nzc ),                      &
                                     prob%A%ptr( : data%dims%c_equality + 1 ), &
                                     prob%C_l, data%K, data%FACTORS,           &
                                     n_more_depen, data%Index_C_more_freed,    &
                                     FDC_control, FDC_inform )
            inform%non_negligible_pivot =                                      &
              MIN( FDC_inform%non_negligible_pivot, inform%non_negligible_pivot )

!  Record output parameters

            inform%alloc_status = FDC_inform%alloc_status
            inform%factorization_status = FDC_inform%factorization_status
            inform%factorization_integer = FDC_inform%factorization_integer
            inform%factorization_real = FDC_inform%factorization_real
            inform%bad_alloc = FDC_inform%bad_alloc
            inform%time%find_dependent =                                       &
              inform%time%find_dependent + FDC_inform%time%total
            inform%time%analyse = inform%time%analyse + FDC_inform%time%analyse
            inform%time%factorize =                                           &
              inform%time%factorize + FDC_inform%time%factorize
            inform%nfacts = 1

            CALL CPU_TIME( dum )
            IF ( control%cpu_time_limit >= zero .AND.                          &
                 dum - time_start > control%cpu_time_limit ) THEN
              inform%status = GALAHAD_error_cpu_limit
              IF ( control%error > 0 .AND. control%print_level > 0 )           &
                WRITE( control%error, 2010 ) inform%status 
              inform%time%total = dum - time_start ; GO TO 800 
            END IF 
            inform%time%preprocess = inform%time%preprocess + dum - time
            IF ( printi .AND. inform%non_negligible_pivot <                    &
                 thousand * control%zero_pivot ) WRITE( control%out, "(        &
           &  /, 1X, 26 ( '*' ), ' WARNING ', 26 ( '*' ), /,                   & 
           &  ' ***  smallest allowed pivot =', ES11.4,' may be too small ***',&
           &  /, ' ***  perhaps increase control%zero_pivot from', ES11.4,     &
           &  '  ***', /, 1X, 26 ( '*' ), ' WARNING ', 26 ( '*' ) )" )         &
               inform%non_negligible_pivot, control%zero_pivot

!  Check for error exits

            IF ( FDC_inform%status /= GALAHAD_ok ) THEN

!  Allocate arrays to hold the matrix vector product

              array_name = 'qpc: data%HX'
              CALL SPACE_resize_array( prob%n, data%HX, inform%status,         &
                     inform%alloc_status, array_name = array_name,             &
                     deallocate_error_fatal = control%deallocate_error_fatal,  &
                     exact_size = control%space_critical,                      &
                     bad_alloc = inform%bad_alloc, out = control%error )
              IF ( inform%status /= GALAHAD_ok ) GO TO 900
        
!  On error exit, compute the current objective function value

              data%HX( : prob%n ) = zero
              CALL QPC_HX( data%dims, prob%n, data%HX( : prob%n ),             &
                            prob%H%ptr( prob%n + 1 ) - 1, prob%H%val,          &
                            prob%H%col, prob%H%ptr, prob%X( : prob%n ), '+' )
              inform%obj = half * DOT_PRODUCT( prob%X( : prob%n ),             &
                                               data%HX( : prob%n ) )           &
                           + DOT_PRODUCT( prob%X( : prob%n ),                  &
                                          prob%G( : prob%n ) ) + prob%f

!  Print details of the error exit

              IF ( control%error > 0 .AND. control%print_level >= 1 ) THEN
                WRITE( control%out, "( ' ' )" )
                WRITE( control%error, 2040 ) FDC_inform%status,                &
                  'FDC_find_dependent'
              END IF
              inform%status = FDC_inform%status
              GO TO 700
            END IF
    
            IF ( control%out > 0 .AND. control%print_level >= 2                &
                 .AND. n_more_depen > 0 )                                      &
              WRITE( control%out, "(/, ' The following ', I0,                  &
           &         ' constraints appear to be dependent', /, ( 8I8 ) )" )    &
                n_more_depen, data%Index_C_more_freed
    
            remap_more_freed = n_more_depen > 0
          ELSE
            remap_more_freed = .FALSE.
          END IF
    
          IF ( remap_more_freed ) THEN

!  Some of the current constraints will be removed by freeing them

            IF ( control%error > 0 .AND. control%print_level >= 1 )            &
              WRITE( control%out, "( /, ' -> ', I0, ' constraints are',        &
             & ' dependent and will be temporarily removed' )" ) n_more_depen

!  Allocate arrays to indicate which constraints have been freed

            array_name = 'qpc: data%C_more_freed'
            CALL SPACE_resize_array( n_more_depen, data%C_more_freed,          &
                   inform%status, inform%alloc_status, array_name = array_name,&
                   deallocate_error_fatal = control%deallocate_error_fatal,    &
                   exact_size = control%space_critical,                        &
                   bad_alloc = inform%bad_alloc, out = control%error )
            IF ( inform%status /= GALAHAD_ok ) GO TO 900

!  Free the constraint bounds as required

            DO i = 1, n_more_depen
              j = data%Index_C_more_freed( i )
              data%C_more_freed( i ) = prob%C_l( j )
              prob%C_l( j ) = - control%infinity
              prob%C_u( j ) = control%infinity
              prob%Y( j ) = zero
            END DO
    
            CALL QPP_initialize( data%QPP_map_more_freed, data%QPP_control )
            data%QPP_control%infinity = control%infinity
            data%QPP_control%treat_zero_bounds_as_general =                    &
              control%QPB_control%treat_zero_bounds_as_general

!  Store the problem dimensions

            data%dims_save_more_freed = data%dims
            a_ne = prob%A%ne 
            h_ne = prob%H%ne 
    
            IF ( printi ) WRITE( control%out,                                  &
               "( /, A, ' problem dimensions before removal of dependencies:', &
     &            /, A,                                                        &
     &         ' n = ', I8, ' m = ', I8, ' a_ne = ', I8, ' h_ne = ', I8 )" )   &
               prefix, prefix, prob%n, prob%m, a_ne, h_ne

!  Perform the preprocessing

            CALL CPU_TIME( time ) 
            CALL QPP_reorder( data%QPP_map_more_freed, data%QPP_control,       &
                              data%QPP_inform, data%dims, prob,                &
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
              IF ( control%out > 0 .AND. control%print_level >= 1 )            &
                WRITE( control%out,"( A, ' status ', I3,' after QPP_reorder')")&
                 prefix, data%QPP_inform%status
              IF ( control%error > 0 .AND. control%print_level > 0 )           &
                WRITE( control%error, 2010 ) inform%status 
              IF ( control%out > 0 .AND. control%print_level > 0 .AND.         &
                   inform%status == GALAHAD_error_upper_entry )                &
                WRITE( control%error, 2240 ) 
              CALL QPP_terminate( data%QPP_map_more_freed, data%QPP_control,   &
                                  data%QPP_inform )
              CALL QPP_terminate( data%QPP_map_fixed, data%QPP_control,        &
                                  data%QPP_inform )
              CALL QPP_terminate( data%QPP_map_freed, data%QPP_control,        &
                                  data%QPP_inform )
              CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
              IF ( data%QPP_inform%status == GALAHAD_error_primal_infeasible ) &
                THEN
                inform%status = GALAHAD_error_primal_infeasible
                GO TO 800 
              ELSE
                GO TO 800 
              END IF 
            END IF 

!  Record revised array lengths

            IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN
              a_ne = prob%m * prob%n
            ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN
              a_ne = prob%A%ptr( prob%m + 1 ) - 1
            ELSE
              a_ne = prob%A%ne 
            END IF

            IF ( SMT_get( prob%H%type ) == 'DIAGONAL' ) THEN
              h_ne = prob%n
            ELSE IF ( SMT_get( prob%H%type ) == 'DENSE' ) THEN
              h_ne = ( prob%n * ( prob%n + 1 ) ) / 2
            ELSE IF ( SMT_get( prob%H%type ) == 'SPARSE_BY_ROWS' ) THEN
              h_ne = prob%H%ptr( prob%n + 1 ) - 1
            ELSE
              h_ne = prob%H%ne 
            END IF
    
            IF ( printi ) WRITE( control%out,                                  &
               "( A, ' problem dimensions after removal of dependencies:',     &
     &            /, A,                                                        &
     &         ' n = ', I8, ' m = ', I8, ' a_ne = ', I8, ' h_ne = ', I8 )" )   &
               prefix, prefix, prob%n, prob%m, a_ne, h_ne

!  Further permute initial working sets if provided

            IF ( control%QPA_control%cold_start == 0 ) THEN
              CALL SORT_inplace_permute( data%QPP_map_more_freed%m,            &
                                         data%QPP_map_more_freed%c_map,        &
                                    IX = C_stat( : data%QPP_map_more_freed%m ) )
              CALL SORT_inplace_permute( data%QPP_map_more_freed%n,            &
                                         data%QPP_map_more_freed%x_map,        &
                                    IX = B_stat( : data%QPP_map_more_freed%n ) )
            END IF
          END IF

!  Experiment!!

!         GO TO 10
        ELSE
          remap_more_freed = .FALSE.
        END IF
      END IF

!  Allocate additional real workspace

      array_name = 'qpc: data%DZ_l'
      CALL SPACE_resize_array( data%dims%x_l_start, data%dims%x_l_end,         &
             data%DZ_l, inform%status,                                         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%DZ_u'
      CALL SPACE_resize_array( data%dims%x_u_start, data%dims%x_u_end,         &
             data%DZ_u, inform%status,                                         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%GRAD'
      CALL SPACE_resize_array( prob%n, data%GRAD, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%X_trial'
      CALL SPACE_resize_array( prob%n, data%X_trial, inform%status,            &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%GRAD_X_phi'
      CALL SPACE_resize_array( prob%n, data%GRAD_X_phi, inform%status,         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%GRAD_C_phi'
      CALL SPACE_resize_array( data%dims%c_l_start, data%dims%c_u_end,         &
             data%GRAD_C_phi, inform%status,                                   &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%S'
      CALL SPACE_resize_array( data%dims%c_e, data%S, inform%status,           &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%A_s'
      CALL SPACE_resize_array( prob%m, data%A_s, inform%status,                &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%H_s'
      CALL SPACE_resize_array( prob%n, data%H_s, inform%status,                &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%Y_last'
      CALL SPACE_resize_array( prob%m, data%Y_last, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%Z_last'
      CALL SPACE_resize_array( prob%n, data%Z_last, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

!  ===============================
!  Approximately solve the problem
!  ===============================

!  overlaps: DY_l => DIST_C_l_trial 
!            DY_u => DIST_C_u_trial
!            DELTA => VECTOR
!            RHS( : c_e ) => R 
!            DZ_l( x_l_start : x_l_end ) => DIST_X_l_trial
!            DZ_u( x_u_start : x_u_end ) => DIST_X_u_trial

      control%QPB_control%LSQP_control%indicator_type =                        &
        control%QPB_control%indicator_type
      control%QPB_control%LSQP_control%indicator_tol_p =                       &
        control%QPB_control%indicator_tol_p
      control%QPB_control%LSQP_control%indicator_tol_pd =                      &
        control%QPB_control%indicator_tol_pd
      control%QPB_control%LSQP_control%indicator_tol_tapia =                   &
        control%QPB_control%indicator_tol_tapia

      IF ( printi ) WRITE( control%out, "( /, A, ' entering QPB ' )" ) prefix

      CALL QPB_solve_main( data%dims, prob%n, prob%m,                          &
                           prob%H%val, prob%H%col, prob%H%ptr,                 &
                           prob%G, prob%f, prob%A%val, prob%A%col,             & 
                           prob%A%ptr, prob%C_l, prob%C_u, prob%X_l,           &
                           prob%X_u, prob%C, prob%X, prob%Y, prob%Z,           &
                           data%RES_x, data%X_trial, data%SOL_y, data%RES_y,   &
                           data%BEST_y, data%SOL, data%RES, data%BEST,         &
                           data%HX, data%GRAD_L, data%DIST_X_l,                &
                           data%DIST_X_u, data%Z_l, data%Z_u, data%BARRIER_X,  &
                           data%Y_l, data%DY_l, data%DIST_C_l, data%Y_u,       &
                           data%DY_u, data%DIST_C_u, data%C, data%BARRIER_C,   &
                           data%SCALE_C, data%DELTA,                           &
                           data%RHS( : data%dims%v_e ),                        &
                           data%GRAD, data%GRAD_X_phi, data%GRAD_C_phi,        &
                           data%DZ_l( data%dims%x_l_start :                    &
                                      data%dims%x_l_end ),                     &
                           data%DZ_u( data%dims%x_u_start :                    &
                                      data%dims%x_u_end ), data%S,             &
                           data%Abycol_val, data%DIAG_X, data%DIAG_C,          &
                           data%IW, data%K_colptr, data%Abycol_ptr,            &
                           data%Abycol_row, data%H_band_ptr,                   &
                           data%K, data%FACTORS, data%CNTL,                    &
                           control%QPB_control, inform%QPB_inform,             &
                           C_last = data%A_s, X_last = data%H_s,               &
                           Y_last = data%Y_last, Z_last = data%Z_last,         &
                           C_stat = C_stat( : prob%m ),                        &
                           B_stat = B_stat( : prob%n ) )
      IF ( printi ) THEN
        IF ( control%QPB_control%out > 0 .AND.                                 &
             control%QPB_control%print_level > 0 )                             &
          WRITE( control%out, "( /, A, ' returned from QPB' )" ) prefix
        WRITE( control%out, "( A, ' on exit from QPB: status = ', I0,          &
       &   ', iterations = ', I0, ', time = ', F0.2 )" ) prefix,               &
          inform%QPB_inform%status, inform%QPB_inform%iter,                    &
          inform%QPB_inform%time%total
        WRITE( control%out, "( A, ' objective value =', ES12.4,                &
       &   /, A, ' # active counstraints: ', I0, ' from ', I0,                 & 
       &   ', bounds: ', I0, ' from ', I0 )" ) prefix, inform%QPB_inform%obj,  &
            prefix, COUNT( C_stat( : prob%m ) /= 0 ),  prob%m,                 &
            COUNT( B_stat( : prob%n ) /= 0 ),  prob%n
      END IF

!  Record output information from QPB

      inform%status = inform%QPB_inform%status
      inform%alloc_status = inform%QPB_inform%alloc_status
      inform%nfacts = inform%nfacts + inform%QPB_inform%nfacts
      inform%nmods = inform%nmods + inform%QPB_inform%nmods
      inform%obj = inform%QPB_inform%obj

!  Record times for components of QPB

      inform%time%total = inform%time%total + inform%QPB_inform%time%total
      inform%time%analyse =                                                    &
        inform%time%analyse + inform%QPB_inform%time%analyse
      inform%time%factorize =                                                  &
        inform%time%factorize + inform%QPB_inform%time%factorize
      inform%time%solve =inform%time%solve + inform%QPB_inform%time%solve
      inform%time%preprocess =                                                 &
        inform%time%preprocess + inform%QPB_inform%time%preprocess

      IF ( printi )                                                            &
        WRITE( control%out, "( A, ' ', I0, ' integer and ', I0,                &
      &  ' real words required for factors in QPB' )" ) prefix,                &
         inform%QPB_inform%factorization_integer,                              &
         inform%QPB_inform%factorization_real

!  Check for error exits

      IF ( inform%status /= GALAHAD_ok .AND.                                   &
           inform%status /= GALAHAD_error_ill_conditioned .AND.                &
           inform%status /= GALAHAD_error_tiny_step .AND.                      &
           inform%status /= GALAHAD_error_max_iterations ) THEN
        GO TO 700
      END IF

!  Check to see if crossover is required

      IF ( control%qpb_or_qpa .AND. inform%status == GALAHAD_ok ) GO TO 700
      IF ( .NOT. control%qpb_or_qpa .AND. control%no_qpa ) GO TO 700

!     IF ( control%no_qpa .AND.                                                &
!          inform%status /= GALAHAD_error_ill_conditioned .AND.                &
!          inform%status /= GALAHAD_error_tiny_step ) GO TO 700

!  Save the solution in case QPA fails

      array_name = 'qpc: data%X_trial'
      CALL SPACE_resize_array( prob%n, data%X_trial, inform%status,            &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%Y_last'
      CALL SPACE_resize_array( prob%m, data%Y_last, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%Z_last'
      CALL SPACE_resize_array( prob%n, data%Z_last, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )

      IF ( inform%status /= GALAHAD_ok ) GO TO 900
      array_name = 'qpc: data%C'
      CALL SPACE_resize_array( prob%m, data%C, inform%status,                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      gotsol = .TRUE.
      data%X_trial( : prob%n ) = prob%X( : prob%n )
      data%Y_last( : prob%m ) = prob%Y( : prob%m )
      data%Z_last( : prob%n ) = prob%Z( : prob%n )
      data%C( : prob%m ) = prob%C( : prob%m )
      best_obj = inform%obj

  150 CONTINUE

!     OPEN( 56 )
!     WRITE( 56, "( '    i st    x - x_l       x_u - x        z ' )" )
!     DO i = 1, prob%n
!       IF ( B_stat( i ) /= 0 ) THEN      
!         WRITE( 56, "( i6, i3, 3ES12.4 )" ) i, B_stat( i ),                    &
!           prob%X( i ) - prob%X_l( i ), prob%X_u( i ) - prob%X( i ), prob%Z( i )
!       END IF
!     END DO

!     WRITE( 56, "( '    i st    c - c_l       c_u - c        y ' )" )
!     DO i = 1, prob%m
!       IF ( C_stat( i ) /= 0 ) THEN      
!         WRITE( 56, "( i6, i3, 3ES12.4 )" ) i, C_stat( i ),                    &
!           prob%C( i ) - prob%C_l( i ), prob%C_u( i ) - prob%C( i ), prob%Y( i )
!       END IF
!     END DO
!     CLOSE( 56 )

!  Store the current constraint statii for later checks

      array_name = 'qpc: data%IW'
      CALL SPACE_resize_array( prob%m + prob%n, data%IW, inform%status,        &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      data%IW( : prob%m ) = C_stat( : prob%m )
      data%IW( prob%m + 1 : prob%m + prob%n ) = B_stat( : prob%n )

!  ===========================================================================
!  Check to see if the constraints in the working set are linearly independent
!  ===========================================================================

!  Allocate workspace arrays

      lbreak = prob%m + data%dims%c_l_end - data%dims%c_u_start +              &
               prob%n - data%dims%x_free + data%dims%x_l_end -                 &
               data%dims%x_u_start + 2

      array_name = 'qpc: data%IBREAK'
      CALL SPACE_resize_array( lbreak, data%IBREAK, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%RES_l'
      CALL SPACE_resize_array( 1, data%dims%c_l_end,                           &
             data%RES_l, inform%status,                                        &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%RES_u'
      CALL SPACE_resize_array( data%dims%c_u_start, prob%m,                    &
             data%RES_u, inform%status,                                        &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%A_norms'
      CALL SPACE_resize_array( prob%m, data%A_norms, inform%status,            &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%H_s'
      CALL SPACE_resize_array( prob%n, data%H_s, inform%status,                &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

!  Compute the initial residuals

      DO i = 1, data%dims%c_u_start - 1
        a_norms = zero ; a_x = zero
        DO ii = prob%A%ptr( i ), prob%A%ptr( i + 1 ) - 1
          a_x = a_x + prob%A%val( ii ) * prob%X( prob%A%col( ii ) )
          a_norms = a_norms + prob%A%val( ii ) ** 2
        END DO
!       write(6,*) 'l', i, a_x, prob%C_l( i )
        data%RES_l( i ) = a_x - prob%C_l( i )
        data%A_norms( i ) = SQRT( a_norms )
      END DO

      DO i = data%dims%c_u_start, data%dims%c_l_end
        a_norms = zero ; a_x = zero
        DO ii = prob%A%ptr( i ), prob%A%ptr( i + 1 ) - 1
          a_x = a_x + prob%A%val( ii ) * prob%X( prob%A%col( ii ) )
          a_norms = a_norms + prob%A%val( ii ) ** 2
        END DO 
!       write(6,*) 'l', i, a_x, prob%C_l( i )
!       write(6,*) 'u', i, a_x, prob%C_u( i )
        data%RES_l( i ) = a_x - prob%C_l( i )
        data%RES_u( i ) = prob%C_u( i ) - a_x
        data%A_norms( i ) = SQRT( a_norms )
      END DO

      DO i = data%dims%c_l_end + 1, prob%m
        a_norms = zero ; a_x = zero
        DO ii = prob%A%ptr( i ), prob%A%ptr( i + 1 ) - 1
          a_x = a_x + prob%A%val( ii ) * prob%X( prob%A%col( ii ) )
          a_norms = a_norms + prob%A%val( ii ) ** 2
        END DO
!       write(6,*) 'u', i, a_x, prob%C_u( i )
        data%RES_u( i ) = prob%C_u( i ) - a_x
        data%A_norms( i ) = SQRT( a_norms )
      END DO

!  If necessary, determine which constraints occur in the reference set

      IF ( control%qpb_or_qpa .OR. control%no_qpb ) THEN

!  cold start from the value set in X; constraints active
!  at X will determine the initial working set.

        IF ( control%QPA_control%cold_start == 1 ) THEN

!  constraints with lower bounds

          DO i = 1, data%dims%c_u_start - 1
            IF ( ABS( data%RES_l( i ) ) <= teneps ) THEN
!             write(6,*) i, data%RES_l( i )
              C_stat( i ) = - 1
            ELSE
              C_stat( i ) = 0
            END IF
          END DO

!  constraints with both lower and upper bounds

          DO i = data%dims%c_u_start, data%dims%c_l_end
            IF ( ABS( data%RES_l( i ) ) <= teneps ) THEN
!             write(6,*) i, data%RES_l( i )
              C_stat( i ) = - 1
            ELSE IF ( ABS( data%RES_u( i ) ) <= teneps ) THEN
!             write(6,*) i, data%RES_u( i )
              C_stat( i ) = 1
            ELSE
              C_stat( i ) = 0
            END IF
          END DO

!  constraints with upper bounds

          DO i = data%dims%c_l_end + 1, prob%m
            IF ( ABS( data%RES_u( i ) ) <= teneps ) THEN
!             write(6,*) i, data%RES_u( i )
              C_stat( i ) = 1
            ELSE
              C_stat( i ) = 0
            END IF
          END DO

!  free variables

          B_stat( : data%dims%x_free ) = 0

!  simple non-negativity

          DO i = data%dims%x_free + 1, data%dims%x_l_start - 1
            IF ( ABS( prob%X( i ) ) <= teneps ) THEN
!             write(6,*) prob%m + i, prob%X( i )
              B_stat( i ) = - 1
            ELSE
              B_stat( i ) = 0
            END IF
          END DO

!  simple bound from below

          DO i = data%dims%x_l_start, data%dims%x_u_start - 1
            IF ( ABS( prob%X( i ) - prob%X_l( i ) ) <= teneps ) THEN
!             write(6,*) prob%m + i, prob%X( i ) - prob%X_l( i )
              B_stat( i ) = - 1
            ELSE
              B_stat( i ) = 0
            END IF
          END DO

!  simple bound from below and above

          DO i = data%dims%x_u_start, data%dims%x_l_end
            IF ( ABS( prob%X( i ) - prob%X_l( i ) ) <= teneps ) THEN
!             write(6,*) prob%m + i, prob%X( i ) - prob%X_l( i )
              B_stat( i ) = - 1
            ELSE IF ( ABS( prob%X( i ) - prob%X_u( i ) ) <= teneps ) THEN
!             write(6,*) prob%m + i, prob%X( i ) - prob%X_u( i )
              B_stat( i ) = 1
            ELSE
              B_stat( i ) = 0
            END IF
          END DO

!  simple bound from above

          DO i = data%dims%x_l_end + 1, data%dims%x_u_end
            IF ( ABS( prob%X( i ) - prob%X_u( i ) ) <= teneps ) THEN
!             write(6,*) prob%m + i, prob%X( i ) - prob%X_u( i )
              B_stat( i ) = 1
            ELSE
              B_stat( i ) = 0
            END IF
          END DO

!  simple non-positivity

          DO i = data%dims%x_u_end + 1, prob%n
            IF ( ABS( prob%X( i ) ) <= teneps ) THEN
!             write(6,*) prob%m + i, prob%X( i )
              B_stat( i ) = 1
            ELSE
              B_stat( i ) = 0
            END IF
          END DO

!  cold start with only equality constraints active

        ELSE IF ( control%QPA_control%cold_start == 3 ) THEN
          B_stat = 0 
          C_stat( : MIN( data%dims%c_equality, prob%n ) ) = 1
          C_stat( MIN( data%dims%c_equality, prob%n ) + 1 : ) = 0

!  cold start with as many active constraints as possible

        ELSE IF ( control%QPA_control%cold_start == 4 ) THEN
          B_stat = 0 ; C_stat = 0
          l = 0 

!  equality constraints

          DO i = 1,  data%dims%c_equality
            IF ( l > prob%n ) EXIT
            C_stat( i ) = - 1
            l = l + 1
          END DO

!  simple bound from below

          DO i = data%dims%x_free + 1, data%dims%x_l_end
            IF ( l > prob%n ) EXIT
            B_stat( i ) = - 1
            l = l + 1
          END DO

!  simple bound from above

          DO i = data%dims%x_l_end + 1, data%dims%x_u_end
            IF ( l > prob%n ) EXIT
            B_stat( i ) = 1
            l = l + 1
          END DO

!  constraints with lower bounds

          DO i = data%dims%c_equality + 1, data%dims%c_l_end
            IF ( l > prob%n ) EXIT
            C_stat( i ) = - 1
            l = l + 1
          END DO

!  constraints with upper bounds

          DO i = data%dims%c_l_end + 1, prob%m
            IF ( l > prob%n ) EXIT
            C_stat( i ) = 1
            l = l + 1
          END DO

!  cold start with no active constraints

        ELSE IF ( control%QPA_control%cold_start /= 0 ) THEN
          B_stat = 0 ; C_stat = 0
        END IF
!       WRITE( control%out, "(' b_stat ', /, ( 10I5 ) )" ) B_stat( : prob%n )
!       WRITE( control%out, "(' c_stat ', /, ( 10I5 ) )" ) C_stat( : prob%m )
!       WRITE( out, "( ' b_stat ', /, ( 10I5 ) )" )                            &
!         B_stat( dims%x_free + 1 : prob%n )

!  select the penalty parameters

        IF( control%rho_g <= zero ) THEN
          prob%rho_g = 2 * prob%m 
        ELSE
          prob%rho_g = control%rho_g
        END IF

        IF( control%rho_b <= zero ) THEN
          prob%rho_b = 2 * prob%n
        ELSE
          prob%rho_b = control%rho_b
        END IF

      ELSE

!  select the penalty parameters

        IF( control%rho_g <= zero ) THEN
          prob%rho_g = prob%m
          prob%rho_g = two * MAX( prob%rho_g, MAXVAL( ABS( prob%Y ) ) )
        ELSE
          prob%rho_g = control%rho_g
        END IF

        IF( control%rho_b <= zero ) THEN
          prob%rho_b = prob%n
          prob%rho_b = two * MAX( prob%rho_b, MAXVAL( ABS( prob%Z ) ) )
        ELSE
          prob%rho_b = control%rho_b
        END IF

        IF ( printi ) WRITE( control%out,                                      &
          "( /, A, ' Warm-starting QPA with rho_g = ', ES8.2,                  &
       &     ', rho_b = ', ES8.2, /, A, ' ', I0, ' constraints and ', I0,      &
       &   ' bounds in the initial working set' )" ) prefix,                   &
          prob%rho_g, prob%rho_b, prefix, COUNT( C_stat( : prob%m ) /= 0 ),    &
          COUNT( B_stat( : prob%n ) /= 0 )

      END IF

      inform%QPA_inform%obj = inform%obj

!  Remove any dependent working constraints

      CALL CPU_TIME( time ) 
      CALL QPA_remove_dependent( prob%n, prob%m, prob%A%val, prob%A%col,       &
                                 prob%A%ptr, data%K, data%FACTORS,             &
                                 data%CNTL, C_stat, B_stat, data%IBREAK,       &
                                 control%QPA_control, inform%QPA_inform,       &
                                 n_depen )
      CALL CPU_TIME( dum ) ; dum = dum - time
      inform%time%preprocess = inform%time%preprocess + dum

!  Check for error exits

      inform%status = inform%QPA_inform%status
      IF ( inform%status /= GALAHAD_ok ) THEN

!  On error exit, compute the current objective function value

        data%H_s( : prob%n ) = zero
        CALL QPC_HX( data%dims, prob%n, data%H_s( : prob%n ),                  &
                     prob%H%ptr( prob%n + 1 ) - 1, prob%H%val,                 &
                     prob%H%col, prob%H%ptr, prob%X( : prob%n ), '+' )
        inform%obj = half * DOT_PRODUCT( prob%X( : prob%n ),                   &
                                         data%H_s( : prob%n ) )                &
                     + DOT_PRODUCT( prob%X( : prob%n ),                        &
                                    prob%G( : prob%n ) ) + prob%f

!  Print details of the error exit

        IF ( control%error > 0 .AND. control%print_level >= 1 ) THEN
          WRITE( control%out, "( ' ' )" )
          IF ( inform%status /= GALAHAD_ok ) WRITE( control%error, 2040 )      &
            inform%status, 'QPA_remove_dependent'
        END IF
        GO TO 700

!       IF ( control%out > 0 .AND. control%print_level >= 2 .AND. n_depen > 0 )&
!       WRITE( control%out, "(/, ' The following ',I7,' constraints appear',   &
!      &       ' to be dependent', /, ( 8I8 ) )" ) n_depen, data%Index_C_freed

      END IF
      IF ( control%out > 0 .AND. control%print_level >= 1 .AND. n_depen > 0 )  &
        WRITE( control%out, "(/, A, ' ', I0,' active constraints appear to',   &
       & ' be dependent and will be ignored ' )" ) prefix, n_depen

!  Continue allocating workspace arrays

!     m_link = MIN( prob%m + prob%n - data%dims%x_free, prob%n )
      m_link = prob%m + prob%n - data%dims%x_free
      k_n_max = prob%n + m_link

!  Allocate real workspace

      array_name = 'qpc: data%BREAKP'
      CALL SPACE_resize_array( lbreak, data%BREAKP, inform%status,             &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%A_s'
      CALL SPACE_resize_array( prob%m, data%A_s, inform%status,                &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%PERT'
      CALL SPACE_resize_array( prob%m + prob%n, data%PERT, inform%status,      &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%GRAD '
      CALL SPACE_resize_array( prob%n, data%GRAD, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%VECTOR'
      CALL SPACE_resize_array( k_n_max, data%VECTOR, inform%status,            &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%RHS'
      CALL SPACE_resize_array( k_n_max + control%QPA_control%max_sc,           &
             data%RHS, inform%status,                                          &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%S'
      CALL SPACE_resize_array( k_n_max + control%QPA_control%max_sc,           &
             data%S, inform%status,                                            &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%B'
      CALL SPACE_resize_array( k_n_max + control%QPA_control%max_sc,           &
             data%B, inform%status,                                            &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%RES'
      CALL SPACE_resize_array( k_n_max + control%QPA_control%max_sc,           &
             data%RES, inform%status,                                          &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%S_perm'
      CALL SPACE_resize_array( k_n_max + control%QPA_control%max_sc,           &
             data%S_perm, inform%status,                                       &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%DX'
      CALL SPACE_resize_array( k_n_max + control%QPA_control%max_sc,           &
             data%DX, inform%status,                                           &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%RES_print'
      CALL SPACE_resize_array( k_n_max + control%QPA_control%max_sc,           &
             data%RES_print, inform%status,                                    &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      IF ( control%QPA_control%precon >= 0 ) THEN
        n_pcg = prob%n
      ELSE
        n_pcg = 0
      END IF

      array_name = 'qpc: data%R_pcg'
      CALL SPACE_resize_array( n_pcg, data%R_pcg, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%X_pcg'
      CALL SPACE_resize_array(n_pcg , data%X_pcg, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%P_pcg'
      CALL SPACE_resize_array( n_pcg, data%P_pcg, inform%status,               &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

!  Allocate integer workspace arrays

      array_name = 'qpc: data%SC'
      CALL SPACE_resize_array(  control%QPA_control%max_sc + 1,                &
             data%SC, inform%status,                                           &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%REF'
      CALL SPACE_resize_array( m_link, data%REF, inform%status,                &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%C_up_or_low'
      CALL SPACE_resize_array( data%dims%c_u_start, data%dims%c_l_end,         &
             data%C_up_or_low, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%X_up_or_low'
      CALL SPACE_resize_array( data%dims%x_u_start, data%dims%x_l_end,         &
             data%X_up_or_low, inform%status,                                  &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%PERM'
      CALL SPACE_resize_array( k_n_max + control%QPA_control%max_sc,           &
             data%PERM, inform%status,                                         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

!  Find the total length of the control%QPA_control%max_sc largest rows

      DO i = 1, prob%m
        data%IBREAK( i ) = prob%A%ptr( i ) - prob%A%ptr( i + 1 )
      END DO
      CALL SORT_heapsort_build( prob%m, data%IBREAK( : prob%m ),     &
                                inform%status )
      lbd = 0
      DO i = 1, MIN( control%QPA_control%max_sc, prob%m )
        ii = prob%m - i + 1
        CALL SORT_heapsort_smallest( ii, data%IBREAK( : ii ), inform%status )
        lbd = lbd - data%IBREAK( ii )
      END DO
      IF ( control%QPA_control%max_sc > prob%m )                               &
        lbd = lbd + control%QPA_control%max_sc - prob%m

!  Allocate arrays

      array_name = 'qpc: data%SCU_mat%BD_col_start'
      CALL SPACE_resize_array( control%QPA_control%max_sc + 1,                 &
             data%SCU_mat%BD_col_start, inform%status,                         &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%SCU_mat%BD_val'
      CALL SPACE_resize_array( lbd, data%SCU_mat%BD_val, inform%status,        &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%SCU_mat%BD_row'
      CALL SPACE_resize_array( lbd, data%SCU_mat%BD_row, inform%status,        &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

      array_name = 'qpc: data%DIAG'
      CALL SPACE_resize_array( 2, K_n_max, data%DIAG, inform%status,           &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) GO TO 900

!  decide on appropriate initial preconditioners and factorizations

      data%auto_prec = control%QPA_control%precon == 0
      data%auto_fact = control%QPA_control%factor == 0

!  If the Hessian has semi-bandwidth smaller than nsemib and the preconditioner 
!  is to be picked automatically, use the full Hessian. Otherwise, use the
!  Hessian of the specified semi-bandwidth.

      IF ( data%auto_prec ) THEN
        data%prec_hist = 2

!  prec_hist indicates which factors are currently being used. Possible values:
!   1 full factors used
!   2 band factors used
!   3 diagonal factors used (as a last resort)

!  Check to see if the Hessian is banded

 dod :  DO type = 1, 6
        
          SELECT CASE( type )
          CASE ( 1 )
        
            hd_start  = 1
            hd_end    = data%dims%h_diag_end_free
            hnd_start = hd_end + 1
            hnd_end   = data%dims%x_free
        
          CASE ( 2 )
        
            hd_start  = data%dims%x_free + 1
            hd_end    = data%dims%h_diag_end_nonneg
            hnd_start = hd_end + 1
            hnd_end   = data%dims%x_l_start - 1
        
          CASE ( 3 )
        
            hd_start  = data%dims%x_l_start
            hd_end    = data%dims%h_diag_end_lower
            hnd_start = hd_end + 1
            hnd_end   = data%dims%x_u_start - 1
        
          CASE ( 4 )
        
            hd_start  = data%dims%x_u_start
            hd_end    = data%dims%h_diag_end_range
            hnd_start = hd_end + 1
            hnd_end   = data%dims%x_l_end
        
          CASE ( 5 )
        
            hd_start  = data%dims%x_l_end + 1
            hd_end    = data%dims%h_diag_end_upper
            hnd_start = hd_end + 1
            hnd_end   = data%dims%x_u_end
        
          CASE ( 6 )
        
            hd_start  = data%dims%x_u_end + 1
            hd_end    = data%dims%h_diag_end_nonpos
            hnd_start = hd_end + 1
            hnd_end   = prob%n
        
          END SELECT
    
!  rows with a diagonal entry
    
          hd_end = MIN( hd_end, prob%n )
          DO i = hd_start, hd_end
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 2
              IF ( ABS( i - prob%H%col( l ) ) > control%QPA_control%nsemib ) THEN
                data%prec_hist = 1
                EXIT dod
              END IF  
            END DO
          END DO
          IF ( hd_end == prob%n ) EXIT
    
!  rows without a diagonal entry
    
          hnd_end = MIN( hnd_end, prob%n )
          DO i = hnd_start, hnd_end
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              IF ( ABS( i - prob%H%col( l ) ) > control%QPA_control%nsemib ) THEN
                data%prec_hist = 1
                EXIT dod
              END IF  
            END DO
          END DO
          IF ( hd_end == prob%n ) EXIT
    
        END DO dod

      END IF

!  =============
!  Now crossover
!  =============

      IF ( printi ) WRITE( control%out, "( /, A, ' entering QPA ' )" ) prefix

      CALL CPU_TIME( time_qpa_start )
      CALL QPA_solve_qp( data%dims, prob%n, prob%m,                            &
                         prob%H%val, prob%H%col, prob%H%ptr,                   &
                         prob%G, prob%f, prob%rho_g, prob%rho_b, prob%A%val,   &
                         prob%A%col, prob%A%ptr, prob%C_l, prob%C_u, prob%X_l, &
                         prob%X_u, prob%X, prob%Y, prob%Z, C_stat, B_stat,     &
                         m_link, K_n_max, lbreak, data%RES_l, data%RES_u,      &
                         data%A_norms, data%H_s, data%BREAKP, data%A_s,        &
                         data%PERT, data%GRAD, data%VECTOR, data%RHS, data%S,  &
                         data%B, data%RES, data%S_perm, data%DX, n_pcg,        &
                         data%R_pcg, data%X_pcg, data%P_pcg, data%Abycol_val,  &
                         data%Abycol_row, data%Abycol_ptr, data%S_val,         &
                         data%S_row, data%S_col, data%S_colptr, data%IBREAK,   &
                         data%SC, data%REF, data%RES_print, data%DIAG,         &
                         data%C_up_or_low, data%X_up_or_low, data%PERM,        &
                         data%FACTORS, data%CNTL,                              &
                         data%AINFO, data%FINFO,                               &
                         data%SCU_mat, data%SCU_info, data%SCU_data, data%K,   &
                         data%seed, time_qpa_start,                            &
                         data%start_print, data%stop_print,                    &
                         data%prec_hist, data%auto_prec, data%auto_fact,       &
                         control%QPA_control%print_level > 0 .AND.             &
                         control%QPA_control%out > 0,                          &
                         control%QPA_control, inform%QPA_inform,               &
                         G_p = G_p, X_p = X_p, Y_p = Y_p, Z_p = Z_p )

      CALL CPU_TIME( time )
      inform%QPA_inform%time%total = time - time_qpa_start
      IF ( printi ) THEN
        IF ( control%QPA_control%out > 0 .AND.                                 &
             control%QPA_control%print_level > 0 )                             &
          WRITE( control%out, "( /, A, ' returned from QPA' )" ) prefix
        WRITE( control%out, "( A, ' on exit from QPA: status = ', I0,          &
       &   ', iterations = ', I0, ', time = ', F0.2 )" ) prefix,               &
          inform%QPA_inform%status, inform%QPA_inform%iter,                    &
          inform%QPA_inform%time%total
        WRITE( control%out, "( A, ' objective value =', ES12.4,                &
       &   /, A, ' # active counstraints: ', I0, ' from ', I0,                 & 
       &   ', bounds: ', I0, ' from ', I0 )" ) prefix, inform%QPA_inform%obj,  &
            prefix, COUNT( C_stat( : prob%m ) /= 0 ),  prob%m,                 &
            COUNT( B_stat( : prob%n ) /= 0 ),  prob%n
      END IF

!  Record output information from QPA

!     IF ( .NOT. inform%QPB_inform%feasible )                                  &
        inform%status = inform%QPA_inform%status
      inform%alloc_status = inform%QPA_inform%alloc_status
      inform%nfacts = inform%nfacts + inform%QPA_inform%nfacts
      inform%nmods = inform%nmods + inform%QPA_inform%nmods
      inform%obj = inform%QPA_inform%obj

!  Record times for components of QPA

      inform%time%total = inform%time%total + inform%QPA_inform%time%total
      inform%time%analyse =                                                    &
        inform%time%analyse + inform%QPA_inform%time%analyse
      inform%time%factorize =                                                  &
        inform%time%factorize + inform%QPA_inform%time%factorize
      inform%time%solve =inform%time%solve + inform%QPA_inform%time%solve
      inform%time%preprocess =                                                 &
        inform%time%preprocess + inform%QPA_inform%time%preprocess

      IF ( printi ) THEN
        WRITE( control%out, "( A, ' ', I0, ' integer and ', I0,                &
      &  ' real words required for factors in QPA' )" ) prefix,                &
         inform%QPA_inform%factorization_integer,                              &
         inform%QPA_inform%factorization_real
      END IF

!  If the crosover fails for any reason, use the active-set predictions from
!  the original QP solve 

      IF ( inform%QPA_inform%status < 0 ) THEN
        DO i = 1, prob%m
          IF ( data%IW( i ) > 0 ) THEN
            C_stat( i ) = 1
          ELSE IF ( data%IW( i ) < 0 ) THEN
            C_stat( i ) = - 1
          ELSE
            C_stat( i ) = 0
          END IF
        END DO
        DO i = 1, prob%n
          IF ( data%IW( prob%m + i ) > 0 ) THEN
            B_stat( i ) = 1
          ELSE IF ( data%IW( prob%m + i ) < 0 ) THEN
            B_stat( i ) = - 1
          ELSE
            B_stat( i ) = 0
          END IF
        END DO

!  If possible, restore the solution from QPB

        IF ( gotsol ) THEN
          prob%X( : prob%n ) = data%X_trial( : prob%n )
          prob%Y( : prob%m ) = data%Y_last( : prob%m )
          prob%Z( : prob%n ) = data%Z_last( : prob%n )
          prob%C( : prob%m ) = data%C( : prob%m )
          inform%obj = best_obj
          inform%status = GALAHAD_ok
        END IF
      ELSE
        inform%p_found = .TRUE.
      END IF

      IF ( control%qpb_or_qpa .OR. .NOT. control%no_qpb ) THEN
        IF ( printd ) THEN
          DO i = 1, prob%m 
            IF ( data%IW( i ) /= C_stat( i ) ) WRITE( control%out,             &
              "( ' C_stat ', I6, ' QPB/LSQP = ', I2, ' vs QPA = ', I2 )" )     & 
                i, data%IW( i ), C_stat( i )
          END DO
          DO i = 1, prob%n
            IF ( data%IW( prob%m + i ) /= B_stat( i ) ) WRITE( control%out,    &
              "( ' B_stat ', I6, ' QPB/LSQP = ', I2, ' vs QPA = ', I2 )" )     &
                i, data%IW( prob%m + i ), B_stat( i )
          END DO
        ELSE
          IF ( printi ) WRITE( control%out,                                    &
            "( A, ' # changes to C_stat, B_stat: QPB/LSQP vs QPA = ', I0,      &
            &  ', ', I0 )" ) prefix,                                           &
             COUNT( C_stat( : prob%m ) /= data%IW( : prob%m ) ),               &
             COUNT( B_stat( : prob%n ) /=                                      &
                   data%IW( prob%m + 1 : prob%m + prob%n ) )
        END IF

        DO i = 1, prob%n
          IF ( data%IW( prob%m + i ) /= B_stat( i ) ) THEN
            IF ( printt ) WRITE( control%out,                                  &
              "( '  B_stat ', I6, ' QPB/LSQP = ', I2, ' vs QPA = ', I2, /,     &
           &     '  dist_x_l, dist_x_u, z = ', 3ES16.8 )" ) i,                 &
              data%IW( prob%m + i ), B_stat( i ), prob%X( i ) - prob%X_l( i ), &
              prob%X_u( i ) - prob%X( i ), prob%Z( i )
            IF ( ABS( prob%X( i ) - prob%X_l( i ) ) <= control%on_bound_tol )  &
              B_stat( i ) = - 1
            IF ( ABS( prob%X( i ) - prob%X_u( i ) ) <= control%on_bound_tol )  &
              B_stat( i ) = 1
          END IF
        END DO
        DO i = 1, prob%m 
          IF ( data%IW( i ) /= C_stat( i ) ) THEN
            IF ( printt ) WRITE( control%out,                                  &
              "( '  C_stat ', I6, ' QPB/LSQP = ', I2, ' vs QPA = ', I2, /,     &
           &     '  dist_c_l, dist_c_u, y = ', 3ES16.8 )" ) i,                 &
              data%IW( i ), C_stat( i ), prob%C( i ) - prob%C_l( i ),          &
              prob%C_u( i ) - prob%C( i ), prob%Y( i )
            IF ( ABS( prob%C( i ) - prob%C_l( i ) ) <= control%on_bound_tol )  &
              C_stat( i ) = - 1
            IF ( ABS( prob%C( i ) - prob%C_u( i ) ) <= control%on_bound_tol )  &
              C_stat( i ) = 1
          END IF
        END DO
      END IF

  700 CONTINUE 

!     IF ( control%yes ) THEN
!      OPEN( 56 )
!      WRITE( 56, "( '    i st        dist_x_l             dist_x_u        z' )")
!      DO i = 1, prob%n
!        WRITE( 56, "( I6, I3, 3ES22.14 )" ) i, B_stat( i ),                    &
!         prob%X( i ) - prob%X_l( i ), prob%X_u( i ) - prob%X( i ), prob%Z( i )
!      END DO

!      WRITE( 56, "( '    i st        dist_c_l             dist_c_u        y' )")
!      DO i = 1, prob%m
!        WRITE( 56, "( I6, I3, 3ES22.14 )" ) i, C_stat( i ),                    &
!         prob%C( i ) - prob%C_l( i ), prob%C_u( i ) - prob%C( i ), prob%Y( i )
!      END DO
!      CLOSE( 56 )
!     END IF

!  If some of the constraints were freed having first been fixed during 
!  the computation, refix them now

      IF ( remap_more_freed ) THEN
        CALL CPU_TIME( time )
        IF ( inform%p_found ) THEN
          IF ( PRESENT( X_p ) ) THEN
            X_p( prob%n + 1 : data%QPP_map_more_freed%n ) = zero
            CALL SORT_inverse_permute( data%QPP_map_more_freed%n,              &
                                       data%QPP_map_more_freed%x_map,          &
                                       X = X_p( : data%QPP_map_more_freed%n ) )
          END IF
          IF ( PRESENT( Y_p ) ) THEN
            Y_p( prob%m + 1 : data%QPP_map_more_freed%m ) = zero
            CALL SORT_inverse_permute( data%QPP_map_more_freed%m,              &
                                       data%QPP_map_more_freed%c_map,          &
                                       X = Y_p( : data%QPP_map_more_freed%m ) )
          END IF
          IF ( PRESENT( Z_p ) ) THEN
            Z_p( prob%n + 1 : data%QPP_map_more_freed%n ) = zero
            CALL SORT_inverse_permute( data%QPP_map_more_freed%n,              &
                                       data%QPP_map_more_freed%x_map,          &
                                       X = Z_p( : data%QPP_map_more_freed%n ) )
          END IF
        END IF
        C_stat( prob%m + 1 : data%QPP_map_more_freed%m ) = 0
        CALL SORT_inverse_permute( data%QPP_map_more_freed%m,                  &
                                   data%QPP_map_more_freed%c_map,              &
                                   IX = C_stat( : data%QPP_map_more_freed%m ) )
        B_stat( prob%n + 1 : data%QPP_map_more_freed%n ) = - 1
!       B_stat( prob%n + 1 : data%QPP_map_more_freed%n ) = 0
        CALL SORT_inverse_permute( data%QPP_map_more_freed%n,                  &
                                   data%QPP_map_more_freed%x_map,              &
                                   IX = B_stat( : data%QPP_map_more_freed%n ) )
        CALL QPP_restore( data%QPP_map_more_freed, data%QPP_inform,            &
                          prob, get_all = .TRUE. )
        CALL QPP_terminate( data%QPP_map_more_freed, data%QPP_control,         &
                            data%QPP_inform )
        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum
        data%dims = data%dims_save_more_freed

!  Fix the temporarily freed constraint bounds

        DO i = 1, n_more_depen
          j = data%Index_c_more_freed( i )
          prob%C_l( j ) = data%C_more_freed( i )
          prob%C_u( j ) = data%C_more_freed( i )
        END DO
      END IF

!  If some of the variables/constraints were fixed during the computation,
!  free them now

  720 CONTINUE
      IF ( remap_fixed ) THEN
        CALL CPU_TIME( time )
        IF ( inform%p_found ) THEN
          IF ( PRESENT( X_p ) ) THEN
            X_p( prob%n + 1 : data%QPP_map_fixed%n ) = zero
            CALL SORT_inverse_permute( data%QPP_map_fixed%n,                   &
                                       data%QPP_map_fixed%x_map,               &
                                       X = X_p( : data%QPP_map_fixed%n ) )
          END IF
          IF ( PRESENT( Y_p ) ) THEN
            Y_p( prob%m + 1 : data%QPP_map_fixed%m ) = zero
            CALL SORT_inverse_permute( data%QPP_map_fixed%m,                   &
                                       data%QPP_map_fixed%c_map,               &
                                       X = Y_p( : data%QPP_map_fixed%m ) )
          END IF
          IF ( PRESENT( Z_p ) ) THEN
            Z_p( prob%n + 1 : data%QPP_map_fixed%n ) = zero
            CALL SORT_inverse_permute( data%QPP_map_fixed%n,                   &
                                       data%QPP_map_fixed%x_map,               &
                                       X = Z_p( : data%QPP_map_fixed%n ) )
          END IF
        END IF
        C_stat( prob%m + 1 : data%QPP_map_fixed%m ) = 0
        CALL SORT_inverse_permute( data%QPP_map_fixed%m,                       &
                                   data%QPP_map_fixed%c_map,                   &
                                   IX = C_stat( : data%QPP_map_fixed%m ) )
        B_stat( prob%n + 1 : data%QPP_map_fixed%n ) = - 1
!       B_stat( prob%n + 1 : data%QPP_map_fixed%n ) = 0
        CALL SORT_inverse_permute( data%QPP_map_fixed%n,                       &
                                   data%QPP_map_fixed%x_map,                   &
                                   IX = B_stat( : data%QPP_map_fixed%n ) )
        CALL QPP_restore( data%QPP_map_fixed, data%QPP_inform,                 &
                          prob, get_all = .TRUE. )
        CALL QPP_terminate( data%QPP_map_fixed, data%QPP_control,              &
                            data%QPP_inform )
        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum
        data%dims = data%dims_save_fixed

!  Release the temporarily fixed problem bounds

        DO i = 1, tiny_x
          j = data%Index_X_fixed( i )
          IF ( j > 0 ) THEN
            prob%X_l( j ) = data%X_fixed( i )
            IF ( B_stat( j ) < 0 ) THEN
              B_stat( j ) = - B_stat( j )
              prob%X( j ) = -  prob%X( j )
            END IF  
          ELSE
            prob%X_u( - j ) = data%X_fixed( i )
          END IF
        END DO
       
!  Do the same for the constraint bounds

        DO i = 1, tiny_c
          j = data%Index_C_fixed( i )
          IF ( j > 0 ) THEN
            prob%C_l( j ) = data%C_fixed( i )
            IF ( C_stat( j ) < 0 ) THEN
              C_stat( j ) = - C_stat( j )
              prob%Y( j ) = -  prob%Y( j )
            END IF  
          ELSE
            prob%C_u( - j ) = data%C_fixed( i )
          END IF
        END DO

      END IF

!  If some of the constraints were freed during the computation, refix them now

      IF ( remap_freed ) THEN
        CALL CPU_TIME( time )
        IF ( inform%p_found ) THEN
          IF ( PRESENT( X_p ) ) THEN
            X_p( prob%n + 1 : data%QPP_map_freed%n ) = zero
            CALL SORT_inverse_permute( data%QPP_map_freed%n,                   &
                                       data%QPP_map_freed%x_map,               &
                                       X = X_p( : data%QPP_map_freed%n ) )
          END IF
          IF ( PRESENT( Y_p ) ) THEN
            Y_p( prob%m + 1 : data%QPP_map_freed%m ) = zero
            CALL SORT_inverse_permute( data%QPP_map_freed%m,                   &
                                       data%QPP_map_freed%c_map,               &
                                       X = Y_p( : data%QPP_map_freed%m ) )
          END IF
          IF ( PRESENT( Z_p ) ) THEN
            Z_p( prob%n + 1 : data%QPP_map_freed%n ) = zero
            CALL SORT_inverse_permute( data%QPP_map_freed%n,                   &
                                       data%QPP_map_freed%x_map,               &
                                       X = Z_p( : data%QPP_map_freed%n ) )
          END IF
        END IF
        C_stat( prob%m + 1 : data%QPP_map_freed%m ) = 0
        CALL SORT_inverse_permute( data%QPP_map_freed%m,                       &
                                   data%QPP_map_freed%c_map,                   &
                                   IX = C_stat( : data%QPP_map_freed%m ) )
!       B_stat( prob%n + 1 : data%QPP_map_freed%n ) = 0
        B_stat( prob%n + 1 : data%QPP_map_freed%n ) = - 1
        CALL SORT_inverse_permute( data%QPP_map_freed%n,                       &
                                   data%QPP_map_freed%x_map,                   &
                                   IX = B_stat( : data%QPP_map_freed%n ) )
        CALL QPP_restore( data%QPP_map_freed, data%QPP_inform,                 &
                          prob, get_all = .TRUE. )
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

  750 CONTINUE 
      data%trans = data%trans - 1
      IF ( data%trans == 0 ) THEN
        CALL CPU_TIME( time )
        IF ( inform%p_found ) THEN
          IF ( PRESENT( X_p ) ) THEN
            X_p( prob%n + 1 : data%QPP_map%n ) = zero
            CALL SORT_inverse_permute( data%QPP_map%n, data%QPP_map%x_map,     &
                                       X = X_p( : data%QPP_map%n ) )
          END IF
          IF ( PRESENT( Y_p ) ) THEN
            Y_p( prob%m + 1 : data%QPP_map%m ) = zero
            CALL SORT_inverse_permute( data%QPP_map%m, data%QPP_map%c_map,     &
                                       X = Y_p( : data%QPP_map%m ) )
          END IF
          IF ( PRESENT( Z_p ) ) THEN
            Z_p( prob%n + 1 : data%QPP_map%n ) = zero
            CALL SORT_inverse_permute( data%QPP_map%n, data%QPP_map%x_map,     &
                                       X = Z_p( : data%QPP_map%n ) )
          END IF
        END IF
        CALL SORT_inverse_permute( data%QPP_map%m, data%QPP_map%c_map,         &
                                   IX = C_stat( : data%QPP_map%m ) )
        B_stat( prob%n + 1 : data%QPP_map%n ) = - 1
        CALL SORT_inverse_permute( data%QPP_map%n, data%QPP_map%x_map,         &
                                   IX = B_stat( : data%QPP_map%n ) )

!  Full restore

        IF ( control%restore_problem >= 2 ) THEN
          CALL QPP_restore( data%QPP_map, data%QPP_inform, prob,               &
                            get_all = .TRUE. )

!  Restore vectors and scalars

        ELSE IF ( control%restore_problem == 1 ) THEN
          CALL QPP_restore( data%QPP_map, data%QPP_inform, prob,               &
                            get_f = .TRUE., get_g = .TRUE.,                    &
                            get_x = .TRUE., get_x_bounds = .TRUE.,             &
                            get_y = .TRUE., get_z = .TRUE.,                    &
                            get_c = .TRUE., get_c_bounds = .TRUE. )

!  Solution recovery

        ELSE
          CALL QPP_restore( data%QPP_map, data%QPP_inform, prob,               &
                             get_x = .TRUE., get_y = .TRUE., get_z = .TRUE.,   &
                             get_c = .TRUE. )
        END IF
        CALL QPP_terminate( data%QPP_map, data%QPP_control, data%QPP_inform )
        CALL CPU_TIME( dum ) ; dum = dum - time
        inform%time%preprocess = inform%time%preprocess + dum
        prob%new_problem_structure = data%new_problem_structure
        data%save_structure = .TRUE.
      END IF

!  Compute total time

      CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
      IF ( control%out > 0 .AND. control%print_level >= 2 )                    &
        WRITE( control%out, 2000 )                                             &
          inform%time%total, inform%time%preprocess,                           &
          inform%time%analyse, inform%time%factorize, inform%time%solve,       &
          inform%QPA_inform%time%total, inform%QPB_inform%time%total,          &
          inform%QPA_inform%time%analyse, inform%QPA_inform%time%factorize,    &
          inform%QPA_inform%time%solve, inform%QPB_inform%time%analyse,        &
          inform%QPB_inform%time%factorize, inform%QPB_inform%time%solve

  800 CONTINUE 
      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' leaving QPC_solve ' )" )

      RETURN  

!  Allocation error

  900 CONTINUE 
      CALL CPU_TIME( time ) ; inform%time%total = time - time_start 
      IF ( control%out > 0 .AND. control%print_level >= 5 )                    &
        WRITE( control%out, "( ' leaving QPC_solve ' )" )

      RETURN  

!  Non-executable statements

 2000 FORMAT( /, ' =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=',               &
              '-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=',                             &
              /, ' =', 5X, 'QPC timing statistics:  total', 0P, F12.2,         &
                 ' preprocess', F12.2, 5X, '=',                                &
              /, ' =', 23X, 'analyse   factorize    solve', 23X, '=',          &
              /, ' =', 18X, 3F11.2, 23X, '=',                                  &
              /, ' =', 12X, 'QPA: total',  F12.2,                              &
                       14X, 'QPB: total',  F12.2, 4X, '=',                     & 
              /, ' =      analyse    factorize     solve',                     &
                 '      analyse    factorize     solve  =',                    &
              /, ' =', 6F12.2, 2x, '=',                                        &
              /, ' =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=',               &
                 '-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=' )
 2010 FORMAT( ' ', /, '   **  Error return ',I3,' from QPC ' ) 
 2040 FORMAT( '   **  Error return ', I6, ' from ', A ) 
 2240 FORMAT( /, '  Warning - an entry from strict upper triangle of H given ' )

!  End of QPC_solve

      END SUBROUTINE QPC_solve

!-*-*-*-*-*-*-   Q P C _ T E R M I N A T E   S U B R O U T I N E   -*-*-*-*-*-*

      SUBROUTINE QPC_terminate( data, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!      ..............................................
!      .                                            .
!      .  Deallocate internal arrays at the end     .
!      .  of the computation                        .
!      .                                            .
!      ..............................................

!  Arguments:
!  =========
!
!   data    see Subroutine QPC_initialize
!   control see Subroutine QPC_initialize
!   inform  see Subroutine QPC_solve

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( QPC_data_type ), INTENT( INOUT ) :: data
      TYPE ( QPC_control_type ), INTENT( IN ) :: control        
      TYPE ( QPC_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      CHARACTER ( LEN = 80 ) :: array_name

!  Deallocate all arrays allocated by QPA and QPB

      CALL QPA_terminate( data, control%QPA_control, inform%QPA_inform )
      IF ( inform%QPA_inform%status /= GALAHAD_ok ) THEN
        inform%status = GALAHAD_error_deallocate
        inform%alloc_status = inform%QPA_inform%alloc_status
!       inform%bad_alloc = inform%QPA_inform%bad_alloc
        IF ( control%deallocate_error_fatal ) RETURN
      END IF

      CALL QPB_terminate( data, control%QPB_control, inform%QPB_inform )
      IF ( inform%QPB_inform%status /= GALAHAD_ok ) THEN
        inform%status = GALAHAD_error_deallocate
        inform%alloc_status = inform%QPB_inform%alloc_status
!       inform%bad_alloc = inform%QPB_inform%bad_alloc
        IF ( control%deallocate_error_fatal ) RETURN
      END IF

!  Deallocate all remaining allocated arrays

      array_name = 'qpc: data%HX'
      CALL SPACE_dealloc_array( data%HX,                                       &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%C_freed'
      CALL SPACE_dealloc_array( data%C_freed,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

!     array_name = 'qpc: prob%X0'
!     CALL SPACE_dealloc_array( prob%X0,                                       &
!        inform%status, inform%alloc_status, array_name = array_name,          &
!        bad_alloc = inform%bad_alloc, out = control%error )
!     IF ( control%deallocate_error_fatal .AND.                                &
!          inform%status /= GALAHAD_ok ) RETURN

!     array_name = 'qpc: prob%WEIGHT'
!     CALL SPACE_dealloc_array( prob%WEIGHT,                                   &
!        inform%status, inform%alloc_status, array_name = array_name,          &
!        bad_alloc = inform%bad_alloc, out = control%error )
!     IF ( control%deallocate_error_fatal .AND.                                &
!          inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%X0'
      CALL SPACE_dealloc_array( data%X0,                                       &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%X_fixed'
      CALL SPACE_dealloc_array( data%X_fixed,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%Index_X_fixed'
      CALL SPACE_dealloc_array( data%Index_X_fixed,                            &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%C_fixed'
      CALL SPACE_dealloc_array( data%C_fixed,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%Index_C_fixed'
      CALL SPACE_dealloc_array( data%Index_C_fixed,                            &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%HX'
      CALL SPACE_dealloc_array( data%HX,                                       &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%C_more_freed'
      CALL SPACE_dealloc_array( data%C_more_freed,                             &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%DZ_l'
      CALL SPACE_dealloc_array( data%DZ_l,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%DZ_u'
      CALL SPACE_dealloc_array( data%DZ_u,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%GRAD'
      CALL SPACE_dealloc_array( data%GRAD,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%X_trial'
      CALL SPACE_dealloc_array( data%X_trial,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%GRAD_X_phi'
      CALL SPACE_dealloc_array( data%GRAD_X_phi,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%GRAD_C_phi'
      CALL SPACE_dealloc_array( data%GRAD_C_phi,                               &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%S'
      CALL SPACE_dealloc_array( data%S,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%A_s'
      CALL SPACE_dealloc_array( data%A_s,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%H_s'
      CALL SPACE_dealloc_array( data%H_s,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%Y_last'
      CALL SPACE_dealloc_array( data%Y_last,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%Z_last'
      CALL SPACE_dealloc_array( data%Z_last,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%IW'
      CALL SPACE_dealloc_array( data%IW,                                       &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%IBREAK'
      CALL SPACE_dealloc_array( data%IBREAK,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%RES_l'
      CALL SPACE_dealloc_array( data%RES_l,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%RES_u'
      CALL SPACE_dealloc_array( data%RES_u,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%A_norms'
      CALL SPACE_dealloc_array( data%A_norms,                                  &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%H_s'
      CALL SPACE_dealloc_array( data%H_s,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%BREAKP'
      CALL SPACE_dealloc_array( data%BREAKP,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%A_s'
      CALL SPACE_dealloc_array( data%A_s,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%PERT'
      CALL SPACE_dealloc_array( data%PERT,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%GRAD'
      CALL SPACE_dealloc_array( data%GRAD,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%VECTOR'
      CALL SPACE_dealloc_array( data%VECTOR,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%RHS'
      CALL SPACE_dealloc_array( data%RHS,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%S'
      CALL SPACE_dealloc_array( data%S,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%B'
      CALL SPACE_dealloc_array( data%B,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%RES'
      CALL SPACE_dealloc_array( data%RES,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%S_perm'
      CALL SPACE_dealloc_array( data%S_perm,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%DX'
      CALL SPACE_dealloc_array( data%DX,                                       &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%RES_print'
      CALL SPACE_dealloc_array( data%RES_print,                                &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%R_pcg'
      CALL SPACE_dealloc_array( data%R_pcg,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%X_pcg'
      CALL SPACE_dealloc_array( data%X_pcg,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%P_pcg'
      CALL SPACE_dealloc_array( data%P_pcg,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%SC'
      CALL SPACE_dealloc_array( data%SC,                                       &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%REF'
      CALL SPACE_dealloc_array( data%REF,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%C_up_or_low'
      CALL SPACE_dealloc_array( data%C_up_or_low,                              &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%X_up_or_low'
      CALL SPACE_dealloc_array( data%X_up_or_low,                              &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%PERM'
      CALL SPACE_dealloc_array( data%PERM,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%SCU_mat%BD_col_start'
      CALL SPACE_dealloc_array( data%SCU_mat%BD_col_start,                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%SCU_mat%BD_val'
      CALL SPACE_dealloc_array( data%SCU_mat%BD_val,                           &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%SCU_mat%BD_row'
      CALL SPACE_dealloc_array( data%SCU_mat%BD_row,                           &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'qpc: data%DIAG'
      CALL SPACE_dealloc_array( data%DIAG,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      RETURN

!  End of subroutine QPC_terminate

      END SUBROUTINE QPC_terminate

!  End of module QPC

   END MODULE GALAHAD_QPC_double
