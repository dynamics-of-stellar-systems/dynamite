! THIS VERSION: GALAHAD 2.2 - 22/04/2008 AT 16:00 GMT.

!-*-*-*-*-*-*-*-*-*- G A L A H A D _ E Q P   M O D U L E -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   development started March 25th 2004
!   originally released GALAHAD Version 2.0. February 16th 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

  MODULE GALAHAD_EQP_double

!      -------------------------------------------------
!     |                                                 |
!     | Solve the equaity-constrained quadratic program |
!     |                                                 |
!     |   minimize     1/2 x(T) H x + g(T) x + f        |
!     |   subject to     A x + c = 0                    |
!     |   and            ||x|| <= Delta                 |
!     |                                                 |
!     | using a projected preconditined CG method       |
!     |                                                 |
!      -------------------------------------------------

      USE GALAHAD_SYMBOLS
!NOT95USE GALAHAD_CPU_time
      USE GALAHAD_SYMBOLS
      USE GALAHAD_SPACE_double
!     USE GALAHAD_SMT_double
      USE GALAHAD_QPT_double
      USE GALAHAD_SBLS_double
      USE GALAHAD_GLTR_double
      USE GALAHAD_SPECFILE_double
   
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: EQP_initialize, EQP_read_specfile, EQP_solve, EQP_terminate,   &
                EQP_solve_main, QPT_problem_type, SMT_type, SMT_put, SMT_get

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

      TYPE, PUBLIC :: EQP_time_type
        REAL :: total, factorize, solve, solve_inter
      END TYPE

      TYPE, PUBLIC :: EQP_control_type
        INTEGER :: error, out, print_level
        INTEGER :: factorization, max_col, indmin, valmin, len_ulsmin, itref_max
        INTEGER :: cg_maxit, preconditioner, new_a, new_h, semi_bandwidth
        REAL ( KIND = wp ) :: pivot_tol, pivot_tol_for_basis, zero_pivot
        REAL ( KIND = wp ) :: inner_fraction_opt, radius
        REAL ( KIND = wp ) :: max_infeasibility_relative
        REAL ( KIND = wp ) :: max_infeasibility_absolute
        REAL ( KIND = wp ) :: inner_stop_relative, inner_stop_absolute
        REAL ( KIND = wp ) :: inner_stop_inter
        LOGICAL :: remove_dependencies, find_basis_by_transpose
        LOGICAL :: space_critical, deallocate_error_fatal
        CHARACTER ( LEN = 30 ) :: prefix
        TYPE ( SBLS_control_type ) :: SBLS_control
        TYPE ( GLTR_control_type ) :: GLTR_control
      END TYPE

      TYPE, PUBLIC :: EQP_data_type
        LOGICAL :: new_c
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: G_f
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: R
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: S
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: VECTOR
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: WORK
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RESID
        TYPE ( SBLS_data_type ) :: SBLS_data
        TYPE ( SMT_type ) :: C
        TYPE ( GLTR_data_type ) :: GLTR_data
      END TYPE

      TYPE, PUBLIC :: EQP_inform_type
        INTEGER :: status, alloc_status, cg_iter, cg_iter_inter
        INTEGER :: factorization_integer, factorization_real
        REAL ( KIND = wp ) :: obj
        TYPE ( EQP_time_type ) :: time
        TYPE ( SBLS_inform_type ) :: SBLS_inform
        TYPE ( GLTR_info_type ) :: GLTR_inform
        CHARACTER ( LEN = 80 ) :: bad_alloc
      END TYPE

!----------------------
!   P a r a m e t e r s
!----------------------

      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: point01 = 0.01_wp
      REAL ( KIND = wp ), PARAMETER :: point1 = 0.1_wp
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: two = 2.0_wp
      REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp
      REAL ( KIND = wp ), PARAMETER :: hundred = 100.0_wp
      REAL ( KIND = wp ), PARAMETER :: epsmch = EPSILON( one )

   CONTAINS

!-*-*-*-*-*-   E Q P _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE EQP_initialize( data, control )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Default control data for EQP. This routine should be called before
!  EQP_solve
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
!   factorization. The factorization to be used.
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
!   len_ulsmin. An initial guess as to the workspace required by ULS
!
!   itref_max. The maximum number of iterative refinements allowed
!
!   cg_maxit. The maximum number of CG iterations allowed. If cg_maxit < 0,
!     this number will be reset to the dimension of the system + 1
!
!   preconditioner. The preconditioner to be used for the CG is defined by 
!    preconditioner. Possible values are
!
!    variable:
!      0  automatic 
!
!    explicit factorization:
!
!      1  no preconditioner, G = I
!      2  full factorization, G = H
!      3  diagonal, G = diag( H )
!      4  banded, G = band( H ) with semi-bandwidth %semi_bandwidth
!     11  G_11 = 0, G_21 = 0, G_22 = H_22
!     12  G_11 = 0, G_21 = H_21, G_22 = H_22
!
!    implicit factorization:
!
!      -1  G_11 = 0, G_21 = 0, G_22 = I
!      -2  G_11 = 0, G_21 = 0, G_22 = H_22
!
!   semi_bandwidth. The semi-bandwidth of a band preconditioner, if appropriate
!
!   new_a. How much of A has changed since the previous factorization.
!    Possible values are
!
!      0  unchanged
!      1  values but not indices have changed
!      2  values and indices have changed 
!
!   new_h. How much of H has changed since the previous factorization.
!    Possible values are
!
!      0  unchanged
!      1  values but not indices have changed
!      2  values and indices have changed
!
!  REAL control parameters:
!
!   pivot_tol. The threshold pivot used by the matrix factorization.
!    See the documentation for SILS for details
!
!   pivot_tol_for_basis. The threshold pivot used by the matrix 
!    factorization when finding the basis.
!    See the documentation for ULS for details
!
!   zero_pivot. Any pivots smaller than zero_pivot in absolute value will 
!    be regarded to be zero when attempting to detect linearly dependent 
!    constraints
!
!   inner_fraction_opt. a search direction which gives at least 
!    inner_fraction_opt times the optimal model decrease will be found
!
!   inner_stop_relative and inner_stop_absolute. The search direction is
!    considered as an acceptable approximation to the minimizer of the
!    model if the gradient of the model in the preconditioning(inverse) 
!    norm is less than 
!     max( inner_stop_relative * initial preconditioning(inverse)
!                                 gradient norm, inner_stop_absolute )
!
!   max_infeasibility_relative and max_infeasibility_absolute. If the 
!     constraints are believed to be rank defficient and the residual
!     at a "typical" feasiblke point is larger than
!      max( max_infeasibility_relative * norm A, max_infeasibility_absolute )
!     the problem will be marked as infeasible
!
!   radius. An upper bound on the permitted step
!
!  LOGICAL control parameters:
!
!   find_basis_by_transpose. If true, implicit factorization preconditioners
!    will be based on a basis of A found by examining A's transpose
!
!   remove_dependencies. If true, the equality constraints will be preprocessed
!    to remove any linear dependencies
!
!   space_critical. If true, every effort will be made to use as little
!    space as possible. This may result in longer computation times
!
!   deallocate_error_fatal. If true, any array/pointer deallocation error
!     will terminate execution. Otherwise, computation will continue
!
!  CHARACTER control parameters:
!
!  prefix (len=30). All output lines will be prefixed by 
!    %prefix(2:LEN(TRIM(%prefix))-1)
!   where %prefix contains the required string enclosed in 
!   quotes, e.g. "string" or 'string'
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      TYPE ( EQP_data_type ), INTENT( OUT ) :: data
      TYPE ( EQP_control_type ), INTENT( OUT ) :: control        

!  Set control parameters

!  Integer parameters

      control%error  = 6
      control%out  = 6
      control%print_level = 0
      control%factorization = 0
      control%max_col = 35
      control%indmin = 1000
      control%valmin = 1000
      control%len_ulsmin = 1000
      control%itref_max = 1
!     control%cg_maxit = - 1
      control%cg_maxit = 200
      control%preconditioner = 0
      control%semi_bandwidth = 5
      control%new_a = 2
      control%new_h = 2

!  Real parameters

      control%radius = SQRT( point1 * HUGE( one ) )
      control%pivot_tol = 0.01_wp
!     control%pivot_tol = epsmch ** 0.75
      control%pivot_tol_for_basis = half
      control%zero_pivot = epsmch ** 0.75
      control%inner_fraction_opt = point1
!     control%inner_stop_relative = zero
      control%inner_stop_relative = point01
      control%inner_stop_absolute = SQRT( epsmch )
      control%max_infeasibility_relative = epsmch ** 0.75
      control%max_infeasibility_absolute = epsmch ** 0.75
      control%inner_stop_inter = 0.01_wp

!  Logical parameters

      control%remove_dependencies = .TRUE.
      control%find_basis_by_transpose = .TRUE.
      control%space_critical = .FALSE.
      control%deallocate_error_fatal = .FALSE.

!  Character parameters

      control%prefix = '""                            '

!  Ensure that the private data arrays have the correct initial status

      CALL SBLS_initialize( data%SBLS_data, control%SBLS_control )
      CALL GLTR_initialize( data%GLTR_data, control%GLTR_control )

!  Reset GLTR and SBLS data for this package

      control%GLTR_control%stop_relative = control%inner_stop_relative
      control%GLTR_control%stop_absolute = control%inner_stop_absolute
      control%GLTR_control%fraction_opt = control%inner_fraction_opt
      control%GLTR_control%unitm = .FALSE.
      control%GLTR_control%error = control%error
      control%GLTR_control%out = control%out
      control%GLTR_control%print_level = control%print_level - 1
      control%GLTR_control%itmax = control%cg_maxit
      control%GLTR_control%boundary = .FALSE.
      control%GLTR_control%space_critical = control%space_critical
      control%GLTR_control%deallocate_error_fatal                              &
        = control%deallocate_error_fatal
      control%GLTR_control%prefix = '" - GLTR:"                     '

      control%SBLS_control%error = control%error
      control%SBLS_control%out = control%out
      control%SBLS_control%print_level = control%print_level - 1
      control%SBLS_control%indmin = control%indmin
      control%SBLS_control%valmin = control%valmin
      control%SBLS_control%len_ulsmin = control%len_ulsmin
      control%SBLS_control%itref_max = control%itref_max
      control%SBLS_control%factorization = control%factorization
      control%SBLS_control%preconditioner = control%preconditioner
      control%SBLS_control%semi_bandwidth = control%semi_bandwidth
      control%SBLS_control%new_a = control%new_a
      control%SBLS_control%new_h = control%new_h
      control%SBLS_control%max_col = control%max_col
      control%SBLS_control%pivot_tol = control%pivot_tol
      control%SBLS_control%pivot_tol_for_basis = control%pivot_tol_for_basis
      control%SBLS_control%zero_pivot = control%zero_pivot
      control%SBLS_control%remove_dependencies = control%remove_dependencies
      control%SBLS_control%find_basis_by_transpose =                           &
        control%find_basis_by_transpose
      control%SBLS_control%deallocate_error_fatal =                            &
        control%deallocate_error_fatal
      control%SBLS_control%prefix = '" - SBLS:"                     '

      data%new_c = .TRUE.
      RETURN  

!  End of EQP_initialize

      END SUBROUTINE EQP_initialize

!-*-*-*-*-   E Q P _ R E A D _ S P E C F I L E  S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE EQP_read_specfile( control, device, alt_specname )

!  Reads the content of a specification file, and performs the assignment of 
!  values associated with given keywords to the corresponding control parameters

!  The defauly values as given by EQP_initialize could (roughly) 
!  have been set as:

!  BEGIN EQP SPECIFICATIONS (DEFAULT)
!   error-printout-device                           6
!   printout-device                                 6
!   print-level                                     0
!   initial-workspace-for-unsymmetric-solver        1000
!   initial-integer-workspace                       1000
!   initial-real-workspace                          1000
!   preconditioner-used                             0
!   semi-bandwidth-for-band-preconditioner          5
!   factorization-used                              0
!   maximum-column-nonzeros-in-schur-complement     35
!   maximum-refinements                             1
!   maximum-number-of-cg-iterations                 200
!   truat-region-radius                             1.0D+19
!   pivot-tolerance-used                            1.0D-12
!   pivot-tolerance-used-for-basis                  0.5
!   zero-pivot-tolerance                            1.0D-12
!   inner-iteration-fraction-optimality-required    0.1
!   inner-iteration-relative-accuracy-required      0.01
!   inner-iteration-absolute-accuracy-required      1.0E-8
!   inner-iteration-intermediate-accuracy-required  1.0D-2
!   max-relative-infeasibility-allowed              1.0E-12
!   max-absolute-infeasibility-allowed              1.0E-12
!   find-basis-by-transpose                         T
!   remove-linear-dependencies                      T
!   space-critical                                  F
!   deallocate-error-fatal                          F
!   output-line-prefix                              ""
!  END EQP SPECIFICATIONS

!  Dummy arguments

      TYPE ( EQP_control_type ), INTENT( INOUT ) :: control        
      INTEGER, INTENT( IN ) :: device
      CHARACTER( LEN = 16 ), OPTIONAL :: alt_specname

!  Programming: Nick Gould and Ph. Toint, January 2002.

!  Local variables

      INTEGER, PARAMETER :: lspec = 35
      CHARACTER( LEN = 16 ), PARAMETER :: specname = 'EQP             '
      TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec

!  Define the keywords

     spec%keyword = ''

!  Integer key-words

      spec(  1 )%keyword = 'error-printout-device'
      spec(  2 )%keyword = 'printout-device'
      spec(  3 )%keyword = 'print-level' 
      spec(  7 )%keyword = 'factorization-used'
      spec(  8 )%keyword = 'maximum-column-nonzeros-in-schur-complement'
      spec(  4 )%keyword = 'initial-workspace-for-unsymmetric-solver'
      spec(  9 )%keyword = 'initial-integer-workspace'
      spec( 10 )%keyword = 'initial-real-workspace'
      spec( 11 )%keyword = 'maximum-refinements'
      spec( 13 )%keyword = 'maximum-number-of-cg-iterations'
      spec( 14 )%keyword = 'preconditioner-used'
      spec( 15 )%keyword = 'semi-bandwidth-for-band-preconditioner'

!  Real key-words

      spec( 19 )%keyword = 'trust-region-radius'
      spec( 23 )%keyword = 'max-relative-infeasibility-allowed'
      spec( 24 )%keyword = 'max-absolute-infeasibility-allowed'
      spec( 26 )%keyword = 'pivot-tolerance-used'
      spec( 27 )%keyword = 'pivot-tolerance-used-for-basis'
      spec( 28 )%keyword = 'zero-pivot-tolerance'
      spec( 29 )%keyword = 'inner-iteration-intermediate-accuracy-required'
      spec( 30 )%keyword = 'inner-iteration-fraction-optimality-required'
      spec( 31 )%keyword = 'inner-iteration-relative-accuracy-required'
      spec( 32 )%keyword = 'inner-iteration-absolute-accuracy-required'

!  Logical key-words

      spec( 33 )%keyword = 'find-basis-by-transpose'
      spec( 34 )%keyword = 'remove-linear-dependencies'
      spec( 16 )%keyword = 'space-critical'
      spec( 35 )%keyword = 'deallocate-error-fatal'

!  Character key-words

!     spec( 36 )%keyword = 'output-line-prefix'

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
      CALL SPECFILE_assign_value( spec( 7 ), control%factorization,            &
                                    control%error )
      CALL SPECFILE_assign_value( spec( 8 ), control%max_col,                  &
                                    control%error )
      CALL SPECFILE_assign_value( spec( 4 ), control%len_ulsmin,               &
                                    control%error )
      CALL SPECFILE_assign_value( spec( 9 ), control%indmin,                   &
                                    control%error )
      CALL SPECFILE_assign_value( spec( 10 ), control%valmin,                  &
                                    control%error )
      CALL SPECFILE_assign_value( spec( 11 ), control%itref_max,               &
                                    control%error )
      CALL SPECFILE_assign_value( spec( 13 ), control%cg_maxit,                &
                                    control%error )
      CALL SPECFILE_assign_value( spec( 14 ), control%preconditioner,          &
                                    control%error )
      CALL SPECFILE_assign_value( spec( 15 ), control%semi_bandwidth,          &
                                    control%error )

!  Set real values

      CALL SPECFILE_assign_value( spec( 19 ), control%radius,                  &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 23 ),                                  &
                                  control%max_infeasibility_relative,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 24 ),                                  &
                                  control%max_infeasibility_absolute,          &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 26 ), control%pivot_tol,               &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 27 ),                                  &
                                  control%pivot_tol_for_basis,                 &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 28 ), control%zero_pivot,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 29 ), control%inner_stop_inter,        &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 30 ), control%inner_fraction_opt,      &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 31 ), control%inner_stop_relative,     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 32 ), control%inner_stop_absolute,     &
                                  control%error )

!  Set logical values

      CALL SPECFILE_assign_value( spec( 33 ),                                  &
                                  control%find_basis_by_transpose,             &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 34 ), control%remove_dependencies,     &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 16 ),                                  &
                                  control%space_critical,                      &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 35 ),                                  &
                                  control%deallocate_error_fatal,              &
                                  control%error )

!  Set charcter values

!     CALL SPECFILE_assign_string( spec( 36 ), control%prefix,                 &
!                                  control%error )

!  Reset GLTR and SBLS data for this package

      control%GLTR_control%stop_relative = control%inner_stop_relative
      control%GLTR_control%stop_absolute = control%inner_stop_absolute
      control%GLTR_control%error = control%error
      control%GLTR_control%out = control%out
      control%GLTR_control%print_level = control%print_level - 1
      control%GLTR_control%itmax = control%cg_maxit

      control%SBLS_control%error = control%error
      control%SBLS_control%out = control%out
      control%SBLS_control%print_level = control%print_level - 1
      control%SBLS_control%indmin = control%indmin
      control%SBLS_control%valmin = control%valmin
      control%SBLS_control%itref_max = control%itref_max
      control%SBLS_control%factorization = control%factorization
      control%SBLS_control%preconditioner = control%preconditioner
      control%SBLS_control%semi_bandwidth = control%semi_bandwidth
      control%SBLS_control%new_a = control%new_a
      control%SBLS_control%new_h = control%new_h
      control%SBLS_control%max_col = control%max_col
      control%SBLS_control%pivot_tol = control%pivot_tol
      control%SBLS_control%pivot_tol_for_basis = control%pivot_tol_for_basis
      control%SBLS_control%zero_pivot = control%zero_pivot
      control%SBLS_control%remove_dependencies = control%remove_dependencies
      control%SBLS_control%find_basis_by_transpose =                           &
        control%find_basis_by_transpose
      control%SBLS_control%space_critical = control%space_critical
      control%SBLS_control%deallocate_error_fatal =                            &
        control%deallocate_error_fatal

!  Read the controls for the preconditioner and iterative solver

      CALL SBLS_read_specfile( control%SBLS_control, device )
      CALL GLTR_read_specfile( control%GLTR_control, device )

      RETURN

      END SUBROUTINE EQP_read_specfile

!-*-*-*-*-*-*-*-*-   E Q P _ S O L V E  S U B R O U T I N E   -*-*-*-*-*-*-*-*-

      SUBROUTINE EQP_solve( prob, data, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Solve the quadratic program
!
!     minimize     q(x) = 1/2 x(T) H x + g(T) x + f
!
!     subject to    A x = c
!
!  where x is a vector of n components ( x_1, .... , x_n ), const is a constant
!  g is an n-vector, H is a symmetric matrix and A is an m by n matrix, 
!  using a projected preconditioned conjugate-gradient method.
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
!    to be solved since the last call to EQP_initialize, and .FALSE. if
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
!   %X is a REAL array of length %n, which must be set by the user
!    to an estimate of the solution x. On successful exit, it will contain
!    the required solution.
!
!   %C is a REAL array of length %m, which must be set by the user
!    to the values of the array c of constant terms of the constraints
!    Ax + c.= 0
!   
!  data is a structure of type EQP_data_type which holds private internal data
!
!  control is a structure of type EQP_control_type that controls the 
!   execution of the subroutine and must be set by the user. Default values for
!   the elements may be set by a call to EQP_initialize. See EQP_initialize 
!   for details
!
!  inform is a structure of type EQP_inform_type that provides 
!    information on exit from EQP_solve. The component status 
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
!          prob%n     >=  1
!          prob%m     >=  0
!          prob%A%type in { 'DENSE', 'SPARSE_BY_ROWS', 'COORDINATE' }
!          prob%H%type in { 'DENSE', 'SPARSE_BY_ROWS', 'COORDINATE', 'DIAGONAL' }
!       has been violated.
!
!    -5 the constraints are likely inconsistent
!
!    -9 an error has occured in SILS_analyse; the status as returned by 
!       AINFO%FLAG is given in the component sils_analyse_status
!
!   -10 an error has occured in SILS_factorize; the status as returned by 
!       FINFO%FLAG is given in the component sils_factorize_status
!
!   -11 an error has occured in SILS_solve; the status as returned by 
!       SINFO%FLAG is given in the component sils_solve_status
!
!   -12 an error has occured in ULS_analyse; the status as returned by 
!       AINFO%FLAG is given in the component uls_analyse_status
!
!   -14 an error has occured in ULS_solve; the status as returned by 
!       SINFO%FLAG is given in the component uls_solve_status
!
!   -15 the computed precondition is insufficient. Try another
!
!   -16 the residuals are large; the factorization may be unsatisfactory
!
!  On exit from EQP_solve, other components of inform give the 
!  following:
!
!     alloc_status = the status of the last attempted allocation/deallocation 
!     bad_alloc = the name of the last array for which (de)allocation failed
!     cg_iter = the total number of conjugate gradient iterations required.
!     factorization_integer = the total integer workspace required for the 
!       factorization.
!     factorization_real = the total real workspace required for the 
!       factorization.
!     obj = the value of the objective function at the best estimate of the 
!       solution determined by EQP_solve.
!     time%total = the total time spent in the package.
!     time%factorize = the time spent factorizing the required matrices.
!     time%solve = the time spent computing the search direction.
!     SBLS_inform = inform components from SBLS
!     GLTR_inform = inform components from GLTR
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      TYPE ( EQP_data_type ), INTENT( INOUT ) :: data
      TYPE ( EQP_control_type ), INTENT( INOUT ) :: control
      TYPE ( EQP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: i, j
      REAL :: time_end

!  prefix for all output 

      CHARACTER ( LEN = LEN( TRIM( control%prefix ) ) - 2 ) :: prefix
      prefix = control%prefix( 2 : LEN( TRIM( control%prefix ) ) - 1 )

!  Set initial values for inform 

      inform%status = GALAHAD_ok 
      inform%alloc_status = 0 ; inform%bad_alloc = ''
      inform%cg_iter = - 1 ; inform%cg_iter_inter = - 1
      inform%obj = zero
      inform%time%factorize = 0.0 ; inform%time%solve = 0.0
      inform%time%solve_inter = 0.0
      CALL CPU_TIME( inform%time%total )

      IF ( prob%n <= 0 .OR. prob%m < 0 .OR.                                    &
           .NOT. QPT_keyword_H( prob%H%type ) .OR.                             &
           .NOT. QPT_keyword_A( prob%A%type ) ) THEN
        inform%status = GALAHAD_error_restrictions
        IF ( control%error > 0 .AND. control%print_level > 0 )                 &
          WRITE( control%error,                                                &
            "( ' ', /, A, ' **  Error return ', I0,' from EQP ' )" )           &
            prefix, inform%status 
        RETURN
      END IF 

!  If required, write out problem 

      IF ( control%out > 0 .AND. control%print_level >= 20 ) THEN
        WRITE( control%out, "( ' n, m = ', 2I8 )" ) prob%n, prob%m
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
        WRITE( control%out, "( ' C = ', /, ( 5ES12.4 ) )" )                    &
          prob%C( : prob%m )
      END IF

!  Call the solver

      CALL EQP_solve_main( prob%n, prob%m, prob%H, prob%G, prob%f,             &
                           prob%A, prob%C, prob%q, prob%X, prob%Y,             &
                           data, control, inform )

      CALL CPU_TIME( time_end )
      inform%time%total = time_end - inform%time%total
      RETURN

!  End of EQP_solve

      END SUBROUTINE EQP_solve

!-*-*-*-*-   E Q P _ S O L V E _ M A I N   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE EQP_solve_main( n, m, H, G, f, A, C, q, X, Y,                 &
                                 data, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Solve the quadratic program
!
!     minimize     q(x) = 1/2 x(T) H x + g(T) x + f
!
!     subject to    (Ax)_i + c_i = 0 , i = 1, .... , m,
!
!  where x is a vector of n components ( x_1, .... , x_n ), f is a constant,
!  g is an n-vector, H is a symmetric matrix, A is an m by n matrix, and
!  c is an m-vector using a projected conjugate gradient method.
!  The subroutine is particularly appropriate when A and H are sparse
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( IN ) :: n, m
      REAL ( KIND = wp ), INTENT( IN ) :: f
      REAL ( KIND = wp ), INTENT( OUT ) :: q
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: G
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: C
      TYPE ( SMT_type ), INTENT( IN ) :: A
      TYPE ( SMT_type ), INTENT( IN ) :: H
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: Y
      TYPE ( EQP_data_type ), INTENT( INOUT ) :: data
      TYPE ( EQP_control_type ), INTENT( INOUT ) :: control
      TYPE ( EQP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: out, i, j, l, a_ne
      LOGICAL :: null_space_prec, explicit_prec, rank_def, solve_inter
      LOGICAL :: printi, printt, printw
      REAL :: time_end
      REAL ( KIND = wp ) :: radius, q_save, model, val, maxres
      REAL ( KIND = wp ) :: pgnorm, stop_inter
      CHARACTER ( LEN = 80 ) :: array_name

!  prefix for all output 

      CHARACTER ( LEN = LEN( TRIM( control%prefix ) ) - 2 ) :: prefix
      prefix = control%prefix( 2 : LEN( TRIM( control%prefix ) ) - 1 )

!  ===========================
!  Control the output printing
!  ===========================

      out = control%out

!  Basic single line of output per iteration

      printi = out > 0 .AND. control%print_level >= 1 

!  As per printi, but with additional timings for various operations

      printt = out > 0 .AND. control%print_level >= 2 

!  As per printt, but with checking of residuals, etc, and also with an 
!  indication of where in the code we are

      printw = out > 0 .AND. control%print_level >= 4

      inform%GLTR_inform%status = 1 
      inform%GLTR_inform%negative_curvature = .TRUE.

      IF ( SMT_get( A%type ) == 'DENSE' ) THEN
        a_ne = m * n
      ELSE IF ( SMT_get( A%type ) == 'SPARSE_BY_ROWS' ) THEN
        a_ne = A%ptr( m + 1 ) - 1
      ELSE
        a_ne = A%ne 
      END IF

!  Set C appropriately

      IF ( data%new_c ) THEN
        data%new_c = .FALSE.
        control%SBLS_control%new_c = 2
        data%C%ne = 0

        array_name = 'eqp: data%C%row'
        CALL SPACE_resize_array( data%C%ne, data%C%row, inform%status,         &
           inform%alloc_status, array_name = array_name,                       &
           deallocate_error_fatal = control%deallocate_error_fatal,            &
           exact_size = control%space_critical,                                &
           bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= GALAHAD_ok ) RETURN

        array_name = 'eqp: data%C%col'
        CALL SPACE_resize_array( data%C%ne, data%C%col, inform%status,         &
           inform%alloc_status, array_name = array_name,                       &
           deallocate_error_fatal = control%deallocate_error_fatal,            &
           exact_size = control%space_critical,                                &
           bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= GALAHAD_ok ) RETURN

        array_name = 'eqp: data%C%val'
        CALL SPACE_resize_array( data%C%ne, data%C%val, inform%status,         &
           inform%alloc_status, array_name = array_name,                       &
           deallocate_error_fatal = control%deallocate_error_fatal,            &
           exact_size = control%space_critical,                                &
           bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= GALAHAD_ok ) RETURN

        array_name = 'eqp: data%C%type'
        CALL SPACE_dealloc_array( data%C%type, inform%status,                  &
           inform%alloc_status, array_name = array_name, out = control%error )
        CALL SMT_put( data%C%type, 'COORDINATE', inform%alloc_status )
      ELSE
        control%SBLS_control%new_c = 0
      END IF

!  ------------------------------------------
!   1. Form and factorize the preconditioner
!  ------------------------------------------

      CALL CPU_TIME( inform%time%factorize )

      control%SBLS_control%error = control%error
      control%SBLS_control%out = control%out
      control%SBLS_control%print_level = control%print_level - 1
      control%SBLS_control%indmin = control%indmin
      control%SBLS_control%valmin = control%valmin
      control%SBLS_control%len_ulsmin = control%len_ulsmin
      control%SBLS_control%itref_max = control%itref_max
      control%SBLS_control%factorization = control%factorization
      control%SBLS_control%preconditioner = control%preconditioner
      control%SBLS_control%semi_bandwidth = control%semi_bandwidth
      control%SBLS_control%new_a = control%new_a
      control%SBLS_control%new_h = control%new_h
      control%SBLS_control%max_col = control%max_col
      control%SBLS_control%pivot_tol = control%pivot_tol
      control%SBLS_control%pivot_tol_for_basis = control%pivot_tol_for_basis
      control%SBLS_control%zero_pivot = control%zero_pivot
      control%SBLS_control%remove_dependencies = control%remove_dependencies
      control%SBLS_control%find_basis_by_transpose =                           &
        control%find_basis_by_transpose
      control%SBLS_control%deallocate_error_fatal =                            &
        control%deallocate_error_fatal
      control%SBLS_control%prefix = '" - SBLS:"                     '

      CALL SBLS_form_and_factorize( n, m, H, A, data%C, data%SBLS_data,        &
                                    control%SBLS_control, inform%SBLS_inform )

      null_space_prec = inform%SBLS_inform%factorization == 3
      explicit_prec = inform%SBLS_inform%preconditioner > 0
      rank_def = inform%SBLS_inform%rank_def

      CALL CPU_TIME( time_end )
      inform%time%factorize = time_end - inform%time%factorize

      IF ( inform%SBLS_inform%status < 0 ) THEN
        inform%status = inform%SBLS_inform%status
        RETURN
      END IF

      IF ( printt ) WRITE( out,                                                &
         "(  A, ' on exit from SBLS: status = ', I0, ', time = ', F0.2 )" )    &
           prefix, inform%SBLS_inform%status, inform%time%factorize

!  ------------------
!   2. Solve the EQP
!  ------------------

      CALL CPU_TIME( inform%time%solve )

      array_name = 'eqp: data%VECTOR'
      CALL SPACE_resize_array( n + m, data%VECTOR, inform%status,              &
         inform%alloc_status, array_name = array_name,                         &
         deallocate_error_fatal = control%deallocate_error_fatal,              &
         exact_size = control%space_critical,                                  &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%RESID'
      CALL SPACE_resize_array( m, data%RESID, inform%status,                   &
         inform%alloc_status, array_name = array_name,                         &
         deallocate_error_fatal = control%deallocate_error_fatal,              &
         exact_size = control%space_critical,                                  &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) RETURN

!  ---------------------------
!   2a. Find a feasible point
!  ---------------------------

!   Find a suitable feasible point x_f (stored in X)

!     data%VECTOR( : n  ) = - G
      data%VECTOR( : n  ) = zero
      data%VECTOR( n + 1 : n + m ) = - C( : m )

      IF ( null_space_prec ) THEN
        CALL SBLS_solve_null_space( data%SBLS_data%nfactors,                   &
                                    control%SBLS_control, inform%SBLS_inform,  &
                                    data%VECTOR )
      ELSE IF ( explicit_prec ) THEN
        CALL SBLS_solve_explicit( n, m, A, data%C, data%SBLS_data%efactors,    &
                                  control%SBLS_control, inform%SBLS_inform,    &
                                  data%VECTOR )
      ELSE
        CALL SBLS_basis_solve( data%SBLS_data%ifactors, control%SBLS_control,  &
                               inform%SBLS_inform, data%VECTOR )
      END IF

      IF ( inform%SBLS_inform%status < 0 ) THEN
        inform%status = inform%SBLS_inform%status
        GO TO 900
      END IF

      X = data%VECTOR( : n  )

!  Compute the constraint residuals

      data%RESID( : m ) = C( : m )
      SELECT CASE ( SMT_get( A%type ) )
      CASE ( 'DENSE' ) 
        l = 0
        DO i = 1, m
          data%RESID( i )                                                     &
            = data%RESID( i ) + DOT_PRODUCT( A%val( l + 1 : l + n ), X )
          l = l + n
        END DO
      CASE ( 'SPARSE_BY_ROWS' )
        DO i = 1, m
          DO l = A%ptr( i ), A%ptr( i + 1 ) - 1
            data%RESID( i ) = data%RESID( i ) + A%val( l ) * X( A%col( l ) )
          END DO
        END DO
      CASE ( 'COORDINATE' )
        DO l = 1, A%ne
          i = A%row( l )
          data%RESID( i ) = data%RESID( i ) + A%val( l ) * X( A%col( l ) )
        END DO
      END SELECT

!  Check to see if the residuals of potentially inconsistemt constraints
!  are satisfactory

      maxres = MAX( MAXVAL( ABS( A%val( : A_ne ) ) ),                         &
                    MAXVAL( ABS( C( : m ) ) ), MAXVAL( ABS( X( : n ) ) ) )
      maxres = MAX( control%max_infeasibility_relative * maxres,              &
                    control%max_infeasibility_absolute )

      IF ( rank_def ) THEN
        IF ( MAXVAL( ABS( data%RESID ) ) <= maxres ) THEN
          IF ( printw ) WRITE( out,                                           &
            "( A, ' residual acceptably small, consistent constraints' )" )   &
              prefix
        ELSE
          IF ( printi ) WRITE( out,                                           &
            "( A, ' residual too large, constraints likely inconsistent' )" ) &
              prefix
          inform%status = GALAHAD_error_primal_infeasible ; GO TO 900
        END IF
      ELSE
        IF ( MAXVAL( ABS( data%RESID ) ) > maxres ) THEN
          IF ( printi ) WRITE( out,                                           &
            "( A, ' residual too large, factorization likely inaccurate' )" ) &
              prefix
          inform%status = GALAHAD_error_ill_conditioned ; GO TO 900
        END IF
      END IF

!  Compute the function and gradient values at x_f

      array_name = 'eqp: data%G_f'
      CALL SPACE_resize_array( n, data%G_f, inform%status,                     &
         inform%alloc_status, array_name = array_name,                         &
         deallocate_error_fatal = control%deallocate_error_fatal,              &
         exact_size = control%space_critical,                                  &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) RETURN

      SELECT CASE ( SMT_get( H%type ) )
      CASE ( 'DIAGONAL' ) 
        DO i = 1, n
          data%G_f( i ) = H%val( i ) * X( i )
        END DO
      CASE ( 'DENSE' ) 
        data%G_f = zero
        l = 0
        DO i = 1, n
          DO j = 1, i
            l = l + 1 ; val = H%val( l )
            data%G_f( i ) = data%G_f( i ) + val * X( j )
            IF ( i /= j ) data%G_f( j ) = data%G_f( j ) + val * X( i )
          END DO
        END DO
      CASE ( 'SPARSE_BY_ROWS' )
        data%G_f = zero
        DO i = 1, n
          DO l = H%ptr( i ), H%ptr( i + 1 ) - 1
            j = H%col( l ) ; val = H%val( l )
            data%G_f( i ) = data%G_f( i ) + val * X( j )
            IF ( i /= j ) data%G_f( j ) = data%G_f( j ) + val * X( i )
          END DO
        END DO
      CASE ( 'COORDINATE' )
        data%G_f = zero
        DO l = 1, H%ne
          i = H%row( l ) ; j = H%col( l ) ; val = H%val( l )
          data%G_f( i ) = data%G_f( i ) + val * X( j )
          IF ( i /= j ) data%G_f( j ) = data%G_f( j ) + val * X( i )
        END DO
      END SELECT

      q_save = f + DOT_PRODUCT( G, X ) + half * DOT_PRODUCT( data%G_f, X )
      data%G_f = G + data%G_f
      control%GLTR_control%f_0 = q_save

      IF ( printi ) WRITE( out,                                                &
         "( /, A, ' constraint & objective (feasibility phase) =',             &
        &  ES11.4, ',', ES12.4 )" ) prefix, MAXVAL( ABS( data%RESID ) ), q_save

!  --------------------------------------------------
!   2b. From the feasible point, use GLTR to minimize 
!       the objective in the null-space of A
!  --------------------------------------------------

!   Compute the correction s from x_f (stored in S)

      array_name = 'eqp: data%R'
      CALL SPACE_resize_array( n + m, data%R, inform%status,                   &
         inform%alloc_status, array_name = array_name,                         &
         deallocate_error_fatal = control%deallocate_error_fatal,              &
         exact_size = control%space_critical,                                  &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%S'
      CALL SPACE_resize_array( n, data%S, inform%status,                       &
         inform%alloc_status, array_name = array_name,                         &
         deallocate_error_fatal = control%deallocate_error_fatal,              &
         exact_size = control%space_critical,                                  &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%WORK'
      CALL SPACE_resize_array( n, data%WORK, inform%status,                    &
         inform%alloc_status, array_name = array_name,                         &
         deallocate_error_fatal = control%deallocate_error_fatal,              &
         exact_size = control%space_critical,                                  &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) RETURN

!  Set initial data
     
      IF ( control%radius > zero ) THEN
        radius = control%radius
      ELSE
        radius = SQRT( point1 * HUGE( one ) )
      END IF 
      q = q_save
      data%R( : n ) = data%G_f
      inform%GLTR_inform%status = 1
      solve_inter = .FALSE.
      inform%cg_iter = 0 ; inform%cg_iter_inter = 0

      control%GLTR_control%stop_relative = control%inner_stop_relative
      control%GLTR_control%stop_absolute = control%inner_stop_absolute
      control%GLTR_control%fraction_opt = control%inner_fraction_opt
      control%GLTR_control%rminvr_zero = hundred * epsmch ** 2
      control%GLTR_control%unitm = .FALSE.
      control%GLTR_control%error = control%error
      control%GLTR_control%out = control%out
      control%GLTR_control%print_level = control%print_level - 1
      control%GLTR_control%itmax = control%cg_maxit
      control%GLTR_control%boundary = .FALSE.
      control%GLTR_control%space_critical = control%space_critical
      control%GLTR_control%deallocate_error_fatal                              &
        = control%deallocate_error_fatal
      control%GLTR_control%prefix = '" - GLTR:"                     '

      DO
        CALL GLTR_solve( n, radius, model, data%S, data%R( : n ),              &
                         data%VECTOR( : n ), data%GLTR_data,                   &
                         control%GLTR_control, inform%GLTR_inform )

!  Check for error returns

!       WRITE(6,"( ' case ', i3  )" ) inform%GLTR_inform%status
        SELECT CASE( inform%GLTR_inform%status )

!  Successful return

        CASE ( GALAHAD_ok )
          EXIT

!  Warnings

        CASE ( GALAHAD_warning_on_boundary, GALAHAD_error_max_iterations )
          IF ( printt ) WRITE( out, "( /,                                      &
         &  A, ' Warning return from GLTR, status = ', I6 )" ) prefix,         &
              inform%GLTR_inform%status
          EXIT
          
!  Allocation errors

           CASE ( GALAHAD_error_allocate )
             inform%status = GALAHAD_error_allocate
             inform%alloc_status = inform%gltr_inform%alloc_status
             inform%bad_alloc = inform%gltr_inform%bad_alloc
             GO TO 900

!  Deallocation errors

           CASE ( GALAHAD_error_deallocate )
             inform%status = GALAHAD_error_deallocate
             inform%alloc_status = inform%gltr_inform%alloc_status
             inform%bad_alloc = inform%gltr_inform%bad_alloc
             GO TO 900

!  Error return

        CASE DEFAULT
          inform%status = inform%gltr_inform%status
          IF ( printt ) WRITE( out, "( /,                                      &
         &  A, ' Error return from GLTR, status = ', I6 )" ) prefix,           &
              inform%GLTR_inform%status
          EXIT

!  Find the preconditioned gradient

        CASE ( 2, 6 )
          IF ( printw ) WRITE( out,                                            &
             "( A, ' ............... precondition  ............... ' )" ) prefix
          data%VECTOR( n + 1 : n + m ) = zero

!     control%SBLS_control%out = 6
!     control%SBLS_control%print_level = 2

          control%SBLS_control%affine = .TRUE.
          CALL SBLS_solve( n, m, A, data%C, data%SBLS_data,                    &
             control%SBLS_control, inform%SBLS_inform, data%VECTOR )

          IF ( inform%SBLS_inform%status < 0 ) THEN
            inform%status = inform%SBLS_inform%status
            GO TO 900
          END IF

          IF ( inform%GLTR_inform%status == 2 ) THEN
            pgnorm = DOT_PRODUCT( data%R( : n ), data%VECTOR( : n ) )
            IF ( ABS( pgnorm ) < ten * EPSILON( one ) ) pgnorm = zero
            pgnorm = SIGN( SQRT( ABS( pgnorm ) ), pgnorm ) 
            IF ( inform%cg_iter == 0 )                                         &
              stop_inter = pgnorm * control%inner_stop_inter
          END IF

!  Compute the residuals

          IF ( printw ) THEN 
            data%RESID = zero
            SELECT CASE ( SMT_get( A%type ) )
            CASE ( 'DENSE' ) 
              l = 0
              DO i = 1, m
                data%RESID( i ) = data%RESID( i )                              &
                  + DOT_PRODUCT( A%val( l + 1 : l + n ), data%VECTOR )
                l = l + n
              END DO
            CASE ( 'SPARSE_BY_ROWS' )
              DO i = 1, m
                DO l = A%ptr( i ), A%ptr( i + 1 ) - 1
                  data%RESID( i )                                              &
                    = data%RESID( i ) + A%val( l ) * data%VECTOR( A%col( l ) )
                END DO
              END DO
            CASE ( 'COORDINATE' )
              DO l = 1, A%ne
                i = A%row( l )
                data%RESID( i )                                                &
                  = data%RESID( i ) + A%val( l ) * data%VECTOR( A%col( l ) )
              END DO
            END SELECT
            WRITE( out,                                                        &
              "( A, ' Constraint residual (optimality phase) ', ES10.4 )" )    &
                prefix, MAXVAL( ABS( data%RESID ) )
          END IF

!  Form the product of VECTOR with H

        CASE ( 3, 7 )

          IF ( inform%GLTR_inform%status == 3 ) THEN
            inform%cg_iter = inform%cg_iter + 1
            IF ( .NOT. solve_inter .AND. pgnorm <= stop_inter ) THEN
              inform%cg_iter_inter = inform%cg_iter
              CALL CPU_TIME( time_end )
              inform%time%solve_inter = time_end - inform%time%solve
              solve_inter = .TRUE.
            END IF
          END IF

          IF ( printw ) WRITE( out,                                            &
            "( A, ' ............ matrix-vector product ..........' )" ) prefix

          data%WORK( : n ) = data%VECTOR( : n )

          SELECT CASE ( SMT_get( H%type ) )
          CASE ( 'DIAGONAL' ) 
            DO i = 1, n
              data%VECTOR( i ) = H%val( i ) * data%WORK( i )
            END DO
          CASE ( 'DENSE' ) 
            data%VECTOR( : n ) = zero
            l = 0
            DO i = 1, n
              DO j = 1, i
                l = l + 1 ; val = H%val( l )
                data%VECTOR( i ) = data%VECTOR( i ) + val * data%WORK( j )
                IF ( i /= j ) data%VECTOR( j )                                 &
                  = data%VECTOR( j ) + val * data%WORK( i )
              END DO
            END DO
          CASE ( 'SPARSE_BY_ROWS' )
            data%VECTOR( : n ) = zero
            DO i = 1, n
              DO l = H%ptr( i ), H%ptr( i + 1 ) - 1
                j = H%col( l ) ; val = H%val( l )
                data%VECTOR( i ) = data%VECTOR( i ) + val * data%WORK( j )
                IF ( i /= j ) data%VECTOR( j )                                 &
                  = data%VECTOR( j ) + val * data%WORK( i )
              END DO
            END DO
          CASE ( 'COORDINATE' )
            data%VECTOR( : n ) = zero
            DO l = 1, H%ne
              i = H%row( l ) ; j = H%col( l ) ; val = H%val( l )
              data%VECTOR( i ) = data%VECTOR( i ) + val * data%WORK( j )
              IF ( i /= j ) data%VECTOR( j )                                   &
                = data%VECTOR( j ) + val * data%WORK( i )
            END DO
          END SELECT

!  Reform the initial residual

        CASE ( 5 )
          
          IF ( printw ) WRITE( out,                                            &
            "( A, ' ................. restarting ................ ' )" ) prefix

          q = q_save
          data%R( : n ) = data%G_f

        END SELECT

      END DO

!  Update the solution

      X = X + data%S

      IF ( .NOT. solve_inter ) THEN
        CALL CPU_TIME( time_end )
        inform%cg_iter_inter = inform%cg_iter
        inform%time%solve_inter = time_end - inform%time%solve
      END IF

!  Compute the residuals

      data%RESID( : m ) = C( : m )
      SELECT CASE ( SMT_get( A%type ) )
      CASE ( 'DENSE' ) 
        l = 0
        DO i = 1, m
          data%RESID( i )                                                      &
            = data%RESID( i ) + DOT_PRODUCT( A%val( l + 1 : l + n ), X )
          l = l + n
        END DO
      CASE ( 'SPARSE_BY_ROWS' )
        DO i = 1, m
          DO l = A%ptr( i ), A%ptr( i + 1 ) - 1
            data%RESID( i ) = data%RESID( i ) + A%val( l ) * X( A%col( l ) )
          END DO
        END DO
      CASE ( 'COORDINATE' )
        DO l = 1, A%ne
          i = A%row( l )
          data%RESID( i ) = data%RESID( i ) + A%val( l ) * X( A%col( l ) )
        END DO
      END SELECT

      inform%obj = model

      CALL CPU_TIME( time_end )
      inform%time%solve = time_end - inform%time%solve

      IF ( printt ) WRITE( out,                                                &
         "(  A, ' on exit from GLTR: status = ', I0, ', CG iterations = ', I0, &
        &   ', time = ', F0.2 )" ) prefix,                                     &
            inform%GLTR_inform%status, inform%cg_iter, inform%time%solve
      IF ( printi ) WRITE( out,                                                &
         "(  A, ' constraint & objective (optimality',                         &
        &   '  phase) =', ES11.4, ',', ES12.4 )" )                             &
            prefix, MAXVAL( ABS( data%RESID ) ), inform%obj

!  Compute the Lagrange multiplier estimates

      Y = data%VECTOR( n + 1 : n + m )

      IF ( printt ) THEN
        SELECT CASE( inform%SBLS_inform%preconditioner )
        CASE( 1 )
          WRITE( out, "( A, ' Preconditioner G = I' )" ) prefix
        CASE( 2 )
          WRITE( out, "( A, ' Preconditioner G = H' )" ) prefix
        CASE( 3 )
          WRITE( out, "( A, ' Preconditioner G = diag(H)' )" ) prefix
        CASE( 4 )
          WRITE( out, "( A, ' Preconditioner G = H_22' )" ) prefix
        CASE( 5 )
          WRITE( out, "( A, ' Preconditioner G = H_22 & H_21' )" ) prefix
        CASE( - 1 )
          WRITE( out, "( A, ' Preconditioner G_22 = I' )" ) prefix
        CASE( - 2 )
          WRITE( out, "( A, ' Preconditioner G_22 = H_22' )" ) prefix
        CASE( - 3 )
          WRITE( out, "( A, ' Preconditioner G_22 = H_22 and G_21 = H_21')")   &
            prefix
        END SELECT
      END IF

      RETURN
 
!  Error returns

  900 CONTINUE
      CALL CPU_TIME( time_end )
      inform%time%solve = time_end - inform%time%solve

      RETURN

!  End of EQP_solve_main

      END SUBROUTINE EQP_solve_main

!-*-*-*-*-*-*-   E Q P _ T E R M I N A T E   S U B R O U T I N E   -*-*-*-*-*-*

      SUBROUTINE EQP_terminate( data, control, inform )

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
!   data    see Subroutine EQP_initialize
!   control see Subroutine EQP_initialize
!   inform  see Subroutine EQP_solve

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( EQP_data_type ), INTENT( INOUT ) :: data
      TYPE ( EQP_control_type ), INTENT( IN ) :: control        
      TYPE ( EQP_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      CHARACTER ( LEN = 80 ) :: array_name

!  Deallocate all arrays allocated by SBLS and GLTR

      CALL SBLS_terminate( data%SBLS_data, control%SBLS_control,               &
                           inform%SBLS_inform )
      CALL GLTR_terminate( data%GLTR_data, control%GLTR_control,               &
                           inform%GLTR_inform )

!  Deallocate all remaining allocated arrays

      array_name = 'eqp: data%C%row'
      CALL SPACE_dealloc_array( data%C%row,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%C%col'
      CALL SPACE_dealloc_array( data%C%col,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%C%val'
      CALL SPACE_dealloc_array( data%C%val,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%C%type'
      CALL SPACE_dealloc_array( data%C%type,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%G_f'
      CALL SPACE_dealloc_array( data%G_f,                                      &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%R'
      CALL SPACE_dealloc_array( data%R,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%S'
      CALL SPACE_dealloc_array( data%S,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%VECTOR'
      CALL SPACE_dealloc_array( data%VECTOR,                                   &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%WORK'
      CALL SPACE_dealloc_array( data%WORK,                                     &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'eqp: data%RESID'
      CALL SPACE_dealloc_array( data%RESID,                                    &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

!  End of subroutine EQP_terminate

      END SUBROUTINE EQP_terminate

!  End of module EQP

   END MODULE GALAHAD_EQP_double
