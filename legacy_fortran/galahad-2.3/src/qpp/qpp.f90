! THIS VERSION: GALAHAD 2.2 - 23/04/2008 AT 08:00 GMT.

!-*-*-*-*-*-*-*-*-*-  G A L A H A D _ Q P P    M O D U L E  -*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released pre GALAHAD Version 1.0. July 29th 1999
!   update released with GALAHAD Version 2.0. February 16th 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_QPP_double

!     --------------------------------------------------------------
!     |                                                            |
!     | Preprocess the data for the quadratic program              |
!     |                                                            |
!     |    minimize     1/2 x(T) H x + g(T) x + f                  |
!     |                                                            |
!     |    subject to     c_l <=  A x  <= c_u                      |
!     |    and            x_l <=   x   <= x_u                      |
!     |                                                            |
!     | to make further manipulations more efficient               |
!     |                                                            |
!     | Optionally instead do this for the parametric problem      |
!     |                                                            |
!     |    minimize   1/2 x(T) H x + g(T) x + f + theta dg(T) x    |
!     |                                                            |
!     |    subject to c_l + theta dc_l <= A x <= c_u + theta dc_u  |
!     |    and        x_l + theta dx_l <=  x  <= x_u + theta dx_u  |
!     |                                                            |
!     | where theta is a scalar parameter                          |
!     |                                                            |
!     --------------------------------------------------------------

      USE GALAHAD_SYMBOLS
      USE GALAHAD_SMT_double
      USE GALAHAD_QPT_double
      USE GALAHAD_SORT_double,                                                 &
        ONLY: SORT_inplace_permute, SORT_inverse_permute, SORT_quicksort

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: SMT_put, SMT_get, QPT_problem_type, QPP_initialize,            &
                QPP_reorder, QPP_apply, QPP_restore, QPP_get_values,           &
                QPP_terminate

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

      TYPE, PUBLIC :: QPP_control_type
        INTEGER :: error
        REAL ( KIND = wp ) :: infinity
        LOGICAL :: treat_zero_bounds_as_general
      END TYPE

      TYPE, PUBLIC :: QPP_inform_type
        INTEGER :: status, alloc_status
      END TYPE

      TYPE, PUBLIC :: QPP_map_type
        INTEGER :: m, n, a_ne, h_ne, h_diag_end_fixed, a_type, h_type
        INTEGER :: m_reordered, n_reordered, a_ne_original, h_ne_original
        LOGICAL :: set, h_perm, a_perm
        CHARACTER, ALLOCATABLE, DIMENSION( : ) :: a_type_original
        CHARACTER, ALLOCATABLE, DIMENSION( : ) :: h_type_original
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: x_map
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: c_map
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IW
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: h_map_inverse
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: a_map_inverse
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ptr_a_fixed
      END TYPE

      TYPE, PUBLIC :: QPP_dims_type
        INTEGER :: nc = -1 , x_s = - 1, x_e = - 1, c_b = - 1, c_s = - 1 
        INTEGER :: c_e = - 1, y_s = - 1, y_i = - 1, y_e = - 1, v_e = - 1
        INTEGER :: x_free = - 1, x_l_start = - 1, x_l_end = - 1
        INTEGER :: x_u_start = - 1, x_u_end = - 1
        INTEGER :: c_equality = - 1, c_l_start = - 1, c_l_end = - 1
        INTEGER :: c_u_start = - 1, c_u_end = - 1
        INTEGER :: h_diag_end_free = - 1, h_diag_end_nonneg = - 1 
        INTEGER :: h_diag_end_lower = - 1, h_diag_end_range = - 1
        INTEGER :: h_diag_end_upper = - 1, h_diag_end_nonpos = - 1
        REAL ( KIND = wp ) :: f = HUGE( 1.0_wp )
      END TYPE

!----------------------
!   P a r a m e t e r s
!----------------------

      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: infinity = HUGE( one )
      INTEGER, PARAMETER :: GALAHAD_error_unordered = - 31
      INTEGER, PARAMETER :: GALAHAD_error_reformat = - 32
      INTEGER, PARAMETER :: GALAHAD_error_ah_unordered = - 33
      INTEGER, PARAMETER :: GALAHAD_error_y_unallocated = - 34
      INTEGER, PARAMETER :: GALAHAD_error_z_unallocated = - 35

   CONTAINS

!-*-*-*-*-*-   Q P P _ i n i t i a l i z e   S U B R O U T I N E  -*-*-*-*-*-

      SUBROUTINE QPP_initialize( map, control )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Default control data for QPP. This routine should be called 
!  before QPP_reorder
! 
!  --------------------------------------------------------------------
!
!  Arguments:
!
!  map      map arrays
!  control  a structure containing control information. Components are -
!
!  INTEGER control parameter:
!
!   error. Error and warning diagnostics occur on stream error 
!   
!  REAL control parameter:
!
!   infinity. Any bound larger or equal to infinity in abolute value 
!    will be considered to be infinite
!
!  LOGICAL control parameter:
!
!    treat_zero_bounds_as_general. If true, any problem bound with the value
!     zero will be treated as if it were a general value
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( QPP_map_type ), INTENT( OUT ) :: map
      TYPE ( QPP_control_type ), INTENT( OUT ) :: control        

!  Set control parameters
!  ======================

!  Integer parameter

      control%error = 6

!  Real parameter

      control%infinity = infinity

!  Logical parameter

      control%treat_zero_bounds_as_general = .FALSE.

      map%set = .FALSE.

      RETURN  

!  End of QPP_initialize

      END SUBROUTINE QPP_initialize

!-*-*-*-*-*-*-*-   Q P P _ r e o r d e r    S U B R O U T I N E  -*-*-*-*-*-*-

      SUBROUTINE QPP_reorder( map, control, inform, dims, prob,               &
                              get_x, get_y, get_z, parametric )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find a reordering of the data for the problem
!
!     minimize          q(x) = 1/2 x(T) H x + g(T) x 
!
!     subject to the bounds  x_l <=  x  <= x_u
!     and constraints        c_l <= A x <= c_u ,
!
!  where x is a vector of n components ( x_1, .... , x_n ), 
!  H is a symmetric matrix, A is an m by n matrix, 
!  and any of the bounds x_l, x_u, c_l, c_u may be infinite.
!
!  The reordered problem has the following properties:
!
!  * the variables are ordered so that their bounds appear in the order
!
!    free                      x
!    non-negativity      0  <= x
!    lower              x_l <= x
!    range              x_l <= x <= x_u
!    upper                     x <= x_u
!    non-positivity            x <=  0
!
!    Fixed variables will be removed. Within each category, the variables 
!    are further ordered so that those with non-zero diagonal Hessian 
!    entries occur before the remainder
!
!  * the constraints are ordered so that their bounds appear in the order
!
!    equality           c_l  = A x
!    lower              c_l <= A x
!    range              c_l <= A x <= c_u
!    upper                     A x <= c_u
!
!    Free constraints will be removed. 

!  * additional constraints may be added, bounds tightened,
!    to reduce the size of the feasible region if this is
!    possible
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!
!  dims is a structure of type QPP_dims_type, whose components hold SCALAR
!   information about the transformed problem on output. 
!
!  prob is a structure of type QPT_type, whose components hold the
!   details of the problem. The following components must be set

!   f is a REAL variable, which must be set by the user to the value of
!    the constant term f in the objective function. On exit, it may have
!    been changed to reflect variables which have been fixed.
!
!   n is an INTEGER variable, which must be set by the user to the
!    number of optimization parameters, n.  RESTRICTION: n >= 1
!                  
!   m is an INTEGER variable, which must be set by the user to the
!    number of general linear constraints, m.  RESTRICTION: m >= 0
!
!   H is a structure of type SMT_type used to hold the LOWER TRIANGULAR part 
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
!   G is a REAL array of length n, which must be set by
!    the user to the value of the gradient, g, of the linear term of the
!    quadratic objective function. The i-th component of G, i = 1, ....,
!    n, should contain the value of g_i.  
!    On exit, G will most likely have been reordered.
!   
!   DG is a REAL array of length n, which, if allocated, must
!    be set by the user to the values of the array dg of the parametric 
!    linear term of the quadratic objective function. 
!    On exit, present DG will most likely have been reordered.
!   
!   A is a structure of type SMT_type used to hold the matrix A. 
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
!   C is a REAL array of length m, which is used to store the values of 
!    A x. On exit, it will have been filled with appropriate values.
!
!   get_x is a LOGICAL variable. See X.
!
!   X is a REAL array of length n, which is used to store the values of 
!   the variables x. If the user has assigned values to X, get_x must be .FALSE.
!   on entry, and X filled with the values of x. In this case, on exit, 
!   X will most likely have been reordered, and any fixed variables moved
!   to their bounds. If the user does not wish to provide values for X, 
!   get_x must be .TRUE. on entry. On exit it will have been filled with 
!   appropriate values.
!
!   X_l, X_u are REAL arrays of length n, which must be set by the user
!    to the values of the arrays x_l and x_u of lower and upper bounds on x.
!    Any bound x_l_i or x_u_i larger than or equal to infinity in absolute value
!    will be regarded as being infinite (see the entry control%infinity).
!    Thus, an infinite lower bound may be specified by setting the appropriate 
!    component of X_l to a value smaller than -infinity, while an infinite 
!    upper bound can be specified by setting the appropriate element of X_u
!    to a value larger than infinity. 
!    On exit, X_l and X_u will most likely have been reordered.
!   
!   DX_l, DX_u are REAL arrays of length n, which, if allocated, must
!    be set by the user to the values of the arrays dx_l and dx_u of parametric 
!    lower and upper bounds on x. 
!    On exit, present DX_l and DX_u will most likely have been reordered.
!   
!   get_z is a LOGICAL variable. See Z
!
!   Z is a REAL array of length n, which are used to store the values
!    of the dual variables (Lagrange multipliers) corresponding to the simple
!    bound constraints x_l <= x and x <= x_u. If the 
!    user has assigned values to Z, get_z must be .FALSE. on entry.
!    In this case, on exit, Z will most likely have been reordered. 
!    If the user does not wish to provide values for Z, get_Z must be 
!    .TRUE.on entry. On exit it will have been filled with appropriate values.
!
!   C_l, C_u are REAL array of length m, which must be set by the user to 
!    the values of the arrays c_l and c_u of lower and upper bounds on A x.
!    Any bound c_l_i or c_u_i larger than or equal to infinity in absolute value
!    will be regarded as being infinite (see the entry control%infinity).
!    Thus, an infinite lower bound may be specified by setting the appropriate 
!    component of C_u to a value smaller than -infinity, while an infinite 
!    upper bound can be specified by setting the appropriate element of BU
!    to a value larger than infinity. 
!    On exit, C_l and C_u will most likely have been reordered.
!  
!   DC_l, DC_u are REAL arrays of length n, which, if allocated, must
!    be set by the user to the values of the arrays dc_l and dc_u of parametric 
!    lower and upper bounds on A x. 
!    On exit, present DC_l and DC_u will most likely have been reordered
!   
!   get_y is a LOGICAL variable. See Y.
!
!   Y is a REAL array of length m, which are used to store the values
!    of the Lagrange multipliers corresponding to the general bound constraints 
!    c_l <= A x and A x <= c_u. If the user has assigned values 
!    to Y, get_y must be .FALSE. on entry. In this case, on exit, 
!    Y will most likely have been reordered. If the user does not 
!    wish to provide values for Y, get_y must be .TRUE. on entry. 
!    On exit it will have been filled with appropriate values.
!
!  map is a structure of type QPP_data_type which holds internal 
!   mapping arrays and related data, and which need not be set by the user
!
!  control is a structure of type QPP_control_type that controls the
!   execution of the subroutine and must be set by the user. Default values for
!   the elements may be set by a call to QPP_initialize. See 
!   QPP_initialize for details
!
!  inform is a structure of type QPP_inform_type that provides 
!    information on exit from QPP_reorder. The component status 
!    has possible values:
!  
!     0 Normal termination
!
!    -1 An allocation error occured; the status is given in the component
!       alloc_status.
!
!    -2 A deallocation error occured; the status is given in the component
!       alloc_status.
!
!    -3 one of the restrictions 
!          prob%n     >=  1
!          prob%m     >=  0
!          prob%A%type in { 'DENSE', 'SPARSE_BY_ROWS', 'COORDINATE' }
!          prob%H%type in { 'DENSE', 'SPARSE_BY_ROWS', 'COORDINATE', 'DIAGONAL' }
!       has been violated
!
!    -5 The constraints are inconsistent
!
!   -23 an entry from the strict upper triangle of H has been input
!
!   -31 an attempt to use QPP_apply/QPP_restore is made prior to a 
!        successful call to QPP_reorder
!
!   -32 the storage format has changed without recalling QPP_reorder
!
!   -33 the array A/H have not been reordered, but the given option requires
!       them to have been
!
!   -34 Neither the array prob%Y nor the pair prob%Y_l and prob%Y_u have
!       been allocated
!
!   -35 Neither the array prob%Z nor the pair prob%Z_l and prob%Z_u have
!       been allocated
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( QPP_map_type ), INTENT( INOUT ) :: map
      TYPE ( QPP_dims_type ), INTENT( OUT ) :: dims
      LOGICAL, INTENT( IN ) :: get_x, get_y, get_z
      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      TYPE ( QPP_control_type ), INTENT( IN ) :: control        
      TYPE ( QPP_inform_type ), INTENT( OUT ) :: inform
      LOGICAL, OPTIONAL :: parametric

!  Local variables

      INTEGER :: i, j, k, l, ll
      INTEGER :: free, nonneg, lower, range, upper, nonpos, fixed, equality
      INTEGER :: a_free, a_lower, a_range, a_upper, a_equality
      INTEGER :: h_free, h_nonneg, h_lower, h_range, h_upper, h_nonpos, h_fixed
      INTEGER :: d_free, o_free, d_nonneg, o_nonneg, d_lower, o_lower, d_range
      INTEGER :: o_range, d_upper, o_upper, d_nonpos, o_nonpos, d_fixed, o_fixed
      REAL ( KIND = wp ) :: xl, xu, cl, cu
      LOGICAL :: reallocate, apy, apyl, apyu, apz, apzl, apzu
      CHARACTER ( LEN = 20 ) :: bad_alloc

!  Check that Y or Y_l/Y_u and Z or Z_l/Z_u has been allocated

      apy = ALLOCATED( prob%Y )
      apyl = ALLOCATED( prob%Y_l ) ; apyu = ALLOCATED( prob%Y_u )
      IF ( .NOT. ( apy .OR. ( apyl .AND. apyu ) ) ) THEN
        inform%status = GALAHAD_error_y_unallocated ; RETURN
      END IF

      apz = ALLOCATED( prob%Z )
      apzl = ALLOCATED( prob%Z_l ) ; apzu = ALLOCATED( prob%Z_u )
      IF ( .NOT. ( apz .OR. ( apzl .AND. apzu ) ) ) THEN
        inform%status = GALAHAD_error_z_unallocated ; RETURN
      END IF

!  Check input parameters

      IF ( prob%n <= 0 .OR. prob%m < 0 .OR.                                    &
           .NOT. QPT_keyword_A( prob%A%type ) ) THEN
        inform%status = GALAHAD_error_restrictions ; RETURN
      ELSE IF ( prob%Hessian_kind < 0 ) THEN
        IF ( .NOT. QPT_keyword_H( prob%H%type ) ) THEN
          inform%status = GALAHAD_error_restrictions ; RETURN
        END IF
      END IF

!  Store original problem dimensions

      map%n = prob%n
      map%m = prob%m
      IF ( prob%Hessian_kind < 0 ) THEN
        CALL QPT_put_H( map%h_type_original, SMT_get( prob%H%type ) )
        IF ( SMT_get( prob%H%type ) == 'COORDINATE' )                          &
           map%h_ne_original = prob%H%ne
      END IF
      CALL QPT_put_A( map%a_type_original, SMT_get( prob%A%type ) )
      IF ( SMT_get( prob%A%type ) == 'COORDINATE' )                            &
         map%a_ne_original = prob%A%ne

      IF ( prob%Hessian_kind < 0 ) THEN

!  See what storage scheme is used for H

        IF ( SMT_get( prob%H%type ) == 'DIAGONAL' ) THEN
          map%h_type = 3
          map%h_ne = prob%n
        ELSE IF ( SMT_get( prob%H%type ) == 'DENSE' ) THEN

!  dense

          map%h_type = 2
          map%h_ne = ( prob%n * ( prob%n + 1 ) ) / 2
        ELSE IF ( SMT_get( prob%H%type ) == 'SPARSE_BY_ROWS' ) THEN

!  row-wise, sparse

          map%h_type = 1
          map%h_ne = prob%H%ptr( prob%n + 1 ) - 1
          DO i = 1, map%n
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              IF ( prob%H%col( l ) > i ) THEN
                inform%status = GALAHAD_error_upper_entry ; RETURN
              END IF
            END DO
          END DO
        ELSE

! co-ordinate, sparse

          map%h_type = 0
          map%h_ne = prob%H%ne
        END IF

      END IF

!  Do the same for A

      IF ( SMT_get( prob%A%type ) == 'DENSE' ) THEN

!  dense

        map%a_type = 2
        map%a_ne = prob%n * prob%m
      ELSE IF ( SMT_get( prob%A%type ) == 'SPARSE_BY_ROWS' ) THEN

!  row-wise, sparse

        map%a_type = 1
        map%a_ne = prob%A%ptr( prob%m + 1 ) - 1
      ELSE

! co-ordinate, sparse

        map%a_type = 0
        map%a_ne = prob%A%ne
      END IF

!  Allocate workspace array

      reallocate = .TRUE.
      IF ( ALLOCATED( map%IW ) ) THEN
        IF ( SIZE( map%IW ) /= MAX( map%n, map%m ) ) THEN ; DEALLOCATE( map%IW )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( map%IW( MAX( map%n, map%m ) ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          bad_alloc = 'map%IW' ; GO TO 900
        END IF
      END IF

!  Check to see which variables have corresponding diagonal Hessian entries.
!  Flag such variables by setting the relevant component of map%IW to 1

      IF ( prob%Hessian_kind < 0 ) THEN
        IF ( map%h_type == 2 .OR. map%h_type == 3 ) THEN
          map%IW( : map%n ) = 1
        ELSE 
          map%IW( : map%n ) = 0
          IF ( map%h_type == 1 ) THEN
            DO i = 1, map%n
              DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
                IF ( i == prob%H%col( l ) ) map%IW( i ) = 1
              END DO
            END DO
          ELSE
            DO i = 1, map%h_ne
              IF ( prob%H%row( i ) == prob%H%col( i ) )                        &
               map%IW( prob%H%row( i ) ) = 1
            END DO
          END IF
        END IF
      END IF

!  =======================================================================
!                         Reorder variables
!  =======================================================================

!  Run through the bounds to see how many fall into each of the 
!  categories:  free(free), non-negativity(nonneg), lower(lower), range(range),
!  upper(upper), non-positivity (nonpos) and fixed (fixed);  of these, 
!  h_free, h_nonneg, h_lower, h_range, h_upper, h_nonpos and h_fixed have
!  diagonal Hessian entries

      free = 0 ; nonneg = 0 ; lower = 0
      range = 0 ; upper = 0 ; nonpos = 0 ; fixed = 0
      h_free = 0 ; h_nonneg = 0 ; h_lower = 0
      h_range = 0 ; h_upper = 0 ; h_nonpos = 0 ; h_fixed = 0

      DO i = 1, map%n
        xl = prob%X_l( i ) ; xu = prob%X_u( i )

!  fixed variable

        IF ( xu == xl ) THEN
          fixed = fixed + 1
          IF ( map%IW( i ) == 1 ) h_fixed = h_fixed + 1
          prob%X( i ) = xl
          IF ( get_z ) THEN
            IF ( apz ) prob%Z( i ) = one
            IF ( apzl ) prob%Z_l( i ) = one
            IF ( apzu ) prob%Z_u( i ) = zero
          END IF
        ELSE IF ( xl <= - control%infinity ) THEN

!  free variable

          IF ( xu >= control%infinity ) THEN
            free = free + 1
            IF ( map%IW( i ) == 1 ) h_free = h_free + 1
            IF ( get_x ) prob%X( i ) = zero
            IF ( get_z ) THEN
              IF ( apz ) prob%Z( i ) = zero
              IF ( apzl ) prob%Z_l( i ) = zero
              IF ( apzu ) prob%Z_u( i ) = zero
            END IF
          ELSE

!  non-positivity

            IF ( xu == zero .AND.                                              &
                .NOT. control% treat_zero_bounds_as_general ) THEN
              nonpos = nonpos + 1  
              IF ( map%IW( i ) == 1 ) h_nonpos = h_nonpos + 1
              IF ( get_x ) prob%X( i ) = - one
              IF ( get_z ) THEN
                IF ( apz ) prob%Z( i ) = - one
                IF ( apzl ) prob%Z_l( i ) = zero
                IF ( apzu ) prob%Z_u( i ) = - one
              END IF

!  upper bounded variable

            ELSE
              upper = upper + 1
              IF ( map%IW( i ) == 1 ) h_upper = h_upper + 1
              IF ( get_x ) prob%X( i ) = xu - one
              IF ( get_z ) THEN
                IF ( apz ) prob%Z( i ) = - one
                IF ( apzl ) prob%Z_l( i ) = zero
                IF ( apzu ) prob%Z_u( i ) = - one
              END IF
            END IF
          END IF
        ELSE
          IF ( xu < control%infinity ) THEN

!  inconsistent bounds

            IF ( xu < xl ) THEN
              inform%status = GALAHAD_error_primal_infeasible
              RETURN

!  range bounded variable

            ELSE
              range = range + 1
              IF ( map%IW( i ) == 1 ) h_range = h_range + 1
              IF ( get_x ) prob%X( i ) = half * ( xl + xu )
              IF ( get_z ) THEN
                IF ( apz ) prob%Z( i ) = zero
                IF ( apzl ) prob%Z_l( i ) = zero
                IF ( apzu ) prob%Z_u( i ) = zero
              END IF
            END IF
          ELSE

!  non-negativity

            IF ( xl == zero .AND.                                              &
                .NOT. control% treat_zero_bounds_as_general ) THEN
              nonneg = nonneg + 1  
              IF ( map%IW( i ) == 1 ) h_nonneg = h_nonneg + 1
              IF ( get_x ) prob%X( i ) = one
              IF ( get_z ) THEN
                IF ( apz ) prob%Z( i ) = one
                IF ( apzl ) prob%Z_l( i ) = one
                IF ( apzu ) prob%Z_u( i ) = zero
              END IF

!  lower bounded variable

            ELSE
              lower = lower + 1
              IF ( map%IW( i ) == 1 ) h_lower = h_lower + 1
              IF ( get_x ) prob%X( i ) = xl + one
              IF ( get_z ) THEN
                IF ( apz ) prob%Z( i ) = one
                IF ( apzl ) prob%Z_l( i ) = one
                IF ( apzu ) prob%Z_u( i ) = zero
              END IF
            END IF
          END IF
        END IF
      END DO

!  now set starting addresses for each division of the variables

      d_free = 0
      o_free = d_free + h_free
      d_nonneg = d_free + free
      o_nonneg = d_nonneg + h_nonneg
      d_lower = d_nonneg + nonneg
      o_lower = d_lower + h_lower
      d_range = d_lower + lower
      o_range = d_range + h_range
      d_upper = d_range + range
      o_upper = d_upper + h_upper
      d_nonpos = d_upper + upper
      o_nonpos = d_nonpos + h_nonpos
      d_fixed = d_nonpos + nonpos
      o_fixed = d_fixed + h_fixed

!  Also set the starting and ending addresses as required 

      dims%h_diag_end_free = o_free
      dims%x_free = d_nonneg
      dims%h_diag_end_nonneg = o_nonneg
      dims%x_l_start = d_lower + 1
      dims%h_diag_end_lower = o_lower
      dims%x_u_start = d_range + 1
      dims%h_diag_end_range = o_range
      dims%x_l_end = d_upper
      dims%h_diag_end_upper = o_upper
      dims%x_u_end = d_nonpos
      dims%h_diag_end_nonpos = o_nonpos
      prob%n = d_fixed
      map%n_reordered = prob%n
      map%h_diag_end_fixed = o_fixed

!  Allocate the inverse mapping array for A

      reallocate = .TRUE.
      IF ( ALLOCATED( map%a_map_inverse ) ) THEN
        IF ( SIZE( map%a_map_inverse ) <= MAX( map%m, map%a_ne ) ) THEN 
          DEALLOCATE( map%a_map_inverse ) ; ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( map%a_map_inverse( MAX( map%m, map%a_ne ) ),                 &
                  STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          bad_alloc = 'map%a_map_inverse' ; GO TO 900
        END IF
      END IF

!  Count how many entries there are in each constraint. Use the array
!  map%a_map_inverse for this purpose

!  Original co-ordinate storage

      IF ( map%a_type == 0 ) THEN
        map%a_map_inverse( : map%m ) = 0
        DO l = 1, map%a_ne
          i = prob%A%row( l )
          map%a_map_inverse( i ) = map%a_map_inverse( i ) + 1
        END DO

!  Original row-wise storage
  
      ELSE IF ( map%a_type == 1 ) THEN
        DO i = 1, map%m
           map%a_map_inverse( i ) = prob%A%ptr( i + 1 ) - prob%A%ptr( i )
        END DO

!  Original dense storage; record the column indices

      ELSE IF ( map%a_type == 2 ) THEN
        map%a_map_inverse( : map%m ) = prob%n
      END IF

!  Run through the constraint bounds to see how many fall into each of the 
!  categories:  free, lower, range, upper and equality

      free = 0 ; lower = 0 ; range = 0 ; upper = 0 ; equality = 0

      DO i = 1, map%m
        cl = prob%C_l( i ) ; cu = prob%C_u( i )

!  equality constraint

        IF ( cu == cl ) THEN
          IF ( map%a_map_inverse( i ) > 0 ) THEN
            equality = equality + 1
            IF ( get_y ) THEN
              IF ( apy ) prob%Y( i ) = one
              IF ( apyl ) prob%Y_l( i ) = one
              IF ( apyu ) prob%Y_u( i ) = zero
            END IF
          ELSE

!  deal with null equality constraint

            IF ( cu == zero ) THEN
              free = free + 1
              IF ( get_y ) THEN
                IF ( apy ) prob%Y( i ) = zero
                IF ( apyl ) prob%Y_l( i ) = zero
                IF ( apyu ) prob%Y_u( i ) = zero
              END IF
            ELSE
              inform%status = GALAHAD_error_primal_infeasible
              RETURN
            END IF
          END IF

        ELSE IF ( cl <= - control%infinity ) THEN

!  free constraint

          IF ( cu >= control%infinity ) THEN
            free = free + 1
            IF ( get_y ) THEN
              IF ( apy ) prob%Y( i ) = zero
              IF ( apyl ) prob%Y_l( i ) = zero
              IF ( apyu ) prob%Y_u( i ) = zero
            END IF
          ELSE

!  upper bounded constraint

            IF ( map%a_map_inverse( i ) > 0 ) THEN
              upper = upper + 1
              IF ( get_y ) THEN
                IF ( apy ) prob%Y( i ) = - one
                IF ( apyl ) prob%Y_l( i ) = zero
                IF ( apyu ) prob%Y_u( i ) = - one
              END IF
            ELSE

!  deal with null upper bounded constraint

              IF ( cu >= zero ) THEN
                free = free + 1
                IF ( get_y ) THEN
                  IF ( apy ) prob%Y( i ) = zero
                  IF ( apyl ) prob%Y_l( i ) = zero
                  IF ( apyu ) prob%Y_u( i ) = zero
                END IF
              ELSE
                inform%status = GALAHAD_error_primal_infeasible
                RETURN
              END IF
            END IF
          END IF
        ELSE
          IF ( cu < control%infinity ) THEN

!  inconsistent constraints

            IF ( cu < cl ) THEN
              inform%status = GALAHAD_error_primal_infeasible
              RETURN

!  range bounded constraint

            ELSE
              IF ( map%a_map_inverse( i ) > 0 ) THEN
                range = range + 1
                IF ( get_y ) THEN
                  IF ( apy ) prob%Y( i ) = one
                  IF ( apyl ) prob%Y_l( i ) = one
                  IF ( apyu ) prob%Y_u( i ) = zero
                END IF
              ELSE

!  deal with null range bounded constraint

                IF ( cl <= zero .AND. cu >= zero ) THEN
                  free = free + 1
                  IF ( get_y ) THEN
                    IF ( apy ) prob%Y( i ) = zero
                    IF ( apyl ) prob%Y_l( i ) = zero
                    IF ( apyu ) prob%Y_u( i ) = zero
                  END IF
                ELSE
                  inform%status = GALAHAD_error_primal_infeasible
                  RETURN
                END IF
              END IF
            END IF
          ELSE

!  lower bounded constraint

            IF ( map%a_map_inverse( i ) > 0 ) THEN
              lower = lower + 1
              IF ( get_y ) THEN
                IF ( apy ) prob%Y( i ) = one
                IF ( apyl ) prob%Y_l( i ) = one
                IF ( apyu ) prob%Y_u( i ) = zero
              END IF
            ELSE

!  deal with null lower bounded constraint

              IF ( cl <= zero ) THEN
                free = free + 1
                IF ( get_y ) THEN
                  IF ( apy ) prob%Y( i ) = zero
                  IF ( apyl ) prob%Y_l( i ) = zero
                  IF ( apyu ) prob%Y_u( i ) = zero
                END IF
              ELSE
                inform%status = GALAHAD_error_primal_infeasible
                RETURN
              END IF
            END IF
          END IF
        END IF
      END DO

!  now set starting addresses for each division of the constraints

      a_equality = 0
      a_lower = a_equality + equality
      a_range = a_lower + lower
      a_upper = a_range + range
      a_free = a_upper + upper

!  Also set the starting and ending addresses as required 

      dims%c_equality = equality
      dims%c_l_start = a_lower + 1
      dims%c_u_start = a_range + 1

      dims%c_l_end = a_upper
      dims%c_u_end = a_free
      prob%m = a_free
      map%m_reordered = prob%m

!  Allocate workspace arrays and permutation map for the variables

      reallocate = .TRUE.
      IF ( ALLOCATED( map%x_map ) ) THEN
        IF ( SIZE( map%x_map ) /= map%n ) THEN ; DEALLOCATE( map%x_map )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( map%x_map( map%n ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          bad_alloc = 'map%x_map' ; GO TO 900
        END IF
      END IF

!  Run through the variable bounds for a second time, this time building 
!  the mapping array

      DO i = 1, map%n
        xl = prob%X_l( i ) ; xu = prob%X_u( i )

!  fixed variable

        IF ( xu == xl ) THEN
          IF ( map%IW( i ) == 1 ) THEN
            d_fixed = d_fixed + 1
            map%x_map( i ) = d_fixed
          ELSE
            o_fixed = o_fixed + 1
            map%x_map( i ) = o_fixed
          END IF
        ELSE IF ( xl <= - control%infinity ) THEN

!  free variable

          IF ( xu >= control%infinity ) THEN
            IF ( map%IW( i ) == 1 ) THEN
              d_free = d_free + 1
              map%x_map( i ) = d_free
            ELSE
              o_free = o_free + 1
              map%x_map( i ) = o_free
            END IF
          ELSE

!  non-positivity

            IF ( xu == zero .AND.                                              &
                .NOT. control% treat_zero_bounds_as_general ) THEN
              IF ( map%IW( i ) == 1 ) THEN
                d_nonpos = d_nonpos + 1
                map%x_map( i ) = d_nonpos
              ELSE
                o_nonpos = o_nonpos + 1
                map%x_map( i ) = o_nonpos
              END IF

!  upper bounded variable

            ELSE
              IF ( map%IW( i ) == 1 ) THEN
                d_upper = d_upper + 1
                map%x_map( i ) = d_upper
              ELSE
                o_upper = o_upper + 1
                map%x_map( i ) = o_upper
              END IF
            END IF
          END IF
        ELSE
          IF ( xu < control%infinity ) THEN

!  range bounded variable

            IF ( map%IW( i ) == 1 ) THEN
              d_range = d_range + 1
              map%x_map( i ) = d_range
            ELSE
              o_range = o_range + 1
              map%x_map( i ) = o_range
            END IF
          ELSE

!  non-negativity

            IF ( xl == zero .AND.                                              &
                .NOT. control% treat_zero_bounds_as_general ) THEN
              IF ( map%IW( i ) == 1 ) THEN
                d_nonneg = d_nonneg + 1
                map%x_map( i ) = d_nonneg
              ELSE
                o_nonneg = o_nonneg + 1
                map%x_map( i ) = o_nonneg
              END IF

!  lower bounded variable

            ELSE
              IF ( map%IW( i ) == 1 ) THEN
                d_lower = d_lower + 1
                map%x_map( i ) = d_lower
              ELSE
                o_lower = o_lower + 1
                map%x_map( i ) = o_lower
              END IF
            END IF
          END IF
        END IF
      END DO

!  Move the gradient and bounds on variables

      IF ( prob%gradient_kind /= 0 .AND. prob%gradient_kind /= 1 )             &
        CALL SORT_inplace_permute( map%n, map%x_map, X = prob%G )
      IF ( PRESENT( parametric ) ) THEN
        IF (  ALLOCATED( prob%DG ) ) THEN
          IF ( SIZE( prob%DG ) == map%n )                                      &
            CALL SORT_inplace_permute( map%n, map%x_map, X = prob%DG )
        END IF
      END IF
      CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X_l )
      CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X_u )
      CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X )
      IF ( apz ) CALL SORT_inplace_permute( map%n, map%x_map, X = prob%Z )
      IF ( apzl ) CALL SORT_inplace_permute( map%n, map%x_map, X = prob%Z_l )
      IF ( apzu ) CALL SORT_inplace_permute( map%n, map%x_map, X = prob%Z_u )
      IF ( PRESENT( parametric ) ) THEN
        IF ( ALLOCATED( prob%DX_l ) ) THEN
          IF ( SIZE( prob%DX_l ) == map%n )                                    &
            CALL SORT_inplace_permute( map%n, map%x_map, X = prob%DX_l )
        END IF
        IF ( ALLOCATED( prob%DX_u ) ) THEN
          IF ( SIZE( prob%DX_u ) == map%n )                                    &
            CALL SORT_inplace_permute( map%n, map%x_map, X = prob%DX_u )
        END IF
      END IF

!  If the problem is a weighted least squares one, permute the weights

      IF ( prob%Hessian_kind > 1 ) THEN
        CALL SORT_inplace_permute( map%n, map%x_map, X = prob%WEIGHT )
        CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X0 )
      ELSE IF ( prob%Hessian_kind > 0 ) THEN
        CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X0 )

!  Permute the rows and columns for a general H. Start by counting how 
!  many entries will be required for each row. map%IW(i) gives the number 
!  of entries in row i

      ELSE IF ( prob%Hessian_kind < 0 ) THEN
        IF ( map%h_type == 3 ) THEN

!  Original diagonal storage; record the column indices

          map%IW( : map%n ) = 1
        ELSE IF ( map%h_type == 2 ) THEN

!  Original dense storage; record the column indices

          DO ll = 1, map%n
            map%IW( ll ) = ll
          END DO
        ELSE 
          map%IW( : map%n ) = 0
          IF ( map%h_type == 0 ) THEN

!  Original co-ordinate storage

            DO l = 1, map%h_ne
              i = map%x_map( prob%H%row( l ) ) ; j = map%x_map( prob%H%col( l ) )
              IF ( i >= j ) THEN
                map%IW( i ) = map%IW( i ) + 1
              ELSE
                map%IW( j ) = map%IW( j ) + 1
              END IF
            END DO
          ELSE

!  Original row-wise storage
  
            DO k = 1, map%n
              i = map%x_map( k )
              DO l = prob%H%ptr( k ), prob%H%ptr( k + 1 ) - 1
                j = map%x_map( prob%H%col( l ) )
                IF ( i >= j ) THEN
                  map%IW( i ) = map%IW( i ) + 1
                ELSE
                  map%IW( j ) = map%IW( j ) + 1
                END IF
              END DO
            END DO
          END IF
        END IF

!  Set starting addresses prior to making the permutation

        j = 1
        DO i = 1, map%n
          k = j
          j = j + map%IW( i )
          map%IW( i ) = k
        END DO

!  Allocate the inverse mapping array for H

        reallocate = .TRUE.
        IF ( ALLOCATED( map%h_map_inverse ) ) THEN
          IF ( SIZE( map%h_map_inverse ) /= map%h_ne ) THEN 
            DEALLOCATE( map%h_map_inverse ) ; ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( map%h_map_inverse( map%h_ne ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            bad_alloc = 'map%h_map_inverse' ; GO TO 900
          END IF
        END IF

        IF ( map%h_type == 0 .OR. map%h_type == 2 .OR. map%h_type == 3 ) THEN

!  Ensure that there is enough space to store the result in a row-wise scheme

          reallocate = .TRUE.
          IF ( ALLOCATED( prob%H%col ) ) THEN
            IF ( SIZE( prob%H%col ) < map%h_ne ) THEN 
              DEALLOCATE( prob%H%col ) ; ELSE ; reallocate = .FALSE.
            END IF
          END IF
          IF ( reallocate ) THEN 
            ALLOCATE( prob%H%col( map%h_ne ), STAT = inform%alloc_status )
            IF ( inform%alloc_status /= 0 ) THEN 
              bad_alloc = 'prob%H%col' ; GO TO 900
            END IF
          END IF

          reallocate = .TRUE.
          IF ( ALLOCATED( prob%H%ptr ) ) THEN
            IF ( SIZE( prob%H%ptr ) < map%n * 1 ) THEN 
              DEALLOCATE( prob%H%ptr ) ; ELSE ; reallocate = .FALSE.
            END IF
          END IF
          IF ( reallocate ) THEN 
            ALLOCATE( prob%H%ptr( map%n + 1 ), STAT = inform%alloc_status )
            IF ( inform%alloc_status /= 0 ) THEN 
              bad_alloc = 'prob%H%ptr' ; GO TO 900
            END IF
          END IF
        END IF

!  Reorder the rows; compute the mapping array for H and renumber its columns

        IF ( map%h_type == 3 ) THEN

!  Original diagonal storage

          DO ll = 1, map%n
            i = map%x_map( ll )
            map%h_map_inverse( map%IW( i ) ) = ll
            prob%H%col( ll ) = i
            map%IW( i ) = map%IW( i ) + 1
          END DO
        ELSE IF ( map%h_type == 2 ) THEN

!  Original dense storage

          l = 0
          DO ll = 1, map%n
            i = map%x_map( ll )
            DO k = 1, ll
              j = map%x_map( k )
              l = l + 1
              IF ( i >= j ) THEN
                map%h_map_inverse( map%IW( i ) ) = l
                prob%H%col( l ) = j
                map%IW( i ) = map%IW( i ) + 1
              ELSE
                map%h_map_inverse( map%IW( j ) ) = l
                prob%H%col( l ) = i
                map%IW( j ) = map%IW( j ) + 1
              END IF
            END DO
          END DO
        ELSE 
          IF ( map%h_type == 0 ) THEN

!  Original co-ordinate storage

            DO l = 1, map%h_ne
              i = map%x_map( prob%H%row( l ) ) ; j = map%x_map( prob%H%col( l ) )
              IF ( i >= j ) THEN
                map%h_map_inverse( map%IW( i ) ) = l
                prob%H%col( l ) = j
                map%IW( i ) = map%IW( i ) + 1
              ELSE
                map%h_map_inverse( map%IW( j ) ) = l
                prob%H%col( l ) = i
                map%IW( j ) = map%IW( j ) + 1
              END IF
            END DO
          ELSE

!  Original row-wise storage
  
            DO k = 1, map%n
              i = map%x_map( k )
              DO l = prob%H%ptr( k ), prob%H%ptr( k + 1 ) - 1
                j = map%x_map( prob%H%col( l ) )
                IF ( i >= j ) THEN
                  map%h_map_inverse( map%IW( i ) ) = l
                  prob%H%col( l ) = j
                  map%IW( i ) = map%IW( i ) + 1
                ELSE
                  map%h_map_inverse( map%IW( j ) ) = l
                  prob%H%col( l ) = i
                  map%IW( j ) = map%IW( j ) + 1
                END IF
              END DO
            END DO
          END IF
        END IF

!  Set the starting addresses for each row in the permuted matrix

        prob%H%ptr( 1 ) = 1
        DO i = 1, map%n
          prob%H%ptr( i + 1 ) = map%IW( i )
        END DO

!  Apply the reordering to H

        CALL SORT_inverse_permute( map%h_ne, map%h_map_inverse,               &
                                   X = prob%H%val, IX = prob%H%col )

!  Reorder the columns so that they appear in increasing order

        CALL QPP_order_rows( map%n, prob%H%val, prob%H%col,                    &
                             prob%H%ptr, map%h_map_inverse )

!  If the original storage scheme was by co-ordinates, record the
!  final row indices

        IF ( map%h_type == 0 ) THEN
          DO i = 1, prob%n
            prob%H%row( prob%H%ptr( i ) : prob%H%ptr( i + 1 ) - 1 ) = i
          END DO
        END IF
      END IF

!  =======================================================================
!                         Reorder constraints
!  =======================================================================

!  Allocate permutation map for the constraints

      reallocate = .TRUE.
      IF ( ALLOCATED( map%c_map ) ) THEN
        IF ( SIZE( map%c_map ) /= map%m ) THEN ; DEALLOCATE( map%c_map )
        ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( map%c_map( map%m ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          bad_alloc = 'map%c_map' ; GO TO 900
        END IF
      END IF

!  Run through the bounds for a second time, this time building the mapping
!  array

      DO i = 1, map%m
        cl = prob%C_l( i ) ; cu = prob%C_u( i )

!  equality constraint

        IF ( cu == cl ) THEN
          IF ( map%a_map_inverse( i ) > 0 ) THEN
            a_equality = a_equality + 1
            map%c_map( i ) = a_equality
          ELSE
            a_free = a_free + 1
            map%c_map( i ) = a_free
          END IF
        ELSE IF ( cl <= - control%infinity ) THEN

!  free constraint

          IF ( cu >= control%infinity ) THEN
            a_free = a_free + 1
            map%c_map( i ) = a_free
          ELSE

!  upper bounded constraint

            IF ( map%a_map_inverse( i ) > 0 ) THEN
              a_upper = a_upper + 1
              map%c_map( i ) = a_upper
            ELSE
              a_free = a_free + 1
              map%c_map( i ) = a_free
            END IF
          END IF
        ELSE

          IF ( cu < control%infinity ) THEN

!  range bounded constraint

            IF ( map%a_map_inverse( i ) > 0 ) THEN
              a_range = a_range + 1
              map%c_map( i ) = a_range
            ELSE
              a_free = a_free + 1
              map%c_map( i ) = a_free
            END IF
          ELSE

!  lower bounded constraint

            IF ( map%a_map_inverse( i ) > 0 ) THEN
              a_lower = a_lower + 1
              map%c_map( i ) = a_lower
            ELSE
              a_free = a_free + 1
              map%c_map( i ) = a_free
            END IF
          END IF
        END IF
      END DO

!  Move the constraint values and bounds

      CALL SORT_inplace_permute( map%m, map%c_map, X = prob%C_l )
      CALL SORT_inplace_permute( map%m, map%c_map, X = prob%C_u )
      IF ( apy ) CALL SORT_inplace_permute( map%m, map%c_map, X = prob%Y )
      IF ( apyl ) CALL SORT_inplace_permute( map%m, map%c_map, X = prob%Y_l )
      IF ( apyu ) CALL SORT_inplace_permute( map%m, map%c_map, X = prob%Y_u )
      IF ( PRESENT( parametric ) ) THEN
        IF ( ALLOCATED( prob%DC_l ) ) THEN
          IF ( SIZE( prob%DC_l ) == map%m )                                    &
            CALL SORT_inplace_permute( map%m, map%c_map, X = prob%DC_l )
        END IF
        IF ( ALLOCATED( prob%DC_u ) ) THEN
          IF ( SIZE( prob%DC_u ) == map%m )                                    &
            CALL SORT_inplace_permute( map%m, map%c_map, X = prob%DC_u )
        END IF
      END IF

!  Now permute the rows and columns of A. Start by counting how many entries
!  will be required for each row. map%IW(i) gives the number of entries in row i
!  Also record the number of nonzeros in each COLUMN corresponding to
!  fixed variables. map%ptr_a_fixed(j) gives the number in permuted column j

      IF ( prob%n < map%n ) THEN
        reallocate = .TRUE.
        IF ( ALLOCATED( map%ptr_a_fixed ) ) THEN
        IF ( LBOUND( map%ptr_a_fixed, 1 ) /= prob%n + 1 .OR.                   &
             UBOUND( map%ptr_a_fixed, 1 ) /= map%n + 1 ) THEN 
          DEALLOCATE( map%ptr_a_fixed )
          ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( map%ptr_a_fixed(  prob%n + 1 : map%n + 1 ),                &
                    STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            bad_alloc = 'map%ptr_a_fixed' ; GO TO 900
          END IF
        END IF
        IF ( inform%alloc_status /= 0 ) GO TO 900
      END IF

      IF ( map%a_type == 2 ) THEN

!  Original dense storage; record the column indices

        map%IW( : map%m ) = prob%n
        IF ( prob%n < map%n ) map%ptr_a_fixed( prob%n + 1 : map%n ) = map%m
      ELSE 
        map%IW( : map%m ) = 0
        IF ( prob%n < map%n ) map%ptr_a_fixed( prob%n + 1 : map%n ) = 0
        IF ( map%a_type == 0 ) THEN

!  Original co-ordinate storage

          DO l = 1, map%a_ne
            j = map%x_map( prob%A%col( l ) )
            IF ( j <= prob%n ) THEN
              i = map%c_map( prob%A%row( l ) )
              map%IW( i ) = map%IW( i ) + 1
            ELSE
              map%ptr_a_fixed( j ) = map%ptr_a_fixed( j ) + 1
            END IF
          END DO
        ELSE

!  Original row-wise storage
  
          DO k = 1, map%m
            i = map%c_map( k )
            DO l = prob%A%ptr( k ), prob%A%ptr( k + 1 ) - 1
              j = map%x_map( prob%A%col( l ) )
              IF ( j <= prob%n ) THEN
                map%IW( i ) = map%IW( i ) + 1
              ELSE
                map%ptr_a_fixed( j ) = map%ptr_a_fixed( j ) + 1
              END IF
            END DO
          END DO
        END IF
      END IF

!  Set starting addresses prior to making the permutation

      j = 1
      DO i = 1, map%m
        k = j
        j = j + map%IW( i )
        map%IW( i ) = k
      END DO

      DO i = prob%n + 1, map%n
        k = j
        j = j + map%ptr_a_fixed( i )
        map%ptr_a_fixed( i ) = k
      END DO      

!  Allocate the inverse mapping array for A

      reallocate = .TRUE.
      IF ( ALLOCATED( map%a_map_inverse ) ) THEN
        IF ( SIZE( map%a_map_inverse ) /= map%a_ne ) THEN 
          DEALLOCATE( map%a_map_inverse ) ; ELSE ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( map%a_map_inverse( map%a_ne ), STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN 
          bad_alloc = 'map%a_map_inverse' ; GO TO 900
        END IF
      END IF

!  Reorder the rows; compute the mapping array for A and renumber its columns
!  NB. Any columns corresponding to FIXED variables, will have been
!  moved to the end of A, and will be stored by COLUMN not by row. In 
!  particular, A_col for these entries gives the row and not the column
!  number

      IF ( map%a_type == 0 .OR. map%a_type == 2 ) THEN

!  Ensure that there is enough space to store the result in a row-wise scheme

        reallocate = .TRUE.
        IF ( ALLOCATED( prob%A%col ) ) THEN
          IF ( SIZE( prob%A%col ) < map%a_ne ) THEN 
            DEALLOCATE( prob%A%col ) ; ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( prob%A%col( map%a_ne ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            bad_alloc = 'prob%A%col' ; GO TO 900
          END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( prob%A%ptr ) ) THEN
          IF ( SIZE( prob%A%ptr ) < map%m * 1 ) THEN 
            DEALLOCATE( prob%A%ptr ) ; ELSE ; reallocate = .FALSE.
          END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( prob%A%ptr( map%m + 1 ), STAT = inform%alloc_status )
          IF ( inform%alloc_status /= 0 ) THEN 
            bad_alloc = 'prob%A%ptr' ; GO TO 900
          END IF
        END IF
      END IF

      IF ( map%a_type == 2 ) THEN

!  Original dense storage

        l = 0
        DO ll = 1, map%m
          DO k = 1, map%n
            l = l + 1
            i = map%c_map( ll ) ; j = map%x_map( k )
            IF ( j <= prob%n ) THEN
              map%a_map_inverse( map%IW( i ) ) = l
              prob%A%col( l ) = j
              map%IW( i ) = map%IW( i ) + 1
            ELSE
              map%a_map_inverse( map%ptr_a_fixed( j ) ) = l
              prob%A%col( l ) = i
              map%ptr_a_fixed( j ) = map%ptr_a_fixed( j ) + 1
            END IF
          END DO
        END DO
      ELSE 
        IF ( map%a_type == 0 ) THEN

!  Original co-ordinate storage

          DO l = 1, map%a_ne
            i = map%c_map( prob%A%row( l ) ) ; j = map%x_map( prob%A%col( l ) )
            IF ( j <= prob%n ) THEN
              map%a_map_inverse( map%IW( i ) ) = l
              prob%A%col( l ) = j
              map%IW( i ) = map%IW( i ) + 1
            ELSE
              map%a_map_inverse( map%ptr_a_fixed( j ) ) = l
              prob%A%col( l ) = i
              map%ptr_a_fixed( j ) = map%ptr_a_fixed( j ) + 1
            END IF
          END DO

        ELSE

!  Original row-wise storage
  
          DO k = 1, map%m
            i = map%c_map( k )
            DO l = prob%A%ptr( k ), prob%A%ptr( k + 1 ) - 1
              j = map%x_map( prob%A%col( l ) )
              IF ( j <= prob%n ) THEN
                map%a_map_inverse( map%IW( i ) ) = l
                prob%A%col( l ) = j
                map%IW( i ) = map%IW( i ) + 1
              ELSE
                map%a_map_inverse( map%ptr_a_fixed( j ) ) = l
                prob%A%col( l ) = i
                map%ptr_a_fixed( j ) = map%ptr_a_fixed( j ) + 1
              END IF
            END DO
          END DO
        END IF
      END IF

!  Set the starting addresses for each row in the permuted matrix

      prob%A%ptr( 1 ) = 1
      DO i = 1, map%m
        prob%A%ptr( i + 1 ) = map%IW( i )
      END DO

      IF ( prob%n < map%n ) THEN
        DO i = map%n, prob%n + 1, - 1
          map%ptr_a_fixed( i + 1 ) = map%ptr_a_fixed( i ) 
        END DO      
        map%ptr_a_fixed( prob%n + 1 ) = prob%A%ptr( map%m + 1 )
      END IF

!  Apply the reordering to A

      CALL SORT_inverse_permute( map%a_ne, map%a_map_inverse, X = prob%A%val,  &
                                 IX = prob%A%col )

!  Reorder the columns so that they appear in increasing order

      CALL QPP_order_rows( map%m, prob%A%val, prob%A%col, prob%A%ptr,          &
                            map%a_map_inverse )

!  If the original storage scheme was by co-ordinates, record the
!  final row indices

      IF ( map%a_type == 0 ) THEN
        DO i = 1, prob%m
          prob%A%row( prob%A%ptr( i ) : prob%A%ptr( i + 1 ) - 1 ) = i
        END DO
      END IF

!  Now evaluate c = A x corresponding to the free variables

      prob%C( : map%m ) = zero
      CALL QPP_AX( map, prob%X, prob%A%val, prob%A%col,                       &
                    prob%A%ptr, map%m, prob%C( : map%m ) )

      IF ( prob%n < map%n ) THEN

!  If a compact gradient format is specified but the Hessian is general,
!  move to a general gradient format to accomodate any fixed variables

        IF ( ( prob%gradient_kind == 0 .OR. prob%gradient_kind == 1 )         &
             .AND. prob%Hessian_kind < 0 ) THEN
          reallocate = .TRUE.
          IF ( ALLOCATED( prob%G ) ) THEN
            IF ( SIZE( prob%G ) < map%n ) THEN 
              DEALLOCATE( prob%G ) ; ELSE ; reallocate = .FALSE.
            END IF
          END IF
          IF ( reallocate ) THEN 
            ALLOCATE( prob%G( map%n ), STAT = inform%alloc_status )
            IF ( inform%alloc_status /= 0 ) THEN 
              bad_alloc = 'prob%G' ; GO TO 900
            END IF
          END IF
          IF ( prob%gradient_kind == 0 ) THEN
            prob%G = zero
          ELSE
            prob%G = one
          END IF
          prob%gradient_kind = 2
        END IF

!  Transform f, g and the bounds on the constraints to account for 
!  fixed variables

        CALL QPP_remove_fixed( map, prob, f = .TRUE., g = .TRUE.,             &
                               c_bounds = .TRUE. )
      END IF

      map%set = .TRUE.
      map%h_perm = .TRUE.
      map%a_perm = .TRUE.
      prob%A%ne = prob%A%ptr( prob%m + 1 ) - 1
!     CALL QPT_put_A( prob%A%type, 'COORDINATE' )
      CALL QPT_put_A( prob%A%type, 'SPARSE_BY_ROWS' )
      IF ( prob%Hessian_kind < 0 ) THEN
        prob%H%ne = prob%H%ptr( prob%n + 1 ) - 1
!       CALL QPT_put_H( prob%H%type, 'COORDINATE' )
        CALL QPT_put_H( prob%H%type, 'SPARSE_BY_ROWS' )
      END IF

      inform%status = GALAHAD_ok
      RETURN  

!  Error returns

  900 CONTINUE
      inform%status = GALAHAD_error_allocate
      IF ( control%error > 0 )                                                 &
        WRITE( control%error, 2900 ) bad_alloc, inform%alloc_status

      RETURN  

!  Non-executable statements

 2900 FORMAT( ' ** Message from -QPP_reorder-', /,                             &
              ' Allocation error, for ', A20, /, ' status = ', I6 ) 

!  End of QPP_reorder

      END SUBROUTINE QPP_reorder

!-*-*-*-*-*-*-*-   Q P P _ a p p l y   S U B R O U T I N E   -*-*-*-*-*-*-*-

      SUBROUTINE QPP_apply( map, inform, prob, get_all,                        &
                            get_all_parametric, get_f, get_g,                  &
                            get_dg, get_x, get_y, get_z, get_x_bounds,         &
                            get_dx_bounds, get_c, get_c_bounds, get_dc_bounds, &
                            get_A, get_H )

!      .........................................
!      .                                       .
!      .  Apply the permutations computed in   .
!      .  QPP_reorder to the data prob as      .
!      .  appropriate                          .
!      .                                       .
!      .........................................

!  Arguments:
!  =========
!
!   prob    see Subroutine QPP_reorder
!   map     see Subroutine QPP_initialize
!
!   inform is a structure of type QPP_inform_type that provides 
!    information on exit from QPP_apply. The component status 
!    has possible values:
!  
!     0 Normal termination.
!
!     1 The mapping arrays have not yet been set. Either QPP_reorder
!       has not yet been called, or the call was unsuccessful.
!
!   get_all       LOGICAL, OPTIONAL. If present, process the entire problem
!   get_all_parametric  LOGICAL, OPTIONAL. If present, process the entire 
!                 problem including parametric parts
!   get_f         LOGICAL, OPTIONAL. If present, process f
!   get_g         LOGICAL, OPTIONAL. If present, process g
!   get_dg        LOGICAL, OPTIONAL. If present, process dg
!   get_c         LOGICAL, OPTIONAL. If present, process c
!   get_x_bounds  LOGICAL, OPTIONAL. If present, process x_l and x_u 
!   get_dx_bounds LOGICAL, OPTIONAL. If present, process dx_l and dx_u 
!   get_c_bounds  LOGICAL, OPTIONAL. If present, process c_l and c_u
!   get_dc_bounds LOGICAL, OPTIONAL. If present, process dc_l and dc_u
!   get_x         LOGICAL, OPTIONAL. If present, process x
!   get_y         LOGICAL, OPTIONAL. If present, process y
!   get_z         LOGICAL, OPTIONAL. If present, process z
!   get_A         LOGICAL, OPTIONAL. If present, process A
!   get_H         LOGICAL, OPTIONAL. If present, process H or weights/x0

!  Dummy arguments

      TYPE ( QPP_map_type ), INTENT( INOUT ) :: map
      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      TYPE ( QPP_inform_type ), INTENT( OUT ) :: inform
      LOGICAL, OPTIONAL, INTENT( IN ) :: get_all, get_g, get_dg,               &
                                         get_c, get_x_bounds, get_c_bounds,    &
                                         get_dx_bounds, get_dc_bounds,         &
                                         get_x, get_y, get_z, get_A, get_H,    &
                                         get_all_parametric, get_f

!  Local variables

      INTEGER :: i, j, k, l, ll
      LOGICAL :: apy, apyl, apyu, apz, apzl, apzu

!  Check that Y or Y_l/Y_u and Z or Z_l/Z_u has been allocated

      apy = ALLOCATED( prob%Y )
      apyl = ALLOCATED( prob%Y_l ) ; apyu = ALLOCATED( prob%Y_u )
      IF ( .NOT. ( apy .OR. ( apyl .AND. apyu ) ) ) THEN
        inform%status = GALAHAD_error_y_unallocated ; RETURN
      END IF

      apz = ALLOCATED( prob%Z )
      apzl = ALLOCATED( prob%Z_l ) ; apzu = ALLOCATED( prob%Z_u )
      IF ( .NOT. ( apz .OR. ( apzl .AND. apzu ) ) ) THEN
        inform%status = GALAHAD_error_z_unallocated ; RETURN
      END IF

!  Check to see that the mapping arrays have been set

      IF ( .NOT. map%set ) THEN
        inform%status = GALAHAD_error_unordered
        RETURN
      END IF

!  Check that the variable and constraint bounds are consistent

      DO i = 1, map%n
        IF ( prob%X_l( i ) > prob%X_u( i ) ) THEN
          inform%status = GALAHAD_error_primal_infeasible
          RETURN
        END IF
      END DO

      DO i = 1, map%m
        IF ( prob%C_l( i ) > prob%C_u( i ) ) THEN
          inform%status = GALAHAD_error_primal_infeasible
          RETURN
        END IF
      END DO

!  Check to see that storage formats have not changed

      IF ( ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.        &
             PRESENT( get_A ) ) .AND. .NOT. map%a_perm ) THEN
        IF ( ( map%a_type == 2 .AND.                                           &
                 SMT_get( prob%A%type ) /= 'DENSE' ) .OR.                      &
             ( map%a_type == 1 .AND.                                           &
                 SMT_get( prob%A%type ) /= 'SPARSE_BY_ROWS' ) .OR.             &
             ( map%a_type == 0 .AND.                                           &
                 ( SMT_get( prob%A%type ) /= 'COORDINATE' .OR.                 &
                   prob%A%ne /= map%a_ne ) ) ) THEN
          inform%status = GALAHAD_error_reformat
          RETURN
        END IF  
      END IF  

      IF ( ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.        &
             PRESENT( get_H ) )                                                &
            .AND. prob%Hessian_kind < 0 .AND. .NOT. map%h_perm ) THEN
        IF ( ( map%h_type == 3 .AND.                                           &
                 SMT_get( prob%H%type ) /= 'DIAGONAL' ) .OR.                   &
             ( map%h_type == 2 .AND.                                           &
                 SMT_get( prob%H%type ) /= 'DENSE' ) .OR.                      &
             ( map%h_type == 1 .AND.                                           &
                 SMT_get( prob%H%type ) /= 'SPARSE_BY_ROWS' ) .OR.             &
             ( map%h_type == 0 .AND.                                           &
                 ( SMT_get( prob%H%type ) /= 'COORDINATE' .OR.                 &
                   prob%H%ne /= map%h_ne ) ) ) THEN
          inform%status = GALAHAD_error_reformat
          RETURN
        END IF  
      END IF  

!  Pick up the correct dimensions

      prob%n = map%n_reordered
      prob%m = map%m_reordered

!  Map A

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_A ) ) THEN

!  The row/column permutations have already been made

        IF ( map%a_perm ) THEN
          CALL SORT_inverse_permute( map%a_ne, map%a_map_inverse,              &
                                     X = prob%A%val )

!  The row/column indices are in their original order

        ELSE

!  Original dense storage; record the column indices

          IF ( map%a_type == 2 ) THEN

!  Compute the number of entries in each row of A, and renumber its columns.
!  NB. Any columns corresponding to FIXED variables, will have been
!  moved to the end of A, and will be stored by COLUMN not by row. In 
!  particular, A_col for these entries gives the row and not the column
!  number

            l = 0
            DO ll = 1, map%m
              i = map%c_map( ll )
              map%IW( ll ) = prob%n
              DO k = 1, map%n
                j = map%x_map( k )
                l = l + 1
                IF ( j <= prob%n ) THEN
                  prob%A%col( l ) = j
                ELSE
                  prob%A%col( l ) = i
                END IF
              END DO
            END DO
          ELSE 
            map%IW( : map%m ) = 0
            IF ( map%a_type == 0 ) THEN

!  Original co-ordinate storage

              DO l = 1, map%a_ne
                i = map%c_map( prob%A%row( l ) )
                j = map%x_map( prob%A%col( l ) )
                IF ( j <= prob%n ) THEN
                  map%IW( i ) = map%IW( i ) + 1
                  prob%A%col( l ) = j
                ELSE
                  prob%A%col( l ) = i
                END IF
              END DO
            ELSE

!  Original row-wise storage
  
              DO k = 1, map%m
                i = map%c_map( k )
                DO l = prob%A%ptr( k ), prob%A%ptr( k + 1 ) - 1
                  j = map%x_map( prob%A%col( l ) )
                  IF ( j <= prob%n ) THEN
                    map%IW( i ) = map%IW( i ) + 1
                    prob%A%col( l ) = j
                  ELSE
                    prob%A%col( l ) = i
                  END IF
                END DO
              END DO
            END IF
          END IF

!  Set the starting addresses for each row in the permuted matrix

          prob%A%ptr( 1 ) = 1
          DO i = 1, map%m
            prob%A%ptr( i + 1 ) = prob%A%ptr( i ) + map%IW( i )
          END DO
    
!  Apply the reordering to A

          CALL SORT_inverse_permute( map%a_ne, map%a_map_inverse,              &
                                     X = prob%A%val, IX = prob%A%col )

!  If the original storage scheme was by co-ordinates, record the
!  final row indices

          IF ( map%a_type == 0 ) THEN
            DO i = 1, prob%m
              prob%A%row( prob%A%ptr( i ) : prob%A%ptr( i + 1 ) - 1 ) = i
            END DO
          END IF
          map%a_perm = .TRUE.
          prob%A%ne = prob%A%ptr( prob%m + 1 ) - 1
        END IF
      END IF

!  Map H or weights/x0

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.        &
           PRESENT( get_H ) ) THEN

!  Weighted-least-distance Hessian

        IF ( prob%Hessian_kind > 1 ) THEN
          CALL SORT_inplace_permute( map%n, map%x_map, X = prob%WEIGHT )
          CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X0 )
        ELSE IF ( prob%Hessian_kind > 0 ) THEN
          CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X0 )
 
!  General Hessian

        ELSE IF ( prob%Hessian_kind < 0 ) THEN  

!  The row/column permutations have already been made

          IF ( map%h_perm ) THEN
            CALL SORT_inverse_permute( map%h_ne, map%h_map_inverse,            &
                                       X = prob%H%val )

!  The row/column indices are in their original order

          ELSE

!  Compute the number of entries in each row of H, and renumber its columns

            IF ( map%h_type == 3 ) THEN

!  Original diagonal storage; record the column indices

              DO l = 1, map%n
                i = map%x_map( l )
                map%IW( i ) = 1
                prob%H%col( l ) = i
              END DO
            ELSE IF ( map%h_type == 2 ) THEN

!  Original dense storage; record the column indices

              l = 0
              DO ll = 1, map%n
                i = map%x_map( ll ) 
                map%IW( ll ) = ll
                DO k = 1, ll
                  l = l + 1
                  prob%H%col( l ) = MIN( i, map%x_map( k ) )
                END DO
              END DO
            ELSE 
              map%IW( : map%n ) = 0
              IF ( map%h_type == 0 ) THEN

!  Original co-ordinate storage

                DO l = 1, map%h_ne
                  i = map%x_map( prob%H%row( l ) )
                  j = map%x_map( prob%H%col( l ) )
                  IF ( i >= j ) THEN
                    map%IW( i ) = map%IW( i ) + 1
                    prob%H%col( l ) = j
                  ELSE
                    map%IW( j ) = map%IW( j ) + 1
                    prob%H%col( l ) = i
                  END IF
                END DO
              ELSE

!  Original row-wise storage
  
                DO k = 1, map%n
                  i = map%x_map( k )
                  DO l = prob%H%ptr( k ), prob%H%ptr( k + 1 ) - 1
                    j = map%x_map( prob%H%col( l ) )
                    IF ( i >= j ) THEN
                      map%IW( i ) = map%IW( i ) + 1
                      prob%H%col( l ) = j
                    ELSE
                      map%IW( j ) = map%IW( j ) + 1
                      prob%H%col( l ) = i
                    END IF
                  END DO
                END DO
              END IF
            END IF

!  Set the starting addresses for each row in the permuted matrix

            prob%H%ptr( 1 ) = 1
            DO i = 1, map%n
              prob%H%ptr( i + 1 ) = prob%H%ptr( i ) + map%IW( i )
            END DO

!  Apply the reordering to H

            CALL SORT_inverse_permute( map%h_ne,  map%h_map_inverse,          &
                                       X = prob%H%val, IX = prob%H%col )

!  If the original storage scheme was by co-ordinates, record the
!  final row indices

            IF ( map%h_type == 0 ) THEN
              DO i = 1, prob%n
                prob%H%row( prob%H%ptr( i ) : prob%H%ptr( i + 1 ) - 1 ) = i
              END DO
            END IF
            map%h_perm = .TRUE.
            prob%H%ne = prob%H%ptr( prob%n + 1 ) - 1
          END IF
        END IF
      END IF

!  Check to see if permuted A and H are available 

      IF ( ( ( PRESENT( get_g ) .OR. PRESENT( get_c_bounds ) )                 &
               .AND. prob%n < map%n ) .AND. .NOT.                              &
           ( map%a_perm .AND. ( map%h_perm .AND. prob%Hessian_kind < 0 ) ) ) THEN
        inform%status = GALAHAD_error_ah_unordered
        RETURN
      END IF
      IF ( ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.        &
               PRESENT( get_c ) ) .AND. .NOT. map%a_perm ) THEN
        inform%status = GALAHAD_error_ah_unordered
        RETURN
      END IF

!  Permute the bounds on x

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.        &
           PRESENT( get_x_bounds ) ) THEN
        CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X_l )
        CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X_u )
      END IF

!  Permute the bounds on dx

      IF ( PRESENT( get_all_parametric ) .OR. PRESENT( get_dx_bounds ) ) THEN
        IF ( ALLOCATED( prob%DX_l ) ) THEN
          IF ( SIZE( prob%DX_l ) == map%n )                                    &
            CALL SORT_inplace_permute( map%n, map%x_map, X = prob%DX_l )
        END IF
        IF ( ALLOCATED( prob%DX_u ) ) THEN
          IF ( SIZE( prob%DX_u ) == map%n )                                    &
            CALL SORT_inplace_permute( map%n, map%x_map, X = prob%DX_u )
        END IF
      END IF

!  Permute the bounds on x

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_x ) ) THEN
        CALL SORT_inplace_permute( map%n, map%x_map, X = prob%X )
      END IF

!  Permute the bounds on y

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_y ) ) THEN
        IF ( apy ) CALL SORT_inplace_permute( map%m, map%c_map, X = prob%Y )
        IF ( apyl ) CALL SORT_inplace_permute( map%m, map%c_map, X = prob%Y_l )
        IF ( apyu ) CALL SORT_inplace_permute( map%m, map%c_map, X = prob%Y_u )
      END IF

!  Permute the bounds on z

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_z ) ) THEN
        IF ( apz ) CALL SORT_inplace_permute( map%n, map%x_map, X = prob%Z )
        IF ( apzl ) CALL SORT_inplace_permute( map%n, map%x_map, X = prob%Z_l )
        IF ( apzu ) CALL SORT_inplace_permute( map%n, map%x_map, X = prob%Z_u )
      END IF

!  Permute the bounds on c

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_c_bounds ) ) THEN
        CALL SORT_inplace_permute( map%m, map%c_map, X = prob%C_l )
        CALL SORT_inplace_permute( map%m, map%c_map, X = prob%C_u )
      END IF

!  Permute the bounds on dc

      IF ( PRESENT( get_all_parametric ) .OR. PRESENT( get_dc_bounds ) ) THEN
        IF ( ALLOCATED( prob%DC_l ) ) THEN
          IF ( SIZE( prob%DC_l ) == map%m )                                    &
            CALL SORT_inplace_permute( map%m, map%c_map, X = prob%DC_l )
        END IF
        IF ( ALLOCATED( prob%DC_u ) ) THEN
          IF ( SIZE( prob%DC_u ) == map%m )                                    &
            CALL SORT_inplace_permute( map%m, map%c_map, X = prob%DC_u )
        END IF
      END IF

!  Permute g

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_g ) ) THEN
        IF ( prob%gradient_kind /= 0 .AND. prob%gradient_kind /= 1 )           &
          CALL SORT_inplace_permute( map%n, map%x_map, X = prob%G )
      END IF

!  Permute dg

      IF ( PRESENT( get_all_parametric ) .OR. PRESENT( get_dg ) ) THEN
        IF ( ALLOCATED( prob%DG ) ) THEN
          IF ( SIZE( prob%DG ) == map%n )                                      &
            CALL SORT_inplace_permute( map%n, map%x_map, X = prob%DG )
        END IF
      END IF

!  Form c = A * x

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_c ) ) THEN
        prob%C( : map%m ) = zero
        CALL QPP_AX( map, prob%X, prob%A%val, prob%A%col, prob%A%ptr,          &
                      map%m, prob%C( : map%m ) )

      END IF

!  Transform f, g and the bounds on the constraints to account for 
!  fixed variables

      IF ( ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) ) .AND.     &
             prob%n < map%n ) THEN
        CALL QPP_remove_fixed( map, prob,                                      &
                               f =.TRUE., g =.TRUE., c_bounds = .TRUE. )
      ELSE IF ( ( PRESENT( get_f ) .OR. PRESENT( get_g ) .OR.                  &
                  PRESENT( get_c_bounds ) ) .AND. prob%n < map%n ) THEN
        CALL QPP_remove_fixed( map, prob,                                      &
                               f = PRESENT( get_f ), g = PRESENT( get_g ),     &
                               c_bounds = PRESENT( get_c_bounds ) )
      END IF

      inform%status = GALAHAD_ok
      RETURN  

!  End of QPP_apply

      END SUBROUTINE QPP_apply

!-*-*-*-*-*-   Q P P _ g e t _ v a l u e s  S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE QPP_get_values( map, inform, prob, X_val, Y_val, Z_val )

!      ................................................
!      .                                              .
!      .  Recover the values of x (primal variables), .
!      .  y (Lagrange multipliers for constraints),   .
!      .  and z (dual variables)                      .
!      .                                              .
!      ................................................

!  Arguments:
!  =========
!
!   map     see Subroutine QPP_initialize
!   prob    see Subroutine QPP_reorder
!   X       REAL, OPTIONAL. If present, returns the value of x
!   Y       REAL, OPTIONAL. If present, returns the value of y
!   Z       REAL, OPTIONAL. If present, returns the value of z
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      TYPE ( QPP_map_type ), INTENT( IN ) :: map
      TYPE ( QPP_inform_type ), INTENT( OUT ) :: inform
      TYPE ( QPT_problem_type ), INTENT( IN ) :: prob

      REAL ( KIND = wp ), OPTIONAL, INTENT( OUT ), DIMENSION( map%n ) :: X_val
      REAL ( KIND = wp ), OPTIONAL, INTENT( OUT ), DIMENSION( map%n ) :: Z_val
      REAL ( KIND = wp ), OPTIONAL, INTENT( OUT ), DIMENSION( map%m ) :: Y_val

!  Check to see that the mapping arrays have been set

      IF ( .NOT. map%set ) THEN
        inform%status = GALAHAD_error_unordered
        RETURN
      END IF

!  Recover the appropriate array(s)

      IF ( PRESENT( X_val ) ) X_val( : map%n ) = prob%X( map%x_map( : map%n ) )
      IF ( PRESENT( Y_val ) ) THEN
        IF ( ALLOCATED( prob%Y ) ) THEN
          Y_val( : map%m ) = prob%Y( map%c_map( : map%m ) )
        ELSE IF ( ALLOCATED( prob%Y_l ) .AND. ALLOCATED( prob%Y_u ) ) THEN
          Y_val( : map%m ) = prob%Y_l( map%c_map( : map%m ) ) +                &
                             prob%Y_u( map%c_map( : map%m ) ) 
        END IF
      END IF
      IF ( PRESENT( Z_val ) ) THEN
        IF ( ALLOCATED( prob%Z ) ) THEN
          Z_val( : map%n ) = prob%Z( map%x_map( : map%n ) )
        ELSE IF ( ALLOCATED( prob%Z_l ) .AND. ALLOCATED( prob%Z_u ) ) THEN
          Z_val( : map%n ) = prob%Z_l( map%x_map( : map%n ) ) +                &
                             prob%Z_u( map%x_map( : map%n ) ) 
        END IF
      END IF
      inform%status = GALAHAD_ok
      RETURN  

!  End of QPP_get_values

      END SUBROUTINE QPP_get_values

!-*-*-*-*-*-*-   Q P P _ r e s t o r e  S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE QPP_restore( map, inform, prob, get_all,                      &
                              get_all_parametric, get_f, get_g,                &
                              get_dg, get_x, get_y, get_z, get_x_bounds,       &
                              get_dx_bounds, get_c, get_c_bounds,              &
                              get_dc_bounds, get_A, get_H )

!      ..............................................
!      .                                            .
!      .  Apply the inverse of the permutations     .
!      .  computed in QPP_reorder to the            .
!      .  data prob as appropriate                  .
!      .                                            .
!      ..............................................

!  Arguments:
!  =========
!
!   prob    see Subroutine QPP_reorder
!   map     see Subroutine QPP_initialize
!
!  inform is a structure of type QPP_inform_type that provides 
!    information on exit from QPP_apply. The component status 
!    has possible values:
!  
!     0 Normal termination.
!
!     1 The mapping arrays have not yet been set. Either QPP_reorder
!      has not yet been called, or the call was unsuccessful.
!
!   get_all       LOGICAL, OPTIONAL. If present, process the entire problem
!   get_all_parametric    LOGICAL, OPTIONAL. If present, process the entire 
!                 problem including parametric parts
!   get_f         LOGICAL, OPTIONAL. If present, process f
!   get_g         LOGICAL, OPTIONAL. If present, process g
!   get_dg        LOGICAL, OPTIONAL. If present, process dg
!   get_x         LOGICAL, OPTIONAL. If present, process x
!   get_y         LOGICAL, OPTIONAL. If present, process y
!   get_z         LOGICAL, OPTIONAL. If present, process z
!   get_x_bounds  LOGICAL, OPTIONAL. If present, process x_l and x_u
!   get_dx_bounds LOGICAL, OPTIONAL. If present, process dx_l and dx_u
!   get_c         LOGICAL, OPTIONAL. If present, process c
!   get_c_bounds  LOGICAL, OPTIONAL. If present, process c_l and c_u
!   get_dc_bounds LOGICAL, OPTIONAL. If present, process dc_l and dc_u
!   get_A         LOGICAL, OPTIONAL. If present, process A
!   get_H         LOGICAL, OPTIONAL. If present, process H or weight/x0

!  Dummy arguments

      TYPE ( QPP_map_type ), INTENT( INOUT ) :: map
      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      TYPE ( QPP_inform_type ), INTENT( OUT ) :: inform
      LOGICAL, OPTIONAL, INTENT( IN ) :: get_all, get_g, get_dg,               &
                                         get_c, get_x_bounds, get_c_bounds,    &
                                         get_dx_bounds, get_dc_bounds,         &
                                         get_x, get_y, get_z, get_A, get_H,    &
                                         get_all_parametric, get_f

!  Local variables

      INTEGER :: i, j, k, l
      LOGICAL :: apy, apyl, apyu, apz, apzl, apzu

!  Check that Y or Y_l/Y_u and Z or Z_l/Z_u has been allocated

      apy = ALLOCATED( prob%Y )
      apyl = ALLOCATED( prob%Y_l ) ; apyu = ALLOCATED( prob%Y_u )
      IF ( .NOT. ( apy .OR. ( apyl .AND. apyu ) ) ) THEN
        inform%status = GALAHAD_error_y_unallocated ; RETURN
      END IF

      apz = ALLOCATED( prob%Z )
      apzl = ALLOCATED( prob%Z_l ) ; apzu = ALLOCATED( prob%Z_u )
      IF ( .NOT. ( apz .OR. ( apzl .AND. apzu ) ) ) THEN
        inform%status = GALAHAD_error_z_unallocated ; RETURN
      END IF

!  Check to see that the mapping arrays have been set

      IF ( .NOT. map%set ) THEN
        inform%status = GALAHAD_error_unordered
        RETURN
      END IF

!  Check to see if permuted A and H are available 

      IF ( ( ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.      &
               PRESENT( get_c ) .OR. PRESENT( get_A ) )                        &
               .AND. .NOT. map%a_perm ) ) THEN
        inform%status = GALAHAD_error_ah_unordered
        RETURN
      END IF
      IF ( ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.        &
             PRESENT( get_H ) )                                                &
            .AND. prob%Hessian_kind < 0 .AND. .NOT. map%h_perm ) THEN
        inform%status = GALAHAD_error_ah_unordered
        RETURN
      END IF
      IF ( ( PRESENT( get_g ) .OR. PRESENT( get_c ) .OR.                       &
             PRESENT( get_c_bounds ) )                                         &
             .AND. map%n_reordered < map%n                                     &
           .AND. .NOT. ( map%a_perm .AND. map%h_perm ) ) THEN
        inform%status = GALAHAD_error_ah_unordered
        RETURN
      END IF

!  If there are fixed variables, compute suitable dual variables

      IF ( ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.        &
             PRESENT( get_z ) ) .AND. map%n_reordered < map%n ) THEN

!  ... initialize them as g

        IF ( prob%gradient_kind == 0 ) THEN
          IF ( apz ) THEN
            prob%Z( map%n_reordered + 1 : map%n ) = zero
          ELSE
            prob%Z_l( map%n_reordered + 1 : map%n ) = zero
          END IF
        ELSE IF ( prob%gradient_kind == 1 ) THEN
          IF ( apz ) THEN
            prob%Z( map%n_reordered + 1 : map%n ) = one
          ELSE
            prob%Z_l( map%n_reordered + 1 : map%n ) = one
          END IF
        ELSE 
          IF ( apz ) THEN
            prob%Z( map%n_reordered + 1 : map%n ) =                            &
              prob%G( map%n_reordered + 1 : map%n )
          ELSE
            prob%Z_l( map%n_reordered + 1 : map%n ) =                          &
              prob%G( map%n_reordered + 1 : map%n )
          END IF
        END IF

!  ... now add suitable rows of H * x

        IF ( prob%Hessian_kind == 1 ) THEN

!   least-distance Hessian

          IF ( apz ) THEN
            prob%Z( prob%n + 1 : map%n ) = prob%Z( prob%n + 1 : map%n ) +      &
              ( prob%X( prob%n + 1 : map%n ) - prob%X0( prob%n + 1 : map%n ) )
          ELSE
            prob%Z_l( prob%n + 1 : map%n ) = prob%Z_l( prob%n + 1 : map%n ) +  &
              ( prob%X( prob%n + 1 : map%n ) - prob%X0( prob%n + 1 : map%n ) )
          END IF
        ELSE IF ( prob%Hessian_kind > 1 ) THEN

!   weighted-least-distance Hessian

          IF ( apz ) THEN
            prob%Z( prob%n + 1 : map%n ) = prob%Z( prob%n + 1 : map%n ) +      &
             ( prob%X( prob%n + 1 : map%n ) - prob%X0( prob%n + 1 : map%n ) )* &
                prob%WEIGHT( prob%n + 1 : map%n ) ** 2 
          ELSE
            prob%Z_l( prob%n + 1 : map%n ) = prob%Z_l( prob%n + 1 : map%n ) +  &
             ( prob%X( prob%n + 1 : map%n ) - prob%X0( prob%n + 1 : map%n ) )* &
                prob%WEIGHT( prob%n + 1 : map%n ) ** 2 
          END IF
        ELSE IF ( prob%Hessian_kind < 0 ) THEN

!   general Hessian

!  ..... rows with a diagonal entry

          DO i = prob%n + 1, map%h_diag_end_fixed
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 2
              j = prob%H%col( l )
              IF ( apz ) THEN
                IF ( j <= prob%n ) THEN
                  prob%Z( i ) = prob%Z( i ) + prob%H%val( l ) * prob%X( j )
                ELSE
                  prob%Z( i ) = prob%Z( i ) + prob%H%val( l ) * prob%X( j )
                  prob%Z( j ) = prob%Z( j ) + prob%H%val( l ) * prob%X( i )
                END IF
              ELSE
                IF ( j <= prob%n ) THEN
                  prob%Z_l( i ) = prob%Z_l( i ) + prob%H%val( l ) * prob%X( j )
                ELSE
                  prob%Z_l( i ) = prob%Z_l( i ) + prob%H%val( l ) * prob%X( j )
                  prob%Z_l( j ) = prob%Z_l( j ) + prob%H%val( l ) * prob%X( i )
                END IF
              END IF
            END DO
            l = prob%H%ptr( i + 1 ) - 1
            IF ( apz ) THEN
              prob%Z( i ) = prob%Z( i ) + prob%H%val( l ) * prob%X( i )
            ELSE
              prob%Z_l( i ) = prob%Z_l( i ) + prob%H%val( l ) * prob%X( i )
            END IF
          END DO

!  ..... rows without a diagonal entry

          DO i = map%h_diag_end_fixed + 1, map%n
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              j = prob%H%col( l )
              IF ( apz ) THEN
                IF ( j <= prob%n ) THEN
                  prob%Z( i ) = prob%Z( i ) + prob%H%val( l ) * prob%X( j )
                ELSE
                  prob%Z( i ) = prob%Z( i ) + prob%H%val( l ) * prob%X( j )
                  prob%Z( j ) = prob%Z( j ) + prob%H%val( l ) * prob%X( i )
                END IF
              ELSE
                IF ( j <= prob%n ) THEN
                  prob%Z_l( i ) = prob%Z_l( i ) + prob%H%val( l ) * prob%X( j )
                ELSE
                  prob%Z_l( i ) = prob%Z_l( i ) + prob%H%val( l ) * prob%X( j )
                  prob%Z_l( j ) = prob%Z_l( j ) + prob%H%val( l ) * prob%X( i )
                END IF
              END IF
            END DO
          END DO
        END IF

!  ... finally subtract suitable rows of A^T * y

        DO i = prob%n + 1, map%n
          DO l = map%ptr_a_fixed( i ), map%ptr_a_fixed( i + 1 ) - 1
            j = prob%A%col( l )
            IF ( apz ) THEN
              IF ( apy ) THEN
                prob%Z( i ) = prob%Z( i ) - prob%A%val( l ) * prob%Y( j )
              ELSE
                prob%Z( i ) = prob%Z( i ) - prob%A%val( l ) *                  &
                  ( prob%Y_l( j ) + prob%Y_u( j ) )
              END IF
            ELSE
              IF ( apy ) THEN
                prob%Z_l( i ) = prob%Z_l( i ) - prob%A%val( l ) * prob%Y( j )
              ELSE
                prob%Z_l( i ) = prob%Z_l( i ) - prob%A%val( l ) *              &
                  ( prob%Y_l( j ) + prob%Y_u( j ) )
              END IF
            END IF
          END DO
        END DO

!  Copy to Z_l and Z_u if required

        IF ( apz ) THEN
          IF ( apzl .AND. apzu ) THEN
            DO i = map%n_reordered + 1, map%n 
              IF ( prob%Z( i ) > 0 ) THEN
                prob%Z_l( i ) = prob%Z( i )
                prob%Z_u( i ) = zero
              ELSE
                prob%Z_u( i ) = prob%Z( i )
                prob%Z_l( i ) = zero
              END IF
            END DO
          END IF
        ELSE
          IF ( apzl .AND. apzu ) THEN
            DO i = map%n_reordered + 1, map%n 
              IF ( prob%Z_l( i ) > 0 ) THEN
                prob%Z_u( i ) = zero
              ELSE
                prob%Z_u( i ) = prob%Z_l( i )
                prob%Z_l( i ) = zero
              END IF
            END DO
          END IF
        END IF
      END IF

!  Form c = A * x

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_c ) ) THEN
        prob%C( : map%m ) = zero
        CALL QPP_AX( map, prob%X, prob%A%val, prob%A%col, prob%A%ptr,          &
                      map%m, prob%C( : map%m ) )
      END IF

!  Transform f, g and the bounds on the constraints to account for 
!  fixed variables

      IF ( ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) )           &
             .AND. map%n_reordered < map%n ) THEN
          CALL QPP_add_fixed( map, prob,                                       &
                              f =.TRUE., g =.TRUE.,                            &
                              c = .TRUE., c_bounds = .TRUE. )
      ELSE IF ( ( PRESENT( get_f ) .OR. PRESENT( get_g ) .OR.                  &
                  PRESENT( get_c ) .OR. PRESENT( get_c_bounds ) )              &
                 .AND. map%n_reordered < map%n ) THEN
        CALL QPP_add_fixed( map, prob,                                         &
                            f = PRESENT( get_f ), g = PRESENT( get_g ),        &
                            c =  PRESENT( get_c ),                             &
                            c_bounds = PRESENT( get_c_bounds ) )
      END IF

!  See if we need to invert the mappings

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_H ) .OR. PRESENT( get_A ) )                            &
        CALL QPP_invert_mapping( map%n, map%x_map, map%IW( : map%n ) )
      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_A ) )                                                  &
        CALL QPP_invert_mapping( map%m, map%c_map, map%IW( : map%m ) )

!  Restore H

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_H ) ) THEN

!  Weighted-least-distance Hessian

        IF ( prob%Hessian_kind > 1 ) THEN
          CALL SORT_inverse_permute( map%n, map%x_map, X = prob%WEIGHT )
          CALL SORT_inverse_permute( map%n, map%x_map, X = prob%X0 )
        ELSE IF ( prob%Hessian_kind > 0 ) THEN
          CALL SORT_inverse_permute( map%n, map%x_map, X = prob%X0 )
 
!  General Hessian

        ELSE IF ( prob%Hessian_kind < 0 ) THEN  

!  Original dense storage

          IF ( map%h_type == 2 .OR. map%h_type == 3 ) THEN
            CALL SORT_inplace_permute( map%h_ne, map%h_map_inverse,            &
                                       X = prob%H%val )

!  Original row-wise storage
   
          ELSE IF ( map%h_type == 1 ) THEN

!  Compute the number of entries in each row of H, and renumber its columns

            map%IW( : map%n ) = 0
            DO k = 1, map%n
              i = map%x_map( k )
              DO l = prob%H%ptr( k ), prob%H%ptr( k + 1 ) - 1
                j = map%x_map( prob%H%col( l ) )
                IF ( i >= j ) THEN
                  map%IW( i ) = map%IW( i ) + 1
                  prob%H%col( l ) = j
                ELSE
                  map%IW( j ) = map%IW( j ) + 1
                  prob%H%col( l ) = i
                END IF
              END DO
            END DO

!  Set the starting addresses for each row in the restored matrix

            prob%H%ptr( 1 ) = 1
            DO i = 1, map%n
              prob%H%ptr( i + 1 ) = prob%H%ptr( i ) + map%IW( i )
            END DO

!  Undo the reordering of H

            CALL SORT_inplace_permute( map%h_ne, map%h_map_inverse,            &
                                       X = prob%H%val, IX = prob%H%col )

!  Original co-ordinate storage

          ELSE

!  Renumber the rows and columns

            DO k = 1, map%n
              i = map%x_map( k )
              DO l = prob%H%ptr( k ), prob%H%ptr( k + 1 ) - 1
                j = map%x_map( prob%H%col( l ) )
                IF ( i >= j ) THEN
                  prob%H%row( l ) = i
                  prob%H%col( l ) = j
                ELSE
                  prob%H%row( l ) = j
                  prob%H%col( l ) = i
                END IF
              END DO
            END DO

!  Undo the reordering of H

            CALL SORT_inplace_permute( map%h_ne, map%h_map_inverse,            &
                X = prob%H%val, IX = prob%H%row, IY = prob%H%col )
          END IF

!  Record that H had been restored

          map%h_perm = .FALSE.
        END IF
      END IF

!  Restore A

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_A ) ) THEN

!  Original dense storage

        IF ( map%a_type == 2 ) THEN
          CALL SORT_inplace_permute(                                           &
              map%a_ne, map%a_map_inverse, X = prob%A%val )

!  Original row-wise storage
  
        ELSE IF ( map%a_type == 1 ) THEN

!  Compute the number of entries in each row of A, and renumber its columns.
!  NB. Any columns corresponding to FIXED variables, will have been
!  moved to the end of A, and will be stored by COLUMN not by row. In 
!  particular, A_col for these entries gives the row and not the column
!  number

          map%IW( : map%m ) = 0
          DO k = 1, map%m
            i = map%c_map( k )
            DO l = prob%A%ptr( k ), prob%A%ptr( k + 1 ) - 1
              map%IW( i ) = map%IW( i ) + 1
              prob%A%col( l ) = map%x_map( prob%A%col( l ) )
            END DO
          END DO

          DO l = map%n_reordered + 1, map%n
            j = map%x_map( l )
            DO k = map%ptr_a_fixed( l ), map%ptr_a_fixed( l + 1 ) - 1
              i =  map%c_map( prob%A%col( k ) )
              map%IW( i ) = map%IW( i ) + 1
              prob%A%col( k ) = j
            END DO
          END DO

!  Set the starting addresses for each row in the permuted matrix

          prob%A%ptr( 1 ) = 1
          DO i = 1, map%m
            prob%A%ptr( i + 1 ) = prob%A%ptr( i ) + map%IW( i )
          END DO
    
!  Undo the reordering of A

          CALL SORT_inplace_permute( map%a_ne, map%a_map_inverse,              &
                                     X = prob%A%val, IX = prob%A%col )

!  Original co-ordinate storage

        ELSE 

!  Renumber the rows and columns

          DO k = 1, map%m
            i = map%c_map( k )
            DO l = prob%A%ptr( k ), prob%A%ptr( k + 1 ) - 1
              prob%A%row( l ) = i
              prob%A%col( l ) = map%x_map( prob%A%col( l ) )
            END DO
          END DO

          DO l = map%n_reordered + 1, map%n
            j = map%x_map( l )
            DO k = map%ptr_a_fixed( l ), map%ptr_a_fixed( l + 1 ) - 1
              prob%A%row( k ) = map%c_map( prob%A%col( k ) )
              prob%A%col( k ) = j
            END DO
          END DO

!  Undo the reordering of A

          CALL SORT_inplace_permute( map%a_ne, map%a_map_inverse,              &
              X = prob%A%val, IX = prob%A%row, IY = prob%A%col )
        END IF

!  Record that A had been restored

        map%a_perm = .FALSE.
      END IF

!  See if we need to reinvert the mappings

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_H ) .OR. PRESENT( get_A ) )                            &
        CALL QPP_invert_mapping( map%n, map%x_map, map%IW( : map%n ) )
      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_A ) )                                                  &
        CALL QPP_invert_mapping( map%m, map%c_map, map%IW( : map%m ) )

!  Restore g

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_g ) ) THEN
        IF ( prob%gradient_kind /= 0 .AND. prob%gradient_kind /= 1 )           &
          CALL SORT_inverse_permute( map%n, map%x_map, X = prob%G )
      END IF

!  Restore dg

      IF ( PRESENT( get_all_parametric ) .OR. PRESENT( get_dg ) ) THEN
        IF ( ALLOCATED( prob%DG ) ) THEN
          IF ( SIZE( prob%DG ) == map%n )                                      &
            CALL SORT_inverse_permute( map%n, map%x_map, X = prob%DG )
        END IF
      END IF

!  Restore c

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_c ) ) THEN
        CALL SORT_inverse_permute( map%m, map%c_map, X = prob%C )
      END IF

!  Restore the bounds on c

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_c_bounds ) ) THEN
        CALL SORT_inverse_permute( map%m, map%c_map, X = prob%C_l )
        CALL SORT_inverse_permute( map%m, map%c_map, X = prob%C_u )
      END IF

!  Restore the bounds on dc

      IF ( PRESENT( get_all_parametric ) .OR. PRESENT( get_dc_bounds ) ) THEN
        IF ( ALLOCATED( prob%DC_l ) ) THEN
          IF ( SIZE( prob%DC_l ) == map%m )                                    &
            CALL SORT_inverse_permute( map%m, map%c_map, X = prob%DC_l )
        END IF
        IF ( ALLOCATED( prob%DC_u ) ) THEN
          IF ( SIZE( prob%DC_u ) == map%m )                                    &
            CALL SORT_inverse_permute( map%m, map%c_map, X = prob%DC_u )
        END IF
      END IF

!  Restore x

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_x ) ) THEN
        CALL SORT_inverse_permute( map%n, map%x_map, X = prob%X )
      END IF

!  Restore the bounds on x

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_x_bounds ) ) THEN
        CALL SORT_inverse_permute( map%n, map%x_map, X = prob%X_l )
        CALL SORT_inverse_permute( map%n, map%x_map, X = prob%X_u )
      END IF

!  Restore the bounds on dx

      IF ( PRESENT( get_all_parametric ) .OR. PRESENT( get_dx_bounds ) ) THEN
        IF ( ALLOCATED( prob%DX_l ) ) THEN
          IF ( SIZE( prob%DX_l ) == map%n )                                    &
            CALL SORT_inverse_permute( map%n, map%x_map, X = prob%DX_l )
        END IF
        IF ( ALLOCATED( prob%DX_u ) ) THEN
          IF ( SIZE( prob%DX_u ) == map%n )                                    &
            CALL SORT_inverse_permute( map%n, map%x_map, X = prob%DX_u )
        END IF
      END IF

!  Restore the dual variables z

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_z ) ) THEN
        IF ( apz ) CALL SORT_inverse_permute( map%n, map%x_map, X = prob%Z )
        IF ( apzl ) CALL SORT_inverse_permute( map%n, map%x_map, X = prob%Z_l )
        IF ( apzu ) CALL SORT_inverse_permute( map%n, map%x_map, X = prob%Z_u )
      END IF

!  Restore the Lagrange multipliers y

      IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
           PRESENT( get_y ) ) THEN
        IF ( apy ) CALL SORT_inverse_permute( map%m, map%c_map, X = prob%Y )
        IF ( apyl ) CALL SORT_inverse_permute( map%m, map%c_map, X = prob%Y_l )
        IF ( apyu ) CALL SORT_inverse_permute( map%m, map%c_map, X = prob%Y_u )
      END IF

!  Pick up the correct dimensions

      prob%n = map%n
      prob%m = map%m
!     IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
!          PRESENT( get_A ) ) prob%A%ne = map%a_ne_original
!     IF ( PRESENT( get_all ) .OR. PRESENT( get_all_parametric ) .OR.          &
!          PRESENT( get_H ) ) prob%H%ne = map%h_ne_original
      IF ( prob%Hessian_kind < 0 ) THEN
        CALL QPT_put_H( prob%H%type, SMT_get( map%h_type_original ) )
        IF ( SMT_get( prob%H%type ) == 'COORDINATE' )                          &
           prob%H%ne = map%h_ne_original
      END IF
      CALL QPT_put_A( prob%A%type, SMT_get( map%a_type_original ) )
      IF ( SMT_get( prob%A%type ) == 'COORDINATE' )                            &
         prob%A%ne = map%a_ne_original

      inform%status = GALAHAD_ok
      RETURN  

!  End of QPP_restore

      END SUBROUTINE QPP_restore

!-*-*-*-*-   Q P P _ f i n a l i z e   S U B R O U T I N E  -*-*-*

      SUBROUTINE QPP_terminate( map, control, inform )

!      ..............................................
!      .                                            .
!      .  Deallocate internal arrays at the end     .
!      .  of the computation                        .
!      .                                            .
!      ..............................................

!  Arguments:
!  =========
!
!   map     see Subroutine QPP_initialize
!   control see Subroutine QPP_initialize
!   inform  see Subroutine QPP_reorder

!  Dummy arguments

      TYPE ( QPP_map_type ), INTENT( INOUT ) :: map
      TYPE ( QPP_control_type ), INTENT( IN ) :: control        
      TYPE ( QPP_inform_type ), INTENT( INOUT ) :: inform

      inform%status = GALAHAD_ok
      map%set = .FALSE.

!  Deallocate all mapping arrays

      IF ( ALLOCATED( map%IW ) ) THEN
        DEALLOCATE( map%IW, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_allocate
          IF ( control%error > 0 )                                             &
              WRITE( control%error, 2900 ) 'map%IW', inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( map%x_map ) ) THEN
        DEALLOCATE( map%x_map, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_allocate
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) 'map%x_map', inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( map%c_map ) ) THEN
        DEALLOCATE( map%c_map, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_allocate
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) 'map%c_map', inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( map%h_map_inverse ) ) THEN
        DEALLOCATE( map%h_map_inverse, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_allocate
          IF ( control%error > 0 )                                             &
           WRITE( control%error, 2900 ) 'map%h_map_inverse', inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( map%a_map_inverse ) ) THEN
        DEALLOCATE( map%a_map_inverse, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_allocate
          IF ( control%error > 0 )                                             &
           WRITE( control%error, 2900 ) 'map%a_map_inverse', inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( map%ptr_a_fixed ) ) THEN
        DEALLOCATE( map%ptr_a_fixed, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_allocate
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) 'map%ptr_a_fixed', inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( map%a_type_original ) ) THEN
        DEALLOCATE( map%a_type_original, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_allocate
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) 'map%a_type_original',                &
              inform%alloc_status
        END IF
      END IF

      IF ( ALLOCATED( map%h_type_original ) ) THEN
        DEALLOCATE( map%h_type_original, STAT = inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_allocate
          IF ( control%error > 0 )                                             &
            WRITE( control%error, 2900 ) 'map%h_type_original',                &
              inform%alloc_status
        END IF
      END IF

!  Non-executable statement

 2900 FORMAT( ' ** Message from -QPP_terminate-', /,                          &
              ' Deallocation error, for data%', A20, /, ' status = ', I6 ) 

!  End of subroutine QPP_terminate

      END SUBROUTINE QPP_terminate

!-*-*-   Q P P _ o r d e r _ r o w s   S U B R O U T I N E  -*-*-*

      SUBROUTINE QPP_order_rows( n_rows, VAL, COL, PTR, MAP_inverse)

!  Reorder the entries in each row so that the column entries appear
!  in order of increasing index; also update (the inverse of) a map from
!  a previous reordering

!  Arguments:
!  =========
!
!   n_rows      number of rows
!   PTR         starting addresses of each row, as well as 1 beyond the last row
!   VAL         values of entries in each row
!   COL         column indices corresponding to values
!   map         a permutation of the numbers 1, .., Ptr( n_rows + 1 ) - 1
!   map_inverse the inverse permutation of the above

!  Dummy arguments

      INTEGER, INTENT( IN ) :: n_rows
      INTEGER, INTENT( IN ), DIMENSION( n_rows + 1 ) :: PTR
      INTEGER, INTENT( INOUT ),                                                &
               DIMENSION( Ptr( n_rows + 1 ) - 1 ) :: COL, MAP_inverse
      REAL ( KIND = wp ), INTENT( INOUT ),                                     &
               DIMENSION( Ptr( n_rows + 1 ) - 1 ) :: VAL

!  Local variables

      INTEGER :: i, current, col_current, inverse, inverse_current
      INTEGER :: previous, next, row_start, row_end, inrow, inform_quicksort
      REAL ( KIND = wp ) :: val_current
      INTEGER, PARAMETER :: do_quicksort = 10

!  loop over the rows

      DO i = 1, n_rows
        row_start = PTR( i )
        row_end = PTR( i + 1 ) - 1

!  For the current row

        inrow = row_end - row_start + 1
        IF ( inrow <= 0 ) CYCLE
        IF ( inrow > do_quicksort ) THEN

!  If the row contains enough entries, use quicksort ...

          DO current = row_start + 1, row_end

!  Skip if the value in the current position is already in order

            IF ( COL( current ) < COL( current - 1 ) ) THEN
              CALL SORT_quicksort( inrow, COL( row_start : row_end ),          &
                                   inform_quicksort,                           &
                                   MAP_inverse( row_start : row_end ),         &
                                   VAL( row_start : row_end ) )
          
              EXIT
            END IF
          END DO
        ELSE

!  ... else use a simple entry shuffle

          DO current = row_start + 1, row_end
            col_current = COL( current )

!  Skip if the value in the current position is already in order

            IF ( col_current < COL( current - 1 ) ) THEN
              val_current = VAL( current )
              inverse_current = MAP_inverse( current )

!  The value in the current position is out of order, but those in
!  positions row_start, current - 1 are now in order

              DO previous = row_start, current - 1
                IF ( col_current < COL( previous ) ) THEN

!  The current value should be inserted at position previous

                  DO next = current, previous + 1, - 1

!  Shift values previous, ... , current + 1 one place to the right

                    COL( next ) = COL( next - 1 )
                    VAL( next ) = VAL( next - 1 )

!  Update the inverse map

                    inverse = MAP_inverse( next - 1 )
                    MAP_inverse( next ) = inverse

                  END DO

!  Insert the current value in its correct position

                  VAL( previous ) = val_current
                  COL( previous ) = col_current

!  Update the inverse map

                  MAP_inverse( previous ) = inverse_current
                  EXIT

                END IF
              END DO
            END IF
          END DO
        END IF
!       write(45, "( ' after  ', /, ( 10I6 ) )" ) COL( row_start : row_end )
!       write(47, "( '        ', /, ( 10I6 ) )" ) COL( row_start : row_end )
      END DO

      RETURN

!  End of subroutine QPP_order_rows

      END SUBROUTINE QPP_order_rows

!-   Q P P _ i n p l a c e _ p e r m u t e  S U B R O U T I N E  -

      SUBROUTINE QPP_inplace_permute( n, MAP, X, IX, IY )

!  Permute the entries of X so that x_i appears in position map_i
!  Do this without resorting to extra vector storage. Optionally,
!  permute the entries of IX and IY so that ixc_i and iy_i 
!  appear in positions map_i

!  Arguments:
!  =========
!
!   n           number of components of x
!   X           the array x
!   MAP         the permutation map
!   IX          the array ix
!   IY          the array iy

!  Dummy arguments

      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( INOUT ), DIMENSION( n ) :: MAP
      INTEGER, INTENT( INOUT ), OPTIONAL, DIMENSION( n ) :: IX, IY
      REAL ( KIND = wp ), INTENT( INOUT ), OPTIONAL, DIMENSION( n ) :: X

!  Local variables

      INTEGER :: i, mi, mi_old, iymi, iymi_old, ixmi, ixmi_old
      REAL ( KIND = wp ) :: xmi, xmi_old

!  For X, IX and IY:

      IF ( PRESENT( IX ) .AND. PRESENT( IY ) ) THEN

!  loop over the entries of X, IX and IY

        DO i = 1, n
          mi = MAP( i )

!  Skip any entry which is already in place

          IF ( mi == i ) THEN
            CYCLE

!  Skip any entry which has already been moved into place, remembering
!  to un-negate the relevant entry in MAP

          ELSE IF ( mi < 0 ) THEN
            MAP( i ) = - mi

!  The i-th entry is not in place. Chase through the list of entries
!  i, MAP( i ), MAP( MAP( i ) ), ... until MAP( ... ( MAP( i ) ) ... ) = i
!  moving entries into place. Negate the relevant entries in MAP so that
!  these entries will not be moved again

          ELSE 
            xmi_old = X( i )
            iymi_old = IY( i )
            ixmi_old = IX( i )
            DO 
              xmi = X( mi )
              iymi = IY( mi )
              ixmi = IX( mi )
              X( mi ) = xmi_old
              IY( mi ) = iymi_old
              IX( mi ) = ixmi_old
              xmi_old = xmi
              iymi_old = iymi
              ixmi_old = ixmi
              mi_old = mi
              mi = MAP( mi_old )
              MAP( mi_old ) = - mi
              IF ( mi == i ) EXIT
            END DO
            X( i ) = xmi_old
            IY( i ) = iymi_old
            IX( i ) = ixmi_old
          END IF
        END DO

!  For X and IX:

      ELSE IF ( PRESENT( IX ) ) THEN

!  loop over the entries of X and IX

        DO i = 1, n
          mi = MAP( i )

!  Skip any entry which is already in place

          IF ( mi == i ) THEN
            CYCLE

!  Skip any entry which has already been moved into place, remembering
!  to un-negate the relevant entry in MAP

          ELSE IF ( mi < 0 ) THEN
            MAP( i ) = - mi

!  The i-th entry is not in place. Chase through the list of entries
!  i, MAP( i ), MAP( MAP( i ) ), ... until MAP( ... ( MAP( i ) ) ... ) = i
!  moving entries into place. Negate the relevant entries in MAP so that
!  these entries will not be moved again

          ELSE 
            xmi_old = X( i )
            ixmi_old = IX( i )
            DO 
              xmi = X( mi )
              ixmi = IX( mi )
              X( mi ) = xmi_old
              IX( mi ) = ixmi_old
              xmi_old = xmi
              ixmi_old = ixmi
              mi_old = mi
              mi = MAP( mi_old )
              MAP( mi_old ) = - mi
              IF ( mi == i ) EXIT
            END DO
            X( i ) = xmi_old
            IX( i ) = ixmi_old
          END IF
        END DO

!  For just X:

      ELSE

!  loop over the entries of X

        DO i = 1, n
          mi = MAP( i )

!  Skip any entry which is already in place

          IF ( mi == i ) THEN
            CYCLE

!  Skip any entry which has already been moved into place, remembering
!  to un-negate the relevant entry in MAP

          ELSE IF ( mi < 0 ) THEN
            MAP( i ) = - mi

!  The i-th entry is not in place. Chase through the list of entries
!  i, MAP( i ), MAP( MAP( i ) ), ... until MAP( ... ( MAP( i ) ) ... ) = i
!  moving entries into place. Negate the relevant entries in MAP so that
!  these entries will not be moved again

          ELSE 
            xmi_old = X( i )
            DO 
              xmi = X( mi )
              X( mi ) = xmi_old
              xmi_old = xmi
              mi_old = mi
              mi = MAP( mi_old )
              MAP( mi_old ) = - mi
              IF ( mi == i ) EXIT
            END DO
            X( i ) = xmi_old
          END IF
        END DO

      END IF

      RETURN

!  End of subroutine SORT_inplace_permute

      END SUBROUTINE QPP_inplace_permute

!-   Q P P _ i n v e r s e _ p e r m u t e  S U B R O U T I N E  -

      SUBROUTINE QPP_inverse_permute( n, MAP_inverse, X, IX )

!  Permute the entries of X so that x(map_i) appears in position i
!  Do this without resorting to extra vector storage. Optionally,
!  permute the entries of IX so that ix(map_i) appears in position i

!  Arguments:
!  =========
!
!   n           number of components of x
!   X           the array x
!   MAP_inverse the permutation map
!   IX          the array IX

!  Dummy arguments

      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( INOUT ), DIMENSION( n ) :: MAP_inverse
      INTEGER, INTENT( INOUT ), OPTIONAL, DIMENSION( n ) :: IX
      REAL ( KIND = wp ), INTENT( INOUT ), OPTIONAL, DIMENSION( n ) :: X

!  Local variables

      INTEGER :: i, mi, mi_old, ixi
      REAL ( KIND = wp ) :: xi

!  For both X and IX:

      IF ( PRESENT( IX ) ) THEN

!  loop over the entries of X and IX

        DO i = 1, n
          mi = MAP_inverse( i )

!  Skip any entry which is already in place

          IF ( mi == i ) THEN
            CYCLE

!  Skip any entry which has already been moved into place, remembering
!  to un-negate the relevant entry in MAP_inverse

          ELSE IF ( mi < 0 ) THEN
            MAP_inverse( i ) = - mi

!  The i-th entry is not in place. Chase through the list of entries
!  i, MAP_inverse( i ), MAP_inverse( MAP_inverse( i ) ), ... until 
!  MAP_inverse( ... ( MAP_inverse( i ) ) ... ) = i, moving entries into place. 
!  Negate the relevant entries in MAP_inverse so that these entries will 
!  not be moved again

          ELSE 
            xi = X( i )
            ixi = IX( i )
            mi_old = i
            DO 
              X( mi_old ) = X( mi )
              IX( mi_old ) = IX( mi )
              mi_old = mi
              mi = MAP_inverse( mi_old )
              MAP_inverse( mi_old ) = - mi
              IF ( mi == i ) EXIT
            END DO
            X( mi_old ) = xi
            IX( mi_old ) = ixi
          END IF
        END DO

!  For just X:

      ELSE

!  loop over the entries of X

        DO i = 1, n
          mi = MAP_inverse( i )

!  Skip any entry which is already in place

          IF ( mi == i ) THEN
            CYCLE

!  Skip any entry which has already been moved into place, remembering
!  to un-negate the relevant entry in MAP_inverse

          ELSE IF ( mi < 0 ) THEN
            MAP_inverse( i ) = - mi

!  The i-th entry is not in place. Chase through the list of entries
!  i, MAP_inverse( i ), MAP_inverse( MAP_inverse( i ) ), ... until 
!  MAP_inverse( ... ( MAP_inverse( i ) ) ... ) = i, moving entries into place. 
!  Negate the relevant entries in MAP_inverse so that these entries will 
!  not be moved again

          ELSE 
            xi = X( i )
            mi_old = i
            DO 
              X( mi_old ) = X( mi )
              mi_old = mi
              mi = MAP_inverse( mi_old )
              MAP_inverse( mi_old ) = - mi
              IF ( mi == i ) EXIT
            END DO
            X( mi_old ) = xi
          END IF
        END DO

      END IF

      RETURN

!  End of subroutine QPP_inverse_permute

      END SUBROUTINE QPP_inverse_permute

!-*-*-*-   Q P P _ i n v e r t _ m a p p i n g  S U B R O U T I N E  -*-*-*-

      SUBROUTINE QPP_invert_mapping( n, MAP, MAP_inverse )

!  Place the inverse mapping to MAP in MAP_inverse, and then use this
!  to overwrite MAP

!  Arguments:
!  =========
!
!   n           number of components of x
!   MAP         the mapping (and, on exit, its inverse)
!   MAP_inverse its inverse

!  Dummy arguments

      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( INOUT ), DIMENSION( n ) :: MAP, MAP_inverse

!  Local variables

      INTEGER :: i

!  Invert the mapping, MAP

      DO i = 1, n
        MAP_inverse( MAP( i ) ) = i
      END DO

!  Copy this back to MAP

      MAP = MAP_inverse

      RETURN

!  End of subroutine QPP_invert_mapping

      END SUBROUTINE QPP_invert_mapping

!-*-   Q P P _ r e m o v e _ f i x e d   S U B R O U T I N E   -*

      SUBROUTINE QPP_remove_fixed( map, prob, f, g, c_bounds )

!      ................................................
!      .                                              .
!      .  Transform f, g and the bounds on the        .
!      .  constraints to account for fixed variables  .
!      .                                              .
!      ................................................

!  Arguments:
!  =========
!
!   map      see Subroutine QPP_initialize
!   prob%    see Subroutine QPP_reorder
!   f        LOGICAL. If true, adjust f to account for fixed variables via
!              f -> f + <g_x,x_x> + 1/2 <x_x,H_{xx} x_x> 
!            where x = (x_r,x_x ) (etc)
!   g        LOGICAL. If true, adjust g to account for fixed variables via
!              g_r -> g_r + H_{xr} x_x
!   c_bounds LOGICAL. If true, adjust c_l and c_u to account for fixed 
!            variables via 
!              c_l <- c_l - A_x x_x and c_u <- c_u - A_x x_x

!  Dummy arguments

      TYPE ( QPP_map_type ), INTENT( IN ) :: map
      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      LOGICAL, OPTIONAL, INTENT( IN ) :: f, g, c_bounds

!  Local variables

      INTEGER :: i, j, l
      REAL ( KIND = wp ) :: x, c
      LOGICAL :: yes_f, yes_g, yes_c_bounds

      IF ( prob%n >= map%n ) RETURN

      IF ( PRESENT ( f ) ) THEN
        yes_f = f
      ELSE
        yes_f = .FALSE.
      END IF

      IF ( PRESENT ( g ) ) THEN
        yes_g = g
      ELSE
        yes_g = .FALSE.
      END IF

      IF ( PRESENT ( c_bounds ) ) THEN
        yes_c_bounds = c_bounds
      ELSE
        yes_c_bounds = .FALSE.
      END IF

      IF ( yes_f ) THEN
        IF ( prob%gradient_kind == 1 ) THEN
          prob%f = prob%f + SUM( prob%X( prob%n + 1 : map%n ) )
        ELSE IF ( prob%gradient_kind /= 0 ) THEN
          prob%f = prob%f + DOT_PRODUCT( prob%G( prob%n + 1 : map%n ),         &
                                         prob%X( prob%n + 1 : map%n ) )
        END IF
      END IF

!  Process f and g together

      IF ( yes_f .AND. yes_g ) THEN

!  Least-distance Hessian
        
        IF ( prob%Hessian_kind == 1 ) THEN
          DO i = prob%n + 1, map%n
            prob%f = prob%f + half * ( prob%X( i ) - prob%X0( i ) ) ** 2 
          END DO

!  Weighted-least-distance Hessian

        ELSE IF ( prob%Hessian_kind > 1 ) THEN
          DO i = prob%n + 1, map%n
            prob%f = prob%f + half *                                           &
              ( prob%WEIGHT( i ) * ( prob%X( i ) - prob%X0( i ) ) ) ** 2 
          END DO

!  General Hessian

        ELSE IF ( prob%Hessian_kind /= 0 ) THEN

!  rows with a diagonal entry

          DO i = prob%n + 1, map%h_diag_end_fixed
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 2
              j = prob%H%col( l )
              IF ( j <= prob%n ) THEN
                prob%G( j ) = prob%G( j ) + prob%H%val( l ) * prob%X( i )
              ELSE
                prob%f = prob%f + prob%X( j ) * prob%H%val( l ) * prob%X( i ) 
              END IF
            END DO
            l = prob%H%ptr( i + 1 ) - 1
            prob%f = prob%f + half * prob%X( i ) * prob%H%val( l ) * prob%X( i ) 
          END DO

!  rows without a diagonal entry

          DO i = map%h_diag_end_fixed + 1, map%n
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              j = prob%H%col( l )
              IF ( j <= prob%n ) THEN
                prob%G( j ) = prob%G( j ) + prob%H%val( l ) * prob%X( i )
              ELSE
                prob%f = prob%f + prob%X( j ) * prob%H%val( l ) * prob%X( i ) 
              END IF
            END DO
          END DO
        END IF

!  Process g separately

      ELSE IF ( yes_g ) THEN
        
        IF ( prob%Hessian_kind < 0 ) THEN

!  rows with a diagonal entry

          DO i = prob%n + 1, map%h_diag_end_fixed
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 2
              j = prob%H%col( l )
              IF ( j <= prob%n ) THEN
                prob%G( j ) = prob%G( j ) + prob%H%val( l ) * prob%X( i )
              END IF
            END DO
          END DO

!  rows without a diagonal entry

          DO i = map%h_diag_end_fixed + 1, map%n
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              j = prob%H%col( l )
              IF ( j <= prob%n ) THEN
                prob%G( j ) = prob%G( j ) + prob%H%val( l ) * prob%X( i )
              END IF
            END DO
          END DO
        END IF

!  Process f separately

      ELSE IF ( yes_f ) THEN

!  Least-distance Hessian
        
        IF ( prob%Hessian_kind == 1 ) THEN
          DO i = prob%n + 1, map%n
            prob%f = prob%f + half * ( prob%X( i ) - prob%X0( i ) ) ** 2 
          END DO

!  Weighted-least-distance Hessian

        ELSE IF ( prob%Hessian_kind > 1 ) THEN
          DO i = prob%n + 1, map%n
            prob%f = prob%f + half *                                           &
              ( prob%WEIGHT( i ) * ( prob%X( i ) - prob%X0( i ) ) ) ** 2 
          END DO

!  General Hessian

        ELSE IF ( prob%Hessian_kind /= 0 ) THEN

!  rows with a diagonal entry

          DO i = prob%n + 1, map%h_diag_end_fixed
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 2
              j = prob%H%col( l )
              IF ( j > prob%n ) THEN
                prob%f = prob%f + prob%X( j ) * prob%H%val( l ) * prob%X( i ) 
              END IF
            END DO
            l = prob%H%ptr( i + 1 ) - 1
            prob%f = prob%f + half * prob%X( i ) * prob%H%val( l ) * prob%X( i ) 
          END DO

!  rows without a diagonal entry

          DO i = map%h_diag_end_fixed + 1, map%n
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              j = prob%H%col( l )
              IF ( j > prob%n ) THEN
                prob%f = prob%f + prob%X( j ) * prob%H%val( l ) * prob%X( i ) 
              END IF
            END DO
          END DO
        END IF
      END IF

!  Process the bounds on c

      IF ( yes_c_bounds ) THEN
        DO j = prob%n + 1, map%n
          x = prob%X( j )
          IF ( x /= zero ) THEN
            DO l = map%ptr_a_fixed( j ), map%ptr_a_fixed( j + 1 ) - 1
              i = prob%A%col( l )   !  NB: this is the row number
              c = prob%A%val( l ) * x
              prob%C_l( i ) = prob%C_l( i ) - c
              prob%C_u( i ) = prob%C_u( i ) - c
            END DO
          END IF
        END DO
      END IF

      RETURN  

!  End of QPP_remove_fixed

      END SUBROUTINE QPP_remove_fixed

!-*-*-   Q P P _ a d d _ f i x e d   S U B R O U T I N E   -*-*-*

      SUBROUTINE QPP_add_fixed( map, prob, f, g, c, c_bounds )

!      ....................................................
!      .                                                  .
!      .  Transform f, g and the bounds on the            .
!      .  constraints when reintroducing fixed variables  .
!      .                                                  .
!      ....................................................

!  Arguments:
!  =========
!
!   map      see Subroutine QPP_initialize
!   prob%    see Subroutine QPP_reorder
!   f        LOGICAL. If true, adjust f to account for fixed variables via
!              f -> f - <g_x,x_x> - 1/2 <x_x,H_{xx} x_x> 
!            where x = (x_r,x_x ) (etc)
!   g        LOGICAL. If true, adjust g to account for fixed variables via
!              g_r -> g_r - H_{xr} x_x
!   c        LOGICAL. If true, adjust c to account for fixed variables via 
!              c <- c + A_x x_x
!   c_bounds LOGICAL. If true, adjust c_l and c_u to account for fixed 
!            variables via 
!              c_l <- c_l + A_x x_x and c_u <- c_u + A_x x_x

!  Dummy arguments

      TYPE ( QPP_map_type ), INTENT( IN ) :: map
      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      LOGICAL, OPTIONAL, INTENT( IN ) :: f, g, c, c_bounds

!  Local variables

      INTEGER :: i, j, l
      REAL ( KIND = wp ) :: x, c_val
      LOGICAL :: yes_f, yes_g, yes_c, yes_c_bounds

      IF ( map%n_reordered >= map%n ) RETURN

      IF ( PRESENT ( f ) ) THEN
        yes_f = f
      ELSE
        yes_f = .FALSE.
      END IF

      IF ( PRESENT ( g ) ) THEN
        yes_g = g
      ELSE
        yes_g = .FALSE.
      END IF

      IF ( PRESENT ( c ) ) THEN
        yes_c = c
      ELSE
        yes_c = .FALSE.
      END IF

      IF ( PRESENT ( c_bounds ) ) THEN
        yes_c_bounds = c_bounds
      ELSE
        yes_c_bounds = .FALSE.
      END IF

!  Process f and g together

      IF ( yes_g .AND. yes_f ) THEN

!  Least-distance Hessian
        
        IF ( prob%Hessian_kind == 1 ) THEN
          DO i = map%n_reordered + 1, map%n
            prob%f = prob%f - half * ( prob%X( i ) - prob%X0( i ) ) ** 2 
          END DO

!  Weighted-least-distance Hessian

        ELSE IF ( prob%Hessian_kind > 1 ) THEN
          DO i = map%n_reordered + 1, map%n
            prob%f = prob%f - half *                                           &
              ( prob%WEIGHT( i ) * ( prob%X( i ) - prob%X0( i ) ) ) ** 2 
          END DO

!  General Hessian

        ELSE IF ( prob%Hessian_kind /= 0 ) THEN

!  rows with a diagonal entry

          DO i = map%n_reordered + 1, map%h_diag_end_fixed
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 2
              j = prob%H%col( l )
              IF ( j <= map%n_reordered ) THEN
                prob%G( j ) = prob%G( j ) - prob%H%val( l ) * prob%X( i )
              ELSE
                prob%f = prob%f - prob%X( j ) * prob%H%val( l ) * prob%X( i ) 
              END IF
            END DO
            l = prob%H%ptr( i + 1 ) - 1
            prob%f = prob%f - half * prob%X( i ) * prob%H%val( l ) * prob%X( i ) 
          END DO

!  rows without a diagonal entry

          DO i = map%h_diag_end_fixed + 1, map%n
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              j = prob%H%col( l )
              IF ( j <= map%n_reordered ) THEN
                prob%G( j ) = prob%G( j ) - prob%H%val( l ) * prob%X( i )
              ELSE
                prob%f = prob%f - prob%X( j ) * prob%H%val( l ) * prob%X( i ) 
              END IF
            END DO
          END DO
        END IF

!  Process g separately

      ELSE IF ( yes_g ) THEN
        
!  rows with a diagonal entry

        DO i = map%n_reordered + 1, map%h_diag_end_fixed
          DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 2
            j = prob%H%col( l )
            IF ( j <= map%n_reordered ) THEN
              prob%G( j ) = prob%G( j ) - prob%H%val( l ) * prob%X( i )
            END IF
          END DO
        END DO

!  rows without a diagonal entry

        DO i = map%h_diag_end_fixed + 1, map%n
          DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
            j = prob%H%col( l )
            IF ( j <= map%n_reordered ) THEN
              prob%G( j ) = prob%G( j ) - prob%H%val( l ) * prob%X( i )
            END IF
          END DO
        END DO

!  Process f separately

      ELSE IF ( yes_f ) THEN

!  Least-distance Hessian
        
        IF ( prob%Hessian_kind == 1 ) THEN
          DO i = map%n_reordered + 1, map%n
            prob%f = prob%f - half * ( prob%X( i ) - prob%X0( i ) ) ** 2 
          END DO

!  Weighted-least-distance Hessian

        ELSE IF ( prob%Hessian_kind > 1 ) THEN
          DO i = map%n_reordered + 1, map%n
            prob%f = prob%f - half *                                           &
              ( prob%WEIGHT( i ) * ( prob%X( i ) - prob%X0( i ) ) ) ** 2 
          END DO

!  General Hessian

        ELSE IF ( prob%Hessian_kind /= 0 ) THEN

!  rows with a diagonal entry

          DO i = map%n_reordered + 1, map%h_diag_end_fixed
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 2
              j = prob%H%col( l )
              IF ( j > map%n_reordered ) THEN
                prob%f = prob%f - prob%X( j ) * prob%H%val( l ) * prob%X( i ) 
              END IF
            END DO
            l = prob%H%ptr( i + 1 ) - 1
            prob%f = prob%f - half * prob%X( i ) * prob%H%val( l ) * prob%X( i ) 
          END DO

!  rows without a diagonal entry

          DO i = map%h_diag_end_fixed + 1, map%n
            DO l = prob%H%ptr( i ), prob%H%ptr( i + 1 ) - 1
              j = prob%H%col( l )
              IF ( j > map%n_reordered ) THEN
                prob%f = prob%f - prob%X( j ) * prob%H%val( l ) * prob%X( i ) 
              END IF
            END DO
          END DO
        END IF
      END IF

      IF ( yes_f ) THEN
        IF ( prob%gradient_kind == 1 ) THEN
          prob%f = prob%f - SUM( prob%X( map%n_reordered + 1 : map%n ) )
        ELSE IF ( prob%gradient_kind /= 0 ) THEN
          prob%f = prob%f - DOT_PRODUCT( prob%G( map%n_reordered + 1 : map%n ),&
                                         prob%X( map%n_reordered + 1 : map%n ) )
        END IF
      END IF

      IF ( yes_c .AND. yes_c_bounds ) THEN

!  Process c and its bounds

        DO j = map%n_reordered + 1, map%n
          x = prob%X( j )
          IF ( x /= zero ) THEN
            DO l = map%ptr_a_fixed( j ), map%ptr_a_fixed( j + 1 ) - 1
              i = prob%A%col( l )   !  NB: this is the row number
              c_val = prob%A%val( l ) * x
              prob%C( i ) = prob%C( i ) + c_val
              prob%C_l( i ) = prob%C_l( i ) + c_val
              prob%C_u( i ) = prob%C_u( i ) + c_val
            END DO
          END IF
        END DO

!  Process the bounds on c

      ELSE IF ( yes_c_bounds ) THEN
        DO j = map%n_reordered + 1, map%n
          x = prob%X( j )
          IF ( x /= zero ) THEN
            DO l = map%ptr_a_fixed( j ), map%ptr_a_fixed( j + 1 ) - 1
              i = prob%A%col( l )   !  NB: this is the row number
              c_val = prob%A%val( l ) * x
              prob%C_l( i ) = prob%C_l( i ) + c_val
              prob%C_u( i ) = prob%C_u( i ) + c_val
            END DO
          END IF
        END DO

!  Process c

      ELSE IF ( yes_c ) THEN
        DO j = map%n_reordered + 1, map%n
          x = prob%X( j )
          IF ( x /= zero ) THEN
            DO l = map%ptr_a_fixed( j ), map%ptr_a_fixed( j + 1 ) - 1
              i = prob%A%col( l )   !  NB: this is the row number
              prob%C( i ) = prob%C( i ) + prob%A%val( l ) * x
            END DO
          END IF
        END DO
      END IF

      RETURN  

!  End of QPP_add_fixed

      END SUBROUTINE QPP_add_fixed

!-*-*-*-*-*-*-*-*-*-*-   Q P P _ A x  S U B R O U T I N E  -*-*-*-*-*-*-*-*-*-

      SUBROUTINE QPP_AX( map, prob_x, prob_A_val, prob_A_col, prob_A_ptr, m, Ax )

!      ..............................................
!      .                                            .
!      .  Perform the operation Ax := Ax + A * x    .
!      .                                            .
!      ..............................................

!  Arguments:
!  =========
!
!   map     see Subroutine QPP_initialize
!   prob_   see Subroutine QPP_reorder
!   m       row dimension of A
!   Ax      the result of adding A * x to Ax
!

!  Dummy arguments

      TYPE ( QPP_map_type ), INTENT( IN ) :: map
      INTEGER, INTENT( IN ), DIMENSION( map%m + 1 ) :: prob_A_ptr
      INTEGER, INTENT( IN ), DIMENSION( map%a_ne ) ::  prob_A_col
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( map%n ) :: prob_x
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( map%a_ne ) :: prob_A_val
      INTEGER, INTENT( IN ) :: m
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: Ax

!  Local variables

      INTEGER :: i, l

      DO i = 1, m
        DO l = prob_A_ptr( i ), prob_A_ptr( i + 1 ) - 1
          Ax( i ) = Ax( i ) + prob_A_val( l ) * prob_x( prob_A_col( l ) )
        END DO
      END DO

      RETURN

!  End of subroutine QPP_Ax

      END SUBROUTINE QPP_Ax

!  End of module QPP

   END MODULE GALAHAD_QPP_double




