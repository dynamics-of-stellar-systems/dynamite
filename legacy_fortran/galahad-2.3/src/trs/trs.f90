! THIS VERSION: GALAHAD 2.2 - 12/06/2008 AT 14:00 GMT.

!-*-*-*-*-*-*-*-  G A L A H A D _ T R S  double  M O D U L E  *-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released GALAHAD Version 2.2. June 3rd, 2008

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_TRS_double

!      -----------------------------------------------
!      |                                             |
!      | Solve the trust-region subproblem           |
!      |                                             |
!      |    minimize     1/2 <x, H x> + <c, x> + f   |
!      |    subject to   ||x||_2 <= radius           |
!      |    or           ||x||_2  = radius           |
!      |                                             |
!      | using a sparse matrix factorization         |
!      |                                             |
!      -----------------------------------------------

      USE GALAHAD_SYMBOLS
      USE GALAHAD_SPACE_double
      USE GALAHAD_RAND_double
      USE GALAHAD_SPECFILE_double
      USE GALAHAD_NORMS_double, ONLY: TWO_NORM
      USE GALAHAD_SILS_double

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: TRS_initialize, TRS_read_specfile, TRS_solve, TRS_terminate,    &
                SMT_type, SMT_put, SMT_get

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

      REAL ( KIND = wp ), PARAMETER :: theta = 0.01_wp
      REAL ( KIND = wp ), PARAMETER :: roots_tol = ten ** ( - 12 )
      REAL ( KIND = wp ), PARAMETER :: eps_get_u = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: eps_easy = ten ** ( - 12 )
      REAL ( KIND = wp ), PARAMETER :: eps_hard = ten ** ( - 12 )
      REAL ( KIND = wp ), PARAMETER :: macheps = EPSILON( one ) 

!--------------------------
!  Derived type definitions
!--------------------------

      TYPE, PUBLIC :: TRS_control_type
        INTEGER :: error, out, print_level, new_h
        REAL ( KIND = wp ) :: initial_multiplier
        LOGICAL :: boundary, equality_problem
        LOGICAL :: space_critical, deallocate_error_fatal
        TYPE ( SILS_control ) :: cholesky
        CHARACTER ( LEN = 30 ) :: prefix
      END TYPE

      TYPE, PUBLIC :: TRS_time_type
        REAL :: total, assemble, analyse, factorize, solve
      END TYPE

      TYPE, PUBLIC :: TRS_inform_type
        INTEGER :: status, alloc_status, factorizations, part_solve
        REAL ( KIND = wp ) :: obj, multiplier
        CHARACTER ( LEN = 80 ) :: bad_alloc
        TYPE ( TRS_time_type ) :: time
        TYPE ( SILS_AINFO ) :: analyse
        TYPE ( SILS_FINFO ) :: factorize
        TYPE ( SILS_SINFO ) :: solve
      END TYPE

      TYPE, PUBLIC :: TRS_data_type
        PRIVATE
        INTEGER :: h_ne
        TYPE ( RAND_seed ) :: seed
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: U, Y, Z
        TYPE ( SMT_type ) :: H_lambda
        TYPE ( SILS_factors ) :: factors
        TYPE ( TRS_control_type ) :: control
      END TYPE

    CONTAINS

!-*-*-*-*-*-*-  T R S _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE TRS_initialize( data, control )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!  .  Set initial values for the TRS control parameters  .
!   
!  Argument:
!  =========
!
!   data     private internal data
!   control  a structure containing control information. Components are -
!            out - output unit
!            error - error output unit
!            print_level - print level > 0 for output
!            initial_multiplier - initial estimate of the Lagrange multipler
!            new_h. How much of H has changed since the previous factorization.
!                   Possible values are
!                   0  unchanged
!                   1  values but not indices have changed
!                   2  values and indices have changed
!            space_critical - .TRUE. if allocated arrays are required to be 
!                            exactly the right size 
!            deallocate_error_fatal - .TRUE. if a deallocation error should
!                                    terminate execution
!            boundary - .TRUE. if the solution is thought to lie on the boundary
!            equality_problem - .TRUE. if the solution is REQUIRED to lie on the
!                             boundary (i.e., the constraint is an equality )
!            prefix - output lines will be prefixed by 
!                        prefix(2:LEN(TRIM(%prefix))-1)
!               where prefix contains the required string enclosed in 
!               quotes, e.g. "string" or 'string'
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!-----------------------------------------------
!   D u m m y   A r g u m e n t
!-----------------------------------------------

      TYPE ( TRS_DATA_TYPE ), INTENT( OUT ) :: data
      TYPE ( TRS_CONTROL_TYPE ), INTENT( OUT ) :: control        

!  Initalize random number seed

      CALL RAND_initialize( data%seed )

!  Set initial values for factorization controls and data

      CALL SILS_initialize( FACTORS = data%factors, CONTROL = control%cholesky )

!  Set initial control parameter values

      control%error = 6
      control%out = 6
      control%print_level = 0
      control%new_h = 2
      control%initial_multiplier = zero
      control%boundary = .FALSE.
      control%equality_problem = .FALSE.
      control%space_critical = .FALSE.
      control%deallocate_error_fatal  = .FALSE.
      control%prefix = '""     '

      RETURN

!  End of subroutine TRS_initialize

      END SUBROUTINE TRS_initialize

!-*-*-*-*-   T R S _ R E A D _ S P E C F I L E  S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE TRS_read_specfile( control, device, alt_specname )

!  Reads the content of a specification file, and performs the assignment of 
!  values associated with given keywords to the corresponding control parameters

!  The defauly values as given by TRS_initialize could (roughly) 
!  have been set as:

!  BEGIN TRS SPECIFICATIONS (DEFAULT)
!   error-printout-device                           6
!   printout-device                                 6
!   print-level                                     0
!   has-h-changed                                   2
!   intitial-multiplier                             0.0
!   solution-is-likely-on-boundary                  F
!   equality-problem                                F
!   space-critical                                  F
!   deallocate-error-fatal                          F
!   output-line-prefix                              ""
!  END TRS SPECIFICATIONS

!  Dummy arguments

      TYPE ( TRS_control_type ), INTENT( INOUT ) :: control        
      INTEGER, INTENT( IN ) :: device
      CHARACTER( LEN = 16 ), OPTIONAL :: alt_specname

!  Programming: Nick Gould and Ph. Toint, January 2002.

!  Local variables

      INTEGER, PARAMETER :: lspec = 16
      CHARACTER( LEN = 16 ), PARAMETER :: specname = 'TRS            '
      TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec

!  Define the keywords

     spec%keyword = ''

!  Integer key-words

      spec(  1 )%keyword = 'error-printout-device'
      spec(  2 )%keyword = 'printout-device'
      spec(  3 )%keyword = 'print-level' 
      spec(  4 )%keyword = 'has-h-changed'

!  Real key-words

      spec(  8 )%keyword = 'intitial-multiplier'

!  Logical key-words

      spec( 11 )%keyword = 'new-structure'
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
      CALL SPECFILE_assign_value( spec( 4 ), control%new_h,                    &
                                  control%error )

!  Set real values

      CALL SPECFILE_assign_value( spec( 8 ), control%initial_multiplier,       &
                                  control%error )

!  Set logical values

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

      END SUBROUTINE TRS_read_specfile

!-*-*-*-*-*-*-*-*-*-*  T R S _ S O L V E   S U B R O U T I N E  -*-*-*-*-*-*-*

      SUBROUTINE TRS_solve( n, radius, f, C, H, X, data, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!  =========
!
!   n - the number of unknowns
!
!   radius - the trust-region radius
!
!   f - the value of constant term for the quadratic function
!
!   C - a vector of values for the linear term c
!
!   H -  a structure of type SMT_type used to hold the LOWER TRIANGULAR part 
!    of the symmetric matrix H. Four storage formats are permitted:
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
!   X - the required solution vector x
!
!   data - private internal data
!
!   control - a structure containing control information. See TRS_initialize
!
!   inform - a structure containing information. Components are
!      status - the output status. Posiible values are
!        0 the solution has been found
!       -1 an array allocation has failed
!       -2 an array deallocation has failed
!       -3 n and/or radius is not positive
!      factorizations - the number of factorizations performed
!      obj - the value of the quadratic function
!      multiplier - the Lagrange multiplier corresponding to the 
!                  trust-region constraint
!      time%total - the total time spent in the package
!      time%assemble - the time spent building H+lambda*I
!      time%analyse - the time spent reordering H+lambda*I prior to
!        factorization
!      time%factorize - the time spent factorizing H+lambda*I
!      time%solve - the time spent solving linear systems inolving H+lambda*I
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ) :: radius
      REAL ( KIND = wp ), INTENT( IN ) :: f
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: C
      TYPE ( SMT_type ), INTENT( IN ) :: H
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X
      TYPE ( TRS_data_type ), INTENT( INOUT ) :: data
      TYPE ( TRS_control_type ), INTENT( IN ) :: control        
      TYPE ( TRS_inform_type ), INTENT( OUT ) :: inform

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, j, l, out, print_level, invit
      REAL ( KIND = wp ) :: lambda, lambda_l, lambda_u, delta_lambda
      REAL ( KIND = wp ) :: alpha, H_f, H_f2, H_inf, utx, distx, val, rayleigh
      REAL ( KIND = wp ) :: x_norm, c_norm, c_norm_over_radius, u_norm, w_norm2
      LOGICAL :: printi, printt, printd, psdef, try_zero, get_initial_u
      CHARACTER ( LEN = 1 ) :: region
      CHARACTER ( LEN = 80 ) :: array_name
      REAL :: time_now, time_start, time

!  prefix for all output 

      CHARACTER ( LEN = LEN( TRIM( control%prefix ) ) - 2 ) :: prefix
      prefix = control%prefix( 2 : LEN( TRIM( control%prefix ) ) - 1 )

      CALL CPU_TIME( time_start )

!  set initial values

      inform%alloc_status = 0 ; inform%bad_alloc = ''
      inform%factorizations = 0
      inform%time%assemble = 0.0 ; inform%time%analyse  = 0.0
      inform%time%factorize = 0.0 ; inform%time%solve = 0.0
      inform%time%total = 0.0

      X = zero ; inform%obj = f ; inform%multiplier = zero
      data%control = control

!  Check for obvious errors

      IF ( n <= 0 ) THEN
        IF ( control%error > 0 .AND. control%print_level > 0 )                  &
          WRITE( control%error,                                                 &
           "( A, ' n = ', I0, ' is too small ' )" ) prefix, n
        inform%status = GALAHAD_error_restrictions
        GO TO 900
      END IF

      IF ( radius <= zero ) THEN
        IF ( control%error > 0 .AND. control%print_level > 0 )                  &
          WRITE( control%error,                                                 &
          "( A, ' The radius ', ES12.4 , ' is not positive ' )" ) prefix, radius
        inform%status = GALAHAD_error_restrictions
        GO TO 900
      END IF

!  record desired output level

      out = control%out
      print_level = control%print_level
      printi = out > 0 .AND. print_level > 0
      printt = out > 0 .AND. print_level > 1
      printd = out > 0 .AND. print_level > 10

!  Choose initial values for the control parameters for the factorization

      IF ( data%control%print_level > 4 ) THEN
        data%control%cholesky%lp = out ; data%control%cholesky%mp = out
        data%control%cholesky%wp = out ; data%control%cholesky%sp = out
        data%control%cholesky%ldiag = 1
      ELSE
        data%control%cholesky%lp = 0 ; data%control%cholesky%mp = 0
        data%control%cholesky%wp = 0 ; data%control%cholesky%sp = 0
        data%control%cholesky%ldiag = 0
      END IF
      data%control%cholesky%pivoting = 2

!  ---------------------------------------------------------------
!  Set up data structure for the matrix H(lambda) = H + lambda * I
!  ---------------------------------------------------------------

!  compute the space required to hold the matrix H

      IF ( data%control%new_h >= 2 ) THEN
        SELECT CASE ( SMT_get( H%type ) )
        CASE ( 'DIAGONAL' ) 
          data%h_ne = n
        CASE ( 'DENSE' ) 
          data%h_ne = ( n * ( n - 1 ) ) / 2
        CASE ( 'SPARSE_BY_ROWS' )
          data%h_ne = H%ptr( n + 1 ) - 1
        CASE ( 'COORDINATE' )
          data%h_ne = H%ne
        END SELECT

!  add space required for lambda * I

        data%H_lambda%n = n
        data%H_lambda%ne = data%h_ne + n

!  allocate space to hold the matrix in co-ordinate form

        array_name = 'spls: H_lambda%row'
        CALL SPACE_resize_array( data%H_lambda%ne, data%H_lambda%row,           &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 910

        array_name = 'spls: H_lambda%col'
        CALL SPACE_resize_array( data%H_lambda%ne, data%H_lambda%col,           &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 910

        array_name = 'spls: H_lambda%val'
        CALL SPACE_resize_array( data%H_lambda%ne, data%H_lambda%val,           &
            inform%status, inform%alloc_status, array_name = array_name,        &
            deallocate_error_fatal = control%deallocate_error_fatal,            &
            exact_size = control%space_critical,                                &
            bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) GO TO 910

        CALL SMT_put( data%H_lambda%type, 'COORDINATE', inform%alloc_status )
        IF ( inform%alloc_status /= 0 ) THEN
          inform%status = GALAHAD_error_allocate
          GO TO 910 
        END IF
      END IF

!  fit the data from H into the coordinate storage scheme provided

      IF ( data%control%new_h >= 2 ) THEN
        SELECT CASE ( SMT_get( H%type ) )
        CASE ( 'DIAGONAL' ) 
          DO i = 1, n
            data%H_lambda%row( i ) = i ; data%H_lambda%col( i ) = i
          END DO
        CASE ( 'DENSE' ) 
          l = 0
          DO i = 1, n
            DO j = 1, i
              l = l + 1
              data%H_lambda%row( l ) = i ; data%H_lambda%col( l ) = j
            END DO
          END DO
        CASE ( 'SPARSE_BY_ROWS' )
          DO i = 1, n
            DO l = H%ptr( i ), H%ptr( i + 1 ) - 1
              data%H_lambda%row( l ) = i 
              data%H_lambda%col( l ) = H%col( l )
            END DO
          END DO
        CASE ( 'COORDINATE' )
          data%H_lambda%row( : data%h_ne ) = H%row( : data%h_ne )
          data%H_lambda%col( : data%h_ne ) = H%col( : data%h_ne )
        END SELECT

!  append the data for lambda * I

        DO i = 1, n
          data%H_lambda%row( data%h_ne + i ) = i
          data%H_lambda%col( data%h_ne + i ) = i
        END DO
      END IF

!  add the numerical values from H

      IF ( data%control%new_h >= 1 )                                           &
        data%H_lambda%val( : data%h_ne ) = H%val( : data%h_ne )

      CALL CPU_TIME( time_now )
      inform%time%assemble = time_now - time_start
      IF ( printt ) WRITE( out, "( /, A,  ' time( assembly ) = ', F10.2 )" )   &
        prefix, inform%time%assemble

!  Perform an analysis of the spasity pattern to identify a good
!  ordering for sparse factorization

      IF ( data%control%new_h >= 2 ) THEN
        CALL CPU_TIME( time )
        CALL SILS_analyse( data%H_lambda, data%factors, data%control%cholesky, &
                           inform%analyse )
        CALL CPU_TIME( time_now )
        inform%time%analyse  = time_now - time
        IF ( printt ) WRITE( out, 2000 ) prefix, inform%time%analyse
           
!  Test that the analysis succeeded

        IF ( inform%analyse%flag < 0 ) THEN
          IF ( printi ) WRITE( out, "( /, A, ' error return from ',            &
         &  'SILS_analyse: status = ', I0 )" ) prefix, inform%analyse%flag
          inform%status = GALAHAD_error_analysis ;  GO TO 900 ; END IF

      END IF

!  =====================
!  Array (re)allocations
!  =====================

!  Allocate Y ands Z

      array_name = 'trs: Y'
      CALL SPACE_resize_array( n, data%Y,                                      &
          inform%status, inform%alloc_status, array_name = array_name,         &
          deallocate_error_fatal = control%deallocate_error_fatal,             &
          exact_size = control%space_critical,                                 &
          bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 910

      array_name = 'trs: Z'
      CALL SPACE_resize_array( n, data%Z,                                      &
          inform%status, inform%alloc_status, array_name = array_name,         &
          deallocate_error_fatal = control%deallocate_error_fatal,             &
          exact_size = control%space_critical,                                 &
          bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= 0 ) GO TO 910

!  Compute the sums of the absolute values of off-diagonal terms of H (in Y),
!  its diagonal terms (in Z) and the square of its Frobenius norm

      data%Y( : n ) = zero ; data%Z( : n ) = zero ; H_f2 = zero
      DO l = 1, data%h_ne
        i = data%H_lambda%row( l ) ; j = data%H_lambda%col( l )
        val = data%H_lambda%val( l )
        IF ( i == j ) THEN
          data%Z( i ) = val
          H_f2 = H_f2 + val ** 2
        ELSE
          data%Y( i ) = data%Y( i ) + ABS( val )
          data%Y( j ) = data%Y( j ) + ABS( val )
          H_f2 = H_f2 + two * val ** 2
        END IF
      END DO

!  compute the Frobenius and infinity norms of H

      H_f = SQRT( H_f2 )
      H_inf = MAXVAL( ABS( data%Z( : n ) ) + data%Y( : n ) )

!  The real line is partitioned into disjoint sets
!     N = { lambda: lambda <= max(0, -lambda_1(H))}
!     L = { lambda: max(0, -lambda_1(H)) < lambda <= lambda_optimal } and
!     G = { lambda: lambda > lambda_optimal }.
!  The aim is to find a lambda in L, as generally then Newton's method 
!  will converge both globally and ultimately quadratically. We also let
!     F = L union G
!
!  We will construct values lambda_l and lambda_u for which
!  lambda_l <= lambda_optimal <= lambda_u, and ensure that all iterates
!  satisfy lambda_l <= lambda <= lambda_u. 

      c_norm = TWO_NORM( C )
      c_norm_over_radius = c_norm / radius
      lambda_l = MAX( zero, - MINVAL( data%Z( : n ) ), c_norm_over_radius -    &
        MIN( MAXVAL( data%Z( : n ) + data%Y( : n ) ), H_f, H_inf ) )
      lambda_u = MAX( zero, c_norm_over_radius +                               &
        MIN( MAXVAL( - data%Z( : n ) + data%Y( : n ) ), H_f, H_inf ) )
      IF ( lambda_l > lambda_u ) THEN
        WRITE( out, "( ' lambda_l = ', ES12.4, ' > lambda_u = ', ES12.4 )" )   &
          lambda_l, lambda_u
        STOP
      END IF

!  =========================
!  Phase 1: find lambda in L
!  =========================

!  assign the initial lambda

      lambda = MAX( data%control%initial_multiplier, lambda_l )
      try_zero = lambda > zero .AND. lambda_l == zero
      region = ' '
      get_initial_u = .TRUE.

      IF ( printi ) THEN
        WRITE( out, "( /, A, ' Phase 1' )" ) prefix
        WRITE( out, 2020 ) prefix
      END IF

      DO
        IF ( printi ) WRITE( out, "( A, A2, 3ES23.15 )" )                       &
          prefix, region, lambda_l, lambda, lambda_u

!  Add lambda * I to H to form H(lambda)

        data%H_lambda%val( data%h_ne + 1 : data%H_lambda%ne ) = lambda

!WRITE( out, "( ' H (co-ordinate) = ' )" )
!WRITE( out, "( ( 2( 2I8, ES12.4 ) ) )" )                     &
!    ( data%H_lambda%row( i ), data%H_lambda%col( i ),        &
!      data%H_lambda%val( i ), i = 1, data%H_lambda%ne )


!  Attempt a Cholesky factorization of H(lambda)

        CALL CPU_TIME( time )
        CALL SILS_factorize( data%H_lambda, data%factors,                       &
                             data%control%cholesky, inform%factorize )
        CALL CPU_TIME( time_now )
        time_now = time_now - time
        inform%time%factorize = inform%time%factorize + time_now
        IF ( printt ) WRITE( out, 2010 ) prefix, time_now
        inform%factorizations = inform%factorizations + 1

!  Test that the factorization succeeded

        IF ( inform%factorize%flag == 0 ) THEN
          psdef = .TRUE.
        ELSE IF ( inform%factorize%flag > 0 .OR.                               &
                  inform%factorize%flag == - 5 .OR.                            &
                  inform%factorize%flag == - 6 ) THEN
          psdef = .FALSE.
        ELSE
          GO TO 920
        END IF
      
!  If H(lambda) is positive definite, solve  H(lambda) x = - c

        IF ( psdef ) THEN 
          CALL CPU_TIME( time )
!write(6,"('C', 3ES12.4)" ) C(:n)
          X = - C
          CALL SILS_solve( data%H_lambda, data%factors, X,                      &
                           data%control%cholesky, inform%solve )
          CALL CPU_TIME( time_now )
          inform%time%solve = inform%time%solve + time_now - time
          x_norm = TWO_NORM( X ) 
! write(6,"('X', 3ES12.4)" ) X(:n)
          IF ( printd ) THEN
            WRITE( out, "( A, 8X, 'lambda', 13X, 'x_norm', 15X, 'radius' )" )   &
              prefix
            WRITE( out, "( A, 3ES20.12 )" ) prefix, lambda, x_norm, radius
          END IF

!  The current estimate gives a good approximation to the required root
!  --------------------------------------------------------------------

          IF ( ABS( x_norm - radius ) <= eps_easy * radius ) THEN
            inform%status = 0
            GO TO 890
          END IF

          IF ( x_norm > radius ) THEN

!  The current lambda lies in L - jump to the Newton loop
!  ------------------------------------------------------
            
            region = 'L'
            EXIT
          ELSE

!  The solution lies in the interior of the trust-region 
!  -----------------------------------------------------

            IF ( lambda == zero ) THEN
              inform%obj = half * DOT_PRODUCT( C, X )
              inform%status = 0
              GO TO 900

!  The current lambda lies in G
!  ----------------------------

            ELSE
              region = 'G'
              IF ( try_zero ) THEN
                lambda = zero
                try_zero = .FALSE.
                CYCLE
              END IF

              lambda_u = lambda

!  Compute ||w|| where ||w||^2 = < x, H^{-1}(lambda) x >

              CALL CPU_TIME( time )
              data%Y( : n ) = X
              CALL SILS_part_solve( data%factors, data%control%cholesky, 'L',   &
                                    data%Y( : n ), inform%part_solve )
              data%Z( : n ) = data%Y( : n )
              CALL SILS_part_solve( data%factors, data%control%cholesky, 'D',   &
                                    data%Z( : n ), inform%part_solve )
              CALL CPU_TIME( time_now )
              inform%time%solve = inform%time%solve + time_now - time
              w_norm2 = DOT_PRODUCT( data%Y( : n ), data%Z( : n ) )

!  Compute the Newton correction, and project this into [lambda_l,lambda_u]

              lambda_l = MAX( lambda_l, lambda +                               &
                 ( ( x_norm - radius ) / radius ) * ( x_norm ** 2 / w_norm2 ) )

!  If it seems as if it may be necessary, build an estmate of the leftmost 
!  eigenvalue and its vector u using inverse iteration

              IF ( lambda_u - lambda_l < eps_get_u ) THEN

!  allocate space to hold u

                IF ( get_initial_u ) THEN
                  get_initial_u = .FALSE.
                  array_name = 'trs: U'
                  CALL SPACE_resize_array( n, data%U,                           &
                    inform%status, inform%alloc_status, array_name = array_name,&
                    deallocate_error_fatal = control%deallocate_error_fatal,    &
                    exact_size = control%space_critical,                        &
                    bad_alloc = inform%bad_alloc, out = control%error )
                  IF ( inform%status /= 0 ) GO TO 910

!  start the inverse iteration with a random vector orthogonal to c

                  DO i = 1, n
                    CALL RAND_random_real( data%seed, .FALSE., data%U( i ) )
                  END DO
                  IF ( c_norm > zero ) THEN
                    alpha = DOT_PRODUCT( C, data%U( : n ) ) / c_norm ** 2
                    data%U( : n ) = data%U( : n ) - alpha * C
                  END IF

!  normalize u

                  u_norm = TWO_NORM( data%U( : n ) )
                  data%U( : n ) = data%U( : n ) / u_norm
                END IF

!  now perform a few iterations of inverse iteration

                invit = 1
                DO i = 1, invit
                  CALL CPU_TIME( time )
                  CALL SILS_solve( data%H_lambda, data%factors, data%U( : n ),  &
                                   data%control%cholesky, inform%solve )
                  CALL CPU_TIME( time_now )
                  inform%time%solve = inform%time%solve + time_now - time
                  u_norm = TWO_NORM( data%U( : n ) )
                  data%U( : n ) = data%U( : n ) / u_norm
                  lambda_l = MAX( lambda_l,  lambda - one / u_norm )
                END DO    

!  Compute the Rayleigh quotient

                data%Y( : n ) = zero
                DO l = 1, data%h_ne
                  i = data%H_lambda%row( l ) ; j = data%H_lambda%col( l )
                  val = data%H_lambda%val( l )
                  IF ( i == j ) THEN
                    data%Y( i ) = data%Y( i ) + val * data%U( i )
                  ELSE
                    data%Y( i ) = data%Y( i ) + val * data%U( j )
                    data%Y( j ) = data%Y( j ) + val * data%U( i )
                  END IF
                END DO
                rayleigh = DOT_PRODUCT( data%U( : n ), data%Y( : n ) )

!  adjust lambda_l to account for the Rayleigh quotient

                lambda_l = MAX( lambda_l,  - rayleigh )

!  compute the next lambda - bias this towards lambda_l

!               lambda = lambda_l + theta * ( lambda_u - lambda_l )
                lambda = lambda_l + ( lambda_u - lambda_l ) ** 1.5
              ELSE
                lambda = MAX( SQRT( lambda_l * lambda_u ),                    &
                              lambda_l + theta * ( lambda_u - lambda_l ) )
              END IF
            END IF
          END IF
        ELSE

!  The current lambda lies in N
!  ----------------------------

          region = 'N'
          try_zero = .FALSE.
          lambda_l = lambda
          lambda = MAX( SQRT( lambda_l * lambda_u ),                           &
                     lambda_l + theta * ( lambda_u - lambda_l ) )
        END IF

        IF ( lambda_u - lambda_l < eps_hard ) THEN

!  The "hard" case
!  ---------------

      IF ( printi ) WRITE( out, "( A, ' Hard case occured' )" ) prefix

!  Compute the step alpha so that X + alpha U lies on the trust-region
!  boundary and gives the smaller value of q

          utx = DOT_PRODUCT( data%U, X ) / radius 
          distx = ( radius - x_norm ) * ( ( radius + x_norm ) / radius ) 
          alpha = sign( distx / ( abs( utx ) +                                 &
                        sqrt( utx ** 2 + distx / radius ) ), utx )

!  Record the optimal values

          X = X + alpha * data%U
          inform%obj = half * ( DOT_PRODUCT( C, X ) - lambda * radius ** 2 )
          inform%status = 0
          GO TO 900
        END IF
        IF ( printt ) WRITE( out, 2020 ) prefix
      END DO

!  ====================
!  Phase 2: lambda in L
!  ====================

      delta_lambda = lambda_u - lambda_l
      IF ( printi ) THEN
        WRITE( out, "( /, A, ' Phase 2' )" ) prefix
        WRITE( out, 2030 ) prefix
      END IF

!  We have found a lambda in L. It is now simply a matter of applying 
!  Newton's method starting from this lambda

      DO
        IF ( printi ) WRITE( out, "( A, 2X, 3ES23.15 )" )                      &
          prefix, lambda, delta_lambda, x_norm - radius 

!  Compute ||w|| where ||w||^2 = < x, H^{-1}(lambda) x >

        CALL CPU_TIME( time )
        data%Y( : n ) = X

        CALL SILS_part_solve( data%factors, data%control%cholesky, 'L',        &
                              data%Y( : n ), inform%part_solve )
        data%Z( : n ) = data%Y( : n )
        CALL SILS_part_solve( data%factors, data%control%cholesky, 'D',        &
                              data%Z( : n ), inform%part_solve )
        CALL CPU_TIME( time_now )
        inform%time%solve = inform%time%solve + time_now - time
        w_norm2 = DOT_PRODUCT( data%Y( : n ), data%Z( : n ) )

!  Compute the Newton correction

        delta_lambda = ( ( x_norm - radius ) / radius ) *                      &
                        ( x_norm ** 2 / w_norm2 )

!  Check that the Newton correction is significant

!       WRITE( out, "( ' lambda, delta ', 2ES22.14 )" ) lambda, delta_lambda

        IF ( ABS( delta_lambda ) < macheps * ABS( lambda ) ) THEN
          inform%status = 0
          GO TO 890
        END IF

!  Compute the new estimate of lambda

        lambda = lambda + delta_lambda

!  Add lambda * I to H to form H(lambda)

        data%H_lambda%val( data%h_ne + 1 : data%H_lambda%ne ) = lambda

!  Attempt a Cholesky factorization of H + lamda * I

        CALL CPU_TIME( time )
        CALL SILS_factorize( data%H_lambda, data%factors,                      &
                             data%control%cholesky, inform%factorize )
        CALL CPU_TIME( time_now )
        time_now = time_now - time
        inform%time%factorize = inform%time%factorize + time_now
        IF ( printt ) WRITE( out, 2010 ) prefix, time_now
        inform%factorizations = inform%factorizations + 1

!  Test that the factorization succeeded

        IF ( inform%factorize%flag == 0 ) THEN
          psdef = .TRUE.
        ELSE IF ( inform%factorize%flag > 0 .OR.                               &
                  inform%factorize%flag == - 5 .OR.                            &
                  inform%factorize%flag == - 6 ) THEN
           psdef = .FALSE.
           write(6,*) ' should not get indefinite H(lambda) in phase 2 '
           STOP
        ELSE
          GO TO 920
        END IF

!  If H(lambda) is positive definite, solve  H(lambda) x = - c

        CALL CPU_TIME( time )
        X = - C
        CALL SILS_solve( data%H_lambda, data%factors, X,                        &
                         data%control%cholesky, inform%solve )
        CALL CPU_TIME( time_now )
        inform%time%solve = inform%time%solve + time_now - time
        x_norm = TWO_NORM( X ) 

        IF ( printd ) WRITE( out, "( A, 3ES20.12 )" )                           &
          prefix, lambda, x_norm, radius

!  Test for convergence

        IF ( x_norm - radius <= eps_easy * radius ) THEN
          inform%status = 0
          GO TO 890
        END IF
        IF ( printt ) WRITE( out, 2030 ) prefix
      END DO

!  Too many iterations

      inform%status = GALAHAD_error_max_iterations

 890  CONTINUE
      inform%obj = half * ( DOT_PRODUCT( C, X ) - lambda * x_norm ** 2 )

!  Exit

 900  CONTINUE
      CALL CPU_TIME( time_now )
      inform%time%solve = time_now - time_start
      inform%multiplier = lambda
      RETURN

!  Allocation error

  910 CONTINUE
      IF ( printi ) WRITE( out, "( ' ', /, A, '   **  Error return ', I0,      &
     & ' from TRS ' )" ) prefix, inform%status 
      CALL CPU_TIME( time_now )
      inform%time%solve = time_now - time_start
      RETURN

!  Factorization failure

  920 CONTINUE
      IF ( printi ) WRITE( out, "( /, A, ' error return from ',                &
     &   'SILS_factorize: status = ', I0 )" ) prefix, inform%factorize%flag
      inform%status = GALAHAD_error_factorization 
      CALL CPU_TIME( time_now )
      inform%time%solve = time_now - time_start
      RETURN

! Non-executable statements

 2000 FORMAT( /, A, ' time( SILS_analyse ) = ', F10.2 )
 2010 FORMAT( /, A, ' time( SILS_factorize ) = ', F10.2 )
 2020 FORMAT( A, '          lambda_l                lambda ',                &
                 '               lambda_u' )
 2030 FORMAT( A, '           lambda                d_lambda',                &
                 '             ||x||-radius' )

!  End of subroutine TRS_solve

      END SUBROUTINE TRS_solve

!-*-*-*-*-*-  T R S _ T E R M I N A T E   S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE TRS_terminate( data, control, inform )

!  ...........................................
!  .                                         .
!  .  Deallocate arrays at end of TRS_solve  .
!  .                                         .
!  ...........................................

!  Arguments:
!  =========
!
!   data    private internal data
!   control see Subroutine TRS_initialize
!   inform    see Subroutine TRS_solve

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      TYPE ( TRS_data_type ), INTENT( INOUT ) :: data
      TYPE ( TRS_control_type ), INTENT( IN ) :: control        
      TYPE ( TRS_inform_type ), INTENT( INOUT ) :: inform

!-----------------------------------------------
!   L o c a l   V a r i a b l e
!-----------------------------------------------

      CHARACTER ( LEN = 80 ) :: array_name

!  Deallocate all internal arrays

      array_name = 'trs: Y'
      CALL SPACE_dealloc_array( data%Y,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

      array_name = 'trs: Z'
      CALL SPACE_dealloc_array( data%Y,                                        &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND. inform%status /= 0 ) RETURN

!  Deallocate all arrays allocated within SILS

      CALL SILS_finalize( data%factors, control%cholesky, inform%alloc_status )
      IF ( inform%alloc_status /= 0 ) THEN
        inform%status = GALAHAD_error_deallocate
        inform%bad_alloc = 'trs: factors'
      END IF

      RETURN

!  End of subroutine TRS_terminate

      END SUBROUTINE TRS_terminate

!-*-*-*-*-*-  End of G A L A H A D _ T R S  double  M O D U L E  *-*-*-*-*-*-

   END MODULE GALAHAD_TRS_double
