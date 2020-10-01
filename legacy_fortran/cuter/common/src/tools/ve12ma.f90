! ( Last modified on 21 Dec 2000 at 23:46:51 )

!-*-*-*-*-*-*-*-*-*-*-*-  VE12_main  *-*-*-*-*-*-*-*-*-*-*-*-*-*

!  Nick Gould
!  Copyright reserved
!  June 12th 1995

PROGRAM VE12_main

  !  Main program for HSL_VE12, an algorithm to solve quadratic programs
  !  using a primal-dual interior-point trust-region iteration

  USE CUTEr_interfaces
  USE READ_input
!S    USE HSL_VE12_single
!D    USE HSL_VE12_double

  IMPLICIT NONE

  !  --------------------------------------------------------------------
  !  Solve the quadratic program
  !
  !     minimize     1/2 x(T) H x + g(T) x
  !
  !     subject to     c_l <= A x <= c_u
  !                    x_l <=  x <= x_u
  !
  !  using VE12
  !
  !  --------------------------------------------------------------------

  !  Parameters

  INTEGER, PARAMETER :: out  = 6
  INTEGER, PARAMETER :: input = 55, inputd = 33
  REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
  REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
  REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp

  !  Scalars

  INTEGER :: n, m, ir, ic, nmax, mmax, la, lh, liw
  INTEGER :: i, j, l, neh, nea, factorization_integer, factorization_real
  INTEGER :: status, mfixed, mdegen, iter, nfacts, nfixed, ndegen, mequal
  INTEGER :: alloc_stat, A_ne, H_ne
  REAL :: CPU( 2 ), CALLS( 7 )
  REAL :: time, timeo, times, timet
  REAL ( KIND = wp ) :: objf, qfval, qfvalt, stopr
  LOGICAL :: fulsol

  CHARACTER ( LEN =  5 ) :: state, solv
  CHARACTER ( LEN = 10 ) :: pname

  !  Arrays

  TYPE ( VE12_data_type ) :: data
  TYPE ( VE12_control_type ) :: VE12_control        
  TYPE ( VE12_inform_type ) :: VE12_inform
  TYPE ( ZD02_problem_type ) :: prob

  !  Allocatable arrays

  CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: VNAME, CNAME
  REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X0, C
  LOGICAL, ALLOCATABLE, DIMENSION( : ) :: EQUATN, LINEAR
  INTEGER, ALLOCATABLE, DIMENSION( : ) :: IW

  !  Open the relevant file.

  OPEN( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD'  )
  REWIND input

  CALL CPU_TIME( time )

  !  Determine the number of variables and constraints

  CALL CDIMEN( input, nmax, mmax )
  nmax = MAX( nmax, 1 ) ; mmax = MAX( mmax, 1 ) ! for fortran 77 CUTEr tools

  !  Allocate suitable arrays

  ALLOCATE( prob%X( nmax ), prob%X_l( nmax ), prob%X_u( nmax ),            &
       prob%G( nmax ), VNAME( nmax ), STAT = alloc_stat )
  IF ( alloc_stat /= 0 ) THEN
     WRITE( out, 2150 ) 'X', alloc_stat ; STOP
  END IF

  ALLOCATE( prob%C_l( mmax ), prob%C_u( mmax ), prob%Y( mmax ), &
       CNAME( mmax ), EQUATN( mmax ), LINEAR( mmax ),                 &
       STAT = alloc_stat )
  IF ( alloc_stat /= 0 ) THEN
     WRITE( out, 2150 ) 'C', alloc_stat ; STOP
  END IF

  !  Set up the data structures necessary to hold the group partially
  !  separable function.

  CALL CSETUP( input, out, n, m, prob%X, prob%X_l, prob%X_u, nmax, EQUATN, &
       LINEAR, prob%Y, prob%C_l, prob%C_u, mmax, .FALSE.,          &
       .FALSE., .FALSE. )

  DEALLOCATE( LINEAR )

  !  Allocate derived types

  ALLOCATE( X0( n ), STAT = alloc_stat )
  IF ( alloc_stat /= 0 ) THEN
     WRITE( out, 2150 ) 'X0', alloc_stat
     STOP
  END IF

  ALLOCATE( prob%C( m ), C( m ), STAT = alloc_stat )
  IF ( alloc_stat /= 0 ) THEN
     WRITE( out, 2150 ) 'C', alloc_stat
     STOP
  END IF

  !  Determine the names of the problem, variables and constraints.

  CALL CNAMES( n, m, pname, VNAME, CNAME )
  WRITE( out, 2020 ) pname 
  WRITE( out, "( /, ' n   = ', i8, ' nmax = ', i8,                         &
       & ' m    = ', i8, ' mmax = ', i8 )" ) n, nmax, m, mmax

  !  Set up the initial estimate of the solution and
  !  right-hand-side of the Kuhn-Tucker system.

  !  Determine the constant terms for the problem functions.

  prob%X( : n ) = MIN( prob%X_u( : n ),                                    &
       MAX( prob%X_l( : n ), prob%X( : n ) ) )

  !  Set X0 to zero to determine the constant terms for the problem functions

  X0 = zero 

  !  Evaluate the constant terms of the objective (objf) and constraint 
  !  functions (C)

  CALL CFN( n, m, X0, objf, m, C( : m ) )
  DO i = 1, m 
     IF ( EQUATN( i ) ) THEN 
        prob%C_l( i ) = prob%C_l( i ) - C( i )
        prob%C_u( i ) = prob%C_l( i )
     ELSE
        prob%C_l( i ) = prob%C_l( i ) - C( i )
        prob%C_u( i ) = prob%C_u( i ) - C( i )
     END IF
  END DO

  !  Determine the number of nonzeros in the Jacobian

  CALL CDIMSJ( la ) ; la = MAX( la, 1 )

  !  Allocate arrays to hold the Jacobian

  ALLOCATE( prob%A_row( la ), prob%A_col( la ), prob%A_val( la ),          &
       STAT = alloc_stat )
  IF ( alloc_stat /= 0 ) THEN
     WRITE( out, 2150 ) 'A', alloc_stat ; STOP
  END IF

  !  Evaluate the linear terms of the constraint functions

  CALL CSGR( n, m, .FALSE., mmax, prob%Y, X0, nea, la, prob%A_val,         &
       prob%A_col, prob%A_row )
  WRITE( out, "( ' nea = ', i8, ' la   = ', i8 )", ADVANCE = 'no' ) &
       nea, la

  DEALLOCATE( X0 )

  !  Exclude zeros; set the linear term for the objective function

  A_ne = 0
  prob%G( : n ) = zero
  DO i = 1, nea
     IF ( prob%A_val( i ) /= zero ) THEN
        IF ( prob%A_row( i ) > 0 ) THEN
           A_ne = A_ne + 1
           prob%A_row( A_ne ) = prob%A_row( i ) 
           prob%A_col( A_ne ) = prob%A_col( i )
           prob%A_val( A_ne ) = prob%A_val( i )
        ELSE
           prob%G( prob%A_col( i ) ) = prob%A_val( i )
        END IF
     END IF
  END DO

  !  Determine the number of nonzeros in the Hessian

  CALL CDIMSH( lh ) ; lh = MAX( lh, 1 )

  !  Allocate arrays to hold the Hessian

  ALLOCATE( prob%H_row( lh ), prob%H_col( lh ), prob%H_val( lh ),          &
       STAT = alloc_stat )
  IF ( alloc_stat /= 0 ) THEN
     WRITE( out, 2150 ) 'H', alloc_stat
     STOP
  END IF

  !  Evaluate the Hessian of the Lagrangian function at the initial point.

  CALL CSH( n, m, prob%X, mmax, prob%Y, neh, lh, prob%H_val, prob%H_row,   &
       prob%H_col )
  WRITE( out, "( ' neh  = ', i8, ' lh   = ', i8 )" ) neh, lh

  !  Remove Hessian out of range

  H_ne = 0
  DO l = 1, neh    
     i = prob%H_row( l ) ; j = prob%H_col( l )
     IF ( i < 1 .OR. i > n .OR. j < 1 .OR. j > n ) CYCLE
     H_ne = H_ne + 1 ; prob%H_val( H_ne ) = prob%H_val( l )
     IF ( i >= j ) THEN
        prob%H_row( H_ne ) = i
        prob%H_col( H_ne ) = j
     ELSE
        prob%H_row( H_ne ) = j
        prob%H_col( H_ne ) = i
     END IF
  END DO

  !  Allocate and initialize dual variables.

  ALLOCATE( prob%Z( nmax ), STAT = alloc_stat )
  IF ( alloc_stat /= 0 ) THEN
     WRITE( out, 2150 ) 'Z', alloc_stat
     STOP
  END IF
  prob%Z( : n ) = one

  liw = MAX( m, n ) + 1
  ALLOCATE( prob%A_ptr( m + 1 ), prob%H_ptr( n + 1 ) )
  ALLOCATE( IW( liw ) )

  !  Transform A to row storage format

  IF ( A_ne /= 0 ) THEN
     CALL REORDER_by_rows( m, n, A_ne, prob%A_row, prob%A_col, A_ne,       &
          prob%A_val, prob%A_ptr, m + 1, IW, liw, i )
  ELSE
     prob%A_ptr = 0
  END IF

  !  Same for H

  IF ( H_ne /= 0 ) THEN
     CALL REORDER_by_rows( n, n, H_ne, prob%H_row, prob%H_col, H_ne,        &
          prob%H_val, prob%H_ptr, n + 1, IW, liw, i )
  ELSE
     prob%H_ptr = 0
  END IF

  !  Deallocate arrays holding matrix row indices

  DEALLOCATE( prob%A_row, prob%H_row )
  DEALLOCATE( IW )
  ALLOCATE( prob%A_row( 0 ), prob%H_row( 0 ) )

  prob%new_problem_structure = .TRUE.

  !  Store the problem dimensions

  prob%n = n
  prob%m = m
  prob%a_ne = - 1
  prob%h_ne = - 1
  prob%f = objf

  !  ------------------- problem set-up complete ----------------------

  CALL CPU_TIME( times )

  !  Open the spec file

  OPEN( inputd, FILE = 'VE12.SPC', FORM = 'FORMATTED', STATUS = 'OLD' )
  REWIND inputd

  !  Update control parameters if required.

  CALL VE12_initialize( data, VE12_control )

  CALL OVERIDE_control( VE12_control%infinity, inputd )
  CALL OVERIDE_control( VE12_control%stop_p , inputd )
  CALL OVERIDE_control( VE12_control%stop_d , inputd )
  CALL OVERIDE_control( VE12_control%stop_c , inputd )
  CALL OVERIDE_control( VE12_control%prfeas, inputd )
  CALL OVERIDE_control( VE12_control%dufeas, inputd )
  CALL OVERIDE_control( VE12_control%muzero, inputd )
  CALL OVERIDE_control( VE12_control%inner_stop_absolute, inputd )
  CALL OVERIDE_control( VE12_control%inner_stop_relative, inputd )
  CALL OVERIDE_control( VE12_control%inner_fraction_opt, inputd )
  CALL OVERIDE_control( VE12_control%pivot_tol, inputd )
  CALL OVERIDE_control( VE12_control%pivot_tol_for_dependencies, inputd )
  CALL OVERIDE_control( VE12_control%initial_radius, inputd )
  CALL OVERIDE_control( VE12_control%VE13_control%reduce_infeas, inputd )

  CALL OVERIDE_control( VE12_control%maxit, inputd )
  CALL OVERIDE_control( VE12_control%print_level, inputd )
  CALL OVERIDE_control( VE12_control%out, inputd )
  CALL OVERIDE_control( VE12_control%error, inputd )
  CALL OVERIDE_control( VE12_control%itrmax, inputd )
  CALL OVERIDE_control( VE12_control%cg_maxit, inputd )
  CALL OVERIDE_control( VE12_control%factor, inputd )
  CALL OVERIDE_control( VE12_control%max_col, inputd )
  CALL OVERIDE_control( VE12_control%precon, inputd )
  CALL OVERIDE_control( VE12_control%nsemib, inputd )
  CALL OVERIDE_control( VE12_control%indmin, inputd )
  CALL OVERIDE_control( VE12_control%valmin, inputd )
  CALL OVERIDE_control( VE12_control%VE13_control%infeas_max, inputd )

  CALL OVERIDE_control( VE12_control%feasol, inputd )
  CALL OVERIDE_control( fulsol, inputd )
  CALL OVERIDE_control( VE12_control%primal, inputd )
  CALL OVERIDE_control( VE12_control%center, inputd )

  VE12_control%restore_problem = 1

  CLOSE( inputd )

  !  Call the optimizer

  qfval = objf 

  CALL CPU_TIME( timeo )

  prob%m = m
  prob%n = n

  solv = ' VE12'
  CALL VE12_solve( prob, data, VE12_control, VE12_inform )
  qfval = VE12_inform%obj

  CALL CREPRT( CALLS, CPU )
  CALL CPU_TIME( timet )

  !  Deallocate arrays from the minimization

  status = VE12_inform%status ; iter = VE12_inform%iter
  nfacts = VE12_inform%nfacts ; stopr = VE12_control%stop_d
  factorization_integer = VE12_inform%factorization_integer 
  factorization_real = VE12_inform%factorization_real
  CALL VE12_terminate( data, VE12_control, VE12_inform )

  !  Print details of the solution obtained

  WRITE( out, 2010 ) status
  IF ( status >= 0 .AND. status <= 5 ) THEN 
     l = 4 ; IF ( fulsol ) l = n 

     !  Print details of the primal and dual variables

     WRITE( out, 2090 ) 
     DO j = 1, 2 
        IF ( j == 1 ) THEN 
           ir = 1 ; ic = MIN( l, n ) 
        ELSE 
           IF ( ic < n - l ) WRITE( out, 2000 ) 
           ir = MAX( ic + 1, n - ic + 1 ) ; ic = n 
        END IF
        DO i = ir, ic 
           state = ' FREE' 
           IF ( ABS( prob%X  ( i ) - prob%X_l( i ) ) < ten * stopr )          &
                state = 'LOWER'
           IF ( ABS( prob%X  ( i ) - prob%X_u( i ) ) < ten * stopr )          &
                state = 'UPPER'
           IF ( ABS( prob%X_l( i ) - prob%X_u( i ) ) <     1.0D-10 )          &
                state = 'FIXED'
           WRITE( out, 2050 ) VNAME( i ), state, prob%X( i ),                 &
                prob%X_l( i ), prob%X_u( i ), prob%Z( i )
        END DO
     END DO

     !  Compute the number of fixed and degenerate variables.

     nfixed = 0 ; ndegen = 0 
     DO i = 1, n 
        IF ( ABS( prob%X( i ) - prob%X_l( i ) ) < stopr ) THEN
           nfixed = nfixed + 1 
           IF ( ABS( prob%Z( i ) ) < ten * stopr ) ndegen = ndegen + 1 
        ELSE IF ( ABS( prob%X( i ) - prob%X_u( i ) ) < ten * stopr ) THEN
           nfixed = nfixed + 1 
           IF ( ABS( prob%Z( i ) ) < ten * stopr ) ndegen = ndegen + 1 
        END IF
     END DO

     !  Print details of the constraints.

     IF ( m > 0 ) THEN 

        WRITE( out, 2040 ) 
        l = 2  ; IF ( fulsol ) l = m 
        DO j = 1, 2 
           IF ( j == 1 ) THEN 
              ir = 1 ; ic = MIN( l, m ) 
           ELSE 
              IF ( ic < m - l ) WRITE( out, 2000 ) 
              ir = MAX( ic + 1, m - ic + 1 ) ; ic = m 
           END IF
           DO i = ir, ic 
              state = ' FREE' 
              IF ( ABS( prob%C( I )   - prob%C_l( i ) ) < ten * stopr )        &
                   state = 'LOWER' 
              IF ( ABS( prob%C( I )   - prob%C_u( i ) ) < ten * stopr )        &
                   state = 'UPPER' 
              IF ( ABS( prob%C_l( i ) - prob%C_u( i ) ) <       stopr )        &
                   state = 'EQUAL' 
              WRITE( out, 2130 ) CNAME( i ), STATE, prob%C( i ),               &
                   prob%C_l( i ), prob%C_u( i ), prob%Y( i ) 
           END DO
        END DO

        !  Compute the number of equality, fixed inequality and degenerate constraints

        mequal = 0 ; mfixed = 0 ; mdegen = 0 
        DO i = 1, m 
           IF ( ABS( prob%C( i ) - prob%C_l( i ) ) < ten * stopr .OR.         &
                ABS( prob%C( i ) - prob%C_u( i ) ) < ten * stopr ) THEN
              IF ( ABS( prob%C_l( i ) - prob%C_u( i ) ) < ten * stopr ) THEN 
                 mequal = mequal + 1 
              ELSE 
                 mfixed = mfixed + 1 
              END IF
              IF ( ABS( prob%Y( i ) ) < stopr ) mdegen = mdegen + 1 
           END IF
        END DO
     END IF
     WRITE( out, 2100 ) n, nfixed, ndegen 
     IF ( m > 0 ) THEN 
        WRITE( out, 2110 ) m, mequal, mdegen 
        IF ( m /= mequal ) WRITE( out, 2120 ) mfixed 
     END IF
     WRITE( out, 2030 ) qfval, iter, nfacts, factorization_integer,         &
          factorization_real 
  END IF

  times = times - time ; timet = timet - timeo
  WRITE( out, 2060 ) times + timet 
  qfvalt = qfval
  WRITE( out, 2070 ) pname 

  !  Compare the variants used so far

  WRITE( out, 2080 ) solv, iter, nfacts, qfvalt, status, times, timet,     &
       times + timet 

  DEALLOCATE( VNAME, CNAME, C )
  CLOSE( input  )

  WRITE ( out, 3000 ) pname, n, m, status, qfval, CPU( 1 ), CPU( 2 ) 

  STOP

  !  Non-executable statements

2000 FORMAT( ' .          .....  ..........',                                &
          '  ..........  ..........  .......... ' ) 
2010 FORMAT( /,' Stopping with inform%status = ', I3, / ) 
2020 FORMAT( /, ' Problem: ', A10 )
2030 FORMAT( /,' Final objective function value ', ES22.14, /,               &
          ' Total number of iterations = ',I6,' Number of factorizations = ', &
          I6, //, I10, ' integer and ', I10, ' real words required',          &
          ' for the factorization' ) 
2040 FORMAT( /,' Constraints : ', //, ' name       state    value  ',        &
          '  Lower bound Upper bound  Multiplier ' ) 
2050 FORMAT( 1X, A10, A6, 4ES12.4 ) 
2060 FORMAT( /, ' Total time = ', 0P, F12.2 ) 
2070 FORMAT( /, ' Problem: ', A10, //,                                       &
          '                                 objective',                   &
          '          < ------ time ----- > ', /,                          &
          ' Method  iterations   factors      value  ',                   &
          '   status setup    VE12   total', /,                           &
          ' ------  ----------   -------    ---------',                   &
          '   ------ -----    ----   -----  ' ) 
2080 FORMAT( A5, 2I10, 6X, ES12.4, I6, 0P, 3F8.2 ) 
2090 FORMAT( /,' Solution : ', //,'                              ',          &
          ' <------ Bounds ------> ', /           &
          ' name       state     value   ',                             &
          '    Lower       Upper      Dual ' ) 
2100 FORMAT( /, ' Of the ', I6, ' variables, ', 2X, I6,                      &
          ' are on bounds, and ', I6, ' are dual degenerate' ) 
2110 FORMAT( ' Of the ',I6,' constraints, ',I6,' are equations, and ',I6,    &
          ' are degenerate' ) 
2120 FORMAT( ' Of the inequality constraints ', I6, ' are on bounds' ) 
2130 FORMAT( 1X, A10, A6, 4ES12.4 ) 
2150 FORMAT( ' Allocation error, variable ', A8, ' status = ', I6 )
2160 FORMAT( ' IOSTAT = ', I6, ' when opening file ', A9, '. Stopping ' )
2180 FORMAT( A10 )
2190 FORMAT( A10, 2I7, 3I6, ES13.4, I6, 0P, F8.2 ) 
3000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //                    &
          ,' Code used               :  VE12',    /                           &
          ,' Problem                 :  ', A10,    /                          &
          ,' # variables             =      ', I10 /                          &
          ,' # constraints           =      ', I10 /                          &
          ' Exit code               =      ', I10 /                          &
          ,' Final f                 = ', E15.7 /                             &
          ,' Set up time             =      ', 0P, F10.2, ' seconds' /        &
          ' Solve time              =      ', 0P, F10.2, ' seconds' //       &
          66('*') / )
  
CONTAINS

  SUBROUTINE REORDER_by_rows( nr, nc, nnz, A_row, A_col, la, A_val,      &
       A_ptr, lptr, IW, liw, iflag )

    !  Reorder a sparse matrix from arbitary coordinate order to row order

    !  Dummy arguments

    INTEGER, INTENT( IN ) :: nr, nc, nnz, la, lptr, liw
    INTEGER, INTENT( OUT ) :: iflag
    REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( la ) :: A_val
    INTEGER, INTENT( INOUT ), DIMENSION( nnz ) :: A_row, A_col
    INTEGER, INTENT( OUT ), DIMENSION( lptr ) :: A_ptr
    INTEGER, INTENT( OUT ), DIMENSION( liw ) :: IW

    !  Local variables

    INTEGER :: i, j, k, k1, k2, l, nzi, ie, iep, je, jep
    INTEGER :: loc, idup, iout, jout, nzout
    INTEGER :: lp = 6, mp = 6
    REAL ( KIND = wp ) :: ae, aep
    SAVE lp, mp

    ! Initialize data

    iflag = 0
    nzout = 0 ; iout = 0 ; jout = 0 ; idup = 0

    !  Check for faulty input data

    IF ( nr < 1 .OR. nc < 1 .OR. nnz < 1 ) THEN
       iflag = - 1
       IF ( LP > 0 ) THEN
          WRITE( lp, 2000 ) iflag
          WRITE( lp, "( 1X, ' nr, nc, or nnz is out of range', /, &
               &           ' nc = ',I6,' nr = ',I6,' nnz = ', I10 )" ) nr, nc, nnz
       END IF
       RETURN
    END IF

    IF ( la < nnz ) THEN
       iflag = - 2
       IF ( lp > 0 ) THEN
          WRITE( lp, 2000 ) iflag
          WRITE( lp, "( 1X, ' increase la from', I8, ' to at least ', I8 )" )&
               la, nnz
       END IF
       RETURN
    END IF

    IF ( liw < MAX( nr, nc ) + 1 ) THEN
       iflag = - 3
       IF ( lp > 0 ) THEN
          WRITE( lp, 2000) iflag
          WRITE( lp, 2020) 'liw ', liw, MAX( nr, nc ) + 1
       END IF
       RETURN
    END IF

    IF ( lptr < nr + 1 ) THEN
       iflag = - 4
       IF ( lp > 0 ) THEN
          WRITE( lp, 2000 ) iflag
          WRITE( lp, 2020 ) 'lptr', lptr, nc + 1
       END IF
       RETURN
    END IF

    !  Record the number of column and row indices out of order in iout and jout
    !  and then remove them from consideration

    DO k = 1,nnz
       i = A_row( k ); j = A_col( k )
       IF ( i > nr .OR. i < 1 ) THEN
          iout = iout + 1
       ELSE IF ( j > nc .OR. j < 1 ) THEN
          jout = jout + 1
       ELSE
          nzout = nzout + 1
          A_row( nzout ) = i ; A_col( nzout ) = j
          A_val( nzout ) = A_val( k )
       END IF
    END DO

    !  Inform the user if there has been faulty data

    IF ( iout > 0 ) THEN
       iflag = 2
       IF ( mp > 0 ) THEN
          WRITE( mp, 2010 ) iflag
          WRITE( mp, "( 1X, I6,' entries input in A_row were out of',        &
               &              ' range and have been ignored by the routine')" ) iout
       END IF
    END IF

    IF ( jout > 0 ) THEN
       iflag = 3
       IF ( mp > 0 ) THEN
          WRITE( mp, 2010 ) iflag
          WRITE( mp, "( 1X, I6,' entries input in A_col were out of',        &
               &              ' range and have been ignored by the routine')" ) jout
       END IF
    END IF

    !  If all the data is faulty, exit

    IF ( iout + jout == nnz ) THEN
       nzout = 0
       RETURN
    END IF

    !  nzout gives the number of nonzero entries following removals. Now sort the 
    !  pattern of a sparse matrix from arbitary order to row order. The
    !  order within each row is unimportant

    !  Record the number of elements in each row in IW

    IW( : nr + 1 ) = 0
    DO k = 1, nzout
       i = A_row( k )
       IW( i ) = IW( i ) + 1
    END DO

    !  Record the positions where each row would begin, in a compressed format 
    !  with the rows in natural order, in A_ptr and IW 

    A_ptr( 1 ) = 1
    DO i = 2, nr + 1
       A_ptr( i ) = IW( i - 1 ) + A_ptr( i - 1 )
       IW( i - 1 ) = A_ptr( i - 1 )
    END DO

    !  Reorder the elements into row order. Fill in each row from the front,
    !  and increase the pointer IW( k ) by 1 as a new entry is placed in row k 

    DO l = 1, nr
       DO k = IW( l ), A_ptr( l + 1 ) - 1
          ie = A_row( k )
          je = A_col( k )
          ae = A_val( k )
          DO j = 1, nzout
             IF ( ie == l ) EXIT
             loc = IW( ie )
             iep = A_row( loc )
             jep = A_col( loc )
             aep = A_val( loc )
             IW( ie ) = loc + 1
             A_row( loc ) = ie
             A_col( loc ) = je
             A_val( loc ) = ae
             ie = iep
             je = jep
             ae = aep
          END DO
          A_row( k ) = ie
          A_col( k ) = je
          A_val( k ) = ae
       END DO
    END DO

    !  Check for duplicates

    nzout = 0
    k1 = 1
    nzi = 0
    IW( : nc ) = 0
    DO i = 1, nr
       k2 = A_ptr( i + 1 ) - 1
       A_ptr( i + 1 ) = A_ptr( i )
       DO K = k1, k2
          j = A_col( k )
          IF ( IW( j ) <= nzi ) THEN
             nzout = nzout + 1
             A_col( nzout ) = j
             A_val( nzout ) = A_val( k )
             A_ptr( i + 1 ) = A_ptr( i + 1 ) + 1
             IW( j ) = nzout

             !  There is a duplicate in row i; sum the values

          ELSE
             idup = idup + 1
             A_val( IW( j ) ) = A_val( IW( j ) ) + A_val( k )
          END IF
       END DO
       k1 = k2 + 1
       nzi = nzout
    END DO
    IF ( idup > 0 ) THEN
       iflag = 1
       IF ( mp > 0 ) THEN
          WRITE( mp, 2010 ) iflag
          WRITE( mp, "( 1X, I6,' Duplicate entries input ' )" ) idup
       END IF
    END IF

    !  Non-executable statements

2000 FORMAT (/,' -- Error return from REORDER_by_rows -- iflag = ', I2 )
2010 FORMAT (/,' -- Warning message from REORDER_by_rows -- iflag = ', I2)
2020 FORMAT (1X, ' increase', A4, ' from', I8,' to at least ', I8 )

    !  End of subroutine REORDER_by_rows
    
  END SUBROUTINE REORDER_by_rows

  !-*-*-*-*-*-*-*-*-*-  C P U _ T I M E  S U B R O U T I N E   -*-*-*-*-*-*-*-*

  SUBROUTINE CPU_TIME( time )

    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Return the cpu time (to be replaced in fortran 95 by a Fortran intrinsic)
    !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    !  Dummy arguments

    REAL, INTENT( OUT ) :: time

    !  Local variables

    REAL :: ZA02AS, dum
    time = ZA02AS( dum )

    !  End of subroutine CPU_TIME

  END SUBROUTINE CPU_TIME


  !  End of program VE12_main

END PROGRAM VE12_main

