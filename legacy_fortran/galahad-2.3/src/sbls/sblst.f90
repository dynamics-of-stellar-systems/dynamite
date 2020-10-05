! THIS VERSION: GALAHAD 2.2 - 24/04/2008 AT 08:00 GMT.
   PROGRAM GALAHAD_SBLS_EXAMPLE
   USE GALAHAD_SBLS_double                            ! double precision version
   IMPLICIT NONE
   INTEGER, PARAMETER :: wp = KIND( 1.0D+0 ) ! set precision
   TYPE ( SMT_type ) :: H, A, C
   TYPE ( SBLS_data_type ) :: data
   TYPE ( SBLS_control_type ) :: control        
   TYPE ( SBLS_inform_type ) :: info
   REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: SOL
   INTEGER :: n, m, h_ne, a_ne, c_ne, prec, preconditioner, factorization
   INTEGER :: data_storage_type, i, tests, status, scratch_out = 56
   INTEGER :: j, l, smt_stat
   REAL ( KIND = wp ) :: val

   IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
   CALL SMT_put( C%type, 'COORDINATE', smt_stat ) ; C%ne = 0
   ALLOCATE( C%val( C%ne ), C%row( C%ne ), C%col( C%ne ) )

   n = 3 ; m = 2 ; h_ne = 4 ; a_ne = 4 
   ALLOCATE( H%ptr( n + 1 ), A%ptr( m + 1 ) )

!  ================
!  error exit tests
!  ================

   WRITE( 6, "( /, ' error exit tests ' )" )
   WRITE( 6, "( /, ' test   status' )" )

   CALL SBLS_initialize( data, control )
   control%print_level = 0
!  control%print_level = 1

!  tests for status = - 3 and - 15

   DO l = 1, 2
     IF ( l == 1 ) THEN
       status = 3
     ELSE
       status = 15
     END IF
     ALLOCATE( H%val( h_ne ), H%col( h_ne ) )
     ALLOCATE( A%val( a_ne ), A%col( a_ne ) )
     IF ( ALLOCATED( H%type ) ) DEALLOCATE( H%type )
     CALL SMT_put( H%type, 'SPARSE_BY_ROWS', smt_stat ) 
     H%val = (/ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp /)
     H%col = (/ 1, 2, 3, 1 /)
     H%ptr = (/ 1, 2, 3, 5 /)
     IF ( ALLOCATED( A%type ) ) DEALLOCATE( A%type )
     CALL SMT_put( A%type, 'SPARSE_BY_ROWS', smt_stat ) 
     A%val = (/ 2.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)
     A%col = (/ 1, 2, 2, 3 /)
     A%ptr = (/ 1, 3, 5 /)

     IF ( status == 3 ) THEN
       n = 0 ; m = - 1
     ELSE
       n = 3 ; m = 2
       A%val = (/ 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)
       A%col = (/ 1, 2, 1, 2 /)
       control%preconditioner = 2
       control%perturb_to_make_definite = .FALSE.
     END IF

     CALL SBLS_form_and_factorize( n, m, H, A, C, data, control, info )
     WRITE( 6, "( I5, I9 )" ) status, info%status
     DEALLOCATE( H%val, H%col )
     DEALLOCATE( A%val, A%col )
   END DO

   CALL SBLS_terminate( data, control, info )
   DEALLOCATE( H%ptr, A%ptr )
   DEALLOCATE( C%val, C%row, C%col )

!  =====================================
!  basic test of various storage formats
!  =====================================

   WRITE( 6, "( /, ' basic tests of storage formats ' )" )

   n = 3 ; m = 2 ; h_ne = 4 ; a_ne = 3 ; c_ne = 3
   ALLOCATE( H%ptr( n + 1 ), A%ptr( m + 1 ), C%ptr( m + 1 ), SOL( n + m ) )

   DO prec = 1, 8
!    WRITE( 6, "( /, 8X, 'storage prec fact  new   status residual     value' )")
     WRITE( 6, "( /, 8X, 'storage prec fact  new   status residual' )" )
     SELECT CASE( prec)
     CASE( 1 : 4 ) ; preconditioner = prec
     CASE( 5 ) ; preconditioner = 11
     CASE( 6 ) ; preconditioner = 12
     CASE( 7 ) ; preconditioner = -1
     CASE( 8 ) ; preconditioner = -2
     END SELECT
     DO factorization = 1, 2
       DO data_storage_type = -3, 0
         CALL SBLS_initialize( data, control )
!        control%print_level = 1
         control%preconditioner = preconditioner
         control%factorization = factorization
         control%get_norm_residual = .TRUE.
         IF ( data_storage_type == 0 ) THEN   ! sparse co-ordinate storage
           IF ( ALLOCATED( H%type ) ) DEALLOCATE( H%type )
           CALL SMT_put( H%type, 'COORDINATE', smt_stat )  ; H%ne = h_ne
           ALLOCATE( H%val( h_ne ), H%row( h_ne ), H%col( h_ne ) )
           H%row = (/ 1, 2, 3, 3 /)
           H%col = (/ 1, 2, 3, 1 /)
           IF ( ALLOCATED( A%type ) ) DEALLOCATE( A%type )
           CALL SMT_put( A%type, 'COORDINATE', smt_stat )  ; A%ne = a_ne
           ALLOCATE( A%val( a_ne ), A%row( a_ne ), A%col( a_ne ) )
           A%row = (/ 1, 1, 2 /)
           A%col = (/ 1, 2, 3 /)
           IF ( preconditioner > 0 ) THEN
             IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
             CALL SMT_put( C%type, 'COORDINATE', smt_stat ) ; C%ne = c_ne
             ALLOCATE( C%val( c_ne ), C%row( c_ne ), C%col( c_ne ) )
             C%row = (/ 1, 2, 2 /)
             C%col = (/ 1, 1, 2 /)
           ELSE
             IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
             CALL SMT_put( C%type, 'COORDINATE', smt_stat ) ; C%ne = 0
             ALLOCATE( C%val( C%ne ), C%row( C%ne ), C%col( C%ne ) )
           END IF
         ELSE IF ( data_storage_type == - 1 ) THEN ! sparse row-wise storage
           IF ( ALLOCATED( H%type ) ) DEALLOCATE( H%type )
           CALL SMT_put( H%type, 'SPARSE_BY_ROWS', smt_stat )
           IF ( ALLOCATED( A%type ) ) DEALLOCATE( A%type )
           CALL SMT_put( A%type, 'SPARSE_BY_ROWS', smt_stat )
           ALLOCATE( H%val( h_ne ), H%col( h_ne ) )
           ALLOCATE( A%val( a_ne ), A%col( a_ne ) )
           H%col = (/ 1, 2, 3, 1 /)
           H%ptr = (/ 1, 2, 3, 5 /)
           A%col = (/ 1, 2, 3 /)
           A%ptr = (/ 1, 3, 4 /)
           IF ( preconditioner > 0 ) THEN
             IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
             CALL SMT_put( C%type, 'SPARSE_BY_ROWS', smt_stat )
             ALLOCATE( C%val( c_ne ), C%col( c_ne ) )
             C%col = (/ 1, 1, 2 /)
             C%ptr = (/ 1, 2, 4 /)
           ELSE
             IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
             CALL SMT_put( C%type, 'COORDINATE', smt_stat ) ; C%ne = 0
             ALLOCATE( C%val( C%ne ), C%row( C%ne ), C%col( C%ne ) )
           END IF
         ELSE IF ( data_storage_type == - 2 ) THEN ! dense storage
           IF ( ALLOCATED( H%type ) ) DEALLOCATE( H%type )
           CALL SMT_put( H%type, 'DENSE', smt_stat )
           IF ( ALLOCATED( A%type ) ) DEALLOCATE( A%type )
           CALL SMT_put( A%type, 'DENSE', smt_stat )
           ALLOCATE( H%val( n * ( n + 1 ) / 2 ) )
           ALLOCATE( A%val( n * m ) )
           IF ( preconditioner > 0 ) THEN
             IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
             CALL SMT_put( C%type, 'DENSE', smt_stat )
             ALLOCATE( C%val( m * ( m + 1 ) / 2 ) )
           ELSE
             IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
             CALL SMT_put( C%type, 'COORDINATE', smt_stat ) ; C%ne = 0
             ALLOCATE( C%val( C%ne ), C%row( C%ne ), C%col( C%ne ) )
           END IF
         ELSE IF ( data_storage_type == - 3 ) THEN ! dense storage
           IF ( ALLOCATED( H%type ) ) DEALLOCATE( H%type )
           CALL SMT_put( H%type, 'DIAGONAL', smt_stat )
           IF ( ALLOCATED( A%type ) ) DEALLOCATE( A%type )
           CALL SMT_put( A%type, 'DENSE', smt_stat )
           ALLOCATE( H%val( n ) )
           ALLOCATE( A%val( n * m ) )
           IF ( preconditioner > 0 ) THEN
             IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
             CALL SMT_put( C%type, 'DIAGONAL', smt_stat )
             ALLOCATE( C%val( m ) )
           ELSE
             IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
             CALL SMT_put( C%type, 'COORDINATE', smt_stat ) ; C%ne = 0
             ALLOCATE( C%val( C%ne ), C%row( C%ne ), C%col( C%ne ) )
           END IF
         END IF

!  test with new and existing data

         DO i = 0, 2
           control%new_a = 2 - i
           control%new_h = 2 - i
           control%new_c = 2 - i
           IF ( data_storage_type == 0 ) THEN     ! sparse co-ordinate storage
             H%val = (/ 1.0_wp, 2.0_wp, 3.0_wp, 1.0_wp /)
             A%val = (/ 2.0_wp, 1.0_wp, 1.0_wp /)
             IF ( preconditioner > 0 ) C%val = (/ 4.0_wp, 1.0_wp, 2.0_wp /)
           ELSE IF ( data_storage_type == - 1 ) THEN  !  sparse row-wise storage
             H%val = (/ 1.0_wp, 2.0_wp, 3.0_wp, 1.0_wp /)
             A%val = (/ 2.0_wp, 1.0_wp, 1.0_wp /)
             IF ( preconditioner > 0 ) C%val = (/ 4.0_wp, 1.0_wp, 2.0_wp /)
           ELSE IF ( data_storage_type == - 2 ) THEN    !  dense storage
             H%val = (/ 1.0_wp, 0.0_wp, 2.0_wp, 1.0_wp, 0.0_wp, 3.0_wp /)
             A%val = (/ 2.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp /)
             IF ( preconditioner > 0 ) C%val = (/ 4.0_wp, 1.0_wp, 2.0_wp /)
           ELSE IF ( data_storage_type == - 3 ) THEN    !  dense storage
             H%val = (/ 1.0_wp, 1.0_wp, 2.0_wp /)
             A%val = (/ 2.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp /)
             IF ( preconditioner > 0 ) C%val = (/ 4.0_wp, 2.0_wp /)
           END IF
           CALL SBLS_form_and_factorize( n, m, H, A, C, data, control, info )
           IF ( info%status < 0 ) THEN
             WRITE( 6, "( A15, 3I5, I9 )" ) SMT_get( H%type ),                 &
               preconditioner, factorization, control%new_h, info%status
             CYCLE
           END IF
           SOL( : n ) = (/ 0.0_wp, 2.0_wp, 0.0_wp /)
           SOL( n + 1 : ) = (/ 2.0_wp, 1.0_wp /)
           CALL SBLS_solve( n, m, A, C, data, control, info, SOL )
           IF ( info%status == 0 ) THEN
!            WRITE( 6, "( A15, 3I5, I9, A9, ES10.2 )" ) SMT_get( H%type ),     &
             WRITE( 6, "( A15, 3I5, I9, A9 )" ) SMT_get( H%type ),             &
               preconditioner, factorization, control%new_h, info%status,      &
               type_residual( info%norm_residual )
!              type_residual( info%norm_residual ), info%norm_residual
           ELSE
             WRITE( 6, "( A15, 3I5, I9 )" ) SMT_get( H%type ),                 &
               preconditioner, factorization, control%new_h, info%status
           END IF
         END DO

         CALL SBLS_terminate( data, control, info )
         IF ( data_storage_type == 0 ) THEN
           DEALLOCATE( H%val, H%row, H%col )
           DEALLOCATE( A%val, A%row, A%col )
         ELSE IF ( data_storage_type == - 1 ) THEN
           DEALLOCATE( H%val, H%col )
           DEALLOCATE( A%val, A%col )
         ELSE IF ( data_storage_type == - 2 ) THEN
           DEALLOCATE( H%val )
           DEALLOCATE( A%val )
         ELSE IF ( data_storage_type == - 3 ) THEN
           DEALLOCATE( H%val )
           DEALLOCATE( A%val )
         END IF
         IF ( preconditioner > 0 ) THEN
           IF ( data_storage_type == 0 ) THEN
             DEALLOCATE( C%val, C%row, C%col )
           ELSE IF ( data_storage_type == - 1 ) THEN
             DEALLOCATE( C%val, C%col )
           ELSE IF ( data_storage_type == - 2 ) THEN
             DEALLOCATE( C%val )
           ELSE IF ( data_storage_type == - 3 ) THEN
             DEALLOCATE( C%val )
           END IF
         ELSE
           DEALLOCATE( C%val, C%row, C%col )
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE( SOL )
   DEALLOCATE( H%ptr, A%ptr )

!  =============================
!  basic test of various options
!  =============================

   WRITE( 6, "( /, ' basic tests of options ' )" )

   n = 2 ; m = 1 ; h_ne = 2 ; a_ne = 2 
   ALLOCATE( H%ptr( n + 1 ), A%ptr( m + 1 ) )
   ALLOCATE( SOL( n + m ) )

   IF ( ALLOCATED( C%type ) ) DEALLOCATE( C%type )
   CALL SMT_put( C%type, 'COORDINATE', smt_stat ) ; C%ne = 0
   ALLOCATE( C%val( C%ne ), C%row( C%ne ), C%col( C%ne ) )

   SOL( : n ) = (/ 0.0_wp, 0.0_wp /)
   SOL( n + 1 : ) = (/ 1.0_wp /)

   ALLOCATE( H%val( h_ne ), H%row( 0 ), H%col( h_ne ) )
   ALLOCATE( A%val( a_ne ), A%row( 0 ), A%col( a_ne ) )
   IF ( ALLOCATED( H%type ) ) DEALLOCATE( H%type )
   CALL SMT_put( H%type, 'SPARSE_BY_ROWS', smt_stat )
   H%col = (/ 1, 2 /)
   H%ptr = (/ 1, 2, 3 /)
   IF ( ALLOCATED( A%type ) ) DEALLOCATE( A%type )
   CALL SMT_put( A%type, 'SPARSE_BY_ROWS', smt_stat )
   A%col = (/ 1, 2 /)
   A%ptr = (/ 1, 3 /)
   CALL SBLS_initialize( data, control )
   control%get_norm_residual = .TRUE.
   
!  test with new and existing data

   WRITE( 6, "( /, ' test   status residual' )" )
   tests = 22
   DO i = 0, tests
     IF ( i == 0 ) THEN
       control%preconditioner = 0
     ELSE IF ( i == 1 ) THEN
       control%preconditioner = 1
     ELSE IF ( i == 2 ) THEN
       control%preconditioner = 2
     ELSE IF ( i == 3 ) THEN
       control%preconditioner = 3
     ELSE IF ( i == 4 ) THEN
       control%preconditioner = 4
     ELSE IF ( i == 5 ) THEN     
       control%factorization = - 1
     ELSE IF ( i == 6 ) THEN     
       control%factorization = 1
     ELSE IF ( i == 7 ) THEN     
       control%max_col = 0
     ELSE IF ( i == 8 ) THEN     
       control%factorization = 2
       control%preconditioner = 0
     ELSE IF ( i == 9 ) THEN
!      control%print_level = 2
       control%preconditioner = 11
     ELSE IF ( i == 10 ) THEN
       control%preconditioner = 12
     ELSE IF ( i == 11 ) THEN
       control%preconditioner = -1
     ELSE IF ( i == 12 ) THEN
       control%preconditioner = -2
     ELSE IF ( i == 13 ) THEN

     ELSE IF ( i == 14 ) THEN
       control%max_col = 5
!      control%primal = .TRUE.       
     ELSE IF ( i == 15 ) THEN
       control%max_col = 75
     ELSE IF ( i == 16 ) THEN
!      control%feasol = .FALSE.
     ELSE IF ( i == 16 ) THEN

     ELSE IF ( i == 17 ) THEN

     ELSE IF ( i == 18 ) THEN

     ELSE IF ( i == 19 ) THEN

     ELSE IF ( i == 20 ) THEN

     ELSE IF ( i == 21 ) THEN

     ELSE IF ( i == 22 ) THEN

     END IF

     H%val = (/ 1.0_wp, 1.0_wp /)
     A%val = (/ 1.0_wp, 1.0_wp /)

!    control%print_level = 4
     CALL SBLS_form_and_factorize( n, m, H, A, C, data, control, info )
     CALL SBLS_solve( n, m, A, C, data, control, info, SOL )
!    write(6,"('x=', 2ES12.4)") X
     IF ( info%status == 0 ) THEN
       WRITE( 6, "( I5, I9, A9 )" )                                            &
         i, info%status, type_residual( info%norm_residual )
     ELSE
       WRITE( 6, "( I5, I9 )" ) i, info%status
     END IF
   END DO
   CALL SBLS_terminate( data, control, info )

   DEALLOCATE( H%val, H%row, H%col )
   DEALLOCATE( A%val, A%row, A%col )
   DEALLOCATE( H%ptr, A%ptr )
   DEALLOCATE( SOL )

!  ============================
!  full test of generic problem
!  ============================

   WRITE( 6, "( /, ' full test of generic problems ' )" )
   WRITE( 6, "( /, ' test   status residual' )" )

   n = 14 ; m = 8 ; h_ne = 28 ; a_ne = 27
   ALLOCATE( SOL( n + m ) )
   ALLOCATE( H%ptr( n + 1 ), A%ptr( m + 1 ) )
   ALLOCATE( H%val( h_ne ), H%row( h_ne ), H%col( h_ne ) )
   ALLOCATE( A%val( a_ne ), A%row( a_ne ), A%col( a_ne ) )
   IF ( ALLOCATED( H%type ) ) DEALLOCATE( H%type )
   CALL SMT_put( H%type, 'COORDINATE', smt_stat )
   IF ( ALLOCATED( A%type ) ) DEALLOCATE( A%type )
   CALL SMT_put( A%type, 'COORDINATE', smt_stat )
   H%ne = h_ne ; A%ne = a_ne 
   SOL( : n ) = (/ 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp,     &
                   0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp /) 
   SOL( n + 1 : ) = (/ 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,         &
                       1.0_wp, 2.0_wp /)
   H%val = (/ 1.0_wp, 1.0_wp, 2.0_wp, 2.0_wp, 3.0_wp, 3.0_wp,                  &
                4.0_wp, 4.0_wp, 5.0_wp, 5.0_wp, 6.0_wp, 6.0_wp,                &
                7.0_wp, 7.0_wp,                                                &
                20.0_wp, 20.0_wp, 20.0_wp, 20.0_wp, 20.0_wp, 20.0_wp, 20.0_wp, &
                20.0_wp, 20.0_wp, 20.0_wp, 20.0_wp, 20.0_wp, 20.0_wp, 20.0_wp /)
   H%row = (/ 1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14,                   &
                1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14 /)
   H%col = (/ 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7,                        &
                1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14 /)
   A%val = (/ 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                          &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                        &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                        &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)
   A%row = (/ 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5,                  &
                6, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8 /)
   A%col = (/ 1, 3, 5, 1, 2, 1, 2, 3, 4, 5, 6, 5, 6, 2, 4, 6,                  &
                8, 10, 12, 8, 9, 8, 10, 11, 12, 13, 14 /)

   CALL SBLS_initialize( data, control )
   control%get_norm_residual = .TRUE.
   control%print_level = 101
   control%itref_max = 3
   control%out = scratch_out
   control%error = scratch_out
!  control%print_level = 1
!  control%out = 6
!  control%error = 6
   OPEN( UNIT = scratch_out, STATUS = 'SCRATCH' )
   CALL SBLS_form_and_factorize( n, m, H, A, C, data, control, info )
   CALL SBLS_solve( n, m, A, C, data, control, info, SOL )
   CLOSE( UNIT = scratch_out )
   IF ( info%status == 0 ) THEN
     WRITE( 6, "( I5, I9, A9 )" )                                              &
       1, info%status, type_residual( info%norm_residual )
   ELSE
     WRITE( 6, "( I5, I9 )" ) 1, info%status
   END IF
   CALL SBLS_terminate( data, control, info )
   DEALLOCATE( H%val, H%row, H%col )
   DEALLOCATE( A%val, A%row, A%col )
   DEALLOCATE( H%ptr, A%ptr )
   DEALLOCATE( SOL )

!  Second problem

   n = 14 ; m = 8 ; h_ne = 14 ; a_ne = 27
   ALLOCATE( SOL( n + m ) )
   ALLOCATE( H%ptr( n + 1 ), A%ptr( m + 1 ) )
   ALLOCATE( H%val( h_ne ), H%row( h_ne ), H%col( h_ne ) )
   ALLOCATE( A%val( a_ne ), A%row( a_ne ), A%col( a_ne ) )
   IF ( ALLOCATED( H%type ) ) DEALLOCATE( H%type )
   CALL SMT_put( H%type, 'COORDINATE', smt_stat )
   IF ( ALLOCATED( A%type ) ) DEALLOCATE( A%type )
   CALL SMT_put( A%type, 'COORDINATE', smt_stat )
   H%ne = h_ne ; A%ne = a_ne 
   SOL( : n ) = (/ 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp,    &
                   0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp /) 
   SOL( n + 1 : ) = (/ 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp,                        &
                       0.0_wp, 0.0_wp, 1.0_wp, 2.0_wp /)
   H%val = (/ 1.0_wp, 1.0_wp, 2.0_wp, 2.0_wp, 3.0_wp, 3.0_wp,                 &
                4.0_wp, 4.0_wp, 5.0_wp, 5.0_wp, 6.0_wp, 6.0_wp,               &
                7.0_wp, 7.0_wp /)
   H%row = (/ 1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14 /)
   H%col = (/ 1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14 /)
   A%val = (/ 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                         &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                       &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                       &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,               &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)
   A%row = (/ 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5,                 &
                6, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8 /)
   A%col = (/ 1, 3, 5, 1, 2, 1, 2, 3, 4, 5, 6, 5, 6, 2, 4, 6,                 &
                8, 10, 12, 8, 9, 8, 10, 11, 12, 13, 14 /)

   CALL SBLS_initialize( data, control )
   control%get_norm_residual = .TRUE.
!  control%print_level = 1
   control%preconditioner = -2 ; control%factorization = 2
   CALL SBLS_form_and_factorize( n, m, H, A, C, data, control, info )
   CALL SBLS_solve( n, m, A, C, data, control, info, SOL )
   IF ( info%status == 0 ) THEN
     WRITE( 6, "( I5, I9, A9 )" )                                              &
       2, info%status, type_residual( info%norm_residual )
   ELSE
     WRITE( 6, "( I5, I9 )" ) 2, info%status
   END IF
   CALL SBLS_terminate( data, control, info )
   DEALLOCATE( H%val, H%row, H%col )
   DEALLOCATE( A%val, A%row, A%col )
   DEALLOCATE( H%ptr, A%ptr )
   DEALLOCATE( SOL )

!  Third problem

!  WRITE( 25, "( /, ' third problem ', / )" )

   n = 14 ; m = 8 ; h_ne = 14 ; a_ne = 26
   ALLOCATE( SOL( n + m ) )
   ALLOCATE( H%ptr( n + 1 ), A%ptr( m + 1 ) )
   ALLOCATE( H%val( h_ne ), H%row( h_ne ), H%col( h_ne ) )
   ALLOCATE( A%val( a_ne ), A%row( a_ne ), A%col( a_ne ) )
   IF ( ALLOCATED( H%type ) ) DEALLOCATE( H%type )
   CALL SMT_put( H%type, 'COORDINATE', smt_stat )
   IF ( ALLOCATED( A%type ) ) DEALLOCATE( A%type )
   CALL SMT_put( A%type, 'COORDINATE', smt_stat )
   H%ne = h_ne ; A%ne = a_ne 
   SOL( : n ) = (/ 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp,    &
            0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 2.0_wp, 0.0_wp, 2.0_wp /) 
   SOL( n + 1 : ) = (/ 0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,                &
                       0.0_wp, 1.0_wp, 2.0_wp /)
   A%val = (/ 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                         &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                       &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                       &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp,                       &
                1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp /)
   A%row = (/ 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5,                    &
                6, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8 /)  
   A%col = (/ 1, 3, 5, 1, 2, 1, 2, 3, 4, 5, 6, 5, 6, 2, 6,                    &
                8, 10, 12, 8, 9, 8, 10, 11, 12, 13, 14 /)
   H%val = (/ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp,         &
                1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp /)
   H%row = (/ 1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14 /)
   H%col = (/ 1, 8, 2, 9, 3, 10, 4, 11, 5, 12, 6, 13, 7, 14 /)
!  SOL = 0.0_wp
!  DO l = 1, a_ne
!    i = A%row( l )
!    j = A%col( l )
!    val = A%val( l )
!    SOL( n + i ) = SOL( n + i ) + val
!    SOL( j ) = SOL( j ) + val
!  END DO
!  DO l = 1, h_ne
!    i = H%row( l )
!    j = H%col( l )
!    val = H%val( l )
!    SOL( i ) = SOL( i ) + val
!    IF ( i /= j ) SOL( j ) = SOL( j ) + val
!  END DO
   CALL SBLS_initialize( data, control )
   control%get_norm_residual = .TRUE.
!  control%print_level = 2 ; control%out = 6 ; control%error = 6
   control%preconditioner = 2 ; control%factorization = 1
   control%min_diagonal = 1.0_wp
   CALL SBLS_form_and_factorize( n, m, H, A, C, data, control, info )
   CALL SBLS_solve( n, m, A, C, data, control, info, SOL )
!     WRITE( 25,"( ' solution ', /, ( 5ES24.16 ) )" ) SOL
   IF ( info%status == 0 ) THEN
     WRITE( 6, "( I5, I9, A9 )" )                                             &
       3, info%status, type_residual( info%norm_residual )
   ELSE
     WRITE( 6, "( I5, I9 )" ) 3, info%status
   END IF
   CALL SBLS_terminate( data, control, info )
   DEALLOCATE( H%val, H%row, H%col )
   DEALLOCATE( A%val, A%row, A%col )
   DEALLOCATE( SOL )

   WRITE( 6, "( /, ' end of tests' )" )

   CONTAINS
     CHARACTER ( len = 8 ) FUNCTION type_residual( residual )
     REAL ( KIND = wp ) :: residual
     REAL, PARAMETER :: ten = 10.0_wp
     REAL, PARAMETER :: tiny = ten ** ( - 12 )
     REAL, PARAMETER :: small = ten ** ( - 8 )
     REAL, PARAMETER :: medium = ten ** ( - 4 )
     REAL, PARAMETER :: large = 1.0_wp
     IF ( ABS( residual ) < tiny ) THEN
       type_residual = '    tiny'
     ELSE IF ( ABS( residual ) < small ) THEN
       type_residual = '   small'
     ELSE IF ( ABS( residual ) < medium ) THEN
       type_residual = '  medium'
     ELSE IF ( ABS( residual ) < large ) THEN
       type_residual = '   large'
     ELSE
       type_residual = '    huge'
     END IF
     RETURN
     END FUNCTION type_residual

   END PROGRAM GALAHAD_SBLS_EXAMPLE


