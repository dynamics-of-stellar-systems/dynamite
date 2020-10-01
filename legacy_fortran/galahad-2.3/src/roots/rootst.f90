! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.
   PROGRAM GALAHAD_ROOTS_test_deck
   USE GALAHAD_ROOTS_double                       ! double precision version
   IMPLICIT NONE
   INTEGER, PARAMETER :: wp = KIND( 1.0D+0 ) ! set precision
   REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp

   INTEGER :: order, type, nroots, status
   REAL ( KIND = wp ) :: tol, A( 0 : 4 ), ROOTS( 4 )
   REAL ( KIND = wp ) :: A1( 0 : 5 ), A2( 0: 2 ), ROOTS2( 1 )

   tol = EPSILON( one ) ** 0.75

   DO order = 2, 4
     DO type = 1, 2
       IF ( order == 2 ) THEN
         A( 0 ) = 2.01_wp
         A( 1 ) = - 3.01_wp
         A( 2 ) = 1.01_wp
         IF ( type == 2 ) A( 2 ) = 0.0_wp
       ELSE IF ( order == 3 ) THEN
         A( 0 ) = - 6.01_wp
         A( 1 ) = 11.01_wp
         A( 2 ) = -6.01_wp
         A( 3 ) = 1.01_wp
         IF ( type == 2 ) A( 3 ) = 0.0_wp
       ELSE
         IF ( type == 1 ) THEN
           A( 0 ) = 24.001_wp
           A( 1 ) = -50.001_wp
           A( 2 ) = 35.001_wp
           A( 3 ) = -10.001_wp
           A( 4 ) = 1.001_wp
         ELSE
           A( 0 ) = 1.00_wp
           A( 1 ) = -4.00_wp
           A( 2 ) = 6.00_wp
           A( 3 ) = -4.00_wp
           A( 4 ) = 1.00_wp
         END IF
       END IF

       IF ( type == 1 ) THEN
         IF ( order == 2 ) THEN
           WRITE( 6, "( /, ' Quadratic ' )" )
           CALL ROOTS_quadratic( A( 0 ), A( 1 ), A( 2 ), tol,                  &
             nroots, ROOTS( 1 ), ROOTS( 2 ) )
         ELSE IF ( order == 3 ) THEN
           WRITE( 6, "( /, ' Cubic ' )" )
           CALL ROOTS_cubic( A( 0 ), A( 1 ), A( 2 ), A( 3 ), tol,              &
             nroots, ROOTS( 1 ), ROOTS( 2 ), ROOTS( 3 ) )
         ELSE
           WRITE( 6, "( /, ' Quartic ' )" )
           CALL ROOTS_quartic( A( 0 ), A( 1 ), A( 2 ), A( 3 ), A( 4 ), tol,    &
             nroots, ROOTS( 1 ), ROOTS( 2 ), ROOTS( 3 ), ROOTS( 4 ) )
         END IF
       ELSE
         IF ( order == 2 ) THEN
           WRITE( 6, "( /, ' Quadratic ' )" )
           CALL ROOTS_solve( A( 0 : order ), tol, nroots,                     &
             ROOTS( : order ), status )
         ELSE IF ( order == 3 ) THEN
           WRITE( 6, "( /, ' Cubic ' )" )
           CALL ROOTS_solve( A( 0 : order ), tol, nroots,                     &
             ROOTS( : order ), status )
         ELSE
           WRITE( 6, "( /, ' Quartic ' )" )
           CALL ROOTS_solve( A( 0 : order ), tol, nroots,                     &
             ROOTS( : order ), status )
         END IF
       END IF
       IF ( nroots == 0 ) THEN
         WRITE( 6, "( ' no real roots ' )" )
       ELSE IF ( nroots == 1 ) THEN
         WRITE( 6, "( ' 1 real root ' )" )
       ELSE IF ( nroots == 2 ) THEN
         WRITE( 6, "( ' 2 real roots ' )" )
       ELSE IF ( nroots == 3 ) THEN
         WRITE( 6, "( ' 3 real roots ' )" )
       ELSE IF ( nroots == 4 ) THEN
         WRITE( 6, "( ' 4 real roots ' )" )
       END IF
       IF ( nroots /= 0 ) WRITE( 6, "( ' roots: ', 4ES12.4 )" ) ROOTS( : nroots )
       
     END DO
   END DO

!  Test for error exits

   WRITE(6,"( /, ' Tests for error exits ' )" )
   A1( 0 : 4 ) = A ; A1( 5 ) = 1.0_wp
   CALL ROOTS_solve( A1, tol, nroots, ROOTS, status )
   WRITE(6,"( ' Test 3: exit status ', I0 )" ) status
   A2 = A(  0 : 2 )
   CALL ROOTS_solve( A2, tol, nroots, ROOTS2, status )
   WRITE(6,"( ' Test 4: exit status ', I0 )" ) status

   STOP

   END PROGRAM GALAHAD_ROOTS_test_deck

