! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.
   USE GALAHAD_ROOTS_double         ! double precision version
   IMPLICIT NONE
   INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )  ! set precision
   REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
   INTEGER :: degree, nroots, status
   REAL ( KIND = wp ) :: tol, A( 0 : 4 ), ROOTS( 4 )
   tol = EPSILON( one ) ** 0.75               ! accuracy requested
   DO degree = 2, 4                           ! polynomials of degree 2 to 4
     IF ( degree == 2 ) THEN
       A( 0 ) = 2.0_wp
       A( 1 ) = - 3.0_wp
       A( 2 ) = 1.0_wp
       WRITE( 6, "( ' Quadratic ' )" )
       CALL ROOTS_solve( A( : degree ), tol, nroots, ROOTS( : degree ), status )
     ELSE IF ( degree == 3 ) THEN
       A( 0 ) = - 6.0_wp
       A( 1 ) = 11.0_wp
       A( 2 ) = -6.0_wp
       A( 3 ) = 1.0_wp
       WRITE( 6, "( /, ' Cubic ' )" )
       CALL ROOTS_solve( A( : degree ), tol, nroots, ROOTS( : degree ), status )
     ELSE
       A( 0 ) = 24.0_wp
       A( 1 ) = -50.0_wp
       A( 2 ) = 35.0_wp
       A( 3 ) = -10.0_wp
       A( 4 ) = 1.0_wp
       WRITE( 6, "( /, ' Quartic ' )" )
       CALL ROOTS_solve( A( : degree ), tol, nroots, ROOTS( : degree ), status )
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
     IF ( nroots /= 0 ) WRITE( 6, "( ' roots: ', 4ES10.2 )" ) ROOTS( : nroots )
   END DO
   END PROGRAM GALAHAD_ROOTS_EXAMPLE

