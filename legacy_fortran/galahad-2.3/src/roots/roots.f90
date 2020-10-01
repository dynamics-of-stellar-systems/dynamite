! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-*-*-  G A L A H A D _ L S Q P    M O D U L E  -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   started life (quadratic roots) in GALAHAD_LSQP ~ 2000
!   released with GALAHAD Version 2.0. April 27th 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_ROOTS_double

!     ---------------------------------------------------
!     |                                                 |
!     |  Find all the real roots of quadratic, cubic    |
!     |  and quartic polynomials with real coefficients |
!     |                                                 |
!     ---------------------------------------------------

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: ROOTS_solve, ROOTS_quadratic, ROOTS_cubic, ROOTS_quartic

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!------------------------------------
!   G e n e r i c   I n t e r f a c e
!------------------------------------

!     INTERFACE ROOTS_solve
!       MODULE PROCEDURE ROOTS_quadratic, ROOTS_cubic, ROOTS_quartic
!     END INTERFACE ROOTS_solve

!----------------------
!   P a r a m e t e r s
!----------------------

      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: two = 2.0_wp
      REAL ( KIND = wp ), PARAMETER :: three = 3.0_wp
      REAL ( KIND = wp ), PARAMETER :: four = 4.0_wp
      REAL ( KIND = wp ), PARAMETER :: point1 = 0.1_wp
      REAL ( KIND = wp ), PARAMETER :: quarter = 0.25_wp
      REAL ( KIND = wp ), PARAMETER :: onethird = one / three
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: twothirds = two / three
!     REAL ( KIND = wp ), PARAMETER :: pi = four * ATAN( 1.0_wp )
!     REAL ( KIND = wp ), PARAMETER :: magic = twothirds * pi
      REAL ( KIND = wp ), PARAMETER :: magic = 2.0943951023931953_wp  !! 2 pi / 3
      REAL ( KIND = wp ), PARAMETER :: infinity = HUGE( one ) * point1

      INTEGER, PARAMETER :: out = 6
      LOGICAL, PARAMETER :: debug = .FALSE.
!     LOGICAL, PARAMETER :: debug = .TRUE.

   CONTAINS

!-*-*-*-*-*-*-*-   R O O T S _ s o l v e   S U B R O U T I N E   -*-*-*-*-*-*-*-

      SUBROUTINE ROOTS_solve( A, tol, nroots, ROOTS, status )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Find all the real roots of a real polynomial sum_i>=0 a(i) x^i of degree <= 4

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: nroots, status
      REAL ( KIND = wp ), INTENT( IN ) :: tol
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( 0 : ) :: A
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: ROOTS 

!  Local variables

      INTEGER :: degree

!  Check input details for errors and consistency

      degree = UBOUND( A, 1 )
      IF ( degree < 0 .OR. degree > 4 ) THEN
        status = - 3 ; RETURN
      ELSE IF ( UBOUND( ROOTS, 1 ) < degree ) THEN
        status = - 4 ; RETURN
      END IF

!  The data appears to be correct

      SELECT CASE( degree )

      CASE ( 0 )

!  polynomials of degree 0

        nroots = 0

      CASE ( 1 )

!  polynomials of degree 1

        IF ( A( 1 ) == zero ) THEN
          IF ( A( 1 ) == zero ) THEN
            nroots = 1
            ROOTS( 1 ) = zero
          ELSE
            nroots = 0
          END IF
        ELSE
          nroots = 1
          ROOTS( 1 ) = - A( 0 ) / A( 1 )
        END IF

      CASE ( 2 )

!  polynomials of degree 2

        CALL ROOTS_quadratic( A( 0 ), A( 1 ), A( 2 ), tol,                     &
          nroots, ROOTS( 1 ), ROOTS( 2 ) ) 
      CASE ( 3 )

!  polynomials of degree 3

        CALL ROOTS_cubic( A( 0 ), A( 1 ), A( 2 ), A( 3 ), tol,                 &
          nroots, ROOTS( 1 ), ROOTS( 2 ), ROOTS( 3 ) ) 
      CASE ( 4 )

!  polynomials of degree 4

        CALL ROOTS_quartic( A( 0 ), A( 1 ), A( 2 ), A( 3 ), A( 4 ), tol,       &
          nroots, ROOTS( 1 ), ROOTS( 2 ), ROOTS( 3 ), ROOTS( 4 ) )
      END SELECT

      status = 0
      RETURN  

!  End of subroutine ROOTS_solve

      END SUBROUTINE ROOTS_solve

!-*-*-*-*-*-   R O O T S _ q u a d r a t i c  S U B R O U T I N E   -*-*-*-*-*-

      SUBROUTINE ROOTS_quadratic( a0, a1, a2, tol, nroots, root1, root2 ) 

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the number and values of real roots of the quadratic equation
! 
!                   a2 * x**2 + a1 * x + a0 = 0
!
!  where a0, a1 and a2 are real 
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: nroots 
      REAL ( KIND = wp ), INTENT( IN ) :: a2, a1, a0, tol 
      REAL ( KIND = wp ), INTENT( OUT ) :: root1, root2

!  Local variables

      REAL ( KIND = wp ) :: rhs, d 

      rhs = tol * a1 * a1 
      IF ( ABS( a0 * a2 ) > rhs ) THEN 
        root2 = a1 * a1 - four * a2 * a0 
        IF ( root2 < zero ) THEN 
          nroots = 0 ; root1 = zero ; root2 = zero 
        ELSE 
          d = - half * ( a1 + SIGN( SQRT( root2 ), a1 ) ) 
          nroots = 2 ; root1 = d / a2 ; root2 = a0 / d 
          IF ( root1 > root2 ) THEN 
            d = root1 ; root1 = root2 ; root2 = d 
          END IF 
        END IF 
      ELSE IF ( a2 * a2 > rhs ) THEN 
        nroots = 2 
        IF ( - a1 / a2 > zero ) THEN 
          root1 = zero ; root2 = - a1 / a2 
        ELSE 
          root1 = - a1 / a2 ; root2 = zero 
        END IF 
      ELSE IF ( a2 == zero .AND. a1 == zero ) THEN 
        IF ( a0 == zero ) THEN
          nroots = 1 ; root1 = zero ; root2 = zero 
        ELSE
          nroots = 0 ; root1 = zero ; root2 = zero 
        END IF
      ELSE IF ( a0 * a0 > rhs ) THEN 
        nroots = 1 ; root1 = - a0 / a1 ; root2 = zero 
      ELSE 
        nroots = 1 ; root1 = zero ; root2 = zero 
      END IF 
      RETURN  

!  End of subroutine ROOTS_quadratic

      END SUBROUTINE ROOTS_quadratic 

!-*-*-*-*-*-*-*-   R O O T S _ c u b i c  S U B R O U T I N E   -*-*-*-*-*-*-*-

      SUBROUTINE ROOTS_cubic( a0, a1, a2, a3, tol, nroots, root1, root2, root3 ) 

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the number and values of real roots of the cubicc equation
! 
!                a3 * x**3 + a2 * x**2 + a1 * x + a0 = 0
!
!  where a0, a1, a2 and a3 are real 
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: nroots
      REAL ( KIND = wp ), INTENT( IN ) :: a3, a2, a1, a0, tol 
      REAL ( KIND = wp ), INTENT( OUT ) :: root1, root2, root3

!  Local variables

      REAL ( KIND = wp ) :: a, b, c, e, f, x, y, z
      REAL ( KIND = wp ) :: c0, c1, c2, b0, b1, p, pprime

!  Check to see if the quartic is actually a cubic

      IF ( a3 == zero ) THEN
        CALL ROOTS_quadratic( a0, a1, a2, tol, nroots, root1, root2 ) 
        root3 = infinity
        RETURN
      END IF

!  Deflate the polnomial if the trailing coefficient is zero

      IF ( a0 == zero ) THEN
        root1 = zero
        CALL ROOTS_quadratic( a1, a2, a3, tol, nroots, root2, root3 ) 
        nroots = nroots + 1
        RETURN
      END IF

!  Use Littlewood's method

!  The cubic is general

      c2 = a2 / ( three * a3 )
      c1 = a1 / ( three * a3 )
      c0 = a0 / a3
      x = c1 - c2 * c2
      y = c0 - c2 * ( x + x + c1 )
      z = four * x ** 3 + y ** 2 
      IF ( z < zero ) THEN
        root1 = - two * SQRT( - x )
        y = y / ( root1 * x )
        x = root1
        y = ATAN2( SQRT( one - y ), SQRT( one + y  ) ) * twothirds
        IF ( c2 < zero ) y = y + magic
        root1 = x * COS( y ) - c2
      ELSE
        a = SQRT( z )
        b = half * ( ABS( y ) + a )
        c = b ** onethird

        IF ( c <= zero ) THEN
          root1 = - c2
          root2 = - c2
          root3 = - c2
          nroots = 3
          GO TO 900
        END IF

        c = c - ( c ** 3 - b ) / ( three * c * c )
        e = c * c + ABS( x )
        f = one/ ( ( x / c ) ** 2 + e )

        IF ( x >= zero ) THEN
          x = e / c
          z = y * f
        ELSE
          x = a * f
          z = SIGN( one, y ) * e / c
        END IF

!   There is one real root

        IF ( z * c2 < zero ) THEN
          nroots = 1
          root2 = half * z - c2
          root3 = half * SQRT( three ) * ABS( x )
          root1 = - c0 / ( root2 * root2 + root3 * root3 )
          GO TO 900
        END IF
        root1 = - z - c2
      END IF

!  Deflate cubic from optimal end

      b0 = - c0 / root1
      IF ( ABS( root1 ** 3 ) <= ABS( c0 ) ) THEN
        b1 = root1 + three * c2
      ELSE
        b1 = ( b0 - three * c1 ) / root1
      END IF

      x = b1 * b1 - four * b0

!  There are three real roots.

      IF ( x >= zero ) THEN
        root3 = - SIGN( half, b1 ) * ( SQRT( x ) + ABS( b1 ) )
        IF ( root3 /= zero ) THEN
          root2 = b0 / root3
        ELSE
          root2 = zero
        END IF

!  Reorder the roots

        IF ( root1 > root2 ) THEN
          a = root2
          root2 = root1
          root1 = a
        END IF

        IF ( root2 > root3 ) THEN
          a = root3
          IF ( root1 > root3 ) THEN
            a = root1
            root1 = root3
          END IF
          root3 = root2
          root2 = a
        END IF
        nroots = 3

!   There is one real root

      ELSE
        root2 = - half * b1
        root3 = half * SQRT( - x )
        nroots = 1
      END IF

  900 CONTINUE

!  Perfom a Newton iteration to ensure that the roots are accurate

      IF ( debug ) THEN
        IF ( nroots == 1 ) THEN
          WRITE( out, "( ' 1 real root ' )" )
        ELSE
          WRITE( out, "( ' 3 real roots ' )" )
        END IF
      END IF

      p = ( ( a3 * root1 + a2 ) * root1 + a1 ) * root1 + a0
      pprime = ( three * a3 * root1 + two * a2 ) * root1 + a1
      IF ( pprime /= zero ) THEN
        IF ( debug ) WRITE( out, 2000 ) 1, root1, p, - p / pprime 
        root1 = root1 - p / pprime
        p = ( ( a3 * root1 + a2 ) * root1 + a1 ) * root1 + a0
      END IF
      IF ( debug ) WRITE( out, 2010 ) 1, root1, p

      IF ( nroots == 3 ) THEN
        p = ( ( a3 * root2 + a2 ) * root2 + a1 ) * root2 + a0
        pprime = ( three * a3 * root2 + two * a2 ) * root2 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 2, root2, p, - p / pprime 
          root2 = root2 - p / pprime
          p = ( ( a3 * root2 + a2 ) * root2 + a1 ) * root2 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 2, root2, p

        p = ( ( a3 * root3 + a2 ) * root3 + a1 ) * root3 + a0
        pprime = ( three * a3 * root3 + two * a2 ) * root3 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 3, root3, p, - p / pprime 
          root3 = root3 - p / pprime
          p = ( ( a3 * root3 + a2 ) * root3 + a1 ) * root3 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 3, root3, p
      END IF

      RETURN

!  Non-executable statements

 2000 FORMAT( ' root ', I1, ': value = ', ES12.4, ' cubic = ', ES12.4,         &
              ' delta = ', ES12.4 )
 2010 FORMAT( ' root ', I1, ': value = ', ES12.4, ' cubic = ', ES12.4 )


!  End of subroutine ROOTS_cubic

      END SUBROUTINE ROOTS_cubic

!-*-*-*-*-*-*-   R O O T S _ q u a r t i c   S U B R O U T I N E   -*-*-*-*-*-*-

      SUBROUTINE ROOTS_quartic( a0, a1, a2, a3, a4, tol, nroots, root1, root2, &
                                root3, root4 ) 

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Find the number and values of real roots of the quartic equation
! 
!        a4 * x**4 + a3 * x**3 + a2 * x**2 + a1 * x + a0 = 0
!
!  where a0, a1, a2, a3 and a4 are real 
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: nroots
      REAL ( KIND = wp ), INTENT( IN ) :: a4, a3, a2, a1, a0, tol 
      REAL ( KIND = wp ), INTENT( OUT ) :: root1, root2, root3, root4

!  Local variables

      INTEGER :: type_roots, nrootsc
      REAL ( KIND = wp ) :: a, alpha, b, beta, c, d, delta, gamma, r
      REAL ( KIND = wp ) :: x1, xm, xmd, xn, xnd
      REAL ( KIND = wp ) :: d3, d2, d1, d0, b4, b3, b2, b1
      REAL ( KIND = wp ) :: rootc1, rootc2, rootc3, p, pprime

!  Check to see if the quartic is actually a cubic

      IF ( a4 == zero ) THEN
        CALL ROOTS_cubic( a0, a1, a2, a3, tol, nroots, root1, root2, root3 ) 
        root4 = infinity
        RETURN
      END IF

!  Use Ferrari's algorithm

!  Initialize 

      nroots = 0
      b1 = a3 / a4
      b2 = a2 / a4
      b3 = a1 / a4
      b4 = a0 / a4
      d3 = one
      d2 =  - b2
      d1 = b1 * b3 - four * b4
      d0 = b4 * ( four * b2 - b1 * b1 ) - b3 * b3

!  Compute the roots of the auxiliary cubic 

      CALL ROOTS_cubic( d0, d1, d2, d3, tol, nrootsc, rootc1, rootc2, rootc3 )
      IF ( nrootsc > 1 ) rootc1 = rootc3
      x1 = b1 * b1 * quarter - b2 + rootc1
      IF ( x1 < zero ) THEN
        xmd = SQRT( - x1 )
        xnd = quarter * ( two * b3 - b1 * rootc1 ) / xmd
        alpha = half * b1 * b1 - rootc1 - b2
        beta = four * xnd - b1 * xmd
        r = SQRT( alpha * alpha + beta * beta )
        gamma = SQRT( half * ( alpha + r ) )
        IF ( gamma == zero ) THEN
          delta = SQRT( - alpha )
        ELSE
          delta = beta * half / gamma
        END IF
        root1 = half * ( - half * b1 + gamma )
        root2 = half * ( xmd + delta )
        root3 = half * ( - half * b1 - gamma )
        root4 = half * ( xmd - delta )
        GO TO 900
      END IF
      IF ( x1 /= zero ) THEN
        xm = SQRT( x1 )
        xn = quarter * ( b1 * rootc1 - two * b3 ) / xm
      ELSE
        xm = zero
        xn = SQRT( quarter * rootc1 * rootc1 - b4 )
      END IF
      alpha = half * b1 * b1 - rootc1 - b2
      beta = four * xn - b1 * xm
      gamma = alpha + beta
      delta = alpha - beta
      a = - half * b1

!  Compute how many real roots there are

      type_roots = 1
      IF ( gamma >= zero ) THEN
        nroots = nroots + 2
        type_roots = 0
        gamma = SQRT( gamma )
      ELSE
        gamma = SQRT( - gamma )
      END IF
      IF ( delta >= zero ) THEN
        nroots = nroots + 2
        delta = SQRT( delta )
      ELSE
        delta = SQRT( - delta )
      END IF
      type_roots = nroots + type_roots

!  Two real roots

      IF ( type_roots == 3 ) THEN
        root1 = half * ( a - xm - delta )
        root2 = half * ( a - xm + delta )
        root3 = half * ( a + xm )
        root4 = half * gamma
        GO TO 900
      ELSE IF ( type_roots /= 4 ) THEN
        IF ( type_roots == 2 ) THEN
          root1 = half * ( a + xm - gamma )
          root2 = half * ( a + xm + gamma )
        ELSE

!  No real roots

          root1 = half * ( a + xm )
          root2 = half * gamma
        END IF
        root3 = half * ( a - xm ) * half
        root4 = half * delta
        GO TO 900
      END IF

!  Four real roots

      b = half * ( a + xm + gamma )
      d = half * ( a - xm + delta )
      c = half * ( a - xm - delta )
      a = half * ( a + xm - gamma )

!  Sort the roots

      root1 = MIN( a, b, c, d )
      root4 = MAX( a, b, c, d )

      IF ( a == root1 ) THEN
        root2 = MIN( b, c, d )
      ELSE IF ( b == root1 ) THEN
        root2 = MIN( a, c, d )
      ELSE IF ( c == root1 ) THEN
        root2 = MIN( a, b, d )
      ELSE
        root2 = MIN( a, b, c )
      END IF

      IF ( a == root4 ) THEN
        root3 = MAX( b, c, d )
      ELSE IF ( b == root4 ) THEN
        root3 = MAX( a, c, d )
      ELSE IF ( c == root4 ) THEN
        root3 = MAX( a, b, d )
      ELSE
        root3 = MAX( a, b, c )
      END IF

  900 CONTINUE

!  Perfom a Newton iteration to ensure that the roots are accurate

      IF ( debug ) THEN
        IF ( nroots == 0 ) THEN
          WRITE( out, "( ' no real roots ' )" )
        ELSE IF ( nroots == 2 ) THEN
          WRITE( out, "( ' 2 real roots ' )" )
        ELSE IF ( nroots == 4 ) THEN
          WRITE( out, "( ' 4 real roots ' )" )
        END IF
      END IF
      IF ( nroots == 0 ) RETURN

      p = ( ( ( a4 * root1 + a3 ) * root1 + a2 ) * root1 + a1 ) * root1 + a0
      pprime = ( ( four * a4 * root1 + three * a3 ) * root1 + two * a2 )       &
                 * root1 + a1
      IF ( pprime /= zero ) THEN
        IF ( debug ) WRITE( out, 2000 ) 1, root1, p, - p / pprime 
        root1 = root1 - p / pprime
        p = ( ( ( a4 * root1 + a3 ) * root1 + a2 ) * root1 + a1 ) * root1 + a0
      END IF
      IF ( debug ) WRITE( out, 2010 ) 1, root1, p

      p = ( ( ( a4 * root2 + a3 ) * root2 + a2 ) * root2 + a1 ) * root2 + a0
      pprime = ( ( four * a4 * root2 + three * a3 ) * root2 + two * a2 )       &
                 * root2 + a1
      IF ( pprime /= zero ) THEN
        IF ( debug ) WRITE( out, 2000 ) 2, root2, p, - p / pprime 
        root2 = root2 - p / pprime
        p = ( ( ( a4 * root2 + a3 ) * root2 + a2 ) * root2 + a1 ) * root2 + a0
      END IF
      IF ( debug ) WRITE( out, 2010 ) 2, root2, p

      IF ( nroots == 4 ) THEN
        p = ( ( ( a4 * root3 + a3 ) * root3 + a2 ) * root3 + a1 ) * root3 + a0
        pprime = ( ( four * a4 * root3 + three * a3 ) * root3 + two * a2 )     &
                   * root3 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 3, root3, p, - p / pprime 
          root3 = root3 - p / pprime
          p = ( ( ( a4 * root3 + a3 ) * root3 + a2 ) * root3 + a1 ) * root3 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 3, root3, p

        p = ( ( ( a4 * root4 + a3 ) * root4 + a2 ) * root4 + a1 ) * root4 + a0
        pprime = ( ( four * a4 * root4 + three * a3 ) * root4 + two * a2 )     &
                   * root4 + a1
        IF ( pprime /= zero ) THEN
          IF ( debug ) WRITE( out, 2000 ) 4, root4, p, - p / pprime 
          root4 = root4 - p / pprime
          p = ( ( ( a4 * root4 + a3 ) * root4 + a2 ) * root4 + a1 ) * root4 + a0
        END IF
        IF ( debug ) WRITE( out, 2010 ) 4, root4, p
      END IF

      RETURN

!  Non-executable statements

 2000 FORMAT( ' root ', I1, ': value = ', ES12.4, ' quartic = ', ES12.4,       &
              ' delta = ', ES12.4 )
 2010 FORMAT( ' root ', I1, ': value = ', ES12.4, ' quartic = ', ES12.4 )

!  End of subroutine ROOTS_quartic

      END SUBROUTINE ROOTS_quartic

!  End of module ROOTS

   END MODULE GALAHAD_ROOTS_double

