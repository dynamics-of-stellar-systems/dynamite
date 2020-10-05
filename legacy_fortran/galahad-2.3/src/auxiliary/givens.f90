! THIS VERSION: GALAHAD 2.1 - 29/10/2007 AT 09:30 GMT.

!-*-*-*-*-*-  G A L A H A D    N O R M S   M O D U L E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould
!
!  History -
!   originally released pre GALAHAD Version 2.0. May 22nd 2004

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

    MODULE GALAHAD_NORMS_double

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: ONE_NORM, TWO_NORM, INFINITY_NORM

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp

!---------------------------------
!   I n t e r f a c e  B l o c k s
!---------------------------------

     INTERFACE NRM2

       FUNCTION SNRM2( n, X, incx )
       REAL :: SNRM2
       INTEGER, INTENT( IN ) :: n, incx
       REAL, INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
       END FUNCTION SNRM2

       FUNCTION DNRM2( n, X, incx )
       DOUBLE PRECISION :: DNRM2
       INTEGER, INTENT( IN ) :: n, incx
       DOUBLE PRECISION, INTENT( IN ), DIMENSION( incx * ( n - 1 ) + 1 ) :: X
       END FUNCTION DNRM2

     END INTERFACE 

    CONTAINS

!-*-*-*-*-  G A L A H A D   O N E  _ N O R M   F U N C T I O N   -*-*-*-*-

       FUNCTION ONE_NORM( X )

!  Compute the l_1 norm of the vector X

!  Dummy arguments

       REAL ( KIND = wp ) :: ONE_NORM
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( : ) :: X

!  Local variable

       INTEGER :: n
       n = SIZE( X )

       IF ( n > 0 ) THEN
         ONE_NORM = SUM( ABS( X ) )
       ELSE
         ONE_NORM = zero
       END IF
       RETURN

!  End of function ONE_NORM

       END FUNCTION ONE_NORM

!-*-*-*-*-  G A L A H A D   T W O  _ N O R M   F U N C T I O N   -*-*-*-*-

       FUNCTION TWO_NORM( X )

!  Compute the l_2 norm of the vector X

!  Dummy arguments

       REAL ( KIND = wp ) :: TWO_NORM
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( : ) :: X

!  Local variable

       INTEGER :: n
       n = SIZE( X )

       IF ( n > 0 ) THEN
         TWO_NORM = NRM2( n, X, 1 )
       ELSE
         TWO_NORM = zero
       END IF
       RETURN

!  End of function TWO_NORM

       END FUNCTION TWO_NORM

!-*-*-*-  G A L A H A D   I N F I N I T Y  _ N O R M   F U N C T I O N   -*-*-*-

       FUNCTION INFINITY_NORM( X )

!  Compute the l_infinity norm of the vector X

!  Dummy arguments

       REAL ( KIND = wp ) :: INFINITY_NORM
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( : ) :: X

!  Local variable

       INTEGER :: n
       n = SIZE( X )

       IF ( n > 0 ) THEN
         INFINITY_NORM = MAXVAL( ABS( X ) )
       ELSE
         INFINITY_NORM = zero
       END IF
       RETURN

!  End of function INFINITY_NORM

       END FUNCTION INFINITY_NORM

!  End of module GALAHAD_NORMS_double

    END MODULE GALAHAD_NORMS_double

