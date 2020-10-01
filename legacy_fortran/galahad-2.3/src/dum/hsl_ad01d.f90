! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.
! Updated 25/06/2002: additional warning information added

!-*-*-*-*-*  L A N C E L O T  -B-  DUMMY AD01_FORWARD  M O D U L E S *-*-*-*

!  Nick Gould, for GALAHAD productions
!  Copyright reserved
!  June 28th 1996

MODULE HSL_AD01_FORWARD_DOUBLE

      IMPLICIT NONE
      PRIVATE
      INTEGER, PARAMETER :: WP = KIND(1D0)
      PUBLIC :: AD01_INITIALIZE
    
!  Dummy HSL_AD01_FORWARD_DOUBLE module

      TYPE, PUBLIC :: AD01_REAL
        PRIVATE
        INTEGER P, CASE
      END TYPE AD01_REAL

CONTAINS

      SUBROUTINE AD01_INITIALIZE(DEGREE,A,VALUE,FULL_THRESHOLD)
        INTEGER, INTENT (IN) :: DEGREE
!       TYPE (AD01_REAL), INTENT (OUT) :: A
        TYPE (AD01_REAL) :: A
        INTEGER, OPTIONAL, INTENT (IN) :: FULL_THRESHOLD
        REAL (WP), INTENT (IN) :: VALUE
    
!  Dummy subroutine available with LANCELOT

        WRITE ( 6, 2000 )
        STOP

!  Non-executable statements

 2000    FORMAT( /, ' We regret that the solution options that you have ', /, &
                    ' chosen are not all freely available with LANCELOT.', //,&
                    ' If you have HSL (formerly the Harwell Subroutine',      &
                    ' Library), this ', /,                                    &
                    ' option may be enabled by replacing the dummy ', /,      &
                    ' module HSL_AD01_FORWARD_DOUBLE with its H...', /,       &
                    ' namesake and dependencies. See', /,                     &
                    '   $GALAHAD/src/makedefs/packages for details.', //,     &
                    ' *** EXECUTION TERMINATING *** ', / )

      END SUBROUTINE AD01_INITIALIZE

END MODULE HSL_AD01_FORWARD_DOUBLE

!  THIS VERSION: 28/06/1996 AT 09:00:00 AM

!-*-*-*-*-*  L A N C E L O T  -B-  DUMMY AD01_BACKWARD  M O D U L E S *-*-*-*

!  Nick Gould, for GALAHAD productions
!  Copyright reserved
!  June 28th 1996

MODULE HSL_AD01_BACKWARD_DOUBLE

      IMPLICIT NONE
      PRIVATE
      INTEGER, PARAMETER :: WP = KIND(1D0)
      PUBLIC :: AD01_INITIALIZE
    
!  Dummy HSL_AD01_BACKWARD_DOUBLE module

      TYPE, PUBLIC :: AD01_REAL
        PRIVATE
        INTEGER P, CASE
      END TYPE AD01_REAL

CONTAINS

      SUBROUTINE AD01_INITIALIZE(DEGREE,A,VALUE,FULL_THRESHOLD)
        INTEGER, INTENT (IN) :: DEGREE
!       TYPE (AD01_REAL), INTENT (OUT) :: A
        TYPE (AD01_REAL) :: A
        INTEGER, OPTIONAL, INTENT (IN) :: FULL_THRESHOLD
        REAL (WP), INTENT (IN) :: VALUE
    
!  Dummy subroutine available with LANCELOT

        WRITE ( 6, 2000 )
        STOP

!  Non-executable statements

 2000    FORMAT( /, ' We regret that the solution options that you have ', /, &
                    ' chosen are not all freely available with LANCELOT.', //,&
                    ' If you have HSL (formerly the Harwell Subroutine',      &
                    ' Library), this ', /,                                    &
                    ' option may be enabled by replacing the dummy ', /,      &
                    ' module HSL_AD01_BACKWARD_DOUBLE with its HSL ',/,       &
                    ' namesake and dependencies. See', /,                     &
                    '   $GALAHAD/src/makedefs/packages for details.', //,     &
                    ' *** EXECUTION TERMINATING *** ', / )

      END SUBROUTINE AD01_INITIALIZE

END MODULE HSL_AD01_BACKWARD_DOUBLE
