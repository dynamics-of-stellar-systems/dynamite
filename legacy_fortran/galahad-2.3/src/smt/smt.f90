! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-*  G A L A H A D _ S M T   M O D U L E  *-*-*-*-*-*-*-*-*

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released pre GALAHAD Version 1.0. December 1st 1997
!   update released with GALAHAD Version 2.0. February 16th 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_SMT_double

!  ==========================
!  Sparse matrix derived type
!  ==========================

     IMPLICIT NONE

!---------------------
!   P r e c i s i o n
!---------------------

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

     PRIVATE
     PUBLIC :: SMT_put, SMT_get

     TYPE, PUBLIC :: SMT_type
       INTEGER :: m, n, ne
       CHARACTER, ALLOCATABLE, DIMENSION( : ) :: id
       CHARACTER, ALLOCATABLE, DIMENSION( : ) :: type
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: row
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: col
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: ptr
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: val
     END TYPE

   CONTAINS

!-*-*-*-*-*-*-*-*-   S M T  _ p u t   S U B R O U T I N E  -*-*-*-*-*-*-*-*-

     SUBROUTINE SMT_put( array, string, stat )

!  Dummy arguments

     CHARACTER, ALLOCATABLE, DIMENSION( : ) :: array
     CHARACTER ( len = * ), INTENT( IN ), OPTIONAL :: string
     INTEGER, INTENT( OUT ), OPTIONAL ::  stat

!  Local variables

     INTEGER :: i, l
     LOGICAL :: ok

     l = 0
     IF ( PRESENT( string ) ) l = LEN_TRIM( string )
     IF ( PRESENT( stat ) ) THEN
!      IF ( ALLOCATED( array ) ) DEALLOCATE( array, STAT = stat )
       ALLOCATE( array( l ), STAT = stat )
       ok = stat == 0
     ELSE
!      IF ( ALLOCATED( array ) ) DEALLOCATE( array )
       ALLOCATE( array( l ) )
       ok = .true.
     END IF
     IF ( ok ) THEN
       DO i = 1, l
         array( i ) = string( i : i )
       END do
     END IF

     RETURN

!  End of SMT_put

     END SUBROUTINE SMT_put

!-*-*-*-*-*-*-*-*-*-*-   S M T  _ g e t   F U N C T I O N  -*-*-*-*-*-*-*-*-*-

     FUNCTION SMT_get( array )

!  Dummy arguments

     CHARACTER, DIMENSION( : ) :: array
     CHARACTER( SIZE( array ) ) :: SMT_get

!  Local variables

     INTEGER :: i

     DO i = 1, SIZE( array )
        SMT_get( i : i ) = array( i )
     END DO

     RETURN

!  End of function SMT_get

     END FUNCTION SMT_get

!  End of module SMT

   END MODULE GALAHAD_SMT_double


