! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.
      PROGRAM GALAHAD_RAND_SPEC
      USE GALAHAD_RAND_double
      IMPLICIT NONE
      INTEGER :: random_integer, seed
      REAL ( kind = KIND( 1.0D+0 ) ) :: random_real
!  Get the current generator word
      CALL RAND_GET_SEED( seed )
      WRITE( 6, "( ' generator word = ', I10 )" ) seed
!  Generate a random real in [-1, 1]
      CALL RAND_RANDOM_REAL( .FALSE., random_real )
      WRITE( 6, "( ' random real = ', F10.2 )" ) random_real
!  Generate another random real
      CALL RAND_RANDOM_REAL( .FALSE., random_real )
      WRITE( 6, "( ' second random real = ', F10.2 )" ) random_real
!  Restore the generator word
      CALL RAND_SET_SEED( seed )
!  Generate a random integer in [1, 100]
      CALL RAND_RANDOM_INTEGER( 100, random_integer ) 
      WRITE( 6, "( ' random integer = ', I3 )" ) random_integer
!  Generate another random integer
      CALL RAND_RANDOM_INTEGER( 100, random_integer )
      WRITE( 6, "( ' second random integer = ', I3 )" ) random_integer
      END PROGRAM GALAHAD_RAND_SPEC
