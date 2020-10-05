C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UDIMEN( INPUT , N )
      INTEGER            INPUT , N
C
C  Compute the basic array dimension for the unconstrained optimization tools.
C
C  Nick Gould, for CGT productions,
C  26th August, 1999.
C
      REWIND INPUT
      READ( INPUT, 1001 ) N
      REWIND INPUT
      RETURN
C
C  Non-executable statements.
C
 1001 FORMAT( I8 )
C
C  End of UDIMEN.
C
      END
