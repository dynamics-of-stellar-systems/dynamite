C     ( Last modified on 10 Sepc 2004 at 16:35:38 )
C  Correction: 10/Sep/2004: undeclared integer variables declared
      SUBROUTINE CDIMEN( INPUT , N , M )
      INTEGER            INPUT , N , M
C
C  Compute the basic array dimensions for the constrained optimization tools.
C
C  Nick Gould, for CGT productions,
C  26th August, 1999.
C
      INTEGER            IALGOR, I , J , NG    , NELNUM, IDUMMY
      INTEGER            NG1   , NEL1  , NSLACK, NOBJGR, IEND
      INTEGER            IARRAY( 10 )
      CHARACTER * 8      PNAME
C
C  Input the problem dimensions.
C
      REWIND INPUT
      READ( INPUT, 1001 ) N , NG, NELNUM
C
C  Input the problem type.
C
      READ( INPUT, 1000 ) IALGOR, PNAME
      IF ( IALGOR .LT. 2 ) THEN
         M = 0
         GO TO 100
      END IF
C
C  Set useful integer values.
C
      NG1    = NG     + 1
      NEL1   = NELNUM + 1
C
C  Print out problem data. input the number of variables, groups,
C  elements and the identity of the objective function group.
C
      IF ( IALGOR .EQ. 2 ) READ( INPUT, 1002 ) NSLACK, NOBJGR
C
C  Input the starting addresses of the elements in each group,
C  of the parameters used for each group and
C  of the nonzeros of the linear element in each group.
C
      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NG1 )
      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NG1 )
      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NG1 )
C
C  Input the starting addresses of the variables and parameters
C  in each element.
C
      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NEL1 )
      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NEL1 )
C
C  Input the group type of each group
C
      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NG )
C
C  Count the number of constraint groups
C
      M = 0
      DO 20 I = 1, NG, 10
         IEND = MIN( I + 9, NG )
         READ( INPUT, 1010 ) ( IARRAY( J - I + 1 ), J = I, IEND )
         DO 10 J = I, IEND
           IF ( IARRAY( J - I + 1 ) .NE. 1 ) M = M + 1
   10    CONTINUE
   20 CONTINUE

  100 CONTINUE
      REWIND INPUT
      RETURN
C
C  Non-executable statements.
C
 1000 FORMAT( I2, A8 )
 1001 FORMAT( 3I8 )
 1002 FORMAT( 2I8 )
 1010 FORMAT( ( 10I8 ) )
C
C  End of CDIMEN.
C
      END
