C     ( Last modified on Mon Feb 25 15:03:37 CST 2002 )
      SUBROUTINE CSCIFG ( N, ICON, X, CI, NNZGCI, LGCI, GCI,
     *                    INDVAR, GRAD )
      INTEGER             N, ICON, NNZGCI, LGCI
      LOGICAL             GRAD
      INTEGER             INDVAR( LGCI )
CS    REAL                CI
CD    DOUBLE PRECISION    CI
CS    REAL                X( N ), GCI( LGCI )
CD    DOUBLE PRECISION    X( N ), GCI( LGCI )
C
C  *********************************************************************
C  *                                                                   *
C  *             OBSOLETE TOOL! Replaced by CCIFSG                     *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER            IOUT
      COMMON / OUTPUT /  IOUT
C
      WRITE( IOUT, 1000 )
      CALL CCIFSG( N, ICON, X, CI, NNZGCI, LGCI, GCI,
     +             INDVAR, GRAD )
      RETURN
C
C  Non executable statement
C
 1000 FORMAT( ' ** SUBROUTINE CSCIFG: this tool is obsolete! ',
     +        ' Please use CCIFSG instead.' )
C
C  end of CSCIFG.
C
      END

