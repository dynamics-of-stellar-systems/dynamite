C     ( Last modified on Mon Feb 25 15:03:37 CST 2002 )
      SUBROUTINE CSCFG ( N, M, X, LC, C, NNZJ, LCJAC, CJAC,
     *                   INDVAR, INDFUN, GRAD )
      INTEGER            N, M, LC, NNZJ, LCJAC
      LOGICAL            GRAD
      INTEGER            INDVAR( LCJAC ), INDFUN( LCJAC )
CS    REAL               X( N ), C( LC ), CJAC  ( LCJAC )
CD    DOUBLE PRECISION   X( N ), C( LC ), CJAC  ( LCJAC )
C
C  *********************************************************************
C  *                                                                   *
C  *             OBSOLETE TOOL! Replaced by CCFSG                     *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER            IOUT
      COMMON / OUTPUT /  IOUT
C
      WRITE( IOUT, 1000 )
      CALL CCFSG ( N, M, X, LC, C, NNZJ, LCJAC, CJAC,
     *             INDVAR, INDFUN, GRAD )
      RETURN
C
C  Non executable statement
C
 1000 FORMAT( ' ** SUBROUTINE CSCFG: this tool is obsolete! ',
     +        ' Please use CCFSG instead.' )
C
C  end of CSCFG.
C
      END

