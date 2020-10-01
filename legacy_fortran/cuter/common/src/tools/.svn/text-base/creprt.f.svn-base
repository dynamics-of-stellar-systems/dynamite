C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CREPRT( CALLS, TIME )
C
C  This routine returns the value of the various counters maintained by the
C  CUTEr tools to the user.  The counters are:

C    CALLS( 1 ): number of calls to the objective function
C    CALLS( 2 ): number of calls to the objective gradient
C    CALLS( 3 ): number of calls to the objective Hessian
C    CALLS( 4 ): number of Hessian times vector products
C    CALLS( 5 ): number of calls to the constraint functions
C    CALLS( 6 ): number of calls to the constraint gradients
C    CALLS( 7 ): number of calls to the constraint Hessians
C
C    TIME( 1 ): CPU time (in seconds) for CSETUP
C    TIME( 2 ): CPU time (in seconds) since the end of CSETUP
C
C  Note that each constraint function is counted separately.
C  Evaluating all the constraints thus results in PNC evaluations, where
C  PNC is the number of constraints in the problem.  Note that PNC does not
C  include repetitions for constraints having full ranges.  Also note that
C
C  N. Gould, D. Orban & Ph. Toint for CUTEr, 2001.
C
      IMPLICIT NONE
C
C  output arguments
C
      REAL               CALLS( 7 )
      REAL               TIME( 2 )
C
C  variables from the PRFCTS common block
C
      INTEGER            NC2OF , NC2OG , NC2OH, NC2CF, NC2CG, NC2CH
      INTEGER            NHVPR , PNC
      REAL               SUTIME, STTIME
      COMMON / PRFCTS /  NC2OF , NC2OG , NC2OH, NC2CF, NC2CG, NC2CH,
     *                   NHVPR , PNC   , SUTIME, STTIME
      SAVE             / PRFCTS /
C
C  local variables
C
      REAL               CPUTIM, DUM
      EXTERNAL           CPUTIM
C
      TIME( 2 )  = CPUTIM( DUM ) - STTIME
      TIME( 1 )  = SUTIME

      CALLS( 1 ) = NC2OF
      CALLS( 2 ) = NC2OG
      CALLS( 3 ) = NC2OH
      CALLS( 4 ) = NHVPR
      IF( PNC .GT. 0 ) THEN
        CALLS( 5 ) = NC2CF / PNC
        CALLS( 6 ) = NC2CG / PNC
        CALLS( 7 ) = NC2CH / PNC
      ELSE
        CALLS( 5 ) = NC2CF
        CALLS( 6 ) = NC2CG
        CALLS( 7 ) = NC2CH
      ENDIF

      RETURN
      END
