C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM GENMA
C
C  Generic package driver (example) for applying package GEN to problems
C  from SIF files.
C
C  Ph. Toint, December 2000 / D. Orban, August 2002
C
      IMPLICIT NONE
      INTEGER           NMAX, MMAX
CTOY  PARAMETER       ( NMAX =   10, MMAX = 10   )
CMED  PARAMETER       ( NMAX =  100, MMAX = 100  )
CBIG  PARAMETER       ( NMAX = 3500, MMAX = 3500 )
CCUS  PARAMETER       ( NMAX =  200, MMAX = 200  )
      INTEGER           INSPEC, INPUT, IOUT, N, M
      PARAMETER       ( INSPEC = 46, INPUT = 47, IOUT = 6 )
      INTEGER           NLIN, NEQ, NBNDS, EXITCODE
CS    REAL              DUMMY
CS    REAL              X(NMAX), BL(NMAX), BU(NMAX)
CS    REAL              V(MMAX), CL(MMAX), CU(MMAX)
      REAL              CPU( 2 ), CALLS( 7 )
CD    DOUBLE PRECISION  DUMMY
CD    DOUBLE PRECISION  X(NMAX), BL(NMAX), BU(NMAX)
CD    DOUBLE PRECISION  V(MMAX), CL(MMAX), CU(MMAX)
      CHARACTER*10      PNAME, VNAMES(NMAX), GNAMES(MMAX)
      LOGICAL           EFIRST, LFIRST, NVFRST
      LOGICAL           EQUATN(MMAX), LINEAR(MMAX)
      LOGICAL           CONSTRAINED
C
C  Open the Spec file for the method (typically called METHOD.SPC)
C
      CALL GENSPC( INSPEC, 'GEN.SPC' )
C
C  Open the relevant problem file.
C
      OPEN ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     .       STATUS = 'OLD' )
      REWIND INPUT
C
C  Get problem dimensions and determine which tools to use
C
      CALL CDIMEN( INPUT, N, M )
      IF( M .EQ. 0 ) THEN
         CONSTRAINED = .FALSE.
      ELSE IF( M .GT. 0 ) THEN
         CONSTRAINED = .TRUE.
      ELSE
         WRITE( 6, '(A)' ) 'Error reading OUTSDIF.d'
         STOP
      ENDIF
C
C  Set up parameters
C
      EFIRST = .TRUE.
      LFIRST = .FALSE.
      NVFRST = .FALSE.
C
C  Set up SIF data from the problem file
C
      IF( CONSTRAINED ) THEN
         CALL CSETUP( INPUT, IOUT, N, M, X, BL, BU, NMAX, EQUATN,
     .        LINEAR, V, CL, CU, MMAX, EFIRST, LFIRST, NVFRST )
      ELSE
         CALL USETUP( INPUT, IOUT, N, X, BL, BU, NMAX )
      ENDIF
C
C  Obtain problem/variables/constraints names.
C
      IF( CONSTRAINED ) THEN
         CALL CNAMES( N, M, PNAME, VNAMES, GNAMES )
      ELSE
         CALL UNAMES( N, PNAME, VNAMES )
      ENDIF
C
C  Obtain info on the problem
C
      NLIN  = 0
      NEQ   = 0
      NBNDS = 0
      IF( CONSTRAINED ) THEN
         CALL GETINFO( N, M, BL, BU, EQUATN, LINEAR, NLIN, NEQ, NBNDS )
      ELSE
         EQUATN( 1 ) = .FALSE.
         LINEAR( 1 ) = .FALSE.
         CALL GETINFO( N, 1, BL, BU, EQUATN, LINEAR, NLIN, NEQ, NBNDS )
      ENDIF
C
C  Call the optimizer.
C
      CALL GEN( DUMMY )
      EXITCODE = 0
C
C  Close the problem file
C
      CLOSE( INPUT  )
C
C  Write the standard statistics (of which some may be irrelevant)
C
C    CALLS( 1 ): number of calls to the objective function
C    CALLS( 2 ): number of calls to the objective gradient
C    CALLS( 3 ): number of calls to the objective Hessian
C    CALLS( 4 ): number of Hessian times vector products
C           --constrained problems only--
C    CALLS( 5 ): number of calls to the constraint functions
C    CALLS( 6 ): number of calls to the constraint gradients
C    CALLS( 7 ): number of calls to the constraint Hessians
C           -----------------------------
C
C    CPU( 1 ) : CPU time (in seconds) for USETUP or CSETUP
C    CPU( 2 ) : CPU time ( in seconds) since the end of USETUP or CSETUP
C
C  Note that each constraint function is counted separately.
C  Evaluating all the constraints thus results in PNC evaluations, where
C  PNC is the number of constraints in the problem.  Note that PNC does not
C  include repetitions for constraints having full ranges.

C  (N, is the dimension of the problem, M is the number of constraints,
C   DUMMY is the final value of the objective function)
C
      IF( CONSTRAINED ) THEN
         CALL CREPRT( CALLS, CPU )
      ELSE
         CALL UREPRT( CALLS, CPU )
      ENDIF
      WRITE ( IOUT, 2000 ) PNAME, N, M, NLIN, NEQ, M-NEQ, NBNDS,
     .     CALLS( 1 ), CALLS( 2 ), CALLS( 3 ), CALLS( 5 ), CALLS( 6 ),
     .     CALLS( 7 )
      WRITE ( IOUT, 2001 ) EXITCODE, DUMMY, CPU( 1 ), CPU( 2 )
C
C  Exit
C
      STOP
C
C  Non-executable statements.
C
C  The following is the complete standard statistics output format: select
C  the items that are relevant to the type of problems solved and adapt the
C  name of the code.
C
C  The only reason for breaking the format in two is for compilers
C  which do not accept more than 19 continuation lines.
C
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used                :  GEN',    /
     *    ,' Variant                  :  name of a variant, if needed',/
     *    ,' Problem                  :  ', A10,    /
     *    ,' # variables              =      ', I10 /
     *    ,' # constraints            =      ', I10 /
     *    ,' # linear constraints     =      ', I10 /
     *    ,' # equality constraints   =      ', I10 /
     *    ,' # inequality constraints =      ', I10 /
     *    ,' # bounds                 =      ', I10 /
     *    ,' # objective functions    =        ', F8.2 /
     *    ,' # objective gradients    =        ', F8.2 /
     *    ,' # objective Hessians     =        ', F8.2 /
     *    ,' # constraints functions  =        ', F8.2 /
     *    ,' # constraints gradients  =        ', F8.2 /
     *    ,' # constraints Hessians   =        ', F8.2 )
 2001 FORMAT(
     *     ' Exit code                =      ', I10 /
     *    ,' Final f                  = ', E15.7 /
     *    ,' Set up time              =      ', 0P, F10.2, ' seconds'/
     *     ' Solve time               =      ', 0P, F10.2, ' seconds'//
     *     66('*') / )
      END

