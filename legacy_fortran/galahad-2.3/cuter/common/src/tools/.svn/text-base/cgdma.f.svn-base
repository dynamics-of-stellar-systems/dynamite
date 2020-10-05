C     ( Last modified on 30 Jul 2004 at 16:31:38 )
      PROGRAM          CGDMA
C
C  CG_DESCENT test driver for problems derived from SIF files.
C
C  Nick Gould, for CGT Productions.
C  July 2004
C
      INTEGER          NMAX  , IOUT  , N , M , INPUT , MSAVE, IOUTCP
      INTEGER          LP    , MP, I , METHOD, ITER  , NF, NG
      INTEGER          STAT  , ICALL , INSPEC, IPRINT( 2 )
CS    REAL             F, TOL, GNORM , BIGINF, ZERO,
CD    DOUBLE PRECISION F, TOL, GNORM , BIGINF, ZERO,
     *                 STPMIN, STPMAX, TLEV
      LOGICAL          BOUNDS
CS    REAL             DELTA, SIGMA, EPSILON, THETA, GAMMA, 
CD    DOUBLE PRECISION DELTA, SIGMA, EPSILON, THETA, GAMMA, 
     *                 RHO, ETA, PSI0, PSI1, PSI2, QUADCU, RSTRTF, 
     *                 MAXITF, FEPS
      INTEGER          NXPAND, NSECNT
      LOGICAL          QUADST, PRNTLV, PRNTFI, STRULE, ERULE, AWOLFE, 
     *                 STEP
CBIG  PARAMETER      ( NMAX = 100000 )
CMED  PARAMETER      ( NMAX =  10000 )
CTOY  PARAMETER      ( NMAX =   1000 )
CCUS  PARAMETER      ( NMAX =  50000 )
      PARAMETER      ( IOUT  = 6 )
      PARAMETER      ( INPUT = 55, INSPEC = 56, IOUTCP = 57 )
CS    PARAMETER      ( BIGINF = 9.0E+19, ZERO = 0.0E0 )
CD    PARAMETER      ( BIGINF = 9.0D+19, ZERO = 0.0D0 )
CS    REAL             X     ( NMAX ), G     ( NMAX ),
CD    DOUBLE PRECISION X     ( NMAX ), G     ( NMAX ),
     *                 D     ( NMAX ), XTEMP  ( NMAX ), GTEMP  ( NMAX )
      CHARACTER * 10   XNAMES( NMAX ), PNAME, CGPARM
      CHARACTER * 15   SPCDAT
      REAL             CPU( 2 ), CALLS( 4 )
      EXTERNAL         EVALF, EVALG
C     
C  Open the Spec file for the method.
C
      SPCDAT = 'CG_DESCENT.SPC'
      OPEN ( INSPEC, FILE = SPCDAT, FORM = 'FORMATTED',
     *      STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
C     TOL         : the required norm of the gradient
C     DELTA       : (Wolfe line search parameter)
C     SIGMA       : (Wolfe line search parameter)
C     EPSILON     : (approximate Wolfe threshold factor)
C     THETA       : (update)
C     GAMMA       : (required decay factor in interval)
C     RHO         : (growth factor in bracket)
C     ETA         : (lower bound for cg's beta_k)
C     PSI0        : (factor used in very initial starting guess)
C     PSI1        : (factor previous step multiplied by in QuadStep)
C     PSI2        : (factor previous step is multipled by for startup)
C     QUADCUTOFF  : (lower bound on rel change in f before QuadStep)
C     RESTART_FAC : (restart cg in restart_fac*n iterations)
C     MAXIT_FAC   : (terminate in maxit_fac*n iterations)
C     FEPS        : (stop when value change <= feps*|f|)
C     NEXPAND     : (number of grow/shrink allowed in bracket)
C     NSECANT     : (number of secant steps allowed in line search)
C     QUADSTEP    : (use initial quad interpolation in line search)
C     PRINTLEVEL  : F (no print) T (intermediate results)
C     PRINTFINAL  : F (no print) T (print error messages, final error)
C     STOPRULE    : F (|grad|_infty <= tol) T (... <= tol*(1+|f|))
C     ERULE       : F (eps_k = eps|f|) T (eps_k = eps)
C     AWOLFE      : F (Wolfe) T (+approx Wolfe) 2 (epsilon pert)
C     STEP        : F (no initial line search guess) T (guess in gnorm)
C
      READ ( INSPEC, 1000 ) TOL, DELTA, SIGMA, EPSILON, THETA, GAMMA, 
     *                      RHO, ETA, PSI0, PSI1, PSI2, QUADCU, RSTRTF, 
     *                      MAXITF, FEPS, NXPAND, NSECNT, QUADST, 
     *                      PRNTLV, PRNTFI, STRULE, ERULE, AWOLFE, STEP
C
C  Close input file.
C
      CLOSE ( INSPEC )
C
C  Create the required data input file cg.parm
C
      CGPARM = 'cg.parm'
      OPEN ( IOUTCP, FILE = CGPARM, FORM = 'FORMATTED',
     *       STATUS = 'UNKNOWN' )
      REWIND IOUTCP
      WRITE( IOUTCP, 2100 ) DELTA, SIGMA, EPSILON, THETA, GAMMA, 
     *                      RHO, ETA, PSI0, PSI1, PSI1
      WRITE( IOUTCP, 2110 ) QUADCU, RSTRTF, MAXITF, FEPS, 
     *                      NXPAND, NSECNT, QUADST, PRNTLV
      WRITE( IOUTCP, 2120 ) PRNTFI, STRULE, ERULE, AWOLFE, STEP
      CLOSE ( IOUTCP )
C
C  Open the relevant file.
C
      OPEN ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
C
C  Check to see if there is sufficient room
C
      CALL UDIMEN( INPUT, N )
      IF ( N .GT. NMAX ) THEN
        WRITE( IOUT, 2040 ) 'X     ', 'NMAX  ', N
        STOP
      END IF
C
C  Set up SIF data.
C
      CALL USETUP( INPUT, IOUT, N, X, XTEMP, GTEMP, NMAX )
C
C  Obtain variable names.
C
      CALL UNAMES( N, PNAME, XNAMES )
C
C  Set up algorithmic input data.
C
      BOUNDS  = .FALSE.
      DO 10 I = 1, N
         IF ( XTEMP( I ) .GT. - BIGINF .OR. GTEMP( I ) .LT. BIGINF )
     *      BOUNDS = .TRUE.
   10 CONTINUE
      IF ( BOUNDS ) WRITE( IOUT, 2030 )
      LP     = IOUT
      MP     = IOUT
C
C  Call the optimizer.
C
      CALL CG_DESCENT( TOL, X, N, EVALF, EVALG, STAT, GNORM, F,
     &                 ITER, NF, NG, D, G, XTEMP, GTEMP )
C
C  Terminal exit.
C
      CALL UREPRT( CALLS, CPU )
      WRITE ( IOUT, 2010 ) F, GNORM
      DO 120 I = 1, N
         WRITE( IOUT, 2020 ) XNAMES( I ), X( I ), G( I )
  120 CONTINUE
      WRITE ( IOUT, 2000 ) PNAME, N, INT( CALLS(1) ), INT( CALLS(2) ),
     *                     STAT, F, CPU(1), CPU(2) 
      CLOSE( INPUT  )
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( 15( D10.3, / ), 2( I10, / ), 6( L10, / ), L10 )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  CG_DESCENT',     /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # objective functions   =      ', I10 /
     *    ,' # objective gradients   =      ', I10 / 
     *     ' Exit code               =      ', I10 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 2010 FORMAT( ' Final objective function value  = ', 1P, D12.4, 
     *        /, ' Final norm of gradient          = ', 1P, D12.4,
     *        //, '                 X         G ' )
 2020 FORMAT(  1X, A10, 1P, 2D12.4 )
 2030 FORMAT(  /, ' ** Warning from CGDMA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
 2040 FORMAT( /, ' ** ERROR from CGDMA. The declared array ', A6, 
     *           ' is too small to hold the problem data.', /, 
     *           ' Increase ', A6, ' in CGDMA to be at least ', I6, 
     *           ' and recompile. Stopping ' )
 2100 FORMAT( 1P, D10.3, 4X, 
     * 'delta        Wolfe line search parameter',
     * /, 1P, D10.3, 4X, 
     * 'sigma        Wolfe line search parameter',
     * /, 1P, D10.3, 4X, 
     * 'epsilon      approximate Wolfe threshold factor',
     * /, 1P, D10.3, 4X, 
     * 'theta        update',
     * /, 1P, D10.3, 4X, 
     * 'gamma        required decay factor in interval',
     * /, 1P, D10.3, 4X, 
     * 'rho          growth factor in bracket',
     * /, 1P, D10.3, 4X, 
     * 'eta          lower bound for cg''s beta_k',
     * /, 1P, D10.3, 4X, 
     * 'psi0         factor used in very initial starting guess',
     * /, 1P, D10.3, 4X, 
     * 'psi1         factor previous step multiplied by in QuadStep',
     * /, 1P, D10.3, 4X, 
     * 'psi2         factor previous step is multipled by for startup')
 2110 FORMAT( 1P, D10.3, 4X, 
     * 'QuadCutOff   lower bound on rel change in f before QuadStep',
     * /, 1P, D10.3, 4X, 
     * 'restart_fac  restart cg in restart_fac*n iterations',
     * /, 1P, D10.3, 4X, 
     * 'maxit_fac    terminate in maxit_fac*n iterations',
     * /, 1P, D10.3, 4X, 
     * 'feps         stop when value change <= feps*|f|',
     * /, I10, 4X, 
     * 'nexpand      number of grow/shrink allowed in bracket',
     * /, I10, 4X, 
     * 'nsecant      number of secant steps allowed in line search',
     * /, L10, 4X, 
     * 'QuadStep     use initial quad interpolation in line search',
     * /, L10, 4X, 
     * 'PrintLevel   F (no print) T (intermediate results)')
 2120 FORMAT(  L10, 4X, 
     * 'PrintFinal   F (no print) T (print error messages, ',
     * 'final error)',
     * /, L10, 4X, 
     * 'StopRule     F (|grad|_infty <= tol) T (... <= tol*(1+|f|))',
     * /, L10, 4X, 
     * 'ERule        F (eps_k = eps|f|) T (eps_k = eps)',
     * /, L10, 4X, 
     * 'AWolfe       F (Wolfe) T (+approx Wolfe) 2 (epsilon pert)',
     * /, L10, 4X, 
     * 'Step         F (no initial line search guess) T (guess in ',
     * 'gnorm')

      END

      SUBROUTINE EVALF( F, X, N )

C  Evaluate the objective function

      INTEGER N
CS    REAL             F
CD    DOUBLE PRECISION F
CS    REAL             X( N )
CD    DOUBLE PRECISION X( N )
      EXTERNAL UFN

      CALL UFN( N, X, F )

      RETURN
      END
      
      SUBROUTINE EVALG( G, X, N )

C  Evaluate the gradiuent of the objective function

      INTEGER N
CS    REAL             X( N ), G( N )
CD    DOUBLE PRECISION X( N ), G( N )
      EXTERNAL UGR

      CALL UGR( N, X, G )

      RETURN
      END
