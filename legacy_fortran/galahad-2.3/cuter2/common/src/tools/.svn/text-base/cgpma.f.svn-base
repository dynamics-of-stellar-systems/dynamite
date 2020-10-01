C     ( Last modified on 30 Jul 2004 at 16:31:38 )
      PROGRAM          CGPMA
C
C  CG+ test driver for problems derived from SIF files.
C
C  Nick Gould, for CGT Productions.
C  July 2004
C
      INTEGER          NMAX  , IOUT  , N , INPUT , MAXIT 
      INTEGER          LP    , MP, I , METHOD, ITER  , NFUN, IREST
      INTEGER          IFLAG , INSPEC, IPRINT( 2 )
CS    REAL             F, EPS, GNORM , BIGINF, ZERO, TLEV
CD    DOUBLE PRECISION F, EPS, GNORM , BIGINF, ZERO, TLEV
      LOGICAL          BOUNDS, FINISH
CBIG  PARAMETER      ( NMAX = 100000 )
CMED  PARAMETER      ( NMAX =  10000 )
CTOY  PARAMETER      ( NMAX =   1000 )
CCUS  PARAMETER      ( NMAX =  50000 )
      PARAMETER      ( IOUT  = 6 )
      PARAMETER      ( INPUT = 55, INSPEC = 56 )
CS    PARAMETER      ( BIGINF = 9.0E+19, ZERO = 0.0E0 )
CD    PARAMETER      ( BIGINF = 9.0D+19, ZERO = 0.0D0 )
CS    REAL             X     ( NMAX ), G     ( NMAX ),
CD    DOUBLE PRECISION X     ( NMAX ), G     ( NMAX ),
     *                 D     (NMAX ), GOLD  ( NMAX ), W     ( NMAX )
      CHARACTER * 10   XNAMES( NMAX ), PNAME, SPCDAT
      INTRINSIC        DABS   , MAX
      COMMON / CGDD /  MP, LP
      COMMON / RUNINF / ITER, NFUN
      REAL             CPU( 2 ), CALLS( 4 )
C     
C  Open the Spec file for the method.
C
      SPCDAT = 'CGPLUS.SPC'
      OPEN ( INSPEC, FILE = SPCDAT, FORM = 'FORMATTED',
     *      STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
C     IPRINT(1): specifies the frequency of output
C     IPRINT(2): specifies the amount of output
C     METHOD   : method used (1=Fletcher-Reeves,2=Polak-Ribiere,3=P-R+)
C     IREST    : no restart (0) or restart every n iterations (1)
C     MAXIT    : maximum number of iterations
C     EPS      : the required norm of the gradient
C
      READ ( INSPEC, 1000 ) IPRINT( 1 ), IPRINT( 2 ), METHOD, IREST,
     *                      MAXIT, EPS
C
C  Close input file.
C
      CLOSE ( INSPEC )
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
      CALL USETUP( INPUT, IOUT, N, X, W, GOLD, NMAX )
C
C  Obtain variable names.
C
      CALL UNAMES( N, PNAME, XNAMES )
C
C  Set up algorithmic input data.
C
      BOUNDS  = .FALSE.
      DO 10 I = 1, N
         IF ( W( I ) .GT. - BIGINF .OR. GOLD( I ) .LT. BIGINF )
     *      BOUNDS = .TRUE.
   10 CONTINUE
      IF ( BOUNDS ) WRITE( IOUT, 2030 )
      LP     = IOUT
      MP     = IOUT
      ITER  = - 1
      IFLAG  = 0
      FINISH = .FALSE.
C
C  Optimization loop
C
   20 CONTINUE
         ITER = ITER + 1
C
C  Evaluate the function and gradient
C
         CALL UOFG( N, X, F, G, .TRUE. )
C
C  Call the optimizer.
C
   30    CONTINUE
         CALL CGFAM( N, X, F, G, D, GOLD, IPRINT, EPS, W,
     *               IFLAG, IREST, METHOD, FINISH )
C
C  Test for termination
C
        IF ( IFLAG .LE. 0 .OR. ITER .GT. MAXIT ) GO TO 50
        IF ( IFLAG .EQ. 1 ) GO TO 20
        IF ( IFLAG .EQ. 2 ) THEN
C
C Termination Test.  The user may replace it by some other test. However, 
C the parameter 'FINISH' must be set to 'TRUE' when the test is satisfied.
C
           TLEV = EPS * ( 1.0D+0 + DABS( F ) )
           DO 40 I = 1, N
              IF( DABS( G( I ) ) .GT. TLEV ) GO TO 30
  40       CONTINUE
           FINISH = .TRUE.
           GO TO 30
        ENDIF

   50 CONTINUE
C
C  Terminal exit.
C
      CALL UREPRT( CALLS, CPU )
      GNORM    = ZERO
      DO 110 I  = 1, N
         GNORM = MAX( GNORM, DABS( G( I ) ) )
  110 CONTINUE
      WRITE ( IOUT, 2010 ) F, GNORM
C      DO 120 I = 1, N
C         WRITE( IOUT, 2020 ) XNAMES( I ), X( I ), G( I )
C  120 CONTINUE
      WRITE ( IOUT, 2000 ) PNAME, N, INT( CALLS(1) ), INT( CALLS(2) ),
     *                     IFLAG, F, CPU(1), CPU(2) 
      CLOSE( INPUT  )
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( 5( I10, / ), D10.3 )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  CG+',     /
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
 2030 FORMAT(  /, ' ** Warning from CGPMA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
 2040 FORMAT( /, ' ** ERROR from CGPMA. The declared array ', A6, 
     *           ' is too small to hold the problem data.', /, 
     *           ' Increase ', A6, ' in CGPMA to be at least ', I6, 
     *           ' and recompile. Stopping ' )
      END
