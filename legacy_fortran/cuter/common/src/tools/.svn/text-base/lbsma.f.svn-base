C     ( Last modified on 30 Jul 2004 at 16:31:38 )
      PROGRAM          LBFMA
C
C  LBFGS test driver for problems derived from SIF files.
C
C  Nick Gould and Ph. Toint, for CGT Productions.
C  July 2004
C
      INTEGER          NMAX  , IOUT  , N , M , INPUT , MSAVE 
      INTEGER          LP    , MP    , LW, I , MAXIT
      INTEGER          IFLAG , ICALL , INSPEC, IPRINT( 2 )
CS    REAL             F, EPS, XTOL, GTOL  , GNORM , BIGINF, SMACHR,
CD    DOUBLE PRECISION F, EPS, XTOL, GTOL  , GNORM , BIGINF, DMACHR,
     *                 ZERO  , STPMIN, STPMAX
      LOGICAL          DIAGCO, BOUNDS
CBIG  PARAMETER      ( NMAX = 100000 )
CMED  PARAMETER      ( NMAX =  10000 )
CTOY  PARAMETER      ( NMAX =   1000 )
CCUS  PARAMETER      ( NMAX =  50000 )
      PARAMETER      ( IOUT  = 6, MSAVE = 7 )
      PARAMETER      ( LW    = NMAX * ( 2 * MSAVE + 1 ) + 2 * MSAVE )
      PARAMETER      ( INPUT = 55, INSPEC = 56 )
CS    PARAMETER      ( BIGINF = 9.0E+19, ZERO = 0.0E0 )
CD    PARAMETER      ( BIGINF = 9.0D+19, ZERO = 0.0D0 )
CS    REAL             X     ( NMAX ), G     ( NMAX )
CS    REAL             DIAG  ( NMAX ), W     ( LW   )
CD    DOUBLE PRECISION X     ( NMAX ), G     ( NMAX )
CD    DOUBLE PRECISION DIAG  ( NMAX ), W     ( LW   )
      CHARACTER * 10   XNAMES( NMAX ), PNAME, SPCDAT
      INTRINSIC        ABS   , MAX
      EXTERNAL         LB2
      COMMON         / LB3 / MP, LP, GTOL, STPMIN, STPMAX
      REAL             CPU( 2 ), CALLS( 4 )
CS    EXTERNAL         SMACHR
CD    EXTERNAL         DMACHR
C     
C  Open the Spec file for the method.
C
      SPCDAT = 'LBFGS.SPC'
      OPEN ( INSPEC, FILE = SPCDAT, FORM = 'FORMATTED',
     *      STATUS = 'OLD' )
      REWIND INSPEC
C                                                                               
C
C  Read input Spec data.
C
C     M        : the number of iterations in the memory
C     IPRINT(1): specifies the frequency of output
C     IPRINT(2): specifies the amount of output
C     MAXIT    : the maximum number of iterations,
C     EPS      : the required norm of the gradient
C
      READ ( INSPEC, 1000 ) M, IPRINT( 1 ), IPRINT( 2 ), MAXIT, 
     *                      EPS
      M = MIN( M, MSAVE )
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
      CALL USETUP( INPUT, IOUT, N, X, W, W( NMAX + 1 ), NMAX )
C
C  Obtain variable names.
C
      CALL UNAMES( N, PNAME, XNAMES )
C
C  Set up algorithmic input data.
C
      BOUNDS  = .FALSE.
      DO 10 I = 1, N
         IF ( W( I ) .GT. - BIGINF .OR. W( NMAX + I ) .LT. BIGINF )
     *      BOUNDS = .TRUE.
   10 CONTINUE
      IF ( BOUNDS ) WRITE( IOUT, 2030 )
      LP     = IOUT
      MP     = IOUT
      ICALL  = 0
      IFLAG  = 0
      DIAGCO = .FALSE.
CS    XTOL = SMACHR( 4 )
CD    XTOL = DMACHR( 4 )
   20 CONTINUE
C
C  Evaluate the function and gradient.
C
      CALL UOFG( N, X, F, G, .TRUE. )
C
C  Call the optimizer.
C
CS    CALL LBFGS( N, M, X, F, G, DIAGCO, DIAG, IPRINT, EPS, XTOL,
CD    CALL LBFGS( N, M, X, F, G, DIAGCO, DIAG, IPRINT, EPS, XTOL,
     *            W, IFLAG )
C
C  Check exit conditions and prepare for re-entry.
C
      IF ( IFLAG. GT. 0 ) THEN
         ICALL = ICALL + 1
         IF ( ICALL .LE. MAXIT ) GO TO 20
      END IF
C
C  Terminal exit.
C
      CALL UREPRT( CALLS, CPU )
      GNORM    = ZERO
      DO 30 I  = 1, N
         GNORM = MAX( GNORM, ABS( G( I ) ) )
   30 CONTINUE
      WRITE ( IOUT, 2010 ) F, GNORM
      DO 40 I = 1, N
         WRITE( IOUT, 2020 ) XNAMES( I ), X( I ), G( I )
   40 CONTINUE
      WRITE ( IOUT, 2000 ) PNAME, N, INT( CALLS(1) ), INT( CALLS(2) ),
     *                     IFLAG, F, CPU(1), CPU(2) 
      CLOSE( INPUT  )
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( 4( I10, / ), D10.3 )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  L-BFGS',    /
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
 2030 FORMAT(  /, ' ** Warning from LBFMA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
 2040 FORMAT( /, ' ** ERROR from LBSMA. The declared array ', A6, 
     *           ' is too small to hold the problem data.', /, 
     *           ' Increase ', A6, ' in LBSMA to be at least ', I6, 
     *           ' and recompile. Stopping ' )
      END

