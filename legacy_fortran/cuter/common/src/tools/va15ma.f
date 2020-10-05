C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM          VA15MA
C
C  VA15 test driver for problems derived from SIF files.
C
C  Nick Gould, for CGT Productions.
C  October 1993.
C  Modified by Ph. Toint.
C
      INTEGER          NMAX  , IOUT  , N , M , INPUT , MSAVE 
      INTEGER          LP    , MP    , INFO  , LW, I , MAXIT
      INTEGER          IFLAG , ICALL , INSPEC, IPRINT( 2 )
CS    REAL             F, EPS, GTOL  , GNORM , BIGINF, ZERO
CD    DOUBLE PRECISION F, EPS, GTOL  , GNORM , BIGINF, ZERO
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
CS    EXTERNAL         VA15C
CD    EXTERNAL         VA15CD
CS    COMMON         / VA15D  / GTOL , MP, LP, INFO
CD    COMMON         / VA15DD / GTOL , MP, LP, INFO
      REAL             CPU( 2 ), CALLS( 4 )
C     
C  Open the Spec file for the method.
C
      SPCDAT = 'VA15.SPC'
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
C     EPS      : the stopping tolerance
C
      READ ( INSPEC, 1000 ) M, IPRINT( 1 ), IPRINT( 2 ), MAXIT, 
     *                      EPS
C
C  Close input file.
C
      CLOSE ( INSPEC )
C
C  Open the relevant file.
C
      OPEN ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND INPUT
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
   20 CONTINUE
C
C  Evaluate the function and gradient.
C
      CALL UOFG( N, X, F, G, .TRUE. )
C
C  Call the optimizer.
C
CS    CALL VA15A ( N, M, X, F, G, DIAGCO, DIAG, IPRINT, EPS, W, IFLAG )
CD    CALL VA15AD( N, M, X, F, G, DIAGCO, DIAG, IPRINT, EPS, W, IFLAG )
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
      WRITE ( IOUT, 2001 ) GNORM
      WRITE ( IOUT, 2010 )
      DO 40 I = 1, N
         WRITE( IOUT, 2020 ) XNAMES( I ), X( I ), G( I )
   40 CONTINUE
      WRITE ( IOUT, 2000 ) PNAME, N, CALLS(1), CALLS(2),
     *                     IFLAG, F, CPU(1), CPU(2) 
      CLOSE( INPUT  )
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( 4( I10, / ), D10.3 )
 2001 FORMAT( /, ' Final norm of G  = ', 1P, D12.4 )
 2010 FORMAT( /, '                 X         G ' )
 2020 FORMAT(  1X, A10, 1P, 2D12.4 )
 2030 FORMAT(  /, ' ** Warning from VA15MA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  VA15',     /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # objective functions   =        ', F8.2 /
     *    ,' # objective gradients   =        ', F8.2 / 
     *     ' Exit code               =      ', I10 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
      END

