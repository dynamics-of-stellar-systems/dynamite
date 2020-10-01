C     ( Last modified on 12 Sepc 2004 at 09:45:12 )
C  Correction: 12/Sep/2004: undeclared integer variable declared
      PROGRAM          UNCMA
C
C  UNCMIN test driver for problems derived from SIF files.
C
C  Nick Gould, for CGT Productions, October 1991.
C  Ph. Toint, December 2000.
C
      INTEGER          NMAX  , IOUT  , N     , METHOD, INPUT , ITRMCD
      INTEGER          IEXP  , MSG   , NDIGIT, ILIM  , IAGFLG, IAHFLG
      INTEGER          INSPEC
CS    REAL             DLT   , GRADTL, STP   , STEPTL, FSCALE, FPLS
CD    DOUBLE PRECISION DLT   , GRADTL, STP   , STEPTL, FSCALE, FPLS
CTOY  PARAMETER      ( NMAX =  10 )
CMED  PARAMETER      ( NMAX = 100 )
CBIG  PARAMETER      ( NMAX = 3500 )
CCUS  PARAMETER      ( NMAX = 200 )
      PARAMETER      ( INPUT = 55, IOUT = 6 , INSPEC = 46)
CS    REAL             X     ( NMAX ), TYPSIZ( NMAX ), XPLS( NMAX )
CS    REAL             A     ( NMAX, NMAX )  , WRK   ( NMAX, 8 )
CS    REAL             GPLS  ( NMAX ), BU    ( NMAX ), BL  ( NMAX )
CD    DOUBLE PRECISION X     ( NMAX ), TYPSIZ( NMAX ), XPLS( NMAX )
CD    DOUBLE PRECISION A     ( NMAX, NMAX )  , WRK   ( NMAX, 8 )
CD    DOUBLE PRECISION GPLS  ( NMAX ), BU    ( NMAX ), BL  ( NMAX )
      CHARACTER*10     XNAMES( NMAX ), PNAME
      EXTERNAL         FCN   , D1FN , D2FN
      INTEGER          I
      LOGICAL          BOUNDS
CS    REAL             XMAX, GNORM, ONE, ZERO, BIGINF, TYPX
CD    DOUBLE PRECISION XMAX, GNORM, ONE, ZERO, BIGINF, TYPX
      INTRINSIC        ABS   , MAX
CS    PARAMETER      ( ONE = 1.0E+0, ZERO = 0.0E+0, BIGINF = 9.0E+19 )
CD    PARAMETER      ( ONE = 1.0D+0, ZERO = 0.0D+0, BIGINF = 9.0D+19 )
      REAL             CPU( 2 ), CALLS( 4 )
C
C  Open the spec file
C
      OPEN( INSPEC, FILE = 'UNCMIN.SPC', FORM = 'FORMATTED',
     *      STATUS = 'OLD' )
C
C  Read input Spec data.
C
C TYPX        <--  TYPICAL SIZE FOR EACH COMPONENT OF X
C FSCALE      <--  ESTIMATE OF SCALE OF MINIMIZATION FUNCTION
C METHOD      <--  ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C                  ( 1 = LINESEARCH, 2 = DOUBLE DOGLEG, 3 = HEBDEN-MORE )
C IEXP        <--  =0 IF MINIMIZATION FUNCTION NOT EXPENSIVE TO EVALUATE
C MSG         <--  MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
C NDIGIT      <--  NUMBER OF GOOD DIGITS IN MINIMIZATION FUNCTION
C ITNLIM      <--  MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C IAGFLG      <--  =0 IF ANALYTIC GRADIENT NOT SUPPLIED
C IAHFLG      <--  =0 IF ANALYTIC HESSIAN NOT SUPPLIED
C DLT         <--  INITIAL TRUST REGION RADIUS
C GRADTL      <--  TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE ENOUGH
C                  TO ZERO TO TERMINATE ALGORITHM
C STEPTL      <--  TOLERANCE AT WHICH SUCCESSIVE ITERATES CONSIDERED
C                  CLOSE ENOUGH TO TERMINATE ALGORITHM
C
      READ ( INSPEC, 1000 ) TYPX, FSCALE, METHOD, IEXP, MSG, NDIGIT,
     +                      ILIM, IAGFLG, IAHFLG, DLT, GRADTL, STEPTL
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
      CALL USETUP( INPUT, IOUT, N, X, BL, BU, NMAX )
C
C  Obtain variable names.
C
      CALL UNAMES( N, PNAME, XNAMES )
C
C  Set up algorithmic input data.
C
      XMAX    = ZERO
      BOUNDS  = .FALSE.
      DO 10 I = 1, N
         TYPSIZ( I ) = TYPX
         XMAX        = MAX( XMAX, ABS( X( I ) ) )
         IF ( BL( I ) .GT. - BIGINF .OR. BU( I ) .LT. BIGINF )
     *      BOUNDS = .TRUE.
   10 CONTINUE
      IF ( BOUNDS ) WRITE( IOUT, 2030 )
      STP    = 1000.0 * MAX( ONE, XMAX )
C
C  Call the optimizer.
C
      CALL OPTIF9( NMAX  , N , X , FCN   , D1FN  , D2FN  , TYPSIZ,
     *             FSCALE, METHOD, IEXP  , MSG   , NDIGIT, ILIM  ,
     *             IAGFLG, IAHFLG, IOUT  , DLT   , GRADTL, STP   ,
     *             STEPTL, XPLS  , FPLS  , GPLS  , ITRMCD, A, WRK )
      CALL UREPRT ( CALLS, CPU )
      GNORM    = ZERO
      DO 20 I  = 1, N
         GNORM = MAX( GNORM, ABS( GPLS( I ) ) )
   20 CONTINUE
      WRITE ( IOUT, 2010 )
      DO 30 I = 1, N
         WRITE( IOUT, 2020 ) XNAMES(I), XPLS( I ), GPLS( I )
   30 CONTINUE
      WRITE ( IOUT, 2000 ) PNAME, N, ( CALLS(I), I = 1, 3 ),
     *                     ITRMCD, FPLS, CPU( 1 ), CPU( 2 )
      CLOSE( INPUT  )
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( 2( D10.3, /), 7( I10, /), 2( D10.3, / ), D10.3 )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  UNCMIN',   /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # objective functions   =        ', F8.2 /
     *    ,' # objective gradients   =        ', F8.2 / 
     *    ,' # objective Hessians    =        ', F8.2 /
     *     ' Exit code               =      ', I10 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 2010 FORMAT( /, '                 X         G ' )
 2020 FORMAT(  A10, 1P, 2D12.4 )
 2030 FORMAT(  /, ' ** Warning from UNCMA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
      END
C
      SUBROUTINE FCN   ( N, X, F )
      INTEGER            N
CS    REAL               F, X( N )
CD    DOUBLE PRECISION   F, X( N )
C
C  Interface for UNCMIN (Schnabel, Koontz and Weiss,
C  ACM Trans. Math. Software, 1982).
C
      EXTERNAL         UFN
      CALL UFN   ( N, X, F )
      RETURN
      END
C
      SUBROUTINE D1FN  ( N, X, G )
      INTEGER            N
CS    REAL               X( N ), G( N )
CD    DOUBLE PRECISION   X( N ), G( N )
C
      EXTERNAL         UGR
      CALL UGR   ( N, X, G )
      END
C
      SUBROUTINE D2FN  ( NR, N, X, H )
      INTEGER            N, NR
CS    REAL               X( N ), H( NR, N )
CD    DOUBLE PRECISION   X( N ), H( NR, N )
C
C  Interface for UNCMIN (Schnabel, Koontz and Weiss,
C  ACM Trans. Math. Software, 1982).
C
      EXTERNAL         UDH
      CALL UDH   ( N, X, NR, H )
      RETURN
      END
