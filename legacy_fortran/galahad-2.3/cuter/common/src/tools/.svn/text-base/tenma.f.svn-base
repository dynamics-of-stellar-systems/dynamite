C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM          TENMA
C
C  TENSOR test driver for problems derived from SIF files.
C
C  Ali Bouaricha, CERFACS (April, 1993), updated by Ph. Toint (March 2001)
C
      INTEGER          NMAX  , IOUT  , N     , METHOD, INPUT , ITNNO
      INTEGER          MSG   , NDIGIT, ILIM  , IAGFLG, IAHFLG, INSPEC
CS    REAL             GRADTL, STEPTL, FSCALE, FPLS  , STEPMX, TYPX
CD    DOUBLE PRECISION GRADTL, STEPTL, FSCALE, FPLS  , STEPMX, TYPX
CTOY  PARAMETER      ( NMAX =  10 )
CMED  PARAMETER      ( NMAX = 100 )
CBIG  PARAMETER      ( NMAX = 500 )
CCUS  PARAMETER      ( NMAX = 200 )
      PARAMETER      ( INPUT = 55, IOUT = 6, INSPEC = 56 )
CS    REAL             X( NMAX ), TYPSIZ( NMAX ), XPLS( NMAX )
CS    REAL             H( NMAX, NMAX ), WRK( NMAX, 8 ), GPLS( NMAX )
CD    DOUBLE PRECISION X( NMAX ), TYPSIZ( NMAX ), XPLS( NMAX )
CD    DOUBLE PRECISION H( NMAX, NMAX ), WRK( NMAX, 8 ), GPLS( NMAX )
      CHARACTER*10     XNAMES( NMAX ), PNAME
      INTEGER          IWRK( NMAX )  , I
      LOGICAL          BOUNDS
CS    REAL             GNORM, BIGINF, ZERO
CD    DOUBLE PRECISION GNORM, BIGINF, ZERO
      INTRINSIC        ABS  , MAX
CS    PARAMETER      ( BIGINF = 9.0E+19, ZERO = 0.0E0 )
CD    PARAMETER      ( BIGINF = 9.0D+19, ZERO = 0.0D0 )
      REAL             CPU( 2 ), CALLS( 4 )
      EXTERNAL         FCN, GRD, HSN
C     
C  Open the Spec file for the method.
C
      OPEN ( INSPEC, FILE = 'TENMIN.SPC', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
      READ ( INSPEC, 1000 ) ILIM, GRADTL, IAGFLG, IAHFLG, FSCALE,
     *                      TYPX, METHOD, MSG
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
      CALL USETUP( INPUT, IOUT, N, X, XPLS, GPLS, NMAX )
C
C  Obtain variable names.
C
      CALL UNAMES( N, PNAME, XNAMES )
C
C  Set up algorithmic input data.
C
CS    NDIGIT = 7
CD    NDIGIT = 15
CS    STEPTL = 0.00001E0
CD    STEPTL = 0.00001D0
      STEPMX = BIGINF
      BOUNDS  = .FALSE.
      DO 10 I = 1, N
         TYPSIZ( I ) = TYPX
         IF ( XPLS( I ) .GT. - BIGINF .OR. GPLS( I ) .LT. BIGINF )
     *      BOUNDS = .TRUE.
   10 CONTINUE
      IF ( BOUNDS ) WRITE( IOUT, 2030 )
C
C  Call the optimizer.
C
      CALL TENSOR( NMAX  , N, X  , FCN   , GRD   , HSN   , TYPSIZ,
     *             FSCALE, GRADTL, STEPTL, ILIM  , STEPMX, IOUT  , 
     *             METHOD, IAGFLG, IAHFLG, NDIGIT, MSG   , XPLS  , 
     *             FPLS  , GPLS  , H     , ITNNO , WRK   , IWRK  )
      CALL UREPRT( CALLS, CPU )
C
      GNORM    = ZERO
      DO 20 I  = 1, N
         GNORM = MAX( GNORM, ABS( GPLS( I ) ) )
   20 CONTINUE
      WRITE ( IOUT, 2010 )
      DO 30 I = 1, N
         WRITE( IOUT, 2020 ) XNAMES(I), XPLS( I ), GPLS( I )
   30 CONTINUE
      WRITE ( IOUT, 2001 ) GNORM
      WRITE ( IOUT, 2000 ) PNAME, N, ( CALLS( I ), I = 1, 3 ),
     *                     ITNNO, FPLS, CPU( 1 ), CPU( 2 )
      CLOSE( INPUT  )
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( I10, /, D10.3, 2(/,I10), 2(/,D10.3), 2(/,I10) )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  TENMIN',   /
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
 2001 FORMAT( / ' FINAL NORMG = ', 1P, D12.4 )
 2010 FORMAT( /, '                 X         G ' )
 2020 FORMAT(  A10, 1P, 2D12.4 )
 2030 FORMAT(  /, ' ** Warning from TENMA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
      END
C
C
C
      SUBROUTINE FCN   ( N, X, F )
      INTEGER            N
CS    REAL               F
CD    DOUBLE PRECISION   F
CS    REAL               X( N )
CD    DOUBLE PRECISION   X( N )
C
C  Interface for TENSOR (Chow, Schnabel, Eskow).
C
C  Ali Bouaricha, CERFACS.
C  April, 1993.
C
      EXTERNAL         UFN
      CALL UFN   ( N, X, F )
      RETURN
C
C  end of FCN.
C
      END
C
C
C
      SUBROUTINE GRD   ( N, X, G )
      INTEGER            N
CS    REAL               X( N ), G( N )
CD    DOUBLE PRECISION   X( N ), G( N )
C
C  Interface for TENSOR (Chow, Schnabel, Eskow).
C
C  Ali Bouaricha, CERFACS.
C  April, 1993.
C
      EXTERNAL         UGR
      CALL UGR   ( N, X, G )
C
C  end of GRD.
C
      END
C
C
C
      SUBROUTINE HSN   ( NR, N, X, H )
      INTEGER            N, NR
CS    REAL               X( N ), H( NR, N )
CD    DOUBLE PRECISION   X( N ), H( NR, N )
C
C  Interface for TENSOR (Chow, Schnabel, Eskow).
C
C  Ali Bouaricha, CERFACS.
C  April, 1993.
C
      EXTERNAL         UDH
      CALL UDH   ( N, X, NR, H )
      RETURN
C
C  end of HSN.
C
      END
