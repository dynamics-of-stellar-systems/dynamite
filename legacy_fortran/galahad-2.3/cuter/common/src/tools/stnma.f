C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM          STNMA
C
C  STENMIN test driver for problems derived from SIF files.
C
C  Ph. Toint, January 1996, for CGT Productions.
C
      INTEGER          NMAX, N, NZ, LIRN, LICN, ILIM, METHOD
      INTEGER          GRDFLG, HSNFLG, NDIGIT, MSG, LWRK, LIWRK
      INTEGER          TERMCD, INFORM, I, IOUT, INPUT, INSPEC
      PARAMETER      ( INPUT =    55, IOUT =      6, INSPEC =   46 )
CTOY  PARAMETER      ( NMAX  =   100, LIRN =    500, LICN =    500 )
CMED  PARAMETER      ( NMAX  =  1000, LIRN =   5000, LICN =   5000 )
CBIG  PARAMETER      ( NMAX  = 20000, LIRN = 500000, LICN = 500000 )
CCUS  PARAMETER      ( NMAX  = 10000, LIRN =  10000, LICN =  10000 )
      PARAMETER      ( LIWRK = 2 * LIRN + 12 * NMAX + 2 )
      PARAMETER      ( LWRK = 7 * NMAX )
      INTEGER          IRN ( LIRN  ), ICN ( LICN ), IWRK( LIWRK )
CD    DOUBLE PRECISION FSCALE, GRADTL, STEPTL, FPLS, STEPMX
CS    REAL             FSCALE, GRADTL, STEPTL, FPLS, STEPMX
CS    REAL             GNORM, BIGINF, ZERO
CD    DOUBLE PRECISION GNORM, BIGINF, ZERO
CS    PARAMETER      ( BIGINF = 9.0E+19, ZERO = 0.0E0 )
CD    PARAMETER      ( BIGINF = 9.0D+19, ZERO = 0.0D0 )
CD    DOUBLE PRECISION X   ( NMAX ),  TYPX( NMAX ), XPLS( NMAX )
CD    DOUBLE PRECISION GPLS( NMAX ),  HESS( LICN ), WRK ( LWRK )
CS    REAL             X   ( NMAX ),  TYPX( NMAX ), XPLS( NMAX )
CS    REAL             GPLS( NMAX ),  HESS( LICN ), WRK ( LWRK )
CS    REAL             BU    ( NMAX ), BL  ( NMAX )
CD    DOUBLE PRECISION BU    ( NMAX ), BL  ( NMAX )
      REAL             CPU( 2 ), CALLS( 4 )
      LOGICAL          BOUNDS
      CHARACTER * 10   XNAMES( NMAX ), PNAME
C
      INTRINSIC        ABS   , MAX
      EXTERNAL         UFN, UGR, USH
C     
C  Open the Spec file for the method.
C
      OPEN ( INSPEC, FILE = 'STNMIN.SPC', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
C     ILIM  :  the maximum number of iterations
C     GRADTL:  the relative gradient stopping tolerance
C     GRDFLG:  the gradient availability and checking flag
C     HSNFLG:  the Hessian availability and checking flag
C     FSCALE:  the typical value of the objective function
C     TYPX  :  the typical value of the problem's variables
C     METHOD:  the method used (0 = Newton, 1 = tensor )
C     NDIGIT:  the number of accurate digits in function values
C     MSG   :  the output specifier
C
      READ ( INSPEC, 1000 ) ILIM, GRADTL, GRDFLG, HSNFLG, FSCALE,
     *                      TYPX( 1 ), METHOD, NDIGIT, MSG
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
      BOUNDS  = .FALSE.
      DO 10 I = 1, N
         TYPX(I) = TYPX( 1 )
         IF ( BL( I ) .GT. - BIGINF .OR. BU( I ) .LT. BIGINF )
     *      BOUNDS = .TRUE.
   10 CONTINUE
      IF ( BOUNDS ) WRITE( IOUT, 2030 )
      INFORM = 0
      STEPTL = GRADTL * GRADTL
      STEPMX = ZERO
C
C  Call the optimizer.
C
CD    CALL STUMCD( N, X, NZ, IRN, LIRN, ICN, LICN, UFN, UGR, 
CS    CALL STUMCS( N, X, NZ, IRN, LIRN, ICN, LICN, UFN, UGR,
     *             USH, TYPX, FSCALE, GRADTL, STEPTL, ILIM, STEPMX,
     *             IOUT, METHOD, GRDFLG, HSNFLG, NDIGIT, MSG, XPLS, 
     *             FPLS, GPLS, HESS, WRK, LWRK, IWRK, LIWRK, TERMCD, 
     *             BL, INFORM )
      CALL UREPRT( CALLS, CPU )
      GNORM    = ZERO
      DO 20 I  = 1, N
         GNORM = MAX( GNORM, ABS( GPLS( I ) ) )
   20 CONTINUE
      WRITE ( IOUT, 2040 ) GNORM
      WRITE ( IOUT, 2010 )
      DO 30 I = 1, N
         WRITE( IOUT, 2020 ) XNAMES(I), XPLS( I ), GPLS( I )
   30 CONTINUE
      WRITE ( IOUT, 2000 ) PNAME, N, CALLS(1), CALLS(2), CALLS(3), 
     *                     TERMCD, FPLS, CPU(1), CPU(2) 
      CLOSE( INPUT  )
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( I10, /, D10.3, 2(/,I10), 2(/,D10.3), 3(/,I10) )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  STENMIN',  /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # objective functions   =        ', F8.2 /
     *    ,' # objective gradients   =        ', F8.2 / 
     *    ,' # objective Hessians    =        ', F8.2 /
     *    ,' Exit code               =      ', I10 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *    ,' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 2010 FORMAT( /, '                 X         G ' )
 2020 FORMAT(  A10, 1P, 2D12.4 )
 2030 FORMAT(  /, ' ** Warning from STNMA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
 2040 FORMAT( ' Final gradient norm = ', D12.4 / )
      END

