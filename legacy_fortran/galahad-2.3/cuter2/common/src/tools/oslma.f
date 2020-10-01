C     ( Last modified on 10 Sepc 2004 at 17:05:38 )
C  Correction: 10/Sep/2004: undeclared integer variables declared
C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM OSLMA
C
C  OSL test driver for problems derived from SIF files.
C
      INTEGER          NMAX  , MMAX  , M , N 
      INTEGER          MCON  , MEQ   , LW    , LIW   
      INTEGER          NNZMAX, NNZJ  , NNZH  , NNZTH
      INTEGER          I , J , IC, IR, ICNH  , INDV ,
     *                 INDF  , IRNH  , MGEQ  , NJAC
      CHARACTER * 10   PNAME
CTOY  PARAMETER      ( NMAX = 100, MMAX = 100 )
CMED  PARAMETER      ( NMAX = 200, MMAX = 200 )
CBIG  PARAMETER      ( NMAX = 500, MMAX = 500 )
CCUS  PARAMETER      ( NMAX = 400, MMAX = 300 )
      PARAMETER      ( NNZMAX = 20 * NMAX                  )
      PARAMETER      ( LW   = 11000 + 10 * NMAX + 10 * MMAX + NNZMAX )
      PARAMETER      ( LIW  = NMAX + 1 + 4 * NNZMAX                  )
      INTEGER          IW    ( LIW )
CS    REAL             X     ( NMAX  ), BL    ( NMAX   ), BU( NMAX ),
CD    DOUBLE PRECISION X     ( NMAX  ), BL    ( NMAX   ), BU( NMAX ),
     *                 G     ( NMAX  ), C     ( MMAX   ),
     *                 CL    ( MMAX  ), CU    ( MMAX   ), W ( LW   )
CS    REAL             OBJF, GET(34)
CD    DOUBLE PRECISION OBJF, GET(34)
      LOGICAL          EQUATN( MMAX  ), LINEAR( MMAX   )
      INTEGER          INPUT, IOUT,  RTCOD
      LOGICAL          DEBUG
      PARAMETER      ( INPUT = 55, IOUT = 6 )
CS    REAL             ZERO
CD    DOUBLE PRECISION ZERO  
CS    PARAMETER      ( ZERO  = 0.0E+0 )
CD    PARAMETER      ( ZERO  = 0.0D+0 )
      REAL             CPU( 2 ), CALLS( 7 )
      DEBUG = .FALSE.
C
C  Open the relevant file.
C
      OPEN ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND INPUT
C
C  Set up the data structures necessary to hold the problem functions.
C
      CALL OSLSE ( INPUT , IOUT  , N , M , MGEQ  , MEQ    , MCON ,
     *              X , BL, BU    , NMAX  , EQUATN, LINEAR,
     *              W , CL, CU    , MMAX  )

C
C  Set up algorithmic input data.
C
      
C
C  Start of main iteration loop.
C
C
C  Set X to zero to determine the constant terms and the linear terms for
C  the problem functions.
C
      OBJF = ZERO
      DO 30 I = 1, N
         X( I ) = ZERO
         W( I ) = ZERO
   30 CONTINUE
C
C  Evaluate the constant terms of the objective and the linear coefficients.
C
      CALL COFG( N ,  X, OBJF, G(1), .TRUE.)
C
C  Partition the integer sparse work vector IW.
C
         INDV  = 0
         INDF  = INDV + NNZMAX
         IRNH  = INDF + NNZMAX
         ICNH  = IRNH + NMAX
C
C  Evaluate the constant and linear terms of the constraint functions
C  in a sparse format.
C
      CALL CSCFG ( N, M, X, MMAX, C, NNZJ, NJAC, G( N + 1 ),
     *             IW( INDV + 1 ), IW( INDF + 1 ), .TRUE. )
C
C  Set up row bounds
C
      DO 40 I = 1, M
         CL(I) = CL(I) - C(I)
         CU(I) = CU(I) - C(I)
   40 CONTINUE
      IF ( DEBUG ) THEN
         WRITE( 6, 2005 ) ( I, X( I ), BL( I ), BU( I ), I = 1, N )
         IF ( MCON .GT. 0 ) WRITE( 6, 2060 ) ( I, C( I ), CL( I ),
     *        CU( I ), EQUATN( I ), LINEAR( I ), I = 1, MCON )
      END IF
C
C  Evaluate the Hessian matrix using sparse format.
C
      CALL CSH ( N , M , X, MMAX, W, NNZH  , NMAX, G( N + NNZJ + 1 ) , 
     *             IW( IRNH + 1 ), IW( ICNH + 1 ) )
C
C  DEBUG output
C
      IF ( DEBUG ) THEN
         WRITE( 6, 2010 ) OBJF
         IF ( M .GT. 0 ) WRITE( 6, 2070 ) ( I, C( I ), I = 1, M )
         WRITE( 6, 2080 )
         WRITE( 6, 2015 ) ( G( I ), I = 1, N )
         DO 11 J = 1, NNZJ
            IC = IW( INDV + J )
            IR =  IW( INDF + J )
            WRITE( 6, 2020 )  IR, IC, G( N + J )
   11    CONTINUE
         DO 12 J = 1, NNZH
            IC = IW( ICNH + J )
            IR =  IW( IRNH + J )
            WRITE( 6, 2030 )  IR, IC, G( N + NNZJ + J )
   12    CONTINUE
      END IF
C
C  Fill-in lower half of Hessian
C
         NNZTH = NNZH
         DO 14 J = 1, NNZTH
            IC = IW( ICNH + J )
            IR =  IW( IRNH + J )
            IF ( IR .NE. IC) THEN
               NNZH = NNZH + 1
               IW( ICNH + NNZH ) = IR
               IW( IRNH + NNZH ) = IC
               G( N + NNZJ +NNZH ) = G( N + NNZJ + J )
            END IF 
   14    CONTINUE

C
C  Describe workspace and allow room for one matrix.
C
      CALL EKKDSCA(RTCOD,W, LW,1)
      IF ( RTCOD .GT. 0 ) CALL CHKRT('EKKDSCA',RTCOD)
      WRITE( 6, 2000 ) N, M, MEQ, MCON, RTCOD
C
C   Describe the model. Minimum of 5 blocks are needed for QP.
C
      CALL EKKDSCM(RTCOD,W,1,5)
      IF (RTCOD.GT.0) CALL CHKRT('EKKDSCM',RTCOD)
C
C   Pass linear model with matrix stored by indices.
C
      CALL EKKLMDL(RTCOD, W, 1, M, N, NNZJ, G(1), CL(1), CU(1),
     *             BL(1), BU(1), IW(INDF+1), IW(INDV+1), G(N+1) )
      IF (RTCOD.GT.0) CALL CHKRT('EKKLMDL',RTCOD)
C
C   Pass quadratic matrix stored by indices.
C
      CALL EKKQMDL(RTCOD, W, 1, NNZH, IW(IRNH+1), IW(ICNH+1), 
     *             G( N + NNZJ + 1 ) )
      IF (RTCOD.GT.0) CALL CHKRT('EKKQMDL',RTCOD)
C
C   Solve the QP using the primal algorithm.
C
      CALL EKKQSLV(RTCOD,W,1,1)
      IF (RTCOD.GT.0) CALL CHKRT('EKKQSLV',RTCOD)
      
C
C   Print the solution.
C
      CALL EKKPRTS(RTCOD,W)
      IF (RTCOD.GT.0) CALL CHKRT('EKKPRTS',RTCOD)

C
C   Compensate for any nonzero constant in the objective
C
      IF (OBJF .NE. ZERO) THEN
        CALL EKKRGET(RTCOD,W,GET,18)
        IF (RTCOD.GT.0) CALL CHKRT('EKKRGET',RTCOD)
        OBJF = GET(18) + OBJF
        WRITE(6,2085) OBJF
      END IF

C
C   Get and print statistics
C
      CALL CREPRT( CALLS, CPU )
C     CALL CNAMES( N, M, PNAME, XNAME, CNAME )
      WRITE ( IOUT, 1500 ) PNAME, N, CALLS(1), CALLS(2), CALLS(3), 
     *                     RTCOD, OBJF, CPU(1), CPU(2)

      STOP
C
C   Non executable statements
C
 1500 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  OSL',  /
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
 2000 FORMAT( /, ' N = ', I6, ' M = ', I6, ' MEQ = ', I6,
     *           ' MCON = ', I6, ' RTCOD = ', I6 )
 2010 FORMAT( /, ' Objective function constant is ', 1P, D22.14 )
 2015 FORMAT( /, ' Objective function linear terms are', /, ( 1P, 
     *        D12.4 ) )
 2020 FORMAT( /, ' Constraint function',I6,' linear terms are with',
     *           ' respect to variable ', 
     *        /,  I6, 1P, D12.4  )
 2030 FORMAT( /, ' Objective function Hessian terms with respect to',
     *           ' variables ', 
     *        /,  I6, 1P, I6, 1P, D12.4  )
 2005 FORMAT( /, ' After projection, the starting point:',
     *        /, '     I      X          BL          BU', /,
     *        ( I6, 1P, 3D12.4 ) )
 2060 FORMAT( /, ' the constraints:', /,
     *         '     I  MULTIPLIER     CL          CU      EQUALITY? ',
     *       '  LINEAR? ', /, ( I6, 1P, 3D12.4, 5X, L1, 10X, L1 ) )
 2070 FORMAT( /, ' the constraint constant values are:',
     *        /, '     I      C ', /, ( I6, 1P, D12.4 ) )
 2080 FORMAT( /, ' Linear terms ' )
 2085 FORMAT( /, ' True objective function value ',
     *           '(OSL ignores constants) ', 1P, D22.14 )
C
C  End of OSLMA.
C
      END
C
C***********************************************************************
C   This subroutine prints the character string RTNAME and the return
C   code RTCOD and stops if RTCOD is large enough to indicate that an
C   OSL error or severe error has occured.
C***********************************************************************
C
      SUBROUTINE CHKRT(RTNAME,RTCOD)
      CHARACTER*7 RTNAME
      INTEGER*4   RTCOD
C
      WRITE(6,9000) RTNAME,RTCOD
      IF (RTCOD.GE.200) STOP 16
      RETURN
9000  FORMAT (1X,'********** ',A7,' return code of ',I4,' **********')
      END
      SUBROUTINE OSLSE( INPUT , IOUT  , N , M , MGEQ  , MEQ    , 
     *                   MCON  , X , BL, BU    , NMAX  , EQUATN, 
     *                   LINEAR, V , CL, CU    , MMAX  )
      INTEGER            INPUT , IOUT  , N , M , MGEQ  , MEQ    , MCON
      INTEGER            NMAX  , MMAX
CS    REAL               X     ( NMAX ), BL    ( NMAX ), BU    ( NMAX )
CS    REAL               V     ( MMAX ), CL    ( MMAX ), CU    ( MMAX )
CD    DOUBLE PRECISION   X     ( NMAX ), BL    ( NMAX ), BU    ( NMAX )
CD    DOUBLE PRECISION   V     ( MMAX ), CL    ( MMAX ), CU    ( MMAX )
      LOGICAL            EQUATN( MMAX ), LINEAR( MMAX )
C
C  Set up the input data for the the OSL minimizer.
C
C
      INTEGER            I
C
C  Set up the data structures necessary to hold the problem functions.
C
      CALL CSETUP ( INPUT , IOUT  , N     , MCON  , X , BL , BU   ,
     *              NMAX  , EQUATN, LINEAR, V , CL, CU     , MMAX,
     *              .TRUE., .FALSE., .FALSE. )
C
C  Count the number of general equality constraints.
C
      M    = MCON
      MGEQ = 0
      DO 20 I = 1, M
         IF ( EQUATN( I ) ) MGEQ = MGEQ + 1
   20 CONTINUE
C     IF ( M .GT. 0 ) WRITE( 6, 2010 ) ( I, V( I ), CL( I ), CU( I ),
C    *                             EQUATN( I ), LINEAR( I ), I = 1, M )
C
      MEQ = MGEQ
      IF ( M .GT. MMAX ) THEN
         IF ( IOUT .GT. 0 )
     *      WRITE( IOUT, 2000 ) 'EQUATN', 'MMAX  ', M - MMAX
         STOP
      END IF
      RETURN
C
C  Non-executable statements.
C
 2000 FORMAT( /, ' ** Program CSETUP: array length ', A6, ' too small.',
     *        /, ' -- Miminimization abandoned.',
     *        /, ' -- Increase the parameter ', A6, ' by at least ', I8,
     *           ' and restart.'  )
C2010 FORMAT( /, ' the constraints:', /,
C    *         '     I  MULTIPLIER     CL          CU      EQUALITY? ',
C    *       '  LINEAR? ', /, ( I6, 1P, 3D12.4, 5X, L1, 10X, L1 ) )
C
C  End of OSLSE.
C
      END
