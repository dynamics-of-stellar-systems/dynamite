C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM VF13MA
C
C  VF13 test driver for problems derived from SIF files.
C
C  Nick Gould, for CGT Productions.
C  July 1991.
C
      INTEGER          NMAX  , MMAX  , INF   , M , N , MAXFUN
      INTEGER          MCON  , LCN   , MEQ   , LW    , LIW   , IPRINT
      INTEGER          I , J , MGEQ
CTOY  PARAMETER ( NMAX = 10, MMAX = 10 )
CMED  PARAMETER ( NMAX = 100, MMAX = 100 )
CBIG  PARAMETER ( NMAX = 500, MMAX = 500 )
CCUS  PARAMETER ( NMAX = 200, MMAX = 200 )
      PARAMETER      ( LCN   = MMAX                      )
      PARAMETER      ( LW   = 5 * NMAX * NMAX / 2 + 43 * NMAX / 2 +
     *                        16 * MMAX + 14                         )
      PARAMETER      ( LIW  = NMAX + 1                               )
      INTEGER          IW    ( LIW )
CS    REAL             X     ( NMAX  ), BL    ( NMAX   ), BU( NMAX ),
CD    DOUBLE PRECISION X     ( NMAX  ), BL    ( NMAX   ), BU( NMAX ),
     *                 G     ( NMAX  ), C     ( MMAX   ),
     *                 CL    ( MMAX  ), CU    ( MMAX   ),
     *                 CN    ( LCN    , MMAX           ), W ( LW   )
      LOGICAL          EQUATN( MMAX  ), LINEAR( MMAX   )
      INTEGER          INPUT , IOUT
      LOGICAL          FIRSTG, DEBUG
      CHARACTER * 10   PNAME, VNAME( NMAX ), CNAME( MMAX )
      PARAMETER      ( INPUT = 55, IOUT = 6 )
CS    REAL             F, ACCREQ, ACC
CD    DOUBLE PRECISION F, ACCREQ, ACC
CS    PARAMETER      ( ACCREQ = 1.0E-7 )
CD    PARAMETER      ( ACCREQ = 1.0D-7 )
      REAL             CPU( 2 ), CALLS( 7 )
C     DEBUG = .TRUE.
      DEBUG = .FALSE.
C
C  Open the relevant file.
C
C    statements for OUTSDIF.d with generic UNIX OPEN statements.
      OPEN ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND INPUT
C
C  Set up the data structures necessary to hold the problem functions.
C
      CALL VF13SE ( INPUT , IOUT  , N , M , MGEQ  , MEQ    , MCON ,
     *              X , BL, BU    , NMAX  , EQUATN, LINEAR,
     *              C , CL, CU    , MMAX  )
      IF ( DEBUG ) THEN
         WRITE( 6, 2030 ) ( I, X( I ), BL( I ), BU( I ), I = 1, N )
         IF ( MCON .GT. 0 ) WRITE( 6, 2060 ) ( I, C( I ), CL( I ),
     *        CU( I ), EQUATN( I ), LINEAR( I ), I = 1, MCON )
      END IF
C
C  Set up algorithmic input data.
C
      MAXFUN = 1000
      IPRINT = 0
      INF    = - 1
      ACC    = ACCREQ
      FIRSTG = .TRUE.
C
C  Start of main iteration loop.
C
   10 CONTINUE
C
C  Evaluate the objective function and constraints.
C
      CALL VF13FN ( N , MGEQ, MEQ, MCON, X , F, MMAX,
     *              C , BL, BU, CL, CU )
      IF ( DEBUG ) THEN
         WRITE( 6, 2010 ) F
         IF ( M .GT. 0 ) WRITE( 6, 2070 ) ( I, C( I ), I = 1, M )
      END IF
C
C  Evaluate the gradient of the objective and constraint functions.
C
      CALL VF13GR ( N , MGEQ  , MEQ   , MCON  , X     , MMAX  , C ,
     *              G, LCN, MMAX  , CN    , BL, BU, CL, CU, FIRSTG )
      IF ( DEBUG ) THEN
         WRITE( 6, 2080 )
         WRITE( 6, 2020 ) ( I, G( I ), I = 1, N )
         DO 11 J = 1, M
            WRITE( 6, 2090 ) J
            WRITE( 6, 2020 ) ( I, CN( I, J ), I = 1, N )
   11    CONTINUE
      END IF
C
C  Perform another iteration of the minimization.
C
CS      CALL VF13A ( N , M , MEQ, X, F , G , C , CN, LCN   , MAXFUN,
CD      CALL VF13AD( N , M , MEQ, X, F , G , C , CN, LCN   , MAXFUN,
     *             ACC   , IPRINT, INF   , W     , LW    , IW     )
      IF ( INF .EQ. 0 ) GO TO 10
      CALL CREPRT ( CALLS, CPU )
      CALL CNAMES( N, M, PNAME, VNAME, CNAME )
      WRITE( 6, 2110 ) F, ( I, VNAME( I ), X( I ), BL( I ), BU( I ), 
     *                      I = 1, N )
      IF ( MCON .GT. 0 ) WRITE( 6, 2120 ) ( I, CNAME( I ), C( I ), 
     *     CL( I ), CU( I ), EQUATN( I ), LINEAR( I ), I = 1, MCON )
      WRITE ( 6, 2000 ) PNAME, N , MCON, CALLS( 1 )  , CALLS( 2 ), 
     *                  CALLS( 5 ), CALLS( 6 ), F, CPU( 1 ), CPU( 2 )
      STOP
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  VF13',     /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # constraints           =      ', I10 /
     *    ,' # objective functions   =        ', F8.2 /
     *    ,' # objective gradients   =        ', F8.2 / 
     *    ,' # constraints functions =        ', F8.2 /
     *    ,' # constraints gradients =        ', F8.2 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     65('*') / )
 2010 FORMAT( /, ' Objective function value is ', 1P, D22.14 )
 2020 FORMAT( /, '     I      GRAD', /, ( I6, 1P, D12.4 ) )
 2030 FORMAT( /, ' After projection, the starting point:',
     *        /, '     I      X          BL          BU', /,
     *        ( I6, 1P, 3D12.4 ) )
 2060 FORMAT( /, ' the constraints:', /,
     *         '     I  MULTIPLIER     CL          CU      EQUALITY? ',
     *       '  LINEAR? ', /, ( I6, 1P, 3D12.4, 5X, L1, 10X, L1 ) )
 2070 FORMAT( /, ' the constraint values are:',
     *        /, '     I      C ', /, ( I6, 1P, D12.4 ) )
 2080 FORMAT( /, ' Objective function ' )
 2090 FORMAT( /, ' Constraint ', I6 )
 2100 FORMAT( /, ' Set up time = ', 0P, F12.2,
     *        /, '  Solve time = ', 0P, F12.2,
     *        /, '  Total time = ', 0P, F12.2, ' seconds' )
 2110 FORMAT( /, ' the objective function value: ', 1P, D12.4, /,
     *        /, ' the variables:', /,
     *        '     I name          value    lower bound upper bound', 
     *        /, ( I6, 1X, A10, 1P, 3D12.4 ) )
 2120 FORMAT( /, ' the constraints:', /,
     *        '     I name          value    lower bound upper bound', 
     *        ' equality?   linear? ', 
     *        /, ( I6, 1X, A10, 1P, 3D12.4, 5X, L1, 10X, L1 ) )
C
C  End of VF13MA.
C
      END
C
      SUBROUTINE VF13SE( INPUT , IOUT  , N , M , MGEQ  , MEQ    , 
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
C  Set up the input data for the the VF13 minimizer.
C
C  Nick Gould, for CGT productions,
C  7th November, 1991.
C
      INTEGER            I
CS    REAL               BIGINF
CD    DOUBLE PRECISION   BIGINF
CS    PARAMETER        ( BIGINF = 9.0E+19 )
CD    PARAMETER        ( BIGINF = 9.0D+19 )

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
      MEQ = MGEQ
C
C  If constraints have both lower and upper bounds, they have to be
C  included twice!
C
      DO 40 I = MGEQ + 1, MCON
         IF ( CL( I ) .GT. - BIGINF .AND.
     *        CU( I ) .LT.   BIGINF ) M = M + 1
   40 CONTINUE
C
C  Include any simple bounds.
C
      DO 50 I = 1, N
         IF ( BL( I ) .EQ.  BU( I ) ) THEN
            MEQ = MEQ + 1
            M   = M   + 1
         ELSE
            IF ( BL( I ) .GT. - BIGINF ) M = M + 1
            IF ( BU( I ) .LT.   BIGINF ) M = M + 1
         END IF
   50 CONTINUE
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
C  End of VF13SE.
C
      END
C
      SUBROUTINE VF13FN( N , MGEQ, MEQ, MCON, X, F, LC,
     *                   C , BL, BU, CL, CU )
      INTEGER            N , MGEQ, MEQ, MCON, LC
CS    REAL               F
CD    DOUBLE PRECISION   F
CS    REAL               X     ( N    ), C     ( LC )
CD    DOUBLE PRECISION   X     ( N    ), C     ( LC )
CS    REAL               BL    ( N    ), BU    ( N    )
CS    REAL               CL    ( LC   ), CU    ( LC   )
CD    DOUBLE PRECISION   BL    ( N    ), BU    ( N    )
CD    DOUBLE PRECISION   CL    ( LC   ), CU    ( LC   )
C
C  Evaluate the objective function and constraints.
C
C  Nick Gould, for CGT productions.
C  November 1991.
C
      INTEGER            I , MT, MFIXED, MFIXVA
CS    REAL               BIGINF
CD    DOUBLE PRECISION   BIGINF
CS    PARAMETER        ( BIGINF = 9.0E+19 )
CD    PARAMETER        ( BIGINF = 9.0D+19 )
      CALL CFN ( N , MCON , X , F , LC, C )
C
C  If there are fixed variables, shift all the inequality constraint values.
C
      MFIXED = MEQ - MGEQ
      MFIXVA = MGEQ
      IF ( MFIXED .GT. 0 ) THEN
         DO 10 I = MCON, MGEQ + 1, - 1
            C( I + MFIXED ) = C( I )
   10    CONTINUE
      END IF
C
C  If constraints have both lower and upper bounds, they have to 
C  be included twice! Reverse the signs of less-than-or-equal-to 
C  constraints.
C
      MT      = MCON + MFIXED
      DO 40 I = MGEQ + 1, MCON
         IF ( CL( I ) .GT. - BIGINF .AND.
     *        CU( I ) .LT.   BIGINF ) THEN
            MT               = MT + 1
            C( MT          ) = CU( I ) - C( I + MFIXED )
            C( I  + MFIXED ) = C( I + MFIXED ) - CL( I )
         ELSE IF ( CL( I ) .GT. - BIGINF ) THEN
            C( I + MFIXED )  = C( I + MFIXED ) - CL( I )
         ELSE IF ( CU( I ) .LT.   BIGINF ) THEN
            C( I + MFIXED )  = CU( I ) - C( I + MFIXED )
         END IF
   40 CONTINUE
C
C  Include any simple bounds, including fixed variables.
C
      DO 50 I = 1, N
         IF ( BL( I ) .EQ.  BU( I ) ) THEN
            MFIXVA      = MFIXVA + 1
            C( MFIXVA ) = X( I ) - BL( I )
         ELSE
            IF ( BL( I ) .GT. - BIGINF ) THEN
               MT       = MT + 1
               C( MT  ) = X( I ) - BL( I )
            END IF
            IF ( BU( I ) .LT.   BIGINF ) THEN
               MT      = MT + 1
               C( MT ) = BU( I ) - X( I )
            END IF
         END IF
   50 CONTINUE
      RETURN
C
C  End of VF13FN.
C
      END
C  THIS VERSION: 03/06/1994 AT 11:38:39 AM.
      SUBROUTINE VF13GR( N , MGEQ  , MEQ  , MCON  , X , LV, V , G , 
     *                   LCN   , MMAX  , CN   , BL, BU, CL, CU, FIRSTG )
      INTEGER            N , MGEQ  , MEQ   , MCON , LV     
      INTEGER            LCN   , MMAX
      LOGICAL            FIRSTG
CS    REAL               X( N ), G( N ), V( LV ), CN( LCN, MMAX )
CS    REAL               BL    ( N    ), BU    ( N    )
CS    REAL               CL    ( MMAX ), CU    ( MMAX )
CD    DOUBLE PRECISION   X( N ), G( N ), V( LV ), CN( LCN, MMAX )
CD    DOUBLE PRECISION   BL    ( N    ), BU    ( N    )
CD    DOUBLE PRECISION   CL    ( MMAX ), CU    ( MMAX )
C
C  Evaluate the gradient of the objective and constraint functions.
C
C  Nick Gould, for CGT productions,
C  November 1991.
C
      INTEGER            I , J , MT    , MFIXED, MFIXVA
CS    REAL               BIGINF, ZERO  , ONE
CD    DOUBLE PRECISION   BIGINF, ZERO  , ONE
CS    PARAMETER        ( ZERO  = 0.0E+0, ONE  = 1.0E+0 )
CD    PARAMETER        ( ZERO  = 0.0D+0, ONE  = 1.0D+0 )
CS    PARAMETER        ( BIGINF = 9.0E+19 )
CD    PARAMETER        ( BIGINF = 9.0D+19 )
C
C  Evaluate the gradient of the objective and constraint functions
C  at the initial point in a dense format.
C
      CALL CGR   ( N , MCON , X, .FALSE., MMAX, V , G ,
     *             .TRUE., LCN, MMAX , CN    )
C
C  If there are fixed variables, shift all the  gradients of the 
C  inequality constraints.
C
      MFIXED = MEQ - MGEQ
      MFIXVA = MGEQ
      IF ( MFIXED .GT. 0 ) THEN
         DO 20 I = MCON, MGEQ + 1, - 1
            DO 10 J = 1, N
               CN( J, I + MFIXED ) = CN( J, I )
   10       CONTINUE
   20    CONTINUE
      END IF
C
C  If constraints have both lower and upper bounds, their gradients
C  have to be included twice! Reverse the signs of less-than-or-equal-
C  -to constraints.
C
      MT      = MCON + MFIXED
      DO 50 I = MGEQ + 1, MCON
         IF ( CL( I ) .GT. - BIGINF .AND.
     *        CU( I ) .LT.   BIGINF ) THEN
            MT      = MT + 1
            DO 30 J = 1, N
               CN( J, MT ) = - CN( J, I + MFIXED )
   30       CONTINUE
         ELSE IF ( CU( I ) .LT.   BIGINF ) THEN
           DO 40 J = 1, N
               CN( J, I + MFIXED ) = - CN( J, I + MFIXED )
   40       CONTINUE
         END IF
   50 CONTINUE
C
C  Include the gradients of any simple bounds, including 
C  fixed variables.
C
      IF ( FIRSTG .OR. MFIXED .GT. 0 ) THEN
         DO 90 I = 1, N
            IF ( BL( I ) .EQ. BU( I ) ) THEN
               MFIXVA  = MFIXVA + 1
               DO 60 J = 1, N
                  CN( J, MFIXVA ) = ZERO
   60          CONTINUE
               CN( I, MFIXVA ) =  ONE
            ELSE
               IF ( FIRSTG ) THEN
                  IF ( BL( I ) .GT. - BIGINF ) THEN
                     MT      = MT + 1
                     DO 70 J = 1, N
                        CN( J, MT ) = ZERO
   70                CONTINUE
                     CN( I, MT ) =  ONE
                  END IF
                  IF ( BU( I ) .LT.   BIGINF ) THEN
                     MT      = MT + 1
                     DO 80 J = 1, N
                        CN( J, MT ) = ZERO
   80                CONTINUE
                     CN( I, MT ) =  - ONE
                  END IF
               END IF
            END IF
   90    CONTINUE
         FIRSTG = .FALSE.
      END IF
      RETURN
C
C  End of VF13GR.
C
      END
