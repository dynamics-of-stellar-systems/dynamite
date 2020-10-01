C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM VE09MA
C
C  --------------------------------------------------------------------
C
C  Solve the quadratic program
C
C     minimize     1/2 x(T) H x + c(T) x
C
C     subject to       A x = b,
C
C                   l <= x <= u
C
C  using the HSL code VE09.
C
C  Nick Gould.
C  December 1991.
C
C  --------------------------------------------------------------------
C
      INTEGER          NMAX, MMAX, NPMMAX, NPLUS, N, M, INSPEC
      INTEGER          INVE09, INPUT, IOUT, I, ITERN, NTOTAL
      INTEGER          ND, NS, NNZJ, NNZH, IND, IOP, MAXIT, LIP
      INTEGER          LV, LD, NDI, NT1, LWS, LIWS, II, IS, NFREE
      INTEGER          NBOTH, NNONEG, NLOWER, L, IR, IC, J, MM
      LOGICAL          CONVEX
      CHARACTER * 5    STATE
      CHARACTER * 10   PNAME
CS    REAL             ONE, ZERO, BIGINF, VL, VU, VX, VM
CD    DOUBLE PRECISION ONE, ZERO, BIGINF, VL, VU, VX, VM
      PARAMETER      ( IOUT  = 6 )
      PARAMETER      ( INPUT = 55, INSPEC = 56 )
CTOY  PARAMETER      ( NMAX = 100, MMAX = 100 )
CTOY  PARAMETER      ( LD = 10000, LWS = 10000, LIWS = 10000 )
CMED  PARAMETER      ( NMAX = 1000, MMAX = 1000 )
CMED  PARAMETER      ( LD = 100000, LWS = 100000, LIWS = 100000 )
CBIG  PARAMETER      ( NMAX = 5000, MMAX = 5000 )
CBIG  PARAMETER      ( LD = 350000, LWS = 200000, LIWS = 350000 )
CCUS  PARAMETER      ( NMAX = 3000, MMAX = 3000 )
CCUS  PARAMETER      ( LD = 200000, LWS = 150000, LIWS = 200000 )
      PARAMETER      ( NPMMAX = NMAX + MMAX, NPLUS = NPMMAX + 3 )
      INTEGER          ISTATE( NPMMAX ), IPERM( NMAX )
      INTEGER          IWS( LIWS ), IRND( LD ), IP( NPLUS )
CS    REAL             D( LD ), BND( 2, NMAX ), TUPPER, T1, T2
CS    REAL             X( NMAX ), X0( NMAX ), P( NMAX ), WS( LWS )
CS    REAL             BLOWER( NPMMAX ), BUPPER( NPMMAX )
CD    DOUBLE PRECISION D( LD ), BND( 2, NMAX ), TUPPER, T1, T2
CD    DOUBLE PRECISION X( NMAX ), X0( NMAX ), P( NMAX ), WS( LWS )
CD    DOUBLE PRECISION BLOWER( NPMMAX ), BUPPER( NPMMAX )
      LOGICAL          EQUATN( MMAX ), LINEAR( MMAX )
      CHARACTER * 10   VNAME( NMAX ), CNAME( MMAX )
CS    PARAMETER      ( ZERO  = 0.0E+0, ONE  = 1.0E+0 )
CD    PARAMETER      ( ZERO  = 0.0D+0, ONE  = 1.0D+0 )
CS    PARAMETER      ( BIGINF = 1.0E+19 )
CD    PARAMETER      ( BIGINF = 1.0D+19 )
      REAL             CPU( 2 ), CALLS( 7 )
      INTEGER          JOUT, JIN, NFRV, ITER, NR, NART
CS    REAL             F, U, GG, GMAX2, XJIN, RHO, SCALE, XART, TOLART
CS    REAL             TOLX, TOLAM, OBJF
CD    DOUBLE PRECISION F, U, GG, GMAX2, XJIN, RHO, SCALE, XART, TOLART
CD    DOUBLE PRECISION TOLX, TOLAM, OBJF
      LOGICAL          FEAS, FSEP, SOC
CS    COMMON /VE09C /  GG, U, F, RHO, GMAX2, XJIN, SCALE, XART, TOLART,
CD    COMMON /VE09CD/  GG, U, F, RHO, GMAX2, XJIN, SCALE, XART, TOLART,
     *                 TOLX, TOLAM, NART, JOUT, JIN, NFRV, ITER,
     *                 NR, FEAS, FSEP, SOC
C     
C  Open the Spec file for the method.
C
      OPEN ( INSPEC, FILE = 'VE09.SPC', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
C     MAXIT : the maximum number of iterations
C     CONVEX: true if the objective function is convex
C
      READ ( INSPEC, 1000 ) MAXIT, CONVEX
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
C  Set up the data structures necessary to hold the group partially
C  separable function.
C
      LV = MIN( LWS, MMAX )
      CALL CSETUP ( INPUT, IOUT, N, M, X0, BLOWER( 1 ),
     *              BUPPER( 1 ), NMAX, EQUATN, LINEAR,
     *              WS, BLOWER( NMAX + 1 ),
     *              BUPPER( NMAX + 1 ), LV, .FALSE., .FALSE., .FALSE. )
C
C  Determine the names of the problem, variables and constraints.
C
      MM = MAX( M, 1 )
      CALL CNAMES ( N, MM, PNAME, VNAME, CNAME )
C
C  Set up the initial estimate of the solution, SOL, and
C  right-hand-side, RHS, of the Kuhn-Tucker system.
C
C  The variables will be permuted, before being passed to
C  the optimizer, so that the first NBOTH variables are
C  bounded on both sides and the last NTOTAL - NFREE + 1 satisfy
C  nonnegativity restrictions (the remaining variables are
C  free).
C
C  On exit from the minimizer, the i-th original variable
C  will have the value X( IPERM( i ) ).
C
      NTOTAL = N
      ND     = 0
      NBOTH  = 0
      NFREE  = 0
C
C  Determine the status of each problem variable.
C
      DO 5 I = 1, N
         IF ( BUPPER( I ) .GT. BIGINF ) THEN
            IF ( BLOWER( I ) .LT. - BIGINF ) THEN
               NFREE = NFREE + 1
               IWS( I ) = 0
            ELSE IF ( BLOWER( I ) .EQ. ZERO ) THEN
               IWS( I ) = 1
            ELSE
               NBOTH = NBOTH + 1
               IWS( I ) = 2
            END IF
         ELSE
            NBOTH = NBOTH + 1
            IWS( I ) = 2
         END IF
    5 CONTINUE
C
C  Determine the status of each slack variable.
C
      DO 10 I = 1, M
         IF ( .NOT. EQUATN( I ) ) THEN
            NTOTAL = NTOTAL + 1
            IF ( ( BLOWER( NMAX + I ) .EQ. ZERO .AND.
     *             BUPPER( NMAX + I ) .GT. BIGINF ) .OR.
     *           ( BUPPER( NMAX + I ) .EQ. ZERO .AND.
     *             BLOWER( NMAX + I ) .LT. - BIGINF ) ) THEN
               IWS( NTOTAL ) = 1
            ELSE
               NBOTH = NBOTH + 1
               IWS( NTOTAL ) = 2
            END IF
         END IF
   10 CONTINUE
      NT1 = NTOTAL + 1
C
C  Permute the variables.
C
      NNONEG = NBOTH  + NFREE
      NLOWER = NNONEG + 1
      NFREE  = NBOTH
      NBOTH  = 0
      DO 15 I = 1, N
         IF ( IWS( I ) .EQ. 0 ) THEN
            NFREE = NFREE + 1
            IPERM( I ) = NFREE
         ELSE IF ( IWS( I ) .EQ. 1 ) THEN
            NNONEG = NNONEG + 1
            IPERM( I ) = NNONEG
         ELSE IF ( IWS( I ) .EQ. 2 ) THEN
            NBOTH = NBOTH + 1
            IPERM( I ) = NBOTH
            BND( 1, NBOTH ) = BLOWER( I )
            BND( 2, NBOTH ) = BUPPER( I )
         END IF
   15 CONTINUE
C
C  Introduce slack variables for inequality constraints.
C  Continue permuting the variables.
C
      IF ( NTOTAL .GT. N ) THEN
         NS = N
         DO 20 I = 1, M
            IF ( .NOT. EQUATN( I ) ) THEN
               ND = ND + 1
               NS = NS + 1
               IF ( IWS( NS ) .EQ. 1 ) THEN
                  NNONEG = NNONEG + 1
                  IPERM( NS ) = NNONEG
                  IF ( BLOWER( NMAX + I ) .EQ. ZERO ) THEN
                     D( ND ) = - ONE
                  ELSE
                     D( ND ) = ONE
                  END IF
               ELSE IF ( IWS( NS ) .EQ. 2 ) THEN
                  NBOTH = NBOTH + 1
                  IPERM( NS ) = NBOTH
                  BND( 1, NBOTH ) = BLOWER( NMAX + I )
                  BND( 2, NBOTH ) = BUPPER( NMAX + I )
                  D( ND ) = - ONE
               END IF
               IRND( ND ) = NTOTAL + I
               IWS( ND ) = IPERM( NS )
            END IF
   20    CONTINUE
      END IF
C
C  Set X to zero to determine the constant terms for
C  the problem functions.
C
      DO 30 I = 1, N
         X( I ) = ZERO
   30 CONTINUE
C
C  Evaluate the constant terms of the objective and constraint functions.
C
      CALL CFN( N , M , X, OBJF, MM, D( ND + 1 ) )
      NS = ND
      DO 40 I = 1, M
         IF ( D( NS + I ) .NE. ZERO ) THEN
            ND = ND + 1
            D( ND ) = - D( NS + I )
            IRND( ND ) = NTOTAL + I
            IWS( ND ) = NT1
         END IF
   40 CONTINUE
C
C  Evaluate the linear terms of the objective and constraint functions
C  in a sparse format.
C
      CALL CSGR ( N, M, .FALSE., LV, WS, X, NNZJ, LD - ND, D( ND + 1 ),
     *            IWS( ND + 1 ), IRND( ND + 1 ) )
      NS = ND
      DO 50 I = 1, NNZJ
         IF ( D( NS + I ) .NE. ZERO ) THEN
            ND = ND + 1
            D( ND ) = D( NS + I )
            IF ( IRND( NS + I ) .GT. 0 ) THEN
               IRND( ND ) = IRND( NS + I ) + NTOTAL
               IWS( ND ) = IPERM( IWS( NS + I ) )
            ELSE
               IRND( ND ) = IWS( NS + I )
               IWS( ND ) = NT1
            END IF
         END IF
   50 CONTINUE
C
C  Evaluate the Hessian of the Lagrangian function at the initial point.
C
      CALL CSH ( N, M, X0, LV, WS, NNZH, LD - ND, D( ND + 1 ),
     *           IRND( ND + 1 ), IWS( ND + 1 ) )
      NS = ND
      DO 60 I = 1, NNZH
         NDI = ND + I
         IF ( D( NDI ) .NE. ZERO ) THEN
            NS = NS + 1
            D( NS ) = D( NDI )
            IRND( NS ) = IPERM( IRND( NDI ) )
            IWS( NS ) = IPERM( IWS( NDI ) )
         END IF
   60 CONTINUE
      NNZH = NS - ND
      DO 70 I = 1, NNZH
         NDI = ND + I
         IF ( IRND( NDI ) .NE. IWS( NDI ) ) THEN
            NS = NS + 1
            D( NS ) = D( NDI )
            IRND( NS ) = IWS( NDI )
            IWS( NS ) = IRND( NDI )
         END IF
   70 CONTINUE
      ND = NS
C
C  Change from a co-ordinate storage scheme to a column-wise
C  scheme.
C
      IND = - 1
      LIP = NT1 + 1
CS    CALL MC49A ( IND, NT1, NTOTAL + M, ND, IRND, IWS, .TRUE., LD, D,
CD    CALL MC49AD( IND, NT1, NTOTAL + M, ND, IRND, IWS, .TRUE., LD, D,
     *             LIP, IP, LIWS - ND, IWS( ND + 1 ), IOP )
      IP( LIP + 1 ) = IP( LIP )
C
C  Prepare to call VE09.
C
      IOP    = 1
      TUPPER = ZERO
      INVE09 = 3
C
C  The major iteration.
C
      DO 300 ITERN = 1, MAXIT
CS       CALL VE09A ( D, LD, IRND, IP, M, NTOTAL, BND, NBOTH, NFREE + 1,
CD       CALL VE09AD( D, LD, IRND, IP, M, NTOTAL, BND, NBOTH, NFREE + 1,
     *                CONVEX, TUPPER, IOP, INVE09, X0, ISTATE,
     *                X, P, T1, T2, IWS, LIWS, WS, LWS )
         IF ( IOP .NE. 3 .AND. IOP .NE. 5 ) GO TO 500
  300 CONTINUE
C
C  End of main iteration loop.
C
  500 CONTINUE
      CALL CREPRT( CALLS, CPU )
C
C  Print details of the solution obtained.
C
      WRITE( 6, 2010 ) IOP
      IF ( IOP .GE. 0 ) THEN
         L = 4
         DO 520 J = 1, 2
            IF ( J .EQ. 1 ) THEN
               IR = 1
               IC = MIN( L, N )
            ELSE
               IF ( IC. LT. N - L ) WRITE( 6, 2000 )
               IR = MAX( IC + 1, N - IC + 1 )
               IC = N
            END IF
            DO 510 II = IR, IC
               I = IPERM( II )
               IS = ISTATE( I )
               IF ( IS .EQ. 0 ) THEN
                  VX = X( I + M )
                  VM = ZERO
                  STATE = ' FREE'
               ELSE
                  IF ( I .LE. NBOTH ) THEN
                     VX = BND( IS, I )
                  ELSE
                     VX = 0.0D+0
                  END IF
                  VM = - X( I + M )
                  IF ( IS .EQ. 1 ) THEN
                     STATE = 'LOWER'
                  ELSE
                     STATE = 'UPPER'
                  END IF
               END IF
               IF ( I .LE. NBOTH ) THEN
                  VL = BND( 1, I )
                  VU = BND( 2, I )
               ELSE
                  VU = BIGINF
                  IF ( I .GE. NLOWER ) THEN
                     VL = ZERO
                  ELSE
                     VL = - BIGINF
                  END IF
               END IF
               WRITE( 6, 2020 ) VNAME( II ), STATE, VX, VL, VU, VM
  510       CONTINUE
  520    CONTINUE
         IF ( M .GT. 0 ) THEN
            DO 610 I = 1, M
               WS( I ) = ZERO
  610       CONTINUE
            DO 620 I = IP( NTOTAL + 1 ), IP( NTOTAL + 2 ) - 1
               IF ( IRND( I ) .GT. NTOTAL ) THEN
                  WS( IRND( I ) - NTOTAL ) =
     *                WS( IRND( I ) - NTOTAL ) - D( I )
               END IF
  620       CONTINUE
            DO 640 II = 1, N
               J = IPERM( II )
               IS = ISTATE( J )
               IF ( IS .EQ. 0 ) THEN
                  VX = X( J + M )
               ELSE
                  IF ( I .LE. NBOTH ) THEN
                     VX = BND( IS, J )
                  ELSE
                     VX = 0.0D+0
                  END IF
               END IF
               DO 630 I = IP( J ), IP( J + 1 ) - 1
                  IF ( IRND( I ) .GT. NTOTAL ) THEN
                     WS( IRND( I ) - NTOTAL ) =
     *                  WS( IRND( I ) - NTOTAL ) + D( I ) * VX
                  END IF
  630          CONTINUE
  640       CONTINUE
            WRITE( 6, 2040 )
            L = 2
            DO 660 J = 1, 2
               IF ( J .EQ. 1 ) THEN
                  IR = 1
                  IC = MIN( L, M )
               ELSE
                  IF ( IC. LT. M - L ) WRITE( 6, 2000 )
                  IR = MAX( IC + 1, M - IC + 1 )
                  IC = M
               END IF
               DO 650 I = IR, IC
                  IF ( EQUATN( I ) ) THEN
                     STATE = 'FIXED'
                  ELSE
                     STATE = ' FREE'
                     IF ( ABS( WS( I ) - BLOWER( NMAX + I ) )
     *                    .LT. 1.0D-12 ) STATE = 'LOWER'
                     IF ( ABS( WS( I ) - BUPPER( NMAX + I ) )
     *                    .LT. 1.0D-12 ) STATE = 'UPPER'
                  END IF
                  WRITE( 6, 2020 ) CNAME( I ), STATE, WS( I ),
     *                       BLOWER( NMAX + I),
     *                       BUPPER( NMAX + I ), X( I )
  650          CONTINUE
  660       CONTINUE
         END IF
         WRITE( 6, 2030 ) F + OBJF, ITERN
      END IF
      WRITE ( 6, 2001 ) PNAME, N, M, (CALLS( I ), I = 1, 3 ),
     *                  (CALLS( I ), I = 5, 7 ),
     *                  IOP, F + OBJF, CPU( 1 ), CPU( 2 )
      CLOSE( INPUT  )
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( I10, /, L10 )
 2000 FORMAT( ' .          .....  ............',
     *        '  ............  ............  ............ ' )
 2001 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  VE09',     /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # constraints           =      ', I10 /
     *    ,' # objective functions   =        ', F8.2 /
     *    ,' # objective gradients   =        ', F8.2 / 
     *    ,' # objective Hessians    =        ', F8.2 /
     *    ,' # constraints functions =        ', F8.2 /
     *    ,' # constraints gradients =        ', F8.2 /
     *    ,' # constraints Hessians  =        ', F8.2 /
     *     ' Exit code               =      ', I10 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 2010 FORMAT( ' Stopping with IOP = ',I3 /,
     *        ' Solution:', //,
     *        ' name       state     value   ',
     *        '   Lower bound   Upper bound  Dual variable ' )
 2020 FORMAT( 1X, A10, A6, 1P, 4D14.6 )
 2030 FORMAT( /, ' Final objective function value ', 1P, D22.14,
     *        /, ' Total number of iterations = ', I6 )
 2040 FORMAT( /, ' Constraints:', //,
     *        ' name       state     value   ',
     *        '   Lower bound   Upper bound    Multiplier ' )
C
C  End of VE09MA.
C
      END
