C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM LA04MA
C
C  --------------------------------------------------------------------
C
C  Solve the linear program
C
C     minimize     1/2 x(T) H x + c(T) x
C
C     subject to       A x = b,
C
C                   l <= x <= u
C
C  using the HSL code LA04.
C
C  Nick Gould.
C  December 1991.
C
C  --------------------------------------------------------------------
C
      INTEGER          NMAX, MMAX, NPMMAX, NPLUS, N, M, INSPEC
      INTEGER          INPUT, IOUT, I, ITERN, NTOTAL, OUTSOL, NV
      INTEGER          NA, NS, NNZJ, IND, JOB, MAXIT, LIP
      INTEGER          LV, LA, NT1, LWS, LIWS, II, NFREE, IORES
      INTEGER          NBOTH, NNONEG, NLOWER, L, IR, IC, J, MM, IOUNIT
      LOGICAL          WRITES, PNAMEE
      CHARACTER * 5    STATE
      CHARACTER * 10   PNAME
      CHARACTER * 14   PNAMES
CS    REAL             ONE, ZERO, BIGINF, VL, VU, VX, VM, OBJF
CD    DOUBLE PRECISION ONE, ZERO, BIGINF, VL, VU, VX, VM, OBJF
      PARAMETER      ( IOUT  = 6 )
      PARAMETER      ( INPUT = 55, INSPEC = 56, OUTSOL = 57 )
CTOY  PARAMETER      ( NMAX = 100, MMAX = 100 )
CTOY  PARAMETER      ( LA = 10000, LWS = 10000, LIWS = 10000 )
CMED  PARAMETER      ( NMAX = 1000, MMAX = 1000 )
CMED  PARAMETER      ( LA = 100000, LWS = 100000, LIWS = 100000 )
CBIG  PARAMETER      ( NMAX =50000, MMAX = 50000 )
CBIG  PARAMETER      ( LA = 1000000, LWS = 2000000, LIWS = 3500000 )
CCUS  PARAMETER      ( NMAX = 3000, MMAX = 3000 )
CCUS  PARAMETER      ( LA = 200000, LWS = 150000, LIWS = 200000 )
      PARAMETER      ( NPMMAX = NMAX + MMAX, NPLUS = NPMMAX + 1 )
      INTEGER          IPERM( NMAX ), INVPRM( NMAX )
      INTEGER          IX( MMAX ), JX( NMAX )
      INTEGER          IWS( LIWS ), IRNA( LA ), IP( NPLUS )
CS    REAL             A( LA ), BND( 2, NMAX )
CS    REAL             X( NPMMAX ), X0( NMAX ), Z( NMAX ), CNTL( 15 )
CS    REAL             BLOWER( NPMMAX ), BUPPER( NPMMAX ), WS( LWS )
CS    REAL             C( NMAX ), B( MMAX ), G( NMAX ), RINFO( 40 )
CD    DOUBLE PRECISION A( LA ), BND( 2, NMAX )
CD    DOUBLE PRECISION X( NPMMAX ), X0( NMAX ), Z( NMAX ), CNTL( 15 )
CD    DOUBLE PRECISION BLOWER( NPMMAX ), BUPPER( NPMMAX ), WS( LWS )
CD    DOUBLE PRECISION C( NMAX ), B( MMAX ), G( NMAX ), RINFO( 40 )
      LOGICAL          EQUATN( MMAX ), LINEAR( MMAX )
      CHARACTER * 10   VNAME( NMAX ), CNAME( MMAX )
CS    PARAMETER      ( ZERO  = 0.0E+0, ONE  = 1.0E+0 )
CD    PARAMETER      ( ZERO  = 0.0D+0, ONE  = 1.0D+0 )
CS    PARAMETER      ( BIGINF = 1.0E+10 )
CD    PARAMETER      ( BIGINF = 1.0D+10 )
      REAL             CPU( 2 ), CALLS( 7 )
C     
C  Open the Spec file for the method.
C
      OPEN ( INSPEC, FILE = 'LA04.SPC', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
C     MAXIT : the maximum number of iterations
C     WRITES: write solution to "problem".sol if true
C
      READ ( INSPEC, 1000 ) MAXIT, IOUNIT, WRITES
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
C   Open the solution file if needed
C
      IF ( WRITES ) THEN
         DO 1 I = 1, 10
           IF ( PNAME( I : I ) .NE. ' ' ) THEN
             PNAMES( I : I ) = PNAME( I : I )
             J = I
           END IF
    1    CONTINUE
         PNAMES( J + 1 : J + 4 ) = '.sol'
         J = J + 4
         INQUIRE( FILE = PNAMES, EXIST = PNAMEE )
         IF ( PNAMEE ) THEN
           OPEN( OUTSOL, FILE = PNAMES, FORM = 'FORMATTED',
     *           STATUS = 'OLD', IOSTAT = IORES )
        ELSE
           OPEN( OUTSOL, FILE = PNAMES, FORM = 'FORMATTED',
     *           STATUS = 'NEW', IOSTAT = IORES )
        END IF
        IF ( IORES .NE. 0 ) THEN 
           WRITE( 6, 2160 ) IORES, PNAMES( : J )
           STOP
        END IF
      END IF
C
C  Set up the initial estimate of the solution, SOL, and
C  right-hand-side, RHS, of the Kuhn-Tucker system.
C
C
C  Set X to zero to determine the constant terms for
C  the problem functions.
C
      DO 30 I = 1, N
         X( I ) = ZERO
         C( I ) = ZERO
   30 CONTINUE
C
C  Evaluate the constant terms of the objective and constraint functions.
C
      CALL CFN( N , M , X, OBJF, MM, B )
      NS = NA
      DO 40 I = 1, M
C       B( I ) = - B( I )
        BLOWER( NMAX + I ) = BLOWER( NMAX + I ) - B( I )
        BUPPER( NMAX + I ) = BUPPER( NMAX + I ) - B( I )
        B( I ) = ZERO
   40 CONTINUE
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
      NA     = 0
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
               NA = NA + 1
               NS = NS + 1
               IF ( IWS( NS ) .EQ. 1 ) THEN
                  NNONEG = NNONEG + 1
                  IPERM( NS ) = NNONEG
                  IF ( BLOWER( NMAX + I ) .EQ. ZERO ) THEN
                     A( NA ) = - ONE
                  ELSE
                     A( NA ) = ONE
                  END IF
               ELSE IF ( IWS( NS ) .EQ. 2 ) THEN
                  NBOTH = NBOTH + 1
                  IPERM( NS ) = NBOTH
                  BND( 1, NBOTH ) = BLOWER( NMAX + I )
                  BND( 2, NBOTH ) = BUPPER( NMAX + I )
                  A( NA ) = - ONE
               END IF
               IRNA( NA ) = I
               IWS( NA ) = IPERM( NS )
            END IF
   20    CONTINUE
      END IF
      NV = NA
      WRITE( 6, "( ' n, m, k, l ', 4I7 )" ) NTOTAL, M, NBOTH, NFREE + 1
      DO 25 I = 1, NTOTAL
         INVPRM( IPERM( I ) ) = I
   25 CONTINUE
C
C  Evaluate the linear terms of the objective and constraint functions
C  in a sparse format.
C
      CALL CSGR ( N, M, .FALSE., LV, WS, X, NNZJ, LA - NA, A( NA + 1 ),
     *            IWS( NA + 1 ), IRNA( NA + 1 ) )
      NS = NA
      DO 50 I = 1, NNZJ
         IF ( A( NS + I ) .NE. ZERO ) THEN
            IF ( IRNA( NS + I ) .GT. 0 ) THEN
               NA = NA + 1
               A( NA ) = A( NS + I )
               IRNA( NA ) = IRNA( NS + I )
               IWS( NA ) = IPERM( IWS( NS + I ) )
C              IWS( NA ) = INVPRM( IWS( NS + I ) )
            ELSE
               C( IPERM( IWS( NS + I ) ) ) = A( NS + I )
C              C( INVPRM( IWS( NS + I ) ) ) = A( NS + I )
            END IF
         END IF
   50 CONTINUE
C     WRITE( 6, "( ( ' row, col, val = ', 2I8, 1PD12.4 ) )" ) 
C    *  ( IRNA( I ), IWS( I ), A( I ), I = 1, NA )
C
C  Change from a co-ordinate storage scheme to a column-wise
C  scheme.
C
      IND = - 1
      LIP = NT1 + 1
CS    CALL MC49A ( IND, NTOTAL, M, NA, IRNA, IWS, .TRUE., LA, A,
CD    CALL MC49AD( IND, NTOTAL, M, NA, IRNA, IWS, .TRUE., LA, A,
     *             LIP, IP, LIWS - NA, IWS( NA + 1 ), JOB )
      IP( LIP + 1 ) = IP( LIP )
C
C  Initialize CNTL
C
CS    CALL LA04I ( CNTL )
CD    CALL LA04ID( CNTL )
      CNTL( 7 ) = IOUNIT
C
C  Prepare to call LA04
C
      JOB = 1
C
C  Loop over a sequence of simplex iterations
C
      DO 300 ITERN = 1, MAXIT
CS       CALL LA04A ( A, LA, IRNA, IP, M, NTOTAL, B, C, BND, NBOTH, 
CD       CALL LA04AD( A, LA, IRNA, IP, M, NTOTAL, B, C, BND, NBOTH, 
     *                NFREE + 1, JOB, CNTL, IX, JX, X, Z, G, RINFO,
     *                WS, LWS, IWS, LIWS )
         IF ( JOB .EQ. 0 ) GO TO 400
         IF ( JOB .LT. 0 ) GO TO 990
  300 CONTINUE
C
C  End of main iteration loop.
C
  400 CONTINUE
      CALL CREPRT( CALLS, CPU )
C
C  Print details of the solution obtained.
C
      WRITE( 6, 2010 ) JOB
      IF ( WRITES ) WRITE( OUTSOL, 2010 ) JOB
      IF ( JOB .GE. 0 ) THEN
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
            DO 510 I = IR, IC
               II = INVPRM( I )
               VM = Z( II )
               IF ( II .LE. NBOTH ) THEN
                  IF ( JX( II ) .LE. 0 ) THEN
                     VX = X( II )
                     STATE = ' FREE'
                  ELSE IF ( JX( II ) .EQ. 1 ) THEN
                     VX = X( II ) + BND( 1, II )
                     STATE = 'LOWER'
                  ELSE
                     VX = X( II ) + BND( 2, II )
                     STATE = 'UPPER'
                  END IF
               ELSE
                  VX = X( II )
                  STATE = ' FREE'
                  IF ( I .GT. NFREE .AND. 
     *                 ABS( VX ) .LT. 1.0D-12 ) STATE = 'LOWER'
               END IF
               IF ( II .LE. NBOTH ) THEN
                  VL = BND( 1, II )
                  VU = BND( 2, II )
               ELSE
                  VU = BIGINF
                  IF ( I .GE. NLOWER ) THEN
                     VL = ZERO
                  ELSE
                     VL = - BIGINF
                  END IF
               END IF
               WRITE( 6, 2020 ) VNAME( I ), STATE, VX, VL, VU, VM
  510       CONTINUE
  520    CONTINUE
         IF ( WRITES ) THEN
            DO 540 I = 1, N
               II = INVPRM( I )
               VM = Z( II )
               IF ( II .LE. NBOTH ) THEN
                  IF ( JX( II ) .LE. 0 ) THEN
                     VX = X( II )
                     STATE = ' FREE'
                  ELSE IF ( JX( II ) .EQ. 1 ) THEN
                     VX = X( II ) + BND( 1, II )
                     STATE = 'LOWER'
                  ELSE
                     VX = X( II ) + BND( 2, II )
                     STATE = 'UPPER'
                  END IF
               ELSE
                  VX = X( II )
                  STATE = ' FREE'
                  IF ( I .GT. NFREE .AND. 
     *                 ABS( VX ) .LT. 1.0D-12 ) STATE = 'LOWER'
               END IF
               IF ( II .LE. NBOTH ) THEN
                  VL = BND( 1, II )
                  VU = BND( 2, II )
               ELSE
                  VU = BIGINF
                  IF ( I .GE. NLOWER ) THEN
                     VL = ZERO
                  ELSE
                     VL = - BIGINF
                  END IF
               END IF
               WRITE( OUTSOL, 2020 ) VNAME( I ), STATE, VX, VL, VU, VM
  540    CONTINUE
         END IF
         IF ( M .GT. 0 ) THEN
C
C  Compute the Lagrange multipliers
C
            DO 570 J = 1, M
               I = IX( J )
               IF ( I .GT. 2 * ( N + M ) ) THEN
                 I = I - 2 * ( N + M ) 
               ELSE IF ( I .GT. ( N + M ) ) THEN
                 I = I - ( N + M ) 
               END IF
               WS( J ) = C( I )
  570       CONTINUE
            J = 7
CS          CALL LA04A ( A, LA, IRNA, IP, M, NTOTAL, B, C, BND, NBOTH, 
CD          CALL LA04AD( A, LA, IRNA, IP, M, NTOTAL, B, C, BND, NBOTH, 
     *                   NFREE + 1, J, CNTL, IX, JX, X, Z, G, RINFO,
     *                   WS, LWS, IWS, LIWS )

            DO 580 J = 1, M
               X( NTOTAL + J ) = WS( J )
  580       CONTINUE
C
C  Now compute the constrainmt residuals
C
            DO 610 I = 1, M
               WS( I ) = ZERO
  610       CONTINUE
            DO 640 J = 1, NTOTAL
               IF ( INVPRM( J ) .LE. N ) THEN
               IF ( J .LE. NBOTH ) THEN
                  IF ( JX( J ) .LE. 0 ) THEN
                     VX = X( J )
                  ELSE IF ( JX( J ) .EQ. 1 ) THEN
                     VX = X( J ) + BND( 1, J )
                  ELSE
                     VX = X( J ) + BND( 2, J )
                     STATE = 'UPPER'
                  END IF
               ELSE
                  VX = X( J )
               END IF
               DO 630 I = IP( J ), IP( J + 1 ) - 1
                  WS( IRNA( I ) ) = WS( IRNA( I ) ) + A( I ) * VX
  630          CONTINUE
               END IF
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
                     STATE = 'EQUAL'
                  ELSE
                     STATE = ' FREE'
                     IF ( ABS( WS( I ) - BLOWER( NMAX + I ) )
     *                    .LT. 1.0D-12 ) STATE = 'LOWER'
                     IF ( ABS( WS( I ) - BUPPER( NMAX + I ) )
     *                    .LT. 1.0D-12 ) STATE = 'UPPER'
                  END IF
                  WRITE( 6, 2020 ) CNAME( I ), STATE, WS( I ),
     *                       BLOWER( NMAX + I),
     *                       BUPPER( NMAX + I ), X( NTOTAL + I )
  650          CONTINUE
  660       CONTINUE
            IF ( WRITES ) THEN
               WRITE( OUTSOL, 2040 )
               DO 680 I = 1, M
                  IF ( EQUATN( I ) ) THEN
                     STATE = 'EQUAL'
                  ELSE
                     STATE = ' FREE'
                     IF ( ABS( WS( I ) - BLOWER( NMAX + I ) )
     *                    .LT. 1.0D-12 ) STATE = 'LOWER'
                     IF ( ABS( WS( I ) - BUPPER( NMAX + I ) )
     *                    .LT. 1.0D-12 ) STATE = 'UPPER'
                  END IF
                  WRITE( OUTSOL, 2020 ) CNAME( I ), STATE, WS( I ),
     *                       BLOWER( NMAX + I),
     *                       BUPPER( NMAX + I ), X( NTOTAL + I )
  680          CONTINUE
            END IF
         END IF
         WRITE( 6, 2030 ) RINFO( 1 ) + OBJF, ITERN
         IF ( WRITES ) WRITE( OUTSOL, 2030 ) RINFO( 1 ) + OBJF, ITERN
      END IF
      WRITE ( 6, 2001 ) PNAME, N, M, (CALLS( I ), I = 1, 3 ),
     *                  (CALLS( I ), I = 5, 7 ),
     *                  JOB, RINFO( 1 ) + OBJF, CPU( 1 ), CPU( 2 )
      CLOSE( INPUT  )
      STOP

  990 CONTINUE
      WRITE( 6, 2050 ) JOB
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( I10, /, I10, /, L10 )
 2000 FORMAT( ' .          .....  ............',
     *        '  ............  ............  ............ ' )
 2001 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  LA04',     /
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
 2010 FORMAT( ' Stopping with JOB = ',I3 /,
     *        ' Solution:', //,
     *        ' name       state     value   ',
     *        '   Lower bound   Upper bound  Dual variable ' )
 2020 FORMAT( 1X, A10, A6, 1P, 4D14.6 )
 2030 FORMAT( /, ' Final objective function value ', 1P, D22.14,
     *        /, ' Total number of iterations = ', I6 )
 2040 FORMAT( /, ' Constraints:', //,
     *        ' name       state     value   ',
     *        '   Lower bound   Upper bound    Multiplier ' )
 2050 FORMAT( /, ' Error return from LA04. Job = ', I5 )
 2160 FORMAT( ' IOSTAT = ', I6, ' when opening file ', A14, 
     *        '. Stopping ' )
C
C  End of LA04MA.
C
      END
