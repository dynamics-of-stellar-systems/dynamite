C     ( Last modified on 6 Sep 2004 at 13:20:00 )
      PROGRAM          TRNMA
C
C  TRON test driver for problems derived from SIF files.
C
C  Nick Gould, for CGT Productions.
C  September 2004
C
      INTEGER          NMAX  , LH, I, IOUT, N , M , INPUT, LB,
     *                 LP, MP, LW, J, MAXIT   , L , NA   , P,
     *                 IFLAG , INSPEC, PMAX, NNZH, NNZH2, LL
CS    REAL             F, EPS, GTOL  , GNORM , ZERO, ONE,
CD    DOUBLE PRECISION F, EPS, GTOL  , GNORM , ZERO, ONE,
     *                 FRTOL, FATOL, FMIN, CGTOL, DELTA, GNORM0
      CHARACTER * 60   TASK
CCUS  PARAMETER      ( NMAX = 100000, LH = 1000000, PMAX = 25 )
CBIG  PARAMETER      ( NMAX = 100000, LH = 1000000, PMAX = 25 )
CMED  PARAMETER      ( NMAX =  10000, LH =  100000, PMAX = 25 )
CTOY  PARAMETER      ( NMAX =   1000, LH =   10000, PMAX = 25 )
      PARAMETER      ( LB = LH, LL = NMAX * PMAX )
      PARAMETER      ( IOUT  = 6 )
      INTEGER          HPTR( NMAX + 1 ), HROW( LH ), 
     *                 LPTR( NMAX + 1 ), LROW( LL ),
     *                 BPTR( NMAX + 1 ), BROW( LB ),
     *                 IFREE( NMAX ), IWA( 3 * NMAX ), ISAVE( 3 )
CS    REAL             X( NMAX ), XL( NMAX ), XU( NMAX ), G( NMAX ), 
CD    DOUBLE PRECISION X( NMAX ), XL( NMAX ), XU( NMAX ), G( NMAX ), 
     *                 XC( NMAX ), S( NMAX ), WA( 7 * NMAX ),
     *                 HVAL( LH ), HDIAG( NMAX ), LVAL( LL ), 
     *                 LDIAG( NMAX ), BVAL( LB ), BDIAG( NMAX ),
     *                 DSAVE( 3 ) 
      PARAMETER      ( INPUT = 55, INSPEC = 56 )
CS    PARAMETER      ( ZERO = 0.0E0, ONE = 1.0E0 )
CD    PARAMETER      ( ZERO = 0.0D0, ONE = 1.0D0 )
      CHARACTER * 10   XNAMES( NMAX ), PNAME, SPCDAT
      REAL             CPU( 2 ), CALLS( 4 )
CS    REAL             SGPNRM2, SNRM2
CD    DOUBLE PRECISION DGPNRM2, DNRM2
CS    EXTERNAL         STRON, SGPNRM2
CD    EXTERNAL         DTRON, DGPNRM2
C     
C  Open the Spec file for the method.
C
      SPCDAT = 'TRON.SPC'
      OPEN ( INSPEC, FILE = SPCDAT, FORM = 'FORMATTED',
     *      STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
C     P        : the number of vectors to hold factorization fill-in
C     MAXIT    : the maximum number of iterations,
C     FATOL    : the absolute error desired in the function
C     FRTOL    : the relative error desired in the function
C     FMIN     : a lower bound for the function
C     CGTOL    : the relative CG decrease required per iteration
C     GTOL     : the required norm of the projected gradient
C
      READ ( INSPEC, 1000 ) M, MAXIT, FATOL, FRTOL, FMIN, CGTOL, GTOL
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
      IF ( P .GT. PMAX ) THEN
        WRITE( IOUT, 2040 ) 'L     ', 'PMAX  ', P
        STOP
      END IF
C
C  Set up SIF data.
C
      CALL USETUP( INPUT, IOUT, N, X, XL, XU, NMAX )
C
C  Check to see if there is sufficient room for the matrices
C
      CALL UDIMSH( NNZH )
      IF ( NNZH .GT. LH ) THEN
        WRITE( IOUT, 2040 ) 'H     ', 'LH    ', NNZH
        STOP
      END IF
C
C  Obtain variable names.
C
      CALL UNAMES( N, PNAME, XNAMES )
C
C  Set up algorithmic input data.
C
      LP     = IOUT
      MP     = IOUT
      IFLAG  = 0
C
C  Optimization loo
C
      TASK = 'START'
   30 CONTINUE
C
C  Evaluate the function, F
C
         IF (TASK( 1: 1 ) .EQ. 'F' .OR. 
     *       TASK( 1: 5 ) .EQ. 'START' ) THEN
            CALL UFN( N, X, F )
         END IF
C
C  Evaluate the gradient, G, and Hessian, H. 
C  NB. The lower triangular part of H is required, 
C  and the row indices are stored in BROW
C
         IF (TASK( 1: 2 ) .EQ. 'GH' .OR. 
     *       TASK( 1: 5 ) .EQ. 'START' ) THEN
            CALL UGRSH( N, X, G, NNZH, LH, HVAL, BROW, HROW  )
C
C  Separate the diagonal of the Hessian from its off diagonal
C
            DO 40 I = 1, N
              HDIAG( I ) = ZERO
   40       CONTINUE
            NNZH2 = NNZH
            NNZH = 0
            DO 50 L = 1, NNZH2
              I = HROW( L ) 
              J = BROW( L )
              IF ( I .EQ. J ) THEN
                 HDIAG( I ) = HVAL( L ) 
              ELSE
                 NNZH = NNZH + 1
                 HROW( NNZH ) = I
                 BROW( NNZH ) = J
                 HVAL( NNZH ) = HVAL( L )
              END IF
   50       CONTINUE
C
C  Move the Hessian from co-ordinate to sparse row format
C
            CALL REORDA( N, NNZH, HROW, BROW, HVAL, HPTR, BPTR )
         END IF
C
C           Initialize the trust region bound.
C
         IF ( TASK( 1: 5 ) .EQ. 'START' ) THEN
            GNORM0 = MAX( DNRM2( N, G, 1 ), ONE )
            DELTA  = GNORM0
         END IF
C
C  Call the optimizer
C
         CALL DTRON( N, X, XL, XU, F, G, HVAL, HDIAG, HPTR, HROW,
     *               FRTOL, FATOL, FMIN, CGTOL, MAXIT, DELTA, TASK,
     *               BVAL, BDIAG, BPTR, BROW,
     *               LVAL, LDIAG, LPTR, LROW,
     *               XC, S, IFREE, ISAVE, DSAVE, WA, IWA )
C
C  Test for convergence
C 
        IF ( TASK( 1: 4 ) .EQ. 'CONV' ) THEN
           IFLAG = 0
         ELSE IF (  TASK( 1: 4 ) .EQ. 'WARN' ) THEN
           IFLAG = 1
         ELSE IF (  TASK( 1: 4 ) .EQ. 'NEWX' ) THEN
           GNORM = DGPNRM2( N, X, XL, XU, G )
           IF ( GNORM .LE. GTOL ) THEN
C          IF ( GNORM .LE. GTOL * GNORM0 ) THEN
              IFLAG = 0
           ELSE
             GO TO 30
           END IF
         ELSE
            GO TO 30
         END IF
C
C  Terminal exit.
C
      CALL UREPRT( CALLS, CPU )
      GNORM = DGPNRM2( N, X, XL, XU, G )
      WRITE ( IOUT, 2010 ) F, GNORM
      DO 120 I = 1, N
         WRITE( IOUT, 2020 ) XNAMES( I ), XL( I ), X( I ), XU( I ), 
     *                       G( I )
  120 CONTINUE
      WRITE ( IOUT, 2000 ) PNAME, N, INT( CALLS(1) ), INT( CALLS(2) ),
     *                     IFLAG, F, CPU(1), CPU(2) 
      CLOSE( INPUT  )
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( 2( I10, / ), 4( D10.3, / ), D10.3 )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  TRON',     /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # objective functions   =      ', I10 /
     *    ,' # objective gradients   =      ', I10 / 
     *     ' Exit code               =      ', I10 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 2010 FORMAT( ' Final objective function value   = ', 1P, D12.4, 
     *        /, ' Final norm of projected gradient = ', 1P, D12.4,
     *        //, '                XL           X        XU', 
     *           '           G ' )
 2020 FORMAT(  1X, A10, 1P, 4D12.4 )
 2040 FORMAT( /, ' ** ERROR from TRNMA. The declared array ', A6, 
     *           ' is too small to hold the problem data.', /, 
     *           ' Increase ', A6, ' in TRNMA to be at least ', I6, 
     *           ' and recompile. Stopping ' )
      END

      SUBROUTINE REORDA( NC, NNZ, IRN, JCN, A, IP, IW )
      INTEGER NC, NNZ
      INTEGER IRN( NNZ  ), JCN( NNZ )
      INTEGER IW( NC + 1 ), IP( NC + 1 )
CS    REAL              A( NNZ )
CD    DOUBLE PRECISION  A( NNZ )

C  Sort a sparse matrix from arbitrary order to column order

C  Nick Gould
C  7th November, 1990

      INTEGER I, J, K, L, IC, NCP1, ITEMP, JTEMP,  LOCAT
CS    REAL             ANEXT , ATEMP
CD    DOUBLE PRECISION ANEXT , ATEMP

C  Initialize the workspace as zero

      NCP1       = NC + 1
      DO 10 J    = 1, NCP1
         IW( J ) = 0
   10 CONTINUE

C  Pass 1. Count the number of elements in each column

      DO 20 K   = 1, NNZ
        J       = JCN( K )
        IW( J ) = IW( J ) + 1
   20 CONTINUE

C  Put the positions where each column begins in
C  a compressed collection into IP and IW

      IP( 1 )       = 1
      DO 30 J       = 2, NCP1
        IP( J )     = IW( J - 1 ) + IP( J - 1 )
        IW( J - 1 ) = IP( J - 1 )
   30 CONTINUE

C  Pass 2. Reorder the elements into column order. 
C          Fill in each column in turn

      DO 70 IC = 1, NC

C  Consider the next unfilled position in column IC

        DO 60 K = IW( IC ), IP( IC + 1 ) - 1

C  The entry should be placed in column J

          I       = IRN( K )
          J       = JCN( K )
          ANEXT   = A( K )
          DO 40 L = 1, NNZ

C  See if the entry is already in place

             IF ( J .EQ. IC ) GO TO 50
             LOCAT = IW( J )

C  As a new entry is placed in column J, increase the pointer 
C  IW( J ) by one

             IW( J  ) = LOCAT + 1

C  Record details of the entry which currently occupies location LOCAT

             ITEMP = IRN( LOCAT )
             JTEMP = JCN( LOCAT )
             ATEMP = A( LOCAT )

C  Move the new entry to its correct place

             IRN( LOCAT ) = I 
             JCN( LOCAT ) = J  
             A( LOCAT )   = ANEXT

C  Make the displaced entry the new entry

             I          = ITEMP
             J          = JTEMP
             ANEXT      = ATEMP
   40     CONTINUE

C  Move the new entry to its correct place 

   50     CONTINUE
          JCN( K ) = J
          IRN( K ) = I
          A( K )   = ANEXT
   60   CONTINUE
   70 CONTINUE

      RETURN

C  End of REORDA

      END
