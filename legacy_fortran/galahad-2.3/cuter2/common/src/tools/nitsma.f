C     ( Last modified on 12 Sepc 2004 at 09:50:12 )
C  Correction: 12/Sep/2004: undeclared integer variable declared
C  ------------------------------------------------------
C
C  Solve a system of nonlinear equations using 
C  Homer Walker's package NITSOL version 0.3 
C  ( ftp://ephraim.wpi.edu/pub/nitsol/nitsol.tar.gz )
C
C  CUTEr interface: Nick Gould
C  June 2003
C
C  ------------------------------------------------------
C
      PROGRAM NITSMA

      INTEGER NMAX, LIPAR, LRPAR, IOUT, INSPEC, LWORK, NWORK, I
      INTEGER NN, MM, INPUT, M
      PARAMETER ( INPUT =  55, IOUT = 6, INSPEC = 46 )
CTOY  PARAMETER ( NMAX = 100 )
CMED  PARAMETER ( NMAX = 1000 )
CBIG  PARAMETER ( NMAX = 50000 )
CCUS  PARAMETER ( NMAX = 3000 )
      PARAMETER ( LWORK = 100 * NMAX )
      PARAMETER ( LIPAR = NMAX + 2, LRPAR = 2 * NMAX )
CS    REAL               ZERO, BIGINF
CD    DOUBLE PRECISION   ZERO, BIGINF
CS    PARAMETER        ( ZERO = 0.0E+0, BIGINF = 1.0E+19 )
CD    PARAMETER        ( ZERO = 0.0D+0, BIGINF = 1.0D+19 )
      INTEGER N, ITERM
CS    REAL             FTOL, STPTOL, FINAL, F
CD    DOUBLE PRECISION FTOL, STPTOL, FINAL, F
      INTEGER NINPUT( 10 ), INFO( 6 ), IPAR( LIPAR )
CS    REAL             X( NMAX ), WORK( LWORK ), RPAR( LRPAR ),
CD    DOUBLE PRECISION X( NMAX ), WORK( LWORK ), RPAR( LRPAR ),
     *                 XFREE( NMAX )
CS    REAL             SDOT, SNRM2
CD    DOUBLE PRECISION DDOT, DNRM2
      LOGICAL          EQUATN( NMAX ), LINEAR( NMAX )
CS    EXTERNAL CUTEF, CUTEJ, CUTEFN, CUTEJN, SDOT, SNRM2
CD    EXTERNAL CUTEF, CUTEJ, CUTEFN, CUTEJN, DDOT, DNRM2

      INTEGER NFREE, NNIMAX, IJACV, IKRYSL, KDMAX, IRPRE
      INTEGER IKSMAX, IRESUP, IFDORD, IBTMAX, IETA, IPSOL
      REAL CPU( 2 ), CALLS( 6 )
      CHARACTER * 10   PNAME, VNAME( NMAX ), CNAME( NMAX )
C
C  NITSOL printing common block (see nitprint.h)
C
      INTEGER IPLVL, IPUNIT
      COMMON / NITPRINT / IPLVL, IPUNIT
C
C  NITSOL info common block (see nitprint.h)
C
      INTEGER INSTEP, NEWSTEP, KRYSTAT
CS    REAL             AVRATE, FCURNRM, ETA
CD    DOUBLE PRECISION AVRATE, FCURNRM, ETA
      COMMON / NITINFO / AVRATE, FCURNRM, ETA, INSTEP, NEWSTEP, KRYSTAT
C
C  common block to tell when to evluate Jacobian again
C
      LOGICAL JKNOWN
      COMMON / NITJEV / JKNOWN
      JKNOWN = .FALSE.
C     
C  Open the Spec file for the method.
C
      OPEN ( INSPEC, FILE = 'NITSOL.SPC', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
C  NNIMAX = maximum number of nonlinear iterations 
C  IJACV  = the method of J*v evaluation 
C          (0 => finite differences, 1 => analytic)
C  IKRYSL = the Krylov solver: 
C          (0 => GMRES 1 => BiCGSTAB, 2 => TFQMR)
C  KDMAX  = max Krylov subspace dimension when GMRES is used 
C  IRPRE = flag for right preconditioning: 
C          ( 0 => no right preconditioning, 1 => right preconditioning)
C  IKSMAX = max allowable number iterations per call to Krylov solver
C  IRESUP = residual update flag when GMRES is used
C           ( 0 => linear combination, 1 => direct evaluation)
C  IFDORD = order of the finite-difference formula (sometimes) used 
C  IBTMAX = maximum allowable number of backtracks per linesearch
C  IETA   = flag determining the forcing term eta
C           ( 0 => abs( ||fcur|| - ||fprev+Jprev*sprev|| )/||fprev||,
C             1 => (||fcur||/||fprev||)**2,
C             2 => gamma*(||fcur||/||fprev||)**alpha 
C             3 => fixed (constant) eta in (0,1) )
C  IPLVL  = printlevel 
C           ( 0 => no printout, 
C             1 => iteration numbers and F-norms,
C             2 => ... + some stats, step norms, and linear model norms
C             3 => ... + some Krylov solver and backtrack information
C             4 => ... + more Krylov solver and backtrack information)
C  IPUNIT = printout unit number
C  IPSOL  = print solution to output channel
C           ( 0 => no, 1=> yes )
C  FTOL   = stopping tolerance on the f-norm
C  STPTOL = stopping tolerance on the steplength
C
      READ( INSPEC, "( I10, 12( /, I10 ), 2( /, 1P, D10.3 ) )" ) 
     *   NNIMAX, IJACV, IKRYSL, KDMAX, IRPRE, IKSMAX, IRESUP, 
     *   IFDORD, IBTMAX, IETA, IPLVL, IPUNIT, IPSOL, FTOL, STPTOL
C
C  Assign values
C
      NINPUT(  1 ) = NNIMAX
      NINPUT(  2 ) = IJACV
      NINPUT(  3 ) = IKRYSL
      NINPUT(  4 ) = KDMAX
      NINPUT(  5 ) = IRPRE
      NINPUT(  6 ) = IKSMAX
      NINPUT(  7 ) = IRESUP
      NINPUT(  8 ) = IFDORD
      NINPUT(  9 ) = IBTMAX
      NINPUT( 10 ) = IETA
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
C  Determine the number of variables and constraints
C
      CALL CDIMEN( INPUT, N, M )
C
C  Check that there is sufficient space
C
      IF ( N .LE. 0 ) THEN
         WRITE( IOUT, "( ' n <= 0 ' )" )
         STOP
      END IF
C
C  Check that there is sufficient space
C
      IF ( N .GT. NMAX ) THEN
         WRITE( IOUT, "( ' Workspace array X of length ', I10, 
     *                   'too small', /, ' Please increase NMAX in',
     *                   ' nitsma.f to have length at least ', I10, 
     *                   ' and recompile ' )" ) NMAX, N
         STOP
      END IF
C
C  Check workspace array is long enough
C
      IF ( NINPUT( 3 ) .EQ. 1 ) THEN
         NWORK = N * ( NINPUT( 4 ) + 5 ) + 
     *        NINPUT( 4 ) * ( NINPUT( 4 ) + 3  )
      ELSE IF ( NINPUT( 3 ) .EQ. 2 ) THEN
         NWORK = 11 * N
      ELSE
         NWORK = 14 * N
      END IF
      IF ( LWORK .LT. NWORK ) THEN
         WRITE( IOUT, "( ' Workspace array WORK of length ', I10, 
     *                   'too small', /, ' Please increase LWORK in',
     *                   ' nitsma.f to have length at least ', I10, 
     *                   ' and recompile ' )" ) LWORK, NWORK
         STOP
      END IF 
C
C  Set up the data structures necessary to hold the group partially
C  separable function.
C
      NN = N
      MM = M
      CALL CSETUP( INPUT, IOUT, N, M, X, WORK( 1 ),
     *             WORK( NN + 1 ), NMAX, EQUATN, LINEAR,
     *             WORK( 2 * NN + 2 * MM + 1 ), 
     *             WORK( 2 * NN + 1 ),
     *             WORK( 2 * NN + MM + 1 ), MM, 
     *             .FALSE., .FALSE., .FALSE. )
C
C  Determine the names of the problem, variables and constraints.
C
      CALL CNAMES( N, M, PNAME, VNAME, CNAME )
C
C  Check that there are no variable bounds
C
      NFREE = 0
      DO 10 I = 1, N
         IF ( WORK( I ) .GT. - BIGINF .OR. 
     *       WORK( N + I ) .LT. BIGINF ) THEN
            IF ( WORK( I ) .NE. WORK( N + I ) ) THEN
               WRITE( IOUT, "( ' Variable ', A10, ' is bounded ',
     *                      /, ' so NITSOL is not appropriate ' )" )  
     *            VNAME( I )
               STOP
            END IF
         ELSE
           NFREE = NFREE + 1
         END IF
   10 CONTINUE
C
C  Check that all constraints are equalities
C
      DO 20 I = 1, M
         IF ( WORK( 2 * N + I ) .NE.  WORK( 2 * N + M + I ) ) THEN
            WRITE( IOUT, "( ' Constraint ', A10, ' is an inequality ',
     *                   /, ' so NITSOL is not appropriate ' )" )  
     *         CNAME( I )
            STOP
         END IF
   20 CONTINUE
C
C  Check that the system is "square"
C
      IF ( NFREE .NE. M ) THEN
         WRITE( IOUT, "( ' n = ', I10, ' /= ', ' m = ', I10,
     *      /, ' so NITSOL is not appropriate ' )" ) NFREE, M    
         STOP
      END IF
C
C  Solve the problem
C
C  No fixed variable case
C
      IF ( NFREE .EQ. N ) THEN
         CALL NITSOL( N, X, CUTEF, CUTEJ, FTOL, STPTOL, NINPUT, 
CS   *                INFO, WORK, RPAR, IPAR, ITERM, SDOT, SNRM2 )
CD   *                INFO, WORK, RPAR, IPAR, ITERM, DDOT, DNRM2 )
C
C  Fixed variable case
C
      ELSE
         IPAR( 1 ) = N
         IPAR( 2 ) = M
         NFREE = 0
         DO 30 I = 1, N
            IF ( WORK( I ) .NE. WORK( N + I ) ) THEN
               NFREE = NFREE + 1
               IPAR( NFREE + 2 ) = I
               XFREE( NFREE ) = X( I )
            ELSE
               X( I ) = WORK( I )
               RPAR( I ) = X( I )
               RPAR( N + I ) = ZERO
            END IF
   30    CONTINUE
         CALL NITSOL( NFREE, XFREE, CUTEFN, CUTEJN, FTOL, STPTOL, 
CS   *           NINPUT, INFO, WORK, RPAR, IPAR, ITERM, SDOT, SNRM2 )
CD   *           NINPUT, INFO, WORK, RPAR, IPAR, ITERM, DDOT, DNRM2 )
         DO 40 I = 1, NFREE
            X( IPAR( I + 2 ) ) = XFREE( I )
   40    CONTINUE
      END IF
C
C  Write results
C
      CALL CREPRT ( CALLS, CPU )
      CALL CFN( N, M, X, F, M, WORK )
CS    FINAL = SNRM2( M, WORK, 1 )
CD    FINAL = DNRM2( M, WORK, 1 )
      WRITE( IOUT, "( /, ' Termination flag iterm:       ', I9,
     *             /, ' Final f-norm:                 ', 1P, E9.3,
     *             /, ' No. function evaluations:     ', I9, 
     *             /, ' No. J*v evaluations:          ', I9,
     *             /, ' No. P(inverse)*v evaluations: ', I9, 
     *             /, ' No. linear iterations:        ', I9, 
     *             /, ' No. nonlinear iterations:     ', i9,
     *             /, ' No. backtracks:               ', i9 )" )
     *      ITERM, FINAL, INFO( 1 ), INFO( 2 ), INFO( 3 ), 
     *      INFO( 4 ), INFO( 5 ), INFO( 6 )
      IF ( IPSOL .GT. 0 ) THEN
        WRITE( IOUT, "( /, ' the variables:', /,
     *          '     I name          value',
     *          /, ( I6, 1X, A10, 1P, D12.4 ) )" )
     *       ( I, VNAME( I ), X( I ), I = 1, N )
        WRITE( IOUT, "( /, ' the constraints:', /,
     *          '     I name          value',
     *          /, ( I6, 1X, A10, 1P, D12.4 ) )" )
     *     (   I, CNAME( I ), WORK( I ), I = 1, M )
      END IF
      WRITE( 6, "( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  NITSOL',     /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables (inc.fixed) =      ', I10 /
     *    ,' # equations             =      ', I10 /
     *    ,' # objective functions   =        ', F8.2 /
     *    ,' # objective gradients   =        ', F8.2 / 
     *    ,' # constraints functions =        ', F8.2 /
     *    ,' # constraints gradients =        ', F8.2 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     65('*') / )" )
     *  PNAME, N , M, CALLS( 1 ), CALLS( 2 ), 
     *  CALLS( 5 ), CALLS( 6 ), FINAL, CPU( 1 ), CPU( 2 )
      CLOSE( INPUT  )
      STOP
      END

      SUBROUTINE CUTEF( N, X, C, RPAR, IPAR, ITRMF )
      INTEGER N, ITRMF
      INTEGER IPAR( * )
CS    REAL             F, X( N ), C( N ), RPAR( * )
CD    DOUBLE PRECISION F, X( N ), C( N ), RPAR( * )
      LOGICAL JKNOWN
      COMMON / NITJEV / JKNOWN
      CALL CFN( N, N, X, F, N, C )
      JKNOWN = .FALSE.
      ITRMF = 0
      RETURN
      END

      SUBROUTINE CUTEJ( N, X, C, IJOB, V, Z, RPAR, IPAR, ITRMJV )
      INTEGER N, IJOB, ITRMJV
      INTEGER IPAR( * )
CS    REAL             X( N ), C( N ), V( N ), Z( N ), RPAR( * )
CD    DOUBLE PRECISION X( N ), C( N ), V( N ), Z( N ), RPAR( * )
      LOGICAL JKNOWN
      COMMON / NITJEV / JKNOWN
      IF ( IJOB .EQ. 0 ) THEN
         CALL CJPROD( N, N, JKNOWN, .FALSE., X, V, N, Z, N )
         JKNOWN = .TRUE.
         ITRMJV = 0
      ELSE
         ITRMJV = 2
      END IF
      RETURN
      END

      SUBROUTINE CUTEFN( NFREE, XFREE, C, RPAR, IPAR, ITRMF )
      INTEGER NFREE, ITRMF
      INTEGER IPAR( * )
CS    REAL             F, XFREE( NFREE ), C( NFREE ), RPAR( * )
CD    DOUBLE PRECISION F, XFREE( NFREE ), C( NFREE ), RPAR( * )
      INTEGER N, M, I
      LOGICAL JKNOWN
      COMMON / NITJEV / JKNOWN
      N = IPAR( 1 )
      M = IPAR( 2 )
      DO 10 I = 1, NFREE
         RPAR( IPAR( I + 2 ) ) = XFREE( I )
   10 CONTINUE
      CALL CFN( N, M, RPAR( 1 ), F, M, C )
      JKNOWN = .FALSE.
      ITRMF = 0
      RETURN
      END

      SUBROUTINE CUTEJN( NFREE, XFREE, C, IJOB, VFREE, Z, RPAR, 
     *                   IPAR, ITRMJV )
      INTEGER NFREE, IJOB, ITRMJV
      INTEGER IPAR( * )
CS    REAL             XFREE( NFREE ), C( NFREE ), VFREE( NFREE ), 
CD    DOUBLE PRECISION XFREE( NFREE ), C( NFREE ), VFREE( NFREE ), 
     *                 Z( NFREE ), RPAR( * )
      INTEGER I, N, M
      LOGICAL JKNOWN
      COMMON / NITJEV / JKNOWN
      IF ( IJOB .EQ. 0 ) THEN
         N = IPAR( 1 )
         M = IPAR( 2 )
         DO 10 I = 1, NFREE
            RPAR( IPAR( I + 2 ) ) = XFREE( I )
            RPAR( N + IPAR( I + 2 ) ) = VFREE( I )
   10    CONTINUE
         CALL CJPROD( N, M, JKNOWN, .FALSE., RPAR( 1 ), RPAR( N + 1 ), 
     *                N, Z, M )
         JKNOWN = .TRUE.
         ITRMJV = 0
      ELSE
         ITRMJV = 2
      END IF
      RETURN
      END

