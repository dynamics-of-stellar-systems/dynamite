C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM VE14MA
C
C  Main program for VE14, an algorithm for solving bound-constrained 
C  QPs via barrier function methods.
C
C  Nick Gould, Rutherford Appleton Laboratory.
C  March 1992.
C  Modified by Ph. Toint.
C
C  Local variable definitions.

      INTEGER          NMAX  , N, NNZH, LWS, LIWS, IOUT , IPRINT
      INTEGER          ITYPEB, I, J, L, IR , IC  , INPUT, NS, LH
      INTEGER          INSPEC
      CHARACTER * 5    STATE
      CHARACTER * 10   PNAME
      PARAMETER      ( IPRINT = 2, IOUT  = 6 )
CS    REAL             ZERO, OBJF, STOPG, QFVAL
CD    DOUBLE PRECISION ZERO, OBJF, STOPG, QFVAL
      PARAMETER      ( INPUT = 55, INSPEC = 56 )
CTOY  PARAMETER      ( NMAX = 100 )
CMED  PARAMETER      ( NMAX = 1000 )
CBIG  PARAMETER      ( NMAX = 20000 )
CCUS  PARAMETER      ( NMAX = 5000 )
CTOY  PARAMETER      ( LH = 10000  , LWS = 10000  , LIWS = 10000 )
CMED  PARAMETER      ( LH = 100000 , LWS = 100000 , LIWS = 100000 )
CBIG  PARAMETER      ( LH = 300000 , LWS = 600000 , LIWS = 600000 )
CCUS  PARAMETER      ( LH = 200000 , LWS = 300000 , LIWS = 300000 )
CS    PARAMETER      ( ZERO   = 0.0E+0 )
CD    PARAMETER      ( ZERO   = 0.0D+0 )
C
C  Array definitions.
C
      INTEGER          ICNTRL( 4 ), INFORM( 7 )
      INTEGER          IRNH( LH ) , ICNH( LH ) , IWS( LIWS )
CS    REAL             X( NMAX ), G( NMAX ), C( NMAX )
CS    REAL             H( LH ), BND( 2, NMAX ), WS( LWS ), CONTRL( 4 )
CD    DOUBLE PRECISION X( NMAX ), G( NMAX ), C( NMAX )
CD    DOUBLE PRECISION H( LH ), BND( 2, NMAX ), WS( LWS ), CONTRL( 4 )
      LOGICAL          LCNTRL( 2 )
      CHARACTER * 10   VNAME( NMAX )
      REAL             CPU( 2 ), CALLS( 4 )
C     
C  Open the Spec file for the method.
C
      OPEN ( INSPEC, FILE = 'VE14.SPC', FORM = 'FORMATTED',
     *      STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
C     ITYPEB indicates the type of barrier used. 
C            ITYPEB = 1 is a normal barrier function, 
C            ITYPEB = 2 is the Jittorntrum and Osborne variation,
C            ITYPEB = 3 is a Lagrangian barrier function.
C     STOPG  gives the stopping tolerance on the projected gradient
C 
      READ ( INSPEC, 1000 ) ITYPEB, STOPG
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
      CALL USETUP ( INPUT, IOUT, N, X, WS( 1 ), WS( NMAX + 1 ), NMAX )
C
C  Determine the names of the problem and variables.
C
      CALL UNAMES( N, PNAME, VNAME )
C
C  Allocate the bounds.
C
CS    CALL SCOPY ( N, WS( 1 ), 1, BND( 1, 1 ), 2 )
CD    CALL DCOPY ( N, WS( 1 ), 1, BND( 1, 1 ), 2 )
CS    CALL SCOPY ( N, WS( NMAX + 1 ), 1, BND( 2, 1 ), 2 )
CD    CALL DCOPY ( N, WS( NMAX + 1 ), 1, BND( 2, 1 ), 2 )
      DO 10 I = 1, N
         WS( I ) = ZERO
   10 CONTINUE   
C
C  Evaluate the constant term, OBJF, of the objective function.
C
      CALL UFN ( N, WS, OBJF )
C
C  Evaluate the linear term of the objective function.
C
      CALL UGR ( N, WS, C )
C
C  Evaluate the Hessian of the objective function at the initial point.
C
      CALL USH ( N, WS, NNZH, LH, H, IRNH, ICNH )
C
C  Remove any zeros. 
C
      NS = 0
      DO 20 I = 1, NNZH
         IF ( H( I ) .NE. ZERO ) THEN
            NS = NS + 1
            IRNH( NS ) = IRNH( I )
            ICNH( NS ) = ICNH( I )
            H( NS ) = H( I )
         END IF
   20 CONTINUE
      NNZH = NS
C
C  ------------------- problem set-up complete ----------------------
C
C
C  Set all default values.
C
CS    CALL VE14I ( ICNTRL, CONTRL, LCNTRL )
CD    CALL VE14ID( ICNTRL, CONTRL, LCNTRL )
C
C  Call the optimizer.
C
CS    CALL VE14A ( ITYPEB, IPRINT, STOPG , N , X , BND   ,
CD    CALL VE14AD( ITYPEB, IPRINT, STOPG , N , X , BND   ,
     *             OBJF  , C     , NNZH  , LH, H , IRNH  , ICNH   ,
     *             G     , IWS   , LIWS  , WS( N + 1 )   , LWS - N, 
     *             QFVAL , ICNTRL, CONTRL, LCNTRL, INFORM )
      CALL UREPRT ( CALLS, CPU )
C
C  Print details of the solution obtained.
C
      WRITE( IOUT, 2030 ) INFORM( 1 )
      IF ( INFORM( 1 ) .GE. 0 .AND. INFORM( 1 ) .LE. 3 ) THEN
         L = 2
         DO 520 J = 1, 2
            IF ( J .EQ. 1 ) THEN
               IR = 1
               IC = MIN( L, N )
            ELSE
               IF ( IC. LT. N - L ) WRITE( IOUT, 2050 )
               IR = MAX( IC + 1, N - IC + 1 )
               IC = N
            END IF
            DO 510 I = IR, IC
               STATE = ' FREE'
               IF ( ABS( X( I ) - BND( 1, I ) ) .LT. 1.0D-12 ) 
     *            STATE = 'LOWER' 
               IF ( ABS( X( I ) - BND( 2, I ) ) .LT. 1.0D-12 ) 
     *            STATE = 'UPPER' 
               IF ( ABS( BND( 1, I ) - BND( 2, I ) ) .LT. 1.0D-12 ) 
     *            STATE = 'FIXED' 
               WRITE( IOUT, 2040 ) VNAME( I ), STATE, X( I ), 
     *                          BND( 1, I ), BND( 2, I ), G( I )
  510       CONTINUE
  520    CONTINUE
      END IF
C
C   Print statistics
C

      WRITE ( IOUT, 2000 ) PNAME, N, ( CALLS( I ), I = 1, 3 ),
     *                     INFORM( 1 ), QFVAL, CPU( 1 ), CPU( 2 )           
      CLOSE( INPUT  )
      STOP
C
C   Non executable statements
C
 1000 FORMAT( I10, /, D10.3 )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  VE14',     /
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
 2030 FORMAT( /, ' Stopping with INFORM( 1 ) = ', I3, //,
     *        ' Solution:', //, 
     *        ' name       state     value   ',
     *        '   Lower bound   Upper bound  Dual variable ' )
 2040 FORMAT( 1X, A10, A6, 1P, 4D14.6 )
 2050 FORMAT( ' .          .....  ............',
     *        '  ............  ............  ............ ' )
C
C  End of VE14MA.
C
      END
