C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM          PDSMA
C
C  PDS test driver for problems derived from SIF files.
C
C  A. R Conn and Ph. Toint for CGT Productions.
C  January 1995, substantially modified September 1996.
C
      INTEGER          NMAX  , SCH, INPUT, INDR, IOUT
      INTEGER          DIM   , RES, IMAX 
CTOY  PARAMETER      ( NMAX  =  10, DIM  =   10, IMAX = 2000 )
CMED  PARAMETER      ( NMAX  = 100, DIM  =   20, IMAX = 4000 )
CBIG  PARAMETER      ( NMAX  = 500, DIM  =   30, IMAX = 6000 )
CCUS  PARAMETER      ( NMAX  = 200, DIM  =   25, IMAX = 5000 )
      PARAMETER      ( INPUT =  55, INDR =   46, RES  =   56 )
      PARAMETER      ( IOUT  =   6, SCH  =   48 )
      INTEGER          CNT   , DEBUG, ERROR , UNIQUE, MAXITR, N
      INTEGER          IFACT , TYPE , RESIZE, SSS   , LIMIT
      PARAMETER      ( LIMIT = ( DIM + 2 ) * IMAX )
      INTEGER          SCHEME( LIMIT ), LIST( LIMIT ), INDEX( LIMIT )
CS    REAL             X     ( NMAX  ), BU  ( NMAX  ), BL   ( NMAX  )
CS    REAL             FACTOR, FBEST , LENGTH, S( DIM * ( DIM + 1 ) )
CS    REAL             SCALE , TOL   , WORK( DIM ) 
CS    REAL             WORK1( -3 : DIM + 1 ), WORK2( -3 : DIM + 1 )
CD    DOUBLE PRECISION X     ( NMAX ), BU   ( NMAX  ), BL   ( NMAX  )
CD    DOUBLE PRECISION FACTOR, FBEST , LENGTH, S( DIM * ( DIM + 1 ) )
CD    DOUBLE PRECISION SCALE , TOL   , WORK( DIM )
CD    DOUBLE PRECISION WORK1( -3 : DIM + 1 ), WORK2( -3 : DIM + 1 )
      CHARACTER*10     XNAMES( NMAX ), PNAME
      INTEGER          I
      LOGICAL          BOUNDS
      REAL             CPU( 2 ), CALLS( 4 )
CS    REAL             BIGINF
CD    DOUBLE PRECISION BIGINF
CS    PARAMETER      ( BIGINF = 9.0E+19 )
CD    PARAMETER      ( BIGINF = 9.0D+19 )
      EXTERNAL         UFN
C
C  Open the Spec file for the method.
C
      OPEN( INDR , FILE = 'PDS.SPC', FORM = 'FORMATTED',
     *      STATUS = 'OLD' )
      REWIND INDR
C
C  Read input Spec data.
C
C   TOL    = the stopping tolerance for the step size
C   MAXITR = the maximum number of iterations allowed
C   TYPE     specifies the type of initial simplex
C            0   User provides the initial simplex
C            1   Automatic generation of a right-angled simplex
C            2   Automatic generation of a regular simplex
C            3   Automatic generation of a scaled right-angled simplex
C   SCALE    If the simplex is automatically generated, scale should
C            contain both the base length and the orientation of the edges
C   DEBUG    should be set to 0, 1, 2, 3 or 4, it controls the amount of 
C            printing.
C            0   no debugging output
C            1   display the iteration count, the best vertex 
C                and its function value
C            2   include the simplex and flag whether or
C                not strict decrease was obtained
C            3   include all vertices constructed and their 
C                finction values
C            4   Include the points that define the search scheme
C   SSS      number of points used in the search strategy
C
C  Set up algorithmic input data.
C
      READ ( INDR, 1000 ) TOL, MAXITR, TYPE, SCALE, DEBUG, SSS
C
C  Close SPEC file
C
      CLOSE( INDR )
C
C  Open the relevant file.
C
      OPEN( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', 
     *      STATUS = 'OLD' )
      REWIND INPUT
C
C  Set up SIF data.
C
      CALL USETUP( INPUT, IOUT, N, X, BL, BU, NMAX )
      CLOSE(INPUT)
C
C  Obtain variable names.
C
      CALL UNAMES( N, PNAME, XNAMES )
C
C  Set up algorithmic input data.
C
      BOUNDS  = .FALSE.
      DO 10 I = 1, N
         IF ( BL( I ) .GT. - BIGINF .OR. BU( I ) .LT. BIGINF )
     *      BOUNDS = .TRUE.
   10 CONTINUE
      IF ( BOUNDS ) WRITE( IOUT, 2030 )
      DO 20 I = 1, N
         S(I) = X(I)
   20 CONTINUE
C
C  Form the scheme
C
      OPEN( SCH, FILE = 'SCHEME', STATUS = 'UNKNOWN', 
     *      FORM = 'UNFORMATTED')
      CALL SEARCH( N, SCH, LIMIT, SCHEME, INDEX, LIST, UNIQUE, IFACT, 
     *             ERROR )
      CLOSE( SCH )
C
      OPEN( RES, FILE = 'RESULT', STATUS = 'UNKNOWN',
     *      FORM = 'FORMATTED')
      FACTOR = FLOAT( IFACT )
      IF ( ERROR .EQ. 0 ) THEN
         WRITE( RES, 100 ) N
         WRITE( RES, 200 ) UNIQUE
         WRITE( RES, 300 ) IFACT
      ELSE
         WRITE( RES, 400 )
      ENDIF
C
C  Re-open the search scheme file.
C
      OPEN( SCH, FILE = 'SCHEME', STATUS = 'OLD',
     *      FORM = 'UNFORMATTED')
      REWIND( SCH )
C
C  Read in the search scheme and determine the "shrink" factor (which
C  depends on the size of the search scheme that has been specified).
C
      CALL GETSS( N, SCH, SSS, SCHEME, FACTOR, RESIZE, ERROR )
C
C  Close the file for the search scheme.
C
      CLOSE( SCH )
      IF ( ERROR .EQ. 0 ) THEN     
C
C  Call the optimizer if no error occured so far
C
         CALL PDS( N, IOUT, TYPE , SCALE , DEBUG , TOL   , MAXITR, 
     *             SSS    , UFN  , FACTOR, SCHEME, RESIZE, S     ,
     *             INDEX  , FBEST, LENGTH, CNT   , WORK  , WORK1 ,
     *             WORK2 )
C
C  Write the results to a file and on the standard output.
C
         WRITE ( IOUT, 2010 )
         DO 30 I = 1, N
            WRITE( IOUT, 2020 ) XNAMES(I), S(I)
   30    CONTINUE
         CALL RESULT( N, CNT, S, FBEST, INDEX, RES )
      ELSE IF ( ERROR .EQ. 1 ) THEN
         WRITE( RES, 500 )
      ELSE IF ( ERROR .EQ. 2 ) THEN
         WRITE( RES, 600 )
      ENDIF
      CLOSE( RES )
C
C  Write results on the standard output
C
      CALL UREPRT( CALLS, CPU )
      WRITE ( IOUT, 2000 ) PNAME, N, CALLS(1), ERROR, FBEST, 
     *                     CPU(1), CPU(2) 
      STOP
C
C  Non-executable statements.
C
 100  FORMAT( 'Successfully completed a search strategy for problems '
     *        ,'of dimension', I6 )
 200  FORMAT( 'The total number of unique points available is', I26 )
 300  FORMAT( 'The factor needed to restore these points to their ',
     *        'real values is', I7 )
 400  FORMAT( 'Returned without a completed search strategy because', /
     *      , 'of internal stack overflow in the QUICKSORT routines.', /
     *      , 'Check the documentation for further details.')
 500  FORMAT( //, ' Search scheme was of the wrong dimension.', /
     *      , ' Exited without calling PDS.' // )
 600  FORMAT( //, ' Insufficient number of points in scheme.', /
     *      , ' Exited without calling PDS.' // )
 1000 FORMAT( D10.3, /, I10, /, I10, /, D10.3, /, I10, /, I10 )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    , ' Code used               :  PDS', /
     *    , ' Problem                 :  ', A10,  /
     *    , ' # variables             =      ', I10 /
     *    , ' # objective functions   =        ', F8.2 /
     *    , ' Exit code               =      ', I10 /
     *    , ' Final f                 = ', E15.7 /
     *    , ' Set up time             =      ', 0P, F10.2, ' seconds' /
     *    , ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *      66('*') / )
 2010 FORMAT( /, '                 X  ' )
 2020 FORMAT(  A10, 1P, D12.4 )
 2030 FORMAT(  /, ' ** Warning from PDSMA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
C
C  End of PDSMA
C
      END
