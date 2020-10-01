C     ( Last modified on 12 Sepc 2004 at 09:35:12 )
C  Correction: 12/Sep/2004: undeclared integer variable declared
      PROGRAM          PRXMA
C
C  PRAXIS test driver for problems derived from SIF files.
C
C  Ph. Toint, for CGT Productions.
C  January 1996.
C
      INTEGER          N, NL, NF, LP, JPRINT, NMAX, ILLCIN, KTM, I
      INTEGER          NFMAX, JRANCH, IOUT, INPUT
      INTEGER          INSPEC, NMX
      PARAMETER      ( NMAX =  20 )
      PARAMETER      ( IOUT  = 6 )
      PARAMETER      ( INPUT = 55, INSPEC = 46 )
      CHARACTER * 10   XNAMES( NMAX ), PNAME
      REAL             CPU( 2 ), CALLS( 4 )
CS    REAL             OBJF, DMACHR
CD    DOUBLE PRECISION OBJF, DMACHR
      EXTERNAL         OBJF, DMACHR

C                                                                               
CS    REAL             V( NMAX, NMAX ), X( NMAX ), D( NMAX ),Q0( NMAX )
CD    DOUBLE PRECISION V( NMAX, NMAX ), X( NMAX ), D( NMAX ),Q0( NMAX )
CS    REAL             Q1( NMAX ), DMIN, EPSMCH, FX, H, QD0, QD1, QF1
CD    DOUBLE PRECISION Q1( NMAX ), DMIN, EPSMCH, FX, H, QD0, QD1, QF1
CS    REAL             SMALL, T, XLDT, XM2, XM4, DSEED, SCBD 
CD    DOUBLE PRECISION SMALL, T, XLDT, XM2, XM4, DSEED, SCBD 
CS    REAL             ONE
CD    DOUBLE PRECISION ONE
CS    PARAMETER      ( ONE = 1.0E0 )
CD    PARAMETER      ( ONE = 1.0D0 )
C                                                                               
      COMMON / CPRAX / V, X, D, Q0, Q1, DMIN, EPSMCH, FX, H, QD0, QD1,
     *                 QF1, SMALL, T, XLDT, XM2, XM4, DSEED, SCBD, N,         
     *                 NL, NF, LP, JPRINT, NMX, ILLCIN, KTM, NFMAX,
     *                 JRANCH
C                                                                               
C  Maximum dimension
C
      NMX = NMAX
C                                                                               
C  LP is the logical unit number for printed output.                            
C                                                                               
      LP = IOUT 
C                                                                               
C  H is an estimate of the distance from the initial point                      
C  to the solution.                                                             
C                                                                               
      H = ONE
C                                                                               
C  EPSMCH is the smallest floating point (real or double precision)             
C  number which, when added to one, gives a result greater than one.            
C                                                                               
      EPSMCH = DMACHR( 1 )
C                                                                               
C  JRANCH = 1 to use BRENT's random,                                            
C  JRANCH = 2 to use function DRANDM.                                           
C                                                                               
      JRANCH = 1                                                                
      CALL RANINI(4.0D0)                                                        
C                                                                               
C  DSEED is an initial seed for DRANDM,                                         
C  a subroutine that generates pseudorandom numbers                             
C  uniformly distributed on (0,1).                                              
C                                                                               
      DSEED = 1234567.0D0
C
C  Open the Spec file for the method.
C
      OPEN ( INSPEC, FILE = 'PRAXIS.SPC', FORM = 'FORMATTED', 
     *       STATUS = 'OLD' )
      REWIND INSPEC
C                                                                               
C
C  Read input Spec data.
C
C     NFMAX :  the maximum number of function calls
C     T     :  the stopping tolerance
C     SCBD  : the upper bound on the scale factors
C     ILLCIN: the "ill-conditioning" flag
C     KTM   :  the maximum number of iterations without improvement
C     JPRINT: the printing specifier
C
      READ ( INSPEC, 1000 ) NFMAX, T, SCBD, ILLCIN, KTM, JPRINT
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
      CALL USETUP( INPUT, IOUT, N, X, Q0, Q1, NMAX )
C
C  Obtain variable names.
C
      CALL UNAMES( N, PNAME, XNAMES )
C
C  Call the optimizer.
C
      CALL PRAXIS( OBJF )
C
C  Exit
C  Note: unfortunately, PRAXIS does not provide any termination code allowing
C        to distinguish whether the routine thinks it has succeeded or not.
C
      CLOSE( INPUT  )
      CALL UREPRT( CALLS, CPU )
      WRITE ( IOUT, 2010 )
      DO 40 I = 1, N
         WRITE( IOUT, 2020 ) XNAMES( I ), X( I )
   40 CONTINUE
      WRITE ( IOUT, 2000 ) PNAME, N, CALLS(1), FX, CPU(1), CPU(2) 
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( I10, 2( /, D10.3), 3( /, I10 ) )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    , ' Code used               :  PRAXIS', /
     *    , ' Problem                 :  ', A10,  /
     *    , ' # variables             =      ', I10 /
     *    , ' # objective functions   =        ', F8.2 /
     *    , ' Final f                 = ', E15.7 /
     *    , ' Set up time             =      ', 0P, F10.2, ' seconds' /
     *      ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *      66('*') / )
 2010 FORMAT( /, '                 X' )
 2020 FORMAT(  1X, A10, 1P, D12.4 )
      END
C
C
C
CS    REAL             FUNCTION OBJF( X, N )
CD    DOUBLE PRECISION FUNCTION OBJF( X, N )
      INTEGER          N
CS    DOUBLE PRECISION X( N )
CD    DOUBLE PRECISION X( N )
      CALL UFN( N, X, OBJF ) 
      RETURN
      END


