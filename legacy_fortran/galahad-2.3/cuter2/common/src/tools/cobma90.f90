!     ( Last modified on 23 Dec 2000 at 22:01:38 )
PROGRAM COBMA
!
!  COBYLA test driver for problems derived from SIF files.
!
!  A. R. Conn and Ph. Toint (based upon Nick Gould's vf13ma.f)
!  January 1995.
!
!  Fortran 90/95 version, D. Orban, November 2006
!
Use CUTEr_precis
Use CUTEr_interfaces

 INTEGER ::         NMAX, MMAX, M, N, MAXFUN, MCON, LW, LIW
 INTEGER ::         IPRINT, I, MGEQ, INDR, INPUT, IOUT, NFIX
! PARAMETER      ( LW   = NMAX * ( 3 * NMAX + 2 * MMAX + 11 )
!     *                        + 4 * MMAX + 6  )
!      PARAMETER      ( LIW  = NMAX + 1 )
 INTEGER, Dimension(:), Allocatable :: IW !( LIW )

 Real( Kind = wp ), Dimension(:), Allocatable :: X, BL, BU, C, CL, CU, W
 Real( Kind = wp ) :: RHOBEG, RHOEND, F
 Real( Kind = 1.0E+0 ), Dimension(2) :: CPU
 Real( Kind = 1.0E+0 ), Dimension(7) :: CALLS

!S    REAL             X( NMAX ), BL( NMAX ), BU( NMAX ),
!D    DOUBLE PRECISION X( NMAX ), BL( NMAX ), BU( NMAX ),
!     *                 C( MMAX ), CL( MMAX ), CU( MMAX ),
!     *                 W( LW ), RHOBEG, RHOEND, F
!      REAL             CPU( 2 ), CALLS( 7 )

 Logical, Dimension(:), Allocatable :: EQUATN, LINEAR
 Character( Len = 10 ), Dimension(:), Allocatable :: VNAME, CNAME
 Character( Len = 10 ) :: PNAME
 Integer, Parameter :: INPUT = 55, INDR = 46, IOUT = 6

!      LOGICAL          EQUATN( MMAX ), LINEAR( MMAX )
!      CHARACTER * 10   PNAME, VNAME( NMAX ), CNAME( MMAX )

!      PARAMETER      ( INPUT = 55, INDR = 46, IOUT = 6 )
!S    REAL             BIGINF
!D    DOUBLE PRECISION BIGINF
!S    PARAMETER      ( BIGINF = 9.0E+19 )
!D    PARAMETER      ( BIGINF = 9.0D+19 )
      COMMON / FCOBFN /  BL, BU, CL, CU, MGEQ, MCON
!
!  Open the relevant file.
!
      OPEN ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', 
     *       STATUS = 'OLD' )
      REWIND INPUT
!
!  Set up the data structures necessary to hold the problem functions.
!
      CALL CSETUP ( INPUT, IOUT, N, MCON, X, BL, BU, NMAX, EQUATN,
     *              LINEAR, C, CL, CU, MMAX, .TRUE., .FALSE., .FALSE. )
      CLOSE ( INPUT )
!
!  Count the number of general equality constraints and ignore them by
!  shifting the remaining constraints at the beginning of the constraints' list
!
      MGEQ = 0
      DO 20 I = 1, MCON
         IF ( EQUATN( I ) ) MGEQ = MGEQ + 1
   20 CONTINUE
      IF ( MGEQ .GT. 0 ) THEN
         WRITE( 6, 3090 ) MGEQ
         DO 30 I = MCON, MGEQ + 1, -1
            CU( I - MGEQ ) = CU( I )
            CL( I - MGEQ ) = CL( I )
   30    CONTINUE
      END IF
      M = MCON - MGEQ
!
!  If constraints have both lower and upper bounds, they have to be
!  included twice!
!
      DO 40 I = 1, MCON - MGEQ
         IF ( CL( I ) .GT. - BIGINF .AND.
     *        CU( I ) .LT.   BIGINF ) M = M + 1
   40 CONTINUE
!
!  Include any simple bounds.
!
      NFIX = 0
      DO 50 I = 1, N
         IF ( BL( I ) .EQ.  BU( I ) ) THEN
            NFIX = NFIX + 1
         ELSE
            IF ( BL( I ) .GT. - BIGINF ) M = M + 1
            IF ( BU( I ) .LT.   BIGINF ) M = M + 1
         END IF
   50 CONTINUE
      IF ( NFIX .GT. 0 ) WRITE( 6, 3020 ) NFIX
!
      IF ( M .GT. MMAX ) THEN
         IF ( IOUT .GT. 0 )
     *      WRITE( IOUT, 3000 ) 'EQUATN', 'MMAX  ', M - MMAX
         STOP
      END IF
!
!  Open the Spec file for the method.
!
      OPEN( INDR, FILE = 'COBYLA.SPC', FORM = 'FORMATTED', 
     *      STATUS = 'OLD')
      REWIND INDR
!
!  Read input Spec data.
!
!   RHOBEG = the size of the simplex initially
!   RHOEND = the size of the simplex at termination
!   MAXFUN = the maximum number of function calls allowed.
!   IPRINT   should be set to 0, 1, 2 or 3, it controls the amount of printing.
!
!  Set up algorithmic input data.
!
      READ ( INDR, 1000 ) RHOBEG, RHOEND, MAXFUN, IPRINT
      CLOSE ( INDR )
!
!  Evaluate the objective function and constraints.
!
      CALL CALCFC( N, X, F, C )
!
!  Perform the minimization.
!
      CALL COBYLA( N, M, X, RHOBEG, RHOEND, IPRINT, MAXFUN, W, IW)
      CALL CREPRT ( CALLS, CPU )
!
      CALL CNAMES( N, MCON, PNAME, VNAME, CNAME )
      CALL CALCFC( N, X, F, C )
      WRITE( 6, 2110 ) ( I, VNAME( I ), X( I ), BL( I ), BU( I ),
     *                      I = 1, N )
      IF ( MCON .GT. 0 ) WRITE( 6, 2120 ) ( I, CNAME( I ), C( I ),
     *     CL( I ), CU( I ), LINEAR( I ), I = 1, MCON )
      WRITE( 6, 2000 ) PNAME, N, MCON, CALLS(1), CALLS(5), F,
     *     CPU(1), CPU(2)
      STOP

!
!  Non-executable statements
!
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  COBYLA ',  /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # constraints           =      ', I10 /
     *    ,' # objective functions   =        ', F8.2 /
     *    ,' # constraints functions =        ', F8.2 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 1000 FORMAT( D12.4, /, D12.4, /,I6, /, I6 )
 2110 FORMAT( /, ' The variables:', /,
     *        '     I name          value    lower bound upper bound',
     *        /, ( I6, 1X, A10, 1P, 3D12.4 ) )
 2120 FORMAT( /, ' The constraints:', /,
     *        '     I name          value    lower bound upper bound',
     *        ' linear? ',
     *        /, ( I6, 1X, A10, 1P, 3D12.4, 5X, L1 ) )
 3000 FORMAT( /,'  ** Program CSETUP: array length ', A6, ' too small.',
     *        /,'  -- Miminimization abandoned.',
     *        /,'  -- Increase the parameter ', A6, ' by at least ', I8,
     *          ' and restart.'  )
 3020 FORMAT( /,'  ** Warning from COBMA. **',
     *        /,'     In the problem as stated , ', I6,
     *          ' variables are fixed: they are changed to free.' )
 3090 FORMAT( /,'  ** Warning from COBMA. **',
     *        /,'     The problem as stated includes ', I6,
     *          ' equality constraints: they are ignored ' )
!
!  End of COBMA.
!
      END
      SUBROUTINE CALCFC( N, X, F, C )
!
!  Evaluates the objective function value in a format compatible with COBYLA,
!  but using the CUTEr tools.
!
!  A. R. Conn and Ph. Toint
!  January 1995.
!
      INTEGER            NMAX, MMAX
!TOY  PARAMETER        ( NMAX =  10, MMAX =  10 )
!MED  PARAMETER        ( NMAX = 100, MMAX = 100 )
!BIG  PARAMETER        ( NMAX = 500, MMAX = 500 )
!CUS  PARAMETER        ( NMAX = 200, MMAX = 100 )
      INTEGER            N, MGEQ, MCON, I, MT
!S    REAL               F, X( N ), C( MMAX ), BL( NMAX ), BU( NMAX )
!D    DOUBLE PRECISION   F, X( N ), C( MMAX ), BL( NMAX ), BU( NMAX )
!S    REAL               CL( MMAX ), CU( MMAX )
!D    DOUBLE PRECISION   CL( MMAX ), CU( MMAX )
      COMMON /FCOBFN/    BL, BU, CL, CU, MGEQ, MCON
!
!S    REAL               BIGINF
!D    DOUBLE PRECISION   BIGINF
!S    PARAMETER        ( BIGINF = 9.0E+19 )
!D    PARAMETER        ( BIGINF = 9.0D+19 )
!
!  Evaluate the objective function and constraints.
!
      CALL CFN ( N, MCON, X, F, MMAX, C )
!
!  If there are equality constraints, ignore them
!  and shift all the inequality constraint values.
!
      IF ( MGEQ .GT. 0 ) THEN
         DO 10 I = MCON, MGEQ + 1, - 1
            C( I - MGEQ ) = C( I )
   10    CONTINUE
      END IF
!
!  If constraints have both lower and upper bounds, they have to
!  be included twice! Reverse the signs of less-than-or-equal-to
!  constraints.
!
      MT = MCON - MGEQ + 1
      DO 40 I = 1, MCON - MGEQ
         IF ( CL( I ) .GT. - BIGINF .AND.
     *        CU( I ) .LT.   BIGINF ) THEN
            C( I )  = CU( I ) - C( I )
            C( MT ) = C( I ) - CL( I )
            MT      = MT + 1
         ELSE IF ( CL( I ) .GT. - BIGINF ) THEN
            C( I )  = C( I ) - CL( I )
         ELSE IF ( CU( I ) .LT.   BIGINF ) THEN
            C( I )  = CU( I ) - C( I )
         END IF
   40 CONTINUE
!
!  Include any simple bounds, including fixed variables.
!
      DO 50 I = 1, N
         IF ( BL( I ) .NE.  BU( I ) ) THEN
            IF ( BL( I ) .GT. - BIGINF ) THEN
               C( MT ) = X( I ) - BL( I )
               MT      = MT + 1
            END IF
            IF ( BU( I ) .LT.   BIGINF ) THEN
               C( MT ) = BU( I ) - X( I )
               MT      = MT + 1
            END IF
         END IF
   50 CONTINUE
      RETURN
!
!  End of CALCFC.
!
      END
