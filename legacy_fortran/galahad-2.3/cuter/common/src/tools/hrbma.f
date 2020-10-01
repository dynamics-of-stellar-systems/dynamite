C  Original release version: 23 Dec 2000 at 22:01:38 GMT
C  This version: 25 Aug 2004 at 14:00:00 GMT
C
C  25/Aug/2004: duplicate diagonal Hessian entries summed and
C    support for sorted RB format added
C
      PROGRAM HRBMA

C  --------------------------------------------------------------------

C  Write out CUTE data in Harwell-Boeing or Rutherford-Boeing Format

C  Nick Gould
C  October 1996

C  --------------------------------------------------------------------

C  Scalar arguments

      INTEGER INPUT, IIN, IOUT, OUTPUT, OUTRHS, N, M, NMAX, MMAX, MATMAX
      INTEGER NA, NE, NH, NJ, NV, PLAST, NROW, NCOL, NNZ, COLMAX
      INTEGER I, J , NTOTAL, LV
CS    REAL             F, BIGINF, PENLTY, ZERO, ONE
CD    DOUBLE PRECISION F, BIGINF, PENLTY, ZERO, ONE
      LOGICAL HB, RB
      CHARACTER * 1 MATYPE, MAFORM, HRB
      CHARACTER * 10 PNAME
      CHARACTER * 12 PRBOUT
      CHARACTER * 17 PRBRHS
      CHARACTER * 80 LINE1
      CHARACTER * 70 LINE2, LINE3
      CHARACTER * 72 LINE4
      CHARACTER * 42 LINE5

C  Parameters

CTOY  PARAMETER ( NMAX = 100 , MMAX = 100, MATMAX = 10000 )
CMED  PARAMETER ( NMAX = 10000, MMAX = 10000, MATMAX = 100000 )
CBIG  PARAMETER ( NMAX = 100000, MMAX = 100000, MATMAX = 1000000 )
CCUS  PARAMETER ( NMAX = 50000, MMAX = 50000, MATMAX = 100000 )
      PARAMETER ( COLMAX = NMAX + MMAX + 1 )
      PARAMETER ( IIN = 5, IOUT = 6, INPUT = 55 )
      PARAMETER ( OUTPUT = 56, OUTRHS = 57 )
CS    PARAMETER ( ZERO = 0.0E+0, ONE = 1.0E+0 )
CD    PARAMETER ( ZERO = 0.0D+0, ONE = 1.0D+0 )
CS    PARAMETER ( BIGINF = 1.0E+19, PENLTY = 1.0E-1 ) 
CD    PARAMETER ( BIGINF = 1.0D+19, PENLTY = 1.0D-1 ) 

C  Array arguments

      INTEGER ROW( MATMAX ), COL( MATMAX ), IP( COLMAX + 1 ), 
     *        IW( COLMAX + 1 )
CS    REAL             X( NMAX ), BL( NMAX ), BU( NMAX ), CL( MMAX ),
CD    DOUBLE PRECISION X( NMAX ), BL( NMAX ), BU( NMAX ), CL( MMAX ),
     *                 CU( MMAX ), V( MMAX ), C( MMAX ), VAL( MATMAX ),
     *                 B( MMAX ), SLACK( MMAX )
      LOGICAL EQUATN( MMAX ), LINEAR( MMAX )
      CHARACTER * 10 VNAMES( NMAX ), GNAMES( MMAX )

C  Open the relevant file

      OPEN( INPUT , FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *      STATUS = 'OLD'  )
      REWIND INPUT

C  Set up the data structures necessary to hold the group partially
C  separable function

      CALL CSETUP( INPUT, IOUT, N, M, X, BL, BU, NMAX, EQUATN,
     *             LINEAR, V, CL, CU, MMAX, .FALSE., .FALSE., .FALSE. ) 

C  Determine whether Harwell-Boeing or Rutherford-Boeing format is required

   10 CONTINUE
      WRITE( IOUT, 2050 )
      READ( IIN, 1000 ) HRB
      IF ( HRB .EQ. 'h' ) HRB = 'H'
      IF ( HRB .EQ. 'r' ) HRB = 'R'
      IF ( HRB .NE. 'H' .AND. HRB .NE. 'R' ) THEN
         WRITE( IOUT, 2060 ) HRB
         GO TO 10
      END IF
      HB = HRB .EQ. 'H'
      RB = .NOT. HB

C  Determine required matrix type

      IF ( M .GT. 0 ) THEN
   20    CONTINUE
         WRITE( IOUT, 2010 )
         READ( IIN, 1000 ) MATYPE
         IF ( MATYPE .EQ. 'a' ) MATYPE = 'A'
         IF ( MATYPE .EQ. 'h' ) MATYPE = 'H'
         IF ( MATYPE .EQ. 'j' ) MATYPE = 'J'
         IF ( MATYPE .EQ. 't' ) MATYPE = 'T'
         IF ( MATYPE .NE. 'A' .AND. MATYPE .NE. 'H' .AND.
     *        MATYPE .NE. 'J' .AND. MATYPE .NE. 'T' ) THEN
            WRITE( IOUT, 2020 ) MATYPE
            GO TO 20
         END IF
         MAFORM = 'A'
      ELSE
         MATYPE = 'H'
      END IF

      IF ( MATYPE .EQ. 'H' ) THEN

C  Determine required matrix format

   30    CONTINUE
         WRITE( IOUT, 2030 )
         READ( IIN, 1000 ) MAFORM
         IF ( MAFORM .EQ. 'a' ) MAFORM = 'A'
         IF ( MAFORM .EQ. 'e' ) MAFORM = 'E'
         IF ( MAFORM .NE. 'A' .AND. MAFORM .NE. 'E' ) THEN
            WRITE( IOUT, 2040 ) MAFORM
            GO TO 30
         END IF
      END IF

C  Determine the names of the problem, variables and constraints

      CALL CNAMES( N, M, PNAME, VNAMES, GNAMES )
      DO 50 PLAST = 8, 1, - 1
         IF ( PNAME( PLAST : PLAST ) .NE. ' ' ) GO TO 60
   50 CONTINUE
   60 CONTINUE

C  Name output file

      PRBOUT = '            '
      IF ( MATYPE .EQ. 'H' ) THEN
         IF ( MAFORM .EQ. 'A' ) THEN
            PRBOUT = PNAME( 1 : PLAST ) // '.hes'
         ELSE
            PRBOUT = PNAME( 1 : PLAST ) // '.ele'
         END IF
      ELSE IF ( MATYPE .EQ. 'A' ) THEN
         PRBOUT = PNAME( 1 : PLAST ) // '.aug'
      ELSE IF ( MATYPE .EQ. 'J' ) THEN
         PRBOUT = PNAME( 1 : PLAST ) // '.jac'
      ELSE
         PRBOUT = PNAME( 1 : PLAST ) // '.jat'
      END IF

C  Open output file

      OPEN( OUTPUT , FILE = PRBOUT, FORM = 'FORMATTED',
     *      STATUS = 'UNKNOWN'  )
      REWIND OUTPUT

C  For RB format, the RHS is stored separately

      IF ( RB ) THEN
         IF ( MATYPE .EQ. 'H' ) THEN
            IF ( MAFORM .EQ. 'A' ) THEN
               PRBRHS = PNAME( 1 : PLAST ) // '.hes.rhsr'
            ELSE
               PRBRHS = PNAME( 1 : PLAST ) // '.ele.rhsr'
            END IF
         ELSE IF ( MATYPE .EQ. 'A' ) THEN
            PRBRHS = PNAME( 1 : PLAST ) // '.aug.rhsr'
         ELSE IF ( MATYPE .EQ. 'J' ) THEN
            PRBRHS = PNAME( 1 : PLAST ) // '.jac.rhsr'
         ELSE
            PRBRHS = PNAME( 1 : PLAST ) // '.jat.rhsr'
         END IF

C  Open output file for RHS

         OPEN( OUTRHS , FILE = PRBRHS, FORM = 'FORMATTED',
     *         STATUS = 'UNKNOWN'  )
         REWIND OUTRHS

      END IF

C  Move X into the feasible region

      DO 110 I = 1, N
         X( I ) = MIN( BU( I ), MAX( BL( I ), X( I ) ) )
  110 CONTINUE

C  Evaluate the constant terms of the objective and constraint functions

      CALL CFN( N, M, X, F, MMAX, B )

      DO 120 I = 1, M
         B( I ) = - B( I )
  120 CONTINUE

c  Use IW to store positions of diagonal entries. Initialize IW as 0

      DO 130 I = 1, N + M
         IW( I ) = 0
  130 CONTINUE

C  Evaluate the linear terms of the objective and constraint functions
C  in a sparse format

      CALL CSGR( N, M, .FALSE., MMAX, V, X, NJ, MATMAX, VAL, COL, ROW )
      
      DO 210 I = 1, M
         SLACK( I ) = ONE
  210 CONTINUE   

C  Remove zero entries

      NA = 0

      DO 220 I = 1, NJ
         IF ( VAL( I ) .NE. ZERO ) THEN
            IF ( ROW( I ) .GT. 0 ) THEN
               NA = NA + 1
               J = ROW( I )  
               ROW( NA ) = J
               COL( NA ) = COL( I )  
               VAL( NA ) = VAL( I )
               SLACK( J ) = MAX( SLACK( J ), ABS( VAL( NA ) ) )
            ELSE
               C( COL( I ) ) = - VAL( I )
            END IF  
         END IF
  220 CONTINUE

C  Determine the status of each slack variable.
C  Introduce slack variables for inequality constraints

      NTOTAL = N
      DO 230 I = 1, M
         IF ( .NOT. EQUATN( I ) ) THEN
            NTOTAL = NTOTAL + 1
            NA = NA + 1
            ROW( NA ) = I
            COL( NA ) = NTOTAL 
            IF ( ( CL( I ) .EQ. ZERO .AND. CU( I ) .GT. BIGINF ) .OR.
     *           ( CU( I ) .EQ. ZERO .AND. CL( I ) .LT. - BIGINF ) ) 
     *             THEN
               IF ( CL( I ) .EQ. ZERO ) THEN
                  VAL( NA ) = - SLACK( I )
               ELSE
                  VAL( NA ) = SLACK( I )
               END IF
            ELSE
               VAL( NA ) = - SLACK( I )
            END IF
         END IF
  230 CONTINUE

      IF ( MATYPE .EQ. 'H' ) THEN
         NNZ = 0
      ELSE
         NNZ = NA
      END IF

      IF ( MATYPE .EQ. 'A' .OR. MATYPE .EQ. 'H' ) THEN

C  The matrix is to be assembled

         IF ( MAFORM .EQ. 'A' ) THEN

C  Evaluate the Hessian of the Lagrangian function at the initial point

            CALL CSH( N, M, X, LV, V, NH, MATMAX - NA, VAL( NA + 1 ),
     *                ROW( NA + 1 ), COL( NA + 1 ) )

C  Remove zero entries

            DO 310 I = NA + 1, NA + NH
               IF ( VAL( I ) .NE. ZERO ) THEN
                  NNZ = NNZ + 1
                  IF ( ROW( I ) .GE. COL( I ) ) THEN
                     ROW( NNZ ) = ROW( I )  
                     COL( NNZ ) = COL( I )  
                  ELSE
                     ROW( NNZ ) = COL( I )  
                     COL( NNZ ) = ROW( I )  
                  END IF  
                  IF ( ROW( I ) .EQ. COL( I ) ) IW( ROW( I ) ) = NNZ
                  VAL( NNZ ) = VAL( I )
               END IF
  310       CONTINUE

C  Include terms representing penalties

            DO 320 I = 1, N
               IF ( BL( I ) .GT. - BIGINF .OR. BU( I ) .LT. BIGINF ) 
     *              THEN
                  J = IW( I )
                  IF ( J .EQ. 0 ) THEN
                     NNZ = NNZ + 1
                     ROW( NNZ ) = I
                     COL( NNZ ) = I
                     IF ( BL( I ) .GT. - BIGINF .AND. BU( I ) .LT. 
     *                   BIGINF ) THEN
                        VAL( NNZ ) = PENLTY + PENLTY
                     ELSE
                        VAL( NNZ ) = PENLTY
                     END IF
                  ELSE
                     IF ( BL( I ) .GT. - BIGINF .AND. BU( I ) .LT. 
     *                   BIGINF ) THEN
                        VAL( J ) = VAL( J ) + PENLTY + PENLTY
                     ELSE
                        VAL( J ) = VAL( J ) + PENLTY
                     END IF
                  END IF
               END IF 
  320       CONTINUE

            NTOTAL = N
            DO 330 I = 1, M
               IF ( .NOT. EQUATN( I ) ) THEN
                  NTOTAL = NTOTAL + 1
                  NNZ = NNZ + 1
                  ROW( NNZ ) = NTOTAL
                  COL( NNZ ) = NTOTAL 
                  IF ( CL( I ) .GT. - BIGINF .AND. CU( I ) .LT. BIGINF ) 
     *                 THEN
                     VAL( NNZ ) = PENLTY + PENLTY
                  ELSE
                     VAL( NNZ ) = PENLTY
                  END IF
               END IF
  330       CONTINUE

C  The matrix is unassembled

         ELSE

C  Evaluate the Hessian of the Lagrangian function at the initial point

            CALL CEH( N, M, X, LV, V, NE, ROW, MATMAX, COLMAX, IP, 
     *                VAL, MATMAX, COL, .TRUE. )

C  Include terms representing penalties

            NV = COL( NE + 1 ) - 1
            DO 340 I = 1, N
               IF ( BL( I ) .GT. - BIGINF .OR. BU( I ) .LT. BIGINF ) 
     *              THEN
                  NE = NE + 1
                  NV = NV + 1
                  ROW( IP( NE ) ) = I
                  IP( NE + 1 ) = IP( NE ) + 1
                  IF ( BL( I ) .GT. - BIGINF .AND. BU( I ) .LT. BIGINF ) 
     *                 THEN
                     VAL( NV ) = PENLTY + PENLTY
                  ELSE
                     VAL( NV ) = PENLTY
                  END IF
               END IF
  340       CONTINUE
            NNZ = IP( NE + 1 ) - 1

         END IF
      END IF

      N = NTOTAL
      WRITE( 6, 2000 ) PNAME, N, M

      LINE1( 1 : 10 ) = PNAME
      DO 510 I = 11, 80, 10
         LINE1( I : I + 9 ) = '          '
  510 CONTINUE
      LINE1( 73 : 80 ) = PNAME( 1 : 8 )

C  Format header cards

      IF ( MATYPE .EQ. 'A' ) THEN

         NROW = N + M
         NCOL = NROW
         LINE1( 12 : 30 ) = 'augmented system - '
         WRITE( UNIT = LINE1( 31: 38 ), FMT = "( I8 )" ) NROW
         LINE1( 39 : 43 ) = ' rows'

         DO 520 I = 1, NA
            ROW( I ) =  ROW( I ) + N 
  520    CONTINUE

      ELSE IF ( MATYPE .EQ. 'H' ) THEN

         NROW = N
         NCOL = NROW
         LINE1( 12 : 21 ) = 'Hessian - '
         WRITE( UNIT = LINE1( 22: 29 ), FMT = "( I8 )" ) NROW
         LINE1( 30 : 34 ) = ' rows'
         IF ( MAFORM .EQ. 'E' ) LINE1( 35 : 51 ) = ' - element format'

      ELSE IF ( MATYPE .EQ. 'J' ) THEN

         NROW = M
         NCOL = N
         LINE1( 12 : 22 ) = 'Jacobian - '
         WRITE( UNIT = LINE1( 23: 30 ), FMT = "( I8 )" ) NROW
         LINE1( 31 : 35 ) = ' rows'
         WRITE( UNIT = LINE1( 36: 43 ), FMT = "( I8 )" ) NCOL
         LINE1( 44 : 51 ) = ' columns'

      ELSE IF ( MATYPE .EQ. 'T' ) THEN

         NROW = N
         NCOL = M
         LINE1( 12 : 32 ) = 'Jacobian transpose - '
         WRITE( UNIT = LINE1( 33: 40 ), FMT = "( I8 )" ) NROW
         LINE1( 41 : 45 ) = ' rows'
         WRITE( UNIT = LINE1( 46: 53 ), FMT = "( I8 )" ) NCOL
         LINE1( 54 : 61 ) = ' columns'

         DO 530 I = 1, NA
            J = COL( I )
            COL( I ) =  ROW( I )
            ROW( I ) = J
  530    CONTINUE

      END IF

C  Transform from co-ordinate to column format

      IF ( MAFORM .EQ. 'A' )
     *   CALL REORDA( NCOL, NNZ, ROW, COL, VAL, IP, IW )

C  Harwell-Boeing format

      IF ( HB ) THEN

C  Continue formatting header cards

         IF ( MAFORM .EQ. 'A' ) THEN
            I = ( NCOL + 10 ) / 10 + ( NNZ  + 9 ) / 10 + 
     *          ( NNZ  +  2 ) / 3  + ( NROW + 2 ) / 3
         ELSE
            I = ( NCOL + 10 ) / 10 + ( NNZ  + 9 ) / 10 + 
     *          ( NV   +  2 ) / 3  + ( NROW + 2 ) / 3
         END IF
         WRITE( UNIT = LINE2( 1: 14 ), FMT = "( I14 )" ) I
         I = ( NCOL + 10 )/10 
         WRITE( UNIT = LINE2( 15: 28 ), FMT = "( I14 )" ) I
         I = ( NNZ + 9 ) / 10
         WRITE( UNIT = LINE2( 29: 42 ), FMT = "( I14 )" ) I
         IF ( MAFORM .EQ. 'A' ) THEN
            I = ( NNZ + 2 ) / 3
         ELSE
            I = ( NV + 2 ) / 3
         END IF
         WRITE( UNIT = LINE2( 43: 56 ), FMT = "( I14 )" ) I
         I = ( NROW + 2 ) / 3
         WRITE( UNIT = LINE2( 57: 70 ), FMT = "( I14 )" ) I
   
         LINE3( 1: 1 ) = 'R'
         IF ( MATYPE .EQ. 'A' .OR. MATYPE .EQ. 'H' ) THEN
            LINE3( 2: 2 ) = 'S'
         ELSE
            LINE3( 2: 2 ) = 'U'
         END IF
         IF ( MAFORM .EQ. 'A' ) THEN
            LINE3( 3: 3 ) = 'A'
         ELSE
            LINE3( 3: 3 ) = 'E'
         END IF
         WRITE( UNIT = LINE3( 4: 14 ), FMT = "( 11X )" ) 
         WRITE( UNIT = LINE3( 15: 28 ), FMT = "( I14 )" ) NROW
         WRITE( UNIT = LINE3( 29 : 42 ), FMT = "( I14 )" ) NCOL
         WRITE( UNIT = LINE3( 43 : 56 ), FMT = "( I14 )" ) NNZ
         IF ( MAFORM .EQ. 'A' ) THEN
            I = 0
         ELSE
            I = NV
         END IF
         WRITE( UNIT = LINE3( 57 : 70 ), FMT = "( I14 )" ) I
   
         WRITE( UNIT = LINE4( 1: 16 ), FMT = "( '(10I8)          ' )" ) 
         WRITE( UNIT = LINE4( 17: 32 ), FMT = "( '(10I8)          ' )" ) 
CS       WRITE( UNIT = LINE4( 33: 52 ), FMT = "( '(1P, 3E24.16)   ' )" ) 
CD       WRITE( UNIT = LINE4( 33: 52 ), FMT = "( '(1P, 3D24.16)   ' )" ) 
CS       WRITE( UNIT = LINE4( 53: 72 ), FMT = "( '(1P, 3E24.16)   ' )" ) 
CD       WRITE( UNIT = LINE4( 53: 72 ), FMT = "( '(1P, 3D24.16)   ' )" ) 
   
         LINE5( 1: 3 ) = 'F  '
         WRITE( UNIT = LINE5( 4: 14 ), FMT = "( 11X )" ) 
         I = 1
         WRITE( UNIT = LINE5( 15: 28 ), FMT = "( I14 )" ) I
         WRITE( UNIT = LINE5( 29: 42 ), FMT = "( I14 )" ) NROW

C  Write out the header

         WRITE( OUTPUT, "( A80, /, A70, /, A70, /, A72, /, A42 )" )
     *          LINE1, LINE2, LINE3, LINE4, LINE5

C  Rutherford-Boeing format

      ELSE

C  Continue formatting header cards

         IF ( MAFORM .EQ. 'A' ) THEN
            I = ( NCOL + 10 ) / 10 + ( NNZ  + 9 ) / 10 + 
     *          ( NNZ  +  2 ) / 3
         ELSE
            I = ( NCOL + 10 ) / 10 + ( NNZ  + 9 ) / 10 + 
     *          ( NV   +  2 ) / 3
         END IF
         WRITE( UNIT = LINE2( 1: 14 ), FMT = "( 1X, I13 )" ) I
         I = ( NCOL + 10 )/10 
         WRITE( UNIT = LINE2( 15: 28 ), FMT = "( 1X, I13 )" ) I
         I = ( NNZ + 9 ) / 10
         WRITE( UNIT = LINE2( 29: 42 ), FMT = "( 1X, I13 )" ) I
         IF ( MAFORM .EQ. 'A' ) THEN
            I = ( NNZ + 2 ) / 3
         ELSE
            I = ( NV + 2 ) / 3
         END IF
         WRITE( UNIT = LINE2( 43: 56 ), FMT = "( 1X, I13 )" ) I
   
         LINE3( 1: 1 ) = 'r'
         IF ( MATYPE .EQ. 'A' .OR. MATYPE .EQ. 'H' ) THEN
            LINE3( 2: 2 ) = 's'
         ELSE
            LINE3( 2: 2 ) = 'u'
         END IF
         IF ( MAFORM .EQ. 'A' ) THEN
            LINE3( 3: 3 ) = 'a'
         ELSE
            LINE3( 3: 3 ) = 'e'
         END IF
         WRITE( UNIT = LINE3( 4: 14 ), FMT = "( 11X )" ) 
         WRITE( UNIT = LINE3( 15: 28 ), FMT = "( 1X, I13 )" ) NROW
         WRITE( UNIT = LINE3( 29 : 42 ), FMT = "( 1X, I13 )" ) NCOL
         WRITE( UNIT = LINE3( 43 : 56 ), FMT = "( 1X, I13 )" ) NNZ
         IF ( MAFORM .EQ. 'A' ) THEN
            I = 0
         ELSE
            I = NV
         END IF
         WRITE( UNIT = LINE3( 57 : 70 ), FMT = "( 1X, I13 )" ) I
   
         WRITE( UNIT = LINE4( 1: 16 ),  FMT = "( '(10I8)          ' )" ) 
         WRITE( UNIT = LINE4( 17: 32 ), FMT = "( '(10I8)          ' )" ) 
CS       WRITE( UNIT = LINE4( 33: 52 ), FMT = "( '(1P, 3E25.16)   ' )" ) 
CD       WRITE( UNIT = LINE4( 33: 52 ), FMT = "( '(1P, 3D25.16)   ' )" ) 
   
C  Write out the header

         WRITE( OUTPUT, "( A80, /, A56, /, A70, /, A52 )" )
     *          LINE1, LINE2( 1 : 56 ), LINE3, LINE4( 1 : 52 )

C  Do the same for the RHS file

         LINE1( 64 : 72 ) = ' - RHS - '

         WRITE( UNIT = LINE2( 1: 5 ), FMT = "( 'rhsrd' )" )
         IF ( MATYPE .EQ. 'A' ) THEN
              WRITE( UNIT = LINE2( 6: 14 ), FMT = "( ' aug sys ' )" )
         ELSE IF ( MATYPE .EQ. 'A' ) THEN
              WRITE( UNIT = LINE2( 6: 14 ), FMT = "( ' hessian ' )" )
         ELSE IF ( MATYPE .EQ. 'J' ) THEN
              WRITE( UNIT = LINE2( 6: 14 ), FMT = "( ' jacobian' )" )
         ELSE IF ( MATYPE .EQ. 'T' ) THEN
              WRITE( UNIT = LINE2( 6: 14 ), FMT = "( ' jactrans' )" )
         END IF

         WRITE( UNIT = LINE2( 15: 16 ), FMT = "( ' r' )" )
         WRITE( UNIT = LINE2( 17: 30 ), FMT = "( 1X, I13 )" ) NROW
         WRITE( UNIT = LINE2( 31: 44 ), FMT = "( 1X, I13 )" ) 1
         WRITE( UNIT = LINE2( 45: 58 ), FMT = "( 1X, I13 )" ) NROW

         WRITE( UNIT = LINE3( 1: 20 ), FMT = "( '(1P, 3D25.16)', 7X )" ) 
         WRITE( UNIT = LINE3( 21: 40 ), FMT = "( 20X )" )
         WRITE( UNIT = LINE3( 41: 60 ), FMT = "( 20X )" )

C  Write out the RHS header

         WRITE( OUTRHS, "( A80, /, A58, /, A60 )" )
     *          LINE1, LINE2( 1 : 58 ), LINE3( 1 : 60 )

      END IF

C  Write out the desired matrix

      IF ( HB ) THEN
         WRITE( OUTPUT, "( (10I8) )" ) ( IP( I ), I = 1, NCOL + 1 )
         WRITE( OUTPUT, "( (10I8) )" ) ( ROW( I ), I = 1, NNZ )
         IF ( MAFORM .EQ. 'A' ) THEN
CS         WRITE( OUTPUT, "( (1P, 3E24.16) )" ) ( VAL( I ), I = 1, NNZ )
CD         WRITE( OUTPUT, "( (1P, 3D24.16) )" ) ( VAL( I ), I = 1, NNZ )
         ELSE
CS         WRITE( OUTPUT, "( (1P, 3E24.16) )" ) ( VAL( I ), I = 1, NV )
CD         WRITE( OUTPUT, "( (1P, 3D24.16) )" ) ( VAL( I ), I = 1, NV )
         END IF
      ELSE
         WRITE( OUTPUT, "( (10I8) )" ) ( IP( I ), I = 1, NCOL + 1 )
         WRITE( OUTPUT, "( (10I8) )" ) ( ROW( I ), I = 1, NNZ )
         IF ( MAFORM .EQ. 'A' ) THEN
CS         WRITE( OUTPUT, "( (1P, 3E25.16) )" ) ( VAL( I ), I = 1, NNZ )
CD         WRITE( OUTPUT, "( (1P, 3D25.16) )" ) ( VAL( I ), I = 1, NNZ )
         ELSE
CS         WRITE( OUTPUT, "( (1P, 3E25.16) )" ) ( VAL( I ), I = 1, NV )
CD         WRITE( OUTPUT, "( (1P, 3D25.16) )" ) ( VAL( I ), I = 1, NV )
         END IF
      END IF

      IF ( HB ) THEN
         IF ( MATYPE .EQ. 'A' ) THEN
CS          WRITE( OUTPUT, "( (1P, 3E24.16) )" ) ( C( I ), I = 1, N ),
CD          WRITE( OUTPUT, "( (1P, 3D24.16) )" ) ( C( I ), I = 1, N ),
     *                                           ( B( I ), I = 1, M )
         ELSE IF ( MATYPE .EQ. 'J' ) THEN
CS          WRITE( OUTPUT, "( (1P, 3E24.16) )" ) ( B( I ), I = 1, M )
CD          WRITE( OUTPUT, "( (1P, 3D24.16) )" ) ( B( I ), I = 1, M )
         ELSE
CS          WRITE( OUTPUT, "( (1P, 3E24.16) )" ) ( C( I ), I = 1, N )
CD          WRITE( OUTPUT, "( (1P, 3D24.16) )" ) ( C( I ), I = 1, N )
         END IF
      ELSE
         IF ( MATYPE .EQ. 'A' ) THEN
CS          WRITE( OUTRHS, "( (1P, 3E25.16) )" ) ( C( I ), I = 1, N ),
CD          WRITE( OUTRHS, "( (1P, 3D25.16) )" ) ( C( I ), I = 1, N ),
     *                                           ( B( I ), I = 1, M )
         ELSE IF ( MATYPE .EQ. 'J' ) THEN
CS          WRITE( OUTRHS, "( (1P, 3E25.16) )" ) ( B( I ), I = 1, M )
CD          WRITE( OUTRHS, "( (1P, 3D25.16) )" ) ( B( I ), I = 1, M )
         ELSE
CS          WRITE( OUTRHS, "( (1P, 3E25.16) )" ) ( C( I ), I = 1, N )
CD          WRITE( OUTRHS, "( (1P, 3D25.16) )" ) ( C( I ), I = 1, N )
         END IF
      END IF

      CLOSE( OUTPUT )
      IF ( RB ) CLOSE( OUTRHS )

      STOP
C
C  Non-executable statements
C
 1000 FORMAT( A1 )
 2000 FORMAT( /, ' Problem name: ', A8,
     *        /, ' Number of variables     = ', I8, 
     *        /, ' Number of equations     = ', I8 )
 2010 FORMAT( /, ' Please state required matrix type. ',
     *        /, ' A (augmented system) or H (Hessian) or',
     *        /, ' J (Jacobian) or T (Transpose of Jacobian) : ' )
 2020 FORMAT( /, ' Matrix type ', A1, 
     *           ' not recognized. Please try again ' )
 2030 FORMAT( /, ' Please state required matrix format. ',
     *        /, ' A (assembled) or E (elemental): ' )
 2040 FORMAT( /, ' Matrix format ', A1, 
     *           ' not recognized. Please try again ' )
 2050 FORMAT( /, ' Please state whether you require'
     *        /, ' Harwell-Boeing or Rutherford-Boeing format. ', /, 
     *           ' H (Harwell-Boeing) or R (Rutherford-Boeing) : ' )
 2060 FORMAT( /, ' Your response ', A1, 
     *           ' not recognized. Please try again ' )

C  End of GETHRB

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

C  Now sort the entries in each column so that their row
C  indices increase

      DO 110 I = 1, NC
         LOCAT = IP( I )
         IC = IP( I + 1 ) - LOCAT
         J = MAX( IC, 1 )
         CALL SORTIN( IC, J, IRN( LOCAT ), A( LOCAT ) )
  110 CONTINUE
      RETURN

C  End of REORDA

      END

      SUBROUTINE SORTIN( N, LIND, IND, VAL )
      INTEGER N, LIND
      INTEGER IND( LIND )
CS    REAL             VAL( LIND )
CD    DOUBLE PRECISION VAL( LIND )

C  Sort N numbers into ascending order. Yes, we should use
C  quicksort but ...well, we are hoping that N won't be too big

      INTEGER I, J, CURMIN, INDMIN
CS    REAL             VALMIN
CD    DOUBLE PRECISION VALMIN

C Find the I-th smallest value

      DO 100 I = 1, N
         CURMIN = I
         INDMIN = IND( CURMIN )
         DO 90 J = I + 1, N
            IF ( IND( J ) .LT. INDMIN ) THEN
               CURMIN = J
               INDMIN = IND( CURMIN )
            END IF
   90    CONTINUE

C If the I-th smallest value is out of place, swap it to poiition CURMIN

         IF ( CURMIN .NE. I ) THEN
            VALMIN = VAL( CURMIN )
            IND( CURMIN ) = IND( I )
            VAL( CURMIN ) = VAL( I )
            IND( I ) = INDMIN
            VAL( I ) = VALMIN
         END IF
  100 CONTINUE
      RETURN

C  End of SORTIN

      END
