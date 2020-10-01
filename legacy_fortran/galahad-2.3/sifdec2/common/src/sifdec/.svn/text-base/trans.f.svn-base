C  THIS VERSION: 11/08/1995 AT 14:44:21 AM.
C  ** VERSION B **

C     INTEGER IOUT, INPUT, OUTPUT, TEMPRY, IAUTO, LNAMES, LDUMMY
C     PARAMETER ( LNAMES = 10000, LDUMMY = LNAMES )
C     CHARACTER * 10 RNAMES( LNAMES ), ADUMMY( LDUMMY )
C     LOGICAL SINGLE
C     IOUT = 21
C     INPUT = 5
C     OUTPUT = 6
C     TEMPRY = 20
C     SINGLE = .FALSE.
C     IAUTO = 1
C     IAD0 = 2
C     CALL TRANS( IOUT, INPUT, OUTPUT, TEMPRY, SINGLE, IAUTO,
C    *            IAD0, RNAMES, LNAMES, ADUMMY, LDUMMY )
C     STOP
C     END

      SUBROUTINE TRANS( IOUT, INPUT, OUTPUT, TEMPRY, SINGLE, IAUTO,
     *                  IAD0, RNAMES, LNAMES, ADUMMY, LDUMMY )
      INTEGER IOUT, INPUT, OUTPUT, TEMPRY, IAUTO, IAD0, LNAMES, LDUMMY
      CHARACTER * 10 RNAMES( LNAMES ), ADUMMY( LDUMMY )
      LOGICAL SINGLE
C
C  Translate a fortran 77 program so that it is capable of
C  accepting AD01 or AD02 reals instead of ordinary reals
C
C  Nick Gould
C  August 11th, 1995
C
      INTEGER I, I1, I2, II, J, JJ, NREAL, NRE, LRE, MRE, NDUMMY
      LOGICAL INTRI, SUB, FUN, NOFIEL, EBRACK
      LOGICAL FREAL, ENDLIN, STARTR
      LOGICAL NUMBA, CHARA
      CHARACTER * 1 CARD
      CHARACTER * 4 AD0
      CHARACTER * 6 LSP, USP
      CHARACTER * 10 LF, UF, BLANK, FIELD
      CHARACTER * 11 LI, UI
      CHARACTER * 13 LS, US
      CHARACTER * 15 UN, ADTYPE
      CHARACTER * 72 NULINE, OLDLIN, BL72
      CHARACTER * 72 RELINE( 20 )
      CHARACTER * 18 LDP, UDP
      PARAMETER ( USP = ' REAL ', UDP = ' DOUBLE PRECISION ' )
      PARAMETER ( UF  = ' FUNCTION ', US  = ' SUBROUTINE ' )
      PARAMETER ( UI  = ' INTRINSIC ' )
      PARAMETER ( LSP = ' real ', LDP = ' double precision ' ) 
      PARAMETER ( LF  = ' function ', LS  = ' subroutine ' )
      PARAMETER ( LI  = ' intrinsic ' )
      PARAMETER ( UN  = '=AD01_UNDEFINED' )
C
C  Create a blank name
C
      DO 5 I = 1, 10
         BLANK( I : I ) = ' '
    5 CONTINUE
      DO 6 I = 1, 72
         BL72( I : I ) = ' '
    6 CONTINUE
C
C  Determine the type of automatic derivative to be used
C
      IF ( IAUTO .EQ. 1 ) THEN
         IF ( SINGLE ) THEN
            ADTYPE = 'FORWARD_SINGLE '
         ELSE
            ADTYPE = 'FORWARD_DOUBLE '
         END IF
      ELSE
         IF ( SINGLE ) THEN
            ADTYPE = 'BACKWARD_SINGLE'
         ELSE
            ADTYPE = 'BACKWARD_DOUBLE'
         END IF
      END IF
C
C  Determine the AD routine to be used
C
      IF ( IAD0 .EQ. 1 ) THEN
         AD0 = 'AD01'
      ELSE
         AD0 = 'AD02'
      END IF
C
C  Initialize logicals
C
      INTRI = .FALSE.
      SUB   = .FALSE.
      FUN   = .FALSE.
      STARTR = .FALSE.
      REWIND( INPUT )
      REWIND( TEMPRY )
   10 CONTINUE
C
C  Read a new line
C
      READ ( INPUT, 1000, END = 600, ERR = 600 ) NULINE
      IF ( SUB .OR. FUN ) THEN
C
C  Ignore comments
C
         IF ( NULINE( 1: 1 ) .EQ. 'C' .OR. 
     *        NULINE( 1: 1 ) .EQ. 'c' ) GO TO 400
C
C  Is the line a continuation?
C
         IF ( NULINE( 6: 6 ) .NE. ' ' ) THEN
            IF ( CARD .EQ. 'I' ) GO TO 10
            II = 7
            IF ( CARD .EQ. 'B' .OR. CARD .EQ. 'E' ) 
     *         CALL GETDUM( NULINE, II, ADUMMY, LDUMMY, NDUMMY, BLANK )
C
C  Find what kind of line this is. Find its first nonzero character
C  Find the start of the name
C
         ELSE
            DO 20 I = 7, 72
               IF ( NULINE( I : I ) .NE. ' ' ) GO TO 30
   20       CONTINUE
            GO TO 10
   30       CONTINUE
            CARD = ' '
            IF ( NULINE( I : I + 6 ) .EQ. 'INTEGER' ) GO TO 400
            IF ( NULINE( I : I + 6 ) .EQ. 'COMPLEX' ) GO TO 400
            IF ( NULINE( I : I + 6 ) .EQ. 'LOGICAL' ) GO TO 400
            IF ( NULINE( I : I + 8 ) .EQ. 'CHARACTER' ) GO TO 400
            IF ( NULINE( I : I + 8 ) .EQ. 'DIMENSION' ) GO TO 400
            IF ( NULINE( I : I + 7 ) .EQ. 'IMPLICIT' ) GO TO 400
            IF ( NULINE( I : I + 10 ) .EQ. 'EQUIVALENCE' ) GO TO 400
            IF ( NULINE( I : I + 7 ) .EQ. 'EXTERNAL' ) THEN
                 CARD = 'E'
                 CALL GETDUM( NULINE, I + 8, ADUMMY, LDUMMY, 
     *                        NDUMMY, BLANK )
                 GO TO 400
            END IF
            IF ( SINGLE ) THEN
               IF ( NULINE( I : I + 15 ) .EQ. 
     *                 'DOUBLE PRECISION' ) GO TO 400
               IF ( NULINE( I : I + 3 ) .EQ. 'REAL' ) THEN
                  CARD = 'R'
                  II = I + 4
                  GO TO 200
               END IF
            ELSE
               IF ( NULINE( I : I + 3 ) .EQ. 'REAL' ) GO TO 400
               IF ( NULINE( I : I + 15 ) .EQ. 
     *                 'DOUBLE PRECISION' ) THEN
                  CARD = 'R'
                  II = I + 16
                  GO TO 200
               END IF
            END IF
            IF ( NULINE( I : I + 3 ) .EQ. 'SAVE' ) THEN
               WRITE( IOUT, 2100 ) AD0
               CARD = 'S'
               II = I + 4
               GO TO 200
            END IF
            IF ( NULINE( I : I + 8 ) .EQ. 'INTRINSIC' ) THEN
               CARD = 'I'
               GO TO 10
            END IF
            IF ( NULINE( I : I + 5 ) .EQ. 'COMMON' ) THEN
               WRITE( IOUT, 2110 ) AD0
               CARD = 'C'
               II = I + 6
               GO TO 200
            END IF
            IF ( NULINE( I : I + 8 ) .EQ. 'PARAMETER' ) THEN
               CARD = 'P'
               II = I + 9
               GO TO 200
            END IF
            IF ( NULINE( I : I + 3 ) .EQ. 'DATA' ) THEN
               CARD = 'D'
               II = I + 4
               GO TO 200
            END IF
C
C  The body of the procedure has been found. Complete the
C  introduction
C 
            REWIND( TEMPRY )
            CARD = 'B'
   60       CONTINUE
            READ ( TEMPRY, 1000, END = 190, ERR = 190 ) OLDLIN
C
C  Write out comments
C
            IF ( OLDLIN( 1: 1 ) .EQ. 'C' .OR. 
     *           OLDLIN( 1: 1 ) .EQ. 'c' ) GO TO 180
C
C  Search for cards defining appropriate real values
C
            IF ( OLDLIN( 6: 6 ) .NE. ' ' ) THEN
               IF ( CARD .NE. 'R' ) GO TO 180
               II = 7
            ELSE
            IF ( CARD .EQ. 'B' ) THEN
               CARD = 'F'
               GO TO 180
            END IF
            IF ( CARD .EQ. 'F' ) WRITE( OUTPUT, 2020 ) AD0, ADTYPE
C
C  Write out the previous set of real values
C       
               IF ( STARTR ) THEN
                  IF ( NRE .GT. 0 ) THEN
                     DO 62 I = 1, MRE - 1 
                        CALL OUTLIN( RELINE( I ), 72, OUTPUT )
   62                CONTINUE   
                     CALL OUTLIN( RELINE( MRE ), LRE, OUTPUT )
                  END IF
                  STARTR = .FALSE.
               END IF
               DO 70 I = 7, 72
                  IF ( OLDLIN( I : I ) .NE. ' ' ) GO TO 80
   70          CONTINUE
               GO TO 60
   80          CONTINUE
               CARD = ' '
               IF ( OLDLIN( I : I + 3 ) .EQ. 'SAVE' ) GO TO 180
               IF ( OLDLIN( I : I + 3 ) .EQ. 'DATA' ) GO TO 180
               IF ( OLDLIN( I : I + 5 ) .EQ. 'COMMON' ) GO TO 180
               IF ( OLDLIN( I : I + 6 ) .EQ. 'INTEGER' ) GO TO 180
               IF ( OLDLIN( I : I + 6 ) .EQ. 'COMPLEX' ) GO TO 180
               IF ( OLDLIN( I : I + 6 ) .EQ. 'LOGICAL' ) GO TO 180
               IF ( OLDLIN( I : I + 7 ) .EQ. 'EXTERNAL' ) GO TO 180
               IF ( OLDLIN( I : I + 7 ) .EQ. 'IMPLICIT' ) GO TO 180
               IF ( OLDLIN( I : I + 8 ) .EQ. 'CHARACTER' ) GO TO 180
               IF ( OLDLIN( I : I + 8 ) .EQ. 'DIMENSION' ) GO TO 180
               IF ( OLDLIN( I : I + 8 ) .EQ. 'PARAMETER' ) GO TO 180
               IF ( OLDLIN( I : I + 10 ) .EQ. 'EQUIVALENCE' ) GO TO 180
               IF ( SINGLE ) THEN
                  IF ( OLDLIN( I : I + 15 ) .EQ. 
     *                 'DOUBLE PRECISION' ) GO TO 180
                  II = I + 4
                  RELINE( 1 ) = BL72
                  RELINE( 1 )( 1 : 11 ) = '      REAL '
                  LRE = 11
                  MRE = 1
               ELSE
                  IF ( OLDLIN( I : I + 3 ) .EQ. 'REAL' ) GO TO 180
                  II = I + 16
                  RELINE( 1 ) = BL72
                  RELINE( 1 )( 1 : 23 ) = '      DOUBLE PRECISION '
                  MRE = 1
                  LRE = 23
               END IF
               NRE = 0
               CARD = 'R'
               STARTR = .TRUE.
            END IF
  110       CONTINUE
            CALL GETFLD( II, I1, I2, ENDLIN, OLDLIN, BLANK, 
     *                   FIELD, NOFIEL, EBRACK, .FALSE. )
            IF ( NOFIEL ) GO TO 60
            J = II - I1
            DO 120 I = 1, NREAL
               IF ( FIELD .EQ. RNAMES( I ) ) THEN
C
C  The parameter will be of type AD01_REAL or AD02_REAL
C
                  DO 115 J = 1, NDUMMY
                     IF ( FIELD .EQ. ADUMMY( J ) ) THEN
                        WRITE( OUTPUT, 2060 ) AD0,
     *                     ( OLDLIN( JJ : JJ ), JJ = I1, II - 1 )
                        GO TO 130
                     END IF
  115             CONTINUE   
                  IF ( IAD0 .EQ. 1 ) THEN
                     WRITE( OUTPUT, 2060 ) AD0,
     *                  ( OLDLIN( JJ : JJ ), JJ = I1, II - 1 ),
     *                  ( UN( JJ : JJ ), JJ = 1, 15 )
                  ELSE
                     WRITE( OUTPUT, 2060 ) AD0,
     *                  ( OLDLIN( JJ : JJ ), JJ = I1, II - 1 )
                  END IF
                  GO TO 130
               END IF
  120       CONTINUE   
C
C  The parameter will be of type REAL
C
            IF ( NRE .GT. 0 ) THEN
               IF ( LRE + 1 .GT. 72 ) THEN
                  MRE = MRE + 1
                  RELINE( MRE ) = BL72
                  RELINE( MRE )( 1 : 8 ) = '     * ,'
                  LRE = 8
               ELSE
                  RELINE( MRE )( LRE + 1 : LRE + 1 ) = ','
                  LRE = LRE + 1
               END IF
            END IF                   
            IF ( LRE + J + 1 .GT. 72 ) THEN
               MRE = MRE + 1
               RELINE( MRE ) = BL72
               RELINE( MRE )( 1 : 7 ) = '     * '
               LRE = 7
            END IF
            RELINE( MRE )( LRE + 1 : LRE + J ) = OLDLIN( I1 : II - 1 )
            LRE = LRE + J
            NRE = NRE + 1
  130       CONTINUE   
            IF ( ENDLIN ) GO TO 60
C
C  Find the next parameter
C
            GO TO 110
  180       CONTINUE
C
C  Output the current line
C
            CALL OUTLIN( OLDLIN, 72, OUTPUT )
            GO TO 60
  190       CONTINUE
C
C  Write out any unfinished set of real values
C       
            IF ( STARTR ) THEN
               IF ( NRE .GT. 0 ) THEN
                  DO 192 I = 1, MRE - 1
                     CALL OUTLIN( RELINE( I ), 72, OUTPUT )
  192             CONTINUE   
                  CALL OUTLIN( RELINE( MRE ), LRE, OUTPUT )
               END IF
               STARTR = .FALSE.
            END IF
C
C  The introduction is complete
C
            IF ( SUB ) THEN
               SUB = .FALSE.
            ELSE
               FUN = .FALSE.
               IF ( FREAL .AND. IAD0 .EQ. 1 ) 
     *            WRITE( OUTPUT, 2030 ) RNAMES( 1 )( 1 : 6 )
               IF ( FREAL .AND. IAD0 .EQ. 2 ) 
     *            WRITE( OUTPUT, 2040 ) AD0, RNAMES( 1 )( 1 : 6 )
            END IF
            REWIND( TEMPRY )
            GO TO 410
         END IF
C
C  Find all variables mentioned on the current string
C
  200    CONTINUE
C
C  Add the variable to the list
C
         IF ( CARD .EQ. 'R' ) THEN
  210       CONTINUE
            CALL GETFLD( II, I1, I2, ENDLIN, NULINE, BLANK, 
     *                   FIELD, NOFIEL, EBRACK, .FALSE. )
            IF ( NOFIEL ) GO TO 400
            DO 220 I = 1, NREAL
               IF ( FIELD .EQ. RNAMES( I ) ) GO TO 230
  220       CONTINUE   
            NREAL = NREAL + 1
            RNAMES( NREAL ) = FIELD
  230       CONTINUE   
            IF ( ENDLIN ) GO TO 400
            GO TO 210
         END IF
C
C  Remove the variable from the list
C
         IF ( CARD .EQ. 'C' .OR. CARD .EQ. 'D' .OR. 
     *        CARD .EQ. 'P' .OR. CARD .EQ. 'S' ) THEN
  310       CONTINUE
            CALL GETFLD( II, I1, I2, ENDLIN, NULINE, BLANK, 
     *                   FIELD, NOFIEL, EBRACK, .FALSE. )
            IF ( NOFIEL ) GO TO 400
            DO 320 I = 1, NREAL
               IF ( FIELD .EQ. RNAMES( I ) ) THEN
                  RNAMES( I ) = RNAMES( NREAL )
                  NREAL = NREAL - 1
                  GO TO 330
               END IF
  320       CONTINUE   
  330       CONTINUE   
            IF ( ENDLIN ) GO TO 400
C
C  For parameter statements, skip the segments after the "="
C
            IF ( CARD .EQ. 'P' ) THEN
               DO 340 I = II, 72
                  IF ( NULINE( I : I ) .EQ. ',' .OR.
     *                 NULINE( I : I ) .EQ. ')' ) THEN
                     II = I + 1
                     GO TO 310
                  END IF
  340          CONTINUE
            END IF
            GO TO 310
         END IF
  400    CONTINUE
         WRITE( TEMPRY, 2000 ) NULINE
         GO TO 10
      END IF
  410 CONTINUE
      IF ( .NOT. ( SUB .OR. FUN ) ) THEN
C
C  Ignore comments
C
         IF ( NULINE( 1: 1 ) .EQ. 'C' .OR. 
     *        NULINE( 1: 1 ) .EQ. 'c' ) GO TO 500
C
C  Remove lines mentioning intrinsic functions
C
         IF ( INTRI ) THEN
            IF ( NULINE( 6: 6 ) .NE. ' ' ) GO TO 10
            INTRI = .FALSE.
         END IF
         DO 420 I = 1, 62
            IF ( NULINE( I : I + 10 ) .EQ. UI .OR.
     *           NULINE( I : I + 10 ) .EQ. LI ) THEN
              INTRI = .TRUE.
              GO TO 10
            END IF         
  420    CONTINUE
C
C  Hunt for the start of a SUBROUTINE
C
         CARD = ' '
         DO 430 I = 1, 60
            IF ( NULINE( I : I + 11 ) .EQ. US .OR.
     *           NULINE( I : I + 11 ) .EQ. LS ) THEN
               II = I + 12
               SUB = .TRUE.
               CARD = 'B'
               NREAL = 0
               NDUMMY = 0
               CALL GETDUM( NULINE, II, ADUMMY, LDUMMY, NDUMMY, BLANK )
               WRITE( TEMPRY, 2000 ) NULINE
               GO TO 10
            END IF         
  430    CONTINUE   
C
C  Hunt for the start of a FUNCTION
C
         DO 440 I = 1, 63
            IF ( NULINE( I : I + 9 ) .EQ. UF .OR.
     *           NULINE( I : I + 9 ) .EQ. LF ) THEN
               II = I + 10
               FUN = .TRUE.
               CARD = 'B'
C
C  Find what kind of function this is
C
               FREAL = .FALSE.
               IF ( SINGLE ) THEN
C
C  Hunt for the string ' REAL '
C
                  DO 431 J = 1, I
                     IF ( NULINE( J : J + 5 ) .EQ. USP .OR.
     *                    NULINE( J : J + 5 ) .EQ. LSP ) THEN
                        FREAL = .TRUE.
                        NULINE( J : J + 5 ) = '      '
                        NULINE( 6 : 6 ) = '*'
                        WRITE( TEMPRY, 2010 ) AD0
                        GO TO 433
                     END IF
  431             CONTINUE   
C
C  Hunt for the string ' DOUBLE PRECISION '
C
               ELSE
                  DO 432 J = 1, I
                     IF ( NULINE( J : J + 17 ) .EQ. UDP .OR.
     *                    NULINE( J : J + 17 ) .EQ. LDP ) THEN
                        IF ( IAD0 .EQ. 1 ) THEN
                           NULINE( J : J + 17 ) = 
     *                        ' TYPE (' // AD0 // '_REAL) '
                        ELSE
                           NULINE( J : J + 17 ) = 
     *                        '                           '
                          AD0 = 'AD02'
                        END IF
                        FREAL = .TRUE.
                        GO TO 433
                     END IF
  432             CONTINUE
               END IF
               WRITE( TEMPRY, 2000 ) NULINE
               NDUMMY = 0
               CALL GETDUM( NULINE, II, ADUMMY, LDUMMY, NDUMMY, BLANK )
               GO TO 10
C
C  The function will be of type AD01_REAL or AD02_REAL. Find its name
C
  433          CONTINUE
               WRITE( TEMPRY, 2000 ) NULINE
C
C  Find the start of the name
C
               DO 434 J = II, 72
                  IF ( NULINE( J : J ) .NE. ' ' ) GO TO 435
  434          CONTINUE
C
C  No name has been found so far. Read the next card
C   
               READ ( INPUT, 1000, END = 600, ERR = 600 ) NULINE
               II = 7
               GO TO 433
C
C  Find the end of the name
C
  435          CONTINUE   
               DO 436 JJ = J + 1, MIN( 72, J + 5 ) 
                  IF ( .NOT. ( CHARA( NULINE( JJ : JJ ) ) .OR.
     *                         NUMBA( NULINE( JJ : JJ ) ) ) ) GO TO 437
  436          CONTINUE
               JJ = MIN( 72, J + 6 )
  437          CONTINUE
               NREAL = 1
               RNAMES( NREAL ) = BLANK
               RNAMES( NREAL )( 1 : JJ - J ) = NULINE( J : JJ - 1 )
               NDUMMY = 0
               CALL GETDUM( NULINE, JJ, ADUMMY, LDUMMY, NDUMMY, BLANK )
               GO TO 10
            END IF         
  440    CONTINUE
C
C  Hunt for the start of a subroutine
C
  500    CONTINUE
         CALL OUTLIN( NULINE, 72, OUTPUT )
         GO TO 10
      END IF
  600 CONTINUE
      RETURN
C
C  Non-executable statements
C
 1000 FORMAT( A72 )
 2000 FORMAT( A72 )
 2010 FORMAT( '      TYPE ( ', A4, '_REAL )' )
 2020 FORMAT( '      USE HSL_', A4, '_', A15 )
C2030 FORMAT( '      CALL AD01_UNDEFINE(', A6, ')' )
 2030 FORMAT( '      ', A6, ' = AD01_UNDEFINED' )
 2040 FORMAT( '      TYPE (', A4, '_REAL) :: ', A6 )
 2060 FORMAT( '      TYPE (', A4, '_REAL) :: ', 46( A1, : ), /,
     *        ( '     *', 66( A1, : ) ) )
 2100 FORMAT( ' ** Warning: a user-supplied external procedure', 
     *        ' SAVEs parameters.', /,
     *        '    This is not allowed by the automatic',
     *        ' differentiation package ', A4, '.' )
 2110 FORMAT( ' ** Warning: a user-supplied external procedure', 
     *        ' saves parameters via common.', /,
     *        '    This is not allowed by the automatic',
     *        ' differentiation package ', A4, '.' )
C
C  End of TRANS
C
      END 

      LOGICAL FUNCTION CHARA( C )
      CHARACTER * 1 C
      INTEGER I
      CHARACTER * 1 LCHARS( 26 )
      CHARACTER * 1 UCHARS( 26 )
      DATA LCHARS / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 
     *              'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 
     *              'u', 'v', 'w', 'x', 'y', 'z' /
      DATA UCHARS / 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 
     *              'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 
     *              'U', 'V', 'W', 'X', 'Y', 'Z' /
      DO 10 I = 1, 26
         IF ( C .EQ. LCHARS( I ) .OR. C .EQ. UCHARS( I ) ) THEN
            CHARA = .TRUE.
            RETURN
         END IF
   10 CONTINUE
      CHARA = .FALSE.
      RETURN
      END

      LOGICAL FUNCTION NUMBA( C )
      CHARACTER * 1 C
      INTEGER I
      CHARACTER * 1 CHARS( 10 )
      DATA CHARS / '0', '1', '2', '3', '4', '5', 
     *             '6', '7', '8', '9'/
      DO 10 I = 1, 10
         IF ( C .EQ. CHARS( I ) ) THEN
            NUMBA = .TRUE.
            RETURN
         END IF
   10 CONTINUE
      NUMBA = .FALSE.
      RETURN
      END

      SUBROUTINE UPPER( C, N )
      INTEGER N
      CHARACTER * ( * ) C
      INTEGER I, J
      CHARACTER * 1 LCHARS( 26 )
      CHARACTER * 1 UCHARS( 26 )
      DATA LCHARS / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 
     *              'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 
     *              'u', 'v', 'w', 'x', 'y', 'z' /
      DATA UCHARS / 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 
     *              'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 
     *              'U', 'V', 'W', 'X', 'Y', 'Z' /
      DO 20 J = 1, N
         DO 10 I = 1, 26
            IF ( C( J : J ) .EQ. LCHARS( I ) ) THEN
               C( J : J ) = UCHARS( I )
               GO TO 20
            END IF
   10    CONTINUE
   20 CONTINUE   
      RETURN
      END

      SUBROUTINE LOWER( C, N )
      INTEGER N
      CHARACTER * ( * ) C
      INTEGER I, J
      CHARACTER * 1 LCHARS( 26 )
      CHARACTER * 1 UCHARS( 26 )
      DATA LCHARS / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 
     *              'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 
     *              'u', 'v', 'w', 'x', 'y', 'z' /
      DATA UCHARS / 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 
     *              'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 
     *              'U', 'V', 'W', 'X', 'Y', 'Z' /
      DO 20 J = 1, N
         DO 10 I = 1, 26
            IF ( C( J : J ) .EQ. UCHARS( I ) ) THEN
               C( J : J ) = LCHARS( I )
               GO TO 20
            END IF
   10    CONTINUE
   20 CONTINUE   
      RETURN
      END

      SUBROUTINE GETFLD( II, I1, I2, ENDLIN, NULINE, BLANK, 
     *                   FIELD, NOFIEL, EBRACK, IGNORB )
      INTEGER II, I1, I2, J
      LOGICAL ENDLIN, NOFIEL, EBRACK, IGNORB
      CHARACTER * 72 NULINE
      CHARACTER * 10 BLANK, FIELD
C
C  Find the first character string in NULINE beyond position II
C
      INTEGER I
      LOGICAL ARRAYS
      LOGICAL NUMBA, CHARA
      EBRACK = .TRUE.
      ENDLIN = .FALSE.
      DO 20 I = II, 72
         IF ( CHARA( NULINE( I : I ) ) ) THEN
            NOFIEL = .FALSE.
            I1 = I
            GO TO 30
         END IF
   20 CONTINUE   
      NOFIEL = .TRUE.
      RETURN
C
C  Next, find its last character
C
   30 CONTINUE   
      DO 40 I = I1 + 1, MIN( 72, I1 + 6 )
         IF ( .NOT. CHARA( NULINE( I : I ) ) .AND.
     *        .NOT. NUMBA( NULINE( I : I ) ) ) THEN
            I2 = I - 1
            II = I
            IF ( IGNORB ) GO TO 70
            GO TO 50
         END IF
   40 CONTINUE   
      I2 = MIN( 72, I1 + 6 )
      ENDLIN = .TRUE.
C
C  Last, check to see if it is an array 
C
   50 CONTINUE
      ARRAYS = .FALSE.
      DO 60 I = I2 + 1, 72
         IF ( .NOT. ARRAYS ) THEN
            IF ( NULINE( I : I ) .NE. ' ' ) THEN
               IF ( NULINE( I : I ) .NE. '(' ) THEN
                  GO TO 70
               ELSE
                  ARRAYS = .TRUE.
               EBRACK = .FALSE.
               END IF
            END IF
         ELSE
            IF ( NULINE( I : I ) .EQ. ')' ) THEN
               EBRACK = .TRUE.
               II = I + 1
               GO TO 70
            END IF
         END IF
   60 CONTINUE   
   70 CONTINUE
C     WRITE( 6, * ) ' LINE, I1, I2, II ', LINE, I1, I2, II
      J = I2 - I1 + 1
      FIELD = BLANK
      FIELD( 1 : J ) = NULINE( I1 : I2 )
      CALL UPPER( FIELD, 10 )
      RETURN
      END


      SUBROUTINE GETDUM( NULINE, II, ADUMMY, LDUMMY, NDUMMY, BLANK )
      INTEGER II, NDUMMY, LDUMMY
      CHARACTER * 10 BLANK
      CHARACTER * 72 NULINE
      CHARACTER * 10 ADUMMY( LDUMMY )
C
C  Determine the list of variables beyond column II on the line NULINE
C
      INTEGER I, JJ, J1, J2
      CHARACTER * 10 FIELD
      LOGICAL NOFIEL, EBRACK, ENDLIN
      JJ = II
   10 CONTINUE
      CALL GETFLD( JJ, J1, J2, ENDLIN, NULINE, BLANK, 
     *             FIELD, NOFIEL, EBRACK, .TRUE. )
      IF ( NOFIEL ) RETURN
      DO 20 I = 1, NDUMMY
         IF ( FIELD .EQ. ADUMMY( I ) ) GO TO 30
   20 CONTINUE   
      NDUMMY = NDUMMY + 1
      ADUMMY( NDUMMY ) = FIELD
   30 CONTINUE   
      IF ( ENDLIN ) RETURN
      GO TO 10
      END

      SUBROUTINE OUTLIN( NULINE, IEND, OUTPUT )
      INTEGER IEND, OUTPUT
      CHARACTER * 72 NULINE
C
C  Write out the current line, cutting off trailing blanks
C
      INTEGER I, I2
      DO 10 I2 = MIN( IEND, 72 ), 1, - 1
         IF ( NULINE( I2 : I2 ) .NE. ' ' ) GO TO 20
   10 CONTINUE
      I2 = 1
   20 CONTINUE
      WRITE( OUTPUT, 2000 ) ( NULINE( I : I ), I = 1, I2 )
      RETURN
 2000 FORMAT( 72( A1 : ) )
      END
