C ** Correction report.
C  THIS VERSION: 04/04/2002 AT 09:30:00 AM.
C     ( Last modified on 15 Mar 2001 at 22:28:00 )
C ** Correction report.
C ** Correction 1. 19/07/93: 13 lines added **
C ** Correction 2. 19/07/93: 2 lines added **
C ** Correction 3. 12/01/94: 1 line corrected **
C ** Correction 4. 12/01/94: 1 line corrected **
C ** Correction 5. 11/08/95: 1 line corrected **
C ** Correction 6. 15/08/95: 2 lines corrected **
C ** Correction 7. 26/02/01: 2 dummy arguments removed **
C ** Correction 8. 04/04/02: 1 line corrected
C ** End of Correction report.
      SUBROUTINE MAKEGR( INPUT , IOUT  , IOUTGR, INFORM, NGRTYP,
     *                   NGRMAX, NLMAX , NINMAX, PNAME , ANAMES,
     *                   RENAME, INNAME, LONAME, MINAME, EXNAME, GTYPES,
     *                   LDEFND, GPNAME, IGPA  , NGPMAX, DEBUG , LENGTH,
     *                   ITABLE, KEY   , INLIST, SINGLE, NULINE, GOTLIN,
     *                   IPRINT )
      INTEGER            INPUT , IOUT  , IOUTGR, INFORM, LENGTH
      INTEGER            NLMAX , NGRTYP, NINMAX, NGPMAX, NGRMAX
      INTEGER            IPRINT
      LOGICAL            GOTLIN, DEBUG , SINGLE
      CHARACTER * 8      PNAME
      CHARACTER * 160    NULINE
      INTEGER            ITABLE( LENGTH ), IGPA( NGRMAX )
      INTEGER            INLIST( LENGTH )
      LOGICAL            LDEFND( NLMAX  )
      CHARACTER * 12     KEY   ( LENGTH )
      CHARACTER * 10     RENAME( NINMAX ), INNAME( NINMAX )
      CHARACTER * 10     GPNAME( NGPMAX ), EXNAME( NINMAX )
      CHARACTER * 10     LONAME( NINMAX ), GTYPES( NGRMAX )
      CHARACTER * 10     ANAMES( NGRMAX ), MINAME( NINMAX )
C
C  MAKE A GROUP FUNCTION EVALUATION SUBROUTINE
C  -------------------------------------------
C  FROM A GPS GROUP FUNCTION DATA FILE.
C  ------------------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
C  -------------------------------------------------------------------
C
C  FUNCTION INDICATOR CARDS.
C  -------------------------
C
C  DEFINITION   PURPOSE.
C  ----------   --------
C  GROUPS       PROBLEM NAME.
C  TEMPORARIES  NAMES OF ADDITIONAL PARAMETERS USED IN FUNCTION DEFS.
C  GLOBALS      GENERAL PARAMETER ASSIGNMENTS.
C  INDIVIDUALS  SET FUNCTION AND DERIVATIVE VALUES FOR EACH GROUP-TYPE.
C  ENDATA       END OF INPUT DATA.
C
C  DATA CARD DESCRIPTION.
C  ----------------------
C
C  SEE 'A PROPOSAL FOR A STANDARD DATA INPUT FORMAT FOR LARGE-SCALE
C       NONLINEAR PROGRAMMING PROBLEMS', SECTION 4,
C       A. R. CONN, N. I. M. GOULD AND PH. L. TOINT, 
C       REPORT CS-89-61, DEPT OF COMPUTER SCIENCE, U. OF WATERLOO,
C       WATERLOO, ONTARIO, N2L3G1, CANADA.
C
C  -------------------------------------------------------------------
C  RETURNS WITH NEGATIVE VALUES OF INFORM INDICATE THAT INSUFFICIENT
C  ARRAY SPACE HAS BEEN ALLOWED, AS FOLLOWS:
C
C    INFORM = - 1  WHEN LENGTH NOT LARGE ENOUGH
C    INFORM = - 2  WHEN MAX( NINNAM, NRENAM, NLONAM, NENAM, NMINAM )
C                  .GT. NINMAX
C
      INTEGER          I, IFIELD, IFREE, ITYPE, INTYPE, IVAR, K1, K2
      INTEGER          NINNAM, NLOOP, NRENAM, NMINAM, NPNAME
      INTEGER          NEXNAM, NINCRS, LINENO, IIRES, K, ILINES, NLINES
      INTEGER          MBLANK, MFIXED, MFREE, MNAME, MTEMP, MGLOB
      INTEGER          MINDIV, MENDAT, MAXNUL, MXRECL, NLONAM
      LOGICAL          DEFNAM, ENDPAR, ENDGEN, FIRSTG
      LOGICAL          SETF, SETG, SETH, STARTP, FIXED
      LOGICAL          ENDF
      PARAMETER        ( IIRES = 20 )
      PARAMETER        ( NINCRS = 2 )
      CHARACTER * 2    FIELD1
      CHARACTER * 6    INCRSE( NINCRS )
      CHARACTER * 8    FIELD2, FIELD3, FIELDI( IIRES )
      CHARACTER * 12   FIELD
      CHARACTER * 41   FIELD7
      EXTERNAL         HASHB , HASHC
C
C  PARAMETER DEFINITIONS.
C
      PARAMETER        ( MXRECL = 160 )
      CHARACTER * 160  BLNKLN
      PARAMETER        ( MBLANK =  1, MFIXED =  2, MFREE  = 3  )
      PARAMETER        ( MNAME  =  4, MTEMP  =  5 )
      PARAMETER        ( MGLOB  =  6, MINDIV =  7, MENDAT =  8 )
      INTEGER          LENIND( MENDAT )
      CHARACTER * 12   INDIC8( MENDAT ), HEADER
      PARAMETER        ( MAXNUL = 20 )
      CHARACTER * 65   NULINA( MAXNUL )
C
C  DATA DECLARATIONS.
C
      DATA INCRSE / 'LENGTH', 'NINMAX' /
      DATA INDIC8( MBLANK ) / '            ' /, LENIND( MBLANK ) / 0  / 
      DATA INDIC8( MFIXED ) / 'FIXED FORMAT' /, LENIND( MFIXED ) / 12 /
      DATA INDIC8( MFREE  ) / 'FREE FORMAT ' /, LENIND( MFREE  ) / 11 /
      DATA INDIC8( MNAME  ) / 'GROUPS      ' /, LENIND( MNAME  ) / 6  /
      DATA INDIC8( MTEMP  ) / 'TEMPORARIES ' /, LENIND( MTEMP  ) / 11 / 
      DATA INDIC8( MGLOB  ) / 'GLOBALS     ' /, LENIND( MGLOB  ) / 7  /
      DATA INDIC8( MINDIV ) / 'INDIVIDUALS ' /, LENIND( MINDIV ) / 11 / 
      DATA INDIC8( MENDAT ) / 'ENDATA      ' /, LENIND( MENDAT ) / 6  /
      DATA FIELDI(  1 ) / 'GROUP   ' /,  FIELDI(  2 ) / 'GVALUE  ' /
      DATA FIELDI(  3 ) / 'LGVALU  ' /,  FIELDI(  4 ) / 'FVALUE  ' /
      DATA FIELDI(  5 ) / 'NCALCG  ' /,  FIELDI(  6 ) / 'ITYPEG  ' /
      DATA FIELDI(  7 ) / 'ICALCG  ' /,  FIELDI(  8 ) / 'DERIVS  ' /
      DATA FIELDI(  9 ) / 'IGRTYP  ' /,  FIELDI( 10 ) / 'IGROUP  ' /
      DATA FIELDI( 11 ) / 'GPVALU  ' /,  FIELDI( 12 ) / 'ISTGPA  ' /
      DATA FIELDI( 13 ) / 'IPSTRT  ' /,  FIELDI( 14 ) / 'JCALCG  ' /
      DATA FIELDI( 15 ) / 'LTYPEG  ' /,  FIELDI( 16 ) / 'LSTGPA  ' /
      DATA FIELDI( 17 ) / 'LCALCG  ' /,  FIELDI( 18 ) / 'LFVALU  ' /
      DATA FIELDI( 19 ) / 'LGPVLU  ' /,  FIELDI( 20 ) / 'IGSTAT  ' /

      IF ( IOUT .GT. 0 ) WRITE( IOUT, 2900 )
C
C  SET INITIAL VALUES FOR INTEGER VARIABLES.
C
      NINNAM = 0
      NRENAM = 0
      NLONAM = 0
      NEXNAM = 0
      NMINAM = 0
      LINENO = 0
      INTYPE = 1
      ILINES = 0
      NLINES = 0
C
C  SET INITIAL VALUES FOR LOGICAL VARIABLES.
C
      DEFNAM = .FALSE.
      ENDPAR = .FALSE.
      SETH   = .FALSE.
      STARTP = .FALSE.
      ENDGEN = .FALSE.
      FIRSTG = .TRUE.
      FIXED  = .TRUE.
C
C  FIND WHICH GROUP-TYPES ARE NONTRIVIAL.
C
      DO 20 ITYPE        = 1, NGRTYP
         LDEFND( ITYPE ) = .FALSE.
   20 CONTINUE
C
C  INSERT THE LIST OF GROUP-TYPE ARGUMENTS INTO THE DICTIONARY.
C
      DO 30 ITYPE = 1, NGRTYP
         FIELD = ANAMES( ITYPE ) // 'PG'
         CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               GO TO 700
            END IF
         ELSE
            NRENAM = NRENAM + 1
            IF ( NRENAM .GT. NINMAX ) THEN
               INFORM = - 2
               GO TO 700
            END IF
            RENAME( NRENAM ) = ANAMES( ITYPE )
         END IF
   30 CONTINUE
C
C  INCLUDE THE NAMES OF THE GROUP PARAMETERS USED
C  IN THIS DICTIONARY.
C
      IF ( NGRTYP .GT. 0 ) THEN
         NPNAME  = IGPA( NGRTYP + 1 ) - 1
         DO 40 I = 1, NPNAME
            FIELD = GPNAME( I ) // 'PG'
            CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
            IF ( IFREE .LE. 0 ) THEN
               IF ( IFREE .EQ. 0 ) THEN
                  INFORM = - 1
                  GO TO 700
               END IF
            ELSE
               NRENAM = NRENAM + 1
               IF ( NRENAM .GT. NINMAX ) THEN
                  INFORM = - 2
                  GO TO 700
               END IF
               RENAME( NRENAM ) = GPNAME( I )
            END IF
   40    CONTINUE
      END IF
C
C  SET A BLANK LINE.
C
      DO 50 I = 1, MXRECL
         BLNKLN( I: I ) = ' '
   50 CONTINUE    
C
C  READ NEXT LINE.
C
  100 CONTINUE
      IF ( ILINES + 1 .GT. NLINES ) THEN
C
C  READ NEXT LINE FROM THE INPUT FILE.
C
         LINENO = LINENO + 1
         IF ( FIXED ) THEN
            IF ( GOTLIN ) THEN
               GOTLIN = .FALSE.
            ELSE
               NULINE = BLNKLN
               READ ( INPUT, 1000, END = 590, ERR = 590 ) NULINE
            END IF
            IF ( IOUT .GT. 0 .AND. DEBUG ) WRITE( IOUT, 2990 )
     *           LINENO, NULINE
         ELSE
            IF ( GOTLIN ) THEN
               GOTLIN = .FALSE.
            ELSE
               NULINE = BLNKLN
               READ ( INPUT, 1010, END = 590, ERR = 590 ) NULINE
            END IF
            IF ( IOUT .GT. 0 .AND. DEBUG ) WRITE( IOUT, 2970 )
     *           LINENO, NULINE
C
C  IF THE CARD IS IN FREE FORMAT, TRANSLATE IT INTO FIXED FORMAT.
C
            CALL  FREEFM( NULINE, MXRECL, MENDAT, INDIC8, LENIND,
     *                    NULINA, MAXNUL, NLINES, .FALSE.,
     *                    INFORM, IOUT )
            IF ( INFORM .GT. 0 ) GO TO 800
            IF ( NLINES .GT. 0 ) THEN
C
C  IF THERE ARE NON-BLANK LINES ON THE FREE FORMAT CARD, READ THE FIRST.
C
               ILINES = 1
               NULINE = BLNKLN
               NULINE = NULINA( ILINES )
               IF ( IOUT .GT. 0 .AND. DEBUG ) WRITE( IOUT, 2980 )
     *              LINENO, ILINES, NULINE
            ELSE
C
C  THERE ARE ONLY BLANK LINES ON THE FREE FORMAT CARD.
C
               GO TO 100
            END IF
         END IF
      ELSE
C
C  READ NEXT LINE FROM THE LAST ENCOUNTERED FREE FORMAT CARD.
C
         ILINES = ILINES + 1
         NULINE = BLNKLN
         NULINE = NULINA( ILINES )
         IF ( IOUT .GT. 0 .AND. DEBUG ) WRITE( IOUT, 2980 )
     *        LINENO, ILINES, NULINE
      END IF
C
C  CONSIDER THE HEADER PART OF THE CARD.
C
      HEADER = NULINE( 1: 12 )
C
C  IGNORE BLANK LINES.
C
      IF ( HEADER .EQ. INDIC8( MBLANK ) ) GO TO 100
      IF ( NULINE( 1: 1 ) .NE. ' ' ) THEN
C
C  IGNORE COMMENT CARDS.
C
         IF ( NULINE( 1: 1 ) .EQ. '*' ) GO TO 100
C
C  CHECK IF WE HAVE ENTERED FIXED-FORMAT INPUT.
C
         IF ( HEADER .EQ. INDIC8( MFIXED ) ) THEN
            FIXED = .TRUE.
            GO TO 100
         END IF
C
C  CHECK IF WE HAVE ENTERED FREE-FORMAT INPUT.
C
         IF ( HEADER .EQ. INDIC8( MFREE ) ) THEN
            FIXED = .FALSE.
            GO TO 100
         END IF
C
C  CHECK THAT THE FIRST ENCOUNTERED INDICATOR CARD IS THE GROUPS CARD.
C
         IF ( .NOT. DEFNAM  ) THEN
            IF ( HEADER .NE. INDIC8( MNAME ) ) THEN
               IF ( NGRTYP .GT. 0 ) GO TO 930
               IF ( IOUT .GT. 0 .AND. IPRINT .NE. 0 ) WRITE( IOUT, 2010)
               GOTLIN = .TRUE.
               GO TO 600
            ELSE
C
C  INDICATOR CARD IS GROUPS.
C  -------------------------
C
               IF ( PNAME  .NE. NULINE( 15: 22 ) ) THEN
                  INFORM = 51
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2510 )
                  GO TO 800
               ELSE
                  DEFNAM = .TRUE.
                  GO TO 100
               END IF
            END IF
         END IF
C
C  AN INDICATOR CARD HAS BEEN FOUND.
C
         DO 110 I = INTYPE, MENDAT
            IF ( HEADER .EQ. INDIC8( I ) ) THEN
               INTYPE = I
               GO TO 120
            END IF
  110    CONTINUE
C
C  THE INDICATOR CARD IS NOT RECOGNISED.
C
         INFORM = 2
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2020 )
         GO TO 800
  120    CONTINUE
C
C  THE PARAMETER VALUES HAVE BEEN COMPLETED. WRITE OUT THE
C  FIRST PART OF THE GENERATED SUBROUTINE.
C
         IF ( INTYPE .GE. MGLOB .AND. .NOT. ENDPAR ) THEN
            ENDPAR = .TRUE.
            NLOOP  = NGRTYP + 1
C
C  INSERT THE LIST OF RESERVED INTEGER/REAL/LOGICAL VARIABLES INTO
C  THE DICTIONARY.
C
            DO 130 I = 1, IIRES
               FIELD = FIELDI( I ) // '  PG'
               CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
               IF ( IFREE .LE. 0 ) THEN
                  IF ( IFREE .EQ. 0 ) THEN
                     INFORM = - 1
                     GO TO 700
                  END IF
                  INFORM = 59
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2590 ) FIELDI( I )
                  GO TO 800
               END IF
  130       CONTINUE
C
C  -------- SET UP SUBROUTINE CALL AND RESERVED PARAMETER DECLARATIONS.
C
            IF ( SINGLE ) THEN
               WRITE( IOUTGR, 3001 ) ( FIELDI( I )( 1 : 6 ), I = 1, 4 ),
     *                    FIELDI( 11 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *                    FIELDI(  6 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *                    FIELDI(  7 )( 1 : 6 ), 
     *                  ( FIELDI(  I )( 1 : 6 ), I = 15, 19 ),
     *                    FIELDI(  8 )( 1 : 6 ), FIELDI( 20 )( 1 : 6 ),
     *                    FIELDI(  3 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *                  ( FIELDI(  I )( 1 : 6 ), I = 15, 20 ),
     *                    FIELDI(  8 )( 1 : 6 ), 
     *                    FIELDI(  6 )( 1 : 6 ), FIELDI( 15 )( 1 : 6 ),
     *                    FIELDI( 12 )( 1 : 6 ), FIELDI( 16 )( 1 : 6 ),
     *                    FIELDI(  7 )( 1 : 6 ), FIELDI( 17 )( 1 : 6 ),
     *                    FIELDI(  2 )( 1 : 6 ), FIELDI(  3 )( 1 : 6 ),
     *                    FIELDI(  4 )( 1 : 6 ), FIELDI( 18 )( 1 : 6 ),
     *                    FIELDI( 11 )( 1 : 6 ), FIELDI( 19 )( 1 : 6 ),
     *                    PNAME
            ELSE
               WRITE( IOUTGR, 3000 ) ( FIELDI( I )( 1 : 6 ), I = 1, 4 ),
     *                    FIELDI( 11 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *                    FIELDI(  6 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *                    FIELDI(  7 )( 1 : 6 ), 
     *                  ( FIELDI(  I )( 1 : 6 ), I = 15, 19 ),
     *                    FIELDI(  8 )( 1 : 6 ), FIELDI( 20 )( 1 : 6 ),
     *                    FIELDI(  3 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *                  ( FIELDI(  I )( 1 : 6 ), I = 15, 20 ),
     *                    FIELDI(  8 )( 1 : 6 ), 
     *                    FIELDI(  6 )( 1 : 6 ), FIELDI( 15 )( 1 : 6 ),
     *                    FIELDI( 12 )( 1 : 6 ), FIELDI( 16 )( 1 : 6 ),
     *                    FIELDI(  7 )( 1 : 6 ), FIELDI( 17 )( 1 : 6 ),
     *                    FIELDI(  2 )( 1 : 6 ), FIELDI(  3 )( 1 : 6 ),
     *                    FIELDI(  4 )( 1 : 6 ), FIELDI( 18 )( 1 : 6 ),
     *                    FIELDI( 11 )( 1 : 6 ), FIELDI( 19 )( 1 : 6 ),
     *                    PNAME
            END IF
C ** Correction 8. 04/04/02: 1 line replaced by 4
            IF ( NGRTYP .EQ. 0 ) THEN
              WRITE( IOUTGR, 3009 ) FIELDI( 20 )( 1 : 6 )
              GO TO 910
            END IF
            WRITE( IOUTGR, 3002 ) FIELDI(  9 )( 1 : 6 ),
     *                    FIELDI( 10 )( 1 : 6 ), FIELDI( 13 )( 1 : 6 ),
     *                    FIELDI( 14 )( 1 : 6 )
C
C --------- INSERT INTEGER DECLARATIONS.
C
            IF ( NINNAM .GT. 0 )
     *         WRITE( IOUTGR, 3010 ) ( INNAME( I ), I = 1, NINNAM )
C
C --------- INSERT REAL DECLARATIONS.
C
            IF ( NRENAM .GT. 0 ) THEN
               IF ( SINGLE ) THEN
                  WRITE( IOUTGR, 3019 ) ( RENAME( I ), I = 1, NRENAM )
               ELSE
                  WRITE( IOUTGR, 3020 ) ( RENAME( I ), I = 1, NRENAM )
               END IF
            END IF
C
C --------- INSERT LOGICAL DECLARATIONS.
C
            IF ( NLONAM .GT. 0 )
     *         WRITE( IOUTGR, 3023 ) ( LONAME( I ), I = 1, NLONAM )
C
C --------- INSERT INTRINSIC DECLARATIONS.
C
            IF ( NMINAM .GT. 0 )
     *         WRITE( IOUTGR, 3021 ) ( MINAME( I ), I = 1, NMINAM )
C
C --------- INSERT EXTERNAL DECLARATIONS.
C
            IF ( NEXNAM .GT. 0 )
     *         WRITE( IOUTGR, 3022 ) ( EXNAME( I ), I = 1, NEXNAM )
            WRITE( IOUTGR, 3009 ) FIELDI( 20 )( 1 : 6 )
         END IF
C
C  THE GENERAL PARAMETER ASSIGNMENTS HAVE BEEN COMPLETED.
C  CONTINUE WITH THE CONSTRUCTION OF THE GENERATED SUBROUTINE.
C
         IF ( INTYPE .GE. MINDIV .AND. .NOT. ENDGEN ) THEN
            ENDGEN = .TRUE.
C
C --------- START LOOP OVER GROUPS.
C
            WRITE( IOUTGR, 3050 ) NLOOP,  FIELDI( 14 )( 1 : 6 ),
     *             FIELDI(  5 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *             FIELDI(  7 )( 1 : 6 ), FIELDI( 14 )( 1 : 6 ),
     *             FIELDI(  9 )( 1 : 6 ),
     *             FIELDI(  6 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *             FIELDI(  9 )( 1 : 6 ), NLOOP,
     *             FIELDI( 13 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *             FIELDI( 10 )( 1 : 6 )
            IF ( NGRTYP .GT. 1 ) THEN
               WRITE( IOUTGR, 3051 ) ( I, I = 1, NGRTYP )
               WRITE( IOUTGR, 3052 ) FIELDI(  9 )( 1 : 6 )
            END IF
         END IF
C
C  INDICATOR CARD IS ENDATA.
C  -------------------------
C
         IF ( INTYPE .EQ. MENDAT ) GO TO 900
         GO TO 100
      ELSE
C
C  CHECK THAT THE FIRST NON COMMENT CARD IS THE GROUPS INDICATOR CARD.
C
         IF ( .NOT. DEFNAM  ) THEN
            IF ( NGRTYP .GT. 0 ) GO TO 930
            IF ( IOUT .GT. 0 .AND. IPRINT .NE. 0 ) WRITE( IOUT, 2010 )
            GOTLIN = .TRUE.
            GO TO 600
         END IF
C
C  A DATA CARD HAS BEEN FOUND.
C  READ THE CHARACTER FIELDS 1, 2, 3 AND 7 FROM THE CARD.
C
         FIELD1 = NULINE(  2:  3 )
         FIELD2 = NULINE(  5: 12 )
         FIELD3 = NULINE( 15: 22 )
         FIELD7 = NULINE( 25: 65 )
C ** Correction 1. 19/07/93: 13 lines added **
C
C  CHECK THAT FIELD3 IS BLANK ON 'A', 'F' AND 'G' CARDS.
C
            IF ( FIELD1( 1: 1 ) .EQ. 'A' .OR.
     *           FIELD1( 1: 1 ) .EQ. 'F' .OR. 
     *           FIELD1( 1: 1 ) .EQ. 'G' .OR. 
     *           FIELD1( 1: 1 ) .EQ. 'H' ) THEN
               IF ( ( FIELD1( 1: 1 ) .NE. 'A' .AND. FIELD2 .NE. 
     *              '       ' ) .OR.FIELD3 .NE. '       ' ) THEN
                  INFORM = 73
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2730 )
                  GO TO 800
               END IF
            END IF
C ** Correction 1. 19/07/93: end of correction **
      END IF
C
C  BRANCH ON THE VALUE OF INTYPE.
C
      GO TO ( 100, 100, 100, 100, 290, 300, 400, 900 ), INTYPE
C
C  INDICATOR CARD IS TEMPORARIES.
C  ------------------------------
C
  290 CONTINUE
C
C  CHECK TO SEE IF THE PARAMETER IS INTEGER, REAL, LOGICAL OR A FUNCTION.
C
      IF ( FIELD1 .NE. 'I ' .AND. FIELD1 .NE. 'R ' .AND.
     *     FIELD1 .NE. 'M ' .AND. FIELD1 .NE. 'F ' .AND.
     *     FIELD1 .NE. 'L' ) THEN
         INFORM = 54
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2540 )
         GO TO 800
      END IF
C
C  IF THE PARAMETER IS A FUNCTION, CHECK TO SEE THAT THE NAME HAS
C  NOT ALREADY BEEN USED.
C
      IF ( FIELD1 .EQ. 'F ' ) THEN
         FIELD = FIELD2 // '  GU'
         CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               GO TO 700
            END IF
         ELSE
            NEXNAM = NEXNAM + 1
            IF ( NEXNAM .GT. NINMAX ) THEN
               INFORM = - 2
               GO TO 700
            END IF
            EXNAME( NEXNAM ) = FIELD2
         END IF
      ELSE
C
C  CHECK TO SEE THAT THE PARAMETER NAME HAS NOT ALREADY BEEN USED.
C
         FIELD = FIELD2 // '  PG'
         CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               GO TO 700
            END IF
         ELSE
            IF ( FIELD1 .EQ. 'R ' ) THEN
               NRENAM = NRENAM + 1
               IF ( NRENAM .GT. NINMAX ) THEN
                  INFORM = - 2
                  GO TO 700
               END IF
               RENAME( NRENAM ) = FIELD2
            ELSE
               IF ( FIELD1 .EQ. 'M ' ) THEN
                  NMINAM = NMINAM + 1
                  IF ( NMINAM .GT. NINMAX ) THEN
                     INFORM = - 2
                     GO TO 700
                  END IF
                  MINAME( NMINAM ) = FIELD2
               ELSE
                  IF ( FIELD1 .EQ. 'L ' ) THEN
                     NLONAM = NLONAM + 1
                     IF ( NLONAM .GT. NINMAX ) THEN
                        INFORM = - 2
                        GO TO 700
                     END IF
                     LONAME( NLONAM ) = FIELD2
                  ELSE
                     NINNAM = NINNAM + 1
                     IF ( NINNAM .GT. NINMAX ) THEN
                        INFORM = - 2
                        GO TO 700
                     END IF
                     INNAME( NINNAM ) = FIELD2
                  END IF
               END IF
            END IF
         END IF
      END IF
      GO TO 100
C
C  INDICATOR CARD IS GLOBAL.
C  -------------------------
C
  300 CONTINUE
      IF ( FIELD1 .EQ. 'A ' .OR. FIELD1 .EQ. 'I ' .OR.
     *     FIELD1 .EQ. 'E ' ) THEN
         STARTP = .TRUE.
C
C  START A PARAMETER ASSIGNMENT. CHECK TO SEE THAT THE PARAMETER HAS
C  BEEN DEFINED.
C
         FIELD = FIELD2 // '  PG'
         CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
         IF ( IFIELD .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               GO TO 700
            END IF
            INFORM = 57
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2570 )
            GO TO 800
         END IF
C
C --------- MAKE GENERAL PARAMETER ASSIGNMENTS.
C
         IF ( FIELD1 .EQ. 'A ' ) THEN
            WRITE( IOUTGR, 3030 ) FIELD2( 1 : 6 ), FIELD7
C
C --------- MAKE CONDITIONAL PARAMETER ASSIGNMENTS.
C
         ELSE   
C
C  CHECK THAT THE LOGICAL VARIABLE HAS BEEN DEFINED.
C
            FIELD = FIELD3 // '  PG'
            CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               IF ( IFREE .EQ. 0 ) THEN
                  INFORM = - 1
                  GO TO 700
               END IF
               INFORM = 57
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2570 )
               GO TO 800
            END IF
            IF ( FIELD1 .EQ. 'I ' ) THEN
               WRITE( IOUTGR, 3031 ) FIELD2( 1 : 6 ),
     *                               FIELD3( 1 : 6 ), FIELD7
            ELSE
               WRITE( IOUTGR, 3032 ) FIELD2( 1 : 6 ),
     *                               FIELD3( 1 : 6 ), FIELD7
            END IF   
         END IF
      ELSE
         IF ( FIELD1( 2: 2 ) .EQ. '+' .AND. STARTP ) THEN
C
C --------- CONTINUE A PARAMETER ASSIGNMENT.
C
            WRITE( IOUTGR, 3040 ) FIELD7
         ELSE
            INFORM = 55
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2550 )
            GO TO 800
         END IF
      END IF
      GO TO 100
C
C  INDICATOR CARD IS INDIVIDUALS.
C  ------------------------------
C
  400 CONTINUE
C
C  CHECK IF A NEW GROUP HAS BEEN ENCOUNTERED.
C
      IF ( FIELD1 .EQ. 'T ' ) THEN
         IF ( FIRSTG ) THEN
C
C  CHECK IF THIS IS THE FIRST GROUP-TYPE.
C
            FIRSTG = .FALSE.
         ELSE
C
C  FINISH OF THE PREVIOUS GROUP, IF ANY.
C
            IF ( .NOT. SETH ) THEN
               INFORM = 63
C ** Correction 5. 11/08/95: 1 line corrected **
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2630 )
C ** Correction 5. 11/08/95: end of correction **
               GO TO 800
            END IF
            IF ( .NOT. SETG ) THEN
               INFORM = 62
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2620 )
               GO TO 800
            END IF
C
C ---------- WIND UP F AND G
C
            IF ( SETF ) THEN
               WRITE( IOUTGR, 3190 )
            ELSE
               INFORM = 61
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
               GO TO 800
            END IF
            IF ( ITYPE .LT. NGRTYP ) WRITE( IOUTGR, 3191 ) NLOOP
         END IF
C
C  FIND ITYPE, THE GROUP-TYPE.
C
         FIELD = FIELD2 // '  GT'
         CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
C
C  THE GROUP-TYPE IS UNKNOWN.
C
         IF ( IFIELD .LE. 0 ) THEN
            INFORM = 19
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2190 )
            GO TO 800
         END IF
C
C --------- FIND TYPE OF CURRENT GROUP.
C
         ITYPE = INLIST( IFIELD )
         WRITE( IOUTGR, 3060 ) FIELD2
         IF ( NGRTYP .GT. 1 ) WRITE( IOUTGR, 3061 ) ITYPE
         WRITE( IOUTGR, 3062 ) ANAMES( ITYPE )( 1 : 6 ),
     *             FIELDI(  4 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
C
C --------- SET GROUP PARAMETERS.
C
         K1       = IGPA( ITYPE )
         K2       = IGPA( ITYPE + 1 ) - 1
         DO 435 K = K1, K2
            IVAR  = K - K1 + 1
            WRITE( IOUTGR, 3063 ) GPNAME( K ), FIELDI( 11 )( 1 : 6 ),
     *          FIELDI( 13 )( 1 : 6 ), IVAR
  435    CONTINUE
         IF ( LDEFND( ITYPE ) ) THEN
            INFORM = 64
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2640 )
            GO TO 800
         ELSE
            LDEFND( ITYPE ) = .TRUE.
         END IF
C
C  INITIALIZE LOGICALS WHICH DETERMINE WHETHER THE DATA HAS BEEN
C  INPUT IN THE CORRECT ORDER.
C
         STARTP = .FALSE.
         SETF   = .FALSE.
         SETG   = .FALSE.
         SETH   = .FALSE.
         ENDF   = .TRUE.
      ELSE
         IF ( FIELD1( 1: 1 ) .EQ. 'A' .OR. FIELD1( 1: 1 ) 
     *        .EQ. 'I' .OR. FIELD1( 1: 1 ) .EQ. 'E' ) THEN
            IF ( SETF ) THEN
               IF ( .NOT. ENDF ) THEN
                  WRITE( IOUTGR, 3120 )
                  ENDF = .TRUE.
               END IF
            END IF
C
C  START A PARAMETER ASSIGNMENT. CHECK TO SEE THAT THE PARAMETER HAS
C  BEEN DEFINED.
C
            IF ( FIELD1( 2: 2 ) .EQ. ' ' ) THEN
               STARTP = .TRUE.
               FIELD = FIELD2 // '  PG'
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .LE. 0 ) THEN
                  INFORM = 57
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2570 )
                  GO TO 800
               END IF
C
C --------- MAKE GROUP-SPECIFIC PARAMETER ASSIGNMENTS.
C
               IF ( FIELD1( 1: 1 ) .EQ. 'A' ) THEN
                  IF ( .NOT. SETF ) THEN
                     WRITE( IOUTGR, 3080 ) FIELD2( 1 : 6 ), FIELD7
                  ELSE
                     WRITE( IOUTGR, 3083 ) FIELD2( 1 : 6 ), FIELD7
                  END IF
C
C --------- MAKE CONDITIONAL PARAMETER ASSIGNMENTS.
C
               ELSE   
C
C  CHECK THAT THE LOGICAL VARIABLE HAS BEEN DEFINED.
C
                  FIELD = FIELD3 // '  PG'
                  CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
                  IF ( IFIELD .LE. 0 ) THEN
                     IF ( IFREE .EQ. 0 ) THEN
                        INFORM = - 1
                        GO TO 700
                     END IF
                     INFORM = 58
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2580 )
                     GO TO 800
                  END IF
                  IF ( FIELD1( 1: 1 ) .EQ. 'I' ) THEN
                     IF ( .NOT. SETF ) THEN
                        WRITE( IOUTGR, 3081 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                     ELSE
                        WRITE( IOUTGR, 3084 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                     END IF
                  ELSE
                     IF ( .NOT. SETF ) THEN
                        WRITE( IOUTGR, 3082 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                     ELSE
                        WRITE( IOUTGR, 3085 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                     END IF
                  END IF   
               END IF
            ELSE
               IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                  IF ( STARTP ) THEN
C
C --------- CONTINUATION OF A PARAMETER ASSIGNMENT.
C
                     IF ( .NOT. SETF ) THEN
                        WRITE( IOUTGR, 3090 ) FIELD7
                     ELSE
                        WRITE( IOUTGR, 3091 ) FIELD7
                     END IF
                  ELSE
                     INFORM = 56
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
                     GO TO 800
                  END IF
               END IF
            END IF
         ELSE
            STARTP = .FALSE.
            IF ( FIELD1( 1: 1 ) .EQ. 'F' ) THEN
C
C  SET THE FUNCTION VALUE.
C
               IF ( FIELD1( 2: 2 ) .EQ. ' ' ) THEN
                  SETF = .TRUE.
                  ENDF = .FALSE.
C
C --------- START G.
C
                  WRITE( IOUTGR, 3100 ) FIELDI(  8 )( 1 : 6 ),
     *            FIELDI( 2 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ), FIELD7
               ELSE
                  IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                     IF ( SETF ) THEN
C
C --------- CONTINUATION OF G.
C
                        WRITE( IOUTGR, 3110 ) FIELD7
                     ELSE
                        INFORM = 56
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
                        GO TO 800
                     END IF
                  END IF
               END IF
            ELSE
               IF ( FIELD1( 1: 1 ) .EQ. 'G' ) THEN
C
C  NO FUNCTION VALUE HAS BEEN SPECIFIED.
C
                  IF ( .NOT. SETF ) THEN
                     INFORM = 61
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
                     GO TO 800
                  END IF
C
C  SET THE FIRST DERIVATIVE VALUE.
C
                  IF ( FIELD1( 2: 2 ) .EQ. ' ' ) THEN
                     IF ( .NOT. SETG ) THEN
                        SETG = .TRUE.
C
C --------- START GDASH.
C
                        IF ( .NOT. ENDF ) THEN
                           WRITE( IOUTGR, 3120 )
                           ENDF = .TRUE.
                        END IF
                     END IF
                     WRITE( IOUTGR, 3130 ) FIELDI( 2 )( 1 : 6 ),
     *                      FIELDI( 10 )( 1 : 6 ), 2, FIELD7
                  ELSE
                     IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                        IF ( SETG ) THEN
C
C --------- CONTINUATION OF GDASH.
C
                           WRITE( IOUTGR, 3140 ) FIELD7
                        ELSE
                           INFORM = 56
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
                           GO TO 800
                        END IF
                     END IF
                  END IF
               ELSE
                  IF ( FIELD1( 1: 1 ) .EQ. 'H' ) THEN
C
C  SET THE SECOND DERIVATIVE VALUE.
C
                     IF ( FIELD1( 2: 2 ) .EQ. ' ' ) THEN
                        IF ( .NOT. SETH ) THEN
C
C  THE FIRST DERIVATIVE HAS NOT BEEN SET.
C
                           IF ( .NOT. SETG ) THEN
                              INFORM = 62
                              IF ( IOUT .GT. 0 ) WRITE( IOUT, 2620 )
                              GO TO 800
                           END IF
                           SETH = .TRUE.
                        END IF
C
C --------- SET G2DASH.
C
                        WRITE( IOUTGR, 3130 ) FIELDI( 2 )( 1 : 6 ),
     *                         FIELDI( 10 )( 1 : 6 ), 3, FIELD7
                     ELSE
                        IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                           IF ( SETH ) THEN
C
C --------- CONTINUATION OF G2DASH.
C
                              WRITE( IOUTGR, 3140 ) FIELD7
                           ELSE
                              INFORM = 56
                              IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
                              GO TO 800
                           END IF
                        END IF
                     END IF
                  ELSE
                     INFORM = 56
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
                     GO TO 800
                  END IF
               END IF
            END IF
         END IF
      END IF
      GO TO 100
C
C  THE END OF THE INPUT FILE HAS BEEN REACHED BEFORE THE ENDATA CARD.
C
  590 CONTINUE 
C
C  IF THE ELEMENTS CARD HAS NOT BEEN ENCOUNTERED, EXIT.
C
      IF ( DEFNAM ) THEN
         INFORM = 52
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2520 )
         RETURN
      END IF
      IF ( NGRTYP .GT. 0 ) GO TO 930
      IF ( IOUT .GT. 0 .AND. IPRINT .NE. 0 ) WRITE( IOUT, 2010 )
C
C  A DUMMY ROUTINE WILL BE SUBSTITUTED.
C
  600 CONTINUE 
C
C  WRITE A DUMMY GROUPS ROUTINE.
C
      IF ( SINGLE ) THEN
         WRITE( IOUTGR, 3001 ) ( FIELDI( I )( 1 : 6 ), I = 1, 4 ),
     *              FIELDI( 11 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *              FIELDI(  6 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *              FIELDI(  7 )( 1 : 6 ), 
     *            ( FIELDI(  I )( 1 : 6 ), I = 15, 19 ),
     *              FIELDI(  8 )( 1 : 6 ), FIELDI( 20 )( 1 : 6 ),
     *              FIELDI(  3 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *            ( FIELDI(  I )( 1 : 6 ), I = 15, 20 ),
     *              FIELDI(  8 )( 1 : 6 ), 
     *              FIELDI(  6 )( 1 : 6 ), FIELDI( 15 )( 1 : 6 ),
     *              FIELDI( 12 )( 1 : 6 ), FIELDI( 16 )( 1 : 6 ),
     *              FIELDI(  7 )( 1 : 6 ), FIELDI( 17 )( 1 : 6 ),
     *              FIELDI(  2 )( 1 : 6 ), FIELDI(  3 )( 1 : 6 ),
     *              FIELDI(  4 )( 1 : 6 ), FIELDI( 18 )( 1 : 6 ),
     *              FIELDI( 11 )( 1 : 6 ), FIELDI( 19 )( 1 : 6 ),
     *              PNAME
      ELSE
         WRITE( IOUTGR, 3000 ) ( FIELDI( I )( 1 : 6 ), I = 1, 4 ),
     *              FIELDI( 11 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *              FIELDI(  6 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *              FIELDI(  7 )( 1 : 6 ), 
     *            ( FIELDI(  I )( 1 : 6 ), I = 15, 19 ),
     *              FIELDI(  8 )( 1 : 6 ), FIELDI( 20 )( 1 : 6 ),
     *              FIELDI(  3 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *            ( FIELDI(  I )( 1 : 6 ), I = 15, 20 ),
     *              FIELDI(  8 )( 1 : 6 ), 
     *              FIELDI(  6 )( 1 : 6 ), FIELDI( 15 )( 1 : 6 ),
     *              FIELDI( 12 )( 1 : 6 ), FIELDI( 16 )( 1 : 6 ),
     *              FIELDI(  7 )( 1 : 6 ), FIELDI( 17 )( 1 : 6 ),
     *              FIELDI(  2 )( 1 : 6 ), FIELDI(  3 )( 1 : 6 ),
     *              FIELDI(  4 )( 1 : 6 ), FIELDI( 18 )( 1 : 6 ),
     *              FIELDI( 11 )( 1 : 6 ), FIELDI( 19 )( 1 : 6 ),
     *              PNAME
      END IF
      WRITE( IOUTGR, 3009 ) FIELDI( 20 )( 1 : 6 )
      WRITE( IOUTGR, 3210 )
      INFORM = 0
      RETURN
C
C  INSUFFICIENT SPACE TO CONTINUE CONSTRUCTION.
C
  700 CONTINUE
      IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 ) INCRSE( - INFORM )
      RETURN
C
C  SUBROUTINE INCOMPLETE.
C
  800 CONTINUE
      IF ( IOUT .GT. 0 ) WRITE( IOUT, 2990 ) LINENO, NULINE
      RETURN
C
C  SUBROUTINE SUCCESSFULLY COMPLETED.
C
  900 CONTINUE
      IF ( .NOT. FIRSTG ) THEN
C
C  FINISH OF THE PREVIOUS GROUP, IF ANY.
C
         IF ( .NOT. SETH ) THEN
            INFORM = 63
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2630 )
            GO TO 800
         END IF
         IF ( .NOT. SETG ) THEN
            INFORM = 62
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2620 )
            GO TO 800
         END IF
C
C ---------- WIND UP F AND G.
C
         IF ( SETF ) THEN
            WRITE( IOUTGR, 3190 )
         ELSE
            INFORM = 61
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
            GO TO 800
         END IF
         IF ( ITYPE .LT. NGRTYP ) WRITE( IOUTGR, 3191 ) NLOOP
      END IF
C
C ---------- END DO LOOP.
C
      WRITE( IOUTGR, 3200 ) NLOOP
  910 CONTINUE
C
C ---------- SUCCESSFUL RUN. WIND UP OUTPUT.
C
      WRITE( IOUTGR, 3210 )
      INFORM = 0
C
C   CHECK THAT ALL ELEMENT TYPES HAVE BEEN DEFINED.
C
  930 CONTINUE
      DO 940 ITYPE = 1, NGRTYP
         IF ( .NOT. LDEFND( ITYPE ) ) THEN
            INFORM = 53
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2530 ) GTYPES( ITYPE )
         END IF
  940 CONTINUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 1000 FORMAT( A72 )
 1010 FORMAT( A160 )
 2000 FORMAT( ' ** Exit from MAKEGR - insufficient space.',
     *        ' Increase size of ', A6 )
 2010 FORMAT( ' ** Exit from MAKEGR - warning.',
     *        ' First card not groups. ', /, '    A dummy',
     *        ' routine will be substituted ' )
 2020 FORMAT( ' ** Exit from MAKEGR - indicator card not recognised ' )
 2190 FORMAT( ' ** Exit from MAKEGR - group type not recognised:',
     *        ' name is ', A10 )
 2510 FORMAT( ' ** Exit from MAKEGR -',
     *        ' Name on card not that specified on input ' )
 2520 FORMAT( ' ** Exit from MAKEGR -',
     *        ' data file incomplete. No ENDATA card ' )
 2530 FORMAT( ' ** Exit from MAKEGR - warning, group type ', A8,
     *        ' undefined ' )
 2540 FORMAT( ' ** Exit from MAKEGR -',
     *        ' unrecognised field 1 in TEMPORARIES section' )
 2550 FORMAT( ' ** Exit from MAKEGR -',
     *        ' unrecognised field 1 in GLOBALS section' )
 2560 FORMAT( ' ** Exit from MAKEGR -',
     *        ' unrecognised field 1 in INDIVIDUAL section' )
 2570 FORMAT( ' ** Exit from MAKEGR -',
     *        ' undefined parameter in GLOBALS section' )
 2580 FORMAT( ' ** Exit from MAKEGR -',
     *        ' undefined parameter in INDIVIDUALS section' )
 2590 FORMAT( ' ** Exit from MAKEGR - repeated parameter name ', A8 )
 2610 FORMAT( ' ** Exit from MAKEGR - function not set '  )
 2620 FORMAT( ' ** Exit from MAKEGR -',
     *        ' one or more first derivatives not set ' )
 2630 FORMAT( ' ** Exit from MAKEGR -',
     *        ' one or more second derivatives not set ' )
 2640 FORMAT( ' ** Exit from MAKEGR - group type already defined ' )
C ** Correction 2. 19/07/93: 2 lines added **
 2730 FORMAT( ' ** Exit from MAKEGR - field 2 or 3 not blank on',
     *        ' A, F, G or H card ' )
C ** Correction 2. 19/07/93: end of correction **
 2900 FORMAT( ' ' )
 2970 FORMAT( ' Line ', I5, 4X, A160 )
 2980 FORMAT( ' Line ', I5, '.', I2, 1X, A65 )
 2990 FORMAT( ' Line ', I5, 4X, A65 )
 3000 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 4( A6, ', ' ), A6, ' )', /,
     *        '      INTEGER ', 3( A6, ', ' ), A6, /,
     *        '      INTEGER ', 3( A6, ', ' ), A6, /,
     *        '      LOGICAL ', A6, /,
     *        '      INTEGER ', A6, '(', A6, '), ', A6, '(', A6, 
     *                       '), ', A6, '(', A6, ')', /,
     *        '      DOUBLE PRECISION ', A6, '(', A6, ',3), ',
     *                                   A6, '(', A6, '), ', A6, 
     *                                '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C' )
 3001 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 4( A6, ', ' ), A6, ' )', /,
     *        '      INTEGER ', 3( A6, ', ' ), A6, /,
     *        '      INTEGER ', 3( A6, ', ' ), A6, /,
     *        '      LOGICAL ', A6, /,
     *        '      INTEGER ', A6, '(', A6, '), ', A6, '(', A6, 
     *                       '), ', A6, '(', A6, ')', /,
     *        '      REAL             ', A6, '(', A6, ',3), ',
     *                                   A6, '(', A6, '), ', A6, 
     *                                '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C' )
 3002 FORMAT( '      INTEGER ', A6, ', ', A6, ', ', A6, ', ', A6 )
 3009 FORMAT( '      ', A6, ' = 0' )
 3010 FORMAT( ( '      INTEGER ', A6, 4( :, ', ', A6 ) ) )
 3019 FORMAT( ( '      REAL             ', A6, 4( :, ', ', A6 ) ) )
 3020 FORMAT( ( '      DOUBLE PRECISION ', A6, 4( :, ', ', A6 ) ) )
 3021 FORMAT( ( '      INTRINSIC ', A6, 4( :, ', ', A6 ) ) )
 3022 FORMAT( ( '      EXTERNAL ', A6, 4( :, ', ', A6 ) ) )
 3023 FORMAT( ( '      LOGICAL ', A6, 4( :, ', ', A6 ) ) )
 3030 FORMAT( '      ', A6, ' = ', A41 )
 3031 FORMAT( '      IF (', A6, ') ', A6, ' = ', A41 )
 3032 FORMAT( '      IF (.NOT.', A6, ') ', A6, ' = ', A41 )
 3040 FORMAT( '     *         ', A41 )
 3050 FORMAT( '      DO ', I5, 1X, A6, ' = 1, ', A6, /,
     *        '       ', A6, ' = ', A6, '(', A6, ')', /,
     *        '       ', A6, ' = ', A6, '(', A6, ')', /,
     *        '       IF ( ', A6, ' .EQ. 0 ) GO TO ', I5, /,
     *        '       ', A6, ' = ', A6, '(', A6, ') - 1' )
C ** Correction 6. 15/08/95: 2 lines corrected **
 3051 FORMAT( '       GO TO (', 8( I5, :, ',' ), /,
     *      ( '     *        ', 8( I5, :, ',' ) ) )
C ** Correction 6. 15/08/95: end of correction **
 3052 FORMAT( '     *        ', 48X, '), ', A6 )
 3060 FORMAT( 'C', /, 'C  GROUP TYPE : ', A8, /, 'C' )
 3061 FORMAT( I5, '  CONTINUE ' )
 3062 FORMAT( '       ', A6, '= ', A6, '(', A6, ')' )
 3063 FORMAT( '       ', A6, '= ', A6, '(', A6, '+', I6, ')' )
 3080 FORMAT( '       ', A6, '= ', A41 )
 3081 FORMAT( '       IF (', A6, ') ', A6, ' = ', A41 )
C ** Correction 3. 12/01/94: 1 line corrected **
 3082 FORMAT( '       IF (.NOT.', A6, ') ', A6, '=', A41 )
C ** Correction 3. 12/01/94: end of correction **
 3083 FORMAT( '        ', A6, '          = ', A41 )
 3084 FORMAT( '        IF (', A6, ') ', A6, ' = ', A41 )
C ** Correction 4. 12/01/94: 1 line corrected **
 3085 FORMAT( '        IF (.NOT.', A6, ')', A6, '=', A41 )
C ** Correction 4. 12/01/94: end of correction **
 3090 FORMAT( '     *          ', A41 )
 3091 FORMAT( '     *           ', A41 )
 3100 FORMAT( '       IF ( .NOT. ', A6, ' ) THEN', /,
     *        '        ', A6, '(', A6, ',1)= ', A41 )
 3110 FORMAT( '     *                  ', A41 )
 3120 FORMAT( '       ELSE' )
 3130 FORMAT( '        ', A6, '(', A6, ',', I1, ')= ', A41 )
 3140 FORMAT( '     *                         ', A41 )
 3190 FORMAT( '       END IF' )
 3191 FORMAT( '       GO TO', I6 )
 3200 FORMAT( I5,  ' CONTINUE' )
 3210 FORMAT( '      RETURN', /,  '      END' )
C
C  END OF MAKEGR.
C
      END
