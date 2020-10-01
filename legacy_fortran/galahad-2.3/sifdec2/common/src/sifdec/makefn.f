C  THIS VERSION: 26/02/2001 AT 09:30:00 AM.
C     ( Last modified on 15 Mar 2001 at 22:28:00 )
C ** Correction report.
C ** Correction -4. 13/02/02: Combine arguments in RANGE
C ** Correction -3. 14/01/02: Add arguments to RANGE/ELFUN (and change names)
C ** Correction -2. 29/12/01: Remove SETTYP and add extra argument to RANGES
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C ** Correction 1. 19/07/93: 12 lines added **
C ** Correction 2. 19/07/93: 2 lines added **
C ** Correction 3. 12/01/94: 1 line corrected **
C ** Correction 4. 12/01/94: 1 line corrected **
C ** Correction 5. 21/02/94: 1 line corrected **
C ** Correction 6. 21/02/94: 1 line corrected, 2 lines added **
C ** Correction 7. 15/08/95: 3 lines corrected **
C ** Correction 8. 15/08/95: 3 lines corrected **
C ** Correction 9. 26/02/01: 2 dummy arguments removed **
C ** End of Correction report.
      SUBROUTINE MAKEFN( INPUT , IOUT  , IOUTFN, IOUTRA, INFORM,
C ** Correction 9. 26/02/01: 2 dummy arguments removed **
     *                   NLMAX , NIMAX , NETMAX, NINMAX, NUMAX ,
     *                   NEL   , NELTYP, PNAME , ENAMES, INAMES, RENAME,
     *                   INNAME, LONAME, MINAME, EXNAME, ETYPES, LDEFND,
     *                   LENGTH, ITABLE, KEY   , IELV  , IINV  , INLIST,
     *                   EPNAME, IEPA  , NEPMAX, DEBUG , IJUMP ,
     *                   U     , SETVEC, NSETVC, SINGLE, NULINE, GOTLIN,
     *                   IPRINT )
      INTEGER            INPUT , IOUT  , IOUTFN, IOUTRA, INFORM
      INTEGER            NLMAX , NIMAX , NETMAX, NELTYP, NINMAX
      INTEGER            NEPMAX, LENGTH, NSETVC, NEL   , NUMAX , IPRINT
      LOGICAL            DEBUG , SINGLE, GOTLIN
      INTEGER            ITABLE( LENGTH ), IJUMP ( NLMAX  )
      INTEGER            IELV  ( NLMAX  ), IINV  ( NLMAX  )
      INTEGER            INLIST( LENGTH )
      INTEGER            IEPA  ( NLMAX  )
      LOGICAL            LDEFND( NLMAX  ), SETVEC( NSETVC )
      DOUBLE PRECISION   U     ( NUMAX  )
      CHARACTER * 12     KEY   ( LENGTH )
      CHARACTER * 8      PNAME
      CHARACTER * 10     EPNAME( NEPMAX ), EXNAME( NINMAX )
      CHARACTER * 10     INAMES( NIMAX  ), RENAME( NINMAX )
      CHARACTER * 10     LONAME( NINMAX ), INNAME( NINMAX )
      CHARACTER * 10     MINAME( NINMAX ), ENAMES( NETMAX )
      CHARACTER * 10     ETYPES( NLMAX  )
      CHARACTER * 160    NULINE
C
C  MAKE A FUNCTION EVALUATION SUBROUTINE AND A RANGE TRANSFORMATION
C  ----------------------------------------------------------------
C  SUBROUTINE FROM A GPS FUNCTION DATA FILE.
C  -----------------------------------------
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
C  ELEMENTS     PROBLEM NAME.
C  TEMPORARIES  NAMES OF ADDITIONAL PARAMETERS USED IN FUNCTION DEFS.
C  GLOBALS      GENERAL PARAMETER ASSIGNMENTS.
C  INDIVIDUALS  DEFINE THE TRANSFORMATION FROM THE ELEMENTAL TO THE
C               INTERNAL VARIABLES FOR ALL ELEMENTS WITH INTERNAL VARS.
C               SET FUNCTION AND DERIVATIVE VALUES AND MAKE
C               ELEMENT SPECIFIC PARAMETER ASSIGNMENTS.
C  ENDATA       END OF INPUT DATA.
C
C  DATA CARD DESCRIPTION.
C  ----------------------
C
C  SEE 'A PROPOSAL FOR A STANDARD DATA INPUT FORMAT FOR LARGE-SCALE
C       NONLINEAR PROGRAMMING PROBLEMS', SECTION 3,
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
C    INFORM = - 11 WHEN NUMAX NOT LARGE ENOUGH
C    INFORM = - 12 WHEN NSETVC NOT LARGE ENOUGH
C
      INTEGER          IFIELD, IFREE, IHVAR, IVAR, INTYPE, NMINAM
      INTEGER          JVAR, K, K1, K2, NH, NHESS, NVARS, ISETTY
      INTEGER          I, NINAME, NINNAM, NINVAR, NLOOP, NRENAM, NEXNAM
      INTEGER          ITYPE, J, IS, JS, NOVALS, NELV, NINV, NN, NENAME
      INTEGER          NPNAME, NINCRS, LINENO, IIRES, ILINES, NLINES
      INTEGER          MBLANK, MFIXED, MFREE, MNAME, MTEMP, MGLOB
      INTEGER          MINDIV, MENDAT, MAXNUL, MXRECL, NLONAM
      DOUBLE PRECISION ZERO, VALUES( 2 )
      LOGICAL          NOINTE, DEFNAM, ENDPAR, ENDGEN, FIRSTL, SETRAN
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      LOGICAL          STARTF, STARTG, STARTH, STARTP, QPROD
      LOGICAL          ENDOFF, ENDOFG, ENDOFH, FIXED, NOMORG
      PARAMETER        ( IIRES = 32 )
C ** Correction 5. 21/02/94: 1 line corrected **
      PARAMETER        ( NINCRS = 12 )
C ** Correction 5. 21/02/94: end of correction **
      CHARACTER * 2    FIELD1
      CHARACTER * 6    INCRSE( NINCRS )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      CHARACTER * 10   FIELD2
      CHARACTER * 8    FIELD3, FIELDS( 2 ), FIELDI( IIRES )
      CHARACTER * 12   FIELD
      CHARACTER * 41   FIELD7
      INTRINSIC        MIN
      EXTERNAL         HASHB , HASHC , GETVAL, OUTRAN
C
C  PARAMETER DEFINITIONS.
C
      PARAMETER        ( MXRECL = 160 )
      CHARACTER * 160    BLNKLN
      PARAMETER        ( MBLANK =  1, MFIXED =  2, MFREE  = 3  )
      PARAMETER        ( MNAME  =  4, MTEMP  =  5 )
      PARAMETER        ( MGLOB  =  6, MINDIV =  7, MENDAT =  8 )
      INTEGER          LENIND( MENDAT )
      CHARACTER * 12   INDIC8( MENDAT ), HEADER
      PARAMETER        ( MAXNUL = 20 )
      CHARACTER * 65   NULINA( MAXNUL )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      CHARACTER * 10   CQSQR, CQPROD
      PARAMETER      ( CQSQR  = '123456789S' )
      PARAMETER      ( CQPROD = '123456789P' )
C
C  DATA DECLARATIONS.
C
C ** Correction 6. 21/02/94: 1 line corrected, 2 lines added **
      DATA INCRSE / 'LENGTH', 'NINMAX', '      ', '      ', '      ',
     *              '      ', '      ', '      ', '      ', '      ',
     *              'NUMAX ', 'NSETVC'                               /
C ** Correction 6. 21/02/94: end of correction **
      DATA INDIC8( MBLANK ) / '            ' /, LENIND( MBLANK ) / 0  / 
      DATA INDIC8( MFIXED ) / 'FIXED FORMAT' /, LENIND( MFIXED ) / 12 /
      DATA INDIC8( MFREE  ) / 'FREE FORMAT ' /, LENIND( MFREE  ) / 11 /
      DATA INDIC8( MNAME  ) / 'ELEMENTS    ' /, LENIND( MNAME  ) / 8  /
      DATA INDIC8( MTEMP  ) / 'TEMPORARIES ' /, LENIND( MTEMP  ) / 11 / 
      DATA INDIC8( MGLOB  ) / 'GLOBALS     ' /, LENIND( MGLOB  ) / 7  /
      DATA INDIC8( MINDIV ) / 'INDIVIDUALS ' /, LENIND( MINDIV ) / 11 / 
      DATA INDIC8( MENDAT ) / 'ENDATA      ' /, LENIND( MENDAT ) / 6  /
      DATA FIELDI(  1 ) / 'ELFUN   ' /,  FIELDI(  2 ) / 'LFUVAL  ' /
      DATA FIELDI(  3 ) / 'FUVALS  ' /,  FIELDI(  4 ) / 'XVALUE  ' /
      DATA FIELDI(  5 ) / 'NCALCF  ' /,  FIELDI(  6 ) / 'ITYPEE  ' /
      DATA FIELDI(  7 ) / 'ISTAEV  ' /,  FIELDI(  8 ) / 'IELVAR  ' /
      DATA FIELDI(  9 ) / 'INTVAR  ' /,  FIELDI( 10 ) / 'ISTADH  ' /
      DATA FIELDI( 11 ) / 'ICALCF  ' /,  FIELDI( 12 ) / 'IFFLAG  ' /
      DATA FIELDI( 13 ) / 'IELEMN  ' /,  FIELDI( 14 ) / 'IELTYP  ' /
      DATA FIELDI( 15 ) / 'IHSTRT  ' /,  FIELDI( 16 ) / 'ILSTRT  ' /
      DATA FIELDI( 17 ) / 'IGSTRT  ' /,  FIELDI( 18 ) / 'EPVALU  ' /
      DATA FIELDI( 19 ) / 'ISTEPA  ' /,  FIELDI( 20 ) / 'IPSTRT  ' /
      DATA FIELDI( 21 ) / 'JCALCF  ' /,  FIELDI( 22 ) / 'LTYPEE  ' /  
      DATA FIELDI( 23 ) / 'LSTAEV  ' /,  FIELDI( 24 ) / 'LELVAR  ' /
      DATA FIELDI( 25 ) / 'LNTVAR  ' /,  FIELDI( 26 ) / 'LSTADH  ' /
      DATA FIELDI( 27 ) / 'LSTEPA  ' /,  FIELDI( 28 ) / 'LCALCF  ' /    
      DATA FIELDI( 29 ) / 'LFVALU  ' /,  FIELDI( 30 ) / 'LXVALU  ' /
      DATA FIELDI( 31 ) / 'LEPVLU  ' /,  FIELDI( 32 ) / 'IFSTAT  ' /
      DATA ZERO         / 0.0D+0 /
      IF ( IOUT .GT. 0 ) WRITE( IOUT, 2900 )
C
C  SET INITIAL VALUES FOR INTEGER VARIABLES.
C
      NINNAM = 0
      NRENAM = 0
      NLONAM = 0
      NMINAM = 0
      NEXNAM = 0
      LINENO = 0
      NLOOP  = NELTYP + 1
      INTYPE = 1
      ILINES = 0
      NLINES = 0
C
C  SET INITIAL VALUES FOR LOGICAL VARIABLES.
C
      DEFNAM = .FALSE.
      ENDPAR = .FALSE.
      STARTH = .FALSE.
      STARTP = .FALSE.
      ENDGEN = .FALSE.
      FIRSTL = .TRUE.
      NOINTE = .FALSE.
      FIXED  = .TRUE.
      GOTLIN = .FALSE.
C
C  CREATE A DICTIONARY OF THE INTERNAL VARIABLE NAMES USED.
C
      NINAME  = IINV( NELTYP + 1 ) - 1
      DO 20 I = 1, NINAME
         FIELD = INAMES( I ) // 'PF'
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
            RENAME( NRENAM ) = INAMES( I )
         END IF
   20 CONTINUE
C
C  INCLUDE THE NAMES OF THE ELEMENTAL VARIABLES USED IN THIS DICTIONARY.
C
      NENAME  = IELV( NELTYP + 1 ) - 1
      DO 30 I = 1, NENAME
         FIELD = ENAMES( I ) // 'PF'
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
            RENAME( NRENAM ) = ENAMES( I )
         END IF
   30 CONTINUE
C
C  INCLUDE THE NAMES OF THE ELEMENTAL PARAMETERS USED
C  IN THIS DICTIONARY.
C
      NPNAME  = IEPA( NELTYP + 1 ) - 1
      DO 40 I = 1, NPNAME
         FIELD = EPNAME( I ) // 'PF'
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
            RENAME( NRENAM ) = EPNAME( I )
         END IF
   40 CONTINUE
C
C  FIND WHICH ELEMENT TYPES HAVE AN INTERNAL REPRESENTATION.
C
      ISETTY             = 0
      DO 50 ITYPE        = 1, NELTYP
         LDEFND( ITYPE ) = .FALSE.
         IF ( IELV( ITYPE + 1 ) - IELV( ITYPE ) .EQ.
     *        IINV( ITYPE + 1 ) - IINV( ITYPE ) ) THEN
            IJUMP( ITYPE ) = 99998
            NOINTE         = .TRUE.
         ELSE
            IJUMP( ITYPE ) = ITYPE
         END IF
   50 CONTINUE
C
C  SET A BLANK LINE.
C
      DO 60 I = 1, MXRECL
         BLNKLN( I: I ) = ' '
   60 CONTINUE    
C
C  READ NEXT LINE.
C
  100 CONTINUE
      IF ( ILINES + 1 .GT. NLINES ) THEN
C
C  READ NEXT LINE FROM THE INPUT FILE.
C
         LINENO = LINENO + 1
         NULINE = BLNKLN
         IF ( FIXED ) THEN
            READ ( INPUT, 1000, END = 590, ERR = 590 ) NULINE
            IF ( IOUT .GT. 0 .AND. DEBUG ) WRITE( IOUT, 2990 )
     *           LINENO, NULINE
         ELSE
            READ ( INPUT, 1010, END = 590, ERR = 590 ) NULINE
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
C  CHECK THAT THE FIRST ENCOUNTERED INDICATOR CARD IS THE ELEMENTS CARD.
C
         IF ( .NOT. DEFNAM  ) THEN
            IF ( HEADER .NE. INDIC8( MNAME ) ) THEN
               IF ( NELTYP .GT. 0 ) GO TO 930
               IF ( IOUT .GT. 0 .AND. IPRINT .NE. 0 ) WRITE( IOUT, 2010)
               GOTLIN = .TRUE.
               GO TO 600
            ELSE
C
C  INDICATOR CARD IS ELEMENTS.
C  ---------------------------
C
               IF ( PNAME  .NE. NULINE( 15: 22 ) ) THEN
                  INFORM = 51
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2510 )
                  GO TO 800
               ELSE
                  DEFNAM = .TRUE.
C
C  -------- SET UP SUBROUTINE CALL FOR RANGE ROUTINE.
C
                  IF ( SINGLE ) THEN
                     WRITE( IOUTRA, 4001 ) PNAME
                  ELSE
                     WRITE( IOUTRA, 4000 ) PNAME
                  END IF
                  IF ( NELTYP .GT. 1 ) THEN
                     WRITE( IOUTRA, 4040 ) ( IJUMP( I ), I = 1, NELTYP )
                     WRITE( IOUTRA, 4050 )
                  END IF
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
C
C  INSERT THE LIST OF RESERVED INTEGER/REAL/LOGICAL VARIABLES INTO
C  THE DICTIONARY.
C
            DO 130 I = 1, IIRES
               FIELD = FIELDI( I ) // '  PF'
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
               WRITE( IOUTFN, 3001 ) FIELDI( 1 )(1:6),
     *   FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6), 
     * ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *   FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     *   FIELDI( 12 )(1:6), FIELDI( 32 )(1:6),
     *   FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     * ( FIELDI(  I )(1:6), I = 22, 32 ),
     *   FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     *   FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *   FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *   FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *   FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *   FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *   FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *   FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *   FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *   FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     *   PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     *   FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     *   FIELDI( 17 )(1:6), FIELDI( 20 )(1:6), FIELDI( 21 )(1:6)
            ELSE   
               WRITE( IOUTFN, 3000 ) FIELDI( 1 )(1:6),
     *   FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     * ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *   FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     *   FIELDI( 12 )(1:6), FIELDI( 32 )(1:6),
     *   FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     * ( FIELDI(  I )(1:6), I = 22, 32 ),
     *   FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     *   FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *   FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *   FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *   FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *   FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *   FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *   FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *   FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *   FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     *   PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     *   FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     *   FIELDI( 17 )(1:6), FIELDI( 20 )(1:6), FIELDI( 21 )(1:6)
            END IF
C
C --------- INSERT INTEGER DECLARATIONS.
C
            IF ( NINNAM .GT. 0 )
     *         WRITE( IOUTFN, 3010 ) ( INNAME( I ), I = 1, NINNAM )
C
C --------- INSERT REAL DECLARATIONS.
C
            IF ( NRENAM .GT. 0 ) THEN
               IF ( SINGLE ) THEN
                  WRITE( IOUTFN, 3019 ) ( RENAME( I ), I = 1, NRENAM )
               ELSE
                  WRITE( IOUTFN, 3020 ) ( RENAME( I ), I = 1, NRENAM )
               END IF
            END IF
C
C --------- INSERT LOGICAL DECLARATIONS.
C
            IF ( NLONAM .GT. 0 )
     *         WRITE( IOUTFN, 3023 ) ( LONAME( I ), I = 1, NLONAM )
C
C --------- INSERT INTRINSIC DECLARATIONS.
C
            IF ( NMINAM .GT. 0 )
     *         WRITE( IOUTFN, 3021 ) ( MINAME( I ), I = 1, NMINAM )
C
C --------- INSERT EXTERNAL DECLARATIONS.
C
            IF ( NEXNAM .GT. 0 )
     *         WRITE( IOUTFN, 3022 ) ( EXNAME( I ), I = 1, NEXNAM )
            WRITE( IOUTFN, 3009 ) FIELDI( 32 )(1:6)
         END IF
C
C  THE GENERAL PARAMETER ASSIGNMENTS HAVE BEEN COMPLETED.
C  CONTINUE WITH THE CONSTRUCTION OF THE GENERATED SUBROUTINE.
C
         IF ( INTYPE .GE. MINDIV .AND. .NOT. ENDGEN ) THEN
            ENDGEN = .TRUE.
C
C --------- START LOOP OVER ELEMENTS.
C
            WRITE( IOUTFN, 3050 ) NLOOP,
     *             FIELDI( 21 )(1:6), FIELDI(  5 )(1:6),
     *             FIELDI( 13 )(1:6), FIELDI( 11 )(1:6),
     *             FIELDI( 21 )(1:6),
     *             FIELDI( 16 )(1:6), FIELDI(  7 )(1:6),
     *             FIELDI( 13 )(1:6), FIELDI( 17 )(1:6),
     *             FIELDI(  9 )(1:6), FIELDI( 13 )(1:6),
     *             FIELDI( 20 )(1:6), FIELDI( 19 )(1:6),
     *             FIELDI( 13 )(1:6),
     *             FIELDI( 12 )(1:6), FIELDI( 15 )(1:6),
     *             FIELDI( 10 )(1:6), FIELDI( 13 )(1:6)
            IF ( NELTYP .GT. 1 ) THEN
               WRITE( IOUTFN, 3051 )
     *            FIELDI( 14 )(1:6), FIELDI(  6 )(1:6),
     *            FIELDI( 13 )(1:6), ( I, I = 1, NELTYP )
               WRITE( IOUTFN, 3052 ) FIELDI( 14 )(1:6)
            END IF
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C
C  MAKE SURE THAT QUADRATIC HESSIAN TERMS ARE INCLUDED.
C
            DO 140 ITYPE = 1, MIN( 2, NELTYP )
C
C  DIAGONAL TERM.
C
               IF ( ETYPES( ITYPE ) .EQ. CQSQR ) THEN
                  WRITE( IOUTFN, 3060 ) ETYPES( ITYPE )
                  IF ( NELTYP .GT. 1 ) WRITE( IOUTFN, 3061 ) ITYPE
                  IF ( SINGLE ) THEN
                    WRITE( IOUTFN, 3053 ) 'E', 'E'
                  ELSE
                    WRITE( IOUTFN, 3053 ) 'D', 'D'
                  END IF
                  LDEFND( ITYPE ) = .TRUE.
                  ISETTY = ISETTY + 1
                  IF ( ISETTY .LT. NELTYP ) WRITE( IOUTFN, 3191 ) NLOOP
C
C  OFF-DIAGONAL TERM.
C
               ELSE IF ( ETYPES( ITYPE ) .EQ. CQPROD ) THEN
                  WRITE( IOUTFN, 3060 ) ETYPES( ITYPE )
                  IF ( NELTYP .GT. 1 ) WRITE( IOUTFN, 3061 ) ITYPE
                  IF ( SINGLE ) THEN
                    WRITE( IOUTFN, 3054 ) 'E', 'E', 'E'
                  ELSE
                    WRITE( IOUTFN, 3054 ) 'D', 'D', 'D'
                  END IF
                  LDEFND( ITYPE ) = .TRUE.
                  ISETTY = ISETTY + 1
                  IF ( ISETTY .LT. NELTYP ) WRITE( IOUTFN, 3191 ) NLOOP
               END IF
  140       CONTINUE   
         END IF
C
C  INDICATOR CARD IS ENDATA.
C  -------------------------
C
         IF ( INTYPE .EQ. MENDAT ) GO TO 900
         GO TO 100
      ELSE
C
C  CHECK THAT THE FIRST NON COMMMENT CARD IS THE ELEMENTS INDICATOR CARD.
C
         IF ( .NOT. DEFNAM  ) THEN
            IF ( NELTYP .GT. 0 ) GO TO 930
            IF ( IOUT .GT. 0 .AND. IPRINT .NE. 0 ) WRITE( IOUT, 2010 )
            GOTLIN = .TRUE.
            GO TO 600
         END IF
C
C  A DATA CARD HAS BEEN FOUND.
C  READ THE CHARACTER FIELDS 1 AND 2 FROM THE CARD.
C
         FIELD1 = NULINE(  2:  3 )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
         FIELD2 = NULINE(  5: 14 )
         IF ( INTYPE .EQ. MINDIV .AND. FIELD1 .EQ. 'R ' ) THEN
C
C  READ THE CHARACTER FIELDS 3 AND 5 FROM THE CARD.
C
            FIELDS( 1 ) = NULINE( 15: 22 )
            FIELDS( 2 ) = NULINE( 40: 47 )
C
C  CHECK TO SEE IF THERE IS ARE ANY NUMERICAL VALUES TO BE READ.
C
            NOVALS = 0
            IF ( FIELDS( 1 ) .NE. '        ' .AND.
     *           NULINE( 15: 15 ) .NE. '[' ) THEN
               NOVALS = 1
               CALL GETVAL( NULINE( 25: 36 ), VALUES( 1 ) )
               IF ( FIELDS( 2 ) .NE. '        ' .AND.
     *           NULINE( 40: 40 ) .NE. '[' ) THEN
                  NOVALS = 2
                  CALL GETVAL( NULINE( 50: 61 ), VALUES( 2 ) )
C
C  REMOVE FIELDS WITH NUMERICAL VALUES OF ZERO.
C
                  IF ( VALUES( 2 ) .EQ. ZERO ) THEN
                     NOVALS = 1
                  END IF
               END IF
               IF ( VALUES( 1 ) .EQ. ZERO ) THEN
                  IF ( NOVALS .EQ. 2 ) THEN
                     VALUES( 1 ) = VALUES( 2 )
                     FIELDS( 1 ) = FIELDS( 2 )
                  END IF
                  NOVALS = NOVALS - 1
               END IF
            END IF
         ELSE
C
C  READ THE CHARACTER FIELDS 3 AND 7 FROM THE CARD.
C
            FIELD3 = NULINE( 15: 22 )
            FIELD7 = NULINE( 25: 65 )
C ** Correction 1. 19/07/93: 12 lines added **
C
C  CHECK THAT FIELD3 IS BLANK ON 'A', 'F' AND 'G' CARDS.
C
            IF ( FIELD1( 1: 1 ) .EQ. 'A' .OR.
     *           FIELD1( 1: 1 ) .EQ. 'F' .OR. 
     *           FIELD1( 1: 1 ) .EQ. 'G' ) THEN
               IF ( FIELD3 .NE. '       ' ) THEN
                  INFORM = 72
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2720 )
                  GO TO 800
               END IF
            END IF
C ** Correction 1. 19/07/93: end of correction **
         END IF
      END IF
C
C  BRANCH ON THE VALUE OF INTYPE.
C
      GO TO ( 100, 100, 100, 100, 200, 300, 400, 900 ), INTYPE
C
C  INDICATOR CARD IS TEMPORARIES.
C  ------------------------------
C
  200 CONTINUE
C
C  CHECK TO SEE IF THE PARAMETER IS INTEGER, REAL, LOGICAL OR A FUNCTION.
C
      IF ( FIELD1 .NE. 'I ' .AND. FIELD1 .NE. 'R ' .AND.
     *     FIELD1 .NE. 'M ' .AND. FIELD1 .NE. 'F ' .AND.
     *     FIELD1 .NE. 'L ' ) THEN
         INFORM = 54
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2540 )
         GO TO 800
      END IF
C
C  IF THE PARAMETER IS A FUNCTION, CHECK TO SEE THAT THE NAME HAS
C  NOT ALREADY BEEN USED.
C
      IF ( FIELD1 .EQ. 'F ' ) THEN
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
         FIELD = FIELD2 // 'FU'
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
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
         FIELD = FIELD2 // 'PF'
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
C  INDICATOR CARD IS GLOBALS.
C  --------------------------
C
  300 CONTINUE
      IF ( FIELD1 .EQ. 'A ' .OR. FIELD1 .EQ. 'I ' .OR.
     *     FIELD1 .EQ. 'E ' ) THEN
         STARTP = .TRUE.
C
C  START A PARAMETER ASSIGNMENT. CHECK TO SEE THAT THE PARAMETER HAS
C  BEEN DEFINED.
C
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
         FIELD = FIELD2 // 'PF'
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
            WRITE( IOUTFN, 3030 ) FIELD2(1:6), FIELD7
C
C --------- MAKE CONDITIONAL PARAMETER ASSIGNMENTS.
C
         ELSE   
C
C  CHECK THAT THE LOGICAL VARIABLE HAS BEEN DEFINED.
C
            FIELD = FIELD3 // '  PF'
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
               WRITE( IOUTFN, 3031 ) FIELD2(1:6),
     *                               FIELD3(1:6), FIELD7
            ELSE
               WRITE( IOUTFN, 3032 ) FIELD2(1:6),
     *                               FIELD3(1:6), FIELD7
            END IF   
         END IF
      ELSE
         IF ( FIELD1( 2: 2 ) .EQ. '+' .AND. STARTP ) THEN
C
C --------- CONTINUE A PARAMETER ASSIGNMENT.
C
            WRITE( IOUTFN, 3040 ) FIELD7
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
C  CHECK IF A NEW ELEMENT HAS BEEN ENCOUNTERED.
C
      IF ( FIELD1 .EQ. 'T ' ) THEN
C
C  CHECK TO SEE IF THE RANGE OF A NEW ELEMENT IS TO BE DEFINED.
C
         IF ( FIRSTL ) THEN
C
C  CHECK IF THIS IS THE FIRST ELEMENT.
C
            FIRSTL = .FALSE.
         ELSE
C
C  FINISH OF THE PREVIOUS ELEMENT, IF ANY.
C
            IF ( STARTH ) THEN
               IF ( .NOT. ENDOFH ) THEN
                  DO 410 IHVAR = 1, NHESS
C
C --------- SET A COMPONENT OF H.
C
                     IF ( .NOT. SETVEC( IHVAR ) ) THEN
                        IF ( SINGLE ) THEN
                           WRITE( IOUTFN, 3162 ) FIELDI(  3 )(1:6),
     *                            FIELDI( 15 )(1:6), IHVAR
                        ELSE
                           WRITE( IOUTFN, 3161 ) FIELDI(  3 )(1:6),
     *                           FIELDI( 15 )(1:6), IHVAR
                        END IF   
                     END IF   
  410             CONTINUE
                  ENDOFH = .TRUE.
               END IF
C
C ---------- WIND UP H.
C
               WRITE( IOUTFN, 3180 )
            END IF
            IF ( STARTG ) THEN
C
C  SET THE REMAINING GRADIENT COMPONENTS TO ZERO.
C
               IF ( .NOT. ENDOFG ) THEN
                  DO 415 IVAR = 1, NINVAR
C
C --------- SET A COMPONENT OF G.
C
                     IF ( .NOT. SETVEC( IVAR ) ) THEN
                        IF ( SINGLE ) THEN
                           WRITE( IOUTFN, 3132 )
     *                        FIELDI(  3 )(1:6),
     *                        FIELDI( 17 )(1:6), IVAR
                        ELSE
                           WRITE( IOUTFN, 3131 )
     *                        FIELDI(  3 )(1:6),
     *                        FIELDI( 17 )(1:6), IVAR
                        END IF   
                     END IF   
  415             CONTINUE
                  ENDOFG = .TRUE.
               END IF   
C
C ---------- WIND UP F AND G
C
            END IF
            IF ( STARTF ) THEN
               WRITE( IOUTFN, 3190 )
            ELSE
               INFORM = 61
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
               GO TO 800
            END IF
            IF ( ISETTY .LT. NELTYP ) WRITE( IOUTFN, 3191 ) NLOOP
         END IF
C
C  FIND ITYPE, THE ELEMENT TYPE.
C
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
         FIELD = FIELD2 // 'ET'
         CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
C
C  THE ELEMENT TYPE IS UNKNOWN.
C
         IF ( IFIELD .LE. 0 ) THEN
            INFORM = 9
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2090 )
            GO TO 800
         END IF
C
C --------- FIND TYPE OF CURRENT ELEMENT.
C
         ITYPE = INLIST( IFIELD )
         WRITE( IOUTFN, 3060 ) FIELD2
         IF ( NELTYP .GT. 1 ) WRITE( IOUTFN, 3061 ) ITYPE
         IF ( LDEFND( ITYPE ) ) THEN
            INFORM = 67
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2670 )
            GO TO 800
         ELSE
            LDEFND( ITYPE ) = .TRUE.
            ISETTY = ISETTY + 1
         END IF
C
C  FIND THE ROW AND COLUMN DIMENSIONS (NINV AND NELV, RESP.) OF THE
C  TRANSFORMATION MATRIX U. U IS STORED IN VECTOR FORM BY COLUMNS.
C
         IS   = IINV( ITYPE ) - 1
         JS   = IELV( ITYPE ) - 1
         NELV = IELV( ITYPE + 1 ) - IELV( ITYPE )
         NINV = IINV( ITYPE + 1 ) - IINV( ITYPE )
         NN   = NINV * NELV
         IF ( NN .GT. NUMAX ) THEN
            INFORM = - 11
            GO TO 700
         END IF
C
C --------- FIND TYPE OF CURRENT ELEMENT.
C
         IF ( NELV .GT. NINV ) WRITE( IOUTRA, 4060 ) FIELD2, ITYPE
C
C  INITIALIZE U AS THE ZERO MATRIX.
C
         DO 420 I  = 1, NN
            U( I ) = ZERO
  420    CONTINUE
         SETRAN = NELV .GT. NINV
C
C --------- SET ELEMENTAL VARIABLES.
C
         K1       = IELV( ITYPE )
         K2       = IELV( ITYPE + 1 ) - 1
         DO 430 K = K1, K2
            IVAR  = K - K1 + 1
            WRITE( IOUTFN, 3070 ) ENAMES( K ), FIELDI(  4 )(1:6),
     *          FIELDI(  8 )(1:6), FIELDI( 16 )(1:6), IVAR
  430    CONTINUE
C
C --------- SET ELEMENTAL PARAMETERS.
C
         K1       = IEPA( ITYPE )
         K2       = IEPA( ITYPE + 1 ) - 1
         DO 435 K = K1, K2
            IVAR  = K - K1 + 1
            WRITE( IOUTFN, 3071 ) EPNAME( K ), FIELDI( 18 )(1:6),
     *          FIELDI( 20 )(1:6), IVAR
  435    CONTINUE
C
C  FIND THE NUMBER OF INTERNAL VARIABLES AND THE REQUIRED SIZE OF
C  THE LOWER TRIANGULAR PORTION OF THE HESSIAN MATRIX.
C
         K1     = IINV( ITYPE )
         K2     = IINV( ITYPE + 1 ) - 1
         NINVAR = K2 - K1 + 1
         NHESS  = NINVAR * ( NINVAR + 1 ) / 2
         NVARS  = 0
         NH     = 0
         STARTP = .FALSE.
         STARTF = .FALSE.
         STARTG = .FALSE.
         STARTH = .FALSE.
         ENDOFF = .TRUE.
         ENDOFG = .TRUE.
         ENDOFH = .TRUE.
         NOMORG = .FALSE.
      ELSE
         IF ( FIELD1 .EQ. 'R ' ) THEN
C
C  THE RANGE TRANSFORMATION MATRIX U IS NOW DEFINED, ENTRY BY ENTRY.
C  DETERMINE WHICH INTERNAL VARIABLE IS GIVEN IN FIELD2.
C
            DO 440 I = 1, NINV
               IF ( FIELD2 .EQ. INAMES( IS + I ) ) GO TO 450
  440       CONTINUE
C
C  THE INTERNAL VARIABLE NAME IS UNRECOGNISED.
C
            INFORM = 65
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2650 )
            GO TO 800
C
C  THE INTERNAL VARIABLE IS THE I-TH IN THE LIST.
C
  450       CONTINUE
C
C  DETERMINE WHICH ELEMENTAL VARIABLE(S) OCCUR IN FIELDS.
C
            IF ( NOVALS .GT. 0 ) THEN
               DO 480 K    = 1, NOVALS
                  DO 460 J = 1, NELV
                     IF ( FIELDS( K ) .EQ. ENAMES( JS + J ) ) GO TO 470
  460             CONTINUE
C
C  THE ELEMENTAL VARIABLE NAME IS UNRECOGNISED.
C
                  INFORM = 66
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2660 )
                  GO TO 800
C
C  THE ELEMENTAL VARIABLE IS THE J-TH IN THE LIST.
C
  470             CONTINUE
C
C  INSERT THE VALUE OF THE NEW NONZERO INTO U.
C
                  U( NINV * ( J - 1 ) + I ) = VALUES( K )
  480          CONTINUE
            END IF
         ELSE
            IF ( FIELD1( 1: 1 ) .EQ. 'A' .OR. FIELD1( 1: 1 ) 
     *           .EQ. 'I' .OR. FIELD1( 1: 1 ) .EQ. 'E' ) THEN
               IF ( STARTF ) THEN
                  IF ( .NOT. ENDOFF ) THEN
                     WRITE( IOUTFN, 3120 )
                     ENDOFF = .TRUE.
                  END IF
               END IF
               IF ( STARTG ) THEN
                  IF ( ENDOFG .AND. .NOT. NOMORG ) THEN
                     WRITE( IOUTFN, 3150 ) FIELDI( 12 )(1:6)
                     NOMORG = .TRUE.
                  END IF
               END IF
C
C  START A PARAMETER ASSIGNMENT.
C
               IF ( FIELD1( 2: 2 ) .EQ. ' ' ) THEN
                  STARTP = .TRUE.
C
C  SET UP THE TRANSFORMATIONS FOR THE ELEMENT.
C
                  IF ( SETRAN ) THEN
                      CALL OUTRAN( NELV, NINV, U, IOUTFN, IOUTRA,
     *                             ENAMES( JS + 1 ), INAMES( IS + 1 ),
     *                             SINGLE )
                      SETRAN = .FALSE.
                  END IF
C
C  CHECK TO SEE THAT THE PARAMETER HAS BEEN DEFINED.
C
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
                  FIELD = FIELD2 // 'PF'
                  CALL HASHC (LENGTH, 12, FIELD, KEY, ITABLE, IFIELD)
                  IF ( IFIELD .LE. 0 ) THEN
                     INFORM = 58
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2580 )
                     GO TO 800
                  END IF
C
C --------- MAKE ELEMENT-SPECIFIC PARAMETER ASSIGNMENTS.
C
                  IF ( FIELD1( 1: 1 ) .EQ. 'A' ) THEN
                     IF ( .NOT. STARTF ) THEN
                        WRITE( IOUTFN, 3080 ) FIELD2(1:6), FIELD7
                     ELSE
                        WRITE( IOUTFN, 3083 ) FIELD2(1:6), FIELD7
                     END IF
C
C --------- MAKE CONDITIONAL PARAMETER ASSIGNMENTS.
C
                  ELSE   
C
C  CHECK THAT THE LOGICAL VARIABLE HAS BEEN DEFINED.
C
                     FIELD = FIELD3 // '  PF'
                     CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE,
     *                            IFIELD )
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
                        IF ( .NOT. STARTF ) THEN
                           WRITE( IOUTFN, 3081 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                        ELSE
                           WRITE( IOUTFN, 3084 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                        END IF
                     ELSE
                     IF ( .NOT. STARTF ) THEN
                           WRITE( IOUTFN, 3082 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                        ELSE
                           WRITE( IOUTFN, 3085 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                        END IF
                     END IF   
                  END IF
               ELSE
                  IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                     IF ( STARTP ) THEN
C
C --------- CONTINUATION OF A PARAMETER ASSIGNMENT.
C
                     IF ( .NOT. STARTF ) THEN
                        WRITE( IOUTFN, 3090 ) FIELD7
                     ELSE
                        WRITE( IOUTFN, 3091 ) FIELD7
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
                     STARTF = .TRUE.
                     ENDOFF = .FALSE.
C
C  SET UP THE TRANSFORMATIONS FOR THE ELEMENT.
C
                     IF ( SETRAN ) THEN
                         CALL OUTRAN( NELV, NINV, U, IOUTFN, IOUTRA,
     *                             ENAMES( JS + 1 ), INAMES( IS + 1 ),
     *                             SINGLE )
                         SETRAN = .FALSE.
                     END IF
C
C --------- START F.
C
                     WRITE( IOUTFN, 3100 ) FIELDI( 12 )(1:6),
     *               FIELDI( 3 )(1:6), FIELDI( 13 )(1:6), FIELD7
                  ELSE
                     IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                        IF ( STARTF ) THEN
C
C --------- CONTINUATION OF F.
C
                           WRITE( IOUTFN, 3110 ) FIELD7
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
                     IF ( .NOT. STARTF ) THEN
                        INFORM = 61
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
                        GO TO 800
                     END IF
C
C  SET THE GRADIENT VALUES.
C
                     IF ( FIELD1( 2: 2 ) .EQ. ' ' ) THEN
                        IF ( .NOT. STARTG ) THEN
                           STARTG = .TRUE.
                           ENDOFG = .FALSE.
C
C  USE THE LOGICAL ARRAY SETVEC TO ENSURE THAT ALL GRADIENTS ARE SET.
C
                           IF ( NINVAR .GT. NSETVC ) THEN
                              INFORM = - 12
                              GO TO 700
                           END IF
                           DO 510 I       = 1, NINVAR
                              SETVEC( I ) = .FALSE.
  510                      CONTINUE
C
C --------- START G.
C
                           IF ( .NOT. ENDOFF ) THEN
                              WRITE( IOUTFN, 3120 )
                              ENDOFF = .TRUE.
                           END IF
                        END IF
C
C  FIND WHICH COMPONENT IS TO BE SET.
C
                        DO 520 K = K1, K2
                           IVAR  = K - K1 + 1
                           IF ( FIELD2 .EQ. INAMES( K ) ) GO TO 525
  520                   CONTINUE
C
C  THE COMPONENT NAME IS UNRECOGNISED.
C
                        INFORM = 60
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2600 )
                        GO TO 800
  525                   CONTINUE
C
C --------- SET A COMPONENT OF G.
C
                        IF ( SETVEC( IVAR ) ) THEN
                           INFORM = 69
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2690 )
                           GO TO 800
                        END IF
                        SETVEC( IVAR ) = .TRUE.
                        NVARS  = NVARS + 1
                        ENDOFG = NVARS .EQ. NINVAR
                        WRITE( IOUTFN, 3130 ) FIELDI(  3 )(1:6),
     *                         FIELDI( 17 )(1:6), IVAR, FIELD7
                     ELSE
                        IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                           IF ( STARTG .AND. .NOT. NOMORG ) THEN
C
C --------- CONTINUATION OF G.
C
                              WRITE( IOUTFN, 3140 ) FIELD7
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
C  SET THE HESSIAN VALUES.
C
                        IF ( FIELD1( 2: 2 ) .EQ. ' ' ) THEN
                           IF ( .NOT. STARTH ) THEN
C
C  SET THE REMAINING GRADIENT COMPONENTS TO ZERO.
C
                              IF ( .NOT. STARTG ) THEN
                                 DO 530 IVAR = 1, NINVAR
C
C --------- SET A COMPONENT OF G.
C
                                    IF ( SINGLE ) THEN
                                       WRITE( IOUTFN, 3132 )
     *                                    FIELDI(  3 )(1:6),
     *                                    FIELDI( 17 )(1:6), IVAR
                                    ELSE
                                       WRITE( IOUTFN, 3131 )
     *                                    FIELDI(  3 )(1:6),
     *                                    FIELDI( 17 )(1:6), IVAR
                                    END IF   
  530                            CONTINUE
                                 STARTG = .TRUE.
                                 ENDOFG = .TRUE.
                              END IF
                              IF ( .NOT. ENDOFG ) THEN
                                 DO 535 IVAR = 1, NINVAR
C
C --------- SET A COMPONENT OF G.
C
                                    IF ( .NOT. SETVEC( IVAR ) ) THEN
                                       IF ( SINGLE ) THEN
                                          WRITE( IOUTFN, 3132 )
     *                                       FIELDI(  3 )(1:6),
     *                                       FIELDI( 17 )(1:6), IVAR
                                       ELSE
                                          WRITE( IOUTFN, 3131 )
     *                                       FIELDI(  3 )(1:6),
     *                                       FIELDI( 17 )(1:6), IVAR
                                       END IF   
                                    END IF   
  535                            CONTINUE
                                 ENDOFG = .TRUE.
                              END IF
                              IF ( .NOT. NOMORG ) THEN
                                 WRITE( IOUTFN, 3150 ) 
     *                                  FIELDI( 12 )(1:6)
                                 NOMORG = .TRUE.
                              END IF
                              STARTH = .TRUE.
                              ENDOFH = .FALSE.
C
C  USE THE LOGICAL ARRAY SETVEC TO ENSURE THAT ALL HESSIANS ARE SET.
C
                              IF ( NHESS .GT. NSETVC ) THEN
C ** Correction 7. 21/02/94: 1 line corrected **
                                 INFORM = - 12
C ** Correction 7. 21/02/94: end of correction **
                                 GO TO 700
                              END IF
                              DO 540 I       = 1, NHESS
                                 SETVEC( I ) = .FALSE.
  540                         CONTINUE
                           END IF
C
C ---------  START H.
C
C
C  FIND WHICH COMPONENT IS TO BE SET.
C
                           DO 550 K = K1, K2
                              IVAR  = K - K1 + 1
                              IF ( FIELD2 .EQ. INAMES( K ) ) GO TO 560
  550                      CONTINUE
C
C  THE COMPONENT NAME FIELD2 IS UNRECOGNISED.
C
                           INFORM = 71
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2710 )
                           GO TO 800
  560                      CONTINUE
                           DO 570 K = K1, K2
                              JVAR  = K - K1 + 1
                              IF ( FIELD3 .EQ. INAMES( K ) ) GO TO 580
  570                      CONTINUE
C
C  THE COMPONENT NAME FIELD3 IS UNRECOGNISED.
C
                           INFORM = 71
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2710 )
                           GO TO 800
  580                      CONTINUE
C
C  FIND THE ADDRESS OF THE COMPONENT OF THE HESSIAN. THE MATRIX IS
C  STORED AS AN UPPER TRIANGLE BY ROWS.
C
                           IF ( IVAR .GT. JVAR ) THEN
                              I    = IVAR
                              IVAR = JVAR
                              JVAR = I
                           END IF
                           IHVAR = IVAR + JVAR * ( JVAR - 1 ) / 2
C
C  ENSURE THAT THE COMPONENT HAS NOT ALREADY BEEN SET.
C
                           IF ( SETVEC( IHVAR ) ) THEN
                              INFORM = 70
                              IF ( IOUT .GT. 0 ) WRITE( IOUT, 2700 )
                              GO TO 800
                           END IF
                           SETVEC( IHVAR ) = .TRUE.
                           NH              = NH + 1
                           ENDOFH          = NH .EQ. NHESS
C
C --------- SET A COMPONENT OF H.
C
                           WRITE( IOUTFN, 3160 ) FIELDI(  3 )(1:6),
     *                            FIELDI( 15 )(1:6), IHVAR, FIELD7
                        ELSE
                           IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                              IF ( STARTH ) THEN
C
C --------- CONTINUATION OF H.
C
                                 WRITE( IOUTFN, 3170 ) FIELD7
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
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      QPROD = .FALSE.
      DO 591 ITYPE = 1, MIN( 2, NELTYP ) 
         IF ( ETYPES( ITYPE ) .NE. CQSQR .AND.
     *        ETYPES( ITYPE ) .NE. CQPROD ) GO TO 930
         IF ( ETYPES( ITYPE ) .EQ. CQPROD ) QPROD = .TRUE.
  591 CONTINUE   
      IF ( IOUT .GT. 0 .AND. IPRINT .NE. 0 ) WRITE( IOUT, 2010 )
C
C  A DUMMY ROUTINE WILL BE SUBSTITUTED.
C
  600 CONTINUE 
C
C  WRITE A DUMMY ELFUNS ROUTINE.
C
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      IF ( NELTYP .EQ. 0 ) THEN
         IF ( SINGLE ) THEN
            WRITE( IOUTFN, 3003 ) FIELDI( 1 )(1:6),
     *     FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     *   ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *     FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     *     FIELDI( 12 )(1:6), FIELDI( 32 )(1:6),
     *     FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     *   ( FIELDI(  I )(1:6), I = 22, 32 ),
     *     FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     *     FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *     FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *     FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *     FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *     FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *     FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *     FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *     FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *     FIELDI( 18 )(1:6), FIELDI( 31 )(1:6), PNAME
         ELSE   
            WRITE( IOUTFN, 3002 ) FIELDI( 1 )(1:6),
     *     FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     *   ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *     FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     *     FIELDI( 12 )(1:6), FIELDI( 32 )(1:6),
     *     FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     *   ( FIELDI(  I )(1:6), I = 22, 32 ),
     *     FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     *     FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *     FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *     FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *     FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *     FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *     FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *     FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *     FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *     FIELDI( 18 )(1:6), FIELDI( 31 )(1:6), PNAME
         END IF
         WRITE( IOUTFN, 3009 ) FIELDI( 32 )(1:6)
         WRITE( IOUTFN, 3201 )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      ELSE
         IF ( SINGLE ) THEN
            WRITE( IOUTFN, 3001 ) FIELDI( 1 )(1:6),
     *     FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     *   ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *     FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     *     FIELDI( 12 )(1:6), FIELDI( 32 )(1:6), FIELDI(  5 )(1:6), 
     *     FIELDI( 12 )(1:6), ( FIELDI(  I )(1:6), I = 22, 32 ),
     *     FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     *     FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *     FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *     FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *     FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *     FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *     FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *     FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *     FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *     FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     *     PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     *     FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     *     FIELDI( 17 )(1:6), FIELDI( 20 )(1:6), FIELDI( 21 )(1:6)
            IF ( QPROD ) THEN
               WRITE( IOUTFN, 3019 ) 'X     ', 'Y     '
            ELSE
               WRITE( IOUTFN, 3019 ) 'X     '
            END IF
         ELSE   
            WRITE( IOUTFN, 3000 ) FIELDI( 1 )(1:6),
     *     FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     *   ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *     FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     *     FIELDI( 12 )(1:6), FIELDI( 32 )(1:6), FIELDI(  5 )(1:6), 
     *     FIELDI( 12 )(1:6), ( FIELDI(  I )(1:6), I = 22, 32 ),
     *     FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     *     FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *     FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *     FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *     FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *     FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *     FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *     FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *     FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *     FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     *     PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     *     FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     *     FIELDI( 17 )(1:6), FIELDI( 20 )(1:6), FIELDI( 21 )(1:6)
            IF ( QPROD ) THEN
               WRITE( IOUTFN, 3020 ) 'X     ', 'Y     '
            ELSE
               WRITE( IOUTFN, 3020 ) 'X     '
            END IF
         END IF
         WRITE( IOUTFN, 3009 ) FIELDI( 32 )(1:6)
         WRITE( IOUTFN, 3050 ) NLOOP,
     *          FIELDI( 21 )(1:6), FIELDI(  5 )(1:6),
     *          FIELDI( 13 )(1:6), FIELDI( 11 )(1:6),
     *          FIELDI( 21 )(1:6),
     *          FIELDI( 16 )(1:6), FIELDI(  7 )(1:6),
     *          FIELDI( 13 )(1:6), FIELDI( 17 )(1:6),
     *          FIELDI(  9 )(1:6), FIELDI( 13 )(1:6),
     *          FIELDI( 20 )(1:6), FIELDI( 19 )(1:6),
     *          FIELDI( 13 )(1:6),
     *          FIELDI( 12 )(1:6), FIELDI( 15 )(1:6),
     *          FIELDI( 10 )(1:6), FIELDI( 13 )(1:6)
         IF ( NELTYP .GT. 1 ) THEN
            WRITE( IOUTFN, 3051 )
     *         FIELDI( 14 )(1:6), FIELDI(  6 )(1:6),
     *         FIELDI( 13 )(1:6), ( I, I = 1, NELTYP )
            WRITE( IOUTFN, 3052 ) FIELDI( 14 )(1:6)
         END IF
C
C  MAKE SURE THAT QUADRATIC HESSIAN TERMS ARE INCLUDED.
C
         DO 640 ITYPE = 1, MIN( 2, NELTYP )
C
C  DIAGONAL TERM.
C
            IF ( ETYPES( ITYPE ) .EQ. CQSQR ) THEN
               WRITE( IOUTFN, 3060 ) ETYPES( ITYPE )
               IF ( NELTYP .GT. 1 ) WRITE( IOUTFN, 3061 ) ITYPE
               IF ( SINGLE ) THEN
                 WRITE( IOUTFN, 3053 ) 'E', 'E'
               ELSE
                 WRITE( IOUTFN, 3053 ) 'D', 'D'
               END IF
               LDEFND( ITYPE ) = .TRUE.
               ISETTY = ISETTY + 1
               IF ( ISETTY .LT. NELTYP ) WRITE( IOUTFN, 3191 ) NLOOP
C
C  OFF-DIAGONAL TERM.
C
            ELSE IF ( ETYPES( ITYPE ) .EQ. CQPROD ) THEN
               WRITE( IOUTFN, 3060 ) ETYPES( ITYPE )
               IF ( NELTYP .GT. 1 ) WRITE( IOUTFN, 3061 ) ITYPE
               IF ( SINGLE ) THEN
                 WRITE( IOUTFN, 3054 ) 'E', 'E', 'E'
               ELSE
                 WRITE( IOUTFN, 3054 ) 'D', 'D', 'D'
               END IF
               LDEFND( ITYPE ) = .TRUE.
               ISETTY = ISETTY + 1
               IF ( ISETTY .LT. NELTYP ) WRITE( IOUTFN, 3191 ) NLOOP
            END IF
  640    CONTINUE   
         WRITE( IOUTFN, 3200 ) NLOOP
      END IF
C
C  WRITE A DUMMY RANGE ROUTINE.
C
      IF ( SINGLE ) THEN
         WRITE( IOUTRA, 4003 ) PNAME
      ELSE
         WRITE( IOUTRA, 4002 ) PNAME
      END IF
      WRITE( IOUTRA, 4080 )
      WRITE( IOUTRA, 4090 )
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
      IF ( .NOT. FIRSTL ) THEN
C
C  FINISH OF THE PREVIOUS ELEMENT, IF ANY.
C
         IF ( STARTH ) THEN
            IF ( .NOT. ENDOFH ) THEN
               DO 910 IHVAR = 1, NHESS
C
C --------- SET A COMPONENT OF H.
C
                  IF ( .NOT. SETVEC( IHVAR ) ) THEN
                     IF ( SINGLE ) THEN
                        WRITE( IOUTFN, 3162 ) FIELDI(  3 )(1:6),
     *                         FIELDI( 15 )(1:6), IHVAR
                     ELSE
                        WRITE( IOUTFN, 3161 ) FIELDI(  3 )(1:6),
     *                         FIELDI( 15 )(1:6), IHVAR
                     END IF   
                  END IF   
  910          CONTINUE
               ENDOFH = .TRUE.
            END IF
C
C ---------- WIND UP H.
C
            WRITE( IOUTFN, 3180 )
         END IF
         IF ( STARTG ) THEN
C
C  SET THE REMAINING GRADIENT COMPONENTS TO ZERO.
C
            IF ( .NOT. ENDOFG ) THEN
               DO 920 IVAR = 1, NINVAR
C
C --------- SET A COMPONENT OF G.
C
                  IF ( .NOT. SETVEC( IVAR ) ) THEN
                     IF ( SINGLE ) THEN
                        WRITE( IOUTFN, 3132 )
     *                     FIELDI(  3 )(1:6),
     *                     FIELDI( 17 )(1:6), IVAR
                     ELSE
                        WRITE( IOUTFN, 3131 )
     *                     FIELDI(  3 )(1:6),
     *                     FIELDI( 17 )(1:6), IVAR
                     END IF   
                  END IF   
  920          CONTINUE
               ENDOFG = .TRUE.
            END IF
C
C ---------- WIND UP F AND G.
C
         END IF
         IF ( STARTF ) THEN
            WRITE( IOUTFN, 3190 )
         ELSE
            INFORM = 61
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
            GO TO 800
         END IF
         IF ( ISETTY .LT. NELTYP ) WRITE( IOUTFN, 3191 ) NLOOP
      END IF
C
C ---------- SUCCESSFUL RUN. WIND UP OUTPUT.
C
      INFORM = 0
      WRITE( IOUTFN, 3200 ) NLOOP
      IF ( NOINTE ) WRITE( IOUTRA, 4070 )
      IF ( NELTYP .EQ. 0 ) WRITE( IOUTRA, 4080 )
      WRITE( IOUTRA, 4090 )
C
C   CHECK THAT ALL ELEMENT TYPES HAVE BEEN DEFINED.
C
  930 CONTINUE
      DO 940 ITYPE = 1, NELTYP
         IF ( .NOT. LDEFND( ITYPE ) ) THEN
            INFORM = 68
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2680 ) ETYPES( ITYPE )
         END IF
  940 CONTINUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 1000 FORMAT( A72 )
 1010 FORMAT( A160 )
 2000 FORMAT( ' ** Exit from MAKEFN - insufficient space.',
     *        ' Increase size of ', A6 )
 2010 FORMAT( ' ** Exit from MAKEFN - warning.',
     *        ' First card not elements. ', /, '    A dummy',
     *        ' routine will be substituted ' )
 2020 FORMAT( ' ** Exit from MAKEFN - indicator card not recognised ' )
 2090 FORMAT( ' ** Exit from MAKEFN - element type not recognised ' )
 2510 FORMAT( ' ** Exit from MAKEFN -',
     *        ' name on card not that specified on input ' )
 2520 FORMAT( ' ** Exit from MAKEFN - data file incomplete.',
     *        ' No ENDATA card ' )
 2540 FORMAT( ' ** Exit from MAKEFN -',
     *        ' unrecognised field 1 in TEMPORARIES section' )
 2550 FORMAT( ' ** Exit from MAKEFN -',
     *        ' unrecognised field 1 in GLOBALS section' )
 2560 FORMAT( ' ** Exit from MAKEFN -',
     *        ' unrecognised field 1 in INDIVIDUALS section' )
 2570 FORMAT( ' ** Exit from MAKEFN -',
     *        ' undefined parameter in GLOBALS section' )
 2580 FORMAT( ' ** Exit from MAKEFN -',
     *        ' undefined parameter in INDIVIDUALS section' )
 2590 FORMAT( ' ** Exit from MAKEFN - repeated parameter name ', A8 )
 2600 FORMAT( ' ** Exit from MAKEFN - unknown component of gradient ' )
 2610 FORMAT( ' ** Exit from MAKEFN - function not set '  )
 2650 FORMAT( ' ** Exit from MAKEFN -',
     *        ' internal variable not recognised ' )
 2660 FORMAT( ' ** Exit from MAKEFN -',
     *        ' elemental variable not recognised ' )
 2670 FORMAT( ' ** Exit from MAKEFN - element type already defined ' )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
 2680 FORMAT( ' ** Exit from MAKEFN - warning, element type ', A10,
     *        ' undefined ' )
 2690 FORMAT( ' ** Exit from MAKEFN -',
     *        ' gradient component already defined ' )
 2700 FORMAT( ' ** Exit from MAKEFN -',
     *        ' Hessian component already defined ' )
 2710 FORMAT( ' ** Exit from MAKEFN - unknown component of Hessian '  )
C ** Correction 2. 19/07/93: 2 lines added **
 2720 FORMAT( ' ** Exit from MAKEFN - field 3 not blank on',
     *        ' A, F or G card ' )
C ** Correction 2. 19/07/93: end of correction **
 2900 FORMAT( ' ' )
 2970 FORMAT( ' Line ', I5, 4X, A160 )
 2980 FORMAT( ' Line ', I5, '.', I2, 1X, A65 )
 2990 FORMAT( ' Line ', I5, 4X, A65 )
 3000 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     *        '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *        '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *        '      INTEGER ', A6, /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', A6, '(', A6, ')', /,
     *        '      DOUBLE PRECISION ', A6, '(', A6, '), ',
     *                                   A6, '(', A6, '), ', 
     *                                   A6, '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *        '      INTEGER ', 5( A6, ', ' ), A6, /,
     *        '      INTEGER ', A6 )
 3001 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     *        '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *        '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *        '      INTEGER ', A6, /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', A6, '(', A6, ')', /,
     *        '      REAL             ', A6, '(', A6, '), ',
     *                                   A6, '(', A6, '), ', 
     *                                   A6, '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *        '      INTEGER ', 5( A6, ', ' ), A6, /,
     *        '      INTEGER ', A6 )
 3002 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     *        '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *        '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *        '      INTEGER ', A6, /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', A6, '(', A6, ')', /,
     *        '      DOUBLE PRECISION ', A6, '(', A6, '), ',
     *                                   A6, '(', A6, '), ', 
     *                                   A6, '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C' )
 3003 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     *        '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *        '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *        '      INTEGER ', A6, /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', A6, '(', A6, ')', /,
     *        '      REAL             ', A6, '(', A6, '), ',
     *                                   A6, '(', A6, '), ', 
     *                                   A6, '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C' )
 3009 FORMAT( '      ', A6, ' = 0' )
 3010 FORMAT( ( '      INTEGER ', A6, :, 4( ', ', A6, : ) ) )
 3019 FORMAT( ( '      REAL             ', A6, :, 4( ', ', A6, : ) ) )
 3020 FORMAT( ( '      DOUBLE PRECISION ', A6, :, 4( ', ', A6, : ) ) )
 3021 FORMAT( ( '      INTRINSIC ', A6, :, 4( ', ', A6, : ) ) )
 3022 FORMAT( ( '      EXTERNAL ', A6, :, 4( ', ', A6, : ) ) )
 3023 FORMAT( ( '      LOGICAL ', A6, 4( :, ', ', A6 ) ) )
 3030 FORMAT( '      ', A6, ' = ', A41 )
 3031 FORMAT( '      IF (', A6, ') ', A6, ' = ', A41 )
 3032 FORMAT( '      IF (.NOT.', A6, ') ', A6, ' = ', A41 )
 3040 FORMAT( '     *         ', A41 )
 3050 FORMAT( '      DO ', I5, 1X, A6, ' = 1, ', A6, /,
     *        '       ', A6, ' = ', A6, '(', A6, ') ', /,
     *        '       ', A6, ' = ', A6, '(', A6, ') - 1', /,
     *        '       ', A6, ' = ', A6, '(', A6, ') - 1', /,
     *        '       ', A6, ' = ', A6, '(', A6, ') - 1', /,
     *        '       IF ( ', A6, ' .EQ. 3 ) ',
     *                A6, ' = ', A6, '(', A6, ') - 1' )
C ** Correction 7. 15/08/95: 3 lines corrected **
 3051 FORMAT( '       ', A6, ' = ', A6, '(', A6, ')', /,
     *        '       GO TO (', 8( I5, :, ',' ), /,
     *      ( '     *        ', 8( I5, :, ',' ) ) )
C ** Correction 7. 15/08/95: end of correction **
 3052 FORMAT( '     *        ', 48X, '), ', A6 )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
 3053 FORMAT( '       X      = XVALUE(IELVAR(ILSTRT+     1))', /,
     *        '       IF ( IFFLAG .EQ. 1 ) THEN', /,
     *        '        FUVALS(IELEMN)= 5.0', A1, '-1 * X * X', /,
     *        '       ELSE', /,
     *        '        FUVALS(IGSTRT+     1)= X', /,
     *        '        IF ( IFFLAG .EQ. 3 ) THEN', /,
     *        '         FUVALS(IHSTRT+     1)= 1.0', A1, '+0', /,
     *        '        END IF', /,
     *        '       END IF' )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
 3054 FORMAT( '       X      = XVALUE(IELVAR(ILSTRT+     1))', /,
     *        '       Y      = XVALUE(IELVAR(ILSTRT+     2))', /,
     *        '       IF ( IFFLAG .EQ. 1 ) THEN', /,
     *        '        FUVALS(IELEMN)= X * Y', /,
     *        '       ELSE', /,
     *        '        FUVALS(IGSTRT+     1)= Y', /,
     *        '        FUVALS(IGSTRT+     2)= X', /,
     *        '        IF ( IFFLAG .EQ. 3 ) THEN', /,
     *        '         FUVALS(IHSTRT+     1)= 0.0', A1, '+0', /,
     *        '         FUVALS(IHSTRT+     2)= 1.0', A1, '+0', /,
     *        '         FUVALS(IHSTRT+     3)= 0.0', A1, '+0', /,
     *        '        END IF', /,
     *        '       END IF' )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
 3060 FORMAT( 'C', /, 'C  ELEMENT TYPE : ', A10, /, 'C' )
 3061 FORMAT( I5, '  CONTINUE' )
 3070 FORMAT( '       ', A6, ' = ', A6, '(', A6, '(', A6, '+', I6, '))')
 3071 FORMAT( '       ', A6, ' = ', A6, '(', A6, '+', I6, ')')
 3080 FORMAT( '       ', A6, ' = ', A41 )
 3081 FORMAT( '       IF (', A6, ') ', A6, ' = ', A41 )
C ** Correction 3. 12/01/94: 1 line corrected **
 3082 FORMAT( '       IF (.NOT.', A6, ') ', A6, '=', A41 )
C ** Correction 3. 12/01/94: end of correction **
 3083 FORMAT( '        ', A6, '       = ', A41 )
 3084 FORMAT( '        IF (', A6, ') ', A6, ' = ', A41 )
C ** Correction 4. 12/01/94: 1 line corrected **
 3085 FORMAT( '        IF (.NOT.', A6, ')', A6, '=', A41 )
C ** Correction 4. 12/01/94: end of correction **
 3090 FORMAT( '     *          ', A41 )
 3091 FORMAT( '     *           ', A41 )
 3100 FORMAT( '       IF ( ', A6, ' .EQ. 1 ) THEN', /,
     *        '        ', A6, '(', A6, ')= ', A41 )
 3110 FORMAT( '     *                  ', A41 )
 3120 FORMAT( '       ELSE' )
 3130 FORMAT( '        ', A6, '(', A6, '+', I6, ')= ', A41 )
 3131 FORMAT( '        ', A6, '(', A6, '+', I6, ')= 0.0D+0' )
 3132 FORMAT( '        ', A6, '(', A6, '+', I6, ')= 0.0E+0' )
 3140 FORMAT( '     *                         ', A41 )
 3150 FORMAT( '        IF ( ', A6, ' .EQ. 3 ) THEN' )
 3160 FORMAT( '         ', A6, '(', A6, '+', I6, ')=', A41 )
 3161 FORMAT( '         ', A6, '(', A6, '+', I6, ')=0.0D+0' )
 3162 FORMAT( '         ', A6, '(', A6, '+', I6, ')=0.0E+0' )
 3170 FORMAT( '     *                         ', A41 )
 3180 FORMAT( '        END IF' )
 3190 FORMAT( '       END IF' )
 3191 FORMAT( '       GO TO', I6 )
 3200 FORMAT( I5,  ' CONTINUE', /, '      RETURN', /,
     *        '      END' )
 3201 FORMAT( '      RETURN', /,
     *        '      END' )
 4000 FORMAT( '      SUBROUTINE RANGE( IELEMN, TRANSP, W1,',
     *        ' W2, NELVAR, NINVAR,', /,
     *        '     *                  ITYPE, LW1, LW2 )', /,
     *        '      INTEGER IELEMN, NELVAR, NINVAR, ITYPE,',
     *        ' LW1, LW2', /,
     *        '      LOGICAL TRANSP', /,
     *        '      DOUBLE PRECISION W1( LW1 ), W2( LW2 )', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *        'C  TRANSP = .FALSE. <=> W2 = U * W1', /,
     *        'C  TRANSP = .TRUE.  <=> W2 = U(TRANSPOSE) * W1', /,
     *        'C', /,
     *        '      INTEGER I' )
 4001 FORMAT( '      SUBROUTINE RANGE( IELEMN, TRANSP, W1,',
     *        ' W2, NELVAR, NINVAR,', /,
     *        '     *                  ITYPE, LW1, LW2 )', /,
     *        '      INTEGER IELEMN, NELVAR, NINVAR, ITYPE,',
     *        ' LW1, LW2', /,
     *        '      LOGICAL TRANSP', /,
     *        '      REAL             W1( LW1 ), W2( LW2 )', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *        'C  TRANSP = .FALSE. <=> W2 = U * W1', /,
     *        'C  TRANSP = .TRUE.  <=> W2 = U(TRANSPOSE) * W1', /,
     *        'C', /,
     *        '      INTEGER I' )
 4002 FORMAT( '      SUBROUTINE RANGE( IELEMN, TRANSP, W1,',
     *        ' W2, NELVAR, NINVAR,', /,
     *        '     *                  ITYPE, LW1, LW2 )', /,
     *        '      INTEGER IELEMN, NELVAR, NINVAR, ITYPE,',
     *        ' LW1, LW2', /,
     *        '      LOGICAL TRANSP', /,
     *        '      DOUBLE PRECISION W1( LW1 ), W2( LW2 )', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *        'C  TRANSP = .FALSE. <=> W2 = U * W1', /,
     *        'C  TRANSP = .TRUE.  <=> W2 = U(TRANSPOSE) * W1', /,
     *        'C' )
 4003 FORMAT( '      SUBROUTINE RANGE( IELEMN, TRANSP, W1,',
     *        ' W2, NELVAR, NINVAR,', /,
     *        '     *                  ITYPE, LW1, LW2 )', /,
     *        '      INTEGER IELEMN, NELVAR, NINVAR, ITYPE,',
     *        ' LW1, LW2', /,
     *        '      LOGICAL TRANSP', /,
     *        '      REAL             W1( LW1 ), W2( LW2 )', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *        'C  TRANSP = .FALSE. <=> W2 = U * W1', /,
     *        'C  TRANSP = .TRUE.  <=> W2 = U(TRANSPOSE) * W1', /,
     *        'C' )
C ** Correction 8. 15/08/95: 3 lines corrected **
 4040 FORMAT( '      GO TO (', 8( I5, :, ',' ), /,
     *      ( '     *       ', 8( I5, :, ',' ) ) )
C ** Correction 8. 15/08/95: end of correction **
 4050 FORMAT( '     *        ', 48X, '), ITYPE' )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
 4060 FORMAT( 'C', /, 'C  ELEMENT TYPE : ', A10, /, 'C', /,
     *        I5, ' CONTINUE', /,
     *        '      IF ( TRANSP ) THEN' )
 4070 FORMAT( 'C', /,
     *        'C  ELEMENTS WITHOUT INTERNAL VARIABLES.', /,
     *        'C', /,
     *        '99998 CONTINUE', /,
     *        '      DO 99999 I = 1, NELVAR', /,
     *        '         W2( I ) = W1( I )', /,
     *        '99999 CONTINUE', /,
     *        '      RETURN' )
 4080 FORMAT( '      RETURN' )
 4090 FORMAT( '      END' )
C
C  END OF MAKEFN.
C
      END
C     THIS VERSION: 21ST JUNE 1990.
      SUBROUTINE OUTRAN( NELV, NINV, U, IOUTFN, IOUTRA, ENAMES, INAMES,
     *                   SINGLE )
      INTEGER          NELV, NINV, IOUTFN, IOUTRA
      LOGICAL          SINGLE
      DOUBLE PRECISION U( NINV, NELV )
      CHARACTER * 10   ENAMES( * ), INAMES( * )
C
C  PRINT OUT THE GATHER AND SCATTER PART OF THE GENERATED RANGE ROUTINE
C  AND THE GATHER PART OF THE GENERATED FUNCTION EVALUATION ROUTINE.
C
      INTEGER          I, J, K
      DOUBLE PRECISION UIJ, ONE, ZERO, EPSMCH, DMACHR
      LOGICAL          ANYNNZ
      CHARACTER * 6    EVNAME, IVNAME
      INTRINSIC        DABS, MOD
      EXTERNAL         DMACHR
      DATA ZERO, ONE / 0.0D+0, 1.0D+0 /
      EPSMCH = DMACHR( 1 )
C
C  PRINT OUT THE SCATTER PART.
C
      DO 20 J    = 1, NELV
         K       = 0         
         ANYNNZ  = .FALSE.
         DO 10 I = 1, NINV
            UIJ  = U( I, J )
C
C  IGNORE ZERO ENTRIES.
C
            IF ( DABS( UIJ ) .LE. EPSMCH ) GO TO 10
            K = K + 1
            IF ( UIJ .GT. ZERO ) THEN
C
C  THE NONZERO IS POSITIVE.
C
               IF ( DABS( UIJ - ONE ) .LE. EPSMCH ) THEN
C
C  SPECIAL CASE IF NONZERO HAS THE VALUE 1.
C
                  IF ( ANYNNZ ) THEN
                     WRITE( IOUTRA, 4030 ) I
                  ELSE
                     WRITE( IOUTRA, 4040 ) J, I
                  END IF
               ELSE
C
C  NONZERO HAS A VALUE OTHER THAN 1.
C
                  IF ( ANYNNZ ) THEN
                     WRITE( IOUTRA, 4050 ) I, UIJ
                  ELSE
                     WRITE( IOUTRA, 4060 ) J, I, UIJ
                  END IF
               END IF
            ELSE
C
C  THE NONZERO IS NEGATIVE.
C
               IF ( DABS( - UIJ - ONE ) .LE. EPSMCH ) THEN
C
C  SPECIAL CASE IF NONZERO HAS THE VALUE - 1.
C
                  IF ( ANYNNZ ) THEN
                     WRITE( IOUTRA, 4070 ) I
                  ELSE
                     WRITE( IOUTRA, 4080 ) J, I
                  END IF
               ELSE
C
C  NONZERO HAS A VALUE OTHER THAN - 1.
C
                  IF ( ANYNNZ ) THEN
                     WRITE( IOUTRA, 4090 ) I, - UIJ
                  ELSE
                     WRITE( IOUTRA, 4100 ) J, I, - UIJ
                  END IF
               END IF
            END IF
            ANYNNZ = .TRUE.
            IF ( MOD( K, 19 ) .EQ. 0 ) WRITE( IOUTRA, 4112 ) J, J
   10    CONTINUE
         IF ( .NOT. ANYNNZ ) THEN
            IF ( SINGLE ) THEN
               WRITE( IOUTRA, 4111 ) J
            ELSE
               WRITE( IOUTRA, 4110 ) J
            END IF
         END IF
   20 CONTINUE
C
C  ----- THE SCATTER HAS BEEN COMPLETED; START THE GATHER.
C
      WRITE( IOUTRA, 4010 )
C
C  PRINT OUT THE GATHER PART.
C
      DO 40 I      = 1, NINV
         K         = 0
         ANYNNZ    = .FALSE.
         IVNAME    = INAMES( I )(1:6)
         DO 30 J   = 1, NELV
            EVNAME = ENAMES( J )(1:6)
            UIJ    = U( I, J )
C
C  IGNORE ZERO ENTRIES.
C
            IF ( DABS( UIJ ) .LE. EPSMCH ) GO TO 30
            K = K + 1
            IF ( UIJ .GT. ZERO ) THEN
C
C  THE NONZERO IS POSITIVE.
C
               IF ( DABS( UIJ - ONE ) .LE. EPSMCH ) THEN
C
C  SPECIAL CASE IF NONZERO HAS THE VALUE 1.
C
                  IF ( ANYNNZ ) THEN
                     WRITE( IOUTFN, 3030 ) EVNAME
                     WRITE( IOUTRA, 4030 ) J
                  ELSE
                     WRITE( IOUTFN, 3040 ) IVNAME, EVNAME
                     WRITE( IOUTRA, 4040 ) I, J
                  END IF
               ELSE
C
C  NONZERO HAS A VALUE OTHER THAN 1.
C
                  IF ( ANYNNZ ) THEN
                     WRITE( IOUTFN, 3050 ) EVNAME, UIJ
                     WRITE( IOUTRA, 4050 ) J, UIJ
                  ELSE
                     WRITE( IOUTFN, 3060 ) IVNAME, EVNAME, UIJ
                     WRITE( IOUTRA, 4060 ) I, J, UIJ
                  END IF
               END IF
             ELSE
C
C  THE NONZERO IS NEGATIVE.
C
               IF ( DABS( - UIJ - ONE ) .LE. EPSMCH ) THEN
C
C  SPECIAL CASE IF NONZERO HAS THE VALUE - 1.
C
                  IF ( ANYNNZ ) THEN
                     WRITE( IOUTFN, 3070 ) EVNAME
                     WRITE( IOUTRA, 4070 ) J
                  ELSE
                     WRITE( IOUTFN, 3080 ) IVNAME, EVNAME
                     WRITE( IOUTRA, 4080 ) I, J
                  END IF
               ELSE
C
C  NONZERO HAS A VALUE OTHER THAN - 1.
C
                  IF ( ANYNNZ ) THEN
                     WRITE( IOUTFN, 3090 ) EVNAME, - UIJ
                     WRITE( IOUTRA, 4090 ) J, - UIJ
                  ELSE
                     WRITE( IOUTFN, 3100 ) IVNAME, EVNAME, - UIJ
                     WRITE( IOUTRA, 4100 ) I, J, - UIJ
                  END IF
               END IF
            END IF
            ANYNNZ = .TRUE.
            IF ( MOD( K, 19 ) .EQ. 0 ) THEN
               WRITE( IOUTFN, 3040 ) IVNAME, IVNAME
               WRITE( IOUTRA, 4112 ) I, I
            END IF
   30    CONTINUE
         IF ( .NOT. ANYNNZ ) THEN
            IF ( SINGLE ) THEN
               WRITE( IOUTFN, 3111 ) IVNAME
               WRITE( IOUTRA, 4111 ) I
            ELSE
               WRITE( IOUTFN, 3110 ) IVNAME
               WRITE( IOUTRA, 4110 ) I
            END IF
         END IF
   40 CONTINUE
C
C  ----- THE GATHER HAS BEEN COMPLETED; WIND UP THE ELEMENT.
C
      WRITE( IOUTRA, 4020 )
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 3030 FORMAT( '     *          + ', A6 )
 3040 FORMAT( '       ', A6, ' =   ', A6 )
 3050 FORMAT( '     *          + ', A6, ' * ', F12.5 )
 3060 FORMAT( '       ', A6, ' =   ', A6, ' * ', F12.5 )
 3070 FORMAT( '     *          - ', A6 )
 3080 FORMAT( '       ', A6, ' = - ', A6 )
 3090 FORMAT( '     *          - ', A6, ' * ', F12.5 )
 3100 FORMAT( '       ', A6, ' = - ', A6, ' * ', F12.5 )
 3110 FORMAT( '       ', A6, ' = 0.0D+0 ' )
 3111 FORMAT( '       ', A6, ' = 0.0E+0 ' )
 4010 FORMAT( '      ELSE' )
 4020 FORMAT( '      END IF', /, '      RETURN' )
 4030 FORMAT( '     *                 + W1(', I6, ' ) ' )
 4040 FORMAT( '         W2(', I6, ' ) =   W1(', I6, ' ) ' )
 4050 FORMAT( '     *                 + W1(', I6, ' ) * ', F12.5 )
 4060 FORMAT( '         W2(', I6, ' ) =   W1(', I6, ' ) * ', F12.5 )
 4070 FORMAT( '     *                 - W1(', I6, ' ) ' )
 4080 FORMAT( '         W2(', I6, ' ) = - W1(', I6, ' ) ' )
 4090 FORMAT( '     *                 - W1(', I6, ' ) * ', F12.5 )
 4100 FORMAT( '         W2(', I6, ' ) = - W1(', I6, ' ) * ', F12.5 )
 4110 FORMAT( '         W2(', I6, ' ) = 0.0D+0 ' )
 4111 FORMAT( '         W2(', I6, ' ) = 0.0E+0 ' )
 4112 FORMAT( '         W2(', I6, ' ) =   W2(', I6, ' ) ' )
C
C  END OF OUTRAN.
C
      END
