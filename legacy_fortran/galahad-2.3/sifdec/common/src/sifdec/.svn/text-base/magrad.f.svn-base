C ** Correction report
C ** Correction 1. 04/04/02: 1 line corrected
C ** End of Correction report.
C  THIS VERSION: 04/04/2002 AT 09:30:00 AM.
C  ** VERSION B **

      SUBROUTINE MAGRAD( INPUT , IOUT  , IOUTGF, IOUTGD, IOUTEM, INFORM, 
     *                   NGRTYP, NGRMAX, NLMAX , NINMAX, 
     *                   PNAME , ANAMES, RENAME, INNAME, LONAME, MINAME, 
     *                   EXNAME, GTYPES, LDEFND, GPNAME, IGPA  , NGPMAX, 
     *                   DEBUG , LENGTH, ITABLE, KEY   , INLIST, SINGLE, 
     *                   NULINE, GOTLIN, IAUTO , IAD0  , IPRINT )
      INTEGER            INPUT , IOUT  , IOUTGF, INFORM, LENGTH
      INTEGER            NLMAX , NGRTYP, NINMAX, NGPMAX, NGRMAX
      INTEGER            IPRINT, IOUTGD, IOUTEM, IAUTO , IAD0
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
C  Make a group function evaluation subroutine, suitable for 
C  ---------------------------------------------------------
C  automatic differentiation from a GPS function data file
C  -------------------------------------------------------
C
C  Nick Gould 04/05/1995
C  For CGT Productions.
C
C  -------------------------------------------------------------------
C
C  Function indicator cards.
C  -------------------------
C
C  Definition   Purpose.
C  ----------   --------
C  GROUPS       Problem name.
C  TEMPORARIES  Names of additional parameters used in function defs.
C  GLOBALS      General parameter assignments.
C  INDIVIDUALS  Set function and derivative values for each group-type.
C  ENDATA       End of input data.
C
C  Data card description.
C  ----------------------
C
C  See 'A Proposal for a standard data input format for large-scale
C       nonlinear programming problems', Section 4,
C       A. R. Conn, N. I. M. Gould and Ph. L. Toint, 
C       Report CS-89-61, Dept of Computer Science, U. of Waterloo,
C       Waterloo, Ontario, N2L3G1, Canada.
C
C  -------------------------------------------------------------------
C  Returns with negative values of INFORM indicate that insufficient
C  array space has been allowed, as follows:
C
C    INFORM = - 1  when LENGTH not large enough
C    INFORM = - 2  when MAX( NINNAM, NRENAM, NLONAM, NENAM, NMINAM )
C                  .GT. NINMAX
C
      INTEGER          I, IFIELD, IFREE, ITYPE, INTYPE, IVAR, K1, K2
      INTEGER          NINNAM, NLOOP, NRENAM, NMINAM, NPNAME, NRENM1
      INTEGER          NEXNAM, NINCRS, LINENO, IIRES, K, ILINES, NLINES
      INTEGER          MBLANK, MFIXED, MFREE, MNAME, MTEMP, MGLOB, NTEM
      INTEGER          MINDIV, MENDAT, MAXNUL, MXRECL, NLONAM, NRENM2
C     INTEGER          NGTNAM
      LOGICAL          DEFNAM, ENDPAR, ENDGEN, FIRSTG
      LOGICAL          SETF, STARTP, FIXED, STARTV
      LOGICAL          ENDF, OUTGF
      PARAMETER        ( IIRES = 21 )
      PARAMETER        ( NINCRS = 2 )
      CHARACTER * 2    FIELD1
      CHARACTER * 4    AD0
      CHARACTER * 6    INCRSE( NINCRS )
      CHARACTER * 8    FIELD2, FIELD3, FIELDI( IIRES )
      CHARACTER * 10   CTEMP
      CHARACTER * 12   FIELD
      CHARACTER * 15   AORB
      CHARACTER * 41   FIELD7
      CHARACTER * 72   CTEM
      EXTERNAL         HASHB , HASHC
C
C  Parameter definitions.
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
C  Data declarations.
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
      DATA FIELDI(  1 ) / 'GROUPF  ' /,  FIELDI(  2 ) / 'GVALUE  ' /
      DATA FIELDI(  3 ) / 'LGVALU  ' /,  FIELDI(  4 ) / 'FVALUE  ' /
      DATA FIELDI(  5 ) / 'NCALCG  ' /,  FIELDI(  6 ) / 'ITYPEG  ' /
      DATA FIELDI(  7 ) / 'ICALCG  ' /,  FIELDI(  8 ) / 'DERIVS  ' /
      DATA FIELDI(  9 ) / 'IGRTYP  ' /,  FIELDI( 10 ) / 'IGROUP  ' /
      DATA FIELDI( 11 ) / 'GPVALU  ' /,  FIELDI( 12 ) / 'ISTGPA  ' /
      DATA FIELDI( 13 ) / 'IPSTRT  ' /,  FIELDI( 14 ) / 'JCALCG  ' /
      DATA FIELDI( 15 ) / 'LTYPEG  ' /,  FIELDI( 16 ) / 'LSTGPA  ' /
      DATA FIELDI( 17 ) / 'LCALCG  ' /,  FIELDI( 18 ) / 'LFVALU  ' /
      DATA FIELDI( 19 ) / 'LGPVLU  ' /,  FIELDI( 20 ) / 'IGSTAT  ' /
      DATA FIELDI( 21 ) / 'GROUP   ' /
C     DATA FIELDI( 21 ) / 'GROUPD  ' /
      IF ( IOUT .GT. 0 ) WRITE( IOUT, 2900 )
C
C  Set initial values for integer variables.
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
      NTEM   = 0
      OUTGF  = IOUTGF .GT. 0
C
C  SET INITIAL VALUES FOR LOGICAL VARIABLES.
C
      DEFNAM = .FALSE.
      ENDPAR = .FALSE.
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
C     NGTNAM = NRENAM
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
C  Insert the list of reserved integer/real/logical variables into
C  the dictionary.
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
C  -------- Set up subroutine call and reserved parameter declarations.
C
            IF ( IAUTO .EQ. 1 ) THEN
               IF ( SINGLE ) THEN
                  AORB = 'FORWARD_SINGLE '
               ELSE
                  AORB = 'FORWARD_DOUBLE '
               END IF
            ELSE
               IF ( SINGLE ) THEN
                  AORB = 'BACKWARD_SINGLE'
               ELSE
                  AORB = 'BACKWARD_DOUBLE'
               END IF
            END IF
            IF ( IAD0 .EQ. 1 ) THEN
               AD0 = 'AD01'
            ELSE
               AD0 = 'AD02'
            END IF
            IF ( SINGLE ) THEN
               IF ( OUTGF ) 
     *         WRITE( IOUTGF, 3001 ) ( FIELDI( I )( 1 : 6 ), I = 1, 4 ),
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
               WRITE( IOUTGD, 3005 ) FIELDI( 21 )( 1 : 6 ), 
     *                  ( FIELDI( I )( 1 : 6 ), I = 2, 4 ),
     *                    FIELDI( 11 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *                    FIELDI(  6 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *                    FIELDI(  7 )( 1 : 6 ), 
     *                  ( FIELDI(  I )( 1 : 6 ), I = 15, 19 ),
     *                    FIELDI(  8 )( 1 : 6 ), FIELDI( 20 )( 1 : 6 ),
     *         AD0, AORB, FIELDI(  3 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
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
               IF ( OUTGF ) 
     *         WRITE( IOUTGF, 3000 ) ( FIELDI( I )( 1 : 6 ), I = 1, 4 ),
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
               WRITE( IOUTGD, 3004 ) FIELDI( 21 )( 1 : 6 ), 
     *                  ( FIELDI( I )( 1 : 6 ), I = 2, 4 ),
     *                    FIELDI( 11 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *                    FIELDI(  6 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *                    FIELDI(  7 )( 1 : 6 ), 
     *                  ( FIELDI(  I )( 1 : 6 ), I = 15, 19 ),
     *                    FIELDI(  8 )( 1 : 6 ), FIELDI( 20 )( 1 : 6 ),
     *         AD0, AORB, FIELDI(  3 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
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
            IF ( IAD0 .EQ. 1 ) THEN
               WRITE( IOUTGD, 3006 )
            ELSE
               WRITE( IOUTGD, 3007 )
            END IF
C ** Correction 1. 04/04/02: 1 line replaced by 5
            IF ( NGRTYP .EQ. 0 ) THEN
               IF ( OUTGF ) WRITE( IOUTGF, 3009 ) FIELDI( 20 )( 1 : 6 )
               WRITE( IOUTGD, 3009 ) FIELDI( 20 )( 1 : 6 )
               GO TO 910
            END IF
            IF ( OUTGF ) WRITE( IOUTGF, 3002 ) FIELDI(  9 )( 1 : 6 ),
     *                   FIELDI( 10 )( 1 : 6 ), FIELDI( 13 )( 1 : 6 ),
     *                   FIELDI( 14 )( 1 : 6 )
            WRITE( IOUTGD, 3002 ) FIELDI(  9 )( 1 : 6 ),
     *                    FIELDI( 10 )( 1 : 6 ), FIELDI( 13 )( 1 : 6 ),
     *                    FIELDI( 14 )( 1 : 6 )
C
C --------- Insert integer declarations.
C
            IF ( NINNAM .GT. 0 .AND. OUTGF )
     *         WRITE( IOUTGF, 3010 ) ( INNAME( I ), I = 1, NINNAM )
            IF ( NINNAM .GT. 0 )
     *         WRITE( IOUTGD, 3010 ) ( INNAME( I ), I = 1, NINNAM )
C
C  Order the real values so that the list of variables which belong
C  to intrinsic or external functions follow those which do not.
C
            IF ( NRENAM .GT. 0 ) THEN
               NRENM1 = 0
               NRENM2 = NRENAM + 1
  140          CONTINUE
               IF ( NRENM1 + 1 .EQ. NRENM2 ) GO TO 180
               DO 150 I = 1, NMINAM
                  IF ( RENAME( NRENM1 + 1 ) .EQ. MINAME( I ) ) GO TO 170
  150          CONTINUE
               DO 160 I = 1, NEXNAM
                  IF ( RENAME( NRENM1 + 1 ) .EQ. EXNAME( I ) ) GO TO 170
  160          CONTINUE
               NRENM1 = NRENM1 + 1
               GO TO 140
  170          CONTINUE
               NRENM2 = NRENM2 - 1
               CTEMP = RENAME( NRENM2 )
               RENAME( NRENM2 ) = RENAME( NRENM1 + 1 )
               RENAME( NRENM1 + 1 ) = CTEMP
               GO TO 140
  180          CONTINUE
C
C --------- Insert real declarations.
C
               IF ( SINGLE ) THEN
                  IF ( OUTGF ) 
     *            WRITE( IOUTGF, 3019 ) ( RENAME( I ), I = 1, NRENAM )
               ELSE
                  IF ( OUTGF ) 
     *            WRITE( IOUTGF, 3020 ) ( RENAME( I ), I = 1, NRENAM )
               END IF
               IF ( IAD0 .EQ. 1 ) THEN
                  IF ( NRENM1 .GT. 0 ) WRITE( IOUTGD, 3018 ) 
     *                 ( RENAME( I ), I = 1, NRENM1 )
               ELSE
                  IF ( NRENM1 .GT. 0 ) WRITE( IOUTGD, 3017 ) 
     *                 ( AD0, RENAME( I ), I = 1, NRENM1 )
               END IF
               IF ( NRENM2 .LE. NRENAM ) WRITE( IOUTGD, 3017 ) 
     *              ( AD0, RENAME( I ), I = NRENM2, NRENAM )
            END IF
C
C --------- Insert logical declarations.
C
            IF ( NLONAM .GT. 0 .AND. OUTGF )
     *         WRITE( IOUTGF, 3023 ) ( LONAME( I ), I = 1, NLONAM )
            IF ( NLONAM .GT. 0 )
     *         WRITE( IOUTGD, 3023 ) ( LONAME( I ), I = 1, NLONAM )
C
C --------- Insert intrinsic declarations.
C
            IF ( NMINAM .GT. 0 .AND. OUTGF )
     *         WRITE( IOUTGF, 3021 ) ( MINAME( I ), I = 1, NMINAM )
C
C --------- Insert external declarations.
C
            IF ( NEXNAM .GT. 0 .AND. OUTGF )
     *         WRITE( IOUTGF, 3022 ) ( EXNAME( I ), I = 1, NEXNAM )
            IF ( NEXNAM .GT. 0 )
     *         WRITE( IOUTGD, 3022 ) ( EXNAME( I ), I = 1, NEXNAM )
            IF ( OUTGF ) WRITE( IOUTGF, 3009 ) FIELDI( 20 )( 1 : 6 )
            WRITE( IOUTGD, 3009 ) FIELDI( 20 )( 1 : 6 )
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
            IF ( OUTGF ) 
     *      WRITE( IOUTGF, 3050 ) NLOOP,  FIELDI( 14 )( 1 : 6 ),
     *             FIELDI(  5 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *             FIELDI(  7 )( 1 : 6 ), FIELDI( 14 )( 1 : 6 ),
     *             FIELDI(  9 )( 1 : 6 ),
     *             FIELDI(  6 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *             FIELDI(  9 )( 1 : 6 ), NLOOP,
     *             FIELDI( 13 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *             FIELDI( 10 )( 1 : 6 )
            IF ( IAD0 .EQ. 2 ) THEN
              WRITE( IOUTGD, 3011 )
C             DO I = 1, NGTNAM
              DO I = 1, NRENM1
                 WRITE( IOUTGD, 3016 ) AD0, RENAME( I )
              END DO
            END IF
            WRITE( IOUTGD, 3050 ) NLOOP,  FIELDI( 14 )( 1 : 6 ),
     *             FIELDI(  5 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *             FIELDI(  7 )( 1 : 6 ), FIELDI( 14 )( 1 : 6 ),
     *             FIELDI(  9 )( 1 : 6 ),
     *             FIELDI(  6 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *             FIELDI(  9 )( 1 : 6 ), NLOOP,
     *             FIELDI( 13 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *             FIELDI( 10 )( 1 : 6 )
            IF ( NGRTYP .GT. 1 ) THEN
               IF ( OUTGF ) WRITE( IOUTGF, 3051 ) ( I, I = 1, NGRTYP )
               WRITE( IOUTGD, 3051 ) ( I, I = 1, NGRTYP )
               IF ( OUTGF ) WRITE( IOUTGF, 3052 ) FIELDI(  9 )( 1 : 6 )
               WRITE( IOUTGD, 3052 ) FIELDI(  9 )( 1 : 6 )
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
            IF ( OUTGF ) WRITE( IOUTGF, 3030 ) FIELD2( 1 : 6 ), FIELD7
            NTEM = NTEM + 1
            WRITE( IOUTEM, 3080 ) FIELD2( 1 : 6 ), FIELD7
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
               IF ( OUTGF ) WRITE( IOUTGF, 3031 ) FIELD2( 1 : 6 ),
     *                               FIELD3( 1 : 6 ), FIELD7
               NTEM = NTEM + 1
               WRITE( IOUTEM, 3081 ) FIELD2( 1 : 6 ),
     *                               FIELD3( 1 : 6 ), FIELD7
            ELSE
               IF ( OUTGF ) WRITE( IOUTGF, 3032 ) FIELD2( 1 : 6 ),
     *                               FIELD3( 1 : 6 ), FIELD7
               NTEM = NTEM + 1
               WRITE( IOUTEM, 3082 ) FIELD2( 1 : 6 ),
     *                               FIELD3( 1 : 6 ), FIELD7
            END IF   
         END IF
      ELSE
         IF ( FIELD1( 2: 2 ) .EQ. '+' .AND. STARTP ) THEN
C
C --------- CONTINUE A PARAMETER ASSIGNMENT.
C
            IF ( OUTGF ) WRITE( IOUTGF, 3040 ) FIELD7
            NTEM = NTEM + 1
            WRITE( IOUTEM, 3040 ) FIELD7
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
C ---------- WIND UP F AND G
C
            IF ( SETF ) THEN
               IF ( .NOT. ENDF ) THEN
                  IF ( IAD0 .EQ. 1 ) THEN
                     WRITE( IOUTGD, 3122 ) FIELDI( 8 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
                  ELSE
                     WRITE( IOUTGD, 3123 ) FIELDI( 8 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
                  END IF
                  ENDF = .TRUE.
               END IF
               WRITE( IOUTGD, 3190 )
            ELSE
               INFORM = 61
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
               GO TO 800
            END IF
            IF ( ITYPE .LT. NGRTYP .AND. OUTGF ) 
     *         WRITE( IOUTGF, 3191 ) NLOOP
            IF ( ITYPE .LT. NGRTYP ) WRITE( IOUTGD, 3191 ) NLOOP
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
         IF ( OUTGF ) WRITE( IOUTGF, 3060 ) FIELD2
         WRITE( IOUTGD, 3060 ) FIELD2
         IF ( NGRTYP .GT. 1 .AND. OUTGF ) WRITE( IOUTGF, 3061 ) ITYPE
         IF ( NGRTYP .GT. 1 ) WRITE( IOUTGD, 3061 ) ITYPE
         IF ( OUTGF ) WRITE( IOUTGF, 3062 ) ANAMES( ITYPE )( 1 : 6 ),
     *             FIELDI(  4 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
         IF ( IAD0 .EQ. 1 ) THEN
            WRITE( IOUTGD, 3064 ) 
     *                FIELDI(  4 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                ANAMES( ITYPE )( 1 : 6 )
         ELSE
            WRITE( IOUTGD, 3065 ) 
     *                FIELDI(  4 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                ANAMES( ITYPE )( 1 : 6 )
         END IF
C
C --------- SET GROUP PARAMETERS.
C
         K1       = IGPA( ITYPE )
         K2       = IGPA( ITYPE + 1 ) - 1
         DO 435 K = K1, K2
            IVAR  = K - K1 + 1
            IF ( OUTGF ) 
     *      WRITE( IOUTGF, 3063 ) GPNAME( K ), FIELDI( 11 )( 1 : 6 ),
     *          FIELDI( 13 )( 1 : 6 ), IVAR
C           IF ( IAD0 .EQ. 2 ) WRITE( IOUTGD, 3015 ) AD0, GPNAME( K )
            WRITE( IOUTGD, 3063 ) GPNAME( K ), FIELDI( 11 )( 1 : 6 ),
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
         ENDF   = .TRUE.
         STARTV = .FALSE.
      ELSE
         IF ( FIELD1( 1: 1 ) .EQ. 'A' .OR. FIELD1( 1: 1 ) 
     *        .EQ. 'I' .OR. FIELD1( 1: 1 ) .EQ. 'E' ) THEN
C
C  Finish off the function assignment
C
            IF ( SETF ) THEN
               IF ( .NOT. ENDF ) THEN
                  IF ( IAD0 .EQ. 1 ) THEN
                     WRITE( IOUTGD, 3122 ) FIELDI( 8 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
                  ELSE
                     WRITE( IOUTGD, 3123 ) FIELDI( 8 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                  FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
                  END IF
                  ENDF = .TRUE.
               END IF
            END IF
C
C  START A PARAMETER ASSIGNMENT. CHECK TO SEE THAT THE PARAMETER HAS
C  BEEN DEFINED.
C
            IF ( .NOT. SETF ) THEN
               IF ( FIELD1( 2: 2 ) .EQ. ' ' ) THEN
                  STARTP = .TRUE.
C
C  Include the global parameters
C
                  IF ( .NOT. STARTV ) THEN
                     REWIND( IOUTEM )
                     DO 483 I = 1, NTEM
                        READ( IOUTEM, 1000 ) CTEM
                        WRITE( IOUTGD, 1000 ) CTEM
  483                CONTINUE   
                     STARTV = .TRUE.
                  END IF
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
                        IF ( OUTGF ) 
     *                  WRITE( IOUTGF, 3080 ) FIELD2( 1 : 6 ), FIELD7
                        WRITE( IOUTGD, 3080 ) FIELD2( 1 : 6 ), FIELD7
                     ELSE
                        IF ( OUTGF ) 
     *                  WRITE( IOUTGF, 3083 ) FIELD2( 1 : 6 ), FIELD7
                        WRITE( IOUTGD, 3083 ) FIELD2( 1 : 6 ), FIELD7
                     END IF
C
C --------- MAKE CONDITIONAL PARAMETER ASSIGNMENTS.
C
                  ELSE   
C
C  CHECK THAT THE LOGICAL VARIABLE HAS BEEN DEFINED.
C
                     FIELD = FIELD3 // '  PG'
                     CALL HASHC( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD)
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
                           IF ( OUTGF ) 
     *                     WRITE( IOUTGF, 3081 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                           WRITE( IOUTGD, 3081 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                        ELSE
                           IF ( OUTGF ) 
     *                     WRITE( IOUTGF, 3084 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                           WRITE( IOUTGD, 3084 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                        END IF
                     ELSE
                        IF ( .NOT. SETF ) THEN
                           IF ( OUTGF ) 
     *                     WRITE( IOUTGF, 3082 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                           WRITE( IOUTGD, 3082 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                        ELSE
                           IF ( OUTGF ) 
     *                     WRITE( IOUTGF, 3085 ) FIELD2( 1 : 6 ),
     *                                        FIELD3( 1 : 6 ), FIELD7
                           WRITE( IOUTGD, 3085 ) FIELD2( 1 : 6 ),
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
                           IF ( OUTGF ) 
     *                     WRITE( IOUTGF, 3090 ) FIELD7
                           WRITE( IOUTGD, 3090 ) FIELD7
                        ELSE
                           IF ( OUTGF ) 
     *                     WRITE( IOUTGF, 3091 ) FIELD7
                           WRITE( IOUTGD, 3091 ) FIELD7
                        END IF
                     ELSE
                        INFORM = 56
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
                        GO TO 800
                     END IF
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
C  Include the global parameters
C
                  IF ( .NOT. STARTV ) THEN
                     REWIND( IOUTEM )
                     DO 484 I = 1, NTEM
                        READ( IOUTEM, 1000 ) CTEM
                        WRITE( IOUTGD, 1000 ) CTEM
  484                CONTINUE   
                     STARTV = .TRUE.
                  END IF
C
C --------- START G.
C
                  IF ( OUTGF ) 
     *            WRITE( IOUTGF, 3100 ) FIELDI(  8 )( 1 : 6 ),
     *            FIELDI( 2 )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ), FIELD7
                  WRITE( IOUTGD, 3101 ) FIELD7
               ELSE
                  IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                     IF ( SETF ) THEN
C
C --------- CONTINUATION OF G.
C
                        IF ( OUTGF ) WRITE( IOUTGF, 3110 ) FIELD7
                        WRITE( IOUTGD, 3110 ) FIELD7
                     ELSE
                        INFORM = 56
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
                        GO TO 800
                     END IF
                  END IF
               END IF
            ELSE IF ( FIELD1( 1: 1 ) .EQ. 'G' .OR.
     *                FIELD1( 1: 1 ) .EQ. 'H' ) THEN
               IF ( SETF ) THEN
                  IF ( .NOT. ENDF ) THEN
                     IF ( IAD0 .EQ. 1 ) THEN
                        WRITE( IOUTGD, 3122 ) FIELDI( 8 )( 1 : 6 ),
     *                     FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                     FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                     FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
                     ELSE
                        WRITE( IOUTGD, 3123 ) FIELDI( 8 )( 1 : 6 ),
     *                     FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                     FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *                     FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
                     END IF
                     ENDF = .TRUE.
                  END IF
               END IF
            ELSE
              INFORM = 56
              IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
              GO TO 800
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
      IF ( IAUTO .EQ. 1 ) THEN
         IF ( SINGLE ) THEN
            AORB = 'FORWARD_SINGLE '
         ELSE
            AORB = 'FORWARD_DOUBLE '
         END IF
      ELSE
         IF ( SINGLE ) THEN
            AORB = 'BACKWARD_SINGLE'
         ELSE
            AORB = 'BACKWARD_DOUBLE'
         END IF
      END IF
      IF ( IAD0 .EQ. 1 ) THEN
         AD0 = 'AD01'
      ELSE
         AD0 = 'AD02'
      END IF
      IF ( SINGLE ) THEN
         IF ( OUTGF ) 
     *   WRITE( IOUTGF, 3001 ) ( FIELDI( I )( 1 : 6 ), I = 1, 4 ),
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
         WRITE( IOUTGD, 3005 ) FIELDI( 21 )( 1 : 6 ), 
     *            ( FIELDI( I )( 1 : 6 ), I = 2, 4 ),
     *              FIELDI( 11 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *              FIELDI(  6 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *              FIELDI(  7 )( 1 : 6 ), 
     *            ( FIELDI(  I )( 1 : 6 ), I = 15, 19 ),
     *              FIELDI(  8 )( 1 : 6 ), FIELDI( 20 )( 1 : 6 ),
     *   AD0, AORB, FIELDI(  3 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
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
         IF ( OUTGF ) 
     *   WRITE( IOUTGF, 3000 ) ( FIELDI( I )( 1 : 6 ), I = 1, 4 ),
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
         WRITE( IOUTGD, 3004 ) FIELDI( 21 )( 1 : 6 ), 
     *            ( FIELDI( I )( 1 : 6 ), I = 2, 4 ),
     *              FIELDI( 11 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
     *              FIELDI(  6 )( 1 : 6 ), FIELDI( 12 )( 1 : 6 ),
     *              FIELDI(  7 )( 1 : 6 ), 
     *            ( FIELDI(  I )( 1 : 6 ), I = 15, 19 ),
     *              FIELDI(  8 )( 1 : 6 ), FIELDI( 20 )( 1 : 6 ),
     *   AD0, AORB, FIELDI(  3 )( 1 : 6 ), FIELDI(  5 )( 1 : 6 ),
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
      IF ( OUTGF ) WRITE( IOUTGF, 3009 ) FIELDI( 20 )( 1 : 6 )
      WRITE( IOUTGD, 3009 ) FIELDI( 20 )( 1 : 6 )
      IF ( OUTGF ) WRITE( IOUTGF, 3210 )
      WRITE( IOUTGD, 3210 )
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
C ---------- WIND UP F AND G.
C
         IF ( SETF ) THEN
            IF ( .NOT. ENDF ) THEN
               IF ( IAD0 .EQ. 1 ) THEN
                  WRITE( IOUTGD, 3122 ) FIELDI( 8 )( 1 : 6 ),
     *               FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *               FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *               FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
               ELSE
                  WRITE( IOUTGD, 3123 ) FIELDI( 8 )( 1 : 6 ),
     *               FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *               FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 ),
     *               FIELDI( 2  )( 1 : 6 ), FIELDI( 10 )( 1 : 6 )
               END IF
               ENDF = .TRUE.
            END IF
            WRITE( IOUTGD, 3190 )
         ELSE
            INFORM = 61
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
            GO TO 800
         END IF
         IF ( ITYPE .LT. NGRTYP .AND. OUTGF ) 
     *      WRITE( IOUTGF, 3191 ) NLOOP
         IF ( ITYPE .LT. NGRTYP ) WRITE( IOUTGD, 3191 ) NLOOP
      END IF
C
C ---------- END DO LOOP.
C
      IF ( OUTGF ) WRITE( IOUTGF, 3200 ) NLOOP
      WRITE( IOUTGD, 3200 ) NLOOP
      IF ( IAD0 .EQ. 2 ) WRITE( IOUTGD, 3192 )
  910 CONTINUE
C
C ---------- SUCCESSFUL RUN. WIND UP OUTPUT.
C
      IF ( OUTGF ) WRITE( IOUTGF, 3210 )
      WRITE( IOUTGD, 3210 )
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
 2730 FORMAT( ' ** Exit from MAKEGR - field 2 or 3 not blank on',
     *        ' A, F, G or H card ' )
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
 3004 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 4( A6, ', ' ), A6, ' )', /,
     *        '      USE HSL_', A4, '_', A15, /,
     *        '      INTEGER ', 3( A6, ', ' ), A6, /,
     *        '      INTEGER ', 3( A6, ', ' ), A6, /,
     *        '      LOGICAL ', A6, /,
     *        '      INTEGER ', A6, '(', A6, '), ', A6, '(', A6, 
     *                       '), ', A6, '(', A6, ')', /,
     *        '      DOUBLE PRECISION ', A6, '(', A6, ',3), ',
     *                                   A6, '(', A6, '), ', A6, 
     *                                '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *        '      INTEGER, POINTER :: H_index( : ) ', /, 
     *        '      DOUBLE PRECISION, POINTER :: H_result( : ) ', /, 
     *        '      DOUBLE PRECISION :: A_int( 1 ) ' )
 3005 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 4( A6, ', ' ), A6, ' )', /,
     *        '      USE HSL_', A4, '_', A15, /,
     *        '      INTEGER ', 3( A6, ', ' ), A6, /,
     *        '      INTEGER ', 3( A6, ', ' ), A6, /,
     *        '      LOGICAL ', A6, /,
     *        '      INTEGER ', A6, '(', A6, '), ', A6, '(', A6, 
     *                       '), ', A6, '(', A6, ')', /,
     *        '      REAL             ', A6, '(', A6, ',3), ',
     *                                   A6, '(', A6, '), ', A6, 
     *                                '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *        '      INTEGER, POINTER :: H_index( : ) ', /, 
     *        '      REAL, POINTER :: H_result( : ) ', /, 
     *        '      REAL :: A_int( 1 ) ' )
 3006 FORMAT( '      TYPE (AD01_REAL) :: G_value = AD01_UNDEFINED', /,
     *        '      TYPE (AD01_REAL) :: A_value( 1 )' )
 3007 FORMAT( '      INTEGER :: ERROR_AD02', /, 
     *        '      TYPE (AD02_REAL) :: G_value', /,
     *        '      TYPE (AD02_REAL) :: A_value( 1 )', /,
     *        '      TYPE (AD02_DATA), POINTER :: DATA_AD02' )
 3009 FORMAT( '      ', A6, ' = 0' )
 3010 FORMAT( ( '      INTEGER ', A6, 4( :, ', ', A6 ) ) )
C3011 FORMAT( '      NULLIFY( DATA_AD02 )', /,
C    *        '      CALL AD02_INITIALIZE(2, G_value,', /,
C    *        '     *                     GVALUE(1,1),', /,
C    *        '     *                     DATA_AD02, 0)' )
 3011 FORMAT( '      CALL AD02_INITIALIZE_DATA(DATA_AD02, ERROR_AD02)' )
C3015 FORMAT( '       CALL ', A4, '_UNDEFINE( ', A6,
C    *        ', DATA_AD02 )' )
 3016 FORMAT( '      CALL ', A4, '_UNDEFINE( ', A6,
     *        ', DATA_AD02 )' )
 3017 FORMAT( ( '      TYPE (', A4, '_REAL) :: ', A6 ) )
 3018 FORMAT( ( '      TYPE (AD01_REAL) :: ', A6, 
     *          ' = AD01_UNDEFINED' ) )
 3019 FORMAT( ( '      REAL             ', A6, 4( :, ', ', A6 ) ) )
 3020 FORMAT( ( '      DOUBLE PRECISION ', A6, 4( :, ', ', A6 ) ) )
 3021 FORMAT( ( '      INTRINSIC ', A6, 4( :, ', ', A6 ) ) )
 3022 FORMAT( ( '      EXTERNAL ', A6, 4( :, ', ', A6 ) ) )
 3023 FORMAT( ( '      LOGICAL ', A6, 4( :, ', ', A6 ) ) )
 3030 FORMAT( '      ', A6, ' = ', A41 )
 3031 FORMAT( '      IF (', A6, ') ', A6, ' = ', A41 )
 3032 FORMAT( '      IF (.NOT.', A6, ') ', /, 
     *        '     *   ', A6, ' = ', A41 )
 3040 FORMAT( '     *         ', A41 )
 3050 FORMAT( '      DO ', I5, 1X, A6, ' = 1, ', A6, /,
     *        '       ', A6, ' = ', A6, '(', A6, ')', /,
     *        '       ', A6, ' = ', A6, '(', A6, ')', /,
     *        '       IF ( ', A6, ' .EQ. 0 ) GO TO ', I5, /,
     *        '       ', A6, ' = ', A6, '(', A6, ') - 1' )
 3051 FORMAT( '       GO TO (', 8( I5, :, ',' ), /,
     *      ( '     *        ', 8( I5, :, ',' ) ) )
 3052 FORMAT( '     *        ', 48X, '), ', A6 )
 3060 FORMAT( 'C', /, 'C  GROUP TYPE : ', A8, /, 'C' )
 3061 FORMAT( I5, '  CONTINUE ' )
 3062 FORMAT( '       ', A6, '= ', A6, '(', A6, ')' )
 3063 FORMAT( '       ', A6, '= ', A6, '(', A6, '+', I6, ')' )
 3064 FORMAT( '       A_int( 1 ) = ', A6, '(', A6, ')', /,
     *        '       CALL AD01_INITIALIZE(2, A_value( : 1 ),', 
     *                ' A_int( : 1 ), 0) ', /, 
     *        '       ', A6, ' = A_value( 1 ) ' ) 
 3065 FORMAT( '       A_int( 1 ) = ', A6, '(', A6, ')', /,
     *        '       CALL AD02_INITIALIZE_COMP(2, A_value( : 1 ),', 
     *                ' A_int( : 1 ),', /,
     *        '     *                      DATA_AD02, ERROR_AD02, 0)', 
     *     /, '       ', A6, ' = A_value( 1 ) ' ) 
 3080 FORMAT( '       ', A6, '= ', A41 )
 3081 FORMAT( '       IF (', A6, ') ', A6, ' = ', A41 )
 3082 FORMAT( '       IF (.NOT.', A6, ') ', A6, ' = ', A41 )
 3083 FORMAT( '        ', A6, '          = ', A41 )
 3084 FORMAT( '        IF (', A6, ') ', A6, ' = ', A41 )
 3085 FORMAT( '        IF (.NOT.', A6, ') ', A6, ' = ', A41 )
 3090 FORMAT( '     *          ', A41 )
 3091 FORMAT( '     *           ', A41 )
 3100 FORMAT( '       IF ( .NOT. ', A6, ' ) ',A6,'(', A6, ',1)= ', A41 )
 3101 FORMAT( '       G_value = ', A41 )
 3110 FORMAT( '     *                  ', A41 )
 3120 FORMAT( '       ELSE' )
 3121 FORMAT( '        WRITE(6,*) '' impossible value IFLAG = '', ', 
     *        'IFFLAG, '' in GROUPF '' ' )
 3122 FORMAT( '       IF ( ', A6, ' ) THEN', /,
     *        '        CALL AD01_HESSIAN(G_value, ', 
     *        A6, '(', A6, ',3), ', A6, '(', A6, ',2))', /,
     *        '       ELSE',/,
     *        '        CALL AD01_VALUE(G_value, ',A6,'(', A6, ',1))' )
 3123 FORMAT( '       IF ( ', A6, ' ) THEN', /,
     *        '        CALL AD02_HESSIAN(G_value, ', 
     *        A6, '(', A6, ',3), ERROR_AD02,', /, 
     *        '     *                    ', A6, '(', A6, ',2))', /,
     *        '       ELSE',/,
     *        '        CALL AD02_VALUE(G_value, ',A6,'(', A6, ',1),',
     *        ' ERROR_AD02)' )
 3130 FORMAT( '        ', A6, '(', A6, ',', I1, ')= ', A41 )
 3140 FORMAT( '     *                         ', A41 )
 3190 FORMAT( '       END IF' )
 3191 FORMAT( '       GO TO', I6 )
 3192 FORMAT( '      CALL AD02_FINALIZE_DATA(DATA_AD02, ERROR_AD02)' )
 3200 FORMAT( I5,  ' CONTINUE' )
 3210 FORMAT( '      RETURN', /,  '      END' )
C
C  END OF MAGRAD.
C
      END
