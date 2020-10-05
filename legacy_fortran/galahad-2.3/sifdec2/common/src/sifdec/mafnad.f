C  THIS VERSION: 04/08/1995 AT 11:40:35 AM.
C  ** VERSION B **

      SUBROUTINE MAFNAD( INPUT , IOUT  , IOUTFF, IOUTFD, IOUTRA, 
     *                   IOUTEM, INFORM, NLMAX , NIMAX , NETMAX, 
     *                   NINMAX, NUMAX , NEL   , NELTYP, PNAME , ENAMES, 
     *                   INAMES, RENAME, INNAME, LONAME, MINAME, EXNAME, 
     *                   ETYPES, LDEFND, LENGTH, ITABLE, KEY   , IELV  , 
     *                   IINV  , INLIST, EPNAME, IEPA  , NEPMAX, DEBUG , 
     *                   IJUMP , U     , SETVEC, NSETVC, SINGLE, 
     *                   NULINE, GOTLIN, IAUTO , IAD0  , IPRINT )
      INTEGER            INPUT , IOUT  , IOUTFF, IOUTRA, INFORM
      INTEGER            NLMAX , NIMAX , NETMAX, NELTYP, NINMAX
      INTEGER            NEPMAX, LENGTH, NSETVC, NEL   , NUMAX , IPRINT
      INTEGER            IOUTFD, IOUTEM, IAUTO , IAD0
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
C  Make a function evaluation subroutine, suitable for automatic
C  -------------------------------------------------------------
C  differentiation, and a range transformation subroutine from a 
C  -------------------------------------------------------------
C  GPS function data file.
C  -----------------------
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
C  ELEMENTS     Problem name.
C  TEMPORARIES  Names of additional parameters used in function defs.
C  GLOBALS      General parameter assignments.
C  INDIVIDUALS  Define the transformation from the elemental to the
C               internal variables for all elements with internal vars.
C               Set function and derivative values and make
C               element specific parameter assignments.
C  ENDATA       End of input data.
C
C  Data card description.
C  ----------------------
C
C  See 'A Proposal for a standard data input format for large-scale
C       nonlinear programming problems', Section 3,
C       A. R. Conn, N. I. M. Gould AND PH. L. Toint, 
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
C    INFORM = - 11 when NUMAX not large enough
C    INFORM = - 12 when nsetvc not large enough
C
C     INTEGER          JHVAR
      INTEGER          IFIELD, IFREE, IHVAR, IVAR, INTYPE, NMINAM
      INTEGER          K, K1, K2, NH, NHESS, NVARS, ISETTY
      INTEGER          I, NINAME, NINNAM, NINVAR, NLOOP, NRENAM, NEXNAM
      INTEGER          ITYPE, J, IS, JS, NOVALS, NELV, NINV, NN, NENAME
      INTEGER          NPNAME, NINCRS, LINENO, IIRES, ILINES, NLINES
      INTEGER          MBLANK, MFIXED, MFREE, MNAME, MTEMP, MGLOB
      INTEGER          MINDIV, MENDAT, MAXNUL, MXRECL, NLONAM
      INTEGER          MAXNEL, MAXNIN, NTEM  , NRENM1, NRENM2
C     INTEGER          NETNAM
      INTEGER          I1, I2, I3, I4, I5, I6
      DOUBLE PRECISION ZERO, VALUES( 2 )
      LOGICAL          NOINTE, DEFNAM, ENDPAR, ENDGEN, FIRSTL, SETRAN
      LOGICAL          STARTF, STARTP, STARTV, QPROD
      LOGICAL          ENDOFF, FIXED , CQSQRT, CQPRDT, OUTFF
      PARAMETER        ( IIRES = 33 )
      PARAMETER        ( NINCRS = 12 )
      CHARACTER * 4    AD0
      CHARACTER * 6    NUNAME
      CHARACTER * 2    FIELD1
      CHARACTER * 6    XVAR  , YVAR  , INCRSE( NINCRS )
      CHARACTER * 10   FIELD2
      CHARACTER * 8    FIELD3, FIELDS( 2 ), FIELDI( IIRES )
      CHARACTER * 10   CTEMP
      CHARACTER * 12   FIELD
      CHARACTER * 15   AORB
      CHARACTER * 41   FIELD7
      CHARACTER * 72   CTEM
      INTRINSIC        MIN
      EXTERNAL         HASHB , HASHC , GETVAL, OUTRN2, NUNAME
C
C  Parameter definitions.
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
      CHARACTER * 10   CQSQR, CQPROD
      PARAMETER      ( CQSQR  = '123456789S' )
      PARAMETER      ( CQPROD = '123456789P' )
C
C  Data declarations.
C
      DATA INCRSE / 'LENGTH', 'NINMAX', '      ', '      ', '      ',
     *              '      ', '      ', '      ', '      ', '      ',
     *              'NUMAX ', 'NSETVC'                               /
      DATA INDIC8( MBLANK ) / '            ' /, LENIND( MBLANK ) / 0  / 
      DATA INDIC8( MFIXED ) / 'FIXED FORMAT' /, LENIND( MFIXED ) / 12 /
      DATA INDIC8( MFREE  ) / 'FREE FORMAT ' /, LENIND( MFREE  ) / 11 /
      DATA INDIC8( MNAME  ) / 'ELEMENTS    ' /, LENIND( MNAME  ) / 8  /
      DATA INDIC8( MTEMP  ) / 'TEMPORARIES ' /, LENIND( MTEMP  ) / 11 / 
      DATA INDIC8( MGLOB  ) / 'GLOBALS     ' /, LENIND( MGLOB  ) / 7  /
      DATA INDIC8( MINDIV ) / 'INDIVIDUALS ' /, LENIND( MINDIV ) / 11 / 
      DATA INDIC8( MENDAT ) / 'ENDATA      ' /, LENIND( MENDAT ) / 6  /
      DATA FIELDI(  1 ) / 'ELFUNF  ' /,  FIELDI(  2 ) / 'LFUVAL  ' /
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
      DATA FIELDI( 33 ) / 'ELFUN   ' /    
C     DATA FIELDI( 33 ) / 'ELFUND  ' /    
      DATA ZERO         / 0.0D+0 /
      IF ( IOUT .GT. 0 ) WRITE( IOUT, 2900 )
C
C  Set initial values for integer variables.
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
      NTEM   = 0
C
C  Set initial values for logical variables.
C
      DEFNAM = .FALSE.
      ENDPAR = .FALSE.
      STARTP = .FALSE.
      ENDGEN = .FALSE.
      FIRSTL = .TRUE.
      NOINTE = .FALSE.
      FIXED  = .TRUE.
      GOTLIN = .FALSE.
      OUTFF = IOUTFF .GT. 0
C
C  Assign unique names to variables from quadratic terms
C
      I1 = 0
      I2 = 1
      I3 = 1
      I4 = 1
      I5 = 1
      I6 = 1
      CQSQRT = .FALSE.
      CQPRDT = .FALSE.
      DO 1 ITYPE = 1, MIN( 2, NELTYP )
         IF ( ETYPES( ITYPE ) .EQ. CQSQR ) CQSQRT = .TRUE.
         IF ( ETYPES( ITYPE ) .EQ. CQPROD ) CQPRDT = .TRUE.
    1 CONTINUE
      IF ( CQPRDT ) THEN
         YVAR = NUNAME( I1, I2, I3, I4, I5, I6, IIRES, 
     *                  NINMAX, NRENAM, NINNAM, NLONAM, 
     *                  NMINAM, NEXNAM, NLMAX , NELTYP, '      ',
     *                  FIELDI, RENAME, INNAME, LONAME, 
     *                  MINAME, EXNAME, ETYPES )
      ELSE
         YVAR = '      '
      END IF
      IF ( CQSQRT .OR.  CQPRDT ) THEN
         XVAR = NUNAME( I1, I2, I3, I4, I5, I6, IIRES, 
     *                  NINMAX, NRENAM, NINNAM, NLONAM, 
     *                  NMINAM, NEXNAM, NLMAX , NELTYP, YVAR  ,
     *                  FIELDI, RENAME, INNAME, LONAME, 
     *                  MINAME, EXNAME, ETYPES )
      ELSE
         XVAR = '      '
      END IF
C
C  Create a dictionary of the internal variable names used.
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
C  Include the names of the elemental variables used in this dictionary.
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
C     NETNAM = NRENAM
C
C  Include the names of the elemental parameters used
C  in this dictionary.
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
C  Find which element types have an internal representation.
C
      MAXNIN = 1
      MAXNEL = 1
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
         MAXNIN = MAX( MAXNIN, IINV( ITYPE + 1 ) - IINV( ITYPE ) )
         MAXNEL = MAX( MAXNEL, IELV( ITYPE + 1 ) - IELV( ITYPE ) )
   50 CONTINUE
C
C  Set a blank line.
C
      DO 60 I = 1, MXRECL
         BLNKLN( I: I ) = ' '
   60 CONTINUE    
C
C  Read next line.
C
  100 CONTINUE
      IF ( ILINES + 1 .GT. NLINES ) THEN
C
C  Read next line from the input file.
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
C  If the card is in free format, translate it into fixed format.
C
            CALL  FREEFM( NULINE, MXRECL, MENDAT, INDIC8, LENIND,
     *                    NULINA, MAXNUL, NLINES, .FALSE.,
     *                    INFORM, IOUT )
            IF ( INFORM .GT. 0 ) GO TO 800
            IF ( NLINES .GT. 0 ) THEN
C
C  If there are non-blank lines on the free format card, read the first.
C
               ILINES = 1
               NULINE = BLNKLN
               NULINE = NULINA( ILINES )
               IF ( IOUT .GT. 0 .AND. DEBUG ) WRITE( IOUT, 2980 )
     *              LINENO, ILINES, NULINE
            ELSE
C
C  There are only blank lines on the free format card.
C
               GO TO 100
            END IF
         END IF
      ELSE
C
C  Read next line from the last encountered free format card.
C
         ILINES = ILINES + 1
         NULINE = BLNKLN
         NULINE = NULINA( ILINES )
         IF ( IOUT .GT. 0 .AND. DEBUG ) WRITE( IOUT, 2980 )
     *        LINENO, ILINES, NULINE
      END IF
C
C  Consider the header part of the card.
C
      HEADER = NULINE( 1: 12 )
C
C  Ignore blank lines.
C
      IF ( HEADER .EQ. INDIC8( MBLANK ) ) GO TO 100
      IF ( NULINE( 1: 1 ) .NE. ' ' ) THEN
C
C  Ignore comment cards.
C
         IF ( NULINE( 1: 1 ) .EQ. '*' ) GO TO 100
C
C  Check if we have entered fixed-format input.
C
         IF ( HEADER .EQ. INDIC8( MFIXED ) ) THEN
            FIXED = .TRUE.
            GO TO 100
         END IF
C
C  Check if we have entered free-format input.
C
         IF ( HEADER .EQ. INDIC8( MFREE ) ) THEN
            FIXED = .FALSE.
            GO TO 100
         END IF
C
C  Check that the first encountered indicator card is the elements card.
C
         IF ( .NOT. DEFNAM  ) THEN
            IF ( HEADER .NE. INDIC8( MNAME ) ) THEN
               IF ( NELTYP .GT. 0 ) GO TO 930
               IF ( IOUT .GT. 0 .AND. IPRINT .NE. 0 ) WRITE( IOUT, 2010)
               GOTLIN = .TRUE.
               GO TO 600
            ELSE
C
C  Indicator card is elements.
C  ---------------------------
C
               IF ( PNAME  .NE. NULINE( 15: 22 ) ) THEN
                  INFORM = 51
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2510 )
                  GO TO 800
               ELSE
                  DEFNAM = .TRUE.
C
C  -------- set up subroutine call for RANGE routine.
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
C  An indicator card has been found.
C
         DO 110 I = INTYPE, MENDAT
            IF ( HEADER .EQ. INDIC8( I ) ) THEN
               INTYPE = I
               GO TO 120
            END IF
  110    CONTINUE
C
C  The indicator card is not recognised.
C
         INFORM = 2
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2020 )
         GO TO 800
  120    CONTINUE
C
C  The parameter values have been completed. Write out the
C  first part of the generated subroutine.
C
         IF ( INTYPE .GE. MGLOB .AND. .NOT. ENDPAR ) THEN
            ENDPAR = .TRUE.
C
C  Insert the list of reserved INTEGER/REAL/LOGICAL variables into
C  the dictionary.
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
               IF ( OUTFF ) WRITE( IOUTFF, 3001 ) FIELDI( 1 )(1:6),
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
               WRITE( IOUTFD, 3005 ) FIELDI( 33 )(1:6),
     *    FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     *  ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *    FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     *    FIELDI( 12 )(1:6), FIELDI( 32 )(1:6), AD0, AORB,
     *    FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     *  ( FIELDI(  I )(1:6), I = 22, 32 ), FIELDI(  6 )(1:6), 
     *    FIELDI( 22 )(1:6), FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *    FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *    FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *    FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *    FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *    FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *    FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *    FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *    FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     *    PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     *    FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     *    FIELDI( 17 )(1:6), FIELDI( 20 )(1:6),
     *    FIELDI( 21 )(1:6), MAXNIN
            ELSE   
               IF ( OUTFF ) WRITE( IOUTFF, 3000 ) FIELDI( 1 )(1:6),
     *    FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     *  ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *    FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     *    FIELDI( 12 )(1:6), FIELDI( 32 )(1:6), FIELDI(  5 )(1:6), 
     *    FIELDI( 12 )(1:6), ( FIELDI(  I )(1:6), I = 22, 32 ),
     *    FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     *    FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *    FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *    FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *    FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *    FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *    FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *    FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *    FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *    FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     *    PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     *    FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     *    FIELDI( 17 )(1:6), FIELDI( 20 )(1:6), FIELDI( 21 )(1:6)
               WRITE( IOUTFD, 3004 ) FIELDI( 33 )(1:6),
     *    FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     *  ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *    FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     *    FIELDI( 12 )(1:6), FIELDI( 32 )(1:6), AD0, AORB,
     *    FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     *  ( FIELDI(  I )(1:6), I = 22, 32 ),
     *    FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), FIELDI(  7 )(1:6), 
     *    FIELDI( 23 )(1:6), FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *    FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *    FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *    FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *    FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *    FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *    FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *    FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     *    PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     *    FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     *    FIELDI( 17 )(1:6), FIELDI( 20 )(1:6),
     *    FIELDI( 21 )(1:6), MAXNIN
            END IF
            IF ( IAD0 .EQ. 1 ) THEN
               WRITE( IOUTFD, 3024 ) MAXNEL, MAXNIN
            ELSE
               WRITE( IOUTFD, 3025 ) MAXNEL, MAXNIN
            END IF
C
C --------- Insert integer declarations.
C
            IF ( NINNAM .GT. 0 .AND. OUTFF )
     *         WRITE( IOUTFF, 3010 ) ( INNAME( I ), I = 1, NINNAM )
            IF ( NINNAM .GT. 0 )
     *         WRITE( IOUTFD, 3010 ) ( INNAME( I ), I = 1, NINNAM )
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
               IF ( OUTFF ) THEN
                  IF ( SINGLE ) THEN
                    WRITE( IOUTFF, 3019 ) ( RENAME( I ), I = 1, NRENAM )
                  ELSE
                    WRITE( IOUTFF, 3020 ) ( RENAME( I ), I = 1, NRENAM )
                  END IF
               END IF
               IF ( IAD0 .EQ. 1 ) THEN
                  IF ( NRENM1 .GT. 0 ) WRITE( IOUTFD, 3018 ) 
     *                 ( RENAME( I ), I = 1, NRENM1 )
               ELSE
                  IF ( NRENM1 .GT. 0 ) WRITE( IOUTFD, 3017 ) 
     *                 ( AD0, RENAME( I ), I = 1, NRENM1 )
               END IF
               IF ( NRENM2 .LE. NRENAM ) WRITE( IOUTFD, 3017 ) 
     *              ( AD0, RENAME( I ), I = NRENM2, NRENAM )
            END IF
C
C --------- Insert logical declarations.
C
            IF ( NLONAM .GT. 0 .AND. OUTFF )
     *         WRITE( IOUTFF, 3023 ) ( LONAME( I ), I = 1, NLONAM )
            IF ( NLONAM .GT. 0 )
     *         WRITE( IOUTFD, 3023 ) ( LONAME( I ), I = 1, NLONAM )
C
C --------- Insert intrinsic declarations.
C
            IF ( NMINAM .GT. 0 .AND. OUTFF )
     *         WRITE( IOUTFF, 3021 ) ( MINAME( I ), I = 1, NMINAM )
C
C --------- Insert external declarations.
C
            IF ( NEXNAM .GT. 0 .AND. OUTFF )
     *         WRITE( IOUTFF, 3022 ) ( EXNAME( I ), I = 1, NEXNAM )
            IF ( NEXNAM .GT. 0 )
     *         WRITE( IOUTFD, 3022 ) ( EXNAME( I ), I = 1, NEXNAM )
C
C --------- Insert variables for quadratic terms (if any)
C
            IF ( XVAR .NE. '      ' ) THEN
               IF ( YVAR .NE. '      ' ) THEN
                  IF ( SINGLE ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3019 ) XVAR, YVAR
                     WRITE( IOUTFD, 3019 ) XVAR, YVAR
                  ELSE
                     IF ( OUTFF ) WRITE( IOUTFF, 3020 ) XVAR, YVAR
                     WRITE( IOUTFD, 3020 ) XVAR, YVAR
                  END IF
               ELSE
                  IF ( SINGLE ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3019 ) XVAR
                     WRITE( IOUTFD, 3019 ) XVAR
                  ELSE
                     IF ( OUTFF ) WRITE( IOUTFF, 3020 ) XVAR
                     WRITE( IOUTFD, 3020 ) XVAR
                  END IF
               END IF
            END IF
            IF ( OUTFF ) WRITE( IOUTFF, 3009 ) FIELDI( 32 )(1:6)
            WRITE( IOUTFD, 3009 ) FIELDI( 32 )(1:6)
         END IF
C
C  The general parameter assignments have been completed.
C  Continue with the construction of the generated subroutine.
C
         IF ( INTYPE .GE. MINDIV .AND. .NOT. ENDGEN ) THEN
            ENDGEN = .TRUE.
C
C --------- Start loop over elements.
C
            IF ( OUTFF ) WRITE( IOUTFF, 3050 ) NLOOP,
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
            IF ( IAD0 .EQ. 1 ) THEN
              WRITE( IOUTFD, 3008 ) 
            ELSE
              WRITE( IOUTFD, 3011 )
C             DO I = 1, NETNAM
              DO I = 1, NRENM1
                 WRITE( IOUTFD, 3016 ) AD0, RENAME( I )
              END DO
            END IF  
            WRITE( IOUTFD, 3050 ) NLOOP,
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
               IF ( OUTFF ) WRITE( IOUTFF, 3051 )
     *            FIELDI( 14 )(1:6), FIELDI(  6 )(1:6),
     *            FIELDI( 13 )(1:6), ( I, I = 1, NELTYP )
               WRITE( IOUTFD, 3051 )
     *            FIELDI( 14 )(1:6), FIELDI(  6 )(1:6),
     *            FIELDI( 13 )(1:6), ( I, I = 1, NELTYP )
               IF ( OUTFF ) WRITE( IOUTFF, 3052 ) FIELDI( 14 )(1:6)
               WRITE( IOUTFD, 3052 ) FIELDI( 14 )(1:6)
            END IF
C
C  Make sure that quadratic Hessian terms are included
C
            DO 190 ITYPE = 1, MIN( 2, NELTYP )
C
C  Diagonal term
C
               IF ( ETYPES( ITYPE ) .EQ. CQSQR ) THEN
                  IF ( OUTFF ) WRITE( IOUTFF, 3060 ) ETYPES( ITYPE )
                  WRITE( IOUTFD, 3060 ) ETYPES( ITYPE )
                  IF ( NELTYP .GT. 1 ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3061 ) ITYPE
                     WRITE( IOUTFD, 3061 ) ITYPE
                  END IF
                  IF ( SINGLE ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3055 ) 'E'
                     WRITE( IOUTFD, 3057 ) 
     *                  XVAR, 'E', XVAR, XVAR, XVAR, 'E'
                  ELSE
                     IF ( OUTFF ) WRITE( IOUTFF, 3055 ) 'D'
                     WRITE( IOUTFD, 3057 )
     *                  XVAR, 'D', XVAR, XVAR, XVAR, 'D'
                  END IF
                  LDEFND( ITYPE ) = .TRUE.
                  ISETTY = ISETTY + 1
                  IF ( ISETTY .LT. NELTYP ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3191 ) NLOOP
                     WRITE( IOUTFD, 3191 ) NLOOP
                  END IF
C
C  Off-diagonal term
C
               ELSE IF ( ETYPES( ITYPE ) .EQ. CQPROD ) THEN
                  IF ( OUTFF ) WRITE( IOUTFF, 3060 ) ETYPES( ITYPE )
                  WRITE( IOUTFD, 3060 ) ETYPES( ITYPE )
                  IF ( NELTYP .GT. 1 ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3061 ) ITYPE
                     WRITE( IOUTFD, 3061 ) ITYPE
                  END IF
                  IF ( SINGLE ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3056 )
                     WRITE( IOUTFD, 3058 ) 
     *                 XVAR, YVAR, XVAR, YVAR, YVAR, XVAR, 'E', 'E', 'E'
                  ELSE
                     IF ( OUTFF ) WRITE( IOUTFF, 3056 )
                     WRITE( IOUTFD, 3058 )
     *                 XVAR, YVAR, XVAR, YVAR, YVAR, XVAR, 'D', 'D', 'D'
                  END IF
                  LDEFND( ITYPE ) = .TRUE.
                  ISETTY = ISETTY + 1
                  IF ( ISETTY .LT. NELTYP ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3191 ) NLOOP
                     WRITE( IOUTFD, 3191 ) NLOOP
                  END IF
               END IF
  190       CONTINUE   
         END IF
C
C  Indicator card is ENDATA.
C  -------------------------
C
         IF ( INTYPE .EQ. MENDAT ) GO TO 900
         GO TO 100
      ELSE
C
C  Check that the first non commment card is the elements indicator card.
C
         IF ( .NOT. DEFNAM  ) THEN
            IF ( NELTYP .GT. 0 ) GO TO 930
            IF ( IOUT .GT. 0 .AND. IPRINT .NE. 0 ) WRITE( IOUT, 2010 )
            GOTLIN = .TRUE.
            GO TO 600
         END IF
C
C  A data card has been found.
C  Read the character fields 1 and 2 from the card.
C
         FIELD1 = NULINE(  2:  3 )
         FIELD2 = NULINE(  5: 14 )
         IF ( INTYPE .EQ. MINDIV .AND. FIELD1 .EQ. 'R ' ) THEN
C
C  Read the character fields 3 and 5 from the card.
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
            IF ( OUTFF ) WRITE( IOUTFF, 3030 ) FIELD2(1:6), FIELD7
            NTEM = NTEM + 1
            WRITE( IOUTEM, 3080 ) FIELD2(1:6), FIELD7
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
               IF ( OUTFF ) WRITE( IOUTFF, 3031 ) FIELD2(1:6),
     *                               FIELD3(1:6), FIELD7
               NTEM = NTEM + 1
               WRITE( IOUTEM, 3081 ) FIELD2(1:6),
     *                               FIELD3(1:6), FIELD7
            ELSE
               IF ( OUTFF ) WRITE( IOUTFF, 3032 ) FIELD2(1:6),
     *                               FIELD3(1:6), FIELD7
               NTEM = NTEM + 1
               WRITE( IOUTEM, 3082 ) FIELD2(1:6),
     *                               FIELD3(1:6), FIELD7
            END IF   
         END IF
      ELSE
         IF ( FIELD1( 2: 2 ) .EQ. '+' .AND. STARTP ) THEN
C
C --------- CONTINUE A PARAMETER ASSIGNMENT.
C
           IF ( OUTFF ) WRITE( IOUTFF, 3040 ) FIELD7
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
            IF ( STARTF ) THEN
               IF ( .NOT. ENDOFF ) THEN
                  IF ( OUTFF ) THEN
                     WRITE( IOUTFF, 3120 )
                     WRITE( IOUTFF, 3121 )
                  END IF
                  IF ( IAD0 .EQ. 1 ) THEN
                     WRITE( IOUTFD, 3122 ) 
     *                  FIELDI( 12 )(1:6), FIELDI( 3  )(1:6), 
     *                  FIELDI( 13 )(1:6), FIELDI( 3  )(1:6), 
     *                  FIELDI( 17 )(1:6), FIELDI( 17 )(1:6), 
     *                  NINVAR
                  ELSE
                     WRITE( IOUTFD, 3123 ) 
     *                  FIELDI( 12 )(1:6), FIELDI( 3  )(1:6), 
     *                  FIELDI( 13 )(1:6), FIELDI( 3  )(1:6), 
     *                  FIELDI( 17 )(1:6), FIELDI( 17 )(1:6), 
     *                  NINVAR
                  END IF
                  WRITE( IOUTFD, 3150 ) FIELDI( 12 )(1:6)
                  IF ( IAD0 .EQ. 1 ) THEN
                    WRITE( IOUTFD, 3151 ) 
                  ELSE
                    WRITE( IOUTFD, 3152 ) 
                  END IF
                  DO 402 JS = 1, NINVAR
                     DO 401 IS = 1, JS
                        IHVAR = ( JS * ( JS - 1 ) ) / 2 + IS
C                       JHVAR = NINVAR * ( IS - 1 ) + JS -
C    *                          ( IS * ( IS - 1 ) ) / 2
                        IF ( IS .EQ. JS ) THEN
                           WRITE( IOUTFD, 3163 ) 
     *                        FIELDI(  3 )(1:6),
     *                        FIELDI( 15 )(1:6), IHVAR, IHVAR
C    *                        FIELDI( 15 )(1:6), IHVAR, JHVAR
                        ELSE
                           WRITE( IOUTFD, 3164 ) 
     *                        FIELDI(  3 )(1:6),
     *                        FIELDI( 15 )(1:6), IHVAR, IHVAR
C    *                        FIELDI( 15 )(1:6), IHVAR, JHVAR
                        END IF
  401                CONTINUE
  402             CONTINUE   
                  ENDOFF = .TRUE.
               END IF
               IF ( OUTFF ) WRITE( IOUTFF, 3190 )
               WRITE( IOUTFD, 3180 )
               WRITE( IOUTFD, 3190 )
            ELSE
               INFORM = 61
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
               GO TO 800
            END IF
            IF ( ISETTY .LT. NELTYP .AND. OUTFF ) 
     *         WRITE( IOUTFF, 3191 ) NLOOP
            IF ( ISETTY .LT. NELTYP ) WRITE( IOUTFD, 3191 ) NLOOP
         END IF
C
C  FIND ITYPE, THE ELEMENT TYPE.
C
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
         IF ( OUTFF ) WRITE( IOUTFF, 3060 ) FIELD2
         WRITE( IOUTFD, 3060 ) FIELD2
         IF ( NELTYP .GT. 1 .AND. OUTFF ) WRITE( IOUTFF, 3061 ) ITYPE
         IF ( NELTYP .GT. 1 ) WRITE( IOUTFD, 3061 ) ITYPE
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
         K1 = IELV( ITYPE )
         K2 = IELV( ITYPE + 1 ) - 1
         IF ( SETRAN ) THEN
            IF ( IAD0 .EQ. 1 ) THEN
               WRITE( IOUTFD, 3230 ) NELV, 
     *                FIELDI(  4 )(1:6), FIELDI(  8 )(1:6), 
     *                FIELDI( 16 )(1:6), FIELDI( 16 )(1:6), NELV
            ELSE
               WRITE( IOUTFD, 3231 ) NELV, 
     *                FIELDI(  4 )(1:6), FIELDI(  8 )(1:6), 
     *                FIELDI( 16 )(1:6), FIELDI( 16 )(1:6), NELV
            END IF
            DO 430 K = K1, K2
               IVAR  = K - K1 + 1
               IF ( OUTFF ) 
     *         WRITE( IOUTFF, 3070 ) ENAMES( K ), FIELDI(  4 )(1:6),
     *             FIELDI(  8 )(1:6), FIELDI( 16 )(1:6), IVAR
               WRITE( IOUTFD, 3220 ) ENAMES( K ), IVAR
  430       CONTINUE
         ELSE
            IF ( IAD0 .EQ. 1 ) THEN
               WRITE( IOUTFD, 3210 ) FIELDI( 12 )(1:6), NINV, 
     *                FIELDI(  4 )(1:6), FIELDI(  8 )(1:6), 
     *                FIELDI( 16 )(1:6), FIELDI( 16 )(1:6), NINV
            ELSE
               WRITE( IOUTFD, 3211 ) FIELDI( 12 )(1:6), NINV, 
     *                FIELDI(  4 )(1:6), FIELDI(  8 )(1:6), 
     *                FIELDI( 16 )(1:6), FIELDI( 16 )(1:6), NINV
            END IF
            DO 431 K = K1, K2
               IVAR  = K - K1 + 1
               IF ( OUTFF ) 
     *         WRITE( IOUTFF, 3070 ) ENAMES( K ), FIELDI(  4 )(1:6),
     *             FIELDI(  8 )(1:6), FIELDI( 16 )(1:6), IVAR
               WRITE( IOUTFD, 3220 ) ENAMES( K ), IVAR
  431       CONTINUE
         END IF
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
         ENDOFF = .TRUE.
         STARTV = .FALSE.
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
C
C  Finish the function assignment
C
               IF ( STARTF .AND. .NOT. ENDOFF ) THEN
                  IF ( OUTFF ) THEN
                     WRITE( IOUTFF, 3120 )
                     WRITE( IOUTFF, 3121 )
                  END IF
                  IF ( IAD0 .EQ. 1 ) THEN
                     WRITE( IOUTFD, 3122 ) 
     *                  FIELDI( 12 )(1:6), FIELDI( 3  )(1:6), 
     *                  FIELDI( 13 )(1:6), FIELDI( 3  )(1:6), 
     *                  FIELDI( 17 )(1:6), FIELDI( 17 )(1:6), 
     *                  NINVAR
                  ELSE
                     WRITE( IOUTFD, 3123 ) 
     *                  FIELDI( 12 )(1:6), FIELDI( 3  )(1:6), 
     *                  FIELDI( 13 )(1:6), FIELDI( 3  )(1:6), 
     *                  FIELDI( 17 )(1:6), FIELDI( 17 )(1:6), 
     *                  NINVAR
                  END IF
                  WRITE( IOUTFD, 3150 ) FIELDI( 12 )(1:6)
                  IF ( IAD0 .EQ. 1 ) THEN
                    WRITE( IOUTFD, 3151 ) 
                  ELSE
                    WRITE( IOUTFD, 3152 ) 
                  END IF
                  DO 482 JS = 1, NINVAR
                     DO 481 IS = 1, JS
                        IHVAR = ( JS * ( JS - 1 ) ) / 2 + IS
C                       JHVAR = NINVAR * ( IS - 1 ) + JS -
C    *                          ( IS * ( IS - 1 ) ) / 2
                        IF ( IS .EQ. JS ) THEN
                           WRITE( IOUTFD, 3163 ) 
     *                        FIELDI(  3 )(1:6),
     *                        FIELDI( 15 )(1:6), IHVAR, IHVAR
C    *                        FIELDI( 15 )(1:6), IHVAR, JHVAR
                        ELSE
                           WRITE( IOUTFD, 3164 ) 
     *                        FIELDI(  3 )(1:6),
     *                        FIELDI( 15 )(1:6), IHVAR, IHVAR
C    *                        FIELDI( 15 )(1:6), IHVAR, JHVAR
                        END IF
  481                CONTINUE
  482             CONTINUE   
                  ENDOFF = .TRUE.
               END IF
               IF ( .NOT. STARTF ) THEN
C
C  START A PARAMETER ASSIGNMENT.
C
                  IF ( FIELD1( 2: 2 ) .EQ. ' ' ) THEN
                     STARTP = .TRUE.
C
C  SET UP THE TRANSFORMATIONS FOR THE ELEMENT.
C
                     IF ( SETRAN ) THEN
                         CALL OUTRN2( NELV, NINV, U, IOUTFF, IOUTFD,
     *                                IOUTRA, ENAMES( JS + 1 ), 
     *                                INAMES( IS + 1 ), SINGLE, AD0 )
                         SETRAN = .FALSE.
                     END IF
C
C --------- SET ELEMENTAL PARAMETERS.
C
                     K1       = IEPA( ITYPE )
                     DO 483 K = K1, IEPA( ITYPE + 1 ) - 1
                        IVAR  = K - K1 + 1
                        IF ( OUTFF ) 
     *                  WRITE( IOUTFF, 3071 ) EPNAME( K ), 
     *                  FIELDI( 18 )(1:6), FIELDI( 20 )(1:6), IVAR
C                       IF ( IAD0 .EQ. 2 ) 
C    *                     WRITE( IOUTFD, 3015 ) AD0, EPNAME( K )
                        WRITE( IOUTFD, 3071 ) EPNAME( K ), 
     *                      FIELDI( 18 )(1:6), FIELDI( 20 )(1:6), IVAR
  483                CONTINUE
C
C  Include the global parameters
C
                     IF ( .NOT. STARTV ) THEN
                        REWIND( IOUTEM )
                        DO 484 I = 1, NTEM
                           READ( IOUTEM, 1000 ) CTEM
                           WRITE( IOUTFD, 1000 ) CTEM
  484                   CONTINUE   
                        STARTV = .TRUE.
                     END IF
C
C  CHECK TO SEE THAT THE PARAMETER HAS BEEN DEFINED.
C
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
                           IF ( OUTFF ) 
     *                     WRITE( IOUTFF, 3080 ) FIELD2(1:6), FIELD7
                           WRITE( IOUTFD, 3080 ) FIELD2(1:6), FIELD7
                        ELSE
                           IF ( OUTFF ) 
     *                     WRITE( IOUTFF, 3083 ) FIELD2(1:6), FIELD7
                           WRITE( IOUTFD, 3083 ) FIELD2(1:6), FIELD7
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
     *                               IFIELD )
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
                              IF ( OUTFF ) 
     *                        WRITE( IOUTFF, 3081 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                              WRITE( IOUTFD, 3081 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                           ELSE
                              IF ( OUTFF ) 
     *                        WRITE( IOUTFF, 3084 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                              WRITE( IOUTFD, 3084 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                           END IF
                        ELSE
                        IF ( .NOT. STARTF ) THEN
                              IF ( OUTFF ) 
     *                        WRITE( IOUTFF, 3082 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                              WRITE( IOUTFD, 3082 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                           ELSE
                              IF ( OUTFF ) 
     *                        WRITE( IOUTFF, 3085 ) FIELD2(1:6),
     *                                           FIELD3(1:6), FIELD7
                              WRITE( IOUTFD, 3085 ) FIELD2(1:6),
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
                              IF ( OUTFF ) WRITE( IOUTFF, 3090 ) FIELD7
                              WRITE( IOUTFD, 3090 ) FIELD7
                           ELSE
                              IF ( OUTFF ) WRITE( IOUTFF, 3091 ) FIELD7
                              WRITE( IOUTFD, 3091 ) FIELD7
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
                     STARTF = .TRUE.
                     ENDOFF = .FALSE.
C
C  SET UP THE TRANSFORMATIONS FOR THE ELEMENT.
C
                     IF ( SETRAN ) THEN
                         CALL OUTRN2( NELV, NINV, U, IOUTFF, IOUTFD, 
     *                                IOUTRA, ENAMES( JS + 1 ), 
     *                                INAMES( IS + 1 ), SINGLE, AD0 )
                         SETRAN = .FALSE.
                     END IF
C
C --------- SET ELEMENTAL PARAMETERS.
C
                     K1       = IEPA( ITYPE )
                     DO 485 K = K1, IEPA( ITYPE + 1 ) - 1
                        IVAR  = K - K1 + 1
                        IF ( OUTFF ) 
     *                  WRITE( IOUTFF, 3071 ) EPNAME( K ), 
     *                  FIELDI( 18 )(1:6), FIELDI( 20 )(1:6), IVAR
C                       IF ( IAD0 .EQ. 2 ) 
C    *                     WRITE( IOUTFD, 3015 ) AD0, EPNAME( K )
                        WRITE( IOUTFD, 3071 ) EPNAME( K ), 
     *                      FIELDI( 18 )(1:6), FIELDI( 20 )(1:6), IVAR
  485                CONTINUE
C
C  Include the global parameters
C
                     IF ( .NOT. STARTV ) THEN
                        REWIND( IOUTEM )
                        DO 486 I = 1, NTEM
                           READ( IOUTEM, 1000 ) CTEM
                           WRITE( IOUTFD, 1000 ) CTEM
  486                      CONTINUE   
                        STARTV = .TRUE.
                     END IF
C
C --------- START F.
C
                     IF ( OUTFF ) 
     *               WRITE( IOUTFF, 3100 ) FIELDI( 12 )(1:6),
     *               FIELDI( 3 )(1:6), FIELDI( 13 )(1:6), FIELD7
                     WRITE( IOUTFD, 3101 ) FIELD7
                  ELSE
                     IF ( FIELD1( 2: 2 ) .EQ. '+' ) THEN
                        IF ( STARTF ) THEN
C
C --------- CONTINUATION OF F.
C
                           IF ( OUTFF ) WRITE( IOUTFF, 3110 ) FIELD7
                           WRITE( IOUTFD, 3110 ) FIELD7
                        ELSE
                           INFORM = 56
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
                           GO TO 800
                        END IF
                     END IF
                  END IF
               ELSE IF ( FIELD1( 1: 1 ) .EQ. 'G' .OR.
     *                   FIELD1( 1: 1 ) .EQ. 'H' ) THEN
                  IF ( STARTF .AND. .NOT. ENDOFF ) THEN
                     IF ( OUTFF ) THEN
                        WRITE( IOUTFF, 3120 )
                        WRITE( IOUTFF, 3121 )
                     END IF
                     IF ( IAD0 .EQ. 1 ) THEN
                        WRITE( IOUTFD, 3122 ) 
     *                     FIELDI( 12 )(1:6), FIELDI( 3  )(1:6), 
     *                     FIELDI( 13 )(1:6), FIELDI( 3  )(1:6), 
     *                     FIELDI( 17 )(1:6), FIELDI( 17 )(1:6), 
     *                     NINVAR
                     ELSE
                        WRITE( IOUTFD, 3123 ) 
     *                     FIELDI( 12 )(1:6), FIELDI( 3  )(1:6), 
     *                     FIELDI( 13 )(1:6), FIELDI( 3  )(1:6), 
     *                     FIELDI( 17 )(1:6), FIELDI( 17 )(1:6), 
     *                     NINVAR
                     END IF
                     WRITE( IOUTFD, 3150 ) FIELDI( 12 )(1:6)
                     IF ( IAD0 .EQ. 1 ) THEN
                       WRITE( IOUTFD, 3151 ) 
                     ELSE
                       WRITE( IOUTFD, 3152 ) 
                     END IF
                     DO 512 JS = 1, NINVAR
                        DO 511 IS = 1, JS
                           IHVAR = ( JS * ( JS - 1 ) ) / 2 + IS
C                          JHVAR = NINVAR * ( IS - 1 ) + JS -
C    *                             ( IS * ( IS - 1 ) ) / 2
                           IF ( IS .EQ. JS ) THEN
                              WRITE( IOUTFD, 3163 ) 
     *                           FIELDI(  3 )(1:6),
     *                           FIELDI( 15 )(1:6), IHVAR, IHVAR
C    *                           FIELDI( 15 )(1:6), IHVAR, JHVAR
                           ELSE
                              WRITE( IOUTFD, 3164 ) 
     *                           FIELDI(  3 )(1:6),
     *                           FIELDI( 15 )(1:6), IHVAR, IHVAR
C    *                           FIELDI( 15 )(1:6), IHVAR, JHVAR
                           END IF
  511                   CONTINUE
  512                CONTINUE   
                     ENDOFF = .TRUE.
                  END IF
               ELSE
                  INFORM = 56
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2560 )
                  GO TO 800
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
C     IF ( NELTYP .GT. 0 ) GO TO 930
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
      IF ( NELTYP .EQ. 0 ) THEN
         IF ( SINGLE ) THEN
            IF ( OUTFF ) WRITE( IOUTFF, 3003 ) FIELDI( 1 )(1:6),
     *     FIELDI(  3 )(1:6), FIELDI(  4 )(1:6),
     *     FIELDI( 18 )(1:6), ( FIELDI(  I )(1:6), I = 5, 10 ),
     *     FIELDI( 19 )(1:6), FIELDI( 11 )(1:6), 
     *   ( FIELDI(  I )(1:6), I = 22, 31 ),
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
            WRITE( IOUTFD, 3007 ) FIELDI( 33 )(1:6),
     *     FIELDI(  3 )(1:6), FIELDI(  4 )(1:6),
     *     FIELDI( 18 )(1:6), ( FIELDI(  I )(1:6), I = 5, 10 ),
     *     FIELDI( 19 )(1:6), FIELDI( 11 )(1:6), 
     *   ( FIELDI(  I )(1:6), I = 22, 31 ), 
     *     FIELDI( 12 )(1:6), FIELDI( 32 )(1:6), AD0, AORB, 
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
            IF ( OUTFF ) WRITE( IOUTFF, 3002 ) FIELDI( 1 )(1:6),
     *    FIELDI(  3 )(1:6), FIELDI(  4 )(1:6),
     *    FIELDI( 18 )(1:6), ( FIELDI(  I )(1:6), I = 5, 10 ),
     *    FIELDI( 19 )(1:6), FIELDI( 11 )(1:6), 
     *  ( FIELDI(  I )(1:6), I = 22, 31 ),
     *    FIELDI( 12 )(1:6), FIELDI( 32 )(1:6),
     *    FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     *  ( FIELDI(  I )(1:6), I = 22, 32 ),
     *    FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     *    FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *    FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *    FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *    FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *    FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *    FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *    FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *    FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *    FIELDI( 18 )(1:6), FIELDI( 31 )(1:6), PNAME
            WRITE( IOUTFD, 3006 ) FIELDI( 33 )(1:6),
     *    FIELDI(  3 )(1:6), FIELDI(  4 )(1:6),
     *    FIELDI( 18 )(1:6), ( FIELDI(  I )(1:6), I = 5, 10 ),
     *    FIELDI( 19 )(1:6), FIELDI( 11 )(1:6), 
     *  ( FIELDI(  I )(1:6), I = 22, 31 ), 
     *    FIELDI( 12 )(1:6), FIELDI( 32 )(1:6), AD0, AORB, 
     *    FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     *  ( FIELDI(  I )(1:6), I = 22, 32 ),
     *    FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     *    FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     *    FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     *    FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     *    FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     *    FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     *    FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     *    FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     *    FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     *    FIELDI( 18 )(1:6), FIELDI( 31 )(1:6), PNAME
         END IF
         IF ( OUTFF ) WRITE( IOUTFF, 3009 ) FIELDI( 32 )(1:6)
         WRITE( IOUTFD, 3009 ) FIELDI( 32 )(1:6)
         IF ( OUTFF ) WRITE( IOUTFF, 3201 )
C        IF ( IAD0 .EQ. 2 ) THEN
C          WRITE( IOUTFD, 3203 )
C        ELSE
            WRITE( IOUTFD, 3201 )
C        END IF
      ELSE
         IF ( SINGLE ) THEN
           IF ( OUTFF ) WRITE( IOUTFF, 3001 ) FIELDI( 1 )(1:6),
     * FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     * ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     *  FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     * FIELDI( 12 )(1:6), FIELDI( 32 )(1:6), FIELDI(  5 )(1:6), 
     * FIELDI( 12 )(1:6), ( FIELDI(  I )(1:6), I = 22, 32 ),
     * FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     * FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     * FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     * FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     * FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     * FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     * FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     * FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     * FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     * FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     * PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     * FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     * FIELDI( 17 )(1:6), FIELDI( 20 )(1:6), FIELDI( 21 )(1:6)
            WRITE( IOUTFD, 3001 ) FIELDI( 33 )(1:6),
     * FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     * ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     * FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     * FIELDI( 12 )(1:6), FIELDI( 32 )(1:6),
     * FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     * ( FIELDI(  I )(1:6), I = 22, 32 ),
     * FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     * FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     * FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     * FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     * FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     * FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     * FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     * FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     * FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     * FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     * PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     * FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     * FIELDI( 17 )(1:6), FIELDI( 20 )(1:6), FIELDI( 21 )(1:6)
            IF ( QPROD ) THEN
               IF ( OUTFF ) WRITE( IOUTFF, 3019 ) 'X     ', 'Y     '
               WRITE( IOUTFD, 3019 ) 'X     ', 'Y     '
            ELSE
               IF ( OUTFF ) WRITE( IOUTFF, 3019 ) 'X     '
               WRITE( IOUTFD, 3019 ) 'X     '
            END IF
         ELSE   
            IF ( OUTFF ) WRITE( IOUTFF, 3000 ) FIELDI( 1 )(1:6),
     * FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     * ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     * FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     * FIELDI( 12 )(1:6), FIELDI( 32 )(1:6),
     * FIELDI(  5 )(1:6), FIELDI( 12 )(1:6),
     * ( FIELDI(  I )(1:6), I = 22, 32 ), FIELDI(  6 )(1:6), 
     * FIELDI( 22 )(1:6), FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     * FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     * FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     * FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     * FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     * FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     * FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     * FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     * FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     * PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     * FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     * FIELDI( 17 )(1:6), FIELDI( 20 )(1:6),
     * FIELDI( 21 )(1:6)
            WRITE( IOUTFD, 3000 ) FIELDI( 33 )(1:6),
     * FIELDI(  3 )(1:6), FIELDI(  4 )(1:6), FIELDI( 18 )(1:6),
     * ( FIELDI(  I )(1:6), I = 5, 10 ), FIELDI( 19 )(1:6), 
     * FIELDI( 11 )(1:6), ( FIELDI(  I )(1:6), I = 22, 31 ),
     * FIELDI( 12 )(1:6), FIELDI( 32 )(1:6), FIELDI(  5 )(1:6), 
     * FIELDI( 12 )(1:6), ( FIELDI(  I )(1:6), I = 22, 32 ),
     * FIELDI(  6 )(1:6), FIELDI( 22 )(1:6), 
     * FIELDI(  7 )(1:6), FIELDI( 23 )(1:6), 
     * FIELDI(  8 )(1:6), FIELDI( 24 )(1:6), 
     * FIELDI(  9 )(1:6), FIELDI( 25 )(1:6), 
     * FIELDI( 10 )(1:6), FIELDI( 26 )(1:6), 
     * FIELDI( 19 )(1:6), FIELDI( 27 )(1:6), 
     * FIELDI( 11 )(1:6), FIELDI( 28 )(1:6), 
     * FIELDI(  3 )(1:6), FIELDI( 29 )(1:6), 
     * FIELDI(  4 )(1:6), FIELDI( 30 )(1:6),
     * FIELDI( 18 )(1:6), FIELDI( 31 )(1:6),
     * PNAME, FIELDI( 13 )(1:6), FIELDI( 14 )(1:6),
     * FIELDI( 15 )(1:6), FIELDI( 16 )(1:6),
     * FIELDI( 17 )(1:6), FIELDI( 20 )(1:6), FIELDI( 21 )(1:6)
            IF ( QPROD ) THEN
               IF ( OUTFF ) WRITE( IOUTFF, 3020 ) 'X     ', 'Y     '
               WRITE( IOUTFD, 3020 ) 'X     ', 'Y     '
            ELSE
               IF ( OUTFF ) WRITE( IOUTFF, 3020 ) 'X     '
               WRITE( IOUTFD, 3020 ) 'X     '
            END IF
         END IF
         IF ( OUTFF ) WRITE( IOUTFF, 3009 ) FIELDI( 32 )(1:6)
         WRITE( IOUTFD, 3009 ) FIELDI( 32 )(1:6)
         IF ( OUTFF ) WRITE( IOUTFF, 3050 ) NLOOP,
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
         WRITE( IOUTFD, 3050 ) NLOOP,
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
            IF ( OUTFF ) WRITE( IOUTFF, 3051 )
     *         FIELDI( 14 )(1:6), FIELDI(  6 )(1:6),
     *         FIELDI( 13 )(1:6), ( I, I = 1, NELTYP )
            WRITE( IOUTFD, 3051 )
     *         FIELDI( 14 )(1:6), FIELDI(  6 )(1:6),
     *         FIELDI( 13 )(1:6), ( I, I = 1, NELTYP )
           IF ( OUTFF ) WRITE( IOUTFF, 3052 ) FIELDI( 14 )(1:6)
            WRITE( IOUTFD, 3052 ) FIELDI( 14 )(1:6)
         END IF
C
C  Make sure that quadratic Hessian terms are included
C
         DO 640 ITYPE = 1, MIN( 2, NELTYP )
C
C  Diagonal term
C
            IF ( ETYPES( ITYPE ) .EQ. CQSQR ) THEN
               IF ( OUTFF ) WRITE( IOUTFF, 3060 ) ETYPES( ITYPE )
               WRITE( IOUTFD, 3060 ) ETYPES( ITYPE )
               IF ( NELTYP .GT. 1 ) WRITE( IOUTFF, 3061 ) ITYPE
               IF ( NELTYP .GT. 1 ) WRITE( IOUTFD, 3061 ) ITYPE
               IF ( SINGLE ) THEN
                 IF ( OUTFF ) WRITE( IOUTFF, 3055 ) 'E'
                 WRITE( IOUTFD, 3053 ) 'E', 'E'
               ELSE
                 IF ( OUTFF ) WRITE( IOUTFF, 3055 ) 'D'
                 WRITE( IOUTFD, 3053 ) 'D', 'D'
               END IF
               LDEFND( ITYPE ) = .TRUE.
               ISETTY = ISETTY + 1
               IF ( ISETTY .LT. NELTYP ) THEN
                  IF ( OUTFF ) WRITE( IOUTFF, 3191 ) NLOOP
                  WRITE( IOUTFD, 3191 ) NLOOP
               END IF
C
C  Off-diagonal term
C
            ELSE IF ( ETYPES( ITYPE ) .EQ. CQPROD ) THEN
               IF ( OUTFF ) WRITE( IOUTFF, 3060 ) ETYPES( ITYPE )
               WRITE( IOUTFD, 3060 ) ETYPES( ITYPE )
               IF ( NELTYP .GT. 1 ) THEN
                  IF ( OUTFF ) WRITE( IOUTFF, 3061 ) ITYPE
                  WRITE( IOUTFD, 3061 ) ITYPE
               END IF
               IF ( SINGLE ) THEN
                 IF ( OUTFF ) WRITE( IOUTFF, 3056 )
                 WRITE( IOUTFD, 3054 ) 'E', 'E', 'E'
               ELSE
                 IF ( OUTFF ) WRITE( IOUTFF, 3056 )
                 WRITE( IOUTFD, 3054 ) 'D', 'D', 'D'
               END IF
               LDEFND( ITYPE ) = .TRUE.
               ISETTY = ISETTY + 1
               IF ( ISETTY .LT. NELTYP ) THEN
                  IF ( OUTFF ) WRITE( IOUTFF, 3191 ) NLOOP
                  WRITE( IOUTFD, 3191 ) NLOOP
               END IF
            END IF
  640    CONTINUE   
         IF ( OUTFF ) WRITE( IOUTFF, 3200 ) NLOOP
         IF ( IAD0 .EQ. 2 ) THEN
           WRITE( IOUTFD, 3202 ) NLOOP
         ELSE
           WRITE( IOUTFD, 3200 ) NLOOP
         END IF
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
         IF ( STARTF ) THEN
            IF ( .NOT. ENDOFF ) THEN
               IF ( OUTFF ) THEN
                  WRITE( IOUTFF, 3120 )
                  WRITE( IOUTFF, 3121 )
               END IF
               IF ( IAD0 .EQ. 1 ) THEN
                  WRITE( IOUTFD, 3122 ) 
     *               FIELDI( 12 )(1:6), FIELDI( 3  )(1:6), 
     *               FIELDI( 13 )(1:6), FIELDI( 3  )(1:6), 
     *               FIELDI( 17 )(1:6), FIELDI( 17 )(1:6), 
     *               NINVAR
               ELSE
                  WRITE( IOUTFD, 3123 ) 
     *               FIELDI( 12 )(1:6), FIELDI( 3  )(1:6), 
     *               FIELDI( 13 )(1:6), FIELDI( 3  )(1:6), 
     *               FIELDI( 17 )(1:6), FIELDI( 17 )(1:6), 
     *               NINVAR
               END IF
               WRITE( IOUTFD, 3150 ) FIELDI( 12 )(1:6)
               IF ( IAD0 .EQ. 1 ) THEN
                 WRITE( IOUTFD, 3151 ) 
               ELSE
                 WRITE( IOUTFD, 3152 ) 
               END IF
               DO 902 JS = 1, NINVAR
                  DO 901 IS = 1, JS
                     IHVAR = ( JS * ( JS - 1 ) ) / 2 + IS
C                    JHVAR = NINVAR * ( IS - 1 ) + JS -
C    *                       ( IS * ( IS - 1 ) ) / 2
                     IF ( IS .EQ. JS ) THEN
                        WRITE( IOUTFD, 3163 ) 
     *                     FIELDI(  3 )(1:6),
     *                     FIELDI( 15 )(1:6), IHVAR, IHVAR
C    *                     FIELDI( 15 )(1:6), IHVAR, JHVAR
                     ELSE
                        WRITE( IOUTFD, 3164 ) 
     *                     FIELDI(  3 )(1:6),
     *                     FIELDI( 15 )(1:6), IHVAR, IHVAR
C    *                     FIELDI( 15 )(1:6), IHVAR, JHVAR
                     END IF
  901             CONTINUE
  902          CONTINUE   
               ENDOFF = .TRUE.
            END IF
           IF ( OUTFF ) WRITE( IOUTFF, 3190 )
            WRITE( IOUTFD, 3180 )
            WRITE( IOUTFD, 3190 )
         ELSE
            INFORM = 61
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2610 )
            GO TO 800
         END IF
         IF ( ISETTY .LT. NELTYP .AND. OUTFF ) 
     *      WRITE( IOUTFF, 3191 ) NLOOP
         IF ( ISETTY .LT. NELTYP ) WRITE( IOUTFD, 3191 ) NLOOP
      END IF
C
C ---------- SUCCESSFUL RUN. WIND UP OUTPUT.
C
      INFORM = 0
      IF ( OUTFF ) WRITE( IOUTFF, 3200 ) NLOOP
      IF ( IAD0 .EQ. 2 ) THEN
        WRITE( IOUTFD, 3202 ) NLOOP
      ELSE
        WRITE( IOUTFD, 3200 ) NLOOP
      END IF
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
 2000 FORMAT( ' ** Exit from MAFNAD - insufficient space.',
     *        ' Increase size of ', A6 )
 2010 FORMAT( ' ** Exit from MAFNAD - warning.',
     *        ' First card not elements. ', /, '    A dummy',
     *        ' routine will be substituted ' )
 2020 FORMAT( ' ** Exit from MAFNAD - indicator card not recognised ' )
 2090 FORMAT( ' ** Exit from MAFNAD - element type not recognised ' )
 2510 FORMAT( ' ** Exit from MAFNAD -',
     *        ' name on card not that specified on input ' )
 2520 FORMAT( ' ** Exit from MAFNAD - data file incomplete.',
     *        ' No ENDATA card ' )
 2540 FORMAT( ' ** Exit from MAFNAD -',
     *        ' unrecognised field 1 in TEMPORARIES section' )
 2550 FORMAT( ' ** Exit from MAFNAD -',
     *        ' unrecognised field 1 in GLOBALS section' )
 2560 FORMAT( ' ** Exit from MAFNAD -',
     *        ' unrecognised field 1 in INDIVIDUALS section' )
 2570 FORMAT( ' ** Exit from MAFNAD -',
     *        ' undefined parameter in GLOBALS section' )
 2580 FORMAT( ' ** Exit from MAFNAD -',
     *        ' undefined parameter in INDIVIDUALS section' )
 2590 FORMAT( ' ** Exit from MAFNAD - repeated parameter name ', A8 )
 2600 FORMAT( ' ** Exit from MAFNAD - unknown component of gradient ' )
 2610 FORMAT( ' ** Exit from MAFNAD - function not set '  )
 2650 FORMAT( ' ** Exit from MAFNAD -',
     *        ' internal variable not recognised ' )
 2660 FORMAT( ' ** Exit from MAFNAD -',
     *        ' elemental variable not recognised ' )
 2670 FORMAT( ' ** Exit from MAFNAD - element type already defined ' )
 2680 FORMAT( ' ** Exit from MAFNAD - warning, element type ', A10,
     *        ' undefined ' )
 2690 FORMAT( ' ** Exit from MAFNAD -',
     *        ' gradient component already defined ' )
 2700 FORMAT( ' ** Exit from MAFNAD -',
     *        ' Hessian component already defined ' )
 2710 FORMAT( ' ** Exit from MAFNAD - unknown component of Hessian '  )
 2720 FORMAT( ' ** Exit from MAFNAD - field 3 not blank on',
     *        ' A, F or G card ' )
 2900 FORMAT( ' ' )
 2970 FORMAT( ' Line ', I5, 4X, A160 )
 2980 FORMAT( ' Line ', I5, '.', I2, 1X, A65 )
 2990 FORMAT( ' Line ', I5, 4X, A65 )
 3000 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *  '     *                   ', 5( A6, ', ' ), /,
     *  '     *                   ', 5( A6, ', ' ), /,
     *  '     *                   ', 5( A6, ', ' ), /,
     *  '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     *  '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *  '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *  '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                       A6, '(', A6, ')', /,
     *  '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                       A6, '(', A6, ')', /,
     *  '      INTEGER ', A6, '(', A6, ')', /,
     *  '      DOUBLE PRECISION ', A6, '(', A6, '), ',
     *                             A6, '(', A6, '), ', 
     *                             A6, '(', A6, ')', /,
     *  'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *  '      INTEGER ', 5( A6, ', ' ), A6, /,
     *  '      INTEGER ', A6 )
 3001 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *  '     *                   ', 5( A6, ', ' ), /,
     *  '     *                   ', 5( A6, ', ' ), /,
     *  '     *                   ', 5( A6, ', ' ), /,
     *  '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     *  '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *  '      INTEGER ', 5( A6, ', ' ),  A6, /,
     *  '      INTEGER ', A6, /,
     *  '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                       A6, '(', A6, ')', /,
     *  '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                       A6, '(', A6, ')', /,
     *  '      INTEGER ', A6, '(', A6, ')', /,
     *  '      REAL             ', A6, '(', A6, '), ',
     *                             A6, '(', A6, '), ', 
     *                             A6, '(', A6, ')', /,
     *  'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     *  '      INTEGER ', 5( A6, ', ' ), A6, /,
     *  '      INTEGER ', A6 )
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
 3004 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     * '      USE HSL_', A4, '_', A15, /, '      INTEGER ', 
     *  5( A6, ', ' ),  A6, /, '      INTEGER ', 5( A6, ', ' ),  A6, /,
     * '      INTEGER ', A6, /, '      INTEGER ', 2( A6, '(', 
     * A6, '), ' ), A6, '(', A6, ')', /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', A6, '(', A6, ')', /,
     *        '      DOUBLE PRECISION ', A6, '(', A6, '), ',
     *                                   A6, '(', A6, '), ', 
     *                                   A6, '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     * '      INTEGER ', 5( A6, ', ' ), A6, /, '      INTEGER ', A6, /, 
     *        '      INTEGER, POINTER :: H_index( : ) ', /, 
     *        '      DOUBLE PRECISION, POINTER :: H_result( : ) ', /, 
     *        '      DOUBLE PRECISION X_int( ', I6, ') ' )
 3005 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     * '      USE HSL_', A4, '_', A15, /, '      INTEGER ', 5( A6, 
     * ', ' ),  A6, /, '      INTEGER ', 5( A6, ', ' ),  A6, /,
     * '      INTEGER ', A6, /, '      INTEGER ', 2( A6, '(', A6, 
     * '), ' ), A6, '(', A6, ')', /,
     *        '      INTEGER ', 2( A6, '(', A6, '), ' ), 
     *                             A6, '(', A6, ')', /,
     *        '      INTEGER ', A6, '(', A6, ')', /,
     *        '      REAL             ', A6, '(', A6, '), ',
     *                                   A6, '(', A6, '), ', 
     *                                   A6, '(', A6, ')', /,
     *        'C', /, 'C  PROBLEM NAME : ', A8, /, 'C', /,
     * '      INTEGER ', 5( A6, ', ' ), A6, /, '      INTEGER ', A6, /,
     *        '      INTEGER, POINTER :: H_index( : ) ', /, 
     *        '      REAL, POINTER :: H_result( : ) ', /, 
     *        '      REAL X_int( ', I6, ') ' )
 3006 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     *        '      USE HSL_', A4, '_', A15, /,
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
 3007 FORMAT( '      SUBROUTINE ', A6, '( ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 5( A6, ', ' ), /,
     *        '     *                   ', 2( A6, ', ' ), A6, ' )', /,
     *        '      USE HSL_', A4, '_', A15, /,
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
 3008 FORMAT( '      X_AD01_int = AD01_UNDEFINED' )
C3008 FORMAT( '      CALL AD01_UNDEFINE( X_AD01_int ) ' )
 3009 FORMAT( '      ', A6, ' = 0' )
 3010 FORMAT( ( '      INTEGER ', A6, :, 4( ', ', A6, : ) ) )
C3011 FORMAT( '      NULLIFY( DATA_AD02 )', /,
C    *        '      CALL AD02_INITIALIZE(IFFLAG-1, X_value(1),', /,
C    *        '     *          XVALUE(IELVAR(ISTAEV(1)+1)),', /,
C    *        '     *                      DATA_AD02, 0)' )
 3011 FORMAT( '      CALL AD02_INITIALIZE_DATA(DATA_AD02, ERROR_AD02)' )
C3015 FORMAT( '       CALL ', A4, '_UNDEFINE( ', A6,
C    *        ', DATA_AD02 )' )
 3016 FORMAT( '      CALL ', A4, '_UNDEFINE( ', A6,
     *        ', DATA_AD02 )' )
 3017 FORMAT( ( '      TYPE (', A4, '_REAL) :: ', A6 ) )
 3018 FORMAT( ( '      TYPE (AD01_REAL) :: ', A6, 
     *          ' = AD01_UNDEFINED' ) )
 3019 FORMAT( ( '      REAL             ', A6, :, 4( ', ', A6, : ) ) )
 3020 FORMAT( ( '      DOUBLE PRECISION ', A6, :, 4( ', ', A6, : ) ) )
 3021 FORMAT( ( '      INTRINSIC ', A6, :, 4( ', ', A6, : ) ) )
 3022 FORMAT( ( '      EXTERNAL ', A6, :, 4( ', ', A6, : ) ) )
 3023 FORMAT( ( '      LOGICAL ', A6, 4( :, ', ', A6 ) ) )
 3024 FORMAT( '      TYPE (AD01_REAL) :: F_value = AD01_UNDEFINED', /,
     *        '      TYPE (AD01_REAL) :: X_value(', I6, '),', 
     *             ' X_AD01_int(',I6, ')' )
 3025 FORMAT( '      INTEGER :: ERROR_AD02', /, 
     *        '      TYPE (AD02_REAL) :: F_value', /,
     *        '      TYPE (AD02_REAL) :: X_value(', I6, '),', 
     *             ' X_AD02_int(',I6, ')', /, 
     *        '      TYPE (AD02_DATA), POINTER :: DATA_AD02' )
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
 3051 FORMAT( '       ', A6, ' = ', A6, '(', A6, ')', /,
     *        '       GO TO (', 8( I5, :, ',' ), /,
     *      ( '     *        ', 8( I5, :, ',' ) ) )
 3052 FORMAT( '     *        ', 48X, '), ', A6 )
 3053 FORMAT( '       X      = XVALUE(IELVAR(ILSTRT+     1))', /,
     *        '       IF ( IFFLAG .EQ. 1 ) THEN', /,
     *        '        FUVALS(IELEMN)= 5.0', A1, '-1 * X * X', /,
     *        '       ELSE', /,
     *        '        FUVALS(IGSTRT+     1)= X', /,
     *        '        IF ( IFFLAG .EQ. 3 ) THEN', /,
     *        '         FUVALS(IHSTRT+     1)= 1.0', A1, '+0', /,
     *        '        END IF', /,
     *        '       END IF' )
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
 3055 FORMAT( '       X      = XVALUE(IELVAR(ILSTRT+     1))', /,
     *        '       IF ( IFFLAG .EQ. 1 ) THEN', /,
     *        '        FUVALS(IELEMN)= 5.0', A1, '-1 * X * X', /,
     *        '       ELSE', /,
     *        '        WRITE(6,*) '' impossible value IFFLAG = '', ', 
     *        'IFFLAG, '' in ELFUNF '' ', /,
     *        '       END IF' )
 3056 FORMAT( '       X      = XVALUE(IELVAR(ILSTRT+     1))', /,
     *        '       Y      = XVALUE(IELVAR(ILSTRT+     2))', /,
     *        '       IF ( IFFLAG .EQ. 1 ) THEN', /,
     *        '        FUVALS(IELEMN)= X * Y', /,
     *        '       ELSE', /,
     *        '        WRITE(6,*) '' impossible value IFFLAG = '', ', 
     *        'IFFLAG, '' in ELFUNF '' ', /,
     *        '       END IF' )
 3057 FORMAT( '       ', A6, ' = XVALUE(IELVAR(ILSTRT+     1))', /,
     *        '       IF ( IFFLAG .EQ. 1 ) THEN', /,
     *        '        FUVALS(IELEMN)= 5.0', A1, '-1 * ', A6, 
     *        ' * ', A6, /,
     *        '       ELSE', /,
     *        '        FUVALS(IGSTRT+     1)= ', A6, /,
     *        '        IF ( IFFLAG .EQ. 3 ) THEN', /,
     *        '         FUVALS(IHSTRT+     1)= 1.0', A1, '+0', /,
     *        '        END IF', /,
     *        '       END IF' )
 3058 FORMAT( '       ', A6, ' = XVALUE(IELVAR(ILSTRT+     1))', /,
     *        '       ', A6, ' = XVALUE(IELVAR(ILSTRT+     2))', /,
     *        '       IF ( IFFLAG .EQ. 1 ) THEN', /,
     *        '        FUVALS(IELEMN)= ', A6, ' * ', A6, /,
     *        '       ELSE', /,
     *        '        FUVALS(IGSTRT+     1)= ', A6, /,
     *        '        FUVALS(IGSTRT+     2)= ', A6, /,
     *        '        IF ( IFFLAG .EQ. 3 ) THEN', /,
     *        '         FUVALS(IHSTRT+     1)= 0.0', A1, '+0', /,
     *        '         FUVALS(IHSTRT+     2)= 1.0', A1, '+0', /,
     *        '         FUVALS(IHSTRT+     3)= 0.0', A1, '+0', /,
     *        '        END IF', /,
     *        '       END IF' )
 3060 FORMAT( 'C', /, 'C  ELEMENT TYPE : ', A10, /, 'C' )
 3061 FORMAT( I5, '  CONTINUE' )
 3070 FORMAT( '       ', A6, ' = ', A6, '(', A6, '(', A6, '+', I6, '))')
 3071 FORMAT( '       ', A6, ' = ', A6, '(', A6, '+', I6, ')')
 3080 FORMAT( '       ', A6, ' = ', A41 )
 3081 FORMAT( '       IF (', A6, ') ', A6, ' = ', A41 )
 3082 FORMAT( '       IF (.NOT.', A6, ') ', A6, '=', A41 )
 3083 FORMAT( '        ', A6, '       = ', A41 )
 3084 FORMAT( '        IF (', A6, ') ', A6, ' = ', A41 )
 3085 FORMAT( '        IF (.NOT.', A6, ')', A6, '=', A41 )
 3090 FORMAT( '     *          ', A41 )
 3091 FORMAT( '     *           ', A41 )
 3100 FORMAT( '       IF ( ', A6, ' .EQ. 1 ) THEN', /,
     *        '        ', A6, '(', A6, ')= ', A41 )
 3101 FORMAT( '       F_value = ', A41 )
 3110 FORMAT( '     *                  ', A41 )
 3120 FORMAT( '       ELSE' )
 3121 FORMAT( '        WRITE(6,*) '' impossible value IFFLAG = '', ', 
     *        'IFFLAG, '' in ELFUNF '' ' )
 3122 FORMAT( '       IF ( ', A6, ' .EQ. 1 ) THEN', /,
     *        '        CALL AD01_VALUE(F_value, ', A6, '(', A6, '))', /,
     *        '       ELSE',/,
     *        '        CALL AD01_GRAD(F_value, ', A6, '(', A6, '+1:',
     *        A6, '+', I6, '))' )
 3123 FORMAT( '       IF ( ', A6, ' .EQ. 1 ) THEN', /,
     *        '        CALL AD02_VALUE(F_value, ', A6, '(', A6, '),', 
     *        ' ERROR_AD02)', /,
     *        '       ELSE',/,
     *        '        CALL AD02_GRAD(F_value, ', A6, '(', A6, '+1:',
     *        A6, '+', I6, '),', /, 
     *        '     *                 ERROR_AD02)' )
 3130 FORMAT( '        ', A6, '(', A6, '+', I6, ')= ', A41 )
 3131 FORMAT( '        ', A6, '(', A6, '+', I6, ')= 0.0D+0' )
 3132 FORMAT( '        ', A6, '(', A6, '+', I6, ')= 0.0E+0' )
 3140 FORMAT( '     *                         ', A41 )
 3150 FORMAT( '        IF ( ', A6, ' .EQ. 3 ) THEN' )
 3151 FORMAT( '         CALL AD01_DERIVS(F_value, 2,',
     *        ' H_index, H_result)' )
 3152 FORMAT( '         CALL AD02_DERIVS(F_value, 2,',
     *        ' H_index, H_result, ERROR_AD02)' )
 3160 FORMAT( '         ', A6, '(', A6, '+', I6, ')=', A41 )
 3161 FORMAT( '         ', A6, '(', A6, '+', I6, ')=0.0D+0' )
 3162 FORMAT( '         ', A6, '(', A6, '+', I6, ')=0.0E+0' )
 3163 FORMAT( '         ', A6, '(', A6, '+', I6, ')=2.0*H_result(', 
     *        I6, ')')
 3164 FORMAT( '         ', A6, '(', A6, '+', I6, ')=H_result(', I6, ')')
 3170 FORMAT( '     *                         ', A41 )
 3180 FORMAT( '        END IF' )
 3190 FORMAT( '       END IF' )
 3191 FORMAT( '       GO TO', I6 )
 3200 FORMAT( I5,  ' CONTINUE', /, '      RETURN', /,
     *        '      END' )
 3201 FORMAT( '      RETURN', /,
     *        '      END' )
 3202 FORMAT( I5,  ' CONTINUE', /, 
     *        '      CALL AD02_FINALIZE_DATA(DATA_AD02, ERROR_AD02)', /,
     *        '      RETURN', /,   '      END' )
 3203 FORMAT( '      CALL AD02_FINALIZE_DATA(DATA_AD02, ERROR_AD02)', /, 
     *        '      RETURN', /, '      END' )
 3210 FORMAT( '       CALL AD01_INITIALIZE(',A6,' - 1, X_value(:', I6, 
     *               '),', /, '     *', 22X, 
     *               A6,'(',A6,'(',A6,'+1:',A6,'+', I6, ')), 0) ' )
C3211 FORMAT( '       CALL AD02_INITIALIZE(',A6,' - 1, X_value(:', I6, 
C    *               '),', /, '     *', 22X, 
C    *               A6,'(',A6,'(',A6,'+1:',A6,'+', I6, ')),', /, 
C    *        '     *                      DATA_AD02, 0)' )
 3211 FORMAT( '       CALL AD02_INITIALIZE_COMP(',A6,' - 1, X_value(:', 
     *               I6, '),', /, '     *', 22X, 
     *               A6,'(',A6,'(',A6,'+1:',A6,'+', I6, ')),', /, 
     *        '     *                      DATA_AD02, ERROR_AD02, 0)' )
 3220 FORMAT( '       ', A6, ' = X_value(', I6, ')' )
 3230 FORMAT( '       CALL AD01_INITIALIZE(0, X_value(:', I6, 
     *               '),', /, '     *', 22X, 
     *               A6,'(',A6,'(',A6,'+1:',A6,'+', I6, ')), 0) ' )
 3231 FORMAT( '       CALL AD02_INITIALIZE_COMP(0, X_value(:', I6, 
     *               '),', /, '     *', 22X, 
     *               A6,'(',A6,'(',A6,'+1:',A6,'+', I6, ')),', /, 
     *        '     *                      DATA_AD02, ERROR_AD02, 0)' )
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
 4040 FORMAT( '      GO TO (', 8( I5, :, ',' ), /,
     *      ( '     *       ', 8( I5, :, ',' ) ) )
 4050 FORMAT( '     *        ', 48X, '), ITYPE' )
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
C  END OF MAFNAD.
C
      END
C     THIS VERSION: 21ST JUNE 1990.

      SUBROUTINE OUTRN2( NELV, NINV, U, IOUTFF, IOUTFD, IOUTRA, 
     *                   ENAMES, INAMES, SINGLE, AD0 )
      INTEGER          NELV, NINV, IOUTFF, IOUTFD, IOUTRA
      LOGICAL          SINGLE
      CHARACTER * 4    AD0
      DOUBLE PRECISION U( NINV, NELV )
      CHARACTER * 10   ENAMES( * ), INAMES( * )
C
C  PRINT OUT THE GATHER AND SCATTER PART OF THE GENERATED RANGE ROUTINE
C  AND THE GATHER PART OF THE GENERATED FUNCTION EVALUATION ROUTINE.
C
      INTEGER          I, J, K
      DOUBLE PRECISION UIJ, ONE, ZERO, EPSMCH, DMACHR
      LOGICAL          ANYNNZ, OUTFF
      CHARACTER * 6    EVNAME, IVNAME
      INTRINSIC        DABS, MOD
      EXTERNAL         DMACHR
      DATA ZERO, ONE / 0.0D+0, 1.0D+0 /
      OUTFF = IOUTFF .GT. 0
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
                     IF ( OUTFF ) WRITE( IOUTFF, 3030 ) EVNAME
                     WRITE( IOUTFD, 3030 ) EVNAME
                     WRITE( IOUTRA, 4030 ) J
                  ELSE
                     IF ( OUTFF ) WRITE( IOUTFF, 3040 ) IVNAME, EVNAME
                     WRITE( IOUTFD, 3041 ) AD0, I, EVNAME
                     WRITE( IOUTRA, 4040 ) I, J
                  END IF
               ELSE
C
C  NONZERO HAS A VALUE OTHER THAN 1.
C
                  IF ( ANYNNZ ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3050 ) EVNAME, UIJ
                     WRITE( IOUTFD, 3050 ) EVNAME, UIJ
                     WRITE( IOUTRA, 4050 ) J, UIJ
                  ELSE
                     IF ( OUTFF ) 
     *                  WRITE( IOUTFF, 3060 ) IVNAME, EVNAME, UIJ
                     WRITE( IOUTFD, 3061 ) AD0, I, EVNAME, UIJ
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
                     IF ( OUTFF ) WRITE( IOUTFF, 3070 ) EVNAME
                     WRITE( IOUTFD, 3070 ) EVNAME
                     WRITE( IOUTRA, 4070 ) J
                  ELSE
                     IF ( OUTFF ) WRITE( IOUTFF, 3080 ) IVNAME, EVNAME
                     WRITE( IOUTFD, 3081 ) AD0, I, EVNAME
                     WRITE( IOUTRA, 4080 ) I, J
                  END IF
               ELSE
C
C  NONZERO HAS A VALUE OTHER THAN - 1.
C
                  IF ( ANYNNZ ) THEN
                     IF ( OUTFF ) WRITE( IOUTFF, 3090 ) EVNAME, - UIJ
                     WRITE( IOUTFD, 3090 ) EVNAME, - UIJ
                     WRITE( IOUTRA, 4090 ) J, - UIJ
                  ELSE
                     IF ( OUTFF ) 
     *                  WRITE( IOUTFF, 3100 ) IVNAME, EVNAME, - UIJ
                     WRITE( IOUTFD, 3101 ) AD0, I, EVNAME, - UIJ
                     WRITE( IOUTRA, 4100 ) I, J, - UIJ
                  END IF
               END IF
            END IF
            ANYNNZ = .TRUE.
            IF ( MOD( K, 19 ) .EQ. 0 ) THEN
               IF ( OUTFF ) WRITE( IOUTFF, 3040 ) IVNAME, IVNAME
               WRITE( IOUTFD, 3041 ) AD0, I, IVNAME
               WRITE( IOUTRA, 4112 ) I, I
            END IF
   30    CONTINUE
         IF ( .NOT. ANYNNZ ) THEN
            IF ( SINGLE ) THEN
              IF ( OUTFF ) WRITE( IOUTFF, 3120 ) IVNAME
               WRITE( IOUTFD, 3121 ) AD0, I
               WRITE( IOUTRA, 4111 ) I
            ELSE
              IF ( OUTFF ) WRITE( IOUTFF, 3110 ) IVNAME
               WRITE( IOUTFD, 3111 ) AD0, I
               WRITE( IOUTRA, 4110 ) I
            END IF
         END IF
         IF ( AD0 .EQ. 'AD01' ) THEN
            WRITE( IOUTFD, 3150 ) I, I
         ELSE
            WRITE( IOUTFD, 3151 ) I, I
         END IF
   40 CONTINUE
C
C  ----- THE GATHER HAS BEEN COMPLETED; WIND UP THE ELEMENT.
C
      WRITE( IOUTRA, 4020 )
      IF ( AD0 .EQ. 'AD01' ) THEN
         WRITE( IOUTFD, 3130 ) NINV, NINV
      ELSE
         WRITE( IOUTFD, 3131 ) NINV, NINV
      END IF
      DO 50 I   = 1, NINV
         IVNAME = INAMES( I )(1:6)
         WRITE( IOUTFD, 3140 ) IVNAME, I
   50 CONTINUE   
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 3030 FORMAT( '     *          + ', A6 )
 3040 FORMAT( '       ', A6, ' =   ', A6 )
 3041 FORMAT( '       X_', A4, '_int(', I6, ') =   ', A6 )
 3050 FORMAT( '     *          + ', A6, ' * ', F12.5 )
 3060 FORMAT( '       ', A6, ' =   ', A6, ' * ', F12.5 )
 3061 FORMAT( '       X_', A4, '_int(', I6, ') =   ', A6, ' * ', F12.5 )
 3070 FORMAT( '     *          - ', A6 )
 3080 FORMAT( '       ', A6, ' = - ', A6 )
 3081 FORMAT( '       X_', A4, '_int(', I6, ') = - ', A6 )
 3090 FORMAT( '     *          - ', A6, ' * ', F12.5 )
 3100 FORMAT( '       ', A6, ' = - ', A6, ' * ', F12.5 )
 3101 FORMAT( '       X_', A4, '_int(', I6, ') = - ', A6, ' * ', F12.5 )
 3110 FORMAT( '       ', A6, ' = 0.0D+0 ' )
 3111 FORMAT( '       X_', A4, '_int(', I6, ') = 0.0D+0 ' )
 3120 FORMAT( '       ', A6, ' = 0.0E+0 ' )
 3121 FORMAT( '       X_', A4, '_int(', I6, ') = 0.0E+0 ' )
 3130 FORMAT( '       CALL AD01_INITIALIZE(IFFLAG - 1, X_value(:', I6, 
     *               '),', /, '     *', 22X, 'X_int(:', I6, '), 0) ')
 3131 FORMAT( '       CALL AD02_INITIALIZE_COMP(IFFLAG - 1, X_value(:', 
     *                I6, '),', /, '     *', 22X, 'X_int(:', I6, '),', 
     *        ' DATA_AD02, ERROR_AD02, 0)' )
 3140 FORMAT( '       ', A6, ' = X_value(', I6, ')' )
 3150 FORMAT( '       CALL AD01_VALUE(X_AD01_int(', I6, 
     *                '), X_int(', I6, '))' )
 3151 FORMAT( '       CALL AD02_VALUE(X_AD02_int(', I6, 
     *                '), X_int(', I6, '), ERROR_AD02)' )
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
C  END OF OUTRN2.
C
      END
      CHARACTER * 6 FUNCTION NUNAME( I1, I2, I3, I4, I5, I6, IIRES, 
     *                               NINMAX, NRENAM, NINNAM, NLONAM, 
     *                               NMINAM, NEXNAM, NLMAX , NELTYP, 
     *                               YNAME , FIELDI, RENAME, INNAME, 
     *                               LONAME, MINAME, EXNAME, ETYPES )
C 
C  Find a name that does not occur in any other list
C
      INTEGER        I1, I2, I3, I4, I5, I6, IIRES, NLMAX, NINMAX
      INTEGER        NRENAM, NINNAM, NLONAM, NMINAM, NEXNAM, NELTYP
      CHARACTER * 6  YNAME
      CHARACTER * 8  FIELDI( IIRES )
      CHARACTER * 10 RENAME( NINMAX ), INNAME( NINMAX )
      CHARACTER * 10 LONAME( NINMAX )
      CHARACTER * 10 MINAME( NINMAX ), EXNAME( NINMAX )
      CHARACTER * 10 ETYPES( NLMAX  )
      CHARACTER * 1  CHARAC( 36 )
      DATA CHARAC / 'Z', 'Y', 'X', 'W', 'V', 'U', 'T', 'S', 'R', 
     *              'Q', 'P', 'O', 'N', 'M', 'L', 'K', 'J', 'I',
     *              'H', 'G', 'F', 'E', 'D', 'C', 'B', 'A', '0',
     *              '9', '8', '7', '6', '5', '4', '3', '2', '1' / 
   10 CONTINUE
C
C  Find the next name in the list
C
      I1 = I1 + 1
      IF ( I1 .EQ. 27 ) THEN
         I1 = 1
         I2 = I2 + 1
         IF ( I2 .EQ. 37 ) THEN
            I2 = 1
            I3 = I3 + 1
            IF ( I3 .EQ. 37 ) THEN
               I3 = 1
               I4 = I4 + 1
               IF ( I4 .EQ. 37 ) THEN
                  I4 = 1
                  I5 = I5 + 1
                  IF ( I5 .EQ. 37 ) THEN
                     I5 = 1
                     I6 = I6 + 1
                     IF ( I6 .EQ. 37 ) THEN
                        write( 6, * ) ' no characters left '
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END IF
      NUNAME = CHARAC( I1 ) // CHARAC( I2 ) // CHARAC( I3 ) // 
     *         CHARAC( I4 ) // CHARAC( I5 ) // CHARAC( I6 )
C
C  See if the name has already been used
C
      DO 110 I = 1, NRENAM
         IF ( RENAME( I )(1:6) .EQ. NUNAME ) GO TO 10
  110 CONTINUE
      DO 120 I = 1, NINNAM
         IF ( INNAME( I )(1:6) .EQ. NUNAME ) GO TO 10
  120 CONTINUE
      DO 130 I = 1, NLONAM
         IF ( LONAME( I )(1:6) .EQ. NUNAME ) GO TO 10
  130 CONTINUE
      DO 140 I = 1, NMINAM
         IF ( MINAME( I )(1:6) .EQ. NUNAME ) GO TO 10
  140 CONTINUE
      DO 150 I = 1, NEXNAM
         IF ( EXNAME( I )(1:6) .EQ. NUNAME ) GO TO 10
  150 CONTINUE
      DO 160 I = 1, NELTYP
         IF ( ETYPES( I )(1:6) .EQ. NUNAME ) GO TO 10
  160 CONTINUE
      DO 170 I = 1, IIRES
         IF ( FIELDI( I )(1:6) .EQ. NUNAME ) GO TO 10
  170 CONTINUE
      IF ( NUNAME .EQ. YNAME ) GO TO 10
      RETURN
C
C  End of NUNAME
C
      END
