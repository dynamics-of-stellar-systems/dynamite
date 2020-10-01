C  THIS VERSION: 30/06/2004 AT 16:45:00 GMT
C ** Correction report.
C ** Correction -5. 07/09/00: Check for non-useful transformations added
C ** Correction -4. 07/09/00: Error return for too small NIMAX added
C ** Correction -3. 21/02/00: Code to handle incomplete/missing data added
C ** Correction -2. 21/02/00: QSECTION added as alias for QUADOBJ
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C ** Correction 0. 20/12/99: Array holding variable types introduced.
C ** Correction 1. 30/11/93: 3 lines interchanged **
C ** Correction 2. 26/02/01: 4 dummy arguments removed from SVAR1 **
C ** Correction 3. 26/02/01: 1 dummy argument removed from SBOUND **
C ** Correction 4. 26/02/01: 1 dummy argument removed from SQHESS **
C ** Correction 5. 26/02/01: 2 dummy arguments removed from SETYPE **
C ** Correction 6. 26/02/01: 2 dummy arguments removed from SGTYPE **
C ** Correction 7. 26/02/01: 2 dummy arguments removed from SOBBND **
C ** Correction 8. 26/02/01: 1 dummy argument removed from PROCAA **
C ** Correction 9. 27/02/01: Character debug output format increased **
C ** Correction 10. 16/04/02: No double assignment of internals for quad terms
C ** Correction 11  30/06/04: Assigned goto statements replaced
C ** End of Correction report.
      SUBROUTINE GPSMPS( LA    , NMAX  , NGMAX , NOMAX , NLMAX , NELMAX,
     *                   NIMAX , NETMAX, NEVMAX, NGRMAX, NSMAX , NEPMAX,
     *                   NGPMAX, NBMAX , NOBMAX, NNZA  , LENGTH, N , NG,
     *                   NOBJ  , NCONST, NRANGE, NBND  , NSTART, NELTYP,
     *                   NGRTYP, NLVARS, NNLVRS, NLISGP, LIWK  ,
     *                   NELNUM, NELING, NARRAY, NINDEX, NEGMAX, NEPVMX,
     *                   NGPVMX, MAXINS, MAXLEV, MAXARA, NRLNDX, NOBBND,
     *                   PNAME , ICOORD, IELING, INLIST, ITABLE, ISTATE,
     *                   IELV  , IINV  , ITYPEE, IDROWS, IELVAR, ISTADG,
     *                   ITYPEG, IEPA  , IGPA  , IWK   , ISTEP , ISTAEV,
C ** Correction 0. 20/12/99: Array holding variable types introduced.
     *                   ISTGP , INDVAL, INSTR , NINSTR, IARRAY, ITYPEV,
     *                   A, BND, VSTART, CSTART, RSCALE, CSCALE, RDROWS, 
     *                   REALVL, DFAULT, RVALUE, VARRAY, EPVALU, BNDFLT,
     *                   GPVALU, GPTEMP, FARRAY, FBOUND, WEIGHT, NAMIIN,
     *                   NAMRIN, GNAMES, VNAMES, BNAMES, ETYPES, INAMES,
     *                   LNAMES, ONAMES, ENAMES, SNAMES, ANAMES, GTYPES,
     *                   EPNAME, GPNAME, OBNAME, ARRAY , CARRAY, KEY   ,
     *                   SINGLE, INPUT , IOUT  , INFORM, DEBUG )
      INTEGER          LA    , NMAX  , NGMAX , NOMAX , NLMAX , NELMAX
      INTEGER          NEVMAX, NGRMAX, NSMAX , NOBMAX, NOBBND, NRLNDX
      INTEGER          NEPMAX, NGPMAX, NBMAX , NNZA  , LENGTH, N , NG
      INTEGER          NCONST, NRANGE, NBND  , NSTART, NELTYP, NGRTYP
      INTEGER          NLVARS, NNLVRS, NELNUM, NELING, NARRAY, NINDEX
      INTEGER          NEGMAX, NEPVMX, NGPVMX, MAXINS, MAXLEV, MAXARA
      INTEGER          NLISGP, INPUT , NETMAX, IOUT  , INFORM, NIMAX
      INTEGER          NOBJ  , LIWK
      LOGICAL          SINGLE, DEBUG
      CHARACTER * 8    PNAME
      INTEGER          INLIST( LENGTH ), ISTAEV( NELMAX )
      INTEGER          ISTATE( NGMAX ), ITABLE ( LENGTH )
      INTEGER          IELV  ( NLMAX  ), IINV  ( NLMAX )
      INTEGER          ITYPEE( NELMAX ), IELING( NEGMAX, 2 )
      INTEGER          IEPA  ( NLMAX  ), IGPA  ( NGRMAX )
      INTEGER          IDROWS( 2, NGMAX ), ITYPEG( NGMAX )
C ** Correction 0. 20/12/99: Array holding variable types introduced.
      INTEGER          IELVAR( NEVMAX ), ITYPEV( NMAX )
      INTEGER          ISTEP ( NELMAX ), ISTGP( NGMAX ), IWK( LIWK )
      INTEGER          INDVAL( NINDEX ), INSTR( 5, MAXINS, MAXLEV )
      INTEGER          NINSTR( MAXLEV ), IARRAY( 5, 3, MAXARA )
      INTEGER          ICOORD( LA, 2 ) , ISTADG( NGMAX )
      DOUBLE PRECISION GPTEMP( NGPVMX )
      DOUBLE PRECISION EPVALU( NEPVMX ), GPVALU( NGPVMX ), DFAULT( NMAX)
      DOUBLE PRECISION A( LA ), BND( 2, NMAX, NBMAX ), REALVL( NRLNDX )
      DOUBLE PRECISION BNDFLT( 2, NBMAX ), CSTART( NGMAX, NSMAX )
      DOUBLE PRECISION RSCALE( NGMAX ), CSCALE( NMAX )
      DOUBLE PRECISION RDROWS( 2, NGMAX ), VSTART( NMAX, NSMAX )
      DOUBLE PRECISION RVALUE( MAXARA, 3 ), VARRAY( 2, MAXARA )
      DOUBLE PRECISION FBOUND( 2, NOBMAX ), WEIGHT( NEGMAX )
      CHARACTER * 2    FARRAY( MAXARA )
      CHARACTER * 10   NAMIIN( NINDEX ), NAMRIN( NRLNDX )
      CHARACTER * 10   BNAMES( NBMAX  )
      CHARACTER * 10   GNAMES( NGMAX  ), VNAMES( NMAX   )
      CHARACTER * 10   ETYPES( NLMAX  ), INAMES( NIMAX  )
      CHARACTER * 10   LNAMES( NELMAX ), OBNAME( NOBMAX )
      CHARACTER * 10   EPNAME( NEPMAX ), GPNAME( NGPMAX )
      CHARACTER * 10   ONAMES( NOMAX  ), ENAMES( NETMAX )
      CHARACTER * 10   SNAMES( NSMAX )
      CHARACTER * 10   ANAMES( NGRMAX ), GTYPES( NGRMAX )
      CHARACTER * 10   ARRAY( 3, MAXARA ), CARRAY( 2, MAXARA )
      CHARACTER * 12   KEY   ( LENGTH )
C
C  .....................................................................
C
C  READ A GPS MPS DATA FILE.
C  -------------------------
C
C  NICK GOULD 12/08/1989
C  FOR CGT PRODUCTIONS.
C
C  MPS INDICATOR CARDS.
C  --------------------
C
C  DEFINITION   PURPOSE.
C  ----------   --------
C  NAME         PROBLEM NAME.
C  ROWS         NAMES OF ROWS (ALIAS GROUP NAMES).
C  COLUMNS      NAMES OF COLUMNS (ALIAS VARIABLE NAMES).
C  RHS          RIGHT-HAND-SIDES (ALIAS CONSTANT TERMS IN GROUPS).
C  RHS'         ALIAS FOR RHS.
C  RANGES       ADDITIONAL BOUNDS ON ROWS.
C  BOUNDS       BOUNDS ON COLUMNS.
C  ENDATA       END OF INPUT DATA.
C
C  ADDITIONAL INDICATOR CARDS.
C  ---------------------------
C
C  DEFINITION   PURPOSE.
C  ----------   --------
C  GROUPS       ALIAS FOR ROWS.
C  CONSTRAINTS  ALIAS FOR ROWS.
C  VARIABLES    ALIAS FOR COLUMNS.
C  CONSTANTS    ALIAS FOR RHS.
C  START_POINT  ESTIMATE OF MINIMIZER.
C  HESSIAN      QUADRATIC TERMS
C  QUADRATIC    ALIAS FOR HESSIAN
C  QUADS        ALIAS FOR HESSIAN
C  QUADOBJ      ALIAS FOR HESSIAN
C  QSECTION     ALIAS FOR HESSIAN
C  ELEMENT_TYPE TYPES OF NONLINEAR ELEMENTS.
C  ELEMENT_USES DEFINITIONS OF NONLINEAR ELEMENTS.
C  GROUP_TYPE   TYPES OF NONTRIVIAL GROUPS.
C  GROUP_USES   DEFINITIONS OF GROUPS.
C
C  DATA CARD DESCRIPTION.
C  ----------------------
C
C  SEE 'A PROPOSAL FOR A STANDARD DATA INPUT FORMAT FOR LARGE-SCALE
C       NONLINEAR PROGRAMMING PROBLEMS', SECTION 2,
C       A. R. CONN, N. I. M. GOULD AND PH. L. TOINT,
C       REPORT CS-89-61, DEPT OF COMPUTER SCIENCE, U. OF WATERLOO,
C       WATERLOO, ONTARIO, N2L3G1, CANADA.
C
C  -------------------------------------------------------------------
C  RETURNS WITH NEGATIVE VALUES OF INFORM INDICATE THAT INSUFFICIENT
C  ARRAY SPACE HAS BEEN ALLOWED, AS FOLLOWS:
C
C    INFORM = - 1  WHEN LENGTH NOT LARGE ENOUGH
C    INFORM = - 2  WHEN NNZA .GT. LA
C    INFORM = - 3  WHEN NELTYP .GE. NLMAX
C    INFORM = - 4  WHEN NGRTYP .GE. NGRMAX
C    INFORM = - 5  WHEN NOBJ .GT. NOMAX
C    INFORM = - 6  WHEN NG .GT. NGMAX
C    INFORM = - 7  WHEN N .GT. NMAX
C    INFORM = - 8  WHEN NSTART .GT. NSMAX
C    INFORM = - 9  WHEN NELNUM .GT. NELMAX
C    INFORM = - 10 WHEN NELING .GT. NEGMAX
C    INFORM = - 11 WHEN NINSTR( 1 OR 2 OR 3 ) .GT. MAXINS
C    INFORM = - 12 WHEN NARRAY .GT. MAXARA
C    INFORM = - 13 WHEN NBND .GT. NBMAX
C    INFORM = - 14 WHEN NELN .GT. NETMAX
C    INFORM = - 15 WHEN NLISEV .GT. NEVMAX
C    INFORM = - 16 WHEN NINN .GT. NIMAX
C    INFORM = - 17 WHEN NLISEP .GT. NEPVMX
C    INFORM = - 18 WHEN NLISGP .GT. NGPVMX
C    INFORM = - 19 WHEN NEPN .GT. NEPMAX
C    INFORM = - 20 WHEN NGPN .GT. NGPMAX
C    INFORM = - 21 WHEN NUSEIN .GT. NINDEX
C    INFORM = - 22 WHEN NUSERE .GT. NRLNDX
C    INFORM = - 23 WHEN NOBBND .GT. NOBMAX
C
C  .....................................................................
C
      INTEGER          I, IP, IS, INTYPE, INTYPO, J, K, K1, K2, L
      INTEGER          IFREE, IFIELD, NOVALS, NVAR, NCOL, NELN, NINN
      INTEGER          MBLANK, MNAME, MROWS, MGROUP, MCNSTR, MOBBND
      INTEGER          MVARS, MCONST, MRHS, MRHSP, MRANGE, MBOUND, MCOLS
      INTEGER          MSTART, METYPE, MGTYPE, MEUSES, MGUSES, MENDAT
      INTEGER          MFREE, MFIXED, MAXNUL, NLINES, ILINES, MXRECL
      INTEGER          NEPN, NGPN, NGRUPE, NELMNT, NUSEIN, NUSERE
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      INTEGER          MQHESS, MQUADO, MQUADR, MQUADS, IPTYPE, ISTYPE
C ** Correction -2. 21/02/00: QSECTION added as alias for QUADOBJ
      INTEGER          MQSECT
      INTEGER          NLISEP, NLISEV, L1, L2, L3, LEVEL2, NDTYPE
      INTEGER          LEVEL3, LEV1, LEV2, LEV3, LEV1S, LEV2S, LEV3S
      INTEGER          LEV1E, LEV2E, LEV3E, LEV1I, LEV2I, LEV3I, LEVL3A
      INTEGER          LEVEL, IJUMP, NINCRS, LINENO, LOOP( 4 )
      DOUBLE PRECISION BIG, ONE, ZERO, VALUE4, VALUE6
      LOGICAL          DEFNAM, INCARD, INREP, DEFAUT, DOLOOP, STRTGU
      LOGICAL          ENDBND, ENDST, SETANA, ENDELT, ENDGRT, DELSET
      LOGICAL          ENDELU, GRP1ST, GRPYET, VARYET, FIXED, DGRSET
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      LOGICAL          QGROUP, QSQR, QPROD
      CHARACTER * 2    FIELD1, COLFIE
      CHARACTER * 10   FIELD2, FIELD3, FIELD5, GRUPE, ELMNT
      CHARACTER * 10   DETYPE, DGTYPE
      CHARACTER * 12   FIELD
      EXTERNAL         SGRP1, SGRP2, SVAR1, SVAR2, SBOUND, SSTART
      EXTERNAL         SETYPE, SEUSES, SGTYPE, SGUSES, SOBBND, GETRIN
      EXTERNAL         GETVAL, GETLIN, GETIIN, PROCAI, PROCAD, PROCAA
      EXTERNAL         HASHA , HASHB , HASHC, REORDA
      INTRINSIC        ABS
C
C  PARAMETER DEFINITIONS.
C
      PARAMETER        ( MXRECL = 160 )
      CHARACTER * 160  NULINE, BLNKLN
      PARAMETER        ( NINCRS = 23 )
      CHARACTER * 6    INCRSE( NINCRS )
      PARAMETER        ( MBLANK =  1, MFIXED =  2, MFREE  = 3  )
      PARAMETER        ( MNAME  =  4,  MROWS  =  5 )
      PARAMETER        ( MGROUP =  6, MCNSTR =  7, MCOLS  =  8 )
      PARAMETER        ( MVARS  =  9, MCONST = 10, MRHS   = 11 )
      PARAMETER        ( MRHSP  = 12, MRANGE = 13, MBOUND = 14 )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      PARAMETER        ( MSTART = 15, MQHESS = 16, MQUADR = 17 )
C ** Correction -2. 21/02/00: QSECTION added as alias for QUADOBJ
      PARAMETER        ( MQUADS = 18, MQUADO = 19, MQSECT = 20 )
      PARAMETER        ( METYPE = 21, MEUSES = 22, MGTYPE = 23 )
      PARAMETER        ( MGUSES = 24, MOBBND = 25, MENDAT = 26 )
C ** Correction  10. 16/04/02: 3 lines added
      CHARACTER * 10   CQSQR, CQPROD
      PARAMETER      ( CQSQR  = '123456789S' )
      PARAMETER      ( CQPROD = '123456789P' )
      INTEGER          LENIND( MENDAT )
      CHARACTER * 12   INDIC8( MENDAT ), HEADER
      PARAMETER        ( MAXNUL = 20 )
      CHARACTER * 65   NULINA( MAXNUL )
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0, BIG = 1.0D+20 )
C
C  DATA DECLARATIONS.
C
      DATA INCRSE / 'LENGTH', 'LA    ', 'NLMAX ', 'NGRMAX', 'NOMAX ',
     *              'NGMAX ', 'NMAX  ', 'NSMAX ', 'NELMAX', 'NEGMAX',
     *              'MAXINS', 'MAXARA', 'NBMAX ', 'NETMAX', 'NEVMAX',
     *              'NIMAX ', 'NEPVMX', 'NGPVMX', 'NEPMAX', 'NGPMAX',
     *              'NINDEX', 'NRLNDX', 'NOBMAX' /
      DATA INDIC8( MBLANK ) / '            ' /, LENIND( MBLANK ) / 0  /
      DATA INDIC8( MFIXED ) / 'FIXED FORMAT' /, LENIND( MFIXED ) / 12 /
      DATA INDIC8( MFREE  ) / 'FREE FORMAT ' /, LENIND( MFREE  ) / 11 /
      DATA INDIC8( MNAME  ) / 'NAME        ' /, LENIND( MNAME  ) / 4  /
      DATA INDIC8( MROWS  ) / 'ROWS        ' /, LENIND( MROWS  ) / 4  /
      DATA INDIC8( MGROUP ) / 'GROUPS      ' /, LENIND( MGROUP ) / 6  /
      DATA INDIC8( MCNSTR ) / 'CONSTRAINTS ' /, LENIND( MCNSTR ) / 11 /
      DATA INDIC8( MCOLS  ) / 'COLUMNS     ' /, LENIND( MCOLS  ) / 7  /
      DATA INDIC8( MVARS  ) / 'VARIABLES   ' /, LENIND( MVARS  ) / 9  /
      DATA INDIC8( MCONST ) / 'CONSTANTS   ' /, LENIND( MCONST ) / 9  /
      DATA INDIC8( MRHS   ) / 'RHS         ' /, LENIND( MRHS   ) / 3  /
      DATA INDIC8( MRHSP  ) / 'RHS''       ' /, LENIND( MRHSP  ) / 4  /
      DATA INDIC8( MRANGE ) / 'RANGES      ' /, LENIND( MRANGE ) / 6  /
      DATA INDIC8( MBOUND ) / 'BOUNDS      ' /, LENIND( MBOUND ) / 6  /
      DATA INDIC8( MSTART ) / 'START POINT ' /, LENIND( MSTART ) / 11 /
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      DATA INDIC8( MQHESS ) / 'HESSIAN     ' /, LENIND( MQHESS ) / 7  /
      DATA INDIC8( MQUADR ) / 'QUADRATIC   ' /, LENIND( MQUADR ) / 9  /
      DATA INDIC8( MQUADS ) / 'QUADS       ' /, LENIND( MQUADS ) / 5  /
      DATA INDIC8( MQUADO ) / 'QUADOBJ     ' /, LENIND( MQUADO ) / 7  /
C ** Correction -2. 21/02/00: QSECTION added as alias for QUADOBJ
      DATA INDIC8( MQSECT ) / 'QSECTION    ' /, LENIND( MQSECT ) / 8  /
      DATA INDIC8( METYPE ) / 'ELEMENT TYPE' /, LENIND( METYPE ) / 12 /
      DATA INDIC8( MEUSES ) / 'ELEMENT USES' /, LENIND( MEUSES ) / 12 /
      DATA INDIC8( MGTYPE ) / 'GROUP TYPE  ' /, LENIND( MGTYPE ) / 10 /
      DATA INDIC8( MGUSES ) / 'GROUP USES  ' /, LENIND( MGUSES ) / 10 /
      DATA INDIC8( MOBBND ) / 'OBJECT BOUND' /, LENIND( MOBBND ) / 12 /
      DATA INDIC8( MENDAT ) / 'ENDATA      ' /, LENIND( MENDAT ) / 6  /
C
C  SET INITIAL VALUES FOR INTEGER VARIABLES.
C
      INTYPE = 1
      INTYPO = 1
      INFORM = 0
      LINENO = 0
      NVAR   = 0
      NNZA   = 0
      NG     = 0
      NBND   = 0
      NSTART = 0
      NOBJ   = 0
      NELTYP = 0
      NGRTYP = 0
      NELNUM = 0
      NGRUPE = 0
      NLISEV = 0
      NLISEP = 0
      NLISGP = 0
      NELING = 0
      NUSEIN = 0
      NUSERE = 0
      NOBBND = 0
      NDTYPE = 0
      NELN   = 0
      NINN   = 0
      NEPN   = 0
      NLVARS = - 1
      NNLVRS = - 1
      NCONST = - 1
      NRANGE = - 1
      LEVEL  = 0
      ILINES = 0
      NLINES = 0
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      IPTYPE = 0
      ISTYPE = 0
C
C  SET INITIAL VALUES FOR LOGICAL VARIABLES.
C
      DEFNAM = .FALSE.
      DOLOOP = .FALSE.
      ENDBND = .FALSE.
      ENDST  = .FALSE.
      ENDELT = .FALSE.
      ENDELU = .FALSE.
      ENDGRT = .FALSE.
      STRTGU = .FALSE.
      GRPYET = .FALSE.
      VARYET = .FALSE.
      DELSET = .FALSE.
      DGRSET = .FALSE.
      GRP1ST = .TRUE.
      FIXED  = .TRUE.
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      QGROUP = .FALSE.
      QSQR   = .FALSE.
      QPROD  = .FALSE.
C
C  SET INITIAL VALUES FOR REAL VARIABLES.
C
      VALUE4 = 0.0D+0
      VALUE6 = 0.0D+0
C
C  SET UP ITABLE DATA.
C
      CALL HASHA ( LENGTH, ITABLE )
C
C  INITIALIZE ROW DATA.
C
      DO 10 I           = 1, NGMAX
         RSCALE( I )    = ONE
   10 CONTINUE
C
C  INITIALIZE COLUMN DATA.
C
      DO 20 I           = 1, NMAX
C ** Correction 0. 20/12/99: Array holding variable types introduced.
         ITYPEV( I )    = 0
         CSCALE( I )    = ONE
         BND( 1, I, 1 ) = ZERO
         BND( 2, I, 1 ) = BIG
   20 CONTINUE
C
C  INITIALIZE DICTIONARY DATA.
C
      DO 30 I        = 1, LENGTH
         INLIST( I ) = 0
   30 CONTINUE
C
C  SET A BLANK LINE.
C
      DO 40 I = 1, MXRECL
         BLNKLN( I: I ) = ' '
   40 CONTINUE
C
C  START OF MAIN LOOP.
C
  100 CONTINUE
      IF ( ILINES + 1 .GT. NLINES ) THEN
C
C  READ NEXT LINE FROM THE INPUT FILE.
C
         LINENO = LINENO + 1
         NULINE = BLNKLN
         IF ( FIXED ) THEN
C ** Correction -3. 21/02/00: Code to handle incomplete/missing data added
            READ ( INPUT, 1000, END = 810, ERR = 810 ) NULINE
            IF ( IOUT .GT. 0 .AND. DEBUG ) WRITE( IOUT, 2990 )
     *           LINENO, NULINE
         ELSE
C ** Correction -3. 21/02/00: Code to handle incomplete/missing data added
            READ ( INPUT, 1010, END = 810, ERR = 810 ) NULINE
            IF ( IOUT .GT. 0 .AND. DEBUG ) WRITE( IOUT, 2970 )
     *           LINENO, NULINE
C
C  IF THE CARD IS IN FREE FORMAT, TRANSLATE IT INTO FIXED FORMAT.
C
            CALL  FREEFM( NULINE, MXRECL, MENDAT, INDIC8, LENIND,
     *                    NULINA, MAXNUL, NLINES, .TRUE., INFORM, IOUT )
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
      IF ( HEADER .EQ. INDIC8( MBLANK ) ) THEN
         IF (  NULINE( 13: 14 ) .EQ. '  ' .AND.
     *         NULINE( 15: 24 ) .EQ. '          ' .AND.
     *         NULINE( 40: 49 ) .EQ. '          ' ) GO TO 100
      END IF
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
C  CHECK THAT THE FIRST ENCOUNTERED INDICATOR CARD IS THE NAME CARD.
C
         IF ( .NOT. DEFNAM  ) THEN
            IF ( HEADER .NE. INDIC8( MNAME ) ) THEN
               INFORM = 1
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2010 )
               GO TO 800
            ELSE
C
C  INDICATOR CARD IS NAME.
C  -----------------------
C
               DEFNAM = .TRUE.
               PNAME  = NULINE( 15: 22 )
               GO TO 100
            END IF
         END IF
         INCARD = .TRUE.
C
C  AN INDICATOR CARD HAS BEEN FOUND.
C
         IF ( .NOT. GRP1ST ) INTYPE = MROWS
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
         IF ( INTYPE .EQ. MGROUP .OR. INTYPE .EQ. MCNSTR) INTYPE = MROWS
         IF ( INTYPE .EQ. MRHS .OR. INTYPE .EQ. MRHSP ) INTYPE = MCONST
         IF ( INTYPE .EQ. MVARS ) INTYPE = MCOLS
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C ** Correction -2. 21/02/00: QSECTION added as alias for QUADOBJ
         IF ( INTYPE .EQ. MQUADR .OR. INTYPE .EQ. MQUADS .OR. 
     *        INTYPE .EQ. MQUADO .OR. INTYPE .EQ. MQSECT ) 
     *      INTYPE = MQHESS
C
C  ENSURE THAT THE GROUPS AND VARIABLES SECTIONS DO NOT GET MIXED UP.
C
         IF ( .NOT. GRP1ST .AND. VARYET .AND. INTYPE .EQ. MCOLS ) THEN
            INFORM = 21
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2210 )
            GO TO 800
         END IF
         IF ( INTYPE .EQ. MROWS ) GRPYET = .TRUE.
         IF ( INTYPE .EQ. MCOLS ) VARYET = .TRUE.
         IF ( VARYET .AND. .NOT. GRPYET ) GRP1ST = .FALSE.
C
C  ENSURE THAT PREVIOUSLY STARTED DO-LOOPS HAVE BEEN FINISHED.
C
         IF ( INTYPE .NE. INTYPO .AND. DOLOOP ) THEN
            INFORM = 38
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2380 )
            GO TO 800
         END IF
         INTYPO = INTYPE
C
C  ALL OF THE LINEAR VARIABLES HAVE BEEN SPECIFIED.
C
         IF ( INTYPE .GE. MCONST ) THEN
            IF ( NLVARS .LT. 0 ) THEN
               NLVARS = NVAR
               N      = NLVARS
            END IF
         END IF
C
C  THE RIGHT-HAND-SIDE VECTORS HAVE BEEN COMPLETED.
C
         IF ( INTYPE .GE. MRANGE ) THEN
            IF ( NCONST .LT. 0 ) NCONST = NVAR - NLVARS
            IF ( NCONST .EQ. 0 ) THEN
               NCONST = 1
               NVAR   = NVAR + 1
            END IF
         END IF
C
C  THE RANGE VECTORS HAVE BEEN COMPLETED.
C
         IF ( INTYPE .GE. MBOUND ) THEN
            IF ( NRANGE .LT. 0 ) NRANGE = NVAR - NCONST - NLVARS
         END IF
C
C  THE BOUND VECTORS HAVE BEEN COMPLETED.
C
         IF ( INTYPE .GE. MSTART ) THEN
            IF ( .NOT. ENDBND ) THEN
               ENDBND = .TRUE.
               IF ( NBND .EQ. 0 ) NBND = 1
            END IF
         END IF
C
C  THE STARTING VECTORS HAVE BEEN COMPLETED.
C
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
         IF ( INTYPE .GE. MQHESS ) THEN
            IF ( .NOT. ENDST ) THEN
               ENDST = .TRUE.
               IF ( NSTART .EQ. 0 ) NSTART = 1
            END IF
         END IF
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C
C  THE QUADRATIC HESSIAN HAS BEEN COMPLETED.
C
         IF ( INTYPE .GE. METYPE ) THEN
         END IF
C
C  THE ELEMENT TYPES HAVE ALL BEEN SPECIFIED.
C
         IF ( INTYPE .GE. MEUSES ) THEN
            IF ( .NOT. ENDELT ) THEN
               ENDELT = .TRUE.
C
C  IF THE LAST ELEMENT HAS NO EXPLICIT INTERNAL REPRESENTATION,
C  USE ITS ELEMENTAL REPRESENTATION.
C
               IF ( NELTYP .GT. 0 ) THEN
C ** Correction  10. 16/04/02: 2 lines added
                  IF ( ETYPES( NELTYP ) .NE. CQSQR .AND.
     *                 ETYPES( NELTYP ) .NE. CQPROD ) THEN
                  IF ( .NOT. INREP ) THEN
                     DO 150 K          = IELV( NELTYP ), NELN
                        NINN           = NINN + 1
C ** Correction -4. 07/09/00: Error return for too small NIMAX added
                        IF ( NINN .GT. NIMAX ) THEN
                           INFORM = - 16
                           GO TO 700
                        END IF
                        INAMES( NINN ) = ENAMES( K )
  150                CONTINUE
C ** Correction -5a. 07/09/00: Check for non-useful transformations added
                  ELSE
                     IF ( NINN - IINV( NELTYP ) .GE.
     *                    NELN - IELV( NELTYP ) ) THEN
                        INFORM = 76
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2760 )
                        GO TO 800
                     END IF
                  END IF
C ** Correction  10. 16/04/02: 1 line added
                  END IF
                  IF ( NELTYP .GE. NLMAX ) THEN
                     INFORM = - 3
                     GO TO 700
                  END IF
               END IF
               IELV( NELTYP + 1 ) = NELN + 1
               IINV( NELTYP + 1 ) = NINN + 1
               IEPA( NELTYP + 1 ) = NEPN + 1
            END IF
         END IF
C
C  THE NONLINEAR ELEMENTS HAVE ALL BEEN SPECIFIED.
C
         IF ( INTYPE .GE. MGTYPE ) THEN
C
C  CHECK IF THERE ARE ANY NONLINEAR VARIABLES.
C
            IF ( NNLVRS .LT. 0 ) NNLVRS = N - NLVARS
            IF ( .NOT. ENDELU ) THEN
               ENDELU = .TRUE.
               IF ( NELNUM .GT. 0 ) THEN
C
C  CHECK THAT THE NONLINEAR ELEMENTS HAVE BEEN COMPLETELY SPECIFIED.
C  FIRST CHECK THE PARAMETER VALUES HAVE BEEN SET.
C
                  DO 153 J = 1, NELNUM
                     K     = ITYPEE( J )
                     IP    = IEPA( K ) - 1
                     K1    = IEPA( K + 1 ) - IEPA( K )
                     K2    = ISTEP( J ) - 1
                     DO 151 I = 1, K1
                        IF ( EPVALU( K2 + I ) .EQ. BIG ) THEN
                           INFORM = 28
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2280 )
     *                        LNAMES( J ), EPNAME( IP + I )
                        END IF
  151                CONTINUE
C
C  NOW CHECK THE ELEMENTAL VARIABLES HAVE BEEN SET.
C
                     IS       = IELV( K ) - 1
                     K1       = IELV( K + 1 ) - IELV( K )
                     K2       = ISTAEV( J ) - 1
                     DO 152 I = 1, K1
                        IF ( IELVAR( K2 + I ) .EQ. 0 ) THEN
                           INFORM = 16
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2160 )
     *                        LNAMES( J ), ENAMES( IS + I )
                        END IF
  152                CONTINUE
  153             CONTINUE
               END IF
               IF ( INFORM .NE. 0 ) RETURN
            END IF
            ISTEP( NELNUM + 1 ) = NLISEP + 1
            ISTAEV( NELNUM + 1 ) = NLISEV + 1
         END IF
C
C  THE GROUP TYPES HAVE ALL BEEN SPECIFIED.
C
         IF ( INTYPE .GE. MGUSES .AND. NGRTYP .GT. 0 ) THEN
C
C  CHECK THAT THE ARGUMENT FOR THE LAST GROUP-TYPE HAS BEEN SET.
C
            IF ( .NOT. ENDGRT ) THEN
               ENDGRT = .TRUE.
               IF ( .NOT. SETANA ) THEN
                  INFORM = 25
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2250 )
                  RETURN
               END IF
               IF ( NGRTYP .GE. NGRMAX ) THEN
                  INFORM = - 4
                  GO TO 700
               END IF
               IGPA( NGRTYP + 1 ) = NGPN + 1
            END IF
         END IF
         IF ( INTYPE .EQ. MENDAT ) THEN
C
C  CHECK THAT THE GROUPS HAVE BEEN COMPLETELY SPECIFIED BY
C  CHECKING THAT THE PARAMETER VALUES HAVE BEEN SET.
C
            NLISGP = 0
            DO 157 J = 1, NG
               K     = ITYPEG( J )
               IF ( K .LT. 0 ) THEN
                  K = - K - 1
                  ITYPEG( J ) = K
                  ISTGP( J ) = NLISGP + 1
                  IF ( K .NE. 0 ) THEN
                    K1 = IGPA( K + 1 ) - IGPA( K )
                    IF ( K1 .GT. 0 ) THEN
                       IP = IGPA( K ) - 1
                       K2 = ISTGP( J ) - 1
                       INFORM = 34
                       DO 155 I = 1, K1
                          NLISGP = NLISGP + 1
                          GPVALU( NLISGP ) = BIG
                          IF ( IOUT .GT. 0 ) WRITE( IOUT, 2340 )
     *                       GNAMES( J ), GPNAME( IP + I )
  155                 CONTINUE
                    END IF 
                  END IF 
               ELSE IF ( K .EQ. 0 ) THEN
                  ISTGP( J ) = NLISGP + 1
               ELSE
                  ISTGP( J ) = NLISGP + 1
                  IP    = IGPA( K ) - 1
                  K1    = IGPA( K + 1 ) - IGPA( K )
                  K2    = ISTGP( J ) - 1
                  DO 156 I = 1, K1
                     NLISGP = NLISGP + 1
                     GPVALU( NLISGP ) = GPTEMP( K2 + I )
                     IF ( GPVALU( NLISGP ) .EQ. BIG ) THEN
                        INFORM = 34
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2340 )
     *                      GNAMES( J ), GPNAME( IP + I )
                     END IF
  156             CONTINUE
               END IF
  157       CONTINUE
            ISTGP( NG + 1 ) = NLISGP + 1
            IF ( INFORM .NE. 0 ) RETURN
C
C  SORT THE LIST OF ELEMENTS FOR EACH GROUP, SO THAT THE
C  ELEMENTS FOR GROUP I PRECEDE THOSE FOR GROUP I + 1, I = 1, NG - 1.
C
            IF ( NELING .GT. 0 ) THEN
               CALL REORDA( NG, NELING, IELING( 1, 1 ), IELING( 1, 2 ),
     *                      WEIGHT, ISTADG, IWK )
            ELSE
               DO 158 I       = 1, NG + 1
                  ISTADG( I ) = 1
  158          CONTINUE
            END IF
         END IF
C
C  INDICATOR CARD IS ENDATA.
C  -------------------------
C
         IF ( INTYPE .EQ. MENDAT ) GO TO 900
         GO TO 100
      END IF
C
C  A DATA CARD HAS BEEN FOUND.
C
      INCARD = .FALSE.
C
C  READ THE CHARACTER FIELDS 1, 2, 3 AND 5 FROM THE NEW DATA LINE.
C
      FIELD1 = NULINE(  2:  3 )
      FIELD2 = NULINE(  5: 14 )
      FIELD3 = NULINE( 15: 24 )
      FIELD5 = NULINE( 40: 49 )
C
C  START OF A DO-LOOP.
C  ===================
C
      IF ( FIELD1 .EQ. 'DO' ) THEN
         IF ( LEVEL .GE. 3 ) THEN
            INFORM = 13
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2130 )
            GO TO 800
         END IF
C
C  THIS IS THE FIRST LEVEL OF THE LOOP.
C
         IF ( LEVEL .EQ. 0 ) THEN
            DOLOOP      = .TRUE.
            NARRAY      = 0
            NINSTR( 1 ) = 0
            NINSTR( 2 ) = 0
            NINSTR( 3 ) = 0
C
C  THIS IS THE SECOND OR THIRD LEVEL OF THE LOOP.
C
         ELSE
            NINSTR( LEVEL ) = NINSTR( LEVEL ) + 1
            IF ( NINSTR( LEVEL ) .GT. MAXINS ) THEN
               INFORM = - 11
               GO TO 700
            END IF
            IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4010 )
     *         LEVEL, NINSTR( LEVEL )
            INSTR( 1, NINSTR( LEVEL ), LEVEL ) = 1
         END IF
C
C  RECORD THE LOCATION OF THE DO-LOOP VARIABLE IN THE ARRAY INLIST.
C
         FIELD = FIELD2( 1 : 10 ) // 'II'
         CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) GO TO 700
            IFREE = - IFREE
         ELSE
            NUSEIN = NUSEIN + 1
            IF ( NUSEIN .GT. NINDEX ) THEN
               INFORM = - 21
               GO TO 700
            END IF
            INLIST( IFREE )  = NUSEIN
            NAMIIN( NUSEIN ) = FIELD( 1: 7 )
         END IF
         IF ( LEVEL .EQ. 0 ) THEN
            LOOP( 1 ) = INLIST( IFREE )
         ELSE
            INSTR( 2, NINSTR( LEVEL ), LEVEL ) = INLIST( IFREE )
         END IF
C
C  RECORD THE STARTING VALUE OF THE DO-LOOP VARIABLE.
C
         FIELD = FIELD3( 1 : 10 ) // 'II'
         CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
         IF ( IFIELD .LE. 0 ) THEN
            INFORM = 3
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1: 10 )
            GO TO 800
         END IF
         IF ( LEVEL .EQ. 0 ) THEN
            LOOP( 2 ) = INLIST( IFIELD )
         ELSE
            INSTR( 3, NINSTR( LEVEL ), LEVEL ) = INLIST( IFIELD )
         END IF
C
C  RECORD THE FINISHING VALUE OF THE DO-LOOP VARIABLE AND
C  SET THE INCREMENTAL VALUE OF THE DO-LOOP VARIABLE TO 1.
C
         FIELD = FIELD5( 1 : 10 ) // 'II'
         CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
         IF ( IFIELD .LE. 0 ) THEN
            INFORM = 3
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
            GO TO 800
         END IF
         IF ( LEVEL .EQ. 0 ) THEN
            LOOP( 3 ) = INLIST( IFIELD )
            LOOP( 4 ) = - 1
         ELSE
            INSTR( 4, NINSTR( LEVEL ), LEVEL ) = INLIST( IFIELD )
            INSTR( 5, NINSTR( LEVEL ), LEVEL ) = - 1
         END IF
         LEVEL = LEVEL + 1
         GO TO 100
      END IF
C
C  A DO-LOOP VARIABLE IS TO HAVE A NON-TRIVIAL INCREMENT.
C
      IF ( FIELD1 .EQ. 'DI' ) THEN
C
C  RECORD THE LOCATION OF THE DO-LOOP VARIABLE IN THE ARRAY INLIST.
C
         IF ( ( LEVEL .EQ. 1 .AND.  FIELD2( 1: 10 ) .EQ.
     *        NAMIIN( LOOP( 1 ) ) ) .OR.( LEVEL .GT. 1 .AND.
     *        FIELD2( 1: 10 ) .EQ. NAMIIN(
     *        INSTR( 2, NINSTR( LEVEL - 1 ), LEVEL - 1 ) ) ) ) THEN
            IF ( DEBUG .AND. IOUT .GT. 0 .AND. LEVEL .GT. 1 )
     *         WRITE( IOUT, 4030 ) LEVEL - 1, NINSTR( LEVEL - 1 )
            FIELD = FIELD3( 1 : 10 ) // 'II'
            CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 3
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
               GO TO 800
            END IF
            IF ( LEVEL .EQ. 1 ) THEN
               LOOP( 4 ) = INLIST( IFIELD )
            ELSE
               INSTR( 5, NINSTR( LEVEL - 1 ), LEVEL - 1 ) =
     *            INLIST( IFIELD )
            END IF
         END IF
         GO TO 100
      END IF
C
C  END OF ONE OR MORE DO-LOOPS.
C  ============================
C
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C     IF ( FIELD1 .EQ. 'OD' .OR. FIELD1 .EQ. 'ND' ) THEN
      IF ( FIELD1 .NE. 'OD' .AND. FIELD1 .NE. 'ND' ) GO TO 341
C
C  TERMINATE THE CURRENT LEVEL OF LOOP.
C
         IF ( FIELD1 .EQ. 'OD' ) THEN
            NINSTR( LEVEL ) = NINSTR( LEVEL ) + 1
            IF ( NINSTR( LEVEL ) .GT. MAXINS ) THEN
               INFORM = - 11
               GO TO 700
            END IF
            INSTR( 1, NINSTR( LEVEL ), LEVEL ) = 2
            IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4020 )
     *         LEVEL, NINSTR( LEVEL )
            LEVEL = LEVEL - 1
         ELSE
            DO 210 I = LEVEL, 1, - 1
               NINSTR( LEVEL ) = NINSTR( LEVEL ) + 1
               IF ( NINSTR( LEVEL ) .GT. MAXINS ) THEN
                  INFORM = - 11
                  GO TO 700
               END IF
               INSTR( 1, NINSTR( LEVEL ), LEVEL ) = 2
               IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4020 )
     *            LEVEL, NINSTR( LEVEL )
               LEVEL = LEVEL - 1
  210       CONTINUE
         END IF
C
C  EXECUTE DO-LOOP INSTRUCTIONS.
C  =============================
C
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C        IF ( LEVEL .EQ. 0 ) THEN
         IF ( LEVEL .NE. 0 ) GO TO 339
C
C  EXECUTE LEVEL-1 DO-LOOP INSTRUCTIONS.
C
            LEV1S = INDVAL( LOOP( 2 ) )
            LEV1E = INDVAL( LOOP( 3 ) )
            IF ( LOOP( 4 ) .LE. 0 ) THEN
               LEV1I = 1
            ELSE
               LEV1I = INDVAL( LOOP( 4 ) )
            END IF
C
C  MOCK DO-LOOP.
C
            LEV1 = LEV1S
  220       CONTINUE
C ** Correction 11  30/06/04: Assigned goto statements replaced. 2 lines replaced
C           IF ( ( LEV1I .GT. 0 .AND. LEV1 .LE. LEV1E ) .OR.
C    *           ( LEV1I .LT. 0 .AND. LEV1 .GE. LEV1E ) ) THEN
            IF ( .NOT. ( ( LEV1I .GT. 0 .AND. LEV1 .LE. LEV1E ) .OR.
     *           ( LEV1I .LT. 0 .AND. LEV1 .GE. LEV1E ) ) ) GO TO 337
               L1 = 0
               L2 = 0
               L3 = 0
C
C  SET THE LOOP VALUE.
C
               INDVAL( LOOP( 1 ) ) = LEV1
               IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 5000 ) 1,
     *            NAMIIN( LOOP( 1 ) ), LEV1
C
C  EXECUTE THE REMAINING LIST OF LEVEL-1 INSTRUCTIONS.
C
  230          CONTINUE
               L1 = L1 + 1
C
C  SEE IF THE LEVEL-1 LOOP IS TO BE TERMINATED.
C
               IF ( INSTR( 1, L1, 1 ) .EQ. 2 ) GO TO 300
C
C  SEE IF A LEVEL-2 LOOP IS TO BE STARTED.
C
C              IF ( INSTR( 1, L1, 1 ) .EQ. 1 ) THEN
               IF ( INSTR( 1, L1, 1 ) .NE. 1 ) GO TO 295
C
C  EXECUTE LEVEL-2 DO-LOOP INSTRUCTIONS.
C
                  LEV2S  = INDVAL( INSTR( 3, L1, 1 ) )
                  LEV2E  = INDVAL( INSTR( 4, L1, 1 ) )
                  IF ( INSTR( 5, L1, 1 ) .LE. 0 ) THEN
                     LEV2I = 1
                  ELSE
                     LEV2I = INDVAL( INSTR( 5, L1, 1 ) )
                  END IF
                  LEVEL2 = L2
                  LEVL3A = L3
C
C  MOCK DO-LOOP.
C
                  LEV2 = LEV2S
  240             CONTINUE
                  L2 = LEVEL2
C ** Correction 1. 30/11/93: 3 lines interchanged **
C ** Correction 11  30/06/04: Assigned goto statements replaced. 2 lines replaced
C                 IF ( ( LEV2I .GT. 0 .AND. LEV2 .LE. LEV2E ) .OR.
C    *                 ( LEV2I .LT. 0 .AND. LEV2 .GE. LEV2E ) ) THEN
                  IF ( .NOT. ( LEV2I .GT. 0 .AND. LEV2 .LE. LEV2E ) .OR.
     *               ( LEV2I .LT. 0 .AND. LEV2 .GE. LEV2E ) ) GO TO 292
                     L3 = LEVL3A
C ** Correction 1. 30/11/93: end of correction **
C
C  SET THE LOOP VALUE.
C
                     INDVAL( INSTR( 2, L1, 1 ) ) = LEV2
                     IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 5000 )
     *                  2, NAMIIN( INSTR( 2, L1, 1 ) ), LEV2
C
C  EXECUTE THE REMAINING LIST OF LEVEL-2 INSTRUCTIONS.
C
  250                CONTINUE
                     L2 = L2 + 1
C
C  SEE IF THE LEVEL-2 LOOP IS TO BE TERMINATED.
C
                     IF ( INSTR( 1, L2, 2 ) .EQ. 2 ) GO TO 290
C
C  SEE IF A LEVEL-3 LOOP IS TO BE STARTED.
C
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C                    IF ( INSTR( 1, L2, 2 ) .EQ. 1 ) THEN
                     IF ( INSTR( 1, L2, 2 ) .NE. 1 ) GO TO 283
C
C  EXECUTE LEVEL-3 DO-LOOP INSTRUCTIONS.
C
                        LEV3S  = INDVAL( INSTR( 3, L2, 2 ) )
                        LEV3E  = INDVAL( INSTR( 4, L2, 2 ) )
                        IF ( INSTR( 5, L2, 2 ) .LE. 0 ) THEN
                           LEV3I = 1
                        ELSE
                           LEV3I = INDVAL( INSTR( 5, L2, 2 ) )
                        END IF
                        LEVEL3 = L3
C
C  MOCK DO-LOOP.
C
                        LEV3 = LEV3S
  260                   CONTINUE
                        L3 = LEVEL3
C ** Correction 11  30/06/04: Assigned goto statements replaced. 2 lines replaced
C                       IF ( ( LEV3I .GT. 0 .AND. LEV3 .LE. LEV3E ) .OR.
C    *                     ( LEV3I .LT. 0 .AND. LEV3 .GE. LEV3E ) ) THEN
                        IF ( .NOT.
     *                     ( LEV3I .GT. 0 .AND. LEV3 .LE. LEV3E ) .OR.
     *                     ( LEV3I .LT. 0 .AND. LEV3 .GE. LEV3E ) ) 
     *                       GO TO 281
C
C  SET THE LOOP VALUE.
C
                           INDVAL( INSTR( 2, L2, 2 ) ) = LEV3
                           IF ( DEBUG .AND. IOUT .GT. 0 )
     *                        WRITE( IOUT, 5000 ) 3,
     *                           NAMIIN( INSTR( 2, L2, 2 ) ), LEV3
C
C  EXECUTE THE REMAINING LIST OF LEVEL-3 INSTRUCTIONS.
C
  270                      CONTINUE
                           L3 = L3 + 1
C
C  SEE IF THE LEVEL-3 LOOP IS TO BE TERMINATED.
C
                           IF ( INSTR( 1, L3, 3 ) .EQ. 2 ) GO TO 280
C
C  EXECUTE LEVEL-3 INDEX INSTRUCTIONS.
C
                           IF ( INSTR( 1, L3, 3 ) .GE. 21 .AND.
     *                          INSTR( 1, L3, 3 ) .LE. 50 ) THEN
                              CALL GETIIN( NINDEX, INDVAL, NRLNDX,
     *                                     REALVL, INSTR( 1, L3, 3 ) )
                              IF ( DEBUG .AND. IOUT .GT. 0 )
     *                           WRITE( IOUT, 5010 ) 3, L3,
     *                           NAMIIN( INSTR( 2, L3, 3 ) ),
     *                           INDVAL( INSTR( 2, L3, 3 ) )
                           END IF
                           IF ( INSTR( 1, L3, 3 ) .GE. 51 .AND.
     *                          INSTR( 1, L3, 3 ) .LE. 99 ) THEN
                              CALL GETRIN( NINDEX, NRLNDX, INDVAL,
     *                                     REALVL, RVALUE( L3, 3 ),
     *                                     INSTR( 1, L3, 3 ), INFORM )
                              IF ( INFORM .GT. 0 ) GO TO 800
                              IF ( DEBUG .AND. IOUT .GT. 0 )
     *                           WRITE( IOUT, 5020 ) 3, L3,
     *                           NAMRIN( INSTR( 2, L3, 3 ) ),
     *                           REALVL( INSTR( 2, L3, 3 ) )
                           END IF
                           IF ( INSTR( 1, L3, 3 ) .GE. 100 ) THEN
                              NARRAY = INSTR( 2, L3, 3 )
                              CALL GETLIN( NINDEX, NRLNDX, INDVAL,
     *                                     IARRAY( 1, 1, NARRAY ),
     *                                     VARRAY( 1, NARRAY ),
     *                                     ARRAY( 1, NARRAY ),
     *                                     CARRAY( 1, NARRAY ),
     *                                     FARRAY( NARRAY ), REALVL,
     *                                     NAMIIN, NOVALS,
     *                                     INSTR( 1, L3, 3 ), FIELD1,
     *                                     FIELD2, FIELD3, VALUE4,
     *                                     FIELD5, VALUE6, IOUT, INFORM,
     *                                     LENGTH, KEY, ITABLE, INLIST )
                              IF ( INFORM .GT. 0 ) GO TO 800
                              IF ( INFORM .LT. 0 ) GO TO 700
                              IF ( DEBUG .AND. IOUT .GT. 0 )
     *                           WRITE( IOUT, 5060 )
     *                               3, L3, FIELD1, FIELD2,
     *                               FIELD3, VALUE4, FIELD5, VALUE6
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C                             ASSIGN 270 TO IJUMP
                              IJUMP = 3
                              GO TO 400
                           END IF
                           GO TO 270
  280                      CONTINUE
                           LEV3 = LEV3 + LEV3I
                           GO TO 260
C                       ELSE
C
C  THE DO-LOOP IS NOT EXECUTED. FIND THE NEXT RELEVANT INSTRUCTION.
C
  281                      CONTINUE
                           L3 = L3 + 1
                           IF ( INSTR( 1, L3, 3 ) .NE. 2 ) GO TO 281
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C                       END IF
C
C  END OF LEVEL-3 DO-LOOP.
C
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C                    END IF
  283                CONTINUE
C
C  EXECUTE LEVEL-2 INDEX INSTRUCTIONS.
C
                     IF ( INSTR( 1, L2, 2 ) .GE. 21 .AND.
     *                    INSTR( 1, L2, 2 ) .LE. 50 ) THEN
                        CALL GETIIN( NINDEX, INDVAL, NRLNDX, REALVL,
     *                               INSTR( 1, L2, 2 ) )
                        IF ( DEBUG .AND. IOUT .GT. 0 )
     *                     WRITE( IOUT, 5010 ) 2, L2,
     *                        NAMIIN( INSTR( 2, L2, 2 ) ),
     *                        INDVAL( INSTR( 2, L2, 2 ) )
                     END IF
                     IF ( INSTR( 1, L2, 2 ) .GE. 51 .AND.
     *                    INSTR( 1, L2, 2 ) .LE. 99 ) THEN
                        CALL GETRIN( NINDEX, NRLNDX, INDVAL,
     *                               REALVL, RVALUE( L2, 2 ),
     *                               INSTR( 1, L2, 2 ), INFORM )
                        IF ( INFORM .GT. 0 ) GO TO 800
                        IF ( DEBUG .AND. IOUT .GT. 0 )
     *                     WRITE( IOUT, 5020 ) 2, L2,
     *                     NAMRIN( INSTR( 2, L2, 2 ) ),
     *                     REALVL( INSTR( 2, L2, 2 ) )
                     END IF
                     IF ( INSTR( 1, L2, 2 ) .GE. 100 ) THEN
                        NARRAY = INSTR( 2, L2, 2 )
                        CALL GETLIN( NINDEX, NRLNDX, INDVAL,
     *                               IARRAY( 1, 1, NARRAY ),
     *                               VARRAY( 1, NARRAY ),
     *                               ARRAY( 1, NARRAY ),
     *                               CARRAY( 1, NARRAY ),
     *                               FARRAY( NARRAY ), REALVL,
     *                               NAMIIN, NOVALS,
     *                               INSTR( 1, L2, 2 ), FIELD1,
     *                               FIELD2, FIELD3, VALUE4,
     *                               FIELD5, VALUE6, IOUT, INFORM,
     *                               LENGTH, KEY, ITABLE, INLIST )
                        IF ( INFORM .GT. 0 ) GO TO 800
                        IF ( INFORM .LT. 0 ) GO TO 700
                        IF ( DEBUG .AND. IOUT .GT. 0 )
     *                     WRITE( IOUT, 5060 )
     *                        2, L2, FIELD1, FIELD2, FIELD3,
     *                        VALUE4, FIELD5, VALUE6
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C                       ASSIGN 250 TO IJUMP
                        IJUMP = 2
                        GO TO 400
                     END IF
                     GO TO 250
  290                CONTINUE
                     LEV2 = LEV2 + LEV2I
                     GO TO 240
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C                 ELSE
C
C  THE DO-LOOP IS NOT EXECUTED. FIND THE NEXT RELEVANT INSTRUCTION.
C
  292                CONTINUE
                     L2 = L2 + 1
                     IF ( INSTR( 1, L2, 2 ) .NE. 2 ) GO TO 292
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C                 END IF
  293             CONTINUE
                  LEVEL2 = L2
C
C  END OF LEVEL-2 DO-LOOP.
C
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C              END IF
  295          CONTINUE
C
C  EXECUTE LEVEL-1 INDEX INSTRUCTIONS.
C
               IF ( INSTR( 1, L1, 1 ) .GE. 21 .AND.
     *              INSTR( 1, L1, 1 ) .LE. 50 ) THEN
                  CALL GETIIN( NINDEX, INDVAL, NRLNDX, REALVL,
     *                         INSTR( 1, L1, 1 ) )
                  IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 5010 )
     *               1, L1, NAMIIN( INSTR( 2, L1, 1 ) ),
     *               INDVAL( INSTR( 2, L1, 1 ) )
               END IF
               IF ( INSTR( 1, L1, 1 ) .GE. 51 .AND.
     *              INSTR( 1, L1, 1 ) .LE. 99 ) THEN
                  CALL GETRIN( NINDEX, NRLNDX, INDVAL,
     *                         REALVL, RVALUE( L1, 1 ),
     *                         INSTR( 1, L1, 1 ), INFORM )
                  IF ( INFORM .GT. 0 ) GO TO 800
                  IF ( DEBUG .AND. IOUT .GT. 0 )
     *               WRITE( IOUT, 5020 ) 1, L1,
     *               NAMRIN( INSTR( 2, L1, 1 ) ),
     *               REALVL( INSTR( 2, L1, 1 ) )
               END IF
               IF ( INSTR( 1, L1, 1 ) .GE. 100 ) THEN
                  NARRAY = INSTR( 2, L1, 1 )
                  CALL GETLIN( NINDEX, NRLNDX, INDVAL,
     *                         IARRAY( 1, 1, NARRAY ),
     *                         VARRAY( 1, NARRAY ), ARRAY( 1, NARRAY ),
     *                         CARRAY( 1, NARRAY ),
     *                         FARRAY( NARRAY ), REALVL,
     *                         NAMIIN, NOVALS,
     *                         INSTR( 1, L1, 1 ), FIELD1,
     *                         FIELD2, FIELD3, VALUE4,
     *                         FIELD5, VALUE6, IOUT, INFORM,
     *                         LENGTH, KEY, ITABLE, INLIST )
                  IF ( INFORM .GT. 0 ) GO TO 800
                  IF ( INFORM .LT. 0 ) GO TO 700
                  IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 5060 )
     *                             1, L1, FIELD1, FIELD2, FIELD3,
     *                             VALUE4, FIELD5, VALUE6
C                 ASSIGN 230 TO IJUMP
                  IJUMP = 1
                  GO TO 400
               END IF
               GO TO 230
  300          CONTINUE
               LEV1 = LEV1 + LEV1I
               GO TO 220
  337       CONTINUE   
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C           END IF
C
C  END OF LEVEL-1 DO-LOOP.
C
            DOLOOP = .FALSE.
  339    CONTINUE
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C        END IF
         GO TO 100
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
C     END IF
  341 CONTINUE
C
C  CONSTRUCT A LIST OF DO-LOOP INSTRUCTIONS: 1) ARITHMETIC INSTRUCTIONS.
C  =====================================================================
C
      IF ( DOLOOP ) THEN
         NINSTR( LEVEL ) = NINSTR( LEVEL ) + 1
         IF ( NINSTR( LEVEL ) .GT. MAXINS ) THEN
            INFORM = - 11
            GO TO 700
         END IF
C
C  AN ARITHMETIC INSTRUCTION IS TO BE PROCESSED.
C
         IF ( FIELD1 .EQ. 'IE' .OR. FIELD1 .EQ. 'IA' .OR.
     *        FIELD1 .EQ. 'IS' .OR. FIELD1 .EQ. 'IM' .OR.
     *        FIELD1 .EQ. 'ID' .OR. FIELD1 .EQ. 'IR' .OR.
     *        FIELD1 .EQ. 'I=' .OR. FIELD1 .EQ. 'I+' .OR.
     *        FIELD1 .EQ. 'I-' .OR. FIELD1 .EQ. 'I*' .OR.
     *        FIELD1 .EQ. 'I/' .OR. FIELD1 .EQ. 'RE' .OR.
     *        FIELD1 .EQ. 'RA' .OR. FIELD1 .EQ. 'RS' .OR.
     *        FIELD1 .EQ. 'RM' .OR. FIELD1 .EQ. 'RD' .OR.
     *        FIELD1 .EQ. 'RI' .OR. FIELD1 .EQ. 'RF' .OR.
     *        FIELD1 .EQ. 'R=' .OR. FIELD1 .EQ. 'R+' .OR.
     *        FIELD1 .EQ. 'R-' .OR. FIELD1 .EQ. 'R*' .OR.
     *        FIELD1 .EQ. 'R/' .OR. FIELD1 .EQ. 'R(' ) THEN
            CALL PROCAI( NINDEX, NRLNDX, LENGTH, NUSEIN, NUSERE,
     *                   INFORM, IOUT, LEVEL, NINSTR( LEVEL ),
     *                   DEBUG, RVALUE( NINSTR( LEVEL ), LEVEL ),
     *                   INLIST, ITABLE, NAMIIN, NAMRIN,
     *                   INSTR( 1, NINSTR( LEVEL ), LEVEL ), KEY,
     *                   FIELD1, FIELD2, FIELD3, FIELD5,
     *                   NULINE( 25: 36 ) )
            IF ( INFORM .GT. 0 ) GO TO 800
            IF ( INFORM .LT. 0 ) GO TO 700
         ELSE
C
C  CONSTRUCT A LIST OF DO-LOOP INSTRUCTIONS: 2) ARRAY DEFINITIONS.
C  ===============================================================
C
            IF ( FIELD1( 1: 1 ) .NE. 'X' .AND. FIELD1( 1: 1 ) .NE. 'Z'
     *           .AND. FIELD1( 1: 1 ) .NE. 'A' ) THEN
               INFORM = 6
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2060 )
               GO TO 800
            END IF
            NARRAY = NARRAY + 1
            IF ( NARRAY .GT. MAXARA ) THEN
               INFORM = - 12
               GO TO 700
            END IF
            CALL PROCAD( NINDEX, NRLNDX, LEVEL, NINSTR( LEVEL ),
     *                   NUSERE, LENGTH, NARRAY, INTYPE, INFORM, IOUT,
     *                   DEBUG, GRP1ST,
     *                   FIELD1, FIELD2, FIELD3, FIELD5,
     *                   NULINE( 25 : 36 ), NULINE( 50 : 61 ),
     *                   INLIST,
     *                   INSTR( 1, NINSTR( LEVEL ), LEVEL ), ITABLE,
     *                   IARRAY( 1, 1, NARRAY ), VARRAY( 1, NARRAY ),
     *                   FARRAY( NARRAY ), NAMIIN, NAMRIN,
     *                   ARRAY( 1, NARRAY ), CARRAY( 1, NARRAY ), KEY )
            IF ( INFORM .GT. 0 ) GO TO 800
            IF ( INFORM .LT. 0 ) GO TO 700
C
C  THE ARRAY DEFINITION IS COMPLETE.
C
         END IF
         GO TO 100
      ELSE
C
C  EXECUTE A NON-DO-LOOP INSTRUCTION.
C  ==================================
C
C  THE INSTRUCTION IS AN ARRAY INSTRUCTION.
C
         IF ( FIELD1 .EQ. 'IE' .OR. FIELD1 .EQ. 'IA' .OR.
     *        FIELD1 .EQ. 'IS' .OR. FIELD1 .EQ. 'IM' .OR.
     *        FIELD1 .EQ. 'ID' .OR. FIELD1 .EQ. 'IR' .OR.
     *        FIELD1 .EQ. 'I=' .OR. FIELD1 .EQ. 'I+' .OR.
     *        FIELD1 .EQ. 'I-' .OR. FIELD1 .EQ. 'I*' .OR.
     *        FIELD1 .EQ. 'I/' .OR. FIELD1 .EQ. 'RE' .OR.
     *        FIELD1 .EQ. 'RA' .OR. FIELD1 .EQ. 'RS' .OR.
     *        FIELD1 .EQ. 'RM' .OR. FIELD1 .EQ. 'RD' .OR.
     *        FIELD1 .EQ. 'RI' .OR. FIELD1 .EQ. 'RF' .OR.
     *        FIELD1 .EQ. 'R=' .OR. FIELD1 .EQ. 'R+' .OR.
     *        FIELD1 .EQ. 'R-' .OR. FIELD1 .EQ. 'R*' .OR.
     *        FIELD1 .EQ. 'R/' .OR. FIELD1 .EQ. 'R(' ) THEN
C
C  1) AN ARITHMETIC INSTRUCTION. DECODE THE INSTRUCTION.
C
            CALL PROCAI( NINDEX, NRLNDX, LENGTH, NUSEIN, NUSERE,
     *                   INFORM, IOUT, 0, 1, DEBUG, RVALUE( 1, 1 ),
     *                   INLIST, ITABLE, NAMIIN, NAMRIN,
     *                   INSTR( 1, 1, 1 ), KEY,
     *                   FIELD1, FIELD2, FIELD3, FIELD5,
     *                   NULINE( 25: 36 ) )
            IF ( INFORM .GT. 0 ) GO TO 800
            IF ( INFORM .LT. 0 ) GO TO 700
C
C  EXECUTE THE INSTRUCTION.
C
            IF ( FIELD1( 1 : 1 ) .EQ. 'I' ) THEN
               CALL GETIIN( NINDEX, INDVAL, NRLNDX, REALVL,
     *                      INSTR( 1, 1, 1 ) )
               IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 5010 ) 0, 1,
     *            NAMIIN( INSTR( 2, 1, 1 ) ), INDVAL( INSTR( 2, 1, 1 ) )
            ELSE
               CALL GETRIN( NINDEX, NRLNDX, INDVAL, REALVL,
     *                      RVALUE( 1, 1 ), INSTR( 1, 1, 1 ), INFORM )
               IF ( INFORM .GT. 0 ) GO TO 800
               IF ( DEBUG .AND. IOUT .GT. 0 )
     *            WRITE( IOUT, 5020 ) 0, 1, NAMRIN( INSTR( 2, 1, 1 ) ),
     *            REALVL( INSTR( 2, 1, 1 ) )
            END IF
            GO TO 100
         ELSE
            IF ( FIELD1( 1 : 1 ) .EQ. 'X' .OR. FIELD1( 1: 1 ) .EQ. 'Z'
     *           .OR. FIELD1( 1 : 1 ) .EQ. 'A' ) THEN
C
C  2) AN ARRAY DEFINITION. DECODE THE INSTRUCTION.
C
               CALL PROCAD( NINDEX, NRLNDX, 0, 1, NUSERE,
     *                      LENGTH, 1, INTYPE, INFORM, IOUT,
     *                      DEBUG, GRP1ST,
     *                      FIELD1, FIELD2, FIELD3, FIELD5,
     *                      NULINE( 25 : 36 ), NULINE( 50 : 61 ),
     *                      INLIST, INSTR( 1, 1, 1 ), ITABLE,
     *                      IARRAY( 1, 1, 1 ), VARRAY( 1, 1 ),
     *                      FARRAY( 1 ), NAMIIN, NAMRIN,
     *                      ARRAY( 1, 1 ), CARRAY( 1, 1 ), KEY )
               IF ( INFORM .GT. 0 ) GO TO 800
               IF ( INFORM .LT. 0 ) GO TO 700
C
C  EXECUTE THE INSTRUCTION.
C
               CALL GETLIN( NINDEX, NRLNDX, INDVAL,
     *                      IARRAY( 1, 1, 1 ), VARRAY( 1, 1 ),
     *                      ARRAY( 1, 1 ), CARRAY( 1, 1 ),
     *                      FARRAY( 1 ), REALVL, NAMIIN, NOVALS,
     *                      INSTR( 1, 1, 1 ), FIELD1,
     *                      FIELD2, FIELD3, VALUE4,
     *                      FIELD5, VALUE6, IOUT, INFORM,
     *                      LENGTH, KEY, ITABLE, INLIST )
               IF ( INFORM .GT. 0 ) GO TO 800
               IF ( INFORM .LT. 0 ) GO TO 700
               IF ( DEBUG .AND. IOUT .GT. 0 )
     *            WRITE( IOUT, 5060 ) 0, 1, FIELD1, FIELD2,
     *                            FIELD3, VALUE4, FIELD5, VALUE6
               GO TO 400
            ELSE
C
C  THE INSTRUCTION IS NOT AN ARRAY INSTRUCTION.
C  CHECK TO SEE IF THERE IS ARE ANY NUMERICAL VALUES TO BE READ.
C
               NOVALS = 0
               IF ( ( FIELD3 .NE. '          ' .AND.
     *                NULINE( 15: 15 ) .NE. '$' )
C    *                .OR. FIELD1( 1: 1 ) .EQ. 'D'
C    *                .OR. FIELD2 .EQ. '''DEFAULT'' '
C    *                .OR. FIELD2 .EQ. '''default'' '
C    *                .OR. FIELD3 .EQ. '''DEFAULT'' '
C    *                .OR. FIELD3 .EQ. '''default'' '
     *                .OR. INTYPE .EQ. MOBBND ) THEN
                  NOVALS = NOVALS + 1
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
                  IF (   INTYPE .LE. MQHESS .OR. INTYPE .EQ. MOBBND .OR.
     *               ( ( INTYPE .EQ. MEUSES .OR. INTYPE .EQ. MGUSES )
     *                   .AND. FIELD1 .EQ. 'P ' ) .OR. ( INTYPE .EQ.
     *                         MGUSES .AND. FIELD1 .EQ. 'E ' ) ) THEN
                     IF ( INTYPE .EQ. MGUSES .AND. FIELD1 .EQ. 'E '.AND.
     *                 NULINE( 25: 36 ) .EQ. '            ' ) THEN
                       VALUE4 = ONE
                     ELSE
                        CALL GETVAL( NULINE( 25: 36 ), VALUE4 )
                     END IF
                     IF ( INTYPE .EQ. MRANGE ) VALUE4 = ABS( VALUE4 )
                  END IF
                  IF ( FIELD5 .NE. '          ' .AND.
     *                 NULINE( 40: 40 ) .NE. '$' ) THEN
                     NOVALS = NOVALS + 1
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
                     IF (  INTYPE .LE. MQHESS .OR.
     *                 ( ( INTYPE .EQ. MEUSES .OR. INTYPE .EQ. MGUSES )
     *                      .AND. FIELD1 .EQ. 'P ' ) .OR. ( INTYPE .EQ.
     *                           MGUSES .AND. FIELD1 .EQ. 'E ' ) ) THEN
                        IF ( INTYPE .EQ. MGUSES .AND.
     *                       FIELD1 .EQ. 'E '   .AND.
     *                       NULINE( 50: 61 ) .EQ. '            ' ) THEN
                          VALUE6 = ONE
                        ELSE
                           CALL GETVAL( NULINE( 50: 61 ), VALUE6 )
                        END IF
                        IF ( INTYPE .EQ. MRANGE ) VALUE6 = ABS( VALUE6 )
                        IF ( INTYPE .LT. MCONST ) THEN
                           IF ( VALUE6 .EQ. ZERO ) NOVALS = 1
                        END IF
                     END IF
C
C  REMOVE FIELDS WITH NUMERICAL VALUES OF ZERO.
C
                  END IF
                  IF ( INTYPE .LT. MCONST ) THEN
                     IF ( VALUE4 .EQ. ZERO ) THEN
                        IF ( NOVALS .EQ. 2 ) THEN
                           VALUE4 = VALUE6
                           FIELD3 = FIELD5
                        END IF
                       NOVALS = NOVALS - 1
                     END IF
                  END IF
                  IF ( FIELD3 .EQ. '''SCALE''   ' .OR.
     *                 FIELD3 .EQ. ' ''SCALE''  ' ) THEN
                     NOVALS = 0
                     VALUE4 = ABS( VALUE4 )
                  END IF
               END IF
               IF ( DEBUG .AND. IOUT .GT. 0 )
     *            WRITE( IOUT, 5060 ) 0, 1, FIELD1, FIELD2,
     *                            FIELD3, VALUE4, FIELD5, VALUE6
            END IF
         END IF
      END IF
  400 CONTINUE
C
C  EXECUTE REAL PARAMETER ARRAY CARD.
C
      IF ( FIELD1 .EQ. 'AE' .OR. FIELD1 .EQ. 'AA' .OR.
     *     FIELD1 .EQ. 'AS' .OR. FIELD1 .EQ. 'AM' .OR.
     *     FIELD1 .EQ. 'AD' .OR. FIELD1 .EQ. 'AI' .OR.
     *     FIELD1 .EQ. 'AF' .OR. FIELD1 .EQ. 'A=' .OR.
     *     FIELD1 .EQ. 'A+' .OR. FIELD1 .EQ. 'A-' .OR.
     *     FIELD1 .EQ. 'A*' .OR. FIELD1 .EQ. 'A/' .OR.
     *     FIELD1 .EQ. 'A(' ) THEN
         CALL    PROCAA( NINDEX, LENGTH, NUSERE, INFORM, IOUT,
C ** Correction 8. 26/02/01: 1 dummy argument removed from PROCAA **
     *                   NRLNDX, INLIST, ITABLE,
     *                   NAMRIN, KEY, INDVAL, REALVL,
     *                   FIELD1, FIELD2, FIELD3, FIELD5, VALUE4 )
         IF ( INFORM .GT. 0 ) GO TO 800
         IF ( INFORM .LT. 0 ) GO TO 700
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
C
C  BRANCH DEPENDING ON THE CURRENT INDICATOR CARD.
C  ===============================================
C
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C ** Correction -2. 21/02/00: QSECTION added as alias for QUADOBJ
      GO TO ( 100, 100, 100, 100, 420, 420, 420, 430, 430, 430,
     *        430, 430, 430, 440, 450, 455, 455, 455, 455, 455, 
     *        460, 470, 480, 490, 500, 900 ), INTYPE
C
C  INDICATOR CARD IS GROUPS/ROWS/CONSTRAINTS.
C  ------------------------------------------
C
  420 CONTINUE
      IF ( GRP1ST ) THEN
         CALL    SGRP1 ( NG, NGMAX, NOMAX, LENGTH, NOBJ, NOVALS,
     *                   INLIST, IDROWS, ISTATE, ITYPEG, ITABLE,
     *                   RDROWS, RSCALE,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   GNAMES, ONAMES, KEY, IOUT, INFORM )
      ELSE
         CALL    SGRP2 ( NG, NGMAX, NOMAX, LA, LENGTH,
     *                   NNZA, NOBJ, NOVALS, INLIST, IDROWS,
     *                   ISTATE, ITYPEG, ITABLE, ICOORD,
     *                   A, RDROWS, RSCALE,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   GNAMES, ONAMES, KEY, IOUT, INFORM )
      END IF
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C
C  INDICATOR CARD IS COLUMNS/VARIABLES, CONSTANTS/RHS/RHS' OR RANGES.
C  ------------------------------------------------------------------
C
  430 CONTINUE
      IF ( INTYPE .EQ. MCOLS  ) COLFIE = 'VA'
      IF ( INTYPE .EQ. MCONST ) COLFIE = 'CO'
      IF ( INTYPE .EQ. MRANGE ) COLFIE = 'RA'
      IF ( GRP1ST .OR. INTYPE .NE. MCOLS ) THEN
         CALL    SVAR2 ( NMAX, NGMAX, LA, LENGTH, NNZA, NVAR, NOVALS,
     *                   NRLNDX, INTYPE .EQ. MCOLS, COLFIE,
C ** Correction 0. 20/12/99: Array holding variable types introduced.
     *                   ICOORD, ISTATE, ITYPEV, INLIST, ITABLE,
     *                   A, CSCALE, REALVL, DFAULT,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   VNAMES, KEY, IOUT, INFORM )
      ELSE
         CALL    SVAR1 ( NMAX, LENGTH, NVAR, COLFIE,
C ** Correction 0. 20/12/99: Array holding variable types introduced.
     *                   ITYPEV, INLIST, ITABLE, CSCALE,
C ** Correction 2. 26/02/01: 4 dummy arguments removed from SVAR1 **
     *                   FIELD2, FIELD3, VALUE4, 
     *                   VNAMES, KEY, INFORM )
      END IF
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C
C  INDICATOR CARD IS BOUNDS.
C  -------------------------
C
  440 CONTINUE
      CALL       SBOUND( NMAX, NBMAX, LENGTH, NLVARS, NBND, NCOL,
     *                   NRLNDX, DEFAUT, INLIST, ITABLE, BND, REALVL,
C ** Correction 3. 26/02/01: 1 dummy argument removed from SBOUND **
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, 
     *                   BNAMES, BNDFLT, KEY, IOUT, INFORM )
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C
C  INDICATOR CARD IS START POINT.
C  ------------------------------
C
  450 CONTINUE
      CALL       SSTART( NMAX, NGMAX, NSMAX, LENGTH, NLVARS, NG, NSTART,
     *                   NCOL, NRLNDX, DEFAUT, INLIST, ITABLE, VSTART,
     *                   CSTART, REALVL, FIELD1, FIELD2, FIELD3, VALUE4,
     *                   FIELD5, VALUE6, SNAMES, NOVALS, KEY,
     *                   IOUT, INFORM )
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C
C  INDICATOR CARD IS QHESS.
C  ------------------------
C
  455 CONTINUE
      CALL       SQHESS( NEGMAX, NGMAX, NLMAX, NELMAX, NEVMAX, NETMAX, 
     *                   NOMAX, NIMAX, LENGTH, NG, NOBJ, NGRUPE, NOVALS, 
     *                   NELN, NINN, NEPN, NELTYP, NLISEP, NLISEV,
     *                   NELNUM, NELING, QGROUP, QSQR, QPROD, IELV, 
     *                   IINV, IEPA, ITYPEG, IELING, ISTAEV, IELVAR,
C ** Correction 4. 26/02/01: 1 dummy argument removed from SQHESS **
     *                   INLIST, ITABLE, ISTATE, ITYPEE, ISTEP,
     *                   IPTYPE, ISTYPE, INREP,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   ENAMES, GNAMES, ETYPES, ONAMES, INAMES, 
     *                   LNAMES, WEIGHT, KEY, IOUT, INFORM )
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C
C  INDICATOR CARD IS ELEMENT TYPE.
C  -------------------------------
C
  460 CONTINUE
      CALL       SETYPE( NLMAX,  NIMAX, NETMAX, NEPMAX, LENGTH,
     *                   NOVALS, NELN, NINN, NEPN, NELTYP,
     *                   INREP, IELV, IINV, IEPA, INLIST, ITABLE,
C ** Correction 5. 26/02/01: 2 dummy arguments removed from SETYPE **
     *                   FIELD1, FIELD2, FIELD3, FIELD5,
     *                   ENAMES, INAMES, EPNAME, ETYPES, KEY,
     *                   IOUT, INFORM )
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C
C  INDICATOR CARD IS ELEMENT USES.
C  --------------------------------
C
  470 CONTINUE
      CALL       SEUSES( NLMAX, NELMAX, NETMAX, NEVMAX, NLISEV, NLISEP,
     *                   NOVALS, NEPMAX, NEPVMX, LENGTH, NELNUM, NELMNT,
     *                   NMAX, N, ELMNT, IELV, IEPA, ITYPEE, IELVAR,
     *                   INLIST, ITABLE, ISTAEV, ISTEP, DELSET, DETYPE,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   EPVALU, ENAMES, LNAMES, EPNAME, VNAMES,
     *                   KEY, IOUT, INFORM )
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C
C  INDICATOR CARD IS GROUP TYPE.
C  -----------------------------
C
  480 CONTINUE
      CALL       SGTYPE( NGRMAX, NGPMAX, NOVALS, LENGTH, NGRTYP, NGPN,
     *                   SETANA, INLIST, IGPA, ITABLE,
C ** Correction 6. 26/02/01: 2 dummy arguments removed from SGTYPE **
     *                   FIELD1, FIELD2, FIELD3, FIELD5, 
     *                   ANAMES, GTYPES, GPNAME, KEY, IOUT, INFORM )
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C
C  INDICATOR CARD IS GROUP USES.
C  -----------------------------
C
  490 CONTINUE
      CALL       SGUSES( NEGMAX, NGPMAX, NGRMAX, NGMAX, NGPVMX,
     *                   LENGTH, NG, NGRUPE, NLISGP, NOVALS, NELING,
     *                   NDTYPE, STRTGU, GRUPE, IGPA, ITYPEG, IELING,
     *                   INLIST, ITABLE, ISTGP, ISTATE, DGRSET, DGTYPE,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   GPTEMP, GNAMES, GPNAME, WEIGHT,
     *                   KEY, IOUT, INFORM )
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C
C  INDICATOR CARD IS OBJECT BOUND.
C  -------------------------------
C
  500 CONTINUE
      CALL       SOBBND( NOBBND, NOBMAX, NRLNDX, LENGTH,
     *                   INLIST, ITABLE, FBOUND, REALVL,
C ** Correction 7. 26/02/01: 2 dummy arguments removed from SOBBND **
     *                   FIELD1, FIELD2, VALUE4, FIELD5,
     *                   OBNAME, KEY   , SINGLE, IOUT  , INFORM )
      IF ( INFORM .EQ. 0 ) THEN
C ** Correction 11  30/06/04: Assigned goto statements replaced. 1 line replaced
         IF ( DOLOOP ) GO TO 600
         GO TO 100
      END IF
      IF ( INFORM .GT. 0 ) GO TO 800
      GO TO 700
C ** Correction 11  30/06/04: Assigned goto statements replaced. 11 lines added
C
C  BRANCH BACK INTO DO LOOPS AT APPROPRIATE POINT
C
  600 CONTINUE
      IF ( IJUMP .EQ. 1 ) THEN
        GO TO 230
      ELSE IF ( IJUMP .EQ. 2 ) THEN
        GO TO 250
      ELSE
        GO TO 270
      END IF
C
C  INSUFFICIENT SPACE.
C
  700 CONTINUE
      IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 ) INCRSE( - INFORM )
      RETURN
C
C  DATA CARD INCONSISTENCY.
C
  800 CONTINUE
      IF ( IOUT .GT. 0 ) THEN
         IF ( DOLOOP ) THEN
            WRITE( IOUT, 2960 ) LINENO, FIELD1, FIELD2, FIELD3,
     *                          VALUE4, FIELD5, VALUE6
         ELSE
            WRITE( IOUT, 2990 ) LINENO, NULINE
         END IF
      END IF
      GO TO 960
C ** Correction -3. 21/02/00: Code to handle incomplete/missing data added
C
C  MISSING/INCOMPLETE DATA CARD.
C
  810 CONTINUE
      IF ( .NOT. DEFNAM ) THEN
         INFORM = 74
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2740 )
      ELSE
         INFORM = 75
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2750 )
      END IF
      GO TO 960
C
C  SUCCESSFUL RETURN.
C
  900 CONTINUE
      INFORM = 0
      IF ( DEBUG .AND. IOUT .GT. 0 ) THEN
         WRITE( IOUT, 3000 ) ( GNAMES( J ), J = 1, NG )
         WRITE( IOUT, 3010 ) ( VNAMES( J ), J = 1, N )
         WRITE( IOUT, 3020 )
     *       ( ICOORD( J, 1 ), ICOORD( J, 2 ), A( J ), J = 1, NNZA )
         WRITE( IOUT, 3030 )
     *    ( ( I, J, BND( 1, J, I ), BND( 2, J, I ), J = 1, NLVARS ),
     *        I =1, NBND )
         WRITE( IOUT, 3100 )
     *    ( ( I, J, VSTART( J, I ), J = 1, NLVARS ), I = 1, NSTART )
         WRITE( IOUT, 3040 ) ( RSCALE( J ), J = 1, NG )
         WRITE( IOUT, 3050 ) ( CSCALE( J ), J = 1, N )
         IF ( NELTYP .GT. 0 )
     *      WRITE( IOUT, 3060 ) ( ETYPES( I ), IELV( I + 1 ) - IELV( I),
     *                            IINV( I + 1 ) - IINV( I ),
     *                            IEPA( I + 1 ) - IEPA( I ),
     *                            I = 1, NELTYP )
         IF ( NGRTYP .GT. 0 )
     *      WRITE( IOUT, 3110 ) ( GTYPES( I ), ANAMES( I ),
     *                            IGPA( I + 1 ) - IGPA( I ),
     *                            I = 1, NGRTYP )
         WRITE( IOUT, 3070 )
         DO 950 I = 1, NG
            K1    = ISTADG( I )
            K2    = ISTADG( I + 1 ) - 1
            IS    = ITYPEG( I )
            IF ( K1 .LE. K2 ) THEN
               IF ( IS .EQ. 0 ) THEN
                  DO 930 K = K1, K2
                     L     = IELING( K, 1 )
                     WRITE( IOUT, 3080 ) GNAMES( I ),
     *               'TRIVIAL   ', LNAMES( L ),
     *               ETYPES( ITYPEE( L ) ), ( VNAMES( IELVAR( J ) ),
     *               J = ISTAEV( L ), ISTAEV( L + 1 ) - 1 )
  930             CONTINUE
               ELSE
                  DO 940 K = K1, K2
                     L     = IELING( K, 1 )
                     WRITE( IOUT, 3080 ) GNAMES( I ),
     *               GTYPES( IS ), LNAMES( L ),
     *               ETYPES( ITYPEE( L ) ), ( VNAMES( IELVAR( J ) ),
     *               J = ISTAEV( L ), ISTAEV( L + 1 ) - 1 )
  940             CONTINUE
               END IF
            ELSE
               IF ( IS .EQ. 0 ) THEN
                  WRITE( IOUT, 3090 ) GNAMES( I ), 'TRIVIAL   '
               ELSE
                  WRITE( IOUT, 3090 ) GNAMES( I ), GTYPES( IS )
               END IF
            END IF
  950    CONTINUE
      END IF
  960 CONTINUE
      IF ( DEBUG .AND. IOUT .GT. 0 ) THEN
         DO 980 I = 1, NUSEIN
            WRITE( IOUT, 4000 ) I, NAMIIN( I ), INDVAL( I )
  980    CONTINUE
         DO 990 I = 1, NUSERE
            WRITE( IOUT, 4100 ) I, NAMRIN( I ), REALVL( I )
  990    CONTINUE
      END IF
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 1000 FORMAT( A65 )
 1010 FORMAT( A160 )
 2000 FORMAT( ' ** Exit from GPSMPS - insufficient space.',
     *        ' Increase size of ', A6 )
 2010 FORMAT( ' ** Exit from GPSMPS - first card is not NAME ' )
 2020 FORMAT( ' ** Exit from GPSMPS - indicator card not recognised ' )
 2030 FORMAT( ' ** Exit from GPSMPS - index parameter name ', A10,
     *        ' not recognised ' )
 2060 FORMAT( ' ** Exit from GPSMPS - non array defn. within do-loop ' )
 2130 FORMAT( ' ** Exit from GPSMPS - do loop level greater than 3 ' )
 2160 FORMAT( ' ** Exit from GPSMPS - element ', A10, ' variable ',
     *        A10, ' not set ' )
 2210 FORMAT( ' ** Exit from GPSMPS -',
     *        ' groups and variables sections mixed')
 2250 FORMAT( ' ** Exit from GPSMPS - no group-type arg. given ' )
 2280 FORMAT( ' ** Exit from GPSMPS - element ', A10, ' parameter ',
     *        A10, ' not set ' )
 2340 FORMAT( ' ** Exit from GPSMPS - group ', A10, ' parameter ',
     *        A10, ' not set ' )
 2380 FORMAT( ' ** Exit from GPSMPS - do loop not completed ' )
C ** Correction -3. 21/02/00: Code to handle incomplete/missing data added
 2740 FORMAT( ' ** Exit from GPSMPS - data file empty ' )
 2750 FORMAT( ' ** Exit from GPSMPS - data file incomplete.',
     *        ' No ENDATA card ' )
C ** Correction -5b. 07/09/00: Check for non-useful transformations added
 2760 FORMAT( ' ** Exit from GPSMPS - #internal vars >= #elementals' )
 2960 FORMAT( /, ' From within do loop ending on line ', I5,
     *        ', current line is ', /,
     *        2X, A2, 1X, A10, A10, 1P, D12.4, 3X, A10, D12.4 )
 2970 FORMAT( ' Line ', I5, 4X, A160 )
 2980 FORMAT( ' Line ', I5, '.', I2, 1X, A65 )
 2990 FORMAT( ' Line ', I5, 4X, A65 )
 3000 FORMAT( /, ' Row names ', /, ' --------- ', /, 8( 1X, A8 ) )
 3010 FORMAT( /, ' Column names ', /, ' ------------', /, 8( 1X, A8 ) )
 3020 FORMAT( /, 3('  Col   Row    Value  '),
     *        /, 3('  ---   ---    -----  '),
     *        /, ( 3( 2I5, 1P, D12.4 ) ) )
 3030 FORMAT( /, 2(' No. var.  Lower bnd   upper bnd '),
     *        /, 2(' --- ----  ---------   --------- '),
     *        /,  ( 2( I3, I6, 1P, 2E12.4 ) ) )
 3040 FORMAT( /, ' Row scaling ', /, ( 1P, D12.4 ) )
 3050 FORMAT( /, ' Column scaling ', /, ( 1P, D12.4 ) )
C ** Correction 9. 27/02/01: Character debug output format increased **
 3060 FORMAT( /, '    Element type No. el. vars. No. in. vars.',
     *           ' No. parameters ',
     *        /, '    ------------ ------------- ------------- ',
     *           ' -------------- ', /, ( 5X, A10, I12, 2( 2X, I12 ) ) )
 3070 FORMAT( /, ' Group      Gr. type    Element   El. type',
     *           '     Variables ',
     *        /, ' -----      --------    -------   --------',
     *           '     --------- ' )
 3080 FORMAT( 1X, A10, 1X, A10, 2X, A10, A10, 4X, 5A10,
     *        /, ( 48X, 5A10 ) )
 3090 FORMAT( 1X, A10, 1X, A10, 2X, '   -    ' )
 3100 FORMAT( /, 2(' No. var.  Start point '),
     *        /, 2(' --- ----  ----------- '),
     *        /,  ( 2( I3, I6, 1P, D12.4 ) ) )
 3110 FORMAT( /, '    Group type   Argument   No. parameters',
     *        /, '    ----------   --------   -------------- ',
     *        /, ( 5X, 2A10, I14 ) )
 4000 FORMAT( ' Int. par. num. ', I5, ' Name = ', A10, ' Value = ', I12)
 4010 FORMAT( ' Level-', I1, ' Instruction ', I4, ' Starting do-loop ' )
 4020 FORMAT( ' Level-', I1, ' Instruction ', I4, ' Ending do-loop ' )
 4030 FORMAT( ' Level-', I1, ' Instruction ', I4,
     *        ' Incrementing do-loop ' )
 4100 FORMAT( ' Real par. num. ', I5, ' Name = ', A10,' Value = ',
     *          1P, D12.4)
 5000 FORMAT( /, ' Level-', I1, ' loop index ', A10, ' = ', I12 )
 5010 FORMAT( ' Level-', I1, ' instruction ', I3,
     *        ' Index ', A10, ' = ', I12 )
 5020 FORMAT( ' Level-', I1, ' instruction ', I3,
     *        ' Index ', A10, ' = ', 1P, D12.4 )
 5060 FORMAT( ' Level-', I1, ' instruction ', I3, ' Set line ', /,
     *        '    FIELD1 = ', A12, ' FIELD2 = ', A10, ' FIELD3 = ',
     *        A10, /, '    VALUE4 = ', 1P, D12.4, ' FIELD5 = ', A10,
     *        ' VALUE6 = ', 1P, D12.4 )
C
C  END OF GPSMPS.
C
      END
