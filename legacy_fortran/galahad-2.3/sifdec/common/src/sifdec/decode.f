C  THIS VERSION: 27/02/2001 AT 21:30:00 PM.
C     ( Last modified on 15 Mar 2001 at 22:28:00 )
C ** Correction report.
C ** Correction -3. 03/11/00: quadratic elements have no internal representation
C ** Correction -2. 07/09/00: Check for non-useful transformations added
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C ** Correction 0. 20/12/99: Array holding variable types introduced.
C ** Correction 1. 28/10/92: 1 line added **
C ** Correction 2. 26/02/01: 4 dummy arguments removed from SVAR1 **
C ** Correction 3. 26/02/01: 1 dummy argument removed from SBOUND **
C ** Correction 4. 26/02/01: 1 dummy argument removed from SQHESS **
C ** Correction 5. 26/02/01: 2 dummy arguments removed from SETYPE **
C ** Correction 6. 26/02/01: 2 dummy arguments removed from SGTYPE **
C ** Correction 7. 26/02/01: 2 dummy arguments removed from SOBBND **
C ** Correction 8. 26/02/01: unused BIG and J removed from SQHESS **
C ** Correction 9. 27/02/01: check to see if there are quadratic elements
C ** Correction 10. 04/04/02: Number of group parameters for trivial group set
C ** End of Correction report.
      SUBROUTINE SGRP1 ( NG, NGMAX, NOMAX, LENGTH, NOBJ, NOVALS,
     *                   INLIST, IDROWS, ISTATE, ITYPEG, ITABLE,
     *                   RDROWS, RSCALE,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   GNAMES, ONAMES, KEY, IOUT, INFORM )
      INTEGER        IOUT, INFORM, LENGTH
      INTEGER        NG, NGMAX, NOMAX
      INTEGER        NOBJ, NOVALS
      DOUBLE PRECISION VALUE4, VALUE6
      CHARACTER * 2  FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        INLIST( LENGTH ), IDROWS( 2, NGMAX )
      INTEGER        ISTATE( NGMAX ), ITYPEG( NGMAX )
      INTEGER        ITABLE ( LENGTH )
      DOUBLE PRECISION           RDROWS( 2, NGMAX ), RSCALE( NGMAX )
      CHARACTER * 10 GNAMES( NGMAX ), ONAMES( NOMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS GROUPS/ROWS/CONSTRAINTS.
C  ------------------------------------------
C  THE GROUPS SECTION APPEARS BEFORE THE VARIABLES SECTION.
C  --------------------------------------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, IFREE, IFIELD, IS, J
      DOUBLE PRECISION             ZERO
      CHARACTER * 12   FIELD
      EXTERNAL         HASHB , HASHC
      PARAMETER      ( ZERO = 0.0D+0 )
C
C  IGNORE 'MARKER' CARDS.
C
      IF ( FIELD3 .EQ. '''MARKER''  ' ) RETURN
C
C  FIND A PLACE TO INSERT THE NEW GROUP NAME IN THE HASH-TABLE.
C
      CALL HASHB ( LENGTH, 12, FIELD2 // 'GR', KEY, ITABLE, IFREE )
      IF ( IFREE .LE. 0 ) THEN
         IF ( IFREE .EQ. 0 ) THEN
            INFORM = - 1
            RETURN
         END IF
C
C  THE GROUP HAS APPEARED BEFORE AS THE J-TH IN THE LIST.
C
         J = INLIST( - IFREE )
      ELSE
C
C  MARK ANY OBJECTIVE FUNCTION ROWS AS SPECIAL.
C
         IF ( FIELD1 .EQ. 'N ' .OR. FIELD1 .EQ. ' N' .OR.
     *        FIELD1 .EQ. 'DN' .OR. FIELD1 .EQ. 'XN' .OR.
     *        FIELD1 .EQ. 'ZN' ) THEN
            NOBJ = NOBJ + 1
            IF( NOBJ .GT. NOMAX ) THEN
               INFORM = - 5
               RETURN
            END IF
            ONAMES ( NOBJ ) = FIELD2
         END IF
C
C  THE GROUP IS THE NG-TH ENCOUNTERED.
C
         NG = NG + 1
         IF ( NG .GE. NGMAX ) THEN
            INFORM = - 6
            RETURN
         END IF
C
C  RECORD THE POSITION OF THE NEW GROUP IN THE TABLE, RECORD ITS
C  NAME AND INITIALISE ITS TYPE AS TRIVIAL.
C
         J               = NG
         INLIST( IFREE ) = NG
         GNAMES( NG )    = FIELD2
         ITYPEG( NG )    = - 1
         IS              = 0
C
C  RECORD THE STATUS, ISTATE, OF THE GROUP. ISTATE( NG ) = :
C  1 THE GROUP IS AN OBJECTIVE FUNCTION TYPE.
C  2 THE GROUP IS AN EQUALITY FUNCTION TYPE.
C  3 THE GROUP IS A LESS-THAN-OR-EQUAL-TO TYPE.
C  4 THE GROUP IS A GREATER-THAN-OR-EQUAL-TO TYPE.
C  5 THE GROUP IS OF D-TYPE AND AN OBJECTIVE FUNCTION TYPE.
C  6 THE GROUP IS OF D-TYPE AND AN EQUALITY FUNCTION TYPE.
C  7 THE GROUP IS OF D-TYPE AND A LESS-THAN-OR-EQUAL-TO TYPE.
C  8 THE GROUP IS OF D-TYPE AND A GREATER-THAN-OR-EQUAL-TO TYPE.
C
         IF ( FIELD1 .EQ. 'N ' .OR. FIELD1 .EQ. ' N' .OR.
     *        FIELD1 .EQ. 'XN' .OR. FIELD1 .EQ. 'ZN' ) IS = 1
         IF ( FIELD1 .EQ. 'E ' .OR. FIELD1 .EQ. ' E' .OR.
     *        FIELD1 .EQ. 'XE' .OR. FIELD1 .EQ. 'ZE' ) IS = 2
         IF ( FIELD1 .EQ. 'L ' .OR. FIELD1 .EQ. ' L' .OR.
     *        FIELD1 .EQ. 'XL' .OR. FIELD1 .EQ. 'ZL' ) IS = 3
         IF ( FIELD1 .EQ. 'G ' .OR. FIELD1 .EQ. ' G' .OR.
     *        FIELD1 .EQ. 'XG' .OR. FIELD1 .EQ. 'ZG' ) IS = 4
         IF ( FIELD1 .EQ. 'DN' ) IS = 5
         IF ( FIELD1 .EQ. 'DE' ) IS = 6
         IF ( FIELD1 .EQ. 'DL' ) IS = 7
         IF ( FIELD1 .EQ. 'DG' ) IS = 8
         IF ( IS .EQ. 0 ) THEN
            INFORM = 10
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2100 ) FIELD1
            RETURN
         ENDIF
         ISTATE( NG ) = IS
      END IF
C
C  INCLUDE GROUP SCALE FACTORS.
C
      IF ( FIELD3 .EQ. '''SCALE''   ' .OR. FIELD3
     *     .EQ. ' ''SCALE''  ' ) RSCALE( J ) = VALUE4
C
C  MARK 'D'-TYPE GROUPS AND RECORD THEIR MULTIPLICATIVE FACTORS.
C  IDROWS(1, ), IDROWS(2, ) GIVE THE NUMBERS OF THE GROUPS REFERRED
C  TO BY THE NEW D-TYPE GROUP AND RDROWS(1, ) AND RDROWS(2, ) GIVE
C  THE MULTIPLICATIVE FACTORS.
C
      IF ( FIELD1( 1: 1 ) .EQ. 'D' ) THEN
         IF ( ISTATE( J ) .LE. 4 ) THEN
            INFORM = 22
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2220 )
            RETURN
         END IF
         DO 10 I = 1, 2
            IF ( I .GT. NOVALS ) THEN
               IDROWS( I, J ) = 1
               RDROWS( I, J ) = ZERO
            ELSE
C
C  CHECK THAT THE GROUPS REFERRED TO WHEN CONSTRUCTING A D-TYPE
C  GROUP ALREADY EXIST.
C
               IF ( I .EQ. 1 ) THEN
                  FIELD = FIELD3 // 'GR'
               ELSE
                  FIELD = FIELD5 // 'GR'
               END IF
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .GT. 0 ) THEN
                  IDROWS( I, J ) = INLIST( IFIELD )
                  IF ( I .EQ. 1 ) THEN
                     RDROWS( I, J ) = VALUE4
                  ELSE
                     RDROWS( I, J ) = VALUE6
                  END IF
               ELSE
                  INFORM = 4
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2040 ) FIELD( 1: 10 )
                  RETURN
               END IF
            END IF
   10    CONTINUE
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2040 FORMAT( ' ** Exit from GPSMPS - group/row name not recognised:',
     *        ' name is ', A10 )
 2100 FORMAT( ' ** Exit from GPSMPS - field 1 ', A2,
     *        '  not recognised in GROUPS section ' )
 2220 FORMAT( ' ** Exit from GPSMPS -',
     *        ' conflicting field 1 on GROUPS card')
      END
C  THIS VERSION: 28/10/1992 AT 02:36:30 PM.
      SUBROUTINE SGRP2 ( NG, NGMAX, NOMAX, LA, LENGTH,
     *                   NNZA, NOBJ, NOVALS, INLIST, IDROWS,
     *                   ISTATE, ITYPEG, ITABLE, ICOORD,
     *                   A, RDROWS, RSCALE,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   GNAMES, ONAMES, KEY, IOUT, INFORM )
      INTEGER        IOUT, INFORM, LA, LENGTH
      INTEGER        NG, NGMAX, NOMAX
      INTEGER        NNZA, NOBJ, NOVALS
      DOUBLE PRECISION VALUE4, VALUE6
      CHARACTER * 2  FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        INLIST( LENGTH ), IDROWS( 2, NGMAX )
      INTEGER        ISTATE( NGMAX ), ITYPEG( NGMAX )
      INTEGER        ITABLE ( LENGTH )
      INTEGER        ICOORD( LA, 2 )
      DOUBLE PRECISION A( LA ), RDROWS( 2, NGMAX ), RSCALE( NGMAX )
      CHARACTER * 10 GNAMES( NGMAX ), ONAMES( NOMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS GROUPS/ROWS/CONSTRAINTS.
C  ------------------------------------------
C  THE GROUPS SECTION APPEARS AFTER THE VARIABLES SECTION.
C  --------------------------------------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, IFREE, IFIELD, IS, J
      DOUBLE PRECISION ZERO
      CHARACTER * 12   FIELD
      EXTERNAL         HASHB , HASHC
      PARAMETER      ( ZERO = 0.0D+0 )
C
C  IGNORE 'MARKER' CARDS.
C
      IF ( FIELD3 .EQ. '''MARKER''  ' ) RETURN
C
C  FIND A PLACE TO INSERT THE NEW GROUP NAME IN THE HASH-TABLE.
C
      CALL HASHB ( LENGTH, 12, FIELD2 // 'GR', KEY, ITABLE, IFREE )
C
C  THE NAME ALREADY EXISTS. IT IS THE J-TH NAME IN THE LIST.
C
      IF ( IFREE .LE. 0 ) THEN
         IF ( IFREE .EQ. 0 ) THEN
            INFORM = - 1
            RETURN
         END IF
         J = INLIST( - IFREE )
      ELSE
C
C  MARK ANY OBJECTIVE FUNCTION ROWS AS SPECIAL.
C
         IF ( FIELD1 .EQ. 'N ' .OR. FIELD1 .EQ. ' N' .OR.
     *        FIELD1 .EQ. 'DN' .OR. FIELD1 .EQ. 'XN' .OR.
     *        FIELD1 .EQ. 'ZN' ) THEN
            NOBJ = NOBJ + 1
            IF( NOBJ .GT. NOMAX ) THEN
               INFORM = - 5
               RETURN
            END IF
            ONAMES ( NOBJ ) = FIELD2
         END IF
C
C  THE GROUP IS THE NG-TH ENCOUNTERED.
C
         NG = NG + 1
         J  = NG
         IF ( NG .GE. NGMAX ) THEN
            INFORM = - 6
            RETURN
         END IF
C
C  RECORD THE POSITION OF THE NEW GROUP IN THE TABLE, RECORD ITS
C  NAME AND INITIALISE ITS TYPE AS TRIVIAL.
C
         INLIST( IFREE ) = NG
         GNAMES( NG )    = FIELD2
         ITYPEG( NG )    = - 1
         IS              = 0
C
C  RECORD THE STATUS, ISTATE, OF THE GROUP. ISTATE( NG ) = :
C  1 THE GROUP IS AN OBJECTIVE FUNCTION TYPE.
C  2 THE GROUP IS AN EQUALITY FUNCTION TYPE.
C  3 THE GROUP IS A LESS-THAN-OR-EQUAL-TO TYPE.
C  4 THE GROUP IS A GREATER-THAN-OR-EQUAL-TO TYPE.
C  5 THE GROUP IS OF D-TYPE AND AN OBJECTIVE FUNCTION TYPE.
C  6 THE GROUP IS OF D-TYPE AND AN EQUALITY FUNCTION TYPE.
C  7 THE GROUP IS OF D-TYPE AND A LESS-THAN-OR-EQUAL-TO TYPE.
C  8 THE GROUP IS OF D-TYPE AND A GREATER-THAN-OR-EQUAL-TO TYPE.
C
         IF ( FIELD1 .EQ. 'N ' .OR. FIELD1 .EQ. ' N' .OR.
     *        FIELD1 .EQ. 'XN' .OR. FIELD1 .EQ. 'ZN' ) IS = 1
         IF ( FIELD1 .EQ. 'E ' .OR. FIELD1 .EQ. ' E' .OR.
     *        FIELD1 .EQ. 'XE' .OR. FIELD1 .EQ. 'ZE' ) IS = 2
         IF ( FIELD1 .EQ. 'L ' .OR. FIELD1 .EQ. ' L' .OR.
     *        FIELD1 .EQ. 'XL' .OR. FIELD1 .EQ. 'ZL' ) IS = 3
         IF ( FIELD1 .EQ. 'G ' .OR. FIELD1 .EQ. ' G' .OR.
     *        FIELD1 .EQ. 'XG' .OR. FIELD1 .EQ. 'ZG' ) IS = 4
         IF ( FIELD1 .EQ. 'DN' ) IS = 5
         IF ( FIELD1 .EQ. 'DE' ) IS = 6
         IF ( FIELD1 .EQ. 'DL' ) IS = 7
         IF ( FIELD1 .EQ. 'DG' ) IS = 8
         IF ( IS .EQ. 0 ) THEN
            INFORM = 10
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2100 ) FIELD1
            RETURN
         ENDIF
         ISTATE( NG ) = IS
      END IF
C
C  INCLUDE GROUP SCALE FACTORS.
C
      IF ( FIELD3 .EQ. '''SCALE''   ' .OR. FIELD3
     *     .EQ. ' ''SCALE''  ' ) THEN
         RSCALE( J ) = VALUE4
         RETURN
      END IF
C
C  MARK 'D'-TYPE GROUPS AND RECORD THEIR MULTIPLICATIVE FACTORS.
C  IDROWS(1, ), IDROWS(2, ) GIVE THE NUMBERS OF THE GROUPS REFERRED
C  TO BY THE NEW D-TYPE GROUP AND RDROWS(1, ) AND RDROWS(2, ) GIVE
C  THE MULTIPLICATIVE FACTORS.
C
      IF ( FIELD1( 1: 1 ) .EQ. 'D' ) THEN
         IF ( ISTATE( J ) .LE. 4 ) THEN
            INFORM = 22
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2220 )
            RETURN
         END IF
         DO 10 I = 1, 2
            IF ( I .GT. NOVALS ) THEN
               IDROWS( I, J ) = 1
               RDROWS( I, J ) = ZERO
            ELSE
C
C  CHECK THAT THE GROUPS REFERRED TO WHEN CONSTRUCTING A D-TYPE
C  GROUP ALREADY EXIST.
C
               IF ( I .EQ. 1 ) THEN
                  FIELD = FIELD3 // 'GR'
               ELSE
                  FIELD = FIELD5 // 'GR'
               END IF
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .GT. 0 ) THEN
                  IDROWS( I, J ) = INLIST( IFIELD )
                  IF ( I .EQ. 1 ) THEN
                     RDROWS( I, J ) = VALUE4
                  ELSE
                     RDROWS( I, J ) = VALUE6
                  END IF
               ELSE
                  INFORM = 4
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2040 ) FIELD( 1: 10 )
                  RETURN
               END IF
            END IF
   10    CONTINUE
      ELSE
         IF ( NOVALS .GT. 0 ) THEN
C
C  CHECK THAT DATA HAS NOT BEEN SPECIFIED FOR A 'D'-GROUP.
C
            IF ( ISTATE( J ) .GE. 5 ) THEN
               INFORM = 8
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2080 )
               RETURN
            END IF
C
C  ENTRIES FOR THE LINEAR ELEMENT FOR GROUP J ARE TO BE SPECIFIED.
C  FIND THE VARIABLE NUMBERS FOR THE INPUT ENTRIES.
C
            DO 20 I = 1, NOVALS
               IF ( I .EQ. 1 ) THEN
                  FIELD = FIELD3 // 'VA'
               ELSE
                  FIELD = FIELD5 // 'VA'
               END IF
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .GT. 0 ) THEN
C
C  THE NNZA-TH NONZERO HAS BEEN SPECIFIED.
C
                  NNZA = NNZA + 1
                  IF ( NNZA .GT. LA ) THEN
                     INFORM = - 2
                    RETURN
                  END IF
C
C  THE NONZERO IS FOR GROUP ICOORD( ,1) AND VARIABLE ICOORD( ,2).
C  ITS VALUE IS GIVEN IN A( ).
C
                  ICOORD( NNZA, 1 ) = J
                  ICOORD( NNZA, 2 ) = INLIST( IFIELD )
                  IF ( I .EQ. 1 ) THEN
                     A( NNZA ) = VALUE4
                  ELSE
                     A( NNZA ) = VALUE6
                  END IF
               ELSE
C
C  THE VARIABLE NAME IS UNKNOWN.
C
                  INFORM = 5
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2050 ) FIELD( 1 : 10 )
                  RETURN
               END IF
   20       CONTINUE
         END IF
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2040 FORMAT( ' ** Exit from GPSMPS - group/row name not recognised:',
     *        ' name is ', A10 )
 2050 FORMAT( ' ** Exit from GPSMPS - column/var name not recognised:',
     *        ' name is ', A10 )
 2080 FORMAT( ' ** Exit from GPSMPS - ''D'' group/row contains data ' )
 2100 FORMAT( ' ** Exit from GPSMPS - field 1 ', A2,
     *        '  not recognised in GROUPS section ' )
 2220 FORMAT( ' ** Exit from GPSMPS -',
     *        ' conflicting field 1 on GROUPS card')
      END
C  THIS VERSION: 20/12/1999 AT 18:00:00 PM.
      SUBROUTINE SVAR1 ( NMAX, LENGTH, NVAR, COLFIE,
C ** Correction 0. 20/12/99: Array holding variable types introduced.
     *                   ITYPEV, INLIST, ITABLE, CSCALE,
C ** Correction 2. 26/02/01: 4 dummy arguments removed from SVAR1 **
     *                   FIELD2, FIELD3, VALUE4, 
     *                   VNAMES, KEY, INFORM )
      INTEGER        INFORM, LENGTH
      INTEGER        NMAX
      INTEGER        NVAR
      DOUBLE PRECISION VALUE4
      CHARACTER * 2  COLFIE
      CHARACTER * 10 FIELD2, FIELD3
      INTEGER        INLIST( LENGTH ),  ITABLE ( LENGTH )
C ** Correction 0. 20/12/99: Array holding variable types introduced.
      INTEGER        ITYPEV( NMAX )
      DOUBLE PRECISION CSCALE( NMAX )
      CHARACTER * 10 VNAMES( NMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS COLUMNS/VARIABLES, CONSTANTS/RHS/RHS' OR RANGES.
C  ------------------------------------------------------------------
C  THE VARIABLES SECTION BEFORE THE GROUPS SECTION.
C  -------------------------------------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          IFREE, J
      EXTERNAL         HASHB
C
C  IGNORE 'MARKER' CARDS.
C
      IF ( FIELD3 .EQ. '''MARKER''  ' ) RETURN
C
C  FIND A PLACE TO INSERT THE NEW VARIABLE NAME IN THE HASH-TABLE
C  IF IT HAS NOT ALREADY BEEN ENTERED.
C
      CALL HASHB ( LENGTH, 12, FIELD2 // COLFIE, KEY, ITABLE, IFREE )
      IF ( IFREE .LE. 0 ) THEN
C
C  THE VARIABLE NAME ALREADY EXISTS. IT IS THE J-TH NAMED VARIABLE.
C
         IF ( IFREE .EQ. 0 ) THEN
            INFORM = - 1
            RETURN
         END IF
         J = INLIST( - IFREE )
      ELSE
C
C  THE VARIABLE NAME IS NEW. THE VARIABLE IS THE NVAR-TH ENCOUNTERED AND
C  IT OCCURS IN POSITION IFREE IN THE TABLE.
C
         NVAR = NVAR + 1
         IF ( NVAR .GT. NMAX ) THEN
            INFORM = - 7
            RETURN
         END IF
         J               = NVAR
         INLIST( IFREE ) = NVAR
         VNAMES( NVAR )  = FIELD2
      END IF
C
C  INCLUDE COLUMN SCALE FACTORS IF THEY ARE ALLOWED.
C
      IF ( FIELD3 .EQ. '''SCALE''   ' .OR.
     *     FIELD3 .EQ. ' ''SCALE''  ' ) THEN
         CSCALE( J ) = VALUE4
      END IF
C ** Correction 0. 20/12/99: Array holding variable types introduced.
C
C  MARK ZERO-ONE AND INTEGER VARIABLES
C
      IF ( FIELD3 .EQ. '''ZERO-ONE''' ) ITYPEV( J ) = 1
      IF ( FIELD3 .EQ. '''INTEGER'' ' ) ITYPEV( J ) = 2
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
      END
C  THIS VERSION: 20/12/1999 AT 18:00:00 PM.
      SUBROUTINE SVAR2 ( NMAX, NGMAX, LA, LENGTH, NNZA, NVAR, NOVALS,
C ** Correction 0. 20/12/99: Array holding variable types introduced.
     *                   NRLNDX, VARSEC, COLFIE, ICOORD, ISTATE, ITYPEV,
     *                   INLIST, ITABLE, A, CSCALE, REALVL, DFAULT,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   VNAMES, KEY, IOUT, INFORM )
      INTEGER        IOUT, INFORM, LA, LENGTH
      INTEGER        NMAX, NGMAX, NNZA, NVAR, NOVALS, NRLNDX
      DOUBLE PRECISION VALUE4, VALUE6
      LOGICAL        VARSEC
      CHARACTER * 2  FIELD1, COLFIE
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        ICOORD( LA, 2 ), ISTATE( NGMAX )
      INTEGER        INLIST( LENGTH ), ITABLE ( LENGTH )
C ** Correction 0. 20/12/99: Array holding variable types introduced.
      INTEGER        ITYPEV( NMAX )
      DOUBLE PRECISION A( LA ), CSCALE( NMAX )
      DOUBLE PRECISION REALVL( NRLNDX ), DFAULT( NMAX )
      CHARACTER * 10 VNAMES( NMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS COLUMNS/VARIABLES, CONSTANTS/RHS/RHS' OR RANGES.
C  ------------------------------------------------------------------
C  THE VARIABLES SECTION APPEARS AFTER THE GROUPS SECTION.
C  -------------------------------------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, IFREE, IFIELD, J
      DOUBLE PRECISION ZERO, BIGINF
      CHARACTER * 12   FIELD
      EXTERNAL         HASHB , HASHC
      PARAMETER      ( ZERO = 0.0D+0, BIGINF = 1.0D+20 )
C
C  IGNORE 'MARKER' CARDS.
C
      IF ( FIELD3 .EQ. '''MARKER''  ' ) RETURN
C
C  FIND A PLACE TO INSERT THE NEW VARIABLE NAME IN THE HASH-TABLE
C  IF IT HAS NOT ALREADY BEEN ENTERED.
C
      CALL HASHB ( LENGTH, 12, FIELD2 // COLFIE, KEY, ITABLE, IFREE )
C
C  THE VARIABLE NAME ALREADY EXISTS. IT IS THE J-TH NAMED VARIABLE.
C
      IF ( IFREE .LE. 0 ) THEN
         IF ( IFREE .EQ. 0 ) THEN
            INFORM = - 1
            RETURN
         END IF
         J = INLIST( - IFREE )
      ELSE
C
C  THE VARIABLE NAME IS NEW. THE VARIABLE IS THE NVAR-TH ENCOUNTERED AND
C  IT OCCURS IN POSITION IFREE IN THE TABLE.
C
         NVAR = NVAR + 1
         IF ( NVAR .GT. NMAX ) THEN
            INFORM = - 7
            RETURN
         END IF
         J               = NVAR
         INLIST( IFREE ) = NVAR
         VNAMES( NVAR )  = FIELD2
         IF ( COLFIE .EQ. 'CO' ) DFAULT( NVAR ) = ZERO
         IF ( COLFIE .EQ. 'RA' ) DFAULT( NVAR ) = BIGINF
      END IF
C
C  INCLUDE COLUMN SCALE FACTORS IF THEY ARE ALLOWED.
C
      IF ( FIELD3 .EQ. '''SCALE''   ' .OR.
     *     FIELD3 .EQ. ' ''SCALE''  ' ) THEN
         IF ( .NOT. VARSEC ) THEN
            INFORM = 7
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2070 )
         ELSE
            CSCALE( J ) = VALUE4
         END IF
         RETURN
      END IF
C ** Correction 0. 20/12/99: Array holding variable types introduced.
C
C  MARK ZERO-ONE AND INTEGER VARIABLES
C
      IF ( FIELD3 .EQ. '''ZERO-ONE''' ) THEN
         IF ( .NOT. VARSEC ) THEN
            INFORM = 7
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2060 )
         ELSE
            ITYPEV( J ) = 1
         END IF
         RETURN
      END IF
      IF ( FIELD3 .EQ. '''INTEGER'' ' ) THEN
         IF ( .NOT. VARSEC ) THEN
            INFORM = 7
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2060 )
         ELSE
            ITYPEV( J ) = 2
         END IF
         RETURN
      END IF
C
C  A NONTRIVIAL DEFAULT VALUE HAS BEEN SPECIFIED FOR A CONSTANT OR
C  RANGE VECTOR.
C
      IF ( FIELD3 .EQ. '''DEFAULT'' ' .AND. .NOT. VARSEC ) THEN
         IF ( FIELD1 .EQ. 'Z ' ) THEN
            CALL HASHC ( LENGTH, 12,  FIELD5( 1 : 10 ) // 'RI',
     *                   KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 3
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD5( 1 : 10 )
               RETURN
            END IF
            VALUE4 = REALVL( INLIST( IFIELD ) )
         END IF
         DFAULT( NVAR ) = VALUE4
      ELSE
C
C  FIND THE GROUP NUMBERS FOR THE INPUT NONZERO(S).
C
         IF ( NOVALS .GT. 0 ) THEN
            DO 10 I = 1, NOVALS
               IF ( I .EQ. 1 ) THEN
                  FIELD = FIELD3 // 'GR'
               ELSE
                  FIELD = FIELD5 // 'GR'
               END IF
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .GT. 0 ) THEN
C
C  CHECK THAT DATA HAS NOT BEEN SPECIFIED FOR A 'D'-GROUP.
C
                  IF ( ISTATE( INLIST( IFIELD ) ) .GE. 5
     *                 .AND. VARSEC ) THEN
                     INFORM = 8
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2080 )
                     RETURN
                  END IF
C
C  THE NNZA-TH NONZERO HAS BEEN SPECIFIED.
C
                  NNZA = NNZA + 1
                  IF ( NNZA .GT. LA ) THEN
                     INFORM = - 2
                     RETURN
                  END IF
C
C  THE NONZERO IS FOR GROUP ICOORD( ,1) AND VARIABLE ICOORD( ,2).
C  ITS VALUE IS GIVEN IN A( ).
C
                  ICOORD( NNZA, 1 ) = INLIST( IFIELD )
                  ICOORD( NNZA, 2 ) = J
                  IF ( I .EQ. 1 ) THEN
                     A( NNZA ) = VALUE4
                  ELSE
                     A( NNZA ) = VALUE6
                  END IF
               ELSE
C
C  THE GROUP NAME IS UNKNOWN.
C
                  INFORM = 4
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2040 ) FIELD( 1: 10 )
                  RETURN
               END IF
   10       CONTINUE
         END IF
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2030 FORMAT( ' ** Exit from GPSMPS - index parameter name ', A10,
     *        ' not recognised ' )
 2040 FORMAT( ' ** Exit from GPSMPS - group/row name not recognised:',
     *        ' name is ', A10 )
C ** Correction 0. 20/12/99: Array holding variable types introduced.
 2060 FORMAT( ' ** Exit from GPSMPS - type given for RHS or RANGES ' )
 2070 FORMAT( ' ** Exit from GPSMPS - scale given for RHS or RANGES ' )
 2080 FORMAT( ' ** Exit from GPSMPS - ''D'' group/row contains data ' )
      END
C  THIS VERSION: 28/10/1992 AT 02:36:30 PM.
      SUBROUTINE SBOUND( NMAX, NBMAX, LENGTH, NLVARS, NBND, NCOL,
     *                   NRLNDX, DEFAUT, INLIST, ITABLE, BND, REALVL,
C ** Correction 3. 26/02/01: 1 dummy argument removed from SBOUND **
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5,
     *                   BNAMES, BNDFLT, KEY, IOUT, INFORM )
      INTEGER        IOUT, INFORM, LENGTH
      INTEGER        NMAX, NBMAX, NLVARS, NBND, NCOL, NRLNDX
      LOGICAL        DEFAUT
      DOUBLE PRECISION VALUE4
      CHARACTER * 2  FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        INLIST( LENGTH ), ITABLE ( LENGTH )
      DOUBLE PRECISION BND( 2, NMAX, NBMAX ), REALVL( NRLNDX )
      DOUBLE PRECISION BNDFLT( 2, NBMAX )
      CHARACTER * 10 BNAMES( NBMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS BOUNDS.
C  -------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, IFIELD
      DOUBLE PRECISION BIG, ZERO
      EXTERNAL         HASHC
      PARAMETER      ( ZERO = 0.0D+0, BIG = 1.0D+20 )
C
C  THE FIRST PAIR OF BOUND VECTORS ARE TO BE ASSIGNED.
C
      IF ( NBND .EQ. 0 ) THEN
         NBND = 1
         IF ( NBND .GT. NBMAX ) THEN
            INFORM = - 13
            RETURN
         END IF
         BNAMES( NBND ) = FIELD2
         DEFAUT         = .TRUE.
         BNDFLT( 1, NBND )    = ZERO
         BNDFLT( 2, NBND )    = BIG
         DO 10 I              = 1, NLVARS
            BND( 1, I, NBND ) = ZERO
            BND( 2, I, NBND ) = BIG
   10    CONTINUE
      END IF
C
C  A NEW PAIR OF BOUND VECTORS ARE TO BE ASSIGNED.
C
      IF ( FIELD2 .NE. BNAMES( NBND ) ) THEN
         NBND = NBND + 1
         IF ( NBND .GT. NBMAX ) THEN
            INFORM = - 13
            RETURN
         END IF
         BNAMES( NBND ) = FIELD2
         DEFAUT         = .TRUE.
         BNDFLT( 1, NBND )    = ZERO
         BNDFLT( 2, NBND )    = BIG
         DO 20 I              = 1, NLVARS
            BND( 1, I, NBND ) = ZERO
            BND( 2, I, NBND ) = BIG
   20    CONTINUE
      END IF
C
C  ENSURE THAT DEFAULT VALUES ARE ASSIGNED FIRST.
C
      IF ( FIELD3 .EQ. '''DEFAULT'' ' ) THEN
         IF ( .NOT. DEFAUT ) THEN
            INFORM = 20
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2200 )
            RETURN
         END IF
         IF ( FIELD1 .EQ. 'ZL' .OR. FIELD1 .EQ. 'ZU' .OR.
     *        FIELD1 .EQ. 'ZX' ) THEN
            CALL HASHC ( LENGTH, 12,  FIELD5( 1 : 10 ) // 'RI',
     *                   KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 3
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD5( 1 : 10 )
               RETURN
            END IF
            VALUE4 = REALVL( INLIST( IFIELD ) )
         END IF
         IF ( FIELD1 .EQ. 'LO' .OR. FIELD1 .EQ. 'UP' .OR.
     *        FIELD1 .EQ. 'FX' .OR. FIELD1 .EQ. 'FR' .OR.
     *        FIELD1 .EQ. 'MI' .OR. FIELD1 .EQ. 'PL' .OR.
     *        FIELD1 .EQ. 'XL' .OR. FIELD1 .EQ. 'XU' .OR.
     *        FIELD1 .EQ. 'XX' .OR. FIELD1 .EQ. 'XR' .OR.
     *        FIELD1 .EQ. 'XM' .OR. FIELD1 .EQ. 'XP' .OR.
     *        FIELD1 .EQ. 'ZL' .OR. FIELD1 .EQ. 'ZU' .OR.
     *        FIELD1 .EQ. 'ZX' ) THEN
C
C  ASSIGN DEFAULT LOWER BOUNDS FOR VARIABLES.
C
            IF ( ( FIELD1( 2 : 2 ) .EQ. 'L' .AND.
     *             FIELD1 .NE. 'PL' ) .OR.
     *           FIELD1( 2 : 2 ) .EQ. 'X' .OR.
     *           FIELD1( 2 : 2 ) .EQ. 'R' .OR.
     *           FIELD1( 2 : 2 ) .EQ. 'M' .OR.
     *           FIELD1 .EQ. 'LO' .OR. FIELD1 .EQ. 'MI' ) THEN
C
C  A FINITE LOWER BOUND IS SPECIFIED.
C
               IF ( FIELD1( 2 : 2 ) .EQ. 'L' .OR.
     *              FIELD1( 2 : 2 ) .EQ. 'X' .OR.
     *              FIELD1 .EQ. 'LO' ) THEN
                  BNDFLT( 1, NBND )    = VALUE4
                  DO 30 I              = 1, NLVARS
                     BND( 1, I, NBND ) = VALUE4
   30             CONTINUE
C
C  AN INFINITE LOWER BOUND IS SPECIFIED.
C
               ELSE
                  BNDFLT( 1, NBND )    = - BIG
                  DO 40 I              = 1, NLVARS
                     BND( 1, I, NBND ) = - BIG
   40             CONTINUE
                  IF ( FIELD1( 2 : 2 ) .EQ. 'M' .OR.
     *                 FIELD1 .EQ. 'MI' ) THEN
                     BNDFLT( 2, NBND )    = ZERO
                     DO 41 I              = 1, NLVARS
                        BND( 2, I, NBND ) = ZERO
   41                CONTINUE
                  END IF
               END IF
            END IF
C
C  ASSIGN DEFAULT UPPER BOUNDS FOR VARIABLES.
C
            IF ( FIELD1( 2 : 2 ) .EQ. 'U' .OR.
     *           FIELD1( 2 : 2 ) .EQ. 'X' .OR.
     *           FIELD1( 2 : 2 ) .EQ. 'R' .OR.
     *           FIELD1( 2 : 2 ) .EQ. 'P' ) THEN
C
C  A FINITE UPPER BOUND IS SPECIFIED.
C
               IF ( FIELD1( 2 : 2 ) .EQ. 'U' .OR.
     *              FIELD1( 2 : 2 ) .EQ. 'X' .OR.
     *              FIELD1 .EQ. 'UP' ) THEN
                  IF ( ( FIELD1( 2 : 2 ) .EQ. 'U' .OR.
     *                   FIELD1 .EQ. 'UP' )
     *                 .AND. VALUE4 .EQ. ZERO
     *                 .AND. BNDFLT( 1, NBND ) .EQ. ZERO
     *                 .AND. BNDFLT( 2, NBND ) .EQ. BIG ) THEN
                     BNDFLT( 1, NBND )    = - BIG
                     DO 51 I              = 1, NLVARS
                        BND( 1, I, NBND ) = - BIG
   51                CONTINUE
                  END IF
                  BNDFLT( 2, NBND )    = VALUE4
                  DO 50 I              = 1, NLVARS
                     BND( 2, I, NBND ) = VALUE4
   50             CONTINUE
C
C  AN INFINITE UPPER BOUND IS SPECIFIED.
C
               ELSE
                  BNDFLT( 2, NBND )    = BIG
                  DO 60 I              = 1, NLVARS
                     BND( 2, I, NBND ) = BIG
   60             CONTINUE
               END IF
            END IF
C
C  FIELD 1 IS NOT RECOGNISED.
C
         ELSE
            INFORM = 10
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2100 ) FIELD1
            RETURN
         END IF
      ELSE
C
C  AN INDIVIDUAL BOUND IS TO BE ASSIGNED.
C
         DEFAUT = .FALSE.
         IF ( FIELD1 .EQ. 'LO' .OR. FIELD1 .EQ. 'XL' .OR.
     *        FIELD1 .EQ. 'UP' .OR. FIELD1 .EQ. 'XU' .OR.
     *        FIELD1 .EQ. 'FX' .OR. FIELD1 .EQ. 'XX' .OR.
     *        FIELD1 .EQ. 'FR' .OR. FIELD1 .EQ. 'XR' .OR.
     *        FIELD1 .EQ. 'MI' .OR. FIELD1 .EQ. 'XM' .OR.
     *        FIELD1 .EQ. 'PL' .OR. FIELD1 .EQ. 'XP' .OR.
     *        FIELD1 .EQ. 'ZL' .OR. FIELD1 .EQ. 'ZU' .OR.
     *        FIELD1 .EQ. 'ZX' ) THEN
C
C  FIND WHICH VARIABLE IS BEING ASSIGNED.
C
            CALL HASHC ( LENGTH, 12, FIELD3//'VA', KEY, ITABLE, IFIELD )
            IF ( IFIELD .GT. 0 ) THEN
               NCOL = INLIST( IFIELD )
C
C  ASSIGN A LOWER BOUND FOR THIS VARIABLE.
C
               IF ( FIELD1 .EQ. 'LO' .OR. FIELD1 .EQ. 'XL' .OR.
     *              FIELD1 .EQ. 'FX' .OR. FIELD1 .EQ. 'XX' .OR.
     *              FIELD1 .EQ. 'FR' .OR. FIELD1 .EQ. 'XR' .OR.
     *              FIELD1 .EQ. 'MI' .OR. FIELD1 .EQ. 'XM' .OR.
     *              FIELD1 .EQ. 'ZL' .OR. FIELD1 .EQ. 'ZX' ) THEN
C
C  A FINITE LOWER BOUND IS SPECIFIED.
C
                  IF ( FIELD1 .EQ. 'LO' .OR. FIELD1 .EQ. 'XL' .OR.
     *                 FIELD1 .EQ. 'ZL' .OR. FIELD1 .EQ. 'ZX' .OR.
     *                 FIELD1 .EQ. 'FX' .OR. FIELD1 .EQ. 'XX' ) THEN
                     BND( 1, NCOL, NBND ) = VALUE4
C
C  AN INFINITE LOWER BOUND IS SPECIFIED.
C
                  ELSE
                     IF ( ( FIELD1 .EQ. 'MI' .OR. FIELD1 .EQ. 'XM' )
     *                      .AND. BND( 1, NCOL, NBND ) .EQ. ZERO
     *                      .AND. BND( 2, NCOL, NBND ) .EQ. BIG )
     *                  BND( 2, NCOL, NBND ) = ZERO
                     BND( 1, NCOL, NBND ) = - BIG
                  END IF
               END IF
C
C  ASSIGN AN UPPER BOUND FOR THE VARIABLE.
C
               IF ( FIELD1 .EQ. 'UP' .OR. FIELD1 .EQ. 'XU' .OR.
     *              FIELD1 .EQ. 'FX' .OR. FIELD1 .EQ. 'XX' .OR.
     *              FIELD1 .EQ. 'FR' .OR. FIELD1 .EQ. 'XR' .OR.
     *              FIELD1 .EQ. 'PL' .OR. FIELD1 .EQ. 'XP' .OR.
     *              FIELD1 .EQ. 'ZU' .OR. FIELD1 .EQ. 'ZX' ) THEN
C
C  A FINITE UPPER BOUND IS SPECIFIED.
C
                  IF ( FIELD1 .EQ. 'UP' .OR. FIELD1 .EQ. 'XU' .OR.
     *                 FIELD1 .EQ. 'ZU' .OR. FIELD1 .EQ. 'ZX' .OR.
     *                 FIELD1 .EQ. 'FX' .OR. FIELD1 .EQ. 'XX' ) THEN
                     IF ( ( FIELD1 .EQ. 'UP' .OR. FIELD1 .EQ. 'XU' .OR.
     *                      FIELD1 .EQ. 'ZU' ) .AND. VALUE4 .EQ. ZERO
     *                      .AND. BND( 1, NCOL, NBND ) .EQ. ZERO
     *                      .AND. BND( 2, NCOL, NBND ) .EQ. BIG )
     *                    BND( 1, NCOL, NBND ) = - BIG
                     BND( 2, NCOL, NBND ) = VALUE4

C
C  AN INFINITE UPPER BOUND IS SPECIFIED.
C
                  ELSE
                     BND( 2, NCOL, NBND ) = BIG
                  END IF
               END IF
            ELSE
               INFORM = 5
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2050 ) FIELD3
               RETURN
            END IF
C
C  FIELD 1 IS NOT RECOGNISED.
C
         ELSE
            INFORM = 10
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2100 ) FIELD1
            RETURN
         END IF
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2030 FORMAT( ' ** Exit from GPSMPS - index parameter name ', A10,
     *        ' not recognised ' )
 2050 FORMAT( ' ** Exit from GPSMPS - column/var name not recognised:',
     *        ' name is ', A10 )
 2100 FORMAT( ' ** Exit from GPSMPS - field 1 ', A2,
     *        '  not recognised in BOUNDS section ' )
 2200 FORMAT( ' ** Exit from GPSMPS - default specified out of order ' )
      END
C  THIS VERSION: 28/10/1992 AT 02:36:30 PM.
      SUBROUTINE SSTART( NMAX, NGMAX, NSMAX, LENGTH, NLVARS, NG,
     *                   NSTART, NCOL, NRLNDX, DEFAUT, INLIST, ITABLE,
     *                   VSTART, CSTART, REALVL,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   SNAMES, NOVALS, KEY, IOUT, INFORM )
      INTEGER        IOUT, INFORM, LENGTH, NOVALS, NG, NGMAX
      INTEGER        NMAX, NSMAX, NLVARS, NSTART, NCOL, NRLNDX
      LOGICAL        DEFAUT
      DOUBLE PRECISION VALUE4, VALUE6
      CHARACTER * 2  FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        INLIST( LENGTH ), ITABLE ( LENGTH )
      DOUBLE PRECISION VSTART( NMAX, NSMAX ), CSTART( NGMAX, NSMAX )
      DOUBLE PRECISION REALVL( NRLNDX )
      CHARACTER * 10 SNAMES( NSMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS START POINT.
C  ------------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, IFIELD
      DOUBLE PRECISION ZERO
      CHARACTER * 10   FIELD
      EXTERNAL         HASHC
      PARAMETER      ( ZERO = 0.0D+0 )
C
C  THE STARTING VECTOR IS  TO BE ASSIGNED.
C
      IF ( NSTART .EQ. 0 ) THEN
         NSTART = 1
         IF ( NSTART .GT. NSMAX ) THEN
            INFORM = - 8
            RETURN
         END IF
         SNAMES( NSTART ) = FIELD2
         DEFAUT           = .TRUE.
         DO 10 I                = 1, NLVARS
            VSTART( I, NSTART ) = ZERO
   10    CONTINUE
         DO 20 I                = 1, NG
            CSTART( I, NSTART ) = ZERO
   20    CONTINUE
      END IF
C
C  A NEW STARTING VECTOR IS TO BE ASSIGNED.
C
      IF ( FIELD2 .NE. SNAMES( NSTART ) ) THEN
         NSTART = NSTART + 1
         IF ( NSTART .GT. NSMAX ) THEN
            INFORM = - 8
            RETURN
         END IF
         SNAMES( NSTART ) = FIELD2
         DEFAUT           = .TRUE.
C
C  ASSIGN A DEFAULT VALUE OF ZERO TO THE VARIABLES.
C
         DO 30 I                = 1, NLVARS
            VSTART( I, NSTART ) = ZERO
   30    CONTINUE
C
C  ASSIGN A DEFAULT VALUE OF ZERO TO THE LAGRANGE MULTIPLIERS.
C
         DO 40 I                = 1, NG
            CSTART( I, NSTART ) = ZERO
   40    CONTINUE
      END IF
C
C  ENSURE THAT DEFAULT VALUES ARE ASSIGNED FIRST.
C
      IF ( FIELD3 .EQ. '''DEFAULT'' ' ) THEN
         IF ( .NOT. DEFAUT ) THEN
            INFORM = 20
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2200 )
            RETURN
         END IF
         IF ( FIELD1( 1: 1 ) .EQ. 'Z' ) THEN
            CALL HASHC ( LENGTH, 12,  FIELD5( 1 : 10 ) // 'RI',
     *                   KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 3
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD5( 1 : 10 )
               RETURN
            END IF
            VALUE4 = REALVL( INLIST( IFIELD ) )
         END IF
C
C  ASSIGN DEFAULT VALUES TO THE STARTING POINT.
C
         IF ( FIELD1 .EQ. '  ' .OR. FIELD1 .EQ. 'V ' .OR.
     *        FIELD1 .EQ. 'X ' .OR. FIELD1 .EQ. 'Z ' .OR.
     *        FIELD1 .EQ. 'XV' .OR. FIELD1 .EQ. 'ZV' ) THEN
            DO 50 I                = 1, NLVARS
               VSTART( I, NSTART ) = VALUE4
   50       CONTINUE
         END IF
C
C  ASSIGN DEFAULT VALUES TO THE LAGRANGE MULTIPLIERS.
C
         IF ( FIELD1 .EQ. '  ' .OR. FIELD1 .EQ. 'M ' .OR.
     *        FIELD1 .EQ. 'X ' .OR. FIELD1 .EQ. 'Z ' .OR.
     *        FIELD1 .EQ. 'XM' .OR. FIELD1 .EQ. 'ZM' ) THEN
            DO 60 I                = 1, NG
               CSTART( I, NSTART ) = VALUE4
   60       CONTINUE
         END IF
      ELSE
C
C  AN INDIVIDUAL STARTING VALUE IS TO BE ASSIGNED.
C
         IF ( FIELD1 .EQ. 'X ' .OR. FIELD1 .EQ. '  ' .OR.
     *        FIELD1 .EQ. 'Z ' ) THEN
            DEFAUT = .FALSE.
C
C  FIND WHICH VALUE IS, OR VALUES ARE, BEING ASSIGNED.
C
            DO 70 I = 1, NOVALS
               IF ( I .EQ. 1 ) THEN
                  FIELD = FIELD3
               ELSE
                  FIELD = FIELD5
               END IF
C
C  SEE IF THE NAME BELONGS TO A VARIABLE.
C
               CALL HASHC ( LENGTH, 12, FIELD // 'VA',
     *                      KEY, ITABLE, IFIELD )
               IF ( IFIELD .GT. 0 ) THEN
                  NCOL = INLIST( IFIELD )
C
C  ASSIGN THE STARTING VALUE FOR THIS VARIABLE.
C
                  IF ( I .EQ. 1 ) THEN
                     VSTART( NCOL, NSTART ) = VALUE4
                  ELSE
                     VSTART( NCOL, NSTART ) = VALUE6
                  END IF
               ELSE
C
C  SEE IF THE NAME BELONGS TO A GROUP.
C
                  CALL HASHC ( LENGTH, 12, FIELD // 'GR',
     *                         KEY, ITABLE, IFIELD )
                  IF ( IFIELD .GT. 0 ) THEN
                     NCOL = INLIST( IFIELD )
C
C  ASSIGN THE STARTING VALUE FOR THE LAGRANGE MULTIPLIER FOR THIS GROUP.
C
                     IF ( I .EQ. 1 ) THEN
                        CSTART( NCOL, NSTART ) = VALUE4
                     ELSE
                        CSTART( NCOL, NSTART ) = VALUE6
                     END IF
                  ELSE
                     INFORM = 5
                     IF ( I .EQ. 1 ) THEN
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2050 ) FIELD3
                     ELSE
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2050 ) FIELD5
                     END IF
                     RETURN
                  END IF
               END IF
   70       CONTINUE
         ELSE
C
C  AN INDIVIDUAL STARTING VALUE FOR A VARIABLE IS TO BE ASSIGNED.
C
            IF ( FIELD1 .EQ. 'V ' .OR. FIELD1 .EQ. 'XV' .OR.
     *           FIELD1 .EQ. 'ZV' ) THEN
               DEFAUT = .FALSE.
C
C  FIND WHICH VALUE IS, OR VALUES ARE, BEING ASSIGNED.
C
               DO 80 I = 1, NOVALS
                  IF ( I .EQ. 1 ) THEN
                     FIELD = FIELD3
                  ELSE
                     FIELD = FIELD5
                  END IF
C
C  SEE IF THE NAME BELONGS TO A VARIABLE.
C
                  CALL HASHC ( LENGTH, 12, FIELD // 'VA',
     *                         KEY, ITABLE, IFIELD )
                  IF ( IFIELD .GT. 0 ) THEN
                     NCOL = INLIST( IFIELD )
C
C  ASSIGN THE STARTING VALUE FOR THIS VARIABLE.
C
                     IF ( I .EQ. 1 ) THEN
                        VSTART( NCOL, NSTART ) = VALUE4
                     ELSE
                        VSTART( NCOL, NSTART ) = VALUE6
                     END IF
                  ELSE
                     INFORM = 4
                     IF ( I .EQ. 1 ) THEN
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2040 ) FIELD3
                     ELSE
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2040 ) FIELD5
                     END IF
                     RETURN
                  END IF
   80          CONTINUE
            ELSE
C
C  AN INDIVIDUAL STARTING LAGRANGE MULTIPLIER VALUE IS TO BE ASSIGNED.
C
               IF ( FIELD1 .EQ. 'M ' .OR. FIELD1 .EQ. 'XM' .OR.
     *              FIELD1 .EQ. 'ZM' ) THEN
                  DEFAUT = .FALSE.
C
C  FIND WHICH VALUE IS, OR VALUES ARE, BEING ASSIGNED.
C
                  DO 90 I = 1, NOVALS
                     IF ( I .EQ. 1 ) THEN
                        FIELD = FIELD3
                     ELSE
                        FIELD = FIELD5
                     END IF
C
C  SEE IF THE NAME BELONGS TO A GROUP.
C
                     CALL HASHC ( LENGTH, 12, FIELD // 'GR',
     *                            KEY, ITABLE, IFIELD )
                     IF ( IFIELD .GT. 0 ) THEN
                        NCOL = INLIST( IFIELD )
C
C  ASSIGN THE STARTING VALUE FOR THE LAGRANGE MULTIPLIER FOR THIS GROUP.
C
                        IF ( I .EQ. 1 ) THEN
                           CSTART( NCOL, NSTART ) = VALUE4
                        ELSE
                           CSTART( NCOL, NSTART ) = VALUE6
                        END IF
                     ELSE
                        INFORM = 5
                        IF ( I .EQ. 1 ) THEN
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2050 ) FIELD3
                         ELSE
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2050 ) FIELD5
                        END IF
                        RETURN
                     END IF
   90             CONTINUE
C
C  FIELD 1 IS NOT RECOGNISED.
C
               ELSE
                  INFORM = 10
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2100 ) FIELD1
                  RETURN
               END IF
            END IF
         END IF
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2030 FORMAT( ' ** Exit from GPSMPS - index parameter name ', A10,
     *        ' not recognised ' )
 2040 FORMAT( ' ** Exit from GPSMPS - group/row name not recognised:',
     *        ' name is ', A10 )
 2050 FORMAT( ' ** Exit from GPSMPS - column/var name not recognised:',
     *        ' name is ', A10 )
 2100 FORMAT( ' ** Exit from GPSMPS - field 1 ', A2,
     *        '  not recognised in START POINT section ' )
 2200 FORMAT( ' ** Exit from GPSMPS - default specified out of order ' )
      END
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C  THIS VERSION: 20/12/1999 AT 18:00:00 PM.
      SUBROUTINE SQHESS( NEGMAX, NGMAX, NLMAX, NELMAX, NEVMAX, NETMAX, 
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
      INTEGER        IOUT, INFORM, LENGTH, IPTYPE, ISTYPE, NIMAX
      INTEGER        NEGMAX, NGMAX, NOBJ, NELN, NINN, NEPN, NELTYP
      INTEGER        NG, NELING, NOVALS, NGRUPE, NLISEP, NLISEV
      INTEGER        NLMAX, NELMAX, NEVMAX, NETMAX, NOMAX, NELNUM
      DOUBLE PRECISION VALUE4, VALUE6
      LOGICAL        QGROUP, QSQR, QPROD, INREP
      CHARACTER * 2  FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        INLIST( LENGTH ), ITABLE ( LENGTH )
      INTEGER        IELV  ( NLMAX ), IINV  ( NLMAX ), IEPA( NLMAX )
      INTEGER        IELING( NEGMAX, 2 )
      INTEGER        ITYPEG( NGMAX ), ITYPEE( NELMAX ), ISTEP( NELMAX )
      INTEGER        ISTATE( NGMAX )
      INTEGER        IELVAR( NEVMAX ), ISTAEV( NELMAX )
      DOUBLE PRECISION WEIGHT( NEGMAX )
      CHARACTER * 10 GNAMES( NGMAX ), ONAMES( NOMAX ), INAMES( NIMAX )
      CHARACTER * 10 ETYPES( NLMAX ), ENAMES( NETMAX ), LNAMES( NELMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS QUADRATIC.
C  ----------------------------
C
C  NICK GOULD 18/12/1999
C  FOR CGT PRODUCTIONS.
C
C ** Correction 8. 26/02/01: unused BIG and J removed from SQHESS **
      INTEGER          I, IFIELD, K, NEVARS
      INTEGER          NCOL1, NCOL2, NTERMS, IFREE
      DOUBLE PRECISION VALUE
      CHARACTER * 12   FIELD
      CHARACTER * 10   CQGROU, CQSQR, CQPROD
      EXTERNAL         HASHB , HASHC
C ** Correction 8a. 26/02/01: 1 line removed **
      PARAMETER      ( CQGROU = '123456789G' )
      PARAMETER      ( CQSQR  = '123456789S' )
      PARAMETER      ( CQPROD = '123456789P' )
C
C  FIND THE FIRST VARIABLE.
C
      CALL HASHC ( LENGTH, 12, FIELD2 // 'VA', KEY, ITABLE, IFIELD )
      IF ( IFIELD .GT. 0 ) THEN
         NCOL1 = INLIST( IFIELD )
      ELSE
         INFORM = 5
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2050 ) FIELD2
         RETURN
      END IF
C
C  FIND THE SECOND VARIABLE.
C
      IF ( FIELD1 .EQ. 'Z ' ) THEN
         NTERMS = 1
      ELSE
         NTERMS = NOVALS
      END IF


      DO 110 I = 1, NTERMS
         IF ( I .EQ. 1 ) THEN
            FIELD = FIELD3 // 'VA'
            VALUE = VALUE4
         ELSE
            FIELD = FIELD5 // 'VA'
            VALUE = VALUE6
         END IF
         IF ( VALUE .NE. 0.0D+0 ) THEN
            CALL HASHC( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .GT. 0 ) THEN
               NCOL2 = INLIST( IFIELD )
            ELSE
               INFORM = 5
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2050 ) FIELD( 1 : 10 )
               RETURN
            END IF
            IF ( .NOT. QGROUP ) THEN
C
C  THIS IS THE FIRST HESSIAN TERM. MAKE IT A NEW GROUP.
C  FIND A PLACE TO INSERT THE NEW GROUP NAME IN THE HASH-TABLE.
C
               CALL HASHB ( LENGTH, 12, CQGROU // 'GR', KEY, 
     *                      ITABLE, IFREE )
               IF ( IFREE .LE. 0 ) THEN
                  INFORM = - 1
                  RETURN
               ELSE
C
C  MARK THIS AS AN OBJECTIVE FUNCTION GROUP.
C
                  NOBJ = NOBJ + 1
                  IF( NOBJ .GT. NOMAX ) THEN
                     INFORM = - 5
                     RETURN
                  END IF
                  ONAMES ( NOBJ ) =  CQGROU
C
C  THE GROUP IS THE NG-TH ENCOUNTERED.
C
                  NG = NG + 1
                  IF ( NG .GE. NGMAX ) THEN
                     INFORM = - 6
                     RETURN
                  END IF
C
C  RECORD THE POSITION OF THE NEW GROUP IN THE TABLE, RECORD ITS
C  NAME AND INITIALISE ITS TYPE AS TRIVIAL.
C
                  NGRUPE          = NG
C ** Correction 8b. 26/02/01: 1 line removed **
                  INLIST( IFREE ) = NG
                  GNAMES( NG )    = CQGROU
                  ITYPEG( NG )    = 0
                  ISTATE( NG )    = 1
                  QGROUP          = .TRUE.
               END IF
            END IF
            IF ( NCOL1 .EQ. NCOL2 ) THEN
C
C  CHECK IF THIS IS THE FIRST OCCURENCE OF A DIAGONAL TERM
C
               IF ( .NOT. QSQR ) THEN
C
C  CHECK IF THIS IS THE FIRST ELEMENT TYPE.
C
                  IF ( NELTYP .EQ. 0 ) THEN
                     ISTYPE = 1
                     NELN   = 1
                     NINN   = 0
                     NEPN   = 0
                     NELTYP = 1
                  ELSE
                     ISTYPE = 2
                     NELTYP = NELTYP + 1
                     IF ( NELTYP .GT. NLMAX ) THEN
                        INFORM = - 3
                        RETURN
                     END IF
                     NELN = NELN + 1
                     IF ( NELN .GT. NETMAX ) THEN
                        INFORM = - 14
                        RETURN
                     END IF
                  END IF
C
C  INPUT THE NAMES OF THE ELEMENT TYPE.
C
                  CALL HASHB ( LENGTH, 12, CQSQR // 'ET', KEY, 
     *                         ITABLE, IFREE )
                  IF ( IFREE .LE. 0 ) THEN
                     IF ( IFREE .EQ. 0 ) THEN
                        INFORM = - 1
                        RETURN
                     END IF
                     INFORM = 18
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2180 )
                     RETURN
                  END IF
                  INLIST( IFREE )  = NELTYP
                  IELV( NELTYP )   = NELN
                  IINV( NELTYP )   = NINN + 1
                  IEPA( NELTYP )   = NEPN + 1
                  ETYPES( NELTYP ) = CQSQR
C
C  INSERT THE ELEMENT VARIABLE.
C
                  ENAMES( NELN )   = 'X         '
                  NINN = NINN + 1
                  IF ( NINN .GT. NIMAX ) THEN
                     INFORM = - 16
                     RETURN
                  END IF
                  INAMES( NINN ) = ENAMES( NELN )
C ** Correction -3. 03/11/00: quadratic elements have no internal representation
                  INREP          = .FALSE.
                  QSQR           = .TRUE.
               END IF
C
C  THE NEW ELEMENT IS THE NELNUM-TH NONLINEAR ELEMENT.
C
               NELNUM = NELNUM + 1
               IF ( NELNUM .GT. NELMAX ) THEN
                  INFORM = - 9
                  RETURN
               END IF
C
C  INSERT THE NAME INTO THE TABLE.
C
               WRITE( UNIT = FIELD, FMT = "( '%', I9, 'EL' )" ) 
     *                123456789 - NELNUM
               CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
C
C  RECORD THE ELEMENTS POSITION IN THE TABLE ALONG WITH ITS NAME.
C
               INLIST( IFREE )  = NELNUM
               LNAMES( NELNUM ) = FIELD( 1 : 10 )
C
C  DETERMINE THE NUMBER OF THE ELEMENT TYPE, K, THE STARTING
C  ADDRESSES FOR THE PARAMETERS AND VARIABLES FOR THE ELEMENT
C  AND THE NUMBER OF PARAMETERS AND VARIABLES INVOLVED.
C
               K                = ISTYPE
               ITYPEE( NELNUM ) = K
               ISTEP ( NELNUM ) = NLISEP + 1
               ISTAEV( NELNUM ) = NLISEV + 1
               NEVARS           = 1
               IF ( NLISEV + NEVARS .GT. NEVMAX ) THEN
                  INFORM = - 15
                  RETURN
               END IF
               IELVAR( NLISEV + 1 ) = NCOL1
               NLISEV = NLISEV + NEVARS
C
C  ASSIGN THE ELEMENT AS NELING IN GROUP NGRUPE
C
               NELING = NELING + 1
               IF ( NELING .GT. NEGMAX ) THEN
                  INFORM = - 10
                  RETURN
               END IF
C
C  THE ELEMENT IS IELING( ,1) AND THE GROUP IS GIVEN BY IELING( ,2).
C
               IELING( NELING, 1 ) = NELNUM
               IELING( NELING, 2 ) = NGRUPE
C
C  THE ELEMENT IS WEIGHTED BY THE CONSTANT WEIGHT().
C
               WEIGHT( NELING ) = VALUE
            ELSE
               IF ( .NOT. QPROD ) THEN
C
C  CHECK IF THIS IS THE FIRST ELEMENT TYPE.
C
                  IF ( NELTYP .EQ. 0 ) THEN
                     IPTYPE = 1
                     NELN   = 1
                     NINN   = 0
                     NEPN   = 0
                     NELTYP = 1
                  ELSE
                     IPTYPE = 2
                     NELTYP = NELTYP + 1
                     IF ( NELTYP .GT. NLMAX ) THEN
                        INFORM = - 3
                        RETURN
                     END IF
                     NELN = NELN + 1
                     IF ( NELN .GT. NETMAX ) THEN
                        INFORM = - 14
                        RETURN
                     END IF
                  END IF
C
C  INPUT THE NAMES OF THE ELEMENT TYPE.
C
                  CALL HASHB ( LENGTH, 12, CQPROD // 'ET', KEY, 
     *                         ITABLE, IFREE )
                  IF ( IFREE .LE. 0 ) THEN
                     IF ( IFREE .EQ. 0 ) THEN
                        INFORM = - 1
                        RETURN
                     END IF
                     INFORM = 18
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2180 )
                     RETURN
                  END IF
                  INLIST( IFREE )  = NELTYP
                  IELV( NELTYP )   = NELN
                  IINV( NELTYP )   = NINN + 1
                  IEPA( NELTYP )   = NEPN + 1
                  ETYPES( NELTYP ) = CQPROD
C
C  INSERT THE ELEMENT VARIABLE.
C
                  ENAMES( NELN )   = 'X         '
                  NINN = NINN + 1
                  IF ( NINN .GT. NIMAX ) THEN
                     INFORM = - 16
                     RETURN
                  END IF
                  INAMES( NINN ) = ENAMES( NELN )
                  NELN = NELN + 1
                  IF ( NELN .GT. NETMAX ) THEN
                     INFORM = - 14
                     RETURN
                  END IF
                  ENAMES( NELN )   = 'Y         '
                  NINN = NINN + 1
                  IF ( NINN .GT. NIMAX ) THEN
                     INFORM = - 16
                     RETURN
                  END IF
                  INAMES( NINN ) = ENAMES( NELN )
C ** Correction -3. 03/11/00: quadratic elements have no internal representation
                  INREP          = .FALSE.
                  QPROD          = .TRUE.
               END IF
C
C  THE NEW ELEMENT IS THE NELNUM-TH NONLINEAR ELEMENT.
C
               NELNUM = NELNUM + 1
               IF ( NELNUM .GT. NELMAX ) THEN
                  INFORM = - 9
                  RETURN
               END IF
C
C  INSERT THE NAME INTO THE TABLE.
C
               WRITE( UNIT = FIELD, FMT = "( '%', I9, 'EL' )" ) 
     *                123456789 - NELNUM
               CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
C
C  RECORD THE ELEMENTS POSITION IN THE TABLE ALONG WITH ITS NAME.
C
               INLIST( IFREE )  = NELNUM
               LNAMES( NELNUM ) = FIELD( 1 : 10 )
C
C  DETERMINE THE NUMBER OF THE ELEMENT TYPE, K, THE STARTING
C  ADDRESSES FOR THE PARAMETERS AND VARIABLES FOR THE ELEMENT
C  AND THE NUMBER OF PARAMETERS AND VARIABLES INVOLVED.
C
               K                = IPTYPE
               ITYPEE( NELNUM ) = K
               ISTEP ( NELNUM ) = NLISEP + 1
               ISTAEV( NELNUM ) = NLISEV + 1
               NEVARS           = 2
               IF ( NLISEV + NEVARS .GT. NEVMAX ) THEN
                  INFORM = - 15
                  RETURN
               END IF
               IELVAR( NLISEV + 1 ) = NCOL1
               IELVAR( NLISEV + 2 ) = NCOL2
               NLISEV = NLISEV + NEVARS
C
C  ASSIGN THE ELEMENT AS NELING IN GROUP NGRUPE
C
               NELING = NELING + 1
               IF ( NELING .GT. NEGMAX ) THEN
                  INFORM = - 10
                  RETURN
               END IF
C
C  THE ELEMENT IS IELING( ,1) AND THE GROUP IS GIVEN BY IELING( ,2).
C
               IELING( NELING, 1 ) = NELNUM
               IELING( NELING, 2 ) = NGRUPE
C
C  THE ELEMENT IS WEIGHTED BY THE CONSTANT WEIGHT().
C
               WEIGHT( NELING ) = VALUE
            END IF
         END IF 
  110 CONTINUE   
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2050 FORMAT( ' ** Exit from GPSMPS - column/var name not recognised:',
     *        ' name is ', A10 )
 2180 FORMAT( ' ** Exit from GPSMPS - duplicate element-type name ' )
C
C  END OF SUBROUTINE SQHESS.
C
      END
C  THIS VERSION: 28/10/1992 AT 02:36:30 PM.
      SUBROUTINE SETYPE( NLMAX, NIMAX, NETMAX, NEPMAX, LENGTH,
     *                   NOVALS, NELN, NINN, NEPN, NELTYP,
     *                   INREP, IELV, IINV, IEPA, INLIST, ITABLE,
C ** Correction 5. 26/02/01: 2 dummy arguments removed from SETYPE **
     *                   FIELD1, FIELD2, FIELD3, FIELD5, 
     *                   ENAMES, INAMES, EPNAME, ETYPES, KEY,
     *                   IOUT, INFORM )
      INTEGER        IOUT, INFORM, LENGTH
      INTEGER        NLMAX, NIMAX, NETMAX, NEPMAX
      INTEGER        NOVALS, NELN, NINN, NELTYP, NEPN
      LOGICAL        INREP
      CHARACTER * 2  FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        INLIST( LENGTH ), ITABLE ( LENGTH )
      INTEGER        IELV  ( NLMAX ), IINV  ( NLMAX ), IEPA( NLMAX )
      CHARACTER * 10 ETYPES( NLMAX ), INAMES( NIMAX ), ENAMES( NETMAX )
      CHARACTER * 10 EPNAME( NEPMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS ELEMENT TYPE.
C  -------------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, IFREE, K
      EXTERNAL         HASHB
C ** Correction 8b. 27/02/01: check to see if there are quadratic elements
      CHARACTER * 10   CQSQR, CQPROD
      PARAMETER      ( CQSQR  = '123456789S' )
      PARAMETER      ( CQPROD = '123456789P' )
C
C  CHECK IF THIS IS THE FIRST ELEMENT TYPE.
C
      IF ( NELTYP .EQ. 0 ) THEN
         CALL HASHB ( LENGTH, 12, FIELD2 // 'ET', KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               RETURN
            END IF
            INFORM = 18
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2180 )
            RETURN
         END IF
         NELTYP           = 1
         NELN             = 0
         NINN             = 0
         NEPN             = 0
         INREP            = .FALSE.
         INLIST( IFREE )  = NELTYP
         IELV( NELTYP )   = NELN + 1
         IINV( NELTYP )   = NINN + 1
         IEPA( NELTYP )   = NEPN + 1
         ETYPES( NELTYP ) = FIELD2
      END IF
C
C  CHECK IF THE COLUMN IS NEW.
C
      IF ( FIELD2 .NE. ETYPES( NELTYP ) ) THEN
C ** Correction 8. 27/02/01: check to see if there are quadratic elements
         IF ( ETYPES( NELTYP ) .NE. CQSQR .AND.
     *        ETYPES( NELTYP ) .NE. CQPROD ) THEN
C
C  IF THE PREVIOUS ELEMENT HAS NO EXPLICIT INTERNAL REPRESENTATION,
C  USE ITS ELEMENTAL REPRESENTATION.
C
            IF ( .NOT. INREP ) THEN
               DO 10 K = IELV( NELTYP ), NELN
                  NINN = NINN + 1
                  IF ( NINN .GT. NIMAX ) THEN
                     INFORM = - 16
                     RETURN
                  END IF
                  INAMES( NINN ) = ENAMES( K )
   10          CONTINUE
C ** Correction -2a. 07/09/00: Check for non-useful transformations added
            ELSE
              IF ( NINN - IINV( NELTYP ) .GE.
     *             NELN - IELV( NELTYP ) ) THEN
                 INFORM = 76
                 IF ( IOUT .GT. 0 ) WRITE( IOUT, 2760 )
                 RETURN
              END IF
            END IF
         END IF
C
C  RECORD THE NAME AND STARTING ADDRESS OF THE NEW ELEMENT TYPE.
C
         CALL HASHB ( LENGTH, 12, FIELD2 // 'ET', KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               RETURN
            END IF
            INFORM = 18
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2180 )
            RETURN
         END IF
         NELTYP = NELTYP + 1
         INREP  = .FALSE.
         IF ( NELTYP .GT. NLMAX ) THEN
            INFORM = - 3
            RETURN
         END IF
         INLIST( IFREE )  = NELTYP
         IELV( NELTYP )   = NELN + 1
         IINV( NELTYP )   = NINN + 1
         IEPA( NELTYP )   = NEPN + 1
         ETYPES( NELTYP ) = FIELD2
      END IF
C
C  INPUT THE NAME OF AN INTERNAL VARIABLE.
C
      IF ( FIELD1 .EQ. 'IV' ) THEN
         IF ( NOVALS .GT. 0 ) THEN
            INREP   = .TRUE.
            DO 30 I = 1, NOVALS
C
C  CHECK THE NAME HAS NOT ALREADY BEEN USED IN THE CURRENT ELEMENT.
C
                DO 20 K = IINV( NELTYP ), NINN
                   IF ( ( I .EQ. 1 .AND. FIELD3 .EQ. INAMES( K ) )
     *                  .OR. ( I .EQ. 2 .AND. FIELD5 .EQ. INAMES( K ) )
     *                ) THEN
                      INFORM = 12
                      IF ( IOUT .GT. 0 ) WRITE( IOUT, 2120 )
                      RETURN
                   END IF
   20           CONTINUE
C
C  THE NAME IS NEW. RECORD IT IN THE ARRAY INAMES.
C
               NINN = NINN + 1
               IF ( NINN .GT. NIMAX ) THEN
                  INFORM = - 16
                  RETURN
               END IF
               IF ( I .EQ. 1 ) THEN
                  INAMES( NINN ) = FIELD3
               ELSE
                  INAMES( NINN ) = FIELD5
               END IF
   30       CONTINUE
         END IF
      ELSE
C
C  INPUT THE NAME OF AN ELEMENTAL VARIABLE.
C
         IF ( FIELD1 .EQ. 'EV' ) THEN
            IF ( NOVALS .GT. 0 ) THEN
               DO 50 I = 1, NOVALS
C
C  CHECK THE NAME HAS NOT ALREADY BEEN USED IN THE CURRENT ELEMENT.
C
                  DO 40 K = IELV( NELTYP ), NELN
                      IF ( ( I .EQ. 1 .AND. FIELD3 .EQ. ENAMES( K ) )
     *                .OR. ( I .EQ. 2 .AND. FIELD5 .EQ. ENAMES( K ) )
     *                ) THEN
                        INFORM = 11
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2110 )
                        RETURN
                     END IF
   40             CONTINUE
C
C  THE NAME IS NEW. RECORD IT IN THE ARRAY ENAMES.
C
                  NELN = NELN + 1
                  IF ( NELN .GT. NETMAX ) THEN
                     INFORM = - 14
                     RETURN
                  END IF
                  IF ( I .EQ. 1 ) THEN
                     ENAMES( NELN ) = FIELD3
                  ELSE
                     ENAMES( NELN ) = FIELD5
                  END IF
   50          CONTINUE
            END IF
         ELSE
C
C  INPUT THE NAME OF AN ELEMENT PARAMETER.
C
            IF ( FIELD1 .EQ. 'EP' ) THEN
               IF ( NOVALS .GT. 0 ) THEN
                  DO 70 I = 1, NOVALS
C
C  CHECK THE NAME HAS NOT ALREADY BEEN USED IN THE CURRENT ELEMENT.
C
                     DO 60 K = IEPA( NELTYP ), NEPN
                         IF ( ( I .EQ. 1 .AND. FIELD3 .EQ. EPNAME( K ) )
     *                   .OR. ( I .EQ. 2 .AND. FIELD5 .EQ. EPNAME( K ) )
     *                   ) THEN
                           INFORM = 23
                           IF ( IOUT .GT. 0 ) WRITE( IOUT, 2230 )
                           RETURN
                        END IF
   60                CONTINUE
C
C  THE NAME IS NEW. RECORD IT IN THE ARRAY EPNAME.
C
                     NEPN = NEPN + 1
                     IF ( NEPN .GT. NEPMAX ) THEN
                        INFORM = - 19
                        RETURN
                     END IF
                     IF ( I .EQ. 1 ) THEN
                        EPNAME( NEPN ) = FIELD3
                     ELSE
                        EPNAME( NEPN ) = FIELD5
                     END IF
   70             CONTINUE
               END IF
C
C  FIELD1 NOT RECOGNISED.
C
            ELSE
               INFORM = 10
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2100 ) FIELD1
               RETURN
            END IF
         END IF
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2100 FORMAT( ' ** Exit from GPSMPS - field 1 ', A2,
     *        '  not recognised in ELEMENT TYPE section ' )
 2110 FORMAT( ' ** Exit from GPSMPS - duplicate elemental var. name ' )
 2120 FORMAT( ' ** Exit from GPSMPS - duplicate internal var. name ' )
 2180 FORMAT( ' ** Exit from GPSMPS - duplicate element-type name ' )
 2230 FORMAT( ' ** Exit from GPSMPS - duplicate elemental param. name ')
C ** Correction -2b. 07/09/00: Check for non-useful transformations added
 2760 FORMAT( ' ** Exit from GPSMPS - #internal vars >= #elementals' )
      END
C  THIS VERSION: 28/10/1992 AT 02:36:30 PM.
      SUBROUTINE SEUSES( NLMAX, NELMAX, NETMAX, NEVMAX, NLISEV, NLISEP,
     *                   NOVALS, NEPMAX, NEPVMX, LENGTH, NELNUM, NELMNT,
     *                   NMAX, N, ELMNT, IELV, IEPA, ITYPEE, IELVAR,
     *                   INLIST, ITABLE, ISTAEV, ISTEP, DELSET, DETYPE,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   EPVALU, ENAMES, LNAMES, EPNAME, VNAMES,
     *                   KEY, IOUT, INFORM )
      INTEGER        IOUT, INFORM, LENGTH
      INTEGER        NMAX, NLMAX, NELMAX, NETMAX, NEVMAX
      INTEGER        NELNUM, NLISEP, NLISEV
      INTEGER        NOVALS, NEPVMX, NEPMAX, NELMNT, N
      DOUBLE PRECISION VALUE4, VALUE6
      LOGICAL        DELSET
      CHARACTER * 2  FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5, ELMNT, DETYPE
      INTEGER        INLIST( LENGTH ), ITABLE ( LENGTH )
      INTEGER        IELV( NLMAX ), IEPA( NLMAX )
      INTEGER        IELVAR( NEVMAX ), ITYPEE( NELMAX )
      INTEGER        ISTAEV( NELMAX ), ISTEP( NELMAX )
      DOUBLE PRECISION EPVALU( NEPVMX )
      CHARACTER * 10 ENAMES( NETMAX ), LNAMES( NELMAX )
      CHARACTER * 10 EPNAME( NEPMAX ), VNAMES( NMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS ELEMENT USES.
C  -------------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, IFREE, IFIELD, IP, IS, J, K, MLISEP, MLISEV
      INTEGER          NEPARS, NEVARS
      DOUBLE PRECISION BIG
      CHARACTER * 12 FIELD
      EXTERNAL         HASHB , HASHC
      PARAMETER      ( BIG = 1.0D+20 )

      INTEGER LENFIELD
C
C  THE CURRENT CARD DEFINES A DEFAULT TYPE.
C
      IF ( FIELD2 .EQ. '''DEFAULT'' ' ) THEN
         IF ( DELSET ) THEN
            INFORM = 26
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2260 )
         END IF
         DELSET = .TRUE.
         DETYPE = FIELD3
      ELSE
C
C  IF THE ELEMENT NAMED IN FIELD2 IS NOT THAT OF THE PREVIOUS CARD,
C  DETERMINE THE CHARACTERISTICS OF THE ELEMENT.
C
         IF ( ELMNT .NE. FIELD2 ) THEN
C
C  LOOK THE NAME UP IN THE DICTIONARY TO SEE IF IT ALREADY EXISTS.
C
            CALL HASHB ( LENGTH, 12, FIELD2 // 'EL', KEY,
     *                   ITABLE, IFREE )
C
C  IF THE ELEMENT NAME IS RECOGNISED, RECOVER ITS CHARACTERISTICS.
C
            IF ( IFREE .LE. 0 ) THEN
               IF ( IFREE .EQ. 0 ) THEN
                  INFORM = - 1
                  RETURN
               END IF
               NELMNT = INLIST( - IFREE )
            ELSE
C
C  THE NEW ELEMENT IS THE NELNUM-TH NONLINEAR ELEMENT.
C
               NELNUM = NELNUM + 1
               IF ( NELNUM .GT. NELMAX ) THEN
                  INFORM = - 9
                  RETURN
               END IF
C
C  RECORD THE ELEMENTS POSITION IN THE TABLE ALONG WITH ITS NAME.
C
               INLIST( IFREE )  = NELNUM
               LNAMES( NELNUM ) = FIELD2
               NELMNT           = NELNUM
            END IF
C
C  RECORD THE NONLINEAR ELEMENT'S NAME.
C
            ELMNT = FIELD2
C
C  IF THE ELEMENT HAS NOT YET BEEN ALLOCATED A TYPE, SET IT.
C
            IF ( FIELD1 .EQ. 'T ' .OR. FIELD1 .EQ. 'XT'
     *           .OR. IFREE .GT. 0 ) THEN
C
C  RECORD THE ELEMENT TYPE.
C
               IF ( FIELD1 .EQ. 'T ' .OR. FIELD1 .EQ. 'XT' ) THEN
                  FIELD = FIELD3 // 'ET'
               ELSE
                  IF ( DELSET ) THEN
                     FIELD = DETYPE // 'ET'
                  ELSE
C
C  THE ELEMENT NAME IS NEW. CHECK THAT IF A DEFAULT ELEMENT TYPE
C  IS REQUIRED, A DEFAULT HAS BEEN SET.
C
                     INFORM = 41
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2410 )
                     RETURN
                  END IF
               END IF
C
C  IF THE GROUP IS NON-TRIVIAL, DETERMINE ITS CHARACTERISTICS.
C
               LENFIELD = LEN( FIELD )
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .LE. 0 ) THEN
                  INFORM = 9
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2090 ) FIELD3
                  RETURN
               END IF
C
C  DETERMINE THE NUMBER OF THE ELEMENT TYPE, K, THE STARTING
C  ADDRESSES FOR THE PARAMETERS AND VARIABLES FOR THE ELEMENT
C  AND THE NUMBER OF PARAMETERS AND VARIABLES INVOLVED.
C
               K                = INLIST( IFIELD )
               ITYPEE( NELNUM ) = K
               ISTEP ( NELNUM ) = NLISEP + 1
               ISTAEV( NELNUM ) = NLISEV + 1
               NEPARS           = IEPA( K + 1 ) - IEPA( K )
               NEVARS           = IELV( K + 1 ) - IELV( K )
               IF ( NLISEV + NEVARS .GT. NEVMAX ) THEN
                  INFORM = - 15
                  RETURN
               END IF
               IF ( NLISEP + NEPARS .GT. NEPVMX ) THEN
                  INFORM = - 17
                  RETURN
               END IF
C
C  INITIALIZE THE SET OF PROBLEM VARIABLES.
C
               DO 10 I                 = 1, NEVARS
                  IELVAR( NLISEV + I ) = 0
   10          CONTINUE
C
C  INITIALIZE THE SET OF PARAMETER VALUES.
C
               DO 20 I                 = 1, NEPARS
                  EPVALU( NLISEP + I ) = BIG
   20          CONTINUE
C
C  FIND THE STARTING ADDRESSES FOR THE LISTS OF THE ELEMENTS
C  PARAMETERS AND VARIABLES.
C
               NLISEP = NLISEP + NEPARS
               NLISEV = NLISEV + NEVARS
               IF ( FIELD1 .EQ. 'T ' .OR. FIELD1 .EQ. 'XT' ) RETURN
            END IF
         END IF
C
C  CHECK THAT THE CARDS ARE IN THE CORRECT ORDER.
C
         IF ( FIELD1 .EQ. 'T ' .OR. FIELD1 .EQ. 'XT' ) THEN
            INFORM = 27
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2270 )
            RETURN
         END IF
C
C  DETERMINE THE NUMBER OF THE ELEMENT TYPE, K, THE STARTING
C  ADDRESSES FOR THE PARAMETERS AND VARIABLES FOR THE ELEMENT
C  AND THE NUMBER OF PARAMETERS AND VARIABLES INVOLVED.
C
         K      = ITYPEE( NELMNT )
         NEPARS = IEPA( K + 1 ) - IEPA( K )
         NEVARS = IELV( K + 1 ) - IELV( K )
         MLISEP = ISTEP( NELMNT ) - 1
         MLISEV = ISTAEV( NELMNT ) - 1
         IP     = IEPA( K ) - 1
         IS     = IELV( K ) - 1
C
C  THE CARD CONTAINS NAMES OF ELEMENTAL VARIABLES.
C
         IF ( FIELD1 .EQ. 'V ' .OR. FIELD1 .EQ. 'ZV' ) THEN
C
C  THE ELEMENTAL VARIABLE IS DEFINED IN FIELD3.
C
            DO 110 I = 1, NEVARS
               IF ( FIELD3 .EQ. ENAMES( IS + I ) ) GO TO 120
  110       CONTINUE
C
C  THE ELEMENTAL VARIABLE NAME IS NOT RECOGNISED.
C
            INFORM = 15
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2150 )
            RETURN
  120       CONTINUE
C
C  CHECK THAT THE VARIABLE HAS NOT ALREADY BEEN SET.
C
            IF ( IELVAR( MLISEV + I ) .NE. 0 ) THEN
               INFORM = 30
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2300 )
               RETURN
            END IF
C
C  SEARCH THE TABLE FOR THE NAME OF THE INPUT VARIABLE.
C
            CALL HASHB ( LENGTH, 12, FIELD5//'VA', KEY, ITABLE, IFREE )
            IF ( IFREE .LE. 0 ) THEN
               IF ( IFREE .EQ. 0 ) THEN
                  INFORM = - 1
                  RETURN
               END IF
C
C  THE VARIABLE HAS APPEARED BEFORE. STORE ITS NUMBER.
C
               IELVAR( MLISEV + I ) = INLIST( - IFREE )
            ELSE
C
C  THE VARIABLE IS COMPLETELY NEW (AND THUS NONLINEAR).
C  IT WILL BE RECORDED AS VARIABLE N.
C
               N = N + 1
               IF ( N .GT. NMAX ) THEN
                  INFORM = - 7
                  RETURN
               END IF
C
C  RECORD THE POSITION OF THE NEW GROUP IN THE TABLE, RECORD ITS
C  NAME, INITIALISE ITS TYPE AS TRIVIAL AND RECORD ITS STATUS AS
C  AN EQUALITY GROUP.
C
               INLIST( IFREE )      = N
               VNAMES( N )          = FIELD5
               IELVAR( MLISEV + I ) = N
            END IF
         ELSE
C
C  THE CARD CONTAINS NAMES AND VALUES OF ELEMENTAL PARAMETERS.
C
            IF ( FIELD1 .EQ. 'P ' .OR. FIELD1 .EQ. 'XP' .OR.
     *           FIELD1 .EQ. 'ZP' ) THEN
               IF ( NOVALS .GT. 0 ) THEN
                  DO 230 J = 1, NOVALS
C
C  CHECK THE NAME HAS NOT ALREADY BEEN USED IN THE CURRENT ELEMENT.
C  THE PARAMETER NAME OCCURS IN FIELD3 OR FIELD5.
C
                     DO 210 I = 1, NEPARS
                        IF (
     *                  ( J .EQ. 1 .AND. FIELD3 .EQ. EPNAME( IP + I ) )
     *             .OR. ( J .EQ. 2 .AND. FIELD5 .EQ. EPNAME( IP + I ) )
     *                       ) GO TO 220
  210                CONTINUE
C
C  THE ELEMENTAL PARAMETER NAME IS NOT RECOGNISED.
C
                     INFORM = 28
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2280 )
     *                  LNAMES( NELNUM ), EPNAME( IP + I )
                     RETURN
C
C  THE ELEMENTAL PARAMETER NAME IS THE I-TH PARAMETER IN THE LIST.
C
  220                CONTINUE
C
C  CHECK THAT THE VALUE HAS NOT ALREADY BEEN SET.
C
                     IF ( EPVALU( MLISEP + I ) .LT. BIG ) THEN
                        INFORM = 29
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2290 )
                        RETURN
                     END IF
C
C  READ THE ASSOCIATED VALUE FROM FIELD4 OR FIELD 6.
C
                     IF ( J .EQ. 1 ) THEN
                        EPVALU( MLISEP + I ) = VALUE4
                     ELSE
                        EPVALU( MLISEP + I ) = VALUE6
                     END IF
  230             CONTINUE
               END IF
C
C  FIELD1 NOT RECOGNISED.
C
            ELSE
               INFORM = 10
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2100 ) FIELD1
               RETURN
            END IF
         END IF
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2090 FORMAT( ' ** Exit from GPSMPS - element type not recognised:',
     *        ' name is ', A10 )
 2100 FORMAT( ' ** Exit from GPSMPS - field 1 ', A2,
     *        '  not recognised in ELEMENT USES section ' )
 2150 FORMAT( ' ** Exit from GPSMPS - element variable unrecognised ' )
 2260 FORMAT( ' ** Exit from GPSMPS - duplicate default element type ' )
 2270 FORMAT( ' ** Exit from GPSMPS - type for element already set ' )
 2280 FORMAT( ' ** Exit from GPSMPS - element ', A10, ' parameter ',
     *        A10, ' unrecognised ' )
 2290 FORMAT( ' ** Exit from GPSMPS - element parameter already set ' )
 2300 FORMAT( ' ** Exit from GPSMPS - element variable already set ' )
 2410 FORMAT( ' ** Exit from GPSMPS - element type unrecognised ' )
      END
C  THIS VERSION: 28/10/1992 AT 02:36:30 PM.
      SUBROUTINE SGTYPE( NGRMAX, NGPMAX, NOVALS, LENGTH, NGRTYP, NGPN,
     *                   SETANA, INLIST, IGPA, ITABLE,
C ** Correction 6. 26/02/01: 2 dummy arguments removed from SGTYPE **
     *                   FIELD1, FIELD2, FIELD3, FIELD5, 
     *                   ANAMES, GTYPES, GPNAME, KEY, IOUT, INFORM )
      INTEGER        IOUT, INFORM, NOVALS, LENGTH
      INTEGER        NGRMAX, NGPMAX
      INTEGER        NGRTYP, NGPN
      LOGICAL        SETANA
C ** Correction 6. 26/02/01: 2 dummy arguments removed from SGTYPE **
      CHARACTER * 2  FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        ITABLE ( LENGTH ), IGPA( NGRMAX )
      INTEGER        INLIST( LENGTH )
      CHARACTER * 10 ANAMES( NGRMAX ), GTYPES( NGRMAX ), GPNAME( NGPMAX)
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS GROUP TYPE.
C  -----------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, IFREE, K
      EXTERNAL         HASHB
C
C  CHECK IF THIS IS THE FIRST GROUP TYPE.
C
      IF ( NGRTYP .EQ. 0 ) THEN
         CALL HASHB ( LENGTH, 12, FIELD2 // 'GT', KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               RETURN
            END IF
            INFORM = 17
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2170 )
            RETURN
         END IF
         NGRTYP           = 1
         NGPN             = 0
         SETANA           = .FALSE.
         INLIST( IFREE )  = NGRTYP
         IGPA( NGRTYP )   = NGPN + 1
         GTYPES( NGRTYP ) = FIELD2
      END IF
C
C  CHECK IF THE GROUP-TYPE IS NEW.
C
      IF ( FIELD2 .NE. GTYPES( NGRTYP ) ) THEN
C
C  CHECK THAT THE ARGUMENT FOR THE PREVIOUS GROUP-TYPE HAS BEEN SET.
C
         IF ( .NOT. SETANA ) THEN
            INFORM = 25
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2250 )
            RETURN
         END IF
C
C  RECORD THE NAME AND STARTING ADDRESS OF THE NEW GROUP TYPE.
C
         CALL HASHB ( LENGTH, 12, FIELD2 // 'GT', KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               RETURN
            END IF
            INFORM = 17
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2170 )
            RETURN
         END IF
         NGRTYP = NGRTYP + 1
         SETANA = .FALSE.
         IF ( NGRTYP .GT. NGRMAX ) THEN
            INFORM = - 4
            RETURN
         END IF
         INLIST( IFREE )  = NGRTYP
         IGPA( NGRTYP )   = NGPN + 1
         GTYPES( NGRTYP ) = FIELD2
      END IF
C
C  INPUT THE NAME OF THE GROUP-TYPE ARGUMENT.
C
      IF ( FIELD1 .EQ. 'GV' ) THEN
         SETANA           = .TRUE.
         ANAMES( NGRTYP ) = FIELD3
      ELSE
C
C  INPUT THE NAME OF AN GROUP PARAMETER.
C
         IF ( FIELD1 .EQ. 'GP' ) THEN
            IF ( NOVALS .GT. 0 ) THEN
               DO 20 I = 1, NOVALS
C
C  CHECK THE NAME HAS NOT ALREADY BEEN USED IN THE CURRENT GROUP.
C
                  DO 10 K = IGPA( NGRTYP ), NGPN
                      IF ( ( I .EQ. 1 .AND. FIELD3 .EQ. GPNAME( K ) )
     *                .OR. ( I .EQ. 2 .AND. FIELD5 .EQ. GPNAME( K ) )
     *                ) THEN
                        INFORM = 24
                        IF ( IOUT .GT. 0 ) WRITE( IOUT, 2240 )
                        RETURN
                     END IF
   10             CONTINUE
C
C  THE NAME IS NEW. RECORD IT IN THE ARRAY GPNAME.
C
                  NGPN = NGPN + 1
                  IF ( NGPN .GT. NGPMAX ) THEN
                     INFORM = - 20
                     RETURN
                  END IF
                  IF ( I .EQ. 1 ) THEN
                     GPNAME( NGPN ) = FIELD3
                  ELSE
                     GPNAME( NGPN ) = FIELD5
                  END IF
   20          CONTINUE
            END IF
C
C  FIELD1 NOT RECOGNISED.
C
         ELSE
            INFORM = 10
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2100 ) FIELD1
            RETURN
         END IF
      END IF
C
C  NON-EXECUTABLE STATEMENTS.
C
 2170 FORMAT( ' ** Exit from GPSMPS - duplicate group-type name ' )
 2100 FORMAT( ' ** Exit from GPSMPS - field 1 ', A2,
     *        '  not recognised in GROUP TYPE section ' )
 2240 FORMAT( ' ** Exit from GPSMPS - duplicate group param. name ' )
 2250 FORMAT( ' ** Exit from GPSMPS - no group-type arg. given ' )
      END
C  THIS VERSION: 20/12/1999 AT 18:00:00 PM.
      SUBROUTINE SGUSES( NEGMAX, NGPMAX, NGRMAX, NGMAX, NGPVMX,
     *                   LENGTH, NG, NGRUPE, NLISGP, NOVALS, NELING,
     *                   NDTYPE, STRTGU, GRUPE, IGPA, ITYPEG, IELING,
     *                   INLIST, ITABLE, ISTGP, ISTATE, DGRSET, DGTYPE,
     *                   FIELD1, FIELD2, FIELD3, VALUE4, FIELD5, VALUE6,
     *                   GPTEMP, GNAMES, GPNAME, WEIGHT,
     *                   KEY, IOUT, INFORM )
      INTEGER        IOUT, INFORM, LENGTH
      INTEGER        NEGMAX, NGPMAX, NGRMAX, NGMAX, NGPVMX
      INTEGER        NG, NLISGP, NELING, NOVALS, NGRUPE, NDTYPE
      DOUBLE PRECISION VALUE4, VALUE6
      LOGICAL        DGRSET, STRTGU
      CHARACTER * 2  FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5, GRUPE, DGTYPE
      INTEGER        INLIST( LENGTH ), ITABLE ( LENGTH )
      INTEGER        IGPA( NGRMAX ), IELING( NEGMAX, 2 )
      INTEGER        ITYPEG( NGMAX )
      INTEGER        ISTGP( NGMAX ), ISTATE( NGMAX )
      DOUBLE PRECISION GPTEMP( NGPVMX ), WEIGHT( NEGMAX )
      CHARACTER * 10 GPNAME( NGPMAX ), GNAMES( NGMAX )
      CHARACTER * 12 KEY( LENGTH )
C
C  INDICATOR CARD IS GROUP USES.
C  -----------------------------
C
C  NICK GOULD 01/08/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, IFIELD, IP, J, K, MLISGP, NGPARS
      DOUBLE PRECISION BIG
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      CHARACTER * 12   FIELD
      CHARACTER * 10   CQGROU
      EXTERNAL         HASHB , HASHC
      INTRINSIC        ABS
      PARAMETER      ( BIG = 1.0D+20 )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      PARAMETER      ( CQGROU = '123456789G' )
C
C  THE CURRENT CARD DEFINES A DEFAULT TYPE.
C
      IF ( FIELD2 .EQ. '''DEFAULT'' ' ) THEN
         IF ( DGRSET ) THEN
            INFORM = 42
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2420 )
            RETURN
         END IF
         DGRSET = .TRUE.
         DGTYPE = FIELD3
C
C  FIND THE NUMBER ALLOCATED TO THE GROUP TYPE.
C
         CALL HASHC ( LENGTH, 12, DGTYPE // 'GT', KEY, ITABLE, IFIELD )
         IF ( IFIELD .LE. 0 ) THEN
            INFORM = 19
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2190 ) FIELD3
            RETURN
         END IF
C
C  RESET THE DEFAULTS FOR EACH OF THE GROUPS ALLOCATED IN PREVIOUS
C  SECTIONS.
C
         NDTYPE = INLIST( IFIELD )
         DO 10 I = 1, NG
            IF ( ITYPEG( I ) .EQ. - 1 ) ITYPEG( I ) = - NDTYPE - 1
   10    CONTINUE
         RETURN
      END IF
C
C  IF THE GROUP NAMED IN FIELD2 IS NOT THAT OF THE PREVIOUS CARD,
C  DETERMINE THE CHARACTERISTICS OF THE GROUP.
C
      IF ( .NOT. STRTGU .OR. GRUPE .NE. FIELD2 ) THEN
         STRTGU = .TRUE.
C
C  LOOK THE NAME UP IN THE DICTIONARY TO SEE IF IT EXISTS.
C
         CALL HASHC ( LENGTH, 12, FIELD2 // 'GR', KEY,
     *                ITABLE, IFIELD )
         IF ( IFIELD .GT. 0 ) THEN
            NGRUPE = INLIST( IFIELD )
            ISTATE( NGRUPE ) = - ABS( ISTATE( NGRUPE ) )
         ELSE
C
C  THE GROUP NAME IS UNKNOWN.
C
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2040 ) FIELD2
            INFORM = 4
            RETURN
         END IF
C
C  RECORD THE GROUP'S NAME.
C
         GRUPE = FIELD2
C
C  RECORD THE GROUP TYPE.
C
         IF ( FIELD1 .EQ. 'T ' .OR. FIELD1 .EQ. 'XT' ) THEN
            FIELD = FIELD3 // 'GT'
C
C  IF THE GROUP IS NON-TRIVIAL, DETERMINE ITS CHARACTERISTICS.
C
            CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 19
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2190 ) FIELD3
               RETURN
            END IF
C
C  DETERMINE THE NUMBER OF THE GROUP TYPE, K, THE STARTING
C  ADDRESSES FOR THE PARAMETERS FOR THE GROUP AND THE
C  NUMBER OF PARAMETERS INVOLVED.
C
            K                = INLIST( IFIELD )
C ** Correction 10. 1 line replaced by 5
            IF ( K .EQ. 0 ) THEN
              NGPARS = 0
            ELSE
              NGPARS = IGPA( K + 1 ) - IGPA( K )
            END IF
            ITYPEG( NGRUPE ) = K
            ISTGP ( NGRUPE ) = NLISGP + 1
            IF ( NLISGP + NGPARS .GT. NGPVMX ) THEN
               INFORM = - 18
               RETURN
            END IF
C
C  INITIALIZE THE SET OF PARAMETER VALUES.
C
            DO 50 I                 = 1, NGPARS
               GPTEMP( NLISGP + I ) = BIG
   50       CONTINUE
C
C  FIND THE STARTING ADDRESSES FOR THE LISTS OF THE GROUP PARAMETERS.
C
            NLISGP = NLISGP + NGPARS
            RETURN
         ELSE
C
C  THE GROUP IS NEW AND OF DEFAULT TYPE. DETERMINE THE STARTING
C  ADDRESSES FOR THE PARAMETERS FOR THE GROUP AND THE
C  NUMBER OF PARAMETERS INVOLVED.
C
            K                = NDTYPE
C ** Correction 10. 1 line replaced by 5
            IF ( K .EQ. 0 ) THEN
              NGPARS = 0
            ELSE
              NGPARS = IGPA( K + 1 ) - IGPA( K )
            END IF
            ITYPEG( NGRUPE ) = K
            ISTGP ( NGRUPE ) = NLISGP + 1
            IF ( NLISGP + NGPARS .GT. NGPVMX ) THEN
               INFORM = - 18
               RETURN
            END IF
C
C  INITIALIZE THE SET OF PARAMETER VALUES.
C
            DO 55 I                 = 1, NGPARS
               GPTEMP( NLISGP + I ) = BIG
   55       CONTINUE
C
C  FIND THE STARTING ADDRESSES FOR THE LISTS OF THE GROUP PARAMETERS.
C
            NLISGP = NLISGP + NGPARS
         END IF
      END IF
C
C  CHECK THAT THE CARDS ARE IN THE CORRECT ORDER.
C
      IF ( FIELD1 .EQ. 'T ' .OR. FIELD1 .EQ. 'XT' ) THEN
         INFORM = 31
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2310 )
         RETURN
      END IF
C
C  THE CARD CONTAINS NAMES OF NONLINEAR ELEMENTS
C
      IF ( FIELD1 .EQ. 'E ' .OR. FIELD1 .EQ. 'XE' .OR.
     *     FIELD1 .EQ. 'ZE' ) THEN
         IF ( NOVALS .GT. 0 ) THEN
C
C  CHECK THE NAME HAS NOT ALREADY BEEN USED IN THE CURRENT ELEMENT.
C  THE PARAMETER NAME OCCURS IN FIELD3 OR FIELD5.
C
            DO 110 I = 1, NOVALS
               IF ( I .EQ. 1 ) THEN
                  FIELD = FIELD3 // 'EL'
               ELSE
                  FIELD = FIELD5 // 'EL'
               END IF
               CALL HASHC( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD)
               IF ( IFIELD .GT. 0 ) THEN
C
C  THE NELING-TH ELEMENT HAS BEEN ASSIGNED.
C
                  NELING = NELING + 1
                  IF ( NELING .GT. NEGMAX ) THEN
                     INFORM = - 10
                     RETURN
                  END IF
C
C  THE ELEMENT IS IELING( ,1) AND THE GROUP IS GIVEN BY IELING( ,2).
C
                  IELING( NELING, 1 ) = INLIST( IFIELD )
                  IELING( NELING, 2 ) = NGRUPE
C
C  THE ELEMENT IS WEIGHTED BY THE CONSTANT WEIGHT().
C
                  IF ( I .EQ. 1 ) THEN
                     WEIGHT( NELING ) = VALUE4
                  ELSE
                     WEIGHT( NELING ) = VALUE6
                  END IF
C
C  THE ELEMENT NAME IS UNKNOWN.
C
               ELSE
                  INFORM = 43
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2430 )
                  RETURN
               END IF
  110       CONTINUE
         END IF
      ELSE
C
C  THE CARD CONTAINS NAMES AND VALUES OF ELEMENTAL PARAMETERS.
C
         IF ( FIELD1 .EQ. 'P ' .OR. FIELD1 .EQ. 'XP' .OR.
     *        FIELD1 .EQ. 'ZP' ) THEN
C
C  DETERMINE THE NUMBER OF THE GROUP TYPE, K, THE STARTING
C  ADDRESSES FOR THE PARAMETERS FOR THE GROUP AND THE
C  NUMBER OF PARAMETERS INVOLVED.
C
            IF ( ITYPEG( NGRUPE ) .LT. 0 ) 
     *        ITYPEG( NGRUPE ) = - ITYPEG( NGRUPE ) - 1
            ITYPEG( NGRUPE ) = ABS( ITYPEG( NGRUPE ) )
            K      = ITYPEG( NGRUPE )
            NGPARS = IGPA( K + 1 ) - IGPA( K )
            IP     = IGPA( K ) - 1
            MLISGP = ISTGP( NGRUPE ) - 1
            IF ( NOVALS .GT. 0 ) THEN
               DO 230 J = 1, NOVALS
C
C  CHECK THE NAME HAS NOT ALREADY BEEN USED IN THE CURRENT ELEMENT.
C  THE PARAMETER NAME OCCURS IN FIELD3 OR FIELD5.
C
                  DO 210 I = 1, NGPARS
                     IF (
     *               ( J .EQ. 1 .AND. FIELD3 .EQ. GPNAME( IP + I ) )
     *          .OR. ( J .EQ. 2 .AND. FIELD5 .EQ. GPNAME( IP + I ) )
     *                    ) GO TO 220
  210             CONTINUE
C
C  THE GROUP PARAMETER NAME IS NOT RECOGNISED.
C
                  INFORM = 33
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2330 )
                  RETURN
C
C  THE GROUP PARAMETER NAME IS THE I-TH PARAMETER IN THE LIST.
C
  220             CONTINUE
C
C  CHECK THAT THE VALUE HAS NOT ALREADY BEEN SET.
C
                  IF ( GPTEMP( MLISGP + I ) .LT. BIG ) THEN
                     INFORM = 32
                     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2320 )
                     RETURN
                  END IF
C
C  READ THE ASSOCIATED VALUE FROM FIELD4 OR FIELD 6.
C
                  IF ( J .EQ. 1 ) THEN
                     GPTEMP( MLISGP + I ) = VALUE4
                  ELSE
                     GPTEMP( MLISGP + I ) = VALUE6
                  END IF
  230          CONTINUE
            END IF
C
C  FIELD1 NOT RECOGNISED.
C
         ELSE
            INFORM = 10
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2100 ) FIELD1
            RETURN
         END IF
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2040 FORMAT( ' ** Exit from GPSMPS - group/row name not recognised:',
     *        ' name is ', A10 )
 2100 FORMAT( ' ** Exit from GPSMPS - field 1 ', A2,
     *        '  not recognised in GROUP USES section' )
 2190 FORMAT( ' ** Exit from GPSMPS - group type not recognised:',
     *        ' name is ', A10 )
 2310 FORMAT( ' ** Exit from GPSMPS - type for group already set ' )
 2320 FORMAT( ' ** Exit from GPSMPS - group parameter already set ' )
 2330 FORMAT( ' ** Exit from GPSMPS - group parameter unrecognised ' )
 2420 FORMAT( ' ** Exit from GPSMPS - default group type already set ' )
 2430 FORMAT( ' ** Exit from GPSMPS - element name not recognised ' )
      END
C  THIS VERSION: 28/10/1992 AT 02:36:30 PM.
      SUBROUTINE SOBBND( NOBBND, NOBMAX, NRLNDX, LENGTH,
     *                   INLIST, ITABLE, FBOUND, REALVL,
C ** Correction 7. 26/02/01: 2 dummy arguments removed from SOBBND **
     *                   FIELD1, FIELD2, VALUE4, FIELD5,
     *                   OBNAME, KEY   , SINGLE, IOUT  , INFORM )
      INTEGER          IOUT, INFORM, LENGTH, NOBBND, NOBMAX, NRLNDX
      DOUBLE PRECISION VALUE4
      LOGICAL          SINGLE
      CHARACTER * 2    FIELD1
      CHARACTER * 10   FIELD2, FIELD5
      INTEGER          INLIST( LENGTH ), ITABLE ( LENGTH )
      DOUBLE PRECISION FBOUND( 2, NOBMAX ), REALVL( NRLNDX )
      CHARACTER * 10   OBNAME( NOBMAX )
      CHARACTER * 12   KEY( LENGTH )
C
C  INDICATOR CARD IS OBJECT BOUND.
C  -------------------------------
C
C  NICK GOULD 12/01/1990
C  FOR CGT PRODUCTIONS.
C
      INTEGER          IFREE, IFIELD, J
      REAL             SMACHR
      DOUBLE PRECISION DMACHR, BIG
      EXTERNAL         HASHB , HASHC, DMACHR
      INTRINSIC        DBLE
      IF ( SINGLE ) THEN
         BIG = 9.0D-1 * DBLE( SMACHR( 5 ) )
      ELSE
         BIG = 9.0D-1 * DMACHR( 5 )
      END IF
C
C  FIND A PLACE TO INSERT THE OBJECTIVE BOUND NAME IN THE HASH-TABLE.
C
      CALL HASHB ( LENGTH, 12, FIELD2 // 'OB', KEY, ITABLE, IFREE )
C
C  THE NAME ALREADY EXISTS. IT IS THE J-TH NAME IN THE LIST.
C
      IF ( IFREE .LE. 0 ) THEN
         IF ( IFREE .EQ. 0 ) THEN
            INFORM = - 1
            RETURN
         END IF
         J = INLIST( - IFREE )
      ELSE
C
C  THE OBJECTIVE FUNCTION BOUND IS THE NOBBND-TH SPECIFIED.
C
         NOBBND = NOBBND + 1
         IF( NOBBND .GT. NOBMAX ) THEN
            INFORM = - 23
            RETURN
         END IF
         J  = NOBBND
C
C  RECORD THE DEFAULT BOUNDS.
C
         FBOUND( 1, NOBBND ) = - BIG
         FBOUND( 2, NOBBND ) = BIG
C
C  RECORD THE POSITION OF THE NEW BOUND IN THE TABLE AND RECORD ITS
C
         INLIST( IFREE )  = NOBBND
         OBNAME( NOBBND ) = FIELD2
      END IF
C
C  RECORD THE BOUND GIVEN.
C
      IF ( FIELD1( 1: 1 ) .EQ. 'Z' ) THEN
         CALL HASHC ( LENGTH, 12,  FIELD5( 1 : 10 ) // 'RI',
     *                KEY, ITABLE, IFIELD )
         IF ( IFIELD .LE. 0 ) THEN
            INFORM = 3
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD5( 1 : 10 )
            RETURN
         END IF
         VALUE4 = REALVL( INLIST( IFIELD ) )
      END IF
      IF ( FIELD1 .EQ. 'XL' .OR. FIELD1 .EQ. 'ZL' .OR.
     *     FIELD1 .EQ. 'LO' ) FBOUND( 1, J ) = VALUE4
      IF ( FIELD1 .EQ. 'XU' .OR. FIELD1 .EQ. 'ZU' .OR.
     *     FIELD1 .EQ. 'UP' ) FBOUND( 1, J ) = VALUE4
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2030 FORMAT( ' ** Exit from GPSMPS - index parameter name ', A10,
     *        ' not recognised ' )
      END









