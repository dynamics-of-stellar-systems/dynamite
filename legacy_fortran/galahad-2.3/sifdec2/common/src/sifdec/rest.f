C  THIS VERSION: 26/02/2001 AT 09:30:00 AM.
C     ( Last modified on 15 Mar 2001 at 22:28:00 )
C ** Correction report.
C ** Correction -4. 20/04/12: QMATRIX added as alias for QUADOBJ
C ** Correction -3. 21/02/00: QSECTION added as alias for QUADOBJ
C ** Correction -2. 21/02/00: Code to process ZERO-ONE and INTEGER cards added.
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C ** Correction 1. 26/02/01: 1 dummy argument removed from PROCAA **
C ** End of Correction report.
      SUBROUTINE PROCAI( NINDEX, NRLNDX, LENGTH, NUSEIN, NUSERE,
     *                   INFORM, IOUT, LEVEL, NINSTR,
     *                   DEBUG, RVALUE, INLIST, ITABLE,
     *                   NAMIIN, NAMRIN, INSTR, KEY,
     *                   FIELD1, FIELD2, FIELD3, FIELD5, FIELD4 )
      INTEGER       NINDEX, NRLNDX, LENGTH, NUSEIN, NUSERE
      INTEGER       INFORM, IOUT, LEVEL, NINSTR
      LOGICAL       DEBUG
      DOUBLE PRECISION RVALUE
      CHARACTER *  2 FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      CHARACTER * 12 FIELD4
      INTEGER       INSTR( 5 ), INLIST( LENGTH ), ITABLE ( LENGTH )
      CHARACTER * 10 NAMIIN( NINDEX ), NAMRIN( NRLNDX )
      CHARACTER * 12 KEY( LENGTH )
C
C  CONSTRUCT A LIST OF DO-LOOP INTEGER AND REAL
C  ARITHMETIC INSTRUCTIONS.
C
C  NICK GOULD, 8/09/89
C  FOR CGT PRODUCTIONS.
C
      INTEGER        NFUNCT, I, IFIELD, IFREE
      CHARACTER * 12 FIELD
      EXTERNAL       HASHB , HASHC, GETINT, GETVAL
      PARAMETER      ( NFUNCT = 14 )
      CHARACTER * 6  FUNCTN( NFUNCT )
      DATA  FUNCTN / 'ABS   ', 'SQRT  ', 'EXP   ', 'LOG   ', 'LOG10 ',
     *               'SIN   ', 'COS   ', 'TAN   ', 'ARCSIN', 'ARCCOS',
     *               'ARCTAN', 'HYPSIN', 'HYPCOS', 'HYPTAN' /
C
C  DECIDE WHAT SORT OF INSTRUCTION IS TO BE PERFORMED.
C  INTEGER INSTRUCTIONS.
C
      IF ( FIELD1 .EQ. 'IE' .OR. FIELD1 .EQ. 'IA' .OR.
     *     FIELD1 .EQ. 'IS' .OR. FIELD1 .EQ. 'IM' .OR.
     *     FIELD1 .EQ. 'ID' .OR. FIELD1 .EQ. 'IR' .OR.
     *     FIELD1 .EQ. 'I=' .OR. FIELD1 .EQ. 'I+' .OR.
     *     FIELD1 .EQ. 'I-' .OR. FIELD1 .EQ. 'I*' .OR.
     *     FIELD1 .EQ. 'I/' ) THEN
         IF ( FIELD1 .EQ. 'IE' ) INSTR( 1 ) = 21
         IF ( FIELD1 .EQ. 'IA' ) INSTR( 1 ) = 22
         IF ( FIELD1 .EQ. 'IS' ) INSTR( 1 ) = 23
         IF ( FIELD1 .EQ. 'IM' ) INSTR( 1 ) = 24
         IF ( FIELD1 .EQ. 'ID' ) INSTR( 1 ) = 25
         IF ( FIELD1 .EQ. 'IR' ) INSTR( 1 ) = 26
         IF ( FIELD1 .EQ. 'I=' ) INSTR( 1 ) = 31
         IF ( FIELD1 .EQ. 'I+' ) INSTR( 1 ) = 32
         IF ( FIELD1 .EQ. 'I-' ) INSTR( 1 ) = 33
         IF ( FIELD1 .EQ. 'I*' ) INSTR( 1 ) = 34
         IF ( FIELD1 .EQ. 'I/' ) INSTR( 1 ) = 35
         IF ( FIELD1 .EQ. 'IE' .OR. FIELD1 .EQ. 'IA' .OR.
     *        FIELD1 .EQ. 'IS' .OR. FIELD1 .EQ. 'IM' .OR.
     *        FIELD1 .EQ. 'ID' ) THEN
C
C  READ THE INTEGER VALUE, IVALUE, FROM FIELD 4.
C
            CALL GETINT( FIELD4, INSTR( 4 ) )
         ELSE
C
C  OBTAIN THE INTEGER VALUE, IVALUE, AS THE VALUE OF THE INDEX IN
C  FIELD 5, FIRST ENSURING THAT THE INDEX EXISTS.
C
            IF ( FIELD1 .NE. 'I=' .AND. FIELD1 .NE. 'IR' ) THEN
               FIELD = FIELD5( 1 : 10 ) // 'II'
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .LE. 0 ) THEN
                  INFORM = 3
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
                  RETURN
               END IF
               INSTR( 4 ) = INLIST( IFIELD )
            END IF
         END IF
C
C  IF A DEFINITION IS TO BE MADE FROM A PREVIOUSLY DEFINED INDEX,
C  ENSURE THAT THE INDEX EXISTS.
C
         IF ( FIELD1 .NE. 'IE' ) THEN
            IF ( FIELD1 .NE. 'IR' ) THEN
               FIELD = FIELD3( 1 : 10 ) // 'II'
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .LE. 0 ) THEN
                  INFORM = 3
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
                  RETURN
               END IF
               INSTR( 3 ) = INLIST( IFIELD )
            ELSE
               FIELD = FIELD3( 1 : 10 ) // 'RI'
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .LE. 0 ) THEN
                  INFORM = 3
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
                  RETURN
               END IF
               INSTR( 3 ) = INLIST( IFIELD )
            END IF
         END IF
C
C  RECORD THE ADDRESS OF THE INDEX WHICH IS TO BE SET.
C
         FIELD = FIELD2( 1 : 10 ) // 'II'
         CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               RETURN
            END IF
            IFREE = - IFREE
         ELSE
            NUSEIN = NUSEIN + 1
            IF ( NUSEIN .GT. NINDEX ) THEN
               INFORM = - 21
               RETURN
            END IF
            INLIST( IFREE )  = NUSEIN
            NAMIIN( NUSEIN ) = FIELD( 1 : 10 )
         END IF
         INSTR( 2 ) = INLIST( IFREE )
C
C  PRINT DETAILS OF THE INSTRUCTION.
C
         IF ( DEBUG .AND. IOUT .GT. 0 ) THEN
            IF ( INSTR( 1 ) .EQ. 21 )
     *         WRITE( IOUT, 4030 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), INSTR( 4 )
            IF ( INSTR( 1 ) .EQ. 22 )
     *         WRITE( IOUT, 4040 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), NAMIIN( INSTR( 3 ) ), INSTR( 4 )
            IF ( INSTR( 1 ) .EQ. 23 )
     *         WRITE( IOUT, 4041 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), NAMIIN( INSTR( 3 ) ), INSTR( 4 )
            IF ( INSTR( 1 ) .EQ. 24 )
     *         WRITE( IOUT, 4050 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), NAMIIN( INSTR( 3 ) ), INSTR( 4 )
            IF ( INSTR( 1 ) .EQ. 25 )
     *         WRITE( IOUT, 4051 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), INSTR( 4 ), NAMIIN( INSTR( 3 ) )
            IF ( INSTR( 1 ) .EQ. 26 )
     *         WRITE( IOUT, 4055 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), NAMRIN( INSTR( 3 ) )
            IF ( INSTR( 1 ) .EQ. 31 )
     *         WRITE( IOUT, 4059 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), NAMIIN( INSTR( 3 ) )
            IF ( INSTR( 1 ) .EQ. 32 )
     *         WRITE( IOUT, 4060 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), NAMIIN( INSTR( 3 ) ),
     *         NAMIIN( INSTR( 4 ) )
            IF ( INSTR( 1 ) .EQ. 33 )
     *         WRITE( IOUT, 4061 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), NAMIIN( INSTR( 4 ) ),
     *         NAMIIN( INSTR( 3 ) )
            IF ( INSTR( 1 ) .EQ. 34 )
     *         WRITE( IOUT, 4070 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), NAMIIN( INSTR( 3 ) ),
     *         NAMIIN( INSTR( 4 ) )
            IF ( INSTR( 1 ) .EQ. 35 )
     *         WRITE( IOUT, 4071 ) LEVEL, NINSTR,
     *         NAMIIN( INSTR( 2 ) ), NAMIIN( INSTR( 3 ) ),
     *         NAMIIN( INSTR( 4 ) )
         END IF
      ELSE
C
C  REAL INSTRUCTIONS.
C
         IF ( FIELD1 .EQ. 'RE' ) INSTR( 1 ) = 51
         IF ( FIELD1 .EQ. 'RA' ) INSTR( 1 ) = 52
         IF ( FIELD1 .EQ. 'RS' ) INSTR( 1 ) = 53
         IF ( FIELD1 .EQ. 'RM' ) INSTR( 1 ) = 54
         IF ( FIELD1 .EQ. 'RD' ) INSTR( 1 ) = 55
         IF ( FIELD1 .EQ. 'RI' ) INSTR( 1 ) = 56
         IF ( FIELD1 .EQ. 'RF' ) INSTR( 1 ) = 57
         IF ( FIELD1 .EQ. 'R=' ) INSTR( 1 ) = 61
         IF ( FIELD1 .EQ. 'R+' ) INSTR( 1 ) = 62
         IF ( FIELD1 .EQ. 'R-' ) INSTR( 1 ) = 63
         IF ( FIELD1 .EQ. 'R*' ) INSTR( 1 ) = 64
         IF ( FIELD1 .EQ. 'R/' ) INSTR( 1 ) = 65
         IF ( FIELD1 .EQ. 'R(' ) INSTR( 1 ) = 67
         IF ( FIELD1 .EQ. 'RE' .OR. FIELD1 .EQ. 'RA' .OR.
     *        FIELD1 .EQ. 'RS' .OR. FIELD1 .EQ. 'RM' .OR.
     *        FIELD1 .EQ. 'RD' .OR. FIELD1 .EQ. 'RF' ) THEN
C
C  READ THE REAL VALUE, RVALUE, FROM FIELD 4.
C
            CALL GETVAL( FIELD4, RVALUE )
         END IF
         IF ( FIELD1 .EQ. 'R+' .OR. FIELD1 .EQ. 'R-' .OR.
     *        FIELD1 .EQ. 'R*' .OR. FIELD1 .EQ. 'R/' .OR.
     *        FIELD1 .EQ. 'R(' ) THEN
C
C  OBTAIN THE REAL VALUE, RVALUE, AS THE VALUE ASSOCIATED WITH THE REAL
C  INDEX IN FIELD 5, FIRST ENSURING THAT THE INDEX EXISTS.
C
            FIELD = FIELD5( 1 : 10 ) // 'RI'
            CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 3
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
               RETURN
            END IF
            INSTR( 4 ) = INLIST( IFIELD )
         END IF
C
C  IF A DEFINITION IS TO BE MADE FROM A PREVIOUSLY DEFINED INDEX,
C  ENSURE THAT THE INDEX EXISTS.
C
         IF ( FIELD1 .EQ. 'RI' ) THEN
            FIELD = FIELD3( 1 : 10 ) // 'II'
            CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 3
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
               RETURN
            END IF
            INSTR( 3 ) = INLIST( IFIELD )
         ELSE
            IF ( FIELD1 .NE. 'RF' .AND. FIELD1 .NE. 'R(' .AND.
     *           FIELD1 .NE. 'RE' ) THEN
               FIELD = FIELD3( 1 : 10 ) // 'RI'
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               IF ( IFIELD .LE. 0 ) THEN
                  INFORM = 3
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
                  RETURN
               END IF
               INSTR( 3 ) = INLIST( IFIELD )
            ELSE
C
C  THE VALUE IS TO BE OBTAINED USING A SPECIAL FUNCTION. DETERMINE
C  WHICH ONE.
C
               IF ( FIELD1 .NE. 'RE' ) THEN
                  DO 10 I = 1, NFUNCT
                     IF ( FIELD3( 1 : 10 ) .EQ. FUNCTN( I ) ) GO TO 20
   10             CONTINUE
                  INFORM = 39
                  IF ( IOUT .GT. 0 ) WRITE( IOUT, 2390 ) FIELD3( 1 : 10)
                  RETURN
   20             CONTINUE
                  INSTR( 3 ) = I
               END IF
            END IF
         END IF
C
C  RECORD THE ADDRESS OF THE INDEX WHICH IS TO BE SET.
C
         FIELD = FIELD2( 1 : 10 ) // 'RI'
         CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
         IF ( IFREE .LE. 0 ) THEN
            IF ( IFREE .EQ. 0 ) THEN
               INFORM = - 1
               RETURN
            END IF
            IFREE = - IFREE
         ELSE
            NUSERE = NUSERE + 1
            IF ( NUSERE .GT. NRLNDX ) THEN
               INFORM = - 22
               RETURN
            END IF
            INLIST( IFREE )  = NUSERE
            NAMRIN( NUSERE ) = FIELD( 1 : 10 )
         END IF
         INSTR( 2 ) = INLIST( IFREE )
C
C  PRINT DETAILS OF THE INSTRUCTION.
C
         IF ( DEBUG .AND. IOUT .GT. 0 ) THEN
            IF ( INSTR( 1 ) .EQ. 51 )
     *         WRITE( IOUT, 4130 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), RVALUE
            IF ( INSTR( 1 ) .EQ. 52 )
     *         WRITE( IOUT, 4140 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), NAMRIN( INSTR( 3 ) ), RVALUE
            IF ( INSTR( 1 ) .EQ. 53 )
     *         WRITE( IOUT, 4141 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), NAMRIN( INSTR( 3 ) ), RVALUE
            IF ( INSTR( 1 ) .EQ. 54 )
     *         WRITE( IOUT, 4150 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), NAMRIN( INSTR( 3 ) ), RVALUE
            IF ( INSTR( 1 ) .EQ. 55 )
     *         WRITE( IOUT, 4151 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), RVALUE, NAMRIN( INSTR( 3 ) )
            IF ( INSTR( 1 ) .EQ. 56 )
     *         WRITE( IOUT, 4180 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), NAMIIN( INSTR( 3 ) )
            IF ( INSTR( 1 ) .EQ. 57 )
     *         WRITE( IOUT, 4110 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), FUNCTN( INSTR( 3 ) ), RVALUE
            IF ( INSTR( 1 ) .EQ. 61 )
     *         WRITE( IOUT, 4159 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), NAMRIN( INSTR( 3 ) )
            IF ( INSTR( 1 ) .EQ. 62 )
     *         WRITE( IOUT, 4160 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), NAMRIN( INSTR( 3 ) ),
     *         NAMRIN( INSTR( 4 ) )
            IF ( INSTR( 1 ) .EQ. 63 )
     *         WRITE( IOUT, 4161 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), NAMRIN( INSTR( 4 ) ),
     *         NAMRIN( INSTR( 3 ) )
            IF ( INSTR( 1 ) .EQ. 64 )
     *         WRITE( IOUT, 4170 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), NAMRIN( INSTR( 3 ) ),
     *         NAMRIN( INSTR( 4 ) )
            IF ( INSTR( 1 ) .EQ. 65 )
     *         WRITE( IOUT, 4171 ) LEVEL, NINSTR,
     *         NAMRIN( INSTR( 2 ) ), NAMRIN( INSTR( 3 ) ),
     *         NAMRIN( INSTR( 4 ) )
            IF ( INSTR( 1 ) .EQ. 67 )
     *         WRITE( IOUT, 4120 ) LEVEL, NINSTR, NAMRIN( INSTR( 2 ) ),
     *         FUNCTN( INSTR( 3 ) ), NAMRIN( INSTR( 4 ) )
         END IF
      END IF
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2030 FORMAT( ' ** Exit from GPSMPS - index parameter name ', A10,
     *        ' not recognised ' )
 2390 FORMAT( ' ** Exit from GPSMPS - specified function name ', A10,
     *        ' not recognised ' )
 4030 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' to the value ', I6 )
 4040 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by adding ', A10, ' to the value ', I6 )
 4041 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by subtracting ', A10, ' from the value ', I6 )
 4050 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by multiplying ', A10, ' by the value ', I6 )
 4051 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by dividing the value ', I6, ' by ', A10 )
 4055 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' to the integer equivalent of ', A10 )
 4059 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' to ', A10 )
 4060 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by adding ', A10, ' to ', A10 )
 4061 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by subtracting ', A10, ' from ', A10 )
 4070 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by multiplying ', A10, ' and ', A10 )
 4071 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by dividing ', A10, ' by ', A10 )
 4110 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' to the value ', A6, '(', 1P, D12.4, ')' )
 4120 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' to the value ', A6, '(', A10, ')' )
 4130 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' to the value ', 1P, D12.4 )
 4140 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by adding ', A10, ' to the value ', 1P, D12.4 )
 4141 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by subtracting ', A10, ' from the value ', 1P, D12.4 )
 4150 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by multiplying ', A10, ' by the value ', 1P, D12.4 )
 4151 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by dividing the value ', 1P, D12.4, ' by ', A10 )
 4159 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' to ', A10 )
 4160 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by adding ', A10, ' to ', A10 )
 4161 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by subtracting ', A10, ' from ', A10 )
 4170 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by multiplying ', A10, ' and ', A10 )
 4171 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' by dividing ', A10, ' by ', A10 )
 4180 FORMAT( ' Level ', I2, ' instruction ', I4, ' set ', A10,
     *        ' to the fl. pt. value of ', A10 )
C
C  END OF PROCAI.
C
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE PROCAD( NINDEX, NRLNDX, LEVEL, NINSTR, NUSERE,
     *                   LENGTH, NARRAY, INTYPE, INFORM, IOUT,
     *                   DEBUG, GRP1ST,
     *                   FIELD1, FIELD2, FIELD3, FIELD5,
     *                   FIELD4, FIELD6,
     *                   INLIST, INSTR, ITABLE, IARRAY,
     *                   VARRAY, FARRAY, NAMIIN, NAMRIN,
     *                   ARRAY, CARRAY, KEY )
      INTEGER        NINDEX, NRLNDX, LEVEL, LENGTH, INFORM, IOUT
      INTEGER        NINSTR, NUSERE, INTYPE, NARRAY
      LOGICAL        DEBUG, GRP1ST
      CHARACTER *  2 FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      CHARACTER * 12 FIELD4, FIELD6
      INTEGER        INSTR( 5 ), IARRAY( 5, 3 )
      INTEGER        INLIST( LENGTH ), ITABLE ( LENGTH )
      DOUBLE PRECISION VARRAY( 2 )
      CHARACTER *  2 FARRAY
      CHARACTER * 10 NAMIIN( NINDEX ), NAMRIN( NRLNDX )
      CHARACTER * 10 ARRAY( 3 ), CARRAY( 2 )
      CHARACTER * 12 KEY( LENGTH )
C
C  CONSTRUCT A LIST OF DO-LOOP ARRAY DEFINITIONS.
C
C  NICK GOULD, 8/09/89
C  FOR CGT PRODUCTIONS.
C
      INTEGER        I, KINDAR, IFREE, MFREE, MFIXED
      INTEGER        MBLANK, MNAME, MROWS, MGROUP, MCNSTR, MCOLS, MVARS
      INTEGER        MCONST, MRHS, MRHSP, MRANGE, MBOUND, MSTART
      INTEGER        METYPE, MGTYPE, MEUSES, MGUSES, MOBBND, MENDAT
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C ** Correction -3. 21/02/00: QSECTION added as alias for QUADOBJ
C ** Correction -4. 20/04/12: QMATRIX added as alias for QUADOBJ
      INTEGER        MQHESS, MQUADO, MQUADR, MQUADS, MQSECT, MQMATR
      EXTERNAL       INTFIE, GETVAL, HASHB
C
C  PARAMETER DEFINITIONS.
C
      PARAMETER        ( MBLANK =  1, MFIXED =  2, MFREE  = 3  )
      PARAMETER        ( MNAME  =  4, MROWS  =  5 )
      PARAMETER        ( MGROUP =  6, MCNSTR =  7, MCOLS  =  8 )
      PARAMETER        ( MVARS  =  9, MCONST = 10, MRHS   = 11 )
      PARAMETER        ( MRHSP  = 12, MRANGE = 13, MBOUND = 14 )
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
      PARAMETER        ( MSTART = 15, MQHESS = 16, MQUADR = 17 )
C ** Correction -3. 21/02/00: QSECTION added as alias for QUADOBJ
      PARAMETER        ( MQUADS = 18, MQUADO = 19, MQSECT = 20 )
C ** Correction -4. 20/04/12: QMATRIX added as alias for QUADOBJ
      PARAMETER        ( MQMATR = 21 )
      PARAMETER        ( METYPE = 22, MEUSES = 23, MGTYPE = 24 )
      PARAMETER        ( MGUSES = 25, MOBBND = 26, MENDAT = 27 )
C
C  DETERMINE HOW MUCH INFORMATION MUST BE SAVED BY DETERMINING
C  THE KIND OF ARRAY DEFINITION BEING MADE.
C
      KINDAR = - 1
C
C  REAL INDEX ARRAY DEFINITIONS.
C
      IF ( FIELD1 .EQ. 'AE' .OR. FIELD1 .EQ. 'AA' .OR.
     *     FIELD1 .EQ. 'AS' .OR. FIELD1 .EQ. 'AM' .OR.
     *     FIELD1 .EQ. 'AD' .OR. FIELD1 .EQ. 'AI' .OR.
     *     FIELD1 .EQ. 'A=' .OR. FIELD1 .EQ. 'A+' .OR.
     *     FIELD1 .EQ. 'A-' .OR. FIELD1 .EQ. 'A*' .OR.
     *     FIELD1 .EQ. 'A/' .OR. FIELD1 .EQ. 'AF' .OR.
     *     FIELD1 .EQ. 'A(' ) THEN
         IF ( FIELD1 .EQ. 'AE' ) KINDAR = 107
         IF ( FIELD1 .EQ. 'AF' ) KINDAR = 108
         IF ( FIELD1 .EQ. 'A(' ) KINDAR = 105
         IF ( FIELD1 .EQ. 'AA' ) KINDAR = 101
         IF ( FIELD1 .EQ. 'AS' ) KINDAR = 101
         IF ( FIELD1 .EQ. 'AM' ) KINDAR = 101
         IF ( FIELD1 .EQ. 'AD' ) KINDAR = 101
         IF ( FIELD1 .EQ. 'A=' ) KINDAR = 103
         IF ( FIELD1 .EQ. 'A+' ) KINDAR = 104
         IF ( FIELD1 .EQ. 'A-' ) KINDAR = 104
         IF ( FIELD1 .EQ. 'A*' ) KINDAR = 104
         IF ( FIELD1 .EQ. 'A/' ) KINDAR = 104
         IF ( FIELD1 .EQ. 'AI' ) KINDAR = 106
      ELSE
C
C  GROUPS SECTION.
C
         IF ( INTYPE .EQ. MROWS  .OR. INTYPE .EQ. MGROUP .OR.
     *        INTYPE .EQ. MCNSTR ) THEN
            IF ( FIELD1( 2: 2 ) .EQ. 'N' .OR.
     *           FIELD1( 2: 2 ) .EQ. 'G' .OR.
     *           FIELD1( 2: 2 ) .EQ. 'L' .OR.
     *           FIELD1( 2: 2 ) .EQ. 'E' ) THEN
               IF ( GRP1ST ) THEN
                  IF ( FIELD3( 1: 7 ) .EQ. '''SCALE''' ) THEN
                     IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
                        KINDAR = 115
                     ELSE
                        KINDAR = 108
                     END IF
                  ELSE
                     KINDAR = 100
                  END IF
               ELSE
                  IF ( FIELD3( 1: 7 ) .EQ. '''SCALE''' ) THEN
                     IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
                        KINDAR = 115
                     ELSE
                        KINDAR = 108
                     END IF
                  ELSE
                     IF ( FIELD3( 1: 10 ) .EQ. '          ' ) THEN
                        KINDAR = 100
                     ELSE
                        IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
                           KINDAR = 113
                        ELSE
                           IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                              KINDAR = 101
                           ELSE
                              KINDAR = 102
                           END IF
                        END IF
                     END IF
                  END IF
               END IF
            END IF
         END IF
C
C  VARIABLES SECTION.
C
         IF ( INTYPE .EQ. MCOLS  .OR. INTYPE .EQ. MVARS ) THEN
            IF ( GRP1ST ) THEN
               IF ( FIELD3( 1: 7 ) .EQ. '''SCALE''' ) THEN
                  IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
                     KINDAR = 115
                  ELSE
                     KINDAR = 108
                  END IF
C ** Correction -2a. 21/02/00: Code to process ZERO-ONE and INTEGER cards added.
               ELSE IF ( FIELD3 .EQ. '''ZERO-ONE''' .OR. 
     *                   FIELD3 .EQ. '''INTEGER'' ' ) THEN
                  KINDAR = 106
               ELSE
                  IF ( FIELD3( 1: 10 ) .EQ. '          ' ) THEN
                     KINDAR = 100
                  ELSE
                     IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
                        KINDAR = 113
                     ELSE
                        IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                           KINDAR = 101
                        ELSE
                           KINDAR = 102
                        END IF
                     END IF
                  END IF
               END IF
            ELSE
               IF ( FIELD3( 1: 7 ) .EQ. '''SCALE''' ) THEN
                  IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
                     KINDAR = 115
                  ELSE
                     KINDAR = 108
                  END IF
C ** Correction -2b. 21/02/00: Code to process ZERO-ONE and INTEGER cards added.
               ELSE IF ( FIELD3 .EQ. '''ZERO-ONE''' .OR. 
     *                   FIELD3 .EQ. '''INTEGER'' ' ) THEN
                  KINDAR = 106
               ELSE
                  KINDAR = 100
               END IF
            END IF
         END IF
C
C  CONSTANTS SECTION.
C
         IF ( INTYPE .EQ. MCONST .OR. INTYPE .EQ. MRHS  .OR.
     *        INTYPE .EQ. MRHSP  ) THEN
            IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
               KINDAR = 116
            ELSE
               IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                  KINDAR = 111
               ELSE
                  KINDAR = 112
               END IF
            END IF
         END IF
C
C  RANGES SECTION.
C
         IF ( INTYPE .EQ. MRANGE ) THEN
            IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
               KINDAR = 116
            ELSE
               IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                  KINDAR = 111
               ELSE
                  KINDAR = 112
               END IF
            END IF
         END IF
C
C  BOUNDS SECTION.
C
         IF ( INTYPE .EQ. MBOUND ) THEN
            IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
               KINDAR = 116
            ELSE
               IF ( FIELD1( 2: 2 ) .EQ. 'R' .OR.
     *              FIELD1( 2: 2 ) .EQ. 'M' .OR.
     *              FIELD1( 2: 2 ) .EQ. 'P' ) KINDAR = 110
               IF ( FIELD1( 2: 2 ) .EQ. 'L' .OR.
     *              FIELD1( 2: 2 ) .EQ. 'U' .OR.
     *              FIELD1( 2: 2 ) .EQ. 'X' ) KINDAR = 111
            END IF
         END IF
C
C  START POINT SECTION.
C
         IF ( INTYPE .EQ. MSTART ) THEN
            IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
               KINDAR = 116
            ELSE
               IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                  KINDAR = 111
               ELSE
                  KINDAR = 112
               END IF
            END IF
         END IF
C ** Correction -1. 20/12/99: Code to process QUADOBJ cards added.
C
C  HESSIAN SECTION.
C
C ** Correction -3. 21/02/00: QSECTION added as alias for QUADOBJ
C ** Correction -4. 20/04/12: QMATRIX added as alias for QUADOBJ
         IF ( INTYPE .EQ. MQUADR .OR. INTYPE .EQ. MQUADS .OR. 
     *        INTYPE .EQ. MQUADO .OR. INTYPE .EQ. MQSECT .OR.
     *        INTYPE .EQ. MQHESS .OR. INTYPE .EQ. MQMATR) THEN
            IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
               KINDAR = 113
            ELSE
               IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                  KINDAR = 101
               ELSE
                  KINDAR = 102
               END IF
            END IF
         END IF
C
C  ELEMENT USES SECTION.
C
         IF ( INTYPE .EQ. MEUSES ) THEN
            IF ( FIELD1( 2: 2 ) .EQ. 'T' ) KINDAR = 106
            IF ( FIELD1( 2: 2 ) .EQ. 'V' ) KINDAR = 105
            IF ( FIELD1( 2: 2 ) .EQ. 'P' ) THEN
               IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
                  KINDAR = 115
               ELSE
                  IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                     KINDAR = 108
                  ELSE
                     KINDAR = 109
                  END IF
               END IF
            END IF
         END IF
C
C  GROUP USES SECTION.
C
         IF ( INTYPE .EQ. MGUSES ) THEN
            IF ( FIELD1( 2: 2 ) .EQ. 'T' ) KINDAR = 106
            IF ( FIELD1( 2: 2 ) .EQ. 'E' ) THEN
               IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
                  KINDAR = 113
               ELSE
                  IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                     KINDAR = 101
                  ELSE
                     KINDAR = 102
                     IF ( FIELD6( 1: 12 ) .EQ. '            ' )
     *                    FIELD6( 1: 3 ) = '1.0'
                  END IF
                  IF ( FIELD4( 1: 12 ) .EQ. '            ' )
     *                 FIELD4( 1: 3 ) = '1.0'
               END IF
            END IF
            IF ( FIELD1( 2: 2 ) .EQ. 'P' ) THEN
               IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
                  KINDAR = 115
               ELSE
                  IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                     KINDAR = 108
                  ELSE
                     KINDAR = 109
                  END IF
               END IF
            END IF
         END IF
C
C  RANGES SECTION.
C
         IF ( INTYPE .EQ. MOBBND ) THEN
            IF ( FIELD1( 1: 1 )  .EQ. 'Z' ) THEN
               KINDAR = 116
            ELSE
               IF ( FIELD5( 1: 10 ) .EQ. '          ' ) THEN
                  KINDAR = 111
               ELSE
                  KINDAR = 112
               END IF
            END IF
         END IF
      END IF
C
C  CHECK THAT THE TYPE OF ARRAY DEFINITION HAS BEEN RECOGNISED.
C
      IF ( KINDAR .LT. 0 ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2140 )
         INFORM = 14
         RETURN
      ELSE
         FARRAY     = FIELD1
         INSTR( 1 ) = KINDAR
         INSTR( 2 ) = NARRAY
      END IF
C
C  AN ARRAY NAME OCCURS IN FIELD 2. INTERPRET THE CONTENTS OF THIS
C  FIELD.
C
      IF ( ( KINDAR .GE. 100 .AND. KINDAR .LE. 109 ) .OR.
     *     ( KINDAR .GE. 113 .AND. KINDAR .LE. 115 ) ) THEN
         CALL INTFIE( LENGTH, 12, KEY, ITABLE, INLIST,
     *                FIELD2, ARRAY( 1 ),
     *                IARRAY( 1, 1 ), IOUT, INFORM )
         IF ( INFORM .NE. 0 ) RETURN
         IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4080 ) LEVEL,
     *                    NINSTR, ARRAY( 1 ),
     *                   ( NAMIIN( IARRAY( 2 + I, 1 ) ),
     *                      I = 1, IARRAY( 2, 1 ) )
C
C  IF THE ARRAY NAME IS JUST A SCALAR NAME, RECORD ITS NAME.
C
         IF ( FIELD1 .EQ. 'AE' .OR. FIELD1 .EQ. 'AA' .OR.
     *        FIELD1 .EQ. 'AS' .OR. FIELD1 .EQ. 'AM' .OR.
     *        FIELD1 .EQ. 'AD' .OR. FIELD1 .EQ. 'AI' .OR.
     *        FIELD1 .EQ. 'A=' .OR. FIELD1 .EQ. 'A+' .OR.
     *        FIELD1 .EQ. 'A-' .OR. FIELD1 .EQ. 'A*' .OR.
     *        FIELD1 .EQ. 'A/' .OR. FIELD1 .EQ. 'AF' .OR.
     *        FIELD1 .EQ. 'A(' ) THEN
            IF ( IARRAY( 2, 1 ) .EQ. 0 ) THEN
               CALL HASHB ( LENGTH, 12, ARRAY( 1 ) // 'RI',
     *                      KEY, ITABLE, IFREE )
               IF ( IFREE .LE. 0 ) THEN
                  IF ( IFREE .EQ. 0 ) THEN
                     INFORM = - 1
                     RETURN
                  END IF
               ELSE
                  NUSERE = NUSERE + 1
                  IF ( NUSERE .GT. NRLNDX ) THEN
                     INFORM = - 22
                     RETURN
                  END IF
                  INLIST( IFREE )  = NUSERE
                  NAMRIN( NUSERE ) = ARRAY( 1 )
               END IF
            END IF
         END IF
      END IF
C
C  AN ARRAY NAME OCCURS IN FIELD 3. INTERPRET THE CONTENTS OF THIS
C  FIELD.
C
      IF ( ( KINDAR .GE. 101 .AND. KINDAR .LE. 104 ) .OR.
     *     ( KINDAR .GE. 110 .AND. KINDAR .LE. 112 ) .OR.
     *       KINDAR .EQ. 113 .OR.  KINDAR .EQ. 116 ) THEN
         CALL INTFIE( LENGTH, 12, KEY, ITABLE, INLIST,
     *                FIELD3, ARRAY( 2 ),
     *                IARRAY( 1, 2 ), IOUT, INFORM )
         IF ( INFORM .NE. 0 ) RETURN
         IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4090 ) LEVEL,
     *                    NINSTR, ARRAY( 2 ),
     *                   ( NAMIIN( IARRAY( 2 + I, 2 ) ),
     *                      I = 1, IARRAY( 2, 2 ) )
      END IF
C
C  AN ARRAY NAME OCCURS IN FIELD 5. INTERPRET THE CONTENTS OF THIS
C  FIELD.
C
      IF ( KINDAR .EQ. 102 .OR.  KINDAR .EQ. 104 .OR.
     *     KINDAR .EQ. 105 .OR.  KINDAR .EQ. 112 .OR.
     *   ( KINDAR .GE. 113 .AND. KINDAR .LE. 116 ) ) THEN
         CALL INTFIE( LENGTH, 12, KEY, ITABLE, INLIST,
     *                FIELD5, ARRAY( 3 ),
     *                IARRAY( 1, 3 ), IOUT, INFORM )
         IF ( INFORM .NE. 0 ) RETURN
         IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4100 ) LEVEL,
     *                    NINSTR, ARRAY( 3 ),
     *                   ( NAMIIN( IARRAY( 2 + I, 3 ) ),
     *                     I = 1, IARRAY( 2, 3 ) )
      END IF
C
C  AN NAME OCCURS IN FIELD 2.
C
      IF ( ( KINDAR .GE. 110 .AND. KINDAR .LE. 112 ) .OR.
     *       KINDAR .EQ. 116 ) THEN
         CARRAY( 1 ) = FIELD2
         IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4110 ) LEVEL,
     *                    NINSTR, CARRAY( 1 )
      END IF
C
C  AN NAME OCCURS IN FIELD 3.
C
      IF ( KINDAR .EQ. 105 .OR. KINDAR .EQ. 106 .OR.
     *     KINDAR .EQ. 108 .OR. KINDAR .EQ. 109 .OR.
     *     KINDAR .EQ. 115 ) THEN
         CARRAY( 1 ) = FIELD3
         IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4120 ) LEVEL,
     *                    NINSTR, CARRAY( 1 )
      END IF
C
C  AN NAME OCCURS IN FIELD 5.
C
      IF ( KINDAR .EQ. 109 ) THEN
         CARRAY( 2 ) = FIELD5
         IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4130 ) LEVEL,
     *                    NINSTR, CARRAY( 2 )
      END IF
C
C  A NUMERICAL VALUE OCCURS IN FIELD 4.
C
      IF ( KINDAR .EQ. 101 .OR. KINDAR .EQ. 102 .OR.
     *     KINDAR .EQ. 107 .OR. KINDAR .EQ. 108 .OR.
     *     KINDAR .EQ. 109 .OR. KINDAR .EQ. 111 .OR.
     *     KINDAR .EQ. 112 ) THEN
         CALL GETVAL( FIELD4, VARRAY( 1 ) )
         IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4140 ) LEVEL,
     *                    NINSTR, VARRAY( 1 )
      END IF
C
C  A NUMERICAL VALUE OCCURS IN FIELD 6.
C
      IF ( KINDAR .EQ. 102 .OR. KINDAR .EQ. 109 .OR.
     *     KINDAR .EQ. 112 ) THEN
         CALL GETVAL( FIELD6, VARRAY( 2 ) )
         IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 4150 ) LEVEL,
     *                    NINSTR, VARRAY( 2 )
      END IF
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2140 FORMAT( ' ** Exit from GPSMPS - type of array defn. unrecognised')
 4080 FORMAT( ' Level ', I2, ' instruction ', I4, ' field 2 array ',
     *          A10, ' indices ', 3( A10, 1X ) )
 4090 FORMAT( ' Level ', I2, ' instruction ', I4, ' field 3 array ',
     *          A10, ' indices ', 3( A10, 1X ) )
 4100 FORMAT( ' Level ', I2, ' instruction ', I4, ' field 5 array ',
     *          A10, ' indices ', 3( A10, 1X ) )
 4110 FORMAT( ' Level ', I2, ' instruction ', I4, ' field 2 name ', A10)
 4120 FORMAT( ' Level ', I2, ' instruction ', I4, ' field 3 name ', A10)
 4130 FORMAT( ' Level ', I2, ' instruction ', I4, ' field 5 name ', A10)
 4140 FORMAT( ' Level ', I2, ' instruction ', I4, ' field 4 value ',
     *        1P, D12.4 )
 4150 FORMAT( ' Level ', I2, ' instruction ', I4, ' field 6 value ',
     *        1P, D12.4 )
C
C  END OF PROCAD.
C
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE PROCAA( NINDEX, LENGTH, NUSERE, INFORM, IOUT, NRLNDX,
C ** Correction 8. 26/02/01: 1 dummy argument removed from PROCAA **
     *                   INLIST, ITABLE,
     *                   NAMRIN, KEY, INDVAL, REALVL,
     *                   FIELD1, FIELD2, FIELD3, FIELD5, RVALUE )
      INTEGER        NINDEX, LENGTH, NUSERE, INFORM, IOUT
      INTEGER        NRLNDX
      DOUBLE PRECISION RVALUE
      CHARACTER *  2 FIELD1
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        INLIST( LENGTH ), ITABLE ( LENGTH )
      INTEGER        INDVAL( NINDEX )
      DOUBLE PRECISION REALVL( NRLNDX )
      CHARACTER * 10 NAMRIN( NRLNDX )
      CHARACTER * 12 KEY( LENGTH )
C
C  CONSTRUCT AND EXECUTE A LIST OF DO-LOOP REAL ARRAY
C  ARITHMETIC INSTRUCTIONS.
C
C  NICK GOULD, 8/09/89
C  FOR CGT PRODUCTIONS.
C
      INTEGER        I, IFIELD, IFREE, INSTR2, INSTR3, INSTR4, NFUNCT
      CHARACTER * 12 FIELD
      INTRINSIC      FLOAT
      EXTERNAL       RINTRN, HASHB , HASHC
      PARAMETER      ( NFUNCT = 14 )
      CHARACTER * 6  FUNCTN( NFUNCT )
      DATA  FUNCTN / 'ABS   ', 'SQRT  ', 'EXP   ', 'LOG   ', 'LOG10 ',
     *               'SIN   ', 'COS   ', 'TAN   ', 'ARCSIN', 'ARCCOS',
     *               'ARCTAN', 'HYPSIN', 'HYPCOS', 'HYPTAN' /
      IF ( FIELD1 .EQ. 'A+' .OR. FIELD1 .EQ. 'A-' .OR.
     *     FIELD1 .EQ. 'A*' .OR. FIELD1 .EQ. 'A/' .OR.
     *     FIELD1 .EQ. 'A(' ) THEN
C
C  OBTAIN THE REAL VALUE, RVALUE, AS THE VALUE ASSOCIATED WITH THE REAL
C  INDEX IN FIELD 5, FIRST ENSURING THAT THE INDEX EXISTS.
C
         FIELD = FIELD5( 1 : 10 ) // 'RI'
         CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
         IF ( IFIELD .LE. 0 ) THEN
            INFORM = 3
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
            RETURN
         END IF
         INSTR4 = INLIST( IFIELD )
      END IF
C
C  IF A DEFINITION IS TO BE MADE FROM A PREVIOUSLY DEFINED INDEX,
C  ENSURE THAT THE INDEX EXISTS.
C
      IF ( FIELD1 .EQ. 'AI' ) THEN
         FIELD = FIELD3( 1 : 10 ) // 'II'
         CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
         IF ( IFIELD .LE. 0 ) THEN
            INFORM = 3
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
            RETURN
         END IF
         INSTR3 = INLIST( IFIELD )
      ELSE
         IF ( FIELD1 .NE. 'AF' .AND. FIELD1 .NE. 'A(' .AND.
     *        FIELD1 .NE. 'AE' ) THEN
            FIELD = FIELD3( 1 : 10 ) // 'RI'
            CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 3
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
               RETURN
            END IF
            INSTR3 = INLIST( IFIELD )
         ELSE
C
C  THE VALUE IS TO BE OBTAINED USING A SPECIAL FUNCTION. DETERMINE
C  WHICH ONE.
C
            IF ( FIELD1 .NE. 'AE' ) THEN
               DO 10 I = 1, NFUNCT
                  IF ( FIELD3( 1 : 10 ) .EQ. FUNCTN( I ) ) GO TO 20
   10          CONTINUE
               INFORM = 39
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2390 ) FIELD3( 1 : 10 )
               RETURN
   20          CONTINUE
               INSTR3 = I
            END IF
         END IF
      END IF
C
C  RECORD THE ADDRESS OF THE INDEX WHICH IS TO BE SET.
C
      FIELD = FIELD2( 1 : 10 ) // 'RI'
      CALL HASHB ( LENGTH, 12, FIELD, KEY, ITABLE, IFREE )
      IF ( IFREE .LE. 0 ) THEN
         IF ( IFREE .EQ. 0 ) THEN
            INFORM = - 1
            RETURN
         END IF
         IFREE = - IFREE
      ELSE
         NUSERE = NUSERE + 1
         IF ( NUSERE .GT. NRLNDX ) THEN
            INFORM = - 22
            RETURN
         END IF
         INLIST( IFREE )  = NUSERE
         NAMRIN( NUSERE ) = FIELD( 1 : 10 )
      END IF
      INSTR2 = INLIST( IFREE )
      IF ( FIELD1 .EQ. 'AE' ) REALVL( INSTR2 ) = RVALUE
      IF ( FIELD1 .EQ. 'AA' ) REALVL( INSTR2 ) =
     *                        RVALUE + REALVL( INSTR3 )
      IF ( FIELD1 .EQ. 'AS' ) REALVL( INSTR2 ) =
     *                        RVALUE - REALVL( INSTR3 )
      IF ( FIELD1 .EQ. 'AM' ) REALVL( INSTR2 ) =
     *                        RVALUE * REALVL( INSTR3 )
      IF ( FIELD1 .EQ. 'AD' ) REALVL( INSTR2 ) =
     *                        RVALUE / REALVL( INSTR3 )
      IF ( FIELD1 .EQ. 'AI' ) REALVL( INSTR2 ) =
     *                        FLOAT( INDVAL( INSTR3 ) )
      IF ( FIELD1 .EQ. 'AF' ) CALL RINTRN( REALVL( INSTR2 ),
     *                                     RVALUE, INSTR3, INFORM )
      IF ( FIELD1 .EQ. 'A=' ) REALVL( INSTR2 ) = REALVL( INSTR3 )
      IF ( FIELD1 .EQ. 'A+' ) REALVL( INSTR2 ) =
     *                        REALVL( INSTR3 ) + REALVL( INSTR4 )
      IF ( FIELD1 .EQ. 'A-' ) REALVL( INSTR2 ) =
     *                        REALVL( INSTR3 ) - REALVL( INSTR4 )
      IF ( FIELD1 .EQ. 'A*' ) REALVL( INSTR2 ) =
     *                        REALVL( INSTR3 ) * REALVL( INSTR4 )
      IF ( FIELD1 .EQ. 'A/' ) REALVL( INSTR2 ) =
     *                        REALVL( INSTR3 ) / REALVL( INSTR4 )
      IF ( FIELD1 .EQ. 'A(' ) CALL RINTRN( REALVL( INSTR2 ),
     *                        REALVL( INSTR4 ), INSTR3, INFORM )
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2030 FORMAT( ' ** Exit from GPSMPS - index parameter name ', A10,
     *        ' not recognised ' )
 2390 FORMAT( ' ** Exit from GPSMPS - specified function name ', A10,
     *        ' not recognised ' )
C
C  END OF PROCAA.
C
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE GETFIE( NINDEX, INDVAL, IARRAY, ARRAY, FIELD, INFORM )
      INTEGER       NINDEX, INFORM
      INTEGER       INDVAL( NINDEX )
      INTEGER       IARRAY( 5 )
      CHARACTER * 10 ARRAY, FIELD
C
C  CONSTRUCT AN EXPANDED ARRAY NAME FROM ITS CONSTITUENT PARTS.
C
C  NICK GOULD 02/08/89
C  FOR CGT PRODUCTIONS.
C
      INTEGER       I, INDCES, IVALUE, J, NDIGIT
      CHARACTER * 9 FIELD9
      INTRINSIC     ABS, FLOAT, LOG10
      J              = IARRAY( 1 )
      FIELD          = ARRAY( 1: J )
      J              = J + 1
      INDCES         = IARRAY( 2 )
      DO 10 I        = 1, INDCES
         IVALUE      = INDVAL( IARRAY( 2 + I ) )
         IF ( IVALUE .EQ. 0 ) THEN
            NDIGIT = 1
         ELSE
            NDIGIT = LOG10( ABS( FLOAT( IVALUE ) ) ) + 1
            IF ( IVALUE .LT. 0 ) NDIGIT = NDIGIT + 1
         END IF
         IF ( ( I .LT. INDCES .AND. J + NDIGIT .GT. 10 ) .OR.
     *        ( I .EQ. INDCES .AND. J + NDIGIT .GT. 11 ) ) THEN
            INFORM = 35
            RETURN
         END IF
         WRITE( UNIT = FIELD9, FMT = 2000 ) IVALUE
         FIELD( J: J + NDIGIT - 1 ) = FIELD9( 10 - NDIGIT: 9 )
         J = J + NDIGIT
         IF ( I .LT. INDCES ) THEN
            FIELD( J: J ) = ','
            J = J + 1
         END IF
  10  CONTINUE
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( I9 )
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE GETIIN( NINDEX, INDVAL, NRLNDX, REALVL, INSTR )
      INTEGER NINDEX, NRLNDX
      INTEGER INDVAL( NINDEX ), INSTR( 4 )
      DOUBLE PRECISION REALVL( NRLNDX )
C
C  EXECUTE INTEGER INDEX INSTRUCTIONS.
C
C  NICK GOULD 02/08/89
C  FOR CGT PRODUCTIONS.
C
      IF ( INSTR( 1 ) .EQ. 21 ) INDVAL( INSTR( 2 ) ) = INSTR( 4 )
      IF ( INSTR( 1 ) .EQ. 22 ) INDVAL( INSTR( 2 ) ) =
     *   INSTR( 4 ) + INDVAL( INSTR( 3 ) )
      IF ( INSTR( 1 ) .EQ. 23 ) INDVAL( INSTR( 2 ) ) =
     *   INSTR( 4 ) - INDVAL( INSTR( 3 ) )
      IF ( INSTR( 1 ) .EQ. 24 ) INDVAL( INSTR( 2 ) ) =
     *   INSTR( 4 ) * INDVAL( INSTR( 3 ) )
      IF ( INSTR( 1 ) .EQ. 25 ) INDVAL( INSTR( 2 ) ) =
     *   INSTR( 4 ) / INDVAL( INSTR( 3 ) )
      IF ( INSTR( 1 ) .EQ. 26 ) INDVAL( INSTR( 2 ) ) =
     *   REALVL( INSTR( 3 ) )
      IF ( INSTR( 1 ) .EQ. 31 ) INDVAL( INSTR( 2 ) ) =
     *   INDVAL( INSTR( 3 ) )
      IF ( INSTR( 1 ) .EQ. 32 ) INDVAL( INSTR( 2 ) ) =
     *   INDVAL( INSTR( 3 ) ) + INDVAL( INSTR( 4 ) )
      IF ( INSTR( 1 ) .EQ. 33 ) INDVAL( INSTR( 2 ) ) =
     *   INDVAL( INSTR( 3 ) ) - INDVAL( INSTR( 4 ) )
      IF ( INSTR( 1 ) .EQ. 34 ) INDVAL( INSTR( 2 ) ) =
     *   INDVAL( INSTR( 3 ) ) * INDVAL( INSTR( 4 ) )
      IF ( INSTR( 1 ) .EQ. 35 ) INDVAL( INSTR( 2 ) ) =
     *   INDVAL( INSTR( 3 ) ) / INDVAL( INSTR( 4 ) )
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE GETRIN( NINDEX, NRLNDX, INDVAL, REALVL,
     *                   RVALUE, INSTR, INFORM )
      INTEGER NINDEX, NRLNDX, INFORM
      DOUBLE PRECISION    RVALUE
      INTEGER INSTR( 4 ), INDVAL( NINDEX )
      DOUBLE PRECISION    REALVL( NRLNDX )
C
C  EXECUTE REAL INDEX INSTRUCTIONS.
C
C  NICK GOULD 02/08/89
C  FOR CGT PRODUCTIONS.
C
      INTEGER I
      INTRINSIC FLOAT
      EXTERNAL RINTRN
      INFORM = 0
      I      = INSTR( 1 ) - 50
      GO TO ( 110, 120, 130, 140, 150, 160, 170, 300, 300, 300,
     *        210, 220, 230, 240, 250, 300, 270, 300, 300 ), I
  110 CONTINUE
      REALVL( INSTR( 2 ) ) = RVALUE
      RETURN
  120 CONTINUE
      REALVL( INSTR( 2 ) ) = RVALUE + REALVL( INSTR( 3 ) )
      RETURN
  130 CONTINUE
      REALVL( INSTR( 2 ) ) = RVALUE - REALVL( INSTR( 3 ) )
      RETURN
  140 CONTINUE
      REALVL( INSTR( 2 ) ) = RVALUE * REALVL( INSTR( 3 ) )
      RETURN
  150 CONTINUE
      REALVL( INSTR( 2 ) ) = RVALUE / REALVL( INSTR( 3 ) )
      RETURN
  160 CONTINUE
      REALVL( INSTR( 2 ) ) = FLOAT( INDVAL( INSTR( 3 ) ) )
      RETURN
  170 CONTINUE
      CALL RINTRN( REALVL( INSTR( 2 ) ), RVALUE, INSTR( 3 ), INFORM )
      RETURN
  210 CONTINUE
      REALVL( INSTR( 2 ) ) = REALVL( INSTR( 3 ) )
      RETURN
  220 CONTINUE
      REALVL( INSTR( 2 ) ) = REALVL( INSTR( 3 ) ) + REALVL( INSTR( 4 ) )
      RETURN
  230 CONTINUE
      REALVL( INSTR( 2 ) ) = REALVL( INSTR( 3 ) ) - REALVL( INSTR( 4 ) )
      RETURN
  240 CONTINUE
      REALVL( INSTR( 2 ) ) = REALVL( INSTR( 3 ) ) * REALVL( INSTR( 4 ) )
      RETURN
  250 CONTINUE
      REALVL( INSTR( 2 ) ) = REALVL( INSTR( 3 ) ) / REALVL( INSTR( 4 ) )
      RETURN
  270 CONTINUE
      CALL RINTRN( REALVL( INSTR( 2 ) ),
     *             REALVL( INSTR( 4 ) ), INSTR( 3 ), INFORM )
      RETURN
  300 CONTINUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE RINTRN( EVALUE, VALUE, IVALUE, INFORM )
      INTEGER IVALUE, INFORM
      DOUBLE PRECISION VALUE, EVALUE
C
C  RETURNS THE VALUE OF THE APPROPRIATE INTRINSIC FUNCTION.
C
C  NICK GOULD 10/11/89
C  FOR CGT PRODUCTIONS.
C
      INTRINSIC ABS, SQRT, EXP, LOG, LOG10, SIN, COS, TAN, ASIN, ACOS,
     *          ATAN, SINH, COSH, TANH, DBLE
      DOUBLE PRECISION DVALUE
      DVALUE = DBLE( VALUE )
      GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120,
     *        130, 140 ), IVALUE
   10 CONTINUE
      EVALUE = ABS( DVALUE )
      RETURN
   20 CONTINUE
      IF ( VALUE .LT. 0.0 ) THEN
         INFORM = 40
         WRITE( 6, 2400 ) VALUE, 'SQRT  '
      ELSE
         EVALUE = SQRT( DVALUE )
      END IF
      RETURN
   30 CONTINUE
      EVALUE = EXP( DVALUE )
      RETURN
   40 CONTINUE
      IF ( VALUE .LE. 0.0 ) THEN
         INFORM = 40
         WRITE( 6, 2400 ) VALUE, 'LOG   '
      ELSE
         EVALUE = LOG( DVALUE )
      END IF
      RETURN
   50 CONTINUE
      IF ( VALUE .LE. 0.0 ) THEN
         INFORM = 40
         WRITE( 6, 2400 ) VALUE, 'LOG10 '
      ELSE
         EVALUE = LOG10( DVALUE )
      END IF
      RETURN
   60 CONTINUE
      EVALUE = SIN( DVALUE )
      RETURN
   70 CONTINUE
      EVALUE = COS( DVALUE )
      RETURN
   80 CONTINUE
      EVALUE = TAN( DVALUE )
      RETURN
   90 CONTINUE
      IF ( ABS( VALUE ) .GT. 1.0 ) THEN
         INFORM = 40
         WRITE( 6, 2400 ) VALUE, 'ASIN  '
      ELSE
         EVALUE = ASIN( DVALUE )
      END IF
      RETURN
  100 CONTINUE
      IF ( ABS( VALUE ) .GT. 1.0 ) THEN
         INFORM = 40
         WRITE( 6, 2400 ) VALUE, 'ACOS  '
      ELSE
         EVALUE = ACOS( DVALUE )
      END IF
      RETURN
  110 CONTINUE
      EVALUE = ATAN( DVALUE )
      RETURN
  120 CONTINUE
      EVALUE = SINH( DVALUE )
      RETURN
  130 CONTINUE
      EVALUE = COSH( DVALUE )
      RETURN
  140 CONTINUE
      EVALUE = TANH( DVALUE )
      RETURN
 2400 FORMAT( ' ** Exit from GPSMPS - argument value ', 1P, D9.1,
     *        ' is illegal for function ', A6 )
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE GETINT( FIELD, IVALUE )
      INTEGER        IVALUE
      CHARACTER * 12 FIELD
C
C  READ THE INTEGER NUMBER IVALUE STORED IN THE CHARACTER FIELD.
C
C  NICK GOULD 02/08/89
C  FOR CGT PRODUCTIONS.
C
      INTEGER        I, J
      CHARACTER * 12 FIELD2
C
C  RIGHT-SHIFT THE FIELD, ELIMINATING BLANKS.
C
      FIELD2  = '            '
            J = 12
      DO 10 I = 12, 1, - 1
         IF ( FIELD( I : I ) .EQ. ' ' ) GO TO 10
         FIELD2( J : J ) = FIELD( I : I )
         J = J - 1
   10 CONTINUE
      READ( UNIT = FIELD2, FMT = 2000 ) IVALUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( I12 )
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE GETLIN( NINDEX, NRLNDX, INDVAL, IARRAY, VARRAY,
     *                   ARRAY, CARRAY, FARRAY, REALVL, NAMIIN, NOVALS,
     *                   KINDAR, FIELD1, FIELD2, FIELD3, VALUE4,
     *                   FIELD5, VALUE6, IOUT, INFORM,
     *                   LENGTH, KEY, ITABLE, INLIST )
      INTEGER        NINDEX, NRLNDX, KINDAR, NOVALS, IOUT, INFORM
      INTEGER        LENGTH
      DOUBLE PRECISION VALUE4, VALUE6
      CHARACTER *  2 FIELD1, FARRAY
      CHARACTER * 10 FIELD2, FIELD3, FIELD5
      INTEGER        INDVAL( NINDEX ), IARRAY( 5, 3 )
      INTEGER        INLIST ( LENGTH ), ITABLE ( LENGTH )
      DOUBLE PRECISION VARRAY( 2 ), REALVL( NRLNDX )
      CHARACTER *  7 NAMIIN( NINDEX )
      CHARACTER * 10 ARRAY( 3 ), CARRAY( 2 )
      CHARACTER * 12 KEY   ( LENGTH )
C
C  TRANSLATE THE CONTENTS OF AN ARRAY CARD TO ITS SCALAR VALUES.
C  THE EXPECTED CONTENTS ARE AT THE CONTROL OF THE PARAMETER KINDAR
C  AS FOLLOWS: (N=NAME ,A=ARRAY NAME, V=NUMERICAL VALUE,
C               R=REAL INDEX ARRAY VALUE)
C
C  KINDAR       FIELD2    FIELD3    FIELD4     FIELD5     FIELD6
C  ------       ------    ------    ------     ------     ------
C   100            A
C   101            A         A         V
C   102            A         A         V          A          V
C   103            A         A
C   104            A         A                    A
C   105            A         N                    A
C   106            A         N
C   107            A                   V
C   108            A         N         V
C   109            A         N         V          N          V
C   110            N         A
C   111            N         A         V
C   112            N         A         V          A          V
C   113            A         A          < ------- R
C   114            A                    < ------- R
C   115            A         N          < ------- R
C   116            N         A          < ------- R
C
C  NICK GOULD 02/08/89
C  FOR CGT PRODUCTIONS.
C
      INTEGER  I, IFIELD
      CHARACTER * 12 FIELD
      EXTERNAL GETFIE, HASHC
      FIELD1 = FARRAY
      NOVALS = 0
C
C  AN ARRAY NAME OCCURS IN FIELD 2. INTERPRET THE CONTENTS OF THIS
C  FIELD.
C
      IF ( ( KINDAR .GE. 100 .AND. KINDAR .LE. 109 ) .OR.
     *     ( KINDAR .GE. 113 .AND. KINDAR .LE. 115 ) ) THEN
         CALL GETFIE( NINDEX, INDVAL, IARRAY( 1, 1 ),
     *                ARRAY( 1 ), FIELD2, INFORM )
         IF ( INFORM .NE. 0 ) THEN
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2350 )
     *         ARRAY( 1 )( 1: IARRAY( 1, 1 ) ),
     *         ( NAMIIN( IARRAY( 2 + I, 1 ) ),
     *           INDVAL( IARRAY( 2 + I, 1 ) ), I = 1, IARRAY( 2, 1 ) )
            INFORM = 35
            RETURN
         END IF
      ELSE
         FIELD2 = '          '
      END IF
C
C  AN ARRAY NAME OCCURS IN FIELD 3. INTERPRET THE CONTENTS OF THIS
C  FIELD.
C
      IF ( ( KINDAR .GE. 101 .AND. KINDAR .LE. 104 ) .OR.
     *     ( KINDAR .GE. 110 .AND. KINDAR .LE. 112 ) .OR.
     *       KINDAR .EQ. 113 .OR.  KINDAR .EQ. 116 ) THEN
         CALL GETFIE( NINDEX, INDVAL, IARRAY( 1, 2 ),
     *                ARRAY( 2 ), FIELD3, INFORM )
         IF ( INFORM .NE. 0 ) THEN
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2350 )
     *         ARRAY( 2 )( 1: IARRAY( 1, 2 ) ),
     *         ( NAMIIN( IARRAY( 2 + I, 2 ) ),
     *           INDVAL( IARRAY( 2 + I, 2 ) ), I = 1, IARRAY( 2, 2 ) )
            INFORM = 35
            RETURN
         END IF
      ELSE
         FIELD3 = '          '
      END IF
      IF ( KINDAR .EQ. 103 ) NOVALS = 1
C
C  AN ARRAY NAME OCCURS IN FIELD 5. INTERPRET THE CONTENTS OF THIS
C  FIELD.
C
      IF ( KINDAR .EQ. 102 .OR.  KINDAR .EQ. 104 .OR.
     *     KINDAR .EQ. 105 .OR.  KINDAR .EQ. 112 .OR.
     *   ( KINDAR .GE. 113 .AND. KINDAR .LE. 116 ) ) THEN
         CALL GETFIE( NINDEX, INDVAL, IARRAY( 1, 3 ),
     *                ARRAY( 3 ), FIELD5, INFORM )
         IF ( INFORM .NE. 0 ) THEN
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2350 )
     *         ARRAY( 3 )( 1: IARRAY( 1, 3 ) ),
     *         ( NAMIIN( IARRAY( 2 + I, 3 ) ),
     *           INDVAL( IARRAY( 2 + I, 3 ) ), I = 1, IARRAY( 2, 3 ) )
            INFORM = 35
            RETURN
         END IF
      ELSE
         FIELD5 = '          '
      END IF
      IF ( KINDAR .EQ. 104 ) NOVALS = 2
C
C  AN NAME OCCURS IN FIELD 2.
C
      IF ( ( KINDAR .GE. 110 .AND. KINDAR .LE. 112 ) .OR.
     *       KINDAR .EQ. 116 ) THEN
         FIELD2 = CARRAY( 1 )
      END IF
C
C  AN NAME OCCURS IN FIELD 3.
C
      IF ( KINDAR .EQ. 105 .OR. KINDAR .EQ. 106 .OR.
     *     KINDAR .EQ. 108 .OR. KINDAR .EQ. 109 .OR.
     *     KINDAR .EQ. 115 ) THEN
         FIELD3 = CARRAY( 1 )
      END IF
C
C  AN NAME OCCURS IN FIELD 5.
C
      IF ( KINDAR .EQ. 109 ) THEN
         FIELD5 = CARRAY( 2 )
      END IF
C
C  A NUMERICAL VALUE OCCURS IN FIELD 4.
C
      IF ( KINDAR .EQ. 101 .OR. KINDAR .EQ. 102 .OR.
     *     KINDAR .EQ. 107 .OR. KINDAR .EQ. 108 .OR.
     *     KINDAR .EQ. 109 .OR. KINDAR .EQ. 111 .OR.
     *     KINDAR .EQ. 112 ) THEN
         VALUE4 = VARRAY( 1 )
         NOVALS = 1
      ELSE
C
C  A REAL INDEX VALUE IS TO BE PLACES IN FIELD 4.
C
         IF ( KINDAR .EQ. 113 .OR. KINDAR .EQ. 114 .OR.
     *        KINDAR .EQ. 115 .OR. KINDAR .EQ. 116 ) THEN
            FIELD = FIELD5( 1 : 10 ) // 'RI'
            CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 3
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1 : 10 )
               RETURN
            END IF
            VALUE4 = REALVL( INLIST( IFIELD ) )
            NOVALS = 1
         ELSE
            VALUE4 = 0.0D+0
         END IF
      END IF
C
C  A NUMERICAL VALUE OCCURS IN FIELD 6.
C
      IF ( KINDAR .EQ. 102 .OR. KINDAR .EQ. 109 .OR.
     *     KINDAR .EQ. 112 ) THEN
         VALUE6 = VARRAY( 2 )
         NOVALS = 2
      ELSE
         VALUE6 = 0.0D+0
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2030 FORMAT( ' ** Exit from GPSMPS - index parameter name ', A10,
     *        ' not recognised ' )
 2350 FORMAT( ' ** Exit from GPSMPS - expanded array name > 10 chars.',
     *        ' Array name = ', A9, /, ( '    Index ', A10,
     *        ' has the value ', I6, : ) )
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE GETVAL( FIELD, VALUE )
      DOUBLE PRECISION VALUE
      CHARACTER * 12 FIELD
C
C  READ THE REAL NUMBER VALUE STORED IN THE CHARACTER FIELD.
C
C  NICK GOULD 02/08/89
C  FOR CGT PRODUCTIONS.
C
      READ( UNIT = FIELD, FMT = 2000 ) VALUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( BN, F12.0 )
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE INTFIE( LENGTH, NCHAR, KEY, ITABLE, INLIST,
     *                   FIELDA, ARRAY, IARRAY, IOUT, INFORM )
      INTEGER        LENGTH, NCHAR, IOUT, INFORM
      CHARACTER * 10 FIELDA, ARRAY
      CHARACTER * 12 KEY   ( LENGTH )
      INTEGER        IARRAY( 5 )
      INTEGER       INLIST ( LENGTH )
      INTEGER       ITABLE ( LENGTH )
C
C  INTERPRET THE CONTENTS OF FIELDA.
C
C  NICK GOULD 02/08/89
C  FOR CGT PRODUCTIONS.
C
      INTEGER        I, IFIELD, INDCES, J
      CHARACTER * 12 FIELD
      EXTERNAL       HASHC
C
C  FIRST FIND THE ARRAY NAME BY SEARCHING THE FIELD FOR THE STRING '('.
C
      DO 10 I = 1, 10
         IF ( FIELDA( I: I ) .EQ. '(' ) GO TO 20
   10 CONTINUE
C     IF ( IOUT .GT. 0 ) WRITE( IOUT, 2370 ) FIELDA
C     INFORM = 37
      IARRAY( 1 ) = 10
      IARRAY( 2 ) = 0
      ARRAY  = FIELDA
      INFORM = 0
      RETURN
C
C  THE STRING '(' OCCURS IN POSITION I.
C
   20 CONTINUE
      J           = I - 1
      IARRAY( 1 ) = J
      ARRAY       = FIELDA( 1: J )
      INDCES      = 0
C
C  NOW FIND THE ARRAY INDICES. SEARCH FOR ONE OF THE STRINGS ')' OR ','.
C
   30 CONTINUE
         J       = I + 1
         DO 40 I = J, 10
            IF ( FIELDA( I: I ) .EQ. ',' .OR.
     *           FIELDA( I: I ) .EQ. ')' ) GO TO 50
   40    CONTINUE
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2370 ) FIELDA
         INFORM = 37
         RETURN
C
C  THE STRING ',' OR ')' OCCURS IN POSITION J.
C
   50    CONTINUE
         IF ( I .NE. J ) THEN
            INDCES = INDCES + 1
            IF ( INDCES .GT. 3 ) THEN
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2360 )
               INFORM = 36
               RETURN
            END IF
            FIELD( 1: 12 )    = '          II'
            FIELD( 1: I - J ) = FIELDA( J: I - 1 )
C
C  CHECK THAT THE ARRAY INDEX EXISTS AND DETERMINE ITS ADDRESS.
C
            CALL HASHC ( LENGTH, NCHAR, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .LE. 0 ) THEN
               INFORM = 3
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 ) FIELD( 1: 10 )
               RETURN
            END IF
            IARRAY( 2 + INDCES ) = INLIST( IFIELD )
         END IF
         IF ( FIELDA( I: I ) .EQ. ',' ) GO TO 30
C
C  THE ARRAY DEFINITION IS COMPLETE.
C
      IARRAY( 2 ) = INDCES
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2030 FORMAT( ' ** Exit from GPSMPS - index parameter name ', A10,
     *        ' not recognised ' )
 2370 FORMAT( ' ** Exit from GPSMPS - incorrect array name', A10,
     *        ' in do-loop ')
 2360 FORMAT( ' ** Exit from GPSMPS - > 3 array name indices ' )
C
C  END OF INTFIE.
C
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      SUBROUTINE FREEFM( NULINE, LELINE, MENDAT, INDIC8, LENIND, NULINA,
     *                   MAXNUL, NLINES, ISSDIF, INFORM, IOUT )
      INTEGER          LELINE, MENDAT, MAXNUL, NLINES, INFORM, IOUT
      LOGICAL          ISSDIF
      INTEGER          LENIND( MENDAT )
      CHARACTER * 1    NULINE( LELINE )
      CHARACTER * 12   INDIC8( MENDAT )
      CHARACTER * 65   NULINA( MAXNUL )
C
C  CONSTRUCT A FIXED FORMAT LINE FROM A FREE FORMAT ONE.
C
C  NICK GOULD, 12/07/90
C  FOR CGT PRODUCTIONS.
C
      INTEGER I, J, K, NFIELD, NSTFIE, LFIELD, LEN, ICARD, LELIP1
      LOGICAL FIELD, NEXTL, WAITNL, SPSAL2
      INTEGER LENFIE( 6, 2 ), ISTFIE( 6, 2 ), NFIE( 2 )
      DATA    LENFIE / 2, 10, 10, 12, 10, 12,
     *                 2, 10, 10, 41, 0, 0 /
      DATA    ISTFIE / 2, 5, 15, 25, 40, 50,
     *                 2, 5, 15, 25, 65, 65 /
      DATA    NFIE / 6, 4  /
C
C  IF ISSDIF IS .TRUE., THE CALL IS MADE FROM SUBROUTINE GPS, WHERE THE
C  CARD LENGTH IS 61. OTHERWISE, THE CARD LENGTH IS 65.
C
      IF ( ISSDIF ) THEN
         ICARD  = 1
      ELSE
         ICARD  = 2
      END IF
      LELIP1 = LELINE + 1
      NFIELD = 0
      NLINES = 0
      FIELD  = .FALSE.
      NEXTL  = .FALSE.
      SPSAL2 = .FALSE.
      WAITNL = .TRUE.
C
C  PROCESS THE NEXT CHARACTER ON THE CARD.
C
      DO 500 I = 1, LELINE
C
C  COPY COMMENTS UNCHANGED.
C
         IF ( WAITNL .AND. NULINE( I ) .EQ. '$' ) THEN
            NLINES = NLINES + 1
            DO 10 J = 1, 65
               NULINA( NLINES )( J: J ) = ' '
   10       CONTINUE
            NULINA( NLINES )( 1: 1 ) = '*'
            DO 20 J = 2, LELIP1 - I
               NULINA( NLINES )( J: J ) = NULINE( I + J - 1 )
   20       CONTINUE
            GO TO 600
         END IF
C
C  IF WE ARE LOOKING FOR AN END OF LINE MARKER, CHECK WHETHER WE HAVE
C  FOUND IT.
C
         IF ( NEXTL ) THEN
            IF ( NULINE( I ) .EQ. ';' ) THEN
               NEXTL = .FALSE.
               WAITNL = .TRUE.
C
C  RESET THE CARD TYPE TO 2 WHEN A SPECIAL CARD HAS BEEN FINISHED.
C
               IF ( SPSAL2 ) THEN
                  SPSAL2 = .FALSE.
                  ICARD  = 2
               END IF
            END IF
            GO TO 500
         END IF
C
C  NEXT CHECK WHETHER WE HAVE FOUND AN END OF LINE MARKER ANYWAY.
C
         IF ( NULINE( I ) .EQ. ';' ) THEN
            WAITNL = .TRUE.
C
C  FINISH OFF THE CURRENT LINE.
C
            J       = ISTFIE( NFIELD, ICARD ) - 1
            DO 30 K = 1, LFIELD
               NULINA( NLINES )( J + K: J + K ) =
     *            NULINE( NSTFIE + K )
   30       CONTINUE
            FIELD = .FALSE.
C
C  RESET THE CARD TYPE TO 2 WHEN A SPECIAL CARD HAS BEEN FINISHED.
C
            IF ( SPSAL2 ) THEN
               SPSAL2 = .FALSE.
               ICARD  = 2
            END IF
            GO TO 500
         END IF
C
C  A FIELD HAS BEEN STARTED.
C
         IF ( FIELD ) THEN
C
C  THE FIELD HAS NOW COME TO AN END.
C
            IF ( ( NULINE( I ) .EQ. ' ' .AND. .NOT.
     *              ( ICARD .EQ. 2 .AND. NFIELD .EQ. 4 ) ) .OR.
     *           NULINE( I ) .EQ. '_' ) THEN
               FIELD = .FALSE.
C
C  STORE THE FIELD IN ITS CORRECT POSITION IN NULINA.
C
               J       = ISTFIE( NFIELD, ICARD ) - 1
               DO 40 K = 1, LFIELD
                  NULINA( NLINES )( J + K: J + K ) =
     *               NULINE( NSTFIE + K )
   40          CONTINUE
C
C  THE FIELD HAS NOW COME TO AN END AND A BLANK FIELD FOLLOWS.
C
               IF ( NULINE( I ) .EQ. '_' ) THEN
                  NFIELD = NFIELD + 1
                  LFIELD = 0
                  IF ( NFIELD .GT. NFIE( ICARD ) ) THEN
                     INFORM = 45
                     WRITE( IOUT, 2450 )
                     RETURN
                  END IF
               END IF
            ELSE
C
C  AN EXTRA CHARACTER HAS BEEN ADDED TO THE FIELD.
C
               LFIELD = LFIELD + 1
C
C  CHECK THAT THE FIELD HAS NOT EXCEEDED ITS ALLOWED SPACE.
C  THIS MAY HAPPEN IF A) THERE IS AN ERROR ON THE CARD, B) THE
C  CARD IS OF TYPE 2 AND ONLY BLANKS REMAIN OR C) THE FIELD IS
C  ACTUALLY AN INDICATOR CARD. CHECK WHICH.
C
               IF ( LENFIE( NFIELD, ICARD ) .LT. LFIELD ) THEN
C
C  IF WE ARE CONSIDERING FIELD 4 WHEN THE CARD IS OF TYPE 2 AND
C  ALL THE SPACE HAS BEEN EXHAUSTED, FINISH THE LINE.
C
                  IF ( ICARD .EQ. 2 .AND. NFIELD .EQ. 4 ) THEN
                     WAITNL  = .TRUE.
                     J       = ISTFIE( NFIELD, ICARD ) - 1
                     DO 50 K = 1, LFIELD
                        NULINA( NLINES )( J + K: J + K ) =
     *                     NULINE( NSTFIE + K )
   50                CONTINUE
                     FIELD = .FALSE.
                     GO TO 500
                  END IF
C
C  THERE IS AN ERROR IN THE CURRENT FIELD.
C
                  IF ( NFIELD .GT. 1 ) THEN
                     INFORM = 44
                     WRITE( IOUT, 2440 ) NFIELD
                     RETURN
                  ELSE
C
C  THE FIRST FIELD MAY BE AN INDICATOR CARD. CHECK.
C
                     DO 70 J    = 2, MENDAT
                        LEN     = LENIND( J )
                        DO 60 K = 1, LEN
                           IF ( NULINE( NSTFIE + K ) .NE.
     *                     INDIC8( J )( K: K ) ) GO TO 70
   60                   CONTINUE
                        GO TO 80
   70                CONTINUE
C
C  THE INDICATOR CARD IS UNKNOWN. EXIT WITH AN ERROR MESSAGE.
C
                     INFORM = 2
                     WRITE( IOUT, 2020 )
                     RETURN
C
C  THE INDICATOR CARD IS RECOGNISED. OUTPUT THIS CARD AS THE NEXT
C  LINE AND AWAIT A FURTHER LINE. (THE TITLE CARD IS AN EXCEPTION AS
C  THE TITLE HAS STILL TO ARRIVE)
C
   80                CONTINUE
                     IF ( J .NE. 4 ) THEN
                        FIELD = .FALSE.
                        NEXTL = .TRUE.
                        NULINA( NLINES )( 1: 12 ) =
     *                       INDIC8( J )( 1: 12 )
                     END IF
                  END IF
               END IF
            END IF
C
C  WE ARE BETWEEN FIELDS.
C
         ELSE
C
C  A NEW FIELD HAS STARTED.
C
            IF ( NULINE( I ) .NE. ' ' ) THEN
C
C  IT IS THE FIRST FIELD.
C
               IF ( WAITNL ) THEN
                  WAITNL = .FALSE.
                  NLINES = NLINES + 1
C
C  INITIALIZE THE NEW LINE, NULINA( NLINES ), AS A BLANK LINE.
C
                  DO 90 J = 1, 65
                     NULINA( NLINES )( J: J ) = ' '
   90             CONTINUE
                  NFIELD = 1
                  IF ( NULINE( I ) .EQ. '_' ) LFIELD = 0
C
C  IF A SPECIAL CARD OCCURS (THAT IS, A CARD ON WHICH THE RANGE
C  TRANSFORMATION IS SPECIFIED WITHIN MAKEFN), MARK IT.
C
                  IF ( NULINE( I ) .EQ. 'R' .AND. ICARD .EQ. 2 ) THEN
                     SPSAL2 = .TRUE.
                     ICARD  = 1
                  END IF
               ELSE
C
C  IT IS NOT THE FIRST FIELD.
C
                  NFIELD = NFIELD + 1
                  IF ( NFIELD .GT. NFIE( ICARD ) ) THEN
                     INFORM = 45
                     WRITE( IOUT, 2450 )
                     RETURN
                  END IF
C
C  IF THE STRING IS IN FIELDS 3 OR 5 AND STARTS WITH A '$', THE
C  REMAINDER OF THE CARD IS CONSIDERED TO BE A COMMENT.
C
                  IF ( ( NFIELD .EQ. 3 .OR. NFIELD .EQ. 5 ) .AND.
     *                 NULINE( I ) .EQ. '$' ) THEN
                     J = ISTFIE( NFIELD, ICARD ) - 1
                     DO 100 K = 1, 66 - I
                        NULINA( NLINES )( J + K: J + K ) =
     *                     NULINE( I + K - 1 )
  100                CONTINUE
                     GO TO 600
                  END IF
               END IF
C
C  SKIP A FIELD IF A '_' IS ENCOUNTERED.
C
               IF ( NULINE( I ) .EQ. '_' ) THEN
                  LFIELD = 0
               ELSE
C
C  SET THE CURRENT LENGTH OF THE FIELD, LFIELD, THE STARTING ADDRESS
C  IN NULINE - 1, NSTFIE, AND THE FIELD NUMBER, NFIELD.
C
                  FIELD = .TRUE.
                  LFIELD = 1
                  NSTFIE = I - 1
               END IF
            END IF
         END IF
  500 CONTINUE
  600 CONTINUE
C
C  FINISH OFF THE LAST LINE.
C
      IF ( FIELD ) THEN
         J        = ISTFIE( NFIELD, ICARD ) - 1
         DO 610 K = 1, LFIELD
            NULINA( NLINES )( J + K: J + K ) =
     *         NULINE( NSTFIE + K )
  610    CONTINUE
      END IF
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2020 FORMAT( ' ** Exit from GPSMPS - indicator card not recognised ' )
 2440 FORMAT( ' ** Exit from GPSMPS - field ', I1, ' on free-form',
     *        ' card too long' )
 2450 FORMAT( ' ** Exit from GPSMPS - too many fields on free-form',
     *        ' card' )
C
C  END OF FREEFM.
C
      END
C  THIS VERSION: 28/04/1992 AT 10:12:34 AM.
      INTEGER FUNCTION ONLY1( NUM )
      INTEGER                 NUM
C
C  RETURNS THE VALUE 1 IF NUM IS ONE AND 2 OTHERWISE
C
C  NICK GOULD, FOR CGT PRODUCTIONS.
C  DECEMBER 7TH, 1990.
C
      ONLY1 = 2
      IF ( NUM .EQ. 1 ) ONLY1 = 1
      RETURN
C
C  END OF ONLY1.
C
      END
