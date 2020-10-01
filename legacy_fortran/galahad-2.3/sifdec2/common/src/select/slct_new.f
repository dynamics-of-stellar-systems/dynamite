C     ( Last modified on 13 Aug 2005 at 23:44:33 )
C-----------------------------------------------------------------------------
C
      PROGRAM SELECT
C
C-----------------------------------------------------------------------------
C
C  The purpose of this program is to interrogate the file containing problem
C  classifications (by default, $MASTSIF/CLASSF.DB) for obtaining a list
C  of problems matching interactively defined characteristics.
C
C  The dialog with the user is on the standard input/output.
C
C  Programming: I. Bongartz, A.R. Conn and Ph. Toint for CGT Productions.
C               Salford fortran by Kristjan Jonasson
C               Revision: Ph. Toint, Aug 2005.
C
C--------- THE FOLLOWING SPECIFICATIONS MAY BE MODIFIED BY THE USER ----------
C
C  Standard input default definition is device 5.  Standard output default
C  definition is device 6. Change them to whatever values are appropriate
C  on your system.
C
      INTEGER       STDIN,     STDOUT
      PARAMETER   ( STDIN = 5, STDOUT = 6 )
C
C  Name of the classification database
C
      CHARACTER*32  NAME
      PARAMETER   ( NAME = 'CLASSF.DB' )
C
C  Device number for reading the classification database
C
      INTEGER       CLSDVC, FLSDVC
      PARAMETER   ( CLSDVC = 55, FLSDVC = 56 )
C
C    default directory for classification file.
C
C  Device number and file name for file containing default directory
C
      INTEGER       DATDVC
      PARAMETER   ( DATDVC = 57 )
      CHARACTER*32  DATNAM
      PARAMETER   ( DATNAM = 'SLCT.DAT' )
C
C---------------- END OF THE USER MODIFIABLE SPECIFICATION ------------------
C
C  Classification constants
C 
      INTEGER       OBJ, CON, REG, DER, INTRST, INTVAR, VARN, VARM,
     *              VFR, V1S, V2S, VFX, C1SL, C2SL, CEQL, C1SG, C2SG,
     *              CEQG, UMAX
      PARAMETER   ( OBJ    =  1, CON    =  2, REG  =  3, DER  = 4, 
     *              INTRST =  5, INTVAR =  6, VARN =  7, VARM = 8,
     *              VFR    =  9, V1S    = 10, V2S  = 11, VFX  = 12, 
     *              C1SL   = 13, C2SL   = 14, CEQL = 15, C1SG = 16, 
     *              C2SG   = 17, CEQG   = 18, UMAX = 999999999 )
C
C  Maximum number of simultaneous targets in search
C
      INTEGER       MAXTRG
      PARAMETER   ( MAXTRG = 19 )
C
C  Variable definitions
C
      CHARACTER*80  LINE
C    names for classification file
      CHARACTER*72  FILEN
C    default directory for classification file
      CHARACTER*256 DFTDIR
C    names for output listing file
C                  Addition by Kristjan Jonasson
      CHARACTER*72  FILES
      CHARACTER*200 PBCLS
      CHARACTER*18  TARGET( MAXTRG )
      CHARACTER*8   LIST(5)
      CHARACTER*1   CHOICE, CHAR, UPPER
      INTEGER       I, IM1, J, L, K, NUM, NMATCH, CONVERT, NBT,  NBI, 
     *              SIZE, LOW(VARN:CEQG), UPP(VARN:CEQG)
      INTEGER       CHOSEN( MAXTRG, VARN:CEQG )
      LOGICAL       REJECT, ANYFIX(VARN:CEQG), MATCH
      INTRINSIC     MIN
C
C  Banner
C
      WRITE ( STDOUT, 1000 )
C
C    in order that this program can give the full path name for the
C    default classification file.
C
C  Open the file containing the default directory name.
C
      OPEN ( UNIT = DATDVC, FILE = DATNAM, STATUS = 'OLD' )
C
C  Read in the default directory
C
      READ ( DATDVC, 8000 ) DFTDIR
      CLOSE ( DATDVC )
      DO 910 I = 1, 256
         IF ( DFTDIR( I : I ) .EQ. ' ' ) THEN
            SIZE = I - 1
            GO TO 920
         END IF
  910 CONTINUE
      SIZE = 256
  920 CONTINUE
      FILEN = DFTDIR( 1 : SIZE ) // '/' // NAME
      SIZE  = LEN( FILEN )
C
C  MAIN LOOP
C
  1   CONTINUE
C
C  Target initialization
C
      DO 4 I = 1, MAXTRG
        DO 3 J = VARN, CEQG
          CHOSEN( I, J ) = -2
 3      CONTINUE
        TARGET(I) = 'XXXXXXXXXXXXXXXXXX'
  4   CONTINUE
      TARGET(1) =   '******************'
      DO 21 I = VFR, VARM
        ANYFIX( I ) = .FALSE.
 21   CONTINUE
C
C  Bounds on the number of variables and constraints are initialized
C  to be inactive.
C
      DO 22 I = VFR, VARM
        LOW( I ) = 0
        UPP( I ) = UMAX
 22   CONTINUE
C
C  Verify the name of the classification database file
C
 77   CONTINUE
      WRITE ( STDOUT, 4950 ) FILEN(1:SIZE)
      WRITE ( STDOUT, 4957 )
      READ  ( STDIN, '( A1 )' ) CHOICE
      IF ( UPPER( CHOICE ) .EQ. 'Y' ) THEN
        WRITE ( STDOUT, 4955 )
        READ  ( STDIN , FMT = '( A )', ERR = 78 ) FILEN
        SIZE = LEN( FILEN )
        GO TO 77
 78     WRITE ( STDOUT, 1201 )
        GO TO 77
      ENDIF
C
C  Writes the current specification
C
 5    CONTINUE
      WRITE ( STDOUT, 5002 )
      NUM = NBT ( OBJ, TARGET, MAXTRG )
      WRITE ( STDOUT, 5003 ) ( TARGET(L)(OBJ:OBJ), L = 1, NUM )
      NUM = NBT ( CON, TARGET, MAXTRG )
      WRITE ( STDOUT, 5004 ) ( TARGET(L)(CON:CON), L = 1, NUM )
      WRITE ( STDOUT, 5005 ) TARGET(1)(REG:REG)
      NUM = NBT ( DER, TARGET, MAXTRG )
      WRITE ( STDOUT, 5006 ) ( TARGET(L)(DER:DER), L = 1, NUM )
      NUM = NBT ( INTRST, TARGET, MAXTRG )
      WRITE ( STDOUT, 5007 ) ( TARGET(L)(INTRST:INTRST), L = 1, NUM )
      WRITE ( STDOUT, 5008 ) TARGET(1)(INTVAR:INTVAR)
C
C     Number of variables
C
      WRITE ( STDOUT, 5009 )
      CHAR = TARGET(1)(VFR:VFR)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5020 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5021 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5022 ) LOW( VFR ), UPP( VFR )
      ELSE
        IF ( ANYFIX(VFR) ) THEN
          WRITE ( STDOUT, 5023 )
        ELSE
          NUM = NBI ( CHOSEN( 1, VFR ), MAXTRG )
          WRITE ( STDOUT, 5024 ) ( CHOSEN( L, VFR ), L = 1, NUM )
        ENDIF
      ENDIF
C
      CHAR = TARGET(1)(V1S:V1S)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5025 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5026 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5027 ) LOW( V1S), UPP( V1S )
      ELSE
        IF ( ANYFIX(V1S) ) THEN
          WRITE ( STDOUT, 5028 )
        ELSE
          NUM = NBI ( CHOSEN( 1, V1S ), MAXTRG )
          WRITE ( STDOUT, 5029 ) ( CHOSEN( L, V1S ), L = 1, NUM )
        ENDIF
      ENDIF
C
      CHAR = TARGET(1)(V2S:V2S)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5030 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5031 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5032 ) LOW( V2S ), UPP( V2S )
      ELSE
        IF ( ANYFIX(V2S) ) THEN
          WRITE ( STDOUT, 5033 )
        ELSE
          NUM = NBI ( CHOSEN( 1, V2S ), MAXTRG )
          WRITE ( STDOUT, 5034 ) ( CHOSEN( L, V2S ), L = 1, NUM )
        ENDIF
      ENDIF
C
      CHAR = TARGET(1)(VFX:VFX)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5035 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5036 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5037 ) LOW( VFX ), UPP( VFX )
      ELSE
        IF ( ANYFIX(VFX) ) THEN
          WRITE ( STDOUT, 5038 )
        ELSE
          NUM = NBI ( CHOSEN( 1, VFX ), MAXTRG )
          WRITE ( STDOUT, 5039 ) ( CHOSEN( L, VFX ), L = 1, NUM )
        ENDIF
      ENDIF
C
      CHAR = TARGET(1)(VARN:VARN)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5040 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5041 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5042 ) LOW( VARN ), UPP( VARN )
      ELSE
        IF ( ANYFIX(VARN) ) THEN
          WRITE ( STDOUT, 5043 )
        ELSE
          NUM = NBI ( CHOSEN( 1, VARN ), MAXTRG )
          WRITE ( STDOUT, 5044 ) ( CHOSEN( L, VARN ), L = 1, NUM )
        ENDIF
      ENDIF
C
C  Numbers of linear constraints
C
      WRITE ( STDOUT, 5050 )
      CHAR = TARGET(1)(C1SL:C1SL)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5025 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5026 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5027 ) LOW( C1SL ), UPP( C1SL )
      ELSE
        IF ( ANYFIX(C1SL) ) THEN
          WRITE ( STDOUT, 5028 )
        ELSE
          NUM = NBI ( CHOSEN( 1, C1SL ), MAXTRG )
          WRITE ( STDOUT, 5029 ) ( CHOSEN( L, C1SL ), L = 1, NUM )
        ENDIF
      ENDIF
C
      CHAR = TARGET(1)(C2SL:C2SL)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5030 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5031 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5032 ) LOW( C2SL ), UPP( C2SL )
      ELSE
        IF ( ANYFIX(C2SL) ) THEN
          WRITE ( STDOUT, 5033 )
        ELSE
          NUM = NBI( CHOSEN( 1, C2SL ), MAXTRG )
          WRITE ( STDOUT, 5034 ) ( CHOSEN( L, C2SL ), L = 1, NUM )
        ENDIF
      ENDIF
C
      CHAR = TARGET(1)(CEQL:CEQL)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5045 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5046 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5047 ) LOW( CEQL ), UPP( CEQL )
      ELSE
        IF ( ANYFIX(CEQL) ) THEN
          WRITE ( STDOUT, 5048 )
        ELSE
          NUM = NBI ( CHOSEN( 1, CEQL ), MAXTRG )
          WRITE ( STDOUT, 5049 ) ( CHOSEN( L, CEQL ), L = 1, NUM )
        ENDIF
      ENDIF
C
C  Numbers of general constraints
C
      WRITE ( STDOUT, 5051 )
      CHAR = TARGET(1)(C1SG:C1SG)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5025 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5026 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5027 ) LOW( C1SG ), UPP( C1SG )
      ELSE
        IF ( ANYFIX(C1SG) ) THEN
          WRITE ( STDOUT, 5028 )
        ELSE
          NUM = NBI ( CHOSEN( 1, C1SG ), MAXTRG )
          WRITE ( STDOUT, 5029 ) ( CHOSEN( L, C1SG ), L = 1, NUM )
        ENDIF
      ENDIF
C
      CHAR = TARGET(1)(C2SG:C2SG)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5030 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5031 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5032 ) LOW( C2SG ), UPP( C2SG )
      ELSE
        IF ( ANYFIX(C2SG) ) THEN
          WRITE ( STDOUT, 5033 )
        ELSE
          NUM = NBI ( CHOSEN( 1, C2SG ), MAXTRG )
          WRITE ( STDOUT, 5034 ) ( CHOSEN( L, C2SG ), L = 1, NUM )
        ENDIF
      ENDIF
C
      CHAR = TARGET(1)(CEQG:CEQG)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5045 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5046 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5047 ) LOW( CEQG ), UPP( CEQG )
      ELSE
        IF ( ANYFIX(CEQG) ) THEN
          WRITE ( STDOUT, 5048 )
        ELSE
          NUM = NBI ( CHOSEN( 1, CEQG ), MAXTRG )
          WRITE ( STDOUT, 5049 ) ( CHOSEN( L, CEQG ), L = 1, NUM )
        ENDIF
      ENDIF
C
      CHAR = TARGET(1)(VARM:VARM)
      IF ( CHAR .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5040 )
      ELSE IF ( CHAR .EQ. '*' ) THEN
        WRITE ( STDOUT, 5041 )
      ELSE IF ( CHAR .EQ. 'I' ) THEN
        WRITE ( STDOUT, 5042 ) LOW( VARM ), UPP( VARM )
      ELSE
        IF ( ANYFIX(VARM)) THEN
          WRITE ( STDOUT, 5043 )
        ELSE
          NUM = NBI ( CHOSEN( 1, VARM ), MAXTRG )
          WRITE ( STDOUT, 5044 ) ( CHOSEN( L, VARM ), L = 1, NUM )
        ENDIF
      ENDIF
C
C  Open the classification file.
C
      OPEN( UNIT = CLSDVC, FILE = FILEN, STATUS = 'OLD' )
C
C  Read in the problem characteristic one wishes to specify
C
      WRITE ( STDOUT, 5000 )
 6    CONTINUE
      WRITE ( STDOUT, 5001 )
      READ  ( STDIN, '( A1 )', ERR = 2 ) CHAR
      CHAR = UPPER( CHAR )
C
C  Get objective function target types
C
      IF ( CHAR .EQ. 'O' ) THEN
 41     CONTINUE
        NUM = 0
        WRITE ( STDOUT, 1002 )
        DO 30 I = 1, MIN( MAXTRG, 6 )
 10       CONTINUE
          IF ( NUM .EQ. 0) WRITE ( STDOUT, 1103 )
          IF ( NUM .GT. 0) WRITE ( STDOUT, 1003 )
          READ  ( STDIN , FMT = '( A1 )', ERR = 20 ) CHOICE
          CHOICE = UPPER( CHOICE )
          IF ( REJECT( CHOICE, OBJ, NUM ) ) GO TO 20
          IF ( I .GT. 1 ) THEN
            IM1 = I - 1
            DO 35 J = 1, IM1
              IF ( TARGET(J)(OBJ:OBJ) .EQ. CHOICE ) THEN
                WRITE ( STDOUT, 1101 )
                GO TO 10
              END IF
 35         CONTINUE
          ENDIF
          TARGET(I)(OBJ:OBJ) = CHOICE
          IF ( CHOICE .NE. ' ' ) NUM = NUM + 1
          IF ( CHOICE .EQ. ' ' .OR. CHOICE .EQ. '*' ) GO TO 40
          GO TO 30
 20       CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 10
 30     CONTINUE
C
C   Verify the type of objective function
C
 40     CONTINUE
        IF ( NUM .EQ. 0 ) THEN
          WRITE ( STDOUT, 2102 )
          GO TO 41
        ELSE
          WRITE ( STDOUT, 2002 ) ( TARGET(K)(OBJ:OBJ), K = 1, NUM )
        ENDIF
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 41
C
C  Get constraints target types
C
      ELSE IF ( CHAR .EQ. 'C' ) THEN
 111    CONTINUE
        NUM = 0
        WRITE ( STDOUT, 1004 )
        DO 130 I = 1, MIN( MAXTRG, 6 )
 110      CONTINUE
          IF ( NUM .EQ. 0) WRITE ( STDOUT, 1105 )
          IF ( NUM .GT. 0) WRITE ( STDOUT, 1005 )
          READ  ( STDIN , FMT = '( A1 )', ERR = 120 ) CHOICE
          CHOICE = UPPER( CHOICE )
          IF ( REJECT( CHOICE, CON, NUM ) ) GO TO 120
          IF ( I .GT. 1 ) THEN
            IM1 = I - 1
            DO 135 J = 1, IM1
              IF ( TARGET(J)(CON:CON) .EQ. CHOICE ) THEN
                WRITE ( STDOUT, 1101 )
                GO TO 110
              END IF
 135        CONTINUE
          ENDIF
          TARGET(I)(CON:CON) = CHOICE
          IF ( CHOICE .NE. ' ' ) NUM = NUM + 1
          IF ( CHOICE .EQ. ' ' .OR. CHOICE .EQ. '*' ) GO TO 140
          GO TO 130
 120      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 110
 130    CONTINUE
C
C   Verify the type of constraints
C
 140    CONTINUE
        IF ( NUM .EQ. 0 ) THEN
          WRITE ( STDOUT, 2104 )
          GO TO 111
        ELSE
          WRITE ( STDOUT, 2004 ) ( TARGET(K)(CON:CON), K = 1, NUM )
        ENDIF
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 111
C
C  Get regularity target types
C
      ELSE IF ( CHAR .EQ. 'R' ) THEN
 211    CONTINUE
        WRITE ( STDOUT, 1006 )
 210    CONTINUE
        WRITE ( STDOUT, 1007 )
        READ  ( STDIN , FMT = '( A1 )', ERR = 220 ) CHOICE
        CHOICE = UPPER( CHOICE )
        IF ( REJECT( CHOICE, REG, 0 ) ) GO TO 220
        TARGET(1)(REG:REG) = CHOICE
        GO TO 240
 220    CONTINUE
        WRITE ( STDOUT, 1001 )
        GO TO 210
 240    CONTINUE
C
C   Verify the problem's regularity
C
        WRITE ( STDOUT, 2005 ) TARGET(1)(REG:REG)
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 211
C
C  Get derivative degree
C
      ELSE IF ( CHAR .EQ. 'D' ) THEN
 311    CONTINUE
        NUM = 0
        WRITE ( STDOUT, 1008 )
        DO 330 I = 1, MIN( MAXTRG, 3 )
 310      CONTINUE
          IF ( NUM .EQ. 0 ) WRITE ( STDOUT, 1109 )
          IF ( NUM .GT. 0 ) WRITE ( STDOUT, 1009 )
          READ  ( STDIN , '( A1 )', ERR = 320 ) CHOICE
          CHOICE = UPPER( CHOICE )
          IF ( REJECT( CHOICE, DER, NUM ) ) GO TO 320
          IF ( I .GT. 1 ) THEN
            IM1 = I - 1
            DO 335 J = 1, IM1
              IF ( TARGET(J)(DER:DER) .EQ. CHOICE ) THEN
                WRITE ( STDOUT, 1101 )
                GO TO 310
              END IF
 335        CONTINUE
          ENDIF
          TARGET(I)(DER:DER) = CHOICE
          IF ( CHOICE .NE. ' ' ) NUM = NUM + 1
          IF ( CHOICE .EQ. ' ' .OR. CHOICE .EQ. '*' ) GO TO 340
          GO TO 330
 320      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 310
 330    CONTINUE
C
C   Verify the degree of available derivatives
C
 340    CONTINUE
        IF ( NUM .EQ. 0 ) THEN
          WRITE ( STDOUT, 2106 )
          GO TO 311
        ELSE
          WRITE ( STDOUT, 2006 ) ( TARGET(K)(DER:DER), K = 1, NUM )
        ENDIF
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 311
C
C  Get problem interest
C
      ELSE IF ( CHAR .EQ. 'I' ) THEN
 411    CONTINUE
        NUM = 0
        WRITE ( STDOUT, 1010 )
        DO 430 I = 1, MIN( MAXTRG, 3 )
 410      CONTINUE
          IF ( NUM .EQ. 0 ) WRITE ( STDOUT, 1111 )
          IF ( NUM .GT. 0 ) WRITE ( STDOUT, 1011 )
          READ  ( STDIN , FMT = '( A1 )', ERR = 420 ) CHOICE
          CHOICE = UPPER( CHOICE )
          IF ( REJECT( CHOICE, INTRST, NUM ) ) GO TO 420
          IF ( I .GT. 1 ) THEN
            IM1 = I - 1
            DO 435 J = 1, IM1
              IF ( TARGET(J)(INTRST:INTRST) .EQ. CHOICE ) THEN
                WRITE ( STDOUT, 1101 )
                GO TO 410
              END IF
 435        CONTINUE
          ENDIF
          TARGET(I)(INTRST:INTRST) = CHOICE
          IF ( CHOICE .NE. ' ' ) NUM = NUM + 1
          IF ( CHOICE .EQ. ' ' .OR. CHOICE .EQ. '*' ) GO TO 440
          GO TO 430
 420      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 410
 430    CONTINUE
C
C   Verify the problem's interest
C
 440    CONTINUE
        IF ( NUM .EQ. 0 ) THEN
          WRITE ( STDOUT, 2107 )
          GO TO 411
        ELSE
          WRITE ( STDOUT, 2007 ) ( TARGET(K)(INTRST:INTRST), K = 1, NUM)
        ENDIF
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 411
C
C  Get the internal variables indicator
C
      ELSE IF ( CHAR .EQ. 'S' ) THEN
 511    CONTINUE
        WRITE ( STDOUT, 1012 )
 510    CONTINUE
        WRITE ( STDOUT, 1013 )
        READ  ( STDIN , '( A1 )', ERR = 520 ) CHOICE
        CHOICE = UPPER( CHOICE )
        IF ( REJECT( CHOICE, INTVAR, 0 ) ) GO TO 520
        TARGET(1)(INTVAR:INTVAR) = CHOICE
        GO TO 540
 520    CONTINUE
        WRITE ( STDOUT, 1001 )
        GO TO 510
C
C   Verify the indicator for explicit internal variables
C
 540    CONTINUE
        IF ( CHOICE .EQ. 'Y' ) THEN
          WRITE ( STDOUT, 2008 )
        ELSE IF ( CHOICE .EQ. 'N' ) THEN
          WRITE ( STDOUT, 2018 )
        ELSE
          WRITE ( STDOUT, 2017 )
        ENDIF
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 511
C
C  Get the number of variables
C
      ELSE IF ( CHAR .EQ. 'F' ) THEN

        CALL GETDIM( VFR, TARGET(1)(VFR:VFR), ANYFIX(VFR), 
     1               CHOSEN(1,VFR), LOW(VFR), UPP(VFR), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'free variables+' ) 
   
      ELSE IF ( CHAR .EQ. 'A' ) THEN

        CALL GETDIM( V1S, TARGET(1)(V1S:V1S), ANYFIX(V1S), 
     1               CHOSEN(1,V1S), LOW(V1S), UPP(V1S), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'variables bounded below or above+' ) 
   
      ELSE IF ( CHAR .EQ. 'B' ) THEN
        CALL GETDIM( V2S, TARGET(1)(V2S:V2S), ANYFIX(V2S), 
     1               CHOSEN(1,V2S), LOW(V2S), UPP(V2S), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'variables bounded below and above+' ) 
   
      ELSE IF ( CHAR .EQ. 'X' ) THEN

        CALL GETDIM( VFX, TARGET(1)(VFX:VFX), ANYFIX(VFX), 
     1               CHOSEN(1,VFX), LOW(VFX), UPP(VFX), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'fixed variables+' ) 
   
      ELSE IF ( CHAR .EQ. 'N' ) THEN

        CALL GETDIM( VARN, TARGET(1)(VARN:VARN), ANYFIX(VARN), 
     1               CHOSEN(1,VARN), LOW(VARN), UPP(VARN), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'variables+' ) 
C
C  Get the number of constraints
C
      ELSE IF ( CHAR .EQ. 'P' ) THEN
        CALL GETDIM( C1SL, TARGET(1)(C1SL:C1SL), ANYFIX(C1SL), 
     1               CHOSEN(1,C1SL), LOW(C1SL), UPP(C1SL), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'one-sided linear constraints+' ) 

      ELSE IF ( CHAR .EQ. 'Q' ) THEN
        CALL GETDIM( C2SL, TARGET(1)(C2SL:C2SL), ANYFIX(C2SL), 
     1               CHOSEN(1,C2SL), LOW(C2SL), UPP(C2SL), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'two-sided linear constraints+' ) 

      ELSE IF ( CHAR .EQ. 'E' ) THEN
        CALL GETDIM( CEQL, TARGET(1)(CEQL:CEQL), ANYFIX(CEQL), 
     1               CHOSEN(1,CEQL), LOW(CEQL), UPP(CEQL), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'linear equality constraints+' ) 

      ELSE IF ( CHAR .EQ. 'V' ) THEN
        CALL GETDIM( C1SG, TARGET(1)(C1SG:C1SG), ANYFIX(C1SG), 
     1               CHOSEN(1,C1SG), LOW(C1SG), UPP(C1SG), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'one-sided general constraints+' )

      ELSE IF ( CHAR .EQ. 'W' ) THEN
        CALL GETDIM( C2SG, TARGET(1)(C2SG:C2SG), ANYFIX(C2SG), 
     1               CHOSEN(1,C2SG), LOW(C2SG), UPP(C2SG), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'two-sided general constraints+' ) 

      ELSE IF ( CHAR .EQ. 'G' ) THEN
        CALL GETDIM( CEQG, TARGET(1)(CEQG:CEQG), ANYFIX(CEQG), 
     1               CHOSEN(1,CEQG), LOW(CEQG), UPP(CEQG), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'general equality constraints+' ) 

      ELSE IF ( CHAR .EQ. 'M' ) THEN
        CALL GETDIM( VARM, TARGET(1)(VARM:VARM), ANYFIX(VARM), 
     1               CHOSEN(1,VARM), LOW(VARM), UPP(VARM), MAXTRG, 
     1               STDIN, STDOUT, UMAX, 
     1               'constraints+' ) 
C
C  All characteristics have been recorded. Quit
C
      ELSE IF ( CHAR .EQ. ' ' ) THEN
        GO TO 7
C
C  Error in the choice of characteristic
C
      ELSE
        GO TO 2
      ENDIF
C
C  Loop for another characteristic
C
      CLOSE ( CLSDVC )
      GO TO 5
C
C  Handle the error in characteristic choice
C
 2    CONTINUE
      WRITE ( STDOUT, 1001 )
      GO TO 6
C
C  The target(s) are now defined.
C
 7    CONTINUE
      WRITE (STDOUT, 4001 )
C
C  Loop on the problem classification and print names of matching ones
C  in lots of 5.  Also accumulate the number of matching problems found.
C
      NMATCH = 0
      L      = 0
 3000 CONTINUE
      READ ( CLSDVC, '( A200 )', END = 3010 ) PBCLS
      IF ( PBCLS(1:1) .NE. '*' ) THEN
        IF ( MATCH( PBCLS, TARGET, MAXTRG, ANYFIX, CHOSEN,
     1              LOW, UPP) ) THEN
          NMATCH = NMATCH + 1
          L = L + 1
          LIST(L) = PBCLS(1:8)
          IF ( L .EQ. 5 ) THEN
            WRITE ( STDOUT, 4002 ) ( LIST(I), I = 1, 5 )
            L = 0
          END IF
        END IF
      ENDIF
      GO TO 3000
C
C  End of the database processing for the main loop.
C
C  Print selected problem names still
C  in waiting list for output.
C
 3010 CONTINUE
      IF ( L .GE. 1 ) WRITE ( STDOUT, 4002 ) ( LIST(I), I = 1, L )
      WRITE ( STDOUT, 4000 ) NMATCH
      CLOSE ( CLSDVC )
      WRITE ( STDOUT, 7001 )
      READ  ( STDIN, '( A1 )' ) CHOICE
      IF ( UPPER( CHOICE ) .EQ. 'Y' ) THEN
        OPEN( UNIT = CLSDVC, FILE = FILEN, STATUS = 'OLD' )
        WRITE ( STDOUT, 4955 )
        READ  ( STDIN , FMT = '( A )' ) FILES
        SIZE = LEN( FILES )
        OPEN( UNIT = FLSDVC, FILE = FILES, STATUS = 'UNKNOWN' )
 3020   CONTINUE
        READ ( CLSDVC, '( A36 )', END = 3050 ) PBCLS
        IF ( PBCLS(1:1) .NE. '*' ) THEN
          IF ( MATCH( PBCLS, TARGET, MAXTRG, ANYFIX, CHOSEN,
     1                LOW, UPP ) ) THEN
            NBLK = 8
            DO 3030 NLBK = 8, 2, -1
              IF ( PBCLS(NBLK:NBLK) .NE. ' ' ) GO TO 3040
 3030       CONTINUE
 3040       CONTINUE
            WRITE ( FLSDVC, '(A)' ) PBCLS(1:NBLK)
          END IF
        END IF
        GO TO 3020
 3050   CLOSE ( FLSDVC )
        CLOSE ( CLSDVC )
      ENDIF
      WRITE ( STDOUT, 7000 )
      READ  ( STDIN, '( A1 )' ) CHOICE
      IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 1
C
C  End of the database processing.
C
      STOP
C
C  Non excutable statements
C
 1000 FORMAT( /'      *************************************************'
     1        /'      *                                               *'
     1        /'      *         CONSTRAINED AND UNCONSTRAINED         *'
     1        /'      *              TESTING ENVIRONMENT              *'
     1        /'      *                                               *'
     1        /'      *                   ( CUTE )                    *'
     1        /'      *                                               *'
     1        /'      *         INTERACTIVE PROBLEM SELECTION         *'
     1        /'      *                                               *'
     1        /'      *                CGT PRODUCTIONS                *'
     1        /'      *                                               *'
     1        /'      *************************************************'
     1        / )
 1101 FORMAT( '   *** THIS SELECTION IS REPETITION.  Please choose ',
     1             'again.'
     1        / )
 1201 FORMAT( '   *** YOU USED MORE THAN 32 CHARACTERS.  Please choose',
     1             ' again.'
     1        / )
 1001 FORMAT( '   *** YOUR ANSWER IS NOT ALLOWED.  Please choose again.'
     1        / )
 4950 FORMAT( /'   Your current classification file is : ',
     1                A )
 4957 FORMAT( /'   Do you wish to change this [<CR> = N] ?',
     1             ' (N/Y)')
 4955 FORMAT(  '   Input the filename you want (up to 32 characters): ')
 5000 FORMAT( /'   CHOOSE A PROBLEM CHARACTERISTIC THAT YOU WANT',
     1         ' TO SPECIFY :'
     1        /'   ---------------------------------------------',
     1         '-------------' )
 5001 FORMAT( /'     O   : Objective type          C : Constraint',
     1                   ' type'
     1        /'     R   : Regularity              I : Problem',
     1                   ' interest'
     1        /'     D   : Degree of available analytic derivatives'
     1        /'     S   : Presence of explicit internal variables'
     1        /'     F   : Number of free variables'
     1        /'     A   : Number of variables bounded below OR above'
     1        /'     B   : Number of variables bounded below AND above'
     1        /'     X   : Number of fixed variables'
     1        /'     N   : Total number of variables'
     1        /'     P   : Number of 1-sided linear inequality ',
     1                     'constraints'
     1        /'     Q   : Number of 2-sided linear inequality ',
     1                     'constraints'
     1        /'     E   : Number of linear equality constraints'
     1        /'     V   : Number of 1-sided general inequality ',
     1                     'constraints'
     1        /'     W   : Number of 2-sided general inequality ',
     1                     'constraints'
     1        /'     G   : Number of general equality constraints'
     1        /'     M   : Total number of constraints'
     1        /'    <CR> : No further characteristic, perform',
     1               ' selection'
     1        /' '
     1        /'   Your choice :' )
 5002 FORMAT( /'   Your current problem selection key is:'
     1        /'     ( * = anything goes )' )
 5003 FORMAT( /'     Objective function type         : ',
     1               10 ( A1, 1X ) )
 5004 FORMAT(  '     Constraints type                : ',
     1               10 ( A1, 1X ) )
 5005 FORMAT(  '     Regularity                      : ',
     1               10 ( A1, 1X ) )
 5006 FORMAT(  '     Degree of available derivatives : ',
     1               10 ( A1, 1X ) )
 5007 FORMAT(  '     Problem interest                : ',
     1               10 ( A1, 1X ) )
 5008 FORMAT(  '     Explicit internal variables     : ',
     1               10 ( A1, 1X ) )
 5009 FORMAT(  '     Number of variables' )
 5020 FORMAT(  '          Free                       : v' )
 5021 FORMAT(  '          Free                       : *' )
 5022 FORMAT(  '          Free                       : in [ ',
     1        I9,', ',I9,' ]' )
 5023 FORMAT(  '          Free                       : ',
     1         'any fixed number ')
 5024 FORMAT(  '          Free                       :', 3( 1X, I9 ),
     1        /'                                      ', 3( 1X, I9 ) )
 5025 FORMAT(  '          Bounded below  OR above    : v' )
 5026 FORMAT(  '          Bounded below  OR above    : *' )
 5027 FORMAT(  '          Bounded below  OR above    : in [ ',
     1        I9,', ',I9,' ]' )
 5028 FORMAT(  '          Bounded below  OR above    : ',
     1         'any fixed number ')
 5029 FORMAT(  '          Bounded below  OR above    :', 3( 1X, I9 ),
     1        /'                                      ', 3( 1X, I9 ) )
 5030 FORMAT(  '          Bounded below AND above    : v' )
 5031 FORMAT(  '          Bounded below AND above    : *' )
 5032 FORMAT(  '          Bounded below AND above    : in [ ',
     1        I9,', ',I9,' ]' )
 5033 FORMAT(  '          Bounded below AND above    : ',
     1         'any fixed number ')
 5034 FORMAT(  '          Bounded below AND above    :', 3( 1X, I9 ),
     1        /'                                      ', 3( 1X, I9 ) )
 5035 FORMAT(  '          Fixed                      : v' )
 5036 FORMAT(  '          Fixed                      : *' )
 5037 FORMAT(  '          Fixed                      : in [ ',
     1        I9,', ',I9,' ]' )
 5038 FORMAT(  '          Fixed                      : ',
     1         'any fixed number ')
 5039 FORMAT(  '          Fixed                      :', 3( 1X, I9 ),
     1        /'                                      ', 3( 1X, I9 ) )
 5040 FORMAT(  '        Total                        : v' )
 5041 FORMAT(  '        Total                        : *' )
 5042 FORMAT(  '        Total                        : in [ ',
     1        I9,', ',I9,' ]' )
 5043 FORMAT(  '        Total                        : ',
     1         'any fixed number ')
 5044 FORMAT(  '        Total                        :', 3( 1X, I9 ),
     1        /'                                      ', 3( 1X, I9 ) )
 5045 FORMAT(  '          Equalities                 : v' )
 5046 FORMAT(  '          Equalities                 : *' )
 5047 FORMAT(  '          Equalities                 : in [ ',
     1        I9,', ',I9,' ]' )
 5048 FORMAT(  '          Equalities                 : ',
     1         'any fixed number ')
 5049 FORMAT(  '          Equalities                 :', 3( 1X, I9 ),
     1        /'                                      ', 3( 1X, I9 ) )
 5050 FORMAT(  '     Number of constraints ',/,'        Linear:' )
 5051 FORMAT(  '        General:' )
 1002 FORMAT( /'   OBJECTIVE FUNCTION TYPE :'
     1        /'   -------------------------' )
 1103 FORMAT( /'     C   : Constant                L : Linear'
     1        /'     Q   : Quadratic               S : Sum of squares'
     1        /'     N   : No objective'
     1        /'     O   : Other (that is none of the above)'
     1        /'    <CR> : Any of the above (*)'
     1        /' '
     1        /'   Your choice :' )
 1003 FORMAT( /'     C   : Constant                L : Linear'
     1        /'     Q   : Quadratic               S : Sum of squares'
     1        /'     N   : No objective'
     1        /'     O   : Other (that is none of the above)'
     1        /'    <CR> : No further type'
     1        /' '
     1        /'   Your choice :' )
 1004 FORMAT( /'   CONSTRAINTS TYPE :'
     1        /'   ------------------' )
 1105 FORMAT( /'     U   : No constraint           X : Fixed variables',
     1                                                 ' only'
     1        /'     B   : Bounds only             N : Linear network'
     1        /'     L   : Linear                  Q : Quadratic'
     1        /'     O   : Other (that is more general than any of the',
     1                     ' above alone)'
     1        /'    <CR> : Any of the above (*)'
     1        /' '
     1        /'   Your choice :' )
 1005 FORMAT( /'     U   : No constraint           X : Fixed variables',
     1                                                 ' only'
     1        /'     B   : Bounds only             N : Linear network'
     1        /'     L   : Linear                  Q : Quadratic'
     1        /'     O   : Other (that is more general than any of the',
     1                     ' above alone)'
     1        /'    <CR> : No further type'
     1        /' '
     1        /'   Your choice :' )
 1006 FORMAT( /'   PROBLEM REGULARITY TYPE :'
     1        /'   -------------------------' )
 1007 FORMAT( /'     R   : Twice continuously differentiable'
     1        /'     I   : Other'
     1        /'    <CR> : Any of the above (*)'
     1        /' '
     1        /'   Your choice :' )
 1008 FORMAT( /'   DEGREE OF AVAILABLE ANALYTICAL DERIVATIVES :'
     1        /'   --------------------------------------------' )
 1109 FORMAT( /'     0   : No analytical der.      1 : Analytical first'
     1        /'     2   : Analytical second'
     1        /'    <CR> : Any of the above (*)'
     1        /' '
     1        /'   Your choice :' )
 1009 FORMAT( /'     0   : No analytical der.      1 : Analytical first'
     1        /'     2   : Analytical second'
     1        /'    <CR> : No further degree'
     1        /' '
     1        /'   Your choice :' )
 1010 FORMAT( /'   PROBLEM INTEREST TYPE :'
     1        /'   -----------------------' )
 1111 FORMAT( /'     A   : Academic                R : Real application'
     1        /'     M   : Modelling'
     1        /'    <CR> : Any of the above (*)'
     1        /' '
     1        /'   Your choice :' )
 1011 FORMAT( /'     A   : Academic                R : Real application'
     1        /'     M   : Modelling'
     1        /'    <CR> : No further type'
     1        /' '
     1        /'   Your choice :' )
 1012 FORMAT( /'   PRESENCE OF EXPLICIT INTERNAL VARIABLES :'
     1        /'   -----------------------------------------' )
 1013 FORMAT( /'     Y   : Yes                     N : No'
     1        /'    <CR> : Any of the above (*)'
     1        /' '
     1        /'   Your choice :' )
 2002 FORMAT( /'   You have specified objective of type(s): ',
     1             10( A1, 1X ) )
 2102 FORMAT( /'   *** Please choose an objective function type.' )
 2003 FORMAT(  '   Do you wish to reconsider your choice [<CR> = N] ?',
     1             ' (N/Y)')
 2004 FORMAT( /'   You have specified constraints of type(s): ',
     1             10( A1, 1X ) )
 2104 FORMAT( /'   *** Please choose a constraint type.' )
 2005 FORMAT( /'   You have specified regularity of type: ', A1 )
 2006 FORMAT( /'   You have specified derivatives of degree: ',
     1             10( A1, 1X ) )
 2106 FORMAT(/'   *** Please choose a degree of available derivatives.')
 2007 FORMAT( /'   You have specified problem interest of type: ',
     1             10( A1, 1X ) )
 2107 FORMAT( /'   *** Please choose a problem interest type.' )
 2008 FORMAT( /'   You have specified problems with explicit',
     1         ' internal variables.' )
 2018 FORMAT( /'   You have specified problems without explicit',
     1         ' internal variables.' )
 2017 FORMAT( /'   You have specified problems with or without',
     1         ' explicit internal variables.' )
 4000 FORMAT( /'   ', I5, ' Problem(s) match(es) the specification.' / )
 4001 FORMAT( /'   MATCHING PROBLEMS :'
     1        /'   -------------------' / )
 4002 FORMAT(  '     ', 5( A8, 3X ) )
C                            Added by Kristjan Jonasson.
 7001 FORMAT( /'   Do you wish to save the problem names to a file',
     1             ' [<CR> = N] ? (N/Y)')
 7000 FORMAT( /'   Do you wish to make another selection [<CR> = N] ?',
     1             ' (N/Y)')
 8000 FORMAT( A256 )
      END
C
C-----------------------------------------------------------------------------
C
      SUBROUTINE GETDIM( IDX, T, ANYIDX, NIDX, LOW, UPP, MAXTRG, 
     1                   STDIN, STDOUT, UMAX, THING ) 
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to obtain the target for one particular
C  problem dimension.
C
C  Programming: Ph. Toint, Aug 2005.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      INTEGER       IDX, MAXTRG, STDOUT, STDIN, UMAX, NIDX( MAXTRG ),
     1              LOW, UPP
      LOGICAL       ANYIDX
      CHARACTER*1   T
      CHARACTER*30  THING
C
C  Other variables
C
      INTEGER       NUM, I, J, IM1, K, LTH, CONVERT
      LOGICAL       REJECT
      CHARACTER*1   CHOICE, UPPER
      CHARACTER*30  UTHING, UNDER
      CHARACTER*80  LINE
C
      DO 100 I = 1, 30
        IF ( THING(I:I) .NE. '+' ) THEN
           LTH = I
           UTHING(I:I) = UPPER( THING(I:I) )
           UNDER(I:I)  = '-'
        ELSE
           GO TO 611
        END IF
 100  CONTINUE
C
C  Get the number of objects.
C
 611  CONTINUE
      WRITE( STDOUT, 1014 ) UTHING(1:LTH), UNDER(1:LTH)
      NUM = 0
 610  CONTINUE
      WRITE ( STDOUT, 1015 ) THING(1:LTH)
      READ  ( STDIN , FMT = '( A1 )', ERR = 620 ) CHOICE
      CHOICE = UPPER( CHOICE )
      IF ( REJECT( CHOICE, IDX, 0 ) ) GO TO 620
      T = CHOICE
      GO TO 640
 620  CONTINUE
      WRITE ( STDOUT, 1001 )
      GO TO 610
C
C  Verify the number type
C
 640  CONTINUE
      IF ( CHOICE .EQ. 'F' ) THEN
        WRITE ( STDOUT, 2009 ) THING(1:LTH)
      ELSE IF ( CHOICE .EQ. 'V' ) THEN
        WRITE ( STDOUT, 2019 ) THING(1:LTH)
      ELSE IF ( CHOICE .EQ. 'I' ) THEN
        WRITE ( STDOUT, 4009 ) THING(1:LTH)
      ELSE
        WRITE ( STDOUT, 2029 ) THING(1:LTH)
      ENDIF
      WRITE ( STDOUT, 2003 )
      READ  ( STDIN, '( A1 )' ) CHOICE
      IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 611
C
C  Get the fixed number
C
      IF ( T .EQ. 'F' ) THEN
 641    CONTINUE
        NUM = 0
        DO 670 I = 1, MAXTRG
          WRITE( STDOUT, 1018 ) UTHING(1:LTH), UNDER(1:LTH)
 650      CONTINUE
          IF ( NUM .EQ. 0) 
     1       WRITE( STDOUT, 1119 ) THING(1:LTH), THING(1:LTH)
          IF ( NUM .GT. 0) 
     1       WRITE( STDOUT, 1019 ) THING(1:LTH), THING(1:LTH)
          READ ( STDIN , '( A80 )', ERR = 660 ) LINE
          CHOICE = UPPER( LINE(1:1) )
          IF ( NUM .EQ. 0 ) 
     1       ANYIDX = CHOICE .EQ. '*' .OR. CHOICE .EQ. ' '
          IF ( NUM .NE. 0 ) ANYIDX  = CHOICE .EQ. '*'
          IF ( CHOICE .EQ. ' ' .OR. ANYIDX ) GO TO 695
          NIDX( I ) = CONVERT( LINE(1:9),9 )
          IF ( I .GT. 1 ) THEN
            IM1 = I - 1
            DO 655 J = 1, IM1
              IF ( NIDX( J ) .EQ. NIDX( I ) ) GO TO 665
 655        CONTINUE
          ENDIF
          IF ( NIDX( I ) .LT. 0 .OR. NIDX( I ) .GT. UMAX ) GO TO 660
          NUM = NUM + 1
          GO TO 670
 660      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 650
 665      CONTINUE
          WRITE ( STDOUT, 1101 )
          GO TO 650
 670    CONTINUE
C
C   Verify the number
C
 695    CONTINUE
        IF ( ANYIDX ) THEN
          WRITE ( STDOUT, 2020 ) THING(1:LTH)
        ELSE
          IF ( NUM .EQ. 0 ) THEN
            WRITE ( STDOUT, 2110 ) THING(1:LTH)
            GO TO 641
          ELSE
            WRITE ( STDOUT, 2010 ) 
     1           THING(1:LTH), ( NIDX( K ), K = 1, NUM )
          ENDIF
        ENDIF
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 641
C
C  Get an interval for the number
C  a) lower bound
C
      ELSE IF ( T .EQ. 'I' ) THEN
 800    CONTINUE
        LOW = 0
        UPP = UMAX
        WRITE ( STDOUT, 4003 ) UTHING(1:LTH), UNDER(1:LTH)
 801    CONTINUE
        WRITE ( STDOUT, 4004 ) THING(1:LTH), THING(1:LTH)
        READ ( STDIN, '( A80 ) ', ERR = 802 ) LINE
        CHOICE = UPPER( LINE(1:1) )
        IF ( CHOICE .NE. ' ' ) THEN
          LOW = CONVERT( LINE(1:9),9 )
          IF ( LOW .LT. 0 .OR. LOW. GT. UMAX ) GO TO 802
        ENDIF
        GO TO 805
 802    CONTINUE
        WRITE ( STDOUT, 1001 )
        GO TO 801
C
C  b) upper bound
C
 805    CONTINUE
        WRITE ( STDOUT, 4005 ) UTHING(1:LTH), UNDER(1:LTH)
 803    CONTINUE
        WRITE ( STDOUT, 4006 ) THING(1:LTH), THING(1:LTH)
        READ ( STDIN, '( A80 ) ', ERR = 804 ) LINE
        CHOICE = UPPER( LINE(1:1) )
        IF ( CHOICE .NE. ' ' ) THEN
          UPP = CONVERT( LINE(1:9),9 )
          IF ( UPP .LT. 0 .OR. UPP. GT. UMAX .OR. LOW. GT. UPP )
     1       GO TO 804
        ENDIF
        GO TO 806
 804    CONTINUE
        WRITE ( STDOUT, 1001 )
        GO TO 803
C
C  Verify the bounds on the number
C
 806    CONTINUE
        WRITE ( STDOUT, 4007 ) THING(1:LTH), LOW, UPP
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 800
      ENDIF
      RETURN
C
C  Non executable statements
C
 1001 FORMAT( '   *** YOUR ANSWER IS NOT ALLOWED.  Please choose again.'
     1        / )
 1101 FORMAT( '   *** THIS SELECTION IS A REPETITION.  Please choose ',
     1             'again.'
     1        / )
 1014 FORMAT( /'   NUMBER OF ',A,':'
     1        /'   -----------',A )
 1015 FORMAT( /'     F   : Fixed                   V : Variable'
     1        /'     I   : In an interval'
     1        /'    <CR> : Any number of ',A,' (*)'
     1        /' '
     1        /'   Your choice :' )
 1018 FORMAT( /'   SELECT A NUMBER OF ', A,':'
     1        /'   --------------------', A )
 1119 FORMAT( /'    (INT) : Select only problems with (INT) ',A10,
     1        /'            (minimum 0, maximum 999999999, multiple ',
     1               'choices are allowed)'
     1        /'    <CR>  : Any fixed number of ',A,' (*)'
     1        /' '
     1        /'   Your choice :' )
 1019 FORMAT( /'    (INT) : Select only problems with (INT) ', A,
     1        /'            (minimum 0, maximum 999999999, multiple ',
     1               'choices are allowed)'
     1        /'     *    : Any fixed number of ', A,
     1        /'    <CR>  : No further selection',
     1        /' '
     1        /'   Your choice :' )
 2003 FORMAT(  '   Do you wish to reconsider your choice [<CR> = N] ?',
     1             ' (N/Y)')
 2009 FORMAT( /'   You have specified a fixed number of ',A,'.' )
 2019 FORMAT( /'   You have specified a variable number of ',A,'.' )
 2029 FORMAT( /'   You have specified any number of ',A,'.' )
 2010 FORMAT( /'   You have specified a number of ',A,
     1         ' in the set: ',
     1        /'       ', 6( I9, 1X ) )
 2020 FORMAT( /'   You have specified any fixed number of ',A,'.' )
 2110 FORMAT( /'   *** Please choose a number of ',A,'.' )
 4003 FORMAT( /'   LOWER BOUND ON THE NUMBER OF ', A,':'
     1        /'   ------------------------------', A )
 4004 FORMAT( /'    (INT) : Problems with at least (INT) ',A,
     1        /'    <CR>  : No lower bound on the number of ', A,
     1        /' '
     1        /'   Your choice : ' )
 4005 FORMAT( /'   UPPER BOUND ON THE NUMBER OF ', A,':'
     1        /'   ------------------------------', A )
 4006 FORMAT( /'    (INT) : Problems with at most (INT) ',A,
     1        /'    <CR>  : No upper bound on the number of ',A,
     1        /' '
     1        /'   Your choice : ' )
 4007 FORMAT( /'   You have specified a number of ', A,
     1        /'       in the interval [ ',I9,', ',I9,' ]' )
 4009 FORMAT( /'   You have specified an interval for the number of',
     1        /'       ',A,'.' )
C
      END
C
C-----------------------------------------------------------------------------
C
      LOGICAL FUNCTION REJECT( CHOICE, ITEM, NUM )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to verify if the character CHOICE
C  is a valid specification for item ITEM, where ITEM is one of
C         1 = objective type,
C         2 = constraints type,
C         3 = problem regularity,
C         4 = degree of available derivatives,
C         5 = problem interest,
C         6 = presence of explicit internal variables,
C         7 = number of variables,
C         8 = number of constraints,
C         9 = number of free variables,
C        10 = number of variables bounded below or above,
C        11 = number of variables bounded below and above,
C        12 = number of fixed variables,
C        13 = number of 1-sided linear constraints,
C        14 = number of 2-sided linear constraints,
C        15 = number of linear equality constraints,
C        16 = number of 1-sided general constraints,
C        17 = number of 2-sided gereral constraints,
C        18 = number of general equality constraints,
C
C  Programming: A. R. Conn and Ph. Toint for CGT Productions.
C
C-----------------------------------------------------------------------------
C
C
C  Classification constants
C
      INTEGER      OBJ, CON, REG, DER, INTRST, INTVAR, VARN, VARM
      PARAMETER  ( OBJ    = 1, CON  = 2, REG  = 3, DER = 4, INTRST = 5,
     1             INTVAR = 6, VARN = 7 )
C
C  Arguments
C
      CHARACTER*1  CHOICE
      INTEGER      ITEM, NUM
C
C  Other variables
C
      LOGICAL      ADMIT
C
C  Control choices
C
      IF ( ( CHOICE .EQ. ' ' .OR. CHOICE .EQ. '*' )
     1                              .AND. NUM .EQ. 0 ) THEN
        CHOICE = '*'
        REJECT = .FALSE.
        RETURN
C
C  Objective function type
C
      ELSE IF ( ITEM .EQ. OBJ ) THEN
        ADMIT = CHOICE .EQ. 'N' .OR. CHOICE .EQ. 'C' .OR.
     1          CHOICE .EQ. 'L' .OR. CHOICE .EQ. 'Q' .OR.
     1          CHOICE .EQ. 'S' .OR. CHOICE .EQ. 'O' .OR.
     1          CHOICE .EQ. ' '
C
C  Constraint type
C
      ELSE IF ( ITEM .EQ. CON ) THEN
        ADMIT = CHOICE .EQ. 'U' .OR. CHOICE .EQ. 'B' .OR.
     1          CHOICE .EQ. 'N' .OR. CHOICE .EQ. 'L' .OR.
     1          CHOICE .EQ. 'Q' .OR. CHOICE .EQ. 'O' .OR.
     1          CHOICE .EQ. 'X' .OR. CHOICE .EQ. ' '
C
C  Problem regularity
C
      ELSE IF ( ITEM .EQ. REG ) THEN
        ADMIT = CHOICE .EQ. 'R' .OR. CHOICE .EQ. 'I'
C
C  Degree of analytical derivatives
C
      ELSE IF ( ITEM .EQ. DER ) THEN
        ADMIT = CHOICE .EQ. '0' .OR. CHOICE .EQ. '1' .OR.
     1          CHOICE .EQ. '2' .OR. CHOICE .EQ. ' '
C
C  Problem interest
C
      ELSE IF ( ITEM .EQ. INTRST ) THEN
        ADMIT = CHOICE .EQ. 'A' .OR. CHOICE .EQ. 'R' .OR.
     1          CHOICE .EQ. 'M' .OR. CHOICE .EQ. ' '
C
C  Presence of explicit internal variables
C
      ELSE IF ( ITEM .EQ. INTVAR ) THEN
         ADMIT = CHOICE .EQ. 'Y' .OR. CHOICE .EQ. 'N'
C
C  Number of variables and constraints
C
      ELSE IF ( ITEM .GE. VARN ) THEN
        ADMIT = CHOICE .EQ. 'F' .OR. CHOICE .EQ. 'V' .OR.
     1          CHOICE .EQ. 'I'
      ENDIF
      REJECT = .NOT. ADMIT
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      LOGICAL FUNCTION MATCH( C, T, MAXTRG, ANYFIX, CHOSEN, LOW, UPP )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to check if the classification C
C  matches one of the targets specified by T, ANYFIX, CHOSEN, LOW and
C  UPP.
C
C  Programming: Ph. Toint, Dec 1991, for CGT Productions. Revised Aug 2005.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      CHARACTER*8  T( MAXTRG )
      CHARACTER*27 C
      LOGICAL      ANYFIX( 18 )
      INTEGER      CHOSEN( MAXTRG, 18 ), MAXTRG, LOW( 18 ), UPP( 18 )
C
C  Other variables
C
      CHARACTER*8    PBNAME
      CHARACTER*1    OBJTYP,  CONTYP, REGTYP, DERLVL, INTRST, INTVAR
      INTEGER        PBN, PBM, PBNFX,  PBN1S,  PBN2S, PBNFR, PBM1SL, 
     1               PBM2SL, PBMEQL,  PBM1SG, PBM2SG, PBMEQG
      INTEGER        I, J, CONVERT, LSTN, FRSTM, LASTM
      LOGICAL        ERROR, MATCHD
      CHARACTER*1    CH
      INTEGER        WRONG, UNKNWN, VARIAB
      PARAMETER    ( WRONG = -10, UNKNWN = -5, VARIAB = -1 )

      INTEGER        OBJ, CON, REG, DER, INTR, INTV, VARN, VARM,  
     *               VFR, V1S, V2S, VFX, C1SL, C2SL, CEQL, C1SG, C2SG,
     *               CEQG
      PARAMETER    ( OBJ  =  1, CON  =  2, REG  =  3, DER  = 4, 
     *               INTR =  5, INTV =  6, VARN =  7, VARM = 8,
     *               VFR  =  9, V1S  = 10, V2S  = 11, VFX  = 12, 
     *               C1SL = 13, C2SL = 14, CEQL = 15, C1SG = 16, 
     *               C2SG = 17, CEQG = 18 )
C
      MATCH = .FALSE.
C
C  Parse the current classification record.
C
c      write(*,*) ' C = ', C
      CALL PARSE ( C, PBNAME, OBJTYP, CONTYP, REGTYP,  
     1             DERLVL, INTRST, INTVAR, PBN, PBM, PBNFR,
     2             PBN1S, PBN2S, PBNFX, PBM1SL, PBM2SL, PBMEQL,
     3             PBM1SG, PBM2SG, PBMEQG, ERROR )
      IF ( ERROR ) RETURN
c      write(*,*) ' PBMEQL = ', PBMEQL
c      write(*,*) ' PBMEQG = ', PBMEQG
C
C  Match objective function type
C
      DO 10 I = 1, MAXTRG
        CH = T(I)(OBJ:OBJ)
        IF ( CH .EQ. ' ' ) GO TO 20
        MATCH = MATCH .OR. CH .EQ. OBJTYP .OR. CH .EQ. '*'
 10   CONTINUE
 20   CONTINUE
      IF ( .NOT. MATCH ) RETURN
C
C  Match constraint type
C
      MATCH = .FALSE.
      DO 30 I = 1, MAXTRG
        CH = T(I)(CON:CON)
        IF ( CH .EQ. ' ' ) GO TO 40
        MATCH = MATCH .OR. CH .EQ. CONTYP .OR. CH .EQ. '*'
 30   CONTINUE
 40   CONTINUE
      IF ( .NOT. MATCH ) RETURN
C
C  Match regularity type
C
      CH = T(1)(REG:REG)
      MATCH = CH .EQ. REGTYP .OR. CH .EQ. '*'
      IF ( .NOT. MATCH ) RETURN
C
C  Match degree of available derivatives
C
      MATCH = .FALSE.
      DO 50 I = 1, MAXTRG
        CH = T(I)(DER:DER)
        IF ( CH .EQ. ' ' ) GO TO 60
        MATCH = MATCH .OR. CH .EQ. DERLVL .OR. CH .EQ. '*'
 50   CONTINUE
 60   CONTINUE
      IF ( .NOT. MATCH ) RETURN
C
C  Match interest of the problem
C
      MATCH = .FALSE.
      DO 70 I = 1, MAXTRG
        CH = T(I)(INTR:INTR)
        IF ( CH .EQ. ' ' ) GO TO 80
        MATCH = MATCH .OR. CH .EQ. INTRST .OR. CH .EQ. '*'
 70   CONTINUE
 80   CONTINUE
      IF ( .NOT. MATCH ) RETURN
C
C  Match for explicit internal variables
C
      CH = T(1)(INTV:INTV)
      MATCH = CH .EQ. INTVAR .OR. CH .EQ. '*'
      IF ( .NOT. MATCH ) RETURN
C
C  Match the numbers of variables
C
      MATCH = MATCHD( PBNFR, T(1)(VFR:VFR), ANYFIX(VFR), LOW(VFR), 
     1                UPP(VFR), CHOSEN(1,VFR), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBN1S, T(1)(V1S:V1S), ANYFIX(V1S), LOW(V1S), 
     1                UPP(V1S), CHOSEN(1,V1S), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBN2S, T(1)(V2S:V2S), ANYFIX(V2S), LOW(V2S), 
     1                UPP(V2S), CHOSEN(1,V2S), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBNFX, T(1)(VFX:VFX), ANYFIX(VFX), LOW(VFX), 
     1                UPP(VFX), CHOSEN(1,VFX), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBN, T(1)(VARN:VARN), ANYFIX(VARN), LOW(VARN), 
     1                UPP(VARN), CHOSEN(1,VARN), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
C  Match the numbers of constraints
C
      MATCH = MATCHD( PBM1SL, T(1)(C1SL:C1SL), ANYFIX(C1SL), LOW(C1SL), 
     1                UPP(C1SL), CHOSEN(1,C1SL), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBM2SL, T(1)(C2SL:C2SL), ANYFIX(C2SL), LOW(C2SL), 
     1                UPP(C2SL), CHOSEN(1,C2SL), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBMEQL, T(1)(CEQL:CEQL), ANYFIX(CEQL), LOW(CEQL), 
     1                UPP(CEQL), CHOSEN(1,CEQL), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBM1SG, T(1)(C1SG:C1SG), ANYFIX(C1SG), LOW(C1SG), 
     1                UPP(C1SG), CHOSEN(1,C1SG), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBM2SG, T(1)(C2SG:C2SG), ANYFIX(C2SG), LOW(C2SG), 
     1                UPP(C2SG), CHOSEN(1,C2SG), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBMEQG, T(1)(CEQG:CEQG), ANYFIX(CEQG), LOW(CEQG), 
     1                UPP(CEQG), CHOSEN(1,CEQG), MAXTRG )
      IF ( .NOT. MATCH )  RETURN
C
      MATCH = MATCHD( PBM, T(1)(VARM:VARM), ANYFIX(VARM), LOW(VARM), 
     1                UPP(VARM), CHOSEN(1,VARM), MAXTRG )

      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      SUBROUTINE PARSE ( CLASSF, PBNAME, OBJTYP, CONTYP, REGTYP,  
     1                   DERLVL, INTRST, INTVAR, PBN, PBM, PBNFR,
     2                   PBN1S, PBN2S, PBNFX, PBM1SL, PBM2SL, PBMEQL,
     3                   PBM1SG, PBM2SG, PBMEQG, ERROR )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to parse one line of the CLASSF.DB
C  database into its components, i.e. a problem name and the problem's
C  caracteristics.  It assumes the separators are either one hyphen or
C  a string of at least one blank.
C
C  Programming: Ph. Toint, Aug 2005.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      CHARACTER*200  CLASSF
      CHARACTER*8    PBNAME
      CHARACTER*1    OBJTYP,  CONTYP, REGTYP, DERLVL, INTRST, INTVAR
      INTEGER        PBN, PBM, PBNFX,  PBN1S,  PBN2S, PBNFR, PBM1SL, 
     1               PBM2SL, PBMEQL,  PBM1SG, PBM2SG, PBMEQG
      LOGICAL        ERROR
C
C  Other variables
C
      INTEGER        IOBJ, IINT, IN, IM, INFX, IN1S, IN2S, INFR, IM1SL,
     1               IM2SL, IMEQL, IM1SG, IM2SG, IMEQG, JOBJ, JINT, JN,
     1               JM, JNFX, JN1S, JN2S, JNFR, JM1SL, JM2SL, JMEQL, 
     1               JM1SG, JM2SG, JMEQG, J, L, MAXL
      INTEGER        CONVERT
      PARAMETER    ( MAXL = 200 )
      INTEGER        WRONG, UNKNWN, VARIAB
      PARAMETER    ( WRONG = -10, UNKNWN = -5, VARIAB = -1 )
C
      ERROR = .FALSE.
C
C  Get the problem name.
C
      PBNAME = CLASSF(1:8)
C
C  Find the start of each field.
C
      CALL NXTSTR( CLASSF, 9, IOBJ, JOBJ )
C      write(*,*) ' IOBJ, JOBJ = ', IOBJ, JOBJ
      IF ( IOBJ .LT. 0 ) GO TO 40
      CALL NXTSTR( CLASSF, JOBJ+2, IINT, JINT )
C      write(*,*) ' IINT, JINT = ', IINT, JINT
      IF ( IINT .LT. 0 ) GO TO 40
      CALL NXTSTR( CLASSF, JINT+2, IN, JN )
C      write(*,*) ' IN, JN = ', IN, JN
      IF ( IN .LT. 0 ) GO TO 40
      CALL NXTSTR( CLASSF, JN+2, IM, JM )
C      write(*,*) ' IM, JM = ', IM, JM
      IF ( IM .LT. 0 ) GO TO 40
      CALL NXTSTR( CLASSF, JM+2, INFR, JNFR )
C      write(*,*) ' INFR,  JNFR = ', INFR, JNFR
      IF ( INFR .GT. 0 ) THEN
         CALL NXTSTR( CLASSF, JNFR+2, IN1S, JN1S )
C        write(*,*) ' IN1S,  JN1S = ', IN1S, JN1S
         IF ( IN1S .LT. 0 ) GO TO 40
         CALL NXTSTR( CLASSF, JN1S+2, IN2S, JN2S )
C        write(*,*) ' IN2S,  JN2S = ', IN2S, JN2S
         IF ( IN2S .LT. 0 ) GO TO 40
         CALL NXTSTR( CLASSF, JN2S+2, INFX, JNFX )
C        write(*,*) ' INFX,  JNFX = ', INFX, JNFX
         IF ( INFX .LT. 0 ) GO TO 40
         CALL NXTSTR( CLASSF, JNFX+2, IM1SL, JM1SL )
C        write(*,*) ' IM1SL, JM1SL = ', IM1SL, JM1SL
         IF ( IM1SL .LT. 0 ) GO TO 40
         CALL NXTSTR( CLASSF, JM1SL+2, IM2SL, JM2SL )
C        write(*,*) ' IM2SL, JM2SL = ', IM2SL, JM2SL
         IF ( IM2SL .LT. 0 ) GO TO 40
         CALL NXTSTR( CLASSF, JM2SL+2, IMEQL, JMEQL )
C        write(*,*) ' IMEQL, JMEQL = ', IMEQL, JMEQL
         IF ( IMEQL .LT. 0 ) GO TO 40
         CALL NXTSTR( CLASSF, JMEQL+2, IM1SG, JM1SG )
C        write(*,*) ' IM1SG, JM1SG = ', IM1SG, JM1SG
         IF ( IM1SG .LT. 0 ) GO TO 40
         CALL NXTSTR( CLASSF, JM1SG+2, IM2SG, JM2SG )
C        write(*,*) ' IM2SG, JM2SG = ', IM2SG, JM2SG
         IF ( IM2SG .LT. 0 ) GO TO 40
         CALL NXTSTR( CLASSF, JM2SG+2, IMEQG, JMEQG )
C        write(*,*) ' IMEQG, JMEQG = ', IMEQG, JMEQG
         IF ( IMEQG .LT. 0 ) GO TO 40
      END IF
C 
C  Decrypt the objective string.
C
      L = IOBJ
      OBJTYP = CLASSF(L:L)
      L = L + 1
      CONTYP = CLASSF(L:L)
      L = L + 1
      REGTYP = CLASSF(L:L)
      L = L + 1
      DERLVL = CLASSF(L:L)
C 
C  Decrypt the interest string.
C
      L = IINT
      INTRST = CLASSF(L:L)
      L = L + 1
      INTVAR = CLASSF(L:L)
C
C  Convert the remaining fields, if present.
C
      PBN = CONVERT( CLASSF(IN:JN), JN-IN+1 )
      IF ( PBN .EQ. WRONG ) GO TO 40
      PBM = CONVERT( CLASSF(IM:JM), JM-IM+1 )
      IF ( PBM .EQ. WRONG ) GO TO 40
      IF ( INFR .GT. 0 ) THEN
         PBNFR  = CONVERT( CLASSF(INFR:JNFR), JNFR-INFR+1 )
         IF ( PBNFR .EQ. WRONG ) GO TO 40
         PBN1S  = CONVERT( CLASSF(IN1S:JN1S), JN1S-IN1S+1 )
         IF ( PBN1S .EQ. WRONG ) GO TO 40
         PBN2S  = CONVERT( CLASSF(IN2S:JN2S), JN2S-IN2S+1 )
         IF ( PBN2S .EQ. WRONG ) GO TO 40
         PBNFX  = CONVERT( CLASSF(INFX:JNFX), JNFX-INFX+1 )
         IF ( PBNFX .EQ. WRONG ) GO TO 40
         PBM1SL = CONVERT( CLASSF(IM1SL:JM1SL), JM1SL-IM1SL+1 )
         IF ( PBM1SL .EQ. WRONG ) GO TO 40
         PBM2SL = CONVERT( CLASSF(IM2SL:JM2SL), JM2SL-IM2SL+1 )
         IF ( PBM2SL .EQ. WRONG ) GO TO 40
         PBMEQL = CONVERT( CLASSF(IMEQL:JMEQL), JMEQL-IMEQL+1 )
         IF ( PBMEQL .EQ. WRONG ) GO TO 40
         PBM1SG = CONVERT( CLASSF(IM1SG:JM1SG), JM1SG-IM1SG+1 )
         IF ( PBM1SG .EQ. WRONG ) GO TO 40
         PBM2SG = CONVERT( CLASSF(IM2SG:JM2SG), JM2SG-IM2SG+1 )
         IF ( PBM2SG .EQ. WRONG ) GO TO 40
         PBMEQG = CONVERT( CLASSF(IMEQG:JMEQG), JMEQG-IMEQG+1 )
         IF ( PBMEQG .EQ. WRONG ) GO TO 40
      ELSE
         PBNFR  = UNKNWN
         PBN1S  = UNKNWN
         PBN2S  = UNKNWN
         PBNFX  = UNKNWN
         PBM1SL = UNKNWN
         PBM2SL = UNKNWN
         PBMEQL = UNKNWN
         PBM1SG = UNKNWN
         PBM2SG = UNKNWN
         PBMEQG = UNKNWN
      END IF
      RETURN
C
C  Error
C
 40   CONTINUE
      ERROR = .TRUE.
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      LOGICAL FUNCTION MATCHD( PVAL, T, ANYFIX, LOW, UPP, 
     *                         CHOSEN, MAXTRG )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to decide if the information 
C  on a problem contained in PVAL matches the pattern specified 
C  by ANYFIX, LOW, UPP and CHOSEN.
C
C  Programming: Ph. Toint, Aug 2005.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      INTEGER      PVAL, LOW, UPP, MAXTRG, CHOSEN( MAXTRG )
      LOGICAL      ANYFIX
      CHARACTER*1  T
C
      INTEGER        VARIAB
      PARAMETER    ( VARIAB = -1 )
C
C  Match the number of constraints
C
      MATCHD = .FALSE.
      IF ( T .EQ. '*' ) THEN
        MATCHD = .TRUE.
      ELSE IF ( T .EQ. 'V' ) THEN
        MATCHD = PVAL .EQ. VARIAB
      ELSE IF ( T .EQ. 'I' ) THEN
C    interval.  A variable number of constraints is no longer considered to
C    match a number in an interval.
        IF ( PVAL .NE. VARIAB ) THEN
          MATCHD = PVAL. GE. LOW .AND. PVAL .LE. UPP
        ENDIF
      ELSE
        IF ( ANYFIX ) THEN
          MATCHD = .TRUE.
        ELSE
C    A variable number of constraints is no longer considered to match a fixed
C    number.
          IF ( PVAL .NE. VARIAB ) THEN
            DO 110 I = 1, MAXTRG
              IF ( CHOSEN(I) .LT. 0 ) GO TO 110
              MATCHD = CHOSEN(I) .EQ. PVAL
              IF ( MATCHD ) GO TO 120
 110        CONTINUE
 120        CONTINUE
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      SUBROUTINE NXTSTR ( STRING, IPOS, ISTART, ISTOP )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to find the position of the next separator
C  ('-' or ' ') in STRING, after the next meaningful substring (i.e. not
C  containing separators).  If no meaningful substring is found, NXTSEP = -1
C  is returned.
C
C  Programming: Ph. Toint, Aug 2005.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      CHARACTER*200 STRING
      INTEGER       IPOS, ISTART, ISTOP
C
C  Other variables
C
      INTEGER       J, MAXL
      PARAMETER   ( MAXL = 200 )
C
C  Remove initial blanks
C
      ISTOP  = -1
      DO 5  ISTART = IPOS, MAXL
         IF ( STRING(ISTART:ISTART) .NE. ' ' ) GO TO 15
 5    CONTINUE
      ISTART = -1
      RETURN
C
C  Find the end of the meaningful string.
C
 15   CONTINUE
      J = ISTART - 1
 100  CONTINUE
      J = J + 1
         IF ( STRING(J:J) .NE. ' ' .AND. STRING(J:J) .NE. '-' ) THEN
            ISTOP = 0
         ELSE
            IF ( ISTOP .EQ. 0 ) THEN
               ISTOP = J - 1
               RETURN
            END IF
         END IF
      IF ( J .LE. MAXL ) GO TO 100
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      INTEGER FUNCTION CONVERT( LINE, LENL )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to convert the nonnegative integer contained
C  in the string LINE into a proper integer.
C  UNKNWN is returned in case of unknown number, 
C  VARIAB for a variable number,
C  WRONG  in case of error.
C
C  Programming: Ph. Toint, Aug 2005, for CGT Productions.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      CHARACTER*200    LINE
      INTEGER          LENL
C
C  Other variables
C
      INTEGER          IS, IE, L
      INTEGER          WRONG, UNKNWN, VARIAB
      PARAMETER      ( WRONG = -10, UNKNWN = -5, VARIAB = -1 )
C
C  Remove trailing blanks
C
      DO 5  IS = 1, LENL
         IF ( LINE(IS:IS) .NE. ' ' ) GO TO 15
 5    CONTINUE
      GO TO 30
 15   CONTINUE
      DO 10 IE= LENL, 1, -1
         IF ( LINE(IE:IE) .NE. ' ' ) GO TO 20
 10   CONTINUE
      GO TO 30
C
C  Read the integer
C
 20   CONTINUE
      L = IE - IS + 1
      IF ( L .EQ. 1 ) THEN
        IF ( LINE(IS:IE) .EQ. '?' ) THEN
           CONVERT = UNKNWN
        ELSE IF ( LINE(IS:IE) .EQ. 'V' ) THEN
           CONVERT = VARIAB
        ELSE
           READ ( LINE(IS:IE), '( I1 )', ERR = 30 ) CONVERT
        END IF
      ELSE IF ( L .EQ. 2 ) THEN
        READ ( LINE(IS:IE), '( I2 )', ERR = 30 ) CONVERT
      ELSE IF ( L .EQ. 3 ) THEN
        READ ( LINE(IS:IE), '( I3 )', ERR = 30 ) CONVERT
      ELSE IF ( L .EQ. 4 ) THEN
        READ ( LINE(IS:IE), '( I4 )', ERR = 30 ) CONVERT
      ELSE IF ( L .EQ. 5 ) THEN
        READ ( LINE(IS:IE), '( I5 )', ERR = 30 ) CONVERT
      ELSE IF ( L .EQ. 6 ) THEN
        READ ( LINE(IS:IE), '( I6 )', ERR = 30 ) CONVERT
      ELSE IF ( L .EQ. 7 ) THEN
        READ ( LINE(IS:IE), '( I7 )', ERR = 30 ) CONVERT
      ELSE IF ( L .EQ. 8 ) THEN
        READ ( LINE(IS:IE), '( I8 )', ERR = 30 ) CONVERT
      ELSE
        READ ( LINE(IS:IE), '( I9 )', ERR = 30 ) CONVERT
      ENDIF
      RETURN
C
C  Error
C
 30   CONTINUE
      CONVERT = WRONG
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      INTEGER FUNCTION NBT ( I, T, MX )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to return the number of non trivial
C  choices of item I specified in the target list T.
C
C  Programming: Ph. Toint, Dec 1991, for CGT Productions.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      CHARACTER*18 T( MX )
      INTEGER      I, MX
C
C  Other variables
C
      INTEGER L
C
      NBT = 1
      DO 10 L = 2, MX
        IF ( T(L)(I:I) .EQ. 'X' ) RETURN
        NBT = NBT + 1
 10   CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      INTEGER FUNCTION NBI ( N, MX )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to return the number of non trivial
C  dimensions specified in the integer dimension vector N.
C
C  Programming: Ph. Toint, Dec 1991, for CGT Productions.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      INTEGER N( MX ), MX
C
C  Other variables
C
      INTEGER L
C
      NBI = 0
      DO 10 L = 1, MX
        IF ( N(L) .LT. 0 ) RETURN
        NBI = NBI + 1
 10   CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      CHARACTER*1 FUNCTION UPPER( CH ) 
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to transform the character CH
C  to upper case, if it is not already.
C
C  Programming: Ph. Toint, Dec 1991, for CGT Productions.
C
C-----------------------------------------------------------------------------
C
      CHARACTER*1  CH
C
      INTRINSIC    ICHAR, CHAR
C
      INTEGER      ICH, LSTART
C
      LSTART = ICHAR( 'a' )
      ICH    = ICHAR( CH )
      IF ( ICH .GE. LSTART .AND. ICH .LE. ICHAR( 'z' ) ) THEN
        UPPER = CHAR( ICHAR( 'A' ) + ICH - LSTART )
      ELSE
        UPPER = CH
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
