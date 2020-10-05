C     ( Last modified on 14 Jan 2001 at 19:03:33 )
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
C
C--------- THE FOLLOWING SPECIFICATIONS MAY BE MODIFIED BY THE USER ----------
C
C  Standard input default definition is device 5.  Standard output default
C  definition is device 6. Change them to whatever values are appropriate
C  on your system.
C
      INTEGER      STDIN,     STDOUT
      PARAMETER  ( STDIN = 5, STDOUT = 6 )
C
C  Name of the classification database
C
      CHARACTER*32 NAME
      PARAMETER  ( NAME = 'CLASSF.DB' )
C
C  Device number for reading the classification database
C
      INTEGER      CLSDVC, FLSDVC
      PARAMETER  ( CLSDVC = 55, FLSDVC = 56 )
C
C    default directory for classification file.
C
C  Device number and file name for file containing default directory
C
      INTEGER      DATDVC
      PARAMETER  ( DATDVC = 57 )
      CHARACTER*32 DATNAM
      PARAMETER  ( DATNAM = 'SLCT.DAT' )
C
C---------------- END OF THE USER MODIFIABLE SPECIFICATION ------------------
C
C  Classification constants
C
      INTEGER      OBJ, CON, REG, DER, INTRST, INTVAR, VARN, VARM
      PARAMETER  ( OBJ    = 1, CON  = 2, REG  = 3, DER = 4, INTRST = 5,
     *             INTVAR = 6, VARN = 7, VARM = 8 )
C
C  Maximum number of simultaneous targets in search
C
      INTEGER      MAXTRG
      PARAMETER  ( MAXTRG = 7 )
C
C  Variable definitions
C
      CHARACTER*80 LINE
C    names for classification file
      CHARACTER*72 FILEN
C    default directory for classification file
      CHARACTER*256 DFTDIR
C    names for output listing file
C                  Addition by Kristjan Jonasson
      CHARACTER*72 FILES
      CHARACTER*28 PBCLS
      CHARACTER*8  TARGET( MAXTRG ), LIST(5)
      CHARACTER*1  CHOICE, CHAR, UPPER
      INTEGER      I, IM1, J, NVAR( MAXTRG ), NCON( MAXTRG ), L, K, NUM,
     *             NMATCH, CONVERT, NBT,  NBI, LN, UN, LM, SIZE, UM
      LOGICAL      REJECT, ANYFNV,  ANYFNC, MATCH
      INTRINSIC    MIN
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
        NVAR(I)   = -2
        NCON(I)   = -2
        TARGET(I) = 'XXXXXXXX'
  4   CONTINUE
      TARGET(1) = '********'
      ANYFNV = .FALSE.
      ANYFNC = .FALSE.
C
C  Bounds on the number of variables and constraints are initialized
C  to be inactive.
C
      LN     = 0
      LM     = 0
      UN     = 99999999
      UM     = 99999999
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
      IF ( TARGET(1)(VARN:VARN) .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5009 )
      ELSE IF ( TARGET(1)(VARN:VARN) .EQ. '*' ) THEN
        WRITE ( STDOUT, 5010 )
      ELSE IF ( TARGET(1)(VARN:VARN) .EQ. 'I' ) THEN
        WRITE ( STDOUT, 4008 ) LN, UN
      ELSE
        IF ( ANYFNV ) THEN
          WRITE ( STDOUT, 5011 )
        ELSE
          NUM = NBI ( NVAR, MAXTRG )
          WRITE ( STDOUT, 5012 ) ( NVAR(L), l = 1, NUM )
        ENDIF
      ENDIF
      IF ( TARGET(1)(VARM:VARM) .EQ. 'V' ) THEN
        WRITE ( STDOUT, 5013 )
      ELSE IF ( TARGET(1)(VARM:VARM) .EQ. '*' ) THEN
        WRITE ( STDOUT, 5014 )
      ELSE IF ( TARGET(1)(VARM:VARM) .EQ. 'I' ) THEN
        WRITE ( STDOUT, 6008 ) LM, UM
      ELSE
        IF ( ANYFNV ) THEN
          WRITE ( STDOUT, 5015 )
        ELSE
          NUM = NBI ( NCON, MAXTRG )
          WRITE ( STDOUT, 5016 ) ( NCON(L), l = 1, NUM )
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
              IF ( TARGET(J)(OBJ:OBJ) .EQ. CHOICE ) GO TO 25
 35         CONTINUE
          ENDIF
          TARGET(I)(OBJ:OBJ) = CHOICE
          IF ( CHOICE .NE. ' ' ) NUM = NUM + 1
          IF ( CHOICE .EQ. ' ' .OR. CHOICE .EQ. '*' ) GO TO 40
          GO TO 30
 20       CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 10
 25       CONTINUE
          WRITE ( STDOUT, 1101 )
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
              IF ( TARGET(J)(CON:CON) .EQ. CHOICE ) GO TO 125
 135        CONTINUE
          ENDIF
          TARGET(I)(CON:CON) = CHOICE
          IF ( CHOICE .NE. ' ' ) NUM = NUM + 1
          IF ( CHOICE .EQ. ' ' .OR. CHOICE .EQ. '*' ) GO TO 140
          GO TO 130
 120      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 110
 125      CONTINUE
          WRITE ( STDOUT, 1101 )
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
              IF ( TARGET(J)(DER:DER) .EQ. CHOICE ) GO TO 325
 335         CONTINUE
          ENDIF
          TARGET(I)(DER:DER) = CHOICE
          IF ( CHOICE .NE. ' ' ) NUM = NUM + 1
          IF ( CHOICE .EQ. ' ' .OR. CHOICE .EQ. '*' ) GO TO 340
          GO TO 330
 320      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 310
 325      CONTINUE
          WRITE ( STDOUT, 1101 )
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
              IF ( TARGET(J)(INTRST:INTRST) .EQ. CHOICE ) GO TO 425
 435        CONTINUE
          ENDIF
          TARGET(I)(INTRST:INTRST) = CHOICE
          IF ( CHOICE .NE. ' ' ) NUM = NUM + 1
          IF ( CHOICE .EQ. ' ' .OR. CHOICE .EQ. '*' ) GO TO 440
          GO TO 430
 420      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 410
 425      CONTINUE
          WRITE ( STDOUT, 1101 )
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
      ELSE IF ( CHAR .EQ. 'N' ) THEN
 611    CONTINUE
        NUM = 0
        WRITE ( STDOUT, 1014 )
 610    CONTINUE
        WRITE ( STDOUT, 1015 )
        READ  ( STDIN , FMT = '( A1 )', ERR = 620 ) CHOICE
        CHOICE = UPPER( CHOICE )
        IF ( REJECT( CHOICE, VARN, 0 ) ) GO TO 620
        TARGET(1)(VARN:VARN) = CHOICE
        GO TO 640
 620    CONTINUE
        WRITE ( STDOUT, 1001 )
        GO TO 610
C
C  Verify the number of variables type
C
 640    CONTINUE
        IF ( CHOICE .EQ. 'F' ) THEN
          WRITE ( STDOUT, 2009 )
        ELSE IF ( CHOICE .EQ. 'V' ) THEN
          WRITE ( STDOUT, 2019 )
        ELSE IF ( CHOICE .EQ. 'I' ) THEN
          WRITE ( STDOUT, 4009 )
        ELSE
          WRITE ( STDOUT, 2029 )
        ENDIF
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 611
C
C  Get the fixed number of variables
C
        IF ( TARGET(1)(VARN:VARN) .EQ. 'F' ) THEN
 641      CONTINUE
          NUM = 0
          DO 670 I = 1, MAXTRG
            WRITE( STDOUT, 1018 )
 650        CONTINUE
            IF ( NUM .EQ. 0) WRITE( STDOUT, 1119 )
            IF ( NUM .GT. 0) WRITE( STDOUT, 1019 )
            READ ( STDIN , '( A80 )', ERR = 660 ) LINE
            CHOICE = UPPER( LINE(1:1) )
            IF ( NUM .EQ. 0 ) ANYFNV =
     1           CHOICE .EQ. '*' .OR. CHOICE .EQ. ' '
            IF ( NUM .NE. 0 ) ANYFNV = CHOICE .EQ. '*'
            IF ( CHOICE .EQ. ' ' .OR. ANYFNV ) GO TO 695
            NVAR( I ) = CONVERT( LINE(1:5) )
            IF ( I .GT. 1 ) THEN
              IM1 = I - 1
              DO 655 J = 1, IM1
                IF ( NVAR(J) .EQ. NVAR(I) ) GO TO 665
 655          CONTINUE
            ENDIF
            IF ( NVAR(I) .LT. 0 .OR. NVAR(I) .GT. 99999999 ) GO TO 660
            NUM = NUM + 1
            GO TO 670
 660        CONTINUE
            WRITE ( STDOUT, 1001 )
            GO TO 650
 665        CONTINUE
            WRITE ( STDOUT, 1101 )
            GO TO 650
 670      CONTINUE
C
C   Verify the number of variables
C
 695      CONTINUE
          IF ( ANYFNV ) THEN
            WRITE ( STDOUT, 2020 )
          ELSE
            IF ( NUM .EQ. 0 ) THEN
              WRITE ( STDOUT, 2110 )
              GO TO 641
            ELSE
              WRITE ( STDOUT, 2010 ) ( NVAR(K), K = 1, NUM )
            ENDIF
          ENDIF
          WRITE ( STDOUT, 2003 )
          READ  ( STDIN, '( A1 )' ) CHOICE
          IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 641
C
C  Get an interval for the number of variables
C  a) lower bound
C
        ELSE IF ( TARGET(1)(VARN:VARN) .EQ. 'I' ) THEN
 800      CONTINUE
          LN = 0
          UN = 99999999
          WRITE ( STDOUT, 4003 )
 801      CONTINUE
          WRITE ( STDOUT, 4004 )
          READ ( STDIN, '( A80 ) ', ERR = 802 ) LINE
          CHOICE = UPPER( LINE(1:1) )
          IF ( CHOICE .NE. ' ' ) THEN
            LN = CONVERT( LINE(1:5) )
            IF ( LN .LT. 0 .OR. LN. GT. 99999999 ) GO TO 802
          ENDIF
          GO TO 805
 802      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 801
C
C  b) upper bound
C
 805      CONTINUE
          WRITE ( STDOUT, 4005 )
 803      CONTINUE
          WRITE ( STDOUT, 4006 )
          READ ( STDIN, '( A80 ) ', ERR = 804 ) LINE
          CHOICE = UPPER( LINE(1:1) )
          IF ( CHOICE .NE. ' ' ) THEN
            UN = CONVERT( LINE(1:5) )
            IF ( UN .LT. 0 .OR. UN. GT. 99999999 .OR. LN. GT. UN )
     1         GO TO 804
          ENDIF
          GO TO 806
 804      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 803
C
C  Verify the bounds on the number of variables
C
 806      CONTINUE
          WRITE ( STDOUT, 4007 ) LN, UN
          WRITE ( STDOUT, 2003 )
          READ  ( STDIN, '( A1 )' ) CHOICE
          IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 800
        ENDIF
C
C  Get the number of constraints
C
      ELSE IF ( CHAR .EQ. 'M' ) THEN
        NUM = 0
 711    CONTINUE
        WRITE ( STDOUT, 1016 )
 710    CONTINUE
        WRITE ( STDOUT, 1017 )
        READ  ( STDIN , FMT = '( A1 )', ERR = 720 ) CHOICE
        CHOICE = UPPER( CHOICE )
        IF ( REJECT( CHOICE, VARM, 0 ) ) GO TO 720
        TARGET(1)(VARM:VARM) = CHOICE
        GO TO 740
 720    CONTINUE
        WRITE ( STDOUT, 1001 )
        GO TO 710
C
C  Verify the constraints number type
C
 740    CONTINUE
        IF ( CHOICE .EQ. 'F' ) THEN
          WRITE ( STDOUT, 2011 )
        ELSE IF ( CHOICE .EQ. 'V' ) THEN
          WRITE ( STDOUT, 2021 )
        ELSE IF ( CHOICE .EQ. 'I' ) THEN
          WRITE ( STDOUT, 6009 )
        ELSE
          WRITE ( STDOUT, 2031 )
        ENDIF
        WRITE ( STDOUT, 2003 )
        READ  ( STDIN, '( A1 )' ) CHOICE
        IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 711
C
C  Get the fixed number of constraints
C
        IF ( TARGET(1)(VARM:VARM) .EQ. 'F' ) THEN
 741      CONTINUE
          NUM = 0
          DO 770 I = 1, MAXTRG
            WRITE( STDOUT, 1020 )
 750        CONTINUE
            IF ( NUM .EQ. 0) WRITE( STDOUT, 1121 )
            IF ( NUM .GT. 0) WRITE( STDOUT, 1021 )
            READ ( STDIN , FMT = '( A80 )', ERR = 760 ) LINE
            CHOICE = UPPER( LINE(1:1) )
            IF ( NUM .EQ. 0 ) ANYFNC =
     1           CHOICE .EQ. '*' .OR. CHOICE .EQ. ' '
            IF ( NUM .NE. 0 ) ANYFNC = CHOICE .EQ. '*'
            IF ( CHOICE .EQ. ' ' .OR. ANYFNC ) GO TO 795
            NCON(I) = CONVERT( LINE(1:5) )
            IF ( I .GT. 1 ) THEN
              IM1 = I - 1
              DO 755 J = 1, IM1
                IF ( NCON(J) .EQ. NCON(I) ) GO TO 765
 755          CONTINUE
            ENDIF
            IF ( NCON(I) .LT. 0 .OR. NCON(I) .GT. 99999999 ) GO TO 760
            NUM = NUM + 1
            GO TO 770
 760        CONTINUE
            WRITE ( STDOUT, 1001 )
            GO TO 750
 765        CONTINUE
            WRITE ( STDOUT, 1101 )
            GO TO 750
 770      CONTINUE
C
C  Verify the number of constraints
C
 795      CONTINUE
          IF ( ANYFNC ) THEN
            WRITE ( STDOUT, 2022 )
          ELSE
            IF ( NUM .EQ. 0 ) THEN
              WRITE ( STDOUT, 2120 )
              GO TO 741
            ELSE
              WRITE ( STDOUT, 2012 ) ( NCON(K), K = 1, NUM )
            ENDIF
          ENDIF
          WRITE ( STDOUT, 2003 )
          READ  ( STDIN, '( A1 )' ) CHOICE
          IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 741
C
C  Get an interval for the number of constraints
C  a) lower bound
C
        ELSE IF ( TARGET(1)(VARM:VARM) .EQ. 'I' ) THEN
 900      CONTINUE
          LM = 0
          UM = 99999999
          WRITE ( STDOUT, 6003 )
 901      CONTINUE
          WRITE ( STDOUT, 6004 )
          READ ( STDIN, '( A80 ) ', ERR = 902 ) LINE
          CHOICE = UPPER( LINE(1:1) )
          IF ( CHOICE .NE. ' ' ) THEN
            LM = CONVERT( LINE(1:5) )
            IF ( LM .LT. 0 .OR. LM. GT. 99999999 ) GO TO 902
          ENDIF
          GO TO 905
 902      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 901
C
C  b) upper bound
C
 905      CONTINUE
          WRITE ( STDOUT, 6005 )
 903      CONTINUE
          WRITE ( STDOUT, 6006 )
          READ ( STDIN, '( A80 ) ', ERR = 904 ) LINE
          CHOICE = UPPER( LINE(1:1) )
          IF ( CHOICE .NE. ' ' ) THEN
            UM = CONVERT( LINE(1:5) )
            IF ( UM .LT. 0 .OR. UM. GT. 99999999 .OR. LM. GT. UM )
     1          GO TO 904
          ENDIF
          GO TO 906
 904      CONTINUE
          WRITE ( STDOUT, 1001 )
          GO TO 903
C
C  Verify the bounds on the number of constraints
C
 906      CONTINUE
          WRITE ( STDOUT, 6007 ) LM, UM
          WRITE ( STDOUT, 2003 )
          READ  ( STDIN, '( A1 )' ) CHOICE
          IF ( UPPER( CHOICE ) .EQ. 'Y' ) GO TO 900
        ENDIF
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
      READ ( CLSDVC, '( A28 )', END = 3010 ) PBCLS
      IF ( MATCH( PBCLS(10:28), TARGET, MAXTRG, ANYFNV, ANYFNC,
     1             NVAR, NCON, LN, UN, LM, UM ) ) THEN
        NMATCH = NMATCH + 1
        L = L + 1
        LIST(L) = PBCLS(1:8)
        IF ( L .EQ. 5 ) THEN
          WRITE ( STDOUT, 4002 ) ( LIST(I), I = 1, 5 )
          L = 0
        ENDIF
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
        READ ( CLSDVC, '( A28 )', END = 3050 ) PBCLS
        IF ( MATCH( PBCLS(10:28), TARGET, MAXTRG, ANYFNV, ANYFNC,
     1               NVAR, NCON, LN, UN, LM, UM ) ) THEN
          NBLK = 8
          DO 3030 NLBK = 8, 2, -1
            IF ( PBCLS(NBLK:NBLK) .NE. ' ' ) GO TO 3040
 3030     CONTINUE
 3040     CONTINUE
          WRITE ( FLSDVC, '(A)' ) PBCLS(1:NBLK)
        ENDIF
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
     1               ' type'
     1        /'     R   : Regularity              I : Problem',
     1               ' interest'
     1        /'     N   : Number of variables     M : Number of',
     1               ' constraints'
     1        /'     D   : Degree of available analytic derivatives'
     1        /'     S   : Presence of explicit internal variables'
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
 5009 FORMAT(  '     Number of variables             : v' )
 5010 FORMAT(  '     Number of variables             : *' )
 5011 FORMAT(  '     Number of variables             : any fixed number
     1               ')
 5012 FORMAT(  '     Number of variables             : ',
     1               10 ( I5, 1X ) )
 5013 FORMAT(  '     Number of constraints           : v' )
 5014 FORMAT(  '     Number of constraints           : *' )
 5015 FORMAT(  '     Number of constraints           : any fixed number
     1               ')
 5016 FORMAT(  '     Number of constraints           : ',
     1               10 ( I5, 1X ) )
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
 1014 FORMAT( /'   NUMBER OF VARIABLES :'
     1        /'   ---------------------' )
 1015 FORMAT( /'     F   : Fixed                   V : Variable'
     1        /'     I   : In an interval'
     1        /'    <CR> : Any number of variables (*)'
     1        /' '
     1        /'   Your choice :' )
 1016 FORMAT( /'   NUMBER OF CONSTRAINTS :'
     1        /'   -----------------------' )
 1017 FORMAT( /'     F   : Fixed                   V : Variable'
     1        /'     I   : In an interval'
     1        /'    <CR> : Any number of constraints (*)'
     1        /' '
     1        /'   Your choice :' )
 1018 FORMAT( /'   SELECT A NUMBER OF VARIABLES:'
     1        /'   -----------------------------' )
 1119 FORMAT( /'    (INT) : Select only problems with (INT) variables'
     1        /'            (minimum 0, maximum 99999999, multiple ',
     1               'choices are allowed)'
     1        /'    <CR>  : Any fixed number of variables (*)'
     1        /' '
     1        /'   Your choice :' )
 1019 FORMAT( /'    (INT) : Select only problems with (INT) variables'
     1        /'            (minimum 0, maximum 99999999, multiple ',
     1               'choices are allowed)'
     1        /'     *    : Any fixed number of variables'
     1        /'    <CR>  : No further selection',
     1        /' '
     1        /'   Your choice :' )
 1020 FORMAT( /'   SELECT A NUMBER OF CONSTRAINTS:'
     1        /'   -------------------------------' )
 1121 FORMAT( /'    (INT) : Select only problems with (INT) constraints'
     1        /'            (minimum 0, maximum 99999999, multiple ',
     1               'choices are allowed)'
     1        /'    <CR>  : Any fixed number of variables (*)'
     1        /' '
     1        /'   Your choice :' )
 1021 FORMAT( /'    (INT) : Select only problems with (INT) constraints'
     1        /'            (minimum 0, maximum 99999999, multiple ',
     1               'choices are allowed)'
     1        /'     *    : Any fixed number of variables '
     1        /'    <CR>  : No further selection',
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
 2009 FORMAT( /'   You have specified a fixed number of variables.' )
 2019 FORMAT( /'   You have specified a variable number of variables.' )
 2029 FORMAT( /'   You have specified any number of variables.' )
 2010 FORMAT( /'   You have specified a number of variables',
     1         ' in the set: ',
     1        /'       ', 10( I5, 1X ) )
 2020 FORMAT( /'   You have specified any fixed number of variables.' )
 2110 FORMAT( /'   *** Please choose a number of variables.' )
 2011 FORMAT( /'   You have specified a fixed number of constraints.' )
 2021 FORMAT(/'   You have specified a variable number of constraints.')
 2031 FORMAT( /'   You have specified any number of constraints.' )
 2012 FORMAT( /'   You have specified a number of constraints',
     1         ' in the set: ',
     1        /'       ', 10( I5, 1X ) )
 2022 FORMAT( /'   You have specified any fixed number of constraints.')
 2120 FORMAT( /'   *** Please choose a number of constraints.' )
 4000 FORMAT( /'   ', I5, ' Problem(s) match(es) the specification.' / )
 4001 FORMAT( /'   MATCHING PROBLEMS :'
     1        /'   -------------------' / )
 4002 FORMAT(  '     ', 5( A8, 3X ) )
 4003 FORMAT( /'   LOWER BOUND ON THE NUMBER OF VARIABLES :'
     1        /'   ----------------------------------------' )
 4004 FORMAT( /'    (INT) : Problems with at least (INT) variables'
     1        /'    <CR>  : No lower bound on the number of variables'
     1        /' '
     1        /'   Your choice : ' )
 4005 FORMAT( /'   UPPER BOUND ON THE NUMBER OF VARIABLES :'
     1        /'   ----------------------------------------' )
 4006 FORMAT( /'    (INT) : Problems with at most (INT) variables'
     1        /'    <CR>  : No upper bound on the number of variables'
     1        /' '
     1        /'   Your choice : ' )
 4007 FORMAT( /'   You have specified a number of variables in [ ',
     1        I5,', ',I5,' ]' )
 4008 FORMAT(  '     Number of variables             : in [ ',
     1        I5,', ',I5,' ]' )
 4009 FORMAT( /'   You have specified an interval for the number of',
     1         ' variables.' )
 6003 FORMAT( /'   LOWER BOUND ON THE NUMBER OF CONSTRAINTS :'
     1        /'   -----------------------------------------' )
 6004 FORMAT( /'    (INT) : Problems with at least (INT) constraints'
     1        /'    <CR>  : No lower bound on the number of constraints'
     1        /' '
     1        /'   Your choice : ' )
 6005 FORMAT( /'   UPPER BOUND ON THE NUMBER OF CONSTRAINTS :'
     1        /'   -----------------------------------------' )
 6006 FORMAT( /'    (INT) : Problems with at most (INT) constraints'
     1        /'    <CR>  : No upper bound on the number of constraints'
     1        /' '
     1        /'   Your choice : ' )
 6007 FORMAT( /'   You have specified a number of constraints in [ ',
     1        I5,', ',I5,' ]' )
 6008 FORMAT(  '     Number of constraints           : in [ ',
     1        I5,', ',I5,' ]' )
 6009 FORMAT( /'   You have specified an interval for the number of',
     1         ' constraints.' )
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
      LOGICAL FUNCTION REJECT( CHOICE, ITEM, NUM)
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
C         8 = number of constraints.
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
     1             INTVAR = 6, VARN = 7, VARM = 8 )
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
      ELSE IF ( ITEM .EQ. VARN .OR. ITEM .EQ. VARM ) THEN
        ADMIT = CHOICE .EQ. 'F' .OR. CHOICE .EQ. 'V' .OR.
     1          CHOICE .EQ. 'I'
      ENDIF
      REJECT = .NOT. ADMIT
      RETURN
      END
C
C-----------------------------------------------------------------------------
C
      LOGICAL FUNCTION MATCH ( C, T, MAXTRG, ANYFNV, ANYFNC, NV, NC,
     1                         LN, UN, LM, UM )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to check if the classification C
C  matches one of the targets specified by T, ANYFNV, ANYFNC, NV, NC, LN
C  UN, LM and UM.
C
C  Programming: Ph. Toint, Dec 1991, for CGT Productions.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      CHARACTER*8  T( MAXTRG )
      CHARACTER*19 C
      LOGICAL      ANYFNV, ANYFNC
      INTEGER      NV( MAXTRG ), NC( MAXTRG ), MAXTRG, LN, UN, LM, UM
C
C  Other variables
C
      INTEGER      I, J, CONVERT
      CHARACTER*1  CH
C
C  Match objective function type
C
      MATCH = .FALSE.
      DO 10 I = 1, MAXTRG
        CH = T(I)(1:1)
        IF ( CH .EQ. ' ' ) GO TO 20
        MATCH = MATCH .OR. CH .EQ. C(1:1) .OR. CH .EQ. '*'
 10   CONTINUE
 20   CONTINUE
      IF ( .NOT. MATCH ) RETURN
C
C  Match constraint type
C
      MATCH = .FALSE.
      DO 30 I = 1, MAXTRG
        CH = T(I)(2:2)
        IF ( CH .EQ. ' ' ) GO TO 40
        MATCH = MATCH .OR. CH .EQ. C(2:2) .OR. CH .EQ. '*'
 30   CONTINUE
 40   CONTINUE
      IF ( .NOT. MATCH ) RETURN
C
C  Match regularity type
C
      CH = T(1)(3:3)
      MATCH = CH .EQ. C(3:3) .OR. CH .EQ. '*'
      IF ( .NOT. MATCH ) RETURN
C
C  Match degree of available derivatives
C
      MATCH = .FALSE.
      DO 50 I = 1, MAXTRG
        CH = T(I)(4:4)
        IF ( CH .EQ. ' ' ) GO TO 60
        MATCH = MATCH .OR. CH .EQ. C(4:4) .OR. CH .EQ. '*'
 50   CONTINUE
 60   CONTINUE
      IF ( .NOT. MATCH ) RETURN
C
C  Match interest of the problem
C
      MATCH = .FALSE.
      DO 70 I = 1, MAXTRG
        CH = T(I)(5:5)
        IF ( CH .EQ. ' ' ) GO TO 80
        MATCH = MATCH .OR. CH .EQ. C(6:6) .OR. CH .EQ. '*'
 70   CONTINUE
 80   CONTINUE
      IF ( .NOT. MATCH ) RETURN
C
C  Match for explicit internal variables
C
      CH = T(1)(6:6)
      MATCH = CH .EQ. C(7:7) .OR. CH .EQ. '*'
      IF ( .NOT. MATCH ) RETURN
C
C  Match the number of variables
C
      MATCH = .FALSE.
      CH = T(1)(7:7)
      IF ( CH .EQ. '*' ) THEN
        MATCH = .TRUE.
      ELSE IF ( CH .EQ. 'V' ) THEN
        MATCH = C(13:13) .EQ. 'V'
      ELSE IF ( CH .EQ. 'I' ) THEN
C    interval.  A variable number of variables is no longer considered to
C    match a number in an interval.
        IF ( C(13:13) .NE. 'V' ) THEN
          I = CONVERT( C(9:13) )
          MATCH = I. GE. LN .AND. I .LE. UN
        ENDIF
      ELSE
        IF ( ANYFNV ) THEN
          MATCH = .TRUE.
        ELSE
C    A variable number of variables is no longer considered to match a fixed
C    number.
          IF ( C(13:13) .NE. 'V' ) THEN
            J = CONVERT( C(9:13) )
            DO 90 I = 1, MAXTRG
              IF ( NV(I) .LT. 0 ) GO TO 90
              MATCH = NV(I) .EQ. J
              IF ( MATCH ) GO TO 100
 90         CONTINUE
 100        CONTINUE
          ENDIF
        ENDIF
      ENDIF
      IF ( .NOT. MATCH ) RETURN
C
C  Match the number of constraints
C
      MATCH = .FALSE.
      CH = T(1)(8:8)
      IF ( CH .EQ. '*' ) THEN
        MATCH = .TRUE.
      ELSE IF ( CH .EQ. 'V' ) THEN
        MATCH = C(19:19) .EQ. 'V'
      ELSE IF ( CH .EQ. 'I' ) THEN
C    interval.  A variable number of constraints is no longer considered to
C    match a number in an interval.
        IF ( C(19:19) .NE. 'V' ) THEN
          I = CONVERT( C(15:19) )
          MATCH = I. GE. LM .AND. I .LE. UM
        ENDIF
      ELSE
        IF ( ANYFNC ) THEN
          MATCH = .TRUE.
        ELSE
C    A variable number of constraints is no longer considered to match a fixed
C    number.
          IF ( C(19:19) .NE. 'V' ) THEN
            J = CONVERT( C(15:19) )
            DO 110 I = 1, MAXTRG
              IF ( NC(I) .LT. 0 ) GO TO 110
              MATCH = NC(I) .EQ. J
              IF ( MATCH ) GO TO 120
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
      INTEGER FUNCTION CONVERT( LINE )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this function is to convert the nonnegative integer contained
C  in the first 5 characters of the string LINE into a proper integer.
C  -1 is returned in case of error.
C
C  Programming: Ph. Toint, Dec 1991, for CGT Productions.
C
C-----------------------------------------------------------------------------
C
C  Arguments
C
      CHARACTER*5    LINE
C
C  Other variables
C
      INTEGER L
C
C  Remove trailing blanks
C
      DO 10 L= 5, 1, -1
        IF ( LINE(L:L) .NE. ' ' ) GO TO 20
 10   CONTINUE
      GO TO 30
C
C  Read the integer
C
 20   CONTINUE
      IF ( L .EQ. 1 ) THEN
        READ ( LINE(1:1), '( I1 )', ERR = 30 ) CONVERT
      ELSE IF ( L .EQ. 2 ) THEN
        READ ( LINE(1:2), '( I2 )', ERR = 30 ) CONVERT
      ELSE IF ( L .EQ. 3 ) THEN
        READ ( LINE(1:3), '( I3 )', ERR = 30 ) CONVERT
      ELSE IF ( L .EQ. 4 ) THEN
        READ ( LINE(1:4), '( I4 )', ERR = 30 ) CONVERT
      ELSE
        READ ( LINE(1:5), '( I5 )', ERR = 30 ) CONVERT
      ENDIF
      RETURN
C
C  Error
C
 30   CONTINUE
      CONVERT = -1
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
      CHARACTER*8 T( MX )
      INTEGER     I, MX
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
