C     ( Last modified on 10 Aug 2005 at 15:54:33 )
      PROGRAM CLSSFY
C----------------------------------------------------------------------------
C
C  The purpose of this program is to add the classification characteristics
C  of the problem (expressed in SIF format and containing an explicit
C  classification to a database file (CLASSF.DB).  This classification must
C  occur in the file before the GROUPS or VARIABLES section of the SDIF,  and
*  must be of the form
C
C  * classification XXXn-XX-m-m
C
C  or
C
C  * CLASSIFICATION XXXn-XX-m-m
C
C  where n and m are integers containing from 1 to 9 digits.
C  The problem is read on standard input  using default device 4.
C  Other input is read on standard input using default device 5
C  and error messages are written on the error output, default device 0.
C
C  The database file is read (the default name is CLASSF.DB but this can be
C  changed by the user), the new problem classification is then inserted in
C  an updated copy of the database (the default name is CLASSF.UDB but this
C  can also be changed by the user).
C
C  Programming:  A.R. Conn and Ph. Toint, Aug 1992, for CGT Productions.
C
C--------- THE FOLLOWING SPECIFICATIONS MAY BE MODIFIED BY THE USER ----------
C
C  The SIF file is read on standard input using default device 4.
C  Standard input default definition is device 5.  Standard output default
C  definition is device 6. Standard error default definition is device 6.
C  Change them to whatever values are appropriate on your system.
C
      INTEGER      NSTDIN,     STDIN,     STDOUT,     STDERR
C                            (Units may not be 0 in Salford FORTRAN.)
C                            Change by Kristjan Jonasson, March 1993.
      PARAMETER  ( NSTDIN = 4, STDIN = 5, STDOUT = 6, STDERR = 6)

C
C  Classification database file name
C
      CHARACTER*32 DBNAME
      PARAMETER  ( DBNAME = 'CLASSF.DB' )
C
C  Input string for the classification database
C
      INTEGER      DBDVC
      PARAMETER  ( DBDVC = 55 )
C
C  Updated classification database file name
C
      CHARACTER*32 UDBNAM
      PARAMETER  ( UDBNAM = 'CLASSF.UDB' )
C
C  Output string for the updated database
C
      INTEGER      UDBDVC
      PARAMETER  ( UDBDVC = 56 )
C
C  Input string for mode of operation and problem name
C
      INTEGER      MPDVC
      PARAMETER  ( MPDVC = 57 )
      CHARACTER*32 MPNAM
      PARAMETER  ( MPNAM = 'CLSF.DAT' )
C
C---------------- END OF THE USER MODIFIABLE SPECIFICATION ------------------
C
C  Variable definitions
C
      CHARACTER*80 LINE
      CHARACTER*32 IFNAME
      CHARACTER*36 CLS, CLSDB
      CHARACTER*8  NAME, NAMEDB
      CHARACTER*1  CHOICE, OBJ, CON, SMTH, INTRST, INTVAR, UPPER,
     1             VARN, VARM
      LOGICAL      FOUND, WRITTN, AUTO, OVERW
      INTEGER      H, ILEN, IDER, K, NCON, NEWH, NVAR, START, TSIZE
C
C  General initializations
C
      FOUND  = .FALSE.
      AUTO   = .FALSE.
      OVERW  = .FALSE.
      NAME   = '?'
      OBJ    = '?'
      CON    = '?'
      SMTH   = '?'
      INTRST = '?'
      INTVAR = '?'
      VARN   = '?'
      VARM   = '?'
      IDER   = 0
      NVAR   = 0
      NCON   = 0
C
C  Open the mode and problem file
C
      OPEN ( UNIT = MPDVC, FILE = MPNAM, STATUS = 'UNKNOWN' )
C
C  Check for automatic execution
C
      READ  ( MPDVC, '( A1 )' ) CHOICE
      IF ( UPPER( CHOICE ) .EQ. 'Y' ) AUTO = .TRUE.
C
C  Input the initial name of the sif file with the classification you want
C  to add
C
 2    CONTINUE
      READ  ( MPDVC , FMT = '( A )', ERR = 4 ) IFNAME
      CLOSE ( MPDVC )
C                            Inspired by Kristjan Jonasson, March 1993.
C Since sed (or basename) does not exist in DOS, the following
C 9 lines allow that IFNAME be specified with or without the .SIF
C extension. This means that CLASSIFY.BAT and CLASSALL.BAT work.
      TSIZE=INDEX(IFNAME,' ')
      IF (TSIZE .GT. 4 .AND. IFNAME(TSIZE-4:TSIZE-1) .EQ. '.SIF') THEN
        TSIZE=TSIZE-4
      ELSE
        IFNAME(TSIZE:TSIZE) = '.'
        IFNAME(TSIZE+1:TSIZE+1) = 'S'
        IFNAME(TSIZE+2:TSIZE+2) = 'I'
        IFNAME(TSIZE+3:TSIZE+3) = 'F'
      ENDIF
      WRITE ( STDOUT, 6000 ) IFNAME(1:TSIZE+3)
      GO TO 100
 4    CONTINUE
      WRITE ( STDOUT, 1201 )
      GO TO 2
C
C  read the next line if not already at the end of the GROUPS section
C
100   CONTINUE
      OPEN ( UNIT = NSTDIN, FILE = IFNAME, STATUS = 'OLD' )
101   CONTINUE
      IF ( FOUND ) GO TO 200
      READ( NSTDIN, '(A80)', END = 99) LINE
C
C  exit when the GROUPS or VARIABLES section is encountered
C
      IF ( LINE(1:6) .EQ. 'GROUPS'     .OR.
     1     LINE(1:7) .EQ. 'COLUMNS'    .OR.
     1     LINE(1:9) .EQ. 'VARIABLES'  .OR.
     1     LINE(1:4) .EQ. 'ROWS'            ) GO TO 99
C
C  remove trailing blanks and compute line length
C
      DO 10 ILEN = 80, 1, -1
        IF ( LINE(ILEN:ILEN) .NE. ' ' ) GO TO 20
10    CONTINUE
      ILEN = 0
20    CONTINUE
C
C  find the problem's name
C
      IF ( LINE(1:4) .EQ. 'NAME' ) NAME = LINE(15:24)
C
C  skip non-comment lines
C
      IF ( LINE(1:1) .EQ. '*' ) THEN
C
C  find the position of the first word in the comment line
C
       DO 30 K = 2, ILEN
          IF ( LINE(K:K) .NE. ' ' ) GO TO 40
 30    CONTINUE
        K = 1
 40     CONTINUE
C
C  is it the classification ?
C
        IF ( K .GT. 1  .AND. LINE(K:K+13) .EQ. 'CLASSIFICATION' .OR.
     1       LINE(K:K+13) .EQ. 'classification' ) THEN
          FOUND = .TRUE.
C
C  find the classification string
C
          DO 50 START = K + 15, ILEN
             IF ( LINE(START:START) .NE. ' ' ) GO TO 60
 50       CONTINUE
          GO TO 300
 60       CONTINUE
C
C  objective classification
C
          OBJ  = LINE(START:START)
          IF ( OBJ .NE. 'N' .AND.
     1         OBJ .NE. 'C' .AND.
     1         OBJ .NE. 'L' .AND.
     1         OBJ .NE. 'Q' .AND.
     1         OBJ .NE. 'S' .AND.
     1         OBJ. NE. 'O'      ) THEN
            WRITE ( STDERR, 1001 ) NAME
          ENDIF
C
C  constraints classification
C
          CON  = LINE(START+1:START+1)
          IF ( CON .NE. 'U' .AND.
     1         CON .NE. 'X' .AND.
     1         CON .NE. 'B' .AND.
     1         CON .NE. 'N' .AND.
     1         CON .NE. 'L' .AND.
     1         CON .NE. 'Q' .AND.
     1         CON. NE. 'O'      ) THEN
            WRITE ( STDERR, 1002 ) NAME
          ENDIF
C
C  smoothness classification
C
          SMTH = LINE(START+2:START+2)
          IF ( SMTH .NE. 'R' .AND. SMTH .NE. 'I' ) THEN
            WRITE ( STDERR, 1003 ) NAME
          ENDIF
C
C  derivatives degree
C
          CALL GETINT( STDERR, LINE, START + 3, H, IDER )
          IF ( IDER. NE. 0 .AND. IDER .NE. 1 .AND. IDER .NE. 2 ) THEN
            WRITE ( STDERR, 1004 ) NAME
          ENDIF
          IF ( H .EQ. 80 ) THEN
            GO TO 300
          ENDIF
C
C  interest of the problem
C
          INTRST = LINE(H+1:H+1)
          IF ( INTRST .NE. 'A' .AND.
     1         INTRST .NE. 'M' .AND.
     1         INTRST .NE. 'R'      ) THEN
            WRITE ( STDERR, 1005 ) NAME
          ENDIF
C
C  internal variables?
C
          INTVAR = LINE(H+2:H+2)
          IF ( INTVAR .NE. 'Y' .AND. INTVAR .NE. 'N' ) THEN
            WRITE ( STDERR, 1006 ) NAME
          ENDIF
C
C  number of variables
C
          VARN = LINE(H+4:H+4)
          IF ( VARN .NE. 'V' ) THEN
            CALL GETINT( STDERR, LINE, H + 4, NEWH, NVAR )
            H = NEWH
            IF ( NVAR .LE. 0 ) THEN
              WRITE ( STDERR, 1007 ) NAME
              NVAR = 0
            ENDIF
            VARN = 'F'
            IF ( NEWH .EQ. 80 ) THEN
              GO TO 300
            ENDIF
          ELSE
            NVAR = 0
            H = H + 5
          ENDIF
C
C  number of constraints
C
          VARM = LINE(H+1:H+1)
          IF ( VARM .NE. 'V' ) THEN
            CALL GETINT( STDERR, LINE, H + 1, NEWH, NCON )
            IF ( NCON .LT. 0 ) THEN
              WRITE ( STDERR, 1008 ) NAME
              NCON = 0
            ENDIF
            VARM = 'F'
          ELSE
            NCON = 0
          ENDIF
        ENDIF
      ENDIF
C
C  loop for next line
C
      GO TO 101
C
C  incomplete classification
C
 300  CONTINUE
      WRITE ( STDERR, 1009) NAME
C
C  end of file or next SIF section encountered: missing classification
C
  99  CONTINUE
      WRITE ( STDERR, 1010 ) NAME
      STOP
C
C  format the new classification
C
 200  CONTINUE
      CLS = '                            '
      CLS(1:8) = NAME
      CLS(10:10) = OBJ
      CLS(11:11) = CON
      CLS(12:12) = SMTH
      WRITE ( CLS(13:13), '( I1 )' ) IDER
      CLS(14:14) = '-'
      CLS(15:15) = INTRST
      CLS(16:16) = INTVAR
      CLS(17:17) = '-'
      IF ( VARN .EQ. 'V' ) THEN
        CLS(26:26) = VARN
      ELSE
        WRITE ( CLS(18:26), '( I9 )' ) NVAR
      ENDIF
      CLS(27:27) = '-'
      IF ( VARM .EQ. 'V' ) THEN
        CLS(36:36) = VARM
      ELSE
        WRITE ( CLS(28:36), '( I9 )' ) NCON
      ENDIF
C
C  open the classification database
C
      OPEN ( UNIT = DBDVC, FILE = DBNAME, STATUS = 'UNKNOWN' )
C
C  open the updated classification database, according to the
C  conventions of the local operating system
C
      OPEN ( UNIT = UDBDVC, FILE = UDBNAM, STATUS = 'UNKNOWN' )
C
C  read the current line in the existing database
C
      WRITTN = .FALSE.
 400  CONTINUE
      READ ( DBDVC, '( A36 )', END = 900 ) CLSDB
      NAMEDB = CLSDB(1:8)
C
C  verify that the name of the problem to be classified is not already
C  present in the database
C
      IF ( NAMEDB .EQ. NAME ) THEN
        WRITE( STDOUT, 2000 ) NAMEDB
        IF ( .NOT. AUTO ) THEN
          IF ( CLSDB .EQ. CLS ) THEN
             WRITE( STDOUT, 2003 ) DBNAME
          ELSE
             WRITE( STDOUT, 2001 )
             READ  ( STDIN, '( A1 )' ) CHOICE
             IF ( UPPER( CHOICE ) .EQ. 'Y' ) THEN
               OVERW = .TRUE.
             ENDIF
          ENDIF
        ENDIF
      ENDIF
C
C  write the new problem classification, if it is not already present
C  in the existing database or if overwrite is specified
C
      IF ( WRITTN .OR. NAMEDB .LT. NAME ) THEN
        WRITE ( UDBDVC, '( A36 )' ) CLSDB
      ELSEIF ( NAMEDB .EQ. NAME ) THEN
        IF ( CLSDB .GT. CLS .AND. .NOT. OVERW) THEN
          WRITE ( UDBDVC, '( A36 )' ) CLS
          WRITTN = .TRUE.
          WRITE ( UDBDVC, '( A36 )' ) CLSDB
        ELSEIF ( CLSDB .LT. CLS .AND. .NOT. OVERW) THEN
          WRITE ( UDBDVC, '( A36 )' ) CLSDB
          WRITE ( UDBDVC, '( A36 )' ) CLS
          WRITTN = .TRUE.
        ELSE
          WRITE ( UDBDVC, '( A36 )' ) CLS
          WRITTN = .TRUE.
        ENDIF
      ELSE
        WRITE ( UDBDVC, '( A36 )' ) CLS
        WRITE ( UDBDVC, '( A36 )' ) CLSDB
        WRITTN = .TRUE.
      ENDIF
C
C  read next line of the existing database
C
      GO TO 400
C
C  close the database files
C
 900  CONTINUE
      CLOSE ( DBDVC )
      IF ( .NOT. WRITTN ) WRITE ( UDBDVC, '( A36 )' ) CLS
      CLOSE ( UDBDVC )
      STOP
C
C  non excutable statements
C
 1001 FORMAT( ' Problem ', A8,
     1        ' has impossible objective classification!' )
 1002 FORMAT( ' Problem ', A8,
     1        ' has impossible constraints classification!' )
 1003 FORMAT( ' Problem ', A8,
     1        ' has impossible smoothnes classification' )
 1004 FORMAT( ' Problem ', A8,
     1        ' has impossible degree of available derivatives!' )
 1005 FORMAT( ' Problem ', A8,
     1        ' has impossible problem interest specification!' )
 1006 FORMAT( ' Problem ', A8,
     1        ' has  internal variable specification!' )
 1007 FORMAT( ' Problem ', A8,
     1        ' has impossible number of variables!' )
 1008 FORMAT( ' Problem ', A8,
     1        ' has impossible number of constraints!' )
 1009 FORMAT( ' Incomplete classification for problem ', A8 )
 1010 FORMAT( ' Classification not found for Problem ', A8 )
 1011 FORMAT( ' ', A8, ' ', 60A1 )
 1201 FORMAT( '   *** YOU USED MORE THAN 28 CHARACTERS.  Please choose',
     1             ' again.'
     1        / )
 2000 FORMAT( ' Duplicate name: ', A8 )
 2001 FORMAT( ' Do you wish to overwrite the existing classification?',
     1        '  [<CR> = N] ? (N/Y)' )
 2003 FORMAT( ' Identical classification record already appears in ',
     1                A32 )
 6000 FORMAT( /' Your current SIF filename file is : ',
     1                A )
C
      END
C
C
C
      SUBROUTINE GETINT ( STDERR, LINE, START, HYPHEN, I )
C
C-----------------------------------------------------------------------------
C
C  The purpose of this subroutine is to read, in the string LINE, an
C  integer number, from the position START to the first following hyphen
C  or EOL. The position of the hyphen or EOL is returned in HYPHEN, while
C  the integer is returned in I.
C
C  Programming: Ph. Toint, Dec 1991, for CGT Productions.
C
C-----------------------------------------------------------------------------
C
      CHARACTER*80 LINE
      INTEGER      START, HYPHEN, I, STDERR, END
C
C   find the position of the hyphen
C
      END = START - 1
      DO 10 HYPHEN = START, 80
        IF ( LINE( HYPHEN:HYPHEN ) .EQ. '-' ) GO TO 20
        IF ( LINE( HYPHEN:HYPHEN ) .NE. ' ' ) END = END + 1
 10   CONTINUE
 20   CONTINUE
C
C   read the integer
C
      IF ( END .EQ. START ) THEN
        READ ( LINE( START:END ), FMT= '( I1 )', ERR = 30 ) I
      ELSE IF ( END .EQ. START + 1 ) THEN
        READ ( LINE( START:END ), FMT= '( I2 )', ERR = 30 ) I
      ELSE IF ( END .EQ. START + 2 ) THEN
        READ ( LINE( START:END ), FMT= '( I3 )', ERR = 30 ) I
      ELSE IF ( END .EQ. START + 3 ) THEN
        READ ( LINE( START:END ), FMT= '( I4 )', ERR = 30 ) I
      ELSE IF ( END .EQ. START + 4 ) THEN
        READ ( LINE( START:END ), FMT= '( I5 )', ERR = 30 ) I
      ELSE IF ( END .EQ. START + 5 ) THEN
        READ ( LINE( START:END ), FMT= '( I6 )', ERR = 30 ) I
      ELSE IF ( END .EQ. START + 6 ) THEN
        READ ( LINE( START:END ), FMT= '( I7 )', ERR = 30 ) I
      ELSE IF ( END .EQ. START + 7 ) THEN
        READ ( LINE( START:END ), FMT= '( I8 )', ERR = 30 ) I
      ELSE 
        READ ( LINE( START:END ), FMT= '( I9 )', ERR = 30 ) I
      ENDIF
      RETURN
C
C   error in reading the integer
C
 30   CONTINUE
      WRITE ( STDERR, 1000 )
      RETURN
C
C   non executable statements
C
 1000 FORMAT(' Error in reading an integer!' )
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
