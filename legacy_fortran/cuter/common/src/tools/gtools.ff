      SUBROUTINE GETVEC( IPTR, X, NMAX, N )
      INTEGER*4       IPTR, NMAX, N
      REAL*8          X( NMAX )
C
C  Copy input argument pointed to by IPTR into the target vector X.
C  Make sure the input argument is a real*8 vector, and that
C  it doesn't exceed the length of the target vector.
C  Put the length of the input vector in N.
C
C  MEX functions
C
      INTEGER*4       MXGETPR
      INTEGER*4       MXGETM, MXGETN
      INTEGER*4       MXISDOUBLE
C
C  Local variables
C
      INTEGER*4       NROWS, NCOLS, M, PR
      CHARACTER*100   MSG
C
      NROWS = MXGETM( IPTR )
      NCOLS = MXGETN( IPTR )
      M     = MIN( NROWS, NCOLS )
      N     = MAX( NROWS, NCOLS )
      IF ( M .NE. 1 ) THEN
         MSG = 'INPUT ARGUMENT MUST BE A VECTOR.'
         CALL MEXERRMSGTXT( MSG )
      END IF
      IF ( N .GT. NMAX ) THEN
         MSG = 'INPUT ARGUMENT IS TOO LONG.'
         CALL MEXERRMSGTXT( MSG )
      END IF
      IF ( MXISDOUBLE( IPTR ) .EQ. 0 ) THEN
         MSG = 'INPUT ARGUMENT MUST BE REAL*8.'
         CALL MEXERRMSGTXT( MSG )
      END IF
      PR = MXGETPR( IPTR )
      CALL MXCOPYPTRTOREAL8( PR, X, N )
      RETURN
      END
C
      SUBROUTINE CNVSPR( NCOL, NNZSH, SH, IRNSH, ICNSH, A, IROW, JCOL )
      INTEGER*4       NCOL, NNZSH
      INTEGER*4       IRNSH( NNZSH ), ICNSH( NNZSH ), IROW( NNZSH )
      INTEGER*4       JCOL( NCOL + 1 )
      REAL*8          SH( NNZSH ), A( NNZSH )
C
C  Convert the CUTE sparse matrix described by SH, IRNSH and ICNSH
C  to the MATLAB sparse matrix described by A, IROW and JCOL.
C
C  Local variables
C
      INTEGER*4       NCOLP1, I, J, K, JJ
C
C  Zero out JCOL and then use it to store the number of nonzero
C  entries in column J of the sparse matrix.
C
      NCOLP1 = NCOL + 1
      DO 100 J = 1, NCOLP1
         JCOL( J ) = 0
  100 CONTINUE
      DO 110 K = 1, NNZSH
         J         = ICNSH( K )
         JCOL( J ) = JCOL( J ) + 1
  110 CONTINUE
      JCOL( NCOLP1 ) = NNZSH + 1
C
C  Now go backwards through JCOL to find the starting index for the
C  nonzero entries in column J.
C
      DO 120 J = NCOL, 1, -1
         JCOL( J ) = JCOL( J + 1 ) - JCOL( J )
  120 CONTINUE
C
C  Copy the entries from SH and IRNSH into the correct entries
C  in A and IROW.  Since MATLAB is written in C, subtract 1 from
C  each of the entries in IROW.  Use JCOL to keep track of what
C  place we're in for each column.
C
      DO 130 K = 1, NNZSH
         I          = IRNSH( K )
         J          = ICNSH( K )
         JJ         = JCOL ( J )
         A   ( JJ ) = SH   ( K )
         IROW( JJ ) = I  - 1
         JCOL( J  ) = JJ + 1
  130 CONTINUE
C
C  Restore the entries in JCOL to the starting index for the nonzero
C  entries in each column.  Since MATLAB is written in C, subtract 1
C  from each of the entries in JCOL.
C
      DO 140 J = NCOL, 2, -1
         JCOL( J ) = JCOL( J - 1 ) - 1
  140 CONTINUE
      JCOL( 1 ) = 0
      JCOL( NCOLP1 ) = JCOL( NCOLP1 ) - 1
      RETURN
      END
C
      SUBROUTINE CNVSAB( NCOL, NNZSH, SH, IRNSH, ICNSH, NNZA, A, IROWA,
     *                   JCOLA, NNZB, B, IROWB, JCOLB )
      INTEGER*4       NCOL, NNZSH, NNZA, NNZB
      INTEGER*4       IRNSH( NNZSH ), ICNSH( NNZSH ), IROWA( NNZSH )
      INTEGER*4       IROWB( NCOL )
      INTEGER*4       JCOLA( NCOL + 1 ), JCOLB( NCOL + 1 )
      REAL*8          SH( NNZSH ), A( NNZSH ), B( NCOL )
C
C  The CUTE sparse matrix described by SH, IRNSH and ICNSH contains
C  both a sparse matrix and a sparse vector.  The k-th nonzero entry
C  in SH belongs to the sparse vector if IRNSH(k)=0, and belongs
C  to row i of the sparse matrix if IRNSH(k)=i>0.
C  This subroutine returns the MATLAB sparse matrix described by
C  NNZA, A, IROWA, JCOLA, and the MATLAB sparse vector described by
C  NNZB, B, IROWB, JCOLB.
C
C  Local variables
C
      INTEGER*4       NCOLP1, I, J, K, JJ
C
C  Zero out JCOLA and JCOLB.  Use JCOLA to store the number of nonzero
C  entries in column J of the sparse matrix, and use JCOLB to store
C  the number of nonzero entries in column J of the sparse vector.
C
      NCOLP1 = NCOL + 1
      DO 100 J = 1, NCOLP1
         JCOLA( J ) = 0
         JCOLB( J ) = 0
  100 CONTINUE
      NNZA = 0
      NNZB = 0
      DO 110 K = 1, NNZSH
         I         = IRNSH( K )
         J         = ICNSH( K )
         IF ( I .GT. 0 ) THEN
            JCOLA( J ) = JCOLA( J ) + 1
            NNZA       = NNZA + 1
         ELSE
            JCOLB( J ) = JCOLB( J ) + 1
            NNZB       = NNZB + 1
         END IF
  110 CONTINUE
      JCOLA( NCOLP1 ) = NNZA + 1
      JCOLB( NCOLP1 ) = NNZB + 1
C
C  Now go backwards through JCOLA and JCOLB to find the starting index
C  for the nonzero entries in column J.
C
      DO 120 J = NCOL, 1, -1
         JCOLA( J ) = JCOLA( J + 1 ) - JCOLA( J )
         JCOLB( J ) = JCOLB( J + 1 ) - JCOLB( J )
  120 CONTINUE
C
C  Copy the entries from SH and IRNSH into the correct entries
C  in A and IROWA, or B if IRNSH(k)=0.  Since MATLAB is written
C  in C, subtract 1 from each of the entries in IROWA.  Set all
C  of the entries in IROWB=0.  Use JCOLA and JCOLB to keep track
C  of what place we're in for each column.
C
      DO 130 K = 1, NNZSH
         I          = IRNSH( K )
         J          = ICNSH( K )
         IF ( I .GT. 0 ) THEN
            JJ          = JCOLA( J )
            A    ( JJ ) = SH   ( K )
            IROWA( JJ ) = I  - 1
            JCOLA( J  ) = JJ + 1
         ELSE
            JJ          = JCOLB( J )
            B    ( JJ ) = SH   ( K )
            IROWB( JJ ) = 0
            JCOLB( J  ) = JJ + 1
         END IF
  130 CONTINUE
C
C  Restore the entries in JCOLA and JCOLB to the starting index for
C  the nonzero entries in each column.  Since MATLAB is written in C,
C  subtract 1 from each of the entries in JCOLA and JCOLB.
C
      DO 140 J = NCOL, 2, -1
         JCOLA( J ) = JCOLA( J - 1 ) - 1
         JCOLB( J ) = JCOLB( J - 1 ) - 1
  140 CONTINUE
      JCOLA( 1 ) = 0
      JCOLB( 1 ) = 0
      JCOLA( NCOLP1 ) = JCOLA( NCOLP1 ) - 1
      JCOLB( NCOLP1 ) = JCOLB( NCOLP1 ) - 1
      RETURN
      END
C
      SUBROUTINE CNVCHR( CNAMES, NWORD, LWORD, WORK )
      CHARACTER*(*) CNAMES( NWORD )
      INTEGER*4     NWORD, LWORD
      REAL*8        WORK( NWORD, LWORD )
C
C  Convert individual character entries in character array CNAMES
C  to the floating point values of their position in the collating
C  sequence (i.e., convert to floating point ASCII values).
C  Put these values in the REAL*8 array WORK.
C
      DO 120 I = 1, NWORD
         DO 100 J = 1, LWORD
            WORK( I, J ) = DFLOAT( ICHAR( CNAMES(I)(J:J) ) )
  100    CONTINUE
  120 CONTINUE
      RETURN
      END
C
C     UPPER - CONVERT STRING TO UPPER CASE EXCEPT WHERE ENCLOSED IN QUOTES
C
C     SOURCE:
C     FORTRAN TOOLS FOR VAX/VMS AND MS-DOS
C     BY RUSSELL K. JONES AND TRACY CRABTREE,
C     JOHN WILEY & SONS, 1988.
C
      CHARACTER*(*) FUNCTION UPPER(STRING)
      CHARACTER*(*) STRING
      CHARACTER*1   NULL, QUOTE
      INTEGER       I
C
      NULL = CHAR(0)
      QUOTE = CHAR(39)
      I = 1
      DO WHILE (STRING(I:I) .NE. NULL)
         IF (STRING(I:I) .EQ. QUOTE) THEN
            UPPER(I:I) = QUOTE
            I = I + 1
            DO WHILE (STRING(I:I) .NE. QUOTE .AND. STRING(I:I) .NE.
     &       NULL)
               UPPER(I:I) = STRING(I:I)
               I = I + 1
            END DO
            UPPER(I:I) = STRING(I:I)
            IF (STRING(I:I) .NE. NULL) I = I + 1
         ELSE
            IF (LLT(STRING(I:I),'a') .OR. LGT(STRING(I:I),'z')) THEN
               UPPER(I:I) = STRING(I:I)
            ELSE
               UPPER(I:I) =
     &          CHAR(ICHAR('A') - ICHAR('a') + ICHAR(STRING(I:I)))
            END IF
            I = I + 1
         END IF
      END DO
      UPPER(I:I) = NULL
      RETURN
      END
C
C     NULLSTR - NULL-TERMINATE AN UNTERMINATED STRING
C
C     SOURCE:
C     FORTRAN TOOLS FOR VAX/VMS AND MS-DOS
C     BY RUSSELL K. JONES AND TRACY CRABTREE,
C     JOHN WILEY & SONS, 1988.
C
      CHARACTER*(*) FUNCTION NULLSTR(STR)
      CHARACTER*(*) STR
      CHARACTER*1   NULL, BLANK
      INTEGER       I
      NULL = CHAR(0)
      BLANK = CHAR(32)
      I = LEN(STR)
      DO WHILE (I .GT. 0 .AND. STR(I:I) .EQ. BLANK)
         I = I - 1
      END DO
      IF (I .EQ. 0) THEN
         NULLSTR = NULL
      ELSE
         NULLSTR = STR(1:I)//NULL
      END IF
      RETURN
      END
C
C     EQUAL - COMPARE TWO STRINGS FOR EQUALITY
C
C     SOURCE:
C     FORTRAN TOOLS FOR VAX/VMS AND MS-DOS
C     BY RUSSELL K. JONES AND TRACY CRABTREE,
C     JOHN WILEY & SONS, 1988.
C
      LOGICAL FUNCTION EQUAL(S,T)
      CHARACTER*(*) S, T
      CHARACTER*1   NULL
      INTEGER       I
      NULL = CHAR(0)
      I = 1
      DO WHILE (S(I:I) .NE. NULL)
         IF (S(I:I) .NE. T(I:I)) THEN
            EQUAL = .FALSE.
            RETURN
         END IF
         I = I + 1
      END DO
      IF (T(I:I) .EQ. NULL) THEN
         EQUAL = .TRUE.
      ELSE
         EQUAL = .FALSE.
      END IF
      RETURN
      END
C
