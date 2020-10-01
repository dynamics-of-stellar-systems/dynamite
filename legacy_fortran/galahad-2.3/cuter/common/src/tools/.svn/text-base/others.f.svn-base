C     ( Last modified on 23 Dec 2000 at 22:01:38 )
C
C
C
CS    SUBROUTINE SGELIM( M , N , IPVT  , JCOL  , A      )
CD    SUBROUTINE DGELIM( M , N , IPVT  , JCOL  , A      )
      INTEGER            M , N
      INTEGER            IPVT  ( M    ), JCOL  ( N     )
CS    REAL               A     ( M     , N     )
CD    DOUBLE PRECISION   A     ( M     , N     )
C
C  Perform the first M steps of Gaussian Elimination with
C  complete pivoting on the M by N ( M <= N) matrix A.
C
C  Nick Gould, 23rd September 1991.
C  For CGT productions.
C
      INTEGER            I , J , K     , IPIVOT, JPIVOT
CS    REAL               APIVOT, ONE   , ATEMP
CD    DOUBLE PRECISION   APIVOT, ONE   , ATEMP
CS    PARAMETER        ( ONE = 1.0E+0 )
CD    PARAMETER        ( ONE = 1.0D+0 )
C
C  Initialize the column indices.
C
      DO 10 J = 1, N
         JCOL( J ) = J
   10 CONTINUE
C
C  Main loop.
C
      DO 100 K = 1, M
C
C  Compute the K-th pivot.
C
         APIVOT = - ONE
         DO 30 J = K, N
            DO 20 I = K, M
               IF ( ABS( A( I, J ) ) .GT. APIVOT ) THEN
                  APIVOT = ABS( A( I, J ) )
                  IPIVOT = I
                  JPIVOT = J
               END IF
   20       CONTINUE
   30    CONTINUE
C
C  Interchange rows I and IPIVOT.
C
         IPVT( K ) = IPIVOT
         IF ( IPIVOT .GT. K ) THEN
            DO 40 J   = K, N
               ATEMP          = A( IPIVOT, J )
               A( IPIVOT, J ) = A( K     , J )
               A( K     , J ) = ATEMP
   40       CONTINUE
         END IF
C
C  Interchange columns J and JPIVOT.
C
         IF ( JPIVOT .GT. K ) THEN
            J              = JCOL( JPIVOT )
            JCOL( JPIVOT ) = JCOL( K      )
            JCOL( K      ) = J
            DO 50 I = 1, M
               ATEMP          = A( I, JPIVOT )
               A( I, JPIVOT ) = A( I, K )
               A( I, K      ) = ATEMP
   50       CONTINUE
         END IF
C
C  Perform the elimination.
C
         APIVOT       = A( K, K )
         DO 70 I      = K + 1, M
            ATEMP     = A( I, K ) / APIVOT
            A( I, K ) = ATEMP
            DO 60 J      = K + 1, N
               A( I, J ) = A( I, J ) - ATEMP * A( K, J )
   60       CONTINUE
   70    CONTINUE
  100 CONTINUE
      RETURN
C
C  End of subroutine GELIM.
C
      END
C
C
C
CS    SUBROUTINE SGESLV( M     , IPVT  , A , X  )
CD    SUBROUTINE DGESLV( M     , IPVT  , A , X  )
      INTEGER            M
      INTEGER            IPVT  ( M    )
CS    REAL               A     ( M     , M     ), X    ( M       )
CD    DOUBLE PRECISION   A     ( M     , M     ), X    ( M       )
C
C  Solve the equations A(T)x = b. The vector b is input in X.
C  The LU factors of P A are input in A; The permutation P is stored
C  in IPVT. The solution x is output in X.
C
C  Nick Gould, 23rd September 1991.
C  For CGT productions.
C
      INTEGER            I , K
CS    REAL               XTEMP, ZERO
CD    DOUBLE PRECISION   XTEMP, ZERO
CS    PARAMETER        ( ZERO = 0.0E+0 )
CD    PARAMETER        ( ZERO = 0.0D+0 )
C
C  Solve U(T)y = b. The vector b is input in X; y is output in X.
C
      DO 20 K  = 1, M
         XTEMP = ZERO
         DO 10 I   = 1, K - 1
            XTEMP  = XTEMP + A( I, K ) * X( I )
   10    CONTINUE
         X( K ) = ( X( K ) - XTEMP ) / A( K, K )
   20 CONTINUE
C
C  Solve L(T) x = y. The vector y is input in X; x is output in X.
C
      DO 40 K  = M - 1, 1, - 1
         XTEMP = ZERO
         DO 30 I   = K + 1, M
            XTEMP  = XTEMP + A( I, K ) * X( I )
   30    CONTINUE
         X( K ) = X( K ) - XTEMP
         I      = IPVT( K )
         IF ( I .NE. K ) THEN
            XTEMP  = X( I )
            X( I ) = X( K )
            X( K ) = XTEMP
         END IF
   40 CONTINUE
      RETURN
C
C  End of subroutine GESLV.
C
      END
C
C
C  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
C     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.
C
CS    SUBROUTINE SSYMMH( MAXSZH, ISYMMH, ISYMMD )
CD    SUBROUTINE DSYMMH( MAXSZH, ISYMMH, ISYMMD )
      INTEGER MAXSZH
      INTEGER ISYMMH( MAXSZH, MAXSZH ), ISYMMD( MAXSZH )
C
C  GIVEN A COLUMNWISE STORAGE SCHEME OF THE UPPER TRIANGLE OF A
C  SYMMETRIC MATRIX OF ORDER MAXSZH, COMPUTE THE POSITION OF THE
C  I,J-TH ENTRY OF THE SYMMETRIC MATRIX IN THIS SCHEME.
C
C  THE VALUE ISYMMH( I, J ) + 1 GIVES THE POSITION OF THE I,J-TH
C  ENTRY OF THE MATRIX IN THE UPPER TRIANGULAR SCHEME.
C
C  NICK GOULD, 10TH OF MAY 1989.
C  FOR CGT PRODUCTIONS.
C
      INTEGER I, J, K
      K       = 0
      DO 20 J = 1, MAXSZH
CDIR$ IVDEP
         DO 10 I           = 1, J - 1
            ISYMMH( I, J ) = K
            ISYMMH( J, I ) = K
            K              = K + 1
   10    CONTINUE
         ISYMMD( J )    = K
         ISYMMH( J, J ) = K
         K              = K + 1
   20 CONTINUE
      RETURN
C
C  END OF SYMMH.
C
      END
C
C
C
CS    SUBROUTINE SSETVL( N, X, INCX, VL )
CD    SUBROUTINE DSETVL( N, X, INCX, VL )
C
C     ******************************************************************
C
      INTEGER          INCX, N
CS    REAL             VL, X( * )
CD    DOUBLE PRECISION VL, X( * )
      INTEGER          IX, M, I, J
      IF ( N .LE. 0 ) RETURN
      IF ( INCX .NE. 1 ) THEN
         IF( INCX .LT. 0 ) THEN
             IX = ( - N + 1 ) * INCX + 1
         ELSE
             IX = 1
         END IF
         J         = IX + ( N - 1 ) * INCX
         DO 100 I  = IX, J, INCX
            X( I ) = VL
  100    CONTINUE
      ELSE
         M = MOD( N, 5 )
         IF ( M .NE. 0 ) THEN
            DO 200 I  = 1, M
               X( I ) = VL
  200       CONTINUE
         END IF
         IF ( N .GE. 5 ) THEN
            IX            = M + 1
            DO 300 I      = IX, N, 5
               X( I )     = VL
               X( I + 1 ) = VL
               X( I + 2 ) = VL
               X( I + 3 ) = VL
               X( I + 4 ) = VL
  300       CONTINUE
         END IF
      END IF
      RETURN
      END
C
C
C
CS    SUBROUTINE SSETVI( NVAR, X, IVAR, VL )
CD    SUBROUTINE DSETVI( NVAR, X, IVAR, VL )
C
C     ******************************************************************
C
      INTEGER          NVAR, IVAR( * ), I
CS    REAL             VL, X( * )
CD    DOUBLE PRECISION VL, X( * )
      IF ( NVAR .LE. 0 ) RETURN
      DO 10 I           = 1, NVAR
         X( IVAR( I ) ) = VL
   10 CONTINUE
      RETURN
      END
