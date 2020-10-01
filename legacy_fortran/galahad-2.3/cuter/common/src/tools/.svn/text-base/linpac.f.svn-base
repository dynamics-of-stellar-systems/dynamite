C     ( Last modified on 23 Dec 2000 at 22:01:38 )
CS    SUBROUTINE SCOPY(N,DX,INCX,DY,INCY)
CD    SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
CS    REAL             DX(*),DY(*)
CD    DOUBLE PRECISION DX(*),DY(*)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END
C
C
C
CS    REAL             FUNCTION SDOT(N,DX,INCX,DY,INCY)
CD    DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
CS    REAL             DX(*),DY(*),DTEMP,ZERO
CD    DOUBLE PRECISION DX(*),DY(*),DTEMP,ZERO
CS    PARAMETER ( ZERO = 0.0E+0 )
CD    PARAMETER ( ZERO = 0.0D+0 )
C
CS    SDOT = ZERO
CD    DDOT = ZERO
      DTEMP = ZERO
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
CS    SDOT = DTEMP
CD    DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 CONTINUE
CS    SDOT = DTEMP
CD    DDOT = DTEMP
      RETURN
      END
C
C
C
CS    REAL             FUNCTION SNRM2 ( N, DX, INCX)
CD    DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER          N, INCX, NEXT
CS    REAL             DX( * ), CUTLO, CUTHI, HITEST, SUM,
CD    DOUBLE PRECISION DX( * ), CUTLO, CUTHI, HITEST, SUM,
     *                 XMAX, ZERO, ONE
      INTRINSIC        ABS, SQRT
      INTEGER          I, J, NN
CS    PARAMETER ( ZERO = 0.0E+0, ONE = 1.0E+0 )
CD    PARAMETER ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
CS    DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
CD    DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
      IF ( N .LE. 0) THEN
CS       SNRM2  = ZERO
CD       DNRM2  = ZERO
      ELSE
         NEXT = 1
         SUM  = ZERO
         NN   = N * INCX
C
C  BEGIN MAIN LOOP
C
         I = 1
   20    CONTINUE
         GO TO ( 30, 50, 70, 110 ), NEXT
   30    CONTINUE
         IF( ABS( DX( I ) ) .GT. CUTLO ) GO TO 85
         NEXT = 2
         XMAX = ZERO
C
C  PHASE 1.  SUM IS ZERO
C
   50    CONTINUE
         IF ( DX( I ) .EQ. ZERO ) GO TO 200
         IF ( ABS( DX( I ) ) .GT. CUTLO ) GO TO 85
C
C  PREPARE FOR PHASE 2.
C
         NEXT = 3
         GO TO 105
C
C  PREPARE FOR PHASE 4.
C
  100    CONTINUE
         I    = J
         NEXT = 4
         SUM  = ( SUM / DX( I ) ) / DX( I )
  105    CONTINUE
         XMAX = ABS( DX( I ) )
         SUM  = SUM + ( DX( I ) / XMAX ) ** 2
         GO TO 200
C
C  PHASE 2.  SUM IS SMALL. SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70    CONTINUE
         IF ( ABS( DX( I ) ) .GT. CUTLO ) THEN
C
C  PREPARE FOR PHASE 3.
C
            SUM = ( SUM * XMAX) * XMAX
            GO TO 85
         END IF
C
C  COMMON CODE FOR PHASES 2 AND 4.
C  IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110    CONTINUE
         IF ( ABS( DX( I ) ) .GT. XMAX ) THEN
            SUM  = ONE + SUM * ( XMAX / DX( I ) ) ** 2
            XMAX = ABS( DX( I ) )
         ELSE
            SUM = SUM + ( DX( I ) / XMAX ) ** 2
         END IF
  200    CONTINUE
         I = I + INCX
         IF ( I .LE. NN ) GO TO 20
C
C  END OF MAIN LOOP.
C
C  COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
CS       SNRM2 = XMAX * SQRT(SUM)
CD       DNRM2 = XMAX * SQRT(SUM)
         GO TO 300
C
C  FOR REAL OR D.P. SET HITEST = CUTHI/N
C
   85    CONTINUE
         HITEST = CUTHI/FLOAT( N )
C
C  PHASE 3. SUM IS MID-RANGE.  NO SCALING.
C
         DO 95 J = I, NN, INCX
            IF( ABS( DX( J ) ) .GE. HITEST ) GO TO 100
            SUM = SUM + DX( J ) ** 2
   95    CONTINUE
CS       SNRM2 = SQRT( SUM )
CD       DNRM2 = SQRT( SUM )
      END IF
  300 CONTINUE
      RETURN
      END
C
C
C
CS    SUBROUTINE SAXPY(N,DA,DX,INCX,DY,INCY)
CD    SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
CS    REAL             DX(*),DY(*),DA,ZERO
CD    DOUBLE PRECISION DX(*),DY(*),DA,ZERO
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0 )
C
      IF(N.LE.0)RETURN
      IF (DA .EQ. ZERO) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
C
C
C
CS    SUBROUTINE SROT (N,DX,INCX,DY,INCY,C,S)
CD    SUBROUTINE DROT (N,DX,INCX,DY,INCY,C,S)
C
C     APPLIES A PLANE ROTATION.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
CS    REAL             DX(*),DY(*),DTEMP,C,S
CD    DOUBLE PRECISION DX(*),DY(*),DTEMP,C,S
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = C*DX(IX) + S*DY(IY)
        DY(IY) = C*DY(IY) - S*DX(IX)
        DX(IX) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        DTEMP = C*DX(I) + S*DY(I)
        DY(I) = C*DY(I) - S*DX(I)
        DX(I) = DTEMP
   30 CONTINUE
      RETURN
      END
CS    SUBROUTINE SROTG(DA,DB,C,S)
CD    SUBROUTINE DROTG(DA,DB,C,S)
C
C     CONSTRUCT GIVENS PLANE ROTATION.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
CS    REAL             DA,DB,C,S,ROE,SCALE,R,Z,ZERO,ONE
CD    DOUBLE PRECISION DA,DB,C,S,ROE,SCALE,R,Z,ZERO,ONE
      INTRINSIC ABS, SQRT, SIGN
CS    PARAMETER ( ZERO = 0.0E+0, ONE = 1.0E+0 )
CD    PARAMETER ( ZERO = 0.0D+0, ONE = 1.0D+0 )
C
      ROE = DB
      IF( ABS(DA) .GT. ABS(DB) ) ROE = DA
      SCALE = ABS(DA) + ABS(DB)
      IF( SCALE .NE. ZERO ) GO TO 10
         C = ONE
         S = ZERO
         R = ZERO
         GO TO 20
   10 R = SCALE*SQRT((DA/SCALE)**2 + (DB/SCALE)**2)
      R = SIGN(ONE,ROE)*R
      C = DA/R
      S = DB/R
   20 Z = ONE
      IF( ABS(DA) .GT. ABS(DB) ) Z = S
      IF( ABS(DB) .GE. ABS(DA) .AND. C .NE. ZERO ) Z = ONE/C
      DA = R
      DB = Z
      RETURN
      END
C
C
C
CS    SUBROUTINE  SSCAL(N,DA,DX,INCX)
CD    SUBROUTINE  DSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
CS    REAL             DA,DX(*)
CD    DOUBLE PRECISION DA,DX(*)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
C
C
C
CS    SUBROUTINE  SSWAP (N,DX,INCX,DY,INCY)
CD    SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
CS    REAL             DX(*),DY(*),DTEMP
CD    DOUBLE PRECISION DX(*),DY(*),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
      END
C
C
C
CS    INTEGER FUNCTION ISAMAX(N,DX,INCX)
CD    INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C     finds the index of element having max. absolute value.
C     jack dongarra, linpack, 3/11/78.
C     modified 3/93 to return if incx .le. 0.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
CS    REAL DX(*),DMAX
CD    DOUBLE PRECISION DX(*),DMAX
      INTEGER I,INCX,IX,N
C
CS    ISAMAX = 0
CD    IDAMAX = 0
      IF( N.LT.1 .OR. INCX.LE.0 ) RETURN
CS    ISAMAX = 1
CD    IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        code for increment not equal to 1
C
      IX = 1
      DMAX = ABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(ABS(DX(IX)).LE.DMAX) GO TO 5
CS       ISAMAX = I
CD       IDAMAX = I
         DMAX = ABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        code for increment equal to 1
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
         IF(ABS(DX(I)).LE.DMAX) GO TO 30
CS       ISAMAX = I
CD       IDAMAX = I
         DMAX = ABS(DX(I))
   30 CONTINUE
      RETURN
      END
C
C
C
CS    REAL FUNCTION SASUM(N,DX,INCX)
CD    DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
C
C     takes the sum of the absolute values.
C     jack dongarra, linpack, 3/11/78.
C     modified 3/93 to return if incx .le. 0.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
CS    REAL DX(*),DTEMP
CD    DOUBLE PRECISION DX(*),DTEMP
      INTEGER I,INCX,M,MP1,N,NINCX
C
CS    SASUM = 0.0E0
CD    DASUM = 0.0D0
      DTEMP = 0.0D0
      IF( N.LE.0 .OR. INCX.LE.0 )RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        code for increment not equal to 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DTEMP = DTEMP + ABS(DX(I))
   10 CONTINUE
CS    SASUM = DTEMP
CD    DASUM = DTEMP
      RETURN
C
C        code for increment equal to 1
C
C
C        clean-up loop
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + ABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + ABS(DX(I)) + ABS(DX(I + 1)) + ABS(DX(I + 2))
     *  + ABS(DX(I + 3)) + ABS(DX(I + 4)) + ABS(DX(I + 5))
   50 CONTINUE
   60 CONTINUE
CS    SASUM = DTEMP
CD    DASUM = DTEMP
      RETURN
      END
