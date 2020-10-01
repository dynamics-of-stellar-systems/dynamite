* *******************************************************************
* COPYRIGHT (c) 1972 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
* Licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* Please note that for an HSL ARCHIVE Licence:
*
* 1. The Package must not be copied for use by any other person.
*    Supply of any part of the library by the Licensee to a third party
*    shall be subject to prior written agreement between AEA
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 8 Dec 1992
C       Toolpack tool decs employed.
C       B21 reference removed in MA20BD.
C       ZERO made PARAMETER.
C       Arg dimensions set to *.
C
      SUBROUTINE MA20AD(Q,D,A,R,S,IQ,M,N,TOLER)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION TOLER
      INTEGER IQ,M,N
      DOUBLE PRECISION A(*),D(*),Q(IQ,*),R(*)
      INTEGER S(*)
      DOUBLE PRECISION B,B21,BIG,MAX,MIN,PIVOT,QOUTIN,SUM
      INTEGER I,IN,J,K,KL,KOUNT,KR,L,M1,M2,N1,N2,OUT
      LOGICAL STAGE,TEST
      DOUBLE PRECISION FD05AD
      EXTERNAL FD05AD
      INTRINSIC DABS
      BIG = FD05AD(5)*.01D0
      M2 = M + 2
      N2 = N + 2
      M1 = M + 1
      N1 = N + 1
      DO 1 J = 1,N
        Q(M2,J) = J
    1 A(J) = 0.
      DO 3 I = 1,M
        Q(I,N2) = N + I
        D(I) = 0.
        IF (Q(I,N1).GE.0) GO TO 3
        DO 2 J = 1,N2
    2   Q(I,J) = -Q(I,J)
    3 CONTINUE
      DO 5 J = 1,N1
        SUM = 0.
        DO 4 I = 1,M
    4   SUM = SUM + Q(I,J)
    5 Q(M1,J) = SUM
      STAGE = .TRUE.
      KOUNT = 0
      KR = 1
      KL = 1
    6 MAX = -1.
      DO 7 J = KR,N
        IF (DABS(Q(M2,J)).GT.N) GO TO 7
        B = DABS(Q(M1,J))
        IF (B.LE.MAX) GO TO 7
        MAX = B
        IN = J
    7 CONTINUE
      IF (Q(M1,IN).GE.0) GO TO 9
      DO 8 I = 1,M2
    8 Q(I,IN) = -Q(I,IN)
    9 K = 0
      DO 10 I = KL,M
        B = Q(I,IN)
        IF (B.LE.TOLER) GO TO 10
        K = K + 1
        R(K) = Q(I,N1)/B
        S(K) = I
        TEST = .TRUE.
   10 CONTINUE
   11 IF (K.GT.0) GO TO 12
      TEST = .FALSE.
      GO TO 14
   12 MIN = BIG
      DO 13 I = 1,K
        IF (R(I).GE.MIN) GO TO 13
        J = I
        MIN = R(I)
        OUT = S(I)
   13 CONTINUE
      R(J) = R(K)
      S(J) = S(K)
      K = K - 1
   14 IF (TEST .OR. .NOT.STAGE) GO TO 16
      DO 15 I = 1,M2
        B = Q(I,KR)
        Q(I,KR) = Q(I,IN)
   15 Q(I,IN) = B
      KR = KR + 1
      GO TO 25
   16 IF (TEST) GO TO 17
      Q(M2,N1) = 2.
      GO TO 34
   17 PIVOT = Q(OUT,IN)
      IF (Q(M1,IN)-PIVOT-PIVOT.LE.TOLER) GO TO 19
      DO 18 J = KR,N1
        B = Q(OUT,J)
        Q(M1,J) = Q(M1,J) - B - B
   18 Q(OUT,J) = -B
      Q(OUT,N2) = -Q(OUT,N2)
      GO TO 11
   19 DO 20 J = KR,N1
        IF (J.EQ.IN) GO TO 20
        Q(OUT,J) = Q(OUT,J)/PIVOT
   20 CONTINUE
      QOUTIN = Q(OUT,IN)
      Q(OUT,IN) = ZERO
      DO 21 J = KR,N1
        IF (J.EQ.IN) GO TO 21
        B21 = -Q(OUT,J)
        DO 22 I = 1,M1
          Q(I,J) = Q(I,J) + B21*Q(I,IN)
   22   CONTINUE
   21 CONTINUE
      Q(OUT,IN) = QOUTIN
      DO 23 I = 1,M1
        IF (I.EQ.OUT) GO TO 23
        Q(I,IN) = -Q(I,IN)/PIVOT
   23 CONTINUE
      Q(OUT,IN) = 1./PIVOT
      B = Q(OUT,N2)
      Q(OUT,N2) = Q(M2,IN)
      Q(M2,IN) = B
      KOUNT = KOUNT + 1
      IF (.NOT.STAGE) GO TO 26
      KL = KL + 1
      DO 24 J = KR,N2
        B = Q(OUT,J)
        Q(OUT,J) = Q(KOUNT,J)
   24 Q(KOUNT,J) = B
   25 IF (KOUNT+KR.NE.N1) GO TO 6
      STAGE = .FALSE.
   26 MAX = -BIG
      DO 28 J = KR,N
        B = Q(M1,J)
        IF (B.GE.0) GO TO 27
        IF (B.GT.-2.) GO TO 28
        B = -B - 2.
   27   IF (B.LE.MAX) GO TO 28
        MAX = B
        IN = J
   28 CONTINUE
      IF (MAX.LE.TOLER) GO TO 30
      IF (Q(M1,IN).GT.0) GO TO 9
      DO 29 I = 1,M2
   29 Q(I,IN) = -Q(I,IN)
      Q(M1,IN) = Q(M1,IN) - 2.
      GO TO 9
   30 L = KL - 1
      DO 32 I = 1,L
        IF (Q(I,N1).GE.0) GO TO 32
        DO 31 J = KR,N2
   31   Q(I,J) = -Q(I,J)
   32 CONTINUE
      Q(M2,N1) = 0.
      IF (KR.NE.1) GO TO 34
      DO 33 J = 1,N
        B = DABS(Q(M1,J))
        IF (B.LE.TOLER .OR. 2D0-B.LE.TOLER) GO TO 34
   33 CONTINUE
      Q(M2,N1) = 1.
   34 DO 37 I = 1,M
        K = Q(I,N2)
        B = Q(I,N1)
        IF (K.GT.0) GO TO 35
        K = -K
        B = -B
   35   IF (I.GE.KL) GO TO 36
        A(K) = B
        GO TO 37
   36   K = K - N
        D(K) = B
   37 CONTINUE
      Q(M2,N2) = KOUNT
      Q(M1,N2) = N1 - KR
      SUM = 0.
      DO 38 I = KL,M
   38 SUM = SUM + Q(I,N1)
      Q(M1,N1) = SUM
      RETURN
      END
      SUBROUTINE MA20BD(Q,D,A,R,S,IQ,M,N,TOLER)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION TOLER
      INTEGER IQ,M,N
      DOUBLE PRECISION A(*),D(*),Q(IQ,*),R(*)
      INTEGER S(M)
      DOUBLE PRECISION A1,A2,B,BIG,MAX,MIN,PIVOT,QOUTIN,SUM
      INTEGER I,IFORCE,IN,IN2,J,K,KOUNT,M1,M2,N1,N2,OUT
      LOGICAL TEST
      DOUBLE PRECISION FD05AD
      EXTERNAL FD05AD
      INTRINSIC DABS,DMAX1,IABS,INT
      BIG = FD05AD(5)*.01D0
      M2 = M + 2
      N2 = N + 2
      M1 = M + 1
      N1 = N + 1
      IFORCE = INT(1.5*N)
      DO 1 J = 1,N
        Q(M2,J) = J
    1 A(J) = 0.
      DO 3 I = 1,M
        Q(I,N2) = N + I
        D(I) = 0.
        IF (Q(I,N1).GE.0) GO TO 3
        DO 2 J = 1,N2
    2   Q(I,J) = -Q(I,J)
    3 CONTINUE
      DO 5 J = 1,N1
        SUM = 0.
        DO 4 I = 1,M
    4   SUM = SUM + Q(I,J)
    5 Q(M1,J) = SUM
      KOUNT = 0
    6 A1 = -BIG
      A2 = -BIG
      DO 9 J = 1,N
        B = Q(M1,J)
        IF (DABS(Q(M2,J)).LE.N) GO TO 8
        IF (B.GE.0.) GO TO 7
        B = -B - 2.
    7   IF (B.LE.A2) GO TO 9
        A2 = B
        IN2 = J
        GO TO 9
    8   IF (B.LE.A1) GO TO 9
        A1 = B
        IN = J
    9 CONTINUE
      IF (KOUNT.LE.IFORCE .AND. A1.GT.TOLER) GO TO 11
      MAX = DMAX1(A1,A2)
      IF (MAX.LE.TOLER) GO TO 24
      IF (MAX.EQ.A1) GO TO 11
      IN = IN2
      IF (Q(M1,IN).GE.0.) GO TO 11
      DO 10 I = 1,M2
   10 Q(I,IN) = -Q(I,IN)
      Q(M1,IN) = Q(M1,IN) - 2.
   11 K = 0
      DO 12 I = 1,M
        B = Q(I,IN)
        IF (B.LE.TOLER) GO TO 12
        K = K + 1
        R(K) = Q(I,N1)/B
        S(K) = I
        TEST = .TRUE.
   12 CONTINUE
   13 IF (K.GT.0) GO TO 14
      TEST = .FALSE.
      GO TO 16
   14 MIN = BIG
      DO 15 I = 1,K
        IF (R(I).GE.MIN) GO TO 15
        J = I
        MIN = R(I)
        OUT = S(I)
   15 CONTINUE
      PIVOT = Q(OUT,IN)
      IF (DABS(Q(OUT,N2)).LE.N) GO TO 19
      R(J) = R(K)
      S(J) = S(K)
      K = K - 1
   16 IF (TEST) GO TO 17
      Q(M2,N1) = 2.
      GO TO 25
   17 IF (Q(M1,IN)-PIVOT-PIVOT.LE.TOLER) GO TO 19
      DO 18 J = 1,N1
        B = Q(OUT,J)
        Q(M1,J) = Q(M1,J) - B - B
   18 Q(OUT,J) = -B
      Q(OUT,N2) = -Q(OUT,N2)
      GO TO 13
   19 DO 20 J = 1,N1
        IF (J.EQ.IN) GO TO 20
        Q(OUT,J) = Q(OUT,J)/PIVOT
   20 CONTINUE
      QOUTIN = Q(OUT,IN)
      Q(OUT,IN) = ZERO
      DO 21 J = 1,N1
        IF (J.EQ.IN) GO TO 21
        B = -Q(OUT,J)
        DO 22 I = 1,M1
          Q(I,J) = Q(I,J) + B*Q(I,IN)
   22   CONTINUE
   21 CONTINUE
      Q(OUT,IN) = QOUTIN
      DO 23 I = 1,M1
        IF (I.EQ.OUT) GO TO 23
        Q(I,IN) = -Q(I,IN)/PIVOT
   23 CONTINUE
      Q(OUT,IN) = 1./PIVOT
      B = Q(OUT,N2)
      Q(OUT,N2) = Q(M2,IN)
      Q(M2,IN) = B
      KOUNT = KOUNT + 1
      GO TO 6
   24 Q(M2,N1) = 1.
   25 DO 28 I = 1,M
        K = Q(I,N2)
        B = Q(I,N1)
        IF (IABS(K).GT.N) GO TO 26
        A(K) = B
        GO TO 28
   26   IF (K.GT.0) GO TO 27
        K = -K
        B = -B
   27   K = K - N
        D(K) = B
   28 CONTINUE
      Q(M2,N2) = KOUNT
      SUM = 0.
      DO 29 I = 1,M
   29 SUM = SUM + DABS(D(I))
      Q(M1,N1) = SUM
      RETURN
      END
