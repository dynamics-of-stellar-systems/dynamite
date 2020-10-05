C *******************************************************************
C COPYRIGHT (c) 1994 Council for the Central Laboratory
*                    of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
C SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
C
C Please note that for an ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date 30 November 1995

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MA51AD(M,N,LA,IRN,KEEP,RANK,ROWS,COLS,W)
      INTEGER M,N,LA,IRN(LA),KEEP(*),RANK,ROWS(M),COLS(N),W(*)
      INTEGER I,IQB(1),IPTRD,IPTRL,IPTRO,IPTRU,J,JB,J1,K1,K2,
     +        KBLOCK,MBLOCK,NB,NBLOCK,NC,NE,NP,NR,RANKB
      EXTERNAL MA51ZD
      IQB(1) = 0
      IPTRL = M + N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      NB = KEEP(KBLOCK+3)
      NE = KEEP(IPTRO+N+1)-1
      K2 = 0
      RANK = 0
      DO 20 JB = 1,NB
         NC = KEEP(NBLOCK+3*JB)
         NR = NC
         IF(NB.EQ.1) NR = M
         K1 = K2 + 1
         K2 = K2 + NC
         NP = KEEP(MBLOCK+3*JB)
         IF (NP.LT.0) THEN
            DO 10 J = K1,K2
               ROWS(J) = 1
               COLS(J) = 1
   10       CONTINUE
            RANK = RANK + NC
         ELSE
            J1 = 1
            IF (JB.GT.1) J1 = KEEP(KBLOCK+3*JB)
            CALL MA51ZD(NR,NC,IQB,NP,LA+1-NE-J1,IRN(NE+J1),
     +                 KEEP(IPTRL+K1),KEEP(IPTRU+K1),RANKB,
     +                 ROWS(K1),COLS(K1),W)
            RANK = RANK + RANKB
         END IF
   20 CONTINUE
      DO 30 I = 1,M
         W(I) = ROWS(KEEP(I))
   30 CONTINUE
      K1 = 1
      K2 = M
      DO 40 I = 1,M
         IF(W(I).EQ.1)THEN
            ROWS(K1) = I
            K1 = K1 + 1
         ELSE
            ROWS(K2) = I
            K2 = K2 - 1
         END IF
   40 CONTINUE
      DO 50 I = 1,N
         W(KEEP(M+I)) = COLS(I)
   50 CONTINUE
      K1 = 1
      K2 = N
      DO 60 I = 1,N
         IF(W(I).EQ.1)THEN
            COLS(K1) = I
            K1 = K1 + 1
         ELSE
            COLS(K2) = I
            K2 = K2 - 1
         END IF
   60 CONTINUE
      END
      SUBROUTINE MA51BD(M,N,IQ,NP,LFACT,IRNF,IPTRL,
     +                  IPTRU,RANK,ROWS,COLS,W)
      INTEGER M,N,IQ(*),NP,LFACT,IRNF(LFACT),IPTRL(N),IPTRU(N)
      INTEGER RANK,ROWS(M),COLS(N),W(*)
      INTEGER I,K1,K2
      EXTERNAL MA51ZD
      CALL MA51ZD(M,N,IQ,NP,LFACT,IRNF,IPTRL,
     +                  IPTRU,RANK,ROWS,COLS,W)
      DO 10 I = 1,M
         W(I) = ROWS(I)
   10 CONTINUE
      K1 = 1
      K2 = M
      DO 20 I = 1,M
         IF(W(I).EQ.1)THEN
            ROWS(K1) = I
            K1 = K1 + 1
         ELSE
            ROWS(K2) = I
            K2 = K2 - 1
         END IF
   20 CONTINUE
      DO 30 I = 1,N
         W(I) = COLS(I)
   30 CONTINUE
      K1 = 1
      K2 = N
      DO 40 I = 1,N
         IF(W(I).EQ.1)THEN
            COLS(K1) = I
            K1 = K1 + 1
         ELSE
            COLS(K2) = I
            K2 = K2 - 1
         END IF
   40 CONTINUE
      END
      SUBROUTINE MA51CD(M,N,LA,A,IRN,KEEP,SGNDET,LOGDET,W)
      INTEGER M,N,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),KEEP(*),SGNDET
      DOUBLE PRECISION LOGDET
      INTEGER W(N)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      INTEGER I,IQB(1),IPTRD,IPTRL,IPTRO,IPTRU,J,JB,J1,K,K1,K2,
     +        KBLOCK,L,MBLOCK,NB,NBLOCK,NC,NE,NP,NR,SGNDT
      DOUBLE PRECISION PIV,LOGDT
      EXTERNAL MA51DD
      IQB(1) = 0
      IPTRL = M + N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      NB = KEEP(KBLOCK+3)
      NE = KEEP(IPTRO+N+1)-1
      LOGDET = ZERO
      SGNDET = 1
      IF (M.NE.N) THEN
          SGNDET = 0
          RETURN
      END IF
      DO 10 I = 1,N
         W(I) = 1
   10 CONTINUE
      DO 40 I = 1,N
         K = KEEP(I)
         IF (W(K).EQ.1) THEN
            DO 30 J = 1,N
               L = KEEP(K)
               W(K) = 0
               IF(K.EQ.I)GO TO 40
               SGNDET = -SGNDET
               K = L
   30       CONTINUE
         END IF
   40 CONTINUE
      DO 60 I = 1,N
         K = KEEP(N+I)
         IF (W(K).EQ.0) THEN
            DO 50 J = 1,N
               L = KEEP(N+K)
               W(K) = 1
               IF(K.EQ.I)GO TO 60
               SGNDET = -SGNDET
               K = L
   50       CONTINUE
         END IF
   60 CONTINUE
      K2 = 0
      DO 80 JB = 1,NB
         NC = KEEP(NBLOCK+3*JB)
         NR = NC
         IF(NB.EQ.1) NR = M
         K1 = K2 + 1
         K2 = K2 + NC
         NP = KEEP(MBLOCK+3*JB)
         IF (NP.LT.0) THEN
            DO 70 J = K1,K2
               PIV = A(KEEP(IPTRD+J))
               IF(PIV.LT.ZERO)THEN
                  SGNDET  = -SGNDET
                  LOGDET  = LOGDET + LOG(-PIV)
               ELSE
                  LOGDET  = LOGDET + LOG(PIV)
               END IF
   70       CONTINUE
         ELSE
            J1 = 1
            IF (JB.GT.1) J1 = KEEP(KBLOCK+3*JB)
            CALL MA51DD(NR,NC,IQB,NP,LA+1-J1-NE,A(NE+J1),IRN(NE+J1),
     +                  KEEP(IPTRL+K1),KEEP(IPTRU+K1),SGNDT,LOGDT,W)
            SGNDET = SGNDET*SGNDT
            LOGDET = LOGDET+LOGDT
            IF (SGNDET.EQ.0) THEN
               LOGDET = 0
               RETURN
            END IF
         END IF
   80 CONTINUE
      END
      SUBROUTINE MA51DD(M,N,IQ,NP,LFACT,FACT,IRNF,IPTRL,
     +                  IPTRU,SGNDET,LOGDET,W)
      INTEGER K,L,M,N,IQ(*),NP,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N),SGNDET
      DOUBLE PRECISION LOGDET
      INTEGER W(N)
      DOUBLE PRECISION A
      INTEGER I,IA1,IF1,J,MF,NF
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      IF1 = IPTRL(N) + 1
      MF = IRNF(2)
      NF = N - NP
      LOGDET = 0
      SGNDET = 1
      IF (MF.GT.0) THEN
         IF (M.NE.N .OR. MF.NE.NF .OR. IRNF(IF1+MF+NF-1).LT.0) THEN
            SGNDET = 0
            RETURN
         END IF
         CALL MA51XD(MF,FACT(IF1),MF,IRNF(IF1+MF),SGNDET,LOGDET)
      END IF
      DO 10 I = 1,NP
         IA1 = IPTRU(I) + 1
         A = FACT(IA1)
         IF(A.LT.ZERO)THEN
            SGNDET  = -SGNDET
            LOGDET  = LOGDET - LOG(-A)
         ELSE
            LOGDET  = LOGDET - LOG(A)
         END IF
         W(I) = IRNF(IA1)
   10 CONTINUE
      DO 20 I = 1,MF
         W(NP+I) = IRNF(IF1+I-1)
   20 CONTINUE
      DO 40 I = 1,N
         K = W(I)
         IF (K.GT.0) THEN
            DO 30 J = 1,N
               L = W(K)
               W(K) = 0
               IF(K.EQ.I)GO TO 40
               SGNDET = -SGNDET
               K = L
   30       CONTINUE
         END IF
   40 CONTINUE
      IF(IQ(1).LE.0)RETURN
      DO 60 I = 1,N
         K = IQ(I)
         IF (W(K).EQ.0) THEN
            DO 50 J = 1,N
               L = IQ(K)
               W(K) = 1
               IF(K.EQ.I)GO TO 60
               SGNDET = -SGNDET
               K = L
   50       CONTINUE
         END IF
   60 CONTINUE
      END
      SUBROUTINE MA51XD(N,A,LDA,IPIV,SGNDET,LOGDET)
      IMPLICIT NONE
      INTEGER LDA,N
      DOUBLE PRECISION A(LDA,N)
      INTEGER IPIV(N),SGNDET
      DOUBLE PRECISION LOGDET
      INTEGER J
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      INTRINSIC LOG
      SGNDET  = 1
      LOGDET  = ZERO
      DO 30 J = 1,N
         IF (IPIV(J).NE.J) SGNDET  = -SGNDET
         IF (A(J,J).LT.ZERO) THEN
            SGNDET  = -SGNDET
            LOGDET  = LOGDET + LOG(-A(J,J))
         ELSE
            LOGDET  = LOGDET + LOG(A(J,J))
         END IF
   30 CONTINUE
      END
      SUBROUTINE MA51YD(M,N,IPIV,RANK,ROWS,COLS)
      INTEGER M, N, IPIV(N), RANK, ROWS(M), COLS(N)
      INTEGER I,J,K
      DO 10 K = 1,N
         IF (IPIV(K).LT.0) GO TO 20
         ROWS(K) = 1
         COLS(K) = 1
   10 CONTINUE
   20 RANK = K - 1
      DO 30 K = RANK+1, M
          ROWS(K) = 0
   30 CONTINUE
      DO 40 K = RANK+1, N
         COLS(K) = 0
   40 CONTINUE
      DO 50 I = RANK + 1,N
         K = -IPIV(I)
         J = COLS(I)
         COLS(I) = COLS(K)
         COLS(K) = J
   50 CONTINUE
      DO 60 I = RANK,1,-1
         K = IPIV(I)
         J = ROWS(I)
         ROWS(I) = ROWS(K)
         ROWS(K) = J
   60 CONTINUE
      END
      SUBROUTINE MA51ZD(M,N,IQ,NP,LFACT,IRNF,IPTRL,
     +                  IPTRU,RANK,ROWS,COLS,W)
      INTEGER M,N,IQ(*),NP,LFACT,IRNF(LFACT),IPTRL(N),IPTRU(N)
      INTEGER RANK,ROWS(M),COLS(N),W(*)
      INTEGER I,IA1,IF1,J,MF,NF
      EXTERNAL MA51YD
      IF1 = IPTRL(N) + 1
      MF = IRNF(2)
      NF = N - NP
      IF (MF.GT.0 .AND. NF.GT.0) THEN
         CALL MA51YD(MF,NF,IRNF(IF1+MF),RANK,W,COLS(NP+1))
      ELSE
         RANK = 0
         DO 10 I = 1,MF
            W(I) = 0
   10    CONTINUE
         DO 20 I = 1,NF
            COLS(NP+I) = 0
   20    CONTINUE
      END IF
      DO 30 I = 1,M
         ROWS(I) = 1
   30 CONTINUE
      DO 40 I = MF,1,-1
         J = IRNF(IF1+I-1)
         ROWS(J) = W(I)
   40 CONTINUE
      RANK = RANK + M - MF
      DO 220 J = NP,1,-1
         IA1 = IPTRU(J)
         IF (IA1.GE.IPTRL(J)) THEN
            COLS(J) = 0
         ELSE
            COLS(J) = 1
         END IF
  220 CONTINUE
      IF (IQ(1).GT.0) THEN
         DO 230 I = 1,N
            W(I) = COLS(I)
  230    CONTINUE
         DO 240 I = 1,N
            COLS(IQ(I)) = W(I)
  240    CONTINUE
      END IF
      END
