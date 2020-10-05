* *******************************************************************
* COPYRIGHT (c) 1976 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Original date 21 Jan 1993
C       Toolpack tool decs employed.
C	Double version of MC13D (name change only)
C 10 August 2001 DOs terminated with CONTINUE
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC13DD(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW)
      INTEGER LICN,N,NUM
      INTEGER IB(N),ICN(LICN),IOR(N),IP(N),IW(N,3),LENR(N)
      EXTERNAL MC13ED
      CALL MC13ED(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW(1,1),IW(1,2),IW(1,3))
      RETURN
      END
      SUBROUTINE MC13ED(N,ICN,LICN,IP,LENR,ARP,IB,NUM,LOWL,NUMB,PREV)
      INTEGER LICN,N,NUM
      INTEGER ARP(N),IB(N),ICN(LICN),IP(N),LENR(N),LOWL(N),NUMB(N),
     +        PREV(N)
      INTEGER DUMMY,I,I1,I2,ICNT,II,ISN,IST,IST1,IV,IW,J,K,LCNT,NNM1,STP
      INTRINSIC MIN
      ICNT = 0
      NUM = 0
      NNM1 = N + N - 1
      DO 20 J = 1,N
        NUMB(J) = 0
        ARP(J) = LENR(J) - 1
   20 CONTINUE
      DO 120 ISN = 1,N
        IF (NUMB(ISN).NE.0) GO TO 120
        IV = ISN
        IST = 1
        LOWL(IV) = 1
        NUMB(IV) = 1
        IB(N) = IV
        DO 110 DUMMY = 1,NNM1
          I1 = ARP(IV)
          IF (I1.LT.0) GO TO 60
          I2 = IP(IV) + LENR(IV) - 1
          I1 = I2 - I1
          DO 50 II = I1,I2
            IW = ICN(II)
            IF (NUMB(IW).EQ.0) GO TO 100
            LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
   50     CONTINUE
          ARP(IV) = -1
   60     IF (LOWL(IV).LT.NUMB(IV)) GO TO 90
          NUM = NUM + 1
          IST1 = N + 1 - IST
          LCNT = ICNT + 1
          DO 70 STP = IST1,N
            IW = IB(STP)
            LOWL(IW) = N + 1
            ICNT = ICNT + 1
            NUMB(IW) = ICNT
            IF (IW.EQ.IV) GO TO 80
   70     CONTINUE
   80     IST = N - STP
          IB(NUM) = LCNT
          IF (IST.NE.0) GO TO 90
          IF (ICNT.LT.N) GO TO 120
          GO TO 130
   90     IW = IV
          IV = PREV(IV)
          LOWL(IV) = MIN(LOWL(IV),LOWL(IW))
          GO TO 110
  100     ARP(IV) = I2 - II - 1
          PREV(IW) = IV
          IV = IW
          IST = IST + 1
          LOWL(IV) = IST
          NUMB(IV) = IST
          K = N + 1 - IST
          IB(K) = IV
  110   CONTINUE
  120 CONTINUE
  130 DO 140 I = 1,N
        II = NUMB(I)
        ARP(II) = I
  140 CONTINUE
      RETURN
      END
