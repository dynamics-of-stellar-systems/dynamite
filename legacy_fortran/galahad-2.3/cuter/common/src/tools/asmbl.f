C     ( Last modified on 23 Dec 2000 at 22:01:38 )
CS    SUBROUTINE SASMBL( N , NG, MAXSEL, NSEMIB, LH    , LIH   , NNZH  ,
CD    SUBROUTINE DASMBL( N , NG, MAXSEL, NSEMIB, LH    , LIH   , NNZH  ,
     *                   NFREE , IFREE , ISTADH, LSTADH, ICNA  , LICNA ,  
     *                   ISTADA, LSTADA, INTVAR, LNTVAR, IELVAR, LELVAR, 
     *                   IELING, LELING, ISTADG, LSTADG, ISTAEV, LSTAEV, 
     *                   ISTAGV, LNSTGV, ISVGRP, LNVGRP, IRNH  , JCNH  ,   
     *                   NXTROW, LNXTRW, IWK   , LIWK  , A , LA, GUVALS, 
     *                   LNGUVL, HUVALS, LNHUVL, GVALS2, GVALS3, GSCALE, 
     *                   ESCALE, LESCAL, H , WK, LWK   , GXEQX , LGXEQX, 
     *                   INTREP, LINTRE, ITYPEE, LITYPE, RANGE , IPRINT, 
     *                   IOUT  , BAND  , MAXSBW, INFORM, NOZERO, FIXSTR)
C
C  ********************************************************************
C
C  Assemble the second derivative matrix of a groups partially
C  separable function in either co-ordinate or band format.
C
C  Nick Gould, February 20TH 1989, for CGT Productions.
C
C  ********************************************************************
C
      INTEGER          N , NG, MAXSEL, NFREE , LNXTRW, NNZH  , LITYPE
      INTEGER          LA, LH, IOUT  , INFORM, NSEMIB, IPRINT, LIH
      INTEGER          LSTADH, LICNA , LSTADA, LNTVAR, LELVAR, LELING
      INTEGER          LSTADG, LSTAEV, LNSTGV, LNVGRP, LIWK,   MAXSBW
      INTEGER          LNGUVL, LNHUVL, LESCAL, LWK,    LGXEQX, LINTRE
      LOGICAL          BAND  , NOZERO, FIXSTR
      INTEGER          IFREE(  N ),      ISTADH( LSTADH ), ICNA( LICNA )  
      INTEGER          ISTADA( LSTADA ), INTVAR( LNTVAR ), JCNH( LIH )
      INTEGER          IELVAR( LELVAR ), IELING( LELING ), IRNH( LIH )
      INTEGER          ISTADG( LSTADG ), ISTAEV( LSTAEV )
      INTEGER          ISTAGV( LNSTGV ), ISVGRP( LNVGRP ), IWK( LIWK )
      INTEGER          ITYPEE( LITYPE )
      INTEGER          NXTROW( 2,        LNXTRW )
      LOGICAL          GXEQX( LGXEQX ),  INTREP( LINTRE )
CS    REAL             A(     LA ),      GVALS2( NG ),     GVALS3( NG ),
CD    DOUBLE PRECISION A(     LA ),      GVALS2( NG ),     GVALS3( NG ),
     *                 GUVALS( LNGUVL ), HUVALS( LNHUVL ), GSCALE( NG ), 
     *                 ESCALE( LESCAL ), H(      LH ),     WK(     LWK )
      EXTERNAL         RANGE 
C
C  LOCAL VARIABLES.
C
      INTEGER          I, II,  IG, J,  JJ, K,  L, IP,  NN,     NNN
      INTEGER          NEWPT,  NIN,    IELL,   IEL,    IHNEXT
      INTEGER          NVAREL, IG1,    LISTVS, LISTVE, IELH
      INTEGER          INEXT,  IJHESS, IROW,   JCOL,   JCOLST, ISTART
CS    REAL             ONE,    ZERO,   WKI,    HESNEW, GDASH,  G2DASH,
CD    DOUBLE PRECISION ONE,    ZERO,   WKI,    HESNEW, GDASH,  G2DASH,
     *                 SCALEE
      CHARACTER * 2    MATRIX( 36, 36 )
CS    EXTERNAL         SSETVL, SSETVI
CD    EXTERNAL         DSETVL, DSETVI
      INTRINSIC        ABS,    MAX,    MIN
C
C  Set constant real parameters.
C
CS    PARAMETER ( ZERO   = 0.0E+0, ONE    = 1.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0, ONE    = 1.0D+0 )
      IF ( .NOT. BAND .AND. NFREE .GT. LNXTRW ) GO TO 620
C
C  Renumber the free variables so that they are variables 1 to NFREE.
C
      DO 10 I     = 1, N
         IWK( I ) = 0
   10 CONTINUE
      DO 20 I = 1, NFREE
         IWK( IFREE( I ) ) = I
C
C  Initialize the link list which points to the row numbers which
c  are used in the columns of the assembled Hessian.
C
C  NXTROW( 1, . ) gives the link list. the list for column J starts
C                 in NXTROW( 1, J ) and ends when NXTROW( 1, K ) = - 1.
C  NXTROW( 2, . ) gives the position in H of the current link.
C
         IF ( .NOT. BAND ) NXTROW( 1, I ) = - 1
   20 CONTINUE
      IF ( IPRINT .GE. 10 ) THEN
         WRITE( IOUT, 2060 ) NFREE, ( IFREE( I ), I = 1, NFREE )
      END IF
C
C  If a band storage scheme is to be used, initialize the entries
C  within the band as zero.
C
      IF ( BAND ) THEN
         MAXSBW  = 0
         IF ( NFREE * ( NSEMIB + 1 ) .GT. LH ) GO TO 610
         DO 30 I = 1, NFREE * ( NSEMIB + 1 )
            H( I ) = ZERO
   30    CONTINUE   
      ELSE
         NNZH  = 0
         NEWPT = NFREE + 1
         IF ( NEWPT .GT. LNXTRW ) GO TO 620
      END IF
C
C  ------------------------------------------------------
C  Form the rank-one second order term for the Ith group.
C  ------------------------------------------------------
C
      DO 200 IG = 1, NG
         IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2070 ) IG
         IF ( GXEQX( IG ) ) GO TO 200
         IF ( .NOT. FIXSTR .AND. GSCALE( IG ) .EQ. ZERO ) GO TO 200
         G2DASH = GSCALE( IG ) * GVALS3( IG )
         IF ( IPRINT .GE. 100 ) WRITE( 6, * ) ' GVALS3(IG) ', GVALS3(IG)
         IF ( NOZERO .AND. G2DASH .EQ. ZERO ) GO TO 200
         IG1    = IG + 1
         LISTVS = ISTAGV( IG )
         LISTVE = ISTAGV( IG1 ) - 1
C
C  Form the gradient of the IG-th group.
C
CS       CALL SSETVI( LISTVE - LISTVS + 1, WK, ISVGRP( LISTVS ), ZERO )
CD       CALL DSETVI( LISTVE - LISTVS + 1, WK, ISVGRP( LISTVS ), ZERO )
C
C  Consider any nonlinear elements for the group.
C
         DO 130 IELL = ISTADG( IG ), ISTADG( IG1 ) - 1
            IEL      = IELING( IELL )
            K        = INTVAR( IEL )
            L        = ISTAEV( IEL )
            NVAREL   = ISTAEV( IEL + 1 ) - L
            SCALEE   = ESCALE( IELL )
            IF ( INTREP( IEL ) ) THEN
C
C  The IEL-th element has an internal representation.
C
               NIN = INTVAR( IEL + 1 ) - K
               CALL RANGE ( IEL, .TRUE., GUVALS( K ),
     *                      WK( N + 1 ), NVAREL, NIN, 
     *                      ITYPEE( IEL ), NIN, NVAREL )
               DO 110 I   = 1, NVAREL
                  J       = IELVAR( L )
                  WK( J ) = WK( J ) + SCALEE * WK( N + I )
                  L       = L + 1
  110          CONTINUE
            ELSE
C
C  The IEL-th element has no internal representation.
C
               DO 120 I   = 1, NVAREL
                  J       = IELVAR( L )
                  WK( J ) = WK( J ) + SCALEE * GUVALS( K )
                  K       = K + 1
                  L       = L + 1
  120          CONTINUE
            END IF
  130    CONTINUE
C
C  Include the contribution from the linear element.
C
         DO 140 K   = ISTADA( IG ), ISTADA( IG1 ) - 1
            J       = ICNA( K )
            WK( J ) = WK( J ) + A( K )
  140    CONTINUE
C
C  The gradient is complete. form the J-th column of the rank-one matrix
C
         DO 190 L = LISTVS, LISTVE
            JJ    = ISVGRP( L )
            J     = IWK( JJ )
            IF ( J .EQ. 0 ) GO TO 190
C
C  Find the entry in row i of this column.
C
            DO 180 K = LISTVS, LISTVE
               II    = ISVGRP( K )
               I     = IWK( II )
               IF ( I .EQ. 0 .OR. I .GT. J ) GO TO 180
C
C  Skip all elements which lie outside a band of width NSEMIB.
C
               IF ( BAND ) MAXSBW = MAX( MAXSBW, J - I )
               IF ( J - I .GT. NSEMIB ) GO TO 180
               HESNEW = WK( II ) * WK( JJ ) * G2DASH
               IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2090 ) I, J, HESNEW
               IF ( NOZERO .AND. HESNEW .EQ. ZERO ) GO TO 180
C
C  Obtain the appropriate storage location in H for the new entry.
C
C
C  Case 1: band matrix storage scheme.
C
               IF ( BAND ) THEN
C
C  The entry belongs on the diagonal.
C
                  IF ( I .EQ. J ) THEN
                     H( I ) = H( I ) + HESNEW
C
C  The entry belongs off the diagonal.
C
                  ELSE
                     H( NFREE + J - I + NSEMIB * ( I - 1 ) ) = 
     *         H( NFREE + J - I + NSEMIB * ( I - 1 ) ) + HESNEW
                  END IF   
C
C  Case 2: coordinate storage scheme.
C
               ELSE
                  ISTART = J
  150             CONTINUE
                  INEXT = NXTROW( 1, ISTART )
                  IF ( INEXT .EQ. - 1 ) THEN
                     IF ( NEWPT .GT. LNXTRW ) GO TO 620
C
C  The (I,J)-th location is empty. Place the new entry in this location
C  and add another link to the list.
C
                     NNZH = NNZH + 1
                     IF ( NNZH .GT. LH .OR. NNZH .GT. LIH ) GO TO 610
                     IRNH( NNZH )        = I
                     JCNH( NNZH )        = J
                     H( NNZH )           = HESNEW
                     NXTROW( 1, ISTART ) = NEWPT
                     NXTROW( 2, ISTART ) = NNZH
                     NXTROW( 1, NEWPT )  = - 1
                     NEWPT               = NEWPT + 1
                  ELSE
C
C  Continue searching the linked list for an entry in row I, column J.
C
                     IF ( IRNH( NXTROW( 2, ISTART ) ) .EQ. I ) THEN
                        IP = NXTROW( 2, ISTART )
                        H( IP ) = H( IP ) + HESNEW
                     ELSE
                        ISTART = INEXT
                        GO TO 150
                     END IF
                  END IF
               END IF
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
C
C  Reset the workspace array to zero.
C
CS    CALL SSETVL( MAXSEL, WK, 1, ZERO )
CD    CALL DSETVL( MAXSEL, WK, 1, ZERO )
C
C  --------------------------------------------------------
C  Add on the low rank first order terms for the I-th group.
C  --------------------------------------------------------
C
      DO 300 IG = 1, NG
         IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2100 ) IG
         IF ( .NOT. FIXSTR .AND. GSCALE( IG ) .EQ. ZERO ) GO TO 300
         IF ( GXEQX( IG ) ) THEN
            GDASH = GSCALE( IG )
         ELSE
            GDASH = GSCALE( IG ) * GVALS2( IG )
            IF ( NOZERO .AND. GDASH .EQ. ZERO ) GO TO 300
         END IF
         IG1 = IG + 1
C
C  See if the group has any nonlinear elements.
C
         DO 290 IELL = ISTADG( IG ), ISTADG( IG1 ) - 1
            IEL      = IELING( IELL )
            LISTVS   = ISTAEV( IEL )
            LISTVE   = ISTAEV( IEL + 1 ) - 1
            NVAREL   = LISTVE - LISTVS + 1
            IELH     = ISTADH( IEL )
            IHNEXT   = IELH
            SCALEE   = ESCALE( IELL )
            DO 250 L = LISTVS, LISTVE
               J     = IWK( IELVAR( L ) )
               IF ( J .NE. 0 ) THEN
C
C  The IEL-th element has an internal representation.
C  Compute the J-th column of the element Hessian matrix.
C
                  IF ( INTREP( IEL ) ) THEN
C
C  Compute the J-th column of the Hessian.
C
                     WK( L - LISTVS + 1 ) = ONE
C
C  Find the internal variables.
C
                     NIN = INTVAR( IEL + 1 ) - INTVAR( IEL )
                     CALL RANGE ( IEL, .FALSE., WK( 1 ),
     *                        WK( MAXSEL + 1 ), NVAREL, NIN,
     *                        ITYPEE( IEL ), NVAREL, NIN )
C
C  Multiply the internal variables by the element Hessian.
C
                     NN = MAXSEL + NIN
CS                   CALL SSETVL( NIN, WK( NN + 1 ), 1, ZERO )
CD                   CALL DSETVL( NIN, WK( NN + 1 ), 1, ZERO )
C
C  Only the upper triangle of the element Hessian is stored.
C
                     JCOLST = IELH - 1
                     DO 220 JCOL = 1, NIN
                        IJHESS   = JCOLST
                        JCOLST   = JCOLST + JCOL
                        WKI      = WK( MAXSEL + JCOL ) * GDASH
                        DO 210 IROW = 1, NIN
                           IF ( IROW .LE. JCOL ) THEN
                              IJHESS = IJHESS + 1
                           ELSE
                              IJHESS = IJHESS + IROW - 1
                           END IF
                           WK( NN + IROW ) = WK( NN + IROW ) +
     *                                       WKI * HUVALS( IJHESS )
  210                   CONTINUE
  220                CONTINUE
C
C  Scatter the product back onto the elemental variables.
C
                     NNN = NN + NIN
                     CALL RANGE ( IEL, .TRUE., WK( NN + 1 ),
     *                            WK( NNN + 1 ), NVAREL, NIN,
     *                            ITYPEE( IEL ), NIN, NVAREL )
                     WK( L - LISTVS + 1 ) = ZERO
                  END IF
C
C  Find the entry in row I of this column.
C
                  DO 240 K = LISTVS, L
                     I     = IWK( IELVAR( K ) )
C
C  Skip all elements which lie outside a band of width NSEMIB.
C
                     IF ( BAND .AND. I .NE. 0 ) MAXSBW = 
     *                    MAX( MAXSBW, ABS( J - I ) )
                     IF ( ABS( I - J ) .LE. NSEMIB .AND. I .NE. 0 ) THEN
C
C  Only the upper triangle of the matrix is stored.
C
                        IF ( I .LE. J ) THEN
                           II = I
                           JJ = J
                        ELSE
                           II = J
                           JJ = I
                        END IF
C
C  Obtain the appropriate storage location in H for the new entry.
C
                        IF ( INTREP( IEL ) ) THEN
                           HESNEW = SCALEE * WK( NNN + K - LISTVS + 1 )
                        ELSE
                           HESNEW = SCALEE * HUVALS( IHNEXT ) * GDASH
                        END IF
                        IF ( IPRINT .GE. 100 )
     *                     WRITE( 6, 2080 ) II, JJ, IEL, HESNEW
C
C  Case 1: band matrix storage scheme.
C
                        IF ( BAND ) THEN
C
C  The entry belongs on the diagonal.
C
                           IF ( II .EQ. JJ ) THEN
                              H( II ) = H( II ) + HESNEW
                              IF ( K .NE. L ) H( II ) = H( II ) + HESNEW
C
C  The entry belongs off the diagonal.
C
                           ELSE
                              H( NFREE + JJ - II + NSEMIB * ( II - 1 ) )
     *                     =  H( NFREE + JJ - II + NSEMIB * ( II - 1 ) )
     *                       + HESNEW
                           END IF   
C
C  Case 2: coordinate storage scheme.
C
                        ELSE
                           IF ( .NOT. NOZERO .OR. HESNEW .NE. ZERO)THEN
                              ISTART = JJ
  230                         CONTINUE
                              INEXT = NXTROW( 1, ISTART )
                              IF ( INEXT .EQ. - 1 ) THEN
                                 IF ( NEWPT .GT. LNXTRW ) GO TO 620
C
C  The (I,J)-th location is empty. Place the new entry in this location
C  and add another link to the list.
C
                                 NNZH = NNZH + 1
                                 IF ( NNZH .GT. LH .OR. 
     *                                NNZH .GT. LIH ) GO TO 610
                                 IRNH( NNZH )        = II
                                 JCNH( NNZH )        = JJ
                                 H( NNZH )           = HESNEW
                                 IF ( K .NE. L .AND. II .EQ. JJ )  
     *                              H( NNZH ) = HESNEW + HESNEW
                                 NXTROW( 1, ISTART ) = NEWPT
                                 NXTROW( 2, ISTART ) = NNZH
                                 NXTROW( 1, NEWPT )  = - 1
                                 NEWPT               = NEWPT + 1
                              ELSE
C
C  Continue searching the linked list for an entry in row I, column J.
C
                              IF ( IRNH( NXTROW( 2, ISTART ) )
     *                             .EQ. II ) THEN
                                    IP      = NXTROW( 2, ISTART )
                                    H( IP ) = H( IP ) + HESNEW
                                    IF ( K .NE. L .AND. II .EQ. JJ )  
     *                                 H( IP ) = H( IP ) + HESNEW
                                 ELSE
                                    ISTART = INEXT
                                    GO TO 230
                                 END IF
                              END IF
                           END IF
                        END IF
                     END IF
                     IHNEXT = IHNEXT + 1
  240             CONTINUE
               ELSE
                  IHNEXT = IHNEXT + L - LISTVS + 1
               END IF
  250       CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C  ----------------------------------------
C  For debugging, print the nonzero values.
C  ----------------------------------------
C
      IF ( IPRINT .GE. 10 ) THEN
         IF ( .NOT. BAND ) WRITE( IOUT, 2000 )
     *      ( IRNH( I ), JCNH( I ), H( I ), I = 1, NNZH )
         IF ( NFREE .LE. 36 ) THEN
            DO 420 I = 1, NFREE
               DO 410 J = 1, NFREE
                  MATRIX( I, J ) = '  '
  410          CONTINUE
  420       CONTINUE
            IF ( BAND ) THEN
               DO 440 I = 1, NFREE
                  IF ( H( I ) .NE. ZERO ) MATRIX( I, I ) = ' *'
                  DO 430 J = 1, MIN( NSEMIB, NFREE - I )
                     K = NFREE + J + NSEMIB * ( I - 1 )
                     IF ( H( K ) .NE. ZERO ) THEN
                        MATRIX( I + J, I ) = ' *'
                        MATRIX( I, I + J ) = ' *'
                     END IF
  430             CONTINUE
  440          CONTINUE
            ELSE
               DO 450 I = 1, NNZH
                  IF ( IRNH( I ) .GT. NFREE ) THEN
                     WRITE( IOUT, * ) ' ENTRY OUT OF BOUNDS IN ASEMBL ',
     *                                ' ROW NUMBER = ', IRNH( I )
                     STOP
                  END IF
                  IF ( JCNH( I ) .GT. NFREE ) THEN
                     WRITE( IOUT, * ) ' ENTRY OUT OF BOUNDS IN ASEMBL ',
     *                                ' COL NUMBER = ', JCNH( I )
                     STOP
                  END IF
                  MATRIX( IRNH( I ), JCNH( I ) ) = ' *'
                  MATRIX( JCNH( I ), IRNH( I ) ) = ' *'
  450          CONTINUE
            END IF
            WRITE( IOUT, 2040 ) ( I, I = 1, NFREE )
            DO 460 I = 1, NFREE
               WRITE( IOUT, 2050 ) I, ( MATRIX( I, J ), J = 1, NFREE )
  460       CONTINUE
         END IF
      END IF
      INFORM = 0
      RETURN
C
C  Unsuccessful returns.
C
  610 CONTINUE
      INFORM = 1
      IF ( LH .LE. LIH ) THEN
         WRITE( IOUT, 2010 ) LH
      ELSE
         WRITE( IOUT, 2030 ) LIH
      END IF
      RETURN
  620 CONTINUE
      INFORM = 2
      WRITE( IOUT, 2020 ) LNXTRW
      RETURN
C
C  Non-executable statements
C
 2000 FORMAT( '    Row  Column    Value        Row  Column    Value ', /
     *        '    ---  ------    -----        ---  ------    ----- ', /
     *        ( 2I6, 1P, D24.16, 2I6, 1P, D24.16 ) )
 2010 FORMAT( ' ** Array dimension LH =', I6, ' too small in ASMBL. ' )
 2020 FORMAT( ' ** Array dimension LNXTRW =',I8, ' too small in ASMBL.')
 2030 FORMAT( ' ** Array dimension LIH =', I6, ' too small in ASMBL.' )
 2040 FORMAT( /, 5X, 36I2 )
 2050 FORMAT( I3, 2X, 36A2 )
 2060 FORMAT( /, I6, ' free variables. They are ', 8I5, /, ( 14I5 ) )
 2070 FORMAT( ' Group ', I5, ' rank-one terms ' )
 2080 FORMAT( ' Row ', I6, ' Column ', I6, ' used from element ', I6,
     *        ' value = ', 1P, D24.16 )
 2090 FORMAT( ' Row ', I6, ' column ', I6, ' used. Value = ',1P,D24.16)
 2100 FORMAT( ' Group ', I5, ' second-order terms ' )
C
C  END OF ASMBL.
C
      END
