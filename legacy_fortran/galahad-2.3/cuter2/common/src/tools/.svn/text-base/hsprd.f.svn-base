C     ( Last modified on 23 Dec 2000 at 22:01:38 )
C
C  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
C     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.
C
CS    SUBROUTINE SHSPRD( N , NN, NG    , NGEL  , NFREE , NVAR1 , NVAR2 , 
CD    SUBROUTINE DHSPRD( N , NN, NG    , NGEL  , NFREE , NVAR1 , NVAR2 , 
     *                   NBPROD, ALLLIN, IVAR  , ISTAEV, LSTAEV, ISTADH, 
     *                   LSTADH, INTVAR, LNTVAR, IELING, LELING, IELVAR, 
     *                   LELVAR, ISTAJC, LNSTJC, ISELTS, LNELTS, ISPTRS, 
     *                   LNPTRS, IGCOLJ, LNGCLJ, ISLGRP, LNLGRP, ISWKSP, 
     *                   LNWKSP, ISVGRP, LNVGRP, ISTAGV, LNSTGV, IVALJR, 
     *                   LNVLJR, ITYPEE, LITYPE, NNONNZ, INONNZ, LNNNON,
     *                   IUSED , LNIUSE, INONZ2, LNNNO2, ISYMMH, MAXSZH,
     *                   P , Q , GVALS2, GVALS3, GRJAC , LGRJAC,
     *                   GSCALE, ESCALE, LESCAL, HUVALS, LHUVAL,
     *                   WK    , LNWK  , WKB   , LNWKB , WKC   , LNWKC , 
     *                   GXEQX , LGXEQX, INTREP, LINTRE, DENSEP, RANGE )
C
C  EVALUATE Q, THE PRODUCT OF THE HESSIAN OF A GROUP  PARTIALLY
C  SEPARABLE FUNCTION WITH THE VECTOR P.
C
C  THE NONZERO COMPONENTS OF P HAVE INDICES IVAR( I ), I = NVAR1, NVAR2.
C  THE NONZERO COMPONENTS OF THE PRODUCT HAVE INDICES INNONZ( I ),
C  I = 1, NNONNZ.
C
C  THE ELEMENTS OF THE ARRAY IUSED MUST BE SET TO ZERO ON ENTRY.
C  THE WORKSPACE ARRAY WK MUST HAVE LENGTH AT LEAST MAX( NG, N + 2 *
C  MAXIMUM NUMBER OF INTERNAL VARIABLES )
C
C  NICK GOULD, 10TH MAY 1989.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          N, NN, NG, NGEL, NFREE, NVAR1, NVAR2, NBPROD
      INTEGER          NNONNZ, MAXSZH, LGXEQX, LINTRE, LITYPE
      INTEGER          LSTAEV, LSTADH, LNTVAR, LELING, LELVAR
      INTEGER          LNSTJC, LNELTS, LNPTRS, LNGCLJ, LNLGRP, LNWKSP
      INTEGER          LNVGRP, LNSTGV, LNVLJR, LNNNON, LNIUSE, LNNNO2
      INTEGER          LGRJAC, LESCAL, LHUVAL, LNWK, LNWKB, LNWKC
      LOGICAL          ALLLIN, DENSEP
      INTEGER          IVAR( N ), ISTAEV( LSTAEV ), ISTADH( LSTADH )
      INTEGER          ISTAJC( LNSTJC ), ISELTS( LNELTS )
      INTEGER          IGCOLJ( LNGCLJ ), ISLGRP( LNLGRP )
      INTEGER          ISWKSP( LNWKSP ), ISPTRS( LNPTRS )
      INTEGER          INTVAR( LNTVAR ), IELING( LELING )
      INTEGER          IELVAR( LELVAR ), ISVGRP( LNVGRP )
      INTEGER          ISTAGV( LNSTGV ), IVALJR( LNVLJR )
      INTEGER          INONNZ( LNNNON ), IUSED( LNIUSE )
      INTEGER          ITYPEE( LITYPE )
      INTEGER          INONZ2( LNNNO2 ), ISYMMH( MAXSZH, MAXSZH )
CS    REAL             P( N ), GVALS2( NG ), GVALS3( NG ),
CD    DOUBLE PRECISION P( N ), GVALS2( NG ), GVALS3( NG ),
     *                 GSCALE( NG ), ESCALE( LESCAL ),
     *                 WK( LNWK ), WKB( LNWKB ), WKC( LNWKC ),
     *                 Q( N ), GRJAC( LGRJAC ), HUVALS( LHUVAL )
      LOGICAL          GXEQX( LGXEQX ), INTREP( LINTRE )
      EXTERNAL         RANGE 
C
C  LOCAL VARIABLES.
C
      INTEGER          I, IEL, IG, IPT, J, IROW, JCOL, IJHESS, LTHVAR
      INTEGER          IELL, II, K, L, LL, NIN, NVAREL, IELHST, NNONZ2
CS    REAL             ZERO, PI, GI
CD    DOUBLE PRECISION ZERO, PI, GI
      LOGICAL          NULLWK
CS    EXTERNAL         SSETVL
CD    EXTERNAL         DSETVL
      INTRINSIC        ABS
C
C  COMMON VARIABLES.
C
CS    REAL             EPSMCH, EPSNEG, TINY, BIG
CD    DOUBLE PRECISION EPSMCH, EPSNEG, TINY, BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY, BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY, BIG
CS    SAVE   / SMACHN /
CD    SAVE   / DMACHN /
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0 )
C
C  ======================= RANK-ONE TERMS ==========================
C
C  IF THE IG-TH GROUP IS NON-TRIVIAL, FORM THE PRODUCT OF P WITH THE
C  SUM OF RANK-ONE FIRST ORDER TERMS, A(TRANS) * GVALS3 * A. A IS
C  STORED BY BOTH ROWS AND COLUMNS. FOR MAXIMUM EFFICIENCY,
C  THE PRODUCT IS FORMED IN DIFFERENT WAYS IF P IS SPARSE OR DENSE.
C
C  -----------------  CASE 1. P IS NOT SPARSE -----------------------
C
      IF ( DENSEP ) THEN
C
C  INITIALIZE WK AND Q AS ZERO.
C
CS       CALL SSETVL( NG, WK, 1, ZERO )
CD       CALL DSETVL( NG, WK, 1, ZERO )
CS       CALL SSETVL( N,  Q,  1, ZERO )
CD       CALL DSETVL( N,  Q,  1, ZERO )
C
C  FORM THE MATRIX-VECTOR PRODUCT WK = A * P, USING THE COLUMN-WISE
C  STORAGE OF A.
C
         DO 20 J       = NVAR1, NVAR2
            I          = IVAR( J )
            PI         = P( I )
CDIR$ IVDEP
            DO 10 K    = ISTAJC( I ), ISTAJC( I + 1 ) - 1
               L       = IGCOLJ( K )
               WK( L ) = WK( L ) + PI * GRJAC( K )
   10       CONTINUE
   20    CONTINUE
C
C  MULTIPLY WK BY THE DIAGONAL MATRIX GVALS3.
C
         DO 30 IG = 1, NG
            IF ( GXEQX( IG ) ) THEN
               WK( IG ) = GSCALE( IG ) * WK( IG )
            ELSE
               WK( IG ) = GSCALE( IG ) * WK( IG ) * GVALS3( IG )
            END IF
   30    CONTINUE
C
C  FORM THE MATRIX-VECTOR PRODUCT Q = A(TRANS) * WK, ONCE
C  AGAIN USING THE COLUMN-WISE STORAGE OF A.
C
         NNONNZ   = 0
         DO 50 J  = 1, NFREE
            I     = IVAR( J )
            PI    = ZERO
CDIR$ IVDEP
            DO  40 K = ISTAJC( I ), ISTAJC( I + 1 ) - 1
               PI    = PI + WK( IGCOLJ( K ) ) * GRJAC( K )
   40       CONTINUE
            Q( I ) = PI
   50    CONTINUE
      ELSE
C
C  ------------------- CASE 2. P IS SPARSE --------------------------
C
         NNONZ2  = 0
C
C  FORM THE MATRIX-VECTOR PRODUCT WK = A * P, USING THE COLUMN-WISE
C  STORAGE OF A. KEEP TRACK OF THE NONZERO COMPONENTS OF WK IN INONZ2.
C  ONLY STORE COMPONENTS CORRESPONDING TO NON TRIVIAL GROUP .
C
         DO 120 J    = NVAR1, NVAR2
            I        = IVAR( J )
            PI       = P( I )
CDIR$ IVDEP
            DO 110 K = ISTAJC( I ), ISTAJC( I + 1 ) - 1
               IG    = IGCOLJ( K )
               IF ( IUSED( IG ) .EQ. 0 ) THEN
                  WK( IG )         = PI * GRJAC( K )
                  IUSED( IG )      = 1
                  NNONZ2           = NNONZ2 + 1
                  INONZ2( NNONZ2 ) = IG
               ELSE
                  WK( IG ) = WK( IG ) + PI * GRJAC( K )
               END IF
  110       CONTINUE
  120    CONTINUE
C
C  RESET IUSED TO ZERO.
C
         DO 130 I = 1, NNONZ2
            IUSED( INONZ2( I ) ) = 0
  130    CONTINUE
C
C  FORM THE MATRIX-VECTOR PRODUCT Q = A(TRANS) * WK, USING THE
C  ROW-WISE STORAGE OF A.
C
         NNONNZ   = 0
         DO 160 J = 1, NNONZ2
            IG    = INONZ2( J )
            IF ( .NOT. GXEQX( IG ) ) THEN
C
C  IF GROUP IG IS NON TRIVIAL, THERE ARE CONTRIBUTIONS FROM ITS
C  RANK-ONE TERM.
C
               PI = GSCALE( IG ) * GVALS3( IG ) * WK( IG )
CDIR$ IVDEP
               DO 150 K = ISTAGV( IG ), ISTAGV( IG + 1 ) - 1
                  L     = ISVGRP( K )
C
C  IF Q HAS A NONZERO IN POSITION L, STORE ITS INDEX IN INONNZ.
C
                  IF ( IUSED( L ) .EQ. 0 ) THEN
                     Q( L )           = PI * GRJAC( IVALJR( K ) )
                     IUSED( L )       = 1
                     NNONNZ           = NNONNZ + 1
                     INONNZ( NNONNZ ) = L
                  ELSE
                     Q( L ) = Q( L ) + PI * GRJAC( IVALJR( K ) )
                  END IF
  150          CONTINUE
            END IF
  160    CONTINUE
      END IF
      IF ( .NOT. ALLLIN ) THEN
C
C  ======================= SECOND-ORDER TERMS =======================
C
C  NOW CONSIDER THE PRODUCT OF P WITH THE SECOND ORDER TERMS (THAT IS.
C  (2ND DERIVATIVES OF THE ELEMENTS). AGAIN, FOR MAXIMUM EFFICIENCY,
C  THE PRODUCT IS FORMED IN DIFFERENT WAYS IF P IS SPARSE OR DENSE.
C
C  --------------------- CASE 1. P IS NOT SPARSE ---------------------
C
         IF ( DENSEP ) THEN
            DO 280 IELL = 1, NGEL
               IEL      = IELING( IELL )
               IG       = ISLGRP( IELL )
               NVAREL   = ISTAEV( IEL + 1 ) - ISTAEV( IEL )
               IF ( GXEQX( IG ) ) THEN
                  GI = GSCALE( IG ) * ESCALE( IELL )
               ELSE
                  GI = GSCALE( IG ) * ESCALE( IELL ) * GVALS2( IG ) 
               END IF
               ISWKSP( IEL ) = NBPROD
               IF ( INTREP( IEL ) ) THEN
C
C  THE IEL-TH ELEMENT HESSIAN HAS AN INTERNAL REPRESENTATION.
C  COPY THE ELEMENTAL VARIABLES INTO WK.
C
                  NULLWK      = .TRUE.
                  LL          = ISTAEV( IEL )
CDIR$ IVDEP
                  DO 210 II   = 1, NVAREL
                     PI       = P( IELVAR( LL ) )
                     WK( II ) = PI
                     IF ( PI .NE. ZERO ) NULLWK = .FALSE.
                     LL       = LL + 1
  210             CONTINUE
C
C  SKIP THE ELEMENT IF WK IS NULL.
C
                  IF ( NULLWK ) GO TO 280
C
C  FIND THE INTERNAL VARIABLES, WKB.
C
                  NIN = INTVAR( IEL + 1 ) - INTVAR( IEL )
                  CALL RANGE ( IEL, .FALSE., WK, WKB, NVAREL, NIN, 
     *                         ITYPEE( IEL ), NVAREL, NIN )
C
C  MULTIPLY THE INTERNAL VARIABLES BY THE ELEMENT HESSIAN AND PUT THE
C  PRODUCT IN WKC. CONSIDER THE FIRST COLUMN OF THE ELEMENT HESSIAN.
C
                  IELHST = ISTADH( IEL )
                  PI     = GI * WKB( 1 )
CDIR$ IVDEP
                  DO 220 IROW = 1, NIN
                     IJHESS   = ISYMMH( 1, IROW ) + IELHST
                     WKC( IROW ) = PI * HUVALS( IJHESS )
  220             CONTINUE
C
C  NOW CONSIDER THE REMAINING COLUMNS OF THE ELEMENT HESSIAN.
C
                  DO 240 JCOL = 2, NIN
                     PI       = GI * WKB( JCOL )
                     IF ( PI .NE. ZERO ) THEN
CDIR$ IVDEP
                        DO 230 IROW = 1, NIN
                           IJHESS  = ISYMMH( JCOL, IROW ) + IELHST
                           WKC( IROW ) = WKC( IROW ) + PI *
     *                                   HUVALS( IJHESS )
  230                   CONTINUE
                     END IF
  240             CONTINUE
C
C  SCATTER THE PRODUCT BACK ONTO THE ELEMENTAL VARIABLES, WK.
C
                  CALL RANGE ( IEL, .TRUE., WKC, WK, NVAREL, NIN, 
     *                         ITYPEE( IEL ), NIN, NVAREL )
C
C  ADD THE SCATTERED PRODUCT TO Q.
C
                  LL        = ISTAEV( IEL )
CDIR$ IVDEP
                  DO 250 II = 1, NVAREL
                     L      = IELVAR( LL )
                     Q( L ) = Q( L ) + WK( II )
                     LL     = LL + 1
  250             CONTINUE
               ELSE
C
C  THE IEL-TH ELEMENT HESSIAN HAS NO INTERNAL REPRESENTATION.
C
                  LTHVAR      = ISTAEV( IEL ) - 1
                  IELHST      = ISTADH( IEL )
                  DO 270 JCOL = 1, NVAREL
                     PI       = GI * P( IELVAR( LTHVAR + JCOL ) )
                     IF ( PI .NE. ZERO ) THEN
CDIR$ IVDEP
                        DO 260 IROW = 1, NVAREL
                           IJHESS   = ISYMMH( JCOL, IROW ) + IELHST
                           L        = IELVAR( LTHVAR + IROW )
                           Q( L )   = Q( L ) + PI * HUVALS( IJHESS )
  260                   CONTINUE
                     END IF
  270             CONTINUE
               END IF
  280       CONTINUE
         ELSE
C
C  -------------------- CASE 2. P IS SPARSE ------------------------
C
            DO 380 J = NVAR1, NVAR2
C
C  CONSIDER EACH NONZERO COMPONENT OF P SEPARATELY.
C
               I   = IVAR( J )
               IPT = ISPTRS( I )
               IF ( IPT .GE. 0 ) THEN
C
C  THE INDEX OF THE I-TH COMPONENT LIES IN THE IEL-TH NONLINEAR ELEMENT.
C
                  IELL = ISELTS( I )
  300             CONTINUE
C
C  CHECK TO ENSURE THAT THE IEL-TH ELEMENT HAS NOT ALREADY BEEN USED.
C
                  IF ( ISWKSP( IELL ) .LT. NBPROD ) THEN
                     ISWKSP( IELL ) = NBPROD
                     IEL    = IELING( IELL )
                     NVAREL = ISTAEV( IEL + 1 ) - ISTAEV( IEL )
                     IG     = ISLGRP( IELL )
                     IF ( GXEQX( IG ) ) THEN
                        GI = GSCALE( IG ) * ESCALE( IELL )
                     ELSE
                        GI = GSCALE( IG ) * ESCALE( IELL ) * GVALS2( IG) 
                     END IF
                     IF ( INTREP( IEL ) ) THEN
C
C  THE IEL-TH ELEMENT HESSIAN HAS AN INTERNAL REPRESENTATION.
C  COPY THE ELEMENTAL VARIABLES INTO WK.
C
                        LL          = ISTAEV( IEL )
CDIR$ IVDEP
                        DO 310 II   = 1, NVAREL
                           WK( II ) = P( IELVAR( LL ) )
                           LL       = LL + 1
  310                   CONTINUE
C
C  FIND THE INTERNAL VARIABLES.
C
                        NIN = INTVAR( IEL + 1 ) - INTVAR( IEL )
                        CALL RANGE ( IEL, .FALSE., WK, WKB, NVAREL, NIN,
     *                              ITYPEE( IEL ), NVAREL, NIN )
C
C  MULTIPLY THE INTERNAL VARIABLES BY THE ELEMENT HESSIAN AND PUT THE 
C  PRODUCT IN WKB. CONSIDER THE FIRST COLUMN OF THE ELEMENT HESSIAN.
C
                        IELHST = ISTADH( IEL )
                        PI     = GI * WKB( 1 )
CDIR$ IVDEP
                        DO 320 IROW = 1, NIN
                           IJHESS   = ISYMMH( 1, IROW ) + IELHST
                           WKC( IROW ) = PI * HUVALS( IJHESS )
  320                   CONTINUE
C
C  NOW CONSIDER THE REMAINING COLUMNS OF THE ELEMENT HESSIAN.
C
                        DO 340 JCOL = 2, NIN
                           PI       = GI * WKB( JCOL )
                           IF ( PI .NE. ZERO ) THEN
CDIR$ IVDEP
                              DO 330 IROW = 1, NIN
                                 IJHESS  = ISYMMH( JCOL, IROW ) + IELHST
                                 WKC( IROW ) = WKC( IROW ) + PI
     *                                     * HUVALS( IJHESS )
  330                         CONTINUE
                           END IF
  340                   CONTINUE
C
C  SCATTER THE PRODUCT BACK ONTO THE ELEMENTAL VARIABLES, WK.
C
                        CALL RANGE ( IEL, .TRUE., WKC, WK, NVAREL, NIN,
     *                               ITYPEE( IEL ), NIN, NVAREL )
C
C  ADD THE SCATTERED PRODUCT TO Q.
C
                        LL        = ISTAEV( IEL )
CDIR$ IVDEP
                        DO 350 II = 1, NVAREL
                           L      = IELVAR( LL )
C
C  IF Q HAS A NONZERO IN POSITION L, STORE ITS INDEX IN INONNZ.
C
                           IF ( ABS( WK( II ) ) .GT. TINY ) THEN
                              IF ( IUSED( L ) .EQ. 0 ) THEN
                                 Q( L )     = WK( II )
                                 IUSED( L ) = 1
                                 NNONNZ     = NNONNZ + 1
                                 INONNZ( NNONNZ ) = L
                              ELSE
                                 Q( L ) = Q( L ) + WK( II )
                              END IF
                           END IF
                           LL = LL + 1
  350                   CONTINUE
                     ELSE
C
C  THE IEL-TH ELEMENT HESSIAN HAS NO INTERNAL REPRESENTATION.
C
                        LTHVAR = ISTAEV( IEL ) - 1
                        IELHST = ISTADH( IEL )
                        DO 370 JCOL = 1, NVAREL
                           PI       = GI * P( IELVAR( LTHVAR + JCOL ) )
                           IF ( PI .NE. ZERO ) THEN
CDIR$ IVDEP
                              DO 360 IROW = 1, NVAREL
                                 IJHESS  = ISYMMH( JCOL, IROW ) + IELHST
C
C  IF Q HAS A NONZERO IN POSITION L, STORE ITS INDEX IN INONNZ.
C
                                 IF ( ABS( HUVALS( IJHESS ) )
     *                              .GT. TINY ) THEN
                                    L = IELVAR( LTHVAR + IROW )
                                    IF ( IUSED( L ) .EQ. 0 ) THEN
                                       Q( L ) = PI * HUVALS( IJHESS )
                                       IUSED( L ) = 1
                                       NNONNZ     = NNONNZ + 1
                                       INONNZ( NNONNZ ) = L
                                    ELSE
                                       Q( L ) = Q( L ) + PI *
     *                                          HUVALS( IJHESS )
                                    END IF
                                 END IF
  360                         CONTINUE
                           END IF
  370                   CONTINUE
                     END IF
                  END IF
C
C  CHECK TO SEE IF THERE ARE ANY FURTHER ELEMENTS WHOSE VARIABLES
C  INCLUDE THE I-TH VARIABLE.
C
                  IF ( IPT .GT. 0 ) THEN
                     IELL = ISELTS( IPT )
                     IPT  = ISPTRS( IPT )
                     GO TO 300
                  END IF
               END IF
  380       CONTINUE
         END IF
      END IF
C
C  ==================== THE PRODUCT IS COMPLETE =======================
C
C  RESET IUSED TO ZERO.
C
      IF ( .NOT. DENSEP ) THEN
         DO 390 I                = 1, NNONNZ
            IUSED( INONNZ( I ) ) = 0
  390    CONTINUE
      END IF
      RETURN
C
C  END OF HSPRD.
C
      END
