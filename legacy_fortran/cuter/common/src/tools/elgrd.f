C     ( Last modified on 23 Dec 2000 at 22:01:38 )
C
C  ** FOR THE CRAY 2, LINES STARTING CDIR$ IVDEP TELL THE COMPILER TO
C     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.
C
CS    SUBROUTINE SELGRD( N , NG, FIRSTG, ICNA  , LICNA , ISTADA, LSTADA,
CD    SUBROUTINE DELGRD( N , NG, FIRSTG, ICNA  , LICNA , ISTADA, LSTADA,
     *                   IELING, LELING, ISTADG, LSTADG, ITYPEE, LITYPE,
     *                   ISTAEV, LSTAEV, IELVAR, LELVAR, INTVAR, LNTVAR, 
     *                   ISVGRP, LNVGRP, ISTAJC, LNSTJC, ISTAGV, LNSTGV, 
     *                   A , LA, GVALS2, LGVALS, GUVALS, LGUVAL, GRAD  ,
     *                   GSCALE, LGSCAL, ESCALE, LESCAL, GRJAC , LNGRJC, 
     *                   WKPVAR, WKEVAR, LNWKEV, GXEQX , LGXEQX, INTREP, 
     *                   LINTRE, RANGE  )
      INTEGER          N , NG, LA    , LITYPE
      INTEGER          LICNA , LSTADA, LELING, LSTADG, LSTAEV, LELVAR
      INTEGER          LNTVAR, LNVGRP, LNSTJC, LNSTGV, LGVALS, LNWKEV
      INTEGER          LGSCAL, LESCAL, LGXEQX, LINTRE, LGUVAL, LNGRJC
      LOGICAL          FIRSTG
      INTEGER          ICNA( LICNA ), ISTADA( LSTADA ), IELING( LELING )
      INTEGER          ISTADG( LSTADG ), ISTAEV( LSTAEV )
      INTEGER          IELVAR( LELVAR ), INTVAR( LNTVAR )
      INTEGER          ISVGRP( LNVGRP ), ISTAJC( LNSTJC )
      INTEGER          ISTAGV( LNSTGV ), ITYPEE( LITYPE )
CS    REAL             A( LA ), GVALS2( LGVALS ), GUVALS( LGUVAL ), 
CD    DOUBLE PRECISION A( LA ), GVALS2( LGVALS ), GUVALS( LGUVAL ), 
     *                 GRAD( N ), GSCALE( LGSCAL ), ESCALE( LESCAL ), 
     *                 GRJAC( LNGRJC ), WKPVAR( N ), WKEVAR( LNWKEV )
      LOGICAL          GXEQX( LGXEQX ), INTREP( LINTRE )
      EXTERNAL         RANGE 
C
C  CALCULATE THE GRADIENT OF EACH GROUP, GRJAC, AND THE GRADIENT OF THE
C  OBJECTIVE FUNCTION, GRAD.
C
C  NICK GOULD, 4TH JULY 1990.
C  FOR CGT PRODUCTIONS.
C
C  LOCAL VARIABLES.
C
      INTEGER          I, IEL, IG, IG1, II,J, JJ, K, L, LL
      INTEGER          NIN, NVAREL, NELOW, NELUP, ISTRGV, IENDGV
CS    REAL             ZERO, GI, SCALEE
CD    DOUBLE PRECISION ZERO, GI, SCALEE
      LOGICAL          NONTRV
CS    EXTERNAL         SSETVL, SSETVI
CD    EXTERNAL         DSETVL, DSETVI
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0 )
C
C  INITIALIZE THE GRADIENT AS ZERO.
C
CS    CALL SSETVL( N, GRAD, 1, ZERO )
CD    CALL DSETVL( N, GRAD, 1, ZERO )
C
C  CONSIDER THE IG-TH GROUP.
C
      DO 160 IG = 1, NG
         IG1    = IG + 1
         ISTRGV = ISTAGV( IG )
         IENDGV = ISTAGV( IG1 ) - 1
         NELOW  = ISTADG( IG )
         NELUP  = ISTADG( IG1 ) - 1
         NONTRV = .NOT. GXEQX( IG )
C
C  COMPUTE THE FIRST DERIVATIVE OF THE GROUP.
C
         IF ( NONTRV ) THEN
            GI = GSCALE( IG ) * GVALS2( IG ) 
         ELSE
            GI = GSCALE( IG )
         END IF 
C
C  THIS IS THE FIRST GRADIENT EVALUATION OR THE GROUP HAS NONLINEAR
C  ELEMENTS.
C
         IF ( FIRSTG .OR. NELOW .LE. NELUP ) THEN
CS          CALL SSETVI( IENDGV - ISTRGV + 1, WKPVAR, ISVGRP( ISTRGV ),
CD          CALL DSETVI( IENDGV - ISTRGV + 1, WKPVAR, ISVGRP( ISTRGV ),
     *                   ZERO )
C
C  LOOP OVER THE GROUP'S NONLINEAR ELEMENTS.
C
            DO 30 II  = NELOW, NELUP
               IEL    = IELING( II )
               K      = INTVAR( IEL )
               L      = ISTAEV( IEL )
               NVAREL = ISTAEV( IEL + 1 ) - L
               SCALEE = ESCALE( II )
               IF ( INTREP( IEL ) ) THEN
C
C  THE IEL-TH ELEMENT HAS AN INTERNAL REPRESENTATION.
C
                  NIN = INTVAR( IEL + 1 ) - K
                  CALL RANGE ( IEL, .TRUE., GUVALS( K ),
     *                         WKEVAR, NVAREL, NIN, 
     *                         ITYPEE( IEL ), NIN, NVAREL )
CDIR$ IVDEP
                  DO 10 I        = 1, NVAREL
                     J           = IELVAR( L )
                     WKPVAR( J ) = WKPVAR( J ) + SCALEE * WKEVAR( I )
                     L           = L + 1
   10             CONTINUE
               ELSE
C
C  THE IEL-TH ELEMENT HAS NO INTERNAL REPRESENTATION.
C
CDIR$ IVDEP
                  DO 20 I        = 1, NVAREL
                     J           = IELVAR( L )
                     WKPVAR( J ) = WKPVAR( J ) + SCALEE * GUVALS( K )
                     K           = K + 1
                     L           = L + 1
   20             CONTINUE
               END IF
   30       CONTINUE
C
C  INCLUDE THE CONTRIBUTION FROM THE LINEAR ELEMENT.
C
CDIR$ IVDEP
            DO 40 K        = ISTADA( IG ), ISTADA( IG1 ) - 1
               J           = ICNA( K )
               WKPVAR( J ) = WKPVAR( J ) + A( K )
   40       CONTINUE
C
C  FIND THE GRADIENT OF THE GROUP.
C
            IF ( NONTRV ) THEN
C
C  THE GROUP IS NON-TRIVIAL.
C
CDIR$ IVDEP
               DO 50 I       = ISTRGV, IENDGV
                  LL         = ISVGRP( I )
                  GRAD( LL ) = GRAD( LL ) + GI * WKPVAR( LL )
C
C  AS THE GROUP IS NON-TRIVIAL, ALSO STORE THE NONZERO ENTRIES OF THE
C  GRADIENT OF THE FUNCTION IN GRJAC.
C
                  JJ           = ISTAJC( LL )
                  GRJAC( JJ )  = WKPVAR( LL )
C
C  INCREMENT THE ADDRESS FOR THE NEXT NONZERO IN THE COLUMN OF
C  THE JACOBIAN FOR VARIABLE LL.
C
                  ISTAJC( LL ) = JJ + 1
   50          CONTINUE
            ELSE
C
C  THE GROUP IS TRIVIAL.
C
CDIR$ IVDEP
               DO 60 I       = ISTRGV, IENDGV
                  LL         = ISVGRP( I )
                  GRAD( LL ) = GRAD( LL ) +  GI * WKPVAR( LL )
   60          CONTINUE
            END IF
C
C  THIS IS NOT THE FIRST GRADIENT EVALUATION AND THERE IS ONLY A LINEAR
C  ELEMENT.
C
         ELSE
C
C  ADD THE GRADIENT OF THE LINEAR ELEMENT TO THE OVERALL GRADIENT
C
CDIR$ IVDEP
            DO 130 K      = ISTADA( IG ), ISTADA( IG1 ) - 1
               LL         = ICNA( K )
               GRAD( LL ) = GRAD( LL ) + GI * A( K )
  130       CONTINUE
C
C  THE GROUP IS NON-TRIVIAL; INCREMENT THE STARTING ADDRESSES FOR
C  THE GROUP  USED BY EACH VARIABLE IN THE (UNCHANGED) LINEAR
C  ELEMENT TO AVOID RESETTING THE NONZEROS IN THE JACOBIAN.
C
            IF ( NONTRV ) THEN
CDIR$ IVDEP
               DO 140 I        = ISTRGV, IENDGV
                  LL           = ISVGRP( I )
                  ISTAJC( LL ) = ISTAJC( LL ) + 1
  140          CONTINUE
            END IF
         END IF
  160 CONTINUE
C
C  RESET THE STARTING ADDRESSES FOR THE LISTS OF GROUP  USING
C  EACH VARIABLE TO THEIR VALUES ON ENTRY.
C
      DO 170 I       = N, 2, - 1
         ISTAJC( I ) = ISTAJC( I - 1 )
  170 CONTINUE
      ISTAJC( 1 ) = 1
      RETURN
C
C  END OF ELGRD.
C
      END
