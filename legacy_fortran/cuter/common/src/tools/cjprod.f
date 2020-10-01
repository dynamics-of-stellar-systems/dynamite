C     ( Last modified on 10 Sepc 2004 at 16:55:38 )
C  Correction: 10/Sep/2004: undeclared integers variables declared
      SUBROUTINE CJPROD( N , M , GOTJ  , JTRANS, X , V , LV, R, LR )
      INTEGER           N , M , LV, LR
      LOGICAL           GOTJ, JTRANS
CS    REAL              X( N ), V( LV ), R( LR )
CD    DOUBLE PRECISION  X( N ), V( LV ), R( LR )
C
C  Compute the matrix-vector product between the Jacobian matrix
C  of the constraints (JTRANS = .FALSE.), or its transpose 
C  (JTRANS = .TRUE.) for the problem, and a given vector P. 
C  The result is placed in R. If GOTJ is .TRUE. the first derivatives 
C  are assumed to have already been computed. If the user is unsure, 
C  set GOTJ = .FALSE. the first time a product is required with the 
C  Jacobian evaluated at X. X is not used if GOTJ = .TRUE.
C
C  Nick Gould, for GOT/CUTEr productions.
C  June, 2003.
C
      INTEGER            LIWK  , LWK   , LFUVAL, LLOGIC, LCHARA
C
C  ---------------------------------------------------------------------
C
C  Parameters whose value might be changed by the user:
C
C  The following parameters define the sizes of problem
C  dependent arrays. These may be changed by the user to
C  suit a particular problem or system configuration.
C
C  The TOOLS will issue error messages if any of these sizes
C  is too small, telling which parameter to increase.
C
C  ---------------------------------------------------------------------
C
C#{sizing}
C
      INTEGER                IWK( LIWK    )
      LOGICAL              LOGI ( LLOGIC  )
      CHARACTER * 10       CHA  ( LCHARA  )
CS    REAL                   WK ( LWK     )
CD    DOUBLE PRECISION       WK ( LWK     )
CS    REAL              FUVALS  ( LFUVAL  )
CD    DOUBLE PRECISION  FUVALS  ( LFUVAL  )
C
C  ---------------------------------------------------------------------
C
C  End of parameters which might be changed by the user.
C
C  ---------------------------------------------------------------------
C
C  Integer variables from the GLOBAL common block.
C
      INTEGER            NG    , NELNUM, NGEL  , NVARS , NNZA  , NGPVLU
      INTEGER            NEPVLU, NG1   , NEL1  , ISTADG, ISTGP , ISTADA
      INTEGER            ISTAEV, ISTEP , ITYPEG, KNDOFC, ITYPEE
      INTEGER            IELING, IELVAR, ICNA  , ISTADH, INTVAR, IVAR
      INTEGER            ICALCF, ITYPEV, IWRK  , A     , B
      INTEGER            U     , GPVALU, EPVALU
      INTEGER            ESCALE, GSCALE, VSCALE, GVALS , XT    , DGRAD
      INTEGER            Q     , WRK   , INTREP, GXEQX , GNAMES, VNAMES
      INTEGER            LO    , CH    , LIWORK, LWORK , NGNG  , FT
      INTEGER            LA, LB, NOBJGR, LU, LELVAR
      INTEGER            LSTAEV, LSTADH, LNTVAR, LCALCF
      INTEGER            LELING, LINTRE, LFT, LGXEQX, LSTADG, LGVALS
      INTEGER            LICNA , LSTADA, LKNDOF, LGPVLU, LEPVLU
      INTEGER            LGSCAL, LESCAL, LVSCAL, LCALCG
C
C  Integer variables from the LOCAL common block.
C
      INTEGER            LFXI  , LGXI  , LHXI  , LGGFX , LDX   , LGRJAC
      INTEGER            LQGRAD, LBREAK, LP    , LXCP  , LX0   , LGX0
      INTEGER            LDELTX, LBND  , LWKSTR, LSPTRS, LSELTS, LINDEX
      INTEGER            LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC
      INTEGER            LNNONZ, LNONZ2, LSYMMD, LSYMMH
      INTEGER            LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND
      INTEGER            LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV
      INTEGER            LNSTJC, LNIUSE, LNFREC, LNNNON, LNNNO2, LNSYMD
      INTEGER            LNSYMH, LNLGRP, LNVGRP, LNGCLJ, LNVLJR, LNQGRD
      INTEGER            LNBRAK, LNP   , LNBND , LNFXI , LNGXI , LNGUVL
      INTEGER            LNHXI , LNHUVL, LNGGFX, LNDX  , LNGRJC, LIWK2
      INTEGER            LWK2  , MAXSIN, NINVAR, MAXSEL
      INTEGER            NTYPE , NSETS , LSTYPE, LSSWTR, LSSIWT, LSIWTR
      INTEGER            LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR
      INTEGER            LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      LOGICAL            ALTRIV, FIRSTG
C
C  Integer variables from the PRFCTS common block.
C
      INTEGER            NC2OF , NC2OG , NC2OH, NC2CF, NC2CG, NC2CH
      INTEGER            NHVPR , PNC
      REAL               SUTIME, STTIME
      COMMON / GLOBAL /  IWK   , WK    , FUVALS, LOGI  ,
     *                   NG    , NELNUM, NGEL  , NVARS , NNZA  , NGPVLU,
     *                   NEPVLU, NG1   , NEL1  , ISTADG, ISTGP , ISTADA,
     *                   ISTAEV, ISTEP , ITYPEG, KNDOFC, ITYPEE,
     *                   IELING, IELVAR, ICNA  , ISTADH, INTVAR, IVAR  ,
     *                   ICALCF, ITYPEV, IWRK  , A     , B     ,
     *                   U     , GPVALU, EPVALU,
     *                   ESCALE, GSCALE, VSCALE, GVALS , XT    , DGRAD ,
     *                   Q     , WRK   , INTREP, GXEQX , GNAMES, VNAMES,
     *                   LO    , CH    , LIWORK, LWORK , NGNG  , FT    ,
     *                   ALTRIV, FIRSTG,
     *                   LA, LB, NOBJGR, LU, LELVAR,
     *                   LSTAEV, LSTADH, LNTVAR, LCALCF,
     *                   LELING, LINTRE, LFT, LGXEQX, LSTADG, LGVALS,
     *                   LICNA , LSTADA, LKNDOF, LGPVLU, LEPVLU,
     *                   LGSCAL, LESCAL, LVSCAL, LCALCG
      COMMON / CHARA /   CHA
      COMMON / LOCAL /   LFXI  , LGXI  , LHXI  , LGGFX , LDX   , LGRJAC,
     *                   LQGRAD, LBREAK, LP    , LXCP  , LX0   , LGX0  ,
     *                   LDELTX, LBND  , LWKSTR, LSPTRS, LSELTS, LINDEX,
     *                   LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC,
     *                   LNNONZ, LNONZ2, LSYMMD, LSYMMH,
     *                   LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND ,
     *                   LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV,
     *                   LNSTJC, LNIUSE, LNFREC, LNNNON, LNNNO2, LNSYMD,
     *                   LNSYMH, LNLGRP, LNVGRP, LNGCLJ, LNVLJR, LNQGRD,
     *                   LNBRAK, LNP   , LNBND , LNFXI , LNGXI , LNGUVL,
     *                   LNHXI , LNHUVL, LNGGFX, LNDX  , LNGRJC, LIWK2 ,
     *                   LWK2  , MAXSIN, NINVAR, MAXSEL, NTYPE ,
     *                   NSETS , LSTYPE, LSSWTR, LSSIWT, LSIWTR,
     *                   LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR,
     *                   LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      INTEGER            IOUT
      COMMON / OUTPUT /  IOUT
      COMMON / PRFCTS /  NC2OF , NC2OG , NC2OH, NC2CF, NC2CG, NC2CH,
     *                   NHVPR , PNC   , SUTIME, STTIME
      INTEGER            NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / OUTPUT /,
     *                 / DIMS   /, / PRFCTS /
C
C  Local variables
C
      INTEGER            I , IG, J , ICON, K, IG1, II
      INTEGER            L, IEL, NVAREL , NIN
      INTEGER            IFSTAT, IGSTAT
CS    REAL               ZERO  , ONE   , FTT, PROD, SCALEE
CD    DOUBLE PRECISION   ZERO  , ONE   , FTT, PROD, SCALEE
CS    PARAMETER        ( ZERO  = 0.0E+0, ONE = 1.0E+0 )
CD    PARAMETER        ( ZERO  = 0.0D+0, ONE = 1.0D+0 )
      EXTERNAL           RANGE 
      IF ( NUMCON .EQ. 0 ) RETURN
C
C  Check input data.
C
      IF ( ( JTRANS .AND. LV .LT. M ) .OR. 
     *     ( .NOT. JTRANS .AND. LV .LT. N ) ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF
      IF ( ( JTRANS .AND. LR .LT. N ) .OR. 
     *     ( .NOT. JTRANS .AND. LR .LT. M ) ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2020 )
         STOP
      END IF
C
C  There are non-trivial group functions.
C
      IF ( .NOT. GOTJ ) THEN
         DO 10 I = 1, MAX( NELNUM, NG )
           IWK( ICALCF + I ) = I
   10    CONTINUE
C
C  Evaluate the element function values.
C
         CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM,
     *                IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ),
     *                IWK( IELVAR + 1 ), IWK( INTVAR + 1 ),
     *                IWK( ISTADH + 1 ), IWK( ISTEP  + 1 ),
     *                IWK( ICALCF + 1 ), 
     *                LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH, 
     *                LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU, 
     *                1, IFSTAT )
C
C  Evaluate the element function values.
C
         CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM,
     *                IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ),
     *                IWK( IELVAR + 1 ), IWK( INTVAR + 1 ),
     *                IWK( ISTADH + 1 ), IWK( ISTEP  + 1 ),
     *                IWK( ICALCF + 1 ), 
     *                LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH, 
     *                LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU, 
     *                3, IFSTAT )
C
C  Compute the group argument values ft.
C
         DO 70 IG = 1, NG
            FTT    = - WK( B + IG )
C
C  Include the contribution from the linear element.
C
            DO 30 J = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
               FTT  = FTT + WK( A + J ) * X( IWK( ICNA + J ) )
   30       CONTINUE
C
C  Include the contributions from the nonlinear elements.
C
            DO 60 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
               FTT  = FTT + WK( ESCALE + J ) * FUVALS( IWK( IELING + J))
   60       CONTINUE
            WK( FT + IG ) = FTT
C
C  Record the derivatives of trivial groups.
C
            IF ( LOGI( GXEQX + IG ) ) THEN
               WK( GVALS + NG + IG ) = ONE
            END IF
   70    CONTINUE
C
C  Evaluate the group derivative values.
C
         IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG,
     *         WK( FT + 1 ), WK ( GPVALU + 1 ), NG,
     *         IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ),
     *         IWK( ICALCF + 1 ),
     *         LCALCG, NG1   , LCALCG, LCALCG, LGPVLU,
     *         .TRUE., IGSTAT )
      END IF
C
C  Ensure that there is sufficient space.
C
      IF ( LWK2 .LT. N ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF
C
C  Form the product r = J(transpose) v
C
      IF ( JTRANS ) THEN
C
C  Initialize R
C
         DO 110 I = 1, N
            R( I ) = ZERO
  110    CONTINUE
C
C  Consider the IG-th group.
C
         DO 190 IG = 1, NG
            ICON   = IWK( KNDOFC + IG )
            IF ( ICON .GT. 0 ) THEN
               IG1 = IG + 1

C
C  Compute the product of v(i) with the (scaled) group derivative
C
               IF ( LOGI( GXEQX + IG ) ) THEN
                  PROD = V( ICON ) * WK( GSCALE + IG )
               ELSE
                  PROD = V( ICON ) * WK( GSCALE + IG )
     *                             * WK( GVALS + NG + IG )
               END IF
C
C  Loop over the group's nonlinear elements.
C
               DO 150 II = IWK( ISTADG + IG  ), 
     *                     IWK( ISTADG + IG1 ) - 1
                  IEL    = IWK( IELING + II  )
                  K      = IWK( INTVAR + IEL )
                  L      = IWK( ISTAEV + IEL )
                  NVAREL = IWK( ISTAEV + IEL + 1 ) - L
                  SCALEE = WK( ESCALE + II ) * PROD
                  IF ( LOGI( INTREP + IEL ) ) THEN
C
C  The IEL-th element has an internal representation.
C
                     NIN = IWK( INTVAR + IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., FUVALS( K ),
     *                            WK( WRK + 1 ), NVAREL, NIN,
     *                            IWK( ITYPEE + IEL ),
     *                            NIN, NVAREL )
CDIR$ IVDEP
                     DO 130 I = 1, NVAREL
                        J = IWK( IELVAR + L )
                        R( J ) = R( J ) + SCALEE * WK( WRK + I )
                        L      = L + 1
  130                CONTINUE
                  ELSE
C
C  The IEL-th element has no internal representation.
C
CDIR$ IVDEP
                     DO 140 I = 1, NVAREL
                        J = IWK( IELVAR + L )
                        R( J ) = R( J ) + SCALEE * FUVALS( K )
                        K      = K + 1
                        L      = L + 1
  140                CONTINUE
                  END IF
  150          CONTINUE
C
C  Include the contribution from the linear element.
C
CDIR$ IVDEP
               DO 160 K = IWK( ISTADA + IG  ),
     *                    IWK( ISTADA + IG1 ) - 1
                  J = IWK( ICNA + K )
                  R( J ) = R( J ) + WK( A + K ) * PROD
  160          CONTINUE
            END IF
  190    CONTINUE
C
C  Form the product r = J v
C
      ELSE
C
C  Consider the IG-th group.
C
         DO 290 IG = 1, NG
            ICON   = IWK( KNDOFC + IG )
            IF ( ICON .GT. 0 ) THEN
               IG1  = IG + 1
               PROD = ZERO
C
C  Compute the first derivative of the group.
C
C  Loop over the group's nonlinear elements.
C
               DO 250 II = IWK( ISTADG + IG  ), 
     *                     IWK( ISTADG + IG1 ) - 1
                  IEL    = IWK( IELING + II  )
                  K      = IWK( INTVAR + IEL )
                  L      = IWK( ISTAEV + IEL )
                  NVAREL = IWK( ISTAEV + IEL + 1 ) - L
                  SCALEE = WK( ESCALE + II )
                  IF ( LOGI( INTREP + IEL ) ) THEN
C
C  The IEL-th element has an internal representation.
C
                     NIN = IWK( INTVAR + IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., FUVALS( K ),
     *                            WK( WRK + 1 ), NVAREL, NIN,
     *                            IWK( ITYPEE + IEL ),
     *                            NIN, NVAREL )
CDIR$ IVDEP
                     DO 230 I = 1, NVAREL
                        PROD  = PROD + V( IWK( IELVAR + L ) ) 
     *                               * SCALEE * WK( WRK + I )
                        L     = L + 1
  230                CONTINUE
                  ELSE
C
C  The IEL-th element has no internal representation.
C
CDIR$ IVDEP
                     DO 240 I = 1, NVAREL
                        PROD  = PROD + V( IWK( IELVAR + L ) ) 
     *                           * SCALEE * FUVALS( K )
                        K     = K + 1
                        L     = L + 1
  240                CONTINUE
                  END IF
  250          CONTINUE
C
C  Include the contribution from the linear element.
C
CDIR$ IVDEP
               DO 260 K = IWK( ISTADA + IG  ),
     *                    IWK( ISTADA + IG1 ) - 1
                  PROD = PROD + V( IWK( ICNA + K ) ) * WK( A + K )
  260          CONTINUE
C
C  Multiply the product by the (scaled) group derivative
C
               IF ( LOGI( GXEQX + IG ) ) THEN
                  R( ICON ) = PROD * WK( GSCALE + IG )
               ELSE
                  R( ICON ) = PROD * WK( GSCALE + IG )
     *                             * WK( GVALS + NG + IG )
               END IF
            END IF
  290    CONTINUE
      END IF
C
C  Update the counters for the report tool.
C
      IF ( .NOT. GOTJ ) THEN
         NC2OG = NC2OG + 1
         NC2CG = NC2CG + PNC
      END IF
      RETURN
C
C Non-executable statements.
C
 2000 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of V ' )
 2020 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of R ' )
C
C  end of CJPROD.
C
      END
