C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CGR   ( N , M , X     , GRLAGF, LV, V ,
     *                   G     , JTRANS, LCJAC1, LCJAC2, CJAC  )
      INTEGER            N , M , LV    , LCJAC1, LCJAC2
      LOGICAL            GRLAGF, JTRANS
CS    REAL               X( N ), G( N ), V( LV ), CJAC( LCJAC1, LCJAC2 )
CD    DOUBLE PRECISION   X( N ), G( N ), V( LV ), CJAC( LCJAC1, LCJAC2 )
C
C  Compute both the gradients of the objective, or Lagrangian, and
C  general constraint functions of a problem initially written in
C  Standard Input Format (SIF).
C
C  G	 is an array which gives the value of the gradient of
C	 the objective function evaluated at X (GRLAGF = .FALSE.)
C        of of the Lagrangian function evaluated at X and V
C        (GRLAGF = .TRUE.),
C
C  CJAC	 is a two-dimensional array of dimension ( LCJAC1, LCJAC2 )
C	 which gives the value of the Jacobian matrix of the
C	 constraint functions, or its transpose, evaluated at X.
C	 If JTRANS is .TRUE., the i,j-th component of the array
C        will contain the i-th derivative of the j-th constraint
C        function. Otherwise, if JTRANS is .FALSE., the i,j-th
C        component of the array will contain the j-th derivative
C        of the i-th constraint function.
C
C  Based on the minimization subroutine LANCELOT/SBMIN
C  by Conn, Gould and Toint.
C
C  Nick Gould, for CGT productions,
C  November 1991.
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
C  integer variables from the GLOBAL common block.
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
C  integer variables from the LOCAL common block.
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
      COMMON / PRFCTS /  NC2OF , NC2OG , NC2OH, NC2CF, NC2CG, NC2CH,
     *                   NHVPR , PNC   , SUTIME, STTIME
      INTEGER            IOUT
      COMMON / OUTPUT /  IOUT
      INTEGER            NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / OUTPUT /,
     *                 / DIMS   /
C
C  local variables.
C
      INTEGER            I , J , IEL, K, IG, II, IG1, L, JJ, LL, ICON
      INTEGER            NIN   , NVAREL, NELOW , NELUP , ISTRGV, IENDGV
      INTEGER            IFSTAT, IGSTAT
      LOGICAL            NONTRV
CS    EXTERNAL           SSETVL, SSETVI, RANGE 
CD    EXTERNAL           DSETVL, DSETVI, RANGE 
CS    REAL               FTT   , ONE   , ZERO  , GI    , SCALEE, GII
CD    DOUBLE PRECISION   FTT   , ONE   , ZERO  , GI    , SCALEE, GII
CS    PARAMETER        ( ZERO  = 0.0E+0, ONE  = 1.0E+0 )
CD    PARAMETER        ( ZERO  = 0.0D+0, ONE  = 1.0D+0 )
C
C  Check input parameters.
C
C                            dimension-checking.
      IF ( NUMCON .GT. 0 ) THEN
         IF ( JTRANS ) THEN
            IF ( LCJAC1 .LT. N .OR. LCJAC2 .LT. M ) THEN
               IF ( LCJAC1 .LT. N ) WRITE( IOUT, 2000 )
               IF ( LCJAC2 .LT. M ) WRITE( IOUT, 2010 )
               STOP
            END IF
         ELSE
            IF ( LCJAC1 .LT. M .OR. LCJAC2 .LT. N ) THEN
               IF ( LCJAC1 .LT. M ) WRITE( IOUT, 2000 )
               IF ( LCJAC2 .LT. N ) WRITE( IOUT, 2010 )
               STOP
            END IF
         END IF
      END IF
C
C  there are non-trivial group functions.
C
      DO 10 I = 1, MAX( NELNUM, NG )
        IWK( ICALCF + I ) = I
   10 CONTINUE
C
C  evaluate the element function values.
C
      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM,
     *             IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ),
     *             IWK( IELVAR + 1 ), IWK( INTVAR + 1 ),
     *             IWK( ISTADH + 1 ), IWK( ISTEP  + 1 ),
     *             IWK( ICALCF + 1 ), 
     *             LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH, 
     *             LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU, 
     *             1, IFSTAT )
C
C  evaluate the element function values.
C
      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM,
     *             IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ),
     *             IWK( IELVAR + 1 ), IWK( INTVAR + 1 ),
     *             IWK( ISTADH + 1 ), IWK( ISTEP  + 1 ),
     *             IWK( ICALCF + 1 ), 
     *             LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH, 
     *             LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU, 
     *             2, IFSTAT )
C
C  compute the group argument values ft.
C
      DO 40 IG = 1, NG
         FTT   = - WK( B + IG )
C
C  include the contribution from the linear element.
C
         DO 20 J = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
            FTT  = FTT + WK( A + J ) * X( IWK( ICNA + J ) )
   20    CONTINUE
C
C  include the contributions from the nonlinear elements.
C
         DO 30 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
            FTT  = FTT + WK( ESCALE + J ) * FUVALS( IWK( IELING + J ) )
   30    CONTINUE
         WK( FT + IG ) = FTT
C
C  Record the derivatives of trivial groups.
C
         IF ( LOGI( GXEQX + IG ) ) WK( GVALS + NG + IG ) = ONE
   40 CONTINUE
C
C  evaluate the group derivative values.
C
      IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG,
     *      WK( FT + 1 ), WK ( GPVALU + 1 ), NG,
     *      IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ),
     *      IWK( ICALCF + 1 ),
     *      LCALCG, NG1   , LCALCG, LCALCG, LGPVLU,
     *      .TRUE., IGSTAT )
C
C  For unconstrained problems, skip construction of gradient 
C  and Jacobian. Call ELGRD instead.
C
      IF ( NUMCON .GT. 0 ) THEN
C
C  Compute the gradient values. Initialize the gradient and
C  Jacobian (or its transpose) as zero.
C
         DO 120 J = 1, N
            G( J ) = ZERO
            DO 110 I = 1, M
               IF ( JTRANS ) THEN
                  CJAC( J, I ) = ZERO
               ELSE
                  CJAC( I, J ) = ZERO
               END IF
  110       CONTINUE
  120    CONTINUE
C
C  Consider the IG-th group.
C
         DO 290 IG = 1, NG
            IG1    = IG + 1
            ICON   = IWK( KNDOFC + IG )
            ISTRGV = IWK( LSTAGV + IG  )
            IENDGV = IWK( LSTAGV + IG1 ) - 1
            NELOW  = IWK( ISTADG + IG  )
            NELUP  = IWK( ISTADG + IG1 ) - 1
            NONTRV = .NOT. LOGI( GXEQX + IG )
C
C  Compute the first derivative of the group.
C
            GI = WK( GSCALE + IG )
            IF ( ICON .EQ. 0 ) THEN
               GII = GI
            ELSE
               IF ( GRLAGF ) GII = GI * V( IWK( KNDOFC + IG  ) )
            END IF
            IF ( NONTRV ) THEN
               GI  = GI  * WK( GVALS + NG + IG )
               IF ( GRLAGF ) GII = GII * WK( GVALS + NG + IG )
            END IF
C
C  This is the first gradient evaluation or the group has nonlinear
C  elements.
C
            IF ( FIRSTG .OR. NELOW .LE. NELUP ) THEN
CS             CALL SSETVI( IENDGV - ISTRGV + 1, WK( WRK + 1 ),
CD             CALL DSETVI( IENDGV - ISTRGV + 1, WK( WRK + 1 ),
     *                      IWK( LSVGRP + ISTRGV ), ZERO )
C
C  Loop over the group's nonlinear elements.
C
               DO 150 II  = NELOW, NELUP
                  IEL    = IWK( IELING + II      )
                  K      = IWK( INTVAR + IEL     )
                  L      = IWK( ISTAEV + IEL     )
                  NVAREL = IWK( ISTAEV + IEL + 1 ) - L
                  SCALEE = WK( ESCALE + II )
                  IF ( LOGI( INTREP + IEL ) ) THEN
C
C  The IEL-th element has an internal representation.
C
                     NIN = IWK( INTVAR + IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., FUVALS( K ),
     *                            WK( WRK + N + 1 ), NVAREL, NIN,
     *                            IWK( ITYPEE + IEL ),
     *                            NIN, NVAREL )
CDIR$ IVDEP
                     DO 130 I         = 1, NVAREL
                        J             = IWK( IELVAR + L )
                        WK( WRK + J ) = WK( WRK + J ) +
     *                                     SCALEE * WK( WRK + N + I )
                        L             = L + 1
  130                CONTINUE
                  ELSE
C
C  The IEL-th element has no internal representation.
C
CDIR$ IVDEP
                     DO 140 I         = 1, NVAREL
                        J             = IWK( IELVAR + L )
                        WK( WRK + J ) = WK( WRK + J ) +
     *                                     SCALEE * FUVALS( K )
                        K             = K + 1
                        L             = L + 1
  140                CONTINUE
                  END IF
  150          CONTINUE
C
C  Include the contribution from the linear element.
C
CDIR$ IVDEP
               DO 160 K         = IWK( ISTADA + IG  ),
     *                            IWK( ISTADA + IG1 ) - 1
                  J             = IWK( ICNA + K )
                  WK( WRK + J ) = WK( WRK + J ) + WK( A + K )
  160          CONTINUE
C
C  Allocate a gradient.
C
CDIR$ IVDEP
               DO 190 I = ISTRGV, IENDGV
                  LL    = IWK( LSVGRP + I )
C
C  The group belongs to the objective function.
C
                  IF ( ICON .EQ. 0 ) THEN
                     G( LL ) = G( LL ) +  GI * WK( WRK + LL )
C
C  The group defines a constraint.
C
                  ELSE
                     IF ( JTRANS ) THEN
                        CJAC( LL, ICON ) =   GI * WK( WRK + LL )
                     ELSE
                        CJAC( ICON, LL ) =   GI * WK( WRK + LL )
                     END IF
                     IF ( GRLAGF )
     *                  G( LL ) = G( LL ) + GII * WK( WRK + LL )
                  END IF
C
C  If the group is non-trivial, also store the nonzero entries of the
C  gradient of the function in GRJAC.
C
                  IF ( NONTRV ) THEN
                     JJ                    = IWK( LSTAJC + LL )
                     FUVALS( LGRJAC + JJ ) = WK ( WRK    + LL )
C
C  Increment the address for the next nonzero in the column of
C  the jacobian for variable LL.
C
                     IWK( LSTAJC + LL ) = JJ + 1
                  END IF
  190          CONTINUE
C
C  This is not the first gradient evaluation and there is only a linear
C  element.
C
            ELSE
C
C  Allocate a gradient.
C
CDIR$ IVDEP
               DO 210 K = IWK( ISTADA + IG ), IWK( ISTADA + IG1 ) - 1
                  LL    = IWK( ICNA + K )
C
C  The group belongs to the objective function.
C
                  IF ( ICON .EQ. 0 ) THEN
                     G( LL ) = G( LL ) +  GI * WK( A + K )
C
C  The group defines a constraint.
C
                  ELSE
                     IF ( JTRANS ) THEN
                        CJAC( LL, ICON )   =            GI * WK( A + K )
                     ELSE
                        CJAC( ICON, LL )   =            GI * WK( A + K )
                     END IF
                     IF ( GRLAGF ) G( LL ) = G( LL ) + GII * WK( A + K )
                  END IF
  210          CONTINUE
C
C  The group is non-trivial; increment the starting addresses for
C  the groups used by each variable in the (unchanged) linear
C  element to avoid resetting the nonzeros in the jacobian.
C
               IF ( NONTRV ) THEN
CDIR$ IVDEP
                  DO 220 I              = ISTRGV, IENDGV
                     LL                 = IWK( LSVGRP + I  )
                     IWK( LSTAJC + LL ) = IWK( LSTAJC + LL ) + 1
  220             CONTINUE
               END IF
            END IF
  290    CONTINUE
C
C  Reset the starting addresses for the lists of groups using
C  each variable to their values on entry.
C
         DO 300 I             = N, 2, - 1
            IWK( LSTAJC + I ) = IWK( LSTAJC + I - 1 )
  300    CONTINUE
         IWK( LSTAJC + 1 ) = 1
      ELSE
C
C  Compute the gradient value.
C
CS       CALL SELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA,
CD       CALL DELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA,
     *                IWK( ISTADA + 1 ), LSTADA, IWK( IELING + 1 ),
     *                LELING, IWK( ISTADG + 1 ), LSTADG,
     *                IWK( ITYPEE + 1 ), LINTRE,
     *                IWK( ISTAEV + 1 ), LSTAEV, IWK( IELVAR + 1 ),
     *                LELVAR, IWK( INTVAR + 1 ), LNTVAR,
     *                IWK( LSVGRP + 1 ),
     *                LNVGRP, IWK( LSTAJC + 1 ), LNSTJC,
     *                IWK( LSTAGV + 1 ), LNSTGV, WK( A + 1 ), LA,
     *                WK( GVALS + NG + 1 ), LGVALS,
     *                FUVALS, LNGUVL, FUVALS( LGGFX + 1 ),
     *                WK( GSCALE + 1 ), LGSCAL,
     *                WK( ESCALE + 1 ), LESCAL, FUVALS( LGRJAC + 1 ),
     *                LNGRJC, WK( WRK + 1 ), WK( WRK + N + 1 ), MAXSEL,
     *                LOGI( GXEQX + 1 ), LGXEQX,
     *                LOGI( INTREP + 1 ), LINTRE, RANGE  )
C
C  Store the gradient value.
C
         DO 400 I = 1, N
            G( I ) = FUVALS( LGGFX + I )
  400    CONTINUE
      END IF
      FIRSTG = .FALSE.
C
C  Update the counters for the report tool.
C
      NC2OG = NC2OG + 1
      NC2CG = NC2CG + PNC
      RETURN
C
C Non-executable statements.
C
 2000 FORMAT( ' ** SUBROUTINE CGR: Increase the leading dimension',
     *        ' of CJAC ' )
 2010 FORMAT( ' ** SUBROUTINE CGR: Increase the second dimension',
     *        ' of CJAC ' )
C
C  end of CGR.
C
      END
