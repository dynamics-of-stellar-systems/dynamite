C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CSGRSH( N , M , X     , GRLAGF, LV, V , NNZJ  ,
     *                   LCJAC , CJAC  , INDVAR, INDFUN, NNZH  ,
     *                   LH, H , IRNH  , ICNH  )
      INTEGER            N , M , LV    , NNZJ  , NNZH  , LCJAC , LH
      LOGICAL            GRLAGF
      INTEGER            INDVAR( LCJAC), INDFUN( LCJAC )
      INTEGER            IRNH  ( LH   ), ICNH  ( LH )
CS    REAL               X     ( N    ), V     ( LV ),
CD    DOUBLE PRECISION   X     ( N    ), V     ( LV ),
     *                   H     ( LH   ), CJAC  ( LCJAC )
C
C  Compute the Hessian matrix of the Lagrangian function of
C  a problem initially written in Standard Input Format (SIF).
C  Also compute the Hessian matrix of the Lagrangian function of
C  the problem
C
C  CJAC	 is an array which gives the values of the nonzeros of the
C	 gradients of the objective, or Lagrangian, and general
C	 constraint functions evaluated  at X and V. The i-th entry
C	 of CJAC gives the value of the derivative with respect to
C	 variable INDVAR(i) of function INDFUN(i). INDFUN(i) = 0
C        indicates the objective function whenever GRLAGF is .FALSE.
C        or the Lagrangian function when GRLAGF is .TRUE., while
C        INDFUN(i) = j > 0 indicates the j-th general constraint
C        function.
C
C H      is an array which gives the values of entries of the
C        upper triangular part of the Hessian matrix of the
C        Lagrangian function, stored in coordinate form, i.e.,
C        the entry H(i) is the derivative with respect to variables
C        with indices X(IRNH(i)) and X(ICNH(i)) for i = 1, ...., NNZH.
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
C  integer variables from the PRFCTS common block.
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
     *                 / DIMS   /, / PRFCTS /
C
C  Local variables
C
      INTEGER            LIWKH , ICON  , LIH   , LGTEMP, IFSTAT
      INTEGER            LNXTRW, LINXTR, INFORM, IENDGV, IGSTAT
      INTEGER            I , J , IEL, K, IG, II, IG1, L, JJ, LL
      INTEGER            NIN   , NVAREL, NELOW , NELUP , ISTRGV
      LOGICAL            NONTRV
CS    EXTERNAL           SSETVL, SSETVI, RANGE 
CD    EXTERNAL           DSETVL, DSETVI, RANGE 
CS    REAL               FTT   , ONE   , ZERO  , GI    , SCALEE, GII
CD    DOUBLE PRECISION   FTT   , ONE   , ZERO  , GI    , SCALEE, GII
CS    PARAMETER        ( ZERO  = 0.0E+0, ONE = 1.0E+0 )
CD    PARAMETER        ( ZERO  = 0.0D+0, ONE = 1.0D+0 )
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
     *             3, IFSTAT )
C
C  compute the group argument values ft.
C
      DO 70 IG = 1, NG
         FTT    = - WK( B + IG )
C
C  include the contribution from the linear element.
C
         DO 30 J = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
            FTT  = FTT + WK( A + J ) * X( IWK( ICNA + J ) )
   30    CONTINUE
C
C  include the contributions from the nonlinear elements.
C
         DO 60 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
            FTT  = FTT + WK( ESCALE + J ) * FUVALS( IWK( IELING + J ) )
   60    CONTINUE
         WK( FT + IG ) = FTT
C
C  Record the derivatives of trivial groups.
C
         IF ( LOGI( GXEQX + IG ) ) THEN
            WK( GVALS +     NG + IG ) = ONE
            WK( GVALS + 2 * NG + IG ) = ZERO
         END IF
   70 CONTINUE
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
C  Define the real work space needed for ELGRD.
C  Ensure that there is sufficient space.
C
      IF ( LWK2 .LT. NG ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF
      IF ( NUMCON .GT. 0 ) THEN
C
C  Change the group weightings to include the contributions from
C  the Lagrange multipliers.
C
         DO 80 IG = 1, NG
            I     = IWK( KNDOFC + IG )
            IF ( I .EQ. 0 ) THEN
               WK( LWKSTR + IG ) = WK( GSCALE + IG )
            ELSE
               WK( LWKSTR + IG ) = WK( GSCALE + IG ) * V( I )
            END IF
   80    CONTINUE
C
C  Compute the gradient values. Initialize the gradient of the
C  objective function as zero.
C
         NNZJ     = 0
         LGTEMP   = WRK + N + MAXSEL
         DO 120 J = 1, N
            WK( LGTEMP + J ) = ZERO
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
            GI  = WK( GSCALE + IG )
            GII = WK( LWKSTR + IG )
            IF ( NONTRV ) THEN
               GI  = GI  * WK( GVALS + NG + IG )
               GII = GII * WK( GVALS + NG + IG )
            END IF
CS          CALL SSETVI( IENDGV - ISTRGV + 1, WK( WRK + 1 ),
CD          CALL DSETVI( IENDGV - ISTRGV + 1, WK( WRK + 1 ),
     *                   IWK( LSVGRP + ISTRGV ), ZERO )
C
C  This is the first gradient evaluation or the group has nonlinear
C  elements.
C
            IF ( FIRSTG .OR. NELOW .LE. NELUP ) THEN
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
                     WK( LGTEMP + LL ) = WK( LGTEMP + LL ) +
     *                                   GI * WK( WRK + LL )
C
C  The group defines a constraint.
C
                  ELSE
                     NNZJ           = NNZJ + 1
                     IF ( NNZJ .LE. LCJAC ) THEN
                        CJAC  ( NNZJ ) = GI * WK( WRK + LL )
                        INDFUN( NNZJ ) = ICON
                        INDVAR( NNZJ ) = LL
                     END IF
                     IF ( GRLAGF )
     *                  WK( LGTEMP + LL ) = WK( LGTEMP + LL ) +
     *                                      GII * WK( WRK + LL )
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
C                             linear element improved. 43 lines replace 40
C
C  Include the contribution from the linear element.
C
CDIR$ IVDEP
               DO 210 K = IWK( ISTADA + IG ),IWK( ISTADA + IG1 ) - 1
                  J             = IWK( ICNA + K )
                  WK( WRK + J ) = WK( WRK + J ) + WK( A + K )
  210          CONTINUE
C
C  Allocate a gradient.
C
CDIR$ IVDEP
               DO 220 I = ISTRGV, IENDGV
                  LL    = IWK( LSVGRP + I )
C
C  The group belongs to the objective function.
C
                  IF ( ICON .EQ. 0 ) THEN
                     WK( LGTEMP + LL ) = WK( LGTEMP + LL ) +
     *                                   GI * WK( WRK + LL )
C
C  The group defines a constraint.
C
                  ELSE
                     NNZJ = NNZJ + 1
                     IF ( NNZJ .LE. LCJAC ) THEN
                        CJAC  ( NNZJ ) = GI * WK( WRK + LL )
                        INDFUN( NNZJ ) = ICON
                        INDVAR( NNZJ ) = LL
                     END IF
                     IF ( GRLAGF )
     *                  WK( LGTEMP + LL ) = WK( LGTEMP + LL ) +
     *                                      GII * WK( WRK + LL )
                  END IF
C
C  Increment the address for the next nonzero in the column of
C  the jacobian for variable LL.
C
                  IF ( NONTRV ) THEN
                     JJ                 = IWK( LSTAJC + LL )
                     IWK( LSTAJC + LL ) = JJ + 1
                  END IF
  220          CONTINUE
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
C
C  Transfer the gradient of the objective function to the sparse
C  storage scheme.
C
         DO 310 I = 1, N
C           IF ( WK( LGTEMP + I ) .NE. ZERO ) THEN
               NNZJ           = NNZJ + 1
               IF ( NNZJ .LE. LCJAC ) THEN
                  CJAC  ( NNZJ ) = WK( LGTEMP + I )
                  INDFUN( NNZJ ) = 0
                  INDVAR( NNZJ ) = I
               END IF
C           END IF
  310    CONTINUE
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
C  Transfer the gradient of the objective function to the sparse
C  storage scheme.
C
         NNZJ = 0
         DO 400 I = 1, N
C           IF ( FUVALS( LGGFX + I ) .NE. ZERO ) THEN
               NNZJ           = NNZJ + 1
               IF ( NNZJ .LE. LCJAC ) THEN
                  CJAC  ( NNZJ ) = FUVALS( LGGFX + I )
                  INDFUN( NNZJ ) = 0
                  INDVAR( NNZJ ) = I
               END IF
C           END IF
  400    CONTINUE
      END IF
      FIRSTG = .FALSE.
C
C  Verify that the Jacobian can fit in the alloted space
C
      IF ( NNZJ .GT. LCJAC ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2020 ) NNZJ - LCJAC 
         STOP
      END IF
C
C  Define the real work space needed for ASMBL.
C  Ensure that there is sufficient space.
C
C                            for unconstrained problems.
      IF ( NUMCON .GT. 0 ) THEN
         IF ( LWK2 .LT. N + 3 * MAXSEL + NG ) THEN
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 )
            STOP
         END IF
      ELSE
         IF ( LWK2 .LT. N + 3 * MAXSEL ) THEN
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 )
            STOP
         END IF
      END IF
C
C  Define the integer work space needed for ASMBL.
C  Ensure that there is sufficient space.
C
      LIWKH  = LIWK2 - N
      LINXTR = LIWKH / 2
      IF ( LINXTR .LT. N ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF
C
C  Set starting addresses for partitions of the integer workspace.
C
      LIH    = LH
      LNXTRW = 0
      DO 320 I           = 1, N
         IWK( IVAR + I ) = I
  320 CONTINUE
C
C  Assemble the Hessian.
C
      IF ( NUMCON .GT. 0 ) THEN
CS       CALL SASMBL( N , NG, MAXSEL, N, LH    , LIH   , NNZH  ,
CD       CALL DASMBL( N , NG, MAXSEL, N, LH    , LIH   , NNZH  ,
     *             N, IWK( IVAR + 1), IWK( ISTADH + 1 ), LSTADH,
     *             IWK( ICNA + 1 )  , LICNA ,
     *             IWK( ISTADA + 1 ), LSTADA, IWK( INTVAR + 1 ), LNTVAR,
     *             IWK( IELVAR + 1 ), LELVAR, IWK( IELING + 1 ), LELING,
     *             IWK( ISTADG + 1 ), LSTADG, IWK( ISTAEV + 1 ), LSTAEV,
     *             IWK( LSTAGV + 1 ), LNSTGV, IWK( LSVGRP + 1 ), LNVGRP,
     *             IRNH, ICNH, IWK( LSEND + LNXTRW + 1 ) , LINXTR,
     *             IWK( LSEND + LIWKH + 1 ) , N,
     *             WK( A + 1 ) , LA, FUVALS, LNGUVL, FUVALS, LNHUVL,
     *             WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ),
     *             WK( LWKSTR + 1 ), WK( ESCALE + 1 ), LESCAL,
     *             H, WK ( LWKSTR + NG + 1 ), LWK2 - NG ,
     *             LOGI( GXEQX + 1 ), LGXEQX, LOGI( INTREP + 1 ),
     *             LINTRE, IWK( ITYPEE + 1 ), LINTRE,
     *             RANGE , 1, IOUT  , .FALSE., I, INFORM,
     *             .FALSE., .TRUE. )
      ELSE
CS       CALL SASMBL( N , NG, MAXSEL, N, LH    , LIH   , NNZH  ,
CD       CALL DASMBL( N , NG, MAXSEL, N, LH    , LIH   , NNZH  ,
     *             N, IWK( IVAR + 1), IWK( ISTADH + 1 ), LSTADH,
     *             IWK( ICNA + 1 )  , LICNA ,
     *             IWK( ISTADA + 1 ), LSTADA, IWK( INTVAR + 1 ), LNTVAR,
     *             IWK( IELVAR + 1 ), LELVAR, IWK( IELING + 1 ), LELING,
     *             IWK( ISTADG + 1 ), LSTADG, IWK( ISTAEV + 1 ), LSTAEV,
     *             IWK( LSTAGV + 1 ), LNSTGV, IWK( LSVGRP + 1 ), LNVGRP,
     *             IRNH, ICNH, IWK( LSEND + LNXTRW + 1 ) , LINXTR,
     *             IWK( LSEND + LIWKH + 1 ) , N,
     *             WK( A + 1 ) , LA, FUVALS, LNGUVL, FUVALS, LNHUVL,
     *             WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ),
     *             WK( GSCALE + 1 ), WK( ESCALE + 1 ), LESCAL,
     *             H, WK ( LWKSTR + 1 ), LWK2 - NG ,
     *             LOGI( GXEQX + 1 ), LGXEQX, LOGI( INTREP + 1 ),
     *             LINTRE, IWK( ITYPEE + 1 ), LINTRE,
     *             RANGE , 1, IOUT  , .FALSE., I, INFORM,
     *             .FALSE., .TRUE. )
      END IF
C
C  Check that there is sufficient integer workspace.
C
      IF ( INFORM .GT. 0 ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF
C
C  Update the counters for the report tool.
C
      NC2CG = NC2CG + PNC
      NC2OG = NC2OG + 1
      NC2OH = NC2OH + 1
      NC2CH = NC2CH + PNC
      RETURN
C
C Non-executable statements.
C
 2000 FORMAT( ' ** SUBROUTINE CSGRSH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CSGRSH: Increase the sizes of IWK and LH')
 2020 FORMAT( /, ' ** SUBROUTINE CSGRSH: array length LCJAC too small.',
     *        /, ' -- Minimization abandoned.' ,
     *        /, ' -- Increase the parameter LCJAC by at least ', I8,
     *           ' and restart.'  )
C
C  end of CSGRSH.
C
      END
