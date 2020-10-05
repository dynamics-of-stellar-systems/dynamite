C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CCFSG ( N, M, X, LC, C, NNZJ, LCJAC, CJAC,
     *                   INDVAR, INDFUN, GRAD )
      INTEGER            N, M, LC, NNZJ, LCJAC
      LOGICAL            GRAD
      INTEGER            INDVAR( LCJAC ), INDFUN( LCJAC )
CS    REAL               X( N ), C( LC ), CJAC  ( LCJAC )
CD    DOUBLE PRECISION   X( N ), C( LC ), CJAC  ( LCJAC )
C
C  Compute the values of the constraint functions and their gradients
C  for constraints initially written in Standard Input Format (SIF).
C  The Jacobian must be stored in a sparse format.
C  (Subroutine CCFG performs the same calculations for a dense Jacobian.)
C
C  CJAC  is an array which gives the values of the nonzeros of the
C        general constraint functions evaluated at X and V.
C        The i-th entry of CJAC gives the value of the derivative
C        with respect to variable INDVAR(i) of constraint function 
C        INDFUN(i) (i.e., INDFUN(i) = j > 0 indicates the j-th
C        general constraint function).
C
C  Based on the subroutines cfn.f and csgr.f by Nick Gould, which are
C  in turn based on the subroutine SBMIN by Conn, Gould and Toint.
C
C  Ingrid Bongartz
C  April 1992.
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
      COMMON / GLOBAL /  IWK   , WK    , FUVALS, LOGI,
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
      INTEGER            NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      INTEGER            IOUT
      COMMON / OUTPUT /  IOUT
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / DIMS /,
     *                 / PRFCTS /, / OUTPUT /
C
C  local variables.
C
      INTEGER            I , J , IEL, K, IG, II, IG1, L, LL, ICON, ICNT
      INTEGER            NIN   , NVAREL, NELOW , NELUP , ISTRGV, IENDGV
      INTEGER            LLO   , LLWRK , IFSTAT, IGSTAT
      LOGICAL            NONTRV
CS    EXTERNAL           SSETVL, SSETVI, RANGE 
CD    EXTERNAL           DSETVL, DSETVI, RANGE 
CS    REAL               FTT   , ONE   , ZERO  , GI    , SCALEE 
CD    DOUBLE PRECISION   FTT   , ONE   , ZERO  , GI    , SCALEE
CS    PARAMETER        ( ZERO  = 0.0E+0, ONE  = 1.0E+0 )
CD    PARAMETER        ( ZERO  = 0.0D+0, ONE  = 1.0D+0 )
C
      IF ( NUMCON .EQ. 0 ) RETURN
C
C  Must identify which elements are included in constraints.
C  Use logical work vector to keep track of elements already included.
C  First ensure there is sufficient room in LOGI.
C
      LLO = GXEQX + NGNG
      LLWRK = LLOGIC - LLO
      IF ( LLWRK .LT. NELNUM ) THEN
          IF ( IOUT .GT. 0 ) WRITE( IOUT, 2010 ) NELNUM - LLWRK 
          STOP
      END IF
      DO 410 I = 1, NELNUM
         LOGI( LLO + I ) = .FALSE.
  410 CONTINUE
C
C  Now identify elements in first M constraint groups.
C
      ICNT = 0
      DO 10 IG = 1, NG
         ICON = IWK( KNDOFC + IG )
         IF ( ICON .GT. 0 .AND. ICON .LE. M ) THEN
            NELOW  = IWK( ISTADG + IG  )
            NELUP  = IWK( ISTADG + IG + 1 ) - 1
            DO 20 II = NELOW, NELUP
               IEL   = IWK( IELING + II )
               IF ( .NOT. LOGI( LLO + IEL ) ) THEN
                  LOGI( LLO + IEL ) = .TRUE.
                  ICNT = ICNT + 1
                  IWK( ICALCF + ICNT ) = IEL
               END IF
   20       CONTINUE
         END IF
   10 CONTINUE
C
C  Evaluate the element function values.
C
      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), ICNT,
     *             IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ),
     *             IWK( IELVAR + 1 ), IWK( INTVAR + 1 ),
     *             IWK( ISTADH + 1 ), IWK( ISTEP  + 1 ),
     *             IWK( ICALCF + 1 ), 
     *             LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH, 
     *             LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU, 
     *             1, IFSTAT )
      IF ( GRAD ) THEN
C
C  Evaluate the element function derivatives.
C
         CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), ICNT,
     *                IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ),
     *                IWK( IELVAR + 1 ), IWK( INTVAR + 1 ),
     *                IWK( ISTADH + 1 ), IWK( ISTEP  + 1 ),
     *                IWK( ICALCF + 1 ), 
     *                LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH, 
     *                LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU, 
     *                2, IFSTAT )
      END IF
C
C  Compute the group argument values ft.
C
      DO 100 IG = 1, NG
         FTT    = ZERO
C
C  Consider only those groups in the constraints.
C
         ICON = IWK( KNDOFC + IG )
         IF ( ICON .GT. 0 .AND. ICON .LE. M ) THEN
            FTT = - WK( B + IG )
C
C  Include contributions from the linear element 
C  only if the variable belongs to the first N variables.
C
            DO 30 I = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
               J = IWK( ICNA + I )
               IF ( J .LE. N )
     *            FTT  = FTT + WK( A + I ) * X( J )
   30       CONTINUE
C
C  Include the contributions from the nonlinear elements.
C
            DO 60 I = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
               FTT  = FTT +
     *                WK( ESCALE + I ) * FUVALS( IWK( IELING + I ) )
   60       CONTINUE
C
C  Record derivatives of trivial groups.
C
            IF ( LOGI( GXEQX + IG ) ) WK( GVALS + NG + IG ) = ONE 
         END IF
         WK( FT + IG ) = FTT
  100 CONTINUE
C
C  Compute the group function values.
C
C  All group functions are trivial.
C
      IF ( ALTRIV ) THEN
CS       CALL SCOPY( NG, WK( FT + 1 ), 1, WK( GVALS + 1 ), 1 )
CD       CALL DCOPY( NG, WK( FT + 1 ), 1, WK( GVALS + 1 ), 1 )
CS       CALL SSETVL( NG, WK( GVALS + NG + 1 ), 1, ONE )
CD       CALL DSETVL( NG, WK( GVALS + NG + 1 ), 1, ONE )
      ELSE
C
C  Evaluate the group function values.
C  Evaluate groups belonging to the first M constraints only.
C
         ICNT = 0
         DO 400 IG = 1, NG
            ICON   = IWK( KNDOFC + IG )
            IF ( ICON .GT. 0 .AND. ICON .LE. M ) THEN
               ICNT = ICNT + 1
               IWK( ICALCF + ICNT ) = IG
            END IF 
  400    CONTINUE
         CALL GROUP ( WK ( GVALS  + 1 ), NG, WK( FT + 1 ),
     *                WK ( GPVALU + 1 ), ICNT,
     *                IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ),
     *                IWK( ICALCF + 1 ),
     *                LCALCG, NG1   , LCALCG, LCALCG, LGPVLU,
     *                .FALSE., IGSTAT )
      END IF
C
C  Compute the constraint function values.
C
      DO 110 IG = 1, NG
         I = IWK( KNDOFC + IG )
         IF ( I .GT. 0 .AND. I .LE. M ) THEN
            IF ( LOGI( GXEQX  + IG ) ) THEN
               C( I ) = WK( GSCALE + IG ) * WK( FT + IG )
            ELSE 
               C( I ) = WK( GSCALE + IG ) * WK( GVALS + IG )
            END IF
         END IF
  110 CONTINUE
      IF ( GRAD ) THEN
C
C  Evaluate the group derivative values.
C
         IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG,
     *         WK( FT + 1 ), WK ( GPVALU + 1 ), ICNT,
     *         IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ),
     *         IWK( ICALCF + 1 ),
     *         LCALCG, NG1   , LCALCG, LCALCG, LGPVLU,
     *         .TRUE., IGSTAT )
C
C  Compute the gradient values.  Initialize the Jacobian as zero.
C
         NNZJ = 0 
         DO 120 J = 1, LCJAC
            CJAC( J ) = ZERO
  120    CONTINUE
C
C  Consider the IG-th group.
C
         DO 290 IG = 1, NG
            ICON   = IWK( KNDOFC + IG )
C
C  Consider only those groups in the first M constraints.
C
            IF ( ICON .EQ. 0 .OR. ICON .GT. M ) GO TO 290
            IG1    = IG + 1
            ISTRGV = IWK( LSTAGV + IG  )
            IENDGV = IWK( LSTAGV + IG1 ) - 1
            NELOW  = IWK( ISTADG + IG  )
            NELUP  = IWK( ISTADG + IG1 ) - 1
            NONTRV = .NOT. LOGI( GXEQX + IG )
C
C  Compute the first derivative of the group.
C
            GI = WK( GSCALE + IG )
            IF ( NONTRV ) GI  = GI  * WK( GVALS + NG + IG )
CS          CALL SSETVI( IENDGV - ISTRGV + 1, WK( WRK + 1 ),
CD          CALL DSETVI( IENDGV - ISTRGV + 1, WK( WRK + 1 ),
     *                   IWK( LSVGRP + ISTRGV ), ZERO )
C
C  The group has nonlinear elements.
C
            IF ( NELOW .LE. NELUP ) THEN
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
C  Include contributions from the first N variables only.
C
                  IF ( LL .LE. N ) THEN
                     NNZJ           = NNZJ + 1
                     IF ( NNZJ .LE. LCJAC ) THEN
                        CJAC  ( NNZJ ) = GI * WK( WRK + LL )
                        INDFUN( NNZJ ) = ICON
                        INDVAR( NNZJ ) = LL
                     END IF
                  END IF
  190          CONTINUE
C
C  The group has only linear elements.
C
            ELSE
C                             linear element improved. 26 lines replace 19
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
C  Include contributions from the first N variables only.
C
                  IF ( LL .LE. N ) THEN
                     NNZJ           = NNZJ + 1
                     IF ( NNZJ .LE. LCJAC ) THEN
                        CJAC  ( NNZJ ) = GI * WK( WRK + LL )
                        INDFUN( NNZJ ) = ICON
                        INDVAR( NNZJ ) = LL
                     END IF
                  END IF
  220          CONTINUE
            END IF
  290    CONTINUE
C
C  Verify that the Jacobian can fit in the allotted space
C
         IF ( NNZJ .GT. LCJAC ) THEN
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 ) NNZJ - LCJAC 
            STOP
         END IF
      END IF
C
C  Update the counters for the report tool.
C
      NC2CF = NC2CF + PNC
      IF ( GRAD ) NC2CG = NC2CG + PNC
      RETURN
C
C  Non-executable statements.
C
 2000 FORMAT( /  ' ** SUBROUTINE CCFSG: array length LCJAC too small.'
     *        /  ' -- Minimization abandoned.' ,
     *        /  ' -- Increase the parameter LCJAC by at least ', I8,
     *           ' and restart.'  )
 2010 FORMAT( /  ' ** SUBROUTINE CCFSG: array length LLOGIC too small.'
     *        /  ' -- Minimization abandoned.' ,
     *        /  ' -- Increase the parameter LLOGIC by at least ', I8,
     *           ' and restart.'  )
C
C  end of CSCFG.
C
      END

