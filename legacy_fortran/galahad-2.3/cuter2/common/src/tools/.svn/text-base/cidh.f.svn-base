C     ( Last modified on 10 Sepc 2004 at 16:45:38 )
C  Correction: 10/Sep/2004: undeclared integer variables declared
      SUBROUTINE CIDH   ( N , X , IPROB, LH1, H )
      INTEGER            N , IPROB, LH1
CS    REAL               X     ( N    ), H( LH1, N )
CD    DOUBLE PRECISION   X     ( N    ), H( LH1, N )
C
C  Compute the Hessian matrix of a specified problem function 
C  (IPROB = 0 is the objective function, while IPROB > 0 is the 
C   IPROB-th constraint) of a problem initially written in 
C  Standard Input Format (SIF).
C
C  H is a two-dimensional array which gives the value of the
C    Hessian matrix of the specified problem function evaluated at
C    X. The i,j-th component of the array will contain
C    the derivative with respect to variables X(i) and X(j).
C
C  Based on the minimization subroutine LANCELOT/SBMIN
C  by Conn, Gould and Toint.
C
C  Nick Gould, for CGT productions,
C  August 1998.
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
      INTEGER            I , J , NCALCF, NCALCG, IFSTAT, IGSTAT
      INTEGER            IG    , LIH   , LIWKH , LNXTRW, LINXTR, INFORM
      INTEGER            K , LH, LIRNH , LJCNH , NNZH  , LWKH
CS    REAL               ZERO  , ONE   , FTT
CD    DOUBLE PRECISION   ZERO  , ONE   , FTT
CS    PARAMETER        ( ZERO  = 0.0E+0, ONE = 1.0E+0 )
CD    PARAMETER        ( ZERO  = 0.0D+0, ONE = 1.0D+0 )
      EXTERNAL           RANGE 
C
C  Check input parameters.
C
      IF ( IPROB .LT. 0 ) THEN
         WRITE( IOUT, 2020 )
         STOP
      END IF
C
C  Find group index IG of constraint IPROB.
C
      IF ( IPROB .GT. 0 ) THEN
         IG = 0
         DO 10 I = 1, NG
            IF ( IWK( KNDOFC + I ) .EQ. IPROB ) THEN
              IG = I
              GO TO 20
            END IF
   10    CONTINUE
   20    CONTINUE
         IF ( IG .EQ. 0 ) THEN
            WRITE( IOUT, 2020 )
            STOP
         END IF
      END IF
C
C  Find which elements are involved in the required problem function.
C  Initialize the list to zero.
C
      DO 30 I = 1, NELNUM
        IWK( ICALCF + I ) = 0
   30 CONTINUE
C
C  If group IG is involved, record its elements.
C
      DO 50 IG = 1, NG
         IF ( IWK( KNDOFC + IG ) .EQ. IPROB ) THEN
            DO 40 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
               IWK( ICALCF + IWK( IELING + J ) ) = 1
   40       CONTINUE
         END IF
   50 CONTINUE
C
C  Only compute the elements which are involved in the required function.
C
      NCALCF = 0
      DO 80 I = 1, NELNUM
         IF ( IWK( ICALCF + I ) .EQ. 1 ) THEN
            NCALCF = NCALCF + 1
            IWK( ICALCF + NCALCF ) = I
         ELSE
C
C  If this is the first ever evaluation, initialize FUVALS.
C
            IF ( FIRSTG ) THEN
               FUVALS( I ) = ZERO
               DO 60 J = IWK( INTVAR + I ), IWK( INTVAR + I + 1 ) - 1
                  FUVALS( J ) = ZERO
   60          CONTINUE
               DO 70 J = IWK( ISTADH + I ), IWK( ISTADH + I + 1 ) - 1
                  FUVALS( J ) = ZERO
   70          CONTINUE
            END IF
         END IF
   80 CONTINUE
C
C  evaluate the element function values.
C
      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NCALCF,
     *             IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ),
     *             IWK( IELVAR + 1 ), IWK( INTVAR + 1 ),
     *             IWK( ISTADH + 1 ), IWK( ISTEP  + 1 ),
     *             IWK( ICALCF + 1 ), 
     *             LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH, 
     *             LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU, 
     *             1, IFSTAT )
C
C  evaluate the element function derivatives.
C
      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NCALCF,
     *             IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ),
     *             IWK( IELVAR + 1 ), IWK( INTVAR + 1 ),
     *             IWK( ISTADH + 1 ), IWK( ISTEP  + 1 ),
     *             IWK( ICALCF + 1 ), 
     *             LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH, 
     *             LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU, 
     *             3, IFSTAT )
C
C  Compute the list of groups involved in the required problem function.
C
      NCALCG = 0
      DO 130 IG = 1, NG
         IF ( IWK( KNDOFC + IG ) .EQ. IPROB ) THEN
            NCALCG = NCALCG + 1
            IWK( ICALCF + NCALCG ) = IG
C
C  compute the group argument values ft.
C
            FTT    = - WK( B + IG )
C
C  include the contribution from the linear element.
C
            DO 110 J = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
               FTT  = FTT + WK( A + J ) * X( IWK( ICNA + J ) )
  110       CONTINUE
C
C  include the contributions from the nonlinear elements.
C
            DO 120 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
               FTT  = FTT + WK( ESCALE + J ) * 
     *                FUVALS( IWK( IELING + J ) )
  120          CONTINUE
            WK( FT + IG ) = FTT
C
C  Record the derivatives of trivial groups.
C
            IF ( LOGI( GXEQX + IG ) ) THEN
               WK( GVALS +     NG + IG ) = ONE
               WK( GVALS + 2 * NG + IG ) = ZERO
            END IF
         ELSE
C
C  If this is the first ever evaluation, initialize GVALS.
C
            IF ( FIRSTG ) THEN
               WK( GVALS +     NG + IG ) = ONE
               WK( GVALS + 2 * NG + IG ) = ZERO
            END IF
         END IF
  130 CONTINUE
C
C  Evaluate the group derivative values.
C
      IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG,
     *      WK( FT + 1 ), WK ( GPVALU + 1 ), NCALCG,
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
         DO 140 IG = 1, NG
            I     = IWK( KNDOFC + IG )
            IF ( I .EQ. IPROB ) THEN
               WK( LWKSTR + IG ) = WK( GSCALE + IG )
            ELSE
               WK( LWKSTR + IG ) = ZERO
            END IF
  140    CONTINUE
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
     *                WK( LWKSTR + 1 ), NG,
     *                WK( ESCALE + 1 ), LESCAL, FUVALS( LGRJAC + 1 ),
     *                LNGRJC, WK( WRK + 1 ), WK( WRK + N + 1 ), MAXSEL,
     *                LOGI( GXEQX + 1 ), LGXEQX,
     *                LOGI( INTREP + 1 ), LINTRE, RANGE  )
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
      END IF
      FIRSTG = .FALSE.
C
C  Define the real work space needed for ASMBL.
C  Ensure that there is sufficient space.
C
      IF ( NUMCON .GT. 0 ) THEN
         LWKH  = LWK2 - N - 3 * MAXSEL - NG
      ELSE
         LWKH  = LWK2 - N - 3 * MAXSEL
      END IF
      IF ( LWKH .LE. 0 ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF
C
C  Define the integer work space needed for ASMBL.
C  Ensure that there is sufficient space.
C
      LIWKH  = LIWK2 - N
      LH     = MIN( LWKH, ( LIWKH - 3 * N ) / 4 )
      LINXTR = LH + N
      IF ( LH .LE. 0 ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF
C
C  Set starting addresses for partitions of the integer workspace.
C
      LIH    = LH
      LIRNH  = 0
      LJCNH  = LIRNH  + LIH
      LNXTRW = LJCNH  + LIH
      DO 90 I            = 1, N
         IWK( IVAR + I ) = I
   90 CONTINUE
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
     *             IWK( LSEND + LIRNH + 1 ) , IWK ( LSEND + LJCNH + 1 ),
     *             IWK( LSEND + LNXTRW + 1 ) , LINXTR,
     *             IWK( LSEND + LIWKH + 1 )  , N ,
     *             WK( A + 1 ) , LA, FUVALS, LNGUVL, FUVALS, LNHUVL,
     *             WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ),
     *             WK( LWKSTR + 1 ), WK( ESCALE + 1 ), LESCAL,
     *             WK( LWKSTR + NG + 1 ), WK( LWKSTR + LWKH + NG + 1 ),
     *             LWK2 - LWKH    , LOGI( GXEQX + 1 ),
     *             LGXEQX, LOGI( INTREP + 1 ), LINTRE,
     *             IWK( ITYPEE + 1 ), LINTRE,
     *             RANGE , 1, IOUT    , .FALSE., I, INFORM, .FALSE.,
     *             .FALSE. )
      ELSE
CS       CALL SASMBL( N , NG, MAXSEL, N, LH    , LIH   , NNZH  ,
CD       CALL DASMBL( N , NG, MAXSEL, N, LH    , LIH   , NNZH  ,
     *             N, IWK( IVAR + 1), IWK( ISTADH + 1 ), LSTADH,
     *             IWK( ICNA + 1 )  , LICNA ,
     *             IWK( ISTADA + 1 ), LSTADA, IWK( INTVAR + 1 ), LNTVAR,
     *             IWK( IELVAR + 1 ), LELVAR, IWK( IELING + 1 ), LELING,
     *             IWK( ISTADG + 1 ), LSTADG, IWK( ISTAEV + 1 ), LSTAEV,
     *             IWK( LSTAGV + 1 ), LNSTGV, IWK( LSVGRP + 1 ), LNVGRP,
     *             IWK( LSEND + LIRNH + 1 ) , IWK ( LSEND + LJCNH + 1 ),
     *             IWK( LSEND + LNXTRW + 1 ) , LINXTR,
     *             IWK( LSEND + LIWKH + 1 )  , N ,
     *             WK( A + 1 ) , LA, FUVALS, LNGUVL, FUVALS, LNHUVL,
     *             WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ),
     *             WK( GSCALE + 1 ), WK( ESCALE + 1 ), LESCAL,
     *             WK( LWKSTR + 1 ), WK( LWKSTR + LWKH + 1 ),
     *             LWK2 - LWKH    , LOGI( GXEQX + 1 ),
     *             LGXEQX, LOGI( INTREP + 1 ), LINTRE,
     *             IWK( ITYPEE + 1 ), LINTRE,
     *             RANGE , 1, IOUT    , .FALSE., I, INFORM, .FALSE.,
     *             .FALSE. )
      END IF
C
C  Check that there is sufficient integer workspace.
C
      IF ( INFORM .GT. 0 ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF
C
C  Initialize the dense matrix.
C
      DO 220 J = 1, N
         DO 210 I = 1, N
            H( I, J ) = ZERO
  210    CONTINUE
  220 CONTINUE
C
C  Transfer the matrix from co-ordinate to dense storage and
C  symmetrize the martix.
C
      IF ( NUMCON .GT. 0 ) THEN
         DO 230 K     = 1, NNZH
            I         = IWK( LSEND + LIRNH + K )
            J         = IWK( LSEND + LJCNH + K )
            H( I, J ) = WK( LWKSTR + NG + K )
            H( J, I ) = WK( LWKSTR + NG + K )
  230    CONTINUE
      ELSE
         DO 240 K     = 1, NNZH
            I         = IWK( LSEND + LIRNH + K )
            J         = IWK( LSEND + LJCNH + K )
            H( I, J ) = WK( LWKSTR + K )
            H( J, I ) = WK( LWKSTR + K )
  240    CONTINUE
      END IF
C
C  Update the counters for the report tool.
C
      IF( IPROB .EQ. 0 ) THEN
         NC2OH = NC2OH + 1
      ELSE 
         NC2CH = NC2CH + 1
      ENDIF 
      RETURN
C
C Non-executable statements.
C
 2000 FORMAT( ' ** SUBROUTINE CIDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CIDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE CIDH: invalid problem index IPROB ' )
C
C  end of CIDH.
C
      END
