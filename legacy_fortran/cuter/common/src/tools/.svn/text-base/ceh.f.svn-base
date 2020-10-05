C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CEH   ( N , M , X , LV, V , NE, IRNHI , LIRNHI, LE    ,
     *                   IPRNHI, HI    , LHI   , IPRHI , BYROWS )
      INTEGER            N , M , LV    , NE , LE, LIRNHI, LHI 
      LOGICAL            BYROWS
      INTEGER            IRNHI ( LIRNHI )
      INTEGER            IPRNHI( LE     ), IPRHI ( LE     )
CS    REAL               X     ( N      ), V     ( LV     ), 
CD    DOUBLE PRECISION   X     ( N      ), V     ( LV     ), 
     *                   HI    ( LHI )
C
C  Compute the Hessian matrix of the Lagrangian function of
C  a problem initially written in Standard Input Format (SIF).
C
C  The matrix is represented in "finite element format", i.e., 
C
C           ne
C      H = sum H_i, 
C          i=1
C
C  where each element H_i involves a small subset of the rows of H.
C  H is stored as a list of the row indices involved in each element
C  and the upper triangle of H_i (stored by rows or columns). 
C  Specifically,
C
C  NE     (integer) number of elements
C  IRNHI  (integer array) a list of the row indices involved which each
C          element. Those for element i directly proceed those for 
C          element i+1, i = 1, ..., NE-1
C  IPRNHI (integer array) pointers to the position in IRNHI of the first 
C          row index in each element. IPRNHI(NE+1) points to the first 
C          empty location in IRPNHI
C  HI     (real array) a list of the nonzeros in the upper triangle of
C          H_i, stored by rows, or by columns, for each element. Those 
C          for element i directly proceed those for element, i+1, 
C          i = 1, ..., NE-1
C  IPRHI  (integer array) pointers to the position in HI of the first 
C          nonzero in each element. IPRHI(NE+1) points to the first 
C          empty location in HI
C  BYROWS (logical) must be set .TRUE. if the upper triangle of each H_i 
C          is to be stored by rows, and .FALSE. if it is to be stored
C          by columns.
C
C  Based on the minimization subroutine LANCELOT/SBMIN
C  by Conn, Gould and Toint.
C
C  Nick Gould, for CGT productions,
C  November 1994.
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
     *                 / DIMS   /, / PRFCTS /
C
C  Local variables
C
      INTEGER            I , J , IFSTAT, IGSTAT
      INTEGER            IG    , LIWKH , INFORM
CS    REAL               ZERO  , ONE   , FTT
CD    DOUBLE PRECISION   ZERO  , ONE   , FTT
CS    PARAMETER        ( ZERO  = 0.0E+0, ONE = 1.0E+0 )
CD    PARAMETER        ( ZERO  = 0.0D+0, ONE = 1.0D+0 )
CS    EXTERNAL           RANGE , ELFUN , GROUP , SELGRD, SASMBE
CD    EXTERNAL           RANGE , ELFUN , GROUP , DELGRD, DASMBE
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
C  evaluate the element function gradients and Hessians.
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
C  Define the real work space needed for ASMBE.
C  Ensure that there is sufficient space.
C
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
C  Define the integer work space needed for ASMBE.
C  Ensure that there is sufficient space.
C
      LIWKH = LIWK2 - N
C
C  Assemble the Hessian.
C

      IF ( NUMCON .GT. 0 ) THEN
CS       CALL SASMBE( N , NG, MAXSEL, 
CD       CALL DASMBE( N , NG, MAXSEL, 
     *                IWK( ISTADH + 1 ), LSTADH,
     *                IWK( ICNA + 1 )  , LICNA ,
     *                IWK( ISTADA + 1 ), LSTADA,
     *                IWK( INTVAR + 1 ), LNTVAR,
     *                IWK( IELVAR + 1 ), LELVAR,
     *                IWK( IELING + 1 ), LELING,
     *                IWK( ISTADG + 1 ), LSTADG,
     *                IWK( ISTAEV + 1 ), LSTAEV,
     *                IWK( LSTAGV + 1 ), LNSTGV,
     *                IWK( LSVGRP + 1 ), LNVGRP,
     *                IWK( LIWKH + 1 ) , LIWK2 - LIWKH  ,
     *                WK( A + 1 ) , LA, FUVALS, LNGUVL, FUVALS, LNHUVL,
     *                WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ),
     *                WK( LWKSTR + 1 ), WK( ESCALE + 1 ), LESCAL,
     *                WK ( LWKSTR + NG + 1 ), LWK2 - NG ,
     *                LOGI( GXEQX + 1 ), LGXEQX, LOGI( INTREP + 1 ),
     *                LINTRE, IWK( ITYPEE + 1 ), LINTRE, RANGE , NE,
     *                IRNHI , LIRNHI, IPRNHI, HI    , LHI   , IPRHI ,
     *                BYROWS, 1, IOUT  , INFORM )
      ELSE
CS       CALL SASMBE( N , NG, MAXSEL, 
CD       CALL DASMBE( N , NG, MAXSEL, 
     *                IWK( ISTADH + 1 ), LSTADH,
     *                IWK( ICNA + 1 )  , LICNA ,
     *                IWK( ISTADA + 1 ), LSTADA,
     *                IWK( INTVAR + 1 ), LNTVAR,
     *                IWK( IELVAR + 1 ), LELVAR,
     *                IWK( IELING + 1 ), LELING,
     *                IWK( ISTADG + 1 ), LSTADG,
     *                IWK( ISTAEV + 1 ), LSTAEV,
     *                IWK( LSTAGV + 1 ), LNSTGV,
     *                IWK( LSVGRP + 1 ), LNVGRP,
     *                IWK( LIWKH + 1 ) , LIWK2 - LIWKH  ,
     *                WK( A + 1 ) , LA, FUVALS, LNGUVL, FUVALS, LNHUVL,
     *                WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ),
     *                WK( GSCALE + 1 ), WK( ESCALE + 1 ), LESCAL,
     *                WK ( LWKSTR + 1 ), LWK2 - NG ,
     *                LOGI( GXEQX + 1 ), LGXEQX, LOGI( INTREP + 1 ),
     *                LINTRE, IWK( ITYPEE + 1 ), LINTRE, RANGE , NE,
     *                IRNHI , LIRNHI, IPRNHI, HI    , LHI   , IPRHI ,
     *                BYROWS, 1, IOUT  , INFORM )
      END IF
C
C  Check that there is room for the elements
C
      IF ( INFORM .GT. 0 ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2020 )
         STOP
      END IF

      write( 6, * ) 'ceh: lirnhi = ', LIRNHI
      write( 6, * ) 'ceh: irnhi = ', (IRNHI(I), I=1,LIRNHI)

C
C  Update the counters for the report tool.
C
      NC2OH = NC2OH + 1
      NC2CH = NC2CH + PNC
      RETURN
C
C Non-executable statements.
C
 2000 FORMAT( ' ** SUBROUTINE CEH: Increase the size of WK ' )
 2020 FORMAT( ' ** SUBROUTINE CEH: Increase the size of',
     *        ' IPNRHI, IPRHI, IRNHI or HI ' )
C
C  end of CEH.
C
      END



