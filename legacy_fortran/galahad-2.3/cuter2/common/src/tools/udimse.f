C     ( Last modified on 10 Sepc 2004 at 16:45:38 )
C  Correction: 10/Sep/2004: undeclared integer variable declared
      SUBROUTINE UDIMSE( NE, NZH, NZIRNH )
      INTEGER            NE, NZH, NZIRNH
C
C  Compute the number of elements and the space required to store
C  the Hessian matrix of the objective function of a problem
C  initially written in Standard Input Format (SIF).
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
C
C  NE     (integer) number of elements
C  NZH    (integer) number of entries needed to store the real values of
C                   H. Specifically, the sum of the number of entries in
C                   the upper triangle of each H_i.
C  NZIRNH (integer) number of entries needed to store the integer entries
C                   of H. Specifically, the sum of the row dimensions of 
C                   each H_i.
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
      INTEGER            IG    , NVARG , IG1
C
C  Initilaize counts
C
      NE       = 0
      NZH      = 0
      NZIRNH   = 0
C
C  Loop over the groups
C
      DO 10 IG = 1, NG
         IG1   = IG + 1
C
C  Only consider nonlinear groups
C
         IF ( IWK( ISTADG + IG ) .LT. IWK( ISTADG + IG1 ) .OR. 
     *        .NOT. LOGI( GXEQX + IG ) ) THEN
            NE = NE + 1
            NVARG = IWK( LSTAGV + IG1 ) - IWK( LSTAGV + IG )
            NZIRNH = NZIRNH + NVARG
            NZH = NZH + ( NVARG * ( NVARG + 1 ) ) / 2
         END IF
   10 CONTINUE
      RETURN
C
C  end of UDIMSE.
C
      END
