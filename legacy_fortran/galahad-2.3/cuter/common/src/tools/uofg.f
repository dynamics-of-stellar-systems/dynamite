C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UOFG  ( N, X, F, G, GRAD )
      INTEGER            N
CS    REAL               F
CD    DOUBLE PRECISION   F
CS    REAL               X( N ), G( N )
CD    DOUBLE PRECISION   X( N ), G( N )
      LOGICAL            GRAD
C
C  Compute the value of the objective function and its gradient
C  for a function initially written in Standard Input Format (SIF).
C
C  G     is an array which gives the value of the gradient of
C        the objective function evaluated at X.
C        G(i) gives the partial derivative of the objective
C        function with respect to variable X(i).
C
C  Based on the subroutines ufn.f and ugr.f by Nick Gould, which are
C  in turn based on the subroutine SBMIN by Conn, Gould and Toint.
C
C  Ingrid Bongartz 
C  February 1993.
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
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / PRFCTS /
C
C  local variables.
C
      INTEGER            I , J , IG    , IFSTAT, IGSTAT
CS    REAL               FTT   , SDOT  , ONE   , ZERO
CD    DOUBLE PRECISION   FTT   , DDOT  , ONE   , ZERO
CS    PARAMETER        ( ZERO  = 0.0E+0, ONE  = 1.0E+0 ) 
CD    PARAMETER        ( ZERO  = 0.0D+0, ONE  = 1.0D+0 ) 
      EXTERNAL           RANGE 
C
C  There are non-trivial group functions.
C
      DO 10 I = 1, MAX( NELNUM, NG )
        IWK( ICALCF + I ) = I
   10 CONTINUE
C
C  Evaluate the element function values.
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
C  Compute the group argument values ft.
C
      DO 100 IG = 1, NG
         FTT     = - WK( B + IG )
C
C  Include the contribution from the linear element 
C  only if the variable belongs to the first N variables.
C
         DO 30 I = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
            J = IWK( ICNA + I ) 
            IF ( J .LE. N ) 
     *         FTT  = FTT + WK( A + I ) * X( J )
   30    CONTINUE
C
C  Include the contributions from the nonlinear elements.
C
         DO 60 I = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
            FTT  = FTT + 
     *             WK( ESCALE + I ) * FUVALS( IWK( IELING + I ) )
   60    CONTINUE
         WK( FT + IG ) = FTT
C
C  Record the derivatives of trivial groups.
C
         IF ( LOGI( GXEQX + IG ) ) WK( GVALS + NG + IG ) = ONE
  100 CONTINUE
C
C  Compute the group function values.
C
C  All group functions are trivial.
C
      IF ( ALTRIV ) THEN
CS       F = SDOT( NG, WK( GSCALE + 1 ), 1, WK( FT + 1 ), 1 )
CD       F = DDOT( NG, WK( GSCALE + 1 ), 1, WK( FT + 1 ), 1 )
CS       CALL SCOPY( NG, WK( FT + 1 ), 1, WK( GVALS + 1 ), 1 )
CD       CALL DCOPY( NG, WK( FT + 1 ), 1, WK( GVALS + 1 ), 1 )
CS       CALL SSETVL( NG, WK( GVALS + NG + 1 ), 1, ONE )
CD       CALL DSETVL( NG, WK( GVALS + NG + 1 ), 1, ONE )
      ELSE
C
C  Evaluate the group function values.
C
         CALL GROUP ( WK ( GVALS  + 1 ), NG, WK( FT + 1 ),
     *                WK ( GPVALU + 1 ), NG,
     *                IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ),
     *                IWK( ICALCF + 1 ),
     *                LCALCG, NG1   , LCALCG, LCALCG, LGPVLU,
     *                .FALSE., IGSTAT )
C
C  Compute the objective function value.
C
         F         = ZERO
         DO 110 IG = 1, NG
            IF ( LOGI( GXEQX + IG ) ) THEN
               F = F + WK( GSCALE + IG ) * WK( FT + IG )
            ELSE
               F = F + WK( GSCALE + IG ) * WK( GVALS + IG )
            END IF
  110    CONTINUE
      END IF
      IF ( GRAD ) THEN
C
C  Evaluate the element function derivatives.
C
         CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM,
     *                IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ),
     *                IWK( IELVAR + 1 ), IWK( INTVAR + 1 ),
     *                IWK( ISTADH + 1 ), IWK( ISTEP  + 1 ),
     *                IWK( ICALCF + 1 ), 
     *                LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH, 
     *                LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU, 
     *                2, IFSTAT )
C
C  Evaluate the group derivative values.
C
         IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG, 
     *         WK( FT + 1 ), WK ( GPVALU + 1 ), NG,
     *         IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ),
     *         IWK( ICALCF + 1 ),
     *         LCALCG, NG1   , LCALCG, LCALCG, LGPVLU,
     *         .TRUE., IGSTAT )
C
C  Compute the gradient values. 
C
CS       CALL SELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA,
CD       CALL DELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA,
     *                IWK( ISTADA + 1 ), LSTADA, IWK( IELING + 1 ),
     *                LELING, IWK( ISTADG + 1 ), LSTADG,
     *                IWK( ITYPEE + 1 ), LINTRE,
     *                IWK( ISTAEV + 1 ), LSTAEV, IWK( IELVAR + 1 ),
     *                LELVAR, IWK( INTVAR + 1 ), LNTVAR, IWK(LSVGRP+1),
     *                LNVGRP, IWK( LSTAJC + 1 ), LNSTJC,
     *                IWK( LSTAGV + 1 ), LNSTGV, WK( A + 1 ), LA,
     *                WK( GVALS + NG + 1 ), LGVALS,
C ** from LGGFX to LNGUVL to stop zero length arrays on VAX.
     *                FUVALS, LNGUVL, FUVALS( LGGFX + 1 ),
     *                WK( GSCALE + 1 ), LGSCAL,
     *                WK( ESCALE + 1 ), LESCAL, FUVALS( LGRJAC + 1 ),
     *                LNGRJC, WK( WRK + 1 ), WK( WRK + N + 1 ), MAXSEL,
     *                LOGI( GXEQX + 1 ), LGXEQX,
     *                LOGI( INTREP + 1 ), LINTRE, RANGE  )
         FIRSTG = .FALSE.
C
C  Store the gradient value.
C
         DO 300 I = 1, N
            G( I ) = FUVALS( LGGFX + I )
  300    CONTINUE
      END IF
C
C  Update the counters for the report tool.
C
      NC2OF = NC2OF + 1
      IF( GRAD ) NC2OG = NC2OG + 1
      RETURN
C
C  end of UOFG.
C
      END
