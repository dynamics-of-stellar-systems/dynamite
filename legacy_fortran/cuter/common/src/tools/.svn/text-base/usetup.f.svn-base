C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE USETUP( INPUT , IOUT  , N , X , BL, BU, NMAX )
      INTEGER            INPUT , IOUT  , N     , NMAX
CS    REAL               X( NMAX ), BL( NMAX ), BU( NMAX )
CD    DOUBLE PRECISION   X( NMAX ), BL( NMAX ), BU( NMAX )
C
C  Set up the input data for the remaining unconstrained
C  optimization tools.
C
C  Nick Gould, for CGT productions,
C  30th October, 1991.
C
C
C  Workspace arrays.
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
      INTEGER            LNNONZ, LNONZ2, LSYMMD, LSYMMH, NELTYP, NGRTYP
      INTEGER            LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND
      INTEGER            LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV
      INTEGER            LNSTJC, LNIUSE, LNFREC, LNNNON, LNNNO2, LNSYMD
      INTEGER            LNSYMH, LNLGRP, LNVGRP, LNGCLJ, LNVLJR, LNQGRD
      INTEGER            LNBRAK, LNP   , LNBND , LNFXI , LNGXI , LNGUVL
      INTEGER            LNHXI , LNHUVL, LNGGFX, LNDX  , LNGRJC, LIWK2
      INTEGER            LWK2  , MAXSIN, NINVAR, MAXSEL, LNIWTR
      INTEGER            NTYPE , NSETS , LSTYPE, LSSWTR, LSSIWT
      INTEGER            LSIWTR, LSWTRA, LNTYPE, LNSWTR, LNSIWT
      INTEGER            LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      LOGICAL            ALTRIV, FIRSTG
C
C  variables from the PRFCTS common block
C
      INTEGER            NC2OF , NC2OG , NC2OH,  NC2CF,  NC2CG,  NC2CH
      INTEGER            NHVPR , PNC
      REAL               SUTIME, STTIME
C
C  the common blocks
C
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
      INTEGER            IOUT2
      COMMON / OUTPUT /  IOUT2
      COMMON / PRFCTS /  NC2OF , NC2OG , NC2OH, NC2CF, NC2CG, NC2CH,
     *                   NHVPR , PNC   , SUTIME, STTIME
      INTEGER            NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / OUTPUT /,
     *                 / DIMS /, / PRFCTS /
C
C  local variables.
C
      INTEGER            IALGOR, IPRINT, INFORM, I, NSLACK
      LOGICAL            FDGRAD, DEBUG
      REAL               DUM,    CPUTIM
CS    REAL               OBFBND( 2 )
CD    DOUBLE PRECISION   OBFBND( 2 )
      CHARACTER * 8      PNAME
      CHARACTER * 10     CHTEMP
      EXTERNAL           RANGE , CPUTIM
      SUTIME = CPUTIM( DUM )
      IOUT2  = IOUT
      DEBUG  = .FALSE.
      DEBUG  = DEBUG .AND. IOUT .GT. 0
      IPRINT = 0
      IF ( DEBUG ) IPRINT = 3
C
C  Input the problem dimensions.
C
      READ( INPUT, 1001 ) N , NG, NELNUM, NGEL  , NVARS , NNZA , NGPVLU,
     *                    NEPVLU, NELTYP, NGRTYP
      IF ( N .LE. 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2020 )
         STOP
      END IF
      IF ( NG .LE. 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2030 )
         STOP
      END IF
      IF ( N .GT. NMAX ) THEN
         CLOSE( INPUT )
         IF ( IOUT .GT. 0 ) THEN
            WRITE( IOUT, 2000 ) 'X   ', 'NMAX  ', N - NMAX
            WRITE( IOUT, 2000 ) 'BL  ', 'NMAX  ', N - NMAX
            WRITE( IOUT, 2000 ) 'BU  ', 'NMAX  ', N - NMAX
         END IF
         STOP
      END IF
C
C  Input the problem type.
C
      READ( INPUT, 1000 ) IALGOR, PNAME
C
C  Set useful integer values.
C
      NG1    = NG     + 1
      NGNG   = NG     + NG
      NEL1   = NELNUM + 1
C
C  Partition the integer workspace.
C
      ISTADG = 0
      ISTGP  = ISTADG + NG1
      ISTADA = ISTGP  + NG1
      ISTAEV = ISTADA + NG1
      ISTEP  = ISTAEV + NEL1
      ITYPEG = ISTEP  + NEL1
      KNDOFC = ITYPEG + NG
      ITYPEE = KNDOFC + NG
      IELING = ITYPEE + NELNUM
      IELVAR = IELING + NGEL
      ICNA   = IELVAR + NVARS
      ISTADH = ICNA   + NNZA
      INTVAR = ISTADH + NEL1
      IVAR   = INTVAR + NEL1
      ICALCF = IVAR   + N
      ITYPEV = ICALCF + MAX( NELNUM, NG )
      IWRK   = ITYPEV + N
      LIWORK = LIWK   - IWRK
C
C  Ensure there is sufficient room.
C
      IF ( LIWORK .LT. 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 )
     *       'IWK   ', 'LIWK  ', - LIWORK
         STOP
      END IF
C
C  Partition the real workspace.
C
      A      = 0
      B      = A      + NNZA
      U      = B      + NG
      GPVALU = U      + NG
      EPVALU = GPVALU + NGPVLU
      ESCALE = EPVALU + NEPVLU
      GSCALE = ESCALE + NGEL
      VSCALE = GSCALE + NG
      GVALS  = VSCALE + N
      XT     = GVALS  + 3 * NG
      DGRAD  = XT     + N
      Q      = DGRAD  + N
      FT     = Q      + N
      WRK    = FT     + NG
      LWORK  = LWK    - WRK
C
C  Ensure there is sufficient room.
C
      IF ( LWORK .LT. 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 )
     *       'WK   ', 'LWK   ', - LWORK
         STOP
      END IF
C
C  Partition the logical workspace.
C
      INTREP = 0
      GXEQX  = INTREP + NELNUM
      LO     = GXEQX  + NGNG
C
C  Ensure there is sufficient room.
C
      IF ( LLOGIC .LT. LO ) THEN
         CLOSE( INPUT )
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 )
     *       'LOGI  ', 'LLOGIC', LO - LLOGIC
         STOP
      END IF
C
C  Partition the character workspace.
C
      GNAMES = 0
      VNAMES = GNAMES + NG
      CH     = VNAMES + N
C
C  Ensure there is sufficient room.
C
      IF ( LCHARA .LT. CH + 1 ) THEN
         CLOSE( INPUT )
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 )
     *       'CHA   ', 'LCHARA', CH + 1 - LCHARA
         STOP
      END IF
C
C  Record the lengths of arrays.
C
      LSTADG = MAX( 1, NG1 )
      LSTADA = MAX( 1, NG1 )
      LSTAEV = MAX( 1, NEL1 )
      LKNDOF = MAX( 1, NG )
      LELING = MAX( 1, NGEL )
      LELVAR = MAX( 1, NVARS )
      LICNA  = MAX( 1, NNZA )
      LSTADH = MAX( 1, NEL1 )
      LNTVAR = MAX( 1, NEL1 )
      LCALCF = MAX( 1, NELNUM, NG )
      LCALCG = MAX( 1, NG )
      LA     = MAX( 1, NNZA )
      LB     = MAX( 1, NG )
      LU     = MAX( 1, NG )
      LESCAL = MAX( 1, NGEL )
      LGSCAL = MAX( 1, NG )
      LVSCAL = MAX( 1, N )
      LFT    = MAX( 1, NG )
      LGVALS = MAX( 1, NG )
      LINTRE = MAX( 1, NELNUM )
      LGXEQX = MAX( 1, NGNG )
      LGPVLU = MAX( 1, NGPVLU )
      LEPVLU = MAX( 1, NEPVLU )
C     LSTGP  = MAX( 1, NG1 )
C     LSTEP  = MAX( 1, NEL1 )
C     LTYPEG = MAX( 1, NG )
C     LTYPEE = MAX( 1, NELNUM )
C     LIVAR  = MAX( 1, N )
C     LBL    = MAX( 1, N )
C     LBU    = MAX( 1, N )
C     LX     = MAX( 1, N )
C     LXT    = MAX( 1, N )
C     LDGRAD = MAX( 1, N )
C     LQ     = MAX( 1, N )

C
C  Print out problem data. input the number of variables, groups,
C  elements and the identity of the objective function group.
C
      IF ( IALGOR .EQ. 2 ) THEN
         READ( INPUT, 1002 ) NSLACK, NOBJGR
      ELSE
         NSLACK = 0
      END IF
      IF ( DEBUG ) WRITE( IOUT, 1100 ) PNAME, N, NG, NELNUM
      CHA( CH + 1 ) = PNAME // '  '
C
C  Input the starting addresses of the elements in each group,
C  of the parameters used for each group and
C  of the nonzeros of the linear element in each group.
C
      READ( INPUT, 1010 ) ( IWK( ISTADG + I ), I = 1, NG1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTADG',
     *   ( IWK( ISTADG + I ), I = 1, NG1 )
      READ( INPUT, 1010 ) ( IWK( ISTGP + I ), I = 1, NG1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTGP ',
     *   ( IWK( ISTGP + I ), I = 1, NG1 )
      READ( INPUT, 1010 ) ( IWK( ISTADA + I ), I = 1, NG1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTADA',
     *   ( IWK( ISTADA + I ), I = 1, NG1 )
C
C  Input the starting addresses of the variables and parameters
C  in each element.
C
      READ( INPUT, 1010 ) ( IWK( ISTAEV + I ), I = 1, NEL1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTAEV',
     *   ( IWK( ISTAEV + I ), I = 1, NEL1 )
      READ( INPUT, 1010 ) ( IWK( ISTEP + I ), I = 1, NEL1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTEP ',
     *   ( IWK( ISTEP + I ), I = 1, NEL1 )
C
C  Input the group type of each group
C
      READ( INPUT, 1010 ) ( IWK( ITYPEG + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ITYPEG',
     *   ( IWK( ITYPEG + I ), I = 1, NG )
      IF ( IALGOR .GE. 2 ) THEN
         READ( INPUT, 1010 )( IWK( KNDOFC + I ), I = 1, NG )
         IF ( DEBUG ) WRITE( IOUT, 1110 ) 'KNDOFC',
     *      ( IWK( KNDOFC + I ), I = 1, NG )
         DO 10 I = 1, NG
            IF ( ABS( IWK( KNDOFC + I ) ) .GE. 2 ) THEN
               CLOSE( INPUT )
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2010 )
               STOP
            END IF
   10    CONTINUE
      END IF
C
C  Input the element type of each element
C
      READ( INPUT, 1010 ) ( IWK( ITYPEE + I ), I = 1, NELNUM )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ITYPEE',
     *   ( IWK( ITYPEE + I ), I = 1, NELNUM )
C
C  Input the number of internal variables for each element.
C
      READ( INPUT, 1010 ) ( IWK( INTVAR + I ), I = 1, NELNUM )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'INTVAR',
     *   ( IWK( INTVAR + I ), I = 1, NELNUM )
C
C  Input the identity of each individual element.
C
      READ( INPUT, 1010 ) ( IWK( IELING + I ), I = 1, NGEL )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'IELING',
     *   ( IWK( IELING + I ), I = 1, NGEL )
C
C  Input the variables in each group's elements.
C
      NVARS = IWK( ISTAEV + NEL1 ) - 1
      READ( INPUT, 1010 ) ( IWK( IELVAR + I ), I = 1, NVARS )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'IELVAR',
     *   ( IWK( IELVAR + I ), I = 1, NVARS )
C
C  Input the column addresses of the nonzeros in each linear element.
C
      READ( INPUT, 1010 ) ( IWK( ICNA + I ), I = 1, NNZA )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ICNA  ',
     *   ( IWK( ICNA + I ), I = 1, NNZA )
C
C  Input the values of the nonzeros in each linear element, the
C  constant term in each group, the lower and upper bounds on
C  the variables and the starting point for the minimization.
C
      READ( INPUT, 1020 ) ( WK( A + I ), I = 1, NNZA )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'A     ',
     *   ( WK( A + I ), I = 1, NNZA )
      READ( INPUT, 1020 ) ( WK( B + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'B     ',
     *   ( WK( B + I ), I = 1, NG )
      IF ( IALGOR .LE. 2 ) THEN
         READ( INPUT, 1020 ) ( BL( I ), I = 1, N )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BL    ',
     *      ( BL( I ), I = 1, N )
         READ( INPUT, 1020 ) ( BU( I ), I = 1, N )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BU    ',
     *      ( BU( I ), I = 1, N )
      ELSE
C
C  Use GVALS and FT as temporary storage for the constraint bounds.
C
         READ( INPUT, 1020 ) ( BL( I ), I = 1, N ),
     *                       ( WK( GVALS + I ), I = 1, NG )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BL    ',
     *      ( BL( I ), I = 1, N ), ( WK( GVALS + I ), I = 1, NG )
         READ( INPUT, 1020 ) ( BU( I ), I = 1, N ),
     *                       ( WK( FT + I ), I = 1, NG )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BU    ',
     *      ( BU( I ), I = 1, N ), ( WK( FT + I ), I = 1, NG )
      END IF
      READ( INPUT, 1020 ) ( X( I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'X     ',
     *   ( X( I ), I = 1, N )
      IF ( IALGOR .GE. 2 ) THEN
         READ( INPUT, 1020 )( WK( U + I ), I = 1, NG )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'U     ',
     *      ( WK( U + I ), I = 1, NG )
      END IF
C
C  Input the parameters in each group.
C
      READ( INPUT, 1020 ) ( WK( GPVALU + I ), I = 1, NGPVLU )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'GPVALU',
     *   ( WK( GPVALU + I ), I = 1, NGPVLU )
C
C  Input the parameters in each individual element.
C
      READ( INPUT, 1020 ) ( WK( EPVALU + I ), I = 1, NEPVLU )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'EPVALU',
     *   ( WK( EPVALU + I ), I = 1, NEPVLU )
C
C  Input the scale factors for the nonlinear elements.
C
      READ( INPUT, 1020 ) ( WK( ESCALE + I ), I = 1, NGEL )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'ESCALE',
     *   ( WK( ESCALE + I ), I = 1, NGEL )
C
C  Input the scale factors for the groups.
C
      READ( INPUT, 1020 ) ( WK( GSCALE + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'GSCALE',
     *   ( WK( GSCALE + I ), I = 1, NG )
C
C  Input the scale factors for the variables.
C
      READ( INPUT, 1020 ) ( WK( VSCALE + I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'VSCALE',
     *   ( WK( VSCALE + I ), I = 1, N )
C
C  Input the lower and upper bounds on the objective function.
C
      READ( INPUT, 1080 ) OBFBND( 1 ), OBFBND( 2 )
      IF ( DEBUG ) WRITE( IOUT, 1180 ) 'OBFBND',
     *    OBFBND( 1 ), OBFBND( 2 )
C
C  Input a logical array which says whether an element has internal
C  varaiables.
C
      READ( INPUT, 1030 ) ( LOGI( INTREP + I ), I = 1, NELNUM )
      IF ( DEBUG ) WRITE( IOUT, 1130 ) 'INTREP',
     *   ( LOGI( INTREP + I ), I = 1, NELNUM )
C
C  Input a logical array which says whether a group is trivial.
C
      READ( INPUT, 1030 ) ( LOGI( GXEQX + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1130 ) 'GXEQX ',
     *   ( LOGI( GXEQX + I ), I = 1, NG )
C
C  Input the names given to the groups and to the variables.
C
      READ( INPUT, 1040 ) ( CHA( GNAMES + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1140 ) 'GNAMES',
     *   ( CHA( GNAMES + I ), I = 1, NG )
      READ( INPUT, 1040 ) ( CHA( VNAMES + I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1140 ) 'VNAMES',
     *   ( CHA( VNAMES + I ), I = 1, N )
C
C  Dummy input for the names given to the element and group types.
C
      READ( INPUT, 1040 ) ( CHTEMP, I = 1, NELTYP )
      READ( INPUT, 1040 ) ( CHTEMP, I = 1, NGRTYP )
C
C  Input the type of each variable.
C
      READ( INPUT, 1010 ) ( IWK( ITYPEV + I ), I = 1, N )
      CLOSE( INPUT )

      NUMVAR = N

C
C  Partition the workspace arrays FUVALS, IWK and WK. Initialize
C  certain portions of IWK.
C
      FIRSTG = .TRUE.
      FDGRAD = .FALSE.
CS    CALL SINITW( N, NG, NELNUM, IWK(IELING+1), LELING, IWK(ISTADG+1),
CD    CALL DINITW( N, NG, NELNUM, IWK(IELING+1), LELING, IWK(ISTADG+1),
     *    LSTADG, IWK(IELVAR+1), LELVAR, IWK(ISTAEV+1), LSTAEV,
     *    IWK(INTVAR+1), LNTVAR, IWK(ISTADH+1), LSTADH,
     *    IWK(ICNA+1), LICNA, IWK(ISTADA+1), LSTADA,
     *    IWK(ITYPEE+1), LINTRE,
     *    LOGI(GXEQX+1), LGXEQX, LOGI(INTREP+1), LINTRE,
     *    LFUVAL, ALTRIV, .TRUE., FDGRAD, LFXI,LGXI,LHXI,LGGFX,
     *    LDX   , LGRJAC, LQGRAD, LBREAK, LP,     LXCP  , LX0   , 
     *    LGX0  , LDELTX, LBND,   LWKSTR, LSPTRS, LSELTS, LINDEX, 
     *    LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC, LNNONZ, LNONZ2, 
     *    LSYMMD, LSYMMH, LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND , 
     *    LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV, LNSTJC, LNIUSE, 
     *    LNFREC, LNNNON, LNNNO2, LNSYMD, LNSYMH, LNLGRP, LNVGRP,
     *    LNGCLJ, LNVLJR, LNQGRD, LNBRAK, LNP,    LNBND ,
     *    LNFXI,  LNGXI,  LNGUVL, LNHXI,  LNHUVL, LNGGFX,
     *    LNDX  , LNGRJC, LIWK2 , LWK2  , MAXSIN, NINVAR,
     *    NTYPE , NSETS , MAXSEL, LSTYPE, LSSWTR, LSSIWT,
     *    LSIWTR, LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR,
     *    LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE, RANGE ,
     *    IWK(IWRK+1),LIWORK,WK(WRK+1),LWORK,IPRINT,IOUT,INFORM )
      IF ( INFORM .NE. 0 ) STOP
C
C  Shift the starting addresses for the real workspace relative to WRK.
C
      LQGRAD = LQGRAD + WRK
      LBREAK = LBREAK + WRK
      LP     = LP     + WRK
      LXCP   = LXCP   + WRK
      LX0    = LX0    + WRK
      LGX0   = LGX0   + WRK
      LDELTX = LDELTX + WRK
      LBND   = LBND   + WRK
      LSWTRA = LSWTRA + WRK
      LWKSTR = LWKSTR + WRK
C
C  Shift the starting addresses for the integer workspace relative
C  to IWRK.
C
      LSPTRS = LSPTRS + IWRK
      LSELTS = LSELTS + IWRK
      LINDEX = LINDEX + IWRK
      LSWKSP = LSWKSP + IWRK
      LSTAGV = LSTAGV + IWRK
      LSTAJC = LSTAJC + IWRK
      LIUSED = LIUSED + IWRK
      LFREEC = LFREEC + IWRK
      LNNONZ = LNNONZ + IWRK
      LNONZ2 = LNONZ2 + IWRK
      LSYMMD = LSYMMD + IWRK
      LSYMMH = LSYMMH + IWRK
      LSLGRP = LSLGRP + IWRK
      LSVGRP = LSVGRP + IWRK
      LGCOLJ = LGCOLJ + IWRK
      LVALJR = LVALJR + IWRK
      LSTYPE = LSTYPE + IWRK
      LSSWTR = LSSWTR + IWRK
      LSSIWT = LSSIWT + IWRK
      LSIWTR = LSIWTR + IWRK
      LSISET = LSISET + IWRK
      LSSVSE = LSSVSE + IWRK
      LSEND  = LSEND  + IWRK
C
C  Initialize the performance counters and variables
C
      NC2OF  = 0
      NC2OG  = 0
      NC2OH  = 0
      NC2CF  = 0
      NC2CG  = 0
      NC2CH  = 0
      NHVPR  = 0
      PNC    = 0
      STTIME = CPUTIM( DUM )
      SUTIME = STTIME - SUTIME
      RETURN
C
C  Non-executable statements.
C
 1000 FORMAT( I2, A8 )
 1001 FORMAT( 10I8 )
 1002 FORMAT( 2I8 )
 1010 FORMAT( ( 10I8 ) )
 1020 FORMAT( ( 1P, 4D16.8 ) )
 1030 FORMAT( ( 72L1 ) )
 1040 FORMAT( ( 8A10 ) )
 1080 FORMAT( 1P, 2D16.8 )
 1100 FORMAT( A8, 3I8 )
 1110 FORMAT( 1X, A6, /, ( 1X, 10I8 ) )
 1120 FORMAT( 1X, A6, /, ( 1X, 1P, 4D16.8 ) )
 1130 FORMAT( 1X, A6, /, ( 1X, 72L1 ) )
 1140 FORMAT( 1X, A6, /, ( 1X, 8A10 ) )
 1180 FORMAT( 1X, A6, /, 1P, 2D16.6 )
 2000 FORMAT( /, ' ** Program USETUP: array length ', A6, ' too small.',
     *        /, ' -- Miminimization abandoned.',
     *        /, ' -- Increase the parameter ', A6, ' by at least ', I8,
     *           ' and restart.'  )
 2010 FORMAT( /, ' ** Program USETUP: the problem includes general',
     *           ' constraints. Execution terminating ' )
 2020 FORMAT( /, ' ** Program USETUP: the problem uses no variables.',
     *           ' Execution terminating ' )
 2030 FORMAT( /, ' ** Program USETUP: the problem is vacuous.',
     *           ' Execution terminating ' )
C
C  End of USETUP.
C
      END
