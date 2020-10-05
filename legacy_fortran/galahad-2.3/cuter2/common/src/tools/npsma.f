C     ( Last modified on 23 Dec 2000 at 22:01:38 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     File npsma.f
C
C     Driver for running NPSOL Version 4.06 on CUTEr problems.
C
C     May 1993. Peihuang Chen
C     modified September 1993. Ingrid Bongartz
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      PROGRAM          NPSMA
C
C  Set up parameters, variables and arrays required by constrained tools.
C
      INTEGER          NMAX  , MMAX
CTOY  PARAMETER      ( NMAX =  55, MMAX =  51 )
CMED  PARAMETER      ( NMAX = 105, MMAX = 101 )
CBIG  PARAMETER      ( NMAX = 2505, MMAX = 501 )
CCUS  PARAMETER      ( NMAX = 1005, MMAX = 301 )
      CHARACTER * 10   PNAME , VNAME( NMAX ), CNAME( MMAX )
CS    REAL             X     ( NMAX  ), BL    ( NMAX   ),
CD    DOUBLE PRECISION X     ( NMAX  ), BL    ( NMAX   ),
     *                 BU    ( NMAX  ), V     ( MMAX   ),
     *                 CL    ( MMAX  ), CU    ( MMAX   )
      LOGICAL          EQUATN( MMAX  ), LINEAR( MMAX   )
      INTEGER          INPUT , IOUT
      PARAMETER      ( INPUT = 55, IOUT = 6 )
C
C  Set up parameters, variables and arrays required by NPSOL.
C
*  =====================================================================
*  Set the declared array dimensions for NPSOL
*
*  LDA    = the declared leading dimension of  A.
*  LDCJ   = the declared leading dimension of  CJAC.
*  LDR    = the declared leading dimension of  R.
*  MAXBND = maximum no. of variables + linear & nonlinear constrnts.
*  LIWORK = the length of the integer work array.
*  LWORK  = the length of the double precision work array.
*
*  Lengths of working arrays
*  according to User's Guide for NPSOL (Version 4.0):
*     LIWORK >= 3*N + NCLIN + 2*NCNLN
*     LWORK  >= 2*N*N + N*NCLIN + 2*N*NCNLN + 20*N + 11*NCLIN + 21*NCNLN
*
*  =====================================================================
      INTEGER          LDA, LDCJ, LDR, LIWORK, LWORK, MAXBND
      PARAMETER      ( LDA = MMAX , LDCJ = MMAX, LDR = NMAX + MMAX,
     *                 LIWORK = 3*( NMAX + MMAX ),
     *                 LWORK  = 2*NMAX*NMAX+3*NMAX*MMAX+20*NMAX+32*MMAX,
     *                 MAXBND = NMAX + MMAX )
      INTEGER          N, NCLIN, NCNLN, INFORM, ITER
      INTEGER          IWORK( LIWORK ), ISTATE( MAXBND )
CS    REAL             F
CD    DOUBLE PRECISION F
CS    REAL             G( NMAX ), R( LDR, NMAX )
CD    DOUBLE PRECISION G( NMAX ), R( LDR, NMAX )
CS    REAL             A( LDA, NMAX ), C( MMAX ), CJAC( LDCJ, NMAX )
CD    DOUBLE PRECISION A( LDA, NMAX ), C( MMAX ), CJAC( LDCJ, NMAX )
CS    REAL             BLOWER( MAXBND ), BUPPER( MAXBND ),
CD    DOUBLE PRECISION BLOWER( MAXBND ), BUPPER( MAXBND ),
     *                 CLAMDA( MAXBND )
CS    REAL             WORK( LWORK )
CD    DOUBLE PRECISION WORK( LWORK )
      EXTERNAL         OBJFUN, CONFUN
C
C  Local variable declarations.
C
      LOGICAL          DEBUG
      INTEGER          IOPTNS, IPRINT, I, M, MM
      REAL             CPU( 2 ), CALLS( 7 )
      CHARACTER*10     CBGBND
C                            finite-difference gradients allowed
      LOGICAL          FDGRAD
      COMMON / FDG   / FDGRAD
C                            and dimension IPADNP and IPSVNP only once
      INTEGER          MXPARM
      PARAMETER       ( MXPARM = 30 )
      INTEGER          IDBGNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4,
     *                 LDBGNP, LFORMH, LVLDER, LVERFY, MSGNP , NLNF  ,
     *                 NLNJ  , NLNX  , NNCNLN, NSAVE , NLOAD , KSAVE ,
     *                 IPADNP, IPSVNP
      COMMON / NPPAR1/ IPSVNP( MXPARM ),
     *                 IDBGNP, ITMXNP, JVRFY1, JVRFY2, JVRFY3, JVRFY4,
     *                 LDBGNP, LFORMH, LVLDER, LVERFY, MSGNP , NLNF  ,
     *                 NLNJ  , NLNX  , NNCNLN, NSAVE , NLOAD , KSAVE ,
     *                 IPADNP(12)
C
C     DEBUG = .TRUE.
      DEBUG = .FALSE.
C
C  Open the relevant file.
C
      OPEN ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND INPUT
C
C  Set up the data structures necessary to hold the problem functions.
C
      CALL NPSSE ( INPUT , IOUT  , NMAX  , MMAX  , MAXBND,
     *             N     , NCLIN , NCNLN , X     , BL    , BU    ,
     *             EQUATN, LINEAR, V     , CL    , CU    ,
     *             BLOWER, BUPPER, CLAMDA, LDA   , A     ,
     *             LDCJ  , CJAC )
C
C  Get the problem name and write some debug messages.
C
      M = NCLIN + NCNLN
C    call to CNAMES.
C    This fix avoids zero-length array difficulties on the VAX.
      MM = MAX( M, 1 )
      CALL CNAMES(N, MM, PNAME, VNAME, CNAME)
      WRITE ( IOUT, 2080 )  PNAME, N, NCLIN, NCNLN
      IF ( DEBUG ) THEN
         WRITE( IOUT, 2030 ) ( I, X( I ), BLOWER( I ),
     *     BUPPER( I ), I = 1, N )
         IF ( NCLIN .GT. 0 ) WRITE( IOUT, 2060 ) ( I, CLAMDA( N + I ),
     *     BLOWER( N + I ), BUPPER( N + I ), EQUATN( NCNLN + I ),
     *     I = 1, NCLIN )
         IF ( NCNLN .GT. 0 ) WRITE( IOUT, 2070 ) ( I, CLAMDA( N + NCLIN
     *     + I ), BLOWER( N + NCLIN + I ), BUPPER( N + NCLIN + I ),
     *     EQUATN( I ), I = 1, NCNLN )
      END IF
C
      CBGBND = '1.0D+15'
C
C  IOPTNS = the unit number for reading the options file.
C  IPRINT = the unit number for writing the output file.
C  Open and then read the options file.
C
      IOPTNS = 4
      IPRINT = 9
      OPEN ( UNIT=IOPTNS, FILE= 'NPSOL.SPC', STATUS='UNKNOWN' )
      CALL NPFILE( IOPTNS, INFORM )
C    gradients are used if requested
      FDGRAD = LVLDER .EQ. 0
      IF ( INFORM .NE. 0 .AND. IOUT .GT. 0 ) THEN
         WRITE ( IOUT, 3000 ) INFORM
         IF ( INFORM .EQ. 1 ) THEN
            WRITE ( IOUT, 3001 )
         ELSE IF ( INFORM .EQ. 2 ) THEN
            WRITE ( IOUT, 3002 )
         ELSE IF ( INFORM .EQ. 3 ) THEN
            WRITE ( IOUT, 3003 )
         ELSE IF ( INFORM .EQ. 4 ) THEN
            WRITE ( IOUT, 3004 )
         END IF
         STOP
      END IF
      CALL NPOPTN( 'Infinite Bound size ='//CBGBND )
C
C  Solve the problem.
C
      CALL NPSOL ( N, NCLIN, NCNLN, LDA, LDCJ, LDR,
     *             A, BLOWER, BUPPER,
     *             CONFUN, OBJFUN,
     *             INFORM, ITER, ISTATE,
     *             C, CJAC, CLAMDA, F, G, R, X,
     *             IWORK, LIWORK, WORK, LWORK )
      CALL CREPRT( CALLS, CPU )
C
C  Print messages about abnormal NPSOL inform codes.
C
      IF ( INFORM .GT. 0 .AND. IOUT .GT. 0 ) THEN
         WRITE ( IOUT, 4000 ) INFORM
         IF ( INFORM .EQ. 1 ) THEN
            WRITE ( IOUT, 4001 )
         ELSE IF ( INFORM .EQ. 2 ) THEN
            WRITE ( IOUT, 4002 )
         ELSE IF ( INFORM .EQ. 3 ) THEN
            WRITE ( IOUT, 4003 )
         ELSE IF ( INFORM .EQ. 4 ) THEN
            WRITE ( IOUT, 4004 )
         ELSE IF ( INFORM .EQ. 6 ) THEN
            WRITE ( IOUT, 4006 )
         ELSE IF ( INFORM .EQ. 7 ) THEN
            WRITE ( IOUT, 4007 )
         ELSE IF ( INFORM .EQ. 9 ) THEN
            WRITE ( IOUT, 4009 )
         END IF
      END IF
C
C  Output final objective function value and timing information.
C
      IF ( IPRINT .GT. 0 )
     *  WRITE ( IPRINT, 2000 ) PNAME, N, M, CALLS(1), 
     *                         CALLS(2), CALLS(5), CALLS(6), INFORM,
     *                         F, CPU(1), CPU(2)
      IF ( IOUT   .GT. 0 )
     *  WRITE ( IOUT  , 2000 ) PNAME, N, M, CALLS(1), 
     *                         CALLS(2), CALLS(5), CALLS(6), INFORM,
     *                         F, CPU(1), CPU(2)
C
      STOP
C
C  Non-executable statements.
C
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  NPSOL',    /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # constraints           =      ', I10 /
     *    ,' # objective functions   =        ', F8.2 /
     *    ,' # objective gradients   =        ', F8.2 / 
     *    ,' # constraints functions =        ', F8.2 /
     *    ,' # constraints gradients =        ', F8.2 /
     *    ,' Exit code               =      ', I10 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *    ,' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 2030 FORMAT( /, ' The starting point:',
     *        /, '     I      X        BLOWER      BUPPER',
     *        /, ( I6, 1P, 3D12.4 ) )
 2060 FORMAT( /, ' The linear constraints:',
     *        /, '     I  MULTIPLIER   BLOWER      BUPPER    EQUALITY?',
     *        /, ( I6, 1P, 3D12.4, 5X, L1 ) )
 2070 FORMAT( /, ' The nonlinear constraints:',
     *        /, '     I  MULTIPLIER   BLOWER      BUPPER    EQUALITY?',
     *        /, ( I6, 1P, 3D12.4, 5X, L1 ) )
 2080 FORMAT( /, ' Run NPSOL on Problem ', A10,
     *        /, ' N = ', I5, ' NCLIN = ', I5, ' NCNLN = ', I5 )
 3000 FORMAT( /, ' NPFILE terminated with  INFORM =', I3 )
 3001 FORMAT(    ' IOPTNS .LT. 0 or IOPTNS .GT. 99 ' )
 3002 FORMAT(    ' BEGIN was found, but end-of-file occurred before',
     *           ' END was found.' )
 3003 FORMAT(    ' End-of-file occurred before BEGIN or ENDRUN were',
     *           ' found.' )
 3004 FORMAT(    ' ENDRUN was found before BEGIN.' )
 4000 FORMAT( /, ' NPSOL  terminated with  INFORM =', I3 )
 4001 FORMAT(    ' Final iterate satisfies first-order Kuhn-Tucker',
     *           ' conditions',
     *        /, ' to accuracy requested, but iterates have not yet',
     *           ' converged.',
     *        /, ' No improvement could be made in merit function.' )
 4002 FORMAT(    ' No feasible point found for linear constraints and',
     *        /, ' bounds.  The problem is infeasible.' )
 4003 FORMAT(    ' No feasible point found for nonlinear constraints.',
     *        /, ' The problem may be infeasible.' )
 4004 FORMAT(    ' Maximum number of iterations reached.' )
 4006 FORMAT(    ' Final iterate does not satisfy Kuhn-Tucker',
     *           ' conditions',
     *        /, ' and no improved point could be found.')
 4007 FORMAT(    ' The provided derivatives of the objective function',
     *        /, ' or nonlinear constraints appear to be incorrect.' )
 4009 FORMAT(    ' An input parameter is invalid.' )
*     end of driver for NPSOL
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE NPSSE ( INPUT , IOUT  , NMAX  , MMAX  , MAXBND,
     *                   N     , NCLIN , NCNLN , X     , BL    , BU    ,
     *                   EQUATN, LINEAR, V     , CL    , CU    ,
     *                   BLOWER, BUPPER, CLAMDA, LDA   , AA    ,
     *                   LDCJ  , CJAC   )
      INTEGER            INPUT , IOUT  , NMAX  , MMAX  , MAXBND,
     *                   N     , NCLIN , NCNLN , LDA   , LDCJ
CS    REAL               X ( NMAX   ), BL( NMAX   ), BU( NMAX   )
CD    DOUBLE PRECISION   X ( NMAX   ), BL( NMAX   ), BU( NMAX   )
CS    REAL               V ( MMAX   ), CL( MMAX   ), CU( MMAX   )
CD    DOUBLE PRECISION   V ( MMAX   ), CL( MMAX   ), CU( MMAX   )
CS    REAL               BLOWER( MAXBND ), BUPPER( MAXBND ),
CD    DOUBLE PRECISION   BLOWER( MAXBND ), BUPPER( MAXBND ),
     *                   CLAMDA( MAXBND )
CS    REAL               AA( LDA , NMAX ), CJAC( LDCJ, NMAX )
CD    DOUBLE PRECISION   AA( LDA , NMAX ), CJAC( LDCJ, NMAX )
      LOGICAL            EQUATN( MMAX   ), LINEAR( MMAX   )
C
C  Set up the input data for NPSOL.
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
      INTEGER                IWK( LIWK    )
CS    REAL                   WK ( LWK     )
CD    DOUBLE PRECISION       WK ( LWK     )
      LOGICAL              LOGI ( LLOGIC  )
      CHARACTER * 10       CHA  ( LCHARA  )
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
      INTEGER            NNOV  , NNJV
      COMMON / NNVARS /  NNOV  , NNJV
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / NNVARS /
C
C  Local variable declarations.
C
      INTEGER            I , IB, IC, II, IG, J , JG, M , MEND
      LOGICAL            EFIRST, LFIRST, NVFRST, LTEMP
CS    REAL               ATEMP
CD    DOUBLE PRECISION   ATEMP
CS    REAL               ZERO
CD    DOUBLE PRECISION   ZERO
CS    PARAMETER        ( ZERO = 0.0E+0 )
CD    PARAMETER        ( ZERO = 0.0D+0 )
C
C  Input problem data using CSETUP.
C
      EFIRST = .FALSE.
      LFIRST = .FALSE.
      NVFRST = .FALSE.
      CALL CSETUP ( INPUT , IOUT  , N , M , X , BL , BU , NMAX,
     *              EQUATN, LINEAR, V , CL , CU , MMAX,
     *              EFIRST, LFIRST, NVFRST )
      CLOSE( INPUT )
C
C  Determine the number of linear and nonlinear constraints.
C
      NCLIN = 0
      NCNLN = 0
      DO 100 IG = 1, NG
         I = IWK( KNDOFC + IG )
         IF ( I .GT. 0 ) THEN
            IWK( LSEND + I ) = IG
            IF ( LINEAR( I ) ) THEN
               NCLIN = NCLIN + 1
            ELSE
               NCNLN = NCNLN + 1
            END IF
         END IF
  100 CONTINUE
      IF ( NCLIN .EQ. 0 .OR. NCNLN .EQ. 0 ) GO TO 130
C
C  Reorder the constraints so that the nonlinear constraints occur before
C  the linear ones.  The constraints are ordered in this way so that CCFG
C  need evaluate the Jacobian for only the first NCNLN constraints.
C
      MEND = M
C
C  Run forward through the constraints until a linear constraint
C  is encountered.
C
      DO 120 I = 1, M
         IF ( I .GT. MEND ) GO TO 130
         IG = IWK( LSEND + I )
         IF ( LINEAR( I ) ) THEN
C
C  Constraint I is linear. Now, run backwards through the
C  constraints until a nonlinear one is encountered.
C
            DO 110 J = MEND, I, - 1
               JG    = IWK( LSEND + J )
               IF ( .NOT. LINEAR( J ) ) THEN
                  MEND = J - 1
C
C  Interchange the data for constraints I and J.
C
                  IWK ( LSEND  +  I ) = JG
                  IWK ( LSEND  +  J ) = IG
                  IWK ( KNDOFC + IG ) = J
                  IWK ( KNDOFC + JG ) = I
                  LTEMP               = LINEAR( I )
                  LINEAR( I )         = LINEAR( J )
                  LINEAR( J )         = LTEMP
                  LTEMP               = EQUATN( I )
                  EQUATN( I )         = EQUATN( J )
                  EQUATN( J )         = LTEMP
                  ATEMP               = V     ( I )
                  V     ( I )         = V     ( J )
                  V     ( J )         = ATEMP
                  ATEMP               = CL    ( I )
                  CL    ( I )         = CL    ( J )
                  CL    ( J )         = ATEMP
                  ATEMP               = CU    ( I )
                  CU    ( I )         = CU    ( J )
                  CU    ( J )         = ATEMP
                  ATEMP               = BL    ( N + I )
                  BL    ( N + I )     = BL    ( N + J )
                  BL    ( N + J )     = ATEMP
                  ATEMP               = BU    ( N + I )
                  BU    ( N + I )     = BU    ( N + J )
                  BU    ( N + J )     = ATEMP
                 GO TO 120
               END IF
  110       CONTINUE
            GO TO 130
         END IF
  120 CONTINUE
  130 CONTINUE
C
C  Set up the lower bound vector BLOWER and upper bound vector BUPPER
C  in the order required by NPSOL.
C  For i=1 to n, set BLOWER (BUPPER) to the lower (upper) bound
C  on the variables.  (CSETUP put these bounds in BL and BU.)
C  For i=n+1 to n+nclin, set BLOWER (BUPPER) to the lower (upper)
C  bounds on the linear constraints.
C  For i=n+nclin+1 to n+nclin+ncnln, set BLOWER (BUPPER) to the lower
C  (upper) bounds on the nonlinear constraints.
C  At the same time, copy the multiplier estimates from V to CLAMDA.
C  CLAMDA has the same ordering as BLOWER and BUPPER.
C
      DO 150 I = 1, N
         BLOWER( I ) = BL( I )
         BUPPER( I ) = BU( I )
         CLAMDA( I ) = ZERO
  150 CONTINUE
      DO 160 I = 1, NCLIN
         IB = N + I
         IC = NCNLN + I
         BLOWER( IB ) = CL( IC )
         BUPPER( IB ) = CU( IC )
         CLAMDA( IB ) = V ( IC )
  160 CONTINUE
      DO 170 I = 1, NCNLN
         IB = N + NCLIN + I
         IC = I
         BLOWER( IB ) = CL( IC )
         BUPPER( IB ) = CU( IC )
         CLAMDA( IB ) = V ( IC )
  170 CONTINUE
C
C  Initialize AA and CJAC to zero.
C
      DO 230 J = 1, N
         DO 210 I = 1, NCLIN
            AA( I, J ) = ZERO
  210    CONTINUE
         DO 220 I = 1, NCNLN
            CJAC( I, J ) = ZERO
  220    CONTINUE
  230 CONTINUE
C
C  Set up the matrix AA, which contains the coefficients of the linear
C  constraints.  Also incorporate nonzero RHS constants of linear
C  constraints into the lower and upper bounds.
C
      DO 250 I = 1, NCLIN
         II = NCNLN + I
         IG = IWK( LSEND + II )
         DO 240 IC = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
            J = IWK( ICNA + IC )
            AA( I, J ) = WK( A + IC ) * WK( GSCALE + IG )
  240    CONTINUE
         IF ( WK( B + IG ) .NE. ZERO ) THEN
            IB = N + I
            IC = B + IG
            BLOWER( IB ) = BLOWER( IB ) + WK( IC )
     *                           * WK( GSCALE + IG )
            BUPPER( IB ) = BUPPER( IB ) + WK( IC )
     *                           * WK( GSCALE + IG  )
            WK    ( IC ) = ZERO
         END IF
  250 CONTINUE
C
C  Incorporate nonzero RHS constants of nonlinear constraints into
C  the lower and upper bounds.
C
      DO 260 I = 1, NCNLN
         IG = IWK( LSEND + I )
         IF ( WK( B + IG ) .NE. ZERO ) THEN
            IB = N + NCLIN + I
            IC = B + IG
            BLOWER( IB ) = BLOWER( IB ) + WK( IC )
     *                           * WK( GSCALE + IG  )
            BUPPER( IB ) = BUPPER( IB ) + WK( IC )
     *                           * WK( GSCALE + IG )
            WK    ( IC ) = ZERO
         END IF
  260 CONTINUE
      RETURN
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  THIS VERSION: 03/06/1994 AT 11:38:26 AM.
      SUBROUTINE OBJFUN( MODE, N, X, F, G, NSTATE )
      INTEGER           MODE, N, NSTATE
CS    REAL              F
CD    DOUBLE PRECISION  F
CS    REAL              X( N ), G( N )
CD    DOUBLE PRECISION  X( N ), G( N )
C
C  Local variables
C
      INTEGER           J
      LOGICAL           GRAD, FDGRAD
      COMMON / FDG /    FDGRAD
C
      IF ( MODE .EQ. 0 ) THEN
         GRAD = .FALSE.
      ELSE
         GRAD = .TRUE.
      END IF
      CALL COFG ( N, X, F, G, GRAD )
C    gradients are used if requested by setting the gradient to
C    appropriate NPSOL values.
      IF ( GRAD .AND. FDGRAD ) THEN
         DO 20 J = 1, N
            G( J ) = -11111.0D+0
   20    CONTINUE
      END IF
      RETURN
C     END OF OBJFUN
      END
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  THIS VERSION: 03/06/1994 AT 11:38:26 AM.
      SUBROUTINE CONFUN( MODE, NCNLN, N, LDCJ,
     *                   NEEDC, X, C, CJAC, NSTATE )
      INTEGER            MODE, NCNLN, N, LDCJ, NSTATE
      INTEGER            NEEDC(*)
CS    REAL               X( N ), C( LDCJ ), CJAC( LDCJ, N )
CD    DOUBLE PRECISION   X( N ), C( LDCJ ), CJAC( LDCJ, N )
C
C  Local variables
C
      INTEGER           I, J
      LOGICAL           GRAD, FDGRAD
      COMMON / FDG /    FDGRAD
C
      IF ( MODE .EQ. 0 ) THEN
         GRAD = .FALSE.
      ELSE
         GRAD = .TRUE.
      END IF
      CALL CCFG ( N, NCNLN, X, LDCJ, C, .FALSE., LDCJ, N, CJAC, GRAD )
C    gradients are used if requested by setting the Jacobian to
C    appropriate NPSOL values.
      IF ( GRAD .AND. FDGRAD ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, NCNLN
               CJAC( I, J ) = -11111.0D+0
   10       CONTINUE
   20    CONTINUE
      END IF
      RETURN
C     END OF CONFUN
      END

