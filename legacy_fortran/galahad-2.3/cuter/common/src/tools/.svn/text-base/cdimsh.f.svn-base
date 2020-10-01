C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CDIMSH( NNZH  )
      INTEGER            NNZH
C
C  Compute the space required to store the Hessian matrix of the 
C  Lagrangian function of a problem initially written in 
C  Standard Input Format (SIF).
C
C  The upper triangle of the Hessian is stored in coordinate form,
C  i.e., the entry H(i) has row index IRNH(i)
C  for i = 1, ...., NNZH.
C
C  Based on the minimization subroutine LANCELOT/SBMIN
C  by Conn, Gould and Toint.
C
C  Nick Gould, for CGT productions,
C  August 1999.
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
      INTEGER            NXTRW1, NXTRW2, LNXTRW, LIRNH, IRNH
      INTEGER            I, II,  IG, J,  JJ, K,  L, IEL, IELL  , INEXT 
      INTEGER            NEWPT,  IG1,    LISTVS, LISTVE, ISTART
C
C  Define the integer work space needed for ASMBI.
C  Ensure that there is sufficient space
C
      NNZH   = 0
      NEWPT  = NUMVAR + 1
      LIRNH  = ( LIWK2 - 2 * NUMVAR ) / 3
      IRNH   = LSEND
      LNXTRW = ( LIWK2 - LIRNH ) / 2
      NXTRW1 = IRNH + LIRNH
      NXTRW2 = NXTRW1 + LNXTRW
      IF ( NEWPT .GT. LNXTRW .OR. LIRNH .LE. 0 ) GO TO 900
C
C  NXTROW( 1, . ) gives the link list. The list for column J starts
C                 in NXTROW( 1, J ) and ends when NXTROW( 1, K ) = - 1.
C  NXTROW( 2, . ) gives the position in H of the current link.
C
C  Initialize the link list which points to the row numbers which
C  are used in the columns of the assembled Hessian
C
      DO 20 I = 1, NUMVAR
         IWK( NXTRW1 + I ) = - 1
   20 CONTINUE
C
C  -------------------------------------------------------
C  Form the rank-one second order term for the IG-th group
C  -------------------------------------------------------
C
      DO 200 IG = 1, NG
         IF ( LOGI( GXEQX + IG ) ) GO TO 200
         IG1    = IG + 1
         LISTVS = IWK( LSTAGV + IG )
         LISTVE = IWK( LSTAGV + IG1 ) - 1
C
C  Form the J-th column of the rank-one matrix
C
         DO 190 L = LISTVS, LISTVE
            J     = IWK( LSVGRP + L )
            IF ( J .EQ. 0 ) GO TO 190
C
C  Find the entry in row I of this column
C
            DO 180 K = LISTVS, LISTVE
               I     = IWK( LSVGRP + K )
               IF ( I .EQ. 0 .OR. I .GT. J ) GO TO 180
C
C  Obtain the appropriate storage location in H for the new entry
C
               ISTART = J
  150          CONTINUE
               INEXT = IWK( NXTRW1 + ISTART )
               IF ( INEXT .EQ. - 1 ) THEN
                  IF ( NEWPT .GT. LNXTRW ) GO TO 900
C
C  The (I,J)-th location is empty. Place the new entry in this location
C  and add another link to the list
C
                  NNZH = NNZH + 1
                  IF ( NNZH .GT. LIRNH ) GO TO 900
                  IWK( IRNH + NNZH )     = I
                  IWK( NXTRW1 + ISTART ) = NEWPT
                  IWK( NXTRW2 + ISTART ) = NNZH
                  IWK( NXTRW1 + NEWPT )  = - 1
                  NEWPT                  = NEWPT + 1
               ELSE
C
C  Continue searching the linked list for an entry in row I, column J
C
                  IF ( IWK( IRNH + IWK( NXTRW2 + ISTART ) ).NE.I ) THEN
                     ISTART = INEXT
                     GO TO 150
                  END IF
               END IF
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
C
C  --------------------------------------------------------
C  Add on the low rank first order terms for the I-th group
C  --------------------------------------------------------
C
      DO 300 IG = 1, NG
         IG1    = IG + 1
C
C  See if the group has any nonlinear elements
C
         DO 290 IELL = IWK( ISTADG + IG ), IWK( ISTADG + IG1 ) - 1
            IEL      = IWK( IELING + IELL )
            LISTVS   = IWK( ISTAEV + IEL )
            LISTVE   = IWK( ISTAEV + IEL + 1 ) - 1
            DO 250 L = LISTVS, LISTVE
               J     = IWK( IELVAR + L )
               IF ( J .NE. 0 ) THEN
C
C  The IEL-th element has an internal representation.
C  Compute the J-th column of the element Hessian matrix
C
C  Find the entry in row I of this column
C
                  DO 240 K = LISTVS, L
                     I     = IWK( IELVAR + K )
                     IF ( I .NE. 0 ) THEN
C
C  Only the upper triangle of the matrix is stored
C
                        IF ( I .LE. J ) THEN
                           II = I
                           JJ = J
                        ELSE
                           II = J
                           JJ = I
                        END IF
C
C  Obtain the appropriate storage location in H for the new entry
C
                        ISTART = JJ
  230                   CONTINUE
                        INEXT = IWK( NXTRW1 + ISTART )
                        IF ( INEXT .EQ. - 1 ) THEN
                           IF ( NEWPT .GT. LNXTRW ) GO TO 900
C
C  The (I,J)-th location is empty. Place the new entry in this location
C  and add another link to the list
C
                           NNZH = NNZH + 1
                           IF ( NNZH .GT. LIRNH ) GO TO 900
                           IWK( IRNH + NNZH )     = II
                           IWK( NXTRW1 + ISTART ) = NEWPT
                           IWK( NXTRW2 + ISTART ) = NNZH
                           IWK( NXTRW1 + NEWPT )  = - 1
                           NEWPT                  = NEWPT + 1
                        ELSE
C
C  Continue searching the linked list for an entry in row I, column J
C
                           IF ( IWK( IRNH + IWK( NXTRW2 + ISTART ) )
     *                       .EQ. II ) THEN
                           ELSE
                              ISTART = INEXT
                              GO TO 230
                           END IF
                        END IF
                     END IF
  240             CONTINUE
               END IF
  250       CONTINUE
  290    CONTINUE
  300 CONTINUE
      RETURN
C
C  Unsuccessful returns.
C
  900 CONTINUE
      WRITE( IOUT, 2000 )
      STOP
C
C Non-executable statements.
C
 2000 FORMAT( ' ** SUBROUTINE CDIMSH: Increase the size of IWK ' )
C
C  end of CDIMSH.
C
      END
