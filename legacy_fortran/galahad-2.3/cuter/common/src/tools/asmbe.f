C     ( Last modified on 23 Dec 2000 at 22:01:38 )
CS    SUBROUTINE SASMBE( N , NG, MAXSEL, ISTADH, LSTADH, ICNA  , LICNA ,
CD    SUBROUTINE DASMBE( N , NG, MAXSEL, ISTADH, LSTADH, ICNA  , LICNA ,
     *                   ISTADA, LSTADA, INTVAR, LNTVAR, IELVAR, LELVAR,
     *                   IELING, LELING, ISTADG, LSTADG, ISTAEV, LSTAEV,
     *                   ISTAGV, LNSTGV, ISVGRP, LNVGRP, IWK   , LIWK  , 
     *                   A , LA, GUVALS, LNGUVL, HUVALS, LNHUVL, GVALS2, 
     *                   GVALS3, GSCALE, ESCALE, LESCAL, WK    , LWK   , 
     *                   GXEQX , LGXEQX, INTREP, LINTRE, ITYPEE, LITYPE,
     *                   RANGE , NE    , IRNHI , LIRNHI, IPRNHI, HI    , 
     *                   LHI   , IPRHI , BYROWS, IPRINT, IOUT  , INFORM)
C
C  ********************************************************************
C
C  Assemble the second derivative matrix of a groups partially
C  separable function into finite-element format
C
C  Nick Gould, November 25th 1994, for CGT Productions.
C
C  ********************************************************************
C
      INTEGER          N , NG, MAXSEL, LIRNHI, LHI   , NE
      INTEGER          LA    , LITYPE, IOUT  , INFORM, IPRINT
      INTEGER          LSTADH, LICNA , LSTADA, LNTVAR, LELVAR, LELING
      INTEGER          LSTADG, LSTAEV, LNSTGV, LNVGRP, LIWK
      INTEGER          LNGUVL, LNHUVL, LESCAL, LWK   , LGXEQX, LINTRE
      LOGICAL          BYROWS
      INTEGER          ISTADH( LSTADH ), ICNA  ( LICNA  )
      INTEGER          ISTADA( LSTADA ), INTVAR( LNTVAR )
      INTEGER          IELVAR( LELVAR ), IELING( LELING )
      INTEGER          ISTADG( LSTADG ), ISTAEV( LSTAEV )
      INTEGER          ISTAGV( LNSTGV ), ISVGRP( LNVGRP )
      INTEGER          IRNHI ( LIRNHI ), IWK   ( LIWK   )
      INTEGER          IPRNHI( NG + 1 ), IPRHI ( NG + 1 )
      INTEGER          ITYPEE( LITYPE )
      LOGICAL          GXEQX( LGXEQX ),  INTREP( LINTRE )
CS    REAL             A(     LA ),      GVALS2( NG ),     GVALS3( NG ),
CD    DOUBLE PRECISION A(     LA ),      GVALS2( NG ),     GVALS3( NG ),
     *                 GUVALS( LNGUVL ), HUVALS( LNHUVL ), GSCALE( NG ),
     *                 ESCALE( LESCAL ), HI    ( LHI    ), WK(     LWK )
      EXTERNAL         RANGE 
C
C  Local variables.
C
      INTEGER          I , II, J, IHI, JJ, K , L     , NN, IG, NNN
      INTEGER          NIN,    IELL,   IEL   , IHNEXT, NVARG, NSIZEE
      INTEGER          NVAREL, IG1,    LISTVS, LISTVE, IELH
      INTEGER          IJHESS, IROW,   JCOL,   JCOLST
CS    REAL             ONE,    ZERO,   WKI,    HESNEW, GDASH,  G2DASH,
CD    DOUBLE PRECISION ONE,    ZERO,   WKI,    HESNEW, GDASH,  G2DASH,
     *                 SCALEE
CS    EXTERNAL         SSETVL, SSETVI
CD    EXTERNAL         DSETVL, DSETVI
C
C  Set constant real parameters.
C
CS    PARAMETER ( ZERO   = 0.0E+0, ONE    = 1.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0, ONE    = 1.0D+0 )
C
C  ------------------------------------------------------
C  Form the rank-one second order term for the i-th group.
C  ------------------------------------------------------
C
      NE          = 0
      IPRNHI( 1 ) = 1
      IPRHI ( 1 ) = 1
      DO 200 IG = 1, NG
         IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2070 ) IG
         IG1    = IG + 1
         IF ( GXEQX( IG ) ) THEN
            G2DASH = ZERO
         ELSE
            G2DASH = GSCALE( IG ) * GVALS3( IG )
            IF ( IPRINT .GE. 100 )WRITE(6,*) ' GVALS3(IG) ', GVALS3(IG)
         END IF
C
C  Ignore linear groups
C
         IF ( ISTADG( IG ) .GE. ISTADG( IG1 ) .AND. GXEQX( IG ) ) 
     *      GO TO 200
C
C  The group is nonlinear
C
         NE     = NE + 1
         LISTVS = ISTAGV( IG )
         LISTVE = ISTAGV( IG1 ) - 1
         NVARG  = LISTVE - LISTVS + 1
C  
C  Set the starting addresses for the integer and real arrays for
C  group IG + 1
C
         IPRNHI( NE + 1 ) = IPRNHI( NE ) + NVARG
         IF ( IPRNHI( NE + 1 ) .GT. LIRNHI ) THEN
            WRITE( IOUT, 2030 ) LIRNHI
            GO TO 610
         END IF
         NSIZEE = ( NVARG * ( NVARG + 1 ) ) / 2
         IPRHI( NE + 1 ) = IPRHI( NE ) + NSIZEE
         IF ( IPRHI( NE + 1 ) .GE. LHI ) THEN
            WRITE( IOUT, 2040 ) LHI, IPRHI( NE + 1 )
            GO TO 610
         END IF
C
C  Record the row indices involved in super-element NE
C
         IF ( GXEQX( IG ) .OR. G2DASH .EQ. ZERO ) THEN 
CS          CALL SSETVL( NSIZEE, HI( IPRHI( NE ) ), 1, ZERO )
CD          CALL DSETVL( NSIZEE, HI( IPRHI( NE ) ), 1, ZERO )
            GO TO 200
         END IF
         K = IPRNHI( NE )
         DO 110 L = LISTVS, LISTVE
            IRNHI( K ) = ISVGRP( L )
            K = K + 1
  110    CONTINUE   
C
C  Form the gradient of the ig-th group.
C
CS       CALL SSETVI( NVARG, WK, ISVGRP( LISTVS ), ZERO )
CD       CALL DSETVI( NVARG, WK, ISVGRP( LISTVS ), ZERO )
C
C  Consider any nonlinear elements for the group.
C
         DO 160 IELL = ISTADG( IG ), ISTADG( IG1 ) - 1
            IEL      = IELING( IELL )
            K        = INTVAR( IEL )
            L        = ISTAEV( IEL )
            NVAREL   = ISTAEV( IEL + 1 ) - L
            SCALEE   = ESCALE( IELL )
            IF ( INTREP( IEL ) ) THEN
C
C  The IEL-th element has an internal representation.
C
               NIN = INTVAR( IEL + 1 ) - K
               CALL RANGE ( IEL, .TRUE., GUVALS( K ),
     *                      WK( N + 1 ), NVAREL, NIN,
     *                      ITYPEE( IEL ), NIN, NVAREL )
               DO 140 I   = 1, NVAREL
                  J       = IELVAR( L )
                  WK( J ) = WK( J ) + SCALEE * WK( N + I )
                  L       = L + 1
  140          CONTINUE
            ELSE
C
C  The IEL-th element has no internal representation.
C
               DO 150 I   = 1, NVAREL
                  J       = IELVAR( L )
                  WK( J ) = WK( J ) + SCALEE * GUVALS( K )
                  K       = K + 1
                  L       = L + 1
  150          CONTINUE
            END IF
  160    CONTINUE
C
C  Include the contribution from the linear element.
C
         DO 170 K   = ISTADA( IG ), ISTADA( IG1 ) - 1
            J       = ICNA( K )
            WK( J ) = WK( J ) + A( K )
  170    CONTINUE
C
C  The gradient is complete. Form the J-th column of the rank-one matrix
C
         DO 190 L = LISTVS, LISTVE
            JJ    = ISVGRP( L )
            J = L - LISTVS + 1
C
C  Find the entry in row I of this column.
C
            DO 180 K = LISTVS, L
               II    = ISVGRP( K )
               I = K - LISTVS + 1
               IF ( BYROWS ) THEN
                  IHI = IPRHI( NE ) - 1 + NVARG * ( I - 1 ) - 
     *                  ( ( I - 1 ) * I ) / 2 + J
               ELSE
                  IHI = IPRHI( NE ) - 1 + I + ( J * ( J - 1 ) ) / 2
               END IF
               HI( IHI ) = WK( II ) * WK( JJ ) * G2DASH
               IF ( IPRINT .GE. 100 ) 
     *            WRITE( IOUT, 2090 ) II, JJ, HI( IHI )
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
C
C  Reset the workspace array to zero.
C
CS    CALL SSETVL( MAXSEL, WK, 1, ZERO )
CD    CALL DSETVL( MAXSEL, WK, 1, ZERO )
C
C  --------------------------------------------------------
C  Add on the low rank first order terms for the I-th group.
C  --------------------------------------------------------
C
      NE = 0
      DO 300 IG = 1, NG
         IG1 = IG + 1
C
C  Once again, ignore linear groups
C
         IF ( ISTADG( IG ) .GE. ISTADG( IG1 ) .AND. GXEQX( IG ) ) 
     *      GO TO 300
C
C  The group is nonlinear
C
         NE = NE + 1
         IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2100 ) IG
         IF ( GXEQX( IG ) ) THEN
            GDASH = GSCALE( IG )
         ELSE
            GDASH = GSCALE( IG ) * GVALS2( IG )
            IF ( IPRINT .GE. 100 )WRITE(6,*) ' GVALS2(IG) ', GVALS2(IG)
         END IF
         IF ( GDASH .EQ. ZERO ) GO TO 300
C
C  Map the problem variables to the elemental variables.
C
         NVARG = IPRNHI( NE + 1 ) - IPRNHI( NE )
         DO 210 I = IPRNHI( NE ), IPRNHI( NE + 1 ) - 1
            IWK( IRNHI( I ) ) = I + 1 - IPRNHI( NE )
  210    CONTINUE   
C
C  See if the group has any nonlinear elements.
C
         DO 290 IELL = ISTADG( IG ), ISTADG( IG1 ) - 1
            IEL      = IELING( IELL )
            LISTVS   = ISTAEV( IEL )
            LISTVE   = ISTAEV( IEL + 1 ) - 1
            NVAREL   = LISTVE - LISTVS + 1
            IELH     = ISTADH( IEL )
            IHNEXT   = IELH
            SCALEE   = ESCALE( IELL )
            DO 250 L = LISTVS, LISTVE
               J     = IWK( IELVAR( L ) )
C
C  The IEL-th element has an internal representation.
C  Compute the J-th column of the element Hessian matrix.
C
               IF ( INTREP( IEL ) ) THEN
C
C  Compute the J-th column of the Hessian.
C
                  WK( L - LISTVS + 1 ) = ONE
C
C  Find the internal variables.
C
                  NIN = INTVAR( IEL + 1 ) - INTVAR( IEL )
                  CALL RANGE ( IEL, .FALSE., WK( 1 ),
     *                     WK( MAXSEL + 1 ), NVAREL, NIN,
     *                     ITYPEE( IEL ), NVAREL, NIN )
C
C  Multiply the internal variables by the element Hessian.
C
                  NN = MAXSEL + NIN
CS                CALL SSETVL( NIN, WK( NN + 1 ), 1, ZERO )
CD                CALL DSETVL( NIN, WK( NN + 1 ), 1, ZERO )
C
C  Only the upper triangle of the element Hessian is stored.
C
                  JCOLST = IELH - 1
                  DO 230 JCOL = 1, NIN
                     IJHESS   = JCOLST
                     JCOLST   = JCOLST + JCOL
                     WKI      = WK( MAXSEL + JCOL ) * GDASH
                     DO 220 IROW = 1, NIN
                        IF ( IROW .LE. JCOL ) THEN
                           IJHESS = IJHESS + 1
                        ELSE
                           IJHESS = IJHESS + IROW - 1
                        END IF
                        WK( NN + IROW ) = WK( NN + IROW ) +
     *                                    WKI * HUVALS( IJHESS )
  220                CONTINUE
  230             CONTINUE
C
C  Scatter the product back onto the elemental variables.
C
                  NNN = NN + NIN
                  CALL RANGE ( IEL, .TRUE., WK( NN + 1 ),
     *                         WK( NNN + 1 ), NVAREL, NIN,
     *                         ITYPEE( IEL ), NIN, NVAREL )
                  WK( L - LISTVS + 1 ) = ZERO
C
C  Find the entry in row I of this column.
C
               END IF
               DO 240 K = LISTVS, L
                  I     = IWK( IELVAR( K ) )
C
C  Only the upper triangle of the matrix is stored.
C
                  IF ( I .LE. J ) THEN
                     II = I
                     JJ = J
                  ELSE
                     II = J
                     JJ = I
                  END IF
C
C  Obtain the appropriate storage location in H for the new entry.
C
                  IF ( INTREP( IEL ) ) THEN
                     HESNEW = SCALEE * WK( NNN + K - LISTVS + 1 )
                  ELSE
                     HESNEW = SCALEE * HUVALS( IHNEXT ) * GDASH
                  END IF
                  IF ( IPRINT .GE. 100 )
     *               WRITE( 6, 2080 ) II, JJ, IEL, HESNEW
                  IF ( BYROWS ) THEN
                     IHI = IPRHI( NE ) - 1 + NVARG * ( II - 1 ) - 
     *                     ( ( II - 1 ) * II ) / 2 + JJ
                  ELSE
                     IHI = IPRHI( NE ) - 1 + II + 
     *                     ( JJ * ( JJ - 1 ) ) / 2
                  END IF
                  HI( IHI ) = HI( IHI ) + HESNEW
                  IF ( K .NE. L .AND. II .EQ. JJ )
     *                 HI( IHI ) = HI( IHI ) + HESNEW
                  IHNEXT = IHNEXT + 1
  240          CONTINUE
  250       CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C  ----------------------------------------
C  For debugging, print the nonzero values.
C  ----------------------------------------
C
      IF ( IPRINT .GE. 10 ) THEN
         DO 410 IG = 1, NE
            WRITE( IOUT, 2000 ) IG
            WRITE( IOUT, 2010 )
     *        ( IRNHI( I ), I = IPRNHI( IG ), IPRNHI( IG + 1 ) - 1 )
            WRITE( IOUT, 2020 )
     *        ( HI( I ), I = IPRHI( IG ), IPRHI( IG + 1 ) - 1 )
  410    CONTINUE   
      END IF
      INFORM = 0
      RETURN
C
C  Unsuccessful returns.
C
  610 CONTINUE
      INFORM = 1
      RETURN
C
C  Non-executable statements
C
 2000 FORMAT( ' Super-element ', I10 )   
 2010 FORMAT( ' Super-element variables     ', 8I7, /, ( 11I7 ) )
 2020 FORMAT( ' Nonzeros   ', 1P, 6D12.4, /, ( 7D12.4 ) )
 2030 FORMAT( ' ** Array dimension LIRNHI = ',I8,' too small in ASMBE.')
 2040 FORMAT( ' ** Array dimension LHI = ', I8,' too small in ASMBE. 
     *     Increase to at least ', I8 )
 2070 FORMAT( ' Group ', I5, ' rank-one terms ' )
 2080 FORMAT( ' Row ', I6, ' Column ', I6, ' used from element ', I6,
     *        ' value = ', 1P, D24.16 )
 2090 FORMAT( ' Row ', I6, ' column ', I6, ' used. Value = ',1P,D24.16)
 2100 FORMAT( ' Group ', I5, ' second-order terms ' )
C
C  End of ASMBE.
C
      END
