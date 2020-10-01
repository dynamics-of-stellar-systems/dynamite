C  THIS VERSION: 04/08/1995 AT 11:38:29 AM.
C  ** VERSION B **
C ** Correction report.
C ** Correction 1. 27/03/2002: allow for non-zero defaults for QPs
      SUBROUTINE INLANC( N     , NLVARS, NG    , NELNUM, NOBJ  , LENGTH, 
     *                   LSTADG, LELVAR, LSTAEV, LNTVAR, LICNA , LSTADA, 
     *                   LA, LB, LBNDS , LINTRE, LIWK  , LWK   , NMAX  , 
     *                   NGMAX , NBMAX , NSMAX , NLMAX , NELMAX, NEGMAX, 
     *                   NOBMAX, NGRMAX, NGPVMX, NEPVMX, NOMAX , NLISGP, 
     *                   NBND  , NNZA  , NCONST, NSTART, NRANGE, NOBJGR, 
     *                   NOBBND, NELTYP, NGRTYP, 
     *                   PNAME , ONEOBJ, NAMEOB, NAMERH, NAMERA,
     *                   NAMEBN, NAMEST, NAMEOF, ISTADG, IELVAR, ISTAEV,
     *                   INTVAR, ICNA  , ISTADA, ICOORD, INLIST, ITABLE,
     *                   ISTATE, IDROWS, IELV  , IINV  , IGPA  , IELING,
     *                   ISTEP , ISTGP , ITYPEE, ITYPEG, ITYPEV, IWK   ,
     *                   A, BND, VSTART, CSTART, RSCALE, CSCALE, 
     *                   RDROWS, DFAULT, WEIGHT, BNDFLT, WK, 
     *                   GPVALU, EPVALU, FBOUND, ABYROW, B , BL, BU, X , 
     *                   CLMULT, ESCALE, GSCALE, VSCALE, INTREP, GXEQX ,
     *                   KEY   , GNAMES, VNAMES, BNAMES, SNAMES, ONAMES,
     *                   ETYPES, GTYPES, OBNAME, IALGOR, IAUTO ,
     *                   IOUT  , IOUTDA, SINGLE, INFORM, DEBUG  )
C
C  CONVERT THE OUTPUT FROM GPSMPS INTO A FORM SUITABLE FOR ANY OF THE
C  ------------------------------------------------------------------
C  LANCELOT PROGRAMS.
C  ------------------
C  SBMIN ( IALGOR = 1 ), AUGLG ( IALGOR = 2 ) OR BARIA ( IALGOR = 3 ).
C  -------------------------------------------------------------------
C
C  NICK GOULD, 16/01/1990
C  FOR CGT PRODUCTIONS.
C
      INTEGER          N , NG, LENGTH, LSTADG, LELVAR, LSTAEV, LNTVAR
      INTEGER          LSTADA, LA, LB, LBNDS , LINTRE, LIWK  , LWK
      INTEGER          NMAX  , NGMAX , NBMAX , NSMAX , NLMAX , NLVARS
      INTEGER          NEGMAX, NGRMAX, NGPVMX, NEPVMX, NLISGP, NOMAX
      INTEGER          NBND  , NNZA  , NCONST, NSTART, NRANGE, NOBJGR
      INTEGER          NELMAX, IALGOR, IOUT  , IOUTDA, INFORM, NOBBND
      INTEGER          NELTYP, NGRTYP, LICNA , NOBMAX, NOBJ  , NELNUM
      INTEGER          IAUTO
      LOGICAL          ONEOBJ, SINGLE, DEBUG
      CHARACTER * 8    PNAME
      CHARACTER *10    NAMEOB, NAMERH, NAMERA, NAMEBN, NAMEST,NAMEOF
      INTEGER          ISTADG( LSTADG ), IELVAR( LELVAR ), ICNA( LICNA )
      INTEGER          ISTAEV( LSTAEV ), INTVAR( LNTVAR )
      INTEGER          ISTADA( LSTADA ), ITABLE( LENGTH )
      INTEGER          ICOORD( LA, 2  ), INLIST( LENGTH )
      INTEGER          ISTATE( NGMAX  ), IDROWS( 2, NGMAX )
      INTEGER          IELV  ( NLMAX  ), IINV  ( NLMAX  )
      INTEGER          ISTEP ( NELMAX ), ISTGP ( NGMAX  )
      INTEGER          IELING( NEGMAX ), IGPA  ( NGRMAX ) 
      INTEGER          ITYPEE( NELMAX ), ITYPEG( NGMAX ), IWK( LIWK )
      INTEGER          ITYPEV( NMAX )
      DOUBLE PRECISION A( LA ),  RDROWS( 2, NGMAX ), WEIGHT( NEGMAX )
      DOUBLE PRECISION BND( 2, NMAX, NBMAX ), VSTART( NMAX, NSMAX )
      DOUBLE PRECISION BNDFLT( 2, NBMAX ), CSTART( NGMAX, NSMAX )
      DOUBLE PRECISION RSCALE( NGMAX ), CSCALE( NMAX )
      DOUBLE PRECISION GPVALU( NGPVMX ), EPVALU( NEPVMX )
      DOUBLE PRECISION DFAULT( NMAX ), FBOUND( 2, NOBMAX ), WK( LWK )
      DOUBLE PRECISION ABYROW( LA ), B( LB ), BL( LBNDS ), BU( LBNDS )
      DOUBLE PRECISION X( NMAX ), GSCALE( NGMAX ), ESCALE( NEGMAX )
      DOUBLE PRECISION VSCALE( NMAX ), CLMULT( NGMAX )
      LOGICAL          INTREP( LINTRE ), GXEQX( NGMAX )
      CHARACTER * 12 KEY( LENGTH  )
      CHARACTER * 10 GNAMES( NGMAX ), BNAMES( NBMAX ), ONAMES( NOMAX )
      CHARACTER * 10 VNAMES( NMAX  ), SNAMES( NSMAX ), OBNAME( NOBMAX )
      CHARACTER * 10 ETYPES( NLMAX ), GTYPES( NGRMAX )
      EXTERNAL        HASHC
C
C  LOCAL VARIABLES.
C
      INTEGER          I, IC, IFIELD, IG, IROW, IS, ITYPE, J, JBND, K
      INTEGER          K1, K2, NELV, NINV, NG1, NEL, NGPV, NGR
      INTEGER          JSTART, NEL1, JCOL, JCONST, JRANGE
      INTEGER          NNZ, NSLACK, JOBBND
      DOUBLE PRECISION AVALUE, BIG, RZERO, RROW
      DOUBLE PRECISION ONE, ZERO, BIGINF, OBFBND( 2 )
C ** Correction 1a. 1 line added
      CHARACTER * 10   CQGROU
      CHARACTER * 12   FIELD
      INTRINSIC        ABS, MAX
      EXTERNAL         REORDA
      PARAMETER      ( ZERO   = 0.0D+0,  ONE = 1.0D+0, RZERO = 0.0D+0 )
      PARAMETER      ( BIGINF = 1.0D+20, BIG = 1.0E+20 )
C ** Correction 1b. 1 line added
      PARAMETER      ( CQGROU = '123456789G' )
      IF ( LIWK .LT.  NG .OR. LWK .LT. MAX( N, NLISGP ) ) THEN
         INFORM = - 1
         IF ( IOUT  .GT. 0 ) WRITE( IOUT, 2000 )
         GO TO 800
      END IF
C
C  DECIDE WHICH OPTIMIZATION METHOD TO USE.
C  ---------------------------------------
C
      IF ( IALGOR .LE. 0 ) THEN
         IALGOR = 1
         DO 10 IG = 1, NG
            IF ( ABS( ISTATE( IG ) ) .GE. 2 ) IALGOR = 2
   10    CONTINUE   
      END IF
      DO 20 IG = 1, NG
         ISTATE( IG ) = ABS( ISTATE( IG ) )
   20 CONTINUE   
C
C  SELECT THE BOUNDS ON THE VARIABLES.
C  -----------------------------------
C
      IF ( NBND .EQ. 1 ) THEN
         JBND = 1
      ELSE
C
C  FIND THE KEY WORD IN THE LIST OF BOUNDS.
C
         DO 110 JBND = 1, NBND
            IF ( NAMEBN .EQ. BNAMES( JBND ) ) GO TO 120
  110    CONTINUE
         INFORM = 47
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2470 ) NAMEBN
         GO TO 800
  120    CONTINUE
      END IF
C
C  THE REQUIRED VECTOR OF BOUNDS IS COLUMN JBND OF BND. COPY THIS
C  VECTOR INTO BL( ) AND BU( ). RECORD THE SCALE FACTORS.
C
      DO 130 J       = 1, NLVARS
         BL( J )     = BND( 1, J, JBND )
         BU( J )     = BND( 2, J, JBND )
         VSCALE( J ) = CSCALE( J )
  130 CONTINUE
C
C  THE BOUNDS ON THE NONLINEAR VARIABLES ARE SET TO DEFAULT VALUES.
C
      DO 140 J   = NLVARS + 1, N
         BL( J ) = BNDFLT( 1, JBND )
         BU( J ) = BNDFLT( 2, JBND )
  140 CONTINUE
C
C  SELECT THE CONSTANT/RHS AND RANGES.
C  -----------------------------------
C
C  FIND THE NAMED CONSTANT (R.H.S.) VECTOR IN THE LIST.
C
      IF ( NCONST .EQ. 1 ) THEN
         JCONST = NLVARS + 1
      ELSE
         FIELD  = NAMERH // 'CO'
         CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
         IF ( IFIELD .GT. 0 ) THEN
            JCONST = INLIST( IFIELD )
         ELSE
            INFORM = 46
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2460 ) NAMERH
            GO TO 800
         END IF
      END IF
C
C  FIND THE NAMED RANGE VECTOR IN THE LIST.
C
      IF ( NRANGE .GT. 0 ) THEN
         IF ( NRANGE .EQ. 1 ) THEN
            JRANGE = NLVARS + NCONST + 1
         ELSE
            FIELD  = NAMERA // 'RA'
            CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
            IF ( IFIELD .GT. 0 ) THEN
               JRANGE = INLIST( IFIELD )
            ELSE
               INFORM = 48
               IF ( IOUT .GT. 0 ) WRITE( IOUT, 2480 ) NAMERA
               GO TO 800
            END IF
         END IF
      ELSE
         JRANGE = 0
      END IF
C
C  INITIALIZE THE VECTOR OF CONSTANTS, B, AS ITS DEFAULT VALUE.
C
      DO 210 I  = 1, NG
C ** Correction 1c. 6 lines added
        IF ( GNAMES ( I ) .EQ. CQGROU ) THEN
         B( I ) = ZERO
         BL( N + I ) = ZERO
         BU( N + I ) = BIGINF
         GSCALE( I ) = ONE
        ELSE
         B( I ) = DFAULT( JCONST )
C
C  INITIALIZE LOWER AND UPPER BOUNDS ON SLACK VARIABLES AS ZERO
C  AND THE DEFAULT RESPECTIVELY.
C
         BL( N + I ) = ZERO
         IF ( NRANGE .EQ. 0 ) THEN
            BU( N + I ) = BIGINF
         ELSE
            BU( N + I ) = DFAULT( JRANGE )
         END IF
C
C  RECORD THE GROUP SCALE FACTORS.
C
         GSCALE( I ) = ONE / RSCALE( I )
C ** Correction 1d. 1 line added
        END IF
  210 CONTINUE
C
C  SWEEP THROUGH THE ENTRIES OF A. LOOK FOR ENTRIES IN COLUMNS
C  JCONST AND JRANGE. SUBSEQUENTLY REMOVE ALL CONSTANT/RHS AND
C  RANGE COLUMNS TO LEAVE ONLY ENTRIES CORRESPONDING TO LINEAR
C  ELEMENTS.
C
      NNZ       = 0
      DO 220 K  = 1, NNZA
         I      = ICOORD( K, 1 )
         J      = ICOORD( K, 2 )
         AVALUE = A( K )
C
C  SEE IF THE ENTRY BELONGS TO THE SELECTED CONSTANT/RHS VECTOR.
C
         IF ( J .EQ. JCONST ) B( I ) = AVALUE
C
C  SEE IF THE ENTRY BELONGS TO THE SELECTED RANGE VECTOR.
C
         IF ( J .EQ. JRANGE ) BU( N + I ) = AVALUE
C
C  CHECK IF THE ENTRY BELONGS TO A LINEAR ELEMENT.
C
         IF ( J .LE. NLVARS ) THEN
            NNZ = NNZ + 1
C
C  RECORD THE COORDINATES AND VALUE OF THE ENTRY FROM THE LINEAR
C  ELEMENT.
C
            ICOORD( NNZ, 1 ) = I
            ICOORD( NNZ, 2 ) = J
            A( NNZ )         = AVALUE
         END IF
  220 CONTINUE
      NNZA = NNZ
C
C  THE MATRIX IS STORED IN COORDINATE FORM. RESORT IT SO THAT
C  IT IS STORED BY ROWS.
C
      IF ( NNZA .GT. 0 ) CALL REORDA( NG, NNZA, ICOORD( 1, 2 ),
     *                                ICOORD( 1, 1 ), A, ISTADA, IWK )
C
C  DECODE THE 'D'-GROUPS/ROWS. SET THE WORKSPACE ARRAY WK TO ZERO.
C
      NNZ        = 0
      NSLACK     = 0
      DO 310 I   = 1, N
         WK( I ) = RZERO
  310 CONTINUE
C
C  GXEQX IS TRUE IF THE GROUP FUNCTION IS TRIVIAL.
C
      DO 315 IG     = 1, NG
        GXEQX( IG ) = ITYPEG( IG ) .EQ. 0
  315 CONTINUE
C
C  SET THE COEFFICIENTS OF THE LINEAR ELEMENTS.
C  --------------------------------------------
C
C  CONSIDER THE GROUPS IN ORDER.
C
      DO 390 IG       = 1, NG
         K1           = ISTADA( IG )
         ISTADA( IG ) = NNZ + 1
         IF ( ISTATE( IG ) .LE. 4 ) THEN
C
C  First pass: determine the nonzeros in the row
C
            DO 320 K = k1, ISTADA( IG + 1 ) - 1
               J     = ICOORD( K, 2 )
               IF ( WK( J ) .NE. RZERO ) THEN
                  WK( J ) = WK( J ) + A( K )
               ELSE
                  WK( J ) = A( K )
               END IF
  320       CONTINUE
C
C  Second pass: only record nonzeros
C
            DO 325 K = K1, ISTADA( IG + 1 ) - 1
               J     = ICOORD( K, 2 )
               IF ( WK( J ) .NE. RZERO ) THEN
                  NNZ           = NNZ + 1
                  ICNA( NNZ )   = J
                  ABYROW( NNZ ) = WK( J )
                  WK( J ) = RZERO
C               ELSE
C                 write(6,*) ' duplicate or zero removed ', J, A( K )
               END IF
  325       CONTINUE
         ELSE
C
C  THE IG-TH GROUP IS A 'D'-GROUP. CONSTRUCT THE NEW GROUP FROM ITS
C  TWO DONORS. CONSIDER EACH DONOR ROW IN TURN. FORM THE NEW ROW IN WK.
C
            DO 340 I    = 1, 2
               IROW     = IDROWS( I, IG )
               RROW     = RDROWS( I, IG )
               DO 330 K = ISTADA( IROW ), ISTADA( IROW + 1 ) - 1
                  IC = ICNA( K )
                  IF ( IC .LE. NLVARS ) WK( IC ) = WK( IC ) +
     *                             ABYROW( K ) * RROW
  330          CONTINUE
  340       CONTINUE
C
C  MOVE THE NEW ROW INTO ABYROW, RESETTING WK TO ZERO AS WE PROCEED.
C
            DO 360 I    = 1, 2
               IROW     = IDROWS( I, IG )
               DO 350 K = ISTADA( IROW ), ISTADA( IROW + 1 ) - 1
                  IC    = ICNA( K )
                  IF ( IC .LE. NLVARS ) THEN
                     IF ( WK( IC ) .NE. RZERO ) THEN
                        NNZ           = NNZ + 1
                        ICNA( NNZ )   = IC
                        ABYROW( NNZ ) = WK( IC )
                        WK( IC )      = RZERO
                     END IF
                  END IF
  350          CONTINUE
  360       CONTINUE
         END IF
C
C  IF THE GROUP IS OF TYPE 'L' OR 'G', INSERT A SLACK VARIABLE IN
C  THE LINEAR ELEMENT.
C
         IS = ISTATE( IG )
         IF ( IS .GT. 4 ) IS = IS - 4
         IF ( IALGOR .LE. 2 ) THEN
            IF ( IS .NE. 1 .AND. IS .NE. 2 ) THEN
               IF ( .NOT. GXEQX( IG ) ) THEN
                 WRITE( IOUT, 1990 )
                 INFORM = 51
                 GO TO 800
               END IF
               NNZ         = NNZ + 1
               NSLACK      = NSLACK  + 1
               JCOL        = N + NSLACK
               ICNA( NNZ ) = JCOL
               IF ( IS .EQ. 3 ) THEN
                  ABYROW( NNZ ) =   ONE
               ELSE
                  ABYROW( NNZ ) = - ONE
               END IF
C
C  GIVE THE SLACK VARIABLE THE SAME NAME AS ITS CORRESPONDING GROUP.
C
               VNAMES( JCOL ) = GNAMES( IG )
C
C  ASSIGN THE CORRECT BOUNDS FOR THE SLACK VARIABLE.
C
               BL( JCOL ) = BL( N + IG )
               BU( JCOL ) = BU( N + IG )
            END IF
         ELSE
            IF ( IS .EQ. 3 ) BL( N + IG ) = - BU( N + IG )
            IF ( IS .EQ. 1 .OR. IS .EQ. 2 .OR. IS .EQ. 3 )
     *                       BU( N + IG ) = ZERO
         END IF
         ISTATE( IG ) = IS
  390 CONTINUE
C
C  RESET THE NUMBER OF VARIABLES TO INCLUDE THE SLACKS.
C
      N             = N + NSLACK
      NNZA          = NNZ + 1
      NG1           = NG + 1
      ISTADA( NG1 ) = NNZA
      IF ( DEBUG .AND. IOUT .GT. 0 ) WRITE( IOUT, 3020 )
     *    ( ( I, ICNA( K ), ABYROW( K ), K = ISTADA( I ),
     *        ISTADA( I + 1 ) - 1 ), I = 1, NG )
C
C  SELECT THE STARTING POINT FOR THE MINIMIZATION.
C  ----------------------------------------------
C
      IF ( NSTART .EQ. 1 ) THEN
         JSTART = 1
      ELSE
C
C  FIND THE KEY WORD IN THE LIST OF STARTING POINTS.
C
         DO 410 JSTART = 1, NSTART
            IF ( NAMEST .EQ. SNAMES( JSTART ) ) GO TO 420
  410    CONTINUE
         INFORM = 49
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2490 ) NAMEST
         GO TO 800
  420    CONTINUE
      END IF
C
C  RECORD THE STARTING POINT.
C
      DO 430 J  = 1, NLVARS
         X( J ) = VSTART( J, JSTART )
  430 CONTINUE
C
C  INITIALIZE ALL SLACK AND NONLINEAR VARIABLES AS ZERO, WITH WEIGHT 1.
C
      DO 440 J       = NLVARS + 1, N
         X( J )      = ZERO
         VSCALE( J ) = ONE
  440 CONTINUE
C
C  RECORD THE LAGRANGE MULTIPLIERS, IF ANY.
C
      IF ( IALGOR .GE. 2 ) THEN
         DO 450 IG       = 1, NG
            CLMULT( IG ) = CSTART( IG, JSTART )
  450    CONTINUE
      END IF
C
C  NONLINEAR ELEMENT INFORMATION.
C  ------------------------------
C
C  THE PARAMETER VALUES FOR THE NONLINEAR ELEMENTS MAY BE UNORDERED.
C  ORDER THE PARAMETERS SO THAT THOSE FOR GROUP I PRECEDE THOSE
C  FROM GROUP I+1. PLACE THE REORDERED SET IN WK.
C
      J         = 1
      DO 520 IG = 1, NG
         K      = ITYPEG( IG )
         ISTGP( IG ) = J
         IF ( K .GT. 0 ) THEN
            K1 = IGPA( K + 1 ) - IGPA( K )
            K2 = ISTGP( IG ) - 1
            DO 510 I    = 1, K1
               WK( J )  = GPVALU( K2 + I )
               J        = J + 1
  510       CONTINUE
         END IF
  520 CONTINUE
      ISTGP( NG1 ) = J
C
C  OVERWRITE GPVALU WITH WK.
C
      DO 530 I = 1, J - 1
         GPVALU( I ) = WK( I )
  530 CONTINUE
C
C  RECORD THE SCALE FACTORS FOR THE NONLINEAR ELEMENTS.
C
      DO 540 I       = 1, ISTADG( NG1 ) - 1
         ESCALE( I ) = WEIGHT( I )
  540 CONTINUE
      DO 550 I = 1, NELNUM
C
C  DETERMINE WHETHER THE NONLINEAR ELEMENTS HAVE INTERNAL
C  REPRESENTATIONS.
C
         ITYPE       = ITYPEE( I )
         NELV        = IELV( ITYPE + 1 ) - IELV( ITYPE )
         NINV        = IINV( ITYPE + 1 ) - IINV( ITYPE )
         INTREP( I ) = NINV .LT. NELV
C
C  STORE THE NUMBER OF INTERNAL VARIABLES FOR EACH ELEMENT.
C
         INTVAR( I ) = NINV
  550 CONTINUE
      IF ( IALGOR .GE. 2 ) THEN
C
C  SELECT THE OBJECTIVE FUNCTION GROUP.
C  ------------------------------------
C
C  FIND THE NAMED OBJECTIVE FUNCTION GROUP IN THE LIST.
C  MARK THE REMAINING OBJECTIVE FUNCTION GROUPS FOR REMOVAL.
C
         IF ( ONEOBJ ) THEN
            NOBJGR    = 0
            DO 610 I  = 1, NOBJ
               FIELD  = ONAMES( I ) // 'GR'
               CALL HASHC ( LENGTH, 12, FIELD, KEY, ITABLE, IFIELD )
               K = INLIST( IFIELD )
               IF ( NAMEOF .EQ. ONAMES( I ) ) THEN
                  NOBJGR = K
                  IF ( IOUT .GT. 0 .AND. DEBUG )
     *                 WRITE( IOUT, 3010 ) ONAMES( I ), K
               ELSE
                  ISTADG( K ) = - ISTADG( K )
                  IF ( IOUT .GT. 0 .AND. DEBUG )
     *                 WRITE( IOUT, 3000 ) ONAMES( I )
               END IF
  610       CONTINUE
C
C  REMOVE REDUNDANT GROUP INFORMATION.
C
             IF ( NOBJ .GT. 1 .OR.
     *          ( NOBJ .EQ. 1 .AND. NOBJGR .EQ. 0 ) ) THEN
               NNZ      = 1
               NEL      = 1
               NGPV     = 1
               NGR      = 0
               DO 650 I = 1, NG
                  IF ( ISTADG( I ) .GT. 0 ) THEN
                     NGR = NGR + 1
                     IF ( I .EQ. NOBJGR ) NOBJGR = NGR
C
C  SHIFT THE GROUP STATUS, NAME, TYPE, CONSTANT,
C  TRIVIALITY INDICATOR, LAGRANGE MULTIPLIER AND WEIGHT.
C
                     IF ( IALGOR .GE. 2 ) ISTATE( NGR ) = ISTATE( I )
                     GNAMES( NGR ) = GNAMES( I )
                     ITYPEG( NGR ) = ITYPEG( I )
                     B( NGR )      = B( I )
                     GXEQX( NGR )  = GXEQX( I )
                     IF ( IALGOR .GE. 2 ) CLMULT( NGR ) = CLMULT( I )
                     GSCALE( NGR ) = GSCALE( I )
C
C  SHIFT THE LIST OF ELEMENTS AND WEIGHTS IN THE I-TH GROUP.
C
                     K1               = ISTADG( I )
                     ISTADG( NGR )    = NEL
                     DO 620 K         = K1, ABS( ISTADG( I + 1 ) ) - 1
                        IELING( NEL ) = IELING( K )
                        ESCALE( NEL ) = ESCALE( K )
                        NEL           = NEL + 1
  620                CONTINUE
C
C  SHIFT THE LIST OF PARAMETERS IN THE I-TH GROUP.
C
                     K1                = ISTGP( I )
                     ISTGP( NGR )      = NGPV
                     DO 630 K          = K1, ISTGP( I + 1 ) - 1
                        GPVALU( NGPV ) = GPVALU( K )
                        NGPV           = NGPV + 1
  630             CONTINUE
C
C  SHIFT THE LIST OF COEFFICIENTS AND POSITIONS OF THE NONZEROS
C  FOR THE LINEAR ELEMENT IN THE I-TH GROUP.
C
                     K1             = ISTADA( I )
                     ISTADA( NGR )  = NNZ
                     DO 640 K       = K1, ISTADA( I + 1 ) - 1
                        A( NNZ )    = A( K )
                        ICNA( NNZ ) = ICNA( K )
                        NNZ         = NNZ + 1
  640                CONTINUE
                  END IF
  650          CONTINUE
               NG               = NGR
               ISTADG( NG + 1 ) = NEL
               ISTGP ( NG + 1 ) = NGPV
               ISTADA( NG + 1 ) = NNZ
            END IF
         END IF
      END IF
C
C  SET THE REQUIRED LOWER AND UPPER BOUNDS ON THE OBJECTIVE FUNCTION.
C
      IF ( NOBBND .EQ. 0 ) THEN
         OBFBND( 1 ) = - BIGINF
         OBFBND( 2 ) =   BIGINF
      ELSE
C
C  FIND THE KEY WORD IN THE LIST OF STARTING POINTS.
C
         DO 660 JOBBND = 1, NOBBND
            IF ( NAMEOB .EQ. OBNAME( JOBBND ) ) GO TO 670
  660    CONTINUE
         INFORM = 50
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2500 ) NAMEOB
         GO TO 800
  670    CONTINUE
         OBFBND( 1 ) = FBOUND( 1, JOBBND )
         OBFBND( 2 ) = FBOUND( 2, JOBBND )
      END IF
C
C  IF NO OUTPUT IS REQUIRED, EXIT.
C
      IF ( IOUTDA .LE. 0 ) GO TO 900
      NEL1 = NELNUM + 1
      WRITE( IOUTDA, 3180 ) N,  NG, NELNUM,     ISTADG( NG1  ) - 1,
     *                      ISTAEV( NEL1 ) - 1, ISTADA( NG1  ) - 1,
     *                      ISTGP ( NG1  ) - 1, ISTEP ( NEL1 ) - 1,
     *                      NELTYP, NGRTYP 
C
C  PRINT OUT PROBLEM DATA. OUTPUT THE NUMBER OF VARIABLES, GROUPS AND
C  ELEMENTS AND, PERHAPS, THE IDENTITY OF THE OBJECTIVE FUNCTION GROUP.
C
      WRITE( IOUTDA, 3100 ) IALGOR, PNAME, IAUTO
      IF ( IALGOR .EQ. 2 ) WRITE( IOUTDA, 3170 ) NSLACK, NOBJGR
C
C  OUTPUT THE STARTING ADDRESSES OF THE ELEMENTS IN EACH GROUP,
C  OF THE PARAMETERS USED FOR EACH GROUP AND
C  OF THE NONZEROS OF THE LINEAR ELEMENT IN EACH GROUP.
C
      WRITE( IOUTDA, 3110 ) ( ISTADG( I ), I = 1, NG1 )
      WRITE( IOUTDA, 3110 ) ( ISTGP ( I ), I = 1, NG1 )
      WRITE( IOUTDA, 3110 ) ( ISTADA( I ), I = 1, NG1 )
C
C  OUTPUT THE STARTING ADDRESSES OF THE VARIABLES AND PARAMETERS
C  IN EACH ELEMENT.
C
      WRITE( IOUTDA, 3110 ) ( ISTAEV( I ), I = 1, NEL1 )
      WRITE( IOUTDA, 3110 ) ( ISTEP( I ), I = 1, NEL1 )
C
C  OUTPUT THE GROUP TYPE OF EACH GROUP AND ITS STATUS.
C
      WRITE( IOUTDA, 3110 ) ( ITYPEG( I ), I = 1, NG )
      IF ( IALGOR .GE. 2 ) WRITE( IOUTDA, 3110 )( ISTATE( I ), I = 1,NG)
C
C  OUTPUT THE ELEMENT TYPE OF EACH ELEMENT.
C
      WRITE( IOUTDA, 3110 ) ( ITYPEE( I ), I = 1, NELNUM )
C
C  OUTPUT THE ELEMENT TYPE OF EACH ELEMENT
C  AND ITS NUMBER OF INTERNAL VARIABLES.
C
      WRITE( IOUTDA, 3110 ) ( INTVAR( I ), I = 1, NELNUM )
C
C  OUTPUT THE IDENTITY OF EACH INDIVIDUAL ELEMENT.
C
      WRITE( IOUTDA, 3110 ) ( IELING( I ), I = 1, ISTADG( NG1 ) - 1 )
C
C  OUTPUT THE VARIABLES IN EACH GROUP'S ELEMENTS.
C
      WRITE( IOUTDA, 3110 ) ( IELVAR( I ), I = 1, ISTAEV( NEL1 ) - 1 )
C
C  OUTPUT THE COLUMN ADDRESSES OF THE NONZEROS IN EACH LINEAR ELEMENT.
C
      NNZA = ISTADA( NG1 ) - 1
      WRITE( IOUTDA, 3110 ) ( ICNA( I ), I = 1, NNZA )
C
C  WRITE SINGLE PRECISION FORMAT
C
      IF ( SINGLE ) THEN
C
C  OUTPUT THE VALUES OF THE NONZEROS IN EACH LINEAR ELEMENT, THE
C  CONSTANT TERM IN EACH GROUP, THE LOWER AND UPPER BOUNDS ON
C  THE VARIABLES AND THE STARTING POINT FOR THE MINIMIZATION
C
         WRITE( IOUTDA, 3121 ) ( ABYROW( I ), I = 1, NNZA )
         WRITE( IOUTDA, 3121 ) ( B( I ), I = 1, NG )
         IF ( IALGOR .LE. 2 ) THEN
            WRITE( IOUTDA, 3121 ) ( BL( I ), I = 1, N )
            WRITE( IOUTDA, 3121 ) ( BU( I ), I = 1, N )
         ELSE
            WRITE( IOUTDA, 3121 ) ( BL( I ), I = 1, N + NG )
            WRITE( IOUTDA, 3121 ) ( BU( I ), I = 1, N + NG )
         END IF
         WRITE( IOUTDA, 3121 ) ( X( I ), I = 1, N )
         IF ( IALGOR .GE. 2 ) WRITE( IOUTDA, 3121 )( CLMULT( I ), 
     *                                               I = 1, NG )
C
C  OUTPUT THE PARAMETERS IN EACH GROUP.
C
         WRITE( IOUTDA, 3121 ) ( GPVALU( I ), I = 1, ISTGP( NG1 ) - 1 )
C
C  OUTPUT THE PARAMETERS IN EACH INDIVIDUAL ELEMENT.
C
         WRITE( IOUTDA, 3121 ) ( EPVALU( I ), I = 1, ISTEP( NEL1 ) - 1 )
C
C  OUTPUT THE SCALE FACTORS FOR THE NONLINEAR ELEMENTS.
C
         WRITE( IOUTDA, 3121 ) ( ESCALE( I ), I = 1, ISTADG( NG1 ) - 1 )
C
C  OUTPUT THE SCALE FACTORS FOR THE GROUPS.
C
         WRITE( IOUTDA, 3121 ) ( GSCALE( I ), I = 1, NG )
C
C  OUTPUT THE SCALE FACTORS FOR THE VARIABLES.
C
         WRITE( IOUTDA, 3121 ) ( VSCALE( I ), I = 1, N )
C
C  OUTPUT THE LOWER AND UPPER BOUNDS ON THE OBJECTIVE FUNCTION.
C
         WRITE( IOUTDA, 3161 ) OBFBND( 1 ), OBFBND( 2 )
C
C  WRITE DOUBLE PRECISION FORMAT
C
      ELSE
C
C  OUTPUT THE VALUES OF THE NONZEROS IN EACH LINEAR ELEMENT, THE
C  CONSTANT TERM IN EACH GROUP, THE LOWER AND UPPER BOUNDS ON
C  THE VARIABLES AND THE STARTING POINT FOR THE MINIMIZATION
C
         WRITE( IOUTDA, 3120 ) ( ABYROW( I ), I = 1, NNZA )
         WRITE( IOUTDA, 3120 ) ( B( I ), I = 1, NG )
         IF ( IALGOR .LE. 2 ) THEN
            WRITE( IOUTDA, 3120 ) ( BL( I ), I = 1, N )
            WRITE( IOUTDA, 3120 ) ( BU( I ), I = 1, N )
         ELSE
            WRITE( IOUTDA, 3120 ) ( BL( I ), I = 1, N + NG )
            WRITE( IOUTDA, 3120 ) ( BU( I ), I = 1, N + NG )
         END IF
         WRITE( IOUTDA, 3120 ) ( X( I ), I = 1, N )
         IF ( IALGOR .GE. 2 ) WRITE( IOUTDA, 3120 )( CLMULT( I ), 
     *                                               I = 1, NG )
C
C  OUTPUT THE PARAMETERS IN EACH GROUP.
C
         WRITE( IOUTDA, 3120 ) ( GPVALU( I ), I = 1, ISTGP( NG1 ) - 1 )
C
C  OUTPUT THE PARAMETERS IN EACH INDIVIDUAL ELEMENT.
C
         WRITE( IOUTDA, 3120 ) ( EPVALU( I ), I = 1, ISTEP( NEL1 ) - 1 )
C
C  OUTPUT THE SCALE FACTORS FOR THE NONLINEAR ELEMENTS.
C
         WRITE( IOUTDA, 3120 ) ( ESCALE( I ), I = 1, ISTADG( NG1 ) - 1 )
C
C  OUTPUT THE SCALE FACTORS FOR THE GROUPS.
C
         WRITE( IOUTDA, 3120 ) ( GSCALE( I ), I = 1, NG )
C
C  OUTPUT THE SCALE FACTORS FOR THE VARIABLES.
C
         WRITE( IOUTDA, 3120 ) ( VSCALE( I ), I = 1, N )
C
C  OUTPUT THE LOWER AND UPPER BOUNDS ON THE OBJECTIVE FUNCTION.
C
         WRITE( IOUTDA, 3160 ) OBFBND( 1 ), OBFBND( 2 )
      END IF
C
C  OUTPUT A LOGICAL ARRAY WHICH SAYS WHETHER AN ELEMENT HAS INTERNAL
C  VARIABLES.
C
      WRITE( IOUTDA, 3130 ) ( INTREP( I ), I = 1, NELNUM )
C
C  OUTPUT A LOGICAL ARRAY WHICH SAYS WHETHER A GROUP IS TRIVIAL.
C
      WRITE( IOUTDA, 3130 ) ( GXEQX( I ), I = 1, NG )
C
C  OUTPUT THE NAMES GIVEN TO THE GROUPS AND TO THE VARIABLES.
C
      WRITE( IOUTDA, 3140 ) ( GNAMES( I ), I = 1, NG )
      WRITE( IOUTDA, 3140 ) ( VNAMES( I ), I = 1, N )
C
C  OUTPUT THE NAMES GIVEN TO THE ELEMENT AND GROUP TYPES.
C
      WRITE( IOUTDA, 3140 ) ( ETYPES( I ), I = 1, NELTYP )
      WRITE( IOUTDA, 3140 ) ( GTYPES( I ), I = 1, NGRTYP )
C
C  OUTPUT THE TYPE OF EACH VARIABLE.
C
      WRITE( IOUTDA, 3110 ) ( ITYPEV( I ), I = 1, N )
      GO TO 900
C
C  INCORRECT DATA SPECIFIED.
C
  800 CONTINUE
      RETURN
C
C  SUCCESSFUL RETURN
C
  900 CONTINUE
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 1990 FORMAT( ' ** Exit from INLANC - ', /, 
     *        ' Although the manual may suggest otherwise,',
     *        ' non-trivial',/,
     *        ' groups are not allowed for inequality constraints',/)
 2000 FORMAT( ' ** Exit from INLANC - increase one or more of NGPVMX,
     *          NELMAX and NGMAX ' )
 2460 FORMAT( ' ** Exit from INLANC - constant name ', A8,
     *        ' not recognised ' )
 2470 FORMAT( ' ** Exit from INLANC - bound name ', A8,
     *        ' not recognised ' )
 2480 FORMAT( ' ** Exit from INLANC - range name ', A8,
     *        ' not recognised ' )
 2490 FORMAT( ' ** Exit from INLANC - start point name ', A8,
     *        ' not recognised ' )
 2500 FORMAT( ' ** Exit from INLANC - obj. bound name ', A8,
     *        ' not recognised ' )
 3000 FORMAT( ' Group ', A10, ' removed as a redundant objective ' )
 3010 FORMAT( /, ' Objective function ', A10, ' is group number ', I8 )
 3020 FORMAT( /, 3('  Row   Col    Value  '),
     *        /, 3('  ---   ---    -----  '),
     *        /, ( 3( 2I5, 1P, D12.4 ) ) )
 3100 FORMAT( I2, A8, I2 )
 3110 FORMAT( ( 10I8 ) )
 3120 FORMAT( ( 1P, 4D16.8 ) )
 3121 FORMAT( ( 1P, 4E16.8 ) )
 3130 FORMAT( ( 72L1 ) )
 3140 FORMAT( ( 8A10 ) )
 3160 FORMAT( 1P, 2D16.8 )
 3161 FORMAT( 1P, 2E16.8 )
 3170 FORMAT( 2I8 )
 3180 FORMAT( 10I8 )
C
C  END OF INLANC.
C
      END
