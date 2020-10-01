C  THIS VERSION: 26/02/2001 AT 09:30:00 AM.
C     ( Last modified on 15 Mar 2001 at 22:28:00 )
C ** Correction report.
C ** Correction -1. 03/03/00: Integer formats increased
C ** Correction 0. 20/12/99: Array holding variable types introduced.
C ** Correction 1. 13/01/94: 1 line corrected **
C ** Correction 2. 26/02/01: 3 dummy arguments removed **
C ** End of Correction report.
C ** Correction 2. 26/02/01: 5 dummy arguments removed **
      SUBROUTINE PRINTP( NMAX, NGMAX, NLMAX,
     *                   NELMAX, NETMAX,
     *                   NEVMAX, NEPMAX, NGRMAX, NEGMAX, NEPVMX,
     *                   NGPVMX, NGPMAX, LSTADA, LICNA, LIWK,
     *                   N, NG, NLVARS, NELNUM,
     *                   ISTATE, ISTADG, IELVAR, ITYPEG, ITYPEE,
     *                   IELV, IINV, IEPA, IGPA,
     *                   ISTADA, ICNA, ISTGP, ISTEP, ISTEV,
C ** Correction 0. 20/12/99: Array holding variable types introduced.
     *                   IELING, ITYPEV, IWK,
     *                   ABYROW, B, BL, BU, X, EPVALU, GPVALU,
     *                   GSCALE, ESCALE, VSCALE,
     *                   PNAME, VNAMES, GNAMES,
     *                   LNAMES, ETYPES, ENAMES,
     *                   ANAMES, EPNAME, GPNAME, GTYPES,
     *                   IOUT, IPRINT )
      INTEGER       IOUT, IPRINT
      INTEGER       NMAX, NGMAX, NLMAX, NELMAX, NETMAX
      INTEGER       NEPMAX, NEVMAX, NGRMAX, LIWK
      INTEGER       N, NG
      INTEGER       NLVARS, NELNUM, LICNA, LSTADA
      INTEGER       NEGMAX, NEPVMX, NGPVMX, NGPMAX
      DOUBLE PRECISION EPVALU( NEPVMX ), GPVALU( NGPVMX )
      INTEGER       IELING( NEGMAX )
      INTEGER       ISTATE( NGMAX ), ISTADA( LSTADA ), ICNA( LICNA )
      INTEGER       IELV  ( NLMAX ), IINV  ( NLMAX )
C ** Correction 0. 20/12/99: Array holding variable types introduced.
      INTEGER       IELVAR( NEVMAX ), ITYPEV( NMAX )
      INTEGER       ISTADG(  NGMAX ), ITYPEG( NGMAX ), ITYPEE( NELMAX )
      INTEGER       IEPA( NLMAX ), IGPA( NGRMAX ), IWK( LIWK )
C ** Correction 1. 13/01/94: 1 line corrected **
      INTEGER       ISTEP( NELMAX ), ISTEV( NELMAX ), ISTGP( NGMAX )
C ** Correction 1. 13/01/94:  end of correction **
      CHARACTER * 8  PNAME
      CHARACTER * 10 GNAMES( NGMAX ), VNAMES( NMAX )
      CHARACTER * 10 ETYPES( NLMAX ), LNAMES( NELMAX )
      CHARACTER * 10 ENAMES( NETMAX )
      CHARACTER * 10 ANAMES( NGRMAX ), GTYPES( NGRMAX )
      CHARACTER * 10 EPNAME( NEPMAX ), GPNAME( NGPMAX )
      DOUBLE PRECISION B( NGMAX ), BL( NMAX ), BU( NMAX ), X( NMAX )
      DOUBLE PRECISION GSCALE( NGMAX ), ESCALE( NEGMAX ), VSCALE( NMAX )
      DOUBLE PRECISION ABYROW( LICNA )
C
C  PRINT DETAILS OF THE PROBLEM PREVIOUSLY SPECIFIED IN AN SIF FILE.
C  ------------------------------------------------------------------
C
C THE LEVEL OF PRINTING PERFORMED IS DETERMINED BY THE VALUE OF IPRINT.
C  POSSIBLE VALUES ARE:
C
C  >= 1, A SIMPLE SUMMARY OF THE PROBLEM NAME, THE NUMBER OF VARIABLES,
C        GROUPS AND ELEMENTS.
C  >= 2, A LIST OF THE VARIABLES USED.
C  >= 3, A BREAKDOWN OF THE GROUPS. A LIST OF THE NONLINEAR ELEMENTS
C        USED, THE TYPES OF EACH GROUP, THE STATUS OF THE GROUP AND
C        A STATEMENT THAT THE GROUP USES OR DOES NOT USE A LINEAR
C        ELEMENT.
C   = 4, FURTHER DETAILS OF EACH GROUP. THE NAME OF THE GROUP-TYPE
C        VARIABLE AND A LIST OF THE VALUES ASSOCIATED WITH EACH
C        PARAMETER.
C  >= 5, DETAILS OF EACH ELEMENT. THE NUMBERS OF ELEMENTAL AND
C        INTERNAL VARIABLES AND THE NUMBERS OF PARAMETERS.
C  >= 6, FURTHER DETAILS OF EACH ELEMENT. THE NAMES OF THE
C        ELEMENTAL VARIABLES TOGETHER WITH THEIR ASSOCIATED
C        PROBLEM VARIABLES, A LIST OF THE VALUES ASSOCIATED
C        WITH EACH PARAMETER AND THE VARIABLES INVOLVED IN
C        THE LINEAR ELEMENT.
C  >= 7, DETAILS OF THE COEFFICIENTS OF THE LINEAR ELEMENTS.
C  >= 8, FULL DETAILS OF THE VARIABLES USED INCLUDING THEIR LOWER
C        AND UPPER BOUNDS AND STARTING VALUES.
C  >= 9, ALL OF THE ABOVE.
C
C  NICK GOULD 28/07/1989
C  FOR CGT PRODUCTIONS.
C
      INTEGER I, IEL, IG, IS, J, K, K1, K2, K3, K4, K5, K6, L
      INTEGER IELTYP
      INTRINSIC ABS
      CHARACTER * 12 STATUS( 4 )
C ** Correction 0. 20/12/99: Array holding variable types introduced.
      CHARACTER * 3 VARTYP( 3 )
      DATA VARTYP / '   ', '0-1', 'int' /
      DATA STATUS / 'an objective', 'an equality ',
     *              'a negativity', 'a positivity' /
     *              
      IF ( IOUT .LE. 0 ) RETURN
      IF ( IPRINT .GE. 1 ) THEN
         DO 10 I     = 1, NELNUM
            IWK( I ) = 0
   10    CONTINUE
         WRITE( IOUT, 2000 ) PNAME, N, N - NLVARS, NG, NELNUM
C
C  LIST OF VARIABLES.
C
         IF ( IPRINT .GE. 2 ) THEN
             WRITE( IOUT, 2010 ) ( VNAMES( J ), J = 1, NLVARS )
             IF ( NLVARS .LT. N )
     *          WRITE( IOUT, 2020 ) ( VNAMES( J ), J = NLVARS + 1, N )
C
C  GROUP DETAILS.
C
            IF ( IPRINT .GE. 3 ) THEN
               DO 500 IG = 1, NG
                  IS     = ITYPEG( IG )
                  IF ( IS .EQ. 0 ) THEN
                     IF ( ISTATE( IG ) .EQ. 1 ) THEN
                       WRITE( IOUT, 2030 ) IG, GNAMES( IG ),
     *                 STATUS( 1 ), 'TRIVIAL   ', GSCALE( IG )
                     ELSE
                       WRITE( IOUT, 2030 ) IG, GNAMES( IG ),
     *                 STATUS( ISTATE( IG ) ),
     *                 'TRIVIAL   ', GSCALE( IG )
                     END IF
                  ELSE
                     IF ( ABS( ISTATE( IG ) ) .EQ. 1 ) THEN
                       WRITE( IOUT, 2030 ) IG, GNAMES( IG ),
     *                 STATUS( 1 ), GTYPES( IS ), GSCALE( IG )
                     ELSE
                       WRITE( IOUT, 2030 ) IG, GNAMES( IG ),
     *                 STATUS( ISTATE( IG ) ),
     *                 GTYPES( IS ), GSCALE( IG )
                     END IF
                  END IF
                  K1 = ISTADG( IG )
                  K2 = ISTADG( IG + 1 ) - 1
                  L  = K2 - K1 + 1
                  IF ( K1 .LE. K2 ) THEN
                     IF ( K1 .EQ. K2 ) THEN
                        WRITE( IOUT, 2060 ) LNAMES( IELING( K1 ) )
                     ELSE
                        WRITE( IOUT, 2070 ) L,
     *                       ( LNAMES( IELING( K ) ), K = K1, K2 )
                     END IF
                  ELSE
                    WRITE( IOUT, 2080 )
                  END IF
C
C  FURTHER GROUP DETAILS.
C
                  IF ( IPRINT .EQ. 4 .OR. IPRINT .GE. 7 ) THEN
                     IF ( IS .GT. 0 ) THEN
                        K3 = IGPA( IS ) - 1
                        K4 = ISTGP( IG ) - 1
                        L  = ISTGP( IG + 1 ) - K4 - 1
                        IF ( IS .GT. 0 ) THEN
                           IF ( L .EQ. 1 ) THEN
                              WRITE( IOUT, 2090 )
     *                          ANAMES( IS ), 'is ', L, '. '
                           ELSE
                              WRITE( IOUT, 2090 )
     *                          ANAMES( IS ), 'are', L, 's.'
                           END IF
                        END IF
                        IF ( L .GT. 0 )
     *                     WRITE( IOUT, 2100 ) ( GPNAME( K3 + I ),
     *                            GPVALU( K4 + I ), I = 1, L )
                     END IF
                  END IF
C
C  ELEMENT DETAILS.
C
                  IF ( IPRINT .GE. 5 ) THEN
                     DO 400 K  = K1, K2
                        IEL    = IELING( K )
                        IELTYP = ITYPEE( IEL )
                        IF ( IWK( IEL ) .EQ. 0 ) THEN
                           WRITE( IOUT, 2110 ) LNAMES( IEL ),
     *                            IEL, ETYPES( IELTYP ),
     *                            IELV( IELTYP + 1 ) - IELV( IELTYP ),
     *                            IINV( IELTYP + 1 ) - IINV( IELTYP ),
     *                            IEPA( IELTYP + 1 ) - IEPA( IELTYP ),
     *                            ESCALE( K )
                           IF ( IPRINT .LT. 6 ) IWK( IEL ) = 1
                        ELSE
                           WRITE( IOUT, 2120 ) LNAMES( IEL ), IEL,
     *                            ESCALE( K )
                        END IF
C
C  FURTHER ELEMENT DETAILS.
C
                        IF ( IPRINT .GE. 6 ) THEN
                           IF ( IWK( IEL ) .EQ. 0 ) THEN
                              IWK( IEL ) = 1
                              K3 = IELV ( IELTYP ) - 1
                              K4 = ISTEV( IEL ) - 1
                              L  = ISTEV( IEL + 1 ) - K4 - 1
                              WRITE( IOUT, 2130 ) ( ENAMES( K3 + I ),
     *                               VNAMES( IELVAR( K4 + I ) ),
     *                                                    I = 1, L )
                              K3 = IEPA( IELTYP )  - 1
                              K4 = ISTEP( IEL ) - 1
                              L  = ISTEP( IEL + 1 ) - K4 - 1
                              IF ( L .GT. 0 )
     *                          WRITE( IOUT, 2150 ) ( EPNAME( K3 + I ),
     *                                 EPVALU( K4 + I ), I = 1, L )
                           END IF
                        END IF
  400                CONTINUE
                  END IF
C
C  LINEAR ELEMENT DETAILS.
C
                  K5 = ISTADA( IG )
                  K6 = ISTADA( IG + 1 ) - 1
                  IF ( IPRINT .EQ. 6 ) THEN
                     IF ( K5 .LE. K6 ) THEN
                        IF ( K5 .EQ. K6 ) THEN
                           WRITE( IOUT, 2040 ) K6 - K5 + 1, '. '
                           IF ( IPRINT .GE. 6 ) WRITE( IOUT, 2140 ) ' ',
     *                       ( VNAMES( ICNA( K ) ), K = K5, K6 )
                        ELSE
                           WRITE( IOUT, 2040 ) K6 - K5 + 1, 's.'
                           IF ( IPRINT .GE. 6 ) WRITE( IOUT, 2140 ) 's',
     *                       ( VNAMES( ICNA( K ) ), K = K5, K6 )
                        END IF
                     ELSE
                        IF ( ABS( B( IG ) ) .LT. 1.0D-12 ) THEN
                           WRITE( IOUT, 2050 )
                        ELSE
                           WRITE( IOUT, 2160 )
                        END IF
                     END IF
                  ELSE
C
C  FURTHER LINEAR ELEMENT DETAILS.
C
                     IF ( K5 .LE. K6 ) THEN
                        IF ( K5 .EQ. K6 ) THEN
                           WRITE( IOUT, 2040 ) K6 - K5 + 1, '. '
                           IF ( IPRINT .GE. 6 ) WRITE( IOUT, 2200 ) ' ',
     *                       ') ', ( VNAMES( ICNA( K ) ),
     *                         ABYROW( K ), K = K5, K6 )
                        ELSE
                           WRITE( IOUT, 2040 ) K6 - K5 + 1, 's.'
                           IF ( IPRINT .GE. 6 ) WRITE( IOUT, 2200 ) 's',
     *                       's)', ( VNAMES( ICNA( K ) ),
     *                         ABYROW( K ), K = K5, K6 )
                        END IF
                        IF ( ABS( B( IG ) ) .LT. 1.0D-12 ) THEN
                           WRITE( IOUT, 2230 )
                        ELSE
                           WRITE( IOUT, 2220 ) B( IG )
                        END IF
                     ELSE
                        IF ( ABS( B( IG ) ) .LT. 1.0D-12 ) THEN
                           WRITE( IOUT, 2050 )
                        ELSE
                           WRITE( IOUT, 2210 ) B( IG )
                        END IF
                     END IF
                  END IF
  500          CONTINUE
            END IF
         END IF
      END IF
      IF ( IPRINT .GE. 8 ) THEN
         WRITE( IOUT, 2170 )
         DO 510 I = 1, N
C ** Correction 0. 20/12/99: Array holding variable types introduced.
            WRITE( IOUT, 2180 ) I, VNAMES( I ), BL( I ), X( I ),
     *                          BU( I ), VSCALE( I ), 
     *                          VARTYP( ITYPEV( I ) + 1 )
  510    CONTINUE
      END IF
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' Problem name ', A8, /,
     *        /, ' There are ', I8, ' VARIABLES of which ', I8,
     *           ' are artificials',
     *        /, ' There are ', I8, ' GROUPS',
     *        /, ' There are ', I8, ' NONLINEAR ELEMENTS ' )
 2010 FORMAT( /, ' Names of problem variables ',
     *        /, ' ----- -- ------- --------- ',
     *        /, 7( 1X, A10 ) )
 2020 FORMAT( /, ' Names of artificial variables ',
     *        /, ' ----- -- ---------- --------- ',
     *        /, 7( 1X, A10 ) )
C ** Correction -1a. 03/03/00: Integer formats increased
 2030 FORMAT( /, ' Group ', I8, ' is named ', A10, /, '  * It is ',
     *        A12, ' group of type ', A10,
     *        /, '  * The group is scaled by the factor ', 1P, D12.4 )
 2040 FORMAT( '  * The group has a LINEAR ELEMENT with ', I6,
     *        ' variable', A2 )
 2050 FORMAT( '  * The group has no LINEAR ELEMENT. ' )
 2060 FORMAT( '  * The group uses a single NONLINEAR',
     *           ' ELEMENT.  This is element ', A10 )
C ** Correction 0. 20/12/99: Array holding variable types introduced.
 2070 FORMAT( '  * The group uses ', I5, ' NONLINEAR',
     *        ' ELEMENTS. These are elements ', A10,
     *         /, ( 3X, 6( 1X, A10 ) ) )
 2080 FORMAT( '  * The group uses no NONLINEAR ELEMENTS. ' )
C ** Correction -1b. 03/03/00: Integer format increased
 2090 FORMAT( '  * The group-type argument is ', A10, ' and there ', A3,
     *           I8, ' parameter', A2 )
 2100 FORMAT( ( '    * Group parameter ', A10, ' has the value ',
     *          1P, D12.4, '.' ) )
C ** Correction -1c. 03/03/00: Integer format increased
 2110 FORMAT(  '  * Group uses nonlinear element ', A10,
     *         ' number ', I5, ' of type ', A10, /,
     *         '    * No. elemental variables =', I8,
     *         '. No. internal variables =', I8, '.',
     *         /, '    * No. parameter values    =', I8, '.', /,
     *         '    * The element is scaled by the factor ', 1P, D12.4 )
 2120 FORMAT( '  * Group uses nonlinear element ', A10,
     *        ' number ', I5, ' described above. ', /,
     *        '    * The element is scaled by the factor ', 1P, D12.4 )
 2130 FORMAT( ( '    * Elemental variable ', A10, '  is assigned',
     *          ' problem variable ' , A10 ) )
 2140 FORMAT( '    * The linear element uses variable', A1, 1X, 3A10,
     *        /, ( 6X, 6A10 ) )
 2150 FORMAT( ( '    * Elemental parameter ', A10, ' has the value ',
     *          1P, D12.4, '.' ) )
 2160 FORMAT( '  * The group has a constant LINEAR ELEMENT. ' )
C ** Correction 0. 20/12/99: Array holding variable types introduced.
 2170 FORMAT( /, '  #  variable name ',
     *           'lower bound start value upper bound scale factor',
     *        ' type', /, '   -  ------------- ',
     *           '----------- ----------- ----------- ------------',
     *        ' ----' )
 2180 FORMAT( I5, 2X, A10, 1X, 1P, 4D12.4, 3X, A3 )
 2200 FORMAT( '    * The linear element uses variable', A1,
     *        ' (with coefficient', A2, 
     *        /, ( 6X, 3( A10, 1X, 1P, D10.2, 1X ) ) )
 2210 FORMAT( '    * The group has a constant LINEAR ELEMENT with',
     *        ' coefficient ', 1P, D10.2 )
 2220 FORMAT( '    * The constant term for the LINEAR ELEMENT has',
     *        ' coefficient ', 1P, D10.2 )
 2230 FORMAT( '    * There is no constant term for the LINEAR ELEMENT' )
C
C  END OF PRINTP.
C
      END



