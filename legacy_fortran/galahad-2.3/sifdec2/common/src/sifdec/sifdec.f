C     ( Last modified on 7 Dec 2001 at 09:05:00 )
      PROGRAM       SIFDEC
C
C  This is the main program for running the SIF decoder for the GALAHAD and
C  CUTEr optimization packages. It calls the driver routine SDLANC which does
C  all the work. The purpose of this main program is to open and close all
C  files, and to care for the proper filenames when possible.
C
C  Nick Gould, for CGT Productions.
C  December 7th, 1990.
C
C
      INTEGER       IINGPS, IOUTDA, IINFN , IOUTFN, IOUTRA, IINGR
      INTEGER       IOUTGR, IOUT  , IOUTEX, NEWL  , I
      INTEGER       IINEX , IPRINT, INFORM, IALGOR, IOUTFF, IOUTFD
      INTEGER       IAUTO , IOUTGF, IOUTGD, IOUTEM, IOUTEA, IAD0
      LOGICAL       NONAME, SINGLE
      EXTERNAL      SDLANC
C
C  ASSIGN THE STANDARD OUTPUT UNIT NUMBERS.
C
      PARAMETER   ( IOUT = 6 )
      CHARACTER * 10 PRB
      CHARACTER * 10 PBNAME
      CHARACTER * 24 PRBDAT, PRBFN , PRBRA , PRBOUT, PRBGR , PRBET
      CHARACTER * 24 PRBEX, PRBFF , PRBFD , PRBGF , PRBGD , PRBEA
C
C  ASSIGN THE REMAINING I/O UNIT NUMBERS.
C
      PARAMETER   ( IOUTDA = 55, IINGPS = 61 )
      PARAMETER   ( IOUTFN = 52, IOUTRA = 53, IOUTGR = 54 )
      PARAMETER   ( IOUTEX = 57, IOUTFD = 59, IOUTEA = 66 )
      PARAMETER   ( IOUTGD = 63, IOUTEM = 67 )
      PARAMETER   ( IOUTFF = 0 , IOUTGF = 0 )
C     PARAMETER   ( IOUTFF = 58, IOUTGD = 63 )
      PARAMETER   ( IINFN  = IINGPS,  IINGR = IINGPS, IINEX = IINGPS )
      PARAMETER   ( NONAME = .FALSE. )
C
C  READ PROBLEM'S NAME, BUILD DEFAULT FILE NAMES AND ASSIGN
C  THE ACTUAL VALUES USED.
C
      OPEN(71, FILE='SIFDECODE.CNF', FORM='FORMATTED', STATUS='OLD')
      READ ( 71, 1000 ) PRB
C      WRITE( IOUT, 2000 ) PRB
      DO 10 I = 1, 10
         IF ( PRB( I: I ) .EQ. ' ' ) THEN
           NEWL = I - 1
           GO TO 20
         END IF
   10 CONTINUE
      NEWL = 10
   20 CONTINUE
C
C  specify the method to be used (1=SBMIN,2=AUGLG,3=BARIA).
C
      READ( 71, 1020 ) IALGOR
C
C  specify whether the problem should be described(<0=DEBUG,0=NO,>0=YES)
C
      READ( 71, 1030 ) IPRINT
C
C  read the actual problem name and use it for initial output
C
      READ( 71, 1000 ) PBNAME
      WRITE( IOUT, 2000 ) PBNAME
C
C  specify whether the derivatives are supplied or are to be computed
C  using automatic differentiation
C
      READ( 71, 1020 ) IAUTO
C
C  specify whether AD01 or AD02 should be used to perform the
C  automatic differentiation
C
      READ( 71, 1020 ) IAD0
C
C  Specify the precision of the output files (single=0,double=1)
C
      READ( 71, 1020 ) I
      CLOSE(71)
      SINGLE = I .EQ. 0
C
C  Open output files (if needed)
C
      PRBDAT = '                       '
      PRBOUT = '                       '
      PRBFN  = '                       '
      PRBRA  = '                       '
      PRBGR  = '                       '
      PRBET  = '                       '
      PRBEX  = '                       '
      PRBDAT = PRB( 1 : NEWL ) // '.SIF'
      PRBOUT = 'OUTSDIF.d'
      PRBFN  = 'ELFUN.f'
      PRBFF  = 'ELFUNF.f'
      PRBFD  = 'ELFUND.f'
      PRBRA  = 'RANGE.f'
      PRBGR  = 'GROUP.f'
      PRBGF  = 'GROUPF.f'
      PRBGD  = 'GROUPD.f'
      PRBET  = 'SETTYP.f'
      PRBEX  = 'EXTER.f'
      PRBEA  = 'EXTERA.f'
C
C  OPEN THE RELEVANT FILES - UNIX SYSTEMS.
C
      OPEN ( IINGPS, FILE = PRBDAT, FORM = 'FORMATTED',
     *       STATUS = 'UNKNOWN' )
      REWIND IINGPS
      OPEN ( IOUTDA, FILE = PRBOUT, FORM = 'FORMATTED',
     *       STATUS = 'UNKNOWN' )
      REWIND IOUTDA
      OPEN ( IOUTRA, FILE = PRBRA,  FORM = 'FORMATTED',
     *       STATUS = 'UNKNOWN' )
      REWIND IOUTRA
      OPEN ( IOUTEX, FILE = PRBEX,  FORM = 'FORMATTED',
     *       STATUS = 'UNKNOWN' )
      REWIND IOUTEX
      IF ( IAUTO .EQ. 0 ) THEN
        OPEN ( IOUTFN, FILE = PRBFN,  FORM = 'FORMATTED',
     *         STATUS = 'UNKNOWN' )
        REWIND IOUTFN
        OPEN ( IOUTGR, FILE = PRBGR,  FORM = 'FORMATTED',
     *         STATUS = 'UNKNOWN' )
        REWIND IOUTGR
      ELSE
        IF ( IOUTFF .GT. 0 ) THEN
          OPEN ( IOUTFF, FILE = PRBFF,  FORM = 'FORMATTED',
     *           STATUS = 'UNKNOWN' )
          REWIND IOUTFF
        END IF
        OPEN ( IOUTFD, FILE = PRBFD,  FORM = 'FORMATTED',
     *         STATUS = 'UNKNOWN' )
        REWIND IOUTFD
        IF ( IOUTGF .GT. 0 ) THEN
          OPEN ( IOUTGF, FILE = PRBGF,  FORM = 'FORMATTED',
     *           STATUS = 'UNKNOWN' )
          REWIND IOUTGF
        END IF
        OPEN ( IOUTGD, FILE = PRBGD,  FORM = 'FORMATTED',
     *         STATUS = 'UNKNOWN' )
        REWIND IOUTGD
        OPEN ( IOUTEA, FILE = PRBEA,  FORM = 'FORMATTED',
     *         STATUS = 'UNKNOWN' )
        REWIND IOUTEA
      END IF
      OPEN( UNIT = IOUTEM )
C
C     DO THE WORK
C
      CALL SDLANC( IINGPS, IOUTDA, IINFN , IOUTFN, IOUTFF, IOUTFD,
     *             IOUTRA, IINGR , IOUTGR, IOUTGF, IOUTGD,
     *             IINEX , IOUTEX, IOUTEM, IOUTEA, IPRINT, IOUT  ,
     *             NONAME, IALGOR, IAUTO , IAD0  , SINGLE, INFORM )
C
C   CLOSE THE OPENED FILES
C
      CLOSE( IINGPS )
      IF ( INFORM .EQ. 0 ) THEN
        CLOSE( IOUTDA )
        CLOSE( IOUTRA )
        CLOSE( IOUTEX )
        IF ( IAUTO .EQ. 0 ) THEN
          CLOSE( IOUTFN )
          CLOSE( IOUTGR )
        ELSE
          IF ( IOUTFF .GT. 0 ) CLOSE( IOUTFF )
          CLOSE( IOUTFD )
          IF ( IOUTGF .GT. 0 ) CLOSE( IOUTGF )
          CLOSE( IOUTGD )
          CLOSE( IOUTEA )
        END IF
C
C  IF AN ERROR HAS BEEN DISCOVERED, DELETE THE OUTPUT FILES.
C
      ELSE
        CLOSE( IOUTDA, STATUS = 'DELETE' )
        CLOSE( IOUTRA, STATUS = 'DELETE' )
        IF ( IAUTO .EQ. 0 ) THEN
          CLOSE( IOUTFN, STATUS = 'DELETE' )
          CLOSE( IOUTGR, STATUS = 'DELETE' )
          CLOSE( IOUTEX, STATUS = 'DELETE' )
        ELSE
          IF ( IOUTFF .GT. 0 ) CLOSE( IOUTFF, STATUS = 'DELETE' )
          CLOSE( IOUTFD, STATUS = 'DELETE' )
          IF ( IOUTGF .GT. 0 ) CLOSE( IOUTGF, STATUS = 'DELETE' )
          CLOSE( IOUTGD, STATUS = 'DELETE' )
          CLOSE( IOUTEA, STATUS = 'DELETE' )
        END IF
      END IF
      CLOSE( IOUTEM, STATUS = 'DELETE' )
      STOP
 1000 FORMAT( A10 )
 1020 FORMAT( I2 )
 1030 FORMAT( I6 )
 2000 FORMAT( /, ' Problem name: ', A10 )
      END
C
      SUBROUTINE SDLANC( IINGPS, IOUTDA, IINFN , IOUTFN, IOUTFF, IOUTFD,
     *                   IOUTRA, IINGR , IOUTGR, IOUTGF, IOUTGD,
     *                   IINEX , IOUTEX, IOUTEM, IOUTEA, IPRINT, IOUT  ,
     *                   NONAME, IALGOR, IAUTO , IAD0  , SINGLE, INFORM)
C
C  DECODE A SIF FILE AND CONVERT THE DATA INTO A FORM SUITABLE FOR
C  INPUT TO SBMIN, AUGLG OR BARIA.
C
C  NICK GOULD, FOR CGT PRODUCTIONS.
C  DECEMBER 7TH, 1990.
C
      INTEGER          IINGPS, IINFN , IINGR , INFORM, NMAX  , NGMAX
      INTEGER          NOBMAX, ONLY1 , NIMAX , LIWK  , LWK   , NCONST
      INTEGER          IOUTDA, IOUTRA, IOUTFN, IOUTGR, LENGTH
      INTEGER          I , IG, ISG   , IINEX , IOUTEX, IOUT  , LA, LB
      INTEGER          NSMAX , NBMAX , NETMAX, NOMAX , NLMAX , NELMAX
      INTEGER          IPRINT, NELNUM, NELING, NEGMAX, NEPVMX, NGPVMX
      INTEGER          NINDEX, MAXINS, MAXLEV, MAXARA, NARRAY, NOBJGR
      INTEGER          IALGOR, NGRMAX, NRLNDX, NEPMAX, NGPMAX, NEVMAX
      INTEGER          NINMAX, NUMAX , NSETVC, LSTADG, LSTADA, LELVAR
      INTEGER          LSTAEV, LNTVAR, LBNDS , LINTRE, LICNA , NREAL
      INTEGER          NLINOB, NNLNOB, NLINEQ, NNLNEQ, NLININ, NNLNIN
      INTEGER          NFREE , NFIXED, NLOWER, NUPPER, NBOTH , NSLACK
      INTEGER          N , NG, NBND  , NELTYP, NLVARS, NOBJ  , NRANGE
      INTEGER          NNZA  , NGRTYP, NSTART, NLISGP, NNLVRS, NOBBND
      INTEGER          IOUTFF, IOUTFD, IOUTGF, IOUTGD, IOUTEM, IAUTO
      INTEGER          IOUTEA, IAD0
      DOUBLE PRECISION BIG   , BLO   , BUP
      LOGICAL          DEBUG , NONAME, SINGLE, ONEOBJ, GOTLIN
      CHARACTER * 10   NAMEOF, NAMERH, NAMERA, NAMEBN, NAMEST, NAMEOB
      CHARACTER * 72   LINEEX
      PARAMETER      ( BIG  = 1.0D+20 )
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
C  ---------------------------------------------------------------------
C
C
C  End of parameters which might be changed by the user.
C
C  ---------------------------------------------------------------------
C
C  DEPENDENCIES ON THE MAXIMUM NUMBER OF NONTRIVIAL GROUP TYPES.
C  NGPMAX IS THE TOTAL NUMBER OF GROUP PARAMETERS.
C
      PARAMETER      ( NGPMAX = NGRMAX )
C
C  DEPENDENCIES ON THE MAXIMUM NUMBER OF GROUPS
C
      PARAMETER      ( NOMAX  = NGMAX )
      PARAMETER      ( NUMAX  = NGMAX )
      PARAMETER      ( LSTADG = NGMAX )
      PARAMETER      ( LSTADA = NGMAX )
      PARAMETER      ( LB     = NGMAX )
      PARAMETER      ( LBNDS  = NMAX + NGMAX )
C
C  DEPENDENCIES ON THE MAXIMUM TOTAL NUMBER OF REAL
C  PARAMETERS ASSOCIATED WITH GROUPS.
C
      PARAMETER      ( LWK    = NGPVMX )
C
C  DEPENDENCIES ON THE MAXIMUM NUMBER OF NONLINEAR ELEMENT TYPES.
C  NETMAX, NIMAX AND NEPMAX ARE THE TOTAL NUMBER OF ELEMENTAL AND
C  INTERNAL VARIABLES AND PARAMETERS RESPECTIVELY.
C
      PARAMETER      ( NETMAX = 5 * NLMAX )
      PARAMETER      ( NIMAX  = 5 * NLMAX )
      PARAMETER      ( NEPMAX = 3 * NLMAX )
C
C  DEPENDENCIES ON THE MAXIMUM NUMBER OF NONLINEAR ELEMENTS.
C
      PARAMETER      ( NEGMAX = NELMAX )
      PARAMETER      ( LSTAEV = NELMAX )
      PARAMETER      ( LNTVAR = NELMAX + 1 )
      PARAMETER      ( LINTRE = NELMAX )
      PARAMETER      ( LIWK   = NELMAX + NGMAX )
C
C  DEPENDENCIES ON THE MAXIMUM TOTAL NUMBER OF ELEMENTAL VARIABLES.
C
      PARAMETER      ( LELVAR = NEVMAX )
C
C  DEPENDENCIES ON THE MAXIMUM NUMBER OF NONZEROS IN LINEAR ELEMENTS.
C
      PARAMETER      ( LICNA  = LA     )
C
C  MAXIMUM NUMBER OF STATEMENTS IN A DO-LOOP.
C
      PARAMETER      ( MAXINS = 200    )
C
C  MAXIMUM NESTING OF DO-LOOPS
C
      PARAMETER      ( MAXLEV = 3      )
C
C  MAXIMUM NUMBER OF ARRAY INSTRUCTIONS.
C
      PARAMETER      ( MAXARA = 150    )
C
C  MAXIMUM SIZE OF DICTIONARY.
C
      PARAMETER      ( LENGTH = NMAX + NGMAX + NELMAX + NINMAX + 1000 )
C
C  ARRAY DEFINITIONS.
C
      INTEGER          INLIST( LENGTH ), ISTAEV( NELMAX )
      INTEGER          ISTATE( NGMAX ), ITABLE ( LENGTH )
      INTEGER          IELV  ( NLMAX  ), IINV  ( NLMAX )
      INTEGER          ITYPEE( NELMAX ), IELING( NEGMAX, 2 )
      INTEGER          IEPA  ( NLMAX  ), IGPA  ( NGRMAX )
      INTEGER          IDROWS( 2, NGMAX ), ITYPEG( NGMAX )
      INTEGER          IELVAR( LELVAR ), INTVAR( LNTVAR )
      INTEGER          ISTADA( LSTADA ), ICNA  ( LICNA  )
      INTEGER          ISTEP ( NELMAX ), ISTGP( NGMAX ), IWK( LIWK )
      INTEGER          INDVAL( NINDEX ), INSTR( 5, MAXINS, MAXLEV )
      INTEGER          NINSTR( MAXLEV ), IARRAY( 5, 3, MAXARA )
      INTEGER          ICOORD( LA, 2 ), ISTADG( NGMAX ), IJUMP( NLMAX )
      INTEGER          ITYPEV( NMAX )
      DOUBLE PRECISION GPTEMP( NGPVMX )
      DOUBLE PRECISION EPVALU( NEPVMX ), GPVALU( NGPVMX ), DFAULT( NMAX)
      DOUBLE PRECISION A( LA ), BND( 2, NMAX, NBMAX ), REALVL( NRLNDX )
      DOUBLE PRECISION BNDFLT( 2, NBMAX ), CSTART( NGMAX, NSMAX )
      DOUBLE PRECISION RSCALE( NGMAX ), CSCALE( NMAX ), WK( LWK )
      DOUBLE PRECISION RDROWS( 2, NGMAX ), VSTART( NMAX, NSMAX )
      DOUBLE PRECISION RVALUE( MAXARA, 3 ), VARRAY( 2, MAXARA )
      DOUBLE PRECISION FBOUND( 2, NOBMAX ), WEIGHT( NEGMAX )
      DOUBLE PRECISION ABYROW( LA ), B( LB ), BL( LBNDS ), BU( LBNDS )
      DOUBLE PRECISION X( NMAX ), U( NUMAX ), ESCALE( NEGMAX )
      DOUBLE PRECISION VSCALE( NMAX ), GSCALE( NGMAX ), CLMULT( NGMAX )
      LOGICAL          INTREP( LINTRE ), LDEFND( NLMAX )
      LOGICAL          SETVEC( NSETVC ), GXEQX( NGMAX )
      CHARACTER * 1    S( 2 )
      CHARACTER * 2    FARRAY( MAXARA )
      CHARACTER * 4    ARE( 2 )
      CHARACTER * 8    PNAME
      CHARACTER * 10   NAMIIN( NINDEX ), NAMRIN( NRLNDX )
      CHARACTER * 10   LONAME( NINMAX ), BNAMES( NBMAX  )
      CHARACTER * 10   GNAMES( NGMAX  ), VNAMES( NMAX   )
      CHARACTER * 10   ETYPES( NLMAX  ), GTYPES( NGRMAX )
      CHARACTER * 10   LNAMES( NELMAX ), OBNAME( NOBMAX )
      CHARACTER * 10   EPNAME( NEPMAX ), GPNAME( NGPMAX )
      CHARACTER * 10   ONAMES( NOMAX  ), ENAMES( NETMAX )
      CHARACTER * 10   EXNAME( NINMAX ), SNAMES( NSMAX )
      CHARACTER * 10   ANAMES( NGRMAX ), INAMES( NIMAX  )
      CHARACTER * 10   MINAME( NINMAX ), RENAME( NINMAX )
      CHARACTER * 10   INNAME( NINMAX )
      CHARACTER * 10   ARRAY( 3, MAXARA ), CARRAY( 2, MAXARA )
      CHARACTER * 12   KEY   ( LENGTH )
      CHARACTER * 160  NULINE
      EXTERNAL         GPSMPS, INLANC, PRINTP, MAKEFN, MAKEGR, ONLY1
      DATA S / ' ', 's' /, ARE / ' is ', 'are ' /
      DATA ONEOBJ / .FALSE. /
      DEBUG  = IPRINT .LT. 0
      IF ( SINGLE ) THEN
         WRITE( IOUT, 2050 )
      ELSE
         WRITE( IOUT, 2060 )
      END IF
      IF ( IPRINT .NE. 0 ) IPRINT = 9
C
C  READ THE GPS MPS DATA.
C
      CALL       GPSMPS( LA    , NMAX  , NGMAX , NOMAX , NLMAX , NELMAX,
     *                   NIMAX , NETMAX, NEVMAX, NGRMAX, NSMAX , NEPMAX,
     *                   NGPMAX, NBMAX , NOBMAX, NNZA  , LENGTH, N , NG,
     *                   NOBJ  , NCONST, NRANGE, NBND  , NSTART, NELTYP,
     *                   NGRTYP, NLVARS, NNLVRS, NLISGP, LIWK  ,
     *                   NELNUM, NELING, NARRAY, NINDEX, NEGMAX, NEPVMX,
     *                   NGPVMX, MAXINS, MAXLEV, MAXARA, NRLNDX, NOBBND,
     *                   PNAME , ICOORD, IELING, INLIST, ITABLE, ISTATE,
     *                   IELV  , IINV  , ITYPEE, IDROWS, IELVAR, ISTADG,
     *                   ITYPEG, IEPA  , IGPA  , IWK   , ISTEP , ISTAEV,
     *                   ISTGP , INDVAL, INSTR , NINSTR, IARRAY, ITYPEV,
     *                   A, BND, VSTART, CSTART, RSCALE, CSCALE, RDROWS,
     *                   REALVL, DFAULT, RVALUE, VARRAY, EPVALU, BNDFLT,
     *                   GPVALU, GPTEMP, FARRAY, FBOUND, WEIGHT, NAMIIN,
     *                   NAMRIN, GNAMES, VNAMES, BNAMES, ETYPES, INAMES,
     *                   LNAMES, ONAMES, ENAMES, SNAMES, ANAMES, GTYPES,
     *                   EPNAME, GPNAME, OBNAME, ARRAY , CARRAY, KEY   ,
     *                   SINGLE, IINGPS, IOUT  , INFORM, DEBUG )
      IF ( INFORM .NE. 0 ) THEN
         WRITE( IOUT, 2010 ) INFORM
         RETURN
      END IF
C
C  ASSIGN THE GROUPS TO CONSTRAINT TYPES AND OBJECTIVES.
C
      NLINOB = 0
      NNLNOB = 0
      NLINEQ = 0
      NNLNEQ = 0
      NLININ = 0
      NNLNIN = 0
      DO 100 IG = 1, NG
         ISG    = ISTATE( IG )
         IF ( ISG .GT. 0 ) THEN
            ISG = MOD( ISG - 1, 4 )
            IF ( ISG .EQ. 0 ) NLINOB = NLINOB + 1
            IF ( ISG .EQ. 1 ) NLINEQ = NLINEQ + 1
            IF ( ISG .GE. 2 ) NLININ = NLININ + 1
         ELSE
            ISG = MOD( ISG + 1, 4 )
            IF ( ISG .EQ.   0 ) NNLNOB = NNLNOB + 1
            IF ( ISG .EQ. - 1 ) NNLNEQ = NNLNEQ + 1
            IF ( ISG .LE. - 2 ) NNLNIN = NNLNIN + 1
         END IF
  100 CONTINUE
C
C  SELECT RHS, RANGES AND BOUNDS.
C
      IF ( NCONST .GT. 0 ) NAMERH = VNAMES( NLVARS + 1 )
      IF ( NRANGE .GT. 0 ) NAMERA = VNAMES( NLVARS + NCONST + 1 )
      IF ( NBND   .GT. 0 ) NAMEBN = BNAMES( 1 )
      IF ( NSTART .GT. 0 ) NAMEST = SNAMES( 1 )
      IF ( NOBJ   .GT. 0 .AND. ONEOBJ ) NAMEOF = ONAMES( 1 )
      IF ( NOBBND .GT. 0 ) NAMEOB = OBNAME( 1 )
      IF ( IPRINT .NE. 0 ) WRITE( IOUT, 2070 ) NCONST, NRANGE, NBND,
     *     NSTART, NOBJ, NOBBND
C
C  CONVERT TO INPUT FOR ONE OF THE LANCELOT PROGRAMS.
C
      CALL       INLANC( N     , NLVARS, NG    , NELNUM, NOBJ  , LENGTH,
     *                   LSTADG, LELVAR, LSTAEV, LNTVAR, LICNA , LSTADA,
     *                   LA, LB, LBNDS , LINTRE, LIWK  , LWK   , NMAX  ,
     *                   NGMAX , NBMAX , NSMAX , NLMAX , NELMAX, NEGMAX,
     *                   NOBMAX, NGRMAX, NGPVMX, NEPVMX, NOMAX , NLISGP,
     *                   NBND  , NNZA  , NCONST, NSTART, NRANGE, NOBJGR,
     *                   NOBBND, NELTYP, NGRTYP, PNAME , ONEOBJ,
     *                   NAMEOB, NAMERH, NAMERA, NAMEBN, NAMEST, NAMEOF,
     *                   ISTADG, IELVAR, ISTAEV, INTVAR, ICNA  , ISTADA,
     *                   ICOORD, INLIST, ITABLE, ISTATE,
     *                   IDROWS, IELV  , IINV  , IGPA  , IELING( 1, 1 ),
     *                   ISTEP , ISTGP , ITYPEE, ITYPEG, ITYPEV, IWK   ,
     *                   A, BND, VSTART, CSTART, RSCALE, CSCALE,
     *                   RDROWS, DFAULT, WEIGHT, BNDFLT, WK    ,
     *                   GPVALU, EPVALU, FBOUND, ABYROW, B , BL, BU, X ,
     *                   CLMULT, ESCALE, GSCALE, VSCALE, INTREP, GXEQX,
     *                   KEY   , GNAMES, VNAMES, BNAMES, SNAMES, ONAMES,
     *                   ETYPES, GTYPES, OBNAME, IALGOR, IAUTO,
     *                   IOUT  , IOUTDA, SINGLE, INFORM, DEBUG )
      IF ( INFORM .NE. 0 ) THEN
         WRITE( IOUT, 2020 ) INFORM
         RETURN
      END IF

C
C  ASSIGN THE VARIABLES TO BOUND TYPES.
C
      NFREE  = 0
      NFIXED = 0
      NLOWER = 0
      NUPPER = 0
      NBOTH  = 0
      IF ( IALGOR .LE. 2 ) THEN
         NSLACK = NLININ + NNLNIN
      ELSE
         NSLACK = 0
      END IF
      NREAL  = N - NSLACK
      DO 110 I = 1, NREAL
         BLO   = BL( I )
         BUP   = BU( I )
         IF ( BLO .LE. - BIG .AND. BUP .GE. BIG ) NFREE  = NFREE  + 1
         IF ( BLO .LE. - BIG .AND. BUP .LT. BIG ) NUPPER = NUPPER + 1
         IF ( BLO .GT. - BIG .AND. BUP .GE. BIG ) NLOWER = NLOWER + 1
         IF ( BLO .GT. - BIG .AND. BUP .LT. BIG ) THEN
            IF ( BLO .EQ. BUP ) THEN
                NFIXED = NFIXED + 1
            ELSE
                NBOTH  = NBOTH  + 1
            END IF
         END IF
  110 CONTINUE
C
C  PRINT PROBLEM SUMMARY.
C
      IF ( NLINOB .GT. 0 ) WRITE( IOUT, 2100 )
     *          NLINOB, S( ONLY1( NLINOB ) )
      IF ( NNLNOB .GT. 0 ) WRITE( IOUT, 2110 )
     *          NNLNOB, S( ONLY1( NNLNOB ) )
      IF ( NLINEQ + NLININ + NNLNEQ + NNLNIN .GT. 0 ) WRITE( IOUT, 2000)
      IF ( NLINEQ .GT. 0 ) WRITE( IOUT, 2120 ) ARE( ONLY1( NLINEQ ) ),
     *          NLINEQ, S( ONLY1( NLINEQ ) )
      IF ( NLININ .GT. 0 ) WRITE( IOUT, 2130 ) ARE( ONLY1( NLININ ) ),
     *          NLININ, S( ONLY1( NLININ ) )
      IF ( NNLNEQ .GT. 0 ) WRITE( IOUT, 2140 ) ARE( ONLY1( NNLNEQ ) ),
     *          NNLNEQ, S( ONLY1( NNLNEQ ) )
      IF ( NNLNIN .GT. 0 ) WRITE( IOUT, 2150 ) ARE( ONLY1( NNLNIN ) ),
     *          NNLNIN, S( ONLY1( NNLNIN ) )
      WRITE( IOUT, 2000 )
      IF ( NFREE  .GT. 0 ) WRITE( IOUT, 2200 ) ARE( ONLY1( NFREE  ) ),
     *          NFREE , S( ONLY1( NFREE  ) )
      IF ( NUPPER .GT. 0 ) WRITE( IOUT, 2210 ) ARE( ONLY1( NUPPER ) ),
     *          NUPPER, S( ONLY1( NUPPER ) )
      IF ( NLOWER .GT. 0 ) WRITE( IOUT, 2220 ) ARE( ONLY1( NLOWER ) ),
     *          NLOWER, S( ONLY1( NLOWER ) )
      IF ( NBOTH  .GT. 0 ) WRITE( IOUT, 2230 ) ARE( ONLY1( NBOTH  ) ),
     *          NBOTH,  S( ONLY1( NBOTH  ) )
      IF ( NFIXED .GT. 0 ) WRITE( IOUT, 2240 ) ARE( ONLY1( NFIXED ) ),
     *          NFIXED, S( ONLY1( NFIXED ) )
      IF ( NSLACK .GT. 0 ) WRITE( IOUT, 2250 ) ARE( ONLY1( NSLACK ) ),
     *          NSLACK, S( ONLY1( NSLACK ) )
      WRITE( IOUTDA, 2080 ) PNAME,
     *          NFREE , NFIXED, NLOWER, NUPPER, NBOTH , NSLACK,
     *          NLINOB, NNLNOB, NLINEQ, NNLNEQ, NLININ, NNLNIN
C
C  PRINT DETAILS OF THE PROBLEM.
C
      CALL       PRINTP( NMAX, NGMAX, NLMAX,
     *                   NELMAX, NETMAX,
     *                   NEVMAX, NEPMAX, NGRMAX, NEGMAX, NEPVMX,
     *                   NGPVMX, NGPMAX, LSTADA, LICNA, LIWK,
     *                   N, NG, NLVARS, NELNUM,
     *                   ISTATE, ISTADG, IELVAR, ITYPEG, ITYPEE,
     *                   IELV, IINV, IEPA, IGPA,
     *                   ISTADA, ICNA, ISTGP, ISTEP, ISTAEV,
     *                   IELING, ITYPEV, IWK, ABYROW, B, BL, BU, X,
     *                   EPVALU, GPVALU, GSCALE, ESCALE, VSCALE,
     *                   PNAME, VNAMES, GNAMES,
     *                   LNAMES, ETYPES, ENAMES,
     *                   ANAMES, EPNAME, GPNAME, GTYPES,
     *                   IOUT, IPRINT )
      IF ( NONAME ) PNAME = '        '
C
C  MAKE SUBROUTINES ELFUN AND RANGE.
C
      IF ( IAUTO .EQ. 0 ) THEN
         CALL MAKEFN( IINFN , IOUT  , IOUTFN, IOUTRA, INFORM,
     *                NLMAX , NIMAX , NETMAX, NINMAX, NUMAX ,
     *                NELNUM, NELTYP, PNAME , ENAMES, INAMES, RENAME,
     *                INNAME, LONAME, MINAME, EXNAME, ETYPES, LDEFND,
     *                LENGTH, ITABLE, KEY   , IELV  , IINV  , INLIST,
     *                EPNAME, IEPA  , NEPMAX, DEBUG , IJUMP ,
     *                U     , SETVEC, NSETVC, SINGLE, NULINE, GOTLIN,
     *                IPRINT )
         IF ( INFORM .NE. 0 ) THEN
            WRITE( IOUT, 2030 ) INFORM
            RETURN
         END IF
C
C  MAKE SUBROUTINES ELFUNF, ELFUND AND RANGE
C
      ELSE
         CALL MAFNAD( IINFN , IOUT  , IOUTFF, IOUTFD, IOUTRA,
     *                IOUTEM, INFORM, NLMAX , NIMAX , NETMAX,
     *                NINMAX, NUMAX , NELNUM, NELTYP, PNAME , ENAMES,
     *                INAMES, RENAME, INNAME, LONAME, MINAME, EXNAME,
     *                ETYPES, LDEFND, LENGTH, ITABLE, KEY   , IELV  ,
     *                IINV  , INLIST, EPNAME, IEPA  , NEPMAX, DEBUG ,
     *                IJUMP , U     , SETVEC, NSETVC, SINGLE,
     *                NULINE, GOTLIN, IAUTO , IAD0  , IPRINT )
         IF ( INFORM .NE. 0 ) THEN
            WRITE( IOUT, 2090 ) INFORM
            RETURN
         END IF
      END IF
C
C  MAKE SUBROUTINE GROUP AND OBTAIN GROUP INFORMATION.
C
      IF ( IAUTO .EQ. 0 ) THEN
         CALL MAKEGR( IINGR , IOUT  , IOUTGR, INFORM, NGRTYP,
     *                NGRMAX, NLMAX , NINMAX, PNAME , ANAMES,
     *                RENAME, INNAME, LONAME, MINAME, EXNAME, GTYPES,
     *                LDEFND, GPNAME, IGPA  , NGPMAX, DEBUG , LENGTH,
     *                ITABLE, KEY   , INLIST, SINGLE, NULINE, GOTLIN,
     *                IPRINT )
         IF ( INFORM .NE. 0 ) THEN
            WRITE( IOUT, 2040 ) INFORM
            RETURN
         END IF
C
C  MAKE SUBROUTINES GROUPF AND GROUPD
C
      ELSE
         CALL MAGRAD( IINGR , IOUT  , IOUTGF, IOUTGD, IOUTEM, INFORM,
     *                NGRTYP, NGRMAX, NLMAX , NINMAX,
     *                PNAME , ANAMES, RENAME, INNAME, LONAME, MINAME,
     *                EXNAME, GTYPES, LDEFND, GPNAME, IGPA  , NGPMAX,
     *                DEBUG , LENGTH, ITABLE, KEY   , INLIST, SINGLE,
     *                NULINE, GOTLIN, IAUTO , IAD0  , IPRINT )
         IF ( INFORM .NE. 0 ) THEN
            WRITE( IOUT, 2160 ) INFORM
            RETURN
         END IF
      END IF
C
C  FINALLY, READ ANY ADDITIONAL PROGRAMS.
C
  500 CONTINUE
      IF ( GOTLIN ) THEN
         LINEEX( 1: 72 ) = NULINE( 1: 72 )
         GOTLIN = .FALSE.
      ELSE
         READ( UNIT = IINEX, FMT = 1000, END = 600, ERR = 600 ) LINEEX
      END IF
C
C  SKIP BLANK LINES.
C
      DO 510 I = 1, 72
         IF ( LINEEX( I: I ) .NE. ' ' ) THEN
            WRITE( IOUTEX, 1000 ) LINEEX
            GO TO 500
         END IF
  510 CONTINUE
      GO TO 500
  600 CONTINUE
C
C  IF REQUIRED, TRANSLATE ANY EXTERNAL FILE TO ACCEPT AUTOMATIC
C  DIFFERENTIATION CONSTRUCTS
C
      IF ( IAUTO .EQ. 1 .OR. IAUTO .EQ. 2 )
     *     CALL TRANS( IOUT, IOUTEX, IOUTEA, IOUTEM, SINGLE, IAUTO,
     *                 IAD0, NAMRIN, NRLNDX, NAMIIN, NINDEX )
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 1000 FORMAT( A72 )
 2000 FORMAT( ' ' )
 2010 FORMAT( /, ' Return from GPSMPS, INFORM = ', I3 )
 2020 FORMAT( /, ' Return from INLANC, INFORM = ', I3 )
 2030 FORMAT( /, ' Return from MAKEFN, INFORM = ', I3 )
 2040 FORMAT( /, ' Return from MAKEGR, INFORM = ', I3 )
 2050 FORMAT( /, ' Single precision version will be formed. ', / )
 2060 FORMAT( /, ' Double precision version will be formed. ', / )
 2070 FORMAT( /, '  NCONST  NRANGE    NBND  NSTART    NOBJ  NOBBND ',
     *        /, 6I8, / )
 2080 FORMAT( A8, 12I8 )
 2090 FORMAT( /, ' Return from MAFNAD, INFORM = ', I3 )
 2100 FORMAT( ' The objective function uses ', I8, ' linear group', A1 )
 2110 FORMAT( ' The objective function uses ', I8,
     *        ' nonlinear group', A1 )
 2120 FORMAT( ' There ', A4, I8, ' linear equality constraint', A1 )
 2130 FORMAT( ' There ', A4, I8, ' linear inequality constraint', A1 )
 2140 FORMAT( ' There ', A4, I8, ' nonlinear equality constraint', A1 )
 2150 FORMAT( ' There ', A4, I8,
     *          ' nonlinear inequality constraint', A1 )
 2160 FORMAT( /, ' Return from MAGRAD, INFORM = ', I3 )
 2200 FORMAT( ' There ', A4, I8, ' free variable', A1 )
 2210 FORMAT( ' There ', A4, I8, ' variable', A1,
     *          ' bounded only from above ' )
 2220 FORMAT( ' There ', A4, I8, ' variable', A1,
     *          ' bounded only from below ' )
 2230 FORMAT( ' There ', A4, I8,
     *        ' variable', A1, ' bounded from below and above ' )
 2240 FORMAT( ' There ', A4, I8, ' fixed variable', A1 )
 2250 FORMAT( ' There ', A4, I8, ' slack variable', A1 )
C
C  END OF SDLANC.
C
      END

C
C END OF PROGRAM SIFDEC
C
