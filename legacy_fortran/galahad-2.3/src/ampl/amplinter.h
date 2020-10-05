
/*
 * ==========================================
 * Main header file for the interface between
 * the Galahad solvers and the Ampl language.
 *
 * D. Orban@ECE,            Chicago 2002-2003
 * ==========================================
 */

#include <math.h>

/* AIX does not append underscore to Fortran subroutine names */
#ifdef _AIX
#define FUNDERSCORE(a)   a
#else
#define FUNDERSCORE(a)   a##_
#endif

#ifdef  Fujitsu_frt
#define MAINENTRY MAIN__
#elif  Lahey_lf95
#define MAINENTRY MAIN__
#else
#define MAINENTRY main
#endif

#include "galahad_ampl.h"

#define MAX(m,n)  ((m)>(n)?(m):(n))
#define MIN(m,n)  ((m)<(n)?(m):(n))

/* Some real constants -- GalahadReal is defined in galahad.h */
#define ZERO             (GalahadReal)0.0
#define ONE              (GalahadReal)1.0
#define TWO              (GalahadReal)2.0
#define THREE            (GalahadReal)3.0
#define FIVE             (GalahadReal)5.0
#define FORTRAN_INFINITY (GalahadReal)pow( 10, 20 )

/*
 * =================================
 *  T Y P E   D E F I N I T I O N S
 * =================================
 */

/*
 * Define Fortran types for integer and double precision
 * The following choices are from f2c.h
 */

typedef long int integer;
typedef long int logical;
#define FORTRAN_FALSE (0)     /* Fortran FALSE */
#define FORTRAN_TRUE  (1)     /* Fortran  TRUE */

/*
 * Define shortcuts for the CUTEr library functions,
 * and try to avoid the trailing underscore.
 *
 */

#define USETUP   FUNDERSCORE(usetup)
#define CSETUP   FUNDERSCORE(csetup)

#define UDIMEN   FUNDERSCORE(udimen)
#define UDIMSH   FUNDERSCORE(udimsh)
#define UDIMSE   FUNDERSCORE(udimse)
#define UVARTY   FUNDERSCORE(uvarty)
#define UNAMES   FUNDERSCORE(unames)
#define UREPRT   FUNDERSCORE(ureprt)

#define CDIMEN   FUNDERSCORE(cdimen)
#define CDIMSJ   FUNDERSCORE(cdimsj)
#define CDIMSH   FUNDERSCORE(cdimsh)
#define CDIMSE   FUNDERSCORE(cdimse)
#define CVARTY   FUNDERSCORE(cvarty)
#define CNAMES   FUNDERSCORE(cnames)
#define CREPRT   FUNDERSCORE(creprt)

#define UFN      FUNDERSCORE(ufn)
#define UGR      FUNDERSCORE(ugr)
#define UOFG     FUNDERSCORE(uofg)
#define UBANDH   FUNDERSCORE(ubandh)
#define UDH      FUNDERSCORE(udh)
#define USH      FUNDERSCORE(ush)
#define UEH      FUNDERSCORE(ueh)
#define UGRDH    FUNDERSCORE(ugrdh)
#define UGRSH    FUNDERSCORE(ugrsh)
#define UGREH    FUNDERSCORE(ugreh)
#define UPROD    FUNDERSCORE(uprod)

#define CFN      FUNDERSCORE(cfn)
#define COFG     FUNDERSCORE(cofg)
#define CCFG     FUNDERSCORE(ccfg)
#define CGR      FUNDERSCORE(cgr)
#define CSGR     FUNDERSCORE(csgr)
#define CCFSG    FUNDERSCORE(ccfsg)
#define CCIFG    FUNDERSCORE(ccifg)
#define CCIFSG   FUNDERSCORE(ccifsg)
#define CGRDH    FUNDERSCORE(cgrdh)
#define CDH      FUNDERSCORE(cdh)
#define CSH      FUNDERSCORE(csh)
#define CEH      FUNDERSCORE(ceh)
#define CIDH     FUNDERSCORE(cidh)
#define CISH     FUNDERSCORE(cish)
#define CSGRSH   FUNDERSCORE(csgrsh)
#define CSGREH   FUNDERSCORE(csgreh)
#define CPROD    FUNDERSCORE(cprod)

#define ELFUN    FUNDERSCORE(elfun)
#define RANGE    FUNDERSCORE(range)
#define GROUP    FUNDERSCORE(group)

#define FORTRAN_OPEN  FUNDERSCORE(fortran_open)
#define FORTRAN_CLOSE FUNDERSCORE(fortran_close)

/*
 * Prototypes for CUTEr FORTRAN routines found in libcuter.a
 * See http://cuter.rl.ac.uk/cuter-www/
 */

/* Setup routines */
void USETUP( integer *funit, integer *iout, integer *n, GalahadReal *x, 
	     GalahadReal *bl, GalahadReal *bu, integer *nmax );

void CSETUP( integer *funit, integer *iout, integer *n, integer *m,
	     GalahadReal *x, GalahadReal *bl, GalahadReal *bu, integer *nmax,
	     logical *equatn, logical *linear, GalahadReal *v, GalahadReal *cl,
	     GalahadReal *cu, integer *mmax, logical *efirst,
	     logical *lfirst, logical *nvfrst );

/* Unconstrained dimensioning and report routines */
void UDIMEN( integer funit, integer *n );
void UDIMSH( integer *nnzh );
void UDIMSE( integer *ne, integer *nzh, integer *nzirnh );
void UVARTY( integer *n, integer *ivarty );
void UNAMES( integer *n, char *pname, char *vnames );
void UREPRT( float *calls, float *time );

/* Constrained dimensioning and report routines */
void CDIMEN( integer *funit, integer *n, integer *m );
void CDIMSJ( integer *nnzj );
void CDIMSH( integer *nnzh );
void CDIMSE( integer *ne, integer *nzh, integer *nzirnh );
void CVARTY( integer *n, integer *ivarty );
void CNAMES( integer *n, integer *m, char *pname, char *vnames,
	      char *gnames );
void CREPRT( float *calls, float *time );


/* Unconstrained optimization routines */
void UFN( integer *n, GalahadReal *x, GalahadReal *f );
void UGR( integer *n, GalahadReal *x, GalahadReal *g );
void UOFG( integer *n, GalahadReal *x, GalahadReal *f, GalahadReal *g,
	    logical *grad );

void UBANDH( integer *n, logical *goth, GalahadReal *x, integer *nsemib,
	      GalahadReal *bandh, integer *lbandh );

void UDH( integer *n, GalahadReal *x, integer *lh1, GalahadReal *h );
void USH( integer *n, GalahadReal *x, integer *nnzh, integer *lh,
	   GalahadReal *h, integer *irnh, integer *icnh );
void UEH( integer *n, GalahadReal *x, integer *ne, integer *irnhi,
	   integer *lirnhi, integer *le, integer *iprnhi, GalahadReal *hi,
	   integer *lhi, integer *iprhi, logical *byrows );

void UGRDH( integer *n, GalahadReal *x, GalahadReal *g, integer *lh1,
	     GalahadReal *h);
void UGRSH( integer *n, GalahadReal *x, GalahadReal *g, integer *nnzh,
	     integer *lh, GalahadReal *h, integer *irnh, integer *icnh );
void UGREH( integer *n, GalahadReal *x, GalahadReal *g, integer *ne,
	     integer *irhni, integer *lirnhi, integer *le, integer *iprnhi,
	     GalahadReal *hi, integer *lhi, integer *iprhi, logical *byrows );

void UPROD( integer *n, logical *goth, GalahadReal *x, GalahadReal *p,
	     GalahadReal *q );

/* Constrained optimization routines */
void CFN(  integer *n, integer *m, GalahadReal *x, GalahadReal *f, integer *lc,
	    GalahadReal *c );
void COFG( integer *n, GalahadReal *x, GalahadReal *f, GalahadReal *g,
	    logical *grad );

void CCFG( integer *n, integer *m, GalahadReal *x, integer *lc,
	    GalahadReal *c, logical *jtrans, integer *lcjac1, integer *lcjac2,
	    GalahadReal *cjac, logical *grad );
void CGR(  integer *n, integer *m, GalahadReal *x, logical *grlagf,
	    integer *lv, GalahadReal *v, GalahadReal *g, logical *jtrans,
	    integer *lcjac1, integer *lcjac2, GalahadReal *cjac );
void CSGR( integer *n, integer *m, logical *grlagf, integer *lv,
	    GalahadReal *v, GalahadReal *x, integer *nnzj, integer *lcjac,
	    GalahadReal *cjac, integer *indvar, integer *indfun );

void CCFSG(  integer *n, integer *m, GalahadReal *x, integer *lc,
	      GalahadReal *c, integer *nnzj, integer *lcjac,
	      GalahadReal *cjac, integer *indvar, integer *indfun,
	      logical *grad );
void CCIFG(  integer *n, integer *i, GalahadReal *x, GalahadReal *ci,
	      GalahadReal *gci, logical *grad );
void CCIFSG( integer *n, integer *i, GalahadReal *x, GalahadReal *ci,
	      integer *nnzsgc, integer *lsgci, GalahadReal *sgci,
	      integer *ivsgci, logical *grad );

void CGRDH( integer *n, integer *m, GalahadReal *x, logical *grlagf,
	     integer *lv, GalahadReal *v, GalahadReal *g, logical *jtrans,
	     integer *lcjac1, integer *lcjac2, GalahadReal *cjac,
	     integer *lh1, GalahadReal *h );

void CDH( integer *n, integer *m, GalahadReal *x, integer *lv, GalahadReal *v,
	   integer *lh1, GalahadReal *h );
void CSH( integer *n, integer *m, GalahadReal *x, integer *lv, GalahadReal *v,
	   integer *nnzh, integer *lh, GalahadReal *h, integer *irnh,
	   integer *icnh );
void CEH( integer *n, integer *m, GalahadReal *x, integer *lv, GalahadReal *v,
	   integer *ne, integer *irnhi, integer *lirhni, integer *le,
	   integer *iprnhi, GalahadReal *hi, integer *lhi, integer *iprhi,
	   logical *byrows );

void CIDH( integer *n, GalahadReal *x, integer *iprob, integer *lh1,
	    GalahadReal *h );
void CISH( integer *n, GalahadReal *x, integer *iprob, integer *nnzh,
	    integer *lh, GalahadReal *h, integer *irnh, integer *icnh );

void CSGRSH( integer *n, integer *m, GalahadReal *x, logical *grlagf,
	      integer *lv, GalahadReal *v, integer *nnzj, integer *lcjac,
	      GalahadReal *cjac, integer *indvar, integer *indfun,
	      integer *nnzh, integer *lh, GalahadReal *h, integer *irnh,
	      integer *icnh );
void CSGREH( integer *n, integer *m, GalahadReal *x, logical *grlagf,
	      integer *lv, GalahadReal *v, integer *nnzj, integer *lcjac,
	      GalahadReal *cjac, integer *indvar, integer *indfun,
	      integer *ne, integer *irnhi, integer *lirhni, integer *le,
	      integer *iprnhi, GalahadReal *hi, integer *lhi, integer *iprhi,
	      logical *byrows );

void CPROD( integer *n, integer *m, logical *goth, GalahadReal *x,
	     integer *lv, GalahadReal *v, GalahadReal *p, GalahadReal *q );

/* For backward compatibility with previous versions of CUTE */
#define CSCFG( n, m, x, lc, c, nnzj, lcjac, cjac, indvar, indfun, grad ) \
        CCFSG( n, m, x, lc, c, nnzj, lcjac, cjac, indvar, indfun, grad )
#define CSCIFG( n, i, x, ci, nnzsgc, lsgci, sgci, ivsgci, grad ) \
        CCIFSG( n, i, x, ci, nnzsgc, lsgci, sgci, ivsgci, grad )

/* FORTRAN auxiliary subroutines to retrieve stream unit numbers */
void FORTRAN_OPEN(  integer *funit, char *fname, integer *ierr );
void FORTRAN_CLOSE( integer *funit, integer *ierr );

/*  Low-level SIF functions required by Lancelot-B  */
/* Arrays appear in uppercase, scalars in lowercase */

/*
void ELFUN( GalahadReal *FUVALS, GalahadReal *XVALUE, GalahadReal *EPVALU,
	    integer *ncalcf, integer *ITYPEE, integer *ISTAEV, integer *IELVAR,
	    integer *INTVAR, integer *ISTADH, integer *ISTEPA, integer *ICALCF,
	    integer *ltypee, integer *lstaev, integer *lelvar, integer *lntvar,
	    integer *lstadh, integer *lstepa, integer *lcalcf, integer *lfvalu,
	    integer *lxvalu, integer *lepvlu, integer *ifflag, integer *ifstat );

void RANGE( integer *ielemn, logical *transp, GalahadReal *W1, GalahadReal *W2, integer *nelvar,
	    integer *ninvar, integer *itype, integer *lw1, integer *lw2 );

void GROUP( GalahadReal *GVALUE, integer *lgvalu, GalahadReal *FVALUE, GalahadReal *GPVALU,
	    integer *ncalcg, integer *ITYPEG, integer *ISTGPA, integer *ICALCG,
	    integer *ltypeg, integer *lstgpa, integer *lcalcg, integer *lfvalu,
	    integer *lgpvlu, logical *derivs, integer *igstat );
*/

/* Functions declared in galahad_ps.c */
/*
void ps_initialize( void );
void ranges( fint *IELEM, fint *TRANSP, GalahadReal *W1, GalahadReal *W2,
	     fint *NELV, fint *NINV);
void derivs(fint j, GalahadReal *g, GalahadReal *h);
void funcvals(fint *ICALCF, int n, GalahadReal *FUVALS, GalahadReal *XT);
void funcgrads(fint *ICALCF, int n, GalahadReal *FUVALS,
	       fint *INTVAR, fint *ISTADH);
void groupvals(fint *ICALCG, int n, GalahadReal *ft, GalahadReal *gv);
void groupgrads(fint *ICALCG, int n, GalahadReal *g1, GalahadReal *g2);
*/
