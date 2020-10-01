/*
 * ======================================================================
 *
 * cuter.h
 * Data type definitions, constants definitions and function prototypes
 * to interface the CUTEr testing environment Fortran library with C
 *
 * This header file is built from different sources and different authors
 * have contributed to it. In any case, many thanks to Andreas Waechter.
 *
 * D. Orban. This version, July 10 2002.
 * ======================================================================
 */

#include <stdlib.h>

#ifndef CUTER_DOT_H_INCLUDED
#define CUTER_DOT_H_INCLUDED
#endif

/*
 * Define name of main() function on a
 * compiler by compiler basis.
 */
#ifdef Isg95
#define MAINENTRY MAIN_
#endif
#ifdef Ispgf
#define MAINENTRY MAIN_
#endif
#ifdef Isifr
#define MAINENTRY MAIN__
#endif

#ifndef MAINENTRY
#define MAINENTRY main
#endif

/*
 * Define Fortran types for integer and double precision
 * The following choices are from f2c.h
 */
/* typedef long int integer; */
typedef int      integer;
typedef float    real;
typedef double   doublereal;
typedef long int logical;
#define FALSE_ (0)     /* Fortran FALSE */
#define TRUE_  (1)     /* Fortran  TRUE */
#define max(a,b) ((a)>(b)?(a):(b))

#define ZERO     0e0
#define ONE      1e0
#define CUTE_INF 1e20    /* 'infinity' in CUTEr interface */
#define FSTRING_LEN 10   /* Length of Fortran strings     */

/* AIX does not append underscore to Fortran subroutine names */
#ifdef _AIX
#define FUNDERSCORE(a)   a
#else
#define FUNDERSCORE(a)   a##_
#endif

typedef struct VarTypes {
	int nbnds, neq, nlin, nrange, nlower, nupper,
		nineq, nineq_lin, nineq_nlin, neq_lin,
		neq_nlin;
} VarTypes;

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

#define CONNAMES FUNDERSCORE(connames)
#define PBNAME   FUNDERSCORE(pbname)
#define VARNAMES FUNDERSCORE(varnames)

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
#define CSH1     FUNDERSCORE(csh1)
#define CEH      FUNDERSCORE(ceh)
#define CIDH     FUNDERSCORE(cidh)
#define CISH     FUNDERSCORE(cish)
#define CSGRSH   FUNDERSCORE(csgrsh)
#define CSGREH   FUNDERSCORE(csgreh)
#define CPROD    FUNDERSCORE(cprod)
#define CPROD1   FUNDERSCORE(cprod1)
#define CJPROD   FUNDERSCORE(cjprod)

#define FORTRAN_OPEN  FUNDERSCORE(fortran_open)
#define FORTRAN_CLOSE FUNDERSCORE(fortran_close)

/*
 * Prototypes for CUTEr FORTRAN routines found in libcuter.a
 * See http://cuter.rl.ac.uk/cuter-www/
 */

/* Setup routines */
void USETUP( integer *funit, integer *iout, integer *n, doublereal *x,
	      doublereal *bl, doublereal *bu, integer *nmax );
void CSETUP( integer *funit, integer *iout, integer *n, integer *m,
	      doublereal *x, doublereal *bl, doublereal *bu, integer *nmax,
	      logical *equatn, logical *linear, doublereal *v, doublereal *cl,
	      doublereal *cu, integer *mmax, logical *efirst,
	      logical *lfirst, logical *nvfrst );

/* Unconstrained dimensioning and report routines */
void UDIMEN( integer *funit, integer *n );
void UDIMSH( integer *nnzh );
void UDIMSE( integer *ne, integer *nzh, integer *nzirnh );
void UVARTY( integer *n, integer *ivarty );
void UNAMES( integer *n, char *pname, char *vnames );
void UREPRT( real *calls, real *time );

/* Constrained dimensioning and report routines */
void CDIMEN( integer *funit, integer *n, integer *m );
void CDIMSJ( integer *nnzj );
void CDIMSH( integer *nnzh );
void CDIMSE( integer *ne, integer *nzh, integer *nzirnh );
void CVARTY( integer *n, integer *ivarty );
void CNAMES( integer *n, integer *m, char *pname, char *vnames,
             char *gnames );
void CREPRT( real *calls, real *time );

void CONNAMES( integer *m, char *gname );
void PBNAME( integer *n, char *pname );
void VARNAMES( integer *n, char *vname );

/* Unconstrained optimization routines */
void UFN( integer *n, doublereal *x, doublereal *f );
void UGR( integer *n, doublereal *x, doublereal *g );
void UOFG( integer *n, doublereal *x, doublereal *f, doublereal *g,
	    logical *grad );

void UBANDH( integer *n, logical *goth, doublereal *x, integer *nsemib,
	      doublereal *bandh, integer *lbandh );

void UDH( integer *n, doublereal *x, integer *lh1, doublereal *h );
void USH( integer *n, doublereal *x, integer *nnzh, integer *lh,
	   doublereal *h, integer *irnh, integer *icnh );
void UEH( integer *n, doublereal *x, integer *ne, integer *irnhi,
	   integer *lirnhi, integer *le, integer *iprnhi, doublereal *hi,
	   integer *lhi, integer *iprhi, logical *byrows );

void UGRDH( integer *n, doublereal *x, doublereal *g, integer *lh1,
	     doublereal *h);
void UGRSH( integer *n, doublereal *x, doublereal *g, integer *nnzh,
	     integer *lh, doublereal *h, integer *irnh, integer *icnh );
void UGREH( integer *n, doublereal *x, doublereal *g, integer *ne,
	     integer *irhni, integer *lirnhi, integer *le, integer *iprnhi,
	     doublereal *hi, integer *lhi, integer *iprhi, logical *byrows );

void UPROD( integer *n, logical *goth, doublereal *x, doublereal *p,
	     doublereal *q );

/* Constrained optimization routines */
void CFN(  integer *n, integer *m, doublereal *x, doublereal *f, integer *lc,
	    doublereal *c );
void COFG( integer *n, doublereal *x, doublereal *f, doublereal *g,
	    logical *grad );

void CCFG( integer *n, integer *m, doublereal *x, integer *lc,
	    doublereal *c, logical *jtrans, integer *lcjac1, integer *lcjac2,
	    doublereal *cjac, logical *grad );
void CGR(  integer *n, integer *m, doublereal *x, logical *grlagf,
	    integer *lv, doublereal *v, doublereal *g, logical *jtrans,
	    integer *lcjac1, integer *lcjac2, doublereal *cjac );
void CSGR( integer *n, integer *m, logical *grlagf, integer *lv,
	    doublereal *v, doublereal *x, integer *nnzj, integer *lcjac,
	    doublereal *cjac, integer *indvar, integer *indfun );

void CCFSG(  integer *n, integer *m, doublereal *x, integer *lc,
	      doublereal *c, integer *nnzj, integer *lcjac,
	      doublereal *cjac, integer *indvar, integer *indfun,
	      logical *grad );
void CCIFG(  integer *n, integer *i, doublereal *x, doublereal *ci,
	      doublereal *gci, logical *grad );
void CCIFSG( integer *n, integer *i, doublereal *x, doublereal *ci,
	      integer *nnzsgc, integer *lsgci, doublereal *sgci,
	      integer *ivsgci, logical *grad );

void CGRDH( integer *n, integer *m, doublereal *x, logical *grlagf,
	     integer *lv, doublereal *v, doublereal *g, logical *jtrans,
	     integer *lcjac1, integer *lcjac2, doublereal *cjac,
	     integer *lh1, doublereal *h );

void CDH( integer *n, integer *m, doublereal *x, integer *lv, doublereal *v,
	   integer *lh1, doublereal *h );
void CSH( integer *n, integer *m, doublereal *x, integer *lv, doublereal *v,
	   integer *nnzh, integer *lh, doublereal *h, integer *irnh,
	   integer *icnh );
void CSH1( integer *n, integer *m, doublereal *x, integer *lv, doublereal *v,
	   integer *nnzh, integer *lh, doublereal *h, integer *irnh,
	   integer *icnh );
void CEH( integer *n, integer *m, doublereal *x, integer *lv, doublereal *v,
	   integer *ne, integer *irnhi, integer *lirhni, integer *le,
	   integer *iprnhi, doublereal *hi, integer *lhi, integer *iprhi,
	   logical *byrows );

void CIDH( integer *n, doublereal *x, integer *iprob, integer *lh1,
	    doublereal *h );
void CISH( integer *n, doublereal *x, integer *iprob, integer *nnzh,
	    integer *lh, doublereal *h, integer *irnh, integer *icnh );

void CSGRSH( integer *n, integer *m, doublereal *x, logical *grlagf,
	      integer *lv, doublereal *v, integer *nnzj, integer *lcjac,
	      doublereal *cjac, integer *indvar, integer *indfun,
	      integer *nnzh, integer *lh, doublereal *h, integer *irnh,
	      integer *icnh );
void CSGREH( integer *n, integer *m, doublereal *x, logical *grlagf,
	      integer *lv, doublereal *v, integer *nnzj, integer *lcjac,
	      doublereal *cjac, integer *indvar, integer *indfun,
	      integer *ne, integer *irnhi, integer *lirhni, integer *le,
	      integer *iprnhi, doublereal *hi, integer *lhi, integer *iprhi,
	      logical *byrows );

void CPROD( integer *n, integer *m, logical *goth, doublereal *x,
	     integer *lv, doublereal *v, doublereal *p, doublereal *q );

void CPROD1( integer *n, integer *m, logical *goth, doublereal *x,
	     integer *lv, doublereal *v, doublereal *p, doublereal *q );

void CJPROD( integer *n, integer *m, logical *gotj, logical *jtrans,
             doublereal *x, doublereal *p, integer *lp, doublereal *r,
             integer *lr );

/* For backward compatibility with previous versions of CUTE */
#define CSCFG( n, m, x, lc, c, nnzj, lcjac, cjac, indvar, indfun, grad ) CCFSG( n, m, x, lc, c, nnzj, lcjac, cjac, indvar, indfun, grad )
#define CSCIFG( n, i, x, ci, nnzsgc, lsgci, sgci, ivsgci, grad ) CCIFSG( n, i, x, ci, nnzsgc, lsgci, sgci, ivsgci, grad )

/* FORTRAN auxiliary subroutines to retrieve stream unit numbers */
void FORTRAN_OPEN(  integer *funit, char *fname, integer *ierr );
void FORTRAN_CLOSE( integer *funit, integer *ierr );

/*
 * Memory allocation shortcuts
 */

void *CUTEr_malloc( void *object, int length, size_t s );
void *CUTEr_calloc( void *object, int length, size_t s );
void *CUTEr_realloc( void *object, int length, size_t s );
void  CUTEr_free( void **object );

#ifndef MALLOC
#define MALLOC(object,length,type)  object = (type *)CUTEr_malloc(object,length,sizeof(type))
#endif
#ifndef CALLOC
#define CALLOC(object,length,type)  object = (type *)CUTEr_calloc(object,length,sizeof(type))
#endif
#ifndef REALLOC
#define REALLOC(object,length,type) object = (type *)CUTEr_realloc(object,length,sizeof(type))
#endif
#ifndef FREE
#define FREE(object) CUTEr_free((void **)(&(object)))
#endif

