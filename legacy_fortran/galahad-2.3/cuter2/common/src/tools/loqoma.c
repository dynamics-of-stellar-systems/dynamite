
/* ====================================================
 * CUTEr interface for LOQO          October 22nd, 2003
 *
 * Originally written by Andreas Wachter, IBM TJ Watson
 * Cosmetically rearranged by   D. Orban,  Northwestern
 * Specs reading part by   Hande Benson, Drexel College
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#define LOQOMA

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "myalloc.h"
#include "loqo.h"
#include "cuter.h"

typedef struct keyword keyword;
struct keyword {
	char    *name;
	int     type;
	int     ivalue;
	double  dvalue;
};

#define KW(a,b,c,d) {a,b,c,d}

static keyword keywds[] = {
	KW("bndpush", 2, 0, -1.0),
	KW("convex", 0, 0, 0.0),
	KW("dense", 1, -1, 0.0),
	KW("dual", 0, 0, 0.0),
	KW("epsdiag", 2, 0, 1.0e-10),
	KW("epsnum", 2, 0, 0.0),
	KW("epssol", 2, 0, 1.0e-6),
	KW("honor_bnds", 0, 0, 0.0),
	KW("honor_bnds_init", 0, 0, 0.0),
	KW("inftol", 2, 0, 1e-6),
	KW("inftol2", 2, 0, 1e+5),
	KW("iterlim", 1, 500, 0.0),
	KW("lincons", 1, 0, 0.0),
	KW("max", 0, 0, 0.0),
	KW("maximize", 0, 0, 0.0),
	KW("maxit", 1, 200, 0.0),
	KW("min", 0, 0, 0.0),
	KW("mindeg", 0, 0, 0.0),
	KW("minimize", 0, 0, 0.0),
    KW("minlocfil", 0, 0, 0.0),
    KW("mufactor", 2, 0, -1.0),
    KW("noreord", 0, 0, 0.0),
    KW("outlev", 1, 0, 0.0),
    KW("pred_corr", 0, 0, 0.0),
    KW("primal", 0, 0, 0.0),
    KW("quadratic", 0, 0, 0.0),
    KW("sdp", 0, 0, 0.0),
    KW("sigfig", 1, 8, 0.0),
    KW("stablty", 2, 0, 1.0),
    KW("steplen", 2, 0, 0.95),
    KW("timlim", 2, 0, -1.0),
    KW("verbose", 1, 0, 0.0)
};

static int bad_opns = 0;
static int spec = 0;

/* Prototypes */

    double objval(  double *x );
    void   objgrad( double *c, double *x );
    void   hessian( double *Q, double *x, double *y );
    void   conval(  double *h, double *x );
    void   congrad( double *A, double *At, double *x );
    int    compare( const void *A, const void *B );
    void   convert_sparse_format( int ishess, int n, int m, int nnz_in,
				  integer *irow_in, integer *jcol_in,
				  int *nnz_out, int **row_out,
				  int **col_out, int **cute2loqo );
	int rd_opns(char *s);
	int rd_specs(char *Spec_name);
	void *binsearch(char **sp);
	void set_opns(LOQO *lp);

/*
 * Global variables used by auxilliary library functions in ccuter.c
 */

    typedef struct {
	int i, j, pos;
    } CUTErentry;

    integer CUTEr_nvar;        /* number of variables */
    integer CUTEr_ncon;        /* number of constraints */
    integer CUTEr_nnzj;        /* number of nonzeros in Jacobian */
    integer CUTEr_nnzh;        /* number of nonzeros in upper triangular part of Hessian
				  of Lagrangian */

/* Counters for number of function and derivative evaluations */
    int count_f = 0;   /* objective function */
    int count_g = 0;   /* gradient of objective function */
    int count_h = 0;   /* Hessian */
    int count_c = 0;   /* constraints */
    int count_a = 0;   /* constraint Jacobian */

/* Variables used in LOQO main driver */

    int *Acuter2loqo, *Atcuter2loqo, qnz, *Qcuter2loqo;
    double *cscale;      /* some of the constraints must be scaled by minus one, */
    /* if they are inequalities with infinite lower bound.  */

/* ============ */
/* Main program */
/* ============ */

    int main( void ) {
	LOQO *lp;                  /* LOQO data structure */
	char *fname = "OUTSDIF.d"; /* CUTEr data file */
	integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
	integer iout = 6;          /* FORTRAN unit number for error output */
	integer ierr;              /* Exit flag from OPEN and CLOSE */

	integer nconp1;

	integer idummy, *indvar, *indfun, *irnh, *icnh;
	doublereal *x, *bl, *bu, *c, *J, *H;
	doublereal *v = NULL, *cl = NULL, *cu = NULL;
	logical *equatn = NULL, *linear = NULL;
	logical efirst = FALSE_, lfirst = FALSE_, nvfrst = FALSE_, grad;
	char *pname, *vnames, *gnames;
	doublereal f, cmax;
	real calls[7], cpu[2];

	int *iA=NULL, *kA=NULL, *iQ=NULL, *kQ=NULL, *iAt=NULL, *kAt=NULL;

	int i, status;

	/* Open problem description file OUTSDIF.d */
	FORTRAN_OPEN( &funit, fname, &ierr );
	if( ierr ) {
	    printf("Error opening file OUTSDIF.d.\nAborting.\n");
	    exit(1);
	}

	/* Determine problem size */
	CDIMEN( &funit, &CUTEr_nvar, &CUTEr_ncon );

	/* Reserve memory for variables, bounds, and multipliers */
	MALLOC( x,      CUTEr_nvar, doublereal );
	MALLOC( bl,     CUTEr_nvar, doublereal );
	MALLOC( bu,     CUTEr_nvar, doublereal );
	MALLOC( equatn, CUTEr_ncon+1, logical    );
	MALLOC( linear, CUTEr_ncon+1, logical    );
	MALLOC( v,      CUTEr_ncon+1, doublereal );
	MALLOC( cl,     CUTEr_ncon+1, doublereal );
	MALLOC( cu,     CUTEr_ncon+1, doublereal );

	/* Seems to be needed for some Solaris C compilers */
	nconp1 = CUTEr_ncon + 1;

	/* Call initialization routine for CUTEr */
	CSETUP( &funit, &iout, &CUTEr_nvar, &CUTEr_ncon, x, bl, bu,
		&CUTEr_nvar, equatn, linear, v, cl, cu, &nconp1,
		&efirst, &lfirst, &nvfrst );

	/* Free unneeded arrays */
	FREE( equatn );
	FREE( linear );

	if( CUTEr_ncon > 0 )
	{
	    /* Determine number of nonzeros in Jacobian */
	    CDIMSJ( &CUTEr_nnzj );
	    /* CUTEr_nnzj -= CUTEr_nvar; */   /* substract dense gradient of objective function */

	    /* Get Jacobian at starting point */
	    MALLOC( c, CUTEr_ncon, doublereal );
	    MALLOC( J, CUTEr_nnzj, doublereal );
	    MALLOC( indvar, CUTEr_nnzj, integer );
	    MALLOC( indfun, CUTEr_nnzj, integer );
	    grad = TRUE_;
	    /* Here, idummy will be set to nnzj again */
	    CCFSG( &CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon, c, &idummy,
		   &CUTEr_nnzj, J, indvar, indfun, &grad );
	    FREE( c );
	    FREE( J );

#ifdef DEBUG
	    for( i=0; i<CUTEr_nnzj; i++){
		printf("i = %d var = %d fun = %d\n",i,indvar[i],indfun[i]);
	    }
#endif

	    /* Substract nonzero count of dense gradient of Lagrangian */
	    CUTEr_nnzj -= CUTEr_nvar;

	    /* Convert Jacobian sparsity structure */
	    convert_sparse_format( 0, CUTEr_ncon, CUTEr_nvar, CUTEr_nnzj,
				   indfun, indvar, &i, &iA,
				   &kA, &Acuter2loqo );
#ifdef DEBUG
	    for( i=0; i<=CUTEr_nvar; i++ )
		printf("i = %d kA = %d\n",i,kA[i]);
	    for( i=0; i<CUTEr_nnzj; i++ )
		printf("i=%d iA=%d map=%d\n",i,iA[i],Acuter2loqo[i]);
#endif

	    /* LOQO's congrad function also requires the mapping for the 
	       transpose of the Jacobian */
	    convert_sparse_format( 0, CUTEr_nvar, CUTEr_ncon, CUTEr_nnzj,
				   indvar, indfun, &i, &iAt,
				   &kAt, &Atcuter2loqo );
#ifdef DEBUG
	    for( i=0; i<=CUTEr_ncon; i++ )
		printf("i = %d kAt = %d\n",i,kAt[i]);
	    for( i=0; i<CUTEr_nnzj; i++ )
		printf("i=%d iAt=%d map=%d\n",i,iAt[i],Atcuter2loqo[i]);
#endif

	    FREE( iAt );
	    FREE( kAt );
	    FREE( indvar );
	    FREE( indfun );
	}
	else
	{
	    CUTEr_nnzj = 0;
	    CALLOC( kA, CUTEr_nvar+1, int );
	}

	/* Determine number of nonzeros in Hessian of Lagrangian */
	CDIMSH( &CUTEr_nnzh );

	if( CUTEr_nnzh > 0 )
	{
	    /* Get Hessian of Lagrangian at starting point */
	    MALLOC( H, CUTEr_nnzh, doublereal );
	    MALLOC( irnh, CUTEr_nnzh, integer );
	    MALLOC( icnh, CUTEr_nnzh, integer );
	    if( CUTEr_ncon == 0 ) idummy = CUTEr_nnzh; /* for unconstrained problems */
	    /* idummy will be set to nnzh again */
	    CSH( &CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon,
		 v, &idummy, &CUTEr_nnzh, H, irnh, icnh );
	    FREE( H );

#ifdef DEBUG
	    for( i=0; i<CUTEr_nnzh; i++){
		printf("i = %d irnh = %d icnh = %d\n",i,irnh[i],icnh[i]);
	    }
#endif

	    /* Convert Hessian sparsity structure */
	    convert_sparse_format( 1, CUTEr_nvar, CUTEr_nvar, CUTEr_nnzh,
				   irnh, icnh, &qnz, &iQ,
				   &kQ, &Qcuter2loqo );
#ifdef DEBUG
	    for( i=0; i<=CUTEr_nvar; i++ )
		printf("i = %d kQ = %d\n",i,kQ[i]);
	    for( i=0; i<qnz; i++ )
		printf("i=%d iQ=%d map=%d\n",i,iQ[i],Qcuter2loqo[i]);
#endif

	    FREE( irnh );
	    FREE( icnh );
	}
	else
	{
	    CALLOC( kQ, CUTEr_nvar+1, int );
	}

	/* now we can forget the initial multipliers */
	FREE( v );

	/* Now we can start to define the problem for LOQO */
	lp = openlp();

	/* Problem size */
	lp->n   = CUTEr_nvar;
	lp->m   = CUTEr_ncon;

	CALLOC( lp->A, CUTEr_nnzj, double );
	lp->nz = CUTEr_nnzj;
	lp->iA = iA;
	lp->kA = kA;

	/* Hessian */
	CALLOC( lp->Q, qnz, double );
	lp->qnz = qnz;
	lp->iQ  = iQ;
	lp->kQ  = kQ;

	/* Bounds on variables */
	MALLOC( lp->l, CUTEr_nvar, double );
	MALLOC( lp->u, CUTEr_nvar, double );
	for(i=0;i<CUTEr_nvar;i++)
	{
	    if( bl[i] <= -CUTE_INF )
		lp->l[i] = -HUGE_VAL;
	    else
		lp->l[i] = bl[i];

	    if( bu[i] >= CUTE_INF )
		lp->u[i] = HUGE_VAL;
	    else
		lp->u[i] = bu[i];
	}

	/* Bounds for the constraints */
	MALLOC( lp->b,  CUTEr_ncon, double );
	MALLOC( lp->r,  CUTEr_ncon, double );
	MALLOC( cscale, CUTEr_ncon, double );
	for( i=0; i<CUTEr_ncon; i++)
	{
	    if( cl[i] <= -CUTE_INF ) /* Scale the constraint by -1 */
	    {
		cscale[i] = -1.;
		lp->b[i] = -1.*cu[i];
		lp->r[i] = HUGE_VAL;
	    }
	    else
	    {
		cscale[i] = 1.;
		lp->b[i] = cl[i];
		if( cu[i] >= CUTE_INF )
		    lp->r[i] = HUGE_VAL;
		else
		    lp->r[i] = cu[i]-cl[i];
	    }
	}

#ifdef DEBUG
	for(i=0; i<CUTEr_ncon; i++) 
	    printf("cscale[%d] = %f r[%d] = %e b[%d] = %e\n", i, 
		   cscale[i], i, lp->r[i], i, lp->b[i]);
#endif

	/* Space for the objective function gradient */
	MALLOC( lp->c, CUTEr_nvar, double );

	/* Read algorithmic parameters from spec file */
	bad_opns = rd_specs( "LOQO.SPC" );
	if( bad_opns ) {
		printf( "Error in LOQO.SPC, line %d \n", bad_opns );
		return 0;
	}
    set_opns( lp );

	nlsetup( lp );

	/* Assign starting point */
	lp->x = x;

	status = solvelp( lp );

	/* Free memory */
	FREE( Qcuter2loqo );
	FREE( cscale );
	FREE( Atcuter2loqo );
	FREE( Acuter2loqo );

	printf("LOQO status: %d\n",status);

	/* Get problem name */
	MALLOC( pname, FSTRING_LEN+1, char );
	MALLOC( vnames, CUTEr_nvar*FSTRING_LEN, char );
	MALLOC( gnames, CUTEr_ncon*FSTRING_LEN, char );
	CNAMES( &CUTEr_nvar, &CUTEr_ncon, pname, vnames, gnames );

	FREE( vnames );
	FREE( gnames );
	pname[FSTRING_LEN] = '\0';
	i = FSTRING_LEN - 1;
	while( i-- > 0 && pname[i] == ' ') {
	    pname[i] = '\0';
	}

	/* Compute final value of objective function and constraint violation */
	MALLOC( c, CUTEr_ncon, doublereal );
	CFN( &CUTEr_nvar, &CUTEr_ncon, lp->x, &f, &CUTEr_ncon, c );
	cmax = 0.;
	for( i=0; i<CUTEr_ncon; i++)
	{
	    if( cl[i] > -CUTE_INF )
		cmax = max(cmax, cl[i] - c[i]);
	    if( cu[i] < CUTE_INF )
		cmax = max(cmax, c[i] - cu[i]);
	}
	FREE( c );
	for( i=0; i<CUTEr_nvar; i++)
	{
	    if( bl[i] > -CUTE_INF )
		cmax = max(cmax, bl[i] - lp->x[i]);
	    if( bu[i] < CUTE_INF )
		cmax = max(cmax, lp->x[i] - bu[i]);
	}

	/* Output some Loqo info */
	printf( "# Iterations\t%-5d\n", lp->iter );
	printf( "  KKT residual\t%-g\n", max( ABS(lp->pres), ABS(lp->dres) ) );
	printf( "  Infeasibility\t%-g\n", cmax );
	printf( "  LOQO CPU  \t%-10.2f\n", lp->elaptime );
	printf( "# Eval f(x) \t%-6d\n", count_f );
	printf( "# Eval g(x) \t%-6d\n", count_g );
	printf( "# Eval c(x) \t%-6d\n", count_c );
	printf( "# Eval J(x) \t%-6d\n", count_a );
	printf( "# Eval H(x) \t%-6d\n", count_h );
										
	/* Get CUTEr statistics */
	CREPRT( calls, cpu );

	printf("\n\n ************************ CUTEr statistics ************************\n\n");
	printf(" Code used               : LOQO\n");
	printf(" Problem                 : %-s\n", pname);
	printf(" # variables             = %-10d\n", CUTEr_nvar);
	printf(" # constraints           = %-10d\n", CUTEr_ncon);
	printf(" # objective functions   = %-15.7g\n", calls[0]);
	printf(" # objective gradients   = %-15.7g\n", calls[1]);
	printf(" # objective Hessians    = %-15.7g\n", calls[2]);
	printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]);
	printf(" # constraints functions = %-15.7g\n", calls[4]);
	printf(" # constraints gradients = %-15.7g\n", calls[5]);
	printf(" # constraints Hessians  = %-15.7g\n", calls[6]);
	printf(" Exit code               = %-10d\n", status);
	printf(" Final f                 = %-23.15g\n", f);
	printf(" Final ||c||_inf         = %-23.15g\n", cmax);
	printf(" Set up time             = %-10.2f seconds\n", cpu[0]);
	printf(" Solve time              = %-10.2f seconds\n", cpu[1]);
	printf(" count_f (LOQO)          = %-8d\n", count_f);
	printf(" count_g (LOQO)          = %-8d\n", count_g);
	printf(" count_c (LOQO)          = %-8d\n", count_c);
	printf(" count_a (LOQO)          = %-8d\n", count_a);
	printf(" count_h (LOQO)          = %-8d\n", count_h);
	printf(" ******************************************************************\n\n");

	FORTRAN_CLOSE( &funit, &ierr );
	if( ierr ) {
	    printf( "Error closing file %s", fname );
	    return 1;
	}
	return 0;
    }

/* -------------------------------------------------- */
/* Interface-specific functions                       */
/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "objval"
    double objval( double *x ) {
	logical grad = FALSE_;
	doublereal *dummy, f;

	count_f++;

	COFG( &CUTEr_nvar, x, &f, dummy, &grad );
	return f;
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "objgrad"
    void objgrad( double *c, double *x ) {
	logical grad = TRUE_;
	doublereal fdummy;

	count_g++;

	COFG( &CUTEr_nvar, x, &fdummy, c, &grad );
	return;
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "hessian"
    void hessian( double *Q, double *x, double *y ) {
	doublereal *h;
	doublereal *v;
	integer idummy, *irnh, *icnh;
	int i;

	count_h++;

	MALLOC( v, CUTEr_ncon, doublereal );
	MALLOC( h, CUTEr_nnzh, doublereal );
	MALLOC( irnh, CUTEr_nnzh, integer );
	MALLOC( icnh, CUTEr_nnzh, integer );

	/* Rescale the multipliers according to cscale - also, note that LOQO
	   uses '-' in the definition of the Lagrangian while CUTEr uses '+' */
	for( i=0; i<CUTEr_ncon; i++ )
	    v[i] = -1.*cscale[i]*y[i];

	CSH( &CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon,
	     v, &idummy, &CUTEr_nnzh, h, irnh, icnh );

	/* Now copy values into correct places in Q */
	for( i=0; i<qnz; i++ )
	    Q[i] = h[Qcuter2loqo[i]];

	FREE( icnh );
	FREE( irnh );
	FREE( h );
	FREE( v );
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "conval"
    void conval( double *h, double *x ) {
	int i;
	logical jtrans = FALSE_;
	logical grad = FALSE_;
	integer izero=0;
	doublereal *cjac;

	count_c++;

	CCFG( &CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon,
	      h, &jtrans, &izero, &izero, cjac, &grad );

	/* Scale according to cscale */
	for( i=0; i<CUTEr_ncon; i++ )
	    h[i] = cscale[i]*h[i];
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "congrad"
    void congrad( double *A, double *At, double *x ) {
	doublereal *c, *cjac;
	integer idummy, *indvar, *indfun;
	logical grad=TRUE_;

	int i;

	count_a++;

	MALLOC( c, CUTEr_ncon, doublereal );
	MALLOC( cjac, CUTEr_nnzj, doublereal );
	MALLOC( indvar, CUTEr_nnzj, integer );
	MALLOC( indfun, CUTEr_nnzj, integer );

	CCFSG( &CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon,
	       c, &idummy, &CUTEr_nnzj, cjac, indvar, indfun, &grad );

	/* Scale entries according to cscale */
	for( i=0; i<CUTEr_nnzj; i++)
	    cjac[i] = cscale[indfun[i]-1]*cjac[i];

#ifdef DEBUG
	for(i=0;i<CUTEr_nnzj;i++)
	    printf("indfun[%d]= %d  cjac[%d] = %e\n",i,indfun[i],i,cjac[i]);
#endif

	/* Copy entries into LOQO's format */
	for( i=0; i<CUTEr_nnzj; i++) {
	    A[i]  = cjac[Acuter2loqo[i]];
	    At[i] = cjac[Atcuter2loqo[i]];
	}

#ifdef DEBUG
	printf("\nJacobian:\n");
	for(i=0;i<CUTEr_nnzj;i++)
	    printf("%d A = %f  At = %f\n",i,A[i],At[i]);
#endif

	FREE( indfun );
	FREE( indvar );
	FREE( cjac );
	FREE( c );
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "compare"
    int compare(const void *A, const void *B) {
	CUTErentry *a, *b;
	a = (CUTErentry *) A;
	b = (CUTErentry *) B;
	if (a->j < b->j) return(-1);
	if (a->j > b->j) return(1);
	if (a->i < b->i) return(-1);
	if (a->i > b->i) return(1);
	printf("Two entries for same element: i = %d j = %d.  Abort.\n",a->i,a->j);
	exit(1);
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "convert_sparse_format"
    void convert_sparse_format(int ishess,      /* 0 for Jac, 1 for Hessian */
			       int n,           /* number of rows */
			       int m,           /* number of columns */
			       int nnz_in,      /* number of nonzeros in input */
			       integer *irow_in,    /* row position of input*/
			       integer *jcol_in,    /* column position of input */
			       int *nnz_out,    /* number of nonzeros in output */
			       int **i_out,     /* sparse format - rows */
			       int **k_out,     /* sparse format - columns */
			       int **cute2loqo  /* mapping for values */
	) {
	int i, j, icuter;
	CUTErentry *A;

	MALLOC(A, 2*nnz_in, CUTErentry);   /* we are generous... */

	/* copy data to sort array */
	*nnz_out = 0;
	if( ishess ) {
	    for( i=0; i<nnz_in; i++) {
		A[*nnz_out].i   = irow_in[i]-1;  /* Input in Fortran numbering */
		A[*nnz_out].j   = jcol_in[i]-1;
		A[*nnz_out].pos = i;
		(*nnz_out)++;
		if( irow_in[i] != jcol_in[i] )
		{
		    A[*nnz_out].j   = irow_in[i]-1;
		    A[*nnz_out].i   = jcol_in[i]-1;
		    A[*nnz_out].pos = i;
		    (*nnz_out)++;
		}
	    }
	} else {
	    for( i=0; i<nnz_in; i++) {
		A[*nnz_out].i   = irow_in[i]-1;  /* Input in Fortran numbering */
		A[*nnz_out].j   = jcol_in[i]-1;
		A[*nnz_out].pos = i;
		(*nnz_out)++;
	    }
	}

#ifdef DEBUG
	printf("Before qsort:\n\n");
	for(i=0;i<*nnz_out; i++) {
	    printf("i = %d  j = %d pos = %d\n",A[i].i,A[i].j,A[i].pos);
	}
#endif

	/* sort the entries */
	qsort(A, *nnz_out, sizeof(CUTErentry), compare);

#ifdef DEBUG
	printf("After qsort:\n\n");
	for(i=0;i<*nnz_out; i++) {
	    printf("i = %d  j = %d pos = %d\n",A[i].i,A[i].j,A[i].pos);
	}
#endif

	/* Copy sorted elements into output arrays */

	MALLOC(*i_out, *nnz_out, int);
	MALLOC(*k_out, m+1, int);
	MALLOC(*cute2loqo, *nnz_out, int);

	icuter = 0;
	for( j=0; j<m; j++ ) {
	    (*k_out)[j] = icuter;
	    while( icuter < *nnz_out && A[icuter].j == j ) {
		(*i_out)[icuter] = A[icuter].i;
		(*cute2loqo)[icuter] = A[icuter].pos;
		icuter++;
	    }
	}
	(*k_out)[m] = icuter;

	/* That should be it */
	FREE( A );
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "rd_opns"
int rd_opns(char *s) {
	char *s1;
	keyword *kw;

	while (*s) {
		while (*s && *s==' ') s++;
		if (!*s) return 1;
		if ((*s >= 'a' && *s <= 'z') || (*s >= 'A' && *s <= 'Z')) {
			kw = (keyword *)binsearch(&s);
			if (bad_opns) return 0;
			while (*s && *s == ' ') s++;
			switch (kw->type) {
				case 0:
					kw->ivalue = 1; /* depending on keyword, this may need to change */
					break;
				case 1:
					if (*s == '=') {
						s++;
						while (*s && *s == ' ') s++;
						if (*s <= '9' && *s >= '0') {
							kw->ivalue = (int)strtol(s1 = s, &s, 10);
						} else {
							return 0;
						}
					} else {
						return 0;
					}
					break;
				case 2:
                    if (*s == '=') {
                        s++;
                        while (*s && *s == ' ') s++;
                        if (*s <= '9' && *s >= '0') {
                            kw->dvalue = strtod(s1 = s, &s);
                        } else {
                            return 0;
                        }
                    } else {
                        return 0;
                    }
                    break;
			}
		} else {
			return 0;
		}
	}
					
	return 1;
}

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "binsearch"
void *binsearch(char **sp) {
	keyword *kw = keywds;
	keyword *kw1;
	int n = 32; /* number of keywords */
	int n1;
	int c1, c2;
	char *s, *s1, *s2;

	s = *sp;
	while (n > 0) {
		kw1 = kw + (n1 = n >> 1);
		s2 = *(char **)kw1;
		for (s1 = s;; s1++) {
			c1 = tolower(*(unsigned char *)s1);
			if (!(c2 = *s2++)) {
				if (c1 <= ' ' || c1 == '=') {
					*sp = s1;
					return kw1;
				}
				break;
			}
			if (c1 != c2) break;
		}
		if (c1 == '=' || c1 < c2) n = n1;
		else {
			n -= n1 + 1;
			kw = kw1 + 1;
		}
	}
	bad_opns++;
	return 0;
}

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "set_opns"
void set_opns(LOQO *lp) {
	keyword *kw = keywds;
	int verbose, itnlim;
	int noreord, md, mlf;
	int primal, dual;
	int max, maximize, min, minimize;

	lp->bndpush = kw->dvalue; kw++;
	lp->convex = kw->ivalue; kw++;
	lp->dense = kw->ivalue; kw++;
	dual = kw->ivalue; kw++;
	lp->epsdiag = kw->dvalue; kw++;
	lp->epsnum = kw->dvalue; kw++;
	lp->epssol = kw->dvalue; kw++;
	lp->honor_bnds = kw->ivalue; kw++;
	lp->honor_bnds_init = kw->ivalue; kw++; 
	lp->inftol = kw->dvalue; kw++;
	lp->inftol2 = kw->dvalue; kw++;
	itnlim = kw->ivalue; kw++;
	lp->lincons = kw->ivalue; kw++; 
	max = kw->ivalue; kw++; 
	maximize = kw->ivalue; kw++; 
	if (itnlim == 200) lp->itnlim = kw->ivalue; else lp->itnlim = itnlim; kw++; 
	min = kw->ivalue; kw++; 
	md = kw->ivalue; kw++;
	minimize = kw->ivalue; kw++; 
	mlf = kw->ivalue; kw++;
	lp->mufactor = kw->dvalue; kw++;
	noreord = kw->ivalue; kw++;
	verbose = kw->ivalue; kw++; 
	lp->pred_corr = kw->ivalue; kw++;
	primal = kw->ivalue; kw++;
	lp->quadratic = kw->ivalue; kw++; 
	lp->sdp = kw->ivalue; kw++; 
	lp->sf_req = kw->ivalue; kw++;
	lp->stablty = kw->dvalue; kw++; 
	lp->steplen = kw->dvalue; kw++;
	if (kw->dvalue != -1.0) lp->timlim = kw->dvalue; kw++; 
	if (verbose == 0) lp->verbose = kw->ivalue; else lp->verbose = verbose; kw++;
	lp->method=1;	
	if (noreord) lp->method = 0;
	else if (mlf) lp->method = 2;
		else if (md) lp->method = 1;

	lp->pdf=0;
	if (primal) lp->pdf = 1;
	if (dual) lp->pdf = 2;

	lp->max=1;
	if (max || maximize) lp->max = -1;
	if (min || minimize) lp->max = 1;

	if (lp->quadratic) lp->lincons=1;
}
 
/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "rd_specs"
int rd_specs(char *Spec_name) {
	char *s1;
	keyword *kw;
	char	*CUTEr_loc;
	char	*Spec_loc;
	char	specline[80];
	char	option[15], val[15], comment[50];
	char	*optionptr, *valptr;
	int 	i, j;
	FILE	*specfile;
	int	s;

	spec = 1;

	/* open the specs file */
	specfile = fopen(Spec_name, "r");
	fgets(specline, 80, specfile);
	s = sscanf(specline, "%s%s%s", option, val, comment);
	i = 0;
	while (s && !feof(specfile)) {
		i++;
		if ((option[0] >= 'a' && option[0] <= 'z') || (option[0] >= 'A' && option[0] <= 'Z')) {
			optionptr = &option[0];
			kw = (keyword *)binsearch(&optionptr);
			if (bad_opns) { fclose(specfile); return i; }
			switch (kw->type) {
				case 0:
					if (val[0] == 'T') kw->ivalue = 1;
					else if (val[0] == 'F') kw->ivalue = 0;
					else { fclose(specfile); return i; }
					break;
				case 1:
					if ((val[0] <= '9' && val[0] >= '0') || (val[0] == '-' && val[1] <= '9' && val[1] >= 0)) {
						valptr = &val[0];
						kw->ivalue = (int)strtol(s1 = valptr, &valptr, 10);
					} else {
						fclose(specfile);
						return i;
					}
					break;
				case 2:
	                        if (val[0] <= '9' && val[0] >= '0') {
					valptr = &val[0];
      	                      kw->dvalue = strtod(s1 = valptr, &valptr);
            	            } else {
					    fclose(specfile);
                        	    return i;
                        	}
		                  break;
			}
		} else {
			if (option[0] != '*') { /* Comment */
				fclose(specfile);
				return i;
			}
		}
		fgets(specline, 80, specfile);
		s = sscanf(specline, "%s%s%s", option, val, comment);
	}
	fclose(specfile);					
	return 0;
}

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
