
/* ====================================================
 * CUTEr interface for generic package     Feb. 3, 2003
 *
 * D. Orban
 *
 * Take a look at $CUTER/common/include/cuter.h     and
 * $CUTER/common/src/tools/loqoma.c  for more examples.
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GENCMA

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cuter.h"

#define GENC    genc
#define GENSPC  genspc
#define GETINFO getinfo

#ifdef Isg95
#define MAINENTRY MAIN_
#else
#define MAINENTRY main
#endif

    doublereal GENC( doublereal );
    void GENSPC( integer, char* );
    void GETINFO( integer, integer, doublereal*, doublereal*,
                  doublereal*, doublereal*, logical*, logical*,
                  VarTypes* );
                 

    integer CUTEr_nvar;        /* number of variables */
    integer CUTEr_ncon;        /* number of constraints */
    integer CUTEr_nnzj;        /* number of nonzeros in Jacobian */
    integer CUTEr_nnzh;        /* number of nonzeros in upper triangular
                                  part of the Hessian of the Lagrangian */

    int MAINENTRY( void ) {

        char *fname = "OUTSDIF.d"; /* CUTEr data file */
        integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
        integer iout = 6;          /* FORTRAN unit number for error output */
        integer ierr;              /* Exit flag from OPEN and CLOSE */

        VarTypes vtypes;

        integer    *indvar = NULL, *indfun = NULL, ncon_dummy, nnzj;
        doublereal *x, *bl, *bu, *dummy1, *dummy2;
        doublereal *v = NULL, *cl = NULL, *cu = NULL, *c = NULL, *cjac = NULL;
        logical    *equatn = NULL, *linear = NULL;
        char       *pname, *vnames, *gnames;
        logical     efirst = FALSE_, lfirst = FALSE_, nvfrst = FALSE_, grad;
        logical     constrained = FALSE_;

        real        calls[7], cpu[2];
        integer     nlin = 0, nbnds = 0, neq = 0;
        doublereal  dummy;
        integer     ExitCode;
        int         i;

        /* Open problem description file OUTSDIF.d */
        ierr = 0;
        FORTRAN_OPEN( &funit, fname, &ierr );
        if( ierr ) {
            printf("Error opening file OUTSDIF.d.\nAborting.\n");
            exit(1);
        }

        /* Determine problem size */
        CDIMEN( &funit, &CUTEr_nvar, &CUTEr_ncon );

        /* Determine whether to call constrained or unconstrained tools */
        if( CUTEr_ncon ) constrained = TRUE_;

        /* Seems to be needed for some Solaris C compilers */
        ncon_dummy = CUTEr_ncon + 1;

        /* Reserve memory for variables, bounds, and multipliers */
        /* and call appropriate initialization routine for CUTEr */
        MALLOC( x,      CUTEr_nvar, doublereal );
        MALLOC( bl,     CUTEr_nvar, doublereal );
        MALLOC( bu,     CUTEr_nvar, doublereal );
        if( constrained ) {
            MALLOC( equatn, CUTEr_ncon+1, logical    );
            MALLOC( linear, CUTEr_ncon+1, logical    );
            MALLOC( v,      CUTEr_ncon+1, doublereal );
            MALLOC( cl,     CUTEr_ncon+1, doublereal );
            MALLOC( cu,     CUTEr_ncon+1, doublereal );
            CSETUP( &funit, &iout, &CUTEr_nvar, &CUTEr_ncon, x, bl, bu,
                    &CUTEr_nvar, equatn, linear, v, cl, cu, &ncon_dummy,
                    &efirst, &lfirst, &nvfrst );
        } else {
            MALLOC( equatn, 1, logical    );
            MALLOC( linear, 1, logical    );
            MALLOC( cl, 1, doublereal );
            MALLOC( cu, 1, doublereal );
            USETUP( &funit, &iout, &CUTEr_nvar, x, bl, bu, &CUTEr_nvar );
        }

        /* Free unneeded arrays */
        FREE( equatn );
        FREE( linear );

        /* Get problem name */
        MALLOC( pname, FSTRING_LEN+1, char );
        MALLOC( vnames, CUTEr_nvar*FSTRING_LEN, char );
        if( constrained ) {
            MALLOC( gnames, CUTEr_ncon*FSTRING_LEN, char );
            CNAMES( &CUTEr_nvar, &CUTEr_ncon, pname, vnames, gnames );
            FREE( gnames );
        } else {
            UNAMES( &CUTEr_nvar, pname, vnames );
        }
        FREE( vnames );

        /* Make sure to null-terminate problem name */
        pname[FSTRING_LEN] = '\0';
        i = FSTRING_LEN - 1;
        while( i-- > 0 && pname[i] == ' ') {
            pname[i] = '\0';
        }

        /* Obtain basic info on problem */
        if( constrained ) {
            
            GETINFO( CUTEr_nvar, CUTEr_ncon, bl, bu, cl, cu, equatn,
                     linear, &vtypes );
            MALLOC( c, CUTEr_ncon+1, doublereal );
            
            CDIMSJ( &nnzj );
            MALLOC( cjac, nnzj, doublereal );
            MALLOC( indvar, nnzj, integer );
            MALLOC( indfun, nnzj, integer );
            grad = TRUE_;

            CCFSG( &CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon, c, &nnzj,
                   &nnzj, cjac, indvar, indfun, &grad );
            
            printf( "%-10s\t%-10s\t%-10s\n", "cl[i]", "c[i]", "cu[i]" );
            for( i = 0; i < CUTEr_ncon; i++ )
                printf( "%-10g\t%-10g\t%-10g\n", cl[i], c[i], cu[i] );

            printf( "\n" );
            printf( "%-10s\t%-10s\t%-10s\n", "fun", "var", "grad cj" );
            for( i = 0; i < nnzj; i++ )
                printf( "%-10i\t%-10i\t%-10g\n",
                        indfun[i], indvar[i], cjac[i] );
            
        } else {
            
            equatn[0] = FALSE_;
            linear[0] = FALSE_;
            GETINFO( CUTEr_nvar, 1, bl, bu, cl, cu, equatn,
                     linear, &vtypes );
            
        }

        /* Call the optimizer */
        dummy = GENC( ONE );
        ExitCode = 0;

        /* Get CUTEr statistics */
        CREPRT( calls, cpu );

        printf("\n\n ************************ CUTEr statistics ************************\n\n");
        printf(" Code used               : GENC\n");
        printf(" Problem                 : %-s\n", pname);
        printf(" # variables             = %-10d\n", CUTEr_nvar);
        printf(" # constraints           = %-10d\n", CUTEr_ncon);
        printf(" # linear constraints    = %-10d\n", vtypes.nlin);
        printf(" # equality constraints  = %-10d\n", vtypes.neq);
        printf(" # inequality constraints= %-10d\n", vtypes.nineq);
        printf(" # bound constraints     = %-10d\n", vtypes.nbnds);
        printf(" # objective functions   = %-15.7g\n", calls[0]);
        printf(" # objective gradients   = %-15.7g\n", calls[1]);
        printf(" # objective Hessians    = %-15.7g\n", calls[2]);
        printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]);
        printf(" # constraints functions = %-15.7g\n", calls[4]);
        printf(" # constraints gradients = %-15.7g\n", calls[5]);
        printf(" # constraints Hessians  = %-15.7g\n", calls[6]);
        printf(" Exit code               = %-10d\n", ExitCode);
        printf(" Final f                 = %-15.7g\n",dummy);
        printf(" Set up time             = %-10.2f seconds\n", cpu[0]);
        printf(" Solve time              = %-10.2f seconds\n", cpu[1]);
        printf(" ******************************************************************\n\n");

        ierr = 0;
        FORTRAN_CLOSE( &funit, &ierr );
        if( ierr ) {
            printf( "Error closing %s on unit %d.\n", fname, funit );
            printf( "Trying not to abort.\n" );
        }

        /* Free workspace */
        FREE( pname );
        FREE( x ); FREE( bl ); FREE( bu );
        FREE( v ); FREE( cl ); FREE( cu );

        return 0;
    }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
