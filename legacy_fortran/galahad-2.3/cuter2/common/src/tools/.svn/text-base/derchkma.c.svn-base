
/* ====================================================
 * CUTEr interface to derivative checker   May 25, 2011
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

#define DERCHKMA

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cuter.h"

#define DERCHK    derchk
#define DERCHKSPC derchkspc
#define GETINFO   getinfo

    doublereal DERCHK(doublereal);
    void DERCHKSPC(integer, char*);
    void GETINFO( integer, integer, doublereal*, doublereal*,
                  doublereal*, doublereal*, logical*, logical*,
                  VarTypes* );

    integer CUTEr_nvar;        /* number of variables */
    integer CUTEr_ncon;        /* number of constraints */
    integer CUTEr_nnzj;        /* number of nonzeros in Jacobian */
    integer CUTEr_nnzh;        /* number of nonzeros in upper triangular
                                  part of the Hessian of the Lagrangian */

    int MAINENTRY(void) {

        char *fname = "OUTSDIF.d"; /* CUTEr data file */
        integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
        integer iout = 6;          /* FORTRAN unit number for error output */
        integer ierr;              /* Exit flag from OPEN and CLOSE */

        VarTypes vtypes;

        integer ncon_dummy;
        doublereal *x, *bl, *bu, *dummy1, *dummy2;
        doublereal *v = NULL, *cl = NULL, *cu = NULL;
        logical *equatn = NULL, *linear = NULL;
        char *pname, *vnames, *gnames, *cptr;
        char **Vnames, **Gnames; /* vnames and gnames as arrays of strings */
        logical efirst = FALSE_, lfirst = FALSE_, nvfrst = FALSE_, grad;
        logical constrained = FALSE_;

        real calls[7], cpu[2];
        integer nlin = 0, nbnds = 0, neq = 0;
        doublereal dummy;
        integer ExitCode;
        int i, j;

        doublereal h, fxp, fxm, approx, derr, xi;
        doublereal *cxp, *cxm, *g;
        int nerr = 0;

        /* Open problem description file OUTSDIF.d */
        ierr = 0;
        FORTRAN_OPEN(&funit, fname, &ierr);
        if (ierr) {
            printf("Error opening file OUTSDIF.d.\nAborting.\n");
            exit(1);
        }

        /* Determine problem size */
        CDIMEN(&funit, &CUTEr_nvar, &CUTEr_ncon);

        /* Determine whether to call constrained or unconstrained tools */
        if (CUTEr_ncon) constrained = TRUE_;

        /* Seems to be needed for some Solaris C compilers */
        ncon_dummy = CUTEr_ncon + 1;

        /* Reserve memory for variables, bounds, and multipliers */
        /* and call appropriate initialization routine for CUTEr */
        MALLOC(x,      CUTEr_nvar, doublereal);
        MALLOC(bl,     CUTEr_nvar, doublereal);
        MALLOC(bu,     CUTEr_nvar, doublereal);
        if (constrained) {
            MALLOC(equatn, CUTEr_ncon+1, logical);
            MALLOC(linear, CUTEr_ncon+1, logical);
            MALLOC(v,      CUTEr_ncon+1, doublereal);
            MALLOC(cl,     CUTEr_ncon+1, doublereal);
            MALLOC(cu,     CUTEr_ncon+1, doublereal);
            CSETUP(&funit, &iout, &CUTEr_nvar, &CUTEr_ncon, x, bl, bu,
                &CUTEr_nvar, equatn, linear, v, cl, cu, &ncon_dummy,
                &efirst, &lfirst, &nvfrst);
        } else {
            MALLOC(equatn, 1, logical);
            MALLOC(linear, 1, logical);
            MALLOC(cl, 1, doublereal);
            MALLOC(cu, 1, doublereal);
            USETUP(&funit, &iout, &CUTEr_nvar, x, bl, bu, &CUTEr_nvar);
        }

        /* Get problem, variables and constraints names */
        MALLOC(pname, FSTRING_LEN+1, char);
            MALLOC(vnames, CUTEr_nvar * FSTRING_LEN, char);  /* For Fortran */
            MALLOC(Vnames, CUTEr_nvar, char*);          /* Array of strings */
            for (i = 0; i < CUTEr_nvar; i++)
            MALLOC(Vnames[i], FSTRING_LEN+1, char);

        if (constrained) {
            MALLOC(gnames, CUTEr_ncon * FSTRING_LEN, char);   /* For Fortran */
            MALLOC(Gnames, CUTEr_ncon, char*);        /* Array of strings */
            for (i = 0; i < CUTEr_ncon; i++)
                MALLOC(Gnames[i], FSTRING_LEN+1, char);
            CNAMES(&CUTEr_nvar, &CUTEr_ncon, pname, vnames, gnames);
        } else {
            UNAMES(&CUTEr_nvar, pname, vnames);
        }

        /* Make sure to null-terminate problem name */
        pname[FSTRING_LEN] = '\0';

        /* Transfer variables and constraint names into arrays of
         * null-terminated strings.
         * If you know of a simpler way to do this portably, let me know!
         */
        for (i = 0; i < CUTEr_nvar; i++) {
          cptr = vnames + i * FSTRING_LEN;
          for (j = 0; j < FSTRING_LEN; j++) {
            Vnames[i][j] = *cptr;
            cptr++;
          }
          Vnames[i][FSTRING_LEN] = '\0';
        }

        for (i = 0; i < CUTEr_ncon; i++) {
          cptr = vnames + i * FSTRING_LEN;
          for (j = 0; j < FSTRING_LEN; j++) {
            Gnames[i][j] = *cptr;
            cptr++;
          }
          Gnames[i][FSTRING_LEN] = '\0';
        }

        /* Fortran strings no longer needed */
        FREE(vnames);
        if (constrained) FREE(gnames);

    /* Obtain basic info on problem */
    if (constrained) {
        GETINFO( CUTEr_nvar, CUTEr_ncon, bl, bu, cl, cu, equatn,
             linear, &vtypes );
    } else {
        equatn[0] = FALSE_;
        linear[0] = FALSE_;
        GETINFO( CUTEr_nvar, 1, bl, bu, cl, cu, equatn,
             linear, &vtypes );
    }

    /* Scramble initial point. */
    for (i = 0; i < CUTEr_nvar; i++)
        x[i] = pow(-1,i) * (2*i+1);

    /* Evaluate gradient at initial point. */
    MALLOC(g, CUTEr_nvar, doublereal);
    UGR(&CUTEr_nvar, x, g);

    /* Check first derivatives of objective and constraints. */
    h = 1.0e-8;
    if (constrained) {
        MALLOC(cxp, CUTEr_ncon+1, doublereal);
        MALLOC(cxm, CUTEr_ncon+1, doublereal);
    }
    fprintf(stderr, "%10s  %22s  %22s  %7s\n",
            "Variable", "AD", "Numerical", "Error");
    for (i = 0; i < CUTEr_nvar; i++) {
        xi = x[i];
        x[i] = xi + h;
        if (constrained)
            CFN(&CUTEr_nvar, &CUTEr_ncon, x, &fxp, &CUTEr_ncon, cxp);
        else
            UFN(&CUTEr_nvar, x, &fxp);
        x[i] = xi - h;
        if (constrained)
            CFN(&CUTEr_nvar, &CUTEr_ncon, x, &fxm, &CUTEr_ncon, cxm);
        else
            UFN(&CUTEr_nvar, x, &fxm);
        x[i] = xi;      /* Restore x */

        /* Check i-th derivative of objective */
        approx = (fxp-fxm)/(2*h);
        derr = fabs(g[i] - approx)/fmax(1,fabs(g[i]));
        if (derr > 1.0e-4) {
            nerr++;
            fprintf(stderr, "%10s  %22.15e  %22.15e  %7.1e\n",
                    Vnames[i], g[i], approx, derr);
        }
    }
    FREE(g);
    if (constrained) {
        FREE(cxp);
        FREE(cxm);
    }

    /* Get CUTEr statistics */
    CREPRT(calls, cpu);

    printf("\n\n ************************ CUTEr statistics");
    printf(" ************************\n\n");
    printf(" Code used               : GENC\n");
    printf(" Problem                 : %-s\n", pname);
    printf(" # variables             = %-10d\n", (int)CUTEr_nvar);
    printf(" # constraints           = %-10d\n", (int)CUTEr_ncon);
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
    printf(" Exit code               = %-10d\n", nerr);
    printf(" Final f                 = %-15.7g\n", fxp);
    printf(" Set up time             = %-10.2f seconds\n", cpu[0]);
    printf(" Solve time              = %-10.2f seconds\n", cpu[1]);
    printf(" *********************************************");
    printf("*********************\n\n");

    ierr = 0;
    FORTRAN_CLOSE(&funit, &ierr);
    if (ierr)
        printf("Error closing %s on unit %d.\n", fname, (int)funit);

    /* Free workspace */
    FREE(pname);
    FREE(x); FREE(bl); FREE(bu);
    FREE(v); FREE(cl); FREE(cu);
    FREE(equatn);
    FREE(linear);

    /* Free memory for variable names */
    for (i = 0; i < CUTEr_nvar; i++) FREE(Vnames[i]);
    FREE(Vnames);

    /* Free memory for constraint names */
    for (i = 0; i < CUTEr_ncon; i++) FREE(Gnames[i]);
    if (constrained) FREE(Gnames);

    return 0;

    }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
