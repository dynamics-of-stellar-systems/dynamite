
/* ====================================================
 * CUTEr interface to KNITRO 7           April 14, 2011
 *
 * D. Orban
 *
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define KNITROMA

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cuter.h"
#include "knitro.h"

  logical somethingTrue = TRUE_;
  logical somethingFalse = FALSE_;
  integer CUTEr_nvar;        /* number of variables */
  integer CUTEr_ncon;        /* number of constraints */
  integer CUTEr_lcjac;       /* length of Jacobian arrays */
  integer CUTEr_nnzj;        /* number of nonzeros in Jacobian */
  integer CUTEr_nnzh;        /* number of nonzeros in upper triangular
                                part of the Hessian of the Lagrangian */
  integer *jacIndexVars, *jacIndexCons, *hessIndexRows, *hessIndexCols;
  doublereal *CUTEr_Jac, *CUTEr_Hess, *Hv, f;

  /* ======================================================================== */
  /* Callback function to evaluate the objective function and constraints */
  /* ======================================================================== */

  int  callbackEvalFC (const int             evalRequestCode,
                       const int             n,
                       const int             m,
                       const int             nnzJ,
                       const int             nnzH,
                       const double * const  x,
                       const double * const  lambda,
                       double * const  obj,
                       double * const  c,
                       double * const  objGrad,
                       double * const  jac,
                       double * const  hessian,
                       double * const  hessVector,
                       void   *        userParams) {

    int i;

    if (evalRequestCode != KTR_RC_EVALFC) {
      fprintf (stderr, "*** callbackEvalFC incorrectly called with code %d\n",
               evalRequestCode);
      return -1;
    }

    if (CUTEr_ncon > 0) {
      CFN(&CUTEr_nvar, &CUTEr_ncon, x, obj, &CUTEr_ncon, c);
      //for (i = 0; i < CUTEr_nvar; i++)
      //    printf(" x = %5d  %22.15e\n", i, x[i]);
      //printf(" f = %22.15e\n", obj[0]);
      //for (i = 0; i < CUTEr_ncon; i++)
      //    printf(" c = %5d  %22.15e\n", i, c[i]);
    } else {
      UFN(&CUTEr_nvar, x, obj);
    }

    return 0;
  }

  /* ======================================================================== */
  /* Callback function to evaluate the objective gradient and Jacobian */
  /* ======================================================================== */

  int  callbackEvalGA (const int             evalRequestCode,
                       const int             n,
                       const int             m,
                       const int             nnzJ,
                       const int             nnzH,
                       const double * const  x,
                       const double * const  lambda,
                       double * const  obj,
                       double * const  c,
                       double * const  objGrad,
                       double * const  jac,
                       double * const  hessian,
                       double * const  hessVector,
                       void   *        userParams) {

    int i;

    if (evalRequestCode != KTR_RC_EVALGA) {
      fprintf (stderr, "*** callbackEvalGA incorrectly called with code %d\n",
               evalRequestCode);
      return -1;
    }

    if (CUTEr_ncon > 0) {

      CSGR(&CUTEr_nvar, &CUTEr_ncon, &somethingFalse, &CUTEr_ncon, lambda,
           x, &CUTEr_nnzj, &CUTEr_lcjac, CUTEr_Jac, jacIndexVars, jacIndexCons);
      for (i = 0; i < CUTEr_nvar; i++) objGrad[i] = 0.0;
      for (i = 0; i < CUTEr_nnzj; i++)
        if (jacIndexCons[i] == 0)
          objGrad[ (int)(jacIndexVars[i])-1] = CUTEr_Jac[i];
        else
          jac[i] = CUTEr_Jac[i];
      //CCFSG(&CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon, c, &CUTEr_nnzj,
      //      &CUTEr_lcjac, jac, jacIndexVars, jacIndexCons, &somethingTrue);
      //COFG(&CUTEr_nvar, x, &f, objGrad, &somethingTrue);

    } else
      UGR(&CUTEr_nvar, x, objGrad);

    // for (i = 0; i < CUTEr_nvar; i++)
    //  printf(" g = %22.15e\n", objGrad[i]);

    return 0;
  }

  /* ======================================================================== */
  /* Callback function to evaluate the Hessian or Hessian-vector product */
  /* ======================================================================== */

  int  callbackEvalHess (const int             evalRequestCode,
                         const int             n,
                         const int             m,
                         const int             nnzJ,
                         const int             nnzH,
                         const double * const  x,
                         const double * const  lambda,
                         double * const  obj,
                         double * const  c,
                         double * const  objGrad,
                         double * const  jac,
                         double * const  hessian,
                         double * const  hessVector,
                         void   *        userParams) {

    int i;

    if (evalRequestCode == KTR_RC_EVALH) {

      if (CUTEr_ncon > 0) {
        CSH(&CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon, lambda, &CUTEr_nnzh,
             &CUTEr_nnzh, hessian, hessIndexRows, hessIndexCols);
      } else {
        USH(&CUTEr_nvar, x, &CUTEr_nnzh,
             &CUTEr_nnzh, hessian, hessIndexRows, hessIndexCols);
      }

      /* Adjust indices to be zero-based. */
      for (i = 0; i < CUTEr_nnzh; i++) {
          hessIndexRows[i] = (int)(hessIndexRows[i] - 1);
          hessIndexCols[i] = (int)(hessIndexCols[i] - 1);
          //printf(" H = %5d  %5d  %22.15e\n", hessIndexRows[i],
          //        hessIndexCols[i], hessian[i]);

      }

    } else if (evalRequestCode == KTR_RC_EVALHV) {

      if (! Hv) MALLOC(Hv, CUTEr_nvar, doublereal);

      if (CUTEr_ncon > 0)
        CPROD(&CUTEr_nvar, &CUTEr_ncon, &somethingTrue, x, &CUTEr_ncon, \
              lambda, hessVector, Hv);
      else
        UPROD(&CUTEr_nvar, &somethingTrue, x, hessVector, Hv);

      for (i = 0; i < CUTEr_nvar; i++) hessVector[i] = Hv[i];

    } else {
      fprintf(stderr, "*** callbackEvalHess incorrectly called with code %d\n",
              evalRequestCode);
      return -1;
    }

    return 0;
  }

  /* ======================================================================== */
  /* Main function */
  /* ======================================================================== */

  int MAINENTRY(void) {

    char *fname = "OUTSDIF.d"; /* CUTEr data file */
    integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
    integer iout = 6;          /* FORTRAN unit number for error output */
    integer ierr;              /* Exit flag from OPEN and CLOSE */

    KTR_context *KnitroData;
    char      szVersion[15 + 1];

    VarTypes vtypes;

    integer ncon_dummy;
    doublereal *x, *bl, *bu, *dummy1, *dummy2;
    doublereal *v = NULL, *cl = NULL, *cu = NULL;
    logical *equatn = NULL, *linear = NULL;
    char *pname, *vnames, *gnames;
    logical efirst = TRUE_, lfirst = TRUE_, nvfrst = FALSE_, grad;
    logical constrained = FALSE_;

    doublereal *c, f;

    real calls[7], cpu[2];
    integer nlin = 0, nbnds = 0, neq = 0;
    doublereal dummy;
    integer ExitCode;
    int nHessOpt, i;

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
      MALLOC(equatn, CUTEr_ncon+1,            logical   );
      MALLOC(linear, CUTEr_ncon+1,            logical   );
      MALLOC(v,      CUTEr_ncon+CUTEr_nvar+1, doublereal);
      MALLOC(cl,     CUTEr_ncon+1,            doublereal);
      MALLOC(cu,     CUTEr_ncon+1,            doublereal);
      CSETUP(&funit, &iout, &CUTEr_nvar, &CUTEr_ncon, x, bl, bu,
              &CUTEr_nvar, equatn, linear, v, cl, cu, &ncon_dummy,
              &efirst, &lfirst, &nvfrst);
    } else {
      MALLOC(equatn, 1,            logical   );
      MALLOC(linear, 1,            logical   );
      MALLOC(cl,     1,            doublereal);
      MALLOC(cu,     1,            doublereal);
      MALLOC(v,      CUTEr_nvar+1, doublereal);
      USETUP(&funit, &iout, &CUTEr_nvar, x, bl, bu, &CUTEr_nvar);
    }

    /* Get problem name */
    MALLOC(pname, FSTRING_LEN+1, char);
    MALLOC(vnames, CUTEr_nvar*FSTRING_LEN, char);
    if (constrained) {
      MALLOC(gnames, CUTEr_ncon*FSTRING_LEN, char);
      CNAMES(&CUTEr_nvar, &CUTEr_ncon, pname, vnames, gnames);
      FREE(gnames);
    } else
      UNAMES(&CUTEr_nvar, pname, vnames);
    FREE(vnames);

    /* Make sure to null-terminate problem name */
    pname[FSTRING_LEN] = '\0';
    i = FSTRING_LEN - 1;
    while(i-- > 0 && pname[i] == ' ') {
      pname[i] = '\0';
    }

    /* Obtain basic info on problem */
    /*
    if (constrained)
      GETINFO(CUTEr_nvar, CUTEr_ncon, bl, bu, cl, cu, equatn, linear, &vtypes);
    else {
      equatn[0] = FALSE_;
      linear[0] = FALSE_;
      GETINFO(CUTEr_nvar, 1, bl, bu, cl, cu, equatn, linear, &vtypes);
    }
    */

    /* Initialize KNITRO context */
    KTR_get_release(15, szVersion);
    KnitroData = KTR_new();
    if (KnitroData == NULL) {
      fprintf(stderr, "Failed to find a valid KNITRO license.\n");
      fprintf(stderr, "%s\n", szVersion);
      exit(1);
    }

    /* Read spec file */
    if (KTR_load_param_file(KnitroData, "knitro.opt")) {
      fprintf(stderr, "Cannot open spec file knitro.opt...");
      fprintf(stderr, " using all default options\n");
    }
    if (KTR_get_int_param_by_name(KnitroData, "hessopt", &nHessOpt))
      exit(-1);

    /* Obtain Jacobian sparsity pattern and initial objective value */
#ifdef KNIT_DEBUG
    fprintf(stderr, "Obtaining Jacobian sparsity pattern...\n");
#endif

    if (constrained) {

      // Constrained problem.

      CDIMSJ(&CUTEr_lcjac);

      MALLOC(jacIndexVars, CUTEr_lcjac + 1, integer   );
      MALLOC(jacIndexCons, CUTEr_lcjac + 1, integer   );
      MALLOC(CUTEr_Jac,    CUTEr_lcjac + 1, doublereal);

      MALLOC(c, CUTEr_ncon, doublereal);

      CCFSG(&CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon, c, &CUTEr_nnzj,
            &CUTEr_lcjac, CUTEr_Jac, jacIndexVars, jacIndexCons, &somethingTrue);
      CFN(&CUTEr_nvar, &CUTEr_ncon, x, &f, &CUTEr_ncon, c);

      FREE(c);

      /* Convert to 0-based indexing */
      for (i = 0; i < CUTEr_nnzj; i++) {
        jacIndexVars[i] = (int)(jacIndexVars[i] - 1);
        jacIndexCons[i] = (int)(jacIndexCons[i] - 1);
        //   printf("jacIndexVars[%d]=%d\n", i, jacIndexVars[i]);
        // printf("jacIndexCons[%d]=%d\n", i, jacIndexCons[i]);
      }

    } else {

      // Unconstrained problem.

      jacIndexVars = NULL;
      jacIndexCons = NULL;
      CUTEr_Jac = NULL;
      UFN(&CUTEr_nvar, x, &f);
    }

    /* Obtain Hessian sparsity pattern */
#ifdef KNIT_DEBUG
    fprintf(stderr, "Obtaining Hessian sparsity pattern...\n");
#endif
    CDIMSH(&CUTEr_nnzh);

    MALLOC(hessIndexRows, CUTEr_nnzh, integer   );
    MALLOC(hessIndexCols, CUTEr_nnzh, integer   );
    MALLOC(CUTEr_Hess,    CUTEr_nnzh, doublereal);

    CSH(&CUTEr_nvar, &CUTEr_ncon, x, &CUTEr_ncon, v, &CUTEr_nnzh, &CUTEr_nnzh,
         CUTEr_Hess, hessIndexRows, hessIndexCols);

    /* Convert to 0-based indexing */
    for (i = 0; i < CUTEr_nnzh; i++) {
      hessIndexRows[i] = (int)(hessIndexRows[i] - 1);
      hessIndexCols[i] = (int)(hessIndexCols[i] - 1);
      // printf("hessIndexRows[%d]=%d\n", i, hessIndexRows[i]);
      // printf("hessIndexCols[%d]=%d\n", i, hessIndexCols[i]);
    }

    /* Register callback functions */
#ifdef KNIT_DEBUG
    fprintf(stderr, "Registering callback functions...\n");
#endif
    if (KTR_set_func_callback(KnitroData, &callbackEvalFC)) {
      fprintf(stderr, "Could not register FC callback function\n");
      goto terminate;
    }
    if (KTR_set_grad_callback(KnitroData, &callbackEvalGA)) {
      fprintf(stderr, "Could not register FC callback function\n");
      goto terminate;
    }
    if ((nHessOpt == KTR_HESSOPT_EXACT) || (nHessOpt == KTR_HESSOPT_PRODUCT))
      if (KTR_set_hess_callback(KnitroData, &callbackEvalHess)) {
        fprintf(stderr, "Could not register Hess callback function\n");
        goto terminate;
      }

    // Convert infinite bounds.
    for (i = 0; i < CUTEr_nvar; i++) {
        if (bl[i] == -CUTE_INF) bl[i] = -KTR_INFBOUND;
        if (bu[i] ==  CUTE_INF) bu[i] =  KTR_INFBOUND;
    }
    for (i = 0; i < CUTEr_ncon; i++) {
        if (cl[i] == -CUTE_INF) cl[i] = -KTR_INFBOUND;
        if (cl[i] ==  CUTE_INF) cu[i] =  KTR_INFBOUND;
    }

#ifdef KNIT_DEBUG
    fprintf(stderr, "Initializing KNITRO data structure...\n");
#endif
    ExitCode = KTR_init_problem(KnitroData,
                                 CUTEr_nvar,
                                 KTR_OBJGOAL_MINIMIZE,
                                 KTR_OBJTYPE_GENERAL,
                                 bl, bu,
                                 (int)CUTEr_ncon,
                                 (int *)linear,
                                 cl, cu,
                                 (int)CUTEr_nnzj,
                                 (int *)jacIndexVars,
                                 (int *)jacIndexCons,
                                 (int)CUTEr_nnzh,
                                 (int *)hessIndexRows,
                                 (int *)hessIndexCols,
                                 x,
                                 NULL);

#ifdef KNIT_DEBUG
    fprintf(stderr, "Exit code from initialization: %d\n", ExitCode);
#endif

    /* Call the optimizer */
    ExitCode = KTR_solve(KnitroData, x, v, 0, &f,
                          NULL, NULL, NULL, NULL, NULL, NULL);

    /* Release KNITRO context */
#ifdef KNIT_DEBUG
    fprintf(stderr, "Terminating...\n");
#endif
    if (KTR_free(&KnitroData))
      fprintf(stderr, "Knitro-CUTEr:: Error freeing Knitro data structure\n");

    /* Get CUTEr statistics */
    CREPRT(calls, cpu);

    printf("\n\n ************************ CUTEr statistics ************************\n\n");
    printf(" Code used               : %-15s\n", szVersion);
    printf(" Problem                 : %-s\n", pname);
    printf(" # variables             = %-10d\n", CUTEr_nvar);
    printf(" # constraints           = %-10d\n", CUTEr_ncon);
    /*
    printf(" # linear constraints    = %-10d\n", vtypes.nlin);
    printf(" # equality constraints  = %-10d\n", vtypes.neq);
    printf(" # inequality constraints= %-10d\n", vtypes.nineq);
    printf(" # bound constraints     = %-10d\n", vtypes.nbnds);
    */
    printf(" # objective functions   = %-15.7g\n", calls[0]);
    printf(" # objective gradients   = %-15.7g\n", calls[1]);
    printf(" # objective Hessians    = %-15.7g\n", calls[2]);
    printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]);
    printf(" # constraints functions = %-15.7g\n", calls[4]);
    printf(" # constraints gradients = %-15.7g\n", calls[5]);
    printf(" # constraints Hessians  = %-15.7g\n", calls[6]);
    printf(" Exit code               = %-10d\n", ExitCode);
    printf(" Final f                 = %-15.7g\n", f);
    printf(" Set up time             = %-10.2f seconds\n", cpu[0]);
    printf(" Solve time              = %-10.2f seconds\n", cpu[1]);
    printf(" ******************************************************************\n\n");

  terminate:

    ierr = 0;
    FORTRAN_CLOSE(&funit, &ierr);
    if (ierr) {
      fprintf(stderr, "Error closing %s on unit %d.\n", fname, funit);
      fprintf(stderr, "Trying not to abort.\n");
    }

    /* Free workspace */
    FREE(pname);
    FREE(x); FREE(bl); FREE(bu);
    FREE(v); FREE(cl); FREE(cu);
    FREE(equatn);
    FREE(linear);

    if (constrained) {
      FREE(CUTEr_Jac);
      FREE(jacIndexVars);
      FREE(jacIndexCons);
    }
    FREE(hessIndexRows);
    FREE(hessIndexCols);
    FREE(CUTEr_Hess);
    if (Hv) FREE(Hv);
    return 0;

  }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif

