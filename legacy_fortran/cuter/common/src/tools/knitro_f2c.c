/*******************************************************/
/* Copyright (c) 2001-2004 by Ziena Optimization, Inc. */
/*               2001-2004 by Northwestern University  */
/* All Rights Reserved                                 */
/*******************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "knitro.h"

/* KNITRO gateway routine for Fortran interface */

 /* This is some glue so we can use the default external symbol
     convention on both Win32 and GNU tools. */
#ifdef _WIN32
# define ktrsolvef_ ktrsolvef
#endif

#ifdef __cplusplus
extern "C" {   /* Prevent C++ compiler from mangling symbols */
#endif

/* This is called from user fortran code. */

void KNITRO_4_API ktrsolvef_(double *f, int *ftype, 
			     int *n, double *x, double *bl, double *bu, double *fgrad, 
			     int *m, double *c, double *cl, double *cu, int *ctype,
			     int *nnzj, double *cjac, int *indvar, int *indfun, 
			     double *lambda, 
			     int *nnzh, double *hess, int *hrow, int *hcol, 
			     double *vector, int *status,
			     int *gradopt, int *hessopt)
{
  static KTR_context *kc;
  int errorflag;

  /* Do some special initialization on first call */
  
  if (*status == KTR_RC_INITIAL) {
    
    /* Create problem instance */
    kc = KTR_new(1);
    
    /* Set user options here */

    /* Uncomment below to read or write options file */
    KTR_load_param_file(kc,"knitro.opt");
    //KTR_save_param_file(kc,"knitro.opt");

    /* Set the options hessopt and gradopt
       specified in the Fortran driver file */

    errorflag = KTR_set_int_param(kc,KTR_PARAM_HESSOPT,*hessopt);
    if (errorflag) {
      fprintf (stderr, "Failure to get param. %d file=%s line=%d\n",
	       errorflag,__FILE__,__LINE__);
      exit(1);
    }
    errorflag = KTR_set_int_param(kc,KTR_PARAM_GRADOPT,*gradopt);
    if (errorflag) {
      fprintf (stderr, "Failure to get param. %d file=%s line=%d\n",
	       errorflag,__FILE__,__LINE__);
      exit(1);
    }
  }

  /* Call KNITRO solver routine */

  *status = KTR_solve(kc, f, *ftype, *n, x, bl, bu, fgrad,
		      *m, c, cl, cu, ctype, 
		      *nnzj, cjac, indvar, indfun, lambda, 
		      *nnzh, hess, hrow, hcol, vector, NULL);

  /* If optimization finished, delete problem instance  */
  
  if (*status <= 0) {
    KTR_free(&kc);
  }

}

#ifdef __cplusplus
    }   /* End of Extern "C" block */
#endif
