/* 
The calling syntax :

   [x, iters, info] = cgmex(n, lA, Adiag, lL, Ldiag, b, maxiter, rtol) ;

*/
#if ARCH == aix-4
#define dpcg_ dpcg
#endif

#include <stdio.h>
#include <stdlib.h> 

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Right hand side */

   int n ;
   double *a, *l, *adiag, *ldiag ;
   int *acol_ptr, *arow_ind, *lcol_ptr, *lrow_ind ;
   double *b, rtol ; 
   int maxiter ;
   
   /* Left hand side */

   double *x, *iters, *info ;

   /* Working arrays */
   
   double *w ;

   /* Local variables */

   int i, nnza, nnzl, info_dpcg, iters_dpcg ; 

   /* Right hand side */

   n = (int) mxGetScalar(prhs[0]);  
   a =  mxGetPr(prhs[1]);  
   arow_ind =  mxGetIr(prhs[1]);  
   acol_ptr =  mxGetJc(prhs[1]);  
   adiag = mxGetPr(prhs[2]) ;      
   l =  mxGetPr(prhs[3]);  
   lrow_ind =  mxGetIr(prhs[3]);  
   lcol_ptr =  mxGetJc(prhs[3]);  
   ldiag = mxGetPr(prhs[4]) ;      
   b =  mxGetPr(prhs[5]);  
   maxiter = (int) mxGetScalar(prhs[6]);  
   rtol = mxGetScalar(prhs[7]);  

   /* Change a and l's i j into fortran index */

   nnza = acol_ptr[n] ;
   nnzl = lcol_ptr[n] ;

   for (i = 0; i < nnza; i++) arow_ind[i] = arow_ind[i] + 1 ;
   for (i = 0; i < nnzl; i++) lrow_ind[i] = lrow_ind[i] + 1 ;
   for (i = 0; i <= n; i++) {
      acol_ptr[i] = acol_ptr[i] + 1 ; 
      lcol_ptr[i] = lcol_ptr[i] + 1 ; 
   }

   /* Left hand side */

   if (plhs[0] != NULL) mxDestroyArray(plhs[0]) ;
   plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL) ;
   x = mxGetPr(plhs[0]) ;  

   if (plhs[1] != NULL) mxDestroyArray(plhs[1]) ;
   plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL) ;
   iters = mxGetPr(plhs[1]) ;  

   if (plhs[2] != NULL) mxDestroyArray(plhs[2]) ;
   plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL) ;
   info = mxGetPr(plhs[2]) ;  

   /* Allocate working memory */

   if (( (w =(double *)calloc(4*n, sizeof(double)))==(double *)NULL )) {
     printf("Not enough memory\n") ;
   } 

   dpcg_(&n,a,adiag,acol_ptr,arow_ind,l,ldiag,
	 lcol_ptr,lrow_ind,b,&rtol,&maxiter,x,
	 &iters_dpcg,&info_dpcg,w,w+n,w+2*n,w+3*n) ;

   *iters = (double) iters_dpcg ;
   *info = (double) info_dpcg ;

   /* Change a's i j into C index */

   for (i = 0; i < nnza; i++) arow_ind[i] = arow_ind[i] - 1 ;
   for (i = 0; i < nnzl; i++) lrow_ind[i] = lrow_ind[i] - 1 ;
   for (i = 0; i <= n; i++) {
      acol_ptr[i] = acol_ptr[i] - 1 ; lcol_ptr[i] = lcol_ptr[i] - 1 ; 
   }


   /* Free memory */

   free(w) ;   
} /* mexFunction */

