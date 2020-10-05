/* 
The calling syntax :

   [L, Ldiag] = icfmex(n, lA, Adiag, p) ;

*/
#if ARCH == aix-4
#define dicfs_ dicfs
#endif

#include <stdio.h>
#include <stdlib.h> 

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int n, p ; 
   double *a, *l;
   double *adiag, *ldiag;
   int *acol_ptr, *arow_ind ;
   int *lcol_ptr, *lrow_ind ;

   /* Working arrays */

   double *w1, *w2 ; 
   int *iw ;

   /* Local variables */

   int i, nnz ;
   double alpha ;

   /* Right hand side */

   n = (int) mxGetScalar(prhs[0]);  
   a =  mxGetPr(prhs[1]);  
   arow_ind =  mxGetIr(prhs[1]);  
   acol_ptr =  mxGetJc(prhs[1]);  
   adiag = mxGetPr(prhs[2]) ;      
   p = (int) mxGetScalar(prhs[3]) ;

   nnz = acol_ptr[n] ;

   /* Left-hand side */

   if (plhs[0] != NULL) mxDestroyArray(plhs[0]) ;
   plhs[0] = mxCreateSparse(n, n, nnz+n*p, mxREAL) ;

   l =  mxGetPr(plhs[0]);  
   lrow_ind =  mxGetIr(plhs[0]);  
   lcol_ptr =  mxGetJc(plhs[0]);  

   if (plhs[1] != NULL) mxDestroyArray(plhs[1]) ;
   plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL) ;
   ldiag = mxGetPr(plhs[1]) ;

   /* Allocate working arrays */

   if (( (w1=(double *)calloc(n, sizeof(double)))==(double *)NULL ) ||
       ( (w2=(double *)calloc(n, sizeof(double)))==(double *)NULL ) ||
       ( (iw=(int *)calloc(3*n, sizeof(int)))==(int *)NULL )) {
     printf("Not enough memory\n") ;
   } 

   /* Change a's i j into fortran index */

   for (i = 0; i < nnz; i++) arow_ind[i] = arow_ind[i] + 1 ;
   for (i = 0; i <= n; i++)  acol_ptr[i] = acol_ptr[i] + 1 ;

   /* Call icf fortran subroutine */

   alpha = 0.0 ;

   dicfs_(&n, &nnz, a, adiag, acol_ptr, arow_ind,
	  l, ldiag, lcol_ptr, lrow_ind, &p, &alpha,
	  iw, w1, w2) ;


   /* Change a's i j into C index */
   for (i = 0; i < nnz; i++) 
      arow_ind[i] = arow_ind[i] - 1 ;
   for (i = 0; i < lcol_ptr[n]; i++) 
      lrow_ind[i] = lrow_ind[i] - 1 ;
   for (i = 0; i <= n; i++) {
      acol_ptr[i] = acol_ptr[i] - 1 ;
      lcol_ptr[i] = lcol_ptr[i] - 1 ;
   }

   /* Free memory */

   free(w1) ;   
   free(w2) ;   
   free(iw) ;   
} /* mexFunction */

