
/* Generic C solver, to be used with
 * generic C driver, gencma.c
 */

#include "stdio.h"

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cuter.h"

    doublereal genc( doublereal dummy ) {

	printf( "\n\tThis is the generic C solver" );
	printf( "\n\thooked to CUTEr." );
	printf( "\n\tThe magic number is 41.9999995555555\n" );
	return 41.9999995555555;

    }

    void genspc( integer funit, char *fname ) {

	integer ierr;

	/* This is a dummy routine to read a spec file.
	   Possibly, this routine contains precision-dependent directives */

	/* Open relevant file */
	FORTRAN_OPEN( &funit, fname, &ierr );
	if( ierr ) {
	    printf( "Error opening spec file %s.\nAborting.\n", fname );
	    exit(1);
	}

	/* ... Do something ... */

	FORTRAN_CLOSE( &funit, &ierr );
	return;

    }

    void getinfo( integer n, integer m, doublereal *bl, doublereal *bu,
				  doublereal *cl, doublereal *cu, logical *equatn, logical *linear,
				  VarTypes *vartypes ) {

	int i;

	vartypes->nlin = 0; vartypes->neq = 0; vartypes->nbnds = 0;
	vartypes->nrange = 0;
	vartypes->nlower = 0; vartypes->nupper = 0; vartypes->nineq = 0;
	vartypes->nineq_lin = 0; vartypes->nineq_nlin = 0;
	vartypes->neq_lin = 0; vartypes->neq_nlin = 0;

	for( i = 0; i < n; i++ )
	    if( bl[i] > -CUTE_INF || bu[i] < CUTE_INF ) vartypes->nbnds++;
	for( i = 0; i < m; i++ ) {
		if( linear[i] ) vartypes->nlin++;
		if( equatn[i] ) {
			vartypes->neq++;
			if( linear[i] )
				vartypes->neq_lin++;
			else
				vartypes->neq_nlin++;
		} else {
			if( cl[i] > -CUTE_INF ) {
				if( cu[i] < CUTE_INF )
					vartypes->nrange++;
				else {
					vartypes->nlower++; vartypes->nineq++;
                }
			} else {
				if( cu[i] < CUTE_INF ) {
					vartypes->nupper++; vartypes->nineq++;
                }
			}
			if( !equatn[i] && linear[i] ) {
				vartypes->nineq_lin++;
			} else {
				vartypes->nineq_nlin++;
			}
		}
	}
	return;
    }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
