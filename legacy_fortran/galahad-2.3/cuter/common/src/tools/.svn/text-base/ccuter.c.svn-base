
/* ================================================== */
/* Auxilliary functions                               */
/* ================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cuter.h"

/* Track origin of error messages */

#define CUTEr_ERRQ(errcode,msg) {                                    \
  printf( "CUTEr C error:: Code = %d, Msg :: %s\n", errcode, msg );  \
  printf( "Error occured in function %s, file %s at line %d\n",      \
          __FUNCT__, __FILE__, __LINE__ );                           \
  exit( errcode );                                                   \
}

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "CUTEr_malloc"
    void *CUTEr_malloc( void *object, int length, size_t s ) {
	object = malloc( length * s );
	if( !object ) {
	    if( !length ) object = malloc( sizeof(double) );
	    if( !object ) CUTEr_ERRQ(-1,"Unable to allocate memory");
	}
	return object;
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "CUTEr_calloc"
    void *CUTEr_calloc( void *object, int length, size_t s ) {
	object = calloc( (size_t)length, s );
	if( !object ) {
	    if( !length ) object = calloc( (size_t)1, s );
	    if( !object ) CUTEr_ERRQ(-2,"Unable to allocate pointer");
	}
	return object;
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "CUTEr_realloc"
    void *CUTEr_realloc( void *object, int length, size_t s ) {
	object = realloc( object, length * s );
	if( !object ) {
	    if( !length ) object = realloc( object, sizeof( double ) );
	    if( !object ) CUTEr_ERRQ(-3,"Unable to reallocate");
	}
	return object;
    }

/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "CUTEr_free"
    void CUTEr_free( void *object ) {
	if( object ) free( object );
    }


/* -------------------------------------------------- */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif

#ifdef __cplusplus
}              /* Closing brace for  extern "C"  block */
#endif
