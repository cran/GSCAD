

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <assert.h>

//#ifdef HAVE_MKL   // obsolete
//#include <mkl_cblas.h>
//#else
//#include "cblas.h"
//#endif
//#ifdef USE_BLAS_LIB
//#include "blas.h"
//#else
//#include "cblas.h"  // dependency upon cblas libraries has been removed in a recent version
//#endif

#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef MATLAB_MEX_FILE
typedef int mwSize;
#endif

#ifndef MAX_THREADS
#define MAX_THREADS 64
#endif

// MIN, MAX macros
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SIGN(a) (((a) < 0) ? -1.0 : 1.0)
#define ABS(a) (((a) < 0) ? -(a) : (a))
// DEBUG macros
#define PRINT_I(name) printf(#name " : %d\n",name);
#define PRINT_F(name) printf(#name " : %g\n",name);
#define PRINT_S(name) printf("%s\n",name);
#define FLAG(a) printf("flag : %d \n",a);

// ALGORITHM constants
#define EPSILON 10e-10
#ifndef INFINITY
#define INFINITY 10e20
#endif
#define EPSILON_OMEGA 0.001
#define TOL_CGRAD 10e-6
#define MAX_ITER_CGRAD 40



#endif
