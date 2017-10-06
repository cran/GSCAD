// The code in this file as well as the included header files are based on
// the SPAMS software. Some modification was made. The original code can be
// found at http://spams-devel.gforge.inria.fr/downloads.html.

// #include <Rcpp.h>
#include <RcppArmadillo.h>
// #include "linalg/linalg.h"
#include "decomp.h"
extern "C" {
#if R_VERSION >= R_Version(2,6,0)
#define VMAXTYPE void *
#else
#define VMAXTYPE char *
#endif
}


#define SWIG_OK                    (0)
#define SWIG_ERROR                 (-1)
#define SWIG_IsOK(r)               (r >= 0)
#define SWIG_ArgError(r)           ((r != SWIG_ERROR) ? r : SWIG_TypeError)

/*  Errors in SWIG */
#define  SWIG_UnknownError    	   -1
#define  SWIG_IOError        	   -2
#define  SWIG_RuntimeError   	   -3
#define  SWIG_IndexError     	   -4
#define  SWIG_TypeError      	   -5
#define  SWIG_DivisionByZero 	   -6
#define  SWIG_OverflowError  	   -7
#define  SWIG_SyntaxError    	   -8
#define  SWIG_ValueError     	   -9
#define  SWIG_SystemError    	   -10
#define  SWIG_AttributeError 	   -11
#define  SWIG_MemoryError    	   -12
#define  SWIG_NullReferenceError   -13

static int R_result_pos = 0;
SEXP appendOutput(SEXP value,SEXP result) {
  R_result_pos++;
  if(LENGTH(result) > R_result_pos)
    SET_VECTOR_ELT(result,R_result_pos,value);
  return result;
}


// SWIGINTERN
int
  SWIG_AsVal_long (SEXP obj, long *val)
  {
    if (val) *val = Rf_asInteger(obj);
    return SWIG_OK;
  }


// SWIGINTERN
int
  SWIG_AsVal_int (SEXP obj, int *val)
  {
    long v;
    int res = SWIG_AsVal_long (obj, &v);
    if (SWIG_IsOK(res)) {
      if ((v < INT_MIN || v > INT_MAX)) {
        return SWIG_OverflowError;
      } else {
        if (val) *val = static_cast< int >(v);
      }
    }
    return res;
  }




/* ********************
 * this is a hack to handle multiple values outrput
 * return value is a vector initailized by swig in typemap out
 * R_result_pos is the index in the vector : it must be set to 0 in typemap(out)
 * typemap(argout) cannot be used without typemap(ou) in a function
 */
/*    end of hack */

//SWIGEXPORT SEXP

// [[Rcpp::export]]
RcppExport SEXP swig_omp(SEXP X, SEXP D, SEXP return_reg_path, SEXP given_L, SEXP L, SEXP given_eps, SEXP eps, SEXP given_Lambda, SEXP Lambda, SEXP numThreads)
{
  SpMatrix< double > *result = 0 ;
  Matrix< double > *arg1 = (Matrix< double > *) 0 ;
  Matrix< double > *arg2 = (Matrix< double > *) 0 ;
  Matrix< double > **arg3 = (Matrix< double > **) 0 ;
  bool arg4 ;
  bool arg5 ;
  Vector< int > *arg6 = (Vector< int > *) 0 ;
  bool arg7 ;
  Vector< double > *arg8 = (Vector< double > *) 0 ;
  bool arg9 ;
  Vector< double > *arg10 = (Vector< double > *) 0 ;
  int arg11 ;
  Matrix< double > *data_temp3 ;
  unsigned int r_nprotect = 0;
  SEXP r_ans = R_NilValue ;
  VMAXTYPE r_vmax = vmaxget() ;
  SEXP R_OutputValues;

  {
    arg3 = &data_temp3;
  }
  {
    /*@SWIG:R_typemaps.i,113,map_matrix@*/
    SEXP rmat=X;
    SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
    if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
    {
      /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
      {
        Rf_error("Expected double dense matrix as argument %d",1); return R_NilValue;
      };
    }
    arg1 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

    /*@SWIG@*/
  }
  {
    /*@SWIG:R_typemaps.i,113,map_matrix@*/
    SEXP rmat=D;
    SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
    if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
    {
      /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
      {
        Rf_error("Expected double dense matrix as argument %d",2); return R_NilValue;
      };
    }
    arg2 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

    /*@SWIG@*/
  }
  arg4 = LOGICAL(return_reg_path)[0] ? true : false;
  arg5 = LOGICAL(given_L)[0] ? true : false;
  {
    SEXP rvec=L;
    if (TYPEOF(rvec) != INTSXP || ! Rf_isVector(rvec))
    {
      /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
      {
        Rf_error("Expected int Vector as argument %d",6); return R_NilValue;
      };
    }

    arg6 = new Vector<int>((int*) INTEGER(rvec), LENGTH(rvec));
  }
  arg7 = LOGICAL(given_eps)[0] ? true : false;
  {
    SEXP rvec=eps;
    if (TYPEOF(rvec) != REALSXP || ! Rf_isVector(rvec))
    {
      /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
      {
        Rf_error("Expected double Vector as argument %d",8); return R_NilValue;
      };
    }

    arg8 = new Vector<double>((double*) REAL(rvec), LENGTH(rvec));
  }
  arg9 = LOGICAL(given_Lambda)[0] ? true : false;
  {
    SEXP rvec=Lambda;
    if (TYPEOF(rvec) != REALSXP || ! Rf_isVector(rvec))
    {
      /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
      {
        Rf_error("Expected double Vector as argument %d",10); return R_NilValue;
      };
    }

    arg10 = new Vector<double>((double*) REAL(rvec), LENGTH(rvec));
  }
  arg11 = static_cast< int >(INTEGER(numThreads)[0]);
  try {
    result = (SpMatrix< double > *)_omp< double >(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11);
  }
  catch(char const *_e) {
    Rf_error("Runtime Error %s",_e);
    return R_NilValue;
  }

  {
    R_result_pos = 0;
    R_len_t m = result->m();
    R_len_t n = result->n();
    R_len_t nzmax = result->nzmax();
    SEXP indptr, indices, vdata, dims, output;
    PROTECT(indptr = Rf_allocVector(INTSXP,n + 1));
    PROTECT(indices = Rf_allocVector(INTSXP,nzmax));
    PROTECT(vdata = Rf_allocVector(REALSXP,nzmax));
    PROTECT(dims = Rf_allocVector(VECSXP,2));
    SET_VECTOR_ELT(dims,0,Rf_ScalarInteger(m));
    SET_VECTOR_ELT(dims,1,Rf_ScalarInteger(n));

    double *xdata = result->v();
    memcpy(REAL(vdata),xdata,nzmax * sizeof(double));
    int *pB = result->pB();
    int *r = result->r();
    memcpy(INTEGER(indices),r,nzmax * sizeof(int));
    memcpy(INTEGER(indptr),pB,(n + 1) * sizeof(int));

    PROTECT(output = Rf_allocVector(VECSXP,4));
    SET_VECTOR_ELT(output,0,indptr);
    SET_VECTOR_ELT(output,1,indices);
    SET_VECTOR_ELT(output,2,vdata);
    SET_VECTOR_ELT(output,3,dims);
    delete result;
    r_ans = output;
    UNPROTECT(5);
  }
  Rf_protect(r_ans);
  Rf_protect(R_OutputValues = Rf_allocVector(VECSXP,2));
  r_nprotect += 2;
  SET_VECTOR_ELT(R_OutputValues, 0, r_ans);
  r_ans = R_OutputValues;
  {
    //# test argout
    if(data_temp3 != NULL) {
      R_len_t m = data_temp3->m();
      R_len_t n = data_temp3->n();
      double *data = data_temp3->rawX();
      SEXP rmat;
      PROTECT(rmat = Rf_allocMatrix(REALSXP,m,n));
      if (!rmat) {
        Rf_error("Cannot alloc R matrix for arg %d",3); return R_NilValue;
      };
      double *rdata = (double *)REAL(rmat);
      memcpy(rdata,data,m * n * sizeof(double));
      delete data_temp3;
      r_ans = appendOutput(rmat, R_OutputValues);;
      UNPROTECT(1);
    }
  }
  {
    delete arg1;
  }
  {
    delete arg2;
  }



  {
    delete arg6;
  }

  {
    delete arg8;
  }

  {
    delete arg10;
  }

  vmaxset(r_vmax);
  if(r_nprotect)  Rf_unprotect(r_nprotect);

  return r_ans;
}


//// Lasso Related ////

//SWIGEXPORT SEXP

// [[Rcpp::export]]
RcppExport SEXP  swig_lassoD ( SEXP X, SEXP D, SEXP return_reg_path, SEXP L, SEXP constraint, SEXP lambda2, SEXP mode, SEXP pos, SEXP ols, SEXP numThreads, SEXP max_length_path, SEXP verbose, SEXP cholevsky)
  {
    SpMatrix< double > *result = 0 ;
    Matrix< double > *arg1 = (Matrix< double > *) 0 ;
    Matrix< double > *arg2 = (Matrix< double > *) 0 ;
    Matrix< double > **arg3 = (Matrix< double > **) 0 ;
    bool arg4 ;
    int arg5 ;
    double arg6 ;
    double arg7 ;
    constraint_type arg8 ;
    bool arg9 ;
    bool arg10 ;
    int arg11 ;
    int arg12 ;
    bool arg13 ;
    bool arg14 ;
    Matrix< double > *data_temp3 ;
    int val8 ;
    int ecode8 = 0 ;
    unsigned int r_nprotect = 0;
    SEXP r_ans = R_NilValue ;
    VMAXTYPE r_vmax = vmaxget() ;
    SEXP R_OutputValues;

    {
      arg3 = &data_temp3;
    }
    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=X;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",1); return R_NilValue;
        };
      }
      arg1 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=D;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",2); return R_NilValue;
        };
      }
      arg2 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    arg4 = LOGICAL(return_reg_path)[0] ? true : false;
    arg5 = static_cast< int >(INTEGER(L)[0]);
    arg6 = static_cast< double >(REAL(constraint)[0]);
    arg7 = static_cast< double >(REAL(lambda2)[0]);
    ecode8 = SWIG_AsVal_int(mode, &val8);
    if (!SWIG_IsOK(ecode8)) {
      // SWIG_exception_fail(SWIG_ArgError(ecode8), "in method '" "lassoD" "', argument " "8"" of type '" "constraint_type""'");
      Rf_error("In method lassoD, argument mode can only be L1COEFFS, L2ERROR, PENALTY"); return R_NilValue;
    }
    arg8 = static_cast< constraint_type >(val8);
    arg9 = LOGICAL(pos)[0] ? true : false;
    arg10 = LOGICAL(ols)[0] ? true : false;
    arg11 = static_cast< int >(INTEGER(numThreads)[0]);
    arg12 = static_cast< int >(INTEGER(max_length_path)[0]);
    arg13 = LOGICAL(verbose)[0] ? true : false;
    arg14 = LOGICAL(cholevsky)[0] ? true : false;
    try {
      result = (SpMatrix< double > *)_lassoD< double >(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14);
    }
    catch(char const *_e) {
      Rf_error("Runtime Error %s",_e);
      return R_NilValue;

    }

    {
      R_result_pos = 0;
      R_len_t m = result->m();
      R_len_t n = result->n();
      R_len_t nzmax = result->nzmax();
      SEXP indptr, indices, vdata, dims, output;
      PROTECT(indptr = Rf_allocVector(INTSXP,n + 1));
      PROTECT(indices = Rf_allocVector(INTSXP,nzmax));
      PROTECT(vdata = Rf_allocVector(REALSXP,nzmax));
      PROTECT(dims = Rf_allocVector(VECSXP,2));
      SET_VECTOR_ELT(dims,0,Rf_ScalarInteger(m));
      SET_VECTOR_ELT(dims,1,Rf_ScalarInteger(n));

      double *xdata = result->v();
      memcpy(REAL(vdata),xdata,nzmax * sizeof(double));
      int *pB = result->pB();
      int *r = result->r();
      memcpy(INTEGER(indices),r,nzmax * sizeof(int));
      memcpy(INTEGER(indptr),pB,(n + 1) * sizeof(int));

      PROTECT(output = Rf_allocVector(VECSXP,4));
      SET_VECTOR_ELT(output,0,indptr);
      SET_VECTOR_ELT(output,1,indices);
      SET_VECTOR_ELT(output,2,vdata);
      SET_VECTOR_ELT(output,3,dims);
      delete result;
      r_ans = output;
      UNPROTECT(5);
    }
    Rf_protect(r_ans);
    Rf_protect(R_OutputValues = Rf_allocVector(VECSXP,2));
    r_nprotect += 2;
    SET_VECTOR_ELT(R_OutputValues, 0, r_ans);
    r_ans = R_OutputValues;
    {
      //# test argout
      if(data_temp3 != NULL) {
        R_len_t m = data_temp3->m();
        R_len_t n = data_temp3->n();
        double *data = data_temp3->rawX();
        SEXP rmat;
        PROTECT(rmat = Rf_allocMatrix(REALSXP,m,n));
        if (!rmat) {
          Rf_error("Cannot alloc R matrix for arg %d",3); return R_NilValue;
        };
        double *rdata = (double *)REAL(rmat);
        memcpy(rdata,data,m * n * sizeof(double));
        delete data_temp3;
        r_ans = appendOutput(rmat, R_OutputValues);;
        UNPROTECT(1);
      }
    }
    {
      delete arg1;
    }
    {
      delete arg2;
    }



    vmaxset(r_vmax);
    if(r_nprotect)  Rf_unprotect(r_nprotect);

    return r_ans;
  }


//SWIGEXPORT SEXP
// [[Rcpp::export]]
RcppExport SEXP  swig_lassoQq ( SEXP X, SEXP Q, SEXP q, SEXP return_reg_path, SEXP L, SEXP constraint, SEXP lambda2, SEXP mode, SEXP pos, SEXP ols, SEXP numThreads, SEXP max_length_path, SEXP verbose, SEXP cholevsky)
  {
    SpMatrix< double > *result = 0 ;
    Matrix< double > *arg1 = (Matrix< double > *) 0 ;
    Matrix< double > *arg2 = (Matrix< double > *) 0 ;
    Matrix< double > *arg3 = (Matrix< double > *) 0 ;
    Matrix< double > **arg4 = (Matrix< double > **) 0 ;
    bool arg5 ;
    int arg6 ;
    double arg7 ;
    double arg8 ;
    constraint_type arg9 ;
    bool arg10 ;
    bool arg11 ;
    int arg12 ;
    int arg13 ;
    bool arg14 ;
    bool arg15 ;
    Matrix< double > *data_temp4 ;
    int val9 ;
    int ecode9 = 0 ;
    unsigned int r_nprotect = 0;
    SEXP r_ans = R_NilValue ;
    VMAXTYPE r_vmax = vmaxget() ;
    SEXP R_OutputValues;

    {
      arg4 = &data_temp4;
    }
    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=X;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",1); return R_NilValue;
        };
      }
      arg1 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=Q;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",2); return R_NilValue;
        };
      }
      arg2 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=q;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",3); return R_NilValue;
        };
      }
      arg3 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    arg5 = LOGICAL(return_reg_path)[0] ? true : false;
    arg6 = static_cast< int >(INTEGER(L)[0]);
    arg7 = static_cast< double >(REAL(constraint)[0]);
    arg8 = static_cast< double >(REAL(lambda2)[0]);
    ecode9 = SWIG_AsVal_int(mode, &val9);
    if (!SWIG_IsOK(ecode9)) {
      // SWIG_exception_fail(SWIG_ArgError(ecode9), "in method '" "lassoQq" "', argument " "9"" of type '" "constraint_type""'");
      Rf_error("In method lassoQq, argument mode can only be L1COEFFS, L2ERROR, PENALTY"); return R_NilValue;

    }
    arg9 = static_cast< constraint_type >(val9);
    arg10 = LOGICAL(pos)[0] ? true : false;
    arg11 = LOGICAL(ols)[0] ? true : false;
    arg12 = static_cast< int >(INTEGER(numThreads)[0]);
    arg13 = static_cast< int >(INTEGER(max_length_path)[0]);
    arg14 = LOGICAL(verbose)[0] ? true : false;
    arg15 = LOGICAL(cholevsky)[0] ? true : false;
    try {
      result = (SpMatrix< double > *)_lassoQq< double >(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15);
    }
    catch(char const *_e) {
      Rf_error("Runtime Error %s",_e);
      return R_NilValue;

    }

    {
      R_result_pos = 0;
      R_len_t m = result->m();
      R_len_t n = result->n();
      R_len_t nzmax = result->nzmax();
      SEXP indptr, indices, vdata, dims, output;
      PROTECT(indptr = Rf_allocVector(INTSXP,n + 1));
      PROTECT(indices = Rf_allocVector(INTSXP,nzmax));
      PROTECT(vdata = Rf_allocVector(REALSXP,nzmax));
      PROTECT(dims = Rf_allocVector(VECSXP,2));
      SET_VECTOR_ELT(dims,0,Rf_ScalarInteger(m));
      SET_VECTOR_ELT(dims,1,Rf_ScalarInteger(n));

      double *xdata = result->v();
      memcpy(REAL(vdata),xdata,nzmax * sizeof(double));
      int *pB = result->pB();
      int *r = result->r();
      memcpy(INTEGER(indices),r,nzmax * sizeof(int));
      memcpy(INTEGER(indptr),pB,(n + 1) * sizeof(int));

      PROTECT(output = Rf_allocVector(VECSXP,4));
      SET_VECTOR_ELT(output,0,indptr);
      SET_VECTOR_ELT(output,1,indices);
      SET_VECTOR_ELT(output,2,vdata);
      SET_VECTOR_ELT(output,3,dims);
      delete result;
      r_ans = output;
      UNPROTECT(5);
    }
    Rf_protect(r_ans);
    Rf_protect(R_OutputValues = Rf_allocVector(VECSXP,2));
    r_nprotect += 2;
    SET_VECTOR_ELT(R_OutputValues, 0, r_ans);
    r_ans = R_OutputValues;
    {
      //# test argout
      if(data_temp4 != NULL) {
        R_len_t m = data_temp4->m();
        R_len_t n = data_temp4->n();
        double *data = data_temp4->rawX();
        SEXP rmat;
        PROTECT(rmat = Rf_allocMatrix(REALSXP,m,n));
        if (!rmat) {
          Rf_error("Cannot alloc R matrix for arg %d",4); return R_NilValue;
        };
        double *rdata = (double *)REAL(rmat);
        memcpy(rdata,data,m * n * sizeof(double));
        delete data_temp4;
        r_ans = appendOutput(rmat, R_OutputValues);;
        UNPROTECT(1);
      }
    }
    {
      delete arg1;
    }
    {
      delete arg2;
    }
    {
      delete arg3;
    }












    vmaxset(r_vmax);
    if(r_nprotect)  Rf_unprotect(r_nprotect);

    return r_ans;
  }


//SWIGEXPORT SEXP

// [[Rcpp::export]]
RcppExport SEXP  swig_lassoMask ( SEXP X, SEXP D, SEXP B, SEXP L, SEXP constraint, SEXP lambda2, SEXP mode, SEXP pos, SEXP numThreads, SEXP verbose)
  {
    SpMatrix< double > *result = 0 ;
    Matrix< double > *arg1 = (Matrix< double > *) 0 ;
    Matrix< double > *arg2 = (Matrix< double > *) 0 ;
    Matrix< bool > *arg3 = (Matrix< bool > *) 0 ;
    int arg4 ;
    double arg5 ;
    double arg6 ;
    constraint_type arg7 ;
    bool arg8 ;
    int arg9 ;
    bool arg10 ;
    int val7 ;
    int ecode7 = 0 ;
    unsigned int r_nprotect = 0;
    SEXP r_ans = R_NilValue ;
    VMAXTYPE r_vmax = vmaxget() ;

    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=X;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",1); return R_NilValue;
        };
      }
      arg1 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=D;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",2); return R_NilValue;
        };
      }
      arg2 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    {
      SEXP rmat=B;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != LGLSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected bool dense matrix as argument %d",3); return R_NilValue;
        };
      }
      arg3 = new Matrix<bool> (Rf_nrows(rmat),Rf_ncols(rmat));
      int *pi = (int *)LOGICAL(rmat);
      bool *po = arg3->rawX();
      for(int i =0;i < Rf_nrows(rmat) * Rf_ncols(rmat);i++)
        *po++ = (bool) *pi++;

    }
    arg4 = static_cast< int >(INTEGER(L)[0]);
    arg5 = static_cast< double >(REAL(constraint)[0]);
    arg6 = static_cast< double >(REAL(lambda2)[0]);
    ecode7 = SWIG_AsVal_int(mode, &val7);
    if (!SWIG_IsOK(ecode7)) {
      // SWIG_exception_fail(SWIG_ArgError(ecode7), "in method '" "lassoMask" "', argument " "7"" of type '" "constraint_type""'");
      Rf_error("In method lassoMask, argument mode can only be L1COEFFS, L2ERROR, PENALTY"); return R_NilValue;

    }
    arg7 = static_cast< constraint_type >(val7);
    arg8 = LOGICAL(pos)[0] ? true : false;
    arg9 = static_cast< int >(INTEGER(numThreads)[0]);
    arg10 = LOGICAL(verbose)[0] ? true : false;
    try {
      result = (SpMatrix< double > *)_lassoMask< double >(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10);
    }
    catch(char const *_e) {
      Rf_error("Runtime Error %s",_e);
      return R_NilValue;

    }

    {
      R_result_pos = 0;
      R_len_t m = result->m();
      R_len_t n = result->n();
      R_len_t nzmax = result->nzmax();
      SEXP indptr, indices, vdata, dims, output;
      PROTECT(indptr = Rf_allocVector(INTSXP,n + 1));
      PROTECT(indices = Rf_allocVector(INTSXP,nzmax));
      PROTECT(vdata = Rf_allocVector(REALSXP,nzmax));
      PROTECT(dims = Rf_allocVector(VECSXP,2));
      SET_VECTOR_ELT(dims,0,Rf_ScalarInteger(m));
      SET_VECTOR_ELT(dims,1,Rf_ScalarInteger(n));

      double *xdata = result->v();
      memcpy(REAL(vdata),xdata,nzmax * sizeof(double));
      int *pB = result->pB();
      int *r = result->r();
      memcpy(INTEGER(indices),r,nzmax * sizeof(int));
      memcpy(INTEGER(indptr),pB,(n + 1) * sizeof(int));

      PROTECT(output = Rf_allocVector(VECSXP,4));
      SET_VECTOR_ELT(output,0,indptr);
      SET_VECTOR_ELT(output,1,indices);
      SET_VECTOR_ELT(output,2,vdata);
      SET_VECTOR_ELT(output,3,dims);
      delete result;
      r_ans = output;
      UNPROTECT(5);
    }
    {
      delete arg1;
    }
    {
      delete arg2;
    }
    {
      delete arg3;
    }







    vmaxset(r_vmax);
    if(r_nprotect)  Rf_unprotect(r_nprotect);

    return r_ans;
  }


//SWIGEXPORT SEXP

// [[Rcpp::export]]
RcppExport SEXP   swig_lassoWeighted ( SEXP X, SEXP D, SEXP W, SEXP L, SEXP constraint, SEXP mode, SEXP pos, SEXP numThreads, SEXP verbose)
  {
    SpMatrix< double > *result = 0 ;
    Matrix< double > *arg1 = (Matrix< double > *) 0 ;
    Matrix< double > *arg2 = (Matrix< double > *) 0 ;
    Matrix< double > *arg3 = (Matrix< double > *) 0 ;
    int arg4 ;
    double arg5 ;
    constraint_type arg6 ;
    bool arg7 ;
    int arg8 ;
    bool arg9 ;
    int val6 ;
    int ecode6 = 0 ;
    unsigned int r_nprotect = 0;
    SEXP r_ans = R_NilValue ;
    VMAXTYPE r_vmax = vmaxget() ;

    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=X;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",1); return R_NilValue;
        };
      }
      arg1 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=D;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",2); return R_NilValue;
        };
      }
      arg2 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    {
      /*@SWIG:R_typemaps.i,113,map_matrix@*/
      SEXP rmat=W;
      SEXP dims = Rf_getAttrib(rmat,Rf_install("dim"));
      if (TYPEOF(rmat) != REALSXP || LENGTH(dims) != 2)
      {
        /*SG_ERROR("Expected Double Vector as argument %d\n", m_rhs_counter);*/
        {
          Rf_error("Expected double dense matrix as argument %d",3); return R_NilValue;
        };
      }
      arg3 = new Matrix<double> ((double *)REAL(rmat),Rf_nrows(rmat),Rf_ncols(rmat));

      /*@SWIG@*/
    }
    arg4 = static_cast< int >(INTEGER(L)[0]);
    arg5 = static_cast< double >(REAL(constraint)[0]);
    ecode6 = SWIG_AsVal_int(mode, &val6);
    if (!SWIG_IsOK(ecode6)) {
      // SWIG_exception_fail(SWIG_ArgError(ecode6), "in method '" "lassoWeighted" "', argument " "6"" of type '" "constraint_type""'");
      Rf_error("In method lassoWeighted, argument mode can only be L1COEFFS, L2ERROR, PENALTY"); return R_NilValue;

    }
    arg6 = static_cast< constraint_type >(val6);
    arg7 = LOGICAL(pos)[0] ? true : false;
    arg8 = static_cast< int >(INTEGER(numThreads)[0]);
    arg9 = LOGICAL(verbose)[0] ? true : false;
    try {
      result = (SpMatrix< double > *)_lassoWeighted< double >(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9);
    }
    catch(char const *_e) {
      Rf_error("Runtime Error %s",_e);
      return R_NilValue;

    }

    {
      R_result_pos = 0;
      R_len_t m = result->m();
      R_len_t n = result->n();
      R_len_t nzmax = result->nzmax();
      SEXP indptr, indices, vdata, dims, output;
      PROTECT(indptr = Rf_allocVector(INTSXP,n + 1));
      PROTECT(indices = Rf_allocVector(INTSXP,nzmax));
      PROTECT(vdata = Rf_allocVector(REALSXP,nzmax));
      PROTECT(dims = Rf_allocVector(VECSXP,2));
      SET_VECTOR_ELT(dims,0,Rf_ScalarInteger(m));
      SET_VECTOR_ELT(dims,1,Rf_ScalarInteger(n));

      double *xdata = result->v();
      memcpy(REAL(vdata),xdata,nzmax * sizeof(double));
      int *pB = result->pB();
      int *r = result->r();
      memcpy(INTEGER(indices),r,nzmax * sizeof(int));
      memcpy(INTEGER(indptr),pB,(n + 1) * sizeof(int));

      PROTECT(output = Rf_allocVector(VECSXP,4));
      SET_VECTOR_ELT(output,0,indptr);
      SET_VECTOR_ELT(output,1,indices);
      SET_VECTOR_ELT(output,2,vdata);
      SET_VECTOR_ELT(output,3,dims);
      delete result;
      r_ans = output;
      UNPROTECT(5);
    }
    {
      delete arg1;
    }
    {
      delete arg2;
    }
    {
      delete arg3;
    }






    vmaxset(r_vmax);
    if(r_nprotect)  Rf_unprotect(r_nprotect);

    return r_ans;
  }


