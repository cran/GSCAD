# The code in this file as well as corresponding cpp code in scr/decomp_swag.cpp
# scr/decomp.h and scr/linalg are based on the SPAMS software.
# Some modifications were made. The original code can be
# found at http://spams-devel.gforge.inria.fr/downloads.html.


`omp.core` = function(X, D, return_reg_path, given_L, L, given_eps, eps, given_Lambda, Lambda, numThreads)
{
  if (inherits(X, "ExternalReference")) X = slot(X,"ref")
  if (inherits(D, "ExternalReference")) D = slot(D,"ref")
  return_reg_path = as.logical(return_reg_path);
  given_L = as.logical(given_L);
  if (inherits(L, "ExternalReference")) L = slot(L,"ref")
  given_eps = as.logical(given_eps);
  if (inherits(eps, "ExternalReference")) eps = slot(eps,"ref")
  given_Lambda = as.logical(given_Lambda);
  if (inherits(Lambda, "ExternalReference")) Lambda = slot(Lambda,"ref")
  numThreads = as.integer(numThreads);

  if(length(numThreads) > 1) {
    warning("using only the first element of numThreads");
  };

  ;ans = .Call('swig_omp', X, D, return_reg_path, given_L, L, given_eps, eps, given_Lambda, Lambda, numThreads, PACKAGE='GSCAD');
  # ans <- new("_p_SpMatrixT_double_t", ref=ans) ;
  ans

}


# Start of lassoD
`lassoD.core` = function(X, D, return_reg_path, L, constraint, lambda2, mode, pos, ols, numThreads, max_length_path, verbose, cholevsky)
{
  if (inherits(X, "ExternalReference")) X = slot(X,"ref")
  if (inherits(D, "ExternalReference")) D = slot(D,"ref")
  return_reg_path = as.logical(return_reg_path);
  L = as.integer(L);

  if(length(L) > 1) {
    warning("using only the first element of L");
  };



  # mode = enumToInteger(mode, "_constraint_type");

  if(length(mode) > 1) {
    warning("using only the first element of mode");
  };

  pos = as.logical(pos);
  ols = as.logical(ols);
  numThreads = as.integer(numThreads);

  if(length(numThreads) > 1) {
    warning("using only the first element of numThreads");
  };

  max_length_path = as.integer(max_length_path);

  if(length(max_length_path) > 1) {
    warning("using only the first element of max_length_path");
  };

  verbose = as.logical(verbose);
  cholevsky = as.logical(cholevsky);
  ;ans = .Call('swig_lassoD', X, D, return_reg_path, L, constraint, lambda2, mode, pos, ols, numThreads, max_length_path, verbose, cholevsky, PACKAGE='GSCAD');
  # ans <- new("_p_SpMatrixT_double_t", ref=ans) ;

  ans

}

# Start of lassoQq

`lassoQq.core` = function(X, Q, q, return_reg_path, L, constraint, lambda2, mode, pos, ols, numThreads, max_length_path, verbose, cholevsky)
{
  if (inherits(X, "ExternalReference")) X = slot(X,"ref")
  if (inherits(Q, "ExternalReference")) Q = slot(Q,"ref")
  if (inherits(q, "ExternalReference")) q = slot(q,"ref")
  return_reg_path = as.logical(return_reg_path);
  L = as.integer(L);

  if(length(L) > 1) {
    warning("using only the first element of L");
  };



  # mode = enumToInteger(mode, "_constraint_type");

  if(length(mode) > 1) {
    warning("using only the first element of mode");
  };

  pos = as.logical(pos);
  ols = as.logical(ols);
  numThreads = as.integer(numThreads);

  if(length(numThreads) > 1) {
    warning("using only the first element of numThreads");
  };

  max_length_path = as.integer(max_length_path);

  if(length(max_length_path) > 1) {
    warning("using only the first element of max_length_path");
  };

  verbose = as.logical(verbose);
  cholevsky = as.logical(cholevsky);
  ;ans = .Call('swig_lassoQq', X, Q, q, return_reg_path, L, constraint, lambda2, mode, pos, ols, numThreads, max_length_path, verbose, cholevsky, PACKAGE='GSCAD');
  # ans <- new("_p_SpMatrixT_double_t", ref=ans) ;

  ans

}

# Start of lassoMask

`lassoMask.core` = function(X, D, B, L, constraint, lambda2, mode, pos, numThreads, verbose)
{
  if (inherits(X, "ExternalReference")) X = slot(X,"ref")
  if (inherits(D, "ExternalReference")) D = slot(D,"ref")
  if (inherits(B, "ExternalReference")) B = slot(B,"ref")
  L = as.integer(L);

  if(length(L) > 1) {
    warning("using only the first element of L");
  };



  # mode = enumToInteger(mode, "_constraint_type");

  if(length(mode) > 1) {
    warning("using only the first element of mode");
  };

  pos = as.logical(pos);
  numThreads = as.integer(numThreads);

  if(length(numThreads) > 1) {
    warning("using only the first element of numThreads");
  };

  verbose = as.logical(verbose);
  ;ans = .Call('swig_lassoMask', X, D, B, L, constraint, lambda2, mode, pos, numThreads, verbose, PACKAGE='GSCAD');
  # ans <- new("_p_SpMatrixT_double_t", ref=ans) ;

  ans

}
# Start of lassoWeighted

`lassoWeighted.core` = function(X, D, W, L, constraint, mode, pos, numThreads, verbose)
{
  if (inherits(X, "ExternalReference")) X = slot(X,"ref")
  if (inherits(D, "ExternalReference")) D = slot(D,"ref")
  if (inherits(W, "ExternalReference")) W = slot(W,"ref")
  L = as.integer(L);

  if(length(L) > 1) {
    warning("using only the first element of L");
  };


  # mode = enumToInteger(mode, "_constraint_type");

  if(length(mode) > 1) {
    warning("using only the first element of mode");
  };

  pos = as.logical(pos);
  numThreads = as.integer(numThreads);

  if(length(numThreads) > 1) {
    warning("using only the first element of numThreads");
  };

  verbose = as.logical(verbose);
  ;ans = .Call('swig_lassoWeighted', X, D, W, L, constraint, mode, pos, numThreads, verbose, PACKAGE='GSCAD');
  # ans <- new("_p_SpMatrixT_double_t", ref=ans) ;

  ans

}



##########################################################

# .verif_enum <- function(arg,ename,msg) {
#   defName = paste(".__E___", ename, sep = "")
#   l = eval.parent(parse(text = sprintf("l <- %s",defName)))
#   if (! (arg %in% names(l)))
#     stop("ERROR : bad enum value ",msg,"\n")
# }

#' @export
.omp <- function(X,D,L = NULL,eps = NULL,lambda1 = NULL,return_reg_path = FALSE, numThreads = -1) {
  path = NULL
  given_L = FALSE
  given_eps = FALSE
  given_lambda1 = FALSE
  if (is.null(L)) {
    L = as.vector(c(0),mode='integer')
  } else {
    given_L = TRUE
    if(length(L) == 1 && ! is.integer(L)) {
      L = as.vector(c(L),mode='integer')
    }
  }
  if (is.null(eps)) {
    eps = as.vector(c(0.),mode='double')
  } else {
    given_eps = TRUE
  }
  if (is.null(lambda1)) {
    lambda1 = as.vector(c(0.),mode='double')
  } else {
    given_lambda1 = TRUE
  }

  #  if(! is.vector(eps)) {
  #    eps = as.vector(c(eps),mode='double')
  #  }
  x = omp.core(X,D,return_reg_path,given_L,L,given_eps,eps,given_lambda1,lambda1, numThreads)
  if(return_reg_path) {
    path = x[[2]]
  }
  indptr = x[[1]][[1]]
  indices = x[[1]][[2]]
  data = x[[1]][[3]]
  shape = x[[1]][[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  if (return_reg_path)
    return (list(alpha,path))
  else
    return(alpha)
}

#' @export
.lasso <- function(X,D= NULL,Q = NULL,q = NULL,return_reg_path = FALSE,L= -1,lambda1= NULL,lambda2= 0.,
                        mode= 2,pos= FALSE,ols= FALSE,numThreads= -1,
                        max_length_path= -1,verbose=FALSE,cholesky= FALSE) {
  #  require('Matrix')
  # Note : 'L' and 'max_length_path' default to -1 so that their effective default values
  # will be set in spams.h

  if (! is.null(Q)) {
    if (is.null(q)) {
      stop("ERROR lasso : q is needed when Q is given\n")
    }
  } else {
    if(is.null(D)) {
      stop("ERROR lasso : you must give D or Q and q\n")
    }
  }
  if(is.null(lambda1)) {
    stop("ERROR lasso : lambda1 must be defined\n")
  }
  # .verif_enum(mode,'constraint_type','mode in Lasso')
  path = NULL
  x = NULL
  if(! is.null(q)) {
    ##    x = do.call(spams_wrap.lassoQq,c(list(X,D,q,0,return_reg_path),params))
    x = lassoQq.core(X,Q,q,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,max_length_path,verbose,cholesky)
    ##    x = .mycall('lassoQq',c('X','D','q',0,'return_reg_path',params))
  } else {
    ##    x = do.call(spams_wrap.lassoD,c(list(X,D,0,return_reg_path),params))
    x = lassoD.core(X,D,return_reg_path,L,lambda1,lambda2,mode,pos,ols,numThreads,max_length_path,verbose,cholesky)
    ##    x = .mycall('lassoD',c('X','D',0,'return_reg_path',params))
  }
  if(return_reg_path) {
    path = x[[2]]
  }
  indptr = x[[1]][[1]]
  indices = x[[1]][[2]]
  data = x[[1]][[3]]
  shape = x[[1]][[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  if (return_reg_path)
    return (list(alpha,path))
  else
    return(alpha)

}
#' @export
.lassoMask <- function(X,D,B,L= -1,lambda1= NULL,lambda2= 0.,
                            mode= 'PENALTY',pos= FALSE,numThreads= -1,
                            verbose=FALSE) {
  #  require('Matrix')
  # Note : 'L' and 'max_length_path' default to -1 so that their effective default values
  # will be set in spams.h

  if(is.null(lambda1)) {
    stop("ERROR lassoMask : lambda1 must be defined\n")
  }
  # .verif_enum(mode,'constraint_type','mode in Lasso')
  x = lassoMask.core(X,D,B,L,lambda1,lambda2,mode,pos,numThreads,verbose)
  indptr = x[[1]]
  indices = x[[2]]
  data = x[[3]]
  shape = x[[4]]
  cat("LASSO : ", length(shape),"\n")
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  return(alpha)
}
#' @export
.lassoWeighted <- function(X,D,W,L= -1,lambda1= NULL,
                                mode= 'PENALTY',pos= FALSE,numThreads= -1,
                                verbose=FALSE) {
  #  require('Matrix')
  # Note : 'L' and 'max_length_path' default to -1 so that their effective default values
  # will be set in spams.h

  if(is.null(lambda1)) {
    stop("ERROR lassoWeighted : lambda1 must be defined\n")
  }
  # .verif_enum(mode,'constraint_type','mode in Lasso')
  x = lassoWeighted.core(X,D,W,L,lambda1,mode,pos,numThreads,verbose)
  indptr = x[[1]]
  indices = x[[2]]
  data = x[[3]]
  shape = x[[4]]
  alpha = sparseMatrix(i = indices, p = indptr, x = data,dims = shape, index1 = FALSE)
  return(alpha)
}
