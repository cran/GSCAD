# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

swig_omp <- function(X, D, return_reg_path, given_L, L, given_eps, eps, given_Lambda, Lambda, numThreads) {
    .Call('_GSCAD_swig_omp', PACKAGE = 'GSCAD', X, D, return_reg_path, given_L, L, given_eps, eps, given_Lambda, Lambda, numThreads)
}

swig_lassoD <- function(X, D, return_reg_path, L, constraint, lambda2, mode, pos, ols, numThreads, max_length_path, verbose, cholevsky) {
    .Call('_GSCAD_swig_lassoD', PACKAGE = 'GSCAD', X, D, return_reg_path, L, constraint, lambda2, mode, pos, ols, numThreads, max_length_path, verbose, cholevsky)
}

swig_lassoQq <- function(X, Q, q, return_reg_path, L, constraint, lambda2, mode, pos, ols, numThreads, max_length_path, verbose, cholevsky) {
    .Call('_GSCAD_swig_lassoQq', PACKAGE = 'GSCAD', X, Q, q, return_reg_path, L, constraint, lambda2, mode, pos, ols, numThreads, max_length_path, verbose, cholevsky)
}

swig_lassoMask <- function(X, D, B, L, constraint, lambda2, mode, pos, numThreads, verbose) {
    .Call('_GSCAD_swig_lassoMask', PACKAGE = 'GSCAD', X, D, B, L, constraint, lambda2, mode, pos, numThreads, verbose)
}

swig_lassoWeighted <- function(X, D, W, L, constraint, mode, pos, numThreads, verbose) {
    .Call('_GSCAD_swig_lassoWeighted', PACKAGE = 'GSCAD', X, D, W, L, constraint, mode, pos, numThreads, verbose)
}

D2optim_mat <- function(A, lambda, c, rho) {
    .Call('_GSCAD_D2optim_mat', PACKAGE = 'GSCAD', A, lambda, c, rho)
}

pgscad <- function(D, c, lambda) {
    .Call('_GSCAD_pgscad', PACKAGE = 'GSCAD', D, c, lambda)
}

imageSplit <- function(A, m, n, stepsize) {
    .Call('_GSCAD_imageSplit', PACKAGE = 'GSCAD', A, m, n, stepsize)
}

imageReCon <- function(B, mm, nn, m, n, stepsize) {
    .Call('_GSCAD_imageReCon', PACKAGE = 'GSCAD', B, mm, nn, m, n, stepsize)
}

