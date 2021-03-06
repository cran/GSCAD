% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSCAD_DL.R
\name{gscad.DL}
\alias{gscad.DL}
\alias{gscad.DLmask}
\title{Learn dictionary under GSCAD regularization}
\usage{
gscad.DL(Y, D0 = NULL, p0, sigma = 0, c = 3.7, lambda = 0.05,
  maxrun = 20, maxrun_ADMM = 20, err_bnd = 1e-06, err_bnd2 = 1e-04,
  rho = 16, cor_bnd = 1, L = 30, LassoMode = 1, LassoLambda = NULL)

gscad.DLmask(Y, Mask, D0 = NULL, p0, sigma = 0, c = 3.7, lambda = 0.05,
  maxrun = 20, maxrun_ADMM = 20, err_bnd = 1e-06, err_bnd2 = 1e-06,
  rho = 16, LassoMode = 1, LassoLambda = NULL)
}
\arguments{
\item{Y}{An m by n matrix in Y=DA. Each column of Y is a sample of size m
(usually a vectorized sqrt(m) by sqrt(m) patches).}

\item{D0}{Initial dictionary. If D0 specified, p0 is not needed, otherwise D0
is evaluated as overcompleted DCT basis using function ODCT(m,p0). Either D0
or p0 needs to be specified.}

\item{p0}{Initial size of the dictionary.}

\item{sigma}{Noise level.}

\item{c, lambda}{Parameters for GSCAD.}

\item{maxrun}{(optional) Maximun number of outer iterations to run. Default is 20.}

\item{maxrun_ADMM}{(optional) Maximun number of iterations to run for updating
dictionary using ADMM. Default is 20.}

\item{err_bnd}{(optional) Stopping criterion for iterations. Default is 1e-6.}

\item{err_bnd2}{(optional) Stopping criterion for updating dictionary
UpDic. Default is 1e-4.}

\item{rho}{(optional) Parameter for ADMM. Default is 16.}

\item{cor_bnd}{(optional) When normalize dictionary, checking if the correlation
of any two atoms are above the cor_bnd, one of the atom is removed. Default is 1.}

\item{L}{(optional) This parameter controls the maximum number of non-zero elements
in each column of sparsecolding A.}

\item{LassoMode}{(optional) At the sparse coding stage, the optimization can be
done in three modes L1COEFFS (0), L2ERROR (1), PENALTY(2). Default is 1.}

\item{LassoLambda}{(optional) Tuning parameter for Lasso}

\item{Mask}{{0,1} matrix of the same size as Y
to indicate the location of corrupted pixels.}
}
\value{
The learned dictionary \code{dictionary}
}
\description{
This function learns the dictionary D under the framework of matrix
factorization \eqn{Y=DA}. Y is a given m by n matrix of n samples.
D is an m by p matrix, where p is unkonw as well, and A ia a p by n matrix.
Both D and A are unkonwn. GSCAD regularization is applied to D and lasso
regularization is applied to A.
}
\details{
See \url{https://arxiv.org/abs/1605.07870}
}
\section{Functions}{
\itemize{
\item \code{gscad.DLmask}: Adding Mask {0,1} matrix of the same size as Y to indicate
the location of corrupted pixel
}}

\examples{
I = lena_crop #use a smaller image as an example
## add noise
sigma=20;m=64
I_noise=AddNoise(I,sigma)
## spliting image into patches
Y_nc = ImageSplit(I_noise,sqrt(m),sqrt(m));
mu=colMeans(Y_nc)
Y=Y_nc-rep(mu,each=nrow(Y_nc))
## learning dictionary
\dontrun{
dictionary=gscad.DL(Y,p0=256,sigma=sigma)
}
I=lena_crop
## corrupt 30\% of the image
out_corrupt=AddHoles(I,0.3)
I_corrupt=out_corrupt$corruptedImage
I_mask=out_corrupt$maskImage
## split image
m=64
Y_nc = ImageSplit(I_corrupt,sqrt(m),sqrt(m));
M = ImageSplit(I_mask,sqrt(m),sqrt(m));
mu=colSums(Y_nc*M)/colSums(M)
Y=Y_nc-M*rep(mu,each=nrow(Y_nc))
mask = matrix(as.logical(M),ncol=ncol(M))
## learn dictionary for inpainting, this function is slow
\dontrun{
dic=gscad.DLmask(Y, mask,  p0=100, sigma=1)
}
}
