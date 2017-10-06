#' Use a given dictionary D to denoise image.
#'
#' The noisy image is split into sqrt(m) by sqrt(m)
#' patches. Each patch is vectorized into a column of matrix Y.
#' Using the given D, the sparse coding A_hat in Y=DA is obtained.
#' Then Y_denoise=DA_hat. The final
#' denoised image is reconstruncted on the denoised patches.
#'
#' See \url{https://arxiv.org/abs/1605.07870}
#'
#' @param I_noise The image to be denoised. In form of matrix.
#' @param sigma Noise level.
#' @param D D is the dictionary used in Y=DA to denoise.
#' @param stepsize (optional) The stepsize when splicting the image. Default is 1
#' @return The denoised image in for of a matrix.
#'
#' @examples
#' I = lena_crop #use a smaller image as an example
#' ## add noise
#' sigma=20;
#' I_noise=AddNoise(I,sigma)
#' ## use ODCT dictionary
#' D0=ODCT(64,100)
#' ## denoise
#' \dontrun{
#' I_clean=denoiseImage(I_noise,D0,sigma)
#' }

#' @export

denoiseImage <- function(I_noise,D,sigma,stepsize=1)
{

  m=nrow(D)

  ## get patches
  Y_nc = ImageSplit(I_noise,sqrt(m),sqrt(m),stepsize);
  mu=colMeans(Y_nc)
  Y=Y_nc-rep(mu,each=nrow(Y_nc))
  n=ncol(Y)

  ## normalize dictionary using L2 norm
  dic_n=D/rep(sqrt(colSums(D^2)),each=nrow(D))

  ## final estimate of the spares coding
  eps = qchisq(0.9,m)*(sigma/255)^2

  alpha_rcn=.omp(Y,dic_n,eps=eps)

  ## reconstruct image by Y_hat=dic_hat * A_hat
  Y_rcn=as.matrix(dic_n %*% alpha_rcn)
  Y_rcn=Y_rcn+rep(mu,each=nrow(Y))
  Image_fit=ImageReCon(Y_rcn,nrow(I_noise),ncol(I_noise),sqrt(m),sqrt(m),stepsize)

  Image_fit
}


#' Use a given dictionary D to denoise image.
#'
#' Denoise Y give D in Y=DA.
#'
#' See \url{https://arxiv.org/abs/1605.07870}

#' @param Y Each column of Y is a vectorized image patch to be denoised.
#' @param sigma Noise level.
#' @param D D is the dictionary used in Y=DA to denoise.
#' @return The denoised matrix Y.

#' @examples
#' I = lena_crop #use a smaller image as an example
#' ## add noise
#' sigma=20;
#' I_noise=AddNoise(I,sigma)
#' ## spliting image into patches
#' m=64;
#' Y_nc = ImageSplit(I_noise,sqrt(m),sqrt(m));
#' mu=colMeans(Y_nc)
#' Y=Y_nc-rep(mu,each=nrow(Y_nc))
#' ## use ODCT dictionary
#' D0=ODCT(64,100)
#' ## denoise
#' \dontrun{
#' Y_denoise=denoise(Y,D0,sigma)
#' }
#' @export
denoise <- function(Y,D,sigma)
{

  m=nrow(D)
  n=ncol(Y)

  ## normalize dictionary using L2 norm
  dic_n=D/rep(sqrt(colSums(D^2)),each=nrow(D))

  ## final estimate of the spares coding
  eps = qchisq(0.9,m)*(sigma/255)^2
  alpha_rcn=.omp(Y,dic_n,eps=eps)

  ## reconstruct image by Y_hat=dic_hat * A_hat
  Y_rcn=as.matrix(dic_n %*% alpha_rcn)
  Y_rcn
}


######################### gscad.Denoise ###############################
#' Use GSCAD to denoise image
#'
#' Use the GSCAD method under dictionary learning framework to learn a proper
#' sized dictionary and denoise image. The noisy image is split into m by m
#' patches and a dictionary with each atom of size m is learned. The final
#' denoised image is reconstruncted on the denoised patches.
#'
#' See \url{https://arxiv.org/abs/1605.07870}
#'
#' @param I_noise The image to be denoised. In form of matrix.
#' @param sigma Noise level.
#' @param D0 Initial dictionary. If D0 specified, m and p0 are not needed, otherwise D0
#' is evaluated as overcompleted DCT basis using function ODCT(m,p0). Either D0
#' or (m,p0) needs to be specified.
#' @param m The size of the small patches to be split.
#' @param p0 Initial size of the dictionary.
#' @param c,lambda Parameters for GSCAD.
#' @param maxrun (optional) Maximun number of outer iterations to run. Default is 20.
#' @param maxrun_ADMM (optional) Maximun number of iterations to run for updating
#' dictionary using ADMM. Default is 20.
#' @param err_bnd (optional) Stopping criterion for iterations. Default is 1e-4.
#' @param err_bnd2 (optional) Stopping criterion for updating dictionary
#' @param rho (optional) Parameter for ADMM. Default is 16.
#' @param cor_bnd (optional) When normalize dictionary, checking if the correlation
#' of any two atoms are above the cor_bnd, one of the atom is removed. Default is 1.
#' @param L (optional) This parameter controls the maximum number of non-zero elements
#' in each column of sparsecolding A.Default is m.

#' @return The learned dictionary \code{dictionary}, its size \code{p}
#' and the denoised image \code{fitted image}.
#'
#' @examples
#' I = lena_crop #use a smaller image as an example
#' ## add noise
#' sigma=20;
#' I_noise=AddNoise(I,sigma)
#' ## denoising using GSCAD
#' \dontrun{
#' out=gscad.denoise(I_noise,sigma,m=64,p0=100)
#' }

#' @export
gscad.denoise <- function(I_noise,sigma,D0=NULL,m,p0,c=3.7, lambda=0.05,
                          maxrun=20, maxrun_ADMM=20, err_bnd=1e-4,err_bnd2=1e-4,
                          rho=16, cor_bnd=1,L=NULL)
{
  # m is the size of the patches eg m=8*8, 16*16
  # p0 initial size of the dictionary
  if(!is.null(D0)){
    m=nrow(D0);
    p0=ncol(p0);
  }


  if(is.null(L)){L=m}

  ## get patches
  Y_nc = ImageSplit(I_noise,sqrt(m),sqrt(m));
  mu=colMeans(Y_nc)
  Y=Y_nc-rep(mu,each=nrow(Y_nc))
  n=ncol(Y)

  ## learn the dictionary

  dic=gscad.DL(Y,D0, p0,sigma,c, lambda,
               maxrun, maxrun_ADMM, err_bnd,err_bnd2,
               rho,cor_bnd,L,
               LassoMode=1,LassoLambda=NULL)
  p=ncol(dic)

  ## normalize dictionary using L2 norm
  dic_n=dic/rep(sqrt(colSums(dic^2)),each=nrow(dic))

  ## final estimate of the spares coding
  eps=qchisq(0.9,m)*(sigma/255)^2
  alpha_rcn=.omp(Y,dic_n,eps=eps)

  ## reconstruct image by Y_hat=dic_hat * A_hat
  Y_rcn=as.matrix(dic_n %*% alpha_rcn)
  Y_rcn=Y_rcn+rep(mu,each=nrow(Y))
  Image_fit=ImageReCon(Y_rcn,nrow(I_noise),ncol(I_noise),sqrt(m),sqrt(m))

  list(dictionary=dic,fitted_image=Image_fit, p=p)
}

