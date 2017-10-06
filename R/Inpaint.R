#' Use a given dictionary D to inpaint image
#'
#' The corrupted image is split into sqrt(m) by sqrt(m)
#' patches. Each patch is vectorized into a column of matrix Y.
#' Mask matrix M is constructed to indicate the location of missing pixels.
#' Using the given D, the sparse coding A_hat in Y=DA is obtained.
#' Then Y_inpaint=DA_hat. Then the final
#' inpainted image is reconstruncted.

#' See \url{https://arxiv.org/abs/1605.07870}
#'
#' @param I_corrupt The image to be inpainted. In form of matrix.
#' @param I_mask {0,1} matrix of the same size as I_corrupt
#' to indicate the location of corrupted pixels.
#' @param D D is the dictionary used in Y=DA to inpaint.
#' @param L (optional) This parameter controls the maximum number of non-zero elements
#' in each column of sparsecolding A.
#' @param eps (optional) A lasso tuning paramter
#' @param sigma Noise level.
#' @param stepsize (optional) The stepsize when splicting the image. Default is 1

#' @return The denoised image in for of a matrix.
#' @examples
#' I=lena_crop
#' ## corrupt 30% of the image
#' out_corrupt=AddHoles(I,0.3)
#' I_corrupt=out_corrupt$corruptedImage
#' I_mask=out_corrupt$maskImage
#' ## use ODCT dictionary
#' D0=ODCT(64,100)
#' ## inpaint
#' out=inpaintImage(I_corrupt,I_mask,D0)

#' @export
inpaintImage <- function(I_corrupt,I_mask,D,L=30,eps=NULL,sigma=0,stepsize=1)
{

  m=nrow(D)

  if(sigma==0){
    if(is.null(eps)) eps=qchisq(0.9,m)*(1/255)^2
  }else{
    eps=qchisq(0.9,m)*(sigma/255)^2
  }

  # get patches
  Y_nc = ImageSplit(I_corrupt,sqrt(m),sqrt(m),stepsize);
  M = ImageSplit(I_mask,sqrt(m),sqrt(m),stepsize);

  mu=colSums(Y_nc*M)/colSums(M)
  Y=Y_nc-M*rep(mu,each=nrow(Y_nc))
  mask = matrix(as.logical(M),ncol=ncol(M))

  ## inpaint
  D=D/rep(sqrt(colSums(D^2)),each=nrow(D))
  alpha_rcn = .lassoMask(Y,D,mask,L=L,lambda1=eps,mode=1)
  Y_rcn=as.matrix(D %*% alpha_rcn)
  Y_rcn=Y_rcn+rep(mu,each=nrow(Y))

  ## reconstruct
  Image_fit=ImageReCon(Y_rcn,nrow(I_corrupt),ncol(I_corrupt),sqrt(m),sqrt(m))
  Image_ip=Image_fit*(1-I_mask)+I_corrupt*I_mask
  Image_ip
}


#' Use a given dictionary D to inpaint image
#'
#' Given D, obtain the sparse coding A_hat in Y=DA.
#' Then Y_inpaint=DA_hat.

#' See \url{https://arxiv.org/abs/1605.07870}
#'
#' @param Y Each column of Y is a vectorized image patch to be denoised.
#' @param Mask {0,1} matrix of the same size as Y
#' to indicate the location of corrupted pixels.
#' @param D D is the dictionary used in Y=DA to inpaint.
#' @param L (optional) This parameter controls the maximum number of non-zero elements
#' in each column of sparsecolding A.
#' @param eps (optional) A lasso tuning paramter
#' @param sigma Noise level.

#' @return The inpainted matrix Y.
#'
#' @examples
#' I=lena_crop
#' ## corrupt 30% of the image
#' out_corrupt=AddHoles(I,0.3)
#' I_corrupt=out_corrupt$corruptedImage
#' I_mask=out_corrupt$maskImage
#' ## split image
#' m=64
#' Y_nc = ImageSplit(I_corrupt,sqrt(m),sqrt(m));
#' M = ImageSplit(I_mask,sqrt(m),sqrt(m));
#' mu=colSums(Y_nc*M)/colSums(M)
#' Y=Y_nc-M*rep(mu,each=nrow(Y_nc))
#' mask = matrix(as.logical(M),ncol=ncol(M))
#' ## use ODCT dictionary
#' D0=ODCT(64,100)
#' ## inpaint
#' Y_inpaint=inpaint(Y,mask,D0)

#' @export
inpaint <- function(Y,Mask,D,L=30,eps=NULL,sigma=0)
{

  m=nrow(D)

  if(sigma==0){
    if(is.null(eps)) eps=qchisq(0.9,m)*(1/255)^2
  }else{
    eps=qchisq(0.9,m)*(sigma/255)^2
  }

  ## inpaint
  D=D/rep(sqrt(colSums(D^2)),each=nrow(D))
  alpha_rcn = .lassoMask(Y,D,Mask,L=L,lambda1=eps,mode=1)
  Y_rcn=as.matrix(D %*% alpha_rcn)

  ## reconstruct
  Y_ip=Y_rcn*(1-Mask)+Y*Mask
  Y_ip
}
