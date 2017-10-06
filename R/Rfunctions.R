#' @importFrom fields image.plot
#' @importFrom fields tim.colors
NULL

############### AddNoise ##################################
#' Add noise to a clean image.
#'
#' Add white noise from N(0,sigma^2) to each pixel of clean image I
#' @param I Clean image.
#' @param sigma Noise level.
#' @return A noisy image of the same size of I.
#' @examples
#' sigma=25; I=lena_crop
#' I_noise=AddNoise(I,sigma)

#' @export
AddNoise=function(I,sigma){
  dim1=dim(I)[1];dim2=dim(I)[2];
  I_noise = I+matrix(rnorm(dim1*dim2,sd=sigma/255),nrow=dim1)
  I_noise
}

############### AddHole ##################################
#' Corrupt a clean image.
#'
#' Corrupt a clean image with randomly selected locations. The corrupted
#' pixcels are set to 0.

#' @param I Clean image.
#' @param pc Percetage of pixels to corrupt
#' @return A corrupted image \code{corruptedImage} and \code{maskImage} to
#' record the corrupted locations.
#'
#' @examples
#' I=lena_crop
#' I_noise=AddNoise(I,30)
#' @export
AddHoles=function(I,pc){
  dim1=dim(I)[1];dim2=dim(I)[2];
  mask=matrix(1,nrow=dim1)

  N=dim1*dim2;
  loc_crp=sample(1:N, floor(N*pc))
  mask=matrix(1,nrow=dim1,ncol=dim2)
  mask[loc_crp]=0

  #
  I_crp=I
  I_crp[loc_crp]=0

  list(corruptedImage=I_crp,maskImage=mask)
}

############### ODCT ##################################

#' Overcomplete Discrete Cosine Transform (ODCT) Generating Function.
#'
#' \code{ODCT} returns a matrix which is used in dictionary learning as an initial dictionary matrix.
#'
#' @param m Number of grid in each base or the size of each atom in the dictionary.
#' @param p Number of bases or the size of the dictionary.
#' @return A m by p matrix. Each column is a DCT base.
#' @examples
#' dic=ODCT(64,100)
#' # generating an initial dictionary of size 100 with each atom/image patch of size 64.
#'
#'
#' @export
ODCT=function(m,p){
  temp=matrix(cos(rep(0:(m-1),p-1)*pi*rep(1:(p-1),each=m)/p),nrow=m)
  temp=temp-rep(colMeans(temp),each=m)
  temp=temp/rep(sqrt(colSums(temp^2)),each=m)
  temp1=rep(1,m)/sqrt(m)
  dic=cbind(temp1,temp)
  dic
}

############### NormMax ##################################
#' Normalize dictionary atoms using \eqn{l_{\infty}} norm
#'
#' Remove zero columns of D. Then let \eqn{D=(D_1,...,D_p)}. The normalization is done by
#' \deqn{D_ik=D_ik/max_k D_ik, for i=1,...,p.}
#' For any D_i, and D_j, if Cor(D_i,D_j)> cor_bnd, one of D_i and D_j is eliminated.
#'
#' @param D The dictionary to be normalized.
#' @param cor_bnd An upper bound of the correlation between atoms.
#' @return The nomalized dictionary.
#' @examples
#' D=matrix(rnorm(50),nrow=5)
#' D=NormMax(D)
#'

#' @export
NormMax=function(D,cor_bnd=0.95){
  ## remove zero columns
  D_colmax=apply(D,2,function(d) max(abs(d)))
  l_Dnonzero=which(D_colmax!=0)
  D=D[,l_Dnonzero]

  ## normalize each remaining column
  n=nrow(D)
  #cMax=max(abs(D))
  cMax=1
  D_scale=pmax(D_colmax[l_Dnonzero],cMax)
  D=D/rep(D_scale,each=n)

  ## remove columns that are highly correlated ( threshold is cor_bnd )
  p=ncol(D)
  Dcor=cor(D)
  col_keep=unlist(lapply(2:p,function(i){
    all(abs(Dcor[i,1:(i-1)])<cor_bnd)
  }))

  col_keep=c(TRUE,col_keep) # keep the first column
  D=D[,col_keep]

  D
}
############### ImageSplit ##################################
#' Splict image into small patches
#'
#' Splict an image matrix A of size M by N into overlapped small patches of size m by n.
#' Total of (M-m+1)(N-n+1) patches.
#'
#' @param A An M by N matrix.
#' @param m,n Size of the small patches (m by n).
#' @param stepsize (optional) The stepsize when splicting the image. Default is 1
#' @return A (mn) by (M-m+1)(N-n+1) matrix, each column of which is a
#'  vectorized small patch.
#'
#' @export
ImageSplit <- function(A, m, n, stepsize=1) {
  # .Call('GSCAD_ImageSplit', PACKAGE = 'GSCAD', A, m, n, stepsize)
  imageSplit(A, m, n, stepsize)
}

############### ImageReCon ##################################

#' Combine patches of image into full image
#'
#' Suppose an image matrix A of size M by N is splited into (M-m+1)(N-n+1)
#' overlapped small patches of size m by n. Given a matrix B of size (mn) by (M-m+1)(N-n+1)
#' , where each column of B is a vectorized small patch, this function try to
#' reconstruct the full image A.
#'
#'
#' @param B An  (mn) by (M-m+1)(N-n+1) matrix.
#' @param m,n Size of the small patches (m by n).
#' @param mm,nn Full image of size (mm by nn).
#' @param stepsize (optional) The stepsize when splicting the image. Default is 1
#' @return A mm by nn matrix
#'
#' @export
ImageReCon <- function(B, mm, nn, m, n, stepsize=1) {
  # .Call('GSCAD_ImageReCon', PACKAGE = 'GSCAD', B, mm, nn, m, n, stepsize)
  imageReCon(B, mm, nn, m, n, stepsize)
}



########################PlotDic########################################
#' Plot atoms of a given dictionary
#'
#' Given a learned dictionary D (m by p) matrix, this function plots all p atoms
#' of D corresponding to p columns. Each column is mapped to a m1 by m2 small
#' pactch.
#'
#' @param D Dictionary of size m by p
#' @param n1,n2 (Optional) Number of patches to be displayed in each column and
#' each row. If specified, both n1 and n2 have to be specified.
#' @param m1,m2 (Optional) Size of the small patch to be mapped. Default is
#' m1=m2=sqrt(m). If specified, both n1 and n2 have to be specified.
#' @param title Title of the plot.
#' @param color (Optional) If TRUE, the plot is generated using tim.colors()
#' option. If FALSE, the plot is generated using gray scale. Default is FALSE.
#' @param background (Optional) Backgroup color. Default is black.

#' @return Plot will be generated.

#' @examples
#' D=matrix(runif(64*100),nrow=64)
#' # PlotDic(D)

#' @export
PlotDic <- function(D,n1=NULL,n2=NULL, m1=NULL,m2=NULL,title=NULL,color=F,background='black')
{
  p=ncol(D)
  m=nrow(D)
  if(is.null(m1) | is.null(m2)) {
    m1=sqrt(m)
    m2=m/m1
  }

  ## n1 n2: number of dictionary atoms in each col and row
  if(is.null(n1) | is.null(n2)){
    n1=floor(sqrt(p))
    n2=ceiling(p/n1)
  }


  ## Matrix to be plot
  if(background=='white')
  {
    c0=max(D)
  }else if(background == 'gray'){
    c0=(min(D)+max(D))/2
  }else{
    c0=min(D)
  }

  D_plot=matrix(c0,nrow=n1*(m1+1)+1,ncol=n2*(m2+1)+1)
  for(i in 1:p){
    loc_i=ceiling(i/n2)
    loc_j=i-(loc_i-1)*n2

    D_plot[((loc_i-1)*(m1+1)+2):((loc_i-1)*(m1+1)+1+m1),
           ((loc_j-1)*(m2+1)+2):((loc_j-1)*(m2+1)+1+m2)]=matrix(D[,i],ncol=m2)

  }
  D_plot=t(apply(D_plot,2,rev))

  #dev.new(height=n1*m1, width=n2*m2)
  if(color==F){
    image.plot(D_plot,col=gray((1:32)/32),axes = FALSE,main=title,
               asp=n1*m1/n2/m2,legend.shrink=0.8, legend.width=0.5,
               legend.mar = 3.5)
  }else{
    image.plot(D_plot,col=tim.colors(),axes = FALSE,main=title,
               asp=n1*m1/n2/m2,legend.shrink=0.8, legend.width=0.5,
               legend.mar = 3.5)
  }
}
