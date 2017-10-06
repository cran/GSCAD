#' @useDynLib GSCAD
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix nnzero
#' @import Rcpp
#' @importFrom grDevices gray
#' @importFrom methods slot
#' @importFrom stats cor qchisq rnorm
NULL



######################### GSCAD_DL ##################################
#' Learn dictionary under GSCAD regularization
#'
#' This function learns the dictionary D under the framework of matrix
#' factorization \eqn{Y=DA}. Y is a given m by n matrix of n samples.
#' D is an m by p matrix, where p is unkonw as well, and A ia a p by n matrix.
#' Both D and A are unkonwn. GSCAD regularization is applied to D and lasso
#' regularization is applied to A.
#'
#'
#' See \url{https://arxiv.org/abs/1605.07870}
#'
#' @param Y An m by n matrix in Y=DA. Each column of Y is a sample of size m
#' (usually a vectorized sqrt(m) by sqrt(m) patches).
#' @param D0 Initial dictionary. If D0 specified, p0 is not needed, otherwise D0
#' is evaluated as overcompleted DCT basis using function ODCT(m,p0). Either D0
#' or p0 needs to be specified.
#' @param p0 Initial size of the dictionary.
#' @param sigma Noise level.
#' @param c,lambda Parameters for GSCAD.
#' @param maxrun (optional) Maximun number of outer iterations to run. Default is 20.
#' @param maxrun_ADMM (optional) Maximun number of iterations to run for updating
#' dictionary using ADMM. Default is 20.
#' @param err_bnd (optional) Stopping criterion for iterations. Default is 1e-6.
#' @param err_bnd2 (optional) Stopping criterion for updating dictionary
#' UpDic. Default is 1e-4.
#' @param rho (optional) Parameter for ADMM. Default is 16.
#' @param cor_bnd (optional) When normalize dictionary, checking if the correlation
#' of any two atoms are above the cor_bnd, one of the atom is removed. Default is 1.
#' @param L (optional) This parameter controls the maximum number of non-zero elements
#' in each column of sparsecolding A.
#' @param LassoMode (optional) At the sparse coding stage, the optimization can be
#' done in three modes L1COEFFS (0), L2ERROR (1), PENALTY(2). Default is 1.
#' @param LassoLambda (optional) Tuning parameter for Lasso
#' @param Mask {0,1} matrix of the same size as Y
#' to indicate the location of corrupted pixels.


#' @return The learned dictionary \code{dictionary}

#' @examples
#' I = lena_crop #use a smaller image as an example
#' ## add noise
#' sigma=20;m=64
#' I_noise=AddNoise(I,sigma)
#' ## spliting image into patches
#' Y_nc = ImageSplit(I_noise,sqrt(m),sqrt(m));
#' mu=colMeans(Y_nc)
#' Y=Y_nc-rep(mu,each=nrow(Y_nc))
#' ## learning dictionary
#' \dontrun{
#' dictionary=gscad.DL(Y,p0=256,sigma=sigma)
#' }

#' @export
gscad.DL<- function(Y,D0=NULL, p0,sigma=0, c=3.7, lambda=0.05,
                    maxrun=20, maxrun_ADMM=20, err_bnd=1e-6,err_bnd2=1e-4,
                    rho=16,cor_bnd=1,L=30,
                    LassoMode=1,LassoLambda=NULL)
{
  # m is the size of the patches eg m=8*8, 16*16
  # p0 initial size of the dictionary

  ## GSCAD initialization
  m=nrow(Y);n=ncol(Y)

  ## check convexity condition
  if(lambda^2 > rho/m |
     (c-1)*{rho*(1+lambda^2)^2-m*lambda^2} < 1+lambda^2){
    print('warning: initial parameters do not satisfy convexity condition ')
  }


  if(is.null(D0)){
    D0=ODCT(m,p0)
  }
  dic=D0

  if(LassoMode == 1)
  {
    if(sigma==0){
      LassoLambda = qchisq(0.9,m)*(1/255)^2
    }else{
      LassoLambda = qchisq(0.9,m)*(sigma/255)^2
    }
  }
  if(LassoMode == 2 & is.null(LassoLambda)){
    LassoLambda = 1.2/sqrt(m)
  }

  alpha=matrix();err_alpha=err_dic=-1;
  nrun=1;conti=T;
  while(nrun<=maxrun & conti==T){
    ## update x ##
    alpha_old=alpha
    alpha <- .lasso(Y,dic,lambda1 = LassoLambda, L=L ,mode = LassoMode)
    # print(paste('nnzero of alpha = ', nnzero(alpha)/ncol(alpha)))



    ## update dictionary ##
    dic_old=dic
    out_dic<-UpDic(Y,alpha,rho=rho,c=c,lambda=lambda,
                   err_bnd=err_bnd2,maxrun_ADMM=maxrun_ADMM)



    ## normalize columns of B
    dic=NormMax(out_dic,cor_bnd)
    #dic=out_dic
    # if(nrun %% 5 ==1){
    #   image(dic,main=paste('natom=',ncol(dic)),col=gray((0:32)/32))
    # }

    ## print intermediate result for the current iteration
    nnzero=sum(colSums(dic^2)!=0)
    print(paste('Iteration ',nrun, ': dictionary size = ', nnzero, sep=''))

    ## check stopping
    if(ncol(dic_old)==ncol(dic) & nrow(alpha_old)==nrow(alpha)){
      p=ncol(dic)
      err_alpha<-sum((alpha-alpha_old)^2)/(n*m)
      err_dic <- sum((dic-dic_old)^2)/(p*m)
      print(paste('err_alpha=',err_alpha,' err_dic=',err_dic,sep=''))
      if(err_alpha<err_bnd & err_dic<err_bnd) conti=F
    }


    # if dictionary size is very small, restart the procedure
    if(ncol(dic)<=1){
      conti=F;
      print('Bad initial value. Restarting the procedure with a new initial value.
            If this keep happening, increase the initial size of the dictionary p0.')
      #dic=D0;
      #nrun=0;
    }
    nrun=nrun+1
  }

  dic
}
######################### UpDic ##################################
#' Update dictionary in GSCAD
#'
#' When the sparce coding \code{A={alpha_1,...,alpha_n}} is given, update
#'  dictionary by solving problem (6) in \url{https://arxiv.org/abs/1605.07870}
#'  using ADMM. #'
#'
#' @param Y The image to be denoised. Inform of matrix.
#' @param alpha Sparse coding.
#' @param rho Parameter for the augmented Lagrangian function.
#' @param c,lambda parameters for GSCAD
#' @param maxrun_ADMM Maximun number of iterations run for ADMM
#' @param err_bnd Stopping criterion for iterations.
#' @param mask {0,1} matrix of the same size as Y
#' to indicate the location of corrupted pixels.

#' @return The updated sparse dictionary .
#' @export
UpDic=function(Y,alpha,rho=1,lambda=0.01,c=3.7,maxrun_ADMM=100,err_bnd=1e-4){
  ## check convexity

  ## initial
  alpha=as.matrix(alpha)
  p=nrow(alpha)
  n=ncol(Y)
  m=nrow(Y)
  D2=matrix(0,nrow=m,ncol=p)
  xi=D2;

  AAt0=tcrossprod(alpha)
  AAt=AAt0
  diag(AAt)=diag(AAt)+rho
  alpha_yt=tcrossprod(alpha,Y)



  nrun_ADMM=1;continue_ADMM=TRUE;
  while(continue_ADMM & nrun_ADMM<=maxrun_ADMM){
    ## update D1
    D1=t(solve(AAt,rho*t(D2-xi)+alpha_yt))
    D1_colmax=apply(D1,2,function(d) max(abs(d)))
    c0=max(D1_colmax)
    if(c0>1){
      c0=1
    }
    # c0=1
    D1_scale=pmax(D1_colmax,c0)
    D1=D1/rep(D1_scale,each=m)

    ## update D2
    D2_old=D2
    # D2=.Call('GSCAD_D2optim_mat', PACKAGE = 'GSCAD', as.matrix(D1+xi), lambda, c, rho)
    D2=D2optim_mat(as.matrix(D1+xi), lambda, c, rho)

    ## update xi
    xi=xi+(D1-D2)

    ## update rho
    err_r=sum((D1-D2)^2)/(m*p)
    err_s=sum((D2-D2_old)^2)/(m*p)
    #
    if(nrun_ADMM<=10 & err_r>10*err_s ){
      rho=2*rho
      AAt=AAt0
      diag(AAt)=diag(AAt)+rho
    }

    if(nrun_ADMM<=10 & err_s>10*err_r ){
      rho_2=rho/2
      if(lambda^2<= rho_2/m &
         (c-1)*{rho_2*(1+lambda^2)^2-m*lambda^2} >= 1+lambda^2){
        rho=rho/2
        AAt=AAt0
        diag(AAt)=diag(AAt)+rho
      }
    }

    if(err_r<err_bnd & err_s<err_bnd) continue_ADMM=FALSE
    nrun_ADMM=nrun_ADMM+1
  }#while
  print(paste('nrADMM', nrun_ADMM,'primal error:' ,err_r, ', dual error:', err_s))

  D2
}

########################## gscad.DLmask#########################################
#' @describeIn gscad.DL Adding Mask {0,1} matrix of the same size as Y to indicate
#' the location of corrupted pixel

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
#' ## learn dictionary for inpainting, this function is slow
#' \dontrun{
#' dic=gscad.DLmask(Y, mask,  p0=100, sigma=1)
#' }
#' @export
gscad.DLmask<- function(Y, Mask, D0=NULL, p0, sigma=0, c=3.7, lambda=0.05,
                        maxrun=20, maxrun_ADMM=20, err_bnd=1e-6, err_bnd2=1e-6,rho=16,
                        LassoMode=1,  LassoLambda=NULL)
{
  # m is the size of the patches eg m=8*8, 16*16
  # p0 initial size of the dictionary
  m=nrow(Y);n=ncol(Y);
  ## GSCAD initialization
  if(is.null(D0)){
    D0=ODCT(m,p0)
  }
  dic=D0
  mask= matrix(as.logical(Mask),ncol=ncol(Mask))

  if(LassoMode == 1)
  {
    if(sigma==0){
      LassoLambda = qchisq(0.9,m)*(1/255)^2
    }else{
      LassoLambda = qchisq(0.9,m)*(sigma/255)^2
    }
  }

  alpha=matrix();
  nrun=1;conti=T;
  while(nrun<=maxrun & conti==T){
    ## update x ##
    alpha_old=alpha
    # alpha = spams.ompMask(Y,dic,mask,L=L)
    alpha <- .lassoMask(Y,dic, mask,
                        lambda1 = LassoLambda,mode = LassoMode)

    ## update dictionary ##
    dic_old=dic
    out_dic<-UpDicMask(Y,alpha,mask,
                       rho=rho,c=c,lambda=lambda,err_bnd=err_bnd,maxrun_ADMM=maxrun_ADMM)
    # normalize columns of B
    dic=NormMax(out_dic,0.95)

    #PlotDic(dic,title=ncol(dic))
    ## print intermediate result for the current iteration
    print(paste('Iteration ',nrun, ': dictionary size = ', ncol(dic), sep=''))

    ## check stopping
    if(ncol(dic_old)==ncol(dic) & nrow(alpha_old)==nrow(alpha)){
      p=ncol(dic)
      err_alpha<-sum((alpha-alpha_old)^2)/n
      err_dic <- sum((dic-dic_old)^2)/p
      print(paste('err_alpha=',err_alpha,', err_dic=',err_dic,sep=''))
      if(err_alpha<err_bnd & err_dic<err_bnd) conti=F
    }


    # if dictionary size is very small, restart the procedure
    if(ncol(dic)<=1){
      conti=F;
      print('Bad initial value. Restarting the procedure with a new initial value.
            If this keep happening, increase the initial size of the dictionary p0.')
      #dic=D0;
      #nrun=0;
    }
    nrun=nrun+1
  }

  dic
}

######################### UpDicMask ##################################
#' @describeIn UpDic Adding mask {0,1} matrix of the same size as Y to indicate
#' the location of corrupted pixel


#' @export
UpDicMask=function(Y,alpha,mask, rho=1,lambda=0.01,c=3.7,maxrun_ADMM=100,err_bnd=1e-4){

  ## initial
  p=nrow(alpha)
  n=ncol(Y)
  m=nrow(Y)
  D1=matrix(0,nrow=m,ncol=p)
  D2=D1;xi=D1;
  alpha=as.matrix(alpha)

  nrun_ADMM=1;continue_ADMM=TRUE;
  while(continue_ADMM & nrun_ADMM<=maxrun_ADMM){
    ## update D1
    # print(paste('ADMM iter=', nrun_ADMM, ', Update D1'))
    for(j in 1:m){
      mask1_loc=(1:n)[mask[j,]]
      A=tcrossprod(alpha[,mask1_loc])
      diag(A)=diag(A)+rho
      # D1[j,]=solve(A,
      #              rho*(D2[j,]-xi[j,])+tcrossprod(alpha[,mask1_loc],Y[j,mask1_loc]))@x
      D1[j,]=solve(A,
                   rho*(D2[j,]-xi[j,])+alpha[,mask1_loc]%*% Y[j,mask1_loc])

    }

    D1_colmax=apply(D1,2,function(d) max(abs(d)))
    c0=max(D1_colmax)
    if(c0>1){
      c0=1
    }
    # c0=1
    D1_scale=pmax(D1_colmax,c0)
    D1=D1/rep(D1_scale,each=m)

    ## update D2
    D2_old=D2
    # D2=.Call('GSCAD_D2optim_mat', PACKAGE = 'GSCAD', as.matrix(D1+xi), lambda, c, rho)
    D2=D2optim_mat(as.matrix(D1+xi), lambda, c, rho)


    ## update xi
    xi=xi+(D1-D2)

    ## update rho
    err_r=sum((D1-D2)^2)/(m*p)
    err_s=sum((D2-D2_old)^2)/(m*p)
    # print(paste('primal error:' ,err_r, ', dual error:', err_s))

    #
    if(nrun_ADMM<=10 & err_r>10*err_s ){
      rho=2*rho
      # print(paste('update rho=', rho))
    }

    if(nrun_ADMM<=10 & err_s>10*err_r ){
      rho_2=rho/2
      if(lambda^2<= rho_2/m &
         (c-1)*{rho_2*(1+lambda^2)^2-m*lambda^2} >= 1+lambda^2){
        rho=rho/2
      }
      # print(paste('update rho=', rho))
    }

    if(err_r<err_bnd & err_s<err_bnd) continue_ADMM=FALSE
    nrun_ADMM=nrun_ADMM+1
  }#while

  D2
}


