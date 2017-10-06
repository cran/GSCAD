
#ifndef myDECOMP_h
#define myDECOMP_h

#include "linalg/linalg.h"
static char low='l';
static char nonUnit='n';
/* **************************
 * Implementation in spams.h
 * **************************/
template <typename T>
SpMatrix<T> *_omp(Matrix<T> *X,Matrix<T> *D,Matrix<T> **path,bool return_reg_path,bool given_L,Vector<int>*L,bool given_eps,Vector<T>*eps,bool given_Lambda,Vector<T>*Lambda,const int numThreads) throw(const char *){
  SpMatrix<T> *alpha = new SpMatrix<T>();
  int n = X->m();
  int nD = D->m();
  int K = D->n();
  if (n != nD)
    throw("omp : incompatible matrix dimensions");
  int sizeL = L->n();
  int sizeE = eps->n();
  int sizeLambda = Lambda->n();
  T *pE = eps->rawX();
  T *pLambda = Lambda->rawX();
  int *pL = L->rawX();
  bool vecL = false;
  bool vecEps = false;
  bool vecLambda = false;
  if (! given_L && ! given_eps && ! given_Lambda)
    throw("omp : You should either provide L, eps or lambda");
  int scalar_L = MIN(n,K);
  if(! given_L)
    pL = &scalar_L;
  else if (sizeL > 1)
    vecL = true;
  if(! given_eps) {
    T scalar_eps = 0.;
    pE = &scalar_eps;
  } else if (sizeE > 1)
    vecEps = true;
  if(! given_Lambda) {
    T scalar_Lambda = 0.;
    pLambda = &scalar_Lambda;
  } else if(sizeLambda > 1)
    vecLambda = true;
  if(return_reg_path) {
    *path = new Matrix<T>(K,scalar_L);
    (*path)->setZeros();
    omp((Matrix<T> &)(*X),(Matrix<T> &)(*D),(SpMatrix<T> &)(*alpha),pL,pE,pLambda,vecL,vecEps,vecLambda,numThreads,*path);
  } else {
    *path = NULL;
    omp((Matrix<T> &)(*X),(Matrix<T> &)(*D),(SpMatrix<T> &)(*alpha),pL,pE,pLambda,vecL,vecEps,vecLambda,numThreads);
  }
  return alpha;
}

template <typename T>
SpMatrix<T> *_ompMask(Matrix<T> *X,Matrix<T> *D,Matrix<bool> *B,Matrix<T> **path,bool return_reg_path,bool given_L,Vector<int>*L,bool given_eps,Vector<T>*eps,bool given_Lambda,Vector<T>*Lambda,const int numThreads) throw(const char *){
  SpMatrix<T> *alpha = new SpMatrix<T>();
  int n = X->m();
  int M = X->n();
  int nD = D->m();
  int K = D->n();
  int nM = B->m();
  int mM = B->n();
  if (n != nD )
    throw("ompMask : incompatible matrix dimensions");
  if (nM != n || mM != M)
    throw("ompMask : Mash has non acceptable dimensions");
  int sizeL = L->n();
  int sizeE = eps->n();
  int sizeLambda = Lambda->n();
  T *pE = eps->rawX();
  T *pLambda = Lambda->rawX();
  int *pL = L->rawX();
  bool vecL = false;
  bool vecEps = false;
  bool vecLambda = false;
  if (! given_L && ! given_eps && ! given_Lambda)
    throw("omp : You should either provide L, eps or lambda");
  int scalar_L = MIN(n,K);
  if(! given_L)
    pL = &scalar_L;
  else if (sizeL > 1)
    vecL = true;
  if(! given_eps) {
    T scalar_eps = 0.;
    pE = &scalar_eps;
  } else if (sizeE > 1)
    vecEps = true;
  if(! given_Lambda) {
    T scalar_Lambda = 0.;
    pLambda = &scalar_Lambda;
  } else if(sizeLambda > 1)
    vecLambda = true;
  if(return_reg_path) {
    *path = new Matrix<T>(K,scalar_L);
    (*path)->setZeros();
    omp_mask((Matrix<T> &)(*X),(Matrix<T> &)(*D),(SpMatrix<T> &)(*alpha),(Matrix<bool> &)(*B),pL,pE,pLambda,vecL,vecEps,vecLambda,numThreads,*path);
  } else {
    *path = NULL;
    omp_mask((Matrix<T> &)(*X),(Matrix<T> &)(*D),(SpMatrix<T> &)(*alpha),(Matrix<bool> &)(*B),pL,pE,pLambda,vecL,vecEps,vecLambda,numThreads);
  }
  return alpha;
}
/* **************************
 * Greedy Forward Selection
 * **************************/

/// Forward Selection (or Orthogonal matching pursuit)
/// Address the problem of:
/// \forall i, \min_{\alpha_i} ||X_i-D\alpha_i||_2^2
///                        s.t. ||\alphai||_0 <= L or
/// \forall i, \min_{\alpha_i} ||\alpha_i||_0
///                        s.t. ||\X_i-D\alpha_i||_2^2 <= epsilon
/// This function is
///   * based on Cholesky decompositions
///   * parallel
///   * optimized for a large number of signals (precompute the Gramm matrix

template <typename T>
void omp(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha,
         const int *L, const T* eps, const T* lambda, const bool vecL = false,
         const bool vecEps = false, const bool Lambda=false, const int numThreads=-1,
         Matrix<T>* path = NULL);

// template <typename T>
// void omp_mask(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha, const Matrix<bool>& mask,
//               const int *L, const T* eps, const T* lambda, const bool vecL = false,
//               const bool vecEps = false, const bool Lambda=false, const int numThreads=-1,
//               Matrix<T>* path = NULL);

/// Auxiliary function of omp
template <typename T>
void coreORMP(Vector<T>& scores, Vector<T>& norm, Vector<T>& tmp,
              Matrix<T>& Un, Matrix<T>& Undn, Matrix<T>& Unds, Matrix<T>& Gs,
              Vector<T>& Rdn, const AbstractMatrix<T>& G, Vector<INTM>& ind,
              Vector<T>& RUn, T& normX, const T* eps, const int* L, const T* lambda,
              T* path = NULL);


/// Auxiliary function of omp
template <typename T>
void coreORMPB(Vector<T>& RtD, const AbstractMatrix<T>& G, Vector<INTM>& ind,
               Vector<T>& coeffs, T& normX, const int L, const T eps, const T lambda = 0);

/* *********************
 * Implementation of OMP
 * *********************/

/// Forward Selection (or Orthogonal matching pursuit)
/// Address the problem of:
/// \forall i, \min_{\alpha_i} ||X_i-D\alpha_i||_2^2
///                        s.t. ||\alphai||_0 <= L or
/// \forall i, \min_{\alpha_i} ||\alpha_i||_0
///                        s.t. ||\X_i-D\alpha_i||_2^2 <= epsilon
/// This function is
///   * efficient (Cholesky-based)
///   * parallel
///   * optimized for a big number of signals (precompute the Gramm matrix

template <typename T>
void omp(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha,
         const int* pL, const T* peps, const T* pLambda,
         const bool vecL, const bool vecEps,
         const bool vecLambda, const int numThreads, Matrix<T>* path) {
    int L;
    if (!vecL) {
        L=*pL;
    } else {
        Vector<int> vL(const_cast<int*>(pL),X.n());
        L=vL.maxval();
    }
    spalpha.clear();
    if (L <= 0) return;
    const INTM M = X.n();
    const INTM K = D.n();
    L = MIN(X.m(),MIN(L,K));
    Matrix<T> vM(L,M);
    Matrix<INTM> rM(L,M);

    ProdMatrix<T> G(D, K < 25000 && M > 10);

    int NUM_THREADS=init_omp(numThreads);

    Vector<T>* scoresT=new Vector<T>[NUM_THREADS];
    Vector<T>* normT=new Vector<T>[NUM_THREADS];
    Vector<T>* tmpT=new Vector<T>[NUM_THREADS];
    Vector<T>* RdnT=new Vector<T>[NUM_THREADS];
    Matrix<T>* UnT=new Matrix<T>[NUM_THREADS];
    Matrix<T>* UndnT=new Matrix<T>[NUM_THREADS];
    Matrix<T>* UndsT=new Matrix<T>[NUM_THREADS];
    Matrix<T>* GsT=new Matrix<T>[NUM_THREADS];
    for (int i = 0; i<NUM_THREADS; ++i) {
        scoresT[i].resize(K);
        normT[i].resize(K);
        tmpT[i].resize(K);
        RdnT[i].resize(K);
        UnT[i].resize(L,L);
        UnT[i].setZeros();
        UndnT[i].resize(K,L);
        UndsT[i].resize(L,L);
        GsT[i].resize(K,L);
    }

    int i;
#pragma omp parallel for private(i)
    for (i = 0; i< M; ++i) {
#ifdef _OPENMP
        int numT=omp_get_thread_num();
#else
        int numT=0;
#endif
        Vector<T> Xi;
        X.refCol(i,Xi);
        T normX = Xi.nrm2sq();

        Vector<INTM> ind;
        rM.refCol(i,ind);
        ind.set(-1);

        Vector<T> RUn;
        vM.refCol(i,RUn);

        Vector<T>& Rdn=RdnT[numT];
        D.multTrans(Xi,Rdn);
        coreORMP(scoresT[numT],normT[numT],tmpT[numT],UnT[numT],UndnT[numT],UndsT[numT],
                 GsT[numT],Rdn,G,ind,RUn, normX, vecEps ? peps+i : peps,
                 vecL ? pL+i : pL, vecLambda ? pLambda+i : pLambda,
                 path && i==0 ? path->rawX() : NULL);
    }

    delete[](scoresT);
    delete[](normT);
    delete[](tmpT);
    delete[](RdnT);
    delete[](UnT);
    delete[](UndnT);
    delete[](UndsT);
    delete[](GsT);

    /// convert the sparse matrix into a proper format
    spalpha.convert(vM,rM,K);
};
// template <typename T>
// void omp_mask(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha, const Matrix<bool>& mask,
//               const int *pL, const T* peps, const T* pLambda, const bool vecL,
//               const bool vecEps, const bool vecLambda, const int numThreads,
//               Matrix<T>* path) {
//     int L;
//     if (!vecL) {
//         L=*pL;
//     } else {
//         Vector<int> vL(const_cast<int*>(pL),X.n());
//         L=vL.maxval();
//     }
//     spalpha.clear();
//     if (L <= 0) return;
//     const int M = X.n();
//     const int K = D.n();
//     L = MIN(X.m(),MIN(L,K));
//     Matrix<T> vM(L,M);
//     Matrix<INTM> rM(L,M);
//
//     ProdMatrix<T> G(D, K < 25000 && M > 10);
//
//     int NUM_THREADS=init_omp(numThreads);
//
//     Vector<T>* scoresT=new Vector<T>[NUM_THREADS];
//     Vector<T>* normT=new Vector<T>[NUM_THREADS];
//     Vector<T>* tmpT=new Vector<T>[NUM_THREADS];
//     Vector<T>* RdnT=new Vector<T>[NUM_THREADS];
//     Matrix<T>* UnT=new Matrix<T>[NUM_THREADS];
//     Matrix<T>* UndnT=new Matrix<T>[NUM_THREADS];
//     Matrix<T>* UndsT=new Matrix<T>[NUM_THREADS];
//     Matrix<T>* GsT=new Matrix<T>[NUM_THREADS];
//     ProdMatrix<T>* GT=new ProdMatrix<T>[NUM_THREADS];
//     Matrix<T>* DmaskT=new Matrix<T>[NUM_THREADS];
//     Vector<T>* XmaskT=new Vector<T>[NUM_THREADS];
//     for (int i = 0; i<NUM_THREADS; ++i) {
//         DmaskT[i].resize(D.m(),D.n());
//         XmaskT[i].resize(X.m());
//         scoresT[i].resize(K);
//         normT[i].resize(K);
//         tmpT[i].resize(K);
//         RdnT[i].resize(K);
//         UnT[i].resize(L,L);
//         UnT[i].setZeros();
//         UndnT[i].resize(K,L);
//         UndsT[i].resize(L,L);
//         GsT[i].resize(K,L);
//     }
//
//     int i;
// #pragma omp parallel for private(i)
//     for (i = 0; i< M; ++i) {
// #ifdef _OPENMP
//         int numT=omp_get_thread_num();
// #else
//         int numT=0;
// #endif
//         Vector<T> Xi;
//         X.refCol(i,Xi);
//
//         Vector<INTM> ind;
//         rM.refCol(i,ind);
//         ind.set(-1);
//
//         Vector<T> RUn;
//         vM.refCol(i,RUn);
//
//         Vector<bool> maski;
//         mask.refCol(i,maski);
//         Vector<T>& Rdn=RdnT[numT];
//         if (maski.allfalse()) continue;
//         if (maski.alltrue()) {
//             D.multTrans(Xi,Rdn);
//             T normX = Xi.nrm2sq();
//             coreORMP(scoresT[numT],normT[numT],tmpT[numT],UnT[numT],UndnT[numT],UndsT[numT],
//                      GsT[numT],Rdn,G,ind,RUn, normX, vecEps ? peps+i : peps,
//                      vecL ? pL+i : pL, vecLambda ? pLambda+i : pLambda,
//                      path && i==0 ? path->rawX() : NULL);
//         } else {
//             D.copyMask(DmaskT[numT],maski);
//             Xi.copyMask(XmaskT[numT],maski);
//             T normX = XmaskT[numT].nrm2sq();
//             DmaskT[numT].multTrans(XmaskT[numT],Rdn);
//             GT[numT].setMatrices(DmaskT[numT],false);
//             GT[numT].addDiag(T(1e-10));
//             T eps_mask= (vecEps ? *(peps+i) : *peps)*XmaskT[numT].n()/Xi.n();
//             coreORMP(scoresT[numT],normT[numT],tmpT[numT],
//                      UnT[numT],UndnT[numT],UndsT[numT],
//                      GsT[numT],Rdn,GT[numT],ind,RUn,
//                      normX, &eps_mask, vecL ? pL+i : pL,
//                      vecLambda ? pLambda+i : pLambda,
//                      path && i==0 ? path->rawX() : NULL);
//
//             DmaskT[numT].setm(D.m());
//             DmaskT[numT].setn(D.n());
//             XmaskT[numT].setn(X.m());
//         }
//     }
//
//     delete[](GT);
//     delete[](XmaskT);
//     delete[](DmaskT);
//     delete[](scoresT);
//     delete[](normT);
//     delete[](tmpT);
//     delete[](RdnT);
//     delete[](UnT);
//     delete[](UndnT);
//     delete[](UndsT);
//     delete[](GsT);
//
//     /// convert the sparse matrix into a proper format
//     spalpha.convert(vM,rM,K);
// };
//
// // /// Auxiliary function of omp
// // template <typename T>
// // void coreORMPB(Vector<T>& RtD, const AbstractMatrix<T>& G, Vector<INTM>& ind,
// //                Vector<T>& coeffs, T& normX, const int L, const T eps, const T lambda) {
// //     const int K = G.n();
// //     Vector<T> scores(K);
// //     Vector<T> norm(K);
// //     Vector<T> tmp(K);
// //     Matrix<T> Un(L,L);
// //     Matrix<T> Undn(K,L);
// //     Matrix<T> Unds(L,L);
// //     Matrix<T> Gs(K,L);
// //     ind.set(-1);
// //     coreORMP(scores,norm,tmp,Un,Undn,Unds,Gs,RtD,G,ind,coeffs,normX,&eps,&L,&lambda);
// // };
//
/// Auxiliary function of omp
template <typename T>
void coreORMP(Vector<T>& scores, Vector<T>& norm, Vector<T>& tmp, Matrix<T>& Un,
              Matrix<T>& Undn, Matrix<T>& Unds, Matrix<T>& Gs, Vector<T>& Rdn,
              const AbstractMatrix<T>& G,
              Vector<INTM>& ind, Vector<T>& RUn,
              T& normX, const T* peps, const int* pL, const T* plambda,
              T* path) {
    const T eps = abs<T>(*peps);
    const int L = MIN(*pL,Gs.n());
    const T lambda=*plambda;
    if ((normX <= eps) || L == 0) return;
    const int K = scores.n();
    scores.copy(Rdn);
    norm.set(T(1.0));
    Un.setZeros();

    // permit unsafe low level access
    T* const prUn = Un.rawX();
    //T* const prUnds = Unds.rawX();
    T* const prUndn = Undn.rawX();
    T* const prGs = Gs.rawX();
    T* const prRUn= RUn.rawX();
    if (path)
        memset(path,0,K*L*sizeof(T));

    int j;
    for (j = 0; j<L; ++j) {
        const int currentInd=scores.fmax();
        if (norm[currentInd] < 1e-8) {
            ind[j]=-1;
            break;
        }
        const T invNorm=T(1.0)/sqrt(norm[currentInd]);
        const T RU=Rdn[currentInd]*invNorm;
        const T delta = RU*RU;
        if (delta < 2*lambda) {
            break;
        }

        RUn[j]=RU;
        normX -= delta;
        ind[j]=currentInd;
        //for (int k = 0; k<j; ++k) prUn[j*L+k]=0.0;
        //prUn[j*L+j]=T(1.0);

        //    for (int k = 0; k<j; ++k) prUnds[k*L+j]=prUndn[k*K+currentInd];
        // MGS algorithm, Update Un
        //      int iter = norm[currentInd] < 0.5 ? 2 : 1;
        //int iter=1;
        //     for (int k = 0; k<iter; ++k) {
        ///       for (int l = 0; l<j; ++l) {
        //         T scal=-cblas_dot<T>(j+1-l,prUn+j*L+l,1,prUnds+l*L+l,1);
        //        T scal = -prUnds[l*L+j];
        //         cblas_axpy<T>(l+1,scal,prUn+l*L,1,prUn+j*L,1);
        //       }
        //    }

        prUn[j*L+j]=-T(1.0);
        cblas_copy<T>(j,prUndn+currentInd,K,prUn+j*L,1);
        cblas_trmv<T>(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,j,prUn,L,prUn+j*L,1);
        cblas_scal<T>(j+1,-invNorm,prUn+j*L,1);

        if (j == L-1 || (normX <= eps)) {
            ++j;
            break;
        }

        if (path) {
            T* last_path=path+(L-1)*K;
            cblas_copy<T>(j+1,prRUn,1,last_path,1);
            cblas_trmv<T>(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
                          j+1,prUn,L,last_path,1);
            for (int k = 0; k<=j; ++k) {
                path[j*K+ind[k]]=last_path[k];
            }
        }

        // update the variables Gs, Undn, Unds, Rdn, norm, scores
        Vector<T> Gsj;
        Gs.refCol(j,Gsj);
        G.copyCol(currentInd,Gsj);
        cblas_gemv<T>(CblasColMajor,CblasNoTrans,K,j+1,T(1.0),prGs,K,prUn+j*L,1,
                      T(0.0),prUndn+j*K,1);
        // prUnds[j*L+j] = prUndn[j*K+currentInd];
        Vector<T> Undnj;
        Undn.refCol(j,Undnj);
        Rdn.add(Undnj,-RUn[j]);
        tmp.sqr(Undnj);
        norm.sub(tmp);
        scores.sqr(Rdn);
        scores.div(norm);
        for (int k = 0; k<=j; ++k) scores[ind[k]]=T();
    }
    // compute the final coefficients
    cblas_trmv<T>(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
                  j,prUn,L,prRUn,1);
    if (path) {
        memset(path+(L-1)*K,0,L*sizeof(T));
        for (int k = 0; k<j; ++k) {
            path[(j-1)*K+ind[k]]=prRUn[k];
        }
    }
};

/////////////////////////// Lasso ////////////////////////////////////////
//// decomp.h ////

enum constraint_type { L1COEFFS, L2ERROR, PENALTY};

/// Implementation of LARS-Lasso for solving
/// \forall i, \min_{\alpha_i} ||X_i-D\alpha_i||_2^2
///                        s.t. ||\alphai||_1 <= constraint or
/// \forall i, \min_{\alpha_i} ||\alpha_i||_1
///                        s.t. ||\X_i-D\alpha_i||_2^2 <= constraint or
/// \forall i, \min_{\alpha_i} constraint*||\alpha_i||_1 + ...
///                        ... ||\X_i-D\alpha_i||_2^2 <= T
/// Optionally, the solution might be positive (boolean pos), and a
/// Least-Square can be solved as a post-processing step.
/// L is a maximum number of coefficients.
/// This function is
///   * efficient (Cholesky-based)
///   * parallel
///   * optimized for a big number of signals (precompute the Gramm matrix
template <typename T>
void lasso(const Matrix<T>& X, const Matrix<T>& D,
           SpMatrix<T>& spalpha,
           int L, const T constraint, const T lambda2 = 0, constraint_type mode = PENALTY,
           const bool pos = false, const bool ols = false, const int numThreads=-1,
           Matrix<T>* path = NULL, const int length_path=-1);

template <typename T>
void lasso(const Data<T>& X, const AbstractMatrix<T>& G, const AbstractMatrix<T>& DtX,
           SpMatrix<T>& spalpha,
           int L, const T constraint, constraint_type mode = PENALTY,
           const bool pos = false, const bool ols = false, const int numThreads=-1,
           Matrix<T>* path = NULL, const int length_path=-1);

/// second implementation using matrix inversion lemma
template <typename T>
void lasso2(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha,
            int L, const T constraint,const T lambda2=0, constraint_type mode = PENALTY, const bool pos = false,
            const int numThreads = -1, Matrix<T>* path = NULL, const int length_path=-1);

template <typename T>
void lasso2(const Data<T>& X, const AbstractMatrix<T>& G, const AbstractMatrix<T>& DtX,
            SpMatrix<T>& spalpha,
            int L, const T constraint, constraint_type mode = PENALTY, const bool pos = false,
            const int numThreads = -1, Matrix<T>* path = NULL, const int length_path=-1);

/// second implementation using matrix inversion lemma
template <typename T>
void lasso_mask(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha, const Matrix<bool>& mask,
                int L, const T constraint,const T lambda2=0, constraint_type mode = PENALTY, const bool pos = false,
                const int numThreads = -1);

/// Auxiliary function for lasso
template <typename T>
void coreLARS(Vector<T>& Rdn, Vector<T>& Xdn, Vector<T>& A,
              Vector<T>& u, Vector<T>& sig,
              Vector<T>& av, Vector<T>& RUn, Matrix<T>& Un,
              Matrix<T>& Unds, Matrix<T>& Gs,
              Matrix<T>& Gsa, Matrix<T>& workT, Matrix<T>& R,
              const AbstractMatrix<T>& G,T& normX,
              Vector<int>& ind,Vector<T>& coeffs,const T constraint,
              const bool ols = false,
              const bool pos =false,
              constraint_type mode = L1COEFFS,
              T* path = NULL, int length_path=-1);

template <typename T>
void coreLARS2(Vector<T>& DtR, const AbstractMatrix<T>& G,
               Matrix<T>& Gs,
               Matrix<T>& Ga,
               Matrix<T>& invGs,
               Vector<T>& u,
               Vector<T>& coeffs,
               Vector<INTM>& ind,
               Matrix<T>& work,
               T& normX,
               const constraint_type mode,
               const T constraint, const bool pos = false,
               T* pr_path = NULL, int length_path = -1);

template <typename T>
void coreLARS2(Vector<T>& DtR, const AbstractMatrix<T>& G,
               Vector<T>& coeffs, T normX,
               const constraint_type mode,
               const T constraint, const bool pos = false);

template <typename T>
void coreLARS2W(Vector<T>& DtR, const AbstractMatrix<T>& G,
                Matrix<T>& Gs,
                Matrix<T>& Ga,
                Matrix<T>& invGs,
                Vector<T>& u,
                Vector<T>& coeffs,
                const Vector<T>& weights,
                Vector<INTM>& ind,
                Matrix<T>& work,
                T& normX,
                const constraint_type mode,
                const T constraint, const bool pos = false);

template <typename T>
void coreLARS2W(Vector<T>& DtR, const AbstractMatrix<T>& G,
                Vector<T>& coeffs, const Vector<T>& weights, T normX,
                const constraint_type mode,
                const T constraint, const bool pos = false);

/// Auxiliary functoni for coreLARS (Cholesky downdate)
template <typename T>
void downDateLasso(int& j,int& minBasis,T& normX,const bool ols,
                   const bool pos, Vector<T>& Rdn, INTM* ind,
                   T* coeffs, Vector<T>& sig, Vector<T>& av,
                   Vector<T>& Xdn, Vector<T>& RUn,Matrix<T>& Unm, Matrix<T>& Gsm,
                   Matrix<T>& Gsam, Matrix<T>& Undsm, Matrix<T>& Rm);


/* **************
* LARS - Lasso
* **************/

/// Implementation of LARS-Lasso for solving
/// \forall i, \min_{\alpha_i} ||X_i-D\alpha_i||_2^2
///                        s.t. ||\alphai||_1 <= constraint or
/// \forall i, \min_{\alpha_i} ||\alpha_i||_1
///                        s.t. ||\X_i-D\alpha_i||_2^2 <= constraint or
/// \forall i, \min_{\alpha_i} constraint*||\alpha_i||_1 + ...
///                        ... ||\X_i-D\alpha_i||_2^2 <= T
/// Optionally, the solution might be positive (boolean pos), and a
/// Least-Square can be solved as a post-processing step.
/// L is a maximum number of coefficients.
/// This function is
///   * efficient (Cholesky-based)
///   * parallel
///   * optimized for a big number of signals (precompute the Gramm matrix

template <typename T>
void lasso(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha,
           int L, const T lambda, const T lambda2, constraint_type mode,
           const bool pos, const bool ols, const int numThreads,
           Matrix<T>* path, const int length_path) {
  ProdMatrix<T> G(D, X.n() > 10 && D.n() < 50000);
  G.addDiag(MAX(lambda2,1e-10));
  ProdMatrix<T> DtX(D,X,false);
  lasso(X,G,DtX,spalpha,L,lambda,mode,pos,ols,numThreads,path,length_path);
}

template <typename T>
void lasso(const Data<T>& X, const AbstractMatrix<T>& G,
           const AbstractMatrix<T>& DtX, SpMatrix<T>& spalpha,
           int L, const T lambda, constraint_type mode,
           const bool pos, const bool ols, const int numThreads,
           Matrix<T>* path, const int length_path) {

  spalpha.clear();
  const INTM M = X.n();
  const INTM K = G.n();
  Matrix<T> vM;
  Matrix<INTM> rM;
  vM.resize(L,M);
  rM.resize(L,M);

  if (L <= 0) return;
  if (path) path->setZeros();

  int NUM_THREADS=init_omp(numThreads);

  //ProdMatrix<T> G(D, K < 25000 && M > 10);

  Vector<T>* RdnT=new Vector<T>[NUM_THREADS];
  Vector<T>* XdnT =new Vector<T>[NUM_THREADS];
  Vector<T>* AT=new Vector<T>[NUM_THREADS];
  Vector<T>* uT=new Vector<T>[NUM_THREADS];
  Vector<T>* sigT=new Vector<T>[NUM_THREADS];
  Vector<T>* avT=new Vector<T>[NUM_THREADS];
  Vector<T>* RUnT = new Vector<T>[NUM_THREADS];
  Matrix<T>* UnT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* RT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* UndsT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* GsT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* GsaT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* workT=new Matrix<T>[NUM_THREADS];
  for (int i = 0; i<NUM_THREADS; ++i) {
    RdnT[i].resize(K);
    if (ols) XdnT[i].resize(K);
    AT[i].resize(K);
    uT[i].resize(L);
    sigT[i].resize(L);
    avT[i].resize(L);
    if (ols) RUnT[i].resize(L);
    UnT[i].resize(L,L);
    UnT[i].setZeros();
    UndsT[i].resize(L,L);
    UndsT[i].setZeros();
    GsT[i].resize(K,L);
    GsaT[i].resize(L,L);
    workT[i].resize(K,2);
    RT[i].resize(L,L);
  }

  Vector<T> norms;
  X.norm_2sq_cols(norms);
  int i;
#pragma omp parallel for private(i)
  for (i = 0; i< M; ++i) {
#ifdef _OPENMP
    int numT=omp_get_thread_num();
#else
    int numT=0;
#endif
    T normX = norms[i];

    Vector<INTM> ind;
    rM.refCol(i,ind);
    Vector<T> coeffs;
    vM.refCol(i,coeffs);
    coeffs.setZeros();

    Vector<T>& Rdn=RdnT[numT];
    DtX.copyCol(i,Rdn);
    coreLARS(Rdn,XdnT[numT], AT[numT], uT[numT], sigT[numT], avT[numT],
             RUnT[numT], UnT[numT], UndsT[numT], GsT[numT], GsaT[numT],
                                                                workT[numT],RT[numT],G,normX, ind,coeffs,lambda,ols,pos,
                                                                mode,path && i==0 ? path->rawX() : NULL, length_path);
  }

  delete[](RdnT);
  delete[](XdnT);
  delete[](AT);
  delete[](uT);
  delete[](sigT);
  delete[](avT);
  delete[](RUnT);
  delete[](UnT);
  delete[](RT);
  delete[](UndsT);
  delete[](GsT);
  delete[](GsaT);
  delete[](workT);

  /// convert the sparse matrix into a proper format
  spalpha.convert(vM,rM,K);
};

/// Auxiliary function for lasso
template <typename T>
void coreLARS(Vector<T>& Rdnv, Vector<T>& Xdnv, Vector<T>& Av,
              Vector<T>& uv, Vector<T>& sigv, Vector<T>& avv, Vector<T>& RUnv,
              Matrix<T>& Unm, Matrix<T>& Undsm, Matrix<T>& Gsm,
              Matrix<T>& Gsam, Matrix<T>& workm, Matrix<T>& Rm,
              const AbstractMatrix<T>& Gm,T& normX,
              Vector<INTM>& indv,Vector<T>& coeffsv,const T constraint,
              const bool ols,const bool pos, constraint_type mode,
              T* path, int length_path) {
  if (mode == L2ERROR && normX < constraint) return;

  const int LL = Gsm.n();
  const int K = Gsm.m();
  const int L = MIN(LL,K);
  if (length_path <= 1) length_path=4*L;
  // permit unsafe fast low level access
  T* const Rdn = Rdnv.rawX();
  T* const Xdn = Xdnv.rawX();
  T* const A = Av.rawX();
  T* const u = uv.rawX();
  T* const sig = sigv.rawX();
  //T* const av = avv.rawX();
  T* const RUn = RUnv.rawX();
  T* const Un = Unm.rawX();
  T* const Unds = Undsm.rawX();
  T* const Gs = Gsm.rawX();
  T* const Gsa = Gsam.rawX();
  T* const work = workm.rawX();
  //T* const G = Gm.rawX();
  //T* const R = Rm.rawX();
  INTM* ind = indv.rawX();
  T* coeffs = coeffsv.rawX();

  coeffsv.setZeros();
  indv.set(-1);

  if (ols) Xdnv.copy(Rdnv);
  int currentInd= pos ? Rdnv.max() : Rdnv.fmax();
  bool newAtom=true;
  T Cmax = 0;
  int iter=1;
  T thrs = 0.0;

  //   INTM* const ind_orig = ind;
  //   T* const coeffs_orig = coeffs;

  int j;
  for (j = 0; j<L; ++j) {
    if (newAtom) {
      ind[j]=currentInd;

      if (pos) {
        Cmax = Rdn[currentInd];
        sig[j]=1.0;
      } else {
        Cmax = abs<T>(Rdn[currentInd]);
        sig[j] = SIGN(Rdn[currentInd]);
      }
      for (int k = 0; k<=j; ++k) Un[j*L+k]=0.0;
      Un[j*L+j]=1.0;
      Gm.extract_rawCol(currentInd,Gs+K*j);
      for (int k = 0; k<j; ++k) Gs[K*j+ind[k]] *= sig[k];
      if (sig[j] < 0) {
        Rdn[currentInd]=-Rdn[currentInd];
        if (ols) Xdn[currentInd]=-Xdn[currentInd];
        cblas_scal<T>(K,sig[j],Gs+K*j,1);
        cblas_scal<T>(j+1,sig[j],Gs+currentInd,K);
      }
      cblas_copy<T>(j+1,Gs+currentInd,K,Gsa+j*L,1);
      for (int k = 0; k<j; ++k) Gsa[k*L+j]=Gsa[j*L+k];

      // <d_j,d_i>
      cblas_copy<T>(j,Gsa+j*L,1,Unds+j,L);
      // <U_j final,d_i>
      cblas_trmv<T>(CblasColMajor,CblasUpper,CblasTrans,CblasNonUnit,
                    j+1,Un,L,Unds+j,L);
      // norm2
      T norm2=Gsa[j*L+j];
      for (int k = 0; k<j; ++k) norm2 -= Unds[k*L+j]*Unds[k*L+j];
      if (norm2 < 1e-15) {
        ind[j]=-1;
        //      cerr << "bad exit" << endl;
        break;
      }

      //   int iter2 = norm2 < 0.5 ? 2 : 1;
      //   for(int k = 0; k<iter2; ++k) {
      //      for (int l = 0; l<j; ++l) {
      //         T scal=-cblas_dot<T>(j+1-l,Un+j*L+l,1,Unds+l*L+l,1);
      //         cblas_axpy<T>(l+1,scal,Un+l*L,1,Un+j*L,1);
      //      }
      //   }
      Un[j*L+j]=-T(1.0);
      cblas_copy<T>(j,Unds+j,L,Un+j*L,1);
      cblas_trmv<T>(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,j,Un,L,Un+j*L,1);

      /// Un is the orthogonalized vectors in the D basis
      T invNorm=1.0/sqrt(norm2);
      cblas_scal<T>(j+1,-invNorm,Un+j*L,1);
      Unds[j*L+j]=cblas_dot<T>(j+1,Un+j*L,1,Gsa+j*L,1);
    }

    for (int k = 0; k<=j; ++k) u[k]=T(1.0);
    cblas_trmv<T>(CblasColMajor,CblasUpper,CblasTrans,CblasNonUnit,
                  j+1,Un,L,u,1);

    T a = T(1.0)/cblas_nrm2<T>(j+1,u,1);

    cblas_trmv<T>(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
                  j+1,Un,L,u,1);
    cblas_scal<T>(j+1,a,u,1);

    cblas_gemv<T>(CblasColMajor,CblasNoTrans,K,j+1,T(1.0),Gs,K,u,1,T(0.0),A,1);

    T potentNorm=0.0;
    if (!ols) {
      for (int k = 0; k<=j; ++k)  potentNorm += Rdn[ind[k]]*u[k];
    }

    if (pos) {
      for (int k = 0; k<K; ++k) {
        T diff = a-A[k];
        work[k]= diff <= 0 ? INFINITY : (Cmax-Rdn[k])/diff;
      }
      for (int k = 0; k<=j; ++k) {
        work[ind[k]]=INFINITY;
      }
      for (int k = 0; k<K; ++k)
        if (work[k] <=0) work[k]=INFINITY;
        currentInd =cblas_iamin<T>(K,work,1);
    } else {
      memset(work,0,2*K*sizeof(T));
      for (int k = 0; k<=j; ++k) {
        const int index=2*ind[k];
        work[index]=INFINITY;
        work[index+1]=INFINITY;
      }
      for (int k = 0; k<K; ++k) {
        const int index=2*k;
        if (!work[index]) {
          const T diff1=a-A[k];
          work[index]= diff1 <= 0 ? INFINITY : (Cmax-Rdn[k])/diff1;
          const T diff2=a+A[k];
          work[index+1]=diff2 <= 0 ? INFINITY : (Cmax+Rdn[k])/diff2;
        }
      }
      currentInd =cblas_iamin<T>(2*K,work,1);
    }
    T gamma=work[currentInd];
    T gammaMin=0;
    int minBasis=0;

    //if (j == L-1) gamma=potentNorm;

    if (mode == PENALTY) {
      gamma=MIN(gamma,(Cmax-constraint)/a);
    }

    //      if (j > 0) {
    vDiv<T>(j+1,coeffs,u,work);
    cblas_scal<T>(j+1,-T(1.0),work,1);
    /// voir pour petites valeurs
    for (int k=0; k<=j; ++k)
      if (coeffs[k]==0 || work[k] <=0) work[k]=INFINITY;
      minBasis=cblas_iamin<T>(j+1,work,1);
      gammaMin=work[minBasis];
      if (gammaMin < gamma) gamma=gammaMin;
      //     }

      if (mode == L1COEFFS) {
        T Tu = 0.0;
        for (int k = 0; k<=j; ++k) Tu += u[k];

        if (Tu > EPSILON)
          gamma= MIN(gamma,(constraint-thrs)/Tu);
        thrs+=gamma*Tu;
      }

      // compute the norm of the residdual

      if (ols == 0) {
        const T t = gamma*gamma - 2*gamma*potentNorm;
        if (t > 0 || isnan(t) || isinf(t)) {
          //      cerr << "bad bad exit" << endl;
          //       cerr << t << endl;
          ind[j]=-1;
          break;
        }
        normX += t;
      } else {
        // plan the last orthogonal projection
        if (newAtom) {
          RUn[j]=0.0;
          for (int k = 0; k<=j; ++k) RUn[j] += Xdn[ind[k]]*
            Un[j*L+k];
          normX -= RUn[j]*RUn[j];
        }
      }

      // Update the coefficients
      cblas_axpy<T>(j+1,gamma,u,1,coeffs,1);

      if (pos) {
        for (int k = 0; k<j+1; ++k)
          if (coeffs[k] < 0) coeffs[k]=0;
      }

      cblas_axpy<T>(K,-gamma,A,1,Rdn,1);
      if (!pos) currentInd/= 2;
      if (path) {
        for (int k = 0; k<=j; ++k)
          path[iter*K+ind[k]]=coeffs[k]*sig[k];
      }

      if (gamma == gammaMin) {
        downDateLasso<T>(j,minBasis,normX,ols,pos,Rdnv,ind,coeffs,sigv,
                         avv,Xdnv, RUnv, Unm, Gsm, Gsam,Undsm,Rm);
        newAtom=false;
        Cmax=abs<T>(Rdn[ind[0]]);
        --j;
      } else {
        newAtom=true;
      }
      ++iter;

      if (mode == PENALTY) {
        thrs=abs<T>(Rdn[ind[0]]);
      }

      if ((j == L-1) ||
          (mode == PENALTY && (thrs - constraint < 1e-15)) ||
          (mode == L1COEFFS && (thrs - constraint > -1e-15)) ||
          (newAtom && mode == L2ERROR && (normX - constraint < 1e-15)) ||
          (normX < 1e-15) ||
          (iter >= length_path)) {
        //       cerr << "exit" << endl;
        //       PRINT_F(thrs)
        //       PRINT_F(constraint)
        //       PRINT_F(normX)
        break;
      }

  }
  if (ols) {
    cblas_copy<T>(j+1,RUn,1,coeffs,1);
    cblas_trmv<T>(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
                  j+1,Un,L,coeffs,1);
  }
  vMul<T>(j+1,coeffs,sig,coeffs);
};

/// Auxiliary functoni for coreLARS (Cholesky downdate)
template <typename T>
inline void downDateLasso(int& j,int& minBasis,T& normX,const bool ols,
                          const bool pos,
                          Vector<T>& Rdnv, INTM* ind,
                          T* coeffs, Vector<T>& sigv, Vector<T>& avv,
                          Vector<T>& Xdnv, Vector<T>& RUnv,Matrix<T>& Unm, Matrix<T>& Gsm,
                          Matrix<T>& Gsam, Matrix<T>& Undsm, Matrix<T>& Rm) {
  const int L = Gsm.n();
  const int K = Gsm.m();
  T* const Rdn = Rdnv.rawX();
  T* const Xdn = Xdnv.rawX();
  T* const sig = sigv.rawX();
  T* const av = avv.rawX();
  T* const RUn = RUnv.rawX();
  T* const Un = Unm.rawX();
  T* const Unds = Undsm.rawX();
  T* const Gs = Gsm.rawX();
  T* const Gsa = Gsam.rawX();
  T* const R = Rm.rawX();

  int indB=ind[minBasis];

  if (!pos && sig[minBasis] < 0) {
    // Update Rdn
    Rdn[indB]=-Rdn[indB];
    if (ols) Xdn[indB]=-Xdn[indB];
  }

  int num=j-minBasis;
  for (int k = 0; k<num*num;++k) R[k]=0.0;
  for (int k = 0; k<num; ++k) R[k*num+k]=1.0;
  // Update Un
  for (int k = minBasis+1; k<=j; ++k) {
    T a = -Un[k*L+minBasis]/Un[minBasis*L+minBasis];
    av[k-minBasis-1] = a;
    cblas_axpy<T>(minBasis,a,Un+minBasis*L,1,Un+k*L,1);
  }
  for (int k = minBasis+1; k<=j; ++k) {
    cblas_copy<T>(minBasis,Un+k*L,1,Un+(k-1)*L,1);
    cblas_copy<T>(num,Un+k*L+minBasis+1,1,Un+(k-1)*L+minBasis,1);
  }
  T alpha=1.0;
  T alphab,gamma;
  for (int k = 0; k<num; ++k) {
    alphab=alpha+av[k]*av[k];
    R[k*num+k]=sqrt(alphab/alpha);
    gamma=av[k]*R[k*num+k]/alphab;
    alpha=alphab;
    cblas_copy<T>(num-k-1,av+k+1,1,R+k*num+k+1,1);
    cblas_scal<T>(num-k-1,gamma,R+k*num+k+1,1);
  }
  if (num > 0) {
    trtri<T>(low,nonUnit,num,R,num);
    cblas_trmm<T>(CblasColMajor,CblasRight,CblasLower,CblasTrans,CblasNonUnit,
                  j,num,T(1.0),R,num,Un+minBasis*L,L);
  }

  // Update Unds
  for (int k = minBasis+1; k<=j; ++k)
    cblas_axpy<T>(j-minBasis,av[k-minBasis-1],Unds+minBasis*L+minBasis+1,1,
                  Unds+k*L+minBasis+1,1);
  for (int k = 0; k<minBasis; ++k)
    for (int l = minBasis+1; l<=j; ++l)
      Unds[k*L+l-1]=Unds[k*L+l];
  for (int k = minBasis+1; k<=j; ++k)
    cblas_copy<T>(j-minBasis,Unds+k*L+minBasis+1,1,Unds+(k-1)*L+minBasis,1);
  if (num > 0)
    cblas_trmm<T>(CblasColMajor,CblasRight,CblasLower,CblasTrans,CblasNonUnit,
                  j-minBasis,num,T(1.0),R,num,Unds+minBasis*L+minBasis,L);
  for (int k = minBasis+1; k<=j; ++k)
    for (int l = 0; l<k; ++l) Unds[k*L+l]=0.0;

  // Update Gs
  for (int k = minBasis+1; k<=j; ++k) {
    cblas_copy<T>(K,Gs+k*K,1,Gs+(k-1)*K,1);
  }
  if (!pos && sig[minBasis] < T(0.0)) cblas_scal<T>(j,T(-1.0),Gs+indB,K);
  // Update Gsa
  for (int k = minBasis+1; k<=j; ++k) {
    cblas_copy<T>(minBasis,Gsa+k*L,1,Gsa+(k-1)*L,1);
    cblas_copy<T>(j-minBasis,Gsa+k*L+minBasis+1,1,Gsa+(k-1)*L+minBasis,1);
  }
  for (int k = 0; k<minBasis; ++k) {
    for (int l = minBasis+1; l<=j; ++l) Gsa[k*L+l-1]=Gsa[k*L+l];
  }

  // Update sig
  for (int k = minBasis+1; k<=j && !pos; ++k) sig[k-1]=sig[k];
  // Update ind
  for (int k = minBasis+1; k<=j; ++k) ind[k-1]=ind[k];
  ind[j]=-1;

  for (int k = minBasis+1; k<=j; ++k) coeffs[k-1]=coeffs[k];
  coeffs[j]=0.0;

  if (ols) {
    // Update RUn and normX
    for (int k = minBasis; k<=j; ++k)
      normX += RUn[k]*RUn[k];
    for (int k = minBasis; k<j; ++k) {
      RUn[k]=0.0;
      for (int l = 0; l<=k; ++l) RUn[k] += Xdn[ind[l]]*
        Un[k*L+l];
      normX -= RUn[k]*RUn[k];
    }
  }

  // Update j
  --j;
}





template <typename T>
void lassoWeight(const Matrix<T>& X, const Matrix<T>& D, const Matrix<T>& weights,
                 SpMatrix<T>& spalpha,
                 int L, const T constraint, constraint_type mode, const bool pos,
                 const int numThreads) {

  spalpha.clear();
  const int M = X.n();
  const int K = D.n();
  Matrix<T> vM;
  Matrix<INTM> rM;
  vM.resize(L,M);
  rM.resize(L,M);

  if (L <= 0) return;

  int NUM_THREADS=init_omp(numThreads);

  //ProdMatrix<T> G(D, K < 25000 && M > 10);
  ProdMatrix<T> G(D, K < 50000);
  //Matrix<T> G;
  //D.XtX(G);
  G.addDiag(1e-10);

  Vector<T>* DtRT=new Vector<T>[NUM_THREADS];
  Vector<T>* uT=new Vector<T>[NUM_THREADS];
  Matrix<T>* GsT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* GaT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* invGsT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* workT=new Matrix<T>[NUM_THREADS];
  for (int i = 0; i<NUM_THREADS; ++i) {
    DtRT[i].resize(K);
    uT[i].resize(K);
    uT[i].setZeros();
    GsT[i].resize(L,L);
    invGsT[i].resize(L,L);
    GaT[i].resize(K,L);
    workT[i].resize(K,3);
    workT[i].setZeros();
  }

  int i;
#pragma omp parallel for private(i)
  for (i = 0; i< M; ++i) {
#ifdef _OPENMP
    int numT=omp_get_thread_num();
#else
    int numT=0;
#endif
    Vector<T> Xi;
    X.refCol(i,Xi);
    T normX = Xi.nrm2sq();

    Vector<INTM> ind;
    rM.refCol(i,ind);
    Vector<T> coeffs;
    vM.refCol(i,coeffs);

    Vector<T>& DtR=DtRT[numT];
    D.multTrans(Xi,DtR);
    Vector<T> we;
    weights.refCol(i,we);

    coreLARS2W(DtR,G,GsT[numT],GaT[numT],invGsT[numT],uT[numT],coeffs,we,
               ind,workT[numT],normX,mode,constraint,pos);
  }

  delete[](DtRT);
  delete[](uT);
  delete[](GsT);
  delete[](GaT);
  delete[](invGsT);
  delete[](workT);

  /// convert the sparse matrix into a proper format
  spalpha.convert(vM,rM,K);
};

/// second implementation using matrix inversion lemma
template <typename T>
void lasso_mask(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha, const Matrix<bool>& mask,
                int L, const T constraint,const T lambda2, constraint_type mode, const bool pos,
                const int numThreads) {
  spalpha.clear();
  const int M = X.n();
  const int K = D.n();
  Matrix<T> vM;
  Matrix<INTM> rM;
  vM.resize(L,M);
  rM.resize(L,M);

  if (L <= 0) return;

  int NUM_THREADS=init_omp(numThreads);

  ProdMatrix<T> G(D,K < 25000 && M > 10);
  G.addDiag(MAX(lambda2,1e-10));

  Vector<T>* DtRT=new Vector<T>[NUM_THREADS];
  Vector<T>* uT=new Vector<T>[NUM_THREADS];
  Vector<T>* XmaskT=new Vector<T>[NUM_THREADS];
  Matrix<T>* GsT=new Matrix<T>[NUM_THREADS];
  ProdMatrix<T>* GT=new ProdMatrix<T>[NUM_THREADS];
  Matrix<T>* DmaskT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* GaT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* invGsT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* workT=new Matrix<T>[NUM_THREADS];
  for (int i = 0; i<NUM_THREADS; ++i) {
    DmaskT[i].resize(D.m(),D.n());
    DtRT[i].resize(K);
    uT[i].resize(K);
    XmaskT[i].resize(X.m());
    uT[i].setZeros();
    GsT[i].resize(L,L);
    invGsT[i].resize(L,L);
    GaT[i].resize(K,L);
    workT[i].resize(K,3);
    workT[i].setZeros();
  }

  int i;
#pragma omp parallel for private(i)
  for (i = 0; i< M; ++i) {
#ifdef _OPENMP
    int numT=omp_get_thread_num();
#else
    int numT=0;
#endif
    Vector<T> Xi;
    X.refCol(i,Xi);
    Vector<bool> maski;
    mask.refCol(i,maski);
    Vector<INTM> ind;
    rM.refCol(i,ind);
    Vector<T> coeffs;
    vM.refCol(i,coeffs);
    Vector<T>& DtR=DtRT[numT];

    if (maski.allfalse()) continue;
    if (maski.alltrue()) {
      T normX = Xi.nrm2sq();
      D.multTrans(Xi,DtR);
      coreLARS2(DtR,G,GsT[numT],GaT[numT],invGsT[numT],uT[numT],coeffs,
                ind,workT[numT],normX,mode,constraint,pos);
    } else {
      D.copyMask(DmaskT[numT],maski);
      Xi.copyMask(XmaskT[numT],maski);
      T constraint_mask = mode == PENALTY || mode == L2ERROR ? constraint*XmaskT[numT].n()/Xi.n() : constraint;
      T normX = XmaskT[numT].nrm2sq();
      DmaskT[numT].multTrans(XmaskT[numT],DtR);
      GT[numT].setMatrices(DmaskT[numT],false);
      GT[numT].addDiag(MAX(lambda2,T(1e-10)));
      coreLARS2(DtR,GT[numT],
                GsT[numT],GaT[numT],invGsT[numT],uT[numT],coeffs,
                ind,workT[numT],normX,mode,constraint_mask,pos);
      DmaskT[numT].setm(D.m());
      DmaskT[numT].setn(D.n());
      XmaskT[numT].setn(X.m());
    }
  }

  delete[](GT);
  delete[](XmaskT);
  delete[](DmaskT);
  delete[](DtRT);
  delete[](uT);
  delete[](GsT);
  delete[](GaT);
  delete[](invGsT);
  delete[](workT);

  /// convert the sparse matrix into a proper format
  spalpha.convert(vM,rM,K);

};

template <typename T>
void lasso2(const Matrix<T>& X, const Matrix<T>& D, SpMatrix<T>& spalpha,
            int L, const T constraint, const T lambda2, constraint_type mode, const bool pos,
            const int numThreads, Matrix<T>* path, int length_path) {
  ProdMatrix<T> G(D,X.n() > 10 && D.n() < 50000);
  ProdMatrix<T> DtX(D,X,false);
  G.addDiag(MAX(lambda2,1e-10));
  lasso2(X,G,DtX,spalpha,L,constraint,mode,pos,numThreads,path, length_path);
}


template <typename T>
void lasso2(const Data<T>& X, const AbstractMatrix<T>& G, const AbstractMatrix<T>& DtX,
            SpMatrix<T>& spalpha,
            int L, const T constraint, constraint_type mode, const bool pos,
            const int numThreads, Matrix<T>* path, int length_path) {
  spalpha.clear();
  const INTM M = X.n();
  const INTM K = G.n();
  Matrix<T> vM;
  Matrix<INTM> rM;
  vM.resize(L,M);
  rM.resize(L,M);

  if (L <= 0) return;
  if (path) path->setZeros();

  int NUM_THREADS=init_omp(numThreads);

  Vector<T>* DtRT=new Vector<T>[NUM_THREADS];
  Vector<T>* uT=new Vector<T>[NUM_THREADS];
  Matrix<T>* GsT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* GaT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* invGsT=new Matrix<T>[NUM_THREADS];
  Matrix<T>* workT=new Matrix<T>[NUM_THREADS];
  for (int i = 0; i<NUM_THREADS; ++i) {
    DtRT[i].resize(K);
    uT[i].resize(K);
    uT[i].setZeros();
    GsT[i].resize(L,L);
    invGsT[i].resize(L,L);
    GaT[i].resize(K,L);
    workT[i].resize(K,3);
    workT[i].setZeros();
  }
  INTM i;
  Vector<T> norms;
  X.norm_2sq_cols(norms);
#pragma omp parallel for private(i)
  for (i = 0; i< M; ++i) {
#ifdef _OPENMP
    int numT=omp_get_thread_num();
#else
    int numT=0;
#endif
    //  Vector<T> Xi;
    //  X.refCol(i,Xi);
    //  T normX = Xi.nrm2sq();
    T normX = norms[i];

    Vector<INTM> ind;
    rM.refCol(i,ind);
    Vector<T> coeffs;
    vM.refCol(i,coeffs);

    Vector<T>& DtR=DtRT[numT];
    DtX.copyCol(i,DtR);
    //D.multTrans(Xi,DtR);
    coreLARS2(DtR,G,GsT[numT],GaT[numT],invGsT[numT],
              uT[numT],coeffs,
              ind,workT[numT],normX,mode,constraint,pos,
              path && i==0 ? path->rawX() : NULL,length_path);
  }

  delete[](DtRT);
  delete[](uT);
  delete[](GsT);
  delete[](GaT);
  delete[](invGsT);
  delete[](workT);

  /// convert the sparse matrix into a proper format
  spalpha.convert(vM,rM,K);
};

template <typename T>
void coreLARS2W(Vector<T>& DtR, const AbstractMatrix<T>& G,
                Vector<T>& coeffs, const Vector<T>& weights, T normX,
                const constraint_type mode,
                const T constraint, const bool pos) {
  const INTM p = G.m();
  const INTM L = p;
  Vector<T> v;
  v.resize(L);
  Vector<INTM> r;
  r.resize(L);
  Vector<T> u;
  u.resize(p);
  Matrix<T> Gs;
  Gs.resize(L,L);
  Matrix<T> invGs;
  invGs.resize(L,L);
  Matrix<T> Ga;
  Ga.resize(p,L);
  Matrix<T> work;
  work.resize(p,3);
  coreLARS2W(DtR,G,Gs,Ga,invGs,u,v,weights,r,work,normX,mode,constraint,pos);
  coeffs.setZeros();
  for (int i = 0; i< L; ++i) {
    if (r[i] < 0) break;
    coeffs[r[i]]=v[i];
  };
};

template <typename T>
void coreLARS2(Vector<T>& DtR, const AbstractMatrix<T>& G,
               Vector<T>& coeffs, T normX,
               const constraint_type mode,
               const T constraint, const bool pos) {
  const INTM p = G.m();
  const INTM L = p;
  Vector<T> v;
  v.resize(L);
  Vector<INTM> r;
  r.resize(L);
  Vector<T> u;
  u.resize(p);
  Matrix<T> Gs;
  Gs.resize(L,L);
  Matrix<T> invGs;
  invGs.resize(L,L);
  Matrix<T> Ga;
  Ga.resize(p,L);
  Matrix<T> work;
  work.resize(p,3);
  coreLARS2(DtR,G,Gs,Ga,invGs,u,v,r,work,normX,mode,constraint,pos);
  coeffs.setZeros();
  for (int i = 0; i< L; ++i) {
    if (r[i] < 0) break;
    coeffs[r[i]]=v[i];
  };
};

/// Auxiliary function for lasso
template <typename T>
void coreLARS2(Vector<T>& DtR, const AbstractMatrix<T>& G,
               Matrix<T>& Gs,
               Matrix<T>& Ga,
               Matrix<T>& invGs,
               Vector<T>& u,
               Vector<T>& coeffs,
               Vector<INTM>& ind,
               Matrix<T>& work,
               T& normX,
               const constraint_type mode,
               const T constraint,
               const bool pos,
               T* path, int length_path) {
  const int LL = Gs.n();
  const int K = G.n();
  const int L = MIN(LL,K);
  if (length_path <= 1) length_path=4*L;

  coeffs.setZeros();
  ind.set(-1);

  T* const pr_Gs = Gs.rawX();
  T* const pr_invGs = invGs.rawX();
  T* const pr_Ga = Ga.rawX();
  T* const pr_work = work.rawX();
  T* const pr_u = u.rawX();
  T* const pr_DtR = DtR.rawX();
  T* const pr_coeffs = coeffs.rawX();
  INTM* const pr_ind = ind.rawX();

  // Find the most correlated element
  int currentInd = pos ? DtR.max() : DtR.fmax();
  if (mode == PENALTY && abs(DtR[currentInd]) < constraint) return;
  if (mode == L2ERROR && normX < constraint) return;
  bool newAtom=true;

  int i;
  int iter=0;
  T thrs = 0;
  for (i = 0; i<L; ++i) {
    ++iter;
    if (newAtom) {
      pr_ind[i]=currentInd;
      //    cerr << "Add " << currentInd << endl;
      G.extract_rawCol(pr_ind[i],pr_Ga+i*K);
      for (int j = 0; j<=i; ++j)
        pr_Gs[i*LL+j]=pr_Ga[i*K+pr_ind[j]];

      // Update inverse of Gs
      if (i == 0) {
        pr_invGs[0]=T(1.0)/pr_Gs[0];
      } else {
        cblas_symv<T>(CblasColMajor,CblasUpper,i,T(1.0),
                      pr_invGs,LL,pr_Gs+i*LL,1,T(0.0),pr_u,1);
        const T schur =
          T(1.0)/(pr_Gs[i*LL+i]-cblas_dot<T>(i,pr_u,1,pr_Gs+i*LL,1));
        pr_invGs[i*LL+i]=schur;
        //            cblas_copy<T>(i,pr_u,1,pr_invGs+i*LL,1);
        memcpy(pr_invGs+i*LL,pr_u,i*sizeof(T));
        cblas_scal<T>(i,-schur,pr_invGs+i*LL,1);
        cblas_syr<T>(CblasColMajor,CblasUpper,i,schur,pr_u,1,
                     pr_invGs,LL);
      }
    }

    // Compute the path direction
    for (int j = 0; j<=i; ++j)
      pr_work[j]= pr_DtR[pr_ind[j]] > 0 ? T(1.0) : T(-1.0);
    cblas_symv<T>(CblasColMajor,CblasUpper,i+1,T(1.0),pr_invGs,LL,
                  pr_work,1,T(0.0),pr_u,1);

    // Compute the step on the path
    T step_max = INFINITY;
    int first_zero = -1;
    for (int j = 0; j<=i; ++j) {
      T ratio = -pr_coeffs[j]/pr_u[j];
      if (ratio > 0 && ratio <= step_max) {
        step_max=ratio;
        first_zero=j;
      }
    }
    //     PRINT_F(step_max)

    T current_correlation = abs<T>(pr_DtR[pr_ind[0]]);
    cblas_gemv<T>(CblasColMajor,CblasNoTrans,K,i+1,T(1.0),pr_Ga,
                  K,pr_u,1,T(0.0),pr_work+2*K,1);
    memcpy(pr_work+K,pr_work+2*K,K*sizeof(T));
    memcpy(pr_work,pr_work+K,K*sizeof(T));
    //      cblas_copy<T>(K,pr_work+2*K,1,pr_work+K,1);
    //     cblas_copy<T>(K,pr_work+2*K,1,pr_work,1);

    for (int j = 0; j<=i; ++j) {
      pr_work[pr_ind[j]]=INFINITY;
      pr_work[pr_ind[j]+K]=INFINITY;
    }
    for (int j = 0; j<K; ++j) {
      pr_work[j] = ((pr_work[j] < INFINITY) && (pr_work[j] > T(-1.0))) ? (pr_DtR[j]+current_correlation)/(T(1.0)+pr_work[j]) : INFINITY;
    }
    //     work.print("work");
    for (int j = 0; j<K; ++j) {
      pr_work[j+K] = ((pr_work[j+K] < INFINITY) && (pr_work[j+K] < T(1.0))) ? (current_correlation-pr_DtR[j])/(T(1.0)-pr_work[j+K]) : INFINITY;
    }
    //     work.print("work");

    if (pos) {
      for (int j = 0; j<K; ++j) {
        pr_work[j]=INFINITY;
      }
    }
    //     work.print("work");
    //     coeffs.print("coeffs");
    int index = cblas_iamin<T>(2*K,pr_work,1);
    T step = pr_work[index];

    // Choose next element
    currentInd = index % K;

    // compute the coefficients of the polynome representing normX^2
    T coeff1 = 0;
    for (int j = 0; j<=i; ++j)
      coeff1 += pr_DtR[pr_ind[j]] > 0 ? pr_u[j] : -pr_u[j];
    T coeff2 = 0;
    for (int j = 0; j<=i; ++j)
      coeff2 += pr_DtR[pr_ind[j]]*pr_u[j];
    T coeff3 = normX-constraint;


    T step_max2;
    if (mode == PENALTY) {
      step_max2 = current_correlation-constraint;
    } else if (mode == L2ERROR) {
      /// L2ERROR
      const T delta = coeff2*coeff2-coeff1*coeff3;
      step_max2 = delta < 0 ? INFINITY : (coeff2-sqrt(delta))/coeff1;
      step_max2 = MIN(current_correlation,step_max2);
    } else {
      /// L1COEFFS
      step_max2 = coeff1 < 0 ? INFINITY : (constraint-thrs)/coeff1;
      step_max2 = MIN(current_correlation,step_max2);
    }
    step = MIN(MIN(step,step_max2),step_max);
    if (step == INFINITY) break; // stop the path

    // Update coefficients
    cblas_axpy<T>(i+1,step,pr_u,1,pr_coeffs,1);

    if (pos) {
      for (int j = 0; j<i+1; ++j)
        if (pr_coeffs[j] < 0) pr_coeffs[j]=0;
    }

    // Update correlations
    cblas_axpy<T>(K,-step,pr_work+2*K,1,pr_DtR,1);

    // Update normX
    normX += coeff1*step*step-2*coeff2*step;

    // Update norm1
    thrs += step*coeff1;

    if (path) {
      for (int k = 0; k<=i; ++k)
        path[iter*K+ind[k]]=pr_coeffs[k];
    }

    // Choose next action

    if (step == step_max) {
      //   cerr << "Remove " << pr_ind[first_zero] << endl;
      /// Downdate, remove first_zero
      /// Downdate Ga, Gs, invGs, ind, coeffs
      for (int j = first_zero; j<i; ++j) {
        cblas_copy<T>(K,pr_Ga+(j+1)*K,1,pr_Ga+j*K,1);
        pr_ind[j]=pr_ind[j+1];
        pr_coeffs[j]=pr_coeffs[j+1];
      }
      pr_ind[i]=-1;
      pr_coeffs[i]=0;
      for (int j = first_zero; j<i; ++j) {
        cblas_copy<T>(first_zero,pr_Gs+(j+1)*LL,1,pr_Gs+j*LL,1);
        cblas_copy<T>(i-first_zero,pr_Gs+(j+1)*LL+first_zero+1,1,
                      pr_Gs+j*LL+first_zero,1);
      }
      const T schur = pr_invGs[first_zero*LL+first_zero];
      cblas_copy<T>(first_zero,pr_invGs+first_zero*LL,1,pr_u,1);
      cblas_copy<T>(i-first_zero,pr_invGs+(first_zero+1)*LL+first_zero,LL,
                    pr_u+first_zero,1);
      for (int j = first_zero; j<i; ++j) {
        cblas_copy<T>(first_zero,pr_invGs+(j+1)*LL,1,pr_invGs+j*LL,1);
        cblas_copy<T>(i-first_zero,pr_invGs+(j+1)*LL+first_zero+1,1,
                      pr_invGs+j*LL+first_zero,1);
      }
      cblas_syr<T>(CblasColMajor,CblasUpper,i,T(-1.0)/schur,
                   pr_u,1,pr_invGs,LL);
      newAtom=false;
      i=i-2;
    } else {
      newAtom=true;
    }
    if ((iter >= length_path-1) || abs(step) < 1e-15 ||
        step == step_max2 || (normX < 1e-15) ||
        (i == (L-1)) ||
        (mode == L2ERROR && normX - constraint < 1e-15) ||
        (mode == L1COEFFS && (constraint-thrs < 1e-15))) {
      break;
    }
  }
}

/// Auxiliary function for lasso
template <typename T>
void coreLARS2W(Vector<T>& DtR, const AbstractMatrix<T>& G,
                Matrix<T>& Gs,
                Matrix<T>& Ga,
                Matrix<T>& invGs,
                Vector<T>& u,
                Vector<T>& coeffs,
                const Vector<T>& weights,
                Vector<INTM>& ind,
                Matrix<T>& work,
                T& normX,
                const constraint_type mode,
                const T constraint,
                const bool pos) {
  const int LL = Gs.n();
  const int K = G.n();
  const int L = MIN(LL,K);
  coeffs.setZeros();
  ind.set(-1);

  T* const pr_Gs = Gs.rawX();
  T* const pr_invGs = invGs.rawX();
  T* const pr_Ga = Ga.rawX();
  //  T* const pr_G = G.rawX();
  T* const pr_work = work.rawX();
  T* const pr_u = u.rawX();
  T* const pr_DtR = DtR.rawX();
  T* const pr_coeffs = coeffs.rawX();
  T* const pr_weights = weights.rawX();
  INTM* const pr_ind = ind.rawX();

  DtR.div(weights);

  // Find the most correlated element
  int currentInd = pos ? DtR.max() : DtR.fmax();
  if (mode == PENALTY && abs(DtR[currentInd]) < constraint) return;
  if (mode == L2ERROR && normX < constraint) return;
  bool newAtom=true;

  int i;
  int iter=0;
  T thrs = 0;
  for (i = 0; i<L; ++i) {
    ++iter;
    if (newAtom) {
      pr_ind[i]=currentInd;
      // Update upper part of Gs and Ga
      G.extract_rawCol(pr_ind[i],pr_Ga+i*K);
      for (int j = 0; j<=i; ++j)
        pr_Gs[i*LL+j]=pr_Ga[i*K+pr_ind[j]];

      // Update inverse of Gs
      if (i == 0) {
        pr_invGs[0]=T(1.0)/pr_Gs[0];
      } else {
        cblas_symv<T>(CblasColMajor,CblasUpper,i,T(1.0),
                      pr_invGs,LL,pr_Gs+i*LL,1,T(0.0),pr_u,1);
        const T schur =
          T(1.0)/(pr_Gs[i*LL+i]-cblas_dot<T>(i,pr_u,1,pr_Gs+i*LL,1));
        pr_invGs[i*LL+i]=schur;
        cblas_copy<T>(i,pr_u,1,pr_invGs+i*LL,1);
        cblas_scal<T>(i,-schur,pr_invGs+i*LL,1);
        cblas_syr<T>(CblasColMajor,CblasUpper,i,schur,pr_u,1,
                     pr_invGs,LL);
      }
    }

    // Compute the path direction
    for (int j = 0; j<=i; ++j)
      pr_work[j]= pr_DtR[pr_ind[j]] > 0 ? weights[pr_ind[j]] : -weights[pr_ind[j]];
    cblas_symv<T>(CblasColMajor,CblasUpper,i+1,T(1.0),pr_invGs,LL,
                  pr_work,1,T(0.0),pr_u,1);

    // Compute the step on the path
    T step_max = INFINITY;
    int first_zero = -1;
    for (int j = 0; j<=i; ++j) {
      T ratio = -pr_coeffs[j]/pr_u[j];
      if (ratio > 0 && ratio <= step_max) {
        step_max=ratio;
        first_zero=j;
      }
    }

    T current_correlation = abs<T>(pr_DtR[pr_ind[0]]);
    cblas_gemv<T>(CblasColMajor,CblasNoTrans,K,i+1,T(1.0),pr_Ga,
                  K,pr_u,1,T(0.0),pr_work+2*K,1);
    vDiv<T>(K,pr_work+2*K,pr_weights,pr_work+2*K);
    cblas_copy<T>(K,pr_work+2*K,1,pr_work+K,1);
    cblas_copy<T>(K,pr_work+2*K,1,pr_work,1);

    for (int j = 0; j<=i; ++j) {
      pr_work[pr_ind[j]]=INFINITY;
      pr_work[pr_ind[j]+K]=INFINITY;
    }
    for (int j = 0; j<K; ++j) {
      pr_work[j] = ((pr_work[j] < INFINITY) && (pr_work[j] > T(-1.0))) ? (pr_DtR[j]+current_correlation)/(T(1.0)+pr_work[j]) : INFINITY;
    }
    for (int j = 0; j<K; ++j) {
      pr_work[j+K] = ((pr_work[j+K] < INFINITY) && (pr_work[j+K] < T(1.0))) ? (current_correlation-pr_DtR[j])/(T(1.0)-pr_work[j+K]) : INFINITY;
    }

    if (pos) {
      for (int j = 0; j<K; ++j) {
        pr_work[j]=INFINITY;
      }
    }
    int index = cblas_iamin<T>(2*K,pr_work,1);
    T step = pr_work[index];
    // Choose next element
    currentInd = index % K;

    // compute the coefficients of the polynome representing normX^2
    T coeff1 = 0;
    for (int j = 0; j<=i; ++j)
      coeff1 += pr_DtR[pr_ind[j]] > 0 ? pr_weights[pr_ind[j]]*pr_u[j] :
      -pr_weights[pr_ind[j]]*pr_u[j];
    T coeff2 = 0;
    for (int j = 0; j<=i; ++j)
      coeff2 += pr_DtR[pr_ind[j]]*pr_u[j]*pr_weights[pr_ind[j]];
    T coeff3 = normX-constraint;

    T step_max2;
    if (mode == PENALTY) {
      step_max2 = current_correlation-constraint;
    } else if (mode == L2ERROR) {
      /// L2ERROR
      const T delta = coeff2*coeff2-coeff1*coeff3;
      step_max2 = delta < 0 ? INFINITY : (coeff2-sqrt(delta))/coeff1;
    } else {
      /// L1COEFFS
      step_max2 = coeff1 < 0 ? INFINITY : (constraint-thrs)/coeff1;
    }
    step = MIN(MIN(step,step_max2),step_max);

    if (step == INFINITY) break; // stop the path

    // Update coefficients
    cblas_axpy<T>(i+1,step,pr_u,1,pr_coeffs,1);

    // Update correlations
    cblas_axpy<T>(K,-step,pr_work+2*K,1,pr_DtR,1);

    // Update normX
    normX += coeff1*step*step-2*coeff2*step;

    // Update norm1
    thrs += step*coeff1;

    if (step == step_max) {
      /// Downdate, remove first_zero
      /// Downdate Ga, Gs, invGs, ind, coeffs
      for (int j = first_zero; j<i; ++j) {
        cblas_copy<T>(K,pr_Ga+(j+1)*K,1,pr_Ga+j*K,1);
        pr_ind[j]=pr_ind[j+1];
        pr_coeffs[j]=pr_coeffs[j+1];
      }
      pr_ind[i]=-1;
      pr_coeffs[i]=0;
      for (int j = first_zero; j<i; ++j) {
        cblas_copy<T>(first_zero,pr_Gs+(j+1)*LL,1,pr_Gs+j*LL,1);
        cblas_copy<T>(i-first_zero,pr_Gs+(j+1)*LL+first_zero+1,1,
                      pr_Gs+j*LL+first_zero,1);
      }
      const T schur = pr_invGs[first_zero*LL+first_zero];
      cblas_copy<T>(first_zero,pr_invGs+first_zero*LL,1,pr_u,1);
      cblas_copy<T>(i-first_zero,pr_invGs+(first_zero+1)*LL+first_zero,LL,
                    pr_u+first_zero,1);
      for (int j = first_zero; j<i; ++j) {
        cblas_copy<T>(first_zero,pr_invGs+(j+1)*LL,1,pr_invGs+j*LL,1);
        cblas_copy<T>(i-first_zero,pr_invGs+(j+1)*LL+first_zero+1,1,
                      pr_invGs+j*LL+first_zero,1);
      }
      cblas_syr<T>(CblasColMajor,CblasUpper,i,T(-1.0)/schur,
                   pr_u,1,pr_invGs,LL);
      newAtom=false;
      i=i-2;
    } else {
      newAtom=true;
    }
    // Choose next action
    if (iter > 4*L || abs(step) < 1e-10 ||
        step == step_max2 || (normX < 1e-10) ||
        (i == (L-1)) ||
        (mode == L2ERROR && normX - constraint < 1e-10) ||
        (mode == L1COEFFS && (constraint-thrs < 1e-10))) {
      break;
    }
  }
}


//// spmas.h
template <typename T>
SpMatrix<T> *_lassoD(Matrix<T> *X, Matrix<T> *D,Matrix<T> **path,bool return_reg_path,
                     int L, const T constraint, const T lambda2, constraint_type mode,
                     const bool pos, const bool ols, const int numThreads,
                     int max_length_path,const bool verbose, bool cholevsky)
  throw(const char *)
  {
    SpMatrix<T> *alpha = new SpMatrix<T>();
    int n = X->m();
    int nD = D->m();
    int K = D->n();
    if (n != nD)
      throw("lasso : incompatible matrix dimensions");
    if(L < 0) L = K;
    if(max_length_path < 0) max_length_path = 4 * L;
    if (L> n && !(mode == PENALTY && isZero(constraint) && !pos && lambda2 > 0)) {
      if (verbose)
        // printf("L is changed to %d\n",n);
      L=n;
    }
    if (L > K) {
      if (verbose)
        // printf("L is changed to %d\n",K);
      L=K;
    }
    if(return_reg_path)
      *path = new Matrix<T>(K,max_length_path);
    else
      *path = NULL;
    if(ols) cholevsky = ols;
    if (cholevsky) {
      lasso((Matrix<T> &)(*X),(Matrix<T> &)(*D),(SpMatrix<T> &)(*alpha),L,constraint,lambda2,mode,pos,ols,numThreads,*path,max_length_path);
    } else {
      lasso2((Matrix<T> &)(*X),(Matrix<T> &)(*D),(SpMatrix<T> &)(*alpha),L,constraint,lambda2,mode,pos,numThreads,*path,max_length_path);
    }
    return alpha;
  }

template <typename T>
SpMatrix<T> *_lassoQq(Matrix<T> *X, Matrix<T> *Q, Matrix<T> *q,Matrix<T> **path,bool return_reg_path,
                      int L, const T constraint, const T lambda2, constraint_type mode,
                      const bool pos, const bool ols, const int numThreads,
                      int max_length_path,const bool verbose, bool cholevsky)
  throw(const char *)
  // lambda2 is ignored
  {
    SpMatrix<T> *alpha = new SpMatrix<T>();
    int n = X->m();
    int M = X->n();
    int K1 = Q->m();
    int K2 = Q->n();
    if(K1 != K2)
      throw("lasso : Q must be square");
    int K = K1;
    int K3 = q->m();
    int M2 = q->n();
    if (K1 != K3 || M != M2)
      throw("lasso : incompatible matrix dimensions");

    if(L < 0) L = K1;
    if(max_length_path < 0) max_length_path = 4 * L;
    if (L> n && !(mode == PENALTY && isZero(constraint) && !pos && lambda2 > 0)) {
      if (verbose)
        // printf("L is changed to %d\n",n);
      L=n;
    }
    if (L > K) {
      if (verbose)
        // printf("L is changed to %d\n",K);
      L=K;
    }
    if(return_reg_path)
      *path = new Matrix<T>(K,max_length_path);
    else
      *path = NULL;
    if(ols) cholevsky = ols;
    if (cholevsky)
      lasso((Data<T> &)(*X),(AbstractMatrix<T> &)(*Q),(AbstractMatrix<T> &)(*q),(SpMatrix<T> &)(*alpha),L,constraint,mode,pos,ols,numThreads,*path,max_length_path);
    else
      lasso2((Data<T> &)(*X),(AbstractMatrix<T> &)(*Q),(AbstractMatrix<T> &)(*q),(SpMatrix<T> &)(*alpha),L,constraint,mode,pos,numThreads,*path,max_length_path);
    return alpha;
  }

template <typename T>
SpMatrix<T> *_lassoMask(Matrix<T> *X, Matrix<T> *D,Matrix<bool> *B,
                        int L, const T constraint, const T lambda2, constraint_type mode,
                        const bool pos, const int numThreads,bool verbose)
  throw(const char *)
  {
    SpMatrix<T> *alpha = new SpMatrix<T>();
    int n = X->m();
    int nD = D->m();
    int K = D->n();
    if (n != nD)
      throw("lassoMask : incompatible matrix dimensions");
    if(L < 0) L = K;
    if (L> n && !(mode == PENALTY && isZero(constraint) && !pos && lambda2 > 0)) {
      if (verbose)
        // printf("L is changed to %d\n",n);
      L=n;
    }
    if (L > K) {
      if (verbose)
        // printf("L is changed to %d\n",K);
      L=K;
    }
    lasso_mask((Matrix<T> &)(*X),(Matrix<T> &)(*D),(SpMatrix<T> &)(*alpha),(Matrix<bool> &)(*B),L,constraint,lambda2,mode,pos,numThreads);
    return alpha;
  }

template <typename T>
SpMatrix<T> *_lassoWeighted(Matrix<T> *X, Matrix<T> *D,Matrix<T> *W,
                            int L, const T constraint, constraint_type mode,
                            const bool pos, const int numThreads,bool verbose)
  throw(const char *)
  {
    SpMatrix<T> *alpha = new SpMatrix<T>();
    int n = X->m();
    int M = X->n();
    int nD = D->m();
    int K = D->n();
    if (n != nD)
      throw("lassoWeighted : incompatible matrix dimensions");
    if(L < 0) L = K;
    if (L> n ) {
      if (verbose)
        // printf("L is changed to %d\n",n);
      L=n;
    }
    if (L > K) {
      if (verbose)
        // printf("L is changed to %d\n",K);
      L=K;
    }
    int KK = W->m();
    int MM = W->n();
    if (K != KK || M != MM)
      throw("lassoWeighted : inconsistent dimensions of matrix W");

    lassoWeight((Matrix<T> &)(*X),(Matrix<T> &)(*D),(Matrix<T> &)(*W),(SpMatrix<T> &)(*alpha),L,constraint,mode,pos,numThreads);
    return alpha;
  }




#endif /* myDECOMP_h */
