#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix imageSplit(NumericMatrix A,int m, int n, int stepsize)
{

  int mm = A.nrow();
  int nn = A.ncol();

  int M = m * n;
  int N =std::ceil( double(mm - m + 1)/stepsize)*std::ceil(double(nn -n + 1)/stepsize);
  // int N=std::ceil(2.3);
  NumericMatrix  B(M,N);


  double *po = B.begin();
  double *pi = A.begin();
  for(int j = 0; j <= nn - n;j+=stepsize) {
    for(int i = 0;i <= mm - m; i+=stepsize) {
      for(int kj = j;kj < j + n;kj++) {
        int kj1 = kj;
        for(int ki = i;ki < i + m;ki++) {
          *po++ = *(pi + ki + kj1 * mm);
        }
      }
    }
  }
  return B;

}



// [[Rcpp::export]]
NumericMatrix imageReCon(NumericMatrix B, int mm, int nn, int m, int n,int stepsize)
{
  // int M = m * n;
  // int N = (mm - m + 1) * (nn -n + 1);

  NumericMatrix MatSum(mm,nn);
  NumericMatrix MatCount(mm,nn);


  double *pb = B.begin();
  double *ps = MatSum.begin();
  double *pc = MatCount.begin();

  for(int j = 0; j <= nn - n;j+=stepsize) {
    for(int i = 0;i <= mm - m; i+=stepsize) {
      for(int kj = j;kj < j + n;kj++) {
        int kj1 = kj;
        for(int ki = i;ki < i + m;ki++) {
          *(ps + ki + kj1 * mm) += *pb++ ;
          *(pc + ki + kj1 * mm) += 1 ;
        }
      }
    }
  }

  NumericMatrix A(mm,nn);
  ps = MatSum.begin();
  pc = MatCount.begin();
  double *pa = A.begin();
  for(int j = 0; j < nn;j++) {
    for(int i = 0;i < mm; i++) {
      *(pa + i + j*mm) = *(ps + i + j*mm) / *(pc +i +j*mm) ;
    }
  }


  return A;

}


