#include <Rcpp.h>
#include <R_ext/Applic.h>
#include <math.h>
#include <cmath>
#include <map>


using namespace Rcpp;
class OptData
{
public:
  NumericVector a;
  double rho;
  double lambda;
  double c;
};


void pSCAD(double* pscad, int m, double* b, double lam, double c) {
  // m is the length of b
  double pMax=(c+1)*pow(lam,2)/2;
  for(int idx=0;idx<m;++idx){
    double bb=std::abs(b[idx]);
    if(bb<=lam)
    {
      pscad[idx]=lam*bb;
    }else if(bb>=c*lam)
    {
      pscad[idx]=pMax;
    }else
    {
      pscad[idx]=-(pow(bb,2)-2*c*bb*lam+pow(lam,2))/(2*(c-1));
    }
  }

}


/////// pSCAD' function  /////////////////

void p1SCAD(double* p1scad, int m, double* b, double lam, double c)
{
  for(int idx=0;idx<m;++idx)
  {
    double bb=std::abs(b[idx]);

    if(bb>c*lam)
    {
      p1scad[idx]=0;
    // }else if(bb==0)
    // {
    //   p1scad[idx]=0;
    }else if(bb<=lam)
    {
      p1scad[idx]=b[idx]<0 ? -lam:lam;
    }else
    {
      p1scad[idx]=b[idx]<0 ? -(c*lam-bb)/(c-1):(c*lam-bb)/(c-1);
    }
  }
}



inline double vecSum(double* vec, int len)
{
  double sum=0;
  for (int idx=0; idx<len; idx++)
  {
    sum+=vec[idx];
  }
  return sum;
}



double fun(int m, double *par, void *ex)
{
  OptData* dat = (OptData*) ex;

  double* pscad = new double[m];
  pSCAD(pscad, m, par, dat->lambda, dat->c);
  double obj=log(1+vecSum(pscad,m));
  for(int idx=0;idx<m; idx++)
    obj+=(dat->rho)*pow((dat->a[idx]-par[idx]),2)/2;

  delete [] pscad;
  return obj;
}

// //' @export
// //[[Rcpp::export]]
// double test_obj(NumericVector a,NumericVector b, double rho, double lambda, double c)
// {
//   OptData dat;
//   dat.a=a;
//   dat.rho=rho;
//   dat.lambda=lambda;
//   dat.c=c;
//
//   double out= fun(a.size(), b.begin(), &dat);
//   return(out);
// }
//




void gradient(int m, double *par, double *gr, void *ex)
{
  OptData* dat = (OptData*) ex;

  double* pscad = new double[m];
  pSCAD(pscad, m, par, dat->lambda, dat->c);

  double* p1scad = new double[m];
  p1SCAD(p1scad, m, par, dat->lambda, dat->c);

  double c1=1/(1+vecSum(pscad,m));

  for(int idx=0;idx<m; idx++)
    gr[idx]=dat->rho*(par[idx]-dat->a[idx])+p1scad[idx]*c1;

  delete[] pscad;
  delete[] p1scad;
}



NumericVector D2optim_vec(NumericVector a, double lambda, double c, double rho)
{
  OptData dat;
  dat.a =Rcpp::abs(a);
  dat.rho = rho;
  dat.lambda = lambda;
  dat.c = c;

  // all parameters
  int m = a.size();
  int lmm = 5;
  NumericVector x(m);

  double* lower = new double[m];
  double* upper = new double[m];
  int* nbd = new int[m];
  for (int idx=0;idx<m;idx++)
  {
    nbd[idx]=2;
    lower[idx]=0;
    upper[idx]=std::abs(a[idx]);
    // if(a[idx]>0)
    // {
    //   lower[idx]=0; upper[idx]=a[idx];
    // }
    // else
    // {
    //   lower[idx]=a[idx];upper[idx]=0;
    // }
  }
  double Fmin;
  int fail;

  double factr=1e+07;
  double pgtol=0;
  int fncount;
  int grcount;
  int maxit=100;
  char msg[1024];
  int trace=0;
  int nREPORT=10; // used when trace is positive

  lbfgsb(m, lmm, x.begin(), lower, upper, nbd,
         &Fmin, fun, gradient, &fail, &dat, factr, pgtol,
         &fncount, &grcount, maxit, msg, trace, nREPORT );


         //correct the numeric error
         double *pscad = new double[m];
         pSCAD(pscad, m, x.begin(), lambda, c);
         for(int idx=0;idx<m; idx++)
         {
           if(std::abs(x[idx])< (1e-4)*lambda)
           {
             x[idx]=0;
           }else if (a[idx]<0)
           {
             x[idx]=-x[idx];
           }

         }

         delete [] pscad;
         delete [] lower;
         delete [] upper;
         delete [] nbd;

         // return List::create(x, Fmin);
         return x;
}

////' @export
//[[Rcpp::export]]
NumericMatrix D2optim_mat(NumericMatrix A, double lambda, double c, double rho)
{
  NumericMatrix B(A.nrow(), A.ncol());
  // NumericVector obj(A.ncol());
  // int m = A.nrow();
  // for (int idx=0; idx<A.ncol(); idx++)
  // {
  //   NumericVector a_nz, nz_idx, b_nz;
  //   NumericVector a=A(_,idx);

    // if a_j less than thred, a_j=0
    // double *pscad=new double[m];
    // pSCAD(pscad, m, a.begin(), lambda, c);
    // double thred=lambda/((1+vecSum(pscad,m))*rho);
    // delete[] pscad;


  //   for(int j=0;j<m;j++)
  //   {
  //     if(std::abs(a(j)) > 0)
  //       {
  //         nz_idx.push_back(j);
  //         a_nz.push_back(a(j));
  //       }
  //   }
  //   b_nz=D2optim_vec(a_nz,lambda,c, rho, &obj(idx));
  //   for(int j=0; j<nz_idx.size(); j++)
  //   {
  //     B(nz_idx(j),idx)=b_nz(j);
  //   }
  // }
    for(int idx=0; idx<A.ncol(); idx++)
    {
      B(_,idx)=D2optim_vec(A(_,idx),lambda,c, rho);

    }

  // return List::create(Named("B")=B,Named("obj")=obj);
  return B;
}

//[[Rcpp::export]]
double pgscad(NumericMatrix D, double c, double lambda)
{
  int m=D.nrow();
  int p=D.ncol();
  double* pscad = new double[m];
  double obj=0;
  for(int idx=0; idx<p; idx++)
  {
    pSCAD(pscad, m, D(_,idx).begin(), lambda, c);
    obj += log(1+vecSum(pscad,m));
  }
  delete[] pscad;
  return(obj);
}


