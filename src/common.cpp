#include "common.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma ;
using namespace Rcpp ;
using namespace std;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
int randWrapper(const int n) { 
  return floor(unif_rand()*n); 
}

//Truncated normal N(y;mu,s)I(a<y<b); where b can be R_PosInf
// Rejection algorithm with a truncated expoential proposal for N(0,1)I(a<x<b) when a is very large: |a| < |b|
double rtexp(double a, double b){
  int stop = false;
  double twoasp = 2*std::pow(a,2);
  double expab = std::exp(-a*(b-a)) - 1;
  double z, e;
  while(!stop){
    R_CheckUserInterrupt();
    z = std::log(1 + unif_rand()*expab);
    e = -std::log(unif_rand());
    stop = (twoasp*e > std::pow(z,2));
  }
  return (a - z/a);
}
double trun_rnorm(const double mu, const double s, double a, double b){
  double xmin = -2.00443204036;                 // Left bound
  double xmax =  3.48672170399;                 // Right bound
  int stop = false;
  double r;
  //if( mu+ELARGE<0 ) return(a);
  //scalling
  if(mu!=0 || s!=1){
    a=(a-mu)/s;
    b=(b-mu)/s;
  }
  // Check if a < b
  if(a>=b){
    Rprintf( "*** B must be greater than A ! ***" ); return(NA_REAL);
  }
  else if(std::abs(a)>std::abs(b)) r = -trun_rnorm(0, 1, -b, -a);
  // If a in the right tail (a > xmax), use rejection algorithm with a truncated exponential proposal  
  else if(a>xmax) r = rtexp(a,b);
  // If a in the left tail (a < xmin), use rejection algorithm with a Gaussian proposal
  else if(a<xmin){
    while(!stop){
      R_CheckUserInterrupt();
      r = norm_rand();
      stop = (r>=a) && (r<=b);
    }
  }  
  // In other cases (xmin < a < xmax)
  else{
    double CDFa = Rf_pnorm5(a, 0, 1.0, true, false);
    double CDFb = Rf_pnorm5(b, 0, 1.0, true, false);
    double u = unif_rand();
    double CDFxi = CDFa + u*(CDFb - CDFa);
    r = Rf_qnorm5(CDFxi, 0, 1, true, false);
  }
  // Scaling
  if(mu!=0 || s!=1)
  r = r*s + mu;
  return r;
}

// generate multivariate normal (mu, sigma)
arma::vec mvrnorm(arma::vec mu, arma::mat sigma) {
  int ncols = mu.size();
  arma::vec Y(ncols);
  for (int i=0; i<ncols; i++){
    Y[i] = norm_rand();
  }
  arma::mat temp = ((arma::chol(sigma)).t())*Y;
  arma::vec res = mu + temp.col(0);
  return( res );
}

// density of multivariate normal (mu, sigma)
double mvdnorm(arma::vec x, arma::vec mu, arma::mat sigma, bool logd=true) { 
  int xdim = x.size();
  arma::mat rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(arma::log(rooti.diag()));
  double constants = -(double)(xdim)*0.5*std::log(2.0*M_PI);
  arma::vec z = rooti*( x - mu) ;    
  double res = constants - 0.5*arma::sum(z%z) + rootisum;     
  if (logd == false) {
    res = std::exp(res);
  }
  return(res);
}

// generate Wishart random matrices
arma::mat rwish(arma::mat Sig, int n) {
  int ncols = Sig.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  arma::mat X = Y * arma::chol(Sig);
  return( X.t()*X );
}

// sample(Nseq, prob=w), where Nseq is n-dim vector, and w.siz()=n
int sample(Rcpp::IntegerVector Nseq, Rcpp::NumericVector w){
  int k = 0;
  double u = unif_rand();;
  double cdf = w[0];
  while(u>cdf){
    cdf += w[k];
    ++k;
  }
  return (Nseq[k]);
}

// calculate qnorm(x) for a vector of x
arma::vec qnormvec(arma::vec x){
  int n = x.size();
  arma::vec res(n);
  for (int i=0; i<n; ++i){
    res[i] = std::min( Rf_qnorm5(x[i], 0, 1, true, false), 8.209536);
  }
  return (res);
}

// calculate pnorm(x) for a vector of x
arma::vec pnormvec(arma::vec x){
  int n = x.size();
  arma::vec res(n);
  for (int i=0; i<n; ++i){
    res[i] = Rf_pnorm5(x[i], 0, 1, true, false);
  }
  return (res);
}

// inverse link function
double ilinkf(double x, int link){
  if(link==1){
    // ilogit
    return ( 1.0/(1.0+std::exp(-x)) );
  }else if(link==2){
    // iprobit
    return(Rf_pnorm5(x, 0, 1.0, true, false));
  }else if(link==3) {
    // icloglog
    return ( 1 - std::exp(-std::exp(x)) );
  }else{
    // iloglog
    return ( std::exp(-std::exp(-x)) );
  }
}

// link function
double linkf(double y, int link){
  if(link==1){
    // logit
    return ( std::log(y/(1.0-y)) );
  }else if(link==2){
    // probit
    return(Rf_qnorm5(y, 0, 1.0, true, false));
  }else if(link==3) {
    // cloglog
    return ( std::log(-std::log(1-y)) );
  }else{
    // loglog
    return ( -std::log(-std::log(y)) );
  }
}

// logit link function
double ilogit(double x){
  return ( 1.0/(1.0+std::exp(-x)) );
}
double logit(double y){
  return ( std::log(y/(1.0-y)) );
}

// probit link functions
double iprobit(double x){
  return(Rf_pnorm5(x, 0, 1.0, true, false));
}
double probit(double y){
  return(Rf_qnorm5(y, 0, 1.0, true, false));
}

// complementary loglog link function
double icloglog(double x){
  return ( 1 - std::exp(-std::exp(x)) );
}
double cloglog(double y){
  return ( std::log(-std::log(1-y)) );
}

// loglog link function
double iloglog(double x){
  return ( std::exp(-std::exp(-x)) );
}
double loglog(double y){
  return ( -std::log(-std::log(y)) );
}

