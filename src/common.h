#ifndef _common_h
#define _common_h

#include <RcppArmadillo.h>

#define ESMALL 1e-15  /* small number */
#define ELARGE 1e+15 /* large number */
#define SYSMIN 1e-305  /* small number */
#define SYSMAX 1e+305 /* large number */
#define LOGSYSMAX 702.28845336318397585 /* large number */
#define LOGSYSMIN -702.28845336318397585 /* large number */
typedef Rcpp::NumericMatrix::iterator mat_iterator;
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
int randWrapper(const int n);

//Truncated normal N(y;mu,s)I(a<y<b)
double rtexp(double a, double b); // Rejection algorithm with a truncated expoential proposal for N(0,1)I(a<x<b) when a is very large: |a| < |b|
double trun_rnorm(const double mu, const double s, double a, double b);

// generate multivariate normal (mu, sigma)
arma::vec mvrnorm(arma::vec mu, arma::mat sigma);

// density of multivariate normal (mu, sigma)
double mvdnorm(arma::vec x, arma::vec mu, arma::mat sigma, bool logd);

// generate Wishart random matrices
arma::mat rwish(arma::mat Sig, int n);

// sample(Nseq, prob=w), where Nseq is n-dim vector, and w.siz()=n
int sample(Rcpp::IntegerVector Nseq, Rcpp::NumericVector w);

// calculate qnorm(x) for a vector of x
arma::vec qnormvec(arma::vec x);

// calculate pnorm(x) for a vector of x
arma::vec pnormvec(arma::vec x);

// inverse link function
double ilinkf(double x, int link);
// link function
double linkf(double y, int link);

// logit link function
double ilogit(double x);
double logit(double y);

// probit link functions
double iprobit(double x);
double probit(double y);

// complementary loglog link function
double icloglog(double x);
double cloglog(double y);

// loglog link function
double iloglog(double x);
double loglog(double y);

#endif
