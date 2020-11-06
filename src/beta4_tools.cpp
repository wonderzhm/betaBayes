#include "beta4_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

// log likelihood given data and parameters 
void beta4_mean_loglik(
    const Rcpp::NumericVector& y, const Rcpp::NumericVector& Xbeta, double phi,
    double th1, double th2, int link, double& ll){
  int n = y.size();
  double mi=0;
  ll = n*(Rf_lgammafn(phi) - (phi-1)*std::log(th2-th1));
  for(int i=0; i<n; ++i){
    mi = ilinkf(Xbeta[i], link);
    ll += (phi*mi-1)*std::log(y[i]-th1) + (phi*(1-mi)-1)*std::log(th2-y[i]);
    ll += -Rf_lgammafn(phi*mi) - Rf_lgammafn(phi*(1-mi));
  }
}
void beta4_mode_loglik(
    const Rcpp::NumericVector& y, const Rcpp::NumericVector& Xbeta, double phi,
    double th1, double th2, int link, double& ll){
  int n = y.size();
  double mi=0;
  ll = n*(Rf_lgammafn(phi+2) - (phi+1)*std::log(th2-th1));
  for(int i=0; i<n; ++i){
    mi = ilinkf(Xbeta[i], link);
    ll += (phi*mi)*std::log(y[i]-th1) + (phi*(1-mi))*std::log(th2-y[i]);
    ll += -Rf_lgammafn(phi*mi+1) - Rf_lgammafn(phi*(1-mi)+1);
  }
}
// Calculate loglikelihood for each obervation i
arma::vec beta4_mean_logliki(
    const Rcpp::NumericVector& y, const Rcpp::NumericVector& Xbeta, double phi, 
    double th1, double th2, int link){
  int n = y.size();
  double mi=0;
  arma::vec res(n);
  for(int i=0; i<n; ++i){
    mi = ilinkf(Xbeta[i], link);
    res[i] = Rf_lgammafn(phi) - (phi-1)*std::log(th2-th1) 
      + (phi*mi-1)*std::log(y[i]-th1) + (phi*(1-mi)-1)*std::log(th2-y[i])
      - Rf_lgammafn(phi*mi) - Rf_lgammafn(phi*(1-mi));
  } 
  return(res);
}
arma::vec beta4_mode_logliki(
    const Rcpp::NumericVector& y, const Rcpp::NumericVector& Xbeta, double phi, 
    double th1, double th2, int link){
  int n = y.size();
  double mi=0;
  arma::vec res(n);
  for(int i=0; i<n; ++i){
    mi = ilinkf(Xbeta[i], link);
    res[i] = Rf_lgammafn(phi+2) - (phi+1)*std::log(th2-th1) 
      + (phi*mi)*std::log(y[i]-th1) + (phi*(1-mi))*std::log(th2-y[i])
      - Rf_lgammafn(phi*mi+1) - Rf_lgammafn(phi*(1-mi)+1);
  } 
  return(res);
}
