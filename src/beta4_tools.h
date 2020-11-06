#ifndef _beta4_tools_h
#define _beta4_tools_h

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "common.h"

// log likelihood given data and parameters 
void beta4_mean_loglik(const Rcpp::NumericVector& y, const Rcpp::NumericVector& Xbeta, double phi, 
                       double th1, double th2, int link, double& ll);
void beta4_mode_loglik(const Rcpp::NumericVector& y, const Rcpp::NumericVector& Xbeta, double phi, 
                       double th1, double th2, int link, double& ll);
// Calculate loglikelihood for each obervation i
arma::vec beta4_mean_logliki(const Rcpp::NumericVector& y, const Rcpp::NumericVector& Xbeta, double phi, 
                             double th1, double th2, int link);
arma::vec beta4_mode_logliki(const Rcpp::NumericVector& y, const Rcpp::NumericVector& Xbeta, double phi, 
                             double th1, double th2, int link);
#endif
