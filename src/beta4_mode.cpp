#include "beta4_tools.h"

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace std;

RcppExport SEXP 
  beta4_mode_reg(
    SEXP nburn_, SEXP nsave_, SEXP nskip_, SEXP ndisplay_, SEXP y_, SEXP X_, 
    SEXP beta_, SEXP phi_, SEXP theta_, SEXP phia0_, SEXP phib0_, 
    SEXP th1a0_, SEXP th1b0_, SEXP th2a0_, SEXP th2b0_, 
    SEXP Vhat_, SEXP beta0_, SEXP S0inv_, SEXP Shat_, SEXP l0_, SEXP adapter_, 
    SEXP link_, SEXP Xpred_
  ){
  BEGIN_RCPP
  
  // Transfer R variables into C++;
  const int nburn = as<int>(nburn_);
  const int nsave = as<int>(nsave_);
  const int nskip = as<int>(nskip_);
  const int ndisplay = as<int>(ndisplay_);
  const Rcpp::NumericVector y(y_); // n by 1;
  const Rcpp::NumericMatrix X(X_); // n by p;
  const int n = y.size();
  const int p = X.ncol();
  const int l0 = as<int>(l0_);
  const double adapter = as<double>(adapter_);
  const int link = as<int>(link_);
  const arma::mat Xpred = as<arma::mat>(Xpred_); // m by p;
  const int m = Xpred.n_rows;

  // prameters to be updated
  Rcpp::NumericVector theta(theta_); // 2 by 1
  Rcpp::NumericVector beta(beta_); // p by 1
  double phi = as<double>(phi_);
  
  // hyperparameters
  // prior of theta
  const double th1a0 = as<double>(th1a0_); 
  const double th1b0 = as<double>(th1b0_); 
  const double th2a0 = as<double>(th2a0_); 
  const double th2b0 = as<double>(th2b0_); 
  const arma::mat Vhat = as<arma::mat>(Vhat_);     // 2 by 2
  // prior of beta
  const arma::vec beta0 = as<arma::vec>(beta0_); // p by 1
  const arma::mat S0inv = as<arma::mat>(S0inv_); // p by p
  const arma::mat Shat = as<arma::mat>(Shat_);   // p by p
  // prior of phi
  const double phia0 = as<double>(phia0_);
  const double phib0 = as<double>(phib0_);
  
  // temp variables
  const int nscan = nburn + (nskip+1)*nsave;
  Rcpp::NumericVector Xbeta(n, 0.0);
  int skiptally=0; 
  int isave=0;
  int distally=0;
  
  // things to save;
  Rcpp::NumericMatrix theta_save(2, nsave);
  Rcpp::NumericMatrix beta_save(p, nsave);
  Rcpp::NumericVector phi_save(nsave);
  Rcpp::NumericMatrix mode_save(m, nsave);
  double rejtheta=0;
  double rejbeta=0;
  double rejphi=0;
  
  // Make arma objects
  arma::mat X_r(const_cast<NumericMatrix&>(X).begin(), n, p, false);
  arma::vec theta_r(theta.begin(), 2, false);
  arma::vec beta_r(beta.begin(), p, false);
  arma::vec Xbeta_r(Xbeta.begin(), n, false);
  
  // Working temp variables
  arma::mat logLik=arma::zeros<arma::mat>(n, nsave);
  arma::vec Dvalues = arma::zeros<arma::vec>(nsave);
  Rcpp::NumericVector sumtheta(2, 0.0); // 2 by 1
  Rcpp::NumericVector sumbeta(p, 0.0); // p by 1
  double sumphi=0;
  double llold, llnew;
  double ratio, uu, nn;
  
  // for theta & beta
  arma::vec thetax = arma::zeros<arma::vec>(2);
  if(th1a0<th1b0) thetax[0] = std::log((theta[0]-th1a0)/(th1b0-theta[0]));
  if(th2a0<th2b0) thetax[1] = std::log((theta[1]-th2a0)/(th2b0-theta[1]));
  int p1 = p + 1;
  arma::vec beth1(p1);
  arma::vec beth1old(p1); 
  arma::vec beth1Barold(p1);
  arma::vec beth1Barnew=arma::zeros<arma::vec>(p1);
  arma::mat beth1Snew=arma::zeros<arma::mat>(p1, p1);
  arma::mat Ip1 = ESMALL*arma::eye(p1, p1);
  arma::mat SVhat1 = arma::zeros<arma::mat>(p1, p1);
  SVhat1.submat(0, 0, p-1, p-1) = Shat;
  beth1.subvec(0, p-1) = beta_r;
  if((th1a0<th1b0)&(th2a0>=th2b0)){
    SVhat1(p,p) = Vhat(0,0);
    beth1[p] = thetax[0];
  } 
  if((th1a0>=th1b0)&(th2a0<th2b0)){
    SVhat1(p,p) = Vhat(1,1);
    beth1[p] = thetax[1];
  } 
  int p2 = p + 2;
  arma::vec beth2(p2);
  arma::vec beth2old(p2); 
  arma::vec beth2Barold(p2);
  arma::vec beth2Barnew=arma::zeros<arma::vec>(p2);
  arma::mat beth2Snew=arma::zeros<arma::mat>(p2, p2);
  arma::mat Ip2 = ESMALL*arma::eye(p2, p2);
  arma::mat SVhat2 = arma::zeros<arma::mat>(p2, p2);
  SVhat2.submat(0, 0, p-1, p-1) = Shat;
  SVhat2.submat(p, p, p+1, p+1) = Vhat;
  beth2.subvec(0, p-1) = beta_r;
  beth2.subvec(p, p+1) = thetax;
  
  // for beta
  arma::vec betaold(p);
  arma::vec beBarold(p);
  arma::vec beBarnew=arma::zeros<arma::vec>(p);
  arma::mat beSnew=arma::zeros<arma::mat>(p,p);
  arma::mat Ip = ESMALL*arma::eye(p,p);
  //double eta=0;
  //double abnormal=0;
  // for logphi
  double phibarnew=0, phibarold=0, phisnew=0, logphi=0, phiold=0, phishat=1.0/(n+0.0);
  
  RNGScope scope;
  
  // Set the Armadillo seed from R's 
  // int seed = (int)Rf_runif(0.0, 10000.0);
  // std::srand(seed);
  
  ////////////////////////////////////////////////////////////////////////
  // Start MCMC
  ////////////////////////////////////////////////////////////////////////
  Xbeta_r = X_r*beta_r;
  logphi = std::log(phi);
  
  for (int iscan=0; iscan<nscan; iscan++){
    R_CheckUserInterrupt();
    //Rprintf( "iscan = %d\n", iscan );
    //Rprintf( "lambda = %f\n", lambda );
    //Rprintf( "phi = %f\n", phi );
    
    ///////////////////////////////////////////////
    // update beta & theta as a block
    //////////////////////////////////////////////
    if((th1a0<th1b0)&(th2a0<th2b0)){
      beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llold);
      llold += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
      llold += thetax[0] - 2*std::log(1+std::exp(thetax[0]));
      llold += thetax[1] - 2*std::log(1+std::exp(thetax[1]));
      beth2old = beth2;
      if(iscan>l0){
        beth2 = mvrnorm(beth2old, beth2Snew);
      }else{
        beth2 = mvrnorm(beth2old, SVhat2);
      }
      beta_r = beth2.subvec(0, p-1);
      thetax = beth2.subvec(p, p+1);
      theta[0] = (th1a0+th1b0*std::exp(thetax[0]))/(1+std::exp(thetax[0]));
      theta[1] = (th2a0+th2b0*std::exp(thetax[1]))/(1+std::exp(thetax[1]));
      Xbeta_r = X_r*beta_r;
      beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llnew);
      llnew += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
      llnew += thetax[0] - 2*std::log(1+std::exp(thetax[0]));
      llnew += thetax[1] - 2*std::log(1+std::exp(thetax[1]));
      ratio = exp(llnew-llold);
      uu = unif_rand();
      if((uu>ratio)){
        beth2=beth2old;
        if(iscan>=nburn) {
          rejbeta+=1.0;
          rejtheta+=1.0;
        }
        beta_r = beth2.subvec(0, p-1);
        thetax = beth2.subvec(p, p+1);
        theta[0] = (th1a0+th1b0*std::exp(thetax[0]))/(1+std::exp(thetax[0]));
        theta[1] = (th2a0+th2b0*std::exp(thetax[1]))/(1+std::exp(thetax[1]));
        Xbeta_r = X_r*beta_r;
      }
      nn = iscan+1;
      beth2Barold = beth2Barnew;
      beth2Barnew = (nn)/(nn+1.0)*beth2Barold + beth2/(nn+1.0);
      beth2Snew = (nn-1.0)/nn*beth2Snew + adapter/(p2+0.0)/nn*(nn*beth2Barold*beth2Barold.t() 
                                                                 - (nn+1.0)*beth2Barnew*beth2Barnew.t()
                                                                 + beth2*beth2.t() + Ip2 );
    }else if((th1a0<th1b0)&(th2a0>=th2b0)){
      beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llold);
      llold += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
      llold += thetax[0] - 2*std::log(1+std::exp(thetax[0]));
      beth1old = beth1;
      if(iscan>l0){
        beth1 = mvrnorm(beth1old, beth1Snew);
      }else{
        beth1 = mvrnorm(beth1old, SVhat1);
      }
      beta_r = beth1.subvec(0, p-1);
      thetax[0] = beth1[p];
      theta[0] = (th1a0+th1b0*std::exp(thetax[0]))/(1+std::exp(thetax[0]));
      Xbeta_r = X_r*beta_r;
      beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llnew);
      llnew += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
      llnew += thetax[0] - 2*std::log(1+std::exp(thetax[0]));
      ratio = exp(llnew-llold);
      uu = unif_rand();
      if((uu>ratio)){
        beth1=beth1old;
        if(iscan>=nburn) {
          rejbeta+=1.0;
          rejtheta+=1.0;
        }
        beta_r = beth1.subvec(0, p-1);
        thetax[0] = beth1[p];
        theta[0] = (th1a0+th1b0*std::exp(thetax[0]))/(1+std::exp(thetax[0]));
        Xbeta_r = X_r*beta_r;
      }
      nn = iscan+1;
      beth1Barold = beth1Barnew;
      beth1Barnew = (nn)/(nn+1.0)*beth1Barold + beth1/(nn+1.0);
      beth1Snew = (nn-1.0)/nn*beth1Snew + adapter/(p1+0.0)/nn*(nn*beth1Barold*beth1Barold.t() 
                                                                 - (nn+1.0)*beth1Barnew*beth1Barnew.t()
                                                                 + beth1*beth1.t() + Ip1 );
    }else if((th1a0>=th1b0)&(th2a0<th2b0)){
      beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llold);
      llold += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
      llold += thetax[1] - 2*std::log(1+std::exp(thetax[1]));
      beth1old = beth1;
      if(iscan>l0){
        beth1 = mvrnorm(beth1old, beth1Snew);
      }else{
        beth1 = mvrnorm(beth1old, SVhat1);
      }
      beta_r = beth1.subvec(0, p-1);
      thetax[1] = beth1[p];
      theta[1] = (th2a0+th2b0*std::exp(thetax[1]))/(1+std::exp(thetax[1]));
      Xbeta_r = X_r*beta_r;
      beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llnew);
      llnew += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
      llnew += thetax[1] - 2*std::log(1+std::exp(thetax[1]));
      ratio = exp(llnew-llold);
      uu = unif_rand();
      if((uu>ratio)){
        beth1=beth1old;
        if(iscan>=nburn) {
          rejbeta+=1.0;
          rejtheta+=1.0;
        }
        beta_r = beth1.subvec(0, p-1);
        thetax[1] = beth1[p];
        theta[1] = (th2a0+th2b0*std::exp(thetax[1]))/(1+std::exp(thetax[1]));
        Xbeta_r = X_r*beta_r;
      }
      nn = iscan+1;
      beth1Barold = beth1Barnew;
      beth1Barnew = (nn)/(nn+1.0)*beth1Barold + beth1/(nn+1.0);
      beth1Snew = (nn-1.0)/nn*beth1Snew + adapter/(p1+0.0)/nn*(nn*beth1Barold*beth1Barold.t() 
                                                                 - (nn+1.0)*beth1Barnew*beth1Barnew.t()
                                                                 + beth1*beth1.t() + Ip1 );
    }else {
      beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llold);
      llold += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
      betaold = beta_r;
      if(iscan>l0){
        beta_r = mvrnorm(betaold, beSnew);
      }else{
        beta_r = mvrnorm(betaold, Shat);
      }
      Xbeta_r = X_r*beta_r;
      beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llnew);
      llnew += -0.5*arma::dot( (beta_r-beta0), (S0inv*(beta_r-beta0)) );
      ratio = exp(llnew-llold);
      uu = unif_rand();
      if((uu>ratio)){
        beta_r=betaold;
        if(iscan>=nburn) rejbeta+=1.0;
        Xbeta_r = X_r*beta_r;
      }
      nn = iscan+1;
      beBarold = beBarnew;
      beBarnew = (nn)/(nn+1.0)*beBarold + beta_r/(nn+1.0);
      beSnew = (nn-1.0)/nn*beSnew + adapter/(p+0.0)/nn*(nn*beBarold*beBarold.t() 
                                                          - (nn+1.0)*beBarnew*beBarnew.t() 
                                                          + beta_r*beta_r.t() + Ip );
    }
    
    ///////////////////////////////////////////////
    // update phi
    //////////////////////////////////////////////
    beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llold);
    llold += phia0*logphi - phib0*phi;
    phiold = logphi;
    if(iscan>l0){
      logphi = Rf_rnorm(phiold, std::sqrt(phisnew));
    }else{
      logphi = Rf_rnorm(phiold, std::sqrt(phishat));
    }
    phi = std::exp(logphi);
    beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, llnew);
    llnew += phia0*logphi - phib0*phi;
    ratio = exp(llnew-llold);
    uu = unif_rand();
    if(uu>ratio){
      logphi = phiold;
      phi = std::exp(logphi);
      if(iscan>=nburn) rejphi+=1.0;
    }
    nn = iscan+1;
    phibarold = phibarnew;
    phibarnew = (nn)/(nn+1.0)*phibarold + logphi/(nn+1.0);
    phisnew = (nn-1.0)/nn*phisnew + adapter/nn*(nn*pow(phibarold,2) - (nn+1.0)*pow(phibarnew,2)
                                                  + pow(logphi,2) + ESMALL );
    
    ///////////////////////////////////////////////
    // Save the sample
    //////////////////////////////////////////////
    if(iscan>=nburn){
      ++skiptally;
      if(skiptally>nskip){
        // calculate loglikelihood
        logLik.col(isave) = beta4_mode_logliki(y, Xbeta, phi, theta[0], theta[1], link);
        // calculate -2loglikelihood
        double tempLik = 0;
        beta4_mode_loglik(y, Xbeta, phi, theta[0], theta[1], link, tempLik);
        Dvalues(isave) = -2.0*tempLik;
        // save regression coefficient
        beta_save(_,isave) = beta;
        sumbeta += beta;
        // save boundaries
        theta_save(_,isave) = theta;
        sumtheta += theta; 
        // precision parameter
        phi_save[isave] = phi;
        sumphi += phi;
        // save mode
        arma::vec Xpredbeta = Xpred*beta_r;
        for(int i=0; i<m; ++i){
          mode_save(i, isave) = ilinkf(Xpredbeta[i], link)*(theta[1]-theta[0])+theta[0];
        }
        
        ++isave;
        ++distally;
        if(distally>=ndisplay){
          Rprintf( "scan = %d\n", isave );
          distally = 0;
        }
        skiptally=0;
      }
    }
  }
  
  // get acceptance rate
  double totscans = nscan-nburn+0.0;
  double ratetheta = 1.0 - rejtheta/totscans;
  double ratebeta = 1.0 - rejbeta/totscans;
  double ratephi = 1.0 - rejphi/totscans;
  
  // get cpo
  arma::mat fpostinv = arma::exp(-logLik);
  arma::vec Linvmean = arma::mean(fpostinv, 1);
  arma::vec cpo = 1.0/Linvmean;
  
  // get DIC
  double meanD = arma::mean(Dvalues);
  Xbeta_r = X_r*as<arma::vec>(sumbeta)/(nsave+0.0);
  double Dmean = 0;
  beta4_mode_loglik(y, Xbeta, sumphi/(nsave+0.0), sumtheta[0]/(nsave+0.0), sumtheta[1]/(nsave+0.0), link, Dmean);
  double pD = meanD + 2.0*Dmean;
  double DIC = meanD + pD; 
  arma::vec DIC_pD(2); DIC_pD[0]=DIC; DIC_pD[1]=pD;
  
  //get stabilized cpo
  arma::mat fpost = arma::exp(logLik);
  arma::mat fpostinv_bd = arma::zeros<mat>(n, nsave);
  fpostinv_bd.each_col() = std::sqrt(nsave+0.0)*Linvmean;
  fpostinv_bd = arma::min(fpostinv, fpostinv_bd);
  arma::vec cpo_stab = arma::mean(fpostinv_bd%fpost, 1)/arma::mean(fpostinv_bd, 1);
  
  // get WAIC
  arma::vec fpostmean=arma::mean(fpost, 1);
  double lpd = arma::sum(arma::log(fpostmean));
  arma::vec logfpostvar=arma::var(logLik, 0, 1);
  double pwaic = arma::sum(logfpostvar);
  double WAIC = -2.0*lpd + 2.0*pwaic;
  arma::vec WAIC_pwaic(2); WAIC_pwaic[0]=WAIC; WAIC_pwaic[1]=pwaic;
  
  return List::create(Named("theta")=theta_save,
                      Named("beta")=beta_save,
                      Named("phi")=phi_save,
                      Named("fitted")=mode_save,
                      Named("cpo")=cpo,
                      Named("cpo_stab")=cpo_stab,
                      Named("DIC_pD")=DIC_pD,
                      Named("WAIC_pwaic")=WAIC_pwaic,
                      Named("ratetheta")=ratetheta,
                      Named("ratebeta")=ratebeta,
                      Named("ratephi")=ratephi);
  END_RCPP
}

// Get Cox-Snell residuals
RcppExport SEXP 
  beta4_mode_cox_snell(
    SEXP y_, SEXP X_, SEXP beta_, SEXP phi_, SEXP theta_, 
    SEXP link_, SEXP tgrid_, SEXP CI_
  ){
  BEGIN_RCPP
  // Transfer R variables into C++;
  const Rcpp::NumericVector y(y_); // n by 1;
  const arma::mat X=as<arma::mat>(X_); // n by p;
  const int n = y.size();
  const arma::mat beta=as<arma::mat>(beta_); // p by nsave;
  const arma::vec phi=as<arma::vec>(phi_); // nsave by 1;
  const arma::mat theta=as<arma::mat>(theta_); // 2 by nsave;
  const int link = as<int>(link_);
  Rcpp::NumericVector tgrid(tgrid_);
  double CI = as<double>(CI_);
  int ngrid = tgrid.size();
  int nsave = beta.n_cols;
  int low = nsave*(1.0-CI)*0.5 - 1;
  int up = nsave*(CI+(1.0-CI)*0.5) - 1;
  
  // things to save;
  arma::mat resid=arma::zeros<arma::mat>(n, nsave);
  arma::mat esth = arma::zeros<arma::mat>(ngrid, nsave);
  
  // get residuals
  for(int isave=0; isave<nsave; ++isave){
    double th1 = theta(0,isave);
    double th2 = theta(1,isave);
    double phi0 = phi[isave];
    double mi = 0;
    double St = 0;
    arma::vec xbeta = X*beta.col(isave);
    for(int i=0; i<n; ++i){
      mi = ilinkf(xbeta[i], link);
      resid(i, isave) = - std::log( Rf_pbeta((y[i]-th1)/(th2-th1), phi0*mi+1,
                                    phi0*(1-mi)+1, false, false) );
    }
    for(int j=0; j<ngrid; ++j){
      St = 0;
      for(int i=0; i<n; ++i){
        if(resid(i,isave)>tgrid[j]) ++St;
      }
      esth(j, isave) = -std::log( St/(n+0.0) );
    }
  }
  arma::vec hhat = arma::mean(esth, 1);
  arma::mat temp = arma::sort(esth,"ascend", 1);
  arma::vec hhatlow = temp.col(low);
  arma::vec hhatup = temp.col(up);
  
  return List::create(Named("resid")=resid,
                      Named("Hhat")=hhat,
                      Named("Hhatlow")=hhatlow,
                      Named("Hhatup")=hhatup,
                      Named("H")=esth);
  END_RCPP
}


