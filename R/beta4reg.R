"beta4reg" <- function(formula, data, na.action, link="logit", model = "mode", 
                       mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500), 
                       prior=NULL, start=NULL, Xpred=NULL){
  #########################################################################################
  # call parameters
  #########################################################################################
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf);
  linkcode = switch(link, logit=1, probit=2, 
                    cloglog=3, loglog=4);
  if(linkcode==1){
    linkf = .logit
    ilinkf = .ilogit
    thvar0 = pi^2/3
  }else if(linkcode==2){
    linkf = .probit
    ilinkf = .iprobit
    thvar0 = 1
  }else if(linkcode==3){
    linkf = .cloglog
    ilinkf = .icloglog
    thvar0 = 1.6
  }else{
    linkf = .loglog
    ilinkf = .iloglog
    thvar0 = 1.6
  }
  
  n <- length(y)
  p <- ncol(X)
  
  #########################################################################################
  # initial MLE analysis and mcmc parameters
  #########################################################################################
  yu = unique(y)
  yuorder = order(yu)
  y.min = yu[yuorder[1]]
  y.max = yu[yuorder[length(yu)]]
  
  #########################################################################################
  # priors
  #########################################################################################
  nburn <- mcmc$nburn;
  nsave <- mcmc$nsave;
  nskip <- mcmc$nskip;
  ndisplay <- mcmc$ndisplay;
  ma0=prior$ma0; if(is.null(ma0)) ma0=1;
  mb0=prior$mb0; if(is.null(mb0)) mb0=1;
  phia0=prior$phia0; if(is.null(phia0)) phia0=.001;
  phib0=prior$phib0; if(is.null(phib0)) phib0=.001;
  beta0=prior$beta0; 
  if(linkcode==1){
    bb0 = digamma(ma0) - digamma(mb0)
    vv0 = trigamma(ma0) + trigamma(mb0)
  }else{
    bb0 = integrate(function(x) linkf(x)*dbeta(x, ma0, mb0), 
                    lower = 0, upper = 1, rel.tol = 1e-6)$value
    vv0 = integrate(function(x) (linkf(x)-bb0)^2*dbeta(x, ma0, mb0), 
                    lower = 0, upper = 1, rel.tol = 1e-6)$value
  }
  if(is.null(beta0)){
    beta0 = rep(0, p)
    beta0[1] = bb0
  } 
  S0 <- prior$S0; 
  if(is.null(S0)){
    g = vv0/p
    S0= g*n*solve(t(X)%*%X)
  } 
  S0inv <- solve(S0)
  th1a0 <- prior$th1a0; if(is.null(th1a0)) th1a0 = y.min - 2*sd(y)
  th1b0 <- prior$th1b0; if(is.null(th1b0)) th1b0 = y.min - .RSMALL
  th2a0 <- prior$th2a0; if(is.null(th2a0)) th2a0 = y.max + .RSMALL 
  th2b0 <- prior$th2b0; if(is.null(th2b0)) th2b0 = y.max + 2*sd(y)
  if(th1b0>=y.min) stop("th1b0 needs to be strictly smaller than min(y)")
  if(th2a0<=y.max) stop("th2a0 needs to be strictly greater than max(y)")
  
  if(is.null(Xpred)) Xpred = X
  
  prior <- list(beta0 = beta0, S0 = S0, ma0=ma0, mb0=mb0, 
                phia0 = phia0, phib0 = phib0, th1a0 = th1a0, th1b0 = th1b0,
                th2a0 = th2a0, th2b0 = th2b0)
  ############ initials #########
  # theta
  delta.th = sd(y)/sqrt(n)
  theta = start$theta; 
  theta.i = start$theta; 
  if(is.null(theta)){
    theta = c(ifelse((th1b0-th1a0)>(2*delta.th), 
                           th1b0-delta.th, (th1a0+th1b0)/2), 
                    ifelse((th2b0-th2a0)>(2*delta.th), 
                           th2a0+delta.th, (th2a0+th2b0)/2) )
    theta.i= c(ifelse((th1b0-th1a0)>(2*delta.th), 
                             th1b0-delta.th, (th1a0+th1b0)/2), 
                      ifelse((th2b0-th2a0)>(2*delta.th), 
                             th2a0+delta.th, (th2a0+th2b0)/2) )
    
  }else{
    if((theta[1]<th1a0)|(theta[1]>th1b0)|(theta[2]<th2a0)|(theta[2]>th2b0)){
      stop("start$theta needs to fall within [th1a0, th1b0] and [th2a0, th2b0]")
    }
  }
  th1.var = (th1b0-th1a0)^2/12
  th2.var = (th2b0-th2a0)^2/12
  th.var = var(y)/n/12
  Vhat0 = diag(rep(thvar0, 2))
  # beta
  ## initial betareg fit
  a0 <- theta.i[1] + 0
  b0 <- theta.i[2] + 0
  fit0 <- betareg::betareg(I((y-a0)/(b0-a0)) ~ X - 1)
  sfit0 <- summary(fit0)
  beta = start$beta; if(is.null(beta)) beta = fit0$coefficients$mean+0
  beta.i = start$beta; if(is.null(beta.i)) beta.i = fit0$coefficients$mean+0
  Shat0 = sfit0$vcov[1:p, 1:p]
  #Shat0 = diag(rep(1/n, p))
  
  # variability
  phi = start$phi; if(is.null(phi)) phi = fit0$coefficients$precision
  phi.i = start$phi; if(is.null(phi.i)) phi.i = fit0$coefficients$precision
  start0 <- list(theta = theta.i, beta = beta.i, phi = phi.i)
  
  #########################################################################################
  # calling the c++ code and # output
  #########################################################################################
  if(model == "mode"){
    model.name <- "Beta mode regression with unknown boundaries:";
    foo <- .Call("beta4_mode_reg", nburn_=nburn, nsave_=nsave, nskip_=nskip, 
                 ndisplay_=ndisplay, y_=y, X_=X, beta_=beta, phi_=phi, 
                 theta_=theta, phia0_=phia0, phib0_=phib0,
                 th1a0_=th1a0, th1b0_=th1b0, th2a0_=th2a0, th2b0_=th2b0,
                 Vhat_=Vhat0, beta0_=beta0, S0inv_=S0inv, Shat_=Shat0, 
                 l0_=round(min(1000,nburn/2)), adapter_=2.38^2, 
                 link_=linkcode, Xpred_=Xpred, PACKAGE = "betaBayes");
  }else if(model == "mean"){
    model.name <- "Beta mean regression with unknown boundaries:";
    foo <- .Call("beta4_mean_reg", nburn_=nburn, nsave_=nsave, nskip_=nskip, 
                 ndisplay_=ndisplay, y_=y, X_=X, beta_=beta, phi_=phi, 
                 theta_=theta, phia0_=phia0, phib0_=phib0,
                 th1a0_=th1a0, th1b0_=th1b0, th2a0_=th2a0, th2b0_=th2b0,
                 Vhat_=Vhat0, beta0_=beta0, S0inv_=S0inv, Shat_=Shat0, 
                 l0_=round(min(1000,nburn/2)), adapter_=2.38^2, 
                 link_=linkcode, Xpred_=Xpred, PACKAGE = "betaBayes");
  }else{
    stop("it only supports mode or mean regressions")
  }
  
  #########################################################################################
  # save state
  #########################################################################################
  #### coefficients
  coeff1 <- c(apply(matrix(foo$beta, p, nsave), 1, mean));
  coeff2 <- c(apply(foo$theta, 1, mean));
  coeff <- c(coeff1, coeff2);
  names(coeff) = c(colnames(X),"theta1", "theta2");
  #### Save to a list
  output <- list(modelname=model.name,
                 model = model,
                 terms = mt,
                 link = link,
                 coefficients=coeff,
                 call=cl,
                 mcmc=mcmc,
                 n=n,
                 p=p,
                 y=y,
                 X = X,
                 beta = matrix(foo$beta, p, nsave),
                 theta = foo$theta,
                 phi = foo$phi,
                 yhat = foo$fitted,
                 prior = prior,
                 start = start0,
                 cpo = foo$cpo_stab,
                 pD = foo$DIC_pD[2], 
                 DIC = foo$DIC_pD[1],
                 pW = foo$WAIC_pwaic[2],
                 WAIC = foo$WAIC_pwaic[1],
                 ratetheta = foo$ratetheta,
                 ratebeta = foo$ratebeta,
                 ratephi = foo$ratephi);
  class(output) <- c("beta4reg")
  output
}

#### print, summary, plot
"print.beta4reg" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior means for regression coefficients and boundaries:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat("\nLPML:", sum(log(x$cpo)))
  cat("\nDIC:", x$DIC)
  cat("\nWAIC:", x$WAIC)
  cat("\nn=",x$n, "\n", sep="")
  invisible(x)
}

"summary.beta4reg" <- function(object, CI.level=0.95, ...) {
  ans <- c(object[c("call", "modelname")])
  
  ### CPO
  ans$cpo <- object$cpo
  
  ### Median information
  mat <- as.matrix(object$beta, nrow=object$p)
  coef.p <- object$coefficients[(1:object$p)];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(quantile(x, probs=c((1-CI.level)/2, 1-(1-CI.level)/2))) )
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$coeff <- coef.table
  
  ### Boundary Information
  mat <- as.matrix(object$theta)
  coef.p <- object$coefficients[object$p+(1:2)];
  coef.m <- apply(mat, 1, median)    
  coef.sd <- apply(mat, 1, sd)
  limm <- apply(mat, 1, function(x) as.vector(quantile(x, probs=c((1-CI.level)/2, 1-(1-CI.level)/2))) )
  coef.l <- limm[1,]
  coef.u <- limm[2,]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$bounds <- coef.table
  
  ### Precision parameter
  mat <- object$phi
  coef.p <- mean(mat); names(coef.p)="phi";    
  coef.m <- median(mat)    
  coef.sd <- sd(mat)
  limm <- as.vector(quantile(mat, probs=c((1-CI.level)/2, 1-(1-CI.level)/2)))
  coef.l <- limm[1]
  coef.u <- limm[2]
  
  coef.table <- cbind(coef.p, coef.m, coef.sd, coef.l , coef.u)
  dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", 
                                                paste(CI.level*100, "%CI-Low", sep=""),
                                                paste(CI.level*100, "%CI-Upp", sep="")))
  ans$phivar <- coef.table

  ## LPML and DIC
  ans$n <- object$n
  ans$p <- object$p
  ans$LPML <- sum(log(object$cpo))
  ans$DIC <- object$DIC
  ans$WAIC <- object$WAIC
  
  ### acceptance rates
  ans$ratetheta = object$ratetheta;
  ans$ratebeta = object$ratebeta;
  ans$ratephi = object$ratephi;
  
  class(ans) <- "summary.beta4reg"
  return(ans)
}

"print.summary.beta4reg"<-function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat(x$modelname,"\nCall:\n", sep = "")
  print(x$call)
  
  cat("\nPosterior inference of regression coefficients\n")
  cat("(Adaptive M-H acceptance rate: ", x$ratebeta, "):\n", sep="")
  print.default(format(x$coeff, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat("\nPosterior inference of boundary parameters\n")
  cat("(Adaptive M-H acceptance rate: ", x$ratetheta, "):\n", sep="")
  print.default(format(x$bounds, digits = digits), print.gap = 2, quote = FALSE)
  
  cat("\nPosterior inference of precision parameter phi\n")
  cat("(Adaptive M-H acceptance rate: ", x$ratephi, "):\n", sep="")
  print.default(format(x$phivar, digits = digits), print.gap = 2, 
                quote = FALSE)
  
  cat("\nLog pseudo marginal likelihood: LPML=", x$LPML, sep="")
  cat("\nDeviance Information Criterion: DIC=", x$DIC, sep="")
  cat("\nWatanabe-Akaike information criterion: WAIC=", x$WAIC, sep="")
  cat("\nNumber of subjects: n=", x$n, "\n", sep="")
  
  invisible(x)
}

"cox.snell.beta4reg" <- function (x, ncurves = 10, CI = 0.95, PLOT = TRUE) {
  res <- list()
  if(is(x,"beta4reg")){
    linkcode = switch(x$link, logit=1, probit=2, 
                      cloglog=3, loglog=4);
    tgrid = seq(0, 7, 0.05)
    if(x$model=="mode"){
      foo <- .Call("beta4_mode_cox_snell", y_=x$y, X_=x$X,
                   beta_=x$beta, phi_=x$phi, theta_=x$theta,
                   link_=linkcode, tgrid_=tgrid, CI_=CI, 
                   PACKAGE = "betaBayes");
    }else if(x$model=="mean"){
      foo <- .Call("beta4_mean_cox_snell", y_=x$y, X_=x$X,
                   beta_=x$beta, phi_=x$phi, theta_=x$theta,
                   link_=linkcode, tgrid_=tgrid, CI_=CI, 
                   PACKAGE = "betaBayes");
    }else{
      stop("This function only supports beta mean and mode models");
    }
    res$tgrid = tgrid
    #res$resid=foo$resid;
    res$Hhat=foo$Hhat
    res$Hhatlow = foo$Hhatlow
    res$Hhatup = foo$Hhatup
    res$H = foo$H
    #res$coverage = mean((foo$hhatlow<tgrid)&(foo$hhatup>tgrid)) 
    
    if(PLOT){
      r.max <- max(tgrid)
      xlim <- c(0, r.max); ylim <- c(0, r.max)
      xx <- seq(0, r.max, 0.01)
      par(cex = 1.5, mar = c(2.1,2.1,1,1), cex.lab = 1.4, cex.axis = 1.1)
      plot(tgrid, foo$hhat, xlim = xlim,
           ylim = ylim, lwd = 3, "l", col="white")
      lines(xx, xx, lty = 1, lwd = 2, col = "darkgrey")
      #lines(tgrid, foo$Hhat, lty = 1, lwd = 3)
      lines(tgrid, foo$Hhatlow, lty = 2, lwd = 3, col=2)
      lines(tgrid, foo$Hhatup, lty = 2, lwd = 3, col=2)
      indx = sample(x$mcmc$nsave, ncurves)
      for(i in 1:ncurves){
        lines(tgrid, foo$H[,indx[i]], lwd=2, lty=3)
      }
    }
  }
  res
  invisible(res)
}

