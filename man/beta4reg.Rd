\name{beta4reg}
\alias{beta4reg}
\alias{print.beta4reg}
\alias{summary.beta4reg}
\alias{print.summary.beta4reg}

\title{Bayesian Beta Regression Models}
\description{
This function fits Bayesian beta regression models. The response distribution can be either the beta with the support on (0,1)  or the four-parameter beta with an unknown final support. The logarithm of the pseudo marginal likelihood (LPML), the deviance information criterion (DIC), and the Watanabe-Akaike information criterion (WAIC) are provided for model comparison.
}
   
\usage{
beta4reg(formula, data, na.action, link="logit", model = "mode",
         mcmc=list(nburn=3000, nsave=2000, nskip=0, ndisplay=500),
         prior=NULL, start=NULL, Xpred=NULL)
}

\arguments{
   \item{formula}{a formula expression of the form \code{y ~ x}.}
   
   \item{data}{a data frame in which to interpret the variables named in the \code{formula} argument.}
   
   \item{na.action}{a missing-data filter function, applied to the \code{model.frame}.}
   
   \item{link}{a character string for the link function. Choices include \code{"logit"}, \code{"probit"}, \code{"loglog"} and \code{"cloglog"}.}
   
   \item{model}{a character string for the regression type. The options include \code{"mean"} for a mean regression, \code{"mode"} for a mode regression.}
   
   \item{mcmc}{a list giving the MCMC parameters. The list must include the following elements: \code{nburn} an integer giving the number of burn-in scans, \code{nskip} an integer giving the thinning interval, \code{nsave} an integer giving the total number of scans to be saved, \code{ndisplay} an integer giving the number of saved scans to be displayed on screen (the function reports on the screen when every \code{ndisplay} iterations have been carried out).}
   
   \item{prior}{a list giving the prior information. The function itself provides all default priors. The following components can be specified here: \code{ma0} and \code{mb0} for the prior of marginal population mode or mean, \code{phia0} and \code{phib0} for the precision parameter, \code{beta0} and \code{S0} for the coefficients beta, \code{th1a0} and \code{th1b0} for the lower bound of the support, \code{th2a0} and \code{th2b0} for the upper bound of the support.}
   
   \item{start}{a list giving the starting values of the parameters. The function itself provides all default choices. The following components can be specified here: \code{beta}, \code{theta}, \code{phi}. } 
   
   \item{Xpred}{A new design matrix at which estimates of the response model or mean are required. The default is the design matrix returned by the argument \code{formula}.}
}


\value{
This class of objects is returned by the \code{beta4reg} function to represent a fitted Bayesian beta regression model. Objects of this class have methods for the functions \code{print} and \code{summary}.

The \code{beta4reg} object is a list containing the following components:
\item{modelname}{
the name of the fitted model
}
\item{terms}{
the \code{\link{terms}} object used
}
\item{link}{
the link function used 
}
\item{model}{
the model fitted: mean or mode
}
\item{coefficients}{
a named vector of coefficients. The last two elements are the estimates of theta1 and theta2 involved in the support of the four-parameter beta distribution. 
}
\item{call}{
the matched call
}
\item{prior}{
the list of hyperparameters used in all priors. 
}
\item{start}{
the list of starting values used for all parameters. 
}
\item{mcmc}{
the list of MCMC parameters used 
}
\item{n}{
the number of row observations used in fitting the model
}
\item{p}{
the number of columns in the model matrix
}
\item{y}{
the response observations
}
\item{X}{
the n by (p+1) orginal design matrix
}
\item{beta}{
the (p+1) by nsave matrix of posterior samples for the coefficients in the \code{linear.predictors}
}
\item{theta}{
the 2 by nsave matrix of posterior samples for theta1 and theta2 involved in the support.
}
\item{phi}{
the vector of posterior samples for the precision parameter. 
}
\item{cpo}{
the length n vector of the stabilized estiamte of CPO; used for calculating LPML
}
\item{pD}{
the effective number of parameters involved in DIC
}
\item{DIC}{
the deviance information criterion (DIC)
}
\item{pW}{
the effective number of parameters involved in WAIC
}
\item{WAIC}{
the Watanabe-Akaike information criterion (WAIC)
}
\item{ratetheta}{
the acceptance rate in the posterior sampling of theta vector involved in the support
}
\item{ratebeta}{
the acceptance rate in the posterior sampling of beta coefficient vector
}
\item{ratephi}{
the acceptance rate in the posterior sampling of precision parameter
}

The use of the \code{summary} function to the object will return new object with the following additional components:
\item{coeff}{
A table that presents the posterior summaries for the regression coefficients
}
\item{bounds}{
A table that presents the posterior summaries for the support boundaries theta1 and theta2
}
\item{phivar}{
A table that presents the posterior summaries for the precision phi.
}

}

\seealso{
\code{\link{cox.snell.beta4reg}}
}

\references{
Zhou, H. and Huang, X. (2020). Bayesian beta regression for bounded responses with unknown supports. Submitted. 
}


\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Xianzheng Huang \email{<huang@stat.sc.edu>}
}

\examples{
rm(list=ls())
library(betaBayes)
library(coda)
library(betareg)

## Data from Ferrari and Cribari-Neto (2004)
data("FoodExpenditure", package = "betareg")

## four-parameter beta mean regression
mcmc=list(nburn=2000, nsave=1000, nskip=4, ndisplay=1000);
# Note larger nburn, nsave and nskip should be used in practice.
prior = list(th1a0 = 0, th2b0 = 1) 
# here the natural bound (0,1) is used to specify the prior
set.seed(100)
res1 <- beta4reg(I(food/income) ~ income + persons, data = FoodExpenditure, 
                link = "logit", model = "mean",
                mcmc = mcmc, prior = prior)
(sfit1 <- summary(res1))
cox.snell.beta4reg(res1) # Cox-Snell plot

## two-parameter beta mean regression with support (0,1)
mcmc=list(nburn=2000, nsave=1000, nskip=4, ndisplay=1000);
# Note larger nburn, nsave and nskip should be used in practice.
prior = list(th1a0 = 0, th1b0 = 0, th2a0 = 1, th2b0 = 1)
# this setting forces the support to be (0,1)
set.seed(100)
res2 <- beta4reg(I(food/income) ~ income + persons, data = FoodExpenditure, 
                link = "logit", model = "mean",
                mcmc = mcmc, prior = prior)
(sfit2 <- summary(res2))
cox.snell.beta4reg(res2) # Cox-Snell plot

}

