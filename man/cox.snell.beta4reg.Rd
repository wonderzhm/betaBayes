\name{cox.snell.beta4reg}
\alias{cox.snell.beta4reg}

\title{Cox-Snell Diagnostic Plot}
\description{
This function provides the Cox-Snell diagnostic plot for fitting  for Bayesian beta regression models.
}
   
\usage{
cox.snell.beta4reg(x, ncurves = 10, CI = 0.95, PLOT = TRUE)
}

\arguments{
   \item{x}{an object obtained from the function \code{\link{beta4reg}}.}
   
   \item{ncurves}{the number of posterior draws.}
   
   \item{CI}{the level of confidence for point-wise credible intervals.}
   
   \item{PLOT}{a logical value indicating whether the Cox-Snell residuals will be plotted.}
}


\value{
The function returns the plot (if \code{PLOT = TRUE}) and a list with the following components:
\item{tgrid}{
the x-axis values with length, say \code{ngrid}
}
\item{Hhat}{
the \code{ngrid} by 1 averaged cumulative hazard values across the \code{nsave} posterior samples
}
\item{Hhatlow}{
the \code{ngrid} by 1 lower bound cumulative hazard values 
}
\item{Hhatup}{
the \code{ngrid} by 1 upper bound cumulative hazard values 
}
\item{H}{
the \code{ngrid} by \code{nsave} cumulative hazard values 
}
}

\seealso{
\code{\link{beta4reg}} 
}

\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Xianzheng Huang \email{<huang@stat.sc.edu>}
}


