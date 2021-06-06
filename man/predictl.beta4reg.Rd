\name{predict.beta4reg}
\alias{predict.beta4reg}

\title{Predict method for beta4 model fits}
\description{
Posterior predicted response values based on beta4 model object
}
   
\usage{
\method{predict}{beta4reg}(object, newx, ...)
}

\arguments{
   \item{object}{an object obtained from the function \code{\link{beta4reg}}.}
   
   \item{newx}{an m by p matrix at which predictions are required. If not specified, the original design matrix will be used.}
   
   \item{...}{further arguments passed to or from other methods.}
}


\value{
The function returns an m by nsave matrix of posterior samples for response predictions at newx. 
}

\seealso{
\code{\link{beta4reg}} 
}

\author{
Haiming Zhou \email{<zhouh@niu.edu>} and Xianzheng Huang \email{<huang@stat.sc.edu>}
}


