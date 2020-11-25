# betaBayes
Bayesian Beta Regression
This R package provides a class of Bayesian beta regression modes for the analysis of continuous data with support restricted to an unknown finite support. The response variable is modeled using a four-parameter beta distribution with the mean or mode parameter depending linearly on covariates through a link function. When the response support is known to be (0,1), the above class of models reduce to traditional (0,1) supported beta regression models. Model choice is carried out via the logarithm of the pseudo marginal likelihood (LPML), the deviance information criterion (DIC), and the Watanabe-Akaike information criterion (WAIC).
To install this R package, type:
library(devtools)
install_github("wonderzhm/betaBayes")
When you are asked "These packages have more recent versions available. It is recommended to update all of them. Which would you like to update?", choose "None". 
