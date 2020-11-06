##
.RSMALL = 1e-15
## logit link functions
.ilogit <- function(x){
  1.0/(1.0+exp(-x))
}
.logit <- function(y){
  log(y/(1.0-y))
}

## probit link functions
.iprobit <- function(x){
  pnorm(x)
}
.probit <- function(y){
  qnorm(y)
}

## complementary log log link functions
.icloglog <- function(x){
  1-exp(-exp(x))
}
.cloglog <- function(y){
  log(-log(1-y))
}

## loglog link functions
.iloglog <- function(x){
  exp(-exp(-x))
}
.loglog <- function(y){
  -log(-log(y))
}
