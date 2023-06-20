## This file created by Hassan Pazira at 16-12-2022
## Updated at 09-09-2023

negloglik.beta.sigma <- function(beta_sig, y, X, Lambda, family){
  if (is.data.frame(X)) X <- data.matrix(X) # if X is data.frame
  nl <- nrow(X)
  p <- ncol(X)
  beta <- beta_sig[1:p]
  v <- beta_sig[length(beta_sig)] # v=log(sigma2) or sigma2=exp(2*v)
  power <- c(0)
  if (family == "binomial") {
    for (j in 1:nl) {
      power <- power + y[j] * as.numeric(X[j,] %*% beta) - log(1 + exp(as.numeric(X[j,] %*% beta)))
    }
    negloglik  <- - power + (t(beta) %*% Lambda %*% beta)/2
  }
  if (family == "gaussian") {
    for (j in 1:nl) {
      power <- power + (y[j] - as.numeric(X[j,] %*% beta))^2
    }
    negloglik  <- power/exp(2*v) + nl * log(exp(2*v)) + exp(2*v) * as.numeric(Lambda[(p+1), (p+1)]) +
      t(beta) %*% Lambda[-(p+1),-(p+1)] %*% beta
  }
  return(as.numeric(negloglik))
}
