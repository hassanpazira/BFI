## This file created by Hassan Pazira at 16-12-2022
negloglik_beta_sigma <- function(beta_sig, y, X, Gamma, family){
  nl <- nrow(X)
  p <- ncol(X)
  beta <- beta_sig[1:p]
  sigma2_e <- beta_sig[length(beta_sig)]
  power <- c(0)
  if (family == "binomial") {
    for (j in 1:nl) {
      power <- power + y[j] * as.numeric(X[j,] %*% beta) - log(1 + exp(as.numeric(X[j,] %*% beta)))
    }
    negloglik  <- - power + (t(beta) %*% Gamma %*% beta)/2
  }
  if (family == "gaussian") {
    for (j in 1:nl) {
      power <- power + (y[j] - as.numeric(X[j,] %*% beta))^2
    }
    negloglik  <- power/sigma2_e + nl * log(sigma2_e) + sigma2_e * as.numeric(Gamma[(p+1), (p+1)]) +
      t(beta) %*% Gamma[-(p+1),-(p+1)] %*% beta
  }
  return(as.numeric(negloglik))
}
