## This file created by Hassan Pazira at 16-12-2022
A_l_maker <- function(y, X, Gamma, family, theta_hat){
  nl <- nrow(X)
  p <- ncol(X)
  if (family == "binomial") {
    len_thet <- length(theta_hat)
    beta_hat <- theta_hat
    power <- matrix(0, len_thet, len_thet)
    for (j in 1:nl) {
      power <- power +  t(X[j,,drop=F]) %*% X[j,,drop=F] *
        exp(as.numeric(X[j,] %*% beta_hat))/(1+exp(as.numeric(X[j,] %*% beta_hat)))^2
    }
    A_l <- power + Gamma
    colnames(A_l) <- colnames(X)
  }
  if (family == "gaussian") {
    len_thet <- length(theta_hat)
    beta_hat <- theta_hat[-len_thet]
    sigma2_e_hat <- theta_hat[len_thet]
    power <- matrix(0, len_thet-1, len_thet-1)
    power2 <- c(0)
    power3 <- c(rep(0, p))
    for (j in 1:nl) {
      power <- power +  t(X[j,,drop=F]) %*% X[j,,drop=F]
      power2 <- power2 + (y[j] - as.numeric(X[j,] %*% beta_hat))^2
      power3 <- power3 + (y[j] - as.numeric(X[j,] %*% beta_hat)) * t(X[j,,drop=F])
    }
    A_l <- matrix(0, len_thet, len_thet)
    A_l[-len_thet, -len_thet] <- power/sigma2_e_hat + Gamma[-(p+1),-(p+1)]
    A_l[len_thet, len_thet] <- 2 * Gamma[len_thet,len_thet] * sigma2_e_hat + 2 * power2 / sigma2_e_hat
    A_l[1:(len_thet-1), len_thet] <- A_l[len_thet, 1:(len_thet-1)] <- -2 * power3 / sigma2_e_hat
    colnames(A_l) <- rownames(A_l) <- c(colnames(X), "sigma2_e")
  }
  return(A_l)
}
