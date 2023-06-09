## This file created by Hassan Pazira at 16-12-2022

## Required Packages:
if(!require(pracma)){install.packages("pracma")} # for the NR methods and the condition number of a matrix 
library(pracma)


#######################
## Required function(s)
#----------------------

Estimators_maker <- function(response, X, family=c("binomial","gaussian"), Gamma){
  # Gamma is the 'inverse' covariance matrix (could be matrix or list of matrices)!
  y <- response
  family <- match.arg(family)
  if (!family %in% c("binomial", "gaussian")) {
    stop("Distributions that can be used are 'binomial' and 'gaussian' in this version of the package!")
  }
  if (is.data.frame(X)) X <- as.matrix(X)
  n <- nrow(X)
  if (any(as.numeric(X[,1])!=rep(1,n))) {
    if (is.null(colnames(X))) {
      colnames(X) <- paste0("X",1:ncol(X))
    }
    X <- cbind(rep(1, n), X)
    colnames(X)[1] <- c("Intercept")
  } else {
    if (is.null(colnames(X))) {
      colnames(X) <- c("Intercept", paste0("X",1:ncol(X)))
    }
  }
  p <- ncol(X)   # number of predictors with intercept!
  if (dim(Gamma)[1]!=dim(Gamma)[2]) stop("'Gamma' should be a square matrix" )
  if (family=="binomial") {
    if (dim(Gamma)[1] != p) stop("'Gamma' in this family should have the 'p' dimension: 
                                 where 'p' is number of coefficients with intercept!" )
  }
  if (family=="gaussian") {
    if (dim(Gamma)[1] != (p+1)) stop(paste("'Gamma' should have the 'p+1' dimension:",
                                           "where 'p' is number of coefficients with intercept,",
                                           "and '1' is for error variance!"))
  }
  if (family=="binomial") {
    if (NCOL(y) == 1) {
      if (length(y)!=n) stop("length of 'y' != sample size 'n'")
      if (is.character(y)) y <- factor(y)
      if (is.factor(y)){
        if (nlevels(y)!=2) stop("only factors with two levels are allowed")
        yy <- as.numeric(y) - 1
      } else {
        if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
        if (any(abs(y - round(y)) > 0.001)) warning("non-integer #successes in a binomial!") #
        yy <- y
      }
    } else if (NCOL(y) == 2) {
      if (dim(y)[1]!=n) stop("sample size != 'n'")
      if(any(abs(y - round(y)) > 0.001)) warning("non-integer counts in a binomial!")
      nn <- y[, 1] + y[, 2]
      y <- ifelse(nn == 0, 0, y[, 1]/nn)
      yy <- y
    } else stop(paste("For the 'binomial' family, y must be",
                      "a vector of 0 and 1\'s or a 2 column",
                      "matrix where col 1 is no. successes",
                      "and col 2 is no. failures"))
  }
  if (family=="gaussian") {
    if (NCOL(y) != 1) stop(paste0("For this family, '",family,"', y must be a vector, 
                                  not a matrix!"))
    if (length(y)!=n) stop("length of 'y' != sample size 'n'")
    if (is.character(y)) stop(paste0("It's not allowed to use a character or a factor as 
                                     a response vector for the ",family," family (-Inf < y < Inf)!"))
    if (any(is.na(y))) stop(paste0("NA is not allowed for the ",family," family!"))
    yy <- y
  }
  y <- yy
  
  ## Optimizations for beta's (and sigma2_e), where theta_hat: b0, b1, ..., bp (and sigma2_e)
  if (family == "binomial") {
    initial_beta_sig  <- c(rep(0, p))
    theta_hat <- try(optim(initial_beta_sig, fn=negloglik_beta_sigma, gr=NULL, y=y, X=X, Gamma=Gamma, 
                           family=family, lower=rep(-Inf,p), upper=rep(Inf,p), method="L-BFGS-B")$par, TRUE)
    if(class(theta_hat)[1] == "try-error") {
      warning("try-error in theta_hat !!!")
    }
    names(theta_hat) <- colnames(X)
  }
  if (family == "gaussian") {
    initial_beta_sig  <- c(rep(0, p), 0.5)
    theta_hat <- try(optim(initial_beta_sig, fn=negloglik_beta_sigma, gr=NULL, y=y, X=X, Gamma=Gamma, 
                           family=family, lower=rep(-Inf,p), upper=rep(Inf,p), method="L-BFGS-B")$par, TRUE)
    if(class(theta_hat)[1] == "try-error") {
      warning("try-error in theta_hat !!!")
    }
    names(theta_hat) <- c(colnames(X), "sigma2_e")
  }
  
  ## A_hat: curvature matrix estimator
  A_hat  <- A_l_maker(y=y, X=X, Gamma=Gamma, family=family, theta_hat=theta_hat)
  ## sd of A_hat
  sd_A <- sqrt(diag(ginv(as.matrix(A_hat))))
  
  output <- list(theta_hat=theta_hat, A_hat=A_hat, sd=sd_A, Gamma=Gamma)
  return(output)
}

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

BFI <- function(theta_hats=NULL, A_hats, Gamma, L=NULL) {
  if (!is.null(theta_hats) & !is.list(theta_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.list(A_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.null(theta_hats) & (length(theta_hats) != length(A_hats))) {
    stop("Length of inputs are not the same.")
  }
  if (is.null(L)) {
    L <- length(A_hats)
  } else {
    if (length(theta_hats) != L) {
      stop("Number of locations and length of inputs are not the same.")
    }
  }
  if (is.matrix(Gamma) | (is.list(Gamma) & length(Gamma)==1)) { # all locations have the same Gamma.
    A_fed <- Reduce("+", A_hats) + (1-L) * Gamma
  } else {
    if (!is.list(Gamma) | (is.list(Gamma) & length(Gamma)!=(L+1))) {
      stop("Gamma should contains of L+1 lists; 1:L for local datastes and last one for combined.")
    }
    A_fed <- Reduce("+", A_hats) + Gamma[[L+1]] - Reduce("+", Gamma[1:L]) 
  }
  sd_fed <- sqrt(diag(ginv(as.matrix(A_fed))))
  if (!is.null(theta_hats)) {
    theta_hat_fed <- ginv(as.matrix(A_fed)) %*% Reduce("+", Map("%*%", A_hats , theta_hats))
  } else theta_hat_fed <- theta_hats
  output <- list(theta_hat=theta_hat_fed, A_hat=A_fed, sd=sd_fed)
  return(output)
}

