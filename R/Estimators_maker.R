## This file created by Hassan Pazira at 16-12-2022
estimators_maker <- function(y, X, family=c("binomial","gaussian"), Gamma){
  # Gamma is the 'inverse' covariance matrix (could be matrix or list of matrices)!
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
    initial_beta_sig  <- rep(0, p+1)
    theta_hat <- try(optim(initial_beta_sig, fn=negloglik_beta_sigma, gr=NULL, y=y, X=X, Gamma=Gamma,
                           family=family, lower=rep(-Inf,p), upper=rep(Inf,p), method="L-BFGS-B")$par, TRUE)
    if(class(theta_hat)[1] == "try-error") {
      warning("try-error in theta_hat !!!")
    }
    # Now, theta_hat returns sigma_e
    theta_hat[length(theta_hat)] <- exp(2 * theta_hat[length(theta_hat)]) # Now, theta_hat returns sigma2_e
    names(theta_hat) <- c(colnames(X), "sigma2_e")
  }

  ## A_hat: curvature matrix estimator
  A_hat  <- A_l_maker(y=y, X=X, Gamma=Gamma, family=family, theta_hat=theta_hat)
  ## sd of A_hat
  sd_A <- sqrt(diag(solve(as.matrix(A_hat))))

  output <- list(theta_hat=theta_hat, A_hat=A_hat, sd=sd_A, Gamma=Gamma)
  return(output)
}
