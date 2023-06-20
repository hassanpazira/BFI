## This file created by Hassan Pazira at 16-12-2022
## Updated at 14-09-2023

MAP.estimation <- function(y, X, family=gaussian, Lambda, intercept=TRUE, initial=NULL,
                           control=list()) {
  # Lambda is the 'inverse' covariance matrix (could be matrix or list of matrices)!
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- do.call(family, args = list(), envir = parent.frame())
  if (is.null(family$family)) stop("'family' not recognized")
  if (!family$family %in% c("binomial", "gaussian")) {
    stop("Distributions that can be used are 'binomial' and 'gaussian' in this version of the package!")
  }
  #if (!missing(intercept) && family == "cox") warning("Cox model has no intercept")
  if (is.matrix(X))
    warning("Note that, if in 'X' there is a 'categorical' covariate with more than 2 levels, \n",
                            "then 'X' must be a 'data.frame' instead of a 'matrix'. ")
  n <- NROW(X)
  if (all(as.numeric(X[,1])==rep(1,n))) {
    X <- X[,-1, drop=F]
  }
  if (any(c("Intercept", "(Intercept)") %in% colnames(X)))
    stop("'intercept' should be the first column of 'X'.")
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X",1:ncol(X))
  } else {
    if (any(c("Intercept", "(Intercept)") %in% colnames(X)))
      stop("'intercept' should be the first column of 'X'.")
  }
  design_matrix <- paste(colnames(X), collapse=" + ")
  formula <- as.formula(paste("y", design_matrix, sep=" ~ "))
  y_X_vars <- model.maker(formula, as.data.frame(X))
  #formula1 <- noquote(paste("y", paste(colnames(y_X_vars$X)[-1], collapse=" + "), sep=" ~ "))
  formula1 <- noquote(paste("y", design_matrix, sep=" ~ "))
  X <- y_X_vars$X
  y <- y_X_vars$y
  if (NROW(y) != NROW(X)) stop("NROW(y) != NROW(X)")
  if (intercept == F) {
    X <- X[,-1, drop=F]
    if (any(c("Intercept", "(Intercept)") %in% colnames(X)))
      stop("'intercept' should be the first column of 'X'.")
  }
  if (any(diag(Lambda) <= 0) & all(Lambda[row(Lambda)!=col(Lambda)] == 0)) # second arg. is not needed
    stop("Diagonal elements of Lambda matrix should be positive (> 0)")
  if (intercept==T) {
    xvar <- apply(X[,-1, drop=F], 2, sd)
  } else {
    xvar <- apply(X, 2, sd)
  }
  xvar[xvar < 10 * .Machine$double.eps] <- 0
  const_vars <- sqrt(xvar) == 0
  if (any(const_vars)) {
    if (intercept==T) {
      stop("Predictor(s) ", paste0(colnames(X[,-1, drop=F])[which(const_vars)],sep=" "),
           "has zero variance! The 'bfi()' function can handle this variable.")
    } else {
      stop("Predictor(s) ", paste0(colnames(X)[which(const_vars)],sep=" "),
           "has zero variance! The 'bfi()' function can handle this variable.")
    }
  }
  p <- ncol(X)         # number of regression parameters/coefficients, i.e., predictors
                       # with intercept if intercept=T, or only predictors if intercept=F.
  if (dim(Lambda)[1]!=dim(Lambda)[2])
    stop("'Lambda' should be a square matrix." )
  if (family$family=="binomial") {
    if (dim(Lambda)[1] != p)
      stop(paste("'Lambda' in this family should have p =",sQuote(p),"dimensions where",
                 "'p' is number of coefficients with","\n",
                 "intercept (if intercept=T) or without intercept (if intercept=F))." ))
    if (NCOL(y) == 1) {
      if (length(y)!=n) stop("length of 'y' != sample size 'n'.")
      if (is.character(y)) y <- factor(y)
      if (is.factor(y)){
        if (nlevels(y)!=2) stop("only factors with two levels are allowed.")
        yy <- as.numeric(y) - 1
      } else {
        if (any(y < 0 | y > 1)) stop("values of 'y' must be 0 <= y <= 1")
        if (any(abs(y - round(y)) > 0.001)) warning("non-integer #successes in a binomial.") #
        yy <- y
      }
    } else if (NCOL(y) == 2) {
      if (dim(y)[1]!=n) stop("sample size != 'n'.")
      if(any(abs(y - round(y)) > 0.001)) warning("non-integer counts in a binomial.")
      nn <- y[, 1] + y[, 2]
      y <- ifelse(nn == 0, 0, y[, 1]/nn)
      yy <- y
    } else stop(paste("For the 'binomial' family, 'y' must be",
                      "a vector of 0 and 1\'s or a 2 column matrix \n",
                      "where col 1 is number successes",
                      "and col 2 is number failures"))
  }
  if (family$family=="gaussian") {
    if (dim(Lambda)[1] != (p+1))
      stop(paste("'Lambda' should have p+1 =",sQuote(p+1),"dimensions where",
                 " p =",sQuote(p),"is number of coefficients with", "\n",
                 " intercept (if intercept=T) or without intercept (if intercept=F), and '1' is for error variance."))
    if (NCOL(y) != 1)
      stop(paste0("For the ",sQuote(family$family)," family, 'y' must be a vector, not a matrix."))
    if (length(y)!=n)
      stop("length of 'y' != sample size 'n'")
    if (is.character(y))
      stop(paste0("It's not allowed to use a character or a factor as a response vector for \n",
                  "the ",sQuote(family$family)," family (-Inf < y < Inf)."))
    if (any(is.na(y)))
      stop(paste0("NA is not allowed for the ",sQuote(family$family)," family."))
    yy <- y
  }
  y <- yy

  ## Optimizations for beta's (and sigma2), where theta_hat is: b0, b1, ..., bp (and sigma2)
  if (family$family == "binomial") {
    if (is.null(initial)) {
      initial_beta_sig  <- c(rep(1, p))
    } else {
      if (length(initial)==p & (inherits(initial, "numeric") | inherits(initial, "integer"))) {
        initial_beta_sig  <- initial
      } else {
        stop(paste("'initial' should be a numerical vector of length",sQuote(p),"."))
      }
    }
    theta_optim <- try(optim(initial_beta_sig, fn=negloglik.beta.sigma, gr=NULL,
                             y=y, X=X, Lambda=Lambda, family=family$family, lower=rep(-Inf,p),
                             upper=rep(Inf,p), method="L-BFGS-B", control=control), TRUE)
    if(!is.null(attr(theta_optim,"class")) | inherits(theta_optim, "try-error")) {
      #conv <- FALSE
      stop("The algorithm did not converge.")
    }
    if (theta_optim$convergence == 0)
      conv <- 0 #paste("The algorithm has converged.")
    if (theta_optim$convergence == 1)
      conv <- 1 #paste("The iteration limit 'maxit' had been reached.")
    if (theta_optim$convergence %in% c(51, 52))
      conv <- noquote(paste(" 2  (Message from the 'L-BFGS-B' method: ''", theta_optim$message, "'')"))
    theta_hat <- theta_optim$par
    value_theta_hat <- theta_optim$value
    names(theta_hat) <- colnames(X)
  }
  if (family$family == "gaussian") {
    if (is.null(initial)) {
      initial_beta_sig  <- c(rep(1, p+1))
    } else {
      if (length(initial)==(p+1) & (inherits(initial, "numeric") | inherits(initial, "integer"))) {
        initial_beta_sig  <- initial
      } else {
        stop(paste("'initial' should be a numerical vector of length", sQuote(p+1),"."))
      }
    }
    theta_optim <- try(optim(initial_beta_sig, fn=negloglik.beta.sigma, gr=NULL,
                             y=y, X=X, Lambda=Lambda, family=family$family, lower=c(rep(-Inf,p-1),0),
                             upper=rep(Inf,p), method="L-BFGS-B", control=control), TRUE)
    if(!is.null(attr(theta_optim,"class")) | inherits(theta_optim, "try-error")) {
      #conv <- FALSE
      stop("The algorithm did not converge.")
    }
    if (theta_optim$convergence == 0)
      conv <- 0 #paste("The algorithm has converged.")
    if (theta_optim$convergence == 1)
      conv <- 1 #paste("The iteration limit 'maxit' had been reached.")
    if (theta_optim$convergence %in% c(51, 52))
      conv <- noquote(paste(" 2  (Message from the 'L-BFGS-B' method: ''", theta_optim$message, "'')"))
    theta_hat <- theta_optim$par
    value_theta_hat <- theta_optim$value # this is before transferring back the sigma2
    # Now, theta_hat returns sigma2
    theta_hat[length(theta_hat)] <- exp(2 * theta_hat[length(theta_hat)]) # Now, theta_hat returns sigma2
    names(theta_hat) <- c(colnames(X), "sigma2")
  }
  ## A_hat: curvature matrix estimator
  A_hat  <- A.l.maker(y=y, X=X, Lambda=Lambda, family=family$family, theta_hat=theta_hat)
  ## sd of A_hat
  sd_A <- sqrt(diag(solve(as.matrix(A_hat))))
  # if (intercept == F) {
  #   theta_hat <- theta_hat
  #   names(theta_hat)[1] <- c("(Intercept)")
  #   if (family$family == "gaussian") {
  #     A_hat0 <- matrix(0, p+2, p+2)
  #   }
  #   if (family$family == "binomial") {
  #     A_hat0 <- matrix(0, p+1, p+1)
  #   }
  #   A_hat0[-1,-1] <- A_hat
  #   colnames(A_hat0) <- rownames(A_hat0) <- c("(Intercept)", colnames(A_hat))
  #   A_hat <- A_hat0
  #   ## sd of A_hat
  #   sd_A <- c(0L, sd_A)
  #   names(sd_A)[1] <- c("(Intercept)")
  #   Lambda <- Lambda_orig
  # }
  names <- names(theta_hat)
  if (is.null(colnames(Lambda))) {
    colnames(Lambda) <- rownames(Lambda) <- names
  }
  output <- list(theta_hat=theta_hat, A_hat=A_hat, sd=sd_A, Lambda=Lambda, formula=formula1,
                 names=names, n=n, np=p, value=value_theta_hat, family=family$family,
                 intercept=intercept, convergence=conv, control=control)
  class(output) <- "bfi"
  return(output)
}
