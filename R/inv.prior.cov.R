## This file created by Hassan Pazira at 20-06-2023
## Updated at 09-07-2023

inv.prior.cov <- function(X, lambda=1, family=gaussian, intercept=TRUE, independ=TRUE, set_seed=NULL) {
  if (any(lambda <= 0)) stop("'lambda' element(s) should be positive (> 0)")
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- do.call(family, args = list(), envir = parent.frame())
  if (is.null(family$family)) stop("'family' not recognized")
  if (!family$family %in% c("binomial", "gaussian")) {
    stop("Distributions that can be used are 'binomial' and 'gaussian' in this version of the package!")
  }
  if (independ==F & !is.null(set_seed)) set.seed(set_seed)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X",1:ncol(X))
  }
  if (is.matrix(X)) warning("Note that, if in 'X' there is a 'categorical' covariate with more than 2 levels,","\n",
                            "then 'X' must be a 'data.frame' instead of a 'matrix'! ")
  y <- rnorm(NROW(X)) # a fake response
  design_matrix <- paste(colnames(X), collapse=" + ")
  formula <- as.formula(paste("y", design_matrix, sep=" ~ "))
  X <- model.maker(formula, as.data.frame(X))$X
  X <- X[, -1, drop = FALSE] # without 'intercept'
  if (intercept==T) {
    name_Lambda <- c("(Intercept)", colnames(X), "sigma2")
  } else {
    name_Lambda <- c(colnames(X), "sigma2")
  }
  np <- NCOL(X)  # without intercept
  if (family$family=="gaussian") {
    if (intercept==T) {
      p <- np+2    # intercept and error variance
    } else {
      p <- np+1    # only error variance (no intercept)
    }
    if (length(lambda) == 1) {
      lambda1 <- lambda2 <- lambda
    }
    if (length(lambda) == 2) {
      lambda1 <- lambda[1]
      lambda2 <- lambda[2]
    }
    if (length(lambda) > 2) {
      if (all(lambda == lambda[1])) {
        lambda1 <- lambda2 <- lambda[1]
      } else {
        if (length(lambda) != p) {
          stop("'lambda' should be a vector of ", sQuote(p), " elements")
        } else {
          lambda1 <- lambda[1:(p-1)]
          lambda2 <- lambda[p]
        }
      }
    }
    if (length(lambda1)==1) {
      lambda1 <- rep(lambda1, p-1)
    } else {
      lambda1 <- lambda1
    }
    if (independ == T) {
      Lambda <- diag(c(lambda1, lambda2), p)
    } else {
      Lambda <- matrix(rnorm(p*p),p)
      #diag(Lambda) <- c(lambda1, lambda2)
    }
    rownames(Lambda) <- colnames(Lambda) <- name_Lambda
  }
  if (family$family=="binomial") {
    if (intercept==T) {
      p <- np+1 # only intercept (and of course no error variance)
    } else {
      p <- np # no intercept (and of course no error variance)
    }
    if (length(lambda) == 1) {
      lambda1 <- lambda
    }
    if (length(lambda) >= 2) {
      if (all(lambda == lambda[1])) {
        lambda1 <- lambda[1]
      } else {
        if (length(lambda) != p) {
          stop("'lambda' should be a vector of ", sQuote(p), " elements")
        } else {
          lambda1 <- lambda
        }
      }
    }
    if (independ == T) {
      Lambda <- diag(ifelse(length(lambda1)==1, rep(lambda1, p), lambda1), p)
    } else {
      Lambda <- matrix(rnorm(p*p),p)
      #diag(Lambda) <- ifelse(length(lambda1)==1, rep(lambda1, p), lambda1)
    }
    rownames(Lambda) <- colnames(Lambda) <- name_Lambda[1:dim(Lambda)[1]]
  }
  return(Lambda)
}
