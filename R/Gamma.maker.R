## This file created by Hassan Pazira at 20-06-2023
## Updated at 22-06-2023
Gamma.maker <- function(X, lambda=1, family=c("binomial","gaussian"), independ=T, set_seed=NULL) {
  family <- match.arg(family)
  if (!family %in% c("binomial", "gaussian")) {
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
  name_Gamma <- c("Intercept", colnames(X), "sigma2_e")
  np <- NCOL(X)  # without intercept
  if (family=="gaussian") {
    p <- np+2    # intercept and error variance
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
      Gamma <- diag(c(lambda1, lambda2), p)
    } else {
      Gamma <- matrix(rnorm(p*p),p)
      diag(Gamma) <- c(lambda1, lambda2)
    }
    rownames(Gamma) <- colnames(Gamma) <- name_Gamma
  }
  if (family=="binomial") {
    p <- np+1 # only intercept
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
      Gamma <- diag(ifelse(length(lambda1)==1, rep(lambda1, p), lambda1), p)
    } else {
      Gamma <- matrix(rnorm(p*p),p)
      diag(Gamma) <- ifelse(length(lambda1)==1, rep(lambda1, p), lambda1)
    }
    rownames(Gamma) <- colnames(Gamma) <- name_Gamma[1:dim(Gamma)[1]]
  }
  return(Gamma)
}
