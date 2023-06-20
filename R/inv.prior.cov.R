## This file created by Hassan Pazira at 20-06-2023
## Updated at 29-09-2023

inv.prior.cov <- function(X, lambda=1, family=gaussian, intercept=TRUE, stratified=FALSE, strat_par=NULL, L=2L) {
  if (any(lambda <= 0)) stop("'lambda' element(s) should be positive (> 0)")
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- do.call(family, args = list(), envir = parent.frame())
  if (is.null(family$family)) stop("'family' not recognized")
  if (!family$family %in% c("binomial", "gaussian")) {
    stop("Distributions that can be used are 'binomial' and 'gaussian' in this version of the package!")
  }
  if (is.matrix(X)) warning("Note that, if in 'X' there is a 'categorical' covariate with more than 2 levels,","\n",
                            "then 'X' must be a 'data.frame' instead of a 'matrix'! ")
  if (all(as.numeric(X[,1])==rep(1,NROW(X)))) {
    X <- X[,-1, drop=F]
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X",1:ncol(X))
  } else {
    if (any(c("Intercept", "(Intercept)") %in% colnames(X)))
      stop("'intercept' should be the first column of 'X'.")
  }
  y <- rnorm(NROW(X)) # a fake response
  design_matrix <- paste(colnames(X), collapse=" + ")
  formula <- as.formula(paste("y", design_matrix, sep=" ~ "))
  X <- model.maker(formula, as.data.frame(X))$X
  if (family$family=="gaussian") {
    name_Lambda <- c(colnames(X), "sigma2")
  }
  if (family$family=="binomial") {
    name_Lambda <- colnames(X)
  }
  X <- X[, -1, drop = FALSE] # without 'intercept'
  np <- NCOL(X)  # without intercept
  if (stratified == F) {
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
      Lambda <- diag(c(lambda1, lambda2), length(c(lambda1, lambda2)))
      if (intercept==T) {
        rownames(Lambda) <- colnames(Lambda) <- name_Lambda
      } else {
        rownames(Lambda) <- colnames(Lambda) <- name_Lambda[-1]
      }
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
      if (length(lambda1)==1) {
        dia_lam <- rep(lambda1, p)
      } else {
        dia_lam <- lambda1
      }
      Lambda <- diag(dia_lam, length(dia_lam))
      if (intercept==T) {
        rownames(Lambda) <- colnames(Lambda) <- name_Lambda
      } else {
        rownames(Lambda) <- colnames(Lambda) <- name_Lambda[-1]
      }
    }
  } else {
    if (is.null(strat_par))
      stop("In stratified analysis, 'strat_par' should not be 'NULL'. Select the stratification parameter(s).")
    if (! 1 %in% strat_par & ! 2 %in% strat_par) {
      stop("'strat_par' should contain '1' and/or '2'.")
    }
    if (family$family == "gaussian") {
      if (!is.numeric(strat_par))
        stop("'strat_par' should be one of the integers: 1 or 2, or a vector of both.")
      if (length(strat_par) > 2)
        stop("For the 'gaussian' family the number of 'strat_par' parameters is two, i.e., 'intercept' and 'sigma2'.")
      if (intercept == F) {
        if (length(strat_par) > 1 | 1 %in% strat_par) {
          stop("Since 'intercept = F', for the 'gaussian' family 'strat_par' should only be the integer 2")
        }
      }
      if (length(strat_par) == 1) {
        if (1 %in% strat_par) {   # only for intercept == T
          p <- L + np + 1   # L intercepts + NCOL(X) + sigma2
          if (length(lambda) == 2) {
            lambda1 <- rep(lambda[1], L+np)
            lambda2 <- c()
            lambda3 <- rep(lambda[2], 1)
          }
          if (length(lambda) == 3) {
            lambda1 <- rep(lambda[1], L)
            lambda2 <- rep(lambda[2], np)
            lambda3 <- rep(lambda[3], 1)
          }
          if (length(lambda) > 3) {
            if (all(lambda == lambda[1])) {
              lambda1 <- rep(lambda[1], p)
              lambda2 <- lambda3 <- c()
            } else {
              if (length(lambda) != p) {
                stop("'lambda' should be a vector of ", sQuote(p), " elements")
              } else {
                lambda1 <- lambda
                lambda2 <- lambda3 <- c()
              }
            }
          }
        } else {
          if (intercept==T) {
            p <- 1 + np + L   # 1 intercept + NCOL(X) + L sigma2's
            if (length(lambda) == 2) {
              lambda1 <- rep(lambda[1], 1+np)
              lambda2 <- c()
              lambda3 <- rep(lambda[2], L)
            }
            if (length(lambda) == 3) {
              lambda1 <- rep(lambda[1], 1)
              lambda2 <- rep(lambda[2], np)
              lambda3 <- rep(lambda[3], L)
            }
            if (length(lambda) > 3) {
              if (all(lambda == lambda[1])) {
                lambda1 <- rep(lambda[1], p)
                lambda2 <- lambda3 <- c()
              } else {
                if (length(lambda) != p) {
                  stop("'lambda' should be a vector of ", sQuote(p), " elements")
                } else {
                  lambda1 <- lambda
                  lambda2 <- lambda3 <- c()
                }
              }
            }
          } else {
            p <- L + np       # NCOL(X) + L sigma2's
            if (length(lambda) == 2) {
              lambda1 <- c()
              lambda2 <- rep(lambda[1], np)
              lambda3 <- rep(lambda[2], L)
            }
            if (length(lambda) == 3) {
              stop("'lambda' can be a vector of ", sQuote(1), ", ", sQuote(2), " or ", sQuote(p), " elements")
            }
            if (length(lambda) > 3) {
              if (all(lambda == lambda[1])) {
                lambda1 <- rep(lambda[1], p)
                lambda2 <- lambda3 <- c()
              } else {
                if (length(lambda) != p) {
                  stop("'lambda' should be a vector of ", sQuote(p), " elements")
                } else {
                  lambda1 <- lambda
                  lambda2 <- lambda3 <- c()
                }
              }
            }
          }
        }
      } else { # length(strat_par) == 2
        p <- 2*L + np        # L intercepts + NCOL(X) + L sigma2's
        if (length(lambda) == 2) {
          lambda1 <- rep(lambda[1], L+np)
          lambda2 <- c()
          lambda3 <- rep(lambda[2], L)
        }
        if (length(lambda) == 3) {
          lambda1 <- rep(lambda[1], L)
          lambda2 <- rep(lambda[2], np)
          lambda3 <- rep(lambda[3], L)
        }
        if (length(lambda) > 3) {
          if (all(lambda == lambda[1])) {
            lambda1 <- rep(lambda[1], p)
            lambda2 <- lambda3 <- c()
          } else {
            if (length(lambda) != p) {
              stop("'lambda' should be a vector of ", sQuote(p), " elements")
            } else {
              lambda1 <- lambda
              lambda2 <- lambda3 <- c()
            }
          }
        }
      }
      if (length(lambda) == 1) {
        lambda1 <- rep(lambda, p)
        lambda2 <- lambda3 <- c()
      }
      Lambda <- diag(c(lambda1, lambda2, lambda3), length(c(lambda1, lambda2, lambda3)))
      if (intercept==T) {
        if (length(strat_par) == 1) {
          if (1 %in% strat_par) {
            rownames(Lambda) <- colnames(Lambda) <- c(paste0(name_Lambda[1], rep("_loc", L), 1:L),
                                                      name_Lambda[-c(1, length(name_Lambda))],
                                                      name_Lambda[length(name_Lambda)])
          }
          if (2 %in% strat_par) {
            rownames(Lambda) <- colnames(Lambda) <- c(name_Lambda[-length(name_Lambda)],
                                                      paste0(name_Lambda[length(name_Lambda)], rep("_loc", L), 1:L))
          }
        } else {
          rownames(Lambda) <- colnames(Lambda) <- c(paste0(name_Lambda[1], rep("_loc", L), 1:L),
                                                    name_Lambda[-c(1, length(name_Lambda))],
                                                    paste0(name_Lambda[length(name_Lambda)], rep("_loc", L), 1:L))
        }
      } else {
        rownames(Lambda) <- colnames(Lambda) <- c(name_Lambda[-c(1, length(name_Lambda))],
                                                  paste0(name_Lambda[length(name_Lambda)], rep("_loc", L), 1:L))
      }
    }
    if (family$family == "binomial") {
      if (intercept == F) {
        stop("Since 'intercept = F', for the 'binomial' family the stratified analysis is not possible!")
      }
      if (!is.numeric(strat_par))
        stop("'strat_par' should be an integer")
      if (length(strat_par) > 1)
        stop("For the 'binomial' family the number of 'strat_par' parameters is one, i.e., 'strat_par = 1'.")
      if (strat_par == 2)
        stop("'strat_par' should only be the integer 1, i.e., 'strat_par = 1'.")
      if (strat_par == 1) { # Must be, as it is the only case!
        p <- L + np   # L intercepts + NCOL(X)
      }
      lambda1 <- c()
      if (length(lambda) == 1) {
        lambda1 <- lambda
      }
      if (length(lambda) == 2) {
        lambda1[1] <- lambda[1]
        lambda1[2] <- lambda[2]
      }
      if (length(lambda) > 2) {
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
      if (length(lambda1)==1) {
        dia_lam <- rep(lambda1, p)
      } else {
        if (length(lambda1)==2) {
          dia_lam <- c(rep(lambda1[1], L), rep(lambda1[2], np))
        } else {
          dia_lam <- lambda1
        }
      }
      Lambda <- diag(dia_lam, length(dia_lam))
      rownames(Lambda) <- colnames(Lambda) <- c(paste0(name_Lambda[1], rep("_loc", L), 1:L), name_Lambda[-1])
    }
  }
  return(Lambda)
}
