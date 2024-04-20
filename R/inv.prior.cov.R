## This file created by Hassan Pazira

#' @export

inv.prior.cov <- function(X, lambda = 1, L = 2L,
                          family = c("gaussian", "binomial", "survival"),
                          treatment = NULL, treat_round = NULL, intercept = TRUE,
                          stratified = FALSE, strat_par = NULL, center_spec = NULL,
                          basehaz = c("weibul", "exp", "gomp", "poly", "pwexp", "unspecified"),
                          max_order = 2, n_intervals = 4) {
  family <- match.arg(family)
  if (any(lambda <= 0)) stop("'lambda' element(s) should be positive (> 0)")
  if (is.null(family)) stop("'family' not recognized")
  if (is.matrix(X)) {
    warning(
      "If in 'X' there is a 'categorical' covariate with more than 2 levels,",
      " then 'X' must be a 'data.frame' instead of a 'matrix'! "
    )
  }
  if ((!is.null(treatment)) & is.null(treat_round)) {
    stop("Which round? treat_round = 'first' or treat_round = 'second'?")
  }
  if ((is.null(treatment)) & !is.null(treat_round)) {
    stop("Both 'treatment' and 'treat_round' must be either both NULL or both not NULL.")
  }
  if ((!is.null(treatment))) { # & is.null(gamma_bfi) & is.null(RCT_propens)
    if (length(treatment) > 1)
      stop("Only one covariate for treatment.")
    if (treat_round=="first") {
      family <- c("binomial")
    }
    if (!is.character(treatment)) stop("Covariate '", treatment,"' should be a character.")
  }
  if (is.null(colnames(X))) {
    # colnames(X) <- paste0("X", seq_len(ncol(X)))
    if (length(X[,-1])==nrow(X)) stop()
    else stop("Colnames of X cannot be NULL.")
  }
  if (family %in% c("binomial", "gaussian")) {
    if (all(as.numeric(X[, 1]) == rep(1, NROW(X)))) {
      X <- X[, -1, drop = FALSE]
    }
    if (any(c("Intercept", "(Intercept)") %in% colnames(X)))
      stop("'Intercept' should be the first column of 'X'.")
  }
  #X <- X[,sort(colnames(X))]
  if ((!is.null(treatment)) ) { #& is.null(gamma_bfi) & is.null(RCT_propens)
    if (!treatment %in% colnames(X)) stop("Treatment should be included in X as a covariate.")
    if (nlevels(factor(X[, treatment])) != 2)
      stop("Covariate '", treatment,"' must have only two categories, e.g., 0 and 1.")
    X_treat_old <- X
    wich_treat <- which(colnames(X) == treatment)
    Z_only_treat <- as.matrix(as.numeric(X_treat_old[, wich_treat]))
    if (treat_round=="second") {
      X <- X[, wich_treat, drop = FALSE]
    } else {
      X <- X[, -wich_treat, drop = FALSE]
    }
  }
  design_matrix <- paste(colnames(X), collapse = " + ")
  formula <- as.formula(paste(" ", design_matrix, sep = " ~ "))
  X <- model.maker(formula, as.data.frame(X), family)$X
  if ((!is.null(treatment)) ) {
    if (treat_round=="second") {
      if (intercept == TRUE) {
        colnames(X) <- c("(Intercept)", treatment)
      } else {
        colnames(X) <- c(treatment)
      }
    }
  }
  if (family == "gaussian") {
    name_Lambda <- c(colnames(X), "sigma2")
  }
  if (family == "binomial") {
    name_Lambda <- colnames(X)
  }
  if (family == "survival") {
    basehaz <- match.arg(basehaz)
    if (basehaz == "unspecified") {
      name_Lambda <- colnames(X)[-1]
    }
    if (basehaz == "exp") {
      name_Lambda <- c(colnames(X), paste("omega",c(1), sep="_"))[-1]
    }
    if (basehaz %in% c("gomp", "weibul")) {
      name_Lambda <- c(colnames(X), paste("omega",c(1:2), sep="_"))[-1]
    }
    if (basehaz=="pwexp") {
      name_Lambda <- c(colnames(X), paste("omega",c(1:(n_intervals)), sep="_"))[-1]
    }
    if (basehaz == "poly") {
      name_Lambda <- c(colnames(X), paste("omega",c(0:(max_order)), sep="_"))[-1]
    }
  }
  X <- X[, -1, drop = FALSE] # without 'intercept'
  np <- NCOL(X) # without intercept
  if (stratified == TRUE & is.null(strat_par) & is.null(center_spec)) {
    stop("Since 'stratified = TRUE', only one of 'strat_par' or 'center_spec' can be NULL.")
  }
  if (stratified == TRUE & !is.null(strat_par) & !is.null(center_spec)) {
    stop("Since 'stratified = TRUE', only one of 'strat_par' or 'center_spec' could be non-NULL.")
  }
  if (stratified == FALSE & (!is.null(strat_par) | !is.null(center_spec))) {
    stop("Since 'stratified = FALSE', both 'strat_par' and 'center_spec' should be NULL.")
  }
  strat_par <- sort(strat_par)
  if (stratified == FALSE) {
    if (family == "gaussian") {
      if (intercept == TRUE) {
        p <- np + 2 # intercept and error variance
      } else {
        p <- np + 1 # only error variance (no intercept)
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
            stop("'lambda' could be a vector of ", sQuote(1),", ",
                 sQuote(2)," or ", sQuote(p), " elements.")
          } else {
            lambda1 <- lambda[seq_len(p - 1)]
            lambda2 <- lambda[p]
          }
        }
      }
      if (length(lambda1) == 1) {
        lambda1 <- rep(lambda1, p - 1)
      } else {
        lambda1 <- lambda1
      }
      Lambda <- diag(c(lambda1, lambda2), length(c(lambda1, lambda2)))
      if (intercept == TRUE) {
        rownames(Lambda) <- colnames(Lambda) <- name_Lambda
      } else {
        rownames(Lambda) <- colnames(Lambda) <- name_Lambda[-1]
      }
    }
    if (family == "binomial") {
      if (intercept == TRUE) {
        p <- np + 1 # only intercept (and of course no error variance)
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
            stop("'lambda' could be a vector of ", sQuote(1),
                 " or ", sQuote(p), " elements.")
          } else {
            lambda1 <- lambda
          }
        }
      }
      if (length(lambda1) == 1) {
        dia_lam <- rep(lambda1, p)
      } else {
        dia_lam <- lambda1
      }
      Lambda <- diag(dia_lam, length(dia_lam))
      if (intercept == TRUE) {
        rownames(Lambda) <- colnames(Lambda) <- name_Lambda
      } else {
        rownames(Lambda) <- colnames(Lambda) <- name_Lambda[-1]
      }
    }
    if (family == "survival") {
      intercept <- FALSE  # In 'survival' there is no intercept !
      if (basehaz == "unspecified") {
        p <- np
      }
      if (basehaz == "exp") {
        p <- np + 1
      }
      if (basehaz %in% c("gomp", "weibul")) {
        p <- np + 2
      }
      if (basehaz=="pwexp") {
        p <- np + n_intervals
      }
      if (basehaz == "poly") {
        p <- np + (max_order + 1)
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
            stop("'lambda' could be a vector of ", sQuote(1),", ",
                 sQuote(2)," or ", sQuote(p), " elements.")
          } else {
            lambda1 <- lambda[seq_len(np)]
            lambda2 <- lambda[(np+1):p]
          }
        }
      }
      if (length(lambda1) == 1) {
        lambda1 <- rep(lambda1, np)
        lambda2 <- rep(lambda2, p - np)
      }
      Lambda <- diag(c(lambda1, lambda2), length(c(lambda1, lambda2)))
      rownames(Lambda) <- colnames(Lambda) <- name_Lambda
    }
  } else {
    if (is.null(center_spec)) {
      if (any(duplicated(strat_par))) stop("There shouldn't be any duplicates in 'strat_par'.")
      if (family == "gaussian") {
        if (!all(strat_par %in% c(1, 2))) {
          stop("'strat_par' should contain '1' and/or '2'.")
        }
        if (!is.numeric(strat_par)) {
          stop("'strat_par' should be one of the integers: 1 or 2, or a vector of both.")
        }
        if (length(strat_par) > 2) {
          stop("For the 'gaussian' family the number of 'strat_par' parameters
               is two, i.e., 'intercept' and 'sigma2'.")
        }
        if (intercept == FALSE) {
          if (length(strat_par) > 1 | 1 %in% strat_par) {
            stop("Since 'intercept = FALSE', for the 'gaussian' family
                 'strat_par' should only be the integer 2")
          }
        }
        if (length(strat_par) == 1) {
          if (1 %in% strat_par) { # only for intercept == TRUE
            p <- L + np + 1 # L intercepts + NCOL(X) + sigma2
            if (length(lambda) == 2) {
              lambda1 <- rep(lambda[1], L + np)
              lambda2 <- NULL
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
                lambda2 <- lambda3 <- NULL
              } else {
                if (length(lambda) != p) {
                  stop("'lambda' should be a vector of ",
                       sQuote(p), " elements")
                } else {
                  lambda1 <- lambda
                  lambda2 <- lambda3 <- NULL
                }
              }
            }
          } else {
            if (intercept == TRUE) {
              p <- 1 + np + L # 1 intercept + NCOL(X) + L sigma2's
              if (length(lambda) == 2) {
                lambda1 <- rep(lambda[1], 1 + np)
                lambda2 <- NULL
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
                  lambda2 <- lambda3 <- NULL
                } else {
                  if (length(lambda) != p) {
                    stop("'lambda' should be a vector of ",
                         sQuote(p), " elements")
                  } else {
                    lambda1 <- lambda
                    lambda2 <- lambda3 <- NULL
                  }
                }
              }
            } else {
              p <- L + np # NCOL(X) + L sigma2's
              if (length(lambda) == 2) {
                lambda1 <- NULL
                lambda2 <- rep(lambda[1], np)
                lambda3 <- rep(lambda[2], L)
              }
              if (length(lambda) == 3) {
                stop("'lambda' can be a vector of ", sQuote(1), ", ", sQuote(2),
                     " or ", sQuote(p), " elements")
              }
              if (length(lambda) > 3) {
                if (all(lambda == lambda[1])) {
                  lambda1 <- rep(lambda[1], p)
                  lambda2 <- lambda3 <- NULL
                } else {
                  if (length(lambda) != p) {
                    stop("'lambda' should be a vector of ", sQuote(p),
                         " elements")
                  } else {
                    lambda1 <- lambda
                    lambda2 <- lambda3 <- NULL
                  }
                }
              }
            }
          }
        } else { # length(strat_par) == 2
          p <- 2 * L + np # L intercepts + NCOL(X) + L sigma2's
          if (length(lambda) == 2) {
            lambda1 <- rep(lambda[1], L + np)
            lambda2 <- NULL
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
              lambda2 <- lambda3 <- NULL
            } else {
              if (length(lambda) != p) {
                stop("'lambda' should be a vector of ", sQuote(p), " elements")
              } else {
                lambda1 <- lambda
                lambda2 <- lambda3 <- NULL
              }
            }
          }
        }
        if (length(lambda) == 1) {
          lambda1 <- rep(lambda, p)
          lambda2 <- lambda3 <- NULL
        }
        Lambda <- diag(c(lambda1, lambda2, lambda3), length(c(lambda1, lambda2, lambda3)))
        if (intercept == TRUE) {
          if (length(strat_par) == 1) {
            if (1 %in% strat_par) {
              rownames(Lambda) <- colnames(Lambda) <- c(
                paste0(name_Lambda[1], rep("_loc", L), seq_len(L)),
                name_Lambda[-c(1, length(name_Lambda))],
                name_Lambda[length(name_Lambda)]
              )
            }
            if (2 %in% strat_par) {
              rownames(Lambda) <- colnames(Lambda) <- c(
                name_Lambda[-length(name_Lambda)],
                paste0(name_Lambda[length(name_Lambda)], rep("_loc", L),
                       seq_len(L))
              )
            }
          } else {
            rownames(Lambda) <- colnames(Lambda) <- c(
              paste0(name_Lambda[1], rep("_loc", L), seq_len(L)),
              name_Lambda[-c(1, length(name_Lambda))],
              paste0(name_Lambda[length(name_Lambda)], rep("_loc", L),
                     seq_len(L))
            )
          }
        } else {
          rownames(Lambda) <- colnames(Lambda) <- c(
            name_Lambda[-c(1, length(name_Lambda))],
            paste0(name_Lambda[length(name_Lambda)], rep("_loc", L), seq_len(L))
          )
        }
      }
      if (family == "binomial") {
        if (intercept == FALSE) {
          stop("Since 'intercept = FALSE', for the 'binomial' family
               the stratified analysis is not possible!")
        }
        if (!is.numeric(strat_par)) {
          stop("'strat_par' should be an integer")
        }
        if (length(strat_par) > 1) {
          stop("For the 'binomial' family the number of 'strat_par' parameters
               is one, i.e., 'strat_par = 1'.")
        }
        if (strat_par == 2) {
          stop("'strat_par' should only be the integer 1, i.e., 'strat_par = 1'.")
        }
        if (strat_par == 1) { # Must be, as it is the only case!
          p <- L + np # L intercepts + NCOL(X)
        }
        lambda1 <- NULL
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
        if (length(lambda1) == 1) {
          dia_lam <- rep(lambda1, p)
        } else {
          if (length(lambda1) == 2) {
            dia_lam <- c(rep(lambda1[1], L), rep(lambda1[2], np))
          } else {
            dia_lam <- lambda1
          }
        }
        Lambda <- diag(dia_lam, length(dia_lam))
        rownames(Lambda) <- colnames(Lambda) <-
          c(paste0(name_Lambda[1], rep("_loc", L), seq_len(L)), name_Lambda[-1])
      }
      if (family == "survival") {
        if (!is.numeric(strat_par)) {
          stop("'strat_par' should contain integers.")
        }
        intercept <- FALSE  # In 'survival' there is no intercept !
        if (basehaz == "unspecified") {
          stop("For the 'unspecified' baseline hazard, 'strat_par' can not be used.")
        }
        if (basehaz == "exp") {
          if (length(strat_par) > 1 | !all(strat_par %in% 1)) {
            stop("For the 'exp' baseline hazard, 'strat_par' should contain only the integer 1.")
          }
          n_omega <- 1
          p <- np + L
        }
        if (basehaz %in% c("gomp", "weibul")) {
          if (length(strat_par) > 2 | !all(strat_par %in% c(1, 2))) {
            stop("For 'gomp' and 'weibul', 'strat_par' can take the values 1, 2, or 1:2.")
          }
          n_omega <- 2
          # p <- np + ifelse(length(strat_par) == 1, L + 1, 2 * L)
          p <- np + length(strat_par) * L + (2 - length(strat_par))
        }
        if (basehaz=="pwexp") {
          if (length(strat_par) > n_intervals | !all(strat_par %in% seq_len(n_intervals))) {
            stop("For 'pwexp', 'strat_par' can be any combination of the values 1 to 'n_intervals'.")
          }
          n_omega <- n_intervals
          p <- np + length(strat_par) * L + (n_intervals - length(strat_par))
        }
        if (basehaz == "poly") {
          if (length(strat_par) > (max_order + 1) | !all(strat_par %in% seq_len(max_order + 1))) {
            stop("For 'poly', 'strat_par' can be any combination of the values 1 to 'n_intervals'.")
          }
          n_omega <- (max_order + 1)
          p <- np + length(strat_par) * L + ((max_order + 1) - length(strat_par))
        }
        if (length(lambda) != p) {
          stop("When 'stratified == TRUE' and 'strat_par != NULL', ",
               "'lambda' should be a vector of ", sQuote(p),
               " elements; length(beta's) + length(strat_par) * L + length(rest of omega's).")
        }
        Lambda <- diag(lambda, length(lambda))
        omegas_names <- NULL
        for (i in seq_len(n_omega)) {
          if (i %in% strat_par) {
            omegas_name <- paste0(name_Lambda[np+i], rep("_loc", L), seq_len(L))
          } else {
            omegas_name <- name_Lambda[np+i]
          }
          omegas_names <- c(omegas_names, omegas_name)
        }
        rownames(Lambda) <- colnames(Lambda) <- c(name_Lambda[1:np], omegas_names)
      }
    } else {
      if (length(center_spec) != L) {
        stop(
          "Either the length of the vector 'center_spec' should be equal to ",
          sQuote(L), " or L should be ", sQuote(length(center_spec))
        )
      }
      if (!is.factor(center_spec)) {
        center_spec <- factor(center_spec)
      }
      lev_cen <- levels(center_spec)
      K <- length(lev_cen)
      if (K < 2) {
        stop("The number of levels of 'center_spec' should be >= 2.")
      }
      if (K == L) {
        stop(
          "The number of levels of 'center_spec' should NOT be equal to ",
          "the number of centers. Otherwise, use 'stratified = TRUE', ",
          "'center_spec = NULL' but set 'strat_par' to 1, 2 or 1:2."
        )
      }
      strat_par <- 1
      # For now, we assume 'strat_par' is 'intercept'.
      # However, the input 'strat_par' is NULL.
      if (length(strat_par) == 1) {
        if (1 %in% strat_par) { # For now, we assume 'strat_par' is 'intercept'.
          strat_par <- 1
          noncore <- cbind(seq_len(K))
          new_noncore <- noncore
        }
        # else {
        #   strat_par <- c(length(theta_hats[[1]]))
        #   noncore <- cbind((last_Lam_dim - K + 1):last_Lam_dim)
        #   new_noncore <- as.matrix(lev_cen)
        # }
      }
      # else {
      #   strat_par <- c(1, length(theta_hats[[1]]))
      #   noncore <- new_noncore <- cbind(lev_cen,
      #                                   (last_Lam_dim - K + 1):last_Lam_dim)
      #   new_noncore[, 2] <- noncore[,1] + K
      # }
      if (family == "gaussian") {
        if (intercept == FALSE) {
          if (length(strat_par) > 1 | 1 %in% strat_par) {
            stop(
              "'intercept' can not ne 'FALSE' for centers, while in the current ",
              "version of the package the center specific covariate is only ",
              "possible for 'intercept'!"
            )
          }
        }
        if (1 %in% strat_par) { # only for intercept == TRUE
          p <- K + np + 1 # K intercepts + NCOL(X) + sigma2
          if (length(lambda) == 2) {
            lambda1 <- rep(lambda[1], K + np)
            lambda2 <- NULL
            lambda3 <- rep(lambda[2], 1)
          }
          if (length(lambda) == 3) {
            lambda1 <- rep(lambda[1], K)
            lambda2 <- rep(lambda[2], np)
            lambda3 <- rep(lambda[3], 1)
          }
          if (length(lambda) > 3) {
            if (all(lambda == lambda[1])) {
              lambda1 <- rep(lambda[1], p)
              lambda2 <- lambda3 <- NULL
            } else {
              if (length(lambda) != p) {
                stop("'lambda' should be a vector of ", sQuote(p), " elements")
              } else {
                lambda1 <- lambda
                lambda2 <- lambda3 <- NULL
              }
            }
          }
        }
        if (length(lambda) == 1) {
          lambda1 <- rep(lambda, p)
          lambda2 <- lambda3 <- NULL
        }
        Lambda <- diag(c(lambda1, lambda2, lambda3),
                       length(c(lambda1, lambda2, lambda3)))
        if (intercept == TRUE) {
          if (length(strat_par) == 1) {
            if (1 %in% strat_par) {
              rownames(Lambda) <- colnames(Lambda) <- c(
                paste0(name_Lambda[1], rep("_", K), lev_cen),
                name_Lambda[-c(1, length(name_Lambda))],
                name_Lambda[length(name_Lambda)]
              )
            }
            if (2 %in% strat_par) {
              # can not be possible for the current version!
              rownames(Lambda) <- colnames(Lambda) <- c(
                name_Lambda[-length(name_Lambda)],
                paste0(name_Lambda[length(name_Lambda)], rep("_", K), lev_cen)
              )
            }
          } else { # can not be possible for the current version!
            rownames(Lambda) <- colnames(Lambda) <- c(
              paste0(name_Lambda[1], rep("_", K), lev_cen),
              name_Lambda[-c(1, length(name_Lambda))],
              paste0(name_Lambda[length(name_Lambda)], rep("_", K), lev_cen)
            )
          }
        } else { # can not be possible for the current version!
          rownames(Lambda) <- colnames(Lambda) <- c(
            name_Lambda[-c(1, length(name_Lambda))],
            paste0(name_Lambda[length(name_Lambda)], rep("_", K), lev_cen)
          )
        }
      }
      if (family == "binomial") {
        if (intercept == FALSE) {
          stop(
            "'intercept = FALSE' for centers, while in the current version of the",
            "package the center specific covariate is only possible for 'intercept'!"
          )
        }
        if (!is.numeric(strat_par)) {
          stop("'strat_par' should be an integer")
        }
        if (length(strat_par) > 1) { # can not be possible for the current version!
          stop("For the 'binomial' family the number of 'strat_par'
               parameters is one, i.e., 'strat_par = 1'.")
        }
        if (strat_par == 2) { # can not be possible for the current version!
          stop("'strat_par' should only be the integer 1, i.e., 'strat_par = 1'.")
        }
        if (strat_par == 1) { # Must be, as it is the only case!
          p <- K + np # K intercepts + NCOL(X)
        }
        lambda1 <- NULL
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
        if (length(lambda1) == 1) {
          dia_lam <- rep(lambda1, p)
        } else {
          if (length(lambda1) == 2) {
            dia_lam <- c(rep(lambda1[1], K), rep(lambda1[2], np))
          } else {
            dia_lam <- lambda1
          }
        }
        Lambda <- diag(dia_lam, length(dia_lam))
        rownames(Lambda) <- colnames(Lambda) <-
          c(paste0(name_Lambda[1], rep("_", K), lev_cen), name_Lambda[-1])
      }
    }
  }
  return(Lambda)
}
