## This file created by Hassan Pazira

#' @export
#' @importFrom stats update.formula


MAP.estimation <- function(y, X, family = c("gaussian", "binomial", "survival"),
                           Lambda, intercept = TRUE,
                           basehaz = c("weibul","exp","gomp","poly","pwexp","unspecified"),
                           treatment = NULL, treat_round = NULL, refer_treat, gamma_bfi = NULL,
                           RCT_propens = NULL, initial = NULL, alpha = 0.1, max_order = 2,
                           n_intervals = 4, min_max_times, center_zero_sample = FALSE,
                           zero_sample_cov, refer_cat, zero_cat, control = list()) {
  # Lambda is the 'inverse' covariance matrix (could be a matrix or list of matrices)!
  if (is.null(control$maxit)) control$maxit <- 100
  family <- match.arg(family)
  # if (is.matrix(X))
  #   warning(
  #     "Note: if in 'X' there is a 'categorical' covariate with more than 2 levels, \n",
  #     "then 'X' MUST be a 'data.frame' instead of a 'matrix'.")
  n <- NROW(X)
  if (NROW(y) != NROW(X)) stop("NROW(y) != NROW(X)")
  if (dim(Lambda)[1] != dim(Lambda)[2])
    stop("'Lambda' should be a square matrix.")
  if (any(diag(Lambda) <= 0) & all(Lambda[row(Lambda) != col(Lambda)] == 0)) {
    # second arg. is not needed
    stop("Diagonal elements of Lambda matrix should be positive (> 0)")
  }
  if ((!is.null(treatment)) & is.null(treat_round)) {
    stop("Which round? treat_round = 'first' or treat_round = 'second'?")
  }
  if ((is.null(treatment)) & !is.null(treat_round)) {
    stop("Both 'treatment' and 'treat_round' must be either both NULL or both not NULL.")
  }
  if ((!is.null(treatment)) & !is.null(treat_round)) {
    if (treat_round=="first" & (!is.null(gamma_bfi) | !is.null(RCT_propens))) {
      stop("When treat_round = 'first', 'gamma_bfi' and 'RCT_propens' must be 'NULL'.")
    }
    if (treat_round=="second" & (is.null(gamma_bfi) & is.null(RCT_propens))) {
      stop("When treat_round = 'second', both 'gamma_bfi' and 'RCT_propens' cannot be 'NULL'.")
    }
  }
  if ((!is.null(treatment)) & !is.null(treat_round)) { # or treat_round=="first"
    if (length(treatment) > 1)
      stop("Only one covariate for treatment.")
    if (!is.character(treatment)) stop("Covariate '", treatment,"' should be a character.")
    if (treat_round=="first") {
      family <- c("binomial")
    }
  }
  wli <- rep(1, n)
  if (family %in% c("binomial", "gaussian")) {
    if (all(as.numeric(X[, 1]) == rep(1, n))) {
      X <- X[, -1, drop = FALSE]
    }
    if (is.null(colnames(X))) {
      #colnames(X) <- paste0("X", seq_len(ncol(X)))
      stop("Column names of X cannot be NULL.")
    } else {
      if (any(c("Intercept", "(Intercept)") %in% colnames(X))) {
        stop("'Intercept' should be the first column of 'X'.")
      }
    }
    #X <- X[,sort(colnames(X))]
    if (!is.null(treatment)) {
      if (!treatment %in% colnames(X)) stop("Treatment should be included in X as a covariate.")
      if (nlevels(factor(X[, treatment])) != 2)
        stop("Covariate '", treatment,"' must have only two categories, e.g., 0 and 1.")
      X_treat_old <- X
      wich_treat <- which(colnames(X) == treatment)
      # Z_only_treat <- as.matrix(as.numeric(X_treat_old[, wich_treat])) #!!!
      Z_only_treat <- X_treat_old[, wich_treat]
      X_no_treat <- X[, -wich_treat, drop = FALSE] # no intercept!
      X <- X_no_treat
      if (treat_round=="first") {
        y <- as.factor(Z_only_treat)
      }
    }
    design_matrix <- paste(colnames(X), collapse = " + ")
    formula <- as.formula(paste("y", design_matrix, sep = " ~ "))
    y_X_vars <- model.maker(formula, as.data.frame(X), family)
    # formula1 <- noquote(paste("y", paste(colnames(y_X_vars$X)[-1],
    #                                      collapse=" + "), sep=" ~ "))
    if (is.null(treatment)) {
      formula1 <- noquote(paste("y", design_matrix, sep = " ~ "))
    } else {
      formula1 <- noquote(paste(treatment, design_matrix, sep = " ~ "))
    }
    X <- y_X_vars$X
    y <- y_X_vars$y
    if (NROW(y) != NROW(X)) stop("NROW(y) != NROW(X)")
    if (intercept == FALSE) {
      X <- X[, -1, drop = FALSE]
      if (any(c("Intercept", "(Intercept)") %in% colnames(X))) {
        stop("'Intercept' should be the first column of 'X'.")
      }
    }
    if (intercept == TRUE) {
      xvar <- apply(X[, -1, drop = FALSE], 2, sd)
    } else {
      xvar <- apply(X, 2, sd)
    }
    xvar[xvar < 10 * .Machine$double.eps] <- 0
    const_vars <- sqrt(xvar) == 0
    if (any(const_vars)) {
      if (intercept == TRUE) {
        warning(
          "Predictor(s) ",
          paste0(colnames(X[, -1, drop = FALSE])[which(const_vars)], sep = " "),
          "has zero variance! The 'bfi()' function can be used to address this variable(s)."
        )
      } else {
        warning(
          "Predictor(s) ", paste0(colnames(X)[which(const_vars)], sep = " "),
          "has zero variance! The 'bfi()' function can be used to address this variable(s)."
        )
      }
    }
    p <- ncol(X) # number of regression parameters/coefficients, i.e.,
    # predictors with intercept if intercept=T, or only predictors if intercept=F.
    if (family == "binomial") {
      if (NCOL(y) == 1) {
        if (length(y) != n) stop("length of 'y' != sample size 'n'.")
        if (is.character(y)) y <- factor(y)
        if (is.factor(y)) {
          if (nlevels(y) != 2) stop("only factors with two levels are allowed.")
          yy <- as.numeric(y) - 1
        } else {
          if (any(y < 0 | y > 1)) stop("values of 'y' must be 0 <= y <= 1")
          if (any(abs(y - round(y)) > 0.001))
            warning("non-integer #successes in a binomial.") #
          yy <- y
        }
      } else if (NCOL(y) == 2) {
        if (dim(y)[1] != n) stop("sample size != 'n'.")
        if (any(abs(y - round(y)) > 0.001))
          warning("non-integer counts in a binomial.")
        nn <- y[, 1] + y[, 2]
        y <- ifelse(nn == 0, 0, y[, 1] / nn)
        yy <- y
      } else {
        stop(paste(
          "For the 'binomial' family, 'y' must be",
          "a vector of 0 and 1\'s or a 2 column matrix \n",
          "where col 1 is number successes",
          "and col 2 is number failures"
        ))
      }
    }
    if (family == "gaussian") {
      if (NCOL(y) != 1) {
        stop(paste0("For the ", sQuote(family),
                    " family, 'y' must be a vector, not a matrix."))
      }
      if (length(y) != n) {
        stop("length of 'y' != sample size 'n'")
      }
      if (is.character(y)) {
        stop(paste0(
          "It's not allowed to use a character or a factor as a response \n",
          "vector for the ", sQuote(family), " family (-Inf < y < Inf)."
        ))
      }
      if (any(is.na(y))) {
        stop(paste0("NA is not allowed for the ", sQuote(family), " family."))
      }
      yy <- y
    }
    y <- yy
    if (!is.null(names(y))) y <- unname(y) #!

    if ((!is.null(treatment)) & !is.null(treat_round) &
        ((!is.null(gamma_bfi)) | (!is.null(RCT_propens)))) { # or treat_round=="second"
      if (treat_round=="second" & (is.null(gamma_bfi)) & (is.null(RCT_propens)))
        stop("In the second round, both 'gamma_bfi' and 'RCT_propens' cannot be 'NULL'.")
      if ((!is.null(gamma_bfi)) & (!is.null(RCT_propens)))
        stop("When 'treatment' is not 'NULL', one of 'RCT_propens' or 'gamma_bfi' should be 'NULL'. ",
             "In the first round, both should be 'NULL'.")
      # if (!is.null(gamma_bfi)) gamma_bfi <- as.data.frame(gamma_bfi)
      if (!treatment %in% colnames(Lambda)) stop("Treatment should be included in 'Lambda'.")

      if (missing(refer_treat)) refer_treat <- levels(as.factor(Z_only_treat))[1]
      Zli <- relevel(as.factor(Z_only_treat), ref = refer_treat)  # Set refer_treat as reference.
      Zli_0_1 <- as.numeric(Zli) - 1  # Convert to 0/1, so that 0 is reference.
      if (length(Zli_0_1) != length(y)) # It is not necessary!
        stop("The length of treatment should be equal to the number of individuals, length(y).")
      X_no_treat <- X
      if (intercept == TRUE) {
        X <- cbind(X[, 1], Zli_0_1)
        colnames(X) <- c("(Intercept)", treatment)
      } else {
        X <- cbind(Zli_0_1)
        colnames(X) <- treatment
      }
      p <- NCOL(X)
      if (!is.null(initial)) {
        if (family == "binomial") initial <- initial[1:p]
        if (family == "gaussian") initial <- initial[c(1:p,length(initial))]
      }
      if (!is.null(gamma_bfi)) {
        if (is.null(names(gamma_bfi))) {
          if (!is.null(colnames(gamma_bfi))) names(gamma_bfi) <- colnames(gamma_bfi)
          else stop("'gamma_bfi' should have names. names(gamma_bfi) should not be 'NULL'.")
        }
        if ("(Intercept)" %in% names(gamma_bfi)) {
          if (length(gamma_bfi) != dim(X_no_treat)[2]) {
            if (intercept == FALSE) stop("In the first round, intercept == TRUE, while in the second it is FALSE.")
            stop("Check the length of 'gamma_bfi'. 'gamma_bfi' should consist of regression coefficients and/or Intercept, excluding nuisance parameters.")
          }
          if (!identical(names(gamma_bfi), colnames(X_no_treat)))
            stop("names(gamma_bfi) and X[, -which(colnames(X) == treatment)] differ in elements, order, or both.")
        } else {
          if (length(gamma_bfi) != dim(X_no_treat)[2]) {
            if (intercept == TRUE) stop("In the first round, intercept == FALSE, while in the second it is TRUE.")
            stop("Check the length of 'gamma_bfi'. 'gamma_bfi' should consist of regression coefficients and/or Intercept, excluding nuisance parameters.")
          }
          if (!identical(names(gamma_bfi), colnames(X_no_treat)))
            stop("names(gamma_bfi) and X[, -which(colnames(X) == treatment)] differ in elements, order, or both.")
        }
        if (!identical(names(gamma_bfi), colnames(X_no_treat)))
          stop("names(gamma_bfi) and X[, -which(colnames(X) == treatment)] differ in elements, order, or both.")
        beta_bfi <- gamma_bfi
        eli <- exp(as.numeric(as.matrix(X_no_treat) %*% as.numeric(t(beta_bfi))))/
            (1 + exp(as.numeric(as.matrix(X_no_treat) %*% as.numeric(t(beta_bfi)))))
      }
      if (!is.null(RCT_propens)) {
        if (length(RCT_propens) != nrow(X))
          stop("Length of 'RCT_propens' should be n_l.")
        eli <- RCT_propens
      }
      if (length(eli) != nrow(X)) # or length(y)
        stop("The length of Propensity Scores should be equal to the number of individuals, length(y).")
      m_1 <- sum(Zli_0_1) # the number of patients in the treatment
      m_2 <- nrow(X) - m_1  # the number of patients in the reference group
      sum_z_y <- sum(Zli_0_1 * y)
      sum_z_y2 <- sum((Zli_0_1 * y)^2)
      sum_z_e <- sum(Zli_0_1 / eli)
      sum_z_y_e <- sum(Zli_0_1 * y / eli)
      sum_1_z_1_e <- sum((1 - Zli_0_1) / (1 - eli))
      sum_1_z_y_1_e <- sum((1 - Zli_0_1) * y / (1 - eli))
      sum_1_z_y <- sum((1 - Zli_0_1) * y)
      for_ATE <- c(m_1, m_2, sum_z_y, sum_z_y2, sum_z_e, sum_z_y_e, sum_1_z_1_e, sum_1_z_y_1_e, sum_1_z_y)
      propensity <- eli
      wli_ATE <- (Zli_0_1/eli) + ((1-Zli_0_1)/(1-eli))
      wli_ATT <- Zli_0_1 + (eli * (1-Zli_0_1)/(1-eli))
      wli <- wli_ATE #!
    } else {
      gamma_bfi <- NULL
      RCT_propens <- NULL
      for_ATE <- NULL
      propensity <- NULL
      if ((!is.null(treatment)) & is.null(gamma_bfi) & is.null(RCT_propens) & missing(refer_treat))
        refer_treat <- NULL
    }
    if (family == "binomial") {
      if (dim(Lambda)[1] != p) {
        stop(paste(
          "'Lambda' in this family should have p =", sQuote(p), "dimensions where",
          "'p' is number of coefficients with", "\n",
          "intercept (if intercept=T) or without intercept (if intercept=F))."
        ))
      }
    }
    if (family == "gaussian") {
      if (dim(Lambda)[1] != (p + 1)) {
        stop(paste(
          "'Lambda' should have p+1 =", sQuote(p + 1), "dimensions where",
          " p =", sQuote(p), "is number of coefficients with", "\n",
          " intercept (if intercept=T) or without intercept (if intercept=F),",
          "and '1' is for error variance."
        ))
      }
    }

    ## Optimizations for beta's (and sigma2),
    ## where theta_hat is: b0, b1, ..., bp (and sigma2)

    if (family == "gaussian") p <- p + 1
    if (is.null(initial)) {
      initial_beta_sig <- c(rep(0, p))
    } else {
      if (length(initial) == p & (inherits(initial, "numeric") |
                                  inherits(initial, "integer"))) {
        initial_beta_sig <- initial
      } else {
        stop(paste("'initial' should be a numerical vector of length",
                   sQuote(p), "."))
      }
    }
    if (length(initial_beta_sig) != dim(Lambda)[1])
      stop("Dim. of 'Lambda' is ", dim(Lambda)[1],"x",dim(Lambda)[2],
           ", while the length of 'theta' is ",length(initial_beta_sig),".")
    theta_optim <- try(optim(initial_beta_sig, fn = negloglik.theta, gr = NULL,
                             y = y, X = X, Lambda = Lambda, family = family,
                             wli = wli, lower = c(rep(-Inf, p)),
                             upper = rep(Inf, p), method = "L-BFGS-B",
                             control = control), TRUE)
    if (!is.null(attr(theta_optim, "class")) |
        inherits(theta_optim, "try-error")) {
      # conv <- FALSE
      stop("The algorithm did not converge. Change 'initial' values.")
    }
    if (theta_optim$convergence == 0) {
      conv <- 0
    } # paste("The algorithm has converged.")
    if (theta_optim$convergence == 1) {
      conv <- 1
    } # paste("The iteration limit 'maxit' had been reached.")
    if (theta_optim$convergence %in% c(51, 52)) {
      conv <- noquote(paste(" 2  (Message from the 'L-BFGS-B' method: ''",
                            theta_optim$message, "'')"))
    }
    theta_hat <- theta_optim$par
    value_theta_hat <- theta_optim$value
    if (family == "binomial") {
      names(theta_hat) <- colnames(X)
    }
    if (family == "gaussian") {
      # this is before transferring back to sigma2
      # Now, theta_hat returns log(sigma2)
      theta_hat[length(theta_hat)] <- exp(theta_hat[length(theta_hat)])
      # Now, theta_hat returns sigma2
      names(theta_hat) <- c(colnames(X), "sigma2")
      p <- p - 1
    }

    ## A_hat: curvature matrix estimator
    A_hat <- A.l.maker(y = y, X = X, Lambda = Lambda, family = family,
                       theta_hat = theta_hat, wli = wli)
    ## sd of A_hat
    sd_A <- sqrt(diag(solve(as.matrix(A_hat))))
  }
  if (family == "survival") {
    intercept <- FALSE
    q_l <- NULL
    basehaz <- match.arg(basehaz)
    if (all(as.numeric(X[, 1]) == rep(1, n))) {
      X <- X[, -1, drop = FALSE]
    }
    if (is.null(colnames(X))) {
      #colnames(X) <- paste0("X",seq_len(ncol(X)))
      stop("Column names of X cannot be NULL.")
    }
    #X <- X[,sort(colnames(X))]
    X_old <- X
    if ((!is.null(treatment))) {
      if (!treatment %in% colnames(X))
        stop("Treatment should be included in X as a covariate.")
      if (nlevels(factor(X[, treatment])) != 2)
        stop("Covariate '", treatment,"' must have only two categories, e.g., 0 and 1.")
      wich_treat <- which(colnames(X) == treatment)
      X_no_treat <- X[, -wich_treat, drop = FALSE] # no intercept!
      X <- X_no_treat # This 'X' or 'X_no_treat' is without dummy variables!
    }
    design_matrix <- paste(colnames(X), collapse=" + ")
    formula <- as.formula(paste(" ", design_matrix, sep=" ~ "))
    X_vars <- model.maker(formula, as.data.frame(X), family="survival")
    if (is.null(treatment)) {
      formula1 <- noquote(paste("y", design_matrix, sep = " ~ "))
    } else {
      formula1 <- noquote(paste("y", treatment, sep = " ~ "))
    }
    X <- X_vars$X[,-1, drop=F]
    if (NROW(y) != NROW(X)) stop("NROW(y) != NROW(X)")
    p <- ncol(X)
    xvar <- apply(X, 2, sd)
    xvar[xvar < 10 * .Machine$double.eps] <- 0
    const_vars <- sqrt(xvar) == 0
    if (any(const_vars)) {
      warning("Predictor(s) ", paste0(colnames(X)[which(const_vars)],sep=" "), "has zero variance!",
           "The function 'bfi()' can be used to accommodate this variable(s).")
    }
    if (NCOL(y) == 2) {
      if (dim(y)[1] != n) stop("number of rows of 'y' != sample size 'n'.")
      if (is.character(y$status)) y$status <- factor(y$status)
      if (is.factor(y$status)) {
        if (nlevels(y$status) != 2) stop("only factors with two levels are allowed.")
        y$status <- as.numeric(y$status) - 1 # if character
      } else {
        if (any(!y$status %in% c(0, 1))) stop("values of 'status' must be 0 and/or 1")
        if (any(abs(y[,2] - round(y[,2])) > 0.001))
          warning("non-integer counts in 'status'.")
      }
      if (any(is.na(y$status)) | any(is.na(y$time)))
        stop(paste0("NA is not allowed for 'time' or 'status'."))
    } else {
      stop(paste(
        "For the 'survival' family, 'y' must be a 2-column matrix, where the first column \n",
        "      is 'time' and the second is 'status' with 0 (censored) and 1 (event)."
      ))
    }

    if (!is.null(treatment)) {
      if (treat_round=="second" & (is.null(gamma_bfi)) & (is.null(RCT_propens)))
        stop("In the second round, both 'gamma_bfi' and 'RCT_propens' cannot be 'NULL'.")
      if ((!is.null(gamma_bfi)) & (!is.null(RCT_propens)))
        stop("When 'treatment' is not 'NULL', one of 'RCT_propens' or 'gamma_bfi' should be 'NULL'. ",
             "In the first round, both should be 'NULL'.")
      if (!is.character(treatment))
        stop("Covariate '", treatment,"' should be a character.")
      # if (!is.null(gamma_bfi)) gamma_bfi <- as.data.frame(gamma_bfi)
      if (!treatment %in% colnames(Lambda))
        stop("Treatment should be included in 'Lambda'.") #!

      # Z_only_treat <- as.matrix(as.numeric(X_old[, wich_treat])) #!!!
      Z_only_treat <- X_old[, wich_treat]
      if (missing(refer_treat)) refer_treat <- levels(as.factor(Z_only_treat))[1]
      Zli <- relevel(as.factor(Z_only_treat), ref = refer_treat)  # Set refer_treat as reference.
      Zli_0_1 <- as.numeric(Zli) - 1  # Convert to 0/1, so that 0 is reference.
      if (length(Zli_0_1) != nrow(y)) # It is not necessary!
        stop("The length of treatment should be equal to the number of individuals, nrow(y).")
      X_no_treat <- X # This 'X' or 'X_no_treat' is with dummy variables!
      X <- as.matrix(Zli_0_1)
      colnames(X) <- treatment
      p <- NCOL(X)
      if (!is.null(gamma_bfi)) {
        if (is.null(names(gamma_bfi))) {
          if (!is.null(colnames(gamma_bfi))) names(gamma_bfi) <- colnames(gamma_bfi)
          else stop("'gamma_bfi' should have names. names(gamma_bfi) should not be 'NULL'.")
        }
        if ("(Intercept)" %in% names(gamma_bfi)) {
          if ((length(gamma_bfi)-1) != dim(X_no_treat)[2]) {
            stop("Check the length of 'gamma_bfi'. 'gamma_bfi' should consist of regression coefficients and/or Intercept, excluding nuisance parameters.")
          }
          if (!identical(names(gamma_bfi)[-1], colnames(X_no_treat)))
            stop("names(gamma_bfi) and X[, -which(colnames(X) == treatment)] differ in elements, order, or both.")
        } else {
          if ((length(gamma_bfi)) != dim(X_no_treat)[2]) {
            stop("Check the length of 'gamma_bfi'. 'gamma_bfi' should consist of regression coefficients and/or Intercept, excluding nuisance parameters.")
          }
          if (!identical(names(gamma_bfi), colnames(X_no_treat)))
            stop("names(gamma_bfi) and X[, -which(colnames(X) == treatment)] differ in elements, order, or both.")
        }
        beta_bfi <- gamma_bfi
        if ("(Intercept)" %in% names(gamma_bfi)) {
          eli <- exp(as.numeric(as.matrix(cbind(1, X_no_treat)) %*% as.numeric(t(beta_bfi))))/
            (1 + exp(as.numeric(as.matrix(cbind(1, X_no_treat)) %*% as.numeric(t(beta_bfi)))))
        } else {
          eli <- exp(as.numeric(as.matrix(X_no_treat) %*% as.numeric(t(beta_bfi))))/
            (1 + exp(as.numeric(as.matrix(X_no_treat) %*% as.numeric(t(beta_bfi)))))
        }
      }
      if (!is.null(RCT_propens)) {
        if (length(RCT_propens) != nrow(X))
          stop("Length of 'RCT_propens' should be n_l.")
        eli <- RCT_propens
      }
      if (length(eli) != nrow(X))
        stop("The length of Propensity Scores should be equal to the number of individuals, nrow(y).")
      for_ATE <- NULL # for survival 'for_ATE' is NULL!
      propensity <- eli
      wli_ATE <- (Zli_0_1/eli) + ((1-Zli_0_1)/(1-eli))
      wli_ATT <- Zli_0_1 + (eli * (1-Zli_0_1)/(1-eli))
      wli <- wli_ATE # maybe 'wli_ATT' in the future!
    } else {
      gamma_bfi <- NULL
      RCT_propens <- NULL
      for_ATE <- NULL
      propensity <- NULL
      if ((!is.null(treatment)) & is.null(gamma_bfi) & is.null(RCT_propens) & missing(refer_treat))
        refer_treat <- NULL
    }

    ## theta_hat, A_hat and sd_A:
    if (basehaz == "poly") {
      if (length(max_order) > 1 | max_order < 1)
        stop("'max_order' should be an integer >= 1.")
      theta_len <- p + max_order + 1 # maximum number
      theta_M_ql <- array(data = NA, dim = c(theta_len, theta_len, max_order+2))
      rownames(theta_M_ql) <- c(colnames(X), paste("omega",c(0:(max_order)), sep="_"))
      if (is.null(initial)) {
        initial_ini <- c(rep(0.0, p+1))
      } else {
        if (length(initial)>=(p+1) &
            (inherits(initial, "numeric") | inherits(initial, "integer"))) {
          initial_ini <- initial[1:(p+1)]
        } else {
          stop(paste("'initial' should be a numerical vector of length >= ", sQuote(p+1),"."))
        }
      }
      # Creating 'q_l'
      q_l <- try(ql.LRT(initial_ini, y=y, X=X, alpha=alpha, max_order=max_order), TRUE)
      if(!is.null(attr(q_l, "class")) | inherits(q_l, "try-error")) {
        q_l <- max_order
        # warning("q_l <- max_order")
      }
      if (q_l > max_order) q_l <- max_order
      conv <- value_theta_hat <- NULL
      if (dim(Lambda)[1] < (p+max_order+1))
        stop("Dimension of 'Lambda' is less than ", (p+max_order+1),".")
      for (ql in q_l:max_order) {
        Gamma_dot <- Lambda[1:(p+ql+1),1:(p+ql+1)] # This Lambda should be p+max_order+1 dim!
        ### Optimizations for theta's
        if (is.null(initial)) {
          initial_beta_omega <- c(rep(0.0, p), rep(-0.5, ql+1))
        } else {
          if (length(initial)==(p+ql+1) &
              (inherits(initial, "numeric") | inherits(initial, "integer"))) {
            initial_beta_omega <- initial
          } else {
            stop(paste("'initial' should be a numerical vector of length", sQuote(p+ql+1),"."))
          }
        }
        if (length(initial_beta_omega) != dim(Gamma_dot)[1])
          stop("Dim. of 'Lambda' is ", dim(Gamma_dot)[1],"x",dim(Gamma_dot)[2],
               ", while the length of 'theta' is ",length(initial_beta_omega),".")
        theta_for_M_ql <- optim.survival(initial_beta_loga_logb=initial_beta_omega,
                                    q_l=ql, tps=tps, y=y, X=X, Lambda=Gamma_dot,
                                    family=family, wli = wli, basehaz=basehaz,
                                    control = control)
        if(!is.null(attr(theta_for_M_ql, "class")) |
           inherits(theta_for_M_ql, "try-error")) {
          stop("The algorithm did not converge. Change 'initial' values.")
        }
        theta_M_ql[1:(p+ql+1),ql+1,1] <- theta_for_M_ql$par
        # A_hat
        A_l_maker_new <- A.l.maker(y=y, X=X, Lambda=Gamma_dot, q_l=ql, tps=tps,
                                   basehaz=basehaz, family=family, wli = wli,
                                   theta_hat=theta_for_M_ql$par)
        theta_M_ql[1:(p+ql+1),1:(p+ql+1),ql+2] <- A_l_maker_new
        if (ql == q_l) {
          A_hat <- A_l_maker_new
          sd_A <- sqrt(diag(solve(as.matrix(A_hat))))
          theta_hat <- theta_for_M_ql$par
          names(theta_hat) <- colnames(A_hat)
        }
        if (theta_for_M_ql$convergence == 0) {
          conv_ql <- 0
        } # paste("The algorithm has converged.")
        if (theta_for_M_ql$convergence == 1) {
          conv_ql <- 1
        } # paste("The iteration limit 'maxit' had been reached.")
        if (theta_for_M_ql$convergence %in% c(51, 52)) {
          conv_ql <- noquote(paste(" 2  (Message from the 'L-BFGS-B' method: ''",
                                theta_for_M_ql$message, "'')"))
        }
        conv <- c(conv, conv_ql)
        value_theta_hat <- c(value_theta_hat, theta_for_M_ql$value)
      }
      # Zero-patient category
      if (center_zero_sample == TRUE) {
        if (missing(zero_sample_cov))
          stop("If there is a categorical covariate with zero sample in only one of",
               "these categories, 'zero_sample_cov' should not be NULL.")
        if (length(zero_sample_cov) > 1)
          stop("For now, only one covariate for each center!")
        if (is.numeric(zero_sample_cov)) {
          zero_sample_cov <- colnames(X_old)[zero_sample_cov]
          # X_old is the original data without 'time' and 'status'.
        }
        if (nlevels(factor(X_old[, zero_sample_cov])) < 2)
          stop("Covariate '", zero_sample_cov,"' must have at least two categories,",
               "one of which is the reference.")
        levs_no_ref_no_zero_poly <- levels(relevel(factor(X_old[, zero_sample_cov]),
                                             ref = refer_cat))
      } else levs_no_ref_no_zero_poly <- NULL
    } else if (basehaz == "unspecified") {
      theta_len <- p

      ## Optimizations for theta's
      if (is.null(initial)) {
        initial_beta_omega <- c(rep(0.0, p), rep(-0.5,theta_len-p))
      } else {
        if (length(initial)==(theta_len) &
            (inherits(initial, "numeric") | inherits(initial, "integer"))) {
          initial_beta_omega <- initial
        } else {
          stop(paste("'initial' should be a numerical vector of length",
                     sQuote(theta_len),"."))
        }
      }
      if (length(initial_beta_omega) != dim(Lambda)[1])
        stop("Dim. of 'Lambda' is ", dim(Lambda)[1],"x",dim(Lambda)[2],
             ", while the length of 'theta' is ",length(initial_beta_omega),".")
      theta_optim <- optim.survival(initial_beta_loga_logb=initial_beta_omega,
                                    q_l=q_l, tps=tps, y=y, X=X, Lambda=Lambda,
                                    family=family, basehaz=basehaz, wli=wli,
                                    control = control)
      if (!is.null(attr(theta_optim, "class")) | inherits(theta_optim, "try-error")) {
        # conv <- FALSE
        stop("The algorithm did not converge. Change 'initial' values.")
      }
      if (theta_optim$convergence == 0) {
        conv <- 0
      } # paste("The algorithm has converged.")
      if (theta_optim$convergence == 1) {
        conv <- 1
      } # paste("The iteration limit 'maxit' had been reached.")
      if (theta_optim$convergence %in% c(51, 52)) {
        conv <- noquote(paste(" 2  (Message from the 'L-BFGS-B' method: ''",
                              theta_optim$message, "'')"))
      }
      theta_hat <- theta_optim$par
      value_theta_hat <- theta_optim$value
      ## A_hat: curvature matrix estimator
      A_hat <- A.l.maker(y = y, X = X, Lambda = Lambda, family = family,
                         theta_hat = theta_hat, q_l = q_l, tps=tps,
                         basehaz = basehaz, wli = wli)
      names(theta_hat) <- colnames(A_hat)
      if (any(diag(solve(as.matrix(A_hat))) < 0))
        stop("In this center, some diagonal elements of the inverse curvature ",
             "matrix are negative! Increase the sample size, e.g., > 50.")
      ## sd_A
      sd_A <- sqrt(diag(solve(as.matrix(A_hat))))

      # Zero-patient category
      if (center_zero_sample == TRUE) { # It needs to be updated for propensity scores !
        if (missing(zero_sample_cov))
          stop("If there is a categorical covariate with zero sample in only one of",
               "these categories, 'zero_sample_cov' should not be NULL.")
        if (length(zero_sample_cov) != 1)
          stop("For now only one covariate for each center!")
        if (is.numeric(zero_sample_cov)) {
          zero_sample_cov <- colnames(X_old)[zero_sample_cov]
          # X_old is the original data without 'time' and 'status'.
        }
        if (missing(refer_cat)) stop("The reference category is missing.")
        if (length(zero_sample_cov) != length(refer_cat))
          stop("'zero_sample_cov' and 'refer_cat' must have the same length.")
        if (!is.character(zero_cat)) stop("'zero_cat' should be character.")
        if (!is.character(refer_cat)) stop("'refer_cat' should be character.")
        if (missing(zero_cat)) stop("The category with zero sample is missing.")
        if (any(zero_cat == refer_cat))
          stop("The category with no patient cannot be used as the reference.")
        # For now only one covariate and one 'zero_cat'!
        for (zz in 1:length(zero_sample_cov)) {
          which_zz <- zero_sample_cov[zz]
          ZZ_old2 <- factor(X_old[, which_zz])
          if (length(levels(ZZ_old2)) < 2)
            stop("Covariate '", which_zz,"' must have at least two categories,",
                 "one of which is the reference.")
          ZZ_old2 <- relevel(ZZ_old2, ref = refer_cat[zz])
          # cat("The categorical covariate '", which_zz,"' had",
          #     length(levels(factor(X_old[,which_zz]))),
          #     "levels (",levels(factor(X_old[,which_zz])),"), but now it has",
          #     length(zero_cat), "additional category (",zero_cat,
          #     ") with no sample in it.","\n")
          lev_zero_cov <- sort(as.character(c(levels(ZZ_old2)[-1], zero_cat)))
          names_after_dammy <- paste(which_zz,lev_zero_cov, sep="")
          which_cat_zero <- which(lev_zero_cov %in% zero_cat)
          names_cat_zero <- names_after_dammy[which_cat_zero]
          if (length(names_cat_zero) != length(zero_cat))
            stop("length(names_cat_zero) != length(zero_cat)")
          # Positions to add the new elements
          for (wi in seq_len(length(zero_cat))) {
            # This 'for' loop should be updated for more than one 'zero_cat'!
            if (wi == 1) which_element <- NULL
            if (which_cat_zero == 1) {
              which_element[wi] <- which(colnames(X) == names_after_dammy[2])
            } else {
              if (length(lev_zero_cov) > 2) {
                if (length(names_after_dammy) == which_cat_zero) {
                  which_element[wi] <- which(colnames(X) ==
                                               names_after_dammy[which_cat_zero-1]) - 1
                } else {
                  which_element[wi] <- which(colnames(X) ==
                                               names_after_dammy[which_cat_zero + 1])
                }
              } else {
                which_element[wi] <- which(colnames(X) == names_after_dammy[1]) + 1
              }
            }
          }
          names_initial_vec_beta <- colnames(X)
          # Add the new elements to the vector
          for (i in seq_len(length(which_element))) { # or seq_along(names_cat_zero)
            names_initial_vec_beta <- append(names_initial_vec_beta, names_cat_zero[i],
                                             after = which_element[i] - 1)
          }
          initial_vec_beta <- rep(0, (ncol(X) + length(names_cat_zero)))
          names(initial_vec_beta) <- names_initial_vec_beta
          initial_vec_beta[-which_element] <- theta_hat[1:p]
          initial_vec_beta <- c(initial_vec_beta,
                                theta_hat[(p+1):length(theta_hat)])
          theta_hat <- initial_vec_beta

          # A_hat
          # create a new 'Gamma_dot' for all parameters!
          Gamma <- Lambda[c(1:p),c(1:p)]
          Gamma_new <- Gamma
          if (ncol(Gamma_new) < which_element[1] &
              abs(ncol(Gamma_new) - which_element[1])==1) {
            for (ii in seq_len(length(which_element))) {
              Gamma_new <- rbind(Gamma_new,
                                 c(Gamma_new[ncol(Gamma_new),1],
                                   Gamma_new[ncol(Gamma_new),-ncol(Gamma_new)]))
            }
            # Add the columns at the specified position to A_hat
            for (jj in seq_len(length(which_element))) {
              Gamma_new <- cbind(Gamma_new,
                                 c(Gamma_new[nrow(Gamma_new),],
                                   Gamma_new[nrow(Gamma_new)-1, ncol(Gamma_new)]))
            }
          } else {
            for (ii in seq_len(length(which_element))) {
              Gamma_new <- rbind(Gamma_new[1:(which_element[ii] - 1), ],
                                 c(Gamma_new[which_element[ii]+1,1:which_element[ii]],
                                   Gamma_new[which_element[ii]+1,which_element[ii]],
                                   Gamma_new[which_element[ii]+1,(which_element[ii]+2):ncol(Gamma_new)]),
                                 Gamma_new[which_element[ii]:nrow(Gamma_new), ])
            }
            # Add the columns at the specified position to A_hat
            for (jj in seq_len(length(which_element))) {
              Gamma_new <- cbind(Gamma_new[, 1:(which_element[jj] - 1)],
                                 c(Gamma[,which_element[jj]],
                                   Gamma[ncol(Gamma),which_element[jj]]),
                                 Gamma_new[, which_element[jj]:ncol(Gamma_new)])
            }
          }
          if (basehaz != "poly")
            Gamma_dot_new <- b.diag(Gamma_new, Lambda[-c(1:p),-c(1:p)])

          # Add the rows at the specified position to A_hat
          new_M_with_rows <- A_hat
          for (ii in seq_len(length(which_element))) {
            new_M_with_rows <- rbind(new_M_with_rows[1:(which_element[ii] - 1), ],
                                     Gamma_dot_new[which_element[ii],-which_element[ii]],
                                     new_M_with_rows[which_element[ii]:nrow(new_M_with_rows), ])
          }
          # Add the columns at the specified position to A_hat
          new_M_with_columns <- new_M_with_rows
          for (jj in seq_len(length(which_element))) {
            new_M_with_columns <- cbind(new_M_with_columns[, 1:(which_element[jj] - 1)],
                                        Gamma_dot_new[,which_element[jj]],
                                        new_M_with_columns[, which_element[jj]:ncol(new_M_with_columns)])
            colnames(new_M_with_columns)[which_element[jj]] <- names_cat_zero[jj]
          }
          A_hat <- new_M_with_columns
          if (any(diag(solve(as.matrix(A_hat)))<0))
            stop("In this center, some diagonal elements of the inverse curvature",
                 "matrix are negative! Increase the sample size, e.g., > 50.")
          sd_A <- sqrt(diag(solve(as.matrix(A_hat))))
        }
        Lambda <- Gamma_dot_new
      }
      # levs_no_ref_no_zero_poly <- NULL
    } else {
      if (basehaz == "exp") {
        theta_len <- p + 1
      }
      if (basehaz %in% c("gomp", "weibul")) {
        theta_len <- p + 2
      }
      if (basehaz=="pwexp") {
        theta_len <- p + n_intervals
        if (missing(min_max_times))
          stop("The minimum of the maximum times is missing.")
        tps <- seq(0, min_max_times + 0.001, length.out = n_intervals + 1)
        which.int.obs <- matrix(data = NA, nrow = nrow(X), ncol = length(tps) - 1)
        dd <- matrix(data = NA, nrow = nrow(X), ncol = length(tps) - 1)
        for(i in 1:nrow(X)){
          for(k in 1:(length(tps) - 1)){
            dd[i, k] <- ifelse(y$time[i] - tps[k] > 0, 1, 0) * ifelse(tps[k + 1] - y$time[i] > 0, 1, 0)
            which.int.obs[i, k] <- dd[i, k] * k
          }
        }
        int.obs <- rowSums(which.int.obs)
        # this vector tells us at which interval each observation is.
        n.obs.int <- colSums(dd) # no. observations in each interval.
        cat("\n", "No. observations in the intervals : ", n.obs.int, "\n", "\n")
      }

      ## Optimizations for theta's
      if (is.null(initial)) {
        initial_beta_omega <- c(rep(0.0, p), rep(-0.5,theta_len-p))
      } else {
        if (length(initial)==(theta_len) &
            (inherits(initial, "numeric") | inherits(initial, "integer"))) {
          initial_beta_omega <- initial
        } else {
          stop(paste("'initial' should be a numerical vector of length",
                     sQuote(theta_len),"."))
        }
      }
      if (length(initial_beta_omega) != dim(Lambda)[1])
        stop("Dim. of 'Lambda' is ", dim(Lambda)[1],"x",dim(Lambda)[2],
             ", while the length of 'theta' is ",length(initial_beta_omega),".")
      theta_optim <- optim.survival(initial_beta_loga_logb=initial_beta_omega,
                               q_l=q_l, tps=tps, y=y, X=X, Lambda=Lambda, family=family,
                               wli = wli, basehaz = basehaz, control = control)
      if (!is.null(attr(theta_optim, "class")) | inherits(theta_optim, "try-error")) {
        # conv <- FALSE
        stop("The algorithm did not converge. Change 'initial' values.")
      }
      if (theta_optim$convergence == 0) {
        conv <- 0
      } # paste("The algorithm has converged.")
      if (theta_optim$convergence == 1) {
        conv <- 1
      } # paste("The iteration limit 'maxit' had been reached.")
      if (theta_optim$convergence %in% c(51, 52)) {
        conv <- noquote(paste(" 2  (Message from the 'L-BFGS-B' method: ''",
                              theta_optim$message, "'')"))
      }
      theta_hat <- theta_optim$par
      value_theta_hat <- theta_optim$value
      ## A_hat: curvature matrix estimator
      A_hat <- A.l.maker(y = y, X = X, Lambda = Lambda, family = family,
                         theta_hat = theta_hat, q_l = q_l, tps = tps,
                         wli = wli, basehaz = basehaz)
      names(theta_hat) <- colnames(A_hat)
      if (any(diag(solve(as.matrix(A_hat))) < 0))
        stop("In this center, some diagonal elements of the inverse curvature ",
             "matrix are negative! Increase the sample size, e.g., > 50.")
      ## sd_A
      sd_A <- sqrt(diag(solve(as.matrix(A_hat))))
      # Zero-patient category
      if (center_zero_sample == TRUE) {
        if (missing(zero_sample_cov))
          stop("If there is a categorical covariate with zero sample in only one of",
               "these categories, 'zero_sample_cov' should not be NULL.")
        if (length(zero_sample_cov) != 1)
          stop("For now only one covariate for each center!")
        if (is.numeric(zero_sample_cov)) {
          zero_sample_cov <- colnames(X_old)[zero_sample_cov]
          # X_old is the original data without 'time' and 'status'.
        }
        if (missing(refer_cat)) stop("The reference category is missing.")
        if (length(zero_sample_cov) != length(refer_cat))
          stop("'zero_sample_cov' and 'refer_cat' must have the same length.")
        if (!is.character(zero_cat)) stop("'zero_cat' should be character.")
        if (!is.character(refer_cat)) stop("'refer_cat' should be character.")
        if (missing(zero_cat)) stop("The category with zero sample is missing.")
        if (any(zero_cat == refer_cat))
          stop("The category with no patient cannot be used as the reference.")
        # For now only one covariate and one 'zero_cat'!
        for (zz in 1:length(zero_sample_cov)) {
          which_zz <- zero_sample_cov[zz]
          ZZ_old2 <- factor(X_old[, which_zz])
          if (length(levels(ZZ_old2)) < 2)
            stop("Covariate '", which_zz,"' must have at least two categories,",
                 "one of which is the reference.")
          ZZ_old2 <- relevel(ZZ_old2, ref = refer_cat[zz])
          # cat("The categorical covariate '", which_zz,"' had",
          #     length(levels(factor(X_old[,which_zz]))),
          #     "levels (",levels(factor(X_old[,which_zz])),"), but now it has",
          #     length(zero_cat), "additional category (",zero_cat,
          #     ") with no sample in it.","\n")
          lev_zero_cov <- sort(as.character(c(levels(ZZ_old2)[-1], zero_cat)))
          names_after_dammy <- paste(which_zz,lev_zero_cov, sep="")
          which_cat_zero <- which(lev_zero_cov %in% zero_cat)
          names_cat_zero <- names_after_dammy[which_cat_zero]
          if (length(names_cat_zero) != length(zero_cat))
            stop("length(names_cat_zero) != length(zero_cat)")
          # Positions to add the new elements
          for (wi in seq_len(length(zero_cat))) {
            # This 'for' loop should be updated for more than one 'zero_cat'!
            if (wi == 1) which_element <- NULL
            if (which_cat_zero == 1) {
              which_element[wi] <- which(colnames(X) == names_after_dammy[2])
            } else {
              if (length(lev_zero_cov) > 2) {
                if (length(names_after_dammy) == which_cat_zero) {
                  which_element[wi] <- which(colnames(X) ==
                                               names_after_dammy[which_cat_zero-1]) - 1
                } else {
                  which_element[wi] <- which(colnames(X) ==
                                               names_after_dammy[which_cat_zero + 1])
                }
              } else {
                which_element[wi] <- which(colnames(X) == names_after_dammy[1]) + 1
              }
            }
          }
          names_initial_vec_beta <- colnames(X)
          # Add the new elements to the vector
          for (i in seq_len(length(which_element))) { # or seq_along(names_cat_zero)
            names_initial_vec_beta <- append(names_initial_vec_beta, names_cat_zero[i],
                                             after = which_element[i] - 1)
          }
          initial_vec_beta <- rep(0, (ncol(X) + length(names_cat_zero)))
          names(initial_vec_beta) <- names_initial_vec_beta
          initial_vec_beta[-which_element] <- theta_hat[1:p]
          initial_vec_beta <- c(initial_vec_beta,
                                theta_hat[(p+1):length(theta_hat)])
          theta_hat <- initial_vec_beta

          # A_hat
          # create a new 'Gamma_dot' for all parameters!
          Gamma <- Lambda[c(1:p),c(1:p)]
          Gamma_new <- Gamma
          if (ncol(Gamma_new) < which_element[1] &
              abs(ncol(Gamma_new) - which_element[1])==1) {
            for (ii in seq_len(length(which_element))) {
              Gamma_new <- rbind(Gamma_new,
                                 c(Gamma_new[ncol(Gamma_new),1],
                                   Gamma_new[ncol(Gamma_new),-ncol(Gamma_new)]))
            }
            # Add the columns at the specified position to A_hat
            for (jj in seq_len(length(which_element))) {
              Gamma_new <- cbind(Gamma_new,
                                 c(Gamma_new[nrow(Gamma_new),],
                                   Gamma_new[nrow(Gamma_new)-1, ncol(Gamma_new)]))
            }
          } else {
            for (ii in seq_len(length(which_element))) {
              Gamma_new <- rbind(Gamma_new[1:(which_element[ii] - 1), ],
                                 c(Gamma_new[which_element[ii]+1,1:which_element[ii]],
                                   Gamma_new[which_element[ii]+1,which_element[ii]],
                                   Gamma_new[which_element[ii]+1,(which_element[ii]+2):ncol(Gamma_new)]),
                                 Gamma_new[which_element[ii]:nrow(Gamma_new), ])
            }
            # Add the columns at the specified position to A_hat
            for (jj in seq_len(length(which_element))) {
              Gamma_new <- cbind(Gamma_new[, 1:(which_element[jj] - 1)],
                                 c(Gamma[,which_element[jj]],
                                   Gamma[ncol(Gamma),which_element[jj]]),
                                 Gamma_new[, which_element[jj]:ncol(Gamma_new)])
            }
          }
          if (basehaz != "poly")
            Gamma_dot_new <- b.diag(Gamma_new, Lambda[-c(1:p),-c(1:p)])

          # Add the rows at the specified position to A_hat
          new_M_with_rows <- A_hat
          for (ii in seq_len(length(which_element))) {
            new_M_with_rows <- rbind(new_M_with_rows[1:(which_element[ii] - 1), ],
                                     Gamma_dot_new[which_element[ii],-which_element[ii]],
                                     new_M_with_rows[which_element[ii]:nrow(new_M_with_rows), ])
          }
          # Add the columns at the specified position to A_hat
          new_M_with_columns <- new_M_with_rows
          for (jj in seq_len(length(which_element))) {
            new_M_with_columns <- cbind(new_M_with_columns[, 1:(which_element[jj] - 1)],
                                        Gamma_dot_new[,which_element[jj]],
                                        new_M_with_columns[, which_element[jj]:ncol(new_M_with_columns)])
            colnames(new_M_with_columns)[which_element[jj]] <- names_cat_zero[jj]
          }
          A_hat <- new_M_with_columns
          if (any(diag(solve(as.matrix(A_hat)))<0))
            stop("In this center, some diagonal elements of the inverse curvature",
                 "matrix are negative! Increase the sample size, e.g., > 50.")
          sd_A <- sqrt(diag(solve(as.matrix(A_hat))))
        }
        Lambda <- Gamma_dot_new
      }
      # levs_no_ref_no_zero_poly <- NULL
    }
  }
  if (is.null(colnames(Lambda))) {
    if (length(theta_hat) != dim(Lambda)[1])
      Lambda <- Lambda[1:length(theta_hat),1:length(theta_hat)]
    colnames(Lambda) <- rownames(Lambda) <- names(theta_hat)
  }
  if (center_zero_sample == FALSE) {
    zero_sample_cov = NULL
    refer_cat = NULL
    zero_cat = NULL
  }
  if (family == "survival") {
    # formula1 <- as.formula(deparse(update.formula(formula1, Survival(time, status) ~ .)))
    # formula1 <- deparse(update.formula(formula1, Survival(time, status) ~ .))
    formula1 <- gsub("\\s+", " ", paste(deparse(update.formula(formula1, Survival(time, status) ~ .)), collapse = ""))
    if (basehaz == "poly") {
      output <- list(
        theta_hat = theta_hat, A_hat = A_hat, sd = sd_A, Lambda = Lambda,
        formula = formula1, names = names(theta_hat), n = n, np = p, q_l = q_l,
        theta_A_poly = theta_M_ql, lev_no_ref_zero = levs_no_ref_no_zero_poly,
        treatment = treatment, zero_sample_cov = zero_sample_cov, refer_cat = refer_cat,
        zero_cat = zero_cat, value = value_theta_hat, family = family,
        basehaz = basehaz, intercept = intercept, convergence = conv,
        control = control
      )
    } else {
      if (!is.null(treatment)) {
        output <- list(
          theta_hat = theta_hat, A_hat = A_hat, sd = sd_A, Lambda = Lambda,
          formula = formula1, names = names(theta_hat), n = n, np = p,
          treatment = treatment, refer_treat = refer_treat, gamma_bfi = gamma_bfi,
          RCT_propens = RCT_propens, propensity = propensity, for_ATE = for_ATE,
          zero_sample_cov = zero_sample_cov, refer_cat = refer_cat, zero_cat = zero_cat,
          value = value_theta_hat, family = family, basehaz = basehaz, intercept = intercept,
          convergence = conv, control = control
        )
      } else {
        output <- list(
          theta_hat = theta_hat, A_hat = A_hat, sd = sd_A, Lambda = Lambda,
          formula = formula1, names = names(theta_hat), n = n, np = p,
          treatment = treatment, zero_sample_cov = zero_sample_cov, refer_cat = refer_cat,
          zero_cat = zero_cat, value = value_theta_hat, family = family,
          basehaz = basehaz, intercept = intercept, convergence = conv,
          control = control
        )
      }
    }
  } else {
    if (!is.null(treatment)) {
      output <- list(
        theta_hat = theta_hat, A_hat = A_hat, sd = sd_A, Lambda = Lambda,
        formula = formula1, names = names(theta_hat), n = n, np = p,
        treatment = treatment, refer_treat = refer_treat, gamma_bfi = gamma_bfi,
        RCT_propens = RCT_propens, propensity = propensity, for_ATE = for_ATE,
        zero_sample_cov = zero_sample_cov, refer_cat = refer_cat, zero_cat = zero_cat,
        value = value_theta_hat, family = family, basehaz = NULL, intercept = intercept,
        convergence = conv, control = control
      )
    } else {
      output <- list(
        theta_hat = theta_hat, A_hat = A_hat, sd = sd_A, Lambda = Lambda,
        formula = formula1, names = names(theta_hat), n = n, np = p,
        treatment = treatment, zero_sample_cov = zero_sample_cov, refer_cat = refer_cat,
        zero_cat = zero_cat, value = value_theta_hat, family = family,
        basehaz = basehaz, intercept = intercept, convergence = conv,
        control = control
      )
    }
  }
  class(output) <- "bfi"
  return(output)
}
