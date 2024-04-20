## This file created by Hassan Pazira

#' @export

A.l.maker <- function(y, X, Lambda, family, theta_hat, q_l, tps, basehaz, wli) {
  nl <- nrow(X)
  p <- ncol(X)
  if (dim(Lambda)[1]!=dim(Lambda)[2])
    stop("'Lambda' should be a square matrix.")
  if (any(diag(Lambda) <= 0))
    stop("Diagonal elements of 'Lambda' matrix should be positive (> 0)")
  if (family == "binomial") {
    len_thet <- length(theta_hat)
    beta_hat <- theta_hat
    power <- matrix(0, len_thet, len_thet)
    for (j in seq_len(nl)) {
      power <- power + t(X[j, , drop = FALSE]) %*% X[j, , drop = FALSE] * wli[j] *
        exp(as.numeric(X[j, ] %*% beta_hat)) /
        (1 + exp(as.numeric(X[j, ] %*% beta_hat)))^2
    }
    A_l <- power + Lambda
    colnames(A_l) <- colnames(X)
  }
  if (family == "gaussian") {
    len_thet <- length(theta_hat)
    beta_hat <- theta_hat[-len_thet]
    sigma2_e_hat <- theta_hat[len_thet]
    power <- matrix(0, len_thet - 1, len_thet - 1)
    power2 <- c(0)
    power3 <- c(rep(0, p))
    for (j in seq_len(nl)) {
      power <- power + t(X[j, , drop = FALSE]) %*% X[j, , drop = FALSE] * wli[j]
      power2 <- power2 + (y[j] - as.numeric(X[j, ] %*% beta_hat))^2 * wli[j]
      power3 <- power3 + (y[j] - as.numeric(X[j, ] %*% beta_hat)) * wli[j] *
        t(X[j, , drop = FALSE])
    }
    A_l <- matrix(0, len_thet, len_thet)
    A_l[-len_thet, -len_thet] <- power / sigma2_e_hat + Lambda[-(p + 1), -(p + 1)]
    A_l[len_thet, len_thet] <- 2 * Lambda[len_thet, len_thet] * sigma2_e_hat +
      2 * power2 / sigma2_e_hat
    A_l[seq_len(len_thet - 1), len_thet] <- A_l[len_thet, seq_len(len_thet - 1)] <-
      -2 * power3 / sigma2_e_hat
    colnames(A_l) <- rownames(A_l) <- c(colnames(X), "sigma2")
  }
  if (family == "survival") {
    if (colnames(y)[1]=="time" & colnames(y)[2]=="status") {
      time <- y$time
      status <- y$status
    } else {
      time <- y[,1]
      status <- y[,2]
    }
    Z <- X
    if (!is.matrix(Z)) Z <- data.matrix(Z) #as.matrix(Z)
    p <- ncol(Z)
    if (basehaz == "poly") {
      if (length(theta_hat) != (p + q_l + 1)) {
        theta_hat <- theta_hat[1:(p + q_l + 1)]
      }
    }
    Gamma <- Lambda[1:p,1:p]
    #Gamma_omeg <- Lambda[-(1:p),-(1:p)]
    beta_hat <- theta_hat[1:p]
    len_theta <- length(theta_hat)
    if (basehaz %in% c("exp", "gomp", "weibul"))
      log_a_hat <- theta_hat[p+1]
    if (basehaz == "gomp" | basehaz == "weibul")
      log_b_hat <- theta_hat[p+2]
    M_l <- matrix(NA, len_theta, len_theta)
    if (basehaz == "pwexp") {
      omega_hat <- theta_hat[-c(1:p)] # it's not a log transformation of omega's
      b <- i.basis(tps, time, ibasis = FALSE) #compute basis splines of degree zero, i.e. piecewise constant basis
      B <- i.basis(tps, time, ibasis = TRUE)
      H_est <- as.numeric(B %*% exp(omega_hat))
      Z_beta_hat = as.numeric(Z %*% beta_hat)
      M_l[1:p, 1:p] <- (t(Z) %*% diag(wli * exp(Z_beta_hat) * H_est) %*% Z)
      M_l[1:p, (p + 1):len_theta] <- (t(Z) %*% diag(wli * exp(Z_beta_hat))) %*% (B %*% diag(exp(omega_hat)))
      M_l[(p + 1):len_theta, 1:p] <- t(M_l[1:p, (p + 1):len_theta])
      M_l[(p + 1):len_theta, (p + 1):len_theta] <- diag(as.numeric((wli * exp(Z_beta_hat)) %*% B) * exp(omega_hat))
      M_l <- M_l/1 + Lambda
    } else if (basehaz == "unspecified") {
      ## This matrix can be obtain from optim(hessian=T) in optim.survival()
      # Initialize the Hessian matrix
      hessian <- matrix(0, nrow = p, ncol = p)
      for (i in which(status == 1)) {  # Loop over events only
        # Risk set: subjects with time >= time[i]
        risk_set <- which(time >= time[i])
        # Check if risk_set is non-empty
        if (length(risk_set) == 0) next
        # Exponentiated linear predictors for the risk set
        exp_eta <- exp(Z[risk_set, ] %*% as.matrix(beta_hat))
        # Calculate the weighted sum of covariates
        weighted_covariates <- sweep(Z[risk_set, , drop = FALSE], 1, wli[risk_set] * exp_eta, "*")
        # Weighted average of Z
        weighted_cov <- colSums(weighted_covariates) / sum(wli[risk_set] * exp_eta)
        # Weighted covariance matrix of Z in the risk set
        weighted_cov_matrix <- t(Z[risk_set, , drop = FALSE]) %*%
          sweep(Z[risk_set, , drop = FALSE], 1, wli[risk_set] * exp_eta, "*") / sum(wli[risk_set] * exp_eta)
        # Outer product of weighted covariate averages
        outer_product <- weighted_cov %o% weighted_cov # == weighted_cov %*% t(weighted_cov)
        # Update the Hessian matrix
        hessian <- hessian - (weighted_cov_matrix - outer_product)
      }
      if (dim(M_l)[1] == dim(hessian)[1]) {
        M_l <- - hessian + Gamma # not Lambda!
      } else {
        stop("Length of 'theta_hat' should be 'p'.") #!
      }
    } else {
      for (l in 1:len_theta) {
        for (m in 1:len_theta) {
          if (basehaz == "exp") {
            power <- 0
            if (l <= p & m <= p) {
              for (j in 1:nl) {
                power <- power + wli[j] * time[j] * Z[j,][l] * Z[j,][m] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat))
              }
              M_l[l,m] <- power # + Lambda[l,m]
            }
            if (l <= p & m == p+1) {
              for (j in 1:nl) {
                power <- power + wli[j] * time[j] * Z[j,][l] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat))
              }
              M_l[l,m] <- M_l[m,l] <- power # + Lambda[l,m]
            }
            if (l == p+1 & m == p+1) {
              for (j in 1:nl) {
                power <- power + wli[j] * time[j] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat))
              }
              M_l[l,m] <- power # + Lambda[l,m]
            }
            M_l <- M_l/1 + Lambda
          }
          if (basehaz == "gomp") {
            power <- 0
            if (l <= p & m <= p) {
              for (j in 1:nl) {
                power <- power + wli[j] * Z[j,][l] * Z[j,][m] * exp(log_a_hat - log_b_hat +
                                                                      as.numeric(Z[j,] %*% beta_hat)) *
                  (exp(time[j] * exp(log_b_hat)) - 1)
              }
              M_l[l,m] <- power # + Lambda[l,m]
            }
            if (l <= p & m == p+1) {
              for (j in 1:nl) {
                power <- power +  wli[j] * (Z[j,][l] * exp(log_a_hat + time[j] * exp(log_b_hat) +
                                                   as.numeric(Z[j,] %*% beta_hat) - log_b_hat) -
                  Z[j,][l] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) - log_b_hat))
              }
              M_l[l,m] <- M_l[m,l] <- power # + Lambda[l,m]
            }
            if (l <= p & m == p+2) {
              for (j in 1:nl) {
                power <- power +  wli[j] * (Z[j,][l] * (time[j] * exp(log_b_hat) -1) *
                  exp(log_a_hat + time[j] * exp(log_b_hat) + as.numeric(Z[j,] %*% beta_hat) - log_b_hat) +
                  Z[j,][l] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) - log_b_hat))
              }
              M_l[l,m] <- M_l[m,l] <- power # + Lambda[l,m]
            }
            if (l == p+1 & m == p+2) {
              for (j in 1:nl) {
                power <- power +  wli[j] * ((time[j] * exp(log_b_hat) -1) *
                  exp(log_a_hat + time[j] * exp(log_b_hat) + as.numeric(Z[j,] %*% beta_hat) - log_b_hat) +
                  exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) - log_b_hat))
              }
              M_l[l,m] <- M_l[m,l] <- power # + Lambda[l,m]
            }
            if (l == p+1 & m == p+1) {
              for (j in 1:nl) {
                power <- power + wli[j] * (exp(log_a_hat + time[j] * exp(log_b_hat) + as.numeric(Z[j,] %*% beta_hat)
                                     - log_b_hat) - exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) - log_b_hat))
              }
              M_l[l,m] <- power # + Lambda[l,m]
            }
            if (l == p+2 & m == p+2) {
              for (j in 1:nl) {
                power <- power + wli[j] * (((time[j] * exp(log_b_hat) )^2 - time[j] * exp(log_b_hat) + 1) *
                  exp(log_a_hat + time[j] * exp(log_b_hat) + as.numeric(Z[j,] %*% beta_hat) - log_b_hat) -
                  status[j] * time[j] * exp(log_b_hat) - exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) - log_b_hat))
              }
              M_l[l,m] <- power # + Lambda[l,m]
            }
            M_l <- M_l/1 + Lambda
          }
          if (basehaz == "weibul") {
            power <- 0
            if (l <= p & m <= p) {
              for (j in 1:nl) {
                power <- power + wli[j] * Z[j,][l] * Z[j,][m] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) +
                                                             exp(log_b_hat) * log(time[j]))
              }
              M_l[l,m] <- power # + Lambda[l,m]
            }
            if (l <= p & m == p+1) {
              for (j in 1:nl) {
                power <- power + wli[j] * Z[j,][l] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) +
                                                  exp(log_b_hat) * log(time[j]))
              }
              M_l[l,m] <- M_l[m,l] <- power # + Lambda[l,m]
            }
            if (l <= p & m == p+2) {
              for (j in 1:nl) {
                power <- power + wli[j] * Z[j,][l] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) + log_b_hat +
                                                  exp(log_b_hat) * log(time[j])) * log(time[j])
              }
              M_l[l,m] <- M_l[m,l] <- power # + Lambda[l,m]
            }
            if (l == p+1 & m == p+2) {
              for (j in 1:nl) {
                power <- power + wli[j] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) + log_b_hat +
                                        exp(log_b_hat) * log(time[j])) * log(time[j])
              }
              M_l[l,m] <- M_l[m,l] <- power # + Lambda[l,m]
            }
            if (l == p+1 & m == p+1) {
              for (j in 1:nl) {
                power <- power + wli[j] * exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) +
                                       exp(log_b_hat) * log(time[j]))
              }
              M_l[l,m] <- power # + Lambda[l,m]
            }
            if (l == p+2 & m == p+2) {
              for (j in 1:nl) {
                power <- power + wli[j] * (exp(log_a_hat + as.numeric(Z[j,] %*% beta_hat) + log_b_hat +
                                       exp(log_b_hat) * log(time[j])) * log(time[j]) *
                  (1 + exp(log_b_hat) * log(time[j])) - exp(log_b_hat) * log(time[j]) * status[j])
              }
              M_l[l,m] <- power # + Lambda[l,m]
            }
            M_l <- M_l/1 + Lambda
          }
          if (basehaz == "poly") {
            M_k1m2_s_lambda_Taylor <- function(s, m, p, q_l, beta_dotomega){
              log_lam_0l_omeg_l <- 0
              for (ql in 0:q_l) {
                log_lam_0l_omeg_l <- log_lam_0l_omeg_l + (beta_dotomega[p+ql+1]) * s^(ql)
              }
              return(s^(m-p-1) * exp(log_lam_0l_omeg_l))
            }
            M_k2m2_s_lambda_Taylor <- function(s, l, m, p, q_l, beta_dotomega){
              log_lam_0l_omeg_l <- 0
              for (ql in 0:q_l) {
                log_lam_0l_omeg_l <- log_lam_0l_omeg_l + (beta_dotomega[p+ql+1]) * s^(ql)
              }
              return(s^(l+m-2*p-2) * exp(log_lam_0l_omeg_l))
            }
            if (l <= p & m <= p) {
              Lam_0l_omeg_l <- sapply(time, FUN = function(x) integrate(lambda.poly, lower=0,
                                                                        upper=x, p=p, q_l=q_l,
                                                                        beta_dotomega=theta_hat)$value)
              power <- as.numeric(crossprod(Z[,l] * Z[,m] , wli * Lam_0l_omeg_l * exp(as.numeric(Z %*% beta_hat))))
              M_l[l,m] <- power # + Gamma[l,m]
            }
            if (l <= p & m >= p+1) {
              s_k1m2_lam_integ <- sapply(
                time, FUN = function(x) integrate(M_k1m2_s_lambda_Taylor, lower=0,
                                                  upper=x, m=m, p=p, q_l=q_l,
                                                  beta_dotomega=theta_hat)$value)
              power <- as.numeric(crossprod(Z[,l] * s_k1m2_lam_integ, wli * exp(as.numeric(Z %*% beta_hat))))
              M_l[l,m] <- M_l[m,l] <- power # + Lambda[l,m]
            }
            if (l >= p+1 & m >= p+1) {
              s_k2m2_lam_integ <- sapply(
                time, FUN = function(x) integrate(M_k2m2_s_lambda_Taylor, lower=0,
                                                  upper=x, l=l, m=m, p=p, q_l=q_l,
                                                  beta_dotomega=theta_hat)$value)
              power <- as.numeric(crossprod(s_k2m2_lam_integ, wli * exp(as.numeric(Z %*% beta_hat))))
              M_l[l,m] <- power # + Lambda[l,m]  # or Gamma_omega[l-p,m-p]
            }
            M_l <- M_l/1 + Lambda
          }
        }
      }
    }
    if (basehaz == "poly") {
      colnames(M_l) <- rownames(M_l) <- c(colnames(Z), paste("omega",c(0:(len_theta-p-1)), sep="_"))
    } else if (basehaz == "unspecified") {
      colnames(M_l) <- rownames(M_l) <- colnames(Z)
    } else {
      colnames(M_l) <- rownames(M_l) <- c(colnames(Z), paste("omega",c(1:(len_theta-p)), sep="_"))
    }
    A_l <- M_l
  }
  return(A_l)
}
