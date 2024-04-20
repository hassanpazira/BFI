## This file created by Hassan Pazira

#' @export

negloglik.theta <- function(beta_omega, y, X, Lambda, family,
                            q_l, tps, bas, ibas, basehaz, wli) {
  if (is.data.frame(X)) X <- data.matrix(X) # if X is data.frame
  nl <- nrow(X)
  p <- ncol(X)
  beta <- beta_omega[seq_len(p)]
  Gamma <- Lambda[1:p,1:p]
  # if (!(family == "survival" & basehaz == "unspecified")) Gamma_omeg <- Lambda[-(1:p),-(1:p)]
  if (family == "binomial") {
    power <- 0
    for (j in seq_len(nl)) {
      power <- power + wli[j] * (y[j] * as.numeric(X[j, ] %*% beta) -
        log(1 + exp(as.numeric(X[j, ] %*% beta))))
    }
    negloglik <- -power + (t(beta) %*% Gamma %*% beta) / 2
  }

  if (family == "gaussian") {
    Gamma_omeg <- Lambda[-(1:p),-(1:p)]
    log_sig2 <- beta_omega[-seq_len(p)] # log_sig2=log(sigma2) or sigma2=exp(log_sig2)
    power <- 0
    for (j in seq_len(nl)) {
      power <- power + wli[j] * (y[j] - as.numeric(X[j, ] %*% beta))^2
    }
    negloglik <- power / exp(log_sig2) + nl * log(exp(log_sig2)) +
      exp(log_sig2) * as.numeric(Gamma_omeg) + t(beta) %*% Gamma %*% beta
  }
  if (family == "survival") {
    if (colnames(y)[1]=="time" & colnames(y)[2]=="status") {
      time <- y$time
      status <- y$status
    } else {
      time <- y[,1]
      status <- y[,2]
    }
    if (basehaz == "poly") {
      if (length(beta_omega) != (p+q_l+1)) {
        stop("The length of 'theta' should be 'p+q_l+1' = ",sQuote(p+q_l+1),".")
      }
    }
    Z_beta <- as.numeric(X %*% beta)
    if (basehaz == "exp") {
      Gamma_omeg <- Lambda[-(1:p),-(1:p)]
      if (length(beta_omega) != (p+1))
        stop("length of initial vector should be 'p+1' = ", sQuote(p+1), ".")
      a <- (beta_omega[p+1]) # is omega_a ><0
      sig2_loga <- solve(Gamma_omeg)[1, 1]
      power <- sum(wli * ((- status) * a - status * Z_beta + time * exp(a + Z_beta)))
      negloglik  <- power + (t(beta_omega) %*% Lambda %*% beta_omega)/2
    }
    if (basehaz == "gomp") {
      Gamma_omeg <- Lambda[-(1:p),-(1:p)]
      if (length(beta_omega) != (p+2))
        stop("length of initial vector should be 'p+2' = ", sQuote(p+2), ".")
      a <- (beta_omega[p+1]) # is omega_a ><0
      b <- (beta_omega[p+2]) # is omega_b ><0
      sig2_loga <- solve(Gamma_omeg)[1, 1]
      sig2_logb <- solve(Gamma_omeg)[2, 2]
      power <- sum(wli * ((status) * a + status * time * exp(b) + status * Z_beta -
                     exp(a - b + time * exp(b) + Z_beta) + exp(a - b + Z_beta)))
      negloglik  <- - power + (t(beta_omega) %*% Lambda %*% beta_omega)/2
    }
    if (basehaz == "weibul") {
      Gamma_omeg <- Lambda[-(1:p),-(1:p)]
      if (length(beta_omega) != (p+2))
        stop("length of initial vector should be 'p+2' = ", sQuote(p+2), ".")
      a <- (beta_omega[p+1]) # is omega_a ><0
      b <- (beta_omega[p+2]) # is omega_b ><0
      sig2_loga <- solve(Gamma_omeg)[1, 1]
      sig2_logb <- solve(Gamma_omeg)[2, 2]
      power <- sum(wli * (status * (a + b + (exp(b) - 1) * log(time) + Z_beta) -
                     exp(a) * time^(exp(b)) * exp(Z_beta)))
      negloglik  <- - power + (t(beta_omega) %*% Lambda %*% beta_omega)/2
    }
    if (basehaz == "poly") {
      if (length(beta_omega) != (p+q_l+1))
        stop("length of initial vector should be 'p+q_l+1' = ", sQuote(p+q_l+1), ".")
      lam_0l_omeg_l <- lambda.poly(time, p=p, q_l=q_l, beta_dotomega=beta_omega)
      if (any(lam_0l_omeg_l==0)) {
        which_0 <- which(lam_0l_omeg_l == 0)
        lam_0l_omeg_l[which_0] <- .Machine$double.eps
      }
      Lam_0l_omeg_l <- sapply(time, FUN = function(x)
        integrate(lambda.poly, lower=0, upper=x, p=p, q_l=q_l,
                  beta_dotomega=beta_omega)$value)
      power <- sum(wli * (status * (log(lam_0l_omeg_l) + Z_beta) -
                     Lam_0l_omeg_l * exp(Z_beta)))
      negloglik <- - power + (t(beta_omega) %*% Lambda %*% beta_omega)/2
    }
    if (basehaz == "pwexp") {
      omega <- beta_omega[-seq_len(p)]
      if (length(beta_omega) != (p+length(tps)-1))
        stop("length of initial vector should be 'p+length(tps)-1' = ",
             sQuote(p+length(tps)-1), ".")
      nF <- as.numeric((wli * status) %*% bas) #compute array of number of failures
      H <- as.numeric(ibas %*% exp(omega))
      lp <- Z_beta
      negloglik <- ((as.numeric(t(wli * exp(lp)) %*% H) - as.numeric(nF %*% omega) -
                       as.numeric((wli * status) %*% lp))/1 +
                      (t(beta_omega) %*% Lambda %*% beta_omega)/2)
    }
    if (basehaz == "unspecified") {
      if (length(beta_omega) != p)
        stop("length of initial vector should be 'p' = ", sQuote(p), ".")
      log_likelihood <- 0
      for (i in which(status == 1)) {  # Only consider observed events
        eta_i <- Z_beta[i]  # Linear predictor for the event subject i
        # Risk set: subjects with time >= time[i]
        risk_set <- which(time >= time[i])
        # Calculate the sum of exponentials for the risk set
        sum_exp_eta <- sum(wli[risk_set] * exp(Z_beta[risk_set]))
        # Update the log partial likelihood
        log_likelihood <- log_likelihood + wli[i] * (eta_i - log(sum_exp_eta))
      }
      negloglik  <- - log_likelihood + (t(beta) %*% Gamma %*% beta) / 2
    }
  }
  return(as.numeric(negloglik))
}
