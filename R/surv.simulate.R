## This file created by Hassan Pazira

#' @export

surv.simulate <- function(L = 1, Z, beta, a, b, u1 = 0, u2 = NULL, cen_rate,
                          gen_data_from = c("exp", "weibul", "gomp"),
                          only_u2 = FALSE, n.rep = 100, Trace = FALSE){
  gen_data_from = match.arg(gen_data_from)
  if (missing(Z)) Z <- NULL
  if (is.list(Z)) {
    if (is.data.frame(Z)) L <- 1
    else L <- length(Z)
  }
  if (!is.null(u2)) cen_rate <- "Not_predefined"
  if (!only_u2) {
    if (is.null(Z)) {
      stop("Z is missing. It should be a list with the length of ", sQuote(L), ".")
    } else {
      if (is.matrix(Z)) Z <- list(Z)
      if (is.data.frame(Z)) {
        L <- 1
        col_names_Z <- colnames(Z)
        Z <- list(as.matrix(Z))
        colnames(Z[[1]]) <- col_names_Z
      }
      if (!is.list(Z)) stop("Z should be a list with the length of ", sQuote(L), ".")
      if (length(Z) != L) stop("The length of Z should be ", sQuote(L), ".")
      for (l in seq_len(L)) {
        if (l == 1) N <- NULL
        N[l] <- nrow(Z[[l]])
        if (is.null(colnames(Z[[l]]))) {
          #colnames(Z[[l]]) <- paste0("Z_",seq_len(ncol(Z[[l]])))
          stop("Colnames of Z cannot be NULL.")
        }
        # Z[[l]] <- Z[[l]][,sort(colnames(Z[[l]]))]
        design_matrix <- paste(colnames(Z[[l]]), collapse=" + ")
        formula <- as.formula(paste(" ", design_matrix, sep=" ~ "))
        X_vars <- model.maker(formula, as.data.frame(Z[[l]]), family="cox")
        Z[[l]] <- X_vars$X[,-1, drop=F]
      }
    }
    p <- ncol(Z[[1]])
    if (length(beta) != p) stop("'beta' should be a vector with ", sQuote(p), " elements!")
  } else {
    p <- length(beta)
  }
  if (missing(cen_rate) & is.null(u2)) stop("'cen_rate' is missing! What is the censoring rate?")
  D <- event_time <- time <- status <- list()
  cen_time <- censor_propor <- NULL
  # if (missing(L)) L <- 1
  # if (missing(u2) | is.null(u2)) u2 <- NULL
  if (a <= 0) {
    if (gen_data_from == "exp") {
      stop("'a' must be > 0. The cumulative hazard form:  a * t , where a > 0.")
    }
    if (gen_data_from == "weibul") {
      stop("'a' must be > 0. The cumulative hazard form:  a * t ^ b , where a > 0 , b > 0.")
    }
    if (gen_data_from == "gomp") {
      stop("'a' must be > 0. The cumulative hazard form:  a * (exp(b * t) - 1) / b , where a > 0 , b > 0.")
    }
  }
  if (gen_data_from != "exp") {
    if (b <= 0) {
      if (gen_data_from == "weibul") {
        stop("'b' must be > 0. The cumulative hazard form:  a * t ^ b , where a > 0 , b > 0.")
      }
      if (gen_data_from == "gomp") {
        stop("'b' must be > 0. The cumulative hazard form:  a * (exp(b * t) - 1) / b , where a > 0 , b > 0.")
      }
    }
  }
  if (cen_rate < 0 & is.null(u2)) stop("'cen_rate' must be >= 0")
  if (cen_rate >= 1 & is.null(u2)) stop("'cen_rate' must be < 1")
  if ((!is.null(u2) | cen_rate==0) & !only_u2) {
    if (cen_rate != 0) {
      # if (missing(u1)) u1 <- 0
      if (u1 < 0) stop("Minimum 'u1' is 0.")
    }
    for (l in 1:L) {
      u <- runif(N[l])
      if (gen_data_from == "exp") {
        event_time[[l]] <- - log(u) * exp(- Z[[l]] %*% (beta)) / a
      }
      if (gen_data_from == "gomp") {
        event_time[[l]] <- log(1 - log(u) * exp(- Z[[l]] %*% (beta)) * b / a) / b
      }
      if (gen_data_from == "weibul") {
        event_time[[l]] <- (- log(u) * exp(- Z[[l]] %*% (beta)) / a)^(1/b)
      }
      if (cen_rate == 0) {
        # follow-up times and event indicators
        time[[l]] <- event_time[[l]]
        status[[l]] <- rep(1, N[l]) # no censoring times or number of events equals to N[l]
        censor_propor[l] <- 0
      } else {
        # censoring times
        cen_time <- runif(N[l], u1, u2) # random censoring C_ij = Uniform(u1, u2)
        # follow-up times and event indicators
        time[[l]] <- pmin(event_time[[l]], cen_time)
        status[[l]] <- as.numeric(event_time[[l]] <= cen_time)
        censor_propor[l] <- (sum(event_time[[l]] > cen_time)/N[l])[1]
      }
      D[[l]] <- data.frame(time = time[[l]], status = status[[l]], Z = Z[[l]])
    }
  } else {
    omega_a <- a
    if (gen_data_from %in% c("weibul", "gomp")) omega_b <- b
    # if (missing(u1)) u1 <- 0  # random censoring C_ij = Uniform(u1, u2)
    if (u1 < 0) stop("Minimum 'u1' is 0.")
    if (gen_data_from %in% c("weibul", "exp")) {
      ## Defining u2 by formula when data is generated from Weibull (or exponential) distribution:
      Inside_integ <- function (s, u1, u2, beta, a, b, sigma2){
        parts <- (pgamma(s*u2^(b), 1/b) - pgamma(s*u1^(b), 1/b)) * gamma(1 + 1/b) *
          dlnorm(s, meanlog = log(a), sdlog = as.numeric(sqrt(sigma2%*%(beta^2)))) / (s^(1/b))
        return(parts)
      }
      cens_prop_fun <- function(x, u1, beta, a, b, sigma2, cen_rate){
        Integ_lambda <- integrate(Inside_integ, 0, Inf, u1=u1, u2=x, beta=beta, a=a, b=b, sigma2=sigma2)$value
        return((1/(x-u1)) * Integ_lambda - cen_rate)
      }
      u2 <- NULL
      if (gen_data_from == "exp") {omega_b2 <- 1}
      if (gen_data_from == "weibul") {omega_b2 <- omega_b}
      u2 <- uniroot(cens_prop_fun, interval=c(u1+0.01, 10^9), u1=u1, beta=beta, a=(omega_a),
                    b=(omega_b2), sigma2=rep(1, p), cen_rate=cen_rate)$root # if b=log(omega_b2) then omega_b2=1
      #if (Trace) cat("Desired u2 is:", u2, "\n")
    } else {
      # gomp (need to be updated!)
      u2 <- 20  # u2 is arbitrary!
      second_time <- last_time <- FALSE
      repeat {
        cp_mat <- NULL
        for (ii in 1:n.rep) {
          cp_l <- NULL
          dat_u2 <- surv.simulate(L=L, Z=Z, beta=beta, a=omega_a, b=omega_b, u1=u1,
                                   u2=u2, cen_rate=cen_rate, gen_data_from=gen_data_from)
          for (ll in 1:L) {
            cp_l <- c(cp_l, dat_u2$censor_propor[[ll]])
          }
          cp_mat <- rbind(cp_mat, cp_l)
        }
        means <- apply(cp_mat, 2, mean)
        if (Trace) cat("censoring proportion for silos are: ", round(means,2),
                       "  with u2 =", u2, "\n")
        means_dif <- mean(means)-cen_rate
        if (abs(means_dif) > 0.01) {
          if (means_dif > 0.01) {
            u2 <- u2 + abs(means_dif)*u2
          } else {
            u2 <- u2 - abs(means_dif)*u2
          }
        } else {
          u2 <- u2
          if (Trace) cat("\n", "Desired u2 is:", u2, "( with u1 =", u1,")", "\n \n")
          break
        }
        if (u2-u1<0.01) {
          u2<-30
          u1<-.5
          if (second_time==TRUE) {
            u1<-.01
            if (last_time==TRUE) stop("Decrease censoring rate!")
          }
          second_time <- TRUE
        }
      }
    }
    if (!only_u2) {
      for (l in 1:L) {
        u <- runif(N[l])
        if (gen_data_from == "exp") {
          event_time[[l]] <- - log(u) * exp(- Z[[l]] %*% (beta)) / a
        }
        if (gen_data_from == "gomp") {
          event_time[[l]] <- log(1 - log(u) * exp(- Z[[l]] %*% (beta)) * b / a) / b
        }
        if (gen_data_from == "weibul") {
          event_time[[l]] <- (- log(u) * exp(- Z[[l]] %*% (beta)) / a)^(1/b)
        }
        if (cen_rate == 0) {
          # follow-up times and event indicators
          time[[l]] <- event_time[[l]]
          status[[l]] <- rep(1, N[l]) # no censoring times or number of events equals to N[l]
          censor_propor[l] <- 0
        } else {
          # censoring times
          cen_time <- runif(N[l], u1, u2)
          # follow-up times and event indicators
          time[[l]] <- pmin(event_time[[l]], cen_time)
          status[[l]] <- as.numeric(event_time[[l]] <= cen_time)
          censor_propor[l] <- (sum(event_time[[l]] > cen_time)/N[l])[1]
        }
        D[[l]] <- data.frame(time = time[[l]], status = status[[l]], Z = Z[[l]])
      }
    }
  }
  if (cen_rate == 0) u1 <- NULL
  if (only_u2) {
    return(list(D=NULL, censor_propor=NULL, u1=u1, u2=u2))
  } else {
    return(list(D=D, censor_propor=censor_propor, u1=u1, u2=u2))
  }
}

