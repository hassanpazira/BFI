## This file created by Hassan Pazira at 16-12-2022
## Updated at 22-06-2023
bfi <- function(theta_hats=NULL, A_hats, Gamma, L=NULL, stratified=F, nuisance=1L, family=c("binomial","gaussian")) {
  if (stratified==T) {
    if (is.null(theta_hats)) stop("In stratified analysis, 'theta_hats' should not be NULL.")
    family <- match.arg(family)
    if (!family %in% c("binomial", "gaussian")) {
      stop("Distributions that can be used are 'binomial' and 'gaussian' in this version of the package!")
    }
    if (family=="gaussian") {
      if (!is.numeric(nuisance)) stop("'nuisance' should be a vector of two elements")
      if (length(nuisance) > 2) stop("For this family the number of nuisance parameters is two, i.e., 'intercept' and 'sigma2_e'.")
      #n_covars <- n_par - 2
    }
    if (family=="binomial") {
      if (!is.numeric(nuisance)) stop("'nuisance' should be an scalar.")
      if (length(nuisance) > 1) stop("For this family the number of nuisance parameters is 1, i.e., 'intercept'.")
      #n_covars <- n_par - 1
    }
  }
  if (!is.null(theta_hats) & !is.list(theta_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.list(A_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.null(theta_hats) & (length(theta_hats) != length(A_hats))) {
    stop("Length of inputs are not equal.")
  }
  if (is.null(L)) {
    L <- length(A_hats)
  } else {
    if (length(A_hats) != L) {
      stop("Number of locations and length of inputs should be the same.")
    }
  }
  for (i in 1:L) {
    if (i==1) {
      n_pars <- ncol(A_hats[[i]])
    } else {
      n_pars[i] <- ncol(A_hats[[i]])
    }
  }
  if (all(n_pars == n_pars[1])) {
    n_par <- n_pars[1]
  } else stop("All matrix in A_hats should have the same columns. Number of parameters in locations must be equal.")
  if (is.matrix(Gamma) | (is.list(Gamma) & length(Gamma)==1)) { # all locations have the same Gamma.
    Gamma_all <- list()
    for (i in 1:(L+1)) {
      Gamma_all[[i]] <- Gamma
    }
  } else {
    if (!is.list(Gamma) | (is.list(Gamma) & length(Gamma)!=(L+1))) {
      stop("Gamma should contains of L+1 lists; 1:L for local datastes and last one for combined.")
    }
    Gamma_all <- Gamma
  }
  if (stratified == F) {
    A_bfi <- Reduce("+", A_hats) + Gamma_all[[L+1]] - Reduce("+", Gamma_all[1:L])
    sd_bfi <- sqrt(diag(solve(as.matrix(A_bfi))))
    if (!is.null(theta_hats)) {
      theta_hat_bfi <- t(solve(as.matrix(A_bfi)) %*% Reduce("+", Map(function(x,y) x%*%y, A_hats , theta_hats)))
    } else {
      theta_hat_bfi <- t(theta_hats)
    }
  } else {
    old_thetas <- old_A_diags <- c()
    for (j in 1:L) {
      old_thetas <- rbind(old_thetas, theta_hats[[j]][nuisance])
      old_A_diags <- rbind(old_A_diags, diag(A_hats[[j]])[nuisance])
      if (j==1) {
        A_hats_a_l  <- list(A_hats[[j]][-nuisance,-nuisance,drop=F])
        Gamma_a_l   <- list(Gamma_all[[j]][-nuisance,-nuisance,drop=F])
        A_hats_b_l  <- list(A_hats[[j]][nuisance,nuisance,drop=F])
        Gamma_b_l   <- list(Gamma_all[[j]][nuisance,nuisance,drop=F])
        Gamma_b_l[L+1] <- list(Gamma_all[[L+1]][nuisance,nuisance,drop=F])
        A_hats_tilda_b_l <- list(A_hats_b_l[[j]] + Gamma_b_l[[L+1]] - Gamma_b_l[[j]])
        A_hats_tilda_b_l_inv  <- list(solve(as.matrix(A_hats_tilda_b_l[[j]])))
        A_hats_ab_l <- list(A_hats[[j]][-nuisance,nuisance,drop=F])
        theta_hat_a_l <- list(theta_hats[[j]][-nuisance])
        theta_hat_b_l <- list(theta_hats[[j]][nuisance])
      } else {
        A_hats_a_l[j]  <- list(A_hats[[j]][-nuisance,-nuisance,drop=F])
        Gamma_a_l[j]   <- list(Gamma_all[[j]][-nuisance,-nuisance,drop=F])
        if (j==L) Gamma_a_l[j+1] <- list(Gamma_all[[j+1]][-nuisance,-nuisance,drop=F])
        A_hats_b_l[j]  <- list(A_hats[[j]][nuisance,nuisance,drop=F])
        Gamma_b_l[j]   <- list(Gamma_all[[j]][nuisance,nuisance,drop=F])
        A_hats_tilda_b_l[j]  <- list(A_hats_b_l[[j]] + Gamma_b_l[[L+1]] - Gamma_b_l[[j]])
        A_hats_tilda_b_l_inv[j]  <- list(solve(as.matrix(A_hats_tilda_b_l[[j]])))
        A_hats_ab_l[j] <- list(A_hats[[j]][-nuisance,nuisance,drop=F])
        theta_hat_a_l <- list(theta_hats[[j]][-nuisance])
        theta_hat_b_l <- list(theta_hats[[j]][nuisance])
      }
    }
    A_a_bfi <- Reduce("+", A_hats_a_l) + Gamma_a_l[[L+1]] - Reduce("+", Gamma_a_l[1:L])
    sd_a_bfi <- sqrt(diag(solve(as.matrix(A_a_bfi))))
    A_hat_ab2 <- Map("%*%", A_hats_ab_l , A_hats_tilda_b_l_inv)
    A_hat_ab3 <- Map("%*%", A_hat_ab2, lapply(1:L, function(x) t(A_hats_ab_l[[x]])))
    theta_hat_a_1 <- A_a_bfi - Reduce("+", A_hat_ab3)
    A_hat_a_l_1 <- Map("-", A_hats_a_l , A_hat_ab3)
    I_d2 <- lapply(1:L, function(x) diag(length(nuisance)))
    A_hat_b2 <- Map("%*%", A_hats_tilda_b_l_inv, A_hats_b_l)
    A_hat_a_l_2 <- Map("-", I_d2, A_hat_b2)
    map1 <- Map("%*%", A_hat_a_l_1 , theta_hat_a_l)
    map02 <- Map("%*%", A_hats_ab_l , A_hat_a_l_2)
    map2 <- Map("%*%", map02, theta_hat_b_l)
    theta_hat_a_2 <- Map("+", map1 , map2)
    theta_hat_a_bfi <- solve(as.matrix(theta_hat_a_1)) %*% Reduce("+", theta_hat_a_2) #!!!
    nam_non_nuis <- names(theta_hats[[1]])[-nuisance]
    vec_theta_alls <- vec_A_alls <- new_names <- c()
    kk <- 0
    for (k in 1:n_par) {
      if (k %in% nuisance) {
        which_k <- which(k==nuisance)
        vec_theta_alls <- c(vec_theta_alls, (old_thetas[,which_k]))
        fake_name <- colnames(old_thetas)[which_k]
        pas_nam <- paste(fake_name,1:L)
        new_names <- c(new_names, pas_nam)
        #vec_A_alls <- c(vec_A_alls, as.numeric(old_A_diags[,which_k]))
      } else {
        kk <- kk + 1
        vec_theta_alls <- c(vec_theta_alls, as.numeric(theta_hat_a_bfi[kk]))
        new_names <- c(new_names, nam_non_nuis[kk])

        #vec_A_alls <- c(vec_A_alls, as.numeric(diag(as.matrix(A_a_bfi))[kk]))
      }
    }
    names(vec_theta_alls) <- new_names
    theta_hat_bfi <- vec_theta_alls
    A_bfi <- A_a_bfi # !
    sd_bfi <- sqrt(diag(solve(as.matrix(A_bfi))))
  }
  output <- list(theta_hat=theta_hat_bfi, A_hat=A_bfi, sd=sd_bfi)
  return(output)
}
