## This file created by Hassan Pazira at 16-12-2022
## Updated at 09-07-2023

bfi <- function(theta_hats=NULL, A_hats, Lambda, const_var = NULL, stratified=FALSE, nuisance=1L) {
  if (stratified==T) {
    if (is.null(theta_hats)) stop("In stratified analysis, 'theta_hats' should not be NULL.")
    if (names(theta_hats[[1]])[length(theta_hats[[1]])] == "sigma2") {
      if (!is.numeric(nuisance)) stop("'nuisance' should be a vector of two scalars")
      if (length(nuisance) > 2) stop("For this family the number of nuisance parameters is two, i.e., 'intercept' and 'sigma2'.")
      #n_covars <- n_par - 2
    }
    if (names(theta_hats[[1]])[length(theta_hats[[1]])] != "sigma2") {
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
  L <- length(A_hats)
  if (!is.null(const_var)) {
    if (length(const_var) != L) {
      stop("length of the vector 'const_var' should be equal to ", sQuote(L))
    }
    const_var_fac <- factor(const_var)
    #length(levels(const_var_fac))
    # ...!?
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
  if (is.matrix(Lambda) | (is.list(Lambda) & length(Lambda)==1)) { # all locations have the same Lambda.
    Lambda_all <- list()
    for (i in 1:(L+1)) {
      Lambda_all[[i]] <- Lambda
    }
  } else {
    if (!is.list(Lambda) | (is.list(Lambda) & length(Lambda)!=(L+1))) {
      stop("Lambda should contains of 'L+1' lists; 1:L for local datastes and last one for combined.")
    }
    Lambda_all <- Lambda
  }
  if (stratified == F) {
    A_bfi <- Reduce("+", A_hats) + Lambda_all[[L+1]] - Reduce("+", Lambda_all[1:L])
    sd_bfi <- sqrt(diag(solve(as.matrix(A_bfi))))
    if (!is.null(theta_hats)) {
      theta_hat_bfi <- t(solve(as.matrix(A_bfi)) %*% Reduce("+", Map(function(x,y) x%*%y, A_hats , theta_hats)))
    } else {
      theta_hat_bfi <- theta_hats
    }
    #output <- list(theta_hat=theta_hat_bfi, A_hat=A_bfi, sd=sd_bfi)
  } else {
    for (j in 1:L) {
      if (j==1) {
        A_hats_a_l <- list(A_hats[[j]][-nuisance,-nuisance,drop=F])
        Lambda_a_l <- list(Lambda_all[[j]][-nuisance,-nuisance,drop=F])
        A_hats_b_l <- list(A_hats[[j]][nuisance,nuisance,drop=F])
        Lambda_b_l <- list(Lambda_all[[j]][nuisance,nuisance,drop=F])
        Lambda_b_l[L+1] <- list(Lambda_all[[L+1]][nuisance,nuisance,drop=F])
        A_hats_tilda_b_l <- list(A_hats_b_l[[j]] + Lambda_b_l[[L+1]] - Lambda_b_l[[j]])
        A_hats_tilda_b_l_inv <- list(solve(as.matrix(A_hats_tilda_b_l[[j]])))
        A_hats_ab_l <- list(A_hats[[j]][-nuisance,nuisance,drop=F])
        theta_hat_a_l <- list(theta_hats[[j]][-nuisance])
        theta_hat_b_l <- list(theta_hats[[j]][nuisance])
      } else {
        A_hats_a_l[j] <- list(A_hats[[j]][-nuisance,-nuisance,drop=F])
        Lambda_a_l[j] <- list(Lambda_all[[j]][-nuisance,-nuisance,drop=F])
        if (j==L) Lambda_a_l[j+1] <- list(Lambda_all[[j+1]][-nuisance,-nuisance,drop=F])
        A_hats_b_l[j] <- list(A_hats[[j]][nuisance,nuisance,drop=F])
        Lambda_b_l[j] <- list(Lambda_all[[j]][nuisance,nuisance,drop=F])
        A_hats_tilda_b_l[j] <- list(A_hats_b_l[[j]] + Lambda_b_l[[L+1]] - Lambda_b_l[[j]])
        A_hats_tilda_b_l_inv[j] <- list(solve(as.matrix(A_hats_tilda_b_l[[j]])))
        A_hats_ab_l[j] <- list(A_hats[[j]][-nuisance,nuisance,drop=F])
        theta_hat_a_l[j] <- list(theta_hats[[j]][-nuisance])
        theta_hat_b_l[j] <- list(theta_hats[[j]][nuisance])
      }
    }
    A_a_bfi <- Reduce("+", A_hats_a_l) + Lambda_a_l[[L+1]] - Reduce("+", Lambda_a_l[1:L])
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
    theta_hat_a_bfi <- solve(as.matrix(theta_hat_a_1)) %*% Reduce("+", theta_hat_a_2)
    part1 <- Map("%*%", A_hats_b_l , theta_hat_b_l)
    for (j in 1:L) {
      if (j==1) {
        theta_hat_a_bfi_list <- list(theta_hat_a_bfi)
      } else {
        theta_hat_a_bfi_list[j] <- list(theta_hat_a_bfi)
      }
    }
    theta_hat_a_l_a_bfi <- Map("-", theta_hat_a_l , theta_hat_a_bfi_list)
    part2 <- Map("%*%", lapply(1:L, function(x) t(A_hats_ab_l[[x]])), theta_hat_a_l_a_bfi)
    part12 <- Map("+", part1 , part2)
    theta_hat_tilda_b_l <- Map("%*%", A_hats_tilda_b_l_inv, part12)
    name_nuis <- names(theta_hats[[1]])[nuisance]
    nam_non_nuis <- names(theta_hats[[1]])[-nuisance]
    theta_nuis <- matrix(unlist(theta_hat_tilda_b_l), ncol = length(nuisance), byrow = TRUE)
    colnames(theta_nuis) <- name_nuis
    vec_theta_alls <- new_names <- c()
    kk <- 0
    for (k in 1:n_par) {
      if (k %in% nuisance) {
        which_k <- which(k==nuisance)
        vec_theta_alls <- c(vec_theta_alls, (theta_nuis[,which_k]))
        #fake_name <- name_nuis[which_k]
        pas_nam <- paste0(name_nuis[which_k], rep("_loc", L), 1:L)
        new_names <- c(new_names, pas_nam)
      } else {
        kk <- kk + 1
        new_names <- c(new_names, nam_non_nuis[kk])
        vec_theta_alls <- c(vec_theta_alls, as.numeric(theta_hat_a_bfi[kk]))
      }
    }
    names(vec_theta_alls) <- new_names
    theta_hat_bfi <- vec_theta_alls
    A_bfi <- A_hats
    names(A_bfi) <- paste0("location_", 1:L)
    for (j in 1:L) {
      A_bfi[[j]][nuisance,nuisance] <- A_hats_tilda_b_l[[j]]
      A_bfi[[j]][-nuisance,-nuisance] <- A_a_bfi
    }
    sd_bfi <- list()
    for (j in 1:L) {
      sd_bfi[[j]] <- sqrt(diag(solve(as.matrix(A_bfi[[j]]))))
    }
    names(sd_bfi) <- paste0("location_", 1:L)
    #output <- list(theta_hat_strat=theta_hat_bfi, A_hat_strat=A_bfi, sd_strat=sd_bfi)
  }
  output <- list(theta_hat=theta_hat_bfi, A_hat=A_bfi, sd=sd_bfi)
  return(output)
}
