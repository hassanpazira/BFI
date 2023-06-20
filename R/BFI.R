## This file created by Hassan Pazira at 16-12-2022
## Updated at 29-09-2023

bfi <- function(theta_hats=NULL, A_hats, Lambda, stratified=FALSE, strat_par=NULL,
                center_spec=NULL, type=c("cat", "con")) {
  if (!is.null(theta_hats) & !is.list(theta_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.list(A_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.null(theta_hats) & (length(theta_hats) != length(A_hats))) {
    stop("Length of inputs are not equal.")
  }
  if (stratified == T & is.null(strat_par) & is.null(center_spec)) {
    stop("Since 'stratified = TRUE', only one of 'strat_par' or 'center_spec' could be NULL.")
  }
  if (stratified == T & !is.null(strat_par) & !is.null(center_spec)) {
    stop("Since 'stratified = TRUE', only one of 'strat_par' or 'center_spec' could be non-NULL.")
  }
  if (stratified == F & (!is.null(strat_par) | !is.null(center_spec))) {
    stop("Since 'stratified = FALSE', both 'strat_par' and 'center_spec' sould be NULL.")
  }
  L <- length(A_hats)
  if (L == 1) stop("he number of locations should be > 1.")
  if (is.data.frame(Lambda)) {
    Lambda <- as.matrix(Lambda)
  }
  if (is.list(Lambda)) {
    col_nam_Lambda <- colnames(Lambda[[1]])
  } else {
    col_nam_Lambda <- colnames(Lambda)
  }
  if (is.null(names(theta_hats[[1]])) & is.null(colnames(A_hats[[1]])) & is.null(col_nam_Lambda)) {
    stop("The elements of the input lists or at least Lambda should be named.")
  } else {
    if (!is.null(names(theta_hats[[1]]))) {
      col_nam <- names(theta_hats[[1]])
    } else {
      if (!is.null(colnames(A_hats[[1]]))) {
        col_nam <- colnames(A_hats[[1]])
      } else {
        col_nam <- col_nam_Lambda
      }
    }
    if (is.null(names(theta_hats[[1]]))) {
      for (i in 1:(L)) {
        names(theta_hats[[i]]) <- col_nam
      }
    }
    if (is.null(colnames(A_hats[[1]]))) {
      for (i in 1:(L)) {
        colnames(A_hats[[i]]) <- rownames(A_hats[[i]]) <- col_nam
      }
    }
    if (is.null(col_nam_Lambda)) {
      if (is.list(Lambda)) {
        for (i in 1:length(Lambda)) {
          colnames(Lambda[[i]]) <- rownames(Lambda[[i]]) <- col_nam
        }
      } else {
        rownames(Lambda) <- colnames(Lambda) <- col_nam
      }
    }
  }
  if (is.matrix(Lambda) | (is.list(Lambda) & length(Lambda)==1)) {
    # all locations and combined data have the same Lambda
    Lambda_all <- list()
    if (is.list(Lambda)) Lambda <- Lambda[[1]]
    for (i in 1:(L+1)) {
      Lambda_all[[i]] <- Lambda
    }
  } else {
    if (is.list(Lambda)) {
      if (length(Lambda) == 2) {
        # all locations have the same Lambda but combined has different Lambda
        Lambda_all <- list()
        for (i in 1:(L)) {
          Lambda_all[[i]] <- Lambda[[1]]
        }
        Lambda_all[[(L+1)]] <- Lambda[[2]]
      } else {
        if (length(Lambda) != (L+1)) {
          stop("Lambda should contain 'L+1' lists; 1:L for local datastes and last one for combined.")
        }
        Lambda_all <- Lambda
      }
    } else {
      stop("Lambda should be a list.")
    }
  }
  if (names(theta_hats[[1]])[length(theta_hats[[1]])] == "sigma2") {
    family <- c("gaussian")
  } else {
    family <- c("binomial")
  }
  if (stratified == T) {
    if (is.null(strat_par))
      stop("In stratified analysis, 'strat_par' should not be 'NULL'. Select the stratification parameter(s).")
    if (! 1 %in% strat_par & ! 2 %in% strat_par) {
      stop("'strat_par' should contain '1' and/or '2'.")
    }
    if (is.null(theta_hats))
      stop("In stratified analysis, 'theta_hats' should not be 'NULL'.")
    if (family == "gaussian") {
      if (!is.numeric(strat_par))
        stop("'strat_par' should be one of the integers: 1 or 2, or a vector of both.")
      if (length(strat_par) > 2)
        stop("For this family ('gaussian') the number of 'strat_par' parameters is two, i.e., 'intercept' and 'sigma2'.")
      if (names(theta_hats[[1]])[1] != "(Intercept)") {
        if (length(strat_par) > 1 | 1 %in% strat_par) {
          stop("Since 'intercept = F', for this family ('gaussian') 'strat_par' should only be the integer 2")
        }
      }
    }
    if (family == "binomial") {
      if (names(theta_hats[[1]])[1] != "(Intercept)") {
        stop("Since 'intercept = F', for this family ('binomial') the stratified analysis is not possible!")
      }
      if (!is.numeric(strat_par))
        stop("'strat_par' should be an integer")
      if (length(strat_par) > 1)
        stop("For this family ('binomial') the number of 'strat_par' parameters is one, i.e., 'intercept'.")
      if (strat_par == 2)
        stop("'strat_par' should only be the integer 1")
    }
  }
  if (stratified==T & is.null(strat_par) & !is.null(center_spec)) {
    if (length(center_spec) != L)
      stop("length of the vector 'center_spec' should be equal to ", sQuote(L))
    center_spec_fact <- factor(center_spec)
    if (length(levels(center_spec_fact)) < 2)
      stop("the number of levels should be >= 2 ")
    #
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
  } else {
    stop("All matrix in A_hats should have the same columns. Number of parameters in locations must be equal.")
  }
  if (stratified == F) {
    if (dim(Lambda_all[[L+1]])[1] != dim(Lambda_all[[1]])[1])
      stop("The last matrix in 'Lambda' should have the same dim. as the other local matrices.")
    A_bfi <- Reduce("+", A_hats) + Lambda_all[[L+1]] - Reduce("+", Lambda_all[1:L])
    sd_bfi <- sqrt(diag(solve(as.matrix(A_bfi))))
    if (!is.null(theta_hats)) {
      theta_hat_bfi <- t(solve(as.matrix(A_bfi)) %*% Reduce("+", Map(function(x,y) x%*%y, A_hats , theta_hats)))
    } else {
      theta_hat_bfi <- theta_hats
    }
    #output <- list(theta_hat=theta_hat_bfi, A_hat=A_bfi, sd=sd_bfi)
  } else {
    if (dim(Lambda_all[[L+1]])[1] == dim(Lambda_all[[1]])[1])
      stop("The last matrix in 'Lambda' should not have the same dim. as the other local matrices.")
    if (length(strat_par) == 1) {
      if (1 %in% strat_par) {
        strat_par <- c(1)
      } else {
        strat_par <- c(length(theta_hats[[1]]))
      }
    } else {
      strat_par <- c(1, length(theta_hats[[1]]))
    }
    for (j in 1:L) {
      if (j==1) {
        A_hats_a_l <- list(A_hats[[j]][-strat_par,-strat_par,drop=F])
        Lambda_a_l <- list(Lambda_all[[j]][-strat_par,-strat_par,drop=F])
        A_hats_b_l <- list(A_hats[[j]][strat_par,strat_par,drop=F])
        Lambda_b_l <- list(Lambda_all[[j]][strat_par,strat_par,drop=F])
        Lambda_b_l[L+1] <- list(Lambda_all[[L+1]][strat_par,strat_par,drop=F])
        A_hats_tilda_b_l <- list(A_hats_b_l[[j]] + Lambda_b_l[[L+1]] - Lambda_b_l[[j]])
        A_hats_tilda_b_l_inv <- list(solve(as.matrix(A_hats_tilda_b_l[[j]])))
        A_hats_ab_l <- list(A_hats[[j]][-strat_par,strat_par,drop=F])
        theta_hat_a_l <- list(theta_hats[[j]][-strat_par])
        theta_hat_b_l <- list(theta_hats[[j]][strat_par])
      } else {
        A_hats_a_l[j] <- list(A_hats[[j]][-strat_par,-strat_par,drop=F])
        Lambda_a_l[j] <- list(Lambda_all[[j]][-strat_par,-strat_par,drop=F])
        if (j==L) Lambda_a_l[j+1] <- list(Lambda_all[[j+1]][-strat_par,-strat_par,drop=F])
        A_hats_b_l[j] <- list(A_hats[[j]][strat_par,strat_par,drop=F])
        Lambda_b_l[j] <- list(Lambda_all[[j]][strat_par,strat_par,drop=F])
        A_hats_tilda_b_l[j] <- list(A_hats_b_l[[j]] + Lambda_b_l[[L+1]] - Lambda_b_l[[j]])
        A_hats_tilda_b_l_inv[j] <- list(solve(as.matrix(A_hats_tilda_b_l[[j]])))
        A_hats_ab_l[j] <- list(A_hats[[j]][-strat_par,strat_par,drop=F])
        theta_hat_a_l[j] <- list(theta_hats[[j]][-strat_par])
        theta_hat_b_l[j] <- list(theta_hats[[j]][strat_par])
      }
    }
    A_a_bfi <- Reduce("+", A_hats_a_l) + Lambda_a_l[[L+1]] - Reduce("+", Lambda_a_l[1:L])
    A_hat_ab2 <- Map("%*%", A_hats_ab_l , A_hats_tilda_b_l_inv)
    A_hat_ab3 <- Map("%*%", A_hat_ab2, lapply(1:L, function(x) t(A_hats_ab_l[[x]])))
    theta_hat_a_1 <- A_a_bfi - Reduce("+", A_hat_ab3)
    A_hat_a_l_1 <- Map("-", A_hats_a_l , A_hat_ab3)
    I_d2 <- lapply(1:L, function(x) diag(length(strat_par)))
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
    name_nuis <- names(theta_hats[[1]])[strat_par]
    nam_non_nuis <- names(theta_hats[[1]])[-strat_par]
    theta_nuis <- matrix(unlist(theta_hat_tilda_b_l), ncol = length(strat_par), byrow = TRUE)
    colnames(theta_nuis) <- name_nuis
    vec_theta_alls <- new_names <- c()
    kk <- 0
    for (k in 1:n_par) {
      if (k %in% strat_par) {
        which_k <- which(k==strat_par)
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
      A_bfi[[j]][strat_par,strat_par] <- A_hats_tilda_b_l[[j]]
      A_bfi[[j]][-strat_par,-strat_par] <- A_a_bfi
    }
    sd_bfi <- list()
    for (j in 1:L) {
      sd_bfi[[j]] <- sqrt(diag(solve(as.matrix(A_bfi[[j]]))))
    }
    names(sd_bfi) <- paste0("location_", 1:L)
    #output <- list(theta_hat_strat=theta_hat_bfi, A_hat_strat=A_bfi, sd_strat=sd_bfi)
  }
  output <- list(theta_hat=theta_hat_bfi, A_hat=A_bfi, sd=sd_bfi, family=family, stratified=stratified)
  class(output) <- "bfi"
  return(output)
}
