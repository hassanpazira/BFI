## This file created by Hassan Pazira at 16-12-2022
## Updated at 29-09-2023

#' @export

bfi <- function(theta_hats = NULL, A_hats, Lambda, stratified = FALSE, strat_par = NULL, center_spec = NULL) {
  if (!is.null(theta_hats) & !is.list(theta_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.list(A_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.null(theta_hats) & (length(theta_hats) != length(A_hats))) {
    stop("Length of inputs are not equal.")
  }
  if (is.null(theta_hats)) {
    stop("argument 'theta_hats' is missing, with no default.")
  }
  if (stratified == TRUE & is.null(strat_par) & is.null(center_spec)) {
    stop("Since 'stratified = TRUE', only one of 'strat_par' or 'center_spec'
         could be NULL.")
  }
  if (stratified == TRUE & !is.null(strat_par) & !is.null(center_spec)) {
    stop("Since 'stratified = TRUE', only one of 'strat_par' or 'center_spec'
         could be non-NULL.")
  }
  if (stratified == FALSE & (!is.null(strat_par) | !is.null(center_spec))) {
    stop("Since 'stratified = FALSE', both 'strat_par' and 'center_spec'
         sould be NULL.")
  }
  L <- length(A_hats)
  if (L == 1) stop("The number of locations should be > 1.")
  if (is.data.frame(Lambda)) {
    Lambda <- as.matrix(Lambda)
  }
  if (is.list(Lambda)) {
    col_nam_Lambda <- colnames(Lambda[[1]])
  } else {
    col_nam_Lambda <- colnames(Lambda)
  }
  if (is.null(names(theta_hats[[1]])) & is.null(colnames(A_hats[[1]])) &
      is.null(col_nam_Lambda)) {
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
      for (i in seq_len(L)) {
        names(theta_hats[[i]]) <- col_nam
      }
    }
    if (is.null(colnames(A_hats[[1]]))) {
      for (i in seq_len(L)) {
        colnames(A_hats[[i]]) <- rownames(A_hats[[i]]) <- col_nam
      }
    }
    if (is.null(col_nam_Lambda)) {
      if (is.list(Lambda)) {
        for (i in seq_len(length(Lambda))) {
          colnames(Lambda[[i]]) <- rownames(Lambda[[i]]) <- col_nam
        }
      } else {
        rownames(Lambda) <- colnames(Lambda) <- col_nam
      }
    }
  }
  if (is.matrix(Lambda) | (is.list(Lambda) & length(Lambda) == 1)) {
    # all locations and combined data have the same Lambda
    Lambda_all <- list()
    if (is.list(Lambda)) Lambda <- Lambda[[1]]
    for (i in seq_len(L + 1)) {
      Lambda_all[[i]] <- Lambda
    }
  } else {
    if (is.list(Lambda)) {
      if (length(Lambda) == 2) {
        # all locations have the same Lambda but combined has different Lambda
        Lambda_all <- list()
        for (i in seq_len(L)) {
          Lambda_all[[i]] <- Lambda[[1]]
        }
        Lambda_all[[(L + 1)]] <- Lambda[[2]]
      } else {
        if (length(Lambda) != (L + 1)) {
          stop("Lambda should contain 'L+1' lists; 1:L for local datastes and
               last one for combined.")
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
  if (stratified == TRUE & !is.null(strat_par)) {
    if ((!1 %in% strat_par) & (!2 %in% strat_par)) {
      stop("'strat_par' should contain '1' and/or '2'.")
    }
    if (is.null(theta_hats)) {
      stop("In stratified analysis, 'theta_hats' should not be 'NULL'.")
    }
    if (family == "gaussian") {
      if (!is.numeric(strat_par)) {
        stop("'strat_par' should be one of the integers: 1 or 2, or a vector of both.")
      }
      if (length(strat_par) > 2) {
        stop("For this family ('gaussian') the number of 'strat_par' parameters
             is two, i.e., 'intercept' and 'sigma2'.")
      }
      if (names(theta_hats[[1]])[1] != "(Intercept)") {
        if (length(strat_par) > 1 | 1 %in% strat_par) {
          stop("Since 'intercept = FALSE', for this family ('gaussian')
               'strat_par' should only be the integer 2")
        }
      }
    }
    if (family == "binomial") {
      if (names(theta_hats[[1]])[1] != "(Intercept)") {
        stop("Since 'intercept = FALSE', for this family ('binomial')
             the stratified analysis is not possible!")
      }
      if (!is.numeric(strat_par)) {
        stop("'strat_par' should be an integer")
      }
      if (length(strat_par) > 1) {
        stop("For this family ('binomial') the number of 'strat_par' parameters
             is one, i.e., 'intercept'.")
      }
      if (strat_par == 2) {
        stop("'strat_par' should only be the integer 1")
      }
    }
  }
  for (i in seq_len(L)) {
    if (i == 1) {
      n_pars <- ncol(A_hats[[i]])
    } else {
      n_pars[i] <- ncol(A_hats[[i]])
    }
  }
  if (all(n_pars == n_pars[1])) {
    n_par <- n_pars[1]
  } else {
    stop("All matrix in A_hats should have the same columns.
         Number of parameters in locations must be equal.")
  }
  last_Lam_dim <- dim(Lambda_all[[L + 1]])[1]
  if (stratified == FALSE) {
    if (last_Lam_dim != dim(Lambda_all[[1]])[1]) {
      stop("The last matrix in 'Lambda' should have the same dim. as the other
           local matrices.")
    }
    A_bfi <- Reduce("+", A_hats) + Lambda_all[[L + 1]] -
      Reduce("+", Lambda_all[seq_len(L)])
    sd_bfi <- sqrt(diag(solve(as.matrix(A_bfi))))
    if (!is.null(theta_hats)) {
      theta_hat_bfi <- t(solve(as.matrix(A_bfi)) %*%
                           Reduce("+", Map(function(x, y) x %*% y,
                                           A_hats, theta_hats)))
    } else {
      theta_hat_bfi <- theta_hats
    }
    # output <- list(theta_hat=theta_hat_bfi, A_hat=A_bfi, sd=sd_bfi)
  } else {
    if (last_Lam_dim == dim(Lambda_all[[1]])[1]) {
      if (all(colnames(Lambda_all[[1]]) == colnames(Lambda_all[[L + 1]]))) {
        stop("The last matrix in 'Lambda' should not have the same dim. as the
             other local matrices.")
      } else {
        if (length(strat_par) == 1) {
          if (1 %in% strat_par) {
            if (colnames(Lambda_all[[1]])[ncol(Lambda_all[[1]])] !=
                colnames(Lambda_all[[L + 1]])[ncol(Lambda_all[[L + 1]])]) {
              stop("Check the elements of 'Lambda'. Is there an 'Sigma2'
                   in the model?")
            }
          } else {
            if (colnames(Lambda_all[[1]])[1] != colnames(Lambda_all[[L + 1]])[1]) {
              stop("Check the elements of 'Lambda'.
                   Is there an 'intercept' in the model?")
            }
          }
        }
      }
    }
    if (is.null(center_spec)) {
      if (length(strat_par) == 1) {
        if (1 %in% strat_par) {
          strat_par <- c(1)
          noncore <- cbind(seq_len(L))
          new_noncore <- noncore
        } else {
          strat_par <- c(length(theta_hats[[1]]))
          noncore <- cbind((last_Lam_dim - L + 1):last_Lam_dim)
          new_noncore <- as.matrix(seq_len(L)) # noncore - (L+1) #
        }
      } else {
        strat_par <- c(1, length(theta_hats[[1]]))
        noncore <- new_noncore <- cbind(seq_len(L),
                                        (last_Lam_dim - L + 1):last_Lam_dim)
        new_noncore[, 2] <- noncore[, 1] + L
      }
      for (j in seq_len(L)) {
        if (j == 1) {
          A_hats_a_l <- list(A_hats[[j]][-strat_par, -strat_par, drop = FALSE ])
          Lambda_a_l <- list(Lambda_all[[j]][-strat_par, -strat_par, drop = FALSE ])
          Lambda_a_l[L + 1] <- list(Lambda_all[[L + 1]][-as.numeric(noncore),
                                                        -as.numeric(noncore), drop = FALSE ])
          A_hats_b_l <- list(A_hats[[j]][strat_par, strat_par, drop = FALSE ])
          Lambda_b_l <- list(Lambda_all[[j]][strat_par, strat_par, drop = FALSE ])
          Lambda_b_l[L + 1] <- list(Lambda_all[[L + 1]][as.numeric(noncore),
                                                        as.numeric(noncore), drop = FALSE ])
          A_hats_tilda_b_l <- list(A_hats_b_l[[j]] +
                                     Lambda_b_l[[L + 1]][new_noncore[j, ],
                                                         new_noncore[j, ], drop = FALSE ] -
            Lambda_b_l[[j]])
          A_hats_tilda_b_l_inv <- list(solve(as.matrix(A_hats_tilda_b_l[[j]])))
          A_hats_ab_l <- list(A_hats[[j]][-strat_par, strat_par, drop = FALSE ])
          theta_hat_a_l <- list(theta_hats[[j]][-strat_par])
          theta_hat_b_l <- list(theta_hats[[j]][strat_par])
        } else {
          A_hats_a_l[j] <- list(A_hats[[j]][-strat_par, -strat_par, drop = FALSE ])
          Lambda_a_l[j] <- list(Lambda_all[[j]][-strat_par, -strat_par, drop = FALSE ])
          A_hats_b_l[j] <- list(A_hats[[j]][strat_par, strat_par, drop = FALSE ])
          Lambda_b_l[j] <- list(Lambda_all[[j]][strat_par, strat_par, drop = FALSE ])
          A_hats_tilda_b_l[j] <- list(A_hats_b_l[[j]] +
                                        Lambda_b_l[[L + 1]][new_noncore[j, ],
                                                            new_noncore[j, ], drop = FALSE ] -
                                        Lambda_b_l[[j]])
          A_hats_tilda_b_l_inv[j] <- list(solve(as.matrix(A_hats_tilda_b_l[[j]])))
          A_hats_ab_l[j] <- list(A_hats[[j]][-strat_par, strat_par, drop = FALSE ])
          theta_hat_a_l[j] <- list(theta_hats[[j]][-strat_par])
          theta_hat_b_l[j] <- list(theta_hats[[j]][strat_par])
        }
      }
      A_a_bfi <- Reduce("+", A_hats_a_l) + Lambda_a_l[[L + 1]] -
        Reduce("+", Lambda_a_l[seq_len(L)])
      A_hat_ab2 <- Map("%*%", A_hats_ab_l, A_hats_tilda_b_l_inv)
      A_hat_ab3 <- Map("%*%", A_hat_ab2,
                       lapply(seq_len(L), function(x) t(A_hats_ab_l[[x]])))
      theta_hat_a_1 <- A_a_bfi - Reduce("+", A_hat_ab3)
      A_hat_a_l_1 <- Map("-", A_hats_a_l, A_hat_ab3)
      I_d2 <- lapply(seq_len(L), function(x) diag(length(strat_par)))
      A_hat_b2 <- Map("%*%", A_hats_tilda_b_l_inv, A_hats_b_l)
      A_hat_a_l_2 <- Map("-", I_d2, A_hat_b2)
      map1 <- Map("%*%", A_hat_a_l_1, theta_hat_a_l)
      map02 <- Map("%*%", A_hats_ab_l, A_hat_a_l_2)
      map2 <- Map("%*%", map02, theta_hat_b_l)
      theta_hat_a_2 <- Map("+", map1, map2)
      theta_hat_a_bfi <- solve(as.matrix(theta_hat_a_1)) %*%
        Reduce("+", theta_hat_a_2)
      part1 <- Map("%*%", A_hats_b_l, theta_hat_b_l)
      for (j in seq_len(L)) {
        if (j == 1) {
          theta_hat_a_bfi_list <- list(theta_hat_a_bfi)
        } else {
          theta_hat_a_bfi_list[j] <- list(theta_hat_a_bfi)
        }
      }
      theta_hat_a_l_a_bfi <- Map("-", theta_hat_a_l, theta_hat_a_bfi_list)
      part2 <- Map("%*%", lapply(seq_len(L), function(x)
        t(A_hats_ab_l[[x]])), theta_hat_a_l_a_bfi)
      part12 <- Map("+", part1, part2)
      theta_hat_tilda_b_l <- Map("%*%", A_hats_tilda_b_l_inv, part12)
      name_nuis <- names(theta_hats[[1]])[strat_par]
      nam_non_nuis <- names(theta_hats[[1]])[-strat_par]
      theta_nuis <- matrix(unlist(theta_hat_tilda_b_l),
                           ncol = length(strat_par), byrow = TRUE)
      colnames(theta_nuis) <- name_nuis
      vec_theta_alls <- new_names <- NULL
      kk <- 0
      for (k in seq_len(n_par)) {
        if (k %in% strat_par) {
          which_k <- which(k == strat_par)
          vec_theta_alls <- c(vec_theta_alls, (theta_nuis[, which_k]))
          # fake_name <- name_nuis[which_k]
          pas_nam <- paste0(name_nuis[which_k], rep("_loc", L), seq_len(L))
          new_names <- c(new_names, pas_nam)
        } else {
          kk <- kk + 1
          new_names <- c(new_names, nam_non_nuis[kk])
          vec_theta_alls <- c(vec_theta_alls, as.numeric(theta_hat_a_bfi[kk]))
        }
      }
      names(vec_theta_alls) <- new_names
      theta_hat_bfi <- vec_theta_alls
      A_bfi <- matrix(0, nrow = last_Lam_dim, ncol = last_Lam_dim)
      colnames(A_bfi) <- rownames(A_bfi) <- new_names
      A_bfi[-noncore, -noncore] <- A_a_bfi
      for (j in seq_len(L)) {
        A_bfi[noncore[j, ], noncore[j, ]] <- A_hats_tilda_b_l[[j]]
        A_bfi[noncore[j, ], -noncore] <- A_bfi[-noncore, noncore[j, ]] <-
          A_hats_ab_l[[j]]
      }
      sd_bfi <- sqrt(diag(solve(as.matrix(A_bfi))))
    } else {
      if (length(center_spec) != L) {
        stop("Length of the vector 'center_spec' should be equal to ", sQuote(L))
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
          "The number of levels of 'center_spec' should NOT be equal to the ",
          "number of centers. Otherwise, use 'stratified = TRUE', ",
          "'center_spec = NULL' but set 'strat_par' to 1, 2 or 1:2."
        )
      }
      strat_par <- c(1)
      # For now, we assume 'strat_par' is 'intercept'.
      # However, the input 'strat_par' is NULL.
      if (length(strat_par) == 1) {
        if (1 %in% strat_par) { # For now, we assume 'strat_par' is 'intercept'.
          strat_par <- c(1)
          noncore <- cbind(seq_len(K))
          new_noncore <- noncore
        } else {
          strat_par <- c(length(theta_hats[[1]]))
          noncore <- cbind((last_Lam_dim - K + 1):last_Lam_dim)
          new_noncore <- as.matrix(lev_cen)
        }
      } else {
        strat_par <- c(1, length(theta_hats[[1]]))
        noncore <- new_noncore <- cbind(lev_cen,
                                        (last_Lam_dim - K + 1):last_Lam_dim)
        new_noncore[, 2] <- noncore[, 1] + K
      }
      if (names(theta_hats[[1]])[1] != "(Intercept)" & (1 %in% strat_par)) {
        stop(
          "It seems 'intercept = FALSE ' for centers, while in the current ",
          "version of the package the center specific covariate is ",
          "only possible for 'intercept'!"
        )
      }
      for (j in seq_len(L)) {
        if (j == 1) {
          A_hats_1a_l <- list(A_hats[[j]][-strat_par, -strat_par, drop = FALSE ])
          Lambda_1a_l <- list(Lambda_all[[j]][-strat_par, -strat_par, drop = FALSE ])
          Lambda_1a_l[L + 1] <-
            list(Lambda_all[[L + 1]][-as.numeric(noncore),
                                     -as.numeric(noncore), drop = FALSE ])
          A_hats_1b_l <- list(A_hats[[j]][strat_par, strat_par, drop = FALSE ])
          Lambda_1b_l <- list(Lambda_all[[j]][strat_par,
                                              strat_par, drop = FALSE ])
          Lambda_1b_l[L + 1] <-
            list(Lambda_all[[L + 1]][as.numeric(noncore),
                                     as.numeric(noncore), drop = FALSE ])
          A_hats_1ab_l <- list(A_hats[[j]][-strat_par,
                                           strat_par, drop = FALSE ])
          theta_hat_1a_l <- list(theta_hats[[j]][-strat_par])
          theta_hat_1b_l <- list(theta_hats[[j]][strat_par])
        } else {
          A_hats_1a_l[j] <- list(A_hats[[j]][-strat_par,
                                             -strat_par, drop = FALSE ])
          Lambda_1a_l[j] <- list(Lambda_all[[j]][-strat_par,
                                                 -strat_par, drop = FALSE ])
          A_hats_1b_l[j] <- list(A_hats[[j]][strat_par,
                                             strat_par, drop = FALSE ])
          Lambda_1b_l[j] <- list(Lambda_all[[j]][strat_par,
                                                 strat_par, drop = FALSE ])
          A_hats_1ab_l[j] <- list(A_hats[[j]][-strat_par,
                                              strat_par, drop = FALSE ])
          theta_hat_1a_l[j] <- list(theta_hats[[j]][-strat_par])
          theta_hat_1b_l[j] <- list(theta_hats[[j]][strat_par])
        }
      }
      if (dim(A_hats_1a_l[[1]])[1] != dim(Lambda_1a_l[[L + 1]])[1]) {
        stop("Check the dimension of the last matrix in Lambda.
             It should be changed.")
      }
      A_1bl_thet1bl <- Map("%*%", A_hats_1b_l, theta_hat_1b_l)
      A_1abl_thet1al <- Map("%*%", lapply(seq_len(L), function(x)
        t(A_hats_1ab_l[[x]])), theta_hat_1a_l)
      for (k in seq_len(K)) {
        z_l_k <- which(center_spec == lev_cen[k])
        if (k == 1) {
          which_lev <- which(lev_cen[k] == lev_cen)
          # should be edited if 'strat_par' is not 'intercept'
          A_hats_1bk <- list(Reduce("+", A_hats_1b_l[z_l_k]) +
                               Lambda_1b_l[[L + 1]][which_lev, which_lev] -
            Reduce("+", Lambda_1b_l[z_l_k])) # !
          A_hats_1bk_inv <- list(solve(as.matrix(A_hats_1bk[[k]])))
          A_hats_1abk <- list(Reduce("+", A_hats_1ab_l[z_l_k]))
          map21 <- list(Reduce("+", A_1bl_thet1bl[z_l_k]))
          map22 <- list(Reduce("+", A_1abl_thet1al[z_l_k]))
        } else {
          which_lev <- which(lev_cen[k] == lev_cen)
          # should be edited if 'strat_par' is not 'intercept'
          A_hats_1bk[k] <- list(Reduce("+", A_hats_1b_l[z_l_k]) +
                                  Lambda_1b_l[[L + 1]][which_lev, which_lev] -
            Reduce("+", Lambda_1b_l[z_l_k])) # !
          A_hats_1bk_inv[k] <- list(solve(as.matrix(A_hats_1bk[[k]])))
          A_hats_1abk[k] <- list(Reduce("+", A_hats_1ab_l[z_l_k]))
          map21[k] <- list(Reduce("+", A_1bl_thet1bl[z_l_k]))
          map22[k] <- list(Reduce("+", A_1abl_thet1al[z_l_k]))
        }
      }
      A_1a_bfi <- Reduce("+", A_hats_1a_l) + Lambda_1a_l[[L + 1]] -
        Reduce("+", Lambda_1a_l[seq_len(L)])
      A_hat_ab2 <- Map("%*%", A_hats_1abk, A_hats_1bk_inv)
      A_hat_ab3 <- Map("%*%", A_hat_ab2,
                       lapply(seq_len(K), function(x) t(A_hats_1abk[[x]])))
      theta_hat_a_1 <- A_1a_bfi - Reduce("+", A_hat_ab3)
      map11 <- Map("%*%", A_hats_1a_l, theta_hat_1a_l)
      map12 <- Map("%*%", A_hats_1ab_l, theta_hat_1b_l)
      map1 <- Reduce("+", map11) + Reduce("+", map12)
      map02 <- Map("+", map21, map22)
      map2 <- Map("%*%", A_hat_ab2, map02)
      theta_hat_a_2 <- map1 - Reduce("+", map2)
      theta_hat_1a_bfi <- solve(as.matrix(theta_hat_a_1)) %*% theta_hat_a_2
      part1 <- Map("%*%", A_hats_1b_l, theta_hat_1b_l)
      for (k in seq_len(K)) {
        if (k == 1) {
          theta_hat_1a_bfi_list <- list(theta_hat_1a_bfi)
        } else {
          theta_hat_1a_bfi_list[k] <- list(theta_hat_1a_bfi)
        }
      }
      part1 <- Map("%*%", lapply(seq_len(K), function(x) t(A_hats_1abk[[x]])),
                   theta_hat_1a_bfi_list)
      part2 <- Map("-", map02, part1)
      theta_hat_1bk_bfi <- Map("%*%", A_hats_1bk_inv, part2)
      name_nuis <- names(theta_hats[[1]])[strat_par]
      nam_non_nuis <- names(theta_hats[[1]])[-strat_par]
      theta_nuis <- matrix(unlist(theta_hat_1bk_bfi), ncol = length(strat_par),
                           byrow = TRUE)
      colnames(theta_nuis) <- name_nuis
      vec_theta_alls <- new_names <- NULL
      kk <- 0
      for (k in seq_len(n_par)) {
        if (k %in% strat_par) {
          which_k <- which(k == strat_par)
          vec_theta_alls <- c(vec_theta_alls, (theta_nuis[, which_k]))
          pas_nam <- paste0(name_nuis[which_k], rep("_", K), lev_cen)
          new_names <- c(new_names, pas_nam)
        } else {
          kk <- kk + 1
          new_names <- c(new_names, nam_non_nuis[kk])
          vec_theta_alls <- c(vec_theta_alls, as.numeric(theta_hat_1a_bfi[kk]))
        }
      }
      names(vec_theta_alls) <- new_names
      theta_hat_bfi <- vec_theta_alls
      if (length(theta_hat_bfi) != last_Lam_dim) {
        stop("Check the dimension of the last matrix in Lambda.
             It should be changed.")
      }
      A_bfi <- matrix(0, nrow = last_Lam_dim, ncol = last_Lam_dim)
      colnames(A_bfi) <- rownames(A_bfi) <- new_names
      A_bfi[-noncore, -noncore] <- A_1a_bfi
      for (k in seq_len(K)) {
        A_bfi[noncore[k, ], noncore[k, ]] <- A_hats_1bk[[k]]
        A_bfi[noncore[k, ], -noncore] <- A_bfi[-noncore, noncore[k, ]] <-
          A_hats_1abk[[k]]
      }
      sd_bfi <- sqrt(diag(solve(as.matrix(A_bfi))))
    }
  }
  output <- list(theta_hat = theta_hat_bfi, A_hat = A_bfi, sd = sd_bfi,
                 family = family, stratified = stratified)
  class(output) <- "bfi"
  return(output)
}
