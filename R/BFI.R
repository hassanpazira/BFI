## This file created by Hassan Pazira

#' @export

bfi <- function(theta_hats = NULL, A_hats, Lambda,
                family = c("gaussian", "binomial", "survival"),
                basehaz = c("weibul","exp","gomp","poly","pwexp","unspecified"),
                stratified = FALSE, strat_par = NULL, center_spec = NULL,
                theta_A_polys = NULL, treat_round = NULL, for_ATE = NULL, p, q_ls,
                center_zero_sample = FALSE, which_cent_zeros, zero_sample_covs,
                refer_cats, zero_cats, lev_no_ref_zeros) {
  family <- match.arg(family)
  basehaz <- match.arg(basehaz)
  if (!is.null(theta_hats) & !is.list(theta_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.null(treat_round)) {
    if (treat_round == "first" & !is.null(for_ATE)) stop("In the first round, 'for_ATE' must be 'NULL'.")
    if (treat_round == "first" & is.null(for_ATE)) family <- c("binomial")
    if (treat_round == "second") {
      if (family == "survival") {
        if (!is.null(for_ATE))
          stop("In the second round, 'for_ATE' must be NULL for survival.")
      } else {
        if (is.null(for_ATE)) stop("In the second round, 'for_ATE' must be a list.")
      }
    }
  } else {
    if (!is.null(for_ATE)) stop("If 'treat_round' is NULL, 'for_ATE' must be NULL.")
  }
  if (family == "survival" & basehaz == "poly" & is.null(theta_A_polys)) {
    stop("When family='survival' and basehaz='poly', 'theta_A_polys' cannot be NULL.")
  }
  if (missing(A_hats) & family != "survival") {
    stop("Argument 'A_hats' cannot be missing.")
  }
  if (family == "survival" & basehaz != "poly" & missing(A_hats) &
      !is.null(theta_A_polys)) {
    stop("'basehaz' should be 'poly'.")
  }
  if (family == "survival" & basehaz != "poly" & missing(A_hats)) {
    stop("Argument 'A_hats' cannot be missing.")
  }
  if (!is.null(theta_A_polys) & !is.list(theta_A_polys)) {
    stop("The input for 'theta_A_polys' should be a list with L elements.")
  }
  if (family == "survival" & basehaz == "poly" & missing(q_ls)) {
    stop("When family='survival' and basehaz='poly', 'q_ls' cannot be missing.")
  }
  if (!is.null(theta_A_polys) & length(theta_A_polys)<2) {
    stop("Length of 'theta_A_polys' should be equal to L (which is > 1).")
  }
  if (is.null(theta_A_polys)) {
    if (!is.list(A_hats)) {
      stop("The input for 'A_hats' should be a list.")
    }
  }
  if (!is.null(theta_hats)) {
    if ((length(theta_hats) != length(A_hats)))
      stop("Length of inputs are not equal.")
  }
  if (is.null(theta_hats) & is.null(theta_A_polys)) {
    stop("Argument 'theta_hats' is missing, with no default.")
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
  if (stratified == TRUE & !is.null(strat_par) & family == "survival" & basehaz == "unspecified") {
    stop("For the 'unspecified' baseline hazard, 'strat_par' can not be used.")
  }
  old_strat_par <- strat_par
  strat_par <- sort(strat_par)
  if (is.data.frame(Lambda)) {
    Lambda <- as.matrix(Lambda)
  }
  if (!(family == "survival" & basehaz == "poly")) {
    L <- length(A_hats)
    for (l in seq_len(L)) {
      if (l == 1) equal_names_theta <- equal_names_A_hat <- list()
      if (is.null(names(theta_hats[[l]]))) stop("Names of 'theta_hat's cannot be NULL.")
      equal_names_theta[[l]] <- names(theta_hats[[l]])
      if (is.null(colnames(A_hats[[l]]))) stop("Colnames of 'A_hat's cannot be NULL.")
      equal_names_A_hat[[l]] <- colnames(A_hats[[l]])
    }
    if (!all(equal_names_theta[[1]] == unlist(equal_names_theta))) {
      if (all(sapply(equal_names_theta, function(x) all(sort(x) == sort(equal_names_theta[[1]]))))) {
        theta_hats <- lapply(theta_hats, function(x) x[order(names(x))])
      } else {
        stop("Names of 'theta_hat's are not the same across centers")
      }
    }
    if (!all(equal_names_A_hat[[1]] == unlist(equal_names_A_hat))) {
      if (all(sapply(equal_names_A_hat, function(x) all(sort(x) == sort(equal_names_A_hat[[1]]))))) {
        A_hats <- lapply(A_hats, function(x) x[order(rownames(x)), order(colnames(x))])
      } else {
        stop("Colnames/Rownames of 'A_hat's are not the same across centers")
      }
    }
  } else {
    L <- length(theta_A_polys)
    for (l in seq_len(L)) {
      if (l == 1) equal_names_theta <- equal_names_A_hat <- list()
      if (is.null(rownames(theta_A_polys[[l]][,,1])))
        stop("Names of 'theta_hat's in 'theta_A_polys' cannot be NULL.")
      equal_names_theta[[l]] <- rownames(theta_A_polys[[l]][,,1])
      if (is.null(rownames(theta_A_polys[[l]][,,2])))
        stop("Colnames of 'A_hat's in 'theta_A_polys' cannot be NULL.")
      equal_names_A_hat[[l]] <- rownames(theta_A_polys[[l]][,,2])
    }
    if (!all(equal_names_theta[[1]] == unlist(equal_names_theta))) {
      if (all(sapply(equal_names_theta, function(x) all(sort(x) == sort(equal_names_theta[[1]]))))) {
        theta_A_polys_changes_theta <- lapply(seq_along(theta_A_polys), function(l) {
          x <- theta_A_polys[[l]]
          nr <- nrow(x[,,1])
          q_lsr <- q_ls[l]
          if ((q_lsr+1) >= nr) stop("q_ls[l] is too large; cannot exclude all rows from sorting.")
          sorted_indices <- order(rownames(x[1:(nr - q_lsr - 1), , 1]))
          x[c(sorted_indices, (nr - q_lsr):nr), , 1]
        })
        for (ll in seq_along(theta_A_polys)) {
          theta_A_polys[[ll]][,,1] <- theta_A_polys_changes_theta[[ll]]
        }
      } else {
        stop("Names of 'theta_hat's are not the same across centers")
      }
    }
    if (!all(equal_names_A_hat[[1]] == unlist(equal_names_A_hat))) {
      if (all(sapply(equal_names_A_hat, function(x) all(sort(x) == sort(equal_names_A_hat[[1]]))))) {
        theta_A_polys_changes_A <- lapply(seq_along(theta_A_polys), function(l) {
          x <- theta_A_polys[[l]]
          nr <- nrow(x[,,1])
          q_lsr <- q_ls[l]
          if ((q_lsr+1) >= nr) stop("q_ls[l] is too large; cannot exclude all rows from sorting.")
          sorted_indices <- order(rownames(x[1:(nr - q_lsr - 1), , 1]))
          x[c(sorted_indices, (nr - q_lsr):nr), c(sorted_indices, (nr - q_lsr):nr), 2:(q_lsr+2)]
        })
        for (ll in seq_along(theta_A_polys)) {
          theta_A_polys[[ll]][,,-1] <- theta_A_polys_changes_A[[ll]]
          # Directly modifying rownames(theta_A_polys[[ll]][,,1]) does NOT work in a 3D array!
          dimnames(theta_A_polys[[ll]])[[1]] <- rownames(theta_A_polys_changes_theta[[ll]])
        }
      } else {
        stop("Colnames/Rownames of 'A_hat's are not the same across centers")
      }
    }
  }
  if (L == 1) stop("The number of locations should be > 1.")
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
          stop("Lambda should contain 'L+1' lists; 1:L for local datastes and last one for combined.")
        }
        Lambda_all <- Lambda
      }
    } else {
      stop("Lambda should be a list.")
    }
  }
  if (family == "survival" & basehaz == "poly") {
    if (length(q_ls) != 1) {
      if (length(q_ls) != L) stop("Length of 'q_ls' should be ", sQuote(L),".")
    }
    q_max <- max(q_ls)
    omega_len <- q_max + 1
    p <- dim(theta_A_polys[[1]])[[1]] - length(grep("omega",rownames(theta_A_polys[[1]])))
    theta_len <- p + omega_len  # for the zero-patients case, 'p' can differ!
    colnames_X <- rownames(theta_A_polys[[1]])[1:p]
    for (i in seq_len(L)) {
      if (ncol(Lambda_all[[i]]) < theta_len)
        stop("Dimension of Lambda should be at least ", sQuote(theta_len),".")
    }
    theta_hats <- A_hats <- list()
    for (l in seq_len(L)) {
      theta_hats[[l]] <- theta_A_polys[[l]][1:theta_len, omega_len, 1]
      A_hats[[l]] <- theta_A_polys[[l]][1:theta_len, 1:theta_len, omega_len+1]
      colnames(A_hats[[l]]) <- rownames(A_hats[[l]])
      if (any(diag(solve(as.matrix(A_hats[[l]]))) < 0)) {
        stop("In the center,", sQuote(l),", some diagonal elements of the inverse",
             "curvature matrix are negative!")
      }
      # Zero-patient category
      if (center_zero_sample == TRUE) {
        if (missing(which_cent_zeros))
          stop("If there is a categorical covariate with zero sample in at least",
               "one of the centers, 'which_cent_zeros' should not be missing.")
        if (missing(zero_sample_covs))
          stop("If there is a categorical covariate with zero sample in only one of",
               "the categories, 'zero_sample_covs' should not be missing.")
        if (!is.numeric(which_cent_zeros)) stop("'which_cent_zeros' should be numeric.")
        if (!is.list(lev_no_ref_zeros)) stop("'lev_no_ref_zeros' should be a list.")
        if (length(which_cent_zeros) != length(zero_sample_covs) |
            length(which_cent_zeros) != length(refer_cats) |
            length(which_cent_zeros) != length(zero_cats) |
            length(which_cent_zeros) != length(lev_no_ref_zeros) )
          stop("Length of 'which_cent_zeros', 'zero_sample_covs', 'refer_cats',",
               "'zero_cats' and 'lev_no_ref_zeros' should be equal.")
        if (missing(refer_cats)) stop("The reference category is missed.")
        if (length(zero_sample_covs) != length(refer_cats))
          stop("'zero_sample_covs' and 'refer_cats' must have the same length.")
        if (!is.character(zero_cats)) stop("'zero_cats' should be character.")
        if (!is.character(refer_cats)) stop("'refer_cats' should be character.")
        if (!is.character(zero_sample_covs)) stop("'zero_sample_covs' should be character.")
        if (missing(zero_cats)) stop("The category with zero sample is missed.")
        if (any(zero_cats == refer_cats))
          stop("The category with no patient cannot be used as the reference.")
        if (max(which_cent_zeros) > L | min(which_cent_zeros) < 1)
          stop("'which_cent_zeros' should be from 1 to L.")
        for (wc in 1:length(which_cent_zeros)) {
          wich_cent <- which_cent_zeros[wc]
          # For each center: only one covariate and one 'zero_cats'!
          lev_zero_cov <- sort(as.character(c(lev_no_ref_zeros[[wc]][-1], zero_cats[wc])))
          names_after_dammy <- paste(zero_sample_covs[wc],lev_zero_cov, sep="")
          which_cat_zero <- which(lev_zero_cov %in% zero_cats[wc])
          names_cat_zero <- names_after_dammy[which_cat_zero]
          if (length(names_cat_zero) != length(zero_cats[wc]))
            stop("length(names_cat_zero) != length(zero_cats[wc])")
          # Positions to add the new elements
          for (wi in seq_len(length(zero_cats[wc]))) {
            # This 'for' loop will be updated for more than one 'zero_cats'!
            if (wi == 1) which_element <- NULL
            if (which_cat_zero == 1) {
              which_element[wi] <- which(colnames_X == names_after_dammy[2])
            } else {
              if (length(lev_zero_cov) > 2) {
                if (length(names_after_dammy) == which_cat_zero) {
                  which_element[wi] <- which(colnames_X ==
                                               names_after_dammy[which_cat_zero-1]) - 1
                } else {
                  which_element[wi] <- which(colnames_X ==
                                               names_after_dammy[which_cat_zero + 1])
                }
              } else {
                which_element[wi] <- which(colnames_X == names_after_dammy[1]) + 1
              }
            }
          }
          names_initial_vec_beta <- colnames_X
          # Add the new elements to the vector
          for (i in seq_len(length(which_element))) { # or seq_along(names_cat_zero)
            names_initial_vec_beta <- append(names_initial_vec_beta, names_cat_zero[i],
                                             after = which_element[i] - 1)
          }
          initial_vec_beta <- rep(0, (p + length(names_cat_zero)))
          names(initial_vec_beta) <- names_initial_vec_beta
          initial_vec_beta[-which_element] <- theta_hats[[wich_cent]][1:p]
          initial_vec_beta <- c(initial_vec_beta,
                                theta_hats[[wich_cent]][(p+1):length(theta_hats[[wich_cent]])])
          theta_hats[[wich_cent]] <- initial_vec_beta

          # A_hats[[wich_cent]]
          # create a new 'Gamma_dot' for all parameters!
          Gamma <- Lambda_all[[wich_cent]][c(1:p),c(1:p)]
          Gamma_new <- Gamma
          if (ncol(Gamma_new) < which_element[1] &
              abs(ncol(Gamma_new) - which_element[1])==1) {
            for (ii in seq_len(length(which_element))) {
              Gamma_new <- rbind(Gamma_new,
                                 c(Gamma_new[ncol(Gamma_new),1],
                                   Gamma_new[ncol(Gamma_new),-ncol(Gamma_new)]))
            }
            # Add the columns at the specified position to A_hats[[wich_cent]]
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
                                   Gamma_new[which_element[ii]+1,
                                             (which_element[ii]+2):ncol(Gamma_new)]),
                                 Gamma_new[which_element[ii]:nrow(Gamma_new), ])
            }
            # Add the columns at the specified position to A_hats[[wich_cent]]
            for (jj in seq_len(length(which_element))) {
              Gamma_new <- cbind(Gamma_new[, 1:(which_element[jj] - 1)],
                                 c(Gamma[,which_element[jj]],
                                   Gamma[ncol(Gamma),which_element[jj]]),
                                 Gamma_new[, which_element[jj]:ncol(Gamma_new)])
            }
          }
          Gamma_dot_new <- b.diag(Gamma_new, Lambda_all[[wich_cent]][-c(1:p),-c(1:p)])

          # Add the rows at the specified position to A_hats[[wich_cent]]
          new_M_with_rows <- A_hats[[wich_cent]]
          for (ii in seq_len(length(which_element))) {
            new_M_with_rows <- rbind(new_M_with_rows[1:(which_element[ii] - 1), ],
                                     Gamma_dot_new[which_element[ii],-which_element[ii]],
                                     new_M_with_rows[which_element[ii]:nrow(new_M_with_rows), ])
          }
          # Add the columns at the specified position to A_hats[[wich_cent]]
          new_M_with_columns <- new_M_with_rows
          for (jj in seq_len(length(which_element))) {
            new_M_with_columns <- cbind(new_M_with_columns[, 1:(which_element[jj] - 1)],
                                        Gamma_dot_new[,which_element[jj]],
                                        new_M_with_columns[, which_element[jj]:ncol(new_M_with_columns)])
            colnames(new_M_with_columns)[which_element[jj]] <- names_cat_zero[jj]
          }
          A_hats[[wich_cent]] <- new_M_with_columns
          Lambda_all[[wich_cent]] <- Gamma_dot_new
        }
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
    stop("All matrices in 'A_hats' should have the same columns.
         Number of parameters in locations must be equal.")
  }
  if (stratified == TRUE & !is.null(strat_par)) {
    if (is.null(theta_hats)) {
      stop("In stratified analysis, 'theta_hats' should not be 'NULL'.")
    }
    if (any(duplicated(strat_par))) stop("There shouldn't be any duplicates in 'strat_par'.")
    if (family == "gaussian") {
      if (!all(strat_par %in% c(1, 2))) {
        stop("'strat_par' should contain '1' and/or '2'.")
      }
      if (!is.numeric(strat_par)) {
        stop("'strat_par' should be one of the integers: 1 or 2, or a vector of both.")
      }
      if (length(strat_par) > 2) {
        stop("For this family ('gaussian'), the number of elements of 'strat_par'
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
        stop("'strat_par' should be an integer.")
      }
      if (length(strat_par) > 1) {
        stop("For this family ('binomial'), the number of elements of 'strat_par' is one, i.e., 'intercept'.")
      }
      if (strat_par == 2) {
        stop("'strat_par' should only be the integer 1")
      }
    }
    if (family == "survival") {
      if (!is.numeric(strat_par)) {
        stop("'strat_par' should contain integers.")
      }
      if (length(strat_par) > (n_par - p)) {
        stop("Maximum length of 'strat_par' can be the number of parameters of the baseline hazard, ",
             (n_par - p))
      }
      if (!all(strat_par %in% seq_len(n_par - p))) {
        stop("All elements of 'strat_par' should be chosen from the integers: ", seq_len(n_par - p))
      }
    }
  }
  last_Lam_dim <- dim(Lambda_all[[L + 1]])[1]
  if (stratified == FALSE) {
    if (last_Lam_dim != dim(Lambda_all[[1]])[1]) {
      stop("The last matrix in 'Lambda' should have equal dim. as the other
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
    if (!is.null(for_ATE)) { # only for "binomial" and 'gaussian'!
      if (!is.list(for_ATE)) {
        stop("The input for 'for_ATE' should be a list.")
      }
      sum_over_L <- Reduce("+", for_ATE)
      m1 <- sum_over_L[1] # treatment
      m2 <- sum_over_L[2] # control
      N <- m1 + m2
      ATE1 <- (sum_over_L[6] - sum_over_L[8]) / N
      ATE2 <- (sum_over_L[6]/sum_over_L[5]) - (sum_over_L[8]/sum_over_L[7])
      ATE <- list(IPTW=ATE1, wIPTW=ATE2)
      # sample_var <- list(
      #   treatment = (sum_over_L[4]/(m1-1)) - (m1/(m1-1))*(sum_over_L[3]/m1)^2,
      #   control = (sum_over_L[4]/(m2-1)) - (m2/(m2-1))*(sum_over_L[3]/m2)^2
      # )
    } else {
      ATE <- NULL
      # sample_var <- NULL
    }
  } else {
    if (last_Lam_dim == dim(Lambda_all[[1]])[1]) {
      if (all(colnames(Lambda_all[[1]]) == colnames(Lambda_all[[L + 1]]))) {
        stop("The last matrix in 'Lambda' should not have the same dim. as the other local matrices.")
      } else {
        if (length(strat_par) == 1 & family != "survival") {
          if (1 %in% strat_par) {
            if (colnames(Lambda_all[[1]])[ncol(Lambda_all[[1]])] !=
                colnames(Lambda_all[[L + 1]])[ncol(Lambda_all[[L + 1]])]) {
              stop("Check the elements of 'Lambda'. Is there an 'Sigma2' in the model?")
            }
          } else {
            if (colnames(Lambda_all[[1]])[1] != colnames(Lambda_all[[L + 1]])[1]) {
              stop("Check the elements of 'Lambda'. Is there an 'intercept' in the model?")
            }
          }
        }
      }
    }
    if (is.null(center_spec)) {
      if (family == "survival") {
        for (i in 1:length(strat_par)) {
          if (i == 1) {
            noncore <- cbind(p + (strat_par[i] - 1) * L + (1:L))
            new_noncore <- cbind(seq_len(L))
          } else {
            noncore <- cbind(noncore, p + (strat_par[i] - 1) * L + (1:L))
            new_noncore <- cbind(new_noncore, new_noncore[, i-1] + L)
          }
        }
        strat_par <- p + strat_par
      } else {
        if (length(strat_par) == 1) {
          if (1 %in% strat_par) {
            strat_par <- c(1)
            noncore <- cbind(seq_len(L))
            new_noncore <- noncore
          } else {
            strat_par <- c(length(theta_hats[[1]]))
            noncore <- cbind((last_Lam_dim - L + 1):last_Lam_dim)
            new_noncore <- as.matrix(seq_len(L)) # == noncore - (L+1)
          }
        } else {
          strat_par <- c(1, length(theta_hats[[1]]))
          noncore <- new_noncore <- cbind(seq_len(L),
                                          (last_Lam_dim - L + 1):last_Lam_dim)
          new_noncore[, 2] <- noncore[, 1] + L
        }
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
                                     Lambda_b_l[[L + 1]][new_noncore[j, ], new_noncore[j, ], drop = FALSE] -
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
                                        Lambda_b_l[[L + 1]][new_noncore[j, ], new_noncore[j, ], drop = FALSE] -
                                        Lambda_b_l[[j]])
          A_hats_tilda_b_l_inv[j] <- list(solve(as.matrix(A_hats_tilda_b_l[[j]])))
          A_hats_ab_l[j] <- list(A_hats[[j]][-strat_par, strat_par, drop = FALSE ])
          theta_hat_a_l[j] <- list(theta_hats[[j]][-strat_par])
          theta_hat_b_l[j] <- list(theta_hats[[j]][strat_par])
        }
      }
      A_a_bfi <- Reduce("+", A_hats_a_l) + Lambda_a_l[[L + 1]] - Reduce("+", Lambda_a_l[seq_len(L)])
      A_hat_ab2 <- Map("%*%", A_hats_ab_l, A_hats_tilda_b_l_inv)
      A_hat_ab3 <- Map("%*%", A_hat_ab2, lapply(seq_len(L), function(x) t(A_hats_ab_l[[x]])))
      theta_hat_a_1 <- A_a_bfi - Reduce("+", A_hat_ab3)
      A_hat_a_l_1 <- Map("-", A_hats_a_l, A_hat_ab3)
      I_d2 <- lapply(seq_len(L), function(x) diag(length(strat_par)))
      A_hat_b2 <- Map("%*%", A_hats_tilda_b_l_inv, A_hats_b_l)
      A_hat_a_l_2 <- Map("-", I_d2, A_hat_b2)
      map1 <- Map("%*%", A_hat_a_l_1, theta_hat_a_l)
      map02 <- Map("%*%", A_hats_ab_l, A_hat_a_l_2)
      map2 <- Map("%*%", map02, theta_hat_b_l)
      theta_hat_a_2 <- Map("+", map1, map2)
      theta_hat_a_bfi <- solve(as.matrix(theta_hat_a_1)) %*% Reduce("+", theta_hat_a_2)
      part1 <- Map("%*%", A_hats_b_l, theta_hat_b_l)
      theta_hat_a_bfi_list <- list(theta_hat_a_bfi)
      for (j in 2:L) theta_hat_a_bfi_list[j] <- list(theta_hat_a_bfi)
      theta_hat_a_l_a_bfi <- Map("-", theta_hat_a_l, theta_hat_a_bfi_list)
      part2 <- Map("%*%", lapply(seq_len(L), function(x)
        t(A_hats_ab_l[[x]])), theta_hat_a_l_a_bfi)
      part12 <- Map("+", part1, part2)
      theta_hat_tilda_b_l <- Map("%*%", A_hats_tilda_b_l_inv, part12)
      name_nuis <- names(theta_hats[[1]])[strat_par]
      nam_non_nuis <- names(theta_hats[[1]])[-strat_par]
      theta_nuis <- matrix(unlist(theta_hat_tilda_b_l), ncol = length(strat_par), byrow = TRUE)
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
        A_bfi[noncore[j, ], -noncore] <- A_bfi[-noncore, noncore[j, ]] <- A_hats_ab_l[[j]]
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
    ATE <- NULL
    # sample_var <- NULL
  }
  if (family != "survival") basehaz <- NULL
  output <- list(theta_hat = theta_hat_bfi, A_hat = A_bfi, sd = sd_bfi, family = family,
                 basehaz = basehaz, stratified = stratified, strat_par = old_strat_par,
                 Ave_Treat = ATE) #S_var = sample_var
  class(output) <- "bfi"
  return(output)
}
