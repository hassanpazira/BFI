# Regression tests for the corrected Gaussian BFI parameterization.
# Put this file in tests/testthat/ in the MODIFIED BFI package.
# The internal parameter is eta = log(sigma2); the public estimate is sigma2.
# Gaussian precision matrices in these tests are block diagonal.

bfi_gaussian_fixture <- function(n_per_center = 85L, seed = 20261023L) {
  set.seed(seed)
  L <- 3L
  n <- L * n_per_center
  X <- data.frame(x1 = stats::rnorm(n))
  center <- rep(seq_len(L), each = n_per_center)
  y <- 0.5 + 0.8 * X$x1 +
    stats::rnorm(n, sd = sqrt(c(1.2, 1.9, 1.5)[center]))
  parts <- split(seq_len(n), center)
  Lambda <- inv.prior.cov(
    X = X, lambda = c(0.1, 0.4), family = "gaussian"
  )
  fits <- lapply(parts, function(i) {
    MAP.estimation(
      y = y[i], X = X[i, , drop = FALSE],
      family = "gaussian", Lambda = Lambda,
      control = list(maxit = 2000)
    )
  })
  list(X = X, y = y, parts = parts, L = L,
       Lambda = Lambda, fits = fits)
}

bfi_gaussian_manual_objective <- function(z, y, Xd, Lambda, w) {
  p <- ncol(Xd)
  beta <- z[seq_len(p)]
  eta <- z[p + 1L]
  residual <- y - drop(Xd %*% beta)
  0.5 * (
    sum(w * residual^2) * exp(-eta) +
      sum(w) * eta +
      drop(crossprod(z, Lambda %*% z))
  )
}

test_that("weighted Gaussian objective and analytic Hessian agree with independent formulas", {
  set.seed(1701)
  n <- 95L
  Xd <- cbind("(Intercept)" = 1, x1 = stats::rnorm(n))
  y <- 0.4 + 0.7 * Xd[, "x1"] + stats::rnorm(n, sd = 1.3)
  w <- seq(0.7, 2.2, length.out = n)
  nm <- c(colnames(Xd), "sigma2")
  Lambda <- diag(c(0.15, 0.3, 0.6))
  dimnames(Lambda) <- list(nm, nm)
  z <- c(0.15, 0.65, log(1.7))

  Q_package <- function(z) {
    as.numeric(negloglik.theta(
      beta_omega = z, y = y, X = Xd, Lambda = Lambda,
      family = "gaussian", wli = w
    ))
  }
  Q_manual <- function(z) {
    bfi_gaussian_manual_objective(z, y, Xd, Lambda, w)
  }

  expect_gt(abs(sum(w) - n), 1)
  expect_equal(Q_package(z), Q_manual(z), tolerance = 1e-11)

  theta_public <- setNames(c(z[1:2], exp(z[3])), nm)
  A_package <- A.l.maker(
    y = y, X = Xd, Lambda = Lambda,
    family = "gaussian", theta_hat = theta_public, wli = w
  )
  A_numeric <- stats::optimHess(par = z, fn = Q_manual)
  relative_error <- max(abs(A_package - A_numeric)) /
    max(abs(A_package))

  expect_lt(relative_error, 1e-4)
  expect_true(isSymmetric(A_package))
  expect_identical(rownames(A_package), nm)
})

test_that("Gaussian MAP reports raw sigma2, but its curvature uses log(sigma2)", {
  dat <- bfi_gaussian_fixture()
  expect_true(all(vapply(dat$fits, function(f) f$convergence == 0, logical(1))))

  for (l in seq_len(dat$L)) {
    fit <- dat$fits[[l]]
    idx <- dat$parts[[l]]
    Xd <- stats::model.matrix(~ x1, data = dat$X[idx, , drop = FALSE])
    nm <- c(colnames(Xd), "sigma2")
    expect_s3_class(fit, "bfi")
    expect_identical(names(fit$theta_hat), nm)
    expect_gt(unname(fit$theta_hat["sigma2"]), 0)
    z <- unname(fit$theta_hat)
    z[length(z)] <- log(z[length(z)])
    Q <- function(z) bfi_gaussian_manual_objective(
      z, dat$y[idx], Xd, dat$Lambda, rep(1, length(idx))
    )
    expect_equal(as.numeric(fit$value), Q(z), tolerance = 1e-8)
    A_numeric <- stats::optimHess(par = z, fn = Q)
    expect_lt(max(abs(fit$A_hat - A_numeric)) /
                max(abs(fit$A_hat)), 1e-4)
  }
})

test_that("non-stratified Gaussian BFI aggregates on eta scale and uses delta method", {
  dat <- bfi_gaussian_fixture()
  local_theta <- lapply(dat$fits, `[[`, "theta_hat")
  local_A <- lapply(dat$fits, `[[`, "A_hat")
  fit <- bfi(
    theta_hats = local_theta, A_hats = local_A,
    Lambda = dat$Lambda, family = "gaussian"
  )

  A_ref <- Reduce("+", local_A) - (dat$L - 1L) * dat$Lambda
  rhs <- Reduce("+", Map(function(A, th) {
    z <- unname(th)
    z[length(z)] <- log(z[length(z)])
    drop(A %*% z)
  }, local_A, local_theta))

  z_ref <- drop(solve(A_ref, rhs))
  theta_ref <- z_ref
  theta_ref[length(theta_ref)] <- exp(theta_ref[length(theta_ref)])
  nm <- names(local_theta[[1]])
  names(theta_ref) <- nm

  expect_equal(unname(coef(fit)), unname(theta_ref), tolerance = 1e-8)
  expect_equal(unname(as.matrix(fit$A_hat)), unname(A_ref), tolerance = 1e-8)
  expect_equal(unname(as.numeric(fit$sd)),
               unname(sqrt(diag(solve(A_ref)))), tolerance = 1e-8)

  d <- rep(1, length(theta_ref))
  d[length(d)] <- theta_ref["sigma2"]
  V_ref <- solve(A_ref) * tcrossprod(d)
  expect_equal(unname(vcov(fit)), unname(V_ref), tolerance = 1e-9)
  expect_identical(names(coef(fit)), nm)
  expect_gt(unname(coef(fit)["sigma2"]), 0)

  s <- summary(fit)
  expect_s3_class(s, "summary.bfi")
  reported <- names(s$estimate)
  expect_equal(unname(as.numeric(s$sd)),
               unname(sqrt(diag(V_ref))[match(reported, nm)]),
               tolerance = 1e-8)
  expect_true(length(capture.output(print(s))) > 0L)
})

test_that("all Gaussian stratified cases match an independently embedded quadratic", {
  dat <- bfi_gaussian_fixture()
  theta_list <- lapply(dat$fits, `[[`, "theta_hat")
  A_list <- lapply(dat$fits, `[[`, "A_hat")

  for (strat_par in list(1L, 2L, c(1L, 2L))) {
    Lambda_global <- inv.prior.cov(
      X = dat$X, lambda = c(0.1, 0.4), family = "gaussian",
      stratified = TRUE, strat_par = strat_par, L = dat$L
    )
    fit <- bfi(
      theta_hats = theta_list, A_hats = A_list,
      Lambda = list(dat$Lambda, Lambda_global),
      family = "gaussian", stratified = TRUE, strat_par = strat_par
    )

    nm <- names(coef(fit))
    expect_setequal(nm, rownames(Lambda_global))
    Lambda_g <- Lambda_global[nm, nm, drop = FALSE]
    A_ref <- Lambda_g
    rhs <- setNames(numeric(length(nm)), nm)

    # S_l maps global parameters to the center-l parameter vector.
    # This formulation is independent of the package's Schur-complement code.
    for (l in seq_len(dat$L)) {
      local_names <- names(theta_list[[l]])
      global_names <- local_names
      if (1L %in% strat_par) {
        global_names[local_names == "(Intercept)"] <-
          paste0("(Intercept)_loc", l)
      }
      if (2L %in% strat_par) {
        global_names[local_names == "sigma2"] <-
          paste0("sigma2_loc", l)
      }

      S <- matrix(0, nrow = length(local_names), ncol = length(nm),
                  dimnames = list(local_names, nm))
      expect_true(all(global_names %in% nm))
      S[cbind(seq_along(local_names), match(global_names, nm))] <- 1
      A_l <- as.matrix(A_list[[l]])
      z_l <- unname(theta_list[[l]])
      z_l[local_names == "sigma2"] <- log(z_l[local_names == "sigma2"])

      A_ref <- A_ref + crossprod(S, (A_l - dat$Lambda) %*% S)
      rhs <- rhs + drop(crossprod(S, A_l %*% z_l))
    }

    z_ref <- drop(solve(A_ref, rhs))
    names(z_ref) <- nm
    sigma_idx <- grepl("^sigma2($|_loc[0-9]+$)", nm)
    theta_ref <- z_ref
    theta_ref[sigma_idx] <- exp(theta_ref[sigma_idx])

    expect_equal(unname(coef(fit)), unname(theta_ref), tolerance = 1e-7)
    expect_equal(unname(as.matrix(fit$A_hat)), unname(A_ref), tolerance = 1e-8)
    expect_true(all(coef(fit)[sigma_idx] > 0))

    d <- rep(1, length(nm))
    d[sigma_idx] <- theta_ref[sigma_idx]
    V_ref <- solve(A_ref) * tcrossprod(d)
    expect_equal(unname(vcov(fit)), unname(V_ref), tolerance = 1e-8)

    s <- summary(fit)
    expect_s3_class(s, "summary.bfi")
    reported <- names(s$estimate)
    expect_equal(unname(as.numeric(s$sd)),
                 unname(sqrt(diag(V_ref))[match(reported, nm)]),
                 tolerance = 1e-8)
    if (any(grepl("^sigma2", reported))) {
      sigma_names <- reported[grepl("^sigma2", reported)]
      se_eta <- sqrt(diag(solve(A_ref)))
      zcrit <- stats::qnorm(0.975)
      ci_manual <- cbind(
        theta_ref[sigma_names] * exp(-zcrit * se_eta[sigma_names]),
        theta_ref[sigma_names] * exp(+zcrit * se_eta[sigma_names])
      )
      expect_equal(unname(as.matrix(s$CI)[sigma_names, , drop = FALSE]),
                   unname(ci_manual), tolerance = 1e-8)
    }
  }
})

test_that("two-round Gaussian IPTW agrees with independently weighted objective and ATE", {
  set.seed(20261024)
  L <- 3L
  n_per <- 85L
  n <- L * n_per
  X <- data.frame(x1 = stats::rnorm(n))
  ps <- stats::plogis(-0.25 + 0.7 * X$x1)
  X$trt <- stats::rbinom(n, 1, ps)
  y <- 0.4 + 1.2 * X$trt + 0.6 * X$x1 +
    stats::rnorm(n, sd = sqrt(1.4))
  parts <- split(seq_len(n), rep(seq_len(L), each = n_per))

  Lambda_ps <- inv.prior.cov(
    X = X, lambda = 0.1, family = "binomial",
    treatment = "trt", treat_round = "first"
  )
  Lambda_out <- inv.prior.cov(
    X = X, lambda = c(0.1, 0.4), family = "gaussian",
    treatment = "trt", treat_round = "second"
  )

  local_ps <- lapply(parts, function(i) {
    MAP.estimation(
      y = y[i], X = X[i, , drop = FALSE], family = "binomial",
      Lambda = Lambda_ps, treatment = "trt", treat_round = "first",
      control = list(maxit = 2000)
    )
  })
  expect_true(all(vapply(local_ps, function(f) f$convergence == 0, logical(1))))
  bfi_ps <- bfi(
    theta_hats = lapply(local_ps, `[[`, "theta_hat"),
    A_hats = lapply(local_ps, `[[`, "A_hat"),
    Lambda = Lambda_ps, family = "binomial", treat_round = "first"
  )
  gamma <- coef(bfi_ps)
  expect_identical(names(gamma), c("(Intercept)", "x1"))

  local_out <- lapply(parts, function(i) {
    MAP.estimation(
      y = y[i], X = X[i, , drop = FALSE], family = "gaussian",
      Lambda = Lambda_out, treatment = "trt", treat_round = "second",
      gamma_bfi = gamma, control = list(maxit = 2000)
    )
  })
  expect_true(all(vapply(local_out, function(f) f$convergence == 0, logical(1))))

  for (l in seq_len(L)) {
    i <- parts[[l]]
    fit <- local_out[[l]]
    e <- as.numeric(fit$propensity)
    Z <- X$trt[i]
    w <- Z / e + (1 - Z) / (1 - e)
    Xd <- stats::model.matrix(~ trt, data = X[i, , drop = FALSE])
    z <- unname(fit$theta_hat)
    z[length(z)] <- log(z[length(z)])
    Q <- function(z) bfi_gaussian_manual_objective(
      z, y[i], Xd, Lambda_out, w
    )
    Q_package <- as.numeric(negloglik.theta(
      beta_omega = z, y = y[i], X = Xd, Lambda = Lambda_out,
      family = "gaussian", wli = w
    ))
    expect_equal(Q_package, Q(z), tolerance = 1e-9)
    expect_equal(as.numeric(fit$value), Q(z), tolerance = 1e-8)
    expect_lt(max(abs(stats::optimHess(z, Q) - fit$A_hat)) /
                max(abs(fit$A_hat)), 1e-4)
  }

  fit_bfi <- bfi(
    theta_hats = lapply(local_out, `[[`, "theta_hat"),
    A_hats = lapply(local_out, `[[`, "A_hat"),
    Lambda = Lambda_out, family = "gaussian",
    treat_round = "second",
    for_ATE = lapply(local_out, `[[`, "for_ATE")
  )
  expect_gt(unname(coef(fit_bfi)["sigma2"]), 0)

  X_ps <- stats::model.matrix(~ x1, data = X)
  e <- drop(stats::plogis(X_ps %*% gamma))
  Z <- X$trt
  ATE_manual <- c(
    IPTW = mean(Z * y / e - (1 - Z) * y / (1 - e)),
    wIPTW = sum(Z * y / e) / sum(Z / e) -
      sum((1 - Z) * y / (1 - e)) / sum((1 - Z) / (1 - e))
  )
  expect_equal(unname(unlist(fit_bfi$Ave_Treat)),
               unname(ATE_manual), tolerance = 1e-9)
})
