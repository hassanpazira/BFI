
test_that("survival formulas reproduce legacy local and federated fits", {

  skip_if_not_installed("survival")

  # Reproducible dataset with time, event and covariates.
  lung <- survival::lung

  d <- data.frame(
    time = lung$time / 365.25,
    status = as.integer(lung$status == 2L),
    x = as.numeric(scale(lung$age)),
    z = as.integer(lung$sex == 2L)
  )

  d <- d[stats::complete.cases(d), , drop = FALSE]

  set.seed(3301)
  d <- d[sample.int(nrow(d)), , drop = FALSE]

  centers <- split(
    d, rep(seq_len(2L), length.out = nrow(d))
  )

  # Aggregate local outputs using the existing BFI API.
  aggregate_survival <- function(items, basehaz,
                                 round = NULL) {

    args <- list(
      Lambda = items[[1L]]$Lambda,
      family = "survival",
      basehaz = basehaz
    )

    if (basehaz == "poly") {
      args$theta_A_polys <- lapply(
        items, function(a) a$fit$theta_A_poly
      )
      args$q_ls <- vapply(
        items, function(a) a$fit$q_l,
        numeric(1)
      )
    } else {
      args$theta_hats <- lapply(
        items, function(a) a$fit$theta_hat
      )
      args$A_hats <- lapply(
        items, function(a) a$fit$A_hat
      )
    }

    if (!is.null(round)) {
      args$treat_round <- round
    }

    do.call(bfi, args)
  }

  compare_results <- function(old, new, fields) {
    for (field in fields) {
      expect_equal(
        old[[field]], new[[field]],
        tolerance = 1e-8,
        info = paste("Mismatch in", field)
      )
    }
  }

  # -------------------------------------------------------
  # Part 1: All six baseline hazards, two centers.
  # -------------------------------------------------------

  f <- survival::Surv(time, status) ~ x + z

  for (haz in c(
    "exp", "weibul", "gomp",
    "poly", "pwexp", "unspecified"
  )) {

    extra_prior <- list()
    extra_fit <- list()

    if (haz == "poly") {
      extra_prior$max_order <- 2L
      extra_fit$max_order <- 2L
    }

    if (haz == "pwexp") {
      extra_prior$n_intervals <- 3L
      extra_fit$n_intervals <- 3L
      extra_fit$min_max_times <- min(
        vapply(
          centers,
          function(a) max(a$time),
          numeric(1)
        )
      )
    }

    old <- new <- vector("list", length(centers))

    for (i in seq_along(centers)) {

      dat <- centers[[i]]

      X <- dat[, c("x", "z"), drop = FALSE]
      Y <- dat[, c("time", "status"), drop = FALSE]

      base_prior <- list(
        family = "survival",
        lambda = c(0.1, 1),
        basehaz = haz
      )

      L_old <- do.call(
        inv.prior.cov,
        c(list(X = X), base_prior, extra_prior)
      )

      L_new <- do.call(
        inv.prior.cov,
        c(
          list(formula = f, data = dat),
          base_prior, extra_prior
        )
      )

      expect_equal(
        L_old, L_new,
        tolerance = 1e-8,
        info = paste(haz, "prior", i)
      )

      fit_options <- c(
        list(family = "survival", basehaz = haz),
        extra_fit
      )

      fit_old <- do.call(
        MAP.estimation,
        c(
          list(y = Y, X = X, Lambda = L_old),
          fit_options
        )
      )

      fit_new <- do.call(
        MAP.estimation,
        c(
          list(
            formula = f, data = dat,
            Lambda = L_new
          ),
          fit_options
        )
      )

      compare_results(
        fit_old, fit_new,
        c("theta_hat", "A_hat", "sd",
          "convergence", "value")
      )

      expect_identical(class(fit_old), class(fit_new))

      if (haz == "poly") {
        compare_results(
          fit_old, fit_new,
          c("q_l", "theta_A_poly")
        )
      }

      old[[i]] <- list(Lambda = L_old, fit = fit_old)
      new[[i]] <- list(Lambda = L_new, fit = fit_new)
    }

    central_old <- aggregate_survival(old, haz)
    central_new <- aggregate_survival(new, haz)

    compare_results(
      central_old, central_new,
      c("theta_hat", "A_hat")
    )

    expect_identical(
      class(central_old), class(central_new)
    )
  }

  # -------------------------------------------------------
  # Part 2: Survival treatment, two rounds and RCT path.
  # -------------------------------------------------------

  treatment_fit <- function(dat, api, round,
                            gamma = NULL, rct = FALSE) {

    # For testing only: binary treatment indicator.
    X <- data.frame(
      x = dat$x,
      treatment = dat$treatment
    )

    Y <- dat[, c("time", "status"), drop = FALSE]

    formula <- survival::Surv(time, status) ~
      x + treatment - 1

    haz <- if (round == "first") {
      "weibul"
    } else {
      "unspecified"
    }

    prior_args <- list(
      family = "survival",
      basehaz = haz,
      lambda = 0.01,
      treatment = "treatment",
      treat_round = round
    )

    fit_args <- list(
      family = "survival",
      basehaz = haz,
      treatment = "treatment",
      treat_round = round
    )

    if (api == "formula") {
      prior_args$formula <- formula
      prior_args$data <- dat
      fit_args$formula <- formula
      fit_args$data <- dat
    } else {
      prior_args$X <- X
      fit_args$y <- Y
      fit_args$X <- X
    }

    if (!is.null(gamma)) {
      fit_args$gamma_bfi <- gamma
    }

    if (rct) {
      fit_args$RCT_propens <- rep(0.5, nrow(dat))
    }

    Lambda <- do.call(inv.prior.cov, prior_args)
    fit_args$Lambda <- Lambda

    list(
      Lambda = Lambda,
      fit = do.call(MAP.estimation, fit_args)
    )
  }

  run_centers <- function(cc, api, round,
                          gamma = NULL, rct = FALSE) {

    lapply(
      cc, treatment_fit,
      api = api,
      round = round,
      gamma = gamma,
      rct = rct
    )
  }

  compare_local <- function(old, new) {
    for (i in seq_along(old)) {
      expect_equal(
        old[[i]]$Lambda,
        new[[i]]$Lambda,
        tolerance = 1e-8
      )

      compare_results(
        old[[i]]$fit, new[[i]]$fit,
        c("theta_hat", "A_hat",
          "propensity", "for_ATE")
      )
    }
  }

  for (type in c("numeric", "factor")) {

    cc <- lapply(centers, function(dat) {

      dat$treatment <- dat$z

      if (type == "factor") {
        dat$treatment <- factor(
          dat$treatment, levels = c(0, 1)
        )
      }

      dat
    })

    # First round: propensity models.
    old1 <- run_centers(cc, "old", "first")
    new1 <- run_centers(cc, "formula", "first")

    compare_local(old1, new1)

    central_old1 <- aggregate_survival(
      old1, "weibul", "first"
    )
    central_new1 <- aggregate_survival(
      new1, "weibul", "first"
    )

    compare_results(
      central_old1, central_new1,
      c("theta_hat", "A_hat")
    )

    # Second round: using federated propensity coefficients.
    old2 <- run_centers(
      cc, "old", "second",
      gamma = central_old1$theta_hat
    )

    new2 <- run_centers(
      cc, "formula", "second",
      gamma = central_new1$theta_hat
    )

    compare_local(old2, new2)

    central_old2 <- aggregate_survival(
      old2, "unspecified", "second"
    )
    central_new2 <- aggregate_survival(
      new2, "unspecified", "second"
    )

    compare_results(
      central_old2, central_new2,
      c("theta_hat", "A_hat")
    )

    # Constant propensity branch (software test, not an RCT).
    old_rct <- run_centers(
      cc, "old", "second", rct = TRUE
    )

    new_rct <- run_centers(
      cc, "formula", "second", rct = TRUE
    )

    compare_local(old_rct, new_rct)

    for (i in seq_along(new_rct)) {
      expect_equal(
        new_rct[[i]]$fit$propensity,
        rep(0.5, nrow(cc[[i]])),
        tolerance = 1e-8
      )
    }

    central_old_rct <- aggregate_survival(
      old_rct, "unspecified", "second"
    )
    central_new_rct <- aggregate_survival(
      new_rct, "unspecified", "second"
    )

    compare_results(
      central_old_rct, central_new_rct,
      c("theta_hat", "A_hat")
    )
  }
})
