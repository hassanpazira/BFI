
test_that("no-intercept treatment formulas match the legacy federated workflow", {

  set.seed(20261027)

  make_center <- function(n = 180L) {
    x <- rnorm(n)
    group <- factor(
      rep(c("A", "B", "C"), length.out = n),
      levels = c("A", "B", "C")
    )

    treatment <- rbinom(
      n, 1,
      plogis(-0.3 + 0.35 * x + 0.2 * (group == "B"))
    )

    data.frame(
      x = x,
      group = group,
      treatment = treatment,
      yg = 1 + 0.5 * x + 0.4 * (group == "C") +
        0.7 * treatment + rnorm(n),
      yb = rbinom(
        n, 1,
        plogis(-0.5 + 0.3 * x +
                 0.25 * (group == "B") +
                 0.55 * treatment)
      )
    )
  }

  centers <- replicate(3L, make_center(), simplify = FALSE)


  local_fit <- function(d, family, api, round,
                        gamma = NULL, rct = FALSE) {

    response <- if (family == "gaussian") "yg" else "yb"

    f <- stats::reformulate(
      c("x", "group", "treatment"),
      response = response,
      intercept = FALSE
    )

    X_old <- d[, c("x", "group", "treatment"), drop = FALSE]

    X_old$treatment <- d$treatment

    # Both APIs explicitly specify no intercept.
    prior_args <- list(
      family = family,
      lambda = 0.01,
      treatment = "treatment",
      treat_round = round,
      intercept = FALSE
    )

    fit_args <- list(
      family = family,
      treatment = "treatment",
      treat_round = round,
      intercept = FALSE,
      control = list(maxit = 1000)
    )

    if (api == "formula") {

      prior_args$formula <- f
      prior_args$data <- d

      fit_args$formula <- f
      fit_args$data <- d

    } else {

      prior_args$X <- X_old

      fit_args$y <- d[[response]]
      fit_args$X <- X_old
    }

    if (!is.null(gamma)) {
      fit_args$gamma_bfi <- gamma
    }

    if (rct) {
      fit_args$RCT_propens <- rep(0.5, nrow(d))
    }

    prior <- do.call(
      inv.prior.cov,
      prior_args
    )

    fit_args$Lambda <- prior

    list(
      Lambda = prior,
      fit = do.call(MAP.estimation, fit_args)
    )
  }


  combine_fits <- function(results, family, round) {
    args <- list(
      theta_hats = lapply(
        results, function(z) z$fit$theta_hat
      ),
      A_hats = lapply(
        results, function(z) z$fit$A_hat
      ),
      Lambda = results[[1L]]$Lambda,
      family = family,
      treat_round = round
    )

    if (round == "second") {
      args$for_ATE <- lapply(
        results, function(z) z$fit$for_ATE
      )
    }

    do.call(bfi, args)
  }

  compare_local <- function(old, new) {
    for (i in seq_along(old)) {
      expect_equal(
        old[[i]]$Lambda, new[[i]]$Lambda,
        tolerance = 1e-8
      )

      for (field in c(
        "theta_hat", "A_hat", "propensity", "for_ATE"
      )) {
        expect_equal(
          old[[i]]$fit[[field]],
          new[[i]]$fit[[field]],
          tolerance = 1e-8
        )
      }

      expect_identical(
        class(old[[i]]$fit),
        class(new[[i]]$fit)
      )
    }
  }

  compare_central <- function(old, new) {
    for (field in c("theta_hat", "A_hat", "Ave_Treat")) {
      expect_equal(
        old[[field]], new[[field]],
        tolerance = 1e-8
      )
    }
  }

  for (fam in c("gaussian", "binomial")) {
    for (type in c("numeric", "factor")) {

      cc <- lapply(centers, function(d) {
        if (type == "factor") {
          d$treatment <- factor(
            d$treatment, levels = c(0, 1)
          )
        }
        d
      })

      # Round 1: local propensity models.
      old1 <- lapply(
        cc, local_fit, family = fam,
        api = "old", round = "first"
      )
      new1 <- lapply(
        cc, local_fit, family = fam,
        api = "formula", round = "first"
      )

      compare_local(old1, new1)

      central_old1 <- combine_fits(old1, fam, "first")
      central_new1 <- combine_fits(new1, fam, "first")
      compare_central(central_old1, central_new1)

      # Round 2: federated propensity coefficients.
      old2 <- lapply(
        cc, local_fit, family = fam,
        api = "old", round = "second",
        gamma = central_old1$theta_hat
      )
      new2 <- lapply(
        cc, local_fit, family = fam,
        api = "formula", round = "second",
        gamma = central_new1$theta_hat
      )

      compare_local(old2, new2)

      central_old2 <- combine_fits(old2, fam, "second")
      central_new2 <- combine_fits(new2, fam, "second")
      compare_central(central_old2, central_new2)

      expect_false(is.null(central_new2$Ave_Treat))
      expect_true(
        all(is.finite(unlist(central_new2$Ave_Treat)))
      )

      # RCT: known treatment probabilities.
      old_rct <- lapply(
        cc, local_fit, family = fam,
        api = "old", round = "second",
        rct = TRUE
      )
      new_rct <- lapply(
        cc, local_fit, family = fam,
        api = "formula", round = "second",
        rct = TRUE
      )

      compare_local(old_rct, new_rct)

      central_old_rct <- combine_fits(
        old_rct, fam, "second"
      )
      central_new_rct <- combine_fits(
        new_rct, fam, "second"
      )

      compare_central(central_old_rct, central_new_rct)

      for (i in seq_along(new_rct)) {
        expect_equal(
          new_rct[[i]]$fit$propensity,
          rep(0.5, nrow(cc[[i]])),
          tolerance = 1e-8
        )
      }
    }
  }
})
