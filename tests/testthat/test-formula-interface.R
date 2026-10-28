
test_that("prior construction accepts formula without response data", {

  set.seed(101)
  d <- data.frame(
    x = rnorm(100),
    z = factor(rep(c("A", "B"), 50))
  )

  L_old <- inv.prior.cov(
    X = d,
    family = "gaussian",
    lambda = 0.01
  )

  L_one <- inv.prior.cov(
    formula = ~ x + z,
    data = d,
    family = "gaussian",
    lambda = 0.01
  )

  L_two <- inv.prior.cov(
    formula = y ~ x + z,
    data = d,
    family = "gaussian",
    lambda = 0.01
  )

  expect_equal(L_old, L_one)
  expect_equal(L_old, L_two)
})


test_that("formula detects intercept correctly", {

  set.seed(102)

  d <- data.frame(
    y = rnorm(100),
    x = rnorm(100)
  )

  helper <- BFI:::.bfi_formula_inputs

  with_intercept <- helper(
    y ~ x, d, family = "gaussian"
  )

  without_intercept <- helper(
    y ~ x - 1, d, family = "gaussian"
  )

  expect_true(with_intercept$intercept)
  expect_false(without_intercept$intercept)

  expect_identical(names(with_intercept$X), "x")
  expect_identical(names(without_intercept$X), "x")
})



test_that("no-intercept categorical formulas match the legacy interface", {

  set.seed(103)

  d <- data.frame(
    y = rnorm(120),
    x = rnorm(120),
    group = factor(
      rep(c("A", "B", "C"), each = 40)
    )
  )

  f <- y ~ x + group - 1

  X_old <- d[, c("x", "group"), drop = FALSE]

  out <- BFI:::.bfi_formula_inputs(
    formula = f,
    data = d,
    family = "gaussian"
  )

  # Preserve the original predictor columns
  expect_false(out$intercept)
  expect_identical(names(out$X), c("x", "group"))
  expect_identical(out$X$group, d$group)

  # Compare prior precision matrices
  L_old <- inv.prior.cov(
    X = X_old,
    family = "gaussian",
    lambda = 0.01,
    intercept = FALSE
  )

  L_formula <- inv.prior.cov(
    formula = f,
    data = d,
    family = "gaussian",
    lambda = 0.01,
    intercept = FALSE
  )

  expect_equal(
    L_formula, L_old,
    tolerance = 1e-8
  )

  # Compare MAP estimates
  fit_old <- MAP.estimation(
    y = d$y,
    X = X_old,
    family = "gaussian",
    Lambda = L_old,
    intercept = FALSE
  )

  fit_formula <- MAP.estimation(
    formula = f,
    data = d,
    family = "gaussian",
    Lambda = L_formula,
    intercept = FALSE
  )

  expect_equal(
    fit_formula$theta_hat,
    fit_old$theta_hat,
    tolerance = 1e-8
  )

  expect_equal(
    fit_formula$A_hat,
    fit_old$A_hat,
    tolerance = 1e-8
  )
})



test_that("treatment column is preserved for no-intercept formulas", {

  set.seed(104)

  d <- data.frame(
    y = rnorm(120),
    x = rnorm(120),
    group = factor(
      rep(c("A", "B", "C"), each = 40)
    ),
    treatment = rep(0:1, 60)
  )

  f <- y ~ x + group + treatment - 1

  for (type in c("numeric", "factor")) {

    dd <- d

    if (type == "factor") {
      dd$treatment <- factor(dd$treatment)
    }

    for (fam in c("gaussian", "binomial")) {

      out <- BFI:::.bfi_formula_inputs(
        formula = f,
        data = dd,
        family = fam,
        treatment = "treatment"
      )

      expect_false(out$intercept)

      expect_true(
        "treatment" %in% names(out$X)
      )

      expect_identical(
        out$X$treatment,
        dd$treatment
      )

      expect_identical(
        names(out$X),
        c("x", "group", "treatment")
      )

      expect_identical(
        out$X$group,
        dd$group
      )
    }
  }
})


test_that("unsupported formula specifications are rejected", {

  set.seed(105)

  d <- data.frame(
    y = rnorm(100),
    x = rnorm(100),
    z = rnorm(100)
  )

  expect_error(
    inv.prior.cov(
      formula = y ~ x * z,
      data = d,
      family = "gaussian"
    ),
    "Use additive formulas"
  )

  expect_error(
    inv.prior.cov(
      formula = y ~ log(x),
      data = d,
      family = "gaussian"
    ),
    "Use additive formulas"
  )

  expect_error(
    inv.prior.cov(
      formula = y ~ x + offset(z),
      data = d,
      family = "gaussian"
    ),
    "Offset"
  )

  expect_error(
    inv.prior.cov(
      formula = y ~ x,
      family = "gaussian"
    ),
    "data.frame"
  )

  expect_error(
    MAP.estimation(
      formula = ~ x + z,
      data = d,
      family = "gaussian",
      Lambda = diag(3)
    ),
    "response is required"
  )
})



test_that("intercept must agree with the formula", {

  set.seed(20261027)

  d <- data.frame(
    x1 = rnorm(120),
    x2 = rnorm(120)
  )

  for (fam in c("gaussian", "binomial")) {

    if (fam == "gaussian") {
      d$y <- 0.5 * d$x1 -
        0.3 * d$x2 + rnorm(nrow(d))
    } else {
      d$y <- rbinom(
        nrow(d), 1,
        plogis(0.5 * d$x1 - 0.3 * d$x2)
      )
    }

    f <- y ~ x1 + x2 - 1
    X <- d[, c("x1", "x2"), drop = FALSE]

    # 1. The default intercept = TRUE must be rejected.
    expect_error(
      inv.prior.cov(
        formula = f,
        data = d,
        family = fam,
        lambda = 0.1
      ),
      "intercept.*conflicts"
    )

    expect_error(
      MAP.estimation(
        formula = f,
        data = d,
        family = fam,
        Lambda = diag(2)
      ),
      "intercept.*conflicts"
    )

    # 2. Explicit intercept = TRUE must also be rejected.
    expect_error(
      inv.prior.cov(
        formula = f,
        data = d,
        family = fam,
        lambda = 0.1,
        intercept = TRUE
      ),
      "intercept.*conflicts"
    )

    expect_error(
      MAP.estimation(
        formula = f,
        data = d,
        family = fam,
        Lambda = diag(2),
        intercept = TRUE
      ),
      "intercept.*conflicts"
    )

    # 3. Explicit intercept = FALSE must succeed.
    Lambda_formula <- inv.prior.cov(
      formula = f,
      data = d,
      family = fam,
      lambda = 0.1,
      intercept = FALSE
    )

    Lambda_old <- inv.prior.cov(
      X = X,
      family = fam,
      lambda = 0.1,
      intercept = FALSE
    )

    expect_equal(
      Lambda_formula,
      Lambda_old,
      tolerance = 1e-8
    )

    fit_formula <- MAP.estimation(
      formula = f,
      data = d,
      family = fam,
      Lambda = Lambda_formula,
      intercept = FALSE,
      control = list(maxit = 1000)
    )

    fit_old <- MAP.estimation(
      y = d$y,
      X = X,
      family = fam,
      Lambda = Lambda_old,
      intercept = FALSE,
      control = list(maxit = 1000)
    )

    expect_equal(
      fit_formula$theta_hat,
      fit_old$theta_hat,
      tolerance = 1e-8
    )

    expect_equal(
      fit_formula$A_hat,
      fit_old$A_hat,
      tolerance = 1e-8
    )

    expect_identical(fit_formula$intercept, FALSE)

    # 4. The reverse conflict must also be rejected.
    expect_error(
      inv.prior.cov(
        formula = y ~ x1 + x2,
        data = d,
        family = fam,
        intercept = FALSE
      ),
      "intercept.*conflicts"
    )

    expect_error(
      MAP.estimation(
        formula = y ~ x1 + x2,
        data = d,
        family = fam,
        Lambda = diag(3),
        intercept = FALSE
      ),
      "intercept.*conflicts"
    )
  }
})



test_that("prior does not require response but MAP rejects missing response", {

  set.seed(107)

  d <- data.frame(
    y = rnorm(100),
    x = rnorm(100)
  )

  d$y[1] <- NA

  expect_no_error(
    inv.prior.cov(
      formula = y ~ x,
      data = d,
      family = "gaussian",
      lambda = 0.01
    )
  )

  expect_error(
    MAP.estimation(
      formula = y ~ x,
      data = d,
      family = "gaussian",
      Lambda = diag(3)
    )
  )

  d$x[1] <- NA

  expect_error(
    inv.prior.cov(
      formula = y ~ x,
      data = d,
      family = "gaussian",
      lambda = 0.01
    )
  )
})


test_that("survival formula produces time and status", {

  skip_if_not_installed("survival")

  set.seed(108)

  d <- data.frame(
    time = rexp(100) + 0.1,
    status = rep(0:1, 50),
    x = rnorm(100)
  )

  out <- BFI:::.bfi_formula_inputs(
    survival::Surv(time, status) ~ x,
    d,
    family = "survival"
  )

  expect_s3_class(out$y, "data.frame")

  expect_identical(
    names(out$y),
    c("time", "status")
  )

  expect_equal(out$y$time, d$time)
  expect_equal(out$y$status, d$status)

  out_matrix <- BFI:::.bfi_formula_inputs(
    cbind(time, status) ~ x,
    d,
    family = "survival"
  )

  expect_equal(out$y$time, out_matrix$y$time)
  expect_equal(out$y$status, out_matrix$y$status)

  expect_error(
    BFI:::.bfi_formula_inputs(
      survival::Surv(time, status, type = "left") ~ x,
      d,
      family = "survival"
    ),
    "Only right-censored"
  )
})
