test_that("Local and federated BFI classes work correctly", {

  set.seed(112358)

  X1 <- data.frame(x1 = rnorm(50))
  X2 <- data.frame(x1 = rnorm(60))

  y1 <- 1 + 0.5 * X1$x1 + rnorm(50)
  y2 <- 1 + 0.5 * X2$x1 + rnorm(60)

  Lambda <- inv.prior.cov(
    X1, lambda = 0.01, family = "gaussian"
  )

  fit1 <- MAP.estimation(
    y1, X1, family = "gaussian", Lambda = Lambda
  )

  fit2 <- MAP.estimation(
    y2, X2, family = "gaussian", Lambda = Lambda
  )

  fit_bfi <- bfi(
    theta_hats = list(fit1$theta_hat, fit2$theta_hat),
    A_hats = list(fit1$A_hat, fit2$A_hat),
    Lambda = Lambda,
    family = "gaussian"
  )

  # Class inheritance
  expect_identical(
    class(fit1),
    c("bfi_local", "bfi")
  )

  expect_identical(
    class(fit_bfi),
    c("bfi_federated", "bfi")
  )

  # Shared S3 methods
  for (fit in list(fit1, fit_bfi)) {

    estimates <- coef(fit)
    V <- vcov(fit)

    expect_true(inherits(fit, "bfi"))
    expect_true(is.numeric(estimates))
    expect_equal(
      dim(V),
      c(length(estimates), length(estimates))
    )

    expect_s3_class(
      summary(fit),
      "summary.bfi"
    )

    expect_output(print(fit), "estimates")
    expect_output(
      print(summary(fit)),
      "Coefficients:"
    )
  }

})
