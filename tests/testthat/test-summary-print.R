test_that("summary.bfi supports standard S3 printing", {

  set.seed(112358)

  X <- data.frame(x1 = rnorm(50))
  y <- rnorm(50)

  Lambda <- inv.prior.cov(
    X, lambda = 0.01, family = "gaussian"
  )

  fit <- MAP.estimation(
    y, X, family = "gaussian", Lambda = Lambda
  )

  # Assignment should not print
  expect_silent(
    s <- summary(fit, cur_mat = TRUE)
  )

  # Correct class and numerical components
  expect_s3_class(s, "summary.bfi")
  expect_true(is.matrix(s$CI))
  expect_true(is.matrix(s$coef_table))

  # S3 print method
  expect_output(print(s), "Coefficients:")
  expect_output(
    print(s),
    "Minus the Curvature Matrix"
  )

})
