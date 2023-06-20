# context("Applying MAP estimate")

test_that(desc = "Test MAP.estimation", {
  p <- 3
  theta <- c(1, rep(2, p), 1.5)
  n <- 5
  X <- data.frame(matrix(rnorm(n * p), n, p))
  mu <- gaussian()$linkinv(theta[1] + as.matrix(X) %*% theta[2:4])
  y <- rnorm(n, mu, sd = sqrt(theta[5]))
  Lambda <- inv.prior.cov(X, lambda = 0.05, family = gaussian)
  fit <- MAP.estimation(y, X, family = gaussian, Lambda)
  expect_identical(class(fit), "bfi")
})
