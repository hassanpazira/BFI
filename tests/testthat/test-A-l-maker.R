test_that("A.l.maker gives correct Cox curvature for unspecified baseline", {

  time <- c(1, 2, 3, 4, 5, 6)
  status <- c(1, 0, 1, 1, 0, 1)

  y <- data.frame(
    time = time,
    status = status
  )

  Z <- cbind(
    x1 = c(-1.0, 0.3, 0.8, -0.4, 1.2, 0.5),
    x2 = c( 0.2,-0.7, 0.4,  1.1,-0.3, 0.9)
  )

  beta <- c(
    x1 = 0.25,
    x2 = -0.15
  )

  Gamma <- diag(c(0.4, 0.7))
  dimnames(Gamma) <- list(colnames(Z), colnames(Z))

  neg_logpost_cox <- function(beta, time, status, Z, wli, Gamma) {

    eta <- drop(Z %*% beta)

    loglik <- 0

    for (i in which(status == 1)) {

      risk_set <- which(time >= time[i])

      loglik <- loglik +
        wli[i] * (
          eta[i] -
            log(sum(wli[risk_set] * exp(eta[risk_set])))
        )
    }

    prior <- 0.5 * drop(
      t(beta) %*% Gamma %*% beta
    )

    -loglik + prior
  }

  check_curvature <- function(wli) {

    A_analytic <- A.l.maker(
      y = y,
      X = Z,
      Lambda = Gamma,
      family = "survival",
      theta_hat = beta,
      q_l = NULL,
      tps = NULL,
      basehaz = "unspecified",
      wli = wli
    )

    A_numeric <- stats::optimHess(
      par = beta,
      fn = neg_logpost_cox,
      time = time,
      status = status,
      Z = Z,
      wli = wli,
      Gamma = Gamma
    )

    expect_equal(
      unname(A_analytic),
      unname(A_numeric),
      tolerance = 1e-5
    )
  }

  # Unweighted case
  check_curvature(rep(1, length(time)))

  # Non-constant weights
  check_curvature(c(0.8, 1.4, 0.6, 1.7, 0.9, 1.2))
})

