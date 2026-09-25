test_that("parametric survival curvature adds prior precision exactly once", {

  A_l_maker <- getFromNamespace("A.l.maker", "BFI")

  ## ---------------------------------------------------------
  ## Deterministic survival data
  ## ---------------------------------------------------------

  X <- data.frame(
    x1 = c(-1.20, -0.65, -0.10, 0.35, 0.90, 1.25),
    x2 = c( 0.40, -0.80,  1.10, 0.55, -0.35, 0.20)
  )

  y <- data.frame(
    time   = c(0.40, 0.80, 1.30, 1.90, 2.50, 3.10),
    status = c(1, 0, 1, 1, 0, 1)
  )

  ## Non-unit weights also exercise the weighted likelihood terms.
  wli <- c(1.00, 0.75, 1.30, 0.90, 1.15, 0.80)


  ## ---------------------------------------------------------
  ## Helper: construct two different positive-definite
  ## prior precision matrices.
  ##
  ## Their difference contains both diagonal and off-diagonal
  ## elements, so the test checks the full matrix, not only
  ## its diagonal.
  ## ---------------------------------------------------------

  make_priors <- function(parameter_names) {

    d <- length(parameter_names)

    Lambda1 <- diag(
      seq(0.20, 0.20 + 0.03 * (d - 1), length.out = d)
    )

    Delta <- diag(
      seq(0.04, 0.04 + 0.01 * (d - 1), length.out = d)
    )

    ## Small symmetric off-diagonal changes
    if (d >= 2) {
      Delta[1, 2] <- Delta[2, 1] <- 0.01
    }

    if (d >= 4) {
      Delta[3, 4] <- Delta[4, 3] <- -0.005
    }

    Lambda2 <- Lambda1 + Delta

    dimnames(Lambda1) <-
      list(parameter_names, parameter_names)

    dimnames(Lambda2) <-
      list(parameter_names, parameter_names)

    list(
      Lambda1 = Lambda1,
      Lambda2 = Lambda2,
      Delta = Delta
    )
  }


  ## ---------------------------------------------------------
  ## Model specifications
  ## ---------------------------------------------------------

  models <- list(

    exp = list(

      theta = c(
        x1 = 0.20,
        x2 = -0.15,
        omega_1 = -0.40
      ),

      q_l = NULL
    ),


    weibul = list(

      theta = c(
        x1 = 0.20,
        x2 = -0.15,
        omega_1 = -0.30,
        omega_2 = 0.20
      ),

      q_l = NULL
    ),


    gomp = list(

      theta = c(
        x1 = 0.20,
        x2 = -0.15,
        omega_1 = -0.50,
        omega_2 = -0.80
      ),

      q_l = NULL
    ),


    poly = list(

      theta = c(
        x1 = 0.20,
        x2 = -0.15,
        omega_0 = -0.50,
        omega_1 = 0.10
      ),

      q_l = 1
    )
  )


  ## ---------------------------------------------------------
  ## Core regression test
  ##
  ## At the SAME theta, data, and weights:
  ##
  ## A(Lambda2) - A(Lambda1)
  ##
  ## must equal
  ##
  ## Lambda2 - Lambda1
  ##
  ## exactly up to numerical integration precision.
  ## ---------------------------------------------------------

  for (basehaz in names(models)) {

    spec <- models[[basehaz]]

    priors <- make_priors(
      names(spec$theta)
    )


    A1 <- A_l_maker(
      y = y,
      X = X,
      Lambda = priors$Lambda1,
      family = "survival",
      theta_hat = spec$theta,
      q_l = spec$q_l,
      tps = NULL,
      basehaz = basehaz,
      wli = wli
    )


    A2 <- A_l_maker(
      y = y,
      X = X,
      Lambda = priors$Lambda2,
      family = "survival",
      theta_hat = spec$theta,
      q_l = spec$q_l,
      tps = NULL,
      basehaz = basehaz,
      wli = wli
    )


    observed_prior_contribution <-
      unname(A2 - A1)

    expected_prior_contribution <-
      unname(
        priors$Lambda2 -
          priors$Lambda1
      )


    expect_equal(
      observed_prior_contribution,
      expected_prior_contribution,
      tolerance = 1e-10,
      info = paste(
        "Prior curvature was not added exactly once for",
        basehaz
      )
    )


    ## Additional basic checks
    expect_true(
      isSymmetric(A1),
      info = paste(
        "Curvature matrix is not symmetric for",
        basehaz
      )
    )

    expect_true(
      isSymmetric(A2),
      info = paste(
        "Curvature matrix is not symmetric for",
        basehaz
      )
    )

    expect_false(
      anyNA(A1),
      info = paste(
        "Curvature matrix contains NA for",
        basehaz
      )
    )

    expect_false(
      anyNA(A2),
      info = paste(
        "Curvature matrix contains NA for",
        basehaz
      )
    )
  }
})
