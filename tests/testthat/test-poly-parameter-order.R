test_that("poly bfi is invariant to parameter order across centers", {

  ## ---------------------------------------------------------
  ## Reference parameter order:
  ## deliberately NOT alphabetical
  ## ---------------------------------------------------------

  ref_names <- c(
    "z",
    "a",
    "omega_0",
    "omega_1"
  )

  q_ls <- c(1, 1)


  ## ---------------------------------------------------------
  ## Local MAP estimates
  ## ---------------------------------------------------------

  theta1 <- c(
    z       =  0.25,
    a       = -0.35,
    omega_0 = -0.60,
    omega_1 =  0.08
  )

  theta2 <- c(
    z       =  0.10,
    a       = -0.20,
    omega_0 = -0.55,
    omega_1 =  0.05
  )


  ## ---------------------------------------------------------
  ## Positive-definite local curvature matrices
  ## ---------------------------------------------------------

  A1 <- matrix(
    c(
      12.0,  0.8,  0.4, -0.2,
      0.8, 15.0, -0.3,  0.5,
      0.4, -0.3, 20.0,  0.7,
      -0.2,  0.5,  0.7,  8.0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ref_names, ref_names)
  )

  A2 <- matrix(
    c(
      10.0, -0.6,  0.3,  0.1,
      -0.6, 13.0,  0.2, -0.4,
      0.3,  0.2, 18.0,  0.6,
      0.1, -0.4,  0.6,  9.0
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ref_names, ref_names)
  )


  ## ---------------------------------------------------------
  ## Different prior precision matrices
  ##
  ## Unequal entries are intentional.
  ## This makes incorrect positional handling detectable.
  ## ---------------------------------------------------------

  Lambda1 <- matrix(
    c(
      0.20,  0.01,  0.00,  0.00,
      0.01,  0.30,  0.00,  0.00,
      0.00,  0.00,  0.40,  0.02,
      0.00,  0.00,  0.02,  0.50
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ref_names, ref_names)
  )

  Lambda2 <- matrix(
    c(
      0.25, -0.01,  0.00,  0.00,
      -0.01,  0.35,  0.00,  0.00,
      0.00,  0.00,  0.45,  0.01,
      0.00,  0.00,  0.01,  0.55
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ref_names, ref_names)
  )

  Lambda_global <- matrix(
    c(
      0.15,  0.005, 0.00,  0.00,
      0.005, 0.22,  0.00,  0.00,
      0.00,  0.00,  0.32,  0.01,
      0.00,  0.00,  0.01,  0.42
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ref_names, ref_names)
  )


  ## ---------------------------------------------------------
  ## Helper:
  ## construct a minimal theta_A_poly array for q_l = 1
  ##
  ## For q_l = 1:
  ##   theta_hat is stored in [, 2, 1]
  ##   A_hat     is stored in [, , 3]
  ## ---------------------------------------------------------

  make_poly_array <- function(theta, A) {

    nm <- names(theta)
    d <- length(theta)

    out <- array(
      NA_real_,
      dim = c(d, d, 3),
      dimnames = list(
        nm,
        nm,
        NULL
      )
    )

    out[, 2, 1] <- theta
    out[, , 3] <- A

    out
  }


  ## =========================================================
  ## REFERENCE ANALYSIS
  ##
  ## Both centers use the parameter order of center 1.
  ## =========================================================

  poly1_ref <- make_poly_array(
    theta1,
    A1
  )

  poly2_ref <- make_poly_array(
    theta2,
    A2
  )

  fit_ref <- bfi(
    Lambda = list(
      Lambda1,
      Lambda2,
      Lambda_global
    ),
    family = "survival",
    theta_A_polys = list(
      poly1_ref,
      poly2_ref
    ),
    basehaz = "poly",
    q_ls = q_ls
  )


  ## =========================================================
  ## PERMUTED ANALYSIS
  ##
  ## Center 2 contains EXACTLY the same information,
  ## but its regression coefficients are ordered differently.
  ## =========================================================

  perm_names <- c(
    "a",
    "z",
    "omega_0",
    "omega_1"
  )

  perm <- match(
    perm_names,
    ref_names
  )

  theta2_perm <-
    theta2[perm]

  A2_perm <-
    A2[
      perm,
      perm,
      drop = FALSE
    ]

  Lambda2_perm <-
    Lambda2[
      perm,
      perm,
      drop = FALSE
    ]

  poly2_perm <- make_poly_array(
    theta2_perm,
    A2_perm
  )

  fit_perm <- bfi(
    Lambda = list(
      Lambda1,
      Lambda2_perm,
      Lambda_global
    ),
    family = "survival",
    theta_A_polys = list(
      poly1_ref,
      poly2_perm
    ),
    basehaz = "poly",
    q_ls = q_ls
  )


  ## ---------------------------------------------------------
  ## 1. Reference order must be preserved
  ## ---------------------------------------------------------

  expect_identical(
    rownames(fit_perm$A_hat),
    ref_names
  )

  expect_identical(
    colnames(fit_perm$A_hat),
    ref_names
  )


  ## ---------------------------------------------------------
  ## 2. BFI curvature must be permutation invariant
  ## ---------------------------------------------------------

  expect_equal(
    fit_perm$A_hat,
    fit_ref$A_hat,
    tolerance = 1e-12
  )


  ## ---------------------------------------------------------
  ## 3. BFI point estimates must be permutation invariant
  ## ---------------------------------------------------------

  expect_equal(
    unname(fit_perm$theta_hat),
    unname(fit_ref$theta_hat),
    tolerance = 1e-12
  )


  ## ---------------------------------------------------------
  ## 4. Posterior SDs must be permutation invariant
  ## ---------------------------------------------------------

  expect_equal(
    unname(fit_perm$sd),
    unname(fit_ref$sd),
    tolerance = 1e-12
  )
})



test_that("poly bfi handles different q_ls and scalar q_ls with permuted parameters", {

  ref_names <- c(
    "z",
    "a",
    "omega_0",
    "omega_1"
  )

  ## ---------------------------------------------------------
  ## q_l differs across centers:
  ##
  ## center 1: q_l = 0
  ## center 2: q_l = 1
  ##
  ## Hence q_max = 1.
  ## ---------------------------------------------------------

  q_ls_vec <- c(0, 1)


  ## ---------------------------------------------------------
  ## q = 1 estimates used by BFI
  ## ---------------------------------------------------------

  theta1_q1 <- c(
    z       =  0.24,
    a       = -0.31,
    omega_0 = -0.58,
    omega_1 =  0.07
  )

  theta2_q1 <- c(
    z       =  0.11,
    a       = -0.22,
    omega_0 = -0.52,
    omega_1 =  0.04
  )


  ## center 1 also has a q = 0 fit
  theta1_q0 <- c(
    z       =  0.20,
    a       = -0.28,
    omega_0 = -0.55
  )


  ## ---------------------------------------------------------
  ## Curvature matrices
  ## ---------------------------------------------------------

  A1_q1 <- matrix(
    c(
      12.0,  0.7,  0.3, -0.2,
      0.7, 14.0, -0.2,  0.4,
      0.3, -0.2, 19.0,  0.6,
      -0.2,  0.4,  0.6,  8.5
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ref_names, ref_names)
  )

  A2_q1 <- matrix(
    c(
      10.5, -0.5,  0.2,  0.1,
      -0.5, 13.5,  0.1, -0.3,
      0.2,  0.1, 17.5,  0.5,
      0.1, -0.3,  0.5,  9.5
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(ref_names, ref_names)
  )

  q0_names <- c(
    "z",
    "a",
    "omega_0"
  )

  A1_q0 <- matrix(
    c(
      11.0,  0.5,  0.2,
      0.5, 13.0, -0.1,
      0.2, -0.1, 18.0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(q0_names, q0_names)
  )


  ## ---------------------------------------------------------
  ## Prior precision matrices
  ## ---------------------------------------------------------

  Lambda1 <- diag(
    c(0.20, 0.30, 0.40, 0.50)
  )

  Lambda2 <- diag(
    c(0.25, 0.35, 0.45, 0.55)
  )

  Lambda_global <- diag(
    c(0.15, 0.22, 0.32, 0.42)
  )

  dimnames(Lambda1) <-
    list(ref_names, ref_names)

  dimnames(Lambda2) <-
    list(ref_names, ref_names)

  dimnames(Lambda_global) <-
    list(ref_names, ref_names)


  ## ---------------------------------------------------------
  ## Build theta_A_poly arrays exactly in the structure
  ## expected for max_order = 1:
  ##
  ## slice 1:
  ##   column 1 = q = 0 theta
  ##   column 2 = q = 1 theta
  ##
  ## slice 2 = q = 0 curvature
  ## slice 3 = q = 1 curvature
  ## ---------------------------------------------------------

  poly1 <- array(
    NA_real_,
    dim = c(4, 4, 3),
    dimnames = list(
      ref_names,
      NULL,
      NULL
    )
  )

  poly1[1:3, 1, 1] <- theta1_q0
  poly1[,    2, 1] <- theta1_q1

  poly1[1:3, 1:3, 2] <- A1_q0
  poly1[,    ,    3] <- A1_q1


  poly2 <- array(
    NA_real_,
    dim = c(4, 4, 3),
    dimnames = list(
      ref_names,
      NULL,
      NULL
    )
  )

  ## Center 2 starts at q_l = 1, so q = 0 entries remain NA.
  poly2[, 2, 1] <- theta2_q1
  poly2[, , 3] <- A2_q1


  ## =========================================================
  ## Reference fit
  ## =========================================================

  fit_ref <- bfi(
    Lambda = list(
      Lambda1,
      Lambda2,
      Lambda_global
    ),
    family = "survival",
    theta_A_polys = list(
      poly1,
      poly2
    ),
    basehaz = "poly",
    q_ls = q_ls_vec
  )


  ## =========================================================
  ## Permute regression parameters in center 2
  ## =========================================================

  perm_names <- c(
    "a",
    "z",
    "omega_0",
    "omega_1"
  )

  perm <- match(
    perm_names,
    ref_names
  )

  poly2_perm <- poly2[
    perm,
    ,
    ,
    drop = FALSE
  ]

  ## For curvature slices, the second parameter dimension
  ## must be permuted consistently as well.
  poly2_perm[, , 2] <-
    poly2[
      perm,
      perm,
      2,
      drop = FALSE
    ][, , 1]

  poly2_perm[, , 3] <-
    poly2[
      perm,
      perm,
      3,
      drop = FALSE
    ][, , 1]

  dimnames(poly2_perm)[[1]] <-
    perm_names


  Lambda2_perm <-
    Lambda2[
      perm_names,
      perm_names,
      drop = FALSE
    ]


  fit_perm <- bfi(
    Lambda = list(
      Lambda1,
      Lambda2_perm,
      Lambda_global
    ),
    family = "survival",
    theta_A_polys = list(
      poly1,
      poly2_perm
    ),
    basehaz = "poly",
    q_ls = q_ls_vec
  )


  ## ---------------------------------------------------------
  ## Different local q_l values must not break permutation
  ## invariance.
  ## ---------------------------------------------------------

  expect_identical(
    rownames(fit_perm$A_hat),
    ref_names
  )

  expect_identical(
    colnames(fit_perm$A_hat),
    ref_names
  )

  expect_equal(
    fit_perm$A_hat,
    fit_ref$A_hat,
    tolerance = 1e-12
  )

  expect_equal(
    unname(fit_perm$theta_hat),
    unname(fit_ref$theta_hat),
    tolerance = 1e-12
  )

  expect_equal(
    unname(fit_perm$sd),
    unname(fit_ref$sd),
    tolerance = 1e-12
  )


  ## =========================================================
  ## Scalar q_ls = q_max should produce the same result.
  ## =========================================================

  fit_scalar <- bfi(
    Lambda = list(
      Lambda1,
      Lambda2_perm,
      Lambda_global
    ),
    family = "survival",
    theta_A_polys = list(
      poly1,
      poly2_perm
    ),
    basehaz = "poly",
    q_ls = 1
  )


  expect_equal(
    fit_scalar$A_hat,
    fit_perm$A_hat,
    tolerance = 1e-12
  )

  expect_equal(
    unname(fit_scalar$theta_hat),
    unname(fit_perm$theta_hat),
    tolerance = 1e-12
  )

  expect_equal(
    unname(fit_scalar$sd),
    unname(fit_perm$sd),
    tolerance = 1e-12
  )
})

