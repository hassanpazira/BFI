# Changelog

## BFI v3.2.0

- Added S3 methods
  [`print.bfi()`](https://hassanpazira.github.io/BFI/reference/bfi-methods.md),
  [`coef.bfi()`](https://hassanpazira.github.io/BFI/reference/bfi-methods.md),
  and
  [`vcov.bfi()`](https://hassanpazira.github.io/BFI/reference/bfi-methods.md)
  for objects of class `"bfi"`.

- Improved the handling and presentation of Gaussian residual-variance
  parameters. Posterior covariance matrices returned by
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) are transformed from the
  `log(sigma2)` scale to the original `sigma2` scale using the delta
  method when applicable.

- Improved
  [`summary.bfi()`](https://hassanpazira.github.io/BFI/reference/summary.bfi.md)
  output for center-specific Gaussian residual variances, including
  delta-method standard deviations and back-transformed credible
  intervals.

- Corrected the curvature calculation for parametric survival models
  with exponential, Weibull, Gompertz, and exponentiated-polynomial
  baseline hazards so that the Gaussian-prior precision matrix is added
  once to the likelihood curvature.

- Corrected the weighted curvature calculation for Cox models with an
  unspecified baseline hazard.

- Improved parameter alignment in
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md). Local
  parameter vectors, curvature matrices, and corresponding prior
  precision matrices with the same named parameters may now have
  different parameter orders across centers; the parameter order of the
  first center is used as the reference order.

- Extended name-based parameter alignment to exponentiated-polynomial
  survival models.

- Corrected handling of scalar `q_ls` values for
  exponentiated-polynomial survival models.

## BFI v3.1.0

CRAN release: 2025-05-24

- The `MAP.estimation` function has been revised, and corresponding
  updates have been made to its documentation and examples.

## BFI v3.0.1

CRAN release: 2025-05-19

- In this release, the package has been updated to support both
  observational and randomized trial data for estimating treatment
  effects using Bayesian Federated ‘Causal’ Inference.

## BFI v2.1.0

- In this release, the package has been updated to support a
  (semi-parametric) Cox model with an ‘unspecified’ baseline hazard,
  maximizing partial log-likelihood.
- The documentation and examples have been updated.

## BFI v2.0.1

CRAN release: 2024-07-04

- This is a bugfix release to resolve a minor bug in the `optim`
  function entries related to the `gaussian` family.

## BFI v2.0.0

- The package now supports survival data analysis, adding comprehensive
  tools for time-to-event modeling alongside existing GLM
  functionalities.
- The documentation and examples have been updated.

## BFI v1.1.4

CRAN release: 2024-04-27

- The package has been updated to address a few NOTEs identified during
  the CRAN submission.

## BFI v1.1.1

- The package was prepared for submission to CRAN with minor
  modifications.

## BFI v1.1.0

- The functions
  [`inv.prior.cov()`](https://hassanpazira.github.io/BFI/reference/inv.prior.cov.md)
  and [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md)
  were adapted when there is a center specific variable.
- The function
  [`summary()`](https://hassanpazira.github.io/BFI/reference/summary.bfi.md)
  was updated to be used in the case of stratification.
- The manual pdf was updated in the case of center specific variable.
- Henceforth the `BFI` package can be called from `Python`.
- A vignette was also added to the package for calling `BFI` from
  `Python`.

## BFI v1.0.0

- In this release, the package website was built using the pkgdown
  package.
- Most of the functions were adapted with extra arguments in the case of
  stratification.
- The manual pdf was also updated.

## BFI v0.6.4

- The manual pdf was updated.
- This is also a bugfix release to resolve one minor bug.

## BFI v0.6.3

- In this release, the package was updated so that the outputs of
  [`inv.prior.cov()`](https://hassanpazira.github.io/BFI/reference/inv.prior.cov.md)
  and
  [`MAP.estimation()`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  (for different centers) have the same dimensions when
  `intercept=FALSE` to be used in
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md).
- Moreover, one argument in
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md), i.e.,
  `const_var`, is added to the package to handle the constant variables.
- This is a bugfix release to resolve one minor bug as well.

## BFI v0.5.2

- In this release, the summary function (a S3 method for class `bfi`)
  was added to the package.
- A package pdf manual was created for the package.

## BFI v0.4.2

- In this release, the package was updated with several arguments, e.g.,
  an intercept should be fitted or not.
- A vignette is added to the package as well.

## BFI v0.3.2

- In this release, the package can carry out the stratified analysis.

## BFI v0.2.2

- This is a bugfix release to resolve one minor bug, and add two
  functions related to building Gamma matrix.
- Moreover, the package now handles the categorical covariates with more
  than two levels.
