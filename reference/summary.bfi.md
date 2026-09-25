# Summarizing BFI Fits

Summary method for an object with class 'bfi' created by the
`MAP.estimation` and `bfi` functions.

## Usage

``` r
# S3 method for class 'bfi'
summary(object,
        cur_mat = FALSE,
        digits = max(3, getOption("digits") - 3),
        ...)
```

## Arguments

- object:

  fitted `bfi` object.

- cur_mat:

  logical; if `TRUE`, minus the curvature matrix around the estimated
  parameters is returned and printed. Default is `FALSE`.

- digits:

  significant digits in printout.

- ...:

  additional arguments affecting the summary produced.

## Details

`summary.bfi()` gives information about the MAP estimates of parameters
of the model. It can be used for the `bfi` objects built by the
`MAP.estimation` and `bfi` functions.

The output of the summary method shows the details of the model, i.e.
formula, family and link function used to specify the generalized linear
model, followed by information about the estimates, standard deviations
and credible intervals. Information about the log-likelihood posterior
and convergence status are also provided.

By default, `summary.bfi` function does not return (minus) the curvature
matrix, but the user can use `cur_mat = TRUE` to print it.

## Value

`summary.bfi` returns an object of class `summary.bfi`, a list with the
following components:

- theta_hat:

  the component from `object`. The last element of this vector is the
  estimate of the dispersion parameter (sigma2) if
  `family = "gaussian"`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  and [`bfi`](https://hassanpazira.github.io/BFI/reference/BFI.md)
  functions.

- A_hat:

  the component from `object`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  and [`bfi`](https://hassanpazira.github.io/BFI/reference/BFI.md)
  functions.

- sd:

  the posterior standard deviations corresponding to the reported
  parameter estimates. For Gaussian models with center-specific residual
  variances, the standard deviations corresponding to `sigma2` are
  transformed from the \\\log(\sigma^2)\\ scale to the original
  \\\sigma^2\\ scale using the delta method. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  and [`bfi`](https://hassanpazira.github.io/BFI/reference/BFI.md)
  functions.

- Lambda:

  the component from `object`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  function.

- formula:

  the component from `object`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  function.

- n:

  the component from `object`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  function.

- np:

  the component from `object`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  function.

- family:

  the component from `object`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  function.

- intercept:

  the component from `object`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  function.

- convergence:

  the component from `object`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  function.

- control:

  the component from `object`. See the
  [`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
  function.

- stratified:

  the component from `object`. See the
  [`bfi`](https://hassanpazira.github.io/BFI/reference/BFI.md) function.

- estimate:

  the estimated regression coefficients. For `gaussian` models the
  common dispersion parameter is excluded, but center-specific residual
  variances obtained with `strat_par = 2` or `strat_par = c(1, 2)` are
  included.

- logLikPost:

  the value of the log-likelihood posterior density evaluated at
  estimates (`theta_hat`).

- link:

  the link function only for GLMs, not for the survival family. By
  default the `gaussian` family with `identity` link function and the
  `binomial` family with `logit` link function are used.

- dispersion:

  the estimated residual variance `sigma2` for Gaussian models with a
  common residual variance. For Gaussian models with center-specific
  residual variances, this component is `NULL`. For the `binomial`
  family, the dispersion parameter is fixed at `1`. For the `survival`
  family, this component is `NULL`.

- CI:

  a 95`%` credible interval for the parameters. For center-specific
  residual variances the interval is computed on the \\\log(\sigma^2)\\
  scale and back-transformed, so it is always positive and not symmetric
  around the estimate. The standard deviation printed for these
  parameters is obtained by the delta method on the \\\sigma^2\\ scale.

## Author

Hassan Pazira  
Maintainer: Hassan Pazira <h.pazira@arq.org>

## See also

[`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
and [`bfi`](https://hassanpazira.github.io/BFI/reference/BFI.md)

## Examples

``` r
#-------------
# y ~ Gaussian
#-------------
# model assumption:
theta <- c(1, 2, 3, 4, 1.5)  # coefficients and sigma2 = 1.5

#----------------
# Data Simulation
#----------------
n      <- 40
X      <- data.frame(x1=rnorm(n),                     # continuous variable
                     x2=sample(1:3, n, replace=TRUE)) # categorical variable
Xx2_1  <- ifelse(X$x2 == '2', 1, 0)
Xx2_2  <- ifelse(X$x2 == '3', 1, 0)
X$x2   <- as.factor(X$x2)
eta    <- theta[1] + theta[2] * X$x1 + theta[3] * Xx2_1 + theta[4] * Xx2_2
mu     <- gaussian()$linkinv(eta)
y      <- rnorm(n, mu, sd = sqrt(theta[5]))

#----------------
# MAP estimations
#----------------
Lambda <- inv.prior.cov(X, lambda = c(0.1, 0.5), family = "gaussian")
fit    <- MAP.estimation(y, X, family = "gaussian", Lambda)
class(fit)
#> [1] "bfi"

#-------------------------
# Summary of MAP estimates
#-------------------------
summary(fit)
#> 
#> Summary of the local model:
#> 
#>    Formula: y ~ x1 + x2 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.6469  0.2409  0.1747   1.1191
#> x1            1.5616  0.1490  1.2696   1.8537
#> x22           3.5545  0.3544  2.8598   4.2491
#> x23           4.2078  0.3406  3.5401   4.8755
#> 
#> Dispersion parameter (sigma2):  0.8246 
#>             log Lik Posterior:  -36.43 
#>                   Convergence:  0 
sumfit <- summary(fit, cur_mat = TRUE)
#> 
#> Summary of the local model:
#> 
#>    Formula: y ~ x1 + x2 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.6469  0.2409  0.1747   1.1191
#> x1            1.5616  0.1490  1.2696   1.8537
#> x22           3.5545  0.3544  2.8598   4.2491
#> x23           4.2078  0.3406  3.5401   4.8755
#> 
#> Dispersion parameter (sigma2):  0.8246 
#>             log Lik Posterior:  -36.43 
#>                   Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>             (Intercept)      x1     x22     x23  sigma2
#> (Intercept)     48.6061 -4.4530 14.5518 16.9771  0.0653
#> x1              -4.4530 45.4708 -1.5784 -1.9001  0.1564
#> x22             14.5518 -1.5784 14.6518  0.0000  0.3558
#> x23             16.9771 -1.9001  0.0000 17.0771  0.4212
#> sigma2           0.0653  0.1564  0.3558  0.4212 20.4120
sumfit$estimate
#> (Intercept)          x1         x22         x23 
#>   0.6468995   1.5616445   3.5544857   4.2077952 
sumfit$logLikPost
#> [1] -36.43131
sumfit$dispersion
#>    sigma2 
#> 0.8246386 
sumfit$CI
#>                 2.5 %   97.5 %
#> (Intercept) 0.1747179 1.119081
#> x1          1.2695673 1.853722
#> x22         2.8598380 4.249133
#> x23         3.5401350 4.875455
class(sumfit)
#> [1] "summary.bfi"
```
