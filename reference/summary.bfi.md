# Summarizing BFI fits

Summary method for an object with class 'bfi' created by the
`MAP.estimation` and `bfi` functions.

## Usage

``` r
# S3 method for class 'bfi'
summary(object,
        cur_mat = FALSE,
        digits = max(3, getOption("digits") - 3),
        ...)

# S3 method for class 'summary.bfi'
print(x,
        digits = x$summary_digits,
        cur_mat = x$summary_cur_mat,
        ...)
```

## Arguments

- object:

  A fitted object of class `"bfi"`.

- x:

  An object of class `"summary.bfi"` returned by `summary.bfi()`.

- cur_mat:

  Logical; if `TRUE`, the minus curvature matrix is also printed when
  the summary object is displayed. The matrix remains available in
  `A_hat` regardless of this option. Default is `FALSE` for
  `summary.bfi()`.

- digits:

  Number of digits used to format the printed results.

- ...:

  Additional arguments, currently unused.

## Details

`summary.bfi()` prepares a structured summary of results returned by
[`MAP.estimation()`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
or [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md).

The resulting object can be stored without producing printed output.
When displayed, `print.summary.bfi()` prints the model family, link
function or baseline hazard, parameter estimates, posterior standard
deviations, and 95% credible intervals.

For local models, the log-posterior value and convergence status are
also displayed.

The minus curvature matrix is retained in `A_hat`. When `cur_mat = TRUE`
is supplied to `summary.bfi()`, this matrix is also printed when the
summary object is displayed.

The print method returns the summary object invisibly.

## Value

`summary.bfi` returns an object of class `"summary.bfi"`, containing the
components described below. `print.summary.bfi()` displays the summary
and returns the supplied object invisibly.

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

- coef_table:

  A matrix containing the reported parameter estimates, posterior
  standard deviations, and 95% credible interval limits.

- summary_digits:

  The number of digits selected for printing the summary.

- summary_cur_mat:

  A logical value indicating whether the curvature matrix should be
  displayed when the summary is printed.

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
#> [1] "bfi_local" "bfi"      

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
#> (Intercept)   0.5084  0.3086 -0.0964   1.1131
#> x1            1.7132  0.1985  1.3241   2.1024
#> x22           3.9894  0.4696  3.0691   4.9098
#> x23           4.6264  0.4396  3.7648   5.4879
#> 
#> Dispersion parameter (sigma2):  1.434 
#>             log Lik Posterior:  -29.44 
#>                   Convergence:  0 
sumfit <- summary(fit, cur_mat = TRUE)
sumfit$estimate
#> (Intercept)          x1         x22         x23 
#>   0.5083817   1.7132473   3.9894438   4.6263612 
sumfit$logLikPost
#> [1] -29.44136
sumfit$dispersion
#>  sigma2 
#> 1.43356 
sumfit$CI
#>                   2.5 %   97.5 %
#> (Intercept) -0.09636638 1.113130
#> x1           1.32412402 2.102371
#> x22          3.06910368 4.909784
#> x23          3.76478601 5.487936
class(sumfit)
#> [1] "summary.bfi"
```
