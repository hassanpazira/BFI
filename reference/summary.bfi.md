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
#> (Intercept)   0.5077  0.3072 -0.0945   1.1099
#> x1            1.7133  0.1977  1.3259   2.1007
#> x22           3.9905  0.4676  3.0741   4.9070
#> x23           4.6274  0.4377  3.7695   5.4853
#> 
#> Dispersion parameter (sigma2):  1.421 
#>             log Lik Posterior:  -59.53 
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
#> (Intercept)   0.5077  0.3072 -0.0945   1.1099
#> x1            1.7133  0.1977  1.3259   2.1007
#> x22           3.9905  0.4676  3.0741   4.9070
#> x23           4.6274  0.4377  3.7695   5.4853
#> 
#> Dispersion parameter (sigma2):  1.421 
#>             log Lik Posterior:  -59.53 
#>                   Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>             (Intercept)      x1    x22    x23  sigma2
#> (Intercept)     28.2460  6.0093 7.7402 9.8511  0.0505
#> x1               6.0093 27.0031 0.8318 2.7736  0.1714
#> x22              7.7402  0.8318 7.8402 0.0000  0.3990
#> x23              9.8511  2.7736 0.0000 9.9511  0.4627
#> sigma2           0.0505  0.1714 0.3990 0.4627 20.7106
sumfit$estimate
#> (Intercept)          x1         x22         x23 
#>    0.507736    1.713305    3.990528    4.627400 
sumfit$logLikPost
#> [1] -59.53309
sumfit$dispersion
#>   sigma2 
#> 1.421159 
sumfit$CI
#>                   2.5 %   97.5 %
#> (Intercept) -0.09445457 1.109927
#> x1           1.32586120 2.100748
#> x22          3.07407378 4.906981
#> x23          3.76947073 5.485330
class(sumfit)
#> [1] "summary.bfi"
```
