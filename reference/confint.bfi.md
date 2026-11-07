# Approximate posterior credible intervals for BFI models

Computes approximate posterior credible intervals using the
curvature-based covariance matrix. For Gaussian residual variances,
intervals are constructed on the log(sigma2) scale and back-transformed.

## Usage

``` r
# S3 method for class 'bfi'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object inheriting from class "bfi".

- parm:

  Parameter names or positions. Defaults to all parameters.

- level:

  Credible level. Default is 0.95.

- ...:

  Additional arguments, currently unused.

## Value

A matrix of approximate posterior credible intervals.

## Examples

``` r
set.seed(123)
X <- data.frame(x1 = rnorm(50))
y <- 1 + 2 * X$x1 + rnorm(50)
Lambda <- inv.prior.cov(X, lambda = 0.1, family = "gaussian")
fit <- MAP.estimation(y, X, family = "gaussian", Lambda = Lambda)
confint(fit)
#>                 2.5 %   97.5 %
#> (Intercept) 0.8975243 1.394278
#> x1          1.6904896 2.232006
#> sigma2      0.5430045 1.187844
confint(fit, parm = "x1", level = 0.90)
#>        5 %     95 %
#> x1 1.73402 2.188475
```
