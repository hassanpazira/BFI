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
