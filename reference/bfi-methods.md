# Methods for local and federated BFI objects

Print, coefficient and variance methods for objects returned by
[`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
and [`bfi`](https://hassanpazira.github.io/BFI/reference/BFI.md).

## Usage

``` r
# S3 method for class 'bfi'
coef(object, ...)

# S3 method for class 'bfi'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'bfi'
vcov(object, ...)
```

## Arguments

- ...:

  Further arguments (currently ignored).

- x, object:

  An object inheriting from `"bfi"`, including `"bfi_local"` and
  `"bfi_federated"` objects.

- digits:

  Number of significant digits used for printing.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
vector of the MAP or BFI estimates.
[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the approximate
posterior covariance matrix.

## Details

[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the approximate
posterior covariance matrix corresponding to the parameterization
returned by [`coef()`](https://rdrr.io/r/stats/coef.html). For Gaussian
models, the inverse curvature matrix is defined on the
\\\log(\sigma^2)\\ scale for residual-variance parameters; these rows
and columns are transformed to the original \\\sigma^2\\ scale using the
multivariate delta method.

For weighted models, the returned covariance matrix is a curvature-based
approximation conditional on the weights used in estimation. It does not
automatically account for uncertainty associated with estimated
propensity scores.

## See also

[`summary.bfi`](https://hassanpazira.github.io/BFI/reference/summary.bfi.md)

## Examples

``` r
X <- data.frame(x1 = rnorm(50))
y <- rnorm(50)
Lambda <- inv.prior.cov(X, lambda = 0.01, family = "gaussian")
fit <- MAP.estimation(y, X, family = "gaussian", Lambda = Lambda)
fit
#> 
#> Local MAP estimates (family: gaussian)
#> 
#> (Intercept)           x1       sigma2  
#>   -0.129876     0.003659     0.872124  
#> 
coef(fit)
#>  (Intercept)           x1       sigma2 
#> -0.129875687  0.003659374  0.872123970 
vcov(fit)
#>               (Intercept)            x1        sigma2
#> (Intercept)  1.744312e-02 -2.688538e-04  7.966818e-07
#> x1          -2.688538e-04  1.963049e-02 -3.800053e-08
#> sigma2       7.966818e-07 -3.800053e-08  3.041351e-02
```
