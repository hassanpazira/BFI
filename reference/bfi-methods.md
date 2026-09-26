# Methods for objects of class `"bfi"`

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

  An object of class `"bfi"`.

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
#>   -0.129876     0.003659     0.871924  
#> 
coef(fit)
#>  (Intercept)           x1       sigma2 
#> -0.129875691  0.003659374  0.871924188 
vcov(fit)
#>               (Intercept)            x1        sigma2
#> (Intercept)  1.743912e-02 -2.687922e-04  7.963622e-07
#> x1          -2.687922e-04  1.962599e-02 -3.798895e-08
#> sigma2       7.963622e-07 -3.798895e-08  3.039947e-02
```
