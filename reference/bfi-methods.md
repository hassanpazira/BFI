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
#> $theta_hat
#> (Intercept)          x1      sigma2 
#>   0.1742710   0.1539692   0.9242471 
#> 
#> $A_hat
#>              (Intercept)           x1       sigma2
#> (Intercept) 54.108085482  3.835566511  0.001765771
#> x1           3.835566511 57.913345104  0.001548889
#> sigma2       0.001765771  0.001548889 25.009268865
#> 
#> $sd
#> (Intercept)          x1      sigma2 
#>   0.1362670   0.1317142   0.1999629 
#> 
#> $Lambda
#>             (Intercept)   x1 sigma2
#> (Intercept)        0.01 0.00   0.00
#> x1                 0.00 0.01   0.00
#> sigma2             0.00 0.00   0.01
#> 
#> $formula
#> [1] y ~ x1
#> 
#> $names
#> [1] "(Intercept)" "x1"          "sigma2"     
#> 
#> $n
#> [1] 50
#> 
#> $np
#> [1] 2
#> 
#> $treatment
#> NULL
#> 
#> $zero_sample_cov
#> NULL
#> 
#> $refer_cat
#> NULL
#> 
#> $zero_cat
#> NULL
#> 
#> $value
#> [1] 46.08029
#> 
#> $family
#> [1] "gaussian"
#> 
#> $basehaz
#> [1] "weibul"      "exp"         "gomp"        "poly"        "pwexp"      
#> [6] "unspecified"
#> 
#> $intercept
#> [1] TRUE
#> 
#> $convergence
#> [1] 0
#> 
#> $control
#> $control$maxit
#> [1] 100
#> 
#> 
#> attr(,"class")
#> [1] "bfi"
coef(fit)
#> NULL
vcov(fit)
#> Error in UseMethod("vcov"): no applicable method for 'vcov' applied to an object of class "bfi"
```
