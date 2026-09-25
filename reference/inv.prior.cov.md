# Creates an inverse covariance matrix for a Gaussian prior

`inv.prior.cov` constructs a diagonal inverse covariance matrix for the
Gaussian prior distribution based on the design matrix of covariates.
This construction accounts for the number of regression parameters,
especially when dealing with categorical covariates. For a linear model,
it also includes an additional row and column to represent the variance
of the measurement error. In the case of a survival model, it considers
the parameters of the baseline hazard function as well.

## Usage

``` r
inv.prior.cov(X,
              lambda = 1,
              L = 2,
              family = c("gaussian", "binomial", "survival"),
              treatment = NULL,
              treat_round = NULL,
              intercept = TRUE,
              stratified = FALSE,
              strat_par = NULL,
              center_spec = NULL,
              basehaz = c("weibul", "exp", "gomp", "poly", "pwexp", "unspecified"),
              max_order = 2,
              n_intervals = 4)
```

## Arguments

- X:

  design matrix of dimension \\n \times p\\, where \\n\\ is the number
  of samples observed, and \\p\\ is the number of predictors/variables
  so excluding the intercept.

- lambda:

  the vector used as the diagonal of the (inverse covariance) matrix
  that will be created by `inv.prior.cov()`. The length of the vector
  depends on the number of columns of `X`, type of the covariates
  (continuous/dichotomous or categorical), `family`, whether an
  `intercept` is included in the model, and whether `stratified`
  analysis is desired. When `stratified = FALSE`, `lambda` could be a
  single positive number (if all values in the vector are equal), a
  vector of two elements (the first is used for regression parameters
  including “intercept” and the second for the residual standard
  deviation “sigma” in the gaussian family or for the baseline hazard
  parameters in the survival case), or a vector of length equal to the
  number of model parameters. However, the length of `lambda` is
  different when `stratified = TRUE`, see ‘Details’ for more
  information. Default is `lambda = 1`, which means all model parameters
  are set to one.

- L:

  the number of centers. This argument is used only when
  `stratified = TRUE`. Default is `L = 2`. See ‘Details’ and ‘Examples’.

- family:

  a description of the error distribution. This is a character string
  naming a family of the model. In the current version of the package,
  the family of model can be `"gaussian"` (with `identity` link
  function), `"binomial"` (with `logit` link function), or `"survival"`.
  Can be abbreviated. By default the `gaussian` family is used. In case
  of a linear regression model, `family = "gaussian"`, there is an extra
  model parameter for the variance of measurement error. While in the
  case of survival model, `family = "survival"`, the number of the model
  parameters depend on the choice of baseline hazard functions, see
  ‘Details’ for more information.

- treatment:

  a character string representing the name or place of the binary
  covariate, respectively. This covariate indicates whether the patient
  got the new treatment (\\z\_{\ell i}=1\\) or the placebo/standard
  treatment (\\z\_{\ell i}=0\\). For both the first and second rounds,
  it should not be '`NULL`'. See ‘Details’.

- treat_round:

  a character string representing the `'first'` or `'second'` round of
  estimating treatment effects. In the first round,
  `treat_round = 'first'`, the local estimates of the coefficients
  (\\\gamma\_{\ell}\\) is estimated. In the second round,
  `treat_round = 'second'`, the `propensity` scores and the statistical
  summaries (`for_ATE`) are calculated.

- intercept:

  logical flag for having an intercept. It is not used when
  `family = "survival"`. By changing the `intercept` the dimension of
  the inverse covariance matrix changes. If `intercept = TRUE` (the
  default), the output matrix created by `inv.prior.cov()` has one row
  and one column related to `intercept`, while if `intercept = FALSE`,
  the resulting matrix does not have the row and column called
  `intercept`.

- stratified:

  logical flag for performing the stratified analysis. If
  `stratified = TRUE`, the parameter(s) selected in the `strat_par`
  argument are allowed to be different across centers to deal with
  heterogeneity across centers. This argument should only be used when
  designing the inverse covariance matrix for the (fictive) combined
  data, i.e., the last matrix for the `Lambda` argument in
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md). If
  `inv.prior.cov()` is used for the analysis in the local centers (to
  build the \\L\\ first matrices for the `Lambda` argument in
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md)), this
  argument should be `FALSE`, even if the BFI analysis is stratified.
  Default is `stratified = FALSE`. See ‘Details’ and ‘Examples’.

- strat_par:

  a integer vector for indicating the stratification parameter(s). It
  can be used to deal with heterogeneity due to center-specific
  parameters. For the `"binomial"` and `"gaussian"` families it is a
  one- or two-element integer vector so that the values \\1\\ and/or
  \\2\\ are/is used to indicate that the “intercept” and/or “sigma2” are
  allowed to vary, respectively. For the `"binomial"` family the length
  of the vector should be one which refers to “intercept”, and the value
  of this element should be \\1\\ (to handel heterogeneity across
  outcome means). For `"gaussian"` this vector can be \\1\\ for
  indicating the “intercept” only (handeling heterogeneity across
  outcome means), \\2\\ for indicating the “sigma2” only (handeling
  heterogeneity due to nuisance parameter), and c(\\1\\, \\2\\) for both
  “intercept” and “sigma2”. When `family = "survival"`, this vector can
  contain any combination from 1 to the maximum number of parameters of
  the baseline function, i.e., \\1\\ for `"exp"`, \\2\\ for `"weibul"`
  and `"gomp"`, `max_order + 1` for `"poly"`, and `n_intervals` for
  `"pwexp"`. This argument is only used when `stratified = TRUE`.
  Default is `strat_par = NULL`. If `stratified = TRUE`, `strat_par` can
  not be `NULL` except when `center_spec` is not `NULL` for handeling
  heterogeneity due to clustering and missing covariates. See ‘Details’
  and ‘Examples’.

- center_spec:

  a vector of \\L\\ elements for representing the center specific
  variable. This argument is used only when `stratified = TRUE` and
  `strat_par = NULL`. Each element represents a specific feature of the
  corresponding center. There must be only one specific value or
  attribute for each center. This vector could be a numeric,
  characteristic or factor vector. Note that, the order of the centers
  in the vector `center_spec` must be the same as in the list of the
  argument `theta_hats` in the function
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md). The
  used data type in the argument `center_spec` must be categorical.
  Default is `center_spec = NULL`. See also ‘Details’ and ‘Examples’.

- basehaz:

  a character string representing one of the available baseline hazard
  functions; `exponential` (`"exp"`), `Weibull` (`"weibul"`, the
  default), `Gompertz` (`"gomp"`), `exponentiated polynomial`
  (`"poly"`), `piecewise constant exponential` (`"pwexp"`), and
  `unspecified baseline hazard` (`"unspecified"`). Can be abbreviated.
  It is only used when `family = "survival"`.

- max_order:

  an integer representing the maximum value of `q_l`, which is the
  order/degree minus 1 of the exponentiated polynomial baseline hazard
  function. This argument is only used when `family = "survival"` and
  `basehaz = "poly"`. Default is 2.

- n_intervals:

  an integer representing the number of intervals in the piecewise
  exponential baseline hazard function. This argument is only used when
  `family = "survival"` and `basehaz = "pwexp"`. Default is 4.

## Details

`inv.prior.cov` creates a diagonal matrix with the vector `lambda` as
its diagonal. The argument `stratified = TRUE` should only be used to
construct a matrix for the prior density in case of stratification in
the fictive combined data. Never be used for the construction of the
matrix for analysis in the centers.

When `stratified = FALSE`, the length of the vector `lambda` depends on
the covariate matrix `X`, `family`, `basehaz`, and whether an
“intercept” is included in the model. For example, if the design matrix
`X` has `p` columns with continuous or dichotomous covariates,
`family = gaussian`, and `intercept = TRUE`, then `lambda` should have
\\p+2\\ elements. In this case, if in `X` there is a categorical
covariate with \\q\>2\\ categories, then the length of `lambda`
increases with \\q-2\\.

All values of `lambda` should be positive as they represent the inverse
of the variance of the Gaussian prior. This argument is considered as
the inverse of the variance of the prior distribution for: \\(\beta_0,
\boldsymbol{\beta})\\ if `family = "binomial"` and `intercept = TRUE`;
\\(\beta_0, \boldsymbol{\beta},\sigma)\\ if `family = "gaussian"` and
`intercept = TRUE`; and \\( \boldsymbol{\beta},\boldsymbol{\omega})\\ if
`family = "survival"`.

If all values in the vector `lambda` equal, one value is enough to be
given as entry. If `lambda` is a scalar, the function `inv.prior.cov`
sets each value at the diagonal equal to `lambda`. When `lambda` is two
dimensional: if `family = "binomial"`, the first and second values are
used for the inverse of the variance of the prior distribution for the
intercept (\\\beta_0\\) and regression parameters
(\\\boldsymbol{\beta}\\), respectively; If `family = "gaussian"`, the
first and second values are used for the inverse of the variance of the
prior distribution for the regression parameters including the intercept
(\\\beta_0, \boldsymbol{\beta}\\) and variance of the measurement error
(\\ \sigma\\), respectively; If `family = "survival"`, the first and
second values are used for the inverse of the variance of the prior
distribution for the regression parameters (\\\boldsymbol{\beta}\\) and
baseline hazard parameters (\\ \omega\\), respectively. But if
`stratified = TRUE` the length of the vector `lambda` must be equal to
the number of parameters in the combined model.

If `intercept = FALSE`, for the `binomial` family the stratified
analysis is not possible therefore `stratified` can not be `TRUE`.

If `stratified = FALSE`, both `strat_par` and `center_spec` must be
`NULL` (the defaults), while if `stratified = TRUE` only one of the two
must be `NULL`.

If `stratified = TRUE` and `family = "survival"`, `strat_par = 1` refers
to \\\omega_0\\ when `basehaz = "poly"`, and to \\\omega_1\\ for other
baseline hazards.

The output of `inv.prior.cov()` can be used in the main functions
[`MAP.estimation()`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
and [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md).

## Value

`inv.prior.cov` returns a diagonal matrix. The dimension of the matrix
depends on the number of columns of `X`, type of the covariates
(continuous/dichotomous or categorical), `intercept`, `family`, and
`basehaz`.

## References

Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker M.A.
(2026). *Bayesian federated inference for survival models*, *Journal of
Applied Statistics*, 53(2): 203-223.
\<https://doi.org/10.1080/02664763.2025.2511932\>

Jonker M.A., Pazira H. and Coolen A.C.C. (2025). *Bayesian Federated
Inference for regression models based on non-shared medical center
data*, *Research Synthesis Methods*, 16(2): 383-423.
\<https://doi.org/10.1017/rsm.2025.6\>

Jonker M.A., Pazira H. and Coolen A.C.C. (2024). *Bayesian federated
inference for estimating statistical models based on non-shared
multicenter data sets*, *Statistics in Medicine*, 43(12): 2421-2438.
\<https://doi.org/10.1002/sim.10072\>

## Author

Hassan Pazira  
Maintainer: Hassan Pazira \<h.pazira@arq.org\>

## See also

[`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)

## Examples

``` r
#----------------
# Data Simulation
#----------------
X <- data.frame(x1=rnorm(50),                     # standard normal variable
                x2=sample(0:2, 50, replace=TRUE), # categorical variable
                x3=sample(0:1, 50, replace=TRUE)) # dichotomous variable
X$x2 <- as.factor(X$x2)
X$x3 <- as.factor(X$x3)

# The (inverse) variance value (lambda=0.05) is assumed to be
# the same for Gaussian prior of all parameters (for non-stratified)

#-------------------------------------------------
# Inverse Covariance Matrix for the Gaussian prior
#-------------------------------------------------
# y ~ Binomial with 'intercept'
inv.prior.cov(X, lambda = 0.05, family = 'binomial')
#>             (Intercept)   x1  x21  x22  x31
#> (Intercept)        0.05 0.00 0.00 0.00 0.00
#> x1                 0.00 0.05 0.00 0.00 0.00
#> x21                0.00 0.00 0.05 0.00 0.00
#> x22                0.00 0.00 0.00 0.05 0.00
#> x31                0.00 0.00 0.00 0.00 0.05
# returns a 5-by-5 matrix

# y ~ Binomial without 'intercept'
inv.prior.cov(X, lambda = 0.05, family = "binomial", intercept = FALSE)
#>       x1  x21  x22  x31
#> x1  0.05 0.00 0.00 0.00
#> x21 0.00 0.05 0.00 0.00
#> x22 0.00 0.00 0.05 0.00
#> x31 0.00 0.00 0.00 0.05
# a 4-by-4 matrix

# y ~ Gaussian with 'intercept'
inv.prior.cov(X, lambda = 0.05, family = 'gaussian')
#>             (Intercept)   x1  x21  x22  x31 sigma2
#> (Intercept)        0.05 0.00 0.00 0.00 0.00   0.00
#> x1                 0.00 0.05 0.00 0.00 0.00   0.00
#> x21                0.00 0.00 0.05 0.00 0.00   0.00
#> x22                0.00 0.00 0.00 0.05 0.00   0.00
#> x31                0.00 0.00 0.00 0.00 0.05   0.00
#> sigma2             0.00 0.00 0.00 0.00 0.00   0.05
# returns a 6-by-6 matrix

# Survival family with 'weibul' baseline hazard
inv.prior.cov(X, lambda = c(0.05, 0.1), family = 'survival')
#>           x1  x21  x22  x31 omega_1 omega_2
#> x1      0.05 0.00 0.00 0.00     0.0     0.0
#> x21     0.00 0.05 0.00 0.00     0.0     0.0
#> x22     0.00 0.00 0.05 0.00     0.0     0.0
#> x31     0.00 0.00 0.00 0.05     0.0     0.0
#> omega_1 0.00 0.00 0.00 0.00     0.1     0.0
#> omega_2 0.00 0.00 0.00 0.00     0.0     0.1
# returns a 6-by-6 matrix

# Survival family with 'pwexp' baseline hazard (4 intervals)
inv.prior.cov(X, lambda = 0.05, family = 'survival', basehaz = "pwexp")
#>           x1  x21  x22  x31 omega_1 omega_2 omega_3 omega_4
#> x1      0.05 0.00 0.00 0.00    0.00    0.00    0.00    0.00
#> x21     0.00 0.05 0.00 0.00    0.00    0.00    0.00    0.00
#> x22     0.00 0.00 0.05 0.00    0.00    0.00    0.00    0.00
#> x31     0.00 0.00 0.00 0.05    0.00    0.00    0.00    0.00
#> omega_1 0.00 0.00 0.00 0.00    0.05    0.00    0.00    0.00
#> omega_2 0.00 0.00 0.00 0.00    0.00    0.05    0.00    0.00
#> omega_3 0.00 0.00 0.00 0.00    0.00    0.00    0.05    0.00
#> omega_4 0.00 0.00 0.00 0.00    0.00    0.00    0.00    0.05
# returns a 8-by-8 matrix

# Survival family with 'poly' baseline hazard
inv.prior.cov(X, lambda = c(0.05, 0.1), family = 'survival', basehaz = "poly")
#>           x1  x21  x22  x31 omega_0 omega_1 omega_2
#> x1      0.05 0.00 0.00 0.00     0.0     0.0     0.0
#> x21     0.00 0.05 0.00 0.00     0.0     0.0     0.0
#> x22     0.00 0.00 0.05 0.00     0.0     0.0     0.0
#> x31     0.00 0.00 0.00 0.05     0.0     0.0     0.0
#> omega_0 0.00 0.00 0.00 0.00     0.1     0.0     0.0
#> omega_1 0.00 0.00 0.00 0.00     0.0     0.1     0.0
#> omega_2 0.00 0.00 0.00 0.00     0.0     0.0     0.1
# returns a 7-by-7 matrix

#--------------------
# Stratified analysis
#--------------------
# y ~ Binomial when 'intercept' varies across 3 centers:
inv.prior.cov(X, lambda = c(.2, 1), family = 'binomial', stratified = TRUE,
              strat_par = 1, L = 3)
#>                  (Intercept)_loc1 (Intercept)_loc2 (Intercept)_loc3 x1 x21 x22
#> (Intercept)_loc1              0.2              0.0              0.0  0   0   0
#> (Intercept)_loc2              0.0              0.2              0.0  0   0   0
#> (Intercept)_loc3              0.0              0.0              0.2  0   0   0
#> x1                            0.0              0.0              0.0  1   0   0
#> x21                           0.0              0.0              0.0  0   1   0
#> x22                           0.0              0.0              0.0  0   0   1
#> x31                           0.0              0.0              0.0  0   0   0
#>                  x31
#> (Intercept)_loc1   0
#> (Intercept)_loc2   0
#> (Intercept)_loc3   0
#> x1                 0
#> x21                0
#> x22                0
#> x31                1

# y ~ Gaussian when 'intercept' and 'sigma2' vary across 2 centers; y ~ Gaussian
inv.prior.cov(X, lambda = c(1, 2, 3), family = "gaussian", stratified = TRUE,
              strat_par = c(1, 2))
#>                  (Intercept)_loc1 (Intercept)_loc2 x1 x21 x22 x31 sigma2_loc1
#> (Intercept)_loc1                1                0  0   0   0   0           0
#> (Intercept)_loc2                0                1  0   0   0   0           0
#> x1                              0                0  2   0   0   0           0
#> x21                             0                0  0   2   0   0           0
#> x22                             0                0  0   0   2   0           0
#> x31                             0                0  0   0   0   2           0
#> sigma2_loc1                     0                0  0   0   0   0           3
#> sigma2_loc2                     0                0  0   0   0   0           0
#>                  sigma2_loc2
#> (Intercept)_loc1           0
#> (Intercept)_loc2           0
#> x1                         0
#> x21                        0
#> x22                        0
#> x31                        0
#> sigma2_loc1                0
#> sigma2_loc2                3

# y ~ Gaussian when 'sigma2' varies across 2 centers (with 'intercept')
inv.prior.cov(X, lambda = c(1, 2, 3), family='gaussian', stratified = TRUE,
              strat_par = 2)
#>             (Intercept) x1 x21 x22 x31 sigma2_loc1 sigma2_loc2
#> (Intercept)           1  0   0   0   0           0           0
#> x1                    0  2   0   0   0           0           0
#> x21                   0  0   2   0   0           0           0
#> x22                   0  0   0   2   0           0           0
#> x31                   0  0   0   0   2           0           0
#> sigma2_loc1           0  0   0   0   0           3           0
#> sigma2_loc2           0  0   0   0   0           0           3

# y ~ Gaussian when 'sigma2' varies across 2 centers (without 'intercept')
inv.prior.cov(X, lambda = c(2, 3), family = "gaussian", intercept = FALSE,
              stratified=TRUE, strat_par = 2)
#>             x1 x21 x22 x31 sigma2_loc1 sigma2_loc2
#> x1           2   0   0   0           0           0
#> x21          0   2   0   0           0           0
#> x22          0   0   2   0           0           0
#> x31          0   0   0   2           0           0
#> sigma2_loc1  0   0   0   0           3           0
#> sigma2_loc2  0   0   0   0           0           3

#--------------------------
# Center specific covariate
#--------------------------
# center specific covariate has K = 2 categories across 4 centers; y ~ Binomial
inv.prior.cov(X, lambda = c(0.1:2), family = 'binomial', stratified = TRUE,
              center_spec = c("Iran","Netherlands","Netherlands","Iran"), L=4)
#>                         (Intercept)_Iran (Intercept)_Netherlands  x1 x21 x22
#> (Intercept)_Iran                     0.1                     0.0 0.0 0.0 0.0
#> (Intercept)_Netherlands              0.0                     0.1 0.0 0.0 0.0
#> x1                                   0.0                     0.0 1.1 0.0 0.0
#> x21                                  0.0                     0.0 0.0 1.1 0.0
#> x22                                  0.0                     0.0 0.0 0.0 1.1
#> x31                                  0.0                     0.0 0.0 0.0 0.0
#>                         x31
#> (Intercept)_Iran        0.0
#> (Intercept)_Netherlands 0.0
#> x1                      0.0
#> x21                     0.0
#> x22                     0.0
#> x31                     1.1

# center specific covariate has K = 3 categories across 5 centers; y ~ Gaussian
inv.prior.cov(X, lambda = c(0.5:3), family = 'gaussian', stratified = TRUE,
              center_spec = c("Medium","Big","Small","Big","Small"), L = 5)
#>                    (Intercept)_Big (Intercept)_Medium (Intercept)_Small  x1 x21
#> (Intercept)_Big                0.5                0.0               0.0 0.0 0.0
#> (Intercept)_Medium             0.0                0.5               0.0 0.0 0.0
#> (Intercept)_Small              0.0                0.0               0.5 0.0 0.0
#> x1                             0.0                0.0               0.0 1.5 0.0
#> x21                            0.0                0.0               0.0 0.0 1.5
#> x22                            0.0                0.0               0.0 0.0 0.0
#> x31                            0.0                0.0               0.0 0.0 0.0
#> sigma2                         0.0                0.0               0.0 0.0 0.0
#>                    x22 x31 sigma2
#> (Intercept)_Big    0.0 0.0    0.0
#> (Intercept)_Medium 0.0 0.0    0.0
#> (Intercept)_Small  0.0 0.0    0.0
#> x1                 0.0 0.0    0.0
#> x21                0.0 0.0    0.0
#> x22                1.5 0.0    0.0
#> x31                0.0 1.5    0.0
#> sigma2             0.0 0.0    2.5

# center specific covariate has K = 4 categories across 5 centers; y ~ Gaussian
inv.prior.cov(X, lambda = 1, family = 'gaussian', stratified = TRUE,
              center_spec = c(3,1:4), L=5)
#>               (Intercept)_1 (Intercept)_2 (Intercept)_3 (Intercept)_4 x1 x21
#> (Intercept)_1             1             0             0             0  0   0
#> (Intercept)_2             0             1             0             0  0   0
#> (Intercept)_3             0             0             1             0  0   0
#> (Intercept)_4             0             0             0             1  0   0
#> x1                        0             0             0             0  1   0
#> x21                       0             0             0             0  0   1
#> x22                       0             0             0             0  0   0
#> x31                       0             0             0             0  0   0
#> sigma2                    0             0             0             0  0   0
#>               x22 x31 sigma2
#> (Intercept)_1   0   0      0
#> (Intercept)_2   0   0      0
#> (Intercept)_3   0   0      0
#> (Intercept)_4   0   0      0
#> x1              0   0      0
#> x21             0   0      0
#> x22             1   0      0
#> x31             0   1      0
#> sigma2          0   0      1
```
