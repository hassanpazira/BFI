# Maximum A Posteriori estimation

`MAP.estimation` function is used (in local centers) to compute Maximum
A Posterior (MAP) estimators of the parameters for Generalized Linear
Models (GLM) and Survival models.

## Usage

``` r
MAP.estimation(y,
               X,
               family = c("gaussian", "binomial", "survival"),
               Lambda,
               intercept = TRUE,
               basehaz = c("weibul", "exp", "gomp", "poly", "pwexp", "unspecified"),
               treatment = NULL,
               treat_round = NULL,
               refer_treat,
               gamma_bfi = NULL,
               RCT_propens = NULL,
               initial = NULL,
               alpha = 0.1,
               max_order = 2,
               n_intervals = 4,
               min_max_times,
               center_zero_sample = FALSE,
               zero_sample_cov,
               refer_cat,
               zero_cat,
               control = list())
```

## Arguments

- y:

  response vector. If the “`binomial`” family is used, this argument is
  a vector with entries 0 (failure) or 1 (success). Alternatively, for
  this family, the response can be a matrix where the first column is
  the number of “successes” and the second column is the number of
  “failures”. For the “`survival`” family, the response is a matrix
  where the first column is the survival time, named “time”, and the
  second column is the censoring indicator, named “status”, with 0
  indicating censoring time and 1 indicating event time.

- X:

  design matrix of dimension \\n\_{\ell} \times p\\, where \\p\\ is the
  number of covariates or predictors and \\n\_{\ell}\\ is the number of
  indeviduals in the local center. If there is a categorical covariate,
  then the function [`factor()`](https://rdrr.io/r/base/factor.html)
  should be used to encode the covariate as a factor. The same
  covariates should be used across centers when their local inference
  results are combined by
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md). The
  order of the covariates may differ across centers, provided that the
  corresponding model parameters are uniquely named.

- family:

  a description of the error distribution. This is a character string
  naming a family of the model. In the current version of the package,
  the family of model can be “`gaussian`” (with `identity` link
  function), “`binomial`” (with `logit` link function), or “`survival`”.
  Can be abbreviated. By default the “`gaussian`” family is used. In
  case of a linear regression model, `family = "gaussian"`, there is an
  extra model parameter for the variance of measurement error. While in
  the case of survival model, `family = "survival"`, the number of the
  model parameters depend on the choice of baseline hazard functions,
  see ‘Details’ for more information.

- Lambda:

  the inverse variance-covariance matrix of the Gaussian distribution
  that is used as prior distribution for the model parameters. The
  dimension of the matrix depends on the number of columns of `X`, type
  of the covariates (continuous / dichotomous or categorical), `family`,
  and whether an `intercept` is included (if applicable). However,
  `Lambda` can be easily created by
  [`inv.prior.cov()`](https://hassanpazira.github.io/BFI/reference/inv.prior.cov.md).
  See
  [`inv.prior.cov`](https://hassanpazira.github.io/BFI/reference/inv.prior.cov.md)
  for more information.

- intercept:

  logical flag for fitting an intercept. If `intercept=TRUE` (the
  default), the intercept is fitted, i.e., it is included in the model,
  and if `intercept = FALSE` it is set to zero, i.e., it's not in the
  model. This argument is not used if `family = "survival"`.

- basehaz:

  a character string representing one of the available baseline hazard
  functions; `exponential` (“`exp`”), `Weibull` (“`weibul`”, the
  default), `Gompertz` (“`gomp`”), `exponentiated polynomial`
  (“`poly`”), `piecewise constant exponential` (“`pwexp`”), and
  `unspecified baseline hazard` (“`unspecified`”). Can be abbreviated.
  It is used only when `family = "survival"`. If local sample size is
  large and the shape of the baseline hazard function is completely
  unknown, the “exponentiated polynomial” and “piecewise exponential”
  hazard functions would be preferred above the lower dimensional
  alternatives. However, if the local samples size is low, one should be
  careful using the “piecewise exponential” hazard function with many
  intervals. If `basehaz = "unspecified"`, it means that a
  (semi-parametric) Cox model is considered, and the parameters
  (regression coefficients) are estimated using the partial
  log-likelihood. If `treatment` is not `NULL`, then `basehaz` must be
  set to `"unspecified"`, as the regression coefficients are estimated
  using the weighted partial log-likelihood.

- treatment:

  a character string representing the name or place of the binary
  covariate, respectively. This covariate indicates whether the patient
  got the new treatment (\\z\_{\ell i}=1\\) or the placebo/standard
  treatment (\\z\_{\ell i}=0\\). The treatment effect is estimated when
  this argument is NOT '`NULL`'. If it is set to '`NULL`' (the default),
  the treatment effect will not be estimated. For both the first and
  second rounds, it should not be '`NULL`'. See ‘Details’.

- treat_round:

  a character string representing the `'first'` or `'second'` round of
  estimating treatment effects. In the first round,
  `treat_round = 'first'`, the local estimates of the coefficients
  (\\\gamma\_{\ell}\\) is estimated. In the second round,
  `treat_round = 'second'`, the treatment effect, `propensity` scores
  and the statistical summaries (`for_ATE`, only for 'binomial' and
  'gaussian' families) are calculated to be sent to the central server
  for estimating the BFI treatment effect (\\\hat \zeta\_{BFI}\\) and
  average treatment effects (ATEs).

- refer_treat:

  a character string representing the reference category of the
  treatment variable. The reference category is considered as \\z\_{\ell
  i}=0\\. This argument is used when `treatment` is not '`NULL`'.
  Default is `refer_treat = levels(X$treatment)[1]`.

- gamma_bfi:

  a vector specifying the BFI estimates of the coefficients received
  from the central server in the first round. It can be defined by the
  output of `MAP.estimation()$theta_hat` obtained from the first round.
  The length of `gamma_bfi` equals the number of regression
  coefficients, including the intercept if `intercept=TRUE`, but
  excluding \\\zeta\\, which represents the treatment effect, as well as
  the nuisance parameter \\\sigma\\ in the `gaussian` family and any
  parameters of the baseline hazard (\\\boldsymbol{\omega}\\) for
  `survival`. This argument is used only when the argument `treatment`
  is not '`NULL`'. If `treatment` is not '`NULL`' but
  `gamma_bfi = NULL`, then the argument `RCT_propens` must not be
  '`NULL`', indicating an RCT study. See ‘Details’.

- RCT_propens:

  a vector specifying the propensity scores, which represent the
  probability of receiving the treatment given the covariates, which are
  known in the randomized studies (RCTs). For example, in a 1:1
  randomized trial, the propensity scores are, by definition, equal to
  1/2 (or 0.5), whereas in an unbalanced randomized trial, e.g., a 2:1
  trial, the propensity scores are now 2/3 and 1/3 for the two arms,
  respectively. The length of `RCT_propens` equals to the number of
  individuals in the local center denoted as \\n\_{\ell}\\. This
  argument is used only when the study is a randomized control trial,
  i.e., the propensity scores are known for this local center. In this
  case, there is only ‘one’ round, and the argument `treatment` must not
  be '`NULL`', whereas `gamma_bfi = NULL`. Indeed, when 'treatment' is
  not '`NULL`', one of the arguments '`RCT_propens`' or '`gamma_bfi`'
  could be '`NULL`'. See ‘Details’.

- initial:

  a vector specifying initial values for the parameters to be optimized
  over. The length of `initial` is equal to the number of model
  parameters and thus, is equal to the number of rows or columns of
  `Lambda`. Since the `'L-BFGS-B'` method is used in the algorithm,
  these values should always be finite. Default is a vector of zeros,
  except for the `survival` family with the `poly` function, where it is
  a vector with the first \\p\\ elements as zeros for coefficients
  (\\\boldsymbol{\beta}\\) and -0.5 for the remaining parameters
  (\\\boldsymbol{\omega}\\). For the `gaussian` family, the last element
  of the `initial` vector could also be considered negative, because the
  optimization for this parameter is carried out on the
  \\\log(\sigma^2)\\ scale.

- alpha:

  a significance level used in the chi-squared distribution (with one
  degree of freedom and 1-\\\alpha\\ representing the upper quantile) to
  conduct a likelihood ratio test for obtaining the order of the
  exponentiated polynomial baseline hazard function. It is only used
  when `family = "survival"` and `basehaz = "poly"`. Default is 0.1. See
  ‘Details’.

- max_order:

  an integer representing the maximum value of `q_l`, which is the
  order/degree minus 1 of the exponentiated polynomial baseline hazard
  function. This argument is only used when `family = "survival"` and
  `basehaz = "poly"`. Default is 2.

- n_intervals:

  an integer representing the number of intervals in the piecewise
  exponential baseline hazard function. This argument is only used when
  `family = "survival"` and `basehaz = "pwexp"`. Default is 4.

- min_max_times:

  a scalar representing the minimum of the maximum event times observed
  in the centers. The value of this argument should be defined by the
  central server (which has access to the maximum event times of all the
  centers) and is only used when `family = "survival"` and
  `basehaz = "pwexp"`.

- center_zero_sample:

  logical flag indicating whether the center has a categorical covariate
  with no observations/individuals in one of the categories. Default is
  `center_zero_sample = FALSE`.

- zero_sample_cov:

  either a character string or an integer representing the categorical
  covariate that has no samples/observations in one of its categories.
  This covariate should have at least two categories, one of which is
  the reference. It is used when `center_zero_sample = TRUE`.

- refer_cat:

  a character string representing the reference category. The category
  with no observations (the argument `zero_cat`) cannot be used as the
  reference in the argument `refer_cat`. It is used when
  `center_zero_sample = TRUE`.

- zero_cat:

  a character string representing the category with no
  samples/observations. It is used when `center_zero_sample = TRUE`.

- control:

  a list of control parameters. See ‘Details’.

## Value

`MAP.estimation` returns a list containing the following components:

- theta_hat:

  the vector corresponding to the maximum a posteriori (MAP) estimates
  of the parameters. For the `gaussian` family, the dispersion parameter
  is estimated on the \\\log(\sigma^2)\\ scale and the last element of
  this vector is back-transformed to \\\sigma^2\\. When `treatment` is
  not NULL and `treat_round = 'first'`, this is the MAP estimates of
  only regression coefficients (\\\boldsymbol{\gamma}\_\ell\\) except
  the treatment effect \\\zeta\_\ell\\. While `treatment` is not NULL
  and `treat_round = 'second'`, this is \\\hat{\zeta}\_\ell\\, the
  weighted MAP estimate of the treatment effect \\ \zeta\_\ell\\ in
  center \\\ell\\;

- A_hat:

  minus the curvature (or Hessian) matrix of the log-posterior around
  the point `theta_hat`. The dimension of the matrix is the same as the
  argument `Lambda`. For the `gaussian` family the row and column
  labelled `sigma2` refers to the \\\log(\sigma^2)\\ scale, on which the
  dispersion parameter is estimated;

- sd:

  the vector of (posterior) standard deviation of the MAP estimates in
  `theta_hat`, that is `sqrt(diag(solve(A_hat)))`. For the `gaussian`
  family the element labelled `sigma2` is the standard deviation on the
  \\\log(\sigma^2)\\ scale;

- Lambda:

  the inverse variance-covariance matrix of the Gaussian distribution
  that is used as prior distribution for the parameters. It's exactly
  the same as the argument `Lambda`;

- formula:

  the formula applied;

- names:

  the names of the model parameters;

- n:

  sample size, \\n\_{\ell}\\;

- np:

  the number of coefficients;

- q_l:

  the order/degree minus 1 of the exponentiated polynomial baseline
  hazard function determined for the current center by the likelihood
  ratio test. This output argument, `q_l`, is only shown when
  `family = "survival"` and `basehaz = "poly"`, and will be used in the
  function
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md);

- theta_A_poly:

  an array where the first component is a matrix with columns
  representing the MAP estimates of the parameters for different
  `q_l`'s, i.e., `q_l`, `q_l`+1, ..., `max_order`. The other components
  are minus the curvature matrices for different `q_l`'s, i.e., `q_l`,
  `q_l`+1, ..., `max_order`. Therefore, the first non-NA curvature
  matrix is equal to the output argument `A_hat`. This output argument,
  `theta_A_poly`, is only shown if `family = "survival"` and
  `basehaz = "poly"`, and will be used in the function
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md);

- lev_no_ref_zer:

  a vector containing the names of the levels of the categorical
  covariate that has no samples/observations in one of its categories.
  The name of the category with no samples and the name of the reference
  category are excluded from this vector. This argument is shown when
  `family = "survival"` and `basehaz = "poly"`, and will be used in the
  function
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md);

- treatment:

  a character string representing the name or place of the binary
  covariate, respectively. If it is set to '`NULL`', the treatment
  effect will not be estimated;

- refer_treat:

  the reference category of the treatment. It is shown when `treatment`
  is not '`NULL`', and can be used in the function
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md);

- gamma_bfi:

  a vector specifying the BFI estimates of the coefficients received
  from the central server in the first round. If `treatment = NULL`,
  then `gamma_bfi` must also be '`NULL`';

- RCT_propens:

  a vector specifying the propensity scores, which represent the
  probability of receiving the treatment given the covariates, which are
  known in the randomized studies (RCTs). If `treatment = NULL`, then
  `RCT_propens` must also be '`NULL`';

- propensity:

  a vector specifying the propensity scores (the probability a patient
  gets the treatment given the characteristics measured at baseline)
  calculated by \\Pr(Z\_\ell = 1 \| X\_\ell)\\;

- for_ATE:

  a vector used in the central server to calculate the average treatment
  effect (ATE). For `family` of `binomial` and `gaussian`, its elements
  are:

  first

  :   \\m\_{\ell 1}\\, the number of patients in the treatment group,
      where \\n\_{\ell} = m\_{\ell 1} + m\_{\ell 2}\\,

  second

  :   \\m\_{\ell 2}\\, the number of patients in the reference group,
      where \\m\_{\ell 2} = n\_{\ell} - m\_{\ell 1}\\,

  third

  :   \\\sum\_{i=1}^{n\_{\ell}} z\_{\ell i} y\_{\ell i}\\,

  fourth

  :   \\\sum\_{i=1}^{n\_{\ell}} (z\_{\ell i} y\_{\ell i})^{2} \\,

  fifth

  :   \\\sum\_{i=1}^{n\_{\ell}} z\_{\ell i} / e\_{\ell i}\\,

  sixth

  :   \\\sum\_{i=1}^{n\_{\ell}} z\_{\ell i} y\_{\ell i} / e\_{\ell i}\\,

  seventh

  :   \\\sum\_{i=1}^{n\_{\ell}} (1 - z\_{\ell i}) / (1 - e\_{\ell i})\\,

  eighth

  :   \\\sum\_{i=1}^{n\_{\ell}} (1 - z\_{\ell i}) y\_{\ell i} / (1 -
      e\_{\ell i})\\,

  ninth

  :   \\\sum\_{i=1}^{n\_{\ell}} (1 - z\_{\ell i}) y\_{\ell i}\\,

  but for `survival`, it's 'NULL';

- zero_sample_cov:

  the categorical covariate that has no samples/observations in one of
  its categories. It is shown when `center_zero_sample = TRUE`, and can
  be used in the function
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md);

- refer_cat:

  the reference category. It is shown when `center_zero_sample = TRUE`,
  and can be used in the function
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md);

- zero_cat:

  the category with no samples/observations. It is shown when
  `center_zero_sample = TRUE`, and can be used in the function
  [`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md);

- value:

  the value of minus the log-likelihood posterior density evaluated at
  `theta_hat`;

- family:

  the `family` used;

- basehaz:

  the baseline hazard function used;

- intercept:

  logical flag used to fit an intercept if `TRUE`, or set to zero if
  `FALSE`;

- convergence:

  an integer value used to encode the warnings and the errors related to
  the algorithm used to fit the model. The values returned are:

  0

  :   algorithm has converged,

  1

  :   maximum number of iterations ('`maxit`') has been reached,

  2

  :   Warning from the 'L-BFGS-B' method. See the message after this
      value;

- control:

  the list of control parameters used to compute the MAP estimates.

## Details

`MAP.estimation` function finds the Maximum A Posteriori (MAP) estimates
of the model parameters by maximizing the log-posterior density with
respect to the parameters, i.e., the estimates equal the values for
which the log-posterior density is maximal (the posterior mode). In
other words, `MAP.estimation()` optimizes the log-posterior density with
respect to the parameter vector to obtain its MAP estimates. In addition
to the model parameters (i.e., coefficients (\\\boldsymbol{\beta}\\) and
variance error (\\\sigma^2_e\\) for `gaussian` or the parameters of the
baseline hazard (\\\boldsymbol{\omega}\\) for `survival`), the curvature
matrix (Hessian of the log-posterior) is estimated around the mode.

For the `gaussian` family the dispersion parameter is estimated on the
\\\log(\sigma^2)\\ scale, and the corresponding row and column of the
curvature matrix refer to that scale.

The `MAP.estimation` function returns an object of class \``bfi`\`.
Therefore,
[`summary()`](https://hassanpazira.github.io/BFI/reference/summary.bfi.md)
can be used for the object returned by `MAP.estimation()`.

For the case where `family = "survival"` and `basehaz = "poly"`, we
assume that in all centers the \\q\_\ell\\'s are equal. However, the
order of the estimated polynomials may vary across the centers so that
each center can have different number of parameters, say \\q\_\ell\\+1.
After obtaining the estimates within the local centers (by using
`MAP.estimation()`) and having all estimates in the central server, we
choose the order of the polynomial approximation for the combined data
to be the maximum of the orders of the local polynomial functions, i.e.,
\\\max \\q_1, \ldots, q_L \\\\, to approximate the global baseline
hazard (exponentiated polynomial) function more accurately. This is
because the higher-order polynomial approximation can capture more
complex features and details in the combined data. Using the
higher-order approximation ensures that we account for the higher-order
moments and features present in the data while maintaining accuracy. As
a result, all potential cases are stored in the `theta_A_poly` argument
to be used in
[`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md) by the
central server. For further information on the `survival` family, refer
to the 'References' section.

The three arguments `'treatment'`, `'treat_round'`, `'refer_treat'`,
`'gamma_bfi'`, and `'RCT_propens'` are related to the estimation of the
treatment effect. For observational and non-randomized studies, the
treatment effect is estimated in two rounds; In the first round,
\\\hat{\boldsymbol{\beta}}\_{\ell}\\ (or
\\\hat{\boldsymbol{\gamma}}\_{\ell}\\) are estimated locally and in the
central server \\\hat{\boldsymbol{\beta}}\_{BFI}\\ (or
\\\hat{\boldsymbol{\gamma}}\_{BFI}\\) is estimated and then is sent to
all local centers for the second round to estimate propensity scores,
weights, treatment effect and ATEs. In the first round, the argument
`treatment` should not be '`NULL`' and `treat_round = "first"`, while
`gamma_bfi = NULL` and `RCT_propens = NULL`. Moreover, in the first
round, the family must be set to `binomial`, however this is handled
automatically. In the second round, local weighted MAP estimate of the
treatment effects and propensity scores are estimated, and along with
some summary statistics are sent to the central server to estimate the
average treatment effects ATEs (in this case `treatment` and `gamma_bfi`
should not be '`NULL`' and `treat_round = "second"`, but
`RCT_propens = NULL`). In contrast, for the randomized control trial
(RCT), the treatment effect can be estimated by only one round as the
propensity scores are known (in this case `treatment` and `RCT_propens`
should not be '`NULL`', but `gamma_bfi = NULL`). NOTE: the argument
`gamma_bfi` should not include estimates of the nuisance parameter
\\\sigma\\ in the `gaussian` family or any parameters of the baseline
hazard (\\\boldsymbol{\omega}\\) and the intercept for `survival`. For
more examples on treatment effect estimation, see the ‘Examples’ section
of [`bfi`](https://hassanpazira.github.io/BFI/reference/BFI.md).

To solve unconstrained and bound-constrained optimization problems, the
`MAP.estimation` function utilizes an optimization algorithm called
Limited-memory Broyden-Fletcher-Goldfarb-Shanno with Bound Constraints
(L-BFGS-B), Byrd et. al. (1995). The L-BFGS-B algorithm is a
limited-memory “quasi-Newton” method that iteratively updates the
parameter estimates by approximating the inverse Hessian matrix using
gradient information from the history of previous iterations. This
approach allows the algorithm to approximate the curvature of the
posterior distribution and efficiently search for the optimal solution,
which makes it computationally efficient for problems with a large
number of variables.

By default, the algorithm uses a relative change in the objective
function as the convergence criterion. When the change in the objective
function between iterations falls below a certain threshold
(\``factr`\`) the algorithm is considered to have converged. The
convergence can be checked with the argument `convergence` in the
output. See ‘Value’.

In case of convergence issue, it may be necessary to investigate and
adjust optimization parameters to facilitate convergence. It can be done
using the `initial` and `control` arguments. By the argument `initial`
the initial points of the interative optimization algorithm can be
changed, and the argument `control` is a list that can supply any of the
following components:

- `maxit`::

  is the maximum number of iterations. Default is 150;

- `factr`::

  controls the convergence of the `'L-BFGS-B'` method. Convergence
  occurs when the reduction in the objective is within this factor of
  the machine tolerance. Default for `factr` is 1e7, which gives a
  tolerance of about 1e-9. The exact tolerance can be checked by
  multiplying this value by `.Machine$double.eps`;

- `pgtol`::

  helps to control the convergence of the `'L-BFGS-B'` method. It is a
  tolerance on the projected gradient in the current search direction,
  i.e., the iteration will stop when the maximum component of the
  projected gradient is less than or equal to `pgtol`, where
  `pgtol`\\\geq 0\\. Default is zero, when the check is suppressed;

- `trace`::

  is a non-negative integer. If positive, tracing information on the
  progress of the optimization is produced. Higher values may produce
  more tracing information: for the method `'L-BFGS-B'` there are six
  levels of tracing. To understand exactly what these do see the source
  code of `optim` function in the
  [stats](https://rdrr.io/r/stats/stats-package.html) package;

- `REPORT`::

  is the frequency of reports for the `'L-BFGS-B'` method if
  `'control$trace'` is positive. Default is every 10 iterations;

- `lmm`::

  is an integer giving the number of `BFGS` updates retained in the
  `'L-BFGS-B'` method. Default is 5.

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

Byrd R.H., Lu P., Nocedal J. and Zhu C. (1995). *A limited memory
algorithm for bound constrained optimization*. SIAM Journal on
Scientific Computing, 16, 1190-1208. \<https://doi.org/10.1137/0916069\>

## Author

Hassan Pazira and Marianne Jonker  
Maintainer: Hassan Pazira \<h.pazira@arq.org\>

## See also

[`bfi`](https://hassanpazira.github.io/BFI/reference/BFI.md),
[`inv.prior.cov`](https://hassanpazira.github.io/BFI/reference/inv.prior.cov.md)
and
[`summary.bfi`](https://hassanpazira.github.io/BFI/reference/summary.bfi.md)

## Examples

``` r
###--------------###
### y ~ Gaussian ###
###--------------###

# Setting a seed for reproducibility
set.seed(11235)

# model parameters: coefficients and sigma2 = 1.5
theta <- c(1, 2, 2, 2, 1.5)

#----------------
# Data Simulation
#----------------
n   <- 30   # sample size
p   <- 3    # number of coefficients without intercept
X   <- data.frame(matrix(rnorm(n * p), n, p)) # continuous variables
# linear predictor:
eta <- theta[1] + theta[2] * X$X1 + theta[3] * X$X2 + theta[4] * X$X3
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu  <- gaussian()$linkinv(eta)
y   <- rnorm(n, mu, sd = sqrt(theta[5]))

# Load the BFI package
library(BFI)

#-----------------------------------------------
# MAP estimations for theta and curvature matrix
#-----------------------------------------------
# MAP estimates with 'intercept'
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = "gaussian")
(fit <- MAP.estimation(y, X, family = "gaussian", Lambda))
#> $theta_hat
#> (Intercept)          X1          X2          X3      sigma2 
#>    1.341258    2.236391    2.071001    2.164002    1.054571 
#> 
#> $A_hat
#>             (Intercept)         X1         X2         X3     sigma2
#> (Intercept)  28.5475747  1.9315730 -7.7385763  6.2804549  0.1341231
#> X1            1.9315730 20.4369223 -2.2127254  4.5673270  0.2236029
#> X2           -7.7385763 -2.2127254 49.3206700  3.1359341  0.2070737
#> X3            6.2804549  4.5673270  3.1359341 24.3657850  0.2163723
#> sigma2        0.1341231  0.2236029  0.2070737  0.2163723 16.0545936
#> 
#> $sd
#> (Intercept)          X1          X2          X3      sigma2 
#>   0.1982976   0.2269483   0.1476575   0.2153269   0.2496131 
#> 
#> $Lambda
#>             (Intercept)  X1  X2  X3 sigma2
#> (Intercept)         0.1 0.0 0.0 0.0      0
#> X1                  0.0 0.1 0.0 0.0      0
#> X2                  0.0 0.0 0.1 0.0      0
#> X3                  0.0 0.0 0.0 0.1      0
#> sigma2              0.0 0.0 0.0 0.0      1
#> 
#> $formula
#> [1] y ~ X1 + X2 + X3
#> 
#> $names
#> [1] "(Intercept)" "X1"          "X2"          "X3"          "sigma2"     
#> 
#> $n
#> [1] 30
#> 
#> $np
#> [1] 4
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
#> [1] 35.28046
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
class(fit)
#> [1] "bfi"
summary(fit, cur_mat = TRUE)
#> 
#> Summary of the local model:
#> 
#>    Formula: y ~ X1 + X2 + X3 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   1.3413  0.1983  0.9526   1.7299
#> X1            2.2364  0.2269  1.7916   2.6812
#> X2            2.0710  0.1477  1.7816   2.3604
#> X3            2.1640  0.2153  1.7420   2.5860
#> 
#> Dispersion parameter (sigma2):  1.055 
#>             log Lik Posterior:  -35.28 
#>                   Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>             (Intercept)      X1      X2      X3  sigma2
#> (Intercept)     28.5476  1.9316 -7.7386  6.2805  0.1341
#> X1               1.9316 20.4369 -2.2127  4.5673  0.2236
#> X2              -7.7386 -2.2127 49.3207  3.1359  0.2071
#> X3               6.2805  4.5673  3.1359 24.3658  0.2164
#> sigma2           0.1341  0.2236  0.2071  0.2164 16.0546

# MAP estimates without 'intercept'
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = 'gaussian',
                        intercept = FALSE)
(fit1 <- MAP.estimation(y, X, family = 'gaussian', Lambda, intercept = FALSE))
#> $theta_hat
#>       X1       X2       X3   sigma2 
#> 2.241637 1.832740 2.525235 2.496099 
#> 
#> $A_hat
#>                X1         X2         X3     sigma2
#> X1      8.6921038 -0.9348497  1.9296405  0.2241621
#> X2     -0.9348497 20.8951379  1.3248942  0.1832670
#> X3      1.9296405  1.3248942 10.3520007  0.2525184
#> sigma2  0.2241621  0.1832670  0.2525184 17.4960821
#> 
#> $sd
#>        X1        X2        X3    sigma2 
#> 0.3478802 0.2205610 0.3192969 0.2391506 
#> 
#> $Lambda
#>         X1  X2  X3 sigma2
#> X1     0.1 0.0 0.0      0
#> X2     0.0 0.1 0.0      0
#> X3     0.0 0.0 0.1      0
#> sigma2 0.0 0.0 0.0      1
#> 
#> $formula
#> [1] y ~ X1 + X2 + X3
#> 
#> $names
#> [1] "X1"     "X2"     "X3"     "sigma2"
#> 
#> $n
#> [1] 30
#> 
#> $np
#> [1] 3
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
#> [1] 63.9101
#> 
#> $family
#> [1] "gaussian"
#> 
#> $basehaz
#> [1] "weibul"      "exp"         "gomp"        "poly"        "pwexp"      
#> [6] "unspecified"
#> 
#> $intercept
#> [1] FALSE
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
summary(fit1, cur_mat = TRUE)
#> 
#> Summary of the local model:
#> 
#>    Formula: y ~ X1 + X2 + X3 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   2.2416  0.3479  1.5598   2.9235
#> X2   1.8327  0.2206  1.4004   2.2650
#> X3   2.5252  0.3193  1.8994   3.1510
#> 
#> Dispersion parameter (sigma2):  2.496 
#>             log Lik Posterior:  -63.91 
#>                   Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>             X1      X2      X3  sigma2
#> X1      8.6921 -0.9348  1.9296  0.2242
#> X2     -0.9348 20.8951  1.3249  0.1833
#> X3      1.9296  1.3249 10.3520  0.2525
#> sigma2  0.2242  0.1833  0.2525 17.4961



###-----------------###
### Survival family ###
###-----------------###

# Setting a seed for reproducibility
set.seed(112358)

#-------------------------
# Simulating Survival data
#-------------------------
n    <- 50
beta <- 1:4
p    <- length(beta)
X    <- data.frame(matrix(rnorm(n * p), n, p)) # continuous (normal) variables

## Simulating survival data from Weibull with a predefined censoring rate of 0.3
y <- surv.simulate(Z = list(X), beta = beta, a = 5, b = exp(1.8), u1 = 0.1,
                   cen_rate = 0.3, gen_data_from = "weibul")$D[[1]][, 1:2]

#---------------------------------------
# MAP estimations with "weibul" function
#---------------------------------------
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = 'survival',
                        basehaz = "weibul")
fit2 <- MAP.estimation(y, X, family = 'survival', Lambda = Lambda,
                       basehaz = "weibul")
fit2
#> $theta_hat
#>        X1        X2        X3        X4   omega_1   omega_2 
#> 0.8931059 2.1928543 3.2725773 4.3943088 1.3619890 1.8453340 
#> 
#> $A_hat
#>                  X1         X2          X3          X4     omega_1    omega_2
#> X1       23.9284583   1.055064   0.2149551    4.366821   11.724899  -48.66091
#> X2        1.0550636  24.509514  -5.7248684   -3.048296   -2.381583  -14.77898
#> X3        0.2149551  -5.724868  32.7991448    1.695595    3.371813  -93.43592
#> X4        4.3668210  -3.048296   1.6955949   33.883288   15.514630 -163.46079
#> omega_1  11.7248989  -2.381583   3.3718129   15.514630   36.638199 -116.92019
#> omega_2 -48.6609087 -14.778978 -93.4359197 -163.460791 -116.920190 1267.08368
#> 
#> $sd
#>         X1         X2         X3         X4    omega_1    omega_2 
#> 0.23091779 0.26922733 0.31532254 0.43683670 0.22009211 0.09129591 
#> 
#> $Lambda
#>          X1  X2  X3  X4 omega_1 omega_2
#> X1      0.1 0.0 0.0 0.0       0       0
#> X2      0.0 0.1 0.0 0.0       0       0
#> X3      0.0 0.0 0.1 0.0       0       0
#> X4      0.0 0.0 0.0 0.1       0       0
#> omega_1 0.0 0.0 0.0 0.0       1       0
#> omega_2 0.0 0.0 0.0 0.0       0       1
#> 
#> $formula
#> [1] "Survival(time, status) ~ X1 + X2 + X3 + X4"
#> 
#> $names
#> [1] "X1"      "X2"      "X3"      "X4"      "omega_1" "omega_2"
#> 
#> $n
#> [1] 50
#> 
#> $np
#> [1] 4
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
#> [1] -29.03443
#> 
#> $family
#> [1] "survival"
#> 
#> $basehaz
#> [1] "weibul"
#> 
#> $intercept
#> [1] FALSE
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
summary(fit2, cur_mat = TRUE)
#> 
#> Summary of the local model:
#> 
#>    Formula: Survival(time, status) ~ X1 + X2 + X3 + X4 
#>     Family: ‘survival’ 
#>   Baseline: ‘weibul’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   0.8931  0.2309  0.4405   1.3457
#> X2   2.1929  0.2692  1.6652   2.7205
#> X3   3.2726  0.3153  2.6546   3.8906
#> X4   4.3943  0.4368  3.5381   5.2505
#> 
#> log Lik Posterior:  29.03 
#>       Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>               X1       X2       X3        X4   omega_1   omega_2
#> X1       23.9285   1.0551   0.2150    4.3668   11.7249  -48.6609
#> X2        1.0551  24.5095  -5.7249   -3.0483   -2.3816  -14.7790
#> X3        0.2150  -5.7249  32.7991    1.6956    3.3718  -93.4359
#> X4        4.3668  -3.0483   1.6956   33.8833   15.5146 -163.4608
#> omega_1  11.7249  -2.3816   3.3718   15.5146   36.6382 -116.9202
#> omega_2 -48.6609 -14.7790 -93.4359 -163.4608 -116.9202 1267.0837

#-------------------------------------
# MAP estimations with "poly" function
#-------------------------------------
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = 'survival',
                        basehaz = 'poly')
fit3 <- MAP.estimation(y, X, family = "survival", Lambda = Lambda,
                       basehaz = "poly")
# Degree of the exponentiated polynomial baseline hazard
fit3$q_l + 1
#> [1] 3
# theta_hat for (beta_1, ..., beta_p, omega_0, ..., omega_{q_l})
fit3$theta_A_poly[,,1][,fit3$q_l+1] # equal to fit3$theta_hat
#>         X1         X2         X3         X4    omega_0    omega_1    omega_2 
#>  0.5099345  0.6924321  1.1425359  1.4261694 -1.5219761  2.0961028  0.5668222 
# A_hat
fit3$theta_A_poly[,,fit3$q_l+2] # equal to fit3$A_hat
#>                 [,1]        [,2]         [,3]       [,4]      [,5]      [,6]
#> X1      28.946657217  0.05866364  0.001076572   8.480412 11.763104  4.384700
#> X2       0.058663644 32.59787625 -8.176657379   5.843179 -2.231511 -6.166270
#> X3       0.001076572 -8.17665738 36.893828905   3.283148  3.584772 -4.062458
#> X4       8.480411839  5.84317926  3.283147626  58.816723 15.811387 -6.021616
#> omega_0 11.763104225 -2.23151073  3.584772096  15.811387 48.521908 17.034602
#> omega_1  4.384699587 -6.16626999 -4.062457601  -6.021616 17.034602 26.901089
#> omega_2  2.918077452 -7.87506193 -5.153277223 -11.037167 17.901089 22.460304
#>               [,7]
#> X1        2.918077
#> X2       -7.875062
#> X3       -5.153277
#> X4      -11.037167
#> omega_0  17.901089
#> omega_1  22.460304
#> omega_2  31.360057
summary(fit3, cur_mat = TRUE)
#> 
#> Summary of the local model:
#> 
#>    Formula: Survival(time, status) ~ X1 + X2 + X3 + X4 
#>     Family: ‘survival’ 
#>   Baseline: ‘poly’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   0.5099  0.1979  0.1220   0.8979
#> X2   0.6924  0.1905  0.3191   1.0657
#> X3   1.1425  0.1783  0.7931   1.4919
#> X4   1.4262  0.1554  1.1216   1.7307
#> 
#> log Lik Posterior:  1.174 
#>       Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>              X1      X2      X3       X4 omega_0 omega_1  omega_2
#> X1      28.9467  0.0587  0.0011   8.4804 11.7631  4.3847   2.9181
#> X2       0.0587 32.5979 -8.1767   5.8432 -2.2315 -6.1663  -7.8751
#> X3       0.0011 -8.1767 36.8938   3.2831  3.5848 -4.0625  -5.1533
#> X4       8.4804  5.8432  3.2831  58.8167 15.8114 -6.0216 -11.0372
#> omega_0 11.7631 -2.2315  3.5848  15.8114 48.5219 17.0346  17.9011
#> omega_1  4.3847 -6.1663 -4.0625  -6.0216 17.0346 26.9011  22.4603
#> omega_2  2.9181 -7.8751 -5.1533 -11.0372 17.9011 22.4603  31.3601

#------------------------------------------------------
# MAP estimations with "pwexp" function with 3 intervals
#-------------------------------------------------------
# Assume we have 4 centers
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = 'survival',
                        basehaz = 'pwexp', n_intervals = 3)
# For this baseline hazard ("pwexp"), we need to know
# maximum survival times of the 3 other centers:
max_times <- c(max(rexp(30)), max(rexp(50)), max(rexp(70)))
# Minimum of the maximum values of the survival times of all 4 centers is:
min_max_times <- min(max(y$time), max_times)
fit4 <- MAP.estimation(y, X, family = "survival", Lambda = Lambda,
                       basehaz = "pwexp", n_intervals = 3,
                       min_max_times=max(y$time))
#> 
#>  No. observations in the intervals :  26 16 8 
#>  
summary(fit4, cur_mat = TRUE)
#> 
#> Summary of the local model:
#> 
#>    Formula: Survival(time, status) ~ X1 + X2 + X3 + X4 
#>     Family: ‘survival’ 
#>   Baseline: ‘pwexp’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   0.5313  0.2316  0.0774   0.9851
#> X2   0.5870  0.2235  0.1490   1.0249
#> X3   1.0772  0.2207  0.6446   1.5098
#> X4   1.2916  0.2002  0.8992   1.6841
#> 
#> log Lik Posterior:  -2.408 
#>       Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>              X1      X2      X3      X4 omega_1 omega_2 omega_3
#> X1      23.6816  0.4579 -1.0144  7.0250  7.6918  3.4509  0.6183
#> X2       0.4579 25.5901 -8.0729  4.6451  1.7484 -1.7688 -2.2004
#> X3      -1.0144 -8.0729 31.9807  4.0945  8.3609 -4.1505 -0.6195
#> X4       7.0250  4.6451  4.0945 48.2644 19.3602 -1.4761 -2.0592
#> omega_1  7.6918  1.7484  8.3609 19.3602 18.1279  0.0000  0.0000
#> omega_2  3.4509 -1.7688 -4.1505 -1.4761  0.0000  9.1151  0.0000
#> omega_3  0.6183 -2.2004 -0.6195 -2.0592  0.0000  0.0000  4.1179


#--------------------------
# Semi-parametric Cox model
#--------------------------
Lambda <- inv.prior.cov(X, lambda = c(0.1), family = 'survival',
                        basehaz = "unspecified")
fit5 <- MAP.estimation(y, X, family = 'survival', Lambda = Lambda,
                       basehaz = "unspecified")
summary(fit5, cur_mat = TRUE)
#> 
#> Summary of the local model:
#> 
#>    Formula: Survival(time, status) ~ X1 + X2 + X3 + X4 
#>     Family: ‘survival’ 
#>   Baseline: ‘unspecified’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   1.0138  0.3045  0.4169   1.6107
#> X2   2.4014  0.4784  1.4638   3.3390
#> X3   3.6767  0.6560  2.3910   4.9625
#> X4   4.5452  0.7908  2.9952   6.0951
#> 
#> log Lik Posterior:  -36.18 
#>       Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>         X1      X2      X3      X4
#> X1 14.1200 -0.1991 -3.8524  0.7846
#> X2 -0.1991 15.3394 -7.7693 -1.5301
#> X3 -3.8524 -7.7693 14.4851 -6.2406
#> X4  0.7846 -1.5301 -6.2406  6.6919

```
