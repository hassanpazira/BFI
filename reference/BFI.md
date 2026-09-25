# Bayesian Federated Inference

`bfi` function can be used (on the central server) to combine inference
results from separate datasets (without combining the data) to
approximate what would have been inferred had the datasets been merged.
This function can handle linear, logistic and survival regression
models.

## Usage

``` r
bfi(theta_hats = NULL,
    A_hats,
    Lambda,
    family = c("gaussian", "binomial", "survival"),
    basehaz = c("weibul", "exp", "gomp", "poly", "pwexp", "unspecified"),
    stratified = FALSE,
    strat_par = NULL,
    center_spec = NULL,
    theta_A_polys = NULL,
    treat_round = NULL,
    for_ATE = NULL,
    p,
    q_ls,
    center_zero_sample = FALSE,
    which_cent_zeros,
    zero_sample_covs,
    refer_cats,
    zero_cats,
    lev_no_ref_zeros)
```

## Arguments

- theta_hats:

  a list of \\L\\ vectors of the maximum a posteriori (MAP) estimates of
  the model parameters in the \\L\\ centers. These vectors must have
  equal dimensions. See ‘Details’.

- A_hats:

  a list of \\L\\ minus curvature matrices for \\L\\ centers. These
  matrices must have equal dimensions. See ‘Details’.

- Lambda:

  a list of \\L+1\\ matrices. The \\k^{\th}\\ matrix is the chosen
  inverse variance-covariance matrix of the Gaussian distribution that
  is used as prior distribution in center \\k\\, where
  \\k=1,2,\ldots,L\\. The last matrix is the chosen variance-covariance
  matrix for the Gaussian prior of the (fictive) combined data set. If
  `stratified = FALSE`, all \\L+1\\ matrices must have equal dimensions.
  While, if `stratified = TRUE`, the first \\L\\ matrices must have
  equal dimensions and the last matrix should have a different (greater)
  dimension than the others. See ‘Details’.

- family:

  a character string representing the family name used for the local
  centers. Can be abbreviated.

- basehaz:

  a character string representing one of the available baseline hazard
  functions; `exponential` (`"exp"`), `Weibull` (`"weibul"`, the
  default), `Gompertz` (`"gomp"`), `exponentiated polynomial`
  (`"poly"`), `piecewise exponential` (`"pwexp"`), and
  `unspecified baseline hazard` (`"unspecified"`). It is only used when
  `family = "survival"`. Can be abbreviated. If
  `basehaz = "unspecified"`, it means that a (semi-parametric) Cox model
  is considered, and the parameters (regression coefficients) are
  estimated using the partial log-likelihood.

- stratified:

  logical flag for performing the stratified analysis. If
  `stratified = TRUE`, the parameter(s) selected in the `strat_par`
  argument are allowed to be different across centers (to deal with
  heterogeneity across centers), except when the argument `center_spec`
  is not `NULL`. Default is `stratified = FALSE`. See ‘Details’ and
  ‘Examples’.

- strat_par:

  an integer vector for indicating the stratification parameter(s). It
  can be used to deal with heterogeneity due to center-specific
  parameters. For the `"binomial"` and `"gaussian"` families it is a
  one- or two-element integer vector so that the values \\1\\ and/or
  \\2\\ are/is used to indicate that the “intercept” and/or “sigma2” are
  allowed to vary, respectively. For the `"binomial"` family the length
  of the vector should be at most one which refers to “intercept”, and
  the value of this element should be \\1\\ (to handle heterogeneity
  across outcome means). For `"gaussian"` this vector can be \\1\\ for
  indicating the “intercept” only (handling heterogeneity across outcome
  means), \\2\\ for indicating the “sigma2” only (handling heterogeneity
  due to nuisance parameter), and \\c(1, 2)\\ for both “intercept” and
  “sigma2”. When `family = "survival"`, this vector can contain any
  combination of values ranging from 1 to the maximum number of
  parameters of the baseline hazard function, i.e., \\1\\ for `"exp"`,
  \\2\\ for `"weibul"` and `"gomp"`, `max_order + 1` for `"poly"`, and
  `n_intervals` for `"pwexp"`. For example, for `"weibul"`, `strat_par`
  could be \\1\\, \\2\\ or \\c(1, 2)\\, where \\1\\ represents
  \\\omega_1\\ and \\2\\ represents \\\omega_2\\. This argument is used
  only when `stratified = TRUE` and `center_spec = NULL`. Default is
  `strat_par = NULL`. See ‘Details’ and ‘Examples’.

- center_spec:

  a vector of \\L\\ elements to account for the heterogeneity across
  centers due to clustering. This argument is used only when
  `stratified = TRUE` and `strat_par = NULL`. Each element represents a
  specific feature of the corresponding center. There must be only one
  specific value or attribute for each center. This vector could be a
  numeric, characteristic or factor vector. Note that, the order of the
  centers in the vector `center_spec` must be the same as in the list of
  the argument `theta_hats`. The used data type in the argument
  `center_spec` must be categorical. Default is `center_spec = NULL`.
  See also ‘Details’ and ‘Examples’.

- theta_A_polys:

  a list with \\L\\ elements so that each element is the array
  `theta_A_poly` (the output of the `MAP.estimation` function,
  `MAP.estimation()$theta_A_poly`) for the corresponding center. This
  argument, `theta_A_polys`, is only used if `family = "survival"` and
  `basehaz = "poly"`. See ‘Details’ and ‘Examples’.

- treat_round:

  a character string representing the `"first"` and `"second"` rounds of
  estimating the treatment effect.

- for_ATE:

  a list of \\L\\ vectors of 9 elements to calculate the average
  treatment effects (ATEs) only for the `binomial` and `gaussian`
  families. These vectors must have equal dimensions. If
  `treat_round = "first"`, then `for_ATE` must be `NULL`. If
  `treat_round = "second"`, then `for_ATE` must be a list for `binomial`
  and `gaussian`, while for `survival`, `for_ATE` must be `NULL`. It
  should be defined using the output of `MAP.estimation()$for_ATE`
  obtained from the second round. See ‘Details’ and ‘Examples’.

- p:

  an integer representing the number of covariates/coefficients. It can
  be found from the output of the `MAP.estimation` function,
  `MAP.estimation()$np`). This argument, `p`, is only used if
  `stratified = TRUE` and `family = "survival"`.

- q_ls:

  a vector with \\L\\ elements in which each element is the order
  (minus 1) of the exponentiated polynomial baseline hazard function for
  the corresponding center, i.e., each element is the value of `q_l`
  (the output of the `MAP.estimation` function, `MAP.estimation()$q_l`).
  This argument, `q_ls`, is only used if `family = "survival"`,
  `family = "survival"` and `basehaz = "poly"`. It can also be a scalar
  which represents the maximum value of the `q_l`'s across the centers.

- center_zero_sample:

  logical flag indicating whether the center has a categorical covariate
  with no observations/individuals in one of the categories. It is used
  to address heterogeneity across centers due to center-specific
  covariates. Default is `center_zero_sample = FALSE`. For more details
  see ‘References’.

- which_cent_zeros:

  an integer vector representing the center(s) which has one categorical
  covariate with no individuals in one of the categories. It is used if
  `center_zero_sample = TRUE`.

- zero_sample_covs:

  a vector in which each element is a character string representing the
  categorical covariate that has no samples/observations in one of its
  categories for the corresponding center. Each element of the vector
  can be obtained from the output of the `MAP.estimation` function for
  the corresponding center, `MAP.estimation()$zero_sample_cov`. It is
  used when `center_zero_sample = TRUE`.

- refer_cats:

  a vector in which each element is a character string representing the
  reference category for the corresponding center. Each element of the
  vector can be obtained from the output of the `MAP.estimation`
  function for the corresponding center, `MAP.estimation()$refer_cat`.
  This vector is used when `center_zero_sample = TRUE`.

- zero_cats:

  a vector in which each element is a character string representing the
  category with no samples/observations for the corresponding center.
  Each element of the vector can be obtained from the output of the
  `MAP.estimation` function for the corresponding center, i.e.,
  `MAP.estimation()$zero_cat`. It is used when
  `center_zero_sample = TRUE`.

- lev_no_ref_zeros:

  a list in which the number of elements equals the length of the
  `which_cent_zeros` argument. Each element of the list is a vector
  containing the names of the levels of the categorical covariate that
  has no samples/observations in one of its categories for the
  corresponding center. However, the name of the category with no
  samples and the name of the reference category are excluded from this
  vector. Each element of the list can be obtained from the output of
  the `MAP.estimation` function, i.e.,
  `MAP.estimation()$lev_no_ref_zero`. This argument is used if
  `center_zero_sample = TRUE`.

## Value

`bfi` returns a list containing the following components:

- theta_hat:

  the vector of estimates obtained by combining the inference results
  from the \\L\\ centers with the `'BFI'` methodology. If an intercept
  was fitted in every center and `stratified = FALSE`, there is only one
  general “intercept” in this vector, while if `stratified = TRUE` and
  `strat_par = 1`, there are \\L\\ different intercepts in the model,
  for each center one. If `treatment` is not '`NULL`', when
  `treat_round = 'first'`, `theta_hat` gives
  \\\hat{\boldsymbol{\gamma}}\_{BFI}\\, and when
  `treat_round = 'second'`, `theta_hat` is the treatment effect \\\hat
  \zeta\_{BFI}\\;

- A_hat:

  minus the curvature (or Hessian) matrix obtained by the `'BFI'` method
  for the combined model. If `stratified = TRUE`, the dimension of the
  matrix is always greater than when `stratified = FALSE`. For the
  `gaussian` family the rows and columns corresponding to the dispersion
  parameter refer to the \\\log(\sigma^2)\\ scale;

- sd:

  the vector of (posterior) standard deviation of the estimates in
  `theta_hat` obtained from the matrix in `A_hat`, i.e., the vector
  equals `sqrt(diag(solve(A_hat)))` which equals the square root of the
  elements at the diagonal of the inverse of the `A_hat` matrix. For the
  `gaussian` family the elements corresponding to the dispersion
  parameter are standard deviations on the \\\log(\sigma^2)\\ scale,
  whereas the estimates in `theta_hat` are reported on the \\\sigma^2\\
  scale.

- family:

  the `family` object used;

- basehaz:

  the baseline hazard function used;

- stratified:

  whether a stratified analysis was done or not;

- strat_par:

  the stratification parameter(s) used;

&nbsp;

- Ave_Treat:

  the estimates of the average treatment effect. Two different
  estimations (IPTW and wIPTW) if the family is `gaussian` or
  `binomial`, and for the `survival` family it is 'NULL'. For more
  details see ‘References’.

## Details

`bfi` function implements the BFI approach described in the papers
Jonker et. al. (2024), Pazira et. al. (2026) and Jonker et. al. (2025)
given in the references. The inference results gathered from different
(\\L\\) centers are combined, and the BFI estimates of the model
parameters and curvature matrix evaluated at that point are returned.

The inference result from each center must be obtained using the
`MAP.estimation` function separately, and then all of these results
(coming from different centers) should be compiled into a list to be
used as an input of `bfi()`. The models in the different centers should
contain the same model parameters. The names of the parameters in each
local `theta_hat` vector and the corresponding row and column names of
each local `A_hat` matrix must refer to the same set of parameters
across centers. The order of these parameters may differ between
centers; in that case, `bfi()` aligns the local parameter vectors,
curvature matrices, and corresponding prior precision matrices by
parameter name before combining the local inference results.

Note that the order of the elements in the lists `theta_hats`, `A_hats`
and `Lambda`, must be the same with respect to the centers, so that in
every list the element at the \\\ell^{\th}\\ position is from the center
\\\ell\\. This should also be the case for the vector `center_spec`.

If for the locations `intercept = FALSE`, the stratified analysis is not
possible anymore for the `binomial` family.

If `stratified = FALSE`, both `strat_par` and `center_spec` must be
`NULL` (the defaults), while if `stratified = TRUE` only one of the two
must be `NULL`.

If `stratified = FALSE` and all the \\L+1\\ matrices in `Lambda` are
equal, it is sufficient to give a (list of) one matrix only. In both
cases of the `stratified` argument (`TRUE` or `FALSE`), if only the
first \\L\\ matrices are equal, the argument `Lambda` can be a list of
two matrices, so that the first matrix represents the chosen
variance-covariance matrix for local centers and the second one is the
chosen matrix for the combined data set. The last matrix of the list in
the argument `Lambda` can be built by the function
[`inv.prior.cov()`](https://hassanpazira.github.io/BFI/reference/inv.prior.cov.md).

If the data type used in the argument `center_spec` is continuous or
categorical with the number of categories equal to the number of
centers, one can use `stratified = TRUE` and `center_spec = NULL`, and
set `strat_par` not to `NULL` (i.e., to \\1\\, \\2\\ or both \\(1,
2)\\). Indeed, in this case, the stratification parameter(s) given in
the argument `strat_par` are assumed to be different across the centers.

When `family = 'survival'` and `basehaz = 'poly'`, the arguments
`theta_hats` and `A_hats` should not be provided. Instead, the
`theta_A_polys` and `q_ls` arguments should be defined using the local
information, specifically `MAP.estimation()$theta_A_poly` and
`MAP.estimation()$q_l`, respectively. See Example 3 in ‘Examples’.

For estimating the treatment effect, in the first round
(`treat_round = "first"`), the argument `for_ATE` must be `NULL` (the
default) and the family must be set to `binomial` (family is handled
automatically.)

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

Hassan Pazira and Marianne Jonker  
Maintainer: Hassan Pazira \<h.pazira@arq.org\>

## See also

[`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)
and
[`inv.prior.cov`](https://hassanpazira.github.io/BFI/reference/inv.prior.cov.md)

## Examples

``` r
#################################################
##  Example 1:  y ~ Binomial  (L = 2 centers)  ##
#################################################

# Setting a seed for reproducibility
set.seed(112358)

#------------------------------------#
# Data Simulation for Local Center 1 #
#------------------------------------#
n1 <- 30                                           # sample size of center 1
X1 <- data.frame(x1=rnorm(n1),                     # continuous variable
                 x2=sample(0:2, n1, replace=TRUE)) # categorical variable
# make dummy variables
X1x2_1 <- ifelse(X1$x2 == '1', 1, 0)
X1x2_2 <- ifelse(X1$x2 == '2', 1, 0)
X1$x2  <- as.factor(X1$x2)
# regression coefficients
beta <- 1:4  # beta[1] is the intercept
# linear predictor:
eta1   <- beta[1] + X1$x1 * beta[2] + X1x2_1 * beta[3] + X1x2_2 * beta[4]
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu1    <- binomial()$linkinv(eta1)
y1     <- rbinom(n1, 1, mu1)

#------------------------------------#
# Data Simulation for Local Center 2 #
#------------------------------------#
n2 <- 50                                           # sample size of center 2
X2 <- data.frame(x1=rnorm(n2),                     # continuous variable
                 x2=sample(0:2, n2, replace=TRUE)) # categorical variable
# make dummy variables:
X2x2_1 <- ifelse(X2$x2 == '1', 1, 0)
X2x2_2 <- ifelse(X2$x2 == '2', 1, 0)
X2$x2  <- as.factor(X2$x2)
# linear predictor:
eta2   <- beta[1] + X2$x1 * beta[2] + X2x2_1 * beta[3] + X2x2_2 * beta[4]
# inverse of the link function:
mu2    <- binomial()$linkinv(eta2)
y2     <- rbinom(n2, 1, mu2)

#---------------------------#
# MAP Estimates at Center 1 #
#---------------------------#
# Assume the same inverse covariance matrix (Lambda) for both centers:
Lambda     <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial')
fit1       <- MAP.estimation(y1, X1, family = 'binomial', Lambda)
theta_hat1 <- fit1$theta_hat # intercept and coefficient estimates
A_hat1     <- fit1$A_hat     # minus the curvature matrix

#---------------------------#
# MAP Estimates at Center 2 #
#---------------------------#
fit2       <- MAP.estimation(y2, X2, family='binomial', Lambda)
theta_hat2 <- fit2$theta_hat
A_hat2     <- fit2$A_hat

#-----------------------#
# BFI at Central Server #
#-----------------------#
theta_hats <- list(theta_hat1, theta_hat2)
A_hats     <- list(A_hat1, A_hat2)
bfi        <- bfi(theta_hats, A_hats, Lambda, family='binomial')
class(bfi)
#> [1] "bfi"
summary(bfi, cur_mat=TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   1.4479  0.8145 -0.1485   3.0444
#> x1            2.4951  0.8978  0.7355   4.2547
#> x21           4.1023  1.6217  0.9239   7.2807
#> x22           2.2199  1.2637 -0.2570   4.6968
#> 
#> Dispersion parameter (sigma2):  1 
#> 
#> Minus the Curvature Matrix: 
#> 
#>             (Intercept)      x1     x21     x22
#> (Intercept)      3.8045 -2.8450  0.7721  0.8694
#> x1              -2.8450  4.0346 -1.2184 -0.5469
#> x21              0.7721 -1.2184  0.7821  0.0000
#> x22              0.8694 -0.5469  0.0000  0.8794

###---------------------###
### Stratified Analysis ###
###---------------------###

# By running the following line an error appears because
# when stratified = TRUE, both 'strat_par' and 'center_spec' can not be NULL:
Just4check1 <- try(bfi(theta_hats, A_hats, Lambda, family = 'binomial',
                   stratified = TRUE), TRUE)
class(Just4check1) # By default, both 'strat_par' and 'center_spec' are NULL!
#> [1] "try-error"

# By running the following line an error appears because when stratified = TRUE,
# last matrix in 'Lambda' should not have the same dim. as the other local matrices:
Just4check2 <- try(bfi(theta_hats, A_hats, Lambda, stratified = TRUE,
                   strat_par = 1), TRUE)
class(Just4check2) # All matices in Lambda have the same dimension!
#> [1] "try-error"

# Stratified analysis when 'intercept' varies across two centers:
newLam <- inv.prior.cov(X1, lambda=c(0.1, 0.3), family = 'binomial',
                        stratified = TRUE, strat_par = 1)
bfi <- bfi(theta_hats, A_hats, list(Lambda, newLam), family = 'binomial',
           stratified=TRUE, strat_par=1)
summary(bfi, cur_mat=TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>                  Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)_loc1   2.3167  2.5315 -2.6448   7.2783
#> (Intercept)_loc2   1.1865  0.7586 -0.3003   2.6733
#> x1                 1.4500  0.7562 -0.0321   2.9322
#> x21                1.9848  1.1901 -0.3479   4.3174
#> x22                1.3621  1.0287 -0.6542   3.3784
#> 
#> Dispersion parameter (sigma2):  1 
#> 
#> Minus the Curvature Matrix: 
#> 
#>                  (Intercept)_loc1 (Intercept)_loc2      x1     x21     x22
#> (Intercept)_loc1           0.1564           0.0000  0.0028  0.0082  0.0134
#> (Intercept)_loc2           0.0000           3.8380 -2.8478  0.7640  0.8561
#> x1                         0.0028          -2.8478  4.3246 -1.2184 -0.5469
#> x21                        0.0082           0.7640 -1.2184  1.0721  0.0000
#> x22                        0.0134           0.8561 -0.5469  0.0000  1.1694


###---------------------###
###  Treatment Effect   ###
###---------------------###

set.seed(112358)

#------------------------------------#
# Data Simulation for Local Center 1 #
#------------------------------------#
n1 <- 30                                           # sample size of center 1
X1 <- data.frame(x1=rnorm(n1),                     # continuous variable
                 treatment=sample(1:2, n1, replace=TRUE)) # categorical variable
X1$treatment  <- as.factor(X1$treatment)

# regression coefficients
beta <- 1:3  # beta[1] is the intercept
# make dummy variable
X1x2_2 <- ifelse(X1$treatment == '2', 1, 0)
# linear predictor:
eta1   <- beta[1] + X1$x1 * beta[2] + X1x2_2 * beta[3]
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu1    <- binomial()$linkinv(eta1)
y1     <- rbinom(n1, 1, mu1)

#------------------------------------#
# Data Simulation for Local Center 2 #
#------------------------------------#
n2 <- 50                                           # sample size of center 2
X2 <- data.frame(x1=rnorm(n2),                     # continuous variable
                 treatment=sample(1:2, n2, replace=TRUE)) # categorical variable
X2$treatment  <- as.factor(X2$treatment)
# make dummy variables:
X2x2_2 <- ifelse(X2$treatment == '2', 1, 0)
# linear predictor:
eta2   <- beta[1] + X2$x1 * beta[2] + X2x2_2 * beta[3]
# inverse of the link function:
mu2    <- binomial()$linkinv(eta2)
y2     <- rbinom(n2, 1, mu2)

# The algorithm works even if the order of the covariates are not
# the same across centers
X2 <- X2[,c("treatment","x1")]

#-----------------------#
#  Observational data   #
#-----------------------#

# For observational data (RWD), we need two rounds for estimating treatment effect:

#-------------#
# First Round #
#-------------#

## Center 1:
Lambda1 <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial',
                         treatment = "treatment", treat_round="first")
fit1_r1 <- MAP.estimation(y1, X1, family = 'binomial', Lambda = Lambda1,
                          treatment = "treatment", treat_round = "first")
# In the first round, the output is without the treatment!
summary(fit1_r1)
#> 
#> Summary of the local model:
#> 
#>    Formula: treatment ~ x1 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)  -0.7226  0.3968 -1.5003   0.0552
#> x1            0.1688  0.3960 -0.6074   0.9450
#> 
#> Dispersion parameter (sigma2):  1 
#>             log Lik Posterior:  -19.01 
#>                   Convergence:  0 

## Center 2:
Lambda2 <- inv.prior.cov(X2, lambda = 0.01, family = 'binomial',
                         treatment = "treatment", treat_round="first")
fit2_r1 <- MAP.estimation(y2, X2, family = 'binomial', Lambda = Lambda2,
                          treatment = "treatment", treat_round = "first")
fit2_r1
#> 
#> Local MAP estimates (family: binomial)
#> 
#> (Intercept)           x1  
#>     -0.2509      -0.1263  
#> 

## Centeral Server:
theta_hats_r1 <- list(fit1_r1$theta_hat, fit2_r1$theta_hat)
A_hats_r1 <- list(fit1_r1$A_hat, fit2_r1$A_hat)
fitbfi_r1 <- bfi(theta_hats_r1, A_hats_r1, Lambda1, family = 'binomial',
                 treat_round = "first")
summary(fitbfi_r1, cur_mat = TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)  -0.3964  0.2299 -0.8470   0.0543
#> x1           -0.0437  0.2464 -0.5266   0.4392
#> 
#> Dispersion parameter (sigma2):  1 
#> 
#> Minus the Curvature Matrix: 
#> 
#>             (Intercept)      x1
#> (Intercept)     18.9205  0.3314
#> x1               0.3314 16.4782

#--------------#
# Second Round #
#--------------#

## Center 1:
Lambda11 <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial',
                          treatment = "treatment", treat_round="second")
fit1_r2 <- MAP.estimation(y1, X1, family = 'binomial', Lambda = Lambda11,
                          treatment = "treatment", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)
# In the second round, the output is only with the treatment!
summary(fit1_r2)
#> 
#> Summary of the local model:
#> 
#>    Formula: treatment ~ x1 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.6013  0.3613 -0.1068   1.3095
#> treatment     0.8123  0.6188 -0.4005   2.0251
#> 
#> Dispersion parameter (sigma2):  1 
#>             log Lik Posterior:  -34.09 
#>                   Convergence:  0 

## Center 2:
Lambda22 <- inv.prior.cov(X2, lambda = 0.01, family = 'binomial',
                         treatment = "treatment", treat_round="second")
fit2_r2 <- MAP.estimation(y2, X2, family = 'binomial', Lambda = Lambda22,
                          treatment = "treatment", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)
fit2_r2$propensity # Propensity Score
#>  [1] 0.4011608 0.4046425 0.4016316 0.4065235 0.4013848 0.3921165 0.4006005
#>  [8] 0.4171208 0.3974666 0.3818413 0.4161945 0.4135265 0.4065453 0.4080270
#> [15] 0.4035998 0.4049694 0.3868497 0.4058360 0.4051789 0.4120526 0.4130087
#> [22] 0.3800897 0.4107172 0.4150832 0.3942670 0.4036748 0.4111251 0.4017086
#> [29] 0.4014951 0.3879404 0.3867900 0.4146683 0.4038661 0.3964183 0.4211174
#> [36] 0.4029404 0.4135190 0.4144179 0.3932297 0.4003255 0.4077919 0.3905918
#> [43] 0.4088652 0.4097657 0.3931901 0.3913643 0.3953093 0.4043939 0.4066814
#> [50] 0.4074400
fit2_r2$for_ATE # will be used in central server
#> [1] 22.00000 28.00000 22.00000 22.00000 54.53770 54.53770 46.87703 34.95445
#> [9] 21.00000
fit2_r2
#> 
#> Local MAP estimates (family: binomial)
#> 
#> (Intercept)    treatment  
#>       1.081        5.769  
#> 

## Centeral Server:
theta_hats_r2 <- list(fit1_r2$theta_hat, fit2_r2$theta_hat)
A_hats_r2 <- list(fit1_r2$A_hat, fit2_r2$A_hat)
for_ATEs <- list(fit1_r2$for_ATE, fit2_r2$for_ATE)
fitbfi_r2 <- bfi(theta_hats_r2, A_hats_r2, Lambda11, family = 'binomial',
                 treat_round = "second", for_ATE = for_ATEs)
fitbfi_r2$S_var
#> NULL
fitbfi_r2$Ave_Treat
#> $IPTW
#> [1] 0.2269916
#> 
#> $wIPTW
#> [1] 0.234374
#> 
summary(fitbfi_r2)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.8558  0.2460  0.3737   1.3379
#> treatment     0.6510  0.5564 -0.4395   1.7416
#> 
#> Dispersion parameter (sigma2):  1 
#> 
#> Average Treatment Effect (ATE):  
#> 
#>          IPTW:  0.227 
#>         wIPTW:  0.2344 

#--------------------#
#  Randomized Trial  #
#--------------------#

# For Randomized Control Trial (RCT), we need only one round (the second round) for
# estimating treatment effect. Because we do not need to estimate propensity score.
# For example, in a 1:1 randomized trial, the propensity scores are, by definition,
# equal to 0.5. Here we use 'RCT_propens', instead of 'gamma_bfi':

## Center 1:
Lambda11 <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial',
                          treatment = "treatment", treat_round="second")
fit1_r2 <- MAP.estimation(y1, X1, family = 'binomial', Lambda = Lambda11,
                          treatment = "treatment", treat_round = "second",
                          RCT_propens = rep(0.5, n1)) # gamma_bfi = NULL
summary(fit1_r2)
#> 
#> Summary of the local model:
#> 
#>    Formula: treatment ~ x1 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.6192  0.3311 -0.0298   1.2682
#> treatment     0.7647  0.6481 -0.5056   2.0350
#> 
#> Dispersion parameter (sigma2):  1 
#>             log Lik Posterior:  -35.91 
#>                   Convergence:  0 

## Center 2:
Lambda22 <- inv.prior.cov(X2, lambda = 0.01, family = 'binomial',
                          treatment = "treatment", treat_round="second")
fit2_r2 <- MAP.estimation(y2, X2, family = 'binomial', Lambda = Lambda22,
                          treatment = "treatment", treat_round = "second",
                          RCT_propens = rep(0.5, n2)) # gamma_bfi = NULL
fit2_r2$for_ATE # will be used in central server
#> [1] 22 28 22 22 44 44 56 42 21
fit2_r2
#> 
#> Local MAP estimates (family: binomial)
#> 
#> (Intercept)    treatment  
#>       1.103        5.568  
#> 

## Centeral Server:
theta_hats_r2 <- list(fit1_r2$theta_hat, fit2_r2$theta_hat)
A_hats_r2 <- list(fit1_r2$A_hat, fit2_r2$A_hat)
for_ATEs <- list(fit1_r2$for_ATE, fit2_r2$for_ATE)
fitbfi_r2 <- bfi(theta_hats_r2, A_hats_r2, Lambda11, family = 'binomial',
                 treat_round = "second", for_ATE = for_ATEs)
fitbfi_r2$S_var
#> NULL
fitbfi_r2$Ave_Treat
#> $IPTW
#> [1] -0.1
#> 
#> $wIPTW
#> [1] 0.2291667
#> 
summary(fitbfi_r2)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.8756  0.2259  0.4328   1.3183
#> treatment     0.6161  0.5971 -0.5542   1.7863
#> 
#> Dispersion parameter (sigma2):  1 
#> 
#> Average Treatment Effect (ATE):  
#> 
#>          IPTW:  -0.1 
#>         wIPTW:  0.2292 


#################################################
##  Example 2:  y ~ Gaussian  (L = 3 centers)  ##
#################################################

# Setting a seed for reproducibility
set.seed(112358)

p     <- 3                     # number of coefficients without 'intercept'
theta <- c(1, rep(2, p), 1.5)  # reg. coef.s ('intercept' is 1) & 'sigma2' = 1.5

#------------------------------------#
# Data Simulation for Local Center 1 #
#------------------------------------#
n1   <- 30                                       # sample size of center 1
X1   <- data.frame(matrix(rnorm(n1 * p), n1, p)) # continuous variables
# linear predictor:
eta1 <- theta[1] + as.matrix(X1) %*% theta[2:4]
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu1  <- gaussian()$linkinv(eta1)
y1   <- rnorm(n1, mu1, sd = sqrt(theta[5]))

#------------------------------------#
# Data Simulation for Local Center 2 #
#------------------------------------#
n2   <- 40                                       # sample size of center 2
X2   <- data.frame(matrix(rnorm(n2 * p), n2, p)) # continuous variables
# linear predictor:
eta2 <- theta[1] + as.matrix(X2) %*% theta[2:4]
# inverse of the link function:
mu2  <- gaussian()$linkinv(eta2)
y2   <- rnorm(n2, mu2, sd = sqrt(theta[5]))

#------------------------------------#
# Data Simulation for Local Center 3 #
#------------------------------------#
n3   <- 50                                       # sample size of center 3
X3   <- data.frame(matrix(rnorm(n3 * p), n3, p)) # continuous variables
# linear predictor:
eta3 <- theta[1] + as.matrix(X3) %*% theta[2:4]
# inverse of the link function:
mu3  <- gaussian()$linkinv(eta3)
y3   <- rnorm(n3, mu3, sd = sqrt(theta[5]))

#---------------------------#
# Inverse Covariance Matrix #
#---------------------------#
# Creating the inverse covariance matrix for the Gaussian prior distribution:
# the same for both centers
Lambda <- inv.prior.cov(X1, lambda = 0.05, family='gaussian')

#---------------------------#
# MAP Estimates at Center 1 #
#---------------------------#
fit1       <- MAP.estimation(y1, X1, family = 'gaussian', Lambda)
theta_hat1 <- fit1$theta_hat # intercept and coefficient estimates
A_hat1     <- fit1$A_hat     # minus the curvature matrix

#---------------------------#
# MAP Estimates at Center 2 #
#---------------------------#
fit2       <- MAP.estimation(y2, X2, family = 'gaussian', Lambda)
theta_hat2 <- fit2$theta_hat
A_hat2     <- fit2$A_hat

#---------------------------#
# MAP Estimates at Center 3 #
#---------------------------#
fit3       <- MAP.estimation(y3, X3, family = 'gaussian', Lambda)
theta_hat3 <- fit3$theta_hat
A_hat3     <- fit3$A_hat

#-----------------------#
# BFI at Central Server #
#-----------------------#
A_hats     <- list(A_hat1, A_hat2, A_hat3)
theta_hats <- list(theta_hat1, theta_hat2, theta_hat3)
bfi        <- bfi(theta_hats, A_hats, Lambda, family = 'gaussian')
summary(bfi, cur_mat=TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.8666  0.1008  0.6691   1.0641
#> X1            1.9381  0.1010  1.7401   2.1360
#> X2            1.9914  0.0969  1.8015   2.1814
#> X3            2.0190  0.1001  1.8228   2.2152
#> 
#> Dispersion parameter (sigma2):  1.178 
#> 
#> Minus the Curvature Matrix: 
#> 
#>             (Intercept)      X1       X2       X3  sigma2
#> (Intercept)    103.1343 -1.6101  17.4699 -13.7999  0.1268
#> X1              -1.6101 98.6759   7.0585  -2.9108  0.2902
#> X2              17.4699  7.0585 109.9640  -1.7586  0.2952
#> X3             -13.7999 -2.9108  -1.7586 101.7741  0.3053
#> sigma2           0.1268  0.2902   0.2952   0.3053 60.0788

###---------------------###
### Stratified Analysis ###
###---------------------###

# Stratified analysis when 'intercept' varies across two centers:
newLam1 <- inv.prior.cov(X1, lambda = c(0.1,0.3), family = 'gaussian',
                         stratified = TRUE, strat_par = 1, L = 3)
# 'newLam1' is used as the prior for combined data and
# 'Lambda' is used as the prior for locals
list_newLam1 <- list(Lambda, newLam1)
bfi1 <- bfi(theta_hats, A_hats, list_newLam1, family = 'gaussian',
            stratified = TRUE, strat_par = 1)
summary(bfi1, cur_mat = TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>                  Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)_loc1   0.7078  0.2132  0.2899   1.1257
#> (Intercept)_loc2   0.9409  0.1681  0.6115   1.2704
#> (Intercept)_loc3   0.8822  0.1539  0.5805   1.1839
#> X1                 1.9483  0.1022  1.7480   2.1485
#> X2                 1.9842  0.0973  1.7935   2.1748
#> X3                 2.0255  0.1016  1.8264   2.2246
#> 
#> Dispersion parameter (sigma2):  1.173 
#> 
#> Minus the Curvature Matrix: 
#> 
#>                  (Intercept)_loc1 (Intercept)_loc2 (Intercept)_loc3      X1
#> (Intercept)_loc1          22.1351           0.0000           0.0000  3.3624
#> (Intercept)_loc2           0.0000          38.6210           0.0000 -7.1784
#> (Intercept)_loc3           0.0000           0.0000          42.6283  2.2058
#> X1                         3.3624          -7.1784           2.2058 98.7259
#> X2                         1.1876           9.8372           6.4451  7.0585
#> X3                        -1.3801         -13.1665           0.7467 -2.9108
#> sigma2                     0.0362           0.0466           0.0440  0.2902
#>                        X2       X3  sigma2
#> (Intercept)_loc1   1.1876  -1.3801  0.0362
#> (Intercept)_loc2   9.8372 -13.1665  0.0466
#> (Intercept)_loc3   6.4451   0.7467  0.0440
#> X1                 7.0585  -2.9108  0.2902
#> X2               110.0140  -1.7586  0.2952
#> X3                -1.7586 101.8241  0.3053
#> sigma2             0.2952   0.3053 60.3288

# Stratified analysis when 'sigma2' varies across two centers:
newLam2 <- inv.prior.cov(X1, lambda = c(0.1,0.3), family = 'gaussian',
                         stratified = TRUE, strat_par = 2, L = 3)
# 'newLam2' is used as the prior for combined data and 'Lambda' is used as
# the prior for locals
list_newLam2 <- list(Lambda, newLam2)
bfi2 <- bfi(theta_hats, A_hats, list_newLam2, family = 'gaussian',
            stratified = TRUE, strat_par=2)
summary(bfi2, cur_mat = TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.8662  0.1008  0.6687   1.0637
#> X1            1.9371  0.1010  1.7392   2.1350
#> X2            1.9907  0.0969  1.8007   2.1806
#> X3            2.0179  0.1001  1.8218   2.2140
#> sigma2_loc1   1.3388  0.3421  0.8114   2.2091
#> sigma2_loc2   1.0257  0.2276  0.6639   1.5847
#> sigma2_loc3   1.1641  0.2314  0.7885   1.7187
#> 
#> For the residual variances, the standard deviation is obtained by the 
#> delta method and the credible interval is computed on the log scale 
#> and back-transformed, so the interval is not symmetric.
#> 
#> Minus the Curvature Matrix: 
#> 
#>             (Intercept)      X1       X2       X3 sigma2_loc1 sigma2_loc2
#> (Intercept)    103.1843 -1.6101  17.4699 -13.7999      0.0362      0.0466
#> X1              -1.6101 98.7259   7.0585  -2.9108      0.0954      0.0989
#> X2              17.4699  7.0585 110.0140  -1.7586      0.0957      0.0989
#> X3             -13.7999 -2.9108  -1.7586 101.8241      0.1052      0.0987
#> sigma2_loc1      0.0362  0.0954   0.0957   0.1052     15.3181      0.0000
#> sigma2_loc2      0.0466  0.0989   0.0989   0.0987      0.0000     20.3019
#> sigma2_loc3      0.0440  0.0958   0.1006   0.1013      0.0000      0.0000
#>             sigma2_loc3
#> (Intercept)      0.0440
#> X1               0.0958
#> X2               0.1006
#> X3               0.1013
#> sigma2_loc1      0.0000
#> sigma2_loc2      0.0000
#> sigma2_loc3     25.3088

# Stratified analysis when 'intercept' and 'sigma2' vary across 2 centers:
newLam3 <- inv.prior.cov(X1, lambda = c(0.1,0.2,0.3), family = 'gaussian',
                         stratified = TRUE, strat_par = c(1, 2), L = 3)
# 'newLam3' is used as the prior for combined data and 'Lambda' is used as
# the prior for locals
list_newLam3 <- list(Lambda, newLam3)
bfi3 <- bfi(theta_hats, A_hats, list_newLam3, family = 'gaussian',
            stratified = TRUE, strat_par = 1:2)
summary(bfi3, cur_mat = TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>                  Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)_loc1   0.7077  0.2132  0.2898   1.1256
#> (Intercept)_loc2   0.9404  0.1681  0.6110   1.2699
#> (Intercept)_loc3   0.8826  0.1539  0.5809   1.1843
#> X1                 1.9463  0.1021  1.7462   2.1464
#> X2                 1.9825  0.0972  1.7919   2.1730
#> X3                 2.0234  0.1015  1.8244   2.2224
#> sigma2_loc1        1.3392  0.3422  0.8116   2.2096
#> sigma2_loc2        1.0255  0.2276  0.6638   1.5844
#> sigma2_loc3        1.1641  0.2314  0.7885   1.7187
#> 
#> For the residual variances, the standard deviation is obtained by the 
#> delta method and the credible interval is computed on the log scale 
#> and back-transformed, so the interval is not symmetric.
#> 
#> Minus the Curvature Matrix: 
#> 
#>                  (Intercept)_loc1 (Intercept)_loc2 (Intercept)_loc3      X1
#> (Intercept)_loc1          22.1351           0.0000           0.0000  3.3624
#> (Intercept)_loc2           0.0000          38.6210           0.0000 -7.1784
#> (Intercept)_loc3           0.0000           0.0000          42.6283  2.2058
#> X1                         3.3624          -7.1784           2.2058 98.8259
#> X2                         1.1876           9.8372           6.4451  7.0585
#> X3                        -1.3801         -13.1665           0.7467 -2.9108
#> sigma2_loc1                0.0362           0.0000           0.0000  0.0954
#> sigma2_loc2                0.0000           0.0466           0.0000  0.0989
#> sigma2_loc3                0.0000           0.0000           0.0440  0.0958
#>                        X2       X3 sigma2_loc1 sigma2_loc2 sigma2_loc3
#> (Intercept)_loc1   1.1876  -1.3801      0.0362      0.0000      0.0000
#> (Intercept)_loc2   9.8372 -13.1665      0.0000      0.0466      0.0000
#> (Intercept)_loc3   6.4451   0.7467      0.0000      0.0000      0.0440
#> X1                 7.0585  -2.9108      0.0954      0.0989      0.0958
#> X2               110.1140  -1.7586      0.0957      0.0989      0.1006
#> X3                -1.7586 101.9241      0.1052      0.0987      0.1013
#> sigma2_loc1        0.0957   0.1052     15.3181      0.0000      0.0000
#> sigma2_loc2        0.0989   0.0987      0.0000     20.3019      0.0000
#> sigma2_loc3        0.1006   0.1013      0.0000      0.0000     25.3088

###----------------------------###
### Center Specific Covariates ###
###----------------------------###

# Assume the first and third centers have the same center-specific covariate value
# of 'High', while this value for the second center is 'Low', i.e.,
# center_spec = c('High','Low','High')
newLam4 <- inv.prior.cov(X1, lambda=c(0.1, 0.2, 0.3), family='gaussian',
                         stratified = TRUE, center_spec = c('High','Low','High'),
                         L = 3)
# 'newLam4' is used as the prior for combined data and 'Lambda' is used as
# the prior for locals
l_newLam4 <- list(Lambda, newLam4)
bfi4 <- bfi(theta_hats, A_hats, l_newLam4, family = 'gaussian',
            stratified = TRUE, center_spec = c('High','Low','High'))
summary(bfi4, cur_mat = TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>                  Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)_High   0.8242  0.1251  0.5789   1.0694
#> (Intercept)_Low    0.9398  0.1681  0.6104   1.2692
#> X1                 1.9435  0.1020  1.7435   2.1434
#> X2                 1.9849  0.0971  1.7945   2.1753
#> X3                 2.0253  0.1015  1.8264   2.2242
#> 
#> Dispersion parameter (sigma2):  1.173 
#> 
#> Minus the Curvature Matrix: 
#> 
#>                  (Intercept)_High (Intercept)_Low      X1       X2       X3
#> (Intercept)_High          64.6634          0.0000  5.5682   7.6327  -0.6334
#> (Intercept)_Low            0.0000         38.6210 -7.1784   9.8372 -13.1665
#> X1                         5.5682         -7.1784 98.8259   7.0585  -2.9108
#> X2                         7.6327          9.8372  7.0585 110.1140  -1.7586
#> X3                        -0.6334        -13.1665 -2.9108  -1.7586 101.9241
#> sigma2                     0.0802          0.0466  0.2902   0.2952   0.3053
#>                   sigma2
#> (Intercept)_High  0.0802
#> (Intercept)_Low   0.0466
#> X1                0.2902
#> X2                0.2952
#> X3                0.3053
#> sigma2           60.3288


###---------------------###
###  Treatment Effect   ###
###---------------------###

set.seed(112358)

#-----------------------------#
# New Data for Local Center 1 #
#-----------------------------#
# Generating new data with 'treatment' variable
# We cansider the first variable (X1$X1) to be the treatment:
X1$X1 <- sample(0:1, n1, replace=TRUE) # categorical variable
eta1  <- theta[1] + as.matrix(X1) %*% theta[2:4]
mu1   <- gaussian()$linkinv(eta1)
y1    <- rnorm(n1, mu1, sd = sqrt(theta[5]))

#-----------------------------#
# New Data for Local Center 2 #
#-----------------------------#
# We cansider the first variable (X2$X1) to be the treatment:
X2$X1 <- sample(0:1, n2, replace=TRUE) # categorical variable
eta2  <- theta[1] + as.matrix(X2) %*% theta[2:4]
mu2   <- gaussian()$linkinv(eta2)
y2    <- rnorm(n2, mu2, sd = sqrt(theta[5]))

#-----------------------------#
# New Data for Local Center 3 #
#-----------------------------#
# We cansider the first variable (X3$X1) to be the treatment:
X3$X1 <- sample(0:1, n3, replace=TRUE) # categorical variable
# linear predictor:
eta3  <- theta[1] + as.matrix(X3) %*% theta[2:4]
# inverse of the link function:
mu3   <- gaussian()$linkinv(eta3)
y3    <- rnorm(n3, mu3, sd = sqrt(theta[5]))

#-----------------------#
#  Observational data   #
#-----------------------#

# For observational data (RWD), we need two rounds for estimating treatment effect:

#-------------#
# First Round #
#-------------#

## Center 1:
Lambda1 <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial',
                         treatment = "X1", treat_round="first")
# When treat_round = "first", the family will automatically set to 'binomial',
# even if family = 'gaussian' or family = 'survival'.
fit1_r1 <- MAP.estimation(y1, X1, family = 'gaussian', Lambda = Lambda1,
                          treatment = "X1", treat_round = "first")
# Althghou family = 'gaussian', the output is based on 'binomial'!
# The output without the treatment (X1) in the first round!
summary(fit1_r1)
#> 
#> Summary of the local model:
#> 
#>    Formula: X1 ~ X2 + X3 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.1108  0.3732 -0.6206   0.8423
#> X2           -0.0297  0.3798 -0.7742   0.7148
#> X3           -0.4071  0.4199 -1.2301   0.4158
#> 
#> Dispersion parameter (sigma2):  1 
#>             log Lik Posterior:  -20.24 
#>                   Convergence:  0 

## Center 2:
Lambda2 <- inv.prior.cov(X2, lambda = 0.01, family = 'gaussian',
                         treatment = "X1", treat_round="first")
fit2_r1 <- MAP.estimation(y2, X2, family = 'gaussian', Lambda = Lambda2,
                          treatment = "X1", treat_round = "first")
fit2_r1
#> 
#> Local MAP estimates (family: binomial)
#> 
#> (Intercept)           X2           X3  
#>      0.3157      -0.5803       0.1921  
#> 

## Center 3:
Lambda3 <- inv.prior.cov(X3, lambda = 0.01, family = 'gaussian',
                         treatment = "X1", treat_round="first")
fit3_r1 <- MAP.estimation(y3, X3, family = 'gaussian', Lambda = Lambda3,
                          treatment = "X1", treat_round = "first")

## Centeral Server:
theta_hats_r1 <- list(fit1_r1$theta_hat, fit2_r1$theta_hat, fit3_r1$theta_hat)
A_hats_r1 <- list(fit1_r1$A_hat, fit2_r1$A_hat, fit3_r1$A_hat)
fitbfi_r1 <- bfi(theta_hats_r1, A_hats_r1, Lambda1, family = 'gaussian',
                 treat_round = "first") # same results with 'binomial'
# The output without the treatment (X1) in the first round!
summary(fitbfi_r1, cur_mat = TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.1594  0.1904 -0.2138   0.5325
#> X2           -0.1777  0.1896 -0.5494   0.1939
#> X3           -0.0442  0.1943 -0.4251   0.3367
#> 
#> Dispersion parameter (sigma2):  1 
#> 
#> Minus the Curvature Matrix: 
#> 
#>             (Intercept)      X2      X3
#> (Intercept)     28.7143  4.6181 -3.0898
#> X2               4.6181 28.6170  0.7706
#> X3              -3.0898  0.7706 26.8660

#--------------#
# Second Round #
#--------------#

## Center 1:
Lambda11 <- inv.prior.cov(X1, lambda = 0.01, family = 'gaussian',
                          treatment = "X1", treat_round="second")
fit1_r2 <- MAP.estimation(y1, X1, family = 'gaussian', Lambda = Lambda11,
                          treatment = "X1", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)
# The output with only the treatment (X1) in the second round!
summary(fit1_r2)
#> 
#> Summary of the local model:
#> 
#>    Formula: X1 ~ X2 + X3 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   1.0775  0.6320 -0.1611   2.3161
#> X1            1.6527  0.8983 -0.1080   3.4133
#> 
#> Dispersion parameter (sigma2):  12.28 
#>             log Lik Posterior:  -105.5 
#>                   Convergence:  0 

## Center 2:
Lambda22 <- inv.prior.cov(X2, lambda = 0.01, family = 'gaussian', treatment = "X1",
                          treat_round="second")
fit2_r2 <- MAP.estimation(y2, X2, family = 'gaussian', Lambda = Lambda22,
                          treatment = "X1", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)

## Center 3:
Lambda33 <- inv.prior.cov(X3, lambda = 0.01, family = 'gaussian', treatment = "X1",
                          treat_round="second")
fit3_r2 <- MAP.estimation(y3, X3, family = 'gaussian', Lambda = Lambda33,
                          treatment = "X1", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)

## Central Server:
theta_hats_r2 <- list(fit1_r2$theta_hat, fit2_r2$theta_hat, fit3_r2$theta_hat)
A_hats_r2 <- list(fit1_r2$A_hat, fit2_r2$A_hat, fit3_r2$A_hat)
for_ATEs <- list(fit1_r2$for_ATE, fit2_r2$for_ATE, fit3_r2$for_ATE)
fitbfi_r2 <- bfi(theta_hats_r2, A_hats_r2, Lambda11, family = 'gaussian',
                 treat_round = "second", for_ATE = for_ATEs)
fitbfi_r2$Ave_Treat
#> $IPTW
#> [1] 2.040845
#> 
#> $wIPTW
#> [1] 2.054435
#> 
fitbfi_r2$S_var
#> NULL
summary(fitbfi_r2)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘gaussian’ 
#>       Link: ‘identity’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.9774  0.3461  0.2990   1.6557
#> X1            2.0672  0.4901  1.1067   3.0277
#> 
#> Dispersion parameter (sigma2):  14.69 
#> 
#> Average Treatment Effect (ATE):  
#> 
#>          IPTW:  2.041 
#>         wIPTW:  2.054 


####################################################
##  Example 3:  Survival family  (L = 2 centers)  ##
####################################################

# Setting a seed for reproducibility
set.seed(112358)

p <- 3
theta <- c(1:4, 5, 6)  # regression coefficients (1:4) & omega's (5:6)

#---------------------------------------------#
# Simulating Survival data for Local Center 1 #
#---------------------------------------------#
n1 <- 30
X1 <- data.frame(matrix(rnorm(n1 * p), n1, p)) # continuous (normal) variables
# Simulating survival data ('time' and 'status') from 'Weibull' with
# a predefined censoring rate of 0.3:
y1 <- surv.simulate(Z = list(X1), beta = theta[1:p], a = theta[5],
                    b = theta[6], u1 = 0.1, cen_rate = 0.3,
                    gen_data_from = "weibul")$D[[1]][, 1:2]

## MAP Estimates at Center 1
Lambda <- inv.prior.cov(X1, lambda = c(0.1, 1), family = "survival",
                        basehaz = "poly")
fit1 <- MAP.estimation(y1, X1, family = 'survival', Lambda = Lambda,
                       basehaz = "poly")
theta_hat1 <- fit1$theta_hat  # coefficient estimates
A_hat1     <- fit1$A_hat      # minus the curvature matrix
summary(fit1, cur_mat=TRUE)
#> 
#> Summary of the local model:
#> 
#>    Formula: Survival(time, status) ~ X1 + X2 + X3 
#>     Family: ‘survival’ 
#>   Baseline: ‘poly’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   0.6589  0.2954  0.0801   1.2378
#> X2   0.8425  0.2474  0.3575   1.3274
#> X3   1.4327  0.2974  0.8498   2.0156
#> 
#> log Lik Posterior:  -1.275 
#>       Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>              X1      X2      X3 omega_0 omega_1 omega_2
#> X1      16.4169 -7.9165 -1.8776  2.0008  0.4276 -0.1580
#> X2      -7.9165 23.3009  1.3864  4.1509  0.4604  0.3718
#> X3      -1.8776  1.3864 21.2546  3.4463 -3.3272 -4.4365
#> omega_0  2.0008  4.1509  3.4463 19.5158  9.5086  7.3123
#> omega_1  0.4276  0.4604 -3.3272  9.5086  8.3123  6.5539
#> omega_2 -0.1580  0.3718 -4.4365  7.3123  6.5539  7.3698
fit1$theta_A_poly # Only when family = "survival" and basehaz ="poly"
#> , , 1
#> 
#>         [,1] [,2]       [,3] [,4] [,5] [,6]
#> X1        NA   NA  0.6589476   NA   NA   NA
#> X2        NA   NA  0.8424763   NA   NA   NA
#> X3        NA   NA  1.4327460   NA   NA   NA
#> omega_0   NA   NA -1.5157831   NA   NA   NA
#> omega_1   NA   NA  1.8741816   NA   NA   NA
#> omega_2   NA   NA  1.7874313   NA   NA   NA
#> 
#> , , 2
#> 
#>         [,1] [,2] [,3] [,4] [,5] [,6]
#> X1        NA   NA   NA   NA   NA   NA
#> X2        NA   NA   NA   NA   NA   NA
#> X3        NA   NA   NA   NA   NA   NA
#> omega_0   NA   NA   NA   NA   NA   NA
#> omega_1   NA   NA   NA   NA   NA   NA
#> omega_2   NA   NA   NA   NA   NA   NA
#> 
#> , , 3
#> 
#>         [,1] [,2] [,3] [,4] [,5] [,6]
#> X1        NA   NA   NA   NA   NA   NA
#> X2        NA   NA   NA   NA   NA   NA
#> X3        NA   NA   NA   NA   NA   NA
#> omega_0   NA   NA   NA   NA   NA   NA
#> omega_1   NA   NA   NA   NA   NA   NA
#> omega_2   NA   NA   NA   NA   NA   NA
#> 
#> , , 4
#> 
#>               [,1]       [,2]      [,3]      [,4]       [,5]       [,6]
#> X1      16.4168965 -7.9165022 -1.877557  2.000794  0.4276326 -0.1580233
#> X2      -7.9165022 23.3009020  1.386382  4.150947  0.4603776  0.3717607
#> X3      -1.8775567  1.3863818 21.254639  3.446257 -3.3271931 -4.4364945
#> omega_0  2.0007944  4.1509475  3.446257 19.515762  9.5086145  7.3122662
#> omega_1  0.4276326  0.4603776 -3.327193  9.508614  8.3122662  6.5539443
#> omega_2 -0.1580233  0.3717607 -4.436494  7.312266  6.5539443  7.3698101
#> 

#---------------------------------------------#
# Simulating Survival data for Local Center 2 #
#---------------------------------------------#
n2 <- 30
X2 <- data.frame(matrix(rnorm(n2 * p), n2, p)) # continuous (normal) variables
# Survival simulated data from 'Weibull' with a predefined censoring rate of 0.3:
y2 <- surv.simulate(Z = list(X2), beta = theta[1:p], a = theta[5],
                    b = theta[6],u1 = 0.1, cen_rate = 0.3,
                    gen_data_from = "weibul")$D[[1]][, 1:2]

## MAP Estimates at Center 2
fit2 <- MAP.estimation(y2, X2, family = 'survival', Lambda = Lambda,
                       basehaz = "poly")
theta_hat2 <- fit2$theta_hat
A_hat2 <- fit2$A_hat
summary(fit2, cur_mat=TRUE)
#> 
#> Summary of the local model:
#> 
#>    Formula: Survival(time, status) ~ X1 + X2 + X3 
#>     Family: ‘survival’ 
#>   Baseline: ‘poly’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   0.5961  0.2632  0.0803   1.1120
#> X2   0.4218  0.1997  0.0303   0.8132
#> X3   1.2563  0.2957  0.6768   1.8359
#> 
#> log Lik Posterior:  -6.197 
#>       Convergence:  0 
#> 
#> Minus the Curvature Matrix: 
#> 
#>              X1      X2      X3 omega_0 omega_1 omega_2
#> X1      16.8039  1.7314  2.0927  0.5517 -3.0791 -6.2397
#> X2       1.7314 28.4288  4.3588  1.2219 -3.6742 -4.7424
#> X3       2.0927  4.3588 22.4192  7.9651 -2.9422 -8.3513
#> omega_0  0.5517  1.2219  7.9651 19.0430  8.5748  8.9323
#> omega_1 -3.0791 -3.6742 -2.9422  8.5748  9.9323 13.7961
#> omega_2 -6.2397 -4.7424 -8.3513  8.9323 13.7961 26.8460

#-----------------------#
# BFI at Central Server #
#-----------------------#
# When family = 'survival' and basehaz = "poly", only 'theta_A_polys'
# should be defined instead of 'theta_hats' and 'A_hats':
theta_A_hats <- list(fit1$theta_A_poly, fit2$theta_A_poly)
qls <- c(fit1$q_l, fit2$q_l)
bfi <- bfi(Lambda = Lambda, family = 'survival', theta_A_polys = theta_A_hats,
           basehaz = "poly", q_ls = qls)
summary(bfi, cur_mat=TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘survival’ 
#>   Baseline: ‘poly’
#> 
#> Coefficients:
#> 
#>         Estimate Std.Dev CI 2.5% CI 97.5%
#> X1        0.6377  0.1918  0.2618   1.0135
#> X2        0.7019  0.1531  0.4018   1.0019
#> X3        1.3883  0.2131  0.9705   1.8060
#> omega_0  -1.6270  0.3536 -2.3201  -0.9339
#> omega_1   3.6309  0.7166  2.2263   5.0354
#> omega_2  -0.4161  0.3617 -1.1250   0.2928
#> 
#> Minus the Curvature Matrix: 
#> 
#>              X1      X2       X3 omega_0 omega_1  omega_2
#> X1      33.1208 -6.1851   0.2152  2.5525 -2.6515  -6.3977
#> X2      -6.1851 51.6297   5.7451  5.3728 -3.2138  -4.3706
#> X3       0.2152  5.7451  43.5738 11.4113 -6.2694 -12.7878
#> omega_0  2.5525  5.3728  11.4113 37.5588 18.0834  16.2445
#> omega_1 -2.6515 -3.2138  -6.2694 18.0834 17.2445  20.3501
#> omega_2 -6.3977 -4.3706 -12.7878 16.2445 20.3501  33.2158


###---------------------###
### Stratified Analysis ###
###---------------------###

# Stratified analysis when first parameter ('omega_0') varies across two centers:
(newLam0 <- inv.prior.cov(X1, lambda = c(rep(1, 3), 0.3, 0.7, rep(2,2)),
                          family = 'survival', stratified = TRUE,
                          basehaz = c("poly"), strat_par = 1, L = 2))
#>              X1 X2 X3 omega_0_loc1 omega_0_loc2 omega_1 omega_2
#> X1            1  0  0          0.0          0.0       0       0
#> X2            0  1  0          0.0          0.0       0       0
#> X3            0  0  1          0.0          0.0       0       0
#> omega_0_loc1  0  0  0          0.3          0.0       0       0
#> omega_0_loc2  0  0  0          0.0          0.7       0       0
#> omega_1       0  0  0          0.0          0.0       2       0
#> omega_2       0  0  0          0.0          0.0       0       2
# 'newLam0' is used as the prior for combined data and 'Lambda' is used as for locals:
list_newLam0 <- list(Lambda, newLam0)
bfi0 <- bfi(Lambda = list_newLam0, family = 'survival', theta_A_polys = theta_A_hats,
            stratified = TRUE, basehaz = c("poly"), p = 3, q_ls = qls, strat_par = 1)
summary(bfi0, cur_mat = TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘survival’ 
#>   Baseline: ‘poly’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   0.5408  0.1871  0.1741   0.9076
#> X2   0.5965  0.1496  0.3032   0.8897
#> X3   1.1933  0.2034  0.7947   1.5919
#> 
#> Minus the Curvature Matrix: 
#> 
#>                   X1      X2       X3 omega_0_loc1 omega_0_loc2 omega_1
#> X1           34.0208 -6.1851   0.2152       2.0008       0.5517 -2.6515
#> X2           -6.1851 52.5297   5.7451       4.1509       1.2219 -3.2138
#> X3            0.2152  5.7451  44.4738       3.4463       7.9651 -6.2694
#> omega_0_loc1  2.0008  4.1509   3.4463      18.8158       0.0000  9.5086
#> omega_0_loc2  0.5517  1.2219   7.9651       0.0000      18.7430  8.5748
#> omega_1      -2.6515 -3.2138  -6.2694       9.5086       8.5748 18.2445
#> omega_2      -6.3977 -4.3706 -12.7878       7.3123       8.9323 20.3501
#>               omega_2
#> X1            -6.3977
#> X2            -4.3706
#> X3           -12.7878
#> omega_0_loc1   7.3123
#> omega_0_loc2   8.9323
#> omega_1       20.3501
#> omega_2       34.2158


# Stratified analysis when the first and second parameters ('omega_0' and 'omega_1')
# vary across two centers:
newLam1 <- inv.prior.cov(X1, lambda = c(rep(1, 3), 0.3, 0.7, 0.5, 0.8, 2),
                         family = 'survival', stratified = TRUE, basehaz = c("poly"),
                         strat_par = c(1, 2), L = 2)
# 'newLam1' is used as the prior for combined data:
list_newLam1 <- list(Lambda, newLam1)
bfi1 <- bfi(Lambda = list_newLam1, family = 'survival', theta_A_polys = theta_A_hats,
            stratified = TRUE, basehaz = c("poly"), p = 3, q_ls = qls,
            strat_par = c(1, 2))
summary(bfi1, cur_mat = TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘survival’ 
#>   Baseline: ‘poly’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   0.5626  0.1880  0.1942   0.9310
#> X2   0.6215  0.1513  0.3250   0.9181
#> X3   1.3566  0.2105  0.9441   1.7692
#> 
#> Minus the Curvature Matrix: 
#> 
#>                   X1      X2       X3 omega_0_loc1 omega_0_loc2 omega_1_loc1
#> X1           34.0208 -6.1851   0.2152       2.0008       0.5517       0.4276
#> X2           -6.1851 52.5297   5.7451       4.1509       1.2219       0.4604
#> X3            0.2152  5.7451  44.4738       3.4463       7.9651      -3.3272
#> omega_0_loc1  2.0008  4.1509   3.4463      18.8158       0.0000       9.5086
#> omega_0_loc2  0.5517  1.2219   7.9651       0.0000      18.7430       0.0000
#> omega_1_loc1  0.4276  0.4604  -3.3272       9.5086       0.0000       7.8123
#> omega_1_loc2 -3.0791 -3.6742  -2.9422       0.0000       8.5748       0.0000
#> omega_2      -6.3977 -4.3706 -12.7878       7.3123       8.9323       6.5539
#>              omega_1_loc2  omega_2
#> X1                -3.0791  -6.3977
#> X2                -3.6742  -4.3706
#> X3                -2.9422 -12.7878
#> omega_0_loc1       0.0000   7.3123
#> omega_0_loc2       8.5748   8.9323
#> omega_1_loc1       0.0000   6.5539
#> omega_1_loc2       9.7323  13.7961
#> omega_2           13.7961  34.2158


###---------------------###
###  Treatment Effect   ###
###---------------------###

set.seed(112358)

#-----------------------------#
# New Data for Local Center 1 #
#-----------------------------#
# Generating new data with 'treatment' variable
# We cansider the first variable (X1$X1) to be the treatment
X1$X1 <- sample(0:1, n1, replace=TRUE) # categorical variable
y1 <- surv.simulate(Z = list(X1), beta = theta[1:p], a = theta[5], b = theta[6],
                    u1 = 0.1, cen_rate = 0.3, gen_data_from = "weibul")$D[[1]][, 1:2]

#-----------------------------#
# New Data for Local Center 2 #
#-----------------------------#
# We cansider the first variable (X2$X1) to be the treatment!
X2$X1 <- sample(0:1, n2, replace=TRUE) # categorical variable
y2 <- surv.simulate(Z = list(X2), beta = theta[1:p], a = theta[5], b = theta[6],
                    u1 = 0.1, cen_rate = 0.3, gen_data_from = "weibul")$D[[1]][, 1:2]

#-------------#
# First Round #
#-------------#

## Center 1:
Lambda1 <- inv.prior.cov(X1, lambda = 0.01, family = 'survival',
                         treatment = "X1", treat_round="first")
# When treat_round = "first", the family will automatically set to 'binomial',
# even if family = 'gaussian' or family = 'survival'.
fit1_r1 <- MAP.estimation(y1, X1, family = 'survival', # 'basehaz' is not needed!
                          Lambda = Lambda1, treatment = "X1", treat_round = "first")
# While family = 'survival', the output is based on 'binomial' with no 'Intercept'!
# The output without the treatment (X1) in the first round!
summary(fit1_r1)
#> 
#> Summary of the local model:
#> 
#>    Formula: X1 ~ X2 + X3 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.1108  0.3732 -0.6206   0.8423
#> X2           -0.0297  0.3798 -0.7742   0.7148
#> X3           -0.4071  0.4199 -1.2301   0.4158
#> 
#> Dispersion parameter (sigma2):  1 
#>             log Lik Posterior:  -20.24 
#>                   Convergence:  0 

## Center 2:
Lambda2 <- inv.prior.cov(X2, lambda = 0.01, family = 'survival',
                         treatment = "X1", treat_round="first")
fit2_r1 <- MAP.estimation(y2, X2, family = 'survival', Lambda = Lambda2,
                          treatment = "X1", treat_round = "first")
fit2_r1
#> 
#> Local MAP estimates (family: binomial)
#> 
#> (Intercept)           X2           X3  
#>     0.27040     -0.03733      0.13470  
#> 

## Central Server:
theta_hats_r1 <- list(fit1_r1$theta_hat, fit2_r1$theta_hat)
A_hats_r1 <- list(fit1_r1$A_hat, fit2_r1$A_hat)
fitbfi_r1 <- bfi(theta_hats_r1, A_hats_r1, Lambda1, family = 'survival',
                 treat_round = "first")
# In the first round output is based on 'binomial', and without
# the intercept and treatment (X1):
summary(fitbfi_r1, cur_mat = TRUE)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘binomial’ 
#>       Link: ‘Logit’
#> 
#> Coefficients:
#> 
#>             Estimate Std.Dev CI 2.5% CI 97.5%
#> (Intercept)   0.1999  0.2627 -0.3149   0.7147
#> X2           -0.0210  0.2423 -0.4960   0.4540
#> X3           -0.0697  0.2585 -0.5763   0.4369
#> 
#> Dispersion parameter (sigma2):  1 
#> 
#> Minus the Curvature Matrix: 
#> 
#>             (Intercept)      X2      X3
#> (Intercept)     14.5585  0.8607 -0.5460
#> X2               0.8607 17.0885  0.3653
#> X3              -0.5460  0.3653 14.9987

#--------------#
# Second Round #
#--------------#

## Center 1:
Lambda11 <- inv.prior.cov(X1, lambda = 0.01, family = 'survival',
                          basehaz = "unspecified", treatment = "X1",
                          treat_round="second")
fit1_r2 <- MAP.estimation(y1, X1, family = 'survival', Lambda = Lambda11,
                          basehaz = "unspecified", treatment = "X1",
                          treat_round = "second", gamma_bfi = fitbfi_r1$theta_hat)
# The output with only the treatment (X1) in the second round!
summary(fit1_r2)
#> 
#> Summary of the local model:
#> 
#>    Formula: Survival(time, status) ~ X1 
#>     Family: ‘survival’ 
#>   Baseline: ‘unspecified’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   0.7284  0.3332  0.0753   1.3815
#> 
#> log Lik Posterior:  -121.7 
#>       Convergence:  0 

## Center 2:
Lambda22 <- inv.prior.cov(X2, lambda = 0.01, family = 'survival',
                          basehaz = "unspecified", treatment = "X1",
                          treat_round="second")
fit2_r2 <- MAP.estimation(y2, X2, family = 'survival', basehaz = "unspecified",
                          Lambda = Lambda22, treatment = "X1",
                          treat_round = "second", gamma_bfi = fitbfi_r1$theta_hat)
fit2_r2
#> 
#> Local MAP estimates (family: survival)
#> 
#>    X1  
#> 0.954  
#> 

## Centeral Server:
theta_hats_r2 <- list(fit1_r2$theta_hat, fit2_r2$theta_hat)
A_hats_r2 <- list(fit1_r2$A_hat, fit2_r2$A_hat)
fitbfi_r2 <- bfi(theta_hats_r2, A_hats_r2, Lambda11, family = 'survival',
                 basehaz = "unspecified", treat_round = "second")
# When family = 'survival', 'for_ATE' is not calculated.
summary(fitbfi_r2)
#> 
#> Summary of the BFI model:
#> 
#>     Family: ‘survival’ 
#>   Baseline: ‘unspecified’
#> 
#> Coefficients:
#> 
#>    Estimate Std.Dev CI 2.5% CI 97.5%
#> X1   0.8272  0.2504  0.3364    1.318

```
