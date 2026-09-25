# Compute the estimated (baseline/cumulative) hazard and (baseline) survival functions

For a given vector of times, `hazards.fun` computes the estimated
baseline hazard, cumulative baseline hazard, hazard, baseline survival,
and survival functions. It can be used for prediction on a new sample.

## Usage

``` r
hazards.fun(time,
            z = NULL,
            p,
            theta_hat,
            basehaz = c("weibul", "exp", "gomp", "poly", "pwexp"),
            q_max,
            timax)
```

## Arguments

- time:

  the vector containing the time values for which the hazard rate is
  computed. If the argument `z` is not `NULL`, then the length of the
  argument `time` should be the number of columns of `z`, which is
  \\p\\.

- z:

  a new observation vector of length \\p\\. If `z = NULL` (the default),
  then the relative risk (\\\boldsymbol{z}^{\top} \boldsymbol{\beta}\\)
  is considered a vector of 1 with length \\n\\.

- p:

  the number of coefficients. It is taken equal to the number of
  elements of the argument `z`, if `z` is not `NULL`.

- theta_hat:

  a vector contains the values of the estimated parameters. The first
  \\p\\ values represent the coefficient parameters
  (\\\boldsymbol{\beta}\\), while the remaining values pertain to the
  parameters of the baseline hazard function (\\\boldsymbol{\omega}\\).

- basehaz:

  a character string representing one of the available baseline hazard
  functions; `exponential` (`"exp"`), `Weibull` (`"weibul"`, the
  default), `Gompertz` (`"gomp"`), `exponentiated polynomial`
  (`"poly"`), and `piecewise exponential` (`"pwexp"`). Can be
  abbreviated.

- q_max:

  a value represents the order of the exponentiated polynomial baseline
  hazard function. This argument should only be used when
  `basehaz = "poly"`. In the case of multiple centers, the maximum value
  of the orders should be used.
  [`ql.LRT()`](https://hassanpazira.github.io/BFI/reference/BFI-internal.md)
  can be used to obtain the order of each center.

- timax:

  a value represents the minimum (or maximum) value of the maximum times
  observed in the different centers. This argument should only be used
  when `basehaz = "pwexp"`.

## Details

`hazards.fun` computes the estimated baseline hazard, cumulative
baseline hazard, hazard, baseline survival, and survival functions at
different time points specified in the argument `time`.

The function `hazards.fun()` can be used for prediction purposes with
new sample. The arguments `time` and `z` should be provided for the new
data.

## Value

`hazards.fun` returns a list containing the following components:

- bhazard:

  the vector of estimates of the baseline hazard function at the time
  points given by the argument `time`;

- cbhazard:

  the vector of estimates of the cumulative baseline hazard function at
  the time points specified in the argument `time`;

- bsurvival:

  the vector of estimates of the baseline survival function at the time
  points given by the argument `time`;

- hazard:

  the vector of estimates of the hazard function at the time points
  given by the argument `time`;

- chazard:

  the vector of estimates of the cumulative hazard function at the time
  points specified in the argument `time`;

- survival:

  the vector of estimates of the survival function at the time points
  given by the argument `time`.

## References

Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker M.A.
(2026). *Bayesian federated inference for survival models*, *Journal of
Applied Statistics*, 53(2): 203-223.
\<https://doi.org/10.1080/02664763.2025.2511932\>

## Author

Hassan Pazira  

## See also

[`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)

## Examples

``` r
# Setting a seed for reproducibility
set.seed(1123)

##-------------------------
## Simulating Survival data
##-------------------------

n <- 40
p <- 7
Original_data <- data.frame(matrix(rnorm((n+1) * p), (n+1), p))
X <- Original_data[1:n,]
X_new <- Original_data[(n+1),]
# Simulating survival data from Exponential distribution
# with a predefined censoring rate of 0.2:
Orig_y <- surv.simulate(Z = Original_data, beta = rep(1,p), a = exp(1),
                        cen_rate = 0.2, gen_data_from = "exp")$D[[1]][,1:2]
y <- Orig_y[1:n,]
y_new <- Orig_y[(n+1),]
time_points <- seq(0, max(y$time), length.out=20)

#------------------------
# Weibull baseline hazard
#------------------------

Lambda <- inv.prior.cov(X, lambda = c(0.5, 1), family = 'survival', basehaz = 'weibul')
fit_weib <- MAP.estimation(y, X, family = 'survival', Lambda = Lambda,
                           basehaz = "weibul")
# reltive risk is 1:
hazards.fun(time = time_points, p = p, theta_hat = fit_weib$theta_hat,
            basehaz = "weibul")
#> $bhazard
#>  [1]      Inf 2.978474 2.838078 2.759040 2.704300 2.662590 2.628988 2.600909
#>  [9] 2.576829 2.555773 2.537084 2.520296 2.505066 2.491138 2.478311 2.466429
#> [17] 2.455365 2.445018 2.435302 2.426148
#> 
#> $cbhazard
#>  [1]  0.000000  1.139229  2.171059  3.165896  4.137445  5.092038  6.033332
#>  [8]  6.963708  7.884839  8.797963  9.704032 10.603800 11.497880 12.386779
#> [15] 13.270923 14.150674 15.026347 15.898213 16.766511 17.631454
#> 
#> $bsurvival
#>  [1] 1.000000e+00 3.200656e-01 1.140567e-01 4.217634e-02 1.596358e-02
#>  [6] 6.145483e-03 2.397493e-03 9.455839e-04 3.764071e-04 1.510404e-04
#> [11] 6.103691e-05 2.482151e-05 1.015159e-05 4.173402e-06 1.723898e-06
#> [16] 7.152208e-07 2.979480e-07 1.245931e-07 5.228735e-08 2.201693e-08
#> 
#> $hazard
#>  [1]      Inf 2.978474 2.838078 2.759040 2.704300 2.662590 2.628988 2.600909
#>  [9] 2.576829 2.555773 2.537084 2.520296 2.505066 2.491138 2.478311 2.466429
#> [17] 2.455365 2.445018 2.435302 2.426148
#> 
#> $chazard
#>  [1]  0.000000  1.139229  2.171059  3.165896  4.137445  5.092038  6.033332
#>  [8]  6.963708  7.884839  8.797963  9.704032 10.603800 11.497880 12.386779
#> [15] 13.270923 14.150674 15.026347 15.898213 16.766511 17.631454
#> 
#> $survival
#>  [1] 1.000000e+00 3.200656e-01 1.140567e-01 4.217634e-02 1.596358e-02
#>  [6] 6.145483e-03 2.397493e-03 9.455839e-04 3.764071e-04 1.510404e-04
#> [11] 6.103691e-05 2.482151e-05 1.015159e-05 4.173402e-06 1.723898e-06
#> [16] 7.152208e-07 2.979480e-07 1.245931e-07 5.228735e-08 2.201693e-08
#> 

#-------------------------
# Gompertz baseline hazard
#-------------------------
fit_gomp <- MAP.estimation(y, X, family = 'survival', Lambda = Lambda,
                           basehaz = "gomp")
# different time points
hazards.fun(time=1:max(y*2), p = p, theta_hat = fit_gomp$theta_hat,
            basehaz = "gomp")
#> $bhazard
#>  [1]  3.689346  4.847122  6.368226  8.366677 10.992274 14.441826 18.973902
#>  [8] 24.928216 32.751088 43.028904 56.532062 74.272728 97.580698
#> 
#> $cbhazard
#>  [1]   3.228712   7.470645  13.043765  20.365818  29.985650  42.624338
#>  [7]  59.229249  81.045046 109.706990 147.363511 196.837250 261.836640
#> [13] 347.233879
#> 
#> $bsurvival
#>  [1]  3.960849e-02  5.695609e-04  2.163541e-06  1.429676e-09  9.492877e-14
#>  [6]  3.079536e-19  1.892623e-26  6.347220e-36  2.263911e-48  1.001937e-64
#> [11]  3.270909e-86 1.931043e-114 1.578493e-151
#> 
#> $hazard
#>  [1]  3.689346  4.847122  6.368226  8.366677 10.992274 14.441826 18.973902
#>  [8] 24.928216 32.751088 43.028904 56.532062 74.272728 97.580698
#> 
#> $chazard
#>  [1]   3.228712   7.470645  13.043765  20.365818  29.985650  42.624338
#>  [7]  59.229249  81.045046 109.706990 147.363511 196.837250 261.836640
#> [13] 347.233879
#> 
#> $survival
#>  [1]  3.960849e-02  5.695609e-04  2.163541e-06  1.429676e-09  9.492877e-14
#>  [6]  3.079536e-19  1.892623e-26  6.347220e-36  2.263911e-48  1.001937e-64
#> [11]  3.270909e-86 1.931043e-114 1.578493e-151
#> 


##----------------------------
## Prediction for a new sample
##----------------------------

## Exponentiated polynomial (poly) baseline hazard:
Lambda <- inv.prior.cov(X, lambda = c(0.5, 1), family = 'survival', basehaz = "poly")
fit_poly <- MAP.estimation(y, X, family = 'survival', Lambda = Lambda,
                           basehaz = "poly")
hazards.fun(time = y_new$time, z = X_new, theta_hat = fit_poly$theta_hat,
            basehaz = "poly", q_max = fit_poly$q_l)
#> $bhazard
#> [1] 3.18832
#> 
#> $cbhazard
#> [1] 0.3483998
#> 
#> $bsurvival
#> [1] 0.7058166
#> 
#> $hazard
#> [1] 67.42799
#> 
#> $chazard
#> [1] 7.368112
#> 
#> $survival
#> [1] 0.0006310584
#> 

## Piecewise Exponential (pwexp) baseline hazard:
Lambda <- inv.prior.cov(X, lambda = c(0.5, 1), family = 'survival', basehaz = "pwexp")
fit_pw <- MAP.estimation(y, X, family='survival', Lambda=Lambda, basehaz="pwexp",
                         min_max_times = max(y))
#> 
#>  No. observations in the intervals :  36 1 0 3 
#>  
hazards.fun(time = y_new$time, z = X_new, theta_hat = fit_pw$theta_hat,
            basehaz = "pwexp", timax = max(y))
#> $bhazard
#> [1] 3.305502
#> 
#> $cbhazard
#> [1] 0.3612048
#> 
#> $bsurvival
#> [1] 0.6968363
#> 
#> $hazard
#> [1] 63.9545
#> 
#> $chazard
#> [1] 6.988551
#> 
#> $survival
#> [1] 0.000922382
#> 
```
