# Generate survival data with predefined censoring rates for proportional hazards models

`surv.simulate` simulates one or multiple (right-censored) survival
datasets for proportional hazards models by simultaneously incorporating
a baseline hazard function from three different survival distributions
(exponential, Weibull and Gompertz), a random censoring time generated
from a uniform distribution with an known/unknown upper limit, and a set
of baseline covariates. When the upper limit of the uniform censoring
time distribution is unknown, `surv.simulate` can be used separately to
obtain the upper limit with a predefined censoring rate.

## Usage

``` r
surv.simulate(L = 1,
              Z,
              beta,
              a,
              b,
              u1 = 0,
              u2,
              cen_rate,
              gen_data_from = c("exp", "weibul", "gomp"),
              only_u2 = FALSE,
              n.rep = 100,
              Trace = FALSE)
```

## Arguments

- L:

  the number of datasets to be generated. Default is `L = 1`.

- Z:

  a list of `L` design matrices of dimension \\n\_\ell \times p\\, where
  \\n\_\ell\\ is the number of samples observed for the
  \\\ell^{\text{th}}\\ dataset and \\p\\ is the number of covariables.
  When `L = 1`, `Z` can be a matrix.

- beta:

  the vector of the (true) coefficients values, with a length of \\p\\
  (the number of covariates).

- a:

  scale parameter, which should be non-negative. See ‘Details’ for the
  form of the cumulative hazard that can be used.

- b:

  shape/location parameter, which should be non-negative. It is not used
  when `gen_data_from = `“exp”. See ‘Details’ for the form of the
  cumulative hazard that can be used.

- u1:

  a known non-negative lower limit of the uniform distribution for
  generating random censoring time. Default is `u1 = 0`. If `cen_rate`
  is not equal to 0, then `u1` does not need to be defined.

- u2:

  an non-negative upper limit of the uniform random censoring time
  distribution. The upper limit can be unknown (`u2 = NULL`, the
  default), or predefined. When this argument is assumed to be unknown,
  `u2 = NULL`, it is calculated by the algorithm within
  `surv.simulate()`. However, if the argument `u2` is known, the
  censoring rate cannot be predefined (meaning there is no control over
  it) and is calculated based on the generated dataset. See ‘Details’
  and ‘References’.

- cen_rate:

  a value representing the proportion of observations in the simulated
  survival data that are censored. The range of this argument is from 0
  to 1. When the upper limit is known, `cen_rate` can nor be predefined.
  If there is no censoring (`cen_rate = 0`), the lower (`u1`) and upper
  (`u2`) limits of the uniform distribution do not need to be specified.

- gen_data_from:

  a description of the distribution from which the time to event is
  generated. This is a character string and can be `exponential`
  (“exp”), `Weibull` (“weibul”), or `Gompertz` (“gomp”). Can be
  abbreviated. By default, the `exponential` distribution is used.

- only_u2:

  logical flag for calculating only the upper limit of the uniform
  censoring time distribution. If `only_u2 = TRUE`, the dataset(s) are
  not generated. If `only_u2 = TRUE`, the arguments `Z` and `u2` do not
  need to be specified, and `cen_rate` should not be set to \\0\\.
  Default is `only_u2 = FALSE`.

- n.rep:

  a scalar specifying the number of iterations. This argument is
  exclusively used in the case of the `Gompertz` distribution. Default
  is 100.

- Trace:

  logical flag indicating whether the output of the desired `u2` and the
  censoring proportion for different datasets should be produced for
  each iteration. It works `gen_data_from = `“gomp”.

## Details

`surv.simulate` function generates \\L\\ simulated right-censored
survival datasets from exponential, Weibull, or Gompertz distributions,
incorporating the covariates, `Z`, distributed according to a
multivariate `normal` distribution, with censoring time generated from a
uniform distribution `Uniform(u1, u2)`, where `u1` is known but `u2` can
be either known or unknown.

`surv.simulate()` can also be used to calculate the unknown upper limit
of the uniform distribution, `u2`, with a predefined censoring rate. To
do this, set `u2 = NULL` and `only_u2 = TRUE`. In this case, the
datasets are not generated; only `u2` is.

`surv.simulate()` uses a root-finding algorithm to select the censoring
parameter that achieves predefined censoring rates in the simulated
survival data.

When `gen_data_from = `“exp”:

- the cumulative baseline hazard function is considered as \\\Lambda_0=a
  t\\,

- the event time for the \\\ell^{\text{th}}\\ dataset, \\T\_\ell\\, is
  computed by \\ - log(u) \\ exp(- Z\_\ell \boldsymbol{\beta}) / a\\,
  where \\u\\ follows a standard uniform distribution;

For `gen_data_from = `“weibul”:

- the cumulative hazard function is as \\\Lambda_0=a t ^ b\\,

- the event time is computed by \\T\_\ell= (- log(u) \\ exp(- Z\_\ell
  \boldsymbol{\beta}) / a)^{1/b}\\, where \\u\\ follows a standard
  uniform distribution;

For `gen_data_from = `“gomp”:

- the cumulative hazard function is as \\\Lambda_0=a (exp(b t) - 1) /
  b\\,

- the event time is computed by \\T\_\ell= \log(1- log(u) \\ exp(-
  Z\_\ell \boldsymbol{\beta}) b / a) / b\\, where \\u\\ follows a
  standard uniform distribution;

Finally the survival time is obtained by \\\tilde{T}\_\ell=\min\\T\_\ell
, C\_\ell \\\\.

The function will be updated for `gen_data_from = `“gomp”.

## Value

`surv.simulate` returns a list containing the following components:

- D:

  a list of \\L\\ data frames, with dimension \\n\_\ell \times (p+2)\\.
  The first and second columns, named `time` and `status`, contain the
  simulated survival time and the censoring indicator, respectively,
  where \\0\\ means censored and \\1\\ means uncensored;

- censor_propor:

  the vector of censoring proportions in the simulated datasets `D`,
  containing \\L\\ elements;

- u1:

  the lower limit of the uniform distribution used to generate random
  censoring times with a predefined censoring rate. Sometimes this
  output is less than the value entered by the user, as it is adjusted
  to achieve the desired amount of censoring rate;

- u2:

  the upper limit of the uniform distribution used to generate random
  censoring times. If `u2 = NULL`, this output will be the estimated
  upper limit necessary to achieve the desired censoring rate across the
  \\L\\ datasets.

## References

Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker M.A.
(2026). *Bayesian federated inference for survival models*, *Journal of
Applied Statistics*, 53(2): 203-223.
\<https://doi.org/10.1080/02664763.2025.2511932\>

## Author

Hassan Pazira  
Maintainer: Hassan Pazira <h.pazira@arq.org>

## See also

[`MAP.estimation`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md)

## Examples

``` r

# Setting a seed for reproducibility
set.seed(1123)

#-------------------------
# Simulating Survival data
#-------------------------
N    <- c(7, 10, 13) # the sample sizes of 3 datasets
beta <- 1:4
p    <- length(beta)
L    <- 3

# Define a function to generate multivariate normal samples
mvrnorm_new <- function(n, mu, Sigma) {
    pp <- length(mu)
    e <- matrix(rnorm(n * pp), nrow = n)
    return(crossprod(t(e), chol(Sigma)) + matrix(mu, n, pp, byrow = TRUE))
}
Z <- list()
for (z in seq_len(L)) {
    Z[[z]] <- mvrnorm_new(n = N[z], mu = rep(0, p),
                          Sigma = diag(rep(1, p),p))
    colnames(Z[[z]]) <- paste0("Z_",seq_len(ncol(Z[[z]])))
}

# One simulated dataset from exponential distribution with no censoring:
surv_data <- surv.simulate(Z = Z[[1]], beta = beta, a = exp(-.9),
                           cen_rate = 0, gen_data_from = "exp")
surv_data
#> $D
#> $D[[1]]
#>           time status       Z.Z_1       Z.Z_2      Z.Z_3      Z.Z_4
#> 1  0.075564227      1  0.70473619  0.94551438 -0.6930954  0.6174432
#> 2  0.732609577      1  1.10710882  0.37164284 -0.1714535 -0.2798152
#> 3 15.148291300      1  0.94005308 -1.14369244 -1.1564422  0.9669025
#> 4  0.003073544      1 -0.91096743  0.31234925  1.9137865  0.3402063
#> 5  2.022124176      1  1.38294723 -0.17384622 -1.2173467  0.3753887
#> 6  0.055636648      1  0.19442801  0.05596198 -0.8097336  1.6147903
#> 7  0.328196224      1 -0.03972739  0.79823340 -0.3073386  0.6104483
#> 
#> 
#> $censor_propor
#> [1] 0
#> 
#> $u1
#> NULL
#> 
#> $u2
#> NULL
#> 
surv_data$D[[1]][,1:2] # The simulated survival data
#>           time status
#> 1  0.075564227      1
#> 2  0.732609577      1
#> 3 15.148291300      1
#> 4  0.003073544      1
#> 5  2.022124176      1
#> 6  0.055636648      1
#> 7  0.328196224      1

# Calculate only 'u2' with a predefined censoring rate of 0.4:
u2_new <- surv.simulate(Z = Z[1:2], beta = beta, a = exp(-.9),
                        b = exp(1.8), u1 = 0.1, only_u2 = TRUE,
                        cen_rate = 0.4, gen_data_from = "weibul")$u2
u2_new
#> [1] 3.252022

# Two simulated datasets with a known 'u2':
# Using 'u2_new' to help control over censoring rate (was chosen 0.4)
surv.simulate(Z = Z[1:2], beta = beta, a = exp(-.9), b = exp(1.8),
              u1 = 0.05, u2 = u2_new, gen_data_from = "weibul")
#> $D
#> $D[[1]]
#>        time status       Z.Z_1       Z.Z_2      Z.Z_3      Z.Z_4
#> 1 0.7399626      1  0.70473619  0.94551438 -0.6930954  0.6174432
#> 2 1.1207633      0  1.10710882  0.37164284 -0.1714535 -0.2798152
#> 3 1.7127411      1  0.94005308 -1.14369244 -1.1564422  0.9669025
#> 4 0.3660419      1 -0.91096743  0.31234925  1.9137865  0.3402063
#> 5 0.8653180      1  1.38294723 -0.17384622 -1.2173467  0.3753887
#> 6 0.3775353      1  0.19442801  0.05596198 -0.8097336  1.6147903
#> 7 0.7019206      1 -0.03972739  0.79823340 -0.3073386  0.6104483
#> 
#> $D[[2]]
#>         time status       Z.Z_1      Z.Z_2      Z.Z_3       Z.Z_4
#> 1  3.1926945      0 -0.23027404 -1.2690632 -0.2673166 -1.55025861
#> 2  0.3414662      1 -0.32916907  1.6265643  0.2278069  0.91480231
#> 3  1.1314558      1  0.49944094 -0.7696502  0.3879469 -0.02710971
#> 4  1.1982244      1 -0.92867906  0.2380289 -1.2894161  0.70696930
#> 5  1.4326297      1 -0.96832355  0.2955621 -0.4037798  0.42578390
#> 6  0.9277940      1 -0.05840795  0.9128333 -0.4899637  0.34036840
#> 7  0.3113212      0  1.07646211 -0.8527642  0.2666801  1.41542485
#> 8  0.1139205      0  1.68615000 -0.5295145  1.3293412  0.59569653
#> 9  0.6654312      1  1.93659046  1.7751477 -1.2878720  0.26822334
#> 10 1.0432595      0  0.51202240 -1.6117240  0.9368176 -0.75419444
#> 
#> 
#> $censor_propor
#> [1] 0.1428571 0.4000000
#> 
#> $u1
#> [1] 0.05
#> 
#> $u2
#> [1] 3.252022
#> 

# Three simulated datasets from 'weibul' with an unknown 'u2':
surv.simulate(Z = Z, beta = beta, a = exp(-1), b = exp(1),
               u1 = 0.01, cen_rate = 0.3, gen_data_from = "weibul")
#> $D
#> $D[[1]]
#>        time status       Z.Z_1       Z.Z_2      Z.Z_3      Z.Z_4
#> 1 0.2823065      1  0.70473619  0.94551438 -0.6930954  0.6174432
#> 2 0.9975437      1  1.10710882  0.37164284 -0.1714535 -0.2798152
#> 3 0.8144260      0  0.94005308 -1.14369244 -1.1564422  0.9669025
#> 4 0.1386992      1 -0.91096743  0.31234925  1.9137865  0.3402063
#> 5 1.6849989      1  1.38294723 -0.17384622 -1.2173467  0.3753887
#> 6 0.1502737      0  0.19442801  0.05596198 -0.8097336  1.6147903
#> 7 0.3864806      1 -0.03972739  0.79823340 -0.3073386  0.6104483
#> 
#> $D[[2]]
#>          time status       Z.Z_1      Z.Z_2      Z.Z_3       Z.Z_4
#> 1  7.75069014      0 -0.23027404 -1.2690632 -0.2673166 -1.55025861
#> 2  0.07093456      1 -0.32916907  1.6265643  0.2278069  0.91480231
#> 3  1.34441157      1  0.49944094 -0.7696502  0.3879469 -0.02710971
#> 4  1.85809906      1 -0.92867906  0.2380289 -1.2894161  0.70696930
#> 5  0.88334276      1 -0.96832355  0.2955621 -0.4037798  0.42578390
#> 6  0.72544536      0 -0.05840795  0.9128333 -0.4899637  0.34036840
#> 7  0.10815548      1  1.07646211 -0.8527642  0.2666801  1.41542485
#> 8  0.10718719      1  1.68615000 -0.5295145  1.3293412  0.59569653
#> 9  0.33155256      1  1.93659046  1.7751477 -1.2878720  0.26822334
#> 10 1.31690683      1  0.51202240 -1.6117240  0.9368176 -0.75419444
#> 
#> $D[[3]]
#>          time status       Z.Z_1       Z.Z_2       Z.Z_3      Z.Z_4
#> 1  1.54786085      1 -1.44153951  0.84543041  1.63435118 -1.3325912
#> 2  0.28011636      1  0.42547690  0.91067540 -0.02140968  0.3038766
#> 3  0.07408815      1  0.63158367  0.08532232  0.64390401  1.2853559
#> 4  0.98235600      0 -0.15201586 -0.91096893  0.32015202 -0.1719336
#> 5  0.09308440      0 -1.48957048  0.09515951  1.29546216 -0.1569018
#> 6  0.20262018      1  0.50347542  0.65439581 -0.57036674  1.1085850
#> 7  0.78550533      0  0.95250295 -0.71764421  0.51880688 -1.5608051
#> 8  0.04926487      1  0.04937487  1.07616067  1.37303835  0.3469306
#> 9  0.38080909      1  0.37153942 -0.04582918 -0.45533515  1.2241016
#> 10 3.21427845      0 -0.90718225 -2.49706485  1.40629078 -0.8608085
#> 11 0.11008080      1  2.53099268 -0.03519585  0.12720354  1.4766462
#> 12 0.15020603      1  1.13004355  1.34211808  0.83315557 -0.4462548
#> 13 6.77948385      1 -0.73382185  0.15624687 -1.35326133  0.2320178
#> 
#> 
#> $censor_propor
#> [1] 0.2857143 0.2000000 0.3076923
#> 
#> $u1
#> [1] 0.01
#> 
#> $u2
#> [1] 9.967712
#> 

# Two simulated datasets from 'gomp' with unknown 'u2' and censoring rate of 0.3:
surv.simulate(Z = Z[2:3], beta = beta, a = exp(1), b = exp(2), u1 = 0.1,
              cen_rate = 0.3, gen_data_from = "gomp", Trace = TRUE)
#> censoring proportion for silos are:  0.01 0.01   with u2 = 20 
#> censoring proportion for silos are:  0.01 0.01   with u2 = 14.19 
#> censoring proportion for silos are:  0.02 0.02   with u2 = 10.10383 
#> censoring proportion for silos are:  0.03 0.02   with u2 = 7.28991 
#> censoring proportion for silos are:  0.04 0.03   with u2 = 5.275091 
#> censoring proportion for silos are:  0.05 0.04   with u2 = 3.895249 
#> censoring proportion for silos are:  0.07 0.07   with u2 = 2.896867 
#> censoring proportion for silos are:  0.09 0.06   with u2 = 2.225017 
#> censoring proportion for silos are:  0.12 0.09   with u2 = 1.72242 
#> censoring proportion for silos are:  0.15 0.12   with u2 = 1.39092 
#> censoring proportion for silos are:  0.17 0.15   with u2 = 1.162114 
#> censoring proportion for silos are:  0.17 0.18   with u2 = 0.9992837 
#> censoring proportion for silos are:  0.18 0.2   with u2 = 0.8723362 
#> censoring proportion for silos are:  0.2 0.21   with u2 = 0.7788285 
#> censoring proportion for silos are:  0.21 0.24   with u2 = 0.7064274 
#> censoring proportion for silos are:  0.22 0.24   with u2 = 0.6536355 
#> censoring proportion for silos are:  0.23 0.25   with u2 = 0.6087609 
#> censoring proportion for silos are:  0.24 0.25   with u2 = 0.5731016 
#> censoring proportion for silos are:  0.26 0.26   with u2 = 0.5429697 
#> censoring proportion for silos are:  0.26 0.27   with u2 = 0.5212091 
#> censoring proportion for silos are:  0.26 0.28   with u2 = 0.5016437 
#> censoring proportion for silos are:  0.26 0.27   with u2 = 0.4858998 
#> censoring proportion for silos are:  0.27 0.28   with u2 = 0.4684635 
#> censoring proportion for silos are:  0.29 0.29   with u2 = 0.4560853 
#> censoring proportion for silos are:  0.28 0.28   with u2 = 0.4499632 
#> censoring proportion for silos are:  0.27 0.28   with u2 = 0.4408082 
#> censoring proportion for silos are:  0.28 0.3   with u2 = 0.4296693 
#> 
#>  Desired u2 is: 0.4296693 ( with u1 = 0.1 ) 
#>  
#> $D
#> $D[[1]]
#>            time status       Z.Z_1      Z.Z_2      Z.Z_3       Z.Z_4
#> 1  3.491216e-01      0 -0.23027404 -1.2690632 -0.2673166 -1.55025861
#> 2  8.252691e-04      1 -0.32916907  1.6265643  0.2278069  0.91480231
#> 3  8.251961e-06      1  0.49944094 -0.7696502  0.3879469 -0.02710971
#> 4  3.457653e-01      1 -0.92867906  0.2380289 -1.2894161  0.70696930
#> 5  1.383742e-01      1 -0.96832355  0.2955621 -0.4037798  0.42578390
#> 6  3.500657e-02      1 -0.05840795  0.9128333 -0.4899637  0.34036840
#> 7  3.473370e-04      1  1.07646211 -0.8527642  0.2666801  1.41542485
#> 8  3.788401e-04      1  1.68615000 -0.5295145  1.3293412  0.59569653
#> 9  1.746997e-02      1  1.93659046  1.7751477 -1.2878720  0.26822334
#> 10 1.887190e-01      0  0.51202240 -1.6117240  0.9368176 -0.75419444
#> 
#> $D[[2]]
#>            time status       Z.Z_1       Z.Z_2       Z.Z_3      Z.Z_4
#> 1  1.908368e-01      1 -1.44153951  0.84543041  1.63435118 -1.3325912
#> 2  5.553927e-03      1  0.42547690  0.91067540 -0.02140968  0.3038766
#> 3  2.554467e-05      1  0.63158367  0.08532232  0.64390401  1.2853559
#> 4  2.509140e-01      1 -0.15201586 -0.91096893  0.32015202 -0.1719336
#> 5  7.807589e-02      1 -1.48957048  0.09515951  1.29546216 -0.1569018
#> 6  8.807375e-03      1  0.50347542  0.65439581 -0.57036674  1.1085850
#> 7  2.237849e-01      0  0.95250295 -0.71764421  0.51880688 -1.5608051
#> 8  2.597067e-04      1  0.04937487  1.07616067  1.37303835  0.3469306
#> 9  7.799071e-03      1  0.37153942 -0.04582918 -0.45533515  1.2241016
#> 10 2.913660e-01      0 -0.90718225 -2.49706485  1.40629078 -0.8608085
#> 11 2.111413e-05      1  2.53099268 -0.03519585  0.12720354  1.4766462
#> 12 5.090627e-03      1  1.13004355  1.34211808  0.83315557 -0.4462548
#> 13 4.275830e-01      0 -0.73382185  0.15624687 -1.35326133  0.2320178
#> 
#> 
#> $censor_propor
#> [1] 0.2000000 0.2307692
#> 
#> $u1
#> [1] 0.1
#> 
#> $u2
#> [1] 0.4296693
#> 
```
