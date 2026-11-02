# An Introduction to BFI in R

## Overview

The R package `BFI` (**B**ayesian **F**ederated **I**nference) provides
several functions to perform Bayesian Federated Inference for two types
of models (`GLM` and `Survival`) using multicenter data without the need
to combine or share them. This tutorial focuses on `GLM` models. Two
commonly used families for `GLM` models, `"binomial"` and `"gaussian"`,
are available for this version of the package. The most commonly used
functions include
[`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md),
[`MAP.estimation()`](https://hassanpazira.github.io/BFI/reference/MAP.estimation.md),
and
[`inv.prior.cov()`](https://hassanpazira.github.io/BFI/reference/inv.prior.cov.md).
In the following, we will see how the `BFI` package can be applied to
real data sets included in the package.

## How to use it?

Before we go on, we first install and load the `BFI` package:

``` r

# Install and load the BFI package from CRAN:
# install.packages("BFI")
library(BFI)
```

By using the following code, we can see that there are two available
data sets in the package: `trauma` and `Nurses`.

``` r

data(package = "BFI")
```

The `trauma` data set can be utilized for the `"binomial"` family and
`Nurses` data set can be used for `"gaussian"` family. To avoid
repetition, we will only use the `trauma` data set. By loading the
package, the data sets included will be loaded and can be inspected as
follows:

``` r

# Get the number of rows and columns
dim(trauma)
```

    ## [1] 371   6

``` r

# To get an idea of the data set, print the first 7 rows
head(trauma, 7)
```

    ##   sex age hospital ISS GCS mortality
    ## 1   1  20        3  24  15         0
    ## 2   0  38        3  34  13         0
    ## 3   0  37        3  50  15         0
    ## 4   0  17        3  43   4         1
    ## 5   0  49        3  29  15         0
    ## 6   0  30        3  22  15         0
    ## 7   1  84        2  66   3         1

This data set consists of data of 371 trauma patients from three
hospitals (peripheral hospital without a neuro-surgical unit,
`'status=1'`, peripheral hospital with a neuro-surgical unit,
`status=2`, and academic medical center, `status=3`).

As we can see, the data set has six columns. The covariates `sex`
(dichotomous), `age` (continuous), `ISS` (Injury Severity Score,
continuous), and `GCS` (Glasgow Coma Scale, continuous), which serve as
the predictors. `mortality` is the response variable, while `hospital`
is a categorical variable which indicates the hospitals involved in the
study. For more information about this data set use

``` r

# Get some info about the data set from the help file
?trauma
```

We will analyze the data with a `logistic` regression model. We apply
fixed, prespecified scaling to the continuous predictors. This avoids
calculating scaling parameters from the combined data, which would not
generally be available in a real federated setting. The transformations
correspond to effects per 10 years of age, 10 ISS points, and 5 GCS
points.

``` r

# Fixed, prespecified scaling.
# No pooled means or standard deviations are required.
trauma$age <- trauma$age / 10
trauma$ISS <- trauma$ISS / 10
trauma$GCS <- trauma$GCS / 5
trauma$hospital <- factor(trauma$hospital)
```

By using the following code we can see there are three hospitals
involved in the study:

``` r

length(levels(trauma$hospital))
```

    ## [1] 3

### MAP estimations at the local centers

Each participating center fits the agreed model using its own data. In
this example, three centers are represented by separate subsets of the
trauma dataset. In an actual federated analysis, each center would
execute the local code at its own institution.

The model family, predictor definitions, parameter names, and prior
specification must be agreed upon before local estimation.

``` r

## Center 1:
data1 <- subset(trauma, hospital == 1)
X1 <- data1[, c("sex", "age", "ISS", "GCS")]
Lambda1 <- inv.prior.cov(X1, lambda=0.01, L=3, family="binomial")
fit1 <- MAP.estimation(formula = mortality ~ sex + age + ISS + GCS,
                       data = data1, family="binomial", Lambda=Lambda1, 
                       control = list(maxit = 1000))
summary(fit1)
```

    ## 
    ## Summary of the local model:
    ## 
    ##    Formula: mortality ~ sex + age + ISS + GCS 
    ##     Family: 'binomial' 
    ##       Link: 'Logit'
    ## 
    ## Coefficients:
    ## 
    ##             Estimate Std.Dev  CI 2.5% CI 97.5%
    ## (Intercept)  -5.6186  2.7222 -10.9540  -0.2831
    ## sex          -5.1018  2.5423 -10.0847  -0.1190
    ## age           0.8850  0.3343   0.2298   1.5401
    ## ISS           1.5738  0.6047   0.3887   2.7590
    ## GCS          -2.4937  0.8329  -4.1262  -0.8613
    ## 
    ## Dispersion parameter (sigma2):  1 
    ##             log Lik Posterior:  -10.86 
    ##                   Convergence:  0

``` r

stopifnot(fit1$convergence == 0)
```

The output `fit1` from the local centers should be sent to the central
server for further analysis. To send the list from R to the central
server (which also uses R), you can save them in a format that R can
easily read, such as an RDS file.

``` r

# Save only the required statistical summaries
saveRDS(fit1, file="center1.rds")
```

``` r

## Center 2:
data2 <- subset(trauma, hospital == 2)
X2 <- data2[, c("sex", "age", "ISS", "GCS")]
Lambda2 <- inv.prior.cov(X2, lambda=0.01, L=3, family="binomial")
fit2 <- MAP.estimation(formula = mortality ~ sex + age + ISS + GCS,
                       data = data2, family="binomial", Lambda=Lambda2, 
                       control = list(maxit = 1000))
summary(fit2)
```

    ## 
    ## Summary of the local model:
    ## 
    ##    Formula: mortality ~ sex + age + ISS + GCS 
    ##     Family: 'binomial' 
    ##       Link: 'Logit'
    ## 
    ## Coefficients:
    ## 
    ##             Estimate Std.Dev CI 2.5% CI 97.5%
    ## (Intercept)  -1.4468  1.3298 -4.0532   1.1596
    ## sex          -0.7884  0.6999 -2.1601   0.5833
    ## age           0.9494  0.2259  0.5066   1.3922
    ## ISS           0.3574  0.2527 -0.1379   0.8527
    ## GCS          -1.9243  0.4221 -2.7516  -1.0969
    ## 
    ## Dispersion parameter (sigma2):  1 
    ##             log Lik Posterior:  -36.86 
    ##                   Convergence:  0

``` r

saveRDS(fit2, file = "center2.rds")

## Center 3:
data3 <- subset(trauma, hospital == 3)
X3 <- data3[, c("sex", "age", "ISS", "GCS")]
Lambda3 <- inv.prior.cov(X3, lambda=0.01, L=3, family="binomial")
fit3 <- MAP.estimation(formula = mortality ~ sex + age + ISS + GCS,
                       data = data3, family="binomial", Lambda=Lambda3, 
                       control = list(maxit = 1000))
summary(fit3)
```

    ## 
    ## Summary of the local model:
    ## 
    ##    Formula: mortality ~ sex + age + ISS + GCS 
    ##     Family: 'binomial' 
    ##       Link: 'Logit'
    ## 
    ## Coefficients:
    ## 
    ##             Estimate Std.Dev CI 2.5% CI 97.5%
    ## (Intercept)  -0.1529  1.0466 -2.2043   1.8985
    ## sex           0.1600  0.5718 -0.9607   1.2807
    ## age           0.6317  0.1433  0.3509   0.9126
    ## ISS           0.2053  0.1822 -0.1518   0.5623
    ## GCS          -2.4847  0.4235 -3.3147  -1.6548
    ## 
    ## Dispersion parameter (sigma2):  1 
    ##             log Lik Posterior:  -50.6 
    ##                   Convergence:  0

``` r

saveRDS(fit3, file = "center3.rds")
```

To see more information about the data set, such as the number of
observations and parameters, we can use the output of the
`MAP.estimation` function as follows:

``` r

# number of samples in center 1
fit1$n
```

    ## [1] 49

``` r

# number of parameters in center 1
fit1$np
```

    ## [1] 5

``` r

# number of samples in center 2
fit2$n
```

    ## [1] 106

``` r

# number of samples in center 3
fit3$n
```

    ## [1] 216

Additionally, before conducting the analysis, we can use the `n.par`
function to retrieve this information.

### BFI estimations

Now, the received RDS files can be loaded in R using the following
lines. They are assumed to be stored in the current working directory,
so they can be loaded using relative file paths.

``` r

# Central server: read the received files
center1 <- readRDS("center1.rds")
center2 <- readRDS("center2.rds")
center3 <- readRDS("center3.rds")
```

On the central server, the
[`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md) function
can be used to obtain the BFI estimations:

``` r

theta_hats <- list(center1$theta_hat, center2$theta_hat, center3$theta_hat)
A_hats     <- list(center1$A_hat, center2$A_hat, center3$A_hat)
Lambda_com <- center1$Lambda
Lambdas    <- list(center1$Lambda, center2$Lambda, center3$Lambda, Lambda_com)
BFI_fits   <- bfi(theta_hats, A_hats, Lambdas, family="binomial")
summary(BFI_fits, cur_mat=TRUE)
```

    ## 
    ## Summary of the BFI model:
    ## 
    ##     Family: 'binomial' 
    ##       Link: 'Logit'
    ## 
    ## Coefficients:
    ## 
    ##             Estimate Std.Dev CI 2.5% CI 97.5%
    ## (Intercept)  -1.0767  0.7422 -2.5313   0.3779
    ## sex          -0.2485  0.4185 -1.0687   0.5717
    ## age           0.5915  0.1059  0.3840   0.7991
    ## ISS           0.3430  0.1344  0.0797   0.6064
    ## GCS          -1.7831  0.2556 -2.2840  -1.2822
    ## 
    ## Dispersion parameter (sigma2):  1 
    ## 
    ## Minus the Curvature Matrix: 
    ## 
    ##             (Intercept)     sex      age      ISS      GCS
    ## (Intercept)     30.4742  7.8008 109.8935 120.0781  46.6773
    ## sex              7.8008  7.8108  30.4465  30.3379  11.5360
    ## age            109.8935 30.4465 531.2711 418.2228 201.0404
    ## ISS            120.0781 30.3379 418.2228 532.7746 173.8729
    ## GCS             46.6773 11.5360 201.0404 173.8729  95.6265

### Comparison

To compare the performance of the BFI methodology, we can combine the
data sets and obtain the MAP estimations based on the combined data:

``` r

# MAP estimates of the combined data:
X_combined  <- data.frame(sex=trauma$sex,
                          age=trauma$age,
                          ISS=trauma$ISS,
                          GCS=trauma$GCS)
Lambda   <- inv.prior.cov(X=X_combined, lambda=0.01, L=3, family="binomial")
fit_comb  <- MAP.estimation(y=trauma$mortality, X=X_combined, family="binomial", Lambda=Lambda) 
summary(fit_comb, cur_mat=TRUE)
```

    ## 
    ## Summary of the local model:
    ## 
    ##    Formula: y ~ sex + age + ISS + GCS 
    ##     Family: 'binomial' 
    ##       Link: 'Logit'
    ## 
    ## Coefficients:
    ## 
    ##             Estimate Std.Dev CI 2.5% CI 97.5%
    ## (Intercept)  -1.0955  0.7088 -2.4846   0.2936
    ## sex          -0.3357  0.3990 -1.1178   0.4464
    ## age           0.6635  0.0996  0.4684   0.8586
    ## ISS           0.3801  0.1288  0.1277   0.6326
    ## GCS          -2.0485  0.2439 -2.5265  -1.5705
    ## 
    ## Dispersion parameter (sigma2):  1 
    ##             log Lik Posterior:  -109.3 
    ##                   Convergence:  0 
    ## 
    ## Minus the Curvature Matrix: 
    ## 
    ##             (Intercept)     sex      age      ISS      GCS
    ## (Intercept)     33.5583  8.4925 124.1472 130.9550  51.5212
    ## sex              8.4925  8.5025  34.0277  32.1266  13.1350
    ## age            124.1472 34.0277 617.1691 464.4325 229.6154
    ## ISS            130.9550 32.1266 464.4325 576.5327 189.4027
    ## GCS             51.5212 13.1350 229.6154 189.4027 106.3280

Now, we can see the difference between the BFI and combined estimates:

``` r

# Squared Errors:
(fit_comb$theta_hat - BFI_fits$theta_hat)^2
```

    ##       (Intercept)         sex         age         ISS        GCS
    ## [1,] 0.0003533092 0.007600807 0.005180229 0.001375924 0.07042184

The differences illustrate the agreement between BFI and pooled-data MAP
estimates in this example. Exact equality is not guaranteed because BFI
relies on a quadratic approximation.

For more details see the following references.

## References

Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker
M.A. (2026). , , 53(2): 203-223.
<https://doi.org/10.1080/02664763.2025.2511932>

Jonker M.A., Pazira H. and Coolen A.C.C. (2025). , , 16(2): 383-423.
<https://doi.org/10.1017/rsm.2025.6>

Jonker M.A., Pazira H. and Coolen A.C.C. (2024). , , 43(12): 2421-2438.
<https://doi.org/10.1002/sim.10072>

## Contact

If you find any errors, have any suggestions, or would like to request
that something be added, please file an issue at [issue
report](https://github.com/hassanpazira/BFI/issues/) or send an email
to: <h.pazira@arq.org>.
