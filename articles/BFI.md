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
install.packages("BFI")
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

We will analyze the data with a `logistic` regression model. First we
standardize the covariates. This is not necessary for the analysis, but
is done for the interpretability of the accuracy of the estimates.

``` r

trauma$age <- scale(trauma$age)
trauma$ISS <- scale(trauma$ISS) 
trauma$GCS <- scale(trauma$GCS) 
trauma$hospital <- as.factor(trauma$hospital)
```

By using the following code we can see there are three hospitals
involved in the study:

``` r

length(levels(trauma$hospital))
```

    ## [1] 3

### MAP estimations

Therefore, the `MAP.estimation` function should be applied to these 3
local data sets separately to obtain the MAP estimations. Note that, in
practice, we do not have access to the combined data, and each center
should perform the analysis independently and send the output to the
central server, as follows:

``` r

# Center 1:
# X1 <- data.frame(sex=trauma$sex[trauma$hospital==1],
#                  age=trauma$age[trauma$hospital==1],
#                  ISS=trauma$ISS[trauma$hospital==1],
#                  GCS=trauma$GCS[trauma$hospital==1])
X1 <- subset(trauma, hospital == 1, select = c(sex, age, ISS, GCS))
Lambda1 <- inv.prior.cov(X1, lambda=0.01, L=3, family="binomial")
fit1 <- MAP.estimation(y=trauma$mortality[trauma$hospital==1], X=X1, family="binomial", Lambda=Lambda1)
summary(fit1)
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
    ##             Estimate Std.Dev  CI 2.5% CI 97.5%
    ## (Intercept)  -2.1226  0.8074  -3.7050  -0.5402
    ## sex          -5.2795  2.6305 -10.4351  -0.1239
    ## age           1.8864  0.7210   0.4732   3.2996
    ## ISS           2.3741  0.9250   0.5611   4.1872
    ## GCS          -2.4522  0.8295  -4.0781  -0.8264
    ## 
    ## Dispersion parameter (sigma2):  1 
    ##             log Lik Posterior:  -10.74 
    ##                   Convergence:  0

``` r

# Center 2:
# X2 <- data.frame(sex=trauma$sex[trauma$hospital==2],
#                  age=trauma$age[trauma$hospital==2],
#                  ISS=trauma$ISS[trauma$hospital==2],
#                  GCS=trauma$GCS[trauma$hospital==2])
X2 <- subset(trauma, hospital == 2, select = c(sex, age, ISS, GCS))
Lambda2 <- inv.prior.cov(X2, lambda=0.01, L=3, family="binomial")
fit2 <- MAP.estimation(y=trauma$mortality[trauma$hospital==2], X=X2, family="binomial", Lambda=Lambda2)
summary(fit2)
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
    ## (Intercept)  -0.7854  0.3789 -1.5280  -0.0427
    ## sex          -0.7850  0.6999 -2.1568   0.5868
    ## age           1.9601  0.4662  1.0465   2.8738
    ## ISS           0.5216  0.3673 -0.1983   1.2415
    ## GCS          -1.8737  0.4115 -2.6803  -1.0672
    ## 
    ## Dispersion parameter (sigma2):  1 
    ##             log Lik Posterior:  -36.87 
    ##                   Convergence:  0

``` r

# Center 3:
# X3 <- data.frame(sex=trauma$sex[trauma$hospital==3],
#                  age=trauma$age[trauma$hospital==3],
#                  ISS=trauma$ISS[trauma$hospital==3],
#                  GCS=trauma$GCS[trauma$hospital==3])
X3 <- subset(trauma, hospital == 3, select = c(sex, age, ISS, GCS))
Lambda3 <- inv.prior.cov(X3, lambda=0.01, L=3, family="binomial")
fit3 <- MAP.estimation(y=trauma$mortality[trauma$hospital==3], X=X3, family="binomial", Lambda=Lambda3)
summary(fit3)
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
    ## (Intercept)  -2.3144  0.4020 -3.1024  -1.5264
    ## sex           0.1580  0.5714 -0.9619   1.2780
    ## age           1.3031  0.2955  0.7239   1.8824
    ## ISS           0.2967  0.2638 -0.2203   0.8138
    ## GCS          -2.4221  0.4130 -3.2316  -1.6125
    ## 
    ## Dispersion parameter (sigma2):  1 
    ##             log Lik Posterior:  -50.63 
    ##                   Convergence:  0

It can be seen that all algorithms have converged (`Convergence: 0`). If
`Convergence: 1` occurs, you can increase the number of iteration in the
optimization process of [`optim()`](https://rdrr.io/r/stats/optim.html)
by adding `control = list(maxit=500)` to the function `MAP.estimation`,
as shown below:

``` r

# Example for Center 3:
fit3 <- MAP.estimation(y=trauma$mortality[trauma$hospital==3], X=X3, family="binomial", Lambda=Lambda3, control = list(maxit=500))
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

The outputs `fit1`,`fit2`, and `fit3` from the local centers should be
sent to the central server for further analysis. To send these lists
from R to the central server (which also uses R), you can save them in a
format that R can easily read, such as an RDS file.

``` r

# Save fit1 as an RDS file
saveRDS(fit1, file="fit1.rds")

# Save fit2 as an RDS file
saveRDS(fit2, file="fit2.rds")

# Save fit3 as an RDS file
saveRDS(fit3, file="fit3.rds")
```

### BFI estimations

Now, the received files can be loaded in R using the following lines:

``` r

# Load the RDS files
fit1 <- readRDS("fit1.rds") # use the relative path to the file
fit2 <- readRDS("fit2.rds") # use the relative path to the file
fit3 <- readRDS("fit3.rds") # use the relative path to the file
```

On the central server, the
[`bfi()`](https://hassanpazira.github.io/BFI/reference/BFI.md) function
can be used to obtain the BFI estimations:

``` r

theta_hats <- list(fit1$theta_hat, fit2$theta_hat, fit3$theta_hat)
A_hats     <- list(fit1$A_hat, fit2$A_hat, fit3$A_hat)
Lambda_com <- inv.prior.cov(X1, lambda=0.01, L=3, family="binomial")
Lambdas    <- list(Lambda1, Lambda2, Lambda3, Lambda_com)
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
    ## (Intercept)  -1.4434  0.2465 -1.9265  -0.9602
    ## sex          -0.2473  0.4187 -1.0680   0.5733
    ## age           1.2189  0.2190  0.7897   1.6481
    ## ISS           0.4939  0.1945  0.1127   0.8751
    ## GCS          -1.7375  0.2491 -2.2258  -1.2492
    ## 
    ## Dispersion parameter (sigma2):  1 
    ## 
    ## Minus the Curvature Matrix: 
    ## 
    ##             (Intercept)     sex     age      ISS      GCS
    ## (Intercept)     30.4330  7.7926  0.7578   8.3797 -16.6608
    ## sex              7.7926  7.8026  1.3061   1.8414  -4.6925
    ## age              0.7578  1.3061 31.5230  -4.8356  15.7494
    ## ISS              8.3797  1.8414 -4.8356  30.7881 -11.7189
    ## GCS            -16.6608 -4.6925 15.7494 -11.7189  34.4508

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
    ## (Intercept)  -1.6255  0.2398 -2.0955  -1.1556
    ## sex          -0.3357  0.3990 -1.1178   0.4463
    ## age           1.3703  0.2056  0.9673   1.7733
    ## ISS           0.5498  0.1862  0.1848   0.9149
    ## GCS          -1.9981  0.2379 -2.4645  -1.5318
    ## 
    ## Dispersion parameter (sigma2):  1 
    ##             log Lik Posterior:  -109.4 
    ##                   Convergence:  0 
    ## 
    ## Minus the Curvature Matrix: 
    ## 
    ##             (Intercept)     sex     age      ISS      GCS
    ## (Intercept)     33.5675  8.4942  2.3840   8.3535 -18.2832
    ## sex              8.4942  8.5042  1.8643   1.4055  -4.5326
    ## age              2.3840  1.8643 37.1516  -6.1748  18.0223
    ## ISS              8.3535  1.4055 -6.1748  33.4408 -12.8543
    ## GCS            -18.2832 -4.5326 18.0223 -12.8543  38.5342

Now, we can see the difference between the BFI and combined estimates:

``` r

# Squared Errors:
(fit_comb$theta_hat - BFI_fits$theta_hat)^2
```

    ##      (Intercept)        sex        age         ISS        GCS
    ## [1,]  0.03319286 0.00781985 0.02292217 0.003132086 0.06791077

which are close to zero, as expected!

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
