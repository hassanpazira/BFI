<!-- badges: start -->
[![Build
Status](https://ci.appveyor.com/api/projects/status/github/hassanpazira/bfi?branch=master&svg=true)](https://ci.appveyor.com/project/hassanpazira/bfi/branch/master)
[![codecov.io](https://codecov.io/github/hassanpazira/bfi/coverage.svg?branch=master)](https://app.codecov.io/github/hassanpazira/bfi/branch/master)
[![MIT
license](https://img.shields.io/badge/license-MIT-brightgreen.svg)]( https://opensource.org/license/mit/)
<!-- badges: end -->

# `BFI` <img src="./man/figures/Last_BFI.jpg" align="right" width="170px"/>

> #### Bayesian Federated Inference

## Description

Due to the limited size of the available data sets especially in rare diseases, it is sometimes challenging to identify the most relevant predictive features using multivariable statistical analysis. This issue may be resolved by combining data from multiple centers into one centralized location without sharing their data with each other, but doing so is difficult in reality because of privacy and security concerns.

To address these challenges, we developed and implemented a Bayesian Federated Inference (BFI) framework for multicenter data. It aims to leverage the statistical power of larger (combined) data sets without requiring all the data to be aggregated in one location. The BFI framework allows each center using their own local data to infer the optimal parameter values as well as additional features of the posterior parameter distribution to be able to gather more information which is not captured by alternative techniques. One of the benefit of BFI over alternative approaches is that, only one inference cycle across the centers is required in BFI.

An R package called `BFI` is created to perform Bayesian Federated Inference. The following instructions will install the development version of the `BFI` package to a computer.


## Install R and RStudio

First, you need to install R and RStudio:

* Install [R](https://www.R-project.org/)
* Install [RStudio Desktop](https://posit.co/download/rstudio-desktop/) (once you have R installed)

For more details about installing R and RStudio, see [this page](https://andreashandel.github.io/MADAcourse/Tools_RandRStudio.html).
If you need help learning R, see [RStudio Education](https://education.rstudio.com/learn/).


## Install `BFI` package

In order to install the `BFI` package directly from Github, you need to have the **devtools** package.

Invoke R or RStudio and then type (in Console)

``` r
if(!require(devtools)) {install.packages("devtools")}
```

and then load it by typing:

``` r
library(devtools)
```

Next, install `BFI` as follows:

``` r
devtools::install_github("hassanpazira/BFI", force = TRUE)
```

The package can now be loaded into R and used by:

``` r
library(BFI)
```

## Update

The latest version of the `BFI`package is `0.6.4`. To check the current version of `BFI` installed in your R library, use:

``` r
packageVersion("BFI")
```

## Details

The `BFI` package provides several functions, the most important of which are the following two main functions:

-   `MAP.estimation()`: should be used by the centers, and the result should be sent to a central server.

-   `bfi()`: should be used by a central server.

To access the R documentation for these functions, for example `bfi()`, enter the following command:

``` r
help(bfi, package = "BFI")   
# or, equivalently, after loading the BFI package 
?bfi
```


## Usage

Let's look at the following example to see how the `BFI` package can be used. For more examples and details look at the `BFI` vignette by typing

``` r
devtools::install_github("hassanpazira/BFI", dependencies = TRUE, build_vignettes = TRUE, force = TRUE)
browseVignettes("BFI")  # to see all vignettes from the BFI package in an HTML browser.
```
or use `vignette("BFI")` to see the `BFI` vignette in the Help tab of RStudio.


Now, we generate two independent (local) data sets from Gaussian distribution, and then apply the package to see how it works. First apply the function `MAP.estimation()` to each local data, and then apply the `bfi()` function to the aggregated results.

``` r
# Load the BFI package
library(BFI)

# model assumption:
beta <- 1:4  # regression coefficients (beta[1] = 1 is the intercept)

#-----------------------------------------------------
# Data Simulation for local center 1 when y ~ Binomial
#-----------------------------------------------------
n1 <- 30                                           # sample size of center 1
X1 <- data.frame(x1=rnorm(n1),                     # continuous variable
                 x2=sample(0:2, n1, replace=TRUE)) # categorical variable
# make dummy variables
X1x2_1 <- ifelse(X1$x2 == '1', 1, 0)
X1x2_2 <- ifelse(X1$x2 == '2', 1, 0)
X1$x2  <- as.factor(X1$x2)
# linear predictor:
eta1   <- beta[1] + X1$x1 * beta[2] + X1x2_1 * beta[3] + X1x2_2 * beta[4]
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu1    <- binomial()$linkinv(eta1)
y1     <- rbinom(n1, 1, mu1)

#-----------------------------------------------------
# Data Simulation for local center 2 when y ~ Binomial
#-----------------------------------------------------
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

# assume the same inverse covariance matrix (Lambda) for both centers:
Lambda <- inv.prior.cov(X1, lambda=0.01, family=binomial)

#--------------------------
# MAP estimates at center 1
#--------------------------
fit1       <- MAP.estimation(y1, X1, family=binomial, Lambda)
theta_hat1 <- fit1$theta_hat # intercept and coefficient estimates for center 1
A_hat1     <- fit1$A_hat     # curvature matrix for center 1

#--------------------------
# MAP estimates at center 2
#--------------------------
fit2       <- MAP.estimation(y2, X2, family=binomial, Lambda)
theta_hat2 <- fit2$theta_hat # intercept and coefficient estimates for center 2
A_hat2     <- fit2$A_hat     # curvature matrix for center 2

#----------------------
# BFI at central center
#----------------------
A_hats     <- list(A_hat1, A_hat2)
theta_hats <- list(theta_hat1, theta_hat2)
bfi_fit    <- bfi(theta_hats, A_hats, Lambda)
class(bfi_fit)
summary(bfi_fit, cur_mat = TRUE)

#--------------------
# stratified analysis
#--------------------
bfi(theta_hats, A_hats, Lambda, stratified = TRUE, strat_par = 1L)

```

## Citation

To cite `BFI` in publications, please use:

``` r
citation("BFI")
```


## Documentation

Here are some of technical papers of the package:

-   [Generalized Linear Models (GLMs)](https://arxiv.org/abs/2302.07677)

-   [Survival Models](https://arxiv.org/abs/2302.07677)


## Contact

If you find any errors, have any suggestions, or would like to request that something be added, please file an issue at [issue report](https://github.com/hassanpazira/BFI/issues/) or send an email to: hassan.pazira@radboudumc.nl.

