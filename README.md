
# `BFI` <img src="./man/figures/BFI_Hexagon_blue2.png" align="right" width="110px"/>

> #### Bayesian Federated Inference

<!-- badges: start -->
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/BFI)](https://cran.r-project.org/package=BFI)
[![CRAN total downloads](https://cranlogs.r-pkg.org/badges/grand-total/BFI)](https://cran.r-project.org/package=BFI)
[![CRAN monthly downloads](https://cranlogs.r-pkg.org/badges/BFI)](https://CRAN.R-project.org/package=BFI)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT license](https://img.shields.io/badge/license-MIT-brightgreen.svg)]( https://opensource.org/license/mit)
[![code size](https://img.shields.io/github/languages/code-size/hassanpazira/BFI.svg)](https://github.com/hassanpazira/BFI)
<!-- badges: end -->

## Description

Due to the limited size of the available data sets, especially in rare diseases, identifying the most relevant predictive features using multivariable statistical analysis can be challenging. This issue may be resolved by combining data from multiple centers into one centralized location without sharing their data with each other, but doing so is difficult in reality because of privacy and security concerns.

To address these challenges, we developed and implemented a Bayesian Federated Inference (BFI) framework for multicenter data. It aims to leverage the statistical power of larger (combined) data sets without requiring all the data to be aggregated in one location. The BFI framework enables each center to use its own local data to infer the optimal parameter values and additional features of the posterior distribution, providing insights not captured by alternative techniques. One of the benefits of BFI over alternative approaches is that only one inference cycle across the centers is required in BFI.

An R package called `BFI` has been created to perform Bayesian Federated Inference. The following instructions will install the development version of the `BFI` package to a computer.
Python and SAS users can also apply the BFI methodology in their respective environments.
For instructions, see [here for Python](https://hassanpazira.github.io/BFI/articles/Python.html) and [here for SAS](https://hassanpazira.github.io/BFI/articles/SAS.html).


## Install R and RStudio

First, you need to install R and RStudio:

* Install [R](https://www.R-project.org/)
* Install [RStudio Desktop](https://posit.co/download/rstudio-desktop/) (once you have R installed)

For more details about installing R and RStudio, see [this page](https://andreashandel.github.io/MADAcourse/content/module-intro-tools/tools-randrstudio.html).
If you need help learning R, see [RStudio Education](https://education.rstudio.com/learn/).


## Install `BFI` package

In order to install the `BFI` package, invoke R or RStudio and then follow one of the following steps:


### Installing from R
To install and load the `BFI` package directly from R, type the following (in Console)
``` r
install.packages("BFI")
library(BFI)
```

### Installing from GitHub
To install the `BFI` package directly from GitHub, you need to have the **devtools** package. So type the following to install and then load it:

``` r
install.packages("devtools")
library(devtools)
```

Next, install `BFI` as follows:

``` r
devtools::install_github("hassanpazira/BFI", dependencies = TRUE, build_vignettes = TRUE, force = TRUE)
```

The package can now be loaded into R and used by:

``` r
library(BFI)
```

## Update

The latest version of the `BFI`package is `3.1.0`. To check the current version of `BFI` installed in your R library, use:

``` r
packageVersion("BFI")
```

## Details

The `BFI` package provides several functions, the most important of which are the following two main functions:

-   `MAP.estimation()`: should be used by the centers, and the result should be sent to a central server.

-   `bfi()`: should be used by a central server.

To access the R documentation for these functions, for example `bfi()`, enter the following command:

``` r
help(bfi, package = "BFI")  # without loading the BFI package
# or, equivalently, after loading the BFI package 
?bfi
```


## Usage

Let's look at the following example to see how the `BFI` package can be used. For more examples and details look at the `BFI` vignette by typing

``` r
browseVignettes("BFI")  # to see all vignettes from the BFI package in an HTML browser.
```
or use `vignette("BFI")`, `vignette("SAS")` or `vignette("Python")` to see the `BFI`, `SAS` or `Python` vignettes separately in the Help tab of RStudio.


Now, we generate two independent (local) data sets from Gaussian distribution, and then apply the package to see how it works. First apply the function `MAP.estimation()` to each local data, and then apply the `bfi()` function to the aggregated results.

``` r
#-------------
# y ~ Gaussian
#-------------
# model assumptions:
p     <- 3                     # number of coefficients without intercept
theta <- c(1, rep(2, p), 1.5)  # regression coefficients (theta[1] is the intercept) and sigma2 = 1.5

#-----------------------------------
# Data simulation for local center 1
#-----------------------------------
n1   <- 30                                       # sample size of center 1
X1   <- data.frame(matrix(rnorm(n1 * p), n1, p)) # continuous variables
# linear predictor:
eta1 <- theta[1] + as.matrix(X1) %*% theta[2:4]
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu1  <- gaussian()$linkinv(eta1)
y1   <- rnorm(n1, mu1, sd = sqrt(theta[5]))

#-----------------------------------
# Data simulation for local center 2
#-----------------------------------
n2   <- 50                                       # sample size of center 2
X2   <- data.frame(matrix(rnorm(n2 * p), n2, p)) # continuous variables
# linear predictor:
eta2 <- theta[1] + as.matrix(X2) %*% theta[2:4]
# inverse of the link function:
mu2  <- gaussian()$linkinv(eta2)
y2   <- rnorm(n2, mu2, sd = sqrt(theta[5]))

#---------------------
# Load the BFI package
#---------------------
library(BFI)

#---------------------------
# Inverse Covariance Matrix
#---------------------------
# Creating the inverse covariance matrix for the Gaussian prior distribution:
Lambda <- inv.prior.cov(X1, lambda=0.05, family='gaussian') # the same for both centers

#--------------------------
# MAP estimates at center 1
#--------------------------
fit1       <- MAP.estimation(y1, X1, family='gaussian', Lambda)
theta_hat1 <- fit1$theta_hat # intercept and coefficient estimates
A_hat1     <- fit1$A_hat     # minus the curvature matrix

#--------------------------
# MAP estimates at center 2
#--------------------------
fit2       <- MAP.estimation(y2, X2, family='gaussian', Lambda)
theta_hat2 <- fit2$theta_hat
A_hat2     <- fit2$A_hat

#----------------------
# BFI at central center
#----------------------
A_hats     <- list(A_hat1, A_hat2)
theta_hats <- list(theta_hat1, theta_hat2)
bfi        <- bfi(theta_hats, A_hats, Lambda)
summary(bfi, cur_mat=TRUE)

#--------------------
# Stratified analysis
#--------------------
# Stratified analysis when 'intercept' varies across two centers:
newLambda1 <- inv.prior.cov(X1, lambda=c(0.1, 0.3), family='gaussian', stratified=TRUE, strat_par = 1)
# 'newLambda1' is used the prior for combined data and 'Lambda' is used the prior for locals
bfi1 <- bfi(theta_hats, A_hats, list(Lambda, newLambda1), stratified=TRUE, strat_par=1)
summary(bfi1, cur_mat=TRUE)

# Stratified analysis when 'sigma2' varies across two centers:
newLambda2 <- inv.prior.cov(X1, lambda=c(0.1, 0.3), family='gaussian', stratified=TRUE, strat_par = 2)
# 'newLambda2' is used the prior for combined data and 'Lambda' is used the prior for locals
bfi2 <- bfi(theta_hats, A_hats, list(Lambda, newLambda2), stratified=TRUE, strat_par=2)
summary(bfi2, cur_mat=TRUE)

# Stratified analysis when 'intercept' and 'sigma2' vary across 2 centers:
newLambda3 <- inv.prior.cov(X1, lambda=c(0.1, 0.2, 0.3), family='gaussian', stratified=TRUE, strat_par = c(1, 2))
# 'newLambda3' is used the prior for combined data and 'Lambda' is used the prior for locals
bfi3 <- bfi(theta_hats, A_hats, list(Lambda, newLambda3), stratified=TRUE, strat_par=1:2)
summary(bfi3, cur_mat=TRUE)

```

## Citation

To cite `BFI` in publications, please use:

``` r
citation("BFI")
```


## Documentation

Here are some of technical papers of the package:

-   [BFI for Generalized Linear Models (GLMs)](https://doi.org/10.1002/sim.10072)

-   [BFI for Survival Models](https://arxiv.org/abs/2404.17464)

-   [BFI in Heterogeneous Populations](https://doi.org/10.1017/rsm.2025.6)


## Contact

If you find any errors, have any suggestions, or would like to request that something be added, please file an issue at [issue report](https://github.com/hassanpazira/BFI/issues/) or send an email to: hassan.pazira@radboudumc.nl.

