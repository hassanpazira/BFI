# BFI <img src="./Last_BFI.jpg" align="right" width="160px"/>

> Bayesian Federated Inference

## Description

Due to the limited size of the available data sets especially in rare diseases, it is sometimes challenging to identify the most relevant predictive features using multivariable statistical analysis. This issue may be resolved by combining data from multiple centers into one centralized location without sharing their data with each other, but doing so is difficult in reality because of privacy and security concerns.

To address these challenges, we developed and implemented a Bayesian Federated Inference (BFI) framework for multicenter data. It aims to leverage the statistical power of larger (combined) data sets without requiring all the data to be aggregated in one location. The BFI framework allows each center using their own local data to infer the optimal parameter values as well as additional features of the posterior parameter distribution to be able to gather more information which is not captured by alternative techniques. One of the benefit of BFI over alternative approaches is that, only one inference cycle across the centers is required in BFI.

An R package called `BFI` is created to perform Bayesian Federated Inference. The following instructions will install the development version of the `BFI` package to a computer.


## Install R and RStudio

First, you need to install R and Rstudio:

* Install [R](http://www.r-project.org/)
* Install [RStudio Desktop](https://posit.co/download/rstudio-desktop) (once you have R installed)

For more details about installing R and RStudio, see [this page](https://andreashandel.github.io/MADAcourse/Tools_RandRStudio.html).
If you need help learning R, see [RStudio Education](https://education.rstudio.com/learn/).


## Install BFI package

In order to install the package directly from Github, you need to have the **devtools** package.

Invoke R or RStudio and then type

``` r
if(!require(devtools)) {install.packages("devtools")}
```

and then load it by typing:

``` r
library(devtools)
```

Next, install `BFI` as follows:

``` r
devtools::install_github("hassanpazira/BFI")
```

The package can now be loaded into R and used by:

``` r
library(BFI)
```

## Update

The latest version of the `BFI`package is `0.2.2`. To check the current version of `BFI` installed in your R library, use:

``` r
packageVersion("BFI")
```

## Details

The current version of the `BFI` package provides two main functions:

-   `MAP.estimation()`: should be used by the centers, and the result should be sent to a central server.

-   `bfi()`: should be used by a central server.

To access the R documentation for these functions, for example `bfi()`, enter the following command:

``` r
help(bfi, package = "BFI")   # or ?bfi after loading the BFI package
```


## Usage

Let's look at the following example to see how the `BFI` package can be used.

We generate two independent (local) data sets from Gaussian distribution, and then apply the package to see how it works. First apply the function `MAP.estimation()` to each local data, and then apply the `bfi()` function to the aggregated results.

``` r
# Load BFI package
library(BFI)

## Create local center 1
set.seed(1123581)
n1       <- 30     # sample size of center 1
p        <- 4      # p (number of coefficients and intercept) is equal for all centers
X1       <- data.frame(matrix(rnorm(n1 * p), n1, p))
eta1     <- 1 + 2 * X1[,1]  # with an intercept: b0=1, b1=2, and b2=b3=...=bp=0
mu1      <- gaussian()$linkinv(eta1)
lambda   <- 0.01
sigma2_e <- 0.75   # nuisance parameter
# Gamma (inverse of covariance matrix for prior) is assumed to be equal across centers
Gamma    <- Gamma.maker(X1, lambda=c(lambda,sigma2_e), family="gaussian")
y1       <- rnorm(n1, mu1, sd=sigma2_e)
# Obtain the all required estimates for the center
fit1     <- MAP.estimation(y1, X1, family="gaussian", Gamma)

## Create local center 2
n2   <- 40     # sample size of center 2
X2   <- data.frame(matrix(rnorm(n2 * p), n2, p))
eta2 <- 1 + 2 * X2[,1] # coefficients are equal for all centers
mu2  <- gaussian()$linkinv(eta2)
y2   <- rnorm(n2, mu2, sd=sigma2_e)
# Obtain the all required estimates for the center
fit2 <- MAP.estimation(y2, X2, family="gaussian", Gamma)

## BFI estimates
A_hats     <- list(fit1$A_hat, fit2$A_hat)
theta_hats <- list(fit1$theta_hat, fit2$theta_hat)
BFI_fit    <- bfi(theta_hats, A_hats, Gamma)

# Coefficients and nuisance estimates
t(theta_hat_bfi <- BFI_fit$theta_hat)
#      (Intercept)       X1         X2         X3         X4  sigma2_e
# [1,]   0.9613874 2.015488 0.09699617 0.05426457 0.04615605 0.5868312

# Curvature matrix estimate
(A_bfi <- BFI_fit$A_hat)

# SD of the estimates
(sd_bfi <- BFI_fit$sd)
```

## Citation

To cite `BFI` in publications, please use:

``` r
citation("BFI")
```


## Documentation

Here are some of technical papers of the package:

-   [Generalized Linear Models (GLMs)](https://arxiv.org/abs/2302.07677)

-   [Survival Models]()


## Contact

If you find any errors, have any suggestions, or would like to request that something be added, please let us know using [issues report](https://github.com/hassanpazira/BFI/issues) or sending an email to: hassan.pazira@radboudumc.nl.

