
# BFI4  <img src="./Last_BFI.jpg" align="right" width="230px"/>
> Bayesian Federated Inference


## What is BFI?

Due to the limited size of the available data sets especially in rare diseases, it is sometimes challenging to identify the most relevant predictive features using multivariable statistical analysis. This issue may be resolved by combining data from multiple centers into one centralized location without sharing their data with each other, but doing so is difficult in reality because of privacy and security concerns.

To address these challenges, we developed and implemented a Bayesian Federated Inference (BFI) framework for multicenter data. It aims to leverage the statistical power of larger (combined) data sets without requiring all the data to be aggregated in one location. The BFI framework allows each center using their own local data to infer the optimal parameter values as well as additional features of the posterior parameter distribution to be able to gather more information which is not captured by alternative techniques.
One of the benefit of BFI over alternative approaches is that, only one inference cycle across the centers is required in BFI. 

An R package called `BFI` is created to perform Federated Bayesian Inference.
The following instructions will install the development version of the `BFI` package to your computer.



## How do I install the package?

In order to install the package directly from Github, you need to have the **devtools** package. 

Invoke R and then type

```r
if(!require(devtools)) {install.packages("devtools")}
```

and then load it by typing:

```r
library(devtools)
```

Next, install `BFI` as follows:

```r
devtools::install_github("hassanpazira/BFI")
```

The package can now be loaded into R and used:

```r
library(BFI)
```


## Documentation:

Here are some of technical papers of the package: 

- [Generalized Linear Models (GLMs)](https://arxiv.org/abs/2302.07677)

- [Survival Models]()


