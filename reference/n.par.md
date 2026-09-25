# The Number of Predictors, Coefficients, and Observations

`n.par` returns the number of regression parameters, covariates and
observations present in X based on the selected family.

## Usage

``` r
n.par(X, family = c("gaussian", "binomial", "survival"))
```

## Arguments

- X:

  design matrix of dimension \\n \times p\\, where \\n\\ is the number
  of samples observed, and \\p\\ is the number of
  predictors/covariables. It could be a matrix or a list of matrices.

- family:

  a description of the error distribution used to specify the model.
  This should be a character string, either “`gaussian`”, “`binomial`”,
  or “`survival`”. Can be abbreviated. By default the `gaussian` family
  is used.

## Details

`orig.names` and `covar.names` are the same if the all covariates in `X`
are continuous. However, if there are at least one categorical variable
in `X` with more than two categories, they are different.

## Value

`n.par` returns a list containing the following components:

- n.reg.par:

  the number of regression parameters;

- n.covar:

  the number of covariates;

- n.sample:

  the number of samples/observations;

- orig.names:

  the original variable names excluding dummy variable names;

- covar.names:

  the variables names, including any dummy variable names (if
  applicable).

## Author

Hassan Pazira  
Maintainer: Hassan Pazira <h.pazira@arq.org>

## Examples

``` r
#--------------------
# family = "gaussian"
#--------------------

X0 <- data.frame(x1 = rnorm(50),                     # standard normal variable
                 x2 = sample(0:2, 50, replace=TRUE), # categorical variable
                 x3 = sample(0:1, 50, replace=TRUE)) # dichotomous variable
n.par(X0) # without dummy variables
#> $n.reg.par
#> [1] 3
#> 
#> $n.covar
#> [1] 3
#> 
#> $n.sample
#> [1] 50
#> 
#> $orig.names
#> [1] "x1" "x2" "x3"
#> 
#> $covar.names
#> [1] "x1" "x2" "x3"
#> 
X0$x2 <- as.factor(X0$x2)
X0$x3 <- as.factor(X0$x3)
n.par(X0)  # with dummy variables
#> $n.reg.par
#> [1] 4
#> 
#> $n.covar
#> [1] 3
#> 
#> $n.sample
#> [1] 50
#> 
#> $orig.names
#> [1] "x1" "x2" "x3"
#> 
#> $covar.names
#> [1] "x1"  "x21" "x22" "x31"
#> 

X1 <- data.frame(Intercept = rep(1,30),
                 x1 = rnorm(30),                     # continuous variable
                 x2 = sample(0:2, 30, replace=TRUE)) # categorical variable
n.par(X1) # without dummy variables
#> $n.reg.par
#> [1] 2
#> 
#> $n.covar
#> [1] 2
#> 
#> $n.sample
#> [1] 30
#> 
#> $orig.names
#> [1] "x1" "x2"
#> 
#> $covar.names
#> [1] "x1" "x2"
#> 
X1$x2  <- as.factor(X1$x2)
n.par(X1) # without dummy variables
#> $n.reg.par
#> [1] 3
#> 
#> $n.covar
#> [1] 2
#> 
#> $n.sample
#> [1] 30
#> 
#> $orig.names
#> [1] "x1" "x2"
#> 
#> $covar.names
#> [1] "x1"  "x21" "x22"
#> 

# a list of two data sets:
X01 <- list(X0, X1)
n.par(X01)
#> $n.reg.par
#> [1] 4 3
#> 
#> $n.covar
#> [1] 3 2
#> 
#> $n.sample
#> [1] 50 30
#> 
#> $orig.names
#> $orig.names[[1]]
#> [1] "x1" "x2" "x3"
#> 
#> $orig.names[[2]]
#> [1] "x1" "x2"
#> 
#> 
#> $covar.names
#> $covar.names[[1]]
#> [1] "x1"  "x21" "x22" "x31"
#> 
#> $covar.names[[2]]
#> [1] "x1"  "x21" "x22"
#> 
#> 
```
