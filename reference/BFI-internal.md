# Internal BFI Functions

Internal BFI functions.

## Value

`A.l.maker()` returns a curvature matrix for each local.

`negloglik.theta()` returns the negative log-likelihood of regression
coefficients and error variance, which can be utilized for optimization
purposes.

`model.maker()` returns the outcome and design matrix, incorporating
dummy variables when categorical covariates are present, e.g., by
expanding factors to a set of dummy variables and expanding interactions
similarly.

`optim.survival()` is a general-purpose optimization algorithm based on
the “L-BFGS-B” method, and returns a vector of estimates for the
coefficients and baseline hazard parameters.

`ql.LRT()` returns an optimal order/degree of the exponentiated
polynomial model, obtained by the likelihood ratio test, for the
'survival' family.

`lambda.poly()` returns the baseline hazard function for the
exponentiated polynomial form.

`i.basis()` creates the integral of the output of basis when
`ibasis = TRUE`. This is also an indicator function between the
\\k^{th}\\ and \\(k+1)^{th}\\ time points when `ibasis = FALSE`.

## Details

These functions are not intended for use by users.

## Author

Hassan Pazira  
Maintainer: Hassan Pazira <h.pazira@arq.org>
