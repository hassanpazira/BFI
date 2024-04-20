# BFI v3.1.0 (Latest)
* The `MAP.estimation` function has been revised, and corresponding updates have been made to its documentation and examples.
* This version marks the final contribution of Hassan Pazira as the package maintainer. The new maintainer for future versions is Marianne A. Jonker <Marianne.Jonker@radboudumc.nl>

# BFI v3.0.1
* In this release, the package has been updated to support both observational and randomized trial data for estimating treatment effects using Bayesian Federated 'Causal' Inference.

# BFI v2.1.0
* In this release, the package has been updated to support a (semi-parametric) Cox model with an 'unspecified' baseline hazard, maximizing partial log-likelihood.
* The documentation and examples have been updated.

# BFI v2.0.1
* This is a bugfix release to resolve a minor bug in the `optim` function entries related to the `gaussian` family.

# BFI v2.0.0
* The package now supports survival data analysis, adding comprehensive tools for time-to-event modeling alongside existing GLM functionalities.
* The documentation and examples have been updated.

# BFI v1.1.4
* The package has been updated to address a few NOTEs identified during the CRAN submission.

# BFI v1.1.1
* The package was prepared for submission to CRAN with minor modifications.

# BFI v1.1.0
* The functions `inv.prior.cov()` and `bfi()` were adapted when there is a center specific variable.
* The function `summary()` was updated to be used in the case of stratification.
* The manual pdf was updated in the case of center specific variable.
* Henceforth the `BFI` package can be called from `Python`. 
* A vignette was also added to the package for calling `BFI` from `Python`.

# BFI v1.0.0
+ In this release, the package website was built using the pkgdown package.
* Most of the functions were adapted with extra arguments in the case of stratification.
* The manual pdf was also updated.

# BFI v0.6.4
* The manual pdf was updated.
* This is also a bugfix release to resolve one minor bug.

# BFI v0.6.3
* In this release, the package was updated so that the outputs of `inv.prior.cov()` and `MAP.estimation()` (for different centers) have the same dimensions when `intercept=FALSE` to be used in `bfi()`.
* Moreover, one argument in `bfi()`, i.e., `const_var`, is added to the package to handle the constant variables.
* This is a bugfix release to resolve one minor bug as well.

# BFI v0.5.2
* In this release, the summary function (a S3 method for class `bfi`) was added to the package. 
* A package pdf manual was created for the package.

# BFI v0.4.2
* In this release, the package was updated with several arguments, e.g., an intercept should be fitted or not.
* A vignette is added to the package as well.

# BFI v0.3.2
* In this release, the package can carry out the stratified analysis.

# BFI v0.2.2
* This is a bugfix release to resolve one minor bug, and add two functions related to building Gamma matrix.
* Moreover, the package now handles the categorical covariates with more than two levels.


