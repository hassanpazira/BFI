# BFI v1.1.1 (Latest)
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
* In this release, the package is updated so that the outputs of `inv.prior.cov()` and `MAP.estimation()` (for different centers) have the same dimensions when `intercept=FALSE` to be used in `bfi()`.
* Moreover, one argument in `bfi()`, i.e., `const_var`, is added to the package to handle the constant variables.
* This is a bugfix release to resolve one minor bug as well.

# BFI v0.5.2
* In this release, the summary function (a S3 method for class `bfi`) was added to the package. 
* A package pdf manual was created for the package.

# BFI v0.4.2
* In this release, the package is updated with several arguments, i.e., an intercept should be fitted or not.
* A vignette is added to the package as well.

# BFI v0.3.2
* In this release, the package can carry out the stratified analysis.

# BFI v0.2.2
* This is a bugfix release to resolve one minor bug, and add two functions related to building Gamma matrix.
* Moreover, the package now handles the categorical covariates with more than two levels.


