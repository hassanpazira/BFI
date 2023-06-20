\name{BFI-package}
\alias{BFI-package}
\docType{package}

\title{Bayesian Federated Inference}

\description{
The Bayesian Federated Inference method combines inference results from different (medical) centers without sharing the data. In this version of the package, the user can fit models specifying Gaussian and Binomial (Logistic) families. The package will be updated with more models soon.
}
\details{
\tabular{ll}{
Package: \tab BFI\cr
Type: \tab Package\cr
Version: \tab 0.6.4\cr
Date: \tab 2023-07-16\cr
License: \tab GPL (>=2)
}

\code{MAP.estimation} and \code{bfi} are the main functions. All other functions are utility functions.

Some examples are provided in the vignettes accompanying this package in order to show how the package can be applied to real data. The vignettes can be found on the package website at \url{https://hassanpazira.github.io/BFI/} or from within R once the package has been installed, e.g. via \code{vignette("BFI", package = "BFI")}.
}
\author{
Hassan Pazira, Marianne A. Jonker, Anthony C.C. Coolen\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}
\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023). \emph{Bayesian Federated Inference for Statistical Models}, \emph{Statistics in Medicine}, Vol. 0(0), 0-0. <https://doi.org/10.48550/arXiv.2302.07677>
}

\keyword{Bayesian}
\keyword{Federated}
\keyword{package}
