\name{BFI-package}
\alias{BFI-package}
\docType{package}
\title{Bayesian Federated Inference}
\description{
The Bayesian Federated Inference method combines inference results from different (medical) centers without sharing the data. In this version of the package, the user can fit models specifying Gaussian, Binomial (Logistic) and Survival families.
}
\details{
\tabular{ll}{
Package: \tab BFI\cr
Type: \tab Package\cr
Version: \tab 3.1.0\cr
Date/Publication: \tab 2024-04-27\cr
License: \tab GPL (>=2)
}

\code{MAP.estimation} and \code{bfi} are the main functions. All other functions are utility functions.

Some examples are provided in the vignettes accompanying this package in order to show how the package can be applied to real data. The vignettes can be found on the package website at \url{https://hassanpazira.github.io/BFI/} or within R once the package has been installed, e.g., via \code{vignette("BFI", package = "BFI")}.
}
\author{
Hassan Pazira, Emanuele Massa, Marianne A. Jonker\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}
\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2024). \emph{Bayesian federated inference for estimating statistical models based on non-shared multicenter data sets}, \emph{Statistics in Medicine}, 43(12): 2421-2438. <https://doi.org/10.1002/sim.10072>

Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker M.A. (2025b). \emph{Bayesian Federated Inference for Survival Models}, \emph{Journal of Applied Statistics (Accepted)}. <https://arxiv.org/abs/2404.17464>

Jonker M.A., Pazira H. and Coolen A.C.C. (2025a). \emph{Bayesian Federated Inference for regression models based on non-shared medical center data}, \emph{Research Synthesis Methods}, 1-41. <https://doi.org/10.1017/rsm.2025.6>
}

\keyword{package}
