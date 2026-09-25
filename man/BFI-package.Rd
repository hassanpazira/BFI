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
Version: \tab 3.2.0\cr
License: \tab MIT + file LICENSE
}

\code{MAP.estimation} and \code{bfi} are the main functions. All other functions are utility functions.

Some examples are provided in the vignettes accompanying this package in order to show how the package can be applied to real data. The vignettes can be found on the package website at \url{https://hassanpazira.github.io/BFI/} or within R once the package has been installed, e.g., via \code{vignette("BFI", package = "BFI")}.
}
\author{
Hassan Pazira, Emanuele Massa, Marianne A. Jonker\cr
Maintainer: \packageMaintainer{BFI}
}
\references{
Pazira H. and Jonker M.A. (2026).
\emph{BFI: An R Package for Bayesian Federated Inference}.
<https://arxiv.org/abs/2609.27977>

Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker M.A. (2026).
\emph{Bayesian federated inference for survival models},
\emph{Journal of Applied Statistics}, 53(2): 203-223.
<https://doi.org/10.1080/02664763.2025.2511932>

Jonker M.A., Pazira H. and Coolen A.C.C. (2025). \emph{Bayesian Federated Inference for regression models based on non-shared medical center data}, \emph{Research Synthesis Methods}, 16(2): 383-423. <https://doi.org/10.1017/rsm.2025.6>

Jonker M.A., Pazira H. and Coolen A.C.C. (2024). \emph{Bayesian federated inference for estimating statistical models based on non-shared multicenter data sets}, \emph{Statistics in Medicine}, 43(12): 2421-2438. <https://doi.org/10.1002/sim.10072>
}

\keyword{package}
