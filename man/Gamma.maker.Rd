\name{Gamma.maker}
\alias{Gamma.maker}
\title{Build the inverse of a covariance matrix for }
\description{
\code{Gamma.maker} function builds a matrix which can be used as the inverse of a covariance matrix for Gaussian prior distribution, to be used in the main functions \code{MAP.estimation} and \code{bfi}.
}
\usage{
Gamma.maker(X, lambda = 1, family = c("binomial", "gaussian"), independ = T, set_seed = NULL)
}
\arguments{
}

\details{
}

\value{
}
\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023) <https://doi.org/10.48550/arXiv.2302.07677>
\emph{Bayesian Federated Inference for Statistical Models},
\emph{Statistics in Medicine}, Vol. 0(0), 0-0.
}
\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}}
\seealso{
\code{\link{MAP.estimation}}.
}
\examples{
#-------------
# y ~ Binomial
#-------------
X      <- data.frame(matrix(rnorm(50 * 4), 50, 4))
lambda <- 0.05
# We assume the same (inverse) covariance matrix for all locations
# and consider independency of prior parameters
(Gamma <- Gamma.maker(X, lambda, family="binomial"))

#-------------
# y ~ Gaussian
#-------------
X         <- data.frame(matrix(rnorm(50 * 3), 50, 3))
lambda    <- 0.01
sigma2_e  <- 0.5
# If we consider dependency for all prior parameters:
(Gamma    <- Gamma.maker(X, lambda=c(lambda, sigma2_e), family="gaussian", independ = F, set_seed=1123))
}
