\name{estimators_maker}
\alias{estimators_maker}
\title{Bayesian Federated Inference}
\description{
\code{estimators_maker} function is used (in local centers) to estimate Maximum A Posterior (MAP) for the parameters of the GLM and Survival models.
}
\usage{
estimators_maker(y, X, family=c("binomial","gaussian"), Gamma)
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
\code{\link{bfi}}.
}
\examples{
set.seed(1123)

#-------------
# y ~ Gaussian
#-------------

n         <- 30
p         <- 4    # number of coefficients
X         <- matrix(rnorm(n * p), n, p)
eta       <- 1 + 2 * X[,1]  ## with an intercept b0=1, b1=2, b2=b3=b4=0
mu        <- gaussian()$linkinv(eta)
lambda    <- 0.01
sigma2_e  <- 1.5
# true_theta <- c(1, 2, 0, 0, 0, sigma2_e)
Gamma     <- diag(c(rep(lambda, p+1), sigma2_e), p+2) # inverse of covariance matrix
y         <- rnorm(n, mu, sd=sigma2_e)

# MAP estimates of the parameters of interest and curvature matrix:
(fit       <- estimators_maker(y, X, family="gaussian", Gamma))
}
