\name{inv.prior.cov}
\alias{inv.prior.cov}
\title{Creates an inverse covariance matrix for a Gaussian prior}
\description{
\code{inv.prior.cov} builds a diagonal inverse covariance matrix for the Gaussian prior distribution based on the design matrix of covariates, that takes into account the number of regression parameters in case of categorical covariates.  In case of a linear model, it also includes the variance of the measurement errors.
}
\usage{
inv.prior.cov(X, lambda = 1, family = gaussian, intercept = TRUE,
              independ = TRUE, set_seed = NULL)
}
\arguments{
\item{X}{design matrix of dimension \eqn{n \times p}, where \eqn{n} is the number of samples observed, and \eqn{p} is the number of model prameters except intercept.}

\item{lambda}{the vector used as the diagonal of the (inverse covariance) matrix created by \code{inv.prior.cov()}. It could be either a single positive number, a vector of two elements, or a vector with length equal to \eqn{p}, \eqn{p+1} or \eqn{p+2} (where \eqn{p} is the number of columns of \code{X} argument) depends on the arguments \code{X}, \code{family}, and \code{intercept}. If \code{lambda} is a scalar, the function \code{inv.prior.cov} takes the same value (equal to that scalar value) for each model parameter. Default is \code{lambda=1}. In the linear regression model the last parameter is assumed to be the variance of measurement error. If \code{lambda} is two dimensional, the first value is used for  the regression parameters and the second for the variance of the measurement error. Used when \code{independ=TRUE}.}

\item{family}{a description of the error distribution and link function used to specify the model. This can be a character string naming a family function or the result of a call to a family function (see \code{\link{family}} for details). In the current version, the family of model can be \code{gaussian} (with \code{identity} link function) and \code{binomial} (with \code{logit} link function). By default the \code{gaussian} family is used. In case of a linear regression model, \code{family=gaussian}, there is an extra model parameter for the variance of measurement error.}

\item{intercept}{logical flag for having an intercept. If \code{intercept=TRUE} (the default), the matrix created by \code{inv.prior.cov} will have one additional row and column coresponding to the intercept, while if \code{intercept=FALSE}, the matrix does not have those row and column.}
}

\details{
\code{inv.prior.cov} creates a diagonal matrix with the vector \code{lambda} as its diagonal.
The output of \code{inv.prior.cov()} is used in the main functions \code{MAP.estimation()} and \code{bfi()}.
}

\value{
\code{inv.prior.cov} returns a diagonal matrix. The dimension of the matrix depends on the arguments \code{X}, \code{family}, and \code{intercept}.
}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023). \emph{Bayesian Federated Inference for Statistical Models},
{Statistics in Medicine}, Vol. 0(0), 0-0. <https://doi.org/10.48550/arXiv.2302.07677>
}

\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{MAP.estimation}}
}

\examples{
#-------------
# y ~ Binomial
#-------------
X <- data.frame(matrix(rnorm(50 * 4), 50, 4))
# We assume the same (inverse) covariance value (0.05) for all 5 parameters
# y ~ Binomial with intercept
inv.prior.cov(X, 0.05, family=binomial)
# y ~ Binomial without intercept
inv.prior.cov(X, 0.05, family="binomial", intercept = FALSE)
# y ~ Gaussian with intercept
inv.prior.cov(X, 0.05, family=gaussian)
}
\keyword{Bayesian}
