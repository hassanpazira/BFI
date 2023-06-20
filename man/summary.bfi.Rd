\name{summary.bfi}
\alias{summary}
\alias{summary.bfi}

\title{Summarizing BFI Fits}

\description{Summary method for an object with class 'bfi' created by the \code{MAP.estimation} function.}

\usage{
\method{summary}{bfi}(object, curmat = FALSE,
        digits = max(3, getOption("digits") - 3))
}

\arguments{
\item{object}{fitted \code{bfi} object.}

\item{curmat}{logical; if \code{TRUE}, the curvature matrix around the estimated parameters is returned and printed. Default is \code{FALSE}.}

\item{digits}{significant digits in printout.}
}

\details{
\code{summary.bfi} gives information about the MAP estimates of parameters of the model. It can be used for \code{bfi} objects built by the \code{MAP.estimation} function.

The output of the summary method shows the details of the model, i.e. formula, family and link function used to specify the generalized linear model, followed by information about the estimates, standard deviations and confidence intervals.
Information about the log-likelihood and convergence status are also provided.

By default, \code{summary.bfi} function does not return the curvature matrix, but the user can use \code{curmat=TRUE} to print it.
}

\value{

\code{summary.bfi} returns an object of class \code{summary.bfi}, a list with the following components:

\item{theta_hat}{the component from \code{object}. The last element of this vector is the estimate of the dispersion parameter (sigma2). See the \code{\link{MAP.estimation}} function.}

\item{A_hat}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{sd}{the component from \code{object}. The last element of this vector is the square root of the estimated dispersion. See the \code{\link{MAP.estimation}} function.}

\item{Lambda}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{formula}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{n}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{np}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{value}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{family}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{convergence}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{control}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{logLik}{the value of the loglikelihood function corresponding to estimates (\code{theta_hat}). This is the minus of the \code{value} component.}

\item{link}{the link function. By default the \code{gaussian} family with \code{identity} link function and the \code{binomial} family with \code{logit} link function are used.}

\item{dispersion}{the estimated variance of the random error, i.e., sigma2. The \code{dispersion} is taken as 1 for the \code{binomial} family.}

\item{se}{the standard error. \code{se} is calculated by dividing \code{sd} (standard deviation) by the square root of the sample size.}

\item{CI}{a 95\code{\%} confidence interval.}

%\item{ResSE}{Residual Standard Error, which is the square root of the dispersion.}
}

\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{MAP.estimation}} and \code{\link{bfi}} functions.
}

\examples{
#-------------
# y ~ Gaussian
#-------------

set.seed(1123581)
n       <- 30
p       <- 4
X       <- data.frame(matrix(rnorm(n * p), n, p))
b       <- 1:2
eta     <- b[1] + X[, 1] * b[2]
mu      <- gaussian()$linkinv(eta)
sigma2e <- 0.5
y       <- rnorm(n, mu, sd=sqrt(sigma2e))
lambda  <- 0.1
Lambda  <- inv.prior.cov(X, lambda=c(lambda,sigma2e), family=gaussian)
fit     <- MAP.estimation(y, X, family=gaussian, Lambda)
class(fit)

summary(fit)
sumfit  <- summary(fit, curmat = T)
sumfit$logLik
sumfit$dispersion
sumfit$CI
class(sumfit)
}
\keyword{models}
\keyword{regression}
