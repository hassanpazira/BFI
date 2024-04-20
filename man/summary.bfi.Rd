\name{summary.bfi}
\alias{summary}
\alias{summary.bfi}

\title{Summarizing BFI Fits}

\description{Summary method for an object with class 'bfi' created by the \code{MAP.estimation} and \code{bfi} functions.}

\usage{
\method{summary}{bfi}(object,
        cur_mat = FALSE,
        digits = max(3, getOption("digits") - 3),
        ...)
}

\arguments{
\item{object}{fitted \code{bfi} object.}

\item{cur_mat}{logical; if \code{TRUE}, minus the curvature matrix around the estimated parameters is returned and printed. Default is \code{FALSE}.}

\item{digits}{significant digits in printout.}

\item{...}{additional arguments affecting the summary produced.}
}

\details{
\code{summary.bfi()} gives information about the MAP estimates of parameters of the model. It can be used for the \code{bfi} objects built by the \code{MAP.estimation} and \code{bfi} functions.

The output of the summary method shows the details of the model, i.e. formula, family and link function used to specify the generalized linear model, followed by information about the estimates, standard deviations and credible intervals.
Information about the log-likelihood posterior and convergence status are also provided.

By default, \code{summary.bfi} function does not return (minus) the curvature matrix, but the user can use \code{cur_mat = TRUE} to print it.
}

\value{

\code{summary.bfi} returns an object of class \code{summary.bfi}, a list with the following components:

\item{theta_hat}{the component from \code{object}. The last element of this vector is the estimate of the dispersion parameter (sigma2) if \code{family = "gaussian"}. See the \code{\link{MAP.estimation}} and \code{\link{bfi}} functions.}

\item{A_hat}{the component from \code{object}. See the \code{\link{MAP.estimation}} and \code{\link{bfi}} functions.}

\item{sd}{the component from \code{object}. If \code{family = "gaussian"}, the last element of this vector is the square root of the estimated dispersion. See the \code{\link{MAP.estimation}} and \code{\link{bfi}} functions.}

\item{Lambda}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{formula}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{n}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{np}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{family}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{intercept}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{convergence}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{control}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{stratified}{the component from \code{object}. See the \code{\link{bfi}} function.}

\item{estimate}{the estimated regression coefficients, i.e., without the estimate \code{sigma2}.}

\item{logLikPost}{the value of the log-likelihood posterior density evaluated at estimates (\code{theta_hat}).}

\item{link}{the link function only for GLMs, not for the survival family. By default the \code{gaussian} family with \code{identity} link function and the \code{binomial} family with \code{logit} link function are used.}

\item{dispersion}{the estimated variance of the random error, i.e., \code{sigma2}. The \code{dispersion} is taken as \code{1} for the \code{binomial} family.}

\item{CI}{a 95\code{\%} credible interval of the MAP estimates of the parameters.}
}

\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{MAP.estimation}} and \code{\link{bfi}}
}

\examples{
#-------------
# y ~ Gaussian
#-------------
# model assumption:
theta <- c(1, 2, 3, 4, 1.5)  # coefficients and sigma2 = 1.5

#----------------
# Data Simulation
#----------------
n      <- 40
X      <- data.frame(x1=rnorm(n),                     # continuous variable
                     x2=sample(1:3, n, replace=TRUE)) # categorical variable
Xx2_1  <- ifelse(X$x2 == '2', 1, 0)
Xx2_2  <- ifelse(X$x2 == '3', 1, 0)
X$x2   <- as.factor(X$x2)
eta    <- theta[1] + theta[2] * X$x1 + theta[3] * Xx2_1 + theta[4] * Xx2_2
mu     <- gaussian()$linkinv(eta)
y      <- rnorm(n, mu, sd = sqrt(theta[5]))

#----------------
# MAP estimations
#----------------
Lambda <- inv.prior.cov(X, lambda = c(0.1, 0.5), family = "gaussian")
fit    <- MAP.estimation(y, X, family = "gaussian", Lambda)
class(fit)

#-------------------------
# Summary of MAP estimates
#-------------------------
summary(fit)
sumfit <- summary(fit, cur_mat = TRUE)
sumfit$estimate
sumfit$logLikPost
sumfit$dispersion
sumfit$CI
class(sumfit)

}
\keyword{regression}
\keyword{models}
\keyword{survival}

