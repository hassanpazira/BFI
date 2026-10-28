\name{summary.bfi}
\alias{summary}
\alias{summary.bfi}
\alias{print.summary.bfi}
\title{Summarizing BFI fits}

\description{Summary method for an object with class 'bfi' created by the \code{MAP.estimation} and \code{bfi} functions.}

\usage{
\method{summary}{bfi}(object,
        cur_mat = FALSE,
        digits = max(3, getOption("digits") - 3),
        ...)

\method{print}{summary.bfi}(x,
        digits = x$summary_digits,
        cur_mat = x$summary_cur_mat,
        ...)
}

\arguments{
\item{object}{A fitted object of class \code{"bfi"}.}

\item{x}{An object of class \code{"summary.bfi"} returned by \code{summary.bfi()}.}

\item{cur_mat}{Logical; if \code{TRUE}, the minus curvature matrix is also printed when the summary object is displayed.
The matrix remains available in \code{A_hat} regardless of this option.
Default is \code{FALSE} for \code{summary.bfi()}.}

\item{digits}{Number of digits used to format the printed results.}

\item{...}{Additional arguments, currently unused.}
}

\details{
\code{summary.bfi()} prepares a structured summary of results
returned by \code{MAP.estimation()} or \code{bfi()}.

The resulting object can be stored without producing printed output.
When displayed, \code{print.summary.bfi()} prints the model family,
link function or baseline hazard, parameter estimates, posterior
standard deviations, and 95\% credible intervals.

For local models, the log-posterior value and convergence status
are also displayed.

The minus curvature matrix is retained in \code{A_hat}.
When \code{cur_mat = TRUE} is supplied to \code{summary.bfi()},
this matrix is also printed when the summary object is displayed.

The print method returns the summary object invisibly.
}

\value{

\code{summary.bfi} returns an object of class \code{"summary.bfi"}, containing the components described below.
\code{print.summary.bfi()} displays the summary and returns
the supplied object invisibly.

\item{theta_hat}{the component from \code{object}. The last element of this vector is the estimate of the dispersion parameter (sigma2) if \code{family = "gaussian"}. See the \code{\link{MAP.estimation}} and \code{\link{bfi}} functions.}

\item{A_hat}{the component from \code{object}. See the \code{\link{MAP.estimation}} and \code{\link{bfi}} functions.}

\item{sd}{the posterior standard deviations corresponding to the reported parameter estimates. For Gaussian models with center-specific residual variances, the standard deviations corresponding to \code{sigma2} are transformed from the \eqn{\log(\sigma^2)} scale to the original \eqn{\sigma^2} scale using the delta method. See the \code{\link{MAP.estimation}} and \code{\link{bfi}} functions.}

\item{Lambda}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{formula}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{n}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{np}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{family}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{intercept}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{convergence}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{control}{the component from \code{object}. See the \code{\link{MAP.estimation}} function.}

\item{stratified}{the component from \code{object}. See the \code{\link{bfi}} function.}

\item{estimate}{the estimated regression coefficients. For \code{gaussian} models the common dispersion parameter is excluded, but center-specific residual variances obtained with \code{strat_par = 2} or \code{strat_par = c(1, 2)} are included.}

\item{logLikPost}{the value of the log-likelihood posterior density evaluated at estimates (\code{theta_hat}).}

\item{link}{the link function only for GLMs, not for the survival family. By default the \code{gaussian} family with \code{identity} link function and the \code{binomial} family with \code{logit} link function are used.}

\item{dispersion}{the estimated residual variance \code{sigma2} for Gaussian models with a common residual variance. For Gaussian models with center-specific residual variances, this component is \code{NULL}. For the \code{binomial} family, the dispersion parameter is fixed at \code{1}. For the \code{survival} family, this component is \code{NULL}.}

\item{CI}{a 95\code{\%} credible interval for the parameters. For center-specific residual variances the interval is computed on the \eqn{\log(\sigma^2)} scale and back-transformed, so it is always positive and not symmetric around the estimate. The standard deviation printed for these parameters is obtained by the delta method on the \eqn{\sigma^2} scale.}

\item{coef_table}{A matrix containing the reported parameter
estimates, posterior standard deviations, and 95\% credible
interval limits.}

\item{summary_digits}{The number of digits selected for
printing the summary.}

\item{summary_cur_mat}{A logical value indicating whether the
curvature matrix should be displayed when the summary is printed.}

}

\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{h.pazira@arq.org}
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

