\name{MAP.estimation}
\alias{MAP.estimation}
\title{Maximum A Posteriori estimation}
\description{
\code{MAP.estimation} function is used (in local centers) to compute Maximum A Posterior (MAP) estimators of the parameters for the GLM and soon for Survival models.
}
\usage{
MAP.estimation(y, X, family = gaussian, Lambda, intercept = TRUE,
               initial = NULL, control = list())
}
\arguments{

\item{y}{response vector. If the \code{binomial} family is used, this argument is a vector with entries 0 (failure) or 1 (success). Alternatively, the response can be a matrix where the first column is the number of \dQuote{successes} and the second column is the number of \dQuote{failures}.}

\item{X}{design matrix of dimension \eqn{n \times p}, where \eqn{p} is the number of covariables or predictors.}

\item{family}{a description of the error distribution and link function used to specify the model. This can be a character string naming a family function or the result of a call to a family function (see \code{\link{family}} for details). In the current version of the package, the family of model can be \code{gaussian} (with \code{identity} link function) and \code{binomial} (with \code{logit} link function). By default the \code{gaussian} family is used. In case of a linear regression model, \code{family = gaussian}, there is an extra model parameter for the variance of measurement error.}

\item{Lambda}{the inverse variance-covariance matrix of the Gaussian distribution that is used as prior distribution for the model parameters. The dimension of the matrix depends on the number of columns of \code{X}, type of the covariates (continuous / dichotomous or categorical), \code{family}, and \code{intercept}. However, \code{Lambda} can be easily created by \code{inv.prior.cov()}.}

\item{intercept}{logical flag for fitting an intercept. If \code{intercept=TRUE} (the default), the intercept is fitted, i.e., it is included in the model, and if \code{intercept=FALSE} it is set to zero, i.e., it's not in the model.}

\item{initial}{a vector specifying initial values for the parameters to be optimized over. The length of \code{initial} is equal to the number of model parameters and thus, is equal to the number of rows or columns of \code{Lambda}.
Since the \code{'L-BFGS-B'} method is used in the algorithm, these values should always be finite. Default is a vector of zeros.}

\item{control}{a list of control parameters. See \sQuote{Details}.}
}

\value{

\code{MAP.estimation} returns a list containing the following components:

\item{theta_hat}{the vector corresponding to the maximum a posteriori (MAP) estimates of the parameters;}

\item{A_hat}{minus the curvature (or Hessian) matrix around the point \code{theta_hat}. The dimension of the matrix is the same as the argument \code{Lambda};}

\item{sd}{the vector of standard deviation of the MAP estimates in \code{theta_hat}, that is \code{sqrt(diag(solve(A_hat)))};}

\item{Lambda}{the inverse variance-covariance matrix of the Gaussian distribution that is used as prior distribution for the parameters. It's exactly the same as the argument \code{Lambda};}

\item{formula}{the formula applied;}

\item{names}{the names of the model parameters;}

\item{n}{sample size;}

\item{np}{the number of coefficients;}

\item{value}{the value of minus the log-likelihood posterior density evaluated at \code{theta_hat};}

\item{family}{the \code{\link{family}} object used.;}

\item{intercept}{logical flag used to fit an intercept if \code{TRUE}, or set to zero if \code{FALSE};}

\item{convergence}{an integer value used to encode the warnings and the errors related
to the algorithm used to fit the model. The values returned are:
\describe{
		\item{0}{algorithm has converged;}

		\item{1}{maximum number of iterations ('\code{maxit}') has been reached;}

		\item{2}{Warning from the 'L-BFGS-B' method. See the message after this value;}
	}
}

\item{control}{the list of control parameters used to compute the MAP estimates.}
}

\details{
\code{MAP.estimation} function finds the Maximum A Posteriori (MAP) estimates of the model parameters by maximizing the log-posterior density with respect to the parameters, i.e., the estimates equal the values for which the log-posterior density is maximal (the posterior mode).
In other words, \code{MAP.estimation()} optimizes the log-posterior density with respect to the parameter vector to obtain its MAP estimates.
In addition to the model parameters, i.e., coefficients (\eqn{{\beta}}'s) and variance error (\eqn{\sigma^2_e}), the curvature matrix (Hessian of the log-posterior) is estimated around the mode.

The \code{MAP.estimation} function returns an object of class `\code{bfi}`. Therefore, \code{summary()} can be used for the object returned by \code{MAP.estimation()}.

To solve unconstrained and bound-constrained optimization problems, the \code{MAP.estimation} function utilizes an optimization algorithm called Limited-memory Broyden-Fletcher-Goldfarb-Shanno with Bound Constraints ({L-BFGS-B}), Byrd et. al. (1995).
The L-BFGS-B algorithm is a limited-memory \dQuote{quasi-Newton} method that iteratively updates the parameter estimates by approximating the inverse Hessian matrix using gradient information from the history of previous iterations. This approach allows the algorithm to approximate the curvature of the posterior distribution and efficiently search for the optimal solution, which makes it computationally efficient for problems with a large number of variables.

By default, the algorithm uses a relative change in the objective function as the convergence criterion. When the change in the objective function between iterations falls below a certain threshold (`\code{factr}`) the algorithm is considered to have converged.
The convergence can be checked with the argument \code{convergence} in the output. See \sQuote{Value}.

In case of convergence issue, it may be necessary to investigate and adjust optimization parameters to facilitate convergence. It can be done using the \code{initial} and \code{control} arguments. By the argument \code{initial} the initial points of the interative optimization algorithm can be changed, and the argument \code{control} is a list that can supply any of the following components:

\describe{

\item{\code{maxit}:}{is the maximum number of iterations. Default is {100};}

\item{\code{factr}:}{controls the convergence of the \code{'L-BFGS-B'} method. Convergence occurs when the reduction in the objective is within this factor of the machine tolerance. Default for \code{factr} is {1e7}, which gives a tolerance of about {1e-9}. The exact tolerance can be checked by multiplying this value by \code{.Machine$double.eps};}

\item{\code{pgtol}:}{helps to control the convergence of the \code{'L-BFGS-B'} method. It is a tolerance on the projected gradient in the current search direction, i.e., the iteration will stop when the maximum component of the projected gradient is less than or equal to \code{pgtol}, where \code{pgtol}\eqn{\geq 0}. Default is zero, when the check is suppressed;}

\item{\code{trace}:}{is a non-negative integer. If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information: for the method \code{'L-BFGS-B'} there are six levels of tracing. To understand exactly what these do see the source code of \code{optim} function in the \pkg{\link{stats}} package;}

\item{\code{REPORT}:}{is the frequency of reports for the \code{'L-BFGS-B'} method if \code{'control$trace'} is positive. Default is every 10 iterations;}

\item{\code{lmm}:}{is an integer giving the number of \code{BFGS} updates retained in the \code{'L-BFGS-B'} method. Default is {5}.}
}
}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023). \emph{Bayesian Federated Inference for Statistical Models}. {Statistics in Medicine}, Vol. 0(0), 0-0. <https://doi.org/10.48550/arXiv.2302.07677>

Byrd R.H., Lu P., Nocedal J. and Zhu C. (1995). \emph{A limited memory algorithm for bound constrained optimization}. {SIAM Journal on Scientific Computing}, 16, 1190-1208. <https://doi.org/10.1137/0916069>
}

\author{
Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{bfi}}, \code{\link{inv.prior.cov}} and \code{\link{summary.bfi}}
}
\examples{
#-------------
# y ~ Gaussian
#-------------
# model assumption:
theta <- c(1, 2, 2, 2, 1.5)  # model parameters: coefficients and sigma2 = 1.5

#----------------
# Data Simulation
#----------------
n   <- 30   # sample size
p   <- 3    # number of coefficients without intercept
X   <- data.frame(matrix(rnorm(n * p), n, p)) # continuous variables
# linear predictor:
eta <- theta[1] + theta[2] * X$X1 + theta[3] * X$X2 + theta[4] * X$X3
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu  <- gaussian()$linkinv(eta)
y   <- rnorm(n, mu, sd = sqrt(theta[5]))

#-----------------------------------------------
# MAP estimations for theta and curvature matrix
#-----------------------------------------------
# MAP estimates with 'intercept'
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = gaussian)
(fit <- MAP.estimation(y, X, family = gaussian, Lambda))
class(fit)
summary(fit, cur_mat = TRUE)

# MAP estimates without 'intercept'
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = gaussian, intercept = FALSE)
(fit1 <- MAP.estimation(y, X, family = gaussian, Lambda, intercept = FALSE))
summary(fit1, cur_mat = TRUE)

}
\keyword{Bayesian}
\keyword{Federated}
