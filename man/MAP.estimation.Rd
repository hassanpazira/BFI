\name{MAP.estimation}
\alias{MAP.estimation}
\title{Maximum A Posteriori estimation}
\description{
\code{MAP.estimation} function is used (in local centers) to estimate Maximum A Posterior (MAP) of the parameters for the GLM and Survival models.
}
\usage{
MAP.estimation(y, X, family = gaussian, Lambda,
               intercept = TRUE, initial = NULL,
               control = list())
}
\arguments{

\item{y}{response vector. If the \code{binomial} family is used, this argument can be a vector with entries 0 (failure) or 1 (success). Alternatively, the response can be a matrix where the first column is the number of \dQuote{successes} and the second column is the number of \dQuote{failures}.}

\item{X}{design matrix of dimension \eqn{n \times p}, where \eqn{p} is the number of covariates (predictors) plus intercept.}

\item{family}{a description of the error distribution and link function used to specify the model. This can be a character string naming a family function or the result of a call to a family function (see \code{\link{family}} for details). By default the \code{gaussian} family (with identity link function) is used.}

\item{Lambda}{the matrix used as the prior of the inverse variance-covariance matrix of Gaussian distribution. It can be created by the function \code{inv.prior.cov()}.}

\item{intercept}{logical flag for fitting an intercept. The intercept is fitted if \code{intercept=TRUE} (the default) or set to zero if \code{intercept=FALSE}. Although the intercept is included in the model by default, it can be penalized.}

\item{initial}{a vector specifying initial values for the parameters (intercept, coefficients and/or error variance) to be optimized over. For the \code{gaussian} family, it should be a \eqn{p+1}-dimensional vector, and for \code{binomial} the length of the vector should be \eqn{p}, where \eqn{p} is the number of covariates when \code{intercept=FALSE} or the number of covariates plus intercept if \code{intercept=TRUE}. Since the \code{'L-BFGS-B'} method is used in the algorithm, these values should always be finite. Default is a vector of zeros.}

\item{control}{a list of control parameters. See \sQuote{Details}.}
}

\value{

\code{MAP.estimation} returns a list containing the following components:

\item{theta_hat}{the \eqn{(p + 1)}-dimensional vector corresponding to the maximum a posteriori (MAP) estimation of the parameters. If \code{intercept=FALSE}, the first element of this vector is set to zero;}

\item{A_hat}{the curvature matrix around the point \code{theta_hat}. If \code{intercept=FALSE}, the first row and column of this matrix are set to zero;}

\item{sd}{the \eqn{(p + 1)}-dimensional vector of standard deviation of estimates in \code{theta_hat}, i.e., \code{sqrt(diag(solve(A_hat)))}. If \code{intercept=FALSE}, the first element of this vector is set to zero;}

\item{Lambda}{the matrix used as the prior of the inverse variance-covariance matrix. If \code{intercept=FALSE}, the first row and column of this matrix are set to zero;}

\item{formula}{the formula of the model;}

\item{n}{sample size;}

\item{np}{the number of coefficients/regression parameters. In other words, number of predictors plus the intercept if \code{intercept=TRUE}, or without intercept if \code{intercept=FALSE};}

\item{value}{the value of the 'negative' loglikelihood function corresponding to \code{theta_hat};}

\item{family}{a description of the error distribution used in the model;}

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
\code{MAP.estimation} function implements

%Method "L-BFGS-B" is that of Byrd et. al. (1995) which allows box constraints, that is each variable can be given a lower and/or upper bound. The initial value must satisfy the constraints. This uses a limited-memory modification of the BFGS quasi-Newton method. If non-trivial bounds are supplied, this method will be selected, with a warning.


%\code{dglars} function implements the differential geometric generalization of the least angle regression method (Efron et al., 2004) proposed in Augugliaro et al. (2013) and Pazira et al. (2017).

%As in \dQuote{\code{\link{glm}}}, the user can specify family and link function using the argument \code{family}. When the \code{binomial} family is used, the responce can be a vector with entries 0/1 (failure/success) or, alternatively, a matrix where the first column is the number of \dQuote{successes} and the second column is the number of \dQuote{failures}. Starting with the version 2.0.0, the model can be specified combining family and link functions as describted in the following table:

%When gaussian, gamma or inverse.gaussian is used to model the error distribution, \code{dglars} returns the vector of the estimates of the dispersion parameter \eqn{\phi}{phi}; by default, the generalized Pearson statistic is used as estimator but the user can use the function \code{phihat} to specify other estimators (see \code{\link{phihat}} for more details).

Note that: although an intercept is included in the model by default (except for \code{family=cox}) using \code{intercept=TRUE}, it can be penalized. It can also be set to zero using \code{intercept=FALSE}.

The argument \code{control} is a list that can supply any of the following components:

\describe{

\item{\code{maxit}:}{is the maximum number of iterations. Default is {1e2};}

\item{\code{factr}:}{controls the convergence of the \code{'L-BFGS-B'} method. Convergence occurs when the reduction in the objective is within this factor of the machine tolerance. Default is {1e7}, that is a tolerance of about {1e-8};}

\item{\code{pgtol}:}{helps control the convergence of the \code{'L-BFGS-B'} method. It is a tolerance on the projected gradient in the current search direction. Default is zero, when the check is suppressed;}

\item{\code{trace}:}{is a non-negative integer. If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information: for the method \code{'L-BFGS-B'} there are six levels of tracing. To understand exactly what these do see the source code of \code{optim} function in the \pkg{\link{stats}} package;}

\item{\code{REPORT}:}{is the frequency of reports for the \code{'L-BFGS-B'} method if \code{'control$trace'} is positive. Default is every 10 iterations;}

\item{\code{lmm}:}{is an integer giving the number of \code{BFGS} updates retained in the \code{'L-BFGS-B'} method. Default is {5}.}
}
}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023). \emph{Bayesian Federated Inference for Statistical Models}, \emph{Statistics in Medicine}, Vol. 0(0), 0-0. <https://doi.org/10.48550/arXiv.2302.07677>
}

\author{
Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{bfi}} and \code{\link{summary.bfi}}.
}
\examples{
#-------------
# y ~ Gaussian
#-------------

set.seed(11235813)
n      <- 30
p      <- 3    # number of coefficients (without intercept)
X      <- data.frame(matrix(rnorm(n * p), n, p))
eta    <- 1 + 2 * X[,1]    # with an intercept
mu     <- gaussian()$linkinv(eta)
sigma2 <- 1.5
# the true theta is c(1, 2, 0, 0, sigma2)
y      <- rnorm(n, mu, sd = sqrt(sigma2))
lambda <- 0.01
# inverse of covariance matrix:
Lambda <- inv.prior.cov(X, lambda = c(lambda,sigma2), family = gaussian)

# MAP estimates of the parameters of interest (including 'intercept') and curvature matrix
(fit <- MAP.estimation(y, X, family = gaussian, Lambda))
class(fit)

# MAP estimates without 'intercept'
Lambda <- inv.prior.cov(X, lambda = c(lambda,sigma2), family = gaussian,  intercept = FALSE)
(fit1 <- MAP.estimation(y, X, family = gaussian, Lambda, intercept = FALSE))

}
\keyword{Bayesian}
\keyword{Federated}
