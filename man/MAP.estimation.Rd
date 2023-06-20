\name{MAP.estimation}
\alias{MAP.estimation}
\title{Maximum A Posteriori estimation}
\description{
\code{MAP.estimation} function is used (in local centers) to estimate Maximum A Posterior (MAP) of the parameters for the GLM and Survival models.
}
\usage{
MAP.estimation(y, X, family = gaussian, Lambda, intercept = TRUE,
               initial = NULL, control = list())
}
\arguments{

\item{y}{response vector. If the \code{binomial} family is used, this argument can be a vector with entries 0 (failure) or 1 (success). Alternatively, the response can be a matrix where the first column is the number of \dQuote{successes} and the second column is the number of \dQuote{failures}.}

\item{X}{design matrix of dimension \eqn{n \times p}, where \eqn{p} is the number of covariates (predictors) plus intercept.}

\item{family}{a description of the error distribution and link function used to specify the model. This can be a character string naming a family function or the result of a call to a family function (see \code{\link{family}} for details). By default the \code{gaussian} family (with identity link function) is used.}

\item{Lambda}{the matrix used as the prior of the inverse variance-covariance matrix of Gaussian distribution. It can be created by the function \code{inv.prior.cov()}.}

\item{intercept}{logical flag for fitting an intercept. The intercept is fitted if \code{intercept=TRUE} (the default) or set to zero if \code{intercept=FALSE}. Although the intercept is included in the model by default, it can be penalized.}

\item{initial}{a vector specifying initial values for the parameters (intercept, coefficients and/or error variance) to be optimized over. For the \code{gaussian} family, it should be a \eqn{p+1}-dimensional vector with a non-negative value for the last element, and for \code{binomial} the length of the vector should be \eqn{p}, where \eqn{p} is the number of covariates plus intercept. Since the \code{'L-BFGS-B'} method is used in the algorithm, these values should always be finite. Default is a vector of zeros.}

\item{control}{a list of control parameters. See \sQuote{Details}.}
}

\value{

\code{MAP.estimation} returns a list containing the following components:

\item{theta_hat}{the \eqn{p}- or (\eqn{p + 1})-dimensional vector corresponding to the maximum a posteriori (MAP) estimation of the parameters, where \eqn{p} is the number of covariates plus intercept. If \code{intercept=FALSE}, the first element of the vector is set to zero;}

\item{A_hat}{the curvature matrix around the point \code{theta_hat}. It's a \eqn{(p\times p)}- or \eqn{(p + 1)\times (p+1)}-dimensional matrix depending on the family used. If \code{intercept=FALSE}, the first row and column of the matrix are set to zero;}

\item{sd}{the \eqn{p}- or (\eqn{p + 1})-dimensional vector of standard deviation of estimates in \code{theta_hat}, i.e., \code{sqrt(diag(solve(A_hat)))}. If \code{intercept=FALSE}, the first element of this vector is set to zero;}

\item{Lambda}{the matrix used as the prior of the inverse variance-covariance matrix. If \code{intercept=FALSE}, the first row and column of this matrix are set to zero;}

\item{formula}{the formula of the model;}

\item{n}{sample size;}

\item{np}{the number of coefficients/regression parameters. In other words, number of predictors plus the intercept if \code{intercept=TRUE}, or without intercept if \code{intercept=FALSE};}

\item{value}{the value of the 'negative' log-likelihood function corresponding to \code{theta_hat};}

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
\code{MAP.estimation} function finds the Maximum A Posteriori (MAP) estimates of model parameters, i.e., the values associated with the peak of the log-posterior distribution (the posterior mode).
In other words, \code{MAP.estimation()} optimizes the log-posterior density with respect to the parameter vector to obtain its MAP estimation.
In addition to the model parameters, i.e., coefficients (\eqn{{\beta}}'s) and variance error (\eqn{\sigma^2_e}), the curvature matrix (Hessian of the log-posterior) is estimated around the mode.

The current version of the \pkg{{BFI}} package supports two families of GLM: \code{gaussian} and \code{binomial}. For the \code{gaussian} family, \code{MAP.estimation} function returns the (\eqn{p + 1})-dimensional vector of the MAP estimates of the coefficients \eqn{(\beta_0, \beta_1,\ldots, \beta_{p-1})} and variance error parameter, respectively.
When the \code{binomial} family is used, it returns the MAP estimates of the \eqn{p} coefficients while the variance error or dispersion is set to one.

Note that: although an intercept is included in the model by default (except for \code{family=cox}) using \code{intercept=TRUE}, it can be penalized. It can also be set to zero using \code{intercept=FALSE}.

Since all elements of the lists in the first three arguments of the \code{bfi} function must have equal dimensions, the function \code{MAP.estimation} creates the output for all centers with similar dimensions even if in some centers the intercept is set to zero \code{intercept=FALSE} (of course for the same \code{family}'s) by setting the corresponding elements to zero.

The \code{MAP.estimation} function returns an object of class `\code{bfi}`. Therefore, \code{summary()} can be used for the object returned by \code{MAP.estimation()}.


To solve unconstrained and bound-constrained optimization problems, the \code{MAP.estimation} function utilizes an optimization algorithm called Limited-memory Broyden-Fletcher-Goldfarb-Shanno with Bound Constraints ({L-BFGS-B}), Byrd et. al. (1995).
The L-BFGS-B algorithm is a limited-memory \dQuote{quasi-Newton} method that iteratively updates the parameter estimates by approximating the inverse Hessian matrix using gradient information from the history of previous iterations. This approach allows the algorithm to approximate the curvature of the posterior distribution and efficiently search for the optimal solution, which makes it computationally efficient for problems with a large number of variables.

By default, the algorithm uses a relative change in the objective function as the convergence criterion. When the change in the objective function between iterations falls below a certain threshold (`\code{factr}`) the algorithm is considered to have converged.

In Jonker et. al. (2023), the Gaussian and non-informative priors were used on the parameters.
In this case, the posterior distribution is highly expected to be unimodal for a Gaussian GLM, which means it has only one minimum and no local optima. This property ensures that the optimization algorithm will converge to the global minimum, guaranteeing convergence.
After fitting the model, the convergence status can be examined by checking the \code{convergence} attribute of the model object using the \code{$convergence} attribute. See \sQuote{Value} for the possible options of this attribute.


If a logistic regression model (\code{family=binomial}) is used, it's generally more likely that the posterior distribution will be unimodal with these two priors, especially if the data provides strong information about the parameters or sample size is large.
Since the starting point has a great impact on the converged point especially for densities which are not unimodal, multiple starting points can be used in the \code{initial} argument to check for unimodality. Indeed, the negative log-posterior density is minimized from an chosen starting point \code{initial}.

However, there could be cases with small data sets (non-informative data), high-dimensionality, or inappropriate optimization settings, which might lead to multimodality or other complex behavior in the posterior distribution, resulting in users experiencing convergence issues. In such cases, it may be necessary to investigate and adjust optimization parameters to facilitate convergence. It can be done using the \code{control} argument.

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
Jonker M.A., Pazira H. and Coolen A.C.C. (2023). \emph{Bayesian Federated Inference for Statistical Models}. {Statistics in Medicine}, Vol. 0(0), 0-0. <https://doi.org/10.48550/arXiv.2302.07677>

Byrd R.H., Lu P., Nocedal J. and Zhu C. (1995). \emph{A limited memory algorithm for bound constrained optimization}. {SIAM Journal on Scientific Computing}, 16, 1190-1208. <https://doi.org/10.1137/0916069>
}

\author{
Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{bfi}} and \code{\link{summary.bfi}}
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
