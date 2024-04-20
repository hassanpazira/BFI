\name{MAP.estimation}
\alias{MAP.estimation}
\title{Maximum A Posteriori estimation}
\description{
\code{MAP.estimation} function is used (in local centers) to compute Maximum A Posterior (MAP) estimators of the parameters for Generalized Linear Models (GLM) and Survival models.
}
\usage{
MAP.estimation(y,
               X,
               family = c("gaussian", "binomial", "survival"),
               Lambda,
               intercept = TRUE,
               initial = NULL,
               basehaz = c("weibul", "exp", "gomp", "poly", "pwexp"),
               alpha = 0.1,
               max_order = 2,
               n_intervals = 4,
               min_max_times,
               center_zero_sample = FALSE,
               zero_sample_cov,
               refer_cat,
               zero_cat,
               control = list())
}
\arguments{

\item{y}{response vector. If the \dQuote{\code{binomial}} family is used, this argument is a vector with entries 0 (failure) or 1 (success). Alternatively, for this family, the response can be a matrix where the first column is the number of \dQuote{successes} and the second column is the number of \dQuote{failures}. For the \dQuote{\code{survival}} family, the response is a matrix where the first column is the survival time, named \dQuote{time}, and the second column is the censoring indicator, named \dQuote{status}, with 0 indicating censoring time and 1 indicating event time.}

\item{X}{design matrix of dimension \eqn{n \times p}, where \eqn{p} is the number of covariates or predictors. Note that the order of the covariates must be the same across the centers; otherwise, the output estimates of bfi() will be incorrect.}

\item{family}{a description of the error distribution. This is a character string naming a family of the model. In the current version of the package, the family of model can be \code{"gaussian"} (with \code{identity} link function), \code{"binomial"} (with \code{logit} link function), or \code{"survival"}. Can be abbreviated. By default the \code{gaussian} family is used. In case of a linear regression model, \code{family = "gaussian"}, there is an extra model parameter for the variance of measurement error. While in the case of survival model, \code{family = "survival"}, the number of the model parameters depend on the choice of baseline hazard functions, see \sQuote{Details} for more information.}

\item{Lambda}{the inverse variance-covariance matrix of the Gaussian distribution that is used as prior distribution for the model parameters. The dimension of the matrix depends on the number of columns of \code{X}, type of the covariates (continuous / dichotomous or categorical), \code{family}, and whether an \code{intercept} is included (if applicable). However, \code{Lambda} can be easily created by \code{inv.prior.cov()}. See \code{\link{inv.prior.cov}} for more information.}

\item{intercept}{logical flag for fitting an intercept. If \code{intercept=TRUE} (the default), the intercept is fitted, i.e., it is included in the model, and if \code{intercept=FALSE} it is set to zero, i.e., it's not in the model. This argument is not used if \code{family = "survival"}.}

\item{initial}{a vector specifying initial values for the parameters to be optimized over. The length of \code{initial} is equal to the number of model parameters and thus, is equal to the number of rows or columns of \code{Lambda}. Since the \code{'L-BFGS-B'} method is used in the algorithm, these values should always be finite. Default is a vector of zeros, except for the \code{survival} family with the \code{poly} function, where it is a vector with the first \eqn{p} elements as zeros for coefficients (\eqn{\boldsymbol{\beta}}) and -0.5 for the remaining parameters (\eqn{\boldsymbol{\omega}}). For the \code{gaussian} family, the last element of the \code{initial} vector could also be considered negative, because the Gaussian prior was applied to \eqn{log(\sigma^2)}.}

\item{basehaz}{a character string representing one of the available baseline hazard functions; \code{exponential} (\code{"exp"}), \code{Weibull} (\code{"weibul"}, the default), \code{Gompertz} (\code{"gomp"}), \code{exponentiated polynomial} (\code{"poly"}), and \code{piecewise constant exponential} (\code{"pwexp"}). Can be abbreviated. It is only used when \code{family = "survival"}. If local sample size is large and the shape of the baseline hazard function is completely unknown, the ``exponentiated polynomial'' and ``piecewise exponential'' hazard functions would be preferred above the lower dimensional alternatives. However, if the local samples size is low, one should be careful using the ``piecewise exponential'' hazard function with many intervals.}

\item{alpha}{a significance level used in the chi-squared distribution (with one degree of freedom and 1-\eqn{\alpha} representing the upper quantile) to conduct a likelihood ratio test for obtaining the order of the exponentiated polynomial baseline hazard function. It is only used when \code{family = "survival"} and \code{basehaz = "poly"}. Default is 0.1. See \sQuote{Details}.}

\item{max_order}{an integer representing the maximum value of \code{q_l}, which is the order/degree minus 1 of the exponentiated polynomial baseline hazard function. This argument is only used when \code{family = "survival"} and \code{basehaz = "poly"}. Default is 2.}

\item{n_intervals}{an integer representing the number of intervals in the piecewise exponential baseline hazard function. This argument is only used when \code{family = "survival"} and \code{basehaz = "pwexp"}. Default is 4.}

\item{min_max_times}{a scalar representing the minimum of the maximum event times observed in the centers. The value of this argument should be defined by the central server (which has access to the maximum event times of all the centers) and is only used when \code{family = "survival"} and \code{basehaz = "pwexp"}.}

\item{center_zero_sample}{logical flag indicating whether the center has a categorical covariate with no observations/individuals in one of the categories. Default is \code{center_zero_sample = FALSE}.}

\item{zero_sample_cov}{either a character string or an integer representing the categorical covariate that has no samples/observations in one of its categories. This covariate should have at least two categories, one of which is the reference. It is used when \code{center_zero_sample = TRUE}.}

\item{refer_cat}{a character string representing the reference category. The category with no observations (the argument \code{zero_cat}) cannot be used as the reference in the argument \code{refer_cat}. It is used when \code{center_zero_sample = TRUE}.}

\item{zero_cat}{a character string representing the category with no samples/observations. It is used when \code{center_zero_sample = TRUE}.}

\item{control}{a list of control parameters. See \sQuote{Details}.}
}

\value{

\code{MAP.estimation} returns a list containing the following components:

\item{theta_hat}{the vector corresponding to the maximum a posteriori (MAP) estimates of the parameters. For the \code{gaussian} family, the last element of this vector is \eqn{\sigma^2};}

\item{A_hat}{minus the curvature (or Hessian) matrix around the point \code{theta_hat}. The dimension of the matrix is the same as the argument \code{Lambda};}

\item{sd}{the vector of standard deviation of the MAP estimates in \code{theta_hat}, that is \code{sqrt(diag(solve(A_hat)))};}

\item{Lambda}{the inverse variance-covariance matrix of the Gaussian distribution that is used as prior distribution for the parameters. It's exactly the same as the argument \code{Lambda};}

\item{formula}{the formula applied;}

\item{names}{the names of the model parameters;}

\item{n}{sample size;}

\item{np}{the number of coefficients;}

\item{q_l}{the order/degree minus 1 of the exponentiated polynomial baseline hazard function determined for the current center by the likelihood ratio test. This output argument, \code{q_l}, is only shown when \code{family = "survival"} and \code{basehaz = "poly"}, and will be used in the function \code{bfi()};}

\item{theta_A_poly}{an array where the first component is a matrix with columns representing the MAP estimates of the parameters for different \code{q_l}'s, i.e., \code{q_l}, \code{q_l}+1, ..., \code{max_order}. The other components are minus the curvature matrices for different \code{q_l}'s, i.e., \code{q_l}, \code{q_l}+1, ..., \code{max_order}. Therefore, the first non-NA curvature matrix is equal to the output argument \code{A_hat}. This output argument, \code{theta_A_poly}, is only shown if \code{family = "survival"} and \code{basehaz = "poly"}, and will be used in the function \code{bfi()};}

\item{lev_no_ref_zer}{a vector containing the names of the levels of the categorical covariate that has no samples/observations in one of its categories. The name of the category with no samples and the name of the reference category are excluded from this vector. This argument is shown when \code{family = "survival"} and \code{basehaz = "poly"}, and will be used in the function \code{bfi()};}

\item{zero_sample_cov}{the categorical covariate that has no samples/observations in one of its categories. It is shown when \code{center_zero_sample = TRUE}, and can be used in the function \code{bfi()};}

\item{refer_cat}{the reference category. It is shown when \code{center_zero_sample = TRUE}, and can be used in the function \code{bfi()};}

\item{zero_cat}{the category with no samples/observations. It is shown when \code{center_zero_sample = TRUE}, and can be used in the function \code{bfi()};}

\item{value}{the value of minus the log-likelihood posterior density evaluated at \code{theta_hat};}

\item{family}{the \code{family} used;}

\item{basehaz}{the baseline hazard function used;}

\item{intercept}{logical flag used to fit an intercept if \code{TRUE}, or set to zero if \code{FALSE};}

\item{convergence}{an integer value used to encode the warnings and the errors related to the algorithm used to fit the model. The values returned are:
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
In addition to the model parameters (i.e., coefficients (\eqn{\boldsymbol{\beta}}) and variance error (\eqn{\sigma^2_e}) for \code{gaussian} or the parameters of the baseline hazard (\eqn{\boldsymbol{\omega}}) for \code{survival}), the curvature matrix (Hessian of the log-posterior) is estimated around the mode.

The \code{MAP.estimation} function returns an object of class `\code{bfi}`. Therefore, \code{summary()} can be used for the object returned by \code{MAP.estimation()}.

For the case where \code{family = "survival"} and \code{basehaz = "poly"}, we assume that in all centers the \eqn{q_\ell}'s are equal.
However, the order of the estimated polynomials may vary across the centers so that each center can have different number of parameters, say \eqn{q_\ell}+1.
%Thus, in each location, the vector parameter is a (\eqn{ p+q_\ell} +1)-dimensional vector \eqn{\boldsymbol{\theta}_\ell^\top=(\boldsymbol{\beta}^\top, \boldsymbol{\omega}_\ell^\top)=(\beta_1, \ldots,\beta_p, \omega_{0}, \omega_{1}, \ldots, \omega_{q_\ell})}.
After obtaining the estimates within the local centers (by using \code{MAP.estimation()}) and having all estimates in the central server, we choose the order of the polynomial approximation for the combined data to be the maximum of the orders of the local polynomial functions, i.e., \eqn{\max \{q_1, \ldots, q_L \}}, to approximate the global baseline hazard (exponentiated polynomial) function more accurately. This is because the higher-order polynomial approximation can capture more complex features and details in the combined data. Using the higher-order approximation ensures that we account for the higher-order moments and features present in the data while maintaining accuracy.
As a result, all potential cases are stored in the \code{theta_A_poly} argument to be used in \code{bfi()} by the central server.
For further information on the \code{survival} family, refer to the 'References' section.

To solve unconstrained and bound-constrained optimization problems, the \code{MAP.estimation} function utilizes an optimization algorithm called Limited-memory Broyden-Fletcher-Goldfarb-Shanno with Bound Constraints (L-BFGS-B), Byrd et. al. (1995).
The L-BFGS-B algorithm is a limited-memory \dQuote{quasi-Newton} method that iteratively updates the parameter estimates by approximating the inverse Hessian matrix using gradient information from the history of previous iterations. This approach allows the algorithm to approximate the curvature of the posterior distribution and efficiently search for the optimal solution, which makes it computationally efficient for problems with a large number of variables.

By default, the algorithm uses a relative change in the objective function as the convergence criterion. When the change in the objective function between iterations falls below a certain threshold (`\code{factr}`) the algorithm is considered to have converged.
The convergence can be checked with the argument \code{convergence} in the output. See \sQuote{Value}.

In case of convergence issue, it may be necessary to investigate and adjust optimization parameters to facilitate convergence. It can be done using the \code{initial} and \code{control} arguments. By the argument \code{initial} the initial points of the interative optimization algorithm can be changed, and the argument \code{control} is a list that can supply any of the following components:

\describe{

\item{\code{maxit}:}{is the maximum number of iterations. Default is 150;}

\item{\code{factr}:}{controls the convergence of the \code{'L-BFGS-B'} method. Convergence occurs when the reduction in the objective is within this factor of the machine tolerance. Default for \code{factr} is 1e7, which gives a tolerance of about 1e-9. The exact tolerance can be checked by multiplying this value by \code{.Machine$double.eps};}

\item{\code{pgtol}:}{helps to control the convergence of the \code{'L-BFGS-B'} method. It is a tolerance on the projected gradient in the current search direction, i.e., the iteration will stop when the maximum component of the projected gradient is less than or equal to \code{pgtol}, where \code{pgtol}\eqn{\geq 0}. Default is zero, when the check is suppressed;}

\item{\code{trace}:}{is a non-negative integer. If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information: for the method \code{'L-BFGS-B'} there are six levels of tracing. To understand exactly what these do see the source code of \code{optim} function in the \pkg{\link{stats}} package;}

\item{\code{REPORT}:}{is the frequency of reports for the \code{'L-BFGS-B'} method if \code{'control$trace'} is positive. Default is every 10 iterations;}

\item{\code{lmm}:}{is an integer giving the number of \code{BFGS} updates retained in the \code{'L-BFGS-B'} method. Default is 5.}
}

}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2024). \emph{Bayesian federated inference for estimating statistical models based on non-shared multicenter data sets}, \emph{Statistics in Medicine}, 43(12): 2421-2438. <https://doi.org/10.1002/sim.10072>

Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker M.A. (2024). \emph{Bayesian Federated Inference for Survival Models}, \emph{arXiv}. <https://arxiv.org/abs/2404.17464>

Jonker M.A., Pazira H. and Coolen A.C.C. (2024b). \emph{Bayesian Federated Inference for regression models with heterogeneous multi-center populations}, \emph{arXiv}. <https://arxiv.org/abs/2402.02898>

Byrd R.H., Lu P., Nocedal J. and Zhu C. (1995). \emph{A limited memory algorithm for bound constrained optimization}. {SIAM Journal on Scientific Computing}, 16, 1190-1208. <https://doi.org/10.1137/0916069>
}

\author{
Hassan Pazira and Marianne Jonker \cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{bfi}}, \code{\link{inv.prior.cov}} and \code{\link{summary.bfi}}
}
\examples{

###--------------###
### y ~ Gaussian ###
###--------------###

# Setting a seed for reproducibility
set.seed(11235)

# model parameters: coefficients and sigma2 = 1.5
theta <- c(1, 2, 2, 2, 1.5)

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

# Load the BFI package
library(BFI)

#-----------------------------------------------
# MAP estimations for theta and curvature matrix
#-----------------------------------------------
# MAP estimates with 'intercept'
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = "gaussian")
(fit <- MAP.estimation(y, X, family = "gaussian", Lambda))
class(fit)
summary(fit, cur_mat = TRUE)

# MAP estimates without 'intercept'
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = 'gaussian', intercept = FALSE)
(fit1 <- MAP.estimation(y, X, family = 'gaussian', Lambda, intercept = FALSE))
summary(fit1, cur_mat = TRUE)


###-----------------###
### Survival family ###
###-----------------###

# Setting a seed for reproducibility
set.seed(112358)

#-------------------------
# Simulating Survival data
#-------------------------
n    <- 40
beta <- 1:4
p    <- length(beta)
X    <- data.frame(matrix(rnorm(n * p), n, p)) # continuous (normal) variables

## Simulating survival data from Weibull with a predefined censoring rate of 0.3
y <- surv.simulate(Z = list(X), beta = beta, a = 5, b = exp(1.8), u1 = 0.1,
                   cen_rate = 0.3, gen_data_from = "weibul")$D[[1]][, 1:2]

##
## MAP estimations with "weibul" function
##
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = 'survival', basehaz = "weibul")
fit2 <- MAP.estimation(y, X, family = 'survival', Lambda = Lambda, basehaz = "weibul")
fit2
fit2$theta_hat

##
## MAP estimations with "poly" function
##
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = 'survival', basehaz = 'poly')
fit3 <- MAP.estimation(y, X, family = "survival", Lambda = Lambda, basehaz = "poly")
# Degree of the exponentiated polynomial baseline hazard
fit3$q_l + 1
# theta_hat for (beta_1, ..., beta_p, omega_0, ..., omega_{q_l})
fit3$theta_A_poly[,,1][,fit3$q_l+1] # equal to fit3$theta_hat
# A_hat
fit3$theta_A_poly[,,fit3$q_l+2] # equal to fit3$A_hat
summary(fit3, cur_mat = TRUE)

##
## MAP estimations with "pwexp" function with 3 intervals
##
# Assume we have 4 centers
Lambda <- inv.prior.cov(X, lambda = c(0.1, 1), family = 'survival',
                        basehaz = 'pwexp', n_intervals = 3)
# For this baseline hazard ("pwexp"), we need to know
# maximum survival times of the 3 other centers:
max_times <- c(max(rexp(30)), max(rexp(50)), max(rexp(70)))
# Minimum of the maximum values of the survival times of all 4 centers is:
min_max_times <- min(max(y$time), max_times)
fit4 <- MAP.estimation(y, X, family = "survival", Lambda = Lambda, basehaz = "pwexp",
                       n_intervals = 3, min_max_times=max(y$time))
summary(fit4, cur_mat = TRUE)

}
\keyword{optimize}
\keyword{regression}
\keyword{models}
\keyword{survival}
\keyword{nonparametric}
\keyword{bayesian}
