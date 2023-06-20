\name{bfi}
\alias{bfi}
\title{Bayesian Federated Inference}
\description{
\code{bfi} function is used (in central server) to estimate the parameters in the GLM and Survival models by BFI method using the aggregated results form the local centers.
}
\usage{
bfi(theta_hats = NULL, A_hats, Lambda, L = NULL,
    const_var = NULL, stratified = FALSE, nuisance = 1L)
}
\arguments{

\item{theta_hats}{a list of \eqn{L} vectors of the maximum a posteriori (MAP) estimates of \eqn{L} centers. These vectors must have equal dimensions. See \sQuote{Details}.}

\item{A_hats}{a list of \eqn{L} curvature matrices for \eqn{L} centers. These matrices must have equal dimensions. See \sQuote{Details}.}

\item{Lambda}{a list of \eqn{L+1} matrices. Each matrix used as the prior of the inverse variance-covariance matrix of Gaussian distribution. The first \eqn{L} matrices are for \eqn{L} local centers and the last one is associated with the combined data. These matrices must have equal dimensions. See \sQuote{Details}.}

\item{L}{number of locations/centers.}

\item{const_var}{a vector of \eqn{L} elements. Each element represents a specific feature of the corresponding center. There must be only one specific value or attribute for each center. This vector could be a numeric, characteristic or factor vector.}

\item{stratified}{logical flag for performing the stratified analysis. If \code{stratified=TRUE}, the parameters selected in the \code{nuisance} argument can be different across centers. Default is \code{stratified=FALSE}. See \sQuote{Details}.}

\item{nuisance}{a vector of one or two elements. It is useable if \code{stratified=TRUE}. For the  \code{binomial} family the length of the vector should be one which referes to 'intercept', and for \code{gaussian} this vector can have maximum two elements which refere to 'intercept' and/or 'sigma2'. See \sQuote{Details}.}
}

\value{

\code{bfi} returns a list containing the following components:

\item{theta_hat}{the \eqn{p}- or \eqn{(p + 1)}-dimensional vector of estimates obtained by the \code{'BFI'} method if \code{stratified=FALSE}. When \code{stratified=TRUE}, in this case dimension of this vector is \eqn{L+p-1} (for \code{binomial} family) or \eqn{2*L + p - 1} (for \code{gaussian} family). See \sQuote{Details};}

\item{A_hat}{the curvature matrix obtained by the \code{'BFI'} method if \code{stratified=FALSE}. In this case, it's a \eqn{p\times p}- or \eqn{(p + 1)\times (p+1)}-dimensional matrix depending on the family used. If \code{stratified=TRUE}, it's a list of \eqn{L} matrices corresponding to each center. This is not the same as \code{A_hats} in \sQuote{Arguments}. See \sQuote{Details};}

\item{sd}{the \eqn{p}- or \eqn{(p + 1)}-dimensional vector of standard deviation of estimates in \code{theta_hat} if \code{stratified=FALSE}, i.e., \code{sqrt(diag(solve(A_hat)))}.
If \code{stratified=TRUE}, it's a list of \eqn{L} vectors corresponding to each center. These vectors are standard deviation of parameter estimates obtained from the matrices in \code{A_hat}, i.e., \code{sqrt(diag(solve(A_hat[[j]])))} where \code{j} refers to a center.}
}

\details{
\code{bfi} function implements .....%the differential geometric generalization of the least angle regression method (Efron et al., 2004) proposed in Augugliaro et al. (2013) and Pazira et al. (2017).

%As in \dQuote{\code{\link{glm}}}, the user can specify family and link function using the argument \code{family}. When the \code{binomial} family is used, the responce can be a vector with entries 0/1 (failure/success) or, alternatively, a matrix where the first column is the number of \dQuote{successes} and the second column is the number of \dQuote{failures}. Starting with the version 2.0.0, the model can be specified combining family and link functions as describted in the following table:

%When gaussian, gamma or inverse.gaussian is used to model the error distribution, \code{dglars} returns the vector of the estimates of the dispersion parameter \eqn{\phi}{phi}; by default, the generalized Pearson statistic is used as estimator but the user can use the function \code{phihat} to specify other estimators (see \code{\link{phihat}} for more details).

Since the elements of the lists in the first three arguments of the \code{bfi} function must have equal dimensions, the function \code{MAP.estimation()} creates the outputs with the same dimensions across the centers even if in some of the centers \code{intercept=FALSE} (of course for the same \code{family}'s) by setting the corresponding elements to zero.
}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023). \emph{Bayesian Federated Inference for Statistical Models}, \emph{Statistics in Medicine}, Vol. 0(0), 0-0. <https://doi.org/10.48550/arXiv.2302.07677>
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
L <- 2    # L is the number of locations

##----------------
## Local center 1:
##----------------
set.seed(1123)
n1         <- 30
p          <- 4    # p (number of coefficients) is the same for all L locations.
X1         <- data.frame(matrix(rnorm(n1 * p), n1, p))
#true_beta <- c(1, 2, 0, 0, 0)
eta1       <- 1 + 2 * X1[,1]  ## with an intercept b0=1, b1=2, b2=b3=...=bp=0
mu1        <- binomial()$linkinv(eta1)
y1         <- rbinom(n1, 1, mu1)
lambda     <- 0.01
# we assume the same inverse covariance matrix for all locations.
Lambda      <- inv.prior.cov(X1, lambda, family=binomial)
fit1       <- MAP.estimation(y1, X1, family=binomial, Lambda)
theta_hat1 <- fit1$theta_hat # intercept and coefficient estimates
A_hat1     <- fit1$A_hat     # curvature matrix

##----------------
## Local center 2:
##----------------
n2         <- 50
p          <- 4    # p is the same for all L locations
X2         <- data.frame(matrix(rnorm(n2 * p), n2, p))
eta2       <- 1 + 2 * X2[,1]  ## with an intercept b0=1, b1=2,b2=b3=...=bp=0
mu2        <- binomial()$linkinv(eta2)
y2         <- rbinom(n2, 1, mu2)
fit2       <- MAP.estimation(y2, X2, family=binomial, Lambda)
theta_hat2 <- fit2$theta_hat
A_hat2     <- fit2$A_hat

##--------------
## Combined data
##--------------
y              <- c(y1, y2)
X              <- rbind(X1, X2)
fit_comb       <- MAP.estimation(y, X, family=binomial, Lambda)
theta_hat_comb <- fit_comb$theta_hat # beta estimates of combined data

###-----------------------------
### Bayesian Federated Inference
###-----------------------------
A_hats <- list(A_hat1, A_hat2)
theta_hats <- list(theta_hat1, theta_hat2)
bfi(theta_hats, A_hats, Lambda)

# theta (intercept and coefficient) estimates by BFI
(theta_hat_bfi <- bfi(theta_hats, A_hats, Lambda)$theta_hat)

# Curvature matrix estimate by BFI
bfi(theta_hats, A_hats=A_hats, Lambda=Lambda)$A_hat
# it is not needed to have the first argument ('theta_hats') for 'A_hats' and 'sd':
bfi(A_hats=A_hats, Lambda=Lambda)$A_hat

# SD of the BFI estimates
bfi(A_hats=A_hats, Lambda=Lambda)$sd

### Difference between BFI estimates and estimates with combined data
theta_hat_bfi - theta_hat_comb

}
\keyword{Bayesian}
\keyword{Federated}
