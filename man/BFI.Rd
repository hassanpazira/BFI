\name{bfi}
\alias{bfi}
\title{Bayesian Federated Inference}
\description{
\code{bfi} function can be used (in central server) to combine inference results from separate data sets (without combining the data) to approximate what would have been inferred had the data sets been merged. For now the function can handle linear and logistic regression models, but more models will be available in the near future.
bfi command
}
\usage{
bfi(theta_hats = NULL, A_hats, Lambda, const_var = NULL,
    stratified = FALSE, nuisance = 1L)
}
\arguments{

\item{theta_hats}{a list of \eqn{L} vectors of the maximum a posteriori (MAP) estimates of the model parameters in the \eqn{L} centers. These vectors must have equal dimensions. See \sQuote{Details}.}

\item{A_hats}{a list of \eqn{L} curvature matrices for \eqn{L} centers. These matrices must have equal dimensions. See \sQuote{Details}.}

\item{Lambda}{a list of \eqn{L+1} matrices. Each matrix used as the prior of the inverse variance-covariance matrix of Gaussian distribution. The first \eqn{L} matrices are for \eqn{L} local centers and the last one is the chosen variance-covariance matrix for the Gaussian prior of the (fictive) combined data set. These matrices must have equal dimensions. If all the \eqn{L+1} matrices are the same, this argument can be one matrix or a list of one matrix. See \sQuote{Details}.}

\item{const_var}{a vector of \eqn{L} elements. Each element represents a specific feature of the corresponding center. There must be only one specific value or attribute for each center. This vector could be a numeric, characteristic or factor vector.}

\item{stratified}{logical flag for performing the stratified analysis. If \code{stratified=TRUE}, the parameters selected in the \code{nuisance} argument can be different across centers. Default is \code{stratified=FALSE}. See \sQuote{Details}.}

\item{nuisance}{a vector of one or two elements. It is usable if \code{stratified=TRUE}. For the  \code{binomial} family the length of the vector should be one which refers to 'intercept', and for \code{gaussian} this vector can have maximum two elements which refer to 'intercept' and/or 'sigma2'. See \sQuote{Details}.}
}

\value{

\code{bfi} returns a list containing the following components:

\item{theta_hat}{the \eqn{p}- or \eqn{(p + 1)}-dimensional vector of estimates obtained by combining the inference results from the \eqn{L} centers with the \code{'BFI'} methodology when \code{stratified=FALSE}, where \eqn{p} is the number of regression parameters including intercept. If \code{stratified=TRUE}, dimension of this vector is \eqn{L+p-1} (for \code{binomial} family, or for \code{gaussian} family when \code{nuisance} has one element) or \eqn{2*L + p - 1} (for \code{gaussian} family when \code{nuisance} has two elements). See \sQuote{Details};}

\item{A_hat}{the curvature matrix obtained by the \code{'BFI'} method when \code{stratified=FALSE}, such that it's a \eqn{(p\times p)}- or \eqn{(p + 1)\times (p+1)}-dimensional matrix depending on the family used, where \eqn{p} is the number of regression parameters including intercept. If \code{stratified=TRUE}, it's a list of \eqn{L} matrices corresponding to each center. This is not the same as \code{A_hats} in \sQuote{Arguments}. See \sQuote{Details};}

\item{sd}{the \eqn{p}- or \eqn{(p + 1)}-dimensional vector of standard deviation of the estimates in \code{theta_hat} if \code{stratified=FALSE}, i.e., the vector equals \code{sqrt(diag(solve(A_hat)))}.
If \code{stratified=TRUE}, it's a list of \eqn{L} vectors; one vector for every center. These vectors are standard deviation of parameter estimates obtained from the matrices in \code{A_hat}, i.e., the vector equals \code{sqrt(diag(solve(A_hat[[j]])))} where \code{j} refers to a center.}
}

\details{
\code{bfi} function implements the BFI approach described in the paper Jonker et. al. (2023) given in the references.
The results gathered from different (\eqn{L}) centers are combined, and the BFI estimates of the model parameters and curvature matrix evaluated at that point are returned.

The result from each center must be obtained using the \code{MAP.estimation} function separately, and then all of these results (coming from different centers) should be compiled into a list to be used as an input of \code{bfi()}.
The models in the different centers should be defined in exactly the same way; among others, exactly the same covariates should be included in the models. The parameter vectors should be defined exactly the same, so that the estimates \code{theta_hat}'s and \code{A_hat}'s (the output of the \code{MAP.estimation} function) are defined in the same way.

Note that, the order of the estimates in the lists \code{theta_hats}, \code{A_hats} and \code{Lambda} must be the same, so that in every list the elements at the \eqn{\ell^{\th}}  position are all from the same center. For example, the first three entries of \code{bfi()} could be as follows:

\code{bfi(theta_hats=theta_list, A_hats=A_list, Lambda= Lambda_list)}

where

\code{theta_list=list(obj}\eqn{_1}\code{$theta_hat, ..., obj}\eqn{_\ell}\code{$theta_hat, ..., obj}\eqn{_L}\code{$theta_hat)},

\code{A_list=list(obj}\eqn{_1}\code{$A_hat,..., obj}\eqn{_\ell}\code{$A_hat,..., obj}\eqn{_L}\code{$A_hat)},

\code{Lambda_list=list(obj}\eqn{_1}\code{$Lambda,..., obj}\eqn{_\ell}\code{$Lambda,..., obj}\eqn{_L}\code{$Lambda)},

and

\code{obj}\eqn{_\ell}\code{=MAP.estimation(y=y}\eqn{_\ell}\code{, X=X}\eqn{_\ell}\code{, Lambda=Lambda}\eqn{_\ell}\code{)} # must be done at the center \eqn{\ell}

where \code{y}\eqn{_\ell}, \code{X}\eqn{_\ell} and \code{Lambda}\eqn{_\ell} are the data and the prior of the inverse variance-covariance matrix for the \eqn{\ell^{\th}} center.


Details related to \code{stratified=TRUE} should be added!!!....

Details related to \code{const_var} should be added!!!....
}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023). \emph{Bayesian Federated Inference for Statistical Models}. {Statistics in Medicine}, Vol. 0(0), 0-0. <https://doi.org/10.48550/arXiv.2302.07677>
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
set.seed(112358)
p <- 4 # number of regression coefficients (including intercept), is the same for all L=2 locations

##---------------
## Local center 1
##---------------
n1   <- 30
X1   <- data.frame(x1=rnorm(n1),                     # standard normal variable
                   x2=sample(0:1, n1, replace=TRUE), # dichotomous variable
                   x3=sample(1:3, n1, replace=TRUE)) # categorical variable
#true_beta <- rep(2, p)           # with an intercept b0=b1=b2=b3=2
eta1 <- 2 + apply(2 * X1, 1, sum) # linear predictor: X%*%\beta=\eta
mu1  <- binomial()$linkinv(eta1)  # inverse of the link function: g^{-1}(\eta)=\mu
y1   <- rbinom(n1, 1, mu1)
# we assume the same inverse covariance matrix for all locations
Lambda     <- inv.prior.cov(X1, 0.01, family=binomial)
fit1       <- MAP.estimation(y1, X1, family=binomial, Lambda)
theta_hat1 <- fit1$theta_hat # intercept and coefficient estimates
A_hat1     <- fit1$A_hat     # curvature matrix

##---------------
## Local center 2
##---------------
n2         <- 50
X2         <- data.frame(x1=rnorm(n2),                     # standard normal variable
                         x2=sample(0:1, n2, replace=TRUE), # dichotomous variable
                         x3=sample(1:3, n2, replace=TRUE)) # categorical variable
eta2       <- 2 + apply(2 * X2, 1, sum) # linear predictor: X%*%\beta=\eta
mu2        <- binomial()$linkinv(eta2)  # inverse of the link function: g^{-1}(\eta)=\mu
y2         <- rbinom(n2, 1, mu2)
fit2       <- MAP.estimation(y2, X2, family=binomial, Lambda)
theta_hat2 <- fit2$theta_hat
A_hat2     <- fit2$A_hat

###-----------------------------
### Bayesian Federated Inference
###-----------------------------
A_hats     <- list(A_hat1, A_hat2)
theta_hats <- list(theta_hat1, theta_hat2)
bfi(theta_hats, A_hats, Lambda)

}
\keyword{Bayesian}
\keyword{Federated}
