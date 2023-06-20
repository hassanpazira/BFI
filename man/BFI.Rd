\name{bfi}
\alias{bfi}
\title{Bayesian Federated Inference}
\description{
\code{bfi} function can be used (in the central server) to combine inference results from separate data sets (without combining the data) to approximate what would have been inferred had the data sets been merged. For now the function can handle linear and logistic regression models, but code for more models will be available in the near future.
bfi command
}
\usage{
bfi(theta_hats = NULL, A_hats, Lambda, stratified = FALSE,
    strat_par = NULL, center_spec = NULL, type=c("cat", "con"))
}
\arguments{

\item{theta_hats}{a list of \eqn{L} vectors of the maximum a posteriori (MAP) estimates of the model parameters in the \eqn{L} centers. These vectors must have equal dimensions. See \sQuote{Details}.}

\item{A_hats}{a list of \eqn{L} (minus) curvature matrices for \eqn{L} centers. These matrices must have equal dimensions. See \sQuote{Details}.}

\item{Lambda}{a list of \eqn{L+1} matrices. The \eqn{k^{\th}} matrix is the chosen inverse variance-covariance matrix of the Gaussian distribution that is used as prior distribution in center \eqn{k}, where \eqn{k=1,2,\ldots,L}. The last matrix is the chosen variance-covariance matrix for the Gaussian prior of the (fictive) combined data set.
If \code{stratified = FALSE}, all \eqn{L+1} matrices must have equal dimensions. While, if \code{stratified = TRUE}, the first \eqn{L} matrices must have equal dimensions and the last matrix should have a different (greater) dimention than the others. See \sQuote{Details}.}

\item{stratified}{logical flag for performing the stratified analysis. If \code{stratified = TRUE}, the parameter(s) selected in the \code{strat_par} argument are allowed to be different across centers, except when the argument \code{center_spec} is not \code{NULL}. Default is \code{stratified = FALSE}. See \sQuote{Details} and \sQuote{Examples}.}

\item{strat_par}{a one- or two-element integer vector for indicating the stratification parameter(s). The values \eqn{1} and/or \eqn{2} are/is used to indicate that the ``intercept'' and/or ``sigma2'' are allowed to vary, respectively. This argument is used only when \code{stratified = TRUE} and \code{center_spec = NULL}. Default is \code{strat_par = NULL}, but if \code{stratified = TRUE}, \code{strat_par} can not be \code{NULL} unless there is a center specific variable. For the \code{binomial} family the length of the vector should be at most one which refers to ``intercept'', and the value of this element should be \eqn{1}. For \code{gaussian} family this vector can be \eqn{1} for indicating the ``intercept'' only, \eqn{2} for indicating the ``sigma2'' only, and c(\eqn{1}, \eqn{2}) for both ``intercept'' and ``sigma2''. See \sQuote{Details} and \sQuote{Examples}.}

\item{center_spec}{a vector of \eqn{L} elements for representing the center specific variable. This argument is used only when \code{stratified = TRUE} and \code{strat_par = NULL}. Each element represents a specific feature of the corresponding center. There must be only one specific value or attribute for each center. This vector could be a numeric, characteristic or factor vector. Note that, the order of the centers in the vector \code{center_spec} must be the same as in the list of the argument \code{theta_hats}. If \code{center_spec = NULL} (the default) and \code{stratified = TRUE}, the stratification parameter(s) given in the argument \code{strat_par} are assumed to be different across the centers. See \sQuote{Details} and \sQuote{Examples}.}

\item{type}{a description of the used data type in the \code{center_spec} argument. Type for 'categorical' (\code{type = "cat"}) or 'continuous' (\code{type = "con"}) data. Can be also abbreviated. This argument is used only when \code{center_spec} is not \code{NULL}.}

}

\value{

\code{bfi} returns a list containing the following components:

\item{theta_hat}{the vector of estimates obtained by combining the inference results from the \eqn{L} centers with the \code{'BFI'} methodology. If an intercept was fitted in every center and \code{stratified = FALSE}, there is only one general ``intercept'' in this vector, while if \code{stratified = TRUE} and \code{strat_par = 1}, there are \eqn{L} different intercepts in the model, for each center one;}

\item{A_hat}{minus the curvature (or Hessian) matrix obtained by the \code{'BFI'} method for the combined model. If \code{stratified = TRUE}, the dimension of the matrix is always greater than when \code{stratified = FALSE}. When \code{stratified = TRUE}, some elements at the off-diagonal matrix can not be determined and are therefore \code{NA};}

\item{sd}{the vector of standard deviation of the estimates in \code{theta_hat} obtained from the matrix in \code{A_hat}, i.e., the vector equals \code{sqrt(diag(solve(A_hat)))} which equals the square root of the elements at the diagonal of the inverse of the \code{A_hat} matrix.}

}

\details{
\code{bfi} function implements the BFI approach described in the paper Jonker et. al. (2023) given in the references.
The inference results gathered from different (\eqn{L}) centers are combined, and the BFI estimates of the model parameters and curvature matrix evaluated at that point are returned.

The inference result from each center must be obtained using the \code{MAP.estimation} function separately, and then all of these results (coming from different centers) should be compiled into a list to be used as an input of \code{bfi()}.
The models in the different centers should be defined in exactly the same way; among others, exactly the same covariates should be included in the models. The parameter vectors should be defined exactly the same, so that the L vectors and matrices in the input lists \code{theta_hat}'s and \code{A_hat}'s are defined in the same way (e.g., the covariates need to be included in the models in the same order).

Note that the order of the elements in the lists \code{theta_hats}, \code{A_hats} and \code{Lambda}, must be the same with respect to the centers, so that in every list the element at the \eqn{\ell^{\th}} position is from the center \eqn{\ell}. This should also be the case for the vector \code{center_spec}.

If \code{intercept = FALSE} for the locations, the stratified analysis is not possible anymore for the \code{binomial} family.

If \code{stratified = FALSE}, both \code{strat_par} and \code{center_spec} must be \code{NULL} (the defaults), while if \code{stratified = TRUE} only one of the two must be \code{NULL}.

If \code{stratified = FALSE} and all the \eqn{L+1} matrices are equal, it is sufficient to give a (list of) one matrix only.
In both cases of the \code{stratified} argument (\code{TRUE} or \code{FALSE}), if only the first \eqn{L} matrices are equal, the argument \code{Lambda} can be a list of two matrices, so that the fist matrix represents the chosen variance-covariance matrix for local centers and the second one is the chosen matrix for the combined data set.
The last matrix of the list in the argument \code{Lambda} can be built by the function \code{inv.prior.cov()}.

}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023). \emph{Bayesian Federated Inference for Statistical Models}. {Statistics in Medicine}, Vol. 0(0), 0-0. <https://doi.org/10.48550/arXiv.2302.07677>
}

\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{MAP.estimation}} and \code{\link{inv.prior.cov}}
}

\examples{
#-------------
# y ~ Binomial
#-------------
# model assumption:
beta <- 1:4  # regression coefficients (beta[1] is the intercept)

#-----------------------------------
# Data Simulation for local center 1
#-----------------------------------
n1 <- 30                                           # sample size of center 1
X1 <- data.frame(x1=rnorm(n1),                     # continuous variable
                 x2=sample(0:2, n1, replace=TRUE)) # categorical variable
# make dummy variables
X1x2_1 <- ifelse(X1$x2 == '1', 1, 0)
X1x2_2 <- ifelse(X1$x2 == '2', 1, 0)
X1$x2  <- as.factor(X1$x2)
# linear predictor:
eta1   <- beta[1] + X1$x1 * beta[2] + X1x2_1 * beta[3] + X1x2_2 * beta[4]
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu1    <- binomial()$linkinv(eta1)
y1     <- rbinom(n1, 1, mu1)

#-----------------------------------
# Data Simulation for local center 2
#-----------------------------------
n2 <- 50                                           # sample size of center 2
X2 <- data.frame(x1=rnorm(n2),                     # continuous variable
                 x2=sample(0:2, n2, replace=TRUE)) # categorical variable
# make dummy variables:
X2x2_1 <- ifelse(X2$x2 == '1', 1, 0)
X2x2_2 <- ifelse(X2$x2 == '2', 1, 0)
X2$x2  <- as.factor(X2$x2)
# linear predictor:
eta2   <- beta[1] + X2$x1 * beta[2] + X2x2_1 * beta[3] + X2x2_2 * beta[4]
# inverse of the link function:
mu2    <- binomial()$linkinv(eta2)
y2     <- rbinom(n2, 1, mu2)

#--------------------------
# MAP estimates at center 1
#--------------------------
# assume the same inverse covariance matrix (Lambda) for both centers:
Lambda     <- inv.prior.cov(X1, lambda=0.01, family=binomial)
fit1       <- MAP.estimation(y1, X1, family=binomial, Lambda)
theta_hat1 <- fit1$theta_hat # intercept and coefficient estimates
A_hat1     <- fit1$A_hat     # minus curvature matrix

#--------------------------
# MAP estimates at center 2
#--------------------------
fit2       <- MAP.estimation(y2, X2, family=binomial, Lambda)
theta_hat2 <- fit2$theta_hat
A_hat2     <- fit2$A_hat

#----------------------
# BFI at central center
#----------------------
A_hats     <- list(A_hat1, A_hat2)
theta_hats <- list(theta_hat1, theta_hat2)
bfi        <- bfi(theta_hats, A_hats, Lambda)
class(bfi)
summary(bfi, cur_mat=TRUE)

#--------------------
# stratified analysis
#--------------------
# By running the following line an error appears because when
# stratified = TRUE, both 'strat_par' and 'center_spec' can not be NULL
# bfi(theta_hats, A_hats, Lambda, stratified=TRUE)
#bfi(theta_hats, A_hats, Lambda, stratified=TRUE, strat_par=1)

}
\keyword{Bayesian}
\keyword{Federated}
