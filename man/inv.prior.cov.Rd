\name{inv.prior.cov}
\alias{inv.prior.cov}
\title{Creates an inverse covariance matrix for a Gaussian prior}
\description{
\code{inv.prior.cov} builds a diagonal inverse covariance matrix for the Gaussian prior distribution based on the design matrix of covariates, that takes into account the number of regression parameters in case of categorical covariates.  In case of a linear model, it also includes a row and column for the variance of the measurement errors.
}
\usage{
inv.prior.cov(X, lambda = 1, L = 2, family = gaussian, intercept = TRUE,
              stratified = FALSE, strat_par = NULL, center_spec = NULL)
}

\arguments{
\item{X}{design matrix of dimension \eqn{n \times p}, where \eqn{n} is the number of samples observed, and \eqn{p} is the number of predictors/variables so excluding the intercept.}

\item{lambda}{the vector used as the diagonal of the (inverse covariance) matrix that will be created by \code{inv.prior.cov()}. The length of the vector depends on the number of columns of \code{X}, type of the covariates (continuous/dichotomous or categorical), \code{family}, whether an \code{intercept} is included in the model, and whether \code{stratified} analysis is desired. When \code{stratified = FALSE}, \code{lambda} could be a single positive number (if all values in the vector are equal), a vector of two elements (the first is used for regression parameters including ``intercept'' and the second for the ``sigma2''), or a vector of length equal to the number of model parameters. However, the length of \code{lambda} is different when \code{stratified = TRUE}, see \sQuote{Details} for more information. Default is \code{lambda = 1}.}

\item{L}{the number of centers. This argument is used only when \code{stratified = TRUE}. Default is \code{L = 2}. See \sQuote{Details} and \sQuote{Examples}.}

\item{family}{a description of the error distribution and link function used to specify the model. This can be a character string naming a family function or the result of a call to a family function (see \code{\link{family}} for details). In the current version, the family of model can be \code{gaussian} (with \code{identity} link function) and \code{binomial} (with \code{logit} link function). By default the \code{gaussian} family is used. In case of a linear regression model, \code{family = gaussian}, there is an extra model parameter for the variance of measurement error.}

\item{intercept}{logical flag for having an intercept. By changing the \code{intercept} the dimension of the inverse covariance matrix changes. If \code{intercept = TRUE} (the default), the output matrix created by \code{inv.prior.cov()} has one row and one column related to \code{intercept}, while if \code{intercept = FALSE}, the resulting matrix does not have the row and column called \code{intercept}.}

\item{stratified}{logical flag for performing the stratified analysis. If \code{stratified = TRUE}, the parameter(s) selected in the \code{strat_par} argument are allowed to be different across centers. This argument should only be used when designing the inverse covariance matrix for the (fictive) combined data, i.e., the last matrix for the \code{Lambda} argument in \code{bfi()}. If \code{inv.prior.cov()} is used for the analysis in the local centers (to built the \eqn{L} first matrices for the \code{Lambda} argument in \code{bfi()}), this argument should be \code{FALSE}, even if the BFI analysis is stratified. Default is \code{stratified = FALSE}. See \sQuote{Details} and \sQuote{Examples}.}

\item{strat_par}{a one- or two-element integer vector for indicating the stratification parameter(s). The values \eqn{1} and/or \eqn{2} are/is used to indicate that the ``intercept'' and/or ``sigma2'' are allowed to vary, respectively. This argument is used only when \code{stratified = TRUE}. Default is \code{strat_par = NULL}, but if \code{stratified = TRUE}, \code{strat_par} can not be \code{NULL}. For the \code{binomial} family the length of the vector should be one which refers to ``intercept'', and the value of this element should be \eqn{1}. For \code{gaussian} this vector can be \eqn{1} for indicating the ``intercept'' only, \eqn{2} for indicating the ``sigma2'' only, and c(\eqn{1}, \eqn{2}) for both ``intercept'' and ``sigma2''. See \sQuote{Examples}.}

\item{center_spec}{a vector of \eqn{L} elements for representing the center specific variable. This argument is used only when \code{stratified = TRUE} and \code{strat_par = NULL}. Each element represents a specific feature of the corresponding center. There must be only one specific value or attribute for each center. This vector could be a numeric, characteristic or factor vector. Note that, the order of the centers in the vector \code{center_spec} must be the same as in the list of the argument \code{theta_hats} in the function \code{bfi()}.
The used data type in the argument \code{center_spec} must be categorical. Default is \code{center_spec = NULL}. See also \sQuote{Details} and \sQuote{Examples}.}

}

\details{
\code{inv.prior.cov} creates a diagonal matrix with the vector \code{lambda} as its diagonal. The argument \code{stratified = TRUE} should only be used to construct a matrix for the prior density in case of stratification in the fictive combined data. Never be used for the construction of the matrix for analysis in the centers.

When \code{stratified = FALSE}, the length of the vector \code{lambda} depends on the covariate matrix \code{X}, \code{family}, and whether an ``intercept'' is included in the model. For example, if the design matrix \code{X} has \code{p} columns with continuous or dichotomous covariates, \code{family = gaussian}, and \code{intercept = TRUE}, then \code{lambda} should have \eqn{p+2} elements. In this case, if in \code{X} there is a categorical covariate with \eqn{q>2} categories, then the length of \code{lambda} increases with \eqn{q-2}. All values of lambda should be non-negative as they represent the inverse of the variance of the Gaussian prior.
Note that, if all values in the vector \code{lambda} equal, one value is enough to be given as entry.
If \code{lambda} is a scalar, the function \code{inv.prior.cov} sets each value at the diagonal equal to \code{lambda}. In the linear regression model the last parameter is assumed to be the inverse of the variance of the prior distribution for the measurement error. If \code{lambda} is two dimensional, the first value is used for the prior of the regression parameters and the second for the inverse of the variance of the prior distribution for the measurement error.

If \code{stratified = TRUE} the length of the vector \code{lambda} should be equal to the number of parameters in the combined model.

If \code{intercept = FALSE}, for the \code{binomial} family the stratified analysis is not possible therefore \code{stratified} can not be \code{TRUE}.

If \code{stratified = FALSE}, both \code{strat_par} and \code{center_spec} must be \code{NULL} (the defaults), while if \code{stratified = TRUE} only one of the two must be \code{NULL}.

The output of \code{inv.prior.cov()} can be used in the main functions \code{MAP.estimation()} and \code{bfi()}.
}

\value{
\code{inv.prior.cov} returns a diagonal matrix. The dimension of the matrix depends on the number of columns of \code{X}, type of the covariates (continuous/dichotomous or categorical), \code{family}, and \code{intercept}.
}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023). \emph{Bayesian Federated Inference for Statistical Models},
{Statistics in Medicine}, Vol. 0(0), 0-0. <https://doi.org/10.48550/arXiv.2302.07677>
}

\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{MAP.estimation}}
}

\examples{
#----------------
# Data Simulation
#----------------
X <- data.frame(x1=rnorm(50),                     # standard normal variable
                x2=sample(0:2, 50, replace=TRUE), # categorical variable
                x3=sample(0:1, 50, replace=TRUE)) # dichotomous variable
X$x2 <- as.factor(X$x2)
X$x3 <- as.factor(X$x3)

#---------------------
# Load the BFI package
#---------------------
library(BFI)

# The (inverse) variance value (lambda=0.05) is assumed to be
# the same for Gaussian prior of all parameters (for non-stratified)

#-------------------------------------------------
# Inverse Covariance Matrix for the Gaussian prior
#-------------------------------------------------
# y ~ Binomial with 'intercept'
inv.prior.cov(X, lambda=0.05, family=binomial) # returns a 5-by-5 matrix

# y ~ Binomial without 'intercept'
inv.prior.cov(X, lambda=0.05, family="binomial", intercept = FALSE) # a 4-by-4 matrix

# y ~ Gaussian with 'intercept'
inv.prior.cov(X, lambda=0.05, family=gaussian) # returns a 6-by-6 matrix

#--------------------
# Stratified analysis
#--------------------
# y ~ Binomial when 'intercept' varies across 3 centers:
inv.prior.cov(X, lambda=c(.2, 1), family=binomial, stratified=TRUE, strat_par = 1, L = 3)

# y ~ Gaussian when 'intercept' and 'sigma2' vary across 2 centers; y ~ Gaussian
inv.prior.cov(X, lambda=c(1, 2, 3), family=gaussian, stratified=TRUE, strat_par = c(1, 2))

# y ~ Gaussian when 'sigma2' varies across 2 centers (with 'intercept')
inv.prior.cov(X, lambda=c(1, 2, 3), family=gaussian, stratified=TRUE, strat_par = 2)

# y ~ Gaussian when 'sigma2' varies across 2 centers (without 'intercept')
inv.prior.cov(X, lambda=c(2, 3), family=gaussian, intercept = FALSE, stratified=TRUE,
              strat_par = 2)

#--------------------------
# Center specific covariate
#--------------------------
# center specific covariate has K=2 categories across 4 centers; y ~ Binomial
inv.prior.cov(X, lambda=c(0.1:2), family=binomial, stratified=TRUE,
              center_spec = c("Iran","Netherlands","Netherlands","Iran"), L=4)

# center specific covariate has K=3 categories across 5 centers; y ~ Gaussian
inv.prior.cov(X, lambda=c(0.5:3), family=gaussian, stratified=TRUE,
              center_spec = c("Medium","Big","Small","Big","Small"), L=5)

# center specific covariate has K=4 categories across 5 centers; y ~ Gaussian
inv.prior.cov(X, lambda=1, family=gaussian, stratified=TRUE, center_spec = c(3,1:4), L=5)

}
\keyword{Bayesian}
