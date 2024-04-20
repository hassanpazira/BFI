\name{inv.prior.cov}
\alias{inv.prior.cov}
\title{Creates an inverse covariance matrix for a Gaussian prior}
\description{
\code{inv.prior.cov} constructs a diagonal inverse covariance matrix for the Gaussian prior distribution based on the design matrix of covariates. This construction accounts for the number of regression parameters, especially when dealing with categorical covariates. For a linear model, it also includes an additional row and column to represent the variance of the measurement error. In the case of a survival model, it considers the parameters of the baseline hazard function as well.
}
\usage{
inv.prior.cov(X,
              lambda = 1,
              L = 2,
              family = c("gaussian", "binomial", "survival"),
              treatment = NULL,
              treat_round = NULL,
              intercept = TRUE,
              stratified = FALSE,
              strat_par = NULL,
              center_spec = NULL,
              basehaz = c("weibul", "exp", "gomp", "poly", "pwexp", "unspecified"),
              max_order = 2,
              n_intervals = 4)
}

\arguments{
\item{X}{design matrix of dimension \eqn{n \times p}, where \eqn{n} is the number of samples observed, and \eqn{p} is the number of predictors/variables so excluding the intercept.}

\item{lambda}{the vector used as the diagonal of the (inverse covariance) matrix that will be created by \code{inv.prior.cov()}. The length of the vector depends on the number of columns of \code{X}, type of the covariates (continuous/dichotomous or categorical), \code{family}, whether an \code{intercept} is included in the model, and whether \code{stratified} analysis is desired. When \code{stratified = FALSE}, \code{lambda} could be a single positive number (if all values in the vector are equal), a vector of two elements (the first is used for regression parameters including ``intercept'' and the second for the ``sigma2'' in the gaussian family or for the baseline hazard parameters in the survival case), or a vector of length equal to the number of model parameters. However, the length of \code{lambda} is different when \code{stratified = TRUE}, see \sQuote{Details} for more information. Default is \code{lambda = 1}, which means all model parameters are set to one.}

\item{L}{the number of centers. This argument is used only when \code{stratified = TRUE}. Default is \code{L = 2}. See \sQuote{Details} and \sQuote{Examples}.}

\item{family}{a description of the error distribution. This is a character string naming a family of the model. In the current version of the package, the family of model can be \code{"gaussian"} (with \code{identity} link function), \code{"binomial"} (with \code{logit} link function), or \code{"survival"}. Can be abbreviated. By default the \code{gaussian} family is used. In case of a linear regression model, \code{family = "gaussian"}, there is an extra model parameter for the variance of measurement error. While in the case of survival model, \code{family = "survival"}, the number of the model parameters depend on the choice of baseline hazard functions, see \sQuote{Details} for more information.}

\item{treatment}{a character string representing the name or place of the binary covariate, respectively. This covariate indicates whether the patient got the new treatment (\eqn{z_{\ell i}=1}) or the placebo/standard treatment (\eqn{z_{\ell i}=0}). For both the first and second rounds, it should not be '\code{NULL}'. See \sQuote{Details}.}

\item{treat_round}{a character string representing the \code{'first'} or \code{'second'} round of estimating treatment effects. In the first round, \code{treat_round = 'first'}, the local estimates of the coefficients (\eqn{\gamma_{\ell}}) is estimated. In the second round, \code{treat_round = 'second'}, the \code{propensity} scores and the statistical summaries (\code{for_ATE}) are calculated.}

\item{intercept}{logical flag for having an intercept. It is not used when \code{family = "survival"}. By changing the \code{intercept} the dimension of the inverse covariance matrix changes. If \code{intercept = TRUE} (the default), the output matrix created by \code{inv.prior.cov()} has one row and one column related to \code{intercept}, while if \code{intercept = FALSE}, the resulting matrix does not have the row and column called \code{intercept}.}

\item{stratified}{logical flag for performing the stratified analysis. If \code{stratified = TRUE}, the parameter(s) selected in the \code{strat_par} argument are allowed to be different across centers to deal with heterogeneity across centers. This argument should only be used when designing the inverse covariance matrix for the (fictive) combined data, i.e., the last matrix for the \code{Lambda} argument in \code{bfi()}. If \code{inv.prior.cov()} is used for the analysis in the local centers (to build the \eqn{L} first matrices for the \code{Lambda} argument in \code{bfi()}), this argument should be \code{FALSE}, even if the BFI analysis is stratified. Default is \code{stratified = FALSE}. See \sQuote{Details} and \sQuote{Examples}.}

\item{strat_par}{a integer vector for indicating the stratification parameter(s). It can be used to deal with heterogeneity due to center-specific parameters. For the \code{"binomial"} and \code{"gaussian"} families it is a one- or two-element integer vector so that the values \eqn{1} and/or \eqn{2} are/is used to indicate that the ``intercept'' and/or ``sigma2'' are allowed to vary, respectively. For the \code{"binomial"} family the length of the vector should be one which refers to ``intercept'', and the value of this element should be \eqn{1} (to handel heterogeneity across outcome means). For \code{"gaussian"} this vector can be \eqn{1} for indicating the ``intercept'' only (handeling heterogeneity across outcome means), \eqn{2} for indicating the ``sigma2'' only (handeling heterogeneity due to nuisance parameter), and c(\eqn{1}, \eqn{2}) for both ``intercept'' and ``sigma2''. When \code{family = "survival"}, this vector can contain any combination from 1 to the maximum number of parameters of the baseline function, i.e., \eqn{1} for \code{"exp"}, \eqn{2} for \code{"weibul"} and \code{"gomp"}, \code{max_order + 1} for \code{"poly"}, and \code{n_intervals} for \code{"pwexp"}.
This argument is only used when \code{stratified = TRUE}. Default is \code{strat_par = NULL}. If \code{stratified = TRUE}, \code{strat_par} can not be \code{NULL} except when \code{center_spec} is not \code{NULL} for handeling heterogeneity due to clustering and missing covariates. See \sQuote{Details} and \sQuote{Examples}.}

\item{center_spec}{a vector of \eqn{L} elements for representing the center specific variable. This argument is used only when \code{stratified = TRUE} and \code{strat_par = NULL}. Each element represents a specific feature of the corresponding center. There must be only one specific value or attribute for each center. This vector could be a numeric, characteristic or factor vector. Note that, the order of the centers in the vector \code{center_spec} must be the same as in the list of the argument \code{theta_hats} in the function \code{bfi()}.
The used data type in the argument \code{center_spec} must be categorical. Default is \code{center_spec = NULL}. See also \sQuote{Details} and \sQuote{Examples}.}


\item{basehaz}{a character string representing one of the available baseline hazard functions; \code{exponential} (\code{"exp"}), \code{Weibull} (\code{"weibul"}, the default), \code{Gompertz} (\code{"gomp"}), \code{exponentiated polynomial} (\code{"poly"}), \code{piecewise constant exponential} (\code{"pwexp"}), and \code{unspecified baseline hazard} (\code{"unspecified"}). Can be abbreviated. It is only used when \code{family = "survival"}.}

\item{max_order}{an integer representing the maximum value of \code{q_l}, which is the order/degree minus 1 of the exponentiated polynomial baseline hazard function. This argument is only used when \code{family = "survival"} and \code{basehaz = "poly"}. Default is 2.}

\item{n_intervals}{an integer representing the number of intervals in the piecewise exponential baseline hazard function. This argument is only used when \code{family = "survival"} and \code{basehaz = "pwexp"}. Default is 4.}

}

\details{
\code{inv.prior.cov} creates a diagonal matrix with the vector \code{lambda} as its diagonal. The argument \code{stratified = TRUE} should only be used to construct a matrix for the prior density in case of stratification in the fictive combined data. Never be used for the construction of the matrix for analysis in the centers.

When \code{stratified = FALSE}, the length of the vector \code{lambda} depends on the covariate matrix \code{X}, \code{family}, \code{basehaz}, and whether an ``intercept'' is included in the model. For example, if the design matrix \code{X} has \code{p} columns with continuous or dichotomous covariates, \code{family = gaussian}, and \code{intercept = TRUE}, then \code{lambda} should have \eqn{p+2} elements. In this case, if in \code{X} there is a categorical covariate with \eqn{q>2} categories, then the length of \code{lambda} increases with \eqn{q-2}.

All values of \code{lambda} should be non-negative as they represent the inverse of the variance of the Gaussian prior. This argument is considered as the inverse of the variance of the prior distribution for: \eqn{(\beta_0, \boldsymbol{\beta})} if \code{family = "binomial"} and \code{intercept = TRUE}; \eqn{(\beta_0, \boldsymbol{\beta},\sigma^2)} if \code{family = "gaussian"} and \code{intercept = TRUE}; and \eqn{( \boldsymbol{\beta},\boldsymbol{\omega})} if \code{family = "survival"}.

If all values in the vector \code{lambda} equal, one value is enough to be given as entry.
If \code{lambda} is a scalar, the function \code{inv.prior.cov} sets each value at the diagonal equal to \code{lambda}.
%In the linear regression model the last parameter is assumed to be the inverse of the variance of the prior distribution for the measurement error.
When \code{lambda} is two dimensional: if \code{family = "binomial"}, the first and second values are used for the inverse of the variance of the prior distribution for the intercept (\eqn{\beta_0}) and regression parameters (\eqn{\boldsymbol{\beta}}), respectively;
If \code{family = "gaussian"}, the first and second values are used for the inverse of the variance of the prior distribution for the regression parameters including the intercept (\eqn{\beta_0, \boldsymbol{\beta}}) and variance of the measurement error (\eqn{ \sigma^2}), respectively;
If \code{family = "survival"}, the first and second values are used for the inverse of the variance of the prior distribution for the regression parameters (\eqn{\boldsymbol{\beta}}) and baseline hazard parameters (\eqn{ \omega}), respectively.
But if \code{stratified = TRUE} the length of the vector \code{lambda} must be equal to the number of parameters in the combined model.

If \code{intercept = FALSE}, for the \code{binomial} family the stratified analysis is not possible therefore \code{stratified} can not be \code{TRUE}.

If \code{stratified = FALSE}, both \code{strat_par} and \code{center_spec} must be \code{NULL} (the defaults), while if \code{stratified = TRUE} only one of the two must be \code{NULL}.

If \code{stratified = TRUE} and \code{family = "survival"}, \code{strat_par = 1} refers to \eqn{\omega_0} when \code{basehaz = "poly"}, and to \eqn{\omega_1} for other baseline hazards.

The output of \code{inv.prior.cov()} can be used in the main functions \code{MAP.estimation()} and \code{bfi()}.
}

\value{
\code{inv.prior.cov} returns a diagonal matrix. The dimension of the matrix depends on the number of columns of \code{X}, type of the covariates (continuous/dichotomous or categorical), \code{intercept}, \code{family}, and \code{basehaz}.
}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2024). \emph{Bayesian federated inference for estimating statistical models based on non-shared multicenter data sets}, \emph{Statistics in Medicine}, 43(12): 2421-2438. <https://doi.org/10.1002/sim.10072>

Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker M.A. (2025b). \emph{Bayesian Federated Inference for Survival Models}, \emph{Journal of Applied Statistics (Accepted)}. <https://arxiv.org/abs/2404.17464>

Jonker M.A., Pazira H. and Coolen A.C.C. (2025a). \emph{Bayesian Federated Inference for regression models based on non-shared medical center data}, \emph{Research Synthesis Methods}, 1-41. <https://doi.org/10.1017/rsm.2025.6>
}

\author{Hassan Pazira and Marianne Jonker \cr
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

# The (inverse) variance value (lambda=0.05) is assumed to be
# the same for Gaussian prior of all parameters (for non-stratified)

#-------------------------------------------------
# Inverse Covariance Matrix for the Gaussian prior
#-------------------------------------------------
# y ~ Binomial with 'intercept'
inv.prior.cov(X, lambda = 0.05, family = 'binomial')
# returns a 5-by-5 matrix

# y ~ Binomial without 'intercept'
inv.prior.cov(X, lambda = 0.05, family = "binomial", intercept = FALSE)
# a 4-by-4 matrix

# y ~ Gaussian with 'intercept'
inv.prior.cov(X, lambda = 0.05, family = 'gaussian')
# returns a 6-by-6 matrix

# Survival family with 'weibul' baseline hazard
inv.prior.cov(X, lambda = c(0.05, 0.1), family = 'survival')
# returns a 6-by-6 matrix

# Survival family with 'pwexp' baseline hazard (4 intervals)
inv.prior.cov(X, lambda = 0.05, family = 'survival', basehaz = "pwexp")
# returns a 8-by-8 matrix

# Survival family with 'poly' baseline hazard
inv.prior.cov(X, lambda = c(0.05, 0.1), family = 'survival', basehaz = "poly")
# returns a 7-by-7 matrix

#--------------------
# Stratified analysis
#--------------------
# y ~ Binomial when 'intercept' varies across 3 centers:
inv.prior.cov(X, lambda = c(.2, 1), family = 'binomial', stratified = TRUE,
              strat_par = 1, L = 3)

# y ~ Gaussian when 'intercept' and 'sigma2' vary across 2 centers; y ~ Gaussian
inv.prior.cov(X, lambda = c(1, 2, 3), family = "gaussian", stratified = TRUE,
              strat_par = c(1, 2))

# y ~ Gaussian when 'sigma2' varies across 2 centers (with 'intercept')
inv.prior.cov(X, lambda = c(1, 2, 3), family='gaussian', stratified = TRUE,
              strat_par = 2)

# y ~ Gaussian when 'sigma2' varies across 2 centers (without 'intercept')
inv.prior.cov(X, lambda = c(2, 3), family = "gaussian", intercept = FALSE,
              stratified=TRUE, strat_par = 2)

#--------------------------
# Center specific covariate
#--------------------------
# center specific covariate has K = 2 categories across 4 centers; y ~ Binomial
inv.prior.cov(X, lambda = c(0.1:2), family = 'binomial', stratified = TRUE,
              center_spec = c("Iran","Netherlands","Netherlands","Iran"), L=4)

# center specific covariate has K = 3 categories across 5 centers; y ~ Gaussian
inv.prior.cov(X, lambda = c(0.5:3), family = 'gaussian', stratified = TRUE,
              center_spec = c("Medium","Big","Small","Big","Small"), L = 5)

# center specific covariate has K = 4 categories across 5 centers; y ~ Gaussian
inv.prior.cov(X, lambda = 1, family = 'gaussian', stratified = TRUE,
              center_spec = c(3,1:4), L=5)

}
\keyword{array}
