\name{bfi}
\alias{bfi}
\title{Bayesian Federated Inference}
\description{
\code{bfi} function can be used (on the central server) to combine inference results from separate datasets (without combining the data) to approximate what would have been inferred had the datasets been merged. This function can handle linear, logistic and survival regression models.
}
\usage{
bfi(theta_hats = NULL,
    A_hats,
    Lambda,
    family = c("gaussian", "binomial", "survival"),
    basehaz = c("weibul", "exp", "gomp", "poly", "pwexp", "unspecified"),
    stratified = FALSE,
    strat_par = NULL,
    center_spec = NULL,
    theta_A_polys = NULL,
    treat_round = NULL,
    for_ATE = NULL,
    p,
    q_ls,
    center_zero_sample = FALSE,
    which_cent_zeros,
    zero_sample_covs,
    refer_cats,
    zero_cats,
    lev_no_ref_zeros)
}
\arguments{

\item{theta_hats}{a list of \eqn{L} vectors of the maximum a posteriori (MAP) estimates of the model parameters in the \eqn{L} centers. These vectors must have equal dimensions. See \sQuote{Details}.}

\item{A_hats}{a list of \eqn{L} minus curvature matrices for \eqn{L} centers. These matrices must have equal dimensions. See \sQuote{Details}.}

\item{Lambda}{a list of \eqn{L+1} matrices. The \eqn{k^{\th}} matrix is the chosen inverse variance-covariance matrix of the Gaussian distribution that is used as prior distribution in center \eqn{k}, where \eqn{k=1,2,\ldots,L}. The last matrix is the chosen variance-covariance matrix for the Gaussian prior of the (fictive) combined data set.
If \code{stratified = FALSE}, all \eqn{L+1} matrices must have equal dimensions. While, if \code{stratified = TRUE}, the first \eqn{L} matrices must have equal dimensions and the last matrix should have a different (greater) dimention than the others. See \sQuote{Details}.}

\item{family}{a character string representing the family name used for the local centers. Can be abbreviated.}

\item{basehaz}{a character string representing one of the available baseline hazard functions; \code{exponential} (\code{"exp"}), \code{Weibull} (\code{"weibul"}, the default), \code{Gompertz} (\code{"gomp"}), \code{exponentiated polynomial} (\code{"poly"}), \code{piecewise exponential} (\code{"pwexp"}), and \code{unspecified baseline hazard} (\code{"unspecified"}). It is only used when \code{family = "survival"}. Can be abbreviated. If \code{basehaz = "unspecified"}, it means that a (semi-parametric) Cox model is considered, and the parameters (regression coefficients) are estimated using the partial log-likelihood.}

\item{stratified}{logical flag for performing the stratified analysis. If \code{stratified = TRUE}, the parameter(s) selected in the \code{strat_par} argument are allowed to be different across centers (to deal with heterogeneity across centers), except when the argument \code{center_spec} is not \code{NULL}. Default is \code{stratified = FALSE}. See \sQuote{Details} and \sQuote{Examples}.}

\item{strat_par}{an integer vector for indicating the stratification parameter(s). It can be used to deal with heterogeneity due to center-specific parameters.
For the \code{"binomial"} and \code{"gaussian"} families it is a one- or two-element integer vector so that the values \eqn{1} and/or \eqn{2} are/is used to indicate that the ``intercept'' and/or ``sigma2'' are allowed to vary, respectively. For the \code{"binomial"} family the length of the vector should be at most one which refers to ``intercept'', and the value of this element should be \eqn{1} (to handel heterogeneity across outcome means). For \code{"gaussian"} this vector can be \eqn{1} for indicating the ``intercept'' only (handeling heterogeneity across outcome means), \eqn{2} for indicating the ``sigma2'' only (handeling heterogeneity due to nuisance parameter), and \eqn{c(1, 2)} for both ``intercept'' and ``sigma2''. When \code{family = "survival"}, this vector can contain any combination of values ranging from 1 to the maximum number of parameters of the baseline hazard function, i.e., \eqn{1} for \code{"exp"}, \eqn{2} for \code{"weibul"} and \code{"gomp"}, \code{max_order + 1} for \code{"poly"}, and \code{n_intervals} for \code{"pwexp"}.
For example, for \code{"weibul"}, \code{strat_par} could be \eqn{1}, \eqn{2} or \eqn{c(1, 2)}, where \eqn{1} represents \eqn{\omega_1} and \eqn{2} represents \eqn{\omega_2}.
This argument is used only when \code{stratified = TRUE} and \code{center_spec = NULL}. Default is \code{strat_par = NULL}. See \sQuote{Details} and \sQuote{Examples}.}

\item{center_spec}{a vector of \eqn{L} elements to account for the heterogeneity across centers due to clustering. This argument is used only when \code{stratified = TRUE} and \code{strat_par = NULL}. Each element represents a specific feature of the corresponding center. There must be only one specific value or attribute for each center. This vector could be a numeric, characteristic or factor vector. Note that, the order of the centers in the vector \code{center_spec} must be the same as in the list of the argument \code{theta_hats}.
The used data type in the argument \code{center_spec} must be categorical. Default is \code{center_spec = NULL}. See also \sQuote{Details} and \sQuote{Examples}.}

\item{theta_A_polys}{a list with \eqn{L} elements so that each element is the array \code{theta_A_ploy} (the output of the \code{MAP.estimation} function, \code{MAP.estimation()$theta_A_ploy}) for the corresponding center. This argument, \code{theta_A_polys}, is only used if \code{family = "survival"} and \code{basehaz = "poly"}. See \sQuote{Details} and \sQuote{Examples}.}

\item{treat_round}{a character string representing the \code{"first"} and \code{"second"} rounds of estimating the treatment effect.}

\item{for_ATE}{a list of \eqn{L} vectors of 9 elements to calculate the average treatment effects (ATEs) only for the \code{binomial} and \code{gaussian} families. These vectors must have equal dimensions. If \code{treat_round = "first"}, then \code{for_ATE} must be \code{NULL}. If \code{treat_round = "second"}, then \code{for_ATE} must be a list for \code{binomial} and \code{gaussian}, while for \code{survival}, \code{for_ATE} must be \code{NULL}. It should be defined using the output of \code{MAP.estimation()$for_ATE} obtained from the first round. See \sQuote{Details} and \sQuote{Examples}.}

\item{p}{an integer representing the number of covariates/coefficients. It can be found from the output of the \code{MAP.estimation} function, \code{MAP.estimation()$np}). This argument, \code{p}, is only used if \code{stratified = TRUE} and \code{family = "survival"}.}

\item{q_ls}{a vector with \eqn{L} elements in which each element is the order (minus 1) of the exponentiated polynomial baseline hazard function for the corresponding center, i.e., each element is the value of \code{q_l} (the output of the \code{MAP.estimation} function, \code{MAP.estimation()$q_l}). This argument, \code{q_ls}, is only used if \code{family = "survival"}, \code{family = "survival"} and \code{basehaz = "poly"}. It can also be a scalar which represents the maximum value of the \code{q_l}'s across the centers.}

\item{center_zero_sample}{logical flag indicating whether the center has a categorical covariate with no observations/individuals in one of the categories. It is used to address heterogeneity across centers due to center-specific covariates. Default is \code{center_zero_sample = FALSE}. For more detailes see \sQuote{References}.}

\item{which_cent_zeros}{an integer vector representing the center(s) which has one categorical covariate with no individuals in one of the categories. It is used if \code{center_zero_sample = TRUE}.}

\item{zero_sample_covs}{a vector in which each element is a character string representing the categorical covariate that has no samples/observations in one of its categories for the corresponding center. Each element of the vector can be obtained from the output of the \code{MAP.estimation} function for the corresponding center, \code{MAP.estimation()$zero_sample_cov}. It is used when \code{center_zero_sample = TRUE}.}

\item{refer_cats}{a vector in which each element is a character string representing the reference category for the corresponding center. Each element of the vector can be obtained from the output of the \code{MAP.estimation} function for the corresponding center, \code{MAP.estimation()$refer_cat}. This vector is used when \code{center_zero_sample = TRUE}.}

\item{zero_cats}{a vector in which each element is a character string representing the category with no samples/observations for the corresponding center. Each element of the vector can be obtained from the output of the \code{MAP.estimation} function for the corresponding center, i.e., \code{MAP.estimation()$zero_cat}. It is used when \code{center_zero_sample = TRUE}.}

\item{lev_no_ref_zeros}{a list in which the number of elements equals the length of the \code{which_cent_zeros} argument. Each element of the list is a vector containing the names of the levels of the categorical covariate that has no samples/observations in one of its categories for the corresponding center. However,  the name of the category with no samples and the name of the reference category are excluded from this vector. Each element of the list can be obtained from the output of the \code{MAP.estimation} function, i.e., \code{MAP.estimation()$lev_no_ref_zero}.
This argument is used if \code{center_zero_sample = TRUE}.}

}

\value{

\code{bfi} returns a list containing the following components:

\item{theta_hat}{the vector of estimates obtained by combining the inference results from the \eqn{L} centers with the \code{'BFI'} methodology. If an intercept was fitted in every center and \code{stratified = FALSE}, there is only one general ``intercept'' in this vector, while if \code{stratified = TRUE} and \code{strat_par = 1}, there are \eqn{L} different intercepts in the model, for each center one. If \code{treatment} is not '\code{NULL}', when \code{treat_round = 'first'}, \code{theta_hat} gives \eqn{\hat{\boldsymbol{\gamma}}_{BFI}}, and when \code{treat_round = 'second'}, \code{theta_hat} is the treatment effect \eqn{\hat \zeta_{BFI}};}

\item{A_hat}{minus the curvature (or Hessian) matrix obtained by the \code{'BFI'} method for the combined model. If \code{stratified = TRUE}, the dimension of the matrix is always greater than when \code{stratified = FALSE};}

\item{sd}{the vector of (posterior) standard deviation of the estimates in \code{theta_hat} obtained from the matrix in \code{A_hat}, i.e., the vector equals \code{sqrt(diag(solve(A_hat)))} which equals the square root of the elements at the diagonal of the inverse of the \code{A_hat} matrix.}

\item{family}{the \code{family} object used;}

\item{basehaz}{the baseline hazard function used;}

\item{stratified}{whether a stratified analysis was done or not;}

\item{strat_par}{the stratification parameter(s) used;}

%\item{S_var}{the sample variance of the treatment and control arms. It will be 'NULL' for \code{survival}.}

\item{Ave_Treat}{the estimates of the average treatment effect. Two diffterent estimations (IPTW and wIPTW) if the family is \code{gaussian} or \code{binomial}, and for the \code{survival} family it is 'NULL'. For more detailes see \sQuote{References}.}
}

\details{
\code{bfi} function implements the BFI approach described in the papers Jonker et. al. (2024a), Pazira et. al. (2024) and Jonker et. al. (2024b) given in the references.
The inference results gathered from different (\eqn{L}) centers are combined, and the BFI estimates of the model parameters and curvature matrix evaluated at that point are returned.

The inference result from each center must be obtained using the \code{MAP.estimation} function separately, and then all of these results (coming from different centers) should be compiled into a list to be used as an input of \code{bfi()}.
The models in the different centers should be defined in exactly the same way; among others, exactly the same covariates should be included in the models. The parameter vectors should be defined exactly the same, so that the \eqn{L} vectors and matrices in the input lists \code{theta_hat}'s and \code{A_hat}'s are defined in the same way (e.g., the covariates need to be included in the models in the same order).

Note that the order of the elements in the lists \code{theta_hats}, \code{A_hats} and \code{Lambda}, must be the same with respect to the centers, so that in every list the element at the \eqn{\ell^{\th}} position is from the center \eqn{\ell}. This should also be the case for the vector \code{center_spec}.

If for the locations \code{intercept = FALSE}, the stratified analysis is not possible anymore for the \code{binomial} family.

If \code{stratified = FALSE}, both \code{strat_par} and \code{center_spec} must be \code{NULL} (the defaults), while if \code{stratified = TRUE} only one of the two must be \code{NULL}.

If \code{stratified = FALSE} and all the \eqn{L+1} matrices in \code{Lambda} are equal, it is sufficient to give a (list of) one matrix only.
In both cases of the \code{stratified} argument (\code{TRUE} or \code{FALSE}), if only the first \eqn{L} matrices are equal, the argument \code{Lambda} can be a list of two matrices, so that the fist matrix represents the chosen variance-covariance matrix for local centers and the second one is the chosen matrix for the combined data set.
The last matrix of the list in the argument \code{Lambda} can be built by the function \code{inv.prior.cov()}.

If the data type used in the argument \code{center_spec} is continuous or categorical with the number of categories equal to the number of centers, one can use \code{stratified = TRUE} and \code{center_spec = NULL}, and set \code{strat_par} not to \code{NULL} (i.e., to \eqn{1}, \eqn{2} or both \eqn{(1, 2)}). Indeed, in this case, the stratification parameter(s) given in the argument \code{strat_par} are assumed to be different across the centers.

When \code{family = 'survival'} and \code{basehaz = 'poly'}, the arguments \code{theta_hats} and \code{A_hats} should not be provided. Instead, the \code{theta_A_polys} and \code{q_ls} arguments should be defined using the local information, specifically \code{MAP.estimation()$theta_A_poly} and \code{MAP.estimation()$q_l}, respectively. See Example 3 in \sQuote{Examples}.

For estimating the treatment effect, in the first round (\code{treat_round = "first"}), the argument \code{for_ATE} must be \code{NULL} (the default) and the family must be set to \code{binomial} (family is handled automatically.)
}

\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2024a). \emph{Bayesian federated inference for estimating statistical models based on non-shared multicenter data sets}, \emph{Statistics in Medicine}, 43(12): 2421-2438. <https://doi.org/10.1002/sim.10072>

Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker M.A. (2025b). \emph{Bayesian Federated Inference for Survival Models}, \emph{Journal of Applied Statistics (Accepted)}. <https://arxiv.org/abs/2404.17464>

Jonker M.A., Pazira H. and Coolen A.C.C. (2025a). \emph{Bayesian Federated Inference for regression models based on non-shared medical center data}, \emph{Research Synthesis Methods}, 1-41. <https://doi.org/10.1017/rsm.2025.6>

%van den Heuvel Z.D., Pazira H. and Jonker M.A. (2025c). \emph{Bayesian Federated Causal Inference from observational data}, \emph{arXiv}. <https://arxiv.org/abs/???>
}

\author{Hassan Pazira and Marianne Jonker \cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{MAP.estimation}} and \code{\link{inv.prior.cov}}
}

\examples{
#################################################
##  Example 1:  y ~ Binomial  (L = 2 centers)  ##
#################################################

# Setting a seed for reproducibility
set.seed(112358)

#------------------------------------#
# Data Simulation for Local Center 1 #
#------------------------------------#
n1 <- 30                                           # sample size of center 1
X1 <- data.frame(x1=rnorm(n1),                     # continuous variable
                 x2=sample(0:2, n1, replace=TRUE)) # categorical variable
# make dummy variables
X1x2_1 <- ifelse(X1$x2 == '1', 1, 0)
X1x2_2 <- ifelse(X1$x2 == '2', 1, 0)
X1$x2  <- as.factor(X1$x2)
# regression coefficients
beta <- 1:4  # beta[1] is the intercept
# linear predictor:
eta1   <- beta[1] + X1$x1 * beta[2] + X1x2_1 * beta[3] + X1x2_2 * beta[4]
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu1    <- binomial()$linkinv(eta1)
y1     <- rbinom(n1, 1, mu1)

#------------------------------------#
# Data Simulation for Local Center 2 #
#------------------------------------#
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

#---------------------------#
# MAP Estimates at Center 1 #
#---------------------------#
# Assume the same inverse covariance matrix (Lambda) for both centers:
Lambda     <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial')
fit1       <- MAP.estimation(y1, X1, family = 'binomial', Lambda)
theta_hat1 <- fit1$theta_hat # intercept and coefficient estimates
A_hat1     <- fit1$A_hat     # minus the curvature matrix

#---------------------------#
# MAP Estimates at Center 2 #
#---------------------------#
fit2       <- MAP.estimation(y2, X2, family='binomial', Lambda)
theta_hat2 <- fit2$theta_hat
A_hat2     <- fit2$A_hat

#-----------------------#
# BFI at Central Server #
#-----------------------#
theta_hats <- list(theta_hat1, theta_hat2)
A_hats     <- list(A_hat1, A_hat2)
bfi        <- bfi(theta_hats, A_hats, Lambda, family='binomial')
class(bfi)
summary(bfi, cur_mat=TRUE)

###---------------------###
### Stratified Analysis ###
###---------------------###

# By running the following line an error appears because
# when stratified = TRUE, both 'strat_par' and 'center_spec' can not be NULL:
Just4check1 <- try(bfi(theta_hats, A_hats, Lambda, family = 'binomial',
                   stratified = TRUE), TRUE)
class(Just4check1) # By default, both 'strat_par' and 'center_spec' are NULL!

# By running the following line an error appears because when stratified = TRUE,
# last matrix in 'Lambda' should not have the same dim. as the other local matrices:
Just4check2 <- try(bfi(theta_hats, A_hats, Lambda, stratified = TRUE,
                   strat_par = 1), TRUE)
class(Just4check2) # All matices in Lambda have the same dimension!

# Stratified analysis when 'intercept' varies across two centers:
newLam <- inv.prior.cov(X1, lambda=c(0.1, 0.3), family = 'binomial',
                        stratified = TRUE, strat_par = 1)
bfi <- bfi(theta_hats, A_hats, list(Lambda, newLam), family = 'binomial',
           stratified=TRUE, strat_par=1)
summary(bfi, cur_mat=TRUE)


###---------------------###
###  Treatment Effect   ###
###---------------------###

set.seed(112358)

#------------------------------------#
# Data Simulation for Local Center 1 #
#------------------------------------#
n1 <- 30                                           # sample size of center 1
X1 <- data.frame(x1=rnorm(n1),                     # continuous variable
                 treatment=sample(1:2, n1, replace=TRUE)) # categorical variable
X1$treatment  <- as.factor(X1$treatment)

# regression coefficients
beta <- 1:3  # beta[1] is the intercept
# make dummy variable
X1x2_2 <- ifelse(X1$treatment == '2', 1, 0)
# linear predictor:
eta1   <- beta[1] + X1$x1 * beta[2] + X1x2_2 * beta[3]
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu1    <- binomial()$linkinv(eta1)
y1     <- rbinom(n1, 1, mu1)

#------------------------------------#
# Data Simulation for Local Center 2 #
#------------------------------------#
n2 <- 50                                           # sample size of center 2
X2 <- data.frame(x1=rnorm(n2),                     # continuous variable
                 treatment=sample(1:2, n2, replace=TRUE)) # categorical variable
X2$treatment  <- as.factor(X2$treatment)
# make dummy variables:
X2x2_2 <- ifelse(X2$treatment == '2', 1, 0)
# linear predictor:
eta2   <- beta[1] + X2$x1 * beta[2] + X2x2_2 * beta[3]
# inverse of the link function:
mu2    <- binomial()$linkinv(eta2)
y2     <- rbinom(n2, 1, mu2)

# The algorithm works even if the order of the covariates are not
# the same across centers
X2 <- X2[,c("treatment","x1")]

#-----------------------#
#  Observational data   #
#-----------------------#

# For observational data (RWD), we need two rounds for estimating treatment effect:

#-------------#
# First Round #
#-------------#

## Center 1:
Lambda1 <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial',
                         treatment = "treatment", treat_round="first")
fit1_r1 <- MAP.estimation(y1, X1, family = 'binomial', Lambda = Lambda1,
                          treatment = "treatment", treat_round = "first")
# In the first round, the output is without the treatment!
summary(fit1_r1)

## Center 2:
Lambda2 <- inv.prior.cov(X2, lambda = 0.01, family = 'binomial',
                         treatment = "treatment", treat_round="first")
fit2_r1 <- MAP.estimation(y2, X2, family = 'binomial', Lambda = Lambda2,
                          treatment = "treatment", treat_round = "first")
fit2_r1

## Centeral Server:
theta_hats_r1 <- list(fit1_r1$theta_hat, fit2_r1$theta_hat)
A_hats_r1 <- list(fit1_r1$A_hat, fit2_r1$A_hat)
fitbfi_r1 <- bfi(theta_hats_r1, A_hats_r1, Lambda1, family = 'binomial',
                 treat_round = "first")
summary(fitbfi_r1, cur_mat = TRUE)

#--------------#
# Second Round #
#--------------#

## Center 1:
Lambda11 <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial',
                          treatment = "treatment", treat_round="second")
fit1_r2 <- MAP.estimation(y1, X1, family = 'binomial', Lambda = Lambda11,
                          treatment = "treatment", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)
# In the second round, the output is only with the treatment!
summary(fit1_r2)

## Center 2:
Lambda22 <- inv.prior.cov(X2, lambda = 0.01, family = 'binomial',
                         treatment = "treatment", treat_round="second")
fit2_r2 <- MAP.estimation(y2, X2, family = 'binomial', Lambda = Lambda22,
                          treatment = "treatment", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)
fit2_r2$propensity # Propensity Score
fit2_r2$for_ATE # will be used in central server
fit2_r2

## Centeral Server:
theta_hats_r2 <- list(fit1_r2$theta_hat, fit2_r2$theta_hat)
A_hats_r2 <- list(fit1_r2$A_hat, fit2_r2$A_hat)
for_ATEs <- list(fit1_r2$for_ATE, fit2_r2$for_ATE)
fitbfi_r2 <- bfi(theta_hats_r2, A_hats_r2, Lambda11, family = 'binomial',
                 treat_round = "second", for_ATE = for_ATEs)
fitbfi_r2$S_var
fitbfi_r2$Ave_Treat
summary(fitbfi_r2)

#--------------------#
#  Randomized Trial  #
#--------------------#

# For Randomized Control Trial (RCT), we need only one round (the second round) for
# estimating treatment effect. Because we do not need to estimate propensity score.
# For example, in a 1:1 randomized trial, the propensity scores are, by definition,
# equal to 0.5. Here we use 'RCT_propens', instead of 'gamma_bfi':

## Center 1:
Lambda11 <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial',
                          treatment = "treatment", treat_round="second")
fit1_r2 <- MAP.estimation(y1, X1, family = 'binomial', Lambda = Lambda11,
                          treatment = "treatment", treat_round = "second",
                          RCT_propens = rep(0.5, n1)) # gamma_bfi = NULL
summary(fit1_r2)

## Center 2:
Lambda22 <- inv.prior.cov(X2, lambda = 0.01, family = 'binomial',
                          treatment = "treatment", treat_round="second")
fit2_r2 <- MAP.estimation(y2, X2, family = 'binomial', Lambda = Lambda22,
                          treatment = "treatment", treat_round = "second",
                          RCT_propens = rep(0.5, n2)) # gamma_bfi = NULL
fit2_r2$for_ATE # will be used in central server
fit2_r2

## Centeral Server:
theta_hats_r2 <- list(fit1_r2$theta_hat, fit2_r2$theta_hat)
A_hats_r2 <- list(fit1_r2$A_hat, fit2_r2$A_hat)
for_ATEs <- list(fit1_r2$for_ATE, fit2_r2$for_ATE)
fitbfi_r2 <- bfi(theta_hats_r2, A_hats_r2, Lambda11, family = 'binomial',
                 treat_round = "second", for_ATE = for_ATEs)
fitbfi_r2$S_var
fitbfi_r2$Ave_Treat
summary(fitbfi_r2)


#################################################
##  Example 2:  y ~ Gaussian  (L = 3 centers)  ##
#################################################

# Setting a seed for reproducibility
set.seed(112358)

p     <- 3                     # number of coefficients without 'intercept'
theta <- c(1, rep(2, p), 1.5)  # reg. coef.s ('intercept' is 1) & 'sigma2' = 1.5

#------------------------------------#
# Data Simulation for Local Center 1 #
#------------------------------------#
n1   <- 30                                       # sample size of center 1
X1   <- data.frame(matrix(rnorm(n1 * p), n1, p)) # continuous variables
# linear predictor:
eta1 <- theta[1] + as.matrix(X1) %*% theta[2:4]
# inverse of the link function ( g^{-1}(\eta) = \mu ):
mu1  <- gaussian()$linkinv(eta1)
y1   <- rnorm(n1, mu1, sd = sqrt(theta[5]))

#------------------------------------#
# Data Simulation for Local Center 2 #
#------------------------------------#
n2   <- 40                                       # sample size of center 2
X2   <- data.frame(matrix(rnorm(n2 * p), n2, p)) # continuous variables
# linear predictor:
eta2 <- theta[1] + as.matrix(X2) %*% theta[2:4]
# inverse of the link function:
mu2  <- gaussian()$linkinv(eta2)
y2   <- rnorm(n2, mu2, sd = sqrt(theta[5]))

#------------------------------------#
# Data Simulation for Local Center 3 #
#------------------------------------#
n3   <- 50                                       # sample size of center 3
X3   <- data.frame(matrix(rnorm(n3 * p), n3, p)) # continuous variables
# linear predictor:
eta3 <- theta[1] + as.matrix(X3) %*% theta[2:4]
# inverse of the link function:
mu3  <- gaussian()$linkinv(eta3)
y3   <- rnorm(n3, mu3, sd = sqrt(theta[5]))

#---------------------------#
# Inverse Covariance Matrix #
#---------------------------#
# Creating the inverse covariance matrix for the Gaussian prior distribution:
# the same for both centers
Lambda <- inv.prior.cov(X1, lambda = 0.05, family='gaussian')

#---------------------------#
# MAP Estimates at Center 1 #
#---------------------------#
fit1       <- MAP.estimation(y1, X1, family = 'gaussian', Lambda)
theta_hat1 <- fit1$theta_hat # intercept and coefficient estimates
A_hat1     <- fit1$A_hat     # minus the curvature matrix

#---------------------------#
# MAP Estimates at Center 2 #
#---------------------------#
fit2       <- MAP.estimation(y2, X2, family = 'gaussian', Lambda)
theta_hat2 <- fit2$theta_hat
A_hat2     <- fit2$A_hat

#---------------------------#
# MAP Estimates at Center 3 #
#---------------------------#
fit3       <- MAP.estimation(y3, X3, family = 'gaussian', Lambda)
theta_hat3 <- fit3$theta_hat
A_hat3     <- fit3$A_hat

#-----------------------#
# BFI at Central Server #
#-----------------------#
A_hats     <- list(A_hat1, A_hat2, A_hat3)
theta_hats <- list(theta_hat1, theta_hat2, theta_hat3)
bfi        <- bfi(theta_hats, A_hats, Lambda, family = 'gaussian')
summary(bfi, cur_mat=TRUE)

###---------------------###
### Stratified Analysis ###
###---------------------###

# Stratified analysis when 'intercept' varies across two centers:
newLam1 <- inv.prior.cov(X1, lambda = c(0.1,0.3), family = 'gaussian',
                         stratified = TRUE, strat_par = 1, L = 3)
# 'newLam1' is used as the prior for combined data and
# 'Lambda' is used as the prior for locals
list_newLam1 <- list(Lambda, newLam1)
bfi1 <- bfi(theta_hats, A_hats, list_newLam1, family = 'gaussian',
            stratified = TRUE, strat_par = 1)
summary(bfi1, cur_mat = TRUE)

# Stratified analysis when 'sigma2' varies across two centers:
newLam2 <- inv.prior.cov(X1, lambda = c(0.1,0.3), family = 'gaussian',
                         stratified = TRUE, strat_par = 2, L = 3)
# 'newLam2' is used as the prior for combined data and 'Lambda' is used as
# the prior for locals
list_newLam2 <- list(Lambda, newLam2)
bfi2 <- bfi(theta_hats, A_hats, list_newLam2, family = 'gaussian',
            stratified = TRUE, strat_par=2)
summary(bfi2, cur_mat = TRUE)

# Stratified analysis when 'intercept' and 'sigma2' vary across 2 centers:
newLam3 <- inv.prior.cov(X1, lambda = c(0.1,0.2,0.3), family = 'gaussian',
                         stratified = TRUE, strat_par = c(1, 2), L = 3)
# 'newLam3' is used as the prior for combined data and 'Lambda' is used as
# the prior for locals
list_newLam3 <- list(Lambda, newLam3)
bfi3 <- bfi(theta_hats, A_hats, list_newLam3, family = 'gaussian',
            stratified = TRUE, strat_par = 1:2)
summary(bfi3, cur_mat = TRUE)

###----------------------------###
### Center Specific Covariates ###
###----------------------------###

# Assume the first and third centers have the same center-specific covariate value
# of 'High', while this value for the second center is 'Low', i.e.,
# center_spec = c('High','Low','High')
newLam4 <- inv.prior.cov(X1, lambda=c(0.1, 0.2, 0.3), family='gaussian',
                         stratified = TRUE, center_spec = c('High','Low','High'),
                         L = 3)
# 'newLam4' is used as the prior for combined data and 'Lambda' is used as
# the prior for locals
l_newLam4 <- list(Lambda, newLam4)
bfi4 <- bfi(theta_hats, A_hats, l_newLam4, family = 'gaussian',
            stratified = TRUE, center_spec = c('High','Low','High'))
summary(bfi4, cur_mat = TRUE)


###---------------------###
###  Treatment Effect   ###
###---------------------###

set.seed(112358)

#-----------------------------#
# New Data for Local Center 1 #
#-----------------------------#
# Generating new data with 'treatment' variable
# We cansider the first variable (X1$X1) to be the treatment:
X1$X1 <- sample(0:1, n1, replace=TRUE) # categorical variable
eta1  <- theta[1] + as.matrix(X1) %*% theta[2:4]
mu1   <- gaussian()$linkinv(eta1)
y1    <- rnorm(n1, mu1, sd = sqrt(theta[5]))

#-----------------------------#
# New Data for Local Center 2 #
#-----------------------------#
# We cansider the first variable (X2$X1) to be the treatment:
X2$X1 <- sample(0:1, n2, replace=TRUE) # categorical variable
eta2  <- theta[1] + as.matrix(X2) %*% theta[2:4]
mu2   <- gaussian()$linkinv(eta2)
y2    <- rnorm(n2, mu2, sd = sqrt(theta[5]))

#-----------------------------#
# New Data for Local Center 3 #
#-----------------------------#
# We cansider the first variable (X3$X1) to be the treatment:
X3$X1 <- sample(0:1, n3, replace=TRUE) # categorical variable
# linear predictor:
eta3  <- theta[1] + as.matrix(X3) %*% theta[2:4]
# inverse of the link function:
mu3   <- gaussian()$linkinv(eta3)
y3    <- rnorm(n3, mu3, sd = sqrt(theta[5]))

#-----------------------#
#  Observational data   #
#-----------------------#

# For observational data (RWD), we need two rounds for estimating treatment effect:

#-------------#
# First Round #
#-------------#

## Center 1:
Lambda1 <- inv.prior.cov(X1, lambda = 0.01, family = 'binomial',
                         treatment = "X1", treat_round="first")
# When treat_round = "first", the family will automatically set to 'binomial',
# even if family = 'gaussian' or family = 'survival'.
fit1_r1 <- MAP.estimation(y1, X1, family = 'gaussian', Lambda = Lambda1,
                          treatment = "X1", treat_round = "first")
# Althghou family = 'gaussian', the output is based on 'binomial'!
# The output without the treatment (X1) in the first round!
summary(fit1_r1)

## Center 2:
Lambda2 <- inv.prior.cov(X2, lambda = 0.01, family = 'gaussian',
                         treatment = "X1", treat_round="first")
fit2_r1 <- MAP.estimation(y2, X2, family = 'gaussian', Lambda = Lambda2,
                          treatment = "X1", treat_round = "first")
fit2_r1

## Center 3:
Lambda3 <- inv.prior.cov(X3, lambda = 0.01, family = 'gaussian',
                         treatment = "X1", treat_round="first")
fit3_r1 <- MAP.estimation(y3, X3, family = 'gaussian', Lambda = Lambda3,
                          treatment = "X1", treat_round = "first")

## Centeral Server:
theta_hats_r1 <- list(fit1_r1$theta_hat, fit2_r1$theta_hat, fit3_r1$theta_hat)
A_hats_r1 <- list(fit1_r1$A_hat, fit2_r1$A_hat, fit3_r1$A_hat)
fitbfi_r1 <- bfi(theta_hats_r1, A_hats_r1, Lambda1, family = 'gaussian',
                 treat_round = "first") # same results with 'binomial'
# The output without the treatment (X1) in the first round!
summary(fitbfi_r1, cur_mat = TRUE)

#--------------#
# Second Round #
#--------------#

## Center 1:
Lambda11 <- inv.prior.cov(X1, lambda = 0.01, family = 'gaussian',
                          treatment = "X1", treat_round="second")
fit1_r2 <- MAP.estimation(y1, X1, family = 'gaussian', Lambda = Lambda11,
                          treatment = "X1", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)
# The output with only the treatment (X1) in the second round!
summary(fit1_r2)

## Center 2:
Lambda22 <- inv.prior.cov(X2, lambda = 0.01, family = 'gaussian', treatment = "X1",
                          treat_round="second")
fit2_r2 <- MAP.estimation(y2, X2, family = 'gaussian', Lambda = Lambda22,
                          treatment = "X1", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)

## Center 3:
Lambda33 <- inv.prior.cov(X3, lambda = 0.01, family = 'gaussian', treatment = "X1",
                          treat_round="second")
fit3_r2 <- MAP.estimation(y3, X3, family = 'gaussian', Lambda = Lambda33,
                          treatment = "X1", treat_round = "second",
                          gamma_bfi = fitbfi_r1$theta_hat)

## Centeral Server:
theta_hats_r2 <- list(fit1_r2$theta_hat, fit2_r2$theta_hat, fit3_r2$theta_hat)
A_hats_r2 <- list(fit1_r2$A_hat, fit2_r2$A_hat, fit3_r2$A_hat)
for_ATEs <- list(fit1_r2$for_ATE, fit2_r2$for_ATE, fit3_r2$for_ATE)
fitbfi_r2 <- bfi(theta_hats_r2, A_hats_r2, Lambda11, family = 'gaussian',
                 treat_round = "second", for_ATE = for_ATEs)
fitbfi_r2$Ave_Treat
fitbfi_r2$S_var
summary(fitbfi_r2)


####################################################
##  Example 3:  Survival family  (L = 2 centers)  ##
####################################################

# Setting a seed for reproducibility
set.seed(112358)

p <- 3
theta <- c(1:4, 5, 6)  # regression coefficients (1:4) & omega's (5:6)

#---------------------------------------------#
# Simulating Survival data for Local Center 1 #
#---------------------------------------------#
n1 <- 30
X1 <- data.frame(matrix(rnorm(n1 * p), n1, p)) # continuous (normal) variables
# Simulating survival data ('time' and 'status') from 'Weibull' with
# a predefined censoring rate of 0.3:
y1 <- surv.simulate(Z = list(X1), beta = theta[1:p], a = theta[5],
                    b = theta[6], u1 = 0.1, cen_rate = 0.3,
                    gen_data_from = "weibul")$D[[1]][, 1:2]

## MAP Estimates at Center 1
Lambda <- inv.prior.cov(X1, lambda = c(0.1, 1), family = "survival",
                        basehaz = "poly")
fit1 <- MAP.estimation(y1, X1, family = 'survival', Lambda = Lambda,
                       basehaz = "poly")
theta_hat1 <- fit1$theta_hat  # coefficient estimates
A_hat1     <- fit1$A_hat      # minus the curvature matrix
summary(fit1, cur_mat=TRUE)
fit1$theta_A_poly # Only when family = "survival" and basehaz ="poly"

#---------------------------------------------#
# Simulating Survival data for Local Center 2 #
#---------------------------------------------#
n2 <- 30
X2 <- data.frame(matrix(rnorm(n2 * p), n2, p)) # continuous (normal) variables
# Survival simulated data from 'Weibull' with a predefined censoring rate of 0.3:
y2 <- surv.simulate(Z = list(X2), beta = theta[1:p], a = theta[5],
                    b = theta[6],u1 = 0.1, cen_rate = 0.3,
                    gen_data_from = "weibul")$D[[1]][, 1:2]

## MAP Estimates at Center 2
fit2 <- MAP.estimation(y2, X2, family = 'survival', Lambda = Lambda,
                       basehaz = "poly")
theta_hat2 <- fit2$theta_hat
A_hat2 <- fit2$A_hat
summary(fit2, cur_mat=TRUE)

#-----------------------#
# BFI at Central Server #
#-----------------------#
# When family = 'survival' and basehaz = "poly", only 'theta_A_polys'
# should be defined instead of 'theta_hats' and 'A_hats':
theta_A_hats <- list(fit1$theta_A_poly, fit2$theta_A_poly)
qls <- c(fit1$q_l, fit2$q_l)
bfi <- bfi(Lambda = Lambda, family = 'survival', theta_A_polys = theta_A_hats,
           basehaz = "poly", q_ls = qls)
summary(bfi, cur_mat=TRUE)


###---------------------###
### Stratified Analysis ###
###---------------------###

# Stratified analysis when first parameter ('omega_0') varies across two centers:
(newLam0 <- inv.prior.cov(X1, lambda = c(rep(1, 3), 0.3, 0.7, rep(2,2)),
                          family = 'survival', stratified = TRUE,
                          basehaz = c("poly"), strat_par = 1, L = 2))
# 'newLam0' is used as the prior for combined data and 'Lambda' is used as for locals:
list_newLam0 <- list(Lambda, newLam0)
bfi0 <- bfi(Lambda = list_newLam0, family = 'survival', theta_A_polys = theta_A_hats,
            stratified = TRUE, basehaz = c("poly"), p = 3, q_ls = qls, strat_par = 1)
summary(bfi0, cur_mat = TRUE)


# Stratified analysis when the first and second parameters ('omega_0' and 'omega_1')
# vary across two centers:
newLam1 <- inv.prior.cov(X1, lambda = c(rep(1, 3), 0.3, 0.7, 0.5, 0.8, 2),
                         family = 'survival', stratified = TRUE, basehaz = c("poly"),
                         strat_par = c(1, 2), L = 2)
# 'newLam1' is used as the prior for combined data:
list_newLam1 <- list(Lambda, newLam1)
bfi1 <- bfi(Lambda = list_newLam1, family = 'survival', theta_A_polys = theta_A_hats,
            stratified = TRUE, basehaz = c("poly"), p = 3, q_ls = qls,
            strat_par = c(1, 2))
summary(bfi1, cur_mat = TRUE)


###---------------------###
###  Treatment Effect   ###
###---------------------###

set.seed(112358)

#-----------------------------#
# New Data for Local Center 1 #
#-----------------------------#
# Generating new data with 'treatment' variable
# We cansider the first variable (X1$X1) to be the treatment
X1$X1 <- sample(0:1, n1, replace=TRUE) # categorical variable
y1 <- surv.simulate(Z = list(X1), beta = theta[1:p], a = theta[5], b = theta[6],
                    u1 = 0.1, cen_rate = 0.3, gen_data_from = "weibul")$D[[1]][, 1:2]

#-----------------------------#
# New Data for Local Center 2 #
#-----------------------------#
# We cansider the first variable (X2$X1) to be the treatment!
X2$X1 <- sample(0:1, n2, replace=TRUE) # categorical variable
y2 <- surv.simulate(Z = list(X2), beta = theta[1:p], a = theta[5], b = theta[6],
                    u1 = 0.1, cen_rate = 0.3, gen_data_from = "weibul")$D[[1]][, 1:2]

#-------------#
# First Round #
#-------------#

## Center 1:
Lambda1 <- inv.prior.cov(X1, lambda = 0.01, family = 'survival',
                         treatment = "X1", treat_round="first")
# When treat_round = "first", the family will automatically set to 'binomial',
# even if family = 'gaussian' or family = 'survival'.
fit1_r1 <- MAP.estimation(y1, X1, family = 'survival', # 'basehaz' is not needed!
                          Lambda = Lambda1, treatment = "X1", treat_round = "first")
# While family = 'survival', the output is based on 'binomial' with no 'Intercept'!
# The output without the treatment (X1) in the first round!
summary(fit1_r1)

## Center 2:
Lambda2 <- inv.prior.cov(X2, lambda = 0.01, family = 'survival',
                         treatment = "X1", treat_round="first")
fit2_r1 <- MAP.estimation(y2, X2, family = 'survival', Lambda = Lambda2,
                          treatment = "X1", treat_round = "first")
fit2_r1

## Centeral Server:
theta_hats_r1 <- list(fit1_r1$theta_hat, fit2_r1$theta_hat)
A_hats_r1 <- list(fit1_r1$A_hat, fit2_r1$A_hat)
fitbfi_r1 <- bfi(theta_hats_r1, A_hats_r1, Lambda1, family = 'survival',
                 treat_round = "first")
# In the first round output is based on 'binomial', and without
# the intercept and treatment (X1):
summary(fitbfi_r1, cur_mat = TRUE)

#--------------#
# Second Round #
#--------------#

## Center 1:
Lambda11 <- inv.prior.cov(X1, lambda = 0.01, family = 'survival',
                          basehaz = "unspecified", treatment = "X1",
                          treat_round="second")
fit1_r2 <- MAP.estimation(y1, X1, family = 'survival', Lambda = Lambda11,
                          basehaz = "unspecified", treatment = "X1",
                          treat_round = "second", gamma_bfi = fitbfi_r1$theta_hat)
# The output with only the treatment (X1) in the second round!
summary(fit1_r2)

## Center 2:
Lambda22 <- inv.prior.cov(X2, lambda = 0.01, family = 'survival',
                          basehaz = "unspecified", treatment = "X1",
                          treat_round="second")
fit2_r2 <- MAP.estimation(y2, X2, family = 'survival', basehaz = "unspecified",
                          Lambda = Lambda22, treatment = "X1",
                          treat_round = "second", gamma_bfi = fitbfi_r1$theta_hat)
fit2_r2

## Centeral Server:
theta_hats_r2 <- list(fit1_r2$theta_hat, fit2_r2$theta_hat)
A_hats_r2 <- list(fit1_r2$A_hat, fit2_r2$A_hat)
fitbfi_r2 <- bfi(theta_hats_r2, A_hats_r2, Lambda11, family = 'survival',
                 basehaz = "unspecified", treat_round = "second")
# When family = 'survival', 'for_ATE' is not calculated.
summary(fitbfi_r2)


}
\keyword{regression}
\keyword{models}
\keyword{survival}
\keyword{nonparametric}
\keyword{bayesian}
\keyword{federated}
