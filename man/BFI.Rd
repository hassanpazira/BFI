\name{bfi}
\alias{bfi}
\title{Bayesian Federated Inference}
\description{
\code{bfi} function is used (in central server) to estimate the parameters of the GLM and Survival models by BFI method using the aggregated results form the local centers.
}
\usage{
bfi(theta_hats = NULL, A_hats, family = c("binomial","gaussian"), Gamma, L = NULL, stratified = F, nuisance = 1L)
}
\arguments{
}

\details{
}

\value{
}
\references{
Jonker M.A., Pazira H. and Coolen A.C.C. (2023) <https://doi.org/10.48550/arXiv.2302.07677>
\emph{Bayesian Federated Inference for Statistical Models},
\emph{Statistics in Medicine}, Vol. 0(0), 0-0.
}
\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}}
\seealso{
\code{\link{estimators_maker}}.
}
\examples{
set.seed(1123)

#-------------
# y ~ Binomial
#-------------

L <- 2    # L is number of locations

###----------------
### Local center 1:
###----------------
n1         <- 30
p          <- 4    # p (number of coefficients) is the same for all L locations.
X1         <- data.frame(matrix(rnorm(n1 * p), n1, p))
#true_beta <- c(1, 2, 0, 0, 0)
eta1       <- 1 + 2 * X1[,1]  ## with an intercept b0=1, b1=2, b2=b3=...=bp=0
mu1        <- binomial()$linkinv(eta1)
y1         <- rbinom(n1, 1, mu1)
lambda     <- 0.01
# we assume the same (inverse) covariance matrix for all locations.
Gamma      <- Gamma.maker(X1, lambda, family="binomial")
fit1       <- MAP.estimation(y1, X1, family="binomial", Gamma)
theta_hat1 <- fit1$theta_hat # beta (intercept and coefficient) estimates
A_hat1     <- fit1$A_hat

###----------------
### Local center 2:
###----------------
n2         <- 50
p          <- 4    # p is the same for all L locations !
X2         <- data.frame(matrix(rnorm(n2 * p), n2, p))
eta2       <- 1 + 2 * X2[,1]  ## with an intercept b0=1, b1=2,b2=b3=...=bp=0
mu2        <- binomial()$linkinv(eta2)
y2         <- rbinom(n2, 1, mu2)
fit2       <- MAP.estimation(y2, X2, family="binomial", Gamma)
theta_hat2 <- fit2$theta_hat # intercept and coefficient estimates
A_hat2     <- fit2$A_hat

###--------------
### Combined data
###--------------
y              <- c(y1, y2)
X              <- rbind(X1, X2)
fit_comb       <- MAP.estimation(y, X, family="binomial", Gamma)
theta_hat_comb <- fit_comb$theta_hat # beta estimates of combined data

###-----------------------------
### Bayesian Federated Inference
###-----------------------------

A_hats <- list(A_hat1, A_hat2)
theta_hats <- list(theta_hat1, theta_hat2)

# theta (intercept and coefficient) estimates by BFI
(theta_hat_bfi <- bfi(theta_hats, A_hats, Gamma)$theta_hat)

# Curvature matrix estimate by BFI
A_bfi  <- bfi(A_hats=A_hats, Gamma=Gamma)$A_hat # == bfi(theta_hats, A_hats=A_hats, Gamma=Gamma)$A_hat

# SD of the BFI estimates
sd_bfi <- bfi(A_hats=A_hats, Gamma=Gamma)$sd # == bfi(theta_hats, A_hats=A_hats, Gamma=Gamma)$sd

### Difference between BFI estimates and estimates with combined data
theta_hat_bfi-theta_hat_comb
}
