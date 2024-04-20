\name{hazards.fun}
\alias{hazards.fun}
\title{Compute the estimated (baseline/cumulative) hazard and (baseline) survival functions}
\description{
For a given vector of times, \code{hazards.fun} computes the estimated baseline hazard, cumulative baseline hazard, hazard, baseline survival, and survival functions. It can be used for prediction on a new sample.
}
\usage{
hazards.fun(time,
            z = NULL,
            p,
            theta_hat,
            basehaz = c("weibul", "exp", "gomp", "poly", "pwexp"),
            q_max,
            timax)
}

\arguments{
\item{time}{the vector containing the time values for which the hazard rate is computed. If the argument \code{z} is not \code{NULL}, then the length of the argument \code{time} should be the number of columns of \code{z}, which is \eqn{p}.}

\item{z}{a new observation vector of length \eqn{p}. If \code{z = NULL} (the default), then the relative risk (\eqn{\boldsymbol{z}^{\top} \boldsymbol{\beta}}) is considered a vector of 1 with length \eqn{n}.}

\item{p}{the number of coefficients. It is taken equal to the number of elements of the argument \code{z}, if \code{z} is not \code{NULL}.}

\item{theta_hat}{a vector contains the values of the estimated parameters. The first \eqn{p} values represent the coefficient parameters (\eqn{\boldsymbol{\beta}}), while the remaining values pertain to the parameters of the baseline hazard function  (\eqn{\boldsymbol{\omega}}).}

\item{basehaz}{a character string representing one of the available baseline hazard functions; \code{exponential} (\code{"exp"}), \code{Weibull} (\code{"weibul"}, the default), \code{Gompertz} (\code{"gomp"}), \code{exponentiated polynomial} (\code{"poly"}), and \code{piecewise exponential} (\code{"pwexp"}). Can be abbreviated.}

\item{q_max}{a value represents the order of the exponentiated polynomial baseline hazard function. This argument should only be used when \code{basehaz = "poly"}. In the case of multiple centers, the maximum value of the orders should be used. \code{ql.LRT()} can be used for obtaining of the order of each center.}

\item{timax}{a value represents the minimum (or maximum) value of the maximum times observed in the different centers. This argument should only be used when \code{basehaz = "pwexp"}.}

}

\details{
\code{hazards.fun} computes the estimated baseline hazard, cumulative baseline hazard, hazard, baseline survival, and survival functions at different time points specified in the argument \code{time}.

The function \code{hazards.fun()} can be used for prediction purposes with new sample. The arguments \code{time} and \code{z} should be provided for the new data.

}

\value{

\code{hazards.fun} returns a list containing the following components:

\item{bhazard}{the vector of estimates of the baseline hazard function at the time points given by the argument \code{time};}

\item{cbhazard}{the vector of estimates of the cumulative baseline hazard function at the time points specified in the argument \code{time};}

\item{bsurvival}{the vector of estimates of the baseline survival function at the time points given by the argument \code{time};}

\item{hazard}{the vector of estimates of the hazard function at the time points given by the argument \code{time};}

\item{chazard}{the vector of estimates of the cumulative hazard function at the time points specified in the argument \code{time};}

\item{survival}{the vector of estimates of the survival function at the time points given by the argument \code{time}.}

}

\references{
Pazira H., Massa E., Weijers J.A.M., Coolen A.C.C. and Jonker M.A. (2025b). \emph{Bayesian Federated Inference for Survival Models}, \emph{Journal of Applied Statistics (Accepted)}. <https://arxiv.org/abs/2404.17464>
}

\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\seealso{
\code{\link{MAP.estimation}}
}

\examples{

# Setting a seed for reproducibility
set.seed(1123)

##-------------------------
## Simulating Survival data
##-------------------------

n <- 40
p <- 7
Original_data <- data.frame(matrix(rnorm((n+1) * p), (n+1), p))
X <- Original_data[1:n,]
X_new <- Original_data[(n+1),]
# Simulating survival data from Exponential distribution
# with a predefined censoring rate of 0.2:
Orig_y <- surv.simulate(Z = Original_data, beta = rep(1,p), a = exp(1),
                        cen_rate = 0.2, gen_data_from = "exp")$D[[1]][,1:2]
y <- Orig_y[1:n,]
y_new <- Orig_y[(n+1),]
time_points <- seq(0, max(y$time), length.out=20)

#------------------------
# Weibull baseline hazard
#------------------------

Lambda <- inv.prior.cov(X, lambda = c(0.5, 1), family = 'survival', basehaz = 'weibul')
fit_weib <- MAP.estimation(y, X, family = 'survival', Lambda = Lambda,
                           basehaz = "weibul")
# reltive risk is 1:
hazards.fun(time = time_points, p = p, theta_hat = fit_weib$theta_hat,
            basehaz = "weibul")

#-------------------------
# Gompertz baseline hazard
#-------------------------
fit_gomp <- MAP.estimation(y, X, family = 'survival', Lambda = Lambda,
                           basehaz = "gomp")
# different time points
hazards.fun(time=1:max(y*2), p = p, theta_hat = fit_gomp$theta_hat,
            basehaz = "gomp")


##----------------------------
## Prediction for a new sample
##----------------------------

## Exponentiated polynomial (poly) baseline hazard:
Lambda <- inv.prior.cov(X, lambda = c(0.5, 1), family = 'survival', basehaz = "poly")
fit_poly <- MAP.estimation(y, X, family = 'survival', Lambda = Lambda,
                           basehaz = "poly")
hazards.fun(time = y_new$time, z = X_new, theta_hat = fit_poly$theta_hat,
            basehaz = "poly", q_max = fit_poly$q_l)

## Piecewise Exponential (pwexp) baseline hazard:
Lambda <- inv.prior.cov(X, lambda = c(0.5, 1), family = 'survival', basehaz = "pwexp")
fit_pw <- MAP.estimation(y, X, family='survival', Lambda=Lambda, basehaz="pwexp",
                         min_max_times = max(y))
hazards.fun(time = y_new$time, z = X_new, theta_hat = fit_pw$theta_hat,
            basehaz = "pwexp", timax = max(y))

}
\keyword{survival}
