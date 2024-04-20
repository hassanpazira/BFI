\name{surv.simulate}
\alias{surv.simulate}
\title{Generate survival data with predefined censoring rates for proportional hazards models}
\description{
\code{surv.simulate} simulates one or multiple (right-censored) survival datasets for proportional hazards models by simultaneously incorporating a baseline hazard function from three different survival distributions (exponential, Weibull and Gompertz), a random censoring time generated from a uniform distribution with an known/unknown upper limit, and a set of baseline covariates.
When the upper limit of the uniform censoring time distribution is unknown, \code{surv.simulate} can be used separately to obtain the upper limit with a predefined censoring rate.
}
\usage{
surv.simulate(L = 1,
              Z,
              beta,
              a,
              b,
              u1 = 0,
              u2,
              cen_rate,
              gen_data_from = c("exp", "weibul", "gomp"),
              only_u2 = FALSE,
              n.rep = 100,
              Trace = FALSE)
}

\arguments{
\item{L}{the number of datasets to be generated. Default is \code{L = 1}.}

\item{Z}{a list of \code{L} design matrices of dimension \eqn{n_\ell \times p}, where \eqn{n_\ell} is the number of samples observed for the \eqn{\ell^{\text{th}}} dataset and \eqn{p} is the number of covariables. When \code{L = 1}, \code{Z} can be a matrix.}

\item{beta}{the vector of the (true) coefficients values, with a length of \eqn{p} (the number of covariates).}

\item{a}{scale parameter, which should be non-negative. See \sQuote{Details} for the form of the cumulative hazard that can be used.}

\item{b}{shape/location parameter, which should be non-negative. It is not used when \code{gen_data_from = }\dQuote{exp}. See \sQuote{Details} for the form of the cumulative hazard that can be used.}

\item{u1}{a known non-negative lower limit of the uniform distribution for generating random censoring time. Default is \code{u1 = 0}. If \code{cen_rate} is not equal to 0, then \code{u1} does not need to be defined.}

\item{u2}{an non-negative upper limit of the uniform random censoring time distribution. The upper limit can be unknown (\code{u2 = NULL}, the default), or predefined. When this argument is assumed to be unknown, \code{u2 = NULL}, it is calculated by the algorithm within \code{surv.simulate()}. However, if the argument \code{u2} is known, the censoring rate cannot be predefined (meaning there is no control over it) and is calculated based on the generated dataset. See \sQuote{Details} and \sQuote{References}.}

\item{cen_rate}{a value representing the proportion of observations in the simulated survival data that are censored. The range of this argument is from 0 to 1. When the upper limit is known, \code{cen_rate} can nor be predefined. If there is no censoring (\code{cen_rate = 0}), the lower (\code{u1}) and upper (\code{u2}) limits of the uniform distribution do not need to be specified.}

\item{gen_data_from}{a description of the distribution from which the time to event is generated. This is a character string and can be \code{exponential} (\dQuote{exp}), \code{Weibull} (\dQuote{weibul}), or \code{Gompertz} (\dQuote{gomp}). Can be abbreviated. By default, the \code{exponential} distribution is used.}

\item{only_u2}{logical flag for calculating only the upper limit of the uniform censoring time distribution. If \code{only_u2 = TRUE}, the dataset(s) are not generated. If \code{only_u2 = TRUE}, the arguments \code{Z} and \code{u2} do not need to be specified, and \code{cen_rate} should not be set to \eqn{0}. Default is \code{only_u2 = FALSE}.}

\item{n.rep}{a scalar specifying the number of iterations. This argument is exclusively used in the case of the \code{Gompertz} distribution. Default is 100.}

\item{Trace}{logical flag indicating whether the output of the desired \code{u2} and the censoring proportion for different datasets should be produced for each iteration. It works \code{gen_data_from = }\dQuote{gomp}.}
}

\details{
\code{surv.simulate} function generates \eqn{L} simulated right-censored survival datasets from exponential, Weibull, or Gompertz distributions, incorporating the covariates, \code{Z}, distributed according to a multivariate \code{normal} distribution, with censoring time generated from a uniform distribution \code{Uniform(u1, u2)}, where \code{u1} is known but \code{u2} can be either known or unknown.

\code{surv.simulate()} can also be used to calculate the unknown upper limit of the uniform distribution, \code{u2}, with a predefined censoring rate. To do this, set \code{u2 = NULL} and \code{only_u2 = TRUE}. In this case, the datasets are not generated; only \code{u2} is.

\code{surv.simulate()} uses a root-finding algorithm to select the censoring parameter that achieves predefined censoring rates in the simulated survival data.

When \code{gen_data_from = }\dQuote{exp}:
\itemize{
  \item the cumulative baseline hazard function is considered as \eqn{\Lambda_0=a t},
  \item the event time for the \eqn{\ell^{\text{th}}} dataset, \eqn{T_\ell}, is computed by \eqn{ - log(u) \ exp(- Z_\ell \boldsymbol{\beta}) / a}, where \eqn{u} follows a standard uniform distribution;
}

For \code{gen_data_from = }\dQuote{weibul}:
\itemize{
  \item the cumulative hazard function is as \eqn{\Lambda_0=a t ^ b},
  \item the event time is computed by \eqn{T_\ell= (- log(u) \ exp(- Z_\ell \boldsymbol{\beta}) / a)^{1/b}}, where \eqn{u} follows a standard uniform distribution;
}

For \code{gen_data_from = }\dQuote{gomp}:
\itemize{
  \item the cumulative hazard function is as \eqn{\Lambda_0=a (exp(b t) - 1) / b},
  \item the event time is computed by \eqn{T_\ell= \log(1- log(u) \ exp(- Z_\ell \boldsymbol{\beta}) b / a) / b}, where \eqn{u} follows a standard uniform distribution;
}


Finally the survival time is obtained by \eqn{\tilde{T}_\ell=\min\{T_\ell , C_\ell \}}.

The function will be updated for \code{gen_data_from = }\dQuote{gomp}.

}

\value{

\code{surv.simulate} returns a list containing the following components:

\item{D}{a list of \eqn{L} data frames, with dimension \eqn{n_\ell \times (p+2)}. The first and second columns, named \code{time} and \code{status}, contain the simulated survival time and the censoring indicator, respectively, where \eqn{0} means censored and \eqn{1} means uncensored;}

\item{censor_propor}{the vector of censoring proportions in the simulated datasets \code{D}, containing \eqn{L} elements;}

\item{u1}{the lower limit of the uniform distribution used to generate random censoring times with a predefined censoring rate. Sometimes this output is less than the value entered by the user, as it is adjusted to achieve the desired amount of censoring rate;}

\item{u2}{the upper limit of the uniform distribution used to generate random censoring times. If \code{u2 = NULL}, this output will be the estimated upper limit necessary to achieve the desired censoring rate across the \eqn{L} datasets.}

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

#-------------------------
# Simulating Survival data
#-------------------------
N    <- c(7, 10, 13) # the sample sizes of 3 datasets
beta <- 1:4
p    <- length(beta)
L    <- 3

# Define a function to generate multivariate normal samples
mvrnorm_new <- function(n, mu, Sigma) {
    pp <- length(mu)
    e <- matrix(rnorm(n * pp), nrow = n)
    return(crossprod(t(e), chol(Sigma)) + matrix(mu, n, pp, byrow = TRUE))
}
Z <- list()
for (z in seq_len(L)) {
    Z[[z]] <- mvrnorm_new(n = N[z], mu = rep(0, p),
                          Sigma = diag(rep(1, p),p))
    colnames(Z[[z]]) <- paste0("Z_",seq_len(ncol(Z[[z]])))
}

# One simulated dataset from exponential distribution with no censoring:
surv_data <- surv.simulate(Z = Z[[1]], beta = beta, a = exp(-.9),
                           cen_rate = 0, gen_data_from = "exp")
surv_data
surv_data$D[[1]][,1:2] # The simulated survival data

# Calculate only 'u2' with a predefined censoring rate of 0.4:
u2_new <- surv.simulate(Z = Z[1:2], beta = beta, a = exp(-.9),
                        b = exp(1.8), u1 = 0.1, only_u2 = TRUE,
                        cen_rate = 0.4, gen_data_from = "weibul")$u2
u2_new

# Two simulated datasets with a known 'u2':
# Using 'u2_new' to help control over censoring rate (was chosen 0.4)
surv.simulate(Z = Z[1:2], beta = beta, a = exp(-.9), b = exp(1.8),
              u1 = 0.05, u2 = u2_new, gen_data_from = "weibul")

# Three simulated datasets from 'weibul' with an unknown 'u2':
surv.simulate(Z = Z, beta = beta, a = exp(-1), b = exp(1),
               u1 = 0.01, cen_rate = 0.3, gen_data_from = "weibul")

# Two simulated datasets from 'gomp' with unknown 'u2' and censoring rate of 0.3:
surv.simulate(Z = Z[2:3], beta = beta, a = exp(1), b = exp(2), u1 = 0.1,
              cen_rate = 0.3, gen_data_from = "gomp", Trace = TRUE)

}
\keyword{datagen}
\keyword{survival}
