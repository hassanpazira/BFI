\name{n.par}
\alias{n.par}
\title{The Number of Predictors, Coefficients, and Observations}
\description{
\code{n.par} returns the number of regression parameters, covariates and observations present in X based on the selected family.
}
\usage{
n.par(X, family = c("gaussian", "binomial", "survival"))
}

\arguments{
\item{X}{design matrix of dimension \eqn{n \times p}, where \eqn{n} is the number of samples observed, and \eqn{p} is the number of predictors/covariables. It could be a matrix or a list of matrices.}

\item{family}{a description of the error distribution used to specify the model. This should be a character string, either ``\code{gaussian}'', ``\code{binomial}'', or ``\code{survival}''. Can be abbreviated. By default the \code{gaussian} family is used.}

}

\details{
\code{orig.names} and \code{covar.names} are the same if the all covariates in \code{X} are continuous. However, if there are at least one categorical variable in \code{X} with more than two categories, they are different.

}

\value{
\code{n.par} returns a list containing the following components:

\item{n.reg.par}{the number of regression parameters;}

\item{n.covar}{the number of covariates;}

\item{n.sample}{the number of samples/observations;}

\item{orig.names}{the original variable names excluding dummy variable names;}

\item{covar.names}{the variables names, including any dummy variable names (if applicable).}

}

\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}

\examples{

#--------------------
# family = "gaussian"
#--------------------

X0 <- data.frame(x1 = rnorm(50),                     # standard normal variable
                 x2 = sample(0:2, 50, replace=TRUE), # categorical variable
                 x3 = sample(0:1, 50, replace=TRUE)) # dichotomous variable
n.par(X0) # without dummy variables
X0$x2 <- as.factor(X0$x2)
X0$x3 <- as.factor(X0$x3)
n.par(X0)  # with dummy variables

X1 <- data.frame(Intercept = rep(1,30),
                 x1 = rnorm(30),                     # continuous variable
                 x2 = sample(0:2, 30, replace=TRUE)) # categorical variable
n.par(X1) # without dummy variables
X1$x2  <- as.factor(X1$x2)
n.par(X1) # without dummy variables

# a list of two data sets:
X01 <- list(X0, X1)
n.par(X01)

}
\keyword{regression}
