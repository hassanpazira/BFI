\name{BFI-internal}
\title{Internal BFI Functions}
\alias{A.l.maker}
\alias{negloglik.theta}
\alias{model.maker}
\alias{optim.survival}
\alias{ql.LRT}
\alias{lambda.poly}
\alias{i.basis}
\description{Internal BFI functions.}
\value{

\code{A.l.maker()} returns a curvature matrix for each local.

\code{negloglik.theta()} returns the negative log-likelihood of regression coefficients and error variance, which can be utilized for optimization purposes.

\code{model.maker()} returns the outcome and design matrix, incorporating dummy variables when categorical covariates are present, e.g., by expanding factors to a set of dummy variables and expanding interactions similarly.

\code{optim.survival()} is a general-purpose optimization algorithm based on the ``L-BFGS-B'' method, and returns a vector of estimates for the coefficients and baseline hazard parameters.

\code{ql.LRT()} returns an optimal order/degree of the exponentiated polynomial model, obtained by the likelihood ratio test, for the 'cox' family.

\code{lambda.poly()} returns the baseline hazard function for the exponentiated polynomial form.

\code{i.basis()} creates the integral of the output of basis when \code{ibasis = TRUE}. This is also an indicator function between the \eqn{k^{th}} and \eqn{(k+1)^{th}} time points when \code{ibasis = FALSE}.

}
\details{These functions are not intended for use by users.}
\author{
Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}
\keyword{internal}
\keyword{optimize}
