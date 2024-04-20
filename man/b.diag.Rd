\name{b.diag}
\alias{b.diag}
\title{Create a Block Diagonal Matrix}
\description{
Construct a block diagonal matrix using multiple given block matrices.
}
\usage{
b.diag(..., fill = 0)
}

\arguments{
\item{...}{individual matrices or one list of matrices.}

\item{fill}{non-block-diagonal elements. Default is \eqn{0}.}
}

\details{
Avoid combining matrices and lists for the \code{...} argument.

\code{b.diag()} covers the arguments of type "character".

If a \emph{sparse} matrix needed, run the following:

\code{library(Matrix); Matrix(b_diag, sparse = TRUE)}

where \code{b_diag} is the matrix returned by \code{b.diag()}.
}

\value{
\code{b.diag()} returns a block diagonal matrix obtained by combining the arguments.
}

\author{Hassan Pazira\cr
Maintainer: Hassan Pazira \email{hassan.pazira@radboudumc.nl}
}


\examples{

b.diag(1, matrix(1:3, 3,4), diag(3:2))

b.diag(matrix(1:6, 2), as.character(2))

lists <- list(1, 2:3, diag(4:6), 7, cbind(8,9:12), 13:15)
b.diag(lists)
identical(b.diag(lists), b.diag(lapply(lists, as.matrix)))

b.diag(replicate(3, matrix(round(rnorm(9)), 3, 3), simplify=FALSE))

}
\keyword{array}
