## This file created by Hassan Pazira at 19-09-2026


#' Methods for objects of class \code{"bfi"}
#'
#' Print, coefficient and variance methods for objects returned by
#' \code{\link{MAP.estimation}} and \code{\link{bfi}}.
#'
#' @param x,object An object of class \code{"bfi"}.
#' @param digits Number of significant digits used for printing.
#' @param ... Further arguments (currently ignored).
#'
#' @details \code{vcov()} returns the approximate posterior covariance
#'   matrix corresponding to the parameterization returned by
#'   \code{coef()}. For Gaussian models, the inverse curvature matrix
#'   is defined on the \eqn{\log(\sigma^2)} scale for residual-variance
#'   parameters; these rows and columns are transformed to the original
#'   \eqn{\sigma^2} scale using the multivariate delta method.
#'
#' @return \code{print()} returns \code{x} invisibly. \code{coef()} returns a
#'   named numeric vector of the MAP or BFI estimates. \code{vcov()} returns
#'   the approximate posterior covariance matrix.
#'
#' @seealso \code{\link{summary.bfi}}
#'
#' @examples
#' X <- data.frame(x1 = rnorm(50))
#' y <- rnorm(50)
#' Lambda <- inv.prior.cov(X, lambda = 0.01, family = "gaussian")
#' fit <- MAP.estimation(y, X, family = "gaussian", Lambda = Lambda)
#' fit
#' coef(fit)
#' vcov(fit)
#'
#' @name bfi-methods
#' @rdname bfi-methods
#' @export
print.bfi <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (class(x)[1] != "bfi") {
    stop("print method is not available for objects with other classes.")
  }
  type <- if (is.null(x$stratified)) "Local MAP estimates" else "BFI estimates"
  cat("\n", type, " (family: ", x$family, ")\n\n", sep = "")
  print.default(format(coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n")
  invisible(x)
}
