## This file created by Hassan Pazira at 19-09-2026


#' @rdname bfi-methods
#' @importFrom stats vcov
#' @export
vcov.bfi <- function(object, ...) {
  if (class(object)[1] != "bfi") {
    stop("vcov method is not available for objects with other classes.")
  }

  cf <- coef(object)
  V <- solve(as.matrix(object$A_hat))

  if (object$family == "gaussian") {
    sigma_idx <- grepl("^sigma2", names(cf))

    if (any(sigma_idx)) {
      d <- rep(1, length(cf))
      d[sigma_idx] <- cf[sigma_idx]

      # Transform covariance matrix from the log(sigma2) scale
      # to the original sigma2 scale using the multivariate delta method.
      V <- V * tcrossprod(d)
    }
  }

  nm <- names(cf)
  dimnames(V) <- list(nm, nm)

  V
}
