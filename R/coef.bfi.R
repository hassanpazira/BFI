## This file created by Hassan Pazira at 19-09-2026


#' @rdname bfi-methods
#' @importFrom stats coef
#' @export
coef.bfi <- function(object, ...) {
  if (class(object)[1] != "bfi") {
    stop("coef method is not available for objects with other classes.")
  }
  th <- object$theta_hat
  if (is.matrix(th)) {
    nm <- if (nrow(th) == 1L) colnames(th) else rownames(th)
  } else {
    nm <- names(th)
  }
  if (is.null(nm)) nm <- colnames(object$A_hat)
  th <- as.numeric(th)
  names(th) <- nm
  th
}