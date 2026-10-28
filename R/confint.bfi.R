#' Approximate posterior credible intervals for BFI models
#'
#' Computes approximate posterior credible intervals using the
#' curvature-based covariance matrix.
#' For Gaussian residual variances, intervals are constructed on
#' the log(sigma2) scale and back-transformed.
#'
#' @param object An object inheriting from class "bfi".
#' @param parm Parameter names or positions. Defaults to all parameters.
#' @param level Credible level. Default is 0.95.
#' @param ... Additional arguments, currently unused.
#'
#' @return A matrix of approximate posterior credible intervals.
#'
#' @importFrom stats confint
#' @export
confint.bfi <- function(object, parm, level = 0.95, ...) {

  if (!inherits(object, "bfi")) {
    stop("object must inherit from class 'bfi'.")
  }

  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must be a number between 0 and 1.")
  }

  est <- coef(object)
  se <- sqrt(diag(vcov(object)))

  z <- stats::qnorm((1 + level) / 2)

  ci <- cbind(
    est - z * se,
    est + z * se
  )

  rownames(ci) <- names(est)
  colnames(ci) <- sprintf(
    "%g %%",
    100 * c((1 - level) / 2, (1 + level) / 2)
  )

  if (object$family == "gaussian") {

    sigma_idx <- grepl("^sigma2", names(est))

    if (any(sigma_idx)) {

      if (any(!is.finite(est[sigma_idx])) ||
          any(est[sigma_idx] <= 0)) {
        stop("Residual variance estimates must be positive.")
      }

      # vcov() is on the original sigma2 scale.
      # Recover the log-scale SD using the delta method.
      se_log <- se[sigma_idx] / est[sigma_idx]

      ci[sigma_idx, 1] <-
        est[sigma_idx] * exp(-z * se_log)

      ci[sigma_idx, 2] <-
        est[sigma_idx] * exp(z * se_log)
    }
  }

  if (!missing(parm)) {

    if (is.numeric(parm)) {
      if (anyNA(parm) || any(!is.finite(parm)) ||
          any(parm != floor(parm)) ||
          any(parm < 1) || any(parm > length(est))) {
        stop("Invalid parameter positions in 'parm'.")
      }
      parm <- names(est)[parm]
    }

    if (!is.character(parm) || anyNA(parm) ||
        !all(parm %in% names(est))) {
      stop("Invalid parameter names in 'parm'.")
    }

    ci <- ci[parm, , drop = FALSE]
  }

  ci
}
