## This file created by Hassan Pazira at 16-12-2022

#' @export
#' @noRd

summary.bfi <- function(object, cur_mat = FALSE,
                        digits = max(3, getOption("digits") - 3), ...) {
  if (!inherits(object, "bfi")) {
    stop("summary method is not available for objects with other classes.")
  }
  if (object$family == c("binomial")) {
    linkf <- noquote("Logit")
    object$estimate <- as.numeric(object$theta_hat)
    object$sd <- object$sd
  }
  if (object$family == c("gaussian")) {
    linkf <- noquote("identity")
    nm_all <- names(drop(object$theta_hat))
    # If is.null(object$stratified)==T, it means the 'object' is from MAP.estimation().
    if (is.null(object$stratified) | ((!is.null(object$stratified)) & (!c(2) %in% object$strat_par))) {
      object$estimate <- as.numeric(object$theta_hat[-length(object$theta_hat)])
      names(object$estimate) <- nm_all[-length(nm_all)]
      object$sd <- object$sd[-length(object$sd)]
    } else {
      object$estimate <- as.numeric(object$theta_hat)
      names(object$estimate) <- nm_all
    }
  }
  if (object$family == c("survival")) {
    linkf <- NULL
    len_omegas <- length(grep("omega", names(object$theta_hat)))
    object$estimate <- as.numeric(object$theta_hat[1:(length(object$theta_hat)-len_omegas)])
    object$sd <- object$sd[1:(length(object$sd)-len_omegas)]
  }

  # Standard deviations and 95% credible intervals
  sd_print <- object$sd
  margin <- qnorm(0.975) * object$sd
  ci <- cbind(object$estimate - margin, object$estimate + margin)

  # For Gaussian models with center-specific sigma2:
  # sigma2 estimates are reported on the original scale,
  # whereas their SDs are obtained on the log(sigma2) scale.
  # For printing, these SDs are transformed to the sigma2 scale using the delta metho,
  # whereas the 95 percent credible intervals are computed on the \eqn{\log(\sigma^2)} scale
  # and then back-transformeی  to the original \eqn{\sigma^2} scale.
  if (object$family == "gaussian") {

    sigma_idx <- grepl("^sigma2", names(object$estimate))

    if (any(sigma_idx)) {
      # 95% CI on the log(sigma2) scale, back-transformed to the sigma2 scale
      ci[sigma_idx, 1] <- object$estimate[sigma_idx] *
        exp(- qnorm(0.975) * object$sd[sigma_idx])

      ci[sigma_idx, 2] <- object$estimate[sigma_idx] *
        exp(+ qnorm(0.975) * object$sd[sigma_idx])

      # Delta-method SD on the sigma2 scale for printing
      sd_print[sigma_idx] <-
        object$estimate[sigma_idx] * object$sd[sigma_idx]
    }
  }

  coef_sd <- cbind(object$estimate, sd_print)
  colnames(ci) <- c("2.5 %", " 97.5 %")
  coef_sd_ci <- cbind(coef_sd, ci)
  colnames(coef_sd_ci) <- c("Estimate", "Std.Dev", "CI 2.5%", "CI 97.5%")

  if (is.null(object$stratified)) {
    object$logLikPost <- -object$value
    object$value <- NULL
  }
  object$link <- linkf
  if (object$family == "gaussian") {
    if (is.null(object$stratified) ||
        (!is.null(object$stratified) && !2 %in% object$strat_par)) {
      object$dispersion <- object$theta_hat[length(object$theta_hat)]
    } else {
      object$dispersion <- NULL
    }
  } else if (object$family == "binomial") {
    object$dispersion <- 1
  } else {
    object$dispersion <- NULL
  }

  object$sd <- sd_print
  object$CI <- ci
  object$coef_table <- coef_sd_ci
  object$summary_digits <- digits
  object$summary_cur_mat <- cur_mat

  class(object) <- "summary.bfi"
  object

}
