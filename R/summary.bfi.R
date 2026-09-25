## This file created by Hassan Pazira at 16-12-2022

#' @export

summary.bfi <- function(object, cur_mat = FALSE,
                        digits = max(3, getOption("digits") - 3), ...) {
  if (class(object)[1] != "bfi") {
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

  if (is.null(object$stratified)) # It means the object is from MAP.estimation()
    cat("\nSummary of the local model:\n\n")
  else
    cat("\nSummary of the BFI model:\n\n")
  if (is.null(object$stratified)) {
    cat("   Formula: ")
    if (object$family != "survival") cat(object$formula, "\n")
    if (object$family == "survival") {
      if (is.character(object$formula)) {
        cat(object$formula,"\n")
      }
    }
  }
  cat("    Family:", sQuote(object$family), "\n")
  if (object$family != "survival") cat("      Link:", sQuote(linkf))
  else cat("  Baseline:", sQuote(object$basehaz))
  cat("\n\nCoefficients:\n\n")
  print(round(coef_sd_ci, digits = digits))
  if (object$family == "gaussian" && any(grepl("^sigma2", names(object$estimate)))) {
    cat("\nFor the residual variances, the standard deviation is obtained by the",
        "\ndelta method and the credible interval is computed on the log scale",
        "\nand back-transformed, so the interval is not symmetric.\n")
  }
  #printCoefmat(coef_sd_ci, digits=digits)
  if (object$family == c("gaussian")) {
    if (is.null(object$stratified) | (!is.null(object$stratified) & (!c(2) %in% object$strat_par))) {
      cat("\nDispersion parameter (sigma2): ",
          format(object$theta_hat[length(object$theta_hat)],
                 digits = digits
          ), "\n")
    }
  }
  if (object$family == c("binomial")) {
    # cat("\nDispersion parameter (sigma2) for",object$family, "family taken to be 1 \n")
    cat("\nDispersion parameter (sigma2): ", 1, "\n")
  }
  if (is.null(object$stratified)) { # It means the object is from MAP.estimation()
    if (object$family == c("survival")) {
      cat("\nlog Lik Posterior: ", format(-object$value, digits = digits), "\n")
      cat("      Convergence: ", format(object$convergence, digits = digits), "\n")
    } else {
      cat("            log Lik Posterior: ", format(-object$value, digits = digits), "\n")
      cat("                  Convergence: ", format(object$convergence, digits = digits), "\n")
    }
    object$logLikPost <- -object$value
    object$value <- NULL
  } else {
    if (object$family %in% c("binomial","gaussian")) {
      if (!is.null(object$Ave_Treat)) {
        cat("\nAverage Treatment Effect (ATE): ", "\n")
        cat("\n         IPTW: ", format(object$Ave_Treat$IPTW, digits = digits), "\n")
        cat("        wIPTW: ", format(object$Ave_Treat$wIPTW, digits = digits), "\n")
      }
      # if (!is.null(object$S_var)) {
      #   cat("\nSample Variance: ", "\n")
      #   cat("\n    Treatment: ", format(object$S_var$treatment, digits = digits), "\n")
      #   cat("      Control: ", format(object$S_var$control, digits = digits), "\n")
      # }
    }
    # if (object$family %in% c("survival")) {
    #   if (object$basehaz != c("unspecified")) {
    #     if (!is.null(object$Ave_Treat)) {
    #       cat("\nAverage Treatment Effect (ATE): ", "\n")
    #       cat("\n         IPTW: ", format(object$Ave_Treat$IPTW, digits = digits), "\n")
    #       cat("        wIPTW: ", format(object$Ave_Treat$wIPTW, digits = digits), "\n")
    #     }
    #     # if (!is.null(object$S_var)) {
    #     #   cat("\nSample Variance: ", "\n")
    #     #   cat("\n    Treatment: ", format(object$S_var$treatment, digits = digits), "\n")
    #     #   cat("      Control: ", format(object$S_var$control, digits = digits), "\n")
    #     # }
    #   }
    # }
  }
  if (cur_mat) {
    # cat("---\n\n")
    cat("\nMinus the Curvature Matrix: \n\n")
    print(round(object$A_hat, digits = digits))
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
  object$sd <- sd_print #!
  object$CI <- ci
  class(object) <- "summary.bfi"
  invisible(object)
}
