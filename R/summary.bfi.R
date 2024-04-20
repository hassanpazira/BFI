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
    object$estimate <- as.numeric(object$theta_hat[-length(object$theta_hat)])
    object$sd <- object$sd[-length(object$sd)]
  }
  if (object$family == c("survival")) {
    linkf <- NULL
    len_omegas <- length(grep("omega",rownames(object$theta_hat)))
    object$estimate <- as.numeric(object$theta_hat[1:(length(object$theta_hat)-len_omegas)])
    object$sd <- object$sd[1:(length(object$sd)-len_omegas)]
  }
  coef_sd <- cbind(object$estimate, object$sd)
  margin <- qnorm(0.975) * object$sd
  ci <- cbind(object$estimate - margin, object$estimate + margin)
  colnames(ci) <- c("2.5 %", " 97.5 %")
  coef_sd_ci <- cbind(coef_sd, ci)
  colnames(coef_sd_ci) <- c("Estimate", "Std.Dev", "CI 2.5%", "CI 97.5%")
  if (is.null(object$stratified)) cat("\nSummary of the local model:\n\n")
  else cat("\nSummary of the BFI model:\n\n")
  if (is.null(object$stratified)) {
    cat("   Formula: ")
    if (object$family != "survival") cat(object$formula, "\n")
    if (object$family == "survival") {
      cat(deparse(update.formula(object$formula, Survival(time, status) ~ .)),"\n")
    }
  }
  cat("    Family:", sQuote(object$family), "\n")
  if (object$family != "survival") cat("      Link:", sQuote(linkf))
  else cat("  Baseline:", sQuote(object$basehaz))
  cat("\n\nCoefficients:\n\n")
  print(round(coef_sd_ci, digits = digits))
  #printCoefmat(coef_sd_ci, digits=digits)
  if (object$family == c("gaussian")) {
    cat("\nDispersion parameter (sigma2): ",
        format(object$theta_hat[length(object$theta_hat)],
               digits = digits
        ), "\n")
  }
  if (object$family == c("binomial")) {
    # cat("\nDispersion parameter (sigma2) for",object$family, "family taken to be 1 \n")
    cat("\nDispersion parameter (sigma2): ", 1, "\n")
  }
  if (is.null(object$stratified)) {
    if (object$family == c("survival")) {
      cat("\nlog Lik Posterior: ", format(-object$value, digits = digits), "\n")
      cat("      Convergence: ", format(object$convergence, digits = digits), "\n")
    } else {
      cat("            log Lik Posterior: ", format(-object$value, digits = digits), "\n")
      cat("                  Convergence: ", format(object$convergence, digits = digits), "\n")
    }
    object$logLikPost <- -object$value
    object$value <- NULL
  }
  if (cur_mat) {
    # cat("---\n\n")
    cat("\nMinus the Curvature Matrix: \n\n")
    print(round(object$A_hat, digits = digits))
  }
  object$link <- linkf
  if (object$family != c("binomial")) {
    object$dispersion <- object$theta_hat[length(object$theta_hat)]
  } else {
    object$dispersion <- 1
  }
  object$CI <- ci
  class(object) <- "summary.bfi"
  invisible(object)
}
