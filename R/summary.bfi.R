summary.bfi <- function(object, curmat = FALSE, digits = max(3, getOption("digits") - 3)) {
  if(class(object)[1] != "bfi")
    stop("summary method is not yet available for objects with other classes.")
  #cat("\n\n===============================================================")
  cat("\nSummary of the model:")
  cat("\n\n    Formula: ")
  cat(object$formula,"\n")
  cat("     Family:", sQuote(object$family),"\n")
  if (object$family == c("binomial")) {
    linkf <- noquote("Logit")
  }
  if (object$family == c("gaussian")) {
    linkf <- noquote("identity")
  }
  cat("       Link:", sQuote(linkf))
  cat("\n\nCoefficients:\n\n")
  coef_only <- as.matrix(object$theta_hat[-length(object$theta_hat)])
  se <- object$sd[-length(object$sd)]/sqrt(object$n)
  coef_se <- cbind(coef_only, se)
  margin <- qnorm(0.975) * object$sd[-length(object$sd)]
  ci <- cbind(coef_only-margin, coef_only+margin)
  colnames(ci) <- c("2.5 %" , " 97.5 %")
  coef_sd_ci <- cbind(coef_se, ci)
  colnames(coef_sd_ci) <- c("Estimate", "Std. Error", "2.5 %" , "97.5 %")
  printCoefmat(coef_sd_ci)
  if(object$family != c("binomial"))
    cat("\nDispersion parameter (sigma2): ", format(object$theta_hat[length(object$theta_hat)], digits = digits), "\n")
  else
    #cat("\nDispersion parameter (sigma2) for",object$family, "family taken to be 1 \n")
    cat("\nDispersion parameter (sigma2): ", 1, "\n")
  cat("                      log Lik: ", format(-object$value, digits = digits), "\n")
  cat("                  Convergence: ", format(object$convergence, digits = digits), "\n\n")
  if (curmat) {
    #cat("---\n\n")
    cat("Curvature Matrix: \n\n")
    print(round(object$A_hat, digits = digits))
  }
  object$logLik <- -object$value
  object$link <- linkf
  object$se <- se
  if(object$family != c("binomial"))
    object$dispersion <- object$theta_hat[length(object$theta_hat)]
  else
    object$dispersion <- 1
  # if(object$family != c("binomial"))
  #   object$ResSE <- sqrt(object$theta_hat[[length(object$theta_hat)]])
  # else
  #   object$dispersion <- 1
  object$CI <- ci
  class(object) <- "summary.bfi"
  invisible(object)
}
