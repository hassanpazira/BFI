## This file created by Hassan Pazira at 16-12-2022
## Updated at 14-09-2023

summary.bfi <- function(object, cur_mat = FALSE, digits = max(3, getOption("digits") - 3), ...) {
  if(class(object)[1] != "bfi")
    stop("summary method is not available for objects with other classes.")
  if(is.null(object$stratified) == FALSE) {
    if(object$stratified == TRUE)
      stop("summary method is not available for stratification analysis.")
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
  coef_sd <- cbind(object$estimate, object$sd)
  margin <- qnorm(0.975) * object$sd
  ci <- cbind(object$estimate - margin, object$estimate + margin)
  colnames(ci) <- c("2.5 %" , " 97.5 %")
  coef_sd_ci <- cbind(coef_sd, ci)
  colnames(coef_sd_ci) <- c("Estimate", "Std.Dev", "CI 2.5%" , "CI 97.5%")
  cat("\nSummary of the model:\n\n")
  if(is.null(object$stratified) == TRUE) {
    cat("   Formula: ")
    cat(object$formula,"\n")
  }
  cat("    Family:", sQuote(object$family),"\n")
  cat("      Link:", sQuote(linkf))
  cat("\n\nCoefficients:\n\n")
  print(round(coef_sd_ci, digits=digits)) # printCoefmat(coef_sd_ci, digits=digits)
  if(object$family != c("binomial")) {
    cat("\nDispersion parameter (sigma2): ", format(object$theta_hat[length(object$theta_hat)],
                                                    digits = digits), "\n")
  } else {
    #cat("\nDispersion parameter (sigma2) for",object$family, "family taken to be 1 \n")
    cat("\nDispersion parameter (sigma2): ", 1, "\n")
  }
  if(is.null(object$stratified) == TRUE) {
    cat("                 log Lik Post: ", format(-object$value, digits = digits), "\n")
    cat("                  Convergence: ", format(object$convergence, digits = digits), "\n")
    object$logLikPost <- -object$value
    object$value <- NULL
  }
  if (cur_mat) {
    #cat("---\n\n")
    cat("\nCurvature Matrix: \n\n")
    print(round(object$A_hat, digits = digits))
  }
  object$link <- linkf
  if(object$family != c("binomial"))
    object$dispersion <- object$theta_hat[length(object$theta_hat)]
  else
    object$dispersion <- 1
  object$CI <- ci
  class(object) <- "summary.bfi"
  invisible(object)
}
