summary.bfi <- function(object, curmat = FALSE, digits = max(3, getOption("digits") - 3), ...) {
  if(class(object)[1] != "bfi")
    stop("summary method is not yet available for objects with other classes.")
  #cat("\n\n===============================================================")
  if (object$intercept == TRUE) {
    cat("\nSummary of the model (with 'intercept'):")
  } else {
    cat("\nSummary of the model (without 'intercept'):")
  }
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
  if (!object$intercept) {
    object$Estimate <- object$theta_hat[-c(1, length(object$theta_hat))]
    object$se <- object$sd[-c(1, length(object$sd))]#/sqrt(object$n)  #!!!
    object$CurMat <- object$A_hat[-1,-1]
  } else {
    object$Estimate <- object$theta_hat[-length(object$theta_hat)]
    object$se <- object$sd[-length(object$sd)]#/sqrt(object$n)  #!!!
    object$CurMat <- object$A_hat
  }
  coef_se <- cbind(object$Estimate, object$se)
  margin <- qnorm(0.975) * object$se #!
  ci <- cbind(object$Estimate-margin, object$Estimate+margin)
  colnames(ci) <- c("2.5 %" , " 97.5 %")
  coef_sd_ci <- cbind(coef_se, ci)
  colnames(coef_sd_ci) <- c("Estimate", "Std.Error", "CI 2.5%" , "CI 97.5%")
  #printCoefmat(coef_sd_ci, digits=digits)
  print(round(coef_sd_ci, digits=digits))
  if(object$family != c("binomial"))
    cat("\nDispersion parameter (sigma2): ", format(object$theta_hat[length(object$theta_hat)],
                                                    digits = digits), "\n")
  else
    #cat("\nDispersion parameter (sigma2) for",object$family, "family taken to be 1 \n")
    cat("\nDispersion parameter (sigma2): ", 1, "\n")
  cat("                      log Lik: ", format(-object$value, digits = digits), "\n")
  cat("                  Convergence: ", format(object$convergence, digits = digits), "\n\n")
  if (curmat) {
    #cat("---\n\n")
    cat("Curvature Matrix: \n\n")
    print(round(object$CurMat, digits = digits))
  }
  object$logLik <- -object$value
  object$link <- linkf
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
