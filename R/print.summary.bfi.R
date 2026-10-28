#' @export
#' @noRd

print.summary.bfi <- function(
    x,
    digits = x$summary_digits,
    cur_mat = x$summary_cur_mat,
    ...) {

  if (!inherits(x, "summary.bfi")) {
    stop("Expected an object of class 'summary.bfi'.")
  }

  if (is.null(x$stratified)) {
    cat("\nSummary of the local model:\n\n")
  } else {
    cat("\nSummary of the BFI model:\n\n")
  }

  if (is.null(x$stratified)) {
    cat("   Formula: ")

    if (x$family != "survival") {
      cat(x$formula, "\n")
    } else if (is.character(x$formula)) {
      cat(x$formula, "\n")
    }
  }

  cat("    Family:", sQuote(x$family), "\n")

  if (x$family != "survival") {
    cat("      Link:", sQuote(x$link))
  } else {
    cat("  Baseline:", sQuote(x$basehaz))
  }

  cat("\n\nCoefficients:\n\n")
  print(round(x$coef_table, digits = digits))

  if (x$family == "gaussian" &&
      any(grepl("^sigma2", names(x$estimate)))) {
    cat(
      "\nFor the residual variances, the standard deviation is obtained by the",
      "\ndelta method and the credible interval is computed on the log scale",
      "\nand back-transformed, so the interval is not symmetric.\n"
    )
  }

  if (x$family == "gaussian" &&
      !is.null(x$dispersion)) {
    cat(
      "\nDispersion parameter (sigma2): ",
      format(x$dispersion, digits = digits), "\n"
    )
  }

  if (x$family == "binomial") {
    cat("\nDispersion parameter (sigma2): ", 1, "\n")
  }

  if (is.null(x$stratified)) {
    if (x$family == "survival") {
      cat(
        "\nlog Lik Posterior: ",
        format(x$logLikPost, digits = digits), "\n"
      )
      cat(
        "      Convergence: ",
        format(x$convergence, digits = digits), "\n"
      )
    } else {
      cat(
        "            log Lik Posterior: ",
        format(x$logLikPost, digits = digits), "\n"
      )
      cat(
        "                  Convergence: ",
        format(x$convergence, digits = digits), "\n"
      )
    }
  } else if (x$family %in% c("binomial", "gaussian") &&
             !is.null(x$Ave_Treat)) {
    cat("\nAverage Treatment Effect (ATE): ", "\n")
    cat(
      "\n         IPTW: ",
      format(x$Ave_Treat$IPTW, digits = digits), "\n"
    )
    cat(
      "        wIPTW: ",
      format(x$Ave_Treat$wIPTW, digits = digits), "\n"
    )
  }

  if (cur_mat) {
    cat("\nMinus the Curvature Matrix: \n\n")
    print(round(x$A_hat, digits = digits))
  }

  invisible(x)
}
