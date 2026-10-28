# Convert additive R formulas to the original BFI y/X interface.
# This helper does not modify the statistical estimation procedure.

.bfi_formula_inputs <- function(formula, data, family,
                                response = TRUE, treatment = NULL) {

  if (!inherits(formula, "formula")) {
    stop("'formula' must be an R formula.")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame when using 'formula'.")
  }

  trm <- stats::terms(formula, data = data)

  if (response && attr(trm, "response") != 1L) {
    stop("A response is required in the formula.")
  }

  if (length(attr(trm, "offset"))) {
    stop("Offset terms are not supported.")
  }

  predictors <- attr(trm, "term.labels")

  if (!length(predictors)) {
    stop("At least one predictor is required.")
  }

  if (any(!predictors %in% names(data))) {
    stop(
      "Use additive formulas with original column names. ",
      "Interactions and transformations must be created ",
      "as separate columns first."
    )
  }

  # inv.prior.cov() needs only predictor columns.
  frame_terms <- if (response) {
    trm
  } else {
    stats::delete.response(trm)
  }

  mf <- stats::model.frame(
    formula = frame_terms,
    data = data,
    na.action = stats::na.fail,
    drop.unused.levels = FALSE
  )

  X <- mf[, predictors, drop = FALSE]

  has_intercept <- attr(trm, "intercept") == 1L

  # Preserve the original predictor columns, including factors.
  # The existing BFI functions apply factor contrasts internally,
  # using a reference category for categorical predictors.
  #
  # This also applies to no-intercept formulas, ensuring
  # compatibility with the original y/X interface.

  if (family %in% c("gaussian", "binomial") &&
      !has_intercept) {

    if (ncol(X) < 1L) {
      stop("The design matrix contains no predictors.")
    }

    # The legacy functions would silently remove a first
    # predictor consisting entirely of ones.
    if (is.numeric(X[[1L]]) &&
        all(X[[1L]] == 1)) {

      stop(
        "The first predictor is a column of ones. ",
        "The existing BFI code would remove it."
      )
    }
  }

  y <- if (response) {
    stats::model.response(mf, type = "any")
  } else {
    NULL
  }

  if (response && family == "survival") {

    if (inherits(y, "Surv")) {

      if (!identical(attr(y, "type"), "right")) {
        stop("Only right-censored survival data are supported.")
      }

      y <- data.frame(
        time = as.numeric(y[, 1L]),
        status = as.numeric(y[, 2L])
      )

    } else if (is.matrix(y) && ncol(y) == 2L) {

      y <- data.frame(
        time = y[, 1L],
        status = y[, 2L]
      )

    } else {

      stop(
        "Use survival::Surv(time, status) ",
        "or cbind(time, status) for survival responses."
      )
    }
  }

  list(
    X = X,
    y = y,
    intercept = has_intercept
  )
}
