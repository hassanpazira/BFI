## This file created by Hassan Pazira

#' @export

## This file created by Hassan Pazira at 22-06-2023
model.maker <- function(formula, data, family) {
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (is.empty.model(mt)) stop("Model matrix is empty")
  X <- model.matrix(mt, mf, contrasts = NULL)
  if (family == "survival") {
    arg <- list(X = X)
  } else {
    if (attr(mt, "intercept") == 0)
      stop("Models without intercept are not allowed in this version!")
    y <- model.response(mf, "any")
    arg <- list(y = y, X = X)
  }
  return(arg)
}
