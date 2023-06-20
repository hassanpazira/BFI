## This file created by Hassan Pazira at 22-06-2023
model.maker <- function (formula, data) {
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula","data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (attr(mt, "intercept") == 0)
    stop("Models without intercept are not allowed in this version of the package")
  if (is.empty.model(mt)) {stop("Model matrix is empty")}
  y <- model.response(mf, "any")
  X <- model.matrix(mt, mf, contrasts=NULL)
  arg <- list(y=y, X=X)
  return(arg)
}
