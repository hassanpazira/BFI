## This file created by Hassan Pazira

#' @export

n.par <- function(X, family = c("gaussian", "binomial", "survival")){
  if (is.character(family))
    family <- match.arg(family)
  else stop("'family' should be a character string.")
  if (!is.list(X) | (is.list(X) & is.data.frame(X))) {
    if (family == "survival") {
      if (is.null(colnames(X)))
        stop("There is no colnames for center ",sQuote(l),".")
      wts <- which(colnames(X) %in% c("time", "status"))
      if (length(wts) > 0) {
        ZZ <- X[,-wts, drop=F]
      } else {
        ZZ <- X
      }
      # ZZ <- ZZ[,sort(colnames(ZZ))]
    } else {
      if (all(as.numeric(X[, 1]) == rep(1, nrow(X)))) {
        X <- X[, -1, drop = FALSE]
      }
      if (is.null(colnames(X))) {
        # colnames(X) <- paste0("X", seq_len(ncol(X)))
        stop("Colnames of X cannot be NULL.")
      } else {
        if (any(c("Intercept", "(Intercept)") %in% colnames(X))) {
          stop("'Intercept' should be the first column of 'X'.")
        }
      }
      ZZ <- X#[,sort(colnames(X))]
    }
    pre_p <- ncol(ZZ)
    pre_colnames <- colnames(ZZ)
    design_matrix <- paste(colnames(ZZ), collapse=" + ")
    formula <- try(as.formula(paste(" ", design_matrix, sep=" ~ ")),TRUE)
    if(inherits(formula, "try-error"))
      stop("It seems the data has already some dummy variables.",
           "Give the original data without dummies.")
    X_vars <- model.maker(formula, as.data.frame(ZZ), family)
    ZZ <- X_vars$X[,-1, drop=F]
    if (is.null(dim(ZZ))) {
      p <- 1
    } else {
      p <- ncol(ZZ)
    }
    return(list(n.reg.par = p, n.covar = pre_p, n.sample = nrow(ZZ),
                orig.names = pre_colnames, covar.names = colnames(ZZ)))
  } else {
    pp <- nn <- pre_p <- NULL
    pre_colnames <- pos_colnames <- list()
    for (l in 1:length(X)) {
      if (family == "survival") {
        if (is.null(colnames(X[[l]])))
          stop("There is no colnames for center ",sQuote(l),".")
        wts <- which(colnames(X[[l]]) %in% c("time", "status"))
        if (length(wts) > 0) {
          ZZ <- X[[l]][,-wts, drop=F]
        } else {
          ZZ <- X[[l]]
        }
        # ZZ <- ZZ[,sort(colnames(ZZ))]
      } else {
        if (all(as.numeric(X[[l]][, 1]) == rep(1, nrow(X[[l]])))) {
          X[[l]] <- X[[l]][, -1, drop = FALSE]
        }
        if (is.null(colnames(X[[l]]))) {
          # colnames(X[[l]]) <- paste0("X", seq_len(ncol(X[[l]])))
          stop("Colnames of X cannot be NULL for center ",sQuote(l),".")
        } else {
          if (any(c("Intercept", "(Intercept)") %in% colnames(X[[l]]))) {
            stop("'Intercept' should be the first column of 'X'.")
          }
        }
        ZZ <- X[[l]]#[,sort(colnames(X[[l]]))]
      }
      pre_p[l] <- ncol(ZZ)
      pre_colnames[[l]] <- colnames(ZZ)
      design_matrix <- paste(colnames(ZZ), collapse=" + ")
      formula <- try(as.formula(paste(" ", design_matrix, sep=" ~ ")),TRUE)
      if(inherits(formula, "try-error")) {
        stop("It seems the data has already some dummy variables. ",
             "Give the original data without dummies.")
      }
      X_vars <- model.maker(formula, as.data.frame(ZZ), family)
      ZZ <- X_vars$X[,-1, drop=F]
      if (is.null(dim(ZZ))) { #!
        pp[l] <- 1
        nn[l] <- length(ZZ)
        pos_colnames[[l]] <- names(ZZ)
      } else {
        pp[l] <- ncol(ZZ)
        nn[l] <- nrow(ZZ)
        pos_colnames[[l]] <- colnames(ZZ)
      }
    }
    return(list(n.reg.par=pp, n.covar=pre_p, n.sample=nn, orig.names=pre_colnames,
                covar.names=pos_colnames))
  }
}
