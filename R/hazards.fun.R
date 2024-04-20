## This file created by Hassan Pazira

#' @export

hazards.fun <- function(time, z = NULL, p, theta_hat,
                        basehaz = c("weibul", "exp", "gomp", "poly", "pwexp"),
                        q_max, timax){
  # the inputs for theta_hat are (beta, omega_a, omega_b, ...) where omega's are \in R (+ and -)
  ti <- time
  basehaz <- match.arg(basehaz)
  if (!is.null(z)) { # New sample
    if (is.data.frame(z)) z <- as.matrix(z)
    if (is.vector(z)) z <- matrix(z,1,length(z))
    if (is.matrix(z) & (ncol(z)==1 & nrow(z)!=1)) z <- t(z)
    if (!missing(p)) {if (p != ncol(z)) stop("p != ncol(z)")}
    if (missing(p)) p <- ncol(z)
    if (nrow(z) != 1) stop("New sample is for one individual (n = 1).")
    beta <- theta_hat[1:p]
    reltive_risk <- exp(as.numeric(z %*% beta))
  } else {
    if (missing(p)) stop("When 'z=NULL', 'p' cannot be missing.")
    reltive_risk <- 1
  }
  if (basehaz == "exp") {
    a <- exp(theta_hat[p+1])
    bhazard <- rep(a, length(ti))
    cbhazard <- a * ti
    hazard <- bhazard * reltive_risk
    chazard <- cbhazard * reltive_risk
    survival <- exp(-chazard)
    bsurvival <- exp(-cbhazard)
  }
  if (basehaz == "gomp") {
    a <- exp(theta_hat[p+1])
    b <- exp(theta_hat[p+2])
    bhazard <- a * exp(b * ti)
    cbhazard <- a * (exp(b * ti) - 1) / b
    hazard <- bhazard * reltive_risk
    chazard <- cbhazard * reltive_risk
    survival <- exp(-chazard)
    bsurvival <- exp(-cbhazard)
  }
  if (basehaz == "weibul") {
    a <- exp(theta_hat[p+1])
    b <- exp(theta_hat[p+2])
    bhazard <- a * b * ti^(b-1)
    cbhazard <- a * ti^b
    hazard <- bhazard * reltive_risk
    chazard <- cbhazard * reltive_risk
    survival <- exp(-chazard)
    bsurvival <- exp(-cbhazard)
  }
  if (basehaz == "poly") {
    if (missing(q_max) | is.null(q_max))
      stop("'q_max' should not be missing or NULL for the 'poly' basehaz.")
    if (length(theta_hat) != (p+q_max+1))
      stop("'length(theta_hat)' is not equal to 'p+q_max+1'.")
    bhazard <- lambda.poly(ti, p, q_max, theta_hat)
    cbhazard <- NULL
    for (t in 1:length(ti)) {
      cbhazard[t] <- integrate(lambda.poly, 0, ti[t], p, q_max, theta_hat)$value
    }
    hazard <- bhazard * reltive_risk
    chazard <- cbhazard * reltive_risk
    survival <- exp(-chazard)
    bsurvival <- exp(-cbhazard)
  }
  if (basehaz == "pwexp") {
    if (missing(timax) | is.null(timax))
      stop("'timax' should not be missing or NULL for the 'pwexp' basehaz.")
    omega <- theta_hat[-c(1:p)]
    tps <- seq(0, timax + 0.001, length.out = length(omega)+1)
    b <- i.basis(tps, ti, ibasis = FALSE) #compute basis splines of degree zero, i.e. piecewise constant basis
    B <- i.basis(tps, ti, ibasis = TRUE) #compute integral of basis splines of degree zero ()
    bhazard <- as.numeric(b %*% exp(omega))
    cbhazard <- as.numeric(B %*% exp(omega))
    hazard <- bhazard * reltive_risk
    chazard <- cbhazard * reltive_risk
    survival <- exp(-chazard)
    bsurvival <- exp(-cbhazard)
  }
  names(bhazard) <- names(hazard) <- NULL
  return(list(bhazard=bhazard, cbhazard=cbhazard, bsurvival=bsurvival,
              hazard=hazard, chazard=chazard, survival=survival))
}
