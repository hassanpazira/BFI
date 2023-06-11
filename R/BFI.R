## This file created by Hassan Pazira at 16-12-2022
bfi <- function(theta_hats=NULL, A_hats, Gamma, L=NULL) {
  if (!is.null(theta_hats) & !is.list(theta_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.list(A_hats)) {
    stop("The input for 'theta_hats' should be a list.")
  }
  if (!is.null(theta_hats) & (length(theta_hats) != length(A_hats))) {
    stop("Length of inputs are not the same.")
  }
  if (is.null(L)) {
    L <- length(A_hats)
  } else {
    if (length(theta_hats) != L) {
      stop("Number of locations and length of inputs are not the same.")
    }
  }
  if (is.matrix(Gamma) | (is.list(Gamma) & length(Gamma)==1)) { # all locations have the same Gamma.
    A_fed <- Reduce("+", A_hats) + (1-L) * Gamma
  } else {
    if (!is.list(Gamma) | (is.list(Gamma) & length(Gamma)!=(L+1))) {
      stop("Gamma should contains of L+1 lists; 1:L for local datastes and last one for combined.")
    }
    A_fed <- Reduce("+", A_hats) + Gamma[[L+1]] - Reduce("+", Gamma[1:L])
  }
  sd_fed <- sqrt(diag(solve(as.matrix(A_fed))))
  if (!is.null(theta_hats)) {
    theta_hat_fed <- solve(as.matrix(A_fed)) %*% Reduce("+", Map("%*%", A_hats , theta_hats))
  } else theta_hat_fed <- theta_hats
  output <- list(theta_hat=theta_hat_fed, A_hat=A_fed, sd=sd_fed)
  return(output)
}
