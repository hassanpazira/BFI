## This file created by Hassan Pazira

#' @export

i.basis <- function(tps, t, ibasis = TRUE) {
  #function that creates the integral of the output of basis
  res <- matrix(0, nrow = length(t), ncol = (length(tps) - 1))
  for (k in 1:(length(tps) - 1)) {
    if (ibasis == TRUE) {
      res[, k] <- pmin(t - tps[k], tps[k + 1] - tps[k]) * as.numeric(t >= tps[k])
    } else {
      #indicator function between tps[k] and tps[k+1]
      res[, k] <- as.numeric(t >= tps[k]) * as.numeric(t < tps[k + 1])
    }
  }
  return(res)
}
