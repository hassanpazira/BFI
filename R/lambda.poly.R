## This file created by Hassan Pazira

#' @export

lambda.poly <- function(s, p, q_l, beta_dotomega){
  if (length(beta_dotomega) != (p+q_l+1)) {
    stop("The length of 'theta' should be 'p+q_l+1' = ",sQuote(p+q_l+1),".")
  }
  log_lam_0l_omeg_l <- 0
  for (ql in 0:q_l) {
    log_lam_0l_omeg_l <- log_lam_0l_omeg_l + (beta_dotomega[p+ql+1]) * s^(ql)
  }
  lam_0l_omeg_l <- exp(log_lam_0l_omeg_l)
  return(lam_0l_omeg_l)
}
