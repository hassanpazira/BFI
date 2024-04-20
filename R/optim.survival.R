## This file created by Hassan Pazira


#' @export

optim.survival <- function(initial_beta_loga_logb, q_l, tps, y, X, Lambda, family,
                           basehaz, wli, control){
  if (basehaz=="pwexp") {
    if (colnames(y)[1]=="time") {
      # compute basis splines of degree zero, i.e. piecewise constant basis
      bas <- i.basis(tps,y$time, ibasis = FALSE)
      # compute integral of basis splines of degree zero
      ibas <- i.basis(tps,y$time, ibasis = TRUE)
    } else {
      bas <- i.basis(tps,y[,1], ibasis = FALSE)
      ibas <- i.basis(tps,y[,1], ibasis = TRUE)
    }
  }
  if (basehaz %in% c("exp", "gomp", "weibul", "unspecified")) {
    beta_loga_logb_hat <- try(optim(initial_beta_loga_logb, fn=negloglik.theta,
                                    gr=NULL, y=y, X=X, Lambda=Lambda, family=family, q_l=q_l,
                                    tps=tps, bas=bas, ibas=ibas, basehaz=basehaz, wli=wli,
                                    method="BFGS", control=control), TRUE)
  } else {
    beta_loga_logb_hat <- try(optim(initial_beta_loga_logb, fn=negloglik.theta,
                                    gr=NULL, y=y, X=X, Lambda=Lambda, family=family, q_l=q_l,
                                    tps=tps, bas=bas, ibas=ibas, basehaz=basehaz, wli=wli,
                                    lower=rep(-Inf, length(initial_beta_loga_logb)),
                                    upper=rep(Inf, length(initial_beta_loga_logb)),
                                    method="L-BFGS-B", control=control), TRUE)
  }
  return(beta_loga_logb_hat)
}
