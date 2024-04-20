

#' @export

optim.bfi <- function(initial_beta_loga_logb, q_l, tps, y, X, Lambda, family, basehaz){
  p <- dim(X)[2] # n. parameters
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
  if (basehaz %in% c("exp", "gomp", "weibul")) {
    beta_loga_logb_hat <- try(optim(initial_beta_loga_logb, fn=negloglik.theta,
                                    gr=NULL, y=y, X=X, Lambda=Lambda, family=family,
                                    q_l=q_l, tps=tps, bas=bas, ibas=ibas, basehaz=basehaz,
                                    method="BFGS"), TRUE)
  } else {
    beta_loga_logb_hat <- try(optim(initial_beta_loga_logb, fn=negloglik.theta,
                                    gr=NULL, y=y, X=X, Lambda=Lambda, family=family,
                                    q_l=q_l, tps=tps, bas=bas, ibas=ibas, basehaz=basehaz,
                                    lower=rep(-Inf,length(initial_beta_loga_logb)),
                                    upper=rep(Inf,length(initial_beta_loga_logb)),
                                    method="L-BFGS-B"), TRUE)
  }
  if(!is.null(attr(beta_loga_logb_hat, "class")) | inherits(beta_loga_logb_hat, "try-error")) {
    beta_loga_logb_hat <- try(optim(initial_beta_loga_logb-0.5, fn=negloglik.theta,
                                    gr=NULL, y=y, X=X, Lambda=Lambda, family=family,
                                    q_l=q_l, tps=tps, bas=bas, ibas=ibas, basehaz=basehaz,
                                    lower=rep(-Inf,length(initial_beta_loga_logb)),
                                    upper=rep(Inf,length(initial_beta_loga_logb)),
                                    method="L-BFGS-B"), TRUE)
    if(!is.null(attr(beta_loga_logb_hat, "class")) | inherits(beta_loga_logb_hat, "try-error")) {
      cat("try-error in beta_loga_logb_hat, now is -.5 ","\n", "\n") # !!! be deleted
      if (basehaz %in% c("exp", "gomp", "weibul")) { # !!! be deleted
        initial_beta_loga_logb <- c(rep(-1, p), rep(1, length(initial_beta_loga_logb)-p))
      } else {
        initial_beta_loga_logb <- rep(1, length(initial_beta_loga_logb))
      }
      beta_loga_logb_hat <- try(optim(initial_beta_loga_logb, fn=negloglik.theta,
                                      gr=NULL, y=y, X=X, Lambda=Lambda, family=family,
                                      q_l=q_l, tps=tps, bas=bas, ibas=ibas, basehaz=basehaz,
                                      lower=rep(-Inf,length(initial_beta_loga_logb)),
                                      upper=rep(Inf,length(initial_beta_loga_logb)),
                                      method="L-BFGS-B"), TRUE)
      if(!is.null(attr(beta_loga_logb_hat, "class")) | inherits(beta_loga_logb_hat, "try-error")) { # !!! be deleted
        cat("try-error in beta_loga_logb_hat, now is -5 & 5 ","\n", "\n") # !!! be deleted
        if (basehaz %in% c("exp", "gomp", "weibul")) { # !!! be deleted
          initial_beta_loga_logb <- c(rep(1, p), rep(-1, length(initial_beta_loga_logb)-p))
        } else {
          initial_beta_loga_logb <- rep(-1, length(initial_beta_loga_logb)) # or 0?
        }
        beta_loga_logb_hat <- try(optim(initial_beta_loga_logb, fn=negloglik.theta,
                                        gr=NULL, y=y, X=X, Lambda=Lambda, family=family,
                                        q_l=q_l, tps=tps, bas=bas, ibas=ibas, basehaz=basehaz,
                                        lower=rep(-Inf,length(initial_beta_loga_logb)),
                                        upper=rep(Inf,length(initial_beta_loga_logb)),
                                        method="L-BFGS-B"), TRUE)
        if(!is.null(attr(beta_loga_logb_hat, "class")) | inherits(beta_loga_logb_hat, "try-error")) {
          stop("Initial values should be changed!")
        }
      }
    }
  }
  return(beta_loga_logb_hat)
}
