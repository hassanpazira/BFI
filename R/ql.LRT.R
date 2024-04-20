## This file created by Hassan Pazira

#' @export

ql.LRT <- function(initial_theta, y, X, alpha = 0.1, max_order = 2){
  # if (is.matrix(X))
  #   warning("Note that, if in 'X' there is a 'categorical' covariate with more than 2, \n",
  #           "  levels, then 'X' must be a 'data.frame' instead of a 'matrix'.")
  p <- ncol(X)
  ql <- 0
  if (length(initial_theta)!=(p+ql+1))
    stop(paste("'initial' should be a numerical vector of length 'p+1'", sQuote(p+ql+1),"."))
  if (NCOL(y)!=2) stop("'y' must have two columns including 'time' and 'status'.")
  if (colnames(y)[1]=="time" & colnames(y)[2]=="status") {
    time <- y$time
    status <- y$status
  } else {
    time <- y[,1]
    status <- y[,2]
  }
  negloglik_l_theta_Taylor_mle <- function(beta_dotomega, time, status, X, q_l){
    p <- ncol(X)
    lam_0l_omeg_l <- lambda.poly(time, p=p, q_l=q_l, beta_dotomega=beta_dotomega)
    if (any(lam_0l_omeg_l==0)) {
      which_0 <- which(lam_0l_omeg_l == 0)
      lam_0l_omeg_l[which_0] <- .Machine$double.eps
    }
    Lam_0l_omeg_l <- sapply(time, FUN = function(x)
      integrate(lambda.poly, lower=0, upper=x, p=p, q_l=q_l,
                beta_dotomega=beta_dotomega)$value)
    power <- sum(status * (log(lam_0l_omeg_l) + as.numeric(X %*% beta_dotomega[1:p])) -
                   Lam_0l_omeg_l * exp(as.numeric(X %*% beta_dotomega[1:p])))
    negloglik  <- - power
    return(as.numeric(negloglik))
  }
  optim.poly.MLE.loc <- function(initial_theta, time, status, X, q_l){
    p <- ncol(X)
    # if (length(initial_theta) != (p+q_l+1))
    #   stop("length of initial vector should be 'p+q_l+1'", sQuote(p+q_l+1),".")
    optim_origin <- optim(initial_theta, fn=negloglik_l_theta_Taylor_mle, gr=NULL,
                          time=time, status=status, X=X, q_l=q_l,
                          lower=rep(-Inf,length(initial_theta)),
                          upper=rep(Inf,length(initial_theta)), method="L-BFGS-B")
    optim_origin_par <- optim_origin$par
    theta_hat <- list(beta_hat=optim_origin_par[1:p],
                      omega_hat = optim_origin_par[(p+1):(p+q_l+1)],
                      value=optim_origin$value)
    return(theta_hat)
  }
  res_ql_0 <- try(optim.poly.MLE.loc(initial_theta=initial_theta,
                                     time=time, status=status, X=X, q_l=ql), TRUE)
  if(!is.null(attr(res_ql_0, "class")) | inherits(res_ql_0, "try-error")) {
    res_ql_0 <- try(optim.poly.MLE.loc(initial_theta=c(rep(0.5, p),-0.5),
                                       time=time, status=status, X=X, q_l=ql), TRUE)
    if(!is.null(attr(res_ql_0, "class")) | inherits(res_ql_0, "try-error")) {
      res_ql_0 <- try(optim.poly.MLE.loc(initial_theta=c(rep(1, p),-1),
                                         time=time, status=status, X=X, q_l=ql), TRUE)
      if(!is.null(attr(res_ql_0, "class")) | inherits(res_ql_0, "try-error")) {
        stop("The algorithm in q_l_LRT() does not converged! The point is: ", res_ql_0)
      }
    }
  }
  L0 <- - res_ql_0$value
  repeat {
    ql <- ql + 1
    #cat("ql is --> ",ql, "\n") # !!!
    res_ql_1 <- try(optim.poly.MLE.loc(initial_theta=c(res_ql_0$beta_hat,res_ql_0$omega_hat, 0),
                                       time=time, status=status, X=X, q_l=ql), TRUE)
    if(!is.null(attr(res_ql_1, "class")) | inherits(res_ql_1, "try-error")) {
      res_ql_1 <- try(optim.poly.MLE.loc(initial_theta=c(res_ql_0$beta_hat,res_ql_0$omega_hat, -0.5),
                                         time=time, status=status, X=X, q_l=ql), TRUE)
      if(!is.null(attr(res_ql_1, "class")) | inherits(res_ql_1, "try-error")) {
        res_ql_1 <- try(optim.poly.MLE.loc(initial_theta=c(res_ql_0$beta_hat,res_ql_0$omega_hat, 1),
                                    time=time, status=status, X=X, q_l=ql), TRUE)
        if(!is.null(attr(res_ql_1, "class")) | inherits(res_ql_1, "try-error")) {
          res_ql_1 <- try(optim.poly.MLE.loc(initial_theta=c(res_ql_0$beta_hat,res_ql_0$omega_hat, -1),
                                      time=time, status=status, X=X, q_l=ql), TRUE)
          if(!is.null(attr(res_ql_1, "class")) | inherits(res_ql_1, "try-error")) {
            res_ql_1 <- try(optim.poly.MLE.loc(initial_theta=rep(0.5, p+ql+1),
                                        time=time, status=status, X=X, q_l=ql), TRUE)
            if(!is.null(attr(res_ql_1, "class")) | inherits(res_ql_1, "try-error")) {
              res_ql_1 <- try(optim.poly.MLE.loc(initial_theta=c(rep(0.5, p+ql),-0.5),
                                          time=time, status=status, X=X, q_l=ql), TRUE)
              if(!is.null(attr(res_ql_1, "class")) | inherits(res_ql_1, "try-error")) {
                res_ql_1 <- try(optim.poly.MLE.loc(initial_theta=c(rep(0, p+ql), -1),
                                            time=time, status=status, X=X, q_l=ql), TRUE)
                if(!is.null(attr(res_ql_1, "class")) | inherits(res_ql_1, "try-error")) {
                  res_ql_1 <- try(optim.poly.MLE.loc(initial_theta=c(rep(-0.5, p+ql),-0.5),
                                              time=time, status=status, X=X, q_l=ql), TRUE)
                  if(!is.null(attr(res_ql_1, "class")) | inherits(res_ql_1, "try-error")) {
                    q_l <- ql - 1
                    warning(res_ql_1,"q_l <- ql in ql.LRT()")
                    break
                  }
                }
              }
            }
          }
        }
      }
    }
    L1 <- - res_ql_1$value
    if (-2 * (L0 - L1) < qchisq(1-alpha, df=1, lower.tail = T) | ql==max_order+1) { # ql==max_order
      q_l <- ql - 1
      #cat("q_l is ----> ",q_l, "\n") # !!!
      break
    }
    L0 <- L1
    res_ql_0 <- res_ql_1
  }
  return(q_l)
}
