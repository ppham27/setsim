Mcustomize <-
  function(y,fit,cov) {
    n <- length(y)
	X <- fit$X
	EST <- as.vector(fit$EST)
    lik <- fit$lik
    H.lik <- function(y,X,EST){
      as.numeric(eval(lik))
    }
    max_loglik <- H.lik(y,X,EST)
    input <- list(y=y,X=X,EST=EST,max_loglik=max_loglik)
    H <- function(epsilon,input,VZ,cstar) {
      H.lik(y=y,X=X,EST=as.numeric(EST+epsilon*VZ))-max_loglik+0.5*cstar
    }
    list(cov=cov,input=input,H=H,H.lik=H.lik)
  }
