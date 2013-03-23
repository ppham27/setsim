Mcustomize <-
  function(y,fit,cov) {
    n <- length(y)
	X <- fit$X
	MLE <- as.vector(fit$MLE)
    lik <- fit$lik
    H.lik <- function(y,X,MLE){
      as.numeric(eval(lik))
    }
    max_loglik <- H.lik(y,X,MLE)
    input <- list(y=y,X=X,MLE=MLE,max_loglik=max_loglik)
    H <- function(epsilon,input,VZ,cstar) {
      H.lik(y=y,X=X,MLE=as.numeric(MLE+epsilon*VZ))-max_loglik+0.5*cstar
    }
    list(cov=cov,input=input,H=H,H.lik=H.lik)
  }
