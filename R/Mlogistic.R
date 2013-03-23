
Mlogistic <-
  function(y,fit,cov) {
	n <- length(y)
	X <- model.matrix(fit)
	MLE <- as.vector(fit$coef)
	pi <- exp(MLE%*%t(X))/(1+exp(MLE%*%t(X)))
	max_loglik<-log(pi)%*%y + log(1-pi)%*%(1-y)
	input <- list(y=y,X=X,MLE=MLE,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar) {
      tmp <- exp(input$X %*% (input$MLE+epsilon*VZ))
      input$y %*% log(tmp/(1+tmp)) + (1-input$y) %*% log(1/(1+tmp)) -
        input$max_loglik+0.5*cstar
    }
	H.lik <- function(y,X,MLE) {
      tmp <- exp(X %*% MLE)
      y %*% log(tmp/(1+tmp)) + (1-y) %*% log(1/(1+tmp))
    }
    if (is.null(cov)) {
      cov <- vcov(fit)
      list(cov=cov,input=input,H=H,H.lik=H.lik)
    }
	else {
      list(input=input,H=H,H.lik=H.lik)
    }
  }
