Mpoisson <-
  function(y,fit,cov) {
	n <- length(y)
	X <- model.matrix(fit)
	EST <- as.vector(fit$coef)
	max_loglik <- y %*% X %*% EST -
      sum(exp(X %*% EST)) -
        sum(lgamma(y+1))
	input <- list(y=y,X=X,EST=EST,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar) {
      y %*% input$X %*% (input$EST + epsilon * VZ) -
        sum(exp(input$X %*% (input$EST + epsilon * VZ))) -
          sum(lgamma(input$y + 1))-input$max_loglik + 0.5 * cstar
    }
    H.lik <- function(y,X,EST){
      y %*% X %*% EST -
        sum(exp(X %*% EST)) -
          sum(lgamma(y+1))
    }
	if (is.null(cov)){
      cov <- vcov(fit)
      list(cov=cov,input=input,H=H,H.lik=H.lik)}
	else {
      list(input=input,H=H,H.lik=H.lik)
    }
  }

