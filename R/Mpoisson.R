Mpoisson <-
  function(y,fit,cov) {
	n <- length(y)
	X <- model.matrix(fit)
	MLE <- as.vector(fit$coef)
	max_loglik <- y %*% X %*% MLE -
      sum(exp(X %*% MLE)) -
        sum(lgamma(y+1))
	input <- list(y=y,X=X,MLE=MLE,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar) {
      y %*% input$X %*% (input$MLE + epsilon * VZ) -
        sum(exp(input$X %*% (input$MLE + epsilon * VZ))) -
          sum(lgamma(input$y + 1))-input$max_loglik + 0.5 * cstar
    }
    H.lik <- function(y,X,MLE){
      y %*% X %*% MLE -
        sum(exp(X %*% MLE)) -
          sum(lgamma(y+1))
    }
	if (is.null(cov)){
      cov <- vcov(fit)
      list(cov=cov,input=input,H=H,H.lik=H.lik)}
	else {
      list(input=input,H=H,H.lik=H.lik)
    }
  }

