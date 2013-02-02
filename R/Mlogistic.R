Mlogistic <-
function(y,fit,cov){
	n <- length(y)
	X <- model.matrix(fit)
	MLE <- as.vector(fit$coef)
	pi <- exp(MLE%*%t(X))/(1+exp(MLE%*%t(X)))
	max_loglik<-log(pi)%*%y + log(1-pi)%*%(1-y)
	input <- list(y=y,X=X,MLE=MLE,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar){
		t(input$y) %*% t(log(exp(t(input$MLE+epsilon*VZ)%*%t(input$X))/(1+exp(t(input$MLE+epsilon*VZ)%*%t(input$X)))))+t(1-input$y) %*% t(log(1/(1+exp(t(input$MLE+epsilon*VZ)%*%t(input$X)))))-input$max_loglik+0.5*cstar}
	H.lik <- function(y,X,MLE){t(y) %*% t(log(exp(MLE%*%t(X))/(1+exp(MLE%*%t(X)))))+t(1-y) %*% t(log(1/(1+exp(MLE%*%t(X)))))}
  if (is.null(cov)){
		cov <- vcov(fit)
		list(cov=cov,input=input,H=H,H.lik=H.lik)}
	else {list(input=input,H=H,H.lik=H.lik)}}

