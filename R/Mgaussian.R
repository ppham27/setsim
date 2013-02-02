Mgaussian <-
function(y,fit,cov){
	n <- length(y)
	X <- model.matrix(fit)
	p <- ncol(X)+1
	MLE <- c(fit$coef,(t(y-(X%*%fit$coef)) %*% (y-(X%*%fit$coef)))/n)
	max_loglik <- as.numeric(logLik(fit))
	input <- list(y=y,X=X,MLE=MLE,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar){
    -(n/2)*log(2*pi*(input$MLE[p]+epsilon*VZ[p,1]))-(t(input$y-input$X %*% (input$MLE[1:p-1]+epsilon*VZ[1:(p-1),1])) %*% (input$y-input$X %*% (input$MLE[1:p-1]+epsilon*VZ[1:(p-1),1])))/(2*(input$MLE[p]+epsilon*VZ[p,1]))-input$max_loglik+0.5*cstar}
	H.lik <- function(y,X,MLE){-(n/2)*log(2*pi*MLE[p])-(t(y-X %*% MLE[1:p-1]) %*% (y-X %*% MLE[1:p-1]))/(2*MLE[p])}
  if (is.null(cov)){
		cov <- rbind(cbind(c(MLE[p])*(solve(t(X)%*%X)),rep(0,p-1)), c(rep(0,p-1),(2*(MLE[p]^2))/n))
		list(cov=cov,input=input,H=H,H.lik=H.lik)}
	else {list(input=input,H=H,H.lik=H.lik)}}

