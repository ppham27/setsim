Mgamma <-
function(y,fit,cov){
	n <- length(y)
	X <- model.matrix(fit)
  mle_nu <- gamma.shape(fit)$alpha     # the maximum likelihood estimate 
  MLE <- c(as.vector(fit$coef),as.vector(mle_nu))
  p <- length(MLE)
  mu <- fit$fitted.values  
	max_loglik <- sum(-y*mle_nu/mu + (mle_nu-1)*log(y) + mle_nu*log(mle_nu/mu)) - n*lgamma(mle_nu)
	input <- list(y=y,X=X,MLE=MLE,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar){
    if (fit$family$link=="inverse"){
      mu <- 1/((input$MLE[1:p-1]+epsilon*VZ[1:p-1])%*%t(X))
      if (all(mu>0) & input$MLE[p]+epsilon*VZ[p]>0){
        sum( -input$y*(input$MLE[p]+epsilon*VZ[p])/mu + ((input$MLE[p]+epsilon*VZ[p])-1)*log(input$y) + (input$MLE[p]+epsilon*VZ[p])*log((input$MLE[p]+epsilon*VZ[p])/mu) ) - n*(lgamma((input$MLE[p]+epsilon*VZ[p])))  - input$max_loglik+0.5*cstar}
      else {NaN}}
    else if (fit$family$link=="log"){
      mu <- exp((input$MLE[1:p-1]+epsilon*VZ[1:p-1])%*%t(X))
      sum( -input$y*(input$MLE[p]+epsilon*VZ[p])/mu + ((input$MLE[p]+epsilon*VZ[p])-1)*log(input$y) + (input$MLE[p]+epsilon*VZ[p])*log((input$MLE[p]+epsilon*VZ[p])/mu) ) - n*(lgamma((input$MLE[p]+epsilon*VZ[p])))  - input$max_loglik+0.5*cstar}}
  H.lik <- function(y,X,MLE){
    if (fit$family$link=="inverse"){
      mu <- 1/(MLE[1:p-1]%*%t(X))
      if (all(mu>0) & MLE[p]>0){
        sum(-y*MLE[p]/mu + (MLE[p]-1)*log(y) + MLE[p]*log(MLE[p]/mu)) - n*(lgamma(MLE[p]))}
      else {NaN}}
    else if (fit$family$link=="log"){
      mu <- exp(MLE[1:p-1]%*%t(X))
      sum(-y*MLE[p]/mu + (MLE[p]-1)*log(y) + MLE[p]*log(MLE[p]/mu)) - n*(lgamma(MLE[p]))}
  }
  if (is.null(cov)){
    info_tmp1<-matrix(0,p-1,p-1)
    info_tmp2<-matrix(0,p-1,1)
    if (fit$family$link=="inverse"){
      for (w in 1:n){
        info_tmp1 <- info_tmp1 + mle_nu*(mu[w]^2)*(X[w,]%*%t(X[w,]))
        info_tmp2 <- info_tmp2 + (y[w]-mu[w])*X[w,]}
      cov <- solve(rbind(cbind(info_tmp1,info_tmp2), c(info_tmp2,(1/((gamma.shape(fit)$SE)^2)))))}
    else if (fit$family$link=="log"){
      for (w in 1:n){
        info_tmp1 <- info_tmp1 +  mle_nu*(y[w]/mu[w])*(X[w,]%*%t(X[w,]))
        info_tmp2 <- info_tmp2 - (y[w]/mu[w]-1)*X[w,]}
      cov <- solve(rbind(cbind(info_tmp1,info_tmp2), c(info_tmp2,(1/((gamma.shape(fit)$SE)^2)))))}
    list(cov=cov,input=input,H=H,H.lik=H.lik)}
  else {list(input=input,H=H,H.lik=H.lik)}}