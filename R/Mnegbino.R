Mnegbino <-
  function(y,fit,cov) {
    n <- length(y)
	X <- model.matrix(fit)
	MLE <- c(as.vector(fit$coef),as.vector(fit$theta))
    p <- length(MLE)
    max_loglik <- as.numeric(logLik(fit))
    mu <- fit$fitted.values
	input <- list(y=y,X=X,MLE=MLE,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar) {
      mu <- exp(input$X %*% (input$MLE[1:p-1]+epsilon*VZ[1:p-1]))
      tmp <- input$MLE[p]+epsilon*VZ[p]
      sum(lgamma(tmp + input$y) -
          lgamma(tmp) - lgamma(input$y + 1) +
          tmp*log(tmp) +
          input$y * log(mu + (input$y == 0)) -
          (tmp + input$y) * log(tmp + mu)) -
            input$max_loglik+0.5*cstar
    }
	H.lik <- function(y,X,MLE) {
      mu <- exp(input$X %*% MLE[1:p - 1])
      sum(lgamma(MLE[p] + y) - lgamma(MLE[p]) - lgamma(y + 1) + 
          MLE[p] * log(MLE[p]) + y * log(mu + (y == 0)) - (MLE[p] + 
                                                 y) * log(MLE[p] + mu))
	}
    if (is.null(cov)) {
      info_tmp1 <- matrix(0,p-1,p-1)
      info_tmp2 <- matrix(0,p-1,1)
      for (w in 1:n) {
        info_tmp1 <- info_tmp1 +
          (mu[w]/(1+mu[w]*(1/MLE[p])) +
           (1/MLE[p])*(y[w] - mu[w])*mu[w]/((1+mu[w]*(1/MLE[p]))^2))*(X[w,]%*%t(X[w,]))
        info_tmp2 <- info_tmp2 - ((y[w]-mu[w])*mu[w]/((MLE[p]+mu[w])^2))*X[w,]
      }
      cov <- solve(rbind(cbind(info_tmp1,info_tmp2),
                         c(info_tmp2,(1/(fit$SE.theta^2)))))
      list(cov=cov,input=input,H=H,H.lik=H.lik)
    } else {
      list(input=input,H=H,H.lik=H.lik)
    }
  }
