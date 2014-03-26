Mnegbino <-
  function(y,fit,cov) {
    n <- length(y)
	X <- model.matrix(fit)
	EST <- c(as.vector(fit$coef),as.vector(fit$theta))
    p <- length(EST)
    max_loglik <- as.numeric(logLik(fit))
    mu <- fit$fitted.values
	input <- list(y=y,X=X,EST=EST,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar) {
      mu <- exp(input$X %*% (input$EST[1:p-1]+epsilon*VZ[1:p-1]))
      tmp <- input$EST[p]+epsilon*VZ[p]
      sum(lgamma(tmp + input$y) -
          lgamma(tmp) - lgamma(input$y + 1) +
          tmp*log(tmp) +
          input$y * log(mu + (input$y == 0)) -
          (tmp + input$y) * log(tmp + mu)) -
            input$max_loglik+0.5*cstar
    }
	H.lik <- function(y,X,EST) {
      mu <- exp(input$X %*% EST[1:p - 1])
      sum(lgamma(EST[p] + y) - lgamma(EST[p]) - lgamma(y + 1) + 
          EST[p] * log(EST[p]) + y * log(mu + (y == 0)) - (EST[p] + 
                                                 y) * log(EST[p] + mu))
	}
    if (is.null(cov)) {
      info_tmp1 <- matrix(0,p-1,p-1)
      info_tmp2 <- matrix(0,p-1,1)
      for (w in 1:n) {
        info_tmp1 <- info_tmp1 +
          (mu[w]/(1+mu[w]*(1/EST[p])) +
           (1/EST[p])*(y[w] - mu[w])*mu[w]/((1+mu[w]*(1/EST[p]))^2))*(X[w,]%*%t(X[w,]))
        info_tmp2 <- info_tmp2 - ((y[w]-mu[w])*mu[w]/((EST[p]+mu[w])^2))*X[w,]
      }
      cov <- solve(rbind(cbind(info_tmp1,info_tmp2),
                         c(info_tmp2,(1/(fit$SE.theta^2)))))
      list(cov=cov,input=input,H=H,H.lik=H.lik)
    } else {
      list(input=input,H=H,H.lik=H.lik)
    }
  }
