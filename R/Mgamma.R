Mgamma <-
  function(y,fit,cov) {
	n <- length(y)
	X <- model.matrix(fit)
    est_nu <- gamma.shape(fit)$alpha     # the maximum likelihood estimate 
    EST <- c(as.vector(fit$coef),as.vector(est_nu))
    p <- length(EST)
    mu <- fit$fitted.values  
	max_loglik <- sum(-y*est_nu/mu + (est_nu-1)*log(y) + est_nu*log(est_nu/mu)) - n*lgamma(est_nu)
	input <- list(y=y,X=X,EST=EST,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar) {
      if (fit$family$link=="inverse") {
        mu <- 1/(X %*% (input$EST[1:p-1]+epsilon*VZ[1:p-1]))
        if (all(mu>0) & input$EST[p]+epsilon*VZ[p]>0) {
          tmp <- input$EST[p]+epsilon*VZ[p]
          sum(-input$y*tmp/mu +
              (tmp-1)*log(input$y) +              
              tmp*log(tmp/mu)) -
                n*lgamma(tmp) -
                  input$max_loglik+0.5*cstar
        } else {
          NaN
        }
      } else if (fit$family$link=="log") {
        mu <- exp(X %*% (input$EST[1:p-1]+epsilon*VZ[1:p-1]))
        tmp <- input$EST[p]+epsilon*VZ[p]
        sum(-input$y*tmp/mu +
            (tmp-1)*log(input$y) +
            tmp*log(tmp/mu)) -
              n*lgamma(tmp)  -
                input$max_loglik+0.5*cstar
      }
    }
    H.lik <- function(y,X,EST) {
      if (fit$family$link=="inverse") {
        mu <- 1/(X %*% EST[1:p-1])
        if (all(mu>0) & EST[p]>0) {
          sum(-y*EST[p]/mu + (EST[p]-1)*log(y) + EST[p]*log(EST[p]/mu)) -
            n*lgamma(EST[p])
        } else {
          NaN
        }
      } else if (fit$family$link=="log") {
        mu <- exp(X %*% EST[1:p-1])
        sum(-y*EST[p]/mu + (EST[p]-1)*log(y) + EST[p]*log(EST[p]/mu)) -
          n*lgamma(EST[p])
      }
    }
    if (is.null(cov)) {
      info_tmp1<-matrix(0,p-1,p-1)
      info_tmp2<-matrix(0,p-1,1)
      if (fit$family$link=="inverse") {
        for (w in 1:n) {
          info_tmp1 <- info_tmp1 + est_nu*(mu[w]^2)*(X[w,]%*%t(X[w,]))
          info_tmp2 <- info_tmp2 + (y[w]-mu[w])*X[w,]
        }
        cov <- solve(rbind(cbind(info_tmp1,info_tmp2),
                           c(info_tmp2,(1/((gamma.shape(fit)$SE)^2)))))
      } else if (fit$family$link=="log") {
        for (w in 1:n){
          info_tmp1 <- info_tmp1 +  est_nu*(y[w]/mu[w])*(X[w,]%*%t(X[w,]))
          info_tmp2 <- info_tmp2 - (y[w]/mu[w]-1)*X[w,]}
        cov <- solve(rbind(cbind(info_tmp1,info_tmp2),
                           c(info_tmp2,(1/((gamma.shape(fit)$SE)^2)))))
      }
      list(cov=cov,input=input,H=H,H.lik=H.lik)
    } else {
      list(input=input,H=H,H.lik=H.lik)
    }
  }
