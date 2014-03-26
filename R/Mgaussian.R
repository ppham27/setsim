Mgaussian <-
  function(y,fit,cov) {
	n <- length(y)
	X <- model.matrix(fit)
	p <- ncol(X)+1
	EST <- c(fit$coef,(t(y-(X%*%fit$coef)) %*% (y-(X%*%fit$coef)))/n)
	max_loglik <- as.numeric(logLik(fit))
	input <- list(y=y,X=X,EST=EST,max_loglik=max_loglik)
	H <- function(epsilon,input,VZ,cstar) {
      tmp <- input$y-input$X %*% (input$EST[1:p-1]+epsilon*VZ[1:(p-1)])
      tmp1 <- input$EST[p]+epsilon*VZ[p]
      -(n/2)*log(2*pi*tmp1)-
        sum(tmp*tmp)/(2*tmp1) -
          input$max_loglik+0.5*cstar
    }
	H.lik <- function(y,X,EST) {
      tmp <- y-X %*% EST[1:p-1]
      -(n/2)*log(2*pi*EST[p])- sum(tmp*tmp)/(2*EST[p])
    }
    if (is.null(cov)) {
      cov <- rbind(cbind(c(EST[p])*(solve(t(X)%*%X)),rep(0,p-1)),
                   c(rep(0,p-1),(2*(EST[p]^2))/n))
      list(cov=cov,input=input,H=H,H.lik=H.lik)
    }
	else {
      list(input=input,H=H,H.lik=H.lik)
    }
  }

