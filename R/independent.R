independent <-
  function(model,y,fit,cov,B,max_B){

    if (missing(cov)){cov <- NULL}
    if (missing(B)){B <- NULL}
    if (missing(max_B)){max_B <- 1000}

    is.whole <- function(value) { floor(value)==value }
    Model <- eval(call(model,y,fit,cov))
    input <- Model$input
    if (is.null(cov)){cov <- Model$cov}
    H <- Model$H
    H.lik <- Model$H.lik
    X <- input$X
    MLE <- input$MLE
    max_loglik <- input$max_loglik

    n <- length(y)
    p <- length(MLE)
    e <- eigen(cov)
    V <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors) #Matrix sqrt of V
    alpha <- 0

    MLE <- matrix(MLE,1,length(MLE))
    colnames(MLE) <- paste("MLE", 1:p)
    soltype <- matrix(0,1,4)
    colnames(soltype) <- c("zero","one","two",">two")
    output <- matrix(numeric(0),0,p+5)
    output <- matrix(output,nrow=nrow(output),ncol=ncol(output),dimnames=list(rep("",nrow(output)),c("z","level","n","epsilon",paste("ray.MLE",1:p),"log_lik")))
    out_wald <- matrix(numeric(0), 0, p + 4)
    out_wald <- matrix(out_wald, nrow = nrow(out_wald), ncol = ncol(out_wald), dimnames = list(rep("", nrow(out_wald)), c("z", "1-alpha", "n", "epsilon", paste("ray.MLE", 1:p))))

    if (!is.null(B)){
      z <- mvrnorm(B,rep(0,p),diag(1,nrow=p,ncol=p)) #Generate B zs from N(0,Ip)
      colnames(z) <- paste("ray.z", 1:p)
      norm <- 0
      for (i in 1:dim(z)[1]){norm[i] <- sqrt(sum(z[i,]^2))}
      z <- z[order(norm),]
      for (i in 1:dim(z)[1]){
        VZ <- V %*% as.vector(t(z[i,]))
        cstar <- t(z[i,])%*%z[i,]
        alpha[i] <- pchisq(cstar,p) #Find alpha(z)
        a <- sqrt(cstar)/sqrt(sum(z[i, ]^2))
        epsilon <- seq(-4*a, 4*a, by=(8*a)/(51-1))
        H.values <- NA
        for (j in 1:length(epsilon)){H.values[j] <- H(epsilon[j],input,VZ,cstar)}
        H.final <- function(...){H.values}
        roots <- nroots(H.final,epsilon)
        range <- roots$range
        sol <- 0
        if (roots$n==0){
          output <- rbind(output,c(i,alpha[i],roots$n,rep(NA,p+1)))
          soltype[,1] <- soltype[,1] + 1}
        else{
          for (k in 1:roots$n){
            sol[k] <- bisect(fn=H,a=range[k,1],b=range[k,2],input=input,VZ=VZ,c_value=cstar)
            output <- rbind(output,c(i,alpha[i],roots$n,sol[k],as.vector(MLE+sol[k]*VZ[,1]),H.lik(y=input$y,X=input$X,MLE=as.vector(input$MLE+sol[k]*VZ[,1]))))}
          if (roots$n==1){soltype[,2] <- soltype[,2] + 1}
          else if (roots$n==2){soltype[,3] <- soltype[,3] + 1}
          else{soltype[,4] <- soltype[,4] + 1}}
        sol_w <- 1
        out_wald <- rbind(out_wald, c(i, cstar, roots$n, -sol_w, as.vector(MLE - sol_w * VZ[,1])))
        out_wald <- rbind(out_wald, c(i, cstar, roots$n, sol_w, as.vector(MLE + sol_w * VZ[,1])))
        if (is.whole(i/200)==TRUE){print(paste("Generated",i,"rays."))}
      }
    }

    else{
      out_z <- matrix(numeric(0), 0, p)
      diff <- 1
      stop_rule <- 0.05
      i <- 1
      while (diff > stop_rule){
        z <- mvrnorm(1,rep(0,p),diag(1,nrow=p,ncol=p))
        out_z <- rbind(out_z,z)
        VZ <- V %*% as.vector(t(z))
        cstar <- t(z)%*%z
        alpha <- pchisq(cstar,p,lower.tail=FALSE) #Find alpha(z)
        a <- sqrt(cstar)/sqrt(sum(z^2))
        epsilon <- seq(-4*a, 4*a, by=(8*a)/(51-1))
        H.values <- NA
        for (j in 1:length(epsilon)){H.values[j] <- H(epsilon[j],input,VZ,cstar)}
        H.final <- function(...){H.values}
        roots <- nroots(H.final,epsilon)
        range <- roots$range
        sol <- 0
        if (roots$n==0){
          output <- rbind(output,c(i,alpha,roots$n,rep(NA,p+1)))
          soltype[,1] <- soltype[,1] + 1}
        else{
          for (k in 1:roots$n){
            sol[k] <- bisect(fn=H,a=range[k,1],b=range[k,2],input=input,VZ=VZ,c_value=cstar)
  			output <- rbind(output,c(i,alpha,roots$n,sol[k],as.vector(MLE+sol[k]*VZ[,1]),H.lik(y=input$y,X=input$X,MLE=as.vector(input$MLE+sol[k]*VZ[,1]))))}
          if (roots$n==1){soltype[,2] <- soltype[,2] + 1}
          else if (roots$n==2){soltype[,3] <- soltype[,3] + 1}
          else{soltype[,4] <- soltype[,4] + 1}}
        if (is.whole(i/200)==TRUE){print(paste("Generated",i,"rays."))}
        i <- i + 1
        sol_w <- 1
        out_wald <- rbind(out_wald, c(i, cstar, roots$n, -sol_w, as.vector(MLE - sol_w * VZ[,1])))
        out_wald <- rbind(out_wald, c(i, cstar, roots$n, sol_w, as.vector(MLE + sol_w * VZ[,1])))
        if (i%%100==0 && i >= 50*p) {
          diff_lower <- abs((apply(out_wald[,5:(5+p-1)],2,min)-(MLE-qnorm(0.975)*sqrt(diag(cov))))/(MLE-qnorm(0.975)*sqrt(diag(cov))))
          diff_upper <- abs((apply(out_wald[,5:(5+p-1)],2,max)-(MLE+qnorm(0.975)*sqrt(diag(cov))))/(MLE+qnorm(0.975)*sqrt(diag(cov))))
          diff <- max(diff_lower, diff_upper)
        }
        else if (i>p*max_B){diff <- stop_rule*0.5}
      }
      z <- out_z
      colnames(z) <- paste("ray.z", 1:p)
    }
    print("Done!")
    wald.num <- matrix(c(MLE-qnorm(0.975)*sqrt(diag(cov)), MLE+qnorm(0.975)*sqrt(diag(cov))),nrow=p,ncol=2,dimnames=list(paste("MLE",1:p),c("Lower Bound","Upper Bound")))
    wald.sim <- matrix(c(apply(out_wald[,5:(p+4)],2,min), apply(out_wald[,5:(p+4)],2,max)),nrow=p,ncol=2,dimnames=list(paste("MLE",1:p),c("Lower Bound","Upper Bound")))
    list("MLE"=MLE,"max_loglik"=max_loglik,"diagnosis"=soltype,"independent.sample"=output,ray.z=z,numWald.interval=wald.num,simWald.interval=wald.sim)}

