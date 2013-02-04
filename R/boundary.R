boundary <-
  function (model,y,fit,target,targetvalue,cov,B,max_B,...) {

    if (missing(cov)) { cov <- NULL }
    if (missing(B)) { B <- NULL }
    if (missing(max_B)) { max_B <- 1000 }

    is.whole <- function(value) { floor(value)==value }
    Model <- eval(call(model, y, fit, cov))
    input <- Model$input
    if (is.null(cov)){cov <- Model$cov}
    H <- Model$H
    MLE <- input$MLE

    n <- length(y)
    p <- length(MLE)
    e <- eigen(cov)
    V <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
    if (target=="level") { critvalue <- qchisq(targetvalue,p) }
    if (target=="ratio") { critvalue <- -2*log(targetvalue) }
    if (target=="customized") { critvalue <- cvalue }

    MLE <- matrix(MLE, 1, length(MLE))
    colnames(MLE) <- paste("MLE", 1:p)
    soltype <- matrix(0, length(targetvalue), 4)
    soltype <- cbind(targetvalue, soltype)
    colnames(soltype) <- list(paste(target),"zero","one","two",">two")
    output <- list()
    out_wald <- list()
    
    if (!is.null(B)) {
      z <- mvrnorm(B, rep(0, p), diag(1, nrow = p, ncol = p))
      colnames(z) <- paste("ray.z", 1:p)
      norm <- 0
      for (i in 1:dim(z)[1]) { norm[i] <- sqrt(sum(z[i, ]^2)) }
      z <- z[order(norm),]
      ## start here
      ## out_wald <- .Call("boundary",
      ##                   model, V, z, critvalue, soltype, targetvalue, input, MLE)
      for (i in 1:dim(z)[1]) {
        VZ <- V %*% as.vector(t(z[i, ]))
        for (s in 1:length(critvalue)){
          c_value <- critvalue[s]
          a <- sqrt(critvalue[s])/sqrt(sum(z[i, ]^2))
          epsilon <- seq(-4*a, 4*a, by=(8*a)/(51-1)) 
          H.values <- NA
          for (j in 1:length(epsilon)) { H.values[j] <- H(epsilon[j],input,VZ,c_value) }
          H.final <- function(...){ H.values }
          roots <- nroots(H.final, epsilon)
          range <- roots$range
          sol <- 0
          if (roots$n == 0){
            soltype[s, 2] <- soltype[s, 2] + 1
          } else {
            sol <- sapply(1:roots$n,
                          function(k) {
                            bisect(fn=H,a=range[k,1],b=range[k,2],
                                   input=input,VZ=VZ,c_value=c_value)
                          })
            for (k in 1:roots$n) {
              output[[length(output)+1]] <-
                c(i,targetvalue[s],roots$n,sol[k],as.vector(MLE+sol[k]*VZ[,1]))
            }
            if (roots$n == 1) {
              soltype[s, 3] <- soltype[s, 3] + 1
            } else if (roots$n == 2) {
              soltype[s,4] <- soltype[s, 4] + 1
            } else {
              soltype[s, 5] <- soltype[s, 5] + 1
            }
          }
        }
        sol_w <- sqrt(qchisq(0.95,1))/sqrt(sum(z[i,]^2))
        out_wald[[length(out_wald)+1]] <-
          c(i, targetvalue[s], roots$n, -sol_w, as.vector(MLE - sol_w * VZ[,1]))
        out_wald[[length(out_wald)+1]] <-
          c(i, targetvalue[s], roots$n, sol_w, as.vector(MLE + sol_w * VZ[,1]))
        if (is.whole(i/200)==TRUE) { print(paste("Generated",i,"rays.")) }
      }
      output <- do.call(rbind, output)
      colnames(output) <- c("z","level","n","epsilon",paste("ray.MLE",1:p))
      out_wald <- do.call(rbind, out_wald)
      colnames(out_wald) <- c("z", "1-alpha", "n", "epsilon",
                              paste("ray.MLE", 1:p))
      ## end here
    } else {
      soltmp_w <- sqrt(qchisq(0.95, 1))
      out_z <- matrix(numeric(0), 0, p)
      diff <- 1
      stop_rule <- 0.05
      i <- 1
      while (diff > stop_rule){
        z <- mvrnorm(1,rep(0,p),diag(1,nrow=p,ncol=p))
        out_z <- rbind(out_z,z)
        VZ <- V %*% as.vector(t(z))
        for (s in 1:length(critvalue)){
          c_value <- critvalue[s]
          a <- sqrt(critvalue[s])/sqrt(sum(z^2))
          epsilon <- seq(-4*a, 4*a, by=(8*a)/(51-1))
          H.values <- NA
          for (j in 1:length(epsilon)) {
            H.values[j] <- H(epsilon[j], input, VZ, c_value)
          }
          H.final <- function(...) { H.values }
          roots <- nroots(H.final, epsilon)
          range <- roots$range
          sol <- 0
          if (roots$n==0) {
            soltype[s, 2] <- soltype[s, 2] + 1
          } else {
            for (k in 1:roots$n){
              sol[k] <-
                bisect(fn=H,a=range[k,1],b=range[k,2],
                       input=input,VZ=VZ,c_value=c_value)
              output <-
                rbind(output,
                      c(i,targetvalue[s],roots$n,sol[k],
                        as.vector(MLE+sol[k]*VZ[, 1])))
            }
            if (roots$n==1){soltype[s,3] <- soltype[s,3]+1}
            else if (roots$n==2){soltype[s,4] <- soltype[s,4]+1}
            else {soltype[s,5] <- soltype[s,5]+1}}}
        if (is.whole(i/200)==TRUE){print(paste("Generated",i,"rays."))}
        i <- i + 1
        sol_w <- soltmp_w/sqrt(sum(z^2))
        out_wald <-
          rbind(out_wald,
                c(i,targetvalue[s],roots$n,-sol_w,as.vector(MLE-sol_w*VZ[,1])))
        out_wald <-
          rbind(out_wald,
                c(i,targetvalue[s],roots$n,sol_w,as.vector(MLE+sol_w*VZ[,1])))
        if (i%%100==0 && i >= 50*p) {
          diff_lower <- abs((apply(out_wald[,5:(5+p-1)],2,min) -
                             (MLE-qnorm(0.975)*sqrt(diag(cov))))/
                            (MLE-qnorm(0.975)*sqrt(diag(cov))))
          diff_upper <- abs((apply(out_wald[,5:(5+p-1)],2,max) -
                             (MLE+qnorm(0.975)*sqrt(diag(cov))))/
                            (MLE+qnorm(0.975)*sqrt(diag(cov))))
          diff <- max(diff_lower, diff_upper)
        }
        else if (i>p*max_B){diff <- stop_rule*0.5}
      }
      z <- out_z
      colnames(z) <- paste("ray.z", 1:p)
    }
    print("Done!")
    conv.waldnum <- matrix(c(MLE-qnorm(0.975)*sqrt(diag(cov)),MLE+qnorm(0.975)*sqrt(diag(cov))),p,2,dimnames=list(paste("MLE",1:p),c("Lower Bound","Upper Bound")))
    conv.waldsim <- matrix(c(apply(out_wald[,5:(p+4)],2,min),apply(out_wald[,5:(p+4)],2,max)),p,2,dimnames=list(paste("MLE",1:p),c("Lower Bound","Upper Bound")))
    if (target=="level" || target=="customized") {
      wald.num <- matrix(NA,length(targetvalue),2*p+1)
      lik.sim <- matrix(NA, length(targetvalue), 2*p+1)
      for (i in 1:length(targetvalue)) {
        wald.num[i,1] <- pchisq(qchisq(targetvalue[i], p),1)
        wald.num[i,2:(2*p+1)] <-
          matrix(c(MLE-sqrt(qchisq(targetvalue[i],p))*sqrt(diag(cov)),
                   MLE+sqrt(qchisq(targetvalue[i],p))*sqrt(diag(cov))),
                 1,2*p)
        tmp <- output[output[,2]==targetvalue[i],]
        lik.sim[i,1] <- pchisq(qchisq(targetvalue[i],p),1)
        lik.sim[i,2:(2*p+1)] <- matrix(c(apply(tmp[,5:(p+4)],2,min),apply(tmp[,5:(p+4)],2,max)),nrow=1,ncol=2*p)
      }
      rownames(wald.num) <- targetvalue
      colnames(wald.num) <- c("1D.Level",paste("Lower.MLE",1:p),paste("Upper.MLE",1:p))
      rownames(lik.sim) <- targetvalue
      colnames(lik.sim) <- c("1D.Level",paste("Lower.MLE",1:p),paste("Upper.MLE",1:p))
      list(MLE = MLE, diagnosis = soltype, boundary.sample = output,
           WaldBoundary.sample = out_wald,
           ray.z = z, numWald.interval = wald.num,
           simLik.interval = lik.sim, convnumWald = conv.waldnum, convsimuWald = conv.waldsim)
    } else if (target=="ratio") {
      list(MLE = MLE, diagnosis = soltype, boundary.sample = output,
           boundaryWald.sample = out_wald, ray.z = z,
           convnumWald = conv.waldnum, convsimuWald = conv.waldsim)
    }
  }
