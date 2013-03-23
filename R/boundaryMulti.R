boundaryMulti <-
  function (model,y,fit,target,targetvalue,cov,B,max_B,...) {

    if (missing(cov)) { cov <- NULL }
    if (missing(B)) { B <- NULL }
    if (missing(max_B)) { max_B <- 1000 }
    is.whole <- function(value) { floor(value)==value }

    Model <- eval(call(model, y, fit, cov))
    MLE <- Model$input$MLE
    if (is.null(cov)){cov <- Model$cov}
    n <- length(y)
    p <- length(MLE)
    e <- eigen(cov)
    V <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
    if (target=="level") { critvalue <- qchisq(targetvalue,p) }
    if (target=="ratio") { critvalue <- -2*log(targetvalue) }
    if (target=="customized") { critvalue <- cvalue }

    boundaryMini <- function(idx, BB, Model) {
      library(MASS)
      input <- Model$input
      H <- function(epsilon,input,VZ,cstar) { Model$H(epsilon,input,VZ,cstar) }
      MLE <- Model$input$MLE
      MLE <- matrix(MLE, 1, length(MLE))
      colnames(MLE) <- paste("MLE", 1:p)
      soltype <- matrix(0, length(targetvalue), 4)
      soltype <- cbind(targetvalue, soltype)
      colnames(soltype) <- list(paste(target),"zero","one","two",">two")    
      
      z <- mvrnorm(BB, rep(0, p), diag(1, nrow = p, ncol = p))
      colnames(z) <- paste("ray.z", 1:p)
      norm <- 0
      for (i in 1:dim(z)[1]) { norm[i] <- sqrt(sum(z[i, ]^2)) }
      z <- z[order(norm),]
      out_c <- .Call("boundaryFixed",
                     H, environment(H),
                     V, z, t(z), critvalue,
                     targetvalue, input, MLE,
                     soltype)
      out_c[[length(out_c)+1]] <- soltype
      output <- matrix(out_c[[1]],byrow=TRUE,ncol=p+4)
      colnames(output) <- c("z","level","n","epsilon",paste("ray.MLE",1:p))
      out_wald <- matrix(out_c[[2]],byrow=TRUE,ncol=p+4)
      colnames(out_wald) <- c("z", "1-alpha", "n", "epsilon",
                              paste("ray.MLE", 1:p))
      return(list(output=output, out_wald=out_wald, soltype=soltype,z=z))
    }
    n.cores <- detectCores()
    ## out_b <- mclapply(1:n.cores, boundaryMini, as.integer(B/n.cores), Model,
    ##                   mc.cores=n.cores)
    cl <- makeCluster(n.cores)
    out_b <- parLapply(cl, 1:n.cores, boundaryMini,
                       as.integer(B/n.cores), Model)
    stopCluster(cl)
    output <- do.call(rbind, lapply(out_b, function(l) {l$output}))
    out_wald <- do.call(rbind, lapply(out_b, function(l) {l$out_wald}))
    soltype <- do.call(rbind, lapply(out_b, function(l) {l$soltype}))
    soltype <- as.matrix(aggregate(soltype[,-1],
                                   by=list(level=soltype[,1]),
                                   sum))
    z <- do.call(rbind, lapply(out_b, function(l) {l$z}))
    print("All Done!")
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
