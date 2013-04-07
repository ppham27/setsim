independent <- function(model,y,fit,cov,B,max_B) {
  if (missing(cov)) { cov <- NULL }
  if (missing(B)) { B <- NULL }
  if (missing(max_B)) { max_B <- 1000 }

  Model <- eval(call(model, y, fit, cov))
  MLE <- Model$input$MLE
  if (is.null(cov)) { cov <- Model$cov }
  n <- length(y)
  p <- length(MLE)
  e <- eigen(cov)
  V <- e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors)
  max_loglik <- Model$input$max_loglik
  
  independentMini <- function(idx, BB, Model) {
    ## returns output, out_wald, soltype, and z
    input <- Model$input
    H <- function(epsilon,input,VZ,cstar) { Model$H(epsilon,input,VZ,cstar) }
    H.lik <- function(y, X, MLE) { Model$H.lik(y, X, MLE) }
    MLE <- Model$input$MLE
    MLE <- matrix(MLE, 1, length(MLE))
    colnames(MLE) <- paste("MLE", 1:p)

    soltype <- matrix(0, 1, 4)
    colnames(soltype) <- list("zero","one","two",">two")
    
    z <- mvrnorm(BB, rep(0, p), diag(1, nrow = p, ncol = p))
    colnames(z) <- paste("ray.z", 1:p)
    norm <- apply(z, 1, function(r) { sqrt(sum(r*r)) })
    z <- z[order(norm),]
    cstar <- apply(z, 1, function(r) { sum(r*r) })
    alpha <- pchisq(cstar, p)
    
    out_c <- .Call("independentFixed",
                   H, environment(H),
                   H.lik, environment(H.lik),
                   V, z, t(z), cstar, alpha,
                   input, MLE,
                   soltype)
    output <- matrix(out_c[[1]],byrow=TRUE,ncol=p+5)
    colnames(output) <-
      c("z","level","n","epsilon",paste("ray.MLE",1:p),"log_lik")
    out_wald <- matrix(out_c[[2]],byrow=TRUE,ncol=p+4)
    colnames(out_wald) <- c("z", "1-alpha", "n", "epsilon",
                            paste("ray.MLE", 1:p))
    return(list(output=output, out_wald=out_wald, soltype=soltype, z=z))
  }
  
  n.cores <- detectCores()
  if (n.cores >= 4) {
    n.splits <- as.integer(n.cores/2)
  } else {
    n.splits <- n.cores
  }
  cl <- makeCluster(n.splits, outfile="")
  clusterCall(cl, function() {
    library(MASS)
  })
  if (!is.null(B)) {
    out_b <- parLapply(cl, 1:n.splits, independentMini,
                       as.integer(B/n.splits), Model)
  } else {
    ## use stop rule
    diff <- 1
    stop_rule <- 0.05
    batch.size <- 200
    i <- 1
    while (diff > stop_rule) {
      out_tmp <- parLapply(cl, 1:n.splits, independentMini,
                           batch.size, Model)
      ## append result
      if (i==1) {
        out_b <- out_tmp
      } else {
        out_b <- c(out_b, out_tmp)
      }
      ## check two cases for stopping
      
      ## 1.) does it meet criteria?
      out_wald <- do.call(rbind, lapply(out_b, function(l) {l$out_wald}))
      diff_lower <- abs((apply(out_wald[,5:(5+p-1)],2,min) -
                         (MLE-qnorm(0.975)*sqrt(diag(cov))))/
                        (MLE-qnorm(0.975)*sqrt(diag(cov))))
      diff_upper <- abs((apply(out_wald[,5:(5+p-1)],2,max) -
                         (MLE+qnorm(0.975)*sqrt(diag(cov))))/
                        (MLE+qnorm(0.975)*sqrt(diag(cov))))
      diff <- max(diff_lower, diff_upper)
      
      ## 2.) has it exceeded the max number of rays?
      if (i*n.splits*batch.size > p*max_B) {
        ## diff will be less than stop rule
        ## effectively ends loop
        diff <- stop_rule*0.5
      }
      i <- i + 1
    }
  }
  stopCluster(cl)
  ## aggregate output of multiple processes
  output <- do.call(rbind, lapply(out_b, function(l) {l$output}))
  out_wald <- do.call(rbind, lapply(out_b, function(l) {l$out_wald}))
  soltype <- do.call(rbind, lapply(out_b, function(l) {l$soltype}))
  soltype <- as.matrix(aggregate(soltype[,-1],
                                 by=list(level=soltype[,1]),
                                 sum))
  z <- do.call(rbind, lapply(out_b, function(l) {l$z}))
  print("All Done!")

  ## process calculations
  wald.num <- matrix(c(MLE-qnorm(0.975)*sqrt(diag(cov)), MLE+qnorm(0.975)*sqrt(diag(cov))),nrow=p,ncol=2,dimnames=list(paste("MLE",1:p),c("Lower Bound","Upper Bound")))
  wald.sim <- matrix(c(apply(out_wald[,5:(p+4)],2,min), apply(out_wald[,5:(p+4)],2,max)),nrow=p,ncol=2,dimnames=list(paste("MLE",1:p),c("Lower Bound","Upper Bound")))

  list("MLE"=MLE,
       "max_loglik"=max_loglik,
       "diagnosis"=soltype,
       "independent.sample"=output,
       ray.z=z,
       numWald.interval=wald.num,
       simWald.interval=wald.sim)
}
