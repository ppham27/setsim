Mlmmc <-
  function(y,fit,cov){
    n <- length(y)
    X <- model.matrix(fit,fit$data)
    groups <- as.character(fit$groups[[1]])
    group <- factor(groups,ordered=FALSE)
    n.group <- length(levels(group))
    EST <- c(diag(getVarCov(fit)),fit$sigma^2)
    p <- length(EST)
    if (is.null(cov)){cov <- fit$apVar}
    covi <- vector("list", n.group)
    rbind(covi)
    Zi <- vector("list", n.group)
    Zi <- split(group,group)
    yi <- vector("list",n.group)
    yi <- split(y,group)
    Xi <- vector("list",n.group)
    Xi <- split(as.data.frame(X),group)
    ri <- vector("list",n.group)
    ri1 <- vector("list",n.group)
    ri2 <- vector("list",n.group)
    lik1 <- 0
    lik2 <- 0
    D <- 0
    H.lik <- function(yi, Xi, EST) {
      for (j in 1:n.group) {
        covi[[j]] <- as.matrix(as.numeric(Zi[[j]])) %*% as.matrix(EST[1:p-1]) %*% t(as.matrix(as.numeric(Zi[[j]])))+diag(EST[p],length(yi[[j]]),length(yi[[j]]))
        D[j] <- log(det(covi[[j]]))
        ri1[[j]] <- t(as.matrix(Xi[[j]])) %*% solve(covi[[j]]) %*% as.matrix(Xi[[j]])
        ri2[[j]] <- t(as.matrix(Xi[[j]])) %*% solve(covi[[j]]) %*% yi[[j]]
      }
      for (j in 1:n.group) {
        ri[[j]] <- yi[[j]]-as.matrix(Xi[[j]]) %*% (solve(Reduce("+",ri1)) %*% Reduce("+",ri2))
        lik1[j] <- t(ri[[j]]) %*% solve(covi[[j]]) %*% ri[[j]]
        lik2[j] <- log(det(t(as.matrix(Xi[[j]])) %*% solve(covi[[j]]) %*% as.matrix(Xi[[j]])))
      }
      lik <- -0.5*(n-p)*log(2*pi)-0.5*sum(D)-0.5*sum(lik1)-0.5*sum(lik2)
      lik
    }
    max_loglik <- fit$logLik
    input <- list(y=yi,X=Xi,EST=EST,max_loglik=max_loglik)
    H <- function(epsilon,input,VZ,cstar){
      H.lik(y=yi,X=Xi,EST=as.numeric(EST+epsilon*VZ))-max_loglik+0.5*cstar
    }
    list(cov=cov,input=input,H=H,H.lik=H.lik)
  }
