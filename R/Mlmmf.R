Mlmmf <-
  function (y, fit, cov) {
    n <- length(y)
    X <- model.matrix(fit, fit$data)
    MLE <- as.vector(fit$coef$fixed)
    groups <- as.character(fit$groups[[1]])
    group <- factor(groups, ordered = FALSE)
    n.group <- length(levels(group))
    yi <- vector("list", n.group)
    yi <- split(y, group)
    Xi <- vector("list", n.group)
    Xi <- split(as.data.frame(X), group)
    if (is.null(cov)){cov <- vcov(fit)}
    covi <- vector("list", n.group)
    inv_covi <- covi
    D <- 0
    like1 <- vector("list", n.group)
    like2 <- vector("list", n.group)
    like3 <- vector("list", n.group)
    like4 <- vector("list", n.group)
    for (j in 1:n.group) {
      covi[[j]] <- getVarCov(fit,individuals=levels(group)[j],type="marginal")
      D[j] <- log(det(covi[[j]][[1]]))
      inv_covi[[j]][[1]] <- solve(covi[[j]][[1]])
      like1[[j]] <- t(yi[[j]]) %*% inv_covi[[j]][[1]] %*% yi[[j]]
      like2[[j]] <- t(yi[[j]]) %*% inv_covi[[j]][[1]] %*% as.matrix(Xi[[j]])
      like3[[j]] <- t(as.matrix(Xi[[j]])) %*% inv_covi[[j]][[1]] %*% yi[[j]]
      like4[[j]] <- t(as.matrix(Xi[[j]])) %*% inv_covi[[j]][[1]] %*% as.matrix(Xi[[j]])
    }
    max_loglik <- fit$logLik
    lik1 <- Reduce("+",like1)
    lik2 <- Reduce("+",like2)
    lik3 <- Reduce("+",like3)
    lik4 <- Reduce("+",like4)
    input <- list(y=yi,X=Xi,MLE=MLE,max_loglik=max_loglik)
    H <- function(epsilon,input,VZ,cstar){
      -0.5*n*log(2*pi)-0.5*sum(D)-0.5*(lik1-lik2 %*% (MLE+epsilon*VZ)-t((MLE+epsilon*VZ)) %*% lik3 + t((MLE+epsilon*VZ)) %*% lik4 %*% ((MLE+epsilon*VZ)))-max_loglik+0.5*cstar
    }
    H.lik <- function(yi,Xi,MLE) {
      -0.5*n*log(2*pi)-0.5*sum(D)-0.5*(lik1-lik2 %*% MLE-t(MLE) %*% lik3 + t(MLE) %*% lik4 %*% MLE)
    }
    list(cov=cov,input=input,H=H,H.lik=H.lik)
  }
