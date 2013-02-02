HLCR <-
function(conf,sample,max_loglik,decimal){
  if (missing(decimal)){decimal <- 4}
  p <- ncol(sample)
  n <- nrow(sample)
  c <- 0.5
  limit <- rep((log(c)+max_loglik),n)
  tol <- sum(sample[,(p)]>=limit)/n
  if (tol<conf){
    while (tol<conf){
      c <- c-10^(-1*decimal)
      limit <- rep((log(c)+max_loglik),n)
      tol <- sum(sample[,p]>=limit)/n
    }
  }
  if (tol>conf){
    while (tol>conf){
      c <- c+10^(-1*decimal)
      limit <- rep((log(c)+max_loglik),n)
      tol <- sum(sample[,p]>=limit)/n
    }
  }
  inside <- sample[(sample[,p]>=limit),]
  colnames(inside) <- c(paste(rep("theta",p-1), 1:(p-1)),"log_lik")
  outside <- sample[(sample[,p]<limit),]
  colnames(outside) <- c(paste(rep("theta",p-1),1:(p-1)),"log_lik")
  list(inside=inside, outside=outside, c=c)
}