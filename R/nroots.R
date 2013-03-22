nroots <-
  function(fn,x) {
    ## Initialize Variables
	fn1 <- function(x){fn(x,...)}
	a <- min(x)
	b <- max(x)
	interval <- length(x)
	sign <- matrix(9,length(x),1)
	sign.lag <- matrix(9,length(x),1)
	n <- 0
	Ans <- 0
	
    ## Search Solutions
	index.root <- which(fn(x)==0)
	index.neg <- which(fn(x)<0)
	index.pos <- which(fn(x)>0)
	sign[index.root] <- 0
	sign[index.neg] <- 1
	sign[index.pos] <- 2
	sign.lag[2:interval] <- sign[1:(interval-1)]
	sign.sum <- sign.lag+sign
	index <- which(sign.sum==3)
	index <- c(index,index.root)
	index <- index[order(index)]
	n <- length(index)
	sol <- matrix(c(x[index-1],x[index],fn(x)[index-1],fn(x)[index]),ncol=4,nrow=n,byrow=FALSE,
	              dimnames=list(rep("Interval",n),c("x.lower","x.upper","y.lower","y.upper")))
	list(n=n,position=index,range=sol)
  }

