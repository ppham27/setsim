benchmark <- function() {
  ## gaussian benchmark
  data(data.gaussian)
  y <- data.gaussian$taste
  x1 <- data.gaussian$h2s
  x2 <- data.gaussian$lactic
  n <- length(y)
  fit <- lm(y~(x1+x2))
  target <- "level"
  targetvalue <- c(0.5,0.9)
  start <- Sys.time()
  out_b <- boundary("Mgaussian",y,fit,target,targetvalue,B=2000, n.splits=1)
  print(paste("Gaussian, 1 process, time taken (s):",Sys.time()-start))
  start <- Sys.time()
  out_b <- boundary("Mgaussian",y,fit,target,targetvalue,B=2000, n.splits=2)
  print(paste("Gaussian, 2 processes, time taken (s):",Sys.time()-start))
  start <- Sys.time()
  out_b <- boundary("Mgaussian",y,fit,target,targetvalue,B=2000, n.splits=4)
  print(paste("Gaussian, 4 processes, time taken (s):",Sys.time()-start))

  ## gamma benchmark
  data(data.loggamma)
  x <- log(data.loggamma$Brain)
  y <- data.loggamma$Body
  fit <- glm(y ~ x, family = Gamma(link = "log"))
  target <- "level"
  targetvalue <- c(0.5,0.9)
  out_b <- boundary("Mgamma",y,fit,target,targetvalue,B=2000, n.splits=1)
  print(paste("Gamma, 1 process, time taken (s):",Sys.time()-start))
  start <- Sys.time()
  out_b <- boundary("Mgamma",y,fit,target,targetvalue,B=2000, n.splits=2)
  print(paste("Gamma, 2 processes, time taken (s):",Sys.time()-start))
  start <- Sys.time()
  out_b <- boundary("Mgamma",y,fit,target,targetvalue,B=2000, n.splits=4)
  print(paste("Gamma, 4 processes, time taken (s):",Sys.time()-start))
  
  ## logistic benchmark
  data(data.logistic)
  y <- data.logistic$disease
  fit <- glm(y ~ (data.logistic$age+data.logistic$sector), family=binomial)
  target <- "level"
  targetvalue <- c(0.5,0.9)
  start <- Sys.time()
  out_b <- boundary("Mlogistic",y,fit,target,targetvalue,B=2000, n.splits=1)
  print(paste("Logistic, 1 process, time taken (s):",Sys.time()-start))
  start <- Sys.time()
  out_b <- boundary("Mlogistic",y,fit,target,targetvalue,B=2000, n.splits=2)
  print(paste("Logistic, 2 processes, time taken (s):",Sys.time()-start))
  start <- Sys.time()
  out_b <- boundary("Mlogistic",y,fit,target,targetvalue,B=2000, n.splits=4)
  print(paste("Logistic, 4 processes, time taken (s):",Sys.time()-start))

  ## poisson benchmark
  data(data.poisson)
  Device <- data.poisson$Device
  DataPlan <- data.poisson$DataPlan
  y <- data.poisson$y
  fit <- glm(y ~ (Device+DataPlan), family=poisson)
  target <- "level"
  targetvalue <- c(0.5,0.9)
  start <- Sys.time()
  out_b <- boundary("Mpoisson",y,fit,target,targetvalue,B=2000, n.splits=1)
  print(paste("Poisson, 1 process, time taken (s):",Sys.time()-start))
  start <- Sys.time()
  out_b <- boundary("Mpoisson",y,fit,target,targetvalue,B=2000, n.splits=2)
  print(paste("Poisson, 2 processes, time taken (s):",Sys.time()-start))
  start <- Sys.time()
  out_b <- boundary("Mpoisson",y,fit,target,targetvalue,B=2000, n.splits=4)
  print(paste("Poisson, 4 processes, time taken (s):",Sys.time()-start))
}
