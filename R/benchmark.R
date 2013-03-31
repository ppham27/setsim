setsimBenchmark <- function(model="Mpoisson") {
  if (missing(model) | model=="Mpoisson") {
    data(data.poisson)
    Device <- data.poisson$Device
    DataPlan <- data.poisson$DataPlan
    y <- data.poisson$y
    fit <- glm(y ~ (Device+DataPlan), family=poisson)
    target <- "level"
    targetvalue <- c(0.5,0.9)    # targeted confidence level
    Rprof(tmp <- tempfile())
    ## Simulation with desired number of rays
    out_b <- boundary("Mpoisson",y,fit,target,targetvalue,B=1000)   
    Rprof()
    summaryRprof(tmp)
  } else if (model=="Mgamma") {
    data(data.inversegamma)
    attach(data.inversegamma)
    x <- log(data.inversegamma$Time)
    y <- data.inversegamma$Plasma
    fit <- glm(y ~ x, family = Gamma(link = "inverse"))
    Rprof(tmp <- tempfile())
    out_i <- independent("Mgamma",y,fit,B=1000)
    Rprof()
    detach(data.inversegamma)
    summaryRprof(tmp)
  } else if (model=="Mgaussian") {
    data(data.gaussian)
    y <- data.gaussian$taste
    x1 <- data.gaussian$h2s
    x2 <- data.gaussian$lactic
    n <- length(y)
    fit <- lm(y~(x1+x2))
    target <- "level"
    targetvalue <- c(0.5,0.9)
    Rprof(tmp <- tempfile())
    out_b <- boundary("Mgaussian",y,fit,target,targetvalue,B=1000)
    Rprof()
    summaryRprof(tmp)
  } else {
    stop(paste("model",model,"was not recognized"))
  }
}
