testC <- function() {
  data(data.poisson)
  Device <- data.poisson$Device
  DataPlan <- data.poisson$DataPlan
  y <- data.poisson$y
  fit <- glm(y ~ (Device+DataPlan), family=poisson)
  target <- "level"
  targetvalue <- c(0.5,0.9)
  Model <- eval(call("Mpoisson", y, fit, NULL))

  .Call("runUnitTests",
        Model$H, environment(Model$H), Model$input)
}
