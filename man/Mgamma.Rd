\name{Mgamma}
\alias{Mgamma}

\title{Gamma Regression Model}

\description{
This function is a private function that returns the basic statistics of a selected model. It is only used in conjunction with \code{\link{boundary}} or \code{\link{independent}} sampling method.
}

\usage{
Mgamma(y, fit, cov)
}

\arguments{
  \item{y}{response variables.}
  \item{fit}{basic statistics after fitting a linear model by class \code{lm}.}
  \item{cov}{a covariance matrix of the parameters. System will use default covariance matrix if it is not specified.}
}

\details{
Both "log" link and "inverse" link are available in the system. The link information is stored in argument "fit" when the "lm" class is created.
}

\examples{
\dontrun{
#############################################################
################-----IndependentSampling-----################
#############################################################
#########Log Link#########
data(data.loggamma)
x <- log(data.loggamma$Brain)
y <- data.loggamma$Body
fit <- glm(y ~ x, family = Gamma(link = "log"))
#########Inverse Link#########
data(data.inversegamma)
attach(data.inversegamma)
x <- log(data.inversegamma$Time)
y <- data.inversegamma$Plasma
fit <- glm(y ~ x, family = Gamma(link = "inverse"))
########################################################
out_i <- independent("Mgamma",y,fit,B=1000)
########################################################
out_i$diag    # out_i$diag is equivalent to out_i$diagnosis
out_i$ind[100:140,]    # out_i$ind is equivalent to out_i$independent.sample
out_i$num    #out_i$num is equivalent to out_i$numWald.interval
out_i$sim    #out_i$sim is equivalent to out_i$simWald.interval
########################################################
par(mfrow=c(2,2),pty = "s")
plot(out_i$ind[,5],out_i$ind[,6],xlab=expression(beta[0]),ylab=expression(beta[1]),cex=0.5)
points(out_i$EST[1],out_i$EST[2],pch=16,col="red",cex=1.5)
plot(out_i$ind[,5],out_i$ind[,7],xlab=expression(beta[0]),ylab=expression(nu),cex=0.5)
points(out_i$EST[1],out_i$EST[3],pch=16,col="red",cex=1.5)
plot(out_i$ind[,6],out_i$ind[,7],xlab=expression(beta[1]),ylab=expression(nu),cex=0.5)
points(out_i$EST[2],out_i$EST[3],pch=16,col="red",cex=1.5)

##########################################################
################-----BoundarySampling-----################
##########################################################
#########Log Link#########
data(data.loggamma)
x <- log(data.loggamma$Brain)
y <- data.loggamma$Body
fit <- glm(y ~ x, family = Gamma(link = "log"))
target <- "level"
targetvalue <- c(0.5,0.9)
#########Inverse Link#########
data(data.inversegamma)
x <- log(data.inversegamma$Time)
y <- data.inversegamma$Plasma
fit <- glm(y ~ x, family = Gamma(link = "inverse"))
target <- "level"
targetvalue <- c(0.5,0.9)
########################################################
out_b <- boundary("Mgamma",y,fit,target,targetvalue,B=1000)
########################################################
out_b$diag    # out_b$diag is equivalent to out_b$diagnosis
out_b$bound[1:20,]    # out_b$bound is equivalent to out_b$boundary.sample
out_b$num    # out_b$num is equivalent to out_b$numWald.interval
out_b$sim    # out_b$sim is equivalent to out_b$simWald.interval
out_b$convnum   # out_b$convnum is equivalent to out_b$convnumWald
out_b$convsim   # out_b$convsim is equivalent to out_b$convsimWald
########################################################
par(mfrow=c(2,2))
plot(out_b$bound[,5],out_b$bound[,6],xlab=expression(beta[0]),ylab=expression(beta[1]),cex=0.5)
points(out_b$EST[1],out_b$EST[2],pch=16,col="red",cex=1.5)
plot(out_b$bound[,5],out_b$bound[,7],xlab=expression(beta[0]),ylab=expression(nu),cex=0.5)
points(out_b$EST[1],out_b$EST[3],pch=16,col="red",cex=1.5)
plot(out_b$bound[,6],out_b$bound[,7],xlab=expression(beta[1]),ylab=expression(nu),cex=0.5)
points(out_b$EST[2],out_b$EST[3],pch=16,col="red",cex=1.5)
}}