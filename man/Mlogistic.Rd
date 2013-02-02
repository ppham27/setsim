\name{Mlogistic}
\alias{Mlogistic}

\title{Logistic Regression Model}

\description{
This function is a private function that returns the basic statistics of a selected model. It is only used in conjunction with \code{\link{boundary}} or \code{\link{independent}} sampling method.
}

\usage{
Mlogistic(y, fit, cov)
}

\arguments{
  \item{y}{response variables.}
  \item{fit}{basic statistics after fitting a linear model by class \code{lm}.}
  \item{cov}{a covariance matrix of the parameters. System will use default covariance matrix if it is not specified.}
}

\examples{
\dontrun{
#############################################################
################-----IndependentSampling-----################
#############################################################
library(MASS)
data(data.logistic)
y <- data.logistic$disease
fit <- glm(y ~ (data.logistic$age+data.logistic$sector), family=binomial)
########################################################
out_i <- independent("Mlogistic",y,fit,B=1000)
########################################################
out_i$diag    # out_i$diag is equivalent to out_i$diagnosis
out_i$ind[100:140,]    # out_i$ind is equivalent to out_i$independent.sample
out_i$num    #out_i$num is equivalent to out_i$numWald.interval
out_i$sim    #out_i$sim is equivalent to out_i$simWald.interval
########################################################
par(mfrow=c(2,2))
plot(out_i$ind[,7]~out_i$ind[,6],xlab=expression(beta[1]),ylab=expression(beta[2]),cex=0.5)
points(out_i$MLE[2],out_i$MLE[3],pch=16,col="red",cex=1.5)
plot(out_i$ind[,6]~out_i$ind[,5],xlab=expression(beta[0]),ylab=expression(beta[1]),cex=0.5)
points(out_i$MLE[1],out_i$MLE[2],pch=16,col="red",cex=1.5)
plot(out_i$ind[,7]~out_i$ind[,5],xlab=expression(beta[0]),ylab=expression(beta[2]),cex=0.5)
points(out_i$MLE[1],out_i$MLE[3],pch=16,col="red",cex=1.5)

##########################################################
################-----BoundarySampling-----################
##########################################################
library(MASS)
data(data.logistic)
y <- data.logistic$disease
fit <- glm(y ~ (data.logistic$age+data.logistic$sector), family=binomial)
target <- "level"
targetvalue <- c(0.5,0.9)
########################################################
out_b <- boundary("Mlogistic",y,fit,target,targetvalue,B=1000)
########################################################
out_b$diag    # out_b$diag is equivalent to out_b$diagnosis
out_b$bound[1:20,]    # out_b$bound is equivalent to out_b$boundary.sample
out_b$num    # out_b$num is equivalent to out_b$numWald.interval
out_b$sim    # out_b$sim is equivalent to out_b$simWald.interval
out_b$convnum   # out_b$convnum is equivalent to out_b$convnumWald
out_b$convsim   # out_b$convsim is equivalent to out_b$convsimWald
########################################################
par(mfrow=c(2,2))
plot(out_b$bound[,7]~out_b$bound[,6],xlab=expression(beta[1]),ylab=expression(beta[2]),cex=0.5)
points(out_b$MLE[2],out_b$MLE[3],pch=16,col="red",cex=1.5)
plot(out_b$bound[,6]~out_b$bound[,5],xlab=expression(beta[0]),ylab=expression(beta[1]),cex=0.5)
points(out_b$MLE[1],out_b$MLE[2],pch=16,col="red",cex=1.5)
plot(out_b$bound[,7]~out_b$bound[,5],xlab=expression(beta[0]),ylab=expression(beta[2]),cex=0.5)
points(out_b$MLE[1],out_b$MLE[3],pch=16,col="red",cex=1.5)
}}