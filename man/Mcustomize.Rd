\name{Mcustomize}
\alias{Mcustomize}

\title{Customized Model}

\description{
This function is a private function that returns the basic statistics of a customized model. It is only used in conjunction with \code{\link{boundary}} or \code{\link{independent}} sampling method.
}

\usage{
Mcustomize(y, fit, cov)
}

\arguments{
  \item{y}{response variables.}
  \item{fit}{a list of maximum likelihood estiamtors for the parameters, design(X) matrix and expression of likelihood function.}
  \item{cov}{a covariance matrix of the parameters.}
}

\details{
\code{fit} must be a list. It can be done by \code{\link{list}} function. The names must exactly be ``MLE'', ``X'' and ``lik''. \code{lik} in \code{fit} must be the \code{\link{expression}} class.
}

\examples{
\dontrun{
##########################################################
################-----BoundarySampling-----################
##########################################################
library(MASS)
data(data.poisson)
Device <- data.poisson$Device
DataPlan <- data.poisson$DataPlan
y <- data.poisson$y
fit.glm <- glm(y ~ (Device+DataPlan), family=poisson)
X <- model.matrix(fit.glm)
MLE <- fit.glm$coef
cov <- vcov(fit.glm)
lik <- expression(MLE \%*\% t(X) \%*\% y - sum(exp(MLE \%*\% t(X)))-sum(lgamma(y+1)))
fit <- list(MLE=MLE,X=X,lik=lik)
target <- "level"
targetvalue <- c(0.5,0.9)
########################################################
out_b <- boundary("Mcustomize",y,fit,target,targetvalue,cov=cov)
########################################################
out_b$diag    # out_b$diag is equivalent to out_b$diagnosis
out_b$bound[1:20,]    # out_b$bound is equivalent to out_b$boundary.sample
out_b$num    # out_b$num is equivalent to out_b$numWald.interval
out_b$sim    # out_b$sim is equivalent to out_b$simWald.interval
out_b$convnum   # out_b$convnum is equivalent to out_b$convnumWald
out_b$convsim   # out_b$convsim is equivalent to out_b$convsimWald
########################################################
par(mfrow=c(2,2))
plot(out_b$bound[,6],out_b$bound[,7],xlab=expression(beta[Def]),ylab=expression(beta[Vic]),cex=0.5)
points(out_b$MLE[2],out_b$MLE[3],pch=16,col="red",cex=1.5)
plot(out_b$bound[,7],out_b$bound[,8],xlab=expression(beta[Vic]),ylab=expression(beta[Pen]),cex=0.5)
points(out_b$MLE[3],out_b$MLE[4],pch=16,col="red",cex=1.5)
plot(out_b$bound[,6],out_b$bound[,8],xlab=expression(beta[Def]),ylab=expression(beta[Def*Vic]),cex=0.5)
points(out_b$MLE[2],out_b$MLE[4],pch=16,col="red",cex=1.5)
plot(out_b$bound[,5],out_b$bound[,8],xlab=expression(beta[Def*Pen]),ylab=expression(beta[Vic*Pen]),cex=0.5)
points(out_b$MLE[1],out_b$MLE[4],pch=16,col="red",cex=1.5)
}}