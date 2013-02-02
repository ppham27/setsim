\name{data.poisson}
\alias{data.poisson}
\docType{data}
\title{Smart Phone Users}

\description{
This is a hypothetical data that gives measurements of 
a type of smart phone (device A and B), unlimited data plan (yes, no) and 
a number of emails sent for a sample of 20 smart phone users.
}

\usage{data(data.poisson)}
\format{
  A data frame with 20 observations on the following 3 variables.
  Note that the first line contains the names of the variables.
  \describe{
    \item{\code{Device}}{a factor with levels \code{A} \code{B}}
    \item{\code{DataPlan}}{a factor with levels \code{N} \code{Y}}
    \item{\code{y}}{a numeric vector}
  }
}

\details{
The goal is to study whether or not 
the number of emails sent from mobile devices 
can be influenced by the difference of the devices 
and an unlimited data plan.
}

\source{
http://yatani.jp/HCIstats/GLM
}

\examples{
data(data.poisson)
## maybe str(data.poisson) ; plot(data.poisson) ...
}

\keyword{datasets}
