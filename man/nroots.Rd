\name{nroots}
\alias{nroots}

\title{Finding the Number of Roots}

\description{
This function finds the number of roots for a given numeric function in a specified domain. It returns the number of roots, positions of the roots and two corresponding intervals containing each root and corresponding function value, respectively. A finer domain discretization will result in more accurate intervals for both a root and corresponding function value. If the user defined function does not intersect x axis, 0 root will be returned.
}

\usage{
nroots(fn, x)
}

\arguments{
  \item{fn}{a numeric function defined by the user. An object which is not a numeric function will not be allowed.}
  \item{x}{a numeric vector indicating the domain of the numeric function.}
}

\value{
  \item{n}{the number of roots in the domain.}
  \item{position}{the position of the root in the domain. The position will change if a different domain is specified.}
  \item{range}{two intervals containing the roots and corresponding values of the functions evaluated at those roots. \code{x.lower} and \code{x.upper} indicates the lower end point and upper end point of the interval on x axis that contains a root. \code{y.lower} and \code{y.upper} indicates the lower end point and upper end point of the interval on y axis that contains the value of the function evaluated at a root.}
}

\examples{
x <- seq(-3*pi,3*pi,length=50)
y <- function(x){sin(x)}
nroots(y,x)

a <- seq(10,50,length=30)
fun <- function(a){log(a)}
nroots(fun,a)
}

\keyword{nroots}