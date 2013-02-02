\name{bisect}
\alias{bisect}

\title{Bisection Method}

\description{
This function finds the approximate roots of a numeric function by bisection algorithm. The root-finding method is one-dimensional.
}

\usage{
bisect(fn, a, b, tol = 1e-05, error = 1e-04, decimal = 5, input, VZ, c_value)
}

\arguments{
  \item{fn}{a numeric function defined by the user. An object which is not a numeric function will not be allowed.}
  \item{a}{numeric value of the lower end point of the interval to be searched.}
  \item{b}{numeric value of the upper end point of the interval to be searched.}
  \item{tol}{tolerance to stop bisection algorithm. The default value is 0.00001.}
  \item{error}{significance of the approximate root. The default value is 0.0001.}
  \item{decimal}{the number of decimal places to round off for the approximate root.}
  \item{input}{a list of response variable, design(X) matrix, maximum likelihood estimators for the parameters and maximized log likelihood value. This argument must be ignored if confidence region visualization algorithm is not used.}
  \item{VZ}{a matrix product by multiplying the matrix square root of a covariance matrix for the maximum likelihood estimator and random vector generated from a standard multivariate normal distribution. This argument must be ignored if confidence region visualization algorithm is not used.}
  \item{c_value}{the critical values to visualize confidence regions. This argument must be ignored if confidence region visualization algorithm is not used.}
}

\value{the approximate root by bisection algorithm. The error of approximation, number of decimal places and other requirements can be adjusted by user to produce a different root.}

\details{
The function requires at least one root to process. Check \code{\link{nroots}} if a roots exists.
}

\seealso{
\code{\link{nroots}} for finding the lower and upper end point of the interval containing a root.}

\section{Warning}{If this function is not used for confidence region visualization algorithm, the numeric function to be calculated must include \{\code{\ldots}\} in its argument. This is due to the last three arguments of this function.}

\examples{
y <- function(x,...){x^2-4}
bisect(y,a=-5,b=0,error=0.1,decimal=5)
bisect(y,a=-5,b=0,error=0.0001,decimal=5)

f <- function(x,...){sin(x)}
x <- seq(-50,50,length=100)
roots <- nroots(f,x)
n <- roots$n
range <- roots$range
sol <- 0
for (i in 1:n){
	sol[i] <- bisect(f,a=range[i,1],b=range[i,2])}
sol
}

\keyword{bisect}