\name{data.inversegamma}
\alias{data.inversegamma}
\docType{data}
\title{Clotting Time}

\description{
Consider the gamma regression model with an inverse link. 
The first column is the response variable, percentage concentration plasma (Plasma) and 
the second column represents the predictor, clotting time (seconds) (Time). 
}

\usage{data(data.inversegamma)}

\format{
  A data frame with 9 observations on the following 2 variables.
  \describe{
    \item{\code{Plasma}}{a numeric vector}
    \item{\code{Time}}{a numeric vector}
  }
}

\source{
chap11.pdf in http://www.stat.ncsu.edu/people/davidian/courses/st732/notes/
}

\examples{
data(data.inversegamma)
## maybe str(data.inversegamma) ; plot(data.inversegamma) ...
}

\keyword{datasets}
