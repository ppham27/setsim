\name{data.negbino}
\alias{data.negbino}
\docType{data}
\title{Attendance Behavior}

\description{
This data is an attendance data on 316 high school juniors from two urban high schools.
The response variable is days absent (daysabs). 
The predictors are the math standardized test scores (math), 
the language arts standardized test scores (langarts), and 
a binary indicator of student gender (male).
}

\usage{data(data.negbino)}
\format{
  A data frame with 317 observations on the following 7 variables.
  Note that the first line contains the names of the variables.
  \describe{
    \item{\code{id}}{a numeric vector}
    \item{\code{school}}{a numeric vector}
    \item{\code{male}}{a numeric vector}
    \item{\code{math}}{a numeric vector}
    \item{\code{langarts}}{a numeric vector}
    \item{\code{daysatt}}{a numeric vector}
    \item{\code{daysabs}}{a numeric vector}
  }
}

\details{
The goal is to study the attendance behavior of high school juniors at two schools. 
Note that the variable names, id and school, represent 316 high schools 
and two urban high schools (1 and 2), respectively. 
}

\source{
http://www.ats.ucla.edu/stat/sas/dae/negbinreg.htm
}

\examples{
data(data.negbino)
## maybe str(data.negbino) ; plot(data.negbino) ...
}

\keyword{datasets}
