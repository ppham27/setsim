\name{data.loglinear}
\alias{data.loglinear}
\docType{data}
\title{Florida Death Penalty}

\description{
This data is concerned with 674 defendants in indictments involving cases
with multiple murders in Florida between 1976 and 1987.
The variables are death penalty verdict (Pen), having the categories, (yes and no),
race of defendant (Def), having the categories, (white and black) and 
race of victims (Vic), having the categories, (white and black). 
}

\usage{data(data.loglinear)}
\format{
  A data frame with 8 observations on the following 4 variables.
  Note that the first line contains the names of the variables.
  \describe{
    \item{\code{count}}{a numeric vector}
    \item{\code{Vic}}{a factor with levels \code{b} \code{w}}
    \item{\code{Def}}{a factor with levels \code{b} \code{w}}
    \item{\code{Pen}}{a factor with levels \code{n} \code{y}}
  }
}

\details{
The goal is to study the effects of racial characteristics on whether individuals 
convicted of homicide receive the death penalty.
}

\source{
Agresti, A. (2002). \emph{Categorical Data Analysis} (2nd edn). New York: John Wiley and Sons.
}

\examples{
data(data.loglinear)
## maybe str(data.loglinear) ; plot(data.loglinear) ...
}

\keyword{datasets}
