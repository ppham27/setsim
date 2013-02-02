\name{data.logistic}
\alias{data.logistic}
\docType{data}
\title{Disease Outbreak}

\description{
This data is concerned with 196 individuals who were interviewed to assess 
if certain specific symptoms associated with a disease spread by a mosquito were 
present. The first 98 cases were included in the current data set. 
The variables are disease status (disease), having the categories, (yes=1, no=0), 
age (age), socioeconomic status of household with three levels, (upper, middle, lower), 
and a city sector, having the categories, (sector 1=0 and section 2=1).
Note that a socioeconomic status variable is represented by two indicator variables (col3 and col4) :
(col3,col4) = (0,0) for upper level, (1,0) for middle level and (0,1) for lower level. 
}

\usage{data(data.logistic)}
\format{
  A data frame with 98 observations on the following 6 variables.
  Note that the first line contains the names of the variables.
  \describe{
    \item{\code{obsnum}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{col3}}{a numeric vector}
    \item{\code{col4}}{a numeric vector}
    \item{\code{sector}}{a numeric vector}
    \item{\code{disease}}{a numeric vector}
  }
}

\details{
The goal is to model the odds of disease using 
age, two numerical indicators of socioeconomic status and a numerical indicator of sector.
}

\source{
Kutner, M. H., Nachtsheim, C. J., Neter, J. and Li, W. (2005) 
\emph{Applied Linear Statistical Models} (5th edition). Chicago: Irwin.
}

\examples{
data(data.logistic)
## maybe str(data.logistic) ; plot(data.logistic) ...
}

\keyword{datasets}
