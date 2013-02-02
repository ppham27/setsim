\name{data.gaussian}
\alias{data.gaussian}
\docType{data}
\title{Mature Cheddar Cheese}

\description{
This data set gives measurements of acetic acid, hydrogen sulfide and lactic acid 
and an index of taste for a sample of mature cheddar cheese.
}

\usage{data(data.gaussian)}
\format{
  A data frame with 30 observations on the following 5 variables.
  Note that the first line contains the names of the variables.
  \describe{
    \item{\code{case}}{a numeric vector}
    \item{\code{taste}}{a numeric vector}
    \item{\code{acetic}}{a numeric vector}
    \item{\code{h2s}}{a numeric vector}
    \item{\code{lactic}}{a numeric vector}
  }
}

\details{
The goal is to investigate the relationships between the taste of matured
cheese (response variable, taste) and three chemicals deciding the taste of cheese:
acetic acid (acetic), hydrogen sulde (h2s), and lactic acid (lactic). 
Note that the first two predictors are log transformed.
}

\source{
Moore, D. S. and McCabe, G. P. (1989) \emph{Introduction to the Practice of Statistics}.
New York: W.H. Freeman. 
}

\examples{
data(data.gaussian)
## maybe str(data.gaussian) ; plot(data.gaussian) ...
}

\keyword{datasets}