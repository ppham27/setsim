\name{data.loggamma}
\alias{data.loggamma}
\docType{data}
\title{Brain Weights}

\description{
This data concerns mean weights of the brains and bodies of 62 species of mammals. 
The first column contains values of species (Species), the second column is brain
weight in grams (Brain), and the third column is body weight in kilograms (Body). 
}

\usage{data(data.loggamma)}
\format{
  A data frame with 62 observations on the following 3 variables.
  \describe{
    \item{\code{Species}}{a string vector for species names}
    \item{\code{Brain}}{a numeric vector}
    \item{\code{Body}}{a numeric vector}
  }
}

\details{
The goal is to investigate the linear relationship between the natural logarithm of brain weight 
and that of body weight, assuming that brain weights
are observations from gamma distributions and their means depend on
logarithms of body weights.
}

\source{
Skovgaard, L. T. (2009). \emph{Confidence Intervals in Generalized Regression Models}. 
Boca Raton, FL : CRC Press. 
}

\examples{
data(data.loggamma)
## maybe str(data.loggamma) ; plot(data.loggamma) ...
}

\keyword{datasets}
