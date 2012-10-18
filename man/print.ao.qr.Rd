\name{print.ao.qr}
\alias{print.ao.qr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing the MLE, its associated log-likelihood and the coefficients of an "ao.qr" object}
\description{
Print method of the function \code{\link{ao.qr}}
}
\usage{
\method{print}{ao.qr} (x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "ao.qr"
}
  \item{...}{
Not used
}
}
\details{
Printing the MLE of the transformation parameter, its associated log-likelihood and the coefficients of the regression
}

\examples{
## labor dataset
data(labor)
data <- labor

## fit the model using ao.qr() with a symmetric transformation 
## for the quantile 0.25 with the method "linprog"
fit.ao.qr <- ao.qr(pain ~ treatment*time, data = data, 
                   kappa = 0.25, estimation = "linprog",
                   transfo = "ao.sym") 
fit.ao.qr
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
