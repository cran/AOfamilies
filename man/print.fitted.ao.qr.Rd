\name{print.fitted.ao.qr}
\alias{print.fitted.ao.qr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing the fitted values of a model fitted using \code{\link{ao.qr}} 
}
\description{
Print method of the function \code{\link{fitted.ao.qr}}
}
\usage{
\method{print}{fitted.ao.qr} (x, ...)
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

\examples{
## labor dataset
data(labor)
data <- labor

## fit the model using ao.qr() with a symmetric transformation  
## for the quantile 0.25 with the method "linprog"
fit.ao.qr <- ao.qr(pain ~ treatment*time, data = data, 
                   kappa = 0.25, estimation = "linprog",
                   transfo = "ao.sym")
					
## extract and print the fitted values
fitted.values <- fitted.ao.qr (fit.ao.qr)
fitted.values
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
