\name{print.summary.ao.qr}
\alias{print.summary.ao.qr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing the summary obtained from an "ao.qr" object
}
\description{
Print method of the function \code{\link{summary.ao.qr}}
}
\usage{
\method{print}{summary.ao.qr} (x, ...)
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
This function prints the call of the function, the quantile considered, the MLE and its associated log-likelihood. It formats the coefficients, standard errors, etc. and additionally gives "significance stars"
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

## print the summary of the model
summary (fit.ao.qr)	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
