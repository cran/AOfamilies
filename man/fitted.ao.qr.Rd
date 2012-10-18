\name{fitted.ao.qr}
\alias{fitted.ao.qr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracting fitted values of an "ao.qr" object 
}
\description{
Function to extract the fitted values of a QR model fitted using \code{\link{ao.qr}}
}
\usage{
\method{fitted}{ao.qr} (object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "ao.qr"
}
  \item{...}{
Not used
}
}

\value{
A list containing a dataframe with 3 columns: the fitted values on the untransformed scale, the lower bound and the higher bound of the estimates
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## labor dataset
data(labor)
data <- labor

## fit the model using ao.qr() with a symmetric transformation  
## for the quantile 0.25 with the method "linprog"
fit.ao.qr <- ao.qr(pain ~ treatment*time, data = data, 
                   kappa = 0.25, estimation = "linprog",
                   transfo = "ao.sym") 
					
## extract the fitted values
fitted.values <- fitted.ao.qr (fit.ao.qr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
