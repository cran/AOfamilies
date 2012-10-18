\name{print.predict.ao.qr}
\alias{print.predict.ao.qr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing predicted values based on an "ao.qr" object
}
\description{
Print method of the function \code{\link{predict.ao.qr}}
}
\usage{
\method{print}{predict.ao.qr} (x, ...)
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
data.subset <- subset (data, time == 90)

## fit the model using ao.qr() with a symmetric transformation (transformation parameter = 0.2) 
## for the quantile 0.25 with the method "linprog"
fit.ao.qr <- ao.qr(pain ~ treatment, data = data.subset, 
                   kappa = 0.25, estimation = "linprog",
                   transfo = "ao.sym") 
					
## perform the prediction for a new data point and print the result
predicted.values <- predict.ao.qr (fit.ao.qr, newdata = data.frame (treatment = 0))
predicted.values
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
