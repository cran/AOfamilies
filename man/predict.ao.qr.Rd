\name{predict.ao.qr}
\alias{predict.ao.qr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting values based on an "ao.qr" object
}
\description{
Function to predict values based on a QR model fitted using \code{\link{ao.qr}}
}
\usage{
\method{predict}{ao.qr} (object, newdata = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "ao.qr"
}
  \item{newdata}{
An optional data frame in which to look for variables with which to predict. If
omitted, the fitted values are used
}
  \item{...}{
Not used
}
}

\value{
A list containing a dataframe with 3 columns: the point estimates on the untransformed scale, the lower bound and the higher bound of the estimates
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## labor dataset
data(labor)
data <- labor
data.subset <- subset (data, time == 90)

## fit the model using ao.qr() with a symmetric transformation  
## for the quantile 0.25 with the method "linprog"
fit.ao.qr <- ao.qr(pain ~ treatment, data = data.subset, 
                   kappa = 0.25, estimation = "linprog",
                   transfo = "ao.sym") 
					
## perform the prediction for a new data point 
predicted.values <- predict.ao.qr (fit.ao.qr, newdata = data.frame (treatment = 0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
