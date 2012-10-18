\name{print.predict.ao.glm}
\alias{print.predict.ao.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing predicted values based on an "ao.glm" object
}
\description{
Print method of the function \code{\link{predict.ao.glm}}
}
\usage{
\method{print}{predict.ao.glm} (x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "ao.glm"
}
  \item{...}{
Not used
}
}

\examples{
## beetle dataset  
dose <- c (1.6907,1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)

totals <- c (59,60,62,56,63,59,62,60)

dead <- c (6,13,18,28,52,53,61,60)

prop_dead <- dead / totals

## fit the model using ao.glm() with an asymmetric transformation 
fit.ao <- ao.glm (prop_dead ~ dose, weights = totals, link = ao.asym)

## perform the prediction for a new data point and print the result 
predicted.values <- predict.ao.glm (fit.ao, newdata = data.frame (dose = c(1.7)))
predicted.values
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
