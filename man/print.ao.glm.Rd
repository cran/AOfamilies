\name{print.ao.glm}
\alias{print.ao.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing the MLE, its associated log-likelihood and the coefficients of an "ao.glm" object
}
\description{
Print method of the function \code{\link{ao.glm}}
}
\usage{
\method{print}{ao.glm} (x, ...)
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
\details{
Printing the MLE of the transformation parameter, its associated log-likelihood and the coefficients of the regression
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## beetle dataset  
dose <- c (1.6907,1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)

totals <- c (59,60,62,56,63,59,62,60)

dead <- c (6,13,18,28,52,53,61,60)

prop_dead <- dead / totals

## fit the model using ao.glm() with an asymmetric transformation and print
## the results
fit.ao <- ao.glm (prop_dead ~ dose, weights = totals, link = ao.asym)
fit.ao
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
