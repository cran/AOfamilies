\name{print.ao.glm.dose}
\alias{print.ao.glm.dose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing the lethal dose(s) estimated with an "ao.glm" object
}
\description{
Print method of the function \code{\link{ao.dose.p}}
}
\usage{
\method{print}{ao.glm.dose} (x, ...)
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
Printing a dataframe giving the prediction at each probability level with standard errors and confidence intervals
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## beetle dataset  
dose <- c (1.6907,1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)

totals <- c (59,60,62,56,63,59,62,60)

dead <- c (6,13,18,28,52,53,61,60)

prop_dead <- dead / totals

## fit the model using ao.glm() with an asymmetric transformation 
fit.ao <- ao.glm (prop_dead ~ dose, weights = totals, link = ao.asym)

## estimate and print the lethal doses at p=0.25, p=0.5 and p=0.75
lethal.doses <- ao.dose.p (fit.ao, p = c(0.25,0.5,0.75))
lethal.doses
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
