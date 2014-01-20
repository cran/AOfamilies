\name{ao.dose.p}
\alias{ao.dose.p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating lethal dose(s) for binomial assay model fitted with an AO transformation as link function
}
\description{
Function to calibrate binomial assays, generalizing the calculation of LD50
}
\usage{
ao.dose.p (object, cf = 1:2, p = 0.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "ao.glm"
}
  \item{cf}{
The terms in the coefficient vector giving the intercept and coefficient of (log-) dose
}
  \item{p}{
Probabilities at which to calculate the lethal dose
}
  \item{\dots}{
Not used
}
}
\value{
A list containing the predictions at each probability level, the associated standard errors and the confidence intervals
}
\references{
Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Springer.
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

## estimate the lethal dose at p=0.25, p=0.5 and p=0.75
lethal.doses <- ao.dose.p (fit.ao, p = c(0.25,0.5,0.75))

}
\seealso{
Function dose.p() of package \pkg{MASS}

\code{\link{ao.glm}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
