\name{plot.ao.glm}
\alias{plot.ao.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ploting diagnostics for an "ao.glm" object
}
\description{
This function corresponds to the function plot.lm() of package \pkg{stats}. Six plots (selectable by which) are currently available: a plot of residuals against fitted values, a Scale-Location plot of \emph{sqrt}(|residuals|) against fitted values, a Normal Q-Q plot, a plot of Cook's distances versus row labels, a plot of residuals against leverages, and a plot of Cook's distances against leverage/(1-leverage).
}
\usage{
\method{plot}{ao.glm} (x, which = 1:4, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "ao.glm"
}
  \item{which}{
The first four plots are provided by default. If a subset of the plots is required, specify a subset of the numbers 1:6
}
  \item{...}{
Further arguments passed to or from other methods
}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ao.glm}}
}
\examples{
## beetle dataset  
dose <- c (1.6907,1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)

totals <- c (59,60,62,56,63,59,62,60)

dead <- c (6,13,18,28,52,53,61,60)

prop_dead <- dead / totals

## fit the model using ao.glm() with an asymmetric transformation 
fit.ao <- ao.glm (prop_dead ~ dose, weights = totals, link = ao.asym)

## plot the first two diagnostics
plot.ao.glm (fit.ao, which = 1:2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
