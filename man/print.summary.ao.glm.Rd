\name{print.summary.ao.glm}
\alias{print.summary.ao.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing the summary obtained from an "ao.glm" object 
}
\description{
Print method of the function \code{\link{summary.ao.glm}}
}
\usage{
\method{print}{summary.ao.glm} (x, ...)
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
This function prints the call of the function, the result of the asymmetry test as defined in Aranda-Ordaz (1981), the MLE and its associated log-likelihood. It also gives the 95 percent confidence interval for the MLE and the special case transformations included in it. It formats the coefficients, standard errors, etc. and additionally gives "significance stars" 
}

\examples{
## beetle dataset  
dose <- c (1.6907,1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)

totals <- c (59,60,62,56,63,59,62,60)

dead <- c (6,13,18,28,52,53,61,60)

prop_dead <- dead / totals

## print the summary of the model 
summary (ao.glm (prop_dead ~ dose, weights = totals, link = ao.asym))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
