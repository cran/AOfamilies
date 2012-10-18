\name{ao.glm.fit}
\alias{ao.glm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Defining the fitting method of \code{\link{ao.glm}}}
\description{
This function corresponds to the fitting method of \code{\link{ao.glm}}
}
\usage{
\method{ao.glm}{fit} (x, y, link, phi, weights, maxit = 500, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Design matrix
}
  \item{y}{
Response variable
}
  \item{link}{
2 possibilities: ao.sym or ao.asym (ao.sym for the symmetric transformation and ao.asym for the asymmetric one)}
  \item{phi}{
An optional vector specifying the value(s) of the transformation parameter to use 
}
  \item{weights}{
An optional vector of weights to be used in the fitting process
}
  \item{maxit}{
Maximum number of iterations used in the IWLS algorithm
}
  \item{...}{
Further arguments passed to or from other methods
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
