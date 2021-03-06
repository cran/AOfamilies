\name{ao.glm.fit}
\alias{ao.glm.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Central fitting function
}
\description{
Central fitting function of ao.glm()
}
\usage{
\method{ao.glm}{fit}(x, y, link, phi, weights, maxit = 500, ...)
%ao.glm.fit(x, y, link, phi, weights, maxit = 500, ...)
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
Link function: ao.sym or ao.asym
}
 \item{phi}{
An optional vector specifying the value(s) of the transformation parameter to
}
  \item{weights}{
An optional vector of weights to be used in the fitting process
}
  \item{maxit}{
Maximum number of iterations
}
  \item{\dots}{
Further arguments passed to or from other methods
}
}
\author{Hakim-Moulay Dehbi
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
