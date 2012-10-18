\name{ao.qr.fit}
\alias{ao.qr.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Defining the fitting method of \code{\link{ao.qr}}
}
\description{
This function corresponds to the fitting method of \code{\link{ao.qr}}
}
\usage{
\method{ao.qr}{fit} (x, y, weights, kappa, phi, se, estimation, epsilon, transfo, R, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Design matrix
}
  \item{y}{
Response variable
}
  \item{weights}{
An optional vector of weights to be used in the fitting process
} 
  \item{kappa}{
Quantile to be estimated. Default value is 0.5
}
  \item{phi}{
Transformation parameter. The default value is a sequence from 0 to 1.6 by step of 0.05
}
  \item{se}{
Method used to compute standard error if estimation is performed with "linprog". Options are "boot" (default), "iid", "nid" and "ker". See package \pkg{quantreg}
}
  \item{estimation}{
2 possibilities: "laplace" (default) or "linprog". Choosing "laplace" means fitting with the package \pkg{lqmm} while "linprog" implies \pkg{quantreg} 
}
  \item{epsilon}{
Tolerance parameter. The default value is 0.001
}
  \item{transfo}{
2 possibilities: "ao.sym" or "ao.asym" ("ao.sym" for the symmetric transformation and "ao.asym" for the asymmetric one)
}
\item{R}{
Number of bootstrap replications
}
  \item{...}{
Further arguments passed to or from other methods
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
