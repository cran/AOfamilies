\name{ao.qr.default}
\alias{ao.qr.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default method}
\description{
Default method of ao.qr
}
\usage{
\method{ao.qr}{default}(x, y, weights, method, kappa, phi, estimation, epsilon, transfo, se, R, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix}
  \item{y}{Response variable}
  \item{weights}{
An optional vector of weights to be used in the fitting process
} 
 \item{method}{The algorithmic method used to compute the fit if linprog is chosen as estimation process (see below).
}
  \item{kappa}{
Quantile to be estimated. Default value is 0.5
}
  \item{phi}{
Transformation parameter. The default value is a sequence from 0 to 1.6 by step of 0.05
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
  \item{se}{
Method used to compute standard error if estimation is performed with "linprog". Options are "boot" (default), "iid", "nid" and "ker". See package \pkg{quantreg}
}
\item{R}{
Number of bootstrap replications
}
  \item{\dots}{
Further arguments passed to or from other methods
}
}
\author{
Hakim-Moulay Dehbi}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
