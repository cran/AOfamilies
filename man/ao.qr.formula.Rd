\name{ao.qr.formula}
\alias{ao.qr.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Formula method
}
\description{
Formula method of ao.qr()}
\usage{
\method{ao.qr}{formula}(formula, data = list(), 
weights = rep(1, length(y)), kappa = 0.5, 
phi = seq(0, 1.5, 0.005), estimation = "laplace", 
epsilon = 0.001, transfo = "ao.sym", 
plotit = "TRUE", method = "br", se = "boot", R = 100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A symbolic description of the model to fit. See Details section for further information
}
  \item{data}{
An optional data frame containing the variables in the model
}
  \item{weights}{
An optional vector of weights to be used in the fitting process
} 
  \item{kappa}{
Quantile to be estimated. Default value is 0.5
}
\item{phi}{
Transformation parameter. The default value is a sequence from 0 to 1.5 by step of 0.005
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
\item{plotit}{
If "TRUE" (default), the plot of the profile log-likelihood is produced when calling the summary method, otherwise not
}
\item{method}{The algorithmic method used to compute the fit if linprog is chosen as estimation process (see below).
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
\details{
Regarding the argument formula: a typical predictor has the form \emph{response ~ terms} where \emph{response} is the (numeric) response vector and \emph{terms} is a series of terms which specifies a linear predictor for \emph{response}. 
}
\value{
An object of class "ao.qr" is a list containing at least the following components:
  \item{MLE}{The MLE of the transformation parameter}
  \item{logLik.MLE}{The log-likelihood of the model using the MLE of the transformation parameter}
  \item{ao}{A vector of AO-transformed values based on the initial bounded response variable using the MLE of the transformation parameter} 
  \item{fit.coef.table}{A table containing the coefficients of the QR model, the standard errors, the t value and the test statistic. If the method "laplace" is used, the confidence bounds of the coefficients are included}
  \item{fit.coef}{The vector of coefficients of the QR model}
  \item{predict.all}{A matrix containing the fitted values and the confidence intervals on the transformed scale}
  \item{fit}{A list containing all the outputs from the fitted QR model. The components are the ones of the functions rq() or lqm() in packages \pkg{quantreg} and \pkg{lqmm} respectively}
}

\note{
To obtain special case transformations, use:
a/for the logit: transfo = "ao.sym" and phi = 0,
b/for the cloglog: transfo = "ao.asym" and phi = 0,
c/for the probit: transfo = "ao.sym" and phi = 0.3955, 
d/for the arcsine: transfo = "ao.sym" and phi = 0.6755

Available methods for an object of class "ao.qr": \code{\link{summary.ao.qr}}, \code{\link{fitted.ao.qr}} and \code{\link{predict.ao.qr}} 
}

\author{
Hakim-Moulay Dehbi}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
