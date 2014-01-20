\name{ao.glm.formula}
\alias{ao.glm.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Formula method
}
\description{
Formula method of ao.glm()
}
\usage{
\method{ao.glm}{formula}(formula, data = list(), link, 
phi = seq(-2, 2, 0.01), weights, plotit = "TRUE", 
plot.spline = "TRUE", ...)
%ao.glm.formula(formula, data = list(), link, phi = seq(-2, 2, 0.01), weights, plotit = "TRUE", plot.spline = "TRUE", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A symbolic description of the model to fit. See Details section for further information
}
  \item{data}{
An optional data frame containing the variables in the model
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
  \item{plotit}{
If "TRUE" (default), the plot of the profile log-likelihood is produced when calling the summary method, otherwise not
}
  \item{plot.spline}{
If "TRUE" (default), it fits a cubic smoothing spline to the profile log-likelihood data
}
  \item{\dots}{
Further arguments passed to or from other methods
}
}
\details{
Regarding the argument formula: a typical predictor has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which speci?es a linear predictor for response. The response should consist of proportions of successes for each binomial trial of the dataset
}
\value{
An object of class "ao.glm" is a list containing at least the following components:
  \item{MLE}{The MLE of the transformation parameter}
  \item{fit.MLE}{A list containing all the outputs from the fitted model. The components of the list are the ones of the function glm() in package \pkg{base}}
  \item{fit.MLE.coef}{The vector of estimated coefficients}
  \item{logLik}{The log-likelihood of the model using the MLE of the transformation parameter}
  \item{valid.lambda}{A vector containing the values of the transformation parameters for which convergence was reached by the IWLS algorithm}
  \item{valid.logLik}{A vector containing the log-likelihood associated with the valid values of the transformation parameter} 
  \item{fitted.values}{The fitted mean values, obtained by transforming the linear predictors by the inverse of the link function}
}

\author{
Hakim-Moulay Dehbi
}
\note{
To obtain special case transformations, use:
a/for the logit: transfo = "ao.sym" and phi = 0,
b/for the cloglog: transfo = "ao.asym" and phi = 0,
c/for the probit: transfo = "ao.sym" and phi = 0.3955, 
d/for the arcsine: transfo = "ao.sym" and phi = 0.6755

Available methods for an object of class "ao.glm": \code{\link{summary.ao.glm}}, \code{\link{fitted.ao.glm}}, \code{\link{predict.ao.glm}}, \code{\link{plot.ao.glm}} and \code{\link{ao.dose.p}} 
}

\references{
Aranda-Ordaz, F. J. (1981). On two families of transformations to additivity for binary response data. Biometrika 68 (2), pp. 357-363
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{summary.ao.glm}}, \code{\link{fitted.ao.glm}}, \code{\link{predict.ao.glm}}, \code{\link{plot.ao.glm}} and \code{\link{ao.dose.p}} 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
