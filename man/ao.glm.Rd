\name{ao.glm}
\alias{ao.glm}
\alias{ao.glm.default}
\alias{ao.glm.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting GLMs with an AO transformation as link function
}
\description{
This function implements the AO symmetric and asymmetric transformations for binary response variables in GLMs
}
\usage{
ao.glm (x, ...)

\method{ao.glm}{default} (x, y, link, phi, weights, ...)

\method{ao.glm}{formula} (formula, data = list(), link, phi = seq(-2, 2, 0.01),
        weights, plotit = "TRUE", plot.spline = "TRUE", ...)

%ao.glm (formula, data, link, phi = seq (-2, 2, 0.01), weights, 
%        plotit = "TRUE", plot.spline = "TRUE", maxit = 500, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Design matrix
}
  \item{y}{
Response variable
}
  \item{formula}{
A symbolic description of the model to fit. See Details section for further information
}
  \item{data}{
An optional data frame containing the variables in the model 
}
  \item{link}{
2 possibilities: ao.sym or ao.asym (ao.sym for the symmetric transformation and ao.asym for the asymmetric one)
}
  \item{phi}{
An optional vector specifying the value(s) of the transformation parameter to use 
}
  \item{weights}{
An optional vector of weights to be used in the fitting process
}
  \item{plotit}{
If "TRUE" (default), the plot of the profile log-likelihood is produced when calling the summary method, otherwise not.
}
 \item{plot.spline}{
If "TRUE" (default), it fits a cubic smoothing spline to the profile log-likelihood data
}
%  \item{maxit}{
%Maximum number of iterations used in the IWLS algorithm
%}
  \item{...}{
Further arguments passed to or from other methods
}
}
\details{
Regarding the argument formula: a typical predictor has the form \emph{response ~ terms} where \emph{response} is the (numeric) response vector and \emph{terms} is a series of terms which specifies a linear predictor for \emph{response}. The \emph{response} should consist of proportions of successes for each binomial trial of the dataset 
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
\examples{
## Example 1: beetle dataset  
dose <- c(1.6907,1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839)

totals <- c(59,60,62,56,63,59,62,60)

dead <- c(6,13,18,28,52,53,61,60)

prop_dead <- dead / totals

# fit the model using ao.glm() with an asymmetric transformation
fit.ao <- ao.glm (prop_dead ~ dose, weights = totals, link = ao.asym)

## Example 2: skin cancer dataset 
age <- c(19.5, 29.5, 39.5, 49.5, 59.5, 69.5, 
         79.5, 89.5, 19.5, 29.5, 39.5, 49.5, 
         59.5, 69.5, 89.5)		   

totals <- c(172675, 123065, 96216,
            92051, 72159, 54722, 32185, 8328,
            181343, 146207, 121374, 111353, 
            83004, 55932, 7583)

cases <- c(1, 16, 30, 71,
           102, 130, 133, 40,
		   4, 38, 119, 221, 
           259, 310, 65)

town <- c(rep(0,8),rep(1,7))

prop.cases <- cases/totals

# fit the model using ao.glm() with a symmetric transformation
fit.ao <- ao.glm (prop.cases ~ age * town, 
                  weights = totals, phi = seq(0,0.99,0.01),
                  link = ao.sym)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
