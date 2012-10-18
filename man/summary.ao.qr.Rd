\name{summary.ao.qr}
\alias{summary.ao.qr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing the results of an "ao.qr" object
}
\description{
Summary method of the function \code{\link{ao.qr}}
}
\usage{
\method{summary}{ao.qr} (object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{
An object of class "ao.qr"
}
  \item{...}{
Not used
}
}

\examples{
## labor dataset
data(labor)
data <- labor

## fit the model using ao.qr() with a symmetric transformation  
## for the quantile 0.25 with the method "linprog"
fit.ao.qr <- ao.qr(pain ~ treatment*time, data = data, 
                   kappa = 0.25, estimation = "linprog", 
				   transfo = "ao.sym")

## summarize the model and assign the summary to a "summary.ao.qr" object 				   
summary.fit.ao.qr <- summary.ao.qr (fit.ao.qr)					
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
