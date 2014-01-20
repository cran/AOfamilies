\name{ao.qr}
\alias{ao.qr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic function
}
\description{
Generic function of ao.qr()
}
\usage{
ao.qr(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
For the UseMethod function
}
  \item{\dots}{
Not used
}
}
\author{
Hakim-Moulay Dehbi
}
\examples{
## labor dataset
data(labor)
data <- labor

## fit the model using ao.qr() with a symmetric transformation 
## for the quantile 0.5 with the method "linprog"
fit.ao.qr <- ao.qr(pain ~ treatment*time, data = data, 
                   kappa = 0.5, estimation = "linprog", 
                   method = "fn", phi = seq(-1,1,0.1), 
                   transfo = "ao.sym") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
