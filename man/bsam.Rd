\name{bsem}
\alias{bsem}
\title{Fit Structural Equation Models}
\description{
Fit a Structural Equation Model (SEM).}
\usage{
bsem(..., ngibbs = 50L)
}
\arguments{
\item{...}{Default blavaan arguments.  See \code{\link{blavaan}}.}
\item{ngibbs}{Number of Gibbs sampling iterations for the structural
  model (inside each Stan sample of the measurement model). Defaults to 50.}
}
\details{
  The \code{bsam} function performs "structural after measurement"
  estimation of a \code{blavaan} model.
  }
\value{
  An object of class \code{\linkS4class{blavaan}}, which inherits from \code{\linkS4class{lavaan}}.
}
\references{
Edgar C. Merkle, Ellen Fitzsimmons, James Uanhoro, & Ben Goodrich (2021). Efficient Bayesian Structural Equation Modeling in Stan. Journal of Statistical
Software, 100(6), 1-22. URL http://www.jstatsoft.org/v100/i06/.
}
\seealso{
\code{\link{blavaan}}
}
\examples{
\dontrun{
## The industrialization and Political Democracy Example
## Bollen (1989), page 332
model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
'

## unique priors for mv intercepts; parallel chains
fit <- bsam(model, data=PoliticalDemocracy,
            dp=dpriors(nu="normal(5,10)"))
summary(fit)
}
}
