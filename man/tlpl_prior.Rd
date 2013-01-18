\name{tlpl_prior}
\alias{tlpl_prior}
\title{A convenience function to make a tlpl prior}
\usage{
  tlpl_prior(X, p.a, p.b, r.a, r.b, nr)
}
\arguments{
  \item{X}{a current state}

  \item{p.a}{vector of alpha parameters for Beta( alpha,
  beta) distributions for probabilities}

  \item{p.b}{vector of beta parameters for Beta( alpha,
  beta) distributions for probabilities}

  \item{r.a}{vector of alpha parameters for Gamma(alpha,
  beta) distributions for rates}

  \item{r.b}{vector of beta parameters for Gamma(alpha,
  beta) distributions for rates}

  \item{nr}{scalar indicating the number of particles}
}
\value{
  a list containing all parameters in a list
}
\description{
  A convenience function to make a tlpl prior
}
\author{
  Jarad Niemi \email{niemi@iastate.edu}
}
\seealso{
  \code{\link{tlpl}}
}

