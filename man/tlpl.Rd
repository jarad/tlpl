\name{tlpl}
\alias{tlpl}
\title{Performs tau-leaped particle learning}
\usage{
  tlpl(data, sckm, swarm = NULL, prior = NULL,
    n.particles = NULL, engine = "C", verbose = 0, ...)
}
\arguments{
  \item{data}{a list with elements y and tau}

  \item{sckm}{a list with a bunch of elements}

  \item{swarm}{a particle swarm}

  \item{prior}{a prior to be used for all particles}

  \item{n.particles}{a number of particles}

  \item{engine}{use 'R' or 'C' code}

  \item{verbose}{level of output}

  \item{...}{additional parameters for resampling}
}
\value{
  a list containing sample for states and sufficient
  statistics
}
\description{
  Performs tau-leaped particle learning
}
\author{
  Jarad Niemi \email{niemi@iastate.edu}
}
\seealso{
  \code{\link{tlpl_quantile}}
}

