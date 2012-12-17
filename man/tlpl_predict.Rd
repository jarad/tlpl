\name{tlpl_predict}
\alias{tlpl_predict}
\title{Produces simulations of data based on a current tlpl filtered distribution}
\usage{
  tlpl_predict(sckm, n, tlpl, engine = "R", verbose)
}
\arguments{
  \item{sckm}{a list defining a stochastic chemical kinetic
  model}

  \item{n}{an integer determining how many forward points
  to simulate}

  \item{tlpl}{a list defining a filtered distribution}

  \item{engine}{use 'R' or 'C' code (only R currently
  implemented)}

  \item{verbose}{level of verbosity curing calculations}
}
\value{
  list defining the states and the number of transitions
}
\description{
  Produces simulations of data based on a current tlpl
  filtered distribution
}
\author{
  Jarad Niemi \email{niemi@iastate.edu}
}
\seealso{
  \code{\link{tlpl}}
}

