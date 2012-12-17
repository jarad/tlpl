\name{tau_leap}
\alias{tau_leap}
\title{Performs tau-leaped simulations}
\usage{
  tau_leap = function(sys, tau=1, while.max=1000, engine="R")
}
\arguments{
  \item{sys}{a reaction system}

  \item{tau}{time step between observations}

  \item{while.max}{number of unsuccessful iterations to allow}

  \item{engine}{use 'R' or 'C' code}
}
\value{
  a list containing the states and the number of transitions through time

}
\description{
  Performs tau-leaped simulations
}
\author{
  Jarad Niemi \email{niemi@iastate.edu}
}
\seealso{
  \code{\link{tlpl}}
}

