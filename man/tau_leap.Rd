\name{tau_leap}
\alias{tau_leap}
\title{Performs tau-leaped simulations}
\usage{
  tau_leap(sys, n = 1, tau = 1, while.max = 1000,
    engine = "R")
}
\arguments{
  \item{sys}{a list defining a stochastic chemical kinetic
  system}

  \item{n}{an integer defining the number of time-points to
  simulate}

  \item{tau}{a positive vector defining the times between
  observations}

  \item{while.max}{at each time point the maximum number of
  simulations to try to ensure non-negativity of all
  species}

  \item{engine}{use 'R' or 'C'}
}
\value{
  a list containing the species counts at each time point
  ('X') which is an (n+1) x s matrix and the number of
  transitions between each time point ('nr') which is an n
  x r matrix
}
\description{
  Performs tau-leaped simulations
}
\author{
  Jarad Niemi \email{niemi@iastate.edu}
}

