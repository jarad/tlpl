\name{sckm}
\alias{sckm}
\title{A convenience function to create a variety of stochastic chemical kinetic models.}
\usage{
  sckm(system, ...)
}
\arguments{
  \item{system}{a character string indicating the type of
  system to create, currently implemented models are sir}

  \item{X}{an optional numeric vector of length 3 that
  determines the initial number of susceptible, infectious,
  and recovered individiuals}

  \item{theta}{an optional numberic vector of length 2 that
  determines the parameter controlling the rate at which
  individuals become infected and how quickly they recover}
}
\value{
  an sckm system
}
\description{
  A convenience function to create a variety of stochastic
  chemical kinetic models.
}
\author{
  Jarad Niemi \email{niemi@iastate.edu}
}

