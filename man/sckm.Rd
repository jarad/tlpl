\name{sckm}
\alias{sckm}
\title{A convenience function to create a variety of stochastic chemical kinetic models.}
\usage{
  sckm(system, ...)
}
\arguments{
  \item{system}{a character string indicating the type of
  system to create, currently implemented models are sir
  and seir}

  \item{...}{other parameters sent on to the system
  construction, e.g. the state of the system X and the rate
  parameters theta}
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

