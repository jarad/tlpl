\name{tlpl_quantile}
\alias{tlpl_quantile}
\title{This function will calculate the desired quantiles for tlpl output.}
\usage{
  tlpl_quantile(tlpl, probs = c(0.025, 0.5, 0.975),
    which = "xpr", verbose = 1)
}
\arguments{
  \item{tlpl}{output from tlpl}

  \item{probs}{numeric vector of probabilities with values
  in [0,1].}

  \item{which}{determines which variables quantiles should
  be produced for}

  \item{verbose}{level of verbosity while running}
}
\description{
  This function will calculate the desired quantiles for
  tlpl output.
}
\seealso{
  tlpl
}

