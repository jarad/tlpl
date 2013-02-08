# Calculates the quantiles for a mixture of beta distributions
mix.betaq = function(a,b,probs,w=1)
{
  n = length(a)
  p = length(probs)
  q = rep(NA,p)
  if (w==1) w = rep(1/n,n)

  for (i in 1:length(probs)) 
  {
    tmp = uniroot(function(x) probs[i]-sum(w*pbeta(x,a,b)), c(0,1))
    q[i] = tmp$root
  }

  return(q)
}

# Calculates the quantiles for a mixture of beta distributions
mix.gammaq = function(a,b,probs,w=1)
{
  n = length(a)
  p = length(probs)
  q = rep(NA,p)
  if (w==1) w = rep(1/n,n)

  for (i in 1:length(probs)) 
  {
    tmp = uniroot(function(x) probs[i]-sum(w*pgamma(x,a,b)), c(0,10))
    q[i] = tmp$root
  }

  return(q)
}


#' This function will calculate the desired quantiles for tlpl output.
#'
#' @param tlpl output from tlpl
#' @param probs numeric vector of probabilities with values in [0,1]. 
#' @param which determines which variables quantiles should be produced for
#' @param verbose level of verbosity while running
#' @return a list with elements 'X.quantiles', 'p.quantiles', and 'r.quantiles' with dimensions s (or r) x p x n where s is the number of states, r is the number of reactions, p=length(probs), and n is the number of observations
#' @seealso tlpl
#' @export tlpl_quantile
#'
tlpl_quantile = function(tlpl, probs=c(.025,.5,.975), which="xpr", verbose=1)
{
  which = tolower(which)
  n = dim(tlpl$X)[3]
  s = dim(tlpl$X)[1]
  p = length(probs)

  do =list(x = grepl("x", which), p = grepl("p", which), r = grepl("r", which))
 
  if (do$p || do$r) { r = dim(tlpl$hyper$prob$a)[1] }
 
  if (do$x) { Xq = array(NA, dim=c(s,p,n)) } else { Xq=NULL }
  if (do$p) { pq = array(NA, dim=c(r,p,n)) } else { pq=NULL }
  if (do$r) { rq = array(NA, dim=c(r,p,n)) } else { rq=NULL }

  for (i in 1:n) 
  {   
    if (verbose) cat(paste("Time point ",i,", ",round(i/n*100), "% completed.\n", sep=''))

    if (do$x) 
    {
      for (j in 1:s) Xq[j,,i] = quantile(tlpl$X[j,,i], probs=probs)
    }

    if (do$p) 
    {
      for (j in 1:r) pq[j,,i] = mix.betaq( tlpl$hyper$prob$a[j,,i], tlpl$hyper$prob$b[j,,i], probs)
    }
        
    if (do$r)
    {
      for (j in 1:r) rq[j,,i] = mix.gammaq(tlpl$hyper$rate$a[j,,i], tlpl$hyper$rate$b[j,,i], probs)
    }
  }

  return(list(X.quantiles=Xq, p.quantiles=pq, r.quantiles=rq))
}

