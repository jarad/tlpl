#' Produces simulations of data based on a current tlpl filtered distribution 
#' 
#' @param sckm a list defining a stochastic chemical kinetic model
#' @param n an integer determining how many forward points to simulate
#' @param tlpl a list defining a filtered distribution
#' @param engine use 'R' or 'C' code (only R currently implemented)
#' @param verbose level of verbosity curing calculations
#' @return list defining the states and the number of transitions
#' @seealso \code{\link{tlpl}}
#' @author Jarad Niemi \email{niemi@@iastate.edu}
#' @export tlpl_predict
#'
tlpl_predict = function(sckm, n, tlpl, engine="R", verbose) {
  stopifnot(engine=="R")

  np = ncol(tlpl$X)
  nr = length(sckm$theta)
  ns = length(sckm$X)
	
  X = array(0, dim=c(ns, np, n+1))
  y = array(0, dim=c(nr, np, n))
    for (i in 1:np) {
      if (verbose>0 && i%%100==0) cat("Particle",i,"(",i/np*100,"%)\n") 
        sckm$X = tlpl$X[,i]
        sckm$theta = rgamma(nr, tlpl$hyper$rate$a[,i],
		                tlpl$hyper$rate$b[,i])
        tl = tau.leap(sckm, n, engine=engine)
        X[,i,] = t(tl$X)
		
        p = rbeta(nr, tlpl$hyper$prob$a[,i], 
		      tlpl$hyper$prob$b[,i])

        for (j in 1:nr) 
        {
          y[j,i,] = rbinom(n, tl$nr[,j], p[j])
        }
    }
  return(list(X=X,y=y))
}

