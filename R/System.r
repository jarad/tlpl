# Contains functions for creating stochastic chemical kinetic models or systems


# Creates a test system for use in testing functions
test.system = function(i=1)
{
  foundi = FALSE
  switch(i,
  {            # Constitutive production
    X = c(0)
    Pre  = matrix(0,1,1)
    Post = matrix(1,1,1)
    s = length(X)
    r = nrow(Pre)
    theta =  rgamma(r, 100, 100)
    stoich = t(Post-Pre)
    mult = rep(1,r)
    foundi = TRUE
  }, 
  {
    X = c(0,1,2)
    Pre = rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1),c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1))
    s = length(X)
    r = nrow(Pre)
    theta =  rgamma(r, 100, 100)
    Post = matrix(1,nrow=r,ncol=s)
    stoich = t(Post-Pre)
    mult = rep(1,r)
    foundi = TRUE
  })

  # Used in place of a default switch case
  if (!foundi) stop("Test system doesn't exist.")

  sys = list(r=r,s=s,Pre=Pre,Post=Post,stoich=stoich,theta=theta,X=X,mult=mult)
  check.system(sys)
 
  return(sys)
}


null.system = function(r=2,s=3)
{
  Pre  = matrix(0,r,s)
  Post = matrix(0,r,s)
  stoich = t(Post-Pre)
  theta = rep(0,r)
  X = rep(0,s)    
  mult = rep(1,r)
  
  sys = list(r=r,s=s,Pre=Pre,Post=Post,stoich=stoich,theta=theta,X=X,mult=mult)
  check.system(sys)
 
  return(sys)
}

# 
random.system = function(seed=1,r=NULL,s=NULL)
{
  set.seed(seed)
  if (is.null(r)) r = rpois(1,2)+1
  if (is.null(s)) s = rpois(1,3)+1
  Pre  = matrix(rpois(r*s, 1),r,s)
  Post = matrix(rpois(r*s, 1),r,s)
  stoich = t(Post-Pre)
  theta = rgamma(r,100,100)
  X = rpois(s,10)+1    
  mult = rep(1,r)

  sys = list(r=r,s=s,Pre=Pre,Post=Post,stoich=stoich,theta=theta,X=X,mult=mult)
  check.system(sys)
 
  return(sys)
}




check.system = function(sys) {
  stopifnot(sys$s == ncol(  sys$Pre),
            sys$s == ncol(  sys$Post),
            sys$s == nrow(  sys$stoich),
            sys$s == length(sys$X),

            sys$r == nrow(  sys$Pre),
            sys$r == nrow(  sys$Post),
            sys$r == ncol(  sys$stoich),
            sys$r == length(sys$theta),
            sys$r == length(sys$mult),

            all(sys$theta>=0),
            all(sys$X    >=0),
            all(sys$Pre  >=0),
            all(sys$Post >=0),
            all(sys$mult >=0))

  if (!all.equal(sys$stoich,t(sys$Post-sys$Pre))) warning("sys$stoich!=t(sys$Post-sys$Pre))")
}



#' A convenience function to create a variety of stochastic chemical kinetic models.
#'
#' @param system a character string indicating the type of system to create, currently implemented models are sir and seir
#' @param ... other parameters sent on to the system construction, e.g. the state of the system X and the rate parameters theta
#' @return an sckm system
#' @author Jarad Niemi \email{niemi@@iastate.edu}
#' @export sckm
sckm = function(system, ...) {
  do.call(tolower(system), list(...))
}


sir = function(X=c(1000,10,0),theta=c(0.5,0.25)) {
  # SIR model:
  # rxn 1: S + I -> 2I
  # rxn 2:     I -> R

               # S I R
  Pre  = rbind(c(1,1,0), # S + I 
               c(0,1,0)) # I

               # S I R
  Post = rbind(c(0,2,0), # 2I
               c(0,0,1)) # R
  stoich = t(Post-Pre)
  mult = c(1/sum(X), 1) 
  states = c("S","I","R")
  reactions = c("S+I->2I","I->R")

  sys = list(r=2, s=3, Pre=Pre, Post=Post, stoich=stoich, theta=theta, X=X, mult=mult, states=states, reactions=reactions)
  check.system(sys)
 
  return(sys)
}

seir = function(X=c(1000,0,10,0), theta=c(.5,.25,.25)) {
  # SEIR model:
  # rxn 1: S + I -> E + I
  # rxn 2:     E -> I
  # rxn 3:     I -> R

              # S E I R
  Pre  = rbind(c(1,0,1,0), # S + I
               c(0,1,0,0), # E 
               c(0,0,1,0)) # I
 
  Post = rbind(c(0,1,1,0), # E + I
               c(0,0,1,0), # I
               c(0,0,0,1)) # R

  stoich = t(Post-Pre)
  mult = c(1/sum(X), 1, 1)
  states = c("S","E","I","R")
  reactions = c("S+I->E+I","E->I","I->R")

  sys = list(r=3, s=4, Pre=Pre, Post=Post, stoich=stoich, theta=theta, X=X, mult=mult, states=states, reactions=reactions)
  check.system(sys)

  return(sys)
}
       

