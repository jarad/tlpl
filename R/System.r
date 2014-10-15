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

#updated seir function, now takes n argument, which is a named vector
#ne is number of exposed states
#ni is number of infected states
seir = function(X=c(16000,0,100,0), theta=c(.5,2/3,2/3), n=c(ne=1,ni=1)) {
 
  ni <- n[["ni"]]
  ne <- n[["ne"]]
  
  
  #make sure X and Theta are the right length, initialize them if not.
  if(length(X)!=(ni+ne+2)){
    X=c(16000,rep(0,ne),100,rep(0,ni))
  }
  
  if(length(theta)!=(ni+ni+ne)){
    theta=theta=c(rep(.5,ni),rep((ne/3),ne),rep(ni/3,ni))
  }
 
  
  #create Pre matrix:
  #will consist of two submatrices
  #first will be all ways of getting infected, i.e., rows
  #with a 1 in the S column and a 1 in an infected state column for different infected states.
  #One row for each infected state.
  m11 <- matrix(1,nrow=ni,ncol=1)
  m12 <- matrix(0, nrow=ni, ncol=ne)
  m13 <- diag(ni)
  m14 <- matrix(0,nrow=ni,ncol=1)
  M1 <- cbind(m11,m12,m13,m14)
  
  #next is a submatrix for trasitioning from either: Ei -> E(i+1), Ei -> I1, Ii -> I(i+1), or Ii -> R
  m21 <- matrix(0,nrow=ne+ni,ncol=1)
  m22 <- diag(ne+ni)
  m23 <- matrix(0,nrow=ne+ni,ncol=1)
  M2 <- cbind(m21,m22,m23)
  
  Pre = rbind(M1, M2)
  
  
  #create Pre matrix:
  #will consist of two submatrices
  #first will be result of infection from each different I state, i.e.
  #rows with a 1 in the E1 column (or I1 column if no exposed state) and a 1 in an infected state column.  
  #One row for each infected state.
  if(ne>0){
    m31 <- matrix(0,nrow=ni,ncol=1)
    m32 <- matrix(1,nrow=ni,ncol=1)
    m33 <- matrix(0, nrow=ni,ncol=ne-1)
    m34 <- diag(ni)
    m35 <- matrix(0, nrow=ni, ncol=1)
    M3 <- cbind(m31,m32,m33,m34,m35)
  }
  if(ne==0){
    m31 <- matrix(0,nrow=ni,ncol=1)
    m32 <- matrix(c(2,rep(1,ni-1)),nrow=ni,ncol=1)
    m33 <- rbind(rep(0,ni-1),diag(ni-1))
    m34 <- matrix(0, nrow=ni, ncol=1)
    M3 <- cbind(m31,m32,m33,m34) 
  }
  
  #next submatrix shows results of moving either E-> E, E->I, I->I, or I-> R
  m41<- matrix(0,nrow=ne+ni,ncol=2)
  m42 <- diag(ne+ni)
  M4 <- cbind(m41,m42)
  
  Post <- rbind(M3,M4)
  
  stoich = t(Post-Pre)
  mult = c(rep(1/sum(X),ni),rep(1,ne+ni))
  
  
  #Create list of different states, and then combine them into different reactions.
  if(ne>0){
    numbers<-c(seq(1,ne),seq(1,ni))
    letters<-c(rep("E",ne), rep("I",ni))
  }
  if(ne==0){
    numbers<-c(seq(1,ni))
    letters<-c(rep("I",ni))
  }
  states<-c("S",paste(letters, numbers,sep=""),"R")
  
  nreactions <- dim(Pre)[1]
  reactions <- rep(0,nreactions)
  for(i in 1:nreactions){
    reactions[i] <- paste(paste( states[which(Pre[i,]!=0)],collapse="+"),"->",paste( states[which(Post[i,]!=0)],collapse="+"),sep="")
  }
  
  sys = list(r=nreactions, s=length(states), Pre=Pre, Post=Post, stoich=stoich, theta=theta, X=X, mult=mult, states=states, reactions=reactions)
  check.system(sys)
  return(sys)
}
