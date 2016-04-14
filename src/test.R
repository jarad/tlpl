library(Rcpp)

sourceCpp("src/resample.cpp")

resample(runif(5),10, 2)


sourceCpp("src/gillespie.cpp")

# SIR
A = cbind(c(-1,1,0), c(0,-1,1))
S = t(A)

update_species(S, c(1,1), c(2,3,0))

