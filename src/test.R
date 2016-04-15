library(Rcpp)

sourceCpp("src/resample.cpp")

resample(runif(5),10, 2)


sourceCpp("src/gillespie.cpp")

# SIR
library(tlpl)
sir = sckm('sir')
sir$lmult = log(sir$mult)

sir$X
update_species(sir, c(1, 1))
sir$X

# 
hazard_part(sir, sir$X)

n = 10:15
k = 1:6

lchoose(n,k)

test(n,k)
