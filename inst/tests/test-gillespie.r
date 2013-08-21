
n.reps = 9

context("gillespie")


test_that("random.system passes check.system",
{
  for (i in 1:n.reps) 
  {
    expect_true({check.system(random.system());TRUE}, info=paste("i=",i))
  }
})



########### tau_leap ###########

context("tau_leap")    
set.seed(proc.time()*1e6)
test_that("tau_leap R and C match", { 
  for (i in 1:n.reps) {
    seed = sample(1e6,1)
    sys = sir()
    n = rpois(1,1)+1
    tau = rgamma(1,100,100)
    expect_equal({set.seed(seed); tau_leap(sys, n, tau, engine="R")},
                 {set.seed(seed); tau_leap(sys, n, tau, engine="C")}, 
                 info=paste("seed= ",seed, ", n=", n, ", tau= ", tau))
  }
})


