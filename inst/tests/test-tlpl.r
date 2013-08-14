context("tlpl")

test_that("tlpl throws errors", {
    expect_error(tlpl())
    expect_error(tlpl(1:6))
    expect_error(tlpl(1:6,1:6))
})


test_that("tlpl completes", {
    sys = test.system(1)
    dat = list()
    dat$y = matrix(1:6, 6, 1)
    dat$tau = 1
    expect_true({tlpl(dat,sys); TRUE})
})


test_that("tlpl: constitutive production", {
    sys = test.system(1)
    dat = list()
    n = 6
    dat$y = matrix(1, n, 1)
    dat$tau = 1
    res = tlpl(dat,sys,nonuniformity="ess",threshold=0)
    np = 2
    expect_equal(res$hyper$prob$a, array(rep(1:(n+1),each=2), dim=c(1,np,n+1)))
})




