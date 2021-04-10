set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

# Good Input
test_that("The function could record the parameters", {
    expect_snapshot(SNN.Construction(simdata, 20, 1/15))
    expect_s4_class(SNN.Construction(simdata, 20, 1/15), "dgCMatrix")
})

# Bad Input: k
test_that("Wrong input for k",{
    expect_error(SNN.Construction(simdata, -1, 1/15), "k should be positive integer.")
    expect_error(SNN.Construction(simdata, 0, 1/15), "k should be positive integer.")
    expect_error(SNN.Construction(simdata, 2.5, 1/15), "k should be positive integer.")
    expect_error(SNN.Construction(simdata, 10000, 1/15), "k should be smaller than cell number.")
})

# Bad Input: threshold
test_that("Wrong input for threshold",{
    expect_error(SNN.Construction(simdata, 20, -1), "The parameter threshold should be within zero to one.")
    expect_error(SNN.Construction(simdata, 20, 2), "The parameter threshold should be within zero to one.")
})
