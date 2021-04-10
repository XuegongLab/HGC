set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

# Good Input
test_that("The function could record the parameters", {
    expect_snapshot(KNN.Construction(simdata, 20))
    expect_s4_class(KNN.Construction(simdata, 20), "dgCMatrix")
})

# Bad Input
test_that("Wrong input for k",{
    expect_error(KNN.Construction(simdata, -1), "k should be positive integer.")
    expect_error(KNN.Construction(simdata, 0), "k should be positive integer.")
    expect_error(KNN.Construction(simdata, 2.5), "k should be positive integer.")
    expect_error(KNN.Construction(simdata, 10000), "k should be smaller than cell number.")
})
