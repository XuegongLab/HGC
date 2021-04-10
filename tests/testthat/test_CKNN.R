set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

# Good Input
test_that("The function could record the parameters", {
    expect_snapshot(CKNN.Construction(simdata, 20, 1))
    expect_s4_class(CKNN.Construction(simdata, 20, 1), "dgCMatrix")
})

# Bad Input: k
test_that("Wrong input for k",{
    expect_error(CKNN.Construction(simdata, -1, 1), "k should be positive integer.")
    expect_error(CKNN.Construction(simdata, 0, 1), "k should be positive integer.")
    expect_error(CKNN.Construction(simdata, 2.5, 1), "k should be positive integer.")
    expect_error(CKNN.Construction(simdata, 10000, 1), "k should be smaller than cell number.")
})

# Bad Input: delta
test_that("Wrong input for delta",{
    expect_error(CKNN.Construction(simdata, 20, -1), "The parameter delta should be positive.")
})
