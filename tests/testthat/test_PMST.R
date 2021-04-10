set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

# Good Input
test_that("The function could record the parameters", {
    expect_snapshot(PMST.Construction(simdata, 20, 0.5))
    expect_s4_class(PMST.Construction(simdata, 20, 0.5), "dgCMatrix")
})

# Bad Input: iter
test_that("Wrong input for iter",{
    expect_error(PMST.Construction(simdata, -1, 0.5), "iter should be positive integer.")
    expect_error(PMST.Construction(simdata, 0, 0.5), "iter should be positive integer.")
    expect_error(PMST.Construction(simdata, 2.5, 0.5), "iter should be positive integer.")
})

# Bad Input: r
test_that("Wrong input for r",{
    expect_error(PMST.Construction(simdata, 20, -0.5), "The parameter r should be positive.")
})
