set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

# Good Input
test_that("The function could record the parameters", {
    expect_snapshot(MST.Construction(simdata))
    expect_s4_class(MST.Construction(simdata), "dgCMatrix")
})
