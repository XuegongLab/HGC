set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

# Good Input
test_that("The function could record the parameters", {
  expect_snapshot(RNN.Construction(simdata, 2.5))
  expect_s4_class(RNN.Construction(simdata, 2.5), "dgCMatrix")
})

# Bad Input
test_that("Wrong input for dist",{
  expect_error(RNN.Construction(simdata, -1), "The threshold distance should be non-negative.")
})
