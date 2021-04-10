set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

sim.SNN <- SNN.Construction(simdata)

# Good Input
test_that("The function could record the parameters", {
    expect_snapshot(HGC.parameter(G = sim.SNN))
})

# Bad Input
test_that("Wrong input for G",{
    expect_error(HGC.parameter(G = 552), "Wrong graph data structure")
    expect_error(HGC.parameter(G = simdata), "Wrong graph data structure")
    expect_error(HGC.parameter(G = sim.SNN[1:50,]), "Wrong graph data structure")
    expect_error(HGC.parameter(G = sim.SNN - 1), "Wrong graph data structure")
})
