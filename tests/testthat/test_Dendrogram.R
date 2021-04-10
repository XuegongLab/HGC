set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

sim.SNN <- SNN.Construction(simdata)

# Good Input
test_that("The function could record the parameters", {
    expect_snapshot(HGC.dendrogram(G = sim.SNN))
    expect_s3_class(HGC.dendrogram(G = sim.SNN), "hclust")
})

# Bad Input
test_that("Wrong input for G",{
  expect_error(HGC.dendrogram(G = 552), "Wrong graph data structure")
  expect_error(HGC.dendrogram(G = simdata), "Wrong graph data structure")
  expect_error(HGC.dendrogram(G = sim.SNN[1:50,]), "Wrong graph data structure")
  expect_error(HGC.dendrogram(G = sim.SNN - 1), "Wrong graph data structure")
})
