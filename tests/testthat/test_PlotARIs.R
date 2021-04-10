set.seed(2000)
ncells <- 100
nPCs <- 20
simdata <- matrix(rnorm(nPCs*ncells), ncol=nPCs, nrow=ncells)

sim.SNN <- SNN.Construction(simdata)
sim.ClusteringTree <- HGC.dendrogram(G = sim.SNN)

# Good Input
sim.labels <- data.frame(Tissue = sample(c("A","B"), 100, replace = TRUE),
                            CellLine = sample(c("AA","BA", "BB"), 100, replace = TRUE))

test_that("The function could show the figure", {
  expect_equal(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 2, k.max = 15, sim.labels), 1)
  expect_equal(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 5, k.max = 15, sim.labels), 1)
  expect_equal(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 2, k.max = 25, sim.labels), 1)
  expect_equal(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 5, k.max = 25, sim.labels), 1)
})

# Bad Input: tree
test_that("Wrong input for tree",{
    expect_error(HGC.PlotARIs(tree = 1, k.min = 2, k.max = 15, sim.labels), "Wrong input data structure for tree.")
    expect_error(HGC.PlotARIs(tree = simdata, k.min = 2, k.max = 15, sim.labels), "Wrong input data structure for tree.")
})

# Bad input: k
test_that("Wrong input for k.min or k.max",{
  expect_error(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = -1, k.max = 15, sim.labels), "k.min and k.max should be positive integer.")
  expect_error(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 1.1, k.max = 15, sim.labels), "k.min and k.max should be positive integer.")
  expect_error(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 2, k.max = 15.1, sim.labels), "k.min and k.max should be positive integer.")
  expect_error(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 2, k.max = 10000, sim.labels), "k.min and k.max should be smaller than cell number.")
  expect_error(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 15, k.max = 2, sim.labels), "k.min should not be larger than k.max.")
})

# Bad input: labels
test_that("Wrong input for labels",{
  expect_error(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 2, k.max = 15, sim.labels[-1,]), "The labels do not match the data.")
  expect_error(HGC.PlotARIs(tree = sim.ClusteringTree, k.min = 2, k.max = 15, rbind(sim.labels,sim.labels)), "The labels do not match the data.")
})
